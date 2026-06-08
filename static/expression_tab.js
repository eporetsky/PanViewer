/**
 * Wheat pangene detail: Expression tab — PlantApp tissue + DEG (plotly).
 * Expects globals from pangene_detail: GID, ACC, refGene (optional).
 */
(function () {
    const CS = 'chinesespring';
    let plotlyPromise = null;
    /** Last successfully loaded gene id (skip duplicate fetch on tab revisit). */
    let lastFetchedGeneId = null;
    /** Selected Chinese Spring gene when multiple exist in the cluster. */
    let activeCsGeneId = null;
    let csGeneTabsWired = false;

    function panelEl() {
        return document.getElementById('expression-panel');
    }

    function accLower(gid) {
        if (typeof ACC === 'undefined' || !gid) return '';
        const a = ACC[gid];
        return a != null ? String(a).toLowerCase().trim() : '';
    }

    function parseServerCsGenes() {
        const p = panelEl();
        if (!p || !p.dataset.csGenes) return [];
        try {
            const arr = JSON.parse(p.dataset.csGenes);
            return Array.isArray(arr) ? arr.filter(Boolean) : [];
        } catch (e) {
            return [];
        }
    }

    /**
     * Union: CS genes in the current alignment matrix plus server-listed pangene members (sorted).
     * Organ/tissue labels come from PlantApp sample metadata (`organ` on each record).
     */
    function findAllChineseSpringGeneIds() {
        const fromServer = parseServerCsGenes();
        const fromAlign = [];
        if (typeof GID !== 'undefined' && Array.isArray(GID)) {
            for (let i = 0; i < GID.length; i++) {
                const g = GID[i];
                if (accLower(g) === CS) fromAlign.push(g);
            }
        }
        const seen = Object.create(null);
        const out = [];
        function add(g) {
            if (!g || seen[g]) return;
            seen[g] = true;
            out.push(g);
        }
        fromAlign.forEach(add);
        fromServer.forEach(add);
        out.sort();
        return out;
    }

    function ensurePlotly() {
        if (window.Plotly) return Promise.resolve();
        if (plotlyPromise) return plotlyPromise;
        plotlyPromise = new Promise(function (resolve, reject) {
            const s = document.createElement('script');
            s.src = 'https://cdn.plot.ly/plotly-2.35.2.min.js';
            s.async = true;
            s.onload = function () { resolve(); };
            s.onerror = function () {
                plotlyPromise = null;
                reject(new Error('Plotly CDN failed'));
            };
            document.head.appendChild(s);
        });
        return plotlyPromise;
    }

    function compactTissueToRecords(tissue) {
        if (!tissue || tissue.format !== 'compact') return [];
        const samples = tissue.samples || [];
        const cpm = tissue.cpm_by_sample || {};
        return samples.map(function (s) {
            const sid = s.sample_acc;
            const v = sid != null ? cpm[sid] : null;
            return Object.assign({}, s, { value: v === undefined ? null : v });
        });
    }

    /** Soft-wrap hover fields: max characters per line (word-aware where possible). */
    var HOVER_WRAP_MAX_CHARS = 100;
    /** Also start a new line after this many words (whichever yields a shorter line first). */
    var HOVER_WRAP_MAX_WORDS = 14;

    /** Safe text inside Plotly hover HTML. */
    function escapeHoverText(s) {
        return String(s || '')
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;');
    }

    /**
     * Split plain text into lines: break before exceeding maxChars, or maxWords per line,
     * whichever comes first on each line. Long single tokens are hard-broken at maxChars.
     */
    function wrapPlainToLines(s, maxChars, maxWords) {
        maxChars = maxChars || HOVER_WRAP_MAX_CHARS;
        maxWords = maxWords == null ? HOVER_WRAP_MAX_WORDS : maxWords;
        s = String(s || '').trim();
        if (!s) return [];
        const words = s.split(/\s+/);
        const lines = [];
        let cur = '';
        let curWordCount = 0;
        function flush() {
            if (cur) {
                lines.push(cur);
                cur = '';
                curWordCount = 0;
            }
        }
        for (let wi = 0; wi < words.length; wi++) {
            const w = words[wi];
            const next = cur ? cur + ' ' + w : w;
            const nextWords = cur ? curWordCount + 1 : 1;
            const overChars = next.length > maxChars;
            const overWords = maxWords > 0 && nextWords > maxWords;
            if (overChars || overWords) {
                if (cur) {
                    flush();
                }
                if (w.length > maxChars) {
                    for (let i = 0; i < w.length; i += maxChars) {
                        lines.push(w.slice(i, i + maxChars));
                    }
                } else {
                    cur = w;
                    curWordCount = 1;
                }
            } else {
                cur = next;
                curWordCount = nextWords;
            }
        }
        if (cur) {
            lines.push(cur);
        }
        return lines;
    }

    /** Escaped HTML for one hover field value with soft wraps as <br>. */
    function wrapHoverField(s) {
        return wrapPlainToLines(s, HOVER_WRAP_MAX_CHARS).map(escapeHoverText).join('<br>');
    }

    /** Optional <b>Title</b> line for hover when `title` exists (PlantApp / API). */
    function titleHoverBlock(title) {
        const t = String(title || '').trim();
        if (!t) return '';
        const body = wrapPlainToLines(t, HOVER_WRAP_MAX_CHARS).map(escapeHoverText).join('<br>');
        return '<b>Title</b>: ' + body + '<br>';
    }

    /** DEG y-axis: anchor at 0 when all values are one-sided (PlantApp logFC). */
    function degYAxisRangeFromValues(yVals) {
        const ys = yVals.filter(function (v) { return Number.isFinite(Number(v)); }).map(Number);
        if (!ys.length) return null;
        let lo = Math.min.apply(null, ys);
        let hi = Math.max.apply(null, ys);
        const span = hi - lo;
        const pad = span > 0 ? span * 0.08 : 0.5;
        const allNonNeg = ys.every(function (v) { return v >= 0; });
        const allNonPos = ys.every(function (v) { return v <= 0; });
        if (allNonNeg && allNonPos) {
            return [-1, 1];
        }
        if (allNonNeg) {
            return [0, hi + pad];
        }
        if (allNonPos) {
            return [lo - pad, 0];
        }
        return [lo - pad, hi + pad];
    }

    function tissuePlotLayout(organs) {
        const n = organs.length;
        const xaxis = {
            title: 'Tissue / organ',
            tickmode: 'array',
            tickvals: n ? organs.map(function (_, i) { return i; }) : [],
            ticktext: n ? organs.slice() : [],
            tickangle: -35,
            automargin: true,
            zeroline: false,
            showgrid: true,
        };
        if (n > 0) {
            xaxis.range = [-0.5, n - 0.5];
        }
        return {
            title: {
                text: 'Tissue-specific expression (aligned to Chinese Spring)',
                font: { size: 14, color: '#212529' },
            },
            paper_bgcolor: '#fff',
            plot_bgcolor: '#fafafa',
            font: { color: '#495057' },
            xaxis: xaxis,
            yaxis: {
                title: 'Expression (CPM)',
                automargin: true,
                zeroline: false,
                showline: false,
                showgrid: true,
            },
            margin: { l: 56, r: 24, t: 56, b: 120 },
            height: 400,
            showlegend: false,
            hovermode: 'closest',
            hoverlabel: { align: 'left' },
        };
    }

    function degPlotLayout(axisCategories) {
        const n = axisCategories.length;
        const xaxis = {
            title: 'Category',
            tickmode: 'array',
            tickvals: n ? axisCategories.map(function (_, i) { return i; }) : [],
            ticktext: n ? axisCategories.slice() : [],
            tickangle: -40,
            automargin: true,
            zeroline: false,
        };
        if (n > 0) {
            xaxis.range = [-0.5, n - 0.5];
        }
        return {
            title: { text: 'Differential expression by category', font: { size: 14, color: '#212529' } },
            paper_bgcolor: '#fff',
            plot_bgcolor: '#fafafa',
            font: { color: '#495057' },
            xaxis: xaxis,
            yaxis: {
                title: { text: 'Log\u2082(FC)' },
                zeroline: true,
                zerolinecolor: 'rgba(45,106,79,0.45)',
                automargin: true,
            },
            margin: { l: 56, r: 24, t: 48, b: 160 },
            height: 420,
            showlegend: false,
            hovermode: 'closest',
            hoverlabel: { align: 'left' },
        };
    }

    function fmtFdr(x) {
        if (x == null || x === '') return '';
        const n = Number(x);
        if (!Number.isFinite(n)) return String(x);
        if (n === 0 || n < 0.01) return n.toExponential(2);
        return n.toFixed(n < 1 ? 3 : 2);
    }

    function buildTissueFigure(records) {
        const organs = [];
        const seen = Object.create(null);
        records.forEach(function (r) {
            const o = r.organ;
            if (o != null && o !== '' && !seen[o]) {
                seen[o] = true;
                organs.push(o);
            }
        });
        const organToI = {};
        organs.forEach(function (o, i) { organToI[o] = i; });

        const x = [];
        const y = [];
        const custom = [];
        records.forEach(function (r) {
            const organ = r.organ;
            const val = r.value;
            if (organ == null || val == null || !Number.isFinite(Number(val))) return;
            const i = organToI[organ];
            if (i === undefined) return;
            x.push(i + (Math.random() - 0.5) * 0.08);
            y.push(Number(val));
            custom.push([
                wrapHoverField(r.experiment_acc || ''),
                wrapHoverField(r.short_title || ''),
                wrapHoverField(r.sample_acc || ''),
                wrapHoverField(r.group || ''),
                wrapHoverField(r.genome || ''),
                wrapHoverField(r.species || ''),
                wrapHoverField(String(organ)),
                wrapHoverField(r.stage || ''),
                wrapHoverField(r.inbred || ''),
                titleHoverBlock(r.title),
            ]);
        });

        const traces = [];
        if (x.length) {
            traces.push({
                type: 'scatter',
                mode: 'markers',
                x: x,
                y: y,
                marker: { size: 7, color: '#2d6a4f', opacity: 0.72, line: { width: 0 } },
                customdata: custom,
                hovertemplate:
                    '<b>Organ</b>: %{customdata[6]}<br>' +
                    '<b>CPM</b>: %{y:.2f}<br>' +
                    '<b>Experiment</b>: %{customdata[0]}<br>' +
                    '<b>Short title</b>: %{customdata[1]}<br>' +
                    '%{customdata[9]}' +
                    '<b>Sample</b>: %{customdata[2]}<br>' +
                    '<b>Group</b>: %{customdata[3]}<br>' +
                    '<b>Stage</b>: %{customdata[7]}<br>' +
                    '<extra></extra>',
            });
        }

        organs.forEach(function (organ, idx) {
            const vals = records
                .filter(function (r) { return r.organ === organ && r.value != null && Number.isFinite(Number(r.value)); })
                .map(function (r) { return Number(r.value); });
            if (!vals.length) return;
            const mean = vals.reduce(function (a, b) { return a + b; }, 0) / vals.length;
            traces.push({
                type: 'scatter',
                mode: 'lines',
                x: [idx - 0.22, idx + 0.22],
                y: [mean, mean],
                line: { color: '#d95f02', width: 2 },
                hoverinfo: 'skip',
                showlegend: false,
            });
        });

        return { data: traces, layout: tissuePlotLayout(organs) };
    }

    function buildDegFigure(deg, experiments) {
        const records = deg && deg.records != null ? deg.records : [];
        const expMap = experiments || {};

        /** Full x-axis category list from API (includes categories with no DEG for this gene). */
        let axisCategories = (deg && Array.isArray(deg.categories) && deg.categories.length)
            ? deg.categories.slice()
            : [];

        const catToI = {};
        axisCategories.forEach(function (c, i) {
            const key = c != null ? String(c) : '';
            if (key) {
                catToI[key] = i;
            }
        });

        const withoutDeg = (deg && Array.isArray(deg.categories_without_deg_for_gene))
            ? deg.categories_without_deg_for_gene
            : [];
        withoutDeg.forEach(function (c) {
            const key = c != null ? String(c) : '';
            if (key && catToI[key] === undefined) {
                catToI[key] = axisCategories.length;
                axisCategories.push(key);
            }
        });

        records.forEach(function (r) {
            const cat = (r.category != null) ? String(r.category) : '';
            if (cat && catToI[cat] === undefined) {
                catToI[cat] = axisCategories.length;
                axisCategories.push(cat);
            }
        });

        if (axisCategories.length === 0 && records.length) {
            const seen = Object.create(null);
            records.forEach(function (r) {
                const c = (r.category != null) ? String(r.category) : '';
                if (c && !seen[c]) {
                    seen[c] = true;
                    axisCategories.push(c);
                }
            });
            axisCategories.forEach(function (c, i) { catToI[c] = i; });
        }

        const x = [];
        const y = [];
        const custom = [];
        records.forEach(function (r) {
            const cat = (r.category != null) ? String(r.category) : '';
            const i = catToI[cat];
            if (i === undefined) return;
            x.push(i + (Math.random() - 0.5) * 0.08);
            y.push(Number(r.logFC));
            const eid = r.experiment_acc || '';
            const et = expMap[eid] || {};
            const titleLine = titleHoverBlock(et.title);
            custom.push([
                wrapHoverField(cat),
                wrapHoverField(eid),
                wrapHoverField(et.short_title || ''),
                wrapHoverField(r.comparison || ''),
                wrapHoverField(r.group1 || ''),
                wrapHoverField(r.group2 || ''),
                wrapHoverField(fmtFdr(r.FDR)),
                titleLine,
            ]);
        });

        const traces = [];
        if (x.length) {
            traces.push({
                type: 'scatter',
                mode: 'markers',
                x: x,
                y: y,
                marker: { size: 8, color: '#2d6a4f', opacity: 0.82, line: { width: 0 } },
                customdata: custom,
                hovertemplate:
                    '<b>Category</b>: %{customdata[0]}<br>' +
                    '<b>log\u2082(FC)</b>: %{y:.2f}<br>' +
                    '<b>FDR</b>: %{customdata[6]}<br>' +
                    '<b>Experiment</b>: %{customdata[1]}<br>' +
                    '<b>Short title</b>: %{customdata[2]}<br>' +
                    '%{customdata[7]}' +
                    '<b>Comparison</b>: %{customdata[3]}<br>' +
                    '<extra></extra>',
            });
        } else if (axisCategories.length) {
            traces.push({
                type: 'scatter',
                mode: 'markers',
                x: axisCategories.map(function (_, i) { return i; }),
                y: axisCategories.map(function () { return 0; }),
                marker: { opacity: 0, size: 3 },
                hoverinfo: 'skip',
                showlegend: false,
            });
        }

        const layout = degPlotLayout(axisCategories.length ? axisCategories : []);
        if (!axisCategories.length) {
            layout.xaxis.tickvals = [0];
            layout.xaxis.ticktext = ['(no categories)'];
            layout.xaxis.range = [-0.5, 0.5];
        }

        const yRange = degYAxisRangeFromValues(y);
        if (yRange) {
            layout.yaxis.range = yRange;
            layout.yaxis.autorange = false;
        } else if (axisCategories.length && !x.length) {
            layout.yaxis.range = [-1, 1];
            layout.yaxis.autorange = false;
        }

        return { data: traces, layout: layout };
    }

    function setBanner(html, kind) {
        const el = document.getElementById('expression-banner');
        if (!el) return;
        el.classList.remove('d-none', 'alert-warning', 'alert-info', 'alert-danger', 'alert-secondary');
        if (!html) {
            el.classList.add('d-none');
            el.innerHTML = '';
            return;
        }
        el.classList.add(kind || 'alert-secondary');
        el.innerHTML = html;
    }

    function apiUrl() {
        const p = panelEl();
        return p && p.dataset.plantappApi ? p.dataset.plantappApi : '';
    }

    async function fetchOmics(geneId) {
        const base = apiUrl();
        if (!base) throw new Error('Missing API URL');
        const url = base + (base.indexOf('?') >= 0 ? '&' : '?') + 'gene_id=' + encodeURIComponent(geneId);
        const res = await fetch(url, { credentials: 'same-origin' });
        if (!res.ok) throw new Error('HTTP ' + res.status);
        return res.json();
    }

    function renderEmptyTissue(msg) {
        const el = document.getElementById('expression-tissue');
        if (!el) return;
        el.innerHTML = '<p class="text-muted small p-3 mb-0">' + msg + '</p>';
    }

    function renderEmptyDeg(msg) {
        const el = document.getElementById('expression-deg');
        if (!el) return;
        el.innerHTML = '<p class="text-muted small p-3 mb-0">' + msg + '</p>';
    }

    function renderCsGeneTabs(ids) {
        const wrap = document.getElementById('expression-gene-tabs-wrap');
        const ul = document.getElementById('expression-gene-tabs');
        if (!wrap || !ul) return;

        if (ids.length <= 1) {
            wrap.classList.add('d-none');
            ul.replaceChildren();
            return;
        }

        wrap.classList.remove('d-none');
        if (!activeCsGeneId || ids.indexOf(activeCsGeneId) < 0) {
            if (typeof refGene !== 'undefined' && refGene && ids.indexOf(refGene) >= 0) {
                activeCsGeneId = refGene;
            } else {
                activeCsGeneId = ids[0];
            }
        }

        ul.replaceChildren();
        ids.forEach(function (gid) {
            const li = document.createElement('li');
            li.className = 'nav-item';
            li.setAttribute('role', 'presentation');
            const btn = document.createElement('button');
            btn.type = 'button';
            btn.className = 'nav-link py-1 px-2 small font-monospace' + (gid === activeCsGeneId ? ' active' : '');
            btn.setAttribute('role', 'tab');
            btn.setAttribute('data-gene-id', gid);
            btn.textContent = gid;
            li.appendChild(btn);
            ul.appendChild(li);
        });
    }

    function ensureCsGeneTabsWired() {
        if (csGeneTabsWired) return;
        const ul = document.getElementById('expression-gene-tabs');
        if (!ul) return;
        csGeneTabsWired = true;
        ul.addEventListener('click', function (ev) {
            const btn = ev.target.closest('button[data-gene-id]');
            if (!btn) return;
            const gid = btn.getAttribute('data-gene-id');
            if (!gid || gid === activeCsGeneId) return;
            activeCsGeneId = gid;
            ul.querySelectorAll('button[data-gene-id]').forEach(function (b) {
                b.classList.toggle('active', b.getAttribute('data-gene-id') === activeCsGeneId);
            });
            lastFetchedGeneId = null;
            void loadAndRender();
        });
    }

    async function loadAndRender() {
        const ids = findAllChineseSpringGeneIds();

        if (ids.length === 0) {
            lastFetchedGeneId = null;
            activeCsGeneId = null;
            renderCsGeneTabs(ids);
            setBanner(
                '<i class="bi bi-info-circle me-1"></i>No <strong>Chinese Spring</strong> gene (<code>ChineseSpring</code> accession) is present in this pan cluster. Expression plots are not available.',
                'alert-info'
            );
            renderEmptyTissue('No Chinese Spring gene in this pangene.');
            renderEmptyDeg('No Chinese Spring gene in this pangene.');
            return;
        }

        if (ids.length === 1) {
            activeCsGeneId = ids[0];
        }

        renderCsGeneTabs(ids);
        ensureCsGeneTabsWired();

        const gid = activeCsGeneId;
        if (!gid) return;

        if (gid === lastFetchedGeneId) return;

        const tissueEl = document.getElementById('expression-tissue');
        const degEl = document.getElementById('expression-deg');
        if (!tissueEl || !degEl) return;

        lastFetchedGeneId = gid;

        setBanner('', null);
        tissueEl.innerHTML = '<p class="text-muted small p-3 mb-0"><span class="spinner-border spinner-border-sm me-2"></span>Loading expression from PlantApp…</p>';
        degEl.innerHTML = '';

        try {
            await ensurePlotly();
        } catch (e) {
            lastFetchedGeneId = null;
            setBanner('<i class="bi bi-exclamation-triangle me-1"></i>Could not load chart library.', 'alert-danger');
            renderEmptyTissue('Chart library failed to load.');
            renderEmptyDeg('');
            return;
        }

        let data;
        try {
            data = await fetchOmics(gid);
        } catch (e) {
            lastFetchedGeneId = null;
            setBanner(
                '<i class="bi bi-wifi-off me-1"></i>Could not reach PanViewer or PlantApp for expression data. Try again later.',
                'alert-warning'
            );
            renderEmptyTissue('Request failed (network or server).');
            renderEmptyDeg('Request failed (network or server).');
            return;
        }

        if (data.unknown_gene) {
            setBanner(
                '<i class="bi bi-search me-1"></i>This gene ID is not in PlantApp’s sequence index, so omics data is unavailable (same as an unknown gene on PlantApp’s gene page).',
                'alert-info'
            );
            renderEmptyTissue('Gene not found in PlantApp.');
            renderEmptyDeg('Gene not found in PlantApp.');
            return;
        }

        if (data.error && !data.tissue && !data.deg) {
            lastFetchedGeneId = null;
            const rl = data.rate_limited
                ? '<i class="bi bi-hourglass-split me-1"></i>Rate limited — retry shortly. '
                : '';
            setBanner(
                rl + '<i class="bi bi-exclamation-triangle me-1"></i>' + String(data.error),
                'alert-warning'
            );
            renderEmptyTissue('PlantApp request failed.');
            renderEmptyDeg('');
            return;
        }

        const resolved = data.resolved_gene_id || gid;
        const plantAppGeneUrl = 'https://www.plantapp.org/gene/?gene_id=' + encodeURIComponent(resolved);
        const openPlantAppBtn =
            '<a class="btn btn-primary btn-sm" href="' +
            plantAppGeneUrl +
            '" target="_blank" rel="noopener noreferrer">' +
            '<i class="bi bi-box-arrow-up-right me-1" aria-hidden="true"></i>Open in PlantApp</a>';
        const rateFrag = data.rate_limited
            ? '<div class="small text-warning mb-2"><i class="bi bi-hourglass-split me-1"></i>Part of the PlantApp request may have been rate-limited (e.g. DEG).</div>'
            : '';
        setBanner(
            rateFrag +
                '<div class="d-flex flex-wrap align-items-center gap-2">' +
                openPlantAppBtn +
                '</div>',
            'alert-secondary'
        );

        const tissuePayload = data.tissue;
        const records = tissuePayload ? compactTissueToRecords(tissuePayload) : [];

        if (!records.length) {
            renderEmptyTissue(
                data.tissue
                    ? 'No tissue-included RNA samples with values for this gene (or all values missing).'
                    : 'No tissue expression payload returned.'
            );
        } else {
            const fig = buildTissueFigure(records);
            tissueEl.innerHTML = '';
            const div = document.createElement('div');
            div.style.width = '100%';
            tissueEl.appendChild(div);
            await Plotly.newPlot(div, fig.data, fig.layout, { responsive: true, displaylogo: false });
        }

        if (data.deg_error) {
            renderEmptyDeg('<strong>DEG</strong>: ' + String(data.deg_error));
        } else {
            const deg = data.deg;
            if (!deg) {
                renderEmptyDeg('No differential expression payload returned.');
            } else if (deg.error) {
                renderEmptyDeg('<strong>DEG</strong>: ' + String(deg.error));
            } else {
                const experiments = deg.experiments || {};
                const fig = buildDegFigure(deg, experiments);
                degEl.innerHTML = '';
                const div = document.createElement('div');
                div.style.width = '100%';
                degEl.appendChild(div);
                if (fig.data.length) {
                    await Plotly.newPlot(div, fig.data, fig.layout, { responsive: true, displaylogo: false });
                } else {
                    div.innerHTML =
                        '<p class="text-muted small p-3 mb-0">No DEG points for this gene. Categories defined for this genome may still appear on the axis when data exists.</p>';
                }
            }
        }
    }

    function onShown() {
        const p = panelEl();
        if (!p) return;
        void loadAndRender();
    }

    function invalidateCache() {
        lastFetchedGeneId = null;
        activeCsGeneId = null;
    }

    window.PBExpressionTab = {
        onShown: onShown,
        invalidateCache: invalidateCache,
        findAllChineseSpringGeneIds: findAllChineseSpringGeneIds,
    };
})();
