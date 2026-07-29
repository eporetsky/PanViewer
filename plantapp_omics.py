"""
Server-side gene omics JSON: local ``expression/expression.db`` when present,
otherwise PlantApp ``GET /api/gene-omics``.

Local and remote payloads share the same contract (compact tissue with
``group_stats``, indexed DEG) consumed by ``static/expression_tab.js``.
"""

from __future__ import annotations

import gzip
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

from expression_local import expression_db_available, fetch_local_omics

DEFAULT_PLANTAPP_BASE = "https://www.plantapp.org"
REQUEST_TIMEOUT_SEC = 120.0
USER_AGENT = "PanViewer/1.0"


def plantapp_base_url() -> str:
    return (os.environ.get("PLANTAPP_BASE_URL") or DEFAULT_PLANTAPP_BASE).rstrip("/")


def _decode_http_body(raw: bytes, resp: Any) -> str:
    hdrs = getattr(resp, "headers", None)
    enc = (hdrs.get("Content-Encoding") or "").lower() if hdrs else ""
    if enc == "gzip":
        try:
            raw = gzip.decompress(raw)
        except OSError:
            pass
    return raw.decode("utf-8", errors="replace")


def _http_get_json(url: str) -> tuple[Any | None, str | None, int | None]:
    """
    Returns (parsed_json, error_message, http_status).
    parsed_json may be dict or list on success.
    """
    req = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SEC) as resp:
            status = getattr(resp, "status", None) or 200
            raw_bytes = resp.read()
            raw = _decode_http_body(raw_bytes, resp)
            try:
                return json.loads(raw), None, status
            except json.JSONDecodeError:
                return None, "Response was not valid JSON", status
    except urllib.error.HTTPError as e:
        try:
            raw_bytes = e.read()
            raw = _decode_http_body(raw_bytes, e)
            body = json.loads(raw)
            if isinstance(body, dict) and body.get("error"):
                return None, str(body["error"]), e.code
        except (json.JSONDecodeError, OSError, ValueError):
            pass
        return None, f"HTTP {e.code}", e.code
    except urllib.error.URLError as e:
        return None, f"Network error: {e.reason!s}", None
    except TimeoutError:
        return None, "Request timed out", None
    except OSError as e:
        return None, str(e), None


def fetch_plantapp_omics(gene_id: str, *, genome: str | None = None) -> dict[str, Any]:
    """
    Fetch tissue + DEG for one gene_id.

    If ``expression/expression.db`` exists, serve from it and never call PlantApp.
    Otherwise proxy PlantApp ``/api/gene-omics``.

    ``genome`` is passed through when set (e.g. ``HvMorex`` for barley on PlantApp).

    Returns a dict safe to jsonify:
      query_gene_id, ok, unknown_gene, error (optional), source (``local``|``plantapp``),
      tissue (optional), deg (optional), resolved_gene_id, resolved_genome
    """
    if expression_db_available():
        return fetch_local_omics(gene_id, genome=genome)

    gid = (gene_id or "").strip()
    out: dict[str, Any] = {
        "query_gene_id": gid,
        "ok": False,
        "unknown_gene": False,
        "tissue": None,
        "deg": None,
        "resolved_gene_id": None,
        "resolved_genome": None,
        "source": "plantapp",
    }
    if not gid:
        out["error"] = "gene_id is required"
        return out

    base = plantapp_base_url()
    omics_params: list[tuple[str, str]] = [
        ("gene_id", gid),
        ("tissue_format", "compact"),
        ("deg_format", "indexed"),
        ("group_stats", "1"),
    ]
    if genome:
        omics_params.append(("genome", genome.strip()))
    omics_url = f"{base}/api/gene-omics?{urllib.parse.urlencode(omics_params)}"

    body, err, status = _http_get_json(omics_url)

    if err:
        out["error"] = err
        if status == 429:
            out["rate_limited"] = True
        return out

    if not isinstance(body, dict):
        out["error"] = "Unexpected omics response"
        return out

    if body == {}:
        out["unknown_gene"] = True
        out["ok"] = True
        return out

    if body.get("error"):
        out["error"] = str(body["error"])
        return out

    out["resolved_gene_id"] = body.get("gene_id") or gid
    out["resolved_genome"] = body.get("genome")
    out["tissue"] = body.get("tissue")
    out["deg"] = body.get("deg")
    out["ok"] = True
    return out
