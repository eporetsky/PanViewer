"""
Server-side fetch of PlantApp gene omics JSON (same contract as plantapp/pages/api.py).

Endpoints (GET, query gene_id only; genome resolved on PlantApp):
  /api/gene-tissue-expression
  /api/differential-expression
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

DEFAULT_PLANTAPP_BASE = "https://www.plantapp.org"
REQUEST_TIMEOUT_SEC = 60.0
USER_AGENT = "PanViewer/1.0"


def plantapp_base_url() -> str:
    return (os.environ.get("PLANTAPP_BASE_URL") or DEFAULT_PLANTAPP_BASE).rstrip("/")


def _http_get_json(url: str) -> tuple[Any | None, str | None, int | None]:
    """
    Returns (parsed_json, error_message, http_status).
    parsed_json may be dict or list on success.
    """
    req = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT},
        method="GET",
    )
    try:
        with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SEC) as resp:
            status = getattr(resp, "status", None) or 200
            raw = resp.read().decode("utf-8", errors="replace")
            try:
                return json.loads(raw), None, status
            except json.JSONDecodeError:
                return None, "Response was not valid JSON", status
    except urllib.error.HTTPError as e:
        try:
            raw = e.read().decode("utf-8", errors="replace")
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


def fetch_plantapp_omics(gene_id: str) -> dict[str, Any]:
    """
    Fetch tissue + DEG for one gene_id from PlantApp.

    Returns a dict safe to jsonify:
      query_gene_id, ok, unknown_gene, error (optional),
      tissue (optional), deg (optional), resolved_gene_id, resolved_genome
    """
    gid = (gene_id or "").strip()
    out: dict[str, Any] = {
        "query_gene_id": gid,
        "ok": False,
        "unknown_gene": False,
        "tissue": None,
        "deg": None,
        "resolved_gene_id": None,
        "resolved_genome": None,
    }
    if not gid:
        out["error"] = "gene_id is required"
        return out

    base = plantapp_base_url()
    q = urllib.parse.urlencode({"gene_id": gid, "format": "compact"})
    tissue_url = f"{base}/api/gene-tissue-expression?{q}"
    deg_url = f"{base}/api/differential-expression?{urllib.parse.urlencode({'gene_id': gid})}"

    t_body, t_err, t_status = _http_get_json(tissue_url)
    if t_err:
        out["error"] = t_err
        if t_status == 429:
            out["rate_limited"] = True
        return out

    if not isinstance(t_body, dict):
        out["error"] = "Unexpected tissue response"
        return out

    if t_body == {}:
        out["unknown_gene"] = True
        out["ok"] = True
        return out

    if t_body.get("error"):
        out["error"] = str(t_body["error"])
        return out

    out["resolved_gene_id"] = t_body.get("gene_id") or gid
    out["resolved_genome"] = t_body.get("genome")
    out["tissue"] = t_body.get("tissue")

    d_body, d_err, d_status = _http_get_json(deg_url)
    if d_err:
        out["deg_error"] = d_err
        if d_status == 429:
            out["rate_limited"] = True
        out["ok"] = True
        return out

    if not isinstance(d_body, dict):
        out["deg_error"] = "Unexpected DEG response"
        out["ok"] = True
        return out

    if d_body == {}:
        # Inconsistent vs tissue; still show tissue if we have it.
        out["deg_error"] = "Empty DEG response (gene may be unknown to that endpoint)."
        out["ok"] = True
        return out

    if d_body.get("error"):
        out["deg_error"] = str(d_body["error"])
        out["ok"] = True
        return out

    out["deg"] = d_body.get("deg")
    out["resolved_gene_id"] = d_body.get("gene_id") or out["resolved_gene_id"]
    out["resolved_genome"] = d_body.get("genome") or out["resolved_genome"]
    out["ok"] = True
    return out
