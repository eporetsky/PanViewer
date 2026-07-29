"""
Load ``database/config.json`` and resolve variant → SQLite paths for PanViewer.

When config.json is missing, each ``database/<stem>.db`` becomes its own species tab
(one variant, id = stem).

Config shape (no ``method`` / ``default_variant`` fields)::

    {
      "species": {
        "wheat": {
          "label": "Wheat",
          "variants": {
            "wheat": {"db": "wheat.db", "label": "Pandagma"}
          }
        }
      }
    }

The first variant listed under a species is the default when the species tab id is used.
"""
from __future__ import annotations

import json
import os
from functools import lru_cache
from typing import Any

DATABASE_CONFIG_NAME = "config.json"


def _database_dir() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "database")


def _config_path() -> str:
    return os.path.join(_database_dir(), DATABASE_CONFIG_NAME)


def _load_config_document(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    if not isinstance(doc, dict):
        raise RuntimeError(f"{path}: expected a JSON object at the top level")
    return doc


def _auto_species_from_db_dir(db_dir: str) -> dict[str, Any]:
    """One species tab per ``<stem>.db`` when config.json is absent."""
    species: dict[str, Any] = {}
    if not os.path.isdir(db_dir):
        return {"species": species}
    for fn in sorted(os.listdir(db_dir)):
        if not fn.endswith(".db"):
            continue
        stem = os.path.splitext(fn)[0].strip().lower()
        if not stem or stem.startswith("."):
            continue
        label = stem[:1].upper() + stem[1:] if stem else stem
        species[stem] = {
            "label": label,
            "variants": {stem: {"db": fn, "label": label}},
        }
    return {"species": species}


@lru_cache(maxsize=1)
def raw_database_config() -> dict[str, Any]:
    db_dir = _database_dir()
    json_path = _config_path()
    if os.path.isfile(json_path):
        return _load_config_document(json_path)
    return _auto_species_from_db_dir(db_dir)


def clear_database_config_cache() -> None:
    raw_database_config.cache_clear()


def species_configs() -> dict[str, dict[str, Any]]:
    """Species tab id → {label, variants: {variant_id: {...}}} (insertion order)."""
    raw = raw_database_config()
    block = raw.get("species")
    if not isinstance(block, dict):
        raise RuntimeError(
            f"{_config_path()}: missing or invalid top-level 'species' mapping"
        )
    out: dict[str, dict[str, Any]] = {}
    for sid, entry in block.items():
        species_id = (sid or "").strip().lower()
        if not species_id:
            continue
        if not isinstance(entry, dict):
            raise RuntimeError(f"{_config_path()}: species '{species_id}' must be a mapping")
        variants_in = entry.get("variants")
        if not isinstance(variants_in, dict) or not variants_in:
            raise RuntimeError(
                f"{_config_path()}: species '{species_id}' needs a non-empty 'variants' map"
            )
        variants: dict[str, dict[str, str]] = {}
        for vid, ventry in variants_in.items():
            variant_id = (vid or "").strip().lower()
            if not variant_id or not isinstance(ventry, dict):
                continue
            db_fn = (ventry.get("db") or f"{variant_id}.db").strip()
            if not db_fn.endswith(".db"):
                db_fn = f"{db_fn}.db"
            vlabel = (ventry.get("label") or "").strip()
            if not vlabel:
                vlabel = variant_id[:1].upper() + variant_id[1:]
            kw_ds = (ventry.get("keyword_dataset_id") or species_id).strip().lower()
            variants[variant_id] = {
                "db": db_fn,
                "label": vlabel,
                "keyword_dataset_id": kw_ds,
            }
        if not variants:
            raise RuntimeError(
                f"{_config_path()}: species '{species_id}' has no valid variants"
            )
        label = (entry.get("label") or "").strip()
        if not label:
            label = species_id[:1].upper() + species_id[1:]
        out[species_id] = {
            "label": label,
            "variants": variants,
        }
    return out


def dataset_configs() -> dict[str, dict[str, Any]]:
    """
    Variant id → runtime config (db_path, species_id, variant_label, …).

    Variants whose DB file is missing are omitted.
    """
    db_dir = _database_dir()
    cfg: dict[str, dict[str, Any]] = {}
    for species_id, spec in species_configs().items():
        for variant_id, v in spec["variants"].items():
            db_path = os.path.join(db_dir, v["db"])
            if not os.path.isfile(db_path):
                continue
            cfg[variant_id] = {
                "db_path": db_path,
                "species_id": species_id,
                "species_label": spec["label"],
                "variant_id": variant_id,
                "variant_label": v["label"],
                "label": spec["label"],
                "keyword_dataset_id": v["keyword_dataset_id"],
                "mode": "pandagma",
                "enable_cluster_picker": False,
                "enable_homeologue_panels": False,
            }
    return cfg


def species_id_for_dataset(dataset_id: str) -> str:
    did = (dataset_id or "").strip().lower()
    cfg = dataset_configs()
    if did in cfg:
        return cfg[did]["species_id"]
    return did


def keyword_dataset_id_for_dataset(dataset_id: str) -> str:
    did = (dataset_id or "").strip().lower()
    cfg = dataset_configs()
    if did in cfg:
        return cfg[did]["keyword_dataset_id"]
    return did


def default_variant_for_species(species_id: str) -> str | None:
    """First installed variant listed under the species in config order."""
    sid = (species_id or "").strip().lower()
    spec = species_configs().get(sid)
    if not spec:
        return None
    cfg = dataset_configs()
    for vid in spec["variants"]:
        if vid in cfg:
            return vid
    return None


def resolve_dataset_id(dataset_id: str) -> str | None:
    """
    Accept species tab id (→ first installed variant) or variant id.
    Returns a variant id with an on-disk DB, or None.
    """
    did = (dataset_id or "").strip().lower()
    cfg = dataset_configs()
    if did in cfg:
        return did
    return default_variant_for_species(did)


def variants_for_species(species_id: str) -> list[dict[str, str]]:
    """Installed variants for a species tab, in config order."""
    sid = (species_id or "").strip().lower()
    spec = species_configs().get(sid)
    if not spec:
        return []
    cfg = dataset_configs()
    out: list[dict[str, str]] = []
    for vid in spec["variants"]:
        if vid not in cfg:
            continue
        out.append(
            {
                "id": vid,
                "label": cfg[vid]["variant_label"],
            }
        )
    return out
