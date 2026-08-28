"""
Complete les metadonnees Site/Formation/Age/GC/SMT/Li/Loc/Obs (et lat/lon,
et en mode specimen, sample/site) d'un fichier .prmag DEJA converti,
depuis une table externe - demande explicite utilisateur, en 2 temps :
d'abord "we can build a routine to complete the sample information later
on" (pendant la discussion sur l'import "as-is" des fichiers legacy), puis
la specification concrete : "can we put a menu like complete sample
information with data in a table".

Deux modes, chacun avec ses propres colonnes de table (`#` optionnel en
tete de la ligne d'en-tete, comme les fichiers complement existants de
convert_legacy_ren.py) :

- Mode SITE (`complete_site_info`) : suppose que specimen/sample/site sont
  DEJA corrects dans le .prmag (typiquement le cas apres un import MagIC
  ou une conversion legacy qui a bien derive le site des 6 premiers
  caracteres du specimen) - la table ne fournit QUE les metadonnees
  complementaires, UNE LIGNE PAR SITE, appliquee a TOUS les specimens de
  ce site :
      #site	lat	lon	formation	age	geologic_classes	geologic_types	lithologies	location	obs

- Mode SPECIMEN (`complete_specimen_info`) : ne suppose RIEN de fiable sur
  sample/site dans le .prmag - la table fournit AUSSI sample/site
  (ecrases si presents), UNE LIGNE PAR SPECIMEN :
      #specimen	sample	site	lat	lon	formation	age	geologic_classes	geologic_types	lithologies	location	obs

Dans les deux modes, une colonne peut etre omise/laissee vide sans
consequence (seules les cles reellement fournies, non vides, sont
appliquees) - completer PROGRESSIVEMENT au fil de ce qui devient connu est
le but explicite de cette routine, pas un remplissage en un seul coup.

Patch le fichier .prmag EN PLACE au niveau du TEXTE BRUT des 4 lignes
d'entete par specimen (pas une reserialisation depuis les objets Pmag deja
parses par testlect.read_prmag_file) : seules les cles explicitement
fournies par la table sont remplacees, tout le reste (volume/mass/
elevation/comment, azimuth/dip/date, bed_dip, method_codes, et TOUTES les
lignes de mesure) reste OCTET PRES identique. Une reserialisation depuis
Pmag perdrait `method_codes` (jamais stocke dans le dataclass Pmag - voir
testlect.py) et risquerait de reformater des valeurs numeriques
differemment de l'original.

Une sauvegarde `.bak` du fichier .prmag est ecrite AVANT toute
modification, une seule fois (jamais ecrasee si elle existe deja - garde
le tout premier etat, avant la toute premiere completion)."""

import os
from typing import Dict, List, Optional, Tuple

_ROCHE_TABLE_FIELDS = (
    "formation", "age", "geologic_classes", "geologic_types", "lithologies", "location", "obs",
)


def _split_line(line: str) -> List[str]:
    return [f.strip() for f in line.rstrip("\n").split("\t")]


def _read_table_rows(path: str, encoding: str) -> Tuple[List[str], List[List[str]]]:
    with open(path, "r", encoding=encoding, errors="replace") as f:
        raw_lines = [l for l in f if l.strip()]
    if not raw_lines:
        return [], []
    header = [h.lstrip("#").strip().lower() for h in _split_line(raw_lines[0])]
    rows = [_split_line(l) for l in raw_lines[1:]]
    return header, rows


def _prmag_kv_line(line: str) -> dict:
    """Meme parsing que testlect._prmag_kv_line (non importe directement -
    module de bas niveau, pas de dependance croisee necessaire ici)."""
    result = {}
    for chunk in line.split("\t"):
        if ":" not in chunk:
            continue
        k, _sep, v = chunk.partition(":")
        result[k.strip()] = v.strip()
    return result


def _rewrite_kv_line(original_line: str, updates: Dict[str, str]) -> str:
    """Reconstruit une ligne 'cle: valeur\\tcle2: valeur2...' en ne
    changeant QUE les cles presentes dans `updates`, toutes les autres
    cles/valeurs de `original_line` restant identiques - y compris leur
    ORDRE et leur formatage d'origine (voir docstring module : jamais
    reconstruite depuis des champs Pmag deja parses)."""
    chunks = original_line.split("\t")
    new_chunks = []
    seen = set()
    for chunk in chunks:
        if ":" not in chunk:
            new_chunks.append(chunk)
            continue
        k, _sep, _v = chunk.partition(":")
        key = k.strip()
        if key in updates:
            new_chunks.append(f"{key}: {updates[key]}")
            seen.add(key)
        else:
            new_chunks.append(chunk)
    for key, value in updates.items():
        if key not in seen:
            new_chunks.append(f"{key}: {value}")
    return "\t".join(new_chunks)


def _apply_table(
    prmag_path: str, table: Dict[str, dict], key_field: str,
) -> Tuple[int, List[str]]:
    """Parcourt le .prmag bloc par specimen (meme detection que
    testlect.read_prmag_file), patche les lignes 'specimen:'/'formation:'
    de chaque bloc dont la cle (`site` ou `specimen` selon `key_field`)
    est trouvee dans `table`. Retourne (nb_specimens_mis_a_jour,
    liste_specimens_sans_correspondance)."""
    with open(prmag_path, "r", encoding="utf-8") as f:
        lines = [raw.rstrip("\n") for raw in f]
    original_text = "\n".join(lines) + "\n"
    n = len(lines)
    i = 0
    n_updated = 0
    unmatched: List[str] = []

    def skip_blank():
        # meme comportement que testlect.read_prmag_file._skip_blank_and_comments
        # (lignes vides ET commentaires "#...") - un premier essai qui ne
        # sautait que les lignes vides desalignait la lecture des le
        # premier bloc (les 5 lignes d'en-tete "#..." du fichier .prmag
        # etaient alors lues a tort comme un faux bloc specimen).
        nonlocal i
        while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#")):
            i += 1

    skip_blank()
    while i < n:
        if i + 4 >= n:
            break
        idx_a = i
        line_a = _prmag_kv_line(lines[i]); i += 1
        i += 1  # line_b - jamais touchee
        i += 1  # line_c - jamais touchee
        idx_d = i
        line_d = _prmag_kv_line(lines[i]); i += 1
        i += 1  # ligne d'en-tete des mesures
        while i < n and lines[i].strip() != "":
            i += 1

        specimen = line_a.get("specimen", "").strip()
        site = line_a.get("site", "").strip()
        lookup_key = specimen if key_field == "specimen" else site
        record = table.get(lookup_key)
        if record is None:
            unmatched.append(specimen or f"(line {idx_a + 1})")
            skip_blank()
            continue

        a_updates = {}
        if record.get("lat"):
            try:
                a_updates["lat"] = f"{float(record['lat']):.5f}"
            except ValueError:
                pass
        if record.get("lon"):
            try:
                a_updates["lon"] = f"{float(record['lon']):.5f}"
            except ValueError:
                pass
        if key_field == "specimen":
            if record.get("sample"):
                a_updates["sample"] = record["sample"]
            if record.get("site"):
                a_updates["site"] = record["site"]
        if a_updates:
            lines[idx_a] = _rewrite_kv_line(lines[idx_a], a_updates)

        d_updates = {field: record[field] for field in _ROCHE_TABLE_FIELDS if record.get(field)}
        if d_updates:
            lines[idx_d] = _rewrite_kv_line(lines[idx_d], d_updates)

        if a_updates or d_updates:
            n_updated += 1
        skip_blank()

    if n_updated:
        backup_path = prmag_path + ".bak"
        if not os.path.exists(backup_path):
            with open(backup_path, "w", encoding="utf-8") as f:
                f.write(original_text)
        with open(prmag_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    return n_updated, unmatched


def complete_site_info(prmag_path: str, table_path: str, encoding: str = "utf-8") -> Tuple[int, List[str]]:
    """Table indexee par SITE (colonne 'site' obligatoire) - voir
    docstring module. Retourne (nb_specimens_mis_a_jour,
    liste_specimens_sans_site_correspondant_dans_la_table)."""
    header, rows = _read_table_rows(table_path, encoding=encoding)
    if "site" not in header:
        raise ValueError("Table file must have a 'site' column header (site-level mode).")
    table: Dict[str, dict] = {}
    for row in rows:
        if len(row) < len(header):
            continue
        record = dict(zip(header, row))
        site = record.get("site", "").strip()
        if site:
            table[site] = record
    return _apply_table(prmag_path, table, key_field="site")


def complete_specimen_info(prmag_path: str, table_path: str, encoding: str = "utf-8") -> Tuple[int, List[str]]:
    """Table indexee par SPECIMEN complet (colonne 'specimen' obligatoire) -
    voir docstring module. Retourne (nb_specimens_mis_a_jour,
    liste_specimens_sans_correspondance_dans_la_table)."""
    header, rows = _read_table_rows(table_path, encoding=encoding)
    if "specimen" not in header:
        raise ValueError("Table file must have a 'specimen' column header (specimen-level mode).")
    table: Dict[str, dict] = {}
    for row in rows:
        if len(row) < len(header):
            continue
        record = dict(zip(header, row))
        specimen = record.get("specimen", "").strip()
        if specimen:
            table[specimen] = record
    return _apply_table(prmag_path, table, key_field="specimen")
