"""
Port Python de Convert2newpmagformat/convert_oldpmag_2ren.f95 : mise a
niveau d'anciens fichiers de donnees (1 ou 2 lignes d'entete par
echantillon) vers le format actuel a 3 lignes (`Id:`/`L:`/`Site: "..." ...`),
en completant les metadonnees roche/site depuis un fichier "complement"
externe plutot que depuis le texte deja present dans le vieux fichier.

Trois cas, correspondant aux 3 sous-routines du Fortran d'origine :

- CAS 1 (`withsimplespecinfo`) : vieux fichier a 1 ligne (juste `Id:`,
  aucune ligne `L:`). Le complement est indexe par SPECIMEN complet (le nom
  de site n'est pas forcement un prefixe du nom de specimen) et fournit
  aussi site/sample/lat/long. Colonnes (separees par espaces/tabs) :
  `specimen site sample lat long age gc smt li loc obs` (PAS de colonne
  `fm` dans cette lecture Fortran - `read(...,err=997) specimen,site,
  sample, rlat,rlon,age,gc,smt,li,loc,obs`, verifie sur le code source :
  fm est absent, ecrit comme chaine vide en sortie).

- CAS 2 (`withsimplesiteinfo`) : vieux fichier a 1 ligne, mais le nom de
  site EST un prefixe du nom de specimen. Le complement est indexe par
  site (prefixe) et fournit aussi lat/long/annee. Colonnes :
  `annee site lat long fm age gc smt li loc obs`. Le nom de site final
  ecrit est `<2 chiffres annee>` + le site du complement (ex. annee 2017
  -> "17"+site), pas le site du complement seul.

- CAS 3 (`withnewformat`) : vieux fichier DEJA a 2 lignes (`Id:`+`L:`,
  lat/long/date deja presentes). Le complement est indexe par site
  (prefixe du specimen) et fournit uniquement fm/age/gc/smt/li/loc/obs
  (pas de lat/long, deja dans la ligne `L:` existante). Colonnes :
  `site fm age gc smt li loc obs`. VERIFIE OCTET PRES contre la paire
  reelle Convert2newpmagformat/old_pmag.txt -> old_pmag.ren (via
  complement.txt) : 2489 echantillons, sortie identique.

Seul le CAS 3 a pu etre verifie contre un exemple reel (paire
avant/apres + fichier complement fournis par l'utilisateur). Les CAS 1 et
2 sont des ports fideles du Fortran mais n'ont pas d'exemple reel
correspondant dans ce projet pour verification - a tester sur un vrai
fichier avant un usage en production.

Fichier complement : ligne d'en-tete OPTIONNELLE (voir _FIELD_ALIASES) -
si la premiere ligne du fichier complement correspond exactement a des
noms de colonnes reconnus (vocabulaire aligne sur MagIC, ex. "site",
"formation", "geologic_age"...), l'ordre REEL des colonnes est utilise
plutot que l'ordre positionnel fige ci-dessus - demande explicite de
l'utilisateur ("plus de flexibilite"). Sans en-tete reconnu, le
comportement (et l'ordre des colonnes) reste identique a l'origine, pour
la retro-compatibilite avec les fichiers complement existants.

Point d'entree recommande : `convert_legacy_auto` (fin de ce module),
demande explicite de l'utilisateur ("single import legacy files (able to
import 1 ligne, 2 lignes and 3 lines header) - the app will automatically
recognize if there is one, two or three lines") - detecte PAR SPECIMEN le
nombre de lignes d'entete deja presentes (evite de faire choisir un cas
au numero) et le mode d'indexation du complement (specimen complet ou
prefixe de site) a partir de son en-tete, au lieu du dispatcher manuel
convert_case1_specimen_info/convert_case2_site_info/convert_case3_new_format
ci-dessous (conserves pour un complement SANS en-tete reconnu, ou un
usage direct)."""

import re
from typing import Dict, List, Optional, TextIO, Tuple

from testlect import parse_measure_line

_ROCHE_FIELDS = ("fm", "age", "gc", "smt", "li", "loc", "obs")

# Vocabulaire cod1/cod2 reconnu ailleurs dans le code (selection.py,
# calcul.py, magic_export.py, datatools.py - releve exhaustif par grep sur
# tous les `cod1 ==`/`cod1 in (...)`) - demande explicite utilisateur
# ("you may encounter data with code1 and code2 that are not recognized...
# list the unrecognized code1 and code2... ask the user to replace the not
# recognized codes in the legacy files") : un vieux fichier saisi/edite a
# la main peut porter une faute de frappe (ex. cod1='O' pour '0', cod2
# minuscule errone) qui passerait silencieusement inapercue jusqu'a fausser
# un calcul bien plus tard (mauvaise classification thermique/AF/
# paleointensite) - mieux vaut la signaler des l'import.
# N/D/T/K/S = thermique (NRM/demag thermique, alias historiques) ; A/F =
# desaimantation AF ; I = acquisition IRM ; R/V/P = pas Thellier (chauffe/
# verification/controle pTRM) ; X/Y/Z = pas ATRM ; L/Q = refroidissement
# lent/rapide.
_RECOGNIZED_COD1 = set("NDTKSAFIRVPXYZLQ")
# Lettres A-Z (identifiant de groupe sequentiel pour R/V/P/X/Y/Z, voir
# selection.py/convert_ren_to_r.py), chiffres 0-9 (releve reel sur des
# fichiers .ren/.txt du projet : '0' tres frequent sur le premier point
# NRM, '1'-'5' aussi presents) + 'p' (point simple D/N/S/T/K) + 'i' (I
# minuscule, IRM) + '+'/'-'/'=' (signe ATRM/GRM) - voir magic_export.py
# et datatools.eliminate_grm. Verifie par grep sur 31 fichiers reels du
# dossier reference/ avant de figer cette liste (une premiere version
# sans les chiffres remontait des milliers de faux positifs sur cod2='0').
_RECOGNIZED_COD2 = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") | {"p", "i", "+", "-", "="}


def _check_measurement_codes(
    specimen: str, line: str, anomalies: List[Tuple[str, int, str, str]],
) -> None:
    """Si `line` se parse comme une ligne de mesure (voir
    testlect.parse_measure_line) dont cod1 et/ou cod2 ne sont pas dans le
    vocabulaire reconnu, ajoute (specimen, etape, cod1, cod2) a
    `anomalies` - liste brute, le regroupement/troncature pour l'affichage
    est a la charge de l'appelant (voir app.py, meme convention que la
    liste `unmatched`)."""
    meas = parse_measure_line(line)
    if meas is None:
        return
    if meas.cod1 not in _RECOGNIZED_COD1 or meas.cod2 not in _RECOGNIZED_COD2:
        anomalies.append((specimen, meas.etape, meas.cod1, meas.cod2))

# Alias reconnus pour une eventuelle ligne d'en-tete du fichier complement -
# vocabulaire aligne sur celui deja utilise pour l'import/export MagIC
# (extract_magic.py : geologic_formations/geologic_age/geologic_classes/
# geologic_types/lithology/locality). Permet au fichier complement de
# declarer ses colonnes DANS N'IMPORTE QUEL ORDRE (voire d'en omettre), au
# lieu de dependre d'un ordre positionnel fige - demande explicite de
# l'utilisateur ("plus de flexibilite"). Si la premiere ligne ne correspond
# a AUCUN de ces alias, elle est traitee comme une donnee (ordre par
# defaut, fidele au Fortran d'origine) - retro-compatible avec les
# fichiers complement existants qui n'ont pas cet en-tete.
_FIELD_ALIASES: Dict[str, set] = {
    "specimen": {"specimen"},
    "site": {"site"},
    "sample": {"sample"},
    "year": {"year", "iyear"},
    "lat": {"lat", "latitude"},
    "lon": {"lon", "long", "longitude"},
    "fm": {"fm", "formation", "geologic_formations", "geologic_formation"},
    "age": {"age", "geologic_age"},
    "agemin": {"agemin", "age_low", "age_min"},
    "agemax": {"agemax", "age_high", "age_max"},
    "ageunit": {"ageunit", "age_unit"},
    "gc": {"gc", "geologic_class", "geologic_classes"},
    "smt": {"smt", "geologic_type", "geologic_types"},
    "li": {"li", "lithology", "lithologies"},
    "loc": {"loc", "locality"},
    "obs": {"obs", "description"},
}


def _split_list_directed(line: str) -> List[str]:
    """Equivalent (best-effort) d'un `read(chaine,*)` Fortran list-directed
    sur une ligne separee par tabulations avec champs textuels entre
    guillemets : coupe sur les tabulations, degage les guillemets de
    chaque champ."""
    fields = [f.strip() for f in line.rstrip("\n").split("\t")]
    return [f[1:-1] if len(f) >= 2 and f[0] == '"' and f[-1] == '"' else f for f in fields]


def _detect_header(fields: List[str], allowed_keys: set) -> Optional[List[str]]:
    """Si CHAQUE champ de `fields` correspond (insensible a la casse) a un
    alias reconnu parmi `allowed_keys` (sans doublon), retourne la liste
    des cles internes dans l'ordre reel des colonnes - sinon None (pas un
    en-tete reconnu, a traiter comme une ligne de donnees)."""
    resolved = []
    for raw in fields:
        key = raw.strip().lower()
        match = next((k for k in allowed_keys if key in _FIELD_ALIASES[k]), None)
        if match is None:
            return None
        resolved.append(match)
    if len(set(resolved)) != len(resolved):
        return None
    return resolved


def _read_complement_lines(path: str, encoding: str) -> List[str]:
    with open(path, "r", encoding=encoding) as f:
        return [line for line in f if line.strip() and not line.startswith("!")]


# Unites d'age reconnues par `checkage` (fichiers_magic.f:2225) - seules
# celles-ci sont recherchees en sous-chaine par le Fortran pour deduire
# `age_unit` a partir du texte libre ; toute autre valeur y reste ignoree
# silencieusement. Ici, validee explicitement plutot que silencieusement
# ignoree - demande explicite de l'utilisateur ("voir code Magic").
_VALID_AGE_UNITS = {
    "Ga", "Ma", "Years AD (+/-)", "ka", "Years BP",
    "Years Cal AD (+/-)", "Years Cal BP",
}


def _validate_age_header(header: List[str]) -> None:
    """Si l'en-tete du fichier complement declare une info d'age, exige
    EXACTEMENT l'une des deux combinaisons du schema MagIC (age/age_low/
    age_high/age_unit - voir l'export site.txt de fichiers_magic.f et
    `checkage`) : `agemin`+`agemax`+`ageunit` (plage) OU `age`+`ageunit`
    (valeur unique) - jamais un `age` sans unite, ni un melange des deux
    formes. Demande explicite de l'utilisateur ("obliger les
    utilisateurs")."""
    present = {k for k in header if k in ("age", "agemin", "agemax", "ageunit")}
    if not present or present == {"age", "ageunit"} or present == {"agemin", "agemax", "ageunit"}:
        return
    raise ValueError(
        "Complement file header: age info must be either "
        "'agemin'+'agemax'+'ageunit' or 'age'+'ageunit' (never 'age' alone, "
        "nor a mix of both forms) - found: " + ", ".join(sorted(present))
    )


def _compose_age(record: Dict[str, str]) -> str:
    """Construit la valeur texte du champ `Age:` (meme convention que les
    fichiers complement existants, ex. '-5000 - -10000 Years Cal AD (+/-)'
    ou '861 ± 122 Years Cal AD (+/-)') a partir des colonnes structurees
    age/agemin/agemax/ageunit declarees en en-tete. `ageunit` est valide
    contre le vocabulaire reconnu par `checkage`."""
    ageunit = record.get("ageunit", "").strip()
    if ageunit and ageunit not in _VALID_AGE_UNITS:
        raise ValueError(
            f"Complement file: unrecognized ageunit {ageunit!r} - must be one of "
            + ", ".join(sorted(_VALID_AGE_UNITS))
        )
    if "agemin" in record or "agemax" in record:
        agemin = record.get("agemin", "").strip()
        agemax = record.get("agemax", "").strip()
        return f"{agemin} - {agemax} {ageunit}".strip()
    if "age" in record:
        age = record.get("age", "").strip()
        return f"{age} {ageunit}".strip()
    return record.get("age", "")


_NOT_SPECIFIED = "Not Specified"


def _write_roche_line(out: TextIO, site: str, sample: str, values: Dict[str, str]) -> None:
    """Ecrit la ligne `Site: "..." Sample: "..." Fm: "..." ...` (format
    101 du Fortran, identique a celui deja lu par testlect.decode_roche).
    Chaque valeur passe par un equivalent de `trim()` Fortran - espaces de
    FIN seulement retires, les espaces de DEBUT (significatifs dans
    complement.txt, ex. ' 20 - 200 ka ') sont preserves - verifie octet
    pres contre Convert2newpmagformat/old_pmag.ren.

    Un champ manquant/vide est ecrit "Not Specified" plutot que "" -
    convention MagIC pour une info absente (voir fichiers_magic.f:526,531,
    `site.lithology="Not Specified"` / `site.type1="Not Specified"`,
    generalisee ici aux 7 champs roche) - demande explicite de
    l'utilisateur. Sans effet sur les fichiers complement existants (aucun
    champ vide dans complement.txt - verifie)."""
    parts = [f'Site: "{site.rstrip()}"', f'Sample: "{sample.rstrip()}"']
    parts += [f'{key.upper() if key in ("gc", "smt") else key.capitalize()}: "{values.get(key, "").rstrip() or _NOT_SPECIFIED}"'
              for key in _ROCHE_FIELDS]
    out.write(" ".join(parts) + "\n")


# ---------------------------------------------------------------------------
# CAS 3 : vieux fichier deja a 2 lignes (Id:+L:), complement indexe par site
# (prefixe du specimen), colonnes site/fm/age/gc/smt/li/loc/obs.
# ---------------------------------------------------------------------------

_CASE3_DEFAULT_ORDER = ["site", "fm", "age", "gc", "smt", "li", "loc", "obs"]
# "age" y est aussi remplacable par agemin+agemax+ageunit dans l'en-tete
# (voir _validate_age_header) - pas dans l'ordre par defaut (fidele au
# Fortran, une seule colonne "age" en texte libre sans en-tete).
_CASE3_ALLOWED_KEYS = set(_CASE3_DEFAULT_ORDER) | {"agemin", "agemax", "ageunit"}


def _load_complement_case3(path: Optional[str], encoding: str = "latin-1") -> List[Tuple[str, Dict[str, str]]]:
    rows: List[Tuple[str, Dict[str, str]]] = []
    if not path:
        return rows
    lines = _read_complement_lines(path, encoding)
    if not lines:
        return rows

    header = _detect_header(_split_list_directed(lines[0]), _CASE3_ALLOWED_KEYS)
    if header is not None:
        _validate_age_header(header)
    order = header if header is not None else _CASE3_DEFAULT_ORDER
    data_lines = lines[1:] if header is not None else lines

    for line in data_lines:
        fields = _split_list_directed(line)
        if len(fields) < len(order):
            continue
        record = dict(zip(order, fields))
        if header is not None:
            record["age"] = _compose_age(record)
        site = record.pop("site", "")
        rows.append((site, record))
    return rows


def convert_case3_new_format(
    old_path: str, complement_path: Optional[str], output_path: str,
    encoding: str = "latin-1", complement_encoding: str = "utf-8",
) -> Tuple[int, List[str]]:
    """Port de `withnewformat` : vieux fichier a 2 lignes -> 3 lignes.
    `complement_encoding` distinct de `encoding` : le fichier complement
    (edite a la main, ex. caracteres "±") est souvent en UTF-8 alors que
    les fichiers de donnees .txt/.ren restent en latin-1 - verifie sur
    complement.txt (caractere "±" present, octets UTF-8 confirmes).

    `complement_path` est OPTIONNEL (None ou "") - demande explicite
    utilisateur ("import these files as they are and later we can
    complete the missing information... most of the tasks within the
    program just need the core correction, bedding and volume or mass...
    most of the missing information in the oldest file are for an export
    to Magic") : la ligne `Id:` du vieux fichier (deja recopiee telle
    quelle, jamais touchee par le complement) porte DEJA in:/az:/dip:/
    str:/v: (correction de carotte, pendage, volume/masse - voir
    testlect.parse_id_line), tout ce dont la plupart des taches du
    programme ont besoin. Le complement ne fournit QUE Site/Fm/Age/GC/
    SMT/Li/Loc/Obs, uniquement utiles pour un export MagIC. Un specimen
    sans complement (ou sans complement du tout) recoit donc quand meme
    une ligne roche complete, avec site="Not Specified" et les 7 champs
    roche a "Not Specified" (meme convention que _write_roche_line pour
    un champ individuel manquant), plutot que d'etre saute (ancien
    comportement : la ligne roche etait purement omise, le specimen etant
    seulement liste dans `unmatched`) - une metadonnee MagIC absente peut
    etre completee plus tard (routine a venir), elle ne doit jamais faire
    perdre une ligne de mesure ou bloquer toute la conversion.
    Retourne (nb_echantillons_convertis, liste_specimens_sans_complement) -
    ce 2e element reste purement informatif (aucun n'est plus jamais
    "perdu")."""
    complement = _load_complement_case3(complement_path, encoding=complement_encoding)

    nb_converted = 0
    unmatched: List[str] = []
    specimen = ""
    sample = ""

    with open(old_path, "r", encoding=encoding) as fin, \
         open(output_path, "w", encoding="utf-8") as fout:
        for raw_line in fin:
            line = raw_line.rstrip("\n")

            if line[:3] == "Id:":
                fout.write(line.rstrip() + "\n")
                specimen = line[3:15].strip()
                sample = specimen[:-1] if len(specimen) > 1 else specimen
                continue

            if line[:2] == "L:":
                fout.write(line[:90].rstrip() + "\n")
                match = next((row for row in complement if specimen.startswith(row[0])), None)
                if match is None:
                    unmatched.append(specimen)
                    site, values = "Not Specified", {}
                else:
                    site, values = match
                _write_roche_line(fout, site, sample, values)
                nb_converted += 1
                continue

            fout.write(line.rstrip() + "\n")

    return nb_converted, unmatched


# ---------------------------------------------------------------------------
# CAS 2 : vieux fichier a 1 ligne, site = prefixe du specimen, complement
# indexe par site et fournit aussi annee/lat/long.
# ---------------------------------------------------------------------------

_CASE2_DEFAULT_ORDER = ["year", "site", "lat", "lon", "fm", "age", "gc", "smt", "li", "loc", "obs"]
_CASE2_ALLOWED_KEYS = set(_CASE2_DEFAULT_ORDER) | {"agemin", "agemax", "ageunit"}


def _load_complement_case2(path: Optional[str], encoding: str = "latin-1") -> List[Tuple[str, int, float, float, Dict[str, str]]]:
    rows = []
    if not path:
        return rows
    lines = _read_complement_lines(path, encoding)
    if not lines:
        return rows

    header = _detect_header(_split_list_directed(lines[0]), _CASE2_ALLOWED_KEYS)
    if header is not None:
        _validate_age_header(header)
    order = header if header is not None else _CASE2_DEFAULT_ORDER
    data_lines = lines[1:] if header is not None else lines

    for line in data_lines:
        fields = _split_list_directed(line)
        if len(fields) < len(order):
            continue
        record = dict(zip(order, fields))
        if header is not None:
            record["age"] = _compose_age(record)
        try:
            iyear = int(float(record.get("year", 0)))
            rlat = float(record.get("lat", 0.0))
            rlon = float(record.get("lon", 0.0))
        except ValueError:
            continue
        site1 = record.get("site", "")
        values = {k: record.get(k, "") for k in _ROCHE_FIELDS}
        rows.append((site1, iyear, rlat, rlon, values))
    return rows


def convert_case2_site_info(
    old_path: str, complement_path: Optional[str], output_path: str,
    encoding: str = "latin-1", complement_encoding: str = "utf-8",
) -> Tuple[int, List[str]]:
    """Port de `withsimplesiteinfo` : vieux fichier a 1 ligne (site =
    prefixe du specimen) -> 3 lignes, `L:`+roche generees depuis le
    complement (annee/lat/long/fm/age/gc/smt/li/loc/obs). Voir
    convert_case3_new_format pour la raison de `complement_encoding`.

    `complement_path` est OPTIONNEL (meme demande utilisateur que
    convert_case3_new_format, voir sa docstring) : la ligne `Id:` du vieux
    fichier (deja recopiee, jamais touchee par le complement) porte deja
    in:/az:/dip:/str:/v: (correction de carotte, pendage, volume/masse).
    Le complement ne fournit QUE l'annee/lat/lon et Site/Fm/Age/GC/SMT/Li/
    Loc/Obs (utiles pour la datation et l'export MagIC, pas pour les
    taches courantes). Un specimen sans complement (ou sans complement du
    tout) recoit quand meme une ligne `L:` (annee/lat/lon a 0, faute de
    mieux) et une ligne roche complete (site="Not Specified", les 7
    champs roche a "Not Specified") plutot que d'etre saute entierement -
    ancien comportement, qui perdait meme le volume/la correction de
    carotte deja presents sur la ligne Id: en cas d'absence de
    complement."""
    complement = _load_complement_case2(complement_path, encoding=complement_encoding)

    nb_converted = 0
    unmatched: List[str] = []
    specimen = ""
    sample = ""

    with open(old_path, "r", encoding=encoding) as fin, \
         open(output_path, "w", encoding="utf-8") as fout:
        for raw_line in fin:
            line = raw_line.rstrip("\n")

            if line[:3] == "Id:":
                fout.write(line[:82].rstrip() + "\n")
                specimen = line[3:15].strip()
                sample = specimen[:-1] if len(specimen) > 1 else specimen

                match = next((row for row in complement if specimen.startswith(row[0])), None)
                if match is None:
                    unmatched.append(specimen)
                    site1, iyear, rlat, rlon, values = "Not Specified", 0, 0.0, 0.0, {}
                else:
                    site1, iyear, rlat, rlon, values = match
                yy = (iyear - 2000) if 1999 < iyear < 2010 else (iyear - 1900)
                site = f"{yy:02d}{site1}" if match is not None else site1

                fout.write(
                    f"L:{rlat:10.5f} G:{rlon:10.5f}  H:   0  T:{iyear:4d} "
                    f" 1  1  0  0   azm:  0.0 azs:  0.0  Or:X12_0_3_90\n"
                )
                _write_roche_line(fout, site1, sample, values)
                nb_converted += 1
                continue

            fout.write(line.rstrip() + "\n")

    return nb_converted, unmatched


# ---------------------------------------------------------------------------
# CAS 1 : vieux fichier a 1 ligne, site PAS forcement prefixe du specimen,
# complement indexe par specimen COMPLET, fournit aussi site/sample/lat/long
# (pas de `fm` - absent de la lecture Fortran d'origine).
# ---------------------------------------------------------------------------

_CASE1_DEFAULT_ORDER = ["specimen", "site", "sample", "lat", "lon", "age", "gc", "smt", "li", "loc", "obs"]
# "fm" absent de l'ordre par defaut (fidele a la lecture Fortran d'origine,
# voir docstring module) mais accepte quand meme comme alias reconnu dans
# une eventuelle ligne d'en-tete - un collegue peut choisir de l'ajouter.
_CASE1_ALLOWED_KEYS = set(_CASE1_DEFAULT_ORDER) | {"fm", "agemin", "agemax", "ageunit"}


def _load_complement_case1(path: Optional[str], encoding: str = "latin-1") -> Dict[str, Tuple[str, str, float, float, Dict[str, str]]]:
    rows: Dict[str, Tuple[str, str, float, float, Dict[str, str]]] = {}
    if not path:
        return rows
    lines = _read_complement_lines(path, encoding)
    if not lines:
        return rows

    header = _detect_header(_split_list_directed(lines[0]), _CASE1_ALLOWED_KEYS)
    if header is not None:
        _validate_age_header(header)
    order = header if header is not None else _CASE1_DEFAULT_ORDER
    data_lines = lines[1:] if header is not None else lines

    for line in data_lines:
        fields = _split_list_directed(line)
        if len(fields) < len(order):
            continue
        record = dict(zip(order, fields))
        if header is not None:
            record["age"] = _compose_age(record)
        specimen_key = record.get("specimen", "")
        if not specimen_key:
            continue
        try:
            rlat = float(record.get("lat", 0.0))
            rlon = float(record.get("lon", 0.0))
        except ValueError:
            continue
        values = {k: record.get(k, "") for k in _ROCHE_FIELDS}
        rows[specimen_key] = (record.get("site", ""), record.get("sample", ""), rlat, rlon, values)
    return rows


def convert_case1_specimen_info(
    old_path: str, complement_path: Optional[str], output_path: str,
    encoding: str = "latin-1", complement_encoding: str = "utf-8",
) -> Tuple[int, List[str]]:
    """Port de `withsimplespecinfo` : vieux fichier a 1 ligne, complement
    indexe par specimen COMPLET (site pas necessairement un prefixe - ce
    cas existe justement pour les collegues dont la convention de nommage
    n'encode pas le site en prefixe du specimen). Voir
    convert_case3_new_format pour la raison de `complement_encoding`.

    `complement_path` est OPTIONNEL (meme demande utilisateur que
    convert_case3_new_format, voir sa docstring) : la ligne `Id:` du vieux
    fichier (deja recopiee, jamais touchee par le complement) porte deja
    in:/az:/dip:/str:/v:. Le complement ne fournit QUE site/sample/lat/
    lon/Fm/Age/GC/SMT/Li/Loc/Obs. Un specimen sans complement (ou sans
    complement du tout) recoit quand meme une ligne `L:` (lat/lon a 0,
    faute de mieux) et une ligne roche complete (site="Not Specified",
    sample derive du specimen comme dans les cas 2/3, les 7 champs roche
    a "Not Specified") plutot que d'etre saute entierement."""
    complement = _load_complement_case1(complement_path, encoding=complement_encoding)

    nb_converted = 0
    unmatched: List[str] = []

    with open(old_path, "r", encoding=encoding) as fin, \
         open(output_path, "w", encoding="utf-8") as fout:
        for raw_line in fin:
            line = raw_line.rstrip("\n")

            if line[:3] == "Id:":
                fout.write(line[:82].rstrip() + "\n")
                specimen = line[3:15].strip()

                match = complement.get(specimen)
                if match is None:
                    unmatched.append(specimen)
                    fallback_sample = specimen[:-1] if len(specimen) > 1 else specimen
                    site, sample, rlat, rlon, values = "Not Specified", fallback_sample, 0.0, 0.0, {}
                else:
                    site, sample, rlat, rlon, values = match

                fout.write(
                    f"L:{rlat:10.5f} G:{rlon:10.5f}  H:   0  T:   0  0  0  0  0"
                    f"   azm:  0.0 azs:  0.0  Or:X12_0_3_90\n"
                )
                _write_roche_line(fout, site, sample, values)
                nb_converted += 1
                continue

            fout.write(line.rstrip() + "\n")

    return nb_converted, unmatched


# ---------------------------------------------------------------------------
# Auto-detection : un seul point d'entree, detecte PAR SPECIMEN (pas par
# fichier entier - un meme fichier peut melanger d'anciens et de recents
# specimens) le nombre de lignes d'entete deja presentes, et - pour le
# complement - le mode d'indexation (specimen complet ou prefixe de site) a
# partir de son en-tete, au lieu de faire choisir un cas au numero -
# demande explicite utilisateur ("single import legacy files (able to
# import 1 ligne, 2 lignes and 3 lines header) - the app will
# automatically recognize if there is one, two or three lines").
#
# Un complement SANS en-tete reconnu (vieille convention positionnelle
# case 1/2/3, ambigue entre case 1 et case 2 - toutes deux a 11 colonnes
# par defaut, voir _CASE1_DEFAULT_ORDER/_CASE2_DEFAULT_ORDER) n'est PAS
# supporte ici : utiliser directement convert_case1_specimen_info/
# convert_case2_site_info/convert_case3_new_format pour un tel fichier.
# ---------------------------------------------------------------------------

_AUTO_ALLOWED_KEYS = set(_FIELD_ALIASES.keys())


def _load_complement_auto(path: Optional[str], encoding: str = "latin-1"):
    """Retourne (mode, table). `mode` est "specimen" (complement indexe
    par specimen complet - son en-tete declare une colonne "specimen"),
    "prefix" (indexe par prefixe de site - en-tete "site" sans
    "specimen"), ou None (pas de complement, fichier vide, ou en-tete non
    reconnu/insuffisant - traite comme "pas de complement" avec un
    avertissement, plutot que d'echouer toute la conversion pour une
    ligne d'en-tete mal formee). `table` est un Dict[str, dict] pour le
    mode "specimen", une List[Tuple[str, dict]] (ordre du fichier, pour
    un matching par prefixe) pour le mode "prefix"."""
    if not path:
        return None, None
    lines = _read_complement_lines(path, encoding)
    if not lines:
        return None, None

    header = _detect_header(_split_list_directed(lines[0]), _AUTO_ALLOWED_KEYS)
    if header is None or ("specimen" not in header and "site" not in header):
        print(
            "⚠️  Complement file has no recognized header declaring "
            "'specimen' or 'site' - automatic case detection needs one "
            "(e.g. a first line 'specimen\\tsite\\tsample\\tlat\\tlon\\t...' "
            "or 'site\\tformation\\tage\\t...'). Proceeding WITHOUT this "
            "complement (all specimens will get \"Not Specified\" geology "
            "fields) - use convert_case1_specimen_info/"
            "convert_case2_site_info/convert_case3_new_format directly for "
            "an older, headerless complement file."
        )
        return None, None
    _validate_age_header(header)

    records = []
    for line in lines[1:]:
        fields = _split_list_directed(line)
        if len(fields) < len(header):
            continue
        record = dict(zip(header, fields))
        record["age"] = _compose_age(record)
        records.append(record)

    if "specimen" in header:
        by_specimen = {r["specimen"]: r for r in records if r.get("specimen")}
        return "specimen", by_specimen
    by_prefix = [(r.get("site", ""), r) for r in records]
    return "prefix", by_prefix


def _lookup_auto(mode: Optional[str], table, specimen: str):
    """Cherche `specimen` dans `table` selon `mode` (voir
    _load_complement_auto). Retourne (found, site, sample_override,
    rlat, rlon, iyear, values) - `sample_override` est None sauf si le
    complement fournit sa propre colonne "sample" (mode "specimen"), a
    defaut le nom d'echantillon est derive du specimen comme dans les cas
    1/2/3 (voir l'appelant).

    Le site retourne (sauf mode "specimen" avec sa propre colonne "site",
    convention non-standard - voir CAS 1) vient des 6 premiers caracteres
    du specimen (annee sur 2 chiffres + code de site, convention Rennes) -
    demande explicite utilisateur ("consider the site number in the six
    first characters of the specimen numbers"), verifie contre
    old_pmag.ren (specimen "12CL10001A" -> site "12CL10", alors que
    complement.txt declarait "12CL100" pour la ligne qui matche par
    prefixe) : la colonne "site" du complement en mode "prefix" ne sert
    donc qu'a la RECHERCHE (quelle ligne de geologie s'applique), jamais a
    la valeur ECRITE."""
    site_from_specimen = specimen[:6] if specimen else "Not Specified"

    if mode == "specimen":
        record = table.get(specimen)
        if record is None:
            return False, site_from_specimen, None, 0.0, 0.0, 0, {}
        rlat = float(record.get("lat", "0") or 0.0)
        rlon = float(record.get("lon", "0") or 0.0)
        site = record.get("site", "").strip() or site_from_specimen
        sample_override = record.get("sample", "").strip() or None
        return True, site, sample_override, rlat, rlon, 0, _roche_values_from_record(record)

    if mode == "prefix":
        match = next(((site1, r) for site1, r in table if specimen.startswith(site1)), None)
        if match is None:
            return False, site_from_specimen, None, 0.0, 0.0, 0, {}
        _site1, record = match
        rlat = float(record.get("lat", "0") or 0.0)
        rlon = float(record.get("lon", "0") or 0.0)
        try:
            iyear = int(float(record.get("year", "0") or 0))
        except ValueError:
            iyear = 0
        return True, site_from_specimen, None, rlat, rlon, iyear, _roche_values_from_record(record)

    return False, site_from_specimen, None, 0.0, 0.0, 0, {}


def _roche_values_from_record(record: Dict[str, str]) -> Dict[str, str]:
    return {k: record.get(k, "") for k in _ROCHE_FIELDS}


def _has_multi_char_suffix(specimen: str) -> bool:
    """Vrai si les 2 DERNIERS caracteres du specimen sont tous deux des
    lettres (ex. "10CL0601AP") - demande explicite utilisateur ("do not
    import samples with multi-char suffices, just provide the list of
    non imported samples") : la convention standard (specimen = sample +
    UNE lettre, ex. "10CL0601A" -> sample "10CL0601") ne tient plus des
    que le suffixe fait 2 caracteres ou plus - `specimen[:-1]` donnerait
    alors un faux "sample" ("10CL0601A" au lieu de "10CL0601" pour
    "10CL0601AP"). Verifie sur donnees reelles (old_pmag.txt,
    Tibet_14_15_Pmag.txt) : 229 specimens sur 4839 suivent ce motif,
    TOUJOURS un suffixe "xP" (specimen pilote/paleointensite, x etant la
    lettre normale du specimen) - jamais une autre forme."""
    return len(specimen) >= 2 and specimen[-1].isalpha() and specimen[-2].isalpha()


def convert_legacy_auto(
    old_path: str, complement_path: Optional[str], output_path: str,
    encoding: str = "latin-1", complement_encoding: str = "utf-8",
) -> Tuple[int, List[str], List[Tuple[str, int, str, str]], List[str]]:
    """Point d'entree unique remplacant le choix manuel Case 1/2/3 -
    voir le commentaire de section ci-dessus. Detecte, PAR SPECIMEN :

    - 'Id:' seul, suivi directement d'une mesure -> 1 ligne (ancien case
      1/2) : `L:` et la ligne roche sont ENTIEREMENT FABRIQUEES.
    - 'Id:'+'L:', suivi directement d'une mesure -> 2 lignes (ancien
      case 3) : seule la ligne roche est FABRIQUEE, Id:/L: sont copiees
      telles quelles (elles portent deja in:/az:/dip:/str:/v: et lat/lon/
      date - jamais touchees par le complement).
    - 'Id:'+'L:'+'Site: ...' -> 3 lignes, deja complet : les 3 lignes
      sont copiees telles quelles, rien a fabriquer.

    Un specimen sans correspondance dans le complement (ou sans
    complement du tout, voir _load_complement_auto) recoit "Not
    Specified" pour le site et les 7 champs roche, lat/lon a 0 - jamais
    saute (meme philosophie que convert_case{1,2,3}_*, demande explicite
    utilisateur : "import these files as they are and later we can
    complete the missing information").

    Chaque ligne de mesure (copiee telle quelle, jamais reinterpretee) est
    aussi verifiee contre _RECOGNIZED_COD1/_RECOGNIZED_COD2 - demande
    explicite utilisateur ("you may encounter data with code1 and code2
    that are not recognized... list the unrecognized code1 and code2...
    ask the user to replace the not recognized codes in the legacy
    files") : un cod1/cod2 hors vocabulaire n'est PAS corrige ici (juste
    copie tel quel, la conversion ne doit pas bloquer dessus), seulement
    remonte dans le 3e element retourne pour que l'appelant (app.py)
    invite l'utilisateur a corriger la faute dans le fichier source.

    Un specimen dont le suffixe fait 2+ caracteres (voir
    _has_multi_char_suffix, ex. "10CL0601AP") N'EST PAS importe (ni ses
    lignes d'entete, ni ses mesures) - demande explicite utilisateur ("do
    not import samples with multi-char suffices, just provide the list of
    non imported samples") : la convention specimen=sample+1 lettre ne
    tient plus, "sample" ne peut pas etre derive de facon fiable pour ce
    cas, mieux vaut l'exclure explicitement (et le lister) que d'ecrire
    une metadonnee fausse.

    Retourne (nb_echantillons_convertis, liste_specimens_sans_complement,
    liste_(specimen,etape,cod1,cod2)_non_reconnus,
    liste_specimens_non_importes_suffixe_multi_caracteres) - tous les
    elements a partir du 2e sont purement informatifs, sauf le 4e qui
    correspond a des specimens reellement EXCLUS du fichier de sortie."""
    mode, table = _load_complement_auto(complement_path, encoding=complement_encoding)

    with open(old_path, "r", encoding=encoding) as fin:
        lines = fin.readlines()

    nb_converted = 0
    unmatched: List[str] = []
    code_anomalies: List[Tuple[str, int, str, str]] = []
    skipped_multi_suffix: List[str] = []
    n = len(lines)
    i = 0
    specimen = ""
    skip_current = False

    with open(output_path, "w", encoding="utf-8") as fout:
        while i < n:
            line = lines[i].rstrip("\n")

            if line[:3] != "Id:":
                if skip_current:
                    i += 1
                    continue
                _check_measurement_codes(specimen, line, code_anomalies)
                fout.write(line.rstrip() + "\n")
                i += 1
                continue

            specimen = line[3:15].strip()
            skip_current = _has_multi_char_suffix(specimen)
            if skip_current:
                skipped_multi_suffix.append(specimen)
                # consomme (sans ecrire) les lignes d'entete de CE
                # specimen, quel que soit leur nombre (1/2/3) - la meme
                # logique de detection que ci-dessous, juste sans sortie.
                next_line = lines[i + 1].rstrip("\n") if i + 1 < n else ""
                if next_line[:2] != "L:":
                    i += 1
                else:
                    i += 2
                    third_line = lines[i].rstrip("\n") if i < n else ""
                    if third_line.lstrip().lower().startswith("site:"):
                        i += 1
                continue

            sample = specimen[:-1] if len(specimen) > 1 else specimen

            next_line = lines[i + 1].rstrip("\n") if i + 1 < n else ""
            has_L = next_line[:2] == "L:"

            if not has_L:
                # 1 ligne : L: et roche entierement fabriquees (ancien case 1/2).
                fout.write(line[:82].rstrip() + "\n")
                found, site, sample_override, rlat, rlon, iyear, values = _lookup_auto(
                    mode, table, specimen)
                if not found:
                    unmatched.append(specimen)
                if iyear:
                    fout.write(
                        f"L:{rlat:10.5f} G:{rlon:10.5f}  H:   0  T:{iyear:4d}"
                        f"  1  1  0  0   azm:  0.0 azs:  0.0  Or:X12_0_3_90\n"
                    )
                else:
                    fout.write(
                        f"L:{rlat:10.5f} G:{rlon:10.5f}  H:   0  T:   0  0  0  0  0"
                        f"   azm:  0.0 azs:  0.0  Or:X12_0_3_90\n"
                    )
                _write_roche_line(fout, site, sample_override or sample, values)
                nb_converted += 1
                i += 1
                continue

            # deja au moins 2 lignes : Id:+L: copiees telles quelles (jamais
            # touchees par le complement - voir docstring).
            fout.write(line.rstrip() + "\n")
            fout.write(next_line[:90].rstrip() + "\n")
            i += 2

            third_line = lines[i].rstrip("\n") if i < n else ""
            if third_line.lstrip().lower().startswith("site:"):
                # deja 3 lignes : rien a fabriquer.
                fout.write(third_line.rstrip() + "\n")
                i += 1
                nb_converted += 1
                continue

            # 2 lignes : seule la roche est fabriquee (ancien case 3). `i`
            # n'avance PAS sur third_line (une mesure) - elle sera ecrite
            # normalement au prochain tour de boucle.
            found, site, sample_override, _rlat, _rlon, _iyear, values = _lookup_auto(
                mode, table, specimen)
            if not found:
                unmatched.append(specimen)
            _write_roche_line(fout, site, sample_override or sample, values)
            nb_converted += 1

    return nb_converted, unmatched, code_anomalies, skipped_multi_suffix
