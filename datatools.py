"""
Utilitaires de reduction/conversion de donnees du menu Graphics de
StarmacOSX_x.f95 (convertthelli, removestep, elimineGRM, convzmoins,
exportthellier) - toutes MUTENT `ech.mesures` en place, comme leurs
equivalents Fortran (aucune ne dessine quoi que ce soit ; elles vivent dans
le menu Graphics de l'original mais sont en realite des operations sur les
donnees, generalement suivies d'un nouveau trace Zijderveld manuel).

`testinduite` (correction du champ induit parasite, calcul.f:3687-3774) a
ete porte puis RETIRE a la demande explicite de l'utilisateur ("on peut la
retirer, personne ne va penser a faire ces tests pour tester le champ
residuel au niveau de la mesure du Cryo 2G") - fonctionnalite jugee sans
usage reel pour ce laboratoire, pas seulement deplacee de menu.
"""

import os
from dataclasses import replace

from selection import SelectedSample, polere
from calcul import _fortran_e


# ---------------------------------------------------------------------------
# convertthelli / thelli : sequence Thellier (D/N/T/K/S/R/P/L/X/Y/Z/V) ->
# sequence NRM/TRM simple (cod1='D', cod2='p') pour le trace Zijderveld -
# plotorthog.f:780-921
# ---------------------------------------------------------------------------

def convert_thellier_to_nrm(ech: SelectedSample) -> None:
    """Port exact du `select case (mes(i).cod1)` de `convertthelli`
    (plotorthog.f:800-887). D/N/T/K/S copies telles quelles (cod1='D',
    cod2='p'). R/P/L supprimees. Pour X/Y/Z avec cod2 '+'/'-' : fusionnee
    (moyenne) avec la ligne PRECEDENTE seulement si celle-ci a le MEME
    cod1 (2e passage du meme axe) ; sinon supprimee - ce n'est jamais la
    ligne SUIVANTE qui declenche la fusion. Pour V : fusionnee TOUJOURS
    avec la ligne precedente, QUEL QUE SOIT son cod1 - typiquement un 'R'
    juste avant (meme si ce 'R' est par ailleurs supprime de la sortie) :
    c'est ainsi que R et V sont combines pour extraire la NRM (TRM
    partielle appliquee/mesuree en deux positions Z/-Z ; R porte le pas
    thermique, V le controle pTRM, leur moyenne annule le biais directionnel
    entre les deux). Toute autre cod1 : supprimee (case default)."""
    mesures = ech.mesures
    result = []

    def merge(prev, cur):
        x = (prev.x + cur.x) / 2.0
        y = (prev.y + cur.y) / 2.0
        z = (prev.z + cur.z) / 2.0
        q = (prev.q + cur.q) // 2
        s = prev.s if (cur.s == 0 or prev.s == 0) else (cur.s + prev.s) / 2.0
        result.append(replace(cur, x=x, y=y, z=z, q=q, s=s, cod1="D", cod2="p"))

    for i, m in enumerate(mesures):
        if m.cod1 in ("D", "N", "T", "K", "S"):
            result.append(replace(m, cod1="D", cod2="p"))
        elif m.cod1 in ("R", "P", "L"):
            continue
        elif m.cod1 in ("X", "Y", "Z"):
            if m.cod2 in ("+", "-") and i > 0 and mesures[i - 1].cod1 == m.cod1:
                merge(mesures[i - 1], m)
        elif m.cod1 == "V":
            if i > 0:
                merge(mesures[i - 1], m)
    ech.mesures = result


# ---------------------------------------------------------------------------
# removestep : supprime toutes les lignes d'un etape donne, renomme les
# cod2 de R/V/P sequentiellement - plotpaleoint2.f:1918-1979
# ---------------------------------------------------------------------------

_CLEANPALEO_CODES = ("X", "Y", "Z", "Q", "L", "F")


def remove_step(ech: SelectedSample, etape: int) -> int:
    """Port de `removestep` (plotpaleoint2.f:1918-1979), reduit a UN
    echantillon (le Fortran exige nbech==1 sur la TOTALITE des donnees
    chargees ; ici applique au seul `ech` passe en argument). D'abord
    `cleanpaleo` sur X/Y/Z/Q/L/F (supprime les lignes GRM/qualite non
    pertinentes pour l'interpretation paleointensite), puis supprime
    toutes les lignes de l'etape donne. Renumerotation en 3 PASSES
    INDEPENDANTES (la 1ere mesure, mes(1), est exclue de toutes) : les
    'R' recoivent A,B,C... (compteur propre), puis - separement - les 'V'
    recoivent A,B,C... (compteur propre, PAS partage avec R - coincide
    seulement quand R et V sont en nombre egal et dans le meme ordre, cas
    usuel mais pas garanti), puis chaque 'P' reprend le cod2 de la ligne
    precedente. Retourne le nombre de lignes supprimees (cleanpaleo +
    etape)."""
    before = len(ech.mesures)
    ech.mesures = [m for m in ech.mesures if m.cod1 not in _CLEANPALEO_CODES]
    ech.mesures = [m for m in ech.mesures if m.etape != etape]
    removed = before - len(ech.mesures)

    mesures = ech.mesures
    ic = 0
    for m in mesures[1:]:
        if m.cod1 == "R":
            m.cod2 = chr(65 + ic) if ic < 26 else "?"
            ic += 1
    ic = 0
    for m in mesures[1:]:
        if m.cod1 == "V":
            m.cod2 = chr(65 + ic) if ic < 26 else "?"
            ic += 1
    for i in range(1, len(mesures)):
        if mesures[i].cod1 == "P":
            mesures[i].cod2 = mesures[i - 1].cod2
    return removed


# ---------------------------------------------------------------------------
# elimineGRM (elimineGRM_DZ / elimineGRM_ZD) : combine les triplets
# cod2='X','Y','Z' consecutifs (apres un pas 'F') en un seul point corrige -
# plotorthog.f:1536-1846
# ---------------------------------------------------------------------------

def eliminate_grm(ech: SelectedSample, method: int = 1) -> int:
    """Port exact de `elimineGRM_DZ` (method=1, plotorthog.f:1692-1846) /
    `elimineGRM_ZD` (method=2, plotorthog.f:1536-1691), selectionnees par
    le dispatcher `elimineGRM` (StarmacOSX_x.f95:975) selon la commune
    `/prefstarmac/grm`. Seules les lignes cod1='F' sont concernees : celles
    dont cod2 est '+'/'-'/'=' ou dont l'etape vaut 0 sont copiees telles
    quelles ; un triplet consecutif cod2='X','Y','Z' est reduit en UN
    point (method=1 : substitution axe-par-axe x du point X, y du point Y,
    z du point Z ; method=2 : moyenne simple des 3 points) ; toute autre
    ligne 'F' isolee (ni copiee ni fusionnable) est supprimee, comme dans
    le Fortran (aucune des deux conditions du `select case` ne s'applique).
    Toutes les lignes d'un autre cod1 sont copiees inchangees. Retourne le
    nombre de triplets reduits."""
    mesures = ech.mesures
    n = len(mesures)
    result = []
    i = 0
    count = 0
    while i < n:
        m = mesures[i]
        if m.cod1 != "F":
            result.append(m)
            i += 1
            continue
        if m.cod2 in ("+", "-", "=") or m.etape == 0:
            result.append(m)
            i += 1
        elif (i + 2 < n and m.cod2 == "X" and mesures[i + 1].cod2 == "Y"
                and mesures[i + 2].cod2 == "Z"):
            mx, my, mz = mesures[i], mesures[i + 1], mesures[i + 2]
            if method == 1:
                x, y, z = mx.x, my.y, mz.z
            else:
                x = (mx.x + my.x + mz.x) / 3.0
                y = (mx.y + my.y + mz.y) / 3.0
                z = (mx.z + my.z + mz.z) / 3.0
            result.append(replace(mx, x=x, y=y, z=z, cod2="="))
            count += 1
            i += 3
        else:
            i += 1  # ligne 'F' isolee : ni copiee ni fusionnee (comme le Fortran)
    ech.mesures = result
    return count


# ---------------------------------------------------------------------------
# convzmoins : inverse y/z et recode en 'R' les mesures cod1='Z' cod2='-' -
# calcul.f:3929-3951 (citee integralement, algorithme trivial)
# ---------------------------------------------------------------------------

def convert_z_minus(ech: SelectedSample) -> int:
    """Retourne le nombre de mesures converties."""
    count = 0
    for m in ech.mesures:
        if m.cod1 == "Z" and m.cod2 == "-":
            m.y = -m.y
            m.z = -m.z
            m.cod2 = "R"
            count += 1
    return count


# ---------------------------------------------------------------------------
# exportthellier : export au format .tdt (ThellierTool) - fichiers.f:1033-1091
# ---------------------------------------------------------------------------

_THELLIER_STEP_OFFSET = {"R": 0.1, "P": 0.2, "V": 0.5}


def export_thellier_tdt(ech: SelectedSample, out_dir: str) -> str:
    """Port de `exportthellier` (fichiers.f:1033-1091 - version compilee du
    projet, PAS fichiers_mod_magic.f qui n'est pas dans le makefile).
    Applique d'abord `cleanpaleo` sur X/Y/Z/Q/L/F (lignes GRM/qualite non
    pertinentes) SUR UNE COPIE des mesures (ne mute pas `ech`, contrairement
    au Fortran qui mute les donnees globales - l'export ne doit pas avoir
    d'effet de bord sur la selection en memoire). Ecrit `<out_dir>/<id>.tdt`.
    Retourne le chemin ecrit."""
    mesures = [m for m in ech.mesures if m.cod1 not in _CLEANPALEO_CODES]
    path = os.path.join(out_dir, f"{ech.id}.tdt")
    lines = ["Thellier-tdt", f"{(ech.com or '')[:2]}\t0.0\t0.0\t0.0\t0.0"]
    vol = ech.vol or 1.0
    for m in mesures:
        mag, dec, inc = polere(m.x, m.y, m.z)
        rxx = mag * 1.0e3 / vol
        step = m.etape + _THELLIER_STEP_OFFSET.get(m.cod1, 0.0)
        lines.append(f"{ech.id}\t{step:.1f}\t{_fortran_e(rxx, 12, 5)}\t{dec:.2f}\t{inc:.2f}")
    with open(path, "w", encoding="iso-8859-1", errors="replace") as f:
        f.write("\n".join(lines) + "\n")
    return path
