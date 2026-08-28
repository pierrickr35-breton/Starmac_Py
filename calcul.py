"""
Equivalent Python d'une premiere tranche du menu Fortran "Calcul" :
ajustement de droites par ACP (subroutines `linear`/`eigen` de linesplans.f,
pilotees interactivement par `ajuslig`) et statistiques de Fisher
(subroutine `fisher`/`cpolar` de calcul.f, utilisees par `fishmes`/`fishres`).

Travaille sur les SelectedSample/Measurement de selection.py.
"""

import math
import os
import random
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

from selection import Measurement, SelectedSample, apply_orientation, polere, select_samples
from testlect import Pmag


# ---------------------------------------------------------------------------
# Equivalent de la structure /Resultats/ (starmac_OSX.inc), limite aux
# champs utilises par un ajustement de droite (cat1='L').
# ---------------------------------------------------------------------------

@dataclass
class FitResult:
    """Un ajustement de droite (equivalent d'un enregistrement `res`/`tr`)."""
    c: int = 0              # identifiant anti-collision (equivalent res.c), rempli a l'archivage
    id: str = ""
    cin: float = 0.0
    caz: float = 0.0
    dip: float = 0.0
    str_: float = 0.0
    cat1: str = "L"
    cat2: str = " "
    orig: str = "n"        # 'o' = ancre a l'origine, 'n' = non ancre
    demag: str = ""
    numcomp: int = 1
    nb: int = 0
    dec: float = 0.0       # declinaison/inclinaison BRUTES (repere echantillon,
    inc: float = 0.0       # non corrigees - comme res.dec/res.inc en Fortran)
    mad: float = 0.0       # maximum angular deviation (deg)
    step_first: int = 0
    step_last: int = 0
    tx: Tuple[float, float] = (0.0, 0.0)  # extremites du segment (pour tracer
    ty: Tuple[float, float] = (0.0, 0.0)  # la droite ajustee sur un Zijderveld)
    tz: Tuple[float, float] = (0.0, 0.0)
    # Champs specifiques aux resultats agreges "mean:" (moyenne de site,
    # equivalent id="mean: <site>", cat1='F', cat2='i' - Fisher). Inutilises
    # pour un resultat normal (L/P/f/s). par2_mean/par3_mean : memes
    # positions de colonnes que par2/par3 (step_first/step_last), mais des
    # statistiques continues (pas des numeros d'etape entiers) - champs
    # separes pour ne pas perdre leur precision decimale au passage par
    # step_first/step_last (entiers, arrondis) lors d'un aller-retour
    # lecture/ecriture. par4/par5 = 2e paire lat/long ou VGP selon le
    # contexte, `liste` = "codes:c1:c2:..." (les `c` des resultats
    # individuels combines dans la moyenne).
    lat: float = 0.0
    rlong: float = 0.0
    par2_mean: float = 0.0
    par3_mean: float = 0.0
    par4: float = 0.0
    par5: float = 0.0
    liste: str = ""


@dataclass
class FisherStats:
    """Statistiques de Fisher (1953) sur un jeu de directions."""
    n: int
    dec: float
    inc: float
    r: float
    k: float
    a95: float
    csd: float


# ---------------------------------------------------------------------------
# ajuslig / linear / eigen : ajustement de droite par ACP (Kirschvink 1980)
# ---------------------------------------------------------------------------

def linear_fit(
    points: List[Tuple[float, float, float]], anchored: bool, mad_threshold: float = 15.0,
) -> Optional[dict]:
    """Equivalent de la subroutine `linear` : ajustement par ACP d'une droite
    passant (si `anchored`) par l'origine, sur les points (x,y,z) fournis
    (repere echantillon, dans l'ordre des etapes).

    Retourne None si le test de linearite ou le MAD (>`mad_threshold` deg,
    15 par defaut pour ajuslig/Zijderveld) echouent (equivalent de
    itestlin<0 en Fortran), sinon un dict avec direction (vecteur propre
    principal, unitaire), mad, nb, et les extremites (point1, point2) du
    segment ajuste pour affichage. `linearpal` (paleointensite, meme
    algorithme) utilise un seuil de 35 - voir fit_arai_direction.
    """
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    if n == 0:
        return None

    anchor = np.zeros(3) if anchored else pts[-1]
    # Le centre utilise pour l'ACP est l'origine si ancre, sinon le centroide
    # (c'est le comportement d'origine de `linear` : sx/sy/sz ne sont
    # accumules que si `ancr .ne. 'o'`).
    center = np.zeros(3) if anchored else pts.mean(axis=0)

    # Test de linearite : distance directe (premier point -> ancre) comparee
    # a la somme des longueurs des segments consecutifs le long du chemin.
    sl = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
    if anchored:
        sl += float(np.linalg.norm(pts[-1]))
        mom = float(np.linalg.norm(pts[0]))
    else:
        mom = float(np.linalg.norm(pts[0] - anchor))
    linearity_ratio = (mom / sl) if sl > 0 else 1.0
    if sl > 0 and linearity_ratio < 0.5:
        return None
    if mom == 0:
        return None

    diffs = pts - center
    cov = diffs.T @ diffs

    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    l1, l2, l3 = eigvals[order]
    direction = eigvecs[:, order[0]]

    if l1 <= 0:
        return None
    mad = math.degrees(math.atan(math.sqrt(abs((l2 + l3) / l1))))
    if mad > mad_threshold:
        return None

    first_delta = pts[0] - anchor
    if np.dot(direction, first_delta) < 0:
        direction = -direction

    dot = float(np.dot(direction, center - anchor))
    point1 = center - dot * first_delta / mom
    point2 = point1 + mom * direction

    return {
        "direction": direction,
        "mad": mad,
        "nb": n,
        "linearity_ratio": linearity_ratio,
        "point1": tuple(point1),
        "point2": tuple(point2),
    }


def fit_line(
    ech: SelectedSample,
    jdeb: int,
    jfin: int,
    anchored: bool = True,
    numcomp: int = 1,
) -> Optional[FitResult]:
    """Equivalent non-interactif de `ajuslig` pour un seul segment : ajuste
    une droite sur ech.mesures[jdeb-1:jfin] (indices 1-based, comme affiches
    par list_measurements). Retourne None si le fit echoue (MAD>15 ou test
    de linearite, comme itestlin<0 en Fortran)."""
    subset = ech.mesures[jdeb - 1:jfin]
    fit = linear_fit([(m.x, m.y, m.z) for m in subset], anchored)
    if fit is None:
        return None

    dx, dy, dz = fit["direction"]
    _, dec, inc = polere(dx, dy, dz)

    return FitResult(
        id=ech.id,
        cin=ech.cin,
        caz=ech.caz,
        dip=ech.dip,
        str_=ech.str_,
        cat1="L",
        cat2=" ",
        orig="o" if anchored else "n",
        demag=subset[-1].cod1,
        numcomp=numcomp,
        nb=fit["nb"],
        dec=dec,
        inc=inc,
        mad=fit["mad"],
        step_first=subset[0].etape,
        step_last=subset[-1].etape,
        tx=(float(fit["point1"][0]), float(fit["point2"][0])),
        ty=(float(fit["point1"][1]), float(fit["point2"][1])),
        tz=(float(fit["point1"][2]), float(fit["point2"][2])),
    )


def _read_int(prompt: str, default: Optional[int] = None) -> Optional[int]:
    chaine = input(prompt).strip()
    if not chaine:
        return default
    try:
        return int(chaine)
    except ValueError:
        return default


def fit_line_interactive(selected: List[SelectedSample]) -> List[FitResult]:
    """Equivalent interactif de `ajuslig` : pour chaque echantillon
    selectionne, affiche la liste des etapes, demande le premier/dernier
    step a utiliser (0 = passer a l'echantillon suivant), l'ancrage a
    l'origine et le numero de composante, puis affiche le resultat."""
    results: List[FitResult] = []

    for ech in selected:
        if len(ech.mesures) < 2:
            continue

        steps = [m.etape for m in ech.mesures]
        print(f"\n{ech.id} :")
        print("  " + "  ".join(f"{i + 1}:{s}" for i, s in enumerate(steps)))

        jdeb = _read_int(" First step (1-based, 0=skip): ", 0)
        if not jdeb:
            continue
        jfin = _read_int(" Final step: ", None)
        if jfin is None or jfin < jdeb or jfin > len(steps):
            continue

        ancr = input(" anchored to the origin (Y/n)? ").strip().lower()
        anchored = ancr != "n"

        numcomp = _read_int(" component number (1-9)? ", 1) or 1

        fit = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=numcomp)
        if fit is None:
            print(" negative linear test (MAD too high or non-linear trend)")
            continue

        print(
            f"  dec={fit.dec:5.1f}  inc={fit.inc:5.1f}  mad={fit.mad:5.1f}  "
            f"nb points: {fit.nb}"
        )
        results.append(fit)

    return results


# ---------------------------------------------------------------------------
# fisher / cpolar : statistiques de Fisher (1953)
# ---------------------------------------------------------------------------

def fisher_mean(directions: List[Tuple[float, float]]) -> FisherStats:
    """Equivalent de `fisher`+`cpolar` : moyenne de Fisher d'une liste de
    directions (dec, inc) en degres."""
    n = len(directions)
    if n == 0:
        raise ValueError("empty list of directions")

    xsum = ysum = zsum = 0.0
    for dec, inc in directions:
        decr, incr = math.radians(dec), math.radians(inc)
        xsum += math.cos(incr) * math.cos(decr)
        ysum += math.cos(incr) * math.sin(decr)
        zsum += math.sin(incr)

    r = math.sqrt(xsum ** 2 + ysum ** 2 + zsum ** 2)
    if r == 0:
        raise ValueError("zero resultant - directions undefined")
    x, y, z = xsum / r, ysum / r, zsum / r

    mean_inc = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    mean_dec = math.degrees(math.atan2(y, x)) % 360.0

    if n > 1 and (n - r) > 0:
        k = (n - 1) / (n - r)
    else:
        k = 1.0

    if n > 1:
        fac = 1.0 / (n - 1)
        p = 20 ** fac - 1
        a = 1 - ((n - r) / r) * p
        a95 = math.degrees(math.acos(max(-1.0, min(1.0, a)))) if abs(a) < 1.0 else 90.0
    else:
        a95 = 180.0

    if n > 1:
        tet = 0.0
        for dec, inc in directions:
            decr, incr = math.radians(dec), math.radians(inc)
            xi = math.cos(incr) * math.cos(decr)
            yi = math.cos(incr) * math.sin(decr)
            zi = math.sin(incr)
            chord = math.sqrt((xi - x) ** 2 + (yi - y) ** 2 + (zi - z) ** 2)
            tet += (2 * math.asin(min(1.0, chord / 2))) ** 2
        csd = math.degrees(math.sqrt(tet / (n - 1)))
    else:
        csd = 0.0

    return FisherStats(n=n, dec=mean_dec, inc=mean_inc, r=r, k=k, a95=a95, csd=csd)


def fisher_from_measurements(selected: List[SelectedSample], orientation: int = 1) -> FisherStats:
    """Equivalent de `fishmes` : moyenne de Fisher calculee directement sur
    TOUTES les mesures de la selection (chaque etape de demagnetisation
    compte comme une direction independante), apres correction d'orientation."""
    directions = []
    for ech in selected:
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            _, dec, inc = polere(xx, yy, zz)
            directions.append((dec, inc))
    return fisher_mean(directions)


def fisher_from_results(results: List[FitResult], orientation: int = 1) -> FisherStats:
    """Equivalent (partiel) de `fishres` : moyenne de Fisher sur les
    ajustements de droite (cat1='L'), en reappliquant la correction
    d'orientation de chaque resultat (son propre cin/caz/dip/str). Les
    ajustements de plan (cat1='P', resolus par intersection de grands
    cercles dans le Fortran via `fishgc`) ne sont pas geres ici."""
    directions = []
    for res in results:
        if res.cat1 != "L":
            continue
        zz = math.sin(math.radians(res.inc))
        xx = math.cos(math.radians(res.inc)) * math.cos(math.radians(res.dec))
        yy = math.cos(math.radians(res.inc)) * math.sin(math.radians(res.dec))
        xx, yy, zz = apply_orientation(xx, yy, zz, res, orientation)
        _, dec, inc = polere(xx, yy, zz)
        directions.append((dec, inc))
    if not directions:
        raise ValueError("no line-type result ('L') in the list")
    return fisher_mean(directions)


# ---------------------------------------------------------------------------
# ajusplans / planar : ajustement de plan (grand cercle) par ACP - meme
# principe que linear_fit, mais le pole du plan est le vecteur propre de la
# PLUS PETITE valeur propre (linear_fit prend la plus grande). Pas de test de
# linearite (n'a pas de sens pour un plan) ; seuil MAD 25 deg (pas 15).
# ---------------------------------------------------------------------------

def planar_fit(points: List[Tuple[float, float, float]], normalize: bool = True) -> Optional[dict]:
    """Equivalent de la subroutine `planar` : ajustement par ACP d'un plan
    (grand cercle) passant par les points fournis. `normalize` : equivalent
    de `norm=='o'` (points ramenes sur la sphere unite avant l'ACP - option
    par defaut de `ajusplans`). Retourne None si MAD>25 (equivalent
    itestplan=-1)."""
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    if n == 0:
        return None
    if normalize:
        norms = np.linalg.norm(pts, axis=1)
        norms[norms == 0] = 1.0
        pts = pts / norms[:, None]

    center = pts.mean(axis=0)
    diffs = pts - center
    cov = diffs.T @ diffs

    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    l1, l2, l3 = eigvals[order]
    if l1 <= 0 or l2 <= 0:
        return None
    pole = eigvecs[:, order[2]]  # plus petite valeur propre = pole du plan

    mad = math.degrees(math.atan(math.sqrt(abs(l3 / l2 + l3 / l1))))
    if mad > 25.0:
        return None

    return {"pole": pole, "mad": mad, "nb": n}


def fit_plane(
    ech: SelectedSample,
    jdeb: int,
    jfin: int,
    normalize: bool = True,
    numcomp: int = 1,
) -> Optional[FitResult]:
    """Equivalent non-interactif de `ajusplans` pour un seul segment (indices
    1-based, comme fit_line). `ancr` est toujours 'o' cote Fortran (hardcode,
    pas de prompt) - reproduit ici (orig="o"), sauf si l'appelant l'ecrase
    explicitement (cas de ajusligredo, ou le fichier redo peut demander 'n')."""
    subset = ech.mesures[jdeb - 1:jfin]
    fit = planar_fit([(m.x, m.y, m.z) for m in subset], normalize=normalize)
    if fit is None:
        return None

    px, py, pz = fit["pole"]
    _, dec, inc = polere(px, py, pz)

    return FitResult(
        id=ech.id,
        cin=ech.cin,
        caz=ech.caz,
        dip=ech.dip,
        str_=ech.str_,
        cat1="P",
        cat2=" ",
        orig="o",
        demag=subset[-1].cod1,
        numcomp=numcomp,
        nb=fit["nb"],
        dec=dec,
        inc=inc,
        mad=fit["mad"],
        step_first=subset[0].etape,
        step_last=subset[-1].etape,
    )


# ---------------------------------------------------------------------------
# ajusfisher : direction moyenne de Fisher sur une plage d'etapes (au lieu
# d'un ajustement PCA) - reutilise fisher_mean() (deja verifiee identique a
# `fisher`/`cpolar`). Repli sur une direction brute unique (label 735 du
# Fortran) pour les echantillons a moins de 3 mesures.
# ---------------------------------------------------------------------------

def fit_single_direction(ech: SelectedSample, index: int = -1) -> FitResult:
    """Equivalent de la branche 'direction unique' de `ajusfisher` (label
    735) : utilise directement le vecteur d'UNE mesure (par defaut la
    derniere, comme le Fortran `jdeb=ech(i).ifin`) comme direction, sans
    moyenne de Fisher."""
    m = ech.mesures[index]
    _, dec, inc = polere(m.x, m.y, m.z)
    return FitResult(
        id=ech.id, cin=ech.cin, caz=ech.caz, dip=ech.dip, str_=ech.str_,
        cat1="s", cat2="s", orig="n", demag=m.cod1, numcomp=1,
        nb=1, dec=dec, inc=inc, mad=0.0,
        step_first=m.etape, step_last=m.etape,
    )


def fit_fisher_direction(
    ech: SelectedSample, jdeb: int, jfin: int, numcomp: int = 1
) -> FitResult:
    """Equivalent de `ajusfisher` (branche moyenne de Fisher) sur
    ech.mesures[jdeb-1:jfin]. `tx` porte (kappa, csd) - meme detournement de
    champ que le Fortran (res.tx(1)=rk, res.tx(2)=stv), pas des coordonnees
    de segment."""
    subset = ech.mesures[jdeb - 1:jfin]
    directions = []
    for m in subset:
        _, dec, inc = polere(m.x, m.y, m.z)
        directions.append((dec, inc))
    stats = fisher_mean(directions)

    return FitResult(
        id=ech.id, cin=ech.cin, caz=ech.caz, dip=ech.dip, str_=ech.str_,
        cat1="f", cat2=" ", orig="n",
        demag=subset[-1].cod1, numcomp=numcomp, nb=stats.n,
        dec=stats.dec, inc=stats.inc, mad=stats.a95,
        step_first=subset[0].etape, step_last=subset[-1].etape,
        tx=(stats.k, stats.csd),
    )


# ---------------------------------------------------------------------------
# ajusligauto : version automatisee de ajuslig - une droite par echantillon
# sur la totalite de ses mesures, meme ancrage pour tous, sauvegarde sans
# confirmation (echantillons a moins de 2 mesures ignores).
# ---------------------------------------------------------------------------

def fit_lines_auto(selected: List[SelectedSample], anchored: bool = True) -> List[FitResult]:
    """Equivalent de `ajusligauto`. `numcomp` vaut 1 si ancre, 2 sinon (meme
    regle que le Fortran : `res.numcomp=2 ; if(ancr.ne."n") res.numcomp=1`)."""
    numcomp = 1 if anchored else 2
    results = []
    for ech in selected:
        if len(ech.mesures) < 2:
            continue
        fit = fit_line(ech, 1, len(ech.mesures), anchored=anchored, numcomp=numcomp)
        if fit is not None:
            results.append(fit)
    return results


# ---------------------------------------------------------------------------
# ajusligredo : rejoue des ajustements (droite ou plan) a partir d'un fichier
# texte "redo" (un ajustement par ligne). Selection d'echantillon par id EXACT
# (equivalent simplifie de `seloneech` - pas de joker '*'/'?' ici, contrai-
# rement a selection.select_samples qui les gere deja pour le menu Selection
# donnees). Plage de steps resolue par VALEUR d'etape (pas position) : jdeb/
# jfin = dernier index dont l'etape est <= tempmin/tempmax (balayage lineaire
# "dernier qui correspond", pas une recherche exacte - fidele au Fortran).
# ---------------------------------------------------------------------------

def _resolve_step_range_by_value(
    ech: SelectedSample, tempmin: int, tempmax: int
) -> Optional[Tuple[int, int]]:
    jdeb = jfin = None
    for idx, m in enumerate(ech.mesures, start=1):
        if m.etape <= tempmin:
            jdeb = idx
        if m.etape <= tempmax:
            jfin = idx
    if jdeb is None or jfin is None:
        return None
    return jdeb, jfin


def fit_from_redo_file(
    pmag_list: List[Pmag], redo_lines: List[str]
) -> List[FitResult]:
    """Equivalent de `ajusligredo`. Chaque ligne : 'sample L|P o|n tempmin
    tempmax numcomp' (ex: 'NWD1-11A L n 400 560 2'). Lignes consecutives
    identiques ignorees (la seconde), comme le Fortran (chaine2==chaine).
    Echantillons a moins de 3 mesures ignores (comme `nbmes<3` cote Fortran).

    `pmag_list` : la TOTALITE des donnees chargees (equivalent self.donnees/
    pmag(:)), PAS une selection prealable - chaque ligne recharge elle-meme
    l'echantillon qu'elle designe via `call seloneech(samnum)` (nbech=0;
    nbmes=0; recherche dans pmag(:) par id, etapmin=0/etapmax=9999/
    demag='*'), independamment de toute selection en cours. Un fichier redo
    peut donc rejouer des ajustements sur des specimens jamais selectionnes
    manuellement au prealable."""
    results: List[FitResult] = []
    last_line: Optional[str] = None

    for raw in redo_lines:
        line = raw.strip()
        if not line or line == last_line:
            last_line = line
            continue
        last_line = line

        parts = line.split()
        if len(parts) < 6:
            continue
        samnum, lp, ancr, tempmin_s, tempmax_s, numcomp_s = parts[:6]
        if lp not in ("L", "P"):
            continue
        matches = select_samples(pmag_list, samnum, verbose=False)
        ech = matches[0] if matches else None
        if ech is None or len(ech.mesures) < 3:
            continue
        try:
            tempmin, tempmax, numcomp = int(tempmin_s), int(tempmax_s), int(numcomp_s)
        except ValueError:
            continue

        rng = _resolve_step_range_by_value(ech, tempmin, tempmax)
        if rng is None:
            continue
        jdeb, jfin = rng
        anchored = ancr.lower() != "n"

        if lp == "L":
            fit = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=numcomp)
        else:
            fit = fit_plane(ech, jdeb, jfin, normalize=True, numcomp=numcomp)
            if fit is not None:
                # contrairement a ajusplans (ancr toujours 'o'), ajusligredo
                # derive orig de l'ancr du fichier meme pour un plan (le
                # calcul PCA lui-meme ignore ancr, seul le champ change).
                fit.orig = "o" if anchored else "n"
        if fit is not None:
            results.append(fit)

    return results


# ---------------------------------------------------------------------------
# Resultats : init / listage (equivalents de initres / lisres, limites aux
# lignes - lisres gere aussi les plans et n'est donc que partiellement porte)
# ---------------------------------------------------------------------------

def init_results() -> List[FitResult]:
    """Equivalent de `initres` : reinitialise la liste des resultats."""
    print("\n-- Init the results list -- ")
    return []


def _correct_dec_inc(res: FitResult, orientation: int) -> Tuple[float, float]:
    """Equivalent de `correct(rd,ri,cdip,caz,dip,str)` (StarUtil.f:66),
    reapplique a la direction BRUTE stockee dans un resultat (res.dec/
    res.inc, repere echantillon d'origine) - duplique de
    stereo.py:_correct_dec_inc pour eviter un import circulaire (stereo.py
    importe deja FitResult depuis ce module).

    Bug corrige (signale par l'utilisateur sur des moyennes de site
    importees de MagIC - sites.txt, deja calculees en In-Situ/apres
    pendage via dir_tilt_correction : "les moyennes de site ne se
    calculent pas en coordonnees echantillons mais en In situ et/ou
    bedding correction") : un resultat "mean:" est DEJA fige dans UNE
    orientation precise (par3_mean : 1=echantillon/2=in-situ/3=apres
    pendage), pas dans le repere echantillon brut comme un resultat L/P -
    Starmac n'a pas les directions specimen individuelles pour le
    reprojeter dans une AUTRE orientation (une moyenne de Fisher ne
    s'inverse pas comme une simple rotation, contrairement a un fit
    L/P individuel - voir recompute_fit_geometry, qui resout le meme
    genre de probleme mais UNIQUEMENT parce qu'il peut retrouver le
    specimen d'origine). Sans ce garde-fou, `apply_orientation` utilisait
    cin=caz=0.0 (jamais renseignes sur une moyenne) et `corfor`
    introduisait une rotation parasite (~90 deg, meme symptome que le bug
    deja corrige pour les resultats individuels charges depuis .pmagres) -
    et meme avec orientation=1 (repere echantillon, un no-op), la valeur
    stockee etait affichee SANS AVERTIR qu'elle est en realite dans une
    AUTRE orientation. Retourne desormais TOUJOURS la valeur stockee
    telle quelle pour un "mean:", quelle que soit `orientation` demandee -
    voir list_results pour l'etiquette (Sa)/(IS)/(TC) AFFICHEE, qui
    reflete `par3_mean` (l'orientation REELLE de la moyenne), pas la
    colonne demandee, pour ne pas laisser croire a une reprojection qui
    n'a pas lieu."""
    if res.id[:5] == "mean:":
        return res.dec, res.inc
    incr, decr = math.radians(res.inc), math.radians(res.dec)
    x = math.cos(incr) * math.cos(decr)
    y = math.cos(incr) * math.sin(decr)
    z = math.sin(incr)
    xx, yy, zz = apply_orientation(x, y, z, res, orientation)
    _, dec, inc = polere(xx, yy, zz)
    return dec, inc


_ORIENT_HEADER = {
    1: "results in sample coordinates",
    2: "results in in-situ coordinates",
    3: "results after tilt correction",
}


_ORIENT_MODE_TAG = {1: "Sa", 2: "IS", 3: "TC"}


def _e95(a95: float, inc: float) -> float:
    """Equivalent exact de la formule e95 de `lisres` (dataselect.f, juste
    avant le format 201) : deduite de a95 et de l'inclinaison de la
    moyenne (paleolatitude via VGP), PAS d'une donnee stockee - calculable
    a l'affichage, comme dans le Fortran."""
    if a95 <= 0:
        return 0.0
    ai = math.radians(inc)
    pl = 1.57095 - math.atan(0.5 * math.tan(ai))
    dm = a95 * math.sin(pl) / math.cos(ai)
    dp = 2 * a95 * (1 / (1 + 3 * math.cos(ai) ** 2))
    return (dm + dp) / 2


def _mean_line_plane_counts(r: FitResult, results: List[FitResult]) -> Optional[Tuple[int, int]]:
    """Equivalent de "L nnn  P nnn" dans `lisres` : compte, parmi les
    resultats specimen PRESENTS dans `results` (meme liste passee a
    list_results) dont le `c` figure dans `r.liste`, combien sont des
    droites (cat1=='L') et des plans (cat1=='P'). Contrairement au
    Fortran (qui stocke ces 2 comptes directement dans tr(i).tx(1)/tx(2)
    a l'archivage), le format .pmagres ne les stocke pas - tx(1) y est
    deja pris par k (voir _format_mean_line) - donc on les retrouve ici
    par recoupement via `liste`. Retourne None si `results` ne contient
    aucun des specimens references (ex: `results` == uniquement les
    moyennes, chargees avec carselect='m') - dans ce cas le compte est
    inconnu plutot que faussement 0/0."""
    ids = [tok for tok in r.liste.replace("codes:", "").split(":") if tok.strip()]
    if not ids:
        return None
    by_c = {str(other.c): other for other in results if other.id[:5] != "mean:"}
    matched = [by_c[i] for i in ids if i in by_c]
    if not matched:
        return None
    nlines = sum(1 for m in matched if m.cat1 == "L")
    nplanes = sum(1 for m in matched if m.cat1 == "P")
    return nlines, nplanes


def _mean_site_strike_dip(r: FitResult, donnees) -> Optional[Tuple[float, float]]:
    """Equivalent de la boucle `do inech=1,nb_ech ... if(tr(i).id(7:12)==
    pmag(inech).id(1:6))` : strike/pendage du site, lus sur N'IMPORTE
    QUEL specimen de `donnees` dont l'id commence par le nom du site
    (le pendage est une propriete du site, identique pour tous ses
    specimens). Retourne None si `donnees` n'est pas fourni ou si aucun
    specimen du site n'y est trouve."""
    if not donnees:
        return None
    site = _mean_site_name(r.id)
    ech = next((s for s in donnees if s.id.startswith(site)), None)
    if ech is None:
        return None
    return ech.str_, ech.dip


def list_results(results: List[FitResult], orientation: int = 1, donnees=None) -> str:
    """Equivalent de `lisres` (dataselect.f:690-820) : tableau des
    ajustements. Comme le Fortran (bloc `select case (iorient)` lignes
    722-784), dec/inc AFFICHES sont recalcules a partir de la direction
    BRUTE stockee (repere echantillon d'origine, jamais modifiee) en
    appliquant la correction d'orientation COURANTE - la meme direction
    stockee s'affiche donc differemment selon iorient (échantillon/in
    situ/apres pendage), au lieu de rester figee dans le repere
    echantillon (voir recompute_fit_geometry pour le pre-requis : cin/
    caz/dip/str_ doivent avoir ete completes depuis le specimen d'origine
    pour un resultat charge depuis .pmagres, sinon cette correction est
    fausse).

    Ligne "mean:" (moyenne de site) : dec/inc affiches NE SONT PAS
    recalcules (voir _correct_dec_inc) - une moyenne est deja figee dans
    l'orientation `par3_mean`, affichee avec SA PROPRE etiquette
    (Sa)/(IS)/(TC) plutot que celle demandee en tete de tableau, pour ne
    pas laisser croire a une reprojection qui n'a pas lieu - demande
    explicite utilisateur sur des moyennes importees de MagIC, deja en
    In-Situ/apres pendage : "les moyennes de site ne se calculent pas en
    coordonnees echantillons mais en In situ et/ou bedding correction".
    Memes autres colonnes que `lisres`, y compris (contrairement a une
    version anterieure de cette fonction) :
    - e95 : calcule ici (`_e95`), pure fonction de a95/inclinaison, comme
      dans le Fortran - pas une donnee stockee.
    - "L nnn  P nnn" (compte lignes/plans composant la moyenne) : deduit
      par recoupement de `r.liste` (les `c` des specimens composants)
      CONTRE `results` lui-meme - fiable seulement si `results` contient
      aussi ces specimens (typiquement `load_results(..., carselect='s')`,
      qui les inclut expres pour cet usage) ; sinon affiche "?" plutot
      qu'un faux "0/0" (voir `_mean_line_plane_counts`).
    - strike/dip du site : lus sur un specimen de `donnees` (parametre
      optionnel, ex. self.donnees) dont l'id commence par le nom du
      site ; affiche "?" si `donnees` n'est pas fourni.
    Le Fortran stocke lignes/plans directement dans tr(i).tx(1)/tx(2) a
    l'archivage - le format .pmagres a repurpose ce slot pour k (voir
    _format_mean_line, k=tx[0]), d'ou le recoupement via `liste` ici au
    lieu d'une simple lecture."""
    tag = _ORIENT_MODE_TAG.get(orientation, "Sa")
    lines = [
        _ORIENT_HEADER.get(orientation, _ORIENT_HEADER[1]),
        f"     Sample          comp  cat  orig  demag   step1  stepn   nb   dec    inc   mad   ({tag})",
    ]
    for i, r in enumerate(results, start=1):
        dec, inc = _correct_dec_inc(r, orientation)
        if r.id[:5] == "mean:":
            counts = _mean_line_plane_counts(r, results)
            lp_txt = f"L{counts[0]:3d}  P{counts[1]:3d}" if counts else "L  ?  P  ?"
            e95 = _e95(r.mad, inc)
            strdip = _mean_site_strike_dip(r, donnees)
            strdip_txt = f"str={strdip[0]:5.1f} dip={strdip[1]:4.1f}" if strdip else "str=?  dip=?"
            own_tag = _ORIENT_MODE_TAG.get(int(r.par3_mean), "?")
            lines.append(
                f"{i:4d}: {r.id:<13s}{r.numcomp:5d}   {r.cat1}{r.cat2}    {r.orig}     {r.demag:<3s}"
                f"  {lp_txt}  {r.nb:4d}  {dec:6.1f} {inc:6.1f} ({own_tag})  a95={r.mad:5.1f} e95={e95:5.1f}"
                f"  k={r.tx[0]:8.1f}  lat={r.lat:9.5f} lon={r.rlong:9.5f}"
                f"  VGP=({r.par4:6.1f},{r.par5:6.1f})  {strdip_txt}"
                f"   [{r.liste}]"
            )
        else:
            lines.append(
                f"{i:4d}: {r.id:<13s}{r.numcomp:5d}   {r.cat1}{r.cat2}    {r.orig}     {r.demag:<3s}"
                f"  {r.step_first:5d}  {r.step_last:5d}  {r.nb:4d}  {dec:6.1f} {inc:6.1f}  {r.mad:5.1f}"
            )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# archivres / selres : sauvegarde et lecture du fichier ".r" (equivalent de
# fichiers.f:archivres et dataselect.f:selres). Format fixe (positions de
# colonnes exactes, pas un simple split) pour rester compatible avec les
# fichiers .r reels produits par Starmac_OSX - structure /Resultats/ complete
# (starmac_OSX.inc lignes 37-47) : c, id, cin, caz, dip, str, cat1, cat2,
# orig, demag, numcomp, nb, dec, inc, par1(mad), par2(step1), par3(stepn),
# tx(2), ty(2), tz(2). Le champ `liste` (cas special "mean:", agregats
# multi-echantillons) n'est pas gere ici - hors perimetre de ce portage.
# ---------------------------------------------------------------------------

def _fortran_e(value: float, width: int = 10, decimals: int = 3) -> str:
    """Equivalent d'un champ Fortran ExW.d (ex E10.3) : mantisse "0.ddd"
    (toujours <1 en valeur absolue), exposant signe sur 2 chiffres - PAS le
    style "d.ddE+xx" de Python/C (`%.3e`). Utilise par datatools.py (export
    detaille) - PLUS par le format .r (qui ne stocke plus tx/ty/tz, voir
    plus bas), garde ici comme utilitaire generique."""
    if value == 0.0:
        mantissa, exp = 0.0, 0
    else:
        sign = -1.0 if value < 0 else 1.0
        av = abs(value)
        exp = math.floor(math.log10(av)) + 1
        mantissa = av / (10.0 ** exp)
        if round(mantissa, decimals) >= 1.0:
            mantissa /= 10.0
            exp += 1
        mantissa *= sign
    sign_char = "-" if mantissa < 0 else " "
    digits = f"{abs(mantissa):.{decimals}f}"[2:]  # "0.123" -> "123"
    exp_sign = "+" if exp >= 0 else "-"
    text = f"{sign_char}0.{digits}E{exp_sign}{abs(exp):02d}"
    return text.rjust(width)[:width]


def results_path_for(data_path: str) -> str:
    """Equivalent de `filr=fil1(1:(nlen-4))//'.r'` : nom du fichier
    resultats derive du fichier de donnees charge (meme dossier, extension
    '.pmagres' - anciennement '.r', renomme a la demande explicite de
    l'utilisateur pour eviter toute ambiguite avec d'autres usages de
    '.r')."""
    base, _ext = os.path.splitext(data_path)
    return base + ".pmagres"


# ---------------------------------------------------------------------------
# Format .pmagres (anciennement .r, renomme a la demande de l'utilisateur)
# reorganise (demande explicite utilisateur, remplace l'ancien
# format Fortran a positions de colonnes fixes/E10.3) : deux sections
# separees dans le MEME fichier - "specimen results" (un ajustement par
# ligne) puis "site mean results" (une moyenne de site par ligne),
# tab-separees, cles auto-descriptives. Tel que specifie par l'utilisateur
# (fichier exemple_Pmag.r) :
# - PAS de tx/ty/tz (extremites du segment ajuste, pour le trace sur un
#   Zijderveld) : "to draw the line on the zijderveld plot we can redo
#   the calculation" - voir recompute_fit_geometry, appelee au chargement
#   plutot que stockee.
# - dec/inc restent en repere ECHANTILLON (pas de correction d'orientation
#   stockee : "as we have the orientation of the samples in the main
#   file, they may not be needed in the result file").
# - `random` (l'ancien `c:NNNNN`) EST conserve - PAS decoratif : sert de
#   cle de reference stable, un specimen pouvant porter plusieurs
#   ajustements candidats (ex. deux droites a des plages d'etapes
#   differentes) et la moyenne de site devant pouvoir dire PRECISEMENT
#   lesquels elle combine (`included_samples` = "codes:<random>:<random>:...").
# - cat2 abandonne (redondant - jamais filtre nulle part, cat1 suffit).
# - Insertion "au bon endroit" (pas toujours en fin de fichier comme
#   avant) : un nouveau resultat specimen est insere a la fin de la
#   SECTION 1 (juste avant la section "site mean", pas apres) ; une
#   nouvelle moyenne va en fin de fichier (section 2, toujours derniere).
#   Les lignes "!" (annotations manuelles de l'utilisateur, ex.
#   "! primary magnetization") et les lignes vides sont preservees telles
#   quelles - demande explicite ("In fortran I was writing at the end of
#   the file" -> ne plus faire ca aveuglement).
# ---------------------------------------------------------------------------

_SPECIMEN_HEADER = "#specimen results in sample coordinates"
_MEAN_HEADER = "#site mean results"

# (libelle, largeur) - espaces-alignees comme le fichier exemple fourni par
# l'utilisateur (exemple_Pmag.r), PAS des tabulations : plus proche de
# l'habitude de lecture Fortran a colonnes fixes (voir aussi .prmag). La
# derniere colonne (random / included_samples) n'est pas paddee (longueur
# variable, rien apres elle sur la ligne).
_SPECIMEN_FIELDS = [
    ("specimen", 12), ("step1", 7), ("step2", 7), ("L/P", 5), ("anc/not", 8),
    ("demag", 7), ("comp", 5), ("n", 6), ("dec", 7), ("inc", 7),
    ("mad", 7), ("random", 8),
]
_MEAN_FIELDS = [
    ("site", 9), ("type", 6), ("n", 5), ("dec", 7), ("inc", 7), ("a95", 7),
    ("k", 9), ("IS/TC", 7), ("lat_site", 12), ("long_site", 12),
    ("VGP_lat", 8), ("VGP_lon", 8), ("included_samples", 0),
]

# "anc"/"not" (anchored/not-anchored) au lieu de 'o'/'n' (francais
# origine/non) dans la colonne "anc/not" du fichier .pmagres - demande
# explicite utilisateur. Le champ interne `FitResult.orig` garde 'o'/'n'
# (utilise ainsi partout ailleurs dans calcul.py/app.py/stereo.py) - seule
# la serialisation fichier change. La lecture reste retro-compatible avec
# d'anciens fichiers .pmagres deja ecrits avec 'o'/'n'.
_ANCHOR_TO_FILE_CODE = {"o": "anc", "n": "not"}
_FILE_CODE_TO_ANCHOR = {"anc": "o", "not": "n", "o": "o", "n": "n"}

# 0/100 (convention MagIC dir_tilt_correction : geographique/apres pendage
# complet) au lieu de 2.0/3.0 dans la colonne "IS/TC" du fichier .pmagres -
# demande explicite utilisateur ("for the mean, can we use 0 and 100 for
# IS and full TC like in magic"). Le champ interne `par3_mean` garde la
# convention Starmac existante (1/2/3, partagee avec `self.orientation`
# ailleurs dans l'appli, voir _ORIENT_MODE_TAG) - seule la serialisation
# fichier change, comme _ANCHOR_TO_FILE_CODE ci-dessus. "1" (coordonnees
# echantillon) n'est normalement plus produit pour une moyenne (voir
# extract_magic.magic_site_means) mais reste gere en lecture/ecriture pour
# ne pas casser d'anciens fichiers .pmagres.
_ORIENT_TO_FILE_CODE = {1.0: "1", 2.0: "0", 3.0: "100"}
_FILE_CODE_TO_ORIENT = {"1": 1.0, "0": 2.0, "100": 3.0}


def _parse_anchor_token(token: str) -> str:
    t = token.strip().lower()
    return _FILE_CODE_TO_ANCHOR.get(t, t[:1] if t else "n")


def _row(fields: List[Tuple[str, int]], values: List[str]) -> str:
    parts = []
    for (_label, width), value in zip(fields, values):
        parts.append(value.ljust(width) if width else value)
    return "".join(parts).rstrip()


def _cols_header(fields: List[Tuple[str, int]]) -> str:
    """En-tete de colonnes, prefixee "#" - le "#" REMPLACE un caractere de
    remplissage (pas un ajout) pour que les colonnes suivantes restent
    alignees avec les lignes de donnees en dessous (sans prefixe) - meme
    convention que exemple_Pmag.r ("#specimen" = 9 caracteres, comme
    "14NQ0101A")."""
    labels = [label for label, _w in fields]
    labels[0] = "#" + labels[0]
    return _row(fields, labels)


def _fmt1(value: float, decimals: int = 1) -> str:
    return f"{float(value):.{decimals}f}"


def _mean_site_name(res_id: str) -> str:
    return res_id[6:] if res_id[:6].lower() == "mean: " else res_id


def _format_specimen_line(res: FitResult) -> str:
    return _row(_SPECIMEN_FIELDS, [
        res.id, _fmt1(res.step_first), _fmt1(res.step_last),
        res.cat1, _ANCHOR_TO_FILE_CODE.get(res.orig, res.orig), res.demag,
        str(res.numcomp), str(res.nb),
        _fmt1(res.dec), _fmt1(res.inc), _fmt1(res.mad), str(res.c),
    ])


def _format_mean_line(res: FitResult) -> str:
    tilt_code = _ORIENT_TO_FILE_CODE.get(res.par3_mean, _fmt1(res.par3_mean, 0))
    return _row(_MEAN_FIELDS, [
        _mean_site_name(res.id), "Fi", str(res.nb),
        _fmt1(res.dec), _fmt1(res.inc), _fmt1(res.mad), _fmt1(res.tx[0]),
        tilt_code, f"{res.lat:.5f}", f"{res.rlong:.5f}",
        _fmt1(res.par4), _fmt1(res.par5), res.liste,
    ])


def _format_result_line(res: FitResult) -> str:
    """Equivalent de `archivres` : choisit le format specimen ou "mean:"
    (si `res.id[:5]=="mean:"`) - meme test que le Fortran d'origine."""
    if res.id[:5] == "mean:":
        return _format_mean_line(res)
    return _format_specimen_line(res)


def _parse_specimen_line(parts: List[str]) -> Optional[FitResult]:
    if len(parts) < 12:
        return None
    try:
        return FitResult(
            id=parts[0].strip(),
            step_first=int(round(float(parts[1]))), step_last=int(round(float(parts[2]))),
            cat1=parts[3].strip(), orig=_parse_anchor_token(parts[4]), demag=parts[5].strip(),
            numcomp=int(parts[6]), nb=int(parts[7]),
            dec=float(parts[8]), inc=float(parts[9]), mad=float(parts[10]),
            c=int(parts[11]),
        )
    except ValueError:
        return None


def _parse_mean_line(parts: List[str]) -> Optional[FitResult]:
    if len(parts) < 13:
        return None
    try:
        tilt_token = parts[7].strip()
        orientation = _FILE_CODE_TO_ORIENT.get(tilt_token)
        if orientation is None:
            orientation = float(tilt_token)  # retro-compat anciens fichiers (1.0/2.0/3.0)
        return FitResult(
            id="mean: " + parts[0].strip(), cat1="F", cat2="i",
            nb=int(parts[2]), dec=float(parts[3]), inc=float(parts[4]), mad=float(parts[5]),
            tx=(float(parts[6]), 0.0), par3_mean=orientation,
            lat=float(parts[8]), rlong=float(parts[9]),
            par4=float(parts[10]), par5=float(parts[11]),
            liste=parts[12].strip(),
        )
    except ValueError:
        return None


def _iter_result_lines(path: str):
    """Lit les deux sections du fichier .r - bascule de la section
    "specimen" vers "mean" des que la ligne d'en-tete "#site mean results"
    est rencontree. Les lignes vides et les annotations "!..." (ajoutees a
    la main par l'utilisateur) sont silencieusement ignorees ici (mais
    preservees telles quelles par archivres, qui ne fait jamais table
    rase du fichier)."""
    if not os.path.exists(path):
        return
    in_mean_section = False
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("!"):
                continue
            if stripped.startswith("#"):
                if "site mean" in stripped.lower():
                    in_mean_section = True
                continue
            parts = stripped.split()
            res = _parse_mean_line(parts) if in_mean_section else _parse_specimen_line(parts)
            if res is not None:
                yield res


def archivres(
    res: FitResult, path: str, existing_ids: Optional[set] = None
) -> Tuple[int, set]:
    """Equivalent de `archivres`, mais insere CHAQUE resultat dans la
    bonne section plutot que de toujours ajouter en fin de fichier brut
    (demande explicite utilisateur : "can we split in two the file and
    write at the right place. In fortran I was writing at the end of the
    file") : un resultat specimen va a la fin de la section 1 (juste
    avant l'en-tete "#site mean results" s'il existe deja, en sautant les
    lignes vides qui la precedent pour ne pas casser l'espacement) ; une
    moyenne de site va toujours en fin de fichier (section 2, creee au
    besoin avec son propre en-tete). Genere un id anti-collision aleatoire
    (equivalent de la boucle `randomnum`/`read` qui reessaie jusqu'a 10000
    fois) - `existing_ids` : ensemble des `c`/`random` deja connus (evite
    de rouvrir le fichier a chaque appel ; passer None la premiere fois
    puis reutiliser l'ensemble retourne). Retourne `(c attribue,
    existing_ids mis a jour)`."""
    if existing_ids is None:
        existing_ids = {r.c for r in _iter_result_lines(path)}

    c = random.randint(0, 99999)
    tries = 0
    while c in existing_ids and tries < 10000:
        c = random.randint(0, 99999)
        tries += 1
    existing_ids.add(c)
    res.c = c

    is_mean = res.id[:5] == "mean:"
    new_line = _format_result_line(res)

    lines: List[str] = []
    if os.path.exists(path):
        with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
            lines = f.read().splitlines()

    mean_header_idx = next(
        (i for i, l in enumerate(lines) if "site mean" in l.strip().lower() and l.strip().startswith("#")),
        None,
    )

    if is_mean:
        if mean_header_idx is None:
            if lines and lines[-1].strip():
                lines.append("")
            lines.append(_MEAN_HEADER)
            lines.append(_cols_header(_MEAN_FIELDS))
        lines.append(new_line)
    else:
        if not lines:
            lines = [_SPECIMEN_HEADER, _cols_header(_SPECIMEN_FIELDS)]
        if mean_header_idx is None:
            lines.append(new_line)
        else:
            insert_idx = mean_header_idx
            while insert_idx > 0 and not lines[insert_idx - 1].strip():
                insert_idx -= 1
            lines.insert(insert_idx, new_line)

    with open(path, "w", encoding="iso-8859-1", errors="replace", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return c, existing_ids


def recompute_fit_geometry(res: FitResult, donnees) -> FitResult:
    """Complete un FitResult CHARGE DEPUIS .pmagres avec ce que ce format
    ne stocke plus, en le recalculant/le relisant depuis le specimen
    d'origine (`donnees`, ex. self.donnees) :

    1) cin/caz/dip/str_ (orientation du specimen - core+pendage) : le
       nouveau format .pmagres ne les stocke PAS ("as we have the
       orientation of the samples in the main file, they may not be
       needed in the result file"), donc un FitResult rechargee les a
       tous a 0.0 (valeurs par defaut du dataclass). BUG REEL corrige ici
       (signale par l'utilisateur : "the core correction is not correct
       for the results in the stereo... it seems that there is a -90
       applied instead of the core correction") : `corfor(x,y,z,cin=0,
       caz=0)` ne renvoie PAS l'identite, c'est une vraie rotation dont
       le residu ressemble a un decalage de 90 - `_correct_dec_inc`
       (calcul.py/stereo.py) appliquait donc une fausse correction des
       qu'un resultat charge depuis le fichier etait affiche en
       orientation 2/3 (in-situ/pendage), alors que le Zijderveld (qui
       lit cin/caz directement sur le SelectedSample vivant, jamais sur
       le FitResult) n'etait pas touche - d'ou l'ecart stereo/Zijderveld
       constate. Applique a TOUT resultat specimen (pas seulement L/P).
    2) tx/ty/tz (segment ajuste, pour le trace sur un Zijderveld) : pas
       stockes non plus ("to draw the line on the zijderveld plot we can
       redo the calculation") - recalcule uniquement pour L/P (seuls
       types avec un segment a tracer).

    Les resultats "mean:" (moyennes de site, pas un seul specimen) sont
    retournes inchanges. Si le specimen est introuvable dans `donnees`,
    les champs concernes restent a leur defaut (comportement inchange)."""
    if res.id[:5] == "mean:":
        return res
    ech = next((s for s in donnees if s.id == res.id), None)
    if ech is None:
        return res

    res = replace(res, cin=ech.cin, caz=ech.caz, dip=ech.dip, str_=ech.str_)

    if res.cat1 not in ("L", "P") or not getattr(ech, "mesures", None):
        return res
    jdeb = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == res.step_first), None)
    jfin = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == res.step_last), None)
    if jdeb is None or jfin is None or jfin < jdeb:
        return res

    anchored = res.orig == "o"
    if res.cat1 == "L":
        recomputed = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=res.numcomp)
    else:
        recomputed = fit_plane(ech, jdeb, jfin, normalize=True, numcomp=res.numcomp)
    if recomputed is None:
        return res
    return replace(res, tx=recomputed.tx, ty=recomputed.ty, tz=recomputed.tz)


# ---------------------------------------------------------------------------
# Conversion d'un ANCIEN fichier .r (positions de colonnes fixes, format
# Fortran d'origine) vers le nouveau .pmagres - demande explicite
# utilisateur ("during the convert legacy file to the new format, is it
# possible to convert the .r file with the results"), typiquement appelee
# en meme temps que la conversion .ren -> .prmag (voir convert_ren_to_r.py).
# Positions figees ici UNIQUEMENT pour lire un vieux fichier existant - le
# lecteur "en direct" (_iter_result_lines) n'utilise plus que le nouveau
# format.
# ---------------------------------------------------------------------------

_LEGACY_R_COLS = {
    "prefix": (0, 2), "c": (2, 7), "id": (9, 21),
    "cin": (23, 28), "caz": (29, 34), "dip": (35, 40), "str_": (41, 46),
    "cat1": (48, 49), "cat2": (49, 50), "orig": (53, 54), "demag": (59, 60),
    "numcomp": (65, 66), "nb": (68, 71), "dec": (72, 77), "inc": (78, 83),
    "par1": (85, 89), "par2": (91, 97), "par3": (98, 104),
    "tx1": (105, 115), "tx2": (116, 126), "ty1": (127, 137),
    "ty2": (138, 148), "tz1": (149, 159), "tz2": (160, 170),
}
_LEGACY_R_COLS_MEAN = {
    "tx1": (105, 115), "tx2": (116, 126),
    "lat": (128, 138), "rlong": (140, 150),
    "par4": (152, 158), "par5": (160, 166),
    "liste": (169, None),
}


def _parse_legacy_result_line(line: str) -> Optional[FitResult]:
    """Un champ "orig" (anchored/not, colonne 53:54) vide/blanc dans
    l'ancien fichier .r est traite comme "n" (non ancre) - demande
    explicite utilisateur ("during the import of legacy results, the
    missing information for anchored or not, means not")."""
    if not line.startswith("c:") or len(line) < _LEGACY_R_COLS["par3"][1]:
        return None

    def seg(key, cols=_LEGACY_R_COLS):
        a, b = cols[key]
        return line[a:b if b is not None else len(line)]

    try:
        prefix_kwargs = dict(
            c=int(seg("c")), id=seg("id").strip(),
            cin=float(seg("cin")), caz=float(seg("caz")),
            dip=float(seg("dip")), str_=float(seg("str_")),
            cat1=seg("cat1"), cat2=seg("cat2"),
            orig=seg("orig").strip() or "n",
            demag=seg("demag").strip(), numcomp=int(seg("numcomp")),
            nb=int(seg("nb")), dec=float(seg("dec")), inc=float(seg("inc")),
            mad=float(seg("par1")),
        )
        par2_raw, par3_raw = float(seg("par2")), float(seg("par3"))
    except ValueError:
        return None

    if prefix_kwargs["id"][:5] == "mean:":
        try:
            return FitResult(
                **prefix_kwargs, par2_mean=par2_raw, par3_mean=par3_raw,
                tx=(float(seg("tx1", _LEGACY_R_COLS_MEAN)), float(seg("tx2", _LEGACY_R_COLS_MEAN))),
                lat=float(seg("lat", _LEGACY_R_COLS_MEAN)), rlong=float(seg("rlong", _LEGACY_R_COLS_MEAN)),
                par4=float(seg("par4", _LEGACY_R_COLS_MEAN)), par5=float(seg("par5", _LEGACY_R_COLS_MEAN)),
                liste=seg("liste", _LEGACY_R_COLS_MEAN).strip(),
            )
        except ValueError:
            return None

    prefix_kwargs["step_first"] = int(round(par2_raw))
    prefix_kwargs["step_last"] = int(round(par3_raw))
    try:
        return FitResult(
            **prefix_kwargs,
            tx=(float(seg("tx1")), float(seg("tx2"))),
            ty=(float(seg("ty1")), float(seg("ty2"))),
            tz=(float(seg("tz1")), float(seg("tz2"))),
        )
    except ValueError:
        return None


def convert_legacy_results_file(old_path: str, new_path: str) -> int:
    """Convertit un ANCIEN fichier .r (colonnes fixes) vers le nouveau
    .pmagres (sections specimen/site mean, colonnes auto-descriptives).
    tx/ty/tz de l'ancien fichier sont IGNORES (pas ecrits dans le nouveau
    format - voir recompute_fit_geometry, qui les recalcule a la demande).
    Retourne le nombre de resultats convertis (0 si `old_path` n'existe
    pas ou ne contient aucune ligne reconnue - `new_path` n'est alors pas
    cree)."""
    if not os.path.exists(old_path):
        return 0

    specimen_lines, mean_lines = [], []
    with open(old_path, "r", encoding="iso-8859-1", errors="replace") as f:
        for raw in f:
            res = _parse_legacy_result_line(raw.rstrip("\n"))
            if res is None:
                continue
            if res.id[:5] == "mean:":
                mean_lines.append(_format_mean_line(res))
            else:
                specimen_lines.append(_format_specimen_line(res))

    if not specimen_lines and not mean_lines:
        return 0

    out_lines: List[str] = []
    if specimen_lines:
        out_lines += [_SPECIMEN_HEADER, _cols_header(_SPECIMEN_FIELDS), *specimen_lines]
    if mean_lines:
        if out_lines:
            out_lines.append("")
        out_lines += [_MEAN_HEADER, _cols_header(_MEAN_FIELDS), *mean_lines]

    with open(new_path, "w", encoding="iso-8859-1", errors="replace", newline="\n") as f:
        f.write("\n".join(out_lines) + "\n")
    return len(specimen_lines) + len(mean_lines)


def _load_data_results(
    path: str, pattern: str, cat1: str, numcomp: Optional[int]
) -> List[FitResult]:
    """Equivalent de la branche `carselect=="d"` (par defaut) de `selres` :
    resultats normaux (L/P/f/s...), les lignes "mean:" sont exclues -
    fidele au Fortran, ou une ligne mean ne matche NI la branche mean
    (carselect!='m') NI la branche normale (verifie `chaineres(10:15).ne.
    "mean: "`)."""
    pattern = (pattern or "*").upper().strip()
    nlen = len(pattern) if pattern != "*" else 0
    results = []
    for res in _iter_result_lines(path):
        if res.id[:5] == "mean:":
            continue
        if pattern != "*" and res.id[:nlen].upper() != pattern[:nlen]:
            continue
        if cat1 != "*" and res.cat1 != cat1:
            continue
        if numcomp is not None and res.numcomp != numcomp:
            continue
        results.append(res)
    return results


def _load_mean_results(path: str, pattern: str, iorient: int) -> List[FitResult]:
    """Equivalent de la branche `carselect=="m"` de `selres` : UNIQUEMENT les
    moyennes de site ("mean: <site>"), filtrees par `res.par3==float(iorient)`
    (le code d'orientation enregistre avec la moyenne) et par `pattern`
    (compare apres le prefixe "mean: ", ajoute automatiquement - equivalent
    de `enteteres="mean: "//entete` puis `numero=enteteres//chaine`)."""
    site = (pattern or "*").strip()
    full_pattern = "*" if site in ("", "*") else f"MEAN: {site.upper()}"
    nlen = len(full_pattern) if full_pattern != "*" else 0

    results = []
    for res in _iter_result_lines(path):
        if res.id[:5] != "mean:":
            continue
        if res.par3_mean != float(iorient):
            continue
        if full_pattern != "*" and res.id[:nlen].upper() != full_pattern[:nlen]:
            continue
        results.append(res)
    return results


def available_mean_orientations(path: str, pattern: str = "*") -> List[int]:
    """Diagnostic (PAS dans le Fortran) : orientations (`par3_mean`, meme
    convention que `self.orientation` - 1/2/3) sous lesquelles il existe au
    moins une moyenne "mean:" correspondant a `pattern`, IGNORE le filtre
    d'orientation - sert a expliquer un "0 resultat" en mode m/s de
    `selres` quand l'orientation courante ne correspond a aucune moyenne
    enregistree (`_load_mean_results` exige une correspondance EXACTE)."""
    site = (pattern or "*").strip()
    full_pattern = "*" if site in ("", "*") else f"MEAN: {site.upper()}"
    nlen = len(full_pattern) if full_pattern != "*" else 0

    found = set()
    for res in _iter_result_lines(path):
        if res.id[:5] != "mean:":
            continue
        if full_pattern != "*" and res.id[:nlen].upper() != full_pattern[:nlen]:
            continue
        if res.par3_mean == int(res.par3_mean):
            found.add(int(res.par3_mean))
    return sorted(found)


def _load_site_results(path: str, pattern: str, iorient: int) -> List[FitResult]:
    """Equivalent de la branche `carselect=="s"` de `selres` (label 444) :
    les moyennes matchees par `_load_mean_results`, PLUS pour chacune les
    resultats individuels references dans son champ `liste`
    ("codes:c1:c2:...", les `c` des resultats combines) - equivalent de
    `decodelisteres` (recherche par `c`, ajout sans autre filtre)."""
    means = _load_mean_results(path, pattern, iorient)
    if not means:
        return means

    results = list(means)
    for mean_res in means:
        wanted = set()
        for token in mean_res.liste.replace("codes:", "").split(":"):
            token = token.strip()
            if token.isdigit():
                wanted.add(int(token))
        if not wanted:
            continue
        for res in _iter_result_lines(path):
            if res.id[:5] != "mean:" and res.c in wanted:
                results.append(res)
    return results


def load_results(
    path: str,
    pattern: str = "*",
    carselect: str = "d",
    cat1: str = "*",
    numcomp: Optional[int] = None,
    iorient: int = 1,
) -> List[FitResult]:
    """Equivalent de `selres` (dataselect.f), les 3 modes `carselect` :

    - 'd' (Data, par defaut) : resultats normaux (L/P/f/s), filtres par
      `pattern` (id), `cat1`, `numcomp`. Les moyennes "mean:" sont exclues.
    - 'm' (Mean) : uniquement les moyennes de site, `pattern` = le nom du
      site SANS le prefixe "mean: " (ajoute automatiquement), filtrees par
      orientation (`iorient`, doit correspondre a `res.par3`). `cat1`/
      `numcomp` ignores (pas prompted par le Fortran en mode Mean).
    - 's' (Site = Data+Mean) : une moyenne matchee (meme filtre que 'm')
      PLUS les resultats individuels qui la composent (son champ `liste`).
    """
    if not os.path.exists(path):
        return []
    carselect = (carselect or "d").strip().lower()
    if carselect == "m":
        return _load_mean_results(path, pattern, iorient)
    if carselect == "s":
        return _load_site_results(path, pattern, iorient)
    return _load_data_results(path, pattern, cat1, numcomp)


# ---------------------------------------------------------------------------
# mdf : step de demi-desaimantation (MDF si AF, MDT si thermique), par deux
# methodes (norme simple, distance restante le long du chemin) - calcul.f:519-678
# ---------------------------------------------------------------------------

@dataclass
class MdfResult:
    id: str
    demagcod: str   # 'M.D.T.' / 'M.D.F.'
    unite: str
    xmdf_path: float   # xmdf2 : methode "distance restante sur le chemin"
    xmdf_norm: float   # xmdf1 : methode "norme simple du vecteur"
    rapport: float


_MDF_CODES = {
    "D": ("M.D.T.", "deg C"), "d": ("M.D.T.", "deg C"),
    "S": ("M.D.T.", "deg C"), "s": ("M.D.T.", "deg C"),
    "F": ("M.D.F.", "oersted"), "f": ("M.D.F.", "oersted"),
    "N": ("A.R.N.", "......."), "n": ("A.R.N.", "......."),
    "I": ("A.R.I.", "......."), "i": ("A.R.i.", "......."),
}


def compute_mdf(ech: SelectedSample) -> Optional[MdfResult]:
    """Equivalent de `mdf` : interpolation lineaire du step ou la remanence
    tombe a 50% de sa valeur initiale, par deux methodes. None pour les
    echantillons ARM/AR (dernier cod1 in N/I/n/i), comme le Fortran (qui
    n'imprime alors rien pour cet echantillon)."""
    mesures = ech.mesures
    n = len(mesures)
    if n < 2:
        return None

    steps = [m.etape for m in mesures]
    xs = [m.x for m in mesures]
    ys = [m.y for m in mesures]
    zs = [m.z for m in mesures]

    # steps[j]==steps[j-1] est possible avec des donnees paleointensite
    # (etapes R/V a la meme temperature) - non prevu par le Fortran (division
    # par zero), on saute simplement le croisement dans ce cas (mdf n'a pas
    # de sens pour ce type de protocole de toute facon).
    rint1 = [math.sqrt(xs[i] ** 2 + ys[i] ** 2 + zs[i] ** 2) for i in range(n)]
    rmax1 = rint1[0]
    xmdf1 = 0.0
    if rmax1 != 0:
        for j in range(1, n):
            if steps[j] == steps[j - 1]:
                continue
            if (rint1[j] / rmax1) < 0.5 and (rint1[j - 1] / rmax1) >= 0.5:
                amdf1 = (rint1[j] - rint1[j - 1]) / ((steps[j] - steps[j - 1]) * rmax1)
                bmdf1 = rint1[j] / rmax1 - steps[j] * amdf1
                if amdf1 != 0:
                    xmdf1 = (0.5 - bmdf1) / amdf1

    rint2 = [0.0] * n
    rint2[n - 1] = rint1[n - 1]
    for j in range(n - 2, -1, -1):
        seg = math.sqrt((xs[j] - xs[j + 1]) ** 2 + (ys[j] - ys[j + 1]) ** 2 + (zs[j] - zs[j + 1]) ** 2)
        rint2[j] = rint2[j + 1] + seg
    rmax2 = rint2[0]
    xmdf2 = 0.0
    if rmax2 != 0:
        for j in range(1, n):
            if steps[j] == steps[j - 1]:
                continue
            if (rint2[j] / rmax2) < 0.5 and (rint2[j - 1] / rmax2) >= 0.5:
                amdf2 = (rint2[j] - rint2[j - 1]) / ((steps[j] - steps[j - 1]) * rmax2)
                bmdf2 = rint2[j] / rmax2 - steps[j] * amdf2
                if amdf2 != 0:
                    xmdf2 = (0.5 - bmdf2) / amdf2

    demagcod, unite = _MDF_CODES.get(mesures[-1].cod1, ("A.R.N.", "......."))
    if demagcod in ("A.R.N.", "A.R.I.", "A.R.i."):
        return None

    rapport = xmdf2 / xmdf1 if xmdf1 != 0 else 0.0
    return MdfResult(id=ech.id, demagcod=demagcod, unite=unite,
                      xmdf_path=xmdf2, xmdf_norm=xmdf1, rapport=rapport)


def list_mdf(selected: List[SelectedSample]) -> str:
    lines = []
    for ech in selected:
        r = compute_mdf(ech)
        if r is None:
            continue
        lines.append(
            f"{r.id} {r.demagcod} norm {r.xmdf_path:7.2f} {r.unite}  "
            f"A.R.N. {r.xmdf_norm:7.2f} {r.unite}  Ratio: {r.rapport:7.2f}"
        )
    return "\n".join(lines) if lines else "(no result - ARM/AR samples excluded, or no 50% crossing)"


# ---------------------------------------------------------------------------
# mdi / mds : moyenne arithmetique/geometrique de l'intensite et de la
# susceptibilite sur toute la selection - calcul.f:683-865
# ---------------------------------------------------------------------------

@dataclass
class MeanIntensityResult:
    id6: str
    n: int
    unit: str
    arith_mean: float
    arith_sd: float
    arith_se: float
    geom_mean: float
    geom_inf: float
    geom_sup: float


@dataclass
class MeanSuscResult:
    n: int
    arith_mean: float
    arith_sd: float
    geom_mean: float
    geom_inf: float
    geom_sup: float


def _arith_geom_stats(values: List[float]) -> tuple:
    n = len(values)
    amoy = sum(values) / n
    somcar = sum((v - amoy) ** 2 for v in values)
    ecart = math.sqrt(somcar / (n - 1)) if n > 1 else 0.0
    erreur = ecart / math.sqrt(n) if n > 1 else 0.0

    logs = [math.log10(v) for v in values if v > 0]
    if len(logs) == n and n > 1:
        amoyg = sum(logs) / n
        somcarg = sum((v - amoyg) ** 2 for v in logs)
        ecartg = math.sqrt(somcarg / (n - 1))
        rgeom, rinf, rsup = 10 ** amoyg, 10 ** (amoyg - ecartg), 10 ** (amoyg + ecartg)
    else:
        rgeom = rinf = rsup = 0.0
    return amoy, ecart, erreur, rgeom, rinf, rsup


def compute_mean_intensity(selected: List[SelectedSample]) -> Optional[MeanIntensityResult]:
    """Equivalent de `mdi` : moyennes arithmetique/geometrique de l'intensite
    (norme du moment) sur TOUTES les mesures de la selection. None si la
    selection melange des normalisations masse/volume (comme le Fortran :
    "on ne peut pas melanger les choux et les carottes"), ou si <=1 mesure."""
    if not selected:
        return None
    norme0 = selected[0].norme
    if any(ech.norme != norme0 for ech in selected[1:]):
        return None

    factor = 1.0e3 if norme0 == "m" else 1.0e6
    unit = "Am2/kg" if norme0 == "m" else "A/m"
    rint = [
        math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2) * factor / (ech.vol or 1.0)
        for ech in selected for m in ech.mesures
    ]
    if len(rint) <= 1:
        return None

    amoy, ecart, erreur, rgeom, rinf, rsup = _arith_geom_stats(rint)
    return MeanIntensityResult(
        id6=selected[0].id[:6], n=len(rint), unit=unit,
        arith_mean=amoy, arith_sd=ecart, arith_se=erreur,
        geom_mean=rgeom, geom_inf=rinf, geom_sup=rsup,
    )


def compute_mean_susceptibility(selected: List[SelectedSample]) -> Optional[MeanSuscResult]:
    """Equivalent de `mds`, appele par `mdi`. None si <=2 mesures avec
    susceptibilite non nulle (comme le Fortran, `if(itot==1) return` /
    `if(itot==2) return`)."""
    suscint = [
        m.s * 1e-4 / (ech.vol or 1.0)
        for ech in selected for m in ech.mesures if m.s != 0.0
    ]
    if len(suscint) <= 2:
        return None
    amoy, ecart, _erreur, rgeom, rinf, rsup = _arith_geom_stats(suscint)
    return MeanSuscResult(n=len(suscint), arith_mean=amoy, arith_sd=ecart,
                           geom_mean=rgeom, geom_inf=rinf, geom_sup=rsup)


def format_mean_intensity(
    mi: MeanIntensityResult, ms: Optional[MeanSuscResult] = None,
    lat: float = 0.0, rlong: float = 0.0,
) -> str:
    lines = [
        f"{mi.id6}  Nb: {mi.n:4d}  Mean Arith. Int.: {mi.arith_mean:.3e}"
        f"        sd: {mi.arith_sd:.3e}   error: {mi.arith_se:.3e}",
        f"{mi.id6}  Nb: {mi.n:4d}  Mean Geom. Int.: {mi.geom_mean:.3e}"
        f"  ± sd.inf: {mi.geom_inf:.3e} ± sd.sup: {mi.geom_sup:.3e}",
    ]
    if ms is not None:
        lines.append(
            f"{mi.id6}   {lat:10.5f}  {rlong:10.5f}   "
            f"{mi.geom_mean:.3e}   {mi.geom_mean - mi.geom_inf:.3e}   {mi.geom_sup - mi.geom_mean:.3e}   "
            f"{ms.geom_mean:.3e}   {ms.geom_mean - ms.geom_inf:.3e}   {ms.geom_sup - ms.geom_mean:.3e}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# meaninc : inclinaison moyenne de McFadden & Reid (1982), estimateur du
# maximum de vraisemblance par iteration de Newton-Raphson - calcul.f:1011-1129
# ---------------------------------------------------------------------------

@dataclass
class MeanInclinationResult:
    n: int
    arith_mean: float
    biased_inc: float     # inclinaison biaisee (Arason)
    unbiased_inc: float   # inclinaison non biaisee (McFadden)
    precision: float
    alpha95: float
    alpha_pos: float
    alpha_neg: float


def compute_mean_inclination(
    selected: List[SelectedSample], orientation: int = 1
) -> Optional[MeanInclinationResult]:
    """Equivalent de `meaninc`. None si <2 mesures ou si l'iteration ne
    converge pas (equivalent implicite d'une boucle Fortran sans limite -
    ici bornee a 1000 iterations par securite)."""
    dirs = []
    for ech in selected:
        for m in ech.mesures:
            xx, yy, zz = apply_orientation(m.x, m.y, m.z, ech, orientation)
            _, _dec, inc = polere(xx, yy, zz)
            dirs.append(inc)

    k = len(dirs)
    if k < 2:
        return None

    eps = 1e-4
    conv = math.pi / 180.0
    som = sum(dirs)
    som1 = som2 = som3 = 0.0
    for d in dirs:
        dd = (90.0 - d) * conv
        som1 += math.cos(dd)
        som2 += math.sin(dd)
        som3 += dd

    amar = som / k
    par = som3 / k

    converged = False
    for _ in range(1000):
        a1 = math.cos(par) ** 2
        a2 = math.sin(par) ** 2
        ff = k * math.cos(par) + (a2 - a1) * som1 - 2 * math.sin(par) * math.cos(par) * som2
        ffp = -k * math.sin(par) + 2 * som1 * math.sin(2 * par) - 2 * som2 * math.cos(2 * par)
        if ffp == 0:
            return None
        xx = par - ff / ffp
        converged = abs(xx - par) < eps
        par = xx
        if converged:
            break
    if not converged:
        return None

    c = math.cos(par) * som1 + math.sin(par) * som2
    s = math.sin(par) * som1 - math.cos(par) * som2
    biais = 180.0 * s / (math.pi * c) if c != 0 else 0.0
    prec = (k - 1.0) / (2.0 * (k - c)) if (k - c) != 0 else 0.0
    ainc = 90.0 - math.degrees(par) + biais
    aai = 90.0 - math.degrees(par)
    resu = ((prec - 1.0) * k + 1.0) / prec if prec != 0 else 0.0
    toto = (0.05 ** (-1.0 / (k - 1.0)) - 1.0) * (k - resu) / resu if resu != 0 else 0.0
    toto = max(0.0, min(2.0, toto))
    alp95 = math.degrees(math.acos(1.0 - toto))

    return MeanInclinationResult(
        n=k, arith_mean=amar, biased_inc=aai, unbiased_inc=ainc,
        precision=prec, alpha95=alp95, alpha_pos=alp95 + biais, alpha_neg=alp95 - biais,
    )


def format_mean_inclination(r: MeanInclinationResult) -> str:
    return (
        f"arithmetic mean: inclination = {r.arith_mean:.2f}\n"
        f"biased inclination (Arason)           : {r.biased_inc:.2f}\n"
        f"unbiased inclination (McFadden)        : {r.unbiased_inc:.2f}\n"
        f"precision parameter                    : {r.precision:.2f}\n"
        f"alpha95 estimate                       : {r.alpha95:.2f}\n\n"
        f"mean inclination: {r.unbiased_inc:.2f} +/- {r.alpha95:.2f}\n"
        f"                    | {r.biased_inc:.2f} +{r.alpha_pos:.2f} -{r.alpha_neg:.2f}"
    )


# ---------------------------------------------------------------------------
# Koenigs : rapport de Koenigsberger Qn = aimantation remanente / aimantation
# induite dans un champ de reference - calcul.f:3869-3926
# ---------------------------------------------------------------------------

@dataclass
class KoenigsbergerRow:
    id: str
    etape: int
    cod1: str
    cod2: str
    mag: float
    ind_mag: float
    ratio: float


def compute_koenigsberger(selected: List[SelectedSample], valk: float) -> List[KoenigsbergerRow]:
    """Equivalent de `Koenigs` : `valk` = champ de reference en microteslas."""
    valc = valk / (4 * math.pi / 10.0)
    rows = []
    for ech in selected:
        vol = ech.vol or 1.0
        for m in ech.mesures:
            if m.s == 0.0:
                continue
            rxx = math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2)
            if ech.norme == "m":
                rxx = rxx * 1.0e3 / vol
                rind = valc * m.s * 1e-7 / vol
            else:
                rxx = rxx * 1.0e6 / vol
                rind = valc * m.s * 10.0 * 1e-5 / vol
            ratio = rxx / rind if rind != 0 else 0.0
            rows.append(KoenigsbergerRow(
                id=ech.id, etape=m.etape, cod1=m.cod1, cod2=m.cod2,
                mag=rxx, ind_mag=rind, ratio=ratio,
            ))
    return rows


def format_koenigsberger(rows: List[KoenigsbergerRow]) -> str:
    lines = ["Sample        step       Mag           K            Koenigsberger ratio"]
    for r in rows:
        lines.append(
            f"{r.id:<12s}  {r.etape:4d}{r.cod1}{r.cod2}  {r.mag:.3e}  {r.ind_mag:.3e}   {r.ratio:7.3f}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# diffmes : difference vectorielle entre mesures consecutives (etape j moins
# etape j+1) - dataselect.f:1156-1236
# ---------------------------------------------------------------------------

def list_diff_measurements(selected: List[SelectedSample], orientation: int = 1) -> str:
    """Equivalent de `diffmes`."""
    lines = []
    ij = 1
    for ech in selected:
        vol = ech.vol or 1.0
        factor = 1.0e3 if ech.norme == "m" else 1.0e6
        for j in range(len(ech.mesures) - 1):
            a, b = ech.mesures[j], ech.mesures[j + 1]
            dx, dy, dz = a.x - b.x, a.y - b.y, a.z - b.z
            xx, yy, zz = apply_orientation(dx, dy, dz, ech, orientation)
            mag, dec, inc = polere(xx, yy, zz)
            mag = mag * factor / vol
            lines.append(
                f"{ij:4d}: {ech.id:<12s}  {a.etape:4d}{a.cod1}{a.cod2}  "
                f"{mag:10.3e} {dec:6.1f} {inc:6.1f}  q={a.q:<4d} {a.ins:<2s} s={a.s:.1f}"
            )
            ij += 1
    return "\n".join(lines) if lines else "(aucune mesure)"


# ---------------------------------------------------------------------------
# viscos : indice de viscosite sur des paires N+/N- - viscos.f
# MUTE selected en place (comme le Fortran).
# ---------------------------------------------------------------------------

def apply_viscosity_test(selected: List[SelectedSample]) -> List[str]:
    """Equivalent de `viscos` : pour chaque echantillon dont les deux
    premieres mesures sont codees N+/N-, remplace la premiere par le
    vecteur MOYEN et la seconde par le vecteur DIFFERENCE (cod2 mis a
    '0'/'v'), .q de la seconde recevant l'indice de viscosite (%). Retourne
    les avertissements (echantillons dont les 2 premieres mesures ne sont
    pas N+/N-, comme les 2 messages `write` du Fortran)."""
    warnings = []
    for ech in selected:
        if len(ech.mesures) < 2:
            continue
        m0, m1 = ech.mesures[0], ech.mesures[1]
        if not (m0.cod1 == "N" and m1.cod1 == "N"):
            warnings.append(f"{ech.id} : attention code different de N et N")
            continue
        if not (m0.cod2 == "+" and m1.cod2 == "-"):
            warnings.append(f"{ech.id} : attention code different de + et -")
            continue

        xm, ym, zm = (m0.x + m1.x) / 2.0, (m0.y + m1.y) / 2.0, (m0.z + m1.z) / 2.0
        xd, yd, zd = (m0.x - m1.x) / 2.0, (m0.y - m1.y) / 2.0, (m0.z - m1.z) / 2.0
        iq1 = (m0.q + m1.q) // 2
        num = math.sqrt(xd ** 2 + yd ** 2 + zd ** 2)
        den = math.sqrt(xm ** 2 + ym ** 2 + zm ** 2)
        visco = 100.0 * num / den if den != 0 else 0.0

        m0.x, m0.y, m0.z = xm, ym, zm
        m0.q = iq1
        m0.cod2 = "0"
        m1.x, m1.y, m1.z = xd, yd, zd
        m1.q = int(visco)
        m1.cod2 = "v"
    return warnings


# ---------------------------------------------------------------------------
# soustra : soustrait le vecteur brut d'une mesure de toutes les autres de
# l'echantillon, puis supprime cette ligne - dataselect.f:1239-1283
# MUTE ech.mesures en place.
# ---------------------------------------------------------------------------

def apply_subtraction(ech: SelectedSample, row_index: int) -> Optional[str]:
    """Equivalent de `soustra`. `row_index` : 1-based (comme les autres
    prompts "numero de ligne"). Retourne un message d'erreur si invalide,
    sinon None."""
    idx = row_index - 1
    if not (0 <= idx < len(ech.mesures)):
        return "invalid line number"
    xs, ys, zs = ech.mesures[idx].x, ech.mesures[idx].y, ech.mesures[idx].z
    for m in ech.mesures:
        m.x -= xs
        m.y -= ys
        m.z -= zs
    del ech.mesures[idx]
    return None


# ---------------------------------------------------------------------------
# holderarm : enregistre les 6 mesures ARM d'un porte-echantillon vide
# (x+,x-,y+,y-,z+,z-), pour une future correction dans `anisot` (non porte) -
# calcul.f:2337-2414
# ---------------------------------------------------------------------------

@dataclass
class ArmHolderBackground:
    x: Tuple[float, float, float, float, float, float]
    y: Tuple[float, float, float, float, float, float]
    z: Tuple[float, float, float, float, float, float]


def record_arm_holder(
    ech: SelectedSample, ixp: int, ixm: int, iyp: int, iym: int, izp: int, izm: int
) -> Optional[ArmHolderBackground]:
    """Equivalent de `holderarm`. Indices 1-based dans ech.mesures. None si
    un indice est invalide."""
    idx = [ixp, ixm, iyp, iym, izp, izm]
    if any(not (1 <= i <= len(ech.mesures)) for i in idx):
        return None
    xs = tuple(ech.mesures[i - 1].x for i in idx)
    ys = tuple(ech.mesures[i - 1].y for i in idx)
    zs = tuple(ech.mesures[i - 1].z for i in idx)
    return ArmHolderBackground(x=xs, y=ys, z=zs)


# ---------------------------------------------------------------------------
# vitref : correction de vitesse de refroidissement (TRM rapide labo vs TRM
# lente en four), a partir d'un motif de mesures L/Q (etape lente) precede
# de R/V (etapes rapides, meme cod2) - calcul.f:3461-3625 (branche "live",
# le reste apres `go to 555` est mort/inutilise cote Fortran).
# ---------------------------------------------------------------------------

@dataclass
class CoolingRateResult:
    id: str
    correction_avant: float
    correction_apres: float
    correction_mean: float
    evolution_pct: float
    evolution_lente_avant_pct: float
    evolution_lente_apres_pct: float


def detect_cooling_rate_rows(ech: SelectedSample) -> Optional[Tuple[int, int, int, int, int]]:
    """Equivalent du mode automatise de `vitref` : cherche un pas 'L' suivi
    d'un pas 'Q', puis dans les 8 mesures precedentes les pas 'R' et 'V' de
    MEME cod2. Retourne (irr,irv,ilr,iravant,irapres) en indices 1-based
    (memes conventions que les autres prompts "numero de ligne"), ou None."""
    mesures = ech.mesures
    for j in range(len(mesures) - 1):
        if mesures[j].cod1 == "L" and mesures[j + 1].cod1 == "Q":
            ilr, irapres = j, j + 1
            irr = irv = iravant = None
            for ki in range(1, 9):
                idx = j - ki
                if idx < 0:
                    break
                if mesures[idx].cod1 == "R" and mesures[idx].cod2 == mesures[j].cod2:
                    irr = idx
                if mesures[idx].cod1 == "V" and mesures[idx].cod2 == mesures[j].cod2:
                    irv = idx
                    iravant = idx
            if irr is not None and irv is not None:
                return irr + 1, irv + 1, ilr + 1, iravant + 1, irapres + 1
    return None


def compute_cooling_rate(
    ech: SelectedSample, irr: int, irv: int, ilr: int, iravant: int, irapres: int
) -> CoolingRateResult:
    """Equivalent du calcul de `vitref` (partie "live" avant le `go to 555`
    mort). Indices 1-based."""
    m = ech.mesures
    vol = ech.vol or 1.0

    def mag(x, y, z):
        mg, _dec, _inc = polere(x, y, z)
        return mg * 1.0e6 / vol

    r, v, l, av, ap = m[irr - 1], m[irv - 1], m[ilr - 1], m[iravant - 1], m[irapres - 1]

    arnx, arny, arnz = (r.x + v.x) / 2.0, (r.y + v.y) / 2.0, (r.z + v.z) / 2.0
    atr_rapide_mag = mag((r.x - v.x) / 2.0, (r.y - v.y) / 2.0, (r.z - v.z) / 2.0)  # noqa: F841 (echo Fortran)
    atrlente = mag(l.x - arnx, l.y - arny, l.z - arnz)
    atrrapavant = mag(av.x - arnx, av.y - arny, av.z - arnz)
    atrrapapres = mag(ap.x - arnx, ap.y - arny, ap.z - arnz)

    corr_avant = atrrapavant / atrlente if atrlente else 0.0
    corr_apres = atrrapapres / atrlente if atrlente else 0.0
    evol = ((atrrapapres / atrrapavant) - 1.0) * 100.0 if atrrapavant else 0.0
    evol_la = ((atrlente / atrrapavant) - 1.0) * 100.0 if atrrapavant else 0.0
    evol_lp = ((atrlente / atrrapapres) - 1.0) * 100.0 if atrrapapres else 0.0

    return CoolingRateResult(
        id=ech.id, correction_avant=corr_avant, correction_apres=corr_apres,
        correction_mean=(corr_avant + corr_apres) / 2.0,
        evolution_pct=evol, evolution_lente_avant_pct=evol_la, evolution_lente_apres_pct=evol_lp,
    )


def format_cooling_rate(r: CoolingRateResult) -> str:
    return (
        f"Sample: {r.id}  % Evolution: {r.evolution_pct:.1f}%"
        f"  slow/Fast before: {r.evolution_lente_avant_pct:.1f}%"
        f"  slow/Fast after: {r.evolution_lente_apres_pct:.1f}%\n"
        f"Sample: {r.id}  cooling correction (1st fast): {r.correction_avant:.3f}\n"
        f"Sample: {r.id}  cooling correction (last fast): {r.correction_apres:.3f}\n"
        f"Sample: {r.id}  mean cooling correction: {r.correction_mean:.3f}"
    )


# ---------------------------------------------------------------------------
# anisot / anisoauto : calcul du tenseur d'anisotropie (ARM/TRM 6 positions)
# - calcul.f:2418-3310 (anisot, saisie manuelle) / 3951-4900+ (anisoauto,
# detection automatique). Ne porte QUE la variante 'A0' (tenseur de base,
# calcule directement sur les 6 mesures reelles) - PAS les 14 variantes
# jackknife (A+/A-/A1-A6/B1-B6, reconstruction de chaque position tour a
# tour pour evaluer la robustesse) ni le bloc de correction d'alteration/
# evolution (testevol) : deferrees a une etape ulterieure (demande
# explicite de l'utilisateur, "Etape 1").
#
# VERIFIE OCTET-PRES contre un vrai fichier .ANI fourni par l'utilisateur
# (Miriam_2025b/SanJuan_Pmag.ANI, echantillon "02A", etape 460) : calcul a
# la main des 6 valeurs k11/k22/k33/k12/k23/k13 a partir des mesures
# reelles du .ren (460X+/460X-/460Y+/460Y-/460RH/460VH) - identiques a la
# ligne 'A0' reelle du .ANI (aucune ligne de base porte-echantillon dans
# ce cas, holder=None).
# ---------------------------------------------------------------------------

@dataclass
class AniTensor:
    id: str
    code2: str  # 'A0'=TRM, 'F0'=ARM, 'N0'=susceptibilite
    k11: float
    k22: float
    k33: float
    k12: float
    k23: float
    k13: float


_ANI_CODE2 = {1: "A0", 2: "F0", 3: "N0"}
_ANI_POSITION_KEYS = ("X+", "X-", "Y+", "Y-", "Z+", "Z-")


def detect_six_positions(
    ech: "SelectedSample", force_zplus_label: Optional[str] = None,
) -> Optional[Dict[str, Measurement]]:
    """Equivalent de la detection automatique des 6 positions dans
    anisoauto (calcul.f:4005-4059) : trouve l'etape des mesures 'X'
    (cod1='X', X+/X-), puis a cette MEME etape les mesures codees 'R' et
    'V' (cod1='R'/'V') servent de substituts pour Z+/Z- quand
    l'echantillon n'a pas de vraies mesures Z+/Z- - le "label" de
    substitution est cod1+cod2 de la mesure trouvee (ex. 'RH'/'VH',
    verifie sur donnees reelles). `force_zplus_label` : remplace le label
    Z+ detecte (ex. 'ZB' - voir _check_trm_evolution) plutot que le
    substitut R normal, pour la reclassification apres detection d'une
    evolution/alteration significative (calcul.f:4107 `Zplus="ZB"`).
    Retourne None si l'etape X est absente ou si les 6 positions ne sont
    pas toutes identifiees (equivalent "decoding incomplete" du Fortran)."""
    item_rv = None
    for m in ech.mesures:
        if m.cod1 == "X":
            item_rv = m.etape
    if item_rv is None:
        return None

    zplus_label = zminus_label = None
    for m in ech.mesures:
        if m.etape != item_rv:
            continue
        if m.cod1 == "R":
            zplus_label = m.cod1 + m.cod2
        elif m.cod1 == "V":
            zminus_label = m.cod1 + m.cod2
    if force_zplus_label:
        zplus_label = force_zplus_label

    positions: Dict[str, Measurement] = {}
    for m in ech.mesures:
        testcode = m.cod1 + m.cod2
        if zplus_label and testcode == zplus_label:
            testcode = "Z+"
        elif zminus_label and testcode == zminus_label:
            testcode = "Z-"
        if testcode in _ANI_POSITION_KEYS:
            positions[testcode] = m

    if len(positions) != 6:
        return None
    return positions


def _position_vector(
    positions: Dict[str, Measurement], key: str, idx: int,
    holder: Optional[ArmHolderBackground],
) -> Tuple[float, float, float]:
    """Composantes (x,y,z) de la position `key` (idx 0..5 dans l'ordre
    X+,X-,Y+,Y-,Z+,Z- - meme ordre que holder.x/y/z), ligne de base
    porte-echantillon deja soustraite si fournie."""
    m = positions[key]
    x, y, z = m.x, m.y, m.z
    if holder is not None:
        x -= holder.x[idx]
        y -= holder.y[idx]
        z -= holder.z[idx]
    return x, y, z


@dataclass
class AnisotropyPositionDiag:
    """Une ligne de diagnostic 'TRM values' (calcul.f:4340-4477) : mesure
    d'une position apres soustraction de la NRM residuelle moyenne
    (nrm_mean), en unites physiques (A/m ou Am2/kg) + dec/inc."""
    key: str
    measurement: Measurement
    intensity: float
    dec: float
    inc: float


def _nrm_mean_and_diag(
    positions: Dict[str, Measurement], holder: Optional[ArmHolderBackground],
    norme: str, vol: float,
) -> Tuple[Tuple[float, float, float], List[AnisotropyPositionDiag], float, Dict[str, Tuple[float, float, float]]]:
    """Equivalent de calcul.f:4206-4479 : moyenne des 3 paires (X+/X-,
    Y+/Y-, Z+/Z-) -> `nrm_mean` (la "NRM residuelle" commune aux 6
    positions - une anisotropie basee sur des TRM PARTIELLES laisse
    souvent une part de NRM non remplacee, qui s'annule dans le tenseur
    final par difference mais fausse visuellement chaque position prise
    isolement - demande explicite de l'utilisateur). Retourne (nrm_mean,
    liste de diagnostics par position, deviation_pct (`deviatTRM`),
    composantes brutes NRM-soustraites par position - utilisees par
    _check_position_inversion)."""
    idx_of = {"X+": 0, "X-": 1, "Y+": 2, "Y-": 3, "Z+": 4, "Z-": 5}
    vecs = {k: _position_vector(positions, k, idx_of[k], holder) for k in _ANI_POSITION_KEYS}

    pair_mean_x = tuple((vecs["X+"][c] + vecs["X-"][c]) / 2 for c in range(3))
    pair_mean_y = tuple((vecs["Y+"][c] + vecs["Y-"][c]) / 2 for c in range(3))
    pair_mean_z = tuple((vecs["Z+"][c] + vecs["Z-"][c]) / 2 for c in range(3))
    nrm_mean = tuple((pair_mean_x[c] + pair_mean_y[c] + pair_mean_z[c]) / 3 for c in range(3))

    def scale(intensity: float) -> float:
        return intensity * 1.0e3 / vol if norme == "m" else intensity * 1.0e6 / vol

    diags: List[AnisotropyPositionDiag] = []
    raw_components: Dict[str, Tuple[float, float, float]] = {}
    sommetrm = 0.0
    trmmean = 0.0
    prev_raw = None
    for key in _ANI_POSITION_KEYS:
        raw = tuple(vecs[key][c] - nrm_mean[c] for c in range(3))
        raw_components[key] = raw
        mag, dec, inc = polere(*raw)
        trmmean += mag
        diags.append(AnisotropyPositionDiag(
            key=key, measurement=positions[key], intensity=scale(mag), dec=dec, inc=inc,
        ))
        if key in ("X-", "Y-", "Z-"):
            sommetrm += math.sqrt(sum((raw[c] + prev_raw[c]) ** 2 for c in range(3)))
        prev_raw = raw
    sommetrm /= 3.0
    trmmean /= 6.0
    deviation_pct = 100.0 * sommetrm / trmmean if trmmean else 0.0

    return nrm_mean, diags, deviation_pct, raw_components


def _check_position_inversion(raw: Dict[str, Tuple[float, float, float]]) -> Optional[str]:
    """Equivalent de calcul.f:4386-4393 (X) et 4420-4427 (Y) : detecte une
    inversion de la position lors de l'acquisition TRM - le "+"-composant
    (NRM-soustrait) d'une position "+" doit etre positif, celui de la
    position "-" negatif (physiquement, ce sont des champs opposes).
    Fidele au Fortran, y compris son ASYMETRIE : le test X compare le
    composant BRUT de X+ a la MAGNITUDE (toujours positive) de X-,
    reduisant de facto la condition a "composant X+ negatif" ; le test Y
    compare les deux composants BRUTS (signes) comme attendu. Ce n'est
    probablement pas intentionnel cote Fortran, mais reproduit tel quel
    par fidelite - PAS de test Z equivalent dans le Fortran d'origine.
    Retourne 'X', 'Y' ou None (aucune inversion detectee)."""
    x_plus_raw = raw["X+"][0]
    x_minus_mag = polere(*raw["X-"])[0]
    if x_plus_raw < 0 and x_minus_mag > 0:
        return "X"
    y_plus_raw = raw["Y+"][1]
    y_minus_raw = raw["Y-"][1]
    if y_plus_raw < 0 and y_minus_raw > 0:
        return "Y"
    return None


def _find_zb(ech: "SelectedSample", item_rv: int) -> Optional[Measurement]:
    """Equivalent de calcul.f:4014 : mesure cod1='Z',cod2='B' a la meme
    etape que les mesures X (item_rv) - une repetition de la position Z,
    utilisee comme controle de type pTRM-check (verifie qu'aucune
    alteration/evolution n'a eu lieu entre les mesures) - PAS un
    equivalent MagIC standard, specifique a ce protocole. None si
    absente."""
    for m in ech.mesures:
        if m.etape == item_rv and m.cod1 == "Z" and m.cod2 == "B":
            return m
    return None


def _check_trm_evolution(
    positions: Dict[str, Measurement], zb: Measurement, holder: Optional[ArmHolderBackground],
) -> float:
    """Equivalent de calcul.f:4095-4108 : compare la mesure de controle ZB
    a la moyenne Z+/Z- (arnz) et a la demi-difference Z+/Z- normale
    (atrz) - `trmevol = atrzb.z / atrz.z` proche de 1.0 si ZB est
    coherente avec la paire Z+/Z- normale (pas d'alteration). Retourne
    trmevol (ratio, PAS encore convertit en %) ; 1.0 si atrz.z est nul
    (evite une division par zero, absente du Fortran d'origine)."""
    zpx, zpy, zpz = _position_vector(positions, "Z+", 4, holder)
    zmx, zmy, zmz = _position_vector(positions, "Z-", 5, holder)
    arnz_z = (zpz + zmz) / 2
    atrz_z = (zpz - zmz) / 2
    zb_x, zb_y, zb_z = zb.x, zb.y, zb.z
    if holder is not None:
        zb_z -= holder.z[4]  # meme position "Z+" pour la ligne de base (ZB est une repetition de Z+)
    atrzb_z = zb_z - arnz_z
    return atrzb_z / atrz_z if atrz_z else 1.0


@dataclass
class AnisotropyComputation:
    """Trace complete du calcul (pas seulement le resultat final) - pour
    pouvoir afficher toute la sequence comme le fait le Fortran (console
    verbeuse), demande explicite de l'utilisateur : le tenseur BRUT avant
    symetrisation (`raw`) et les diagnostics par position sont le seul
    moyen de verifier qu'un echantillon etait bien oriente lors de
    l'acquisition de la TRM."""
    tensor: AniTensor  # 'A0', symetrise
    raw: Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]
    positions: Dict[str, Measurement]
    holder_used: bool
    nrm_mean: Tuple[float, float, float]
    position_diags: List[AnisotropyPositionDiag]
    deviation_pct: float  # `deviatTRM`
    swapped_axes: List[str]  # ex. ['X'] si une inversion a ete corrigee
    trm_evolution_pct: Optional[float]  # `TRMevol`, None si pas de ZB trouvee
    zb_used: bool  # Zplus remplace par ZB (evolution > seuil + use_zb_on_evolution)


def compute_anisotropy_tensor(
    ech: "SelectedSample", holder: Optional[ArmHolderBackground] = None,
    use_zb_on_evolution: bool = False,
    positions: Optional[Dict[str, Measurement]] = None,
) -> Optional[AnisotropyComputation]:
    """Tenseur 'A0' (equivalent de la branche de base d'anisot/anisoauto,
    calcul.f:4005-4595) : detecte les 6 positions (detect_six_positions) -
    ou reprend celles fournies via `positions` (saisie manuelle, "n" a
    "automated recognition of the 6 steps ?", calcul.f:2450-2459 - meme
    calcul ensuite, seule l'identification des positions differe),
    soustrait la ligne de base porte-echantillon si fournie, calcule la
    NRM residuelle moyenne et les diagnostics par position (une
    anisotropie basee sur des TRM PARTIELLES laisse souvent une part de
    NRM non remplacee - demande explicite de l'utilisateur), detecte et
    corrige une inversion X+/X- ou Y+/Y- (specimen mal oriente pendant
    l'acquisition - calcul.f:4386-4427), verifie une eventuelle
    alteration via une mesure de controle ZB si presente
    (_check_trm_evolution ; substitue Z+ par ZB si l'evolution depasse
    5% ET `use_zb_on_evolution`), puis calcule le tenseur BRUT
    (antisymetrique, `raw`) et le tenseur symetrique final. Retourne None
    si les 6 positions ne sont pas identifiees.

    `use_zb_on_evolution` correspond au choix demande UNE FOIS pour tout
    le lot dans le Fortran ("using ZB instead of R for a > + or - 5%
    evolution ? y/N"), pas par echantillon."""
    if positions is None:
        positions = detect_six_positions(ech)
        if positions is None:
            return None
    else:
        positions = dict(positions)  # copie - peut etre mutee par le swap X/Y ci-dessous

    swapped_axes: List[str] = []
    for _ in range(4):  # borne de securite - converge en 1-2 iterations en pratique
        nrm_mean, diags, deviation_pct, raw_components = _nrm_mean_and_diag(
            positions, holder, ech.norme, ech.vol)
        axis = _check_position_inversion(raw_components)
        if axis is None:
            break
        if axis == "X":
            positions["X+"], positions["X-"] = positions["X-"], positions["X+"]
        else:
            positions["Y+"], positions["Y-"] = positions["Y-"], positions["Y+"]
        swapped_axes.append(axis)

    trm_evolution_pct: Optional[float] = None
    zb_used = False
    item_rv = positions["X+"].etape
    zb = _find_zb(ech, item_rv)
    if zb is not None:
        trmevol = _check_trm_evolution(positions, zb, holder)
        trm_evolution_pct = (trmevol - 1.0) * 100.0
        if use_zb_on_evolution and (trmevol > 1.05 or trmevol < 0.95):
            reclassified = detect_six_positions(ech, force_zplus_label=zb.cod1 + zb.cod2)
            if reclassified is not None:
                positions = reclassified
                zb_used = True
                nrm_mean, diags, deviation_pct, raw_components = _nrm_mean_and_diag(
                    positions, holder, ech.norme, ech.vol)

    xpx, xpy, xpz = _position_vector(positions, "X+", 0, holder)
    xmx, xmy, xmz = _position_vector(positions, "X-", 1, holder)
    ypx, ypy, ypz = _position_vector(positions, "Y+", 2, holder)
    ymx, ymy, ymz = _position_vector(positions, "Y-", 3, holder)
    zpx, zpy, zpz = _position_vector(positions, "Z+", 4, holder)
    zmx, zmy, zmz = _position_vector(positions, "Z-", 5, holder)

    ani011, ani012, ani013 = (xpx - xmx) / 2, (xpy - xmy) / 2, (xpz - xmz) / 2
    ani021, ani022, ani023 = (ypx - ymx) / 2, (ypy - ymy) / 2, (ypz - ymz) / 2
    ani031, ani032, ani033 = (zpx - zmx) / 2, (zpy - zmy) / 2, (zpz - zmz) / 2

    tensor = AniTensor(
        id=ech.id, code2="A0",
        k11=ani011, k22=ani022, k33=ani033,
        k12=(ani012 + ani021) / 2,
        k23=(ani032 + ani023) / 2,
        k13=(ani013 + ani031) / 2,
    )
    raw = ((ani011, ani012, ani013), (ani021, ani022, ani023), (ani031, ani032, ani033))
    return AnisotropyComputation(
        tensor=tensor, raw=raw, positions=positions, holder_used=holder is not None,
        nrm_mean=nrm_mean, position_diags=diags, deviation_pct=deviation_pct,
        swapped_axes=swapped_axes, trm_evolution_pct=trm_evolution_pct, zb_used=zb_used,
    )


def write_ani_tensor(
    path: str, ech: "SelectedSample", tensor: AniTensor,
    positions: Optional[Dict[str, Measurement]] = None,
    trm_evolution_pct: float = 0.0, deviation_pct: float = 0.0,
) -> None:
    """Ajoute une ligne au fichier .ANI (calcul.f:4592-4595/4565-4566,
    format list-directed Fortran - separateurs espace, relu par
    `read_ani_tensor` via `line.split()`). L'etape et les labels Z+/Z-
    (pour le texte informatif) viennent de `positions` (voir
    AnisotropyComputation.positions - passer ce qui a deja ete detecte
    evite un second appel a detect_six_positions) ; si non fourni,
    redetecte depuis `ech`. `trm_evolution_pct`/`deviation_pct` : voir
    AnisotropyComputation (0.0 par defaut si absents, ex. pas de mesure
    ZB trouvee)."""
    if positions is None:
        positions = detect_six_positions(ech)
    etape = positions["X+"].etape if positions else 0
    zplus_label = positions["Z+"].cod1 + positions["Z+"].cod2 if positions else ""
    zminus_label = positions["Z-"].cod1 + positions["Z-"].cod2 if positions else ""
    info = (
        f'"TRM evo: {trm_evolution_pct:5.1f} deviation:{deviation_pct:5.1f}'
        f'  steps: X+ X- Y+ Y- {zplus_label} {zminus_label}"'
    )
    with open(path, "a", encoding="utf-8") as f:
        f.write(
            f"D  {ech.id} {ech.cin} {ech.caz} {ech.dip} {ech.str_} {etape} "
            f"{tensor.code2} {tensor.k11:.6E} {tensor.k22:.6E} {tensor.k33:.6E} "
            f"{tensor.k12:.6E} {tensor.k23:.6E} {tensor.k13:.6E} 1.00000 {info}\n"
        )


# ---------------------------------------------------------------------------
# inverseani / inverse : applique l'inverse d'un tenseur d'anisotropie
# (lu dans un fichier .ANI) a toutes les mesures d'un echantillon -
# plotpaleoint2.f:1803-1918 + calcul.f:3311-3327 (`inverse`)
# ---------------------------------------------------------------------------

def inverse_symmetric_3x3(
    a: float, b: float, c: float, d: float, e: float, f: float,
    x: float, y: float, z: float,
) -> Tuple[float, float, float]:
    """Equivalent de `inverse` : applique l'inverse de la matrice symetrique
    [[a,b,c],[b,d,e],[c,e,f]] au vecteur (x,y,z). Retourne (x,y,z) inchange
    si le determinant est nul (matrice non inversible)."""
    det = a * d * f - a * e * e - b * b * f + 2 * b * c * e - c * c * d
    if det == 0:
        return x, y, z
    a11, a12, a13 = (d * f - e * e) / det, -(b * f - c * e) / det, (b * e - c * d) / det
    a22, a23 = (a * f - c * c) / det, -(a * e - b * c) / det
    a33 = (a * d - b * b) / det
    x1 = a11 * x + a12 * y + a13 * z
    y1 = a12 * x + a22 * y + a23 * z
    z1 = a13 * x + a23 * y + a33 * z
    return x1, y1, z1


def compute_anicor_factor(
    tensor: AniTensor, direction_specimen_frame: Tuple[float, float, float],
) -> float:
    """Equivalent de `anicor` (plotpaleoint2.f:1686-1785) : facteur de
    correction d'anisotropie `fcor` a appliquer a l'intensite (H*fcor)
    d'apres le tenseur TRM 'A0' et la direction ANCREE (repere specimen,
    PAS encore tournee selon l'orientation - `xnrm_nocor` du Fortran).

    NON PARFAITEMENT VERIFIE : sur l'echantillon reel 06A (tenseur .ANI
    reel + direction ancree calculee), fcor=0.892 est obtenu contre 0.887
    dans le transcript Fortran reel (~0.5% d'ecart) - vraisemblablement
    une petite difference de precision numerique dans le calcul du
    vecteur propre principal (ACP) entre numpy et l'algorithme Fortran
    d'origine, la formule elle-meme est transcrite telle quelle depuis le
    source. A surveiller si des resultats bien plus divergents sont
    observes sur d'autres echantillons."""
    rnorm = (tensor.k11 + tensor.k22 + tensor.k33) / 3.0
    if rnorm == 0:
        return 0.0
    a, b, c = tensor.k11 / rnorm, tensor.k12 / rnorm, tensor.k13 / rnorm
    d, e, f = tensor.k22 / rnorm, tensor.k23 / rnorm, tensor.k33 / rnorm
    xx, yy, zz = direction_specimen_frame
    x1, y1, z1 = inverse_symmetric_3x3(a, b, c, d, e, f, xx, yy, zz)
    norm1 = math.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    if norm1 == 0:
        return 0.0
    x1, y1, z1 = x1 / norm1, y1 / norm1, z1 / norm1
    rhlab = math.sqrt(c * c + e * e + f * f)
    x2 = a * x1 + b * y1 + c * z1
    y2 = d * y1 + e * z1
    z2 = f * z1
    ranc = math.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    return rhlab / ranc if ranc else 0.0


def read_ani_tensor(path: str, sample_id: str, code2: str) -> Optional[AniTensor]:
    """Equivalent de la recherche dans le fichier .ANI de `inverseani`
    (lecture list-directed Fortran : champs separes par des espaces -
    truc,id,cin,caz,dip,str,etape,code2,k11,k22,k33,k12,k23,k13,s, 15 champs).
    VERIFIE contre un vrai fichier .ANI (Miriam_2025b/SanJuan_Pmag.ANI) -
    parsing confirme correct sur des lignes reelles produites par
    anisoauto (voir compute_anisotropy_tensor)."""
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="iso-8859-1", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 14:
                continue
            if parts[1] != sample_id or parts[7] != code2:
                continue
            try:
                k11, k22, k33, k12, k23, k13 = (float(p) for p in parts[8:14])
            except ValueError:
                continue
            return AniTensor(id=parts[1], code2=parts[7], k11=k11, k22=k22, k33=k33,
                              k12=k12, k23=k23, k13=k13)
    return None


def apply_inverse_anisotropy(ech: SelectedSample, tensor: AniTensor) -> None:
    """Equivalent de `inverseani` : normalise le tenseur par sa trace/3,
    puis applique l'inverse a CHAQUE mesure de `ech` - MUTE ech.mesures en
    place."""
    rnorm = (tensor.k11 + tensor.k22 + tensor.k33) / 3.0
    if rnorm == 0:
        return
    a, b, c = tensor.k11 / rnorm, tensor.k12 / rnorm, tensor.k13 / rnorm
    d, e, f = tensor.k22 / rnorm, tensor.k23 / rnorm, tensor.k33 / rnorm
    for m in ech.mesures:
        m.x, m.y, m.z = inverse_symmetric_3x3(a, b, c, d, e, f, m.x, m.y, m.z)
