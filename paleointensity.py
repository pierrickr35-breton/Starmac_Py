"""
Diagramme d'Arai (paleointensite, protocole Thellier a double chauffe) :
port de la sous-routine `paleoin` (plotpaleoint2.f, Starmac_AWE_v22), limite
au coeur du diagramme NRM/TRM (points, checks pTRM, droite ajustee, bloc de
texte NRM/Hlab/H/Q) - PAS le panneau Zijderveld/stereo modifie ni les
corrections anisotropie/vitesse de refroidissement (hors perimetre de cette
premiere passe, portees plus tard si besoin).

Deux conventions de codes portees, auto-detectees par
`detect_method_and_hlab()` a partir des cod1 presents dans les mesures du
specimen (demande explicite utilisateur : "the method is IZZI if there
are S codes and Thellier if there are R and V codes") :

- 'THELLIER' (bloc `reponse==1` "Thellier ordering", lignes ~464-629) -
  verifiee empiriquement (echantillon 10CL1406A de newFormat_file.ren,
  codes N/R/V/P reels) :
  - cod1='N' : NRM initiale (etape 0, avant toute chauffe) - sert
    uniquement de normalisateur `arno` (magnitude de la toute premiere
    mesure).
  - cod1='R' : etape a champ nul (remanence apres chauffe, avant
    application du champ labo) - alimente `yp` (NRM restant, axe
    vertical).
  - cod1='V' : etape en champ (apres application du champ labo,
    acquisition de pTRM) - alimente `xp` (TRM acquis, axe horizontal),
    par difference vectorielle avec le 'R' de MEME temperature.
  - cod1='P' : mesure de verification pTRM.
- 'IZZI' (bloc "coe version" reponse==2, lignes 292-461, convention
  N/S/R/P - S=zero-field, R=in-field) - portee via `_compute_arai_coe`,
  testee sur des donnees IZZI reelles (import MagIC, ex. specimen
  HP01-01) mais PAS verifiee octet-pres contre un transcript Fortran
  (aucun n'existe pour ce protocole) - voir docstring de
  `_compute_arai_coe` pour le detail des ecarts de convention avec le
  bloc THELLIER (pas de /2, signe de winc oppose, filtre cleanpaleo
  partiel).

LIMITE REELLE DECOUVERTE EN TESTANT SUR DONNEES REELLES (specimen
HP01-01, qui porte une experience ATRM AVANT l'IZZI dans la meme liste
de mesures - voir convert_magic_to_r.py) : `arno` (normalisateur) est
TOUJOURS calcule depuis `mesures[0]` (equivalent Fortran `x(ideb)` avec
ideb=1 fixe) - si le premier pas du specimen appartient a une AUTRE
experience (ex. ATRM) plutot qu'au N initial du protocole de
paleointensite lui-meme, `arno` (et donc TOUTE la normalisation xp/yp)
est fausse. Ce n'est pas un bug de portage : le Fortran d'origine a la
MEME hypothese (un specimen = une seule experience de paleointensite) -
mais elle ne tient plus pour un specimen MagIC combinant plusieurs
experiences. PAS corrige ici (changerait un comportement fidele au
Fortran sans instruction explicite) - a traiter si besoin, probablement
en filtrant `ech.mesures` sur le sous-ensemble pertinent avant l'appel a
compute_arai plutot qu'en modifiant l'algorithme lui-meme.

API : `compute_arai()` fait le travail de `paleoin` equivalent aux lignes
464-629 (THELLIER) ou 292-461 (IZZI), `fit_arai_line()` l'ajustement
(lignes 1153-1248), `draw_arai()` le trace (lignes 819-1331) sur un `ctx`
compatible PlotContext/SVGWriter (memes conventions que zijderveld.py/
stereo.py).
"""

import math
from dataclasses import dataclass, field, replace as _dc_replace
from typing import List, Optional, Tuple

import numpy as np
from matplotlib.figure import Figure

from selection import SelectedSample, Measurement, apply_orientation, polere
from plotlib import PlotContext
from calcul import linear_fit
from zijderveld import draw_zijderveld
from stereo import draw_stereo_net, superc

# cartest() : deux alphabets possibles pour decoder cod2 en indice de step
# (slot) 1..25 - itest_type bascule sur '1' si un cod2 vaut litteralement
# '1' quelque part dans les mesures de l'echantillon (sinon reste a 0/defaut).
_CARTEST_DEFAULT = "123456789ABCDEFGHIJKLMNOP"  # itest_type==0 (defaut)
_CARTEST_ALT = "ABCDEFGHIJKLMNOPQRSTUVWXY"       # itest_type==1


def _cartest_table(mesures: List[Measurement]) -> str:
    return _CARTEST_ALT if any(m.cod2 == "1" for m in mesures) else _CARTEST_DEFAULT


def _cartest(ctype: str, table: str) -> int:
    """1-based, -1 si le caractere n'est pas dans la table (equivalent k=-1)."""
    idx = table.find(ctype)
    return idx + 1 if idx >= 0 else -1


def _step_of(m: Measurement) -> float:
    """Valeur de pas a AFFICHER/UTILISER en calcul (AraiPoint.temp,
    PtrmCheck.temp) : `step_value` (precis, voir testlect.Measurement) si
    disponible (.prmag), sinon `float(etape)` (comportement historique,
    .ren) - demande explicite utilisateur : garder `etape` (entier)
    inchange partout ailleurs dans l'application, n'utiliser la valeur
    precise que la ou elle sert reellement (paleointensite)."""
    return m.step_value if m.step_value is not None else float(m.etape)


@dataclass
class AraiPoint:
    """Un pas de temperature du diagramme d'Arai (equivalent d'un slot k)."""
    k: int
    temp: float = 0.0  # step_value (precis) si dispo, sinon float(etape) - voir _step_of
    xp: float = 0.0     # TRM acquis (axe horizontal), normalise par arno
    yp: float = 0.0     # NRM restant (axe vertical), normalise par arno
    dec: float = 0.0    # direction du pTRM acquis (V-R)
    winc: float = 0.0
    decl: float = 0.0   # direction du NRM restant (R, moyenne R+V ligne 505)
    aincl: float = 0.0
    nrm_vec: Tuple[float, float, float] = (0.0, 0.0, 0.0)  # vecteur brut (R+V)/2, pour l'ACP anchoree/libre


@dataclass
class PtrmCheck:
    """Une verification pTRM (jusqu'a 2 par slot k, ntest=1 ou 2). `xt` =
    pTRM re-mesuree (verification vs la V d'origine du slot cible) ; `yt`
    = position NRM pour le trace (moyenne R+V a l'etape du check, PAS
    xtptrm - ne pas confondre) ; `xtptrm` = pTRM "normale" mesuree a
    l'etape du check lui-meme (difference V-R a CETTE etape, pas une
    moyenne) - reference pour `ecart = xt-xtptrm` (table texte, PAS
    utilisee pour le trace graphique)."""
    k: int
    ntest: int
    xt: float = 0.0
    yt: float = 0.0
    xtptrm: float = 0.0
    temp: float = 0.0  # step_value (precis) si dispo, sinon float(etape) - voir _step_of


@dataclass
class AraiFit:
    """Ajustement lineaire + statistiques de qualite (lignes 1153-1275)."""
    a: float
    b: float
    n1: int
    n2: int
    sigma: float = 0.0
    ccr: float = 0.0
    f: float = 0.0
    g: float = 0.0
    qq: float = 0.0
    hlab: float = 0.0
    h: float = 0.0


def apply_thelli(mesures: List[Measurement]) -> List[Measurement]:
    """Equivalent de `thelli` (plotorthog.f:697-778), appelee par
    `zijderplot` UNIQUEMENT quand des codes R/V sont presents (ligne
    1054-1056 : `itrr` puis `call thelli`) - PAS `zijder2` (le Zijderveld
    plein-page standalone, qui affiche chaque mesure brute telle quelle).
    Demande explicite utilisateur, scopee au panneau Zijderveld de "View
    Paleoint Results" ("in the original plot, only the NRM is shown by
    combining R and V steps... on the zijderveld plot... within the panel
    including an arai, zijderveld and stereo").

    Regle (transcrite du Fortran) : R et P sont entierement ECARTES (jamais
    affiches) ; chaque V est remplace par UN SEUL point = moyenne de V et
    de la mesure BRUTE precedente (le R du meme palier, dans l'ordre
    normal Thellier R-puis-V) ; meme mecanisme pour des paires X/X, Y/Y,
    Z/Z consecutives (ATRM) - la premiere de la paire est ecartee, la
    seconde remplacee par la moyenne des deux ; D/N/T/K/S passent
    inchanges ; tout code non reconnu est ecarte. Si mesures ne contient
    ni 'R' ni 'V', retourne la liste TELLE QUELLE (comportement Fortran :
    thelli n'est appelee que si itrr==1)."""
    if not any(m.cod1 in ("R", "V") for m in mesures):
        return mesures

    out: List[Measurement] = []
    for i, m in enumerate(mesures):
        cod1 = m.cod1
        if cod1 in ("D", "N", "T", "K", "S"):
            out.append(m)
        elif cod1 in ("R", "P"):
            continue
        elif cod1 == "V":
            prev = mesures[i - 1] if i > 0 else m
            out.append(_dc_replace(
                m, cod1="D",
                x=(prev.x + m.x) / 2.0, y=(prev.y + m.y) / 2.0, z=(prev.z + m.z) / 2.0,
            ))
        elif cod1 in ("X", "Y", "Z") and i > 0 and mesures[i - 1].cod1 == cod1:
            prev = mesures[i - 1]
            out.append(_dc_replace(
                m, cod1="D",
                x=(prev.x + m.x) / 2.0, y=(prev.y + m.y) / 2.0, z=(prev.z + m.z) / 2.0,
            ))
        # sinon (X/Y/Z isole ou 1er d'une paire, ou code non reconnu) : ecarte
    return out


def detect_method_and_hlab(mesures: List[Measurement], com: str = "") -> Tuple[str, float]:
    """Detecte automatiquement le protocole ('IZZI' ou 'THELLIER') et le
    champ labo (Hlab, uT) a partir des mesures elles-memes plutot que du
    champ `com:` (fragile, position fixe, doit etre tape a la main) -
    demande explicite utilisateur ("the method is IZZI if there are S
    codes and Thellier if there are R and V codes... extract the field
    value from treat_dc_field").

    Methode : 'IZZI' des qu'un cod1=='S' est present (convention "coe
    version" du Fortran, S=zero-field/R=in-field) ; 'THELLIER' si R ET V
    sont presents sans S (convention "Thellier ordering", R=zero-field/
    V=in-field, deja portee dans compute_arai). Par defaut 'THELLIER' si
    aucun des deux motifs n'est detecte (comportement historique).

    Hlab : premiere valeur `treat_dc_field` non-None trouvee parmi les
    mesures (uT, deja dans cette unite - voir testlect.read_prmag_file) ;
    a defaut (fichiers .ren, ce champ n'existe pas), repli sur l'ancien
    mecanisme `parse_com_field(com)["ichamp"]`."""
    cod1_set = {m.cod1 for m in mesures}
    if "S" in cod1_set:
        method = "IZZI"
    elif "R" in cod1_set and "V" in cod1_set:
        method = "THELLIER"
    else:
        method = "THELLIER"

    hlab = 0.0
    for m in mesures:
        if m.treat_dc_field is not None and m.treat_dc_field != 0.0:
            hlab = m.treat_dc_field
            break
    if hlab == 0.0:
        hlab = float(parse_com_field(com)["ichamp"])
    return method, hlab


def _compute_arai_coe(
    mesures: List[Measurement], itpr: bool, arno: float, table: str,
    xs: List[float], ys: List[float], zs: List[float],
) -> Tuple[dict, List[PtrmCheck]]:
    """Equivalent du bloc "coe version" (reponse==2, plotpaleoint2.f
    lignes 292-461) : S=zero-field (point isole, PAS moyenne avec le pas
    suivant - contrairement au 'R' de la convention Thellier), R=in-field
    (difference vectorielle avec le 'S' de MEME etape). AUCUNE division
    par 2 sur xp/yp ici (contrairement a compute_arai) : le Fortran source
    utilise `/arno` seul dans ce bloc, pas `/(arno*2)` - verifie a la
    lecture (lignes 352,392-393 vs 174,252 du bloc reponse==1).

    winc EST DE SIGNE OPPOSE a la convention Thellier ici (pas de '-' -
    ligne 399-400 du Fortran, a comparer a la ligne 558-559 du bloc
    reponse==1 qui neglige avec '-') - difference reelle du source,
    preservee telle quelle, pas une incoherence de ce portage.

    NON PORTE (contrairement au reste) : le filtre cleanpaleo Q/L/F
    (seuls X/Y/Z sont ecartes dans ce bloc du Fortran, ligne 339 - Q/L/F
    ne le sont PAS, a la difference du bloc reponse==1 ou un bug reel de
    ce type avait ete corrige sur donnees reelles ; ici, faute de donnees
    reelles avec Q/L/F en IZZI pour verifier, le comportement du Fortran
    est reproduit tel quel).

    Verification des checks pTRM (P) : convention reponse==2 uniquement
    (xt = difference entre le pas P et le premier 'R' de MEME slot k
    trouve dans la sequence) - reponse==3 (alternative basee sur le pas
    juste avant le check) existe dans le Fortran mais n'est PAS portee ici
    (aucune donnee reelle pour distinguer laquelle s'applique - a
    verifier si un desaccord est observe sur un vrai jeu de donnees
    IZZI)."""
    points: dict = {}
    checks: List[PtrmCheck] = []

    def get_point(k: int) -> AraiPoint:
        if k not in points:
            points[k] = AraiPoint(k=k)
        return points[k]

    def vecdiff(i: int, j: int) -> float:
        if itpr:
            return abs(zs[i] - zs[j]) / arno
        return math.sqrt(
            (xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2 + (zs[i] - zs[j]) ** 2
        ) / arno

    n = len(mesures)
    for i in range(1, n):
        cod1 = mesures[i].cod1
        if cod1 in ("X", "Y", "Z"):  # voir docstring : Q/L/F non filtres ici, fidele au Fortran
            continue

        if cod1 == "S":
            k = _cartest(mesures[i].cod2, table)
            if k < 0:
                continue
            p = get_point(k)
            p.yp = math.sqrt(xs[i] ** 2 + ys[i] ** 2 + zs[i] ** 2) / arno
            p.nrm_vec = (xs[i], ys[i], zs[i])
            p.decl = math.degrees(math.atan2(ys[i], xs[i])) % 360.0
            p.aincl = math.degrees(math.atan2(zs[i], math.hypot(xs[i], ys[i])))
            continue

        if cod1 == "P":
            k = _cartest(mesures[i].cod2, table)
            if k < 0:
                continue
            ntest = 1
            if i > 0 and mesures[i - 1].cod1 == "P":
                if _cartest(mesures[i - 1].cod2, table) == k:
                    ntest = 2
            check = PtrmCheck(k=k, ntest=ntest, temp=_step_of(mesures[i]))
            found_r = False
            found_s = False
            for ji in range(1, n):
                if mesures[ji].cod1 == "R":
                    ki = _cartest(mesures[ji].cod2, table)
                    if ki != k:
                        continue
                    check.xt = vecdiff(ji, i)
                    found_r = True
                    break
                if found_s:
                    continue
                if mesures[ji].cod1 != "S":
                    continue
                if mesures[ji].etape == mesures[i].etape:
                    check.yt = math.sqrt(xs[ji] ** 2 + ys[ji] ** 2 + zs[ji] ** 2) / arno
                    check.xtptrm = vecdiff(i - 1, i)
                    found_s = True
            if found_r or found_s:
                checks.append(check)
            continue

        # sinon ('R') : etape en champ - xp/dec/winc par difference avec le
        # 'S' de MEME temperature (recherche en arriere puis en avant,
        # meme mecanisme que le bloc reponse==1 mais roles R/S inverses).
        k = _cartest(mesures[i].cod2, table)
        if k < 0:
            continue
        j = i - 1
        etap = mesures[i].etape
        found = False
        while j >= 0:
            if mesures[j].cod1 == "S":
                found = True
                break
            j -= 1
            if j >= 0 and mesures[j].etape != etap:
                j = -1
                break
        if not found:
            j = i + 1
            while j < n:
                if mesures[j].cod1 == "S":
                    found = True
                    break
                j += 1
                if j < n and mesures[j].etape != etap:
                    j = n
                    break
        if not found:
            continue
        p = get_point(k)
        p.xp = vecdiff(i, j)
        dx, dy, dz = xs[i] - xs[j], ys[i] - ys[j], zs[i] - zs[j]
        p.dec = math.degrees(math.atan2(dy, dx)) % 360.0
        p.winc = math.degrees(math.atan2(dz, math.hypot(dx, dy)))  # PAS de '-' ici, voir docstring
        p.temp = _step_of(mesures[i])

    return points, checks


def compute_arai(
    ech: SelectedSample, itpr: bool = False, method: Optional[str] = None,
    orientation: int = 1,
) -> Tuple[List[AraiPoint], List[PtrmCheck], float]:
    """Equivalent du tri/normalisation de `paleoin` : construit les points
    du diagramme d'Arai et les checks pTRM a partir des mesures brutes de
    `ech`. Deux conventions, auto-detectees (voir detect_method_and_hlab)
    si `method` n'est pas force explicitement :
    - 'THELLIER' (bloc reponse==1, lignes 464-629, convention N/R/V/P) -
      la SEULE historiquement portee/verifiee (echantillon 10CL1406A).
    - 'IZZI' (bloc "coe version", reponse==2, lignes 292-461, convention
      N/S/R/P) - portee via _compute_arai_coe, demande explicite
      utilisateur ("the method is IZZI if there are S codes").

    `itpr` : equivalent de `itpr==1` - si True, les differences vectorielles
    (xp, xt) n'utilisent que la composante Z (au lieu de la norme complete).

    `orientation` (1=echantillon/2=in-situ/3=pendage, meme convention que
    partout ailleurs) : applique a x/y/z AVANT tout calcul - une rotation
    orthonormale (corfor/corpen) preserve les normes, donc xp/yp/arno/H/Q/
    l'ajustement restent numeriquement identiques quelle que soit
    l'orientation ; seuls decl/aincl/dec/winc (directions) en dependent.
    AUPARAVANT toujours calcule en repere ECHANTILLON quel que soit
    `orientation` (bug reel - `fit_arai_direction`, elle, appliquait deja
    correctement l'orientation demandee, mais pas ce chemin-ci) : affectait
    `gamma` (90-winc, affiche dans les resultats) et le mini-stereo NRM/TRM
    du panneau de revisite - demande explicite utilisateur ("the plots are
    forced in Sample coordinates in this figures... within the panel
    including an arai, zijderveld and stereo").

    Retourne (points, checks, arno). `arno` = magnitude de la toute premiere
    mesure (NRM initiale, equivalent x(ideb)/y(ideb)/z(ideb))."""
    mesures = ech.mesures
    if len(mesures) < 2:
        return [], [], 0.0

    factor = 1.0e6 / (ech.vol or 1.0)  # equivalent x(jj)=mes.x*1e6/vol (toujours 1e6, pas de test norme ici)
    xs, ys, zs = [], [], []
    for m in mesures:
        x, y, z = apply_orientation(m.x * factor, m.y * factor, m.z * factor, ech, orientation)
        xs.append(x); ys.append(y); zs.append(z)

    # `arno` doit venir de la VRAIE NRM initiale (cod1='N'), pas
    # forcement mesures[0] : un specimen peut porter une experience
    # d'anisotropie de TRM (ATRM, D/X/Y/Z) AVANT le protocole de
    # paleointensite lui-meme - ce n'est pas une donnee desordonnee, c'est
    # une sequence reelle (l'ATRM travaille sur une TRM artificielle de
    # labo, pas la NRM naturelle - precision explicite de l'utilisateur).
    # Repli sur mesures[0] si aucun 'N' n'est trouve (comportement
    # d'origine, specimens sans ATRM en tete).
    ideb = next((i for i, m in enumerate(mesures) if m.cod1 == "N"), 0)
    arno = math.sqrt(xs[ideb] ** 2 + ys[ideb] ** 2 + zs[ideb] ** 2)
    if arno == 0.0:
        return [], [], 0.0

    table = _cartest_table(mesures)

    if method is None:
        method, _hlab = detect_method_and_hlab(mesures, ech.com)
    if method == "IZZI":
        points_dict, checks = _compute_arai_coe(mesures, itpr, arno, table, xs, ys, zs)
        result_points = [points_dict[k] for k in sorted(points_dict)]
        return result_points, checks, arno

    points: dict = {}  # k -> AraiPoint
    checks: List[PtrmCheck] = []

    def get_point(k: int) -> AraiPoint:
        if k not in points:
            points[k] = AraiPoint(k=k)
        return points[k]

    def vecdiff(i: int, j: int) -> float:
        if itpr:
            return abs(zs[i] - zs[j]) / (arno * 2.0)
        return math.sqrt(
            (xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2 + (zs[i] - zs[j]) ** 2
        ) / (arno * 2.0)

    n = len(mesures)
    for i in range(1, n):  # equivalent do i=2,il (1-based) -> range(1,n) (0-based, skip mesures[0]=NRM initiale)
        cod1 = mesures[i].cod1
        # equivalent de `cleanpaleo` (plotpaleoint2.f, appele en tete de
        # `paleoin` - "step removed from list: X/Y/Z/Q/L/F") : X/Y/Z sont
        # les pas d'anisotropie de TRM, Q/L le refroidissement rapide/lent,
        # F un pas auxiliaire - aucun n'est un pas Thellier valide. Bug reel
        # corrige : Q/L n'etaient pas filtres ici, donc tombaient dans la
        # branche 'V' par defaut ci-dessous et pouvaient ECRASER le bon
        # point si leur cod2 coincidait avec celui d'un vrai V (verifie sur
        # l'echantillon reel 06A, etape 530 : cod2='J' partage par V/L/Q,
        # 'Q' etant traite en dernier ecrasait xp/dec/winc du vrai 'V').
        if cod1 in ("X", "Y", "Z", "Q", "L", "F"):
            continue

        if cod1 == "R":
            # yp/decl/aincl = moyenne du point R et du point SUIVANT (le V
            # associe a la meme temperature) - transcrit tel quel (ligne
            # ~505-514), meme si la justification physique exacte de cette
            # moyenne n'est pas documentee dans le Fortran.
            k = _cartest(mesures[i].cod2, table)
            if k < 0 or i + 1 >= n:
                continue
            p = get_point(k)
            sx, sy, sz = xs[i] + xs[i + 1], ys[i] + ys[i + 1], zs[i] + zs[i + 1]
            p.yp = math.sqrt(sx * sx + sy * sy + sz * sz) / (arno * 2.0)
            p.nrm_vec = (sx / 2.0, sy / 2.0, sz / 2.0)
            p.decl = math.degrees(math.atan2(sy, sx)) % 360.0
            p.aincl = math.degrees(math.atan2(sz, math.hypot(sx, sy)))
            continue

        if cod1 == "P":
            # verification pTRM : k = slot du cod2 du POINT DE CONTROLE
            # lui-meme (peut correspondre a un slot LATER/different de la
            # temperature reellement re-testee - transcrit tel quel, lignes
            # 33130-33190). ntest=2 si le point precedent est aussi 'P' avec
            # le meme k (deux checks au meme "endroit" du graphique).
            k = _cartest(mesures[i].cod2, table)
            if k < 0:
                continue
            ntest = 1
            if i > 0 and mesures[i - 1].cod1 == "P":
                k_prev = _cartest(mesures[i - 1].cod2, table)
                if k_prev == k:
                    ntest = 2
            check = PtrmCheck(k=k, ntest=ntest, temp=_step_of(mesures[i]))
            found_v = False
            found_r = False
            for ji in range(1, n):
                if not found_v and mesures[ji].cod1 == "V":
                    ki = _cartest(mesures[ji].cod2, table)
                    if ki == k:
                        check.xt = vecdiff(ji, i)
                        found_v = True
                        continue
                if not found_v and not found_r and mesures[ji].cod1 == "R":
                    if mesures[ji].etape == mesures[i].etape and ji + 1 < n:
                        sx = xs[ji] + xs[ji + 1]
                        sy = ys[ji] + ys[ji + 1]
                        sz = zs[ji] + zs[ji + 1]
                        check.yt = math.sqrt(sx * sx + sy * sy + sz * sz) / (arno * 2.0)
                        # xtptrm : pTRM normale (difference V-R, PAS la
                        # somme comme yt) mesuree a l'etape du check -
                        # reference pour ecart=xt-xtptrm (plotpaleoint2.f:
                        # 619-624) - verifie octet-pres sur 06A.
                        check.xtptrm = vecdiff(ji + 1, ji)
                        found_r = True
                if found_v:
                    break
            if found_v:
                checks.append(check)
            continue

        # sinon ('V') : etape en champ - xp/dec/winc par difference avec le
        # 'R' de MEME temperature (recherche en arriere puis en avant,
        # lignes 3370 : ib==1 recule, ib==2 avance).
        k = _cartest(mesures[i].cod2, table)
        if k < 0:
            continue
        j = i - 1
        etap = mesures[i].etape
        found = False
        while j >= 0:
            if mesures[j].cod1 == "R":
                found = True
                break
            j -= 1
            if j >= 0 and mesures[j].etape != etap:
                j = -1
                break
        if not found:
            j = i + 1
            while j < n:
                if mesures[j].cod1 == "R":
                    found = True
                    break
                j += 1
                if j < n and mesures[j].etape != etap:
                    j = n
                    break
        if not found:
            continue
        p = get_point(k)
        p.xp = vecdiff(i, j)
        dx, dy, dz = xs[i] - xs[j], ys[i] - ys[j], zs[i] - zs[j]
        p.dec = math.degrees(math.atan2(dy, dx)) % 360.0
        p.winc = -math.degrees(math.atan2(dz, math.hypot(dx, dy)))
        p.temp = _step_of(mesures[i])

    result_points = [points[k] for k in sorted(points)]
    return result_points, checks, arno


def fit_arai_line(points: List[AraiPoint], n1: int, n2: int, hlab: float = 0.0) -> AraiFit:
    """Equivalent du calcul de pente/intercept + statistiques de qualite
    (lignes 1153-1248). `n1`,`n2` : position 1-based dans `points` (PAS le
    slot k brut - meme convention que `jdeb`/`jfin` dans ouvrir_ajuslig_dialog)."""
    sub = points[n1 - 1:n2]
    m = len(sub)
    if m < 2:
        raise ValueError("At least 2 points are required to fit a line.")
    xm = sum(p.xp for p in sub) / m
    ym = sum(p.yp for p in sub) / m
    x2 = sum((p.xp - xm) ** 2 for p in sub)
    y2 = sum((p.yp - ym) ** 2 for p in sub)
    xy = sum((p.xp - xm) * (p.yp - ym) for p in sub)
    b = -math.sqrt(y2 / x2) if x2 > 0 else 0.0
    a = ym - b * xm
    ccr = xy / math.sqrt(x2 * y2) if x2 > 0 and y2 > 0 else 0.0

    fit = AraiFit(a=a, b=b, n1=n1, n2=n2, hlab=hlab, h=abs(hlab * b), ccr=ccr)
    if hlab == 0.0 or m < 3 or x2 <= 0:
        return fit

    sigma = math.sqrt(max(0.0, (2 * y2 - 2 * b * xy) / ((m - 2) * x2)))
    zp = []
    for p in sub:
        z = b * p.xp + a
        z = p.yp + (z - p.yp) / 2.0 if z >= p.yp else p.yp - (p.yp - z) / 2.0
        zp.append(z)
    denom = zp[0] - zp[-1]
    f = denom / a if a != 0 else 0.0
    # bug reel corrige : le Fortran (plotpaleoint2.f:1220-1223) divise par
    # `denom` DEUX fois (une fois par terme dans la boucle, une fois de
    # plus a la fin) - equivalent a diviser par denom**2 une seule fois,
    # pas par denom une seule fois comme le faisait cette ligne avant -
    # verifie octet-pres contre un vrai transcript Fortran (echantillon
    # 06A, g=0.7752 attendu, 0.8615 obtenu avec l'ancienne formule).
    g_sum = sum((zp[i + 1] - zp[i]) ** 2 for i in range(len(zp) - 1))
    g = 1.0 - (g_sum / (denom * denom)) if denom != 0 else 0.0
    qq = abs(b) * f * g / sigma if sigma != 0 else 0.0

    fit.sigma = sigma
    fit.f = f
    fit.g = g
    fit.qq = qq
    return fit


@dataclass
class AraiDirection:
    """Direction ancree ET libre (ACP, `linearpal`, lignes 990-1151) +
    DANG (angle entre les deux) - PAS le meme ajustement que fit_arai_line
    (qui donne l'intensite, pas une direction)."""
    anchored_dec: Optional[float] = None
    anchored_inc: Optional[float] = None
    anchored_mad: Optional[float] = None
    anchored_specimen_frame: Optional[Tuple[float, float, float]] = None  # avant rotation orientation - pour compute_anicor_factor
    free_dec: Optional[float] = None
    free_inc: Optional[float] = None
    free_mad: Optional[float] = None
    nb: int = 0
    dang: Optional[float] = None


def fit_arai_direction(
    points: List[AraiPoint], n1: int, n2: int,
    ech: SelectedSample, orientation: int = 1,
) -> AraiDirection:
    """Equivalent de `linearpal` applique deux fois (ancree 'o' puis libre
    'n', lignes 990-1151) sur les vecteurs NRM bruts (nrm_vec) des points
    n1..n2 - meme algorithme ACP que `calcul.linear_fit` (deja utilise
    pour ajuslig/Zijderveld), avec le seuil MAD propre a la
    paleointensite (35 deg, pas 15 - `pmad>35` cote Fortran contre
    `itestlin<0`/MAD>15 pour ajuslig). La direction est tournee selon
    `orientation` (1/2/3, meme convention que le reste de l'appli) avant
    de calculer dec/inc. `dang` = angle entre les deux directions (Lisa
    Tauxe DANG)."""
    sub = points[n1 - 1:n2]
    result = AraiDirection(nb=len(sub))
    if len(sub) < 2:
        return result

    raw_pts = [p.nrm_vec for p in sub]
    anchored = linear_fit(raw_pts, anchored=True, mad_threshold=35.0)
    free = linear_fit(raw_pts, anchored=False, mad_threshold=35.0)

    dir_anchored = dir_free = None
    if anchored is not None:
        result.anchored_specimen_frame = tuple(anchored["direction"])
        dx, dy, dz = apply_orientation(*anchored["direction"], ech, orientation)
        _, dec, inc = polere(dx, dy, dz)
        result.anchored_dec, result.anchored_inc, result.anchored_mad = dec, inc, anchored["mad"]
        dir_anchored = (dx, dy, dz)
    if free is not None:
        dx, dy, dz = apply_orientation(*free["direction"], ech, orientation)
        _, dec, inc = polere(dx, dy, dz)
        result.free_dec, result.free_inc, result.free_mad = dec, inc, free["mad"]
        dir_free = (dx, dy, dz)

    if dir_anchored is not None and dir_free is not None:
        na = math.sqrt(sum(c * c for c in dir_anchored))
        nf = math.sqrt(sum(c * c for c in dir_free))
        if na > 0 and nf > 0:
            cosang = sum(dir_anchored[i] * dir_free[i] for i in range(3)) / (na * nf)
            result.dang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))

    return result


@dataclass
class CrmResult:
    """Tableau de contamination CRM (plotpaleoint2.f:1152-1247, bloc
    '544'). `values` : k (position 1-based dans `points`, entre n1 et n2)
    -> valeur crm(k). `crmmax`/`dtrm`/`rcrm` : diagnostics finaux (rcrm =
    '% rcrm' de la ligne de resultats)."""
    values: dict
    crmmax: float = 0.0
    dtrm: float = 0.0
    rcrm: float = 0.0


def compute_crm(points: List[AraiPoint], n1: int, n2: int, arno: float) -> CrmResult:
    """Equivalent du calcul `crm(k)` (plotpaleoint2.f:1152-1247) : test de
    contamination CRM par comparaison de l'inclinaison NRM restante au pas
    k contre celle du premier pas du fit (n1). Transcrit tel quel, y
    compris la constante '3.14152927' (PAS 3.1415927 = pi) utilisee au
    denominateur - coquille reelle du source Fortran, sans consequence
    notable sur la valeur numerique (ecart ~1e-8) mais preservee pour
    fidelite. `n1`,`n2` : position 1-based dans `points`, meme convention
    que `fit_arai_line`."""
    aincl_n1 = points[n1 - 1].aincl
    values: dict = {}
    crmmax = 0.0
    for k in range(n1, n2 + 1):
        p = points[k - 1]
        num = math.sin((aincl_n1 - p.aincl) * 3.1415927 / 180.0)
        den = math.sin((aincl_n1 - 90.0) * 3.14152927 / 180.0)
        val = (num / den) * arno * p.yp if den != 0 else 0.0
        values[k] = val
        crmmax = max(crmmax, abs(val))
    dtrm = points[n2 - 1].xp - points[n1 - 1].xp
    crmmax_norm = crmmax / arno if arno else 0.0
    rcrm = (crmmax_norm / dtrm) * 100.0 if dtrm else 0.0
    return CrmResult(values=values, crmmax=crmmax_norm, dtrm=dtrm, rcrm=rcrm)


@dataclass
class CurvatureResult:
    """Resultat de `AraiCurvature` (adjustcircle.f95) : ajustement d'un
    cercle (Taubin puis Levenberg-Marquardt geometrique, Chernov & Lesort
    2005) sur les points (xp,yp) normalises par leur maximum, comme test
    de courbure de Paterson (2011)."""
    k: float
    sse: float
    taubin_a: float
    taubin_b: float
    taubin_r: float
    lma_a: float
    lma_b: float
    lma_r: float


def _taubin_svd(xy: np.ndarray) -> Tuple[float, float, float]:
    """Equivalent de `TaubinSVD` (adjustcircle.f95:67-140) : ajustement
    algebrique de cercle par Taubin (1991), via SVD."""
    centroid = xy.mean(axis=0)
    x = xy[:, 0] - centroid[0]
    y = xy[:, 1] - centroid[1]
    z = x * x + y * y
    zmean = float(z.mean())
    z0 = (z - zmean) / (2.0 * math.sqrt(zmean))
    zxy = np.column_stack([z0, x, y])
    _, _, vt = np.linalg.svd(zxy, full_matrices=False)
    a0, a1, a2 = vt[2, 0], vt[2, 1], vt[2, 2]
    a0n = a0 / (2.0 * math.sqrt(zmean))
    a3 = -zmean * a0n
    ax = -a1 / (a0n * 2.0) + centroid[0]
    ay = -a2 / (a0n * 2.0) + centroid[1]
    br = math.sqrt(a1 * a1 + a2 * a2 - 4.0 * a0n * a3) / abs(a0n) / 2.0
    return float(ax), float(ay), float(br)


def _var_circle(xy: np.ndarray, a: float, b: float, r: float) -> float:
    """Equivalent de `VarCircle` (adjustcircle.f95:405-432)."""
    n = len(xy)
    d = np.sqrt((xy[:, 0] - a) ** 2 + (xy[:, 1] - b) ** 2) - r
    return float(np.sum(d * d) / (n - 3))


def _lma_circle(xy: np.ndarray, par_ini: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """Equivalent de `LMA` (adjustcircle.f95:163-404) : raffinement
    geometrique (Chernov & Lesort 2005) de l'estimation initiale de
    Taubin, transcrit ligne a ligne (memes noms de variables, meme
    structure GOTO traduite en boucles/`break`/`continue`). Ecarts
    numeriques de l'ordre de 1e-3 a 1e-4 attendus face au Fortran d'origine
    (real*4 simple precision contre float64 ici, sur ~50 iterations d'un
    ajustement non lineaire sensible) - verifie sur l'echantillon reel 06A
    (k=-0.7734 ici contre -0.7736 reel, SSE identique a la precision
    affichee)."""
    n = len(xy)
    factor_up = 10.0
    factor_down = 0.04
    lam = 0.01
    epsilon1 = 0.000001
    two_pi = 2 * math.pi
    iter_max = 50
    adjust_max = 20
    xshift = 0.0
    yshift = 0.0
    dx_shift = 0.1
    dy_shift = 0.0

    a0, b0, r0 = par_ini
    a_new = a0 + xshift
    b_new = b0 + yshift
    anew = 1.0 / (2.0 * r0)
    aabb = a_new * a_new + b_new * b_new
    fnew = (aabb - r0 * r0) * anew
    tnew = math.acos(max(-1.0, min(1.0, -a_new / math.sqrt(aabb))))
    if b_new > 0:
        tnew = two_pi - tnew

    var_new = _var_circle(xy, a0, b0, r0)
    finish = 0
    aold = fold = told = var_old = None

    for _outer in range(iter_max):
        aold, fold, told, var_old = anew, fnew, tnew, var_new
        h = math.sqrt(1 + 4 * aold * fold)
        a_old = -h * math.cos(told) / (2 * aold) - xshift
        b_old = -h * math.sin(told) / (2 * aold) - yshift
        r_old = 1.0 / abs(2 * aold)

        d = math.sqrt(1 + 4 * aold * fold)
        ct, st = math.cos(told), math.sin(told)
        h11 = h12 = h13 = h22 = h23 = h33 = 0.0
        f1 = f2 = f3 = 0.0
        for i in range(n):
            xi = xy[i, 0] + xshift
            yi = xy[i, 1] + yshift
            zi = xi * xi + yi * yi
            ui = xi * ct + yi * st
            vi = -xi * st + yi * ct
            adf = aold * zi + d * ui + fold
            sq = math.sqrt(4 * aold * adf + 1)
            den = sq + 1
            gi = 2 * adf / den
            fact = 2 / den * (1 - aold * gi / sq)
            dgdai = fact * (zi + 2 * fold * ui / d) - gi * gi / sq
            dgdfi = fact * (2 * aold * ui / d + 1)
            dgdti = fact * d * vi
            h11 += dgdai * dgdai
            h12 += dgdai * dgdfi
            h13 += dgdai * dgdti
            h22 += dgdfi * dgdfi
            h23 += dgdfi * dgdti
            h33 += dgdti * dgdti
            f1 += gi * dgdai
            f2 += gi * dgdfi
            f3 += gi * dgdti

        for _adjust in range(adjust_max):
            g11 = math.sqrt(h11 + lam)
            g12 = h12 / g11
            g13 = h13 / g11
            g22 = math.sqrt(h22 + lam - g12 * g12)
            g23 = (h23 - g12 * g13) / g22
            g33 = math.sqrt(h33 + lam - g13 * g13 - g23 * g23)
            d1 = f1 / g11
            d2 = (f2 - g12 * d1) / g22
            d3 = (f3 - g13 * d1 - g23 * d2) / g33
            dt = d3 / g33
            df = (d2 - g23 * dt) / g22
            da = (d1 - g12 * df - g13 * dt) / g11

            anew = aold - da
            fnew = fold - df
            tnew = told - dt
            xxt = 1 + 4 * anew * fnew
            if xxt < epsilon1 and lam > 1.0:
                xshift += dx_shift
                yshift += dy_shift
                h = math.sqrt(1 + 4 * aold * fold)
                a_temp = -h * math.cos(told) / (2 * aold) + dx_shift
                b_temp = -h * math.sin(told) / (2 * aold) + dy_shift
                r_temp = 1.0 / abs(2 * aold)
                anew = 1.0 / (2 * r_temp)
                aabb = a_temp * a_temp + b_temp * b_temp
                fnew = (aabb - r_temp * r_temp) * anew
                tnew = math.acos(max(-1.0, min(1.0, -a_temp / math.sqrt(aabb))))
                if b_temp > 0:
                    tnew = two_pi - tnew
                var_new = var_old
                break

            if 1 + 4 * anew * fnew < epsilon1:
                lam *= factor_up
                continue

            d = math.sqrt(1 + 4 * anew * fnew)
            ct, st = math.cos(tnew), math.sin(tnew)
            gg = 0.0
            for i in range(n):
                xi = xy[i, 0] + xshift
                yi = xy[i, 1] + yshift
                zi = xi * xi + yi * yi
                ui = xi * ct + yi * st
                adf = anew * zi + d * ui + fnew
                sq = math.sqrt(4 * anew * adf + 1)
                den = sq + 1
                gi = 2 * adf / den
                gg += gi * gi
            var_new = gg / (n - 3)
            h = math.sqrt(1 + 4 * anew * fnew)
            a_new = -h * math.cos(tnew) / (2 * anew) - xshift
            b_new = -h * math.sin(tnew) / (2 * anew) - yshift
            r_new = 1.0 / abs(2 * anew)

            if var_new <= var_old:
                progress = (abs(a_new - a_old) + abs(b_new - b_old) + abs(r_new - r_old)) / (r_new + r_old)
                if progress < epsilon1:
                    aold, fold, told, var_old = anew, fnew, tnew, var_new
                    finish = 1
                    break
                lam *= factor_down
                break
            else:
                lam *= factor_up
                continue

        if finish == 1:
            break

    h = math.sqrt(1 + 4 * aold * fold)
    par_a = -h * math.cos(told) / (2 * aold) - xshift
    par_b = -h * math.sin(told) / (2 * aold) - yshift
    par_r = 1.0 / abs(2 * aold)
    return par_a, par_b, par_r


def arai_curvature(xs: List[float], ys: List[float]) -> CurvatureResult:
    """Equivalent de `AraiCurvature` (adjustcircle.f95:1-65) : normalise
    (xs,ys) par leur maximum respectif, ajuste un cercle (Taubin puis
    LMA), et deduit le signe de la courbure `k=1/r` en comparant le centre
    du cercle (repere NORMALISE) a la moyenne des points BRUTS (pas
    normalisee) - transcrit tel quel, y compris cette incoherence
    d'echelle du Fortran d'origine (sans consequence pratique observee :
    xp/yp sont deja des fractions proches de 1)."""
    xs_arr = np.asarray(xs, dtype=float)
    ys_arr = np.asarray(ys, dtype=float)
    max_x = float(xs_arr.max())
    max_y = float(ys_arr.max())
    meanx = float(xs_arr.mean())
    meany = float(ys_arr.mean())
    xn = xs_arr / max_x
    yn = ys_arr / max_y
    xy = np.column_stack([xn, yn])

    taubin = _taubin_svd(xy)
    best_a, best_b, best_r = _lma_circle(xy, taubin)

    if best_a <= meanx and best_b <= meany:
        rk = -1.0 / best_r
    else:
        rk = 1.0 / best_r
    d = np.sqrt((xn - best_a) ** 2 + (yn - best_b) ** 2) - best_r
    sse = float(np.sum(d * d))

    return CurvatureResult(
        k=rk, sse=sse,
        taubin_a=taubin[0], taubin_b=taubin[1], taubin_r=taubin[2],
        lma_a=best_a, lma_b=best_b, lma_r=best_r,
    )


_METHOD_TO_REPONSE = {"T": 1, "C": 2, "I": 3}


def parse_com_field(com: str) -> dict:
    """Equivalent du parsing du champ `com:` (plotpaleoint2.f:97-122,
    format Fortran `(i2,a1,a1,1x,a2,1x,i2,i2)`) : les 12 caracteres du
    champ com encodent `<champ 2 chiffres><methode><air/vide> <four><mois>
    <annee>`, ex. 'com:50' -> champ=50µT (methode/air_vide/four/date
    absents). `methode` ('T'=Thellier,'C'=Coe,'I'=Tauxe/IZZI, sinon/vide
    -> Thellier par defaut - seul Thellier est porte, voir docstring
    module) selectionne le protocole ET remplace le choix de
    `Parameters Paleointensity`. `air_vide`/`four`/mois/annee sont
    PUREMENT INFORMATIFS cote Fortran (aucune correction n'en depend,
    confirme a la lecture du source). Tolerant : champ vide/court/mal
    forme -> valeurs par defaut, ne leve jamais d'exception."""
    text = (com or "").ljust(12)
    result = {
        "has_field": bool(text[:2].strip()),
        "ichamp": 0, "method": "", "reponse": 1,
        "air_vide": "", "four": "", "imois": 0, "iannee": 0,
    }
    try:
        result["ichamp"] = int(text[0:2].strip() or 0)
    except ValueError:
        pass
    method = text[2:3].strip().upper()
    result["method"] = method
    result["reponse"] = _METHOD_TO_REPONSE.get(method, 1)
    result["air_vide"] = text[3:4].strip().upper()
    result["four"] = text[5:7].strip().upper()
    try:
        result["imois"] = int(text[8:10].strip() or 0)
    except ValueError:
        pass
    try:
        result["iannee"] = int(text[10:12].strip() or 0)
    except ValueError:
        pass
    return result


_SYMBO1_MAP = {3: 5, 5: 8, 2: 10, 11: 14}


def draw_arai(
    ctx,
    ech: SelectedSample,
    points: List[AraiPoint],
    checks: List[PtrmCheck],
    arno: float,
    fit: Optional[AraiFit] = None,
    dimens: float = 10.0,
    step_labels: Optional[List[int]] = None,
) -> None:
    """Port du trace du diagramme d'Arai (paleoin, lignes 819-1331). `ctx`
    deja positionne par l'appelant sur l'origine du panneau (equivalent d'un
    `ctx.plot(0,0,-3)` prealable, comme build_zijderveld_figure/
    build_stereo_figure). `dimens` : pas encore calibre contre un export reel
    (contrairement a dimzij=15/dimstereo=21) - valeur de depart plausible.

    Les 4 helpers `plott`/`numbe`/`symbo2`/`symbo1` reproduisent EXACTEMENT
    les wrappers Fortran du meme nom (lignes 1360-1395 de plotpaleoint2.f) :
    ils negent/permutent (x,y) avant d'appeler les primitives reelles - c'est
    ce mecanisme qui fait que le diagramme sort avec NRM vertical et TRM
    horizontal sur la page."""

    def plott(x, y, ip):
        ctx.plot(-y, x, ip)

    def numbe(x, y, dimm, value, angle, ndec):
        ctx.number(-y, x, dimm, value, 0.0, ndec)

    def symbo2(x, y, dimm, text, angle, nchar):
        ctx.plottxt(-y, x, dimm, text, 0.0, nchar=nchar)

    def symbo1(x, y, dimm, ip, angle, ipen):
        ctx.symbol(-y, x, dimm, _SYMBO1_MAP.get(ip, ip), ipen)

    if not points:
        return

    dimplot = dimens
    dimtext = dimens

    xmax = max(p.xp for p in points)
    ymax = max(p.yp for p in points)
    ymin = min(p.yp for p in points)

    ctx.newpen(1)
    ctx.thickn(1.0)
    plott(5.5, -1.2, -3)
    if ymax > 1.1:
        dimplot = dimplot / ymax

    # axe NRM (vertical) - lignes 832-848
    plott(-dimplot / 20.0, 0.0, 3)
    symbo1(dimplot, 0.0, dimtext / 20.0, 3, 0.0, -2)
    plott(dimplot, 0.0, 3)
    if ymax <= 1.0:
        plott(dimplot * 1.1, 0.0, 2)
        symbo2(dimplot * 1.2, 0.15 * dimplot, dimtext / 20.0, "NRM", 270.0, 3)
    else:
        plott((ymax + 0.1) * dimplot, 0.0, 2)
        symbo2((ymax + 0.2) * dimplot, 0.15 * dimplot, dimtext / 20.0, "NRM", 270.0, 3)
    symbo2(dimplot * 0.975, dimplot / 10.0, dimtext / 20.0, "1", 270.0, 1)
    plott(0.0, dimplot / 20.0, 3)

    # axe TRM (horizontal) - 3 branches selon xmax/ymin, lignes 849-928
    if xmax > 1.6:
        phi = xmax / 1.3
        plott(0.0, -(xmax + 0.1) * dimplot / phi, 2)
        symbo2(-dimplot / 10.0, -(xmax - 0.1) * dimplot / phi, dimtext / 20.0, "TRM", 270.0, 3)
        iphi = int(xmax)
        for i in range(1, iphi + 1):
            symbo1(0.0, -i * dimplot / phi, dimtext / 20.0, 3, 0.0, -1)
            numbe(-dimplot / 10.0, (-i * dimplot / phi) + dimplot / 40.0, dimtext / 20.0, float(i), 270.0, -1)
    elif ymin < 0.2 and (xmax / (1.0 - ymin)) < 0.6:
        phi = max(xmax, 0.65)
        plott(0.0, -(phi + 0.6) * dimplot / phi, 2)
        symbo2(-dimplot / 10.0, -(phi - 0.1) * dimplot / phi, dimtext / 20.0, "TRM", 270.0, 3)
        symbo1(0.0, -dimplot / phi, dimtext / 20.0, 3, 0.0, -1)
        numbe(-dimplot / 10.0, (-dimplot / phi) + dimplot / 40.0, dimtext / 20.0, 1.0, 270.0, -1)
    else:
        phi = 1.0
        symbo1(0.0, -dimplot, dimtext / 20.0, 3, 270.0, -2)
        if xmax >= 1.0:
            plott(0.0, -(xmax + 0.1) * dimplot, 2)
            symbo2(-dimplot / 10.0, -(xmax + 0.1) * dimplot, dimtext / 20.0, "TRM", 270.0, 3)
        else:
            plott(0.0, -dimplot * 1.1, 2)
            symbo2(-dimplot / 10.0, -dimplot * 1.1, dimtext / 20.0, "TRM", 270.0, 3)
        symbo2(-dimplot / 10.0, -dimplot * 0.99, dimtext / 20.0, "1", 270.0, 1)

    # points + checks pTRM - lignes 899-915
    checks_by_k: dict = {}
    for c in checks:
        checks_by_k.setdefault(c.k, []).append(c)
    for p in points:
        u = p.yp * dimplot
        v = -p.xp * dimplot / phi
        ctx.newpen(3)
        symbo1(u, v, dimtext / 25.0, 5, 270.0, -1)
        ctx.newpen(5)
        for c in checks_by_k.get(p.k, [])[:2]:
            a_ = c.yt * dimplot
            b_ = -c.xt * dimplot / phi
            if c.ntest == 1:
                plott(u, v, 3)
            else:
                plott(u, b_, 3)
            plott(u, b_, 2)
            symbo1(a_, b_, dimtext / 25.0, 2, 270.0, -2)

    # id echantillon - ligne 916-919
    ctx.newpen(1)
    symbo2(-dimplot / 5.0, -dimplot / 4.0, dimtext / 20.0, ech.id, 270.0, 12)

    # texte NRM - lignes 921-929 (nchar = longueur REELLE du texte forme,
    # equivalent `nlen=len(trim(text))` - pas une valeur figee, sous peine
    # de tronquer le texte comme le faisait une premiere version buguee)
    if ech.norme == "v":
        nrm_text = f"NRM ={arno:6.2f} A/m"
    else:
        nrm_text = f"NRM ={arno * 0.001:6.4f} Am2/kg"
    symbo2(dimplot * 1.2, -dimplot / 2.0, dimtext / 25.0, nrm_text, 270.0, len(nrm_text))

    # droite ajustee + points re-marques - lignes 1185-1197
    if fit is not None:
        ctx.newpen(3)
        plott(fit.a * dimplot, 0.0, 3)
        if fit.b != 0.0 and (-fit.a / fit.b) >= xmax:
            plott(0.0, -(xmax * dimplot / phi), 3)
            plott(0.0, -(xmax + 0.1) * dimplot / phi, 2)
        plott(fit.a * dimplot, 0.0, 3)
        if fit.b != 0.0:
            plott(0.0, fit.a * dimplot / fit.b / phi, 2)
        for p in points[fit.n1 - 1: fit.n2]:
            u = p.yp * dimplot
            v = -p.xp * dimplot / phi
            symbo1(u, v, dimtext / 25.0, 11, 0.0, -1)

        # bloc Hlab/Q/H - lignes 1228-1275
        if fit.hlab != 0.0:
            ctx.newpen(1)
            mu = "µ"
            hlab_text = f"Hlab = {fit.hlab:5.1f} {mu}T"
            q_text = f"Q    = {fit.qq:5.1f}"
            h_text = f"H    = {fit.h:5.1f} {mu}T"
            symbo2(dimplot * 1.1, -dimplot, dimtext / 25.0, hlab_text, 270.0, len(hlab_text))
            symbo2(dimplot * 0.9, -dimplot, dimtext / 25.0, q_text, 270.0, len(q_text))
            symbo2(dimplot * 1.0, -dimplot, dimtext / 25.0, h_text, 270.0, len(h_text))

        # labels d'etapes (temperature) - lignes 1320-1331
        if step_labels:
            angu = 270.0
            if fit.b != 0.0:
                angu = angu + 90.0 - math.degrees(math.atan(abs(fit.b * phi)))
            for pos in step_labels:
                p = points[pos - 1]
                u = p.yp * dimplot
                v = -p.xp * dimplot / phi
                numbe(u, v - dimplot / 30.0, dimtext / 30.0, float(p.temp), angu, -1)

    ctx.newpen(1)


def build_arai_figure(
    ech: SelectedSample,
    points: List[AraiPoint],
    checks: List[PtrmCheck],
    arno: float,
    fit: Optional[AraiFit] = None,
    dimens: float = 10.0,
    step_labels: Optional[List[int]] = None,
    fig: Optional[Figure] = None,
) -> Figure:
    """Enveloppe matplotlib de `draw_arai`, meme pattern que
    `build_zijderveld_figure`/`build_stereo_figure` (zijderveld.py/stereo.py)."""
    if fig is None:
        fig = Figure(figsize=(6.0, 8.0), dpi=100)
    else:
        fig.clear()
    ax = fig.add_subplot(111)
    ctx = PlotContext(ax)
    ctx.clear()
    ctx.plot(0.0, 0.0, -3)

    draw_arai(ctx, ech, points, checks, arno, fit=fit, dimens=dimens, step_labels=step_labels)

    ax.relim()
    ax.autoscale_view()
    fig.tight_layout()
    return fig


def _draw_arai_stereo_directions(
    ctx: PlotContext, points: List[AraiPoint], r: float,
    la: float = -90.0, phi: float = 0.0, iproj: int = 1, point_size: float = 0.22,
) -> None:
    """Equivalent PARTIEL de `stereoNRMTRM` (plotster.f:746-880) : trace la
    direction NRM restante (rdec/rinc famille 1, symboles 8/14) et la
    direction TRM acquis (famille 2, symboles 9/15) de chaque point du
    diagramme d'Arai sur un reseau deja trace (`draw_stereo_net`, meme
    convention `superc`/symboles que `draw_stereo_measurements`). La 3e
    famille du Fortran (direction TRM incrementale point-a-point, calculee
    via `polere` sur la difference vectorielle entre points consecutifs)
    n'est PAS portee ici - diagnostic secondaire, hors perimetre pour
    l'instant."""
    ctx.thickn(0.5)
    ctx.newpen(3)
    for p in points:
        u, v, ifl = superc(la, phi, p.decl, p.aincl, iproj)
        ctx.symbol(v * r, u * r, point_size, 8 if ifl == 5 else 14, -1)
    ctx.newpen(4)
    for p in points:
        u, v, ifl = superc(la, phi, p.dec, p.winc, iproj)
        ctx.symbol(v * r, u * r, point_size, 9 if ifl == 5 else 15, -1)
    ctx.newpen(1)


def build_paleoint_review_figure(
    ech: SelectedSample,
    points: List[AraiPoint],
    checks: List[PtrmCheck],
    arno: float,
    fit: Optional[AraiFit] = None,
    orientation: int = 1,
    dimens: float = 7.5,
    fig: Optional[Figure] = None,
) -> Figure:
    """Equivalent de la mise en page combinee de `visi_paleoin`
    (visi_Paleoint.f, `boite(1)`/`boite(2)`/`boite(3)`) utilisee par "View
    Paleoint Results" pour revisiter rapidement un traitement deja
    effectue : diagramme d'Arai en haut (pleine largeur, boite 1),
    Zijderveld en bas a gauche (boite 2), stereo NRM restant/TRM acquis en
    bas a droite (boite 3, `stereoNRMTRM` partiel - voir
    `_draw_arai_stereo_directions`). PAS le panneau Arai seul de
    `afficher_arai`/`build_arai_figure` (mode d'edition interactive).

    Dimensions calees sur le Fortran (demande explicite utilisateur : "the
    scale of each graph is small... check the original Fortran") :
    - page `call plots(19.5,28.,fname)` (visi_Paleoint.f:436) -> figsize en
      pouces (19.5/2.54, 28/2.54) - le figsize precedent (7.5,9.5) etait
      nettement moins haut (24.1cm) que l'original (28cm), une cause
      plausible de "plots trop petits" une fois mis a l'echelle dans le
      panneau graphique de l'appli.
    - `dimens` par defaut 7.5 (pas 9.0) : valeur d'initialisation reelle du
      common /paleointensite/ (plotpaleoint2.f:1595/1604/1607 - "if
      (dimens==0.0) dimens=7.5"), pas une valeur inventee.
    - height_ratios [390,370] entre boite(1) (paleointensite, y:10-400,
      hauteur 390) et boite(2)/(3) (Zijderveld/stereo, y:420-790, hauteur
      370) - graphicsAWE.f95:1201-1247, remplace l'ancien [1.15,1.0]
      approximatif.
    - `dimzij=9.0` explicitement passe au panneau Zijderveld (visi_Paleoint.f
      lignes 253/277 : "dimzij=9.0", DIFFERENT du dimzij=15.0 de la vue
      Zijderveld pleine page - voir zijderveld.py:draw_zijderveld).

    NON RESOLU (signale, pas corrige silencieusement) : le Fortran de cette
    vue combinee appelle `zijderplot` (visi_Paleoint.f:439), PAS `zijder2`
    (la subroutine sur laquelle draw_zijderveld/zijder2 a ete calibre et
    verifie pixel-pres pour la vue pleine page - voir memoire projet
    project_starmac_fortran_gotchas.md) - deux subroutines Zijderveld
    distinctes existent dans plotorthog.f avec un rendu different
    (monochrome pour zijderplot vs rouge/vert pour zijder2). Ce panneau
    combine reutilise donc draw_zijderveld (donc la logique zijder2) par
    approximation, pas par fidelite verifiee a `zijderplot` - a
    l'utilisateur de dire si c'est acceptable ou si zijderplot merite son
    propre portage."""
    if fig is None:
        fig = Figure(figsize=(19.5 / 2.54, 28.0 / 2.54), dpi=100)
    else:
        fig.clear()
    gs = fig.add_gridspec(2, 2, height_ratios=[390, 370])

    ax_arai = fig.add_subplot(gs[0, :])
    ctx_arai = PlotContext(ax_arai)
    ctx_arai.clear()
    ctx_arai.plot(0.0, 0.0, -3)
    draw_arai(ctx_arai, ech, points, checks, arno, fit=fit, dimens=dimens)
    ax_arai.relim()
    ax_arai.autoscale_view()
    ax_arai.set_title(ech.id)

    ax_zij = fig.add_subplot(gs[1, 0])
    ctx_zij = PlotContext(ax_zij)
    ctx_zij.clear()
    ctx_zij.plot(0.0, 0.0, -3)
    # thelli (voir apply_thelli) : R/V combines en un seul point NRM par
    # palier, R/P ecartes - fidele a zijderplot (visi_Paleoint.f), PAS a
    # zijder2 (Zijderveld standalone, qui reste inchange - demande
    # explicite utilisateur, scopee a ce panneau).
    ech_zij = _dc_replace(ech, mesures=apply_thelli(ech.mesures))
    draw_zijderveld(ctx_zij, ech_zij, orientation=orientation, show_stereo=False, dimzij=9.0)
    ax_zij.relim()
    ax_zij.autoscale_view()

    ax_stereo = fig.add_subplot(gs[1, 1])
    ctx_stereo = PlotContext(ax_stereo)
    ctx_stereo.clear()
    ctx_stereo.plot(0.0, 0.0, -3)
    dimster = 13.0
    point_size = (0.18 * dimster) / 10.0
    r = draw_stereo_net(ctx_stereo, orientation, dimster=dimster, show_orient_label=False)
    _draw_arai_stereo_directions(ctx_stereo, points, r, point_size=point_size)
    ax_stereo.relim()
    ax_stereo.autoscale_view()

    # PAS tight_layout() (marges genereuses pensees pour des graduations/
    # labels d'axe - inutiles ici, les 3 panneaux sont dessines avec
    # ax.axis("off"), voir plotlib.PlotContext.clear) : des marges quasi
    # nulles explicites laissent chaque panneau (contraint en aspect
    # 'equal' par PlotContext, adjustable='datalim') occuper le maximum de
    # sa boite plutot que de rester petit avec de la marge inutilisee
    # autour - demande explicite utilisateur ("each of the three could be
    # larger within the page").
    fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.02, hspace=0.06, wspace=0.03)
    return fig
