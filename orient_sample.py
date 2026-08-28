"""
Declinaison magnetique locale a partir de l'azimuth solaire et de
l'azimuth magnetique deja presents dans les fichiers (orient_paleomag.f :
`sundec`/`angleorient`/`orient_sample`) - demande explicite utilisateur
("dans les fichiers il y a l'information sur l'azimuth solaire et
l'azimuth magnetique et cela permet de calculer la declinaison locale.
Reprendre du code Fortran les calculs").

Principe (technique classique de "sun compass" en paleomagnetisme) : sur
le terrain, on releve l'azimuth d'un MEME repere (encoche/ligne sur la
carotte) avec DEUX methodes independantes - une boussole magnetique
(azmag) et une visee solaire (azsun, a une heure connue). L'azimuth
GEOGRAPHIQUE VRAI de ce repere (azgeo) se calcule uniquement a partir de
la position du soleil (latitude/longitude/date/heure - AUCUN champ
magnetique impliqué) ; comparer azgeo a azmag donne alors la declinaison
magnetique LOCALE reellement mesuree sur site, independante de tout
modele geomagnetique global et potentiellement plus fiable localement
(anomalie magnetique, derive de l'instrument...) que la declinaison
IGRF (calculee separement ici a titre de comparaison/verification, pas
utilisee pour corriger les donnees).

IGRF : le Fortran d'origine (`igrfosx`/`igrf13syn`, IGRFstarmac.f) code
en dur les coefficients de Gauss IGRF-13 directement dans le source -
demande explicite utilisateur ("mettre les coefficients IGRF dans un
fichier externe pour les mises a jour tous les 5 ans") : ici, REUTILISE
`pmagpy.pmag.doigrf` plutot que de porter `igrf13syn` avec ses propres
DATA figees. pmagpy charge deja ses coefficients depuis un fichier TEXTE
EXTERNE (`pmagpy/field_models/igrf14coeffs.txt`, format standard NOAA/
IAGA, un jeu de coefficients par epoque tous les 5 ans + un modele de
variation seculaire pour l'epoque courante) - remplacer ce fichier par la
prochaine generation (IGRF-15, attendue ~2030) suffira a mettre a jour le
modele SANS toucher au code, exactement la demande de l'utilisateur. Deja
IGRF-14 (derniere generation, 2025) ici, plus recent que le IGRF-13 code
en dur dans le Fortran d'origine.
"""

import math
from typing import Tuple

from pmagpy import pmag

_F = math.pi / 180.0


def _sun_declination_and_eot(iy: int, m: int, d: float) -> Tuple[float, float]:
    """Port exact de `sundec` (orient_paleomag.f:442-488, "mise au point
    par Haraldur Audusson et Pierrick Roperch, Corvallis, aout 1988") :
    ephemeride solaire basse precision (declinaison + equation du temps)
    a partir du jour julien. `d` (jour du mois) peut etre FRACTIONNAIRE
    (heure/minute d'observation deja incorporees par l'appelant, voir
    compute_local_declination) - la formule du jour julien l'accepte
    directement (`tjd=...+d`).

    Retourne (declinaison_solaire_deg, equation_du_temps_deg)."""
    tjd = float(367 * iy - (7 * (iy + (m + 9) // 12)) // 4 + (275 * m) // 9) + 1721.014e3 + d
    t = tjd - 2451.545e3
    tt = (t / 36525.0) + 1.0

    eps = 23.452294 - (0.0130 * tt) + 0.0025 * math.cos((125.0435 - 0.0529538076 * t) * _F)
    e = 0.01675104 - (0.00004180 * tt)
    w = 281.220844 + (1.719175 * tt)
    rl = 280.466 + 0.9856473516 * t
    g = rl - w
    ec = (2.0 * e * math.sin(g * _F) + (5.0 / 4.0) * e * e * math.sin(2.0 * g * _F)) / _F
    dele = -0.0047 * math.sin((125.0435 - 0.0529538076 * t) * _F)
    rland = rl + ec + dele

    sun = math.sin(rland * _F) * math.sin(eps * _F)
    decli = math.degrees(math.asin(sun))

    a1 = math.tan(_F * eps / 2.0) ** 2
    a2 = 4.0 * e * a1
    a3 = 0.5 * a1 ** 2
    re = (a1 * math.sin(2.0 * rl * _F) + a2 * math.sin((rl - w) * _F) * math.cos(2.0 * rl * _F)
          - a3 * math.sin(4.0 * rl * _F)) / _F
    avs = re - ec
    return decli, avs


def _spherical_azimuth(p: float, a: float, b: float) -> Tuple[float, float]:
    """Port exact de `angleorient` (orient_paleomag.f:489-505) : resout le
    triangle spherique Pole-Zenith-Astre (angle au sommet `p` = angle
    horaire, cotes `a`/`b` = codeclinaison/colatitude) pour la distance
    zenithale `del` et l'azimuth `ang` de l'astre (ici le soleil)."""
    cal = math.cos(a * _F) * math.cos(b * _F)
    sal = math.sin(a * _F) * math.sin(b * _F) * math.cos(p * _F)
    cs = cal + sal
    cc = min(cs * cs, 1.0)
    y = math.sqrt(1.0 - cc)
    dele = math.atan2(y, cs)
    delta = math.sin(dele) * math.sin(b * _F)
    if delta == 0.0:
        delta = 0.000001
    cs = (math.cos(a * _F) - (math.cos(dele) * math.cos(b * _F))) / delta
    cc = min(cs * cs, 1.0)
    y = math.sqrt(1.0 - cc)
    ang = math.degrees(math.atan2(y, cs))
    return math.degrees(dele), ang


def igrf_declination(lat: float, rlong: float, year: int, month: int, day: float) -> float:
    """Declinaison IGRF (pmagpy.pmag.doigrf, voir docstring module) a la
    position/date donnee - meme decimalisation de la date que le Fortran
    (`year+(month+day/31.0)/12.0`, "pour simplifier on prend 31 jours par
    mois"), altitude 0 (le Fortran d'origine n'en passait pas non plus a
    `igrfosx`)."""
    date_decimal = year + (month + day / 31.0) / 12.0
    x, y, _z, _f = pmag.doigrf(rlong, lat, 0.0, date_decimal)
    return math.degrees(math.atan2(y, x))


def compute_local_declination(
    lat: float, rlong: float, year: int, month: int, day: float,
    hour: float, minute: float, ioutil: int, azmag: float, azsun: float,
) -> Tuple[float, float]:
    """Port exact de `orient_sample` (orient_paleomag.f:596-640).

    `ioutil` : 1 = "aiguille verticale sur platine tournante" (outil
    d'orientation commencant par 'A', ex. "A12_0_3_90"), 2 = "equerre
    pivotante sur platine fixe" (tout le reste) - meme convention que
    `infoech` (dataselect2.f:1029-1030 : `ioutil=2 ; if
    outilorient(1:1)=="A") ioutil=1`), determine laquelle des 2 formules
    `azgeo` s'applique.

    Retourne (decli_igrf, declin) - decli_igrf est le modele global
    (comparaison/verification uniquement), declin est la declinaison
    LOCALE reellement mesuree sur le terrain (voir docstring module) - 0.0
    si aucune visee solaire n'a ete faite (azsun/hour/minute tous nuls,
    meme convention que le Fortran)."""
    decli_igrf = igrf_declination(lat, rlong, year, month, day)

    day_frac = day + (hour + (minute / 60.0)) / 24.0 - 0.5
    dele, av = _sun_declination_and_eot(year, month, day_frac)

    p = rlong + ((hour - 12) + (minute / 60.0)) * 15.0 + av
    a = 90.0 - dele
    b = 90.0 - lat
    if p > 360.0:
        p -= 360.0
    _dell, az = _spherical_azimuth(p, a, b)
    if 0.0 <= p < 180.0:
        az = 360.0 - az
    az += 180.0
    if az > 360.0:
        az -= 360.0

    if ioutil == 2:
        azgeo = 360.0 - azsun + az
    else:
        azgeo = azsun + az
    if azgeo >= 360.0:
        azgeo -= 360.0
    if azgeo < 0.0:
        azgeo += 360.0

    declin = azgeo - azmag
    if declin < -90.0:
        declin += 360.0
    if azsun == 0.0 and hour == 0.0 and minute == 0.0:
        declin = 0.0

    return decli_igrf, declin
