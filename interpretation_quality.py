"""
Evalue la qualite d'interpretations DEJA enregistrees (.pmagres, resultats
de type 'L'/'P') - demande explicite utilisateur ("une routine qui
permettrait d'evaluer les interpretations"), premiere des deux routines
demandees (l'autre, l'interpretation automatique par PCA en fenetre
glissante, est un chantier separe, pas traite ici).

Principe : reutilise EXCLUSIVEMENT les fonctions de fit deja validees
(calcul.linear_fit/fit_line/fit_plane) - aucune nouvelle formule
statistique introduite. La valeur ajoutee vient de calculer, pour chaque
resultat, des diagnostics que l'enregistrement seul ne porte pas :

- MAD RECALCULE en direct depuis les mesures vivantes (donnees), pas
  relu du fichier - detecte une interpretation devenue "perimee" si les
  mesures ont change depuis l'archivage (ex. une mesure corrigee/
  supprimee) sans que le resultat ait ete refait.
- Pour une droite (cat1='L') : DANG-like - angle entre le fit ANCRE et le
  fit LIBRE sur le MEME intervalle de pas. Standard en paleointensite
  (voir paleointensity.fit_arai_direction, meme principe) mais PAS
  habituellement calcule pour une interpretation directionnelle
  ordinaire - ajoute ici. Ne degrade la note QUE si le fit archive est
  lui-meme ANCRE (res.orig=='o') et que son alternatif libre diverge
  beaucoup (l'ancrage force est peut-etre inapproprie) - demande
  explicite utilisateur ("when the component is not anchored to the
  origin, the component is OK and might be a true component ; then the
  DANG will be very high even for an excellent component") : pour un fit
  deja LIBRE, un grand DANG signifie seulement que la composante ne passe
  pas par l'origine (attendu pour une surimpression/composante secondaire
  authentique), reste une simple note informative, jamais une degradation.
- linearity_ratio (calcul.linear_fit) : deja calcule en interne pour
  rejeter un fit manifestement courbe (seuil 0.5), mais jusqu'ici jamais
  expose - remonte ici comme diagnostic continu plutot que pass/fail.
- nrm_fraction : proportion du VDS (vector difference sum, meme
  convention que pmag.dovds - voir _nrm_fraction) qu'occupe l'intervalle
  interprete - demande explicite utilisateur ("la proportion de NRM prise
  en compte dans l'interpretation"). Seuils RECALIBRES (demande explicite
  utilisateur, apres le passage a une base VECTORIELLE ci-dessous qui
  donne des valeurs structurellement plus basses qu'une simple difference
  de magnitude) : >=40% excellent, 15-40% bon, 5-15% limite, <5% pauvre -
  4 paliers au lieu de 3 (ajout de "poor", deja utilise ailleurs pour le
  MAD, voir _GRADE_ORDER). Distingue une composante VRAIMENT dominante/caracteristique
  d'une composante mineure ne representant qu'une petite fraction du
  signal total, meme avec un excellent MAD (un MAD parfait sur 5% du NRM
  ne dit rien de la fiabilite de l'interpretation prise dans son
  ensemble). Base sur la difference VECTORIELLE entre pas consecutifs
  (pas la simple difference de magnitude des extremites) - demande
  explicite utilisateur ("this site does show well-defined secondary
  components with normal polarity and a primary magnetization with
  reverse polarity. You need to test the intensity of magnetization as
  vectorial difference but with the absolute intensity of the
  magnetization") : un specimen multi-composantes peut voir la magnitude
  BRUTE du vecteur total croitre puis decroitre de facon non monotone
  (une composante secondaire disparait pendant qu'une primaire de
  direction differente domine de plus en plus le vecteur somme) - la
  premiere version (difference des magnitudes aux extremites) pouvait
  alors donner une fraction NEGATIVE pour une composante haute-
  temperature parfaitement valide.

Nettoyage prealable des mesures (demande explicite utilisateur : "be sure
to remove the GRM (option 1) automatically, remove the IRM or thermal
demag of IRM") - voir _cleaned_measurements :
- runs IRM/thermique-d'IRM exclus, memes points que le Zijderveld affiche
  (reutilise auto_interpretation._zijderveld_measurements).
- `datatools.eliminate_grm(method=1)` (= `elimineGRM_DZ`, deja porte,
  jusqu'ici seulement accessible via le menu manuel "Suppress GRM")
  applique sur une COPIE JETABLE - jamais sur donnees/ech.mesures, qui
  restent inchanges - et SEULEMENT si le specimen porte reellement un
  triplet GRM-tumbler (cod1='F', cod2 'X'/'Y'/'Z' consecutifs) : sinon
  cette fonction (fidele au Fortran) supprimerait a tort les pas 'F'
  ISOLES d'un AF classique en mode tumbler simple (cod2='T', jamais un
  triplet) - comportement correct pour son usage prevu (echantillons
  GRM), destructeur hors de ce cas.

Ponderation ("a ponderer avec les autres parametres", demande explicite
utilisateur) : la note finale est LA PLUS FAIBLE des deux notes
(MAD-based, nrm_fraction-based) - un maillon faible suffit a degrader la
confiance globale, meme si l'autre critere est excellent.

Seuils de note (MAD) : REPRENNENT ceux deja utilises par l'appli au
moment du fit (calcul.py : 15 deg pour une droite/ajuslig, 25 pour un
plan/ajusplans) - un resultat enregistre a donc TOUJOURS ete <= ce seuil
AU MOMENT DE L'ARCHIVAGE ; la note ici (excellent/bon/limite) subdivise
CET INTERVALLE DEJA ACCEPTE pour distinguer les bons fits des fits
tout juste passables, et "perime" signale un MAD qui a depuis DEPASSE
ce seuil (donnees modifiees). Le seuil DANG (15 deg) est un choix de
CE module, pas une convention MagIC/Fisher etablie - a ajuster si
l'usage montre qu'il est trop/pas assez strict.
"""

import copy
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from calcul import FitResult, linear_fit, planar_fit
from selection import zijderveld_measurements
from datatools import eliminate_grm

_LINE_MAD_THRESHOLD = 15.0
_PLANE_MAD_THRESHOLD = 25.0
_DANG_THRESHOLD = 15.0  # voir docstring module : choix de ce module, pas une convention etablie
# Seuils NRM_fraction (FRAC vectoriel, voir _nrm_fraction) - recalibres par
# l'utilisateur apres le passage a une base vectorielle (valeurs
# structurellement plus basses qu'avec l'ancienne difference de magnitude).
_NRM_FRACTION_EXCELLENT = 0.40
_NRM_FRACTION_GOOD = 0.15
_NRM_FRACTION_MARGINAL = 0.05


@dataclass
class QualityReport:
    id: str
    cat1: str
    numcomp: int
    orig: str
    step_first: float
    step_last: float
    nb: int
    mad_saved: float          # tel que lu dans le .pmagres
    mad_live: Optional[float]  # recalcule depuis les mesures vivantes (None si non recalculable)
    linearity_ratio: Optional[float]
    dang: Optional[float]      # droites seulement - angle ancre/libre
    nrm_fraction: Optional[float]  # proportion du NRM initial decrite par l'intervalle
    grade: str                 # "excellent" / "good" / "marginal" / "stale" / "unrecomputable"
    notes: List[str] = field(default_factory=list)


_GRADE_ORDER = {"stale": 0, "unrecomputable": 1, "poor": 1, "marginal": 2, "good": 3, "excellent": 4}


def _cleaned_ech(ech):
    """Copie JETABLE de `ech` (voir docstring module) avec `.mesures`
    remplace par la sequence Zijderveld (IRM/thermique-d'IRM exclus) puis,
    si le specimen porte reellement un triplet GRM-tumbler, apres passage
    de `datatools.eliminate_grm(method=1)`. `ech`/`donnees` d'origine ne
    sont jamais modifies - seule cette copie l'est."""
    mesures = zijderveld_measurements(ech)
    has_grm_triplet = any(m.cod1 == "F" and m.cod2 == "X" for m in mesures)
    proxy = copy.copy(ech)
    proxy.mesures = list(mesures)
    if has_grm_triplet:
        eliminate_grm(proxy, method=1)
    return proxy


def _find_jdeb_jfin(ech, res: FitResult) -> Optional[tuple]:
    """Meme recherche que calcul.recompute_fit_geometry (indices 1-based
    par correspondance sur `etape`)."""
    jdeb = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == res.step_first), None)
    jfin = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == res.step_last), None)
    if jdeb is None or jfin is None or jfin < jdeb:
        return None
    return jdeb, jfin


def _vector_diff(a, b) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _nrm_fraction(ech, jdeb: int, jfin: int) -> Optional[float]:
    """Proportion du VDS (vector difference sum, meme convention que
    pmag.dovds - somme des |Mi-Mi+1| consecutifs sur TOUTE la sequence de
    l'echantillon, plus la magnitude du DERNIER point comme "difference a
    zero" implicite) qu'occupe l'intervalle [jdeb,jfin].

    PORTE ICI depuis une premiere version basee sur la simple DIFFERENCE
    DE MAGNITUDE des deux extremites ((|M(jdeb)|-|M(jfin)|)/|M(0)|) -
    demande explicite utilisateur ("this site does show well-defined
    secondary components with normal polarity and a primary magnetization
    with reverse polarity. You need to test the intensity of
    magnetization as vectorial difference but with the absolute intensity
    of the magnetization") : quand un specimen porte plusieurs composantes
    de directions differentes (cas normal, pas une anomalie), la magnitude
    BRUTE du vecteur total peut croitre puis decroitre de facon non
    monotone au fil de la desaimantation (une composante secondaire dans
    une direction disparait pendant qu'une composante primaire dans une
    AUTRE direction reste, voire domine de plus en plus le vecteur
    somme) - la version precedente pouvait alors donner une fraction
    NEGATIVE pour une composante haute-temperature parfaitement valide,
    simplement parce que |M| a la fin de son intervalle depassait |M| au
    debut. La difference VECTORIELLE (norme du vecteur DIFFERENCE entre
    deux pas consecutifs, somme le long de l'intervalle) mesure combien de
    signal a reellement ete retire par CETTE composante, independamment de
    ce qui restait d'une autre composante - toujours positive, jamais de
    signe contre-intuitif. None si le specimen a moins de 2 mesures ou si
    le VDS total est nul (specimen sans signal)."""
    mesures = ech.mesures
    if not mesures or len(mesures) < 2:
        return None
    vds_total = sum(_vector_diff(mesures[i], mesures[i + 1]) for i in range(len(mesures) - 1))
    m_last = mesures[-1]
    vds_total += math.sqrt(m_last.x ** 2 + m_last.y ** 2 + m_last.z ** 2)
    if vds_total == 0:
        return None
    interval_vd = sum(_vector_diff(mesures[i], mesures[i + 1]) for i in range(jdeb - 1, jfin - 1))
    return interval_vd / vds_total


def grade_mad(mad: float, threshold: float) -> str:
    """excellent/good/marginal/poor selon MAD - meme subdivision que
    `evaluate_result` (seuil/3, 2*seuil/3, seuil), extraite ici pour etre
    reutilisee par auto_interpretation.py (pas de notion de "stale"/
    "unrecomputable" ici : une suggestion n'a pas de valeur sauvegardee a
    comparer, contrairement a un resultat archive)."""
    if mad > threshold:
        return "poor"
    if mad > threshold * 2 / 3:
        return "marginal"
    if mad > threshold / 3:
        return "good"
    return "excellent"


def _nrm_fraction_grade(fraction: Optional[float]) -> Optional[str]:
    if fraction is None:
        return None
    if fraction >= _NRM_FRACTION_EXCELLENT:
        return "excellent"
    if fraction >= _NRM_FRACTION_GOOD:
        return "good"
    if fraction >= _NRM_FRACTION_MARGINAL:
        return "marginal"
    return "poor"


def evaluate_result(res: FitResult, donnees) -> Optional[QualityReport]:
    """None si le resultat n'est pas de type 'L'/'P' (ex. une moyenne de
    site "mean:") ou si le specimen/l'intervalle ne sont plus retrouvables
    (specimen supprime depuis, .prmag rechargee sans ce specimen...)."""
    if res.cat1 not in ("L", "P"):
        return None
    ech = next((s for s in donnees if s.id == res.id), None)
    if ech is None or not getattr(ech, "mesures", None):
        return None
    ech = _cleaned_ech(ech)  # voir docstring module - copie jetable, IRM/GRM traites

    notes: List[str] = []
    ji = _find_jdeb_jfin(ech, res)
    if ji is None:
        return QualityReport(
            id=res.id, cat1=res.cat1, numcomp=res.numcomp, orig=res.orig,
            step_first=res.step_first, step_last=res.step_last, nb=res.nb,
            mad_saved=res.mad, mad_live=None, linearity_ratio=None, dang=None,
            nrm_fraction=None, grade="unrecomputable",
            notes=["step_first/step_last no longer found in the specimen's measurements"],
        )
    jdeb, jfin = ji
    points = [(m.x, m.y, m.z) for m in ech.mesures[jdeb - 1:jfin]]
    nrm_fraction = _nrm_fraction(ech, jdeb, jfin)

    if res.cat1 == "L":
        anchored = res.orig == "o"
        live = linear_fit(points, anchored=anchored, mad_threshold=180.0)  # pas de rejet ici, on veut TOUJOURS le chiffre
        other = linear_fit(points, anchored=not anchored, mad_threshold=180.0)
        mad_live = live["mad"] if live else None
        linearity_ratio = live["linearity_ratio"] if live else None
        dang = None
        if live is not None and other is not None:
            d1, d2 = live["direction"], other["direction"]
            cosang = max(-1.0, min(1.0, float(sum(a * b for a, b in zip(d1, d2)))))
            dang = math.degrees(math.acos(cosang))
        elif other is None:
            notes.append(f"complementary {'free' if anchored else 'anchored'} fit failed - no DANG-like check available")
        threshold = _LINE_MAD_THRESHOLD
    else:  # 'P'
        live = planar_fit(points, normalize=True)
        mad_live = live["mad"] if live else None
        linearity_ratio = None
        dang = None
        threshold = _PLANE_MAD_THRESHOLD

    if mad_live is None:
        grade = "unrecomputable"
        notes.append("fit could not be recomputed from live measurements (degenerate points?)")
    elif mad_live > threshold:
        grade = "stale"
        notes.append(f"MAD now {mad_live:.1f} > {threshold:.0f} (was {res.mad:.1f} at save time) - measurements changed since archiving?")
    else:
        grade = grade_mad(mad_live, threshold)

    if dang is not None and dang > _DANG_THRESHOLD:
        # Demande explicite utilisateur ("when the component is not
        # anchored to the origin, the component is OK and might be a true
        # component ; then the DANG will be very high even for an
        # excellent component") : DANG mesure l'ecart entre le fit
        # ARCHIVE et son alternatif complementaire, PAS un defaut du fit
        # archive lui-meme. Pour un fit LIBRE (res.orig != 'o'), un grand
        # DANG signifie seulement que la droite ne pointe pas vers
        # l'origine - normal pour une composante secondaire/surimpression
        # authentique qui ne decroit pas a zero, pas un signe de mauvaise
        # qualite. Seul un fit ANCRE dont l'alternatif libre diverge
        # beaucoup est reellement suspect (l'ancrage force peut-etre un
        # passage par l'origine que les donnees ne soutiennent pas) - la
        # degradation de note ne s'applique donc qu'a ce cas.
        if res.orig == "o":
            notes.append(f"anchored/free directions disagree by {dang:.1f} deg - anchoring to the origin may not be appropriate here")
            if grade in ("excellent", "good"):
                grade = "marginal"
        else:
            notes.append(f"anchored/free directions disagree by {dang:.1f} deg - expected for a free (non-anchored) component, not a quality issue")

    if res.nb < 4:
        notes.append(f"only {res.nb} points - statistically fragile even with good MAD")

    # Ponderation demandee par l'utilisateur ("a ponderer avec les autres
    # parametres") : la note finale retient LA PLUS FAIBLE des deux notes
    # (MAD-based ci-dessus, nrm_fraction-based) - un maillon faible degrade
    # la confiance globale meme si l'autre critere est bon. N'intervient
    # pas sur un grade deja "stale"/"unrecomputable" (probleme de donnees,
    # pas de qualite de fit - une bonne fraction de NRM ne compense pas).
    frac_grade = _nrm_fraction_grade(nrm_fraction)
    if frac_grade is not None and grade in ("excellent", "good", "marginal"):
        if _GRADE_ORDER[frac_grade] < _GRADE_ORDER[grade]:
            notes.append(
                f"only {nrm_fraction * 100:.0f}% of initial NRM described by this interval "
                f"(MAD alone would grade this '{grade}')"
            )
            grade = frac_grade
    elif nrm_fraction is None:
        notes.append("could not compute NRM fraction (zero-intensity first measurement?)")

    return QualityReport(
        id=res.id, cat1=res.cat1, numcomp=res.numcomp, orig=res.orig,
        step_first=res.step_first, step_last=res.step_last, nb=res.nb,
        mad_saved=res.mad, mad_live=mad_live, linearity_ratio=linearity_ratio,
        dang=dang, nrm_fraction=nrm_fraction, grade=grade, notes=notes,
    )


def evaluate_results(results: List[FitResult], donnees) -> List[QualityReport]:
    reports = []
    for res in results:
        r = evaluate_result(res, donnees)
        if r is not None:
            reports.append(r)
    return reports


def format_quality_report(reports: List[QualityReport], worst_first: bool = True) -> str:
    if not reports:
        return "No line/plane result could be evaluated (empty list, or all are site means)."
    ordered = sorted(reports, key=lambda r: _GRADE_ORDER.get(r.grade, 99)) if worst_first else reports

    counts: Dict[str, int] = {}
    for r in reports:
        counts[r.grade] = counts.get(r.grade, 0) + 1
    summary = ", ".join(f"{g}: {counts[g]}" for g in ("stale", "unrecomputable", "marginal", "good", "excellent") if g in counts)

    lines = [f"Interpretation quality audit - {len(reports)} result(s) - {summary}", ""]
    for r in ordered:
        dang_txt = f"  DANG(anc/free)={r.dang:.1f}" if r.dang is not None else ""
        lin_txt = f"  lin_ratio={r.linearity_ratio:.2f}" if r.linearity_ratio is not None else ""
        nrm_txt = f"  NRM_frac={r.nrm_fraction * 100:.0f}%" if r.nrm_fraction is not None else ""
        mad_txt = f"MAD={r.mad_live:.1f}" if r.mad_live is not None else "MAD=n/a"
        lines.append(
            f"[{r.grade:>14s}] {r.id:<13s} {r.cat1}{r.numcomp} {r.orig}  "
            f"steps {r.step_first:.0f}-{r.step_last:.0f}  n={r.nb}  {mad_txt} (saved {r.mad_saved:.1f})"
            f"{nrm_txt}{lin_txt}{dang_txt}"
        )
        for note in r.notes:
            lines.append(f"                 -> {note}")
    return "\n".join(lines)
