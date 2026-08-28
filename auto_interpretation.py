"""
Interpretation AUTOMATIQUE (suggestion) de diagrammes de desaimantation -
demande explicite utilisateur ("une routine intelligente qui ferait des
interpretations des diagrammes de desaimantation"), deuxieme des deux
routines demandees (la premiere, l'evaluation d'interpretations deja
faites, est dans interpretation_quality.py).

Principe (choisi explicitement plutot qu'un modele IA/ML - discute avec
l'utilisateur : "je recommande une approche heuristique... plus
verifiable, coherent avec l'esprit du projet") : PCA en extension gloutonne
(greedy), reutilisant EXCLUSIVEMENT calcul.linear_fit (aucune nouvelle
formule statistique) :

1) Composante PRINCIPALE (primary) : part des DERNIERS points (haute
   temperature/champ fort - la ou la composante caracteristique est
   generalement isolee) et etend la fenetre VERS L'ARRIERE tant que le
   MAD reste acceptable et n'augmente pas brusquement (signe qu'une autre
   composante a ete atteinte). Ancree a l'origine par defaut (hypothese
   standard pour une composante caracteristique qui doit decroitre vers
   zero). Si l'ancrage echoue DES le dernier point (magnetisation
   parasite en toute fin de sequence - demande explicite utilisateur :
   "often, the last steps at high temperature is a spurious
   magnetization that departs from zero... it is best to discard the
   last points and keep the vector anchored to the origin"), retire
   d'abord 1 a 3 points du cote haute temperature et reessaie ANCRE
   (voir _primary_search) - repli en libre seulement si meme cela
   echoue.
2) Composante SECONDAIRE (secondary) : meme extension gloutonne mais
   depuis le PREMIER point, VERS L'AVANT, limitee aux points non deja
   couverts par la composante principale. Libre par defaut (une
   surimpression secondaire ne vise pas necessairement l'origine), repli
   en ancree si le libre echoue d'emblee. Absente si aucun point ne
   reste disponible.

Convention "ancre = primary" (demande explicite utilisateur : "by
default, if the component is anchored to the origin, use primary") : si
la recherche secondaire (avant) finit ancree alors que la principale
(arriere) ne l'est pas, les etiquettes sont echangees - l'ancrage a
l'origine prime sur le sens de recherche pour decider quelle composante
est "primary".

Les memes points que le Zijderveld affiche (zijderveld.draw_zijderveld) -
runs IRM/thermique-d'IRM exclus (selection.split_experiments/
experiment_kind) - pour que la composante suggeree corresponde a ce que
l'utilisateur voit reellement dans ce diagramme.

Notation reutilisee de interpretation_quality.py (grade_mad, seuils
MAD/fraction NRM) pour une echelle de confiance COHERENTE entre "evaluer
une interpretation existante" et "noter une suggestion automatique".

LIMITE ASSUMEE (avertissement donne des la premiere discussion, a
rappeler a l'utilisateur) : fonctionne bien sur une decroissance propre a
une seule composante par segment ; reste faillible sur des composantes
qui se chevauchent (transition progressive plutot que franche). Une
suggestion est TOUJOURS a valider par l'utilisateur, jamais archivee
automatiquement - `propose_components` ne fait AUCUNE ecriture dans
self.results/.pmagres, c'est a l'appelant (app.py) de proposer
l'archivage via les mecanismes existants (fit_line + archivres) si
l'utilisateur accepte une suggestion.
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional

from calcul import FitResult, linear_fit, planar_fit
from interpretation_quality import (
    grade_mad, _nrm_fraction_grade, _GRADE_ORDER, _LINE_MAD_THRESHOLD,
    _PLANE_MAD_THRESHOLD, _vector_diff,
)
from selection import zijderveld_measurements, polere

MIN_POINTS = 4
MAX_MAD_JUMP = 5.0  # deg - au-dela, on considere avoir atteint une autre composante


@dataclass
class ComponentSuggestion:
    label: str            # "primary" / "secondary"
    step_first: float
    step_last: float
    anchored: bool
    dec: float
    inc: float
    mad: float
    nb: int
    nrm_fraction: Optional[float]
    grade: str
    notes: List[str] = field(default_factory=list)
    kind: str = "line"    # "line" ou "plane" (grand cercle) - voir _primary_search


_zijderveld_measurements = zijderveld_measurements  # voir selection.py (partage avec interpretation_quality.py)


def _nrm_fraction(mesures: list, start: int, end: int) -> Optional[float]:
    """Meme calcul VECTORIEL (VDS, voir interpretation_quality._nrm_fraction
    et _vector_diff) que le module d'evaluation - notation reutilisee pour
    une echelle de confiance COHERENTE entre suggestion automatique et
    interpretation deja archivee (voir docstring module). Reimplemente ici
    (plutot que d'appeler directement interpretation_quality._nrm_fraction)
    sur la liste de mesures deja filtree (indices 0-based [start,end]
    inclus) plutot que sur ech.mesures/des indices 1-based - evite
    d'exiger que la selection Zijderveld soit un prefixe strict de
    ech.mesures (vrai aujourd'hui mais pas une garantie a figer ici)."""
    if not mesures or len(mesures) < 2:
        return None
    vds_total = sum(_vector_diff(mesures[i], mesures[i + 1]) for i in range(len(mesures) - 1))
    m_last = mesures[-1]
    vds_total += math.sqrt(m_last.x ** 2 + m_last.y ** 2 + m_last.z ** 2)
    if vds_total == 0:
        return None
    interval_vd = sum(_vector_diff(mesures[i], mesures[i + 1]) for i in range(start, end))
    return interval_vd / vds_total


def _greedy_backward(points_xyz: list, anchored: bool, min_start: int = 0) -> Optional[dict]:
    """Meilleure fenetre [start, n-1] trouvee en etendant vers l'arriere
    depuis les derniers MIN_POINTS points - voir docstring module, etape
    1. `min_start` : ne descend pas en dessous de cet indice (voir
    propose_components - exclut le NRM brut, index 0). Retourne None si
    meme la fenetre minimale echoue."""
    n = len(points_xyz)
    if n - min_start < MIN_POINTS:
        return None
    best = None
    prev_mad = None
    for start in range(n - MIN_POINTS, min_start - 1, -1):
        fit = linear_fit(points_xyz[start:n], anchored=anchored, mad_threshold=180.0)
        if fit is None:
            break
        mad = fit["mad"]
        if mad > _LINE_MAD_THRESHOLD:
            break
        if prev_mad is not None and mad - prev_mad > MAX_MAD_JUMP:
            break
        best = (start, n - 1, fit)
        prev_mad = mad
    return best


def _greedy_backward_plane(points_xyz: list, min_start: int = 0) -> Optional[tuple]:
    """Variante PLAN de _greedy_backward - demande explicite utilisateur,
    apres verification sur des donnees reelles (site 14NQ04, fichier
    Tibet) que _primary_search (LIGNE seule) manque ou force un mauvais
    ajustement sur des specimens ou l'interpretation manuelle utilisait un
    plan (grand cercle) pour la composante haute temperature/champ fort -
    cas frequent quand cette composante n'est pas isolee (chevauchement
    avec une composante non totalement retiree). `planar_fit` n'a pas de
    notion d'ancrage (toujours "libre" au sens geometrique - un plan
    contraint une direction sans en choisir une seule ; convention
    Fortran/calcul.fit_plane : orig='o' pour un plan, quoi qu'il en soit -
    voir _make_plane_suggestion) et rejette lui-meme un MAD>25 (voir
    calcul.planar_fit) - pas besoin de reverifier le seuil ici comme pour
    une ligne."""
    n = len(points_xyz)
    if n - min_start < MIN_POINTS:
        return None
    best = None
    prev_mad = None
    for start in range(n - MIN_POINTS, min_start - 1, -1):
        fit = planar_fit(points_xyz[start:n], normalize=True)
        if fit is None:
            break
        mad = fit["mad"]
        if prev_mad is not None and mad - prev_mad > MAX_MAD_JUMP:
            break
        best = (start, n - 1, fit)
        prev_mad = mad
    return best


def _greedy_forward(points_xyz: list, anchored: bool, min_start: int, limit: int) -> Optional[dict]:
    """Symetrique de _greedy_backward : etend vers l'avant depuis les
    premiers MIN_POINTS points a partir de `min_start` (voir
    propose_components - exclut le NRM brut, index 0), sans depasser
    l'indice `limit` (exclu - les points deja pris par la composante
    principale)."""
    if limit - min_start < MIN_POINTS:
        return None
    best = None
    prev_mad = None
    for end in range(min_start + MIN_POINTS - 1, limit):
        fit = linear_fit(points_xyz[min_start:end + 1], anchored=anchored, mad_threshold=180.0)
        if fit is None:
            break
        mad = fit["mad"]
        if mad > _LINE_MAD_THRESHOLD:
            break
        if prev_mad is not None and mad - prev_mad > MAX_MAD_JUMP:
            break
        best = (min_start, end, fit)
        prev_mad = mad
    return best


def _primary_search(points_xyz: list, min_start: int, max_end_trim: int = 3):
    """Recherche de la composante PRINCIPALE - demande explicite
    utilisateur ("often, the last steps at high temperature is a spurious
    magnetization that departs from zero, if the main component was
    going to the origin, it is best to discard the last points and keep
    the vector anchored to the origin") : essaie D'ABORD un fit ANCRE se
    terminant au tout DERNIER point ; si ce dernier point (ou les
    derniers) est une magnetisation parasite qui fait echouer l'ancrage
    des le depart (MAD/linearite), retire 1 a `max_end_trim` points du
    cote haute temperature/champ fort et reessaie ANCRE avant de se
    replier sur un fit LIBRE - preference EXPLICITE pour "ancre mais plus
    court" plutot que "libre mais complet".

    PUIS compare au meilleur ajustement de PLAN (grand cercle) sur la
    meme fenetre haute temperature/champ fort - demande explicite
    utilisateur, apres verification sur des specimens reels (site 14NQ04,
    Tibet) ou l'interpretation manuelle utilisait un plan pour cette
    composante alors que la recherche LIGNE seule echouait entierement ou
    ne trouvait qu'un ajustement mediocre. Le plan n'est retenu que
    lorsqu'il fait clairement mieux qu'une ligne inexistante ou de qualite
    "marginal"/"poor" (voir grade_mad) - une ligne "good"/"excellent"
    reste toujours preferee (une direction unique est plus utile qu'un
    grand cercle des qu'elle est bien resolue). PORTEE VOLONTAIREMENT
    LIMITEE au niveau du specimen individuel - la combinaison de plusieurs
    plans/lignes de PLUSIEURS specimens en une direction de site est un
    probleme distinct et non trivial (l'intersection de grands cercles
    peut converger sur la direction d'une composante secondaire bien
    groupee plutot que sur la primaire - demande explicite utilisateur,
    laisse pour un travail separe), pas traite ici.

    Retourne (start, end, fit, anchored, n_trimmed, kind) ou None -
    kind = "line" ou "plane"."""
    n = len(points_xyz)
    line_result = None
    for end_trim in range(0, max_end_trim + 1):
        end = n - 1 - end_trim
        if end - min_start + 1 < MIN_POINTS:
            break
        candidate = _greedy_backward(points_xyz[:end + 1], anchored=True, min_start=min_start)
        if candidate is not None:
            start, cend, fit = candidate
            line_result = (start, cend, fit, True, end_trim)
            break
    if line_result is None:
        candidate = _greedy_backward(points_xyz, anchored=False, min_start=min_start)
        if candidate is not None:
            start, end, fit = candidate
            line_result = (start, end, fit, False, 0)

    line_grade = grade_mad(line_result[2]["mad"], _LINE_MAD_THRESHOLD) if line_result is not None else None
    if line_result is not None and line_grade in ("excellent", "good"):
        start, end, fit, anchored, n_trimmed = line_result
        return start, end, fit, anchored, n_trimmed, "line"

    plane_result = _greedy_backward_plane(points_xyz, min_start=min_start)
    if plane_result is not None:
        start, end, fit = plane_result
        plane_grade = grade_mad(fit["mad"], _PLANE_MAD_THRESHOLD)
        if plane_grade in ("excellent", "good"):
            return start, end, fit, True, 0, "plane"

    if line_result is not None:
        start, end, fit, anchored, n_trimmed = line_result
        return start, end, fit, anchored, n_trimmed, "line"
    if plane_result is not None:
        start, end, fit = plane_result
        return start, end, fit, True, 0, "plane"
    return None


def _make_suggestion(
    label: str, mesures: list, start: int, end: int, anchored: bool, fit: dict,
    n_trimmed: int = 0, kind: str = "line",
) -> ComponentSuggestion:
    """`fit` est le dict retourne par calcul.linear_fit (cle "direction")
    pour kind="line", ou calcul.planar_fit (cle "pole") pour kind="plane" -
    demande explicite utilisateur ("check the logic implemented in the
    manual interpretation of site 14NQ04") : dec/inc d'un plan sont ceux
    du POLE du grand cercle (meme convention que calcul.fit_plane), PAS
    une direction caracteristique directement comparable a celle d'une
    ligne - a afficher/interpreter comme tel."""
    dx, dy, dz = fit["pole"] if kind == "plane" else fit["direction"]
    _, dec, inc = polere(dx, dy, dz)
    nrm_fraction = _nrm_fraction(mesures, start, end)
    threshold = _PLANE_MAD_THRESHOLD if kind == "plane" else _LINE_MAD_THRESHOLD
    grade = grade_mad(fit["mad"], threshold)
    notes = []
    if kind == "plane":
        notes.append("plane fit (great circle) - dec/inc shown are the POLE, not a characteristic direction")
    if n_trimmed:
        notes.append(
            f"discarded {n_trimmed} trailing high-treatment point(s) that departed from "
            f"the origin-ward trend, to keep the fit anchored"
        )
    frac_grade = _nrm_fraction_grade(nrm_fraction)
    if frac_grade is not None and _GRADE_ORDER[frac_grade] < _GRADE_ORDER.get(grade, 99):
        notes.append(f"only {nrm_fraction * 100:.0f}% of initial NRM described by this interval")
        grade = frac_grade
    if fit["nb"] < 4:
        notes.append(f"only {fit['nb']} points - statistically fragile")
    return ComponentSuggestion(
        label=label, step_first=mesures[start].etape, step_last=mesures[end].etape,
        anchored=anchored, dec=dec, inc=inc, mad=fit["mad"], nb=fit["nb"],
        nrm_fraction=nrm_fraction, grade=grade, notes=notes, kind=kind,
    )


def propose_components(ech) -> List[ComponentSuggestion]:
    """Propose 0, 1 ou 2 composantes (primary, et secondary si des points
    restent disponibles avant elle) pour ce specimen. Ne modifie rien,
    n'ecrit rien - voir docstring module.

    Le NRM BRUT (TOUTE PREMIERE mesure de la sequence, non traitee) est
    TOUJOURS exclu des deux recherches, INCONDITIONNELLEMENT (pas
    seulement si cod1=='N' - demande explicite utilisateur : "by default,
    it is best to not include the NRM0... start the interpretation at
    best at the first demag step") - convention paleomagnetique standard,
    PAS un critere statistique : verifie sur donnees reelles
    (old_pmag.ren) que le MAD change a peine (0.54 -> 0.45 deg) selon
    qu'on l'inclue ou non, alors que l'interpretation humaine l'exclut
    systematiquement (le NRM brut porte souvent une composante visqueuse
    non representative, retiree des le premier palier de traitement,
    meme s'il "tombe" pres de la droite ajustee par coincidence
    statistique)."""
    mesures = _zijderveld_measurements(ech)
    if len(mesures) < MIN_POINTS:
        return []
    points_xyz = [(m.x, m.y, m.z) for m in mesures]
    min_start = 1 if len(mesures) > MIN_POINTS else 0

    primary_result = _primary_search(points_xyz, min_start=min_start)

    primary_start = len(mesures)
    primary_suggestion = None
    if primary_result is not None:
        start, end, fit, anchored, n_trimmed, kind = primary_result
        primary_start = start
        primary_suggestion = _make_suggestion("primary", mesures, start, end, anchored, fit, n_trimmed, kind)

    secondary_fit = _greedy_forward(points_xyz, anchored=False, min_start=min_start, limit=primary_start)
    secondary_anchored = False
    if secondary_fit is None:
        secondary_fit = _greedy_forward(points_xyz, anchored=True, min_start=min_start, limit=primary_start)
        secondary_anchored = True
    secondary_suggestion = None
    if secondary_fit is not None:
        start, end, fit = secondary_fit
        secondary_suggestion = _make_suggestion("secondary", mesures, start, end, secondary_anchored, fit)

    # "by default, if the component is anchored to the origin, use
    # primary" (demande explicite utilisateur) : une composante ancree a
    # l'origine EST par convention la caracteristique/terminale - si la
    # recherche "secondaire" (avant) finit ancree alors que la
    # "principale" (arriere) ne l'est pas, on echange simplement les
    # etiquettes plutot que le sens de recherche (start/end restent
    # corrects, seul le label et le numcomp implicite changent).
    suggestions: List[ComponentSuggestion] = [s for s in (primary_suggestion, secondary_suggestion) if s is not None]
    if (
        primary_suggestion is not None and secondary_suggestion is not None
        and secondary_suggestion.anchored and not primary_suggestion.anchored
    ):
        primary_suggestion.label, secondary_suggestion.label = "secondary", "primary"
        suggestions = [secondary_suggestion, primary_suggestion]

    return suggestions


def format_suggestions(ech_id: str, suggestions: List[ComponentSuggestion]) -> str:
    """Rendu compact - une ligne courte par suggestion (demande explicite
    utilisateur : "the line is too long and difficult to read"), dec/inc
    omis ici (visibles directement sur le Zijderveld superpose)."""
    if not suggestions:
        return f"{ech_id}: no component could be automatically proposed (too few usable points, or too noisy)."
    lines = [f"{ech_id} - auto-interpretation suggestions:"]
    for s in suggestions:
        nrm_txt = f"  NRM={s.nrm_fraction * 100:.0f}%" if s.nrm_fraction is not None else ""
        anc = "plane" if s.kind == "plane" else ("anc" if s.anchored else "free")
        lines.append(
            f"  {s.label:<9s} {s.step_first:.0f}-{s.step_last:.0f}  {anc}  n={s.nb}"
            f"  MAD={s.mad:.1f}{nrm_txt}  [{s.grade}]"
        )
        for note in s.notes:
            lines.append(f"      -> {note}")
    return "\n".join(lines)
