"""
Port de `exportpmagren` ("export detailed Rennes") et `exporttolatex`
("export Latex"), toutes deux dans dataselect.f. Dans le Fortran,
`exportpmagren` appelle `exporttolatex` automatiquement a la fin (un seul
menu declenche les deux exports a la suite) ; ici elles sont deliberement
SEPAREES en deux fonctions/menus independants, comme demande.

Hors perimetre (documente, pas un oubli) :
- Declinaison IGRF (`decli_igrf`/`declin`, via `orient_sample`/
  IGRFstarmac.f) : IGRF n'est pas porte ailleurs dans ce projet
  (selection.py le note deja) - ecrit "n.d" comme le fait le Fortran
  lui-meme pour ses propres cas de donnees manquantes.
- Le bloc `includegraphics{zijder-<id>.pdf}` de exporttolatex (suppose
  des PDF de Zijderveld deja generes sur disque par un export SVG/PDF
  separe, hors du perimetre de cette fonction) : non reproduit.

Les champs Site/Sample/Fm/Age/GC/SMT/Li/Loc viennent directement des
champs magic_* deja decodes depuis la ligne roche (testlect.decode_roche)
- pas besoin de la reparser ici.
"""

from typing import Dict, List, Optional

import math

from selection import SelectedSample, Measurement, polere, corfor, corpen
from calcul import FitResult

_HEADER_COMMENT = [
    "! Paleomagnetic laboratory - Geosciences Rennes",
    "! Equipments 2G magnetometer - Molspin spinner - Agico JR6",
    "! 2G one measurement between two zeros  : code C1",
    "! 2G four measurements between two zeros : code C4",
    "! Molspin spinner 6 positions : code Mo",
    "! Jr6a spinner automated mode : code Ja",
    "! Jr6/Jr5 spinner 2 positions : code J2",
    "! MMTD Furnace D+ and D- = sample orientation changed in the furnace between two steps along +Z and -Z",
    "! AF 2G Online three axis F+ = sequence coils  X,Z,Y, F- = sequence coils  Y,Z,X",
    "! FX = along X; FY = along Y, FZ = along Z;  FG = combine FX,FY and FZ to remove GRM",
    "! Paleointensities with the Coe/Tauxe or Thellier method",
    "! R: indicates field along Z axis; V: indicates field along Z axis; P: Ptrm check",
]


def _sample_display_name(ech: SelectedSample) -> str:
    """Nom d'echantillon = specimen id sans le dernier caractere A-E s'il
    en a un (convention specimen = sample + lettre de sous-carotte)."""
    sid = ech.id.strip()
    if sid and sid[-1] in "ABCDE":
        return sid[:-1]
    return sid


def _dc_field_string(m: Measurement, prev: List[Measurement], rfield: float, thellier: bool) -> str:
    """Equivalent du bloc de 8 `if` construisant `dc_field` (dataselect.f,
    juste avant l'ecriture de chaque ligne de mesure)."""
    def fmt(val: float, template_small: str, template_big: str) -> str:
        return (template_small if abs(val) < 100.0 else template_big) % int(val)

    if m.cod1 == "R":
        return fmt(rfield, "  0:0:%d ", " 0:0:%d")
    if m.cod1 == "V":
        return fmt(-rfield, " 0:0:%d ", " 0:0:%d")
    if m.cod1 == "P":
        if thellier:
            return fmt(rfield, "  0:0:%d ", " 0:0:%d")
        if prev and prev[-1].cod1 == "S":
            return fmt(rfield, "  0:0:%d ", " 0:0:%d")
        return "  0:0:0  "
    if m.cod1 == "Z" and m.cod2 == "+":
        return fmt(rfield, "  0:0:%d ", " 0:0:%d")
    if m.cod1 == "Z" and m.cod2 == "-":
        return fmt(-rfield, " 0:0:%d ", " 0:0:%d")
    if m.cod1 == "Y" and m.cod2 == "+":
        return fmt(rfield, "  0:%d:0 ", " 0:%d:0 ")
    if m.cod1 == "Y" and m.cod2 == "-":
        return fmt(-rfield, " 0:%d:0 ", " 0:%d:0")
    if m.cod1 == "X" and m.cod2 == "+":
        return fmt(rfield, " %d:0:0  ", " %d:0:0 ")
    if m.cod1 == "X" and m.cod2 == "-":
        return fmt(-rfield, " %d:0:0 ", "%d:0:0 ")
    return "  0:0:0  "


def _sample_param_lines(ech: SelectedSample, depthsam: Optional[float]) -> List[str]:
    """Bloc de parametres (str2..str20 du Fortran)."""
    lines = [" --------------  Parameters sample & data   ---------------- ", ""]
    lines.append(f"Site     :  {ech.magic_site}")
    lines.append(f"Sample   :  {_sample_display_name(ech)}")
    lines.append(f"Specimen :  {ech.id}")
    if ech.norme == "v":
        lines.append(f"Volume   :{ech.vol:8.3f}     masse :   n.d ")
    else:
        lines.append(f"Volume   :   n.d      masse :  {ech.vol:8.3f}")

    if depthsam is not None:
        lines.append(f"Lat :{ech.lat:10.5f}    Long :{ech.rlong:12.5f}"
                      f"  height magnetostratigraphy :{depthsam:7.1f}")
    else:
        lines.append(f"Lat :{ech.lat:10.5f}    Long :{ech.rlong:12.5f}"
                      f"     Elevation :{ech.altitude:7.1f}")

    lines.append(f"Sampling date     :  {int(ech.year):4d}   {int(ech.month):2d}   {int(ech.day):2d}")
    lines.append(f"Sampling time UTM :  {int(ech.hour):2d}  {int(ech.minute):2d}")

    lines.append(f"azimuth mag :{ech.azmag:6.1f}  IGRF  Declination :   n.d")
    if ech.azsun == 0.0 and ech.hour == 0.0:
        lines.append(f"azimuth sun :{ech.azsun:6.1f}  Local Declination :    n.d ")
    else:
        lines.append(f"azimuth sun :{ech.azsun:6.1f}  Local Declination :   n.d")

    lines.append('Orientation :  "use AGICO code A12_0_3_90"')
    lines.append(f"core azimuth   :  {ech.caz:6.1f}")
    lines.append(f"core dip       :  {ech.cin:6.1f}")
    lines.append(f"Strike bedding :  {ech.str_:6.1f}")
    lines.append(f"Dip bedding    :  {ech.dip:6.1f}")

    lines.append(f'Formation   :   "{ech.magic_fm}"')
    lines.append(f'Age         :   "{ech.magic_age}"')
    lines.append(f'Geology     :   "{ech.magic_gc} : {ech.magic_smt} : {ech.magic_li}"')
    lines.append(f'Locality    :   "{ech.magic_loc}"')
    lines.append(f'Observation :   "{ech.magic_obs}"')

    rfield = 0.0
    try:
        rfield = float(ech.com[:2])
    except (ValueError, TypeError):
        pass
    if rfield != 0.0:
        lines.append(f"dc applied magnetic field : {rfield:4.1f}  µT")
    else:
        lines.append(f"dc applied magnetic field : {rfield:4.1f}  µT")
    return lines


def _measurement_table_lines(ech: SelectedSample) -> List[str]:
    """Equivalent du tableau `Step code dc_field Mag ... Dsc Isc Dis Iis
    Dtc Itc q ins K` (dataselect.f)."""
    if ech.norme == "m":
        lines = [" Step code  dc_Field    Mag(Am2)     Am2/kg    Dsc   Isc     "
                 "Dis   Iis     Dtc   Itc    q  Mag     K"]
    else:
        lines = [" Step code  dc_field    Mag(Am2)      A/m      Dsc   Isc     "
                 "Dis   Iis     Dtc   Itc    q  Mag     K"]

    rfield = 0.0
    try:
        rfield = float(ech.com[:2])
    except (ValueError, TypeError):
        pass
    thellier = any(m.cod1 == "V" for m in ech.mesures)

    for j, m in enumerate(ech.mesures):
        dc_field = _dc_field_string(m, ech.mesures[:j], rfield, thellier)

        mag1, dec1, inc1 = polere(m.x, m.y, m.z)
        x2, y2, z2 = corfor(m.x, m.y, m.z, ech.cin, ech.caz)
        _mag2, dec2, inc2 = polere(x2, y2, z2)
        x3, y3, z3 = corpen(x2, y2, z2, ech.dip, ech.str_)
        _mag3, dec3, inc3 = polere(x3, y3, z3)

        step_mag = m.etape / 10.0 if m.cod1 == "F" else float(m.etape)
        if ech.norme == "m":
            rxx = mag1 * 1.0e3 / ech.vol if ech.vol else 0.0
            k = m.s * 1.0e-7 / ech.vol if ech.vol else 0.0
        else:
            rxx = mag1 * 1.0e6 / ech.vol if ech.vol else 0.0
            k = m.s * 10.0 * 1.0e-5 / ech.vol if ech.vol else 0.0

        lines.append(
            f"{step_mag:6.1f} {m.cod1}{m.cod2}  {dc_field}  {mag1:10.3E}  {rxx:10.3E}"
            f"  {dec1:6.1f}{inc1:6.1f}  {dec2:6.1f}{inc2:6.1f}  {dec3:6.1f}{inc3:6.1f}"
            f"{m.q:4d}  {m.ins:<2s}  {k:10.3E}"
        )
    return lines


def _fit_result_lines(ech: SelectedSample, results: List[FitResult]) -> List[str]:
    """Equivalent de `lisreslatex` : pour chaque FitResult de cet
    echantillon, une ligne D_is/I_is (in-situ) et D_tc/I_tc (apres
    correction de pendage) - utilise UNIQUEMENT par l'export Latex, pas
    par l'export texte detaille (fidele au Fortran)."""
    matches = [r for r in results if r.id.strip() == ech.id.strip()]
    if not matches:
        return []
    lines = ["", "---- ChRM Best line  or best plane ----", "",
             "   sample    L_P   ori   D  Ncomp  TC    D_is  I_is   D_tc  I_tc    mad     T1    T2"]
    for r in matches:
        incr, decr = r.inc, r.dec
        x = math.cos(math.radians(incr)) * math.cos(math.radians(decr))
        y = math.cos(math.radians(incr)) * math.sin(math.radians(decr))
        z = math.sin(math.radians(incr))
        x2, y2, z2 = corfor(x, y, z, r.cin, r.caz)
        _m2, dec2, inc2 = polere(x2, y2, z2)
        x3, y3, z3 = corpen(x2, y2, z2, r.dip, r.str_)
        _m3, dec3, inc3 = polere(x3, y3, z3)
        if r.demag == "F":
            etapmin, etapmax = r.step_first // 10, r.step_last // 10
        else:
            etapmin, etapmax = r.step_first, r.step_last
        lines.append(
            f"   {r.id:<12s} {r.cat1}   {r.orig}     {r.demag:<3s}  {r.numcomp:1d}  "
            f"{r.nb:3d}   {dec2:5.1f} {inc2:5.1f}   {dec3:5.1f} {inc3:5.1f}  "
            f"{r.mad:4.1f}  {etapmin:5d} {etapmax:5d}"
        )
    return lines


# ---------------------------------------------------------------------------
# export detailed Rennes (exportpmagren)
# ---------------------------------------------------------------------------

def export_detailed_txt(
    samples: List[SelectedSample],
    location: str,
    out_path: str,
    heights: Optional[Dict[str, float]] = None,
) -> None:
    lines = list(_HEADER_COMMENT)
    lines += ["", "", f"Location :{location}", ""]

    for ech in samples:
        if len(ech.mesures) < 2:
            continue
        depthsam = heights.get(ech.id.strip()) if heights else None
        lines.append("")
        lines.append(chr(12) + " --------------  Parameters sample & data   ---------------- ")
        lines.append("")
        lines.extend(_sample_param_lines(ech, depthsam)[1:])  # sans re-repeter le titre
        lines.append("")
        lines.extend(_measurement_table_lines(ech))

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# export Latex (exporttolatex)
# ---------------------------------------------------------------------------

_LATEX_PREAMBLE = r"""\documentclass[a4paper,14pt, titlepage, twoside]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[francais]{babel}
\usepackage{amsmath}
\usepackage{amssymb,amsfonts,textcomp}
\usepackage{color}
\usepackage[colorlinks=true,linkcolor=blue,urlcolor=blue,bookmarks=true]{hyperref}
\usepackage{multicol}
\usepackage{graphicx}
\usepackage[table]{xcolor}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage{datetime}
\usepackage{lastpage}
\usepackage{alltt}
\usepackage{fancyhdr}
\usepackage{adjustbox}
\usepackage{fancyvrb}
\usepackage{geometry}
\geometry{hmargin=1.5cm,vmargin=1.5cm}
\pagestyle{fancy}
\renewcommand{\footrulewidth}{1pt}
\fancyfoot[R]{\small{page~\thepage~sur~\pageref{LastPage}}}
\fancyfoot[c]{\small{Paleomagnetism laboratory -- Geosciences Rennes INSU-CNRS Univ. Rennes1}}
\renewcommand{\headrulewidth}{1pt}
\fancyhead[C]{}
\fancyhead[L]{}
\fancyhead[R]{}
\usepackage{fancybox}
"""


def _latex_escape(s: str) -> str:
    for a, b in (("\\", r"\textbackslash{}"), ("_", r"\_"), ("&", r"\&"),
                 ("%", r"\%"), ("#", r"\#"), ("$", r"\$")):
        s = s.replace(a, b)
    return s


def export_latex(
    samples: List[SelectedSample],
    location: str,
    out_path: str,
    results: Optional[List[FitResult]] = None,
    heights: Optional[Dict[str, float]] = None,
) -> None:
    results = results or []
    parts = [_LATEX_PREAMBLE, r"\begin{document}", ""]

    parts.append(r"\fontsize{28}{28}")
    parts.append(r"\begin{center}\selectfont{")
    parts.append("")
    parts.append(f"Country and study area : {_latex_escape(location)}")
    parts.append("}")
    parts.append(r"\end{center}")

    parts.append(r"\fontsize{8}{10}")
    parts.append(r"\Large List of samples \\")

    sitetest = "xxxxxx"
    for ech in samples:
        if len(ech.mesures) < 2:
            continue
        site6 = ech.id[:6]
        if site6 != sitetest:
            parts.append(r"\\")
            parts.append(rf"\Large Site: {_latex_escape(site6)} \large \\")
            sitetest = site6
        parts.append(rf"\hyperlink{{{ech.id}}}{{ {ech.id}    }}")

    parts.append(r"\fontsize{10}{12}")
    parts.append(r"\begin{verbatim}")
    parts.extend(_HEADER_COMMENT)
    parts.append(r"\end{verbatim}")

    for ech in samples:
        if len(ech.mesures) < 2:
            continue
        depthsam = heights.get(ech.id.strip()) if heights else None
        parts.append(r"\pagebreak")
        parts.append(rf"\hypertarget{{{ech.id}}}{{ }}")
        parts.append(r"\begin{verbatim}")
        parts.extend(_sample_param_lines(ech, depthsam))
        parts.append("")
        parts.extend(_measurement_table_lines(ech))
        parts.extend(_fit_result_lines(ech, results))
        parts.append(r"\end{verbatim}")

    parts.append(r"\end{document}")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(parts) + "\n")
