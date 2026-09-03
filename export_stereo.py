"""
Export de resultats Starmac vers des formats DIRECTEMENT compatibles avec
StereoUtils_Py - demande explicite utilisateur ("is it possible to add the
export results to Stereo_Py as in the original Fortran with in addition a
specific file for the poles", puis correction : "in the export to Stereo_Py
I did not explain that it was with a format compatible to the project in
Stereo_Py").

Deux formats DIFFERENTS, pour deux menus DIFFERENTS de StereoUtils_Py (pas
le meme fichier ni la meme convention de colonnes) :

- `export_results_to_stereo` : produit un fichier au format "Project"
  (menu Project > Load Project..., voir `stereo_project.py`/
  `stereo_project.load_project` dans ce depot) - lignes POSITIONNELLES (PAS
  d'en-tete, PAS d'alias de colonne) :
      layer  id  type(d/m/g)  dec  inc  alpha95  iplot(2/3)  symbol  rgb  size
  UNE entree 'd' par resultat specimen (cat1 in L/P/f/s), groupee par
  LAYER = site (6 premiers caracteres du specimen, meme convention que
  magic_export._site_mean_row/convert_ren_to_r) - permet de recalculer une
  moyenne de Fisher PAR SITE directement dans StereoUtils_Py via son propre
  menu "Fisher Project" (qui moyenne les entrees consecutives d'un meme
  layer). UNE entree 'm' par moyenne de site ("mean:", cat1=='F'), dans un
  layer SEPARE ("Site_Means", jamais le layer du site correspondant) : les
  melanger romprait silencieusement un recalcul "Fisher Project" ulterieur
  (`fisher_project` ne filtre PAS par type - une moyenne deja calculee se
  retrouverait comptee comme un point de donnees supplementaire, biaisant
  le resultat). Une seule orientation par fichier (le format "Project" n'a
  pas de place pour une 2e paire dec/inc TC comme le Fortran d'origine
  `exportres` - voir parametre `orientation`).

- `export_poles_to_stereo` : PAS le format "Project" (poles = coordonnees
  GLOBALES lat/lon d'un VGP, pas des directions locales - la stereonet
  "Project" n'a pas de sens pour un pole). Format a en-tete auto-descriptif
  du menu Data input (colonnes reperees par ALIAS, voir
  `stereo_selection.FIELD_ALIASES`), colonnes "dec"/"inc" portant en
  realite VGP_lon/VGP_lat - EXACTEMENT la convention utilisee par
  `plot_vgps_on_map` (app.py de StereoUtils_Py : "les dec/inc en memoire
  representent directement longitude/latitude du VGP, meme convention que
  le Fortran") : charger ce fichier via le menu "Data input" standard,
  PUIS "Plot VGPs on Map" (Pmag_Python), fonctionne directement, sans
  passer par le format "Project". `site`/`nb`/`a95`/`dp`/`dm`/`p95`
  restent disponibles en colonnes supplementaires (non-alias) pour
  reference."""

from typing import List

from calcul import FitResult, _correct_dec_inc, dp_dm_from_a95


def export_results_to_stereo(results: List[FitResult], out_path: str, orientation: int = 2) -> int:
    """Ecrit un fichier au format "Project" de StereoUtils_Py (voir
    docstring module). `orientation` : 2=in-situ (par defaut, meme
    convention que StarmacApp.orientation), 3=apres pendage complet -
    UNE SEULE des deux, contrairement a l'export "exportres" Fortran
    d'origine qui ecrivait toujours les deux (le format "Project" n'a
    qu'une paire dec/inc par entree). Retourne le nombre de lignes
    ecrites (0 si aucun resultat exploitable)."""
    lines = []
    for r in results:
        if r.id[:5] == "mean:" or r.cat1 not in ("L", "P", "f", "s"):
            continue
        site = r.id[:6]
        dec, inc = _correct_dec_inc(r, orientation)
        lines.append(
            f"{site:<14s} {r.id:<12s} d  {dec:6.1f}  {inc:6.1f}  {0.0:6.1f}  "
            f" 3  c  0_0_0       {0.30:6.2f}"
        )
    for r in results:
        if r.id[:5] != "mean:":
            continue
        site = r.id[6:].strip()
        lines.append(
            f"{'Site_Means':<14s} {site:<12s} m  {r.dec:6.1f}  {r.inc:6.1f}  {r.mad:6.1f}  "
            f" 3  e  0_0_0       {0.55:6.2f}"
        )
    if not lines:
        return 0
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return len(lines)


def export_poles_to_stereo(results: List[FitResult], out_path: str) -> int:
    """Ecrit un fichier "dec"/"inc" (VGP_lon/VGP_lat) directement lisible
    par "Data input" -> "Plot VGPs on Map" de StereoUtils_Py (voir
    docstring module - PAS le format "Project"). Retourne le nombre de
    lignes ecrites (0 si `results` ne contient aucune moyenne de site
    "mean:")."""
    rows = []
    for r in results:
        if r.id[:5] != "mean:":
            continue
        site = r.id[6:].strip()
        vgp_dp, vgp_dm = (
            (r.vgp_dp, r.vgp_dm) if (r.vgp_dp or r.vgp_dm)
            else dp_dm_from_a95(r.mad, r.inc)
        )
        p95 = (vgp_dp + vgp_dm) / 2.0
        rows.append([
            f"{r.par5:.1f}", f"{r.par4:.1f}",  # dec=VGP_lon, inc=VGP_lat
            site, str(r.nb), f"{r.mad:.1f}",
            f"{vgp_dp:.1f}", f"{vgp_dm:.1f}", f"{p95:.1f}",
        ])
    if not rows:
        return 0
    header = ["dec", "inc", "site", "nb", "a95", "dp", "dm", "p95"]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("#" + "\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    return len(rows)
