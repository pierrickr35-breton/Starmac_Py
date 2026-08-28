import io
import os
import subprocess
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from testlect import read_ren_file, read_prmag_file
from convert_legacy_ren import convert_legacy_auto
from complete_sample_info import complete_site_info, complete_specimen_info
from import_new_data import parse_jr6_file, parse_legacy_new_measurements, archive_new_measurements
from convert_ren_to_r import convert_file as convert_ren_to_r_file
from convert_magic_to_r import convert_magic_file
from convert_utrecht_to_r import convert_files as convert_utrecht_files
from irm import build_irm_figure, has_irm_data
from selection import (
    select_samples,
    select_samples_by_site,
    delete_measurements,
    init_selection,
    list_measurements,
    list_xyz,
    list_measurements_vrm,
    list_measurements_depth,
    sample_info,
)
from interpretation_quality import evaluate_result, evaluate_results, format_quality_report
from auto_interpretation import propose_components, format_suggestions
from calcul import (
    FitResult,
    fit_line,
    fit_plane,
    fit_fisher_direction,
    fit_single_direction,
    fit_lines_auto,
    fit_from_redo_file,
    fisher_from_measurements,
    fisher_from_results,
    list_results,
    init_results,
    results_path_for,
    archivres,
    load_results,
    recompute_fit_geometry,
    available_mean_orientations,
    list_mdf,
    compute_mean_intensity,
    compute_mean_susceptibility,
    format_mean_intensity,
    compute_mean_inclination,
    format_mean_inclination,
    compute_koenigsberger,
    format_koenigsberger,
    list_diff_measurements,
    apply_viscosity_test,
    apply_subtraction,
    record_arm_holder,
    detect_cooling_rate_rows,
    compute_cooling_rate,
    format_cooling_rate,
    read_ani_tensor,
    apply_inverse_anisotropy,
    compute_anisotropy_tensor,
    write_ani_tensor,
    compute_anicor_factor,
    _ANI_CODE2,
    _correct_dec_inc,
)
from zijderveld import build_zijderveld_figure, draw_zijderveld
from stereo import build_stereo_figure, build_stereo_results_figure
from xygraph import build_xygraph_figure, has_mixed_demag
from susceptibility import build_susceptibility_figure
from svgwriter import SVGWriter
from paleointensity import (
    compute_arai,
    fit_arai_line,
    fit_arai_direction,
    parse_com_field,
    detect_method_and_hlab,
    compute_crm,
    arai_curvature,
    build_arai_figure,
    build_paleoint_review_figure,
    draw_arai,
)
from paleointensity_magic import compute_magic_paleointensity, format_magic_paleointensity
from datatools import (
    convert_thellier_to_nrm,
    remove_step,
    eliminate_grm,
    convert_z_minus,
    export_thellier_tdt,
)
from magic_export import export_to_magic
from detailed_export import export_detailed_txt, export_latex

SVG_DEBUG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "svg_debug")

# Raccourcis clavier repris de starmac_OSX.f95 (AWE_addMenu(...,"CTRL+...")),
# adaptes suite a un diagnostic reel (diag_clavier.py) sur clavier AZERTY :
# Option (ALT) est le modificateur d'ACCENT sur les claviers europeens -
# Cmd+Option+Shift+E genere le keysym "Ecircumflex" (Ê), jamais "E". Tout
# raccourci combinant Option + une lettre est donc casse par construction
# sur AZERTY (confirme, pas specifique aux chiffres comme suppose au debut).
# On evite donc completement Option : CTRL+ALT+x (Fortran) est traduit ici
# en Cmd+Ctrl+<lettre> plutot qu'en Cmd+Option+Shift+<lettre>, en choisissant
# des lettres qui evitent les raccourcis systeme macOS connus (Cmd+Ctrl+Q =
# verrouiller l'ecran, Cmd+Ctrl+F = plein ecran).
#   - "CTRL+l" (lismes)          -> Cmd+L        (CTRL -> Command)
#   - "CTRL+META+z" (plotzijder) -> Cmd+Ctrl+Z    (META -> Control)
#   - "CTRL+ALT+1" (selce)       -> Cmd+Ctrl+A    (ALT -> Control aussi, pas Option)
# Note de syntaxe Tk (piege verifie empiriquement) : des que Shift est un
# modificateur explicite, Tk exige la lettre du keysym en MAJUSCULE pour que
# le binding se declenche - non applicable ici puisqu'on n'utilise plus Shift
# du tout dans les combinaisons a 3 modificateurs.
# valeur : (accelerateur affiche dans le menu, sequence de bind Tk)
SHORTCUTS = {
    "importpc":   ("Cmd+O",      "<Command-o>"),
    "starend":    ("Cmd+Ctrl+X", "<Command-Control-x>"),
    "selmes":     ("Cmd+E",      "<Command-e>"),
    "selentete":  ("Cmd+Ctrl+E", "<Command-Control-e>"),
    "effmes":     ("Cmd+D",      "<Command-d>"),
    "initmes":    ("Cmd+I",      "<Command-i>"),
    "lismes":     ("Cmd+L",      "<Command-l>"),
    "infoech":    ("Cmd+J",      "<Command-j>"),
    "selce":      ("Cmd+Ctrl+A", "<Command-Control-a>"),
    "selis":      ("Cmd+Ctrl+B", "<Command-Control-b>"),
    "selcp":      ("Cmd+Ctrl+T", "<Command-Control-t>"),
    "ajuslig":    ("Cmd+B",      "<Command-b>"),
    "fishmes":    ("Cmd+F",      "<Command-f>"),
    "fishres":    ("Cmd+Ctrl+F", "<Command-Control-f>"),
    "lisres":     ("Cmd+Ctrl+L", "<Command-Control-l>"),
    "selres":     ("Cmd+Ctrl+R", "<Command-Control-r>"),
    "initres":    ("Cmd+Ctrl+I", "<Command-Control-i>"),
    "plotzijder": ("Cmd+Ctrl+Z", "<Command-Control-z>"),
    "xygraph":    ("Cmd+Ctrl+Y", "<Command-Control-y>"),
}

ORIENTATION_SHORTCUT_NAMES = {1: "selce", 2: "selis", 3: "selcp"}

ORIENTATIONS = {
    "Sample (CE)": 1,
    "In situ (IS)": 2,
    "Tilt cor. (CP)": 3,
}


class StarmacApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Starmac_AWE_4.0 - Paleomagnetism")
        self.root.geometry("1100x700")

        self.donnees = []  # Stockage des données d'échantillons (List[Pmag])
        self.selection = []  # Dernière sélection (List[SelectedSample])
        self.results = []  # Ajustements de droite (List[FitResult], equivalent tr/nbres)
        self.results_path = None  # equivalent filr - fichier .r, fixe au chargement des donnees
        self._archived_ids = None  # cache des `c` deja utilises dans results_path (voir _save_result)
        self._arm_holder_background = None  # equivalent xholarm/yholarm/zholarm (holderarm), pour Anisotropy
        self.entete = ""  # Préfixe de sélection (equivalent selentete)
        self.orientation = tk.IntVar(value=2)  # equivalent iorient (CE/IS/CP) - defaut In situ
        self._current_graphic = None  # ("zijderveld", sample_id) / ("stereo"|"xygraph"|"susceptibility"|"arai", None)
        self._arai_state = None  # (ech, points, checks, arno, fit_ou_None, hlab) pour kind=="arai"
        self._paleoint_review_state = None  # (ech, points, checks, arno, fit) pour kind=="paleoint_review"

        self._setup_menu()
        self._setup_shortcuts()

        self.paned_window = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True)

        # PANNEAU GAUCHE : Zone Graphique (matplotlib embarqué)
        self.graph_frame = ttk.Frame(self.paned_window, width=550)
        self.fig = Figure(figsize=(5.2, 5.2), dpi=100)
        self.canvas_fig = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas_fig.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.paned_window.add(self.graph_frame, weight=1)

        # PANNEAU DROIT : Console / Liste de données
        self.text_frame = ttk.Frame(self.paned_window, width=550)
        self.text_area = tk.Text(
            self.text_frame,
            bg="#ffffff",
            fg="#000000",
            insertbackground="black",
            font=("Courier", 14),
            wrap="none",
        )
        text_yscroll = ttk.Scrollbar(self.text_frame, orient=tk.VERTICAL, command=self.text_area.yview)
        text_xscroll = ttk.Scrollbar(self.text_frame, orient=tk.HORIZONTAL, command=self.text_area.xview)
        self.text_area.configure(yscrollcommand=text_yscroll.set, xscrollcommand=text_xscroll.set)
        self.text_area.tag_configure("prompt", foreground="#c0392b")
        self.text_area.grid(row=0, column=0, sticky="nsew")
        text_yscroll.grid(row=0, column=1, sticky="ns")
        text_xscroll.grid(row=1, column=0, sticky="ew")
        self.text_frame.rowconfigure(0, weight=1)
        self.text_frame.columnconfigure(0, weight=1)
        self.paned_window.add(self.text_frame, weight=1)

        self.load_demo_data()

        # Force l'activation de l'appli au niveau macOS. lift()/focus_force()
        # (niveau Tk) et le bascule topmost (niveau window manager) n'ont pas
        # suffi : la fenetre s'affiche et reagit a la souris mais l'appli ne
        # devient jamais le process "frontmost" pour le clavier tant qu'on
        # n'a pas clique sur la barre de menus. On demande donc directement
        # a System Events d'activer ce process par PID (osascript) - c'est
        # le mecanisme d'activation macOS reel, hors de Tk. Peut demander
        # une autorisation "Automatisation" la premiere fois (a accepter).
        self.root.after(200, self._activate_window)

    def _activate_window(self):
        self.root.attributes("-topmost", True)
        self.root.after(50, lambda: self.root.attributes("-topmost", False))
        try:
            subprocess.run(
                [
                    "osascript", "-e",
                    f'tell application "System Events" to set frontmost of '
                    f'(first process whose unix id is {os.getpid()}) to true',
                ],
                check=False,
                capture_output=True,
                timeout=2,
            )
        except (subprocess.SubprocessError, OSError):
            pass

    @staticmethod
    def _labeled(text, shortcut_name):
        """Libelle avec le raccourci entre parentheses, SANS utiliser le
        parametre `accelerator=` de Tk : sur Aqua, ce dernier semble faire
        que macOS intercepte la combinaison au niveau du menu natif (la
        barre de menus s'active/se surligne) sans jamais invoquer la
        commande Tcl associee - ce qui court-circuite bind_all. Le texte
        seul, lui, est purement cosmetique."""
        return f"{text}    ({SHORTCUTS[shortcut_name][0]})"

    def _setup_menu(self):
        menubar = tk.Menu(self.root)

        # Menu PmagFile (titres repris de StarmacOSX_x.f95, bloc menu anglais
        # lignes 491-608 - AWE_addMenu(HelpMenuUnit,"PmagFile"/"Pmag data"/
        # "Graphics"/"Calcul", ...) - pas le bloc francais "Fichiers"/
        # "Selection donnees"/"Graphiques" utilise plus haut dans le fichier)
        # Reorganise a la demande explicite utilisateur ("est-ce aussi
        # possible de reorganiser le menu Pmag file") - regroupe par
        # fonction plutot qu'en une seule liste plate : ouverture, puis
        # conversions VERS le format .prmag (le format natif courant),
        # puis interoperabilite MagIC (import ET export - aucune de ces
        # commandes ne charge quoi que ce soit dans la session, meme
        # "Import..." : toutes ecrivent un nouveau fichier sur disque),
        # puis exports/rapports.
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label=self._labeled("Open Pmag file...", "importpc"),
                               command=self.ouvrir_fichier_ren)
        file_menu.add_separator()
        file_menu.add_command(label="Import Starmac legacy files...", command=self.ouvrir_convert_legacy_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="Import MagIC contribution to .prmag format...",
                               command=self.ouvrir_convert_magic_to_r_dialog)
        file_menu.add_command(label="Import Utrecht/PMAG2 .col to .prmag format...",
                               command=self.ouvrir_convert_utrecht_to_r_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="Complete sample information...", command=self.ouvrir_complete_sample_info_dialog)
        file_menu.add_command(label="Archive new lab data...", command=self.ouvrir_archive_new_data_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="export to Magic from Starmac_Py", command=self.ouvrir_export_magic_dialog)
        file_menu.add_separator()
        file_menu.add_command(label="export pmag content as a text file", command=self.ouvrir_export_detailed_dialog)
        file_menu.add_command(label="export pmag content as a LaTeX file", command=self.ouvrir_export_latex_dialog)
        file_menu.add_separator()
        file_menu.add_command(label=self._labeled("Quit Starmac", "starend"), command=self.root.quit)
        menubar.add_cascade(label="PmagFile", menu=file_menu)

        # Menu Pmag data
        data_menu = tk.Menu(menubar, tearoff=0)
        data_menu.add_command(label=self._labeled("Select data...", "selmes"),
                               command=self.ouvrir_selection_dialog)
        data_menu.add_command(label="Select site...", command=self.ouvrir_selection_site_dialog)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("Select header...", "selentete"),
                               command=self.ouvrir_entete_dialog)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("Delete some data...", "effmes"),
                               command=self.ouvrir_effmes_dialog)
        data_menu.add_command(label=self._labeled("Init list", "initmes"),
                               command=self.reinitialiser_selection)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("List data", "lismes"), command=self.lister_mesures)
        data_menu.add_command(label="List in XYZ", command=self.lister_xyz)
        data_menu.add_command(label="List data VRM", command=self.lister_vrm)
        data_menu.add_command(label="List and depth...", command=self.ouvrir_lismesdepth_dialog)
        data_menu.add_separator()
        data_menu.add_command(label=self._labeled("Info samples", "infoech"),
                               command=self.afficher_info_echantillons)
        data_menu.add_separator()

        # Radiobuttons directement dans le menu (pas de sous-menu imbrique) :
        # meme raisonnement que ci-dessus, en plus du probleme deja identifie
        # avec les sous-menus en cascade.
        for label, value in ORIENTATIONS.items():
            data_menu.add_radiobutton(
                label=self._labeled(label, ORIENTATION_SHORTCUT_NAMES[value]),
                variable=self.orientation, value=value,
                command=lambda v=value: self._set_orientation(v),
            )

        menubar.add_cascade(label="Pmag data", menu=data_menu)

        # Menu Results (regroupe ce qui etait eparpille entre Pmag
        # data/Calcul/Graphics - demande explicite utilisateur : "adding a
        # menu Results and another Paleointensity"). "data+interpretation"
        # y va (pas dans Paleointensity) - demande explicite utilisateur :
        # "data+interpretation concerne les vecteurs et va dans results".
        results_menu = tk.Menu(menubar, tearoff=0)
        results_menu.add_command(label=self._labeled("Select results...", "selres"),
                                  command=self.ouvrir_selres_dialog)
        results_menu.add_command(label=self._labeled("List results", "lisres"), command=self.lister_resultats)
        results_menu.add_command(label=self._labeled("Init results", "initres"),
                                  command=self.reinitialiser_resultats)
        results_menu.add_separator()
        results_menu.add_command(label=self._labeled("best lines...", "ajuslig"),
                                  command=self.ouvrir_ajuslig_dialog)
        results_menu.add_command(label="best lines auto", command=self.ouvrir_ajusligauto_dialog)
        results_menu.add_command(label="best planes...", command=self.ouvrir_ajusplans_dialog)
        results_menu.add_command(label="best dir Fisher...", command=self.ouvrir_ajusfisher_dialog)
        results_menu.add_command(label="Best fit from redo file...", command=self.ouvrir_ajusligredo_dialog)
        results_menu.add_command(label="Auto-interpret (suggest components)...", command=self.ouvrir_autointerpretation_dialog)
        results_menu.add_separator()
        results_menu.add_command(label="Evaluate interpretations...", command=self.evaluer_interpretations)
        results_menu.add_command(label=self._labeled("Fisher results", "fishres"), command=self.fisher_resultats)
        results_menu.add_separator()
        results_menu.add_command(label="Stereo Results", command=self.afficher_stereo_results)
        results_menu.add_command(label="data+interpretation", command=self.afficher_visres)
        menubar.add_cascade(label="Results", menu=results_menu)

        # Menu Paleointensity (demande explicite utilisateur : "view
        # paleoint va dans le menu paleointensite"). GRM va dans Calcul
        # ("Mettre GRM dans calcul"), Test Induite retire du menu
        # ("supprimer test induite" - fonction conservee dans le code,
        # simplement plus exposee par un intitule de menu).
        paleoint_menu = tk.Menu(menubar, tearoff=0)
        paleoint_menu.add_command(label="Paleointensity...", command=self.afficher_arai)
        paleoint_menu.add_command(label="View Paleoint Results...", command=self.ouvrir_openfilepint_dialog)
        paleoint_menu.add_separator()
        paleoint_menu.add_command(label="Thellier >> NRM", command=self.ouvrir_convertthelli_dialog)
        paleoint_menu.add_command(label="Remove step", command=self.ouvrir_removestep_dialog)
        paleoint_menu.add_command(label="ConvertZ- 2G", command=self.ouvrir_convzmoins_dialog)
        paleoint_menu.add_command(label="export to ThellierTool...", command=self.ouvrir_exportthellier_dialog)
        paleoint_menu.add_separator()
        paleoint_menu.add_command(label="Cooling rate...", command=self.ouvrir_cooling_rate_dialog)
        menubar.add_cascade(label="Paleointensity", menu=paleoint_menu)

        # Menu Calcul (categorie gardee en francais dans la source elle-meme,
        # cf ligne 580 du bloc anglais - seuls les items sont traduits).
        # Items non encore implementes (au-dela d'ajuslig/fisher/resultats) :
        # stubs relies a self._not_implemented, memes intitules que la source,
        # SAUF "Orientation drill cores"/"Orient Core with LowTemp" (lies a
        # l'option Drillcore, retiree de l'appli).
        calcul_menu = tk.Menu(menubar, tearoff=0)
        calcul_menu.add_command(label=self._labeled("Fisher data", "fishmes"), command=self.fisher_mesures)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="MdF-MdT", command=self.afficher_mdf)
        calcul_menu.add_command(label="Mean Intensity", command=self.afficher_mean_intensity)
        calcul_menu.add_command(label="Koenigsberger ratio...", command=self.ouvrir_koenigsberger_dialog)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Mean Inclination", command=self.afficher_mean_inclination)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Test viscosity", command=self.appliquer_viscosity_test)
        calcul_menu.add_command(label="Diff measurements n/n-1", command=self.afficher_diff_measurements)
        calcul_menu.add_command(label="Subtraction...", command=self.ouvrir_subtraction_dialog)
        calcul_menu.add_command(label="Autoinverse", command=lambda: self._not_implemented("Autoinverse"))
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Anisotropy", command=self.ouvrir_anisotropy_dialog)
        calcul_menu.add_command(label="Holder_ARM...", command=self.ouvrir_holderarm_dialog)
        calcul_menu.add_command(label="Inverse_ANI_correction...", command=self.ouvrir_inverseani_dialog)
        calcul_menu.add_separator()
        calcul_menu.add_command(label="Suppress GRM", command=self.ouvrir_elimine_grm_dialog)
        menubar.add_cascade(label="Calcul", menu=calcul_menu)

        # Menu Graphics
        graph_menu = tk.Menu(menubar, tearoff=0)
        graph_menu.add_command(label=self._labeled("Zijderveld", "plotzijder"), command=self.afficher_zijderveld)
        graph_menu.add_command(label="Stereo data", command=self.afficher_stereo)
        graph_menu.add_command(label=self._labeled("XYgraph", "xygraph"), command=self.afficher_xygraph)
        graph_menu.add_command(label="Susceptibility", command=self.afficher_susceptibilite)
        graph_menu.add_command(label="Plot IRM", command=self.afficher_irm)
        graph_menu.add_separator()
        graph_menu.add_command(label="Clear Screen", command=self.clear_screen)
        graph_menu.add_separator()
        graph_menu.add_command(label="Export SVG...", command=self.exporter_svg)
        menubar.add_cascade(label="Graphics", menu=graph_menu)

        self.root.config(menu=menubar)

    def _setup_shortcuts(self):
        """Raccourcis clavier globaux (voir SHORTCUTS). Tous utilisent Command
        comme modificateur de base (jamais Control seul), donc aucun conflit
        avec les raccourcis d'edition Emacs integres a Tk pour les widgets
        Entry/Text (Control+E, Control+D... eux sont en Control seul) : pas
        besoin d'ignorer les raccourcis selon le widget qui a le focus."""
        bindings = {
            "importpc": self.ouvrir_fichier_ren,
            "starend": self.root.quit,
            "selmes": self.ouvrir_selection_dialog,
            "selentete": self.ouvrir_entete_dialog,
            "effmes": self.ouvrir_effmes_dialog,
            "initmes": self.reinitialiser_selection,
            "lismes": self.lister_mesures,
            "infoech": self.afficher_info_echantillons,
            "selce": lambda: self._set_orientation(1),
            "selis": lambda: self._set_orientation(2),
            "selcp": lambda: self._set_orientation(3),
            "ajuslig": self.ouvrir_ajuslig_dialog,
            "fishmes": self.fisher_mesures,
            "fishres": self.fisher_resultats,
            "lisres": self.lister_resultats,
            "selres": self.ouvrir_selres_dialog,
            "initres": self.reinitialiser_resultats,
            "plotzijder": self.afficher_zijderveld,
            "xygraph": self.afficher_xygraph,
        }
        for name, callback in bindings.items():
            _, sequence = SHORTCUTS[name]
            self.root.bind_all(sequence, lambda event, cb=callback: cb())

    def _set_orientation(self, value):
        self.orientation.set(value)
        self._refresh_current_graphic()

    def _write_svg_debug(self, ech, fits_for_sample):
        """Pour test/comparaison : ecrit le rendu SVGWriter (fidele a
        svginit/svgplot) du Zijderveld affiche, AVANT le dessin matplotlib -
        pas encore le pipeline d'affichage definitif, juste un moyen rapide
        de verifier le SVG genere contre un export Starmac_OSX original."""
        try:
            os.makedirs(SVG_DEBUG_DIR, exist_ok=True)
            # page 19x28cm, origine (90,600)px : valeurs reelles de `plots()`
            # (graphicsAWE.f95) - draw_zijderveld part de cette origine de
            # PAGE pour appliquer lui-meme le decalage initial (-2,-5.5) puis
            # l'origine du NEV (u,v), comme le fait `zijder2`.
            writer = SVGWriter(width_cm=19.0, height_cm=28.0)
            writer.set_origin_px(90.0, 600.0)
            draw_zijderveld(writer, ech, orientation=self.orientation.get(), fits=fits_for_sample)
            writer.plotnd()
            path = os.path.join(SVG_DEBUG_DIR, f"zijder-{ech.id}.svg")
            writer.save(path)
        except Exception as e:
            print(f"[svg_debug] echec ecriture SVG pour {ech.id} : {e}")

    def _refresh_current_graphic(self):
        if self._current_graphic is None:
            return
        kind, extra = self._current_graphic
        self._clear_figure()

        if kind == "zijderveld":
            ech = next((s for s in self.selection if s.id == extra), None)
            if ech is None:
                return
            fits_for_sample = [r for r in self.results if r.id == ech.id]
            self._write_svg_debug(ech, fits_for_sample)
            self.fig.set_size_inches(5.5, 8.5, forward=True)
            build_zijderveld_figure(ech, orientation=self.orientation.get(), fits=fits_for_sample, fig=self.fig)
            # liste des donnees de l'echantillon affiche dans la fenetre
            # texte (PAS sur le graphique) - demande explicite de
            # l'utilisateur, meme format que "List data" (Cmd+L).
            buffer = io.StringIO()
            list_measurements([ech], orientation=self.orientation.get(), out=buffer)
            self._afficher(buffer.getvalue())
        elif kind == "stereo":
            if not self.selection:
                return
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_figure(self.selection, orientation=self.orientation.get(), fig=self.fig)
        elif kind == "stereo_results":
            if not self.results:
                return
            self.fig.set_size_inches(5.5, 5.5, forward=True)
            build_stereo_results_figure(
                self.results, orientation=self.orientation.get(),
                nbech=len(self.selection), fig=self.fig,
            )
        elif kind == "xygraph":
            if not self.selection:
                return
            # empile verticalement (pas cote a cote) : le panneau graphique
            # a une largeur fixe, on agrandit la HAUTEUR, pas la largeur.
            height = 8.5 if has_mixed_demag(self.selection) else 4.5
            self.fig.set_size_inches(6.0, height, forward=True)
            build_xygraph_figure(self.selection, fig=self.fig)  # appelle deja fig.tight_layout()
            self._redraw_canvas()
            return
        elif kind == "susceptibility":
            if not self.selection:
                return
            self.fig.set_size_inches(6.0, 4.5, forward=True)
            build_susceptibility_figure(self.selection, fig=self.fig)  # appelle deja fig.tight_layout()
            self._redraw_canvas()
            return
        elif kind == "irm":
            if not self.selection:
                return
            self.fig.set_size_inches(6.0, 8.0, forward=True)
            build_irm_figure(self.selection, fig=self.fig)  # appelle deja fig.tight_layout()
            self._redraw_canvas()
            return
        elif kind == "arai":
            if self._arai_state is None:
                return
            ech, points, checks, arno, fit, _hlab = self._arai_state
            self.fig.set_size_inches(6.0, 8.0, forward=True)
            build_arai_figure(ech, points, checks, arno, fit=fit, fig=self.fig)
        elif kind == "paleoint_review":
            if self._paleoint_review_state is None:
                return
            ech, points, checks, arno, fit = self._paleoint_review_state
            # page Fortran (visi_Paleoint.f: call plots(19.5,28.,fname)) -
            # remplace l'ancien (7.5,9.5) qui ecrasait le figsize deja
            # corrige a l'interieur de build_paleoint_review_figure
            # (celui-ci ne s'applique que si fig=None, or self.fig existe
            # deja ici) - demande explicite utilisateur ("each of the
            # three could be larger within the page").
            self.fig.set_size_inches(19.5 / 2.54, 28.0 / 2.54, forward=True)
            build_paleoint_review_figure(
                ech, points, checks, arno, fit=fit,
                # orientation FORCEE a 1 (echantillon), PAS
                # self.orientation.get() : le Fortran (plotpaleoint2.f et
                # visi_Paleoint.f) force "iorient=1" sans condition en
                # tete de cette routine, quel que soit le reglage
                # d'orientation par ailleurs dans l'appli - demande
                # explicite utilisateur ("the combined plots should be in
                # sample coordinates").
                orientation=1, fig=self.fig,
            )
            # PAS _fit_figure_to_data() : ecrit pour un panneau SIMPLE
            # (ax.axes[0] uniquement) - sur cette figure a 3 panneaux, elle
            # redimensionnerait toute la figure sur le seul rapport
            # largeur/hauteur du panneau Arai (le premier), ecrasant la
            # taille/mise en page deja correctes de build_paleoint_review_
            # figure (page Fortran 19.5x28cm, subplots_adjust) - meme
            # constat que xygraph/susceptibility/irm ci-dessus, qui
            # bypassent deja _fit_figure_to_data pour la meme raison -
            # demande explicite utilisateur ("each of the three could be
            # larger within the page").
            self._redraw_canvas()
            return
        else:
            return

        self._fit_figure_to_data()
        self._redraw_canvas()

    def _redraw_canvas(self):
        """Agrandit le CONTENANT (volet gauche du PanedWindow + fenetre)
        pour accueillir la Figure demandee, plafonne a l'espace ecran
        reellement disponible, force un vrai cycle de redimensionnement -
        puis MESURE la taille reellement obtenue par le widget canvas et
        cale la Figure dessus, plutot que l'inverse.

        Lecon tiree de plusieurs iterations infructueuses : matplotlib
        traite le WIDGET Tk comme source de verite pour la taille de la
        Figure (`FigureCanvasTkAgg` resynchronise la Figure dessus des
        qu'un vrai evenement <Configure> se produit) - essayer d'imposer
        une taille de Figure malgre ca revient a se battre contre le
        framework : la moindre resynchronisation ulterieure (y compris
        celle qu'on declenche nous-meme en "secouant" la fenetre pour
        forcer Tk a se mettre a jour) ecrase la valeur programmee. La bonne
        approche est donc : dimensionner le CONTENANT selon ce qu'on veut
        afficher, puis lire sa taille REELLE une fois stabilisee, et
        aligner la Figure dessus - jamais l'inverse."""
        width_in, height_in = self.fig.get_size_inches()
        dpi = self.fig.dpi

        max_width_px = max(400, self.root.winfo_screenwidth() - 100)
        max_height_px = max(300, self.root.winfo_screenheight() - 100)
        avail_fig_w_in = (max_width_px - 20) / dpi
        avail_fig_h_in = (max_height_px - 90) / dpi
        if width_in > avail_fig_w_in or height_in > avail_fig_h_in:
            scale = min(avail_fig_w_in / width_in, avail_fig_h_in / height_in, 1.0)
            width_in, height_in = width_in * scale, height_in * scale

        width_px = int(width_in * dpi) + 20
        height_px = int(height_in * dpi) + 90  # marge titre/menus/barre d'etat

        try:
            current_w = self.paned_window.sashpos(0)
            if width_px > current_w:
                self.paned_window.sashpos(0, width_px)
        except tk.TclError:
            pass

        if height_px > self.root.winfo_height():
            self.root.geometry(f"{self.root.winfo_width()}x{height_px}")

        # laisse Tk appliquer sashpos()/geometry() avant de lire la taille
        # resultante - sinon winfo_width/height() ci-dessous renverraient
        # encore l'ancienne valeur (la demande est asynchrone).
        self.root.update_idletasks()

        # "secoue" la fenetre pour forcer un vrai cycle <Configure> sur
        # tous les widgets enfants (constate : ni update_idletasks(), ni
        # sashpos()/geometry() seuls ne suffisent - seul un VRAI
        # redimensionnement le fait, d'ou "toucher la fenetre" qui
        # corrige toujours le probleme manuellement).
        w = self.root.winfo_width()
        h = self.root.winfo_height()
        self.root.geometry(f"{w + 1}x{h + 1}")
        self.root.update_idletasks()
        self.root.geometry(f"{w}x{h}")
        self.root.update_idletasks()

        # Aligne la Figure sur la taille REELLE obtenue par le widget
        # (source de verite), pas sur ce qu'on voulait au depart.
        canvas_w = self.canvas_fig.get_tk_widget().winfo_width()
        canvas_h = self.canvas_fig.get_tk_widget().winfo_height()
        if canvas_w > 1 and canvas_h > 1:
            self.fig.set_size_inches(canvas_w / dpi, canvas_h / dpi, forward=False)

        self.canvas_fig.draw()

    def _fit_figure_to_data(self, width_in: float = 5.5, max_height_in: float = 14.0):
        """Cale la taille de la Figure sur le vrai rapport largeur/hauteur des
        donnees tracees (mesure sur `ax.get_xlim()/get_ylim()` apres l'appel a
        `build_*_figure`, qui a deja fait `relim()+autoscale_view()`).

        Necessaire car `PlotContext.clear()` impose `ax.set_aspect('equal')` :
        si la Figure n'a PAS le meme rapport largeur/hauteur que les donnees,
        matplotlib retrecit la zone de tracee (`adjustable='box'`) pour
        respecter l'aspect - le dessin se retrouve reduit a une bande etroite
        au milieu d'une figure trop large/haute, avec de grandes marges vides
        (constate pour le Zijderveld+stereo, dont l'etendue verticale relle
        est bien plus grande que ce qu'une figure 5.5x8.5 fixe peut refleter)."""
        if not self.fig.axes:
            return
        ax = self.fig.axes[0]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        w = xlim[1] - xlim[0]
        h = ylim[1] - ylim[0]
        if w <= 0 or h <= 0:
            return
        height_in = min(width_in * (h / w), max_height_in)
        self.fig.set_size_inches(width_in, height_in, forward=True)
        self.fig.tight_layout()

    # ------------------------------------------------------------------
    # Fichiers
    # ------------------------------------------------------------------

    def ouvrir_fichier_ren(self):
        """Ouvre une boîte de dialogue pour sélectionner un fichier .ren
        ou .prmag et charge ses données.

        Accepte aussi .txt : anciens formats de fichier (ex. Domeyko - une
        seule ligne 'Id:', pas de 'L:'/roche ; Corbières - 'L:' + roche
        collee sans guillemets en fin de ligne plutot que sur une ligne a
        part) que `read_ren_file` sait maintenant lire tout aussi bien.

        .prmag (nouveau format, voir convert_ren_to_r.py/testlect.
        read_prmag_file) est lu par un parseur different mais produit la
        MEME structure List[Pmag] - tout le reste de l'application
        continue de fonctionner sans changement, quel que soit le format
        ouvert - demande explicite utilisateur ("the open file still
        look for .ren")."""
        fichier_path = filedialog.askopenfilename(
            title="Select a .ren, .prmag or .txt data file",
            filetypes=[
                ("REN/PRMAG/TXT files", "*.ren *.prmag *.txt"),
                ("PRMAG files", "*.prmag"),
                ("REN/TXT files", "*.ren *.txt"),
                ("All files", "*.*"),
            ],
        )
        if fichier_path:
            self._load_data_file(fichier_path)

    def _load_data_file(self, fichier_path: str, announce: bool = True) -> bool:
        """Charge `fichier_path` (.ren/.prmag/.txt) dans self.donnees/
        self.results_path et l'affiche - factorise depuis
        ouvrir_fichier_ren pour etre reutilisable a la fin des dialogues
        de conversion (Import Starmac legacy files/MagIC/Utrecht) -
        demande explicite utilisateur ("ce serait bien d'ouvrir les
        fichiers convertis a la fin des conversions"). `announce=False`
        omet le messagebox "Success" (utilise apres une conversion, dont
        le propre messagebox "Conversion complete" a deja informe
        l'utilisateur - eviter 2 popups consecutifs pour la meme action).
        Retourne True si le chargement a reussi."""
        try:
            if fichier_path.lower().endswith(".prmag"):
                self.donnees = read_prmag_file(fichier_path)
            else:
                self.donnees = read_ren_file(fichier_path)
            self.selection = []
            # equivalent filr=fil1(1:(nlen-4))//'.r' : fichier resultats
            # derive du fichier de donnees (meme dossier, extension .r)
            self.results_path = results_path_for(fichier_path)
            self._archived_ids = None  # relit le fichier .r au prochain archivage

            # Réinitialisation et affichage dans la zone de texte
            self.text_area.delete("1.0", tk.END)
            self.text_area.insert(
                tk.END,
                f"File loaded: {os.path.basename(fichier_path)}\n"
                f"Number of samples: {len(self.donnees)}\n"
                "----------------------------------------\n\n"
            )

            for sample in self.donnees:
                self.text_area.insert(
                    tk.END,
                    f"ID: {sample.id} | Measurements: {sample.nbmes} | Inc: {sample.cin} | Az: {sample.caz}\n"
                )

            if announce:
                messagebox.showinfo(
                    "Success", f"{len(self.donnees)} sample(s) loaded!"
                )
            return True
        except Exception as e:
            messagebox.showerror("Error", f"Could not read the file:\n{e}")
            return False

    def ouvrir_convert_legacy_dialog(self):
        """Equivalent GUI de Convert2newpmagformat/convert_oldpmag_2ren.f95 :
        met a niveau un vieux fichier de donnees (1, 2 ou 3 lignes d'entete
        par echantillon, eventuellement MELANGEES dans le meme fichier)
        vers le format actuel a 3 lignes, en completant Site/Sample/Fm/
        Age/GC/SMT/Li/Loc/Obs depuis un fichier "complement" externe
        OPTIONNEL. Un seul point d'entree, `convert_legacy_auto` - demande
        explicite utilisateur ("single import legacy files (able to
        import 1 ligne, 2 lignes and 3 lines header) - the app will
        automatically recognize if there is one, two or three lines") :
        remplace l'ancien choix manuel "Case 1/2/3" (les 3 sous-routines
        dediees convert_case1_specimen_info/convert_case2_site_info/
        convert_case3_new_format restent disponibles directement dans
        convert_legacy_ren.py pour un complement SANS en-tete reconnu,
        seul cas que l'auto-detection ne sait pas lever - voir sa
        docstring). Seul l'ancien "cas 3" (fichier deja a 2 lignes) a pu
        etre verifie octet pres contre un exemple reel fourni par
        l'utilisateur (Convert2newpmagformat/old_pmag.txt + complement.txt
        -> old_pmag.ren).

        Le complement est OPTIONNEL - demande explicite utilisateur
        ("import these files as they are and later we can complete the
        missing information... most of the tasks within the program just
        need the core correction, bedding and volume or mass... most of
        the missing information in the oldest file are for an export to
        Magic") : la ligne `Id:` du vieux fichier porte deja in:/az:/dip:/
        str:/v: (correction de carotte, pendage, volume/masse), jamais
        touchee par le complement - celui-ci ne sert qu'a renseigner
        Site/Fm/Age/GC/SMT/Li/Loc/Obs pour un futur export MagIC. Sans
        complement (ou pour un specimen sans correspondance), la
        conversion procede quand meme, avec "Not Specified" pour les
        champs geologie - a completer plus tard par une routine dediee
        plutot que de bloquer l'import maintenant.

        Enchaine AUSSI la conversion .ren -> .prmag (convert_ren_to_r.
        convert_file) - demande explicite utilisateur ("convert all types
        of legacy files to the .prmag format") : jusqu'ici il fallait
        relancer manuellement "Convert .ren to new format .r..." sur le
        fichier "..._converted.ren" tout juste produit. Le .ren
        intermediaire est CONSERVE sur disque (pas juste un fichier de
        travail jete) : c'est le point d'entree naturel pour completer a
        la main les champs "Not Specified"/les codes non reconnus avant
        de reconvertir en .prmag, sans avoir a repartir du fichier legacy
        d'origine."""
        self.text_area.insert(tk.END, "\n--- Import legacy files (Escape to cancel) ---\n", "prompt")

        old_path = filedialog.askopenfilename(
            title="Select the old data file (1, 2 or 3 header lines)",
            filetypes=[("Text/REN", "*.txt *.ren"), ("All files", "*.*")],
        )
        if not old_path:
            return

        complement_path = None
        if messagebox.askyesno(
            "Complement file",
            "Provide a complement file (Site/Formation/Age/... for a future MagIC export)?\n\n"
            "Not required: core correction, bedding and volume/mass are already in the data "
            "file and will be imported either way.\n\n"
            "The complement file must have a recognized column header (e.g. a first line "
            "'specimen\\tsite\\tsample\\tlat\\tlon\\t...' or 'site\\tformation\\tage\\t...').",
        ):
            complement_path = filedialog.askopenfilename(
                title="Select the complement file",
                filetypes=[("Text", "*.txt"), ("All files", "*.*")],
            )

        base, _ext = os.path.splitext(old_path)
        ren_path = base + "_converted.ren"
        prmag_path = base + "_converted.prmag"

        try:
            nb, unmatched, code_anomalies, skipped_multi_suffix = convert_legacy_auto(
                old_path, complement_path, ren_path)
        except Exception as e:
            messagebox.showerror("Error", f"Conversion failed:\n{e}")
            return

        legacy_results_path = base + ".r"
        try:
            nb_prmag, nb_results = convert_ren_to_r_file(
                ren_path, prmag_path, legacy_results_path=legacy_results_path)
        except Exception as e:
            messagebox.showerror("Error", f"Conversion to .prmag failed:\n{e}")
            return

        msg = f"Converted: {nb} sample(s) -> {ren_path}\n"
        msg += f"Converted: {nb_prmag} sample(s) -> {prmag_path}\n"
        if nb_results:
            msg += f"Converted: {nb_results} result(s) -> {results_path_for(prmag_path)}\n"
        missing_info_path = None
        if unmatched:
            # Liste COMPLETE des echantillons (sample, pas specimen - un
            # "95VN1901A"+"95VN1901B" ne comptent qu'une fois) sans
            # correspondance dans le complement, ecrite dans un fichier
            # texte a part - demande explicite utilisateur ("provide in a
            # text file the list of sample or site... to let the user
            # have a complete list of samples with missing information") :
            # la console tronque a 20 exemples, insuffisant pour aller
            # completer un fichier complement ou une table
            # "Complete sample information" sur un import a des centaines
            # de specimens sans correspondance.
            samples = sorted({s[:-1] if len(s) > 1 else s for s in unmatched})
            missing_info_path = base + "_missing_info.txt"
            with open(missing_info_path, "w", encoding="utf-8") as f:
                f.write("\n".join(samples) + "\n")
            shown = ", ".join(unmatched[:20])
            more = f", ... ({len(unmatched) - 20} more)" if len(unmatched) > 20 else ""
            msg += (f"No complement match ({len(unmatched)} specimen(s), {len(samples)} sample(s)), "
                    f"geology fields set to \"Not Specified\" (to complete later).\n"
                    f"Full sample list written to: {missing_info_path}\n"
                    f"  e.g. {shown}{more}\n")
        if code_anomalies:
            by_pair = {}
            for specimen, etape, cod1, cod2 in code_anomalies:
                by_pair.setdefault((cod1, cod2), []).append(f"{specimen} (step {etape})")
            msg += (f"\nUnrecognized cod1/cod2 ({len(code_anomalies)} measurement(s), "
                    f"{len(by_pair)} distinct code(s)) - please check and replace these in "
                    f"the source file:\n")
            for (cod1, cod2), examples in sorted(by_pair.items(), key=lambda kv: -len(kv[1])):
                shown = ", ".join(examples[:5])
                more = f", ... ({len(examples) - 5} more)" if len(examples) > 5 else ""
                msg += f"  '{cod1}{cod2}' x{len(examples)}: {shown}{more}\n"
        skipped_path = None
        if skipped_multi_suffix:
            # Specimens EXCLUS du fichier de sortie (pas juste "a
            # completer plus tard" comme unmatched/code_anomalies ci-
            # dessus) - demande explicite utilisateur ("do not import
            # samples with multi-char suffices, just provide the list of
            # non imported samples") : liste COMPLETE ecrite a part, meme
            # raison que missing_info_path (console tronquee a 20).
            skipped_path = base + "_not_imported.txt"
            with open(skipped_path, "w", encoding="utf-8") as f:
                f.write("\n".join(sorted(skipped_multi_suffix)) + "\n")
            shown = ", ".join(skipped_multi_suffix[:20])
            more = f", ... ({len(skipped_multi_suffix) - 20} more)" if len(skipped_multi_suffix) > 20 else ""
            msg += (f"\nNOT imported ({len(skipped_multi_suffix)} specimen(s) with a multi-character "
                    f"suffix, e.g. \"10CL0601AP\" - sample name can't be reliably derived):\n"
                    f"Full list written to: {skipped_path}\n"
                    f"  e.g. {shown}{more}\n")
        # Ouvre directement le .prmag converti - demande explicite
        # utilisateur ("ce serait bien d'ouvrir les fichiers convertis a
        # la fin des conversions"). AVANT self._afficher(msg) : le
        # chargement efface la zone de texte (voir _load_data_file), le
        # rapport de conversion (unmatched/code_anomalies/skipped) doit
        # etre affiche APRES, pas remplace par la simple liste de
        # specimens charges.
        self._load_data_file(prmag_path, announce=False)
        self._afficher(msg)
        info = f"{nb} sample(s) converted.\nOutput: {ren_path}\n{prmag_path}"
        if unmatched:
            info += (f"\n{len(unmatched)} specimen(s) had no complement match - geology fields "
                      f"left \"Not Specified\".\nFull sample list: {missing_info_path}")
        if code_anomalies:
            info += (f"\n{len(code_anomalies)} measurement(s) with an unrecognized cod1/cod2 "
                      f"- please check and fix them in the source file (see console).")
        if skipped_multi_suffix:
            info += (f"\n{len(skipped_multi_suffix)} specimen(s) NOT imported (multi-character suffix)."
                      f"\nFull list: {skipped_path}")
        messagebox.showinfo("Conversion complete", info)

    def ouvrir_complete_sample_info_dialog(self):
        """Complete les metadonnees Site/Formation/Age/GC/SMT/Li/Loc/Obs
        (et lat/lon, et en mode specimen sample/site) d'un .prmag DEJA
        converti, depuis une table externe - demande explicite
        utilisateur ("can we put a menu like complete sample information
        with data in a table"), suite du "we can build a routine to
        complete the sample information later on" evoque lors du passage
        a un import legacy "as-is" (voir complete_sample_info.py pour le
        detail exact des deux modes et le format de table attendu).
        Patche le .prmag EN PLACE (une sauvegarde .bak est ecrite avant
        toute modification, une seule fois)."""
        self.text_area.insert(tk.END, "\n--- Complete sample information (Escape to cancel) ---\n", "prompt")
        mode = self._console_input(
            "Table indexed by site (assumes specimen/sample/site already correct) (s) "
            "/ by specimen (also fills sample/site) (p): ", "s")
        if mode is None:
            return
        mode = (mode.strip().lower() or "s")[:1]
        if mode not in ("s", "p"):
            messagebox.showerror("Error", "Must be 's' (site) or 'p' (specimen).")
            return

        prmag_path = filedialog.askopenfilename(
            title="Select the .prmag file to complete",
            filetypes=[("prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not prmag_path:
            return
        table_path = filedialog.askopenfilename(
            title="Select the information table (tab-separated)",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not table_path:
            return

        func = complete_site_info if mode == "s" else complete_specimen_info
        try:
            n_updated, unmatched = func(prmag_path, table_path)
        except Exception as e:
            messagebox.showerror("Error", f"Completion failed:\n{e}")
            return

        msg = f"{os.path.basename(prmag_path)}: {n_updated} specimen(s) updated.\n"
        if unmatched:
            shown = ", ".join(unmatched[:20])
            more = f", ... ({len(unmatched) - 20} more)" if len(unmatched) > 20 else ""
            msg += f"No table match ({len(unmatched)}): {shown}{more}\n"
        self._afficher(msg)
        messagebox.showinfo(
            "Completion done",
            f"{n_updated} specimen(s) updated in {prmag_path}\n"
            f"(backup: {prmag_path}.bak)"
            + (f"\n{len(unmatched)} specimen(s) had no match in the table." if unmatched else "")
            + "\n\nRe-open this file to see the changes (the currently loaded data is not refreshed automatically).",
        )

    def ouvrir_archive_new_data_dialog(self):
        """Archive de nouvelles mesures (acquises APRES la creation du
        .prmag, sur le meme instrument de labo) dans un .prmag DEJA
        existant - demande explicite utilisateur ("After the creation of
        the .prmag file, it is possible that some new data will be
        acquire on the magnetometer in the lab... I need the possibility
        to archive new data acquired in the legacy files. I need also to
        upload those acquire with the JR6 magnetometer") - voir
        import_new_data.py pour le detail exact du format JR6 et de la
        logique de fusion/deduplication (port de
        reference/ImportJR6/ImportJR6data.f95 et importinpmagren.f).

        Un specimen de la nouvelle acquisition SANS correspondance dans
        le .prmag est ECARTE et journalise, PAS cree - demande explicite
        utilisateur ("we assume that data will be archived only for
        specimens already defined in the .prmag file"), a la difference
        du Fortran d'origine. Patche le .prmag EN PLACE (sauvegarde .bak
        avant toute modification, une seule fois) - accepte plusieurs
        fichiers source a la fois (fusionnes avant l'archivage)."""
        self.text_area.insert(tk.END, "\n--- Archive new lab data (Escape to cancel) ---\n", "prompt")
        source = self._console_input(
            "Source format: legacy Rennes file(s) (r) / JR6 file(s) (j): ", "r")
        if source is None:
            return
        source = (source.strip().lower() or "r")[:1]
        if source not in ("r", "j"):
            messagebox.showerror("Error", "Must be 'r' (legacy Rennes) or 'j' (JR6).")
            return

        prmag_path = filedialog.askopenfilename(
            title="Select the .prmag file to archive new data into",
            filetypes=[("prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not prmag_path:
            return

        if source == "r":
            source_paths = filedialog.askopenfilenames(
                title="Select one or more legacy Rennes file(s) with the new data",
                filetypes=[("REN/Text", "*.ren *.txt"), ("All files", "*.*")],
            )
        else:
            source_paths = filedialog.askopenfilenames(
                title="Select one or more JR6 file(s) with the new data",
                filetypes=[("JR6", "*.jr6 *.txt"), ("All files", "*.*")],
            )
        if not source_paths:
            return

        parser = parse_legacy_new_measurements if source == "r" else parse_jr6_file
        new_by_specimen = {}
        try:
            for path in source_paths:
                for specimen, mesures in parser(path).items():
                    new_by_specimen.setdefault(specimen, []).extend(mesures)
        except Exception as e:
            messagebox.showerror("Error", f"Could not read the source file(s):\n{e}")
            return
        if not new_by_specimen:
            messagebox.showwarning("No data", "No usable measurement found in the selected file(s).")
            return

        try:
            n_specimens, n_added, n_dup, unmatched = archive_new_measurements(prmag_path, new_by_specimen)
        except Exception as e:
            messagebox.showerror("Error", f"Archiving failed:\n{e}")
            return

        msg = (
            f"{os.path.basename(prmag_path)}: {n_added} measurement(s) archived "
            f"across {n_specimens} specimen(s) ({n_dup} already present, skipped).\n"
        )
        if unmatched:
            shown = ", ".join(unmatched[:20])
            more = f", ... ({len(unmatched) - 20} more)" if len(unmatched) > 20 else ""
            msg += (f"NOT archived - specimen not found in this .prmag ({len(unmatched)}): "
                    f"{shown}{more}\n")
        self._afficher(msg)
        messagebox.showinfo(
            "Archiving done",
            f"{n_added} measurement(s) archived in {prmag_path}\n"
            f"(backup: {prmag_path}.bak)"
            + (f"\n{len(unmatched)} specimen(s) had no match - not archived." if unmatched else "")
            + "\n\nRe-open this file to see the changes (the currently loaded data is not refreshed automatically).",
        )

    def ouvrir_convert_magic_to_r_dialog(self):
        """Convertit un fichier de contribution MagIC COMBINE (comme pour
        Import MagIC contribution file..., toutes les tables locations/
        sites/samples/specimens/measurements dans un seul .txt) vers le
        nouveau format .r (voir convert_magic_to_r.py) - reconstruit un
        cod1/cod2 a la Rennes a partir de method_codes/treat_*, detecte
        les protocoles hors perimetre Starmac (AMS, hysteresis, MPMS,
        susceptibilite vs champ/frequence/temperature) et les ecarte sans
        faire echouer la conversion, plutot que les demagnetiser/mesurer
        de paleointensite. Le rapport de fin de conversion (protocoles
        ecartes, method_codes non reconnus) est affiche dans la console
        et dans la boite de dialogue.

        Genere AUSSI le .pmagres compagnon (specimens.txt deja interprete
        + sites.txt deja calcule) - demande explicite utilisateur ("is
        the pmagres file also generated during the magic import") - voir
        convert_magic_to_r._convert_magic_results pour le detail."""
        path = filedialog.askopenfilename(
            title="Select the MagIC contribution file (.txt)",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        base, _ext = os.path.splitext(path)
        output_path = base + ".prmag"
        try:
            nb, report, nb_results, nb_means = convert_magic_file(path, output_path)
        except Exception as e:
            messagebox.showerror("Error", f"Conversion failed:\n{e}")
            return
        results_msg = (
            f"Converted: {nb_results} result(s) and {nb_means} site mean(s) -> "
            f"{results_path_for(output_path)}\n"
            if (nb_results or nb_means) else ""
        )
        msg = f"{report}\nConverted: {nb} specimen(s) -> {output_path}\n{results_msg}"
        # Ouvre directement le .prmag converti - demande explicite
        # utilisateur ("ce serait bien d'ouvrir les fichiers convertis a
        # la fin des conversions"). AVANT self._afficher(msg), meme
        # raison que ouvrir_convert_legacy_dialog.
        self._load_data_file(output_path, announce=False)
        self._afficher(msg)
        messagebox.showinfo("Conversion complete", msg)

    def ouvrir_convert_utrecht_to_r_dialog(self):
        """Convertit un ou plusieurs fichiers .col d'Utrecht
        (paleomagnetism.org/PMAG2, JSON malgre l'extension) vers un SEUL
        fichier .prmag combine, PLUS un seul .pmagres compagnon pour les
        interpretations deja calculees par leur logiciel (voir
        convert_utrecht_to_r.py : "import the interpretation in
        Pmagres", "import all the col files in a single .prmag file").
        Site = nom du fichier .col source (sans extension) - ces fichiers
        n'ont pas de notion de site propre, ni extractible du nom de
        specimen ; approximation assumee ("assume that the site is the
        name of the col file"), a raffiner manuellement au besoin.
        Orientation de carotte (azimuth/dip) convertie via la regle
        donnee par l'utilisateur (caz=coreAzimuth+90, cin=90-coreDip) et
        verifiee numeriquement sur les 156 interpretations reelles
        fournies (voir docstring de convert_utrecht_to_r.py) - les vues
        in-situ/apres pendage sont fiables pour ces echantillons."""
        paths = filedialog.askopenfilenames(
            title="Select one or more Utrecht/PMAG2 .col files (JSON)",
            filetypes=[("Utrecht collection", "*.col"), ("All files", "*.*")],
        )
        if not paths:
            return
        default_name = (
            os.path.splitext(os.path.basename(paths[0]))[0] + ".prmag" if len(paths) == 1
            else "utrecht_combined.prmag"
        )
        output_path = filedialog.asksaveasfilename(
            title="Save combined .prmag as",
            defaultextension=".prmag",
            initialfile=default_name,
            initialdir=os.path.dirname(paths[0]),
            filetypes=[("Starmac prmag", "*.prmag"), ("All files", "*.*")],
        )
        if not output_path:
            return
        try:
            nb, nb_results = convert_utrecht_files(list(paths), output_path)
        except Exception as e:
            messagebox.showerror("Error", f"Conversion failed:\n{e}")
            return
        sites = ", ".join(sorted({os.path.splitext(os.path.basename(p))[0] for p in paths}))
        msg = (
            f"Converted: {nb} specimen(s) from {len(paths)} file(s) -> {output_path}\n"
            f"Converted: {nb_results} interpretation(s) -> {results_path_for(output_path)}\n"
            f"\nSite = source .col file name: {sites}\n"
            "Refine manually in Starmac if a file bundles more than one real site.\n"
        )
        # Ouvre directement le .prmag converti - demande explicite
        # utilisateur ("ce serait bien d'ouvrir les fichiers convertis a
        # la fin des conversions"). AVANT self._afficher(msg), meme
        # raison que ouvrir_convert_legacy_dialog.
        self._load_data_file(output_path, announce=False)
        self._afficher(msg)
        messagebox.showinfo("Conversion complete", msg)

    def ouvrir_export_magic_dialog(self):
        """Equivalent GUI de `export2magic` ("export Rennes to Magic",
        fichiers_magic.f) - mode classique uniquement (sites/samples/
        specimens/measurements.txt + locations.txt), sans paleointensite
        ni magnetostratigraphie (voir magic_export.py pour le detail des
        ecarts assumes par rapport au Fortran : tri explicite par site,
        les 2 bugs Fortran reperes sont corriges, formatage numerique
        simple). Les champs Site/Sample/Fm/Age/GC/SMT/Li/Loc utilises
        viennent de la ligne « roche » decodee a l'ouverture du fichier
        (testlect.decode_roche) - un echantillon sans site MagIC decode
        n'aura simplement pas de ligne dans sites.txt."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- export Rennes to Magic (Escape to cancel) ---\n", "prompt")
        lab_analysts = self._console_input(
            "Analysts' names (in quotes, separated by ':'): ", "")
        if lab_analysts is None:
            return

        continent = self._console_input("Continent/ocean (optional): ", "")
        if continent is None:
            return
        country = self._console_input("Country (optional): ", "")
        if country is None:
            return
        region = self._console_input("Region (optional): ", "")
        if region is None:
            return

        out_dir = filedialog.askdirectory(title="MagIC output folder (sites/samples/specimens/measurements.txt)")
        if not out_dir:
            return

        try:
            result = export_to_magic(
                self.selection, self.results, out_dir,
                lab_analysts=lab_analysts,
                continent_ocean=continent, country=country, region=region,
            )
        except OSError as e:
            messagebox.showerror("Error", f"MagIC export failed:\n{e}")
            return

        summary = "\n".join(
            f"  {name}.txt: {result.counts[name]} line(s)" for name in result.paths)
        self._afficher(f"MagIC export finished in {out_dir}:\n{summary}\n")

    def _prompt_magstrat_heights(self):
        """Equivalent du prompt "with magnetostratigraphic data (y/N)"
        commun a exportpmagren et exporttolatex - lit un fichier "id
        depth" (positif vers le haut) si l'utilisateur repond oui.
        Retourne None si annule (Échap), sinon un dict {id: depth} (vide
        si l'utilisateur repond non)."""
        cmagstrat = self._console_input("With magnetostratigraphic data (y/N): ", "N")
        if cmagstrat is None:
            return None
        if cmagstrat.strip().lower() != "y":
            return {}
        path = filedialog.askopenfilename(
            title="Sample/height file (positive upward)",
            filetypes=[("Text", "*.txt *.dat"), ("All files", "*.*")],
        )
        if not path:
            return {}
        heights = {}
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            heights[parts[0].strip()] = float(parts[1])
                        except ValueError:
                            continue
        except OSError as e:
            messagebox.showwarning("Error", f"Could not read {path}:\n{e}")
        return heights

    def ouvrir_export_detailed_dialog(self):
        """Equivalent GUI de `exportpmagren` ("export detailed Rennes",
        dataselect.f) : fichier texte unique avec, pour chaque échantillon
        sélectionné, un bloc de paramètres complet puis le tableau de
        mesures (Dsc/Isc, Dis/Iis, Dtc/Itc = dec/inc en repère échantillon/
        in-situ/après pendage). Séparé de « export Latex » (le Fortran les
        enchaîne automatiquement, ici ce sont deux menus indépendants).
        La déclinaison IGRF n'est pas calculée (IGRF non porté) - affichée
        "n.d" comme le fait déjà le Fortran pour ses propres cas de
        données manquantes."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self.text_area.insert(tk.END, "\n--- export detailed Rennes (Escape to cancel) ---\n", "prompt")
        location = self._console_input("Main study location (country, region): ", "")
        if location is None:
            return
        heights = self._prompt_magstrat_heights()
        if heights is None:
            return
        out_path = filedialog.asksaveasfilename(
            title="Detailed .txt file", defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not out_path:
            return
        try:
            export_detailed_txt(self.selection, location, out_path, heights=heights or None)
        except OSError as e:
            messagebox.showerror("Error", f"Export failed:\n{e}")
            return
        self._afficher(f"Detailed export written: {out_path}\n")

    def ouvrir_export_latex_dialog(self):
        """Equivalent GUI de `exporttolatex` (dataselect.f) : document
        LaTeX (index de sites/échantillons hyperliés + mêmes blocs que
        « export detailed Rennes », plus l'ajustement ChRM le cas échéant)
        - séparé de « export detailed Rennes » (voir docstring de cette
        dernière). L'inclusion des PDF Zijderveld par échantillon
        (`includegraphics{zijder-<id>.pdf}` du Fortran, qui suppose ces
        PDF déjà générés sur disque) n'est pas reproduite."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self.text_area.insert(tk.END, "\n--- export Latex (Escape to cancel) ---\n", "prompt")
        location = self._console_input("Main study location (country, region): ", "")
        if location is None:
            return
        heights = self._prompt_magstrat_heights()
        if heights is None:
            return
        out_path = filedialog.asksaveasfilename(
            title="LaTeX file", defaultextension=".tex",
            filetypes=[("LaTeX", "*.tex"), ("All files", "*.*")],
        )
        if not out_path:
            return
        try:
            export_latex(self.selection, location, out_path, results=self.results, heights=heights or None)
        except OSError as e:
            messagebox.showerror("Error", f"Export failed:\n{e}")
            return
        self._afficher(f"LaTeX export written: {out_path}\n")

    # ------------------------------------------------------------------
    # Graphiques : Zijderveld (equivalent de `plotzijder`/`zijderplot`)
    # ------------------------------------------------------------------

    def load_demo_data(self):
        self.text_area.insert(tk.END, "Ready to import data...\n")

    def _clear_figure(self):
        """Efface juste le CONTENU de la Figure (axes), SANS appeler
        _redraw_canvas - car celle-ci ajuste la geometrie de la fenetre/du
        panneau en fonction de `self.fig.get_size_inches()`, qui a ce stade
        vaut encore la taille de l'ANCIEN graphique affiche, pas celle du
        nouveau qui va etre construit juste apres (`set_size_inches` suit
        toujours cet appel). Appeler `_redraw_canvas()` ici ferait un
        premier ajustement avec la MAUVAISE taille, en course avec le bon
        ajustement suivant - c'etait la cause du glitch d'affichage (taille
        de Figure corrompue, ex. [6.2, 7.7] au lieu de [6.0, 8.5]) qui
        persistait meme apres avoir corrige `_redraw_canvas` elle-meme."""
        self.fig.clear()

    def clear_screen(self):
        self._clear_figure()
        self._redraw_canvas()
        self._current_graphic = None
        self._arai_state = None
        self._paleoint_review_state = None

    def exporter_svg(self):
        """Equivalent fonctionnel de laser1/laser2/... (menu Graphics_SVG) :
        contrairement au Fortran, pas besoin d'un code d'export separe -
        matplotlib ecrit un SVG directement depuis la Figure affichee."""
        if not self.fig.axes:
            messagebox.showwarning("No graphic", "Display a graphic first (e.g. Zijderveld).")
            return

        sample_id = self._current_graphic[1] if self._current_graphic and self._current_graphic[1] else None
        default_name = f"{sample_id}.svg" if sample_id else "graphic.svg"
        path = filedialog.asksaveasfilename(
            title="Export as SVG",
            defaultextension=".svg",
            initialfile=default_name,
            filetypes=[("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            self.fig.savefig(path, format="svg")
        except Exception as e:
            messagebox.showerror("Error", f"SVG export failed:\n{e}")
            return
        messagebox.showinfo("Export successful", f"Graphic exported:\n{path}")

    def afficher_zijderveld(self):
        """Equivalent GUI de la boucle `do i=1,nbech ... call zijder2(...)`
        (linesplans.f) : quand plusieurs echantillons sont selectionnes, on
        les parcourt SEQUENTIELLEMENT, un `read(*,*)` (ici `_console_input`)
        servant de PAUSE entre chaque - pas pour lire une valeur utile, juste
        pour laisser le temps de regarder le diagramme avant de passer au
        suivant (Echap pour arreter la sequence en cours de route)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        samples = [s for s in self.selection if s.mesures]
        if not samples:
            messagebox.showwarning("No measurements", "No selected sample has any measurement.")
            return

        for i, ech in enumerate(samples):
            self._current_graphic = ("zijderveld", ech.id)
            self._refresh_current_graphic()
            if i == len(samples) - 1:
                break
            reste = len(samples) - i - 1
            answer = self._console_input(
                f"[{ech.id}] Enter for next sample ({reste} remaining, Escape to stop): "
            )
            if answer is None:
                break

    def afficher_stereo(self):
        """Equivalent de `stermes` (stereoplot(0))."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self._current_graphic = ("stereo", None)
        self._refresh_current_graphic()

    def afficher_stereo_results(self):
        """Equivalent de `sterres` (stereoplot(1)) : directions moyennes de
        Fisher avec cone de confiance alpha95 (cat1='F'), et grands cercles/
        points pour les ajustements de droite/plan (cat1 'L'/'f'/'P') de
        self.results (a charger au prealable via Ajustement... ou
        Select results...)."""
        if not self.results:
            messagebox.showwarning(
                "No results",
                "No results in memory - run a fit or "
                "load some via Pmag data > Select results...")
            return
        self._current_graphic = ("stereo_results", None)
        self._refresh_current_graphic()

    def afficher_xygraph(self):
        """Equivalent de `xygraph` (courbe de désaimantation)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self._current_graphic = ("xygraph", None)
        self._refresh_current_graphic()

    def afficher_susceptibilite(self):
        """Equivalent de `suscep`."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self._current_graphic = ("susceptibility", None)
        self._refresh_current_graphic()

    def afficher_irm(self):
        """Acquisition d'aimantation remanente isotherme (IRM, mesures
        cod1='I') - pas dans le Fortran, port des scripts GMT de
        l'utilisateur (Scripts_IRM_GMT), demande explicite."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        if not has_irm_data(self.selection):
            messagebox.showwarning(
                "No IRM data",
                "No selected sample has an IRM measurement (step code 'I').",
            )
            return
        self._current_graphic = ("irm", None)
        self._refresh_current_graphic()

    def afficher_arai(self):
        """Equivalent GUI de `paleoin` (diagramme d'Arai/Thellier) - reprend
        la sequence interactive du Fortran (plotpaleoint2.f) presque telle
        quelle : champ com:, tableau des pas NRM/TRM, verifications pTRM,
        courbure de Paterson (Taubin+LMA, adjustcircle.f95, deux appels
        comme l'original), directions ancree/libre + DANG, tableau CRM
        (%rcrm), statistiques de Coe (b/sigma/ccr/f/g/q/h), correction
        d'anisotropie (tenseur .ANI 'A0') et gamma. Verifie octet-pres
        contre un transcript Fortran reel (echantillon 06A,
        SanJuan_Pmag.ren, ni,nj=4,10) - y compris la coquille reelle du
        source (constante '3.14152927' au lieu de pi dans `compute_crm`),
        preservee pour fidelite.

        Volontairement HORS PERIMETRE (non portes) : ARN a 2/3 composantes
        (icomposante 2/3 - traites comme 1), rf1/rf2 et le bloc Prevot Ql
        (code mort cote Fortran - la lecture de susceptibilite y est
        commentee, `sus`/`ql` valent toujours 0)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- Arai diagram / paleointensity (Escape to cancel) ---\n", "prompt")

        ani_path = None
        if self.results_path:
            ani_path = os.path.splitext(self.results_path)[0] + ".ANI"
            if not os.path.exists(ani_path):
                ani_path = os.path.splitext(results_path_for(self.results_path))[0] + ".ANI"
            if not os.path.exists(ani_path):
                ani_path = None

        traite = 0
        for ech in self.selection:
            if len(ech.mesures) < 2:
                continue
            # orientation=1 (echantillon) FORCE, pas self.orientation.get()
            # - Fortran (plotpaleoint2.f/visi_Paleoint.f) : "iorient=1"
            # sans condition en tete de routine, ignore le reglage
            # d'orientation de l'appli - demande explicite utilisateur
            # ("the combined plots should be in sample coordinates").
            points, checks, arno = compute_arai(ech, orientation=1)
            if len(points) < 2:
                continue  # pas de codes N/R/V/P exploitables pour cet echantillon
            traite += 1

            # champ/methode auto-detectes depuis les mesures elles-memes
            # (treat_dc_field pour .prmag, cod1 S/R/V) plutot que du seul
            # commentaire - demande explicite utilisateur ("extract the
            # field value from treat_dc_field... the method is IZZI if
            # there are S codes and Thellier if there are R and V
            # codes"). Repli sur l'ancien parsing du commentaire si aucune
            # mesure ne porte de treat_dc_field (.ren).
            method_auto, ichamp_auto = detect_method_and_hlab(ech.mesures, ech.com)
            confirm = self._console_input(
                f" method={method_auto}, lab field detected = {ichamp_auto:g} uT "
                f"(comment=[{ech.com}]) - is this OK Y/n : ", "Y")
            if confirm is None:
                return
            ncom = ichamp_auto != 0 and confirm.strip().lower() != "n"
            ichamp = ichamp_auto if ncom else 0

            icomp_s = self._console_input(
                " ARN one (1) two (2) or several components (3) : ", "1")
            if icomp_s is None:
                return
            try:
                icomposante = int(icomp_s.strip() or "1")
            except ValueError:
                icomposante = 1
            if icomposante not in (2, 3):
                icomposante = 1
            if icomposante != 1:
                self._afficher(
                    "(ARN 2/3 components not yet ported - treated as 1 component.)\n")

            plot_inv = self._console_input(" Plot TRM in - and + : y/N", "N")
            if plot_inv is None:
                return

            if ichamp != 0:
                hlab = float(ichamp)
                parsed = parse_com_field(ech.com)  # imois/iannee : purement informatif
                self._afficher(
                    f" ech:{ech.id:<12}  field of {ichamp}t "
                    f"month:{parsed['imois']}year:{parsed['iannee']}\n"
                )
            else:
                hlab_s = self._console_input(
                    "  intensity of the laboratory field (in microteslas): ", "0")
                if hlab_s is None:
                    return
                try:
                    hlab = float(hlab_s or 0.0)
                except ValueError:
                    messagebox.showerror("Error", "Hlab must be a number.")
                    continue

            volumasse = 1.0e-3 if ech.norme == "m" else 1.0

            out = io.StringIO()
            out.write(" paleointensity diagram NRM-TRM\n")
            out.write(f"  numero: {ech.id}\n\n")
            out.write("    n  temp.    arn     dec    inc      atr      mom.atr    dec    inc\n\n")
            for i, p in enumerate(points, start=1):
                mom = p.xp * arno * volumasse
                out.write(
                    f"  {i:3d}  {p.temp:4.0f}   {p.yp:9.3E}  {p.decl:7.2f} {p.aincl:7.2f}   "
                    f"{p.xp:9.3E}   {mom:9.3E}   {p.dec:7.1f} {p.winc:7.1f}\n"
                )
            if checks:
                pt_by_k = {p.k: p for p in points}
                out.write("\n   results of the ptrm checks\n")
                out.write("      temp      mom        %atr   ecart/atrt  ecart/atrp\n")
                for c in checks:
                    target = pt_by_k[c.k].temp if c.k in pt_by_k else c.temp
                    mom = c.xt * arno
                    ecart = c.xt - c.xtptrm
                    ratio = ecart / c.xtptrm if c.xtptrm else 0.0
                    out.write(
                        f"   {target:4.0f} {c.temp:4.0f}   {mom:9.3E}   {c.xt:6.2f}   "
                        f"{ecart:6.2f}      {ratio:6.2f}\n"
                    )
            self._afficher(out.getvalue())

            self._arai_state = (ech, points, checks, arno, None, hlab)
            self._current_graphic = ("arai", ech.id)
            self._refresh_current_graphic()

            while True:
                range_s = self._console_input(
                    f" calculation of the slope between steps ni,nj (1-{len(points)}, 0 = skip): ",
                    f"1 {len(points)}")
                if range_s is None:
                    return
                parts = range_s.split()
                if len(parts) == 1 and parts[0].strip() == "0":
                    break
                try:
                    if len(parts) >= 2:
                        n1, n2 = int(parts[0]), int(parts[1])
                    else:
                        n1, n2 = 1, len(points)
                except ValueError:
                    messagebox.showerror("Error", "Enter two integers (ni nj).")
                    continue
                if not (1 <= n1 <= n2 <= len(points)) or n2 - n1 < 1:
                    messagebox.showerror("Error", "Invalid point range (at least 2 points).")
                    continue

                fit = fit_arai_line(points, n1, n2, hlab=hlab)
                direction = fit_arai_direction(points, n1, n2, ech, orientation=1)  # SC forcee, voir compute_arai ci-dessus
                gamma = 90.0 - points[n2 - 1].winc

                fcor = hcorani = None
                if ani_path and direction.anchored_specimen_frame is not None:
                    tensor = read_ani_tensor(ani_path, ech.id, "A0")
                    if tensor is not None:
                        fcor = compute_anicor_factor(tensor, direction.anchored_specimen_frame)
                        hcorani = fit.h * fcor

                # courbure de Paterson (adjustcircle.f95, AraiCurvature) :
                # appel 1 sur (0,1)+points[1..n2] (utilise dans la ligne de
                # resultats finale), appel 2 sur points[n1..n2] uniquement
                # (affiche mais pas repris dans la ligne finale - transcrit
                # tel quel, meme choix que le Fortran).
                xpat1 = [0.0] + [p.xp for p in points[:n2]]
                ypat1 = [1.0] + [p.yp for p in points[:n2]]
                curv0 = arai_curvature(xpat1, ypat1)
                xpat2 = [p.xp for p in points[n1 - 1:n2]]
                ypat2 = [p.yp for p in points[n1 - 1:n2]]
                curv1 = arai_curvature(xpat2, ypat2)

                crm = compute_crm(points, n1, n2, arno) if n2 - n1 >= 1 else None

                self._arai_state = (ech, points, checks, arno, fit, hlab)
                self._refresh_current_graphic()

                res = io.StringIO()
                res.write(f"\n calculation of the slope between steps ni,nj :{n1} {n2}\n\n")
                res.write(" -----------------------------------\n")
                res.write(
                    f" circle Parameter for initialization after Taubin: {curv0.taubin_a:9.6f} "
                    f"{curv0.taubin_b:9.6f} {curv0.taubin_r:8.5f}\n"
                )
                res.write(
                    f" circle Parameter (a,b,r) after LMA: {curv0.lma_a:9.6f} "
                    f"{curv0.lma_b:9.6f} {curv0.lma_r:8.5f}\n"
                )
                res.write(" curvature calculated from point (0,1) to n2\n")
                res.write(f" parameter k (1/r) : {curv0.k:7.4f}  error SSE :{curv0.sse:7.5f}\n")
                res.write("\n -----------------------------------\n")
                res.write(
                    f" circle Parameter for initialization after Taubin: {curv1.taubin_a:9.6f} "
                    f"{curv1.taubin_b:9.6f} {curv1.taubin_r:8.5f}\n"
                )
                res.write(
                    f" circle Parameter (a,b,r) after LMA: {curv1.lma_a:9.6f} "
                    f"{curv1.lma_b:9.6f} {curv1.lma_r:8.5f}\n"
                )
                res.write(" curvature calculated from point n1 to n2\n")
                res.write(f" parameter k (1/r) : {curv1.k:7.4f}  error SSE :{curv1.sse:7.5f}\n")
                res.write("\n -----------------------------------\n\n")
                if direction.anchored_dec is not None:
                    res.write(
                        f" anchored direction: dec={direction.anchored_dec:6.1f}  "
                        f"inc={direction.anchored_inc:6.1f}  mad={direction.anchored_mad:5.1f}  "
                        f"nb points: {direction.nb}\n"
                    )
                if direction.free_dec is not None:
                    res.write(
                        f"free direction:  dec={direction.free_dec:6.1f}  "
                        f"inc={direction.free_inc:6.1f}  mad={direction.free_mad:5.1f}  "
                        f"nb points: {direction.nb}\n"
                    )
                if direction.dang is not None:
                    res.write(f"\n Lisa Tauxe DANG {direction.dang:6.1f}\n")
                res.write("(rf1/rf2 not yet ported.)\n")
                if crm is not None:
                    res.write("\n temp   mom.crm       temp   mom.2,6       temp   mom.crm\n")
                    ks = list(range(n1, n2 + 1))
                    k = n1
                    while k <= n2:
                        kkj = min(k + 2, n2)
                        cols = []
                        for j in range(k, kkj + 1):
                            cols.append(f"{points[j-1].temp:4.0f}   {crm.values[j]:10.4E}")
                        res.write("  " + "     ".join(cols) + "\n")
                        k += 2
                res.write(f"\n slope b= {fit.b:8.4f}\n")
                if fit.sigma:
                    res.write(f"\n sigma = {fit.sigma:7.4f}\n")
                    res.write(f"\n linear correlation coefficient: {fit.ccr:9.5f}\n")
                if crm is not None:
                    res.write(f"\n crm max ={crm.crmmax:8.4f}    deltatrm={crm.dtrm:8.4f}\n")
                res.write(
                    f"\n{ech.id:<12}   f={fit.f:7.4f}    g={fit.g:7.4f}    q={fit.qq:7.4f}"
                    f"    ccr={fit.ccr:8.5f}  h={fit.h:8.3f}"
                    + (f"     % rcrm={crm.rcrm:7.3f}\n" if crm is not None else "\n")
                )
                res.write(
                    "\nNum          t1   t2   N    f       g       q     mad   dang   Hlab"
                    "     b       sb     sb/b     ccr      H      fcor  Hcorani  gamma     k     k_sse\n"
                )
                res.write(
                    f"{ech.id:<12} {points[n1-1].temp:4.0f} {points[n2-1].temp:4.0f}  "
                    f"{n2 - n1 + 1:2d}  {fit.f:6.3f}  {fit.g:6.3f}  {fit.qq:6.3f}  "
                    f"{(direction.free_mad or 0.0):5.1f}  {(direction.dang or 0.0):5.1f}  "
                    f"{hlab:5.1f}  {fit.b:7.4f}  {fit.sigma:6.4f}  "
                    f"{(fit.sigma / fit.b if fit.b else 0.0):7.4f}  {fit.ccr:8.5f}  "
                    f"{fit.h:6.2f}  "
                    + (f"{fcor:6.3f}  {hcorani:6.2f}" if fcor is not None else "   -       -  ")
                    + f"  {gamma:5.1f}  {curv0.k:7.4f}  {curv0.sse:7.5f}\n"
                )
                self._afficher(res.getvalue())

                decision = self._console_input(
                    "Redo (r) / next sample (Enter): ", "")
                if decision is None:
                    return
                if decision.strip().lower() == "r":
                    continue
                break

        if traite == 0:
            messagebox.showinfo(
                "No diagram",
                "No selected sample has usable paleointensity measurements (N/R/V/P codes).",
            )

    def ouvrir_openfilepint_dialog(self):
        """Equivalent GUI de `openfilepint`/`visi_paleoin` ("View Paleoint
        Results", visi_Paleoint.f) : PAS un traitement en lot - l'objectif
        (comme le Fortran) est de revisiter RAPIDEMENT un traitement deja
        effectue, en laissant l'utilisateur choisir l'echantillon dans une
        liste numerotee (fichier "echantillon Tmin Tmax [taux refroid.]",
        une ligne par determination deja faite), rejouee en boucle jusqu'a
        'q' - exactement le motif `openfilepint` (menu numerote, prompt
        "Type the line number...(q to quit)") plutot qu'un traitement de
        toutes les lignes d'un coup.

        Le graphique associe est la mise en page combinee de
        `visi_paleoin` (`boite(1)/(2)/(3)`) : diagramme d'Arai en haut,
        Zijderveld en bas a gauche, stereo NRM/TRM en bas a droite
        (`build_paleoint_review_figure`) - PAS le panneau Arai seul.

        Comme le Fortran (qui cherche l'echantillon dans `pmag(:)` - la
        TOTALITE des donnees chargees, pas une selection prealable),
        chaque echantillon est retrouve dans self.donnees (meme raison que
        afficher_visres, cf. son docstring). Le champ `com:` fournit Hlab
        silencieusement (pas de prompt de confirmation, ni de choix ARN
        composante/signe TRM - tous fixes cote Fortran dans cette
        routine, contrairement a `afficher_arai`/`paleoin`) ; le taux de
        refroidissement optionnel (4e colonne) n'est qu'un multiplicateur
        final sur H (`rHfinal=H*corcool`), PAS la correction complete de
        `vitref`/Cooling rate."""
        if not self.donnees:
            messagebox.showwarning("No data", "Load a .ren file first.")
            return
        list_path = filedialog.askopenfilename(
            title="List file (sample Tmin Tmax [cooling rate])",
            filetypes=[("Text", "*.txt *.lst *.dat"), ("All files", "*.*")],
        )
        if not list_path:
            return
        try:
            with open(list_path, "r", encoding="iso-8859-1", errors="replace") as f:
                lines = [l for l in f.read().splitlines() if l.strip() and not l.strip().startswith("!")]
        except OSError as e:
            messagebox.showerror("Error", f"Could not read {list_path}:\n{e}")
            return

        entries = []  # (sample_id, tmin, tmax, cooling)
        for line in lines:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                sample_id = parts[0]
                tmin, tmax = int(float(parts[1])), int(float(parts[2]))
                cooling = float(parts[3]) if len(parts) > 3 else 0.0
            except ValueError:
                continue
            entries.append((sample_id, tmin, tmax, cooling))
        if not entries:
            messagebox.showwarning("Empty list", f"No usable entry found in {list_path}.")
            return

        ani_path = None
        if self.results_path:
            ani_path = os.path.splitext(self.results_path)[0] + ".ANI"
            if not os.path.exists(ani_path):
                ani_path = os.path.splitext(results_path_for(self.results_path))[0] + ".ANI"
            if not os.path.exists(ani_path):
                ani_path = None

        self.text_area.insert(
            tk.END, "\n--- Rapid view of previous Paleointensity determinations (Escape to quit) ---\n", "prompt")

        last_iligne = None  # dernier numero de ligne affiche - permet d'avancer par simple Retour
        while True:
            menu = io.StringIO()
            for i in range(0, len(entries), 3):
                row = entries[i:i + 3]
                menu.write("   ".join(
                    f"{i + j + 1:3d}  {e[0]:<12}" for j, e in enumerate(row)) + "\n")
            self._afficher(menu.getvalue())

            choice = self._console_input(
                "Type the line number to select the sample (Enter = next, q to quit) : ", "")
            if choice is None:
                return
            choice = choice.strip()
            if choice.lower().startswith("q"):
                break
            if not choice:
                # Retour simple : avance sur la ligne suivante (demande
                # explicite utilisateur - "a simple return just
                # incremented the selection... You could also start to 11
                # and the next one was 12") - seulement s'il y a deja eu
                # un choix explicite ; sans historique, un Retour a vide
                # quitte comme avant (rien a partir de quoi avancer).
                if last_iligne is None:
                    break
                iligne = last_iligne + 1
                if iligne > len(entries):
                    self._afficher("(no more entries)\n")
                    break
            else:
                try:
                    iligne = int(choice)
                except ValueError:
                    continue
                if not (1 <= iligne <= len(entries)):
                    continue
            last_iligne = iligne

            sample_id, tmin, tmax, cooling = entries[iligne - 1]
            matches = select_samples(
                self.donnees, sample_id, step_min=0, step_max=2000,
                demag1="*", demag2="*", verbose=False)
            ech = matches[0] if matches else None
            if ech is None:
                self._afficher(f"{sample_id} sample not found\n")
                continue
            # rend l'echantillon revisite ACTIF (self.selection), comme
            # toutes les autres routines qui affichent un echantillon
            # (afficher_zijderveld/afficher_arai...) - demande explicite
            # utilisateur ("does not select the sample") : cette routine
            # ne le faisait pas, seul self._current_graphic/
            # _paleoint_review_state (pilotant uniquement le graphique)
            # etaient mis a jour, donc les autres menus (Zijderveld,
            # Delete some data, export SVG...) ne voyaient pas
            # l'echantillon comme selectionne apres l'avoir revisite ici.
            self.selection = [ech]

            # orientation=1 (echantillon) FORCE, pas self.orientation.get()
            # - Fortran (plotpaleoint2.f/visi_Paleoint.f) : "iorient=1"
            # sans condition en tete de routine, ignore le reglage
            # d'orientation de l'appli - demande explicite utilisateur
            # ("the combined plots should be in sample coordinates").
            points, checks, arno = compute_arai(ech, orientation=1)
            if len(points) < 2:
                self._afficher(f"{sample_id}: no usable paleointensity measurements (N/R/V/P).\n")
                continue

            n1 = next((i + 1 for i, p in enumerate(points) if p.temp == tmin), None)
            n2 = next((i + 1 for i, p in enumerate(points) if p.temp == tmax), None)
            if n1 is None or n2 is None or n2 - n1 < 1:
                self._afficher(
                    f"{tmin} {tmax} please check the temperature interval for: {ech.id}\n")
                continue

            _method_auto, ichamp_auto = detect_method_and_hlab(ech.mesures, ech.com)
            if ichamp_auto:
                hlab = float(ichamp_auto)
            else:
                hlab_s = self._console_input(
                    "  intensity of the laboratory field (in microteslas): ", "0")
                if hlab_s is None:
                    return
                try:
                    hlab = float(hlab_s or 0.0)
                except ValueError:
                    hlab = 0.0

            fit = fit_arai_line(points, n1, n2, hlab=hlab)
            direction = fit_arai_direction(points, n1, n2, ech, orientation=1)  # SC forcee, voir compute_arai ci-dessus
            gamma = 90.0 - points[n2 - 1].winc

            fcor = hcorani = None
            if ani_path and direction.anchored_specimen_frame is not None:
                tensor = read_ani_tensor(ani_path, ech.id, "A0")
                if tensor is not None:
                    fcor = compute_anicor_factor(tensor, direction.anchored_specimen_frame)
                    hcorani = fit.h * fcor

            h_final = hcorani if hcorani is not None else fit.h
            if cooling:
                h_final = h_final * cooling

            crm = compute_crm(points, n1, n2, arno)
            xpat1 = [0.0] + [p.xp for p in points[:n2]]
            ypat1 = [1.0] + [p.yp for p in points[:n2]]
            curv0 = arai_curvature(xpat1, ypat1)

            self._paleoint_review_state = (ech, points, checks, arno, fit)
            self._current_graphic = ("paleoint_review", ech.id)
            self._refresh_current_graphic()

            res = io.StringIO()
            res.write(f"\n{ech.id} : Tmin={tmin} Tmax={tmax}\n")
            if direction.anchored_dec is not None:
                res.write(
                    f" anchored direction: dec={direction.anchored_dec:6.1f}  "
                    f"inc={direction.anchored_inc:6.1f}  mad={direction.anchored_mad:5.1f}\n"
                )
            res.write(
                f" f={fit.f:.3f}  g={fit.g:.3f}  q={fit.qq:.3f}  ccr={fit.ccr:.3f}  "
                f"h={fit.h:.2f}  % rcrm={crm.rcrm:.1f}\n"
            )
            res.write(
                f" Hlab={hlab:.1f}µT"
                + (f"  fcor={fcor:.3f}  Hcorani={hcorani:.1f}µT" if fcor is not None else "")
                + (f"  fcorCool={cooling:.3f}  HcorCool={h_final:.1f}µT" if cooling else "")
                + f"  gamma={gamma:.1f}  k={curv0.k:.4f}\n"
            )
            self._afficher(res.getvalue())

            # Second traitement, PARALLELE et INDEPENDANT du natif Starmac
            # ci-dessus : appelle le code PmagPy/MagIC reel (pmag.PintPars)
            # sur le MEME intervalle Tmin/Tmax - demande explicite
            # utilisateur ("a second parallel processing of paleointensity.
            # the one from Magic"). Echec attendu et non bloquant pour une
            # partie reelle des specimens (protocoles hors IZZI/Thellier
            # standard, voir paleointensity_magic.py) - la boucle de revue
            # continue, avec juste une ligne d'explication.
            try:
                magic_result = compute_magic_paleointensity(ech, step_first=tmin, step_last=tmax)
                self._afficher(format_magic_paleointensity(magic_result))
            except Exception as e:
                self._afficher(f"MagIC/PmagPy paleointensity: not computed for {ech.id} ({e})\n")

    def afficher_visres(self):
        """Equivalent GUI de `visres` ("data+interpretation",
        plotorthog.f:8-149) : boucle sur CHAQUE résultat de self.results
        INDIVIDUELLEMENT (pas regroupé par échantillon) - pour cat1=='L'
        affiche le Zijderveld de l'échantillon avec CE seul ajustement (les
        autres ajustements du même échantillon ne sont PAS superposés,
        comme le Fortran qui réduit `tr` à un seul élément avant d'appeler
        `zijder`) ; pour cat1 in ('P','f') affiche Stereo Results pour ce
        seul résultat. Les résultats 'F' (moyennes) sont ignorés, comme
        dans le Fortran (aucune des 3 branches if ne les traite). Même
        motif pause/décision (Entrée = suivant, Échap = arrêter) que les
        autres boucles (ajuslig, Arai...).

        Comme le Fortran (qui reconstruit l'échantillon DIRECTEMENT depuis
        `pmag(:)`, la totalité des données chargées, avec etapmin=0/
        etapmax=2000/demag1='*'/demag2='*' figés dans le code - PAS depuis
        une sélection préalablement filtrée), l'échantillon de chaque
        résultat est retrouvé dans self.donnees et non dans self.selection :
        sinon, un résultat chargé via « Select results... » (selres) sans
        avoir aussi fait « Select samples... » sur le même échantillon ne
        trouvait jamais ses données et data+interpretation semblait ne
        rien voir malgré des résultats bien sélectionnés.

        Affiche aussi l'évaluation de qualité de CE résultat (voir
        interpretation_quality.evaluate_result) - demande explicite
        utilisateur ("dans l'option data+interpretation, est-ce possible
        de lister l'évaluation de l'échantillon?")."""
        if not self.results:
            messagebox.showwarning(
                "No results", "No results saved (see Calcul > Ajustement...).")
            return
        if not self.donnees:
            messagebox.showwarning(
                "No data", "Load a .ren file first.")
            return

        self.text_area.insert(tk.END, "\n--- Data + interpretation (Escape to stop) ---\n", "prompt")
        traite = 0
        for r in self.results:
            if r.cat1 not in ("L", "P", "f"):
                continue
            matches = select_samples(
                self.donnees, r.id, step_min=0, step_max=2000, demag1="*", demag2="*", verbose=False)
            ech = matches[0] if matches else None
            if ech is None:
                continue
            traite += 1
            self._clear_figure()

            if r.cat1 == "L":
                self.fig.set_size_inches(5.5, 8.5, forward=True)
                build_zijderveld_figure(ech, orientation=self.orientation.get(), fits=[r], fig=self.fig)
            else:
                self.fig.set_size_inches(5.5, 5.5, forward=True)
                build_stereo_results_figure(
                    [r], orientation=self.orientation.get(), nbech=1, fig=self.fig)
            self._fit_figure_to_data()
            self._redraw_canvas()

            dec, inc = _correct_dec_inc(r, self.orientation.get())
            self._afficher(f"{r.id} ({r.cat1}): dec={dec:.1f}  inc={inc:.1f}\n")
            # Evaluation de l'interpretation affichee - demande explicite
            # utilisateur ("dans l'option data+interpretation, est-ce
            # possible de lister l'evaluation de l'echantillon?"). None
            # pour un resultat 'f' (moyenne Fisher, pas un fit ligne/plan -
            # voir interpretation_quality.evaluate_result), silencieusement
            # pas de ligne d'evaluation dans ce cas, comme pour "mean:".
            report = evaluate_result(r, self.donnees)
            if report is not None:
                self._afficher("\n" + format_quality_report([report]) + "\n")
            decision = self._console_input("Next (Enter) / stop (Escape): ", "")
            if decision is None:
                return

        if traite == 0:
            messagebox.showinfo(
                "No displayable result",
                "No line/plane type result (« mean: » averages are not "
                "displayed by « data+interpretation », matching the Fortran).",
            )

    def ouvrir_convertthelli_dialog(self):
        """Equivalent GUI de `convertthelli` ("Thellier >> NRM") : convertit
        la séquence Thellier de chaque échantillon sélectionné en séquence
        NRM/demag simple (cod1='D'), utilisable comme un Zijderveld classique."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        if not messagebox.askyesno(
            "Thellier >> NRM",
            f"Convert the Thellier sequence of {len(self.selection)} sample(s) "
            "into a simple NRM/demag sequence - irreversible for this session. Continue?",
        ):
            return
        for ech in self.selection:
            before = len(ech.mesures)
            convert_thellier_to_nrm(ech)
            self._afficher(f"{ech.id}: {before} -> {len(ech.mesures)} measurement(s).\n")

    def ouvrir_removestep_dialog(self):
        """Equivalent GUI de `removestep` ("Remove step") : supprime toutes
        les lignes d'un palier de démagnétisation pour un échantillon."""
        if len(self.selection) != 1:
            messagebox.showwarning("Invalid selection", "Select a single sample.")
            return
        ech = self.selection[0]
        steps_text = "  ".join(
            f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
        self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
        etape_s = self._console_input("Step (temperature/field) to remove: ", "")
        if etape_s is None:
            return
        try:
            etape = int(etape_s)
        except ValueError:
            messagebox.showerror("Error", "Must be an integer.")
            return
        if not messagebox.askyesno(
            "Remove step",
            f"Remove all measurements of {ech.id} at step {etape} - irreversible "
            "for this session. Continue?",
        ):
            return
        removed = remove_step(ech, etape)
        self._afficher(f"{ech.id}: {removed} measurement(s) removed ({len(ech.mesures)} remaining).\n")

    def ouvrir_elimine_grm_dialog(self):
        """Equivalent GUI de `elimineGRM` ("Suppress GRM") : réduit les
        triplets de mesures GRM (cod1='F', cod2 'X','Y','Z' consécutifs)
        en un seul point corrigé, méthode au choix (1 = substitution
        axe-par-axe, x du point X / y du point Y / z du point Z ; 2 =
        moyenne simple des 3 points)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        method_s = self._console_input("Method 1 (axis-by-axis substitution) or 2 (average): ", "1")
        if method_s is None:
            return
        try:
            method = int(method_s)
        except ValueError:
            method = 1
        if method not in (1, 2):
            method = 1
        if not messagebox.askyesno(
            "Suppress GRM",
            f"Reduce the GRM triplets (X/Y/Z) of {len(self.selection)} sample(s) "
            "with method " + ("axis-by-axis substitution" if method == 1 else "average")
            + " - irreversible for this session. Continue?",
        ):
            return
        for ech in self.selection:
            n = eliminate_grm(ech, method=method)
            self._afficher(f"{ech.id}: {n} triplet(s) reduced ({len(ech.mesures)} measurements remaining).\n")

    def ouvrir_convzmoins_dialog(self):
        """Equivalent GUI de `convzmoins` ("ConvertZ- 2G") : inverse y/z et
        recode en 'R' les mesures cod1='Z' cod2='-'."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        if not messagebox.askyesno(
            "ConvertZ- 2G",
            f"Convert Z- measurements to R for {len(self.selection)} sample(s) "
            "- irreversible for this session. Continue?",
        ):
            return
        total = 0
        for ech in self.selection:
            n = convert_z_minus(ech)
            total += n
            if n:
                self._afficher(f"{ech.id}: {n} Z- measurement(s) converted.\n")
        if not total:
            self._afficher("No Z- measurement found in the selection.\n")

    def ouvrir_exportthellier_dialog(self):
        """Equivalent GUI de `exportthellier` ("export to ThellierTool") :
        écrit un fichier .tdt par échantillon sélectionné, dans un dossier
        choisi par l'utilisateur."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        out_dir = os.path.dirname(self.results_path) if self.results_path else os.getcwd()
        chosen = filedialog.askdirectory(initialdir=out_dir, title="ThellierTool export folder (.tdt)")
        if not chosen:
            return
        written = []
        for ech in self.selection:
            try:
                path = export_thellier_tdt(ech, chosen)
                written.append(path)
            except OSError as e:
                messagebox.showerror("Error", f"Export failed for {ech.id}:\n{e}")
                return
        self._afficher("Exported .tdt files:\n" + "\n".join(written) + "\n")

    # ------------------------------------------------------------------
    # Aides communes
    # ------------------------------------------------------------------

    def _not_implemented(self, name):
        """Entree de menu presente dans le Fortran d'origine (StarmacOSX_x.f95,
        menu Calcul) mais pas encore portee en Python."""
        messagebox.showinfo("Not yet implemented", f"« {name} » has not been ported yet.")

    def _archive_only(self, fit):
        """Equivalent de `call archivres` seul (sans toucher self.results) -
        pour les flux qui ajoutent `fit` a self.results en amont a titre
        d'apercu (ex. superposition sur le Zijderveld avant confirmation),
        et n'archivent qu'une fois la sauvegarde confirmee par l'utilisateur."""
        if not self.results_path:
            return
        try:
            _, self._archived_ids = archivres(fit, self.results_path, self._archived_ids)
        except OSError as e:
            messagebox.showwarning(
                "Archiving failed",
                f"The result is in memory but could not be written to "
                f"{self.results_path}:\n{e}",
            )

    def _save_result(self, fit):
        """Equivalent de `tr(nbres+1)=res; nbres=nbres+1; call archivres` :
        ajoute `fit` a self.results ET l'archive immediatement dans le
        fichier .r (equivalent filr) - le Fortran appelle archivres des
        qu'un resultat est sauvegarde, pas seulement en fin de session."""
        self.results.append(fit)
        self._archive_only(fit)

    def _save_results(self, fits):
        for fit in fits:
            self._save_result(fit)

    def _afficher(self, text):
        """Ajoute `text` a la suite du contenu existant (n'efface plus la
        console a chaque appel de menu), et fait defiler jusqu'au nouveau
        contenu."""
        if self.text_area.get("1.0", "end-1c").strip():
            self.text_area.insert(tk.END, "\n" + "-" * 60 + "\n")
        self.text_area.insert(tk.END, text)
        self.text_area.see(tk.END)

    @staticmethod
    def _parse_demag(code):
        code = (code or "").strip()
        if not code or code == "*":
            return "*", "*"
        if len(code) == 1:
            return code[0], "*"
        return code[0], code[1]

    def _console_input(self, prompt, default=""):
        """Lit une ligne tapee au clavier DANS la fenetre texte, comme une
        vraie console (equivalent du `write(*,...)` + `read(*,*)` du Fortran
        d'origine) : affiche `prompt`, laisse l'utilisateur taper, puis rend
        la main des qu'il appuie sur Entree. Renvoie None si annule (Echap).

        `default` n'est PAS pre-insere dans le texte (essaye au debut, mais
        source de confusion : si l'utilisateur tapait sans d'abord effacer
        la valeur pre-remplie, elle s'affichait concatenee avec sa saisie,
        ex. "dm" pour un defaut "d" suivi d'un "m" tape) - applique
        uniquement si la ligne validee est vide (Entree seule).

        Bloque l'appelant sans geler l'interface : `wait_variable` continue
        de faire tourner la boucle d'evenements Tk pendant l'attente (pas de
        thread ni de callback a gerer cote appelant).

        Le prompt commence TOUJOURS sur une ligne neuve - demande explicite
        utilisateur ("there are still long lines text output and expecting
        an answer at the end the line... display the question (in red) in
        a new line") : sans ce saut de ligne, un prompt arrivant juste
        apres un texte deja affiche SANS retour a la ligne final (ex. le
        dernier chiffre d'un rapport) se retrouvait concatene a sa suite
        sur la meme ligne, illisible. La couleur rouge existe deja (tag
        "prompt", voir tag_configure) - seul le retour a la ligne manquait
        ici, point d'entree UNIQUE de tous les prompts interactifs de
        l'appli (corrige une fois pour toutes plutot qu'a chaque site
        d'appel individuellement)."""
        last_char = self.text_area.get("end-2c", "end-1c")
        if last_char not in ("", "\n"):
            self.text_area.insert(tk.END, "\n")
        self.text_area.insert(tk.END, prompt, "prompt")
        start_index = self.text_area.index("end-1c")
        self.text_area.mark_set(tk.INSERT, tk.END)
        self.text_area.see(tk.END)
        self.text_area.focus_set()

        outcome = {"value": None}
        done = tk.BooleanVar(value=False)

        def on_return(event):
            typed = self.text_area.get(start_index, "end-1c")
            outcome["value"] = typed if typed.strip() else default
            self.text_area.insert(tk.END, "\n")
            done.set(True)
            return "break"

        def on_escape(event):
            self.text_area.delete(start_index, "end-1c")
            outcome["value"] = None
            self.text_area.insert(tk.END, "\n")
            done.set(True)
            return "break"

        ret_id = self.text_area.bind("<Return>", on_return)
        esc_id = self.text_area.bind("<Escape>", on_escape)
        self.text_area.wait_variable(done)
        self.text_area.unbind("<Return>", ret_id)
        self.text_area.unbind("<Escape>", esc_id)
        self.text_area.see(tk.END)
        return outcome["value"]

    def _lister(self, func):
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        buffer = io.StringIO()
        func(self.selection, orientation=self.orientation.get(), out=buffer)
        self._afficher(buffer.getvalue())

    # ------------------------------------------------------------------
    # Selection donnees : selmes / selentete / effmes / initmes
    # ------------------------------------------------------------------

    def ouvrir_selection_dialog(self):
        """Equivalent GUI de `selmes` - saisie au clavier DANS la fenetre
        texte (console), pas de fenetre popup - questions posees les unes
        apres les autres comme le faisait le Fortran d'origine."""
        if not self.donnees:
            messagebox.showwarning(
                "No data",
                "Load a .ren file or import a MagIC folder first.",
            )
            return

        self.text_area.insert(tk.END, "\n--- Data selection (Escape to cancel) ---\n", "prompt")
        pattern = self._read_prefixed_pattern("Sample number (* = all, ? = wildcard): ")
        if pattern is None:
            return
        step_min_s = self._console_input("Step min: ", "0")
        if step_min_s is None:
            return
        step_max_s = self._console_input("Step max: ", "9999")
        if step_max_s is None:
            return
        demag_s = self._console_input("Demag code (e.g. N0, T, AF, * = all): ", "*")
        if demag_s is None:
            return
        try:
            step_min = int(step_min_s or 0)
            step_max = int(step_max_s or 9999)
        except ValueError:
            messagebox.showerror("Error", "Step min / Step max must be integers.")
            return

        demag1, demag2 = self._parse_demag(demag_s)
        new_matches = select_samples(
            self.donnees,
            pattern=pattern or "*",
            step_min=step_min,
            step_max=step_max,
            demag1=demag1,
            demag2=demag2,
            verbose=False,
        )
        # equivalent de selmes (dataselect.f) : accumule sur la selection
        # existante (nbech n'est jamais remis a zero dans selmes lui-meme,
        # seul initmes/"Init list" le fait) - ne remplace pas self.selection.
        self.selection = self.selection + new_matches
        nb_mesures = sum(len(s.mesures) for s in self.selection)
        self.text_area.insert(
            tk.END, f"Selection: +{len(new_matches)} sample(s) - "
                    f"total {len(self.selection)} sample(s), {nb_mesures} measurement(s)\n")
        self.text_area.see(tk.END)

    def _pick_from_list(self, title, items, allow_all=True):
        """Fenetre modale avec une Listbox pour choisir UN element a la
        souris (double-clic ou bouton Select) - PAS dans le Fortran
        (console texte uniquement), demande explicite de l'utilisateur
        pour remplacer la saisie au clavier du nom de site. Retourne
        l'element choisi, '*' si `allow_all` et l'entree "* (all)" est
        choisie, ou None si annule (bouton Cancel / fermeture)."""
        top = tk.Toplevel(self.root)
        top.title(title)
        top.transient(self.root)
        top.grab_set()

        listbox = tk.Listbox(top, activestyle="dotbox", exportselection=False)
        if allow_all:
            listbox.insert(tk.END, "* (all)")
        for item in items:
            listbox.insert(tk.END, item)
        listbox.selection_set(0)
        listbox.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        result = {"value": None}

        def confirm(_event=None):
            sel = listbox.curselection()
            if sel:
                text = listbox.get(sel[0])
                result["value"] = "*" if allow_all and text == "* (all)" else text
            top.destroy()

        def cancel():
            top.destroy()

        listbox.bind("<Double-Button-1>", confirm)
        listbox.bind("<Return>", confirm)

        btn_frame = ttk.Frame(top)
        btn_frame.pack(fill=tk.X, padx=8, pady=(0, 8))
        ttk.Button(btn_frame, text="Select", command=confirm).pack(side=tk.RIGHT, padx=(4, 0))
        ttk.Button(btn_frame, text="Cancel", command=cancel).pack(side=tk.RIGHT)

        top.protocol("WM_DELETE_WINDOW", cancel)
        listbox.focus_set()
        self.root.wait_window(top)
        return result["value"]

    def ouvrir_selection_site_dialog(self):
        """Pas dans le Fortran (aucune selection par site n'existe dans le
        menu d'origine) : selectionne tous les echantillons d'un site
        MagIC (champ `magic_site`, decode depuis la ligne roche a
        l'ouverture du fichier - voir testlect.decode_roche). Accumule sur
        self.selection comme `selmes` (voir ouvrir_selection_dialog),
        plutot que de la remplacer."""
        if not self.donnees:
            messagebox.showwarning(
                "No data",
                "Load a .ren file or import a MagIC file first.",
            )
            return

        sites = sorted({p.magic_site.strip() for p in self.donnees if p.magic_site.strip()})
        if not sites:
            messagebox.showwarning(
                "No site", "No sample has a decoded MagIC site (see the « roche » line).")
            return

        self.text_area.insert(tk.END, "\n--- Selection by site (Escape to cancel) ---\n", "prompt")
        site = self._pick_from_list("Select site", sites)
        if site is None:
            return
        step_min_s = self._console_input("Step min: ", "0")
        if step_min_s is None:
            return
        step_max_s = self._console_input("Step max: ", "9999")
        if step_max_s is None:
            return
        demag_s = self._console_input("Demag code (e.g. N0, T, AF, * = all): ", "*")
        if demag_s is None:
            return
        try:
            step_min = int(step_min_s or 0)
            step_max = int(step_max_s or 9999)
        except ValueError:
            messagebox.showerror("Error", "Step min / Step max must be integers.")
            return

        demag1, demag2 = self._parse_demag(demag_s)
        new_matches = select_samples_by_site(
            self.donnees, site=site or "*", step_min=step_min, step_max=step_max,
            demag1=demag1, demag2=demag2, verbose=False,
        )
        self.selection = self.selection + new_matches
        nb_mesures = sum(len(s.mesures) for s in self.selection)
        self._afficher(
            f"Site selection « {site} »: +{len(new_matches)} sample(s) - "
            f"total {len(self.selection)} sample(s), {nb_mesures} measurement(s)"
        )

    def ouvrir_entete_dialog(self):
        """Equivalent GUI de `selentete` : une saisie VIDE efface le
        header (lentete=0 dans le Fortran - dataselect.f:1295, blanc ==
        pas de prefixe actif), plutot que de re-appliquer l'ancienne
        valeur. Bug reel corrige : `self.entete` etait utilise comme
        valeur par defaut, rendant le header impossible a reinitialiser
        (Entree seul renvoyait toujours l'ancienne valeur)."""
        current = f" (current: '{self.entete}')" if self.entete else ""
        entete = self._console_input(f"Header{current}, blank to clear: ", "")
        if entete is None:
            return
        self.entete = entete.upper()[:12]
        self._afficher(
            f"Header set: '{self.entete}'\n" if self.entete else "Header cleared.\n"
        )

    def _read_prefixed_pattern(self, prompt):
        """Lit un motif echantillon en appliquant le header actif comme
        VRAI prefixe, fidele a `numero=entete(1:lentete)//chaine(1:len(chaine))`
        (dataselect.f, `selmes`/`selres` branche 'd') : le header est
        concatene devant TOUT ce qui est tape, y compris une saisie vide -
        dans ce cas le motif final est le header seul, pas '*'. Avant ce
        correctif, `self.entete` n'etait utilise que comme valeur par
        defaut (uniquement si la saisie etait vide), donc taper un motif
        n'appliquait jamais le prefixe : bug rapporte par l'utilisateur
        (header 'TRO' sans effet, obligeant a taper 'TRO14' en entier)."""
        label = f"{prompt}{self.entete}" if self.entete else prompt
        raw = self._console_input(label, "")
        if raw is None:
            return None
        if self.entete:
            return self.entete + raw
        return raw if raw.strip() else "*"

    def ouvrir_effmes_dialog(self):
        """Equivalent GUI de `effmes`."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- Delete measurements (Escape to cancel) ---\n", "prompt")
        pattern = self._console_input("Sample to erase (* = all): ", "*")
        if pattern is None:
            return
        step_min_s = self._console_input("Step min: ", "0")
        if step_min_s is None:
            return
        step_max_s = self._console_input("Step max: ", "9000")
        if step_max_s is None:
            return
        demag_s = self._console_input("Demag code (* = all): ", "*")
        if demag_s is None:
            return
        occurrence = self._console_input(
            "Occurrence to delete (* = all matches, 1 = first only, 2 = second only...): ", "*")
        if occurrence is None:
            return
        try:
            step_min = int(step_min_s or 0)
            step_max = int(step_max_s or 9000)
        except ValueError:
            messagebox.showerror("Error", "Step min / Step max must be integers.")
            return

        demag1, demag2 = self._parse_demag(demag_s)
        self.selection = delete_measurements(
            self.selection,
            pattern=pattern or "*",
            step_min=step_min,
            step_max=step_max,
            demag1=demag1,
            demag2=demag2,
            occurrence=occurrence or "*",
            verbose=False,
        )
        nb_mesures = sum(len(s.mesures) for s in self.selection)
        self.text_area.insert(
            tk.END, f"Selection after deletion: {len(self.selection)} sample(s), {nb_mesures} measurement(s)\n")
        self.text_area.see(tk.END)

    def reinitialiser_selection(self):
        """Equivalent GUI de `initmes`."""
        self.selection = init_selection()
        self._afficher("Selection reset - no sample selected.\n")

    # ------------------------------------------------------------------
    # Selection donnees : lismes / listeXYZ / lismesVRM / lismesdepth / infoech
    # ------------------------------------------------------------------

    def lister_mesures(self):
        self._lister(list_measurements)

    def lister_xyz(self):
        self._lister(list_xyz)

    def lister_vrm(self):
        self._lister(list_measurements_vrm)

    def ouvrir_lismesdepth_dialog(self):
        """Equivalent GUI de `lismesdepth`."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        depth_path = filedialog.askopenfilename(
            title="Depth file (2 columns: sample_name  depth)",
            filetypes=[("Text", "*.txt"), ("All files", "*.*")],
        )
        if not depth_path:
            return

        depths = {}
        try:
            with open(depth_path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            depths[parts[0]] = float(parts[1])
                        except ValueError:
                            continue
        except Exception as e:
            messagebox.showerror("Error", f"Could not read the depth file:\n{e}")
            return

        self.text_area.insert(tk.END, "\n--- Expected site direction (Escape to cancel) ---\n", "prompt")
        dec_s = self._console_input("Expected declination (D): ", "0.0")
        if dec_s is None:
            return
        inc_s = self._console_input("Expected inclination (I): ", "0.0")
        if inc_s is None:
            return
        try:
            expected_dec = float(dec_s or 0.0)
            expected_inc = float(inc_s or 0.0)
        except ValueError:
            messagebox.showerror("Error", "Declination / Inclination must be numbers.")
            return

        buffer = io.StringIO()
        list_measurements_depth(
            self.selection, depths,
            expected_dec=expected_dec, expected_inc=expected_inc,
            orientation=self.orientation.get(), out=buffer,
        )
        self._afficher(buffer.getvalue())

    def afficher_info_echantillons(self):
        """Equivalent GUI de `infoech`."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self._afficher(sample_info(self.selection))

    # ------------------------------------------------------------------
    # Calcul : ajuslig / fishmes / fishres / lisres / initres
    # ------------------------------------------------------------------

    def ouvrir_ajuslig_dialog(self):
        """Equivalent GUI de `ajuslig` (ajustement de droite par ACP,
        linesplans.f) : boucle sur TOUS les echantillons de la selection,
        comme le `do i=1,nbech` d'origine - zijderveld affiche pour chaque
        echantillon, premier step a 0 pour passer au suivant sans ajuster,
        puis apres chaque ajustement : sauver/refaire/suivant (Y/r/n), au
        lieu de traiter un seul echantillon choisi a l'avance."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        samples = [s for s in self.selection if len(s.mesures) >= 2]
        if not samples:
            messagebox.showwarning(
                "Not enough measurements", "No selected sample has at least 2 measurements.")
            return

        for ech in samples:
            self._current_graphic = ("zijderveld", ech.id)
            self._refresh_current_graphic()

            steps_text = "  ".join(
                f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
            self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")

            while True:  # boucle "refaire" (r)
                jdeb_s = self._console_input(
                    f"[{ech.id}] First step (0 = skip this sample): ", "1")
                if jdeb_s is None:
                    return
                try:
                    jdeb = int(jdeb_s)
                except ValueError:
                    messagebox.showerror("Error", "First step must be an integer.")
                    continue
                if jdeb == 0:
                    break  # echantillon suivant, sans ajustement

                jfin_s = self._console_input("Last step: ", str(len(ech.mesures)))
                if jfin_s is None:
                    return
                ancr_s = self._console_input("Anchored to origin (Y/n): ", "Y")
                if ancr_s is None:
                    return
                numcomp_s = self._console_input("Component number (1-9): ", "1")
                if numcomp_s is None:
                    return
                try:
                    jfin = int(jfin_s)
                    numcomp = int(numcomp_s or 1)
                except ValueError:
                    messagebox.showerror("Error", "Last step and component number must be integers.")
                    continue
                if not (1 <= jdeb <= jfin <= len(ech.mesures)):
                    messagebox.showerror("Error", "Invalid step range.")
                    continue
                anchored = ancr_s.strip().lower() != "n"

                fit = fit_line(ech, jdeb, jfin, anchored=anchored, numcomp=numcomp)
                if fit is None:
                    messagebox.showwarning(
                        "Fit rejected",
                        "MAD > 15° or non-linear trend (linearity test failed).",
                    )
                    continue

                dec, inc = _correct_dec_inc(fit, self.orientation.get())
                self._afficher(
                    f"{ech.id}: dec={dec:.1f}  inc={inc:.1f}  mad={fit.mad:.1f}  "
                    f"nb points={fit.nb}  ({'anchored' if anchored else 'not anchored'})\n"
                )
                # affiche le zijderveld AVEC cet ajustement superpose (pas
                # encore sauvegarde - ajoute temporairement a self.results,
                # que _refresh_current_graphic utilise deja pour filtrer les
                # ajustements a superposer par echantillon)
                self.results.append(fit)
                self._current_graphic = ("zijderveld", ech.id)
                self._refresh_current_graphic()

                decision = self._console_input(
                    "Save (Y) / redo (r) / next sample (n): ", "Y")
                if decision is None:
                    self.results.remove(fit)
                    return
                decision = decision.strip().lower()
                if decision == "r":
                    self.results.remove(fit)
                    continue  # redo pour le meme echantillon
                if decision != "n":
                    self._archive_only(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                else:
                    self.results.remove(fit)
                break  # echantillon suivant

    def ouvrir_autointerpretation_dialog(self):
        """Pas dans le Fortran (aucun equivalent) - demande explicite
        utilisateur ("une routine intelligente qui ferait des
        interpretations des diagrammes de desaimantation"), deuxieme des
        deux routines demandees (voir auto_interpretation.py pour
        l'algorithme - PCA en extension gloutonne, PAS un modele IA/ML).

        Boucle sur la selection, propose jusqu'a 2 composantes par
        echantillon (primary/secondary). Les DEUX sont converties en
        FitResult et superposees ENSEMBLE sur le Zijderveld AVANT toute
        question de sauvegarde - demande explicite utilisateur ("is it
        possible to visualize the interpretation on the zijderveld before
        to accept it") - puis seule(s) celle(s) choisie(s) sont
        conservees/archivees, les autres retirees de self.results.
        AUCUNE suggestion n'est archivee automatiquement (limite assumee
        et documentee dans auto_interpretation.py : fiable sur une
        decroissance a composante unique, faillible sur des composantes
        qui se chevauchent - chaque suggestion reste a valider)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        self.text_area.insert(tk.END, "\n--- Auto-interpretation (suggestions, Escape to stop) ---\n", "prompt")
        for ech in self.selection:
            if len(ech.mesures) < 2:
                continue
            suggestions = propose_components(ech)
            if not suggestions:
                self._afficher(format_suggestions(ech.id, suggestions))
                continue

            # Construit un FitResult PROVISOIRE par suggestion et les
            # superpose TOUS ensemble sur le Zijderveld avant de demander
            # quoi que ce soit (numcomp=1/2 pour que les deux se
            # distinguent visuellement, meme convention que ajuslig).
            candidates = {}  # label -> FitResult
            for s in suggestions:
                jdeb = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == s.step_first), None)
                jfin = next((i + 1 for i, m in enumerate(ech.mesures) if m.etape == s.step_last), None)
                if jdeb is None or jfin is None:
                    continue
                numcomp = 1 if s.label == "primary" else 2
                if s.kind == "plane":
                    fit = fit_plane(ech, jdeb, jfin, numcomp=numcomp)
                else:
                    fit = fit_line(ech, jdeb, jfin, anchored=s.anchored, numcomp=numcomp)
                if fit is not None:
                    candidates[s.label] = fit

            self._afficher(format_suggestions(ech.id, suggestions))
            if not candidates:
                self._afficher(f"{ech.id}: suggested interval(s) no longer fit (data changed?).\n")
                continue

            self.results.extend(candidates.values())
            self._current_graphic = ("zijderveld", ech.id)
            self._refresh_current_graphic()

            choice = self._console_input(
                "Save which (p=primary, s=secondary, b=both, Enter=discard all): ", "")
            if choice is None:
                for fit in candidates.values():
                    self.results.remove(fit)
                return
            choice = choice.strip().lower()
            keep_labels = set()
            if choice == "b":
                keep_labels = set(candidates.keys())
            elif choice == "p" and "primary" in candidates:
                keep_labels = {"primary"}
            elif choice == "s" and "secondary" in candidates:
                keep_labels = {"secondary"}

            for label, fit in candidates.items():
                if label in keep_labels:
                    self._archive_only(fit)
                else:
                    self.results.remove(fit)
            if keep_labels:
                self._afficher(f"Saved: {', '.join(sorted(keep_labels))}  (total results: {len(self.results)})\n")

    def ouvrir_ajusplans_dialog(self):
        """Equivalent GUI de `ajusplans` (linesplans.f) : ajustement de plan
        (grand cercle) par ACP, meme boucle pause/decision que ajuslig. Pas
        de prompt d'ancrage (le Fortran ancre toujours a l'origine pour un
        plan, aucun choix propose a l'utilisateur)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        samples = [s for s in self.selection if len(s.mesures) >= 2]
        if not samples:
            messagebox.showwarning(
                "Not enough measurements", "No selected sample has at least 2 measurements.")
            return

        for ech in samples:
            self._current_graphic = ("zijderveld", ech.id)
            self._refresh_current_graphic()

            steps_text = "  ".join(
                f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
            self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")

            while True:  # boucle "refaire"
                jdeb_s = self._console_input(
                    f"[{ech.id}] First step (0 = skip this sample): ", "1")
                if jdeb_s is None:
                    return
                try:
                    jdeb = int(jdeb_s)
                except ValueError:
                    messagebox.showerror("Error", "Must be an integer.")
                    continue
                if jdeb == 0:
                    break

                jfin_s = self._console_input("Last step: ", str(len(ech.mesures)))
                if jfin_s is None:
                    return
                try:
                    jfin = int(jfin_s)
                except ValueError:
                    messagebox.showerror("Error", "Must be an integer.")
                    continue
                if not (1 <= jdeb <= jfin <= len(ech.mesures)):
                    messagebox.showerror("Error", "Invalid step range.")
                    continue

                norm_s = self._console_input("Normalize (Y/n): ", "Y")
                if norm_s is None:
                    return
                normalize = norm_s.strip().lower() != "n"

                numcomp_s = self._console_input("Component number (1-9): ", "1")
                if numcomp_s is None:
                    return
                try:
                    numcomp = int(numcomp_s or 1)
                except ValueError:
                    messagebox.showerror("Error", "Must be an integer.")
                    continue

                fit = fit_plane(ech, jdeb, jfin, normalize=normalize, numcomp=numcomp)
                if fit is None:
                    messagebox.showwarning("Fit rejected", "MAD > 25° (poorly defined plane).")
                    continue

                dec, inc = _correct_dec_inc(fit, self.orientation.get())
                self._afficher(
                    f"{ech.id} : pole dec={dec:.1f}  inc={inc:.1f}  mad={fit.mad:.1f}  "
                    f"nb points={fit.nb}\n"
                )

                decision = self._console_input(
                    "Save (Y) / redo (r) / next sample (n): ", "Y")
                if decision is None:
                    return
                decision = decision.strip().lower()
                if decision == "r":
                    continue
                if decision != "n":
                    self._save_result(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                break

    def ouvrir_ajusfisher_dialog(self):
        """Equivalent GUI de `ajusfisher` (linesplans.f) : moyenne de Fisher
        sur une plage d'etapes. Repli automatique sur une direction unique
        (branche 735 du Fortran) pour les echantillons a moins de 3 mesures,
        ou si l'utilisateur choisit un seul step (premier == dernier)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        for ech in self.selection:
            if not ech.mesures:
                continue

            if len(ech.mesures) < 3:
                fit = fit_single_direction(ech)
                dec, inc = _correct_dec_inc(fit, self.orientation.get())
                self._afficher(
                    f"{ech.id}: single direction (fewer than 3 measurements) dec={dec:.1f}  "
                    f"inc={inc:.1f}\n"
                )
                decision = self._console_input(
                    f"[{ech.id}] Save the single direction (Y/n): ", "Y")
                if decision is None:
                    return
                if decision.strip().lower() != "n":
                    self._save_result(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                continue

            steps_text = "  ".join(
                f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
            self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")

            while True:
                jdeb_s = self._console_input(
                    f"[{ech.id}] First step (0 = skip this sample): ", "1")
                if jdeb_s is None:
                    return
                try:
                    jdeb = int(jdeb_s)
                except ValueError:
                    messagebox.showerror("Error", "Must be an integer.")
                    continue
                if jdeb == 0:
                    break

                jfin_s = self._console_input("Last step: ", str(len(ech.mesures)))
                if jfin_s is None:
                    return
                try:
                    jfin = int(jfin_s)
                except ValueError:
                    messagebox.showerror("Error", "Must be an integer.")
                    continue
                if jfin == 0:
                    break
                if not (1 <= jdeb <= jfin <= len(ech.mesures)):
                    messagebox.showerror("Error", "Invalid step range.")
                    continue

                numcomp_s = self._console_input("Component number (1-9): ", "1")
                if numcomp_s is None:
                    return
                try:
                    numcomp = int(numcomp_s or 1)
                except ValueError:
                    messagebox.showerror("Error", "Must be an integer.")
                    continue

                if jfin == jdeb:
                    fit = fit_single_direction(ech, index=jdeb - 1)
                    dec, inc = _correct_dec_inc(fit, self.orientation.get())
                    self._afficher(f"{ech.id}: single direction dec={dec:.1f}  inc={inc:.1f}\n")
                else:
                    fit = fit_fisher_direction(ech, jdeb, jfin, numcomp=numcomp)
                    dec, inc = _correct_dec_inc(fit, self.orientation.get())
                    self._afficher(
                        f"{ech.id} : dec={dec:.1f}  inc={inc:.1f}  a95={fit.mad:.1f}  "
                        f"k={fit.tx[0]:.1f}  nb={fit.nb}\n"
                    )

                decision = self._console_input(
                    "Save (Y) / redo (r) / next sample (n): ", "Y")
                if decision is None:
                    return
                decision = decision.strip().lower()
                if decision == "r":
                    continue
                if decision != "n":
                    self._save_result(fit)
                    self._afficher(f"Result saved: {len(self.results)}\n")
                break

    def ouvrir_ajusligauto_dialog(self):
        """Equivalent GUI de `ajusligauto` : ajustement de droite automatique
        sur la totalite des mesures de chaque echantillon selectionne, meme
        ancrage pour tous, sauvegarde sans confirmation."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return

        ancr_s = self._console_input(
            "Automatic fit anchored to origin (Y/n): ", "Y")
        if ancr_s is None:
            return
        anchored = ancr_s.strip().lower() != "n"

        fits = fit_lines_auto(self.selection, anchored=anchored)
        self._save_results(fits)
        self._afficher(
            f"Automatic fit: {len(fits)} line(s) fitted and saved "
            f"({'anchored' if anchored else 'not anchored'}).\n"
        )

    def ouvrir_ajusligredo_dialog(self):
        """Equivalent GUI de `ajusligredo` : rejoue des ajustements (droite
        ou plan) a partir d'un fichier texte 'redo' (voir docstring de
        calcul.fit_from_redo_file pour le format des lignes). Comme le
        Fortran (`call seloneech(samnum)` a chaque ligne), chaque
        specimen designe dans le fichier redo est recherche et charge
        directement depuis self.donnees - PAS depuis self.selection, qui
        n'a donc pas besoin d'etre preparee au prealable."""
        if not self.donnees:
            messagebox.showwarning("No data", "Load a .ren file first.")
            return

        path = filedialog.askopenfilename(
            title="Redo file (sample L|P o|n tempmin tempmax numcomp)",
            filetypes=[("Text files", "*.txt *.redo"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                lines = f.readlines()
        except Exception as e:
            messagebox.showerror("Error", f"Could not read the file:\n{e}")
            return

        fits = fit_from_redo_file(self.donnees, lines)
        self._save_results(fits)
        self._afficher(
            f"Redo file « {os.path.basename(path)} »: {len(fits)} fit(s) replayed "
            f"and saved.\n"
        )

    def afficher_mdf(self):
        """Equivalent GUI de `mdf` (calcul.f:519-678)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self._afficher(list_mdf(self.selection))

    def afficher_mean_intensity(self):
        """Equivalent GUI de `mdi`/`mds` (calcul.f:683-865)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        mi = compute_mean_intensity(self.selection)
        if mi is None:
            messagebox.showerror(
                "Error",
                "Could not compute (mass/volume mix in the selection, "
                "or not enough measurements).",
            )
            return
        ms = compute_mean_susceptibility(self.selection)
        ech0 = self.selection[0]
        self._afficher(format_mean_intensity(mi, ms, lat=ech0.lat, rlong=ech0.rlong))

    def ouvrir_koenigsberger_dialog(self):
        """Equivalent GUI de `Koenigs` (calcul.f:3869-3926)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        valk_s = self._console_input("Reference field for the Koenigsberger ratio (µT): ", "40")
        if valk_s is None:
            return
        try:
            valk = float(valk_s)
        except ValueError:
            messagebox.showerror("Error", "Must be a number.")
            return
        if valk == 0.0:
            messagebox.showerror("Error", "The reference field cannot be zero.")
            return
        rows = compute_koenigsberger(self.selection, valk)
        self._afficher(
            format_koenigsberger(rows) if rows
            else "(no measurement with non-zero susceptibility in the selection)"
        )

    def afficher_mean_inclination(self):
        """Equivalent GUI de `meaninc` (calcul.f:1011-1129, estimateur de
        McFadden & Reid 1982)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        r = compute_mean_inclination(self.selection, orientation=self.orientation.get())
        if r is None:
            messagebox.showerror("Error", "Not enough measurements, or the iteration does not converge.")
            return
        self._afficher(format_mean_inclination(r))

    def afficher_diff_measurements(self):
        """Equivalent GUI de `diffmes` (dataselect.f:1156-1236)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        self._afficher(list_diff_measurements(self.selection, orientation=self.orientation.get()))

    def appliquer_viscosity_test(self):
        """Equivalent GUI de `viscos` (viscos.f) : MUTE les 2 premières
        mesures de chaque échantillon (N+/N-) en (vecteur moyen, vecteur
        différence) - operation irréversible sur la sélection en mémoire,
        confirmation demandée."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        if not messagebox.askyesno(
            "Viscosity test",
            "This replaces the first 2 measurements (N+/N-) of each selected "
            "sample with a mean vector and a difference vector - irreversible "
            "for this session. Continue?",
        ):
            return
        warnings = apply_viscosity_test(self.selection)
        msg = "Viscosity test applied."
        if warnings:
            msg += "\n" + "\n".join(warnings)
        self._afficher(msg + "\n")

    def ouvrir_subtraction_dialog(self):
        """Equivalent GUI de `soustra` (dataselect.f:1239-1283) : soustrait
        le vecteur d'une mesure de toutes les autres du même échantillon,
        puis supprime cette ligne - un seul échantillon a la fois, comme le
        Fortran."""
        if len(self.selection) != 1:
            messagebox.showwarning(
                "Invalid selection", "Select a single sample for the subtraction.")
            return
        ech = self.selection[0]
        steps_text = "  ".join(
            f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
        self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
        row_s = self._console_input("Line number to subtract: ", "")
        if row_s is None:
            return
        try:
            row = int(row_s)
        except ValueError:
            messagebox.showerror("Error", "Must be an integer.")
            return
        if not messagebox.askyesno(
            "Subtraction",
            f"Subtract line {row} of {ech.id} from all other measurements, "
            f"and remove it - irreversible for this session. Continue?",
        ):
            return
        err = apply_subtraction(ech, row)
        if err:
            messagebox.showerror("Error", err)
            return
        self._afficher(f"{ech.id}: line {row} subtracted and removed ({len(ech.mesures)} measurements remaining).\n")

    def ouvrir_holderarm_dialog(self):
        """Equivalent GUI de `holderarm` (calcul.f:2337-2414) : enregistre
        les 6 mesures ARM d'un porte-échantillon vide - utile pour une
        future correction dans Anisotropy (non encore porté)."""
        if len(self.selection) != 1:
            messagebox.showwarning(
                "Invalid selection",
                "Select the single sample corresponding to the empty holder.")
            return
        ech = self.selection[0]
        steps_text = "  ".join(
            f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
        self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
        idx_s = self._console_input("Line numbers x+ x- y+ y- z+ z-: ", "")
        if idx_s is None:
            return
        try:
            ixp, ixm, iyp, iym, izp, izm = (int(v) for v in idx_s.split())
        except ValueError:
            messagebox.showerror("Error", "Exactly 6 space-separated integers are required.")
            return
        bg = record_arm_holder(ech, ixp, ixm, iyp, iym, izp, izm)
        if bg is None:
            messagebox.showerror("Error", "One of the line numbers is invalid.")
            return
        self._arm_holder_background = bg
        self._afficher(
            "ARM holder recorded:\n"
            + "\n".join(f"{i + 1}: {bg.x[i]:.3e}  {bg.y[i]:.3e}  {bg.z[i]:.3e}" for i in range(6))
            + "\n"
        )

    def ouvrir_anisotropy_dialog(self):
        """Equivalent GUI de `anisoauto` (calcul.f:3951+), tenseur 'A0'
        uniquement (voir calcul.compute_anisotropy_tensor) : detection
        automatique des 6 positions X+/X-/Y+/Y-/Z+/Z- (Z+/Z- pouvant etre
        substituees par des mesures R/V a la meme etape que X, ex. 'RH'/
        'VH'), correction optionnelle par la ligne de base porte-
        echantillon (Holder_ARM...), ecriture dans le fichier .ANI (meme
        emplacement/format que celui deja lu par Inverse_ANI_correction...).
        Les 14 variantes jackknife (A+/A-/A1-A6/B1-B6) ne sont PAS portees
        (etape ulterieure).

        Affiche TOUTE la sequence de calcul dans la fenetre texte, comme
        le fait la console du Fortran d'origine (pas seulement le resultat
        final) - demande explicite de l'utilisateur : les donnees
        d'anisotropie sont souvent basees sur des TRM PARTIELLES (une part
        de la NRM n'a pas ete remplacee), et le tenseur BRUT avant
        symetrisation + les diagnostics par position sont le seul moyen de
        verifier qu'un echantillon etait bien oriente lors de l'acquisition
        de la TRM."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        if not self.results_path:
            messagebox.showwarning(
                "No file", "Load a .ren file first (the .ANI file is derived from it).")
            return
        ani_path = os.path.splitext(self.results_path)[0] + ".ANI"

        self.text_area.insert(tk.END, "\n--- Anisotropy (Escape to cancel) ---\n", "prompt")
        auto_s = self._console_input("Automated recognition of the 6 steps? y/N: ", "N")
        if auto_s is None:
            return
        use_auto = auto_s.strip().lower() == "y"

        use_zb = False
        if use_auto:
            use_zb_s = self._console_input(
                "Using ZB instead of R for a > +/- 5% evolution? y/N: ", "N")
            if use_zb_s is None:
                return
            use_zb = use_zb_s.strip().lower() == "y"

        done, skipped = [], []
        for ech in self.selection:
            lines = [f"\n--- {ech.id}: {'automated ' if use_auto else ''}calcul of TRM or ARM anisotropy with 6 positions ---"]
            lines.append(
                "Remanent magnetization list, normalized by mass (Am2/kg)"
                if ech.norme == "m" else
                "Remanent magnetization list (A/m)"
            )
            buffer = io.StringIO()
            list_measurements([ech], orientation=self.orientation.get(), out=buffer)
            lines.append(buffer.getvalue().rstrip("\n"))
            self._afficher("\n".join(lines) + "\n")
            lines = []

            manual_positions = None
            if not use_auto:
                # equivalent du repli manuel de `anisot` (carreco != "y") :
                # pas de substitution R/V->Z+/Z- dans cette branche du
                # Fortran (zplus/zminus n'y sont jamais lus), seulement les
                # codes litteraux X+/X-/Y+/Y-/Z+/Z-.
                idx_s = self._console_input(
                    f"[{ech.id}] Line numbers x+ x- y+ y- z+ z- (Escape to skip): ", "")
                if idx_s is None:
                    skipped.append(ech.id)
                    continue
                try:
                    idx = [int(v) for v in idx_s.split()]
                    if len(idx) != 6 or any(not (1 <= i <= len(ech.mesures)) for i in idx):
                        raise ValueError
                except ValueError:
                    messagebox.showerror("Error", "6 valid line numbers are required.")
                    skipped.append(ech.id)
                    continue
                manual_positions = dict(zip(
                    ("X+", "X-", "Y+", "Y-", "Z+", "Z-"),
                    (ech.mesures[i - 1] for i in idx),
                ))

            result = compute_anisotropy_tensor(
                ech, holder=self._arm_holder_background, use_zb_on_evolution=use_zb,
                positions=manual_positions,
            )
            if result is None:
                skipped.append(ech.id)
                lines.append("Decoding incomplete: could not identify all 6 positions (X+/X-/Y+/Y-/Z+/Z-).")
                self._afficher("\n".join(lines) + "\n")
                continue

            if result.holder_used:
                bg = self._arm_holder_background
                lines.append("Holder ARM background subtracted (X+,X-,Y+,Y-,Z+,Z-):")
                lines.extend(
                    f"  {i + 1}: {bg.x[i]:.3E}  {bg.y[i]:.3E}  {bg.z[i]:.3E}" for i in range(6)
                )

            if result.swapped_axes:
                lines.extend(
                    f"!! inversion {ax}+ and {ax}- for sample: {ech.id}" for ax in result.swapped_axes
                )

            unit = "Am2/kg" if ech.norme == "m" else "A/m"
            lines.append(
                f"Mean NRM (residual, subtracted below): "
                f"x={result.nrm_mean[0]:.3E}  y={result.nrm_mean[1]:.3E}  z={result.nrm_mean[2]:.3E}"
            )
            lines.append("Mean NRM subtracted - TRM values:")
            lines.extend(
                f"  {d.key}: {d.measurement.etape}{d.measurement.cod1}{d.measurement.cod2}  "
                f"{d.intensity:10.3E} {unit}  dec={d.dec:6.1f}  inc={d.inc:6.1f}"
                for d in result.position_diags
            )
            lines.append(f"Deviation (pair asymmetry / mean TRM intensity): {result.deviation_pct:5.1f} %")

            if result.trm_evolution_pct is not None:
                lines.append(f"ZB check found - TRM evolution: {result.trm_evolution_pct:5.1f} %")
                if result.zb_used:
                    lines.append("  -> evolution exceeds 5%: Z+ replaced by the ZB control measurement.")
            else:
                lines.append("No ZB control measurement found (no pTRM-check available).")

            lines.append("-------- Tensor (raw, before symmetrization) --------")
            lines.extend(f"  {r[0]:10.3E}  {r[1]:10.3E}  {r[2]:10.3E}" for r in result.raw)
            lines.append("(forcing symmetric)")

            t = result.tensor
            lines.append("-------- Tensor 'A0' (symmetric) --------")
            lines.append(f"  {t.k11:10.3E}  {t.k12:10.3E}  {t.k13:10.3E}")
            lines.append(f"  {t.k12:10.3E}  {t.k22:10.3E}  {t.k23:10.3E}")
            lines.append(f"  {t.k13:10.3E}  {t.k23:10.3E}  {t.k33:10.3E}")

            write_ani_tensor(
                ani_path, ech, t, positions=result.positions,
                trm_evolution_pct=result.trm_evolution_pct or 0.0,
                deviation_pct=result.deviation_pct,
            )
            done.append(ech.id)
            self._afficher("\n".join(lines) + "\n")

        summary = f"Tensor 'A0' written for {len(done)} sample(s) -> {ani_path}\n"
        if skipped:
            summary += f"Skipped (6 positions not identified): {', '.join(skipped)}\n"
        self._afficher(summary)

    def ouvrir_inverseani_dialog(self):
        """Equivalent GUI de `inverseani` (plotpaleoint2.f:1803-1918) :
        applique l'inverse d'un tenseur d'anisotropie (lu dans un fichier
        .ANI, dérivé du fichier de données comme le fichier .r) à toutes
        les mesures de l'échantillon sélectionné."""
        if len(self.selection) != 1:
            messagebox.showwarning(
                "Invalid selection", "Select a single sample.")
            return
        ech = self.selection[0]
        if getattr(ech, "flaganiso", False):
            messagebox.showwarning("Already corrected", f"{ech.id} has already been corrected for anisotropy.")
            return
        if not self.results_path:
            messagebox.showwarning(
                "No file", "Load a .ren file first (the .ANI file is derived from it).")
            return
        ani_path = os.path.splitext(self.results_path)[0] + ".ANI"
        if not os.path.exists(ani_path):
            ani_path = os.path.splitext(results_path_for(self.results_path))[0] + ".ANI"
        choice = self._console_input(
            "Correction with 1: TRM tensor (A0)  2: ARM tensor (F0)  3: susceptibility (N0): ", "1")
        if choice is None:
            return
        try:
            ichoice = int(choice)
        except ValueError:
            ichoice = 1
        if ichoice not in (1, 2, 3):
            ichoice = 1
        code2 = _ANI_CODE2[ichoice]

        tensor = read_ani_tensor(ani_path, ech.id, code2)
        if tensor is None:
            messagebox.showerror(
                "Error", f"No '{code2}' tensor for {ech.id} in {ani_path}.")
            return
        apply_inverse_anisotropy(ech, tensor)
        ech.flaganiso = True
        self._afficher(f"{ech.id}: measurements corrected (inverse of the {code2} tensor).\n")

    def ouvrir_cooling_rate_dialog(self):
        """Equivalent GUI de `vitref` (calcul.f:3461-3625, branche "live") :
        boucle sur la sélection, détection automatique du motif L/Q + R/V
        (sinon saisie manuelle des 5 numéros de ligne)."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        auto_s = self._console_input("Automated procedure (Y/n): ", "Y")
        if auto_s is None:
            return
        auto = auto_s.strip().lower() != "n"

        for ech in self.selection:
            rows = detect_cooling_rate_rows(ech) if auto else None
            if rows is None:
                steps_text = "  ".join(
                    f"{i + 1}:{m.etape}{m.cod1}{m.cod2}" for i, m in enumerate(ech.mesures))
                self._afficher(f"{ech.id} - available steps:\n{steps_text}\n")
                if auto:
                    self._afficher(f"{ech.id}: L/Q + R/V pattern not found automatically.\n")
                rows_s = self._console_input(
                    f"[{ech.id}] Line numbers ATR+fast ATR-fast ATR+slow before/after loop "
                    "(empty = skip): ", "")
                if rows_s is None:
                    return
                if not rows_s.strip():
                    continue
                try:
                    rows = tuple(int(v) for v in rows_s.split())
                    if len(rows) != 5:
                        raise ValueError
                except ValueError:
                    messagebox.showerror("Error", "Exactly 5 integers are required.")
                    continue

            result = compute_cooling_rate(ech, *rows)
            self._afficher(format_cooling_rate(result) + "\n")

    @staticmethod
    def _format_fisher(stats, title):
        return (
            f"{title}\n"
            "----------------------------------------\n"
            f"n = {stats.n}\n"
            f"Dec = {stats.dec:.1f}   Inc = {stats.inc:.1f}\n"
            f"R = {stats.r:.3f}   k = {stats.k:.1f}\n"
            f"a95 = {stats.a95:.1f}   csd = {stats.csd:.1f}\n"
        )

    def fisher_mesures(self):
        """Equivalent GUI de `fishmes`."""
        if not self.selection:
            messagebox.showwarning("No selection", "Select some samples first.")
            return
        try:
            stats = fisher_from_measurements(self.selection, orientation=self.orientation.get())
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return
        self._afficher(self._format_fisher(stats, "Fisher on the selected measurements"))

    def fisher_resultats(self):
        """Equivalent GUI de `fishres` (limité aux résultats de type ligne)."""
        if not self.results:
            messagebox.showwarning("No results", "Run one or more line fits first.")
            return
        try:
            stats = fisher_from_results(self.results, orientation=self.orientation.get())
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return
        self._afficher(self._format_fisher(stats, "Fisher on the results (lines)"))

    def lister_resultats(self):
        """Equivalent GUI de `lisres` (limité aux résultats de type ligne).
        Passe l'orientation courante : les dec/inc affichés sont recalculés
        dans ce repère à chaque appel (voir list_results), pas figés dans
        le repère échantillon d'origine."""
        if not self.results:
            messagebox.showwarning("No results", "Run one or more line fits first.")
            return
        self._afficher(list_results(self.results, orientation=self.orientation.get(),
                                     donnees=self.donnees))

    def evaluer_interpretations(self):
        """Pas dans le Fortran (aucun equivalent) - demande explicite
        utilisateur ("une routine qui permettrait d'evaluer les
        interpretations"). Audite self.results (memes resultats deja
        charges/produits que List results) : recalcule chaque fit ligne/
        plan depuis les mesures vivantes (self.donnees) et note la
        qualite (voir interpretation_quality.py pour le detail des
        criteres - MAD recalcule, angle ancre/libre pour les droites,
        ratio de linearite). Les moyennes de site ("mean:") sont ignorees
        (pas un fit ligne/plan)."""
        if not self.results:
            messagebox.showwarning("No results", "Run one or more line/plane fits first, or load some via Select results...")
            return
        reports = evaluate_results(self.results, self.donnees)
        self._afficher(format_quality_report(reports) + "\n")

    def reinitialiser_resultats(self):
        """Equivalent GUI de `initres`."""
        self.results = init_results()
        self._afficher("Results reset - no result saved.\n")

    def ouvrir_selres_dialog(self):
        """Equivalent GUI de `selres` (dataselect.f) : charge des resultats
        DEPUIS le fichier .r (equivalent filr) dans self.results. Trois
        modes, comme le Fortran (`carselect` : "Data (default), mean (m) ou
        site [data+mean] (s)") :
        - Data : resultats normaux (L/P/f/s), filtres par echantillon/type/
          composante - les moyennes "mean:" sont exclues.
        - Mean : uniquement les moyennes de site, filtrees par orientation
          courante (equivalent `res.par3==float(iorient)`) - pas de filtre
          type/composante (non demandes par le Fortran dans ce mode).
        - Site : la/les moyenne(s) matchee(s) PLUS les resultats individuels
          qui la composent (equivalent `decodelisteres`)."""
        if not self.results_path or not os.path.exists(self.results_path):
            messagebox.showwarning(
                "No results file",
                "No .r file found for the loaded data "
                "(run a fit first, or load a .ren file)."
                if not self.results_path else
                f"{self.results_path} does not exist yet (no result archived).",
            )
            return

        self.text_area.insert(tk.END, "\n--- Select results (Escape to cancel) ---\n", "prompt")
        carselect = self._console_input(
            "Data (default), mean (m) or site [data+mean] (s): ", "d")
        if carselect is None:
            return
        carselect = (carselect.strip().lower() or "d")[:1]

        if carselect in ("m", "s"):
            site = self._console_input("Site (name, without « mean: », * = all): ", "*")
            if site is None:
                return
            loaded = load_results(
                self.results_path, pattern=site or "*", carselect=carselect,
                iorient=self.orientation.get(),
            )
            if not loaded:
                available = available_mean_orientations(self.results_path, site or "*")
                if available and self.orientation.get() not in available:
                    noms = {1: "Sample (CE)", 2: "In situ (IS)", 3: "Tilt cor. (CP)"}
                    self._afficher(
                        "No mean found for the current orientation "
                        f"« {noms.get(self.orientation.get(), self.orientation.get())} » - "
                        "but some exist for: "
                        + ", ".join(noms.get(o, str(o)) for o in available)
                        + ". Change the orientation (Pmag data menu) then try again.\n"
                    )
        else:
            pattern = self._read_prefixed_pattern("Sample (* = all): ")
            if pattern is None:
                return
            cat1 = self._console_input("Type (L/P/f/s, * = all): ", "*")
            if cat1 is None:
                return
            numcomp_s = self._console_input("Component number (empty = all): ", "")
            if numcomp_s is None:
                return
            numcomp = None
            if numcomp_s.strip():
                try:
                    numcomp = int(numcomp_s)
                except ValueError:
                    messagebox.showerror("Error", "The component number must be an integer.")
                    return
            loaded = load_results(self.results_path, pattern=pattern or "*",
                                   carselect="d", cat1=cat1 or "*", numcomp=numcomp)

        # tx/ty/tz (segment ajuste, pour le trace sur un Zijderveld) ne sont
        # plus stockes dans le fichier .r (voir calcul.recompute_fit_geometry)
        # - recalcules ici a partir des mesures brutes (self.donnees), une
        # seule fois au chargement plutot qu'a chaque affichage - demande
        # explicite utilisateur ("to draw the line on the zijderveld plot we
        # can redo the calculation").
        self.results = [recompute_fit_geometry(r, self.donnees) for r in loaded]
        self._afficher(f"nb selected results: {len(loaded)}\n")


if __name__ == "__main__":
    root = tk.Tk()
    app = StarmacApp(root)
    root.mainloop()
