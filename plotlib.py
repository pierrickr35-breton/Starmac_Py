"""
Port Python de la bibliotheque de trace Fortran (graphicsAWE.f95 : `plot`,
`symbol`/`label2`, `newpen`, `thickn`, `plottxt`, `number`, `circl2`,
`clear`...). graphic.f95 (QuickDraw, annees 80) documente l'intention
d'origine ; graphicsAWE.f95 est la version reellement compilee dans
Starmac_AWE_v22 (AWE_Canvas/Qt5) et sert de reference pour ce portage.

Objectif : conserver la MEME API imperative a base de "crayon" (position
courante, origine deplacable, coordonnees en cm) pour pouvoir traduire les
routines de trace Fortran restantes (stereo, susceptibilite, XY,
paleointensite...) quasi ligne a ligne, plutot que redevoir la mise en page
de chacune en idiomes matplotlib. Le rendu final passe par une Figure/Axes
matplotlib (independante de la resolution : l'unite de travail est le pouce,
sans le facteur `scrdpi` du Fortran qui n'a de sens que pour un canvas en
pixels physiques).

Non porte (delibbrement hors perimetre, cf. zijderveld.py) : l'export
PostScript/SVG (`postscript`, `plotsvg`, `laser*`) et la mise en page
multi-panneaux fixe de `boite()` (specifique a la fenetre AWE d'origine) -
chaque graphique aura ici sa propre Figure matplotlib independante ;
l'export SVG passe par `fig.savefig(path, format="svg")`, que matplotlib
sait deja faire correctement a partir des memes appels PlotContext.
"""

import math
from typing import List, Optional, Tuple

from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, PathPatch, Polygon
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

_FONT = FontProperties(family="Helvetica")
# Hauteur de capitale (cap-height) d'une TextPath a size=1, mesuree une fois
# empiriquement (police par defaut si Helvetica absente) - sert a convertir
# `height_cm` (hauteur de caractere voulue) en facteur d'echelle du glyphe.
_CAP_HEIGHT = TextPath((0, 0), "H", size=1, prop=_FONT).get_extents().height
# Largeur (en unites de cap-height) d'un chiffre, mesuree une fois de la
# meme facon que _CAP_HEIGHT - sert a decaler du texte de "la largeur d'un
# caractere" pour une hauteur de police donnee (ex. label d'etape Zijderveld).
_DIGIT_WIDTH_RATIO = TextPath((0, 0), "0", size=1, prop=_FONT).get_extents().width / _CAP_HEIGHT


def char_width_cm(height_cm: float) -> float:
    """Largeur approximative (cm) d'un caractere chiffre rendu par
    `plottxt`/`number` a la hauteur `height_cm`."""
    return _DIGIT_WIDTH_RATIO * height_cm


# newpen() : palette 1-9 (graphicsAWE.f95)
_PEN_COLORS = {
    1: "black", 2: "white", 3: "red", 4: "green", 5: "blue",
    6: "cyan", 7: "violet", 8: "yellow", 9: "deeppink",
}

# label2() : types de symboles (voir commentaire d'origine dans graphicsAWE.f95)
#   1/5: plus  2: point  3: barre verticale  4: barre horizontale
#   6: croix  7: asterisque  8: cercle  9: carre  10: triangle  11: peck
#   12: losange  13: etoile  14: cercle plein  15: carre plein
#   16: triangle plein  17: peck plein  18: losange plein  19: etoile pleine
#   20: cercle blanc (rempli de blanc)


class PlotContext:
    """Equivalent de plot()/symbol()/newpen()/thickn()/plottxt()/number()/
    circl2()/clear() : dessine sur une Axes matplotlib fournie par l'appelant.
    Coordonnees en cm, comme dans le Fortran d'origine (`plot(x1,y1,mode)`).
    """

    def __init__(self, ax: Axes):
        self.ax = ax
        self.origin_x = 0.0
        self.origin_y = 0.0
        self.x = 0.0
        self.y = 0.0
        self.rotation = 0.0  # degres (equivalent scrrot)
        self.pen_color = "black"
        self.fill_color = "black"
        self.line_width = 1.0

    # ------------------------------------------------------------------
    # Etat du crayon : newpen / newpencol / thickn / clear
    # ------------------------------------------------------------------

    def newpen(self, ipen: int) -> None:
        """Equivalent de `newpen(ipen)` : palette de couleurs 1-9."""
        self.pen_color = _PEN_COLORS.get(ipen, "black")
        self.fill_color = self.pen_color

    def newpencol(self, sr: int, sg: int, sb: int, fr: int, fg: int, fb: int) -> None:
        """Equivalent de `newpencol(...)` : couleurs RGB (0-255) explicites
        pour le trait (stroke) et le remplissage (fill)."""

        def hexcol(r: int, g: int, b: int) -> str:
            if r < 0 or g < 0 or b < 0 or r > 255 or g > 255 or b > 255:
                return "black"
            return f"#{r:02x}{g:02x}{b:02x}"

        self.pen_color = hexcol(sr, sg, sb)
        self.fill_color = hexcol(fr, fg, fb)

    def thickn(self, sizel: float) -> None:
        """Equivalent de `thickn(sizel)`."""
        self.line_width = sizel * 2

    def clear(self) -> None:
        """Equivalent de `clear()` (sans le fond `boite()`, specifique a la
        fenetre AWE d'origine - chaque Figure matplotlib gere son propre
        fond)."""
        self.ax.clear()
        self.ax.set_aspect("equal", adjustable="datalim")
        self.ax.axis("off")
        self.x = self.y = self.origin_x = self.origin_y = 0.0

    # ------------------------------------------------------------------
    # Transformation de coordonnees (equivalent du bloc commun a plot/
    # plottxt/circl2 : rotation scrrot, origine scrhor/scrvor, cm -> pouces)
    # ------------------------------------------------------------------

    def _transform(self, x_cm: float, y_cm: float) -> Tuple[float, float]:
        x_in, y_in = x_cm / 2.54, y_cm / 2.54
        rad = math.radians(self.rotation)
        cs, sn = math.cos(rad), math.sin(rad)
        # Le Fortran calcule des pixels ecran (y croissant vers le BAS), d'ou
        # son "scrvor - ..." pour compenser. Les coordonnees de donnees
        # matplotlib croissent deja vers le HAUT nativement : pas de flip ici,
        # sous peine de doublement l'inversion (bug trouve via stereo.py, ou
        # N/S se retrouvaient inverses).
        return (
            self.origin_x + (x_in * cs - y_in * sn),
            self.origin_y + (x_in * sn + y_in * cs),
        )

    # ------------------------------------------------------------------
    # plot : deplacement/trace du crayon
    # ------------------------------------------------------------------

    def plot(self, x_cm: float, y_cm: float, mode: int) -> None:
        """Equivalent de `plot(x1,y1,ipen)`.

        mode -4 : revient a l'origine ABSOLUE (0,0) de l'Axes - equivalent
                  de `case(-4): scrhor=90.0; scrvor=600.0` (page origin fixe
                  dans le Fortran ; ici simplement (0,0), l'Axes matplotlib
                  n'ayant pas de page physique). x_cm/y_cm ignores.
        mode -3 : deplace le crayon ET redefinit l'origine sur ce point
                  (les coordonnees suivantes en mode 2/3 seront relatives
                  a ce nouveau point - equivalent de `call plot(x,y,-3)`
                  utilise en debut de routine pour placer un diagramme).
        mode  3 : deplace le crayon, sans tracer (leve de crayon).
        mode  2 : trace une ligne du crayon courant jusqu'a (x,y) (baisse
                  de crayon), avec la couleur/epaisseur courantes.
        autre   : ignore (mode 999 = fin de trace, sans effet ici).
        """
        if mode == -4:
            self.origin_x = self.origin_y = self.x = self.y = 0.0
            return
        tx, ty = self._transform(x_cm, y_cm)
        if mode == -3:
            self.origin_x, self.origin_y = tx, ty
            self.x, self.y = tx, ty
        elif mode == 3:
            self.x, self.y = tx, ty
        elif mode == 2:
            self.ax.add_line(Line2D(
                [self.x, tx], [self.y, ty],
                color=self.pen_color, linewidth=self.line_width,
                solid_capstyle="round",
            ))
            self.x, self.y = tx, ty

    # ------------------------------------------------------------------
    # symbol / label2 : marqueurs
    # ------------------------------------------------------------------

    def symbol(self, x_cm: float, y_cm: float, h_cm: float, ityp: int, nt: int) -> None:
        """Equivalent de `symbol(x,y,h,ityp,nt)`. `nt<=-2` trace d'abord une
        ligne depuis la position courante (relie les points consecutifs
        d'une serie) ; `ityp` selectionne la forme du marqueur (voir
        commentaire de module)."""
        if nt <= -2:
            self.plot(x_cm, y_cm, 2)
        self._label2(x_cm, y_cm, h_cm, ityp)

    def pick_point(self, x_cm: float, y_cm: float, kind: str, data) -> None:
        """Enregistre un point interactif INVISIBLE (scatter alpha=0,
        picker actif) a la meme position qu'un `symbol()`/`plot()` deja
        trace (memes origine/rotation courantes via `_transform`) - pour
        le clic-pour-info cote app.py (`pick_event`), sans jamais changer
        le rendu visible - demande explicite utilisateur ("un collegue
        m'a demande si il etait possible de cliquer sur des donnees d'un
        graphique pour obtenir des informations, par exemple sur un
        zijderveld la temperature, ou sur un stereo results le nom du
        specimen"). N'affecte JAMAIS l'export SVG : svgwriter.py est un
        pipeline de rendu totalement independant (ecrit du SVG brut,
        n'importe pas ce module), see docstring de module ci-dessus.
        `kind` identifie le type de point pour le gestionnaire de clic
        (ex. "zijderveld_step", "stereo_specimen") ; `data` est
        l'information a afficher (objet Measurement, id specimen...),
        recuperee via `event.artist._starmac_pick_data` sans avoir a
        redecoder les coordonnees affichees."""
        tx, ty = self._transform(x_cm, y_cm)
        artist = self.ax.scatter([tx], [ty], s=250, alpha=0, picker=True, zorder=1000)
        artist._starmac_pick_kind = kind
        artist._starmac_pick_data = data

    def _label2(self, x_cm: float, y_cm: float, size_cm: float, itype: int) -> None:
        r = size_cm / 2.0
        if itype in (1, 5):  # plus
            self.plot(x_cm, y_cm - r, 3); self.plot(x_cm, y_cm + r, 2)
            self.plot(x_cm - r, y_cm, 3); self.plot(x_cm + r, y_cm, 2)
        elif itype == 2:  # point (petit disque plein - une ligne de longueur
            # nulle ne produit aucun marqueur visible en matplotlib malgre
            # solid_capstyle="round", bug reel corrige ici - meme fix que
            # AMS_Py/plotlib.py, meme bibliotheque partagee)
            self.circl2(x_cm, y_cm, r * 0.3, 1, 1)
        elif itype == 3:  # barre verticale
            self.plot(x_cm, y_cm - r, 3); self.plot(x_cm, y_cm + r, 2)
        elif itype == 4:  # barre horizontale
            self.plot(x_cm - r, y_cm, 3); self.plot(x_cm + r, y_cm, 2)
        elif itype == 6:  # croix (X)
            self.plot(x_cm - r, y_cm - r, 3); self.plot(x_cm + r, y_cm + r, 2)
            self.plot(x_cm - r, y_cm + r, 3); self.plot(x_cm + r, y_cm - r, 2)
        elif itype == 7:  # asterisque
            self.plot(x_cm, y_cm - r, 3); self.plot(x_cm, y_cm + r, 2)
            self.plot(x_cm - r, y_cm, 3); self.plot(x_cm + r, y_cm, 2)
            self.plot(x_cm - r * .7, y_cm - r * .7, 3); self.plot(x_cm + r * .7, y_cm + r * .7, 2)
            self.plot(x_cm - r * .7, y_cm + r * .7, 3); self.plot(x_cm + r * .7, y_cm - r * .7, 2)
        elif itype == 8:  # cercle
            self.circl2(x_cm, y_cm, r, 1, 0)
        elif itype == 14:  # cercle plein
            self.circl2(x_cm, y_cm, r, 1, 1)
        elif itype == 20:  # cercle blanc
            self.circl2(x_cm, y_cm, r, 1, 1, fill_override="white")
        elif itype in (9, 15):  # carre (plein si 15)
            self._polygon(
                [(x_cm - r, y_cm - r), (x_cm + r, y_cm - r),
                 (x_cm + r, y_cm + r), (x_cm - r, y_cm + r)],
                filled=itype == 15,
            )
        elif itype in (10, 16):  # triangle (plein si 16)
            self._polygon(
                [(x_cm - r, y_cm - r), (x_cm, y_cm + r), (x_cm + r, y_cm - r)],
                filled=itype == 16,
            )
        elif itype in (11, 17):  # peck / triangle inverse (plein si 17)
            self._polygon(
                [(x_cm - r, y_cm + r), (x_cm, y_cm - r), (x_cm + r, y_cm + r)],
                filled=itype == 17,
            )
        elif itype in (12, 18):  # losange (plein si 18)
            self._polygon(
                [(x_cm, y_cm - r), (x_cm + r, y_cm), (x_cm, y_cm + r), (x_cm - r, y_cm)],
                filled=itype == 18,
            )
        elif itype in (13, 19):  # etoile (pleine si 19)
            self._polygon(self._star_points(x_cm, y_cm, size_cm), filled=itype == 19)
        self.plot(x_cm, y_cm, 3)

    def _polygon(self, pts_cm: List[Tuple[float, float]], filled: bool) -> None:
        pts = [self._transform(x, y) for x, y in pts_cm]
        self.ax.add_patch(Polygon(
            pts, closed=True,
            facecolor=self.fill_color if filled else "none",
            edgecolor=self.pen_color, linewidth=self.line_width,
        ))

    @staticmethod
    def _star_points(x: float, y: float, size: float) -> List[Tuple[float, float]]:
        dr = math.pi / 180.0
        r = size / (1.0 + math.cos(36 * dr))
        a = r * math.sin(18 * dr) / math.sin(126 * dr)
        pts = []
        for i in range(1, 6):
            sa = a * math.sin((i * 72 - 36) * dr)
            ca = a * math.cos((i * 72 - 36) * dr)
            sr = r * math.sin(i * 72 * dr)
            cr = r * math.cos(i * 72 * dr)
            pts.append((x + sa, y + ca))
            pts.append((x + sr, y + cr))
        return pts

    def circl2(
        self, xc_cm: float, yc_cm: float, radius_cm: float,
        istroke: int, ifill: int, fill_override: Optional[str] = None,
    ) -> None:
        """Equivalent de `circl2(xc,yc,radius,istroke,ifill)`. ifill=0 ou 2 :
        non rempli ; sinon rempli avec la couleur de remplissage courante
        (ou `fill_override`)."""
        cx, cy = self._transform(xc_cm, yc_cm)
        r_in = radius_cm / 2.54
        face = "none" if ifill in (0, 2) else (fill_override or self.fill_color)
        self.ax.add_patch(Circle(
            (cx, cy), r_in, facecolor=face,
            edgecolor=self.pen_color, linewidth=self.line_width,
        ))

    # ------------------------------------------------------------------
    # Texte : plottxt / number
    # ------------------------------------------------------------------

    def plottxt(
        self, x_cm: float, y_cm: float, height_cm: float, text: str,
        angle: float = 0.0, nchar: Optional[int] = None,
    ) -> None:
        """Equivalent de `plottxt(x,y,height,txtstr,angle,nchar)`.

        Rendu en glyphes vectoriels (TextPath) positionnes/mis a l'echelle
        en coordonnees DATA (comme svgwriter.SVGWriter.plottxt), et non via
        `ax.text(fontsize=...)` (points ecran, fixes) : sinon le texte ne
        suit pas le zoom/redimensionnement de la figure - contrairement aux
        traits et cercles, deja en coordonnees data via `_transform`.

        Pas de decalage `y+1.25*height` : ce decalage n'existe, cote Fortran,
        que pour le rendu ECRAN (AWE_canvasDrawText) - la branche SVG (celle
        dont on part ici) recalcule x2,y2 a partir des coordonnees BRUTES,
        sans ce decalage (verifie au pixel pres contre un export reel)."""
        if nchar is not None:
            text = text[:nchar]
        if not text:
            return
        tx, ty = self._transform(x_cm, y_cm)
        scale = (height_cm / 2.54) / _CAP_HEIGHT
        path = TextPath((0, 0), text, size=1, prop=_FONT)
        transform = Affine2D().scale(scale).rotate_deg(angle).translate(tx, ty) + self.ax.transData
        self.ax.add_patch(PathPatch(
            path, transform=transform, facecolor=self.pen_color,
            edgecolor="none", linewidth=0,
        ))

    def number(
        self, x_cm: float, y_cm: float, height_cm: float, realnb: float,
        angle: float = 0.0, ndec: int = 0,
    ) -> None:
        """Equivalent de `number(x,y,height,realnb,angle,ndec)`."""
        self.plottxt(x_cm, y_cm, height_cm, f"{realnb:.{max(ndec, 0)}f}", angle)

    # ------------------------------------------------------------------
    # No-ops (equivalents purement lies a l'export SVG/PostScript, hors
    # perimetre de ce portage - presents pour que le code transcrit depuis
    # le Fortran n'ait pas besoin d'etre modifie sur ces appels)
    # ------------------------------------------------------------------

    def group(self, start: bool) -> None:
        pass

    def plotnd(self) -> None:
        pass
