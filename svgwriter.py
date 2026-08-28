"""
Port fidele des routines d'ecriture SVG de graphicsAWE.f95 : `svginit`,
`svgplot`, `chckpgsvg`, et les branches SVG de `plottxt`/`circl2`/`newpen`/
`thickn`. Verifie contre des exports SVG reels de Starmac_OSX (memes
balises, memes formats numeriques, meme encodage iso-8859-1).

Point cle repris du Fortran (rate dans une premiere tentative qui ecrivait
un <line> par segment) : `svgplot` ACCUMULE les points consecutifs d'un
trace continu (mode 2 apres mode 2) dans UNE seule balise <polyline
points="x1,y1 x2,y2 ...">, et ne la ferme que lorsque le crayon est leve
(mode 3) ou en fin de trace (mode 999) - via `chckpgsvg`. C'est le meme
`SVGWriter` qui joue le role de `PlotContext` (memes methodes
plot/symbol/circl2/plottxt/newpen/thickn/clear), utilisable directement par
stereo.py/zijderveld.py a la place d'un PlotContext matplotlib.

Perimetre : les symboles bases sur `plot`/`circl2` (croix, plus, cercle
plein/ouvert - types 5,8,13,14,20, ceux utilises par Zijderveld/Stereo)
sont fideles. Les symboles polygonaux (carre/triangle/losange/etoile,
types 9-19, non utilises par les graphiques actuellement portes) passent
par un `<polygon>` generique plutot qu'une transcription de
polysquare/polytriangle/polystar (non lues en detail - a faire si un
graphique futur en a besoin).
"""

import math
from typing import List, Optional, Tuple

# newpen() : palette 1-9 (memes couleurs que graphicsAWE.f95, en hexa SANS '#')
_PEN_COLORS = {
    1: "000000", 2: "FFFFFF", 3: "FF0000", 4: "008000", 5: "0000FF",
    6: "00FFFF", 7: "EE82EE", 8: "FFFF00", 9: "FF1493",
}


class SVGWriter:
    """Equivalent SVG-only de PlotContext (memes methodes publiques), portant
    directement `svginit`/`svgplot`/`chckpgsvg`/les branches SVG de
    `circl2`/`plottxt`. Coordonnees en cm, comme le Fortran d'origine."""

    SCRDPI = 72.0  # points/pouce (constante Fortran, common /calscr/)

    def __init__(self, width_cm: float = 20.0, height_cm: float = 20.0):
        self.width_cm = width_cm
        self.height_cm = height_cm
        self.origin_x = 0.0  # pixels (equivalent scrhor)
        self.origin_y = 0.0  # pixels (equivalent scrvor)
        self._page_origin_x = 0.0  # equivalent de la constante 90.0 (case -4)
        self._page_origin_y = 0.0  # equivalent de la constante 600.0 (case -4)
        self.x = 0.0
        self.y = 0.0
        self.rotation = 0.0
        self.pen_color = "000000"    # hex SANS '#' (strokecol Fortran)
        self.fill_color = "000000"   # polyfil Fortran
        self.line_width = 2.0
        self.ifilpoly = 0

        self._lines: List[str] = []
        self._lninpg = 0  # 0 = pas de polyline ouverte, 2 = polyline en cours
        self._pending_points: List[Tuple[float, float]] = []
        self._pending_attrs: Tuple[str, str, float] = ("none", "000000", 1.0)
        self.xltt = 0.0
        self.yltt = 0.0

    # ------------------------------------------------------------------
    # Etat du crayon : newpen / newpencol / thickn / clear
    # ------------------------------------------------------------------

    def newpen(self, ipen: int) -> None:
        col = _PEN_COLORS.get(ipen, "000000")
        self.pen_color = col
        self.fill_color = col

    def newpencol(self, sr: int, sg: int, sb: int, fr: int, fg: int, fb: int) -> None:
        def hexcol(r: int, g: int, b: int) -> str:
            if r < 0 or g < 0 or b < 0 or r > 255 or g > 255 or b > 255:
                return "000000"
            return f"{r:02X}{g:02X}{b:02X}"

        self.pen_color = hexcol(sr, sg, sb)
        self.fill_color = hexcol(fr, fg, fb)

    def thickn(self, sizel: float) -> None:
        self.line_width = sizel * 2

    def clear(self) -> None:
        """Equivalent de `clear()` : reinitialise le document et l'origine."""
        self._lines = []
        self._lninpg = 0
        self._pending_points = []
        self.x = self.y = self.origin_x = self.origin_y = 0.0
        self.xltt = self.yltt = 0.0
        self._page_origin_x = self._page_origin_y = 0.0

    def set_origin_px(self, x_px: float, y_px: float) -> None:
        """Equivalent de l'initialisation scrhor/scrvor dans `plots()` : fixe
        l'origine ABSOLUE de la page (memorisee pour `plot(...,-4)`, qui y
        revient - equivalent de `case(-4): scrhor=90.0; scrvor=600.0`)."""
        self.origin_x = self.x = self.xltt = x_px
        self.origin_y = self.y = self.yltt = y_px
        self._page_origin_x = x_px
        self._page_origin_y = y_px

    # ------------------------------------------------------------------
    # Transformation de coordonnees (identique a plot()/circl2()/plottxt())
    # ------------------------------------------------------------------

    def _transform(self, x_cm: float, y_cm: float) -> Tuple[float, float]:
        x_in, y_in = x_cm / 2.54, y_cm / 2.54
        rad = math.radians(self.rotation)
        cs, sn = math.cos(rad), math.sin(rad)
        return (
            self.origin_x + (x_in * cs - y_in * sn) * self.SCRDPI,
            self.origin_y - (x_in * sn + y_in * cs) * self.SCRDPI,
        )

    # ------------------------------------------------------------------
    # plot : port de `svgplot` (accumulation de <polyline>)
    # ------------------------------------------------------------------

    def plot(self, x_cm: float, y_cm: float, mode: int) -> None:
        """Equivalent de `plot(x1,y1,ipen)` / `svgplot`.

        mode -4 : revient a l'origine ABSOLUE de la page (equivalent de
                  `case(-4): scrhor=90.0; scrvor=600.0` - constantes fixees
                  ici par le dernier `set_origin_px()`/`clear()`), sans
                  tracer. x_cm/y_cm ignores (le Fortran les ignore aussi).
        mode -3 : redefinit l'origine sur ce point (sans tracer) - relatif a
                  l'origine COURANTE (compose avec elle, ne l'ecrase pas).
        mode  3 : leve le crayon (ferme la polyline en cours, s'il y en a une).
        mode  2 : baisse le crayon - accumule le point dans la polyline
                  ouverte (ou en ouvre une nouvelle si aucune n'est active).
        """
        if mode == -4:
            self.origin_x, self.origin_y = self._page_origin_x, self._page_origin_y
            self.x, self.y = self.origin_x, self.origin_y
            self.xltt, self.yltt = self.x, self.y
            return
        tx, ty = self._transform(x_cm, y_cm)
        if mode == -3:
            self.origin_x, self.origin_y = tx, ty
            self.x, self.y = tx, ty
        elif mode == 3:
            self._close_polyline()
            self.x, self.y = tx, ty
        elif mode == 2:
            if self._lninpg == 0:
                fill = f"#{self.fill_color}" if self.ifilpoly == 1 else "none"
                self._pending_attrs = (fill, self.pen_color, max(self.line_width, 0.5))
                self._pending_points = [(self.xltt, self.yltt)]
                self._lninpg = 2
            self._pending_points.append((tx, ty))
            self.x, self.y = tx, ty
        self.xltt, self.yltt = self.x, self.y

    def _close_polyline(self) -> None:
        if self._lninpg == 2:
            fill, stroke, width = self._pending_attrs
            pts = " ".join(f"{x:.3f},{y:.3f}" for x, y in self._pending_points)
            self._lines.append(
                f'<polyline fill="{fill}" stroke="#{stroke}" '
                f'stroke-width="{width:.3f}" points="{pts}"/>'
            )
            self._lninpg = 0
            self._pending_points = []

    # ------------------------------------------------------------------
    # symbol / label2 : marqueurs (memes decompositions que PlotContext)
    # ------------------------------------------------------------------

    def symbol(self, x_cm: float, y_cm: float, h_cm: float, ityp: int, nt: int) -> None:
        """Equivalent de `symbol(x,y,h,ityp,nt)`."""
        if nt <= -2:
            self.plot(x_cm, y_cm, 2)
        self._label2(x_cm, y_cm, h_cm, ityp)

    def _label2(self, x_cm: float, y_cm: float, size_cm: float, itype: int) -> None:
        r = size_cm / 2.0
        if itype in (1, 5):  # plus
            self.plot(x_cm, y_cm - r, 3); self.plot(x_cm, y_cm + r, 2)
            self.plot(x_cm - r, y_cm, 3); self.plot(x_cm + r, y_cm, 2)
        elif itype == 2:  # point
            self.plot(x_cm, y_cm, 3); self.plot(x_cm, y_cm, 2)
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
            self.circl2(x_cm, y_cm, r, 1, 1, fill_override="FFFFFF")
        elif itype in (9, 15):  # carre (plein si 15) - <polygon> generique
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
        fill = f"#{self.fill_color}" if filled else "none"
        pts_str = " ".join(f"{x:.3f},{y:.3f}" for x, y in pts)
        self._lines.append(
            f'<polygon fill="{fill}" stroke="#{self.pen_color}" '
            f'stroke-width="{max(self.line_width, 0.5):.3f}" points="{pts_str}"/>'
        )

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
        """Equivalent de la branche SVG de `circl2` : balise <circle> directe
        (pas d'accumulation - contrairement a plot(), un seul point)."""
        self._close_polyline()
        cx, cy = self._transform(xc_cm, yc_cm)
        r_px = (radius_cm / 2.54) * self.SCRDPI
        fill_hex = fill_override or self.fill_color
        fill = f"#{fill_hex}" if ifill not in (0, 2) else "none"
        self._lines.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r_px:.2f}" '
            f'fill="{fill}" stroke="#{self.pen_color}" '
            f'stroke-width="{max(self.line_width, 0.5):.2f}"/>'
        )

    # ------------------------------------------------------------------
    # Texte : port de la branche SVG de `plottxt` / `number`
    # ------------------------------------------------------------------

    def plottxt(
        self, x_cm: float, y_cm: float, height_cm: float, text: str,
        angle: float = 0.0, nchar: Optional[int] = None,
    ) -> None:
        """Equivalent de la branche SVG de `plottxt`. Le Fortran d'origine
        ne gere QUE angle=0 ou angle=90 (matrices ecrites en dur) - un
        autre angle n'ecrit rien, comportement reproduit ici.

        Le decalage `y+1.25*height` du Fortran n'existe QUE pour le rendu
        ECRAN (AWE_canvasDrawText) - la branche SVG recalcule x2,y2 a partir
        des coordonnees BRUTES (x/2.54, y/2.54, sans ce decalage). Verifie
        au pixel pres contre zijder-11CL7801A.svg (bloc scale/id/AF/SC)."""
        self._close_polyline()
        if nchar is not None:
            text = text[:nchar]
        tx, ty = self._transform(x_cm, y_cm)
        size_px = height_cm * self.SCRDPI / 2.54 * 5.0 / 3.0
        esc = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

        if angle == 0:
            matrix = f"1 0 0 1 {tx:.2f} {ty:.2f}"
        elif angle == 90:
            matrix = f"0 -1 1 0 {tx:.2f} {ty:.2f}"
        else:
            return

        self._lines.append(
            f'<text transform="matrix({matrix})" fill="#{self.pen_color}" '
            f'font-family="\'Helvetica\'" font-size="{size_px:.2f}">{esc}  </text>'
        )

    def number(
        self, x_cm: float, y_cm: float, height_cm: float, realnb: float,
        angle: float = 0.0, ndec: int = 0,
    ) -> None:
        """Equivalent de `number(x,y,height,realnb,angle,ndec)`."""
        self.plottxt(x_cm, y_cm, height_cm, f"{realnb:.{max(ndec, 0)}f}", angle)

    # ------------------------------------------------------------------
    # No-ops (equivalents purement lies au mode ecran, sans effet ici)
    # ------------------------------------------------------------------

    def group(self, start: bool) -> None:
        pass

    def plotnd(self) -> None:
        """Equivalent partiel du mode 999 (fin de trace) : ferme toute
        polyline en cours. L'ecriture des balises de fin de document
        (</g></svg>) est geree par `to_string()`/`save()`."""
        self._close_polyline()

    # ------------------------------------------------------------------
    # Document : port de `svginit` (en-tete) + fin de document
    # ------------------------------------------------------------------

    def to_string(self) -> str:
        self._close_polyline()
        header = (
            '<?xml version="1.0" encoding="iso-8859-1"?>\n'
            '<svg version="1.2" baseProfile="tiny" xmlns="http://www.w3.org/2000/svg"\n'
            'xmlns:xlink="http://www.w3.org/1999/xlink"\n'
            f' x="0px" y="0px" width="{self.width_cm:.1f}cm" height="{self.height_cm:.1f}cm" xml:space="preserve">\n'
            '<g id="layerpmag0">\n'
        )
        body = "\n".join(self._lines)
        footer = "\n</g>\n</svg>\n"
        return header + body + footer

    def save(self, path: str) -> None:
        with open(path, "w", encoding="iso-8859-1", errors="xmlcharrefreplace") as f:
            f.write(self.to_string())
