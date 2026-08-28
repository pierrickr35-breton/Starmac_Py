from dataclasses import dataclass, field
import re
from typing import List


@dataclass
class PmagData:
  id: str = ""
  cin: safe_float = 0.0
  caz: safe_float = 0.0
  dip: safe_float = 0.0
  str_val: safe_float = 0.0
  norme: str = ""
  vol: safe_float = 0.0
  com: str = ""

  lat: safe_float = 999.0
  rlong: safe_float = 0.0
  altitude: safe_float = 0.0
  year: safe_float = 0.0
  month: safe_float = 0.0
  day: safe_float = 0.0
  hour: safe_float = 0.0
  minute: safe_float = 0.0
  azmag: safe_float = 0.0
  azsun: safe_float = 0.0
  outilorient: str = "A"
  roche: str = ""

  nbmes: int = 0
  etape: List[str] = field(default_factory=list)
  cod1: List[str] = field(default_factory=list)
  cod2: List[str] = field(default_factory=list)
  x: List[float] = field(default_factory=list)
  y: List[float] = field(default_factory=list)
  z: List[float] = field(default_factory=list)
  q: List[float] = field(default_factory=list)
  ins: List[str] = field(default_factory=list)
  s: List[float] = field(default_factory=list)
  xech: List[float] = field(default_factory=list)
  yech: List[float] = field(default_factory=list)
  zech: List[float] = field(default_factory=list)
  heuremes: List[str] = field(default_factory=list)


def import_texte(filename: str) -> List[PmagData]:
  pmag_list: List[PmagData] = []
  current_pmag: PmagData = None

  with open(filename, "r", encoding="latin-1") as f:
    lines = f.readlines()

  idx = 0
  while idx < len(lines):
    line = lines[idx].rstrip("\r\n")

    # Ignorer les lignes vides
    if not line.strip():
      idx += 1
      continue

    # Détection du début d'un échantillon (Identique à `truc=='I'` en Fortran)
    if line.strip().startswith("Id:"):
      current_pmag = PmagData()

      # Extraction des paramètres d'en-tête (Id, in, az, dip, str, v/m, com)
      for match in re.finditer(r"(Id|in|az|dip|str|v|m|com):", line):
        key = match.group(1)
        start = match.end()
        val_str = line[start : start + 12].strip()

        if key == "Id":
          current_pmag.id = val_str
        elif key == "in":
          current_pmag.cin = float(val_str) if val_str else 0.0
        elif key == "az":
          current_pmag.caz = float(val_str) if val_str else 0.0
        elif key == "dip":
          current_pmag.dip = float(val_str) if val_str else 0.0
        elif key == "str":
          current_pmag.str_val = float(val_str) if val_str else 0.0
        elif key in ("v", "m"):
          current_pmag.norme = key
          current_pmag.vol = float(val_str) if val_str else 0.0
        elif key == "com":
          current_pmag.com = line[start : start + 15].strip()

      pmag_list.append(current_pmag)
      idx += 1

      # Traitement de la ligne de coordonnées/orientation (L: G: H: T: ...)
      if idx < len(lines) and lines[idx].strip().startswith("L:"):
        line_coords = lines[idx].strip()

        # Extraction par clés regex
        lat_m = re.search(r"L:\s*([\d.-]+)", line_coords)
        lon_m = re.search(r"G:\s*([\d.-]+)", line_coords)
        alt_m = re.search(r"H:\s*([\d.-]+)", line_coords)
        azm_m = re.search(r"azm:\s*(\S+)", line_coords)
        azs_m = re.search(r"azs:\s*(\S+)", line_coords)
        or_m = re.search(r"Or:\s*(\S+)", line_coords)

        if lat_m:
          current_pmag.lat = float(lat_m.group(1))
        if lon_m:
          current_pmag.rlong = float(lon_m.group(1))
        if alt_m:
          current_pmag.altitude = float(alt_m.group(1))
        if azm_m and azm_m.group(1) != "n.d":
          current_pmag.azmag = float(azm_m.group(1))
        if azs_m and azs_m.group(1) != "n.d":
          current_pmag.azsun = float(azs_m.group(1))
        if or_m:
          current_pmag.outilorient = or_m.group(1)

        idx += 1

      # Ligne éventuelle avec métadonnées Site/Sample/Roche
      if idx < len(lines) and lines[idx].strip().startswith("Site:"):
        current_pmag.roche = lines[idx].strip()
        idx += 1

      continue

    # Lecture des blocs de mesures
    if current_pmag is not None:
      parts = line.split()
      if len(parts) >= 7:
        # Équivalent de read(chaine, 201 ou 2011) en Fortran
        step_code = parts[0]
        current_pmag.etape.append(step_code[:4])
        current_pmag.cod1.append(step_code[4:5] if len(step_code) > 4 else "")
        current_pmag.cod2.append(step_code[5:6] if len(step_code) > 5 else "")

        current_pmag.x.append(float(parts[1]))
        current_pmag.y.append(float(parts[2]))
        current_pmag.z.append(float(parts[3]))
        current_pmag.q.append(float(parts[4]))
        current_pmag.ins.append(parts[5])
        current_pmag.s.append(float(parts[6]))

        # Format long (xech, yech, zech, heuremes)
        if len(parts) >= 11:
          current_pmag.xech.append(float(parts[7]))
          current_pmag.yech.append(float(parts[8]))
          current_pmag.zech.append(float(parts[9]))
          current_pmag.heuremes.append(" ".join(parts[10:]))
        else:
          current_pmag.xech.append(0.0)
          current_pmag.yech.append(0.0)
          current_pmag.zech.append(0.0)
          current_pmag.heuremes.append("2000/01/01 12:00:00")

        current_pmag.nbmes += 1

    idx += 1

  return pmag_list