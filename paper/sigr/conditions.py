"""The four test conditions: cell, temperature, SoC window and the
per-condition settings. Everything else is shared (parameters.py)."""
from dataclasses import dataclass

# Basis for C-rates and EFC [A.h]. The tests define both on the cells' 2.5 Ah
# nominal capacity (C/3 = 0.833 A, EFC = throughput / (2 x 2.5 Ah)). The paper
# runs use the model's own capacity, 2.5947 Ah (C/3 = 0.865 A, ~3.8% higher),
# as in the published fit; set both to 2.5 for the test-protocol basis.
MODEL_CAPACITY_AH = 2.5947
RATE_BASIS_AH = MODEL_CAPACITY_AH
EFC_BASIS_AH = MODEL_CAPACITY_AH

UPPER_CUTOFF_V = 4.195
LOWER_CUTOFF_V = 2.5
RPT_RATE = "C/20"


@dataclass(frozen=True)
class Condition:
    key: str  # file-name stem, e.g. "25deg_full_soc"
    cell: str  # cell ID, also the data-file prefix
    data_folder: str  # folder under the data directory
    temperature_K: float
    window: str  # label: "full SoC", "15-95% SoC", "20-80% SoC"
    F0: float  # pore-buffering transmitted-fraction plateau
    width: float  # pore-buffering transition width (the one fitted per-condition value)
    knee_efc_raw: float  # the cell's capacity-knee EFC (raw test EFC, before the RPT1 shift)

    @property
    def label(self):
        return f"{self.cell} ({self.temperature_K - 273.15:.0f} degC, {self.window})"

    @property
    def partial(self):
        return self.window != "full SoC"


CONDITIONS = {
    c.key: c
    for c in [
        Condition("25deg_full_soc", "CELL064", "CELL064_25C_full", 298.15, "full SoC", 0.7, 5e-3, 152.0),
        Condition("45deg_full_soc", "CELL017", "CELL017_45C", 318.15, "full SoC", 0.8, 3e-4, 318.853946),
        Condition("25deg_15-95_soc", "CELL009", "CELL009_25C_15-95", 298.15, "15-95% SoC", 0.7, 1e-3, 298.1722359745961),
        Condition("25deg_20-80_soc", "CELL026", "CELL026_25C_20-80", 298.15, "20-80% SoC", 0.7, 3e-4, 301.41566051398945),
    ]
}


def ageing_cycle(cond, rate_basis_ah=RATE_BASIS_AH):
    """One ageing cycle as pybamm experiment steps (C/3 = rate_basis_ah / 3).
    15-95%: voltage-limited both ends. 20-80%: CC-CV to 3.969 V on top and a
    fixed 1.500 Ah discharge (never updated; 2.5 V safety floor)."""
    if cond.key == "25deg_15-95_soc":
        return ("Discharge at C/3 until 3.14 V", "Charge at C/3 until 4.15 V", "Hold at 4.15 V until 14 mA")
    if cond.key == "25deg_20-80_soc":
        seconds = 1.5 / (rate_basis_ah / 3) * 3600
        return (f"Discharge at C/3 for {seconds:.1f} seconds or until 2.5 V",
                "Charge at C/3 until 3.969 V", "Hold at 3.969 V until 14 mA")
    return (f"Discharge at C/3 until {LOWER_CUTOFF_V} V", f"Charge at C/3 until {UPPER_CUTOFF_V} V",
            f"Hold at {UPPER_CUTOFF_V} V until C/50")
