"""The Python port must reproduce the original Fortran output (.plt files)."""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import schemes  # noqa: E402

TOLERANCE = 1e-5


def max_difference(scheme: str, folder: str) -> float:
    snapshots = schemes.run(scheme)
    worst = 0.0
    for k in range(schemes.OUTPUTS):
        fortran = schemes.read_plt(str(ROOT / "results" / folder / f"{folder}_{k + 1}.plt"))
        worst = max(worst, float(np.abs(snapshots[k] - fortran).max()))
    return worst


def test_upwind_matches_fortran_output():
    assert max_difference("upwind", "upwind") < TOLERANCE


def test_committed_maccormack_output_is_the_backward_variant():
    assert max_difference("maccormack-run", "maccormack") < TOLERANCE


def test_textbook_maccormack_is_more_accurate_than_upwind():
    for row in schemes.error_table()[:3]:
        assert row["maccormack_rms"] < row["upwind_rms"]


def test_courant_number_is_stable():
    assert 0 < schemes.COURANT <= 1


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok -", name)
    print("max |python - fortran|: upwind", max_difference("upwind", "upwind"),
          "| maccormack-run", max_difference("maccormack-run", "maccormack"))
