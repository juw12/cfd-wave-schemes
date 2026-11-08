"""1-D linear wave (advection) equation u_t + c u_x = 0 solved two ways.

Python port of the two Fortran studies in this repo (fortran/upwind.f90 and
fortran/maccormack.f90). Same grid, same initial pulse, same time stepping and
the same output schedule, in float32 like the Fortran `real` type.

Three schemes:
  upwind           first-order upwind (matches results/upwind/*.plt exactly)
  maccormack       textbook MacCormack: forward predictor, backward corrector
                   (what fortran/maccormack.f90 says)
  maccormack-run   backward predictor + backward corrector. This is what the
                   committed results/maccormack/*.plt were actually produced
                   with (reproduced to 1e-7); it is only first-order accurate.

Run:  python schemes.py            -> prints an error table, writes comparison.png
"""
from __future__ import annotations

import numpy as np

DX = 0.01            # grid spacing
NODES = 301          # grid points, domain x = 0 .. 3.0
DT = 0.01            # time step
SPEED = 0.09         # wave speed c
COURANT = SPEED * DT / DX   # C = c dt / dx = 0.09
PULSE_NODES = 11     # initial half-sine pulse on x = 0 .. 0.1
OUTPUTS = 10         # the Fortran writes 10 snapshots


def initial_pulse() -> np.ndarray:
    u = np.zeros(NODES, dtype=np.float32)
    i = np.arange(PULSE_NODES, dtype=np.float32)
    u[:PULSE_NODES] = np.sin(np.float32(10 * 3.14159265359) * i * np.float32(DX))
    return u


def upwind_step(u: np.ndarray) -> np.ndarray:
    """First-order upwind: u_i <- u_i - C (u_i - u_{i-1}), interior nodes only."""
    c = np.float32(COURANT)
    new = u.copy()
    new[1:-1] = u[1:-1] - c * (u[1:-1] - u[:-2])
    return new


def maccormack_step(u: np.ndarray, previous: np.ndarray) -> np.ndarray:
    """Explicit MacCormack predictor (forward) + corrector (backward).

    `previous` is the last corrector result; the Fortran keeps node 1 from it
    because the corrector loop starts at node 2.
    """
    c = np.float32(COURANT)
    predicted = np.zeros_like(u)
    predicted[:-1] = u[:-1] - c * (u[1:] - u[:-1])
    new = previous.copy()
    new[1:] = np.float32(0.5) * ((u[1:] + predicted[1:]) - c * (predicted[1:] - predicted[:-1]))
    return new


def maccormack_run_step(u: np.ndarray, previous: np.ndarray) -> np.ndarray:
    """The variant behind the committed MacCormack .plt files: the predictor
    uses a backward difference on interior nodes (i = 2 .. xmax-1, up(xmax) = 0),
    so both half-steps lean upwind. Reproduces all ten .plt snapshots to 1e-7."""
    c = np.float32(COURANT)
    predicted = np.zeros_like(u)
    predicted[1:-1] = u[1:-1] - c * (u[1:-1] - u[:-2])
    new = previous.copy()
    new[1:] = np.float32(0.5) * ((u[1:] + predicted[1:]) - c * (predicted[1:] - predicted[:-1]))
    return new


STEPPERS = {
    "upwind": lambda u, previous: upwind_step(u),
    "maccormack": maccormack_step,
    "maccormack-run": maccormack_run_step,
}


def snapshot_steps() -> list[int]:
    """Cumulative step counts at each output.

    The Fortran runs 100*j more steps before writing snapshot j, so the
    totals are 100, 300, 600, ... 5500 steps (t = 1, 3, 6, ... 55 s).
    """
    totals = []
    steps = 0
    for j in range(1, OUTPUTS + 1):
        steps += 100 * j
        totals.append(steps)
    return totals


def run(scheme: str) -> list[np.ndarray]:
    step = STEPPERS[scheme]
    u = initial_pulse()
    after = np.zeros(NODES, dtype=np.float32)
    snapshots = []
    for j in range(1, OUTPUTS + 1):
        for _ in range(100 * j):
            after = step(u, after)
            u = after
        snapshots.append(after.copy())
    return snapshots


def exact(steps: int) -> np.ndarray:
    """The initial pulse moved right by c t (it leaves the 3.0 domain after ~32 s)."""
    x = np.arange(NODES) * DX - SPEED * steps * DT
    u = np.where((x >= 0) & (x <= (PULSE_NODES - 1) * DX), np.sin(10 * np.pi * x), 0.0)
    return u


def read_plt(path: str) -> np.ndarray:
    data = np.loadtxt(path)
    return data[:, 1].astype(np.float32)


def error_table() -> list[dict]:
    rows = []
    results = {name: run(name) for name in STEPPERS}
    for k, steps in enumerate(snapshot_steps()):
        reference = exact(steps)
        row = {"snapshot": k + 1, "t_s": round(steps * DT, 2)}
        for name, snapshots in results.items():
            row[f"{name}_peak"] = float(snapshots[k].max())
            row[f"{name}_rms"] = float(np.sqrt(np.mean((snapshots[k] - reference) ** 2)))
        rows.append(row)
    return rows


def plot(path: str = "comparison.png") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results = {name: run(name) for name in STEPPERS}
    styles = {
        "upwind": ("#0891b2", "upwind (1st order)"),
        "maccormack-run": ("#d97706", "MacCormack as run (backward predictor)"),
        "maccormack": ("#db2777", "MacCormack textbook (2nd order)"),
    }
    x = np.arange(NODES) * DX
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for ax, k in zip(axes, (0, 2, 4)):
        steps = snapshot_steps()[k]
        ax.plot(x, exact(steps), color="#94a3b8", lw=2, label="exact")
        for name, (color, label) in styles.items():
            ax.plot(x, results[name][k], color=color, label=label)
        ax.set_ylabel("u")
        ax.set_title(f"t = {steps * DT:g} s", loc="left", fontsize=10)
        ax.set_xlim(0, 1.6)
        ax.grid(alpha=0.3)
    axes[0].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("x")
    fig.suptitle(f"1-D wave equation, C = {COURANT:.2f}")
    fig.tight_layout()
    fig.savefig(path, dpi=120)


if __name__ == "__main__":
    print(f"Courant number C = {COURANT:.2f}")
    print("RMS error against the exact solution (peak in brackets)")
    print(f"{'#':>2} {'t (s)':>6} {'upwind':>16} {'MacC as run':>16} {'MacC textbook':>16}")
    for row in error_table():
        cells = [f"{row[f'{n}_rms']:.4f} ({row[f'{n}_peak']:.2f})" for n in ("upwind", "maccormack-run", "maccormack")]
        print(f"{row['snapshot']:>2} {row['t_s']:>6g} " + " ".join(f"{c:>16}" for c in cells))
    plot()
    print("wrote comparison.png")
