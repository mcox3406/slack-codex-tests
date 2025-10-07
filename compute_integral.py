"""Evaluate Fresnel-type integrals ∫₀^∞ sin(a x²) dx for several scale factors."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Sequence

import importlib.util

if importlib.util.find_spec("scipy") is None:
    raise ModuleNotFoundError(
        "SciPy is required for this script. Install it with `pip install scipy`."
    )

from scipy import integrate
import matplotlib.pyplot as plt


@dataclass
class Result:
    a: float
    analytic: float
    scipy_value: float
    scipy_error: float
    adaptive_value: float
    tail_estimate: float


def analytic_value(a: float) -> float:
    """Closed form √(π/(8a)) for a > 0."""

    return math.sqrt(math.pi / (8.0 * a))


def simpson(f: Callable[[float], float], a: float, b: float) -> float:
    c = (a + b) / 2.0
    return (b - a) / 6.0 * (f(a) + 4.0 * f(c) + f(b))


def adaptive_simpson(
    f: Callable[[float], float],
    a: float,
    b: float,
    eps: float,
    max_depth: int = 15,
    whole: float | None = None,
    depth: int = 0,
) -> float:
    if whole is None:
        whole = simpson(f, a, b)

    c = (a + b) / 2.0
    left = simpson(f, a, c)
    right = simpson(f, c, b)
    delta = left + right - whole

    if depth >= max_depth or abs(delta) <= 15.0 * eps:
        return left + right + delta / 15.0

    return (
        adaptive_simpson(f, a, c, eps / 2.0, max_depth, left, depth + 1)
        + adaptive_simpson(f, c, b, eps / 2.0, max_depth, right, depth + 1)
    )


def sinx2_tail(cutoff: float) -> float:
    """Asymptotic tail ∫₍cutoff₎^∞ sin(x²) dx using 7-term expansion."""

    r2 = cutoff * cutoff
    return (
        math.cos(r2) / (2.0 * cutoff)
        - math.sin(r2) / (4.0 * cutoff**2)
        - 3.0 * math.cos(r2) / (8.0 * cutoff**3)
        + 15.0 * math.sin(r2) / (16.0 * cutoff**4)
        + 105.0 * math.cos(r2) / (32.0 * cutoff**5)
        - 945.0 * math.sin(r2) / (64.0 * cutoff**6)
        - 10395.0 * math.cos(r2) / (128.0 * cutoff**7)
    )


def adaptive_fresnel(a: float, cutoff: float = 120.0, tol: float = 1e-9) -> tuple[float, float]:
    """Adaptive Simpson evaluation with asymptotic tail correction."""

    integrand = lambda x: math.sin(a * x * x)
    base = adaptive_simpson(integrand, 0.0, cutoff, tol)
    scaled_cutoff = math.sqrt(a) * cutoff
    tail = sinx2_tail(scaled_cutoff) / math.sqrt(a)
    return base + tail, abs(tail)


def scipy_fresnel(a: float) -> tuple[float, float]:
    """SciPy quad integration using oscillatory weighting."""

    def base(x: float) -> float:
        if x == 0.0:
            return 0.0
        return 0.5 / math.sqrt(x)

    value, error = integrate.quad(
        base,
        0.0,
        math.inf,
        weight="sin",
        wvar=a,
        limit=500,
    )
    return value, error


def evaluate(a_values: Sequence[float]) -> List[Result]:
    results: List[Result] = []
    for a in a_values:
        analytic = analytic_value(a)
        scipy_value, scipy_error = scipy_fresnel(a)
        adaptive_value, tail_estimate = adaptive_fresnel(a)
        results.append(
            Result(
                a=a,
                analytic=analytic,
                scipy_value=scipy_value,
                scipy_error=scipy_error,
                adaptive_value=adaptive_value,
                tail_estimate=tail_estimate,
            )
        )
    return results


def render_plot(results: Sequence[Result], output_path: Path) -> None:
    a_values = [r.a for r in results]
    analytic = [r.analytic for r in results]
    scipy_values = [r.scipy_value for r in results]
    adaptive_values = [r.adaptive_value for r in results]

    plt.figure(figsize=(6.0, 4.0))
    plt.plot(a_values, analytic, marker="o", label="Analytic √(π/(8a))")
    plt.plot(a_values, scipy_values, marker="s", label="SciPy quad")
    plt.plot(a_values, adaptive_values, marker="^", label="Adaptive Simpson")
    plt.xlabel("a")
    plt.ylabel(r"$I(a) = \int_0^\infty \sin(a x^2)\,dx$")
    plt.title("Fresnel Integral Comparison")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close()


def print_results(results: Sequence[Result]) -> None:
    header = f"{'a':>6}  {'Analytic':>12}  {'SciPy':>12}  {'|Δ| (SciPy)':>12}  {'Adaptive':>12}  {'|Δ| (Adaptive)':>15}"
    print(header)
    print("-" * len(header))
    for r in results:
        scipy_diff = abs(r.analytic - r.scipy_value)
        adaptive_diff = abs(r.analytic - r.adaptive_value)
        print(
            f"{r.a:6.2f}  {r.analytic:12.9f}  {r.scipy_value:12.9f}  {scipy_diff:12.3e}  "
            f"{r.adaptive_value:12.9f}  {adaptive_diff:15.3e}"
        )


def main() -> None:
    a_values = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
    results = evaluate(a_values)
    print_results(results)

    output_path = Path(__file__).resolve().parent / "assets" / "fresnel_comparison.png"
    render_plot(results, output_path)
    print(f"\nSaved comparison plot to {output_path.relative_to(Path.cwd())}")


if __name__ == "__main__":
    main()
