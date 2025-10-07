"""Compute the Fresnel integral ∫₀^∞ sin(x²) dx numerically.

The integral can be mapped to the Euler gamma function by substituting
``y = x²``:

    ∫₀^∞ sin(x²) dx = ½ ∫₀^∞ y^{-1/2} sin(y) dy.

For ``0 < s < 1`` the identity ``∫₀^∞ y^{s-1} sin(y) dy = Γ(s) sin(π s / 2)``
holds, and setting ``s = 1/2`` yields the closed-form value ``√(π/8)``.  The
script evaluates this expression numerically and also confirms it by directly
numerically integrating ``sin(x²)`` over a large finite interval with an
adaptive Simpson rule plus a vanishing tail estimate.
"""

from __future__ import annotations

import math
from typing import Callable, Tuple


def analytic_value() -> float:
    """Return √(π/8) via the gamma-function identity."""

    return 0.5 * math.gamma(0.5) * math.sin(math.pi / 4.0)


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


def numerical_check(cutoff: float = 120.0, tol: float = 1e-9) -> Tuple[float, float]:
    """Integrate numerically on ``[0, cutoff]`` and estimate the tail."""

    f = lambda x: math.sin(x * x)
    base = adaptive_simpson(f, 0.0, cutoff, tol)

    # Tail estimate from Abramowitz & Stegun 7.3.26 with the first seven terms
    # of the asymptotic expansion.
    r2 = cutoff * cutoff
    tail = (
        math.cos(r2) / (2.0 * cutoff)
        - math.sin(r2) / (4.0 * cutoff**2)
        - 3.0 * math.cos(r2) / (8.0 * cutoff**3)
        + 15.0 * math.sin(r2) / (16.0 * cutoff**4)
        + 105.0 * math.cos(r2) / (32.0 * cutoff**5)
        - 945.0 * math.sin(r2) / (64.0 * cutoff**6)
        - 10395.0 * math.cos(r2) / (128.0 * cutoff**7)
    )
    return base + tail, abs(tail)


def compute_integral() -> Tuple[float, float, float]:
    closed_form = analytic_value()
    numeric, tail_bound = numerical_check()
    return closed_form, numeric, tail_bound


def main() -> None:
    closed_form, numeric, tail_bound = compute_integral()
    print(f"Closed-form evaluation √(π/8)     = {closed_form:.12f}")
    print(f"Adaptive Simpson + tail estimate = {numeric:.12f}")
    print(f"Tail magnitude upper bound        ≈ {tail_bound:.1e}")
    print(f"Absolute difference                = {abs(closed_form - numeric):.3e}")


if __name__ == "__main__":
    main()
