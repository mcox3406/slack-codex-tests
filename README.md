# Fresnel Integral Sweep

Run `python compute_integral.py` to tabulate $I(a)=\int_0^\infty \sin(a x^2)\,dx$ for several scales, comparing the closed form with SciPy's oscillatory quadrature and an adaptive Simpson cross-check.
The script saves the comparison plot locally as `assets/fresnel_comparison.png`.
