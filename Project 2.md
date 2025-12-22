# Tunneling in Exotic Quantum Potentials

## 1. Motivation

Most introductory quantum mechanics courses focus on simple potentials such as the square well or harmonic oscillator. However, many physical systems feature **exotic potentials** that are highly non-polynomial, asymmetric, long-ranged, or even fractal-like. Examples include:

- Cusp-like potentials  
  $$
  V(x) = V_0 |x|^{\alpha}, \quad 0 < \alpha < 1
  $$

- Exponential wells  
  $$
  V(x) = -V_0 e^{-|x|/a}
  $$

- Soft barriers of the form  
  $$
  V(x) = \frac{V_0}{1 + x^2}
  $$

- Fractal or piecewise rough potentials

- Potentials defined by Bessel or Airy functions

The tunneling properties of such potentials are significantly richer than in textbook systems. This project explores how tunneling depends on the structural features of exotic potentials and investigates the breakdown of the WKB approximation in nonstandard regimes.

---

## 2. Research Problem

Let $ V(x) $ be an exotic one-dimensional potential supporting at least one bound state. Relevant examples include:

- **Cusp potential**  
  $$
  V(x) = V_0 |x|^{1/2}
  $$

- **Exponential well**  
  $$
  V(x) = -V_0 e^{-|x|/a}
  $$

- **Soft-barrier potential**  
  $$
  V(x) = \frac{V_0}{1 + x^2}
  $$

- **Fractal well**  
  A piecewise or iterated function with small-scale rough structure (optional)

One primary potential will be chosen, with an optional second for comparison.

### (a) Bound States and Wavefunction Structure

Compute numerically the bound states $ E_n $ and wavefunctions $ \psi_n(x) $ using finite-difference or shooting methods. Analyze:

- How the exotic shape influences localization  
- The dependence of forbidden-region penetration on $ n $  
- The presence of tails or irregular oscillations induced by the potential’s structure

---

### (b) WKB Tunneling in Exotic Potentials

Apply the WKB approximation to estimate the tunneling probability:

$$
T_n^{\text{WKB}} \approx \exp\left[
-\frac{2}{\hbar}
\int_{x_1}^{x_2}
\sqrt{2m(V(x) - E_n)} \, dx
\right]
$$

Explore:

- How the barrier’s curvature (or lack of smoothness) affects WKB accuracy  
- Cases where WKB fails (cusp at $ x = 0 $, long-range tails)  
- Comparison between tunneling from excited versus ground states  

---

### (c) External Field and Modified Barriers

Apply a static electric field:

$$
V(x) \rightarrow V(x) - Fx, \quad F > 0
$$

Investigate:

- How exotic barriers deform under tilting  
- Onset of over-the-barrier ionization  
- Conditions under which the barrier becomes ultrasoft or singular  

---

### (d) Time-Dependent Schrödinger Dynamics

Solve the time-dependent Schrödinger equation numerically for initially bound states:

$$
i\hbar \frac{\partial \psi(x,t)}{\partial t}
=
\left[
-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
+ V(x) - Fx
\right]\psi(x,t)
$$

Extract:

- Escape rates $ \Gamma_n $  
- Probability flux through the barrier  
- Deviations from WKB predictions due to exotic potential features  
- Possible appearance of resonances or quasi-bound states  

---

### (e) Research Questions

- How do exotic features of the potential influence tunneling scaling laws?  
- What structural features make WKB fail?  
- How do long-range potentials differ from compact-support barriers?  
- Do fractal or rough potentials induce anomalous tunneling behavior?  
- What numerical challenges arise from singular or cusp-like potentials?  

---

## 3. Expected Outcomes

- Implement numerical solvers for the stationary and time-dependent Schrödinger equation  
- Compute tunneling probabilities for nonstandard potentials  
- Evaluate the applicability and limits of the WKB approximation  
- Analyze how exotic geometries modify tunneling behavior  

---

## Suggested References

- L. D. Landau and E. M. Lifshitz, *Quantum Mechanics*  
- C. Cohen-Tannoudji et al., *Quantum Mechanics*  
- D. J. Griffiths and D. F. Schroeter, *Introduction to Quantum Mechanics*
