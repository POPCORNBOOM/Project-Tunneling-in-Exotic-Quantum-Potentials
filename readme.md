
# Quantum Tunneling in Cusp Potentials: Bound States, Field-Induced Ionization, and Time-Dependent Dynamics

[[PDF version]](/Quantum%20Tunneling%20in%20Cusp%20Potentials%20Bound%20States,%20Field-Induced%20Ionization,%20and%20Time-Dependent%20Dynamics.pdf)

**Yigu Wang**  
*January 6, 2026*

## Abstract

I investigate quantum tunneling phenomena in cusp potentials of the form $V(x) = V_0|x|^\alpha$ with $0 < \alpha < 1$. The singular nature of the cusp at $x=0$ introduces unique features in both bound state structures and tunneling dynamics. I employ finite difference methods to solve for eigenspectra, analyze forbidden region penetration, and examine field-induced barrier modification under uniform external fields. Using the WKB approximation and direct time-dependent Schrödinger equation (TDSE) simulations via split-operator methods, I characterize escape rates and transmission coefficients. My results reveal critical field strengths for over-the-barrier ionization, the formation of ultrasoft barriers, and significant deviations from standard WKB predictions in the presence of singular potentials.

## 1. Introduction

Quantum tunneling is a fundamental phenomenon in quantum mechanics, manifesting in diverse physical systems ranging from atomic ionization to nuclear decay. The study of tunneling through exotic potential barriers provides insights into non-classical behavior and serves as a testing ground for approximation methods such as the Wentzel-Kramers-Brillouin (WKB) approach.

Most introductory quantum mechanics courses focus on simple potentials such as the square well or harmonic oscillator. However, many physical systems feature **exotic potentials** that are highly non-polynomial, asymmetric, long-ranged, or exhibit singular features. The tunneling properties of such potentials are significantly richer than in textbook systems and often lead to the breakdown of standard approximation schemes.

In this work, I focus on the cusp potential:

$$
V(x) = V_0 |x|^{\alpha}, \quad 0 < \alpha < 1
$$

This potential exhibits a sharp corner at $x=0$, introducing a singularity that significantly affects wavefunction behavior. Unlike smooth potentials, the infinite derivative at the origin creates unique challenges for both analytical approximations and numerical simulations. The cusp potential serves as an important test case for understanding how structural features of exotic potentials influence tunneling behavior and the validity limits of semiclassical methods.

My investigation is structured as follows: In Section 2, I solve for bound states and analyze wavefunction characteristics. Section 3 examines potential barrier deformation under external fields and identifies critical field strengths for ionization. Section 4 presents time-dependent dynamics and escape rate calculations. I conclude with a discussion of my findings and their implications for understanding tunneling in exotic potentials.

## 2. Bound States and Wavefunction Structure

### 2.1 Methodology

I solve the one-dimensional time-independent Schrödinger equation:

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi
$$

using finite difference discretization on a spatial grid of length $L$ with $N$ points. This converts the differential equation into a sparse matrix eigenvalue problem, which I solve using `scipy.sparse.linalg.eigsh`. The finite difference method is particularly well-suited for handling the cusp singularity at $x=0$, as it does not require explicit evaluation of derivatives at that point.

### 2.2 Potential Profile and Eigenspectrum

Figure 1 shows the cusp potential along with the first several energy eigenstates. The wavefunctions are normalized, and their probability densities reveal the spatial localization of each quantum state. The singular point at $x=0$ is clearly visible as a sharp vertex in the potential profile.

![Cusp potential with eigenstates](figures/part_a_plot_potential_and_states.png)  
*Figure 1: Cusp potential $V(x) = V_0|x|^\alpha$ with the first several bound state energies + scaled eigenfunctions and corresponding wavefunctions. The singular point at $x=0$ is clearly visible.*

### 2.3 Ground State Evolution with Potential Depth

To understand how the potential shape influences eigenstate properties, I performed a parameter scan over the potential depth $V_0$. As $V_0$ increases, the ground state becomes increasingly localized near the origin. This is quantified by two metrics:

- **Spatial spread** $\sigma = \sqrt{\langle x^2 \rangle - \langle x \rangle^2}$: decreases with increasing $V_0$, indicating stronger confinement
- **Inverse participation ratio** (IPR) $= \int |\psi(x)|^4 dx$: increases with $V_0$, quantifying the degree of localization

![Ground state wavefunction evolution](figures/ground_state_evolution.gif)  
*Figure 2: Ground state wavefunction evolution as potential depth $V_0$ increases. Four snapshots show the progressive contraction of the wavefunction toward the origin as the potential well deepens. This demonstrates how the exotic shape of the cusp potential influences localization properties.*

### 2.4 Forbidden Region Probability

For each eigenstate with energy $E_n$, I define the classically forbidden region as the set of points where $V(x) > E_n$. The forbidden region probability is:

$$
P_{\text{forbidden}}^{(n)} = \int_{V(x) > E_n} |\psi_n(x)|^2 \, dx
$$

This quantity measures the extent of quantum tunneling into classically inaccessible regions.

Figure 3 shows that higher energy states have *smaller* forbidden region probabilities. This counterintuitive result arises because the cusp potential grows rapidly away from the origin (as $|x|^\alpha$ with $\alpha < 1$). Higher energy states, despite having more energy, remain primarily confined within the wide potential well where $V(x) < E_n$. The effective width of the allowed region increases faster than the tunneling tails extend into the forbidden region.

![Forbidden region probability](figures/part_a_forbidden_region_probability_vs_n.png)  
*Figure 3: Forbidden region probability versus quantum number $n$. Higher energy states exhibit reduced tunneling into classically forbidden regions, a feature specific to the exotic shape of the cusp potential.*

### 2.5 Penetration Depth Analysis

I define the penetration depth as the distance from the classical turning point at which the probability density decays to $1/e$ of its value at the turning point. This provides a characteristic length scale for wavefunction penetration into the forbidden region.

Figure 4 demonstrates that penetration depth increases with quantum number $n$, suggesting that higher energy states can extend further into the forbidden region despite having lower total forbidden probabilities.

This behavior can be understood from the effective shape of the barrier near the turning point: as the bound-state energy increases, the potential becomes locally flatter in this region (since $d^2V/dx^2 \propto |x|^{\alpha-2}$ decreases for larger $|x|$). As a result, the exponential decay of the wavefunction is weaker, leading to a larger characteristic penetration depth even though the integrated forbidden-region probability decreases.

![Penetration depth](figures/part_a_forbidden_penetration_depth_vs_n.png)  
*Figure 4: Penetration depth versus quantum number $n$. The characteristic decay length increases for higher excited states due to the local flattening of the potential barrier at larger distances.*

The spatial decay rate in the forbidden region is shown in Figure 5. As $n$ increases, the effective decay becomes slower, consistent with the flattening of the potential barrier at higher energies. This trend is distinct from standard polynomial potentials and reflects the exotic character of the cusp potential.

![Decay rate](figures/part_a_forbidden_penetration_decay_rate.png)  
*Figure 5: Decay rate of wavefunction amplitude in the forbidden region. Higher states exhibit slower decay due to the reduced effective barrier, demonstrating how the potential's structure influences tunneling characteristics.*

## 3. Field-Induced Barrier Deformation and Ionization

### 3.1 Tilted Potential Under External Field

I apply a uniform external electric field $F$, modifying the potential to:

$$
V_F(x) = V_0|x|^\alpha - Fx
$$

This linear term tilts the potential, lowering the barrier on one side (here, $x > 0$) and potentially enabling over-the-barrier escape. The field-modified potential allows us to investigate how exotic barriers deform under external perturbations and the onset of field ionization.

### 3.2 Barrier Evolution Under Increasing Field

Figure 6 illustrates the progressive barrier suppression as the field strength increases. The singular point at $x=0$ persists regardless of field strength—this is an intrinsic property of the cusp potential that remains even under external field modification.

![Barrier deformation](figures/part_c_tilted_barrier_animation.gif)  
*Figure 6: Potential barrier deformation under increasing external field $F$. Four snapshots illustrate the progressive barrier suppression. The barrier height decreases until over-the-barrier ionization becomes possible. Note that the singular point at $x=0$ persists throughout, demonstrating the robustness of this exotic feature.*

### 3.3 WKB Transmission Coefficient

I compute the WKB transmission probability:

$$
T_{\text{WKB}} = \exp\left(-\frac{2}{\hbar}S\right), \quad S = \int_{x_1}^{x_2} \sqrt{2m(V_F(x) - E_n)} \, dx
$$

where $x_1$ and $x_2$ are the classical turning points. The WKB action $S$ quantifies the barrier penetration difficulty and provides a semiclassical estimate of tunneling probability.

Figure 7 shows $\ln T$ versus field strength. As $F$ increases, the transmission coefficient rises exponentially due to barrier suppression. At a critical field $F \approx 1.87$, the barrier vanishes and $T \to 1$, marking the transition to over-the-barrier ionization.

![Transmission coefficient (log scale)](figures/part_c_transmission_vs_field_ln.png)  
*Figure 7: Logarithm of WKB transmission coefficient versus external field strength. The vertical asymptote indicates the onset of over-the-barrier ionization, where the barrier maximum drops below the bound state energy.*

Figure 8 presents the same data on a linear scale, clearly showing the transition to complete transmission in the ultrasoft barrier regime.

![Transmission coefficient (linear scale)](figures/part_c_transmission_vs_field_linear.png)  
*Figure 8: WKB transmission coefficient (linear scale) versus field strength. The rapid increase marks the ultrasoft barrier regime preceding ionization, where the barrier height becomes comparable to or less than the thermal energy scale.*

### 3.4 Critical Field for Ionization

The critical field $F_{\text{ionization}}$ is determined by the condition that the barrier maximum equals the eigenstate energy. For the cusp potential, I find the barrier maximum by solving $dV_F/dx = 0$:

$$
\alpha V_0 x_{\text{top}}^{\alpha-1} = F \quad \Rightarrow \quad x_{\text{top}} = \left(\frac{F}{\alpha V_0}\right)^{1/(\alpha-1)}
$$

Setting $V_F(x_{\text{top}}) = E_n$ yields the critical field:

$$
F_{\text{ionization}} = \alpha \left( \frac{(1-\alpha)^{(1-\alpha)} V_0}{E_n^{(1-\alpha)}} \right)^{1/\alpha}
$$

This analytical formula provides a precise prediction for the field strength required to ionize a given bound state.

### 3.5 State-Dependent Ionization Thresholds

Figure 9 demonstrates that higher energy states require lower fields for ionization. This reflects the physical intuition that particles with greater energy need less external assistance to overcome the barrier. The dependence is monotonic and follows the scaling predicted by the critical field formula.

![Critical field vs state](figures/part_c_ionization_field_vs_state.png)  
*Figure 9: Critical ionization field versus quantum state number $n$. Higher excited states ionize at progressively weaker fields, as they are already closer in energy to the barrier top. This trend is characteristic of exotic potentials with rapidly growing tails.*

### 3.6 Ultrasoft and Singular Barriers

**Ultrasoft barrier formation:** As $F \to F_{\text{ionization}}$, the barrier height $\Delta V = V_{\text{top}} - E_n$ approaches zero. In this regime:

- The WKB action integral $S \to 0$
- Transmission coefficient $T_{\text{WKB}} \to 1$
- The barrier becomes "ultrasoft" in the sense that its effective height is negligible
- The WKB approximation becomes questionable, as it assumes a well-defined barrier

The ultrasoft regime is characterized by extremely rapid variation of tunneling probability with field strength, making it a regime of particular interest for field-ionization experiments.

**Singular barrier characteristics:** The singularity at $x=0$ is an intrinsic feature of the cusp potential with $\alpha < 1$:

$$
\left.\frac{dV}{dx}\right|_{x \to 0} = \alpha V_0 \, \text{sgn}(x) \, |x|^{\alpha-1} \to \infty
$$

This divergent derivative persists even in the presence of the linear field term. Consequently:

- When turning points approach $x=0$, numerical WKB integration becomes unstable
- Second-order WKB corrections (involving $d^2V/dx^2$) diverge at the singular point
- The connection formulas used in WKB theory may break down near the singularity
- Special care is needed in numerical implementations to avoid spurious results

The ultrasoft and singular features are independent phenomena: ultrasoft barriers arise from field tuning (controllable and reversible), while singularity is inherent to the potential form (unavoidable structural feature).

## 4. Time-Dependent Schrödinger Equation Dynamics

### 4.1 Numerical Method

I solve the time-dependent Schrödinger equation:

$$
i\hbar \frac{\partial \psi(x,t)}{\partial t} = \left[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V_F(x)\right]\psi(x,t)
$$

using the split-operator method with fast Fourier transforms (FFT). This method alternates between position-space evolution (for the potential term) and momentum-space evolution (for the kinetic term), achieving spectral accuracy for smooth wavefunctions.

To prevent artificial reflections from spatial boundaries, I employ complex absorbing potentials (CAP) near the edges of the computational domain:

$$
V_{\text{CAP}}(x) = -i W(x), \quad W(x) = 
\begin{cases}
W_0 \left(\frac{|x| - x_{\text{abs}}}{L/2 - x_{\text{abs}}}\right)^2 & \text{if } |x| > x_{\text{abs}} \\
0 & \text{otherwise}
\end{cases}
$$

The CAP absorbs outgoing probability flux, simulating an open quantum system.

### 4.2 Initial Conditions and Tilted Potential

I initialize the system in an eigenstate $\psi_n$ of the *unperturbed* potential (without field) and then evolve under the tilted potential $V_F(x)$. This mimics a sudden switching-on of the external field, a scenario relevant to field-ionization experiments.

Figure 10 shows the initial configuration with the tilted potential and several bound state energies.

![Tilted potential](figures/part_d_tilted_potential_and_bound_states.png)  
*Figure 10: Tilted potential $V_F(x) = V_0|x|^\alpha - Fx$ with initial bound state wavefunctions. The barrier on the right side is suppressed by the external field, enabling tunneling escape. The initial states are eigenstates of the unperturbed potential.*

### 4.3 Survival Probability and Escape Rate

I define the in-well survival probability as:

$$
P_{\text{inwell}}(t) = \int_{x \leq x_{\text{top}}} |\psi(x,t)|^2 \, dx
$$

where $x_{\text{top}}$ is the position of the barrier maximum. This quantity measures the probability that the particle remains confined within the potential well region.

Figure 11 shows the time evolution of $P_{\text{inwell}}(t)$ for the $n=3$ state. The exponential decay indicates a constant escape rate, characteristic of quantum tunneling from a metastable state.

![Survival probability](figures/part_d_survival_vs_time.png)  
*Figure 11: In-well probability versus time for the fourth energy level ($n=3$). The exponential decay indicates a constant escape rate, validating the use of a simple exponential model for tunneling dynamics.*

### 4.4 Exponential Decay Fitting

Assuming exponential decay $P_{\text{inwell}}(t) \approx A e^{-\Gamma t}$, I fit the latter half of the time series (to avoid transient effects from the initial non-eigenstate evolution) to:

$$
\ln P_{\text{inwell}}(t) = \ln A - \Gamma t
$$

Linear regression yields the escape rate $\Gamma$. For the parameters used ($F = 1.5$, state $n=3$), I obtain $\Gamma \approx 2.20 \times 10^{-4}$ (in atomic units).

![Log survival fit](figures/part_d_log_survival_linear_fit.png)  
*Figure 12: Logarithm of in-well probability versus time with linear fit. The slope gives the escape rate $\Gamma$. The excellent linear fit confirms the exponential decay assumption in the quasi-steady-state regime.*

### 4.5 Probability Flux Analysis

The probability current density is:

$$
j(x,t) = \frac{\hbar}{m} \, \text{Im}\left[\psi^*(x,t) \frac{\partial \psi(x,t)}{\partial x}\right]
$$

This quantity represents the local flow of probability and can be used to verify the escape rate independently.

I evaluate $j(x_{\text{top}}, t)$ at the barrier top and compare the ratio $j/P_{\text{inwell}}$ over time. In the steady-state regime, this ratio should approach the escape rate $\Gamma$, providing a consistency check.

![Flux ratio](figures/part_d_flux_over_survival.png)  
*Figure 13: Ratio of probability flux at barrier top to in-well probability. The stabilization of this ratio confirms the exponential decay assumption and provides an independent estimate of the escape rate consistent with the survival probability analysis.*

### 4.6 Time Evolution Visualization

Figure 14 provides a visual representation of the wavefunction dynamics, showing the real part, imaginary part, and probability density at four time snapshots. The visualization confirms several key features:

- The wavefunction leaks continuously over the barrier, with probability flowing from left to right
- The CAP effectively absorbs outgoing flux without spurious reflections
- The in-well region remains relatively undisturbed by boundary effects
- The phase evolution (visible in the real and imaginary parts) is consistent with quantum tunneling dynamics

![Time evolution snapshots](figures/part_d_time_evolution_animation.gif)  
*Figure 14: Time evolution of the wavefunction at four time snapshots. Each panel shows the real part (green), imaginary part (purple), and probability density (blue) of $\psi(x,t)$. The tilted potential (black) and energy level (red dashed line) are indicated for reference. The absorbing cap region is shown in orange shading. The wavefunction continuously leaks over the barrier and is absorbed at the boundaries, with minimal reflection artifacts.*

### 4.7 Comparison with WKB Predictions

I computed escape rates for all bound states below the barrier top and plotted $\ln \Gamma$ versus the WKB action $S$. According to the WKB approximation, the escape rate should be proportional to the tunneling probability:

$$
\Gamma \propto T_{\text{WKB}} = \exp\left(-\frac{2}{\hbar}S\right)
$$

implying a linear relationship:

$$
\ln \Gamma \approx -\frac{2}{\hbar}S + \text{const}
$$

with slope $-2/\hbar = -2$ in atomic units.

![lnΓ vs S](figures/part_d_lnGamma_vs_S.png)  
*Figure 15: Logarithm of TDSE-derived escape rate versus WKB action. The linear fit yields $\ln\Gamma = -0.0903S - 7.9411$ with correlation coefficient $r = -0.86$. The significant deviation from the expected slope of $-2$ indicates breakdown of the simple WKB approximation.*

The fitted slope is approximately $-0.0903$, significantly deviating from the theoretical value of $-2$. This large discrepancy can be attributed to several factors:

- **Singular potential structure:** The cusp at $x=0$ violates the smoothness assumptions underlying WKB theory. The infinite derivative leads to breakdown of the asymptotic expansions used in deriving WKB formulas.
- **Quantum reflection:** Near the singularity, quantum reflection effects become significant. These are not captured by the simple WKB transmission coefficient, which assumes purely classical motion away from turning points.
- **Low-lying states:** For states with small quantum numbers, the semiclassical approximation is inherently less accurate. The de Broglie wavelength is comparable to the potential length scale, violating the WKB regime condition.

This analysis demonstrates that exotic potentials with singular features require beyond-WKB treatments for quantitative accuracy. While WKB provides useful qualitative insights, careful numerical simulations are essential for precise predictions.

## 5. Discussion and Conclusions

I have conducted a comprehensive investigation of quantum tunneling in cusp potentials, combining analytical WKB methods with numerical time-dependent simulations. This work addresses key research questions about how exotic features of quantum potentials influence tunneling behavior:

### 5.1 Key Findings

1. **Bound state properties and exotic potential features:** The singular point at $x=0$ leads to distinctive wavefunction characteristics not seen in smooth potentials. Higher excited states exhibit larger penetration depths but paradoxically lower forbidden region probabilities due to the rapid growth of the potential. This demonstrates how exotic potential geometries can produce counterintuitive quantum behavior.

2. **Field-induced ionization and barrier modification:** External fields systematically suppress the barrier, leading to ultrasoft barrier formation near the critical field $F_{\text{ionization}}$. Higher energy states require progressively weaker fields for ionization, following the analytical prediction. The singular point at $x=0$ persists under field modification, representing a robust structural feature.

3. **WKB accuracy and breakdown mechanisms:** While WKB provides qualitative insights into tunneling trends, quantitative predictions deviate significantly from TDSE results. The fitted slope of $-0.0903$ versus the theoretical $-2$ indicates that the cusp singularity fundamentally undermines the semiclassical approximation. This answers the research question: *What structural features make WKB fail?* The answer is that non-analytic features (infinite derivatives, cusps) violate the smoothness assumptions essential to WKB theory.

4. **Numerical techniques for exotic potentials:** The split-operator FFT method with complex absorbing potentials proves effective for studying escape dynamics, even in the presence of singularities. The convergence of flux-based and survival-based escape rates validates the numerical implementation. However, special care is needed near singular points to avoid numerical instabilities.

5. **Tunneling scaling laws:** The relationship between escape rate and WKB action shows significant deviation from the standard exponential scaling. This suggests that exotic potentials may exhibit anomalous tunneling behavior, potentially leading to novel phenomena in field-ionization experiments or quantum transport through non-standard barriers.

### 5.2 Comparison with Standard Potentials

Unlike polynomial potentials (square wells, harmonic oscillators, quartic barriers), the cusp potential exhibits several distinctive features:

- Non-polynomial growth leads to counterintuitive scaling of forbidden-region probabilities
- The singular point creates localized quantum reflection effects
- Field-induced barrier modification produces more dramatic changes in tunneling probability
- WKB accuracy degrades more severely than in smooth potential cases

These findings suggest that the broader class of exotic potentials (exponential wells, soft barriers, fractal potentials) may harbor rich tunneling physics worthy of further investigation.

### 5.3 Future Directions

Several promising directions for future work emerge from this study:

- **Higher-order WKB corrections:** Implementing uniform approximations or complex-plane techniques near the singularity may improve agreement with numerical results.
- **Multi-dimensional generalizations:** Extending the cusp potential to higher dimensions could reveal new features.
- **Quantum reflection analysis:** A detailed study of reflection coefficients near the cusp could quantify the deviation from pure tunneling behavior.
- **Comparison with other exotic potentials:** Systematic comparison with exponential wells ($V \sim e^{-|x|/a}$), soft barriers ($V \sim 1/(1+x^2)$), or fractal potentials could identify universal features versus potential-specific behavior.

### 5.4 Broader Implications

The cusp potential serves as an important test case for quantum tunneling theory, highlighting the limitations of semiclassical approximations and the rich physics arising from potential singularities. The significant deviations from WKB predictions underscore the need for careful numerical validation when applying semiclassical methods to exotic systems.

From a fundamental perspective, this work demonstrates that quantum mechanics is exquisitely sensitive to the detailed structure of the potential energy landscape. Features that might seem like minor technical complications (a cusp instead of a smooth minimum) can profoundly alter observable properties like escape rates and ionization thresholds.

## Acknowledgments

I would like to thank Prof. Ciappina for his helpful comments and suggestions provided via email, which greatly improved the clarity and quality of this work.

I acknowledge the course *Quantum Physics II* for providing the opportunity to explore this concrete problem in quantum tunneling through an independent research project.

The numerical simulations and data analysis were performed using custom Python codes built on `scipy`, `numpy`, and `matplotlib`. I am grateful to the developers of these invaluable computational tools.

The full source code, data, and animations are available at [https://github.com/POPCORNBOOM/Project-Tunneling-in-Exotic-Quantum-Potentials](https://github.com/POPCORNBOOM/Project-Tunneling-in-Exotic-Quantum-Potentials).

## References

1. L.D. Landau and E.M. Lifshitz, *Quantum Mechanics: Non-Relativistic Theory*, 3rd ed., Pergamon Press (1977).
2. C. Cohen-Tannoudji, B. Diu, and F. Laloë, *Quantum Mechanics*, Wiley-VCH (1977).
3. D.J. Griffiths and D.F. Schroeter, *Introduction to Quantum Mechanics*, 3rd ed., Cambridge University Press (2018).