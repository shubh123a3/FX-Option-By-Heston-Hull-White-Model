# FX Option Pricing using Heston-Hull-White Model

This project simulates FX option pricing under stochastic volatility (Heston model) and stochastic domestic/foreign interest rates (Hull-White model). Both Monte Carlo simulation and characteristic function (CF) methods are implemented with detailed mathematical explanations.

## ✨ Project Highlights

- Heston Model for stochastic volatility  
- Hull-White Model for interest rates  
- Monte Carlo with Euler discretization  
- Full characteristic function implementation  
- Clean, vectorized Python code  
- Plotting volatility surface, forward paths, and more  

---

## 🧠 Mathematical Theory and Formulas

> Use GitHub with MathJax enabled for LaTeX rendering.

### 1. Euler Discretization

```math
\Delta t = \frac{T}{N_{\text{steps}}}
```

```math
\Delta W_k = \sqrt{\Delta t} \cdot Z_k
```

---

### 2. Variance Process (Heston)

```math
V_{t+\Delta t} = V_t + \kappa (v_{\text{bar}} - V_t)\Delta t + \gamma \sqrt{V_t} \Delta W_v + \gamma \rho_{vrd} \eta_d B_d(t,T) \sqrt{V_t} \Delta t
```

```math
V_t \geq 0
```

---

### 3. FX Process

```math
FX_{t+\Delta t} = FX_t \left(1 + \sqrt{V_t} \Delta W_x - \eta_d B_d(t,T) \Delta W_{rd} + \eta_f B_f(t,T) \Delta W_{rf} \right)
```

---

### 4. Discount Factor (Hull-White)

```math
B_d(t,T) = \frac{1 - e^{-\lambda_d(T-t)}}{\lambda_d}
```

```math
B_f(t,T) = \frac{1 - e^{-\lambda_f(T-t)}}{\lambda_f}
```

---

### 5. Variance Expectation (CIR Process)

```math
E[\sqrt{V(t)}] = \sqrt{2c(t)} \cdot \frac{\Gamma\left(\frac{1+\delta}{2}\right)}{\Gamma\left(\frac{\delta}{2}\right)} \cdot {}_1F_1\left(-\frac{1}{2}, \frac{\delta}{2}, -\frac{\bar{\kappa}(t)}{2}\right)
```

Where:

```math
\delta = \frac{4\kappa v_{\text{bar}}}{\gamma^2}
```

```math
c(t) = \frac{\gamma^2}{4\kappa}(1 - e^{-\kappa t})
```

```math
\bar{\kappa}(t) = \frac{4 \kappa v_0 e^{-\kappa t}}{\gamma^2(1 - e^{-\kappa t})}
```

---

### 6. Characteristic Function (Heston Component)

```math
D_1 = \sqrt{(\kappa - \gamma \rho_{xv} i u)^2 + \gamma^2 (u^2 + iu)}
```

```math
g = \frac{\kappa - \gamma \rho_{xv} i u - D_1}{\kappa - \gamma \rho_{xv} i u + D_1}
```

```math
C(u,\tau) = \frac{1 - e^{-D_1 \tau}}{\gamma^2(1 - g e^{-D_1 \tau})} (\kappa - \gamma \rho_{xv} i u - D_1)
```

---

### 7. Characteristic Function (Full)

Let:

```math
A = I_1(u) + I_2(u)
```

Then:

```math
\phi(u) = \exp\left(A + v_0 C(u,\tau)\right)
```

Where:

```math
I_1(u) = \int_0^\tau \left[\text{temp1}(z) + \text{temp2}(z,u) + \text{temp3}(z,u)\right] C(u,z)\,dz
```

```math
I_2(u) = (u^2 + iu) \int_0^\tau \zeta(\tau - z)\,dz
```

---

## 🧪 Installation

```bash
git clone https://github.com/shubh123a3/FX-Option-By-Heston-Hull-White-Model
cd FX-Option-By-Heston-Hull-White-Model
pip install -r requirements.txt
```

---

## 📈 Sample Outputs

- Forward FX paths (Monte Carlo)  
- Volatility surface  
- Characteristic function plots  
- CIR variance paths  

---

## 📚 References

- John C. Hull – *Options, Futures and Other Derivatives*  
- Steven Shreve – *Stochastic Calculus for Finance*  
- Heston (1993) – Closed-form solution for options with stochastic volatility  
- Hull & White (1990) – Pricing interest-rate derivatives  

---

## 🤝 Contributing

Feel free to fork, open issues, or submit PRs to improve or optimize the simulation.

---

## 🧠 Author

**Shubh Shrishrimal** – Quant Finance Student | CS + Finance Enthusiast  
 or check out more projects at [GitHub](https://github.com/shubh123a3)
