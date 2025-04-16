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

### 1. Euler Discretization

Time step:
$$
\Delta t = \frac{T}{N_{\text{steps}}}
$$

Brownian motion increments:
$$
\Delta W_k = \sqrt{\Delta t} \cdot Z_k
$$

---

### 2. Variance Process (Heston)

Updated as:
$$
V_{t+\Delta t} = V_t + \kappa (v_{\text{bar}} - V_t)\Delta t + \gamma \sqrt{V_t} \Delta W_v + \gamma \rho_{vrd} \eta_d B_d(t,T) \sqrt{V_t} \Delta t
$$

Ensure:
$$
V_t \geq 0
$$

---

### 3. FX Process

Updated as:
$$
FX_{t+\Delta t} = FX_t \left(1 + \sqrt{V_t} \Delta W_x - \eta_d B_d(t,T) \Delta W_{rd} + \eta_f B_f(t,T) \Delta W_{rf} \right)
$$

---

### 4. Discount Factor (Hull-White)

Domestic:
$$
B_d(t,T) = \frac{1 - e^{-\lambda_d(T-t)}}{\lambda_d}
$$

Foreign:
$$
B_f(t,T) = \frac{1 - e^{-\lambda_f(T-t)}}{\lambda_f}
$$

---

### 5. Variance Expectation (CIR Process)

Expected square root of variance:
$$
E[\sqrt{V(t)}] = \sqrt{2c(t)} \cdot \frac{\Gamma\left(\frac{1+\delta}{2}\right)}{\Gamma\left(\frac{\delta}{2}\right)} \cdot {}_1F_1\left(-\frac{1}{2}, \frac{\delta}{2}, -\frac{\bar{\kappa}(t)}{2}\right)
$$

Where:
- $ \delta = \frac{4\kappa v_{\text{bar}}}{\gamma^2} $
- $ c(t) = \frac{\gamma^2}{4\kappa}(1 - e^{-\kappa t}) $
- $ \bar{\kappa}(t) = \frac{4 \kappa v_0 e^{-\kappa t}}{\gamma^2(1 - e^{-\kappa t})} $

---

### 6. Characteristic Function (Heston Component)

Define:
$$
D_1 = \sqrt{(\kappa - \gamma \rho_{xv} i u)^2 + \gamma^2 (u^2 + iu)}
$$

Then:
$$
g = \frac{\kappa - \gamma \rho_{xv} i u - D_1}{\kappa - \gamma \rho_{xv} i u + D_1}
$$

C-function:
$$
C(u,\tau) = \frac{1 - e^{-D_1 \tau}}{\gamma^2(1 - g e^{-D_1 \tau})} (\kappa - \gamma \rho_{xv} i u - D_1)
$$

---

### 7. Characteristic Function (Full)

Let:
$$
A = I_1(u) + I_2(u)
$$

Then:
$$
\phi(u) = \exp\left(A + v_0 C(u,\tau)\right)
$$

With:
$$
I_1(u) = \int_0^\tau \left[\text{temp1}(z) + \text{temp2}(z,u) + \text{temp3}(z,u)\right] C(u,z)\,dz
$$

And:
$$
I_2(u) = (u^2 + iu) \int_0^\tau \zeta(\tau - z)\,dz
$$

---

## 🧪 Installation

```bash
git clone https://github.com/shubh123a3/FX-Option-By-Heston-Hull-White-Model
cd FX-Option-By-Heston-Hull-White-Model
pip install -r requirements.txt
