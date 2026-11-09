# SOT-MRAM Physics & OOMMF Implementation Reference

## Critical Current Density Formulas

### SOT (Spin-Orbit Torque)

For in-plane SOT with an in-plane bias field H_x:

```
Jc = (2 e μ0 Ms t_FL H_k) / (ℏ θ_SH)
```

where:
- `θ_SH`: spin Hall angle of the heavy metal (W: ~0.3, Ta: ~0.15, Pt: ~0.3)
- `H_k`: effective anisotropy field = 2K_u/M_s - M_s (PMA)
- `t_FL`: free layer thickness
- `μ0`: vacuum permeability

With field-like torque ratio ξ_FL = H_FL/H_DL, the effective critical field is modified:

```
H_k_eff = H_k * sqrt(1 + ξ_FL²)   (approximate)
```

### Field-Free SOT Switching

Requires symmetry breaking. Common approaches:
1. **STT assist**: small perpendicular current breaks symmetry
2. **Shape anisotropy**: elongated free layer
3. **Exchange bias**: antiferromagnetic coupling
4. **Wedged layer**: thickness gradient

For STT-assisted SOT:
```
J_SOT,c = Jc_SOT * (1 - J_STT/J_STT,c)   (approximate)
```

## SOT Torque Components

### Damping-Like Torque (DLT)

```
τ_DL = (γ ℏ θ_SH J / (2 e μ0 M_s t_FL)) * (σ × m) × m
```

Equivalent field: `H_DL = H_DL_mag * (σ × m)` where `H_DL_mag = ℏ θ_SH J / (2 e μ0 M_s t_FL)`

### Field-Like Torque (FLT)

```
τ_FL = (γ ℏ θ_SH J ξ_FL / (2 e μ0 M_s t_FL)) * σ × m
```

Equivalent field: `H_FL = ξ_FL * H_DL_mag * σ`

## OOMMF Implementation Options

### Option 1: Equivalent External Field (no extension needed)

Add SOT as `Oxs_UZeeman` or `Oxs_ScriptVectorField` with time-dependence.
Limitation: H_DL depends on m (σ × m), requires script field or multi-stage.

### Option 2: Custom SOT Extension

Known OOMMF SOT extensions:
- `oommf-sot` (various GitHub repos)
- Southampton SOT modules
- Combined STT+SOT via `Oxs_SpinXferEvolve` modifications

### Option 3: Ubermag (oommfc)

- Native `mm.Slonczewski` for STT
- SOT via `mm.Zeeman(H=time_function)` with m-dependent field
- Auto-generates MIF, calls OOMMF, reads OVF via `micromagneticdata`

## Jc Extraction Methods

### Method 1: Coarse Scan + Bisection
1. Log-spaced J scan (10-20 points)
2. Identify threshold J where mz crosses zero
3. Geometric bisection (5 iterations → ~0.1% accuracy)

### Method 2: Pulse Width Dependence
```
Jc(τ) = Jc0 * (1 + τ0/τ)   (thermal activation model)
```
Extract Jc0 (zero-temperature limit) by fitting Jc vs 1/τ.

### Method 3: Analytical Macrospin Estimate (quick screening)
Use the formula above to estimate Jc before running full micromagnetics.

## Switching Criteria

| Criterion | Definition | Use case |
|---|---|---|
| `mz_zero` | mz_final < 0 | Standard, robust |
| `mz_threshold` | mz_final < threshold | Partial switching |
| `energy` | Final energy < initial | Metastable states |

## Typical Parameter Ranges

| Parameter | Range | Notes |
|---|---|---|
| Ms | 0.8–1.5e6 A/m | CoFeB |
| Ku | 0.5–2e6 J/m³ | PMA |
| alpha | 0.005–0.03 | Damping |
| θ_SH | 0.1–0.5 | Heavy metal |
| ξ_FL | 0.05–0.5 | Field-like ratio |
| t_FL | 0.8–2 nm | Free layer |
| Jc | 1e11–5e12 A/m² | SOT critical |
