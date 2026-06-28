# COMFA-OA Thermal Comfort Index for Older Adults

COMFA-OA is a thermal comfort estimation tool designed for older adults. The model outputs the thermal comfort level as an **energy budget** (W/m²), with thresholds listed in Table¹ below.

An optional equivalent temperature metric — **CET (COMFA-OA Equivalent Temperature)** — represents the air temperature that would produce the same energy budget under standardized indoor-like meteorological conditions².

---

## About This Repository

This repository provides the code implementation of the COMFA-OA (Outdoor Thermal Comfort Model for Older Adults), as described in:

Li, X., Zhang, Y., Sang, H., Lee, C., Sullivan, W. C., Maddock, J., Li, D., & Brown, R. D. (2025). *A novel thermal comfort model for older adults – development and validation of the COMFA-OA model.*  
Building and Environment, 113758.  
https://doi.org/10.1016/j.buildenv.2025.113758

---

## Features

- Extends the original COMFA model with physiological adjustments for adults aged 55+.
- Required microclimate inputs: **Ta**, **Tmrt**, **Ws**, **Rh**.
- Required human parameters: **Height**, **Weight**, **clo**, **METs**, **Age**, **Sex**, **bodyPosture**.
- Incorporates major energy balance components: radiation, convection, evaporation, conduction.

---

## Calculation Options

### Option 1 — Run the calculator via code

- `1_COMFA_OA_Rmarkdown`
- `1_COMFA_OA_Python`

### Option 2 — Use the web calculator

[https://comfa.shinyapps.io/10_shinnyapp/](https://comfa.shinyapps.io/COMFA_OA_calculator/)

---

## Input Template (2_COMFA_OA_example_input.csv)

The COMFA-OA model requires microclimate, physiological, and demographic inputs to compute outdoor thermal comfort for older adults.  
The following variables appear in this repository’s example dataset.

### Variable Description

| Variable | Description |
|----------|-------------|
| Ta | Air temperature (°C) |
| Rh | Relative humidity (%) |
| clo | Clothing insulation (clo) |
| Tmrt | Mean radiant temperature (°C) |
| METs | Standard metabolic rate (met) |
| Age | Participant age (years) |
| sex | Biological sex (1 = male, 2 = female) |
| Ws | Wind speed (m/s) |
| Height | Height (m) |
| Weight | Weight (kg) |
| bodyPosture | Participant posture (e.g., sitting, standing) |


### Example Input Rows

| Ta    | Rh   | clo  | Tmrt  | METs | Age | sex | Ws     | Height | Weight | bodyPosture |
|-------|------|------|-------|------|-----|-----|--------|--------|--------|-------------|
| 21.05 | 43.1 | 0.73 | 20.37 | 1    | 86  | 1   | 0.955  | 1.8796 | 95.25  | sitting     |
| 31.34 | 60.5 | 0.34 | 73.86 | 1    | 86  | 1   | 1.451  | 1.8796 | 95.25  | sitting     |
| 35.63 | 43.7 | 0.33 | 48.54 | 1    | 56  | 2   | 1.191  | 1.6256 | 63.50  | sitting     |


---
@article{Li2026_COMFAOA,
  title   = {A novel thermal comfort model for older adults – development and validation of the COMFA-OA model},
  author  = {Li, Xiaoyu and Zhang, Yue and Sang, Huiyan and Lee, Chanam and Sullivan, William C. and Maddock, Jay E. and Li, Dongying and Brown, Robert D.},
  journal = {Building and Environment},
  volume  = {287},
  number  = {Part A},
  pages   = {113758},
  year    = {2026},
  issn    = {0360-1323},
  doi     = {10.1016/j.buildenv.2025.113758},
  url     = {https://www.sciencedirect.com/science/article/pii/S0360132325012284}
}

---

## Notes

### ¹ COMFA-OA Threshold Table

| Predicted TSV | COMFA-OA (W/m²) | Thermal Stress Category |
|---------------|----------------|--------------------------|
| -3 (Cold) | < -192 | High cold stress |
| -2 (Cool) | -192 to -121 | Medium cold stress |
| -1 (Slightly cool) | -120 to -50 | Slight cold stress |
| 0 (Neutral) | -49 to 22 | Thermally neutral |
| 1 (Slightly warm) | 23 to 93 | Slight heat stress |
| 2 (Warm) | 94 to 165 | Medium heat stress |
| 3 (Hot) | > 165 | High heat stress |

---

### ² Standardized Conditions for CET

**Meteorological conditions**
- Ta = Tmrt  
- Relative humidity = 50%  
- Wind speed = 0.3 m/s (~0.5 m/s at 10 m)

**Human parameters**
- Age = 65 years  
- Height = 1.75 m  
- Weight = 75 kg  
- Sex = Male  
- Standard metabolic rate = 2.3 METs
- Clothes insulation = 0.9 clo  
