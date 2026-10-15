# twoSCH4BferLib — a CH₄/O₂ combustion library for OpenFOAM 14

A custom OpenFOAM library providing a **power-law transport package**
(`twoSCH4Bfer`) and an **equivalence-ratio-dependent Arrhenius reaction rate**
(`ArrheniusReactionRateExtend`) for two-stage methane combustion.

Built for **openfoam.org 14** (ported from ESI/OpenCFD 2506).

---

## What it does

This library adds two runtime-selectable components to OpenFOAM:

### 1. Transport model `twoSCH4Bfer`

Power-law dynamic viscosity with a constant Prandtl number (or constant thermal
conductivity):

$$\mu(T) = \mu_0 \left(\frac{T}{T_0}\right)^{\alpha}$$

- constant **Prandtl number** mode: $\kappa = \dfrac{C_p\,\mu}{Pr}$
- constant **thermal conductivity** mode: $\kappa = \kappa_0$

### 2. Reaction rate `ArrheniusReactionRateExtend`

Modified Arrhenius rate with an equivalence-ratio multiplier:

$$k = A\,T^{\beta}\,\exp\!\left(-\frac{T_a}{T}\right)\,f(\varphi)$$

- The multiplier $f(\varphi)$ is selected by the **reaction name** (`R1` or `R2`),
  i.e. each stage of the two-step mechanism has its own correction.
- The equivalence ratio $\varphi$ is read once from `constant/phiEq`.

This is intended for a **two-stage methane mechanism** of the form:

```
R1:  CH4 + O2 → CO + H2O      (fuel consumption, phi-corrected)
R2:  CO + O2 → CO2            (CO oxidation, phi-corrected)
```

---

## Features

- Full equation-of-state coverage (19 thermo/EoS combinations):
  - coefficient-mixing EoS: `perfectGas`, `PengRobinsonGas`,
    `incompressiblePerfectGas`, `rPolynomial`, `rhoConst`
  - single-component EoS: `perfectFluid`, `adiabaticPerfectFluid`, `icoPolynomial`
- Thermo: `janaf`, `hConst`, `eConst`, `hPolynomialThermo`
- Energy: `sensibleEnthalpy`, `sensibleInternalEnergy`
- Reaction types (runtime-selectable):
  - `irreversibleArrheniusTwoSCH4Bfer`
  - `reversibleArrheniusTwoSCH4Bfer`
  - `nonEquilibriumReversibleArrheniusTwoSCH4Bfer`
- Chemistry model registration (`Standard`) and chemistry reduction methods
  (`none`, `DRG`, `DRGEP`, `DAC`, `EFA`, `PFA`).

---

## Directory layout

```
OpenFOAM_org/
└── src/
    └── twoSCH4Bfer/
        ├── twoSCH4Bfer.H / .C / .I.H           # transport package
        ├── ArrheniusReactionRateExtend.H/.I.H  # reaction rate
        ├── forTwoSCH4Bfer.H                    # macro cascades (EoS/thermo)
        ├── twoSCH4BferThermos.C                # thermo registration
        ├── makeReactionsTwoS.C                 # reaction registration
        ├── makeChemistryModelsTwoS.C           # chemistry model registration
        ├── makeChemistryReductionMethodsTwoS.C # reduction methods
        ├── README.md                           # detailed docs
        └── Make/
            ├── files
            └── options
```

---

## Requirements

- **OpenFOAM 14** (openfoam.org), built with `g++` (C++14).
- The library is compiled into `$FOAM_USER_LIBBIN/twoSCH4BferLib.so`.

> **Environment note:** if your `~/.bashrc` sources a *different* OpenFOAM
> version (e.g. `v2506`), you must build/run in a clean environment that sources
> the org 14 `etc/bashrc` (see below).

---

## Installation (compile)

```bash
# Use a clean environment pointing to openfoam.org 14
env -i HOME="$HOME" USER="$USER" TERM="$TERM" \
  PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
  bash --noprofile --norc -c \
  'source /opt/openfoam14/etc/bashrc && \
   cd src/twoSCH4Bfer && \
   wmake libso'
```

The resulting library is written to:

```
$FOAM_USER_LIBBIN/twoSCH4BferLib.so
```

To rebuild cleanly:

```bash
wclean && wmake libso
```

---

## Usage

### 1. Load the library

In `system/controlDict`:

```cpp
libs
(
    "$FOAM_USER_LIBBIN/twoSCH4BferLib.so"
);
```

### 2. Select the thermo package

In `constant/physicalProperties`:

```cpp
thermoType
{
    type            hePsiThermo;      // or heRhoThermo
    mixture         multicomponentMixture;
    transport       twoSCH4Bfer;
    thermo          janaf;            // hConst / eConst / janaf
    energy          sensibleEnthalpy; // or sensibleInternalEnergy
    equationOfState perfectGas;
    specie          specie;
}

defaultSpecie N2;

#include "thermo.compressibleGas"
```

Species transport properties (inside `thermo.compressibleGas`):

```cpp
transport
{
    mu          1.8405e-5;   // mu0  [Pa.s]
    T0          300;         // T0   [K]
    muExponent  0.6759;      // alpha []
    Pr          0.7;         // OR kappa [W/m/K] (not both)
}
```

### 3. Define reactions

In `constant/reactions`:

```cpp
reactions
{
    R1
    {
        type        irreversibleArrheniusTwoSCH4Bfer;
        reaction    "CH4^0.5 + 1.5O2^0.65 = CO + 2H2O";
        A           1.8677e9;
        beta        0;
        Ta          17866;
    }
    R2
    {
        type        reversibleArrheniusTwoSCH4Bfer;
        reaction    "CO^1 + 0.5O2^0.5 = CO2";
        A           1.408e6;
        beta        0.8;
        Ta          6039;
    }
}
```

### 4. Provide the equivalence ratio

In `constant/phiEq`:

```cpp
value   0.8;
```

### 5. Lewis number = 1 (optional, default in org 14)

In `constant/thermophysicalTransport`:

```cpp
laminar
{
    model           unityLewisFourier;   // Le = 1: D = kappa/Cp
}
```

### 6. Run

```bash
foamRun
```

(org 14 uses the modular solver `foamRun` with the `multicomponentFluid` module;
`system/controlDict` must set `solver multicomponentFluid;`.)

> **Do not use spaces in the case path** — OpenFOAM aborts on
> `fileName::stripInvalid()`.

---

## Minimal example

A ready-to-adapt case is provided under:

```
OpenFOAM_run/OpenFOAM_14/phi-0.8-v2606-clear/
```

---

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `Unknown mixture type multicomponentMixture` | Library not loaded — check `libs (...)` in `controlDict`, and that the `.so` was built for org 14. |
| `undefined symbol: chemistryReductionMethod<...>::dictionaryConstructorTablePtr_` | Old build — rebuild (reduction methods are now included). |
| `Min sum of specie mass fractions ... < 0.999` | org 14 sums mass fractions from the **2nd** species; put a product (≈ 0) first in the `species` list. |
| `Cannot open .../constant/phiEq` | Create `constant/phiEq` with entry `value <phi>;`. |
| `fileName::stripInvalid()` | Case path contains a space — move it to a space-free path. |

---

## License

Ported from ESI/OpenCFD (com) 2506 to openfoam.org 14.
Reaction mechanism based on the two-stage model of Bui-Pham (1992), extended
with an equivalence-ratio-dependent rate correction.
