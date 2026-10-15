/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Macro cascades which instantiate the twoSCH4Bfer transport package for the
    full set of gas/liquid equation-of-state + thermo combinations provided by
    the ESI/OpenCFD (com) library:

        - janaf  + perfectGas / PengRobinsonGas / incompressiblePerfectGas
        - hConst + perfectGas / incompressiblePerfectGas / perfectFluid /
                   rPolynomial / adiabaticPerfectFluid / rhoConst
        - eConst + perfectGas / incompressiblePerfectGas / perfectFluid /
                   adiabaticPerfectFluid / rhoConst
        - hPolynomialThermo<icoPolynomial<specie,8>,8>
          (sensibleEnthalpy + sensibleInternalEnergy)

    OpenFOAM 14 marks perfectFluid, adiabaticPerfectFluid and icoPolynomial as
    `noCoefficientMixing`, so they must not be used with
    coefficientMulticomponentMixture.  Two cascades are therefore provided:

        forTwoSCH4BferCoeff  - coefficient-mixing compatible EoS only.
                               Use with coefficientMulticomponentMixture,
                               reactions and chemistry models.

        forTwoSCH4Bfer       - the full set (coefficient + single-component
                               EoS). Use with singleComponentMixture.

    The icoPolynomial / hPolynomialThermo combinations use two-parameter
    templates, which the org `forThermo` macro cannot express, so they are
    declared separately by `defineTwoSCH4BferIcoPoly8Types` (call it exactly
    once per translation unit before `forTwoSCH4Bfer`).

\*---------------------------------------------------------------------------*/

#ifndef forTwoSCH4Bfer_H
#define forTwoSCH4Bfer_H

#include "specie.H"

#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "incompressiblePerfectGas.H"
#include "perfectFluid.H"
#include "rPolynomial.H"
#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "icoPolynomial.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "hPolynomialThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "thermo.H"

#include "forThermo.H"

#include "twoSCH4Bfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Declare the two icoPolynomial-based thermo-physics typedefs.  This must be
// invoked exactly once per translation unit (the generated typedef names are
// fixed, unlike the names produced by forThermo which are unique per
// combination).
#define defineTwoSCH4BferIcoPoly8Types                                         \
                                                                               \
    typedef                                                                    \
        twoSCH4Bfer                                                            \
        <                                                                      \
            species::thermo                                                    \
            <                                                                  \
                hPolynomialThermo                                              \
                <                                                              \
                    icoPolynomial<specie, 8>,                                  \
                    8                                                          \
                >,                                                             \
                sensibleEnthalpy                                               \
            >                                                                  \
        >                                                                      \
        twoSCH4BfericoPoly8HThermoPhysics;                                     \
                                                                               \
    typedef                                                                    \
        twoSCH4Bfer                                                            \
        <                                                                      \
            species::thermo                                                    \
            <                                                                  \
                hPolynomialThermo                                              \
                <                                                              \
                    icoPolynomial<specie, 8>,                                  \
                    8                                                          \
                >,                                                             \
                sensibleInternalEnergy                                         \
            >                                                                  \
        >                                                                      \
        twoSCH4BfericoPoly8EThermoPhysics


// Coefficient-mixing compatible equation-of-state combinations (13 types).
// These may be wrapped in coefficientMulticomponentMixture.
#define forTwoSCH4BferCoeff(Macro, Args...)                                    \
                                                                               \
    /* janaf + sensibleEnthalpy */                                             \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, janafThermo,                            \
        perfectGas, specie, Macro, Args                                        \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, janafThermo,                            \
        PengRobinsonGas, specie, Macro, Args                                   \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, janafThermo,                            \
        incompressiblePerfectGas, specie, Macro, Args                          \
    );                                                                         \
                                                                               \
    /* janaf + sensibleInternalEnergy */                                       \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, janafThermo,                      \
        perfectGas, specie, Macro, Args                                        \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, janafThermo,                      \
        PengRobinsonGas, specie, Macro, Args                                   \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, janafThermo,                      \
        incompressiblePerfectGas, specie, Macro, Args                          \
    );                                                                         \
                                                                               \
    /* hConst + sensibleEnthalpy */                                            \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, hConstThermo,                           \
        perfectGas, specie, Macro, Args                                        \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, hConstThermo,                           \
        incompressiblePerfectGas, specie, Macro, Args                          \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, hConstThermo,                           \
        rPolynomial, specie, Macro, Args                                       \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, hConstThermo,                           \
        rhoConst, specie, Macro, Args                                          \
    );                                                                         \
                                                                               \
    /* eConst + sensibleInternalEnergy */                                      \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, eConstThermo,                     \
        perfectGas, specie, Macro, Args                                        \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, eConstThermo,                     \
        incompressiblePerfectGas, specie, Macro, Args                          \
    );                                                                         \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, eConstThermo,                     \
        rhoConst, specie, Macro, Args                                          \
    )


// Single-component only (noCoefficientMixing) equation-of-state combinations
// (6 types). Do NOT wrap these in coefficientMulticomponentMixture.
#define forTwoSCH4BferSingle(Macro, Args...)                                   \
                                                                               \
    /* hConst + perfectFluid (sensibleEnthalpy) */                             \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, hConstThermo,                           \
        perfectFluid, specie, Macro, Args                                      \
    );                                                                         \
                                                                               \
    /* eConst + perfectFluid (sensibleInternalEnergy) */                       \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, eConstThermo,                     \
        perfectFluid, specie, Macro, Args                                      \
    );                                                                         \
                                                                               \
    /* hConst + adiabaticPerfectFluid (sensibleEnthalpy) */                    \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleEnthalpy, hConstThermo,                           \
        adiabaticPerfectFluid, specie, Macro, Args                             \
    );                                                                         \
                                                                               \
    /* eConst + adiabaticPerfectFluid (sensibleInternalEnergy) */              \
    forThermo                                                                  \
    (                                                                          \
        twoSCH4Bfer, sensibleInternalEnergy, eConstThermo,                     \
        adiabaticPerfectFluid, specie, Macro, Args                             \
    );                                                                         \
                                                                               \
    /* hPolynomialThermo<icoPolynomial<specie,8>,8> (two-parameter templates) */ \
    Macro(Args, twoSCH4BfericoPoly8HThermoPhysics);                            \
    Macro(Args, twoSCH4BfericoPoly8EThermoPhysics)


// Full set (coefficient + single-component, 19 types). Use with
// singleComponentMixture only.
#define forTwoSCH4Bfer(Macro, Args...)                                         \
    forTwoSCH4BferCoeff(Macro, Args);                                          \
    forTwoSCH4BferSingle(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
