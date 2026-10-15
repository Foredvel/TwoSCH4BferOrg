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
    Instantiates the chemistry reduction-method runtime-selection tables for
    all twoSCH4Bfer thermo types, and registers the available reduction
    methods (none, DRG, DRGEP, DAC, EFA, PFA).

    This is required because the org 14 `Standard` chemistry model (used by
    `makeChemistryModelsTwoS.C`) is templated on the thermo type and carries a
    `chemistryReductionMethod<ThermoType>` runtime-selection table which would
    otherwise remain an undefined symbol and prevent the library from loading.

\*---------------------------------------------------------------------------*/

#include "chemistryReductionMethod.H"

#include "noChemistryReduction.H"
#include "DRG.H"
#include "DRGEP.H"
#include "DAC.H"
#include "EFA.H"
#include "PFA.H"

#include "forTwoSCH4Bfer.H"
#include "makeChemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Define the reduction-method runtime-selection tables
    forTwoSCH4BferCoeff(defineChemistryReductionMethod, nullArg);

    // Register the available reduction methods
    forTwoSCH4BferCoeff(makeChemistryReductionMethod, none);
    forTwoSCH4BferCoeff(makeChemistryReductionMethod, DRG);
    forTwoSCH4BferCoeff(makeChemistryReductionMethod, DRGEP);
    forTwoSCH4BferCoeff(makeChemistryReductionMethod, DAC);
    forTwoSCH4BferCoeff(makeChemistryReductionMethod, EFA);
    forTwoSCH4BferCoeff(makeChemistryReductionMethod, PFA);
}

// ************************************************************************* //
