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
    Instantiates the twoSCH4Bfer transport package as a selectable
    multicomponent thermo (both psi- and rho-based), so that it can be used
    with:

        thermoType
        {
            type            heRhoThermo;      // or hePsiThermo
            mixture         multicomponentMixture;
            transport       twoSCH4Bfer;
            thermo          janaf;            // hConst / eConst / janaf
            energy          sensibleEnthalpy; // or sensibleInternalEnergy
            equationOfState perfectGas;
            specie          specie;
        }

\*---------------------------------------------------------------------------*/

#include "psiMulticomponentThermo.H"
#include "rhoFluidMulticomponentThermo.H"

#include "coefficientMulticomponentMixture.H"
#include "singleComponentMixture.H"

#include "forTwoSCH4Bfer.H"

#include "makeFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTwoSCH4BferIcoPoly8Types;

    // psi-based (compressible, p-based)
    forTwoSCH4BferCoeff
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        psiMulticomponentThermo,
        coefficientMulticomponentMixture
    );
    forTwoSCH4Bfer
    (
        makeFluidMulticomponentThermo,
        psiMulticomponentThermo,
        singleComponentMixture
    );

    // rho-based (compressible, rho-based)
    forTwoSCH4BferCoeff
    (
        makeFluidMulticomponentThermos,
        rhoFluidThermo,
        rhoFluidMulticomponentThermo,
        coefficientMulticomponentMixture
    );
    forTwoSCH4Bfer
    (
        makeFluidMulticomponentThermo,
        rhoFluidMulticomponentThermo,
        singleComponentMixture
    );
}

// ************************************************************************* //
