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
    Registers the ArrheniusReactionRateExtend reaction rate for all
    twoSCH4Bfer thermo types, so that reactions can use:

        reactions
        {
            R1
            {
                type    irreversibleArrheniusExtend;
                reaction "CH4 + ...";
                A       ...;
                beta    ...;
                Ta      ...;
            }
        }

\*---------------------------------------------------------------------------*/

#include "forTwoSCH4Bfer.H"
#include "makeReaction.H"

#include "ArrheniusReactionRateExtend.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forTwoSCH4BferCoeff(defineReaction, nullArg);
    forTwoSCH4BferCoeff(makeIRNReactions, ArrheniusReactionRateExtend);
}

// ************************************************************************* //
