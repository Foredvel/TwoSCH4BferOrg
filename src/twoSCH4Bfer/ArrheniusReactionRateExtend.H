/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

Class
    Foam::ArrheniusReactionRateExtend

Description
    Arrhenius reaction rate with an equivalence-ratio dependent multiplier:

        k = A * T^beta * exp(-Ta/T) * f(phi)

    The multiplier f(phi) is selected by the reaction name ("R1" or "R2") and
    the equivalence ratio phi is read once from constant/phiEq.

    Ported from the ESI/OpenCFD (com) version to the openfoam.org (org)
    code base.

SourceFiles
    ArrheniusReactionRateExtendI.H

\*---------------------------------------------------------------------------*/

#ifndef ArrheniusReactionRateExtend_H
#define ArrheniusReactionRateExtend_H

#include "speciesTable.H"
#include "scalarField.H"
#include "typeInfo.H"
#include "IFstream.H"
#include "dictionary.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class ArrheniusReactionRateExtend;

Ostream& operator<<(Ostream&, const ArrheniusReactionRateExtend&);


/*---------------------------------------------------------------------------*\
               Class ArrheniusReactionRateExtend Declaration
\*---------------------------------------------------------------------------*/

class ArrheniusReactionRateExtend
{
    // Private Data

        scalar beta_;
        scalar A_;
        scalar Ta_;
        scalar phiEx_;
        word name_;


public:

    // Constructors

        //- Construct from components
        inline ArrheniusReactionRateExtend
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const word& name = word::null
        );

        //- Construct from dictionary
        inline ArrheniusReactionRateExtend
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );


    // Member Functions

        //- Read the equivalence ratio from constant/phiEq (cached)
        static scalar readPhi();

        //- Return the type name
        static word type()
        {
            return "ArrheniusTwoSCH4Bfer";
        }

        //- Pre-evaluation hook
        inline void preEvaluate() const;

        //- Post-evaluation hook
        inline void postEvaluate() const;

        //- Return the rate
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- The derivative of the rate w.r.t. temperature
        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Is the rate a function of concentration?
        inline bool hasDdc() const;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& ddc
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const ArrheniusReactionRateExtend&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ArrheniusReactionRateExtendI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
