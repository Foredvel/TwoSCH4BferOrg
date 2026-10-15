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
    Foam::twoSCH4Bfer

Description
    Transport package with power-law dynamic viscosity:

        mu(T) = mu0*(T/T0)^alpha

    and either a constant Prandtl number or a constant thermal conductivity.
    Ported from the ESI/OpenCFD (com) version to the openfoam.org (org)
    code base.

    Usage
    \table
        Property    | Description
        mu          | Reference dynamic viscosity [Pa.s] at T0
        T0          | Reference temperature [K]
        muExponent  | Power-law exponent alpha []
        Pr          | Constant Prandtl number []  (either Pr or kappa)
        kappa       | Constant thermal conductivity [W/m/K] (either Pr or kappa)
    \endtable

SourceFiles
    twoSCH4BferI.H
    twoSCH4Bfer.C

\*---------------------------------------------------------------------------*/

#ifndef twoSCH4Bfer_H
#define twoSCH4Bfer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class twoSCH4Bfer;

template<class Thermo>
inline twoSCH4Bfer<Thermo> operator+
(
    const twoSCH4Bfer<Thermo>&,
    const twoSCH4Bfer<Thermo>&
);

template<class Thermo>
inline twoSCH4Bfer<Thermo> operator*
(
    const scalar,
    const twoSCH4Bfer<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const twoSCH4Bfer<Thermo>&
);


/*---------------------------------------------------------------------------*\
                     Class twoSCH4Bfer Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class twoSCH4Bfer
:
    public Thermo
{
    // Private Data

        //- Reference dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reference temperature [K]
        scalar T0_;

        //- Power-law exponent []
        scalar alpha_;

        //- Is the Prandtl Number constant?
        bool constPr_;

        //- Reciprocal constant Prandtl Number []
        scalar rPr_;

        //- Constant thermal conductivity [W/m/K]
        scalar kappa_;


public:

    // Constructors

        //- Construct from components
        inline twoSCH4Bfer
        (
            const Thermo& t,
            const scalar mu,
            const scalar T0,
            const bool constPr,
            const scalar alpha,
            const scalar rPr,
            const scalar kappa
        );

        //- Construct as named copy
        inline twoSCH4Bfer(const word&, const twoSCH4Bfer&);

        //- Construct from name and dictionary
        twoSCH4Bfer(const word& name, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<twoSCH4Bfer> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "twoSCH4Bfer<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const twoSCH4Bfer&);
        inline void operator*=(const scalar);


    // Friend operators

        friend twoSCH4Bfer operator+ <Thermo>
        (
            const twoSCH4Bfer&,
            const twoSCH4Bfer&
        );

        friend twoSCH4Bfer operator* <Thermo>
        (
            const scalar,
            const twoSCH4Bfer&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const twoSCH4Bfer&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoSCH4BferI.H"

#ifdef NoRepository
    #include "twoSCH4Bfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
