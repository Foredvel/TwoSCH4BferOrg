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

\*---------------------------------------------------------------------------*/

#include "ArrheniusReactionRateExtend.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ArrheniusReactionRateExtend::ArrheniusReactionRateExtend
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const word& name
)
:
    beta_(beta),
    A_(A),
    Ta_(Ta),
    phiEx_(readPhi()),
    name_(name)
{}


inline Foam::ArrheniusReactionRateExtend::ArrheniusReactionRateExtend
(
    const speciesTable&,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    beta_(dict.lookup<scalar>("beta", dimless)),
    A_(dict.lookup<scalar>("A", dims/pow(dimTemperature, beta_))),
    Ta_
    (
        dict.found("Ta") || !dict.found("Ea")
      ? dict.lookup<scalar>("Ta", dimTemperature)
      : dict.lookup<scalar>("Ea", dimEnergy/dimMoles)
       /constant::physicoChemical::RR.value()
    ),
    phiEx_(readPhi()),
    name_(dict.dictName())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ArrheniusReactionRateExtend::readPhi()
{
    // Read once per run, then served from cache
    static const scalar phi = []()
    {
        const fileName file
        (
            fileName(getEnv("FOAM_CASE"))/"constant"/"phiEq"
        );

        IFstream is(file);

        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open " << file << nl
                << "Create it with entry: value <phi>;" << nl
                << exit(FatalError);
        }

        const dictionary d(is);
        return d.lookup<scalar>("value");
    }();

    return phi;
}


inline void Foam::ArrheniusReactionRateExtend::preEvaluate() const
{}


inline void Foam::ArrheniusReactionRateExtend::postEvaluate() const
{}


inline Foam::scalar Foam::ArrheniusReactionRateExtend::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/T);
    }

    // Equivalence-ratio dependent multiplier
    if (name_ == "R1")
    {
        const scalar p_0_1 = 1.1,  s_0_1 = 0.09, b_1 = 0.37;
        const scalar p_1_1 = 1.13, s_1_1 = 0.03, c_1 = 6.7;
        const scalar p_2_1 = 1.6,  s_2_1 = 0.22;

        scalar f1 =
            2.0
          / (
                1.0
              + tanh((p_0_1 - phiEx_)/s_0_1)
              + b_1*(1.0 + tanh((phiEx_ - p_1_1)/s_1_1))
              + c_1*(1.0 + tanh((phiEx_ - p_2_1)/s_2_1))
            );

        ak *= f1;
    }
    else if (name_ == "R2")
    {
        const scalar p_0_2 = 0.95, s_0_2 = 0.08, b_2 = 2.5e-5;
        const scalar p_1_2 = 1.3,  s_1_2 = 0.04, c_2 = 0.0087;
        const scalar p_2_2 = 1.2,  s_2_2 = 0.04;
        const scalar p_3_2 = 1.2,  s_3_2 = 0.05;

        scalar f2 =
              0.5*(1.0 + tanh((p_0_2 - phiEx_)/s_0_2))
            + 0.5*b_2*(1.0 + tanh((phiEx_ - p_1_2)/s_1_2))
            + 0.5*c_2
             *(1.0 + tanh((phiEx_ - p_2_2)/s_2_2))
             *(1.0 + tanh((p_3_2 - phiEx_)/s_3_2));

        ak *= f2;
    }

    return ak;
}


inline Foam::scalar Foam::ArrheniusReactionRateExtend::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/T);
    }

    return ak*(beta_ + Ta_/T)/T;
}


inline bool Foam::ArrheniusReactionRateExtend::hasDdc() const
{
    return false;
}


inline void Foam::ArrheniusReactionRateExtend::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    ddc = 0;
}


inline void Foam::ArrheniusReactionRateExtend::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ArrheniusReactionRateExtend& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
