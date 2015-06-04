/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    calUwater 

Description
    Calculates and writes the velocity field U of water phase only.

    This tools is designed for interFoam alike solvers. In some applications,
    we are only interested in the velocity of water phase.

    The -noWrite option just outputs the max/min values without writing
    the field.


\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject alpha1Header
    (
        "alpha.water",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk() && alpha1Header.headerOk())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);

        Info<< "    Reading alpha1" << endl;
        volScalarField alpha1(alpha1Header, mesh);

        Info<< "    Calculating Uwater" << endl;
        volVectorField Uwater
        (
            IOobject
            (
                "Uwater",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            pos(alpha1-0.5)*U
        );

        volScalarField magUwater
        (
            IOobject
            (
                "magUwater",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mag(Uwater)
        );

        Info<< "Uwater max/min : "
            << max(magUwater).value() << " "
            << min(magUwater).value() << endl;

        if (writeResults)
        {
            Uwater.write();
        }
    }
    else
    {
        Info<< "    No U or No alpha.water" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
