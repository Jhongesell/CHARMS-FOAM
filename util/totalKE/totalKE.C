/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    totalKE 

Description
    Calculates and writes the scalar field of total kinetic energy at
    each time.

    In this case, totalKE = \int_V 0.5(U&U) dV. It does not consider the 
    models (RANS or LES). So be careful in interpreting the meaning of totalKE.

Author
    Xiaofeng Liu, Ph.D., P.E.
    Penn State University

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    OFstream totalKE_File
    (
        "totalKE.dat"
    );

    runTime.setTime(timeDirs.last(), timeDirs.size()-1);

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            Info<< "    Calculating total kinetic energy" << endl;

            dimensionedScalar totalK = sum(0.5*(U&U)*mesh.V()); 

            Info<< "    total kinetic energy = " << totalK.value() << endl;

            //save to file
            if(Pstream::parRun())
            {
                if (Pstream::myProcNo() == Pstream::masterNo())
                {
                     totalKE_File << runTime.timeName() 
                                  << " "  << totalK.value() << endl;
                }
            }
            else
            {
                totalKE_File << runTime.timeName() 
                             << " "  << totalK.value() << endl;
            }
        }
        else
        {
            Info<< "    No U" << endl;
        }

        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
