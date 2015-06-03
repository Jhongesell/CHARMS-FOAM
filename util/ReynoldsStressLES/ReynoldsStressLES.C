/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    ReynoldsStressLES

Description
    Calculates and writes the full Reynolds stress R for the current time step.

    It is only for LES simulations. 

    For RANS simulutionls, use the generic R tool.

Author
    Xiaofeng Liu, Ph.D., P.E.
    Penn State University

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject UPrime2MeanHeader
        (
            "UPrime2Mean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOobject RMeanHeader
        (
            "RMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );


        if (UPrime2MeanHeader.headerOk() && RMeanHeader.headerOk())
        {
            Info<< "Reading field " << UPrime2MeanHeader.name() << nl << endl;
            volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);

            Info<< "Reading field " << RMeanHeader.name() << nl << endl;
            volSymmTensorField RMean(RMeanHeader, mesh);

            volSymmTensorField ReynoldsStress
            (
                IOobject
                (
                    "ReynoldsStress",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                UPrime2Mean + RMean
            );

            Info<< "    Writing Reynolds stresses" << endl;
            ReynoldsStress.write();
        }
        else
        {
            Info<< "    no " << RMeanHeader.name() << " field" << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
