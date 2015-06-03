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
    R_anisotropy

Description
    Calculates and writes the Reynolds stress R and then calculate
    the anisotropy  for the current time step.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

#include "fluidThermo.H"
#include "compressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressibleR
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> model
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    Info<< "Getting R and k fields" << nl << endl;

    volSymmTensorField& R = model->R()();

    volScalarField& k = model->k()();

    //normalize R field
    volSymmTensorField B = R/(2.0*k) - 1.0/3.0*I;

    //get the eigen values of the normalized B field
    volVectorField B_eigenValues
    (
        IOobject
        (
            "eigenValues",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        eigenValues(B)
    );

    //calculate the second and third invariants (the first = 0)
    volScalarField IIs
    (
        IOobject
        (
            "IIs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        B_eigenValues.component(vector::X)*
        B_eigenValues.component(vector::Y)+
        B_eigenValues.component(vector::Y)*
        B_eigenValues.component(vector::Z)+
        B_eigenValues.component(vector::X)*
        B_eigenValues.component(vector::Z)
    );
    
    volScalarField IIIs
    (
        IOobject
        (
            "IIIs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        B_eigenValues.component(vector::X)*
        B_eigenValues.component(vector::Y)*
        B_eigenValues.component(vector::Z)
    );

    Info<< "    Writing IIs and IIIs" << endl;
    IIs.write();
    IIIs.write();
}



int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "compressible",
        "calculate compressible R"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const bool compressible = args.optionFound("compressible");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject RHeader
        (
            "R",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOobject kHeader
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (RHeader.headerOk() && kHeader.headerOk())
        {
            Info<< "Reading R and k fields " << nl << endl;
            volSymmTensorField R(RHeader, mesh);
            volScalarField k(kHeader, mesh);

            //normalize R field
            volSymmTensorField B = R/(2.0*k) - 1.0/3.0*I;

            //get the eigen values of the normalized B field
            volVectorField B_eigenValues
            (
               IOobject
               (
                  "eigenValues",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
               ),
               eigenValues(B)
            );

            //calculate the second and third invariants (the first = 0)
            volScalarField IIs
            (
               IOobject
               (
                  "IIs",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
               ),
               B_eigenValues.component(vector::X)*
               B_eigenValues.component(vector::Y)+
               B_eigenValues.component(vector::Y)*
               B_eigenValues.component(vector::Z)+
               B_eigenValues.component(vector::X)*
               B_eigenValues.component(vector::Z)
            );
    
            volScalarField IIIs
            (
               IOobject
               (
                   "IIIs",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               B_eigenValues.component(vector::X)*
               B_eigenValues.component(vector::Y)*
               B_eigenValues.component(vector::Z)
            );

            IIs.write();
            IIIs.write();

            //calculate eta and xi fields
            volScalarField eta
            (
               IOobject
               (
                   "eta",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               sqrt(-IIs/3.0)
            );

            volScalarField xi
            (
               IOobject
               (
                   "xi",
                   runTime.timeName(),
                   mesh,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               sign(IIIs)*pow(mag(IIIs)/2.0, 1.0/3.0)
            );
         
            Info<< "    Writing eta and xi" << endl;
            eta.write();
            xi.write();
        }
        else
        {
            Info<< "    no U or k field" << endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
