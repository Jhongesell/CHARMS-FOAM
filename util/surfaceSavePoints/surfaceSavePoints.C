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
    surfaceSavePoints

Description
    Dump the points of the surface to a XYZ file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "quaternion.H"
#include "mathematicalConstants.H"

#include "MeshedSurfaces.H"

using namespace Foam;
using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform (scale/rotate) a surface. "
        "Like transformPoints but for surfaces."
    );
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("output surfaceFile");

    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];

    Info<< "Reading surf from " << surfFileName << " ..." << nl
        << "Writing surf to " << outFileName << " ..." << endl;

    meshedSurface surf1(surfFileName);

    pointField points(surf1.points());

    OFstream os(outFileName);

    forAll(points, pointI)
    {
       os << points[pointI].x() << " " 
          << points[pointI].y() << " "
          << points[pointI].z() << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
