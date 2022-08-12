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
    caclPressDerOnSlices

Description
    Calculates and prints average pressures over slices along some given
    coordinate. In user specified box
     

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sampledPlane.H"
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    // Info <<timeDirs.size()<<endl;
    runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);
    Info<< "Time = " << runTime.timeName() << endl;
    mesh.readUpdate();
    
    // -- read molar fraction field
    volScalarField CO 
    (
        IOobject
            (
            "CO", 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // -- read rho field
    volScalarField rho
    (
        IOobject
            (
            "rho", 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // -- read rho field
    volScalarField reactingCellZone
    (
        IOobject
            (
            "reactingCellZone", 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // -- read molar mass
    // -- thermophysicalProperties dictionary
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
        "thermophysicalProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    // -- molar mass
    scalar molMR(readScalar(thermophysicalProperties.subDict("mixture").subDict("specie").lookup("molWeight")));
    
    // -- integration (sum c_s*V in cellZone)
    scalar integral(0);

    // // -- loop over all cells in mesh    
    forAll(mesh.cells(), celli)
    {

    //     // integral += mesh.V()[celli]*CO[celli]*rho[celli]/molMR;
    //     // integral += mesh.V()[celli]*CO[celli]*rho[celli]/molMR;
        integral += mesh.V()[celli]*reactingCellZone[celli]*CO[celli]*rho[celli]/(molMR*1e-3);
        // Info<<reactingCellZone[celli]<<endl;
    }

    Info << "Integral reaction source = " << integral << endl;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
