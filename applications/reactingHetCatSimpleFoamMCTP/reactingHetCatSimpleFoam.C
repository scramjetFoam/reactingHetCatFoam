/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    diffReactConvSysThreeWayMacro

Description
    Steady-state non-isothermal solver for heterogeneously catalyzed reactive fluid flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermophysicalTransportModel.H"
#include "simpleControl.H"
// #include "pressureControl.H"
#include "pressureReference.H"
// #include "fvOptions.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "IOporosityModelList.H"
// #include "DEff.H"
// #include "Switch.H"
// #include "MaxwellStefanMy/MaxwellStefan.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
   
    // #include "setRootCase.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "readTransportPropertiesMy.H"
    #include "postProcess.H"
    #include "createFields.H"
    //~ #include "createPorousZones.H"
    #include "createZones.H"
    // #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl; 
        
        {
            scalar nItconc(readScalar(simple.dict().lookup("nConcCorrectors")));
            for (label concIt = 0; concIt < nItconc; concIt++)
            { 
                Info << "\nSolving continuity equation for each specie, iteration "<< concIt+1 << "/" << nItconc << endl;
                // #include "concEq.H"
                #include "concEq.H"
                // #include "concEqMass.H"
            }

            Info << "\nSolving momentum equation"<<endl;
             
            #include "UEqn.H"

            scalar nItTemp(readScalar(simple.dict().lookup("nTempCorrectors")));

            for (label concT=0; concT < nItTemp; concT++)
            {   
                Info << "\nSolving enthalpy balance, iteration "<< concT+1 << "/" << nItTemp << endl;
                #include "EEqn.H"
            }
            
            Info << "\nCorrection of pressure"<<endl;
            #include "pEqn.H"   
            // #include "rhoEqn.H"   
        }              
                
        turbulence->correct();
        thermophysicalTransport->correct();

        runTime.write();
        
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
