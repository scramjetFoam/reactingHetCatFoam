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
#include "argList.H"
// #include "immiscibleIncompressibleTwoPhaseMixture.H"

#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    
    if (argc < 2)                                                       //I need a fieldName - where to correct the BC
    {
        FatalError
            << "Not enough arguments" << nl
            << exit(FatalError);
    }
    
    const label nCuts = strtod(argv[1],NULL);                          //get the number of cuts to make
    argList::addOption
    (
        "cutDir",
        "vector",
        "create cutting planes normal to the given direction <vector> - eg, '(1 0 0)'"
    );
    
    argList::addOption
    (
        "boxMin",
        "vector",
        "specify bottom left corner of the box of interest <vector> - eg, '(0 0 0)'"
    );
    
    argList::addOption
    (
        "boxMax",
        "vector",
        "specify top right corner of the box of interest <vector> - eg, '(1 1 1)'"
    );
    
    argList::validArgs.append("nCuts");
    
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    
    vector cutDir(0,0,1);
    if (args.optionReadIfPresent("cutDir", cutDir))
    {
        Info<< "Cutting plane is normal to: " << cutDir << endl;
    }
    else
    {
        Info<< "Default cutting plane is normal to: " << cutDir << endl;
    }
    
    vector boxMin(0,0,0);
    if (args.optionReadIfPresent("boxMin", boxMin))
    {
        Info<< "Defined box bottom left corner: " << boxMin << endl;
    }
    else
    {
        Info<< "Default box bottom left corner: " << boxMin << endl;
    }
    
    vector boxMax(1,1,1);
    if (args.optionReadIfPresent("boxMax", boxMax))
    {
        Info<< "Defined box top right corner: " << boxMax << endl;
    }
    else
    {
        Info<< "Default box top right corner: " << boxMax << endl;
    }
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();
        
        // read the velocity field
	    volVectorField U 
        (
            IOobject
            (
                "U", 
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField coat
        (
            IOobject
            (
            //	"coatWallFrac", 
                // "coat", 
                "reactionZone", 
                runTime.timeName(0),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
            
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
        
        // prepare the cutting surfaces
        vector dirMin(0,0,0);
        vector dirMax(1,1,1);
        scalar cutCmpt(0);
        vector planeCenter(0,0,0);                                      //prepare the vector for the plane center
        
        // find the bouding box of the mesh (in the given direction)
        for (direction i=0; i<vector::nComponents; i++)
        {
            dirMin[i] = min(mesh.C().component(i)).value();
            dirMax[i] = max(mesh.C().component(i)).value();
            
            if (cutDir[i]!=0)
            {
                cutCmpt = i;
                planeCenter[i]=dirMin[i];
            }
            else
            {
                planeCenter[i]=(dirMax[i] + dirMin[i])/2.0;
            }
            
        }
        
        Info<< "Mesh bounding box: min: " << dirMin << " max: " << dirMax << endl;
        
        scalar stepSize((dirMax[cutCmpt]-dirMin[cutCmpt])/(nCuts+SMALL));
        vector currPos(planeCenter);
        
	    // scalar koefArea(0);
        Info<< "Creating flux averages over the selected planes\n" << endl;
    
	    Info << "z\tkonvCO\tCO\tarea" << endl;
        for (direction i=0; i<nCuts; i++)
        {
            currPos = planeCenter + i*stepSize*cutDir;                  //update the current position
            point planeLocation(currPos);                               //prepare the point with current position
            plane smpPl1(planeLocation,cutDir);                       //auxiliary plane
            word nazev= "smpl";
            Foam::sampledSurfaces::plane smpPl(nazev,mesh,smpPl1);                       //auxiliary plane
            // sampledPlane smpPl("smpPl",mesh,extPlane);                  //this might not work
            // ("smpPl",mesh,extPlane);                  //this might not work
            
            smpPl.update();
            
            Info << smpPl.faces() << endl;
            //~ Info << smpPl.meshCells() << endl;
            
	        scalar mdot(0);
	        scalar mdotCO(0);
            scalar areaWeig(0);
	        scalar nCutCells(0);
	        // scalar areaReal((boxMax[0]-boxMin[0])*(boxMax[1]-boxMin[1]));
            
//             forAll(smpPl.meshCells(),smpCellI)
//             {
//                 if (mesh.C()[smpPl.meshCells()[smpCellI]].component(cutCmpt) < boxMin.component(cutCmpt) or mesh.C()[smpPl.meshCells()[smpCellI]].component(cutCmpt) > boxMax.component(cutCmpt))
//                 {
//                     // avPressBox = -1.0;
//                     nCutCells  = SMALL;
//                     Info << "!! Skipped Plane !!" << endl;
//                     break;
//                 }
//                 bool includeCell(true);
//                 for (direction j=0; j<vector::nComponents; j++)
//                 {
//                     if (mesh.C()[smpPl.meshCells()[smpCellI]].component(j) < boxMin.component(j) or mesh.C()[smpPl.meshCells()[smpCellI]].component(j) > boxMax.component(j))
//                     {
//                         //Info << "cell at loc" << mesh.C()[smpPl.meshCells()[smpCellI]] << "not included"<<endl;
// 			            includeCell = false;
//                         break;
//                     }
//                 }
//                 if (includeCell)
//                 {
//                     // if (field == "U"){
//                     if (coat[smpPl.meshCells()[smpCellI]] == 0)
//                     {
//                         scalar SCell((Foam::pow(mesh.V()[smpPl.meshCells()[smpCellI]],2.0/3.0)));
//                         mdot += U[smpPl.meshCells()[smpCellI]][2]*SCell;
//                         mdotCO += U[smpPl.meshCells()[smpCellI]][2]*SCell*CO[smpPl.meshCells()[smpCellI]];
//                         areaWeig += SCell;
//                         nCutCells  += 1.0;
//                     }
// 			        // }
// 		            // else{
//                     // scalar SCell((Foam::pow(mesh.V()[smpPl.meshCells()[smpCellI]],2.0/3.0)));
//                     // areaWeig += SCell;
//                     // nCutCells  += 1.0;
//                     // }
//                     //scalar SCell(Foam::pow(mesh.V()[smpPl.meshCells()[smpCellI]],2.0/3.0));
//                     //avPressBoxWeig += p[smpPl.meshCells()[smpCellI]]*SCell;
//                     //avPressBoxWeig += avPressBoxWeig*SCell;
                    
//                     //~ Info << mesh.C()[smpPl.meshCells()[smpCellI]] << endl;
//                 }
//             }                                    
//             //Info << "Plane #     : " << i+1 << endl;
//             //Info << "Plane pos.  : " << currPos <<endl;
// //            Info << currPos <<"\t"<< avPressBoxWeig/areaWeig<<endl;
//             //Info << "avgPressure : " << avPressBox/(nCutCells) << " m2s-2" << endl;
//             if (areaWeig == 0)
//             {
//             	Info << "Zero area" << endl;
//             }
//             else
//             {
//                 // if (koefArea == 0){
//                 // koefArea = areaReal/areaWeig;
//                 Info << currPos[2] <<"\t"<< (0.002 - mdotCO/mdot)/0.002*100<<"\t"<<mdotCO/mdot<<"\t"<<areaWeig <<endl;
//                 // koefArea = 0;
//                 // }
//                 // else{
//                 //     Info << currPos[2] <<"\t"<< avPressBoxWeig/areaWeig<<"\t"<<avTBoxWeig/areaWeig<<"\t"<< avYCO/areaWeig <<"\t"<< Vdot*koefArea<<"\t"<<VdotCO*koefArea<<"\t"<<koefArea <<endl;
//                 // }
// 		    }

		}
    }
        
    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
