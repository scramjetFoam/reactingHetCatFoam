Info<< "Reading the concentrations \n" << endl;

// prepare list to store the concentrations in
PtrList<volScalarField> MFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    Info<< "Reading the concentration of the: " << currentSp << endl;

    MFrList.set
    (
        countI,
        new volScalarField
        (
            IOobject
            (
                currentSp,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}


// prepare list to store volScalarField  with the zone inside
// NOTE: 1 - zone is active here
//       0 - zone is not active here

PtrList<volScalarField> zoneLst(zoneNames.size());
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    Info<< "Reading the zone " << currentZone << " and saving it as volScalarField." << endl;

    // -- load cellZone as 0 volScalarField
    zoneLst.set
    (
        countI,
        new volScalarField
        (
            IOobject
            (
                currentZone,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar (currentZone,dimless,0.0),
            "zeroGradient"
        )
    );
}

// -- set to 1, where the cellZone is
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];

    labelList cellZoneIDs_ = mesh.cellZones().findIndices(currentZone);

    const labelList& cells = mesh.cellZones()[cellZoneIDs_[0]];
    forAll(cells, celli)
    {
        zoneLst[countI][cells[celli]] = 1;
    }
}

// -- free channels --> free channels = 1 - sum all zones 
volScalarField sumZones
(
    IOobject
    (
        "sumZones",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sumZones",dimless,0.0)
); 

forAll(zoneLst,countI)
{
    sumZones += zoneLst[countI];
}
volScalarField freeChannel = dimensionedScalar("jedna",dimless,1.0)-sumZones;

// -- porous media effective thermal diffusivity
volScalarField kappaEff
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaEff",kappaA.dimensions(),0.0)
);

// -- porous media heat capacity 
volScalarField CpEff
(
    IOobject
    (
        "CpEff",
        runTime.timeName(),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("CpEff",dimEnergy/dimMass/dimTemperature,0.0)
);

// -- porous media efective mass density 
volScalarField rhoEff
(
    IOobject
    (
        "rhoEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoEff",dimMass/dimVolume,0.0)
);

// -- heat transfer coefficient (W/m2/K)
volScalarField kHEff
(
    IOobject
    (
        "kHEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kHEff",(dimEnergy/dimTime/dimLength/dimLength/dimTemperature),0.0)
);

// -- porosity
volScalarField epsEff
(
    IOobject
    (
        "epsEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("epsEff",dimless,0.0)
);

// -- specific area
volScalarField aZoneEff
(
    IOobject
    (
        "aZoneEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("aZoneEff",dimless/dimLength,0.0)
);

// -- specific area corrected to neigbouring cells
// -- NOTE: prepared for implementation of the interface heat source
// volScalarField aEff
// (
//     IOobject
//     (
//         "aEff",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("aEff",dimless/dimLength,0.0)
// );

// -- efective kappa values, gas kappa in free channel, media kappa in zones     
forAll(zoneLst,countI)
{
    epsEff += zoneLst[countI]*(eps[countI]);
    // kappaEff += zoneLst[countI]*(kappaZ[countI]);
    kappaEff += zoneLst[countI]*(kappaZ[countI]);
    CpEff += zoneLst[countI]*(CpZ[countI]);
    rhoEff += zoneLst[countI]*(rhoZ[countI]);
    kHEff += zoneLst[countI]*(kH[countI]);
    aZoneEff += zoneLst[countI]*(aZone[countI]);
} 
epsEff += freeChannel;
// volScalarField kappaS = kappaEff;

// -- find cells on the surface 
// prepared for implementation of the interface heat source
// counter (number of neigbouring cells with different )
// forAll(mesh.cells(),celli)
// {
//     // forAll(mesh.cells()[celli].faceNeighbour(),facei)
//     // {
//     if(aZoneEff[celli] > 0)
//     {
//         aEff[celli] += aZoneEff[celli];
//     }
//     else
//     {
//         // Info<<mesh.cellCells()[celli]<<endl;
//         forAll(mesh.cellCells()[celli],neigbi)
//         {
//             int neighIndex = mesh.cellCells()[celli][neigbi];
//             if (aZoneEff[neighIndex] > 0)
//             {
//                 aEff[celli] += aZoneEff
//             }
//         }
//     }
// }
// aEff.write();

// -- heat transfer coefficient (W/m3/K)
volScalarField kHEffA = kHEff*aZoneEff;

rhoEff.correctBoundaryConditions();
kappaEff.correctBoundaryConditions();
CpEff.correctBoundaryConditions();
// kHEff.correctBoundaryConditions();
// aZoneEff.correctBoundaryConditions();
// kHEffA.correctBoundaryConditions();


// -- reading field U
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
	"U",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

// -- dummy volScalarField to resolve dimension problem in conc.H fvOptions inconsistent source
volScalarField dummyConcEqn
(
    IOobject
    (
        "dummyConcEqn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("dummyConcEqn",dimMoles/dimVolume,1.0)
    // dimensionedScalar ("dummyConcEqn",dimless,1.0)
);

// -- dummy volScalarField to resolve dimension problem in conc.H fvOptions inconsistent source
//NOTE: I dont know why this is necessary, dimensions are correct
//but openFOAM doesnt think so :(
volScalarField dummyEEqn
(
    IOobject
    (
        "dummyEEqn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("dummyEEqn",dimMass/dimLength*Foam::pow(dimTime,-2)/dimTemperature,1.0)
);

// -- T field from mesh
// Info<< "Reading field T\n" << endl;
// volScalarField T
// (
//     IOobject
//     (
//         "T",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );

// -- T field for porousMedia from mesh
Info<< "Reading field TSolid\n" << endl;
volScalarField TSolid
(
    IOobject
    (
        "TSolid",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);


fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);


#include "compressibleCreatePhi.H"

pressureReference pressureReference
(
    p,
    simple.dict(),
    thermo.incompressible()
);

mesh.schemes().setFluxRequired(p.name());


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

//-- surface field for kappaG interpolation
surfaceScalarField MGradT
(
    IOobject
    (
        "MGradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::snGrad(epsEff)
);


Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"




