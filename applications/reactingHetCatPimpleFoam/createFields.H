// -- dynamical stuff
#include "createRDeltaT.H"

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
            dimensionedScalar (currentZone,dimless,0.0)
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
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sumZones",dimless,0.0)
); 

forAll(zoneLst,countI)
{
    sumZones += zoneLst[countI];
}
volScalarField freeChannel = dimensionedScalar("jedna",dimless,1.0)-sumZones;


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

// -- T field for porousMedia from mesh
Info<< "Reading field TporMedia\n" << endl;
volScalarField TporMedia
(
    IOobject
    (
        "TporMedia",
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
    // p/univR/T*molMList[molMList.size()-1]
);


#include "compressibleCreatePhi.H"

pressureReference pressureReference
(
    p,
    pimple.dict(),
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


Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);
Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"




