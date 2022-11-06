// -- dynamical stuff
#include "createRDeltaT.H"

Info<< "Reading the concentrations \n" << endl;

// -- creating field Mg -- molar mass of gas
Info<< "Reading field Mg\n" << endl;
volScalarField Mg
(
    IOobject
    (
        "Mg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Mg",dimMass/dimMoles,0.0),
    "zeroGradient"
);

// prepare list to store the concentrations in
PtrList<volScalarField> MFrList(specieNames.size());
PtrList<volVectorField> jIList(specieNames.size());

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
    Mg += molMList[countI]*MFrList[countI];

    // -- creating velocity enforced by diffusion for specie
    word jiDiffcurrentSp = "jiDiff" + specieNames[countI];
    jIList.set
    (
        countI,
        new volVectorField
        (
            IOobject
            (
                jiDiffcurrentSp,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector ("jiDiff",dimLength/dimTime*dimMoles/dimVolume, vector(0.0,0.0,0.0))
            // dimensionedVector ("jiDiff",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
        )
    );
}
volVectorField sumJDif
(
    IOobject
    (
        "sumJDif",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    // dimensionedVector ("UDiMiSum",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
    // dimensionedVector ("UDiMiSum",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
    // dimensionedVector ("jiDiff",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
    dimensionedVector ("sumJDif",dimLength/dimTime*dimMoles/dimVolume, vector(0.0,0.0,0.0))
    // dimensionedScalar ("dummyConcEqn",dimless,1.0)
);


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
    
    zoneLst[countI].correctBoundaryConditions();
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
    dimensionedScalar("sumZones",dimless,0.0),
    "zeroGradient"
); 


// -- porous media thermal diffusivity
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

// -- porosity
volVectorField permEff
(
    IOobject
    (
        "permEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("permEff",dimless/dimLength/dimLength,vector (0, 0, 0)),
    "zeroGradient"
);
// -- prepare the effective diffusion coefficient
volScalarField DEff
(
    IOobject
    (
        "DEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("DEff",dimLength*dimLength/dimTime,0.0),
    "zeroGradient"
    // (1.0-max(isOnWallCoat,isInWallCoat))*DifMix +
    // isOnWallCoat*(porEpsC/(tortC*(1.0/DifMix+1.0/DKn)))+
    // isInWallCoat*(porEpsS/(tortS*(1.0/DifMix+1.0/DKn))) 
);

// -- porous media heat capacity 
volScalarField CpS
(
IOobject
(
    "CpS",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
),
mesh,
dimensionedScalar("CpS",dimEnergy/dimMass/dimTemperature,0.0)
);

// -- porous media density 
volScalarField rhoS
(
IOobject
(
    "rhoS",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
),
mesh,
dimensionedScalar("rhoS",dimMass/dimVolume,0.0)
);

// -- heat transfer coefficient
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

// -- prepare sum of the molar fractions without reference
volScalarField sumWithoutRef
(
    IOobject
    (
        "sumWithoutRef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sumWithoutRef",dimless,0.0)
    // (1.0-max(isOnWallCoat,isInWallCoat))*DifMix +
    // isOnWallCoat*(porEpsC/(tortC*(1.0/DifMix+1.0/DKn)))+
    // isInWallCoat*(porEpsS/(tortS*(1.0/DifMix+1.0/DKn))) 
);



// -- efective kappa values, gas kappa in free channel, media kappa in zones
// kappaEff += freeChannel*kappa;         
forAll(zoneLst,countI)
{
    epsEff += zoneLst[countI]*(eps[countI]);
    // kappaEff += zoneLst[countI]*(kappaZ[countI]);
    kappaEff += zoneLst[countI]*(kappaZ[countI]);
    CpS += zoneLst[countI]*(CpZ[countI]);
    rhoS += zoneLst[countI]*(rhoZ[countI]);
    kHEff += zoneLst[countI]*(kH[countI]);
    aZoneEff += zoneLst[countI]*(aZone[countI]);
    permEff += zoneLst[countI]*perm[countI];
    sumZones += zoneLst[countI];
} 

volScalarField freeChannel = dimensionedScalar("jedna",dimless,1.0)-sumZones;
epsEff += freeChannel;
// volScalarField kappaS = kappaEff;

volScalarField kHEffA = kHEff*aZoneEff;


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

// -- reading field UD -- velocity enforced by diffusion
Info<< "Reading field UD\n" << endl;
volVectorField UD
(
    IOobject
    (
        "UD",
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
    // thermo.rho()
    // p/univR/T*molMList[molMList.size()-1]
    p/univR/thermo.T()*Mg
);

surfaceScalarField phiD
(
    IOobject
    (
        "phiD",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho*UD) & mesh.Sf()
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




