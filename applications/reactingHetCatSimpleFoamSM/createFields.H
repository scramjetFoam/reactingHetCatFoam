Info<< "Reading the concentrations \n" << endl;
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
    dimensionedScalar("Mg",dimMass/dimMoles,0.0)
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
    // Info<<molMList[countI]<<max(MFrList[countI])<<min(MFrList[countI])<<endl;

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

scalar n = specieNames.size();
scalar dimMat = n-1;

// -- matrix with the binary Ficks diffusion coefficients

Info << "Creating matrix with binary diffusivities without temperature and pressure dependence" << endl;
PtrList<dimensionedScalar> DijList((n*n));
for(int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        Info << "Seting matrix index (" << i << ", " << j << ")." <<endl;
        word nazevTu = "Dij" + word(i) + word(j);
        int index = (i*(n) + j);
        DijList.set
        (
            index,
            new dimensionedScalar(
                "DAB",
                Foam::pow(dimLength,2)/dimTime/Foam::pow(dimTemperature,1.75)*dimPressure,
                ((0.014)/(Foam::sqrt(2.0/(1e-3/molMList[i] + 1e-3/molMList[j]))*Foam::pow(Foam::pow(volDiffLst[i],1.0/3.0)+Foam::pow(volDiffLst[j],1.0/3.0),2))).value()
            )
        );
    }
}

Info << "Preparing matrix to store inv of B matrix" << endl;

// -- matrix with the inverse of the B matrix for S-M computations
PtrList<volScalarField> BinvList(dimMat*dimMat);
for(int i = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        if (j >= i)
        {
            Info << "Seting matrix index (" << i << ", " << j << ")." <<endl;
            word nazevTu = "Binv" + word(i) + word(j);
            int index = (i*(dimMat) + j);
            BinvList.set
            (
                index,
                new volScalarField
                (
                    IOobject
                    (
                        nazevTu,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    // NOTE: check dimensions
                    dimensionedScalar (nazevTu,dimLength*dimLength/dimTime, 0)
                    // dimensionedVector ("jiDiff",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
                )
            );
        }
    }
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

// -- porosity
volScalarField epsEff
(
IOobject
(
    "epsEff",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
),
mesh,
dimensionedScalar("epsEff",dimless,0.0)
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
    dimensionedVector("permEff",dimless/dimLength/dimLength,vector (0, 0, 0))
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
    dimensionedScalar("DEff",dimLength*dimLength/dimTime,0.0)
    // (1.0-max(isOnWallCoat,isInWallCoat))*DifMix +
    // isOnWallCoat*(porEpsC/(tortC*(1.0/DifMix+1.0/DKn)))+
    // isInWallCoat*(porEpsS/(tortS*(1.0/DifMix+1.0/DKn))) 
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
    epsEff += zoneLst[countI]*(eps[countI]);
    permEff += zoneLst[countI]*perm[countI];
}
volScalarField freeChannel = dimensionedScalar("jedna",dimless,1.0)-sumZones;
epsEff += freeChannel;


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


// // -- creating velocity enforced by diffusion for specie
// volVectorField sumFick
// (
//     IOobject
//     (
//         "sumFick",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     // dimensionedVector ("UDiMiSum",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
//     // dimensionedVector ("UDiMiSum",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
//     // dimensionedVector ("jiDiff",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
//     dimensionedVector ("sumFick",dimLength/dimTime*dimMoles/dimVolume, vector(0.0,0.0,0.0))
//     // dimensionedScalar ("dummyConcEqn",dimless,1.0)
// );

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

//-- surface field for kappaG interpolation

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
    p/univR/thermo.T()*Mg
);


#include "compressibleCreatePhi.H"

Info<< "Reading/calculating face flux field from diffusion flux phiD\n" << endl;
 
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
        IOobject::NO_WRITE
    ),
    // linearInterpolate(U*rho*Cp) & mesh.Sf()
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




