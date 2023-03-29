// -- debug level, so far I have 3 debug levels, see ifs below
int debug = 1;

if (simple.dict().found("debug"))
{
    debug = readScalar(simple.dict().lookup("debug"));
}

Info << "Debug level: " << debug <<endl;

IOobject::writeOption writeOption1 = IOobject::NO_WRITE;
IOobject::writeOption writeOption2 = IOobject::NO_WRITE;
if (debug > 0)
{
    writeOption1 = IOobject::AUTO_WRITE;
}
if (debug == 2)
{
    writeOption2 = IOobject::AUTO_WRITE;
}

#include "createFieldsFromMesh.H"

#include "createFieldsFromTransProp.H"

#include "createEffProp.H"


// -- dummy volScalarFields to resolve dimension problem in conc.H fvOptions inconsistent source
//NOTE: I dont know why this is necessary, dimensions are correct
//but openFOAM doesnt think so :(
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
    // dimensionedScalar ("dummyConcEqn",dimMoles/dimVolume,1.0)
    // dimensionedScalar ("dummyConcEqn",dimMoles*dimMoles/dimMass/dimMass/dimVolume,1.0)
    dimensionedScalar ("dummyConcEqn",dimMass/dimVolume,1.0)
);

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
    // mesh
);


#include "compressibleCreatePhi.H"

Info<< "Reading/calculating face flux field from diffusion flux phiD\n" << endl;

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

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// -- matrix with the binary Fick diffusion coefficients
Info << "Creating matrix with binary diffusivities without temperature and pressure dependence" << endl;
PtrList<dimensionedScalar> DijList(((n)*(n)));
for(int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
    {
        Info << "Seting matrix index (" << i << ", " << j << ")." <<endl;
        word nazevTu = "D" + name(i) + name(j);
        int index = (i*(n) + j);
        DijList.set
        (
            index,
            new dimensionedScalar(
                nazevTu,
                Foam::pow(dimLength,2)/dimTime/Foam::pow(dimTemperature,1.75)*dimPressure,
                ((0.014)/(Foam::sqrt(2.0/(1e-3/molMList[i] + 1e-3/molMList[j]))*Foam::pow(Foam::pow(volDiffLst[i],1.0/3.0)+Foam::pow(volDiffLst[j],1.0/3.0),2))).value()
            )
        );
    }
}

Info << "Creating matrix of mass diffusive fluxes relative to mass average velocity." << endl;
PtrList<volVectorField> JiList(n-1);
for(int i = 0; i < n-1; i++)
{
    Info << "Seting diffusive flux of " << specieNames[i] << endl;
    word nazevTu = "J" + specieNames[i];
    JiList.set
    (
        i,
        new volVectorField
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
            dimensionedVector (nazevTu,dimMass/dimTime/dimLength/dimLength, vector(0,0,0))
        )
    );
}

Info << "Creating matrix of mass diffusive fluxes relative to molar average velocity." << endl;
PtrList<volVectorField> JUiList(n-1);
for(int i = 0; i < n-1; i++)
{
    Info << "Seting diffusive flux of " << specieNames[i] << endl;
    word nazevTu = "JU" + specieNames[i];
    JUiList.set
    (
        i,
        new volVectorField
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
            dimensionedVector (nazevTu,dimMass/dimTime/dimLength/dimLength, vector(0,0,0))
        )
    );
}

// -- matrix with the inverse of the B matrix for S-M computations
Info << "Preparing matrix to store inv of B matrix" << endl;
PtrList<volScalarField> BinvList((dimMat)*(dimMat));
for(int i = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        // if (j >= i)
        // {
            Info << "Seting matrix index (" << i << ", " << j << ")." <<endl;
            word nazevTu = "Binv" + name(i) + name(j);
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
                    dimensionedScalar (nazevTu,dimLength*dimLength/dimTime, 0),
                    "zeroGradient"
                    // dimensionedVector ("jiDiff",dimLength/dimTime*dimMass/dimVolume, vector(0.0,0.0,0.0))
                )
            );
        // }
    }
}

