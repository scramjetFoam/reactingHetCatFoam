// -- debug level, so far I have 3 debug levels, see ifs below
int debug = 0;

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
        writeOption1
    ),
    // thermo.rho()
    p/univR/thermo.T()*Mg
);

// volScalarField cTot
// (
//     IOobject
//     (
//         "cTot",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         writeOption1
//     ),
//     // thermo.rho()
//     p/univR/thermo.T()
// );


#include "compressibleCreatePhi.H"

Info<< "Reading/calculating face flux field from diffusion flux phiD\n" << endl;
 
// surfaceScalarField phiD
// (
//     IOobject
//     (
//         "phiD",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         writeOption2
//     ),
//     linearInterpolate(rho*UD) & mesh.Sf()
// );

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

surfaceScalarField phiD
(
    IOobject
    (
        "phiD",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        writeOption2
    ),
    linearInterpolate(rho*UD) & mesh.Sf()
);



