// -- porosity
volScalarField epsEff
(
    IOobject
    (
        "epsEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption1
    ),
    mesh,
    dimensionedScalar("epsEff",dimless,0.0),
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
        writeOption1
    ),
    mesh,
    dimensionedScalar("DEff",dimLength*dimLength/dimTime,0.0),
    "zeroGradient"
);
DEff.correctBoundaryConditions();

// -- free channels --> free channels = 1 - sum all zones 
volScalarField sumZones
(
    IOobject
    (
        "sumZones",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    dimensionedScalar("sumZones",dimless,0.0),
    "zeroGradient"
); 

// -- prepare the effective electrochemical diffusive coefficient
volScalarField sigmaEffEl
(
    IOobject
    (
        "sigmaEffEl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption1
    ),
    mesh,
    dimensionedScalar("sigmaEffEl",sigmaScInfEl[0].dimensions(),0.0),
    "zeroGradient"
);
sigmaEffEl.correctBoundaryConditions();

// -- prepare the effective electrochemical diffusive coefficient
volScalarField sigmaEffIon
(
    IOobject
    (
        "sigmaEffIon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption1
    ),
    mesh,
    dimensionedScalar("sigmaEffIon",sigmaScInfIon[0].dimensions(),0.0),
    "zeroGradient"
);
sigmaEffIon.correctBoundaryConditions();


forAll(zoneLst,countI)
{
    sumZones += zoneLst[countI];
    epsEff += zoneLst[countI]*(eps[countI]);
}
volScalarField freeChannel = dimensionedScalar("jedna",dimless,1.0)-sumZones;
epsEff += freeChannel;

// -- surface field for kappaG interpolation
surfaceScalarField epsGrad
(
    IOobject
    (
        "epsGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    // linearInterpolate(U*rho*Cp) & mesh.Sf()
    fvc::snGrad(epsEff)
);

// // -- prepare Knudsen diffusivity
volScalarField DKn
(
    IOobject
    (
        "DKn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DKn",dimLength*dimLength/dimTime,0.0)
);



