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

// -- porosity
volVectorField permEff
(
    IOobject
    (
        "permEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption1
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
        writeOption1
    ),
    mesh,
    dimensionedScalar("DEff",dimLength*dimLength/dimTime,0.0),
    "zeroGradient"
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
        writeOption1
    ),
    mesh,
    dimensionedScalar("sumZones",dimless,0.0),
    "zeroGradient"
); 

forAll(zoneLst,countI)
{
    sumZones += zoneLst[countI];
    epsEff += zoneLst[countI]*(eps[countI]);
    permEff += zoneLst[countI]*perm[countI];
}
volScalarField freeChannel = dimensionedScalar("jedna",dimless,1.0)-sumZones;
epsEff += freeChannel;

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