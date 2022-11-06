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

forAll(zoneLst,countI)
{
    sumZones += zoneLst[countI];
    epsEff += zoneLst[countI]*(eps[countI]);
    permEff += zoneLst[countI]*perm[countI];
}

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