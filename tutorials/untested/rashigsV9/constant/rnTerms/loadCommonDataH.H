// get the time connected stuff
const Time& time = mesh_.time();

//load the dictionary with the reaction data
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",      // dictionary name
        time.constant(),            // dict is found in "constant"
        mesh_,                      // registry for the dict
        IOobject::MUST_READ,        // must exist, otherwise failure
        IOobject::NO_WRITE          // dict is only read by the solver
    )
);

// read the list of species from the transport properties
wordList specieNames
(
    transportProperties.lookup("species")
);

labelList activeReacts
(
    transportProperties.subDict("reactProps").lookup("activeReacts")
);

// load heat capacity coefficients
PtrList<dimensionedScalar> ALst(specieNames.size()+1);			// list of coefficients for heat capacity
PtrList<dimensionedScalar> BLst(specieNames.size()+1);			// list of coefficients for heat capacity
PtrList<dimensionedScalar> CLst(specieNames.size()+1);			// list of coefficients for heat capacity
PtrList<dimensionedScalar> DLst(specieNames.size()+1);			// list of coefficients for heat capacity

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    scalar A(readScalar(transportProperties.subDict(currentSp).lookup("A")));
    scalar B(readScalar(transportProperties.subDict(currentSp).lookup("B")));
    scalar C(readScalar(transportProperties.subDict(currentSp).lookup("C")));
    scalar D(readScalar(transportProperties.subDict(currentSp).lookup("D")));
    ALst.set
    (
        countI,
        new dimensionedScalar("A" + currentSp, dimEnergy /dimMoles/Foam::pow(dimTemperature,1),A)
    );
    BLst.set
    (
        countI,
        new dimensionedScalar("B"+ currentSp, dimEnergy /dimMoles/Foam::pow(dimTemperature,2),B)
    );
    CLst.set
    (
        countI,
        new dimensionedScalar("C"+ currentSp, dimEnergy /dimMoles/Foam::pow(dimTemperature,3),C)
    );
    DLst.set
    (
        countI,
        new dimensionedScalar("D"+ currentSp, dimEnergy /dimMoles/Foam::pow(dimTemperature,4),D)
    );
}
scalar A(readScalar(transportProperties.subDict(specieNames[0]).lookup("AVR")));
scalar B(readScalar(transportProperties.subDict(specieNames[0]).lookup("BVR")));
scalar C(readScalar(transportProperties.subDict(specieNames[0]).lookup("CVR")));
scalar D(readScalar(transportProperties.subDict(specieNames[0]).lookup("DVR")));
ALst.set
(
    specieNames.size(),
    new dimensionedScalar("AVR", dimEnergy /dimMoles/Foam::pow(dimTemperature,1),A)
);
BLst.set
(
    specieNames.size(),
    new dimensionedScalar("BVR", dimEnergy /dimMoles/Foam::pow(dimTemperature,2),B)
);
CLst.set
(
    specieNames.size(),
    new dimensionedScalar("CVR", dimEnergy /dimMoles/Foam::pow(dimTemperature,3),C)
);
DLst.set
(
    specieNames.size(),
    new dimensionedScalar("DVR", dimEnergy /dimMoles/Foam::pow(dimTemperature,4),D)
);

// load the general environment properties
dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//gas constant
dimensionedScalar pAtm(transportProperties.subDict("genProps").lookup("pAtm"));//atmpspheric pressure
//~ dimensionedScalar cTemp(transportProperties.subDict("genProps").lookup("cTemp"));//current temperature

const volScalarField& T = mesh_.lookupObject<volScalarField>("T"); 
const volScalarField& p = mesh_.lookupObject<volScalarField>("p"); 

// calculate the total concentration
//~ dimensionedScalar cTot = pAtm/univR/cTemp;                              //total concentration
//~ dimensionedScalar rCTot = 1.0/cTot;                                  //reversed total concentration
volScalarField cTot(p/univR/T);

// load data for clipping the gradients (problem with potential infinite values)
scalar  clipGradMin(readScalar(transportProperties.subDict("reactProps").subDict("numProps").lookup("clipGradMin")));
scalar  clipGradMax(readScalar(transportProperties.subDict("reactProps").subDict("numProps").lookup("clipGradMax")));

// initialize the fields for reaction term and its derivative w.r.t
// current specie
volScalarField R
(
    IOobject
    (
        "R",
        time.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh_,
    dimensionedScalar("R",dimEnergy/dimTime/dimVolume,0.0)
);
volScalarField dRdT
(
    IOobject
    (
        "drdT",
        time.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh_,
    dimensionedScalar("drdT",dimEnergy/dimTime/dimVolume/dimTemperature,0.0)
);
        
volScalarField coatWallFrac
(
    IOobject
    (
        "coatWallFrac",
        time.timeName(0),
        // time.timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_
);
// initialize the fields for reaction heat source
// volScalarField rQ
// (
//     IOobject
//     (
//         "rQ",
//         time.timeName(),
//         mesh_,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh_,
//     dimensionedScalar("rQ",dimEnergy/dimTime/dimVolume,0.0)
// );


