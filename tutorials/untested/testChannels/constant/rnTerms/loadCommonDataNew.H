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

// load the general environment properties
dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//gas constant
dimensionedScalar pAtm(transportProperties.subDict("genProps").lookup("pAtm"));//atmpspheric pressure
dimensionedScalar cTemp(transportProperties.subDict("genProps").lookup("cTemp"));//current temperature

const volScalarField& T = mesh_.lookupObject<volScalarField>("T"); 
// const volScalarField& p = mesh_.lookupObject<volScalarField>("p"); 
const volScalarField& rhoTu = mesh_.lookupObject<volScalarField>("rho"); 

// calculate the total concentration
// dimensionedScalar cTot = pAtm/univR/cTemp;                              //total concentration
//~ dimensionedScalar rCTot = 1.0/cTot;                                  //reversed total concentration
// volScalarField cTot(pAtm/univR/T);
dimensionedScalar  MTu("MTu",dimMass*Foam::pow(dimMoles,-1),28.96e-3);
volScalarField cTot(rhoTu/MTu);

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
    dimensionedScalar("R",dimMoles/dimTime/dimVolume,0.0)
);
volScalarField gradR
(
    IOobject
    (
        "gradR",
        time.timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh_,
    dimensionedScalar("gradR",dimMoles/dimTime/dimVolume,0.0)
    // dimensionedScalar("gradR",dimMoles/dimTime,0.0)
);
            
