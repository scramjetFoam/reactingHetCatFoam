// -- get the time connected stuff
const Time& time = mesh().time();

// -- load the dictionary with the reaction data
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",      // dictionary name
        time.constant(),            // dict is found in "constant"
        mesh(),                      // registry for the dict
        IOobject::MUST_READ,        // must exist, otherwise failure
        IOobject::NO_WRITE          // dict is only read by the solver
    )
);

// -- load the thermophysicalProperties dictionary
IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",      // dictionary name
        time.constant(),            // dict is found in "constant"
        mesh(),                      // registry for the dict
        IOobject::MUST_READ,        // must exist, otherwise failure
        IOobject::NO_WRITE          // dict is only read by the solver
    )
);

// -- load the reactiveProperties dictionary
IOdictionary reactiveProperties
(
    IOobject
    (
        "reactiveProperties",      // dictionary name
        time.constant(),            // dict is found in "constant"
        mesh(),                      // registry for the dict
        IOobject::MUST_READ,        // must exist, otherwise failure
        IOobject::NO_WRITE          // dict is only read by the solver
    )
);

// -- read the list of species from the transport properties
wordList specieNames
(
    transportProperties.lookup("species")
);

// -- read the list with different zones
wordList zoneNames
(
    transportProperties.lookup("zones")
);

// -- list with coating fraction in zone
PtrList<dimensionedScalar> coatWallFrac(zoneNames.size());			
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    coatWallFrac.set
    (
        countI, new dimensionedScalar(transportProperties.subDict(currentZone).lookup("coatFrac"))
    );
}

// -- active reactions
labelList activeReacts
(
    reactiveProperties.subDict("reactProps").lookup("activeReacts")
);

// -- load the general environment properties
dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//gas constant
dimensionedScalar atm("atm", dimpressure, 101325);//gas constant    //atmospheric pressure

// -- reference molar mass 
scalar molMR(readScalar(thermophysicalProperties.subDict("mixture").subDict("specie").lookup("molWeight")));
dimensionedScalar molMRef("molMRef",dimMass/dimMoles,molMR*1e-3);

// -- load temperature and pressure fields
const volScalarField& T = mesh().lookupObject<volScalarField>("T"); 
const volScalarField& p = mesh().lookupObject<volScalarField>("p");  
const volScalarField& rhoTu = mesh().lookupObject<volScalarField>("rho");  

// -- calculate the total concentration
volScalarField cTot(p/T/univR);

// load data for clipping the gradients (problem with potential infinite values)
scalar  clipGradMin(readScalar(reactiveProperties.subDict("reactProps").subDict("numProps").lookup("clipGradMin")));
scalar  clipGradMax(readScalar(reactiveProperties.subDict("reactProps").subDict("numProps").lookup("clipGradMax")));



            
