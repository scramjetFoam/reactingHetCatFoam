Info<< "Reading transport properties\n" << endl; 

// -- transportProperties dictionary
IOdictionary transportProperties
(
    IOobject
    (
	"transportProperties",    // dictionary name
	runTime.constant(),     // dict is found in "constant"
	mesh,                   // registry for the dict
	IOobject::MUST_READ,    // must exist, otherwise failure
	IOobject::NO_WRITE      // dict is only read by the solver
    )
);

// -- thermophysicalProperties dictionary
IOdictionary thermophysicalProperties
(
    IOobject
    (
	"thermophysicalProperties",    // dictionary name
	runTime.constant(),     // dict is found in "constant"
	mesh,                   // registry for the dict
	IOobject::MUST_READ,    // must exist, otherwise failure
	IOobject::NO_WRITE      // dict is only read by the solver
    )
);

//prepare the lists for materials
wordList specieNames
(
    transportProperties.lookup("species")
);

wordList zoneNames
(
    transportProperties.lookup("zones")
);


Info<< "Loaded species are\n" << specieNames << endl;
Info<< "Loaded zones are\n" << zoneNames << endl;

// -- prepare the output variables - volumetric diffusion coeficients for individual species
PtrList<dimensionedScalar> volDiffLst(specieNames.size()+1);		// list of volumetric diffusion coeficients + reference
PtrList<dimensionedScalar> molMList(specieNames.size()+1);		    // list of molatr masses + reference
    
//read the specie properties for each of them
forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    Info<< "Processing the specie: " << currentSp << endl;
    
    Info<< "Preparing volDiffLst and molMList for diffusivity computations" << endl;
    scalar sigmaV(readScalar(transportProperties.subDict(currentSp).lookup("sigmaV")));
    scalar molM(readScalar(transportProperties.subDict(currentSp).lookup("molM")));
    
    volDiffLst.set
    (
        countI, new dimensionedScalar("vollDiff",dimensionSet(0, 3, 0, 0, -1, 0, 0),sigmaV)
    );
    
    molMList.set
    (
        countI, new dimensionedScalar("molM", dimMass /dimMoles,molM)
    );

    Info<< "sigma V("<< currentSp <<") = "<< sigmaV << endl;
    Info<< "molar mass("<< currentSp <<") = "<< molM << endl;
}

// -- zone properties
PtrList<dimensionedScalar> eps(zoneNames.size());		// list of porosities 
PtrList<dimensionedScalar> tort(zoneNames.size());		// list of tortuosities 
PtrList<dimensionedScalar> kappaZ(zoneNames.size());		// list of effective heat conduitivities 
PtrList<dimensionedScalar> dp(zoneNames.size());		// list of pore diameter
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    eps.set
    (
        countI, new dimensionedScalar(transportProperties.subDict(currentZone).lookup("porEps"))
    );
    tort.set
    (
        countI, new dimensionedScalar(transportProperties.subDict(currentZone).lookup("tort"))
    );
    kappaZ.set
    (
        countI, new dimensionedScalar(transportProperties.subDict(currentZone).lookup("kappa"))
    );
    dp.set
    (
        countI, new dimensionedScalar(transportProperties.subDict(currentZone).lookup("dP"))
    );

}

// -- set reference sigma and molar mass
scalar sigmaVR(readScalar(transportProperties.subDict("ref").lookup("sigmaV")));
volDiffLst.set
(
    specieNames.size(), new dimensionedScalar("vollDiff",dimensionSet(0, 3, 0, 0, -1, 0, 0),sigmaVR)
);

scalar molMR(readScalar(thermophysicalProperties.subDict("mixture").subDict("specie").lookup("molWeight")));
molMList.set
(
    specieNames.size(), new dimensionedScalar("molM", dimMass /dimMoles,molMR*1e-3)
);

Info<< "sigma V(reference) = "<< sigmaVR << endl;
Info<< "molar mass(reference) = "<< molMR << endl;

dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//R constant

dimensionedScalar kappaA(transportProperties.subDict("ref").lookup("kappaA"));//R constant
dimensionedScalar kappaB(transportProperties.subDict("ref").lookup("kappaB"));//R constant
dimensionedScalar kappaC(transportProperties.subDict("ref").lookup("kappaC"));//R constant
dimensionedScalar kappaD(transportProperties.subDict("ref").lookup("kappaD"));//R constant
    
    
    
    
