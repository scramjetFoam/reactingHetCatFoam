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
// PtrList<dimensionedScalar> ALst(specieNames.size()+1);			    // list of coefficients for heat capacity
// PtrList<dimensionedScalar> BLst(specieNames.size()+1);			    // list of coefficients for heat capacity
// PtrList<dimensionedScalar> CLst(specieNames.size()+1);			    // list of coefficients for heat capacity
// PtrList<dimensionedScalar> DLst(specieNames.size()+1);			    // list of coefficients for heat capacity

    
//read the specie properties for each of them
forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    Info<< "Processing the specie: " << currentSp << endl;
    
    Info<< "Preparing volDiffLst and molMList for diffusivity computations" << endl;
    scalar sigmaV(readScalar(transportProperties.subDict(currentSp).lookup("sigmaV")));
    scalar molM(readScalar(transportProperties.subDict(currentSp).lookup("molM")));
    // scalar A(readScalar(transportProperties.subDict(currentSp).lookup("A")));
    // scalar B(readScalar(transportProperties.subDict(currentSp).lookup("B")));
    // scalar C(readScalar(transportProperties.subDict(currentSp).lookup("C")));
    // scalar D(readScalar(transportProperties.subDict(currentSp).lookup("D")));
    
    volDiffLst.set
    (
        countI, new dimensionedScalar("vollDiff",dimensionSet(0, 3, 0, 0, -1, 0, 0),sigmaV)
    );
    
    molMList.set
    (
        countI, new dimensionedScalar("molM", dimMass /dimMoles,molM)
    );
    // ALst.set
    // (
	//     countI, new dimensionedScalar("ALst", dimEnergy /dimMoles/Foam::pow(dimTemperature,1),A)
    // );
    // BLst.set
    // (
	//     countI, new dimensionedScalar("BLst", dimEnergy /dimMoles/Foam::pow(dimTemperature,2),B)
    // );
    // CLst.set
    // (
	//     countI, new dimensionedScalar("CLst", dimEnergy /dimMoles/Foam::pow(dimTemperature,3),C)
    // );
    // DLst.set
    // (
	//     countI, new dimensionedScalar("DLst", dimEnergy /dimMoles/Foam::pow(dimTemperature,4),D)
    // );
    Info<< "sigma V("<< currentSp <<") = "<< sigmaV << endl;
    Info<< "molar mass("<< currentSp <<") = "<< molM << endl;
    // Info<< "cp("<< currentSp <<") = "<< A<<" + "<< B<< "* T +" << C << " *T**2 + "<<D<<"*T**3"<< endl;
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

// -- set the reference heat capacity coeficients
// scalar A(readScalar(transportProperties.subDict("ref").lookup("A")));
// scalar B(readScalar(transportProperties.subDict("ref").lookup("B")));
// scalar C(readScalar(transportProperties.subDict("ref").lookup("C")));
// scalar D(readScalar(transportProperties.subDict("ref").lookup("D")));
// ALst.set
// (
//     specieNames.size(),
//     new dimensionedScalar("ALst", dimEnergy /dimMoles/Foam::pow(dimTemperature,1),A)
// );
// BLst.set
// (
//     specieNames.size(),
//     new dimensionedScalar("BLst", dimEnergy /dimMoles/Foam::pow(dimTemperature,2),B)
// );
// CLst.set
// (
//     specieNames.size(),
//     new dimensionedScalar("CLst", dimEnergy /dimMoles/Foam::pow(dimTemperature,3),C)
// );
// DLst.set
// (
//     specieNames.size(),
//     new dimensionedScalar("DLst", dimEnergy /dimMoles/Foam::pow(dimTemperature,4),D)
// );

Info<< "sigma V(reference) = "<< sigmaVR << endl;
Info<< "molar mass(reference) = "<< molMR << endl;
// Info<< "cp(reference) = "<< A<<" + "<< B<< "* T +" << C << " *T**2 + "<<D<<"*T**3"<< endl;


// load the coating properties -- overwriten by general zones
// dimensionedScalar porEpsC(transportProperties.subDict("coatingProps").lookup("porEps"));
// // dimensionedScalar constrC(transportProperties.subDict("coatingProps").lookup("constr"));
// dimensionedScalar tortC(transportProperties.subDict("coatingProps").lookup("tort"));
// dimensionedScalar kappaTC(transportProperties.subDict("coatingProps").lookup("kappaT"));
// dimensionedScalar CpC(transportProperties.subDict("coatingProps").lookup("Cp"));
// dimensionedScalar rhoC(transportProperties.subDict("coatingProps").lookup("rho"));
// dimensionedScalar dpC(transportProperties.subDict("coatingProps").lookup("dP"));

// load the substrate properties
// dimensionedScalar porEpsS(transportProperties.subDict("substrateProps").lookup("porEps"));
// // dimensionedScalar constrS(transportProperties.subDict("substrateProps").lookup("constr"));
// dimensionedScalar tortS(transportProperties.subDict("substrateProps").lookup("tort"));
// dimensionedScalar kappaTS(transportProperties.subDict("substrateProps").lookup("kappaT"));
// dimensionedScalar CpS(transportProperties.subDict("substrateProps").lookup("Cp"));
// dimensionedScalar rhoS(transportProperties.subDict("substrateProps").lookup("rho"));
// dimensionedScalar dpS(transportProperties.subDict("substrateProps").lookup("dP"));

// load the on-wall and in-wall properties
// scalar  coatingVolFracOW(readScalar(transportProperties.subDict("onWallProperties").lookup("coatingVolFrac")));
// scalar  coatingVolFracIW(readScalar(transportProperties.subDict("inWallProperties").lookup("coatingVolFrac")));

// scalar nItconc(readScalar(transportProperties.subDict("nIterations").lookup("conc")));
// scalar nIttemp(readScalar(transportProperties.subDict("nIterations").lookup("temp")));
    
//~ dimensionedScalar cTemp(transportProperties.subDict("genProps").lookup("cTemp"));//current temperature
// dimensionedScalar pAtm(transportProperties.subDict("genProps").lookup("pAtm"));//pressure
dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//R constant
// dimensionedScalar rhoCons(transportProperties.subDict("genProps").lookup("rhoCons"));//R constant
//~ dimensionedScalar Cp(transportProperties.subDict("genProps").lookup("Cp"));//R constant
// dimensionedScalar kappa(transportProperties.subDict("genProps").lookup("kappa"));//R constant
dimensionedScalar kappaA(transportProperties.subDict("ref").lookup("kappaA"));//R constant
dimensionedScalar kappaB(transportProperties.subDict("ref").lookup("kappaB"));//R constant
dimensionedScalar kappaC(transportProperties.subDict("ref").lookup("kappaC"));//R constant
dimensionedScalar kappaD(transportProperties.subDict("ref").lookup("kappaD"));//R constant
    
    
    
    
