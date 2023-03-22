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

// -- porosityProperties dictionary
IOdictionary porosityProperties
(
    IOobject
    (
	"porosityProperties",    // dictionary name
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

scalar n = specieNames.size();
scalar dimMat = n-1;

Info<< "Loaded species are\n" << specieNames << endl;
Info<< "Loaded zones are\n" << zoneNames << endl;

// -- type of the assumed diffusion
// ---- 1) Fick -- Fick diffusion without convective fluxes forced by diffussion
// ---- 2) FickConv -- Fick diffusion with convective fluxes forced by diffussion
// ---- 3) SMConv -- Stefan-Maxwell with convective fluxes forced by diffussion

// -- prepare the output variables - volumetric diffusion coeficients for individual species
PtrList<dimensionedScalar> volDiffLst(specieNames.size());		// list of volumetric diffusion coeficients + reference
PtrList<dimensionedScalar> molMList(specieNames.size());		    // list of molatr masses + reference
    
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
dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//R constant

dimensionedScalar kappaA(transportProperties.subDict("ref").lookup("kappaA"));//R constant
dimensionedScalar kappaB(transportProperties.subDict("ref").lookup("kappaB"));//R constant
dimensionedScalar kappaC(transportProperties.subDict("ref").lookup("kappaC"));//R constant
dimensionedScalar kappaD(transportProperties.subDict("ref").lookup("kappaD"));//R constant
    
    
// -- prepare free difusion coeficients of the current specie against all other species 
// -- here without T and p dependence
PtrList<dimensionedScalar> DABLst(specieNames.size());
forAll(specieNames,lstI)
{
    // diffusion coefficients according to Fuller without temperature and pressure dependence
    DABLst.set
    (
        lstI,
        new dimensionedScalar(
            "DAB",
            Foam::pow(dimLength,2)/dimTime/Foam::pow(dimTemperature,1.75)*dimPressure,
            ((0.014)/(Foam::sqrt(2.0/(1e-3/molMList[lstI] + 1e-3/molMList[specieNames.size()-1]))*Foam::pow(Foam::pow(volDiffLst[lstI],1.0/3.0)+Foam::pow(volDiffLst[specieNames.size()-1],1.0/3.0),2))).value()
        )
    );
}   
        
