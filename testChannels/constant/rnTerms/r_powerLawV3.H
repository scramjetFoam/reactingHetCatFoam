// -- read the corresponding row of stoechiometric matrix
List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

// -- load the current specie
const volScalarField&  currSp = mesh().lookupObject<volScalarField>(currSpName);

// -- load data for numerics from transportProperties
scalar  minCForGradND(readScalar(transportProperties.subDict(currSpName).lookup("minForGrad")));
dimensionedScalar minCForGrad("minCForGrad",currSp.dimensions(),minCForGradND);
volScalarField cForGrad(max(currSp,minCForGrad));

// -- initialize the fields for reaction term and its derivative w.r.t
volScalarField R
(
    IOobject
    (
        "R",
        time.timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh(),
    dimensionedScalar("R",dimMoles/dimTime/dimVolume,0.0)
);
volScalarField gradR
(
    IOobject
    (
        "gradR",
        time.timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh(),
    dimensionedScalar("gradR",dimMoles/dimTime/dimVolume,0.0)
);

// -- for all reactions in nuLst (trasportProperties specie defined)
forAll(nuLst,rnI)                                                      
{    
    if (mag(nuLst[rnI]) > SMALL and activeReacts[rnI])                  //tango only in active reactions
    {   
        // -- name of the reaction
        char rName [50];
        label nChars = sprintf(rName,"reaction%02d",rnI);
        
        // -- load the kinetics info about reaction
        dimensionedScalar k0  = transportProperties.subDict("reactProps").subDict(rName).lookup("k0");
        dimensionedScalar EA  = transportProperties.subDict("reactProps").subDict(rName).lookup("EA");
	    dimensionedScalar kin  = transportProperties.subDict("reactProps").subDict(rName).lookup("kin");
        dimensionedScalar Ein  = transportProperties.subDict("reactProps").subDict(rName).lookup("Ein");
        
        Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = "<< EA.value()<<", kin = "<<kin.value()<<" and Ein = "<< Ein.value()<<")" << endl;
            
        // -- kinetics from article 3D reconstruction ...
        volScalarField k(k0*Foam::exp(-EA/(univR*T)));
        volScalarField Kin(kin*Foam::exp(Ein/(T)));

        // -- load molar fraction fields
        const volScalarField& CO = mesh().lookupObject<volScalarField>("CO"); 
        const volScalarField& O2 = mesh().lookupObject<volScalarField>("O2");
    
        // -- calculate reaction term
        R	= k*nuLst[rnI]*CO*O2*Foam::pow(1+Kin*CO,-2);
        if (currSpName == "CO"){
            // Info << "CO: nu="<<nuLst[rnI] << endl;
            gradR = nuLst[rnI]*(k*O2*(1-Kin*CO)*Foam::pow(1+Kin*CO,-3));
        }
        else if (currSpName == "O2")
        {
            // Info << "O2: nu="<<nuLst[rnI] <<endl;
            gradR   = nuLst[rnI]*(k*CO*Foam::pow(1+Kin*CO,-2));
        }
        else
        {
            // Info << "CO2: nu="<<nuLst[rnI] <<endl;
            gradR   = 0.0*nuLst[rnI]*(k*(CO)*(O2));
        }
    }
    else
    {
        if (activeReacts[rnI])
        {
            Info << currSpName << " is not active in reaction " << rnI << endl;
        }
        else
        {
            Info << "reaction " << rnI << " switched-off" << endl;
        }
    }
}

// bounding the gradR field by maximal and minimal values            
dimensionedScalar minGradR("minGrad",gradR.dimensions(), clipGradMin);
dimensionedScalar maxGradR("maxGrad",gradR.dimensions(), clipGradMax);
gradR = min(max(gradR,minGradR),maxGradR);

// Info<<"min: "<< min(currSp);
// Info<<"max: "<< max(currSp) <<endl;

// eqn += isCoat*coatingVolFrac*fvm::Sp(dummyGrad*gradR,currSp);
// Info << max(gradR)<<min(gradR) << endl;
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    const volScalarField& currZoneVol = mesh().lookupObject<volScalarField>(currentZone); 
    eqn += currZoneVol*coatWallFrac[countI]*fvm::Sp(gradR,currSp);
    eqn += currZoneVol*coatWallFrac[countI]*(R - gradR*currSp);

    // explicit
    // eqn += currZoneVol*coatWallFrac[countI]*(R);
} 


 
// Note: Taylor expansion of kinetics
