// -- read the corresponding row of stoechiometric matrix
List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

// -- load the current specie
const volScalarField&  currSp = mesh().lookupObject<volScalarField>(currSpName);

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

        int kinSetupSwitch(readScalar(reactiveProperties.subDict("reactProps").subDict(rName).lookup("kinSetupSwitch")));

        if (kinSetupSwitch == 1){
            // -- load the kinetics info about reaction
            dimensionedScalar k0  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("k0");
            dimensionedScalar EA  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("EA");
            dimensionedScalar kin  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("kin");
            dimensionedScalar Ein  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("Ein");
            
            Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = "<< EA.value()<<", kin = "<<kin.value()<<" and Ein = "<< Ein.value()<<")." << endl;
                
            // -- Arhenius law
            volScalarField k(k0*Foam::exp(-EA/(univR*TSolid)));
            volScalarField Kin(kin*Foam::exp(Ein/(TSolid)));

            // -- load molar fraction fields
            const volScalarField& CO = mesh().lookupObject<volScalarField>("CO"); 
            const volScalarField& O2 = mesh().lookupObject<volScalarField>("O2");
        
            // -- calculate the reaction term
            R += k*nuLst[rnI]*CO*O2*Foam::pow(1+Kin*CO,-2);

            // -- calculate the gradient for source linearization
            if (currSpName == "CO"){
                gradR += nuLst[rnI]*(k*O2*(1-Kin*CO)*Foam::pow(1+Kin*CO,-3));
            }
            else if (currSpName == "O2")
            {
                gradR += nuLst[rnI]*(k*CO*Foam::pow(1+Kin*CO,-2));
            }
            else
            {
                gradR += 0.0*nuLst[rnI]*(k*(CO)*(O2));
            }
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

// -- bounding the gradR field by maximal and minimal values            
dimensionedScalar minGradR("minGrad",gradR.dimensions(), clipGradMin);
dimensionedScalar maxGradR("maxGrad",gradR.dimensions(), clipGradMax);
gradR = min(max(gradR,minGradR),maxGradR);

// -- reaction occuring inside porous media
// -- linearized
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    const volScalarField& currZoneVol = mesh().lookupObject<volScalarField>(currentZone); 
    eqn += currZoneVol*coatWallFrac[countI]*fvm::Sp(gradR,currSp);
    eqn += currZoneVol*coatWallFrac[countI]*(R - gradR*currSp);
} 
