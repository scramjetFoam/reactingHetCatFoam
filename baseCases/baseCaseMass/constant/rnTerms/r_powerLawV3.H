// -- read the corresponding row of stoechiometric matrix
List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

word currSpNameMass = currSpName + "Mass";
// -- load the current specie
const volScalarField&  currSpMass = mesh().lookupObject<volScalarField>(currSpNameMass);
const volScalarField&  currSp = mesh().lookupObject<volScalarField>(currSpName);

scalar molM(readScalar(transportProperties.subDict(currSpName).lookup("molM")));
dimensionedScalar molMDS("molMDS",dimMass/dimMoles,molM);

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
    // dimensionedScalar("R",dimMoles/dimTime/dimVolume,0.0)
    dimensionedScalar("R",dimMass/dimTime/dimVolume,0.0)
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
    // dimensionedScalar("gradR",dimMoles/dimTime/dimVolume,0.0)
    dimensionedScalar("gradR",dimMass/dimTime/dimVolume,0.0)
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

        // -- load the kinetics info about reaction
        dimensionedScalar k0  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("k0");
        dimensionedScalar EA  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("EA");
        
        // -- NOTETH:: NOT rewritten !
        if (kinSetupSwitch == 1)
        {
            // -- inhibition parameters
            dimensionedScalar kin  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("kin");
            dimensionedScalar Ein  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("Ein");
            
            Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = "<< EA.value()<<", kin = "<<kin.value()<<" and Ein = "<< Ein.value()<<")." << endl;
                
            // -- Arhenius law
            volScalarField k(k0*Foam::exp(-EA/(univR*T)));
            volScalarField Kin(kin*Foam::exp(Ein/(T)));

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
        else if(kinSetupSwitch == 2)
        {   
            // -- Arhenius law
            volScalarField k(k0*Foam::exp(-EA/(univR*T)));

            Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = " << EA.value() << ", min(k) = " << min(k).value() << ", max(k) = " << max(k).value() << ")." << endl;
                    
            //counting exponent for total concetration
            float cTotExp = 0.0;
            forAll(specieNames,countI){
                //this can be probably done nicer, but I used the code above  
                List<scalar> alphaLstCurSp(transportProperties.subDict(specieNames[countI]).lookup("alphaVec"));    
                cTotExp += alphaLstCurSp[rnI];
            }
            
            // calculating real kinetic constant of reaction with arrhenius kinetics and scaling with concentration for right dimension
            //~ dimensionedScalar kReact = k0*Foam::exp(-EA/(univR*cTemp))*Foam::pow(cTot,cTotExp);
            volScalarField kReact = k*Foam::pow(cTot,cTotExp)*molMDS;
            // volScalarField kReact = k*Foam::pow(rhoTu,cTotExp);
            
            volScalarField auxR(0.0*R);
            volScalarField auxGradR(0.0*gradR);
        
            auxR        += kReact * nuLst[rnI];
            auxGradR    += kReact * min(0.0,nuLst[rnI]);
                    
            Info << "reactants in reaction " << rnI << ": ";
            forAll(specieNames,countI)                                      //forAll other species
            {
                word testedSpName = specieNames[countI];
                List<scalar> nuLstTested(transportProperties.subDict(testedSpName).lookup("nuVec"));

                word testedSpNameMass = testedSpName + "Mass";
                
                if (nuLstTested[rnI] < -SMALL)                              // if testedSpecie is reactant
                {
                    Info << testedSpName << " ";
                    
                    List<scalar> alphaLst(transportProperties.subDict(testedSpName).lookup("alphaVec"));
                    const volScalarField& testedSp = mesh().lookupObject<volScalarField>(testedSpName); 
                
                    auxR *= Foam::pow(testedSp, alphaLst[rnI]);
                    
                    // if (testedSpName != currSpName)
                    // {
                    //     auxGradR *= Foam::pow(testedSp,alphaLst[rnI]);
                    // }
                    // else
                    // {
                    //     auxGradR *= alphaLst[rnI]*Foam::pow(testedSp,alphaLst[rnI]-1.0);
                    // }
                }
            }
            Info << endl;

            R     += auxR;
            gradR += auxGradR;
        }
        else if (kinSetupSwitch == 3)
        {
            dimensionedScalar atm("atm", dimPressure, 101325); //gas constant    //atmospheric pressure
            volScalarField k1R(k0*Foam::exp(-EA/(univR*T)));
            dimensionedScalar k2  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("k2");
            dimensionedScalar E2  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("E2");
            dimensionedScalar KEq  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("KEq");
            dimensionedScalar cCuCl  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("cCuCl");
            volScalarField k2R(k2*Foam::exp(-E2/(univR*T)));
            
            const volScalarField& ethylene = mesh().lookupObject<volScalarField>("ethylene"); 
            const volScalarField& O2 = mesh().lookupObject<volScalarField>("O2"); 
            
            volScalarField citatel( k1R * cCuCl * ethylene * O2 * p / atm );
            volScalarField jmenovatel( Foam::pow(ethylene + k1R/k2R * O2, 2) * (1 + KEq * p / atm * O2) ); 
            R += nuLst[rnI] * citatel / jmenovatel * molM;
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
    eqn += currZoneVol*coatWallFrac[countI]*fvm::Sp(gradR,currSpMass);
    eqn += currZoneVol*coatWallFrac[countI]*(R - gradR*currSpMass);
    // eqn += currZoneVol*coatWallFrac[countI]*(R);
    // eqn += molM*currZoneVol*coatWallFrac[countI]*fvm::Sp(gradR,currSp);
    // eqn += molM*currZoneVol*coatWallFrac[countI]*(R - gradR*currSp);

    // eqn += molM*currZoneVol*coatWallFrac[countI]*(R);
} 
