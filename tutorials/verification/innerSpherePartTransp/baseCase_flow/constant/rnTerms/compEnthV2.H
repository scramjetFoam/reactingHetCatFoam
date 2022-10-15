// enthalpy source and its temperature derivation
volScalarField Rent
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
    dimensionedScalar("R",dimEnergy/dimTime/dimVolume,0.0)
);
volScalarField dRdT
(
    IOobject
    (
        "drdT",
        time.timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh(),
    dimensionedScalar("drdT",dimEnergy/dimTime/dimVolume/dimTemperature,0.0)
);
forAll(specieNames,lstI)
{
      // -- get current specie name
      word currSpName = specieNames[lstI];                            
      // -- read the corresponding row of stoechiometric matrix
      List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

      // -- for all reactions in the nuLst
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
                dimensionedScalar sHr = reactiveProperties.subDict("reactProps").subDict(rName).lookup("sHr");
                if (kinSetupSwitch == 1)
                {
                    // -- inhibition parameters
                    dimensionedScalar kin  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("kin");
                    dimensionedScalar Ein  = reactiveProperties.subDict("reactProps").subDict(rName).lookup("Ein");
                    
                    Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = "<< EA.value()<<", kin = "<<kin.value()<<" and Ein = "<< Ein.value()<<"), sHr"<<sHr.value()<< endl;
                        
                    // -- Arhenius law
                    volScalarField k(k0*Foam::exp(-EA/(univR*T)));
                    volScalarField Kin(kin*Foam::exp(Ein/(T)));

                    // -- load molar fraction fields
                    const volScalarField& CO = mesh().lookupObject<volScalarField>("CO"); 
                    const volScalarField& O2 = mesh().lookupObject<volScalarField>("O2");

                    // -- heat source and its linearized part
                    Rent = -(sHr*k*CO*O2*Foam::pow(1+Kin*CO,-2));
                    dRdT = -(sHr*(k0*CO*O2*(CO*EA*Foam::exp(-(Ein*univR+EA)/(univR*T))*kin-2*Foam::exp(-(Ein*univR+EA)/(univR*T))*CO*kin*Ein*univR+EA*Foam::exp(-EA/(univR*T)))/(univR*Foam::pow(T,2)*Foam::pow(1+kin*Foam::exp(-Ein/T)*CO,3))));
                }

                else if(kinSetupSwitch == 2)
                {                   
                    // -- Arhenius law
                    volScalarField k(k0*Foam::exp(-EA/(univR*T)));

                    Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = " << EA.value() << ", min(k) = " << min(k).value() << ", max(k) = " << max(k).value() << ")." << endl;
                            
                    //counting exponent for total concetration
                    //NOTE: I start with -1.0 because of units of kinetic constants 
                    //      in transport properties (maybe there can be problem with
                    //      non-integer alpha)
                    float cTotExp = 0.0;
                    forAll(specieNames,countI){
                        //this can be probably done nicer, but I used the code above  
                        List<scalar> alphaLstCurSp(transportProperties.subDict(specieNames[countI]).lookup("alphaVec"));    
                        cTotExp += alphaLstCurSp[rnI];
                    }
                    
                    // calculating real kinetic constant of reaction with arrhenius kinetics and scaling with concentration for right dimension
                    //~ dimensionedScalar kReact = k0*Foam::exp(-EA/(univR*cTemp))*Foam::pow(cTot,cTotExp);
                    volScalarField kReact = k*Foam::pow(cTot,cTotExp);
                    
                    volScalarField auxR(0.0*Rent);

                    // NOTETH: I am not using linearization of the source term at the time
                    // volScalarField auxGradR(0.0*dRdT);
                
                    auxR        += nuLst[rnI]*kReact*sHr;
                    // auxGradR    += min(0.0,nuLst[rnI])*kReact*sHr;   
                            
                    Info << "reactants in reaction " << rnI << ": ";
                    forAll(specieNames,countI)                                      //forAll other species
                    {
                        word testedSpName = specieNames[countI];
                        List<scalar> nuLstTested(transportProperties.subDict(testedSpName).lookup("nuVec"));
                        
                        if (nuLstTested[rnI] < -SMALL)                              // if testedSpecie is reactant
                        {
                            Info << testedSpName << " ";
                            
                            List<scalar> alphaLst(transportProperties.subDict(testedSpName).lookup("alphaVec"));
                            const volScalarField& testedSp = mesh().lookupObject<volScalarField>(testedSpName); 
                        
                            auxR *= Foam::pow(testedSp,alphaLst[rnI]);
                            
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

                    Rent     += auxR;
                    // gradR += auxGradR;
                }
                // turn of the current reaction
                activeReacts[rnI] = 0.0;

                Info << "reaction " << rnI << " was accounted for and turned off" << nChars<< endl;
            }
            else
            {   
                  Info << currSpName << " is not active in reaction " << rnI<< " or reaction is not active"<<endl;
            }
      }
}

// bounding the gradR field by maximal and minimal values            
// dimensionedScalar minGradR("minGrad",dRdT.dimensions(), clipGradMin);
// dimensionedScalar maxGradR("maxGrad",dRdT.dimensions(), clipGradMax);
// dRdT = min(max(dRdT,minGradR),maxGradR);
// Info <<"rQ - enthalphy source: "<< min(rQ) <<max(rQ)<< endl;         
// eqn += rQ;
// eqn += coatWallFrac*fvm::Sp(dRdT,T);
// eqn += (coatWallFrac*R - coatWallFrac*dRdT*T);
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    const volScalarField& currZoneVol = mesh().lookupObject<volScalarField>(currentZone); 
    // eqn += currZoneVol*coatWallFrac[countI]*fvm::Sp(dRdT,T);
    // eqn += currZoneVol*coatWallFrac[countI]*(Rent - dRdT*T);

    eqn += currZoneVol*coatWallFrac[countI]*(Rent);
} 
// Note: I add this source as purely explicit
