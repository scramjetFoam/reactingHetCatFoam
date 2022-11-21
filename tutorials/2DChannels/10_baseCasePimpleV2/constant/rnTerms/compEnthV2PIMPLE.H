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
                // Info << currSpName << " is active in reaction " << rnI << endl;

                char rName [50];
                label nChars = sprintf(rName,"reaction%02d",rnI);
                
                // -- load data from transport propertis
                dimensionedScalar k0 = reactiveProperties.subDict("reactProps").subDict(rName).lookup("k0");
                dimensionedScalar EA = reactiveProperties.subDict("reactProps").subDict(rName).lookup("EA");
                dimensionedScalar kin = reactiveProperties.subDict("reactProps").subDict(rName).lookup("kin");
                dimensionedScalar Ein = reactiveProperties.subDict("reactProps").subDict(rName).lookup("Ein");
                dimensionedScalar sHr = reactiveProperties.subDict("reactProps").subDict(rName).lookup("sHr");
                Info << currSpName << " is active in reaction " << rnI << " "<<rName<<" "<< " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = "<< EA.value()<<", kin = "<<kin.value()<<" and Ein = "<< Ein.value()<<"), sHr"<<sHr.value()<< endl;

                // -- calculate reaction rate constants 
                volScalarField k(k0*Foam::exp(-EA/(univR*TSolid)));
                volScalarField Kin(kin*Foam::exp(Ein/(TSolid)));

                // --concentration fields
                const volScalarField& CO = mesh().lookupObject<volScalarField>("CO"); 
                const volScalarField& O2 = mesh().lookupObject<volScalarField>("O2"); 

                // -- heat source and its linearized part
                Rent = -(sHr*k*CO*O2*Foam::pow(1+Kin*CO,-2));
                dRdT = -(sHr*(k0*CO*O2*(CO*EA*Foam::exp(-(Ein*univR+EA)/(univR*TSolid))*kin-2*Foam::exp(-(Ein*univR+EA)/(univR*TSolid))*CO*kin*Ein*univR+EA*Foam::exp(-EA/(univR*TSolid)))/(univR*Foam::pow(TSolid,2)*Foam::pow(1+kin*Foam::exp(-Ein/TSolid)*CO,3))));

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
dimensionedScalar minGradR("minGrad",dRdT.dimensions(), clipGradMin);
dimensionedScalar maxGradR("maxGrad",dRdT.dimensions(), clipGradMax);
dRdT = min(max(dRdT,minGradR),maxGradR);
// Info <<"rQ - enthalphy source: "<< min(rQ) <<max(rQ)<< endl;         
// eqn += rQ;
// eqn += coatWallFrac*fvm::Sp(dRdT,T);
// eqn += (coatWallFrac*R - coatWallFrac*dRdT*T);
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    const volScalarField& currZoneVol = mesh().lookupObject<volScalarField>(currentZone); 
    // eqn += currZoneVol*coatingVolFrac[countI]*fvm::Sp(gradR,currSp);
    // eqn += currZoneVol*coatingVolFrac[countI]*(R - gradR*currSp);
    
    eqn += currZoneVol*coatWallFrac[countI]*fvm::Sp(dRdT,TSolid);
    eqn += currZoneVol*coatWallFrac[countI]*(Rent - dRdT*TSolid);

    // eqn += currZoneVol*coatWallFrac[countI]*(Rent);
} 
// Note: I add this source as purely explicit
