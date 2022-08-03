forAll(specieNames,lstI)
{
      word currSpName = specieNames[lstI];                            //get current specie name
    
      // read the corresponding row of stoechiometric matrix
      List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

      forAll(nuLst,rnI)                                                       //for all reactions
      {    
            if (mag(nuLst[rnI]) > SMALL and activeReacts[rnI])                  //tango only in active reactions
            {
                  Info << currSpName << " is active in reaction " << rnI << endl;

                  char rName [50];
                  label nChars = sprintf(rName,"reaction%02d",rnI);
                  
                  //load data from transport propertis
                  dimensionedScalar k0 = transportProperties.subDict("reactProps").subDict(rName).lookup("k0");
                  dimensionedScalar EA = transportProperties.subDict("reactProps").subDict(rName).lookup("EA");
                  dimensionedScalar kin = transportProperties.subDict("reactProps").subDict(rName).lookup("kin");
                  dimensionedScalar Ein = transportProperties.subDict("reactProps").subDict(rName).lookup("Ein");

                  // load the standard reaction heat at temperature 298K
                  dimensionedScalar sHr298 = transportProperties.subDict("reactProps").subDict(rName).lookup("sHr");
                  List<scalar> nuLstReact(transportProperties.subDict("reactProps").subDict(rName).lookup("nuReact"));
                  volScalarField sHr = sHatT(sHr298, ALst,BLst,CLst, DLst, nuLstReact, specieNames, mesh_);

                  Info << min(sHr)<<max(sHr)<<endl;

                  // calculate reaction rate constant 
                  volScalarField k(k0*Foam::exp(-EA/(univR*T)));
                  volScalarField Kin(kin*Foam::exp(Ein/(T)));

                  //concentration fields
                  const volScalarField& CO = mesh_.lookupObject<volScalarField>("CO"); 
                  const volScalarField& O2 = mesh_.lookupObject<volScalarField>("O2"); 

                  //reaction rate
                  R = k*CO*O2*Foam::pow(1+Kin*CO,-2);   

                  //enthaphy source
                  rQ -= isCoat*coatingVolFrac*R*sHr;


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

Info <<"rQ - enthalphy source: "<< min(rQ) <<max(rQ)<< endl;         
eqn += rQ;
// Note: I add this source as purely explicit
