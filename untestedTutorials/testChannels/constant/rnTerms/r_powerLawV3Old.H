// read the corresponding row of stoechiometric matrix
List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

// load the current specie
const volScalarField&  currSp = mesh_.lookupObject<volScalarField>(currSpName);

// load data for numerics form transportProperties
scalar  minCForGradND(readScalar(transportProperties.subDict(currSpName).lookup("minForGrad")));
dimensionedScalar minCForGrad("minCForGrad",currSp.dimensions(),minCForGradND);

volScalarField cForGrad(max(currSp,minCForGrad));


forAll(nuLst,rnI)                                                       //for all reactions
{    
    if (mag(nuLst[rnI]) > SMALL and activeReacts[rnI])                  //tango only in active reactions
    {
        //~ Pout << currSpName << " is active in reaction " << rnI;
        //~ Pout << " => loading kinetics parameters (k0 and EA)" << endl;
        
        char rName [50];
        label nChars = sprintf(rName,"reaction%02d",rnI);
        
        dimensionedScalar k0  = transportProperties.subDict("reactProps").subDict(rName).lookup("k0");
        dimensionedScalar EA  = transportProperties.subDict("reactProps").subDict(rName).lookup("EA");
	dimensionedScalar kin  = transportProperties.subDict("reactProps").subDict(rName).lookup("kin");
        dimensionedScalar Ein  = transportProperties.subDict("reactProps").subDict(rName).lookup("Ein");
        
	Info << currSpName << " is active in reaction " << rnI<<" "<<nChars;
	//~ Pout << currSpName << " is active in reaction " << rnI;
	Info << " => loading kinetics parameters (k0 = "<<k0.value()<<", EA = "<< EA.value()<<", kin = "<<kin.value()<<" and Ein = "<< Ein.value()<<")" << endl;
        //counting exponent for total concetration
        //NOTE: I start with -1.0 because of units of kinetic constants 
        //      in transport properties (maybe there can be problem with
        //      non-integer alpha)
        //~ float cTotExp = -1.0;
        //~ forAll(specieNames,countI){
            //~ //this can be probably done nicer, but I used the code above  
            //~ List<scalar> alphaLstCurSp(transportProperties.subDict(specieNames[countI]).lookup("alphaVec"));    
            //~ cTotExp += alphaLstCurSp[rnI];
        //~ }
	    
	//kinetics from article 3D reconstruction ...
	//~ k0 = k0*44.8*10/cTemp.value();	
	//~ dimensionedScalar k = k0*Foam::exp(-EA/(univR*cTemp));
	//~ dimensionedScalar Kin = kin*Foam::exp(Ein/(cTemp));
	volScalarField k(k0*Foam::exp(-EA/(univR*T)));
	volScalarField Kin(kin*Foam::exp(Ein/(T)));
    
	
        // calculating real kinetic constant of reaction with arrhenius kinetics and scaling with concentration for right dimension

	const volScalarField& CO = mesh_.lookupObject<volScalarField>("CO"); 
	const volScalarField& O2 = mesh_.lookupObject<volScalarField>("O2"); 

	
	//creating dimensionScalar with value 1 and dimensions of the yCO
	//~ dimensionedScalar oneForSum = dimensionedScalar("oneForSum",dimless/dimTime,1.0);

        //volScalarField kReact = k*CO*O2*(Foam::pow((1+Kin*CO),-2));
        
    
	R	= k*nuLst[rnI]*CO*O2*Foam::pow(1+Kin*CO,-2)*Foam::pow(cTot,-1);
	if (currSpName == "CO"){
		Info << "CO: nu="<<nuLst[rnI];
		//gradR	= nuLst[rnI]*((-2.0*Kin*k*CO*O2*(Foam::pow(Kin*CO+1,-3)))+(k*O2*(Foam::pow(Kin*CO+1,-2))))*Foam::pow(cTot,-1);
		gradR = nuLst[rnI]*(k*O2*(1-Kin*CO)*Foam::pow(1+Kin*CO,-3))*Foam::pow(cTot,-1);
		Info << "gradRmin= "<<min(gradR)<<"gradRmax"<<max(gradR)<<endl; 
   	}
	else
	{
		Info << "O2: nu="<<nuLst[rnI];
		gradR   = nuLst[rnI]*(k*CO*Foam::pow(1+Kin*CO,-2))*Foam::pow(cTot,-1);
		Info << "gradRmin= "<<min(gradR)<<"gradRmax"<<max(gradR)<<endl; 
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
            
eqn += isCoat*coatingVolFrac*fvm::Sp(gradR,currSp);
eqn += (isCoat*coatingVolFrac*R - isCoat*coatingVolFrac*gradR*currSp);
//eqn += isCoat*coatingVolFrac*R

Info<<"min: "<< min(currSp);
Info<<"max: "<< max(currSp) <<endl;
 
// Note: Taylor expansion of kinetics
