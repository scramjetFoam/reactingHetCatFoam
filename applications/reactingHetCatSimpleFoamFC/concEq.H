// -- update sum of the molar fractions without reference
sumWithoutRef = 0*sumWithoutRef;
if (diffusion == "FickConv")
{
    sumJDif = 0*sumJDif;
}

// -- solve the concentration equation for all chemical species
forAll(specieNames,lstI)
{
    // -- get current name of the specie that is solved for
    word currSpName = specieNames[lstI];                            	

    Info<<"\nWorking on the specie: "<<currSpName<<endl;

    volScalarField& currSp	= MFrList[lstI];                	//current concentration field
    
    if (lstI != specieNames.size()-1)
    {
        // -- DIFUSION COEFICIENT CALCULATION 
        // -- difussivity of the specie in the reference gas see 22_Hlavaty_MSc

        // -- prepare free difusion coeficients of the current specie against all other species 
        // -- NOTE: at the moment we use only diffusion of the specie in the reference gas (nitrogen) -- last member of the list
        PtrList<dimensionedScalar> DABLst(specieNames.size());
        forAll(DABLst,lstII)
        {
            // diffusion coefficients according to Fuller without temperature and pressure dependence
            DABLst.set
            (
                lstII,
                new dimensionedScalar(
                    "DAB",
                    Foam::pow(dimLength,2)/dimTime/Foam::pow(dimTemperature,1.75)*dimPressure,
                    ((0.014)/(Foam::sqrt(2.0/(1e-3/molMList[lstII] + 1e-3/molMList[lstI]))*Foam::pow(Foam::pow(volDiffLst[lstI],1.0/3.0)+Foam::pow(volDiffLst[lstII],1.0/3.0),2))).value()
                )
            );
        }

        // -- Diffusion coeficient of the mixture -- default reference gas diffusion + add temperature and pressure dependence (see Fuller)
        volScalarField DifMix
        (
            DABLst[DABLst.size()-1]*Foam::pow(thermo.T(),1.75)/thermo.p()
        );

        // // -- prepare Knudsen diffusivity
        volScalarField DKn
        (
            IOobject
            (
                "DKn",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("DKn",DifMix.dimensions(),0.0)
        );

        DKn += freeChannel*DifMix;         // free-channel (not used, but cant be 0)
        forAll(zoneLst,countI)
        {
            DKn += zoneLst[countI]*(1.0/3.0*sqrt((thermo.T()*8*univR)/(3.141592653*molMList[lstI])))*dp[countI];
        }

        DEff = freeChannel*DifMix;         // free-channel
        forAll(zoneLst,countI)
        {
            DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
        } 

        // -- if diffusion coeficient is setup previously dont compute it
        if (transportProperties.subDict(currSpName).found("D"))
        {
            dimensionedScalar D(transportProperties.subDict(currSpName).lookup("D"));//R constant
            DEff = 0*DEff + D;
        }

        // -- inform about calculated properties
        if (debug == 1){
            Info << "Min rho: " << min(rho).value() << ", max rho: " << max(rho).value() << endl;
            Info << "Min DEff[" << specieNames[lstI] << "]: " << min(DEff).value() << ", max DEff["<<specieNames[lstI]<<"]: "<< max(DEff).value() << endl;
            Info << "Min DFree[" << specieNames[lstI] << "]: " << min(DifMix).value() << ", max DFree["<<specieNames[lstI]<<"]: "<< max(DifMix).value() << endl;
            Info << "Min y[" << specieNames[lstI] << "]: " << min(currSp).value() << ", max y[" << specieNames[lstI] << "]: " << max(currSp).value() << endl;
        }
        if (debug <= 2 && debug != 0)
        {
            Info << "Min Mg: " << min(Mg).value() << ", max Mg: " << max(Mg).value() <<" Mgconst "<<molMList[0]<< endl;
        }

        fvScalarMatrix concEq
        (
            // -- convection
            fvm::div((linearInterpolate(rho*(U)/Mg) & mesh.Sf()), currSp)
            ==
            // -- diffusion
            fvm::laplacian(DEff*thermo.p()/univR/thermo.T(), currSp)                              
            + fvc::laplacian(DEff*currSp, thermo.p()/univR/thermo.T())       
            // -- source
            + fvModels.source(dummyConcEqn, currSp)
        );
        if (diffusion == "FickConv")
        {
            concEq = fvScalarMatrix
                (
                    fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
                    ==
                    fvm::laplacian(DEff*thermo.p()/univR/thermo.T(), currSp)                              
                    + fvc::laplacian(DEff*currSp, thermo.p()/univR/thermo.T())     
                    + fvModels.source(dummyConcEqn, currSp)
                );
        }
            
        // fvOptions.constrain(concEq);
        fvConstraints.constrain(concEq);

        concEq.relax();

        concEq.solve();

        fvConstraints.constrain(currSp);

        currSp.correctBoundaryConditions();
        
        forAll(currSp,cellI)
        { 
            currSp[cellI] = max(currSp[cellI],0); 
        }

        sumWithoutRef += currSp;
        if (diffusion == "FickConv")
        {
            sumJDif += (DEff*(fvc::grad(thermo.p()/univR/thermo.T()*currSp)));
        }
    }
    else{
        currSp = dimensionedScalar("jedna",dimless,1) - sumWithoutRef;
    }
}

if (diffusion == "FickConv")
{
    // -- update molar mass of the gas and velocity from diffusion
    UD = -sumJDif/MFrList[MFrList.size()-1]*Mg/rho;

    if (debug != 0)
    {
            Info << "Min UD: " << min(UD).value() << ", max UD: " << max(UD).value() << endl;
    }

    // -- update phi
    // phi = linearInterpolate(rho*U) & mesh.Sf();
    phiD = linearInterpolate(rho*UD) & mesh.Sf();

    UD.correctBoundaryConditions();
}

