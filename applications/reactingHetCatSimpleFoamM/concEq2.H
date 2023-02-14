// -- update sum of the mass fractions without reference
sumWithoutRef = 0 * sumWithoutRef;

// -- calculate the diffusive flux
// -- update sum of the diffusive velocity for boundary condition (all but the reference)
UD = UD * 0;

// -- solve the concentration equation for all but one chemical species
forAll(specieNames,lstI)
{
    // -- get current name of the specie that is solved for
    const word currSpName = specieNames[lstI];                            	

    Info << "\nWorking on the specie: " << currSpName << endl;

    // -- mass fraction of the current specie              	
    volScalarField& currSpMass	= MassFrList[lstI];             	
    
    // -- if not the last specie
    if (lstI != specieNames.size()-1)
    {
        // -- DIFUSION COEFICIENT CALCULATION 
        // -- difussivity of the specie in the reference gas see 22_Hlavaty_MSc

        // -- Diffusion coeficient of the mixture -- default reference gas diffusion + add temperature and pressure dependence (see Fuller)
        volScalarField DifMix
        (
            DABLst[lstI]*Foam::pow(thermo.T(),1.75)/thermo.p()
        );

        // -- calculate the Knudsen difussivity
        DKn = 0*DifMix;                    
        DKn += freeChannel*DifMix;         // free-channel (not used, but cant be 0)
        forAll(zoneLst,countI)
        {
            DKn += zoneLst[countI]*(1.0/3.0*sqrt((thermo.T()*8*univR)/(3.141592653*molMList[lstI])))*dp[countI];
        }

        // -- effective diffusion coeficient
        DEff = freeChannel*DifMix;         // free-channel
        forAll(zoneLst,countI)
        {
            DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
        } 

        // -- if diffusion coeficient is set in transportProperties dont compute it
        if (transportProperties.subDict(currSpName).found("D"))
        {
            dimensionedScalar D(transportProperties.subDict(currSpName).lookup("D"));//R constant
            DEff = 0*DEff + D;
        }

        // -- inform about calculated properties
        if (debug > 0){
            Info << "Min rho: " << min(rho).value() << ", max rho: " << max(rho).value() << endl;
            Info << "Min DEff[" << specieNames[lstI] << "]: " << min(DEff).value() << ", max DEff["<<specieNames[lstI]<<"]: "<< max(DEff).value() << endl;
            Info << "Min DFree[" << specieNames[lstI] << "]: " << min(DifMix).value() << ", max DFree["<<specieNames[lstI]<<"]: "<< max(DifMix).value() << endl;
        }
        else if (debug == 2)
        {
            Info << "Min Mg: " << min(Mg).value() << ", max Mg: " << max(Mg).value() <<" Mgconst "<<molMList[0]<< endl;
        }

        fvScalarMatrix concEq
        (
            // -- specie mass convection
            fvm::div(phi, currSpMass) 

            ==

            // -- specie mass diffusion                          
            fvm::laplacian(DEff*rho, currSpMass) 
            + fvm::div(fvc::interpolate(DEff*rho/Mg*fvc::grad(Mg)) & mesh.Sf(), currSpMass)    

            // -- source
            + fvModels.source(dummyConcEqn, currSpMass)
        );
            
        // fvOptions.constrain(concEq);
        fvConstraints.constrain(concEq);

        concEq.relax();

        concEq.solve();

        fvConstraints.constrain(currSpMass);

        currSpMass.correctBoundaryConditions();

        sumWithoutRef += currSpMass;

        UD -= ( DEff / Mg * fvc::grad(currSpMass * Mg));
        // UDS -= fvc::interpolate( DEff * fvc::grad(currSpMass * Mg) / Mg );
    }
    else
    {
        currSpMass = 1 - sumWithoutRef;

        UD = UD / currSpMass;
    }
}

// volScalarField MgOld = Mg;

sumWiMi = 0 * sumWiMi;
forAll(specieNames,lstI)
{
    sumWiMi += MassFrList[lstI] / molMList[lstI];
}
Mg = 1/sumWiMi;

// UD = UD * MgOld / Mg;

// -- calculate the molar fractions, update UD from overall diffusive flux
forAll(specieNames,lstI)
{
    MFrList[lstI] = MassFrList[lstI] / molMList[lstI] * Mg;
}

// -- inform about mass and molar fractions
if (debug > 0)
{
    forAll(specieNames,lstI)
    {
            Info << "Min y[" << specieNames[lstI] << "]: " << min(MFrList[lstI]).value() << ", max y[" << specieNames[lstI] << "]: " << max(MFrList[lstI]).value() << endl;
            Info << "Min w[" << specieNames[lstI] << "]: " << min(MassFrList[lstI]).value() << ", max w[" << specieNames[lstI] << "]: " << max(MassFrList[lstI]).value() << endl;
    }
}
