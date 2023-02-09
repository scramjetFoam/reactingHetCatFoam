// -- update sum of the mass fractions without reference
sumWithoutRef = 0*sumWithoutRef;

// -- calculate the diffusive flux
// -- update sum of the diffusive fluxes
// volVectorField mDOld = mD;
// sumMDi = sumMDi * 0;

// volVectorField UDOld = UD;
UD = UD * 0;

// -- solve the concentration equation for all but one chemical species
forAll(specieNames,lstI)
{
    // -- get current name of the specie that is solved for
    const word currSpName = specieNames[lstI];                            	

    Info << "\nWorking on the specie: " << currSpName << endl;

    // -- mass fraction of the current specie              	
    volScalarField& currSpMass	= MassFrList[lstI];                	
    volScalarField& currSpMoles	= MFrList[lstI];                	
    
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
            // fvm::div(fvc::interpolate(UDOld * rho) & mesh.Sf(), currSpMass)
            // fvm::div(phi + phiCons + phiCorr, currSpMass) 
            ==
            // -- specie mass diffusion                           
            // fvm::laplacian(DEff*rho, currSpMass) 
            // // + fvc::laplacian(DEff*rho/Mg*currSpMass, Mg)   
            // + fvm::div(linearInterpolate(rho*DEff/Mg*fvc::grad(Mg)) & mesh.Sf(), currSpMass)  

            fvm::laplacian(DEff*rho, currSpMass) 
            // + fvc::laplacian(DEff*currSpMass, rho)   
            + fvm::div(fvc::interpolate(DEff*rho/Mg*fvc::grad(Mg)) & mesh.Sf(), currSpMass)    
            // + fvm::div(linearInterpolate(DEff*fvc::grad(rho)) & mesh.Sf(), currSpMass)    

            // fvm::laplacian(DEff*rho, currSpMass) 
            //  fvc::laplacian(DEff*rho/Mg*molMList[lstI], currSpMoles) 
            // -- source
            + fvModels.source(dummyConcEqn, currSpMass)
            //  fvModels.source(dummyConcEqn, currSpMass)
        );
            
        // fvOptions.constrain(concEq);
        fvConstraints.constrain(concEq);

        concEq.relax();

        concEq.solve();

        fvConstraints.constrain(currSpMass);

        currSpMass.correctBoundaryConditions();

        sumWithoutRef += currSpMass;
        // mD += rho * U * currSpMass;
        // sumMDi += fvc::laplacian( DEff * rho / Mg, currSpMass * Mg);
        // mD -= (DEff*rho/Mg*fvc::grad(currSpMass*Mg));
        // mD -= (DEff*rho/Mg*fvc::grad(currSpMoles)*molMList[lstI]);

        UD -= (DEff/Mg*fvc::grad(currSpMass*Mg));
        // mD -= (DEff/Mg*fvc::grad(currSpMass*Mg));
        // sumWiMi += currSpMass / molMList[lstI];
        // Info << "mD = " << max(mD).value() << endl;

    }
    else
    {
        // mD += rho * U * currSpMass;
        
        currSpMass = 1 - sumWithoutRef;

        UD = UD / currSpMass;

        // mD = mD / currSpMass;
        
        // mD = mD / currSpMass;
        
        // UD -= (DEff/Mg*fvc::grad(currSpMass*Mg));
        // Info << "UD = " << max(UD).value() << endl;
        // currSpMass = (1 / Mg - sumWiMi) * molMList[lstI];
        // sumJMDiff += (DEff*rho*fvc::grad(currSpMass));
    }
}

// surfaceScalarField phiOld = (linearInterpolate(rho*U) & mesh.Sf());

for (label corMg = 0; corMg < 1; corMg++)
{
    volScalarField sumWiMi = 0/Mg;
    forAll(specieNames,lstI)
    {
        sumWiMi += MassFrList[lstI] / molMList[lstI];
    }
    Mg = 1/sumWiMi;

    // -- calculate the molar fractions, update UD from overall diffusive flux
    sumWithoutRef = sumWithoutRef * 0;
    forAll(specieNames,lstI)
    {
        MFrList[lstI] = MassFrList[lstI] / molMList[lstI] * Mg;
        // if (lstI != specieNames.size() - 1){
        //     UD -= (DEff/Mg*fvc::grad(MassFrList[lstI]*Mg));
        // }
        // // if (lstI != specieNames.size() - 1){
        // //     // mD -= (DEff*rho/Mg*fvc::grad(MassFrList[lstI]*Mg));
        // //     // mDF -= fvc::interpolate(DEff*rho/Mg * fvc::grad(MassFrList[lstI]*Mg)) & mesh.Sf();
        // //     // UDF -= fvc::interpolate(DEff/Mg) * fvc::snGrad(MassFrList[lstI]*Mg);
        // //     UD -= (DEff/Mg*fvc::grad(MassFrList[lstI]*Mg));
        // // }
        // //     // dimensionedScalar D(transportProperties.subDict(specieNames[lstI]).lookup("D"));//R constant
        // //     // DEff = 0*DEff + D;
        // // // }
        // else
        // {
        // //     // mD = mD / MassFrList[lstI];
        //     UD = UD / MassFrList[lstI];
        // //     // UDF = UDF / fvc::interpolate(MassFrList[lstI]);
        // }
    }

    // UDFS = UDFSOld * UDF;

    // surfaceScalarField rhoS = fvc::interpolate(rho);

    // forAll(UDFS, faceI)
    // {
        // UDFS[faceI][0] = UDF[faceI];
    // }

    // Info << max(UDFS.boundaryField()) << endl;

    // UDFS = mDF / fvc::interpolate(rho);
    
    // forAll(specieNames,lstI)
    // {
    //     sumJMDiff += (DEff*rho/Mg*fvc::grad(MassFrList[lstI]*Mg));
    // }

    // -- calculate the gas molar mass 
    // sumWiMi = 0/Mg;
    // Mg = Mg*0;

    // forAll(specieNames,lstI)
    // {
    //     // sumWiMi += MassFrList[lstI] / molMList[lstI];
    //     Mg += MFrList[lstI] * molMList[lstI];
    // }

    // sumWithoutRef = sumWithoutRef*0;
    // forAll(specieNames,lstI)
    // {
    //     // if (lstI != specieNames.size() - 1){
    //         MassFrList[lstI] = MFrList[lstI] * molMList[lstI] / Mg;
    //     // }
    //     // else
    //     // {
    //     //     MassFrList[lstI] = 1 - sumWithoutRef;
    //     // }
    // }

    // sumWiMi = 0/Mg;
    // forAll(specieNames,lstI)
    // {
    //     sumWiMi += MassFrList[lstI] / molMList[lstI];
    // }
    // Mg = 1/sumWiMi;
}
// volScalarField rhoStar = p / univR / thermo.T() * (Mg - Mgconst);

// phiCons = ((linearInterpolate(rhoStar * U) & mesh.Sf()));
// volVectorField UStar = - (fvc::div(phi) + fvc::div(phiCons) + fvc::div(((linearInterpolate(UStar) & mesh.Sf()))) * rho) / fvc::div(((linearInterpolate(rho)) & mesh.Sf()));
// phiCorr = fvc::div(linearInterpolate(rho * UStar) & mesh.Sf());
// Mg = 1/sumWiMi;

// volScalarField rhoOld = rho;

// rho = p/univR/thermo.T()*Mg;
// thermo.rho() = rho;

// p = p * rho / rhoOld;

// rho = rhoOld;


// phi = linearInterpolate(rho*U) & mesh.Sf();

// UD = sumJMDiff / Mg;
// phiD = linearInterpolate(rho*UD) & mesh.Sf();

// -- inform about mass and molar fractions
if (debug > 0)
{
    forAll(specieNames,lstI)
    {
            Info << "Min y[" << specieNames[lstI] << "]: " << min(MFrList[lstI]).value() << ", max y[" << specieNames[lstI] << "]: " << max(MFrList[lstI]).value() << endl;
            Info << "Min w[" << specieNames[lstI] << "]: " << min(MassFrList[lstI]).value() << ", max w[" << specieNames[lstI] << "]: " << max(MassFrList[lstI]).value() << endl;
    }
}
