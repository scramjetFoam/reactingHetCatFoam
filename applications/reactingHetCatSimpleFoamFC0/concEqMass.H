// -- update sum of the molar fractions without reference
sumWithoutRef = 0*sumWithoutRef;
// if (diffusion == "FickConv")
// {
//     sumJDif = 0*sumJDif;
// }

volScalarField MgNew = 0*Mg;
PtrList<volScalarField> MassFrList(specieNames.size());

// -- solve the concentration equation for all chemical species
forAll(specieNames,lstI)
{
    // -- get current name of the specie that is solved for
    word currSpName = specieNames[lstI];                            	

    Info<<"\nWorking on the specie: "<<currSpName<<endl;

    // volScalarField& currSp	= MFrList[lstI];                	//current concentration field
    volScalarField currSpMass	= MFrList[lstI] * molMList[lstI] / Mg;                	//mass fraction
    currSpMass.correctBoundaryConditions();
    
    // if (lstI != specieNames.size()-1)
    // {
        // -- DIFUSION COEFICIENT CALCULATION 
        // -- difussivity of the specie in the reference gas see 22_Hlavaty_MSc

        // -- Diffusion coeficient of the mixture -- default reference gas diffusion + add temperature and pressure dependence (see Fuller)
        volScalarField DifMix
        (
            DABLst[lstI]*Foam::pow(thermo.T(),1.75)/thermo.p()
        );

        DKn = 0*DifMix;                    // reset knudsen
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
        if (debug > 0){
            Info << "Min rho: " << min(rho).value() << ", max rho: " << max(rho).value() << endl;
            Info << "Min DEff[" << specieNames[lstI] << "]: " << min(DEff).value() << ", max DEff["<<specieNames[lstI]<<"]: "<< max(DEff).value() << endl;
            Info << "Min DFree[" << specieNames[lstI] << "]: " << min(DifMix).value() << ", max DFree["<<specieNames[lstI]<<"]: "<< max(DifMix).value() << endl;
            // Info << "Min y[" << specieNames[lstI] << "]: " << min(currSp).value() << ", max y[" << specieNames[lstI] << "]: " << max(currSp).value() << endl;
        }
        if (debug <= 2 && debug != 0)
        {
            Info << "Min Mg: " << min(Mg).value() << ", max Mg: " << max(Mg).value() <<" Mgconst "<<molMList[0]<< endl;
        }

        fvScalarMatrix concEq
        (
            // -- convection
            fvm::div(phi, currSpMass)
            ==
            // -- diffusion
            fvm::laplacian(DEff*rho, currSpMass)                              
            // + fvc::laplacian(DEff*currSp, cTot)       
            // -- source
            + fvModels.source(dummyConcEqn*molMList[lstI], currSpMass)
        );

        if (diffusion == "FickConv")
        {
            // concEq = fvScalarMatrix
            //     (
            //         // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            //         fvm::div(phi/fvc::interpolate(Mg), currSp)
            //         // fvm::div(phiD/fvc::interpolate(Mg), currSp)
            //         ==
            //         fvm::laplacian(DEff*rho/Mg, currSp)                               
            //         // + fvc::laplacian(DEff*currSp, thermo.p()/univR/thermo.T())     
            //         + fvModels.source(dummyConcEqn, currSp)
            //     );
        }
            
        // fvOptions.constrain(concEq);
        fvConstraints.constrain(concEq);

        concEq.relax();

        concEq.solve();

        fvConstraints.constrain(currSpMass);

        currSpMass.correctBoundaryConditions();

        MgNew += molMList[lstI] / currSpMass;

        MassFrList[lstI] = currSpMass;
        
        // forAll(currSp,cellI)
        // { 
        //     currSp[cellI] = min(max(currSp[cellI],SMALL),1-SMALL); 
        //     // if(currSp[cellI] < 0){
        //         // Info <<currSp[cellI]<<endl;
        //     // } 
        // }

        // sumWithoutRef += currSp;
        // if (diffusion == "FickConv")
        // {
        //     sumJDif += (DEff*(fvc::grad(thermo.p()/univR/thermo.T()*currSp)));
        // }
    // }
    // else{
    //     // currSp = dimensionedScalar("jedna",dimless,1) - sumWithoutRef;
    //     sumWithoutRef += currSp;
    // }
}

Mg = MgNew;
// forAll(specieNames,lstI)
// {
//     Mg += MFrList[lstI] * molMList[lstI];
// }

// Info << "Min y[sumWithoutRef]: " << min(sumWithoutRef).value() << ", max y[sumWithoutRef]: " << max(sumWithoutRef).value() << endl;
forAll(specieNames,lstI)
{
    // MFrList[lstI] = MFrList[lstI]/;
    MFrList[lstI] = MassFrList[lstI] * Mg / molMList[lstI];
    Info << "Min y[" << specieNames[lstI] << "]: " << min(MFrList[lstI]).value() << ", max y[" << specieNames[lstI] << "]: " << max(MFrList[lstI]).value() << endl;
}

// thermo.p() = p*sumWithoutRef;
        
        // currSp = (1-sumWithoutRef)/(currSp);

if (diffusion == "FickConv")
{
    // -- update molar mass of the gas and velocity from diffusion
    // UD = -sumJDif/MFrList[MFrList.size()-1]*Mg/rho;

    // if (debug != 0)
    // {
    //         Info << "Min UD: " << min(UD).value() << ", max UD: " << max(UD).value() << endl;
    // }

    // UD.correctBoundaryConditions();

    // // -- update phi
    // // phi = linearInterpolate(rho*U) & mesh.Sf();
    // phiD = linearInterpolate(rho*UD) & mesh.Sf();
}


