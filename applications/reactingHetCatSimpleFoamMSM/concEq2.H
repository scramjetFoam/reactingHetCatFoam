// -- update sum of the mass fractions without reference
sumWithoutRef = 0 * sumWithoutRef;

// -- calculate the diffusive flux
// -- update sum of the diffusive velocity for boundary condition (all but the reference)
UD = UD * 0;

// -- update Stefan-Maxwell diffusion coeficients
// -- for all cells
forAll(rho,cellI)
{
    // Info << "Start preparing inv(B)" << endl;
    // -- in each cell prepare matrix to invert
    scalarSquareMatrix B(dimMat);
    // Info << B<<endl;
    // B(0,0) =3;
    // B(1,0) =1;
    // B(0,1) =2;
    // Info << "B" << B <<endl;


    // NOTETH :: important!!! ---> porosity, tortuosity and knudsen not included yet

    // -- creation of the B matrix -- see 1993_TaylorKrisna_MultiComp...
    for (int i = 0; i < dimMat; i++)
    {
        for (int j = 0; j < dimMat; j++)
        {
            // -- diagonal
            if (j == i)
            {
                scalar suma = 0;
                for (int k = 0; k < n; k++)
                {
                    if (k != i)
                    {
                        // suma += (MFrList[k][cellI]/(DijList[i*(n)+k]*Foam::pow(thermo.T()[cellI],1.75)/thermo.p()[cellI])).value();
                        // suma += (MFrList[k][cellI]/DEffFnc(eps,tort,zoneLst,freeChannel,DijList,thermo,molMList,dp,volDiffLst,specieNames, transportProperties,univR,index,i)[cellI]).value();
                        DEffFnc(eps,tort,zoneLst,freeChannel,DEff,DKn,DijList,thermo,molMList,dp,volDiffLst,specieNames,transportProperties,univR,i,k,n);
                        suma += (MFrList[k][cellI]/DEff[cellI]);
                        // suma += (MFrList[k][cellI]).value();
                        // suma += DEffFnc(eps,tort,zoneLst,freeChannel,DijList,thermo,molMList,dp,volDiffLst, specieNames,transportProperties,univR,index,i)[cellI].value();
                    }
                }
                DEffFnc(eps,tort,zoneLst,freeChannel,DEff,DKn,DijList,thermo,molMList,dp,volDiffLst,specieNames,transportProperties,univR,i,n-1,n);
                B(i,i) = (MFrList[i][cellI]/(DEff[cellI])) + suma;
            }

            // -- non-diagonal
            else
            {
                // B(i,j) = B(j,i);
                // B(i,j)
                // int index = i*(n)+j; 
                DEffFnc(eps,tort,zoneLst,freeChannel,DEff,DKn,DijList,thermo,molMList,dp,volDiffLst,specieNames,transportProperties,univR,i,j,n);
                B(i,j) = -(MFrList[i][cellI] * (1.0/(DEff[cellI])));
                // index = i*(n)+(n-1); 
                DEffFnc(eps,tort,zoneLst,freeChannel,DEff,DKn,DijList,thermo,molMList,dp,volDiffLst,specieNames,transportProperties,univR,i,n-1,n);
                B(i,j) += MFrList[i][cellI] * 1.0/(DEff[cellI]);
            }
            // Info <<i << " " <<j << " " << n << " " << B(i,j) << endl;
        }
    } 
    // Info<<B<<endl;

    // -- inversion of the matrix -- maybe too much complicated -- rethink later
    scalarSquareMatrix Binv = scalarSquareMatrix(SVDinv(B));

    // Info << "B" << B << "Binv" << Binv <<endl;
    // Info << Binv(1,1) <<endl;

    // scalarSquareMatrix Binv;
    // B.decompose();
    // B.inv(Binv);

    for (int i = 0; i < dimMat; i++)
    {
        for (int j = 0; j < dimMat; j++)
        {
            // if (j >= i)
            // {
                // if (Binv(i,j) < 1e-10){
                //     Info <<Binv<< i<<j<<endl;
                // }
                BinvList[i*(dimMat)+j][cellI] = Binv (i,j);
            // }
        }
    }
}

for (int i = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        BinvList[i*(dimMat)+j].correctBoundaryConditions();
    }
}


// -- print invB matrix
Info << "min(invBMat)" <<endl;
for (int i  = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        // if (j >= i)
        // {
            Info << " " << min(BinvList[i*(dimMat)+j]).value();
        // }
        // else
        // {
        //      Info << " " << min(BinvList[j*(dimMat)+i]).value();
        // }
    }
    Info <<endl;
}

Info << "max(invBMat)" <<endl;
for (int i  = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        // if (j >= i)
        // {
            Info << " " << max(BinvList[i*(dimMat)+j]).value();
        // }
        // else
        // {
        //      Info << " " << max(BinvList[j*(dimMat)+i]).value();
        // }
    }
    Info <<endl;
}

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
        // volScalarField DifMix
        // (
        //     DABLst[lstI]*Foam::pow(thermo.T(),1.75)/thermo.p()
        // );

        // // -- calculate the Knudsen difussivity
        // DKn = 0*DifMix;                    
        // DKn += freeChannel*DifMix;         // free-channel (not used, but cant be 0)
        // forAll(zoneLst,countI)
        // {
        //     DKn += zoneLst[countI]*(1.0/3.0*sqrt((thermo.T()*8*univR)/(3.141592653*molMList[lstI])))*dp[countI];
        // }

        // // -- if diffusion coeficient is set in transportProperties dont compute it
        // if (transportProperties.subDict(currSpName).found("D"))
        // {
        //     dimensionedScalar D(transportProperties.subDict(currSpName).lookup("D"));//R constant
        //     DifMix = 0*DifMix + D;
        // }

        // // -- effective diffusion coeficient
        // DEff = freeChannel*DifMix;         // free-channel
        // forAll(zoneLst,countI)
        // {
        //     DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
        // } 

        // -- inform about calculated properties
        if (debug > 0){
            Info << "Min rho: " << min(rho).value() << ", max rho: " << max(rho).value() << endl;
            // Info << "Min DEff[" << specieNames[lstI] << "]: " << min(DEff).value() << ", max DEff["<<specieNames[lstI]<<"]: "<< max(DEff).value() << endl;
            Info << "Min binv[" << specieNames[lstI] << "]: " << min(BinvList[lstI*dimMat+lstI]).value() << ", max binv["<<specieNames[lstI]<<"]: "<< max(BinvList[lstI*dimMat+lstI]).value() << endl;
            // Info << "Min DFree[" << specieNames[lstI] << "]: " << min(DifMix).value() << ", max DFree["<<specieNames[lstI]<<"]: "<< max(DifMix).value() << endl;
        }
        else if (debug == 2)
        {
            Info << "Min Mg: " << min(Mg).value() << ", max Mg: " << max(Mg).value() <<" Mgconst "<<molMList[0]<< endl;
        }

        // -- calculate contribution of other gases to the equation
        volVectorField explJi (molMList[lstI] * rho / Mg * BinvList[0] * fvc::grad(MFrList[0]) * 0);
        volVectorField explJiF = explJi * 0;
        for (int j = 0; j < dimMat; j++)
        {
            if (j != lstI)
            {
                Info<<"Adding contribution from " << specieNames[j] << endl;
                explJi += molMList[lstI] * rho / Mg * BinvList[lstI*(dimMat)+j] * fvc::grad(MFrList[j]);
                explJiF -= molMList[lstI] * BinvList[j*(dimMat)+j] * rho / Mg * fvc::grad(MFrList[j]) * (1 - ((MassFrList[n-1] * MFrList[j]) / (MFrList[n-1] * MassFrList[j])));  
            }
            // else if (j < lstI)
            // {
            //     Info<<"Adding contribution from " << specieNames[j] << endl;
            //     explJi += molMList[lstI] * rho / Mg * BinvList[j*(dimMat)+lstI] * fvc::grad(MFrList[j]);
            // }
        }

        fvScalarMatrix concEq
        (
            // -- specie mass convection
            fvm::div(phi, currSpMass) 

            ==

            // -- specie mass diffusion                          
            // fvm::laplacian(DEff*rho, currSpMass) 
            // + fvm::div(fvc::interpolate(DEff*rho/Mg*fvc::grad(Mg)) & mesh.Sf(), currSpMass)  
            // -- Stefan Maxwell diffusion
            // fvm::laplacian(rho*BinvList[lstI*(dimMat)+lstI], currSpMass)
            // + fvm::div(fvc::interpolate(rho * BinvList[lstI*(dimMat)+lstI] / Mg * fvc::grad(Mg)) & mesh.Sf(), currSpMass)
            // + fvc::div(explJi) 
            
            // fvc::div(JiList[lstI])  
    	    
            fvm::laplacian(BinvList[lstI*(dimMat)+lstI]              * rho * (1 - currSpMass * (1 - (MassFrList[n-1] * MFrList[lstI]) / (MFrList[n-1] * MassFrList[lstI]))),currSpMass)
            + fvm::div(fvc::interpolate(BinvList[lstI*(dimMat)+lstI] * rho * (1 - currSpMass * (1 - (MassFrList[n-1] * MFrList[lstI]) / (MFrList[n-1] * MassFrList[lstI]))) / Mg * fvc::grad(Mg) ) & mesh.Sf(), currSpMass)
            // + fvc::div(explJiF)  
            + fvm::div(fvc::interpolate(explJiF) & mesh.Sf(), currSpMass)  

            // -- source
            + fvModels.source(dummyConcEqn, currSpMass)
        );

        // for (int j = 0; j < dimMat; j++)
        // {
        //     concEq -= adawda;
        // }

        // fvOptions.constrain(concEq);
        fvConstraints.constrain(concEq);

        concEq.relax();

        concEq.solve();

        fvConstraints.constrain(currSpMass);

        currSpMass.correctBoundaryConditions();

        sumWithoutRef += currSpMass;
        sumWithoutRef.correctBoundaryConditions();

        // UD -= ( DEff / Mg * fvc::grad(currSpMass * Mg)) / molMList[lstI].value();
        // UD -= ( BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(currSpMass * Mg) ) / molMList[lstI].value();
        // UD -= ( BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(MFrList[lstI]) ) * molMList[lstI] + explJi / rho;
        // UD -= ( BinvList[lstI*dimMat+lstI] * fvc::grad(MFrList[lstI]) ) + explJi / rho * Mg / molMList[lstI];
        // UD -= ( BinvList[lstI*dimMat+lstI] * fvc::grad(MFrList[lstI]) / (Mg/dimensionedScalar("jedna",dimMass/dimMoles,1)) ) + explJi / rho / molMList[lstI] *dimensionedScalar("jedna",dimMass/dimMoles,1);
        // UD -= ( BinvList[lstI*dimMat+lstI] * fvc::grad(MFrList[lstI]) ) + explJi / rho * Mg / molMList[lstI];
        // UD -= ( BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(currSpMass * Mg) + explJi / rho) / molMList[lstI].value();
        // UD -= ( BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(currSpMass * Mg) + explJi / rho);
        // UD.correctBoundaryConditions();

        // JiList[lstI] = ( rho * U * currSpMass + rho * BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(currSpMass * Mg) + explJi ) / molMList[lstI];
        // JiList[lstI] = ( rho * U * currSpMass - rho * BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(currSpMass * Mg) - explJi ) / molMList[lstI];
        // UDS -= fvc::interpolate( DEff * fvc::grad(currSpMass * Mg) / Mg );
        // JUiList[lstI] = ( rho * U * currSpMass - rho * BinvList[lstI*dimMat+lstI] / Mg * fvc::grad(currSpMass * Mg) - explJi );
        JUiList[lstI] = -( BinvList[lstI*(dimMat)+lstI] * rho * (1 - currSpMass * (1 - (MassFrList[n-1] * MFrList[lstI]) / (MFrList[n-1] * MassFrList[lstI]))) / Mg * fvc::grad(Mg * currSpMass)  + explJiF * currSpMass);
        JUiList[lstI].correctBoundaryConditions();
    }
    else
    {
        currSpMass = 1 - sumWithoutRef;
        currSpMass.correctBoundaryConditions();

        // UD = UD / currSpMass;
        // UD = UD / currSpMass * molMList[lstI].value();
        // UD = UD / MFrList[lstI];
        // UD = UD / MFrList[lstI] ;
        volVectorField J = JiList[0] * 0;
        for (int k = 0; k < n-1; k++)
        {
            J += JUiList[k] ;
        }
        UD = J / currSpMass / rho;
        UD.correctBoundaryConditions();
    }
}

// volScalarField MgOld = Mg;

sumWiMi = 0 * sumWiMi;
forAll(specieNames,lstI)
{
    sumWiMi += MassFrList[lstI] / molMList[lstI];
    sumWiMi.correctBoundaryConditions();
}
Mg = 1/sumWiMi;
Mg.correctBoundaryConditions();

for (int i = 0; i < n-1; i++)
{
    JiList[i] = 0*JiList[i];
    for (int k = 0; k < dimMat; k++)
    {
        // -- diagonal
        if (k == i)
        {
            JiList[i] += JUiList[k];
        }
        JiList[i] -= JUiList[k]*MassFrList[i]*(1 - (MassFrList[n-1] * MFrList[k]) / (MFrList[n-1] * MassFrList[k]));
    }
}

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
