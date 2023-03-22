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
                        suma += (MFrList[k][cellI]/(DijList[i*n+k]*Foam::pow(thermo.T()[cellI],1.75)/thermo.p()[cellI])).value();
                    }
                }
                B(i,i) = (MFrList[i][cellI]/(DijList[i*n+(n-1)]*Foam::pow(thermo.T()[cellI],1.75)/thermo.p()[cellI])).value() + suma;
            }

            // -- uppper triangular
            else if (j > i)
            {
                B(i,j) = - (MFrList[i][cellI] * (1/(DijList[i*n+j]*Foam::pow(thermo.T()[cellI],1.75)/thermo.p()[cellI]) - 1/(DijList[i*n+(n-1)]*Foam::pow(thermo.T()[cellI],1.75)/thermo.p()[cellI]))).value();
            }

            // -- lower triangular
            else
            {
                B(i,j) = B(j,i);
            }
        }
    } 

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
            if (j >= i)
            {
                BinvList[i*dimMat+j][cellI] = Binv (i,j);
            }
        }
    }
}

// -- print invB matrix
Info << "min(invBMat)" <<endl;
for (int i  = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        if (j >= i)
        {
            Info << " " << min(BinvList[i*dimMat+j]).value();
        }
        else
        {
             Info << " " << min(BinvList[j*dimMat+i]).value();
        }
    }
    Info <<endl;
}

Info << "max(invBMat)" <<endl;
for (int i  = 0; i < dimMat; i++)
{
    for (int j = 0; j < dimMat; j++)
    {
        if (j >= i)
        {
            Info << " " << max(BinvList[i*dimMat+j]).value();
        }
        else
        {
             Info << " " << max(BinvList[j*dimMat+i]).value();
        }
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

        // -- if diffusion coeficient is set in transportProperties dont compute it
        if (transportProperties.subDict(currSpName).found("D"))
        {
            dimensionedScalar D(transportProperties.subDict(currSpName).lookup("D"));//R constant
            DifMix = 0*DifMix + D;
        }

        // -- effective diffusion coeficient
        DEff = freeChannel*DifMix;         // free-channel
        forAll(zoneLst,countI)
        {
            DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
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

        UD -= ( DEff / Mg * fvc::grad(currSpMass * Mg)) / molMList[lstI].value();
        // UDS -= fvc::interpolate( DEff * fvc::grad(currSpMass * Mg) / Mg );
    }
    else
    {
        currSpMass = 1 - sumWithoutRef;

        UD = UD / currSpMass * molMList[lstI].value();
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
