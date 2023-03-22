// -- update diffusion velocity of the so far species
// jiDiff = 0*jiDiff;
// sumFick = 0*sumFick;

// -- update sum of the molar fractions without reference
sumWithoutRef = 0*sumWithoutRef;
sumJDif = 0*sumJDif;

// -- total molar concentration
volScalarField cTot(rho/Mg);

// -- update Stefan-Maxwell diffusion coeficients
// -- for all cells
forAll(rho,cellI)
{
    // Info << "Start preparing inv(B)" << endl;
    // -- in each cell prepare matrix to invert
    scalarSquareMatrix B(dimMat);

    // NOTETH :: important!!! ---> porosity, tortuosity and knudsen not included yet

    // -- creation of the B matrix
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

// -- solve the concentration equation for all chemical species
forAll(specieNames,i)
{
    // -- get current name of the specie that is solved for
    word currSpName = specieNames[i];                            	

    Info << "\nWorking on the specie: " << currSpName << endl;

    volScalarField& currSp	= MFrList[i];                	//current concentration field

    if (i != specieNames.size()-1)
    {
        // -- prepare Knudsen diffusivity
        // volScalarField DKn
        // (
        //     IOobject
        //     (
        //         "DKn",
        //         runTime.timeName(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("DKn",DifMix.dimensions(),0.0)
        // );

        // DKn += freeChannel*DifMix;         // free-channel (not used, but cant be 0)
        // forAll(zoneLst,countI)
        // {
        //     DKn += zoneLst[countI]*(1.0/3.0*sqrt((thermo.T()*8*univR)/(3.141592653*molMList[lstI])))*dp[countI];
        // }

        // DEff = freeChannel*DifMix;         // free-channel
        // forAll(zoneLst,countI)
        // {
        //     DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
        // } 

        // -- if diffusion coeficient is setup previously dont compute it
        // if (transportProperties.subDict(currSpName).found("D"))
        // {
        //     dimensionedScalar D(transportProperties.subDict(currSpName).lookup("D"));//R constant
        //     DEff = 0*DEff + D;
        // }

        // volScalarField DEff = DABLstEff_Fick()

        // DEff.correctBoundaryConditions();

        // --  calculation of total concentration
        // volScalarField cTot(thermo.rho()/Mg);
        // volScalarField cTot(rho/molMList[0]);

        // -- inform about calculated properties
        Info << "Min cTot: " << min(cTot).value() << ", max cTot: " << max(cTot).value() << endl;
        Info << "Min Mg: " << min(Mg).value() << ", max Mg: " << max(Mg).value() <<" Mgconst "<<molMList[0]<< endl;
        // Info << "Min DEff[" << specieNames[i] << "]: " << min(DEff).value() << ", max DEff["<<specieNames[i]<<"]: "<< max(DEff).value() << endl;
        // Info << "Min DFree[" << specieNames[i] << "]: " << min(DifMix).value() << ", max DFree["<<specieNames[i]<<"]: "<< max(DifMix).value() << endl;
        Info << "Min y[" << currSpName << "]: " << min(currSp).value() << ", max y[" << currSpName << "]: " << max(currSp).value() << endl;
        Info << "Min B[" << currSpName << "]: " << min(BinvList[i*dimMat+i]) <<  "Max B[" << currSpName << "]: " << max(BinvList[i*dimMat+i]) << endl;


        // -- calculate explicit contribution to the equation
        volVectorField explJi (-currSp * BinvList[i*dimMat+i] * fvc::grad(cTot));
        for (int j = 0; j < dimMat; j++)
        {
            if (j > i)
            {
                explJi -= BinvList[i*dimMat+j] * fvc::grad(cTot * MFrList[j]);
            }
            else if (j < i)
            {
                explJi -= BinvList[j*dimMat+i] * fvc::grad(cTot * MFrList[j]);
            }
        }
    
        // -- balance of the specie
        fvScalarMatrix concEq
        (
            // fvm::ddt(cTot,currSp)                                               
            // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(thermo.rho()*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div(phi/molMList[0], currSp)
            // fvm::div(phi/linearInterpolate(Mg), currSp)
            // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(rho*(U)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(rho*(U)/Mg) & mesh.Sf()), currSp)
            // fvm::div((linearInterpolate(rho*(U+UD)/Mg) & mesh.Sf()), currSp)
            // fvm::div(phi/molMList[molMList.size()-1], currSp)
            // fvm::div(linearInterpolate(cTot*U) & mesh.Sf(), currSp)
            // fvm::div((linearInterpolate(thermo.rho()*(U)) & mesh.Sf())/Mg, currSp)
            + fvc::div(linearInterpolate(explJi) & mesh.Sf())
            - fvm::laplacian(cTot*BinvList[i*dimMat+i],currSp)   
            ==
            // fvm::laplacian(DEff*cTot, currSp)                              
            // + fvc::laplacian(DEff*currSp, cTot)    
            // + fvModels.source(dummyConcEqn, currSp)                                      // source term defined in fvOptions, linearized -- part implicit, part explicit
            fvModels.source(dummyConcEqn, currSp)                                      // source term defined in fvModels, linearized -- part implicit, part explicit
        );
            
        // fvOptions.constrain(concEq);
        fvConstraints.constrain(concEq);

        concEq.relax();

        concEq.solve();

        fvConstraints.constrain(currSp);

        currSp.correctBoundaryConditions();
        
        forAll(currSp,cellI)
        { 
        //	 currSp[cellI] = min(max(currSp[cellI],0),1);                    //clip unreal molar fractions
            currSp[cellI] = max(currSp[cellI],0); 
        }

        sumWithoutRef += currSp;
        jIList[i] = - (DEff*(fvc::grad(cTot*currSp)));
        sumJDif += (DEff*(fvc::grad(cTot*currSp)));
    }
    else{
        currSp = dimensionedScalar("jedna",dimless,1) - sumWithoutRef;
        jIList[i] = - sumJDif;
    }
}

// -- update molar mass of the gas and velocity from diffusion
// UD = -sumJDif/MFrList[MFrList.size()-1]*Mg/rho;

Info << "Min UD: " << min(UD).value() << ", max UD: " << max(UD).value() << endl;

// -- update phi
// phi = linearInterpolate(rho*U) & mesh.Sf();
phiD = linearInterpolate(rho*UD) & mesh.Sf();

UD.correctBoundaryConditions();

