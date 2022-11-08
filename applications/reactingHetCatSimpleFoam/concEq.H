// -- solve the concentration equation for all chemical species
forAll(specieNames,lstI)
{
    // -- get current name of the specie that is solved for
    word currSpName = specieNames[lstI];                            	

    Info<<"\nWorking on the specie: "<<currSpName<<endl;

    volScalarField& currSp	= MFrList[lstI];                	//current concentration field

    // -- DIFUSION COEFICIENT CALCULATION 
    // -- difussivity of the specie in the reference gas see 22_Hlavaty_MSc

    // -- prepare free difusion coeficients of the current specie against all other species 
    // -- NOTE: at the moment we use only diffusion of the specie in the reference gas (nitrogen) -- last member of the list
    PtrList<dimensionedScalar> DABLst(specieNames.size()+1);
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

    // -- if diffusion coeficient is setup previously dont compute it
    if (transportProperties.subDict(currSpName).found("D"))
    {
        dimensionedScalar D(transportProperties.subDict(currSpName).lookup("D"));//R constant
        DifMix = 0*DifMix + D;
        Info << "Difussivity loaded from transport properties" <<endl;
    }

    // -- prepare Knudsen diffusivity
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

    DEff.correctBoundaryConditions();

    // --  calculation of total concentration
    volScalarField cTot(thermo.rho()/molMList[molMList.size()-1]);
    // volScalarField cTot(rho/molMList[molMList.size()-1]);

    // -- inform about calculated properties
    Info << "Min cTot: " << min(cTot).value() << ", max cTot: " << max(cTot).value() << endl;
    Info << "Min DEff[" << specieNames[lstI] << "]: " << min(DEff).value() << ", max DEff["<<specieNames[lstI]<<"]: "<< max(DEff).value() << endl;
    Info << "Min DFree[" << specieNames[lstI] << "]: " << min(DifMix).value() << ", max DFree["<<specieNames[lstI]<<"]: "<< max(DifMix).value() << endl;
    Info << "Min y[" << specieNames[lstI] << "]: " << min(currSp).value() << ", max y[" << specieNames[lstI] << "]: " << max(currSp).value() << endl;


    // -- calculate convective molar flux at the face 
    // surfaceScalarField phiUCTot
    // (
    //     IOobject
    //     (
    //         "phiU",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     linearInterpolate(U*cTot) & mesh.Sf()
    // );
    

    // -- balance of the specie
    fvScalarMatrix concEq
    (
        // fvm::ddt(cTot,currSp)                                               
        fvm::div(phi/molMList[molMList.size()-1], currSp)
        ==
        fvm::laplacian(DEff*cTot, currSp)                              
        + fvc::laplacian(DEff*currSp, cTot)       
        // + fvModels.source(dummyConcEqn, currSp)                                      // source term defined in fvOptions, linearized -- part implicit, part explicit
        + fvModels.source(dummyConcEqn, currSp)                                      // source term defined in fvModels, linearized -- part implicit, part explicit
    );

    gradCCO = mag(fvc::snGrad(DEff*MFrList[0]*thermo.p()/thermo.T()/univR));
    gradYCO = mag(fvc::snGrad(DEff*MFrList[0]));
    cCOS = mag(fvc::interpolate(DEff*MFrList[0]*thermo.p()/thermo.T()/univR));
    yCOS = mag(fvc::interpolate(DEff*MFrList[0]));
	    
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
}