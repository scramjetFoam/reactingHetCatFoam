{
    // -- load temperature from thermo package
    T = thermo.T();
    
    // -- calculate the effective thermal diffusivity
    volScalarField kappa (kappaA+kappaB*T+kappaC*Foam::pow(T,2)+kappaD*Foam::pow(T,3));
    volScalarField kappaEff
    (
      IOobject
      (
          "kappaEff",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("kappaEff",kappa.dimensions(),0.0)
    );

    // -- efective kappa values, gas kappa in free channel, media kappa in zones
    kappaEff += freeChannel*kappa;         
    forAll(zoneLst,countI)
    {
      kappaEff += zoneLst[countI]*(kappaZ[countI]);
    } 
    
    // -- calculate heat capacity of the gas (reference at the moment)
    // NOTE: old, now using OpenFOAM thermo
    // volScalarField Cp
    // (
    //   (ALst[ALst.size()-1] + BLst[BLst.size()-1] * T + CLst[CLst.size()-1] * Foam::pow(T,2) + DLst[DLst.size()-1] * Foam::pow(T,3)) / molMList[molMList.size()-1]
    // );
      	
    // NOTE: discuss how to treat CpEff and rhoEff (probably will be important for dynamics)
    // -- calculate efective heat capacity -- not used at the moment
    // volScalarField CpEff
    // (
    //   (1.0-max(isOnWallCoat,isInWallCoat))*Cp  
    //   + isInWallCoat*(porEpsS.value()*Cp+(1.0-porEpsS.value())*CpS)
    //   + isOnWallCoat*(porEpsC.value()*Cp+(1.0-porEpsC.value())*CpC)
    // );
    
    //efective rho
    // volScalarField rhoEff
    // (
    //   (1.0-max(isOnWallCoat,isInWallCoat))*rho
    //   // (1.0-max(isOnWallCoat,isInWallCoat))*thermo.rho()
    //   + isOnWallCoat*(rhoC)
    //   + isInWallCoat*(rhoS)
    // );    
    
    
    // dissipation rate  -- not used at the moment (neglectable)
    // volScalarField K (0.5*magSqr(U));
    
    // -- heat flux surfaceField
    surfaceScalarField phiRhoCpMy
    (
        IOobject
        (
            "phiRhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        // linearInterpolate(U*rho*Cp) & mesh.Sf()
        linearInterpolate(U*thermo.rho()*thermo.Cp()) & mesh.Sf()
    );

    // -- minimal and maximal values of calculated fields
    Info << "Min (kappaEff): " << min(kappaEff).value() << ", max (kappaEff): " << max(kappaEff).value() << "." << endl;
    Info << "Min (thermo.Cp): " << min(thermo.Cp()).value() << ", max (Cp): " << max(thermo.Cp()).value() << "." << endl;
    Info << "Min (thermo.rho): " << min(thermo.rho()).value() << ", max (thermo.rho): " << max(thermo.rho()).value() << "." << endl;
    
    // enthalpy balance
    fvScalarMatrix EEqn
    (
        fvm::div(phiRhoCpMy, T)
        ==
        fvm::laplacian(kappaEff, T)
        // + fvOptions(dummyEEqn,T)
        // + fvModels.source(dummyEEqn, T)
        + fvModels.source(dummyEEqn, T)
    );
    
    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    // fvOptions.correct(thermo.T());
    fvConstraints.constrain(T);

    thermo.he() = T * thermo.Cp();

    thermo.correct();

    Info << "Min (thermo.T): " << min(thermo.T()).value() << ", max (thermo.T): " << max(thermo.T()).value() << "." << endl;
    
    //~ radiation->correct();
}
