{
  // temperature of the gas
  volScalarField& T = thermo.T();
  // volScalarField& he = thermo.he();
    
  // -- calculate the thermal diffusivity of gas
  volScalarField kappa (kappaA+kappaB*T+kappaC*Foam::pow(T,2)+kappaD*Foam::pow(T,3));

  
  // -- calculate heat capacity of the gas (reference at the moment)
  // NOTE: old, now using OpenFOAM thermo
  // volScalarField Cp
  // (
  //   (ALst[ALst.size()-1] + BLst[BLst.size()-1] * T + CLst[CLst.size()-1] * Foam::pow(T,2) + DLst[DLst.size()-1] * Foam::pow(T,3)) / molMList[molMList.size()-1]
  // );
      
  // NOTE: discuss how to treat CpEff and rhoEff (probably will be important for dynamics)
  
  // dissipation rate  -- not used at the moment (neglectable)
  // volScalarField K (0.5*magSqr(U));
    
  //efective rho * Cp
  // volScalarField rhoCpEff(thermo.rho()*thermo.Cp());

  // -- minimal and maximal values of calculated fields
  Info << "Min (kappaEff): " << min(kappaEff).value() << ", max (kappaEff): " << max(kappaEff).value() << "." << endl;
  Info << "Min (thermo.Cp): " << min(thermo.Cp()).value() << ", max (Cp): " << max(thermo.Cp()).value() << "." << endl;
  Info << "Min (thermo.rho): " << min(thermo.rho()).value() << ", max (thermo.rho): " << max(thermo.rho()).value() << "." << endl;
  Info << "Min rhoEff: " << min(rhoEff).value() << ", max rhoEff: " << max(rhoEff).value() << "." << endl;
  Info << "Min CpEff: " << min(CpEff).value() << ", max CpEff: " << max(CpEff).value() << "." << endl;
  
  //enthalpy balance for porous media
    fvScalarMatrix EEqnPM
  (
      CpEff*fvm::ddt(rhoEff,TporMedia)
      ==
      //-- ensure zero gradient on the porous media
      sumZones*fvm::laplacian(kappaEff, TporMedia) +
      - freeChannel*fvm::laplacian(kappaEff, TporMedia) +
      //-- heat source 
       sumZones*fvModels.source(dummyEEqn, TporMedia)

      // - sumZones*fvm::Sp(kHEff,T)
      + sumZones*fvc::Sp(kHEff,T)
      - sumZones*fvm::Sp(kHEff,TporMedia)
      // - kHEff*(TporMedia-T)
    //   fvm::ddt(rho, he) 
    //   + fvm::div(phi, he)
    //   // -dpdt
    // ==
    //   fvm::laplacian(kappaEff/thermo.Cp(), he)
    //   + fvModels.source(rho, he)
  );
  EEqnPM.relax();

  fvConstraints.constrain(EEqnPM);

  EEqnPM.solve();

  // TporMedia = sumZones*TporMedia;

  fvConstraints.constrain(TporMedia);


  Info << "Solving temperature of gas" << endl;

  // enthalpy balance for gas
  fvScalarMatrix EEqn
  (
      // CpEff*fvm::ddt(rhoEff, T) +
      // thermo.Cp()*fvm::ddt(thermo.rho(), T) +
      // CpEff*rhoEff*fvm::ddt(T) + fvm::Sp(fvc::ddt(CpEff,rhoEff),T)+
      thermo.Cp()*rho*fvm::ddt(T) + fvm::Sp(fvc::ddt(thermo.Cp(),rho),T)+
      // CpEff*fvm::ddt(rhoEff, T) +
      fvc::ddt(rho, K) + fvc::div(phi, K) +
      fvm::div(linearInterpolate(thermo.Cp())*phi, T)
      -dpdt
      ==
      fvm::laplacian(kappa, T)
      // + fvModels.source(dummyEEqn, T)
      - sumZones*fvm::Sp(kHEff,T)
      + sumZones*fvc::Sp(kHEff,TporMedia)
      // + sumZones*kHEff*TporMedia

    //   fvm::ddt(rho, he) 
    //   + fvm::div(phi, he)
    //   // -dpdt
    // ==
    //   fvm::laplacian(kappaEff/thermo.Cp(), he)
    //   + fvModels.source(rho, he)
  );
  
  EEqn.relax();

  fvConstraints.constrain(EEqn);

  EEqn.solve();

  // fvOptions.correct(thermo.T());
  fvConstraints.constrain(T);
  // fvConstraints.constrain(he);

  // thermo.he() = T * thermo.Cp();
  thermo.he() = T * thermo.Cp();

  thermo.correct();
  Info << "Min (TporMedia): " << min(TporMedia).value() << ", max (TporMedia): " << max(TporMedia).value() << "." << endl;

  Info << "Min (thermo.T): " << min(thermo.T()).value() << ", max (thermo.T): " << max(thermo.T()).value() << "." << endl;
  
  //~ radiation->correct();
}
