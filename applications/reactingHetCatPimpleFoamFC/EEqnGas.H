{
  // -- temperature of the gas
  volScalarField& T = thermo.T();
  // MGradT = fvc::snGrad(thermo.T());

  // -- calculate the thermal diffusivity of gas
  volScalarField kappa (kappaA+kappaB*T+kappaC*Foam::pow(T,2)+kappaD*Foam::pow(T,3));

  // -- minimal and maximal values of calculated fields
  Info << "Min (kappaEff): " << min(kappaEff).value() << ", max (kappaEff): " << max(kappaEff).value() << "." << endl;
  Info << "Min (epsEff): " << min(epsEff).value() << ", max (kappaEff): " << max(epsEff).value() << "." << endl;
  Info << "Min (thermo.Cp): " << min(thermo.Cp()).value() << ", max (Cp): " << max(thermo.Cp()).value() << "." << endl;
  Info << "Min (thermo.rho): " << min(thermo.rho()).value() << ", max (thermo.rho): " << max(thermo.rho()).value() << "." << endl;
  Info << "Min rhoEff: " << min(rhoS).value() << ", max rhoEff: " << max(rhoS).value() << "." << endl;
  Info << "Min CpEff: " << min(CpS).value() << ", max CpEff: " << max(CpS).value() << "." << endl;

  // -- enthalpy balance for gas
  fvScalarMatrix EEqn
  (
    // -- accumulation
    epsEff*thermo.Cp()*rho*fvm::ddt(T) + epsEff*fvm::Sp(fvc::ddt(thermo.Cp(),rho),T)+
    // -- convection
    fvm::div(linearInterpolate(thermo.Cp())*phi, T)
    ==
    // -- diffusion
    fvm::laplacian(kappa*epsEff, T)


    // -- heat source from reaction 
    +fvModels.source(dummyEEqn, T)
    // -- heat source from solid media
    // - sumZones*fvm::Sp(kHEffA,T)
    // + sumZones*fvc::Sp(kHEffA,TSolid)
    // NOTE TH: 
    // implement surface heat source from porous media at solid interface
  );

  EEqn.relax();

  fvConstraints.constrain(EEqn);

  EEqn.solve();

  fvConstraints.constrain(T);

  thermo.he() = T * thermo.Cp();

  thermo.correct();

  // Info << "Min (TporMedia): " << min(TSolid).value() << ", max (TporMedia): " << max(TSolid).value() << "." << endl;

  Info << "Min (thermo.T): " << min(thermo.T()).value() << ", max (thermo.T): " << max(thermo.T()).value() << "." << endl;
}