{
  // -- temperature of the gas
  volScalarField& T = thermo.T();
  // MGradT = fvc::snGrad(thermo.T());

  // -- calculate the thermal diffusivity of gas
  volScalarField kappa (kappaA+kappaB*T+kappaC*Foam::pow(T,2)+kappaD*Foam::pow(T,3));

  // -- enthalpy balance for gas
  fvScalarMatrix EEqn
  (
    // -- accumulation
    //epsEff*thermo.Cp()*rho*fvm::ddt(T) + epsEff*fvm::Sp(fvc::ddt(thermo.Cp(),rho),T)+
    thermo.Cp()*(fvm::ddt(rho,T)) 
    // -- convection
    +fvm::div(linearInterpolate(thermo.Cp())*phi, T)
    //-dpdt
    //+fvm::div((linearInterpolate(rho*thermo.Cp()*(U)) & mesh.Sf()), T)	
    ==
    // -- diffusion
    fvm::laplacian(kappa*epsEff, T)


//    // -- heat source from reaction 
//    +fvModels.source(dummyEEqn, T)
    // -- heat source from solid media
     - sumZones*fvm::Sp(kHEffA,T)
     + sumZones*fvc::Sp(kHEffA,TSolid)
    // NOTE TH: 
    // implement surface heat source from porous media at solid interface
  );

  EEqn.relax();

  fvConstraints.constrain(EEqn);

  EEqn.solve();

  fvConstraints.constrain(T);

  thermo.he() = T * thermo.Cp();

  thermo.correct();

  Info << "Min (TporMedia): " << min(TSolid).value() << ", max (TporMedia): " << max(TSolid).value() << "." << endl;

  Info << "Min (thermo.T): " << min(thermo.T()).value() << ", max (thermo.T): " << max(thermo.T()).value() << "." << endl;
}
