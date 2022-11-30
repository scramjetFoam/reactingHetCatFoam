{

    // -- minimal and maximal values of calculated fields
    if (debug ==1){
    Info << "Min (kappaEff): " << min(kappaEff).value() << ", max (kappaEff): " << max(kappaEff).value() << "." << endl;
    Info << "Min (epsEff): " << min(epsEff).value() << ", max (epsEff): " << max(epsEff).value() << "." << endl;
    Info << "Min (thermo.Cp): " << min(thermo.Cp()).value() << ", max (Cp): " << max(thermo.Cp()).value() << "." << endl;
    Info << "Min (thermo.rho): " << min(thermo.rho()).value() << ", max (thermo.rho): " << max(thermo.rho()).value() << "." << endl;
    Info << "Min rhoEff: " << min(rhoS).value() << ", max rhoEff: " << max(rhoS).value() << "." << endl;
    Info << "Min CpEff: " << min(CpS).value() << ", max CpEff: " << max(CpS).value() << "." << endl;
    }
    // -- temperature of the gas
    volScalarField& T = thermo.T();
    // MGradT = fvc::snGrad(thermo.T());

    // -- enthalpy balance for porous media
    fvScalarMatrix EEqnS
    (
      // -- accumulation
      (1-epsEff)*CpS*rhoS*fvm::ddt(TSolid)
      // -- diffusion
      - fvm::laplacian((1-epsEff)*kappaEff, TSolid) 
      ==
      // -- heat source from reaction 
      fvModels.source(dummyEEqn, TSolid)+
      // -- heat source into gas
      sumZones*fvc::Sp(kHEffA,T)
      - sumZones*fvm::Sp(kHEffA,TSolid)
    );

    EEqnS.relax();

    fvConstraints.constrain(EEqnS);

    EEqnS.solve();

    fvConstraints.constrain(TSolid);
}
