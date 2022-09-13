{
    // -- temperature of the gas
    volScalarField& T = thermo.T();
    // MGradT = fvc::snGrad(thermo.T());

    // -- enthalpy balance for porous media
    fvScalarMatrix EEqnS
    (
      // -- diffusion of heat (kappaS = kappaEff/(1-epsEff), term = kappaS*(1-epsEff))
      - fvm::laplacian(kappaEff, TSolid) 
      ==
      // -- heat source from reaction 
      fvModels.source(dummyEEqn, TSolid)
      // -- heat source into gas
      + sumZones*fvc::Sp(kHEffA,T)
      - sumZones*fvm::Sp(kHEffA,TSolid)
    );

    EEqnS.relax();

    fvConstraints.constrain(EEqnS);

    EEqnS.solve();

    fvConstraints.constrain(TSolid);
}
