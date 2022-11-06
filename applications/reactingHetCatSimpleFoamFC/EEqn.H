{
    // -- load temperature from thermo package
    volScalarField& T = thermo.T();

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
    // kappaEff += freeChannel*kappa;         
    forAll(zoneLst,countI)
    {
        kappaEff += zoneLst[countI]*(kappaZ[countI]);
    }
    volScalarField kappaPM = kappaEff + kappa*epsEff;


    if (debug == 1){
        // -- minimal and maximal values of calculated fields
        Info << "Min (kappaEff): " << min(kappaEff).value() << ", max (kappaEff): " << max(kappaEff).value() << "." << endl;
        Info << "Min (kappaPokus): " << min(kappaPM).value() << ", max (kappaPokus): " << max(kappaPM).value() << "." << endl;
        Info << "Min (thermo.Cp): " << min(thermo.Cp()).value() << ", max (Cp): " << max(thermo.Cp()).value() << "." << endl;
        Info << "Min (thermo.rho): " << min(thermo.rho()).value() << ", max (thermo.rho): " << max(thermo.rho()).value() << "." << endl;
    }

    // enthalpy balance of the porous media temperature
    fvScalarMatrix EEqn
    (
        // -- convection
        fvm::div(phi*linearInterpolate(thermo.Cp()), T)
        ==
        // -- diffusion
        fvm::laplacian(kappaPM, T)
        // -- heat source
        + fvModels.source(dummyEEqn, T)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(T);

    thermo.he() = T * thermo.Cp();

    thermo.correct();

    Info << "Min (thermo.T): " << min(thermo.T()).value() << ", max (thermo.T): " << max(thermo.T()).value() << "." << endl;
}