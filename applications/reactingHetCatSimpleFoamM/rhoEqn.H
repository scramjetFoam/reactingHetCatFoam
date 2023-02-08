{
    fvScalarMatrix rhoEq
        (
            // -- specie mass convection
            fvm::div(linearInterpolate(U) & mesh.Sf(), rho)
            ==
            fvModels.source(rho)
        );
            
        // fvOptions.constrain(concEq);
        fvConstraints.constrain(rhoEq);

        rhoEq.relax();

        rhoEq.solve();

        fvConstraints.constrain(rho);

        rho.correctBoundaryConditions();

        p = rho * univR * thermo.T() / Mg;

        phi = linearInterpolate(U * rho) & mesh.Sf();
}
