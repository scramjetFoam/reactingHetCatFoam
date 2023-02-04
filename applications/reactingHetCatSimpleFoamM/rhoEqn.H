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

        fvConstraints.constrain(currSpMass);

        currSpMass.correctBoundaryConditions();

        sumWithoutRef += currSpMass;
}
