fvScalarMatrix elEq
(
    // -- specie mass convection
    fvm::laplacian(sigmaEff, elPhi)   

    // ==

    // -- source
    // fvModels.source(elPhi)
);
            
// fvOptions.constrain(concEq);
fvConstraints.constrain(elEq);

elEq.relax();

elEq.solve();

fvConstraints.constrain(elPhi);

elPhi.correctBoundaryConditions();


        
