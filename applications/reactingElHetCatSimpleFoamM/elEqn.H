sigmaEffEl = sigmaEffEl*0;
sigmaEffIon = sigmaEffIon*0;

forAll(zoneLst,countI)
{
    sigmaEffEl += zoneLst[countI]*sigmaScInfEl[countI] * Foam::exp(-sigmaScEAEl[countI]/thermo.T());
    sigmaEffIon += zoneLst[countI]*sigmaScInfIon[countI] * Foam::exp(-sigmaScEAIon[countI]/thermo.T());
}


// sigmaEffEl += 1*freeChannel*dimensionedScalar("je",sigmaScInfEl[0].dimensions(),1) + 1*(1-freeChannel)*dimensionedScalar("je",sigmaScInfEl[0].dimensions(),1);
// sigmaEffIon += 1*freeChannel*dimensionedScalar("je",sigmaScInfEl[0].dimensions(),1) + 1*(1-freeChannel)*dimensionedScalar("je",sigmaScInfEl[0].dimensions(),1);

sigmaEffEl.correctBoundaryConditions();
sigmaEffIon.correctBoundaryConditions();


fvScalarMatrix elEq
(
    // -- specie mass convection
    - fvm::laplacian(sigmaEffEl, elPhi)   

    // ==

    // -- source
    // fvModels.source(dummyElEqn,elPhi)
);

            
// fvOptions.constrain(concEq);


elEq.relax();

fvConstraints.constrain(elEq);

elEq.solve();

fvConstraints.constrain(elPhi);


// forAll(elPhi,cellI)
// { 
//     elPhi[cellI] = max((elPhi[cellI]),SMALL); 
// }

elPhi.correctBoundaryConditions();


fvScalarMatrix elIonEq
(
    // -- specie mass convection
    - fvm::laplacian(sigmaEffIon, elIon)   

    // ==

    // -- source
    // fvModels.source(dummyElEqn,elIon)
);
            
// fvOptions.constrain(concEq);
elIonEq.relax();

fvConstraints.constrain(elIon);

elIonEq.solve();

fvConstraints.constrain(elIon);

elIon.correctBoundaryConditions();     
   

// forAll(elIon,cellI)
// { 
//     elIon[cellI] = max((elIon[cellI]),SMALL); 
// }