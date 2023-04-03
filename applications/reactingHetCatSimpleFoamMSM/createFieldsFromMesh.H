Info<< "Reading the concentrations \n" << endl;

// // prepare list to store mass fractions in
PtrList<volScalarField> MassFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI] + "Mass";
    Info<< "Reading the molar fraction of the: " << currentSp << endl;

    MassFrList.set
    (
        countI,
        new volScalarField
        (
            IOobject
            (
                currentSp,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

// -- creating field Mg -- molar mass of gas
Info<< "Reading field Mg\n" << endl;
volScalarField Mg
(
    IOobject
    (
        "Mg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption1
    ),
    mesh,
    dimensionedScalar("Mg",dimMass/dimMoles,0),
    "zeroGradient"
);

volScalarField sumWiMi
(
    IOobject
    (
        "sumWiMi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    dimensionedScalar("sumWiMi",dimMoles/dimMass,0),
    "zeroGradient"
);

forAll(specieNames,lstI)
{
    sumWiMi += MassFrList[lstI] / molMList[lstI];
    sumWiMi.correctBoundaryConditions();
}
Mg = 1 / sumWiMi;
Mg.correctBoundaryConditions();

// -- prepare list to store molar fractions in (read from BC)
PtrList<volScalarField> MFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSpMol = specieNames[countI];
    Info<< "Computing the molar fraction of the: " << currentSpMol << endl;
    // Info<< "MFrList[countI]: " << min(MFrList[countI]).value() << endl;

    volScalarField currentSpMolar = MassFrList[countI] * Mg / molMList[countI];

    MFrList.set
    (
        countI,
        new volScalarField
        (
            IOobject
            (
                currentSpMol,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            currentSpMolar
        )
    );
}

// // loop over all boundaries
// forAll(MFrList[0].boundaryField(), bfI)
// {
//     if (MFrList[0].boundaryField()[bfI].type() == "fixedValue")
//     {
//         Mg.boundaryFieldRef().set(bfI , fvPatchField<scalar>::New("fixedValue", mesh.boundary()[bfI], Mg));
//         const polyPatch& patch = mesh.boundaryMesh()[bfI];
//         forAll(patch , faceI)
//         {
//             Mg.boundaryFieldRef()[bfI][faceI] = (MgOld.boundaryFieldRef()[bfI][faceI]);
//         }
//     }
//     Info << "I set " << Mg.boundaryField()[bfI] << endl;
// }

// prepare list to store mass fractions in
// PtrList<volScalarField> MassFrList(specieNames.size());

// forAll(specieNames,countI)
// {
//     word currentSp = specieNames[countI] + "Mass";
//     Info<< "Preparing the mass fraction of the: " << currentSp << endl;

//     volScalarField currentSpMass = MFrList[countI] / Mg * molMList[countI];

//     // // loop over all boundaries
//     forAll(currentSpMass.boundaryField(), bfI)
//     {
//         if (MFrList[countI].boundaryField()[bfI].type() == "fixedValue")
//         {
//             currentSpMass.boundaryFieldRef().set(bfI , fvPatchField<scalar>::New("fixedValue", mesh.boundary()[bfI], currentSpMass));
//             // Mg.boundaryFieldRef().set(bfI , fvPatchField<scalar>::New("fixedValue", mesh.boundary()[bfI], Mg));

//             const polyPatch& patch = mesh.boundaryMesh()[bfI];

//             forAll(patch , faceI)
//             {
//             //     // Mg.boundaryFieldRef()[bfI][faceI] = (Mg.boundaryFieldRef()[bfI][faceI]);
//                 currentSpMass.boundaryFieldRef()[bfI][faceI] = (MFrList[countI].boundaryFieldRef()[bfI][faceI] / Mg.boundaryFieldRef()[bfI][faceI] * molMList[countI]).value();
//             }
//         }
//         else if (MFrList[countI].boundaryField()[bfI].type() == "fixedValue")
//         {
//             currentSpMass.boundaryFieldRef().set(bfI , fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[bfI], currentSpMass));
//         }
//         Info << "I set " << currentSpMass.boundaryField()[bfI] << endl;
//     }
    
//     MassFrList.set
//     (
//         countI,
//         new volScalarField
//         (
//             IOobject
//             (
//                 currentSp,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             currentSpMass
//             // mesh
//         )
//     );
//     MassFrList[countI].write();
// }

// -- reading field U
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// -- T field for porousMedia from mesh
// Info<< "Reading field TSolid\n" << endl;
// volScalarField TSolid
// (
//     IOobject
//     (
//         "TSolid",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );





