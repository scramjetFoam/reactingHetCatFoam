Info<< "Reading the concentrations \n" << endl;

// // prepare list to store mass fractions in
// PtrList<volScalarField> MassFrList(specieNames.size());

// forAll(specieNames,countI)
// {
//     word currentSp = specieNames[countI] << "Mass";
//     Info<< "Reading the molar fraction of the: " << currentSp << endl;

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
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//         )
//     );
// }

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
    dimensionedScalar("Mg",dimMass/dimMoles,1),
    "zeroGradient"
);

// prepare list to store mass fractions in
PtrList<volScalarField> MassFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI] + "Mass";
    Info<< "Preparing the mass fraction of the: " << currentSp << endl;

    MassFrList.set
    (
        countI,
        new volScalarField
        // volScalarField
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
            // MFrList[countI] * molMList[countI] / Mg,
            // "fixedValue"
        )
    );
}

volScalarField sumWiMi = 0/Mg;
forAll(specieNames,lstI)
{
    sumWiMi += MassFrList[lstI] / molMList[lstI];
}
Mg = 1/sumWiMi;

// prepare list to store molar fractions in
PtrList<volScalarField> MFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    Info<< "Reading the molar fraction of the: " << currentSp << endl;

    MFrList.set
    (
        countI,
        new volScalarField
        (
            IOobject
            (
                currentSp,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            // mesh
            MassFrList[countI] * Mg / molMList [countI]
        )
    );
}

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

// -- reading field UD -- velocity enforced by diffusion
// Info<< "Reading field UD\n" << endl;
// volVectorField UD
// (
//     IOobject
//     (
//         "UD",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );




