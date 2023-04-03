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
    dimensionedScalar("Mg",dimMass/dimMoles,0)
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
    dimensionedScalar("sumWiMi",dimMoles/dimMass,0)
);

forAll(specieNames,lstI)
{
    sumWiMi += MassFrList[lstI] / molMList[lstI];
}
Mg = 1 / sumWiMi;

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

// -- creating field UD -- velocity enforced by diffusion
volVectorField UD
(
    IOobject
    (
        "UD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
    dimensionedVector ("UD",dimLength/dimTime, vector(0.0,0.0,0.0))
    // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
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

// -- T field from mesh
// Info<< "Reading field T\n" << endl;
// volScalarField T
// (
//     IOobject
//     (
//         "T",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );


