Info<< "Reading the concentrations \n" << endl;

// prepare list to store the concentrations in
PtrList<volScalarField> MFrList(specieNames.size());

forAll(specieNames,countI)
{
    word currentSp = specieNames[countI];
    Info<< "Reading the concentration of the: " << currentSp << endl;

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
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
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

// -- reading field UD -- velocity enforced by diffusion
Info<< "Reading field UD\n" << endl;
volVectorField UD
(
    IOobject
    (
        "UD",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


