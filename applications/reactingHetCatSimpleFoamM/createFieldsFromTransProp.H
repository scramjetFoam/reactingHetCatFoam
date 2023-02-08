// -- creating field with sum of the diffusive fluxes of all species but reference
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

// surfaceScalarField phiD = linearInterpolate(UD * )

// surfaceVectorField UDRef
// (
//     IOobject
//     (
//         "UDRef",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         writeOption2
//     ),
//     mesh,
//     // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
//     dimensionedVector ("UDRef", dimLength/dimTime, vector(0.0,0.0,0.0))
//     // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
// );

volVectorField mD
(
    IOobject
    (
        "mD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
    dimensionedVector ("mD",dimMass/dimVolume*dimLength/dimTime, vector(0.0,0.0,0.0))
    // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
);

surfaceScalarField UDF
(
    IOobject
    (
        "UDF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
    dimensionedScalar ("UDF",dimLength/dimTime, 0)
    // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
);

surfaceVectorField UDFS
(
    IOobject
    (
        "UDFS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
    dimensionedVector ("UDFS",dimLength/dimTime, vector(1.0,0.0,0.0))
    // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
);

surfaceVectorField UDFSOld
(
    IOobject
    (
        "UDFSOld",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
    dimensionedVector ("UDFSOld",dimless, vector(1.0,0.0,0.0))
    // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
);

volVectorField mDBC
(
    IOobject
    (
        "mDBC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
    dimensionedVector ("mDBC",dimMass/dimVolume*dimLength/dimTime, vector(0.0,0.0,0.0)),
    "fixedValue"
    // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
);

// volScalarField sumMDi
// (
//     IOobject
//     (
//         "sumMDi",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         writeOption2
//     ),
//     mesh,
//     // dimensionedVector ("sumJMDiff",dimLength/dimTime, vector(0.0,0.0,0.0))
//     dimensionedScalar ("sumMDi",dimMass/dimVolume/dimTime, 0)
//     // dimensionedVector ("sumJMDiff",dimLength/dimTime*dimMass/dimMoles, vector(0.0,0.0,0.0))
// );

// // -- creating field with sum of the diffusive molar fluxes of all species but reference
// volVectorField sumJDiff
// (
//     IOobject
//     (
//         "sumJDiff",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         writeOption2
//     ),
//     mesh,
//     dimensionedVector ("sumJDiff",dimLength/dimTime*dimMoles/dimVolume, vector(0.0,0.0,0.0))
// );

// -- prepare sum of the molar fractions without reference
volScalarField sumWithoutRef
(
    IOobject
    (
        "sumWithoutRef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    dimensionedScalar("sumWithoutRef",dimless,1.0)
);

// -- zones volScalarFields -- from transport properties or from file
PtrList<volScalarField> zoneLst(zoneNames.size());
forAll(zoneNames,countI)
{
    word currentZone = zoneNames[countI];
    Info<< "Reading the zone " << currentZone << " and saving it as volScalarField." << endl;

    bool loadFromFile = false;

    if (transportProperties.subDict(currentZone).found("loadFromFile"))
    {
        transportProperties.subDict(currentZone).lookup("loadFromFile") >> loadFromFile;
    }

    // -- load from file
    if (loadFromFile)
    {
        zoneLst.set
        (
            countI,
            new volScalarField
            (
                IOobject
                (
                    currentZone,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        Info << "Zone " << currentZone << " has been loaded from file." << endl;
    }
    // -- or create with zeros and then set according to cellZone
    else
    {   
        // -- create with zeros
        zoneLst.set
        (
            countI,
            new volScalarField
            (
                IOobject
                (
                    currentZone,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    writeOption2
                ),
                mesh,
                dimensionedScalar (currentZone,dimless,0.0),
                "zeroGradient"
            )
        );

        // -- set 1 where the cellZone is
        word currentZone = zoneNames[countI];
        labelList cellZoneIDs_ = mesh.cellZones().findIndices(currentZone);
        const labelList& cells = mesh.cellZones()[cellZoneIDs_[0]];
        forAll(cells, celli)
        {
            zoneLst[countI][cells[celli]] = 1;
        }
        zoneLst[countI].correctBoundaryConditions();
        Info << "Zone " << currentZone << " has been created from cellZone." << endl;
    }
}



