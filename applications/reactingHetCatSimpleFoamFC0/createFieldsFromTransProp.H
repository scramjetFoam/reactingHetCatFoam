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
    dimensionedScalar("Mg",dimMass/dimMoles,0.0),
    "zeroGradient"
);

forAll(specieNames,countI)
{
    Mg += molMList[countI]*MFrList[countI];
}

// -- creating field with sum of the diffusive fluxes of all species but reference
volVectorField sumJDif
(
    IOobject
    (
        "sumJDif",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        writeOption2
    ),
    mesh,
    dimensionedVector ("sumJDif",dimLength/dimTime*dimMoles/dimVolume, vector(0.0,0.0,0.0))
);

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
    dimensionedScalar("sumWithoutRef",dimless,0.0)
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



