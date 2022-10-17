// -- function to compure Steffan Maxwell diffusion coefficients 

PtrList<dimensionedScalar> DABLstEff_Fick(
                                            PtrList<dimensionedScalar>& yi,
                                            volScalarField& eps,
                                            fluidThermo& thermo,
                                            PtrList<dimensionedScalar>& molMList,
                                            PtrList<dimensionedScalar>& volDiffLst,
                                            wordList& specieNames,
                                            int& lstI,
                                            volScalarField& DKn,
                                        )
{
    PtrList<dimensionedScalar> DABLst(specieNames.size());
    forAll(DABLst,lstII)
    {
        // -- diffusion coefficients according to Fuller without temperature and pressure dependence
        DABLst.set
        (
            lstII,
            new dimensionedScalar(
                "DAB",
                Foam::pow(dimLength,2)/dimTime/Foam::pow(dimTemperature,1.75)*dimPressure,
                ((0.014)/(Foam::sqrt(2.0/(1e-3/molMList[lstII] + 1e-3/molMList[lstI]))*Foam::pow(Foam::pow(volDiffLst[lstI],1.0/3.0)+Foam::pow(volDiffLst[lstII],1.0/3.0),2))).value()
            )
        );
    }   
    // -- Diffusion coeficient of the mixture -- default reference gas diffusion (last specie) + add temperature and pressure dependence (see Fuller)
    volScalarField DifMix
    (
        DABLst[DABLst.size()-1]*Foam::pow(thermo.T(),1.75)/thermo.p()
    );

    DKn += freeChannel*DifMix;         // free-channel (not used, but cant be 0)
    forAll(zoneLst,countI)
    {
        DKn += zoneLst[countI]*(1.0/3.0*sqrt((thermo.T()*8*univR)/(3.141592653*molMList[lstI])))*dp[countI];
    }

    DEff = freeChannel*DifMix;         // free-channel
    forAll(zoneLst,countI)
    {
        DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
    } 

    return DEff;
}