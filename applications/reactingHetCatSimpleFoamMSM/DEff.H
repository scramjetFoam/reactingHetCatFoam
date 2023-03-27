// -- function to compure Fick effective diffusion coefficients 

void DEffFnc
(
    PtrList<dimensionedScalar>& eps,
    PtrList<dimensionedScalar>& tort,
    PtrList<volScalarField>& zoneLst,
    volScalarField& freeChannel,
    volScalarField& DEff,
    volScalarField& DKn,
    PtrList<dimensionedScalar>& DijList,
    fluidThermo& thermo,
    PtrList<dimensionedScalar>& molMList,
    PtrList<dimensionedScalar>& dp,
    PtrList<dimensionedScalar>& volDiffLst,
    wordList& specieNames,
    IOdictionary& transportProperties,
    dimensionedScalar& univR,
    int i,
    int j,
    int n
)
{
    volScalarField DifMix
    (
        DijList[i*(n)+j]*Foam::pow(thermo.T(),1.75)/thermo.p()
    );
    // -- calculate the Knudsen difussivity
    DKn = 0*DifMix;                    
    DKn += freeChannel*DifMix;         // free-channel (not used, but cant be 0)
    forAll(zoneLst,countI)
    {
        DKn += zoneLst[countI]*(1.0/3.0*sqrt((thermo.T()*8*univR)/(3.141592653*molMList[i])))*dp[countI];
    }

    word DTu = "D" + (name(i)) + (name(j));
    // Info << DTu <<endl;

    // -- if diffusion coeficient is set in transportProperties dont compute it
    if (transportProperties.subDict(specieNames[0]).found(DTu))
    {
        // Info << "Taking diff here." << endl;
        dimensionedScalar D(transportProperties.subDict(specieNames[0]).lookup(DTu));//R constant
        DifMix = 0*DifMix + D;
    }
    // else{
    //     Info << "Diff not taken"<<endl;
    // }

    // -- effective diffusion coeficient
    DEff = freeChannel*DifMix;         // free-channel
    // volScalarField DEff = epsEff/tortEff*DifMix;         // free-channel
    forAll(zoneLst,countI)
    {
        DEff += zoneLst[countI]*(eps[countI]/(tort[countI]*(1.0/DifMix+1.0/DKn)));  // porous media
    } 
    // return DEff;
}