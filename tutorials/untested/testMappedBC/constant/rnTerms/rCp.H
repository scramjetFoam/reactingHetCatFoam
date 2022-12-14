// function to compure standard reaction enthalphy at given temperature
volScalarField sHatT(dimensionedScalar sHr, 
			  PtrList<dimensionedScalar>& A, 
			  PtrList<dimensionedScalar>& B, 
			  PtrList<dimensionedScalar>& C, 
			  PtrList<dimensionedScalar>& D, 
			  List<scalar>& nuLst,
			  wordList specieNames,
			  const fvMesh& mesh_)
{
	const Time& time = mesh_.time();
	const volScalarField& T = mesh_.lookupObject<volScalarField>("T"); 
	dimensionedScalar Tref("Tref",dimTemperature,298.0);
	volScalarField suma
	(
		IOobject
		(
			"suma",
			time.timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("zero",dimEnergy/dimMoles,0.0)
	);
	forAll(specieNames,lstI)
	{	
		//~ Info << "min suma "<< min(suma) <<" max(suma) "<< max(suma)<<endl;
		//~ Info << "A "<< A[lstI] <<min(T)<<max(T)<<"nuLst"<<nuLst<<endl;
		suma = suma + nuLst[lstI] * (A[lstI]*(T-Tref)  +
					B[lstI]/2.0*(Foam::pow(T,2) - Foam::pow(Tref,2))+
					C[lstI]/3.0*(Foam::pow(T,3) - Foam::pow(Tref,3))+
					D[lstI]/4.0*(Foam::pow(T,4) - Foam::pow(Tref,4)));
		//~ Info << "min suma "<< min(suma) <<" max(suma) "<< max(suma)<<endl;
	}
	//~ Info << "min suma "<< min(suma) <<" max(suma) "<< max(suma)<<endl;
	//~ Info << "A "<< A <<endl;
	//~ Info << "B "<< B <<endl;
	//~ Info << "C "<< C <<endl;
	//~ Info << "D "<< D <<endl;
	volScalarField sHatTR
	(
		sHr + suma
		
	);
	return sHatTR;
}
