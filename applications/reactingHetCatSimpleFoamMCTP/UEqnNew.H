    // Construct the Momentum equation
    
    MRF.correctBoundaryVelocity(U);

    // volVectorField explRHS(-fvc::div(phi,UD));
    // // // explRHS -= fvc::div(phiD,U);
    // explRHS -= fvc::div(phiD,UD);
    // explRHS += fvc::div((rho*turbulence->nuEff())*dev2(T(fvc::grad(UD))));
    // explRHS += fvc::laplacian(turbulence->nuEff()*rho,UD);

    // // // NOTETH: check this
    // explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,UD);


    // explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,U);
    // explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,UD);
    // explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,U);

    // volVectorField explRHS = - rho * turbulence->nuEff() * cmptMultiply(permEff,U);

    // volTensorField Cd = rho * turbulence->nuEff() * permEff;
    // volScalarField isoCd = 
    // explRHS = explRHS*0;
    // Info << permEff <<endl;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rhoConst, U)
      + turbulence->divDevTau(U)
    //   + fvm::div(phiD,U)
     ==
        // diffusion velocity
        // - fvc::div(phi,UD)
        // // // - fvm::div(fvc::flux(UD),U)
        // - fvm::div(linearInterpolate(rho*UD) & mesh.Sf(),U)+
        // // // - fvc::div(fvc::interpolate(rho*UD),UD)
        // - fvc::div(linearInterpolate(rho*UD) & mesh.Sf(),UD)
        // // // + fvc::div(turbulence->nuEff()*rho*dev2(T(fvc::grad(UD))))+  
        // // // + fvc::div((rho*turbulence->nuEff())*dev2(T(fvc::grad(UD)))) 
        // + fvc::div((rho*turbulence->nuEff())*dev2(T(fvc::grad(UD)))) 
        // + fvc::laplacian(turbulence->nuEff()*rho,UD)+

        //
        //  explRHS*mesh.V()/dimensionedScalar("jedna",dimVolume,1) +
        //  explRHS +
         fvModels.source(rhoConst, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    Info << max(turbulence->nuEff())<<endl;
    UEqn.relax();

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal
    
    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> tTU = tensor(I)*UEqn.A();
        pZones.addResistance(UEqn, tTU.ref());
        trTU = inv(tTU());
        trTU.ref().rename("rAU");

        fvConstraints.constrain(UEqn);

        volVectorField gradp(fvc::grad(p));

        for (int UCorr=0; UCorr<nUCorr; UCorr++)
        {
            U = trTU() & (UEqn.H() - gradp);
        }
        U.correctBoundaryConditions();

        fvConstraints.constrain(U);
    }
    
    else
    {
        pZones.addResistance(UEqn);

        // // Info << pZones.D() <<endl;

        // // const labelList& cells = mesh.cells();
        // // const cellList& cells = mesh.cells();            
        // const cellList& cells = mesh.cells();                 // Cell to face      
        // const scalarField& V = mesh.V();
        // // const volScalarField& mu = rho * turbulence->nuEff();
        // const volScalarField& mu = turbulence->nuEff();
        // const volVectorField& UN = UEqn.psi();

        // // volTensorField Cd =  tensor(vector(1,0,0),vector(0,1,0),vector(0,0,1)) & (rho * turbulence->nuEff() * permEff);

        // // volScalarField isoCd = Cd & dimensionedVector("1",dimless,vector(1,1,1));


        // forAll(cells,celli)
        // {
        // //     // const label celli = cells[i];

        // //     // const vector Cd = Cd[celli];
        // //     // const scalar isoCd = tr(Cd);
            
        //     // tensor Cd = tensor(vector(mu[celli] * permEff[celli][0],0.,0.),vector(0.,mu[celli] * permEff[celli][0],0.),vector(0.,0.,mu[celli] * permEff[celli][0])) ;
        //     tensor Cd = tensor(vector(mu[celli] * permEff[celli][0],0,0),
        //                        vector(0,mu[celli] * permEff[celli][1],0),
        //                        vector(0,0,mu[celli] * permEff[celli][2]));

        //     scalar isoCd = tr(Cd);

        //     UEqn.diag()[celli] += V[celli]*isoCd;
        //     UEqn.source()[celli] -= V[celli]*((Cd - I*isoCd) & UN[celli]);

        //     // UEqn.diag()[celli] += V[celli]*U[celli]*mu[celli]*permEff[celli][0];
        //     // UEqn.source()[celli] -= V[celli]*((Cd - I*isoCd) & UD[celli]);
        //     // UEqn.source()[celli] -= V[celli]*((Cd) & U[celli]);
        // }

        // pZones.addResistance(UEqn,(turbulence->nuEff*rho*pZones.D_ + (rho *mag(UD))/(2)*pZones.f())*UD);


        fvConstraints.constrain(UEqn);

        solve(UEqn == -fvc::grad(p));

        U.correctBoundaryConditions();
        // solve(UEqn == -fvc::grad(p) + explRHS);

        fvConstraints.constrain(U);

        trAU = 1.0/UEqn.A();
        trAU.ref().rename("rAU");
    }
    
