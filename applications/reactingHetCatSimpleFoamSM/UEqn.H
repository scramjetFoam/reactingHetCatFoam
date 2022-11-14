    // Construct the Momentum equation
    
    MRF.correctBoundaryVelocity(U);

    volVectorField explRHS(-fvc::div(phi,UD));
    // explRHS -= fvc::div(phiD,U);
    explRHS -= fvc::div(phiD,UD);
    explRHS += fvc::div((rho*turbulence->nuEff())*dev2(T(fvc::grad(UD))));
    explRHS += fvc::laplacian(turbulence->nuEff()*rho,UD);

    // NOTETH: check this
    explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,UD);
    // explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,UD);
    // explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,U);

    // Info << permEff <<endl;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevTau(U)
      + fvm::div(phiD,U)
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
         explRHS +
         fvModels.source(rho, U)
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

        // pZones.addResistance(UEqn,(turbulence->nuEff*rho*pZones.D_ + (rho *mag(UD))/(2)*pZones.f())*UD);

        fvConstraints.constrain(UEqn);

        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);

        trAU = 1.0/UEqn.A();
        trAU.ref().rename("rAU");
    }
    
