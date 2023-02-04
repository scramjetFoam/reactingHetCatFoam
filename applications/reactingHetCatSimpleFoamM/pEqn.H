{
    // thermo.psi() = Mg/univR/thermo.T();
    // const volScalarField& psi = thermo.psi();

    // Info<< max(psi)<<endl;

    // const volScalarField& psi = rho/p;
    volScalarField psi = rho/p;
    // volScalarField psi = Mg/univR/thermo.T();


    tmp<volVectorField> tHbyA;
    if (pressureImplicitPorosity)
    {
        tHbyA = constrainHbyA(trTU()&UEqn.H(), U, p);
    }
    else
    {
        tHbyA = constrainHbyA(trAU()*UEqn.H(), U, p);
    }
    volVectorField& HbyA = tHbyA.ref();

    tUEqn.clear();

    bool closedVolume = false;

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    );
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    closedVolume = adjustPhi(phiHbyA, U, p);

    // volVectorField RHSbyA(trAU()*explRHS);

    // surfaceScalarField phiRHSbyA
    // (
    //     "phiRHSbyA",
    //     fvc::interpolate(rho)*fvc::flux(RHSbyA)
    // );

    while (simple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;
        if (pressureImplicitPorosity)
        {
            tpEqn =
            (
                fvm::laplacian(rho*trTU(), p)
              + fvModels.source(psi, p, rho.name())
             ==
                // - fvc::div(phiD) +
                fvc::div(phiHbyA)
                // + fvc::div(phiRHSbyA)
            );
        }

        else
        {
            tpEqn =
            (
                fvm::laplacian(rho*trAU(), p)
              + fvModels.source(psi, p, rho.name())
             ==
                // - fvc::div(linearInterpolate(rho*UD)&mesh.Sf()) +
                // - fvc::div(phiD) +
                fvc::div(phiHbyA)
                // - fvc::div(phiD)
                // + fvc::div(phiRHSbyA)
            );
        }

        fvScalarMatrix& pEqn = tpEqn.ref();

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();
        if (simple.finalNonOrthogonalIter())
        {
            // phi = phiHbyA - pEqn.flux() + phiRHSbyA;
            // phi = phiHbyA - pEqn.flux() + phiRHSbyA;
            // phi = phiHbyA - pEqn.flux() + phiD;
            phi = phiHbyA - pEqn.flux();// - phiD;
            // phiD = linearInterpolate(rho*UD) & mesh.Sf();
            // phi = phiHbyA - pEqn.flux() + fvc::flux(trAU()*explRHS);
            // phi = phiHbyA - pEqn.flux();
        }
    }

    #include "incompressible/continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    if (pressureImplicitPorosity)
    {
        U = HbyA - (trTU() & fvc::grad(p)); //+ trAU()*explRHS;
        // U = HbyA - (trTU() & fvc::grad(p)) + trAU()*explRHS;
    }
    else
    {
        // U = HbyA - trAU()*fvc::grad(p) + trAU()*explRHS;
        // U = HbyA - trAU()*fvc::grad(p) + trAU()*explRHS;
        U = HbyA - trAU()*fvc::grad(p);
        // U = HbyA - trAU()*fvc::grad(p) + RHSbyA;
        // U = HbyA - trAU()*fvc::grad(p) + trAU()*explRHS;
    }

    U.correctBoundaryConditions();
    fvConstraints.constrain(U);

    fvConstraints.constrain(p);

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    // if (closedVolume)
    // {
    //     p += (initialMass - fvc::domainIntegrate(psi*p))
    //         /fvc::domainIntegrate(psi);
    //     p.correctBoundaryConditions();
    // }

    // rho = thermo.rho();
    // rho = p/univR/T*molMList[molMList.size()-1];
    
    // -- update gas molar mass

    // Info <<max(rho)<<endl;

    // -- update rho
    // rho.relax();
    rho = p/univR/thermo.T()*Mg;
    // thermo.rho() = rho;

    // // -- update cTot
    // cTot = p/univR/thermo.T();
    // cTot.correctBoundaryConditions();

    // Info <<max(rho)<<endl;
    
    // rho.relax();
}
