// Construct the Momentum equation

MRF.correctBoundaryVelocity(U);

volVectorField explRHS(-fvc::div(phi,UD));
// // explRHS -= fvc::div(phiD,U);
// explRHS -= fvc::ddt(rho*UD) - fvc::ddt(phiD,UD);
explRHS -= fvc::ddt(UD)*rho - fvc::ddt(rho)*UD;
explRHS -= fvc::div(phiD,UD);
explRHS += fvc::div((rho*turbulence->nuEff())*dev2(T(fvc::grad(UD))));
explRHS += fvc::laplacian(turbulence->nuEff()*rho,UD);

// // NOTETH: check this
explRHS -= rho * turbulence->nuEff() * cmptMultiply(permEff,UD);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
    + MRF.DDt(rho, U)
    + turbulence->divDevTau(U)
    ==
    fvModels.source(rho, U)
    +explRHS
);
fvVectorMatrix& UEqn = tUEqn.ref();

Info << max(turbulence->nuEff())<<endl;
UEqn.relax();

fvConstraints.constrain(UEqn);

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

    fvConstraints.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);

    trAU = 1.0/UEqn.A();
    trAU.ref().rename("rAU");
}
