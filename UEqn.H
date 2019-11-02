// Solve the momentum equation

// Interpolates the source term to the faces
fvVectorMatrix UEqn
    (
        fvm::ddt(U)
        + fvm::div(phi, U)
        - fvm::laplacian(nu, U)
    );


if (piso.momentumPredictor())
{
    solve
    (
        UEqn
        ==
        - fvc::grad(p)
        + fvc::reconstruct
        (
            -fvc::interpolate(beta*(T-TRef))*(g&mesh.Sf())
        )
    );
}
