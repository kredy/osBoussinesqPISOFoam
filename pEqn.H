// OpenFOAM stuff to solve the equations
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

// source term as a scalar field at the faces
surfaceScalarField phig
(
    "phig",
    - rAUf*fvc::interpolate(beta*(T-TRef))
    * (g & mesh.Sf())
);

// ddtCorr reduces decoupling for unstructured meshes
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
    + rAUf*fvc::ddtCorr(U,phi)
);

/*Explicitly adjust velocities at zero gradient pressure boundaries to
  satisfy mass conservation */
adjustPhi(phiHbyA, U, p);

phiHbyA = phiHbyA + phig;
constrainPressure(p, U, phiHbyA, rAUf);


// Non-orthogonal correction
while (piso.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

    pEqn.setReference(pRefCell, pRefValue);
    
    pEqn.solve();
    
    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }    
}


// Correct velocity
U = HbyA + rAU*(fvc::reconstruct(phig/rAUf) - fvc::grad(p));

// Correct boundary conditions
U.correctBoundaryConditions();


