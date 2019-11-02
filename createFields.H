Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


// Read Prandtl number
dimensionedScalar Pr
(
    "Pr",
    dimless,
    transportProperties.lookup("Pr")
    
);

// Read Rayleigh number
dimensionedScalar Ra
(
    "Ra",
    dimless,
    transportProperties.lookup("Ra")
);

// Read reference temperature
dimensionedScalar TRef
(
    "TRef",
    dimTemperature,
    transportProperties.lookup("TRef")
);


// Need better way to impliment the following two scalars
// Read temperature on the cold wall
dimensionedScalar Tcold
(
    "Tcold",
    dimTemperature,
    transportProperties.lookup("Tcold")

);

// Read temperature on the hot wall
dimensionedScalar Thot
(
    "Thot",
    dimTemperature,
    transportProperties.lookup("Thot")
);

// Temporary variable to match dimensions of beta
dimensionedScalar betatmp
(
    "betatmp",
    dimless/dimTemperature,
    scalar(1.0)
);

// Temporary variable to match dimensions of alpha and nu
dimensionedScalar alphatmp
(
    "alphatmp",
    dimViscosity,
    scalar(1.0)
);

// Temporary variable to match dimensions of T
dimensionedScalar Ttmp
(
    "Ttmp",
    dimTemperature,
    scalar(1.0)
);


// Calculate beta
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    1.0
);
beta = Pr * betatmp;

// Calculate alpha
dimensionedScalar alpha
(
    "alpha",
    dimViscosity,
    1.0
);
alpha = alphatmp * (1/sqrt(Ra));

// Calculate nu
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    1.0
);
nu = alpha * Pr;

// Read vector and scalar fields
Info<< "Reading thermophysicalProperties\n" << endl;


Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

T = ((T - Tcold)/(Thot-Tcold)) * Ttmp;

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
     ),
    mesh
);


// Calculate the velocity flux
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
   (
       "phi",
       runTime.timeName(),
       mesh,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

// Read acceleration due to gravity
#include "readGravitationalAcceleration.H"


/* To account for infinite domain with zero gradient
   pressure boundary conditions on all the walls */
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
