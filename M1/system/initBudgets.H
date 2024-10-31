if (!mesh().foundObject<volSymmTensorField>("uiuj")) 
{
    Info<<"creating uiuj"<<endl;
    volSymmTensorField* uiuj = new volSymmTensorField
    (
        IOobject
        (   
            "uiuj",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedSymmTensor(pow(dimVelocity, 2), Zero)
    );
    uiuj->store();
}
if (!mesh().foundObject<volTensorField>("djui")) 
{
    Info<<"creating djui"<<endl;
    volTensorField* djui = new volTensorField
    (
        IOobject
        (   
            "djui",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedTensor(dimVelocity/dimLength, Zero)
    );
    djui->store();
}
if (!mesh().foundObject<volScalarField>("djuidjui")) 
{
    Info<<"creating djuidjui"<<endl;
    volScalarField* djuidjui = new volScalarField
    (
        IOobject
        (   
            "djuidjui",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(pow(dimVelocity/dimLength, 2), Zero)
    );
    djuidjui->store();
}
if (!mesh().foundObject<volVectorField>("Cui")) 
{
    Info<<"creating Cui"<<endl;
    volVectorField* Cui = new volVectorField
    (
        IOobject
        (   
            "Cui",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedVector(dimAcceleration, Zero)
    );
    Cui->store();
}
if (!mesh().foundObject<volScalarField>("uiCui")) 
{
    Info<<"creating uiCui"<<endl;
    volScalarField* uiCui = new volScalarField
    (
        IOobject
        (   
            "uiCui",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimAcceleration*dimVelocity, Zero)
    );
    uiCui->store();
}
if (!mesh().foundObject<volVectorField>("Lui")) 
{
    Info<<"creating Lui"<<endl;
    volVectorField* Lui = new volVectorField
    (
        IOobject
        (   
            "Lui",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedVector(dimAcceleration/dimViscosity, Zero)
    );
    Lui->store();
}
if (!mesh().foundObject<volScalarField>("uiLui")) 
{
    Info<<"creating uiLui"<<endl;
    volScalarField* uiLui = new volScalarField
    (
        IOobject
        (   
            "uiLui",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimAcceleration*dimVelocity/dimViscosity, Zero)
    );
    uiLui->store();
}
if (!mesh().foundObject<volVectorField>("Gcp")) 
{
    Info<<"creating Gcp"<<endl;
    volVectorField* Gcp = new volVectorField
    (
        IOobject
        (   
            "Gcp",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedVector(dimAcceleration, Zero)
    );
    Gcp->store();
}
if (!mesh().foundObject<volScalarField>("uiGcp")) 
{
    Info<<"creating uiGcp"<<endl;
    volScalarField* uiGcp = new volScalarField
    (
        IOobject
        (   
            "uiGcp",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimAcceleration*dimVelocity, Zero)
    );
    uiGcp->store();
}
if (!mesh().foundObject<surfaceScalarField>("onef")) 
{
    Info<<"creating onef"<<endl;
    surfaceScalarField* onef = new surfaceScalarField
    (
        IOobject
        (   
            "onef",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimless, One)
    );
    onef->store();
}
