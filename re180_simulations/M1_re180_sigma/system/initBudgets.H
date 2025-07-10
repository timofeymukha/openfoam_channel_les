if (!mesh().foundObject<volScalarField>("nutsijsij")) 
{
    Info<<"creating nutsijsij"<<endl;
    auto* nutsijsij = new volScalarField
    (
        IOobject
        (   
            "nutsijsij",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimViscosity*pow(dimVelocity/dimLength, 2), Zero)
    );
    nutsijsij->store();
}
if (!mesh().foundObject<volScalarField>("sijsij")) 
{
    Info<<"creating sijsij"<<endl;
    auto* sijsij = new volScalarField
    (
        IOobject
        (   
            "sijsij",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(pow(dimVelocity/dimLength, 2), Zero)
    );
    sijsij->store();
}
if (!mesh().foundObject<volTensorField>("sij")) 
{
    Info<<"creating sij"<<endl;
    auto * sij = new volTensorField
    (
        IOobject
        (   
            "sij",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedTensor(dimVelocity / dimLength, Zero)
    );
    sij->store();
}
if (!mesh().foundObject<volTensorField>("nutsij")) 
{
    Info<<"creating nutsij"<<endl;
    auto * nutsij = new volTensorField
    (
        IOobject
        (   
            "nutsij",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedTensor(dimViscosity * dimVelocity / dimLength, Zero)
    );
    nutsij->store();
}
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

if (!mesh().foundObject<volVectorField>("Lui2")) 
{
    Info<<"creating Lui2"<<endl;
    volVectorField* Lui2 = new volVectorField
    (
        IOobject
        (   
            "Lui2",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedVector(dimAcceleration/dimViscosity, Zero)
    );
    Lui2->store();
}
if (!mesh().foundObject<volVectorField>("LuiSgs")) 
{
    Info<<"creating LuiSgs"<<endl;
    volVectorField* LuiSgs = new volVectorField
    (
        IOobject
        (   
            "LuiSgs",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedVector(dimViscosity*dimVelocity/dimLength/dimLength, Zero)
    );
    LuiSgs->store();
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
if (!mesh().foundObject<volScalarField>("uiLui2")) 
{
    Info<<"creating uiLui2"<<endl;
    volScalarField* uiLui2 = new volScalarField
    (
        IOobject
        (   
            "uiLui2",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimAcceleration*dimVelocity/dimViscosity, Zero)
    );
    uiLui2->store();
}
if (!mesh().foundObject<volScalarField>("uiLuiSgs")) 
{
    Info<<"creating uiLuiSgs"<<endl;
    volScalarField* uiLuiSgs = new volScalarField
    (
        IOobject
        (   
            "uiLuiSgs",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(dimAcceleration / dimTime * dimLength, Zero)
    );
    uiLuiSgs->store();
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
if (!mesh().foundObject<volScalarField>("nut")) 
{
    Info<<"creating nut"<<endl;
    volScalarField* nut = new volScalarField
    (
        IOobject
        (   
            "nut",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(), 
        dimensionedScalar(pow(dimVelocity, 2)/dimTime, Zero)
    );
    nut->store();
}
