/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "turbulentPotential.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Hello

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPotential, 0);
addToRunTimeSelectionTable(RASModel, turbulentPotential, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> turbulentPotential::Ts() const
{ 
	if(tslimiter_ == "true")
	{
        return max(k_/(epsilon_ + epsilonSmall_), 6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)));
	}
	
    return ((k_+k0_)/(epsilon_ + epsilonSmall_));
}


tmp<volScalarField> turbulentPotential::Ls() const
{
	
	volScalarField trueL = pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
	
	if(lslimiter_ == "true")
	{
		return cL1_*max(pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_),cL2_*pow(pow3(nu())/(epsilon_ + epsilonSmall_),0.25));
	}
	
	//Info << "Max trueL: " << gMax(trueL) << " Min trueL: " << gMin(trueL) << endl;
	//Info << "Dims: " << trueL.dimensions() << endl;
	
	return pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotential::turbulentPotential
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.45
        )
    ),
    cEp2con_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.83
        )
    ),
    cEp3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cD1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            0.88
        )
    ),
    cD3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD3",
       	    coeffDict_,
            0.5
        )
    ),
    cD4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD4",
       	    coeffDict_,
            1.12
        )
    ),
    cVv1_
    (
     	dimensionedScalar::lookupOrAddToDict
        (
            "cVv1",
            coeffDict_,
            0.0
        )
    ),
    cTv1_
    (
     	dimensionedScalar::lookupOrAddToDict
        (
            "cTv1",
            coeffDict_,
            0.0
        )
    ),
    cP1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.6
        )
    ),
    cP3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.12
        )
    ),
    cP4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
        )
    ),
    cL1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL1",
            coeffDict_,
            0.36
        )
    ),
    cL2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL2",
            coeffDict_,
            85.0
        )
    ),
    eC1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC1",
            coeffDict_,
            1.4
        )
    ),
    eC2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC2",
            coeffDict_,
            0.3
        )
    ),
    eC3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC3",
            coeffDict_,
            0.0
        )
    ),
    eC4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC4",
            coeffDict_,
            0.0
        )
    ),
    eC5_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "eC5",
            coeffDict_,
            0.0
        )
    ),
    cMu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),
    betaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaK",
            coeffDict_,
            0.09
        )
    ),
    cT_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.0033
        )
    ),
    cPr_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPr",
            coeffDict_,
            0.33
        )
    ),
	cPw_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPw",
            coeffDict_,
            25.0
        )
    ),
    cEhmM_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmM",
            coeffDict_,
            10.0
        )
    ),
    cEhmP_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmP",
            coeffDict_,
            0.67
        )
    ),
    cEhmPK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmPK",
            coeffDict_,
            0.09
        )
    ),
    cEhmPK2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmPK",
            coeffDict_,
            0.67
        )
    ),
    cEhR_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhR",
            coeffDict_,
            1.0
        )
    ),
	cNF_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            1.0
        )
    ),
	rS_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "rS",
            coeffDict_,
            1e-10
        )
    ),
	pMix_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "pMix",
            coeffDict_,
            0.4
        )
    ),
	cPrK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrK",
            coeffDict_,
            0.6
        )
    ),
	cPrP_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrP",
            coeffDict_,
            1.0
        )
    ),
	rPr_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "rPr",
            coeffDict_,
            0.5
        )
    ),
	nutRatMax_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutRatMax",
            coeffDict_,
            1.0e5
        )
    ),
    sigmaKInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaKInit",
            coeffDict_,
            1.0
        )
    ),
    sigmaEpsInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEpsInit",
            coeffDict_,
            0.833
        )
    ),
    sigmaEpsVisc_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEpsVisc",
            coeffDict_,
            1.0
        )
    ),
    sigmaPhiInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPhiInit",
            coeffDict_,
            0.33
        )
    ),
    sigmaPsiInit_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPsiInit",
            coeffDict_,
            1.0
        )
    ),
	prodType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"prodType",
			coeffDict_,
			1.0
		)
	),
    nutScale_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutScale",
            coeffDict_,
            1.0
        )
    ),
    nutBlend_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutBlend",
            coeffDict_,
            1.0
        )
    ),
    psiNuFrac_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "psiNuFrac",
            coeffDict_,
            1.0
        )
    ),
    ellipticSwitch_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ellipticSwitch",
            coeffDict_,
            0.0
        )
    ),

   nutType_
   (
       coeffDict_.lookup("nutType")
   ),

   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solvePsi_
   (
       coeffDict_.lookup("solvePsi")
   ),

   solvePhi_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqnSigmaK_
   (
       coeffDict_.lookup("eqnSigmaK")
   ),

   eqnSigmaEps_
   (
       coeffDict_.lookup("eqnSigmaEps")
   ),

   eqnSigmaPhi_
   (
       coeffDict_.lookup("eqnSigmaPhi")
   ),

   eqnSigmaPsi_
   (
       coeffDict_.lookup("eqnSigmaPsi")
   ),

   eqncEp1_
   (
       coeffDict_.lookup("eqncEp1")
   ),
   
   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),
   
   timeScaleEps_
   (
       coeffDict_.lookup("timeScaleEps")
   ),
   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   lslimiter_
   (
       coeffDict_.lookup("lslimiter")
   ),
   eqncMu_
   (
       coeffDict_.lookup("eqncMu")
   ),
   phiType_
   (
       coeffDict_.lookup("phiType")
   ),
   y_
   (
   mesh_
   ),
    
	k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(k_))
    ),
    
	epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_/max(nut_))
    ),
    
	tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	tpphiSqrt_
    (
        IOobject
        (
            "tpphiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_))
    ),
    
	vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::curl(U_))
    ),
    
	tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::grad(U_))
    ),
	
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (epsilon_/(k_ + k0_))
    ),
    
	kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(k_))
    ),
	
	alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (1.0/(1.0 + 1.5*tpphi_))
    ),
	
	phiSqrt_
    (
        IOobject
        (
            "phiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_*k_))
    ),
    
	gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(kSqrt_))
    ),
    
	sigmaK_
    (
        IOobject
        (
            "sigmaK",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaKInit_)
    ),
    
	sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaEpsInit_)
    ),
    
	sigmaPhi_
    (
        IOobject
        (
            "sigmaPhi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaPhiInit_)
    ),
    
	sigmaPsi_
    (
        IOobject
        (
            "sigmaPsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        (sigmaPsiInit_)
    ),
    
	cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    
	tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((2*nut_*magSqr(symm(fvc::grad(U_)))/k_))
    ),
    
	cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (2.0*(0.5+0.5*((tpProd_*k_)/epsilon_)))
    ),
    
	dimRat_
    (
        IOobject
        (
            "dimRat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (psiReal() & psiReal())/(k_*phiReal())
    ),
    
	gradTpphi_
    (
        IOobject
        (
            "gradTpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tpphi_))
    ),
    
	gradTppsi_
    (
        IOobject
        (
            "gradTppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tppsi_))
    ),
    
	tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqr(tppsi_ & vorticity_))
    ),
    
	tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    )
{

    Info<< "Made it past constructors " << endl;

    //*************************************//	
    // Eddy viscosity - diffusion only
    //*************************************//
    if(solveNut_ == "true")
    {
		nut_ = cMu_*k_*tpphi_*Ts();	
        nut_ = min(nut_,nutRatMax_*nu());        
        nut_.correctBoundaryConditions();
        bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));       
    }
	
	
    //*************************************//	
    // Epsilon-hat
    //*************************************//
    
	if(eqnEpsHat_ == "mod")
	{
        epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "pk")
	{
        epsHat_ = epsilon_/(k_ + (cEhmPK_*nu()*mag(gradk_)/(tpphi_*kSqrt_ + sqrt(k0_))));
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "pk2")
	{
        epsHat_ = epsilon_/(k_ + (cEhmPK2_*nu()*mag(gradk_)/(phiSqrt_ + sqrt(k0_)))); 
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "phi")
	{		
		volVectorField gradphiSqrt("gradphiSqrt", fvc::grad(sqrt(tpphi_*k_))) ;
        epsHat_ = ((epsilon_)/(1.0 + cEhmP_*nu()*mag(gradphiSqrt)/(tpphi_*k_)))/(k_+k0_);
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "dif")
	{
        epsHat_ = (epsilon_ - 2.0*nu()*sqr(mag(gradkSqrt_)))/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "dify")
	{
        epsHat_ = (epsilon_ - 2.0*nu()*k_/sqr(y_))/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "rough")
	{
        epsHat_ = epsilon_/(k_+k0_);
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else
	{
        Info<< "No EpsHat Model Chosen" <<endl;
	    epsHat_ = (epsilon_)/(k_ + cEhmM_*nu()*mag(fvc::grad(kSqrt_)));
	    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}

	
	
	
    Info<< "solveK is: " <<solveK_ <<endl;
    Info<< "solveEps is: " <<solveEps_ <<endl;
    Info<< "solvePhi is: " <<solvePhi_ <<endl;
    Info<< "solvePsi is: " <<solvePsi_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotential::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_))
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotential::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> turbulentPotential::divDevReff(volVectorField& U) const
{
    return
    (
       fvc::grad(phiReal())
     + fvc::curl(psiReal())
     + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
     - fvm::laplacian(nuEff(), U)
    );
}


bool turbulentPotential::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
		cP4_.readIfPresent(coeffDict());
		cL1_.readIfPresent(coeffDict());
		cL2_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		eC1_.readIfPresent(coeffDict());
		eC2_.readIfPresent(coeffDict());
		eC3_.readIfPresent(coeffDict());
		eC4_.readIfPresent(coeffDict());
		eC5_.readIfPresent(coeffDict());
		cEhmP_.readIfPresent(coeffDict());
		cEhmM_.readIfPresent(coeffDict());
		cEhmPK_.readIfPresent(coeffDict());
		cEhmPK2_.readIfPresent(coeffDict());
		cEhR_.readIfPresent(coeffDict());
		cPr_.readIfPresent(coeffDict());
		cPrK_.readIfPresent(coeffDict());
		cPrP_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
		cD3_.readIfPresent(coeffDict());
		cD4_.readIfPresent(coeffDict());
        cVv1_.readIfPresent(coeffDict());
        cTv1_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		sigmaKInit_.readIfPresent(coeffDict());
        sigmaEpsInit_.readIfPresent(coeffDict());
        sigmaEpsVisc_.readIfPresent(coeffDict());
        sigmaPhiInit_.readIfPresent(coeffDict());
		sigmaPsiInit_.readIfPresent(coeffDict());
		nutBlend_.readIfPresent(coeffDict());
		nutScale_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void turbulentPotential::correct()
{

    //**********************************************//	
    // Bounding values not already defined by model
    //**********************************************//
	
    const dimensionedScalar eH0("minEpsHat", epsHat_.dimensions(), ROOTVSMALL);
	const dimensionedScalar nut0("minNut", nut_.dimensions(), ROOTVSMALL);
	const dimensionedScalar tph0("minTpphi", tpphi_.dimensions(), ROOTVSMALL);
	const dimensionedScalar L0("lMin", dimensionSet(0,1,0,0,0,0,0), ROOTVSMALL);

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, k0_);
        bound(epsilon_, epsilonSmall_);
		bound(tpphi_,tph0);
		bound(nut_,nut0);
    }
	
	
    RASModel::correct();

	
    if (!turbulence_)
    {
        return;
    }
	
	
    //*************************************//	
    // Timestep - for use in elliptic switch
    //*************************************//
	
    dimensionedScalar cTime = U_.mesh().time().value();
	word sMMdebug = runTime_.controlDict().lookup("showMaxMin");


    //*************************************//	
    // Vorticity and Gradient
    //*************************************//
    
	vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);
	


    //*************************************//	
    // Length and Time Scales
    //*************************************//	
	
	const volScalarField L("Length",Ls());
	const volScalarField L2("Lsqr",sqr(L));
	const volScalarField T("Time",Ts());	
		

	
	//*************************************//	
    // Misc Terms
    //*************************************//

	const volVectorField gradPhi_("gradPhi", fvc::grad(phiReal()));		
	const volScalarField gradgradPhi_("gradgradPhi", fvc::laplacian(DphiEff(),phiReal()));

	tpphiSqrt_ = sqrt(tpphi_ + ROOTVSMALL);
	const volVectorField gradTpphiSqrt("gradTpphiSqrt",fvc::grad(tpphiSqrt_));

	kSqrt_ = sqrt(mag(k_)+k0_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), sqrt(ROOTVSMALL)));

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
	
    gradTpphi_ = fvc::grad(tpphi_);
    phiSqrt_ = sqrt(tpphi_*k_ + k0_);
	const volVectorField gradPhiSqrt_("gradPhiSqrt",fvc::grad(phiSqrt_));	
	
	
    //*************************************//	
    // K Production
    //*************************************//
 
	const volScalarField S2 = 2*magSqr(dev(symm(uGrad_)));
	const volScalarField magS = sqrt(S2);
	volScalarField G("RASModel::G", nut_*S2); 
	volScalarField GdK("GdK", G/(k_ + k0_));
	const volScalarField Gnut("Gnut", nut_*S2);
	

	if(prodType_.value() == 1.0){
		Info<< "Using strain production term" <<endl;
		tpProd_ = GdK;
	} else if(prodType_.value() == 2.0){
		Info<< "Using mixed 3 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + pMix_*(1.0-alpha_)*cPrK_*alpha_*magS + (1.0 - pMix_)*(1.0 - alpha_)*cPrP_*tpphi_*magS;
		G = tpProd_*k_;
		GdK = tpProd_;	
    } else if(prodType_.value() == 3.0){
		Info<< "Using mixed 4 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + (1.0-alpha_)*GdK;
		G = tpProd_*k_;
		GdK = tpProd_;	
    } else if(prodType_.value() == 4.0){
		Info<< "Using mixed 5 production term" <<endl;
		G = alpha_*mag((tppsi_*k_) & vorticity_) + 0.27*(1.0-alpha_)*(k_ - 1.5*nut_*mag(gradPhiSqrt_))*magS;
		tpProd_ = G/(k_+k0_);
		GdK = tpProd_;			
	} else{
		Info<< "Using psi-vorticity production term" <<endl;
		tpProd_ = tppsi_ & vorticity_;
		G = tpProd_*k_;
		GdK = tpProd_;		
	}
    
	tpProdSqr_ = sqr(tpProd_);
	tpProd3d_ = mag(psiReal() ^ vorticity_);
	
	const volScalarField pOD = G/epsilon_; 
	
	
    //*************************************//
    // Non-constant Constants 
    //*************************************//
    
    if(eqnSigmaK_ == "true")
    {
	    sigmaK_ = 0.33+ 0.67*pOD;
	}

    if(eqnSigmaEps_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.33*pOD; 
    }

    if(eqnSigmaPhi_ == "true")
    {
	    sigmaPhi_ = eC4_ + eC3_*pOD;
	}

    if(eqnSigmaPsi_ == "true")
    {
	    sigmaEps_ = 0.67 + 0.4*pOD;
    }

    if(eqncEp2_ == "true")
    {
        cEp2_ = cEp2con_ - 0.16*exp(-0.25*sqr(k_)/(nu()*(epsilon_ + epsilonSmall_)));
    }
    else
    {
        cEp2_ =  cEp2con_;
    }

	

	//*************************************//	
    // Update Alpha
    //*************************************//
    
	alpha_ = 1.0/(1.0 + 1.5*tpphi_);		

	

	//*************************************//
    // Calculate eddy viscosity
    //*************************************//
    
    if(solveNut_ == "true")
    {
		nut_ = cMu_*k_*tpphi_*T;	  
        nut_ = min(nut_,nutRatMax_*nu());  
		nut_.correctBoundaryConditions();
        bound(nut_,nut0);
    }	
	
	
	
	//*************************************//	
    // Epsilon-hat
    //*************************************//
    
	if(eqnEpsHat_ == "mod")
	{
        epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
        bound(epsHat_,eH0);
	}
	else if(eqnEpsHat_ == "pk")
	{
        epsHat_ = epsilon_/(k_ + (2.0*cEhmPK_*nu()*mag(gradkSqrt_)/(tpphi_ + SMALL)));
        bound(epsHat_,eH0);
	}
	else if(eqnEpsHat_ == "dif")
	{
        epsHat_ = (epsilon_ - 2.0*nu()*sqr(mag(gradkSqrt_)))/(k_ + k0_);
        bound(epsHat_,eH0);
	}
	else if(eqnEpsHat_ == "dify")
	{
		volScalarField eH("eH", (epsilon_ - 2.0*nu()*k_/sqr(y_)));
        epsHat_ = (epsilon_ - 2.0*nu()*k_/sqr(y_))/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "rough")
	{
        epsHat_ = 1.0/Ts();
        bound(epsHat_,eH0);
	}
	else
	{
        Info<< "No EpsHat Model Chosen - using mod" <<endl;
	    epsHat_ = (1.0/(1.0 + (cEhmM_*nu()*mag(gradkSqrt_)/(k_+k0_))))/Ts();
	    bound(epsHat_,eH0);
	}
	
	

	
    //*************************************//
    //Dissipation equation
    //*************************************//
    volScalarField cEp1eqn("cEp1eqn",(cEp1_*(tpphi_/tpphi_)));
	
    if(eqncEp1_ == "true")
    {	
		cEp1eqn = min(1.6*(tpphi_/tpphi_),(cEp1_-0.1)*(1.0+0.137*alpha_));
		Info<< "Using cEps1 Equation" <<endl;
	}

    tmp<fvScalarMatrix> epsEqn   
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       cEp1eqn*G*epsHat_
     - fvm::Sp(cEp2_*epsHat_,epsilon_)
    );

    if(solveEps_ == "true")
    {
    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,epsilonSmall_);
    }
	
	
	
	
    //*************************************//
    // Turbulent kinetic energy equation
    //*************************************//
    
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/(k_+k0_),k_)
    );


    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,k0_);
    }

	
	
	
	
    //*************************************//
    // Phi/K equation 
    //*************************************//


    cP1eqn_ = cP1_*(0.33 + 0.67*((tpProd_*k_)/(epsilon_ + epsilonSmall_)));
	

    tmp<fvScalarMatrix> tpphiEqn
    (
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
      + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_)
      ==
	  // Pressure Strain Slow
	    cP1_*(2.0*alpha_-1.0)*epsHat_*tpphi_
	  // Pressure Strain Fast
	  + cP2_*tpphi_*GdK
	  + cP2_*cD2_*(1.0-alpha_)*tpphi_*GdK  
	  //+ cP3_*epsHat_*tpphi_
	  // Prod from K eqn
      - fvm::Sp(GdK,tpphi_)
	  // Dissipation 
      - fvm::Sp((2.0*alpha_-1.0)*epsHat_,tpphi_)
	  // Pressure diffusion  
	  - fvm::Sp((0.5 + cD4_*sigmaPhi_)*(1.0/(1.0 + cPw_/(reTau())))*(gradPhiSqrt_ & gradPhiSqrt_)*T,tpphi_)
	  // Extra diffusion terms
      + (cVv1_*nu())*(gradk_ & gradTpphi_)/(k_+k0_)
	  - fvm::SuSp((cTv1_*nut_)*(gradk_ & gradTpphi_)/(tpphi_*k_ + k0_),tpphi_)
	  // Transition
      + cT_*tpProd_*sqrt((((nu()/1000.0)+nut_)/nu()))
    );

    if(solvePhi_ == "true")
    {
    tpphiEqn().relax();
    solve(tpphiEqn);
    bound(tpphi_,tph0);
    }
	



	
    //*************************************//   
    // Psi Specific Constants
    //*************************************//
	const volScalarField psiProd("psiProd", (tppsi_ & vorticity_)); 
	
	
    //*************************************//   
    // Psi Equation
    //*************************************//
    
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)

      ==

	  // Production
	    (1.0-cP2_)*tpphi_*vorticity_
      - fvm::Sp(tpProd_,tppsi_)
	  + epsHat_*tppsi_
	  
	  // Fast Pressure strain
	  //+ (1.0 - cP2_)*(tppsi_ & vorticity_)*tppsi_
      //- fvm::Sp(cP2_*(2.0*alpha_+0.5)*tpProd_,tppsi_) 
	  - fvm::Sp(cP2_*tpProd_,tppsi_)
	  
	  // Slow Pressure Strain + Dissipation
      - fvm::Sp(cP1_*(1.0-alpha_)*epsHat_,tppsi_)
      
	  // Dissipation
	  - fvm::Sp(cD1_*alpha_*epsHat_,tppsi_)
	  //+ cD3_*(2*alpha_-1.0)*tpphi_*vorticity_
	  
	  // Gradients
      + (cTv1_*nut_)*(gradk_ & gradTppsi_)/(k_+k0_)

	  // Transition Term
      + cT_*sqrt((((nu()/100.0)+nut_)/nu()))*vorticity_
    );

    if(solvePsi_ == "true") 
    {
    tppsiEqn().relax();
    solve(tppsiEqn);
    }
	
    // Re-calculate psi/k gradient
    gradTppsi_ = fvc::grad(tppsi_);
	




	
    //*************************************//   
    // Output some max values
    //*************************************//
	
	if(sMMdebug == "true")
	{
    
    volScalarField phiActual("phiActual",tpphi_*k_);
	volScalarField psiActual("psiZ",tppsi_.component(2)*k_);
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	
	Info << "Max cEp1: " << max(cEp1eqn) << " Min cEp1: " << min(cEp1eqn) << endl; 
    Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Epsilon: " << gMax(epsilon_) <<endl;
    Info<< "Max Phi: " << gMax(phiActual) << " Max Psi: " << gMax(psiActual) << " Max G: " << gMax(G) << " Max Gnut: " << gMax(Gnut) <<endl;
    Info<< "Max uTauSquared: " << gMax(uTauSquared) << " Max vorticity: " << gMax(vorticity_) << endl;
	
	}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
