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
    cPphi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPphi",
            coeffDict_,
            2.0
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
            0.5
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
   prodType_
   (
       coeffDict_.lookup("prodType")
   ),
   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   eqncMu_
   (
       coeffDict_.lookup("eqncMu")
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
    
	tpphisqrt_
    (
        IOobject
        (
            "tpphi",
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
            IOobject::AUTO_WRITE
        ),
        (sqrt(k_))
    ),
	
	phiSqrt_
    (
        IOobject
        (
            "phiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    )
{

    Info<< "Made it past constructors " << endl;

    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
			if(nutType_ == "strain"){
				nut_ = 0.09*k_*Ts();
			}else{
				nut_ = cMu_*k_*tpphi_*Ts();	
			}
        }
        
        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*k_*tpphi_/epsHat_;
        }

        nut_ = min(nut_,nutRatMax_*nu());        
        nut_.correctBoundaryConditions();
        bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));       
    }
	
	
    //*************************************//	
    // Epsilon-tilda-hat
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
        cMu_.readIfPresent(coeffDict());
		cPphi_.readIfPresent(coeffDict());
		cEhmP_.readIfPresent(coeffDict());
		cEhmM_.readIfPresent(coeffDict());
		cEhmPK_.readIfPresent(coeffDict());
		cEhmPK2_.readIfPresent(coeffDict());
		cEhR_.readIfPresent(coeffDict());
		cPr_.readIfPresent(coeffDict());
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


    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, dimensionedScalar("minK", k_.dimensions(), SMALL));
        bound(epsilon_, dimensionedScalar("minEps", epsilon_.dimensions(), SMALL));
		bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), SMALL));
    }
	
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }
	
       

	   
    //*************************************//	
    // Vorticity
    //*************************************//
    
	vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);
	
	
	
	
	//*************************************//	
    // Alpha 
    //*************************************//
	
	volScalarField alpha_("alpha", 1.0/(1.0 + 1.5*tpphi_));
	

	volVectorField gradPhiSqrt_("gradPhiSqrt",fvc::grad(phiSqrt_));
	
	volVectorField gradPhi_("gradPhi", fvc::grad(phiReal()));		
	volScalarField gradgradPhi_("gradgradPhi", fvc::laplacian(DphiEff(),phiReal()));


	
	
    //*************************************//	
    // Production
    //*************************************//
 
	volScalarField S2 = magSqr(symm(uGrad_));
	volScalarField G("RASModel::G", nut_*2*S2); 
	volScalarField GdK("GdK", G/(k_ + k0_));

	if(prodType_ == "strain"){
		Info<< "Using strain production term" <<endl;
		volScalarField S2 = magSqr(symm(uGrad_));
		G = nut_*2*S2;
		tpProd_ = G/(k_ + k0_);
		GdK = G/(k_ + k0_);
	} else if(prodType_ == "mixed"){
		Info<< "Using mixed production term" <<endl; 
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + (1.0-alpha_)*cPr_*tpphi_*mag(symm(fvc::grad(U_)));
		G = tpProd_*k_;
		GdK = tpProd_;	
	} else if(prodType_ == "mixed2"){
		Info<< "Using mixed 2 production term" <<endl;
		tpProd_ = pMix_*(2*alpha_ - 1.0)*mag(tppsi_ & vorticity_) + (2.0 - pMix_)*cPr_*alpha_*tpphi_*mag(symm(fvc::grad(U_)));
		G = tpProd_*k_;
		GdK = tpProd_;
	} else if(prodType_ == "mixed3"){
		Info<< "Using mixed 3 production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + 0.33*(1.0-alpha_)*0.41*alpha_*sqrt(2.0)*mag(symm(fvc::grad(U_))) + 0.67*(1.0 - alpha_)*tpphi_*sqrt(2.0)*mag(symm(fvc::grad(U_)));
		G = tpProd_*k_;
		GdK = tpProd_;

        Info << "Max difference m3-psV: " << max(G - ((tppsi_ & vorticity_)*k_)) << endl;	

		volScalarField p1("p1",alpha_*mag(tppsi_ & vorticity_));
		volScalarField p1("p2",0.33*(1.0-alpha_)*0.41*alpha_*sqrt(2.0)*mag(symm(fvc::grad(U_))));
		volScalarField p1("p3",0.67*(1.0-alpha_)*tpphi_*sqrt(2.0)*mag(symm(fvc::grad(U_))));
		
		Info << "Min 1: " << gMin(p1) << endl; 
		Info << "Min 2: " << gMin(p2) << endl;
		Info << "Min 3: " << gMin(p3) << endl;
		Info << "Min G: " << gMin(G) << endl;
	} else if(prodType_ == "rough"){
		Info<< "Using rough production term" <<endl;
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + rPr_*(2*alpha_-1.0)*mag(tppsi_ & vorticity_) + (1.0-alpha_)*cPr_*alpha_*tpphi_*mag(symm(fvc::grad(U_)));
		G = tpProd_*k_;
		GdK = tpProd_;	
	} else{
		Info<< "Using psi-vorticity production term" <<endl;
		tpProd_ = tppsi_ & vorticity_;
		G = tpProd_*k_;
		GdK = tpProd_;		
	}
    
	tpProdSqr_ = sqr(tpProd_);
	tpProd3d_ = mag(psiReal() ^ vorticity_);
	

	
	
    //*************************************//
    // Non-constant Constants 
    //*************************************//
    
    if(eqnSigmaK_ == "true")
    {
	    sigmaK_ = 0.67 + 0.33*(tpProd_/epsHat_);
	}

    if(eqnSigmaEps_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.5*(tpProd_/epsHat_);
    }

    if(eqnSigmaPhi_ == "true")
    {
	    sigmaPhi_ = 0.21 + 0.12*(tpProd_/epsHat_);
	}

    if(eqnSigmaPsi_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.4*(tpProd_/epsHat_);
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
    //Dissipation equation
    //*************************************//
    
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       cEp1_*G*epsHat_ 
     - fvm::Sp(cEp2_*epsHat_,epsilon_)
     + cEp3_*tpProd3d_*epsHat_
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
	// Update K-related fields
    //*************************************//
    	
    kSqrt_ = sqrt(mag(k_)+k0_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), sqrt(SMALL)));
    //kSqrt_.correctBoundaryConditions();

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
    
    //Info<< "Made it past Ksqrt" <<endl;
	
	//label patchID1 = mesh_.boundaryMesh().findPatchID("FOIL_LEAD"); 
	//Info<< tppsi_.boundaryField()[patchID1] << endl;

	//label patchID2 = mesh_.boundaryMesh().findPatchID("WALL_BOTTOM"); 
    //Info<< k_.boundaryField()[patchID2] << endl;	
	

	
	
	//*************************************//	
    // Epsilon-tilda-hat
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
        epsHat_ = (1.0/(1.0 + cEhmP_*nu()*mag(gradphiSqrt)/(tpphi_*k_)))/Ts();
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "dif")
	{
        epsHat_ = (epsilon_ - 2.0*nu()*sqr(mag(gradkSqrt_)))/(k_ + k0_);
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else if(eqnEpsHat_ == "rough")
	{
        epsHat_ = 1.0/Ts();
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}
	else
	{
        Info<< "No EpsHat Model Chosen - using mod" <<endl;
	    epsHat_ = (1.0/(1.0 + (cEhmM_*nu()*mag(gradkSqrt_)/(k_+k0_))))/Ts();
	    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	}


	
	
    //*************************************//
    // Phi/K equation
    //*************************************//


    cP1eqn_ = cPphi_*(0.33 + 0.67*((tpProd_*k_)/(epsilon_ + epsilonSmall_)));
	volScalarField ruuModel("ruuModel",1.767*alpha_*k_);
	
	
    tmp<fvScalarMatrix> tpphiEqn
    (
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
      + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_)
      ==
      // cP2_*(1.0 - alpha_)*epsHat_*(ruuModel/(k_+k0_))
	    cPphi_*(2.0*alpha_-1.0)*epsHat_*tpphi_
	  + cD1_*(1.0-alpha_)*GdK*tpphi_
	  //+ cD1_*alpha_*epsHat_*tpphi_
	  + cP2_*GdK*tpphi_
	  // Prod from K eqn
      - fvm::Sp(GdK,tpphi_)
	  // Dissipation 
      - fvm::Sp((2.0*alpha_-1.0)*(epsHat_),tpphi_)
	  // Pressure diffusion 
	  //+ alpha_*tpphi_*((tppsi_ - nut_*vorticity_/k_) & gradPhi_)/(0.5*0.09*kSqrt_ + sqrt(k0_))	
	  //- fvm::Sp(cD1_*(psiReal() & gradPhi_)/(k_*kSqrt_ + k0_*sqrt(k0_)),tpphi_) 
	  - fvm::Sp((0.5 + cD4_*sigmaPhi_)*(gradPhiSqrt_ & gradPhiSqrt_)/epsHat_,tpphi_)
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
    bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), SMALL));
    }

	// Re-calculate phi/k gradient
    gradTpphi_ = fvc::grad(tpphi_);
    phiSqrt_ = sqrt(tpphi_*k_);


	
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

	  // Production and pressure strain
      - fvm::Sp((1.0 - cP2_)*tpProd_,tppsi_)
      + (1.0 - cP2_)*tpphi_*vorticity_
      - fvm::Sp(cD2_*alpha_*tpProd_,tppsi_)
      - fvm::Sp(cP1_*(1.0-alpha_)*epsHat_,tppsi_)
      
	  // Dissipation
	  + (1.0 - alpha_)*(epsHat_)*tppsi_
	  + cD3_*(2*alpha_-1.0)*tpphi_*vorticity_
	  
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
    // Calculate eddy viscosity
    //*************************************//
    
    if(solveNut_ == "true")
    {
        if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            if(nutType_ == "strain"){
				nut_ = 0.09*k_*Ts();
			}else if(nutType_ == "ruu"){
				nut_ = 0.33*alpha_*(tpphi_*k_)*Ts();
			}else if(nutType_ == "fmu"){
				nut_ = cMu_*(0.79 + 0.21*(G/epsilon_))*k_*tpphi_*Ts();	
			}else{
				nut_ = cMu_*k_*tpphi_*Ts();	
			}
        }
        
        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*k_*tpphi_/epsHat_;           
        }
               
        nut_ = min(nut_,nutRatMax_*nu());
		nut_.correctBoundaryConditions();
        bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));
    }
	


	
	
    //*************************************//   
    // Output some max values
    //*************************************//
    
    volScalarField phiActual("phiActual",tpphi_*k_);
	volScalarField psiActual("psiZ",tppsi_.component(2)*k_);
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	
	if(runTime_.outputTime())
	{   
		 volScalarField phiPstrain("phiPstrain", cP2_*(1.0 - alpha_)*epsHat_*(ruuModel/(k_+k0_)) + cP2_*GdK*tpphi_);
		 phiPstrain.write();

		 volScalarField phiDiss1("phiDiss1", -1.0*GdK*tpphi_ );
		 phiDiss1.write();
		
		 volScalarField phiDiss2("phiDiss2", -1.0*(2.0*alpha_-1.0)*epsHat_*tpphi_ );
		 phiDiss2.write();
		
		//volScalarField phiPdiff("phiPdiff", alpha_*tpphi_*((tppsi_ - nut_*vorticity_/k_) & gradPhi_)/(0.5*0.09*kSqrt_ + sqrt(k0_)));
		//phiPdiff.write();
		
		// volScalarField phiGradterm("phiGradterm", (cVv1_*nu())*(gradkSqrt_ & gradTpphi_)/(kSqrt_ + sqrt(k0_)) );
		// phiGradterm.write();
		
		volScalarField phiViscTransport("phiViscTransport", fvc::laplacian(nu(), tpphi_) + (2.0*nu())*(gradk_ & gradTpphi_)/(k_+k0_) );
		phiViscTransport.write();	

		volScalarField phiTurbTransport("phiTurbTransport", fvc::laplacian(sigmaPhi_*nut_, tpphi_) ); 
		phiTurbTransport.write();	

        volScalarField psiProd("psiProd", mag(tppsi_ & vorticity_)*k_);
        psiProd.write();	
		
        //volScalarField sProd("sProd", (1.0-alpha_)*cPr_*alpha_*tpphi_*mag(symm(fvc::grad(U_)))*k_ );
        //sProd.write();	
		
		phiActual.write();
		
		psiActual.write();
		
		alpha_.write();
		
		ruuModel.write();
		
		G.write();
		

        		
	}
	
    Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Epsilon: " << gMax(epsilon_) <<endl;
    Info<< "Max Phi: " << gMax(phiActual) << " Max Psi: " << gMax(psiActual) << " Max Production: " << gMax(G) <<endl;
    Info<< "Max 3D Production: " << gMax(tpProd3d_) << " Max uTauSquared: " << gMax(uTauSquared) << " Max vorticity: " << gMax(vorticity_) << endl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
