/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.1.2.1  2002/10/10 14:40:31  hristov
Updating VirtualMC to v3-09-02

Revision 1.1  2002/10/07 11:25:28  gamez
First version, generation of cosmic muons on the surface


*/

/////////////////////////////////////////////////////////////////////////////
//
//  Contain parametrizations to generate atmospheric muons, and also
//  to generate single muons and muon bundles at surface level.
//
//Begin_Html
/*
<img src="picts/AliGenCRTClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Enrique.Gamez.Flores@cern.ch">Enrique Gamez</a>.
</font>
<pre>
*/
//End_Html
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream.h>

#include "AliRun.h"
#include "AliConst.h"
#include "AliPDG.h"
#include "AliMCProcess.h"

#include "AliGenCRT.h"

ClassImp(AliGenCRT)

//_____________________________________________________________________________
AliGenCRT::AliGenCRT() : AliGenerator(-1)
{
  //
  // Default ctor.
  //
  fIpart = 0;

  fXwidth=0;
  fNx=1;
  fZwidth=0;
  fNz=1;
  fMuonGrid = kFALSE;


  // Set the origin above the vertex, on the surface.
  fOrigin[0] = 0.;
  fOrigin[1] = 0.;
  fOrigin[2] = 0.;

  fZenithMin = 0.; 
  fZenithMax = 0.;

  fAzimuthMin = 0.;
  fAzimuthMax = 0.;

  fPResolution = 1.;
  fPRange = 0.;

  fMomentumDist = 0;
  fZenithDist = 0;

  fPDist = 0;

  fAp = 0;
  fUnfoldedMomentumDist = 0;

  fCRModeName = 0;
}

//_____________________________________________________________________________
AliGenCRT::AliGenCRT(Int_t npart) 
  : AliGenerator(npart)
{
  //
  // Standard ctor.
  //
  fName = "CRT";
  fTitle = "Cosmic Muons generator";

  // Generate Muon- by default
  fIpart = kMuonMinus;

  fXwidth=0;
  fNx=1;
  fZwidth=0;
  fNz=1;
  fMuonGrid = kFALSE;

  // Set the origin above the vertex, on the surface.
  fOrigin[0] = 0.;
  fOrigin[1] = AliCRTConstants::fgDepth; // At the surface by default.
  fOrigin[2] = 0.;

  fZenithMin = 0.; // Generate veritcals by default.
  fZenithMax = 0.;

  fAzimuthMin = 0.;
  fAzimuthMax = 0.;

  fPResolution = 1.; // 1 GeV by default.
  fPRange = 0.;      // 0 GeV by default.

  fMomentumDist = 0;
  fZenithDist = 0;

  fPDist = 0;

  fAp = 0;
  fUnfoldedMomentumDist = 0;

  fCRModeName = 0;
}

//_____________________________________________________________________________
AliGenCRT::AliGenCRT(const AliGenCRT& gen) : AliGenerator()
{
  //
  // Copy ctor.
  //
  gen.Copy(*this);
}

//_____________________________________________________________________________
AliGenCRT& AliGenCRT::operator= (const AliGenCRT& gen)
{
  //
  // Asingment ctor.
  //
  gen.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
AliGenCRT::~AliGenCRT()
{
  //
  // Default dtor.
  //
  if ( fAp )           delete fAp;
  if ( fMomentumDist ) delete fMomentumDist;
  if ( fUnfoldedMomentumDist ) delete fUnfoldedMomentumDist;
  if ( fZenithDist )   delete fZenithDist;
  if ( fPDist )        delete fPDist;
}

//_____________________________________________________________________________
void AliGenCRT::Generate()
{
  //
  // Generate on one trigger
  //

  Float_t polar[3]= {0,0,0}; // Polarization parameters
  //
  Float_t origin[3];
  Float_t p[3];
  Int_t i, j, nt;
  Float_t pmom, pt;
  Float_t zenith, azimuth;
  //
  Float_t random[6];


  // Check if you need to use a distribution function for the momentum
  if ( fMomentumDist ) {
    pmom = - this->GetMomentum(); // Always downwards.

  } else {
    pmom = -fPMin;

  }
  
  zenith = this->GetZenithAngle(pmom);  // In degrees
  
  // Take the azimuth random.
  Rndm(random, 2);
  azimuth = (fAzimuthMin + (fAzimuthMax-fAzimuthMin)*random[0]); // In degrees
  
  origin[0] = AliCRTConstants::fgDepth*TMath::Tan(zenith*kDegrad)*
    TMath::Sin(azimuth*kDegrad);
  origin[1] = AliCRTConstants::fgDepth;
  origin[2] = AliCRTConstants::fgDepth*TMath::Tan(zenith*kDegrad)*
    TMath::Cos(azimuth*kDegrad);
  
  if ( fVertexSmear == kPerEvent ) {
    Rndm(random,6);
    for (j=0;j<3;j++) {
      if ( j == 1 ) {
	// Don't smear the vertical position.
	origin[j] = AliCRTConstants::fgDepth;
      } else {
	origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	  TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
      }
    }
  }


  if ( fCRMode == kSingleMuons ) {

    //
    // Generate kinematics a la AliGenBox
    //

    for(i=0;i<fNpart;i++) {
      Rndm(random,3);

      if(TestBit(kMomentumRange)) {
	pmom = -( fPMin + random[1]*(fPMax - fPMin) ); // always downwards.
	pt=pmom*TMath::Sin(zenith*kDegrad);
      } else {

      	pt= -(fPtMin+random[1]*(fPtMax-fPtMin)); // always downwards.
	pmom=pt/TMath::Sin(zenith*kDegrad);
      }
      
      p[0] = pt*TMath::Sin(azimuth*kDegrad);
      p[1] = pmom*TMath::Cos(zenith*kDegrad);
      p[2] = pt*TMath::Cos(azimuth*kDegrad);
      
      if(fVertexSmear==kPerTrack) {
	Rndm(random,6);
	for (j=0;j<3;j++) {
	  if ( j == 1 ) {
	    origin[j] = AliCRTConstants::fgDepth;
	  } else {
	    origin[j]=fOrigin[j]+fOsigma[j]*
	      TMath::Cos(2*random[2*j]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	  }
	}
      }

    }
    
    // End of generation a la AliGenBox
    
  } else if ( fCRMode == kMuonBundle ) {

    //
    // Generate kinematics a la AliGenScan
    //
    Float_t dx,dz;
    //
    if (fNx > 0) {
      dx=fXwidth/fNx;
    } else {
      dx=1e10;
    }
    
    if (fNz > 0) {
      dz=fZwidth/fNz;
    } else {
      dz=1e10;
    }

    for (Int_t ix=0; ix<fNx; ix++) {
      for (Int_t iz=0; iz<fNz; iz++){
	Rndm(random,6);
	origin[0]+=ix*dx+2*(random[0]-0.5)*fOsigma[0];
	origin[2]+=iz*dz+2*(random[1]-0.5)*fOsigma[2];	     

	pmom = -(fPMin + random[3] *(fPMax - fPMax) ); // Always downward.
	p[0] = TMath::Sin(zenith*kDegrad)*TMath::Sin(azimuth*kDegrad)*pmom;
	p[1] = TMath::Cos(zenith*kDegrad)*pmom;
	p[2] = TMath::Sin(zenith*kDegrad)*TMath::Cos(azimuth*kDegrad)*pmom;
	
      }
    }

    // End of generation a la AliGenScan

  } else if ( fCRMode == kMuonFlux ) {
    
    // Always generate the things downwards.
    p[0] = TMath::Sin(zenith*kDegrad)*TMath::Sin(azimuth*kDegrad)*pmom;
    p[1] = TMath::Cos(zenith*kDegrad)*pmom;
    p[2] = TMath::Sin(zenith*kDegrad)*TMath::Cos(azimuth*kDegrad)*pmom;

  } else {
    Error("Generate", "Generation Mode unknown!\n");
  }

  // Put the track on the stack.
  SetTrack(fTrackIt,-1,fIpart,p,origin,polar,0,kPPrimary,nt);


}

//_____________________________________________________________________________
void AliGenCRT::Init()
{
  //
  // Initialize some internal methods.
  //

  // Determine some specific data members.
  fPRange = TMath::Abs(fPMax-fPMin);

  if ( fCRMode == kSingleMuons ) {
    fCRModeName = new TString("Single Muons");
    // Initialisation, check consistency of selected ranges
    if(TestBit(kPtRange)&&TestBit(kMomentumRange)) 
      Fatal("Init","You should not set the momentum range and the pt range!");
    
    if((!TestBit(kPtRange))&&(!TestBit(kMomentumRange))) 
      Fatal("Init","You should set either the momentum or the pt range!");
    
  } else if ( fCRMode == kMuonBundle ) {
    fCRModeName = new TString("Muon Bundles");

  } else if ( fCRMode == kMuonFlux ) {
    fCRModeName = new TString("Muon Fluxes");
    // Initialize the ditribution functions.
    this->InitMomentumGeneration();
    this->InitZenithalAngleGeneration();
    
  } else {
    Fatal("Generate", "Generation Mode unknown!\n");

  }

}

//____________________________________________________________________________
void AliGenCRT::SetGridRange(Int_t nx,Float_t xwidth, Int_t nz, Float_t zwidth)
{
  //
  // Define the grid
  // This data shuold be used for Muon bundles generation.
  //
  fXwidth=xwidth;
  fNx=nx;
  fZwidth=zwidth;
  fNz=nz;

  // Print a message  about the use, if the Mode has not been set, or
  // it has to a different Mode.
  if ( fCRMode != kMuonBundle ) {
    Warning("SetRange","You have been specified a grid to generate muon bundles, but seems that you haven't choose this generation mode, or you have already select a different one");
    fMuonGrid = kTRUE;
  }
}

//____________________________________________________________________________
void AliGenCRT::InitApWeightFactors()
{
  //
  // This factors will be  to correct the zenithal angle distribution
  // acording the momentum

  //
  // Fill the array for the flux zenith angle dependence.
  // at the index 0 of fAp[] will be the "factor" if we have a muon
  // of 0 GeV.
  Float_t a[6] = {-1.61, -1.50, -1.28, -0.94, -0.61, -0.22};
  Float_t p[6] = { 0., 10., 30., 100., 300., 1000.};

  // Get the information from the momentum
  Int_t pEnd  = TMath::Abs(TMath::Nint(fPMax/fPResolution)) + 1;
  // Initialize the Array of floats to hold the a(p) factors.
  fAp = new TArrayF(pEnd);
  
  Int_t index = 0;

  for (Int_t i = 0; i < pEnd; i++ ) {
    Float_t currentP = ((Float_t)i)*fPResolution;
    if ( currentP < p[1] )                          index = 0;
    else if ( currentP >= p[1] && currentP < p[2] ) index = 1;
    else if ( currentP >= p[2] && currentP < p[3] ) index = 2;
    else if ( currentP >= p[3] && currentP < p[4] ) index = 3;
    else if ( currentP >= p[4] )                    index = 4;

    Float_t ap = (currentP -p[index])*(a[index+1] - a[index])/
                 (p[index+1] - p[index]) + a[index];
    fAp->AddAt(ap, i);
  }

}

//___________________________________________________________________________
void AliGenCRT::InitMomentumGeneration()
{
  //
  // Initialize a funtion to generate the momentum randomly
  // acording this function.
  //

  // Check if we nned to initialize the function
  if ( fPMin != fPMax ) {

    // Check also if the function have been defined yet.
    if ( !fMomentumDist ) {

      // If not, use this function
      const char* y      = "log10(x)";
      
      const char* h1Coef = "[0]*( %s*%s*%s/2 - (5*%s*%s/2) + 3*%s )";
      const char* h2Coef = "[1]*( (-2*%s*%s*%s/3) + (3*%s*%s) - 10*%s/3 + 1 )";
      const char* h3Coef = "[2]*( %s*%s*%s/6 - %s*%s/2 + %s/3 )";
      const char* s2Coef = "[3]*( %s*%s*%s/3 - 2*%s*%s + 11*%s/3 - 2 )";
      
      const char* h = "%s + %s + %s + %s";
      const char* flux = "pow(10., %s)";
      const char* normalizedFlux = "0.86*x*x*x*pow(10., %s)";
      const char* paramNames[4] = {"H1", "H2", "H3", "S1"};
      
      char buffer1[1024];
      char buffer2[1024];
      char buffer3[1024];
      char buffer4[1024];
      char buffer5[1024];
      char buffer6[1024];
      char buffer7[1024];

      sprintf(buffer1, h1Coef, y, y, y, y, y, y);
      sprintf(buffer2, h2Coef, y, y, y, y, y, y);
      sprintf(buffer3, h3Coef, y, y, y, y, y, y);
      sprintf(buffer4, s2Coef, y, y, y, y, y, y);
      
      sprintf(buffer5, h, buffer1, buffer2, buffer3, buffer4);
      
      sprintf(buffer6, flux, buffer5);
      
      fMomentumDist = new TF1("fMomentumDist", buffer6, fPMin, fPMax);
      sprintf(buffer7, normalizedFlux, buffer5);
      fUnfoldedMomentumDist = new TF1("fUnfoldedMomentumDist", buffer7, fPMin, fPMax);
      for (Int_t i = 0; i < 4; i++ ) {
	fMomentumDist->SetParName(i, paramNames[i]);
	fUnfoldedMomentumDist->SetParName(i, paramNames[i]);
      }
      
      fMomentumDist->SetParameter(0, 0.133);
      fMomentumDist->SetParameter(1, -2.521);
      fMomentumDist->SetParameter(2, -5.78);
      fMomentumDist->SetParameter(3, -2.11);

      fUnfoldedMomentumDist->SetParameter(0, 0.133);
      fUnfoldedMomentumDist->SetParameter(1, -2.521);
      fUnfoldedMomentumDist->SetParameter(2, -5.78);
      fUnfoldedMomentumDist->SetParameter(3, -2.11);
      
    }

  }

}

//____________________________________________________________________________
void AliGenCRT::InitZenithalAngleGeneration()
{
  //
  // Initalize a distribution function for the zenith angle.
  // This angle will be obtained randomly acording this function.
  // The generated angles  will been in degrees.

  // Check if we need to create the function.
  if ( fZenithMin != fZenithMax ) {

    // Check also if another function have been defined.
    if ( !fZenithDist ) {
      
      // initialize the momentum dependent coefficients, a(p) 
      InitApWeightFactors();
      
      // Define the standard function.
      char* zenithalDisributionFunction = "1 + [0]*(1 - cos(x*3.14159265358979312/180))";

      Int_t pEnd  = TMath::Abs(TMath::Nint(fPMax/fPResolution)) + 1;
      fPDist = new TClonesArray("TF1", pEnd);
      TClonesArray &angle = *fPDist;
      for ( Int_t i = 0; i < pEnd; i++ ) {
	// Fill the distribution
	TF1* zenith = new(angle[i]) TF1("zenith",zenithalDisributionFunction, fZenithMin, fZenithMax);

	// Fix the momentum dependent coefficients
	zenith->SetParName(0, "a(p)");
	zenith->SetParameter(0, fAp->At(i));

      }

    } 

  }

}

//____________________________________________________________________________
const Float_t AliGenCRT::GetZenithAngle(Float_t mom)
{

  Float_t zenith = 0.;
  // Check if you need to generate a constant zenith angle.
  if ( !fZenithDist ) {
    // Check if you have defined an array of momentum functions
    if ( fPDist ) {
      Int_t pIndex = TMath::Abs(TMath::Nint(mom));
      TF1* zenithAngle = (TF1*)fPDist->UncheckedAt(pIndex);
      zenith = zenithAngle->GetRandom();
      return zenith;
    } else {

      if ( fCRMode != kMuonFlux ) {
	// If you aren't generating muons obeying any ditribution
	// only generate a flat zenith angle, acording the input settings
	Float_t random[2];
	Rndm(random, 2);
	zenith = fZenithMin + (fZenithMax - fZenithMin)*random[0];

      } else {
	// Even if you are generating muons acording some distribution,
	// but you don't want to ...
	zenith = fZenithMin;

      }

    }
  } else {
    zenith = fZenithDist->GetRandom();
  }
     
  return zenith;
  
}

//____________________________________________________________________________
