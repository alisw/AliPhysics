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

/* $Id$ */

/////////////////////////////////////////////////////////////////////////////
//
//  Contain parametrizations to generate atmospheric muons, and also
//  to generate single muons and muon bundles at surface level.
//
//Begin_Html
/*
<img src="picts/AliGenACORDEClass.gif">
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

#include "AliGenACORDE.h"

#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>

#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliGenACORDE)

//_____________________________________________________________________________
AliGenACORDE::AliGenACORDE()
  : AliGenerator(),
    fIpart(0),
    fCRMode(kSingleMuons),
    fCRModeName(0),
    fXwidth(0),
    fNx(1),
    fZwidth(0),
    fNz(1),
    fMuonGrid(kFALSE),
    fZenithMin(0),
    fZenithMax(0),
    fAzimuthMin(0),
    fAzimuthMax(0),
    fPRange(0),
    fPResolution(1),
    fAp(0),
    fMomentumDist(0),
    fUnfoldedMomentumDist(0),
    fZenithDist(0),
    fPDist(0)
{
  //
  // Default ctor.
  //
}

//_____________________________________________________________________________
AliGenACORDE::AliGenACORDE(Int_t npart) 
  : AliGenerator(npart),
    fIpart(kMuonMinus),
    fCRMode(kSingleMuons),
    fCRModeName(0),
    fXwidth(0),
    fNx(1),
    fZwidth(0),
    fNz(1),
    fMuonGrid(kFALSE),
    fZenithMin(0),
    fZenithMax(0),
    fAzimuthMin(0),
    fAzimuthMax(0),
    fPRange(0),
    fPResolution(1),
    fAp(0),
    fMomentumDist(0),
    fUnfoldedMomentumDist(0),
    fZenithDist(0),
    fPDist(0)
{
  //
  // Standard ctor.
  //
  fName = "ACORDE";
  fTitle = "Cosmic Muons generator";

  // Set the origin above the vertex, on the surface.
  fOrigin[0] = 0.;
  fOrigin[1] = AliACORDEConstants::Instance()->Depth(); // At the surface by default.
  fOrigin[2] = 0.;
}

//_____________________________________________________________________________
AliGenACORDE::~AliGenACORDE()
{
  //
  // Default dtor.
  //
  if ( fPDist ) {fPDist->Delete(); delete fPDist; fPDist = 0;}
  if ( fUnfoldedMomentumDist ) delete fUnfoldedMomentumDist;
  if ( fMomentumDist ) delete fMomentumDist;
  if ( fAp )           delete fAp;
  if ( fCRModeName )   delete fCRModeName;
}

//_____________________________________________________________________________
void AliGenACORDE::Generate()
{
  //
  // Generate on one trigger
  // Call the respective method inside the loop for the number
  // of tracks per trigger.

  for (Int_t i = 0; i < fNpart; i++ ) {

    if ( fCRMode == kMuonBundle ) {
      this->GenerateOneMuonBundle();

    } else if ( fCRMode == kSingleMuons ) {
      this->GenerateOneSingleMuon(kTRUE);

    } else {
      // Generate only single muons following the parametrizations
      // for atmospheric muons.
      this->GenerateOneSingleMuon();

    }

  }
}

//_____________________________________________________________________________
void AliGenACORDE::Init()
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
void AliGenACORDE::GenerateOneSingleMuon(Bool_t withFlatMomentum)
{
  //
  // Generate One Single Muon
  // This method will generate only one muon.
  // The momentum will be randomly flat distributed if
  // the paremeter "withFlatMomentum" is set to kTRUE,
  // otherwise the momentum will generate acordingly the parametrization
  // given by 
  // and adpted from Bruno Alessandro's implementation with the
  // CERNLIB to AliRoot.
  // The "withFlatMomentum" parameter also will be used to generate
  // the muons with a flat Zenithal angle distribution.
  // Do the smearing here, so that means per track.

  Float_t polar[3]= {0,0,0}; // Polarization parameters
  Float_t origin[3];
  Int_t nt;
  Float_t p[3];
  Float_t pmom, pt;
  Float_t random[6];

  // Take the azimuth random.
  Rndm(random, 2);
  Float_t azimuth = fAzimuthMin + (fAzimuthMax-fAzimuthMin)*random[0];
  Float_t zenith = fZenithMin + (fZenithMax - fZenithMin)*random[1];

  if ( withFlatMomentum ) {
    Rndm(random,3);
    if(TestBit(kMomentumRange)) {
      pmom = -( fPMin + random[0]*(fPMax - fPMin) ); // always downwards.
      pt = pmom*TMath::Sin(zenith*kDegrad);
    } else {
      pt = -( fPtMin + random[1]*(fPtMax - fPtMin)); // always downwards.
      pmom = pt/TMath::Sin(zenith*kDegrad);
    }

  } else {
    if ( fMomentumDist ) {
      pmom = -this->GetMomentum(); // Always downwards.
    } else {
      pmom = -fPMin;
    }
    zenith = this->GetZenithAngle(pmom);  // In degrees
    pt = pmom*TMath::Sin(zenith*kDegrad);
  }

  p[0] = pt*TMath::Sin(azimuth*kDegrad);
  p[1] = pmom*TMath::Cos(zenith*kDegrad);
  p[2] = pt*TMath::Cos(azimuth*kDegrad);

  // Finaly the origin, with the smearing
  Rndm(random,6);
  origin[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(zenith*kDegrad)*
    TMath::Sin(azimuth*kDegrad)
    + fOsigma[0]* TMath::Cos(2*random[0]*TMath::Pi())*TMath::Sqrt(-2*TMath::Log(random[1]));

  origin[1] = AliACORDEConstants::Instance()->Depth();

  origin[2] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(zenith*kDegrad)*
    TMath::Cos(azimuth*kDegrad)
    + fOsigma[2]* TMath::Cos(2*random[2]*TMath::Pi())*TMath::Sqrt(-2*TMath::Log(random[3]));

  // Put the track on the stack.
  PushTrack(fTrackIt,-1,fIpart,p,origin,polar,0,kPPrimary,nt);

}

//____________________________________________________________________________
void AliGenACORDE::GenerateOneMuonBundle()
{
  //
  // Generate One Muon Bundle method
  // This method will generate a bunch of muons following the
  // procedure of the AliGenScan class.
  // These muons will be generated with flat momentum.

  Float_t polar[3]= {0,0,0}; // Polarization parameters
  Float_t origin[3];
  Float_t p[3];
  Int_t nt;
  Float_t pmom;
  Float_t random[6];

  Rndm(random, 3);
  Float_t zenith = fZenithMin + (fZenithMax - fZenithMin)*random[1];
  Float_t azimuth = fAzimuthMin + (fAzimuthMax-fAzimuthMin)*random[2];
  //Float_t zenith = 10;
  //Float_t azimuth = 30;

  // Generate the kinematics a la AliGenScan (Andreas Morchs)
  Float_t dx, dz;
  if ( fNx > 0 ) {
    dx = fXwidth/fNx;
  } else {
    dx = 1e10;
    //dx = 100.;
  }

  if ( fNz > 0 ) {
    dz = fZwidth/fNz;
  } else {
    dz = 1e10;
    //dz = 100.;
  }

  origin[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(zenith*kDegrad)*
              TMath::Sin(azimuth*kDegrad);
  //origin[0] = 0.;
  origin[1] = AliACORDEConstants::Instance()->Depth();
  //origin[1] = 900;
  origin[2] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(zenith*kDegrad)*
              TMath::Cos(azimuth*kDegrad);
    //origin[2] = 0.;

  for (Int_t ix = 0; ix < fNx; ix++ ) {
    for (Int_t iz = 0; iz < fNz; iz++ ) {
      Rndm(random,6);
      origin[0]+=ix*dx+2*(random[1]-0.5)*fOsigma[0];
      origin[2]+=iz*dz+2*(random[2]-0.5)*fOsigma[2];
      if ( random[4] < 0.5 ) {
        origin[0] = -1*origin[0];
      }
      if ( random[5] < 0.5 ) {
        origin[2] = -1*origin[2];
      }

      pmom = -(fPMin + random[3] *(fPMax - fPMax) ); // Always downwards
      p[0] = TMath::Sin(zenith*kDegrad)*TMath::Sin(azimuth*kDegrad)*pmom;
      p[1] = TMath::Cos(zenith*kDegrad)*pmom;
      p[2] = TMath::Sin(zenith*kDegrad)*TMath::Cos(azimuth*kDegrad)*pmom;

      PushTrack(fTrackIt, -1, fIpart, p, origin, polar, 0, kPPrimary, nt);
    }

  }

}

//____________________________________________________________________________
void AliGenACORDE::SetGridRange(Int_t nx,Float_t xwidth, Int_t nz, Float_t zwidth)
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
void AliGenACORDE::InitApWeightFactors()
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
void AliGenACORDE::InitMomentumGeneration()
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
void AliGenACORDE::InitZenithalAngleGeneration()
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
      this->InitApWeightFactors();

      Int_t pEnd  = TMath::Abs(TMath::Nint(fPMax/fPResolution)) + 1;
      char name[26];
      char title[52];
      fPDist = new TClonesArray("TH1F", pEnd);
      TClonesArray &mom = *fPDist;
      TH1F* zenith = 0;
      Float_t weight = 0;
      for ( Int_t i = 0; i < pEnd; i++ ) {
	// Fill the distribution
	sprintf(name, "zenith%d", i+1);
	sprintf(title, "Zenith distribution, p=%f", fPMin+(Float_t)i);
	zenith = new(mom[i]) TH1F(name, title, TMath::Abs(TMath::Nint(fZenithMax-fZenithMin)), TMath::Cos(fZenithMax*TMath::Pi()/180), TMath::Cos(fZenithMin*TMath::Pi()/180));

	// Make a loop for the angle and fill the histogram for the weight
	Int_t steps = 1000;
	Float_t value = 0;
	for (Int_t j = 0; j < steps; j++ ) {
	  value = TMath::Cos(fZenithMin*TMath::Pi()/180) + (Float_t)j * ( TMath::Cos(fZenithMax*TMath::Pi()/180) - TMath::Cos(fZenithMin*TMath::Pi()/180))/1000;
	  weight = 1 + fAp->At(i)*(1 - value);
	  zenith->Fill(value, weight);
	}

      }

    } 

  }

}

//____________________________________________________________________________
Float_t AliGenACORDE::GetZenithAngle(Float_t mom) const
{

  Float_t zenith = 0.;
  // Check if you need to generate a constant zenith angle.
  if ( !fZenithDist ) {
    // Check if you have defined an array of momentum functions
    if ( fPDist ) {
      Int_t pIndex = TMath::Abs(TMath::Nint(mom));
      TH1F* cosZenithAngle = (TH1F*)fPDist->UncheckedAt(pIndex);
      Float_t tmpzenith = TMath::ACos(cosZenithAngle->GetRandom());
      // Correct the value
      zenith = kRaddeg*tmpzenith;
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

//_____________________________________________________________________________
Float_t AliGenACORDE::GetMomentum() const
{
  //
  //
  //
  return fMomentumDist->GetRandom();
}
