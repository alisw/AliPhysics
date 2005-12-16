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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Simple TRD Monte Carlo class                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <stdlib.h>

#include <TLorentzVector.h>
 
#include "AliRun.h"
#include "AliTRDgeometry.h"
#include "AliTRDparameter.h"
#include "AliTRDsimpleMC.h"
#include "AliTRDv1.h"
#include "AliTRDparameter.h"
#include "AliMC.h"
 
ClassImp(AliTRDsimpleMC)
 
//_____________________________________________________________________________
AliTRDsimpleMC::AliTRDsimpleMC()
{                       
  //
  // AliTRDsimpleMC default constructor
  //
       
  fMaxStep       = 0.0;
  fNStep         = 0;
  fTrack         = 0;
  fTrackPx       = 0.0;
  fTrackPy       = 0.0;
  fTrackPz       = 0.0;
  fTrackPtot     = 0.0;          
  fTrackEtot     = 0.0;
  fTrackX        = 0.0;          
  fTrackY        = 0.0;          
  fTrackZ        = 0.0;       
  fTrackStep     = 0.0;
  fTrackPid      = 0;
  fTrackCharge   = 0.0;
  fTrackMass     = 0.0;
  fTrackEntering = kFALSE;   

  fTRD           = NULL;
  fPar           = NULL;
                                        
}                                                                               

//_____________________________________________________________________________
AliTRDsimpleMC::AliTRDsimpleMC(const char *name, const char *title)
               :TVirtualMC(name,title,kFALSE)
{                       
  //
  // AliTRDsimpleMC default constructor
  //
       
  fMaxStep       = 0.0;
  fNStep         = 0;
  fTrack         = 0;
  fTrackPx       = 0.0;
  fTrackPy       = 0.0;
  fTrackPz       = 0.0;
  fTrackPtot     = 0.0;          
  fTrackEtot     = 0.0;
  fTrackX        = 0.0;          
  fTrackY        = 0.0;          
  fTrackZ        = 0.0;       
  fTrackStep     = 0.0;
  fTrackPid      = 0;
  fTrackCharge   = 0.0;
  fTrackMass     = 0.0;
  fTrackEntering = kFALSE;   

  fTRD           = NULL;
  fPar           = NULL;
                                        
}                                                                               
 
//_____________________________________________________________________________
AliTRDsimpleMC::AliTRDsimpleMC(const AliTRDsimpleMC &m):TVirtualMC(m)
{
  //
  // AliTRDsimpleMC copy constructor
  //
 
  ((AliTRDsimpleMC &) m).Copy(*this);
 
}
 
//_____________________________________________________________________________
AliTRDsimpleMC::~AliTRDsimpleMC()
{
  //
  // AliTRDsimpleMC destructor
  //
 
}                                                                               
 
//_____________________________________________________________________________
AliTRDsimpleMC &AliTRDsimpleMC::operator=(const AliTRDsimpleMC &m)
{
  //
  // Assignment operator
  //
 
  if (this != &m) ((AliTRDsimpleMC &) m).Copy(*this);
  return *this;
 
}
 
//_____________________________________________________________________________
void AliTRDsimpleMC::Copy(TObject &m) const
{
  //
  // Copy function
  //                             
                 
  ((AliTRDsimpleMC &) m).fMaxStep       = fMaxStep;
  ((AliTRDsimpleMC &) m).fNStep         = fNStep;
  ((AliTRDsimpleMC &) m).fTrack         = fTrack;
  ((AliTRDsimpleMC &) m).fTrackPx       = fTrackPx;
  ((AliTRDsimpleMC &) m).fTrackPy       = fTrackPy;
  ((AliTRDsimpleMC &) m).fTrackPz       = fTrackPz;
  ((AliTRDsimpleMC &) m).fTrackPtot     = fTrackPtot;
  ((AliTRDsimpleMC &) m).fTrackEtot     = fTrackEtot;
  ((AliTRDsimpleMC &) m).fTrackX        = fTrackX;
  ((AliTRDsimpleMC &) m).fTrackY        = fTrackY;
  ((AliTRDsimpleMC &) m).fTrackZ        = fTrackZ;
  ((AliTRDsimpleMC &) m).fTrackStep     = fTrackStep;
  ((AliTRDsimpleMC &) m).fTrackPid      = fTrackPid;
  ((AliTRDsimpleMC &) m).fTrackCharge   = fTrackCharge;
  ((AliTRDsimpleMC &) m).fTrackMass     = fTrackMass;
  ((AliTRDsimpleMC &) m).fTrackEntering = fTrackEntering;

}
                                                                                
//_____________________________________________________________________________
void AliTRDsimpleMC::NewTrack(Int_t iTrack, Int_t pdg
                             , Double_t px, Double_t py, Double_t pz)
{
  //
  // Starts a new track.
  // 

  if (!fPar) {
    fPar = new AliTRDparameter("TRDparameter","Standard TRD parameter");
  }

  if (!fTRD) {
    fTRD = (AliTRDv1 *) gAlice->GetDetector("TRD");   
    fX0  = AliTRDgeometry::GetTime0(0) - AliTRDgeometry::DrThick(); 
  }

  fTRD->ResetHits();

  fTrack         = iTrack;
  fMaxStep       = 0.0;
  fNStep         = 0;
  fTrackStep     = 0.0;
  fTrackPid      = pdg;
  fTrackEntering = kTRUE;

  switch (pdg) {
  case kPdgElectron:
    fTrackMass   =  5.11e-4;
    fTrackCharge = -1.0;
    break;
  case kPdgPion:
    fTrackMass   =  0.13957;
    fTrackCharge =  1.0;
    break;
  default:
    printf("<AliTRDsimpleMC::NewTrack> PDG code %d not implemented\n",pdg);
    break;
  };

  Double_t pTot2 = px*px + py*py + pz*pz;
  fTrackPtot = TMath::Sqrt(pTot2);
  fTrackEtot = TMath::Sqrt(pTot2 + fTrackMass*fTrackMass); 
  fTrackPx   = px;
  fTrackPy   = py;
  fTrackPz   = pz;
  
  fTrackX    = fX0;
  fTrackY    = 0.0;
  fTrackZ    = 0.0;

  gAlice->GetMCApp()->SetCurrentTrack(0);

}
                                                                                
//_____________________________________________________________________________
void AliTRDsimpleMC::ProcessEvent()
{
  //
  // Process one single track:
  //   - Determines the step size.
  //   - Calculates the track position
  //   - Calls TRD step manager.
  //

  // The stepsize from an exponential distribution
  fTrackStep = gRandom->Exp(fMaxStep);  

  if ((fTrackEntering) && (fNStep > 0)) {
    fTrackEntering = kFALSE;
  }
  fNStep++;

  // New track position
  Double_t d  = fTrackStep / fTrackPtot;
  fTrackX += fTrackPx * d;
  fTrackY += fTrackPy * d;
  fTrackZ += fTrackPz * d;

  // Call the TRD step manager
  fTRD->StepManager();  

}

//_____________________________________________________________________________
void AliTRDsimpleMC::TrackPosition(TLorentzVector& position) const
{
  //
  // Track Position
  //

  position[0] = fTrackX;
  position[1] = fTrackY;
  position[2] = fTrackZ;

}

//_____________________________________________________________________________
void AliTRDsimpleMC::TrackPosition(Double_t &x, Double_t &y, Double_t &z) const
{
  //
  // Track Position
  //

  x = fTrackX;
  y = fTrackY;
  z = fTrackZ;

}

//_____________________________________________________________________________
void AliTRDsimpleMC::TrackMomentum(TLorentzVector& momentum) const
{
  //
  // Track Momentum
  //

  momentum[0] = fTrackPx;
  momentum[1] = fTrackPy;
  momentum[2] = fTrackPz;
  momentum[3] = fTrackEtot;

}

//_____________________________________________________________________________
void AliTRDsimpleMC::TrackMomentum(Double_t &px, Double_t &py,
				   Double_t &pz, Double_t &et) const
{
  //
  // Track Momentum
  //

  px = fTrackPx;
  py = fTrackPy;
  pz = fTrackPz;
  et = fTrackEtot;

}

//_____________________________________________________________________________
Int_t AliTRDsimpleMC::VolId(const Text_t* volName) const
{
  //
  // Returns the volume IDs:
  //   1 = drift region
  //   2 = amplification region
  //   3 = drift chambers
  // 
 
  Int_t volId = -1;

  if      (strcmp(volName,"UJ00") == 0) {
    volId = kVolDrRg;
  }
  else if (strcmp(volName,"UK00") == 0) {
    volId = kVolAmRg;
  }
  else if (strcmp(volName,"UC00") == 0) {
    volId = kVolDrCh;
  }

  return volId;

}

//_____________________________________________________________________________
Int_t AliTRDsimpleMC::CurrentVolID(Int_t& copyNo) const
{
  //
  // Check for the current volume
  //

  Int_t volId = -1;

  copyNo = 0;

  // Drift region
  if      ((fTrackX-fX0) <  AliTRDgeometry::DrThick()) {
    volId = kVolDrRg;
  }
  else if ((fTrackX-fX0) < (AliTRDgeometry::DrThick() +
                            AliTRDgeometry::AmThick())) {
    volId = kVolAmRg;
  }

  return volId;

}

//_____________________________________________________________________________
const char *AliTRDsimpleMC::CurrentVolName() const
{
  //
  // Check for the current volume
  //

  const char *volName = "UA00";

  // Drift region
  if      ((fTrackX-fX0) <  AliTRDgeometry::DrThick()) {
    volName = "UJ00";
  }
  else if ((fTrackX-fX0) < (AliTRDgeometry::DrThick() +
                            AliTRDgeometry::AmThick())) {
    volName = "UK00";
  }

  return volName;

}

//_____________________________________________________________________________
Int_t AliTRDsimpleMC::CurrentVolOffID(Int_t , Int_t &copyNo) const
{
  //
  // Check for the current volume
  //

  copyNo = 1;
  return kVolDrCh;

}
