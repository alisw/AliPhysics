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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Implementation of the TRD PID class                                    //
//                                                                        //
// Assigns the electron and pion likelihoods to each ESD track.           //
// The function MakePID(AliESDEvent *event) calculates the probability    //
// of having dedx and a maximum timbin at a given                         //
// momentum (mom) and particle type k                                     //
// from the precalculated distributions.                                  //
//                                                                        //
// Authors :                                                              //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de> (orig. version) //
//   Alex Bercuci (a.bercuci@gsi.de)                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTracker.h"
#include "AliRun.h"

#include "AliTRDpidESD.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDtrack.h"
#include "Cal/AliTRDCalPID.h"

ClassImp(AliTRDpidESD)

Bool_t  AliTRDpidESD::fgCheckTrackStatus = kTRUE;
Bool_t  AliTRDpidESD::fgCheckKinkStatus  = kFALSE;
Int_t   AliTRDpidESD::fgMinPlane         = 0;

//_____________________________________________________________________________
AliTRDpidESD::AliTRDpidESD()
  :TObject()
  ,fTrack(0x0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDpidESD::AliTRDpidESD(const AliTRDpidESD &p)
  :TObject(p)
  ,fTrack(0x0)
{
  //
  // AliTRDpidESD copy constructor
  //

  ((AliTRDpidESD &) p).Copy(*this);

}

//_____________________________________________________________________________
AliTRDpidESD::~AliTRDpidESD()
{
  //
  // Destructor
  //

  if(fTrack) delete fTrack;

}

//_____________________________________________________________________________
AliTRDpidESD &AliTRDpidESD::operator=(const AliTRDpidESD &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDpidESD &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDpidESD::Copy(TObject &p) const
{
  //
  // Copy function
  //

  ((AliTRDpidESD &) p).fgCheckTrackStatus         = fgCheckTrackStatus;
  ((AliTRDpidESD &) p).fgCheckKinkStatus          = fgCheckKinkStatus;
  ((AliTRDpidESD &) p).fgMinPlane                 = fgMinPlane;
  ((AliTRDpidESD &) p).fTrack                     = 0x0;
  
}

//_____________________________________________________________________________
Int_t AliTRDpidESD::MakePID(AliESDEvent *event)
{
  //
  // This function calculates the PID probabilities based on TRD signals
  //
  // The method produces probabilities based on the charge
  // and the position of the maximum time bin in each layer.
  // The dE/dx information can be used as global charge or 2 to 3
  // slices. Check AliTRDCalPID and AliTRDCalPIDRefMaker for the actual
  // implementation.
  //
  // Author
  // Alex Bercuci (A.Bercuci@gsi.de) 2nd May 2007
  //

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliErrorGeneral("AliTRDpidESD::MakePID()"
      ,"No access to calibration data");
    return -1;
  }
  
//   AliTRDrecoParam *rec = AliTRDReconstructor::RecoParam();
//   if (!rec) {
//     AliErrorGeneral("AliTRDpidESD::MakePID()", "No TRD reco param.");
//     return 0x0;
//   }

  // Retrieve the CDB container class with the probability distributions
  const AliTRDCalPID *pd = calibration->GetPIDObject(AliTRDpidUtil::kLQ);
  if (!pd) {
    AliErrorGeneral("AliTRDpidESD::MakePID()"
      ,"No access to AliTRDCalPID");
    return -1;
  }

  // Loop through all ESD tracks
  Double_t p[10];
  AliESDtrack *t = 0x0;
  Float_t dedx[AliTRDtrack::kNslice], dEdx;
  Int_t   timebin;
  Float_t mom, length;
  Int_t   nPlanePID;
  for (Int_t i=0; i<event->GetNumberOfTracks(); i++) {
    t = event->GetTrack(i);
    
    // Check track
    if(!CheckTrack(t)) continue;

            
    // Skip tracks which have no TRD signal at all
    if (t->GetTRDsignal() == 0.) continue;
  
    // Loop over detector layers
    mom          = 0.;
    length       = 0.;
    nPlanePID    = 0;
    for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) 
                  p[iSpecies] = 1./AliPID::kSPECIES;

    for (Int_t iLayer = 0; iLayer < AliTRDgeometry::kNlayer; iLayer++) {
      // read data for track segment
      for(int iSlice=0; iSlice<AliTRDtrack::kNslice; iSlice++)
        dedx[iSlice] = t->GetTRDslice(iLayer, iSlice);
      dEdx    = t->GetTRDslice(iLayer, -1);
      timebin = t->GetTRDTimBin(iLayer);

      // check data
      if ((dEdx <=  0.) || (timebin <= -1.)) continue;

      // retrive kinematic info for this track segment
      if(!RecalculateTrackSegmentKine(t, iLayer, mom, length)){
        // information is not fully reliable especialy for length
        // estimation. To be used in the future. 
      }
      
      // this track segment has fulfilled all requierments
      nPlanePID++;

      // Get the probabilities for the different particle species
      for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) {
        p[iSpecies] *= pd->GetProbability(iSpecies, mom, dedx, length, iLayer);
        //p[iSpecies] *= pd->GetProbabilityT(iSpecies, mom, timebin);
      }
    }
    if(nPlanePID == 0) continue;
    
    // normalize probabilities
    Double_t probTotal = 0.;
    for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) probTotal   += p[iSpecies];
    if(probTotal <= 0.){
      AliWarning(Form("The total probability (%e) over all species <= 0 in ESD track %d."
                                      , probTotal, i));
      AliWarning("This may be caused by some error in reference data.");
      AliWarning("Calculation continues but results might be corrupted.");
      continue;
    }
    for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) p[iSpecies] /= probTotal;

    // book PID to the track
    t->SetTRDpid(p);
    t->SetTRDntracklets(nPlanePID<<3);
  }
  
  return 0;

}

//_____________________________________________________________________________
Bool_t AliTRDpidESD::CheckTrack(AliESDtrack *t)
{
  //
  // Check if track is eligible for PID calculations
  //
  
  // Check the ESD track status
  if (fgCheckTrackStatus) {
    if (((t->GetStatus() & AliESDtrack::kTRDout  ) == 0) &&
      ((t->GetStatus() & AliESDtrack::kTRDrefit) == 0)) return kFALSE;
  }

  // Check for ESD kink tracks
  if (fgCheckKinkStatus && (t->GetKinkIndex(0) != 0)) return kFALSE;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpidESD::RecalculateTrackSegmentKine(AliESDtrack *esd
                                              , Int_t layer
                                              , Float_t &mom
                                              , Float_t &length)
{
  //
  // Retrive momentum "mom" and track "length" in TRD chamber from plane
  // "plan" according to information stored in AliESDtrack "t".
  //
  // Origin
  // Alex Bercuci (A.Bercuci@gsi.de)	
  //

  const Float_t kAmHalfWidth = AliTRDgeometry::AmThick() / 2.;
        const Float_t kDrWidth     = AliTRDgeometry::DrThick();
  const Float_t kTime0       = AliTRDgeometry::GetTime0(layer);

  // set initial length value to chamber height 
  length = 2 * kAmHalfWidth + kDrWidth;
    
  // retrive track's outer param
  const AliExternalTrackParam *op = esd->GetOuterParam();
  if(!op){
    mom    = esd->GetP();
    return kFALSE;
  }

  AliExternalTrackParam *param = 0x0;
  if(!fTrack){
    fTrack = new AliExternalTrackParam(*op);
    param = fTrack;
  } else param = new(fTrack) AliExternalTrackParam(*op);
  
  // retrive the magnetic field
  Double_t xyz0[3];
  op->GetXYZ(xyz0);
  Float_t field = AliTracker::GetBz(xyz0); // Bz in kG at point xyz0
  Double_t s, t;

  // propagate to chamber entrance
  if(!param->PropagateTo(kTime0-kAmHalfWidth-kDrWidth, field)){
    mom    = op->GetP();
    s      = op->GetSnp();
    t      = op->GetTgl();
          if (s < 1.) length /= TMath::Sqrt((1. - s*s) / (1. + t*t));
    return kFALSE;
  }
  mom        = param->GetP();
  s = param->GetSnp();
  t = param->GetTgl();
  if (s < 1.) length    /= TMath::Sqrt((1. - s*s) / (1. + t*t));

  // check if track is crossing tracking sector by propagating to chamber exit- maybe is too much :)
  Double_t alpha = param->GetAlpha();
  if(!param->PropagateTo(kTime0+kAmHalfWidth, field)) return kFALSE;
    
  // mark track segments which are crossing SM boundaries along chamber
  if(TMath::Abs(alpha-param->GetAlpha())>.01) return kFALSE;
  
  return kTRUE;

}
