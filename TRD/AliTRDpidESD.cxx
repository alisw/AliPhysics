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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Implementation of the TRD PID class                                    //
//                                                                        //
// Assigns the electron and pion likelihoods to each ESD track.           //
// The function MakePID(AliESD *event) calculates the probability         //
// of having dedx and a maximum timbin at a given                         //
// momentum (mom) and particle type k                                     //
// from the precalculated distributions.                                  //
//                                                                        //
// Original version:                                                      //
// Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliESD.h"
#include "AliESDtrack.h"

#include "AliTRDpidESD.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "Cal/AliTRDCalPIDLQ.h"

ClassImp(AliTRDpidESD)

  Bool_t AliTRDpidESD::fCheckTrackStatus = kTRUE;
  Bool_t AliTRDpidESD::fCheckKinkStatus  = kFALSE;
  Int_t AliTRDpidESD::fMinPlane         = 0;

//_____________________________________________________________________________
AliTRDpidESD::AliTRDpidESD():TObject()
{
  //
  // Default constructor
  //


}

//_____________________________________________________________________________
AliTRDpidESD::AliTRDpidESD(const AliTRDpidESD &p):TObject(p)
{
  //
  // AliTRDpidESD copy constructor
  //

  ((AliTRDpidESD &) p).Copy(*this);

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

 ((AliTRDpidESD &) p).fCheckTrackStatus          = fCheckTrackStatus;
 ((AliTRDpidESD &) p).fCheckKinkStatus           = fCheckKinkStatus;
 ((AliTRDpidESD &) p).fMinPlane                  = fMinPlane;

}

//_____________________________________________________________________________
Int_t AliTRDpidESD::MakePID(AliESD *event)
{
  //
  // This function calculates the PID probabilities based on TRD signals
  //
  // So far this method produces probabilities based on the total charge
  // in each layer and the position of the maximum time bin in each layer.
  // In a final version this should also exploit the charge measurement in
  // the different slices of a given layer.
  //  

  Double_t p[10];
  Int_t    nSpecies  = AliPID::kSPECIES;
  Int_t    nPlanePID = 0;
  Double_t mom       = 0.0;
  Double_t probTotal = 0.0;

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliErrorGeneral("AliTRDpidESD::MakePID"
                   ,"No access to calibration data\n");
    return -1;
  }  

  // Retrieve the CDB container class with the probability distributions
  const AliTRDCalPIDLQ *pd = calibration->GetPIDLQObject();
  if (!pd) {
    AliErrorGeneral("AliTRDpidESD::MakePID"
                   ,"No access to AliTRDCalPIDLQ\n");
    return -1;
  }

  // Loop through all ESD tracks
  Int_t ntrk = event->GetNumberOfTracks();
  for (Int_t i = 0; i < ntrk; i++) {

    AliESDtrack *t = event->GetTrack(i);

    // Check the ESD track status
    if (fCheckTrackStatus) {
      if (((t->GetStatus() & AliESDtrack::kTRDout  ) == 0) &&
	  ((t->GetStatus() & AliESDtrack::kTRDrefit) == 0)) {
        continue;
      }
    }

    // Check for ESD kink tracks
    if (fCheckKinkStatus) {
      if (t->GetKinkIndex(0) != 0) {
        continue;
      }
    }

    // Skip tracks that have no TRD signal at all
    if (t->GetTRDsignal() == 0) {
      continue;
    }

    mom       = t->GetP();
    probTotal = 0.0;
    nPlanePID = 0;
    for (Int_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      p[iSpecies] = 1.0;
    }

    // Check the different detector layers
    for (Int_t iPlan = 0; iPlan < AliTRDgeometry::kNplan; iPlan++) {

      // Use the total charge in a given plane
      Double_t dedx    = t->GetTRDsignals(iPlan,-1);
      Int_t    timebin = t->GetTRDTimBin(iPlan);
      if ((dedx    >  0.0) && 
          (timebin > -1.0)) {

        nPlanePID++;

        // Get the probabilities for the different particle species
        for (Int_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

	  p[iSpecies] *= pd->GetProbability(iSpecies,mom,dedx);
	  p[iSpecies] *= pd->GetProbabilityT(iSpecies,mom,timebin);
  	  p[iSpecies] *= 100.0; // ??????????????

          probTotal   += p[iSpecies];

	}

      }

    } 

    for (Int_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      if ((probTotal >       0.0) &&
          (nPlanePID > fMinPlane)) {
        p[iSpecies] /= probTotal;
      }
      else {
        p[iSpecies] = 1.0;
      }
    }

    t->SetTRDpid(p);

  }

  return 0;

}
