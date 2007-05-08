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
// The function MakePID(AliESD *event) calculates the probability         //
// of having dedx and a maximum timbin at a given                         //
// momentum (mom) and particle type k                                     //
// from the precalculated distributions.                                  //
//                                                                        //
// Authors :                                                              //
// Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de> (Original version)//
// Alex Bercuci (a.bercuci@gsi.de)                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliESD.h"
#include "AliESDtrack.h"

#include "AliTRDpidESD.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliRun.h"
#include "AliTRDtrack.h"
#include "Cal/AliTRDCalPIDLQ.h"

ClassImp(AliTRDpidESD)

  Bool_t AliTRDpidESD::fCheckTrackStatus = kTRUE;
  Bool_t AliTRDpidESD::fCheckKinkStatus  = kFALSE;
  Int_t AliTRDpidESD::fMinPlane          = 0;

//_____________________________________________________________________________
AliTRDpidESD::AliTRDpidESD()
  :TObject()
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDpidESD::AliTRDpidESD(const AliTRDpidESD &p)
  :TObject(p)
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
  // The method produces probabilities based on the charge
  // and the position of the maximum time bin in each layer.
  // The dE/dx information can be used as global charge or 2 to 3
  // slices. Check AliTRDCalPIDLQ and AliTRDCalPIDLQRef for the actual
  // implementation.
  //
  // Author
  // Alex Bercuci (A.Bercuci@gsi.de) 2nd May 2007

	AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
	if (!calibration) {
		AliErrorGeneral("AliTRDpidESD::MakePID()"
			,"No access to calibration data");
		return -1;
	}
	
	// Retrieve the CDB container class with the probability distributions
	const AliTRDCalPIDLQ *pd = calibration->GetPIDLQObject();
	if (!pd) {
		AliErrorGeneral("AliTRDpidESD::MakePID()"
			,"No access to AliTRDCalPIDLQ");
		return -1;
	}


	// Loop through all ESD tracks
	Double_t p[10];
	AliESDtrack *t = 0x0;
	Double_t dedx[AliTRDtrack::kNslice], dEdx;
	Int_t    timebin;
	Float_t mom, length, probTotal;
	Int_t nPlanePID;
	for (Int_t i=0; i<event->GetNumberOfTracks(); i++) {
		t = event->GetTrack(i);

		// Check track
		if(!CheckTrack(t)) continue;
						
		// Skip tracks which have no TRD signal at all
		if (t->GetTRDsignal() == 0.) continue;
	
		// Loop over detector layers
		mom       = 0.; //t->GetP();
		length    = 0.;
		probTotal = 0.;
		nPlanePID    = 0;
		for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) p[iSpecies] = 1.;
		for (Int_t iPlan = 0; iPlan < AliTRDgeometry::kNplan; iPlan++) {
			// read data for track segment
			for(int iSlice=0; iSlice<AliTRDtrack::kNslice; iSlice++)
				dedx[iSlice] = t->GetTRDsignals(iPlan, iSlice);
			dEdx    = t->GetTRDsignals(iPlan, -1);
			timebin = t->GetTRDTimBin(iPlan);

			// check data
			if ((dEdx <=  0.) || (timebin <= -1.)) continue;

			// retrive kinematic info for this track segment
			if(!GetTrackSegmentKine(t, iPlan, mom, length)) continue;
			
			// this track segment has fulfilled all requierments
			nPlanePID++;
			
			// Get the probabilities for the different particle species
			for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) {
				p[iSpecies] *= pd->GetProbability(iSpecies, mom, dedx, length);
				p[iSpecies] *= pd->GetProbabilityT(iSpecies, mom, timebin);
				probTotal   += p[iSpecies];
			}
		}

		// normalize probabilities
		if(probTotal > 0.)
			for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
				if(nPlanePID > fMinPlane) p[iSpecies] /= probTotal;
				else p[iSpecies] = 1.0;


		// book PID to the track
		t->SetTRDpid(p);
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
	if (fCheckTrackStatus) {
		if (((t->GetStatus() & AliESDtrack::kTRDout  ) == 0) &&
			((t->GetStatus() & AliESDtrack::kTRDrefit) == 0)) return kFALSE;
	}

	// Check for ESD kink tracks
	if (fCheckKinkStatus && (t->GetKinkIndex(0) != 0)) return kFALSE;

	return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDpidESD::GetTrackSegmentKine(AliESDtrack *t, Int_t plan, Float_t &mom, Float_t &length)
{
  //
  // Retrive momentum "mom" and track "length" in TRD chamber from plane
  // "plan" according to information stored in AliESDtrack "t".
  // 

	if(!gAlice){
		AliErrorGeneral("AliTRDpidESD::GetTrackSegmentKine()"
		,"No gAlice object to retrive TRDgeometry and Magnetic fied  - this has to be removed in the future.");
		return kFALSE;
	}
	
	// Retrieve TRD geometry -> Maybe there is a better way to do this
	Bool_t kSelfGeom = kFALSE;
	AliTRDgeometry *TRDgeom =0x0;
	if(gAlice) TRDgeom = AliTRDgeometry::GetGeometry(gAlice->GetRunLoader());
	if(!TRDgeom){
		AliWarningGeneral("AliTRDpidESD::GetTrackSegmentKine()", "Cannot load TRD geometry from gAlice! Build a new one.\n");
		TRDgeom = new AliTRDgeometry();
		kSelfGeom = kTRUE;
	}
	const Float_t kAmHalfWidth = AliTRDgeometry::AmThick() / 2.;
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
	

	// retrive the magnetic field
	Double_t xyz0[3] = { 0., 0., 0.}, xyz1[3];
	Double_t b[3], alpha;
	gAlice->Field(xyz0,b);      // b[] is in kilo Gauss
	Float_t field = b[2] * 0.1; // Tesla

		
	// find momentum at chamber entrance and track length in chamber
	AliExternalTrackParam *param = (plan<3) ? new AliExternalTrackParam(*t->GetInnerParam()) : new AliExternalTrackParam(*t->GetOuterParam());

	param->PropagateTo(TRDgeom->GetTime0(plan)+kAmHalfWidth, field);
	param->GetXYZ(xyz0);
	alpha = param->GetAlpha();
	param->PropagateTo(TRDgeom->GetTime0(plan)-kAmHalfWidth-kDrWidth, field);
	// eliminate track segments which are crossing SM boundaries along chamber
	if(fabs(alpha-param->GetAlpha())>.01){
		delete param;
		if(kSelfGeom) delete TRDgeom;
		return kFALSE;
	}
	param->GetXYZ(xyz1);
	length = sqrt(
		(xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0])+
		(xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1])+
		(xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2])
	);
	param->GetPxPyPz(xyz1);
	mom = sqrt(xyz1[0]*xyz1[0] + xyz1[1]*xyz1[1] + xyz1[2]*xyz1[2]);
	delete param;
	if(kSelfGeom) delete TRDgeom;

	return kTRUE;
}

