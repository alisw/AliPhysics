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

#include "AliESDtrack.h"
#include "AliTracker.h"

#include "AliTRDtrackV1.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDrecoParam.h"

ClassImp(AliTRDtrackV1)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//  Local TRD Kalman track                                                   //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1() 
  :AliTRDtrack()
{
  //
  // Default constructor
  //
	
	for(int ip=0; ip<6; ip++){
		fTrackletIndex[ip] = -1;
		fTracklet[ip].Reset();
	}
}

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(const AliESDtrack &t) 
  :AliTRDtrack(t)
{
  //
  // Constructor from AliESDtrack
  //

	//AliInfo(Form("alpha %f", GetAlpha()));
	t.GetTRDtracklets(&fTrackletIndex[0]);
	for(int ip=0; ip<6; ip++) fTracklet[ip].Reset();
}

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(const AliTRDtrackV1 &ref) 
  :AliTRDtrack(ref)
{
  //
  // Copy constructor
  //

	for(int ip=0; ip<6; ip++){ 
		fTrackletIndex[ip] = ref.fTrackletIndex[ip];
		fTracklet[ip]      = ref.fTracklet[ip];
	}
}

//_______________________________________________________________
// AliTRDtrackV1::~AliTRDtrackV1()
// {
// 	
// }
	
//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(AliTRDseedV1 *trklts, const Double_t p[5], const Double_t cov[15]
             , Double_t x, Double_t alpha)
  :AliTRDtrack()
{
  //
  // The stand alone tracking constructor
  // TEMPORARY !!!!!!!!!!!
  // to check :
  // 1. covariance matrix
  // 2. dQdl calculation
  //

  Double_t cnv   = 1.0 / (GetBz() * kB2C);

  Double_t pp[5] = { p[0]    
                   , p[1]
                   , x*p[4] - p[2]
                   , p[3]
                   , p[4]*cnv      };

  Double_t c22 = x*x*cov[14] - 2*x*cov[12] + cov[ 5];
  Double_t c32 =   x*cov[13] -     cov[ 8];
  Double_t c20 =   x*cov[10] -     cov[ 3];
  Double_t c21 =   x*cov[11] -     cov[ 4];
  Double_t c42 =   x*cov[14] -     cov[12];

  Double_t cc[15] = { cov[ 0]
                    , cov[ 1],     cov[ 2]
                    , c20,         c21,         c22
                    , cov[ 6],     cov[ 7],     c32,     cov[ 9]
                    , cov[10]*cnv, cov[11]*cnv, c42*cnv, cov[13]*cnv, cov[14]*cnv*cnv };
  
	Set(x,alpha,pp,cc);
  Double_t s = GetSnp();
  Double_t t = GetTgl();

	Int_t ncl = 0, nclPlane; AliTRDcluster *c = 0x0;
	for(int iplane=0; iplane<6; iplane++){
		fTrackletIndex[iplane] = -1;
		fTracklet[iplane] = trklts[iplane];
		nclPlane = 0;
		for(int ic = 0; ic<AliTRDseed::knTimebins; ic++){
			if(!fTracklet[iplane].IsUsable(ic)) continue;
			if(!(c = fTracklet[iplane].GetClusters(ic))) continue;
			
			fIndex[ncl]      = fTracklet[iplane].GetIndexes(ic);
			Double_t q = TMath::Abs(c->GetQ());
			fClusters[ncl]   = c;
			// temporary !!!
			// This is not correctly. Has to be updated in the AliTRDtrackerV1::FollowBackProlonagation()
			fdQdl[ncl]       = q * (s*s < 1.) ? TMath::Sqrt((1-s*s)/(1+t*t)) : 1.;
			ncl++;
			nclPlane++;
		}
		//printf("%d N clusters plane %d [%d %d].\n", iplane, nclPlane, fTracklet[iplane].GetN2(), trklts[iplane].GetN());
	}
	//printf("N clusters in AliTRDtrackV1 %d.\n", ncl);
	SetNumberOfClusters(/*ncl*/);
}

//_______________________________________________________________
Bool_t AliTRDtrackV1::CookLabel(Float_t wrong)
{
	// set MC label for this tracklet

  Int_t s[kMAXCLUSTERSPERTRACK][2];
  for (Int_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    s[i][0] = -1;
    s[i][1] =  0;
  }

  Bool_t labelAdded;
	Int_t label;
	AliTRDcluster *c    = 0x0;
  for (Int_t ip = 0; ip < AliESDtrack::kNPlane; ip++) {
		if(fTrackletIndex[ip] < 0) continue;
		for (Int_t ic = 0; ic < AliTRDseed::knTimebins; ic++) {
			if(!(c = fTracklet[ip].GetClusters(ic))) continue;
			for (Int_t k = 0; k < 3; k++) { 
				label      = c->GetLabel(k);
				labelAdded = kFALSE; 
				Int_t j = 0;
				if (label >= 0) {
					while ((!labelAdded) && (j < kMAXCLUSTERSPERTRACK)) {
						if ((s[j][0] == label) || 
								(s[j][1] ==     0)) {
							s[j][0] = label; 
							s[j][1]++; 
							labelAdded = kTRUE;
						}
						j++;
					}
				}
			}
		}
	}

  Int_t max = 0;
  label = -123456789;
  for (Int_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    if (s[i][1] <= max) continue;
		max   = s[i][1]; 
		label = s[i][0];
  }

  if ((1. - Float_t(max)/GetNumberOfClusters()) > wrong) label = -label;

  SetLabel(label); 
	
	return kTRUE;
}

//_______________________________________________________________
Bool_t AliTRDtrackV1::CookPID()
{
  //
  // Cook the PID information
  //

// CookdEdx();  // truncated mean ... do we still need it ?

// CookdEdxTimBin(seed->GetID());
	
  // Sets the a priori probabilities
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) {
    fPID[ispec] = 1.0 / AliPID::kSPECIES;	
  }
	
	// steer PID calculation @ tracklet level
	Double_t *prob = 0x0;
	fPIDquality = 0;
	for(int itrklt=0; itrklt<AliESDtrack::kNPlane; itrklt++){
    //for (Int_t iSlice = 0; iSlice < AliESDtrack::kNSlice; iSlice++) fdEdxPlane[itrklt][iSlice] = -1.;

		if(fTrackletIndex[itrklt]<0) continue;
		if(!fTracklet[itrklt].IsOK()) continue;
		if(!(prob = fTracklet[itrklt].GetProbability())) return kFALSE;
		
		Int_t nspec = 0; // quality check of tracklet dEdx
		for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
			if(prob[ispec] < 0.) continue;
			fPID[ispec] *= prob[ispec];
			nspec++;
		}
		if(!nspec) continue;
		
		fPIDquality++;
	}
  
  // no tracklet found for PID calculations
  if(!fPIDquality) return kTRUE;

	// slot for PID calculation @ track level
	

  // normalize probabilities
  Double_t probTotal = 0.0;
  for (Int_t is = 0; is < AliPID::kSPECIES; is++) probTotal += fPID[is];


  if (probTotal <= 0.0) {
    AliWarning("The total probability over all species <= 0. This may be caused by some error in the reference data.");
    return kFALSE;
  }

  for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) fPID[iSpecies] /= probTotal;

  return kTRUE;
}

//_______________________________________________________________
Float_t AliTRDtrackV1::GetMomentum(Int_t plane) const
{
  //
  // Get the momentum at a given plane
  //

	return plane >=0 && plane < 6 && fTrackletIndex[plane] >= 0 ? fTracklet[plane].GetMomentum() : -1.;
}

//_______________________________________________________________
Double_t AliTRDtrackV1::GetPredictedChi2(const AliTRDseedV1 *trklt) const
{
  //
  // Get the predicted chi2
  //

  Double_t x      = trklt->GetX0();
  Double_t p[2]   = { trklt->GetYat(x)
                    , trklt->GetZat(x) };
  Double_t cov[3];
	trklt->GetCovAt(x, cov);

  return AliExternalTrackParam::GetPredictedChi2(p, cov);

}

//_______________________________________________________________
Bool_t AliTRDtrackV1::IsOwner() const
{
  //
  // Check whether track owns the tracklets
  //

	for (Int_t ip = 0; ip < AliESDtrack::kNPlane; ip++) {
		if(fTrackletIndex[ip] < 0) continue;
		if(!fTracklet[ip].IsOwner()) return kFALSE;
	}
	return kTRUE;
}
	
//_____________________________________________________________________________
void AliTRDtrackV1::MakeBackupTrack()
{
  //
  // Creates a backup track
  //

  if (fBackupTrack) {
    fBackupTrack->~AliTRDtrack();
    new(fBackupTrack) AliTRDtrack((AliTRDtrack&)(*this));
  }
  fBackupTrack = new AliTRDtrack((AliTRDtrack&)(*this));
}


//___________________________________________________________
void AliTRDtrackV1::SetNumberOfClusters() 
{
// Calculate the number of clusters attached to this track
	
	Int_t ncls = 0;
	for(int ip=0; ip<6; ip++){
		if(fTrackletIndex[ip] >= 0) ncls += fTracklet[ip].GetN();
	}
	AliKalmanTrack::SetNumberOfClusters(ncls);	
}

	
//_______________________________________________________________
void AliTRDtrackV1::SetOwner(Bool_t own)
{
  //
  // Toggle ownership of tracklets
  //

	for (Int_t ip = 0; ip < AliESDtrack::kNPlane; ip++) {
		if(fTrackletIndex[ip] < 0) continue;
		//AliInfo(Form("p[%d] index[%d]", ip, fTrackletIndex[ip]));
		fTracklet[ip].SetOwner(own);
	}
}

//_______________________________________________________________
void AliTRDtrackV1::SetTracklet(AliTRDseedV1 *trklt, Int_t plane, Int_t index)
{
  //
  // Set the tracklets
  //

	if(plane < 0 || plane >= AliESDtrack::kNPlane) return;
	fTracklet[plane]      = (*trklt);
	fTrackletIndex[plane] = index;
}

//_______________________________________________________________
Bool_t  AliTRDtrackV1::Update(AliTRDseedV1 *trklt, Double_t chisq)
{
  //
  // Update track parameters
  //

  Double_t x      = trklt->GetX0();
  Double_t p[2]   = { trklt->GetYat(x)
                    , trklt->GetZat(x) };
  Double_t cov[3];
	trklt->GetCovAt(x, cov);
	
// 	Print();
// 	AliInfo(Form("cov[%f %f %f]", cov[0], cov[1], cov[2]));

  if(!AliExternalTrackParam::Update(p, cov)) return kFALSE;
	//Print();

	AliTRDcluster *c = 0x0;
	Int_t ic = 0; while(!(c = trklt->GetClusters(ic))) ic++;
	AliTracker::FillResiduals(this, p, cov, c->GetVolumeId());


  // Register info to track
//   Int_t n      = GetNumberOfClusters();
//   fIndex[n]    = index;
//   fClusters[n] = c;
  SetNumberOfClusters(/*GetNumberOfClusters()+trklt->GetN()*/);
  SetChi2(GetChi2() + chisq);
	
	// update tracklet
	trklt->SetMomentum(GetP());
	trklt->SetSnp(GetSnp());
	trklt->SetTgl(GetTgl());
	return kTRUE;
}

//_______________________________________________________________
void AliTRDtrackV1::UpdateESDtrack(AliESDtrack *track)
{
  //
  // Update the ESD track
  //
	
	// copy dEdx to ESD
	Float_t *dedx = 0x0;
	for (Int_t ip = 0; ip < AliESDtrack::kNPlane; ip++) {
		if(fTrackletIndex[ip] < 0) continue;
		fTracklet[ip].CookdEdx(AliESDtrack::kNSlice);
		dedx = fTracklet[ip].GetdEdx();
		for (Int_t js = 0; js < AliESDtrack::kNSlice; js++) track->SetTRDsignals(dedx[js], ip, js);
		//track->SetTRDTimBin(fTimBinPlane[i], i);
	}

	// copy PID to ESD
	track->SetTRDpid(fPID);
	track->SetTRDpidQuality(fPIDquality);
}
