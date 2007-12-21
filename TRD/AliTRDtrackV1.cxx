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
//  TRD track                                                                //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDtrackV1.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDrecoParam.h"

#include "AliESDtrack.h"

ClassImp(AliTRDtrackV1)


//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1() :
	AliTRDtrack()
	,fRecoParam(0x0)
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
AliTRDtrackV1::AliTRDtrackV1(const AliESDtrack &t) :
	AliTRDtrack(t)
	,fRecoParam(0x0)
{
  //
  // Standard constructor
  //

	//AliInfo(Form("alpha %f", GetAlpha()));
	t.GetTRDtracklets(&fTrackletIndex[0]);
	for(int ip=0; ip<6; ip++) fTracklet[ip].Reset();
}

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(const AliTRDtrackV1 &t) :
	AliTRDtrack(t)
	,fRecoParam(0x0)
{
  //
  // Copy constructor
  //

}

//_______________________________________________________________
// AliTRDtrackV1::~AliTRDtrackV1()
// {
// 	
// }

	
//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(AliTRDseedV1 *trklts, const Double_t p[5], const Double_t cov[15], Double_t x, Double_t alpha) :
	AliTRDtrack()
	,fRecoParam(0x0)
{
// The stand alone tracking constructor
// TEMPORARY !!!!!!!!!!!
// to check :
// 1. covariance matrix
// 2. dQdl calculation


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
	SetNumberOfClusters(ncl);
}

//_______________________________________________________________
Bool_t AliTRDtrackV1::CookPID()
{
// CookdEdx();  // truncated mean ... do we still need it ?

// CookdEdxTimBin(seed->GetID());
	for(int itrklt=0; itrklt<6; itrklt++){
    for (Int_t iSlice = 0; iSlice < AliESDtrack::kNSlice; iSlice++) fdEdxPlane[itrklt][iSlice] = -1.;

		if(fTrackletIndex[itrklt]<0) continue;
		fTracklet[itrklt].CookdEdx(fdEdxPlane[itrklt]);
	}

	// retrive calibration db
  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No access to calibration data");
    return kFALSE;
  }

  // Retrieve the CDB container class with the parametric detector response
  const AliTRDCalPID *pd = calibration->GetPIDObject(0/*fRecoParam->GetPIDMethod()*/);
  if (!pd) {
    AliError("No access to AliTRDCalPID object");
    return kFALSE;
  }

// CookPID(pidQ);


}


//_______________________________________________________________
Double_t AliTRDtrackV1::GetPredictedChi2(const AliTRDseedV1 *trklt) const
{
  //
  // Returns the predicted chi2
  //

  Double_t x      = trklt->GetX0();
  Double_t p[2]   = { trklt->GetYat(x)
                    , trklt->GetZat(x) };
  Double_t cov[3];
	trklt->GetCovAt(x, cov);

  return AliExternalTrackParam::GetPredictedChi2(p, cov);

}

//_______________________________________________________________
void AliTRDtrackV1::SetTracklet(AliTRDseedV1 *trklt, Int_t plane, Int_t index)
{
  //
  // Sets the tracklets
  //

	if(plane < 0 || plane >=6) return;
	fTracklet[plane]      = (*trklt);
	fTrackletIndex[plane] = index;
}

//_______________________________________________________________
Bool_t  AliTRDtrackV1::Update(const  AliTRDseedV1 *trklt, Double_t chisq)
{
  //
  // Update the track
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

  // Register info to track
//   Int_t n      = GetNumberOfClusters();
//   fIndex[n]    = index;
//   fClusters[n] = c;
  SetNumberOfClusters(GetNumberOfClusters()+trklt->GetN());
  SetChi2(GetChi2() + chisq);

	return kTRUE;
}

//_______________________________________________________________
void AliTRDtrackV1::UpdateESDdEdx(AliESDtrack *track)
{
  //
  // Update the dedx info
  //

	for (Int_t i = 0; i < AliESDtrack::kNPlane; i++) {
		for (Int_t j = 0; j < AliESDtrack::kNSlice; j++) {
			track->SetTRDsignals(fdEdxPlane[i][j], i, j);
		}
		track->SetTRDTimBin(fTimBinPlane[i], i);
	}
}
