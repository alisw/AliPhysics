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

// AOD track cuts for B2
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TMath.h>
#include <TString.h>
#include <TVector3.h>
#include <AliPID.h>
#include <AliAODEvent.h>
#include <AliAODMCParticle.h>
#include <AliAODTrack.h>
#include <AliAODVertex.h>
#include "AliLnAODtrackCuts.h"

ClassImp(AliLnAODtrackCuts)

AliLnAODtrackCuts::AliLnAODtrackCuts()
: TObject()
, fTrackSel("its_tpc_dca")
, fMaxDCAxy(1)
, fMaxDCAz(2)
, fMaxNSigma(3)
, fMaxEta(0.8)
, fTPCXRows(0)
, fMinTPCnClsOrXRows(70)
, fTOFmatch(0)
{
//
// constructor
}

AliLnAODtrackCuts::~AliLnAODtrackCuts()
{
//
// virtual destructor
}

void AliLnAODtrackCuts::SetSelectionCriteria(const TString& trksel)
{
//
// set track selection criteria
//
	fTrackSel = trksel;
	fTrackSel.ToLower();
	fTOFmatch = (fTrackSel.Contains("tof")) ? kTRUE : kFALSE;
}

Bool_t AliLnAODtrackCuts::IsWithinGeoAcceptance(const AliAODMCParticle* prt) const
{
//
// is particle within the geometrical acceptance?
//
	if( TMath::Abs(prt->Eta()) < fMaxEta ) return kTRUE;
	
	return kFALSE;
}

Bool_t AliLnAODtrackCuts::IsWithinGeoAcceptance(const AliAODTrack* trk) const
{
//
// is track within the geometrical acceptance?
//
//
	if( TMath::Abs(trk->Eta()) < fMaxEta ) return kTRUE;
	
	return kFALSE;
}

Bool_t AliLnAODtrackCuts::IsWithinGeoAcceptance(Double_t p[3]) const
{
//
// is track within the geometrical acceptance?
//
//
	TVector3 trk(p[0],p[1],p[2]);
	if( TMath::Abs(trk.Eta()) < fMaxEta ) return kTRUE;
	
	return kFALSE;
}

Bool_t AliLnAODtrackCuts::IsKinkDaughter(const AliAODTrack* trk) const
{
//
// is this track a kink daughter?
//
//
	AliAODVertex* vtx = trk->GetProdVertex();
	
	if(vtx == 0 ) return kFALSE;
	if(vtx->GetType() == AliAODVertex::kKink) return kTRUE;
	
	return kFALSE;
}

Double_t AliLnAODtrackCuts::GetNTPCXRowsOverFindable(const AliAODTrack* trk) const
{
//
// number of TPC crossed rows over findable
//
	if(trk->GetTPCNclsF() <= 0) return 1;
	
	return static_cast<Double_t>(trk->GetTPCNCrossedRows())/static_cast<Double_t>(trk->GetTPCNclsF());
}

Bool_t AliLnAODtrackCuts::AcceptItsTpcNSigma(const AliAODTrack* trk, Double_t b[2], Double_t bCov[3]) const
{
//
// Check ITS-TPC nsigma base cuts
//
	if(!trk->IsOn(AliAODTrack::kITSrefit)) return kFALSE;
	if(trk->GetITSNcls()<2) return kFALSE;
	//if(trk->GetITSchi2PerCluster()>36) return kFALSE;
	
	if(!trk->TestFilterBit(AliAODTrack::kTrkTPCOnly)) return kFALSE;
	if(!trk->IsOn(AliAODTrack::kTPCrefit) ) return kFALSE;
	
	if(this->IsKinkDaughter(trk)) return kFALSE;
	
	if(fTPCXRows)
	{
		if(trk->GetTPCNCrossedRows()<fMinTPCnClsOrXRows) return kFALSE;
		if(this->GetNTPCXRowsOverFindable(trk)< 0.8 ) return kFALSE;
	}
	else
	{
		if(trk->GetTPCNcls()<fMinTPCnClsOrXRows) return kFALSE;
	}
	
	//if(trk->GetTPCchi2Global()>36) return kFALSE;
	
	if(this->GetNSigmaToVertex(b, bCov) > fMaxNSigma) return kFALSE;
	
	return kTRUE;
}

Bool_t AliLnAODtrackCuts::AcceptItsTpcDCA(const AliAODTrack* trk, Double_t b[2]) const
{
//
// Check ITS-TPC-DCA base cuts
//
	if(!trk->IsOn(AliAODTrack::kITSrefit) ) return kFALSE;
	if(trk->GetITSNcls()<2) return kFALSE;
	//if(trk->GetITSchi2PerCluster()>36) return kFALSE;
	
	if(!trk->IsOn(AliAODTrack::kTPCrefit) ) return kFALSE;
	if(!trk->TestFilterBit(AliAODTrack::kTrkTPCOnly)) return kFALSE;
	
	if(this->IsKinkDaughter(trk)) return kFALSE;
	
	if(fTPCXRows)
	{
		if(trk->GetTPCNCrossedRows()<fMinTPCnClsOrXRows) return kFALSE;
		if(this->GetNTPCXRowsOverFindable(trk)< 0.8 ) return kFALSE;
	}
	else
	{
		if(trk->GetTPCNcls()<fMinTPCnClsOrXRows) return kFALSE;
	}
	
	//if(trk->GetTPCchi2Global()>36) return kFALSE;
	
	if(TMath::Abs(b[0]) > fMaxDCAxy) return kFALSE;
	if(TMath::Abs(b[1]) > fMaxDCAz) return kFALSE;
	
	return kTRUE;
}

Bool_t AliLnAODtrackCuts::AcceptItsTpcStdCut(const AliAODTrack* trk, Double_t b[2]) const
{
//
// standard cuts with very loose DCA
//
	if(!trk->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
	
	if(trk->GetTPCNcls()<70) return kFALSE;
	
	if(TMath::Abs(b[0]) > fMaxDCAxy) return kFALSE;
	if(TMath::Abs(b[1]) > fMaxDCAz) return kFALSE;
	
	return kTRUE;
}

Bool_t AliLnAODtrackCuts::AcceptTOF(const AliAODTrack* trk) const
{
//
// check TOF match signal
//
	if( !trk->IsOn(AliAODTrack::kTOFout) || !trk->IsOn(AliAODTrack::kTIME)) return kFALSE;
	//if( trk->GetIntegratedLength() < 350) return kFALSE;
	if( this->GetIntegratedLength(trk) < 350) return kFALSE;
	if( trk->GetTOFsignal() < 1e-6) return kFALSE;
	
	return kTRUE;
}

Bool_t AliLnAODtrackCuts::AcceptTrack(const AliAODTrack* trk, Double_t b[2], Double_t bCov[3]) const
{
//
// check if the tracks fulfill the track selection criteria
// from a predefined set of track cuts
//
	if(fTrackSel == "its_tpc_dca")
	{
		if(!this->AcceptItsTpcDCA(trk, b)) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_dca_spd1")
	{
		if(!this->AcceptItsTpcDCA(trk, b)) return kFALSE;
		if(!trk->HasPointOnITSLayer(0) ) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_dca_spd" || fTrackSel == "its_tpc_dca_spd2")
	{
		if(!this->AcceptItsTpcDCA(trk, b)) return kFALSE;
		if(!trk->HasPointOnITSLayer(0) && !trk->HasPointOnITSLayer(1) ) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_tof_nsigma")
	{
		if(!this->AcceptItsTpcNSigma(trk, b, bCov)) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_tof_nsigma_spd" || fTrackSel == "its_tpc_tof_nsigma_spd2")
	{
		if(!this->AcceptItsTpcNSigma(trk, b, bCov)) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
		if(!trk->HasPointOnITSLayer(0) && !trk->HasPointOnITSLayer(1) ) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_tof_dca")
	{
		if(!this->AcceptItsTpcDCA(trk, b)) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_tof_dca_spd1")
	{
		if(!this->AcceptItsTpcDCA(trk, b)) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
		if(!trk->HasPointOnITSLayer(0) ) return kFALSE;
	}
	else if(fTrackSel == "its_tpc_tof_dca_spd" || fTrackSel == "its_tpc_tof_dca_spd2")
	{
		if(!this->AcceptItsTpcDCA(trk, b)) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
		if(!trk->HasPointOnITSLayer(0) && !trk->HasPointOnITSLayer(1) ) return kFALSE;
	}
	else if(fTrackSel == "std_its_tpc_dca")
	{
		if(!this->AcceptItsTpcStdCut(trk, b)) return kFALSE;
	}
	else if(fTrackSel == "std_its_tpc_tof_dca")
	{
		if(!this->AcceptItsTpcStdCut(trk, b)) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
	}
	else if(fTrackSel == "std_its_tpc_2010" || fTrackSel == "std_its_tpc_2011")
	{
		if(!trk->TestFilterBit(AliAODTrack::kTrkGlobal)) return kFALSE;
		if(trk->GetTPCNcls()<70) return kFALSE;
	}
	else if(fTrackSel == "std_its_tpc_tof_2010" || fTrackSel == "std_its_tpc_tof_2011")
	{
		if(!trk->TestFilterBit(AliAODTrack::kTrkGlobal)) return kFALSE;
		if(trk->GetTPCNcls()<70) return kFALSE;
		if(!this->AcceptTOF(trk)) return kFALSE;
	}
	else
	{
		// default to standard cuts with tight DCA cut
		return trk->TestFilterBit(AliAODTrack::kTrkGlobal);
	}
	
	return kTRUE;
}

Double_t AliLnAODtrackCuts::GetIntegratedLength(const AliAODTrack* trk, Int_t pid) const
{
//
// track length workaround (cm)
//
	Double_t times[10];
	trk->GetIntegratedTimes(times);
	
	Double_t m = AliPID::ParticleMass(AliPID::kElectron);
	Double_t t = times[0];
	Double_t p = (pid > AliPID::kTriton) ? 2.*trk->P() : trk->P();
	Double_t c = 2.99792458e-2; // cm/ps
	Double_t beta = p/TMath::Sqrt(p*p+m*m);
	
	return c*t*beta;
}

Double_t AliLnAODtrackCuts::GetNSigmaToVertex(Double_t b[2], Double_t bCov[3]) const
{
//
// Number of sigma to the vertex (adapted from AliESDtrackCuts.cxx)
//
  //Double_t b[2];
  Double_t bRes[2];
  //Double_t bCov[3];

  //if(!this->GetImpactParameters(trk, b, b, bCov)) return 1.e+6;
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebugClass(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // -----------------------------------
  // How to get to a n-sigma cut?
  //
  // The accumulated statistics from 0 to d is
  //
  // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
  // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
  //
  // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-d**2)/2)
  // Can this be expressed in a different way?

  if (bRes[0] == 0 || bRes[1] ==0)
    return -1;

  Double_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));

  // work around precision problem
  // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
  // 1e-15 corresponds to nsigma ~ 7.7
  if (TMath::Exp(-d * d / 2) < 1e-15)
    return 1000;

  Double_t nSigma = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  return nSigma;
}
