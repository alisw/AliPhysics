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
//
// This class applies the ITSsa cuts at the AOD level.
// Needed for MuonCalo pass where the FilterBit information was not properly saved.
// It contains also some quality cuts which can be modifed by user.
//
// Author: Igor Lakomov <Igor.Lakomov@cern.ch>
//

#include "AliAODITSsaTrackCuts.h"

AliAODITSsaTrackCuts::AliAODITSsaTrackCuts() : fMinNClustersITS(0), fMaxChi2PerClustersITS(0), fdcaxycut(0), fdcazcut(0), fPrimaryVertex(0)
{
//constructor
}

AliAODITSsaTrackCuts::~AliAODITSsaTrackCuts()
{
//destructor
  delete fdcaxycut;
  delete fdcazcut;
  // Do not delete, not owner  delete fPrimaryVertex;
}

Bool_t AliAODITSsaTrackCuts::AcceptTrack(const AliAODTrack* aodTrack)
{
  if (!fPrimaryVertex) {
    AliFatal("PrimaryVertex is not set! Please, use AliAODITSsaTrackCuts::ExtractAndSetPrimaryVertex(AliVEvent *event)\n");
//    return kFALSE;
  }

  if (aodTrack->IsMuonTrack()) return kFALSE; //reject Muon duplicates

  Int_t nClustersITS = 0;
  nClustersITS = aodTrack->GetITSNcls();
  if (nClustersITS<fMinNClustersITS) return kFALSE; //cut on minimum number of ITS clusters

  Float_t chi2PerClusterITS = -1;
  chi2PerClusterITS = aodTrack->GetITSchi2()/Float_t(nClustersITS);
  if (chi2PerClusterITS>fMaxChi2PerClustersITS) return kFALSE; //cut on max chi2 per ITS cluster

  if ( !( aodTrack->HasPointOnITSLayer(AliESDtrackCuts::kSPD*2) || aodTrack->HasPointOnITSLayer(AliESDtrackCuts::kSPD*2+1) ) ) return kFALSE; //at least one point in the SPD

  UInt_t status = aodTrack->GetStatus();
  if ((status&AliESDtrack::kITSrefit)==0) return kFALSE;

  if ((status & AliESDtrack::kITSin) == 0 || (status & AliESDtrack::kTPCin)) return kFALSE;  else if(!(status & AliESDtrack::kITSpureSA)) return kFALSE;

  Double_t pt = aodTrack->Pt();
  Double_t fDCAXY = CalculateDCAXY(aodTrack);
  if (fDCAXY>fdcaxycut->Eval(pt)) return kFALSE; // pt-dependent cut on DCAXY

  Double_t fDCAZ = CalculateDCAZ(aodTrack);    
  if (TMath::Abs(fDCAZ)>fdcazcut->Eval(pt)) return kFALSE; // pt-dependent cut on DCAZ

  return kTRUE; //default return
}

Double_t AliAODITSsaTrackCuts::CalculateDCAXY(const AliAODTrack* aodTrack)
{
  if (!fPrimaryVertex) {
    AliFatal("PrimaryVertex is not set! Please, use AliAODITSsaTrackCuts::ExtractAndSetPrimaryVertex(AliVEvent *event)\n");
//    return -1;
  }
  Double_t pos[3], v[3];
  fPrimaryVertex->GetXYZ(v);
  aodTrack->GetXYZ(pos);
  Double_t vDCAglobalx  = pos[0] - v[0];
  Double_t vDCAglobaly  = pos[1] - v[1];
  return TMath::Sqrt(vDCAglobalx*vDCAglobalx + vDCAglobaly*vDCAglobaly);
}

Double_t AliAODITSsaTrackCuts::CalculateDCAZ(const AliAODTrack* aodTrack)
{
  if (!fPrimaryVertex) {
    AliFatal("PrimaryVertex is not set! Please, use AliAODITSsaTrackCuts::ExtractAndSetPrimaryVertex(AliVEvent *event)\n");
//    return -1;
  }
  Double_t pos[3], vz;
  vz = fPrimaryVertex->GetZ();
  aodTrack->GetXYZ(pos);
  return pos[2] - vz;
}

void AliAODITSsaTrackCuts::SetDefaultDCAXYptdepCut2015()
{
  if(fdcaxycut){
    delete fdcaxycut;
  }
  fdcaxycut = new TFormula("fdcaxycut","0.0231+0.0315/x^1.3"); // 7*(0.0033+0.0045/pt^1.3)
}

void AliAODITSsaTrackCuts::SetUserDCAXYptdepCut(const char *formula)
{
  if(fdcaxycut){
    delete fdcaxycut;
  }
  fdcaxycut = new TFormula("fdcaxycut",formula);
}

void AliAODITSsaTrackCuts::SetDefaultDCAZptdepCut2015()
{
  if(fdcazcut){
    delete fdcazcut;
  }
  fdcazcut = new TFormula("fdcazcut","1");
}

void AliAODITSsaTrackCuts::SetUserDCAZptdepCut(const char *formula)
{
  if(fdcazcut){
    delete fdcazcut;
  }
  fdcazcut = new TFormula("fdcazcut",formula);
}

AliAODITSsaTrackCuts* AliAODITSsaTrackCuts::GetStandardAODITSsaTrackCuts2015()
{
  AliAODITSsaTrackCuts* itssatrackcuts = new AliAODITSsaTrackCuts();
  itssatrackcuts->SetDefaultDCAXYptdepCut2015();
  itssatrackcuts->SetDefaultDCAZptdepCut2015();
  itssatrackcuts->SetMinNClustersITS(4);
  itssatrackcuts->SetMaxChi2PerClustersITS(2.5);

  return itssatrackcuts;
}
