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
//
// Debug track
// the tree is represented as reduced events
// 
// Authors:
//   M.Fasel <M.Fasel@gsi.de>
//

#include <cstring>

#include "AliHFEreducedTrack.h"

ClassImp(AliHFEreducedTrack)


//_______________________________________
AliHFEreducedTrack::AliHFEreducedTrack():
TObject(),
  fSignedPt(0.),
  fP(0.),
  fEta(0.),
  fPhi(0.),
  fTPCmomentum(0.),
  fFilterBit(20),
  fTrackID(0),
  fMCSignedPt(0.),
  fMCP(0.),
  fMCEta(0.),
  fMCPhi(0.),
  fMCPDG(0),
  fMCMotherPdg(0),
  fMCSignal(kFALSE),
  fMCSource(5),
  fTrackStatus(8),
  fNclustersITS(0),
  fNclustersTPC(0),
  fNclustersTRD(0),
  fITSclusterMap(6),
  fITSstatusMap(6),
  fNclustersTPCPID(0),
  fNclustersTPCAll(0),
  fTPCcrossedRows(0),
  fTPCsharedClusters(0),
  fTPCclusterRatio(0.),
  fTPCclusterRatioAll(0.),
  fTRDtrackletsPID(0),
  fTRDnslices(0),
  fTRDlayer(6),
  fTRDchi2(0.),
  fTPCdEdx(0.),
  fTPCdEdxCorrected(0.),
  fTPCsigmaEl(-1000.),
  fTPCsigmaElCorrected(-1000.),
  fTOFsigmaEl(-1000.),
  fTOFmismatchProb(0.),
  fEoverP(0.),
  fEMCALsigmaEl(-1000.),
  fV0PID(kV0undef)
{
  // 
  // Default Constuctor
  //
  memset(fMCProdVtx, 0, sizeof(Double_t)*3);
  memset(fShowerShape, 0, sizeof(Double_t)*4);
  memset(fDCA, 0, sizeof(Float_t)*2);
}

//_______________________________________
AliHFEreducedTrack::AliHFEreducedTrack(const AliHFEreducedTrack &ref):
  TObject(ref),
  fSignedPt(ref.fSignedPt),
  fP(ref.fP),
  fEta(ref.fEta),
  fPhi(ref.fPhi),
  fTPCmomentum(ref.fTPCmomentum),
  fFilterBit(ref.fFilterBit),
  fTrackID(ref.fTrackID),
  fMCSignedPt(ref.fMCSignedPt),
  fMCP(ref.fMCP),
  fMCEta(ref.fMCEta),
  fMCPhi(ref.fMCPhi),
  fMCPDG(ref.fMCPDG),
  fMCMotherPdg(ref.fMCMotherPdg),
  fMCSignal(ref.fMCSignal),
  fMCSource(ref.fMCSource),
  fTrackStatus(ref.fTrackStatus),
  fNclustersITS(ref.fNclustersITS),
  fNclustersTPC(ref.fNclustersTPC),
  fNclustersTRD(ref.fNclustersTRD),
  fITSclusterMap(ref.fITSclusterMap),
  fITSstatusMap(ref.fITSstatusMap),
  fNclustersTPCPID(ref.fNclustersTPCPID),
  fNclustersTPCAll(ref.fNclustersTPCAll),
  fTPCcrossedRows(ref.fTPCcrossedRows),
  fTPCsharedClusters(ref.fTPCsharedClusters),
  fTPCclusterRatio(ref.fTPCclusterRatio),
  fTPCclusterRatioAll(ref.fTPCclusterRatioAll),
  fTRDtrackletsPID(ref.fTRDtrackletsPID),
  fTRDnslices(ref.fTRDnslices),
  fTRDlayer(ref.fTRDlayer),
  fTRDchi2(ref.fTRDchi2),
  fTPCdEdx(ref.fTPCdEdx),
  fTPCdEdxCorrected(ref.fTPCdEdxCorrected),
  fTPCsigmaEl(ref.fTPCsigmaEl),
  fTPCsigmaElCorrected(ref.fTPCsigmaElCorrected),
  fTOFsigmaEl(ref.fTOFsigmaEl),
  fTOFmismatchProb(ref.fTOFmismatchProb),
  fEoverP(ref.fEoverP),
  fEMCALsigmaEl(ref.fEMCALsigmaEl),
  fV0PID(ref.fV0PID)
{
  // 
  // Copy Constuctor
  //
  memcpy(fMCProdVtx, ref.fMCProdVtx, sizeof(Double_t) *3);
  memcpy(fShowerShape, ref.fShowerShape, sizeof(Double_t)*4);
  memcpy(fDCA, ref.fDCA, sizeof(Float_t)*2);
}

//_______________________________________
AliHFEreducedTrack &AliHFEreducedTrack::operator=(const AliHFEreducedTrack &ref){
  //
  // Assignment Operator
  //
  if(&ref != this){
    TObject::operator=(ref);
    fSignedPt = ref.fSignedPt;
    fP = ref.fP;
    fEta = ref.fEta;
    fPhi = ref.fPhi;
    fTPCmomentum = ref.fTPCmomentum;
    fFilterBit = ref.fFilterBit;
    fTrackID = ref.fTrackID;
    fMCSignedPt = ref.fMCSignedPt;
    fMCP = ref.fMCP;
    fMCEta = ref.fMCEta;
    fMCPhi = ref.fMCPhi;
    fMCPDG = ref.fMCPDG;
    fMCMotherPdg = ref.fMCMotherPdg;
    fMCSignal = ref.fMCSignal;
    fMCSource = ref.fMCSource;
    memcpy(fMCProdVtx, ref.fMCProdVtx, sizeof(Double_t) *3);
    fTrackStatus =ref.fTrackStatus;
    fNclustersITS = ref.fNclustersITS;
    fNclustersTPC = ref.fNclustersTPC;
    fNclustersTRD = ref.fNclustersTRD;
    fITSclusterMap = ref.fITSclusterMap;
    fITSstatusMap = ref.fITSstatusMap;
    fNclustersTPCPID = ref.fNclustersTPCPID;
    fNclustersTPCAll = ref.fNclustersTPCAll;
    fTPCcrossedRows = ref.fTPCcrossedRows;
    fTPCsharedClusters = ref.fTPCsharedClusters;
    fTPCclusterRatio = ref.fTPCclusterRatio;
    fTPCclusterRatioAll = ref.fTPCclusterRatioAll;
    fTRDtrackletsPID = ref.fTRDtrackletsPID;
    fTRDnslices = ref.fTRDnslices;
    fTRDlayer = ref.fTRDlayer;
    fTRDchi2 = ref.fTRDchi2;
    fTPCdEdx = ref.fTPCdEdx;
    fTPCdEdxCorrected = ref.fTPCdEdxCorrected;
    fTPCsigmaEl = ref.fTPCsigmaEl;
    fTPCsigmaElCorrected = ref.fTPCsigmaElCorrected;
    fTOFsigmaEl = ref.fTOFsigmaEl;
    fTOFmismatchProb = ref.fTOFmismatchProb;
    fEoverP = ref.fEoverP;
    fEMCALsigmaEl = ref.fEMCALsigmaEl;
    fV0PID = ref.fV0PID;
    memcpy(fShowerShape, ref.fShowerShape, sizeof(Double_t)*4);
    memcpy(fDCA, ref.fDCA, sizeof(Float_t)*2);
  }
  return *this;
}
