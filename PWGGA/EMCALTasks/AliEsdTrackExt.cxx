// $Id$
//
// Modified track class to be able to store cached quantities.
//
// Author: C.Loizides

#include "AliEsdTrackExt.h"
#include <TVector3.h>
#include "AliESDEvent.h"

//_________________________________________________________________________________________________
AliEsdTrackExt::AliEsdTrackExt() 
  : AliESDtrack(), fEmcEta(-10), fEmcPhi(-10), fNCrossedRows(-10), fChi2TPCConstrainedVsGlobal(-10)
{
  // Default constructor.
}

//_________________________________________________________________________________________________
AliEsdTrackExt::AliEsdTrackExt(const AliESDtrack &t)
  : AliESDtrack(t), fEmcEta(-10), fEmcPhi(-10), fNCrossedRows(-10), fChi2TPCConstrainedVsGlobal(-10)
{
  // Constructor.

  const AliExternalTrackParam *outp = GetOuterParam();
  if (outp&&IsEMCAL()) {
    Double_t trkPos[3] = {0.,0.,0.};
    if (outp->GetXYZ(trkPos)) {
      TVector3 vec(trkPos[0],trkPos[1],trkPos[2]);
      Double_t veta = vec.Eta();
      Double_t vphi = vec.Phi();
      if(vphi<0)
        vphi += 2*TMath::Pi();
      fEmcEta = veta;
      fEmcPhi = vphi;
    }
  }

  fNCrossedRows = GetTPCCrossedRows();

  const AliESDVertex* vertex = 0;
  vertex = fESDEvent->GetPrimaryVertexTracks();
  if (!vertex || !vertex->GetStatus()) 
    vertex = fESDEvent->GetPrimaryVertexSPD();
  if (vertex) {
    fChi2TPCConstrainedVsGlobal            = GetChi2TPCConstrainedVsGlobal(vertex);
    fCacheChi2TPCConstrainedVsGlobalVertex = vertex;
  }
}

//_________________________________________________________________________________________________
void AliEsdTrackExt::DeleteParams()
{ 
  // Delete the unneeded params.

  delete fIp;          fIp          = 0;
  delete fTPCInner;    fTPCInner    = 0;
  delete fHMPIDp;      fHMPIDp      = 0;
  delete fFriendTrack; fFriendTrack = 0;
  if (!IsEMCAL()) {
    delete fOp;        fOp          = 0;
  }
}

//_________________________________________________________________________________________________
void AliEsdTrackExt::MakeMiniTrack(Bool_t dall, Bool_t dcon, Bool_t dtrp, Bool_t dmap, 
                                   Bool_t dits, Bool_t dtpc, Bool_t dtrd, Bool_t dtof, 
                                   Bool_t dhmp)
{ 
  // Make mini track depending on what should be reset.

  if (fCp && !dcon) 
    fCp->ResetCovariance(1);

  if (dtrp) {
    if (dcon) {
      delete fCp;
      fCp = 0;
    }
    delete fIp;          fIp          = 0;
    delete fTPCInner;    fTPCInner    = 0;
    delete fOp;          fOp          = 0;
    delete fHMPIDp;      fHMPIDp      = 0;
    delete fFriendTrack; fFriendTrack = 0;
  }

  if (dmap) {
    fTPCFitMap.Clear();
    fTPCClusterMap.Clear();
    fTPCSharedMap.Clear();
  }

  // Reset ITS track related information
  if (dits) {
    if (dall) {
      fITSchi2       = 0;
      fITSncls       = 0;       
      fITSClusterMap = 0;
      fITSSharedMap  = 0;
      fITSsignal     = 0;     
    }     
    fITSLabel      = 0;  
    for (Int_t i=0;i<4;++i) fITSdEdxSamples[i] = 0.;
    for (Int_t i=0;i<AliPID::kSPECIES;++i) fITSr[i] = 0; 
    for (Int_t i=0;i<12;++i) fITSModule[i] = -1;
  }

  // Reset TPC related track information
  if (dtpc) {
    if (dall) {
      fTPCchi2       = 0;
      fTPCchi2Iter1  = 0;
      fTPCncls       = 0;
      fTPCnclsF      = 0;
      fTPCnclsIter1  = 0;
      fTPCnclsFIter1 = 0;
      fTPCFitMap     = 0;
      fTPCClusterMap = 0;
      fTPCSharedMap  = 0;
      fTPCsignal     = 0;
      fTPCsignalS    = 0;
      fTPCsignalN    = 0;
    }
    fTPCLabel      = 0;
    fdTPC          = 0;
    fzTPC          = 0;
    fCddTPC        = 0;
    fCdzTPC        = 0;
    fCzzTPC        = 0;
    fCchi2TPC      = 0;
    for (Int_t i=0;i<AliPID::kSPECIES;++i) fTPCr[i] = 0; 
    for (Int_t i=0;i<4;++i) fTPCPoints[i] = 0;
    for (Int_t i=0; i<3;++i) fKinkIndexes[i] = 0;
    for (Int_t i=0; i<3;++i) fV0Indexes[i] = 0;
    delete fTPCdEdxInfo; fTPCdEdxInfo = 0;
  }

  // Reset TRD related track information
  if (dtrd) {
    fTRDchi2       = 0;
    fTRDncls       = 0;
    fTRDncls0      = 0;
    fTRDsignal     = 0;
    fTRDLabel      = 0;
    fTRDQuality    = 0;
    fTRDntracklets = 0;
    fTRDslices     = 0;
    fTRDBudget     = 0;
    for (Int_t i=0;i<kTRDnPlanes;++i) fTRDTimBin[i] = 0;
    for (Int_t i=0;i<AliPID::kSPECIES;++i) fTRDr[i] = 0; 
    fTRDnSlices    = 0;
    if(fTRDnSlices)
      delete[] fTRDslices;
  }

  // Reset TOF related track information
  if (dtof) {
    if (dall) {
      fTOFchi2       = 0;        
      fTOFindex      = -1;       
      fTOFsignal     = 99999;      
      fTOFCalChannel = -1;
      fTOFsignalToT  = 99999;
      fTOFsignalRaw  = 99999;
      fTOFsignalDz   = 999;
      fTOFsignalDx   = 999;
      fTOFdeltaBC    = 999;
      fTOFl0l1       = 999;
    }
    for (Int_t i=0;i<AliPID::kSPECIES;++i) fTOFr[i] = 0;
    for (Int_t i=0;i<3;++i) fTOFLabel[i] = -1;
    for (Int_t i=0;i<10;++i) fTOFInfo[i] = 0;
    for (Int_t i=0;i<AliPID::kSPECIES;++i) fTrackTime[i] = 0;
  }

  // Reset HMPID related track information
  if (dhmp) {
    fHMPIDchi2     = 0;     
    fHMPIDqn       = 0;     
    fHMPIDcluIdx   = -1;     
    fHMPIDsignal   = 0;     
    fHMPIDtrkTheta = 0;     
    fHMPIDtrkPhi   = 0;      
    fHMPIDtrkX     = 0;     
    fHMPIDtrkY     = 0;      
    fHMPIDmipX     = 0;
    fHMPIDmipY     = 0;
    for (Int_t i=0;i<AliPID::kSPECIES;++i) fHMPIDr[i] = 0;
  }
}

//_________________________________________________________________________________________________
void AliEsdTrackExt::Setup()
{
  // Setup cache with stored variables.

  fCacheNCrossedRows = fNCrossedRows;

  const AliESDVertex* vertex = 0;
  vertex = fESDEvent->GetPrimaryVertexTracks();
  if (!vertex || !vertex->GetStatus()) 
    vertex = fESDEvent->GetPrimaryVertexSPD();
  if (vertex) {
    fCacheChi2TPCConstrainedVsGlobal       = fChi2TPCConstrainedVsGlobal;
    fCacheChi2TPCConstrainedVsGlobalVertex = vertex;
  }
}
