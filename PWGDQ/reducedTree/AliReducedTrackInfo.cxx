/*
***********************************************************
  Implementation of AliReducedTrackInfo
  Contact: iarsene@cern.ch
  2015/04/08
  *********************************************************
*/

#ifndef ALIREDUCEDTRACKINFO_H
#include "AliReducedTrackInfo.h"
#endif

#include <TMath.h>
#include "AliReducedBaseTrack.h"

ClassImp(AliReducedTrackInfo)

//_______________________________________________________________________________
AliReducedTrackInfo::AliReducedTrackInfo() :
  AliReducedBaseTrack(),
  fTrackId(0),
  fStatus(0),
  fTPCPhi(0.0),
  fTPCPt(0.0),
  fTPCEta(0.0),
  fMomentumInner(0.0),
  fDCA(),
  fTrackLength(0.0),
  fMassForTracking(0.0),
  fChi2TPCConstrainedVsGlobal(0.0),
  fHelixCenter(),
  fHelixRadius(0.0),
  fITSclusterMap(0),
  fITSSharedClusterMap(0),
  fITSsignal(0.0),
  fITSnSig(),
  fITSchi2(0.0),
  fTPCNcls(0),
  fTPCCrossedRows(0),
  fTPCNclsF(0),
  fTPCNclsShared(0),
  fTPCClusterMap(0),
  fTPCsignal(0),
  fTPCsignalN(0),
  fTPCnSig(),
  fTPCchi2(0.0),
  fTPCActiveLength(0.),
  fTPCGeomLength(0.),
  fTOFbeta(0.0),
  fTOFtime(0.0),
  fTOFdx(0.0),
  fTOFdz(0.0),
  fTOFmismatchProbab(0.0),
  fTOFchi2(0.0),
  fTOFnSig(),
  fTOFdeltaBC(0),
  fTRDntracklets(),
  fTRDpid(),
  fTRDpidLQ2D(),
  fCaloClusterId(-999),
  fTrackParam(),
  fCovMatrix(),
  fMCMom(),
  fMCFreezeout(),
  fMCLabels(),
  fMCPdg(),
  fMCGeneratorIndex(-1)

{
  //
  // Constructor
  //
  fDCA[0] = 0.0; fDCA[1]=0.0;
  for(Int_t i=0; i<4; ++i) {fTPCnSig[i]=-999.; fTOFnSig[i]=-999.; fITSnSig[i]=-999.; }
  fHelixCenter[0] = 0.0; fHelixCenter[1] = 0.0;
  fTRDntracklets[0]=0; fTRDntracklets[1]=0;
  fTRDpid[0]=-999.; fTRDpid[1]=-999.;
  fTRDpidLQ2D[0] = -999.; fTRDpidLQ2D[1] = -999.;
  for(Int_t i=0;i<3;++i) {fMCMom[i]=0.; fMCFreezeout[i]=0.;}
  for(Int_t i=0;i<4;++i) {fMCLabels[i]=-9999; fMCPdg[i]=-9999;}
  for(Int_t i=0;i<6;++i) {fTrackParam[i]=0.;}
  for(Int_t i=0;i<21;++i) {fCovMatrix[i]=0.;}
}


//_______________________________________________________________________________
AliReducedTrackInfo::~AliReducedTrackInfo()
{
  //
  // De-Constructor
  //
}
