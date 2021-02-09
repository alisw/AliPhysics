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
  fStatus(0),
  fTPCPhi(0.0),
  fTPCPt(0.0),
  fTPCEta(0.0),
  fMomentumInner(0.0),
  fMomentumOnCalo(0.0),
  fPhiOnCalo(-999.),
  fEtaOnCalo(-999.),
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
  fTPCdEdxInfoQmax(),
  fTPCdEdxInfoQtot(),
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
  fTRDGTUtracklets(0),
  fTRDGTUlayermask(0),
  fTRDGTUpt(0.),
  fTRDGTUsagitta(2.0),
  fTRDGTUPID(0.),
  fMatchedEMCalClusterEnergy(-999),
  fEMCALnSigEle(-999),
  fCaloClusterId(-999),
  fTrackParam(),
  fCovMatrix(),
  fMCMom(),
  fMCFreezeout(),
  fMCLabels(),
  fMCPdg(),
  fHFProc(0),
  fMCGeneratorIndex(-1)

{
  //
  // Constructor
  //
  fDCA[0] = 0.0; fDCA[1]=0.0;
  fTPCDCA[0] = 0.0; fTPCDCA[1]=0.0;
  for(Int_t i=0; i<4; ++i) {fTPCnSig[i]=-999.; fTOFnSig[i]=-999.; fITSnSig[i]=-999.; }
  for(Int_t i=0; i<4; ++i) {fTPCdEdxInfoQmax[i]=-999.; fTPCdEdxInfoQtot[i]=-999.;}
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
AliReducedTrackInfo::AliReducedTrackInfo(const AliReducedTrackInfo &c) :
  AliReducedBaseTrack(c),
  fStatus(c.fStatus),
  fTPCPhi(c.fTPCPhi),
  fTPCPt(c.fTPCPt),
  fTPCEta(c.fTPCEta),
  fMomentumInner(c.fMomentumInner),
  fMomentumOnCalo(c.fMomentumOnCalo),
  fPhiOnCalo(c.fPhiOnCalo),
  fEtaOnCalo(c.fEtaOnCalo),
  fTrackLength(c.fTrackLength),
  fMassForTracking(c.fMassForTracking),
  fChi2TPCConstrainedVsGlobal(c.fChi2TPCConstrainedVsGlobal),
  fHelixRadius(c.fHelixRadius),
  fITSclusterMap(c.fITSclusterMap),
  fITSSharedClusterMap(c.fITSSharedClusterMap),
  fITSsignal(c.fITSsignal),
  fITSchi2(c.fITSchi2),
  fTPCNcls(c.fTPCNcls),
  fTPCCrossedRows(c.fTPCCrossedRows),
  fTPCNclsF(c.fTPCNclsF),
  fTPCNclsShared(c.fTPCNclsShared),
  fTPCClusterMap(c.fTPCClusterMap),
  fTPCsignal(c.fTPCsignal),
  fTPCsignalN(c.fTPCsignalN),
  fTPCchi2(c.fTPCchi2),
  fTPCActiveLength(c.fTPCActiveLength),
  fTPCGeomLength(c.fTPCGeomLength),
  fTOFbeta(c.fTOFbeta),
  fTOFtime(c.fTOFtime),
  fTOFdx(c.fTOFdx),
  fTOFdz(c.fTOFdz),
  fTOFmismatchProbab(c.fTOFmismatchProbab),
  fTOFchi2(c.fTOFchi2),
  fTOFdeltaBC(c.fTOFdeltaBC),
  fTRDGTUtracklets(c.fTRDGTUtracklets),
  fTRDGTUlayermask(c.fTRDGTUlayermask),
  fTRDGTUpt(c.fTRDGTUpt),
  fTRDGTUsagitta(c.fTRDGTUsagitta),
  fTRDGTUPID(c.fTRDGTUPID),
  fMatchedEMCalClusterEnergy(c.fMatchedEMCalClusterEnergy),
  fEMCALnSigEle(c.fEMCALnSigEle),
  fCaloClusterId(c.fCaloClusterId),
  fMCGeneratorIndex(c.fMCGeneratorIndex)
{
   //
   // copy constructor
   //
   fDCA[0] = c.fDCA[0]; fDCA[1]=c.fDCA[1];
   fTPCDCA[0] = c.fTPCDCA[0]; fTPCDCA[1]=c.fTPCDCA[1];
   fHelixCenter[0] = c.fHelixCenter[0]; fHelixCenter[1] = c.fHelixCenter[1];
   fTRDntracklets[0]=c.fTRDntracklets[0]; fTRDntracklets[1]=c.fTRDntracklets[1];
   fTRDpid[0]=c.fTRDpid[0]; fTRDpid[1]=c.fTRDpid[1];
   fTRDpidLQ2D[0] = c.fTRDpidLQ2D[0]; fTRDpidLQ2D[1] = c.fTRDpidLQ2D[1];
   for(Int_t i=0; i<4; ++i) {fTPCnSig[i]=c.fTPCnSig[i]; fTOFnSig[i]=c.fTOFnSig[i]; fITSnSig[i]=c.fITSnSig[i];}
   for(Int_t i=0; i<4; ++i) {fTPCdEdxInfoQmax[i]=c.fTPCdEdxInfoQmax[i]; fTPCdEdxInfoQtot[i]=c.fTPCdEdxInfoQtot[i];}
   for(Int_t i=0;i<6;++i) {fTrackParam[i]=c.fTrackParam[i];}
   for(Int_t i=0;i<21;++i) {fCovMatrix[i]=c.fCovMatrix[i];}
   for(Int_t i=0;i<3;++i) {fMCMom[i]=c.fMCMom[i]; fMCFreezeout[i]=c.fMCFreezeout[i];}
   for(Int_t i=0;i<4;++i) {fMCLabels[i]=c.fMCLabels[i]; fMCPdg[i]=c.fMCPdg[i];}
}

//_______________________________________________________________________________
AliReducedTrackInfo::~AliReducedTrackInfo()
{
  //
  // De-Constructor
  //
}
