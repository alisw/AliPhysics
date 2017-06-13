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
 **********************************************************i****************/

/* $Id$ */ 

// AliFlowTrackCuts:
// Data selection for flow framework
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)
// mods:   Redmer A. Bertens (rbertens@cern.ch)
//
// This class gurantees consistency of cut methods, trackparameter
// selection (global tracks, TPC only, etc..) and parameter mixing
// in the flow framework. Transparently handles different input types:
// ESD, MC, AOD.
// This class works in 2 steps: first the requested track parameters are
// constructed (to be set by SetParamType() ), then cuts are applied.
// the constructed track can be requested AFTER checking the cuts by
// calling GetTrack(), in this case the cut object stays in control,
// caller does not have to delete the track.
// Additionally caller can request an AliFlowTrack object to be constructed
// according the parameter mixing scenario requested by SetParamMix().
// AliFlowTrack is made using MakeFlowTrack() method, its an 'object factory'
// so caller needs to take care of the freshly created object.

#include "TFile.h"
#include <limits.h>
#include <float.h>
#include <TMatrix.h>
#include "TMCProcess.h"
#include "TParticle.h"
#include "TH2F.h"
#include "AliStack.h"
#include "TBrowser.h"
#include "TArrayD.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliVVZERO.h"
#include "AliMCParticle.h"
#include "AliESDkink.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"   // XZhang 20120604
#include "AliMultiplicity.h"
#include "AliMultSelection.h" // available from November 2015
#include "AliAODTrack.h"
#include "AliAODTracklets.h"   // XZhang 20120615
#include "AliFlowTrackSimple.h"
#include "AliFlowTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliLog.h"
#include "AliESDpid.h"
#include "AliESDPmdTrack.h"
#include "AliESDUtils.h" //TODO
#include "AliFlowBayesianPID.h"
#include "AliFlowCandidateTrack.h"
#include "AliKFParticle.h"
#include "AliESDVZERO.h"
#include "AliFlowCommonConstants.h"
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "TF2.h"


using std::cout;
using std::endl;

ClassImp(AliFlowTrackCuts)

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts():
  AliFlowTrackSimpleCuts(),
  fAliESDtrackCuts(NULL),
  fMuonTrackCuts(NULL),  // XZhang 20120604
  fQA(NULL),
  fCutMC(kFALSE),
  fCutMChasTrackReferences(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fCutMCfirstMotherPID(kFALSE),
  fMCfirstMotherPID(0),
  fIgnoreSignInMCPID(kFALSE),
  fCutMCisPrimary(kFALSE),
  fRequireTransportBitForPrimaries(kTRUE),
  fMCisPrimary(kFALSE),
  fRequireCharge(kFALSE),
  fFakesAreOK(kTRUE),
  fCutSPDtrackletDeltaPhi(kFALSE),
  fSPDtrackletDeltaPhiMax(FLT_MAX),
  fSPDtrackletDeltaPhiMin(-FLT_MAX),
  fIgnoreTPCzRange(kFALSE),
  fIgnoreTPCzRangeMax(FLT_MAX),
  fIgnoreTPCzRangeMin(-FLT_MAX),
  fCutChi2PerClusterTPC(kFALSE),
  fMaxChi2PerClusterTPC(FLT_MAX),
  fMinChi2PerClusterTPC(-FLT_MAX),
  fCutFracSharedTPCCluster(kFALSE),
  fMaxFracSharedTPCCluster(FLT_MAX),
  fCutCrossedTPCRows(kFALSE),
  fMinNCrossedRows(0),
  fMinCrossedRowsOverFindableClusters(2.),
  fCutGoldenChi2(kFALSE),
  fMaxGoldenChi2(FLT_MAX),
  fRequireTOFSignal(kFALSE),
  fCutNClustersTPC(kFALSE),
  fNClustersTPCMax(INT_MAX),
  fNClustersTPCMin(INT_MIN),  
  fCutNClustersITS(kFALSE),
  fNClustersITSMax(INT_MAX),
  fNClustersITSMin(INT_MIN),
  fCutChi2PerClusterITS(kFALSE),
  fCutITSClusterGlobal(kFALSE),
  fMaxChi2PerClusterITS(FLT_MAX),
  fUseAODFilterBit(kTRUE),
  fAODFilterBit(1),
  fCutDCAToVertexXY(kFALSE),
  fCutDCAToVertexXYPtDepAOD(kFALSE),
  fCutDCAToVertexXYAOD(kFALSE),
  fMaxDCAxyAOD(FLT_MAX),
  fCutDCAToVertexZAOD(kFALSE),
  fMaxDCAzAOD(FLT_MAX),
  fCutDCAToVertexZ(kFALSE),
  fCutMinimalTPCdedx(kFALSE),
  fMinimalTPCdedx(0.),
  fCutTPCSecbound(kFALSE),
  fCutTPCSecboundMinpt(0.2),
  fCutTPCSecboundVar(kFALSE),
  fPhiCutLow(NULL),
  fPhiCutHigh(NULL),
  fLinearizeVZEROresponse(kFALSE),
  fCentralityPercentileMin(0.),
  fCentralityPercentileMax(5.),
  fPurityLevel(0.8),
  fCutPmdDet(kFALSE),
  fPmdDet(0),
  fCutPmdAdc(kFALSE),
  fPmdAdc(0.),
  fCutPmdNcell(kFALSE),
  fPmdNcell(0.),  
  fMinKinkAngle(TMath::DegToRad()*2.),
  fMinKinkRadius(130.),
  fMaxKinkRadius(200.),
  fMinKinkQt(.05),
  fMaxKinkQt(.5),
  fMinKinkInvMassKmu(0.),
  fMaxKinkInvMassKmu(0.6),
  fForceTPCstandalone(kFALSE),
  fRequireKinkDaughters(kFALSE),
  fParamType(kGlobal),
  fParamMix(kPure),
  fKink(NULL),
  fV0(NULL),
  fTrack(NULL),
  fTrackMass(0.),
  fTrackPt(0.),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(1.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL),
  fTPCtrack(),
  fESDpid(),
  fBayesianResponse(NULL),
  fPIDsource(kTOFpid),
  fTPCpidCuts(NULL),
  fTOFpidCuts(NULL),
  fParticleID(AliPID::kUnknown),
  fParticleProbability(.9),
  fAllowTOFmismatchFlag(kFALSE),
  fRequireStrictTOFTPCagreement(kFALSE),
  fCutRejectElectronsWithTPCpid(kFALSE),
  fProbBayes(0.0),
  fCurrCentr(0.0),
  fPtTOFPIDoff(0.),
  fVZEROgainEqualization(NULL),
  fVZEROgainEqualizationCen(NULL),
  fApplyRecentering(kFALSE),
  fVZEROgainEqualizationPerRing(kFALSE),
  fDivSigma(kTRUE),
  fChi2A(0x0),
  fChi2C(0x0),
  fChi3A(0x0),
  fChi3C(0x0),
  fPIDResponse(NULL),
  fNsigmaCut2(9),
  fPurityFunctionsFile(0),
  fPurityFunctionsList(0),
  fCutITSclusterShared(kFALSE),
  fMaxITSclusterShared(0),
  fCutITSChi2(kFALSE),
  fMaxITSChi2(0),
  fRun(0)
{
  //io constructor 
  SetPriors(); //init arrays
  // New PID procedure (Bayesian Combined PID)
  // allocating here is necessary because we don't 
  // stream this member
  // TODO: fix streaming problems AliFlowBayesianPID
  fBayesianResponse = new AliFlowBayesianPID();
  fBayesianResponse->SetNewTrackParam();
  for(Int_t i(0); i < 4; i++) {
      fVZEROApol[i] = 0;
      fVZEROCpol[i] = 0;
  }
  for(Int_t i(0); i < 8; i++) fUseVZERORing[i] = kTRUE;
    
  for(Int_t i(0) ; i < 180; i++) {
     fPurityFunction[i]=NULL;
  }

  fPhiCutLow = new TF1("fPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);
  fPhiCutHigh = new TF1("fPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts(const char* name):
  AliFlowTrackSimpleCuts(name),
  fAliESDtrackCuts(NULL),
  fMuonTrackCuts(NULL),  // XZhang 20120604
  fQA(NULL),
  fCutMC(kFALSE),
  fCutMChasTrackReferences(kFALSE),
  fCutMCprocessType(kFALSE),
  fMCprocessType(kPNoProcess),
  fCutMCPID(kFALSE),
  fMCPID(0),
  fCutMCfirstMotherPID(kFALSE),
  fMCfirstMotherPID(0),
  fIgnoreSignInMCPID(kFALSE),
  fCutMCisPrimary(kFALSE),
  fRequireTransportBitForPrimaries(kTRUE),
  fMCisPrimary(kFALSE),
  fRequireCharge(kFALSE),
  fFakesAreOK(kTRUE),
  fCutSPDtrackletDeltaPhi(kFALSE),
  fSPDtrackletDeltaPhiMax(FLT_MAX),
  fSPDtrackletDeltaPhiMin(-FLT_MAX),
  fIgnoreTPCzRange(kFALSE),
  fIgnoreTPCzRangeMax(FLT_MAX),
  fIgnoreTPCzRangeMin(-FLT_MAX),
  fCutChi2PerClusterTPC(kFALSE),
  fMaxChi2PerClusterTPC(FLT_MAX),
  fMinChi2PerClusterTPC(-FLT_MAX),
  fCutFracSharedTPCCluster(kFALSE),
  fMaxFracSharedTPCCluster(FLT_MAX),
  fCutCrossedTPCRows(kFALSE),
  fMinNCrossedRows(0),
  fMinCrossedRowsOverFindableClusters(2.),
  fCutGoldenChi2(kFALSE),
  fMaxGoldenChi2(FLT_MAX),
  fRequireTOFSignal(kFALSE),
  fCutNClustersTPC(kFALSE),
  fNClustersTPCMax(INT_MAX),
  fNClustersTPCMin(INT_MIN),  
  fCutNClustersITS(kFALSE),
  fNClustersITSMax(INT_MAX),
  fNClustersITSMin(INT_MIN),
  fCutChi2PerClusterITS(kFALSE),
  fCutITSClusterGlobal(kFALSE),
  fMaxChi2PerClusterITS(FLT_MAX),
  fUseAODFilterBit(kTRUE),
  fAODFilterBit(1),
  fCutDCAToVertexXY(kFALSE),
  fCutDCAToVertexXYPtDepAOD(kFALSE),
  fCutDCAToVertexXYAOD(kFALSE),
  fMaxDCAxyAOD(FLT_MAX),
  fCutDCAToVertexZAOD(kFALSE),
  fMaxDCAzAOD(FLT_MAX),
  fCutDCAToVertexZ(kFALSE),
  fCutMinimalTPCdedx(kFALSE),
  fMinimalTPCdedx(0.),
  fCutTPCSecbound(kFALSE),
  fCutTPCSecboundMinpt(0.2),
  fCutTPCSecboundVar(kFALSE),
  fPhiCutLow(NULL),
  fPhiCutHigh(NULL),
  fLinearizeVZEROresponse(kFALSE),
  fCentralityPercentileMin(0.),
  fCentralityPercentileMax(5.),
  fPurityLevel(0.8),
  fCutPmdDet(kFALSE),
  fPmdDet(0),
  fCutPmdAdc(kFALSE),
  fPmdAdc(0.),
  fCutPmdNcell(kFALSE),
  fPmdNcell(0.),  
  fMinKinkAngle(TMath::DegToRad()*2.),
  fMinKinkRadius(130.),
  fMaxKinkRadius(200.),
  fMinKinkQt(0.05),
  fMaxKinkQt(0.5),
  fMinKinkInvMassKmu(0.0),
  fMaxKinkInvMassKmu(0.6),
  fForceTPCstandalone(kFALSE),
  fRequireKinkDaughters(kFALSE),
  fParamType(kGlobal),
  fParamMix(kPure),
  fKink(NULL),
  fV0(NULL),
  fTrack(NULL),
  fTrackMass(0.),
  fTrackPt(0.),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(1.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL),
  fTPCtrack(),
  fESDpid(),
  fBayesianResponse(NULL),
  fPIDsource(kTOFpid),
  fTPCpidCuts(NULL),
  fTOFpidCuts(NULL),
  fParticleID(AliPID::kUnknown),
  fParticleProbability(.9),
  fAllowTOFmismatchFlag(kFALSE),
  fRequireStrictTOFTPCagreement(kFALSE),
  fCutRejectElectronsWithTPCpid(kFALSE),
  fProbBayes(0.0),
  fCurrCentr(0.0),
  fPtTOFPIDoff(0.),
  fVZEROgainEqualization(NULL),
  fVZEROgainEqualizationCen(NULL),
  fApplyRecentering(kFALSE),
  fVZEROgainEqualizationPerRing(kFALSE),
  fDivSigma(kTRUE),
  fChi2A(0x0),
  fChi2C(0x0),
  fChi3A(0x0),
  fChi3C(0x0),
  fPIDResponse(NULL),
  fNsigmaCut2(9),
  fPurityFunctionsFile(0),
  fPurityFunctionsList(0),
  fCutITSclusterShared(kFALSE),
  fMaxITSclusterShared(0),
  fCutITSChi2(kFALSE),
  fMaxITSChi2(0),
  fRun(0)
{
  //constructor
  SetTitle("AliFlowTrackCuts");
  fESDpid.GetTPCResponse().SetBetheBlochParameters( 0.0283086,
                                                    2.63394e+01,
                                                    5.04114e-11,
                                                    2.12543e+00,
                                                    4.88663e+00 );
  SetPriors(); //init arrays
  // New PID procedure (Bayesian Combined PID)
  fBayesianResponse = new AliFlowBayesianPID();
  fBayesianResponse->SetNewTrackParam();
  for(Int_t i(0); i < 4; i++) {
      fVZEROApol[i] = 0;
      fVZEROCpol[i] = 0;
  }
  for(Int_t i(0); i < 8; i++) fUseVZERORing[i] = kTRUE;
    
  for(int i=0;i<180;i++){
     fPurityFunction[i]=NULL;
  }

  fPhiCutLow = new TF1("fPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);
  fPhiCutHigh = new TF1("fPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::AliFlowTrackCuts(const AliFlowTrackCuts& that):
  AliFlowTrackSimpleCuts(that),
  fAliESDtrackCuts(NULL),
  fMuonTrackCuts(NULL),  // XZhang 20120604
  fQA(NULL),
  fCutMC(that.fCutMC),
  fCutMChasTrackReferences(that.fCutMChasTrackReferences),
  fCutMCprocessType(that.fCutMCprocessType),
  fMCprocessType(that.fMCprocessType),
  fCutMCPID(that.fCutMCPID),
  fMCPID(that.fMCPID),
  fCutMCfirstMotherPID(that.fCutMCfirstMotherPID),
  fMCfirstMotherPID(that.fMCfirstMotherPID),
  fIgnoreSignInMCPID(that.fIgnoreSignInMCPID),
  fCutMCisPrimary(that.fCutMCisPrimary),
  fRequireTransportBitForPrimaries(that.fRequireTransportBitForPrimaries),
  fMCisPrimary(that.fMCisPrimary),
  fRequireCharge(that.fRequireCharge),
  fFakesAreOK(that.fFakesAreOK),
  fCutSPDtrackletDeltaPhi(that.fCutSPDtrackletDeltaPhi),
  fSPDtrackletDeltaPhiMax(that.fSPDtrackletDeltaPhiMax),
  fSPDtrackletDeltaPhiMin(that.fSPDtrackletDeltaPhiMin),
  fIgnoreTPCzRange(that.fIgnoreTPCzRange),
  fIgnoreTPCzRangeMax(that.fIgnoreTPCzRangeMax),
  fIgnoreTPCzRangeMin(that.fIgnoreTPCzRangeMin),
  fCutChi2PerClusterTPC(that.fCutChi2PerClusterTPC),
  fMaxChi2PerClusterTPC(that.fMaxChi2PerClusterTPC),
  fMinChi2PerClusterTPC(that.fMinChi2PerClusterTPC),
  fCutFracSharedTPCCluster(that.fCutFracSharedTPCCluster),
  fMaxFracSharedTPCCluster(that.fMaxFracSharedTPCCluster),
  fCutCrossedTPCRows(that.fCutCrossedTPCRows),
  fMinNCrossedRows(that.fMinNCrossedRows),
  fMinCrossedRowsOverFindableClusters(that.fMinCrossedRowsOverFindableClusters),
  fCutGoldenChi2(that.fCutGoldenChi2),
  fMaxGoldenChi2(that.fMaxGoldenChi2),
  fRequireTOFSignal(that.fRequireTOFSignal),
  fCutNClustersTPC(that.fCutNClustersTPC),
  fNClustersTPCMax(that.fNClustersTPCMax),
  fNClustersTPCMin(that.fNClustersTPCMin),
  fCutNClustersITS(that.fCutNClustersITS),
  fNClustersITSMax(that.fNClustersITSMax),
  fNClustersITSMin(that.fNClustersITSMin),
  fCutChi2PerClusterITS(that.fCutChi2PerClusterITS),
  fCutITSClusterGlobal(that.fCutITSClusterGlobal),
  fMaxChi2PerClusterITS(that.fMaxChi2PerClusterITS),
  fUseAODFilterBit(that.fUseAODFilterBit),
  fAODFilterBit(that.fAODFilterBit),
  fCutDCAToVertexXY(that.fCutDCAToVertexXY),
  fCutDCAToVertexXYPtDepAOD(that.fCutDCAToVertexXYPtDepAOD),
  fCutDCAToVertexXYAOD(that.fCutDCAToVertexXYAOD),
  fMaxDCAxyAOD(that.fMaxDCAxyAOD),
  fCutDCAToVertexZAOD(that.fCutDCAToVertexZAOD),
  fMaxDCAzAOD(that.fMaxDCAzAOD),
  fCutDCAToVertexZ(that.fCutDCAToVertexZ),
  fCutMinimalTPCdedx(that.fCutMinimalTPCdedx),
  fMinimalTPCdedx(that.fMinimalTPCdedx),
  fCutTPCSecbound(that.fCutTPCSecbound),
  fCutTPCSecboundMinpt(that.fCutTPCSecboundMinpt),
  fCutTPCSecboundVar(that.fCutTPCSecboundVar),
  fPhiCutLow(NULL),
  fPhiCutHigh(NULL),
  fLinearizeVZEROresponse(that.fLinearizeVZEROresponse),
  fCentralityPercentileMin(that.fCentralityPercentileMin),
  fCentralityPercentileMax(that.fCentralityPercentileMax),
  fPurityLevel(that.fPurityLevel),
  fCutPmdDet(that.fCutPmdDet),
  fPmdDet(that.fPmdDet),
  fCutPmdAdc(that.fCutPmdAdc),
  fPmdAdc(that.fPmdAdc),
  fCutPmdNcell(that.fCutPmdNcell),
  fPmdNcell(that.fPmdNcell),  
  fMinKinkAngle(that.fMinKinkAngle),
  fMinKinkRadius(that.fMinKinkRadius),
  fMaxKinkRadius(that.fMaxKinkRadius),
  fMinKinkQt(that.fMinKinkQt),
  fMaxKinkQt(that.fMaxKinkQt),
  fMinKinkInvMassKmu(that.fMinKinkInvMassKmu),
  fMaxKinkInvMassKmu(that.fMaxKinkInvMassKmu),
  fForceTPCstandalone(that.fForceTPCstandalone),
  fRequireKinkDaughters(that.fRequireKinkDaughters),
  fParamType(that.fParamType),
  fParamMix(that.fParamMix),
  fKink(NULL),
  fV0(NULL),
  fTrack(NULL),
  fTrackMass(0.),
  fTrackPt(0.),
  fTrackPhi(0.),
  fTrackEta(0.),
  fTrackWeight(1.),
  fTrackLabel(INT_MIN),
  fMCevent(NULL),
  fMCparticle(NULL),
  fEvent(NULL),
  fTPCtrack(),
  fESDpid(that.fESDpid),
  fBayesianResponse(NULL),
  fPIDsource(that.fPIDsource),
  fTPCpidCuts(NULL),
  fTOFpidCuts(NULL),
  fParticleID(that.fParticleID),
  fParticleProbability(that.fParticleProbability),
  fAllowTOFmismatchFlag(that.fAllowTOFmismatchFlag),
  fRequireStrictTOFTPCagreement(that.fRequireStrictTOFTPCagreement),
  fCutRejectElectronsWithTPCpid(that.fCutRejectElectronsWithTPCpid),
  fProbBayes(0.0),
  fCurrCentr(0.0),
  fPtTOFPIDoff(that.fPtTOFPIDoff),
  fVZEROgainEqualization(NULL),
  fVZEROgainEqualizationCen(NULL),
  fApplyRecentering(that.fApplyRecentering),
  fVZEROgainEqualizationPerRing(that.fVZEROgainEqualizationPerRing),
  fDivSigma(that.fDivSigma),
  fChi2A(0x0),
  fChi2C(0x0),
  fChi3A(0x0),
  fChi3C(0x0),
  fPIDResponse(that.fPIDResponse),
  fNsigmaCut2(that.fNsigmaCut2),
  fPurityFunctionsFile(that.fPurityFunctionsFile),
  fPurityFunctionsList(that.fPurityFunctionsList),
  fCutITSclusterShared(kFALSE),
  fMaxITSclusterShared(0),
  fCutITSChi2(kFALSE),
  fMaxITSChi2(0),
  fRun(0)
{
  //copy constructor
  if (that.fTPCpidCuts) fTPCpidCuts = new TMatrixF(*(that.fTPCpidCuts));
  if (that.fTOFpidCuts) fTOFpidCuts = new TMatrixF(*(that.fTOFpidCuts));
  if (that.fAliESDtrackCuts) fAliESDtrackCuts = new AliESDtrackCuts(*(that.fAliESDtrackCuts));
  if (that.fMuonTrackCuts)   fMuonTrackCuts   = new AliMuonTrackCuts(*(that.fMuonTrackCuts));  // XZhang 20120604
  SetPriors(); //init arrays
  if (that.fQA) DefineHistograms();
  // would be neater via copy ctor of TArrayD
  fChi2A = new TArrayD(that.fChi2A->GetSize(), that.fChi2A->GetArray());
  fChi2C = new TArrayD(that.fChi2C->GetSize(), that.fChi2C->GetArray());
  fChi3A = new TArrayD(that.fChi3A->GetSize(), that.fChi3A->GetArray());
  fChi3C = new TArrayD(that.fChi3C->GetSize(), that.fChi3C->GetArray());
  // New PID procedure (Bayesian Combined PID)
  fBayesianResponse = new AliFlowBayesianPID();
  fBayesianResponse->SetNewTrackParam();
    
  // VZERO gain calibration
  // no reason to init fVZEROgainEqualizationPerRing, will be initialized on node if necessary
  // pointer is set to NULL in initialization list of this constructor
//  if (that.fVZEROgainEqualization) fVZEROgainEqualization = new TH1(*(that.fVZEROgainEqualization));
  for(Int_t i(0); i < 4; i++) { // no use to copy these guys since they're only initialized on worked node
      fVZEROApol[i] = 0.;
      fVZEROCpol[i] = 0.;
  }
  for(Int_t i(0); i < 8; i++) fUseVZERORing[i] = that.fUseVZERORing[i];
  
  for(Int_t i(0); i < 180; i++) {
     fPurityFunction[i]=that.fPurityFunction[i];
  }
  
  fPhiCutLow = new TF1("fPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);
  fPhiCutHigh = new TF1("fPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);
}

//-----------------------------------------------------------------------
AliFlowTrackCuts& AliFlowTrackCuts::operator=(const AliFlowTrackCuts& that)
{
  //assignment
  if (this==&that) return *this;

  AliFlowTrackSimpleCuts::operator=(that);
  //the following may seem excessive but if AliESDtrackCuts properly does copy and clone
  //this approach is better memory-fragmentation-wise in some cases
  if (that.fAliESDtrackCuts && fAliESDtrackCuts) *fAliESDtrackCuts=*(that.fAliESDtrackCuts);
  if (that.fAliESDtrackCuts && !fAliESDtrackCuts) fAliESDtrackCuts=new AliESDtrackCuts(*(that.fAliESDtrackCuts));
  if (!that.fAliESDtrackCuts) delete fAliESDtrackCuts; fAliESDtrackCuts=NULL;

  if ( that.fMuonTrackCuts &&  fMuonTrackCuts) *fMuonTrackCuts = *(that.fMuonTrackCuts);                   // XZhang 20120604
  if ( that.fMuonTrackCuts && !fMuonTrackCuts)  fMuonTrackCuts = new AliMuonTrackCuts(*(that.fMuonTrackCuts));  // XZhang 20120604
  if (!that.fMuonTrackCuts) delete fMuonTrackCuts; fMuonTrackCuts = NULL;                                  // XZhang 20120604
  if (!that.fVZEROgainEqualization) delete fVZEROgainEqualization; fVZEROgainEqualization = NULL;
  //these guys we don't need to copy, just reinit
  if (that.fQA) {fQA->Delete(); delete fQA; fQA=NULL; DefineHistograms();} 
  fCutMC=that.fCutMC;
  fCutMChasTrackReferences=that.fCutMChasTrackReferences;
  fCutMCprocessType=that.fCutMCprocessType;
  fMCprocessType=that.fMCprocessType;
  fCutMCPID=that.fCutMCPID;
  fMCPID=that.fMCPID;
  fCutMCfirstMotherPID=that.fCutMCfirstMotherPID;
  fMCfirstMotherPID=that.fMCfirstMotherPID;
  fIgnoreSignInMCPID=that.fIgnoreSignInMCPID,
  fCutMCisPrimary=that.fCutMCisPrimary;
  fRequireTransportBitForPrimaries=that.fRequireTransportBitForPrimaries;
  fMCisPrimary=that.fMCisPrimary;
  fRequireCharge=that.fRequireCharge;
  fFakesAreOK=that.fFakesAreOK;
  fCutSPDtrackletDeltaPhi=that.fCutSPDtrackletDeltaPhi;
  fSPDtrackletDeltaPhiMax=that.fSPDtrackletDeltaPhiMax;
  fSPDtrackletDeltaPhiMin=that.fSPDtrackletDeltaPhiMin;
  fIgnoreTPCzRange=that.fIgnoreTPCzRange;
  fIgnoreTPCzRangeMax=that.fIgnoreTPCzRangeMax;
  fIgnoreTPCzRangeMin=that.fIgnoreTPCzRangeMin;
  fCutChi2PerClusterTPC=that.fCutChi2PerClusterTPC;
  fMaxChi2PerClusterTPC=that.fMaxChi2PerClusterTPC;
  fMinChi2PerClusterTPC=that.fMinChi2PerClusterTPC;
  fCutFracSharedTPCCluster=that.fCutFracSharedTPCCluster;
  fMaxFracSharedTPCCluster=that.fMaxFracSharedTPCCluster;
  fCutCrossedTPCRows=that.fCutCrossedTPCRows;
  fMinNCrossedRows=that.fMinNCrossedRows;
  fMinCrossedRowsOverFindableClusters=that.fMinCrossedRowsOverFindableClusters;
  fCutGoldenChi2=that.fCutGoldenChi2;
  fMaxGoldenChi2=that.fMaxGoldenChi2;
  fRequireTOFSignal=that.fRequireTOFSignal;
  fCutNClustersTPC=that.fCutNClustersTPC;
  fNClustersTPCMax=that.fNClustersTPCMax;
  fNClustersTPCMin=that.fNClustersTPCMin;  
  fCutNClustersITS=that.fCutNClustersITS;
  fNClustersITSMax=that.fNClustersITSMax;
  fNClustersITSMin=that.fNClustersITSMin;
  fCutChi2PerClusterITS=that.fCutChi2PerClusterITS;
  fCutITSClusterGlobal=that.fCutITSClusterGlobal;
  fMaxChi2PerClusterITS=that.fMaxChi2PerClusterITS;
  fUseAODFilterBit=that.fUseAODFilterBit;
  fAODFilterBit=that.fAODFilterBit;
  fCutDCAToVertexXY=that.fCutDCAToVertexXY;
  fCutDCAToVertexXYPtDepAOD=that.fCutDCAToVertexXYPtDepAOD;
  fCutDCAToVertexXYAOD=that.fCutDCAToVertexXYAOD;
  fMaxDCAxyAOD=that.fMaxDCAxyAOD;
  fCutDCAToVertexZAOD=that.fCutDCAToVertexZAOD;
  fMaxDCAzAOD=that.fMaxDCAzAOD;
  fCutDCAToVertexZ=that.fCutDCAToVertexZ;
  fCutMinimalTPCdedx=that.fCutMinimalTPCdedx;
  fMinimalTPCdedx=that.fMinimalTPCdedx;
  fCutTPCSecbound=that.fCutTPCSecbound;
  fCutTPCSecboundMinpt=that.fCutTPCSecboundMinpt;
  fCutTPCSecboundVar=that.fCutTPCSecboundVar;
  fPhiCutLow=NULL;
  fPhiCutHigh=NULL;
  fPtTOFPIDoff=that.fPtTOFPIDoff;
  fLinearizeVZEROresponse=that.fLinearizeVZEROresponse;
  fCentralityPercentileMin=that.fCentralityPercentileMin;
  fCentralityPercentileMax=that.fCentralityPercentileMax;
  fPurityLevel=that.fPurityLevel;
  fCutPmdDet=that.fCutPmdDet;
  fPmdDet=that.fPmdDet;
  fCutPmdAdc=that.fCutPmdAdc;
  fPmdAdc=that.fPmdAdc;
  fCutPmdNcell=that.fCutPmdNcell;
  fPmdNcell=that.fPmdNcell;
  fMinKinkAngle=that.fMinKinkAngle;
  fMinKinkRadius=that.fMinKinkRadius;
  fMaxKinkRadius=that.fMaxKinkRadius;
  fMinKinkQt=that.fMinKinkQt;
  fMaxKinkQt=that.fMaxKinkQt;
  fMinKinkInvMassKmu=that.fMinKinkInvMassKmu;
  fMaxKinkInvMassKmu=that.fMaxKinkInvMassKmu;
  fForceTPCstandalone=that.fForceTPCstandalone;
  fRequireKinkDaughters=that.fRequireKinkDaughters;
  
  fParamType=that.fParamType;
  fParamMix=that.fParamMix;

  fKink=NULL;
  fV0=NULL;
  fTrack=NULL;
  fTrackMass=0.;
  fTrackPt=0.;
  fTrackEta=0.;
  fTrackPhi=0.;
  fTrackWeight=1.;
  fTrackLabel=INT_MIN;
  fMCevent=NULL;
  fMCparticle=NULL;
  fEvent=NULL;

  fESDpid = that.fESDpid;
  fPIDsource = that.fPIDsource;

  delete fTPCpidCuts;
  delete fTOFpidCuts;
  if (that.fTPCpidCuts) fTPCpidCuts = new TMatrixF(*(that.fTPCpidCuts));
  if (that.fTOFpidCuts) fTOFpidCuts = new TMatrixF(*(that.fTOFpidCuts));

  fParticleID=that.fParticleID;
  fParticleProbability=that.fParticleProbability;
  fAllowTOFmismatchFlag=that.fAllowTOFmismatchFlag;
  fRequireStrictTOFTPCagreement=that.fRequireStrictTOFTPCagreement;
  fCutRejectElectronsWithTPCpid=that.fCutRejectElectronsWithTPCpid;
  fProbBayes = that.fProbBayes;
  fCurrCentr = that.fCurrCentr;
    
  fApplyRecentering = that.fApplyRecentering;
  fVZEROgainEqualizationPerRing = that.fVZEROgainEqualizationPerRing;
  fDivSigma = that.fDivSigma;
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
  if (that.fVZEROgainEqualization) fVZEROgainEqualization = new TH1(*(that.fVZEROgainEqualization));
  if (that.fVZEROgainEqualizationCen) fVZEROgainEqualizationCen = new TH2(*(that.fVZEROgainEqualizationCen));
#else
  //PH Lets try Clone, however the result might be wrong
  if (that.fVZEROgainEqualization) fVZEROgainEqualization = (TH1*)that.fVZEROgainEqualization->Clone();
  if (that.fVZEROgainEqualizationCen) fVZEROgainEqualizationCen = (TH2*)that.fVZEROgainEqualizationCen->Clone();
#endif
    
  for(Int_t i(0); i < 4; i++) { // no use to copy these guys since they're only initialized on worked node
      fVZEROApol[i] = that.fVZEROApol[i];
      fVZEROCpol[i] = that.fVZEROCpol[i];
  }
  for(Int_t i(0); i < 8; i++) fUseVZERORing[i] = that.fUseVZERORing[i];
  // would be neater via copy ctor of TArrayD
  fChi2A = new TArrayD(that.fChi2A->GetSize(), that.fChi2A->GetArray());
  fChi2C = new TArrayD(that.fChi2C->GetSize(), that.fChi2C->GetArray());
  fChi3A = new TArrayD(that.fChi3A->GetSize(), that.fChi3A->GetArray());
  fChi3C = new TArrayD(that.fChi3C->GetSize(), that.fChi3C->GetArray());

  fPIDResponse = that.fPIDResponse;
  fNsigmaCut2 = that.fNsigmaCut2;
 
  fRun = that.fRun;
  
  fPhiCutLow = new TF1("fPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);
  fPhiCutHigh = new TF1("fPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);
  
  return *this;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts::~AliFlowTrackCuts()
{
  //dtor
  delete fAliESDtrackCuts;
  delete fTPCpidCuts;
  delete fTOFpidCuts;
  if (fMuonTrackCuts) delete fMuonTrackCuts;  // XZhang 20120604
  if (fQA) { fQA->SetOwner(); fQA->Delete(); delete fQA; }
  if (fVZEROgainEqualization) {
      delete fVZEROgainEqualization;
      fVZEROgainEqualization = 0x0;
  }
  if (fVZEROgainEqualizationCen) {
      delete fVZEROgainEqualizationCen;
      fVZEROgainEqualizationCen = 0x0;
  }
  if(fChi2A) {
      delete fChi2A;
      fChi2A = 0x0;
  }
  if(fChi3A) {
      delete fChi3A;
      fChi3A = 0x0;
  }
  if(fChi2C) {
      delete fChi2C;
      fChi2C = 0x0;
  }
  if(fChi3C) {
      delete fChi3C;
      fChi3C = 0x0;
  }
  if(fPurityFunctionsFile) {
    fPurityFunctionsFile->Close();
    fPurityFunctionsFile = 0x0;
  }
  for(int i=0;i<180;i++){
    if(fPurityFunction[i]) {
        delete fPurityFunction[i];
        fPurityFunction[i] = 0;
    }
  }
  if(fPhiCutLow) delete fPhiCutLow;
  if(fPhiCutHigh) delete fPhiCutHigh;
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::SetEvent(AliVEvent* event, AliMCEvent* mcEvent)
{
  //set the event
  Clear();
  fEvent=event;
  fMCevent=mcEvent;
 
 //set run number
 if(fEvent->GetRunNumber()) fRun = fEvent->GetRunNumber();

  // Get PID response
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man){
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if(inputHandler) fPIDResponse=inputHandler->GetPIDResponse();
  }

  //do the magic for ESD
  AliESDEvent* myESD = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent* myAOD = dynamic_cast<AliAODEvent*>(event);
  if (fCutPID && myESD)
  {
    //TODO: maybe call it only for the TOF options?
    // Added by F. Noferini for TOF PID
    // old procedure now implemented inside fBayesianResponse
    //    fESDpid.MakePID(myESD,kFALSE);
    // new procedure
    fBayesianResponse->SetDetResponse(myESD, fCurrCentr,AliESDpid::kTOF_T0,kTRUE); // centrality = PbPb centrality class (0-100%) or -1 for pp collisions
    fESDpid.SetTOFResponse(myESD,AliESDpid::kTOF_T0);
    // End F. Noferini added part
  }
  if (fCutPID && myAOD){
    fBayesianResponse->SetDetResponse(myAOD, fCurrCentr,AliESDpid::kTOF_T0); // centrality = PbPb centrality class (0-100%) or -1 for pp collisions
    if(myAOD->GetTOFHeader()){
      fESDpid.SetTOFResponse(myAOD,AliESDpid::kTOF_T0);
    }   
    else{ // corrected on the fly track by track if tof header is not present (old AODs)
    }
    // End F. Noferini added part
  }
  
  if(fPIDsource==kTOFbayesian) fBayesianResponse->SetDetAND(1);
  else if(fPIDsource==kTPCbayesian) fBayesianResponse->ResetDetOR(1);
    
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::SetCutMC( Bool_t b )
{
  //will we be cutting on MC information?
  fCutMC=b;

  //if we cut on MC info then also the Bethe Bloch should be the one tuned for MC
  if (fCutMC)
  {
    fESDpid.GetTPCResponse().SetBetheBlochParameters( 2.15898e+00/50.,
                                                       1.75295e+01,
                                                       3.40030e-09,
                                                       1.96178e+00,
                                                       3.91720e+00);
  }
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsSelected(TObject* obj, Int_t id)
{
  //check cuts
  AliVParticle* vparticle = dynamic_cast<AliVParticle*>(obj);
//if (vparticle) return PassesCuts(vparticle);                // XZhang 20120604
  if (vparticle) {                                            // XZhang 20120604
    if (fParamType==kMUON) return PassesMuonCuts(vparticle);  // XZhang 20120604
    return PassesCuts(vparticle);                             // XZhang 20120604
  }                                                           // XZhang 20120604

  AliFlowTrackSimple* flowtrack = dynamic_cast<AliFlowTrackSimple*>(obj);
  if (flowtrack) return PassesCuts(flowtrack);
  AliMultiplicity* tracklets = dynamic_cast<AliMultiplicity*>(obj);
  if (tracklets) return PassesCuts(tracklets,id);
  AliAODTracklets* trkletAOD = dynamic_cast<AliAODTracklets*>(obj);  // XZhang 20120615
  if (trkletAOD) return PassesCuts(trkletAOD,id);                    // XZhang 20120615
  AliESDPmdTrack* pmdtrack = dynamic_cast<AliESDPmdTrack*>(obj);
  if (pmdtrack) return PassesPMDcuts(pmdtrack);
  AliVVZERO* vvzero = dynamic_cast<AliVVZERO*>(obj);    // downcast to base class
  if (vvzero) return PassesVZEROcuts(id);
  AliESDkink* kink = dynamic_cast<AliESDkink*>(obj);
  if (kink) return PassesCuts(kink);
  //AliESDv0* v0 = dynamic_cast<AliESDv0*>(obj);
  //if (v0) return PassesCuts(v0);
 
  return kFALSE;  //default when passed wrong type of object
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsSelectedMCtruth(TObject* obj, Int_t id)
{
  //check cuts
  AliVParticle* vparticle = dynamic_cast<AliVParticle*>(obj);
  if (vparticle) 
  {
    return PassesMCcuts(fMCevent,(fFakesAreOK)?TMath::Abs(vparticle->GetLabel()):vparticle->GetLabel());
  }
  AliMultiplicity* tracklets = dynamic_cast<AliMultiplicity*>(obj);
  if (tracklets)
  {
    Int_t label0 = tracklets->GetLabel(id,0);
    Int_t label1 = tracklets->GetLabel(id,1);
    Int_t label = (fFakesAreOK)?TMath::Abs(label0):label0;
    if (label0!=label1) label = -666;
    return PassesMCcuts(fMCevent,label);
  }
  return kFALSE;  //default when passed wrong type of object
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(const AliFlowTrackSimple* track)
{
  //check cuts on a flowtracksimple

  //clean up from last iteration
  ClearTrack();
  return AliFlowTrackSimpleCuts::PassesCuts(track);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(const AliMultiplicity* tracklet, Int_t id)
{
  //check cuts on a tracklets

  if (id<0) return kFALSE;

  //clean up from last iteration, and init label
  ClearTrack();

  fTrackPhi = tracklet->GetPhi(id);
  fTrackEta = tracklet->GetEta(id);
  fTrackWeight = 1.0;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) return kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) return kFALSE;}

  //check MC info if available
  //if the 2 clusters have different label track cannot be good
  //and should therefore not pass the mc cuts
  Int_t label0 = tracklet->GetLabel(id,0);
  Int_t label1 = tracklet->GetLabel(id,1);
  //if possible get label and mcparticle
  fTrackLabel = (label0==label1)?tracklet->GetLabel(id,1):-1;
  if (!fFakesAreOK && fTrackLabel<0) return kFALSE;
  if (fTrackLabel>=0 && fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  //check MC cuts
  if (fCutMC && !PassesMCcuts()) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(const AliAODTracklets* tracklet, Int_t id)
{
  // XZhang 20120615
  //check cuts on a tracklets

  if (id<0) return kFALSE;

  //clean up from last iteration, and init label
  ClearTrack();

  fTrackPhi = tracklet->GetPhi(id);
//fTrackEta = tracklet->GetEta(id);
  fTrackEta = -1.*TMath::Log(TMath::Tan(tracklet->GetTheta(id)/2.));
  fTrackWeight = 1.0;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) return kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) return kFALSE;}

  //check MC info if available
  //if the 2 clusters have different label track cannot be good
  //and should therefore not pass the mc cuts
  Int_t label0 = tracklet->GetLabel(id,0);
  Int_t label1 = tracklet->GetLabel(id,1);
  //if possible get label and mcparticle
  fTrackLabel = (label0==label1)?tracklet->GetLabel(id,1):-1;
  if (!fFakesAreOK && fTrackLabel<0) return kFALSE;
  if (fTrackLabel>=0 && fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  //check MC cuts
  if (fCutMC && !PassesMCcuts()) return kFALSE;
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMCcuts(AliMCEvent* mcEvent, Int_t label)
{
  //check the MC info
  if (!mcEvent) return kFALSE;
  if (label<0) return kFALSE;//otherwise AliCMevent prints a warning before returning NULL
  AliMCParticle* mcparticle = static_cast<AliMCParticle*>(mcEvent->GetTrack(label));
  if (!mcparticle) {AliError("no MC track"); return kFALSE;}

  if (fCutMCisPrimary)
  {
    if (IsPhysicalPrimary(mcEvent,label,fRequireTransportBitForPrimaries) != fMCisPrimary) return kFALSE;
  }
  if (fCutMCPID)
  {
    Int_t pdgCode = mcparticle->PdgCode();
    if (fIgnoreSignInMCPID) 
    {
      if (TMath::Abs(fMCPID) != TMath::Abs(pdgCode)) return kFALSE;
    }
    else 
    {
      if (fMCPID != pdgCode) return kFALSE;
    }
  }
  if (fCutMCfirstMotherPID)
  {

    TParticle* tparticle=mcparticle->Particle();
    Int_t firstMotherLabel = 0;
    if (tparticle) { firstMotherLabel = tparticle->GetFirstMother(); }
    AliVParticle* firstMotherParticle = mcEvent->GetTrack(firstMotherLabel);
    Int_t pdgcodeFirstMother = 0;
    if (firstMotherParticle) { pdgcodeFirstMother = firstMotherParticle->PdgCode(); }
    if (pdgcodeFirstMother != fMCfirstMotherPID) return kFALSE;
  }
  if ( fCutMCprocessType )
  {
    TParticle* particle = mcparticle->Particle();
    Int_t processID = particle->GetUniqueID();
    if (processID != fMCprocessType ) return kFALSE;
  }
  if (fCutMChasTrackReferences)
  {
    if (mcparticle->GetNumberOfTrackReferences()<1) return kFALSE;
  }
  return kTRUE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMCcuts()
{
  //check MC info
  if (!fMCevent) return kFALSE;
  if (fTrackLabel<0) return kFALSE;//otherwise AliCMevent prints a warning before returning NULL
  fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  return PassesMCcuts(fMCevent,fTrackLabel);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(const AliESDv0* v0)
{
  //check if the v0 passes cuts

  //clean up from last iteration
  ClearTrack();
  
  fV0 = const_cast<AliESDv0*>(v0);

  Bool_t pass = kTRUE;

  //first, extract/prepare the v0
  if (!v0->GetOnFlyStatus()) return kFALSE; //skip offline V0 finder
  const AliExternalTrackParam *negHelix=v0->GetParamN();
  const AliExternalTrackParam *posHelix=v0->GetParamP();
  AliVParticle *v0tracks[2];
  v0tracks[0] = fEvent->GetTrack(v0->GetNindex());
  v0tracks[1] = fEvent->GetTrack(v0->GetPindex());
  if( v0tracks[1]->Charge() < 0)
  {
    v0tracks[1] = fEvent->GetTrack(v0->GetNindex());
    v0tracks[0] = fEvent->GetTrack(v0->GetPindex());
    negHelix=v0->GetParamP();
    posHelix=v0->GetParamN();
  }

  int KalmanPidPairs[4][2] =
  {
    {-11,11}, // e-,e+ (gamma)
    {-211,2212}, // pi- p (lambda)
    {-2212,211}, // anti-p pi+ (anti-lambda)
    {-211,211} // pi-  pi+ (Kshort)
    //   {-321,321} // K- K+ (phi)
  };
    
  Int_t id = 3;
  //refit using a mass hypothesis
  AliKFParticle v0trackKFneg(*(negHelix),KalmanPidPairs[id][0]);
  AliKFParticle v0trackKFpos(*(posHelix),KalmanPidPairs[id][1]);
  AliKFParticle v0particleRefit;
  v0particleRefit += v0trackKFneg;
  v0particleRefit += v0trackKFpos;
  Double_t invMassErr= -999;
  v0particleRefit.GetMass(fTrackMass,invMassErr);
  //Double_t openAngle = v0trackKFneg.GetAngle(v0trackKFpos);
  fTrackEta = v0particleRefit.GetEta();
  fTrackPt = v0particleRefit.GetPt();
  fTrackPhi = TMath::Pi()+v0particleRefit.GetPhi();
  ////find out which mass bin and put the number in fPOItype
  //Int_t massBins = AliFlowCommonConstants::GetMaster()->GetNbinsMass();
  //Double_t massMin = AliFlowCommonConstants::GetMaster()->GetMassMin();
  //Double_t massMax = AliFlowCommonConstants::GetMaster()->GetMassMax();
  //fPOItype = 1 + int(massBins*(fTrackMass-massMin)/(massMax-massMin));
  
  /////////////////////////////////////////////////////////////////////////////
  //apply cuts
  if ( v0tracks[0]->Charge() == v0tracks[1]->Charge() ) pass=kFALSE;
  if ( v0tracks[0]->Pt()<0.15 || v0tracks[1]->Pt()<0.15 ) pass=kFALSE;

  return pass;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(const AliESDkink* kink)
{
  //check if the kink passes cuts
  
  //clean up from last iteration
  ClearTrack();
  fKink=const_cast<AliESDkink*>(kink);

  Bool_t pass = kTRUE;

  Float_t kinkAngle = kink->GetAngle(2);
  if (kinkAngle<fMinKinkAngle) pass = kFALSE;
  Double_t kinkRadius = kink->GetR();
  if (kinkRadius<fMinKinkRadius || kinkRadius>fMaxKinkRadius) pass = kFALSE;

  //Qt
  const TVector3 motherMfromKink(kink->GetMotherP());
  const TVector3 daughterMfromKink(kink->GetDaughterP());
  Float_t qt=kink->GetQt();
  if ( qt < fMinKinkQt || qt > fMaxKinkQt) pass = kFALSE;

  //invariant mass
  Float_t energyDaughterMu = TMath::Sqrt( daughterMfromKink.Mag()*daughterMfromKink.Mag()+
                                          0.105658*0.105658 );
  Float_t p1XM = motherMfromKink.Px();
  Float_t p1YM = motherMfromKink.Py();
  Float_t p1ZM = motherMfromKink.Pz();
  Float_t p2XM = daughterMfromKink.Px();
  Float_t p2YM = daughterMfromKink.Py();
  Float_t p2ZM = daughterMfromKink.Pz();
  Float_t p3Daughter = TMath::Sqrt( ((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+
                                    ((p1ZM-p2ZM)*(p1ZM-p2ZM)) );
  Double_t invariantMassKmu = TMath::Sqrt( (energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-
                                           motherMfromKink.Mag()*motherMfromKink.Mag() );
  if (invariantMassKmu > fMaxKinkInvMassKmu) pass=kFALSE;
  if (invariantMassKmu < fMinKinkInvMassKmu) pass=kFALSE;
  fTrackMass=invariantMassKmu;

  if (fQA)
  {
    QAbefore(13)->Fill(qt);
    if (pass) QAafter(13)->Fill(qt);
    QAbefore(14)->Fill(invariantMassKmu);
    if (pass) QAafter(14)->Fill(invariantMassKmu);
    const Double_t* kinkPosition = kink->GetPosition();
    QAbefore(15)->Fill(kinkPosition[0],kinkPosition[1]);
    if (pass) QAafter(15)->Fill(kinkPosition[0],kinkPosition[1]);
    QAbefore(16)->Fill(motherMfromKink.Mag(),kinkAngle*TMath::RadToDeg());
    if (pass) QAafter(16)->Fill(motherMfromKink.Mag(),kinkAngle*TMath::RadToDeg());
  }

  //mother track cuts
  Int_t indexKinkMother = kink->GetIndex(0);
  AliESDtrack* motherTrack = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(indexKinkMother));
  if (!motherTrack) return kFALSE;
  if (!PassesCuts(motherTrack)) pass = kFALSE;
  
  return pass;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesCuts(AliVParticle* vparticle)
{
  //check cuts for an ESD vparticle

  ////////////////////////////////////////////////////////////////
  //  start by preparing the track parameters to cut on //////////
  ////////////////////////////////////////////////////////////////
  //clean up from last iteration
  ClearTrack();
  Bool_t pass=kTRUE;

  //get the label and the mc particle
  fTrackLabel = (fFakesAreOK)?TMath::Abs(vparticle->GetLabel()):vparticle->GetLabel();
  if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  else fMCparticle=NULL;

  Bool_t isMCparticle = kFALSE; //some things are different for MC particles, check!
  AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*>(vparticle);
  AliAODTrack* aodTrack = NULL;
  if (esdTrack)
  {
    //for an ESD track we do some magic sometimes like constructing TPC only parameters
    //or doing some hybrid, handle that here
    HandleESDtrack(esdTrack);
  }
  else
  {
    HandleVParticle(vparticle);
    //now check if produced particle is MC
    isMCparticle = (dynamic_cast<AliMCParticle*>(fTrack))!=NULL;
    aodTrack = dynamic_cast<AliAODTrack*>(vparticle); //keep the additional dynamic cast out of the way for ESDs
  }
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  if (!fTrack) return kFALSE;
  //because it may be different from global, not needed for aodTrack because we dont do anything funky there
  if (esdTrack) esdTrack = static_cast<AliESDtrack*>(fTrack);
  
  //check the common cuts for the current particle fTrack (MC,AOD,ESD)
  Double_t pt = fTrack->Pt();
  Double_t p = fTrack->P();
  if (!fFakesAreOK) {if (fTrackLabel<0) pass=kFALSE;}
  if (fCutPt) {if (pt < fPtMin || pt >= fPtMax ) pass=kFALSE;}
  if (fCutEta) {if (fTrack->Eta() < fEtaMin || fTrack->Eta() >= fEtaMax ) pass=kFALSE;}
  if (fCutPhi) {if (fTrack->Phi() < fPhiMin || fTrack->Phi() >= fPhiMax ) pass=kFALSE;}
  if (fRequireCharge) {if (fTrack->Charge() == 0) pass=kFALSE;}
  if (fCutCharge && !isMCparticle) {if (fTrack->Charge() != fCharge) pass=kFALSE;}
  if (fCutCharge && isMCparticle)
  { 
    //in case of an MC particle the charge is stored in units of 1/3|e| 
    Int_t charge = TMath::Nint(fTrack->Charge()/3.0); //mc particles have charge in units of 1/3e
    if (charge!=fCharge) pass=kFALSE;
  }
  //if(fCutPID) {if (fTrack->PID() != fPID) pass=kFALSE;}

  //when additionally MC info is required
  if (fCutMC && !PassesMCcuts()) pass=kFALSE;

  //the case of ESD or AOD
  if (esdTrack) { if (!PassesESDcuts(esdTrack)) { pass=kFALSE; } }
  if (aodTrack) { if (!PassesAODcuts(aodTrack,pass)) { pass=kFALSE; } }

  if (fQA)
  {
    if (fMCparticle)
    {
      TParticle* tparticle=fMCparticle->Particle();
      Int_t processID = tparticle->GetUniqueID();
      Int_t firstMotherLabel = tparticle->GetFirstMother();
      Bool_t primary = IsPhysicalPrimary();
      //TLorentzVector v;
      //mcparticle->Particle()->ProductionVertex(v);
      //Double_t prodvtxX = v.X();
      //Double_t prodvtxY = v.Y();

      Int_t pdgcode = fMCparticle->PdgCode();
      AliVParticle* firstMotherParticle = fMCevent->GetTrack(firstMotherLabel);
      Int_t pdgcodeFirstMother = 0;
      if (firstMotherParticle) {pdgcodeFirstMother = firstMotherParticle->PdgCode();}
      
      Float_t pdg = 0;
      switch (TMath::Abs(pdgcode))
      {
        case 11:
          pdg = AliPID::kElectron + 0.5; break;
        case 13:
          pdg = AliPID::kMuon + 0.5; break;
        case 211:
          pdg = AliPID::kPion + 0.5; break;
        case 321:
          pdg = AliPID::kKaon + 0.5; break;
        case 2212:
          pdg = AliPID::kProton + 0.5; break;
        default:
          pdg = AliPID::kUnknown + 0.5; break;
      }
      pdg = TMath::Sign(pdg,static_cast<Float_t>(pdgcode));

      Float_t geantCode = 0;
      switch (pdgcodeFirstMother)
      {
        case 22: //#gamma
          geantCode=1;
          break;
        case -11: //e^{+}
          geantCode=2;
          break;
        case 11: //e^{-}
          geantCode=3;
          break;
        case 12: case 14: case 16: //#nu
          geantCode=4;
          break;
        case -13:
          geantCode=5; //#mu^{+}
          break;
        case 13:
          geantCode=6; //#mu^{-}
          break;
        case 111: //#pi^{0}
          geantCode=7;
          break;
        case 211: //#pi^{+}
          geantCode=8;
          break;
        case -211: //#pi^{-}
          geantCode=9;
          break;
        case 130: //K^{0}_{L}
          geantCode=10;
          break;
        case 321: //K^{+}
          geantCode=11;
          break;
        case -321: //K^{-}
          geantCode=12;
          break;
        case 2112:
          geantCode = 13; //n
          break;
        case 2212:
          geantCode = 14; //p
          break;
        case -2212:
          geantCode=15; //#bar{p}
          break;
        case 310:
          geantCode=16; //K^{0}_{S}
          break;
        case 221:
          geantCode=17; //#eta
          break;
        case 3122:
          geantCode=18; //#Lambda
          break;
        case 3222:
          geantCode=19; //#Sigma^{+}
          break;
        case 3212:
          geantCode=20; //#Sigma^{0}
          break;
        case 3112:
          geantCode=21; //#Sigma^{-}
          break;
        case 3322:
          geantCode=22; //#Theta^{0}
          break;
        case 3312:
          geantCode=23; //#Theta^{-}
          break;
        case 3332:
          geantCode=24; //#Omega^{-}
          break;
        case -2112:
          geantCode=25; //#bar{n}
          break;
        case -3122:
          geantCode=26; //#bar{#Lambda}
          break;
        case -3112:
          geantCode=27; //#bar{#Sigma}^{-}
          break;
        case -3212:
          geantCode=28; //#bar{#Sigma}^{0}
          break;
        case -3222:
          geantCode=29; //#bar{#Sigma}^{+}
          break;
        case -3322:
          geantCode=30; //#bar{#Theta}^{0}
          break;
        case -3312:
          geantCode=31; //#bar{#Theta}^{+}
          break;
        case -3332:
          geantCode=32; //#bar{#Omega}^{+}
          break;
        case -15:
          geantCode=33; //#tau^{+}
          break;
        case 15:
          geantCode=34; //#tau^{-}
          break;
        case 411:
          geantCode=35; //D^{+}
          break;
        case -411:
          geantCode=36; //D^{-}
          break;
        case 421:
          geantCode=37; //D^{0}
          break;
        case -421:
          geantCode=38; //#bar{D}^{0}
          break;
        case 431:
          geantCode=39; //D_{s}^{+}
          break;
        case -431:
          geantCode=40; //#bar{D_{s}}^{-}
          break;
        case 4122:
          geantCode=41; //#Lambda_{c}^{+}
          break;
        case 24:
          geantCode=42; //W^{+}
          break;
        case -24:
          geantCode=43; //W^{-}
          break;
        case 23:
          geantCode=44; //Z^{0}
          break;
        default:
          geantCode = 100;
          break;
      }
      
      QAbefore(2)->Fill(p,pdg);
      QAbefore(3)->Fill(p,primary?0.5:-0.5);
      QAbefore(4)->Fill(p,static_cast<Float_t>(processID));
      QAbefore(7)->Fill(p,geantCode+0.5);
      if (pass) QAafter(2)->Fill(p,pdg);
      if (pass) QAafter(3)->Fill(p,primary?0.5:-0.5);
      if (pass) QAafter(4)->Fill(p,static_cast<Float_t>(processID));
      if (pass) QAafter(7)->Fill(p,geantCode);
    }
  }

  //true by default, if we didn't set any cuts
  return pass;
}

//_______________________________________________________________________
Int_t AliFlowTrackCuts::GetITStype(const AliAODTrack* track) const
{
  // set ITStype (test purpose only)
  Int_t ITStype=0;
  Bool_t bHITs[6] = {kFALSE};
  Int_t hitcnt=0;
  Bool_t isSPD=kFALSE, isSDD=kFALSE, isSSD=kFALSE;
  for (Int_t i=0; i<6; i++) {
    if(track->HasPointOnITSLayer(i)) {
      if(i<2) isSPD=kTRUE;
      if(i>2 && i<4) isSDD=kTRUE;
      if(i>4) isSSD=kTRUE;
      bHITs[i]=kTRUE;
      hitcnt++;
    }
  }
  if(hitcnt==1) {
    for (Int_t i=0; i<6; i++) {
      if(bHITs[i]) ITStype = i+1;
    }
  }
  if(hitcnt>=2) {
    if(isSPD  & !isSDD & !isSSD) ITStype = 7;
    if(isSPD  &  isSDD & !isSSD) ITStype = 8;
    if(isSPD  & !isSDD &  isSSD) ITStype = 9;
    if(!isSPD & isSDD  & !isSSD) ITStype = 10;
    if(!isSPD & isSDD  &  isSSD) ITStype = 11;
    if(!isSPD & !isSDD &  isSSD) ITStype = 12;
    if(isSPD  &  isSDD &  isSSD) {
      if(bHITs[0]  && !bHITs[1]) ITStype = 13;
      if(!bHITs[0] &&  bHITs[1]) ITStype = 14;
      if(bHITs[0]  &&  bHITs[1]) ITStype = 15;
    }
  }
  return ITStype;
}

//_______________________________________________________________________
Bool_t AliFlowTrackCuts::PassesAODcuts(const AliAODTrack* track, Bool_t passedFid)
{
  //check cuts for AOD
  Bool_t pass = passedFid;
//  AliAODPid *pidObj = track->GetDetPid();

  if (fCutNClustersTPC)
  {
    Int_t ntpccls = track->GetTPCNcls();
    if (ntpccls < fNClustersTPCMin || ntpccls > fNClustersTPCMax) pass=kFALSE;
  }

  if (fCutNClustersITS)
  {
    Int_t nitscls = track->GetITSNcls();
    if (nitscls < fNClustersITSMin || nitscls > fNClustersITSMax) pass=kFALSE;
  }
  
  if (fCutChi2PerClusterTPC)
  {
    Double_t chi2tpc = track->Chi2perNDF();
    if (chi2tpc < fMinChi2PerClusterTPC || chi2tpc > fMaxChi2PerClusterTPC) pass=kFALSE;
  }
  
  if (fCutChi2PerClusterITS)
  {
    Int_t nitscls = track->GetITSNcls();
    if(nitscls>0) {
      Double_t chi2TIS = track->GetITSchi2()/track->GetITSNcls();
      if (chi2TIS >= fMaxChi2PerClusterITS) pass=kFALSE;
    }
  }
  
  if (fCutITSClusterGlobal)
  {
    Bool_t bSPDCl = kFALSE;
    for (Int_t i=0; i<2; i++) {
      if(track->HasPointOnITSLayer(i)) bSPDCl = kTRUE;
    }
    Bool_t bSDDCl = track->HasPointOnITSLayer(2);
    Bool_t temppass = kTRUE;
    if (bSPDCl || (!bSPDCl && bSDDCl)) temppass = kTRUE;
    else temppass = kFALSE;
    if(!temppass) pass=kFALSE;
  }
  
  if (fCutFracSharedTPCCluster)
  {
    Int_t ntpccls = track->GetTPCncls();
    if(ntpccls>0) {
      Int_t ntpcclsS = track->GetTPCnclsS();
      Double_t fshtpccls = 1.*ntpcclsS/ntpccls;
      if (fshtpccls > fMaxFracSharedTPCCluster) pass=kFALSE;
    }
  }
  
  if (fCutCrossedTPCRows)
  {
    Int_t nCrossedRows = track->GetTPCNCrossedRows();
    if (nCrossedRows <= fMinNCrossedRows) pass=kFALSE;
    Float_t CrossedRowsOverFindableClusters = track->GetTPCFoundFraction();
    if (CrossedRowsOverFindableClusters < fMinCrossedRowsOverFindableClusters) pass=kFALSE;
  }
  
  if (fCutGoldenChi2)
  {
    Double_t GoldenChi2 = track->GetChi2TPCConstrainedVsGlobal();
    if (GoldenChi2 >= fMaxGoldenChi2) pass=kFALSE;
  }
  
  if (fRequireTOFSignal)
  {
    if(TMath::Abs(track->GetTOFsignalDz())>10.)  pass=kFALSE;
    if(track->GetTOFsignal() < 12000.)  pass=kFALSE;
    if(track->GetTOFsignal() > 25000.)  pass=kFALSE;
  }
  
  if (GetRequireTPCRefit() && !(track->GetStatus() & AliESDtrack::kTPCrefit) ) pass=kFALSE;
  if (GetRequireITSRefit() && !(track->GetStatus() & AliESDtrack::kITSrefit) ) pass=kFALSE;
  
  if (fUseAODFilterBit && !track->TestFilterBit(fAODFilterBit)) pass=kFALSE;
  Double_t DCAxy = track->DCA();
  Double_t DCAz = track->ZAtDCA();
  if(fCutDCAToVertexXYAOD || fCutDCAToVertexZAOD || fCutDCAToVertexXYPtDepAOD) {
    if (std::abs((Int_t)DCAxy)==999 || std::abs((Int_t)DCAz)==999) {
      // re-evaluate the dca as it seems to not be natively present
      // allowed only for tracks inside the beam pipe
      Double_t pos[3] = {-99., -99., -99.};
      track->GetPosition(pos);
      if(pos[0]*pos[0]+pos[1]*pos[1] <= 3.*3.) {
        AliAODTrack copy(*track);       // stack copy
        Double_t b[2] = {-99., -99.};
        Double_t bCov[3] = {-99., -99., -99.};
        if(copy.PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 100., b, bCov)) {
          DCAxy = b[0];
          DCAz = b[1];
        }
      }
    }
    if (fCutDCAToVertexXYAOD) {
      if (TMath::Abs(DCAxy)>fMaxDCAxyAOD) pass=kFALSE;
    }
    if (fCutDCAToVertexXYPtDepAOD) {
      Double_t MaxDCAPtDep = 2.4;
      if (fAODFilterBit!=128) MaxDCAPtDep = 0.0182+0.0350/pow(track->Pt(),1.01);
      else MaxDCAPtDep = 0.4+0.2/pow(track->Pt(),0.3);
      if (TMath::Abs(DCAxy)>MaxDCAPtDep) pass=kFALSE;
    }
    if (fCutDCAToVertexZAOD) {
      if (TMath::Abs(DCAz)>fMaxDCAzAOD) pass=kFALSE;
    }
  }
  Double_t dedx = track->GetTPCsignal();
  if(fCutMinimalTPCdedx) {
    if (dedx < fMinimalTPCdedx) pass=kFALSE;
  }
  Double_t time[9];
  track->GetIntegratedTimes(time);
  
  if(fCutTPCSecbound && track->Pt()>fCutTPCSecboundMinpt) {
    Double_t xyz[3]={-9999.,-9999.,-9999.};
    const double r = 84.; // TPC IROC radius
    if (!track->GetXYZatR(r, fEvent->GetMagneticField(), xyz, NULL)) pass=kFALSE;
    Double_t cra = TMath::ATan2(xyz[1],xyz[0]); // crossing angle
    Double_t dpe = 3.*TMath::TwoPi()/360.;// excluded region (\pm 3 degree)
    for(Int_t nb=-9; nb<=9; nb++) {
      Double_t bnp = nb*TMath::Pi()/9.; // TPC boundaries azimuthal angle
      if(cra<bnp+dpe && cra>bnp-dpe) pass=kFALSE;
    }
  }
  if(fCutTPCSecboundVar) {
    Double_t phimod = track->Phi();
    if(fEvent->GetMagneticField() < 0)   // for negative polarity field
      phimod = TMath::TwoPi() - phimod;
    if(track->Charge() < 0) // for negative charge
      phimod = TMath::TwoPi() - phimod;
    if(phimod < 0)
      cout << "Warning!!!!! phi < 0: " << phimod << endl;
    
    phimod += TMath::Pi()/18.0; // to center gap in the middle
    phimod = fmod(phimod, TMath::Pi()/9.0);
    if(phimod < fPhiCutHigh->Eval(track->Pt()) && phimod > fPhiCutLow->Eval(track->Pt()))
      pass=kFALSE; // reject track
  }

  if (fCutPID && (fParticleID!=AliPID::kUnknown)) //if kUnknown don't cut on PID
  {
    if (!PassesAODpidCut(track)) pass=kFALSE;
  }

  if (fQA) {
    // changed 04062014 used to be filled before possible PID cut
    Double_t momTPC = track->GetTPCmomentum();
    QAbefore( 0)->Fill(momTPC,GetBeta(track, kTRUE));
    if(pass) QAafter( 0)->Fill(momTPC, GetBeta(track, kTRUE));
    QAbefore( 1)->Fill(momTPC,dedx);
    QAbefore( 5)->Fill(track->Pt(),DCAxy);
    QAbefore( 6)->Fill(track->Pt(),DCAz);
    if (pass) QAafter( 1)->Fill(momTPC,dedx);
    if (pass) QAafter( 5)->Fill(track->Pt(),DCAxy);
    if (pass) QAafter( 6)->Fill(track->Pt(),DCAz);
    QAbefore( 8)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kElectron]));
    if (pass) QAafter(  8)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kElectron]));
    QAbefore( 9)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kMuon]));
    if (pass) QAafter(  9)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kMuon]));
    QAbefore(10)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kPion]));
    if (pass) QAafter( 10)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kPion]));
    QAbefore(11)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kKaon]));
    if (pass) QAafter( 11)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kKaon]));
    QAbefore(12)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kProton]));
    if (pass) QAafter( 12)->Fill(track->P(),(track->GetTOFsignal()-time[AliPID::kProton]));
    QAbefore(18)->Fill(track->P(),track->Chi2perNDF());
    if (pass) QAafter( 18)->Fill(track->P(),track->Chi2perNDF());
    Int_t ntpccls = track->GetTPCncls();
    if(ntpccls>0) {
      Int_t ntpcclsS = track->GetTPCnclsS();
      Double_t fshtpccls = 1.*ntpcclsS/ntpccls;
      QAbefore(19)->Fill(track->P(),fshtpccls);
      if (pass) QAafter(19)->Fill(track->P(),fshtpccls);
    }
    QAbefore(20)->Fill(track->P(),ntpccls);
    if (pass) QAafter(20)->Fill(track->P(),ntpccls);
    Int_t nitscls = track->GetITSNcls();
    if(nitscls>0) {
      Double_t chi2TIS = track->GetITSchi2()/track->GetITSNcls();
      QAbefore(21)->Fill(track->P(),chi2TIS);
      if (pass) QAafter(21)->Fill(track->P(),chi2TIS);
    }
  }

  return pass;
}

//_______________________________________________________________________
Bool_t AliFlowTrackCuts::PassesESDcuts(AliESDtrack* track)
{
  //check cuts on ESD tracks
  Bool_t pass=kTRUE;
  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;
  track->GetImpactParameters(dcaxy,dcaz);
  const AliExternalTrackParam* pout = track->GetOuterParam();
  const AliExternalTrackParam* pin = track->GetInnerParam();
  
  if (fIgnoreTPCzRange)
  {
    if (pin&&pout)
    {
      Double_t zin = pin->GetZ();
      Double_t zout = pout->GetZ();
      if (zin*zout<0) pass=kFALSE;   //reject if cross the membrane
      if (zin < fIgnoreTPCzRangeMin || zin > fIgnoreTPCzRangeMax) pass=kFALSE;
      if (zout < fIgnoreTPCzRangeMin || zout > fIgnoreTPCzRangeMax) pass=kFALSE;
    }
  }
 
  Int_t ntpccls = ( fParamType==kTPCstandalone )?
                    track->GetTPCNclsIter1():track->GetTPCNcls();    
  if (fCutChi2PerClusterTPC)
  {
    Float_t tpcchi2 = (fParamType==kTPCstandalone)?
                       track->GetTPCchi2Iter1():track->GetTPCchi2();
    tpcchi2 = (ntpccls>0)?tpcchi2/ntpccls:-FLT_MAX;
    if (tpcchi2<fMinChi2PerClusterTPC || tpcchi2 >=fMaxChi2PerClusterTPC)
      pass=kFALSE;
  }

  if (fCutMinimalTPCdedx) 
  {
    if (track->GetTPCsignal() < fMinimalTPCdedx) pass=kFALSE;
  }

  if (fCutNClustersTPC)
  {
    if (ntpccls < fNClustersTPCMin || ntpccls > fNClustersTPCMax) pass=kFALSE;
  }

  Int_t nitscls = track->GetNcls(0);
  if (fCutNClustersITS)
  {
    if (nitscls < fNClustersITSMin || nitscls > fNClustersITSMax) pass=kFALSE;
  }

  //some stuff is still handled by AliESDtrackCuts class - delegate
  if (fAliESDtrackCuts)
  {
    if (!fAliESDtrackCuts->IsSelected(track)) pass=kFALSE;
  }
 
  //PID part with pid QA
  Double_t beta = GetBeta(track);
  Double_t dedx = Getdedx(track);
  if (fQA)
  {
    if (pass) QAbefore(0)->Fill(track->GetP(),beta);
    if (pass) QAbefore(1)->Fill(pin->GetP(),dedx);
  }
  if (fCutPID && (fParticleID!=AliPID::kUnknown)) //if kUnknown don't cut on PID
  {
    if (!PassesESDpidCut(track)) pass=kFALSE;
  }
  if (fCutRejectElectronsWithTPCpid)
  {
    //reject electrons using standard TPC pid
    //TODO this should be rethought....
    Double_t pidTPC[AliPID::kSPECIES];
    track->GetTPCpid(pidTPC);
    if (pidTPC[AliPID::kElectron]<fParticleProbability) pass=kFALSE;
  }
  if(fCutITSclusterShared) {  
    Int_t counterForSharedCluster=MaxSharedITSClusterCuts(track);
    if(counterForSharedCluster >= fMaxITSclusterShared) pass=kFALSE;
  }
  
  if(fCutITSChi2) {  
    Double_t chi2perClusterITS = MaxChi2perITSClusterCuts(track);
    if(chi2perClusterITS >= fMaxITSChi2) pass=kFALSE;
  }
    
  if (fQA)
  {
    if (pass) QAafter(0)->Fill(track->GetP(),beta);
    if (pass) QAafter(1)->Fill(pin->GetP(),dedx);
  }
  //end pid part with qa
  
  if (fQA)
  {
    Double_t pt = track->Pt();
    QAbefore(5)->Fill(pt,dcaxy);
    QAbefore(6)->Fill(pt,dcaz);
    if (pass) QAafter(5)->Fill(pt,dcaxy);
    if (pass) QAafter(6)->Fill(pt,dcaz);
    QAbefore(17)->Fill(Float_t(track->GetKinkIndex(0)));
    if (pass) QAafter(17)->Fill(Float_t(track->GetKinkIndex(0)));
  }

  return pass;
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::HandleVParticle(AliVParticle* track)
{
  //handle the general case
  switch (fParamType)
  {
    default:
      fTrack = track;
      break;
  }
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::HandleESDtrack(AliESDtrack* track)
{
  //handle esd track
  switch (fParamType)
  {
    case kGlobal:
      fTrack = track;
      break;
    case kTPCstandalone:
      if (!track->FillTPCOnlyTrack(fTPCtrack)) 
      {
        fTrack=NULL;
        fMCparticle=NULL;
        fTrackLabel=-998;
        return;
      }
      fTrack = &fTPCtrack;
      //recalculate the label and mc particle, they may differ as TPClabel != global label
      fTrackLabel = (fFakesAreOK)?TMath::Abs(fTrack->GetLabel()):fTrack->GetLabel();
      if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
      else fMCparticle=NULL;
      break;
    default:
      if (fForceTPCstandalone)
      {
        if (!track->FillTPCOnlyTrack(fTPCtrack))
        {
          fTrack=NULL;
          fMCparticle=NULL;
          fTrackLabel=-998;
          return;
        }
        fTrack = &fTPCtrack;
        //recalculate the label and mc particle, they may differ as TPClabel != global label
        fTrackLabel = (fFakesAreOK)?TMath::Abs(fTrack->GetLabel()):fTrack->GetLabel();
        if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
        else fMCparticle=NULL;
      }
      else 
        fTrack = track;
      break;
  }
}

//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::Count(AliVEvent* event)
{
  //calculate the number of track in given event.
  //if argument is NULL(default) take the event attached
  //by SetEvent()
    
  Int_t multiplicity = 0;
  if (!event)
  {
    for (Int_t i=0; i<GetNumberOfInputObjects(); i++)
    {
      if (IsSelected(GetInputObject(i))) multiplicity++;
    }
  }
  else
  {
    for (Int_t i=0; i<event->GetNumberOfTracks(); i++)
    {
      if (IsSelected(event->GetTrack(i))) multiplicity++;
    }
  }
  return multiplicity;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetAODTrackCutsForFilterBit(UInt_t bit, TString suffix)
{
  // object which in its default form only cuts on filterbit (for AOD analysis)
  // NOTE : the cut object is defined in such a way that ONLY filterbit is tested, 
  // all other paramters are NOT checked 
  // the exception is TPCdecx which is always cut for reasons of backward compatibility
  // but by setting the threshold to larger than -99999999 effectively there is no check
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts(Form("AOD fitlerbit %i, %s", (int)bit, suffix.Data()));
  cuts->SetMinimalTPCdedx(-999999999);
  cuts->SetAODfilterBit(bit);
  cuts->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  return cuts;
}  
//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts()
{
  //returns the lhc10h vzero track cuts, this function
  //is left here for backward compatibility
  //if a run is recognized as 11h, the calibration method will
  //switch to 11h calbiration, which means that the cut 
  //object is updated but not replaced.
  //calibratin is only available for PbPb runs
  return GetStandardVZEROOnlyTrackCuts2010();
}
//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2010()
{
  //get standard VZERO cuts
  //DISCLAIMER: LHC10h VZERO calibration consists (by default) of two steps
  //1) re-weigting of signal
  //2) re-centering of q-vectors
  //step 2 is available only for n==2 and n==3, for the higher harmonics the user
  //is repsonsible for making sure the q-sub distributions are (sufficiently) flat
  //or a sensible NUA procedure is applied !
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard vzero flow cuts");
  cuts->SetParamType(AliFlowTrackCuts::kVZERO);
  cuts->SetEtaRange( -10, +10 );
  cuts->SetEtaGap(-1., 1.);
  cuts->SetPhiMin( 0 );
  cuts->SetPhiMax( TMath::TwoPi() );
  // options for the reweighting
  cuts->SetVZEROgainEqualizationPerRing(kFALSE);
  cuts->SetApplyRecentering(kTRUE);
  // to exclude a ring , do e.g. 
  // cuts->SetUseVZERORing(7, kFALSE);
  // excluding a ring will break the re-centering as re-centering relies on a 
  // database file which tuned to receiving info from all rings
  return cuts;
}
//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts2011()
{
  //get standard VZERO cuts for 2011 data
  //in this case, the vzero segments will be weighted by
  //VZEROEqMultiplicity, 
  //if recentering is enableded, the sub-q vectors
  //will be taken from the event header, so make sure to run 
  //the VZERO event plane selection task before this task !
  //DISCLAIMER: recentering is only available for n==2
  //for the higher harmonics the user
  //is repsonsible for making sure the q-sub distributions are (sufficiently) flat
  //or a sensible NUA procedure is applied !
  //recentering replaces the already evaluated q-vectors, so 
  //when chosen, additional settings (e.g. excluding rings) 
  //have no effect. recentering is true by default
  //
  //NOTE user is responsible for running the vzero event plane
  //selection task in advance, e.g. add to your launcher macro
  //
  //  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  //  AddTaskVZEROEPSelection();
  //
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard vzero flow cuts 2011");
  cuts->SetParamType(kVZERO);
  cuts->SetEtaRange( -10, +10 );
  cuts->SetEtaGap(-1., 1.);
  cuts->SetPhiMin( 0 );
  cuts->SetPhiMax( TMath::TwoPi() );
  cuts->SetApplyRecentering(kTRUE);
  cuts->SetVZEROgainEqualizationPerRing(kFALSE);
 return cuts;
}
//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetBetaVZEROOnlyTrackCuts()
{
  // test patch for VZERO track cuts
  // april 20 2015. 
  // DISCLAIMER: BETA TEST OF CODE, DO NOT USE FOR PHYSICS RESULTS !!!
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("beta version of vzero ");
  cuts->SetParamType(AliFlowTrackCuts::kBetaVZERO);
  cuts->SetEtaRange( -10, +10 );
  cuts->SetEtaGap(-1., 1.);
  cuts->SetPhiMin( 0 );
  cuts->SetPhiMax( TMath::TwoPi() );
  cuts->SetApplyRecentering(kTRUE);
  // idea of the cuts is that calibration is done per ring
  // and that it is transparent for different data taking periods
  return cuts;
}


//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardGlobalTrackCuts2010()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard Global tracks");
  cuts->SetParamType(kGlobal);
  cuts->SetPtRange(0.2,5.);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinChi2PerClusterTPC(0.1);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetMinNClustersITS(2);
  cuts->SetRequireITSRefit(kTRUE);
  cuts->SetRequireTPCRefit(kTRUE);
  cuts->SetMaxDCAToVertexXY(0.3);
  cuts->SetMaxDCAToVertexZ(0.3);
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMinimalTPCdedx(10.);

  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard TPC standalone 2010");
  cuts->SetParamType(kTPCstandalone);
  cuts->SetPtRange(0.2,5.);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinChi2PerClusterTPC(0.2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetMaxDCAToVertexXY(3.0);
  cuts->SetMaxDCAToVertexZ(3.0);
  cuts->SetDCAToVertex2D(kTRUE);
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMinimalTPCdedx(10.);

  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts()
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard TPC standalone");
  cuts->SetParamType(kTPCstandalone);
  cuts->SetPtRange(0.2,5.);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->SetMinNClustersTPC(70);
  cuts->SetMinChi2PerClusterTPC(0.2);
  cuts->SetMaxChi2PerClusterTPC(4.0);
  cuts->SetMaxDCAToVertexXY(3.0);
  cuts->SetMaxDCAToVertexZ(3.0);
  cuts->SetDCAToVertex2D(kTRUE);
  cuts->SetAcceptKinkDaughters(kFALSE);
  cuts->SetMinimalTPCdedx(10.);

  return cuts;
}

//-----------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardITSTPCTrackCuts2009(Bool_t selPrimaries)
{
  //get standard cuts
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard global track cuts 2009");
  delete cuts->fAliESDtrackCuts;
  cuts->fAliESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selPrimaries);
  cuts->SetParamType(kGlobal);
  return cuts;
}

//-----------------------------------------------------------------------------
AliFlowTrackCuts* AliFlowTrackCuts::GetStandardMuonTrackCuts(Bool_t isMC, Int_t passN)
{
// XZhang 20120604
  AliFlowTrackCuts* cuts = new AliFlowTrackCuts("standard muon track cuts");
  cuts->SetParamType(kMUON);
  cuts->SetStandardMuonTrackCuts();
  cuts->SetIsMuonMC(isMC);
  cuts->SetMuonPassNumber(passN);
  return cuts;
}

//-----------------------------------------------------------------------
//AliFlowTrack* AliFlowTrackCuts::FillFlowTrackV0(TObjArray* trackCollection, Int_t trackIndex) const
//{
//  //fill flow track from a reconstructed V0 (topological)
//  AliFlowTrack* flowtrack = static_cast<AliFlowTrack*>(trackCollection->At(trackIndex));
//  if (flowtrack)
//  {
//    flowtrack->Clear();
//  }
//  else 
//  {
//    flowtrack = new AliFlowCandidateTrack();
//    trackCollection->AddAtAndExpand(flowtrack,trackIndex);
//  }
//
//  TParticle *tmpTParticle=NULL;
//  AliMCParticle* tmpAliMCParticle=NULL;
//  AliExternalTrackParam* externalParams=NULL;
//  AliESDtrack* esdtrack=NULL;
//  switch(fParamMix)
//  {
//    case kPure:
//      flowtrack->Set(fTrack);
//      break;
//    case kTrackWithMCkine:
//      flowtrack->Set(fMCparticle);
//      break;
//    case kTrackWithMCPID:
//      flowtrack->Set(fTrack);
//      //flowtrack->setPID(...) from mc, when implemented
//      break;
//    case kTrackWithMCpt:
//      if (!fMCparticle) return NULL;
//      flowtrack->Set(fTrack);
//      flowtrack->SetPt(fMCparticle->Pt());
//      break;
//    case kTrackWithPtFromFirstMother:
//      if (!fMCparticle) return NULL;
//      flowtrack->Set(fTrack);
//      tmpTParticle = fMCparticle->Particle();
//      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
//      flowtrack->SetPt(tmpAliMCParticle->Pt());
//      break;
//    case kTrackWithTPCInnerParams:
//      esdtrack = dynamic_cast<AliESDtrack*>(fTrack);
//      if (!esdtrack) return NULL;
//      externalParams = const_cast<AliExternalTrackParam*>(esdtrack->GetTPCInnerParam());
//      if (!externalParams) return NULL;
//      flowtrack->Set(externalParams);
//      break;
//    default:
//      flowtrack->Set(fTrack);
//      break;
//  }
//  if (fParamType==kMC) 
//  {
//    flowtrack->SetSource(AliFlowTrack::kFromMC);
//    flowtrack->SetID(fTrackLabel);
//  }
//  else if (dynamic_cast<AliAODTrack*>(fTrack))
//  {
//    if (fParamType==kMUON)                            // XZhang 20120604
//      flowtrack->SetSource(AliFlowTrack::kFromMUON);  // XZhang 20120604
//    else                                              // XZhang 20120604
//      flowtrack->SetSource(AliFlowTrack::kFromAOD);   // XZhang 20120604
//    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
//  }
//  else if (dynamic_cast<AliMCParticle*>(fTrack)) 
//  {
//    flowtrack->SetSource(AliFlowTrack::kFromMC);
//    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
//  }
//
//  if (fV0)
//  {
//    //add daughter indices
//  }
//
//  return flowtrack;
//}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::FillFlowTrackKink(TObjArray* trackCollection, Int_t trackIndex) const
{
  //fill flow track from AliVParticle (ESD,AOD,MC)
  AliFlowTrack* flowtrack = static_cast<AliFlowTrack*>(trackCollection->At(trackIndex));
  if (flowtrack)
  {
    flowtrack->Clear();
  }
  else 
  {
    flowtrack = new AliFlowCandidateTrack();
    trackCollection->AddAtAndExpand(flowtrack,trackIndex);
  }

  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  AliExternalTrackParam* externalParams=NULL;
  AliESDtrack* esdtrack=NULL;
  switch(fParamMix)
  {
    case kPure:
      flowtrack->Set(fTrack);
      break;
    case kTrackWithMCkine:
      flowtrack->Set(fMCparticle);
      break;
    case kTrackWithMCPID:
      flowtrack->Set(fTrack);
      //flowtrack->setPID(...) from mc, when implemented
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return NULL;
      flowtrack->Set(fTrack);
      flowtrack->SetPt(fMCparticle->Pt());
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return NULL;
      flowtrack->Set(fTrack);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      break;
    case kTrackWithTPCInnerParams:
      esdtrack = dynamic_cast<AliESDtrack*>(fTrack);
      if (!esdtrack) return NULL;
      externalParams = const_cast<AliExternalTrackParam*>(esdtrack->GetTPCInnerParam());
      if (!externalParams) return NULL;
      flowtrack->Set(externalParams);
      break;
    default:
      flowtrack->Set(fTrack);
      break;
  }
  if (fParamType==kMC) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(fTrackLabel);
  }
  else if (dynamic_cast<AliESDtrack*>(fTrack))
  {
    flowtrack->SetSource(AliFlowTrack::kFromESD);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliESDMuonTrack*>(fTrack))                                  // XZhang 20120604
  {                                                                                 // XZhang 20120604
    flowtrack->SetSource(AliFlowTrack::kFromMUON);                                  // XZhang 20120604
    flowtrack->SetID((Int_t)static_cast<AliESDMuonTrack*>(fTrack)->GetUniqueID());  // XZhang 20120604
  }                                                                                 // XZhang 20120604
  else if (dynamic_cast<AliAODTrack*>(fTrack))
  {
    if (fParamType==kMUON)                            // XZhang 20120604
      flowtrack->SetSource(AliFlowTrack::kFromMUON);  // XZhang 20120604
    else                                              // XZhang 20120604
      flowtrack->SetSource(AliFlowTrack::kFromAOD);   // XZhang 20120604
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliMCParticle*>(fTrack)) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }

  if (fKink)
  {
    Int_t indexMother = fKink->GetIndex(0);
    Int_t indexDaughter = fKink->GetIndex(1);
    flowtrack->SetID(indexMother);
    flowtrack->AddDaughter(indexDaughter);
    flowtrack->SetMass(1.);
    flowtrack->SetSource(AliFlowTrack::kFromKink);
  }

  flowtrack->SetMass(fTrackMass);

  return flowtrack;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::FillFlowTrackGeneric(TObjArray* trackCollection, Int_t trackIndex) const
{
  //fill a flow track from tracklet,vzero,pmd,...
  AliFlowTrack* flowtrack = static_cast<AliFlowTrack*>(trackCollection->At(trackIndex));
  if (flowtrack)
  {
    flowtrack->Clear();
  }
  else 
  {
    flowtrack = new AliFlowTrack();
    trackCollection->AddAtAndExpand(flowtrack,trackIndex);
  }
  
  if (FillFlowTrackGeneric(flowtrack)) return flowtrack;
  else 
  {
    trackCollection->RemoveAt(trackIndex);
    delete flowtrack;
    return NULL;
  }
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::FillFlowTrackGeneric(AliFlowTrack* flowtrack) const
{
  //fill a flow track from tracklet,vzero,pmd,...
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  switch (fParamMix)
  {
    case kPure:
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt(fTrackPt);
      flowtrack->SetMass(fTrackMass);
      break;
    case kTrackWithMCkine:
      if (!fMCparticle) return kFALSE;
      flowtrack->SetPhi( fMCparticle->Phi() );
      flowtrack->SetEta( fMCparticle->Eta() );
      flowtrack->SetPt( fMCparticle->Pt() );
      flowtrack->SetMass(fTrackMass);
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return kFALSE;
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt(fMCparticle->Pt());
      flowtrack->SetMass(fTrackMass);
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return kFALSE;
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      flowtrack->SetMass(fTrackMass);
      break;
    default:
      flowtrack->SetPhi(fTrackPhi);
      flowtrack->SetEta(fTrackEta);
      flowtrack->SetWeight(fTrackWeight);
      flowtrack->SetPt(fTrackPt);
      flowtrack->SetMass(fTrackMass);
      break;
  }
  flowtrack->SetSource(AliFlowTrack::kFromTracklet);
  return kTRUE;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::FillFlowTrackVParticle(TObjArray* trackCollection, Int_t trackIndex) const
{
  //fill flow track from AliVParticle (ESD,AOD,MC)
  
  AliFlowTrack* flowtrack = static_cast<AliFlowTrack*>(trackCollection->At(trackIndex));
  if (flowtrack)
  {
    flowtrack->Clear();
  }
  else 
  {
    flowtrack = new AliFlowTrack();
    trackCollection->AddAtAndExpand(flowtrack,trackIndex);
  }

  if (FillFlowTrackVParticle(flowtrack)) return flowtrack;
  else
  {
    trackCollection->RemoveAt(trackIndex);
    delete flowtrack;
    return NULL;
  }
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::FillFlowTrackVParticle(AliFlowTrack* flowtrack) const
{
  //fill flow track from AliVParticle (ESD,AOD,MC)
  if (!fTrack) return kFALSE;
  TParticle *tmpTParticle=NULL;
  AliMCParticle* tmpAliMCParticle=NULL;
  AliExternalTrackParam* externalParams=NULL;
  AliESDtrack* esdtrack=NULL;
  switch(fParamMix)
  {
    case kPure:
      flowtrack->Set(fTrack);
      break;
    case kTrackWithMCkine:
      flowtrack->Set(fMCparticle);
      break;
    case kTrackWithMCPID:
      flowtrack->Set(fTrack);
      //flowtrack->setPID(...) from mc, when implemented
      break;
    case kTrackWithMCpt:
      if (!fMCparticle) return kFALSE;
      flowtrack->Set(fTrack);
      flowtrack->SetPt(fMCparticle->Pt());
      break;
    case kTrackWithPtFromFirstMother:
      if (!fMCparticle) return kFALSE;
      flowtrack->Set(fTrack);
      tmpTParticle = fMCparticle->Particle();
      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
      flowtrack->SetPt(tmpAliMCParticle->Pt());
      break;
    case kTrackWithTPCInnerParams:
      esdtrack = dynamic_cast<AliESDtrack*>(fTrack);
      if (!esdtrack) return kFALSE;
      externalParams = const_cast<AliExternalTrackParam*>(esdtrack->GetTPCInnerParam());
      if (!externalParams) return kFALSE;
      flowtrack->Set(externalParams);
      break;
    default:
      flowtrack->Set(fTrack);
      break;
  }
  if (fParamType==kMC) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(fTrackLabel);
  }
  else if (dynamic_cast<AliESDtrack*>(fTrack))
  {
    flowtrack->SetSource(AliFlowTrack::kFromESD);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  else if (dynamic_cast<AliESDMuonTrack*>(fTrack))                                  // XZhang 20120604
  {                                                                                 // XZhang 20120604
    flowtrack->SetSource(AliFlowTrack::kFromMUON);                                  // XZhang 20120604
    flowtrack->SetID((Int_t)static_cast<AliESDMuonTrack*>(fTrack)->GetUniqueID());  // XZhang 20120604
  }                                                                                 // XZhang 20120604
  else if (dynamic_cast<AliAODTrack*>(fTrack))
  {
    if (fParamType==kMUON)                            // XZhang 20120604
      flowtrack->SetSource(AliFlowTrack::kFromMUON);  // XZhang 20120604
    else                                              // XZhang 20120604
      flowtrack->SetSource(AliFlowTrack::kFromAOD);   // XZhang 20120604
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(fTrack);
    flowtrack->SetITStype(GetITStype(aodTrack));
  }
  else if (dynamic_cast<AliMCParticle*>(fTrack)) 
  {
    flowtrack->SetSource(AliFlowTrack::kFromMC);
    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
  }
  flowtrack->SetMass(fTrackMass);
  return kTRUE;
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrackCuts::FillFlowTrack(TObjArray* trackCollection, Int_t trackIndex) const
{
  //fill a flow track constructed from whatever we applied cuts on
  //return true on success
  switch (fParamType)
  {
    case kSPDtracklet:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kPMD:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kVZERO:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kBetaVZERO:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kDeltaVZERO:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kKappaVZERO:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kHotfixHI:
      return FillFlowTrackGeneric(trackCollection, trackIndex);
    case kKink:
      return FillFlowTrackKink(trackCollection, trackIndex);
    //case kV0:
    //  return FillFlowTrackV0(trackCollection, trackIndex);
    default:
      return FillFlowTrackVParticle(trackCollection, trackIndex);
  }
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::FillFlowTrack(AliFlowTrack* track) const
{
  //fill a flow track constructed from whatever we applied cuts on
  //return true on success
  switch (fParamType)
  {
    case kSPDtracklet:
      return FillFlowTrackGeneric(track);
    case kPMD:
      return FillFlowTrackGeneric(track);
    case kVZERO:
      return FillFlowTrackGeneric(track);
    case kBetaVZERO:
      return FillFlowTrackGeneric(track);
    case kDeltaVZERO:
      return FillFlowTrackGeneric(track);
    case kKappaVZERO:
      return FillFlowTrackGeneric(track);
    case kHotfixHI:
      return FillFlowTrackGeneric(track);
    default:
      return FillFlowTrackVParticle(track);
  }
}

////-----------------------------------------------------------------------
//AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackSPDtracklet() const
//{
//  //make a flow track from tracklet
//  AliFlowTrack* flowtrack=NULL;
//  TParticle *tmpTParticle=NULL;
//  AliMCParticle* tmpAliMCParticle=NULL;
//  switch (fParamMix)
//  {
//    case kPure:
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      break;
//    case kTrackWithMCkine:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi( fMCparticle->Phi() );
//      flowtrack->SetEta( fMCparticle->Eta() );
//      flowtrack->SetPt( fMCparticle->Pt() );
//      break;
//    case kTrackWithMCpt:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      flowtrack->SetPt(fMCparticle->Pt());
//      break;
//    case kTrackWithPtFromFirstMother:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      tmpTParticle = fMCparticle->Particle();
//      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
//      flowtrack->SetPt(tmpAliMCParticle->Pt());
//      break;
//    default:
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      break;
//  }
//  flowtrack->SetSource(AliFlowTrack::kFromTracklet);
//  return flowtrack;
//}
//
////-----------------------------------------------------------------------
//AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackVParticle() const
//{
//  //make flow track from AliVParticle (ESD,AOD,MC)
//  if (!fTrack) return NULL;
//  AliFlowTrack* flowtrack=NULL;
//  TParticle *tmpTParticle=NULL;
//  AliMCParticle* tmpAliMCParticle=NULL;
//  AliExternalTrackParam* externalParams=NULL;
//  AliESDtrack* esdtrack=NULL;
//  switch(fParamMix)
//  {
//    case kPure:
//      flowtrack = new AliFlowTrack(fTrack);
//      break;
//    case kTrackWithMCkine:
//      flowtrack = new AliFlowTrack(fMCparticle);
//      break;
//    case kTrackWithMCPID:
//      flowtrack = new AliFlowTrack(fTrack);
//      //flowtrack->setPID(...) from mc, when implemented
//      break;
//    case kTrackWithMCpt:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack(fTrack);
//      flowtrack->SetPt(fMCparticle->Pt());
//      break;
//    case kTrackWithPtFromFirstMother:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack(fTrack);
//      tmpTParticle = fMCparticle->Particle();
//      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
//      flowtrack->SetPt(tmpAliMCParticle->Pt());
//      break;
//    case kTrackWithTPCInnerParams:
//      esdtrack = dynamic_cast<AliESDtrack*>(fTrack);
//      if (!esdtrack) return NULL;
//      externalParams = const_cast<AliExternalTrackParam*>(esdtrack->GetTPCInnerParam());
//      if (!externalParams) return NULL;
//      flowtrack = new AliFlowTrack(externalParams);
//      break;
//    default:
//      flowtrack = new AliFlowTrack(fTrack);
//      break;
//  }
//  if (fParamType==kMC) 
//  {
//    flowtrack->SetSource(AliFlowTrack::kFromMC);
//    flowtrack->SetID(fTrackLabel);
//  }
//  else if (dynamic_cast<AliESDtrack*>(fTrack))
//  {
//    flowtrack->SetSource(AliFlowTrack::kFromESD);
//    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
//  }
//  else if (dynamic_cast<AliAODTrack*>(fTrack)) 
//  {
//    flowtrack->SetSource(AliFlowTrack::kFromAOD);
//    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
//  }
//  else if (dynamic_cast<AliMCParticle*>(fTrack)) 
//  {
//    flowtrack->SetSource(AliFlowTrack::kFromMC);
//    flowtrack->SetID(static_cast<AliVTrack*>(fTrack)->GetID());
//  }
//  return flowtrack;
//}
//
////-----------------------------------------------------------------------
//AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackPMDtrack() const
//{
//  //make a flow track from PMD track
//  AliFlowTrack* flowtrack=NULL;
//  TParticle *tmpTParticle=NULL;
//  AliMCParticle* tmpAliMCParticle=NULL;
//  switch (fParamMix)
//  {
//    case kPure:
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      break;
//    case kTrackWithMCkine:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi( fMCparticle->Phi() );
//      flowtrack->SetEta( fMCparticle->Eta() );
//      flowtrack->SetWeight(fTrackWeight);
//      flowtrack->SetPt( fMCparticle->Pt() );
//      break;
//    case kTrackWithMCpt:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      flowtrack->SetPt(fMCparticle->Pt());
//      break;
//    case kTrackWithPtFromFirstMother:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      tmpTParticle = fMCparticle->Particle();
//      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
//      flowtrack->SetPt(tmpAliMCParticle->Pt());
//      break;
//    default:
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      break;
//  }
//
//  flowtrack->SetSource(AliFlowTrack::kFromPMD);
//  return flowtrack;
//}
//
////-----------------------------------------------------------------------
//AliFlowTrack* AliFlowTrackCuts::MakeFlowTrackVZERO() const
//{
//  //make a flow track from VZERO
//  AliFlowTrack* flowtrack=NULL;
//  TParticle *tmpTParticle=NULL;
//  AliMCParticle* tmpAliMCParticle=NULL;
//  switch (fParamMix)
//  {
//    case kPure:
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      break;
//    case kTrackWithMCkine:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi( fMCparticle->Phi() );
//      flowtrack->SetEta( fMCparticle->Eta() );
//      flowtrack->SetWeight(fTrackWeight);
//      flowtrack->SetPt( fMCparticle->Pt() );
//      break;
//    case kTrackWithMCpt:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      flowtrack->SetPt(fMCparticle->Pt());
//      break;
//    case kTrackWithPtFromFirstMother:
//      if (!fMCparticle) return NULL;
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      tmpTParticle = fMCparticle->Particle();
//      tmpAliMCParticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(tmpTParticle->GetFirstMother()));
//      flowtrack->SetPt(tmpAliMCParticle->Pt());
//      break;
//    default:
//      flowtrack = new AliFlowTrack();
//      flowtrack->SetPhi(fTrackPhi);
//      flowtrack->SetEta(fTrackEta);
//      flowtrack->SetWeight(fTrackWeight);
//      break;
//  }
//
//  flowtrack->SetSource(AliFlowTrack::kFromVZERO);
//  return flowtrack;
//}
//
////-----------------------------------------------------------------------
//AliFlowTrack* AliFlowTrackCuts::MakeFlowTrack() const
//{
//  //get a flow track constructed from whatever we applied cuts on
//  //caller is resposible for deletion
//  //if construction fails return NULL
//  //TODO: for tracklets, PMD and VZERO we probably need just one method,
//  //something like MakeFlowTrackGeneric(), wait with this until
//  //requirements quirks are known.
//  switch (fParamType)
//  {
//    case kSPDtracklet:
//      return MakeFlowTrackSPDtracklet();
//    case kPMD:
//      return MakeFlowTrackPMDtrack();
//    case kVZERO:
//      return MakeFlowTrackVZERO();
//    default:
//      return MakeFlowTrackVParticle();
//  }
//}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsPhysicalPrimary() const
{
  //check if current particle is a physical primary
  if (!fMCevent) return kFALSE;
  if (fTrackLabel<0) return kFALSE;
  return IsPhysicalPrimary(fMCevent, fTrackLabel, fRequireTransportBitForPrimaries);
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::IsPhysicalPrimary(AliMCEvent* mcEvent, Int_t label, Bool_t requiretransported)
{
  //check if current particle is a physical primary
  Bool_t physprim=mcEvent->IsPhysicalPrimary(label);
  AliMCParticle* track = static_cast<AliMCParticle*>(mcEvent->GetTrack(label));
  if (!track) return kFALSE;
  TParticle* particle = track->Particle();
  Bool_t transported = particle->TestBit(kTransportBit);
  //printf("label: %i prim: %s, transp: %s, pass: %s\n",label, (physprim)?"YES":"NO ",(transported)?"YES":"NO ",
        //(physprim && (transported || !requiretransported))?"YES":"NO"  );
  return (physprim && (transported || !requiretransported));
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::DefineHistograms()
{
  //define qa histograms
  if (fQA) return;
  
  const Int_t kNbinsP=200;
  Double_t binsP[kNbinsP+1];
  binsP[0]=0.0;
  for(int i=1; i<kNbinsP+1; i++)
  {
    //if(binsP[i-1]+0.05<1.01)
    //  binsP[i]=binsP[i-1]+0.05;
    //else
    binsP[i]=binsP[i-1]+0.05;
  }

  const Int_t nBinsDCA=1000;
  Double_t binsDCA[nBinsDCA+1];
  for(int i=0; i<nBinsDCA+1; i++) {binsDCA[i]=0.01*i-5.;}
  //for(int i=1; i<41; i++) {binsDCA[i+100]=0.1*i+1.0;}

  Bool_t adddirstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fQA=new TList(); fQA->SetOwner();
  fQA->SetName(Form("%s QA",GetName()));
  TList* before = new TList(); before->SetOwner();
  before->SetName("before");
  TList* after = new TList(); after->SetOwner();
  after->SetName("after");
  fQA->Add(before);
  fQA->Add(after);
  before->Add(new TH2F("TOFbeta",";p [GeV/c];#beta",kNbinsP,binsP,1000,0.4,1.1)); //0
  after->Add(new TH2F("TOFbeta",";p [GeV/c];#beta",kNbinsP,binsP,1000,0.4,1.1)); //0
  before->Add(new TH2F("TPCdedx",";p [GeV/c];dEdx",kNbinsP,binsP,500,0,500)); //1
  after->Add(new TH2F("TPCdedx",";p [GeV/c];dEdx",kNbinsP,binsP,500,0,500)); //1
  before->Add(new TH2F("MC pid",";p[GeV/c];species",kNbinsP,binsP,10,-5, 5)); //2
  after->Add(new TH2F("MC pid",";p[GeV/c];species",kNbinsP,binsP,10,-5, 5)); //2
  //primary
  TH2F* hb = new TH2F("MC primary",";p[GeV/c];primary",kNbinsP,binsP,2,-1,1);
  TH2F* ha = new TH2F("MC primary",";p[GeV/c];primary",kNbinsP,binsP,2,-1,1);
  TAxis* axis = NULL;
  axis = hb->GetYaxis(); 
  axis->SetBinLabel(1,"secondary");
  axis->SetBinLabel(2,"primary");
  axis = ha->GetYaxis(); 
  axis->SetBinLabel(1,"secondary");
  axis->SetBinLabel(2,"primary");
  before->Add(hb); //3
  after->Add(ha); //3
  //production process
  hb = new TH2F("MC production process",";p[GeV/c];",kNbinsP,binsP,kMaxMCProcess,
                      -0.5, kMaxMCProcess-0.5);
  ha = new TH2F("MC production process",";p[GeV/c];",kNbinsP,binsP,kMaxMCProcess,
                      -0.5, kMaxMCProcess-0.5);
  axis = hb->GetYaxis();
  for (Int_t i=0; i<kMaxMCProcess; i++)
  {
    axis->SetBinLabel(i+1,TMCProcessName[i]);
  }
  axis = ha->GetYaxis();
  for (Int_t i=0; i<kMaxMCProcess; i++)
  {
    axis->SetBinLabel(i+1,TMCProcessName[i]);
  }
  before->Add(hb); //4
  after->Add(ha); //4
  //DCA
  before->Add(new TH2F("DCAxy",";p_{t}[GeV/c];DCAxy[cm]", 100, 0., 10., nBinsDCA, binsDCA));//5
  after->Add(new TH2F("DCAxy",";p_{t}[GeV/c];DCAxy[cm]", 100, 0., 10., nBinsDCA, binsDCA));//5
  before->Add(new TH2F("DCAz",";p_{t}[GeV/c];DCAz[cm]", 100, 0., 10., nBinsDCA, binsDCA));//6
  after->Add(new TH2F("DCAz",";p_{t}[GeV/c];DCAz[cm]", 100, 0., 10., nBinsDCA, binsDCA));//6
  //first mother
  hb = new TH2F("MC first mother",";p[GeV/c];",kNbinsP,binsP,44,1.,45.);
  ha = new TH2F("MC first mother",";p[GeV/c];",kNbinsP,binsP,44,1.,45.);
  hb->GetYaxis()->SetBinLabel(1, "#gamma");
  ha->GetYaxis()->SetBinLabel(1, "#gamma");
  hb->GetYaxis()->SetBinLabel(2, "e^{+}");
  ha->GetYaxis()->SetBinLabel(2, "e^{+}");
  hb->GetYaxis()->SetBinLabel(3, "e^{-}");
  ha->GetYaxis()->SetBinLabel(3, "e^{-}");
  hb->GetYaxis()->SetBinLabel(4, "#nu");
  ha->GetYaxis()->SetBinLabel(4, "#nu");
  hb->GetYaxis()->SetBinLabel(5, "#mu^{+}");
  ha->GetYaxis()->SetBinLabel(5, "#mu^{+}");
  hb->GetYaxis()->SetBinLabel(6, "#mu^{-}");
  ha->GetYaxis()->SetBinLabel(6, "#mu^{-}");
  hb->GetYaxis()->SetBinLabel(7, "#pi^{0}");
  ha->GetYaxis()->SetBinLabel(7, "#pi^{0}");
  hb->GetYaxis()->SetBinLabel(8, "#pi^{+}");
  ha->GetYaxis()->SetBinLabel(8, "#pi^{+}");
  hb->GetYaxis()->SetBinLabel(9, "#pi^{-}");
  ha->GetYaxis()->SetBinLabel(9, "#pi^{-}");
  hb->GetYaxis()->SetBinLabel(10, "K^{0}_{L}");
  ha->GetYaxis()->SetBinLabel(10, "K^{0}_{L}");
  hb->GetYaxis()->SetBinLabel(11, "K^{+}");
  ha->GetYaxis()->SetBinLabel(11, "K^{+}");
  hb->GetYaxis()->SetBinLabel(12, "K^{-}");
  ha->GetYaxis()->SetBinLabel(12, "K^{-}");
  hb->GetYaxis()->SetBinLabel( 13, "n");
  ha->GetYaxis()->SetBinLabel( 13, "n");
  hb->GetYaxis()->SetBinLabel( 14, "p");
  ha->GetYaxis()->SetBinLabel( 14, "p");
  hb->GetYaxis()->SetBinLabel(15, "#bar{p}");
  ha->GetYaxis()->SetBinLabel(15, "#bar{p}");
  hb->GetYaxis()->SetBinLabel(16, "K^{0}_{S}");
  ha->GetYaxis()->SetBinLabel(16, "K^{0}_{S}");
  hb->GetYaxis()->SetBinLabel(17, "#eta");
  ha->GetYaxis()->SetBinLabel(17, "#eta");
  hb->GetYaxis()->SetBinLabel(18, "#Lambda");
  ha->GetYaxis()->SetBinLabel(18, "#Lambda");
  hb->GetYaxis()->SetBinLabel(19, "#Sigma^{+}");
  ha->GetYaxis()->SetBinLabel(19, "#Sigma^{+}");
  hb->GetYaxis()->SetBinLabel(20, "#Sigma^{0}");
  ha->GetYaxis()->SetBinLabel(20, "#Sigma^{0}");
  hb->GetYaxis()->SetBinLabel(21, "#Sigma^{-}");
  ha->GetYaxis()->SetBinLabel(21, "#Sigma^{-}");
  hb->GetYaxis()->SetBinLabel(22, "#Theta^{0}");
  ha->GetYaxis()->SetBinLabel(22, "#Theta^{0}");
  hb->GetYaxis()->SetBinLabel(23, "#Theta^{-}");
  ha->GetYaxis()->SetBinLabel(23, "#Theta^{-}");
  hb->GetYaxis()->SetBinLabel(24, "#Omega^{-}");
  ha->GetYaxis()->SetBinLabel(24, "#Omega^{-}");
  hb->GetYaxis()->SetBinLabel(25, "#bar{n}");
  ha->GetYaxis()->SetBinLabel(25, "#bar{n}");
  hb->GetYaxis()->SetBinLabel(26, "#bar{#Lambda}");
  ha->GetYaxis()->SetBinLabel(26, "#bar{#Lambda}");
  hb->GetYaxis()->SetBinLabel(27, "#bar{#Sigma}^{-}");
  ha->GetYaxis()->SetBinLabel(27, "#bar{#Sigma}^{-}");
  hb->GetYaxis()->SetBinLabel(28, "#bar{#Sigma}^{0}");
  ha->GetYaxis()->SetBinLabel(28, "#bar{#Sigma}^{0}");
  hb->GetYaxis()->SetBinLabel(29, "#bar{#Sigma}^{+}");
  ha->GetYaxis()->SetBinLabel(29, "#bar{#Sigma}^{+}");
  hb->GetYaxis()->SetBinLabel(30, "#bar{#Theta}^{0}");
  ha->GetYaxis()->SetBinLabel(30, "#bar{#Theta}^{0}");
  hb->GetYaxis()->SetBinLabel(31, "#bar{#Theta}^{+}");
  ha->GetYaxis()->SetBinLabel(31, "#bar{#Theta}^{+}");
  hb->GetYaxis()->SetBinLabel(32, "#bar{#Omega}^{+}");
  ha->GetYaxis()->SetBinLabel(32, "#bar{#Omega}^{+}");
  hb->GetYaxis()->SetBinLabel(33, "#tau^{+}");
  ha->GetYaxis()->SetBinLabel(33, "#tau^{+}");
  hb->GetYaxis()->SetBinLabel(34, "#tau^{-}");
  ha->GetYaxis()->SetBinLabel(34, "#tau^{-}");
  hb->GetYaxis()->SetBinLabel(35, "D^{+}");
  ha->GetYaxis()->SetBinLabel(35, "D^{+}");
  hb->GetYaxis()->SetBinLabel(36, "D^{-}");
  ha->GetYaxis()->SetBinLabel(36, "D^{-}");
  hb->GetYaxis()->SetBinLabel(37, "D^{0}");
  ha->GetYaxis()->SetBinLabel(37, "D^{0}");
  hb->GetYaxis()->SetBinLabel(38, "#bar{D}^{0}");
  ha->GetYaxis()->SetBinLabel(38, "#bar{D}^{0}");
  hb->GetYaxis()->SetBinLabel(39, "D_{s}^{+}");
  ha->GetYaxis()->SetBinLabel(39, "D_{s}^{+}");
  hb->GetYaxis()->SetBinLabel(40, "#bar{D_{s}}^{-}");
  ha->GetYaxis()->SetBinLabel(40, "#bar{D_{s}}^{-}");
  hb->GetYaxis()->SetBinLabel(41, "#Lambda_{c}^{+}");
  ha->GetYaxis()->SetBinLabel(41, "#Lambda_{c}^{+}");
  hb->GetYaxis()->SetBinLabel(42, "W^{+}");
  ha->GetYaxis()->SetBinLabel(42, "W^{+}");
  hb->GetYaxis()->SetBinLabel(43, "W^{-}");
  ha->GetYaxis()->SetBinLabel(43, "W^{-}");
  hb->GetYaxis()->SetBinLabel(44, "Z^{0}");
  ha->GetYaxis()->SetBinLabel(44, "Z^{0}");
  before->Add(hb);//7
  after->Add(ha);//7

  before->Add(new TH2F("TOFkElectron",";p_{t}[GeV/c];TOF signal - IT", kNbinsP,binsP,1000,-2e4, 2e4));//8
  after->Add( new TH2F("TOFkElectron",";p_{t}[GeV/c];TOF signal - IT", kNbinsP,binsP,1000,-2e4, 2e4));//8

  before->Add(new TH2F("TOFkMuon",";p_{t}[GeV/c];TOF signal - IT",     kNbinsP,binsP,1000,-2e4, 2e4));//9
  after->Add( new TH2F("TOFkMuon",";p_{t}[GeV/c];TOF signal - IT",     kNbinsP,binsP,1000,-2e4, 2e4));//9

  before->Add(new TH2F("TOFkPion",";p_{t}[GeV/c];TOF signal - IT",     kNbinsP,binsP,1000,-2e4, 2e4));//10
  after->Add( new TH2F("TOFkPion",";p_{t}[GeV/c];TOF signal - IT",     kNbinsP,binsP,1000,-2e4, 2e4));//10

  before->Add(new TH2F("TOFkKaon",";p_{t}[GeV/c];TOF signal - IT",     kNbinsP,binsP,1000,-2e4, 2e4));//11
  after->Add( new TH2F("TOFkKaon",";p_{t}[GeV/c];TOF signal - IT",     kNbinsP,binsP,1000,-2e4, 2e4));//11

  before->Add(new TH2F("TOFkProton",";p_{t}[GeV/c];TOF signal - IT",   kNbinsP,binsP,1000,-2e4, 2e4));//12
  after->Add( new TH2F("TOFkProton",";p_{t}[GeV/c];TOF signal - IT",   kNbinsP,binsP,1000,-2e4, 2e4));//12

  //kink stuff
  before->Add(new TH1F("KinkQt",";q_{t}[GeV/c];counts", 200, 0., 0.3));//13
  after->Add(new TH1F("KinkQt",";q_{t}[GeV/c];counts", 200, 0., 0.3));//13

  before->Add(new TH1F("KinkMinv",";m_{inv}(#mu#nu)[GeV/c^{2}];counts;Kink M_{inv}", 200, 0., 0.7));//14
  after->Add(new TH1F("KinkMinv",";m_{inv}(#mu#nu)[GeV/c^{2}];counts;Kink M_{inv}", 200, 0., 0.7));//14

  before->Add(new TH2F("KinkVertex",";x[cm];y[cm];Kink vertex position",250,-250.,250., 250,-250.,250.));//15
  after->Add(new TH2F("KinkVertex",";x[cm];y[cm];Kink vertex position",250,-250.,250., 250,-250.,250.));//15

  before->Add(new TH2F("KinkAngleMp",";p_{mother}[GeV/c];Kink decay angle [deg];", 100,0.,6., 100,0.,80.));//16
  after->Add(new TH2F("KinkAngleMp",";p_{mother}[GeV/c];Kink decay angle [deg];", 100,0.,6., 100,0.,80.));//16

  before->Add(new TH1F("KinkIndex",";Kink index;counts", 2, 0., 2.));//17
  after->Add(new TH1F("KinkIndex",";Kink index;counts", 2, 0., 2.));//17
  
  before->Add(new TH2F("Chi2PerClusterTPC",";p[GeV/c];species",kNbinsP,binsP,50,0.,5.)); //18
  after->Add(new TH2F("Chi2PerClusterTPC",";p[GeV/c];species",kNbinsP,binsP,50,0.,5.)); //18
  
  before->Add(new TH2F("FractionSharedClusterTPC",";p[GeV/c];species",kNbinsP,binsP,50,0.,1.)); //19
  after->Add(new TH2F("FractionSharedClusterTPC",";p[GeV/c];species",kNbinsP,binsP,50,0.,1.)); //19
  
  before->Add(new TH2F("NClustersTPC",";p[GeV/c];species",kNbinsP,binsP,150,50.,200.)); //20
  after->Add(new TH2F("NClustersTPC",";p[GeV/c];species",kNbinsP,binsP,150,50.,200.)); //20
  
  before->Add(new TH2F("Chi2PerClusterITS",";p[GeV/c];species",kNbinsP,binsP,100,0.,50.)); //21
  after->Add(new TH2F("Chi2PerClusterITS",";p[GeV/c];species",kNbinsP,binsP,100,0.,50.)); //21

  TH1::AddDirectory(adddirstatus);
}

//-----------------------------------------------------------------------
Int_t AliFlowTrackCuts::GetNumberOfInputObjects() const
{
  //get the number of tracks in the input event according source
  //selection (ESD tracks, tracklets, MC particles etc.)
  AliESDEvent* esd=NULL;
  AliAODEvent* aod=NULL;  // XZhang 20120615
  switch (fParamType)
  {
    case kSPDtracklet:
      if (!fEvent) return 0;                                           // XZhang 20120615
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      aod = dynamic_cast<AliAODEvent*>(fEvent);                        // XZhang 20120615
//    if (!esd) return 0;                                              // XZhang 20120615
//    return esd->GetMultiplicity()->GetNumberOfTracklets();           // XZhang 20120615
      if (esd) return esd->GetMultiplicity()->GetNumberOfTracklets();  // XZhang 20120615
      if (aod) return aod->GetTracklets()->GetNumberOfTracklets();     // XZhang 20120615
    case kMC:
      if (!fMCevent) return 0;
      return fMCevent->GetNumberOfTracks();
    case kPMD:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return 0;
      return esd->GetNumberOfPmdTracks();
    case kVZERO:
      return fgkNumberOfVZEROtracks;
    case kBetaVZERO:
      return fgkNumberOfVZEROtracks;
    case kDeltaVZERO:
      return fgkNumberOfVZEROtracks;
    case kKappaVZERO:
      return fgkNumberOfVZEROtracks;
    case kHotfixHI:
      return fgkNumberOfVZEROtracks;
    case kMUON:                                      // XZhang 20120604
      if (!fEvent) return 0;                         // XZhang 20120604
      esd = dynamic_cast<AliESDEvent*>(fEvent);      // XZhang 20120604
      if (esd) return esd->GetNumberOfMuonTracks();  // XZhang 20120604
      return fEvent->GetNumberOfTracks();  // if AOD // XZhang 20120604
    case kKink:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return 0;
      return esd->GetNumberOfKinks();
    case kV0:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return 0;
      return esd->GetNumberOfV0s();
    default:
      if (!fEvent) return 0;
      return fEvent->GetNumberOfTracks();
  }
  return 0;
}

//-----------------------------------------------------------------------
TObject* AliFlowTrackCuts::GetInputObject(Int_t i)
{
  //get the input object according the data source selection:
  //(esd tracks, traclets, mc particles,etc...)
  AliESDEvent* esd=NULL;
  AliAODEvent* aod=NULL;  // XZhang 20120615
  switch (fParamType)
  {
    case kSPDtracklet:
      if (!fEvent) return NULL;                                              // XZhang 20120615
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      aod = dynamic_cast<AliAODEvent*>(fEvent);                              // XZhang 20120615
//    if (!esd) return NULL;                                                 // XZhang 20120615
//    return const_cast<AliMultiplicity*>(esd->GetMultiplicity());           // XZhang 20120615
      if (esd) return const_cast<AliMultiplicity*>(esd->GetMultiplicity());  // XZhang 20120615
      if (aod) return const_cast<AliAODTracklets*>(aod->GetTracklets());     // XZhang 20120615
    case kMC:
      if (!fMCevent) return NULL;
      return fMCevent->GetTrack(i);
    case kPMD:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return esd->GetPmdTrack(i);
    case kVZERO:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) //contributed by G.Ortona
      {
        aod = dynamic_cast<AliAODEvent*>(fEvent);
        if(!aod)return NULL;
        return aod->GetVZEROData();
      }
      return esd->GetVZEROData();
    case kBetaVZERO:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) //contributed by G.Ortona
      {
        aod = dynamic_cast<AliAODEvent*>(fEvent);
        if(!aod)return NULL;
        return aod->GetVZEROData();
      }
      return esd->GetVZEROData();
   case kDeltaVZERO:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) //contributed by G.Ortona
      {
        aod = dynamic_cast<AliAODEvent*>(fEvent);
        if(!aod)return NULL;
        return aod->GetVZEROData();
      }
      return esd->GetVZEROData();
   case kKappaVZERO:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) //contributed by G.Ortona
      {
        aod = dynamic_cast<AliAODEvent*>(fEvent);
        if(!aod)return NULL;
        return aod->GetVZEROData();
      }
      return esd->GetVZEROData();
   case kHotfixHI:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) //contributed by G.Ortona
      {
        aod = dynamic_cast<AliAODEvent*>(fEvent);
        if(!aod)return NULL;
        return aod->GetVZEROData();
      }
      return esd->GetVZEROData();
    case kMUON:                                  // XZhang 20120604
      if (!fEvent) return NULL;                  // XZhang 20120604
      esd = dynamic_cast<AliESDEvent*>(fEvent);  // XZhang 20120604
      if (esd) return esd->GetMuonTrack(i);      // XZhang 20120604
      return fEvent->GetTrack(i);  // if AOD     // XZhang 20120604
    case kKink:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return esd->GetKink(i);
    case kV0:
      esd = dynamic_cast<AliESDEvent*>(fEvent);
      if (!esd) return NULL;
      return esd->GetV0(i);
    default:
      if (!fEvent) return NULL;
      return fEvent->GetTrack(i);
  }
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::Clear(Option_t*)
{
  //clean up
  fMCevent=NULL;
  fEvent=NULL;
  ClearTrack();
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::ClearTrack(Option_t*)
{
  //clean up last track
  fKink=NULL;
  fV0=NULL;
  fTrack=NULL;
  fMCparticle=NULL;
  fTrackLabel=-997;
  fTrackWeight=1.0;
  fTrackEta=0.0;
  fTrackPhi=0.0;
  fTrackPt=0.0;
  fPOItype=1;
  fTrackMass=0.;
}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesAODpidCut(const AliAODTrack* track )
{
     if(!track->GetAODEvent()->GetTOFHeader()){
          AliAODPid *pidObj = track->GetDetPid();
          if (!pidObj) fESDpid.GetTOFResponse().SetTimeResolution(84.);
          else{
            Double_t sigmaTOFPidInAOD[10];
            pidObj->GetTOFpidResolution(sigmaTOFPidInAOD);
            if(sigmaTOFPidInAOD[0] > 84.){
              fESDpid.GetTOFResponse().SetTimeResolution(sigmaTOFPidInAOD[0]); // use the electron TOF PID sigma as time resolution (including the T0 used)
          }
        }
     }

 //check if passes the selected pid cut for ESDs
  Bool_t pass = kTRUE;
  switch (fPIDsource)
  {
   case kTOFbeta:
      if (!PassesTOFbetaCut(track)) pass=kFALSE;
      break;
  case kTOFbayesian:
      if (!PassesTOFbayesianCut(track)) pass=kFALSE;
      break;
  case kTPCbayesian:
      if (!PassesTPCbayesianCut(track)) pass=kFALSE;
      break;
  case kTPCTOFNsigma:
      if (!PassesTPCTOFNsigmaCut(track)) pass = kFALSE;
      break;
  case kTPCTOFNsigmaPurity:
      if(!PassesTPCTOFNsigmaPurityCut(track)) pass = kFALSE;
      break;
  case kTPCTPCTOFNsigma:
	  if(!PassesTPCTPCTOFNsigmaCut(track)) pass = kFALSE;
	  break;
  default:
    return kTRUE;
    break;
 }
  return pass;

}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesESDpidCut(const AliESDtrack* track )
{
  //check if passes the selected pid cut for ESDs
  Bool_t pass = kTRUE; 
  switch (fPIDsource)    
  {
    case kTPCpid:
      if (!PassesTPCpidCut(track)) pass=kFALSE;
      break;
    case kTPCdedx:
      if (!PassesTPCdedxCut(track)) pass=kFALSE;
      break;
    case kTOFpid:
      if (!PassesTOFpidCut(track)) pass=kFALSE;
      break;
    case kTOFbeta:
      if (!PassesTOFbetaCut(track)) pass=kFALSE;
      break;
    case kTOFbetaSimple:
      if (!PassesTOFbetaSimpleCut(track)) pass=kFALSE;
      break;
    case kTPCbayesian:
      if (!PassesTPCbayesianCut(track)) pass=kFALSE;
      break;
      // part added by F. Noferini
    case kTOFbayesian:
      if (!PassesTOFbayesianCut(track)) pass=kFALSE;
      break;
      // end part added by F. Noferini

      //part added by Natasha
    case kTPCNuclei:
      if (!PassesNucleiSelection(track)) pass=kFALSE;
      break;
      //end part added by Natasha
      
    case kTPCTOFNsigma:
      if (!PassesTPCTOFNsigmaCut(track)) pass = kFALSE;
      break;
    default:
      printf("AliFlowTrackCuts::PassesCuts() this should never be called!\n");
      pass=kFALSE;
      break;
  }
  return pass;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbetaSimpleCut(const AliESDtrack* track )
{
  //check if passes PID cut using timing in TOF
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFpid) &&
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365);
                    
  if (!fAllowTOFmismatchFlag) {if ((track->GetStatus() & AliESDtrack::kTOFmismatch)) return kFALSE;}

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;

  if (!goodtrack) return kFALSE;
  
  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->GetP();
  Float_t l = track->GetIntegratedLength();  
  Float_t trackT0 = fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t timeTOF = track->GetTOFsignal()- trackT0; 
  Float_t beta = l/timeTOF/c;
  Double_t integratedTimes[9] = {-1.0,-1.0,-1.0,-1.0,-1.0, -1.0, -1.0, -1.0, -1.0};
  track->GetIntegratedTimes(integratedTimes);
  Float_t betaHypothesis[5] = {0.0,0.0,0.0,0.0,0.0};
  Float_t s[5] = {0.0,0.0,0.0,0.0,0.0};
  for (Int_t i=0;i<5;i++)
  {
    betaHypothesis[i] = l/integratedTimes[i]/c;
    s[i] = beta-betaHypothesis[i];
  }

  switch (fParticleID)
  {
    case AliPID::kPion:
      return ( (s[2]<0.015) && (s[2]>-0.015) &&
               (s[3]>0.025) &&
               (s[4]>0.03) );
    case AliPID::kKaon:
      return ( (s[3]<0.015) && (s[3]>-0.015) &&
               (s[2]<-0.03) &&
               (s[4]>0.03) );
    case AliPID::kProton:
      return ( (s[4]<0.015) && (s[4]>-0.015) &&
               (s[3]<-0.025) &&
               (s[2]<-0.025) );
    default:
      return kFALSE;
  }
  return kFALSE;
}

//-----------------------------------------------------------------------
Float_t AliFlowTrackCuts::GetBeta(const AliVTrack* track, Bool_t QAmode)
{
  //get beta
  Double_t integratedTimes[9] = {-1.0,-1.0,-1.0,-1.0,-1.0, -1.0, -1.0, -1.0, -1.0};
  track->GetIntegratedTimes(integratedTimes);

  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->P();
  Float_t l = integratedTimes[0]*c;  
  Float_t trackT0 = fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t timeTOF = track->GetTOFsignal()- trackT0; 
  if(QAmode && timeTOF <= 0) return -999;   // avoid division by zero when filling 'before' qa histograms
  return l/timeTOF/c;
}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbetaCut(const AliAODTrack* track )
{
  //check if track passes pid selection with an asymmetric TOF beta cut
  if (!fTOFpidCuts)
  {
    //printf("no TOFpidCuts\n");
    return kFALSE;
  }

  //check if passes PID cut using timing in TOF
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFpid) &&
                     (track->GetTOFsignal() > 12000) &&
                     (track->GetTOFsignal() < 100000);

  if (!goodtrack) return kFALSE;

  const Float_t c = 2.99792457999999984e-02;
  Double_t integratedTimes[9] = {-1.0,-1.0,-1.0,-1.0,-1.0, -1.0, -1.0, -1.0, -1.0};
  track->GetIntegratedTimes(integratedTimes);
  Float_t l = integratedTimes[0]*c;

  goodtrack = goodtrack && (l > 365);

  if (!goodtrack) return kFALSE;

  if (!fAllowTOFmismatchFlag) {if ((track->GetStatus() & AliESDtrack::kTOFmismatch)) return kFALSE;}

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;


  Float_t beta = GetBeta(track);

  //construct the pid index because it's not AliPID::EParticleType
  Int_t pid = 0;
    cout<<"TOFbeta: fParticleID = "<<fParticleID<<endl;
  switch (fParticleID)
  {
    case AliPID::kPion:
      pid=2;
      break;
    case AliPID::kKaon:
      pid=3;
      break;
    case AliPID::kProton:
      pid=4;
      break;
    default:
      return kFALSE;
  }

  //signal to cut on
  Float_t p = track->P();
  Float_t betahypothesis = l/integratedTimes[pid]/c;
  Float_t betadiff = beta-betahypothesis;

  Float_t* arr = fTOFpidCuts->GetMatrixArray();
  Int_t col = TMath::BinarySearch(fTOFpidCuts->GetNcols(),arr,static_cast<Float_t>(p));
  if (col<0) return kFALSE;
  Float_t min = (*fTOFpidCuts)(1,col);
  Float_t max = (*fTOFpidCuts)(2,col);

  Bool_t pass = (betadiff>min && betadiff<max);

  return pass;
}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbetaCut(const AliESDtrack* track )
{
  //check if track passes pid selection with an asymmetric TOF beta cut
  if (!fTOFpidCuts)
  {
    //printf("no TOFpidCuts\n");
    return kFALSE;
  }

  //check if passes PID cut using timing in TOF
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFpid) && 
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365);

  if (!fAllowTOFmismatchFlag) {if ((track->GetStatus() & AliESDtrack::kTOFmismatch)) return kFALSE;}

  if (!goodtrack) return kFALSE;

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;
  
  Float_t beta = GetBeta(track);
  Double_t integratedTimes[9] = {-1.0,-1.0,-1.0,-1.0,-1.0, -1.0, -1.0, -1.0, -1.0};
  track->GetIntegratedTimes(integratedTimes);

  //construct the pid index because it's not AliPID::EParticleType
  Int_t pid = 0;
  switch (fParticleID)
  {
    case AliPID::kPion:
      pid=2;
      break;
    case AliPID::kKaon:
      pid=3;
      break;
    case AliPID::kProton:
      pid=4;
      break;
    default:
      return kFALSE;
  }

  //signal to cut on
  const Float_t c = 2.99792457999999984e-02;  
  Float_t l = track->GetIntegratedLength();  
  Float_t p = track->GetP();  
  Float_t betahypothesis = l/integratedTimes[pid]/c;
  Float_t betadiff = beta-betahypothesis;

  Float_t* arr = fTOFpidCuts->GetMatrixArray();
  Int_t col = TMath::BinarySearch(fTOFpidCuts->GetNcols(),arr,static_cast<Float_t>(p));
  if (col<0) return kFALSE;
  Float_t min = (*fTOFpidCuts)(1,col);
  Float_t max = (*fTOFpidCuts)(2,col);

  Bool_t pass = (betadiff>min && betadiff<max);
  
  return pass;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFpidCut(const AliESDtrack* track) const
{
  //check if passes PID cut using default TOF pid
  Double_t pidTOF[AliPID::kSPECIES];
  track->GetTOFpid(pidTOF);
  if (pidTOF[fParticleID]>=fParticleProbability) return kTRUE;
  return kFALSE;
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCpidCut(const AliESDtrack* track) const
{
  //check if passes PID cut using default TPC pid
  Double_t pidTPC[AliPID::kSPECIES];
  track->GetTPCpid(pidTPC);
  Double_t probablity = 0.;
  switch (fParticleID)
  {
    case AliPID::kPion:
      probablity = pidTPC[AliPID::kPion] + pidTPC[AliPID::kMuon];
      break;
    default:
      probablity = pidTPC[fParticleID];
  }
  if (probablity >= fParticleProbability) return kTRUE;
  return kFALSE;
}

//-----------------------------------------------------------------------
Float_t AliFlowTrackCuts::Getdedx(const AliESDtrack* track) const
{
  //get TPC dedx
  return track->GetTPCsignal();
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCdedxCut(const AliESDtrack* track)
{
  //check if passes PID cut using dedx signal in the TPC
  if (!fTPCpidCuts)
  {
    //printf("no TPCpidCuts\n");
    return kFALSE;
  }

  const AliExternalTrackParam* tpcparam = track->GetInnerParam(); //tpc only params at the inner wall
  if (!tpcparam) return kFALSE;
  Double_t p = tpcparam->GetP();
  Float_t sigExp = fESDpid.GetTPCResponse().GetExpectedSignal(p, fParticleID);
  Float_t sigTPC = track->GetTPCsignal();
  Float_t s = (sigTPC-sigExp)/sigExp;

  Float_t* arr = fTPCpidCuts->GetMatrixArray();
  Int_t arrSize = fTPCpidCuts->GetNcols();
  Int_t col = TMath::BinarySearch( arrSize, arr, static_cast<Float_t>(p));
  if (col<0) return kFALSE;
  Float_t min = (*fTPCpidCuts)(1,col);
  Float_t max = (*fTPCpidCuts)(2,col);

  //printf("------------TPC pid cut %s\n",(s>min && s<max)?"PASS":"FAIL");
  return (s>min && s<max);
}

//-----------------------------------------------------------------------
void AliFlowTrackCuts::InitPIDcuts()
{
  //init matrices with PID cuts
  TMatrixF* t = NULL;
  if (!fTPCpidCuts)
  {
    if (fParticleID==AliPID::kPion)
    {
      t = new TMatrixF(3,15);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.4;  (*t)(2,0)  =   0.0;
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.4;  (*t)(2,1)  =   0.1;
      (*t)(0,2)  = 0.30;  (*t)(1,2)  = -0.4;  (*t)(2,2)  =  0.2;
      (*t)(0,3)  = 0.35;  (*t)(1,3)  = -0.4;  (*t)(2,3)  =  0.2;
      (*t)(0,4)  = 0.40;  (*t)(1,4)  = -0.4;  (*t)(2,4)  =   0.3;
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.4;  (*t)(2,5)  =   0.3;
      (*t)(0,6)  = 0.50;  (*t)(1,6)  = -0.4;  (*t)(2,6)  =  0.25;
      (*t)(0,7)  = 0.55;  (*t)(1,7)  = -0.4;  (*t)(2,7)  =  0.15;
      (*t)(0,8)  = 0.60;  (*t)(1,8)  = -0.4;  (*t)(2,8)  =   0.1;
      (*t)(0,9)  = 0.65;  (*t)(1,9)  = -0.4;  (*t)(2,9)  =  0.05;
      (*t)(0,10)  = 0.70;  (*t)(1,10)  = -0.4;  (*t)(2,10)  =     0;
      (*t)(0,11)  = 0.75;  (*t)(1,11)  = -0.4;  (*t)(2,11)  =     0;
      (*t)(0,12)  = 0.80;  (*t)(1,12)  = -0.4;  (*t)(2,12)  = -0.05;
      (*t)(0,13)  = 0.85;  (*t)(1,13)  = -0.4;  (*t)(2,13)  = -0.1;
      (*t)(0,14)  = 0.90;  (*t)(1,14)  = 0;     (*t)(2,14)  =     0;
    }
    else
    if (fParticleID==AliPID::kKaon)
    {
      t = new TMatrixF(3,12);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.2;  (*t)(2,0)  = 0.2; 
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.2;  (*t)(2,1)  = 0.2;
      (*t)(0,2)  = 0.30;  (*t)(1,2)  = -0.2;  (*t)(2,2)  = 0.2;
      (*t)(0,3)  = 0.35;  (*t)(1,3)  = -0.2;  (*t)(2,3)  = 0.2;
      (*t)(0,4)  = 0.40;  (*t)(1,4)  = -0.1;  (*t)(2,4)  = 0.2;
      (*t)(0,5)  = 0.45;  (*t)(1,5)  = -0.1;  (*t)(2,5)  = 0.2;
      (*t)(0,6)  = 0.50;  (*t)(1,6)  =-0.05;  (*t)(2,6)  = 0.2;
      (*t)(0,7)  = 0.55;  (*t)(1,7)  = -0.1;  (*t)(2,7)  = 0.1;
      (*t)(0,8)  = 0.60;  (*t)(1,8)  =-0.05;  (*t)(2,8)  = 0.1;
      (*t)(0,9)  = 0.65;  (*t)(1,9)  =    0;  (*t)(2,9)  = 0.15;
      (*t)(0,10)  = 0.70;  (*t)(1,10)  = 0.05;  (*t)(2,10)  = 0.2;
      (*t)(0,11)  = 0.75;  (*t)(1,11)  =    0;  (*t)(2,11)  = 0;
    }
    else
    if (fParticleID==AliPID::kProton)
    {
      t = new TMatrixF(3,9);
      (*t)(0,0)  = 0.20;  (*t)(1,0)  = -0.1;  (*t)(2,0)  =  0.1; 
      (*t)(0,1)  = 0.25;  (*t)(1,1)  = -0.2;  (*t)(2,1)  =  0.2; 
      (*t)(0,2)  = 0.80;  (*t)(1,2)  = -0.1;  (*t)(2,2)  =  0.2; 
      (*t)(0,3)  = 0.85;  (*t)(1,3)  =-0.05;  (*t)(2,3)  =  0.2; 
      (*t)(0,4)  = 0.90;  (*t)(1,4)  =-0.05;  (*t)(2,4)  = 0.25; 
      (*t)(0,5)  = 0.95;  (*t)(1,5)  =-0.05;  (*t)(2,5)  = 0.25; 
      (*t)(0,6)  = 1.00;  (*t)(1,6)  = -0.1;  (*t)(2,6)  = 0.25; 
      (*t)(0,7)  = 1.10;  (*t)(1,7)  =-0.05;  (*t)(2,7)  =  0.3; 
      (*t)(0,8) = 1.20;   (*t)(1,8)  =    0;  (*t)(2,8) =    0;
    }
    delete fTPCpidCuts;
    fTPCpidCuts=t;
  }
  t = NULL;
  if (!fTOFpidCuts)
  {
    if (fParticleID==AliPID::kPion)
    {
      //TOF pions, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(1,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(1,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(1,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(1,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(1,4)  = -0.030;  (*t)(2,4)  =   0.030;
      (*t)(0,5)  = 0.250;  (*t)(1,5)  = -0.036;  (*t)(2,5)  =   0.032;
      (*t)(0,6)  = 0.300;  (*t)(1,6)  = -0.038;  (*t)(2,6)  =   0.032;
      (*t)(0,7)  = 0.350;  (*t)(1,7)  = -0.034;  (*t)(2,7)  =   0.032;
      (*t)(0,8)  = 0.400;  (*t)(1,8)  = -0.032;  (*t)(2,8)  =   0.020;
      (*t)(0,9)  = 0.450;  (*t)(1,9)  = -0.030;  (*t)(2,9)  =   0.020;
      (*t)(0,10)  = 0.500;  (*t)(1,10)  = -0.030;  (*t)(2,10)  =   0.020;
      (*t)(0,11)  = 0.550;  (*t)(1,11)  = -0.030;  (*t)(2,11)  =   0.020;
      (*t)(0,12)  = 0.600;  (*t)(1,12)  = -0.030;  (*t)(2,12)  =   0.020;
      (*t)(0,13)  = 0.650;  (*t)(1,13)  = -0.030;  (*t)(2,13)  =   0.020;
      (*t)(0,14)  = 0.700;  (*t)(1,14)  = -0.030;  (*t)(2,14)  =   0.020;
      (*t)(0,15)  = 0.750;  (*t)(1,15)  = -0.030;  (*t)(2,15)  =   0.020;
      (*t)(0,16)  = 0.800;  (*t)(1,16)  = -0.030;  (*t)(2,16)  =   0.020;
      (*t)(0,17)  = 0.850;  (*t)(1,17)  = -0.030;  (*t)(2,17)  =   0.020;
      (*t)(0,18)  = 0.900;  (*t)(1,18)  = -0.030;  (*t)(2,18)  =   0.020;
      (*t)(0,19)  = 0.950;  (*t)(1,19)  = -0.028;  (*t)(2,19)  =   0.028;
      (*t)(0,20)  = 1.000;  (*t)(1,20)  = -0.028;  (*t)(2,20)  =   0.028;
      (*t)(0,21)  = 1.100;  (*t)(1,21)  = -0.028;  (*t)(2,21)  =   0.028;
      (*t)(0,22)  = 1.200;  (*t)(1,22)  = -0.026;  (*t)(2,22)  =   0.028;
      (*t)(0,23)  = 1.300;  (*t)(1,23)  = -0.024;  (*t)(2,23)  =   0.028;
      (*t)(0,24)  = 1.400;  (*t)(1,24)  = -0.020;  (*t)(2,24)  =   0.028;
      (*t)(0,25)  = 1.500;  (*t)(1,25)  = -0.018;  (*t)(2,25)  =   0.028;
      (*t)(0,26)  = 1.600;  (*t)(1,26)  = -0.016;  (*t)(2,26)  =   0.028;
      (*t)(0,27)  = 1.700;  (*t)(1,27)  = -0.014;  (*t)(2,27)  =   0.028;
      (*t)(0,28)  = 1.800;  (*t)(1,28)  = -0.012;  (*t)(2,28)  =   0.026;
      (*t)(0,29)  = 1.900;  (*t)(1,29)  = -0.010;  (*t)(2,29)  =   0.026;
      (*t)(0,30)  = 2.000;  (*t)(1,30)  = -0.008;  (*t)(2,30)  =   0.026;
      (*t)(0,31)  = 2.100;  (*t)(1,31)  = -0.008;  (*t)(2,31)  =   0.024;
      (*t)(0,32)  = 2.200;  (*t)(1,32)  = -0.006;  (*t)(2,32)  =   0.024;
      (*t)(0,33)  = 2.300;  (*t)(1,33)  = -0.004;  (*t)(2,33)  =   0.024;
      (*t)(0,34)  = 2.400;  (*t)(1,34)  = -0.004;  (*t)(2,34)  =   0.024;
      (*t)(0,35)  = 2.500;  (*t)(1,35)  = -0.002;  (*t)(2,35)  =   0.024;
      (*t)(0,36)  = 2.600;  (*t)(1,36)  = -0.002;  (*t)(2,36)  =   0.024;
      (*t)(0,37)  = 2.700;  (*t)(1,37)  = 0.000;  (*t)(2,37)  =   0.024;
      (*t)(0,38)  = 2.800;  (*t)(1,38)  = 0.000;  (*t)(2,38)  =   0.026;
      (*t)(0,39)  = 2.900;  (*t)(1,39)  = 0.000;  (*t)(2,39)  =   0.024;
      (*t)(0,40)  = 3.000;  (*t)(1,40)  = 0.002;  (*t)(2,40)  =   0.026;
      (*t)(0,41)  = 3.100;  (*t)(1,41)  = 0.002;  (*t)(2,41)  =   0.026;
      (*t)(0,42)  = 3.200;  (*t)(1,42)  = 0.002;  (*t)(2,42)  =   0.026;
      (*t)(0,43)  = 3.300;  (*t)(1,43)  = 0.002;  (*t)(2,43)  =   0.026;
      (*t)(0,44)  = 3.400;  (*t)(1,44)  = 0.002;  (*t)(2,44)  =   0.026;
      (*t)(0,45)  = 3.500;  (*t)(1,45)  = 0.002;  (*t)(2,45)  =   0.026;
      (*t)(0,46)  = 3.600;  (*t)(1,46)  = 0.002;  (*t)(2,46)  =   0.026;
      (*t)(0,47)  = 3.700;  (*t)(1,47)  = 0.002;  (*t)(2,47)  =   0.026;
      (*t)(0,48)  = 3.800;  (*t)(1,48)  = 0.002;  (*t)(2,48)  =   0.026;
      (*t)(0,49)  = 3.900;  (*t)(1,49)  = 0.004;  (*t)(2,49)  =   0.024;
      (*t)(0,50)  = 4.000;  (*t)(1,50)  = 0.004;  (*t)(2,50)  =   0.026;
      (*t)(0,51)  = 4.100;  (*t)(1,51)  = 0.004;  (*t)(2,51)  =   0.026;
      (*t)(0,52)  = 4.200;  (*t)(1,52)  = 0.004;  (*t)(2,52)  =   0.024;
      (*t)(0,53)  = 4.300;  (*t)(1,53)  = 0.006;  (*t)(2,53)  =   0.024;
      (*t)(0,54)  = 4.400;  (*t)(1,54)  = 0.000;  (*t)(2,54)  =   0.000;
      (*t)(0,55)  = 4.500;  (*t)(1,55)  = 0.000;  (*t)(2,55)  =   0.000;
      (*t)(0,56)  = 4.600;  (*t)(1,56)  = 0.000;  (*t)(2,56)  =   0.000;
      (*t)(0,57)  = 4.700;  (*t)(1,57)  = 0.000;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(1,58)  = 0.000;  (*t)(2,58)  =   0.000;
      (*t)(0,59)  = 4.900;  (*t)(1,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(1,60)  = 0.000;  (*t)(2,60)  =   0.000;
    }
    else
    if (fParticleID==AliPID::kProton)
    {
      //TOF protons, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(1,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(1,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(1,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(1,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(1,4)  = -0.07;  (*t)(2,4)  =   0.07;
      (*t)(0,5)  = 0.200;  (*t)(1,5)  = -0.07;  (*t)(2,5)  =   0.07;
      (*t)(0,6)  = 0.200;  (*t)(1,6)  = -0.07;  (*t)(2,6)  =   0.07;
      (*t)(0,7)  = 0.200;  (*t)(1,7)  = -0.07;  (*t)(2,7)  =   0.07;
      (*t)(0,8)  = 0.200;  (*t)(1,8)  = -0.07;  (*t)(2,8)  =   0.07;
      (*t)(0,9)  = 0.200;  (*t)(1,9)  = -0.07;  (*t)(2,9)  =   0.07;
      (*t)(0,10)  = 0.200;  (*t)(1,10)  = -0.07;  (*t)(2,10)  =   0.07;
      (*t)(0,11)  = 0.200;  (*t)(1,11)  = -0.07;  (*t)(2,11)  =   0.07;
      (*t)(0,12)  = 0.200;  (*t)(1,12)  = -0.07;  (*t)(2,12)  =   0.07;
      (*t)(0,13)  = 0.200;  (*t)(1,13)  = -0.07;  (*t)(2,13)  =   0.07;
      (*t)(0,14)  = 0.200;  (*t)(1,14)  = -0.07;  (*t)(2,14)  =   0.07;
      (*t)(0,15)  = 0.200;  (*t)(1,15)  = -0.07;  (*t)(2,15)  =   0.07;
      (*t)(0,16)  = 0.200;  (*t)(1,16)  = -0.07;  (*t)(2,16)  =   0.07;
      (*t)(0,17)  = 0.850;  (*t)(1,17)  = -0.070;  (*t)(2,17)  =   0.070;
      (*t)(0,18)  = 0.900;  (*t)(1,18)  = -0.072;  (*t)(2,18)  =   0.072;
      (*t)(0,19)  = 0.950;  (*t)(1,19)  = -0.072;  (*t)(2,19)  =   0.072;
      (*t)(0,20)  = 1.000;  (*t)(1,20)  = -0.074;  (*t)(2,20)  =   0.074;
      (*t)(0,21)  = 1.100;  (*t)(1,21)  = -0.032;  (*t)(2,21)  =   0.032;
      (*t)(0,22)  = 1.200;  (*t)(1,22)  = -0.026;  (*t)(2,22)  =   0.026;
      (*t)(0,23)  = 1.300;  (*t)(1,23)  = -0.026;  (*t)(2,23)  =   0.026;
      (*t)(0,24)  = 1.400;  (*t)(1,24)  = -0.024;  (*t)(2,24)  =   0.024;
      (*t)(0,25)  = 1.500;  (*t)(1,25)  = -0.024;  (*t)(2,25)  =   0.024;
      (*t)(0,26)  = 1.600;  (*t)(1,26)  = -0.026;  (*t)(2,26)  =   0.026;
      (*t)(0,27)  = 1.700;  (*t)(1,27)  = -0.026;  (*t)(2,27)  =   0.026;
      (*t)(0,28)  = 1.800;  (*t)(1,28)  = -0.026;  (*t)(2,28)  =   0.026;
      (*t)(0,29)  = 1.900;  (*t)(1,29)  = -0.026;  (*t)(2,29)  =   0.026;
      (*t)(0,30)  = 2.000;  (*t)(1,30)  = -0.026;  (*t)(2,30)  =   0.026;
      (*t)(0,31)  = 2.100;  (*t)(1,31)  = -0.026;  (*t)(2,31)  =   0.026;
      (*t)(0,32)  = 2.200;  (*t)(1,32)  = -0.026;  (*t)(2,32)  =   0.024;
      (*t)(0,33)  = 2.300;  (*t)(1,33)  = -0.028;  (*t)(2,33)  =   0.022;
      (*t)(0,34)  = 2.400;  (*t)(1,34)  = -0.028;  (*t)(2,34)  =   0.020;
      (*t)(0,35)  = 2.500;  (*t)(1,35)  = -0.028;  (*t)(2,35)  =   0.018;
      (*t)(0,36)  = 2.600;  (*t)(1,36)  = -0.028;  (*t)(2,36)  =   0.016;
      (*t)(0,37)  = 2.700;  (*t)(1,37)  = -0.028;  (*t)(2,37)  =   0.016;
      (*t)(0,38)  = 2.800;  (*t)(1,38)  = -0.030;  (*t)(2,38)  =   0.014;
      (*t)(0,39)  = 2.900;  (*t)(1,39)  = -0.030;  (*t)(2,39)  =   0.012;
      (*t)(0,40)  = 3.000;  (*t)(1,40)  = -0.030;  (*t)(2,40)  =   0.012;
      (*t)(0,41)  = 3.100;  (*t)(1,41)  = -0.030;  (*t)(2,41)  =   0.010;
      (*t)(0,42)  = 3.200;  (*t)(1,42)  = -0.030;  (*t)(2,42)  =   0.010;
      (*t)(0,43)  = 3.300;  (*t)(1,43)  = -0.030;  (*t)(2,43)  =   0.010;
      (*t)(0,44)  = 3.400;  (*t)(1,44)  = -0.030;  (*t)(2,44)  =   0.008;
      (*t)(0,45)  = 3.500;  (*t)(1,45)  = -0.030;  (*t)(2,45)  =   0.008;
      (*t)(0,46)  = 3.600;  (*t)(1,46)  = -0.030;  (*t)(2,46)  =   0.008;
      (*t)(0,47)  = 3.700;  (*t)(1,47)  = -0.030;  (*t)(2,47)  =   0.006;
      (*t)(0,48)  = 3.800;  (*t)(1,48)  = -0.030;  (*t)(2,48)  =   0.006;
      (*t)(0,49)  = 3.900;  (*t)(1,49)  = -0.030;  (*t)(2,49)  =   0.006;
      (*t)(0,50)  = 4.000;  (*t)(1,50)  = -0.028;  (*t)(2,50)  =   0.004;
      (*t)(0,51)  = 4.100;  (*t)(1,51)  = -0.030;  (*t)(2,51)  =   0.004;
      (*t)(0,52)  = 4.200;  (*t)(1,52)  = -0.030;  (*t)(2,52)  =   0.004;
      (*t)(0,53)  = 4.300;  (*t)(1,53)  = -0.028;  (*t)(2,53)  =   0.002;
      (*t)(0,54)  = 4.400;  (*t)(1,54)  = -0.030;  (*t)(2,54)  =   0.002;
      (*t)(0,55)  = 4.500;  (*t)(1,55)  = -0.028;  (*t)(2,55)  =   0.002;
      (*t)(0,56)  = 4.600;  (*t)(1,56)  = -0.028;  (*t)(2,56)  =   0.002;
      (*t)(0,57)  = 4.700;  (*t)(1,57)  = -0.028;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(1,58)  = -0.028;  (*t)(2,58)  =   0.002;
      (*t)(0,59)  = 4.900;  (*t)(1,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(1,60)  = 0.000;  (*t)(2,60)  =   0.000; 
    }
    else
    if (fParticleID==AliPID::kKaon)
    {
      //TOF kaons, 0.9 purity
      t = new TMatrixF(3,61);
      (*t)(0,0)  = 0.000;  (*t)(1,0)  = 0.000;  (*t)(2,0)  =   0.000;
      (*t)(0,1)  = 0.050;  (*t)(1,1)  = 0.000;  (*t)(2,1)  =   0.000;
      (*t)(0,2)  = 0.100;  (*t)(1,2)  = 0.000;  (*t)(2,2)  =   0.000;
      (*t)(0,3)  = 0.150;  (*t)(1,3)  = 0.000;  (*t)(2,3)  =   0.000;
      (*t)(0,4)  = 0.200;  (*t)(1,4)  = -0.05;  (*t)(2,4)  =   0.05;
      (*t)(0,5)  = 0.200;  (*t)(1,5)  = -0.05;  (*t)(2,5)  =   0.05;
      (*t)(0,6)  = 0.200;  (*t)(1,6)  = -0.05;  (*t)(2,6)  =   0.05;
      (*t)(0,7)  = 0.200;  (*t)(1,7)  = -0.05;  (*t)(2,7)  =   0.05;
      (*t)(0,8)  = 0.200;  (*t)(1,8)  = -0.05;  (*t)(2,8)  =   0.05;
      (*t)(0,9)  = 0.200;  (*t)(1,9)  = -0.05;  (*t)(2,9)  =   0.05;
      (*t)(0,10)  = 0.200;  (*t)(1,10)  = -0.05;  (*t)(2,10)  =   0.05;
      (*t)(0,11)  = 0.550;  (*t)(1,11)  = -0.026;  (*t)(2,11)  =   0.026;
      (*t)(0,12)  = 0.600;  (*t)(1,12)  = -0.026;  (*t)(2,12)  =   0.026;
      (*t)(0,13)  = 0.650;  (*t)(1,13)  = -0.026;  (*t)(2,13)  =   0.026;
      (*t)(0,14)  = 0.700;  (*t)(1,14)  = -0.026;  (*t)(2,14)  =   0.026;
      (*t)(0,15)  = 0.750;  (*t)(1,15)  = -0.026;  (*t)(2,15)  =   0.026;
      (*t)(0,16)  = 0.800;  (*t)(1,16)  = -0.026;  (*t)(2,16)  =   0.026;
      (*t)(0,17)  = 0.850;  (*t)(1,17)  = -0.024;  (*t)(2,17)  =   0.024;
      (*t)(0,18)  = 0.900;  (*t)(1,18)  = -0.024;  (*t)(2,18)  =   0.024;
      (*t)(0,19)  = 0.950;  (*t)(1,19)  = -0.024;  (*t)(2,19)  =   0.024;
      (*t)(0,20)  = 1.000;  (*t)(1,20)  = -0.024;  (*t)(2,20)  =   0.024;
      (*t)(0,21)  = 1.100;  (*t)(1,21)  = -0.024;  (*t)(2,21)  =   0.024;
      (*t)(0,22)  = 1.200;  (*t)(1,22)  = -0.024;  (*t)(2,22)  =   0.022;
      (*t)(0,23)  = 1.300;  (*t)(1,23)  = -0.024;  (*t)(2,23)  =   0.020;
      (*t)(0,24)  = 1.400;  (*t)(1,24)  = -0.026;  (*t)(2,24)  =   0.016;
      (*t)(0,25)  = 1.500;  (*t)(1,25)  = -0.028;  (*t)(2,25)  =   0.014;
      (*t)(0,26)  = 1.600;  (*t)(1,26)  = -0.028;  (*t)(2,26)  =   0.012;
      (*t)(0,27)  = 1.700;  (*t)(1,27)  = -0.028;  (*t)(2,27)  =   0.010;
      (*t)(0,28)  = 1.800;  (*t)(1,28)  = -0.028;  (*t)(2,28)  =   0.010;
      (*t)(0,29)  = 1.900;  (*t)(1,29)  = -0.028;  (*t)(2,29)  =   0.008;
      (*t)(0,30)  = 2.000;  (*t)(1,30)  = -0.028;  (*t)(2,30)  =   0.006;
      (*t)(0,31)  = 2.100;  (*t)(1,31)  = -0.026;  (*t)(2,31)  =   0.006;
      (*t)(0,32)  = 2.200;  (*t)(1,32)  = -0.024;  (*t)(2,32)  =   0.004;
      (*t)(0,33)  = 2.300;  (*t)(1,33)  = -0.020;  (*t)(2,33)  =   0.002;
      (*t)(0,34)  = 2.400;  (*t)(1,34)  = -0.020;  (*t)(2,34)  =   0.002;
      (*t)(0,35)  = 2.500;  (*t)(1,35)  = -0.018;  (*t)(2,35)  =   0.000;
      (*t)(0,36)  = 2.600;  (*t)(1,36)  = -0.016;  (*t)(2,36)  =   0.000;
      (*t)(0,37)  = 2.700;  (*t)(1,37)  = -0.014;  (*t)(2,37)  =   -0.002;
      (*t)(0,38)  = 2.800;  (*t)(1,38)  = -0.014;  (*t)(2,38)  =   -0.004;
      (*t)(0,39)  = 2.900;  (*t)(1,39)  = -0.012;  (*t)(2,39)  =   -0.004;
      (*t)(0,40)  = 3.000;  (*t)(1,40)  = -0.010;  (*t)(2,40)  =   -0.006;
      (*t)(0,41)  = 3.100;  (*t)(1,41)  = 0.000;  (*t)(2,41)  =   0.000;
      (*t)(0,42)  = 3.200;  (*t)(1,42)  = 0.000;  (*t)(2,42)  =   0.000;
      (*t)(0,43)  = 3.300;  (*t)(1,43)  = 0.000;  (*t)(2,43)  =   0.000;
      (*t)(0,44)  = 3.400;  (*t)(1,44)  = 0.000;  (*t)(2,44)  =   0.000;
      (*t)(0,45)  = 3.500;  (*t)(1,45)  = 0.000;  (*t)(2,45)  =   0.000;
      (*t)(0,46)  = 3.600;  (*t)(1,46)  = 0.000;  (*t)(2,46)  =   0.000;
      (*t)(0,47)  = 3.700;  (*t)(1,47)  = 0.000;  (*t)(2,47)  =   0.000;
      (*t)(0,48)  = 3.800;  (*t)(1,48)  = 0.000;  (*t)(2,48)  =   0.000;
      (*t)(0,49)  = 3.900;  (*t)(1,49)  = 0.000;  (*t)(2,49)  =   0.000;
      (*t)(0,50)  = 4.000;  (*t)(1,50)  = 0.000;  (*t)(2,50)  =   0.000;
      (*t)(0,51)  = 4.100;  (*t)(1,51)  = 0.000;  (*t)(2,51)  =   0.000;
      (*t)(0,52)  = 4.200;  (*t)(1,52)  = 0.000;  (*t)(2,52)  =   0.000;
      (*t)(0,53)  = 4.300;  (*t)(1,53)  = 0.000;  (*t)(2,53)  =   0.000;
      (*t)(0,54)  = 4.400;  (*t)(1,54)  = 0.000;  (*t)(2,54)  =   0.000;
      (*t)(0,55)  = 4.500;  (*t)(1,55)  = 0.000;  (*t)(2,55)  =   0.000;
      (*t)(0,56)  = 4.600;  (*t)(1,56)  = 0.000;  (*t)(2,56)  =   0.000;
      (*t)(0,57)  = 4.700;  (*t)(1,57)  = 0.000;  (*t)(2,57)  =   0.000;
      (*t)(0,58)  = 4.800;  (*t)(1,58)  = 0.000;  (*t)(2,58)  =   0.000;
      (*t)(0,59)  = 4.900;  (*t)(1,59)  = 0.000;  (*t)(2,59)  =   0.000;
      (*t)(0,60)  = 5.900;  (*t)(1,60)  = 0.000;  (*t)(2,60)  =   0.000;
    }
    delete fTOFpidCuts;
    fTOFpidCuts=t;
  }
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCbayesianCut(const AliAODTrack* track)
{
  fBayesianResponse->ComputeProb(track,track->GetAODEvent()); // fCurrCentr is needed for mismatch fraction
  Float_t *probabilities = fBayesianResponse->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mism$

  Int_t kTPC = fBayesianResponse->GetCurrentMask(0); // is TPC on

  if(! kTPC) return kFALSE;

  fProbBayes = 0.0;

switch (fParticleID)
  {
    case AliPID::kPion:
      fProbBayes = probabilities[2];
      break;
    case AliPID::kKaon:
       fProbBayes = probabilities[3];
     break;
    case AliPID::kProton:
       fProbBayes = probabilities[4];
      break;
    case AliPID::kElectron:
       fProbBayes = probabilities[0];
     break;
    case AliPID::kMuon:
       fProbBayes = probabilities[1];
     break;
    case AliPID::kDeuteron:
       fProbBayes = probabilities[5];
     break;
    case AliPID::kTriton:
       fProbBayes = probabilities[6];
     break;
    case AliPID::kHe3:
       fProbBayes = probabilities[7];
     break;
   default:
      return kFALSE;
  }

 if(fProbBayes > fParticleProbability){
    if(!fCutCharge)
      return kTRUE;
    else if (fCutCharge && fCharge * track->Charge() > 0)
      return kTRUE;
  }
  return kFALSE;
}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCbayesianCut(const AliESDtrack* track)
{
  //cut on TPC bayesian pid
  //if (! fAllowTOFmismatchFlag) {if (track->GetStatus() & AliESDtrack::kTOFmismatch) return kFALSE;}

  //Bool_t statusMatchingHard = TPCTOFagree(track);
  //if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
  //     return kFALSE;
  fBayesianResponse->ComputeProb(track,fCurrCentr); // fCurrCentr is needed for mismatch fraction
  Int_t kTPC = fBayesianResponse->GetCurrentMask(0); // is TPC on
  //Int_t kTOF = fBayesianResponse->GetCurrentMask(1); // is TOF on

  if(! kTPC) return kFALSE;

  //  Bool_t statusMatchingHard = 1;
  //  Float_t mismProb = 0;
  //  if(kTOF){
  //    statusMatchingHard = TPCTOFagree(track);
  //    mismProb = fBayesianResponse->GetTOFMismProb();
  //  }
  //  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
  //       return kFALSE;

  Float_t *probabilities = fBayesianResponse->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mismatch parameterization

  fProbBayes = 0.0;

  switch (fParticleID)
  {
    case AliPID::kPion:
      fProbBayes = probabilities[2];
      break;
    case AliPID::kKaon:
      fProbBayes = probabilities[3];
     break;
    case AliPID::kProton:
      fProbBayes = probabilities[4];
      break;
    case AliPID::kElectron:
       fProbBayes = probabilities[0];
     break;
    case AliPID::kMuon:
       fProbBayes = probabilities[1];
     break;
    case AliPID::kDeuteron:
       fProbBayes = probabilities[5];
     break;
    case AliPID::kTriton:
       fProbBayes = probabilities[6];
     break;
    case AliPID::kHe3:
       fProbBayes = probabilities[7];
     break;
    default:
      return kFALSE;
  }

  if(fProbBayes > fParticleProbability)
    {
      if(!fCutCharge)
	return kTRUE;
      else if (fCutCharge && fCharge * track->GetSign() > 0)
	return kTRUE;
    }
  return kFALSE;
}
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTOFbayesianCut(const AliAODTrack* track)
{  
  //check is track passes bayesian combined TOF+TPC pid cut
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFout) &&
                     (track->GetStatus() & AliESDtrack::kTIME) &&
                     (track->GetTOFsignal() > 12000) &&
                     (track->GetTOFsignal() < 100000);

  if (! goodtrack)   
       return kFALSE;

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;

 fBayesianResponse->ComputeProb(track,track->GetAODEvent()); // fCurrCentr is needed for mismatch fraction
  Float_t *probabilities = fBayesianResponse->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mism$

  Float_t mismProb = fBayesianResponse->GetTOFMismProb(); // mismatch Bayesian probabilities

  fProbBayes = 0.0;

 switch (fParticleID)
  {
    case AliPID::kPion:
      fProbBayes = probabilities[2];
      break;
    case AliPID::kKaon:
       fProbBayes = probabilities[3];
     break;
    case AliPID::kProton:
       fProbBayes = probabilities[4];
      break;
    case AliPID::kElectron:
       fProbBayes = probabilities[0];
     break;
    case AliPID::kMuon:
       fProbBayes = probabilities[1];
     break;
    case AliPID::kDeuteron:
       fProbBayes = probabilities[5];
     break;
    case AliPID::kTriton:
       fProbBayes = probabilities[6];
     break;
    case AliPID::kHe3:
       fProbBayes = probabilities[7];
     break;
   default:
      return kFALSE;
  }

 if(fProbBayes > fParticleProbability && mismProb < 0.5){
    if(!fCutCharge)
      return kTRUE;
    else if (fCutCharge && fCharge * track->Charge() > 0)
      return kTRUE;
  }
  return kFALSE;

}
//-----------------------------------------------------------------------
// part added by F. Noferini (some methods)
Bool_t AliFlowTrackCuts::PassesTOFbayesianCut(const AliESDtrack* track)
{
  //check is track passes bayesian combined TOF+TPC pid cut
  Bool_t goodtrack = (track->GetStatus() & AliESDtrack::kTOFout) && 
                     (track->GetStatus() & AliESDtrack::kTIME) &&
                     (track->GetTOFsignal() > 12000) && 
                     (track->GetTOFsignal() < 100000) && 
                     (track->GetIntegratedLength() > 365);

  if (! goodtrack)
       return kFALSE;

  if (! fAllowTOFmismatchFlag) {if (track->GetStatus() & AliESDtrack::kTOFmismatch) return kFALSE;}

  Bool_t statusMatchingHard = TPCTOFagree(track);
  if (fRequireStrictTOFTPCagreement && (!statusMatchingHard))
       return kFALSE;

  fBayesianResponse->ComputeProb(track,fCurrCentr); // fCurrCentr is needed for mismatch fraction
  Float_t *probabilities = fBayesianResponse->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mismatch parameterization

  Float_t mismProb = fBayesianResponse->GetTOFMismProb(); // mismatch Bayesian probabilities

  fProbBayes = 0.0;

  switch (fParticleID)
  {
    case AliPID::kPion:
      fProbBayes = probabilities[2];
      break;
    case AliPID::kKaon:
       fProbBayes = probabilities[3];
     break;
    case AliPID::kProton:
       fProbBayes = probabilities[4];
      break;
    case AliPID::kElectron:
       fProbBayes = probabilities[0];
     break;
    case AliPID::kMuon:
       fProbBayes = probabilities[1];
     break;
    case AliPID::kDeuteron:
       fProbBayes = probabilities[5];
     break;
    case AliPID::kTriton:
       fProbBayes = probabilities[6];
     break;
    case AliPID::kHe3:
       fProbBayes = probabilities[7];
     break;
   default:
      return kFALSE;
  }

  //  printf("pt = %f -- all prob = [%4.2f,%4.2f,%4.2f,%4.2f,%4.2f] -- prob = %f\n",track->Pt(),fProbBayes[0],fProbBayes[1],fProbBayes[2],fProbBayes[3],fProbBayes[4],prob);
  if(fProbBayes > fParticleProbability && mismProb < 0.5){
    if(!fCutCharge)
      return kTRUE;
    else if (fCutCharge && fCharge * track->GetSign() > 0)
      return kTRUE;
  }
  return kFALSE;
}


//-----------------------------------------------------------------------
 // part added by Natasha
Bool_t AliFlowTrackCuts::PassesNucleiSelection(const AliESDtrack* track)
{
  //pid selection for heavy nuclei
  Bool_t select=kFALSE;

  //if (!track) continue; 

  if (!track->GetInnerParam()) 
    return kFALSE;    //break;

  const AliExternalTrackParam* tpcTrack = track->GetInnerParam();

  Double_t ptotTPC = tpcTrack->GetP();
  Double_t sigTPC = track->GetTPCsignal();
  Double_t dEdxBBA = 0.;
  Double_t dSigma = 0.; 

  switch (fParticleID)
  {
    case AliPID::kDeuteron:
      //pid=10;
      dEdxBBA = AliExternalTrackParam::BetheBlochAleph(ptotTPC/1.8756,
          4.60e+00,
          8.9684e+00,
          1.640e-05,
          2.35e+00,
          2.35e+00);
      dSigma = (sigTPC -  dEdxBBA)/dEdxBBA;

      if( ptotTPC<=1.1 && (dSigma < (0.5 - (0.1818*ptotTPC)) ) && (dSigma > ( (0.218*ptotTPC - 0.4) ) )  )
      {select=kTRUE;}
      break;

    case AliPID::kTriton:
      //pid=11;
      select=kFALSE;
      break;

    case AliPID::kHe3:
      //pid=12;
      // ----- Pass 2 -------
      dEdxBBA = 4.0 * AliExternalTrackParam::BetheBlochAleph( (2.*ptotTPC)/2.8084,
          1.74962,
          27.4992,
          4.00313e-15,
          2.42485,
          8.31768);
      dSigma = (sigTPC -  dEdxBBA)/dEdxBBA;
      if(ptotTPC<=5.0 && (dSigma >= (-0.03968*ptotTPC - 0.1)) && (dSigma <= (0.31 - 0.0217*ptotTPC)))
      {select=kTRUE;}
      break;

    case AliPID::kAlpha:
      //pid=13;
      select=kFALSE;
      break;

    default:
      return kFALSE;
  }       

  return select;
}
// end part added by Natasha
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCTOFNsigmaCut(const AliAODTrack* track) 
{
    // do a simple combined cut on the n sigma from tpc and tof
    // with information of the pid response object (needs pid response task)
    // stub, not implemented yet
    if(!fPIDResponse) return kFALSE;
    if(!track) return kFALSE;

    // check TOF status
    if ((track->GetStatus()&AliVTrack::kTOFout)==0) return kFALSE;
    if ((track->GetStatus()&AliVTrack::kTIME)==0) return kFALSE;

    // check TPC status
    if(track->GetTPCsignal() < 10) return kFALSE;

    Float_t nsigmaTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,fParticleID);
    Float_t nsigmaTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,fParticleID);

    Float_t nsigma2 = nsigmaTPC*nsigmaTPC + nsigmaTOF*nsigmaTOF;

    return (nsigma2 < fNsigmaCut2);

}
//-----------------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCTOFNsigmaCut(const AliESDtrack* track)
{
    // do a simple combined cut on the n sigma from tpc and tof
    // with information of the pid response object (needs pid response task)
    // stub, not implemented yet
    if(!fPIDResponse) return kFALSE;
    if(!track) return kFALSE;

    // check TOF status
    if ((track->GetStatus()&AliVTrack::kTOFout)==0) return kFALSE;
    if ((track->GetStatus()&AliVTrack::kTIME)==0) return kFALSE;

    // check TPC status
    if(track->GetTPCsignal() < 10) return kFALSE;

    Float_t nsigmaTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,fParticleID);
    Float_t nsigmaTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,fParticleID);

    Float_t nsigma2 = nsigmaTPC*nsigmaTPC + nsigmaTOF*nsigmaTOF;

    return (nsigma2 < fNsigmaCut2);

}
//-----------------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCTOFNsigmaPurityCut(const AliAODTrack* track)
{
    // do a combined cut on the n sigma from tpc and tof based on purity of the identified particle. (Particle is selected if Purity>fPurityLevel & nsigma<3)
    // with information of the pid response object (needs pid response task)
    
    if(!fPIDResponse) return kFALSE;
    if(!track) return kFALSE;
    
    Int_t p_int = -999;
    Double_t pInterval=0;
    for(int i=0;i<60;i++){
        pInterval=0.1*i;
        if(track->P()>pInterval && track->P()<pInterval+0.1){p_int = i;}
    }
    /*
     Double_t LowPtPIDTPCnsigLow_Pion[2] = {-3,-3};
     Double_t LowPtPIDTPCnsigLow_Kaon[2] = {-3,-2};
     Double_t LowPtPIDTPCnsigHigh_Pion[2] ={2.4,3};
     Double_t LowPtPIDTPCnsigHigh_Kaon[2] ={3,2.2};
     */
    
    Float_t nsigmaTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,fParticleID);
    Float_t nsigmaTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,fParticleID);
    
    int index = (fParticleID-2)*60 + p_int;
    if ( (track->IsOn(AliAODTrack::kITSin))){
        if(p_int<2) return kFALSE;
        
        if(!fPurityFunction[index]){ cout<<"fPurityFunction[index] does not exist"<<endl; return kFALSE;}
        if(p_int>1){
            if((track->IsOn(AliAODTrack::kTOFpid))){
                if(fPurityFunction[index]->Eval(nsigmaTOF,nsigmaTPC)>fPurityLevel){
                    if(TMath::Sqrt(TMath::Power(nsigmaTPC,2)+TMath::Power(nsigmaTOF,2))<3){
                        return kTRUE;
                    }
                }
            }
        }
        /*if(p_int<4){
         if(fParticleID==AliPID::kKaon && nsigmaTPC>LowPtPIDTPCnsigLow_Kaon[p_int-2] && nsigmaTPC<LowPtPIDTPCnsigHigh_Kaon[p_int-2]){return kTRUE;}
         if(fParticleID==AliPID::kPion && nsigmaTPC>LowPtPIDTPCnsigLow_Pion[p_int-2] && nsigmaTPC<LowPtPIDTPCnsigHigh_Pion[p_int-2]){return kTRUE;}
         if(fParticleID==AliPID::kProton && nsigmaTPC>-3 && nsigmaTPC<3){return kTRUE;}
         }*/
    }
    return kFALSE;
}
//-----------------------------------------------------------------------
void AliFlowTrackCuts::SetPriors(Float_t centrCur){
 //set priors for the bayesian pid selection
  fCurrCentr = centrCur;

  fBinLimitPID[0] = 0.300000;
  fBinLimitPID[1] = 0.400000;
  fBinLimitPID[2] = 0.500000;
  fBinLimitPID[3] = 0.600000;
  fBinLimitPID[4] = 0.700000;
  fBinLimitPID[5] = 0.800000;
  fBinLimitPID[6] = 0.900000;
  fBinLimitPID[7] = 1.000000;
  fBinLimitPID[8] = 1.200000;
  fBinLimitPID[9] = 1.400000;
  fBinLimitPID[10] = 1.600000;
  fBinLimitPID[11] = 1.800000;
  fBinLimitPID[12] = 2.000000;
  fBinLimitPID[13] = 2.200000;
  fBinLimitPID[14] = 2.400000;
  fBinLimitPID[15] = 2.600000;
  fBinLimitPID[16] = 2.800000;
  fBinLimitPID[17] = 3.000000;
 
  // 0-10%
  if(centrCur < 10){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0168;
      fC[1][4] = 0.0122;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0272;
      fC[2][4] = 0.0070;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0562;
      fC[3][4] = 0.0258;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0861;
      fC[4][4] = 0.0496;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1168;
      fC[5][4] = 0.0740;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1476;
      fC[6][4] = 0.0998;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1810;
      fC[7][4] = 0.1296;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2240;
      fC[8][4] = 0.1827;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2812;
      fC[9][4] = 0.2699;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3328;
      fC[10][4] = 0.3714;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3780;
      fC[11][4] = 0.4810;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.4125;
      fC[12][4] = 0.5771;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4486;
      fC[13][4] = 0.6799;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4840;
      fC[14][4] = 0.7668;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4971;
      fC[15][4] = 0.8288;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4956;
      fC[16][4] = 0.8653;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.5173;
      fC[17][4] = 0.9059;   
  }
  // 10-20%
  else if(centrCur < 20){
     fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0132;
      fC[1][4] = 0.0088;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0283;
      fC[2][4] = 0.0068;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0577;
      fC[3][4] = 0.0279;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0884;
      fC[4][4] = 0.0534;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1179;
      fC[5][4] = 0.0794;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1480;
      fC[6][4] = 0.1058;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1807;
      fC[7][4] = 0.1366;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2219;
      fC[8][4] = 0.1891;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2804;
      fC[9][4] = 0.2730;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3283;
      fC[10][4] = 0.3660;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3710;
      fC[11][4] = 0.4647;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.4093;
      fC[12][4] = 0.5566;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4302;
      fC[13][4] = 0.6410;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4649;
      fC[14][4] = 0.7055;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4523;
      fC[15][4] = 0.7440;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4591;
      fC[16][4] = 0.7799;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4804;
      fC[17][4] = 0.8218;
  }
  // 20-30%
  else if(centrCur < 30){
     fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0102;
      fC[1][4] = 0.0064;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0292;
      fC[2][4] = 0.0066;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0597;
      fC[3][4] = 0.0296;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0900;
      fC[4][4] = 0.0589;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1199;
      fC[5][4] = 0.0859;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1505;
      fC[6][4] = 0.1141;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1805;
      fC[7][4] = 0.1454;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2221;
      fC[8][4] = 0.2004;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2796;
      fC[9][4] = 0.2838;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3271;
      fC[10][4] = 0.3682;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3648;
      fC[11][4] = 0.4509;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3988;
      fC[12][4] = 0.5339;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4315;
      fC[13][4] = 0.5995;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4548;
      fC[14][4] = 0.6612;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4744;
      fC[15][4] = 0.7060;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4899;
      fC[16][4] = 0.7388;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4411;
      fC[17][4] = 0.7293;
  }
  // 30-40%
  else if(centrCur < 40){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0102;
      fC[1][4] = 0.0048;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0306;
      fC[2][4] = 0.0079;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0617;
      fC[3][4] = 0.0338;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0920;
      fC[4][4] = 0.0652;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1211;
      fC[5][4] = 0.0955;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1496;
      fC[6][4] = 0.1242;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1807;
      fC[7][4] = 0.1576;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2195;
      fC[8][4] = 0.2097;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2732;
      fC[9][4] = 0.2884;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3204;
      fC[10][4] = 0.3679;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3564;
      fC[11][4] = 0.4449;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3791;
      fC[12][4] = 0.5052;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.4062;
      fC[13][4] = 0.5647;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.4234;
      fC[14][4] = 0.6203;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4441;
      fC[15][4] = 0.6381;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4629;
      fC[16][4] = 0.6496;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4293;
      fC[17][4] = 0.6491;
  }
  // 40-50%
  else if(centrCur < 50){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0093;
      fC[1][4] = 0.0057;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0319;
      fC[2][4] = 0.0075;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0639;
      fC[3][4] = 0.0371;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0939;
      fC[4][4] = 0.0725;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1224;
      fC[5][4] = 0.1045;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1520;
      fC[6][4] = 0.1387;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1783;
      fC[7][4] = 0.1711;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2202;
      fC[8][4] = 0.2269;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2672;
      fC[9][4] = 0.2955;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3191;
      fC[10][4] = 0.3676;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3434;
      fC[11][4] = 0.4321;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3692;
      fC[12][4] = 0.4879;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3993;
      fC[13][4] = 0.5377;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3818;
      fC[14][4] = 0.5547;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4003;
      fC[15][4] = 0.5484;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4281;
      fC[16][4] = 0.5383;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3960;
      fC[17][4] = 0.5374;
  }
  // 50-60%
  else if(centrCur < 60){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0076;
      fC[1][4] = 0.0032;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0329;
      fC[2][4] = 0.0085;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0653;
      fC[3][4] = 0.0423;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0923;
      fC[4][4] = 0.0813;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1219;
      fC[5][4] = 0.1161;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1519;
      fC[6][4] = 0.1520;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1763;
      fC[7][4] = 0.1858;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2178;
      fC[8][4] = 0.2385;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2618;
      fC[9][4] = 0.3070;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.3067;
      fC[10][4] = 0.3625;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3336;
      fC[11][4] = 0.4188;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3706;
      fC[12][4] = 0.4511;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3765;
      fC[13][4] = 0.4729;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3942;
      fC[14][4] = 0.4855;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.4051;
      fC[15][4] = 0.4762;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.3843;
      fC[16][4] = 0.4763;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.4237;
      fC[17][4] = 0.4773;
  }
  // 60-70%
  else if(centrCur < 70){
         fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0071;
      fC[1][4] = 0.0012;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0336;
      fC[2][4] = 0.0097;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0662;
      fC[3][4] = 0.0460;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0954;
      fC[4][4] = 0.0902;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1181;
      fC[5][4] = 0.1306;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1481;
      fC[6][4] = 0.1662;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1765;
      fC[7][4] = 0.1963;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2155;
      fC[8][4] = 0.2433;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2580;
      fC[9][4] = 0.3022;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.2872;
      fC[10][4] = 0.3481;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3170;
      fC[11][4] = 0.3847;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3454;
      fC[12][4] = 0.4258;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3580;
      fC[13][4] = 0.4299;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3903;
      fC[14][4] = 0.4326;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.3690;
      fC[15][4] = 0.4491;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.4716;
      fC[16][4] = 0.4298;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3875;
      fC[17][4] = 0.4083;
  }
  // 70-80%
  else if(centrCur < 80){
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0075;
      fC[1][4] = 0.0007;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0313;
      fC[2][4] = 0.0124;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0640;
      fC[3][4] = 0.0539;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0923;
      fC[4][4] = 0.0992;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1202;
      fC[5][4] = 0.1417;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1413;
      fC[6][4] = 0.1729;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1705;
      fC[7][4] = 0.1999;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.2103;
      fC[8][4] = 0.2472;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2373;
      fC[9][4] = 0.2916;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.2824;
      fC[10][4] = 0.3323;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.3046;
      fC[11][4] = 0.3576;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3585;
      fC[12][4] = 0.4003;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3461;
      fC[13][4] = 0.3982;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3362;
      fC[14][4] = 0.3776;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.3071;
      fC[15][4] = 0.3500;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.2914;
      fC[16][4] = 0.3937;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3727;
      fC[17][4] = 0.3877;
  }
  // 80-100%
  else{
      fC[0][0] = 0.005;
      fC[0][1] = 0.005;
      fC[0][2] = 1.0000;
      fC[0][3] = 0.0010;
      fC[0][4] = 0.0010;

      fC[1][0] = 0.005;
      fC[1][1] = 0.005;
      fC[1][2] = 1.0000;
      fC[1][3] = 0.0060;
      fC[1][4] = 0.0035;

      fC[2][0] = 0.005;
      fC[2][1] = 0.005;
      fC[2][2] = 1.0000;
      fC[2][3] = 0.0323;
      fC[2][4] = 0.0113;

      fC[3][0] = 0.005;
      fC[3][1] = 0.005;
      fC[3][2] = 1.0000;
      fC[3][3] = 0.0609;
      fC[3][4] = 0.0653;

      fC[4][0] = 0.005;
      fC[4][1] = 0.005;
      fC[4][2] = 1.0000;
      fC[4][3] = 0.0922;
      fC[4][4] = 0.1076;

      fC[5][0] = 0.005;
      fC[5][1] = 0.005;
      fC[5][2] = 1.0000;
      fC[5][3] = 0.1096;
      fC[5][4] = 0.1328;

      fC[6][0] = 0.005;
      fC[6][1] = 0.005;
      fC[6][2] = 1.0000;
      fC[6][3] = 0.1495;
      fC[6][4] = 0.1779;

      fC[7][0] = 0.005;
      fC[7][1] = 0.005;
      fC[7][2] = 1.0000;
      fC[7][3] = 0.1519;
      fC[7][4] = 0.1989;

      fC[8][0] = 0.005;
      fC[8][1] = 0.005;
      fC[8][2] = 1.0000;
      fC[8][3] = 0.1817;
      fC[8][4] = 0.2472;

      fC[9][0] = 0.005;
      fC[9][1] = 0.005;
      fC[9][2] = 1.0000;
      fC[9][3] = 0.2429;
      fC[9][4] = 0.2684;

      fC[10][0] = 0.005;
      fC[10][1] = 0.005;
      fC[10][2] = 1.0000;
      fC[10][3] = 0.2760;
      fC[10][4] = 0.3098;

      fC[11][0] = 0.005;
      fC[11][1] = 0.005;
      fC[11][2] = 1.0000;
      fC[11][3] = 0.2673;
      fC[11][4] = 0.3198;

      fC[12][0] = 0.005;
      fC[12][1] = 0.005;
      fC[12][2] = 1.0000;
      fC[12][3] = 0.3165;
      fC[12][4] = 0.3564;

      fC[13][0] = 0.005;
      fC[13][1] = 0.005;
      fC[13][2] = 1.0000;
      fC[13][3] = 0.3526;
      fC[13][4] = 0.3011;

      fC[14][0] = 0.005;
      fC[14][1] = 0.005;
      fC[14][2] = 1.0000;
      fC[14][3] = 0.3788;
      fC[14][4] = 0.3011;

      fC[15][0] = 0.005;
      fC[15][1] = 0.005;
      fC[15][2] = 1.0000;
      fC[15][3] = 0.3788;
      fC[15][4] = 0.3011;

      fC[16][0] = 0.005;
      fC[16][1] = 0.005;
      fC[16][2] = 1.0000;
      fC[16][3] = 0.3788;
      fC[16][4] = 0.3011;

      fC[17][0] = 0.005;
      fC[17][1] = 0.005;
      fC[17][2] = 1.0000;
      fC[17][3] = 0.3788;
      fC[17][4] = 0.3011;
  }
  
  for(Int_t i=18;i<fgkPIDptBin;i++){
    fBinLimitPID[i] = 3.0 + 0.2 * (i-18);
    fC[i][0] = fC[17][0];
    fC[i][1] = fC[17][1];
    fC[i][2] = fC[17][2];
    fC[i][3] = fC[17][3];
    fC[i][4] = fC[17][4];
  }  
}

//---------------------------------------------------------------//
Bool_t AliFlowTrackCuts::TPCTOFagree(const AliVTrack *track)
{
  //check pid agreement between TPC and TOF
  Bool_t status = kFALSE;

  const Float_t c = 2.99792457999999984e-02;

  Float_t mass[5] = {5.10998909999999971e-04,1.05658000000000002e-01,1.39570000000000000e-01,4.93676999999999977e-01,9.38271999999999995e-01};
  

  Double_t exptimes[9];
  track->GetIntegratedTimes(exptimes);
  
  Float_t dedx = track->GetTPCsignal();

  Float_t p = track->P();
  Float_t time = track->GetTOFsignal()- fESDpid.GetTOFResponse().GetStartTime(p);
  Float_t tl = exptimes[0]*c; // to work both for ESD and AOD

  Float_t betagammares =  fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);

  Float_t betagamma1 = tl/(time-5 *betagammares) * 33.3564095198152043;

//  printf("betagamma1 = %f\n",betagamma1);

  if(betagamma1 < 0.1) betagamma1 = 0.1;

  if(betagamma1 < 0.99999) betagamma1 /= TMath::Sqrt(1-betagamma1*betagamma1);
  else betagamma1 = 100;

  Float_t betagamma2 = tl/(time+5 *betagammares) * 33.3564095198152043;
//  printf("betagamma2 = %f\n",betagamma2);

  if(betagamma2 < 0.1) betagamma2 = 0.1;

  if(betagamma2 < 0.99999) betagamma2 /= TMath::Sqrt(1-betagamma2*betagamma2);
  else betagamma2 = 100;


  Float_t momtpc=track->GetTPCmomentum();
 
  for(Int_t i=0;i < 5;i++){
    Float_t resolutionTOF =  fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[i], mass[i]);
    if(TMath::Abs(exptimes[i] - time) < 5 * resolutionTOF){
      Float_t dedxExp = 0;
      if(i==0) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kElectron);
      else if(i==1) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kMuon);
      else if(i==2) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kPion);
      else if(i==3) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kKaon);
      else if(i==4) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kProton);

      Float_t resolutionTPC = 2;
      if(i==0) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kElectron); 
      else if(i==1) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kMuon);
      else if(i==2) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kPion);
      else if(i==3) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kKaon);
      else if(i==4) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);

      if(TMath::Abs(dedx - dedxExp) < 3 * resolutionTPC){
	status = kTRUE;
      }
    }
  }

  Float_t bb1 =  fESDpid.GetTPCResponse().Bethe(betagamma1);
  Float_t bb2 =  fESDpid.GetTPCResponse().Bethe(betagamma2);
  Float_t bbM =  fESDpid.GetTPCResponse().Bethe((betagamma1+betagamma2)*0.5);


  //  status = kFALSE;
  // for nuclei
  Float_t resolutionTOFpr =   fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);
  Float_t resolutionTPCpr =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);
  if(TMath::Abs(dedx-bb1) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  else if(TMath::Abs(dedx-bb2) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  else if(TMath::Abs(dedx-bbM) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  
  return status;
}
// end part added by F. Noferini
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesTPCTPCTOFNsigmaCut(const AliAODTrack* track){

  //Added by Bernhard Hohlweger (bernhard.hohlweger@cern.ch)
  //This method uses only the TPC dEdx up to a certain Pt for particle identification, above this Pt value the combined information from TOF and TPC is used.
  //Default values tuned on 2010h.
  //Pt threshold can be set by calling AliFlowTrackCuts::SetPtTOFPIDoff(Double_t pt).

  Bool_t pass = kTRUE;

  if(!fPIDResponse) pass = kFALSE;
  if(!track) pass = kFALSE;

  // check TPC status
  if(track->GetTPCsignal() < 10) pass = kFALSE;

  if(fPtTOFPIDoff == 0.){
    if(fParticleID == AliPID::kProton){
      fPtTOFPIDoff = 0.75;
    }else if(fParticleID == AliPID::kPion){
      fPtTOFPIDoff = 0.55;
    }else if(fParticleID == AliPID::kKaon){
      fPtTOFPIDoff = 0.45;
    }else{
      //by default just use the standard TPCTOFNsigmaCut
      fPtTOFPIDoff = 0.;
    }
  }
  if(pass){
    Double_t Pt = track->Pt();
    Float_t nsigmaTPC = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,track,fParticleID);
    Float_t nsigma2 = 999.;
    if(Pt < fPtTOFPIDoff){
      nsigma2 = nsigmaTPC*nsigmaTPC;
    }else{
      // check TOF status
      if (((track->GetStatus()&AliVTrack::kTOFout)==0)&&((track->GetStatus()&AliVTrack::kTIME)==0)){
        pass = kFALSE;
      }else{
        Float_t nsigmaTOF = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,track,fParticleID);
        nsigma2 = nsigmaTPC*nsigmaTPC + nsigmaTOF*nsigmaTOF;
      }
    }
    pass = pass && (nsigma2 < fNsigmaCut2);
  }

  return pass;

}
//-----------------------------------------------------------------------
const char* AliFlowTrackCuts::PIDsourceName(PIDsource s)
{
  //get the name of the particle id source
  switch (s)
  {
    case kTPCdedx:
      return "TPCdedx";
    case kTPCbayesian:
      return "TPCbayesian";
    case kTOFbeta:
      return "TOFbeta";
    case kTPCpid:
      return "TPCpid";
    case kTOFpid:
      return "TOFpid";
    case kTOFbayesian:
      return "TOFbayesianPID";
    case kTOFbetaSimple:
      return "TOFbetaSimple";
    case kTPCNuclei:
      return "TPCnuclei";
    case kTPCTOFNsigma:
      return "TPCTOFNsigma";
    case kTPCTOFNsigmaPurity:
      return "TPCTOFNsigmaPurity";
    default:
      return "NOPID";
  }
}

//-----------------------------------------------------------------------
const char* AliFlowTrackCuts::GetParamTypeName(trackParameterType type) 
{
  //return the name of the selected parameter type
  switch (type)
  {
    case kMC:
      return "MC";
    case kGlobal:
      return "Global";
    case kTPCstandalone:
      return "TPCstandalone";
    case kSPDtracklet:
      return "SPDtracklets";
    case kPMD:
      return "PMD";
    case kVZERO:
      return "VZERO";
    case kBetaVZERO:
      return "BetaVZERO";
    case kDeltaVZERO:
      return "DeltaVZERO";
    case kKappaVZERO:
      return "KappaVZERO";
    case kHotfixHI:
      return "HotfixHI";
    case kMUON:       // XZhang 20120604
      return "MUON";  // XZhang 20120604
    case kKink:
      return "Kink";
    case kV0:
      return "V0";
    default:
      return "unknown";
  }
}

//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesPMDcuts(const AliESDPmdTrack* track )
{
  //check PMD specific cuts
  //clean up from last iteration, and init label
  Int_t   det   = track->GetDetector();
  //Int_t   smn   = track->GetSmn();
  Float_t clsX  = track->GetClusterX();
  Float_t clsY  = track->GetClusterY();
  Float_t clsZ  = track->GetClusterZ();
  Float_t ncell = track->GetClusterCells();
  Float_t adc   = track->GetClusterADC();

  fTrack = NULL;
  fMCparticle=NULL;
  fTrackLabel=-996;

  fTrackEta = GetPmdEta(clsX,clsY,clsZ);
  fTrackPhi = GetPmdPhi(clsX,clsY);
  fTrackWeight = 1.0;

  Bool_t pass=kTRUE;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) pass = kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) pass = kFALSE;}
  if (fCutPmdDet) {if(det != fPmdDet) pass = kFALSE;}
  if (fCutPmdAdc) {if(adc < fPmdAdc) pass = kFALSE;}
  if (fCutPmdNcell) {if(ncell < fPmdNcell) pass = kFALSE;}

  return pass;
}
  
//-----------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesVZEROcuts(Int_t id)
{
  //check VZERO cuts
  if (id<0) return kFALSE;

  //clean up from last iter
  ClearTrack();

  fTrackPhi = TMath::PiOver4()*(0.5+id%8);

  // 10102013 weighting vzero tiles - rbertens@cern.ch
  if(!fVZEROgainEqualization) {
      // if for some reason the equalization is not initialized (e.g. 2011 data)
      // the fVZEROxpol[] weights are used to enable or disable vzero rings
    if(id<32) {   // v0c side
      fTrackEta = -3.45+0.5*(id/8);
      if(id < 8) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[0];
      else if (id < 16 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[1];
      else if (id < 24 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[2];
      else if (id < 32 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[3];
    } else {       // v0a side
      fTrackEta = +4.8-0.6*((id/8)-4);
      if( id < 40) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[0];
      else if ( id < 48 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[1];
      else if ( id < 56 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[2];
      else if ( id < 64 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[3];
    }
  } else { // the equalization is initialized
     // note that disabled rings have already been excluded on calibration level in 
     // AliFlowEvent (so for a disabled ring, fVZEROxpol is zero
    if(id<32) {    // v0c side
      fTrackEta = -3.45+0.5*(id/8);
      if(id < 8) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[0]/fVZEROgainEqualization->GetBinContent(1+id);
      else if (id < 16 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[1]/fVZEROgainEqualization->GetBinContent(1+id);
      else if (id < 24 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[2]/fVZEROgainEqualization->GetBinContent(1+id);
      else if (id < 32 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[3]/fVZEROgainEqualization->GetBinContent(1+id);
    } else {       // v0a side
      fTrackEta = +4.8-0.6*((id/8)-4);
      if( id < 40) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[0]/fVZEROgainEqualization->GetBinContent(1+id);
      else if ( id < 48 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[1]/fVZEROgainEqualization->GetBinContent(1+id);
      else if ( id < 56 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[2]/fVZEROgainEqualization->GetBinContent(1+id);
      else if ( id < 64 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[3]/fVZEROgainEqualization->GetBinContent(1+id);
    }
    // printf ( " tile %i and weight %.2f \n", id, fTrackWeight);
  }
 
  // 29052015 weighting vzero tiles - jacopo.margutti@cern.ch
  if(fVZEROgainEqualizationCen) {
    //Added by Bernhard Hohlweger - bhohlweg@cern.ch
	Double_t EventCentrality = -999;
	if(fEvent->GetRunNumber() < 209122){
	  //For Run1 Data the Old Centrality Percentile Method is available whereas for Run2 a new method was implemented
	  //Cut was done for the first run of the LHC15a period
	  EventCentrality = fEvent->GetCentrality()->GetCentralityPercentile("V0M");
	  //09062016 Bernhard Hohlweger
	  Double_t CorrectionFactor = fVZEROgainEqualizationCen->GetBinContent(fVZEROgainEqualizationCen->FindBin(id,EventCentrality));
	  // the fVZEROxpol[] weights are used to enable or disable vzero rings
	  if(id<32) {   // v0c side
	    fTrackEta = -3.45+0.5*(id/8);
		  if(id < 8) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[0];
		  else if (id < 16 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[1];
		  else if (id < 24 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[2];
		  else if (id < 32 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROCpol[3];
		} else {       // v0a side
		  fTrackEta = +4.8-0.6*((id/8)-4);
		  if( id < 40) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[0];
		  else if ( id < 48 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[1];
		  else if ( id < 56 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[2];
		  else if ( id < 64 ) fTrackWeight = fEvent->GetVZEROEqMultiplicity(id)*fVZEROApol[3];
		}
		if(CorrectionFactor) {
		  fTrackWeight /= CorrectionFactor;
	    }
	}else{//for Run2 Data
	  Double_t CorrectionFactor = -1.;
	  AliMultSelection *MultSelection = 0x0;
	  MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
	  if( !MultSelection) {
	    //If you get this warning (and EventCentrality -999) please check that the AliMultSelectionTask actually ran (before your task)
		AliFatal("AliMultSelection not found, did you Run AliMultSelectionTask? \n");
	  }else{
		EventCentrality = MultSelection->GetMultiplicityPercentile("V0M");
		if(EventCentrality < 0 && EventCentrality > 100){
		  AliWarning("No Correction Available for this Centrality \n");
		}else{
		  CorrectionFactor = fVZEROgainEqualizationCen->GetBinContent(fVZEROgainEqualizationCen->FindBin(id,EventCentrality));
		}
	  }
	  if(id<32) {    // v0c side
	    fTrackEta = -3.45+0.5*(id/8);
		if(id < 8) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[0];
		else if (id < 16 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[1];
		else if (id < 24 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[2];
		else if (id < 32 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROCpol[3];
	  } else {       // v0a side
		fTrackEta = +4.8-0.6*((id/8)-4);
		if( id < 40) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[0];
		else if ( id < 48 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[1];
		else if ( id < 56 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[2];
		else if ( id < 64 ) fTrackWeight = fEvent->GetVZEROData()->GetMultiplicity(id)*fVZEROApol[3];
	  }
	  if(CorrectionFactor>0) {
		fTrackWeight /= CorrectionFactor;
	  }
	}
  }//end of if(fVZEROgainEqualizationCen)

  if (fLinearizeVZEROresponse && id < 64)
  {
    //this is only needed in pass1 of LHC10h
    Float_t multVZERO[fgkNumberOfVZEROtracks];
    Float_t dummy=0.0;
    AliESDUtils::GetCorrV0((AliESDEvent*)fEvent,dummy,multVZERO);
    fTrackWeight = multVZERO[id];
  }

  Bool_t pass=kTRUE;
  if (fCutEta) {if (  fTrackEta < fEtaMin || fTrackEta >= fEtaMax ) pass = kFALSE;}
  if (fCutPhi) {if ( fTrackPhi < fPhiMin || fTrackPhi >= fPhiMax ) pass = kFALSE;}

  return pass;
}

//-----------------------------------------------------------------------------
Bool_t AliFlowTrackCuts::PassesMuonCuts(AliVParticle* vparticle)
{
// XZhang 20120604
  ClearTrack();
  fTrackLabel = (fFakesAreOK)?TMath::Abs(vparticle->GetLabel()):vparticle->GetLabel();
  if (fMCevent) fMCparticle = static_cast<AliMCParticle*>(fMCevent->GetTrack(fTrackLabel));
  else fMCparticle=NULL;

  AliESDMuonTrack *esdTrack = dynamic_cast<AliESDMuonTrack*>(vparticle);
  AliAODTrack     *aodTrack = dynamic_cast<AliAODTrack*>(vparticle);
  if ((!esdTrack) && (!aodTrack)) return kFALSE;
  if (!fMuonTrackCuts->IsSelected(vparticle)) return kFALSE;
  HandleVParticle(vparticle);   if (!fTrack)  return kFALSE;
  return kTRUE;
}

//----------------------------------------------------------------------------//
Double_t AliFlowTrackCuts::GetPmdEta(Float_t xPos, Float_t yPos, Float_t zPos)
{
  //get PMD "track" eta coordinate
  Float_t rpxpy, theta, eta;
  rpxpy  = TMath::Sqrt(xPos*xPos + yPos*yPos);
  theta  = TMath::ATan2(rpxpy,zPos);
  eta    = -TMath::Log(TMath::Tan(0.5*theta));
  return eta;
}

//--------------------------------------------------------------------------//
Double_t AliFlowTrackCuts::GetPmdPhi(Float_t xPos, Float_t yPos)
{
  //Get PMD "track" phi coordinate
  Float_t pybypx, phi = 0., phi1;
  if(xPos==0)
    {
      if(yPos>0) phi = 90.;
      if(yPos<0) phi = 270.;
    }
  if(xPos != 0)
    {
      pybypx = yPos/xPos;
      if(pybypx < 0) pybypx = - pybypx;
      phi1 = TMath::ATan(pybypx)*180./3.14159;
      
      if(xPos > 0 && yPos > 0) phi = phi1;        // 1st Quadrant
      if(xPos < 0 && yPos > 0) phi = 180 - phi1;  // 2nd Quadrant
      if(xPos < 0 && yPos < 0) phi = 180 + phi1;  // 3rd Quadrant
      if(xPos > 0 && yPos < 0) phi = 360 - phi1;  // 4th Quadrant
      
    }
  phi = phi*3.14159/180.;
  return   phi;
}

//---------------------------------------------------------------//
void AliFlowTrackCuts::Browse(TBrowser* b)
{
  //some browsing capabilities
  if (fQA) b->Add(fQA);
}

//---------------------------------------------------------------//
Long64_t AliFlowTrackCuts::Merge(TCollection* list)
{
  //merge
  if (!fQA || !list) return 0;
  if (list->IsEmpty()) return 0;
  AliFlowTrackCuts* obj=NULL;
  TList tmplist;
  TIter next(list);
  while ( (obj = dynamic_cast<AliFlowTrackCuts*>(next())) )
  {
    if (obj==this) continue;
    tmplist.Add(obj->GetQA());
  }
  return fQA->Merge(&tmplist);
}
//________________________________________________________________//
void AliFlowTrackCuts::SetTPCTOFNsigmaPIDPurityFunctions(Float_t PurityLevel)
{
    fPurityLevel = PurityLevel;
    cout<<"fPurityLevel = "<<fPurityLevel<<endl;
    
    fPurityFunctionsFile = TFile::Open(Form("alien:///alice/cern.ch/user/n/nmohamma/PurityFunctions_%i-%icent.root",fCentralityPercentileMin,fCentralityPercentileMax));
    if((!fPurityFunctionsFile) || (!fPurityFunctionsFile->IsOpen())) {
        printf("The purity functions file does not exist");
        return;
    }
    fPurityFunctionsList=(TDirectory*)fPurityFunctionsFile->Get("Filterbit1");
    if(!fPurityFunctionsList){printf("The purity functions list is empty"); return;}
    
    TString species[3] = {"pion","kaon","proton"};
    TList *Species_functions[3];
    Int_t ispecie = 0;
    for(ispecie = 0; ispecie < 3; ispecie++) {
        Species_functions[ispecie] = (TList*)fPurityFunctionsList->Get(species[ispecie]);
        if(!Species_functions[ispecie]) {
            cout<<"Purity functions for species: "<<species[ispecie]<<" not found!!!"<<endl;
            return;
        }
    }
    
    for(int i=0;i<180;i++){
        int ispecie = i/60;
        int iPbin = i%60;
        fPurityFunction[i] = (TF2*)Species_functions[ispecie]->FindObject(Form("PurityFunction_%d%d",iPbin,iPbin+1));
        if(!fPurityFunction[i]){printf("Purity function does not exist"); return;}
    }
}
//__________________
Int_t AliFlowTrackCuts::MaxSharedITSClusterCuts(AliESDtrack* track){
    
    Int_t counterForSharedCluster = 0;
    for(int i =0;i<6;i++){
        Bool_t sharedITSCluster = track->HasSharedPointOnITSLayer(i);
//        Bool_t PointOnITSLayer = track->HasPointOnITSLayer(i);
        //   cout << "has a point???    " << PointOnITSLayer << "     is it shared or not??? "  << sharedITSCluster << endl;
        if(sharedITSCluster == 1) counterForSharedCluster++;
    }
    //if(counterForSharedCluster >= maxITSclusterShared) pass=kFALSE;
    
    return counterForSharedCluster;

}
//___________________
Double_t AliFlowTrackCuts::MaxChi2perITSClusterCuts(AliESDtrack* track){
    
    // cout << "  Total Shared Cluster ==   "  <<counterForSharedCluster<< endl;
    if(track->GetNcls(0) == 0) return 999.;
    Double_t chi2perClusterITS = track->GetITSchi2()/track->GetNcls(0);
    //  cout << "  chi2perClusterITS ==   "  <<chi2perClusterITS <<endl;
    
    
    //if(chi2perClusterITS >= maxITSChi2) pass=kFALSE;

    return chi2perClusterITS;
}




