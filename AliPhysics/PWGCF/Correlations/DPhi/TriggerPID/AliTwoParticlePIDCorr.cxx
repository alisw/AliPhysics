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

#include "AliTwoParticlePIDCorr.h"
#include "AliVParticle.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TList.h"
#include "TFile.h"
#include "TGrid.h"
#include "TExMap.h"
#include "AliCentrality.h"
#include "Riostream.h"

#include "AliAnalysisDataSlot.h"
 #include "AliAnalysisDataContainer.h"

#include "AliTHn.h"    
#include "AliCFContainer.h"
#include "THn.h"
#include "THnSparse.h"
#include "TBits.h"
#include <TSpline.h>
#include <AliPID.h>
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include <AliPIDResponse.h>
#include "AliPIDCombined.h"   

#include <AliInputEventHandler.h>
#include "AliAODInputHandler.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "THnSparse.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "TParticle.h"
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliOADBContainer.h"

#include "AliEventPoolManager.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"

#include <random>
using namespace AliPIDNameSpace;
using namespace std;

ClassImp(AliTwoParticlePIDCorr)
ClassImp(LRCParticlePID)


const char * kPIDTypeName[]={"TPC","TOF","TPC-TOF"} ;
const char * kDetectorName[]={"ITS","TPC","TOF"} ;
const char * kParticleSpeciesName[]={"Pions","Kaons","Protons","Undefined"} ;
//Source code::dphicorrelations,VnV0, TaskBFpsi, AliHelperPID, 

//________________________________________________________________________
AliTwoParticlePIDCorr::AliTwoParticlePIDCorr() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fOutput(0),
   fOutputList(0),
  fList(0),
  fCentralityMethod("V0A"),
  fPPVsMult(kFALSE),
 fPileUp_zvtx_INEL_evsel(kTRUE),//as this is by default true in AlippVsMultUtils class for proper event selection
  fSampleType("pPb"),
 fRequestEventPlane(kFALSE),
 fRequestEventPlanemixing(kFALSE),
  fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default)
  trkVtx(0),
  zvtx(0),
  fFilterBit(768),
  fTrackStatus(0),
  fSharedClusterCut(-1),
 fSharedTPCmapCut(-1),
 fSharedfraction_Pair_cut(-1),
  fVertextype(1),
 skipParticlesAbove(0),
  fzvtxcut(10.0),
  fVxMax_MC(0.3),
  fVyMax_MC(0.3),
  fVzMax_MC(10.),
  ffilltrigassoUNID(kFALSE),
  ffilltrigUNIDassoID(kFALSE),
  ffilltrigIDassoUNID(kTRUE),
  ffilltrigIDassoID(kFALSE),
  ffilltrigIDassoIDMCTRUTH(kFALSE),
  fMaxNofMixingTracks(50000),
  fPtOrderMCTruth(kTRUE),
 fPtOrderDataReco(kTRUE),
  fWeightPerEvent(kFALSE),
  fTriggerSpeciesSelection(kFALSE),
  fAssociatedSpeciesSelection(kFALSE),
 fRandomizeReactionPlane(kFALSE),
  fTriggerSpecies(SpPion),
  fAssociatedSpecies(SpPion),
  fCustomBinning(""),
  fBinningString(""),
  fSelectHighestPtTrig(kFALSE),
  fcontainPIDtrig(kTRUE),
  fcontainPIDasso(kFALSE),
  SetChargeAxis(0),
  frejectPileUp(kFALSE),
  fCheckFirstEventInChunk(kFALSE),
  fminPt(0.2),
  fmaxPt(20.0),
  fmineta(-0.8),
  fmaxeta(0.8),
  fselectprimaryTruth(kTRUE),
  fonlyprimarydatareco(kFALSE),
  fdcacutvalue(3.0),
  ffillhistQAReco(kFALSE),
  ffillhistQATruth(kFALSE),
  kTrackVariablesPair(0),
  fminPtTrig(0),
  fmaxPtTrig(0),
  fminPtComboeff(2.0),
  fmaxPtComboeff(4.0), 
  fminPtAsso(0),
  fmaxPtAsso(0),
 fmincentmult(0),
 fmaxcentmult(0), 
 fPriHistShare(0),
  fhistcentrality(0),
 fhistImpactParm(0),
 fhistImpactParmvsMult(0x0),
 fNchNpartCorr(0x0),
  fEventCounter(0),
  fEtaSpectrasso(0),
  fphiSpectraasso(0),
  MCtruthpt(0),
  MCtrutheta(0),
  MCtruthphi(0),
  MCtruthpionpt(0),
  MCtruthpioneta(0),
  MCtruthpionphi(0),
  MCtruthkaonpt(0),
  MCtruthkaoneta(0),
  MCtruthkaonphi(0),
  MCtruthprotonpt(0),
  MCtruthprotoneta(0),
  MCtruthprotonphi(0),
  fPioncont(0),
  fKaoncont(0),
  fProtoncont(0),
  fUNIDcont(0),
  fEventno(0),
  fEventnobaryon(0),
  fEventnomeson(0),
 fhistJetTrigestimate(0),
fTwoTrackDistancePtdip(0x0),
fTwoTrackDistancePtdipmix(0x0),
  fCentralityCorrelation(0x0),
   fCentralityCorrelationMC(0x0),
 fHistVZEROAGainEqualizationMap(0),
  fHistVZEROCGainEqualizationMap(0),
 fHistVZEROChannelGainEqualizationMap(0),
fCentralityWeights(0),
 fHistCentStats(0x0),
 fHistRefmult(0x0),
 fHistEQVZEROvsTPCmultiplicity(0x0),
    fHistEQVZEROAvsTPCmultiplicity(0x0),
    fHistEQVZEROCvsTPCmultiplicity(0x0),
    fHistVZEROCvsEQVZEROCmultiplicity(0x0),
    fHistVZEROAvsEQVZEROAmultiplicity(0x0),
    fHistVZEROCvsVZEROAmultiplicity(0x0),
    fHistEQVZEROCvsEQVZEROAmultiplicity(0x0),
    fHistVZEROSignal(0x0),
fHistEventPlaneTruth(0x0),
fHistPsiMinusPhi(0x0),
fEventPlanePID(0x0),
evplaneMC(999.),
 fgPsi2v0a(999.),
    fgPsi2v0c(999.),
    fgPsi2tpc(999.),
    fgPsi3v0a(999.),
    fgPsi3v0c(999.),
    fgPsi3tpc(999.),
    fgPsi2v0aMC(999.),
    fgPsi2v0cMC(999.),
    fgPsi2tpcMC(999.),
    fgPsi3v0aMC(999.),
    fgPsi3v0cMC(999.),
    fgPsi3tpcMC(999.),
 gReactionPlane(999.),
  fV2(kTRUE),
 fV3(kFALSE),
 fIsAfter2011(kTRUE),
 kNCent_ds(-1),
  fRun(-1),
  fNcluster(70),
 fEPdet("V0A"),  
 fMultV0(NULL),
  fV0Cpol(100),
  fV0Apol(100),
 fHResTPCv0A2(NULL),
fHResTPCv0C2(NULL),
fHResv0Cv0A2(NULL),
fHResTPCv0A3(NULL),
fHResTPCv0C3(NULL),
fHResv0Cv0A3(NULL),
 fHResMA2(NULL),
fHResMC2(NULL),
fHResAC2(NULL),
fHResMA3(NULL),
fHResMC3(NULL),
fHResAC3(NULL),
fPhiRPTPC(NULL),
fPhiRPTPCv3(NULL),
fPhiRPv0A(NULL),
fPhiRPv0C(NULL),
fPhiRPv0Av3(NULL),
fPhiRPv0Cv3(NULL),
 fControlConvResoncances(0),
  fHistoTPCdEdx(0x0),
  fHistoTOFbeta(0x0),
  fTPCTOFPion3d(0),
  fTPCTOFKaon3d(0),
  fTPCTOFProton3d(0),
  fPionPt(0),
  fPionEta(0),
  fPionPhi(0),
  fKaonPt(0),
  fKaonEta(0),
  fKaonPhi(0),
  fProtonPt(0),
  fProtonEta(0),
  fProtonPhi(0),
 kShortSigmahisto(0),
LambdaSigmahisto(0),
fHistdEdxVsPTPCbeforePIDelectron(NULL),
  fHistNSigmaTPCvsPtbeforePIDelectron(NULL),
  fHistdEdxVsPTPCafterPIDelectron(NULL),
  fHistNSigmaTPCvsPtafterPIDelectron(NULL),
  fCorrelatonTruthPrimary(0),
  fCorrelatonTruthPrimarymix(0),
  fTHnCorrUNID(0),
  fTHnCorrUNIDmix(0),
  fTHnCorrID(0),
  fTHnCorrIDmix(0),
  fTHnCorrIDUNID(0),
  fTHnCorrIDUNIDmix(0),
  fTHnTrigcount(0),
  fTHnTrigcountMCTruthPrim(0),
  fPoolMgr(0x0),
  fArrayMC(0),
  fAnalysisType("AOD"), 
  fefffilename(""),
 ffilenamesigmaV0(""),
 ftwoTrackEfficiencyCutDataReco(kTRUE),
fTwoTrackCutMinRadius(0.8),
fTwoTrackCutMaxRadius(2.5),
  twoTrackEfficiencyCutValue(0.02),
  fPID(NULL),
 fPIDCombined(NULL),
 eventno(0),
  fPtTOFPIDmin(0.5),
  fPtTOFPIDmax(6.0),
  fRequestTOFPID(kTRUE),
  fPIDType(NSigmaTPCTOF),
 fFIllPIDQAHistos(kTRUE),
  fNSigmaPID(3),
  fBayesCut(0.8),
 fdiffPIDcutvalues(kFALSE),
 fPt1(2.0),//set Pt ranges for diff. pid cut values < fPtTOFPIDMax
 fPt2(3.5),
 fPt3(4.0),
 fPIDCutval1(0.0),
 fPIDCutval2(0.0),
 fPIDCutval3(0.0),
 fPIDCutval4(0.0),
 fHighPtKaonNSigmaPID(-1),
 fHighPtKaonSigma(3.5),
  fUseExclusiveNSigma(kFALSE), 
  fRemoveTracksT0Fill(kFALSE),
fSelectCharge(0),
fTriggerSelectCharge(0),
fAssociatedSelectCharge(0),
fTriggerRestrictEta(-1),
fEtaOrdering(kFALSE),
fCutConversions(kFALSE),
fCutResonances(kFALSE),
fRejectResonanceDaughters(-1),
  fOnlyOneEtaSide(0),
fInjectedSignals(kFALSE),
  fRemoveWeakDecays(kFALSE),
fRemoveDuplicates(kFALSE),
  fapplyTrigefficiency(kFALSE),
  fapplyAssoefficiency(kFALSE),
  ffillefficiency(kFALSE),
  fmesoneffrequired(kFALSE),
  fkaonprotoneffrequired(kFALSE),
 fRunShufflingbalance(kFALSE),
 //for electron rejection at single rack level
   fElectronRejectionTRUTH(kFALSE),
  fElectronRejection(kFALSE),
  fElectronOnlyRejection(kFALSE),
  fElectronRejectionNSigma(-1.),
  fElectronRejectionMinPt(0.),
  fElectronRejectionMaxPt(1000.), 
fAnalysisUtils(0x0),
 fPPVsMultUtils(0x0),
 fDCAXYCut(0),
 fV0TrigCorr(kFALSE),
  ffillofflineV0(kFALSE),
 fUsev0DaughterPID(kFALSE),
  fMinPtDaughter(1.0),// v0 related cut starts here
  fMaxPtDaughter(4.0),
  fDCAToPrimVtx(0.1),
  fMaxDCADaughter(1.0),
  fMinCPA(0.998),
  lMax(100),
  fHistRawPtCentInvK0s(0x0),
  fHistRawPtCentInvLambda(0x0),
  fHistRawPtCentInvAntiLambda(0x0),
  fHistFinalPtCentInvK0s(0x0),
  fHistFinalPtCentInvLambda(0x0),
  fHistFinalPtCentInvAntiLambda(0x0),
 fCtauCut3D(kTRUE),
  NCtau(3.0),
fCutctauK0s(2.68),
  fCutctauLambda(7.89),
  fCutctauAntiLambda(7.89),
  fRapCutK0s(0.7),
  fRapCutLambda(0.7),
fDaugNClsTPC(70),
 fFracTPCcls(0),
 fCutDaughterPtV0(kFALSE),
 TPCSectoredgecut(kTRUE),//high Pt rel. rise TPCPid related variables
 fNclsusedfordEdXdtr(60),
 fPiondeltacutmin(1.0),
 fPiondeltacutmax(7.0),  
 fProtondeltacutmin(-21),//-25 for systematics
 fProtondeltacutmax(-14),
 deltapion_val(-999999)

{
 
  for ( Int_t i = 0; i < 16; i++) { 
    fHistQA[i] = NULL;
  }

 for ( Int_t i = 0; i < 6; i++ ){
    fTrackHistEfficiency[i] = NULL;
    effcorection[i]=NULL;
    //effmap[i]=NULL;
  }
 for ( Int_t i = 0; i < 2; i++ ){
   fTwoTrackDistancePt[i]=NULL;
   fTwoTrackDistancePtmix[i]=NULL;
}

 for(Int_t ipart=0;ipart<NSpecies;ipart++)
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++)
      fnsigmas[ipart][ipid]=999.;

 for(Int_t ipart=0;ipart<NSpecies;ipart++) {fHasDoubleCounting[ipart]=kFALSE;}

  for(Int_t i = 0; i != 2; ++i)
    for(Int_t j = 0; j != 2; ++j)
      for(Int_t iC = 0; iC < 9; iC++){
	fMeanQ[iC][i][j] = 0.;
	fWidthQ[iC][i][j] = 1.;
        fMeanQv3[iC][i][j] = 0.;
	fWidthQv3[iC][i][j] = 1.;
    }

  }
//________________________________________________________________________
AliTwoParticlePIDCorr::AliTwoParticlePIDCorr(const char *name) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
 fOutput(0),
   fOutputList(0),
   fList(0),
 fCentralityMethod("V0A"),
  fPPVsMult(kFALSE),
 fPileUp_zvtx_INEL_evsel(kTRUE),//as this is by default true in AlippVsMultUtils class for proper event selection
  fSampleType("pPb"),
 fRequestEventPlane(kFALSE),
 fRequestEventPlanemixing(kFALSE),
  fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default)
  trkVtx(0),
  zvtx(0),
  fFilterBit(768),
  fTrackStatus(0),
  fSharedClusterCut(-1),
  fSharedTPCmapCut(-1),
 fSharedfraction_Pair_cut(-1),
  fVertextype(1),
   skipParticlesAbove(0),
  fzvtxcut(10.0),
  fVxMax_MC(0.3),
  fVyMax_MC(0.3),
  fVzMax_MC(10.),
  ffilltrigassoUNID(kFALSE),
  ffilltrigUNIDassoID(kFALSE),
  ffilltrigIDassoUNID(kTRUE),
  ffilltrigIDassoID(kFALSE),
  ffilltrigIDassoIDMCTRUTH(kFALSE),
  fMaxNofMixingTracks(50000),
  fPtOrderMCTruth(kTRUE),
  fPtOrderDataReco(kTRUE),
  fWeightPerEvent(kFALSE),
  fTriggerSpeciesSelection(kFALSE),
  fAssociatedSpeciesSelection(kFALSE),
   fRandomizeReactionPlane(kFALSE),
  fTriggerSpecies(SpPion),
  fAssociatedSpecies(SpPion),
  fCustomBinning(""),
  fBinningString(""),
  fSelectHighestPtTrig(kFALSE),
  fcontainPIDtrig(kTRUE),
  fcontainPIDasso(kFALSE),
  SetChargeAxis(0),
  frejectPileUp(kFALSE),
  fCheckFirstEventInChunk(kFALSE),
  fminPt(0.2),
  fmaxPt(20.0),
  fmineta(-0.8),
  fmaxeta(0.8),
  fselectprimaryTruth(kTRUE),
  fonlyprimarydatareco(kFALSE),
  fdcacutvalue(3.0),
  ffillhistQAReco(kFALSE),
  ffillhistQATruth(kFALSE),
 kTrackVariablesPair(0),
  fminPtTrig(0),
  fmaxPtTrig(0),
  fminPtComboeff(2.0),
  fmaxPtComboeff(4.0), 
  fminPtAsso(0),
  fmaxPtAsso(0),
   fmincentmult(0),
   fmaxcentmult(0),
   fPriHistShare(0),
  fhistcentrality(0),
   fhistImpactParm(0),
   fhistImpactParmvsMult(0x0),
   fNchNpartCorr(0x0),
  fEventCounter(0),
  fEtaSpectrasso(0),
  fphiSpectraasso(0),
  MCtruthpt(0),
  MCtrutheta(0),
  MCtruthphi(0),
  MCtruthpionpt(0),
  MCtruthpioneta(0),
  MCtruthpionphi(0),
  MCtruthkaonpt(0),
  MCtruthkaoneta(0),
  MCtruthkaonphi(0),
  MCtruthprotonpt(0),
  MCtruthprotoneta(0),
  MCtruthprotonphi(0),
  fPioncont(0),
  fKaoncont(0),
  fProtoncont(0),
   fUNIDcont(0),
  fEventno(0),
  fEventnobaryon(0),
  fEventnomeson(0),
  fhistJetTrigestimate(0),
fTwoTrackDistancePtdip(0x0),
fTwoTrackDistancePtdipmix(0x0),
  fCentralityCorrelation(0x0),
     fCentralityCorrelationMC(0x0),
 fHistVZEROAGainEqualizationMap(0),
  fHistVZEROCGainEqualizationMap(0),
   fHistVZEROChannelGainEqualizationMap(0),
fCentralityWeights(0),
  fHistCentStats(0x0),
  fHistRefmult(0x0),
    fHistEQVZEROvsTPCmultiplicity(0x0),
    fHistEQVZEROAvsTPCmultiplicity(0x0),
    fHistEQVZEROCvsTPCmultiplicity(0x0),
    fHistVZEROCvsEQVZEROCmultiplicity(0x0),
    fHistVZEROAvsEQVZEROAmultiplicity(0x0),
    fHistVZEROCvsVZEROAmultiplicity(0x0),
    fHistEQVZEROCvsEQVZEROAmultiplicity(0x0),
    fHistVZEROSignal(0x0),
fHistEventPlaneTruth(0x0),
   fHistPsiMinusPhi(0x0),
fEventPlanePID(0x0),
evplaneMC(999.),
 fgPsi2v0a(999.),
    fgPsi2v0c(999.),
    fgPsi2tpc(999.),
    fgPsi3v0a(999.),
    fgPsi3v0c(999.),
    fgPsi3tpc(999.),
    fgPsi2v0aMC(999.),
    fgPsi2v0cMC(999.),
    fgPsi2tpcMC(999.),
    fgPsi3v0aMC(999.),
    fgPsi3v0cMC(999.),
    fgPsi3tpcMC(999.),
   gReactionPlane(999.),
 fV2(kTRUE),
 fV3(kFALSE),
 fIsAfter2011(kTRUE),
   kNCent_ds(-1),
  fRun(-1),
  fNcluster(70),
   fEPdet("V0A"),  
 fMultV0(NULL),
  fV0Cpol(100),
  fV0Apol(100),
 fHResTPCv0A2(NULL),
fHResTPCv0C2(NULL),
fHResv0Cv0A2(NULL),
fHResTPCv0A3(NULL),
fHResTPCv0C3(NULL),
fHResv0Cv0A3(NULL),
 fHResMA2(NULL),
fHResMC2(NULL),
fHResAC2(NULL),
fHResMA3(NULL),
fHResMC3(NULL),
fHResAC3(NULL),
fPhiRPTPC(NULL),
fPhiRPTPCv3(NULL),
fPhiRPv0A(NULL),
fPhiRPv0C(NULL),
fPhiRPv0Av3(NULL),
fPhiRPv0Cv3(NULL),
  fControlConvResoncances(0), 
  fHistoTPCdEdx(0x0),
  fHistoTOFbeta(0x0),
  fTPCTOFPion3d(0),
  fTPCTOFKaon3d(0),
  fTPCTOFProton3d(0),
  fPionPt(0),
  fPionEta(0),
  fPionPhi(0),
  fKaonPt(0),
  fKaonEta(0),
  fKaonPhi(0),
  fProtonPt(0),
  fProtonEta(0),
  fProtonPhi(0),
  kShortSigmahisto(0),
LambdaSigmahisto(0),
     fHistdEdxVsPTPCbeforePIDelectron(NULL),
  fHistNSigmaTPCvsPtbeforePIDelectron(NULL),
  fHistdEdxVsPTPCafterPIDelectron(NULL),
  fHistNSigmaTPCvsPtafterPIDelectron(NULL),
  fCorrelatonTruthPrimary(0),
 fCorrelatonTruthPrimarymix(0),
  fTHnCorrUNID(0),
  fTHnCorrUNIDmix(0),
  fTHnCorrID(0),
  fTHnCorrIDmix(0),
  fTHnCorrIDUNID(0),
  fTHnCorrIDUNIDmix(0),
  fTHnTrigcount(0),
  fTHnTrigcountMCTruthPrim(0),
  fPoolMgr(0x0),
  fArrayMC(0),
  fAnalysisType("AOD"),
  fefffilename(""),
  ffilenamesigmaV0(""),
  ftwoTrackEfficiencyCutDataReco(kTRUE),
fTwoTrackCutMinRadius(0.8),
fTwoTrackCutMaxRadius(2.5),
  twoTrackEfficiencyCutValue(0.02),
  fPID(NULL),
  fPIDCombined(NULL),
  eventno(0),
 fPtTOFPIDmin(0.5),
  fPtTOFPIDmax(6.0),
  fRequestTOFPID(kTRUE),
  fPIDType(NSigmaTPCTOF),
  fFIllPIDQAHistos(kTRUE),
  fNSigmaPID(3),
  fBayesCut(0.8),
 fdiffPIDcutvalues(kFALSE),
 fPt1(2.0),
 fPt2(3.5),
 fPt3(4.0),
 fPIDCutval1(0.0),
 fPIDCutval2(0.0),
 fPIDCutval3(0.0),
 fPIDCutval4(0.0),
fHighPtKaonNSigmaPID(-1),
 fHighPtKaonSigma(3.5),
  fUseExclusiveNSigma(kFALSE),   
  fRemoveTracksT0Fill(kFALSE),
fSelectCharge(0),
fTriggerSelectCharge(0),
fAssociatedSelectCharge(0),
fTriggerRestrictEta(-1),
fEtaOrdering(kFALSE),
fCutConversions(kFALSE),
fCutResonances(kFALSE),
fRejectResonanceDaughters(-1),
  fOnlyOneEtaSide(0),
fInjectedSignals(kFALSE),
  fRemoveWeakDecays(kFALSE),
fRemoveDuplicates(kFALSE),
  fapplyTrigefficiency(kFALSE),
  fapplyAssoefficiency(kFALSE),
  ffillefficiency(kFALSE),
 fmesoneffrequired(kFALSE),
 fkaonprotoneffrequired(kFALSE),
    fRunShufflingbalance(kFALSE),
//for electron rejection at single rack level
      fElectronRejectionTRUTH(kFALSE),
  fElectronRejection(kFALSE),
  fElectronOnlyRejection(kFALSE),
  fElectronRejectionNSigma(-1.),
  fElectronRejectionMinPt(0.),
  fElectronRejectionMaxPt(1000.),
   fAnalysisUtils(0x0),
   fPPVsMultUtils(0x0),
   fDCAXYCut(0),
   fV0TrigCorr(kFALSE),
   ffillofflineV0(kFALSE),
 fUsev0DaughterPID(kFALSE),
  fMinPtDaughter(1.0),// v0 related cut starts here
  fMaxPtDaughter(4.0),
  fDCAToPrimVtx(0.1),
  fMaxDCADaughter(1.0),
  fMinCPA(0.998),
  lMax(100),
  fHistRawPtCentInvK0s(0x0),
  fHistRawPtCentInvLambda(0x0),
  fHistRawPtCentInvAntiLambda(0x0),
  fHistFinalPtCentInvK0s(0x0),
  fHistFinalPtCentInvLambda(0x0),
  fHistFinalPtCentInvAntiLambda(0x0),
  fCtauCut3D(kTRUE),
  NCtau(3.0),
fCutctauK0s(2.68),
  fCutctauLambda(7.89),
  fCutctauAntiLambda(7.89),
  fRapCutK0s(0.7),
  fRapCutLambda(0.7),
fDaugNClsTPC(70),
   fFracTPCcls(0),
   fCutDaughterPtV0(kFALSE),
 TPCSectoredgecut(kTRUE),//high Pt rel. rise TPCPid related variables
   fNclsusedfordEdXdtr(60),
   fPiondeltacutmin(1.0),
 fPiondeltacutmax(7.0),  
 fProtondeltacutmin(-21),//-25 for systematics
   fProtondeltacutmax(-14),
    deltapion_val(-999999)

    
  // The last in the above list should not have a comma after it
     
{
  
   for ( Int_t i = 0; i < 16; i++) { 
    fHistQA[i] = NULL;
  }
 
for ( Int_t i = 0; i < 6; i++ ){
    fTrackHistEfficiency[i] = NULL;
    effcorection[i]=NULL;
    //effmap[i]=NULL;
  }

for ( Int_t i = 0; i < 2; i++ ){
   fTwoTrackDistancePt[i]=NULL;
   fTwoTrackDistancePtmix[i]=NULL;
}

 for(Int_t ipart=0;ipart<NSpecies;ipart++)
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++)
      fnsigmas[ipart][ipid]=999.;

   for(Int_t ipart=0;ipart<NSpecies;ipart++) {fHasDoubleCounting[ipart]=kFALSE;}

  for(Int_t i = 0; i != 2; ++i)
    for(Int_t j = 0; j != 2; ++j)
      for(Int_t iC = 0; iC < 9; iC++){
	fMeanQ[iC][i][j] = 0.;
	fWidthQ[iC][i][j] = 1.;
        fMeanQv3[iC][i][j] = 0.;
	fWidthQv3[iC][i][j] = 1.;
    }
  
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
     DefineInput(0, TChain::Class());

  DefineOutput(1, TList::Class());                                        // for output list
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());

}

//________________________________________________________________________
AliTwoParticlePIDCorr::~AliTwoParticlePIDCorr()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;

  }

if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;

  }

if(fRequestEventPlane){
if (fList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fList;
 }
 }

  if (fPID) delete fPID;
  if (fPIDCombined) delete fPIDCombined;

  }
//________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////////////////////////

TH2F* AliTwoParticlePIDCorr::GetHistogram2D(const char * name){
  // returns histo named name
  return (TH2F*) fOutputList->FindObject(name);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Float_t AliTwoParticlePIDCorr::PhiRange(Float_t DPhi)

{
	//
	// Puts the argument in the range [-pi/2,3 pi/2].
	//
	
	if (DPhi < -TMath::Pi()/2) DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	

	return DPhi;
	
}
//________________________________________________________________________
void AliTwoParticlePIDCorr::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

// global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nPsiTOF = 10;  
  const Int_t nCentrBin = 9;  


  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!  

  fOutputList = new TList;
  fOutputList->SetOwner();
  fOutputList->SetName("PIDQAList");

  if(fRequestEventPlane){
  fList = new TList;
  fList->SetOwner();
  fList->SetName("EPQAList");
  }
  fEventCounter = new TH1F("fEventCounter","EventCounter", 19, 0.5,19.5);
  fEventCounter->GetXaxis()->SetBinLabel(1,"Event Accesed");
  fEventCounter->GetXaxis()->SetBinLabel(3,"After PileUP Cut");//only for Data
  fEventCounter->GetXaxis()->SetBinLabel(5,"Have A Vertex");
  fEventCounter->GetXaxis()->SetBinLabel(7,"After vertex Cut");
  fEventCounter->GetXaxis()->SetBinLabel(9,"Getting centrality");
  fEventCounter->GetXaxis()->SetBinLabel(11,"After centrality flattening");
  fEventCounter->GetXaxis()->SetBinLabel(13,"Within 0-100% centrality");
  fEventCounter->GetXaxis()->SetBinLabel(15,"Event Analyzed");
  //fEventCounter->GetXaxis()->SetBinLabel(8,"Event Analysis finished");
  fOutput->Add(fEventCounter);
  
fEtaSpectrasso=new TH2F("fEtaSpectraasso","fEtaSpectraasso",180,-0.9,0.9,200,0.0,fmaxPt);
fOutput->Add(fEtaSpectrasso);

fphiSpectraasso=new TH2F("fphiSpectraasso","fphiSpectraasso",72,0,2*TMath::Pi(),200,0.0,fmaxPt);
fOutput->Add(fphiSpectraasso);

 if(fSampleType=="pPb" || fSampleType=="PbPb" || fPPVsMult==kTRUE || fCentralityMethod == "MC_b"){
   if  (fAnalysisType =="MCAOD" ||  fAnalysisType =="MC"){
     fCentralityCorrelationMC = new TH2D("fCentralityCorrelationMC", ";centrality_ImpactParam;multiplicity", 101, 0, 101, 20000, 0,40000);
     fOutput->Add(fCentralityCorrelationMC);
   }

   if  (fAnalysisType =="MCAOD" ||  fAnalysisType =="AOD"){
     fCentralityCorrelation = new TH2D("fCentralityCorrelation", ";centrality_ImpactParam;multiplicity", 101, 0, 101, 20000, 0,40000);
     fOutput->Add(fCentralityCorrelation);
   }
      
 }

if(fCentralityMethod=="V0M" || fCentralityMethod=="V0A" || fCentralityMethod=="V0C" || fCentralityMethod=="CL1" || fCentralityMethod=="ZNA" || fCentralityMethod=="V0AEq" || fCentralityMethod=="V0CEq" || fCentralityMethod=="V0MEq")
  {
 TString gCentName[8] = {"V0A","V0C","V0M","V0AEq","V0CEq","V0MEq","CL1","ZNA"};
  fHistCentStats = new TH2F("fHistCentStats",
                             "Centrality statistics;;Cent percentile",
			    8,-0.5,7.5,220,-5,105);
  for(Int_t i = 1; i <= 8; i++){
    fHistCentStats->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
    //fHistCentStatsUsed->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  }
  fOutput->Add(fHistCentStats);
  }

if(fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMult==kFALSE))
  {
fhistcentrality=new TH1F("fhistcentrality","referencemultiplicity",3001,-0.5,30000.5);
fOutput->Add(fhistcentrality);
  }
 else{
fhistcentrality=new TH1F("fhistcentrality","centrality",220,-5,105);
fOutput->Add(fhistcentrality);
 }
 if(fCentralityMethod=="MC_b"){
fhistImpactParm=new TH1F("fhistImpactParm","Impact_Parameter",300,0,300);
fOutput->Add(fhistImpactParm);
fhistImpactParmvsMult=new TH2F("fhistImpactParmvsMult","fhistImpactParmvsMult",300,0,300,4001,-0.5,40000.5);
fOutput->Add(fhistImpactParmvsMult);
 }

 if(fAnalysisType =="MCAOD"){//for "MC" this creates a lot of problems which were not solved so removed from "MC" case
 if(fSampleType=="pPb" || fSampleType=="PbPb"){
fNchNpartCorr=new TH2F("fNchNpartCorr","fNchNpartCorr",500,0.0,500.0,4001,-0.5,40000.5);
fNchNpartCorr->GetXaxis()->SetTitle("Npart (a.u.)");
fNchNpartCorr->GetYaxis()->SetTitle("Nch(a.u.)");
fOutput->Add(fNchNpartCorr);
   }
 }

 TString gmultName[4] = {"V0A_MANUAL","V0C_MANUAL","V0M_MANUAL","TRACKS_MANUAL"};
  fHistRefmult = new TH2F("fHistRefmult",
                             "Reference multiplicity",
			    4,-0.5,3.5,10000,0,20000);
  for(Int_t i = 1; i <= 4; i++){
    fHistRefmult->GetXaxis()->SetBinLabel(i,gmultName[i-1].Data());
  }
  fOutput->Add(fHistRefmult);

 if(fCentralityMethod == "V0A_MANUAL" || fCentralityMethod == "V0M_MANUAL" || fCentralityMethod == "V0C_MANUAL" ){
 //TPC vs EQVZERO multiplicity
    fHistEQVZEROvsTPCmultiplicity = new TH2F("fHistEQVZEROvsTPCmultiplicity","EqVZERO vs TPC multiplicity",10001,-0.5,10000.5,4001,-0.5,4000.5);
    fHistEQVZEROvsTPCmultiplicity->GetXaxis()->SetTitle("EqVZERO multiplicity (a.u.)");
    fHistEQVZEROvsTPCmultiplicity->GetYaxis()->SetTitle("TPC multiplicity (a.u.)");
    fOutput->Add(fHistEQVZEROvsTPCmultiplicity);


    fHistEQVZEROAvsTPCmultiplicity = new TH2F("fHistEQVZEROAvsTPCmultiplicity","EqVZERO_A vs TPC multiplicity",10001,-0.5,10000.5,4001,-0.5,4000.5);
    fHistEQVZEROAvsTPCmultiplicity->GetXaxis()->SetTitle("EqVZERO_A multiplicity (a.u.)");
    fHistEQVZEROAvsTPCmultiplicity->GetYaxis()->SetTitle("TPC multiplicity (a.u.)");
    fOutput->Add(fHistEQVZEROAvsTPCmultiplicity);


    fHistEQVZEROCvsTPCmultiplicity = new TH2F("fHistEQVZEROCvsTPCmultiplicity","EqVZERO_C vs TPC multiplicity",10001,-0.5,10000.5,4001,-0.5,4000.5);
    fHistEQVZEROCvsTPCmultiplicity->GetXaxis()->SetTitle("EqVZERO_C multiplicity (a.u.)");
    fHistEQVZEROCvsTPCmultiplicity->GetYaxis()->SetTitle("TPC multiplicity (a.u.)");
    fOutput->Add(fHistEQVZEROCvsTPCmultiplicity);

 //EQVZERO vs VZERO multiplicity
  fHistVZEROCvsEQVZEROCmultiplicity = new TH2F("fHistVZEROCvsEQVZEROCmultiplicity","EqVZERO_C vs VZERO_C multiplicity",10001,-0.5,10000.5,10001,-0.5,10000.5);
    fHistVZEROCvsEQVZEROCmultiplicity->GetXaxis()->SetTitle("VZERO_C multiplicity (a.u.)");
    fHistVZEROCvsEQVZEROCmultiplicity->GetYaxis()->SetTitle("EqVZERO_C multiplicity (a.u.)");
    fOutput->Add(fHistVZEROCvsEQVZEROCmultiplicity);


fHistVZEROAvsEQVZEROAmultiplicity = new TH2F("fHistVZEROAvsEQVZEROAmultiplicity","EqVZERO_A vs VZERO_A multiplicity",10001,-0.5,10000.5,10001,-0.5,10000.5);
    fHistVZEROAvsEQVZEROAmultiplicity->GetXaxis()->SetTitle("VZERO_A multiplicity (a.u.)");
    fHistVZEROAvsEQVZEROAmultiplicity->GetYaxis()->SetTitle("EqVZERO_A multiplicity (a.u.)");
    fOutput->Add(fHistVZEROAvsEQVZEROAmultiplicity);


  //VZEROC vs VZEROA multiplicity
fHistVZEROCvsVZEROAmultiplicity = new TH2F("fHistVZEROCvsVZEROAmultiplicity","VZERO_C vs VZERO_A multiplicity",10001,-0.5,10000.5,10001,-0.5,10000.5);
    fHistVZEROCvsVZEROAmultiplicity->GetXaxis()->SetTitle("VZERO_C multiplicity (a.u.)");
    fHistVZEROCvsVZEROAmultiplicity->GetYaxis()->SetTitle("VZERO_A multiplicity (a.u.)");
    fOutput->Add(fHistVZEROCvsVZEROAmultiplicity);



  //EQVZEROC vs EQVZEROA multiplicity
fHistEQVZEROCvsEQVZEROAmultiplicity = new TH2F("fHistEQVZEROCvsEQVZEROAmultiplicity","EqVZERO_C vs EqVZERO_A multiplicity",10001,-0.5,10000.5,10001,-0.5,10000.5);
    fHistEQVZEROCvsEQVZEROAmultiplicity->GetXaxis()->SetTitle("EqVZERO_C multiplicity (a.u.)");
    fHistEQVZEROCvsEQVZEROAmultiplicity->GetYaxis()->SetTitle("EqVZERO_A multiplicity (a.u.)");
    fOutput->Add(fHistEQVZEROCvsEQVZEROAmultiplicity);

 fHistVZEROSignal = new TH2F("fHistVZEROSignal","VZERO signal vs VZERO channel;VZERO channel; Signal (a.u.)",64,0.5,64.5,3001,-0.5,30000.5);
  fOutput->Add(fHistVZEROSignal);
 }


 if(fRequestEventPlane){
//Event plane
 
  fHistPsiMinusPhi = new TH2D("fHistPsiMinusPhi","",4,-0.5,3.5,100,0,2.*TMath::Pi());
  fList->Add(fHistPsiMinusPhi);

  fEventPlanePID = new TH3F("fEventPlanePID",";centrality;eventplane;PID",20,0.0,100.0,4,-0.5,3.5,4,-0.5,3.5);
  fList->Add(fEventPlanePID);

 }
 
if(fCutConversions || fCutResonances)
    {
fControlConvResoncances = new TH2F("fControlConvResoncances", ";id;delta mass", 3, -0.5, 2.5, 100, -0.1, 0.1);
 fOutput->Add(fControlConvResoncances);
    }

fHistoTPCdEdx = new TH2F("fHistoTPCdEdx", ";p_{T} (GeV/c);dE/dx (au.)",200,0.0,fmaxPt,500, 0., 500.);
fOutputList->Add(fHistoTPCdEdx);
fHistoTOFbeta = new TH2F(Form("fHistoTOFbeta"), ";p_{T} (GeV/c);v/c",200, 0.0, fmaxPt, 500, 0.1, 1.1);
  fOutputList->Add(fHistoTOFbeta);
  
   fTPCTOFPion3d=new TH3F ("fTPCTOFpion3d", "fTPCTOFpion3d",200,0.0,fmaxPt, 120,-60.,60.,120,-60.,60);
   fOutputList->Add(fTPCTOFPion3d);
  
   fTPCTOFKaon3d=new TH3F ("fTPCTOFKaon3d", "fTPCTOFKaon3d",200,0.0,fmaxPt, 120,-60.,60.,120,-60.,60);
   fOutputList->Add(fTPCTOFKaon3d);

   fTPCTOFProton3d=new TH3F ("fTPCTOFProton3d", "fTPCTOFProton3d",200,0.0,fmaxPt, 120,-60.,60.,120,-60.,60);
   fOutputList->Add(fTPCTOFProton3d);

if(ffillhistQAReco)
    {
    fPionPt = new TH1F("fPionPt","p_{T} distribution",200,0.0,fmaxPt);
 fOutputList->Add(fPionPt);
    fPionEta= new TH1F("fPionEta","#eta distribution",360,-1.8,1.8);
 fOutputList->Add(fPionEta);
    fPionPhi = new TH1F("fPionPhi","#phi distribution",340,0,6.8);
 fOutputList->Add(fPionPhi);
  
    fKaonPt = new TH1F("fKaonPt","p_{T} distribution",200,0.0,fmaxPt);
 fOutputList->Add(fKaonPt);
    fKaonEta= new TH1F("fKaonEta","#eta distribution",360,-1.8,1.8);
 fOutputList->Add(fKaonEta);
    fKaonPhi = new TH1F("fKaonPhi","#phi distribution",340,0,6.8);
 fOutputList->Add(fKaonPhi);
  
    fProtonPt = new TH1F("fProtonPt","p_{T} distribution",200,0.0,fmaxPt);
 fOutputList->Add(fProtonPt);
    fProtonEta= new TH1F("fProtonEta","#eta distribution",360,-1.8,1.8);
 fOutputList->Add(fProtonEta);
    fProtonPhi= new TH1F("fProtonPhi","#phi distribution",340,0,6.8);
 fOutputList->Add(fProtonPhi);
    }

  fHistQA[0] = new TH1F("fHistQAvx", "Histo Vx All ", 50, -5., 5.);
  fHistQA[1] = new TH1F("fHistQAvy", "Histo Vy All", 50, -5., 5.);
  fHistQA[2] = new TH1F("fHistQAvz", "Histo Vz All", 50, -25., 25.);  
  fHistQA[3] = new TH1F("fHistQAvxA", "Histo Vx  After Cut ", 50, -5., 5.);
  fHistQA[4] = new TH1F("fHistQAvyA", "Histo Vy After Cut", 50, -5., 5.);
  fHistQA[5] = new TH1F("fHistQAvzA", "Histo Vz After Cut", 50, -25., 25.);
  fHistQA[6] = new TH1F("fHistQADcaXyC", "Histo DCAxy after cut", 50, -5., 5.);
  fHistQA[7] = new TH1F("fHistQADcaZC", "Histo DCAz after cut", 50, -5., 5.);   
  fHistQA[8] = new TH1F("fHistQAPt","p_{T} distribution",200,0.0,fmaxPt);
  fHistQA[9] = new TH1F("fHistQAEta","#eta distribution",360,-1.8,1.8);
  fHistQA[10] = new TH1F("fHistQAPhi","#phi distribution",340,0,6.8);
  fHistQA[11] = new TH1F("fHistQANCls","Number of TPC cluster",200,0,200);
  fHistQA[13] = new TH1F("fHistQAChi2","Chi2 per NDF",100,0,10);
 fHistQA[12] = new TH1F("fHistQANCls1","Number of TPC cluster1",200,0,200);
 fHistQA[14] = new TH1F("nCrossedRowsTPC","Number of TPC ccrossed rows",200,0,200);
 fHistQA[15] = new TH1F("ratioCrossedRowsOverFindableClustersTPC","Number of TPC ccrossed rows find clusters",200,0,2);
    
for(Int_t i = 0; i < 16; i++)
    {
      fOutput->Add(fHistQA[i]);
    }

    fPriHistShare = new TH1F ("fPriHistShare","Shared clusters, primaries;#shared clusters;counts",160,0,160);
    fOutput->Add(fPriHistShare);

   Int_t eventplaneaxis=0;

   if (fRequestEventPlane) eventplaneaxis=1;

   kTrackVariablesPair=6+SetChargeAxis+eventplaneaxis;

   if(fcontainPIDtrig && !fcontainPIDasso) kTrackVariablesPair=7+SetChargeAxis+eventplaneaxis;
 
 if(!fcontainPIDtrig && fcontainPIDasso) kTrackVariablesPair=7+SetChargeAxis+eventplaneaxis;
 
 if(fcontainPIDtrig && fcontainPIDasso) kTrackVariablesPair=8+SetChargeAxis+eventplaneaxis;
 
 
// two particle histograms
  Int_t anaSteps   = 1;       // analysis steps
  const char* title = "d^{2}N_{ch}/d#varphid#eta";

  Int_t iBinPair[kTrackVariablesPair];         // binning for track variables
  Double_t* dBinsPair[kTrackVariablesPair];    // bins for track variables  
  TString* axisTitlePair;  // axis titles for track variables
  axisTitlePair=new TString[kTrackVariablesPair];

 TString defaultBinningStr;
  defaultBinningStr =   "eta: -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\n"
    "p_t_assoc: 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 8.0,10.0\n"
    "p_t_leading_course: 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0,10.0\n"
    "p_t_eff:0.0,0.2,0.3,0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0,5.5, 6.0, 7.0, 8.0,9.0,10.0,12.0,14.0,16.0,20.0\n"
    "vertex: -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10\n"
  "delta_phi: -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389\n" // this binning starts at -pi/2 and is modulo 3 
	"delta_eta: -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,2.1, 2.2, 2.3, 2.4\n"
      "multiplicity: 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1\n"
      "multiplicity_mixing: 0., 1., 2., 3., 4., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1\n";

 if(fRequestEventPlane){
  defaultBinningStr += "eventPlane: -0.5,0.5,1.5,2.5,3.5\n"; // Event Plane Bins (Psi: -0.5->0.5 (in plane), 0.5->1.5 (intermediate), 1.5->2.5 (out of plane), 2.5->3.5 (rest))
  }
 if(fRequestEventPlanemixing){
defaultBinningStr += "eventPlanemixing: 0.0*TMath::DegToRad(), 30.0*TMath::DegToRad(), 60.0*TMath::DegToRad(), 90.0*TMath::DegToRad(), 120.0*TMath::DegToRad(),150.0*TMath::DegToRad(),180.1*TMath::DegToRad()\n";
 }
  if(fcontainPIDtrig){
    /*
      if(fV0TrigCorr){//Invariant mass axis instead of Pid binning axis
	
       	if(fLambda) defaultBinningStr += "InvariantMass:1.065,1.066,1.067,1.068,1.069,1.07,1.071,1.072,1.073,1.074,1.075,1.076,1.077,1.078,1.079,1.08,1.081,1.082,1.083,1.084,1.085,1.086,1.087,1.088,1.089,1.09,1.091,1.092,1.093,1.094,1.095,1.096,1.097,1.098,1.099,1.1,1.101,1.102,1.103,1.104,1.105,1.106,1.107,1.108,1.109,1.11,1.111,1.112,1.113,1.114,1.115,1.116,1.117,1.118,1.119,1.12,1.121,1.122,1.123,1.124,1.125,1.126,1.127,1.128,1.129,1.13,1.131,1.132,1.133,1.134,1.135,1.136,1.137,1.138,1.139,1.14,1.141,1.142,1.143,1.144,1.145,1.146,1.147,1.148,1.149,1.15,1.151,1.152,1.153,1.154,1.155,1.156,1.157,1.158,1.159,1.16,1.161,1.162,1.163,1.164,1.165\n";
	
	
	if(fkShort) defaultBinningStr += "InvariantMass:0.398,0.4,0.402,0.404,0.406,0.408,0.41,0.412,0.414,0.416,0.418,0.42,0.422,0.424,0.426,0.428,0.43,0.432,0.434,0.436,0.438,0.44,0.442,0.444,0.446,0.448,0.45,0.452,0.454,0.456,0.458,0.46,0.462,0.464,0.466,0.468,0.47,0.472,0.474,0.476,0.478,0.48,0.482,0.484,0.486,0.488,0.49,0.492,0.494,0.496,0.498,0.5,0.502,0.504,0.506,0.508,0.51,0.512,0.514,0.516,0.518,0.52,0.522,0.524,0.526,0.528,0.53,0.532,0.534,0.536,0.538,0.54,0.542,0.544,0.546,0.548,0.55,0.552,0.554,0.556,0.558,0.56,0.562,0.564,0.566,0.568,0.57,0.572,0.574,0.576,0.578,0.58,0.582,0.584,0.586,0.588,0.59,0.592,0.594,0.596,0.598\n";
	
  }
    */
    defaultBinningStr += "PIDTrig: -0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5\n"; // course
  }
  
  if(fcontainPIDasso){
      defaultBinningStr += "PIDAsso: -0.5,0.5,1.5,2.5,3.5\n"; // course
  }
 
  if(SetChargeAxis==2){
      defaultBinningStr += "TrigCharge: -2.0,0.0,2.0\n"; // course
      defaultBinningStr += "AssoCharge: -2.0,0.0,2.0\n"; // course
  }

  /*  if(fV0TrigCorr){//default values are PbPb(cent:0,5,10,20,40,60,100)

      defaultBinningStr += "kK0sa0: 3.63508e-03, 3.66158e-03, 3.60326e-03, 3.50492e-03, 3.45697e-3, 3.29805e-03\n"; // course
      defaultBinningStr += "kK0sa1: 6.28389e-04, 5.76389e-04, 5.66960e-04, 5.56158e-04, 5.37069e-4, 5.78176e-04\n"; // course
      defaultBinningStr += "kK0sa2: 6.28389e-04, 5.76389e-04, 5.66960e-04, 5.56158e-04, 5.37069e-4, 5.78176e-04\n"; // course
      defaultBinningStr += "kLambdaa0: 1.27024e-03, 1.26385e-03, 1.27434e-03, 1.27346e-03, 1.25620e-03, 1.26722e-03\n"; // course
      defaultBinningStr += "kLambdaa1: 2.221058e-4,  2.12986e-04, 1.97638e-04, 1.72958e-04, 1.74241e-04, 1.98560e-04\n"; // course
      defaultBinningStr += "kLambdaa2: 2.221058e-4,  2.12986e-04, 1.97638e-04, 1.72958e-04, 1.74241e-04, 1.98560e-04\n"; // course
      defaultBinningStr += "kAntiLambdaa0: 1.27024e-03, 1.26385e-03, 1.27434e-03, 1.27346e-03, 1.25620e-03, 1.26722e-03\n"; // course
      defaultBinningStr += "kAntiLambdaa1: 2.221058e-4,  2.12986e-04, 1.97638e-04, 1.72958e-04, 1.74241e-04, 1.98560e-04\n"; // course
      defaultBinningStr += "kAntiLambdaa2: 2.221058e-4,  2.12986e-04, 1.97638e-04, 1.72958e-04, 1.74241e-04, 1.98560e-04\n"; // course

      }*/
 // =========================================================
  // Customization (adopted from AliUEHistograms)
  // =========================================================

  TObjArray* lines = defaultBinningStr.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    TString tag = line(0, line.Index(":")+1);
    if (!fCustomBinning.BeginsWith(tag) && !fCustomBinning.Contains(TString("\n") + tag))
      fBinningString += line + "\n";
    else
      AliInfo(Form("Using custom binning for %s", tag.Data()));
  }
  delete lines;
  fBinningString += fCustomBinning;
  
  AliInfo(Form("Used AliTHn Binning:\n%s",fBinningString.Data()));

 //  =========================================================
  // Now set the bins
  // =========================================================

    dBinsPair[0]       = GetBinning(fBinningString, "multiplicity", iBinPair[0]);
    axisTitlePair[0]   = "multiplicity";

    dBinsPair[1]     = GetBinning(fBinningString, "vertex", iBinPair[1]);
    axisTitlePair[1]  = "v_{Z} (cm)"; 

    dBinsPair[2]     = GetBinning(fBinningString, "p_t_leading_course", iBinPair[2]);
    axisTitlePair[2]    = "p_{T,trig.} (GeV/c)"; 

    dBinsPair[3]     = GetBinning(fBinningString, "p_t_assoc", iBinPair[3]);
    axisTitlePair[3]    = "p_{T,assoc.} (GeV/c)";

    dBinsPair[4]       = GetBinning(fBinningString, "delta_eta", iBinPair[4]);
    axisTitlePair[4]   = "#Delta#eta"; 

    dBinsPair[5]       = GetBinning(fBinningString, "delta_phi", iBinPair[5]);
    axisTitlePair[5]   = "#Delta#varphi (rad)";  

    Int_t dim_val=6;

    if(fRequestEventPlane){
    dBinsPair[dim_val]       = GetBinning(fBinningString, "eventPlane", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "#varphi - #Psi_{2} (a.u.)";
    dim_val=7;
    }

    if(!fcontainPIDtrig && !fcontainPIDasso && SetChargeAxis==2){
    dBinsPair[dim_val]       = GetBinning(fBinningString, "TrigCharge", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "TrigCharge";

    dBinsPair[dim_val+1]       = GetBinning(fBinningString, "AssoCharge", iBinPair[dim_val+1]);
    axisTitlePair[dim_val+1]   = "AssoCharge";
    }

 if(fcontainPIDtrig && !fcontainPIDasso){
   /*
   if(fV0TrigCorr){//Invariant mass axis instead of Pid binning axis
    dBinsPair[dim_val]       = GetBinning(fBinningString, "InvariantMass", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "InvariantMass"; 
   }
   else{
    dBinsPair[dim_val]       = GetBinning(fBinningString, "PIDTrig", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "PIDTrig"; 
   }
   */
    dBinsPair[dim_val]       = GetBinning(fBinningString, "PIDTrig", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "PIDTrig";
    
    if(SetChargeAxis==2){
    dBinsPair[dim_val+1]       = GetBinning(fBinningString, "TrigCharge", iBinPair[dim_val+1]);
    axisTitlePair[dim_val+1]   = "TrigCharge";

    dBinsPair[dim_val+2]       = GetBinning(fBinningString, "AssoCharge", iBinPair[dim_val+2]);
    axisTitlePair[dim_val+2]   = "AssoCharge";
    }
 }

 if(!fcontainPIDtrig && fcontainPIDasso){
    dBinsPair[dim_val]       = GetBinning(fBinningString, "PIDAsso", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "PIDAsso"; 

 if(SetChargeAxis==2){
    dBinsPair[dim_val+1]       = GetBinning(fBinningString, "TrigCharge", iBinPair[dim_val+1]);
    axisTitlePair[dim_val+1]   = "TrigCharge";

    dBinsPair[dim_val+2]       = GetBinning(fBinningString, "AssoCharge", iBinPair[dim_val+2]);
    axisTitlePair[dim_val+2]   = "AssoCharge";
    }
 }

if(fcontainPIDtrig && fcontainPIDasso){
  /*
   if(fV0TrigCorr){//Invariant mass axis instead of Pid binning axis
    dBinsPair[dim_val]       = GetBinning(fBinningString, "InvariantMass", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "InvariantMass"; 
   }
   else{
    dBinsPair[dim_val]       = GetBinning(fBinningString, "PIDTrig", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "PIDTrig"; 
   }
*/

    dBinsPair[dim_val]       = GetBinning(fBinningString, "PIDTrig", iBinPair[dim_val]);
    axisTitlePair[dim_val]   = "PIDTrig";
    
    dBinsPair[dim_val+1]       = GetBinning(fBinningString, "PIDAsso", iBinPair[dim_val+1]);
    axisTitlePair[dim_val+1]   = "PIDAsso";

    if(SetChargeAxis==2){
    dBinsPair[dim_val+2]       = GetBinning(fBinningString, "TrigCharge", iBinPair[dim_val+2]);
    axisTitlePair[dim_val+2]   = "TrigCharge";

    dBinsPair[dim_val+3]       = GetBinning(fBinningString, "AssoCharge", iBinPair[dim_val+3]);
    axisTitlePair[dim_val+3]   = "AssoCharge";
    }
 }
	
	Int_t nEtaBin = -1;
 	Double_t* EtaBin = GetBinning(fBinningString, "eta", nEtaBin);
	
        Int_t nPteffbin = -1;
 	Double_t* Pteff = GetBinning(fBinningString, "p_t_eff", nPteffbin);

        Int_t multmixbin = -1;
 	Double_t* multmix = GetBinning(fBinningString, "multiplicity_mixing", multmixbin);


	//Set the limits from custom binning
	fminPtTrig=dBinsPair[2][0];
        fmaxPtTrig=dBinsPair[2][iBinPair[2]];
        fminPtAsso=dBinsPair[3][0];
        fmaxPtAsso=dBinsPair[3][iBinPair[3]];
        fmincentmult=dBinsPair[0][0];
        fmaxcentmult=dBinsPair[0][iBinPair[0]];

 //centrality binning for lambda , kshort a0,a1 factor//*********************************hardcoded values in array, be carefu

	/*Int_t Cent_ds  = -1;
        kBinCent_ds=GetBinning(fBinningString, "multiplicity", Cent_ds);
	kNCent_ds=Cent_ds;*/

	/*if(fV0TrigCorr){
	   
        Int_t ka0  = -1;
        kK0s_a0=GetBinning(fBinningString, "kK0sa0", ka0);

	Int_t ka1  = -1;
        kK0s_a1=GetBinning(fBinningString, "kK0sa1", ka1);
	
	Int_t ka2  = -1;
        kK0s_a2=GetBinning(fBinningString, "kK0sa2", ka2);

	Int_t La0  = -1;
        kLambda_a0=GetBinning(fBinningString, "kLambdaa0", La0);

	Int_t La1  = -1;
        kLambda_a1=GetBinning(fBinningString, "kLambdaa1", La1);

        Int_t La2  = -1;
        kLambda_a2=GetBinning(fBinningString, "kLambdaa2", La2);
	
	Int_t ALa0  = -1;
        kAntiLambda_a0=GetBinning(fBinningString, "kAntiLambdaa0", ALa0);

        Int_t ALa1  = -1;
        kAntiLambda_a1=GetBinning(fBinningString, "kAntiLambdaa1", ALa1);

	 Int_t ALa2  = -1;
        kAntiLambda_a2=GetBinning(fBinningString, "kAntiLambdaa2", ALa2);

	 }*/
//****************************************************************************************************************//

	//event pool manager
Int_t MaxNofEvents=1000;
const Int_t NofVrtxBins=10+(1+10)*2;
Double_t ZvrtxBins[NofVrtxBins+1]={ -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10, 
				       90,  92,  94,  96,  98, 100, 102, 104, 106, 108, 110, 
				    190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210}; 


if(fRequestEventPlanemixing){
    // Event plane angle (Psi) bins for event mixing
  
    Int_t nPsiBins=-1; 
    Double_t* psibins = GetBinning(fBinningString, "eventPlanemixing", nPsiBins);
    fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,multmixbin,multmix,NofVrtxBins,ZvrtxBins, nPsiBins, psibins);
    if(psibins)  delete [] psibins; 
				    }

 else{
 const Int_t  nPsiBinsd=1;
 Double_t psibinsd[nPsiBinsd+1]={0.0, 2000.0};
fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,multmixbin,multmix,NofVrtxBins,ZvrtxBins, nPsiBinsd, psibinsd);

 }  
fPoolMgr->SetTargetValues(fMaxNofMixingTracks, 0.1, 5);
 
   if(!fPoolMgr){
      AliError("Event Mixing required, but Pool Manager not initialized...");
      return;
    }

	//fminPtComboeff=fminPtTrig;***then this value will be fixed ,even Setter can't change it's value
	//fmaxPtComboeff=fmaxPtTrig;
//THnSparses for calculation of efficiency

 if((fAnalysisType =="MCAOD") && ffillefficiency) {
TString Histrename;
  Int_t effbin[4];
  effbin[0]=iBinPair[0];
  effbin[1]=iBinPair[1];
  effbin[2]=nPteffbin;
  effbin[3]=nEtaBin;
  Int_t effsteps=5;//for each species type::primMCParticles(0),primRecoTracksMatched(1),allRecoTracksMatched(2),primRecoTracksMatchedPID(3),allRecoTracksMatchedPID(4)
for(Int_t jj=0;jj<6;jj++)//PID type binning
    {
     if(jj==5) effsteps=3;//for unidentified particles
  Histrename="fTrackHistEfficiency";Histrename+=jj;
  fTrackHistEfficiency[jj] = new AliTHn(Histrename.Data(), "Tracking efficiency", effsteps, 4, effbin);
  fTrackHistEfficiency[jj]->SetBinLimits(0, dBinsPair[0]);
  fTrackHistEfficiency[jj]->SetVarTitle(0, "Centrality");
  fTrackHistEfficiency[jj]->SetBinLimits(1, dBinsPair[1]);
  fTrackHistEfficiency[jj]->SetVarTitle(1, "zvtx");
  fTrackHistEfficiency[jj]->SetBinLimits(2, Pteff);
  fTrackHistEfficiency[jj]->SetVarTitle(2, "p_{T} (GeV/c)");
  fTrackHistEfficiency[jj]->SetBinLimits(3, EtaBin);
  fTrackHistEfficiency[jj]->SetVarTitle(3, "#eta");
  fOutput->Add(fTrackHistEfficiency[jj]);
    }
 }

//AliThns for Correlation plots(data &  MC)
 
     if(ffilltrigassoUNID)
       {
    fTHnCorrUNID = new AliTHn("fTHnCorrUNID", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fTHnCorrUNID->SetBinLimits(j, dBinsPair[j]);
    fTHnCorrUNID->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fTHnCorrUNID);

 fTHnCorrUNIDmix = new AliTHn("fTHnCorrUNIDmix", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fTHnCorrUNIDmix->SetBinLimits(j, dBinsPair[j]);
    fTHnCorrUNIDmix->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fTHnCorrUNIDmix);
       }

     if(ffilltrigIDassoID)
       {
fTHnCorrID = new AliTHn("fTHnCorrID", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fTHnCorrID->SetBinLimits(j, dBinsPair[j]);
    fTHnCorrID->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fTHnCorrID);

fTHnCorrIDmix = new AliTHn("fTHnCorrIDmix", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fTHnCorrIDmix->SetBinLimits(j, dBinsPair[j]);
    fTHnCorrIDmix->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fTHnCorrIDmix);
       }

     if(ffilltrigUNIDassoID || ffilltrigIDassoUNID)//***********a bit tricky, be careful
       {
fTHnCorrIDUNID = new AliTHn("fTHnCorrIDUNID", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fTHnCorrIDUNID->SetBinLimits(j, dBinsPair[j]);
    fTHnCorrIDUNID->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fTHnCorrIDUNID);


fTHnCorrIDUNIDmix = new AliTHn("fTHnCorrIDUNIDmix", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fTHnCorrIDUNIDmix->SetBinLimits(j, dBinsPair[j]);
    fTHnCorrIDUNIDmix->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fTHnCorrIDUNIDmix);
       }



  //ThnSparse for Correlation plots(truth MC)
     if(ffilltrigIDassoIDMCTRUTH) {//remember that in this case uidentified means other than pions, kaons, protons

fCorrelatonTruthPrimary = new AliTHn("fCorrelatonTruthPrimary", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fCorrelatonTruthPrimary->SetBinLimits(j, dBinsPair[j]);
    fCorrelatonTruthPrimary->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fCorrelatonTruthPrimary);


fCorrelatonTruthPrimarymix = new AliTHn("fCorrelatonTruthPrimarymix", title, anaSteps, kTrackVariablesPair, iBinPair);
for (Int_t j=0; j<kTrackVariablesPair; j++) {
    fCorrelatonTruthPrimarymix->SetBinLimits(j, dBinsPair[j]);
    fCorrelatonTruthPrimarymix->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutput->Add(fCorrelatonTruthPrimarymix);     
 }

    //binning for trigger no. counting

     Int_t ChargeAxis=0;
     if(SetChargeAxis==2) ChargeAxis=1;

	Int_t* fBinst;
	Int_t dims=3+ChargeAxis+eventplaneaxis;
	if(fcontainPIDtrig) dims=4+ChargeAxis+eventplaneaxis;
        fBinst= new Int_t[dims];
   Double_t* dBinsTrig[dims];    // bins for track variables  
   TString* axisTitleTrig;  // axis titles for track variables
   axisTitleTrig=new TString[dims];

	for(Int_t i=0; i<3;i++)
	  {
	    fBinst[i]=iBinPair[i];
	    dBinsTrig[i]=dBinsPair[i];
	    axisTitleTrig[i]=axisTitlePair[i];
	  }
	Int_t dim_val_trig=3;
    if(fRequestEventPlane){
      fBinst[dim_val_trig]=iBinPair[6];//if fRequestEventPlane=TRUE, dim_val already becomes 7.
      dBinsTrig[dim_val_trig]=dBinsPair[6];
      axisTitleTrig[dim_val_trig]=axisTitlePair[6];
      dim_val_trig=4;
    }

if(!fcontainPIDtrig && !fcontainPIDasso && ChargeAxis==1){
fBinst[dim_val_trig]=iBinPair[dim_val];
dBinsTrig[dim_val_trig]=dBinsPair[dim_val];
axisTitleTrig[dim_val_trig]=axisTitlePair[dim_val];
    }

if(fcontainPIDtrig && !fcontainPIDasso){
fBinst[dim_val_trig]=iBinPair[dim_val];
dBinsTrig[dim_val_trig]=dBinsPair[dim_val];
axisTitleTrig[dim_val_trig]=axisTitlePair[dim_val]; 
    if(ChargeAxis==1){
fBinst[dim_val_trig+1]=iBinPair[dim_val+1];
dBinsTrig[dim_val_trig+1]=dBinsPair[dim_val+1];
axisTitleTrig[dim_val_trig+1]=axisTitlePair[dim_val+1];
    }
 }

 if(!fcontainPIDtrig && fcontainPIDasso){
 if(ChargeAxis==1){
    fBinst[dim_val_trig]=iBinPair[dim_val+1];
dBinsTrig[dim_val_trig]=dBinsPair[dim_val+1];
axisTitleTrig[dim_val_trig]=axisTitlePair[dim_val+1];
    }
 }

if(fcontainPIDtrig && fcontainPIDasso){
  fBinst[dim_val_trig]=iBinPair[dim_val];
dBinsTrig[dim_val_trig]=dBinsPair[dim_val];
axisTitleTrig[dim_val_trig]=axisTitlePair[dim_val]; 
    if(ChargeAxis==1){
fBinst[dim_val_trig+1]=iBinPair[dim_val+2];
dBinsTrig[dim_val_trig+1]=dBinsPair[dim_val+2];
axisTitleTrig[dim_val_trig+1]=axisTitlePair[dim_val+2];
    }
    }
 
  //AliTHns for trigger counting(data & reco MC)
  if(ffilltrigassoUNID || ffilltrigUNIDassoID || ffilltrigIDassoUNID || ffilltrigIDassoID)
	  {
	    fTHnTrigcount = new  AliTHn("fTHnTrigcount", "fTHnTrigcount", 2, dims, fBinst); //2 steps;;;;0->same event;;;;;1->mixed event
   for(Int_t i=0; i<dims;i++){
    fTHnTrigcount->SetBinLimits(i, dBinsTrig[i]);
    fTHnTrigcount->SetVarTitle(i, axisTitleTrig[i]);
  } 
  fOutput->Add(fTHnTrigcount);
	  }
  
  if(ffilltrigIDassoIDMCTRUTH) {
  //AliTHns for trigger counting(truth MC)
  fTHnTrigcountMCTruthPrim = new  AliTHn("fTHnTrigcountMCTruthPrim", "fTHnTrigcountMCTruthPrim", 2, dims, fBinst); //2 steps;;;;0->same event;;;;;1->mixed event
 for(Int_t i=0; i<dims;i++){
    fTHnTrigcountMCTruthPrim->SetBinLimits(i, dBinsTrig[i]);
    fTHnTrigcountMCTruthPrim->SetVarTitle(i, axisTitleTrig[i]);
  } 
  fOutput->Add(fTHnTrigcountMCTruthPrim);
 }

if(fAnalysisType=="MCAOD" || fAnalysisType=="MC"){
  if(ffillhistQATruth)
    {
  MCtruthpt=new TH1F ("MCtruthpt","ptdistributiontruthprim",200,0.0,fmaxPt);
  fOutputList->Add(MCtruthpt);

  MCtrutheta=new TH1F ("MCtrutheta","etadistributiontruthprim",360,-1.8,1.8);
  fOutputList->Add(MCtrutheta);

  MCtruthphi=new TH1F ("MCtruthphi","phidisttruthprim",340,0,6.8);
  fOutputList->Add(MCtruthphi);

  MCtruthpionpt=new TH1F ("MCtruthpionpt","MCtruthpionpt",200,0.0,fmaxPt);
  fOutputList->Add(MCtruthpionpt);

  MCtruthpioneta=new TH1F ("MCtruthpioneta","MCtruthpioneta",360,-1.8,1.8);
  fOutputList->Add(MCtruthpioneta);

  MCtruthpionphi=new TH1F ("MCtruthpionphi","MCtruthpionphi",340,0,6.8);
  fOutputList->Add(MCtruthpionphi);

  MCtruthkaonpt=new TH1F ("MCtruthkaonpt","MCtruthkaonpt",200,0.0,fmaxPt);
  fOutputList->Add(MCtruthkaonpt);

  MCtruthkaoneta=new TH1F ("MCtruthkaoneta","MCtruthkaoneta",360,-1.8,1.8);
  fOutputList->Add(MCtruthkaoneta);

  MCtruthkaonphi=new TH1F ("MCtruthkaonphi","MCtruthkaonphi",340,0,6.8);
  fOutputList->Add(MCtruthkaonphi);

  MCtruthprotonpt=new TH1F ("MCtruthprotonpt","MCtruthprotonpt",200,0.0,fmaxPt);
  fOutputList->Add(MCtruthprotonpt);

  MCtruthprotoneta=new TH1F ("MCtruthprotoneta","MCtruthprotoneta",360,-1.8,1.8);
  fOutputList->Add(MCtruthprotoneta);

  MCtruthprotonphi=new TH1F ("MCtruthprotonphi","MCtruthprotonphi",340,0,6.8);
  fOutputList->Add(MCtruthprotonphi);
    }
 fPioncont=new TH2F("fPioncont", "fPioncont",10,-0.5,9.5,200,0.0,fmaxPt);
  fOutputList->Add(fPioncont);

 fKaoncont=new TH2F("fKaoncont","fKaoncont",10,-0.5,9.5,200,0.0,fmaxPt);
  fOutputList->Add(fKaoncont);

 fProtoncont=new TH2F("fProtoncont","fProtoncont",10,-0.5,9.5,200,0.0,fmaxPt);
  fOutputList->Add(fProtoncont);

fUNIDcont=new TH2F("fUNIDcont","fUNIDcont",10,-0.5,9.5,200,0.0,fmaxPt);
  fOutputList->Add(fUNIDcont);
  }

fEventno=new TH2F("fEventno","fEventno",iBinPair[0], dBinsPair[0],iBinPair[1],dBinsPair[1]);
 fEventno->GetXaxis()->SetTitle("Centrality");
 fEventno->GetYaxis()->SetTitle("Z_Vtx");
fOutput->Add(fEventno);
fEventnobaryon=new TH2F("fEventnobaryon","fEventnobaryon",iBinPair[0], dBinsPair[0],iBinPair[1],dBinsPair[1]);
 fEventnobaryon->GetXaxis()->SetTitle("Centrality");
 fEventnobaryon->GetYaxis()->SetTitle("Z_Vtx");
fOutput->Add(fEventnobaryon);
fEventnomeson=new TH2F("fEventnomeson","fEventnomeson",iBinPair[0], dBinsPair[0],iBinPair[1],dBinsPair[1]);
 fEventnomeson->GetXaxis()->SetTitle("Centrality");
 fEventnomeson->GetYaxis()->SetTitle("Z_Vtx");
fOutput->Add(fEventnomeson);

fhistJetTrigestimate=new TH2F("fhistJetTrigestimate","fhistJetTrigestimate",iBinPair[0],dBinsPair[0],6,-0.5,5.5);
fOutput->Add(fhistJetTrigestimate);

   fTwoTrackDistancePtdip = new TH3F("fTwoTrackDistancePtdip", ";#Delta#eta;#Delta#varphi;#Delta p_{T}", 36, -1.8, 1.8, 72,-TMath::Pi()/2, 3*TMath::Pi()/2, 40, 0, 10);
  fOutput->Add(fTwoTrackDistancePtdip);

fTwoTrackDistancePtdipmix = new TH3F("fTwoTrackDistancePtdipmix", ";#Delta#eta;#Delta#varphi;#Delta p_{T}", 36, -1.8, 1.8, 72,-TMath::Pi()/2, 3*TMath::Pi()/2, 40, 0, 10);
  fOutput->Add(fTwoTrackDistancePtdipmix);

  TString Histttrname;
for(Int_t jj=0;jj<2;jj++)// PID type binning
    {
  Histttrname="fTwoTrackDistancePt";Histttrname+=jj;
  fTwoTrackDistancePt[jj] = new TH3F(Histttrname.Data(), ";#Delta#eta;#Delta#varphi^{*}_{min};#Delta p_{T}", 100, -0.15, 0.15, 100, -0.05, 0.05, 20, 0, 10);
  fOutput->Add(fTwoTrackDistancePt[jj]);

 Histttrname="fTwoTrackDistancePtmix";Histttrname+=jj;
  fTwoTrackDistancePtmix[jj] = new TH3F(Histttrname.Data(), ";#Delta#eta;#Delta#varphi^{*}_{min};#Delta p_{T}", 100, -0.15, 0.15, 100, -0.05, 0.05, 20, 0, 10);
  fOutput->Add(fTwoTrackDistancePtmix[jj]);
    }
//Mixing
//DefineEventPool();


  if(fV0TrigCorr){

    // TFile *fsifile1= TFile::Open(ffilenamesigmaV0);

     
    if (TString(ffilenamesigmaV0).BeginsWith("alien:"))
    TGrid::Connect("alien:");
    TFile *fsifile1=TFile::Open(ffilenamesigmaV0);
 
kShortSigmahisto = (TH2F*)fsifile1->Get("v0MassSigmaCent_KS0");
LambdaSigmahisto = (TH2F*)fsifile1->Get("v0MassSigmaCent_LAM");

kShortSigmahisto->SetDirectory(0);
LambdaSigmahisto->SetDirectory(0);

 fsifile1->Close();

  }

  if(fapplyTrigefficiency || fapplyAssoefficiency)
   {
     const Int_t nDimt = 4;//       cent zvtx  pt   eta
     Int_t fBinsCht[nDimt] = {iBinPair[0], iBinPair[1], nPteffbin ,nEtaBin};//*************change it
     Double_t fMinCht[nDimt] = { dBinsPair[0][0],dBinsPair[1][0], Pteff[0], EtaBin[0] };
     Double_t fMaxCht[nDimt] = {dBinsPair[0][iBinPair[0]], dBinsPair[1][iBinPair[1]], Pteff[nPteffbin], EtaBin[nEtaBin]};

  TString Histrexname;
for(Int_t jj=0;jj<6;jj++)// PID type binning
    {
  Histrexname="effcorection";Histrexname+=jj;
  effcorection[jj] = new THnSparseF(Histrexname.Data(),"cent:zvtx::Pt:eta", nDimt, fBinsCht, fMinCht, fMaxCht);
  effcorection[jj]->Sumw2(); 
  effcorection[jj]->GetAxis(0)->Set(iBinPair[0], dBinsPair[0]);
  effcorection[jj]->GetAxis(0)->SetTitle("Centrality");
  effcorection[jj]->GetAxis(1)->Set( iBinPair[1],dBinsPair[1]);
  effcorection[jj]->GetAxis(1)->SetTitle("zvtx"); 
  effcorection[jj]->GetAxis(2)->Set(nPteffbin, Pteff);  
  effcorection[jj]->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  effcorection[jj]->GetAxis(3)->Set( nEtaBin,EtaBin);
  effcorection[jj]->GetAxis(3)->SetTitle("#eta");
  fOutput->Add(effcorection[jj]);
    }
// TFile *fsifile = new TFile(fefffilename,"READ");

 if (TString(fefffilename).BeginsWith("alien:"))
    TGrid::Connect("alien:");
 TFile *fileT=TFile::Open(fefffilename);
 TString Nameg;
for(Int_t jj=0;jj<6;jj++)//type binning
    {
Nameg="effmap";Nameg+=jj;
//effcorection[jj] = (THnSparseF*)fsifile->Get(Nameg.Data());
effcorection[jj] = (THnSparseF*)fileT->Get(Nameg.Data());

//effcorection[jj]->SetDirectory(0);//****************************not present in case oh THnF
    }
//fsifile->Close();
fileT->Close();

   }

 delete [] EtaBin; 
 delete [] Pteff; 
 delete [] multmix; 

 //******************************************************************V0 plots*********************************************//
 if(fV0TrigCorr){
	//histos for v0
	//Double_t BinsV0[]={1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.5,8.0};
	//const Int_t nbinsV0 =sizeof(BinsV0)/sizeof(Double_t)-1;

 
   if(ffillofflineV0){
   fHistRawPtCentInvK0s= new TH3F("fHistRawPtCentInvK0s", "K^{0}_{s}: mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",100,0.398,0.598,200,0.0,fmaxPt,100,0.0,100.);
fOutput->Add(fHistRawPtCentInvK0s);


 fHistRawPtCentInvLambda= new TH3F("fHistRawPtCentInvLambda", "#Lambda: mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",100,1.065,1.165,200,0.0,fmaxPt,100,0.0,100.);
fOutput->Add(fHistRawPtCentInvLambda);


 fHistRawPtCentInvAntiLambda= new TH3F("fHistRawPtCentInvAntiLambda", "#bar{#Lambda} : mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",100,1.065,1.165,200,0.0,fmaxPt,100,0.0,100.);
fOutput->Add(fHistRawPtCentInvAntiLambda);


 fHistFinalPtCentInvK0s= new TH3F("fHistFinalPtCentInvK0s", "K^{0}_{s}: mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",100,0.398,0.598,200,0.0,fmaxPt,100,0.0,100.);
fOutput->Add(fHistFinalPtCentInvK0s);


fHistFinalPtCentInvLambda= new TH3F("fHistFinalPtCentInvLambda", "#Lambda: mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",100,1.065,1.165,200,0.0,fmaxPt,100,0.0,100.);
fOutput->Add(fHistFinalPtCentInvLambda);


 fHistFinalPtCentInvAntiLambda= new TH3F("fHistFinalPtCentInvAntiLambda", "#bar{#Lambda} : mass vs #it{p}_{T};Mass (GeV/#it{c}^2);#it{p}_{T} (GeV/#it{c});centrality",100,1.065,1.165,200,0.0,fmaxPt,100,0.0,100.);
fOutput->Add(fHistFinalPtCentInvAntiLambda);
   }
 }

  //*************************************************************EP plots***********************************************//
  if(fRequestEventPlane){
  // TProfile for resolutions 3 subevents (V0A, V0C, TPC)
  // v2
  fHResTPCv0A2 = new TProfile("hResTPCv0A2","",nCentrBin,0,nCentrBin);
  fHResTPCv0C2 = new TProfile("hResTPCv0C2","",nCentrBin,0,nCentrBin);
  fHResv0Cv0A2 = new TProfile("hResv0Cv0A2","",nCentrBin,0,nCentrBin);

  fList->Add(fHResTPCv0A2);
  fList->Add(fHResTPCv0C2);
  fList->Add(fHResv0Cv0A2);

  // v3
  fHResTPCv0A3 = new TProfile("hResTPCv0A3","",nCentrBin,0,nCentrBin);
  fHResTPCv0C3 = new TProfile("hResTPCv0C3","",nCentrBin,0,nCentrBin);
  fHResv0Cv0A3 = new TProfile("hResv0Cv0A3","",nCentrBin,0,nCentrBin);

  fList->Add(fHResTPCv0A3);
  fList->Add(fHResTPCv0C3);
  fList->Add(fHResv0Cv0A3);

  // MC as in the dataEP resolution (but using MC tracks)
  if(fAnalysisType == "MCAOD"  && fV2){
    fHResMA2 = new TProfile("hResMA2","",nCentrBin,0,nCentrBin);
    fHResMC2 = new TProfile("hResMC2","",nCentrBin,0,nCentrBin);
    fHResAC2 = new TProfile("hResAC2","",nCentrBin,0,nCentrBin);
    fList->Add(fHResMA2); 
    fList->Add(fHResMC2); 
    fList->Add(fHResAC2); 
  }
  if(fAnalysisType == "MCAOD" && fV3){
    fHResMA3 = new TProfile("hResMA3","",nCentrBin,0,nCentrBin);
    fHResMC3 = new TProfile("hResMC3","",nCentrBin,0,nCentrBin);
    fHResAC3 = new TProfile("hResAC3","",nCentrBin,0,nCentrBin);
    fList->Add(fHResMA3); 
    fList->Add(fHResMC3); 
    fList->Add(fHResAC3); 
  }


  // V0A and V0C event plane distributions
  //v2 
  fPhiRPTPC = new TH2F("fPhiRPTPCv2","#phi distribution of EP TPC;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/2,TMath::Pi()/2);
  fPhiRPTPCv3 = new TH2F("fPhiRPTPCv3","#phi distribution of EP TPC;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/3,TMath::Pi()/3);
  fList->Add(fPhiRPTPC);
  fList->Add(fPhiRPTPCv3);

  fPhiRPv0A = new TH2F("fPhiRPv0Av2","#phi distribution of EP VZERO-A;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/2,TMath::Pi()/2);
  fPhiRPv0C = new TH2F("fPhiRPv0Cv2","#phi distribution of EP VZERO-C;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/2,TMath::Pi()/2);
  fList->Add(fPhiRPv0A);
  fList->Add(fPhiRPv0C);

  //v3
  fPhiRPv0Av3 = new TH2F("fPhiRPv0Av3","#phi distribution of EP VZERO-A;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/3,TMath::Pi()/3);
  fPhiRPv0Cv3 = new TH2F("fPhiRPv0Cv3","#phi distribution of EP VZERO-C;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/3,TMath::Pi()/3);
  fList->Add(fPhiRPv0Av3);
  fList->Add(fPhiRPv0Cv3);

  fHistEventPlaneTruth = new TH2F("fHistEventPlaneTruth","#phi distribution of EP MCTRUTHheader;centrality;#phi (rad)",nCentrBin,0,nCentrBin,nPsiTOF,-TMath::Pi()/2,TMath::Pi()/2);
  fList->Add(fHistEventPlaneTruth);

  }
    
  //*****************************************************PIDQA histos*****************************************************//

// for electron rejection only TPC nsigma histograms
  if(fAnalysisType=="MCAOD" || fAnalysisType=="AOD"){
  if(fElectronRejection) {
 
    fHistdEdxVsPTPCbeforePIDelectron = new TH2D ("dEdxVsPTPCbeforeelectron","dEdxVsPTPCbeforeelectron", 1000, -10.0, 10.0, 1000, 0, 1000); 
    fOutputList->Add(fHistdEdxVsPTPCbeforePIDelectron);
    
    fHistNSigmaTPCvsPtbeforePIDelectron = new TH2D ("NSigmaTPCvsPtbeforeelectron","NSigmaTPCvsPtbeforeelectron", 1000, -10, 10, 1000, 0, 500); 
    fOutputList->Add(fHistNSigmaTPCvsPtbeforePIDelectron);
    
    fHistdEdxVsPTPCafterPIDelectron = new TH2D ("dEdxVsPTPCafterelectron","dEdxVsPTPCafterelectron", 1000, -10, 10, 1000, 0, 1000); 
    fOutputList->Add(fHistdEdxVsPTPCafterPIDelectron);

    fHistNSigmaTPCvsPtafterPIDelectron = new TH2D ("NSigmaTPCvsPtafterelectron","NSigmaTPCvsPtafterelectron", 1000, -10, 10, 1000, 0, 500); 
    fOutputList->Add(fHistNSigmaTPCvsPtafterPIDelectron); 
  }
}
 
  //nsigma plot
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-30;
      Double_t maxy=30;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=50;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigma_%d_%d",ipart,ipid),Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0.0,fmaxPt,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //nsigmaRec plot
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-15;
      Double_t maxy=15;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaRec_%d_%d",ipart,ipid),
				  Form("n#sigma for reconstructed %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0.0,fmaxPt,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }

  //BayesRec plot
  if(fPIDType==Bayes){//use bayesianPID
    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();//****************************************Need to know about it

  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    Double_t miny=0.;
    Double_t maxy=1;
    TH2F *fHistoBayes=new TH2F(Form("BayesRec_%d",ipart),
			       Form("probability for reconstructed %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,500,miny,maxy);
    fHistoBayes->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayes->GetYaxis()->SetTitle(Form("Bayes prob %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayes);


   TH2F *fHistoBayesTPC=new TH2F(Form("probBayes_TPC_%d",ipart),
			       Form("probability for Tracks as %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,500,miny,maxy);
    fHistoBayesTPC->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayesTPC->GetYaxis()->SetTitle(Form("Bayes prob TPC %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayesTPC);

  TH2F *fHistoBayesTOF=new TH2F(Form("probBayes_TOF_%d",ipart),
			       Form("probability for Tracks as %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,500,miny,maxy);
    fHistoBayesTOF->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayesTOF->GetYaxis()->SetTitle(Form("Bayes prob TOF %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayesTOF);

 TH2F *fHistoBayesTPCTOF=new TH2F(Form("probBayes_TPCTOF_%d",ipart),
			       Form("probability for Tracks as  %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,500,miny,maxy);
    fHistoBayesTPCTOF->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayesTPCTOF->GetYaxis()->SetTitle(Form("Bayes prob TPCTOF %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayesTPCTOF);
  }
  }

  //nsigma separation power plot 
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
 Double_t miny=0;
 Double_t maxy=10;
   TH2F *Pi_Ka_sep=new TH2F(Form("Pi_Ka_sep_%d",ipid),
			       Form("Pi_Ka separation in %s",kPIDTypeName[ipid]),200,0.0,fmaxPt,200,miny,maxy);
    Pi_Ka_sep->GetXaxis()->SetTitle("P_{T} (GeV/C)");
    Pi_Ka_sep->GetYaxis()->SetTitle(Form("expected seaparation(n#sigma) in %s",kPIDTypeName[ipid]));
    fOutputList->Add(Pi_Ka_sep);

   TH2F *Pi_Pr_sep=new TH2F(Form("Pi_Pr_sep_%d",ipid),
			       Form("Pi_Pr separation in %s",kPIDTypeName[ipid]),200,0.0,fmaxPt,200,miny,maxy);
    Pi_Pr_sep->GetXaxis()->SetTitle("P_{T} (GeV/C)");
    Pi_Pr_sep->GetYaxis()->SetTitle(Form("expected seaparation(n#sigma) in %s",kPIDTypeName[ipid]));
    fOutputList->Add(Pi_Pr_sep);

    TH2F *Ka_Pr_sep=new TH2F(Form("Ka_Pr_sep_%d",ipid),
			       Form("Ka_Pr separation in %s",kPIDTypeName[ipid]),200,0.0,fmaxPt,200,miny,maxy);
    Ka_Pr_sep->GetXaxis()->SetTitle("P_{T} (GeV/C)");
    Ka_Pr_sep->GetYaxis()->SetTitle(Form("expected seaparation(n#sigma) in %s",kPIDTypeName[ipid]));
    fOutputList->Add(Ka_Pr_sep);
    }
    //deltapion histos
    //before deltapion cut
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
      Double_t miny=-40;
      Double_t maxy=20;
      if(ipart!=SpPion) continue;//only around pion's mean position;proton and pion's delta in one histo
      TH2F *fHistodelta=new TH2F(Form("deltapion_%d",ipart),
				  Form("deltapion %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,600,miny,maxy);
      fHistodelta->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistodelta->GetYaxis()->SetTitle(Form("deltapion %s",kParticleSpeciesName[ipart]));
      fOutputList->Add(fHistodelta);
  }

     //After deltapion cut
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
      Double_t miny=-40;
      Double_t maxy=20;
      if(ipart!=SpPion) continue;//only around pion's mean position;protn and pion's delta in one histo
      TH2F *fHistodeltaRec=new TH2F(Form("deltapionRec_%d",ipart),
				  Form("deltapionRec %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,600,miny,maxy);
      fHistodeltaRec->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistodeltaRec->GetYaxis()->SetTitle(Form("deltapionRec %s",kParticleSpeciesName[ipart]));
      fOutputList->Add(fHistodeltaRec);
  }

    //TRUTHMC PID deltapion
    if (fAnalysisType == "MCAOD"){
      for(Int_t ipart=0;ipart<NSpecies;ipart++){
      Double_t miny=-40;
      Double_t maxy=20;
      TH2F *fHistodeltaMC=new TH2F(Form("deltapionMC_%d",ipart),
				  Form("deltapionMC %s",kParticleSpeciesName[ipart]),200,0.0,fmaxPt,600,miny,maxy);
      fHistodeltaMC->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistodeltaMC->GetYaxis()->SetTitle(Form("deltapionMC %s",kParticleSpeciesName[ipart]));
      fOutputList->Add(fHistodeltaMC);
  }
}

  //nsigmaDC plot
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaDC_%d_%d",ipart,ipid),
				  Form("n#sigma for double counting %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0.0,fmaxPt,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //nsigmaMC plot
 if (fAnalysisType == "MCAOD"){
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-30;
      Double_t maxy=30;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=50;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaMC_%d_%d",ipart,ipid),
				  Form("n#sigma for MC %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0.0,fmaxPt,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  }
  //PID signal plot
  for(Int_t idet=0;idet<fNDetectors;idet++){
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
      Double_t maxy=500;
      if(idet==fTOF)maxy=1.1;
      TH2F *fHistoPID=new TH2F(Form("PID_%d_%d",idet,ipart),Form("%s signal - %s",kDetectorName[idet],kParticleSpeciesName[ipart]),200,0.0,fmaxPt,500,-maxy,maxy);
      fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
      fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kDetectorName[idet]));
      fOutputList->Add(fHistoPID);
    }
  }
  //PID signal plot, before PID cut
  for(Int_t idet=0;idet<fNDetectors;idet++){
    Double_t maxy=500;
    if(idet==fTOF)maxy=1.1;
    TH2F *fHistoPID=new TH2F(Form("PIDAll_%d",idet),Form("%s signal",kDetectorName[idet]),200,0.0,fmaxPt,500,-maxy,maxy);
    fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
    fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kDetectorName[idet]));
    fOutputList->Add(fHistoPID);
  }

  PostData(1, fOutput);              // Post data for ALL output slots >0 here, to get at least an empty histogram
  PostData(2, fOutputList);
  if(fRequestEventPlane) PostData(3, fList);
  AliInfo("Finished setting up the Output");

   TH1::AddDirectory(oldStatus);
}
//-------------------------------------------------------------------------------
void AliTwoParticlePIDCorr::UserExec( Option_t * ){
 
  if(fAnalysisType == "AOD") {
    
    doAODevent();
    
  }//AOD--analysis-----

  else if(fAnalysisType == "MCAOD" || fAnalysisType == "MC") {
  
    doMCAODevent();

  }
  
  else return;
  
}
//-------------------------------------------------------------------------
void AliTwoParticlePIDCorr::doMCAODevent() 
{

  // get the event (for generator level: MCEvent())
  AliVEvent* event = NULL;

  if(fAnalysisType == "MC")  event = dynamic_cast<AliVEvent*>(MCEvent());

  else{
    event = dynamic_cast<AliVEvent*>(InputEvent());     
  }
  if(!event) {
    AliError("eventMain not available");
    return;
  }

  Double_t Inv_mass=0.0;//has no meaning for pions, kaons and protons(just set 0.0) to fill the LRCParticlePID position
    evplaneMC=999.;
    fgPsi2v0aMC=999.;
    fgPsi2v0cMC=999.;
    fgPsi2tpcMC=999.;
    fgPsi3v0aMC=999.;
    fgPsi3v0cMC=999.;
    fgPsi3tpcMC=999.;
    gReactionPlane = 999.;

 // get centrality object and check quality(valid for p-Pb and Pb-Pb; coming soon for pp 7 TeV)
  Double_t cent_v0=-1.0;
  Double_t effcent=1.0;
  Double_t refmultReco =0.0;
   Double_t nooftrackstruth=0.0;//in case of pp this will give the multiplicity(for truth case) after the track loop(only for unidentified particles that pass  kinematic cuts)


if(fAnalysisType == "MC"){
   
    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(event);

 if(!gMCEvent) {
      AliError("mcEvent not available");
      return ;
    }
// count all events(physics triggered)   
  fEventCounter->Fill(1);

	AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(gMCEvent->GenEventHeader());
	if(!header) return;  
	  TArrayF gVertexArray;
	  header->PrimaryVertex(gVertexArray);
          Float_t zVtxmc =gVertexArray.At(2);
	  //cout<<"*****************************************************************************************************hi I am here"<<endl;

 
	  cent_v0=GetAcceptedEventMultiplicity((AliVEvent*)gMCEvent,kFALSE); //b value; 2nd argument has no meaning
 
 if(cent_v0<0.) return;//mainly returns impact parameter
//within proper centrality range having positive value upto the maximum range mentioned in the multiplicity binning of AddTask macro 
  fEventCounter->Fill(13);

 //get the event plane in case of PbPb
   if(fRequestEventPlane){
     gReactionPlane=GetEventPlane((AliVEvent*)gMCEvent,kTRUE,cent_v0);//get the truth event plane,middle argument has no meaning in this case
   if(gReactionPlane==999.) return;
 }

   /*   
   TObjArray* tracksMCtruth_t=new TObjArray;//for truth MC particles with PID,here unidentified means any particle other than pion, kaon or proton(Basicaly Spundefined of AliHelperPID)******WARNING::different from data and reco MC
 tracksMCtruth_t->SetOwner(kTRUE);  
   */

TObjArray* tracksMCtruth=new TObjArray;//for truth MC particles with PID,here unidentified means any particle other than pion, kaon or proton(Basicaly Spundefined of AliHelperPID)******WARNING::different from data and reco MC
 tracksMCtruth->SetOwner(kTRUE);  //***********************************IMPORTANT!

for (Int_t iTracks = 0; iTracks < gMCEvent->GetNumberOfPrimaries(); iTracks++) {
	AliMCParticle* partMC = dynamic_cast<AliMCParticle *>(gMCEvent->GetTrack(iTracks));
	if (!partMC) {
	  AliError(Form("Could not receive particle %d", iTracks));
	  continue;
	}
//exclude non stable particles
	if(fselectprimaryTruth && !(gMCEvent->IsPhysicalPrimary(iTracks))) continue;

//consider only charged particles
    if(partMC->Charge() == 0) continue;


//give only kinematic cuts at the generator level  
 if (partMC->Eta() < fmineta || partMC->Eta() > fmaxeta) continue;
 if (partMC->Pt() < fminPt ||  partMC->Pt() > fmaxPt) continue;

 if(!partMC) continue;//for safety

          TParticle *particle = partMC->Particle();
	  if(!particle) continue;
           Int_t particletypeTruth=-999;
	  
	  Int_t pdgtruth = particle->GetPdgCode();

	  //electron rejection(not necessary if we consider physical primary tracks-basiccaly to avoid gamma conversions)
	  if(fElectronRejectionTRUTH){
	    if (TMath::Abs(pdgtruth)==11) continue;
	  }

 //To determine multiplicity in case of PP
 nooftrackstruth++;
 //cout<<"**************************************"<<TMath::Abs(partMC->GetLabel())<<endl;
//only physical primary(all/unidentified)  
if(ffillhistQATruth)
    {
 MCtruthpt->Fill(partMC->Pt());
 MCtrutheta->Fill(partMC->Eta());
 MCtruthphi->Fill(partMC->Phi());
    }
 if (TMath::Abs(pdgtruth)==211)
   {
 particletypeTruth=SpPion;
if(ffillhistQATruth)
    {
 MCtruthpionpt->Fill(partMC->Pt());
 MCtruthpioneta->Fill(partMC->Eta());
 MCtruthpionphi->Fill(partMC->Phi());
    }
      }
 if (TMath::Abs(pdgtruth)==321)
   {
 particletypeTruth=SpKaon;
if(ffillhistQATruth)
    {
 MCtruthkaonpt->Fill(partMC->Pt());
 MCtruthkaoneta->Fill(partMC->Eta());
 MCtruthkaonphi->Fill(partMC->Phi());
  }
    }
if(TMath::Abs(pdgtruth)==2212)
  {
 particletypeTruth=SpProton;
if(ffillhistQATruth)
    {
 MCtruthprotonpt->Fill(partMC->Pt());
 MCtruthprotoneta->Fill(partMC->Eta());
 MCtruthprotonphi->Fill(partMC->Phi());
    }
     }
 if(TMath::Abs(pdgtruth)!=211 && TMath::Abs(pdgtruth)!=321 && TMath::Abs(pdgtruth)!=2212)  particletypeTruth=unidentified;//*********************WARNING:: situation is different from reco MC and data case(here we don't have SpUndefined particles,because here unidentified=SpUndefined)

    if(fRequestEventPlane){
      FillPIDEventPlane(cent_v0,particletypeTruth,partMC->Phi(),gReactionPlane);
    }
    /*
//Exclude resonances
	if(fExcludeResonancesInMC) {
	  TParticle *particle = track->Particle();
	  if(!particle) continue;
	  
	  Bool_t kExcludeParticle = kFALSE;
	  Int_t gMotherIndex = particle->GetFirstMother();
	  if(gMotherIndex != -1) {
	    AliMCParticle* motherTrack = dynamic_cast<AliMCParticle *>(event->GetTrack(gMotherIndex));
	    if(motherTrack) {
	      TParticle *motherParticle = motherTrack->Particle();
	      if(motherParticle) {
		Int_t pdgCodeOfMother = motherParticle->GetPdgCode();
		//if((pdgCodeOfMother == 113)||(pdgCodeOfMother == 213)||(pdgCodeOfMother == 221)||(pdgCodeOfMother == 223)||(pdgCodeOfMother == 331)||(pdgCodeOfMother == 333)) {
		}
		if(pdgCodeOfMother == 113  // rho0
		   || pdgCodeOfMother == 213 || pdgCodeOfMother == -213 // rho+
		   // || pdgCodeOfMother == 221  // eta
		   // || pdgCodeOfMother == 331  // eta'
		   // || pdgCodeOfMother == 223  // omega
		   // || pdgCodeOfMother == 333  // phi
		   || pdgCodeOfMother == 311  || pdgCodeOfMother == -311 // K0
		   // || pdgCodeOfMother == 313  || pdgCodeOfMother == -313 // K0*
		   // || pdgCodeOfMother == 323  || pdgCodeOfMother == -323 // K+*
		   || pdgCodeOfMother == 3122 || pdgCodeOfMother == -3122 // Lambda
		   || pdgCodeOfMother == 111  // pi0 Dalitz
		   ) {
		  kExcludeParticle = kTRUE;
		}
	      }
	    }
	  }
	  
	  //Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	  if(kExcludeParticle) continue;
	}

	//Exclude electrons with PDG
	if(fExcludeElectronsInMC) {
	  
	  TParticle *particle = track->Particle();
	  
	  if (particle){ 
	    if(TMath::Abs(particle->GetPdgCode()) == 11) continue;
	  }
	}
    */

 Float_t effmatrixtruth=1.0;//In Truth MC, no case of efficiency correction so it should be always 1.0
if((partMC->Pt()>=fminPtAsso && partMC->Pt()<=fmaxPtAsso) || (partMC->Pt()>=fminPtTrig && partMC->Pt()<=fmaxPtTrig))//to reduce memory consumption in pool
  {
    Short_t chargeval=0;
    if(partMC->Charge()>0)   chargeval=1;
    if(partMC->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
    const TBits *clustermap=0;
    const TBits *sharemap=0;
    LRCParticlePID* copy6 = new LRCParticlePID(particletypeTruth,Inv_mass,chargeval,partMC->Pt(),partMC->Eta(), partMC->Phi(),effmatrixtruth,clustermap,sharemap);
//copy6->SetUniqueID(eventno * 100000 + TMath::Abs(partMC->GetLabel()));
 copy6->SetUniqueID(eventno * 100000 + (Int_t)nooftrackstruth);
 tracksMCtruth->Add(copy6);//************** TObjArray used for truth correlation function calculation
  }
 }//track loop ends

 if(nooftrackstruth>0.0){
if (fSampleType=="pPb" || fSampleType=="PbPb" || fPPVsMult==kTRUE || fCentralityMethod == "MC_b") fCentralityCorrelationMC->Fill(cent_v0, nooftrackstruth);//only with unidentified tracks(i.e before PID selection);;;;;can be used to remove centrality outliers??????
/*
AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());

 if(collGeometry){
Float_t NpartProj= collGeometry-> ProjectileParticipants();
Float_t NpartTarg = collGeometry->TargetParticipants();
 Float_t Npart= (NpartProj + NpartTarg);
fNchNpartCorr->Fill(Npart,nooftrackstruth);//Creating some problem while being written in the TList in outputcontainer:::NOT SOLVED(ONLY in case of"MC", Fine for "MCAOD")
 }
*/
 }
/*
if(fRunShufflingbalance){
    tracksMCtruth = GetShuffledTracks(tracksMCtruth_t);
  }
 else tracksMCtruth=CloneAndReduceTrackList(tracksMCtruth_t);
*/
 
  if (fRandomizeReactionPlane)//only for TRuth MC??
  {
    Double_t centralityDigits = cent_v0*1000. - (Int_t)(cent_v0*1000.);
    Double_t angle = TMath::TwoPi() * centralityDigits;
    AliInfo(Form("Shifting phi of all tracks by %f (digits %f)", angle, centralityDigits));
    ShiftTracks(tracksMCtruth, angle);  
  }
 

 Float_t weghtval=1.0;
 Float_t bSign = 0;

 if(nooftrackstruth>0.0) 
  {
//no. of events analyzed
  fEventCounter->Fill(15);
 if(ffilltrigIDassoIDMCTRUTH){
 //Fill Correlations for MC truth particles(same event)
if(tracksMCtruth && tracksMCtruth->GetEntriesFast()>0)//hadron triggered correlation
  Fillcorrelation(gReactionPlane,tracksMCtruth,0,cent_v0,zVtxmc,weghtval,kFALSE,bSign,fPtOrderMCTruth,kFALSE,kFALSE,"trigIDassoIDMCTRUTH");//mixcase=kFALSE for same event case

//start mixing
 AliEventPool* pool2 = fPoolMgr->GetEventPool(cent_v0, zVtxmc+200, gReactionPlane);
if (pool2 && pool2->IsReady())
  {//start mixing only when pool->IsReady
if(tracksMCtruth && tracksMCtruth->GetEntriesFast()>0)
  {//proceed only when no. of trigger particles >0 in current event
    Float_t nmix=(Float_t)pool2->GetCurrentNEvents();  
for (Int_t jMix=0; jMix<pool2->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks6 = pool2->GetEvent(jMix);
  if(!bgTracks6) continue;
  Fillcorrelation(gReactionPlane,tracksMCtruth,bgTracks6,cent_v0,zVtxmc,nmix,(jMix == 0),bSign,fPtOrderMCTruth,kFALSE,kTRUE,"trigIDassoIDMCTRUTH");//mixcase=kTRUE for mixing case
  
   }// pool event loop ends mixing case
 }//if(trackstrig && trackstrig->GetEntriesFast()>0) condition ends mixing case
} //if pool->IsReady() condition ends mixing case

 //still in main event loop

 if(tracksMCtruth){
if(pool2)  pool2->UpdatePool(CloneAndReduceTrackList(tracksMCtruth));//ownership of tracksasso is with pool now, don't delete it
 }
  }
  }
 //still in main event loop

 //if(tracksMCtruth_t) delete tracksMCtruth_t;

if(tracksMCtruth) delete tracksMCtruth;

 }//MC type

//"MC" type analysis is finished but still in event loop

 else{//if(fAnalysisType=="MCAOD")

  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }
 
// count all events(physics triggered)   
  fEventCounter->Fill(1);

   

//check the PIDResponse handler
     fPID = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPID) AliFatal("This Task needs the PID response attached to the inputHandler");
    //if (!fPID) return;
// get mag. field required for twotrack efficiency cut
 Float_t bSign = 0;
 bSign = (aod->GetMagneticField() > 0) ? 1 : -1;

 //check for TClonesArray(truth track MC information)
fArrayMC = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fArrayMC) {
    AliFatal("Error: MC particles branch not found!\n");
    return;
  }
  
  //check for AliAODMCHeader(truth event MC information)
  AliAODMCHeader *header=NULL;
  header=(AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
  if(!header) {
    printf("MC header branch not found!\n");
    return;
  }
 
//Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
Float_t zVtxmc =header->GetVtxZ();
 if(TMath::Abs(zVtxmc)>fzvtxcut) return;

 // For productions with injected signals, figure out above which label to skip particles/tracks

 if (fInjectedSignals)
  {
    AliGenEventHeader* eventHeader = 0;
    Int_t headers = 0;

// AOD
      if (!header)
      AliFatal("fInjectedSignals set but no MC header found");
      
      headers = header->GetNCocktailHeaders();
      eventHeader = header->GetCocktailHeader(0);

 if (!eventHeader)
    {
      // We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
      // (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
      AliError("First event header not found. Skipping this event.");
      //fHistos->FillEvent(centrality, AliUEHist::kCFStepAnaTopology);
      return;
    }
skipParticlesAbove = eventHeader->NProduced();
    AliInfo(Form("Injected signals in this event (%d headers). Keeping events of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove));
  }

 if (fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMult==kFALSE))
   {
 //make the event selection with reco vertex cut and centrality cut and return the value of the centrality
     Double_t refmultTruth = GetAcceptedEventMultiplicity((AliVEvent*)aod,kTRUE);  //incase of ref multiplicity it will return the truth MC ref mullt value; need to determine the ref mult value separately for reco Mc case; in case of centrality this is final and fine
     refmultReco = GetAcceptedEventMultiplicity((AliVEvent*)aod,kFALSE); 
     if(refmultTruth<=0 || refmultReco<=0) return;
     cent_v0=refmultTruth;
   }
 else {
 cent_v0=GetAcceptedEventMultiplicity((AliVEvent*)aod,kFALSE); //centrality value; 2nd argument has no meaning
 }


 if(cent_v0<0.) return;
 effcent=cent_v0;// This will be required for efficiency THn filling(specially in case of pp)


//count selected events having centrality betn 0-100%
 fEventCounter->Fill(13);

  //get the event plane in case of PbPb
   if(fRequestEventPlane){
     gReactionPlane=GetEventPlane((AliVEvent*)aod,kTRUE,cent_v0);//get the truth event plane
   if(gReactionPlane==999.) return;
 }
   /*
TObjArray* tracksMCtruth_t=new TObjArray;//for truth MC particles with PID,here unidentified means any particle other than pion, kaon or proton(Basicaly Spundefined of AliHelperPID)******WARNING::different from data and reco MC
 tracksMCtruth_t->SetOwner(kTRUE);  
   */
   
TObjArray* tracksMCtruth=new TObjArray;//for truth MC particles with PID,here unidentified means any particle other than pion, kaon or proton(Basicaly Spundefined of AliHelperPID)******WARNING::different from data and reco MC
 tracksMCtruth->SetOwner(kTRUE);  //***********************************IMPORTANT!

  eventno++;

  //There is a small difference between zvtx and zVtxmc?????? 
  //cout<<"***********************************************zVtxmc="<<zVtxmc<<endl;
  //cout<<"***********************************************zvtx="<<zvtx<<endl;
 
//now process the truth particles(for both efficiency & correlation function)
Int_t nMCTrack = fArrayMC->GetEntriesFast();
  
for (Int_t iMC = 0; iMC < nMCTrack; iMC++) 
{      //MC truth track loop starts
    
AliAODMCParticle *partMC = (AliAODMCParticle*) fArrayMC->At(iMC);
    
    if(!partMC){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",iMC));
      continue;
    }

//consider only charged particles
    if(partMC->Charge() == 0) continue;

//consider only primary particles; neglect all secondary particles including from weak decays
 if(fselectprimaryTruth && !partMC->IsPhysicalPrimary()) continue;


//remove injected signals(primaries above <maxLabel>)
 if (fInjectedSignals && partMC->GetLabel() >= skipParticlesAbove) continue;

//remove duplicates
  Bool_t isduplicate=kFALSE;
 if (fRemoveDuplicates)
   { 
 for (Int_t j=iMC+1; j<nMCTrack; ++j) 
   {//2nd trutuh loop starts
AliAODMCParticle *partMC2 = (AliAODMCParticle*) fArrayMC->At(j);
   if(!partMC2){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",j));
      continue;
    }    
 if (partMC->GetLabel() == partMC2->GetLabel())
   {
isduplicate=kTRUE;
 break;  
   }    
   }//2nd truth loop ends
   }
 if(fRemoveDuplicates && isduplicate) continue;//remove duplicates

//give only kinematic cuts at the truth level  
 if (partMC->Eta() < fmineta || partMC->Eta() > fmaxeta) continue;
 if (partMC->Pt() < fminPt ||  partMC->Pt() > fmaxPt) continue;

 if(!partMC) continue;//for safety

  //get particle ID
Int_t pdgtruth=((AliAODMCParticle*)partMC)->GetPdgCode();
Int_t particletypeTruth=-999;

   //electron rejection(not necessary if we consider physical primary tracks-basiccaly to avoid gamma conversions)
   if(fElectronRejectionTRUTH){
	    if (TMath::Abs(pdgtruth)==11) continue;
	  }
   
 //To determine multiplicity in case of PP
 nooftrackstruth++;
 //cout<<"**************************************"<<TMath::Abs(partMC->GetLabel())<<endl;
//only physical primary(all/unidentified)  
if(ffillhistQATruth)
    {
 MCtruthpt->Fill(partMC->Pt());
 MCtrutheta->Fill(partMC->Eta());
 MCtruthphi->Fill(partMC->Phi());
    }

 if (TMath::Abs(pdgtruth)==211)
   {
 particletypeTruth=SpPion;
if(ffillhistQATruth)
    {

 MCtruthpionpt->Fill(partMC->Pt());
 MCtruthpioneta->Fill(partMC->Eta());
 MCtruthpionphi->Fill(partMC->Phi());
    }
      }
 if (TMath::Abs(pdgtruth)==321)
   {
 particletypeTruth=SpKaon;
if(ffillhistQATruth)
    {
 MCtruthkaonpt->Fill(partMC->Pt());
 MCtruthkaoneta->Fill(partMC->Eta());
 MCtruthkaonphi->Fill(partMC->Phi());
  }
    }
if(TMath::Abs(pdgtruth)==2212)
  {
 particletypeTruth=SpProton;
if(ffillhistQATruth)
    {
 MCtruthprotonpt->Fill(partMC->Pt());
 MCtruthprotoneta->Fill(partMC->Eta());
 MCtruthprotonphi->Fill(partMC->Phi());
    }
     }
 if(TMath::Abs(pdgtruth)!=211 && TMath::Abs(pdgtruth)!=321 && TMath::Abs(pdgtruth)!=2212)  particletypeTruth=unidentified;//*********************WARNING:: situation is different from reco MC and data case(here we don't have SpUndefined particles,because here unidentified=SpUndefined)

    if(fRequestEventPlane){
      FillPIDEventPlane(cent_v0,particletypeTruth,partMC->Phi(),gReactionPlane);
    }

 // -- Fill THnSparse for efficiency and contamination calculation
    if (fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMult==kFALSE)) effcent=15.0;//integrated over multiplicity(so put any fixed value for each track so that practically means there is only one bin in multiplicity i.e. multiplicity intregated out )**************Important

 Double_t primmctruth[4] = {effcent, zVtxmc,partMC->Pt(), partMC->Eta()};
 if(ffillefficiency)
  {
    fTrackHistEfficiency[5]->Fill(primmctruth,0);//for all primary truth particles(4)
    if (TMath::Abs(pdgtruth)==211 || TMath::Abs(pdgtruth)==321) fTrackHistEfficiency[3]->Fill(primmctruth,0);//for  primary truth mesons(3)
    if (TMath::Abs(pdgtruth)==2212 || TMath::Abs(pdgtruth)==321) fTrackHistEfficiency[4]->Fill(primmctruth,0);//for  primary truth kaons+protons(4)
    if (TMath::Abs(pdgtruth)==211)  fTrackHistEfficiency[0]->Fill(primmctruth,0);//for pions
    if (TMath::Abs(pdgtruth)==321)  fTrackHistEfficiency[1]->Fill(primmctruth,0);//for kaons
    if (TMath::Abs(pdgtruth)==2212)  fTrackHistEfficiency[2]->Fill(primmctruth,0);//for protons
  }
   
 Float_t effmatrixtruth=1.0;//In Truth MC, no case of efficiency correction so it should be always 1.0
if((partMC->Pt()>=fminPtAsso && partMC->Pt()<=fmaxPtAsso) || (partMC->Pt()>=fminPtTrig && partMC->Pt()<=fmaxPtTrig))//to reduce memory consumption in pool
  {
    Short_t chargeval=0;
    if(partMC->Charge()>0)   chargeval=1;
    if(partMC->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
    const TBits *clustermap=0;
    const TBits *sharemap=0;
    LRCParticlePID* copy6 = new LRCParticlePID(particletypeTruth,Inv_mass,chargeval,partMC->Pt(),partMC->Eta(), partMC->Phi(),effmatrixtruth,clustermap,sharemap);
//copy6->SetUniqueID(eventno * 100000 + TMath::Abs(partMC->GetLabel()));
 copy6->SetUniqueID(eventno * 100000 + (Int_t)nooftrackstruth);
 tracksMCtruth->Add(copy6);//************** TObjArray used for truth correlation function calculation
  }
  }//MC truth track loop ends

//*********************still in event loop

/*
 if(fRunShufflingbalance){
    tracksMCtruth = GetShuffledTracks(tracksMCtruth_t);
  }
 else tracksMCtruth=CloneAndReduceTrackList(tracksMCtruth_t);
*/
 
 
   if (fRandomizeReactionPlane)//only for TRuth MC??
  {
    Double_t centralityDigits = cent_v0*1000. - (Int_t)(cent_v0*1000.);
    Double_t angle = TMath::TwoPi() * centralityDigits;
    AliInfo(Form("Shifting phi of all tracks by %f (digits %f)", angle, centralityDigits));
    ShiftTracks(tracksMCtruth, angle);  
  }


   
if(nooftrackstruth>0.0){
if (fSampleType=="pPb" || fSampleType=="PbPb" || fPPVsMult==kTRUE || fCentralityMethod == "MC_b") fCentralityCorrelationMC->Fill(cent_v0, nooftrackstruth);//only with unidentified tracks(i.e before PID selection);;;;;can be used to remove centrality outliers??????
 
 if(fSampleType=="pPb" || fSampleType=="PbPb"){ //GetCocktailHeader(0) is creating problem for Phojet(AOD60), anyhow for pp it has no use  
   AliGenEventHeader* eventHeader = header->GetCocktailHeader(0);  // get first MC header from either ESD/AOD (including cocktail header if available)
      if (eventHeader)
      {	     
AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*>(eventHeader);
 if(collGeometry){
Float_t NpartProj= collGeometry-> ProjectileParticipants();
Float_t NpartTarg = collGeometry->TargetParticipants();
Float_t Npart= (NpartProj + NpartTarg);
fNchNpartCorr->Fill(Npart,nooftrackstruth);
         }
      }
 }
 }
 
 Float_t weghtval=1.0;

if(nooftrackstruth>0.0 && ffilltrigIDassoIDMCTRUTH)
  {
 //Fill Correlations for MC truth particles(same event)
if(tracksMCtruth && tracksMCtruth->GetEntriesFast()>0)//hadron triggered correlation
  Fillcorrelation(gReactionPlane,tracksMCtruth,0,cent_v0,zVtxmc,weghtval,kFALSE,bSign,fPtOrderMCTruth,kFALSE,kFALSE,"trigIDassoIDMCTRUTH");//mixcase=kFALSE for same event case

//start mixing
 AliEventPool* pool2 = fPoolMgr->GetEventPool(cent_v0, zVtxmc+200, gReactionPlane);
if (pool2 && pool2->IsReady())
  {//start mixing only when pool->IsReady
if(tracksMCtruth && tracksMCtruth->GetEntriesFast()>0)
  {//proceed only when no. of trigger particles >0 in current event
    Float_t nmix=(Float_t)pool2->GetCurrentNEvents();  
for (Int_t jMix=0; jMix<pool2->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks6 = pool2->GetEvent(jMix);
  if(!bgTracks6) continue;
  Fillcorrelation(gReactionPlane,tracksMCtruth,bgTracks6,cent_v0,zVtxmc,nmix,(jMix == 0),bSign,fPtOrderMCTruth,kFALSE,kTRUE,"trigIDassoIDMCTRUTH");//mixcase=kTRUE for mixing case
  
   }// pool event loop ends mixing case
 }//if(trackstrig && trackstrig->GetEntriesFast()>0) condition ends mixing case
} //if pool->IsReady() condition ends mixing case

 //still in main event loop

 if(tracksMCtruth){
if(pool2)  pool2->UpdatePool(CloneAndReduceTrackList(tracksMCtruth));//ownership of tracksasso is with pool now, don't delete it
 }
  }

 //still in main event loop

//if(tracksMCtruth_t) delete tracksMCtruth_t;

if(tracksMCtruth) delete tracksMCtruth;

//now deal with reco tracks


//Float_t bSign1=aod->GetHeader()->GetMagneticField() ;//used for reconstructed track dca cut
   Float_t bSign1=aod->GetMagneticField() ;//used for reconstructed track dca cut


//detrmine the ref mult in case of Reco(not required if we get centrality info from AliCentrality)
 if (fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMult==kFALSE)) cent_v0=refmultReco;
 effcent=cent_v0;// This will be required for efficiency THn filling(specially in case of pp)

 if(fRequestEventPlane){
   gReactionPlane = GetEventPlane((AliVEvent*)aod,kFALSE,cent_v0);//get the reconstructed event plane
    if(gReactionPlane==999.) return;
 }


  TExMap *trackMap = new TExMap();

// --- track loop for mapping matrix
 if(fFilterBit==128)
   {
 for (Int_t itrk = 0; itrk < aod->GetNumberOfTracks(); itrk++) 
{ //track loop starts for TObjArray(containing track and event information) filling; used for correlation function calculation 
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrk));
  if (!track) continue;
  Int_t tracktype=ClassifyTrack(track,trkVtx,bSign1,kFALSE);//don't fill the histos here
  if(tracktype!=1) continue; 

  if(!track) continue;//for safety

   Int_t gid = track->GetID();
   trackMap->Add(gid,itrk);
 }//track looop ends
   }


 /*
   TObjArray* tracksUNID_t = new TObjArray;
   tracksUNID_t->SetOwner(kTRUE);   
 */
 
   TObjArray* tracksUNID = new TObjArray;
   tracksUNID->SetOwner(kTRUE);

   TObjArray* tracksID = new TObjArray;
   tracksID->SetOwner(kTRUE);


//get the selected v0 particle TObjArray
   TObjArray* tracksIDV0 = new TObjArray;
   tracksIDV0->SetOwner(kTRUE);

   Double_t trackscount=0.0;
// loop over reconstructed tracks 
  for (Int_t itrk = 0; itrk < aod->GetNumberOfTracks(); itrk++) 
{ // reconstructed track loop starts
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrk));
  if (!track) continue;
 //get the corresponding MC track at the truth level (doing reco matching)
  AliAODMCParticle* recomatched = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(track->GetLabel())));
  if(!recomatched) continue;//if a reco track doesn't have corresponding truth track at generated level is a fake track(label==0), ignore it

//remove injected signals 
 if(fInjectedSignals)
   {
    AliAODMCParticle* mother = recomatched;

      while (!mother->IsPhysicalPrimary())
      {// find the primary mother;the first stable mother is searched and checked if it is <= <maxLabel>
	if (mother->GetMother() < 0)
	{
	  mother = 0;
	  break;
	}
	  
   mother =(AliAODMCParticle*) fArrayMC->At(((AliAODMCParticle*)mother)->GetMother());
	if (!mother)
	  break;
      }
 if (!mother)
    {
      Printf("WARNING: No mother found for particle %d:", recomatched->GetLabel());
      continue;
    }
 if (mother->GetLabel() >= skipParticlesAbove) continue;//remove injected signals(primaries above <maxLabel>)
   }//remove injected signals

 if (fRemoveWeakDecays && ((AliAODMCParticle*) recomatched)->IsSecondaryFromWeakDecay()) continue;//remove weak decays
	
  Bool_t isduplicate2=kFALSE;
if (fRemoveDuplicates)
   {
  for (Int_t j =itrk+1; j < aod->GetNumberOfTracks(); j++) 
    {//2nd loop starts
 AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(aod->GetTrack(j));
 if (!track2) continue;
 AliAODMCParticle* recomatched2 = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(track2->GetLabel())));
if(!recomatched2) continue;

if (track->GetLabel() == track2->GetLabel())
   {
isduplicate2=kTRUE;
 break;  
   }
    }//2nd loop ends
   }
 if(fRemoveDuplicates && isduplicate2) continue;//remove duplicates
     
  fHistQA[11]->Fill(track->GetTPCNcls());
  Int_t tracktype=ClassifyTrack(track,trkVtx,bSign1,kTRUE);//dcacut=kFALSE,onlyprimary=kFALSE

 if(tracktype==0) continue; 
 if(tracktype==1)//tracks "not" passed AliAODTrack::kPrimary at reconstructed level & have proper TPC PID response(?)
{
  if(!track) continue;//for safety
 //accepted all(primaries+secondary) reconstructed tracks(pt 0.2 to 10.0,,eta -0.8 to 0.8)

  if (fElectronRejection){
	//Fill QA after the PID
	fHistdEdxVsPTPCafterPIDelectron -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
	fHistNSigmaTPCvsPtafterPIDelectron -> Fill(track->P()*track->Charge(),TMath::Abs(fPID->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron))); }

  
 if(fV0TrigCorr){ 
if(IsTrackFromV0(aod,track)) continue;// remove auto correlation
  }

  AliAODTrack *PIDtrack=track;//for PID purpose, mainly important for TPC only tracks

  if(fFilterBit==128){
Int_t gid1 = track->GetID();
//if(gid1>=0) PIDtrack = track;
 PIDtrack =(AliAODTrack*) aod->GetTrack(trackMap->GetValue(-1-gid1));
if(!PIDtrack) continue;//for safety; so that each of the TPC only tracks have corresponding global track along with it
  }

  trackscount++;

//check for eta , phi holes
 fEtaSpectrasso->Fill(track->Eta(),track->Pt());
 fphiSpectraasso->Fill(track->Phi(),track->Pt());

  Int_t particletypeMC=-9999;

//tag all particles as unidentified
 particletypeMC=unidentified;

 Float_t effmatrix=1.;

// -- Fill THnSparse for efficiency calculation
 if (fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMult==kFALSE)) effcent=15.0;//integrated over multiplicity(so put any fixed value for each track so that practically means there is only one bin in multiplicity i.e. multiplicity intregated out )**************Important
 //NOTE:: this will be used for fillinfg THnSparse of efficiency & also to get the the track by track eff. factor on the fly(only in case of pp)

 //Clone & Reduce track list(TObjArray) for unidentified particles
if((track->Pt()>=fminPtAsso && track->Pt()<=fmaxPtAsso) || (track->Pt()>=fminPtTrig && track->Pt()<=fmaxPtTrig))//to reduce memory consumption in pool
  {
     Short_t chargeval=0;
    if(track->Charge()>0)   chargeval=1;
    if(track->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
 if (fapplyTrigefficiency || fapplyAssoefficiency)//get the trackingefficiency x contamination factor for unidentified particles
   effmatrix=GetTrackbyTrackeffvalue(track,effcent,zvtx,particletypeMC);
 LRCParticlePID* copy = new LRCParticlePID(particletypeMC,Inv_mass,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix,track->GetTPCClusterMapPtr(),track->GetTPCSharedMapPtr());
   copy->SetUniqueID(eventno * 100000 +(Int_t)trackscount);
   tracksUNID->Add(copy);//track information Storage for UNID correlation function(tracks that pass the filterbit & kinematic cuts only)
  }

/*
   if(fRunShufflingbalance){
    tracksUNID = GetShuffledTracks(tracksUNID_t);
  }
   else tracksUNID=CloneAndReduceTrackList(tracksUNID_t);
*/

//get the pdg code of the corresponding truth particle
 Int_t pdgCode = ((AliAODMCParticle*)recomatched)->GetPdgCode();

 Double_t allrecomatchedpid[4] = {effcent, zVtxmc,recomatched->Pt(), recomatched->Eta()};
 if(ffillefficiency) {
fTrackHistEfficiency[5]->Fill(allrecomatchedpid,2);//for allreco matched
 if(TMath::Abs(pdgCode)==211 ||  TMath::Abs(pdgCode)==321)   fTrackHistEfficiency[3]->Fill(allrecomatchedpid,2);//for mesons
 if(TMath::Abs(pdgCode)==321 ||  TMath::Abs(pdgCode)==2212)   fTrackHistEfficiency[4]->Fill(allrecomatchedpid,2);//for kaons+protons
 if(TMath::Abs(pdgCode)==211)  fTrackHistEfficiency[0]->Fill(allrecomatchedpid,2);//for pions  
 if(TMath::Abs(pdgCode)==321)  fTrackHistEfficiency[1]->Fill(allrecomatchedpid,2);//for kaons
 if(TMath::Abs(pdgCode)==2212) fTrackHistEfficiency[2]->Fill(allrecomatchedpid,2);//for protons

 if (((AliAODMCParticle*)recomatched)->IsPhysicalPrimary()) {
 fTrackHistEfficiency[5]->Fill(allrecomatchedpid,1);//for primreco matched
 if(TMath::Abs(pdgCode)==211 ||  TMath::Abs(pdgCode)==321)   fTrackHistEfficiency[3]->Fill(allrecomatchedpid,1);//for mesons
 if(TMath::Abs(pdgCode)==321 ||  TMath::Abs(pdgCode)==2212)   fTrackHistEfficiency[4]->Fill(allrecomatchedpid,1);//for kaons+protons
 if( TMath::Abs(pdgCode)==211)  fTrackHistEfficiency[0]->Fill(allrecomatchedpid,1);//for pions  
 if( TMath::Abs(pdgCode)==321)  fTrackHistEfficiency[1]->Fill(allrecomatchedpid,1);//for kaons
 if( TMath::Abs(pdgCode)==2212) fTrackHistEfficiency[2]->Fill(allrecomatchedpid,1);//for protons
 }
 }

 //now start the particle identification process:)

// DCA XY cut to check contaminations from lambda using dca sigma cut variation with filterbit 16, only for identified particles 
     if (fDCAXYCut)
	{
	  if (!trkVtx) continue;
	  
	  Double_t pos[2];
	  Double_t covar[3];
	  AliAODTrack* clone =(AliAODTrack*) track->Clone();
	  Bool_t success = clone->PropagateToDCA(trkVtx, bSign1, fdcacutvalue, pos, covar);
	  delete clone;
	  if (!success)
	    continue;

// 	  Printf("%f", ((AliAODTrack*)part)->DCA());
// 	  Printf("%f", pos[0]);
	  if (TMath::Abs(pos[0]) > fDCAXYCut->Eval(track->Pt())) continue;
	}
 
Float_t dEdx = PIDtrack->GetTPCsignal();
 fHistoTPCdEdx->Fill(track->Pt(), dEdx);

 if(HasTOFPID(PIDtrack))
{
Double_t beta = GetBeta(PIDtrack);
fHistoTOFbeta->Fill(track->Pt(), beta);
 }

 //remove the tracks which don't have proper TOF response-otherwise the misIDentification rate values will be wrong
if(fRequestTOFPID && track->Pt()>fPtTOFPIDmin && track->Pt()<fPtTOFPIDmax && (!HasTOFPID(PIDtrack)) ) continue;

 //above fPtTOFPIDmax PID is done using rel. rise in TPC only
 if(track->Pt()>=fPtTOFPIDmax){
   if (!HasTPCPID(track)) continue;
   Bool_t TPCsectoredge=TPCCutMIGeo(track,aod);
   if(TPCSectoredgecut) {if(!TPCsectoredge) continue;}
 }

//do track identification(nsigma method)
  particletypeMC=GetParticle(PIDtrack,fFIllPIDQAHistos);//******************************problem is here

  //Plot NSIGMA distr of MC particles with thei PDG code (TRUTH PID)
switch(TMath::Abs(pdgCode)){
  case 2212:
    if(fFIllPIDQAHistos){
      if(track->Pt()<fPtTOFPIDmax){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(PIDtrack)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",SpProton,ipid));
	h->Fill(track->Pt(),fnsigmas[SpProton][ipid]);
      }
    }
else{
        TH2F *h1=GetHistogram2D(Form("NSigmaMC_%d_%d",SpProton,NSigmaTPC));
	h1->Fill(track->Pt(),fnsigmas[SpProton][NSigmaTPC]);
  
	TH2F *h=GetHistogram2D(Form("deltapionMC_%d",SpProton));
  	h->Fill(track->Pt(),deltapion_val);
 }      
    }
    break;
  case 321:
    if(fFIllPIDQAHistos){
      if(track->Pt()<fPtTOFPIDmax){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(PIDtrack)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",SpKaon,ipid));
	h->Fill(track->Pt(),fnsigmas[SpKaon][ipid]);
      }
   }
else{
        TH2F *h2=GetHistogram2D(Form("NSigmaMC_%d_%d",SpKaon,NSigmaTPC));
	h2->Fill(track->Pt(),fnsigmas[SpKaon][NSigmaTPC]);
	
	TH2F *h=GetHistogram2D(Form("deltapionMC_%d",SpKaon));
  	h->Fill(track->Pt(),deltapion_val);
 }  
    }
    break;
  case 211:
    if(fFIllPIDQAHistos){
            if(track->Pt()<fPtTOFPIDmax){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(PIDtrack)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",SpPion,ipid));
	h->Fill(track->Pt(),fnsigmas[SpPion][ipid]);
      }
   }
 else{

        TH2F *h3=GetHistogram2D(Form("NSigmaMC_%d_%d",SpPion,NSigmaTPC));
	h3->Fill(track->Pt(),fnsigmas[SpPion][NSigmaTPC]);
	
	TH2F *h=GetHistogram2D(Form("deltapionMC_%d",SpPion));
  	h->Fill(track->Pt(),deltapion_val);
 } 
    }
    break;
  }
  
//2-d TPCTOF map(for each Pt interval)
  if(HasTOFPID(PIDtrack)){
 fTPCTOFPion3d->Fill(track->Pt(),fnsigmas[SpPion][NSigmaTOF],fnsigmas[SpPion][NSigmaTPC]);
 fTPCTOFKaon3d->Fill(track->Pt(),fnsigmas[SpKaon][NSigmaTOF],fnsigmas[SpKaon][NSigmaTPC]);
 fTPCTOFProton3d->Fill(track->Pt(),fnsigmas[SpProton][NSigmaTOF],fnsigmas[SpProton][NSigmaTPC]); 
  }

 //Pt, Eta , Phi distribution of the reconstructed identified particles
if(ffillhistQAReco)
    {
if (particletypeMC==SpPion)
  {
    fPionPt->Fill(track->Pt());
    fPionEta->Fill(track->Eta());
    fPionPhi->Fill(track->Phi());
  }
if (particletypeMC==SpKaon)
  {
    fKaonPt->Fill(track->Pt());
    fKaonEta->Fill(track->Eta());
    fKaonPhi->Fill(track->Phi());
  }
if (particletypeMC==SpProton)
  {
    fProtonPt->Fill(track->Pt());
    fProtonEta->Fill(track->Eta());
    fProtonPhi->Fill(track->Phi());
  }
    }

 //fill tracking efficiency
 if(ffillefficiency)
   {
 if(particletypeMC==SpPion || particletypeMC==SpKaon)
   {
     if(TMath::Abs(pdgCode)==211 ||  TMath::Abs(pdgCode)==321) {
       fTrackHistEfficiency[3]->Fill(allrecomatchedpid,4);//for mesons
 if (((AliAODMCParticle*)recomatched)->IsPhysicalPrimary())  fTrackHistEfficiency[3]->Fill(allrecomatchedpid,3);//for mesons
     }
   }
 if(particletypeMC==SpKaon || particletypeMC==SpProton)
   {
     if(TMath::Abs(pdgCode)==321 ||  TMath::Abs(pdgCode)==2212) {
       fTrackHistEfficiency[4]->Fill(allrecomatchedpid,4);//for kaons+protons
 if (((AliAODMCParticle*)recomatched)->IsPhysicalPrimary())  fTrackHistEfficiency[4]->Fill(allrecomatchedpid,3);
     }
   }
 if(particletypeMC==SpPion && TMath::Abs(pdgCode)==211)  {
   fTrackHistEfficiency[0]->Fill(allrecomatchedpid,4);//for pions
 if (((AliAODMCParticle*)recomatched)->IsPhysicalPrimary())  fTrackHistEfficiency[0]->Fill(allrecomatchedpid,3);
 } 
 if(particletypeMC==SpKaon && TMath::Abs(pdgCode)==321) {
   fTrackHistEfficiency[1]->Fill(allrecomatchedpid,4);//for kaons
if (((AliAODMCParticle*)recomatched)->IsPhysicalPrimary())  fTrackHistEfficiency[1]->Fill(allrecomatchedpid,3);
 }
 if(particletypeMC==SpProton && TMath::Abs(pdgCode)==2212){
   fTrackHistEfficiency[2]->Fill(allrecomatchedpid,4);//for protons
if (((AliAODMCParticle*)recomatched)->IsPhysicalPrimary())  fTrackHistEfficiency[2]->Fill(allrecomatchedpid,3);
 }
   }


 //for misidentification fraction calculation(do it with tuneonPID)
 if(particletypeMC==SpPion )
   {
     if(TMath::Abs(pdgCode)==211) fPioncont->Fill(1.,track->Pt());
     if(TMath::Abs(pdgCode)==321) fPioncont->Fill(3.,track->Pt());
     if(TMath::Abs(pdgCode)==2212) fPioncont->Fill(5.,track->Pt());
     if(TMath::Abs(pdgCode)!=211 && TMath::Abs(pdgCode)!=321 && TMath::Abs(pdgCode)!=2212) fPioncont->Fill(7.,track->Pt());
   }
if(particletypeMC==SpKaon )
   {
     if(TMath::Abs(pdgCode)==211) fKaoncont->Fill(1.,track->Pt());
     if(TMath::Abs(pdgCode)==321) fKaoncont->Fill(3.,track->Pt());
     if(TMath::Abs(pdgCode)==2212) fKaoncont->Fill(5.,track->Pt());
     if(TMath::Abs(pdgCode)!=211 && TMath::Abs(pdgCode)!=321 && TMath::Abs(pdgCode)!=2212) fKaoncont->Fill(7.,track->Pt());
   }
 if(particletypeMC==SpProton )
   {
     if(TMath::Abs(pdgCode)==211) fProtoncont->Fill(1.,track->Pt());
     if(TMath::Abs(pdgCode)==321) fProtoncont->Fill(3.,track->Pt());
     if(TMath::Abs(pdgCode)==2212) fProtoncont->Fill(5.,track->Pt());
     if(TMath::Abs(pdgCode)!=211 && TMath::Abs(pdgCode)!=321 && TMath::Abs(pdgCode)!=2212) fProtoncont->Fill(7.,track->Pt());
   }
 if(particletypeMC==SpUndefined )//these undefined are not due to absence of proper TOF response, rather due to the PID method only
   {
     if(TMath::Abs(pdgCode)==211) fUNIDcont->Fill(1.,track->Pt());
     if(TMath::Abs(pdgCode)==321) fUNIDcont->Fill(3.,track->Pt());
     if(TMath::Abs(pdgCode)==2212) fUNIDcont->Fill(5.,track->Pt());
     if(TMath::Abs(pdgCode)!=211 && TMath::Abs(pdgCode)!=321 && TMath::Abs(pdgCode)!=2212) fUNIDcont->Fill(7.,track->Pt());
   }

 if(particletypeMC==SpUndefined) continue;


    if(fRequestEventPlane){
      FillPIDEventPlane(cent_v0,particletypeMC,track->Phi(),gReactionPlane);
    }

if((track->Pt()>=fminPtAsso && track->Pt()<=fmaxPtAsso) || (track->Pt()>=fminPtTrig && track->Pt()<=fmaxPtTrig))//to reduce memory consumption in pool
  {
    Short_t chargeval=0;
    if(track->Charge()>0)   chargeval=1;
    if(track->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
if (fapplyTrigefficiency || fapplyAssoefficiency)
    effmatrix=GetTrackbyTrackeffvalue(track,effcent,zvtx,particletypeMC);//get the tracking eff x TOF matching eff x PID eff x contamination factor for identified particles 
 LRCParticlePID* copy1 = new LRCParticlePID(particletypeMC,Inv_mass,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix,track->GetTPCClusterMapPtr(),track->GetTPCSharedMapPtr());
    copy1->SetUniqueID(eventno * 100000 + (Int_t)trackscount);
    tracksID->Add(copy1);
  }
  }// if(tracktype==1) condition structure ands

}//reco track loop ends

  //*************************************************************still in event loop
 
 if(trackMap) delete trackMap;

if(trackscount>0.0)
  { 
//fill the centrality/multiplicity distribution of the selected events
 fhistcentrality->Fill(cent_v0);//*********************************WARNING::binning of cent_v0 is different for pp and pPb/PbPb case

 if (fSampleType=="pPb" || fSampleType=="PbPb" || fPPVsMult==kTRUE || fCentralityMethod == "MC_b") fCentralityCorrelation->Fill(cent_v0, trackscount);//only with unidentified tracks(i.e before PID selection);;;;;can be used to remove centrality outliers??????

//***************************************event no. counting
Bool_t isbaryontrig=kFALSE;
Bool_t ismesontrig=kFALSE;
if(tracksUNID && tracksUNID->GetEntriesFast()>0) fEventno->Fill(cent_v0,zvtx);

if(tracksID && tracksID->GetEntriesFast()>0)
  {
for(Int_t i=0;i<tracksID->GetEntriesFast();i++)
    {  //trigger loop starts
      LRCParticlePID *trig=(LRCParticlePID*)(tracksID->UncheckedAt(i));
      if(!trig) continue;
      if(trig->Pt()<fminPtTrig || trig->Pt()>fmaxPtTrig) continue;
      Int_t particlepidtrig=trig->getparticle(); //either 1 or 2
      if(particlepidtrig==SpProton) isbaryontrig=kTRUE;
      if(particlepidtrig==SpPion) ismesontrig=kTRUE;
    }//trig loop ends
 if (isbaryontrig) fEventnobaryon->Fill(cent_v0,zvtx); 
 if (ismesontrig) fEventnomeson->Fill(cent_v0,zvtx);
  }


 if(fV0TrigCorr){
 tracksIDV0=GetV0Particles((AliVEvent*) aod, cent_v0);
 if(tracksIDV0->GetEntriesFast()<=0) return;
 }
 //same event delte-eta, delta-phi plot
if(tracksUNID && tracksUNID->GetEntriesFast()>0)//hadron triggered correlation
  {//same event calculation starts
    if(ffilltrigassoUNID) Fillcorrelation(gReactionPlane,tracksUNID,0,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigassoUNID");//mixcase=kFALSE (hadron-hadron correlation)
    if(tracksID && tracksID->GetEntriesFast()>0 && ffilltrigUNIDassoID)  Fillcorrelation(gReactionPlane,tracksUNID,tracksID,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigUNIDassoID");//mixcase=kFALSE (hadron-ID correlation)
  }

if(tracksID && tracksID->GetEntriesFast()>0)//ID triggered correlation
  {//same event calculation starts
    if(tracksUNID && tracksUNID->GetEntriesFast()>0 && ffilltrigIDassoUNID) {
      if(fV0TrigCorr) Fillcorrelation(gReactionPlane,tracksIDV0,tracksUNID,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoUNID");//mixcase=kFALSE (ID-hadron correlation)
      
 else Fillcorrelation(gReactionPlane,tracksID,tracksUNID,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoUNID");//mixcase=kFALSE (ID-hadron correlation)
    }
    if(ffilltrigIDassoID)   Fillcorrelation(gReactionPlane,tracksID,0,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoID");//mixcase=kFALSE (ID-ID correlation)
  }

//still in  main event loop
//start mixing
 if(ffilltrigassoUNID || ffilltrigIDassoUNID){//mixing with unidentified particles
  AliEventPool* pool = fPoolMgr->GetEventPool(cent_v0, zvtx,gReactionPlane);//In the pool there is tracksUNID(i.e associateds are unidentified)
if (pool && pool->IsReady())
  {//start mixing only when pool->IsReady
    Float_t nmix1=(Float_t)pool->GetCurrentNEvents();  
 for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks = pool->GetEvent(jMix);
  if(!bgTracks) continue;
  if(ffilltrigassoUNID && tracksUNID && tracksUNID->GetEntriesFast()>0)//*******************************hadron trggered mixing
    Fillcorrelation(gReactionPlane,tracksUNID,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigassoUNID");//mixcase=kTRUE
 if(ffilltrigIDassoUNID && tracksID && tracksID->GetEntriesFast()>0)//***********************************ID trggered mixing
   {
     if(fV0TrigCorr) Fillcorrelation(gReactionPlane,tracksIDV0,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoUNID");//mixcase=kTRUE 
     
 else  Fillcorrelation(gReactionPlane,tracksID,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoUNID");//mixcase=kTRUE 
   }
   }// pool event loop ends mixing case

} //if pool->IsReady() condition ends mixing case
 if(tracksUNID) {
if(pool)
  pool->UpdatePool(CloneAndReduceTrackList(tracksUNID));
 }
 }//mixing with unidentified particles

 if(ffilltrigUNIDassoID || ffilltrigIDassoID){//mixing with identified particles
  AliEventPool* pool1 = fPoolMgr->GetEventPool(cent_v0, zvtx+100,gReactionPlane);//In the pool1 there is tracksID(i.e associateds are identified)
if (pool1 && pool1->IsReady())
  {//start mixing only when pool->IsReady
  Float_t nmix2=(Float_t)pool1->GetCurrentNEvents();  
for (Int_t jMix=0; jMix<pool1->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks2 = pool1->GetEvent(jMix);
  if(!bgTracks2) continue;
if(ffilltrigUNIDassoID && tracksUNID && tracksUNID->GetEntriesFast()>0)
  Fillcorrelation(gReactionPlane,tracksUNID,bgTracks2,cent_v0,zvtx,nmix2,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigUNIDassoID");//mixcase=kTRUE  
if(ffilltrigIDassoID && tracksID && tracksID->GetEntriesFast()>0)
  Fillcorrelation(gReactionPlane,tracksID,bgTracks2,cent_v0,zvtx,nmix2,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoID");//mixcase=kTRUE

   }// pool event loop ends mixing case
} //if pool1->IsReady() condition ends mixing case

if(tracksID) {
if(pool1) 
  pool1->UpdatePool(CloneAndReduceTrackList(tracksID));//ownership of tracksasso is with pool now, don't delete it(tracksUNID is with pool)
 }
 }//mixing with identified particles

  //no. of events analyzed
fEventCounter->Fill(15);
  }

//if(tracksUNID_t) delete tracksUNID_t;
if(tracksUNID)  delete tracksUNID;

if(tracksID) delete tracksID;
 
if(tracksIDV0) delete tracksIDV0;
  


}//AOD || MCAOD condition

//still in the main event loop

}
//________________________________________________________________________
void AliTwoParticlePIDCorr::doAODevent() 
{
  //get AOD
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }
  // AliAODHeader *header=(AliAODHeader*)InputEvent()->GetHeader();
  //printf("Run No: %d, ESD File: %s, Event in ESD File: %d\n",header->GetRunNumber(), header->GetESDFileName().Data(),header->GetEventNumberESDFile());
  
  
  //TString firedTriggerClasses=aod->GetFiredTriggerClasses();
  //if(firedTriggerClasses.Contains("CSH1-B")){
 
 //Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fSelectBit);
 //if(!isSelected) return;
   //cout<<"*************************************************"<<firedTriggerClasses<<endl;
  Double_t Inv_mass=0.0;

// count all events   
  fEventCounter->Fill(1);

   fPID = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPID) AliFatal("This Task needs the PID response attached to the inputHandler");
    //if (!fPID) return;//this should be available with each event even if we don't do PID selection

    fgPsi2v0a=999.;
    fgPsi2v0c=999.;
    fgPsi2tpc=999.;
    fgPsi3v0a=999.;
    fgPsi3v0c=999.;
    fgPsi3tpc=999.;
    gReactionPlane = 999.;
    
  Double_t cent_v0=   -999.;
  Double_t effcent=1.0;
  Float_t bSign = 0.;
  Double_t trackscount=0;//counts particles passed filterbit cuts and kinematic cuts used in this analysis


 bSign = (aod->GetMagneticField() > 0) ? 1 : -1;//for two track efficiency cut in correlation function calculation
 Float_t bSign1=aod->GetMagneticField() ;//for dca cut in ClassifyTrack(), i.e in track loop


// check event cuts and fill event histograms and return the centrality or reference multiplicity value
 if((cent_v0 = GetAcceptedEventMultiplicity((AliVEvent*)aod,kFALSE)) < 0){ 
    return;
  }
 effcent=cent_v0;//required for efficiency correction case********Extremely Important
  //get the event plane in case of PbPb
    if(fRequestEventPlane){
      gReactionPlane = GetEventPlane((AliVEvent*)aod,kFALSE,cent_v0);
    if(gReactionPlane==999.) return;
    }    
    

TExMap *trackMap = new TExMap();
// --- track loop for mapping matrix
 if(fFilterBit==128)
   {
 for (Int_t itrk = 0; itrk < aod->GetNumberOfTracks(); itrk++) 
{ //track loop starts for TObjArray(containing track and event information) filling; used for correlation function calculation 
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrk));
  if (!track) continue;
  Int_t tracktype=ClassifyTrack(track,trkVtx,bSign1,kFALSE);//don't fill the histos here
  if(tracktype!=1) continue; 

  if(!track) continue;//for safety

   Int_t gid = track->GetID();
   trackMap->Add(gid,itrk);
 }//track looop ends
   }

 /*
   TObjArray*  tracksUNID_t= new TObjArray;//track info before doing PID
   tracksUNID_t->SetOwner(kTRUE);  // IMPORTANT!
 */
 
   TObjArray*  tracksUNID= new TObjArray;//track info before doing PID(required for reshuffling of charges for balance function calculation)
   tracksUNID->SetOwner(kTRUE);  // IMPORTANT!
  
   TObjArray* tracksID= new TObjArray;//only pions, kaons,protons i.e. after doing the PID selection
   tracksID->SetOwner(kTRUE);  // IMPORTANT!
 
   //get the selected v0 particle TObjArray
   TObjArray*  tracksIDV0= new TObjArray;
   tracksIDV0->SetOwner(kTRUE);  // IMPORTANT!
    
    eventno++;

    Bool_t fTrigPtmin1=kFALSE;
    Bool_t fTrigPtmin2=kFALSE;
    Bool_t fTrigPtJet=kFALSE;

 for (Int_t itrk = 0; itrk < aod->GetNumberOfTracks(); itrk++) 
{ //track loop starts for TObjArray(containing track and event information) filling; used for correlation function calculation 
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(aod->GetTrack(itrk));
  if (!track) continue;
  fHistQA[11]->Fill(track->GetTPCNcls());
  Int_t particletype=-9999;//required for PID filling
  Int_t tracktype=ClassifyTrack(track,trkVtx,bSign1,kTRUE);//dcacut=kFALSE,onlyprimary=kFALSE
  if(tracktype!=1) continue; 

  if(!track) continue;//for safety

  if (fElectronRejection){
  	//Fill QA after the PID
	fHistdEdxVsPTPCafterPIDelectron -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
	fHistNSigmaTPCvsPtafterPIDelectron -> Fill(track->P()*track->Charge(),TMath::Abs(fPID->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron)));} 

  if(fV0TrigCorr){ 
if(IsTrackFromV0(aod,track)) continue;// remove auto correlation
  }
AliAODTrack *PIDtrack=track;//for PID purpose, mainly important for TPC only tracks

  if(fFilterBit==128){
Int_t gid1 = track->GetID();
//if(gid1>=0) PIDtrack = track;
 PIDtrack =(AliAODTrack*) aod->GetTrack(trackMap->GetValue(-1-gid1));
if(!PIDtrack) continue;//for safety; so that each of the TPC only tracks have corresponding global track along with it
  }

//check for eta , phi holes
 fEtaSpectrasso->Fill(track->Eta(),track->Pt());
 fphiSpectraasso->Fill(track->Phi(),track->Pt());

 trackscount++;
 
 //if no applyefficiency , set the eff factor=1.0
 Float_t effmatrix=1.0;

 //tag all particles as unidentified that passed filterbit & kinematic cuts 
 particletype=unidentified;

 //To count the no. of tracks having an accepted track in a certain PT(e.g. Jet Pt) range
 if(track->Pt()>=fminPtTrig) fTrigPtmin1=kTRUE;
 if(track->Pt()>=(fminPtTrig+0.5)) fTrigPtmin2=kTRUE;
 if(track->Pt()>=fmaxPtTrig) fTrigPtJet=kTRUE;


if (fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMult==kFALSE)) effcent=15.0;//integrated over multiplicity [i.e each track has multiplicity 15.0](so put any fixed value for each track so that practically means there is only one bin in multiplicityi.e multiplicity intregated out )**************Important for efficiency related issues


 //to reduce memory consumption in pool
  if((track->Pt()>=fminPtAsso && track->Pt()<=fmaxPtAsso) || (track->Pt()>=fminPtTrig && track->Pt()<=fmaxPtTrig)) 
  {
 //Clone & Reduce track list(TObjArray) for unidentified particles
    Short_t chargeval=0;
    if(track->Charge()>0)   chargeval=1;
    if(track->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
 if (fapplyTrigefficiency || fapplyAssoefficiency)//get the trackingefficiency x contamination factor for unidentified particles
   effmatrix=GetTrackbyTrackeffvalue(track,effcent,zvtx,particletype);
 LRCParticlePID* copy = new LRCParticlePID(particletype,Inv_mass,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix,track->GetTPCClusterMapPtr(),track->GetTPCSharedMapPtr());
  copy->SetUniqueID(eventno * 100000 + (Int_t)trackscount);
  tracksUNID->Add(copy);//track information Storage for UNID correlation function(tracks that pass the filterbit & kinematic cuts only)
  }
  
  /*
   if(fRunShufflingbalance){
    tracksUNID = GetShuffledTracks(tracksUNID_t);
  }
   else tracksUNID=CloneAndReduceTrackList(tracksUNID_t);
  */
//now start the particle identificaion process:)

  // DCA XY cut to check contaminations from lambda using dca sigma cut variation with filterbit 16, only for identified particles 
     if (fDCAXYCut)
	{
	  if (!trkVtx) continue;
	  
	  Double_t pos[2];
	  Double_t covar[3];
	  AliAODTrack* clone =(AliAODTrack*) track->Clone();
	  Bool_t success = clone->PropagateToDCA(trkVtx, bSign1, fdcacutvalue, pos, covar);
	  delete clone;
	  if (!success)
	    continue;

// 	  Printf("%f", ((AliAODTrack*)part)->DCA());
// 	  Printf("%f", pos[0]);
	  if (TMath::Abs(pos[0]) > fDCAXYCut->Eval(track->Pt())) continue;
	}

//track passing filterbit 768 have proper TPC response,or need to be checked explicitly before doing PID????

  Float_t dEdx = PIDtrack->GetTPCsignal();
  fHistoTPCdEdx->Fill(track->Pt(), dEdx);

  //fill beta vs Pt plots only for tracks having proper TOF response(much less tracks compared to the no. that pass the filterbit & kinematic cuts)
 if(HasTOFPID(PIDtrack))
{
  Double_t beta = GetBeta(PIDtrack);
  fHistoTOFbeta->Fill(track->Pt(), beta);
 }
  
 //remove the tracks which don't have proper TOF response-otherwise the misIDentification rate values will be wrong(in MC)
if(fRequestTOFPID && track->Pt()>fPtTOFPIDmin && track->Pt()<fPtTOFPIDmax && (!HasTOFPID(PIDtrack)) ) continue;

//above fPtTOFPIDmax PID is done using rel. rise in TPC only
 if(track->Pt()>=fPtTOFPIDmax){
   if (!HasTPCPID(track)) continue;
   Bool_t TPCsectoredge=TPCCutMIGeo(track,aod);
   if(TPCSectoredgecut) {if(!TPCsectoredge) continue;}
 }
 
//track identification(using nsigma method)
     particletype=GetParticle(PIDtrack,fFIllPIDQAHistos);//*******************************change may be required(It should return only pion,kaon, proton and Spundefined; NOT unidentifed***************be careful)

//2-d TPCTOF map(for each Pt interval)
  if(HasTOFPID(PIDtrack)){
 fTPCTOFPion3d->Fill(track->Pt(),fnsigmas[SpPion][NSigmaTOF],fnsigmas[SpPion][NSigmaTPC]);
 fTPCTOFKaon3d->Fill(track->Pt(),fnsigmas[SpKaon][NSigmaTOF],fnsigmas[SpKaon][NSigmaTPC]);
 fTPCTOFProton3d->Fill(track->Pt(),fnsigmas[SpProton][NSigmaTOF],fnsigmas[SpProton][NSigmaTPC]); 
  }

//ignore the Spundefined particles as they also contain pion, kaon, proton outside the nsigma cut(also if tracks don't have proper TOF PID in a certain Pt interval) & these tracks are actually counted when we do the the efficiency correction, so considering them as unidentified particles & doing the efficiency correction(i.e defining unidentified=pion+Kaon+proton+SpUndefined is right only without efficiency correction) for them will be two times wrong!!!!! 
  if (particletype==SpUndefined) continue;//this condition creating a modulated structure in delphi projection in mixed event case(only when we are dealing with identified particles i.e. tracksID)!!!!!!!!!!!


    if(fRequestEventPlane){
      FillPIDEventPlane(cent_v0,particletype,track->Phi(),gReactionPlane);
    }
    
 //Pt, Eta , Phi distribution of the reconstructed identified particles
if(ffillhistQAReco)
    {
if (particletype==SpPion)
  {
    fPionPt->Fill(track->Pt());
    fPionEta->Fill(track->Eta());
    fPionPhi->Fill(track->Phi());
  }
if (particletype==SpKaon)
  {
    fKaonPt->Fill(track->Pt());
    fKaonEta->Fill(track->Eta());
    fKaonPhi->Fill(track->Phi());
  }
if (particletype==SpProton)
  {
    fProtonPt->Fill(track->Pt());
    fProtonEta->Fill(track->Eta());
    fProtonPhi->Fill(track->Phi());
  }
    }
 
 if((track->Pt()>=fminPtAsso && track->Pt()<=fmaxPtAsso) || (track->Pt()>=fminPtTrig && track->Pt()<=fmaxPtTrig)) 
  {
    Short_t chargeval=0;
    if(track->Charge()>0)   chargeval=1;
    if(track->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
if (fapplyTrigefficiency || fapplyAssoefficiency)
  effmatrix=GetTrackbyTrackeffvalue(track,effcent,zvtx,particletype);//get the tracking eff x TOF matching eff x PID eff x contamination factor for identified particles; Bool_t mesoneffrequired=kFALSE
 LRCParticlePID* copy1 = new LRCParticlePID(particletype,Inv_mass,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix,track->GetTPCClusterMapPtr(),track->GetTPCSharedMapPtr());
    copy1->SetUniqueID(eventno * 100000 + (Int_t)trackscount);
    tracksID->Add(copy1);
  }
} //track loop ends but still in event loop

 if(trackMap) delete trackMap;

 
if(trackscount<1.0){
  if(tracksUNID) delete tracksUNID;
  if(tracksID) delete tracksID;
  return;
 }

 if (fTrigPtmin1) fhistJetTrigestimate->Fill(cent_v0,0.0);
 if (fTrigPtmin2) fhistJetTrigestimate->Fill(cent_v0,2.0);
 if (fTrigPtJet) fhistJetTrigestimate->Fill(cent_v0,4.0);

 Float_t weightval=1.0;

  
//fill the centrality/multiplicity distribution of the selected events
 fhistcentrality->Fill(cent_v0);//*********************************WARNING::binning of cent_v0 is different for pp and pPb/PbPb case

if(fSampleType=="pPb" || fSampleType=="PbPb" || fPPVsMult==kTRUE) fCentralityCorrelation->Fill(cent_v0, trackscount);//only with unidentified tracks(i.e before PID selection);;;;;can be used to remove centrality outliers??????

//count selected events having centrality betn 0-100%
 fEventCounter->Fill(13);

//***************************************event no. counting
Bool_t isbaryontrig=kFALSE;
Bool_t ismesontrig=kFALSE;
if(tracksUNID && tracksUNID->GetEntriesFast()>0) fEventno->Fill(cent_v0,zvtx);

if(tracksID && tracksID->GetEntriesFast()>0)
  {
for(Int_t i=0;i<tracksID->GetEntriesFast();i++)
    {  //trigger loop starts
      LRCParticlePID *trig=(LRCParticlePID*)(tracksID->UncheckedAt(i));
      if(!trig) continue;
      if(trig->Pt()<fminPtTrig || trig->Pt()>fmaxPtTrig) continue;
      Int_t particlepidtrig=trig->getparticle(); //either 1 or 2
      if(particlepidtrig==SpProton) isbaryontrig=kTRUE;
      if(particlepidtrig==SpPion) ismesontrig=kTRUE;
    }//trig loop ends
 if (isbaryontrig) fEventnobaryon->Fill(cent_v0,zvtx); 
 if (ismesontrig) fEventnomeson->Fill(cent_v0,zvtx);
  }

//Get the TObjArray of V0 particles
 if(fV0TrigCorr){
 tracksIDV0=GetV0Particles((AliVEvent*) aod,cent_v0);
 if(tracksIDV0->GetEntriesFast()<=0) return;
 }

//same event delta-eta-deltaphi plot 
if(tracksUNID && tracksUNID->GetEntriesFast()>0)//hadron triggered correlation
  {//same event calculation starts
    if(ffilltrigassoUNID) Fillcorrelation(gReactionPlane,tracksUNID,0,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigassoUNID");//mixcase=kFALSE (hadron-hadron correlation)
    if(tracksID && tracksID->GetEntriesFast()>0 && ffilltrigUNIDassoID) Fillcorrelation(gReactionPlane,tracksUNID,tracksID,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigUNIDassoID");//mixcase=kFALSE (hadron-ID correlation)
    
  }

if(tracksID && tracksID->GetEntriesFast()>0)//ID triggered correlation
  {//same event calculation starts
    if(tracksUNID && tracksUNID->GetEntriesFast()>0 && ffilltrigIDassoUNID){
if(fV0TrigCorr)  Fillcorrelation(gReactionPlane,tracksIDV0,tracksUNID,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoUNID");//mixcase=kFALSE (ID-hadron correlation)
else  Fillcorrelation(gReactionPlane,tracksID,tracksUNID,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoUNID");//mixcase=kFALSE (ID-hadron correlation)
    }
    if(ffilltrigIDassoID)   Fillcorrelation(gReactionPlane,tracksID,0,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoID");//mixcase=kFALSE (ID-ID correlation)
  }

//still in  main event loop


//start mixing
 if(ffilltrigassoUNID || ffilltrigIDassoUNID){//mixing with unidentified particles
AliEventPool* pool = fPoolMgr->GetEventPool(cent_v0, zvtx,gReactionPlane);//In the pool there is tracksUNID(i.e associateds are unidentified)
if (pool && pool->IsReady())
  {//start mixing only when pool->IsReady
  Float_t nmix1=(Float_t)pool->GetCurrentNEvents();  
 for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks = pool->GetEvent(jMix);
  if(!bgTracks) continue;
  if(ffilltrigassoUNID && tracksUNID && tracksUNID->GetEntriesFast()>0)//*******************************hadron trggered mixing
    Fillcorrelation(gReactionPlane,tracksUNID,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigassoUNID");//mixcase=kTRUE
 if(ffilltrigIDassoUNID && tracksID && tracksID->GetEntriesFast()>0)//***********************************ID trggered mixing
   {
if(fV0TrigCorr)   Fillcorrelation(gReactionPlane,tracksIDV0,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoUNID");//mixcase=kTRUE
else  Fillcorrelation(gReactionPlane,tracksID,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoUNID");//mixcase=kTRUE
   } 
   }// pool event loop ends mixing case
} //if pool->IsReady() condition ends mixing case
 if(tracksUNID) {
if(pool)
  pool->UpdatePool(CloneAndReduceTrackList(tracksUNID));
 }
 }//mixing with unidentified particles


 if(ffilltrigUNIDassoID || ffilltrigIDassoID){//mixing with identified particles
 AliEventPool* pool1 = fPoolMgr->GetEventPool(cent_v0, zvtx+100,gReactionPlane);//In the pool1 there is tracksID(i.e associateds are identified)
if (pool1 && pool1->IsReady())
  {//start mixing only when pool->IsReady
  Float_t nmix2=(Float_t)pool1->GetCurrentNEvents();  
for (Int_t jMix=0; jMix<pool1->GetCurrentNEvents(); jMix++) 
  { //pool event loop start
 TObjArray* bgTracks2 = pool1->GetEvent(jMix);
  if(!bgTracks2) continue;
if(ffilltrigUNIDassoID && tracksUNID && tracksUNID->GetEntriesFast()>0)
  Fillcorrelation(gReactionPlane,tracksUNID,bgTracks2,cent_v0,zvtx,nmix2,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigUNIDassoID");//mixcase=kTRUE  
if(ffilltrigIDassoID && tracksID && tracksID->GetEntriesFast()>0)
  Fillcorrelation(gReactionPlane,tracksID,bgTracks2,cent_v0,zvtx,nmix2,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoID");//mixcase=kTRUE

   }// pool event loop ends mixing case
} //if pool1->IsReady() condition ends mixing case

if(tracksID) {
if(pool1) 
  pool1->UpdatePool(CloneAndReduceTrackList(tracksID));//ownership of tracksasso is with pool now, don't delete it(tracksUNID is with pool)
 }
 }//mixing with identified particles


  //no. of events analyzed
fEventCounter->Fill(15);
 
//still in main event loop
//if(tracksUNID_t) delete tracksUNID_t;

if(tracksUNID)  delete tracksUNID;

if(tracksID) delete tracksID;

if(tracksIDV0) delete tracksIDV0;
  

} // *************************event loop ends******************************************
//________________________________________________________________________
TObjArray* AliTwoParticlePIDCorr::GetShuffledTracks(TObjArray *tracks){//taken from TaskBFpsi
  // Clones TObjArray and returns it with tracks after shuffling the charges

  TObjArray* tracksShuffled = new TObjArray;
  tracksShuffled->SetOwner(kTRUE);

  vector<Short_t> *chargeVector = new vector<Short_t>;   //original charge of accepted tracks 

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    LRCParticlePID* track = (LRCParticlePID*) tracks->UncheckedAt(i);
    chargeVector->push_back(track->Charge());
  }  
 
  std::random_device rd;
  std::default_random_engine engine{rd()};
  std::shuffle(chargeVector->begin(), chargeVector->end(),engine);
  
  for(Int_t i = 0; i < tracks->GetEntriesFast(); i++){
    LRCParticlePID* track = (LRCParticlePID*) tracks->UncheckedAt(i);

    //==============================correction(At the moment we are assuming that the correction factoer is harge independent)
    //Double_t correction = GetTrackbyTrackCorrectionMatrix(track->Eta(), track->Phi(),track->Pt(), chargeVector->at(i), gCentrality);
    //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);
    //tracksShuffled->Add(new AliBFBasicParticle(track->Eta(), track->Phi(), track->Pt(),chargeVector->at(i), correction));


     LRCParticlePID* copy100 = new LRCParticlePID(track->getparticle(),track->GetInvMass(),chargeVector->at(i), track->Pt(),track->Eta(), track->Phi(), track->geteffcorrectionval(),track->GetTPCPadMap(),track->GetTPCSharedMap());
     copy100->SetUniqueID(track->GetUniqueID());
     tracksShuffled->Add(copy100);

  }

  delete chargeVector;
   
  return tracksShuffled;
}

//________________________________________________________________________

TObjArray* AliTwoParticlePIDCorr::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliDPhiBasicParticle which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    LRCParticlePID* particle = (LRCParticlePID*) tracks->UncheckedAt(i);
    LRCParticlePID* copy100 = new LRCParticlePID(particle->getparticle(),particle->GetInvMass(),particle->Charge(), particle->Pt(),particle->Eta(), particle->Phi(), particle->geteffcorrectionval(),particle->GetTPCPadMap(),particle->GetTPCSharedMap());
    copy100->SetUniqueID(particle->GetUniqueID());
    tracksClone->Add(copy100);
  }
  
  return tracksClone;
}

//--------------------------------------------------------------------------------
void AliTwoParticlePIDCorr::Fillcorrelation(Float_t ReactionPlane,TObjArray *trackstrig,TObjArray *tracksasso,Double_t cent,Float_t vtx,Float_t weight,Bool_t firstTime,Float_t bSign,Bool_t fPtOrder,Bool_t twoTrackEfficiencyCut,Bool_t mixcase,TString fillup)
{

  //before calling this function check that either trackstrig & tracksasso are available 

 // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* input = (tracksasso) ? tracksasso : trackstrig;
  TArrayF eta(input->GetEntriesFast());
  for (Int_t i=0; i<input->GetEntriesFast(); i++)
    eta[i] = ((LRCParticlePID*) input->UncheckedAt(i))->Eta();

  //if(trackstrig)
    Int_t jmax=trackstrig->GetEntriesFast();
  if(tracksasso)
     jmax=tracksasso->GetEntriesFast();

// identify K, Lambda candidates and flag those particles
    // a TObject bit is used for this
const UInt_t kResonanceDaughterFlag = 1 << 14;
    if (fRejectResonanceDaughters > 0)
    {
      Double_t resonanceMass = -1;
      Double_t massDaughter1 = -1;
      Double_t massDaughter2 = -1;
      const Double_t interval = 0.02;
 switch (fRejectResonanceDaughters)
      {
	case 1: resonanceMass = 0.9; massDaughter1 = 0.1396; massDaughter2 = 0.9383; break; // method test
	case 2: resonanceMass = 0.4976; massDaughter1 = 0.1396; massDaughter2 = massDaughter1; break; // k0
	case 3: resonanceMass = 1.115; massDaughter1 = 0.1396; massDaughter2 = 0.9383; break; // lambda
	default: AliFatal(Form("Invalid setting %d", fRejectResonanceDaughters));
      }      

for (Int_t i=0; i<trackstrig->GetEntriesFast(); i++)
	trackstrig->UncheckedAt(i)->ResetBit(kResonanceDaughterFlag);
 for (Int_t i=0; tracksasso->GetEntriesFast(); i++)
	  tracksasso->UncheckedAt(i)->ResetBit(kResonanceDaughterFlag);

 for (Int_t i=0; i<trackstrig->GetEntriesFast(); i++)
      {
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
for (Int_t j=0; tracksasso->GetEntriesFast(); j++)
	{
        LRCParticlePID *asso=(LRCParticlePID*)(tracksasso->UncheckedAt(j));

 // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event)
if (trig->IsEqual(asso)) continue;

if (trig->Charge() * asso->Charge() > 0) continue;

 Float_t mass = GetInvMassSquaredCheap(trig->Pt(), trig->Eta(), trig->Phi(), asso->Pt(), asso->Eta(), asso->Phi(), massDaughter1, massDaughter2);
     
if (TMath::Abs(mass - resonanceMass*resonanceMass) < interval*5)
	  {
	    mass = GetInvMassSquared(trig->Pt(), trig->Eta(), trig->Phi(), asso->Pt(), asso->Eta(), asso->Phi(), massDaughter1, massDaughter2);

	    if (mass > (resonanceMass-interval)*(resonanceMass-interval) && mass < (resonanceMass+interval)*(resonanceMass+interval))
	    {
	      trig->SetBit(kResonanceDaughterFlag);
	      asso->SetBit(kResonanceDaughterFlag);
	      
// 	      Printf("Flagged %d %d %f", i, j, TMath::Sqrt(mass));
	    }
	  }
	}
      }
    }

      //Select the highest Pt trigger particle in an event (within a given Pt trigger range)

    Float_t TriggerPtMin=fminPtTrig;
    Float_t TriggerPtMax=fmaxPtTrig;
    Int_t HighestPtTriggerIndx=-99999;
    TH1* triggerWeighting = 0;

if(fSelectHighestPtTrig || fWeightPerEvent)//**************add this data member to the constructor
      {
if (fWeightPerEvent)
    {
      TAxis* axis=0;
   if(ffilltrigassoUNID || ffilltrigUNIDassoID || ffilltrigIDassoUNID || ffilltrigIDassoID) axis = fTHnTrigcount->GetGrid(0)->GetGrid()->GetAxis(2);                                          
  if((fAnalysisType =="MCAOD") && ffilltrigIDassoIDMCTRUTH)    axis = fTHnTrigcountMCTruthPrim->GetGrid(0)->GetGrid()->GetAxis(2);
      triggerWeighting = new TH1F("triggerWeighting", "", axis->GetNbins(), axis->GetXbins()->GetArray());
    }
for(Int_t i=0;i<trackstrig->GetEntriesFast();i++)
    {  //trigger loop starts(highest Pt trigger selection)
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
      if(!trig) continue;
      Float_t trigpt=trig->Pt();
    //to avoid overflow qnd underflow
      if(trigpt<fminPtTrig || trigpt>fmaxPtTrig) continue;
      Int_t particlepidtrig=trig->getparticle();
      if(particlepidtrig==999) continue;
      if(fTriggerSpeciesSelection){ if (particlepidtrig!=fTriggerSpecies) continue;}

      Float_t trigeta=trig->Eta();

      // some optimization
 if (fTriggerRestrictEta > 0 && TMath::Abs(trigeta) > fTriggerRestrictEta)
	continue;

if (fOnlyOneEtaSide != 0)
      {
	if (fOnlyOneEtaSide * trigeta < 0)
	  continue;
      }
  if (fTriggerSelectCharge != 0)
	if (trig->Charge() * fTriggerSelectCharge < 0)
	  continue;
	
      if (fRejectResonanceDaughters > 0)
	if (trig->TestBit(kResonanceDaughterFlag)) continue;

      if(fSelectHighestPtTrig){
 if(trigpt>TriggerPtMin && trigpt<=TriggerPtMax)
          {         
	  HighestPtTriggerIndx=(Int_t)trig->GetUniqueID();
          TriggerPtMin=trigpt;
          }
      }

if (fWeightPerEvent)  triggerWeighting->Fill(trigpt);

    }//trigger loop ends(highest Pt trigger selection)

      }//******************(fSelectHighestPtTrig || fWeightPerEvent) condition ends


 //two particle correlation filling
for(Int_t i=0;i<trackstrig->GetEntriesFast();i++)
    {  //trigger loop starts
      LRCParticlePID *trig=(LRCParticlePID*)(trackstrig->UncheckedAt(i));
      if(!trig) continue;
      Float_t trigpt=trig->Pt();
    //to avoid overflow qnd underflow
      if(trigpt<fminPtTrig || trigpt>fmaxPtTrig) continue;
      /*
      Double_t ParticlePID_InvMass=0.0;
      if(fV0TrigCorr) ParticlePID_InvMass=trig->GetInvMass();
      else{
      Int_t particlepidtrig=trig->getparticle();
      ParticlePID_InvMass=(Double_t) particlepidtrig;
      if(fTriggerSpeciesSelection){ if (particlepidtrig!=fTriggerSpecies) continue;}
      }
*/

      Int_t particlepidtrig=trig->getparticle();
      if(particlepidtrig==999) continue;
       if(fTriggerSpeciesSelection){ if (particlepidtrig!=fTriggerSpecies) continue;}//***********************************forks,lam.Alam their PID numbers have no meaning, only their Inv_mass will be stored
      
      Double_t ParticlePID_InvMass=(Double_t) particlepidtrig;
     
      Float_t trigeta=trig->Eta();

      // some optimization
 if (fTriggerRestrictEta > 0 && TMath::Abs(trigeta) > fTriggerRestrictEta)
	continue;

if (fOnlyOneEtaSide != 0)
      {
	if (fOnlyOneEtaSide * trigeta < 0)
	  continue;
      }
  if (fTriggerSelectCharge != 0)
	if (trig->Charge() * fTriggerSelectCharge < 0)
	  continue;
	
      if (fRejectResonanceDaughters > 0)
	if (trig->TestBit(kResonanceDaughterFlag)) continue;

      if(fSelectHighestPtTrig && HighestPtTriggerIndx!=-99999) {
	if(trig->GetUniqueID()!=(UInt_t)HighestPtTriggerIndx) continue;
      }

      Float_t trigphi=trig->Phi();
      Float_t trackefftrig=1.0;
      if(fapplyTrigefficiency) trackefftrig=trig->geteffcorrectionval();

    // Event plane (determine psi bin)
    Double_t gPsiMinusPhi    =   0.;
    Double_t gPsiMinusPhiBin = -10.;
if(fRequestEventPlane){
    gPsiMinusPhi   = TMath::Abs(trigphi - ReactionPlane);
    //in-plane(Note thet event plane angle has to be defined within 0 to 180 degree(do not use this if event ), otherwise the definition of in plane and out plane particles is wrong)
    if((gPsiMinusPhi <= 7.5*TMath::DegToRad())||
      (gPsiMinusPhi >= 352.5*TMath::DegToRad())||
       ((172.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 0.0;
    /*
 if((gPsiMinusPhi <= 7.5*TMath::DegToRad())||
       ((172.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 0.0;
    */
    //intermediate
    else if(((37.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5*TMath::DegToRad()))||
	    ((127.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5*TMath::DegToRad()))||
	    ((217.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5*TMath::DegToRad()))||
	    ((307.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 1.0;
    //out of plane
    else if(((82.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5*TMath::DegToRad()))||
	    ((262.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 2.0;
    //everything else
    else 
      gPsiMinusPhiBin = 3.0;
    
    fHistPsiMinusPhi->Fill(gPsiMinusPhiBin,gPsiMinusPhi);
 }

      //cout<<"*******************trackefftrig="<<trackefftrig<<endl;
	Double_t* trigval;
	Int_t dim=3;
	Int_t eventplaneAxis=0;
        if(fRequestEventPlane) eventplaneAxis=1;
	if(fcontainPIDtrig && SetChargeAxis==0) dim=4+eventplaneAxis;
	if(!fcontainPIDtrig && SetChargeAxis==2) dim=4+eventplaneAxis;
	if(fcontainPIDtrig && SetChargeAxis==2) dim=5+eventplaneAxis;
        trigval= new Double_t[dim];
      trigval[0] = cent;
      trigval[1] = vtx;
      trigval[2] = trigpt;

      if(fRequestEventPlane){
      trigval[3] = gPsiMinusPhiBin;
      if(fcontainPIDtrig && SetChargeAxis==0) trigval[4] = ParticlePID_InvMass;
      if(!fcontainPIDtrig && SetChargeAxis==2) trigval[4] = trig->Charge();
      if(fcontainPIDtrig && SetChargeAxis==2) {
      trigval[4] = ParticlePID_InvMass;
      trigval[5] = trig->Charge();
       }
      }

  if(!fRequestEventPlane){
      if(fcontainPIDtrig && SetChargeAxis==0) trigval[3] = ParticlePID_InvMass;
      if(!fcontainPIDtrig && SetChargeAxis==2) trigval[3] = trig->Charge();
      if(fcontainPIDtrig && SetChargeAxis==2) {
      trigval[3] = ParticlePID_InvMass;
      trigval[4] = trig->Charge();
       }
      }

 

	if (fWeightPerEvent)
	{
	  // leads effectively to a filling of one entry per filled trigger particle pT bin
	  Int_t weightBin = triggerWeighting->GetXaxis()->FindBin(trigval[2]);
// 	  Printf("Using weight %f", triggerWeighting->GetBinContent(weightBin));
	  trackefftrig *= triggerWeighting->GetBinContent(weightBin);
	}


      //trigger particle counting for both same and mixed event case;;;;;step=0->same event case;;;;;step=1->mixed event case;;;;;;
if(ffilltrigassoUNID==kTRUE && ffilltrigUNIDassoID==kTRUE){
      if(fillup=="trigassoUNID" ) {
if(mixcase==kFALSE)   fTHnTrigcount->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcount->Fill(trigval,1,1.0/trackefftrig); 
      }
    }
 if(ffilltrigassoUNID==kTRUE && ffilltrigUNIDassoID==kFALSE){
   if(fillup=="trigassoUNID" )  
     {
if(mixcase==kFALSE)   fTHnTrigcount->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcount->Fill(trigval,1,1.0/trackefftrig); 
     }
    }
if(ffilltrigassoUNID==kFALSE && ffilltrigUNIDassoID==kTRUE){
  if(fillup=="trigUNIDassoID")  
    {
if(mixcase==kFALSE)   fTHnTrigcount->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcount->Fill(trigval,1,1.0/trackefftrig); 
    }
    }
 //ensure that trigIDassoID , trigassoUNID, trigIDassoUNID & trigUNIDassoID  case FillCorrelation called only once in the event loop for same event correlation function calculation, otherwise there will be multiple counting of pion, kaon,proton,unidentified
if(ffilltrigIDassoUNID==kTRUE && ffilltrigIDassoID==kTRUE){
  if(fillup=="trigIDassoID")  
    {
if(mixcase==kFALSE)   fTHnTrigcount->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcount->Fill(trigval,1,1.0/trackefftrig); 
    }
    }
 if(ffilltrigIDassoUNID==kTRUE && ffilltrigIDassoID==kFALSE){
   if(fillup=="trigIDassoUNID" ) 
     {
if(mixcase==kFALSE)   fTHnTrigcount->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcount->Fill(trigval,1,1.0/trackefftrig); 
     } 
    }
if(ffilltrigIDassoUNID==kFALSE && ffilltrigIDassoID==kTRUE){
  if(fillup=="trigIDassoID")  
    {
if(mixcase==kFALSE)   fTHnTrigcount->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcount->Fill(trigval,1,1.0/trackefftrig); 
    }
    }

 if(fillup=="trigIDassoIDMCTRUTH") { //In truth MC case "Unidentified" means any particle other than pion,kaon or proton and no efficiency correction(default value 1.0)************************be careful!!!! 
if(mixcase==kFALSE)   fTHnTrigcountMCTruthPrim->Fill(trigval,0,1.0/trackefftrig); 
if(mixcase==kTRUE && firstTime)   fTHnTrigcountMCTruthPrim->Fill(trigval,1,1.0/trackefftrig); 
  }

    //asso loop starts within trigger loop
   for(Int_t j=0;j<jmax;j++)
             {
    LRCParticlePID *asso=0;
    if(!tracksasso)
    asso=(LRCParticlePID*)(trackstrig->UncheckedAt(j));
    else
    asso=(LRCParticlePID*)(tracksasso->UncheckedAt(j));

    if(!asso) continue;

    //to avoid overflow and underflow
 if(asso->Pt()<fminPtAsso || asso->Pt()>fmaxPtAsso) continue;//***********************Important

 //if(fmaxPtAsso==fminPtTrig) {if(asso->Pt()==fminPtTrig) continue;}//******************Think about it!

  if(!tracksasso && j==i) continue;

   // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event,i.e. both Trig and asso TObjArray belongs to the same Pt range but say Trig is Unidentified but asso is identified then the serial no. wise particles are not same and and j==i doesn't aplly)
   // if (tracksasso && trig->IsEqual(asso))  continue;

  if (tracksasso && (trig->GetUniqueID()==asso->GetUniqueID())) continue;

  //Pt ordering acts as a nested loop which avoids double counting and a pair counts only once in case of unidentified particles
 if (fPtOrder)
 if (asso->Pt() >= trig->Pt()) continue;

  Int_t particlepidasso=asso->getparticle(); 
  if(fAssociatedSpeciesSelection){ if (particlepidasso!=fAssociatedSpecies) continue;}
	    

if (fAssociatedSelectCharge != 0)
if (asso->Charge() * fAssociatedSelectCharge < 0) continue;
	    
 if (fSelectCharge > 0)
        {
          // skip like sign
          if (fSelectCharge == 1 && asso->Charge() * trig->Charge() > 0)
            continue;
            
          // skip unlike sign
          if (fSelectCharge == 2 && asso->Charge() * trig->Charge() < 0)
            continue;
        }

if (fEtaOrdering)
	{
	  if (trigeta < 0 && asso->Eta() < trigeta)
	    continue;
	  if (trigeta > 0 && asso->Eta() > trigeta)
	    continue;
	}

if (fRejectResonanceDaughters > 0)
	  if (asso->TestBit(kResonanceDaughterFlag))
	  {
// 	    Printf("Skipped j=%d", j);
	    continue;
	  }

	// conversions
	if (fCutConversions && asso->Charge() * trig->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j], asso->Phi(), 0.510e-3, 0.510e-3);
	  
	  if (mass < 0.1)
	  {
	    mass = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(), eta[j], asso->Phi(), 0.510e-3, 0.510e-3);
	    
	    fControlConvResoncances->Fill(0.0, mass);

	    if (mass < 0.04*0.04) 
	      continue;
	  }
	}

	// K0s
	if (fCutResonances && asso->Charge() * trig->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(), eta[j], asso->Phi(), 0.1396, 0.1396);
	  
	  const Float_t kK0smass = 0.4976;
	  
	  if (TMath::Abs(mass - kK0smass*kK0smass) < 0.1)
	  {
	    mass = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j], asso->Phi(), 0.1396, 0.1396);
	    
	    fControlConvResoncances->Fill(1, mass - kK0smass*kK0smass);

	    if (mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))
	      continue;
	  }
	}
	
	// Lambda
	if (fCutResonances && asso->Charge() * trig->Charge() < 0)
	{
	  Float_t mass1 = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(), eta[j], asso->Phi(), 0.1396, 0.9383);
	  Float_t mass2 = GetInvMassSquaredCheap(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j] , asso->Phi(), 0.9383, 0.1396);
	  
	  const Float_t kLambdaMass = 1.115;

	  if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < 0.1)
	  {
	    mass1 = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j], asso->Phi(), 0.1396, 0.9383);

	    fControlConvResoncances->Fill(2, mass1 - kLambdaMass*kLambdaMass);
	    
	    if (mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
	      continue;
	  }
	  if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < 0.1)
	  {
	    mass2 = GetInvMassSquared(trig->Pt(), trigeta, trig->Phi(), asso->Pt(),eta[j] , asso->Phi(), 0.9383, 0.1396);

	    fControlConvResoncances->Fill(2, mass2 - kLambdaMass*kLambdaMass);

	    if (mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
	      continue;
	  }
	}

	if (twoTrackEfficiencyCut)
	{
	  // the variables & cuthave been developed by the HBT group 
	  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
	  Float_t phi1 = trig->Phi();
	  Float_t pt1 = trig->Pt();
	  Float_t charge1 = trig->Charge();
	  Float_t phi2 = asso->Phi();
	  Float_t pt2 = asso->Pt();
	  Float_t charge2 = asso->Charge();

	  Float_t deta= trigeta - eta[j]; 
    
 // optimization
	  if (TMath::Abs(deta) < twoTrackEfficiencyCutValue * 2.5 * 3)
	  {

  // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, fTwoTrackCutMinRadius, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, fTwoTrackCutMaxRadius, bSign);

 const Float_t kLimit = twoTrackEfficiencyCutValue * 3;

	    Float_t dphistarminabs = 1e5;
	    Float_t dphistarmin = 1e5;

 if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=fTwoTrackCutMinRadius; rad<=fTwoTrackCutMaxRadius; rad+=0.01) 
	      {
		Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);

		Float_t dphistarabs = TMath::Abs(dphistar);

	if (dphistarabs < dphistarminabs)
		{
		  dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      if(mixcase==kFALSE)  fTwoTrackDistancePt[0]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));//for same event
	      if(mixcase==kTRUE)  fTwoTrackDistancePtmix[0]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));//for mixed event

if (dphistarminabs < twoTrackEfficiencyCutValue && TMath::Abs(deta) < twoTrackEfficiencyCutValue)
	      {
// 		Printf("Removed track pair %d %d with %f %f %f %f %f %f %f %f %f", i, j, deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
		continue;
	      }
             if(mixcase==kFALSE) fTwoTrackDistancePt[1]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));//for same event
             if(mixcase==kTRUE) fTwoTrackDistancePtmix[1]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));//for mixed event

	    }
	  }
	}

	//pair sharedfraction cut(only between the trig and asso track)
   if(fillup!="trigIDassoIDMCTRUTH")//******************************************NOT for TRUTH MC particles
  {
    if(fSharedfraction_Pair_cut>=0){
	Bool_t passsharedfractionpaircut=CalculateSharedFraction(trig->GetTPCPadMap(),asso->GetTPCPadMap(),trig->GetTPCSharedMap(),asso->GetTPCSharedMap());
	if(!passsharedfractionpaircut) continue;
    }
  }
	Float_t weightperevent=weight;
        Float_t trackeffasso=1.0;
        if(fapplyAssoefficiency) trackeffasso=asso->geteffcorrectionval();
	//cout<<"*******************trackeffasso="<<trackeffasso<<endl;
        Float_t deleta=trigeta-eta[j];
	Float_t delphi=PhiRange(trigphi-asso->Phi()); 

	Float_t delpt=trigpt-asso->Pt();
	//fill it with/without two track efficiency cut	   
	if(mixcase==kFALSE)  fTwoTrackDistancePtdip->Fill(deleta, delphi, TMath::Abs(delpt));//for same event
	if(mixcase==kTRUE)  fTwoTrackDistancePtdipmix->Fill(deleta, delphi, TMath::Abs(delpt));//for mixed event

 //here get the two particle efficiency correction factor
	Float_t effweight=trackefftrig*trackeffasso*weightperevent;
	// if(mixcase==kFALSE)	cout<<"*******************effweight="<<effweight<<endl;
	Double_t* vars;
	Int_t dimused=kTrackVariablesPair+eventplaneAxis;
        vars= new Double_t[dimused];
	vars[0]=cent;
	vars[1]=vtx;
	vars[2]=trigpt;
	vars[3]=asso->Pt();
	vars[4]=deleta;
	vars[5]=delphi;

	Int_t dimension=6;
        if(fRequestEventPlane) 
	{
       vars[6]=gPsiMinusPhiBin;
       dimension=7;
	}

if(!fcontainPIDtrig && !fcontainPIDasso && SetChargeAxis==2){
        vars[dimension]=trig->Charge();
	vars[dimension+1]=asso->Charge();
 }
if(fcontainPIDtrig && !fcontainPIDasso){
        vars[dimension]=ParticlePID_InvMass;
if(SetChargeAxis==2){
        vars[dimension+1]=trig->Charge();
	vars[dimension+2]=asso->Charge();
 }
	}
if(!fcontainPIDtrig && fcontainPIDasso){
        vars[dimension]=particlepidasso;
if(SetChargeAxis==2){
        vars[dimension+1]=trig->Charge();
	vars[dimension+2]=asso->Charge();
   }
 }
 if(fcontainPIDtrig && fcontainPIDasso){
	vars[dimension]=ParticlePID_InvMass;
	vars[dimension+1]=particlepidasso;
if(SetChargeAxis==2){
        vars[dimension+2]=trig->Charge();
	vars[dimension+3]=asso->Charge();
   }
 }

	if (fWeightPerEvent)
	{
	  Int_t weightBin = triggerWeighting->GetXaxis()->FindBin(vars[2]);
// 	  Printf("Using weight %f", triggerWeighting->GetBinContent(weightBin));
	 effweight *= triggerWeighting->GetBinContent(weightBin);
	}
    

	//Fill Correlation ThnSparses
    if(fillup=="trigassoUNID")
      {
    if(mixcase==kFALSE)  fTHnCorrUNID->Fill(vars,0,1.0/effweight);
    if(mixcase==kTRUE)   fTHnCorrUNIDmix->Fill(vars,0,1.0/effweight);
      }
    if(fillup=="trigIDassoID")
      {
	if(mixcase==kFALSE)  fTHnCorrID->Fill(vars,0,1.0/effweight);
	if(mixcase==kTRUE)  fTHnCorrIDmix->Fill(vars,0,1.0/effweight);
      }
    if(fillup=="trigIDassoIDMCTRUTH")//******************************************for TRUTH MC particles
      {//in this case effweight should be 1 always 
	if(mixcase==kFALSE)  fCorrelatonTruthPrimary->Fill(vars,0,1.0/effweight); 
	if(mixcase==kTRUE) fCorrelatonTruthPrimarymix->Fill(vars,0,1.0/effweight);
    }   
    if(fillup=="trigIDassoUNID" || fillup=="trigUNIDassoID")//****************************be careful
      {
	if(mixcase==kFALSE)  fTHnCorrIDUNID->Fill(vars,0,1.0/effweight);
	if(mixcase==kTRUE)   fTHnCorrIDUNIDmix->Fill(vars,0,1.0/effweight);
       }
	
delete[] vars;
   }//asso loop ends 
delete[] trigval;
 }//trigger loop ends 

 if (triggerWeighting)
    {
      delete triggerWeighting;
      triggerWeighting = 0;
    }
}

//------------------------------------------------------------------------------------------------
Bool_t AliTwoParticlePIDCorr:: CalculateSharedFraction(const TBits *triggerPadMap,const TBits *assocPadMap,const TBits *triggerShareMap,const TBits *assocShareMap)
{//source code-AliFemtoShareQualityPairCut.cxx
Double_t nofhits=0;
Double_t nofsharedhits=0;

for(UInt_t imap=0;imap< (triggerPadMap->GetNbits() );imap++)
{
//if they are in same pad
//cout<<triggerPadMap->TestBitNumber(imap)<<"    "<< assocPadMap->TestBitNumber(imap)<<endl;
if (triggerPadMap->TestBitNumber(imap) &&
      assocPadMap->TestBitNumber(imap))
{
//if they share
//cout<<triggerShareMap->TestBitNumber(imap)<<"   "<<assocShareMap->TestBitNumber(imap)<<endl;
if (triggerShareMap->TestBitNumber(imap) &&
      assocShareMap->TestBitNumber(imap))
{
  //cout<<triggerShareMap->TestBitNumber(imap)<<"   "<<assocShareMap->TestBitNumber(imap)<<endl;
nofhits+=2;
nofsharedhits+=2;
}



//not shared
    else {
     
      nofhits+=2;
    }


}
//different pad

//cout<< (triggerPadMap->TestBitNumber(imap) || assocPadMap->TestBitNumber(imap))<<endl;
else if (triggerPadMap->TestBitNumber(imap) ||
      assocPadMap->TestBitNumber(imap)) {
    // One track has a hit, the other does not
   
    nofhits++;
    //cout<<"No hits :"<<nofhits<<endl;
   
      }



}

Double_t SharedFraction=0.0;
if(nofhits>0) SharedFraction=(nofsharedhits/nofhits);

//cout<<"Fraction shared hits :"<<SharedFraction<<endl;

if(SharedFraction>fSharedfraction_Pair_cut) return kFALSE;

return kTRUE;

}

//________________________________________________________________________________________________
Float_t AliTwoParticlePIDCorr::GetTrackbyTrackeffvalue(AliAODTrack* track,Double_t cent,Float_t evzvtx, Int_t parpid)
{
  //This function is called only when applyefficiency=kTRUE; also ensure that "track" is present before calling that function
 Int_t effVars[4];
 Float_t effvalue=1.; 

  if(parpid==unidentified)
            {
	    effVars[0] = effcorection[5]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[5]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[5]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[5]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[5]->GetBinContent(effVars);
	    }
if(parpid==SpPion || parpid==SpKaon)
            {
	      if(fmesoneffrequired && !fkaonprotoneffrequired && track->Pt()>=fminPtComboeff && track->Pt()<=fmaxPtComboeff)
		{
	    effVars[0] = effcorection[3]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[3]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[3]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[3]->GetAxis(3)->FindBin(track->Eta());
            effvalue=effcorection[3]->GetBinContent(effVars);
		}
	      else{
 if(parpid==SpPion)
            {
	    effVars[0] = effcorection[0]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[0]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[0]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[0]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[0]->GetBinContent(effVars);
	    }
	    
 if(parpid==SpKaon)
            {
	    effVars[0] = effcorection[1]->GetAxis(0)->FindBin(cent);
	    effVars[1] =  effcorection[1]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] =  effcorection[1]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] =  effcorection[1]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[1]->GetBinContent(effVars);
	    }
	      }
	    }	
	     
 if(parpid==SpProton)
            {
	    effVars[0] =  effcorection[2]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[2]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[2]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[2]->GetAxis(3)->FindBin(track->Eta()); 
            effvalue=effcorection[2]->GetBinContent(effVars);
	    }

 if(fkaonprotoneffrequired && !fmesoneffrequired && track->Pt()>=fminPtComboeff && track->Pt()<=fmaxPtComboeff)
		{
  if(parpid==SpProton || parpid==SpKaon)
            {
	    effVars[0] = effcorection[4]->GetAxis(0)->FindBin(cent);
	    effVars[1] = effcorection[4]->GetAxis(1)->FindBin(evzvtx); 
	    effVars[2] = effcorection[4]->GetAxis(2)->FindBin(track->Pt()); 
	    effVars[3] = effcorection[4]->GetAxis(3)->FindBin(track->Eta());
            effvalue=effcorection[4]->GetBinContent(effVars);
	    }
		}	    
	    // 	  Printf("%d %d %d %d %f", effVars[0], effVars[1], effVars[2], effVars[3], fEfficiencyCorrectionAssociated->GetBinContent(effVars));
     if(effvalue<=0.) effvalue=1.;

     return effvalue; 

}
//---------------------------------------------------------------------------------



Float_t AliTwoParticlePIDCorr::GetV0_MeanSigma_CentPt(Double_t cent, Float_t V0Pt, Int_t parpid)
{

 Int_t effVars[2];
 Float_t effvalue=1.; 

  if(parpid==0)
            {
	      effVars[0] = kShortSigmahisto->GetXaxis()->FindBin(cent);
	      effVars[1] = kShortSigmahisto->GetYaxis()->FindBin(V0Pt); 
            effvalue=kShortSigmahisto->GetBinContent(effVars[0],effVars[1]);
	    }

if(parpid==1)
            {
	      effVars[0] = LambdaSigmahisto->GetXaxis()->FindBin(cent);
	      effVars[1] = LambdaSigmahisto->GetYaxis()->FindBin(V0Pt); 
	      effvalue=LambdaSigmahisto->GetBinContent(effVars[0],effVars[1]);
	    }

 if(effvalue<=0.) {

   if  (fSampleType=="PbPb") {
     if(parpid==0)   effvalue=0.007;//average no. , no meaning , just for protection
     if(parpid==1)   effvalue=0.0024;//average no. , no meaning , just for protection
   }

   if  (fSampleType=="pPb") {
     if(parpid==0)   effvalue=0.0055;//average no. , no meaning , just for protection//Check and change for pPb
     if(parpid==1)   effvalue=0.0021;//average no. , no meaning , just for protection
   }
 }
       return effvalue; 


}


//____________________________________________________________________________________________
Int_t AliTwoParticlePIDCorr::ClassifyTrack(AliAODTrack* track,AliAODVertex* vertex,Float_t magfield, Bool_t fill)
{  
 
  if(!track) return 0;
  Bool_t trackOK = track->TestFilterBit(fFilterBit);
  if(!trackOK) return 0;
  if (fTrackStatus != 0 && !CheckTrack(track)) return 0;
  //select only primary traks(for data & reco MC tracks) 
  if(fonlyprimarydatareco && track->GetType()!=AliAODTrack::kPrimary) return 0;
  if(track->Charge()==0) return 0;
  if (fill) fHistQA[12]->Fill(track->GetTPCNcls());  
  Float_t dxy, dz;		  
  dxy = track->DCA();
  dz = track->ZAtDCA();
  if (fill) fHistQA[6]->Fill(dxy);
  if (fill) fHistQA[7]->Fill(dz);
  Float_t chi2ndf = track->Chi2perNDF();
  if (fill) fHistQA[13]->Fill(chi2ndf);  
  // Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  Float_t nCrossedRowsTPC = track->GetTPCNCrossedRows();
  if (fill) fHistQA[14]->Fill(nCrossedRowsTPC); 
  //Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  if (track->GetTPCNclsF()>0) {
   Float_t  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
   if (fill) fHistQA[15]->Fill(ratioCrossedRowsOverFindableClustersTPC);
    }
//accepted tracks  
     Float_t pt=track->Pt();
     if(pt< fminPt || pt> fmaxPt) return 0;
     if(TMath::Abs(track->Eta())> fmaxeta) return 0;
     if(track->Phi()<0. || track->Phi()>2*TMath::Pi()) return 0;
     //if (!HasTPCPID(track)) return 0;//trigger & associated particles must have TPC PID,Is it required???
     
 //===========================PID (so far only for electron rejection)===============================//		    
     if(fElectronRejection) {//adapted from TaskBFpsi

	// get the electron nsigma
	Double_t nSigma = TMath::Abs(fPID->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron));
	
	//Fill QA before the PID
	if (fill){
	fHistdEdxVsPTPCbeforePIDelectron -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
	fHistNSigmaTPCvsPtbeforePIDelectron -> Fill(track->P()*track->Charge(),nSigma);} 
	//end of QA-before pid
	
	// check only for given momentum range
	if( pt > fElectronRejectionMinPt && pt < fElectronRejectionMaxPt ){
	  	  
	  //look only at electron nsigma
	  if(fElectronOnlyRejection){
	    
	    //Make the decision based on the n-sigma of electrons
	    if(nSigma < fElectronRejectionNSigma) return 0;
	  }
	  else{
	    
	    Double_t nSigmaPions   = TMath::Abs(fPID->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion));
	    Double_t nSigmaKaons   = TMath::Abs(fPID->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon));
	    Double_t nSigmaProtons = TMath::Abs(fPID->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kProton));
	    
	    //Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
	    if(nSigma < fElectronRejectionNSigma
	       && nSigmaPions   > fElectronRejectionNSigma
	       && nSigmaKaons   > fElectronRejectionNSigma
	       && nSigmaProtons > fElectronRejectionNSigma ) return 0;
	  }
	}
	
      }
      //===========================end of PID (so far only for electron rejection)===============================//
     
// DCA XY cut to check contaminations from lambda using dca sigma cut variation with filterbit 16, only for identified particles
//But if we put the cut here then it will be for all (identified and unidentified) particles. BUT this cut is mainly used for identified particles.
     
     /*if (if(track->Pt() > fminPtTrig) && fDCAXYCut)
	{
	  if (!vertex)
	    return 0;
	  
	  Double_t pos[2];
	  Double_t covar[3];
	  AliAODTrack* clone =(AliAODTrack*) track->Clone();
	  Bool_t success = clone->PropagateToDCA(vertex, magfield, fdcacutvalue, pos, covar);
	  delete clone;
	  if (!success)
	    return 0;

// 	  Printf("%f", ((AliAODTrack*)part)->DCA());
// 	  Printf("%f", pos[0]);
	  if (TMath::Abs(pos[0]) > fDCAXYCut->Eval(track->Pt()))
	    return 0;
	}*/

	if (fSharedClusterCut >= 0)
	{
	  Double_t frac = Double_t(((AliAODTrack*)track)->GetTPCnclsS()) / Double_t(((AliAODTrack*)track)->GetTPCncls());
	  if (frac > fSharedClusterCut)
	    return 0;
	}

   // Rejects tracks with shared clusters after filling a control histogram
   // This overload is used for primaries

     // Get the shared maps
      const TBits sharedMap = track->GetTPCSharedMap();
     // Fill a control histogram
      fPriHistShare->Fill(sharedMap.CountBits());

    // Reject shared clusters
       if (fSharedTPCmapCut >= 0)
	{     
      if((sharedMap.CountBits()) >= 1)  return 0;// Bad track, has too many shared clusters!
	}

     if (fill) fHistQA[8]->Fill(pt);
     if (fill) fHistQA[9]->Fill(track->Eta());
     if (fill) fHistQA[10]->Fill(track->Phi());
     return 1;
  }
  //________________________________________________________________________________
void AliTwoParticlePIDCorr::CalculateNSigmas(AliAODTrack *track, Bool_t FIllQAHistos) 
{
//This function is called within the func GetParticle() for accepted tracks only i.e.after call of Classifytrack() & for those tracks which have proper TPC PID response . combined nsigma(circular) cut only for particles having pt upto  4.0 Gev/c and beyond that use the asymmetric nsigma cut around pion's mean position in TPC ( while filling the  TObjArray for trig & asso )
Float_t pt=track->Pt();

//plot the separation power

Double_t bethe[AliPID::kSPECIES]={0.};
Double_t sigma_TPC[AliPID::kSPECIES]={0.}; 

 Double_t Pi_Ka_sep[NSigmaPIDType+1]={0.};
 Double_t Pi_Pr_sep[NSigmaPIDType+1]={0.};
 Double_t Ka_Pr_sep[NSigmaPIDType+1]={0.};


    Double_t ptpc = track->GetTPCmomentum();
    Int_t dEdxN = track->GetTPCsignalN();
 for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
       bethe[ipart] = fPID->GetTPCResponse().GetExpectedSignal(ptpc, (AliPID::EParticleType)ipart);
      //Double_t diff = dEdx - bethe;
       sigma_TPC[ipart] = fPID->GetTPCResponse().GetExpectedSigma(ptpc, dEdxN, (AliPID::EParticleType)ipart);
      //nSigma[ipart] = diff / sigma;
    }
 Pi_Ka_sep[NSigmaTPC]=TMath::Abs(bethe[AliPID::kPion]-bethe[AliPID::kKaon])/((sigma_TPC[AliPID::kPion]+sigma_TPC[AliPID::kKaon])/2.0);
 Pi_Pr_sep[NSigmaTPC]=TMath::Abs(bethe[AliPID::kPion]-bethe[AliPID::kProton])/((sigma_TPC[AliPID::kPion]+sigma_TPC[AliPID::kProton])/2.0);
 Ka_Pr_sep[NSigmaTPC]=TMath::Abs(bethe[AliPID::kKaon]-bethe[AliPID::kProton])/((sigma_TPC[AliPID::kKaon]+sigma_TPC[AliPID::kProton])/2.0);


Double_t sigma_TOF[AliPID::kSPECIES]={0.}; 

if(HasTOFPID(track) && pt>fPtTOFPIDmin)
   {
 Double_t timei[AliPID::kSPECIES];
 track->GetIntegratedTimes(timei);
 for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {  sigma_TOF[ipart]= fPID->GetTOFResponse().GetExpectedSigma(track->P(), timei[ipart], AliPID::ParticleMass(ipart));}
 Pi_Ka_sep[NSigmaTOF]=TMath::Abs(timei[AliPID::kPion]-timei[AliPID::kKaon])/((sigma_TOF[AliPID::kPion]+sigma_TOF[AliPID::kKaon])/2.0);
 Pi_Pr_sep[NSigmaTOF]=TMath::Abs(timei[AliPID::kPion]-timei[AliPID::kProton])/((sigma_TOF[AliPID::kPion]+sigma_TOF[AliPID::kProton])/2.0);
 Ka_Pr_sep[NSigmaTOF]=TMath::Abs(timei[AliPID::kKaon]-timei[AliPID::kProton])/((sigma_TOF[AliPID::kKaon]+sigma_TOF[AliPID::kProton])/2.0);

  Pi_Ka_sep[NSigmaTPCTOF]=TMath::Abs(Pi_Ka_sep[NSigmaTPC]*Pi_Ka_sep[NSigmaTPC]+Pi_Ka_sep[NSigmaTOF]*Pi_Ka_sep[NSigmaTOF]);
  Pi_Pr_sep[NSigmaTPCTOF]=TMath::Abs(Pi_Pr_sep[NSigmaTPC]*Pi_Pr_sep[NSigmaTPC]+Pi_Pr_sep[NSigmaTOF]*Pi_Pr_sep[NSigmaTOF]);
  Ka_Pr_sep[NSigmaTPCTOF]=TMath::Abs(Ka_Pr_sep[NSigmaTPC]*Ka_Pr_sep[NSigmaTPC]+Ka_Pr_sep[NSigmaTOF]*Ka_Pr_sep[NSigmaTOF]);
   }


//fill separation power histograms
 for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
   if(ipid==0){
	TH2F *h=GetHistogram2D(Form("Pi_Ka_sep_%d",ipid));
	h->Fill(track->Pt(),Pi_Ka_sep[ipid]);
        TH2F *h1=GetHistogram2D(Form("Pi_Pr_sep_%d",ipid));
	h1->Fill(track->Pt(),Pi_Pr_sep[ipid]);
        TH2F *h2=GetHistogram2D(Form("Ka_Pr_sep_%d",ipid));
	h2->Fill(track->Pt(),Ka_Pr_sep[ipid]);
   }
   if(HasTOFPID(track) && pt>fPtTOFPIDmin && ipid!=0){
       TH2F *h=GetHistogram2D(Form("Pi_Ka_sep_%d",ipid));
	h->Fill(track->Pt(),Pi_Ka_sep[ipid]);
        TH2F *h1=GetHistogram2D(Form("Pi_Pr_sep_%d",ipid));
	h1->Fill(track->Pt(),Pi_Pr_sep[ipid]);
        TH2F *h2=GetHistogram2D(Form("Ka_Pr_sep_%d",ipid));
	h2->Fill(track->Pt(),Ka_Pr_sep[ipid]);
   }
 }


//it is assumed that every track that passed the filterbit have proper TPC response(!!)
Float_t nsigmaTPCkPion =fPID->NumberOfSigmasTPC(track, AliPID::kPion);
Float_t nsigmaTPCkKaon =fPID->NumberOfSigmasTPC(track, AliPID::kKaon);
Float_t nsigmaTPCkProton =fPID->NumberOfSigmasTPC(track, AliPID::kProton);

Float_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;
Float_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;

 if(HasTOFPID(track) && pt>fPtTOFPIDmin && pt<fPtTOFPIDmax)
   {

nsigmaTOFkPion =fPID->NumberOfSigmasTOF(track, AliPID::kPion);
nsigmaTOFkKaon =fPID->NumberOfSigmasTOF(track, AliPID::kKaon);
nsigmaTOFkProton =fPID->NumberOfSigmasTOF(track, AliPID::kProton);
//---combined
nsigmaTPCTOFkPion   = TMath::Sqrt(nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion);
nsigmaTPCTOFkKaon   = TMath::Sqrt(nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon);
nsigmaTPCTOFkProton = TMath::Sqrt(nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton);


   }
else{
    // --- combined
    // if TOF is missing and below fPtTOFPID only the TPC information is used
    nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
    nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
    nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);

  }

//set data member fnsigmas
  fnsigmas[SpPion][NSigmaTPC]=nsigmaTPCkPion;
  fnsigmas[SpKaon][NSigmaTPC]=nsigmaTPCkKaon;
  fnsigmas[SpProton][NSigmaTPC]=nsigmaTPCkProton;

  //for all tracks below fPtTOFPIDmin  and also for tracks above fPtTOFPIDmin without proper TOF response these TOF nsigma values will be 999.
  fnsigmas[SpPion][NSigmaTOF]=nsigmaTOFkPion;
  fnsigmas[SpKaon][NSigmaTOF]=nsigmaTOFkKaon;
  fnsigmas[SpProton][NSigmaTOF]=nsigmaTOFkProton;

 //for all tracks below fPtTOFPIDmin  and also for tracks above fPtTOFPIDmin without proper TOF response these TPCTOF nsigma values will be TMath::Abs(TPC only nsigma)
  fnsigmas[SpPion][NSigmaTPCTOF]=nsigmaTPCTOFkPion;
  fnsigmas[SpKaon][NSigmaTPCTOF]=nsigmaTPCTOFkKaon;
  fnsigmas[SpProton][NSigmaTPCTOF]=nsigmaTPCTOFkProton;

 if(FIllQAHistos){
    //Fill NSigma SeparationPlot
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigma_%d_%d",ipart,ipid));
	h->Fill(track->Pt(),fnsigmas[ipart][ipid]);
      }
    }
  }

}
//----------------------------------------------------------------------------
Int_t AliTwoParticlePIDCorr::FindMinNSigma(AliAODTrack *track,Bool_t FillQAHistos) 
{
  //this function is always called after calling the function CalculateNSigmas(AliAODTrack *track)
if(fRequestTOFPID && track->Pt()>fPtTOFPIDmin && track->Pt()<fPtTOFPIDmax && (!HasTOFPID(track)) )return SpUndefined;//so any track having Pt>0.5 withot having proper TOF response will be defined as SpUndefined
//get the identity of the particle with the minimum Nsigma
  Float_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fPIDType){
  case NSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPC]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPC])  ;
    break;
  case NSigmaTOF:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTOF])  ;
    break;
  case NSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF])  ;
    break;
  case Bayes://the nsigma in the bayesian is used to clean with a very large n-sigma value
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF])  ;
    break;
  }


if(fdiffPIDcutvalues){
  if(track->Pt()<=fPt1) fNSigmaPID=fPIDCutval1;
  if(track->Pt()>fPt1 && track->Pt()<=fPt2) fNSigmaPID=fPIDCutval2;
  if(track->Pt()>fPt2 && track->Pt()<=fPt3) fNSigmaPID=fPIDCutval3;
  if(track->Pt()>fPt3) fNSigmaPID=fPIDCutval4;
  }

 // guess the particle based on the smaller nsigma (within fNSigmaPID)
  if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) return SpUndefined;//it is the default value for the three

  if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon < fNSigmaPID)){
    if((fHighPtKaonNSigmaPID>0) && (track->Pt()>fHighPtKaonSigma) && (nsigmaKaon > fHighPtKaonNSigmaPID)) return SpUndefined;//different nsigma cut for kaons above a particular Pt range(within the TPC-TOF PID range)
if(FillQAHistos){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",SpKaon,ipid));
	h->Fill(track->Pt(),fnsigmas[SpKaon][ipid]);
      }
    }
 return SpKaon;
  }
  if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion < fNSigmaPID)) {
 if(FillQAHistos){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",SpPion,ipid));
	h->Fill(track->Pt(),fnsigmas[SpPion][ipid]);
      }
    }
return SpPion;
  }
  if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fNSigmaPID)) {
if(FillQAHistos){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",SpProton,ipid));
	h->Fill(track->Pt(),fnsigmas[SpProton][ipid]);
      }
 }
return SpProton;
  }

// else, return undefined
  return SpUndefined;
  
 
}

//------------------------------------------------------------------------------------------
Bool_t* AliTwoParticlePIDCorr::GetDoubleCounting(AliAODTrack * trk,Bool_t FIllQAHistos){ 
  //this function is always called after calling the function CalculateNSigmas(AliAODTrack *track)

  //if a particle has double counting set fHasDoubleCounting[ipart]=kTRUE
  //fill DC histos
  for(Int_t ipart=0;ipart<NSpecies;ipart++)fHasDoubleCounting[ipart]=kFALSE;//array with kTRUE for second (or third) identity of the track
  
  Int_t MinNSigma=FindMinNSigma(trk,kFALSE);//not filling the NSigmaRec histos
  
  
  if(MinNSigma==SpUndefined)return fHasDoubleCounting;//in case of undefined no Double counting
  
  Float_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fPIDType) {
  case NSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPC]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPC])  ;
    break;
  case NSigmaTOF:
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTOF])  ;
    break;
  case NSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF])  ;
    break;
  case Bayes://the nsigma in the bayesian is used to clean with a very large n-sigma value
    nsigmaProton  =  TMath::Abs(fnsigmas[SpProton][NSigmaTPCTOF]);
    nsigmaKaon	  =  TMath::Abs(fnsigmas[SpKaon][NSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[SpPion][NSigmaTPCTOF])  ;
    break;
  }

  // Actually the tracks in the overlapping region(in TPC-TOF nSigma plane) will be ignored

  if(nsigmaPion<fNSigmaPID && MinNSigma!=SpPion)fHasDoubleCounting[SpPion]=kTRUE;
  if(nsigmaKaon<fNSigmaPID && MinNSigma!=SpKaon)fHasDoubleCounting[SpKaon]=kTRUE;
  if(nsigmaProton<fNSigmaPID && MinNSigma!=SpProton)fHasDoubleCounting[SpProton]=kTRUE;
     
  

if(FIllQAHistos){
    //fill NSigma distr for double counting
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
      if(fHasDoubleCounting[ipart]){//this may be kTRUE only for particles having Pt<=4.0 GeV/C, so this histo contains all the particles having Pt<=4.0 GeV/C in the nsigma overlapping region in TPC/TPC-TOF plane 
	for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	  if((ipid!=NSigmaTPC) && (!HasTOFPID(trk)))continue;//not filling TOF and combined if no TOF PID
	  TH2F *h=GetHistogram2D(Form("NSigmaDC_%d_%d",ipart,ipid));
	  h->Fill(trk->Pt(),fnsigmas[ipart][ipid]);
	}
      }
    }
  }
 
 
  return fHasDoubleCounting;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t* AliTwoParticlePIDCorr::GetAllCompatibleIdentitiesNSigma(AliAODTrack * trk,Bool_t FIllQAHistos){ 
 //mainly intended to check the probability of the PID of the tracks which are in the overlapping nSigma regions and near about the middle position from the   mean position of two ID particle
  Bool_t *IDs=GetDoubleCounting(trk,FIllQAHistos);
  IDs[FindMinNSigma(trk,FIllQAHistos)]=kTRUE;
  return IDs;
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////

UInt_t AliTwoParticlePIDCorr::CalcPIDCombined(AliAODTrack *track, Int_t detMask, Double_t* prob) const{
  //
  // Bayesian PID calculation
  //
  for(Int_t i=0;i<AliPID::kSPECIES;i++)
    {
      prob[i]=0.;
    }
  fPIDCombined->SetDetectorMask(detMask);
  
  return fPIDCombined->ComputeProbabilities((AliAODTrack*)track, fPID, prob);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTwoParticlePIDCorr::GetIDBayes(AliAODTrack * trk, Bool_t FIllQAHistos){ 
  
  Bool_t *IDs=GetAllCompatibleIdentitiesNSigma(trk,FIllQAHistos);


  //Filling of Probability histos
	Double_t probTPC[AliPID::kSPECIES]={0.};
	Double_t probTOF[AliPID::kSPECIES]={0.};
	Double_t probTPCTOF[AliPID::kSPECIES]={0.};

	UInt_t detUsedTPC = 0;
	UInt_t detUsedTOF = 0;
	UInt_t detUsedTPCTOF = 0;

 //get the TPC probability
          fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	  detUsedTPC = fPIDCombined->ComputeProbabilities(trk, fPID, probTPC);
if(detUsedTPC == AliPIDResponse::kDetTPC)
  {
for(Int_t ipart=0;ipart<NSpecies;ipart++){

	TH2F *h=GetHistogram2D(Form("probBayes_TPC_%d",ipart));
	if(ipart==0)	h->Fill(trk->Pt(),probTPC[AliPID::kPion]);
	if(ipart==1)	h->Fill(trk->Pt(),probTPC[AliPID::kKaon]);
	if(ipart==2)	h->Fill(trk->Pt(),probTPC[AliPID::kProton]);
 }
  }

  //get the TOF probability
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	  detUsedTOF = fPIDCombined->ComputeProbabilities(trk, fPID, probTOF);
if(detUsedTOF == AliPIDResponse::kDetTOF)
  {
for(Int_t ipart=0;ipart<NSpecies;ipart++){
	TH2F *h=GetHistogram2D(Form("probBayes_TOF_%d",ipart));
	if(ipart==0)	h->Fill(trk->Pt(),probTOF[AliPID::kPion]);
	if(ipart==1)	h->Fill(trk->Pt(),probTOF[AliPID::kKaon]);
	if(ipart==2)	h->Fill(trk->Pt(),probTOF[AliPID::kProton]);
 }
  }

 //get the TPC-TOF probability
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  detUsedTPCTOF = fPIDCombined->ComputeProbabilities(trk, fPID, probTPCTOF);
if(detUsedTPCTOF == (AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC))
  {
for(Int_t ipart=0;ipart<NSpecies;ipart++){
	TH2F *h=GetHistogram2D(Form("probBayes_TPCTOF_%d",ipart));
	if(ipart==0)	h->Fill(trk->Pt(),probTPCTOF[AliPID::kPion]);
	if(ipart==1)	h->Fill(trk->Pt(),probTPCTOF[AliPID::kKaon]);
	if(ipart==2)	h->Fill(trk->Pt(),probTPCTOF[AliPID::kProton]); 
}
  }

//bayesian PID is basically calculated using TPC-TOF cpmbined info after a particular Pt by default
  Double_t probBayes[AliPID::kSPECIES];
  
  UInt_t detUsed= 0;
  if(trk->Pt()>fPtTOFPIDmin && HasTOFPID(trk)){//use TOF information
    detUsed = CalcPIDCombined(trk, AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC, probBayes);
    if(detUsed != (AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC))return SpUndefined;//check that TPC and TOF are used
  }else{
    detUsed = CalcPIDCombined(trk,AliPIDResponse::kDetTPC, probBayes);
    if(detUsed != AliPIDResponse::kDetTPC)return SpUndefined;//check that TPC is used
  }
  
  //the probability has to be normalized to one, we check it
  Double_t sump=0.;
  for(Int_t ipart=0;ipart<AliPID::kSPECIES;ipart++)sump+=probBayes[ipart];
  if(sump<.99 && sump>1.01){//FIXME precision problem in the sum, workaround
    AliFatal("Bayesian probability not normalized to one");
  }

  if(fdiffPIDcutvalues){
  if(trk->Pt()<=fPt1) fBayesCut=fPIDCutval1;
  if(trk->Pt()>fPt1 && trk->Pt()<=fPt2) fBayesCut=fPIDCutval2;
  if(trk->Pt()>fPt2 && trk->Pt()<=fPt3) fBayesCut=fPIDCutval3;
  if(trk->Pt()>fPt3) fBayesCut=fPIDCutval4;
  }

  
  //probabilities are normalized to one, if the cut is above .5 there is no problem
  if(probBayes[AliPID::kPion]>fBayesCut && IDs[SpPion]==1){
    TH2F *h=GetHistogram2D(Form("BayesRec_%d",SpPion));
    h->Fill(trk->Pt(),probBayes[AliPID::kPion]);
    return SpPion;
  }
  else if(probBayes[AliPID::kKaon]>fBayesCut && IDs[SpKaon]==1){
    TH2F *h=GetHistogram2D(Form("BayesRec_%d",SpKaon));
    h->Fill(trk->Pt(),probBayes[AliPID::kKaon]);
    return SpKaon;
  }
  else if(probBayes[AliPID::kProton]>fBayesCut && IDs[SpProton]==1){
    TH2F *h=GetHistogram2D(Form("BayesRec_%d",SpProton));
    h->Fill(trk->Pt(),probBayes[AliPID::kProton]);
    return SpProton;
  }
  else{
    return SpUndefined;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTwoParticlePIDCorr::GetParticle(AliAODTrack * trk, Bool_t FIllQAHistos){ 
  //return the specie according to the minimum nsigma value
  //no double counting, this has to be evaluated using CheckDoubleCounting()
  
  Int_t ID=SpUndefined;

  //calculate nsigmas (used also by the bayesian)
  CalculateNSigmas(trk,FIllQAHistos);//fill the data member fnsigmas with the nsigmas value [ipart][iPID]


 //Do PID eith TPC+TOF(combined TPC+TOFNSigma/TPCNSigma for low Pt or bayesian)
  if(trk->Pt()<fPtTOFPIDmax){//use any pidmethod mentioned in the enum in .h
  if(fPIDType==Bayes){//use bayesianPID
    
    if(!fPIDCombined) {
      AliFatal("PIDCombined object not found");
    }
    
    ID = GetIDBayes(trk,FIllQAHistos);
    
  }else{ //use nsigma PID  

   ID=FindMinNSigma(trk,FIllQAHistos);
if(fUseExclusiveNSigma){ //if one particle has double counting and exclusive nsigma is requested ID = kSpUndefined
      Bool_t *HasDC;
      HasDC=GetDoubleCounting(trk,FIllQAHistos);
      for(Int_t ipart=0;ipart<NSpecies;ipart++){
	if(HasDC[ipart]==kTRUE)  ID = SpUndefined;
      }
    }
  }
}
  else{//use the deltapion method with rel. rise in TPC only for trk->Pt()>=fPtTOFPIDmax
    Double_t delta_pion=fPID->GetSignalDelta(AliPIDResponse::kTPC,trk, AliPID::kPion, kFALSE);
     deltapion_val=delta_pion;//global variable, can be used anywhere in the code after calling the function GetParticle()
     	//before deltapion cut
    for(Int_t ipart=0;ipart<NSpecies;ipart++){
	  if(ipart!=SpPion) continue;//only around pion's mean position;same histo for Pr and Pi
     TH2F *h=GetHistogram2D(Form("deltapion_%d",ipart));
    h->Fill(trk->Pt(),delta_pion);
    }
    //Now apply deltapion cut
    if (delta_pion>fPiondeltacutmin && delta_pion<fPiondeltacutmax) {
	TH2F *h1=GetHistogram2D(Form("NSigmaRec_%d_%d",SpPion,NSigmaTPC));
	h1->Fill(trk->Pt(),fPID->NumberOfSigmasTPC(trk, AliPID::kPion));

	TH2F *h2=GetHistogram2D(Form("deltapionRec_%d",SpPion));
        h2->Fill(trk->Pt(),delta_pion);
        ID = SpPion;
    }
     //else if (delta_pion>-2 && delta_pion<0) ID = SpKaon;//no meaning, we can't identify kaons with this method in the high Pt(rel. rise in TPC) range
    else if (delta_pion>fProtondeltacutmin && delta_pion<fProtondeltacutmax){
        TH2F *h3=GetHistogram2D(Form("NSigmaRec_%d_%d",SpProton,NSigmaTPC));
	h3->Fill(trk->Pt(),fPID->NumberOfSigmasTPC(trk, AliPID::kProton));

	TH2F *h4=GetHistogram2D(Form("deltapionRec_%d",SpPion));
        h4->Fill(trk->Pt(),delta_pion);
      ID = SpProton;
    }
    else ID = SpUndefined;
   
   // The TPC PID is based on the observable ∆π = dE/dx − ⟨dE/dx⟩π. Particles with 1 < ∆π <7and−21<∆π <−14(−25<∆π <−14)forpass1(pass2)wereidentified as pions and protons, respectively
  }//deltapion cut for trk->Pt()>fPtTOFPIDmax

 if(FIllQAHistos){
//Fill PID signal plot
  if(ID != SpUndefined){
    for(Int_t idet=0;idet<fNDetectors;idet++){
      TH2F *h=GetHistogram2D(Form("PID_%d_%d",idet,ID));
      if(idet==fITS)h->Fill(trk->P(),trk->GetITSsignal()*trk->Charge());
      if(idet==fTPC)h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
      if(idet==fTOF && HasTOFPID(trk))h->Fill(trk->P(),GetBeta(trk)*trk->Charge());
    }
  }
  //Fill PID signal plot without cuts
  for(Int_t idet=0;idet<fNDetectors;idet++){
    TH2F *h=GetHistogram2D(Form("PIDAll_%d",idet));
    if(idet==fITS)h->Fill(trk->P(),trk->GetITSsignal()*trk->Charge());
    if(idet==fTPC)h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
    if(idet==fTOF && HasTOFPID(trk))h->Fill(trk->P(),GetBeta(trk)*trk->Charge());
  }
 }
  return ID;
}

//-------------------------------------------------------------------------------------
Bool_t
AliTwoParticlePIDCorr::HasTPCPID(AliAODTrack *track) const
{  
  // check PID signal 
   AliPIDResponse::EDetPidStatus statustpc = fPID->CheckPIDStatus(AliPIDResponse::kTPC,track);
   if(statustpc!=AliPIDResponse::kDetPidOk) return kFALSE;
   //ULong_t status=track->GetStatus();
   //if  (!( (status & AliAODTrack::kTPCpid  ) == AliAODTrack::kTPCpid  )) return kFALSE;//remove light nuclei
   //if (track->GetTPCsignal() <= 0.) return kFALSE;
    if(track->GetTPCsignalN() < fNclsusedfordEdXdtr) return kFALSE;//tracks with TPCsignalN< 60 have questionable dEdx,cutting on TPCsignalN > 70 or > 60 shouldn't make too much difference in statistics,also  it is IMO safe to use TPC also for MIPs.
   
  return kTRUE;  
}
//___________________________________________________________

Bool_t
AliTwoParticlePIDCorr::HasTOFPID(AliAODTrack *track) const
{
  // check TOF matched track 
  //ULong_t status=track->GetStatus();
  //if  (!( (status & AliAODTrack::kITSin  ) == AliAODTrack::kITSin  )) return kFALSE;
 AliPIDResponse::EDetPidStatus statustof = fPID->CheckPIDStatus(AliPIDResponse::kTOF,track);
 if(statustof!= AliPIDResponse::kDetPidOk) return kFALSE;
  if(track->Pt()<=fPtTOFPIDmin || track->Pt()>=fPtTOFPIDmax) return kFALSE;
 //if(!((status & AliAODTrack::kTOFpid  ) == AliAODTrack::kTOFpid  )) return kFALSE;
 //Float_t probMis = fPIDresponse->GetTOFMismatchProbability(track);
 // if (probMis > 0.01) return kFALSE;
if(fRemoveTracksT0Fill)
    {
Int_t startTimeMask = fPID->GetTOFResponse().GetStartTimeMask(track->P());
      if (startTimeMask < 0)return kFALSE; 
    }
  return kTRUE;
}

//________________________________________________________________________
Double_t AliTwoParticlePIDCorr :: GetBeta(AliAODTrack *track)
{
  //it is called only when TOF PID is available
//TOF beta calculation
  Double_t tofTime=track->GetTOFsignal();
  
  Double_t c=TMath::C()*1.E-9;// m/ns
  Float_t startTime = fPID->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());//in ps
  Double_t length= fPID->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
  tofTime -= startTime;      // subtract startTime to the signal
  Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
  tof=tof*c;
  return length/tof;


  /*
  Double_t p = track->P();
  Double_t time=track->GetTOFsignal()-fPID->GetTOFResponse().GetStartTime(p);
  Double_t timei[5];
  track->GetIntegratedTimes(timei);
  return timei[0]/time;
  */
}
//------------------------------------------------------------------------------------------------------

Float_t AliTwoParticlePIDCorr::GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared
  // same can be achieved, but with more computing time with
  /*TLorentzVector photon, p1, p2;
  p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
  p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
  photon = p1+p2;
  photon.M()*/
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
    tantheta1 = 2 * TMath::Exp(-eta1) / ( 1 - TMath::Exp(-2*eta1));
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
    tantheta2 = 2 * TMath::Exp(-eta2) / ( 1 - TMath::Exp(-2*eta2));
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2 ) ) );
  
  return mass2;
}
//---------------------------------------------------------------------------------

Float_t AliTwoParticlePIDCorr::GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared approximately
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = 1.0-eta2+eta2*eta2/2-eta2*eta2*eta2/6+eta2*eta2*eta2*eta2/24;
    tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  // fold onto 0...pi
  Float_t deltaPhi = TMath::Abs(phi1 - phi2);
  while (deltaPhi > TMath::TwoPi())
    deltaPhi -= TMath::TwoPi();
  if (deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  
  Float_t cosDeltaPhi = 0;
  if (deltaPhi < TMath::Pi()/3)
    cosDeltaPhi = 1.0 - deltaPhi*deltaPhi/2 + deltaPhi*deltaPhi*deltaPhi*deltaPhi/24;
  else if (deltaPhi < 2*TMath::Pi()/3)
    cosDeltaPhi = -(deltaPhi - TMath::Pi()/2) + 1.0/6 * TMath::Power((deltaPhi - TMath::Pi()/2), 3);
  else
    cosDeltaPhi = -1.0 + 1.0/2.0*(deltaPhi - TMath::Pi())*(deltaPhi - TMath::Pi()) - 1.0/24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( cosDeltaPhi + 1.0 / tantheta1 / tantheta2 ) ) );
  
//   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
  
  return mass2;
}
//--------------------------------------------------------------------------------
Float_t  AliTwoParticlePIDCorr::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{ 
  //
  // calculates dphistar
  //
  
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}

//------------------------------------------------------------------------
Double_t* AliTwoParticlePIDCorr::GetBinning(const char* configuration, const char* tag, Int_t& nBins)
{
  // This method is a copy from AliUEHist::GetBinning
  // takes the binning from <configuration> identified by <tag>
  // configuration syntax example:
  // eta: 2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4
  // phi: .....
  //
  // returns bin edges which have to be deleted by the caller
  
  TString config(configuration);
  TObjArray* lines = config.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    if (line.BeginsWith(TString(tag) + ":"))
    {
      line.Remove(0, strlen(tag) + 1);
      line.ReplaceAll(" ", "");
      TObjArray* binning = line.Tokenize(",");
      Double_t* bins = new Double_t[binning->GetEntriesFast()];
      for (Int_t j=0; j<binning->GetEntriesFast(); j++)
	bins[j] = TString(binning->At(j)->GetName()).Atof();
      
      nBins = binning->GetEntriesFast() - 1;

      delete binning;
      delete lines;
      return bins;
    }
  }
  
  delete lines;
  AliFatal(Form("Tag %s not found in %s", tag, configuration));
  return 0;
}

//____________________________________________________________________

Bool_t AliTwoParticlePIDCorr::CheckTrack(AliAODTrack * part)
{
  // check if the track status flags are set
  
  UInt_t status=((AliAODTrack*)part)->GetStatus();
  if ((status & fTrackStatus) == fTrackStatus)
    return kTRUE;
  return kFALSE;
}
//________________________________________________________________________

Bool_t AliTwoParticlePIDCorr::AcceptEventCentralityWeight(Double_t centrality)
{
  // rejects "randomly" events such that the centrality gets flat
  // uses fCentralityWeights histogram

  // TODO code taken and adapted from AliRDHFCuts; waiting for general class AliCentralityFlattening
  
  Double_t weight = fCentralityWeights->GetBinContent(fCentralityWeights->FindBin(centrality));
  Double_t centralityDigits = centrality*100. - (Int_t)(centrality*100.);
  
  Bool_t result = kFALSE;
  if (centralityDigits < weight) 
    result = kTRUE;
  
  AliInfo(Form("Centrality: %f; Digits: %f; Weight: %f; Result: %d", centrality, centralityDigits, weight, result));
  
  return result;
}

//____________________________________________________________________
void AliTwoParticlePIDCorr::ShiftTracks(TObjArray* tracks, Double_t angle)
{
  // shifts the phi angle of all tracks by angle
  // 0 <= angle <= 2pi
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); ++i) 
  {
   LRCParticlePID *part=(LRCParticlePID*)(tracks->UncheckedAt(i));

    Double_t newAngle = part->Phi() + angle; 
    if (newAngle >= TMath::TwoPi())
      newAngle -= TMath::TwoPi();
    
    part->SetPhi(newAngle);
  }
}


//________________________________________________________________________
void  AliTwoParticlePIDCorr::SetVZEROCalibrationFile(const char* filename,const char* lhcPeriod) {
  //Function to setup the VZERO gain equalization
    //============Get the equilization map============//
  TFile *calibrationFile = TFile::Open(filename);
  if((!calibrationFile)||(!calibrationFile->IsOpen())) {
    Printf("No calibration file found!!!");
    return;
  }

  TList *list = dynamic_cast<TList *>(calibrationFile->Get(lhcPeriod));
  if(!list) {
    Printf("Calibration TList not found!!!");
    return;
  }

  fHistVZEROAGainEqualizationMap = dynamic_cast<TH1F *>(list->FindObject("gHistVZEROAGainEqualizationMap"));
  if(!fHistVZEROAGainEqualizationMap) {
    Printf("VZERO-A calibration object not found!!!");
    return;
  }
  fHistVZEROCGainEqualizationMap = dynamic_cast<TH1F *>(list->FindObject("gHistVZEROCGainEqualizationMap"));
  if(!fHistVZEROCGainEqualizationMap) {
    Printf("VZERO-C calibration object not found!!!");
    return;
  }

  fHistVZEROChannelGainEqualizationMap = dynamic_cast<TH2F *>(list->FindObject("gHistVZEROChannelGainEqualizationMap"));
  if(!fHistVZEROChannelGainEqualizationMap) {
    Printf("VZERO channel calibration object not found!!!");
    return;
  }
}

//________________________________________________________________________
Double_t AliTwoParticlePIDCorr::GetChannelEqualizationFactor(Int_t run,Int_t channel) {
  //
  if(!fHistVZEROAGainEqualizationMap) return 1.0;

  for(Int_t iBinX = 1; iBinX <= fHistVZEROChannelGainEqualizationMap->GetNbinsX(); iBinX++) {
    Int_t gRunNumber = atoi(fHistVZEROChannelGainEqualizationMap->GetXaxis()->GetBinLabel(iBinX));
    if(gRunNumber == run)
      return fHistVZEROChannelGainEqualizationMap->GetBinContent(iBinX,channel+1);
  }

  return 1.0;
}

//________________________________________________________________________
Double_t AliTwoParticlePIDCorr::GetEqualizationFactor(Int_t run, const char* side) {
  //
  if(!fHistVZEROAGainEqualizationMap) return 1.0;

  TString gVZEROSide = side;
  for(Int_t iBinX = 1; iBinX < fHistVZEROAGainEqualizationMap->GetNbinsX(); iBinX++) {
    Int_t gRunNumber = atoi(fHistVZEROAGainEqualizationMap->GetXaxis()->GetBinLabel(iBinX));
    //cout<<"Looking for run "<<run<<" - current run: "<<gRunNumber<<endl;
    if(gRunNumber == run) {
      if(gVZEROSide == "A") 
	return fHistVZEROAGainEqualizationMap->GetBinContent(iBinX);
      else if(gVZEROSide == "C") 
	return fHistVZEROCGainEqualizationMap->GetBinContent(iBinX);
    }
  }

  return 1.0;
}
//________________________________________________________________________
Double_t AliTwoParticlePIDCorr::GetReferenceMultiplicityVZEROFromAOD(AliVEvent *mainevent){
  //Function that returns the reference multiplicity from AODs (data or reco MC, Not for Truth)
  //Different ref. mult. implemented: V0M, V0A, V0C, TPC
  if(!mainevent) return -1;

      AliAODEvent* event = dynamic_cast<AliAODEvent*>(mainevent);

  Double_t gRefMultiplicity = 0., gRefMultiplicityTPC = 0.;
  Double_t gRefMultiplicityVZERO = 0., gRefMultiplicityVZEROA = 0., gRefMultiplicityVZEROC = 0.;

  AliAODHeader *header = dynamic_cast<AliAODHeader *>(event->GetHeader());
  if(!header) {
    Printf("ERROR: AOD header not available");
    return -999;
  }
  Int_t gRunNumber = header->GetRunNumber();
 Float_t bSign1=header->GetMagneticField() ;//for dca cut in ClassifyTrack(), i.e in track loop


 for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) 
{ //track loop starts for TObjArray(containing track and event information) filling; used for correlation function calculation 
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(event->GetTrack(itrk));
  if (!track) continue;
  Int_t tracktype=ClassifyTrack(track,trkVtx,bSign1,kFALSE);//don't fill the histos here
  if(tracktype!=1) continue; 

  if(!track) continue;//for safety

    gRefMultiplicityTPC += 1.0;

 }//track looop ends

 if(fCentralityMethod == "V0A_MANUAL" || fCentralityMethod == "V0M_MANUAL" || fCentralityMethod == "V0C_MANUAL" ){
  //VZERO segmentation in two detectors (0-31: VZERO-C, 32-63: VZERO-A)
  for(Int_t iChannel = 0; iChannel < 64; iChannel++) {
    fHistVZEROSignal->Fill(iChannel,event->GetVZEROEqMultiplicity(iChannel));
    
    if(iChannel < 32) 
      gRefMultiplicityVZEROC += event->GetVZEROEqMultiplicity(iChannel);
    else if(iChannel >= 32) 
      gRefMultiplicityVZEROA += event->GetVZEROEqMultiplicity(iChannel);
  }//loop over PMTs
  
  //Equalization of gain
  Double_t gFactorA = GetEqualizationFactor(gRunNumber,"A");
  if(gFactorA != 0)
    gRefMultiplicityVZEROA /= gFactorA;
  Double_t gFactorC = GetEqualizationFactor(gRunNumber,"C");
  if(gFactorC != 0)
    gRefMultiplicityVZEROC /= gFactorC;
  if((gFactorA != 0)&&(gFactorC != 0)) 
    gRefMultiplicityVZERO = (gRefMultiplicityVZEROA/gFactorA)+(gRefMultiplicityVZEROC/gFactorC);

      
  //EQVZERO vs TPC multiplicity
  fHistEQVZEROvsTPCmultiplicity->Fill(gRefMultiplicityVZERO,gRefMultiplicityTPC);
  fHistEQVZEROAvsTPCmultiplicity->Fill(gRefMultiplicityVZEROA,gRefMultiplicityTPC);
  fHistEQVZEROCvsTPCmultiplicity->Fill(gRefMultiplicityVZEROC,gRefMultiplicityTPC);

  //EQVZERO vs VZERO multiplicity
  fHistVZEROCvsEQVZEROCmultiplicity->Fill(event->GetVZEROData()->GetMTotV0C(),gRefMultiplicityVZEROC);
  fHistVZEROAvsEQVZEROAmultiplicity->Fill(event->GetVZEROData()->GetMTotV0A(),gRefMultiplicityVZEROA);

  //VZEROC vs VZEROA multiplicity
  fHistVZEROCvsVZEROAmultiplicity->Fill(event->GetVZEROData()->GetMTotV0C(),event->GetVZEROData()->GetMTotV0A());

  //EQVZEROC vs EQVZEROA multiplicity
  fHistEQVZEROCvsEQVZEROAmultiplicity->Fill(gRefMultiplicityVZEROC,gRefMultiplicityVZEROA);
 }
    fHistRefmult->Fill(3.,gRefMultiplicityTPC);
    fHistRefmult->Fill(2.,gRefMultiplicityVZERO);
    fHistRefmult->Fill(0.,gRefMultiplicityVZEROA);
    fHistRefmult->Fill(1.,gRefMultiplicityVZEROC);


 if(fCentralityMethod == "TRACKS_MANUAL") gRefMultiplicity = gRefMultiplicityTPC;
   
 else if(fCentralityMethod == "V0M_MANUAL") gRefMultiplicity = gRefMultiplicityVZERO;

 else if(fCentralityMethod == "V0A_MANUAL")  gRefMultiplicity = gRefMultiplicityVZEROA;
   
 else if(fCentralityMethod == "V0C_MANUAL") gRefMultiplicity = gRefMultiplicityVZEROC;
 
 else gRefMultiplicity = gRefMultiplicityTPC;
 
  return gRefMultiplicity;
}

//-------------------------------------------------------------------------------------------------------
Double_t AliTwoParticlePIDCorr::GetRefMultiOrCentrality(AliVEvent *mainevent, Bool_t truth){

  if(!mainevent) return -1;
  // get centrality object and check quality
  Double_t cent_v0=-1;
  Bool_t shift_to_TRACKS_MANUAL=kFALSE;//in case of wrong setting automatic shift to Tracks_Manual method

  Double_t gRefMultiplicityTPC_Truth = 0.;
  Double_t gRefMultiplicityVZERO_Truth = 0., gRefMultiplicityVZEROA_Truth = 0., gRefMultiplicityVZEROC_Truth = 0.;

  if(fAnalysisType == "AOD"|| fAnalysisType == "MCAOD") { //centrality in AOD header  //++++++++++++++
      AliAODEvent* event = dynamic_cast<AliAODEvent*>(mainevent);

if(fCentralityMethod=="V0M" || fCentralityMethod=="V0A" || fCentralityMethod=="V0C" || fCentralityMethod=="CL1" || fCentralityMethod=="ZNA" || fCentralityMethod=="V0AEq" || fCentralityMethod=="V0CEq" || fCentralityMethod=="V0MEq")//for PbPb, pPb, pp7TeV(still to be introduced)//data or RecoMC and also for TRUTH
    {
                   
 if(fSampleType=="pp_7" && fPPVsMult==kTRUE)
   {//for pp 7 TeV case only using Alianalysisutils class
     //AliAnalysisUtils *t1;
     //fAnalysisUtils=t1;

      if(!fPPVsMultUtils)
      fPPVsMultUtils=new AliPPVsMultUtils();
             
     if(fPPVsMultUtils){
       
       cent_v0 = fPPVsMultUtils->GetMultiplicityPercentile(mainevent,fCentralityMethod.Data(),fPileUp_zvtx_INEL_evsel);
       
       fHistCentStats->Fill(0.,fPPVsMultUtils->GetMultiplicityPercentile((AliVEvent*)event,"V0A",fPileUp_zvtx_INEL_evsel));
       fHistCentStats->Fill(1.,fPPVsMultUtils->GetMultiplicityPercentile((AliVEvent*)event,"V0C",fPileUp_zvtx_INEL_evsel));
       fHistCentStats->Fill(2.,fPPVsMultUtils->GetMultiplicityPercentile((AliVEvent*)event,"V0M",fPileUp_zvtx_INEL_evsel));

 // This getter(GetMultiplicityPercentile) automatically includes event selection by default and will return negative
 // values for the following types of events:
 //
 // --- Events that don't have at least one tracklet
 // --- Events without reconstructed SPD vertex
 // --- Events with a PV falling outside |z|<10cm
 // --- Events that are tagged as pileup with IsPileupFromSPDInMultBins
  /*
  fHistCentStats->Fill(3.,fAnalysisUtils->GetMultiplicityPercentile((AliVEvent*)event,"V0AEq"));//only available for LHC10d at present (Quantile info)
  fHistCentStats->Fill(4.,fAnalysisUtils->GetMultiplicityPercentile((AliVEvent*)event,"V0CEq"));//only available for LHC10d at present (Quantile info)
  fHistCentStats->Fill(5.,fAnalysisUtils->GetMultiplicityPercentile((AliVEvent*)event,"V0MEq"));//only available for LHC10d at present (Quantile info)
  */
     }
    else cent_v0 = -1;
   
      }
           
else if(fSampleType=="pPb" || fSampleType=="PbPb")
  {
  AliCentrality *centralityObj=0;
  AliAODHeader *header = (AliAODHeader*) event->GetHeader();
  if(!header) return -1;
  centralityObj = header->GetCentralityP();
  // if (centrality->GetQuality() != 0) return ;
  if(centralityObj){
  fHistCentStats->Fill(0.,centralityObj->GetCentralityPercentile("V0A"));
  fHistCentStats->Fill(1.,centralityObj->GetCentralityPercentile("V0C"));
  fHistCentStats->Fill(2.,centralityObj->GetCentralityPercentile("V0M"));
  fHistCentStats->Fill(3.,centralityObj->GetCentralityPercentile("V0AEq"));//only available for LHC10d at present (Quantile info)
  fHistCentStats->Fill(4.,centralityObj->GetCentralityPercentile("V0CEq"));//only available for LHC10d at present (Quantile info)
  fHistCentStats->Fill(5.,centralityObj->GetCentralityPercentile("V0MEq"));//only available for LHC10d at present (Quantile info)

  fHistCentStats->Fill(6.,centralityObj->GetCentralityPercentile("CL1"));
  fHistCentStats->Fill(7.,centralityObj->GetCentralityPercentile("ZNA")); 

   cent_v0 = centralityObj->GetCentralityPercentile(fCentralityMethod);
    }
  else cent_v0= -1;    
  }
 else shift_to_TRACKS_MANUAL=kTRUE;    

    }//centralitymethod condition

 else if(fCentralityMethod=="V0M_MANUAL" || fCentralityMethod=="V0A_MANUAL" || fCentralityMethod=="V0C_MANUAL" || fCentralityMethod=="TRACKS_MANUAL" || shift_to_TRACKS_MANUAL)//data or RecoMc and also for TRUTH
   {


     if(fSampleType=="pp_7") {
        

	 cent_v0=AliPPVsMultUtils::GetStandardReferenceMultiplicity(mainevent); 
       // fHistRefMult08 -> Fill ( AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ) ); 

	 if(cent_v0>0) fHistRefmult->Fill(1.,cent_v0);


     }


   else{//for pPb and PbPb only
     if(!truth){//for data or RecoMC
       cent_v0 = GetReferenceMultiplicityVZEROFromAOD((AliVEvent*)event);
   }//for data or RecoMC

    if(truth && (fAnalysisType == "MCAOD")){//condition for TRUTH case
//check for TClonesArray(truth track MC information)
fArrayMC = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fArrayMC) {
    //AliFatal("Error: MC particles branch not found!\n");
    return -1;
  }
//now process the truth particles(for both efficiency & correlation function)
Int_t nMCTrack = fArrayMC->GetEntriesFast();
  
for (Int_t iMC = 0; iMC < nMCTrack; iMC++) 
{//MC truth track loop starts
    
AliAODMCParticle *partMC = (AliAODMCParticle*) fArrayMC->At(iMC);
    
    if(!partMC){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",iMC));
      continue;
    }

//consider only charged particles
    if(partMC->Charge() == 0) continue;

//consider only primary particles; neglect all secondary particles including from weak decays
 if(fselectprimaryTruth && !partMC->IsPhysicalPrimary()) continue;


//remove injected signals(primaries above <maxLabel>)
 if (fInjectedSignals && partMC->GetLabel() >= skipParticlesAbove) continue;

//remove duplicates
  Bool_t isduplicate=kFALSE;
 if (fRemoveDuplicates)
   { 
 for (Int_t j=iMC+1; j<nMCTrack; ++j) 
   {//2nd trutuh loop starts
AliAODMCParticle *partMC2 = (AliAODMCParticle*) fArrayMC->At(j);
   if(!partMC2){
      AliError(Form("ERROR: Could not retrieve AODMCtrack %d",j));
      continue;
    }    
 if (partMC->GetLabel() == partMC2->GetLabel())
   {
isduplicate=kTRUE;
 break;  
   }    
   }//2nd truth loop ends
   }
 if(fRemoveDuplicates && isduplicate) continue;//remove duplicates


       //  if (fCentralityMethod=="V0M_MANUAL") 
	if((partMC->Eta() < 5.1 && partMC->Eta() > 2.8) || (partMC->Eta() > -3.7 && partMC->Eta() < -1.7)) gRefMultiplicityVZERO_Truth+=1;
	//   else if (fCentralityMethod=="V0A_MANUAL") {
	if(partMC->Eta() < 5.1 && partMC->Eta() > 2.8)  gRefMultiplicityVZEROA_Truth+=1;
	// else if (fCentralityMethod=="V0C_MANUAL") {
	if(partMC->Eta() < -1.7 && partMC->Eta() > -3.7)  gRefMultiplicityVZEROC_Truth+=1;
	//else if (fCentralityMethod=="TRACKS_MANUAL") {
        if (partMC->Eta() > fmineta && partMC->Eta() < fmaxeta) {
	  if (partMC->Pt() > fminPt &&  partMC->Pt() < fmaxPt) gRefMultiplicityTPC_Truth+=1;
           }
     
 }//truth track loop ends

 fHistRefmult->Fill(3.,gRefMultiplicityTPC_Truth);
 fHistRefmult->Fill(2.,gRefMultiplicityVZERO_Truth); 
 fHistRefmult->Fill(0.,gRefMultiplicityVZEROA_Truth);
 fHistRefmult->Fill(1.,gRefMultiplicityVZEROC_Truth);

 if(fCentralityMethod == "TRACKS_MANUAL")   cent_v0=gRefMultiplicityTPC_Truth;

 else if(fCentralityMethod == "V0M_MANUAL")  cent_v0=gRefMultiplicityVZERO_Truth;

 else if(fCentralityMethod == "V0A_MANUAL")  cent_v0=gRefMultiplicityVZEROA_Truth;

 else if(fCentralityMethod == "V0C_MANUAL")  cent_v0=gRefMultiplicityVZEROC_Truth;
 
 else      cent_v0=gRefMultiplicityTPC_Truth;

    }//condition for TRUTH case

   }//condition for else(sampletype other than pp_7 ends)

   }//end of MANUAL method

 else if ((fAnalysisType == "MCAOD") && (fCentralityMethod == "MC_b"))//TRUTH MC in AOD production(Impaact parm is not used in data or RecoMC case)
    {
    AliAODMCHeader* header = (AliAODMCHeader*) event->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!header)
    return -1;
    
      AliGenEventHeader* eventHeader = header->GetCocktailHeader(0);  // get first MC header from either ESD/AOD (including cocktail header if available)
      if (!eventHeader)
      {
	// We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
	// (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
	AliError("Event header not found. Skipping this event.");
	return -1;
      }
      
      AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*> (eventHeader);
     
      
      if (collGeometry) {
  cent_v0 = collGeometry->ImpactParameter();
  fhistImpactParm->Fill(cent_v0);
      }
      else cent_v0=-1.;
    }//end of Impact parameter method

//else return -1
 else cent_v0=-1.;
}//AOD OR MCAOD condition


else if(fAnalysisType == "MC"){
    Double_t gImpactParameter = -1.;
    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(mainevent);
    if(gMCEvent){
      AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());      
      if(headerH){
	gImpactParameter = headerH->ImpactParameter();

 for(Int_t iParticle = 0; iParticle < gMCEvent->GetNumberOfPrimaries(); iParticle++) {
      AliMCParticle* track = dynamic_cast<AliMCParticle *>(gMCEvent->GetTrack(iParticle));
      if (!track) {
	AliError(Form("Could not receive particle %d", iParticle));
	continue;
      }
      
      //exclude non stable particles
      if(fselectprimaryTruth && !(gMCEvent->IsPhysicalPrimary(iParticle))) continue;

      if(track->Charge() == 0) continue;

 //  if (fCentralityMethod=="V0M_MANUAL") 
	if((track->Eta() < 5.1 && track->Eta() > 2.8) || (track->Eta() > -3.7 && track->Eta() < -1.7)) gRefMultiplicityVZERO_Truth+=1;
	//   else if (fCentralityMethod=="V0A_MANUAL") {
	if(track->Eta() < 5.1 && track->Eta() > 2.8)  gRefMultiplicityVZEROA_Truth+=1;
	// else if (fCentralityMethod=="V0C_MANUAL") {
	if(track->Eta() < -1.7 && track->Eta() > -3.7)  gRefMultiplicityVZEROC_Truth+=1;
	//else if (fCentralityMethod=="TRACKS_MANUAL") {
        if (track->Eta() > fmineta && track->Eta() < fmaxeta) {
	  if (track->Pt() > fminPt &&  track->Pt() < fmaxPt) gRefMultiplicityTPC_Truth+=1;}

     }//loop over primaries

 fHistRefmult->Fill(3.,gRefMultiplicityTPC_Truth);
 fHistRefmult->Fill(2.,gRefMultiplicityVZERO_Truth); 
 fHistRefmult->Fill(0.,gRefMultiplicityVZEROA_Truth);
 fHistRefmult->Fill(1.,gRefMultiplicityVZEROC_Truth);
 if (fCentralityMethod == "MC_b"){
       cent_v0=gImpactParameter;
       fhistImpactParm->Fill(gImpactParameter);
       fhistImpactParmvsMult->Fill(gImpactParameter,gRefMultiplicityTPC_Truth);
 }

 else if(fCentralityMethod == "TRACKS_MANUAL")   cent_v0=gRefMultiplicityTPC_Truth;

 else if(fCentralityMethod == "V0M_MANUAL")  cent_v0=gRefMultiplicityVZERO_Truth;

 else if(fCentralityMethod == "V0A_MANUAL")  cent_v0=gRefMultiplicityVZEROA_Truth;

 else if(fCentralityMethod == "V0C_MANUAL")  cent_v0=gRefMultiplicityVZEROC_Truth;
 
 else      cent_v0=gImpactParameter;//default value is the impact parameter
    }//MC event header
    }//MC event cast
    else   cent_v0 = -1.;
  }//MC condition

  else{
    cent_v0 = -1.;
  }
 return cent_v0;
}
//-----------------------------------------------------------------------------------------
Double_t AliTwoParticlePIDCorr::GetAcceptedEventMultiplicity(AliVEvent *event,Bool_t truth){
  //do the event selection(zvtx, pileup, centrality/multiplicity cut) and then return the value of the centrality of that event
  if(!event) return -1;

  Float_t gRefMultiplicity = -1.;

   //***********************************SOURCE CODE-TASKBFPsi

   // Event Vertex MC
    if(fAnalysisType == "MC") {
    AliMCEvent *mcevent = dynamic_cast<AliMCEvent*>(event);
      if(!mcevent) {
       AliError("mcEvent not available");
      return -1.;
      }
      
      if(mcevent){
	AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(mcevent->GenEventHeader());
	if(header){  
	  TArrayF gVertexArray;
	  header->PrimaryVertex(gVertexArray);
 //count events having a proper vertex
          fEventCounter->Fill(5);	  

fHistQA[0]->Fill((gVertexArray.At(0)));fHistQA[1]->Fill((gVertexArray.At(1)));fHistQA[2]->Fill((gVertexArray.At(2))); //for trkVtx only before vertex cut |zvtx|<10 cm

       if(TMath::Abs(gVertexArray.At(0)) < fVxMax_MC) {
	    if(TMath::Abs(gVertexArray.At(1)) < fVyMax_MC) {
	      if(TMath::Abs(gVertexArray.At(2)) < fVzMax_MC) {
//count events after vertex cut
                 fEventCounter->Fill(7);
 fHistQA[3]->Fill((gVertexArray.At(0)));fHistQA[4]->Fill((gVertexArray.At(1)));fHistQA[5]->Fill((gVertexArray.At(2)));//after vertex cut,for trkVtx only

		// get the reference multiplicty or centrality
 gRefMultiplicity = GetRefMultiOrCentrality((AliVEvent*)mcevent,kFALSE);//2nd argument has no meaning

               if(gRefMultiplicity<0) return -1;

 // take events only within the  multiplicity class mentioned in the custom binning
  if(gRefMultiplicity < fmincentmult || gRefMultiplicity > fmaxcentmult) return -1;

//count events having proper centrality/ref multiplicity
                    fEventCounter->Fill(9);

	      }//Vz cut
	    }//Vy cut
	  }//Vx cut	   
      }//MC event header
    }//MC event object
   }//MC

    else  if(fAnalysisType == "MCAOD" || fAnalysisType == "AOD"){// if(fAnalysisType == "MCAOD" || fAnalysisType == "AOD"
  //vertex selection(is it fine for PP?)
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  if(fSampleType=="pp_7"){
    
    if(!fPPVsMultUtils) fPPVsMultUtils=new AliPPVsMultUtils();
    Bool_t eventselected=kFALSE;
    eventselected = fPPVsMultUtils->IsEventSelected(aod); 

      if(eventselected){

    gRefMultiplicity = GetRefMultiOrCentrality((AliVEvent*)aod,kFALSE);
    //count events having proper centrality/ref multiplicity
    if(gRefMultiplicity>0) fEventCounter->Fill(9);
      }
      
      else gRefMultiplicity=-1;
    }


  else{//for pPb and PbPb
 // check first event in chunk (is not needed for new reconstructions)
  if(fCheckFirstEventInChunk){
    AliAnalysisUtils ut;
    if(ut.IsFirstEventInChunk(aod)) 
      return -1.;
  }

 if(frejectPileUp){
    AliAnalysisUtils ut;
    ut.SetUseMVPlpSelection(kTRUE);
    ut.SetUseOutOfBunchPileUp(kTRUE);
    if(ut.IsPileUpEvent(aod))
      return -1.;
  }

//count events after pileup selection
   fEventCounter->Fill(3);

 if (fVertextype==1){//for pPb basically if(!fAnalysisUtils->IsVertexSelected2013pA(aod)) return; 
   trkVtx = aod->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return -1;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return -1;
   zvtx = trkVtx->GetZ();
  const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
  if (!spdVtx || spdVtx->GetNContributors()<=0) return -1;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return -1;
   if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return -1;
  }
  else if(fVertextype==2) {//for pp and pb-pb case ,used in AliAnalysisTaskPhiCorelations.cxx
	Int_t nVertex = aod->GetNumberOfVertices();
  	if( nVertex > 0 ) { 
     trkVtx = (AliAODVertex*)aod->GetPrimaryVertex();
		Int_t nTracksPrim = trkVtx->GetNContributors();
                 zvtx = trkVtx->GetZ();
  		//if (fDebug > 1)AliInfo(Form(" Vertex in = %f with %d particles by  %s data ...",zVertex,nTracksPrim,vertex->GetName()));
  		// Reject TPC only vertex
		TString name(trkVtx->GetName());
		if (name.CompareTo("PrimaryVertex") && name.CompareTo("SPDVertex"))return -1;

		// Select a quality vertex by number of tracks?
  		if( nTracksPrim < fnTracksVertex ) {
		  //if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
  			return -1;
  			}
  		// TODO remove vertexer Z events with dispersion > 0.02: Doesn't work for AOD at present
                //if (strcmp(vertex->GetTitle(), "AliVertexerZ") == 0 && vertex->GetDispersion() > 0.02)
                //  return kFALSE;
		//	if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED...");
	}
	else return -1;

  }
 else if(fVertextype==0){//default case
  trkVtx =(AliAODVertex*) aod->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return -1;//proper number of contributors
  zvtx = trkVtx->GetZ();
  Double32_t fCov[6];
  trkVtx->GetCovarianceMatrix(fCov);
  if(fCov[5] == 0) return -1;//proper vertex resolution
  }
  else {
   AliInfo("Wrong Vertextype set for Primary-vertex Selection: event REJECTED ...");
   return -1;//as there is no proper sample type
  }

fHistQA[0]->Fill((trkVtx->GetX()));fHistQA[1]->Fill((trkVtx->GetY()));fHistQA[2]->Fill((trkVtx->GetZ())); //for trkVtx only before vertex cut |zvtx|<10 cm

//count events having a proper vertex
   fEventCounter->Fill(5);

 if (TMath::Abs(zvtx) > fzvtxcut) return -1;

//count events after vertex cut
  fEventCounter->Fill(7);


 //if(!fAnalysisUtils->IsVertexSelected2013pA(aod)) return;
  
 fHistQA[3]->Fill((trkVtx->GetX()));fHistQA[4]->Fill((trkVtx->GetY()));fHistQA[5]->Fill((trkVtx->GetZ()));//after vertex cut,for trkVtx only

 //get the centrality or multiplicity
 if(truth)  {gRefMultiplicity = GetRefMultiOrCentrality((AliVEvent*)aod,kTRUE);}//kTRUE-->for Truth case(only meaningful in case of ref multiplicity,in case of centrality it has no meaning)

 else {gRefMultiplicity = GetRefMultiOrCentrality((AliVEvent*)aod,kFALSE);}//kFALSE-->for data and RecoMc case(only meaningful in case of ref multiplicity,in case of centrality it has no meaning)

  if(gRefMultiplicity<0) return -1;

 // take events only within the  multiplicity class mentioned in the custom binning
  if(gRefMultiplicity < fmincentmult || gRefMultiplicity > fmaxcentmult) return -1;

//count events having proper centrality/ref multiplicity
  fEventCounter->Fill(9);


// centrality weighting (optional for 2011 if central and semicentral triggers are used);only for data and recoMC
 if (fCentralityWeights && !AcceptEventCentralityWeight(gRefMultiplicity))//**********************
  {
    AliInfo(Form("Rejecting event because of centrality weighting: %f", gRefMultiplicity));
    return -1;
  }

//count events after rejection due to centrality weighting
  fEventCounter->Fill(11);

  }//pPb and PbPb condition ends
    }//AOD or MCAOD ends
    else gRefMultiplicity=-1;

  return gRefMultiplicity;

}
//--------------------------------------------------------------------------------------------------------
Float_t AliTwoParticlePIDCorr::GetEventPlane(AliVEvent *mainevent,Bool_t truth, Double_t v0Centr)
{
  Float_t eventplane=999.;
  // Get the event plane
  if(!mainevent) return 999.;


 //MC: from reaction plane
  if(fAnalysisType == "MC"){
    if(!mainevent) {
      AliError("mcEvent not available");
      return 999.;
    }

    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(mainevent);
    if(gMCEvent){
      AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());    
      if (headerH) {
 Int_t iC = -1;  
    // Impact parameter bins(it is only for Pb-Pb)
    if(v0Centr < 3.50) iC = 0;
    else if(v0Centr < 4.94) iC = 1;
    else if(v0Centr < 6.98) iC = 2;
    else if(v0Centr < 8.55) iC = 3;
    else if(v0Centr < 9.88) iC = 4;
    else if(v0Centr < 11.04) iC = 5;
    else if(v0Centr < 12.09) iC = 6;
    else if(v0Centr < 13.05) iC = 7;
    else iC = 8;

	eventplane = headerH->ReactionPlaneAngle();
	if(eventplane > TMath::Pi()/2 && eventplane <=  TMath::Pi()*3/2) eventplane-=TMath::Pi(); 
	 if(eventplane > TMath::Pi()*3/2) eventplane-=2*TMath::Pi();
         fHistEventPlaneTruth->Fill(iC,eventplane);
	//gReactionPlane *= TMath::RadToDeg();
      }//MC header
    }//MC event cast
  }//MC
  
  else  if(fAnalysisType == "MCAOD" || fAnalysisType == "AOD") {
 //reset Q vector info	

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(mainevent);


    Int_t run = event->GetRunNumber();

    if(run != fRun){
	// Load the calibrations run dependent
      if(! fIsAfter2011) OpenInfoCalbration(run);
      fRun=run;
    }


  Int_t iC = -1;  
  if (v0Centr > 80) return 999.; // analysis only for 0-80% centrality classes
 // centrality bins
    if(v0Centr < 5) iC = 0;
    else if(v0Centr < 10) iC = 1;
    else if(v0Centr < 20) iC = 2;
    else if(v0Centr < 30) iC = 3;
    else if(v0Centr < 40) iC = 4;
    else if(v0Centr < 50) iC = 5;
    else if(v0Centr < 60) iC = 6;
    else if(v0Centr < 70) iC = 7;
    else iC = 8;


     Int_t iCcal = iC;

 //reset Q vector info	 
    Double_t Qxa2 = 0, Qya2 = 0;
    Double_t Qxc2 = 0, Qyc2 = 0;
    Double_t Qxa3 = 0, Qya3 = 0;
    Double_t Qxc3 = 0, Qyc3 = 0;


  //MC: from reaction plane
 if(truth)
{
    AliAODMCHeader* header = (AliAODMCHeader*) event->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (header){
      evplaneMC = header->GetReactionPlaneAngle();//[0, 360]
        //make it within [-pi/2,pi/2] to make it general
	if(evplaneMC > TMath::Pi()/2 && evplaneMC <=  TMath::Pi()*3/2) evplaneMC-=TMath::Pi(); 
	 if(evplaneMC > TMath::Pi()*3/2) evplaneMC-=2*TMath::Pi();
         fHistEventPlaneTruth->Fill(iC,evplaneMC);
	/*
        AliGenEventHeader* eventHeader = header->GetCocktailHeader(0);  // get first MC header from either ESD/AOD (including cocktail header if available)
      if (eventHeader)
      {
	      
	AliCollisionGeometry* collGeometry = dynamic_cast<AliCollisionGeometry*> (eventHeader);     
      
	if (collGeometry){//get the reaction plane from MC header   
	  gReactionPlane = collGeometry->ReactionPlaneAngle();//[0,180]
 }
      }
	*/   
     //taken from vnv0 code(get the TPC, V0A, V0C event plane using truth tracks)
       TClonesArray *mcArray = NULL;
	mcArray = (TClonesArray*)event->GetList()->FindObject(AliAODMCParticle::StdBranchName());
	if(mcArray){
	  Float_t QxMCv2[3] = {0,0,0};
	  Float_t QyMCv2[3] = {0,0,0};
	  Float_t QxMCv3[3] = {0,0,0};
	  Float_t QyMCv3[3] = {0,0,0};
	  Float_t EvPlaneMCV2[3] = {0,0,0};
	  Float_t EvPlaneMCV3[3] = {0,0,0};
	  Float_t etaMin[3] = {2.8,-3.6,-0.8}; // A-side, C-side M-barrel
	  Float_t etaMax[3] = {4.88,-1.8,0.8};

	  // analysis on MC tracks
	  Int_t nMCtrack = mcArray->GetEntries() ;

	  // EP computation with MC tracks
	  for(Int_t iT=0;iT < nMCtrack;iT++){
	    AliAODMCParticle *mctr = (AliAODMCParticle*) mcArray->At(iT);
	    if(!mctr || !(mctr->IsPrimary()) || !(mctr->Charge()) || mctr->Pt() < 0.2) continue;
	    
	    Float_t eta = mctr->Eta();
  for(Int_t iD=0;iD<3;iD++){
	      if(eta > etaMin[iD] && eta < etaMax[iD]){
		Float_t phi = mctr->Phi();
		QxMCv2[iD] += TMath::Cos(2*phi);
		QyMCv2[iD] += TMath::Sin(2*phi);
		QxMCv3[iD] += TMath::Cos(3*phi);
		QyMCv3[iD] += TMath::Sin(3*phi);
	      }
	    }
	  }

	    EvPlaneMCV2[0] = TMath::ATan2(QyMCv2[0],QxMCv2[0])/2.;
	    EvPlaneMCV2[1] = TMath::ATan2(QyMCv2[1],QxMCv2[1])/2.;
	    EvPlaneMCV2[2] = TMath::ATan2(QyMCv2[2],QxMCv2[2])/2.;
	    fHResMA2->Fill(Double_t(iC), TMath::Cos(2*(EvPlaneMCV2[2]-EvPlaneMCV2[0])));
	    fHResMC2->Fill(Double_t(iC), TMath::Cos(2*(EvPlaneMCV2[2]-EvPlaneMCV2[1])));
	    fHResAC2->Fill(Double_t(iC), TMath::Cos(2*(EvPlaneMCV2[0]-EvPlaneMCV2[1])));
            fgPsi2v0aMC = EvPlaneMCV2[0];
            fgPsi2v0cMC = EvPlaneMCV2[1];
            fgPsi2tpcMC = EvPlaneMCV2[2];
	  

	    EvPlaneMCV3[0] = TMath::ATan2(QyMCv3[0],QxMCv3[0])/3.;
	    EvPlaneMCV3[1] = TMath::ATan2(QyMCv3[1],QxMCv3[1])/3.;
	    EvPlaneMCV3[2] = TMath::ATan2(QyMCv3[2],QxMCv3[2])/3.;
	    fHResMA3->Fill(Double_t(iC), TMath::Cos(3*(EvPlaneMCV3[2]-EvPlaneMCV3[0])));
	    fHResMC3->Fill(Double_t(iC), TMath::Cos(3*(EvPlaneMCV3[2]-EvPlaneMCV3[1])));
	    fHResAC3->Fill(Double_t(iC), TMath::Cos(3*(EvPlaneMCV3[0]-EvPlaneMCV3[1])));
            fgPsi3v0aMC = EvPlaneMCV3[0];
            fgPsi3v0cMC = EvPlaneMCV3[1];
            fgPsi3tpcMC = EvPlaneMCV3[2];
	  
	}    
    }

    }
 else{
    Int_t nAODTracks = event->GetNumberOfTracks();

// TPC EP needed for resolution studies (TPC subevent)
   //AliEventplane * ep = (fAOD->GetHeader())->GetEventplaneP();
   //Double_t psiTPC = ep->GetEventplane("Q", fAOD, 2); // in range of [0, pi]
    Double_t Qx2 = 0, Qy2 = 0;
    Double_t Qx3 = 0, Qy3 = 0;

    for(Int_t iT = 0; iT < nAODTracks; iT++) {
      
      AliAODTrack* aodTrack =(AliAODTrack*) event->GetTrack(iT);
      
      if (!aodTrack){
	continue;
      }
      
      Bool_t trkFlag = aodTrack->TestFilterBit(1);

      if ((TMath::Abs(aodTrack->Eta()) > 0.8) || (aodTrack->Pt() < 0.2) || (aodTrack->GetTPCNcls() < fNcluster)  || !trkFlag) 
	continue;
	
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};


      AliAODTrack param(*aodTrack);
      if (!param.PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov)){
	continue;
      }
	    
      if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4))
	continue;
      
      Qx2 += TMath::Cos(2*aodTrack->Phi()); 
      Qy2 += TMath::Sin(2*aodTrack->Phi());
      Qx3 += TMath::Cos(3*aodTrack->Phi()); 
      Qy3 += TMath::Sin(3*aodTrack->Phi());
      
    }
    
   Float_t evPlAng2 = TMath::ATan2(Qy2, Qx2)/2.;
   Float_t evPlAng3 = TMath::ATan2(Qy3, Qx3)/3.;

    fgPsi2tpc = evPlAng2;
    fgPsi3tpc = evPlAng3;

     fPhiRPTPC->Fill(iC,evPlAng2);
     fPhiRPTPCv3->Fill(iC,evPlAng3);



//V0 info    
    AliAODVZERO* aodV0 = event->GetVZEROData();

    for (Int_t iv0 = 0; iv0 < 64; iv0++) {
      Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
      Float_t multv0 = aodV0->GetMultiplicity(iv0);

      if(! fIsAfter2011){
	if(fAnalysisType == "AOD"){//not for reco MC tracks, only for real data
	  if (iv0 < 32){ // V0C
	    Qxc2 += TMath::Cos(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qyc2 += TMath::Sin(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qxc3 += TMath::Cos(3*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qyc3 += TMath::Sin(3*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	  } else {       // V0A
	    Qxa2 += TMath::Cos(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qya2 += TMath::Sin(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qxa3 += TMath::Cos(3*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qya3 += TMath::Sin(3*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	  }
	}
	else{
	  if (iv0 < 32){ // V0C
	    Qxc2 += TMath::Cos(2*phiV0) * multv0;//*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qyc2 += TMath::Sin(2*phiV0) * multv0;//*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qxc3 += TMath::Cos(3*phiV0) * multv0;//*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	    Qyc3 += TMath::Sin(3*phiV0) * multv0;//*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	  } else {       // V0A
	    Qxa2 += TMath::Cos(2*phiV0) * multv0;//*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qya2 += TMath::Sin(2*phiV0) * multv0;//*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qxa3 += TMath::Cos(3*phiV0) * multv0;//*fV0Apol/fMultV0->GetBinContent(iv0+1);
	    Qya3 += TMath::Sin(3*phiV0) * multv0;//*fV0Apol/fMultV0->GetBinContent(iv0+1);
	  }
	}
      }
    }
   //grab for each centrality the proper histo with the Qx and Qy to do the recentering
    Double_t Qxamean2 = fMeanQ[iCcal][1][0];
    Double_t Qxarms2  = fWidthQ[iCcal][1][0];
    Double_t Qyamean2 = fMeanQ[iCcal][1][1];
    Double_t Qyarms2  = fWidthQ[iCcal][1][1];
    Double_t Qxamean3 = fMeanQv3[iCcal][1][0];
    Double_t Qxarms3  = fWidthQv3[iCcal][1][0];
    Double_t Qyamean3 = fMeanQv3[iCcal][1][1];
    Double_t Qyarms3  = fWidthQv3[iCcal][1][1];
    
    Double_t Qxcmean2 = fMeanQ[iCcal][0][0];
    Double_t Qxcrms2  = fWidthQ[iCcal][0][0];
    Double_t Qycmean2 = fMeanQ[iCcal][0][1];
    Double_t Qycrms2  = fWidthQ[iCcal][0][1];	
    Double_t Qxcmean3 = fMeanQv3[iCcal][0][0];
    Double_t Qxcrms3  = fWidthQv3[iCcal][0][0];
    Double_t Qycmean3 = fMeanQv3[iCcal][0][1];
    Double_t Qycrms3  = fWidthQv3[iCcal][0][1];	
    
    Double_t QxaCor2 = (Qxa2 - Qxamean2)/Qxarms2;
    Double_t QyaCor2 = (Qya2 - Qyamean2)/Qyarms2;
    Double_t QxcCor2 = (Qxc2 - Qxcmean2)/Qxcrms2;
    Double_t QycCor2 = (Qyc2 - Qycmean2)/Qycrms2;
    Double_t QxaCor3 = (Qxa3 - Qxamean3)/Qxarms3;
    Double_t QyaCor3 = (Qya3 - Qyamean3)/Qyarms3;
    Double_t QxcCor3 = (Qxc3 - Qxcmean3)/Qxcrms3;
    Double_t QycCor3 = (Qyc3 - Qycmean3)/Qycrms3;
    /*
    //to calculate 2nd order event plane with v0M
 Double_t QxCor2 = (Qxa2 - Qxamean2 + Qxc2 - Qxcmean2)
    /TMath::Sqrt(Qxarms2*Qxarms2 + Qxcrms2*Qxcrms2);
  Double_t QyCor2 = (Qya2 - Qyamean2 + Qyc2 - Qycmean2)
    /TMath::Sqrt(Qyarms2*Qyarms2 + Qycrms2*Qycrms2);

  //here the calculated event plane is within -Pi to +Pi(delete it , no use here , only for definition)
  Double_t psiV0A =(TMath::Pi() + TMath::ATan2(-QyaCor2, -QxaCor2))/2.;
  Double_t psiV0C = (TMath::Pi() + TMath::ATan2(-QycCor2, -QxcCor2))/2.;
  Double_t psiVZero = (TMath::Pi() + TMath::ATan2(-QyCor2, -QxCor2))/2.;

    */

    Float_t evPlAngV0ACor2=999.;
    Float_t evPlAngV0CCor2=999.;
    Float_t evPlAngV0ACor3=999.;
    Float_t evPlAngV0CCor3=999.;

   if(! fIsAfter2011){
      if(fAnalysisType == "AOD"){
	evPlAngV0ACor2 = TMath::ATan2(QyaCor2, QxaCor2)/2.;
	evPlAngV0CCor2 = TMath::ATan2(QycCor2, QxcCor2)/2.;
	evPlAngV0ACor3 = TMath::ATan2(QyaCor3, QxaCor3)/3.;
	evPlAngV0CCor3 = TMath::ATan2(QycCor3, QxcCor3)/3.;
      }
      else{
	evPlAngV0ACor2 = TMath::ATan2(Qya2, Qxa2)/2.;
	evPlAngV0CCor2 = TMath::ATan2(Qyc2, Qxc2)/2.;
	evPlAngV0ACor3 = TMath::ATan2(Qya3, Qxa3)/3.;
	evPlAngV0CCor3 = TMath::ATan2(Qyc3, Qxc3)/3.;
      }
    }
    else{
      AliEventplane *ep =  event->GetEventplane();
      evPlAngV0ACor2 = ep->GetEventplane("V0A", event, 2);
      evPlAngV0CCor2 = ep->GetEventplane("V0C", event, 2);
      evPlAngV0ACor3 = ep->GetEventplane("V0A", event, 3);
      evPlAngV0CCor3 = ep->GetEventplane("V0C", event, 3);
    }

    fgPsi2v0a = evPlAngV0ACor2;
    fgPsi2v0c = evPlAngV0CCor2;
    fgPsi3v0a = evPlAngV0ACor3;
    fgPsi3v0c = evPlAngV0CCor3;

 // Fill EP distribution histograms evPlAng
    
     fPhiRPv0A->Fill(iC,evPlAngV0ACor2);
     fPhiRPv0C->Fill(iC,evPlAngV0CCor2);
    
     fPhiRPv0Av3->Fill(iC,evPlAngV0ACor3);
     fPhiRPv0Cv3->Fill(iC,evPlAngV0CCor3);

    // Fill histograms needed for resolution evaluation
    fHResTPCv0A2->Fill(Double_t(iC), TMath::Cos(2*(evPlAng2 - evPlAngV0ACor2)));
    fHResTPCv0C2->Fill(Double_t(iC), TMath::Cos(2*(evPlAng2 - evPlAngV0CCor2)));
    fHResv0Cv0A2->Fill(Double_t(iC), TMath::Cos(2*(evPlAngV0ACor2 - evPlAngV0CCor2)));
    
    fHResTPCv0A3->Fill(Double_t(iC), TMath::Cos(3*(evPlAng3 - evPlAngV0ACor3)));
    fHResTPCv0C3->Fill(Double_t(iC), TMath::Cos(3*(evPlAng3 - evPlAngV0CCor3)));
    fHResv0Cv0A3->Fill(Double_t(iC), TMath::Cos(3*(evPlAngV0ACor3 - evPlAngV0CCor3)));


    /*   
 Float_t gVZEROEventPlane    = -10.;
  Float_t gReactionPlane      = -10.;
  Double_t qxTot = 0.0, qyTot = 0.0;

    AliEventplane *ep = event->GetEventplane();
    if(ep){ 
      gVZEROEventPlane = ep->CalculateVZEROEventPlane(event,10,2,qxTot,qyTot);
      if(gVZEROEventPlane < 0.) gVZEROEventPlane += TMath::Pi();
      //gReactionPlane = gVZEROEventPlane*TMath::RadToDeg();
      gReactionPlane = gVZEROEventPlane;
    }
    */
  }//AOD,ESD,ESDMC
 //return gReactionPlane;

 //make the final 2nd order event plane within 0 to Pi
     //using data and reco tracks only
      if(fgPsi2v0a!=999. && fgPsi2v0a < 0.) fgPsi2v0a += TMath::Pi();
      if(fgPsi2v0c!=999. && fgPsi2v0c < 0.) fgPsi2v0c += TMath::Pi();
      if(fgPsi2tpc!=999. && fgPsi2tpc < 0.) fgPsi2tpc += TMath::Pi();
      //using truth tracks only
      if(evplaneMC!=999. && evplaneMC < 0.) evplaneMC += TMath::Pi();
      if(fgPsi2v0aMC!=999. && fgPsi2v0aMC < 0.) fgPsi2v0aMC += TMath::Pi();
      if(fgPsi2v0cMC!=999. && fgPsi2v0cMC < 0.) fgPsi2v0cMC += TMath::Pi();
      if(fgPsi2tpcMC!=999. && fgPsi2tpcMC < 0.) fgPsi2tpcMC += TMath::Pi();
      //for the time being leave the 3rd order event planes within -pi/3 t0 +pi/3

      if(truth){//for truth MC
	if(fV2 && fEPdet=="header")eventplane=evplaneMC;
	if(fV2 && fEPdet=="V0A")eventplane=fgPsi2v0aMC;
	if(fV2 && fEPdet=="V0C")eventplane=fgPsi2v0cMC;
	if(fV2 && fEPdet=="TPC")eventplane=fgPsi2tpcMC;

	if(fV3 && fEPdet=="V0A")eventplane=fgPsi3v0aMC;
	if(fV3 && fEPdet=="V0C")eventplane=fgPsi3v0cMC;
	if(fV3 && fEPdet=="TPC")eventplane=fgPsi3tpcMC;
}
      else{//for data and recoMC
	if(fV2 && fEPdet=="V0A")eventplane=fgPsi2v0a;
	if(fV2 && fEPdet=="V0C")eventplane=fgPsi2v0c;
	if(fV2 && fEPdet=="TPC")eventplane=fgPsi2tpc;

	if(fV3 && fEPdet=="V0A")eventplane=fgPsi3v0a;
	if(fV3 && fEPdet=="V0C")eventplane=fgPsi3v0c;
	if(fV3 && fEPdet=="TPC")eventplane=fgPsi3tpc;

      }
     

  }//AOD/MCAOD condition

  else eventplane=999.;
 
  return eventplane;

}
//------------------------------------------------------------------------------------------------------------------
void AliTwoParticlePIDCorr::OpenInfoCalbration(Int_t run){
    TString oadbfilename = "$ALICE_ROOT/OADB/PWGCF/VZERO/VZEROcalibEP.root";
    TFile *foadb = TFile::Open(oadbfilename.Data());

    if(!foadb){
	printf("OADB file %s cannot be opened\n",oadbfilename.Data());
	return;
    }

    AliOADBContainer *cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorr");
    if(!cont){
	printf("OADB object hMultV0BefCorr is not available in the file\n");
	return;	
    }

    if(!(cont->GetObject(run))){
	printf("OADB object hMultV0BefCorr is not available for run %i (used run 137366)\n",run);
	run = 137366;
    }
    fMultV0 = ((TH2F *) cont->GetObject(run))->ProfileX();

    TF1 *fpol0 = new TF1("fpol0","pol0"); 
    fMultV0->Fit(fpol0,"","",0,31);
    fV0Cpol = fpol0->GetParameter(0);
    fMultV0->Fit(fpol0,"","",32,64);
    fV0Apol = fpol0->GetParameter(0);

    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < 9;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya2_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
		fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

		//for v3
		if(iside==0 && icoord==0)
		  snprintf(namecont,100,"hQxc3_%i",i);
		else if(iside==1 && icoord==0)
		  snprintf(namecont,100,"hQxa3_%i",i);
		else if(iside==0 && icoord==1)
		  snprintf(namecont,100,"hQyc3_%i",i);
		else if(iside==1 && icoord==1)
		  snprintf(namecont,100,"hQya3_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    printf("OADB object %s is not available in the file\n",namecont);
		    return;	
		}
		
		if(!(cont->GetObject(run))){
		    printf("OADB object %s is not available for run %i (used run 137366)\n",namecont,run);
		    run = 137366;
		}
		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

     	    }
	}
    }
}
//____________________________________________________________________
void AliTwoParticlePIDCorr::FillPIDEventPlane(Double_t centrality,Int_t par,Float_t trigphi,Float_t fReactionPlane) 
{

 // Event plane (determine psi bin)
    Double_t gPsiMinusPhi    =   0.;
    Double_t gPsiMinusPhiBin = -10.;
if(fRequestEventPlane){
    gPsiMinusPhi   = TMath::Abs(trigphi - fReactionPlane);
    //in-plane
    if((gPsiMinusPhi <= 7.5*TMath::DegToRad())||
      (gPsiMinusPhi >= 352.5*TMath::DegToRad())||
       ((172.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 0.0;
    //intermediate
    else if(((37.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5*TMath::DegToRad()))||
	    ((127.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5*TMath::DegToRad()))||
	    ((217.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5*TMath::DegToRad()))||
	    ((307.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 1.0;
    //out of plane
    else if(((82.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5*TMath::DegToRad()))||
	    ((262.5*TMath::DegToRad() <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5*TMath::DegToRad())))
      gPsiMinusPhiBin = 2.0;
    //everything else
    else 
      gPsiMinusPhiBin = 3.0;

    fEventPlanePID->Fill(centrality,gPsiMinusPhiBin,(Float_t)par); 
 }
}




//___________________________________________________________________________________________

/* Int_t AliTwoParticlePIDCorr::GetCentBin(Double_t cent)
{
  Int_t bin = -1;
  for(Int_t i=0;i<kNCent_ds;i++)
    if ( (cent>=kBinCent_ds[i]) && (cent<kBinCent_ds[i+1]) )
      bin = i;

  return bin;

  }*/
//____________________________________________________________________________________________________

TObjArray* AliTwoParticlePIDCorr::GetV0Particles(AliVEvent* event,Double_t Centrality)
{//Only consider V0 particles having V0Pt>=2.0 GeV/C , REMEMBER that to chenge in case of change in Pt selection...Also note the Pt range in efficiency histogram

  AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(event);

  // Int_t curCentBin = GetCentBin(Centrality);//Now not used, previously used for centrality bin wise storage of ao, a1, a2- pol2 fit values of invariant mass fit sigmas

 //function to select v0's from AODs
  trkVtx=fAOD->GetPrimaryVertex();
  Float_t xv=trkVtx->GetX(), yv=trkVtx->GetY(), zv=trkVtx->GetZ();
  Int_t nV0sTot = fAOD->GetNumberOfV0s();

        TObjArray * selectedV0s = new TObjArray;
	selectedV0s->SetOwner(kTRUE);

 for (Int_t iV0 = 0; iV0 < nV0sTot; iV0++) 
    {
    
    AliAODv0 *v0=fAOD->GetV0(iV0);
    if (!v0) continue;
    if(!CheckStatusv0(v0)) continue;


     Float_t v0Pt=TMath::Sqrt(v0->Pt2V0());

    if (v0Pt< fminPtTrig  || v0Pt> fmaxPtTrig) continue;//*********************************IMPORTANT(depending on Pt range need for physics, Also note efficiency histo range)

    //(track->Pt()>=fminPtTrig && track->Pt()<=fmaxPtTrig)//to reduce memory consumption in Pool

    AliAODTrack *ptrack=(AliAODTrack*)v0->GetDaughter(0);
    AliAODTrack *ntrack=(AliAODTrack*)v0->GetDaughter(1);

    Bool_t cutK0sPID=kFALSE;
    Bool_t cutLambdaPID=kFALSE;
    Bool_t cutAntiLambdaPID=kFALSE;

    if(fUsev0DaughterPID)
{
	//use fHelperPID check PID of the daughter tracks
	//v0 daughter PID may be helpful in distangling k0S and (Anti)Lamda
                                                                                          
        Int_t PIDptrack = GetParticle(ptrack,kFALSE);
        Int_t PIDntrack = GetParticle(ntrack ,kFALSE);

        if(PIDptrack ==0 &&  PIDntrack == 0) cutK0sPID=kTRUE;

        if(PIDptrack==2 && PIDntrack ==0) cutLambdaPID=kTRUE;

        if (PIDptrack ==0 && PIDntrack == 2) cutAntiLambdaPID=kTRUE;

      }

 // effective mass calculations for each hypothesis(without daughter PID)
    Double_t InvMassK0s = v0->MassK0Short();
    Double_t InvMassAntiLambda = v0->MassAntiLambda();
    Double_t InvMassLambda = v0->MassLambda();

  
    
    Float_t v0Eta=v0->Eta();
    Float_t v0Phi=v0->Phi();


    
    //This is simply raw v0 without any specialised cut
    if(ffillofflineV0){
    fHistRawPtCentInvK0s->Fill(InvMassK0s,v0Pt,Centrality);
    fHistRawPtCentInvLambda->Fill(InvMassLambda,v0Pt,Centrality);
    fHistRawPtCentInvAntiLambda->Fill(InvMassAntiLambda,v0Pt,Centrality);
    }

  // Decay vertex
    Double_t xyz[3];   
    v0->GetSecondaryVtx(xyz);
    Float_t dx,dy,dz;
     dx=xyz[0]-xv, dy=xyz[1]-yv, dz=xyz[2]-zv;

    Float_t v0DecayRadius=TMath::Sqrt(dx*dx + dy*dy);
    Float_t v0DecayLength=TMath::Sqrt(dx*dx + dy*dy + dz*dz);
    // VO's main characteristics to check the reconstruction cuts
    // Float_t DcaV0Daughters    = v0->DcaV0Daughters();
    Float_t V0cosPointAngle   = v0->CosPointingAngle(trkVtx);
    // Float_t DcaPosToPrimVertex = v0->DcaPosToPrimVertex();
    //Float_t DcaNegToPrimVertex = v0->DcaNegToPrimVertex();   
    //Float_t Dcav0PVz   = v0->DcaV0ToPrimVertex(); 
    Float_t v0Pz=v0->Pz();
    Float_t v0P= TMath::Sqrt(v0Pt*v0Pt + v0Pz*v0Pz);

    Float_t ctauLambda =999.;
    Float_t ctauAntiLambda = 999.;
    Float_t ctauK0s = 999.;
 if(fCtauCut3D)
      {
	if(v0P > 0){
	 ctauLambda = (v0DecayLength*1.1157)/v0P;
	 ctauAntiLambda = (v0DecayLength*1.1157)/v0P;
	 ctauK0s = (v0DecayLength*0.4977)/v0P;
	}
      }


 else{
   if(v0Pt > 0.0){
          ctauLambda = (v0DecayRadius*1.1157)/v0Pt;
	 ctauAntiLambda = (v0DecayRadius*1.1157)/v0Pt;
	 ctauK0s = (v0DecayRadius*0.4977)/v0Pt;
	}
 }
    
 Bool_t ctauCutK0s= ctauK0s < (NCtau*fCutctauK0s) ; //ctauK0s 2.68 cm, mean life time of K0s is 8.95 x10^(-11)
 Bool_t ctauCutLambda = ctauLambda    < (NCtau*fCutctauLambda); //ctauLambda 7.89 cm ,mean life is 2.6 x10 ^(-10) ***** 3xctau is the accepted limit
 Bool_t ctauCutAntiLambda= ctauAntiLambda < (NCtau*fCutctauAntiLambda);

    Bool_t RapCutK0s = v0->RapK0Short() < fRapCutK0s;
    Bool_t RapCutLambda = v0->RapLambda() < fRapCutLambda; 
    Bool_t RapCutAntiLambda = v0->Y(-3122) < fRapCutLambda;

    Bool_t CPACut= V0cosPointAngle > fMinCPA; //cosine of pointing angle with v0 should be greater than 0.998

    //Now we put a loose mass cut which will be tightened later 
    Bool_t MassCutLooseK0s=(TMath::Abs(InvMassK0s - 0.497614) < 0.1);
    Bool_t MassCutLooseLambda=(TMath::Abs(InvMassLambda - 1.115683) < 0.1); // cut is same for Anti-Lambda
    Bool_t MassCutLooseAntiLambda=(TMath::Abs(InvMassAntiLambda - 1.115683) < 0.1); // cut is same for Anti-Lambda

    // cout<<"MassCutLooseK0s="<<MassCutLooseK0s<<"   "<<"MassCutLooseLambda="<<MassCutLooseLambda<<"   "<<"MassCutLooseAntiLambda="<<MassCutLooseAntiLambda<<endl;
 //Special Cut for Kshort arementeros podalanski plot
    Bool_t ArmenterosCut =kFALSE;
    if(ctauCutK0s && RapCutK0s && CPACut && MassCutLooseK0s)
      {

    Float_t lAlphaV0      =  v0->AlphaV0();
    Float_t lPtArmV0      =  v0->PtArmV0();

     ArmenterosCut = lPtArmV0 > TMath::Abs(0.2*lAlphaV0);

      }

    Bool_t IskShortOk=(ctauCutK0s && RapCutK0s && CPACut && MassCutLooseK0s && ArmenterosCut);

    Bool_t IsLambdaOk=(ctauCutLambda && RapCutLambda && CPACut && MassCutLooseLambda);

    Bool_t IsAntiLambdaOk=(ctauCutAntiLambda && RapCutAntiLambda && CPACut && MassCutLooseAntiLambda);

//Difference on Lambda and Anti-Lambda can be made through daughter PID

    Int_t particletype=999;

    // Disentangle the V0 candidate
    Double_t massK0s = 0., sK0s = 0.;
    Double_t massLambda = 0., sL = 0.;
    Double_t sAL = 0.;

    //if(IskShortOk && (IsLambdaOk || IsAntiLambdaOk)) cout<<"***********************************************This shouldn't happen"<<endl;
    
	if( IskShortOk || cutK0sPID )
	{
	  if(ffillofflineV0) fHistFinalPtCentInvK0s->Fill(InvMassK0s,v0Pt,Centrality);

	 massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();


	 sK0s=GetV0_MeanSigma_CentPt(Centrality,v0Pt,0);//0=kShort

	 // cout<<"InvMass="<<"Pt"<<v0Pt<<"*******************************************************sK0s="<<sK0s<<endl;
	 
         if (TMath::Abs(massK0s-InvMassK0s) < 3*sK0s) particletype=SpKs0;
    
	 if((massK0s-InvMassK0s) >= -8*sK0s && (massK0s-InvMassK0s)<= -5*sK0s) particletype=SpKs0_LS_Bckg; 
	 if((massK0s-InvMassK0s) <= 8*sK0s && (massK0s-InvMassK0s)>= 5*sK0s) particletype=SpKs0_RS_Bckg; 


    Short_t chargeval=0;
    Float_t effmatrix=1.0;
    LRCParticlePID* copy1 = new LRCParticlePID(particletype,InvMassK0s,chargeval,v0Pt,v0Eta, v0Phi,effmatrix,ptrack->GetTPCSharedMapPtr(),ntrack->GetTPCSharedMapPtr());
    copy1->SetUniqueID(eventno * 200000 + (Int_t)iV0);
    selectedV0s->Add(copy1);
  
	}
	//cout<<InvMassK0s<<"      "<<InvMassLambda<<"         "<<InvMassAntiLambda<<"     "<<particletype<<endl;


	if(IsLambdaOk ||  cutLambdaPID || IsAntiLambdaOk ||  cutAntiLambdaPID)
	{
	if(ffillofflineV0) fHistFinalPtCentInvLambda->Fill(InvMassLambda,v0Pt,Centrality);
//Add in the LRCParticle and give Lambda a tag 5

       massLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();


         sL=GetV0_MeanSigma_CentPt(Centrality,v0Pt,1);//1=Lambda

	 // cout<<"Pt"<<v0Pt<<"*******************************************************sL="<<sL<<endl;

	 
         if (TMath::Abs(massLambda-InvMassLambda) < 3*sL) particletype=SpLam;
    
	 if((massLambda-InvMassLambda) >= -8*sL && (massLambda-InvMassLambda)<= -5*sL) particletype=SpLam_LS_Bckg; 
	 if((massLambda-InvMassLambda) <= 8*sL && (massLambda-InvMassLambda)>= 5*sL) particletype=SpLam_RS_Bckg; 
	
	
    Short_t chargeval=0;
    Float_t effmatrix=1.0;
    LRCParticlePID* copy1 = new LRCParticlePID(particletype,InvMassLambda,chargeval,v0Pt,v0Eta, v0Phi,effmatrix,ptrack->GetTPCSharedMapPtr(),ntrack->GetTPCSharedMapPtr());
    copy1->SetUniqueID(eventno * 200000 + (Int_t)iV0);
    selectedV0s->Add(copy1);
	}
	//cout<<InvMassK0s<<"      "<<InvMassLambda<<"         "<<InvMassAntiLambda<<"     "<<particletype<<endl;

	/*
	if(IsAntiLambdaOk ||  cutAntiLambdaPID)
	{
	if(ffillofflineV0) fHistFinalPtCentInvAntiLambda->Fill(InvMassAntiLambda,v0Pt,Centrality);
//Add in the LRCParticle and give Lambda a tag 6

       massLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

        
                	  sAL = kAntiLambda_a0[curCentBin] + kAntiLambda_a1[curCentBin]*v0Pt;
       if(fpol2)	  sAL = kAntiLambda_a0[curCentBin] + kAntiLambda_a1[curCentBin]*v0Pt + kAntiLambda_a2[curCentBin]*(v0Pt*v0Pt);

	
        Bool_t ALamSignal = (TMath::Abs(massLambda-InvMassAntiLambda) < 3*sAL);
        Bool_t ALamBckg =( TMath::Abs(massLambda-InvMassAntiLambda + 6.5*sAL) < 1.5*sAL || TMath::Abs(massLambda-InvMassAntiLambda - 6.5*sAL) < 1.5*sAL  );
	
	if(ALamSignal) particletype=SpALam;
	if(ALamBckg)   particletype=SpALamBckg;
	
	Short_t chargeval=0;
    Float_t effmatrix=1.0;
    LRCParticlePID* copy1 = new LRCParticlePID(particletype,InvMassAntiLambda,chargeval,v0Pt,v0Eta, v0Phi,effmatrix,ptrack->GetTPCSharedMapPtr(),ntrack->GetTPCSharedMapPtr());
    copy1->SetUniqueID(eventno * 200000 + (Int_t)iV0);
    selectedV0s->Add(copy1);
	}
	*/

    }//v0 loop

  return selectedV0s;    
}

//___________________________________________________________________
  Bool_t AliTwoParticlePIDCorr :: CheckStatusv0Daughter(AliAODTrack *t1 ,AliAODTrack *t2)
  {
  if (!t1->IsOn(AliAODTrack::kTPCrefit) || !t2->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  // Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
if(t1->GetTPCClusterInfo(2,1)<fDaugNClsTPC || t2->GetTPCClusterInfo(2,1)<fDaugNClsTPC) return kFALSE ;

// ---------------- Fraction of TPC Shared Cluster 
     Float_t fracPosDaugTPCSharedMap = GetFractionTPCSharedCls(t1);
     Float_t fracNegDaugTPCSharedMap = GetFractionTPCSharedCls(t2);

 if(  (fracPosDaugTPCSharedMap > fFracTPCcls) || (fracNegDaugTPCSharedMap > fFracTPCcls) )
	return kFALSE;
  
  return kTRUE; 

  }
//___________________________________________________________________________________________
   
 Float_t AliTwoParticlePIDCorr :: GetFractionTPCSharedCls( AliAODTrack *track)
{
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for primaries
 
  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();

  return 1.*sharedMap.CountBits()/track->GetTPCNclsF();
  
}
//______________________________________________________________________
  Bool_t AliTwoParticlePIDCorr :: CheckStatusv0(AliAODv0 *v1)
  {

    // Offline reconstructed V0 only
    if (v1->GetOnFlyStatus()) return kFALSE;

     AliAODTrack *ptrack=(AliAODTrack *)v1->GetDaughter(0);
     AliAODTrack *ntrack=(AliAODTrack *)v1->GetDaughter(1);

    if(!ptrack || !ntrack) return kFALSE;

    if(ptrack->Charge()==-1 || ntrack->Charge()==1) return kFALSE; //remove wrongly identified charge pairs 

    if(ptrack->Charge()==0 || ntrack->Charge()==0) return kFALSE; //remove uncharged pairs

    if(ptrack->Charge() == ntrack->Charge()) return kFALSE; //remove like sign pairs 

    if(!CheckStatusv0Daughter(ptrack,ntrack)) return kFALSE;//daughters need to pass some basic cuts    

    if(TMath::Abs(ptrack->Eta()) > fmaxeta || TMath::Abs(ntrack->Eta()) > fmaxeta) return kFALSE; // remove daughters beyond eta bound |0.8|
    if(fCutDaughterPtV0){
    if(ptrack->Pt() < fMinPtDaughter || ntrack->Pt() < fMinPtDaughter) return kFALSE; // remove daughter tracks below minmum p |1.0 GeV/c|

    if(ptrack->Pt() > fMaxPtDaughter || ntrack->Pt() > fMaxPtDaughter) return kFALSE; // remove daughter tracks above maximum p ** to make it compatiable with AliHelperPID**|4.0 GeV/C|
    }

    // Daughters: Impact parameter of daughter to prim vtx
    Float_t xy = v1->DcaPosToPrimVertex();
    if (TMath::Abs(xy)<fDCAToPrimVtx) return kFALSE; //0.1 cm
    xy = v1->DcaNegToPrimVertex();
    if (TMath::Abs(xy)<fDCAToPrimVtx) return kFALSE; //0.1 cm

    // Daughters: DCA
    Float_t dca = v1->DcaV0Daughters();
    if (dca>fMaxDCADaughter) return kFALSE; //1.0 cm
    
    // V0: Cosine of the pointing angle
    Float_t cpa=v1->CosPointingAngle(trkVtx); //0.998
    if (cpa<fMinCPA) return kFALSE;

    // V0: Fiducial volume
    Double_t xyz[3]; v1->GetSecondaryVtx(xyz);
    Float_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
    if (r2<5.*5.) return kFALSE;
    if (r2>lMax*lMax) return kFALSE; //lmax=100 cm

    return kTRUE;


  }
//__________________________________________________________________________________________
Bool_t AliTwoParticlePIDCorr::IsTrackFromV0(AliAODEvent* fAOD,AliAODTrack* track)
{
//to check whether a daughter being taken as associated
	Int_t assoID = track->GetID();

	for(int i=0; i<fAOD->GetNumberOfV0s(); i++){ // loop over V0s
		AliAODv0* aodV0 = fAOD->GetV0(i);
		
		AliAODTrack *trackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
        	AliAODTrack *trackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
			
		
		Int_t negID = trackNeg->GetID();
		Int_t posID = trackPos->GetID();
		
		if ((TMath::Abs(negID)+1)==(TMath::Abs(assoID))){ return kTRUE;}
		if ((TMath::Abs(posID)+1)==(TMath::Abs(assoID))){ return kTRUE;}
		//----------------------------------
	}
	return kFALSE;
}

//________________________________________________________________________
Bool_t AliTwoParticlePIDCorr::TPCCutMIGeo( AliAODTrack* track,  AliAODEvent* evt)
{
  //
  // TPC Cut MIGeo
  //

  if (!track || !evt)
    return kFALSE;

Double_t fgCutGeo = 1.;// Cut variable for TPCCutMIGeo concerning geometry   
Double_t fgCutNcr = 0.85;// Cut variable for TPCCutMIGeo concerning num crossed rows 
Double_t fgCutNcl = 0.7;// Cut variable for TPCCutMIGeo concerning num clusters  

//UShort_t fgCutPureNcl = 60;

  const Short_t sign = track->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[100];

  track->GetXYZ(xyz);
  track->GetPxPyPz(pxpypz);

  AliExternalTrackParam* par = new AliExternalTrackParam(xyz, pxpypz, cv, sign);
  const AliESDtrack dummy;

  const Double_t magField = evt->GetMagneticField();
  Double_t varGeom = dummy.GetLengthInActiveZone(par, 3, 236, magField, 0, 0);
  Double_t varNcr  = track->GetTPCClusterInfo(3, 1);
  Double_t varNcls = track->GetTPCsignalN();

  const Double_t varEval = 130. - 5. * TMath::Abs(1. / track->Pt());
  Bool_t cutGeom   = varGeom > fgCutGeo * varEval;
  Bool_t cutNcr    = varNcr  > fgCutNcr * varEval;
  Bool_t cutNcls   = varNcls > fgCutNcl * varEval;
  // Bool_t cutNcls1 =(track->GetTPCsignalN() >= fgCutPureNcl);//Is it required??
  
  Bool_t kout = cutGeom && cutNcr && cutNcls;
  
  delete par;
  
  return kout;
}


//________________________________________________________________________


//----------------------------------------------------------
void AliTwoParticlePIDCorr::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
  
  
}
//------------------------------------------------------------------ 

