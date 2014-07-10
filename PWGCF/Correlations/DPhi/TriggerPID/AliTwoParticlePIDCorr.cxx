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

#include "AliCentrality.h"
#include "Riostream.h"

#include "AliTHn.h"    
#include "AliCFContainer.h"
#include "THn.h"
#include "THnSparse.h"

#include <TSpline.h>
#include <AliPID.h>
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include <AliPIDResponse.h>
#include "AliPIDCombined.h"   

#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include "AliAODInputHandler.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"

#include "THnSparse.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "TParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliOADBContainer.h"

#include "AliEventPoolManager.h"
#include "AliAnalysisUtils.h"
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
  fSampleType("pPb"),
 fRequestEventPlane(kFALSE),
  fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default)
  trkVtx(0),
  zvtx(0),
  fFilterBit(768),
  fTrackStatus(0),
  fSharedClusterCut(-1),
  fVertextype(1),
 skipParticlesAbove(0),
  fzvtxcut(10.0),
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
  fminPt(0.2),
  fmaxPt(20.0),
  fmineta(-0.8),
  fmaxeta(0.8),
  fselectprimaryTruth(kTRUE),
  fonlyprimarydatareco(kFALSE),
  fdcacut(kFALSE),
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
  fhistcentrality(0),
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
 ftwoTrackEfficiencyCutDataReco(kTRUE),
  twoTrackEfficiencyCutValue(0.02),
  fPID(NULL),
 fPIDCombined(NULL),
 eventno(0),
  fPtTOFPIDmin(0.5),
  fPtTOFPIDmax(4.0),
  fRequestTOFPID(kTRUE),
  fPIDType(NSigmaTPCTOF),
 fFIllPIDQAHistos(kTRUE),
  fNSigmaPID(3),
  fBayesCut(0.8),
 fdiffPIDcutvalues(kFALSE),
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
fAnalysisUtils(0x0),
  fDCAXYCut(0)     

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
  fSampleType("pPb"),
 fRequestEventPlane(kFALSE),
  fnTracksVertex(1),  // QA tracks pointing to principal vertex (= 3 default)
  trkVtx(0),
  zvtx(0),
  fFilterBit(768),
  fTrackStatus(0),
  fSharedClusterCut(-1),
  fVertextype(1),
   skipParticlesAbove(0),
  fzvtxcut(10.0),
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
  fminPt(0.2),
  fmaxPt(20.0),
  fmineta(-0.8),
  fmaxeta(0.8),
  fselectprimaryTruth(kTRUE),
  fonlyprimarydatareco(kFALSE),
   fdcacut(kFALSE),
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
  fhistcentrality(0),
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
  ftwoTrackEfficiencyCutDataReco(kTRUE),
  twoTrackEfficiencyCutValue(0.02),
  fPID(NULL),
  fPIDCombined(NULL),
  eventno(0),
 fPtTOFPIDmin(0.5),
  fPtTOFPIDmax(4.0),
  fRequestTOFPID(kTRUE),
  fPIDType(NSigmaTPCTOF),
  fFIllPIDQAHistos(kTRUE),
  fNSigmaPID(3),
  fBayesCut(0.8),
 fdiffPIDcutvalues(kFALSE),
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
   fAnalysisUtils(0x0),
  fDCAXYCut(0)    
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
 
  DefineOutput(1, TList::Class());                                        // for output list
  DefineOutput(2, TList::Class());

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
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPID = inputHandler->GetPIDResponse();

  //AliAnalysisUtils *fUtils = new AliAnalysisUtils();

//get the efficiency correction map

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

  fList = new TList;
  fList->SetOwner();
  fList->SetName("EPQAList");
  
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
  
fEtaSpectrasso=new TH2F("fEtaSpectraasso","fEtaSpectraasso",180,-0.9,0.9,100,0.,20. );
fOutput->Add(fEtaSpectrasso);

fphiSpectraasso=new TH2F("fphiSpectraasso","fphiSpectraasso",72,0,2*TMath::Pi(),100,0.,20.);
fOutput->Add(fphiSpectraasso);

 if(fSampleType=="pPb" || fSampleType=="PbPb"){ fCentralityCorrelation = new TH2D("fCentralityCorrelation", ";centrality;multiplicity", 101, 0, 101, 20000, 0,40000);
      fOutput->Add(fCentralityCorrelation);
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

if(fCentralityMethod.EndsWith("_MANUAL"))
  {
fhistcentrality=new TH1F("fhistcentrality","referencemultiplicity",30001,-0.5,30000.5);
fOutput->Add(fhistcentrality);
  }
 else{
fhistcentrality=new TH1F("fhistcentrality","centrality",220,-5,105);
fOutput->Add(fhistcentrality);
 }

if(fCentralityMethod.EndsWith("_MANUAL"))
  {
TString gmultName[4] = {"V0A_MANUAL","V0C_MANUAL","V0M_MANUAL","TRACKS_MANUAL"};
  fHistRefmult = new TH2F("fHistRefmult",
                             "Reference multiplicity",
			    4,-0.5,3.5,10000,0,20000);
  for(Int_t i = 1; i <= 4; i++){
    fHistRefmult->GetXaxis()->SetBinLabel(i,gmultName[i-1].Data());
    //fHistCentStatsUsed->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  }
  fOutput->Add(fHistRefmult);


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

fHistoTPCdEdx = new TH2F("fHistoTPCdEdx", ";p_{T} (GeV/c);dE/dx (au.)",200,0.0,10.0,500, 0., 500.);
fOutputList->Add(fHistoTPCdEdx);
fHistoTOFbeta = new TH2F(Form("fHistoTOFbeta"), ";p_{T} (GeV/c);v/c",100, 0., fmaxPt, 500, 0.1, 1.1);
  fOutputList->Add(fHistoTOFbeta);
  
   fTPCTOFPion3d=new TH3F ("fTPCTOFpion3d", "fTPCTOFpion3d",100,0., 10., 120,-60.,60.,120,-60.,60);
   fOutputList->Add(fTPCTOFPion3d);
  
   fTPCTOFKaon3d=new TH3F ("fTPCTOFKaon3d", "fTPCTOFKaon3d",100,0., 10., 120,-60.,60.,120,-60.,60);
   fOutputList->Add(fTPCTOFKaon3d);

   fTPCTOFProton3d=new TH3F ("fTPCTOFProton3d", "fTPCTOFProton3d",100,0., 10., 120,-60.,60.,120,-60.,60);
   fOutputList->Add(fTPCTOFProton3d);

if(ffillhistQAReco)
    {
    fPionPt = new TH1F("fHistQAPionPt","p_{T} distribution",200,0.,10.);
 fOutputList->Add(fPionPt);
    fPionEta= new TH1F("fHistQAPionEta","#eta distribution",360,-1.8,1.8);
 fOutputList->Add(fPionEta);
    fPionPhi = new TH1F("fHistQAPionPhi","#phi distribution",340,0,6.8);
 fOutputList->Add(fPionPhi);
  
    fKaonPt = new TH1F("fHistQAKaonPt","p_{T} distribution",200,0.,10.);
 fOutputList->Add(fKaonPt);
    fKaonEta= new TH1F("fHistQAKaonEta","#eta distribution",360,-1.8,1.8);
 fOutputList->Add(fKaonEta);
    fKaonPhi = new TH1F("fHistQAKaonPhi","#phi distribution",340,0,6.8);
 fOutputList->Add(fKaonPhi);
  
    fProtonPt = new TH1F("fHistQAProtonPt","p_{T} distribution",200,0.,10.);
 fOutputList->Add(fProtonPt);
    fProtonEta= new TH1F("fHistQAProtonEta","#eta distribution",360,-1.8,1.8);
 fOutputList->Add(fProtonEta);
    fProtonPhi= new TH1F("fHistQAProtonPhi","#phi distribution",340,0,6.8);
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
  fHistQA[8] = new TH1F("fHistQAPt","p_{T} distribution",200,0.,10.);
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
    "p_t_eff:0.0,0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0,5.5, 6.0, 7.0, 8.0,9.0,10.0\n"
    "vertex: -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10\n"
  "delta_phi: -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389\n" // this binning starts at -pi/2 and is modulo 3 
	"delta_eta: -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,2.1, 2.2, 2.3, 2.4\n"
      "multiplicity: 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1\n";

 if(fRequestEventPlane){
   defaultBinningStr += "eventPlane: -0.5,0.5,1.5,2.5,3.5\n"; // Event Plane Bins (Psi: -0.5->0.5 (in plane), 0.5->1.5 (intermediate), 1.5->2.5 (out of plane), 2.5->3.5 (rest))
  }

  if(fcontainPIDtrig){
      defaultBinningStr += "PIDTrig: -0.5,0.5,1.5,2.5,3.5\n"; // course
  }
  if(fcontainPIDasso){
      defaultBinningStr += "PIDAsso: -0.5,0.5,1.5,2.5,3.5\n"; // course
  }
 
  if(SetChargeAxis==2){
      defaultBinningStr += "TrigCharge: -2.0,0.0,2.0\n"; // course
      defaultBinningStr += "AssoCharge: -2.0,0.0,2.0\n"; // course
  }
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
    axisTitlePair[0]   = " multiplicity";

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


	fminPtTrig=dBinsPair[2][0];
        fmaxPtTrig=dBinsPair[2][iBinPair[2]];
        fminPtAsso=dBinsPair[3][0];
        fmaxPtAsso=dBinsPair[3][iBinPair[3]];
        fmincentmult=dBinsPair[0][0];
        fmaxcentmult=dBinsPair[0][iBinPair[0]];

	//event pool manager
Int_t MaxNofEvents=1000;
const Int_t NofVrtxBins=10+(1+10)*2;
Double_t ZvrtxBins[NofVrtxBins+1]={ -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10, 
				       90,  92,  94,  96,  98, 100, 102, 104, 106, 108, 110, 
				    190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210}; 

if(fCentralityMethod.EndsWith("_MANUAL"))
   {
 const Int_t NofCentBins=9;
 Double_t CentralityBins[NofCentBins+1]={0.,9.,14.,19.,26.,34.,44.,58.,80.,500.};//Is This binning is fine for pp, or we don't require them....
if(fRequestEventPlane){
    // Event plane angle (Psi) bins
  /*
    Double_t* psibins = NULL;
    Int_t nPsiBins; 
    psibins = GetBinning(fBinningString, "eventPlane", nPsiBins);
  */
 const Int_t  nPsiBins=6;
 Double_t psibins[nPsiBins+1]={0.0*TMath::DegToRad(), 30.0*TMath::DegToRad(), 60.0*TMath::DegToRad(), 90.0*TMath::DegToRad(), 120.0*TMath::DegToRad(),150.0*TMath::DegToRad(),180.1*TMath::DegToRad()};
fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins, nPsiBins, psibins);
// if(psibins)  delete [] psibins; 
				    }

 else{
 const Int_t  nPsiBinsd=1;
 Double_t psibinsd[nPsiBinsd+1]={0.0, 2000.0};
fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins, nPsiBinsd, psibinsd);

// fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins);
 }  
fPoolMgr->SetTargetValues(fMaxNofMixingTracks, 0.1, 5);

   }
 else
   {
 const Int_t  NofCentBins=15;
Double_t CentralityBins[NofCentBins+1]={0., 1., 2., 3., 4., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1 };
 if(fRequestEventPlane){
    // Event plane angle (Psi) bins
   /*
    Double_t* psibins = NULL;
    Int_t nPsiBins; 
    psibins = GetBinning(fBinningString, "eventPlane", nPsiBins);
   */
 const Int_t  nPsiBins=6;
 Double_t psibins[nPsiBins+1]={0.0*TMath::DegToRad(), 30.0*TMath::DegToRad(), 60.0*TMath::DegToRad(), 90.0*TMath::DegToRad(), 120.0*TMath::DegToRad(),150.0*TMath::DegToRad(),180.1*TMath::DegToRad()};
fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins, nPsiBins, psibins);
// if(psibins)  delete [] psibins; 
				    }

 else{
const Int_t  nPsiBinsd=1;
 Double_t psibinsd[nPsiBinsd+1]={0.0, 2000.0};
fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins, nPsiBinsd, psibinsd);

//fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins);
 }  
fPoolMgr->SetTargetValues(fMaxNofMixingTracks, 0.1, 5);
   }

 
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
     if((fAnalysisType == "MCAOD") && ffilltrigIDassoIDMCTRUTH) {//remember that in this case uidentified means other than pions, kaons, protons

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
 
  //ThSparse for trigger counting(data & reco MC)
  if(ffilltrigassoUNID || ffilltrigUNIDassoID || ffilltrigIDassoUNID || ffilltrigIDassoID)
	  {
	    fTHnTrigcount = new  AliTHn("fTHnTrigcount", "fTHnTrigcount", 2, dims, fBinst); //2 steps;;;;0->same event;;;;;1->mixed event
   for(Int_t i=0; i<dims;i++){
    fTHnTrigcount->SetBinLimits(i, dBinsTrig[i]);
    fTHnTrigcount->SetVarTitle(i, axisTitleTrig[i]);
  } 
  fOutput->Add(fTHnTrigcount);
	  }
  
  if((fAnalysisType =="MCAOD") && ffilltrigIDassoIDMCTRUTH) {
  //AliTHns for trigger counting(truth MC)
  fTHnTrigcountMCTruthPrim = new  AliTHn("fTHnTrigcountMCTruthPrim", "fTHnTrigcountMCTruthPrim", 2, dims, fBinst); //2 steps;;;;0->same event;;;;;1->mixed event
 for(Int_t i=0; i<dims;i++){
    fTHnTrigcount->SetBinLimits(i, dBinsTrig[i]);
    fTHnTrigcount->SetVarTitle(i, axisTitleTrig[i]);
  } 
  fOutput->Add(fTHnTrigcountMCTruthPrim);
 }

if(fAnalysisType=="MCAOD"){
  if(ffillhistQATruth)
    {
  MCtruthpt=new TH1F ("MCtruthpt","ptdistributiontruthprim",100,0.,10.);
  fOutputList->Add(MCtruthpt);

  MCtrutheta=new TH1F ("MCtrutheta","etadistributiontruthprim",360,-1.8,1.8);
  fOutputList->Add(MCtrutheta);

  MCtruthphi=new TH1F ("MCtruthphi","phidisttruthprim",340,0,6.8);
  fOutputList->Add(MCtruthphi);

  MCtruthpionpt=new TH1F ("MCtruthpionpt","MCtruthpionpt",100,0.,10.);
  fOutputList->Add(MCtruthpionpt);

  MCtruthpioneta=new TH1F ("MCtruthpioneta","MCtruthpioneta",360,-1.8,1.8);
  fOutputList->Add(MCtruthpioneta);

  MCtruthpionphi=new TH1F ("MCtruthpionphi","MCtruthpionphi",340,0,6.8);
  fOutputList->Add(MCtruthpionphi);

  MCtruthkaonpt=new TH1F ("MCtruthkaonpt","MCtruthkaonpt",100,0.,10.);
  fOutputList->Add(MCtruthkaonpt);

  MCtruthkaoneta=new TH1F ("MCtruthkaoneta","MCtruthkaoneta",360,-1.8,1.8);
  fOutputList->Add(MCtruthkaoneta);

  MCtruthkaonphi=new TH1F ("MCtruthkaonphi","MCtruthkaonphi",340,0,6.8);
  fOutputList->Add(MCtruthkaonphi);

  MCtruthprotonpt=new TH1F ("MCtruthprotonpt","MCtruthprotonpt",100,0.,10.);
  fOutputList->Add(MCtruthprotonpt);

  MCtruthprotoneta=new TH1F ("MCtruthprotoneta","MCtruthprotoneta",360,-1.8,1.8);
  fOutputList->Add(MCtruthprotoneta);

  MCtruthprotonphi=new TH1F ("MCtruthprotonphi","MCtruthprotonphi",340,0,6.8);
  fOutputList->Add(MCtruthprotonphi);
    }
 fPioncont=new TH2F("fPioncont", "fPioncont",10,-0.5,9.5,100,0.,10.);
  fOutputList->Add(fPioncont);

 fKaoncont=new TH2F("fKaoncont","fKaoncont",10,-0.5,9.5,100,0.,10.);
  fOutputList->Add(fKaoncont);

 fProtoncont=new TH2F("fProtoncont","fProtoncont",10,-0.5,9.5,100,0.,10.);
  fOutputList->Add(fProtoncont);

fUNIDcont=new TH2F("fUNIDcont","fUNIDcont",10,-0.5,9.5,100,0.,10.);
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

 
  //nsigma plot
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-30;
      Double_t maxy=30;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=50;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigma_%d_%d",ipart,ipid),Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
  }
  
  //nsigmaRec plot
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaRec_%d_%d",ipart,ipid),
				  Form("n#sigma for reconstructed %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
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
			       Form("probability for reconstructed %s",kParticleSpeciesName[ipart]),200,0,10,500,miny,maxy);
    fHistoBayes->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayes->GetYaxis()->SetTitle(Form("Bayes prob %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayes);


   TH2F *fHistoBayesTPC=new TH2F(Form("probBayes_TPC_%d",ipart),
			       Form("probability for Tracks as %s",kParticleSpeciesName[ipart]),200,0,10,500,miny,maxy);
    fHistoBayesTPC->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayesTPC->GetYaxis()->SetTitle(Form("Bayes prob TPC %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayesTPC);

  TH2F *fHistoBayesTOF=new TH2F(Form("probBayes_TOF_%d",ipart),
			       Form("probability for Tracks as %s",kParticleSpeciesName[ipart]),200,0,10,500,miny,maxy);
    fHistoBayesTOF->GetXaxis()->SetTitle("P_{T} (GeV / c)");
    fHistoBayesTOF->GetYaxis()->SetTitle(Form("Bayes prob TOF %s",kParticleSpeciesName[ipart]));
    fOutputList->Add(fHistoBayesTOF);

 TH2F *fHistoBayesTPCTOF=new TH2F(Form("probBayes_TPCTOF_%d",ipart),
			       Form("probability for Tracks as  %s",kParticleSpeciesName[ipart]),200,0,10,500,miny,maxy);
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
			       Form("Pi_Ka separation in %s",kPIDTypeName[ipid]),50,0,10,200,miny,maxy);
    Pi_Ka_sep->GetXaxis()->SetTitle("P_{T} (GeV/C)");
    Pi_Ka_sep->GetYaxis()->SetTitle(Form("expected seaparation(n#sigma) in %s",kPIDTypeName[ipid]));
    fOutputList->Add(Pi_Ka_sep);

   TH2F *Pi_Pr_sep=new TH2F(Form("Pi_Pr_sep_%d",ipid),
			       Form("Pi_Pr separation in %s",kPIDTypeName[ipid]),50,0,10,200,miny,maxy);
    Pi_Pr_sep->GetXaxis()->SetTitle("P_{T} (GeV/C)");
    Pi_Pr_sep->GetYaxis()->SetTitle(Form("expected seaparation(n#sigma) in %s",kPIDTypeName[ipid]));
    fOutputList->Add(Pi_Pr_sep);

    TH2F *Ka_Pr_sep=new TH2F(Form("Ka_Pr_sep_%d",ipid),
			       Form("Ka_Pr separation in %s",kPIDTypeName[ipid]),50,0,10,200,miny,maxy);
    Ka_Pr_sep->GetXaxis()->SetTitle("P_{T} (GeV/C)");
    Ka_Pr_sep->GetYaxis()->SetTitle(Form("expected seaparation(n#sigma) in %s",kPIDTypeName[ipid]));
    fOutputList->Add(Ka_Pr_sep);
    }

  //nsigmaDC plot
  if(fUseExclusiveNSigma) {
  for(Int_t ipart=0;ipart<NSpecies;ipart++){
    for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==NSigmaTPCTOF){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaDC_%d_%d",ipart,ipid),
				  Form("n#sigma for double counting %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]));
      fOutputList->Add(fHistoNSigma);
    }
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
				  Form("n#sigma for MC %s %s",kParticleSpeciesName[ipart],kPIDTypeName[ipid]),200,0,10,500,miny,maxy);
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
      TH2F *fHistoPID=new TH2F(Form("PID_%d_%d",idet,ipart),Form("%s signal - %s",kDetectorName[idet],kParticleSpeciesName[ipart]),200,0,10,500,-maxy,maxy);
      fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
      fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kDetectorName[idet]));
      fOutputList->Add(fHistoPID);
    }
  }
  //PID signal plot, before PID cut
  for(Int_t idet=0;idet<fNDetectors;idet++){
    Double_t maxy=500;
    if(idet==fTOF)maxy=1.1;
    TH2F *fHistoPID=new TH2F(Form("PIDAll_%d",idet),Form("%s signal",kDetectorName[idet]),200,0,10,500,-maxy,maxy);
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

  else if(fAnalysisType == "MCAOD") {
  
    doMCAODevent();
    
  }
  
  else return;
  
}
//-------------------------------------------------------------------------
void AliTwoParticlePIDCorr::doMCAODevent() 
{
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    AliError("Cannot get the AOD event");
    return;
  }
 
// count all events(physics triggered)   
  fEventCounter->Fill(1);

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

//check the PIDResponse handler
  if (!fPID) return;

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

 if (fSampleType=="pp" && fCentralityMethod.EndsWith("_MANUAL"))
   {
 //make the event selection with reco vertex cut and centrality cut and return the value of the centrality
     Double_t refmultTruth = GetAcceptedEventMultiplicity(aod,kTRUE);  //incase of ref multiplicity it will return the truth MC ref mullt value; need to determine the ref mult value separately for reco Mc case; in case of centrality this is final and fine
     refmultReco = GetAcceptedEventMultiplicity(aod,kFALSE); 
     if(refmultTruth<=0 || refmultReco<=0) return;
     cent_v0=refmultTruth;
   }
 else {
 cent_v0=GetAcceptedEventMultiplicity(aod,kTRUE); //centrality value; 2nd argument has no meaning
 if(cent_v0<0.) return;
 }

 effcent=cent_v0;// This will be required for efficiency THn filling(specially in case of pp)

  //get the event plane in case of PbPb
   if(fRequestEventPlane){
   gReactionPlane=GetEventPlane(aod,kTRUE,cent_v0);//get the truth event plane
   if(gReactionPlane==999.) return;
 }
  

   Double_t nooftrackstruth=0.0;//in case of pp this will give the multiplicity(for truth case) after the track loop(only for unidentified particles that pass  kinematic cuts)

   //TObjArray* tracksMCtruth=0;
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
 //get particle ID
Int_t pdgtruth=((AliAODMCParticle*)partMC)->GetPdgCode();
Int_t particletypeTruth=-999;
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
 if (fSampleType=="pp" && fCentralityMethod.EndsWith("_MANUAL")) effcent=15.0;//integrated over multiplicity(so put any fixed value for each track so that practically means there is only one bin in multiplicity i.e. multiplicity intregated out )**************Important

 Double_t primmctruth[4] = {effcent, zVtxmc,partMC->Pt(), partMC->Eta()};
 if(ffillefficiency)
  {
    fTrackHistEfficiency[5]->Fill(primmctruth,0);//for all primary truth particles(4)
    if (TMath::Abs(pdgtruth)==211 || TMath::Abs(pdgtruth)==321) fTrackHistEfficiency[3]->Fill(primmctruth,0);//for  primary truth mesons(3)
    if (TMath::Abs(pdgtruth)==2212 || TMath::Abs(pdgtruth)==321) fTrackHistEfficiency[4]->Fill(primmctruth,0);//for  primary truth kaons+protons(4)
    if (TMath::Abs(pdgtruth)==211)  fTrackHistEfficiency[0]->Fill(primmctruth,0);//for pions
    if (TMath::Abs(pdgtruth)==321)  fTrackHistEfficiency[1]->Fill(primmctruth,0);//for kaons
    if(TMath::Abs(pdgtruth)==2212)  fTrackHistEfficiency[2]->Fill(primmctruth,0);//for protons
  }
   
 Float_t effmatrixtruth=1.0;//In Truth MC, no case of efficiency correction so it should be always 1.0
if((partMC->Pt()>=fminPtAsso && partMC->Pt()<=fmaxPtAsso) || (partMC->Pt()>=fminPtTrig && partMC->Pt()<=fmaxPtTrig))//to reduce memory consumption in pool
  {
    Short_t chargeval=0;
    if(partMC->Charge()>0)   chargeval=1;
    if(partMC->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
LRCParticlePID* copy6 = new LRCParticlePID(particletypeTruth,chargeval,partMC->Pt(),partMC->Eta(), partMC->Phi(),effmatrixtruth);
//copy6->SetUniqueID(eventno * 100000 + TMath::Abs(partMC->GetLabel()));
 copy6->SetUniqueID(eventno * 100000 + (Int_t)nooftrackstruth);
 tracksMCtruth->Add(copy6);//************** TObjArray used for truth correlation function calculation
  }
  }//MC truth track loop ends

//*********************still in event loop

 if (fSampleType=="PbPb"){
   if (fRandomizeReactionPlane)//only for TRuth MC??
  {
    Double_t centralityDigits = cent_v0*1000. - (Int_t)(cent_v0*1000.);
    Double_t angle = TMath::TwoPi() * centralityDigits;
    AliInfo(Form("Shifting phi of all tracks by %f (digits %f)", angle, centralityDigits));
    ShiftTracks(tracksMCtruth, angle);  
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

if(tracksMCtruth) delete tracksMCtruth;

//now deal with reco tracks

//detrmine the ref mult in case of Reco(not required if we get centrality info from AliCentrality)
 if (fSampleType=="pp" && fCentralityMethod.EndsWith("_MANUAL")) cent_v0=refmultReco;
 effcent=cent_v0;// This will be required for efficiency THn filling(specially in case of pp)

 if(fRequestEventPlane){
   gReactionPlane = GetEventPlane(aod,kFALSE,cent_v0);//get the reconstructed event plane
    if(gReactionPlane==999.) return;
 }
  
   //TObjArray* tracksUNID=0;
   TObjArray* tracksUNID = new TObjArray;
   tracksUNID->SetOwner(kTRUE);

   //TObjArray* tracksID=0;
   TObjArray* tracksID = new TObjArray;
   tracksID->SetOwner(kTRUE);


   Float_t bSign1=aod->GetHeader()->GetMagneticField() ;//used for reconstructed track dca cut


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
  trackscount++;

//check for eta , phi holes
 fEtaSpectrasso->Fill(track->Eta(),track->Pt());
 fphiSpectraasso->Fill(track->Phi(),track->Pt());

  Int_t particletypeMC=-9999;

//tag all particles as unidentified
 particletypeMC=unidentified;

 Float_t effmatrix=1.;

// -- Fill THnSparse for efficiency calculation
 if (fSampleType=="pp" && fCentralityMethod.EndsWith("_MANUAL")) effcent=15.0;//integrated over multiplicity(so put any fixed value for each track so that practically means there is only one bin in multiplicity i.e. multiplicity intregated out )**************Important
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
   LRCParticlePID* copy = new LRCParticlePID(particletypeMC,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix);
   copy->SetUniqueID(eventno * 100000 +(Int_t)trackscount);
   tracksUNID->Add(copy);//track information Storage for UNID correlation function(tracks that pass the filterbit & kinematic cuts only)
  }

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

Float_t dEdx = track->GetTPCsignal();
 fHistoTPCdEdx->Fill(track->Pt(), dEdx);

 if(HasTOFPID(track))
{
Double_t beta = GetBeta(track);
fHistoTOFbeta->Fill(track->Pt(), beta);
 }

//do track identification(nsigma method)
  particletypeMC=GetParticle(track,fFIllPIDQAHistos);//******************************problem is here

switch(TMath::Abs(pdgCode)){
  case 2212:
    if(fFIllPIDQAHistos){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",SpProton,ipid));
	h->Fill(track->Pt(),fnsigmas[SpProton][ipid]);
      }
    }
    break;
  case 321:
    if(fFIllPIDQAHistos){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",SpKaon,ipid));
	h->Fill(track->Pt(),fnsigmas[SpKaon][ipid]);
      }
    }
    break;
  case 211:
    if(fFIllPIDQAHistos){
      for(Int_t ipid=0;ipid<=NSigmaPIDType;ipid++){
	if((ipid!=NSigmaTPC) && (!HasTOFPID(track)))continue;//not filling TOF and combined if no TOF PID
	TH2F *h=GetHistogram2D(Form("NSigmaMC_%d_%d",SpPion,ipid));
	h->Fill(track->Pt(),fnsigmas[SpPion][ipid]);
      }
    }
    break;
  }


//2-d TPCTOF map(for each Pt interval)
  if(HasTOFPID(track)){
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
 if(particletypeMC==SpUndefined )
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

if((track->Pt()>=fminPtAsso && track->Pt()<=fmaxPtAsso) || (track->Pt()>=fminPtTrig && track->Pt()<=fmaxPtTrig))//to reduce memory consumption in pool
  {
    Short_t chargeval=0;
    if(track->Charge()>0)   chargeval=1;
    if(track->Charge()<0)   chargeval=-1;
    if(chargeval==0) continue;
if (fapplyTrigefficiency || fapplyAssoefficiency)
    effmatrix=GetTrackbyTrackeffvalue(track,effcent,zvtx,particletypeMC);//get the tracking eff x TOF matching eff x PID eff x contamination factor for identified particles 
    LRCParticlePID* copy1 = new LRCParticlePID(particletypeMC,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix);
    copy1->SetUniqueID(eventno * 100000 + (Int_t)trackscount);
    tracksID->Add(copy1);
  }
  }// if(tracktype==1) condition structure ands

}//reco track loop ends

  //*************************************************************still in event loop
 

if(trackscount>0.0)
  { 
//fill the centrality/multiplicity distribution of the selected events
 fhistcentrality->Fill(cent_v0);//*********************************WARNING::binning of cent_v0 is different for pp and pPb/PbPb case

 if (fSampleType=="pPb" || fSampleType=="PbPb") fCentralityCorrelation->Fill(cent_v0, trackscount);//only with unidentified tracks(i.e before PID selection);;;;;can be used to remove centrality outliers??????

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

 //same event delte-eta, delta-phi plot
if(tracksUNID && tracksUNID->GetEntriesFast()>0)//hadron triggered correlation
  {//same event calculation starts
    if(ffilltrigassoUNID) Fillcorrelation(gReactionPlane,tracksUNID,0,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigassoUNID");//mixcase=kFALSE (hadron-hadron correlation)
    if(tracksID && tracksID->GetEntriesFast()>0 && ffilltrigUNIDassoID)  Fillcorrelation(gReactionPlane,tracksUNID,tracksID,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigUNIDassoID");//mixcase=kFALSE (hadron-ID correlation)
  }

if(tracksID && tracksID->GetEntriesFast()>0)//ID triggered correlation
  {//same event calculation starts
    if(tracksUNID && tracksUNID->GetEntriesFast()>0 && ffilltrigIDassoUNID)  Fillcorrelation(gReactionPlane,tracksID,tracksUNID,cent_v0,zvtx,weghtval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoUNID");//mixcase=kFALSE (ID-hadron correlation)
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
   Fillcorrelation(gReactionPlane,tracksID,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoUNID");//mixcase=kTRUE 
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

if(tracksUNID)  delete tracksUNID;

if(tracksID) delete tracksID;


PostData(1, fOutput);

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

// count all events   
  fEventCounter->Fill(1);

if (!fPID) return;//this should be available with each event even if we don't do PID selection

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
 Float_t bSign1=aod->GetHeader()->GetMagneticField() ;//for dca cut in ClassifyTrack(), i.e in track loop


// check event cuts and fill event histograms and return the centrality or reference multiplicity value
 if((cent_v0 = GetAcceptedEventMultiplicity(aod,kFALSE)) < 0){ 
    return;
  }
  
  //get the event plane in case of PbPb
    if(fRequestEventPlane){
      gReactionPlane = GetEventPlane(aod,kFALSE,cent_v0);
    if(gReactionPlane==999.) return;
    }    
  
   TObjArray*  tracksUNID= new TObjArray;//track info before doing PID
   tracksUNID->SetOwner(kTRUE);  // IMPORTANT!

   TObjArray* tracksID= new TObjArray;//only pions, kaons,protons i.e. after doing the PID selection
   tracksID->SetOwner(kTRUE);  // IMPORTANT!
 
    
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


 if (fSampleType=="pp") effcent=15.0;//integrated over multiplicity [i.e each track has multiplicity 15.0](so put any fixed value for each track so that practically means there is only one bin in multiplicityi.e multiplicity intregated out )**************Important for efficiency related issues


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
 LRCParticlePID* copy = new LRCParticlePID(particletype,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix);
  copy->SetUniqueID(eventno * 100000 + (Int_t)trackscount);
  tracksUNID->Add(copy);//track information Storage for UNID correlation function(tracks that pass the filterbit & kinematic cuts only)
  }

//now start the particle identificaion process:) 

//track passing filterbit 768 have proper TPC response,or need to be checked explicitly before doing PID????

  Float_t dEdx = track->GetTPCsignal();
  fHistoTPCdEdx->Fill(track->Pt(), dEdx);

  //fill beta vs Pt plots only for tracks having proper TOF response(much less tracks compared to the no. that pass the filterbit & kinematic cuts)
 if(HasTOFPID(track))
{
  Double_t beta = GetBeta(track);
  fHistoTOFbeta->Fill(track->Pt(), beta);
 }
  

//track identification(using nsigma method)
     particletype=GetParticle(track,fFIllPIDQAHistos);//*******************************change may be required(It should return only pion,kaon, proton and Spundefined; NOT unidentifed***************be careful)


//2-d TPCTOF map(for each Pt interval)
  if(HasTOFPID(track)){
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
 LRCParticlePID* copy1 = new LRCParticlePID(particletype,chargeval,track->Pt(),track->Eta(), track->Phi(),effmatrix);
    copy1->SetUniqueID(eventno * 100000 + (Int_t)trackscount);
    tracksID->Add(copy1);
  }
} //track loop ends but still in event loop

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

if(fSampleType=="pPb" || fSampleType=="PbPb") fCentralityCorrelation->Fill(cent_v0, trackscount);//only with unidentified tracks(i.e before PID selection);;;;;can be used to remove centrality outliers??????

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


//same event delta-eta-deltaphi plot 

if(tracksUNID && tracksUNID->GetEntriesFast()>0)//hadron triggered correlation
  {//same event calculation starts
    if(ffilltrigassoUNID) Fillcorrelation(gReactionPlane,tracksUNID,0,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigassoUNID");//mixcase=kFALSE (hadron-hadron correlation)
    if(tracksID && tracksID->GetEntriesFast()>0 && ffilltrigUNIDassoID)  Fillcorrelation(gReactionPlane,tracksUNID,tracksID,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigUNIDassoID");//mixcase=kFALSE (hadron-ID correlation)
  }

if(tracksID && tracksID->GetEntriesFast()>0)//ID triggered correlation
  {//same event calculation starts
    if(tracksUNID && tracksUNID->GetEntriesFast()>0 && ffilltrigIDassoUNID)  Fillcorrelation(gReactionPlane,tracksID,tracksUNID,cent_v0,zvtx,weightval,kFALSE,bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kFALSE,"trigIDassoUNID");//mixcase=kFALSE (ID-hadron correlation)
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
   Fillcorrelation(gReactionPlane,tracksID,bgTracks,cent_v0,zvtx,nmix1,(jMix == 0),bSign,fPtOrderDataReco,ftwoTrackEfficiencyCutDataReco,kTRUE,"trigIDassoUNID");//mixcase=kTRUE 
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


if(tracksUNID)  delete tracksUNID;

if(tracksID) delete tracksID;


PostData(1, fOutput);

} // *************************event loop ends******************************************//_______________________________________________________________________
TObjArray* AliTwoParticlePIDCorr::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliDPhiBasicParticle which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    LRCParticlePID* particle = (LRCParticlePID*) tracks->UncheckedAt(i);
    LRCParticlePID* copy100 = new LRCParticlePID(particle->getparticle(),particle->Charge(), particle->Pt(),particle->Eta(), particle->Phi(), particle->geteffcorrectionval());
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
      Int_t particlepidtrig=trig->getparticle();
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
      if(fcontainPIDtrig && SetChargeAxis==0) trigval[4] = particlepidtrig;
      if(!fcontainPIDtrig && SetChargeAxis==2) trigval[4] = trig->Charge();
      if(fcontainPIDtrig && SetChargeAxis==2) {
      trigval[4] = particlepidtrig;
      trigval[5] = trig->Charge();
       }
      }

  if(!fRequestEventPlane){
      if(fcontainPIDtrig && SetChargeAxis==0) trigval[3] = particlepidtrig;
      if(!fcontainPIDtrig && SetChargeAxis==2) trigval[3] = trig->Charge();
      if(fcontainPIDtrig && SetChargeAxis==2) {
      trigval[3] = particlepidtrig;
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

   // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event,i.e. both Trig and asso TObjArray belongs to the same Pi range but say Trig is Unidentified but asso is identified then the serial no. wise particles are not same and and j==i doesn't aplly)
   // if (tracksasso && trig->IsEqual(asso))  continue;

  if (tracksasso && (trig->GetUniqueID()==asso->GetUniqueID())) continue;
          
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
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);

 const Float_t kLimit = twoTrackEfficiencyCutValue * 3;

	    Float_t dphistarminabs = 1e5;
	    Float_t dphistarmin = 1e5;

 if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
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
        vars[dimension]=particlepidtrig;
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
	vars[dimension]=particlepidtrig;
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

//--------------------------------------------------------------------------------
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
     if(effvalue==0.) effvalue=1.;

     return effvalue; 

}
//---------------------------------------------------------------------------------

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
// DCA XY
	if (fdcacut && fDCAXYCut)
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
	}

	if (fSharedClusterCut >= 0)
	{
	  Double_t frac = Double_t(((AliAODTrack*)track)->GetTPCnclsS()) / Double_t(((AliAODTrack*)track)->GetTPCncls());
	  if (frac > fSharedClusterCut)
	    return 0;
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

 if(HasTOFPID(track) && pt>fPtTOFPIDmin)
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
if(fRequestTOFPID && track->Pt()>fPtTOFPIDmin && (!HasTOFPID(track)) )return SpUndefined;//so any track having Pt>0.6 withot having proper TOF response will be defined as SpUndefined
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
  if(track->Pt()<=4) fNSigmaPID=fPIDCutval1;
  if(track->Pt()>4 && track->Pt()<=6) fNSigmaPID=fPIDCutval2;
  if(track->Pt()>6 && track->Pt()<=8) fNSigmaPID=fPIDCutval3;
  if(track->Pt()>8) fNSigmaPID=fPIDCutval4;
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

  
  Double_t probBayes[AliPID::kSPECIES];
  
  UInt_t detUsed= 0;
  if(HasTOFPID(trk) && trk->Pt()>fPtTOFPIDmin){//use TOF information
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
  if(trk->Pt()<=4) fBayesCut=fPIDCutval1;
  if(trk->Pt()>4 && trk->Pt()<=6) fBayesCut=fPIDCutval2;
  if(trk->Pt()>6 && trk->Pt()<=8) fBayesCut=fPIDCutval3;
  if(trk->Pt()>8) fBayesCut=fPIDCutval4;
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

  CalculateNSigmas(trk,FIllQAHistos);//fill the data member fnsigmas with the nsigmas value [ipart][iPID]


 //Do PID
  if(fPIDType==Bayes){//use bayesianPID
    
    if(!fPIDCombined) {
      AliFatal("PIDCombined object has to be set in the steering macro");
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
   // if(track->GetTPCsignalN() < 60) return kFALSE;//tracks with TPCsignalN< 60 have questionable dEdx,cutting on TPCsignalN > 70 or > 60 shouldn't make too much difference in statistics,also  it is IMO safe to use TPC also for MIPs.
   
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
  if(track->Pt()<=fPtTOFPIDmin) return kFALSE;
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
//_________________________________________________________________________
/*
void AliTwoParticlePIDCorr ::DefineEventPool()
{
Int_t MaxNofEvents=1000;
const Int_t NofVrtxBins=10+(1+10)*2;
Double_t ZvrtxBins[NofVrtxBins+1]={ -10,   -8,  -6,  -4,  -2,   0,   2,   4,   6,   8,  10, 
				       90,  92,  94,  96,  98, 100, 102, 104, 106, 108, 110, 
				      190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210 

//default values are for centrality
Int_t  NofCentBins=15;
Double_t CentralityBins[NofCentBins+1]={0., 1., 2., 3., 4., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.1 };

 if(fCentralityMethod.EndsWith("_MANUAL"))
   {
 Int_t  NofCentBins=9;
 CentralityBins[NofCentBins+1]={0.,9.,14.,19.,26.,34.,44.,58.,80.,500.};//Is This binning is fine for pp, or we don't require them....
   }
fPoolMgr = new AliEventPoolManager(MaxNofEvents,fMaxNofMixingTracks,NofCentBins,CentralityBins,NofVrtxBins,ZvrtxBins);



  
fPoolMgr->SetTargetValues(fMaxNofMixingTracks, 0.1, 5);

//if(!fPoolMgr) return kFALSE;
//return kTRUE;

}
*/
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
Double_t AliTwoParticlePIDCorr::GetReferenceMultiplicityVZEROFromAOD(AliAODEvent *event){
  //Function that returns the reference multiplicity from AODs (data or reco MC)
  //Different ref. mult. implemented: V0M, V0A, V0C, TPC
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


  if(fCentralityMethod == "TRACKS_MANUAL") 
    gRefMultiplicity = gRefMultiplicityTPC;
  else if(fCentralityMethod == "V0M_MANUAL")
    gRefMultiplicity = gRefMultiplicityVZERO;
  else if(fCentralityMethod == "V0A_MANUAL")
    gRefMultiplicity = gRefMultiplicityVZEROA;
  else if(fCentralityMethod == "V0C_MANUAL")
    gRefMultiplicity = gRefMultiplicityVZEROC;

      //ref mult QA
      fHistRefmult->Fill(0.,gRefMultiplicityVZEROA);
      fHistRefmult->Fill(1.,gRefMultiplicityVZEROC);
      fHistRefmult->Fill(2.,gRefMultiplicityVZERO);
      fHistRefmult->Fill(3.,gRefMultiplicityTPC);

  
  return gRefMultiplicity;
}

//-------------------------------------------------------------------------------------------------------
Double_t AliTwoParticlePIDCorr::GetRefMultiOrCentrality(AliAODEvent *event, Bool_t truth){

  if(!event) return -1;
  // get centrality object and check quality
  Double_t cent_v0=-1;
  Double_t nooftrackstruth=0;

if(fCentralityMethod=="V0M" || fCentralityMethod=="V0A" || fCentralityMethod=="V0C" || fCentralityMethod=="CL1" || fCentralityMethod=="ZNA" || fCentralityMethod=="V0AEq" || fCentralityMethod=="V0CEq" || fCentralityMethod=="V0MEq")//for PbPb, pPb, pp7TeV(still to be introduced)//data or RecoMC and also for TRUTH
    {
  AliCentrality *centralityObj=0;
  AliAODHeader *header = (AliAODHeader*) event->GetHeader();
  if(!header) return -1;
  centralityObj = header->GetCentralityP();
  // if (centrality->GetQuality() != 0) return ;

  if(centralityObj)
  {
  fHistCentStats->Fill(0.,centralityObj->GetCentralityPercentile("V0A"));
  fHistCentStats->Fill(1.,centralityObj->GetCentralityPercentile("V0C"));
  fHistCentStats->Fill(2.,centralityObj->GetCentralityPercentile("V0M"));
if(fSampleType=="pp")   fHistCentStats->Fill(3.,centralityObj->GetCentralityPercentile("V0AEq"));//only available for LHC10d at present (Quantile info)
if(fSampleType=="pp")   fHistCentStats->Fill(4.,centralityObj->GetCentralityPercentile("V0CEq"));//only available for LHC10d at present (Quantile info)
if(fSampleType=="pp")   fHistCentStats->Fill(5.,centralityObj->GetCentralityPercentile("V0MEq"));//only available for LHC10d at present (Quantile info)

if(fSampleType=="pPb" || fSampleType=="PbPb")      fHistCentStats->Fill(6.,centralityObj->GetCentralityPercentile("CL1"));
if(fSampleType=="pPb" || fSampleType=="PbPb")      fHistCentStats->Fill(7.,centralityObj->GetCentralityPercentile("ZNA")); 

      cent_v0 = centralityObj->GetCentralityPercentile(fCentralityMethod);
  }
  else cent_v0= -1;    
    }//centralitymethod condition

 else if(fCentralityMethod=="V0M_MANUAL" || fCentralityMethod=="V0A_MANUAL" || fCentralityMethod=="V0C_MANUAL" || fCentralityMethod=="TRACKS_MANUAL")//data or RecoMc and also for TRUTH
   {
     if(!truth){//for data or RecoMC
    cent_v0 = GetReferenceMultiplicityVZEROFromAOD(event);
   }//for data or RecoMC

    if(truth){//condition for TRUTH case
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


      if (fCentralityMethod=="V0M_MANUAL") {
	if(partMC->Eta() > 5.1 || partMC->Eta() < 2.8)    continue;
	if (partMC->Eta() < -3.7 || partMC->Eta() > -1.7) continue;
}
      else if (fCentralityMethod=="V0A_MANUAL") {
	if(partMC->Eta() > 5.1 || partMC->Eta() < 2.8)  continue;}
      else if (fCentralityMethod=="V0C_MANUAL") {
	if(partMC->Eta() > -1.7 || partMC->Eta() < -3.7)  continue;}
      else if (fCentralityMethod=="TRACKS_MANUAL") {
        if (partMC->Eta() < fmineta || partMC->Eta() > fmaxeta) continue;
        if (partMC->Pt() < fminPt ||  partMC->Pt() > fmaxPt) continue;
           }
      else{//basically returns the tracks manual case
//give only kinematic cuts at the truth level  
       if (partMC->Eta() < fmineta || partMC->Eta() > fmaxeta) continue;
       if (partMC->Pt() < fminPt ||  partMC->Pt() > fmaxPt) continue;
      }

 //To determine multiplicity in case of PP
 nooftrackstruth+= 1;;

 }//truth track loop ends
 cent_v0=nooftrackstruth;

    }//condition for TRUTH case

   }//end of MANUAL method

 else if ((fAnalysisType == "MCAOD") && (fCentralityMethod == "MC_b"))//TRUTH MC
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
     
      
     if (collGeometry)   cent_v0 = collGeometry->ImpactParameter();
      else cent_v0=-1.;
    }//end of Impact parameter method

//else return -1
 else cent_v0=-1.;

 return cent_v0;
}
//-----------------------------------------------------------------------------------------
Double_t AliTwoParticlePIDCorr::GetAcceptedEventMultiplicity(AliAODEvent *aod,Bool_t truth){

  if(!aod) return -1;

  Float_t gRefMultiplicity = -1.;

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

  //vertex selection(is it fine for PP?)
 if ( fVertextype==1){//for pPb basically if(!fAnalysisUtils->IsVertexSelected2013pA(aod)) return; 
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
  else if(fVertextype==2) {//for pp and pb-pb case , taken from Jan's code
	Int_t nVertex = aod->GetNumberOfVertices();
  	if( nVertex > 0 ) { 
     trkVtx = aod->GetPrimaryVertex();
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
  trkVtx = aod->GetPrimaryVertex();
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
 if(truth)  {gRefMultiplicity = GetRefMultiOrCentrality(aod,kTRUE);}//kTRUE-->for Truth case(only meaningful in case of ref multiplicity)

 else {gRefMultiplicity = GetRefMultiOrCentrality(aod,kFALSE);}//kFALSE-->for data and RecoMc case(only meaningful in case of ref multiplicity)

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

  return gRefMultiplicity;

}
//--------------------------------------------------------------------------------------------------------
Float_t AliTwoParticlePIDCorr::GetEventPlane(AliAODEvent *event,Bool_t truth, Double_t v0Centr)
{
  // Get the event plane
 //reset Q vector info	

    Int_t run = event->GetRunNumber();

    if(run != fRun){
	// Load the calibrations run dependent
      if(! fIsAfter2011) OpenInfoCalbration(run);
      fRun=run;
    }


  Int_t iC = -1;  
if (v0Centr < 80){ // analysis only for 0-80% centrality classes
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
	else if(evplaneMC > TMath::Pi()*3/2) evplaneMC-=2*TMath::Pi();
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
      
      AliAODTrack* aodTrack = event->GetTrack(iT);
      
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
	if(fV2 && fEPdet=="header")gReactionPlane=evplaneMC;
	if(fV2 && fEPdet=="V0A")gReactionPlane=fgPsi2v0aMC;
	if(fV2 && fEPdet=="V0C")gReactionPlane=fgPsi2v0cMC;
	if(fV2 && fEPdet=="TPC")gReactionPlane=fgPsi2tpcMC;

	if(fV3 && fEPdet=="V0A")gReactionPlane=fgPsi3v0aMC;
	if(fV3 && fEPdet=="V0C")gReactionPlane=fgPsi3v0cMC;
	if(fV3 && fEPdet=="TPC")gReactionPlane=fgPsi3tpcMC;
}
      else{//for data and recoMC
	if(fV2 && fEPdet=="V0A")gReactionPlane=fgPsi2v0a;
	if(fV2 && fEPdet=="V0C")gReactionPlane=fgPsi2v0c;
	if(fV2 && fEPdet=="TPC")gReactionPlane=fgPsi2tpc;

	if(fV3 && fEPdet=="V0A")gReactionPlane=fgPsi3v0a;
	if(fV3 && fEPdet=="V0C")gReactionPlane=fgPsi3v0c;
	if(fV3 && fEPdet=="TPC")gReactionPlane=fgPsi3tpc;

      }
     
 } //centrality cut condition

return gReactionPlane;
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

//----------------------------------------------------------
void AliTwoParticlePIDCorr::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
  
  
}
//------------------------------------------------------------------ 
