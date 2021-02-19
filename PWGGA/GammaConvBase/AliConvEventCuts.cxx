/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Authors: Friederike Bock, Daniel Muehlheim                              *
* A. Marin with help of Evgeny: Addition of double Gap events (Feb2019)   *
* Version 1.0                                                             *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliConvEventCuts.h"

#include <memory>
#include <TSystem.h>
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TObjString.h"
#include "AliMCEvent.h"
#include "AliAODConversionPhoton.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"
#include "AliVCaloCells.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCaloTriggerMimicHelper.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerDecisionContainer.h"



class iostream;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliConvEventCuts)
/// \endcond


const char* AliConvEventCuts::fgkCutNames[AliConvEventCuts::kNCuts] = {
  "HeavyIon",                     //0
  "CentralityMin",                //1
  "CentralityMax",                //2
  "SelectSpecialTrigger",         //3
  "SelectSpecialSubTriggerClass", //4
  "RemovePileUp",                 //5
  "RejectExtraSignals",           //6
  "VertexCut",                    //7
};


//________________________________________________________________________
AliConvEventCuts::AliConvEventCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fHeaderList(NULL),
  fDoLightOutput(0),
  fEventQuality(-1),
  fGeomEMCAL(NULL),
  fAODMCTrackArray(NULL),
  fV0Reader(NULL),
  fIsHeavyIon(0),
  fDetectorCentrality(-1),
  fModCentralityClass(0),
  fEnableVertexCut(kTRUE),
  fMaxVertexZ(10),
  fCentralityMin(0),
  fCentralityMax(0),
  fMultiplicityMethod(0),
  fSpecialTrigger(0),
  fSpecialSubTrigger(0),
  fRemovePileUp(kFALSE),
  fRemovePileUpSPD(kFALSE),
  fUseSphericity(0),
  fUseSphericityTrue(kFALSE),
  fPastFutureRejectionLow(0),
  fPastFutureRejectionHigh(0),
  fDoPileUpRejectV0MTPCout(0),
  fFPileUpRejectV0MTPCout(0),
  fRemovePileUpSDDSSDTPC(0),
  fFPileUpRejectSDDSSDTPC(0),
  fRejectExtraSignals(0),
  fOfflineTriggerMask(0),
  fHasV0AND(kTRUE),
  fIsSDDFired(kTRUE),
  fRandom(0),
  fnHeaders(0),
  fNotRejectedStart(NULL),
  fNotRejectedEnd(NULL),
  fGeneratorNames(NULL),
  fPeriodEnum(kNoPeriod),
  fEnergyEnum(kUnset),
  fTimeRangeCut(),
  fCutString(NULL),
  fCutStringRead(""),
  fUtils(NULL),
  fEtaShift(0.0),
  fDoEtaShift(kFALSE),
  fUseJetFinderForOutlier(kFALSE),
  fUseFilePathForPthard(kFALSE),
  fUseAdditionalOutlierRejection(kFALSE),
  fDoCentralityFlat(0),
  fPathWeightsFlatCent(""),
  fNameHistoNotFlatCentrality(""),
  fDoReweightHistoMCPi0(kFALSE),
  fDoReweightHistoMCEta(kFALSE),
  fDoReweightHistoMCK0s(kFALSE),
  fPathTrFReweighting(""),
  fNameHistoReweightingPi0(""),
  fNameHistoReweightingEta(""),
  fNameHistoReweightingK0s(""),
  fNameFitDataPi0(""),
  fNameFitDataEta(""),
  fNameFitDataK0s(""),
  fDoReweightHistoMCGamma(kFALSE),
  fPathTrFGammaReweighting(""),
  fNameHistoReweightingGamma(""),
  fNameDataHistoReweightingGamma(""),
  fLabelNamePileupCutTPC(""),
  fHistoEventCuts(NULL),
  fHistoPastFutureBits(NULL),
  hCentrality(NULL),
  hCentralityNotFlat(NULL),
  //hCentralityVsNumberOfPrimaryTracks(NULL),
  hVertexZ(NULL),
  hNPileupVertices(NULL),
  hPileupVertexToPrimZ(NULL),
  hPileupVertexToPrimZSPDPileup(NULL),
  hPileupVertexToPrimZTrackletvsHits(NULL),
  hEventPlaneAngle(NULL),
  fEventPlaneAngle(0),
  hTriggerClass(NULL),
  hTriggerClassSelected(NULL),
  hTriggerClassesCorrelated(NULL),
  hReweightMCHistPi0(NULL),
  hReweightMCHistEta(NULL),
  hReweightMCHistK0s(NULL),
  fFitDataPi0(NULL),
  fFitDataEta(NULL),
  fFitDataK0s(NULL),
  hReweightMCHistGamma(NULL),
  hReweightDataHistGamma(NULL),
  fAddedSignalPDGCode(0),
  fPreSelCut(kFALSE),
  fTriggerSelectedManually(kFALSE),
  fSpecialTriggerName(""),
  fSpecialSubTriggerName(""),
  fSpecialSubTriggerNameAdditional(""),
  fNSpecialSubTriggerOptions(0),
  hSPDClusterTrackletBackgroundBefore(NULL),
  hSPDClusterTrackletBackground(NULL),
  hV0MultVsNumberTPCoutTracks(NULL),
  hTPCSDDSSDClusters(NULL),
  fV0ReaderName(""),
  CaloTriggerHelperName(""),
  fCorrTaskSetting(""),
  fCaloTriggers(NULL),
  fTriggerPatchInfo(NULL),
  fMainTriggerPatchEMCAL(NULL),
  fCaloTriggersName(""),
  fCaloTriggerPatchInfoName(""),
  fTriggersEMCAL(0),
  fTriggersEMCALSelected(-1),
  fEMCALTrigInitialized(kFALSE),
  fHistoTriggThresh(NULL),
  fRunNumberTriggerOADB(-1),
  fSecProdBoundary(1.0),
  fMaxPtJetMC(0),
  fMinFacPtHard(-1),
  fMaxFacPtHard(2.5),
  fMaxFacPtHardSingleParticle(1.5),
  fMimicTrigger(kFALSE),
  fPathTriggerMimicSpecialInput(""),
  fRejectTriggerOverlap(kFALSE),
  fDoMultiplicityWeighting(kFALSE),
  fPathReweightingMult(""),
  fNameHistoReweightingMultData(""),
  fNameHistoReweightingMultMC(""),
  hReweightMultData(NULL),
  hReweightMultMC(NULL),
  fPHOSTrigger(kPHOSAny),
  fDebugLevel(0)
{
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());

  fUtils = new AliAnalysisUtils();
  //if you do not want to apply the cut on the distance between the SPD and TRK vertex:
  //fUtils->SetCutOnZVertexSPD(kFALSE);


}

//________________________________________________________________________
AliConvEventCuts::AliConvEventCuts(const AliConvEventCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fHeaderList(ref.fHeaderList),
  fDoLightOutput(ref.fDoLightOutput),
  fEventQuality(ref.fEventQuality),
  fGeomEMCAL(ref.fGeomEMCAL),
  fAODMCTrackArray(ref.fAODMCTrackArray),
  fV0Reader(NULL),
  fIsHeavyIon(ref.fIsHeavyIon),
  fDetectorCentrality(ref.fDetectorCentrality),
  fModCentralityClass(ref.fModCentralityClass),
  fEnableVertexCut(ref.fEnableVertexCut),
  fMaxVertexZ(ref.fMaxVertexZ),
  fCentralityMin(ref.fCentralityMin),
  fCentralityMax(ref.fCentralityMax),
  fMultiplicityMethod(ref.fMultiplicityMethod),
  fSpecialTrigger(ref.fSpecialTrigger),
  fSpecialSubTrigger(ref.fSpecialSubTrigger),
  fRemovePileUp(ref.fRemovePileUp),
  fRemovePileUpSPD(ref.fRemovePileUpSPD),
  fUseSphericity(ref.fUseSphericity),
  fUseSphericityTrue(ref.fUseSphericityTrue),
  fPastFutureRejectionLow(ref.fPastFutureRejectionLow),
  fPastFutureRejectionHigh(ref.fPastFutureRejectionHigh),
  fDoPileUpRejectV0MTPCout(ref.fDoPileUpRejectV0MTPCout),
  fFPileUpRejectV0MTPCout(ref.fFPileUpRejectV0MTPCout),
  fRemovePileUpSDDSSDTPC(ref.fRemovePileUpSDDSSDTPC),
  fFPileUpRejectSDDSSDTPC(ref.fFPileUpRejectSDDSSDTPC),
  fRejectExtraSignals(ref.fRejectExtraSignals),
  fOfflineTriggerMask(ref.fOfflineTriggerMask),
  fHasV0AND(ref.fHasV0AND),
  fIsSDDFired(ref.fIsSDDFired),
  fRandom(ref.fRandom),
  fnHeaders(ref.fnHeaders),
  fNotRejectedStart(NULL),
  fNotRejectedEnd(NULL),
  fGeneratorNames(ref.fGeneratorNames),
  fPeriodEnum(ref.fPeriodEnum),
  fEnergyEnum(kUnset),
  fTimeRangeCut(),
  fCutString(NULL),
  fCutStringRead(""),
  fUtils(NULL),
  fEtaShift(ref.fEtaShift),
  fDoEtaShift(ref.fDoEtaShift),
  fUseJetFinderForOutlier(ref.fUseJetFinderForOutlier),
  fUseFilePathForPthard(ref.fUseFilePathForPthard),
  fUseAdditionalOutlierRejection(ref.fUseAdditionalOutlierRejection),
  fDoCentralityFlat(ref.fDoCentralityFlat),
  fPathWeightsFlatCent(ref.fPathWeightsFlatCent),
  fNameHistoNotFlatCentrality(ref.fNameHistoNotFlatCentrality),
  fDoReweightHistoMCPi0(ref.fDoReweightHistoMCPi0),
  fDoReweightHistoMCEta(ref.fDoReweightHistoMCEta),
  fDoReweightHistoMCK0s(ref.fDoReweightHistoMCK0s),
  fPathTrFReweighting(ref.fPathTrFReweighting),
  fNameHistoReweightingPi0(ref.fNameHistoReweightingPi0),
  fNameHistoReweightingEta(ref.fNameHistoReweightingEta),
  fNameHistoReweightingK0s(ref.fNameHistoReweightingK0s),
  fNameFitDataPi0(ref.fNameFitDataPi0),
  fNameFitDataEta(ref.fNameFitDataEta),
  fNameFitDataK0s(ref.fNameFitDataK0s),
  fDoReweightHistoMCGamma(ref.fDoReweightHistoMCGamma),
  fPathTrFGammaReweighting(ref.fPathTrFGammaReweighting),
  fNameHistoReweightingGamma(ref.fNameHistoReweightingGamma),
  fNameDataHistoReweightingGamma(ref.fNameDataHistoReweightingGamma),
  fLabelNamePileupCutTPC(ref.fLabelNamePileupCutTPC),
  fHistoEventCuts(NULL),
  fHistoPastFutureBits(NULL),
  hCentrality(ref.hCentrality),
  hCentralityNotFlat(ref.hCentralityNotFlat),
  //hCentralityVsNumberOfPrimaryTracks(ref.hCentralityVsNumberOfPrimaryTracks),
  hVertexZ(ref.hVertexZ),
  hNPileupVertices(ref.hNPileupVertices),
  hPileupVertexToPrimZ(ref.hPileupVertexToPrimZ),
  hPileupVertexToPrimZSPDPileup(ref.hPileupVertexToPrimZSPDPileup),
  hPileupVertexToPrimZTrackletvsHits(ref.hPileupVertexToPrimZTrackletvsHits),
  hEventPlaneAngle(ref.hEventPlaneAngle),
  fEventPlaneAngle(ref.fEventPlaneAngle),
  hTriggerClass(NULL),
  hTriggerClassSelected(NULL),
  hTriggerClassesCorrelated(NULL),
  hReweightMCHistPi0(ref.hReweightMCHistPi0),
  hReweightMCHistEta(ref.hReweightMCHistEta),
  hReweightMCHistK0s(ref.hReweightMCHistK0s),
  fFitDataPi0(ref.fFitDataPi0),
  fFitDataEta(ref.fFitDataEta),
  fFitDataK0s(ref.fFitDataK0s),
  hReweightMCHistGamma(ref.hReweightMCHistGamma),
  hReweightDataHistGamma(ref.hReweightDataHistGamma),
  fAddedSignalPDGCode(ref.fAddedSignalPDGCode),
  fPreSelCut(ref.fPreSelCut),
  fTriggerSelectedManually(ref.fTriggerSelectedManually),
  fSpecialTriggerName(ref.fSpecialTriggerName),
  fSpecialSubTriggerName(ref.fSpecialSubTriggerName),
  fSpecialSubTriggerNameAdditional(ref.fSpecialSubTriggerNameAdditional),
  fNSpecialSubTriggerOptions(ref.fNSpecialSubTriggerOptions),
  hSPDClusterTrackletBackgroundBefore(NULL),
  hSPDClusterTrackletBackground(NULL),
  hV0MultVsNumberTPCoutTracks(NULL),
  hTPCSDDSSDClusters(NULL),
  fV0ReaderName(ref.fV0ReaderName),
  CaloTriggerHelperName(ref.CaloTriggerHelperName),
  fCorrTaskSetting(ref.fCorrTaskSetting),
  fCaloTriggers(NULL),
  fTriggerPatchInfo(NULL),
  fMainTriggerPatchEMCAL(NULL),
  fCaloTriggersName(ref.fCaloTriggersName),
  fCaloTriggerPatchInfoName(ref.fCaloTriggerPatchInfoName),
  fTriggersEMCAL(ref.fTriggersEMCAL),
  fTriggersEMCALSelected(ref.fTriggersEMCALSelected),
  fEMCALTrigInitialized(kFALSE),
  fHistoTriggThresh(ref.fHistoTriggThresh),
  fRunNumberTriggerOADB(ref.fRunNumberTriggerOADB),
  fSecProdBoundary(ref.fSecProdBoundary),
  fMaxPtJetMC(ref.fMaxPtJetMC),
  fMinFacPtHard(ref.fMinFacPtHard),
  fMaxFacPtHard(ref.fMaxFacPtHard),
  fMaxFacPtHardSingleParticle(ref.fMaxFacPtHardSingleParticle),
  fMimicTrigger(ref.fMimicTrigger),
  fPathTriggerMimicSpecialInput(ref.fPathTriggerMimicSpecialInput),
  fRejectTriggerOverlap(ref.fRejectTriggerOverlap),
  fDoMultiplicityWeighting(ref.fDoMultiplicityWeighting),
  fPathReweightingMult(ref.fPathReweightingMult),
  fNameHistoReweightingMultData(ref.fNameHistoReweightingMultData),
  fNameHistoReweightingMultMC(ref.fNameHistoReweightingMultMC),
  hReweightMultData(ref.hReweightMultData),
  hReweightMultMC(ref.hReweightMultMC),
  fPHOSTrigger(kPHOSAny),
  fDebugLevel(ref.fDebugLevel)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());
  fUtils = new AliAnalysisUtils();
  // dont copy histograms (if you like histograms, call InitCutHistograms())


}


//________________________________________________________________________
AliConvEventCuts::~AliConvEventCuts() {
  // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  //    delete fHistograms;
  // fHistograms = NULL;
  if(fCutString != NULL){
      delete fCutString;
      fCutString = NULL;
  }
  if(fNotRejectedStart){
      delete[] fNotRejectedStart;
      fNotRejectedStart = NULL;
  }
  if(fNotRejectedEnd){
      delete[] fNotRejectedEnd;
      fNotRejectedEnd = NULL;
  }
  if(fGeneratorNames){
      delete[] fGeneratorNames;
      fGeneratorNames = NULL;
  }
  if(fUtils){
    delete fUtils;
    fUtils = NULL;
  }
  if(fAODMCTrackArray){
    delete[] fAODMCTrackArray;
    fAODMCTrackArray = 0x0;
  }

}

//________________________________________________________________________
void AliConvEventCuts::InitCutHistograms(TString name, Bool_t preCut){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }
  if(fDoLightOutput==2) {
      AliInfo("Minimal output chosen");
      return;
  }

  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("ConvEventCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  if (hReweightMCHistPi0){
    hReweightMCHistPi0->SetName("MCInputForWeightingPi0");
    fHistograms->Add(hReweightMCHistPi0);
  }
  if (hReweightMCHistEta){
    hReweightMCHistEta->SetName("MCInputForWeightingEta");
    fHistograms->Add(hReweightMCHistEta);
  }
  if (hReweightMCHistK0s){
    hReweightMCHistK0s->SetName("MCInputForWeightingK0s");
    fHistograms->Add(hReweightMCHistK0s);
  }

  if (hReweightMCHistGamma){
    hReweightMCHistGamma->SetName("MCInputForWeightingGamma");
    fHistograms->Add(hReweightMCHistGamma);
  }


  if (hReweightMultData){
    hReweightMultData->SetName(Form("hReweightMultData_%s",GetCutNumber().Data()));
    fHistograms->Add(hReweightMultData);
  }
  if (hReweightMultMC){
    hReweightMultMC->SetName(Form("hReweightMultMC_%s",GetCutNumber().Data()));
    fHistograms->Add(hReweightMultMC);
  }

  if(!fDoLightOutput){
    if (fIsHeavyIon == 1){
      hSPDClusterTrackletBackgroundBefore = new TH2F(Form("SPD tracklets vs SPD clusters %s before Pileup Cut",GetCutNumber().Data()),"SPD tracklets vs SPD clusters", 200, 0, 6000, 200, 0, 20000);
      fHistograms->Add(hSPDClusterTrackletBackgroundBefore);
      hSPDClusterTrackletBackground = new TH2F(Form("SPD tracklets vs SPD clusters %s",GetCutNumber().Data()),"SPD tracklets vs SPD clusters", 200, 0, 6000, 200, 0, 20000);
      fHistograms->Add(hSPDClusterTrackletBackground);
    } else{
      hSPDClusterTrackletBackgroundBefore = new TH2F(Form("SPD tracklets vs SPD clusters %s before Pileup Cut",GetCutNumber().Data()),"SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
      fHistograms->Add(hSPDClusterTrackletBackgroundBefore);
      hSPDClusterTrackletBackground = new TH2F(Form("SPD tracklets vs SPD clusters %s",GetCutNumber().Data()),"SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
      fHistograms->Add(hSPDClusterTrackletBackground);
    }

    if (fRemovePileUpSDDSSDTPC){
      hTPCSDDSSDClusters = new TH2F(Form("SDD+SSD clusters vs TPC clusters %s",GetCutNumber().Data()),"SDD+SSD clusters vs TPC clusters", 500, 0., 6.e+6, 500, 0., 5.e+4);
      fHistograms->Add(hTPCSDDSSDClusters);
    }

    if (fDoPileUpRejectV0MTPCout){
      if(fIsHeavyIon == 1)
	    hV0MultVsNumberTPCoutTracks    = new TH2F("V0Mult vs TPCout Tracks", "V0Mult vs TPCout Tracks", 500, 0, 15000, 500, 0, 40000);
      else if(fIsHeavyIon == 2)
	    hV0MultVsNumberTPCoutTracks    = new TH2F("V0Mult vs TPCout Tracks", "V0Mult vs TPCout Tracks", 500, 0, 1000, 500, 0, 2500);
      else
	    hV0MultVsNumberTPCoutTracks    = new TH2F("V0Mult vs TPCout Tracks", "V0Mult vs TPCout Tracks", 200, 0, 400, 500, 0, 1500);
      fHistograms->Add(hV0MultVsNumberTPCoutTracks);
    }
  }

  // if(fIsHeavyIon > 0){ // commented as mult. dep. analyses in pp started
  if( fModCentralityClass == 20){ // high mult 0.1%
    const Int_t centBins = 145;
    Double_t arrCent[centBins + 1];
    for(Int_t i = 0; i < centBins + 1; i++){
      if(i < 50) arrCent[i] = i*0.1;
      else if( i < centBins) arrCent[i] = 5 + (i-50);
      else arrCent[i] = 100;
    }
    hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",centBins,arrCent);
  } else if( fModCentralityClass == 21){// high mult 0.01%
    const Int_t centBins = 235;
    Double_t arrCent[centBins + 1];
    for(Int_t i = 0; i < centBins + 1; i++){
      if(i < 100) arrCent[i] = i*0.01;
      else if(i < 140) arrCent[i] = 1 + 0.1*(i-100);
      else if( i < centBins) arrCent[i] = 5 + (i-140);
      else arrCent[i] = 100;
    }
    hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",centBins,arrCent);
  } else {
    hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",210,0,105);
  }
    fHistograms->Add(hCentrality);
  // }

  //hCentralityVsNumberOfPrimaryTracks=new TH2F(Form("Centrality vs Primary Tracks %s",GetCutNumber().Data()),"Centrality vs Primary Tracks ",400,0,100,4000,0,4000);
  //fHistograms->Add(hCentralityVsNumberOfPrimaryTracks); commented on 3.3.2015 because it's in the main Task

  hVertexZ              = new TH1F(Form("VertexZ %s",GetCutNumber().Data()),"VertexZ",1000,-50,50);
  fHistograms->Add(hVertexZ);

  hNPileupVertices      = new TH1F(Form("NPileupVertices %s",GetCutNumber().Data()),"NPileupVertices",30,-0.5,29.5);
  fHistograms->Add(hNPileupVertices);

  hPileupVertexToPrimZ  = new TH1F(Form("PileupVertexDistance %s",GetCutNumber().Data()),"PileupVertexDistance",600,-15,15);
  fHistograms->Add(hPileupVertexToPrimZ);
  hPileupVertexToPrimZSPDPileup  = new TH1F(Form("PileupVertexDistance_SPDPileup %s",GetCutNumber().Data()),"PileupVertexDistance_SPDPileup",600,-15,15);
  fHistograms->Add(hPileupVertexToPrimZSPDPileup);
  hPileupVertexToPrimZTrackletvsHits  = new TH1F(Form("PileupVertexDistance_TrackletvsHits %s",GetCutNumber().Data()),"PileupVertexDistance_TrackletvsHits",600,-15,15);
  fHistograms->Add(hPileupVertexToPrimZTrackletvsHits);

  if(fIsHeavyIon == 1){
    hEventPlaneAngle = new TH1F(Form("EventPlaneAngle %s",GetCutNumber().Data()),"EventPlaneAngle",60, 0, TMath::Pi());
    fHistograms->Add(hEventPlaneAngle);
  }
  fHistoPastFutureBits=new TH1F(Form("PastFutureBits %s",GetCutNumber().Data()),"Past Future Bits",180,-90*25,90*25);
  fHistograms->Add(fHistoPastFutureBits);

  // Event Cuts and Info
  if(preCut){
    fHistoEventCuts=new TH1F(Form("ESD_EventCuts %s",GetCutNumber().Data()),"Event Cuts",8,-0.5,7.5);
    fHistoEventCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoEventCuts->GetXaxis()->SetBinLabel(2,"OfflineTrigger");
    fHistoEventCuts->GetXaxis()->SetBinLabel(3,"nvtxcontr");
    fHistoEventCuts->GetXaxis()->SetBinLabel(4,"VertexZ");
    fHistoEventCuts->GetXaxis()->SetBinLabel(5,"pileup");
    fHistoEventCuts->GetXaxis()->SetBinLabel(6,"centrsel");
    fHistoEventCuts->GetXaxis()->SetBinLabel(7,"OOB-pileup");
    fHistoEventCuts->GetXaxis()->SetBinLabel(8,"out");
    fHistograms->Add(fHistoEventCuts);

    hTriggerClass= new TH1F(Form("OfflineTrigger %s",GetCutNumber().Data()),"OfflineTrigger",37,-0.5,36.5);
    hTriggerClass->GetXaxis()->SetBinLabel( 1,"kMB");
    hTriggerClass->GetXaxis()->SetBinLabel( 2,"kINT7");
    hTriggerClass->GetXaxis()->SetBinLabel( 3,"kMUON");
    hTriggerClass->GetXaxis()->SetBinLabel( 4,"kHighMult");
    hTriggerClass->GetXaxis()->SetBinLabel( 5,"kEMC1");
    hTriggerClass->GetXaxis()->SetBinLabel( 6,"kCINT5");
    hTriggerClass->GetXaxis()->SetBinLabel( 7,"kCMUS5/kMUSPB");
    hTriggerClass->GetXaxis()->SetBinLabel( 8,"kMUSH7/kMUSHPB");
    hTriggerClass->GetXaxis()->SetBinLabel( 9,"kMUL7/kMuonLikePB");
    hTriggerClass->GetXaxis()->SetBinLabel(10,"kMUU7/kMuonUnlikePB");
    hTriggerClass->GetXaxis()->SetBinLabel(11,"kEMC7/kEMC8");
    hTriggerClass->GetXaxis()->SetBinLabel(12,"kMUS7");
    hTriggerClass->GetXaxis()->SetBinLabel(13,"kPHI1");
    hTriggerClass->GetXaxis()->SetBinLabel(14,"kPHI7/kPHI8/kPHOSPb");
    hTriggerClass->GetXaxis()->SetBinLabel(15,"kEMCEJE");
    hTriggerClass->GetXaxis()->SetBinLabel(16,"kEMCEGA");
    hTriggerClass->GetXaxis()->SetBinLabel(17,"kCentral");
    hTriggerClass->GetXaxis()->SetBinLabel(18,"kSemiCentral");
    hTriggerClass->GetXaxis()->SetBinLabel(19,"kDG5");
    hTriggerClass->GetXaxis()->SetBinLabel(20,"kZED");
    hTriggerClass->GetXaxis()->SetBinLabel(21,"kSPI7/kSPI");
    hTriggerClass->GetXaxis()->SetBinLabel(22,"kINT8");
    hTriggerClass->GetXaxis()->SetBinLabel(23,"kMuonSingleLowPt8");
    hTriggerClass->GetXaxis()->SetBinLabel(24,"kMuonSingleHighPt8");
    hTriggerClass->GetXaxis()->SetBinLabel(25,"kMuonLikeLowPt8");
    hTriggerClass->GetXaxis()->SetBinLabel(26,"kMuonUnlikeLowPt8");
    hTriggerClass->GetXaxis()->SetBinLabel(27,"kMuonUnlikeLowPt0");
    hTriggerClass->GetXaxis()->SetBinLabel(28,"kUserDefined");
    hTriggerClass->GetXaxis()->SetBinLabel(29,"kTRD");
    hTriggerClass->GetXaxis()->SetBinLabel(30,"kFastOnly");
    hTriggerClass->GetXaxis()->SetBinLabel(31,"kAnyINT");
    hTriggerClass->GetXaxis()->SetBinLabel(32,"kAny");
    hTriggerClass->GetXaxis()->SetBinLabel(33,"V0AND");
    hTriggerClass->GetXaxis()->SetBinLabel(34,"NOT kFastOnly");
    hTriggerClass->GetXaxis()->SetBinLabel(35,"kCaloOnly");
    hTriggerClass->GetXaxis()->SetBinLabel(36,"failed Physics Selection");
    hTriggerClass->GetXaxis()->SetBinLabel(37,"mimickedTrigger");
    fHistograms->Add(hTriggerClass);
  }
  if(!preCut){
    hTriggerClassSelected= new TH1F(Form("OfflineTriggerSelected %s",GetCutNumber().Data()),"OfflineTriggerSelected",35,-0.5,34.5);
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 1,"kMB");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 2,"kINT7");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 3,"kMUON");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 4,"kHighMult");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 5,"kEMC1");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 6,"kCINT5");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 7,"kCMUS5/kMUSPB");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 8,"kMUSH7/kMUSHPB");
    hTriggerClassSelected->GetXaxis()->SetBinLabel( 9,"kMUL7/kMuonLikePB");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(10,"kMUU7/kMuonUnlikePB");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(11,"kEMC7/kEMC8");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(12,"kMUS7");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(13,"kPHI1");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(14,"kPHI7/kPHI8/kPHOSPb");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(15,"kEMCEJE");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(16,"kEMCEGA");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(17,"kCentral");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(18,"kSemiCentral");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(19,"kDG5");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(20,"kZED");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(21,"kSPI7/kSPI");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(22,"kINT8");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(23,"kMuonSingleLowPt8");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(24,"kMuonSingleHighPt8");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(25,"kMuonLikeLowPt8");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(26,"kMuonUnlikeLowPt8");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(27,"kMuonUnlikeLowPt0");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(28,"kUserDefined");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(29,"kTRD");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(30,"kFastOnly");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(31,"kAnyINT");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(32,"kAny");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(33,"V0AND");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(34,"NOT kFastOnly");
    hTriggerClassSelected->GetXaxis()->SetBinLabel(35,"mimickedTrigger");
    fHistograms->Add(hTriggerClassSelected);

    if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9 || fSpecialTrigger == 10){
      hTriggerClassesCorrelated= new TH1F(Form("TriggerCorrelations %s",GetCutNumber().Data()),"Triggers Correlated with EMCal triggers",17,-0.5,16.5);
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 1,"kMB");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 2,"kINT7");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 3,"kEMC1");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 4,"kEMC7");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 5,"kEMCEJE");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 6,"kEMCEJ1");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 7,"kEMCEJ2");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 8,"kEMCEGA");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 9,"kEMCEG1");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 10,"kEMCEG2");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 11,"kDMC7");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 12,"kDMCDJE");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 13,"kDMCDJ1");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 14,"kDMCDJ2");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 15,"kDMCDGA");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 16,"kDMCDG1");
      hTriggerClassesCorrelated->GetXaxis()->SetBinLabel( 17,"kDMCDG2");
      fHistograms->Add(hTriggerClassesCorrelated);
    }

  }
  TH1::AddDirectory(kTRUE);
}

///________________________________________________________________________
Bool_t AliConvEventCuts::EventIsSelected(AliVEvent *event, AliMCEvent *mcEvent){
  // Process Event Selection

  Int_t cutindex=0;
  if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
  cutindex++;

  // Check for MC event
  Bool_t isMC = kFALSE;
  if(mcEvent && event->IsA()==AliESDEvent::Class()){
    // Check if MC event is correctly loaded
    AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcHandler){
      fEventQuality = 2;
      return kFALSE;
    }
    if (!mcHandler->InitOk() ){
      fEventQuality = 2;
      return kFALSE;
    }
    if (!mcHandler->TreeK() ){
      fEventQuality = 2;
      return kFALSE;
    }
    // TrackRefs.root is currently not used by anyone
    // and was excluded from future MC productions
    // if (!mcHandler->TreeTR() ) {
    //   fEventQuality = 2;
    //   return kFALSE;
    // }
    isMC = kTRUE;
  }



  // Event Trigger
  //    cout << "before event trigger" << endl;
  if(!IsTriggerSelected(event, isMC )){
    if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
    fEventQuality = 3;
    return kFALSE;
  }
  cutindex++;

  if(event->IsA()==AliESDEvent::Class()){
    AliTriggerAnalysis fTriggerAnalysis;// = new AliTriggerAnalysis;
    fHasV0AND = fTriggerAnalysis.IsOfflineTriggerFired((AliESDEvent*)event, AliTriggerAnalysis::kV0AND);
    if(fHasV0AND&&hTriggerClass)hTriggerClass->Fill(32);
  }
  //   cout << "event number " << ((AliESDEvent*)event)->GetEventNumberInFile() << " entered"<< endl;

  // Number of Contributors Cut
  if (fEnableVertexCut){
    if(GetNumberOfContributorsVtx(event)<=0) {
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      fEventQuality = 5;
      return kFALSE;
    }
  }
  cutindex++;

  // Z Vertex Position Cut
  if (fEnableVertexCut){
    if(!VertexZCut(event)){
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      fEventQuality = 4;
      return kFALSE;
    }
  }
  cutindex++;

  // SPD clusters vs tracklets to check for pileup/background
  Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
  Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
  if(hSPDClusterTrackletBackgroundBefore) hSPDClusterTrackletBackgroundBefore->Fill(nTracklets, (nClustersLayer0 + nClustersLayer1));


  Double_t distZMax     = 0;
  if(event->IsA()==AliESDEvent::Class()){
    Int_t nPileVert = ((AliESDEvent*)event)->GetNumberOfPileupVerticesSPD();
    if (hNPileupVertices) hNPileupVertices->Fill(nPileVert);
    if (nPileVert > 0){
      for(Int_t i=0; i<nPileVert;i++){
        const AliESDVertex* pv  = ((AliESDEvent*)event)->GetPileupVertexSPD(i);
        Int_t nc2               = pv->GetNContributors();
        if(nc2>=3){
          Double_t z1     = ((AliESDEvent*)event)->GetPrimaryVertexSPD()->GetZ();
          Double_t z2     = pv->GetZ();
          Double_t distZ  = z2-z1;
          if (TMath::Abs(distZMax) <  TMath::Abs(distZ) ){
            distZMax      = distZ;
          }
        }
      }
      if (hPileupVertexToPrimZ) hPileupVertexToPrimZ->Fill(distZMax);
    }
  }

  // Pile Up Rejection
  if (fIsHeavyIon == 2){
    if(GetUseNewMultiplicityFramework()){// for Run2 pPb
      if(fUtils->IsPileUpMV(event)){
        if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
        fEventQuality = 6;
        return kFALSE;
      }
    } else{
      if(fUtils->IsFirstEventInChunk(event)){
        if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
        fEventQuality = 6;
        return kFALSE;
      }
      if(fRemovePileUpSPD){
        if(fUtils->IsPileUpEvent(event)){
          if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
          if (hPileupVertexToPrimZSPDPileup) hPileupVertexToPrimZSPDPileup->Fill(distZMax);
          fEventQuality = 6;
          return kFALSE;
        }
        if (fUtils->IsSPDClusterVsTrackletBG(event)){
          if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
          if (hPileupVertexToPrimZTrackletvsHits) hPileupVertexToPrimZTrackletvsHits->Fill(distZMax);
          fEventQuality = 11;
          return kFALSE;
        }
      }
    }
  } else if(fRemovePileUpSPD){
    if(event->IsPileupFromSPD(3,0.8,3.,2.,5.) ){
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      if (hPileupVertexToPrimZSPDPileup) hPileupVertexToPrimZSPDPileup->Fill(distZMax);
      fEventQuality = 6;
      return kFALSE;
    }
    if (fUtils->IsSPDClusterVsTrackletBG(event)){
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      if (hPileupVertexToPrimZTrackletvsHits) hPileupVertexToPrimZTrackletvsHits->Fill(distZMax);
      fEventQuality = 11;
      return kFALSE;
    }
  }
  cutindex++;

  // Centrality Selection
  if(!IsCentralitySelected(event,mcEvent)){
    if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
    fEventQuality = 1;
    return kFALSE;
  }
  cutindex++;

  if(fRemovePileUp && IsOutOfBunchPileupPastFuture(event)){
    if(fHistoEventCuts) fHistoEventCuts->Fill(cutindex);
    fEventQuality = 12;
    return kFALSE;
  }
  cutindex++;
  // Fill Event Histograms
  if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
  if(hCentrality)hCentrality->Fill(GetCentrality(event));
  if(hVertexZ)hVertexZ->Fill(event->GetPrimaryVertex()->GetZ());
  //   if(hCentralityVsNumberOfPrimaryTracks)
  //      hCentralityVsNumberOfPrimaryTracks->Fill(GetCentrality(event),
  //                                               ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()
  //                                                ->GetTask(fV0ReaderName.Data()))->GetNumberOfPrimaryTracks());

  if(fIsHeavyIon == 1){
    AliEventplane *EventPlane = event->GetEventplane();
    fEventPlaneAngle = EventPlane->GetEventplane("V0",event,2);
    if(hEventPlaneAngle)hEventPlaneAngle->Fill(TMath::Abs(fEventPlaneAngle));
  }
  if(hSPDClusterTrackletBackground) hSPDClusterTrackletBackground->Fill(nTracklets, (nClustersLayer0 + nClustersLayer1));

  // for data from LHC18r apply timeRange cut
  if((fPeriodEnum==kLHC18qr) && (!isMC)) {

    // no need to check here for new runNumber, InitFromRunNumber does this internally
    fTimeRangeCut.InitFromRunNumber(event->GetRunNumber());
    if(fTimeRangeCut.CutEvent(event)){
      // since timecut is necessary because of TPC problems
      fEventQuality = 9;
      return kFALSE;
    }
  }

  fEventQuality = 0;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConvEventCuts::UpdateCutString() {
  ///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
void AliConvEventCuts::LoadWeightingFlatCentralityFromFile() {

  AliInfo("Entering loading of weights for centrality flattening");
  TFile *w = TFile::Open(fPathWeightsFlatCent.Data());
  if(!w){
    AliError(Form("file for centrality flattening %s not found",fPathWeightsFlatCent.Data()));
    return;
  }

  if (fNameHistoNotFlatCentrality.CompareTo("") != 0 && (fDoCentralityFlat > 0)){
    cout << "I have to find: " <<  fNameHistoNotFlatCentrality.Data() << endl;
    TH1D *hCentralityNotFlattemp = (TH1D*)w->Get(fNameHistoNotFlatCentrality.Data());
    hCentralityNotFlat = new TH1D(*hCentralityNotFlattemp);
    if (hCentralityNotFlat) AliInfo(Form("%s has been loaded from %s", fNameHistoNotFlatCentrality.Data(),fPathWeightsFlatCent.Data() ));
    else AliWarning(Form("%s not found in %s", fNameHistoNotFlatCentrality.Data() ,fPathWeightsFlatCent.Data()));
    hCentralityNotFlat->SetDirectory(0);
  }

  w->Close();
  delete w;
}

///________________________________________________________________________
void AliConvEventCuts::LoadWeightingMultiplicityFromFile() {

  AliInfo("Entering loading of weights for multiplicity weighting");
  TFile *w = TFile::Open(fPathReweightingMult.Data());
  if(!w){
    AliError(Form("file for multiplicity reweighting %s not found",fPathReweightingMult.Data()));
    return;
  }

  if (fNameHistoReweightingMultData.CompareTo("") != 0 && (fDoMultiplicityWeighting > 0)){
    cout << "I have to find: " <<  fNameHistoReweightingMultData.Data() << endl;
    TH1D *hReweightMultDatatemp = (TH1D*)w->Get(fNameHistoReweightingMultData.Data());
    if(hReweightMultDatatemp){
      hReweightMultData = new TH1D(*hReweightMultDatatemp);
      hReweightMultData->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingMultData.Data(),fPathReweightingMult.Data() ));
    } else  AliError(Form("%s was not contained in %s", fNameHistoReweightingMultData.Data(),fPathReweightingMult.Data() ));
  }
  if (fNameHistoReweightingMultMC.CompareTo("") != 0 && (fDoMultiplicityWeighting > 0)){
    cout << "I have to find: " <<  fNameHistoReweightingMultMC.Data() << endl;
    TH1D *hReweightMultMCtemp = (TH1D*)w->Get(fNameHistoReweightingMultMC.Data());
    if(hReweightMultMCtemp){
      hReweightMultMC = new TH1D(*hReweightMultMCtemp);
      hReweightMultMC->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingMultMC.Data(),fPathReweightingMult.Data() ));
    } else  AliError(Form("%s was not contained in %s", fNameHistoReweightingMultMC.Data(),fPathReweightingMult.Data() ));
  }

  w->Close();
  delete w;
}


///________________________________________________________________________
void AliConvEventCuts::LoadReweightingHistosMCFromFile() {

  AliInfo("Entering loading of histograms for weighting");
  TFile *f = TFile::Open(fPathTrFReweighting.Data());
  if(!f){
    AliError(Form("file for weighting %s not found",fPathTrFReweighting.Data()));
    return;
  }
  if (fNameHistoReweightingPi0.CompareTo("") != 0 && fDoReweightHistoMCPi0 ){
    cout << "I have to find: " <<  fNameHistoReweightingPi0.Data() << endl;
    TH1D *hReweightMCHistPi0temp = (TH1D*)f->Get(fNameHistoReweightingPi0.Data());
    if(hReweightMCHistPi0temp){
      hReweightMCHistPi0 = new TH1D(*hReweightMCHistPi0temp);
      hReweightMCHistPi0->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingPi0.Data(),fPathTrFReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s", fNameHistoReweightingPi0.Data() ,fPathTrFReweighting.Data()));
  }
  if (fNameFitDataPi0.CompareTo("") != 0 && fDoReweightHistoMCPi0 ){
    cout << "I have to find: " <<  fNameFitDataPi0.Data() << endl;
    TF1 *fFitDataPi0temp = (TF1*)f->Get(fNameFitDataPi0.Data());
    if(fFitDataPi0temp){
      fFitDataPi0 = new TF1(*fFitDataPi0temp);
      AliInfo(Form("%s has been loaded from %s", fNameFitDataPi0.Data(),fPathTrFReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s",fPathTrFReweighting.Data(), fNameFitDataPi0.Data() ));
  }

  if (fNameHistoReweightingEta.CompareTo("") != 0 && fDoReweightHistoMCEta){
    cout << "I have to find: " <<  fNameHistoReweightingEta.Data() << endl;
    TH1D *hReweightMCHistEtatemp = (TH1D*)f->Get(fNameHistoReweightingEta.Data());
    if(hReweightMCHistEtatemp){
      hReweightMCHistEta = new TH1D(*hReweightMCHistEtatemp);
      hReweightMCHistEta->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingEta.Data(),fPathTrFReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s", fNameHistoReweightingEta.Data(),fPathTrFReweighting.Data() ));
  }

  if (fNameFitDataEta.CompareTo("") != 0 && fDoReweightHistoMCEta){
    cout << "I have to find: " <<  fNameFitDataEta.Data() << endl;
    TF1 *fFitDataEtatemp = (TF1*)f->Get(fNameFitDataEta.Data());
    if(fFitDataEtatemp){
      fFitDataEta = new TF1(*fFitDataEtatemp);
      AliInfo(Form("%s has been loaded from %s", fNameFitDataEta.Data(),fPathTrFReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s", fNameFitDataEta.Data(),fPathTrFReweighting.Data() ));

  }
  if (fNameHistoReweightingK0s.CompareTo("") != 0 && fDoReweightHistoMCK0s){
    cout << "I have to find: " <<  fNameHistoReweightingK0s.Data() << endl;
    TH1D *hReweightMCHistK0stemp = (TH1D*)f->Get(fNameHistoReweightingK0s.Data());
    if(hReweightMCHistK0stemp){
      hReweightMCHistK0s = new TH1D(*hReweightMCHistK0stemp);
      hReweightMCHistK0s->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingK0s.Data(),fPathTrFReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s", fNameHistoReweightingK0s.Data(),fPathTrFReweighting.Data() ));
  }

  if (fNameFitDataK0s.CompareTo("") != 0 && fDoReweightHistoMCK0s){
    cout << "I have to find: " <<  fNameFitDataK0s.Data() << endl;
    TF1 *fFitDataK0stemp = (TF1*)f->Get(fNameFitDataK0s.Data());
    if(fFitDataK0stemp){
      fFitDataK0s = new TF1(*fFitDataK0stemp);
      AliInfo(Form("%s has been loaded from %s", fNameFitDataK0s.Data(),fPathTrFReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s", fNameFitDataK0s.Data(),fPathTrFReweighting.Data() ));
  }
  f->Close();
  delete f;
}

///________________________________________________________________________
void AliConvEventCuts::LoadGammaPtReweightingHistosMCFromFile() {

  AliInfo("Entering loading of histograms for gamma pT weighting");
  TFile *f = TFile::Open(fPathTrFGammaReweighting.Data());
  if(!f){
    AliError(Form("file for gamma pT  weighting %s not found",fPathTrFGammaReweighting.Data()));
    return;
  }
  if (fNameHistoReweightingGamma.CompareTo("") != 0 && fDoReweightHistoMCGamma ){
    cout << "I have to find: " <<  fNameHistoReweightingGamma.Data() << endl;
    TH1D *hReweightMCHistGammatemp = (TH1D*)f->Get(fNameHistoReweightingGamma.Data());
    if(hReweightMCHistGammatemp){
      hReweightMCHistGamma = new TH1D(*hReweightMCHistGammatemp);
      hReweightMCHistGamma->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingGamma.Data(),fPathTrFGammaReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s", fNameHistoReweightingGamma.Data() ,fPathTrFGammaReweighting.Data()));
  }
  if (fNameDataHistoReweightingGamma.CompareTo("") != 0 && fDoReweightHistoMCGamma ){
    cout << "I have to find: " <<  fNameDataHistoReweightingGamma.Data() << endl;
    TH1D *hReweightDataHistGammatemp = (TH1D*)f->Get(fNameDataHistoReweightingGamma.Data());
    if(hReweightDataHistGammatemp){
      hReweightDataHistGamma = new TH1D(*hReweightDataHistGammatemp);
      hReweightDataHistGamma->SetDirectory(0);
      AliInfo(Form("%s has been loaded from %s", fNameDataHistoReweightingGamma.Data(),fPathTrFGammaReweighting.Data() ));
    } else AliWarning(Form("%s not found in %s",fPathTrFGammaReweighting.Data(), fNameDataHistoReweightingGamma.Data() ));
  }


  f->Close();
  delete f;
}


///________________________________________________________________________
Bool_t AliConvEventCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());

  // Initialize Cuts from a given Cut string
  if(fDoCentralityFlat > 0){
    AliInfo("Centrality flattening was enabled");
    LoadWeightingFlatCentralityFromFile();
  }

  if (fDoMultiplicityWeighting){
    AliInfo("Multiplicity weighting was enabled");
    LoadWeightingMultiplicityFromFile();
  }

  if(fDoReweightHistoMCPi0 || fDoReweightHistoMCEta || fDoReweightHistoMCK0s) {
    AliInfo("Particle Weighting was enabled");
    LoadReweightingHistosMCFromFile();
  }

  if(fDoReweightHistoMCGamma) {
    AliInfo("Gamma pT Weighting was enabled");
    LoadGammaPtReweightingHistosMCFromFile();
  }


  AliInfo(Form("Set Event Cut Number: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
    AliError("Cut selection is not alphanumeric");
    return kFALSE;
  }

  if (fV0ReaderName.CompareTo("") == 0){
    fV0ReaderName = "V0ReaderV1";
  }
  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  PrintCutsWithValues();

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConvEventCuts::SetCut(cutIds cutID, const Int_t value) {
  // /Set individual cut ID
  // "HeavyIon",                     //0
  // "CentralityMin",                //1
  // "CentralityMax",                //2
  // "SelectSpecialTrigger",         //3
  // "SelectSpecialSubTriggerClass", //4
  // "RemovePileUp",                 //5
  // "RejectExtraSignals",           //6
  // "VertexCut",                    //7

  switch (cutID) {
  case kisHeavyIon:
    if( SetIsHeavyIon(value)) {
      fCuts[kisHeavyIon] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kCentralityMin:
    if( SetCentralityMin(value)) {
      fCuts[kCentralityMin] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kCentralityMax:
    if( SetCentralityMax(value)) {
      fCuts[kCentralityMax] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kSelectSpecialTriggerAlias:
    if( SetSelectSpecialTrigger(value)) {
      fCuts[kSelectSpecialTriggerAlias] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kSelectSubTriggerClass:
    if( SetSelectSubTriggerClass(value)) {
      fCuts[kSelectSubTriggerClass] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kremovePileUp:
    if( SetRemovePileUp(value)) {
      fCuts[kremovePileUp] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kExtraSignals:
    if( SetRejectExtraSignalsCut(value)) {
      fCuts[kExtraSignals] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kVertex:
    if( SetVertexCut(value)) {
      fCuts[kVertex] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kNCuts:
    AliError("Cut id out of range");
    return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;
}

///________________________________________________________________________
void AliConvEventCuts::PrintCuts() {
  // Print out current Cut Selection
  for(Int_t ic = 0; ic < kNCuts; ic++) {
    printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
  }
}

void AliConvEventCuts::PrintCutsWithValues() {
  // Print out current Cut Selection with value
  printf("\nEvent cutnumber \n");
  for(Int_t ic = 0; ic < kNCuts; ic++) {
    printf("%d",fCuts[ic]);
  }
  printf("\n\n");
  printf("EnergyVar-enum: '%i', PeriodVar-enum: '%i' \n", fEnergyEnum, fPeriodEnum );
  if (fIsHeavyIon == 0) {
    printf("Running in pp mode \n");
    if (fSpecialTrigger == 0){
      if(fSpecialTriggerName.Contains("INT7")){
        printf("\t only events triggered by V0AND will be analysed \n");
      }else if(fSpecialTriggerName.Contains("INT8")){
        printf("\t only events triggered by T0AND will be analysed \n");
      }else if(!fSpecialTriggerName.IsNull()){
        printf("\t only events triggered by %s will be analysed \n", fSpecialTriggerName.Data());
      }else{
        if (fSpecialSubTrigger == 0){
          printf("\t only events triggered by V0OR will be analysed \n");
        } else if (fSpecialSubTrigger == 1){
          printf("\t only events where SDD was present will be analysed \n");
        }
      }
    } else if (fSpecialTrigger == 1){
      if (fSpecialSubTrigger == 0){
        printf("\t only events triggered by V0AND will be analysed \n");
      } else if(fSpecialSubTrigger == 1){
        printf("\t only events where SDD was present will be analysed and triggered by VOAND\n");
      }
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps");
    } else if (fSpecialTrigger > 1){
      printf("\t only events triggered by %s %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data(), fSpecialSubTriggerNameAdditional.Data());
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps\n\n");
    }
    if ( !(fCentralityMin == 0 && fCentralityMax == 0) && !(fCentralityMax < fCentralityMin) ){
      printf("\t Multiplicity cut %d - %d \n", fCentralityMin, fCentralityMax);
    }
  } else if (fIsHeavyIon == 1){
    printf("Running in PbPb mode \n");
    if (fDetectorCentrality == 0){
      printf("\t centrality selection based on V0M \n");
    } else if (fDetectorCentrality == 1){
      printf("\t centrality selection based on Cl1 \n");
    } else if (fDetectorCentrality == 2){
      printf("\t centrality selection based on ZNA \n");
    }
    if (fModCentralityClass == 0){
      printf("\t %d - %d \n", fCentralityMin*10, fCentralityMax*10);
    } else if ( fModCentralityClass == 1){
      printf("\t %d - %d \n", fCentralityMin*5, fCentralityMax*5);
    } else if ( fModCentralityClass == 2){
      printf("\t %d - %d \n", fCentralityMin, fCentralityMax);
    } else if ( fModCentralityClass == 20){
      printf("\t %f - %f \n", fCentralityMin*0.1, fCentralityMax*0.1);
    } else if ( fModCentralityClass == 21){
      printf("\t %f - %f \n", fCentralityMin*0.01, fCentralityMax*0.01);
    } else if (fModCentralityClass == 3){
      printf("\t %d - %d, with Track mult in MC as data \n", fCentralityMin*10, fCentralityMax*10);
    } else if ( fModCentralityClass == 4){
      printf("\t %d - %d, with Track mult in MC as data \n", fCentralityMin*5, fCentralityMax*5);
    } else if (fModCentralityClass == 5){
      printf("\t %d - %d, with overlapping Track mult in MC as data \n", fCentralityMin*10, fCentralityMax*10);
    } else if ( fModCentralityClass == 6){
      printf("\t %d - %d, with overlapping Track mult in MC as data \n", fCentralityMin*5, fCentralityMax*5);
    }
    if (fSpecialTrigger == 0){
      printf("\t only events triggered by kMB, kCentral, kSemiCentral will be analysed \n");
    } else if (fSpecialTrigger > 1){
      printf("\t only events triggered by %s %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data(),fSpecialSubTriggerNameAdditional.Data());
      printf("\n\t        SpecialTrigger is:  %s\n", fSpecialTriggerName.Data());
      printf("\t        SpecialSubTrigger is: %s\n", fSpecialSubTriggerName.Data());
      if(fNSpecialSubTriggerOptions==2)
      printf("\t        SpecialSubTrigger2 is: %s\n", fSpecialSubTriggerNameAdditional.Data());
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps\n\n");
    }
  } else if (fIsHeavyIon == 2){
    printf("Running in pPb mode \n");
    if (fDetectorCentrality == 0){
      printf("\t centrality selection based on V0A \n");
    } else if (fDetectorCentrality == 1){
      printf("\t centrality selection based on Cl1 \n");
    } else if (fDetectorCentrality == 2){
      printf("\t centrality selection based on ZNA \n");
    }
    if (fModCentralityClass == 0){
      printf("\t %d - %d \n", fCentralityMin*10, fCentralityMax*10);
    }  else if ( fModCentralityClass == 1){
      printf("\t %d - %d \n", fCentralityMin*5, fCentralityMax*5);
    }  else if ( fModCentralityClass == 2){
      printf("\t %d - %d \n", fCentralityMin, fCentralityMax);
    }
    if (fSpecialTrigger == 0){
      printf("\t only events triggered by kINT7 will be analysed \n");
    } else if (fSpecialTrigger > 1){
      printf("\t only events triggered by %s %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data(),fSpecialSubTriggerNameAdditional.Data());
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps\n\n");
    }
  }
  if (fEnableVertexCut) printf("\t Vertex cut with |Z_{vtx}| <%2.2f \n",fMaxVertexZ);
    else printf("\t No vertex cut \n");

  if (fRemovePileUp ==1 ) {
     printf("\t Doing pile up removal  \n");
     if (fRemovePileUpSPD ==1 ){
       printf("\t Doing pile up removal using SPD \n");
     }
     if (fDoPileUpRejectV0MTPCout ==1 ){
       printf("\t Doing extra pile up removal V0M vs TPCout  \n");
     }
     if (fRemovePileUpSDDSSDTPC){
       printf("\t Doing extra pile up removal using SDD+SSD vs TPC clusters\n");
     }
     if (fPastFutureRejectionLow !=0 && fPastFutureRejectionHigh !=0 ){
       printf("\t Doing extra past-future pile up removal\n");
     }
  }

  printf("MC event cuts: \n");
  if (fRejectExtraSignals == 0) printf("\t no rejection was applied \n");
    else if (fRejectExtraSignals == 1) printf("\t only MB header will be inspected \n");
    else if (fRejectExtraSignals == 4) printf("\t special handling for Jets embedded in MB events \n");
    else if (fRejectExtraSignals > 1) printf("\t special header have been selected \n");
  printf("\t minimum factor between jet and pt hard = %2.2f \n", fMinFacPtHard);
  printf("\t maximum factor between jet and pt hard = %2.2f \n", fMaxFacPtHard);
  printf("\t maximum factor between pi0 or eta pt and pt hard = %2.2f \n", fMaxFacPtHardSingleParticle);
}

///________________________________________________________________________
Bool_t AliConvEventCuts::SetIsHeavyIon(Int_t isHeavyIon)
{   // Set Cut
  switch(isHeavyIon){
  case 0: // pp
    fIsHeavyIon=0;
    break;
  case 1: // V0M PbPb & XeXe
    // steps of 10%
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=0;
    break;
  case 2: // CL1 PbPb & XeXe
    // steps of 10%
    fIsHeavyIon=1;
    fDetectorCentrality=1;
    fModCentralityClass=0;
    break;
  case 3: // V0M PbPb & XeXe
    // steps of 5%
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=1;
    break;
  case 4: // V0M PbPb & XeXe & primary track mult for MC different track array
    // steps of 10%
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=5;
    break;
  case 5: // V0M PbPb & XeXe & primary track mult for MC
    // steps of 10%
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=3;
    break;
  case 6: // V0M PbPb & XeXe & primary track mult for MC
    // steps of 5%
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=4;
    break;
  case 7: // V0M PbPb & XeXe & primary track mult for MC different track array
    // steps of 5%
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=6;
    break;
  case 8: // pPb V0A
    // steps of 10%
    fIsHeavyIon=2;
    fDetectorCentrality=0;
    fModCentralityClass=0;
    break;
  case 9: // pPb CL1
    // steps of 10%
    fIsHeavyIon=2;
    fDetectorCentrality=1;
    fModCentralityClass=0;
    break;
  case 10: // a: pPb V0A
    // steps of 5%
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    fIsHeavyIon=2;
    fDetectorCentrality=0;
    fModCentralityClass=1;
    break;
  case 11: // b: pPb CL1
    // steps of 5%
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    fIsHeavyIon=2;
    fDetectorCentrality=1;
    fModCentralityClass=1;
    break;
  case 12: // c: pPb V0A
    // steps of 1%
    // 0 -0%, 1-1%, 2-2%, 3-3%, 4-4%, 5-5%, 6-6%, 7-7%, 8-8%, 9-9%, a-10%, b-11%, c-12%, d-13%, e-14%, f-15%, g-16%, h-17%, i-18%, j-19%, k-20%
    fIsHeavyIon=2;
    fDetectorCentrality=0;
    fModCentralityClass=2;
    break;
  case 13: // d: pPb CL1
    // steps of 1%
    // 0 -0%, 1-1%, 2-2%, 3-3%, 4-4%, 5-5%, 6-6%, 7-7%, 8-8%, 9-9%, a-10%, b-11%, c-12%, d-13%, e-14%, f-15%, g-16%, h-17%, i-18%, j-19%, k-20%
    fIsHeavyIon=2;
    fDetectorCentrality=1;
    fModCentralityClass=2;
    break;
  case 14: // e: pPb ZNA
    // steps of 10%
    fIsHeavyIon=2;
    fDetectorCentrality=2;
    fModCentralityClass=0;
    break;
  case 15: // f: pPb ZNA
    // steps of 5%
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    fIsHeavyIon=2;
    fDetectorCentrality=2;
    fModCentralityClass=1;
    break;
  case 16: // g: pPb ZNA
    // steps of 1%
    // 0 -0%, 1-1%, 2-2%, 3-3%, 4-4%, 5-5%, 6-6%, 7-7%, 8-8%, 9-9%, a-10%, b-11%, c-12%, d-13%, e-14%, f-15%, g-16%, h-17%, i-18%, j-19%, k-20%
    fIsHeavyIon=2;
    fDetectorCentrality=2;
    fModCentralityClass=2;
    break;
  case 17: // h: pp -> Sphericity cuts
    fIsHeavyIon=0;
    fUseSphericity=1;
    break;
  case 18: // i: pp -> Sphericity cuts & mult < 20
    fIsHeavyIon=0;
    fUseSphericity=2;
    break;
  case 19: // j: pp -> Sphericity cuts & mult > 20
    fIsHeavyIon=0;
    fUseSphericity=3;
    break;
  case 20: // k: pp -> Sphericity cuts & axis in EMCal
    fIsHeavyIon=0;
    fUseSphericity=4;
    break;
  case 21: // l: pp -> Sphericity cuts & axis not in EMCal
    fIsHeavyIon=0;
    fUseSphericity=5;
    break;
  case 22: // m: pp -> Multiplicity V0M in 1% bins
    fIsHeavyIon=0;
    fDetectorCentrality=0;
    fModCentralityClass=2;
    break;
  case 23: // n: pp -> Multiplicity V0M in 10% bins
    fIsHeavyIon=0;
    fDetectorCentrality=0;
    fModCentralityClass=0;
    break;
  case 24: // o: pp -> Multiplicity CL1 in 1% bins
    fIsHeavyIon=0;
    fDetectorCentrality=3;
    fModCentralityClass=2;
    break;
  case 25: // p: pp -> Multiplicity CL1 in 10% bins
    fIsHeavyIon=0;
    fDetectorCentrality=3;
    fModCentralityClass=0;
    break;
  case 26: // q: pp -> Multiplicity V0M in 0.1% bins
    fIsHeavyIon=0;
    fDetectorCentrality=0;
    fModCentralityClass=20;
    break;
  case 27: // r: pp -> Multiplicity V0M in 0.01% bins
    fIsHeavyIon=0;
    fDetectorCentrality=0;
    fModCentralityClass=21;
    break;
  case 28: // s: pp -> Multiplicity CL1 in 0.01% bins
    fIsHeavyIon=0;
    fDetectorCentrality=3;
    fModCentralityClass=21;
    break;
  case 29: // t: UPC
    fIsHeavyIon=0;
    break;
  default:
    AliError(Form("SetHeavyIon not defined %d",isHeavyIon));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliConvEventCuts::SetCentralityMin(Int_t minCentrality)
{
  // Set Cut
  if(minCentrality<0||minCentrality>20){
    AliError(Form("minCentrality not defined %d",minCentrality));
    return kFALSE;
  }

  fCentralityMin=minCentrality;
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliConvEventCuts::SetCentralityMax(Int_t maxCentrality)
{
  // Set Cut
  if(maxCentrality<0||maxCentrality>20){
    AliError(Form("maxCentrality not defined %d",maxCentrality));
    return kFALSE;
  }
  fCentralityMax=maxCentrality;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConvEventCuts::SetSelectSpecialTrigger(Int_t selectSpecialTrigger)
{
  // Set Cut
  switch(selectSpecialTrigger){
  case 0:
    fSpecialTrigger=0; // V0OR
    break;
  case 1:
    fSpecialTrigger=1; // V0AND
    break;
  // case 2:
  //   fSpecialTrigger=2; //
  //   break;
  case 3:
    fSpecialTrigger=3; //specific centrality trigger selection
    fOfflineTriggerMask=AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
    fTriggerSelectedManually = kTRUE;
    fSpecialTriggerName="AliVEvent::kCentral/kSemiCentral/kINT7";
    break;
  case 4:
    fSpecialTrigger=4; // trigger alias kTRD
    fOfflineTriggerMask=AliVEvent::kTRD;
    fTriggerSelectedManually = kTRUE;
    fSpecialTriggerName="AliVEvent::kTRD";
    break;
  case 5:
    fSpecialTrigger=5; // trigger alias kEMC
    fOfflineTriggerMask=AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 ;
    fTriggerSelectedManually = kTRUE;
    fTriggersEMCALSelected= 0;
    SETBIT(fTriggersEMCALSelected, kL0);
    fSpecialTriggerName="AliVEvent::kEMC7/kEMC8/kEMC1";
    break;
  case 6:
    fSpecialTrigger=6; // trigger alias kPHI
    fOfflineTriggerMask=AliVEvent::kPHI7 | AliVEvent::kPHI1 | AliVEvent::kPHI8 | AliVEvent::kPHOSPb;
    fTriggerSelectedManually = kTRUE;
    fSpecialTriggerName="AliVEvent::kPHI7/kPHI1/kPHI8/kPHOSPb";
    break;
  case 7:
    fSpecialTrigger=7; // trigger alias kHighMult
    fOfflineTriggerMask=AliVEvent::kHighMult;
    fTriggerSelectedManually = kTRUE;
    fSpecialTriggerName="AliVEvent::kHighMult";
    break;
  case 8:
    fSpecialTrigger=8; // trigger alias kEMCEGA
    fOfflineTriggerMask=AliVEvent::kEMCEGA;
    fTriggerSelectedManually = kTRUE;
    fTriggersEMCALSelected= 0;
    SETBIT(fTriggersEMCALSelected, kG2);
    fSpecialTriggerName="AliVEvent::kEMCEGA";
    break;
  case 9:
    fSpecialTrigger=9; // trigger alias kEMCEJE
    fOfflineTriggerMask=AliVEvent::kEMCEJE;
    fTriggerSelectedManually = kTRUE;
    fTriggersEMCALSelected= 0;
    SETBIT(fTriggersEMCALSelected, kJ2);
    fSpecialTriggerName="AliVEvent::kEMCEJE";
    break;
  case 10: //CALO and CALOFAST
    fSpecialTrigger=10; // trigger alias kEMC
    fOfflineTriggerMask=AliVEvent::kCaloOnly;
    fTriggerSelectedManually = kTRUE;
    fTriggersEMCALSelected= 0;
    fSpecialTriggerName="AliVEvent::kCaloOnly";
    break;
  case 11: // Double gap (DG) events
    fSpecialTrigger=11; // DG events
    fOfflineTriggerMask=0;  // kAny cannot be used for DG events
    fTriggerSelectedManually = kTRUE;
    fSpecialTriggerName="";
    break;
  case 12: // Ultra Peripheral Collision (UPC) events
    fSpecialTrigger=12; // UPC events
    fOfflineTriggerMask=0;  // kAny cannot be used for UPC events
    fTriggerSelectedManually = kTRUE;
    fSpecialTriggerName="";
    break;
  case 13: // d; software trigger applied on minimum bias data clusters
    fSpecialTrigger=13; // software trigger on data
    fSpecialTriggerName="";
    break;
  case 14: // e; software trigger applied on EMCal EG1 trigger data clusters
    fSpecialTrigger=14; // software trigger on data
    fOfflineTriggerMask=AliVEvent::kEMCEGA;
    fTriggerSelectedManually = kTRUE;
    fTriggersEMCALSelected= 0;
    SETBIT(fTriggersEMCALSelected, kG2);
    fSpecialTriggerName="AliVEvent::kEMCEGA";
    break;
  default:
    AliError(Form("Warning: Special Trigger %d Not known",selectSpecialTrigger));
    return 0;
  }
  return 1;
}

///________________________________________________________________________
Bool_t AliConvEventCuts::SetSelectSubTriggerClass(Int_t selectSpecialSubTriggerClass)
{
  // Set Cut
  if (fSpecialTrigger == 0){ //OR
    switch(selectSpecialSubTriggerClass){
    case 0://with VZERO
      fSpecialTrigger=0;
      fSpecialSubTrigger=0;
      // AliInfo("Info: Nothing to be done");
      break;
    case 3: //V0OR with SDD requested (will only work with LHC11a dataset)
      fSpecialSubTrigger=1;
      // cout << "V0OR with SDD requested" << endl;
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 1){ //AND with different detectors
    switch(selectSpecialSubTriggerClass){
    case 0: case 4:  //with VZERO general implementation of V0AND (periods LHC11c onwards)
      fSpecialTrigger=0;
      fSpecialSubTrigger=0;
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fSpecialTriggerName="AliVEvent::kINT7";
    break;
    case 1: //with TZERO
      fSpecialTrigger=0;
      fSpecialSubTrigger=0;
      fOfflineTriggerMask=AliVEvent::kINT8;
      fTriggerSelectedManually = kTRUE;
      fSpecialTriggerName="AliVEvent::kINT8";
      break;
    case 2: //with VZERO (will only work with LHC11a dataset)
      fSpecialTrigger=1;
      fSpecialSubTrigger=0;
      // AliInfo("Info: Nothing to be done");
      break;
    case 3: //V0AND with SDD requested (will only work with LHC11a dataset)
      fSpecialTrigger=1;
      fSpecialSubTrigger=1;
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 3){ // Selecting kCentral and kSemiCentral from trigger classes, not aliases
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // kCentral - no vertex restriction
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHN";
      // cout << "kCentralOpen" << endl;
      break;
    case 2: // kCentral - T00 +- 10 cm
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CCENT";
      // cout << "kCentralVertex" << endl;
      break;
    case 3: // kCentral - both
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHN|CCENT|CSEMI|CVLN";
      // cout << "kCentral both" << endl;
      break;
    case 4: // kSemiCentral - no vertex restriction
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVLN";
      // cout << "kSemiCentralOpen" << endl;
      break;
    case 5: // kSemiCentral - T00 +- 10 cm
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSEMI";
      // cout << "kSemiCentralVertex" << endl;
      break;
    case 6: // kSemiCentral - both
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSEMI%CVLN";
      // cout << "kSemiCentral both" << endl;
      break;
    case 7: // kMB
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPBI1_|CPBI1-";
      // cout << "kMB 1" << endl;
      break;
    case 8: // kMB
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPBI2_|CPBI2-";
      // cout << "kMB 2" << endl;
      break;
    case 9: // kMB
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPBI2_@CPBI2-@CPBI2_@CPBI2-";
      // cout << "kMB both" << endl;
      break;
    case 10: // 0V0M
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="C0V0M";
      break;
    case 11: // 0V0L
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="C0V0L";
      break;
    case 12: // 0VHM
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="C0VHM";
      break;
    case 13: // VOL7
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CV0L7";
      break;
    case 14: // 0STC
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="C0STC";
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 4){ // Subdivision of TRD trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // 7WUHSH - V0AND with single electron in TRD & EMCAL
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7WUHEE";
      break;
    case 2: // 8WUHSH - T0AND with single electron in TRD & EMCAL
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8WUHEE";
      break;
    case 3: // 7WUHSE - V0AND with single high pt electron in TRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7WUHSE";
      break;
    case 4: // 8WUHSE - T0AND with single high pt electron in TRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8WUHSE";
      break;
    case 5: // 7WUHJE - V0AND with jet in TRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7WUHJT";
      break;
    case 6: // 8WUHJE - T0AND with jet in TRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8WUHJT";
      break;
    case 7: // 7WUHQU - V0AND with dielectron pair in TRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7WUHQU";
      break;
    case 8: // 8WUHQU - T0AND with dielectron pair in TRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8WUHQU";
      break;
    case 9: // INT7HSE - V0AND with single high pt electron in TRD run 2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="INT7HSE";
      break;
    case 10: // INT7HQU - V0AND with dielectron  in TRD run 2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="INT7HQU";
      break;
    case 11: // INT7HJT - V0AND with jet in TRD run 2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="INT7HJT";
      break;
    case 12: // INT7HNU - V0AND with nuclei in TRD run 2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="INT7HNU";
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 5){ // Subdivision of kEMC trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // CEMC1 - V0OR and EMCAL fired
      fOfflineTriggerMask=AliVEvent::kEMC1;
      fSpecialTriggerName="AliVEvent::kEMC1";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CEMC1";
      break;
    case 2: // CEMC7 - V0AND and EMCAL fired
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kEMC7;
      fSpecialTriggerName="AliVEvent::kEMC7";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CEMC7";
      break;
    case 3: // CEMC8  - T0OR and EMCAL fired
      fOfflineTriggerMask=AliVEvent::kEMC8;
      fSpecialTriggerName="AliVEvent::kEMC8";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CEMC8";
      break;
    case 4: // CDMC1 - V0OR and DCAL fired
      fOfflineTriggerMask=AliVEvent::kEMC1;
      fSpecialTriggerName="AliVEvent::kEMC1";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CDMC1";
      break;
    case 5: // CDMC7 - V0AND and DCAL fired
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kEMC7;
      fSpecialTriggerName="AliVEvent::kEMC7";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CDMC7";
      break;
    case 6: // CDMC8  - T0OR and DCAL fired
      fOfflineTriggerMask=AliVEvent::kEMC8;
      fSpecialTriggerName="AliVEvent::kEMC8";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CDMC8";
      break;
    case 7: // MC7 - V0AND and EMCAL OR DCAL fired
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kEMC7;
      fSpecialTriggerName="AliVEvent::kEMC7";
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="CEMC7";
      fSpecialSubTriggerNameAdditional="CDMC7";
      break;
    case 8: // MC8  - T0OR and EMCAL OR DCAL fired
      fOfflineTriggerMask=AliVEvent::kEMC8;
      fSpecialTriggerName="AliVEvent::kEMC8";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="CEMC8";
      fSpecialSubTriggerNameAdditional="CDMC8";
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  }else if (fSpecialTrigger == 6){ // Subdivision of kPHI trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // CEMC1 - V0OR and PHOS fired
      fOfflineTriggerMask=AliVEvent::kPHI1;
      fSpecialTriggerName="AliVEvent::kPHI1";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPHI1";
      break;
    case 2: // CEMC7 - V0AND and PHOS fired
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kPHI7;
      fSpecialTriggerName="AliVEvent::kPHI7";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPHI7";
      break;
    case 3: // CEMC8  - T0OR and PHOS fired
      fOfflineTriggerMask=AliVEvent::kPHI8;
      fSpecialTriggerName="AliVEvent::kPHI8";
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPHI8";
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  }else if (fSpecialTrigger == 7){ // Subdivision of kHighMult trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // CSHM1 - V0OR and high mult fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSHM1";
      break;
    case 2: // CSHM7 - V0AND and high mult fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSHM7";
      break;
    case 3: // CSHM8  - T0OR and high mult fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSHM8";
      break;
    case 4: // V0 high mult trigger
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kAny;
      fSpecialTriggerName="V0Mult";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHMV0M-B-";
      break;
    case 5: // SPD high mult trigger
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kAny;
      fSpecialTriggerName="SPMult";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHMSH2-B-";
      break;
    case 6: // V0 high mult trigger with pileup condition on
      fSpecialSubTrigger=1;
      fOfflineTriggerMask=AliVEvent::kAny;
      fSpecialTriggerName="V0Mult";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHMV0M-B-SPD2";
      break;

    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  }else if (fSpecialTrigger == 8){ // Subdivision of kEMCEGA trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // 7EGA - CINT7 EGA
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EGA";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 2: // 8EGA - CINT8 EGA
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EGA";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 3: // 7EG1 - CINT7 EG1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 4: // 8EG1 - CINT8 EG1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 5: // 7EG2 - CINT7 EG2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 6: // 8EG2 - CINT8 EG2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 7: // 7DGA - CINT7 DGA
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DGA";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 8: // 8DGA - CINT8 DGA
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DGA";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 9: // 7DG1 - CINT7 DG1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 10: //a) 8DG1 - CINT8 DG1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 11: //b) 7DG2 - CINT7 DG2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 12: //c) 8DG2 - CINT8 DG2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 13: //d) Gamma Low EMC and DMC
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EG1";
      fSpecialSubTriggerNameAdditional="7DG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 14: //e) Gamma Low EMC and DMC
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EG2";
      fSpecialSubTriggerNameAdditional="7DG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 15: //f) Gamma Low EMC and DMC
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="8EG1";
      fSpecialSubTriggerNameAdditional="8DG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 16: //g) Gamma Low EMC and DMC
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="8EG2";
      fSpecialSubTriggerNameAdditional="8DG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;

    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 9){ // Subdivision of kEMCEJE trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0;
      fSpecialSubTriggerName="";
      // AliInfo("Info: Nothing to be done");
      break;
    case 1: // 7EJE - CINT7 EJE
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EJE";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
    case 2: // 8EJE - CINT8 EJE
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EJE";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
    case 3: // 7EJ1 - CINT7 EJ1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EJ1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ1);
      break;
    case 4: // 8EJ1 - CINT8 EJ1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EJ1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ1);
      break;
    case 5: // 7EJ2 - CINT7 EJ2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EJ2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
    case 6: // 8EJ2 - CINT8 EJ2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EJ2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
   case 7: // 7DJ1 - CINT7 DJ1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DJ1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ1);
      break;
    case 8: // 8DJ1 - CINT8 DJ1
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DJ1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ1);
      break;
    case 9: // 7DJ2 - CINT7 DJ2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DJ2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
    case 10: // 8DJ2 - CINT8 DJ2
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DJ2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
    case 11: // high Jet trigger EMC+DMC
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EJ1";
      fSpecialSubTriggerNameAdditional="7DJ1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ1);
      break;
    case 12: // low Jet trigger EMC+DMC
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EJ2";
      fSpecialSubTriggerNameAdditional="7DJ2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kJ2);
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 10){ // Subdivision of kEMC trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CINT7";
      fSpecialTriggerName="AliVEvent::kCaloOnly/INT7";
      break;
    case 1: // CEMC7 - V0AND and EMCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CEMC7-";
      fSpecialTriggerName="AliVEvent::kCaloOnly/EMC7";
      break;
    case 2: // CEMC7EG2 - V0AND and EMCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EG2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EG2";
      break;
    case 3: // CEMC7EG1  - V0AND and EMCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EG1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EG1";
      break;
    case 4: // CEMC7EJ2 - V0AND and EMCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EJ2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EJ2";
      break;
    case 5: // CEMC7EJ1 - V0AND and EMCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EJ1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EJ1";
      break;
    case 6: // CDMC7 - V0AND and DCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CDMC7-";
      fSpecialTriggerName="AliVEvent::kCaloOnly/DMC7";
      break;
    case 7: // CDMC7DG2 - V0AND and DCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DG2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7DG2";
      break;
    case 8: // CDMC7DG1  - V0AND and DCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DG1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7DG1";
      break;
    case 9: // CDMC7DJ2 - V0AND and DCAL fired
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DJ2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7DJ2";
      break;
    case 10: // DEMC7DJ1 - V0AND and DCAL fired - a
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DJ1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/7DJ1";
      break;
    case 11: // CEMC8 - V0AND and EMCAL fired - b
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CEMC8-";
      fSpecialTriggerName="AliVEvent::kCaloOnly/EMC8";
      break;
    case 12: // CEMC8EG2 - V0AND and EMCAL fired - c
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EG2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8EG2";
      break;
    case 13: // CEMC8EG1  - V0AND and EMCAL fired - d
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EG1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8EG1";
      break;
    case 14: // CEMC8EJ2 - V0AND and EMCAL fired - e
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EJ2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8EJ2";
      break;
    case 15: // CEMC8EJ1 - V0AND and EMCAL fired - f
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8EJ1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8EJ1";
      break;
    case 16: // CDMC8 - V0AND and DCAL fired - g
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CDMC8-";
      fSpecialTriggerName="AliVEvent::kCaloOnly/DMC8";
      break;
    case 17: // CDMC8DG2 - V0AND and DCAL fired - h
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DG2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8DG2";
      break;
    case 18: // CDMC8DG1  - V0AND and DCAL fired - i
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DG1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8DG1";
      break;
    case 19: // CDMC8DJ2 - V0AND and DCAL fired - j
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DJ2";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8DJ2";
      break;
    case 20: // DEMC8DJ1 - V0AND and DCAL fired - k
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DJ1";
      fSpecialTriggerName="AliVEvent::kCaloOnly/8DJ1";
      break;
    case 21: // Gamma Low EMC and DMC - l
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EG1";
      fSpecialSubTriggerName="7EG1";
      fSpecialSubTriggerNameAdditional="7DG1";
      break;
    case 22: // Gamma Low EMC and DMC - m
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EG2";
      fSpecialSubTriggerName="7EG2";
      fSpecialSubTriggerNameAdditional="7DG2";
      break;
    case 23: // high Jet trigger EMC+DMC - n
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EJ1";
      fSpecialSubTriggerName="7EJ1";
      fSpecialSubTriggerNameAdditional="7DJ1";
      break;
    case 24: // low Jet trigger EMC+DMC - o
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialTriggerName="AliVEvent::kCaloOnly/7EJ2";
      fSpecialSubTriggerName="7EJ2";
      fSpecialSubTriggerNameAdditional="7DJ2";
      break;
    case 25: // CPHI7 - V0AND and PHOS fired - p
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPHI7-";
      fSpecialTriggerName="AliVEvent::kCaloOnly/CPHI7";
      break;
    case 26: // V0AND and EMCAL OR DCAL fired - q
      fSpecialSubTrigger=1;
      fSpecialTriggerName="AliVEvent::kCaloOnly/AliVEvent::kEMC7";
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="CEMC7";
      fSpecialSubTriggerNameAdditional="CDMC7";
      break;
    case 27: // T0OR and EMCAL OR DCAL fired - r
      fSpecialSubTrigger=1;
      fSpecialTriggerName="AliVEvent::kCaloOnly/AliVEvent::kEMC8";
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="CEMC8";
      fSpecialSubTriggerNameAdditional="CDMC8";
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 11){ // selection of double gap events
    switch(selectSpecialSubTriggerClass){
    case 0: //CCUP25-B-SPD1-CENTNOTRD
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CCUP25-B-SPD1-CENTNOTRD";
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 12){ // selection of UPC events
    switch(selectSpecialSubTriggerClass){
    case 0: //CCUP8
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CCUP8";
      break;
    case 1: //CCUP9
      fSpecialSubTrigger=2;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CCUP9";
      break;
    default:
      AliError(Form("Warning: Special Subtrigger Class %d Not known",selectSpecialSubTriggerClass));
      return 0;
    }
  } else if (fSpecialTrigger == 13){ // software trigger on min bias data
    switch(selectSpecialSubTriggerClass){
    case 0: // mimick of MC7 - V0AND and EMCAL OR DCAL fired
      fSpecialSubTrigger=0;
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CEMC7_sw";
      fSpecialSubTriggerNameAdditional="CDMC7_sw";
      break;
    case 1: // mimick of L1 low - V0AND and EMCAL OR DCAL fired
      fSpecialSubTrigger=0;
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EG2_sw";
      fSpecialSubTriggerNameAdditional="7DG2_sw";
      // fTriggersEMCALSelected= 0;
      // fTriggersEMCALSelected= 0;
      // SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 2: // mimick of L1 high - V0AND and EMCAL OR DCAL fired
      fSpecialSubTrigger=0;
      fOfflineTriggerMask=AliVEvent::kINT7;
      fTriggerSelectedManually = kTRUE;
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7EG1_sw";
      fSpecialSubTriggerNameAdditional="7DG1_sw";
      // fTriggersEMCALSelected= 0;
      // SETBIT(fTriggersEMCALSelected, kG1);
      break;
    default:
      AliError(Form("Warning: Special Subtrigger Class %d Not known",selectSpecialSubTriggerClass));
      return 0;
    }
  } else if (fSpecialTrigger == 14){ // software trigger on EMCal triggered data
    switch(selectSpecialSubTriggerClass){
    case 0: //e0) Gamma High EMC and DMC software trigger on GammaHigh data
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EG1_EGA_sw";
      fSpecialSubTriggerNameAdditional="7DG1_EGA_sw";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 1: //e1) Gamma Low EMC and DMC  on Gamma Low data
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EG2_EGA_sw";
      fSpecialSubTriggerNameAdditional="7DG2_EGA_sw";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 2: //e2) Gamma High EMC and DMC software trigger on GammaLow data
      fSpecialSubTrigger=1;
      fNSpecialSubTriggerOptions=2;
      fSpecialSubTriggerName="7EG1_EGA_sw";
      fSpecialSubTriggerNameAdditional="7DG1_EGA_sw";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    default:
      AliError(Form("Warning: Special Subtrigger Class %d Not known",selectSpecialSubTriggerClass));
      return 0;
    }
  }
  return 1;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::SetMultiplicityMethod(Int_t multiplicityMethod)
{
  // Set Cut
  fMultiplicityMethod=multiplicityMethod;

  // 0 Photon Multiplicity
  // 1 TPC Track multiplicity
  // 2 V0 Mult
  // 3 SPD Mult

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::SetRemovePileUp(Int_t removePileUp)
{// Set Cut
  switch(removePileUp){
  case 0:
    fRemovePileUp           = kFALSE;
    break;
  case 1:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    break;
  case 2:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fPastFutureRejectionLow =-89;
    fPastFutureRejectionHigh= 89;
    break;
  case 3:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fPastFutureRejectionLow = -4;
    fPastFutureRejectionHigh=  7;
    break;
  case 4:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fPastFutureRejectionLow = -10;
    fPastFutureRejectionHigh=  13;
    break;
  case 5:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fPastFutureRejectionLow = -40;
    fPastFutureRejectionHigh=  43;
    break;
  case 6:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fDoPileUpRejectV0MTPCout = kTRUE;
    fFPileUpRejectV0MTPCout = new TF1("fFPileUpRejectV0MTPCout","[0] + [1]*x",0.,10000.);
    fFPileUpRejectV0MTPCout->SetParameter(0,0.);
    fFPileUpRejectV0MTPCout->SetParameter(1,0.);
    if (fIsHeavyIon==1){
      if(fPeriodEnum == kLHC15o){
         fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,5.0);
         break;
      }else{
         fFPileUpRejectV0MTPCout->SetParameter(0,-1500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,3.0);
         break;
      }
    } else  if(fIsHeavyIon == 2){
       fFPileUpRejectV0MTPCout->SetParameter(0,-200.);
       fFPileUpRejectV0MTPCout->SetParameter(1,2.0);
       break;
    }else{
       fFPileUpRejectV0MTPCout->SetParameter(0,-300.);
       fFPileUpRejectV0MTPCout->SetParameter(1,4.0);
       break;
    }
   break;
  case 7:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fDoPileUpRejectV0MTPCout = kTRUE;
    fFPileUpRejectV0MTPCout = new TF1("fFPileUpRejectV0MTPCout","[0] + [1]*x",0.,10000.);
    fFPileUpRejectV0MTPCout->SetParameter(0,0.);
    fFPileUpRejectV0MTPCout->SetParameter(1,0.);
    if (fIsHeavyIon==1){
      if(fPeriodEnum == kLHC15o){
         fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,5.0);
         break;
      }else{
         fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,3.0);
         break;
      }
    } else  if(fIsHeavyIon == 2){
       fFPileUpRejectV0MTPCout->SetParameter(0,-300.);
       fFPileUpRejectV0MTPCout->SetParameter(1,1.5);
       break;
    }else{
       fFPileUpRejectV0MTPCout->SetParameter(0,-300.);
       fFPileUpRejectV0MTPCout->SetParameter(1,3.0);
       break;
    }
    break;
  case 8:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fDoPileUpRejectV0MTPCout = kTRUE;
    fFPileUpRejectV0MTPCout = new TF1("fFPileUpRejectV0MTPCout","[0] + [1]*x",0.,10000.);
    fFPileUpRejectV0MTPCout->SetParameter(0,0.);
    fFPileUpRejectV0MTPCout->SetParameter(1,0.);
    if (fIsHeavyIon==1){
      if(fPeriodEnum == kLHC15o){
         fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,5.0);
         break;
      }else{
         fFPileUpRejectV0MTPCout->SetParameter(0,-1500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,3.0);
         break;
      }
    } else  if(fIsHeavyIon == 2){
       fFPileUpRejectV0MTPCout->SetParameter(0,-200.);
       fFPileUpRejectV0MTPCout->SetParameter(1,1.5);
       break;
    }else{
       fFPileUpRejectV0MTPCout->SetParameter(0,-300.);
       fFPileUpRejectV0MTPCout->SetParameter(1,4.0);
       break;
    }
   break;
  case 9:
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fPastFutureRejectionLow =-89;
    fPastFutureRejectionHigh= 89;
    fDoPileUpRejectV0MTPCout = kTRUE;
    fFPileUpRejectV0MTPCout = new TF1("fFPileUpRejectV0MTPCout","[0] + [1]*x",0.,10000.);
    fFPileUpRejectV0MTPCout->SetParameter(0,0.);
    fFPileUpRejectV0MTPCout->SetParameter(1,0.);
    if (fIsHeavyIon==1){
      if(fPeriodEnum == kLHC15o){
         fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,5.0);
         break;
      }else{
         fFPileUpRejectV0MTPCout->SetParameter(0,-1500.);
         fFPileUpRejectV0MTPCout->SetParameter(1,3.0);
         break;
      }
    } else  if(fIsHeavyIon == 2){
       fFPileUpRejectV0MTPCout->SetParameter(0,-200.);
       fFPileUpRejectV0MTPCout->SetParameter(1,1.5);
       break;
    }else{
       fFPileUpRejectV0MTPCout->SetParameter(0,-300.);
       fFPileUpRejectV0MTPCout->SetParameter(1,4.0);
       break;
    }
   break;
 case 10: // a           for Pb-Pb
    fRemovePileUp     = kTRUE;
    fRemovePileUpSPD  = kTRUE;
    if(fPeriodEnum == kLHC18qr){
      fUtils->SetASPDCvsTCut(750.);
      fUtils->SetBSPDCvsTCut(4.);
    } else {
      fUtils->SetASPDCvsTCut(200.);
      fUtils->SetBSPDCvsTCut(7.);
    }
    fDoPileUpRejectV0MTPCout = kTRUE;
    fFPileUpRejectV0MTPCout = new TF1("fFPileUpRejectV0MTPCout","[0] + [1]*x",0.,10000.);
    if(fPeriodEnum == kLHC18qr){
      fFPileUpRejectV0MTPCout->SetParameter(0,-2000.);
      fFPileUpRejectV0MTPCout->SetParameter(1,6.0);
    } else {
      fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
      fFPileUpRejectV0MTPCout->SetParameter(1,5.0);
    }
    break;
 case 11: // b            for Pb-Pb
    fRemovePileUp     = kTRUE;
    fRemovePileUpSPD  = kFALSE;
    fDoPileUpRejectV0MTPCout = kTRUE;
    fFPileUpRejectV0MTPCout = new TF1("fFPileUpRejectV0MTPCout","[0] + [1]*x",0.,10000.);
    fFPileUpRejectV0MTPCout->SetParameter(0,-2500.);
    fFPileUpRejectV0MTPCout->SetParameter(1,5.0);
    break;
 case 12: // c
    fRemovePileUp           = kTRUE;
    fRemovePileUpSPD        = kTRUE;
    fUtils->SetASPDCvsTCut(200.);
    fUtils->SetBSPDCvsTCut(7.);
    break;
  case 13: // d         for Pb-Pb LHC18qr
    fRemovePileUp = kTRUE;
    fRemovePileUpSDDSSDTPC = kTRUE;

    fFPileUpRejectSDDSSDTPC = new TF1("fFPileUpRejectSDDSSDTPC", "[0]+[1]*x+[2]*x*x", 0., 1.e+7);
    fFPileUpRejectSDDSSDTPC->SetParameters(-3000., 0.0099,9.426e-10);

    break;
  default:
    AliError("RemovePileUpCut not defined");
    return kFALSE;
  }
  if (fDoPileUpRejectV0MTPCout) fLabelNamePileupCutTPC = "Pileup V0M-TPCout Tracks";
  if (fFPileUpRejectSDDSSDTPC)  fLabelNamePileupCutTPC = "Pileup SDD+SSD-TPC clusters";

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::SetRejectExtraSignalsCut(Int_t extraSignal) {

  switch(extraSignal){
  case 0:
    fRejectExtraSignals = 0;
    break; // No Rejection
  case 1:
    fRejectExtraSignals = 1;
    break; // MinBias Header
  case 2:
    fRejectExtraSignals = 2;
    break; // User String Array
  case 3:
    fRejectExtraSignals = 3;
    break; // Rejection for Gamma Correction only
  case 4:
    fRejectExtraSignals = 4;
    break; // Special handling of Jet weights for Jets embedded in MB events
  default:
    AliError(Form("Extra Signal Rejection not defined %d",extraSignal));
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::SetVertexCut(Int_t vertexCut) {

  switch(vertexCut){
  case 0: // no Vertex required // NOT fully working yet
    fEnableVertexCut   = kFALSE;
    fMaxVertexZ     = 1000;
    break;
  case 1: // vertex within +-15 cm
    fEnableVertexCut   = kTRUE;
    fMaxVertexZ     = 15;
    break;
  case 2: // vertex within +-12.5 cm
    fEnableVertexCut   = kTRUE;
    fMaxVertexZ     = 12.5;
    break;
  case 3: // vertex within +-10 cm
    fEnableVertexCut   = kTRUE;
    fMaxVertexZ     = 10.0;
    break;
  case 4: // vertex within +-7.5 cm
    fEnableVertexCut   = kTRUE;
    fMaxVertexZ     = 7.5;
    break;
  case 5: // vertex within +-5 cm
    fEnableVertexCut   = kTRUE;
    fMaxVertexZ     = 5.;
    break;
  case 6: // vertex within +-2.5 cm
    fEnableVertexCut   = kTRUE;
    fMaxVertexZ     = 2.5;
    break;
  default:
    AliError(Form("Vertex Cut not defined %d",vertexCut));
    return kFALSE;
  }
  return kTRUE;
}

//-------------------------------------------------------------
Bool_t AliConvEventCuts::GetUseNewMultiplicityFramework(){
  switch (fPeriodEnum){
    // pp 5TeV
    case kLHC15n :
    case kLHC17pq :
    // pp 5TeV MC
    case kLHC15l1a2 :
    case kLHC15l1b2 :
    case kLHC16h8a :
    case kLHC16h8b :
    case kLHC16k3a2 :
    case kLHC16k3a :
    case kLHC16k5a :
    case kLHC16k5b :
    case kLHC17e2 :
    case kLHC18j3 :
    case kLHC16h3 :
    case kLHC17l3b :
    case kLHC18j2 :
    case kLHC17l4b :
    case kLHC18b8 :
    case kLHC18b10 :
    case kLHC18l2 :
    case kLHC17P1PHO :
    // pp 13 TeV
    case kLHC15fm :
    case kLHC16NomB :
    case kLHC16LowB :
    case kLHC17NomB :
    case kLHC17LowB :
    case kLHC18NomB :
    case kLHC18LowB :
    // pp 13 TeV MC
    case kLHC15g3a3 :
    case kLHC15g3a :
    case kLHC15g3c2 :
    case kLHC15g3c3 :
    case kLHC15g3 :
    case kLHC16a2a :
    case kLHC16a2b :
    case kLHC16a2c :
    case kLHC15P2EPos :
    case kLHC15P2Pyt8 :
    case kLHC15k5a :
    case kLHC15k5b :
    case kLHC15k5c :
    case kLHC16P1Pyt8 :
    case kLHC16P1Pyt8LowB:
    case kLHC16P1EPOS :
    case kLHC16P1PHO :
    case kLHC16P1JJ :
    case kLHC16P1JJLowB:
    case kLHC17h8a :
    case kLHC17h8b :
    case kLHC17h8c :
    case kLHC17c3b1 :
    case kLHC17c3a1 :
    case kLHC17c3b2 :
    case kLHC17c3a2 :
    case kLHC17i3a1 :
    case kLHC17i3b1 :
    case kLHC17i3b2 :
    case kLHC17i3c1 :
    case kLHC17i3c2 :
    case kLHC20b1b1 :
    case kLHC20b1b2 :
    case kLHC20b1c1 :
    case kLHC20b1c2 :
    case kLHC17P1Pyt8NomB :
    case kLHC17P1Pyt6NomB :
    case kLHC17P1PHONomB13TeV :
    case kLHC17P1Pyt8LowB :
    case kLHC17j5a :
    case kLHC17j5b :
    case kLHC17j5c :
    case kLHC17P1JJ :
    case kLHC17P1JJLowB :
    case kLHC18l6b1 :
    case kLHC18l6b2 :
    case kLHC18l6c1 :
    case kLHC18l6c2 :
    case kLHC18P1JJ :
    case kLHC18P1Pyt8NomB :
    case kLHC18P1Pyt8LowB :
    case kLHC19i3b1 :
    case kLHC19i3b2 :
    case kLHC19i3c1 :
    case kLHC19i3c2 :
    // pPb 5 TeV
    case kLHC13bc :
    case kLHC13de :
    case kLHC13f :
    case kLHC16qt :
    // pPb 5 TeV MC
    case kLHC13b2_efix :
    case kLHC13e7 :
    case kLHC14b2 :
    case kLHC18j5 :
    case kLHC13b4_fix :
    case kLHC13b4_plus :
    case kLHC19a4 :
    case kLHC16c3a :
    case kLHC16c3b :
    case kLHC16c3c :
    case kLHC17g6a1 :
    case kLHC17g6a2 :
    case kLHC17g6a3 :
    case kLHC17f2a :
    case kLHC17f2b :
    case kLHC18f3 :
    case kLHC17g8a :
    case kLHC17d2a :
    case kLHC17d2b :
    // pPb 8 TeV
    case kLHC16r :
    case kLHC16s :
    // pPb 8 TeV MC
    case kLHC17a3a :
    case kLHC17a3b :
    case kLHC17a4a :
    case kLHC17a4b :
    case kLHC18f3bc :
    case kLHC17f3 :
    case kLHC17f3a :
    case kLHC17f3b :
    case kLHC17f4 :
    case kLHC17f4a :
    case kLHC17f4b :
    case kLHC17g6b2a :
    case kLHC17g6b2b :
    case kLHC17g6b3a :
    case kLHC17g6b3b :
    case kLHC16rP1JJ :
    case kLHC16sP1JJ :
    case kLHC16rsGJ :
    // Xe-Xe 5.44 TeV
    case kLHC17n :
    // Xe-Xe 5.44 TeV MC
    case kLHC17XeXeHi :
    // PbPb 5 TeV
    case kLHC15o :
    case kLHC18qr :
    // PbPb 5 TeV MC
    case kLHC15k1a1 :
    case kLHC15k1a2 :
    case kLHC15k1a3 :
    case kLHC16j7 :
    case kLHC16g2 :
    case kLHC16g3 :
    case kLHC16h4 :
    case kLHC16i1a :
    case kLHC16i1b :
    case kLHC16i1c :
    case kLHC16i2a :
    case kLHC16i2b :
    case kLHC16i2c :
    case kLHC16i3a :
    case kLHC16i3b :
    case kLHC16i3c :
    case kLHC16h2a :
    case kLHC16h2b :
    case kLHC16h2c :
    case kLHC16k3b :
    case kLHC16k3b2 :
    case kLHC18b11a :
    case kLHC18b11b :
    case kLHC18b11c :
    case kLHC18e1 :
    case kLHC18e1a :
    case kLHC18e1b :
    case kLHC18e1c :
    case kLHC18l8a :
    case kLHC18l8b :
    case kLHC18l8c :
    case kLHC19h2a :
    case kLHC19h2b :
    case kLHC19h2c :
    case kLHC19h3 :
    case kLHC20e3a :
    case kLHC20e3b :
    case kLHC20e3c :
    case kLHC20g10 :
      return kTRUE;
      break;
    default :
      return kFALSE;
      break;
  }
  return kFALSE;
}

//-------------------------------------------------------------
Float_t AliConvEventCuts::GetCentrality(AliVEvent *event)
{   // Get Event Centrality

  AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
  Int_t runnumber = event->GetRunNumber();
  if(esdEvent){
    if(GetUseNewMultiplicityFramework()){
      AliMultSelection *MultSelection = (AliMultSelection*)event->FindListObject("MultSelection");
      if(!MultSelection){
        AliWarning ("AliMultSelection object not found !");
        return -1;
      }else{
        if(fDetectorCentrality==0){
          if(fIsHeavyIon==2){
            if (runnumber > 266329 && runnumber < 267139)
              return MultSelection->GetMultiplicityPercentile("V0C");// default for Pbp
            else
              return MultSelection->GetMultiplicityPercentile("V0A");// default for pPb
          } else {
            return MultSelection->GetMultiplicityPercentile("V0M");// default
          }
        } else if(fDetectorCentrality==1){
          return MultSelection->GetMultiplicityPercentile("CL1",kTRUE);
        } else if(fDetectorCentrality==2){
          if (runnumber > 266329 && runnumber < 267139)
            return MultSelection->GetMultiplicityPercentile("ZNC",kTRUE);
          else
            return MultSelection->GetMultiplicityPercentile("ZNA",kTRUE);
        } else if(fDetectorCentrality==3){
          return MultSelection->GetMultiplicityPercentile("SPDTracklets",kTRUE);
        }
      }
    }else{
      AliCentrality *fESDCentrality = (AliCentrality*)esdEvent->GetCentrality();
      if(fDetectorCentrality==0){
        if(fIsHeavyIon==2){
          if (runnumber > 196432 && runnumber < 197389)
            return fESDCentrality->GetCentralityPercentile("V0C"); // default for Pbp
          else
            return fESDCentrality->GetCentralityPercentile("V0A"); // default for pPb
        } else {
          return fESDCentrality->GetCentralityPercentile("V0M"); // default
        }
      } else if(fDetectorCentrality==1){
        return fESDCentrality->GetCentralityPercentile("CL1");
      } else if(fDetectorCentrality==2){
        if (runnumber > 196432 && runnumber < 197389)
          return fESDCentrality->GetCentralityPercentile("ZNC");
        else
          return fESDCentrality->GetCentralityPercentile("ZNA");
      }
    }
  }
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);
  if(aodEvent){
    if(GetUseNewMultiplicityFramework()){
      AliMultSelection *MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
      if(!MultSelection){
        AliWarning ("AliMultSelection object not found !");
        return -1;
            } else{
        if(fDetectorCentrality==0){
          if(fIsHeavyIon==2){
            if (runnumber > 266329 && runnumber < 267139)
              return MultSelection->GetMultiplicityPercentile("V0C");// default for Pbp
            else
              return MultSelection->GetMultiplicityPercentile("V0A");// default for pPb
          } else if ( fPeriodEnum == kLHC18qr ) {
            return MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
          } else {
            return MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
          }
        }else if(fDetectorCentrality==1){
          return MultSelection->GetMultiplicityPercentile("CL1",kTRUE);
        }else if(fDetectorCentrality==2) {
          if (runnumber > 266329 && runnumber < 267139)
            return MultSelection->GetMultiplicityPercentile("ZNC",kTRUE);
          else
            return MultSelection->GetMultiplicityPercentile("ZNA",kTRUE);
        } else if(fDetectorCentrality==3){
          return MultSelection->GetMultiplicityPercentile("SPDTracklets",kFALSE);
        }
    }
    }else{
      if(aodEvent->GetHeader()){return ((AliVAODHeader*)aodEvent->GetHeader())->GetCentrality();}
    }
  }
  return -1;
}

//_____________________________________________________________________________________
Bool_t AliConvEventCuts::IsCentralitySelected(AliVEvent *event, AliMCEvent *mcEvent)
{
  // Centrality Selection
  if(!fIsHeavyIon){
    if ((fCentralityMin == 0 && fCentralityMax == 0) || (fCentralityMin > fCentralityMax) || (fUseSphericity > 0) ){
      return kTRUE;
    } else if (fDetectorCentrality == -1){
      Int_t primaryTracksPP[9] = { 0,   2,   5,    10,   15,
                                  30,  50,  100,  1000
                                  };
      Int_t nprimaryTracks = GetV0Reader()->GetNumberOfPrimaryTracks();
      if ( nprimaryTracks >= primaryTracksPP[fCentralityMin] && nprimaryTracks < primaryTracksPP[fCentralityMax]){
        return kTRUE;
      } else {
        return kFALSE;
      }
      return kFALSE;
    }
  }
  if(fCentralityMin == fCentralityMax ) return kTRUE;//0-100%
  else if ( fCentralityMax==0) fCentralityMax=10; //CentralityRange = fCentralityMin-10*multfactor
  Double_t centrality=GetCentrality(event);
  if(centrality<0 && !mcEvent)return kFALSE;
  Double_t addMarginZNA = (fDetectorCentrality==2 && fCentralityMax==10) ? 2.0 : 0.0; // For ZNA multiplicity goes up to 101%

  Int_t centralityC=0;
  if (fModCentralityClass == 0){
    centralityC= Int_t(centrality);
    if(centralityC >= (fCentralityMin*10) && centralityC < (fCentralityMax*10 + addMarginZNA))
      return kTRUE;
    else return kFALSE;
  }
  else if (fModCentralityClass == 1){
    centralityC= Int_t(centrality);
    if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
      return kTRUE;
    } else return kFALSE;
  }
  else if (fModCentralityClass == 2){
    centralityC= Int_t(centrality);
    if(centralityC >= fCentralityMin && centralityC < fCentralityMax){
      return kTRUE;
    } else return kFALSE;
  }
  else if (fModCentralityClass == 20){  // pp 13 TeV 0.1% mult classes
    centralityC= Int_t(centrality*10);
    if(centralityC >= fCentralityMin && centralityC < fCentralityMax){
      return kTRUE;
    } else return kFALSE;
  }
  else if (fModCentralityClass == 21){  // pp 13 TeV 0.01% mult classes
    centralityC= Int_t(centrality*100);
    if(centralityC >= fCentralityMin && centralityC < fCentralityMax){
      return kTRUE;
    } else return kFALSE;
  }
  Int_t nprimaryTracks = GetV0Reader()->GetNumberOfPrimaryTracks();
  Int_t PrimaryTracks10[11][2] =
    {
      {9999,9999}, //  0 //1550 changed to 9999 on 9 Dec
      {1210, 928}, // 10
      { 817, 658}, // 20
      { 536, 435}, // 30
      { 337, 276}, // 40
      { 197, 162}, // 50
      { 106, 100}, // 60
      {  51,  44}, // 70
      {  21,  18}, // 80
      {   0,   0},  // 90
      {   0,   0}// 100 // only max accessible
    };
  Int_t PrimaryTracksLHC11h10[11][2] =
    {
      {9999,9999}, //  0 //1550 changed to 9999 on 9 Dec
      { 985, 928}, // 10
      { 661, 658}, // 20
      { 434, 435}, // 30
      { 275, 276}, // 40
      { 173, 162}, // 50
      { 100, 100}, // 60
      {  42,  44}, // 70
      {  19,  18}, // 80
      {   0,   0},  // 90
      {   0,   0}// 100 // only max accessible
    };
  Int_t PrimaryTracksLHC15o10[11][2] =
    {
      {2500, 2700}, // 0-10% cent class max # of tracks: max value of the data distribution
      {1498, 1498}, // 0-10% cent class min # of tracks
      {1012, 1012}, // 10-20
      { 669,  669}, // 20-30
      { 423,  423}, // 30-40
      { 251,  251}, // 40-50
      { 136,  136}, // 50-60
      {  67,   67}, // 60-70
      {  28,   28}, // 70-80
      {   0,    0}, // 80-90% cent class min # of tracks
      {   0,    0}  // not used
    };
  Int_t PrimaryTracksLHC10h5[21][2] =
    {
      {9999,9999},  // 0 ///1550 changed to 9999 on 9 Dec
      {1485,1168},  // 5
      {1210, 928},  // 10
      { 995, 795},  // 15
      { 817, 658},  // 20
      { 666, 538},  // 25
      { 536, 435},  // 30
      { 428, 350},  // 35
      { 337, 276},  // 40
      { 260, 214},  // 45
      { 197, 162},  // 50
      { 147, 125},  // 55
      { 106, 100},  // 60
      {  75,  63},  // 65
      {  51,  44},  // 70
      {  34,  29},  // 75
      {  21,  18},  // 80
      {  13,  11},  // 85
      {   6,   6},  // 90
      {   3,   3},  // 95
      {   0,   0}   // 100 only max accessible
    };
    Int_t PrimaryTracksLHC11h5[21][2] =
    {
      {9999,9999},  // 0 ///1550 changed to 9999 on 9 Dec
      {1166,1168},  // 5
      { 953, 928},  // 10
      { 805, 795},  // 15
      { 655, 658},  // 20
      { 535, 538},  // 25
      { 435, 435},  // 30
      { 349, 350},  // 35
      { 275, 276},  // 40
      { 214, 214},  // 45
      { 165, 162},  // 50
      { 127, 125},  // 55
      {  93, 100},  // 60
      {  64,  63},  // 65
      {  44,  44},  // 70
      {  30,  29},  // 75
      {  18,  18},  // 80
      {  11,  11},  // 85
      {   6,   6},  // 90
      {   3,   3},  // 95
      {   0,   0}   // 100 only max accessible
    };
    Int_t PrimaryTracksLHC15o5[21][2] =
    {
      { 2500, 2700},  // 0-5% cent class max # of tracks: max value of the data distribution
      { 1827, 1827},  // 0-5% cent class min # of tracks
      { 1498, 1498},  // 5-10
      { 1234, 1234},  // 10-15
      { 1012, 1012},  // 15-20
      {  827,  827},  // 20-25
      {  669,  669},  // 25-30
      {  536,  536},  // 30-35
      {  423,  423},  // 35-40
      {  329,  329},  // 40-45
      {  251,  251},  // 45-50
      {  188,  188},  // 50-55
      {  136,  136},  // 55-60
      {   97,   97},  // 60-65
      {   67,   67},  // 65-70
      {   44,   44},  // 70-75
      {   28,   28},  // 75-80
      {   17,   17},  // 80-85
      {   10,   10},  // 85-90
      {    5,    5},  // 90-95 cent class minimum # of tracks
      {    0,    0}   // 95-100
    };
    Int_t PrimaryTracksLHC17n10[11][2] =
    {
        {9999,9999}, //  0 // 1500 max in hist but set to real max
        { 800, 800}, // 10 // guess
        { 628, 628}, // 20
        { 350, 350}, // 30 // guess
        { 268, 268}, // 40
        { 200, 200}, // 50 // guess
        { 100, 100}, // 60 // guess
        {  51,  44}, // 70 // guess
        {  21,  18}, // 80 // guess
        {   0,   0}, // 90 // guess
        {   0,   0}  // 100 // only max accessible
    };

  Int_t column = 0;
  if(event->IsA()==AliESDEvent::Class()) column = 0;
  if(event->IsA()==AliAODEvent::Class()) column = 1;

  if (fModCentralityClass == 3){
    if(mcEvent){
      // setting specific arry for LHC11h for MC track mult
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h10[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h10[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      // setting specific arry for LHC17n for MC track mult
      } else if(fPeriodEnum == kLHC17XeXeHi ){
        if(nprimaryTracks > PrimaryTracksLHC17n10[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC17n10[fCentralityMin][column])
            return kTRUE;
        else return kFALSE;
      // settings for LHC15o MCs
      } else if( fPeriodEnum == kLHC18e1 || fPeriodEnum == kLHC18e1a || fPeriodEnum == kLHC18e1b || fPeriodEnum == kLHC18e1c || fPeriodEnum == kLHC16h4 || fPeriodEnum == kLHC16i1a ||
                 fPeriodEnum == kLHC16i1b || fPeriodEnum == kLHC16i1c || fPeriodEnum == kLHC16i2a || fPeriodEnum == kLHC16i2b || fPeriodEnum == kLHC16i2c || fPeriodEnum == kLHC16i3a ||
                 fPeriodEnum == kLHC16i3b || fPeriodEnum == kLHC16i3c){
        centralityC= Int_t(centrality/10);
        if(centralityC >= fCentralityMin && centralityC < fCentralityMax){
            if(fCentralityMin==0 && nprimaryTracks >= PrimaryTracksLHC15o10[0][column]) return kFALSE;
            else return kTRUE;
        } else return kFALSE;
      // setting specific arry for LHC10h for MC track mult
      } else {
        if(nprimaryTracks > PrimaryTracks10[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks10[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      }
    } else {
      centralityC= Int_t(centrality/10);
      if(centralityC >= fCentralityMin && centralityC < fCentralityMax)
        return kTRUE;
      else return kFALSE;
    }
  } else if (fModCentralityClass ==4){
    if(mcEvent){
      // setting specific arry for LHC11h for MC track mult
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h5[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h5[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      // settings for LHC15o MCs
      } else if(fPeriodEnum == kLHC18e1 || fPeriodEnum == kLHC18e1a || fPeriodEnum == kLHC18e1b || fPeriodEnum == kLHC18e1c || fPeriodEnum == kLHC16h4 ||
                fPeriodEnum == kLHC16i1a || fPeriodEnum == kLHC16i1b || fPeriodEnum == kLHC16i1c || fPeriodEnum == kLHC16i2a || fPeriodEnum == kLHC16i2b ||
                fPeriodEnum == kLHC16i2c || fPeriodEnum == kLHC16i3a || fPeriodEnum == kLHC16i3b || fPeriodEnum == kLHC16i3c){
        centralityC = Int_t(centrality);
        if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
            if(fCentralityMin==0 && nprimaryTracks >= PrimaryTracksLHC15o5[0][column]) return kFALSE;
            else return kTRUE;
        } else return kFALSE;
        // setting specific arry for LHC10h for MC track mult
      } else {
        if(nprimaryTracks > PrimaryTracksLHC10h5[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC10h5[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      }
    }
    else{
      centralityC= Int_t(centrality);
      if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
        return kTRUE;
      } else return kFALSE;
    }
  }

  Int_t PrimaryTracksLHC11h10AltMin[11][2] =
  {
    {1550,1550}, //  0  - 0
    { 800, 800}, // 10  - 1
    { 600, 600}, // 20  - 2
    { 400, 400}, // 30  - 3
    { 240, 240}, // 40  - 4
    { 130, 130}, // 50  - 5
    {  90,  90}, // 60  - 6
    {  35,  35}, // 70  - 7
    {  15,  15}, // 80  - 8
    {   5,   5}, // 90  - 9
    {   0,   0}  // 100 // only max accessible
  };
  Int_t PrimaryTracksLHC11h10AltMax[11][2] =
  {
    {1550,1550}, //  0 //1550 changed to 9999 on 9 Dec
    {1000,1000}, // 10
    { 700, 700}, // 20
    { 480, 480}, // 30
    { 300, 300}, // 40
    { 200, 200}, // 50
    { 120, 120}, // 60
    {  50,  50}, // 70
    {  22,  22}, // 80
    {  10,  10}, // 90
    {   0,   0}  // 100 // only max accessible
  };
  Int_t PrimaryTracksLHC11h5AltMin[21][2] =
  {
    {1550,1550},  // 0
    {1000,1000},  // 5
    { 800, 800},  // 10
    { 700, 700},  // 15
    { 600, 600},  // 20
    { 500, 500},  // 25
    { 400, 400},  // 30
    { 300, 300},  // 35
    { 240, 240},  // 40
    { 180, 180},  // 45
    { 130, 130},  // 50
    { 127, 125},  // 55
    {  90,  90},  // 60
    {  55,  55},  // 65
    {  35,  35},  // 70
    {  25,  25},  // 75
    {  15,  15},  // 80
    {  11,  11},  // 85
    {   5,   5},  // 90
    {   0,   0},  // 95
    {   0,   0}   // 100 only max accessible
  };
  Int_t PrimaryTracksLHC11h5AltMax[21][2] =
  {
    {1550,1550},  // 0
    {1250,1250},  // 5
    {1000,1000},  // 10
    { 805, 795},  // 15
    { 700, 700},  // 20
    { 585, 585},  // 25
    { 480, 480},  // 30
    { 380, 380},  // 35
    { 300, 300},  // 40
    { 235, 235},  // 45
    { 200, 200},  // 50
    { 140, 140},  // 55
    { 120, 120},  // 60
    {  70,  70},  // 65
    {  50,  50},  // 70
    {  35,  25},  // 75
    {  22,  22},  // 80
    {  15,  15},  // 85
    {  10,  10},  // 90
    {   5,   5},  // 95
    {   0,   0}   // 100 only max accessible
  };

  if (fModCentralityClass == 5){
    if(mcEvent){
      // setting specific arry for LHC11h for MC track mult
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h10AltMin[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h10AltMax[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      // default return
      } else {
        return kFALSE;
      }
    } else {
      centralityC= Int_t(centrality/10);
      if(centralityC >= fCentralityMin && centralityC < fCentralityMax)
        return kTRUE;
      else return kFALSE;
    }
  }
  else if (fModCentralityClass ==6){
    if(mcEvent){
      // setting specific arry for LHC11h for MC track mult
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h5AltMin[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h5AltMax[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      } else {
        return kFALSE;
      }
    } else{
      centralityC= Int_t(centrality);
      if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
        return kTRUE;
      } else return kFALSE;
    }
  }

  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::VertexZCut(AliVEvent *event){
  // Cut on z position of primary vertex
  Double_t fVertexZ=event->GetPrimaryVertex()->GetZ();
  Double_t fVertexZSPD = 0;
  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(event);
  if(fESDEvent){
    fVertexZSPD = fESDEvent->GetPrimaryVertexSPD()->GetZ();
  }
  AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(event);
  if(fAODEvent){
    fVertexZSPD = fAODEvent->GetPrimaryVertexSPD()->GetZ();
  }

  if(TMath::Abs(fVertexZ)>fMaxVertexZ)return kFALSE;


  if (fPeriodEnum == kLHC11h){
    if (TMath::Abs(fVertexZ-fVertexZSPD) > 0.1) return kFALSE;
  }
  if (fIsHeavyIon == 2){
    if(!fUtils->IsVertexSelected2013pA(event)) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::IsOutOfBunchPileupPastFuture(AliVEvent *event)
{
  if(fPastFutureRejectionLow==0 && fPastFutureRejectionHigh==0)
    return kFALSE;
  TBits fIR1 =  event->GetHeader()->GetIRInt1InteractionMap();         // IR1 contains V0 information (VIR)
  TBits fIR2 =  event->GetHeader()->GetIRInt2InteractionMap();         // IR2 contains T0 information
  UShort_t bunchCrossings = event->GetBunchCrossNumber();
  if(fHistoPastFutureBits){
    for(Int_t i = 0; i<180;i++){
      if(fIR1.TestBitNumber(i))
        fHistoPastFutureBits->Fill((i*25)-90*25);
    }
  }

  Bool_t isOutOfBunchPileup = 0;
  Int_t pf1 = fPastFutureRejectionLow +bunchCrossings%4;
  Int_t pf2 = fPastFutureRejectionHigh+bunchCrossings%4;
  if(pf1 < -89) pf1 = -89;
  if(pf2 > 89)  pf2 =  89;
  Int_t pf2maxForT0 = pf2;
  Int_t ir1skip     = 0;
  for (Int_t i=pf1;i<=pf2;i++) {
    if (i==0) continue;
    if (i<=pf2maxForT0) isOutOfBunchPileup|=fIR2.TestBitNumber(90+i); // T0-based clean-up
    if (i>0 && i<=ir1skip) continue; // skip next 2 for old IR definitions
    isOutOfBunchPileup|=fIR1.TestBitNumber(90+i); // V0-based clean-up
  }
  return isOutOfBunchPileup;
}

//________________________________________________________________________
AliV0ReaderV1* AliConvEventCuts::GetV0Reader(){

  if(!fV0Reader){
    fV0Reader = dynamic_cast<AliV0ReaderV1*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()));

    if(!fV0Reader){
      AliError("V0Reader could not be obtained, returning nullptr");
    }
  }
  return fV0Reader;
}

//________________________________________________________________________
Double_t AliConvEventCuts::GetV0Multiplicity(AliVEvent *event) const {

  if (!event){
    AliError("event is a nullptr.");
    return -1.;
  }

  if (fIsHeavyIon==2){
      if(event->GetRunNumber()>266400 && event->GetRunNumber()<267140)
        return  event->GetVZEROData()->GetMTotV0C(); // for Pbp
      else
        return event->GetVZEROData()->GetMTotV0A(); // for pPb
  }else{
    return event->GetVZEROData()->GetMTotV0A() + event->GetVZEROData()->GetMTotV0C() ;
  }
}

//________________________________________________________________________
Int_t AliConvEventCuts::GetNumberOfTPCClusters(AliVEvent *event) const {

  if      (dynamic_cast<AliAODEvent*>(event)) return dynamic_cast<AliAODEvent*>(event)->GetNumberOfTPCClusters();
  else if (dynamic_cast<AliESDEvent*>(event)) return dynamic_cast<AliESDEvent*>(event)->GetNumberOfTPCClusters();
  else {
    AliError("event is a nullptr");
    return -1;
  }
}

//________________________________________________________________________
Bool_t AliConvEventCuts::IsPileUpV0MTPCout(AliVEvent *event){

  if (!fFPileUpRejectV0MTPCout){
    AliError("fFPileUpRejectV0MTPCout is a nullptr.");
    return kTRUE;
  }
  Int_t nTracksTPCout = GetV0Reader()->GetNumberOfTPCoutTracks();
  Double_t multV0M    = GetV0Multiplicity(event);

  if (multV0M < fFPileUpRejectV0MTPCout->Eval(nTracksTPCout)){
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::IsPileUpSDDSSDTPC(AliVEvent *event)
{
  if (!fFPileUpRejectSDDSSDTPC){
    AliError("fFPileUpRejectSDDSSDTPC is a nullptr.");
    return kTRUE;
  }
  Int_t nCluSDDSSD = GetV0Reader()->GetSumSDDSSDClusters(event);
  Int_t nCluTPC    = GetNumberOfTPCClusters(event);

  if (nCluSDDSSD <= fFPileUpRejectSDDSSDTPC->Eval(nCluTPC)){
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
void AliConvEventCuts::FillTPCPileUpHistograms(AliVEvent *event){

  if (hV0MultVsNumberTPCoutTracks){
    hV0MultVsNumberTPCoutTracks->Fill(GetV0Reader()->GetNumberOfTPCoutTracks(), GetV0Multiplicity(event));
  }

  if (hTPCSDDSSDClusters){
    hTPCSDDSSDClusters->Fill(GetNumberOfTPCClusters(event), GetV0Reader()->GetSumSDDSSDClusters(event));
  }
}

//________________________________________________________________________
Int_t AliConvEventCuts::GetNumberOfContributorsVtx(AliVEvent *event){
  // returns number of contributors to the vertex

  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(event);
  if(fESDEvent){
    if (fESDEvent->GetPrimaryVertex() != NULL){
      if(fESDEvent->GetPrimaryVertex()->GetNContributors()>0) {
      // cout << "accepted global" << fESDEvent->GetEventNumberInFile() << " with NCont: " << fESDEvent->GetPrimaryVertex()->GetNContributors() << endl;
        return fESDEvent->GetPrimaryVertex()->GetNContributors();
      }
    }

    if(fESDEvent->GetPrimaryVertexSPD() !=NULL){
      if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
      // cout << "accepted SPD" << fESDEvent->GetEventNumberInFile() << " with NCont: " << fESDEvent->GetPrimaryVertexSPD()->GetNContributors() << endl;
        return fESDEvent->GetPrimaryVertexSPD()->GetNContributors();
      }else {
        AliWarning(Form("Number of contributors from bad vertex type:: %s",fESDEvent->GetPrimaryVertex()->GetName()));
            //  cout << "rejected " << fESDEvent->GetEventNumberInFile() << endl;
        return 0;
      }
    }
  }

  AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(event);
  if(fAODEvent){
    if (fAODEvent->GetPrimaryVertex() != NULL){
      if(fAODEvent->GetPrimaryVertex()->GetNContributors()>0) {
        return fAODEvent->GetPrimaryVertex()->GetNContributors();
      }
    }
    if(fAODEvent->GetPrimaryVertexSPD() !=NULL){
      if(fAODEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
        return fAODEvent->GetPrimaryVertexSPD()->GetNContributors();
      } else {
        AliWarning(Form("Number of contributors from bad vertex type:: %s",fAODEvent->GetPrimaryVertex()->GetName()));
        return 0;
      }
    }
  }
  // cout << "rejected " << fESDEvent->GetEventNumberInFile() << endl;
  return 0;
}

//________________________________________________________________________
// Analysing Jet-Jet MC's
//________________________________________________________________________
Bool_t AliConvEventCuts::IsJetJetMCEventAccepted(AliMCEvent *mcEvent, Double_t& weight, Float_t& pthard, AliVEvent* event, Double_t maxJetPt ){
  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;
  weight                               = -1;
  fMaxPtJetMC                          = 0;

  if (  fPeriodEnum != kLHC18b8 && fPeriodEnum != kLHC18b10 && fPeriodEnum != kLHC18l2 &&           // LHC17pq pp 5TeV JetJet MC's
        fPeriodEnum != kLHC19i3b1 && fPeriodEnum != kLHC19i3c1  &&                                  // LHC18 JetJet MC decay EMCal triggered
        fPeriodEnum != kLHC19i3b2 && fPeriodEnum != kLHC19i3c2  &&                                  // LHC18 JetJet MC decay DCal/PHOS triggered
        fPeriodEnum != kLHC18l6b1 && fPeriodEnum != kLHC18l6c1  &&                                  // LHC17 JetJet MC decay EMCal triggered
        fPeriodEnum != kLHC18l6b2 && fPeriodEnum != kLHC18l6c2  &&                                  // LHC17 JetJet MC decay DCal/PHOS triggered
        fPeriodEnum != kLHC17i3a1 &&                                                                // LHC16ijklop GammaJet MC EMCal triggered
        fPeriodEnum != kLHC17i3b1 && fPeriodEnum != kLHC17i3c1 &&                                   // LHC16ijklop JetJet MC EMCal triggered
        fPeriodEnum != kLHC17i3b2 && fPeriodEnum != kLHC17i3c2 &&                                   // LHC16ijklop JetJet MC DCal/PHOS triggered
        fPeriodEnum != kLHC20b1b1 && fPeriodEnum != kLHC20b1c1 &&                                   // LHC16ijklop JetJet MC EMCal triggered new production
        fPeriodEnum != kLHC20b1b2 && fPeriodEnum != kLHC20b1c2 &&                                   // LHC16ijklop JetJet MC DCal/PHOS triggered new production
        fPeriodEnum != kLHC17g8a &&                                                                 // LHC16qt pPb 5TeV JetJet MC's
        fPeriodEnum != kLHC16rP1JJ &&  fPeriodEnum != kLHC16sP1JJ &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC16rsGJ &&                                                                // LHC16sr pPb 8TeV GammaJet MC's
        fPeriodEnum != kLHC17g6b2a &&  fPeriodEnum != kLHC17g6b2b &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC17g6b3a &&  fPeriodEnum != kLHC17g6b3b &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC16P1JJ && fPeriodEnum != kLHC16P1JJLowB &&                               // LHC16X Jet Jet MC's
        fPeriodEnum != kLHC17P1JJ && fPeriodEnum != kLHC17P1JJLowB &&                               // LHC17X Jet Jet MC's
        fPeriodEnum != kLHC18P1JJ &&                                                                // LHC18X Jet Jet MC's
        fPeriodEnum != kLHC16h3  &&                                                                 // LHC15n Jet Jet MC's
        fPeriodEnum != kLHC15a3a && fPeriodEnum != kLHC15a3a_plus && fPeriodEnum != kLHC15a3b &&    // LHC13g Jet Jet MC's
        fPeriodEnum != kLHC15g1a && fPeriodEnum != kLHC15g1b &&                                     // LHC11a Jet Jet MC's
        fPeriodEnum != kLHC13b4_fix && fPeriodEnum != kLHC13b4_plus &&                              // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC19a4 &&                                                                  // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC16c3a && fPeriodEnum != kLHC16c3b && fPeriodEnum != kLHC16c3c &&         // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC17g6a1 && fPeriodEnum != kLHC17g6a2 && fPeriodEnum != kLHC17g6a3 &&      // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC12P2JJ && fPeriodEnum != kLHC17g5b && fPeriodEnum != kLHC17g5c &&        // LHC12 JetJet MC
        fPeriodEnum != kLHC17g5a1 && fPeriodEnum != kLHC17g5a2 &&                                    // LHC12 GammaJet MC
        fPeriodEnum != kLHC14k1a  &&  fPeriodEnum != kLHC14k1b                                      // LHC11 JetJet MC
     ){

    weight = 1;
    return kTRUE;
  }

  if(mcEvent){
    cHeader           = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
    if(cHeader) headerFound   = kTRUE;
  }else{
    //no mcEvent available -> not running on MC
    weight = 1;
    return kTRUE;
  }
  if(headerFound){
    TList *genHeaders         = 0x0;
    if(cHeader) genHeaders    = cHeader->GetHeaders();
    AliGenEventHeader* gh     = 0;
    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName = gh->GetName();
      if (GeneratorName.CompareTo("AliGenPythiaEventHeader") == 0 || GeneratorName.Contains("Pythia8Jets")){
        Bool_t eventAccepted = kTRUE;
        Int_t nTriggerJets = dynamic_cast<AliGenPythiaEventHeader*>(gh)->NTriggerJets();
        Float_t ptHard = dynamic_cast<AliGenPythiaEventHeader*>(gh)->GetPtHard();
        Float_t tmpjet[]={0,0,0,0};
        if(maxJetPt==-1){
          for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
            dynamic_cast<AliGenPythiaEventHeader*>(gh)->TriggerJet(ijet, tmpjet);
            TParticle jet(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
            //Compare jet pT and pt Hard
            if(jet.Pt() > fMaxFacPtHard * ptHard)
              eventAccepted= kFALSE;
            //set highest jet pT
            if (jet.Pt() > fMaxPtJetMC) fMaxPtJetMC = jet.Pt();
          }
        } else {
          fMaxPtJetMC = maxJetPt;
          if(maxJetPt > (fMaxFacPtHard * ptHard)){
            eventAccepted= kFALSE;
          }
        }
        // if minimum jet pT compared to pT hard is required, reject event based on it
        if(fMaxPtJetMC < fMinFacPtHard * ptHard)
          eventAccepted= kFALSE;

        if (mcEvent){
          for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
            AliMCParticle* particle = (AliMCParticle*) mcEvent->GetTrack(i);
            if (!particle) continue;
            // if (TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221){
                if (particle->Pt() > fMaxFacPtHardSingleParticle*ptHard && TMath::Abs(particle->PdgCode()) > 21){
                  eventAccepted= kFALSE;
                }
            // }
          }
        }

        Int_t pthardbin = -1;
        if(fUseFilePathForPthard) pthardbin = GetPtHardBinFromPath(GetV0Reader()->GetCurrentFileName(),event);

        if ( fPeriodEnum == kLHC16P1JJLowB || fPeriodEnum == kLHC17P1JJLowB ) {
          Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                            21, 28, 36, 45, 57,
                                            70, 85, 99, 115, 132,
                                            150, 169, 190, 212, 235,
                                            1000000};
          Double_t weightsBins[20]      = { 43.7553,  13.5848, 6.788, 2.67826, 0.975255,
                                            0.39069, 0.127342, 0.0465597, 0.0206539, 0.00750243,
                                            0.00319118, 0.00122291, 0.000641232, 0.000321437, 0.000168273,
                                            9.17033e-05, 5.34755e-05, 3.01354e-05, 1.74518e-05, 2.8004e-05};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }

        } else  if ( fPeriodEnum == kLHC16P1JJ || fPeriodEnum == kLHC17P1JJ || fPeriodEnum == kLHC18P1JJ ){
          Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                            21, 28, 36, 45, 57,
                                            70, 85, 99, 115, 132,
                                            150, 169, 190, 212, 235,
                                            1000000};
          Double_t weightsBins[20]      = { 43.8654,  13.6215, 6.79856, 2.67526, 0.978794,
                                            0.390797,  0.127769, 0.0465714, 0.0206173, 0.00750282,
                                            0.00318773,  0.00122533, 0.000644385, 0.000321225,  0.00016846,
                                            9.18305e-05, 5.33507e-05, 3.00677e-05, 1.74608e-05, 2.80823e-05};

          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }

        } else  if ( fPeriodEnum == kLHC16h3 ){
          Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                            21, 28, 36, 45, 57,
                                            70, 85, 99, 115, 132,
                                            150, 169, 190, 212, 235,
                                            1000000};
          Double_t weightsBins[20]      = { 16.0869, 4.61169, 2.14976, 0.782544, 0.264854,
                                            9.7619E-02, 2.92747E-02, 9.89515E-03, 4.05152E-03, 1.35393E-03,
                                            5.29864E-04, 1.88317E-04, 9.23E-05, 4.29E-05, 2.09E-05,
                                            1.06E-05, 5.76E-06, 3.00E-06, 1.62E-06, 2.10E-06 };
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }


        } else if ( fPeriodEnum == kLHC14k1a ){
          Double_t ptHardBinRanges[7]  = { 5,  7,  9, 12, 16,
                                            21, 1000000};
          Double_t weightsBins[6]      = { 2.327372e-02, 1.783327e-02, 1.678043e-02, 1.167544e-02, 7.066289e-03, 8.714857e-03};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC14k1b ){
          Double_t ptHardBinRanges[8]  = { 10,  14,  19, 26, 35,
                                            48, 66, 1000000};
          Double_t weightsBins[7]      = {  6.174824e-04 ,  6.557521e-04 ,  6.472503e-04 ,  4.857432e-04 ,  3.402152e-04 ,  1.873434e-04 ,  1.376054e-04 };
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 7) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC15a3b || fPeriodEnum == kLHC15g1b ){
          Double_t ptHardBinRanges[13]  = { 5,  7,  9, 12, 16,
                                            21, 28, 36, 45, 57,
                                            70, 85, 1000};
          Double_t weightsBins[12]      = { 7.858393e-03, 4.718691e-03, 4.077575e-03, 2.814527e-03, 1.669625e-03,
                                            1.007535e-03, 4.536554e-04, 2.111041e-04, 1.094840e-04, 4.404973e-05,
                                            1.933238e-05, 1.562895e-05};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 12) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC15g1a ){
          Double_t ptHardBinRanges[20]  = { 5,      11,     21,     36,      57,
                                            84,     117,    152,    191, 1000000,
                                            5,       7,      9,     12,      16,
                                            21,      28,     36,     45,      57 };
          Double_t weightsBins[19]      = { 4.43629 , 0.49523, 0.0394921, 0.00383174, 0.000446559,
                                            6.37374e-05, 1.03134e-05, 2.27012e-06, 7.59281e-07, 0,
                                            2.62906, 1.12884, 0.656873, 0.262822,  0.0876732,
                                            0.0307759, 0.0087083, 0.0027664, 0.00106203};

          Int_t bin = 0;
          Int_t binFromFile = GetV0Reader()->GetPtHardFromFile();
          if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 19) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC15a3a || fPeriodEnum == kLHC15a3a_plus ) {
          Double_t ptHardBinRanges[20]  = { 5,      11,     21,     36,      57,
                                            84,     117,    152,    191, 1000000,
                                            5,       7,      9,     12,      16,
                                            21,      28,     36,     45,      57 };
          // LHC15a3a
          Double_t weightsBins[19]      = { 4.43897 , 0.495766, 0.039486, 0.00383011, 0.000447104,
                                            6.37277e-05, 1.03166e-05, 2.26971e-06, 7.59023e-07, 0,
                                            2.63331, 1.12815, 0.657034, 0.262756,  0.0877227,
                                            0.0307638, 0.00870635, 0.00276658, 0.00106229};
          Int_t bin = 0;
          Int_t binFromFile = GetV0Reader()->GetPtHardFromFile();
          if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 19) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC18b10 ){
          Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
          Double_t weightsBins[6]      = { 0.000103227,      1.31851e-05,     1.94129e-06,     3.26392e-07,    6.51801e-08,
                                            2.21123e-08 };

          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC18l2 ){
          Double_t ptHardBinRanges[13]  = { 5, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 10000};
          Double_t weightsBins[12]      = { 0.000103227, 0.000221821, 0.00039017, 0.000503916, 0.000498568,
                                            0.000463416, 0.000304633, 0.00018642, 0.000127846, 6.50537e-05,
                                            3.63397e-05, 1.70296e-05 };

          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 12) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)) weight = weightsBins[pthardbin-1];
        } else if ( fPeriodEnum == kLHC12P2JJ ){
          Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                            21, 28, 36, 45, 57,
                                            70, 85, 99, 115, 132,
                                            150, 169, 190, 212, 235,
                                            1000000};
          Double_t weightsBins[20]      = { 28.3084, 8.43277, 4.07753, 1.54359, 0.543318,
                                            0.208394, 0.0652349, 0.0186904, 0.00834528, 0.00301414,
                                            0.00125939, 0.000474403, 0.000244052, 0.00011924, 6.09838e-05,
                                            3.24148e-05, 1.84314e-05, 1.00926e-05, 5.68632e-06, 8.38092e-06};
            Int_t bin = 0;
            while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
            if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC17g5a1 ||fPeriodEnum == kLHC17g5a2 ){
          Double_t ptHardBinRanges[7]  = { 5,  11,  21, 36, 57,
                                            84, 1000000};
          Double_t weightsBins[6]      = {1.43277693e-04, 1.88432909e-05, 2.86381554e-06, 5.18655146e-07, 1.08660208e-07, 3.89688100e-08};
            Int_t bin = 0;
            while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
            if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }

        } else if ( fPeriodEnum == kLHC16c3a ){ // ALIROOT-5901
          Double_t ptHardBinRanges[6]   = {  7, 9, 12, 16, 21, 1000};
          Double_t weightsBins[5]       = {  6.731200e-03, 7.995602e-03, 6.778717e-03, 4.643571e-03, 6.014497e-03};
          Int_t bin = 0;
          pthardbin++;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 5) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if (fPeriodEnum == kLHC16c3b ){ // ALIROOT-5901
          Double_t ptHardBinRanges[7]   = {  14, 19, 26, 35, 48, 66, 1000};
          Double_t weightsBins[6]       = {  6.07559700e-03, 3.93844754e-03, 2.00223387e-03, 9.85233093e-04, 3.89161623e-04, 1.86561078e-04};
          Int_t bin = 0;
          pthardbin++;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if (fPeriodEnum == kLHC16c3c ){ // ALIROOT-5901
          Double_t ptHardBinRanges[8]   = {  0, 5, 11, 21, 36, 57, 84, 1000};
          Double_t weightsBins[7]       = {  0.00151999, 0.000100346, 1.27688e-05, 1.82388e-06, 3.08506e-07, 6.00308e-08, 1.88414e-08}; //preliminary estimates
          Int_t bin = 0;
          pthardbin++;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 7) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC17g6a1 ){ //ALIROOT-7271 - GJ
          Double_t ptHardBinRanges[7]   = {  5, 11, 21, 36, 57, 84, 1000};
          Double_t weightsBins[6]       = {  6.731200e-03, 6.731200e-03, 7.995602e-03, 6.778717e-03, 4.643571e-03, 6.014497e-03};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        } else if (fPeriodEnum == kLHC17g6a2 ){ //ALIROOT-7271 - JJ low EMC trigg
          Double_t ptHardBinRanges[7]   = {  5, 7, 9, 12, 16, 21, 1000 };
          Double_t weightsBins[6]       = {  1.31493822e-02, 9.76893424e-03, 9.52744251e-03, 7.12499226e-03, 4.53465002e-03, 5.66990089e-03};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        } else if (fPeriodEnum == kLHC17g6a3 ){ //ALIROOT-7271 - JJ high EMC trigg
          Double_t ptHardBinRanges[9]   = {  8, 10, 14, 19, 26, 35, 48, 66, 1000};
          Double_t weightsBins[8]       = {  2.17156487e-04, 4.81510348e-04, 5.17290689e-04, 4.92575600e-04, 3.54346740e-04, 2.37913195e-04, 1.25136140e-04, 8.48443967e-05};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        } else if ( fPeriodEnum == kLHC13b4_fix || fPeriodEnum == kLHC13b4_plus ){
          Double_t ptHardBinRanges[11]  = { 5,     11,   21,   36,   57,
                                            84,    117,   152,  191,    234,
                                            1000};
          Double_t weightsBins[10]      = { 2.24185e-6 , 2.48463e-7, 2.23171e-8, 2.43667e-9, 3.29934e-10,
                                            5.34592e-11, 1.00937e-11, 2.6493e-12, 8.53912e-13, 5.43077e-13};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 10) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC19a4 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
          Double_t weightsBins[20]      = {17.8908, 4.99797, 2.31865, 0.824485, 0.276159, 0.100522, 0.0299966, 0.00994219, 0.00405425, 0.00134633, 0.00052021, 0.000182129, 8.93822e-05, 4.07073e-05, 1.97037e-05, 9.97165e-06, 5.33217e-06, 2.79808e-06, 1.47593e-06, 1.85158e-06};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC17g8a ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
          Double_t weightsBins[20]      = {1.565930E+01, 4.598350E+00, 2.081240E+00, 7.744650E-01, 2.644240E-01, 1.002330E-01, 2.979190E-02, 9.696490E-03, 3.950930E-03, 1.333040E-03, 5.210630E-04, 1.927180E-04, 9.235930E-05, 4.346820E-05, 2.120660E-05, 1.073260E-05, 5.701210E-06, 3.047490E-06, 1.664780E-06, 2.123400E-06};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC16rP1JJ ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
          Double_t weightsBins[20]      = {2.723740E+01, 8.127540E+00, 3.934400E+00, 1.492720E+00, 5.268010E-01, 2.033790E-01, 6.361520E-02, 2.256080E-02, 9.638840E-03, 3.372890E-03, 1.381980E-03, 5.121390E-04, 2.613120E-04, 1.260940E-04, 6.393150E-05, 3.386080E-05, 1.926040E-05, 1.046950E-05, 5.895950E-06, 8.658420E-06};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC16sP1JJ ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
          Double_t weightsBins[20]      = {2.716550E+01, 8.121430E+00, 3.932100E+00, 1.492830E+00, 5.272190E-01, 2.023090E-01, 6.371860E-02, 2.245360E-02, 9.590340E-03, 3.369300E-03, 1.384470E-03, 5.119390E-04, 2.606910E-04, 1.259110E-04, 6.408650E-05, 3.396290E-05, 1.917340E-05, 1.044610E-05, 5.882680E-06, 8.672390E-06};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC16rsGJ ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
          Double_t weightsBins[6]      = {1.54182815e-04, 2.03902180e-05, 3.23031581e-06, 5.63579758e-07,1.16970192e-07,4.27062629e-08};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC17i3a1 ){ // weights obtained from ga_pp_mc_aod train 912
           Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
           Double_t weightsBins[6]      = { 0.0002181, 3.13684e-05, 5.01515e-06, 9.50662e-07, 2.08186e-07, 7.96555e-08};
           Int_t bin = 0;
           while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
           if (bin < 6) weight = weightsBins[bin];
           if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC17i3b1 ){ // preliminary weights obtained from local running
           Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
           Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
           Int_t bin = 0;
           while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
           if (bin < 8) weight = weightsBins[bin];
           if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC17i3b2 ){ // preliminary weights obtained from local running
           Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
           Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
           Int_t bin = 0;
           while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
           if (bin < 8) weight = weightsBins[bin];
           if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC17i3c1 ){ // preliminary weights obtained from local running
           Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
           Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
           Int_t bin = 0;
           while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
           if (bin < 8) weight = weightsBins[bin];
           if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC17i3c2 ){ // preliminary weights obtained from local running
           Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
           Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
           Int_t bin = 0;
           while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
           if (bin < 8) weight = weightsBins[bin];
           if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC20b1b1 ){ // preliminary weights obtained from local running
           Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 22, 10000};
           Double_t weightsBins[6]      = { 0.0448228, 0.0388829, 0.0366336, 0.0278004, 0.0175832, 0.0241481};
           Int_t bin = 0;
           while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
           if (bin < 6) weight = weightsBins[bin];
           if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC20b1b2 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 22, 10000};
          Double_t weightsBins[6]      = { 0.0318163, 0.0280099, 0.0271878, 0.0208271, 0.0130291, 0.0181411};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
             weight = weightsBins[pthardbin-1];
             if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
           }
        } else if ( fPeriodEnum == kLHC20b1c1 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
          Double_t weightsBins[8]      = { 0.000771445, 0.00170422, 0.00183904, 0.00185169, 0.0014099, 0.00103708, 0.000604256, 0.000494368};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC20b1c2 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
          Double_t weightsBins[8]      = { 0.000528778, 0.00121486, 0.00132257, 0.00132752, 0.00102868, 0.000741674, 0.000437791, 0.000350834};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC18l6c1 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
          Double_t weightsBins[8]      = { 0.000802655, 0.00171876 ,0.00185997, 0.00184664 , 0.00140638  , 0.00102373  , 0.000599643 , 0.000496501 };
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
      } else if ( fPeriodEnum == kLHC18l6c2 ){ // preliminary weights obtained from local running
         Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
         Double_t weightsBins[8]      = { 0.000583836, 0.00124525 ,0.00133701, 0.00132045 , 0.0010136  , 0.000738478  , 0.00043112 , 0.00035756 };
         Int_t bin = 0;
         while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
         if (bin < 8) weight = weightsBins[bin];
         if(fUseFilePathForPthard && (pthardbin > -1)){
           weight = weightsBins[pthardbin-1];
           if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
         }
      } else if ( fPeriodEnum == kLHC18l6b1 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
          Double_t weightsBins[6]      = { 0.0479994 , 0.039711  ,0.039082 , 0.0287247 , 0.0182263  , 0.0248411};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
       } else if ( fPeriodEnum == kLHC18l6b2 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
          Double_t weightsBins[6]      = { 0.0343221 , 0.0287842  ,0.0282014 , 0.0207503 , 0.0132192  , 0.0180455};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
       } else if ( fPeriodEnum == kLHC19i3b1 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
          Double_t weightsBins[6]      = { 0.0479994 , 0.039711  ,0.039082 , 0.0287247 , 0.0182263  , 0.0248411};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC19i3b2 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
          Double_t weightsBins[6]      = { 0.0343118 , 0.0287709  ,0.028205 , 0.0207509 , 0.0132104  , 0.0180347};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC19i3c1 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
          Double_t weightsBins[8]      = { 0.000810919, 0.0017328 ,0.00185665, 0.00183507 , 0.00140593  , 0.00102492  , 0.000597323 , 0.000494353 };
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC19i3c2 ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
          Double_t weightsBins[8]      = { 0.000583836, 0.00124525 ,0.00133701, 0.00132045 , 0.0010136  , 0.000738478  , 0.00043112 , 0.00035756 };
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC17g5b ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[7]  = { 5, 7, 9, 12, 16, 21, 10000};
          Double_t weightsBins[6]  = { 0.0235115, 0.0186173, 0.0182963, 0.0135686, 0.00874758, 0.0115796};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else if ( fPeriodEnum == kLHC17g5c ){ // preliminary weights obtained from local running
          Double_t ptHardBinRanges[9]  = { 8, 10, 14, 19, 26, 35, 48, 66, 10000};
          Double_t weightsBins[8]  = { 0.00039861, 0.000868597, 0.000939155, 0.000917328, 0.00068357, 0.000475878, 0.000266069, 0.000203764};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
      } else if ( fPeriodEnum == kLHC17g6b2a ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = { 5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]  = { 0.0267755, 0.0206222, 0.0199551, 0.0149312, 0.00936012, 0.0124922};
        Int_t bin = 0;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 6) weight = weightsBins[bin];
        if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
      } else if ( fPeriodEnum == kLHC17g6b2b ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = { 5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]  = { 0.0186043, 0.0144098, 0.0145382, 0.0106785, 0.00688777, 0.00923482};
        Int_t bin = 0;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 6) weight = weightsBins[bin];
        if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
      } else if ( fPeriodEnum == kLHC17g6b3a ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = { 8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]  = { 0.000460503, 0.000942402, 0.00101341, 0.00103905, 0.000734619, 0.000520764, 0.000297075, 0.000217233};
        Int_t bin = 0;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 8) weight = weightsBins[bin];
        if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
      } else if ( fPeriodEnum == kLHC17g6b3b ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = { 8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]  = { 0.000305143, 0.000672797, 0.000737654, 0.0007214, 0.000543925, 0.000377044, 0.000214699, 0.00015541};
        Int_t bin = 0;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 8) weight = weightsBins[bin];
        if(fUseFilePathForPthard && (pthardbin > -1)){
            weight = weightsBins[pthardbin-1];
            if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
          }
        } else {
          weight = 1;
        }
        if(isnan(weight) || weight == 0){
          TString fFileNameBroken =  ((TString)GetV0Reader()->GetCurrentFileName()).Data();
          TString debugMessage=Form("JJ weight = %.05f for file: %s \n event Nr in File: %i \n CutNr: %s", weight, fFileNameBroken.Data(), mcEvent->GetEventNumberInFile(), fCutStringRead.Data());
          AliFatal(debugMessage.Data());
        }
        if (weight == -1) return kFALSE;
        else return eventAccepted;

      }
    }
  } else {
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    TString eventHeaderName     = eventHeader->ClassName();
    Bool_t eventAccepted = kFALSE;
    if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0 || eventHeaderName.Contains("Pythia8Jets")){
      eventAccepted = kTRUE;
    }else { //special case for pythia8jets embedded in EPOSLHC for AODs
      if(event->IsA()==AliAODEvent::Class()){
        AliAODMCHeader *mch = NULL;
        AliAODEvent * aod = dynamic_cast<AliAODEvent*> (event);
        if(aod){
          mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
          if ( mch ){
            Int_t nGenerators = mch->GetNCocktailHeaders();
            if ( nGenerators > 0  ){
              for(Int_t igen = 0; igen < nGenerators; igen++)
              {
                AliGenEventHeader * eventHeaderGen = mch->GetCocktailHeader(igen) ;
                TString name = eventHeaderGen->GetName();
                if (name.CompareTo("AliGenPythiaEventHeader") == 0 || name.Contains("Pythia8Jets") || name.Contains("Pythia8GammaJet")){
                  eventAccepted = kTRUE;
                  eventHeader = eventHeaderGen;
                }
              }
            }
          }
        }
      }
    }
    if(eventAccepted){
      Int_t nTriggerJets =  dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->NTriggerJets();
      Float_t ptHard = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->GetPtHard();
      pthard = ptHard;
      Float_t tmpjet[]={0,0,0,0};
      if(maxJetPt==-1){
        for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
          dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->TriggerJet(ijet, tmpjet);
          TParticle jet(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
          //Compare jet pT and pt Hard
          if(jet.Pt() > fMaxFacPtHard * ptHard){
            eventAccepted= kFALSE;
          }
          //set highest jet pT
          if (jet.Pt() > fMaxPtJetMC){
            fMaxPtJetMC = jet.Pt();
          }
        }
      } else {
        fMaxPtJetMC = maxJetPt;
        if(maxJetPt > (fMaxFacPtHard * ptHard)){
          eventAccepted= kFALSE;
        }
      }
      // if minimum jet pT compared to pT hard is required, reject event based on it
      if(fMaxPtJetMC < fMinFacPtHard * ptHard){
        eventAccepted= kFALSE;
      }

      if (mcEvent){
        for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
          // TParticle* particle = (TParticle *)mcEvent->Particle(i);
          AliMCParticle* particle = (AliMCParticle*) mcEvent->GetTrack(i);
          if (!particle) continue;
          // if (TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221){
              if (particle->Pt() > fMaxFacPtHardSingleParticle*ptHard && TMath::Abs(particle->PdgCode()) > 21){
                eventAccepted= kFALSE;
              }
          // }
        }
      }
      Int_t pthardbin = -1;
      if(fUseFilePathForPthard) pthardbin = GetPtHardBinFromPath(GetV0Reader()->GetCurrentFileName(),event);

      if ( fPeriodEnum == kLHC16P1JJLowB || fPeriodEnum == kLHC16P1JJLowB) {
        Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                          21, 28, 36, 45, 57,
                                          70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
        Double_t weightsBins[20]      = { 43.7553,  13.5848, 6.788, 2.67826, 0.975255,
                                          0.39069, 0.127342, 0.0465597, 0.0206539, 0.00750243,
                                          0.00319118, 0.00122291, 0.000641232, 0.000321437, 0.000168273,
                                          9.17033e-05, 5.34755e-05, 3.01354e-05, 1.74518e-05, 2.8004e-05};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else  if ( fPeriodEnum == kLHC16P1JJ || fPeriodEnum == kLHC17P1JJ || fPeriodEnum == kLHC18P1JJ){
        Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                          21, 28, 36, 45, 57,
                                          70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
        Double_t weightsBins[20]      = { 43.8654,  13.6215, 6.79856, 2.67526, 0.978794,
                                          0.390797,  0.127769, 0.0465714, 0.0206173, 0.00750282,
                                          0.00318773,  0.00122533, 0.000644385, 0.000321225,  0.00016846,
                                          9.18305e-05, 5.33507e-05, 3.00677e-05, 1.74608e-05, 2.80823e-05};

        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else  if ( fPeriodEnum == kLHC16h3 ){
        Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                          21, 28, 36, 45, 57,
                                          70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
        Double_t weightsBins[20]      = { 16.0869, 4.61169, 2.14976, 0.782544, 0.264854,
                                          9.7619E-02, 2.92747E-02, 9.89515E-03, 4.05152E-03, 1.35393E-03,
                                          5.29864E-04, 1.88317E-04, 9.23E-05, 4.29E-05, 2.09E-05,
                                          1.06E-05, 5.76E-06, 3.00E-06, 1.62E-06, 2.10E-06 };
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC15a3b || fPeriodEnum == kLHC15g1b ){
        Double_t ptHardBinRanges[13]  = { 5,  7,  9, 12, 16,
                                          21, 28, 36, 45, 57,
                                          70, 85, 1000};
        Double_t weightsBins[12]      = { 7.858393e-03, 4.718691e-03, 4.077575e-03, 2.814527e-03, 1.669625e-03,
                                          1.007535e-03, 4.536554e-04, 2.111041e-04, 1.094840e-04, 4.404973e-05,
                                          1.933238e-05, 1.562895e-05};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 12) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC15g1a ){
        Double_t ptHardBinRanges[20]  = { 5,      11,     21,     36,      57,
                                          84,     117,    152,    191, 1000000,
                                          5,       7,      9,     12,      16,
                                          21,      28,     36,     45,      57 };
        Double_t weightsBins[19]      = { 4.43629 , 0.49523, 0.0394921, 0.00383174, 0.000446559,
                                          6.37374e-05, 1.03134e-05, 2.27012e-06, 7.59281e-07, 0,
                                          2.62906, 1.12884, 0.656873, 0.262822,  0.0876732,
                                          0.0307759, 0.0087083, 0.0027664, 0.00106203};

        Int_t bin = 0;
        Int_t binFromFile = GetV0Reader()->GetPtHardFromFile();
        if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 19) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC15a3a || fPeriodEnum == kLHC15a3a_plus ) {
        Double_t ptHardBinRanges[20]  = { 5,      11,     21,     36,      57,
                                          84,     117,    152,    191, 1000000,
                                          5,       7,      9,     12,      16,
                                          21,      28,     36,     45,      57 };
        // LHC15a3a
        Double_t weightsBins[19]      = { 4.43897 , 0.495766, 0.039486, 0.00383011, 0.000447104,
                                          6.37277e-05, 1.03166e-05, 2.26971e-06, 7.59023e-07, 0,
                                          2.63331, 1.12815, 0.657034, 0.262756,  0.0877227,
                                          0.0307638, 0.00870635, 0.00276658, 0.00106229};
        Int_t bin = 0;
        Int_t binFromFile = GetV0Reader()->GetPtHardFromFile();
        if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 19) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC12P2JJ ){
        Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                          21, 28, 36, 45, 57,
                                          70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
        Double_t weightsBins[20]      = { 28.3084, 8.43277, 4.07753, 1.54359, 0.543318,
                                          0.208394, 0.0652349, 0.0186904, 0.00834528, 0.00301414,
                                          0.00125939, 0.000474403, 0.000244052, 0.00011924, 6.09838e-05,
                                          3.24148e-05, 1.84314e-05, 1.00926e-05, 5.68632e-06, 8.38092e-06};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }

      } else if ( fPeriodEnum == kLHC17g5a1 || fPeriodEnum == kLHC17g5a2 ){
        Double_t ptHardBinRanges[7]  = { 5,  11,  21, 36, 57,
                                            84, 1000000};
        Double_t weightsBins[6]      = {1.43277693e-04, 1.88432909e-05, 2.86381554e-06, 5.18655146e-07, 1.08660208e-07, 3.89688100e-08};
            Int_t bin = 0;
            if(ptHard >= ptHardBinRanges[0]){
              while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
              if (bin < 6) weight = weightsBins[bin];
              if(fUseFilePathForPthard && (pthardbin > -1)){
                weight = weightsBins[pthardbin-1];
                if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
              }
            } else{ // smaller than lowest border
               weight = weightsBins[0];
               if(fUseFilePathForPthard && (pthardbin > -1)){
                weight = weightsBins[pthardbin-1];
                if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
               }
            }
      } else if ( fPeriodEnum == kLHC16c3a ){ //ALIROOT-5901
        Double_t ptHardBinRanges[6]   = {  7, 9, 12, 16, 21, 1000};
        Double_t weightsBins[5]       = {  6.73298726e-03, 8.00549934e-03, 6.77989565e-03, 4.64169953e-03, 6.01322269e-03};
        Int_t bin = 0;
        pthardbin++;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 5) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if (fPeriodEnum == kLHC16c3b ){ //ALIROOT-5901
        Double_t ptHardBinRanges[7]   = {  14, 19, 26, 35, 48, 66, 1000};
        Double_t weightsBins[6]       = {  6.07559700e-03, 3.93844754e-03, 2.00223387e-03, 9.85233093e-04, 3.89161623e-04, 1.86561078e-04};
        Int_t bin = 0;
        pthardbin++;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if (fPeriodEnum == kLHC16c3c ){ //ALIROOT-5901
        Double_t ptHardBinRanges[8]   = {  0, 5, 11, 21, 36, 57, 84, 1000};
        Double_t weightsBins[7]       = {  0.00151999, 0.000100346, 1.27688e-05, 1.82388e-06, 3.08506e-07, 6.00308e-08, 1.88414e-08}; //preliminary estimates
        Int_t bin = 0;
        pthardbin++;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 7) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g6a1 ){ //ALIROOT-7271 - GJ
        Double_t ptHardBinRanges[7]   = {  5, 11, 21, 36, 57, 84, 1000};
        Double_t weightsBins[6]       = {  6.731200e-03, 6.731200e-03, 7.995602e-03, 6.778717e-03, 4.643571e-03, 6.014497e-03};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if (fPeriodEnum == kLHC17g6a2 ){ //ALIROOT-7271 - JJ low EMC trigg
        Double_t ptHardBinRanges[7]   = {  5, 7, 9, 12, 16, 21, 1000 };
        Double_t weightsBins[6]       = {  1.31493822e-02, 9.76893424e-03, 9.52744251e-03, 7.12499226e-03, 4.53465002e-03, 5.66990089e-03};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if (fPeriodEnum == kLHC17g6a3 ){ //ALIROOT-7271 - JJ high EMC trigg
        Double_t ptHardBinRanges[9]   = {  8, 10, 14, 19, 26, 35, 48, 66, 1000};
        Double_t weightsBins[8]       = {  2.17156487e-04, 4.81510348e-04, 5.17290689e-04, 4.92575600e-04, 3.54346740e-04, 2.37913195e-04, 1.25136140e-04, 8.48443967e-05};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC14k1a ){
        Double_t ptHardBinRanges[7]  = { 5,  7,  9, 12, 16,
                                          21, 1000000};
        Double_t weightsBins[6]      = { 2.327372e-02, 1.783327e-02, 1.678043e-02, 1.167544e-02, 7.066289e-03, 8.714857e-03};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC14k1b ){
        Double_t ptHardBinRanges[8]  = { 10,  14,  19, 26, 35,
                                          48, 66, 1000000};
        Double_t weightsBins[7]      = {  6.174824e-04 ,  6.557521e-04 ,  6.472503e-04 ,  4.857432e-04 ,  3.402152e-04 ,  1.873434e-04 ,  1.376054e-04 };
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 7) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC13b4_fix || fPeriodEnum == kLHC13b4_plus ){
        Double_t ptHardBinRanges[11]  = { 5,     11,   21,   36,   57,
                                          84,    117,   152,  191,    234,
                                          1000};
        Double_t weightsBins[10]      = { 2.24185e-6 , 2.48463e-7, 2.23171e-8, 2.43667e-9, 3.29934e-10,
                                          5.34592e-11, 1.00937e-11, 2.6493e-12, 8.53912e-13, 5.43077e-13};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 10) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC19a4 ){
        Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
        Double_t weightsBins[20]      = {17.8908, 4.99797, 2.31865, 0.824485, 0.276159, 0.100522, 0.0299966, 0.00994219, 0.00405425, 0.00134633, 0.00052021, 0.000182129, 8.93822e-05, 4.07073e-05, 1.97037e-05, 9.97165e-06, 5.33217e-06, 2.79808e-06, 1.47593e-06, 1.85158e-06};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
    } else if ( fPeriodEnum == kLHC17g8a ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
        Double_t weightsBins[20]      = {1.565930E+01, 4.598350E+00, 2.081240E+00, 7.744650E-01, 2.644240E-01, 1.002330E-01, 2.979190E-02, 9.696490E-03, 3.950930E-03, 1.333040E-03, 5.210630E-04, 1.927180E-04, 9.235930E-05, 4.346820E-05, 2.120660E-05, 1.073260E-05, 5.701210E-06, 3.047490E-06, 1.664780E-06, 2.123400E-06};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
    } else if ( fPeriodEnum == kLHC16rP1JJ ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
        Double_t weightsBins[20]      = {2.723740E+01, 8.127540E+00, 3.934400E+00, 1.492720E+00, 5.268010E-01, 2.033790E-01, 6.361520E-02, 2.256080E-02, 9.638840E-03, 3.372890E-03, 1.381980E-03, 5.121390E-04, 2.613120E-04, 1.260940E-04, 6.393150E-05, 3.386080E-05, 1.926040E-05, 1.046950E-05, 5.895950E-06, 8.658420E-06};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
    } else if ( fPeriodEnum == kLHC16sP1JJ ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
        Double_t weightsBins[20]      = {2.716550E+01, 8.121430E+00, 3.932100E+00, 1.492830E+00, 5.272190E-01, 2.023090E-01, 6.371860E-02, 2.245360E-02, 9.590340E-03, 3.369300E-03, 1.384470E-03, 5.119390E-04, 2.606910E-04, 1.259110E-04, 6.408650E-05, 3.396290E-05, 1.917340E-05, 1.044610E-05, 5.882680E-06, 8.672390E-06};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
    } else if ( fPeriodEnum == kLHC16rsGJ ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
        Double_t weightsBins[6]      = {1.54182815e-04, 2.03902180e-05, 3.23031581e-06, 5.63579758e-07,1.16970192e-07,4.27062629e-08};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17i3a1 ){ // weights obtained from ga_pp_mc_aod train 912
        Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
        Double_t weightsBins[6]      = { 0.0002181, 3.13684e-05, 5.01515e-06, 9.50662e-07, 2.08186e-07, 7.96555e-08};
         Int_t bin = 0;
         if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
         }
      } else if ( fPeriodEnum == kLHC17i3b1 ){ // preliminary weights obtained from local running
         Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
         Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
         Int_t bin = 0;
         if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17i3b2 ){ // preliminary weights obtained from local running
         Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
         Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
         Int_t bin = 0;
         if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17i3c1 ){ // preliminary weights obtained from local running
         Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
         Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
         Int_t bin = 0;
         if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17i3c2 ){ // preliminary weights obtained from local running
         Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
         Double_t weightsBins[8]      = { 0.000813592, 0.00172074, 0.00187963, 0.00184331, 0.00142672, 0.0010083, 0.000599846, 0.000499877};
         Int_t bin = 0;
         if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC20b1b1 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 22, 10000};
        Double_t weightsBins[6]      = { 0.0448228, 0.0388829, 0.0366336, 0.0278004, 0.0175832, 0.0241481};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC20b1b2 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 22, 10000};
        Double_t weightsBins[6]      = { 0.0318163, 0.0280099, 0.0271878, 0.0208271, 0.0130291, 0.0181411};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC20b1c1 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]      = { 0.000771445, 0.00170422, 0.00183904, 0.00185169, 0.0014099, 0.00103708, 0.000604256, 0.000494368};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC20b1c2 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]      = { 0.000528778, 0.00121486, 0.00132257, 0.00132752, 0.00102868, 0.000741674, 0.000437791, 0.000350834};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18b8 ){
        Double_t ptHardBinRanges[21]  = { 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000};
        Double_t weightsBins[20]      = { 16.1083,      4.60917,     2.15196,     0.782021,    0.26541,
                                           0.0978374,   0.0294286,   0.00989457,  0.0040615,   0.00135787,
                                           0.000531766, 0.000188772, 9.23331e-05, 4.30245e-05, 2.10196e-05,
                                           1.06695e-05, 5.78742e-06, 3.02897e-06, 1.62702e-06, 2.12118e-06 };

        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18b10 ){
        Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
        Double_t weightsBins[6]      = { 0.000103227,      1.31851e-05,     1.94129e-06,     3.26392e-07,    6.51801e-08,
                                           2.21123e-08 };

        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18l2 ){
        Double_t ptHardBinRanges[7]  = { 5, 11, 21, 36, 57, 84, 10000};
        Double_t weightsBins[6]      = { 0.000103227,      1.31851e-05,     1.94129e-06,     3.26392e-07,    6.51801e-08,
                                           2.21123e-08 };

        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18l6c1 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]      = { 0.000802655, 0.00171876 ,0.00185997, 0.00184664 , 0.00140638  , 0.00102373  , 0.000599643 , 0.000496501 };
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18l6c2 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]      = { 0.000583836, 0.00124525 ,0.00133701, 0.00132045 , 0.0010136  , 0.000738478  , 0.00043112 , 0.00035756 };
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18l6b1 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]      = { 0.0479994 , 0.039711  ,0.039082 , 0.0287247 , 0.0182263  , 0.0248411};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC18l6b2 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]      = { 0.0343221 , 0.0287842  ,0.0282014 , 0.0207503 , 0.0132192  , 0.0180455};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC19i3b1 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]      = { 0.0479994 , 0.039711  ,0.039082 , 0.0287247 , 0.0182263  , 0.0248411};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC19i3b2 ){ // preliminary weights obtained from local running
         Double_t ptHardBinRanges[7]  = {  5, 7, 9, 12, 16, 21, 10000};
         Double_t weightsBins[6]      = { 0.0343118 , 0.0287709  ,0.028205 , 0.0207509 , 0.0132104  , 0.0180347};
         Int_t bin = 0;
         if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC19i3c1 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]      = { 0.000810919, 0.0017328 ,0.00185665, 0.00183507 , 0.00140593  , 0.00102492  , 0.000597323 , 0.000494353 };
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC19i3c2 ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = {  8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]      = { 0.000583836, 0.00124525 ,0.00133701, 0.00132045 , 0.0010136  , 0.000738478  , 0.00043112 , 0.00035756 };
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g5b ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = { 5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]  = { 0.0235115, 0.0186173, 0.0182963, 0.0135686, 0.00874758, 0.0115796};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g5c ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = { 8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]  = { 0.00039861, 0.000868597, 0.000939155, 0.000917328, 0.00068357, 0.000475878, 0.000266069, 0.000203764};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g6b2a ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = { 5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]  = { 0.0267755, 0.0206222, 0.0199551, 0.0149312, 0.00936012, 0.0124922};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g6b2b ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[7]  = { 5, 7, 9, 12, 16, 21, 10000};
        Double_t weightsBins[6]  = { 0.0186043, 0.0144098, 0.0145382, 0.0106785, 0.00688777, 0.00923482};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g6b3a ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = { 8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]  = { 0.000460503, 0.000942402, 0.00101341, 0.00103905, 0.000734619, 0.000520764, 0.000297075, 0.000217233};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
      } else if ( fPeriodEnum == kLHC17g6b3b ){ // preliminary weights obtained from local running
        Double_t ptHardBinRanges[9]  = { 8, 10, 14, 19, 26, 35, 48, 66, 10000};
        Double_t weightsBins[8]  = { 0.000305143, 0.000672797, 0.000737654, 0.0007214, 0.000543925, 0.000377044, 0.000214699, 0.00015541};
        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 8) weight = weightsBins[bin];
          if(fUseFilePathForPthard && (pthardbin > -1)){
              weight = weightsBins[pthardbin-1];
              if(fUseAdditionalOutlierRejection && ((ptHard < ptHardBinRanges[pthardbin-1]) || (ptHard > ptHardBinRanges[pthardbin]))) eventAccepted= kFALSE;
            }
        }
    } else {
      weight = 1;
    }

    if(isnan(weight) || weight == 0){
      TString fFileNameBroken =  ((TString)GetV0Reader()->GetCurrentFileName()).Data();
      TString debugMessage=Form("JJ weight = %.05f for file: %s \n event Nr in File: %i \n CutNr: %s", weight, fFileNameBroken.Data(), mcEvent->GetEventNumberInFile(), fCutStringRead.Data());
      AliFatal(debugMessage.Data());
    }
    if (weight == -1) return kFALSE;
    else return eventAccepted;

    } else {
      return kFALSE;
    }
  }

  return kFALSE;
}

//________________________________________________________________________
// Analysing Jet-Jet MC's
//________________________________________________________________________
void AliConvEventCuts::GetXSectionAndNTrials(AliMCEvent *mcEvent, Float_t &XSection, Float_t &NTrials, AliVEvent* event){

  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;
  if (  fPeriodEnum != kLHC18b8 && fPeriodEnum != kLHC18b10 && fPeriodEnum != kLHC18l2 &&           // LHC17pq pp 5TeV JetJet MC's
        fPeriodEnum != kLHC17i3a1 &&                                                                // LHC16ijklop GammaJet MC EMCal triggered
        fPeriodEnum != kLHC17i3b1 &&                                                                // LHC16ijklop JetJet MC (with decay photon in EMC acc.)
        fPeriodEnum != kLHC17i3c1 &&                                                                // LHC16ijklop JetJet MC (with decay photon in EMC acc.)
        fPeriodEnum != kLHC17i3b2 &&                                                                // LHC16ijklop JetJet MC (with decay photon in DCal/PHOS acc.)
        fPeriodEnum != kLHC17i3c2 &&                                                                // LHC16ijklop JetJet MC (with decay photon in DCal/PHOS acc.)
        fPeriodEnum != kLHC20b1b1 &&                                                                // LHC16ijklop JetJet MC (with decay photon in EMC acc.) new prod.
        fPeriodEnum != kLHC20b1c1 &&                                                                // LHC16ijklop JetJet MC (with decay photon in EMC acc.) new prod.
        fPeriodEnum != kLHC20b1b2 &&                                                                // LHC16ijklop JetJet MC (with decay photon in DCal/PHOS acc.) new prod.
        fPeriodEnum != kLHC20b1c2 &&                                                                // LHC16ijklop JetJet MC (with decay photon in DCal/PHOS acc.) new prod.
        fPeriodEnum != kLHC17g8a &&                                                                 // LHC16qt pPb 5TeV JetJet MC's
        fPeriodEnum != kLHC16rP1JJ &&  fPeriodEnum != kLHC16sP1JJ &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC16rsGJ &&                                                                // LHC16sr pPb 8TeV GammaJet MC's
        fPeriodEnum != kLHC17g6b2a &&  fPeriodEnum != kLHC17g6b2b &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC17g6b3a &&  fPeriodEnum != kLHC17g6b3b &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC16P1JJ && fPeriodEnum != kLHC16P1JJLowB &&                               // LHC16X Jet Jet MC's
        fPeriodEnum != kLHC17P1JJ && fPeriodEnum != kLHC17P1JJLowB &&                               // LHC17X Jet Jet MC's
        fPeriodEnum != kLHC18P1JJ &&                                                                // LHC18X Jet Jet MC's
        fPeriodEnum != kLHC16h3 &&                                                                  // LHC15n Jet Jet MC's
        fPeriodEnum != kLHC15a3a && fPeriodEnum != kLHC15a3a_plus && fPeriodEnum != kLHC15a3b &&    // LHC13g Jet Jet MC's
        fPeriodEnum != kLHC15g1a && fPeriodEnum != kLHC15g1b &&                                     // LHC11a Jet Jet MC's
        fPeriodEnum != kLHC13b4_fix && fPeriodEnum != kLHC13b4_plus &&                              // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC19a4 &&                                                                  // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC16c3a && fPeriodEnum != kLHC16c3b && fPeriodEnum != kLHC16c3c &&         // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC17g6a1 && fPeriodEnum != kLHC17g6a2 && fPeriodEnum != kLHC17g6a3 &&      // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC12P2JJ  &&  fPeriodEnum != kLHC17g5b  &&  fPeriodEnum != kLHC17g5c  &&   // LHC12 JetJet MC
        fPeriodEnum != kLHC17g5a1  &&  fPeriodEnum != kLHC17g5a2 &&                                 // LHC12 GammaJet MC
        fPeriodEnum != kLHC18b11a  &&                                                               // LHC18 GammaJet MC anchored to LHC15o
        fPeriodEnum != kLHC18b11b  &&                                                               // LHC18 GammaJet MC anchored to LHC15o
        fPeriodEnum != kLHC18b11c  &&                                                               // LHC18 GammaJet MC anchored to LHC15o
        fPeriodEnum != kLHC18l6b1  &&                                                               // LHC17 GammaJet MC anchored to LHC17 (with decay photon in EMC acc.)
        fPeriodEnum != kLHC18l6c1  &&                                                               // LHC17 GammaJet MC anchored to LHC17 (with decay photon in EMC acc.)
        fPeriodEnum != kLHC18l6b2  &&                                                               // LHC17 GammaJet MC anchored to LHC17 (with decay photon in DCal/PHOS acc.)
        fPeriodEnum != kLHC18l6c2  &&                                                               // LHC17 GammaJet MC anchored to LHC17 (with decay photon in DCal/PHOS acc.)
        fPeriodEnum != kLHC19i3b1  &&                                                               // LHC18 JetJet MC anchored to LHC18 (with decay photon in EMC acc.)
        fPeriodEnum != kLHC19i3c1  &&                                                               // LHC18 JetJet MC anchored to LHC18 (with decay photon in EMC acc.)
        fPeriodEnum != kLHC19i3b2  &&                                                               // LHC18 JetJet MC anchored to LHC18 (with decay photon in DCal/PHOS acc.)
        fPeriodEnum != kLHC19i3c2  &&                                                               // LHC18 JetJet MC anchored to LHC18 (with decay photon in DCal/PHOS acc.)
        fPeriodEnum != kLHC14k1a  &&  fPeriodEnum != kLHC14k1b                                      // LHC11 JetJet MC
     ){
    NTrials = -1;
    XSection = -1;
    return;
  }

  if(mcEvent){
    cHeader                   = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
    if(cHeader) headerFound   = kTRUE;
  }else{
    //no mcEvent available -> not running on MC
    NTrials = -1;
    XSection = -1;
    return;
  }

  if(headerFound){
    TList *genHeaders         = 0x0;
    if(cHeader) genHeaders    = cHeader->GetHeaders();
    AliGenEventHeader* gh     = 0;
    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName   = gh->GetName();
      if (GeneratorName.CompareTo("AliGenPythiaEventHeader") == 0 || GeneratorName.Contains("Pythia8Jets")){
        AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(gh);
        NTrials = gPythia->Trials();
        XSection = gPythia->GetXsection();
        return;
      }
    }
  } else {
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    if(eventHeader){
      TString eventHeaderName     = eventHeader->ClassName();
      if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0 || eventHeaderName.Contains("Pythia8Jets")){
        AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader);
        NTrials = gPythia->Trials();
        XSection = gPythia->GetXsection();
        return;
      }
    }
  }
  // the following part is necessary for pythia8jets embedded in EPOS for AODs
  if(event){
    if(event->IsA()==AliAODEvent::Class()){
      AliAODMCHeader *mch = NULL;
      AliAODEvent * aod = dynamic_cast<AliAODEvent*> (event);
      if(aod){
        mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
        if ( mch ){
          Int_t nGenerators = mch->GetNCocktailHeaders();
          if ( nGenerators > 0  ){
            for(Int_t igen = 0; igen < nGenerators; igen++){
              AliGenEventHeader * eventHeaderGen = mch->GetCocktailHeader(igen) ;
              TString name = eventHeaderGen->GetName();
              if (name.CompareTo("AliGenPythiaEventHeader") == 0 || name.Contains("Pythia8Jets") || name.Contains("Pythia8GammaJet")){
                AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(eventHeaderGen);
                NTrials = gPythia->Trials();
                XSection = gPythia->GetXsection();
                return;
              }
            }
          }
        }
      }
    }
  }

  NTrials = -1;
  XSection = -1;
  return;
}


//________________________________________________________________________
// Analysing Jet-Jet MC's
//________________________________________________________________________
Float_t AliConvEventCuts::GetPtHard(AliMCEvent *mcEvent, AliVEvent* event){
  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;

  if (  fPeriodEnum != kLHC18b8 && fPeriodEnum != kLHC18b10 && fPeriodEnum != kLHC18l2 &&           // LHC17pq pp 5TeV JetJet MC's
        fPeriodEnum != kLHC19i3b1 && fPeriodEnum != kLHC19i3c1 &&                                   // LHC18 JetJet MC with decay photons in EMCal acc.
        fPeriodEnum != kLHC19i3b2 && fPeriodEnum != kLHC19i3c2 &&                                   // LHC18 JetJet MC with decay photons in DCal/PHOS acc.
        fPeriodEnum != kLHC18l6b1 && fPeriodEnum != kLHC18l6c1 &&                                   // LHC17 JetJet MC with decay photons in EMCal acc.
        fPeriodEnum != kLHC18l6b2 && fPeriodEnum != kLHC18l6c2 &&                                   // LHC17 JetJet MC with decay photons in DCal/PHOS acc.
        fPeriodEnum != kLHC17i3b1 && fPeriodEnum != kLHC17i3c1 &&                                   // LHC16 JetJet MC with decay photons in EMCal acc.
        fPeriodEnum != kLHC17i3b2 && fPeriodEnum != kLHC17i3c2 &&                                   // LHC16 JetJet MC with decay photons in DCal/PHOS acc.
        fPeriodEnum != kLHC20b1b1 && fPeriodEnum != kLHC20b1c1 &&                                   // LHC16 JetJet MC with decay photons in EMCal acc. new prod.
        fPeriodEnum != kLHC20b1b2 && fPeriodEnum != kLHC20b1c2 &&                                   // LHC16 JetJet MC with decay photons in DCal/PHOS acc. new prod.
        fPeriodEnum != kLHC17g8a &&                                                                 // LHC16qt pPb 5TeV JetJet MC's
        fPeriodEnum != kLHC16rP1JJ &&  fPeriodEnum != kLHC16sP1JJ &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC16rsGJ &&                                                                // LHC16sr pPb 8TeV GammaJet MC's
        fPeriodEnum != kLHC17g6b2a &&  fPeriodEnum != kLHC17g6b2b &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC17g6b3a &&  fPeriodEnum != kLHC17g6b3b &&                                // LHC16sr pPb 8TeV JetJet MC's
        fPeriodEnum != kLHC16P1JJ && fPeriodEnum != kLHC16P1JJLowB &&                               // LHC16X Jet Jet MC's
        fPeriodEnum != kLHC17P1JJ && fPeriodEnum != kLHC17P1JJLowB &&                               // LHC17X Jet Jet MC's
        fPeriodEnum != kLHC18P1JJ &&                                                                // LHC18X Jet Jet MC's
        fPeriodEnum != kLHC16h3 &&                                                                  // LHC15n Jet Jet MC's
        fPeriodEnum != kLHC15a3a && fPeriodEnum != kLHC15a3a_plus && fPeriodEnum != kLHC15a3b &&    // LHC13g Jet Jet MC's
        fPeriodEnum != kLHC15g1a && fPeriodEnum != kLHC15g1b &&                                     // LHC11a Jet Jet MC's
        fPeriodEnum != kLHC13b4_fix && fPeriodEnum != kLHC13b4_plus &&                              // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC19a4 &&                                                                  // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC16c3a && fPeriodEnum != kLHC16c3b && fPeriodEnum != kLHC16c3c &&         // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC17g6a1 && fPeriodEnum != kLHC17g6a2 && fPeriodEnum != kLHC17g6a3 &&      // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC12P2JJ  && fPeriodEnum != kLHC17g5b  && fPeriodEnum != kLHC17g5c  &&     // LHC12 JetJet MC
        fPeriodEnum != kLHC17g5a1  && fPeriodEnum != kLHC17g5a2  &&     // LHC12 GammaJet MC
        fPeriodEnum != kLHC14k1a  &&  fPeriodEnum != kLHC14k1b                                      // LHC11 JetJet MC
    ) return -1;

  if(mcEvent){
    cHeader           = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
    if(cHeader) headerFound   = kTRUE;
  }else{
    //no mcEvent available -> not running on MC
    return -1;
  }

  if(headerFound){
    TList *genHeaders         = 0x0;
    if(cHeader) genHeaders    = cHeader->GetHeaders();
    AliGenEventHeader* gh     = 0;
    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh             = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName   = gh->GetName();
      if (GeneratorName.CompareTo("AliGenPythiaEventHeader") == 0 || GeneratorName.Contains("Pythia8Jets")){
        return dynamic_cast<AliGenPythiaEventHeader*>(gh)->GetPtHard();
      }
    }
  } else {
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    if(eventHeader){
      TString eventHeaderName     = eventHeader->ClassName();
      if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0 || eventHeaderName.Contains("Pythia8Jets")){
        return dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->GetPtHard();
      }
    }

    if(event->IsA()==AliAODEvent::Class()){
      AliAODMCHeader *mch = NULL;
      AliAODEvent * aod = dynamic_cast<AliAODEvent*> (event);
      if(aod){
        mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
        if ( mch ){
          Int_t nGenerators = mch->GetNCocktailHeaders();
          if ( nGenerators > 0  ){
            for(Int_t igen = 0; igen < nGenerators; igen++)
            {
              AliGenEventHeader * eventHeaderGen = mch->GetCocktailHeader(igen) ;
              TString name = eventHeaderGen->GetName();
              if (name.CompareTo("AliGenPythiaEventHeader") == 0 || name.Contains("Pythia8Jets")){
                return dynamic_cast<AliGenPythiaEventHeader*>(eventHeaderGen)->GetPtHard();
              }
            }
          }
        }
      }
    }
  }

  return -1;
}


Int_t AliConvEventCuts::GetPtHardBinFromPath(const char* currFile,AliVEvent *event)
{
  TString file(currFile);
  Int_t pthard = -1;

  // Determine archive type
  TString archivetype;
  std::unique_ptr<TObjArray> walk(file.Tokenize("/"));
  for(auto const &t : *walk){
    TString &tok = static_cast<TObjString *>(t)->String();
    if(tok.Contains(".zip")){
      archivetype = tok;
      Int_t pos = archivetype.Index(".zip");
      archivetype.Replace(pos, archivetype.Length() - pos, "");
    }
  }
  if(archivetype.Length()){
    AliDebugStream(1) << "Auto-detected archive type " << archivetype << std::endl;
    Ssiz_t pos1 = file.Index(archivetype,archivetype.Length(),0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebugStream(1) << "File name: " << file << std::endl;

  // Build virtual file name
  // Support for train tests
  TString virtualFileName;
  if(file.Contains("__alice")){
    TString tmp(file);
    Int_t pos = tmp.Index("__alice");
    tmp.Replace(0, pos, "");
    tmp.ReplaceAll("__", "/");
    // cut out tag for archive and root file
    // this needs a determin
    std::unique_ptr<TObjArray> toks(tmp.Tokenize("/"));
    TString tag = "_" + archivetype;
    for(auto const &t : *toks){
      TString &path = static_cast<TObjString *>(t)->String();
      if(path.Contains(tag)){
        Int_t posTag = path.Index(tag);
        path.Replace(posTag, path.Length() - posTag, "");
      }
      virtualFileName += "/" + path;
    }
  } else {
    virtualFileName = file;
  }

  AliDebugStream(1) << "Physical file name " << file << ", virtual file name " << virtualFileName << std::endl;

  // Get the pt hard bin
  TString strPthard(virtualFileName);

  // Idea: match different informations
  // + Year clearly 2000+
  // + Run number can be match to the one in the event
  // + If we know it is not year or run number, it must be the pt-hard bin if we start from the beginning
  // The procedure is only valid for the current implementations and unable to detect non-pt-hard bins
  // It will also fail in case of arbitrary file names

  bool binfound = false;
  std::unique_ptr<TObjArray> tokens(strPthard.Tokenize("/"));
  for(auto const &t : *tokens) {
    TString &tok = static_cast<TObjString *>(t)->String();
    if(tok.IsDec()){
      Int_t number = tok.Atoi();
      if(number > 2000 && number < 3000){
        // Year
        continue;
      } else if(number == event->GetRunNumber()){
        // Run number
        continue;
      } else {
        if(!binfound){
          // the first number that is not one of the two must be the pt-hard bin
          binfound = true;
          pthard = number;
          break;
        }
      }
    }
  }
  if(!binfound) {
    AliErrorStream() << "Could not extract file number from path " << strPthard << std::endl;
  } else {
    AliDebugStream(1) << "Auto-detecting pt-hard bin " << pthard << std::endl;
  }

  AliDebugStream(1) << "File: " << file << std::endl;

  return pthard;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::MimicTrigger(AliVEvent *event, Bool_t isMC ){

  // abort if mimicing not enabled
  if (!fMimicTrigger) return kTRUE;
  if(!(fSpecialTrigger == 5 || fSpecialTrigger == 6 || fSpecialTrigger == 8 || fSpecialTrigger == 10 || fSpecialTrigger == 13  || fSpecialTrigger == 14 )) return kTRUE;   // not the correct trigger for mimcking

  // Trigger mimicking based on decision by the AliAnalysisTaskEmcalTriggerSelection for L1 triggers
  // To get the correct values one has to select the correct dataset
  if(fMimicTrigger == 2){
    auto triggercont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(event->FindListObject("EmcalTriggerDecision"));
    if(!triggercont){
      AliFatal("Trigger decision container not found in event - not possible to select EMCAL triggers");
    } else {
      if( (fSpecialTrigger == 8 || fSpecialTrigger == 10 ) && (fSpecialSubTriggerName.CompareTo("7EG2")==0 ||fSpecialSubTriggerName.CompareTo("8EG2")==0) ){
        if( (triggercont->IsEventSelected("EG2")) || (triggercont->IsEventSelected("DG2")) ) return kTRUE;
      } else if( (fSpecialTrigger == 8 || fSpecialTrigger == 10 ) && (fSpecialSubTriggerName.CompareTo("7EGA")==0 || fSpecialSubTriggerName.CompareTo("8EGA")==0 || fSpecialSubTriggerName.CompareTo("7EG1")==0 ||fSpecialSubTriggerName.CompareTo("8EG1")==0 ) ){
        if( (triggercont->IsEventSelected("EG1")) || (triggercont->IsEventSelected("DG1")) ) return kTRUE;
      } else if( fSpecialTrigger == 5 || fSpecialTrigger == 10 ){
        if( triggercont->IsEventSelected("EMCL0") || triggercont->IsEventSelected("DMCL0") ) return kTRUE;
      } else {
        return kTRUE; // In case no suitable fSpecialTrigger was selected
      }
    }
    return kFALSE;
  }
  if ((fMimicTrigger == 3)||(fMimicTrigger == 4)){
    if (fSpecialTrigger == 6){
      AliCaloTriggerMimicHelper* tempMimickHelper = 0x0;
      tempMimickHelper = (AliCaloTriggerMimicHelper*) (AliAnalysisManager::GetAnalysisManager()->GetTask(CaloTriggerHelperName.Data()));
      if (tempMimickHelper){
        return tempMimickHelper->GetEventChosenByTrigger();
      } else {
        AliFatal(Form("AliCaloTriggerMimicHelper tempMimickHelper was not found for fSpecialTrigger == %d and fMimicTrigger == %d (name:%s)", fSpecialTrigger, fMimicTrigger,CaloTriggerHelperName.Data()));
      }
    } else {
        AliFatal(Form("fSpecialTrigger == %d was not implemented in MimicTrigger case fMimicTrigger == %d", fSpecialTrigger, fMimicTrigger));
    }
  }


    //Get the clusters
    TClonesArray * arrClustersMimic = NULL;
    Int_t nclus = 0;
    if(!fCorrTaskSetting.CompareTo("")){
      nclus = event->GetNumberOfCaloClusters();
    } else {
      arrClustersMimic = dynamic_cast<TClonesArray*>(event->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
      if(!arrClustersMimic)
        AliFatal(Form("%sClustersBranch was not found in AliConvEventCuts! Check the correction framework settings!",fCorrTaskSetting.Data()));
      nclus = arrClustersMimic->GetEntries();
    }
    // return if no Clusters in the event
    if(nclus == 0)  return kFALSE;

    // Trigger mimicking based on cluster energy
    // thresholds are loaded from the OADB (OADB/PWGGA/EMCalTriggerMimicOADB.root)
    // Case1: if a cluster has energy above threshold -> accept event
    if(fMimicTrigger == 1){

      // Loading trigger thresholds from OADB
      // Only load histo if is not loaded already or if the runnumber has changed!
      Int_t runnumber = event->GetRunNumber();
      if(!fHistoTriggThresh || fRunNumberTriggerOADB != runnumber){
        fRunNumberTriggerOADB = runnumber;
        std::unique_ptr<AliOADBContainer> contfileTriggThresh(new AliOADBContainer(""));

        if(!fPathTriggerMimicSpecialInput.CompareTo("")){ // load from standard OADB on EOS
          TFile *fileTriggThresh=TFile::Open(AliDataFile::GetFileNameOADB("PWGGA/EMCalTriggerMimicOADB.root").data(),"read");
          if (!fileTriggThresh || fileTriggThresh->IsZombie())
          {
            AliFatal("OADB/PWGGA/EMCalTriggerMimicOADB.root was not found");
          }
          if (fileTriggThresh) delete fileTriggThresh;
          contfileTriggThresh->InitFromFile(AliDataFile::GetFileNameOADB("PWGGA/EMCalTriggerMimicOADB.root").data(),"AliEMCalTriggerMimic");
          if(!contfileTriggThresh){
            AliFatal("AliOADBContainer could not be loaded from PWGGA/EMCalTriggerMimicOADB.root");
          } else{
            contfileTriggThresh->SetOwner(kTRUE);
          }
        } else { // load from special OADB file from AliEn
          TFile *fileTriggThresh=TFile::Open(AliDataFile::GetFileNameOADB(((char*)Form("PWGGA/%s",fPathTriggerMimicSpecialInput.Data()))).data(),"read");
          if (!fileTriggThresh || fileTriggThresh->IsZombie())
          {
            AliFatal(Form("%s was not found",fPathTriggerMimicSpecialInput.Data()));
          }
          if (fileTriggThresh) delete fileTriggThresh;
          contfileTriggThresh->InitFromFile(AliDataFile::GetFileNameOADB(((char*)Form("PWGGA/%s",fPathTriggerMimicSpecialInput.Data()))).data(),"AliEMCalTriggerMimic");
          if(!contfileTriggThresh){
            AliFatal(Form("AliOADBContainer could not be loaded from %s",fPathTriggerMimicSpecialInput.Data()));
          } else{
            contfileTriggThresh->SetOwner(kTRUE);
          }
        }

        TObjArray *arrayTriggThresh=(TObjArray*)contfileTriggThresh->GetObject(runnumber);
        if (!arrayTriggThresh)
        {
          AliFatal(Form("No Trigger threshold found for run number: %d", runnumber));
        }
        // EMCal L0 trigger
        if( (fSpecialTrigger == 8 || fSpecialTrigger == 10 ) && (fSpecialSubTriggerName.CompareTo("7EGA")==0 || fSpecialSubTriggerName.CompareTo("8EGA")==0 ||
            fSpecialSubTriggerName.CompareTo("7EG1")==0 ||fSpecialSubTriggerName.CompareTo("8EG1")==0 ) ) fHistoTriggThresh  = (TH1S*)arrayTriggThresh->FindObject("EMCalL1G1");
        // EMCal L1 G1 trigger
        else if((fSpecialTrigger == 8 || fSpecialTrigger == 10 ) && (fSpecialSubTriggerName.CompareTo("7EG2")==0 ||fSpecialSubTriggerName.CompareTo("8EG2")==0) ) fHistoTriggThresh  = (TH1S*)arrayTriggThresh->FindObject("EMCalL1G2");
        // PHOS L0 trigger
        else if(fSpecialTrigger == 5 || fSpecialTrigger == 10) fHistoTriggThresh  = (TH1S*)arrayTriggThresh->FindObject("EMCalL0");
        // EMCal L1 G2 trigger
        else if((fSpecialTrigger == 6) && (fSpecialSubTriggerName.CompareTo("CPHI7")==0 )) fHistoTriggThresh  = (TH1S*)arrayTriggThresh->FindObject("PHOSL0");
        // return true if mimicking for fSpecialTrigger is not defined
        else return kTRUE;

        if(!fHistoTriggThresh){
          AliFatal(Form("No histogram for trigger threshold found for run number: %d", runnumber));
        }
      }


    // Get individual threshold for every Supermodule (if no Supermodulewise was defined threshold is the same for all SMs)
    Float_t fTriggThresh[20] = {0};
    fRandom.SetSeed(0);

    // load EMCal geometry if needed
    if(fSpecialTrigger != 6){
      if(!fGeomEMCAL) fGeomEMCAL = AliEMCALGeometry::GetInstance();
      if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
    }


  // Loop over EMCal clusters
    for(Int_t i = 0; i < nclus; i++){
      AliVCluster* clus = NULL;
      std::unique_ptr<AliVCluster> tmpcluster;  // takes care about deleting clusters constructed with new
      if(event->IsA()==AliESDEvent::Class()){
        if(arrClustersMimic){
          tmpcluster = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMimic->At(i)));
          clus = tmpcluster.get();
        } else
          clus = event->GetCaloCluster(i);
      } else if(event->IsA()==AliAODEvent::Class()){
        if(arrClustersMimic) {
          tmpcluster = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMimic->At(i)));
          clus = tmpcluster.get();
        }
        else
          clus = event->GetCaloCluster(i);
      }
      if (!clus) {
        continue;
      }
       if ((fSpecialTrigger!=6 && !clus->IsEMCAL()) || (fSpecialTrigger==6 && !clus->IsPHOS())) {
        continue;
      }
      if (clus->GetM02()<0.1) {
        continue;
      }
      if (clus->GetNCells()<2) {
        continue;
      }
      Int_t iSuperModule = 0;
      if (clus->IsEMCAL()){
        // Get the supermodule from cluster position
        if(fHistoTriggThresh->GetNbinsX() > 1){
          Float_t clusPos[3]={0,0,0};
          clus->GetPosition(clusPos);
          TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
          fGeomEMCAL->SuperModuleNumberFromEtaPhi(clusterVector.Eta(),clusterVector.Phi(),iSuperModule);
          if(iSuperModule >= fHistoTriggThresh->GetNbinsX() ){
            AliFatal("Supermodule nr. does not match with input histogramm");
          }
        }
      }
      if(fTriggThresh[iSuperModule] < 0.1) fTriggThresh[iSuperModule] = fRandom.Gaus(fHistoTriggThresh->GetBinContent(iSuperModule + 1)*0.01, fHistoTriggThresh->GetBinError(iSuperModule + 1)*0.01);
      if (clus->E() > fTriggThresh[iSuperModule]){
        return kTRUE;
      }
    }
  }

  // trigger mimicking applied on data clusters (software trigger)
  if(fMimicTrigger == 5){
    // find out the correct threshold for the trigger
    Double_t minClusE = 0.;
    if( (fSpecialTrigger == 13) && (fSpecialSubTriggerName.CompareTo("CEMC7_sw")==0 ))
      minClusE = 2.5;
    else if( (fSpecialTrigger == 13) && (fSpecialSubTriggerName.CompareTo("7EG2_sw")==0 ))
      minClusE = 3.5;
    else if( (fSpecialTrigger == 13) && (fSpecialSubTriggerName.CompareTo("7EG1_sw")==0 ))
      minClusE = 8.8;
    else if( (fSpecialTrigger == 14) && (fSpecialSubTriggerName.CompareTo("7EG2_EGA_sw")==0 ))
      minClusE = 4.0;
    else if( (fSpecialTrigger == 14) && (fSpecialSubTriggerName.CompareTo("7EG1_EGA_sw")==0 ))
      minClusE = 9.5;

  // Loop over clusters
    for(Int_t i = 0; i < nclus; i++){
      AliVCluster* clus = NULL;
      std::unique_ptr<AliVCluster> tmpcluster;  // takes care about deleting clusters constructed with new
      if(event->IsA()==AliESDEvent::Class()){
        if(arrClustersMimic){
          tmpcluster = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMimic->At(i)));
          clus = tmpcluster.get();
        } else
          clus = event->GetCaloCluster(i);
      } else if(event->IsA()==AliAODEvent::Class()){
        if(arrClustersMimic) {
          tmpcluster = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMimic->At(i)));
          clus = tmpcluster.get();
        }
        else
          clus = event->GetCaloCluster(i);
      }
      if (!clus) {
        continue;
      }
      if ( !clus->IsEMCAL()) {
        continue;
      }
      if (clus->GetM02()<0.1) {
        continue;
      }
      if (clus->GetNCells()<2) {
        continue;
      }
      if (clus->GetIsExotic()) {
        continue;
      }
      if (clus->IsEMCAL()){
        if(clus->E() > minClusE){
          return kTRUE;
        }
      }
    }
  }

  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::IsTriggerSelected(AliVEvent *event, Bool_t isMC)
{

  AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());


  UInt_t isSelected = AliVEvent::kAny;

  if (fInputHandler==NULL) return kFALSE;
  if( fInputHandler->GetEventSelection() || event->IsA()==AliAODEvent::Class()) {

    TString firedTrigClass = event->GetFiredTriggerClasses();
    // if no trigger has been selected manually, select kAny in case of presel (also important for AOD filtering!)
    // in other cases select standards depending on system
    if (!fTriggerSelectedManually){
      if (fPreSelCut) fOfflineTriggerMask = AliVEvent::kAny;
      else {
        if (fIsHeavyIon == 1){
          if( fPeriodEnum == kLHC18qr ){
            fOfflineTriggerMask = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
          } else {
            fOfflineTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
          }
        } else if (fIsHeavyIon == 2){
            fOfflineTriggerMask = AliVEvent::kINT7;
        } else {
            fOfflineTriggerMask = AliVEvent::kMB;
        }
      }
    }

    // in case of MC switch to kAny if no MB/INT7/INT8 has been selected
    if(isMC){
      if( fIsHeavyIon == 0){
        if( fOfflineTriggerMask != AliVEvent::kMB && fOfflineTriggerMask != AliVEvent::kINT7 && fOfflineTriggerMask != AliVEvent::kINT8 ){
          fOfflineTriggerMask = AliVEvent::kAny;
        }
      }else{
        fOfflineTriggerMask = AliVEvent::kAny;
      }
    }

    // DG event selection; special condition
    if ( (fSpecialTrigger == 11)  && fTriggerSelectedManually &&  fSpecialSubTriggerName.CompareTo("CCUP25-B-SPD1-CENTNOTRD") == 0 ) {
      if (firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 1;
    }

    // UPC event selection; special condition
    if ( (fSpecialTrigger == 12)  && fTriggerSelectedManually &&  ( fSpecialSubTriggerName.CompareTo("CCUP8") == 0 || fSpecialSubTriggerName.CompareTo("CCUP9") == 0 ) ) {
      if (firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 1;
    }

    if (fOfflineTriggerMask){
      isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
      if (isSelected && !fPreSelCut){
        // cout << firedTrigClass.Data() << endl;
        // cout << "Special trigger: "<< fSpecialTrigger << " initialized " << fEMCALTrigInitialized << endl;
        // if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9){ // EMCAL triggers
        //   if (!fEMCALTrigInitialized ) InitializeEMCALTrigger(event);
        //   fTriggersEMCAL= GetTriggerList();
        // }
        if (fSpecialSubTrigger>0 && !isMC){
          if(fNSpecialSubTriggerOptions==2){ // in case two special triggers are available
            if (fSpecialTrigger == 13 || fSpecialTrigger == 14) {
              if(!firedTrigClass.Contains(((TString)fSpecialSubTriggerName(0,4)).Data()) && !firedTrigClass.Contains(((TString)fSpecialSubTriggerNameAdditional(0,4)).Data())){
                isSelected = 0;
              }
            } else {
              if(!firedTrigClass.Contains(fSpecialSubTriggerName.Data()) && !firedTrigClass.Contains(fSpecialSubTriggerNameAdditional.Data())) isSelected = 0;
            }
          } else { // standard case for just one trigger
            if (!firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 0;
          }
          if (fRejectTriggerOverlap){
            // trigger rejection EMC1,7,8
            if (fSpecialTrigger == 5){
              if(fNSpecialSubTriggerOptions==2){
                // trigger rejection for EMC and DMC triggers together
                if (fSpecialSubTriggerName.CompareTo("CEMC7") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("CDMC7") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CEMC8") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("CDMC8") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                }
              } else {
                // separate rejection for EMC and DMC triggers
                if( fSpecialSubTriggerName.CompareTo("CEMC7") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CEMC1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kMB) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CEMC8") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CDMC7") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CDMC1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kMB) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CDMC8") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                }
              }
            }
            // gamma triggers -> no overlap with L0 and MB trigger required
            if (fSpecialTrigger == 6){
                 if( fSpecialSubTriggerName.CompareTo("CPHI7") == 0){
                     if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                 } else if( fSpecialSubTriggerName.CompareTo("CPHI8") == 0){
                     if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                 }
            }
            if (fSpecialTrigger == 8){
              // trigger rejection EGA
              if( fSpecialSubTriggerName.CompareTo("7EGA") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              } else if (fSpecialSubTriggerName.CompareTo("8EGA") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              } else if (fSpecialSubTriggerName.CompareTo("7DGA") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              } else if (fSpecialSubTriggerName.CompareTo("8DGA") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              }
              // trigger rejection EG1 & EG2
              // EG1 is the trigger with the highest threshold
              if(fNSpecialSubTriggerOptions==2){
                // trigger rejection for EMC and DMC triggers together
                if ((fSpecialSubTriggerName.CompareTo("7EG1") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG1") == 0)
                    || (fSpecialSubTriggerName.CompareTo("7EG1_EGA_sw") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG1_EGA_sw") == 0)
                  ){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                } else if ((fSpecialSubTriggerName.CompareTo("8EG1") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("8DG1") == 0)
                    || (fSpecialSubTriggerName.CompareTo("7EG1_EGA_sw") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG1_EGA_sw") == 0)
                  ){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("8DG2"))  isSelected = 0;
                } else if ((fSpecialSubTriggerName.CompareTo("7EG2") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG2") == 0)
                    || (fSpecialSubTriggerName.CompareTo("7EG2_EGA_sw") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG2_EGA_sw") == 0)
                  ){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                } else if ((fSpecialSubTriggerName.CompareTo("8EG2") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("8DG2") == 0)
                    || (fSpecialSubTriggerName.CompareTo("7EG2_EGA_sw") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG2_EGA_sw") == 0)
                  ){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                }
              } else {
                // separate rejection for EMC and DMC triggers
                if (fSpecialSubTriggerName.CompareTo("7EG1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EG1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7EG2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EG2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DG1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DG1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8DG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DG2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DG2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                }
              }
            }
            // jet triggers -> no overlap with gamma trigger and lower triggers required
            if (fSpecialTrigger == 9){
              if(fNSpecialSubTriggerOptions==2){
                // trigger rejection for EMC and DMC triggers together
                if (fSpecialSubTriggerName.CompareTo("7EJ1") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DJ1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EJ2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DJ2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7EJ2") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DJ2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG1"))  isSelected = 0;
                }
              } else {
                // separate rejection for EMC and DMC triggers
                if( fSpecialSubTriggerName.CompareTo("7EJE") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EGA"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EJE") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8EGA"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7EJ1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EJ2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EJ1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("8EG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("8EJ2"))  isSelected = 0;
                } else   if (fSpecialSubTriggerName.CompareTo("7EJ2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
                } else   if (fSpecialSubTriggerName.CompareTo("8EJ2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("8EG1"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DJ1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7DG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DJ2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DJ1") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8DG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("8DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("8DJ2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DJ2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG1"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DG2") == 0){
                  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                  if (firedTrigClass.Contains("8DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("8DG1"))  isSelected = 0;
                }
              }
            }
            if (fSpecialTrigger == 10 && (fInputHandler->IsEventSelected() & AliVEvent::kCaloOnly) ){
              if(fNSpecialSubTriggerOptions==2){
                // trigger rejection for EMC and DMC triggers together
                if (fSpecialSubTriggerName.CompareTo("7EJ1") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DJ1") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EJ2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DJ2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7EJ2") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DJ2") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG1"))  isSelected = 0;
                }
                // trigger rejection for EMC and DMC triggers together
                if (fSpecialSubTriggerName.CompareTo("7EG1") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG1") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EG1") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("8DG1") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC8-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC8-")) isSelected = 0;
                  if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
                  if (firedTrigClass.Contains("8DG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7EG2") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("7DG2") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EG2") == 0 && fSpecialSubTriggerNameAdditional.CompareTo("8DG2") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC8-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC8-")) isSelected = 0;
                }
              } else {
                // trigger rejection L0 triggers
                if (fSpecialSubTriggerName.CompareTo("CEMC7-") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CEMC1-") == 0){
                  if (firedTrigClass.Contains("INT1-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CEMC8-") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CDMC7-") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CDMC1-") == 0){
                  if (firedTrigClass.Contains("INT1-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("CDMC8-") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                }
                // trigger rejection EGA
                if (fSpecialSubTriggerName.CompareTo("7EGA") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EGA") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC8-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DGA") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DGA") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC8-")) isSelected = 0;
                }
                // trigger rejection L1 triggers
                if(fSpecialSubTriggerName.CompareTo("7EG1") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EG1") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC8-")) isSelected = 0;
                  if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7EG2") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8EG2") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("EMC8-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DG1") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC7-")) isSelected = 0;
                  if (firedTrigClass.Contains("7DG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DG1") == 0){
                  if (firedTrigClass.Contains("INT8-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC8-")) isSelected = 0;
                  if (firedTrigClass.Contains("8DG2"))  isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("7DG2") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC7-")) isSelected = 0;
                } else if (fSpecialSubTriggerName.CompareTo("8DG2") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                  if (firedTrigClass.Contains("DMC8-")) isSelected = 0;
                }
                // trigger rejection PHOS triggers
                if (fSpecialSubTriggerName.CompareTo("CPHI7-") == 0){
                  if (firedTrigClass.Contains("INT7-")) isSelected = 0;
                }
              }
            }
          }
          if (isSelected != 0 ){
            // cout << "I am here" << " :" << fSpecialSubTriggerName.Data() <<endl;
            if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9 ){
              if (hTriggerClassesCorrelated){
                if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClassesCorrelated->Fill(0);
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)hTriggerClassesCorrelated->Fill(1);
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)hTriggerClassesCorrelated->Fill(2);
                if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC7) && firedTrigClass.Contains("EMC"))hTriggerClassesCorrelated->Fill(3);
                if (firedTrigClass.Contains("7EJE") || firedTrigClass.Contains("8EJE")) hTriggerClassesCorrelated->Fill(4);
                if (firedTrigClass.Contains("7EJ1") || firedTrigClass.Contains("8EJ1")) hTriggerClassesCorrelated->Fill(5);
                if (firedTrigClass.Contains("7EJ2") || firedTrigClass.Contains("8EJ2")) hTriggerClassesCorrelated->Fill(6);
                if (firedTrigClass.Contains("7EGA") || firedTrigClass.Contains("8EGA")) hTriggerClassesCorrelated->Fill(7);
                if (firedTrigClass.Contains("7EG1") || firedTrigClass.Contains("8EG1")) hTriggerClassesCorrelated->Fill(8);
                if (firedTrigClass.Contains("7EG2") || firedTrigClass.Contains("8EG2")) hTriggerClassesCorrelated->Fill(9);
                if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC7) && firedTrigClass.Contains("DMC"))hTriggerClassesCorrelated->Fill(10);
                if (firedTrigClass.Contains("7DJE") || firedTrigClass.Contains("8DJE")) hTriggerClassesCorrelated->Fill(11);
                if (firedTrigClass.Contains("7DJ1") || firedTrigClass.Contains("8DJ1")) hTriggerClassesCorrelated->Fill(12);
                if (firedTrigClass.Contains("7DJ2") || firedTrigClass.Contains("8DJ2")) hTriggerClassesCorrelated->Fill(13);
                if (firedTrigClass.Contains("7DGA") || firedTrigClass.Contains("8DGA")) hTriggerClassesCorrelated->Fill(14);
                if (firedTrigClass.Contains("7DG1") || firedTrigClass.Contains("8DG1")) hTriggerClassesCorrelated->Fill(15);
                if (firedTrigClass.Contains("7DG2") || firedTrigClass.Contains("8DG2")) hTriggerClassesCorrelated->Fill(16);
              }
            } else if ( fSpecialTrigger == 10 ){
              if (hTriggerClassesCorrelated){
                if (fInputHandler->IsEventSelected() & AliVEvent::kCaloOnly){
                  hTriggerClassesCorrelated->Fill(0);
                  if (firedTrigClass.Contains("INT7-"))hTriggerClassesCorrelated->Fill(1);
                  if (firedTrigClass.Contains("EMC1-"))hTriggerClassesCorrelated->Fill(2);
                  if (firedTrigClass.Contains("EMC7-")|| firedTrigClass.Contains("EMC8-"))hTriggerClassesCorrelated->Fill(3);
                  if (firedTrigClass.Contains("7EJE") || firedTrigClass.Contains("8EJE")) hTriggerClassesCorrelated->Fill(4);
                  if (firedTrigClass.Contains("7EJ1") || firedTrigClass.Contains("8EJ1")) hTriggerClassesCorrelated->Fill(5);
                  if (firedTrigClass.Contains("7EJ2") || firedTrigClass.Contains("8EJ2")) hTriggerClassesCorrelated->Fill(6);
                  if (firedTrigClass.Contains("7EGA") || firedTrigClass.Contains("8EGA")) hTriggerClassesCorrelated->Fill(7);
                  if (firedTrigClass.Contains("7EG1") || firedTrigClass.Contains("8EG1")) hTriggerClassesCorrelated->Fill(8);
                  if (firedTrigClass.Contains("7EG2") || firedTrigClass.Contains("8EG2")) hTriggerClassesCorrelated->Fill(9);
                  if (firedTrigClass.Contains("DMC7-")|| firedTrigClass.Contains("DMC8-"))hTriggerClassesCorrelated->Fill(10);
                  if (firedTrigClass.Contains("7DJE") || firedTrigClass.Contains("8DJE")) hTriggerClassesCorrelated->Fill(11);
                  if (firedTrigClass.Contains("7DJ1") || firedTrigClass.Contains("8DJ1")) hTriggerClassesCorrelated->Fill(12);
                  if (firedTrigClass.Contains("7DJ2") || firedTrigClass.Contains("8DJ2")) hTriggerClassesCorrelated->Fill(13);
                  if (firedTrigClass.Contains("7DGA") || firedTrigClass.Contains("8DGA")) hTriggerClassesCorrelated->Fill(14);
                  if (firedTrigClass.Contains("7DG1") || firedTrigClass.Contains("8DG1")) hTriggerClassesCorrelated->Fill(15);
                  if (firedTrigClass.Contains("7DG2") || firedTrigClass.Contains("8DG2")) hTriggerClassesCorrelated->Fill(16);
                }
              }
            }
          }

        } else if (isMC){
          if (fSpecialTrigger == 5 || fSpecialTrigger == 6 || fSpecialTrigger == 8 || fSpecialTrigger == 9){ // EMCAL triggers
            // isSelected = 0;
            // if (fTriggersEMCAL > 0)cout << "Special Trigger " << fSpecialTrigger << " triggers: " << fTriggersEMCAL << "    selected triggers: " << fTriggersEMCALSelected << " run number: " <<event->GetRunNumber()<<endl;
            // if (fTriggersEMCAL&fTriggersEMCALSelected){
            //   cout << "accepted ++++++++++++++++++++" << endl;
              isSelected = 1;
            // }
          }
        }
        //if for specific centrality trigger selection
        if(fSpecialSubTrigger == 1){
          if(fSpecialSubTriggerName.Contains("|")  && GetCentrality(event) <= 10.){
            TObjArray *ClassesList = fSpecialSubTriggerName.Tokenize("|");
            for (Int_t i=0; i<ClassesList->GetEntriesFast();++i){
              TObjString *NameClass = (TObjString*)ClassesList->At(i);
              if (firedTrigClass.Contains(NameClass->GetString())) isSelected = 1;
            }
          } else if(fSpecialSubTriggerName.Contains("%")){
            TObjArray *ClassesList = fSpecialSubTriggerName.Tokenize("%");
            for (Int_t i=0; i<ClassesList->GetEntriesFast();++i){
              TObjString *NameClass = (TObjString*)ClassesList->At(i);
              if (firedTrigClass.Contains(NameClass->GetString())) isSelected = 1;
            }
          } else if(fSpecialSubTriggerName.Contains("@")){
            TObjArray *ClassesList = fSpecialSubTriggerName.Tokenize("@");
            for (Int_t i=0; i<ClassesList->GetEntriesFast();++i){
              TObjString *NameClass = (TObjString*)ClassesList->At(i);
              if (firedTrigClass.Contains(NameClass->GetString())) isSelected = 1;
            }
          } else if(fSpecialSubTriggerName.Contains("&")){ //logic AND of two classes
            TObjArray *ClassesList = fSpecialSubTriggerName.Tokenize("&");
            TString CheckClass = "";
            for (Int_t i=0; i<ClassesList->GetEntriesFast(); i++){
              TObjString *NameClass = (TObjString*)ClassesList->At(i);
              if (firedTrigClass.Contains(NameClass->GetString())) CheckClass+="1";
              else CheckClass+="0";
            }
            if(CheckClass.Contains("0")) isSelected = 0;
          }
          else if(firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 1;
        }
      }
    }
    //******************************************************//
    // uncomment the following lines for trigger debugging
    //   if (!fPreSelCut){
    //     if (isSelected) cout << fSpecialTrigger << "\t"<<  fSpecialTriggerName.Data() <<  "\t"<<  fSpecialSubTrigger << "\t" << fSpecialSubTriggerName.Data()  << endl;
    //   } else {
    //      cout << "Preselection: " << firedTrigClass.Data() << endl;
    //   }
    //******************************************************//
  }
  fIsSDDFired = !(fInputHandler->IsEventSelected() & AliVEvent::kFastOnly);

  Bool_t mimickedTrigger = kTRUE;
  if (fMimicTrigger) mimickedTrigger = MimicTrigger(event, isMC);
  // cout << "mimicked decision \t" << mimickedTrigger << "expect decision? "<< fMimicTrigger<< endl;
  // Fill Histogram
  if(hTriggerClass){
    if (fIsSDDFired) hTriggerClass->Fill(34);
    if (mimickedTrigger){
      if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClass->Fill(0);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)hTriggerClass->Fill(1);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUON)hTriggerClass->Fill(2);
      if (fInputHandler->IsEventSelected() & AliVEvent::kHighMult)hTriggerClass->Fill(3);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)hTriggerClass->Fill(4);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCINT5)hTriggerClass->Fill(5);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCMUS5)hTriggerClass->Fill(6);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMUSPB)hTriggerClass->Fill(6);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSH7)hTriggerClass->Fill(7);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMUSHPB)hTriggerClass->Fill(7);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUL7)hTriggerClass->Fill(8);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikePB)hTriggerClass->Fill(8);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUU7)hTriggerClass->Fill(9);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikePB)hTriggerClass->Fill(9);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7)hTriggerClass->Fill(10);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kEMC8)hTriggerClass->Fill(10);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUS7)hTriggerClass->Fill(11);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI1)hTriggerClass->Fill(12);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI7)hTriggerClass->Fill(13);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kPHI8)hTriggerClass->Fill(13);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kPHOSPb)hTriggerClass->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)hTriggerClass->Fill(14);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA)hTriggerClass->Fill(15);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)hTriggerClass->Fill(16);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)hTriggerClass->Fill(17);
      if (fInputHandler->IsEventSelected() & AliVEvent::kDG5)hTriggerClass->Fill(18);
      if (fInputHandler->IsEventSelected() & AliVEvent::kZED)hTriggerClass->Fill(19);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSPI7)hTriggerClass->Fill(20);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kSPI)hTriggerClass->Fill(20);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT8)hTriggerClass->Fill(21);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt8)hTriggerClass->Fill(22);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt8)hTriggerClass->Fill(23);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt8)hTriggerClass->Fill(24);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8)hTriggerClass->Fill(25);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt0)hTriggerClass->Fill(26);
      if (fInputHandler->IsEventSelected() & AliVEvent::kUserDefined)hTriggerClass->Fill(27);
      if (fInputHandler->IsEventSelected() & AliVEvent::kTRD)hTriggerClass->Fill(28);
      if (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly)hTriggerClass->Fill(29);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCaloOnly)hTriggerClass->Fill(30);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAnyINT)hTriggerClass->Fill(31);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAny)hTriggerClass->Fill(32);
      if (!fInputHandler->IsEventSelected()) hTriggerClass->Fill(35);
    }
    if (mimickedTrigger && fMimicTrigger) hTriggerClass->Fill(36);
  }

  if(hTriggerClassSelected && isSelected){

    if (mimickedTrigger){
      if (!fIsSDDFired) hTriggerClassSelected->Fill(33);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClassSelected->Fill(0);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)hTriggerClassSelected->Fill(1);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUON)hTriggerClassSelected->Fill(2);
      if (fInputHandler->IsEventSelected() & AliVEvent::kHighMult)hTriggerClassSelected->Fill(3);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)hTriggerClassSelected->Fill(4);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCINT5)hTriggerClassSelected->Fill(5);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCMUS5)hTriggerClassSelected->Fill(6);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMUSPB)hTriggerClassSelected->Fill(6);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSH7)hTriggerClassSelected->Fill(7);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMUSHPB)hTriggerClassSelected->Fill(7);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUL7)hTriggerClassSelected->Fill(8);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikePB)hTriggerClassSelected->Fill(8);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUU7)hTriggerClassSelected->Fill(9);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikePB)hTriggerClassSelected->Fill(9);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7)hTriggerClassSelected->Fill(10);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kEMC8)hTriggerClassSelected->Fill(10);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUS7)hTriggerClassSelected->Fill(11);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI1)hTriggerClassSelected->Fill(12);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI7)hTriggerClassSelected->Fill(13);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kPHI8)hTriggerClassSelected->Fill(13);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kPHOSPb)hTriggerClassSelected->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)hTriggerClassSelected->Fill(14);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA)hTriggerClassSelected->Fill(15);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)hTriggerClassSelected->Fill(16);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)hTriggerClassSelected->Fill(17);
      if (fInputHandler->IsEventSelected() & AliVEvent::kDG5)hTriggerClassSelected->Fill(18);
      if (fInputHandler->IsEventSelected() & AliVEvent::kZED)hTriggerClassSelected->Fill(19);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSPI7)hTriggerClassSelected->Fill(20);
    //       if (fInputHandler->IsEventSelected() & AliVEvent::kSPI)hTriggerClassSelected->Fill(20);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT8)hTriggerClassSelected->Fill(21);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt8)hTriggerClassSelected->Fill(22);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt8)hTriggerClassSelected->Fill(23);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt8)hTriggerClassSelected->Fill(24);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8)hTriggerClassSelected->Fill(25);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt0)hTriggerClassSelected->Fill(26);
      if (fInputHandler->IsEventSelected() & AliVEvent::kUserDefined)hTriggerClassSelected->Fill(27);
      if (fInputHandler->IsEventSelected() & AliVEvent::kTRD)hTriggerClassSelected->Fill(28);
      if (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly)hTriggerClassSelected->Fill(29);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAnyINT)hTriggerClassSelected->Fill(30);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAny)hTriggerClassSelected->Fill(31);
    }
    if (mimickedTrigger && fMimicTrigger) hTriggerClassSelected->Fill(34);
  }

  if(!isSelected)return kFALSE;
  if (fMimicTrigger)
    if (!mimickedTrigger ) return kFALSE;
  return kTRUE;

}

//________________________________________________________________________
TString AliConvEventCuts::GetCutNumber(){
  // returns TString with current cut number
  return fCutStringRead;
}

// todo: refactoring seems worthwhile
//________________________________________________________________________
void AliConvEventCuts::GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *event){

  if(fNotRejectedStart){
    delete[] fNotRejectedStart;
    fNotRejectedStart         = NULL;
  }
  if(fNotRejectedEnd){
    delete[] fNotRejectedEnd;
    fNotRejectedEnd         = NULL;
  }
  if(fGeneratorNames){
    delete[] fGeneratorNames;
    fGeneratorNames         = NULL;
  }

  if(rejection == 0) return; // No Rejection

  AliGenCocktailEventHeader *cHeader   = 0x0;
  AliAODMCHeader *cHeaderAOD       = 0x0;
  Bool_t headerFound           = kFALSE;
  AliMCEvent *fMCEvent           = 0x0;
  TClonesArray *fMCEventAOD       = 0x0;
  if(event->IsA()==AliMCEvent::Class()){
    if(dynamic_cast<AliMCEvent*>(event)){
      cHeader               = dynamic_cast<AliGenCocktailEventHeader*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());
      fMCEvent              = dynamic_cast<AliMCEvent*>(event);
      if(cHeader) headerFound   = kTRUE;
    }
  }
  if(event->IsA()==AliAODEvent::Class()){ // event is a AODEvent in case of AOD
    cHeaderAOD              = dynamic_cast<AliAODMCHeader*>(event->FindListObject(AliAODMCHeader::StdBranchName()));
    fMCEventAOD             = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if(cHeaderAOD) headerFound     = kTRUE;
  }

  if (fDebugLevel > 0 ) cout << "event starts here" << endl;
  if(headerFound){
    TList *genHeaders         = 0x0;
    if(cHeader) genHeaders    = cHeader->GetHeaders();
    if(cHeaderAOD){
      genHeaders              = cHeaderAOD->GetCocktailHeaders();
      if(genHeaders->GetEntries()==1){
        SetRejectExtraSignalsCut(0);
        return;
      }
    }
    AliGenEventHeader* gh     = 0;
    fnHeaders                 = 0;
    Int_t firstindexA         = 0;
    Int_t lastindexA          = -1;
    if(rejection == 1 || rejection == 3) fnHeaders = 1; // MinBiasHeader
    if(rejection == 2 || rejection == 4){ // TList of Headers Names
      for(Int_t i = 0; i<genHeaders->GetEntries();i++){
        gh                    = (AliGenEventHeader*)genHeaders->At(i);
        TString GeneratorName = gh->GetName();
        lastindexA            = lastindexA + gh->NProduced();
        if (fDebugLevel > 0 ) cout << i << "\t" << GeneratorName.Data() << endl;
        for(Int_t j = 0; j<HeaderList->GetEntries();j++){
          TString GeneratorInList   = ((TObjString*)HeaderList->At(j))->GetString();
          if (fDebugLevel > 0 )  cout << GeneratorInList.Data() << endl;
          if(GeneratorInList.Contains(GeneratorName) ){
            if (fDebugLevel > 0 ) cout << "accepted" << endl;
            if (GeneratorInList.BeginsWith("PARAM") || GeneratorInList.CompareTo("BOX") == 0 ){
              if(fMCEvent){
                if (fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c ){
                  if (fDebugLevel > 2 )cout << "number of produced particle: " <<  gh->NProduced() << endl;
                  if (fDebugLevel > 2 )cout << "pdg-code of first particle: " <<  fMCEvent->Particle(firstindexA)->GetPdgCode() << endl;
                  if (fMCEvent->Particle(firstindexA)->GetPdgCode() == fAddedSignalPDGCode ) {
                    if (gh->NProduced() > 10 && fMCEvent->Particle(firstindexA+10)->GetPdgCode() == fAddedSignalPDGCode && GeneratorInList.CompareTo("BOX") == 0){
                      if (fDebugLevel > 0 ) cout << "cond 1: "<< fnHeaders << endl;
                      fnHeaders++;
                      continue;
                    } else if (gh->NProduced() == 3 && GeneratorInList.Contains("PARAM_EMC") && (i == 3 || i == 5) ){
                      if (fDebugLevel > 2 ) cout << "accepted EMC header "<< endl;
                      if (fDebugLevel > 0 ) cout << "cond 1: "<< fnHeaders << endl;
                      fnHeaders++;
                      continue;
                    } else if (gh->NProduced() > 2 && GeneratorInList.Contains("PARAM_PHOS") && (i == 4 || i == 6) ){
                      if (fDebugLevel > 2 ) cout << "accepted PHOS header "<< endl;
                      if (fDebugLevel > 0 ) cout << "cond 1: "<< fnHeaders << endl;
                      fnHeaders++;
                      continue;

                    }
                    continue;
                  }
                } else {
                  if (fDebugLevel > 0 ) cout << "cond 2: " << fnHeaders << endl;
                  fnHeaders++;
                  continue;
                }
              }
              if ( fMCEventAOD){
                AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fMCEventAOD->At(firstindexA));
                if (aodMCParticle && (fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c) ){
                  if (  aodMCParticle->GetPdgCode() == fAddedSignalPDGCode ){
                    if (gh->NProduced() > 10 && GeneratorInList.CompareTo("BOX") == 0){
                      AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(fMCEventAOD->At(firstindexA+10));
                      if (  aodMCParticle2->GetPdgCode() == fAddedSignalPDGCode ){
                        if (fDebugLevel > 0 ) cout << "cond 1: " << fnHeaders << endl;
                        fnHeaders++;
                        continue;
                      }
                    } else if (gh->NProduced() == 3 && GeneratorInList.Contains("PARAM_EMC") && (i == 3 || i == 5) ){
                      if (fDebugLevel > 2 ) cout << "accepted EMC header "<< endl;
                      if (fDebugLevel > 0 ) cout << "cond 1: "<< fnHeaders << endl;
                      fnHeaders++;
                      continue;
                    } else if (gh->NProduced() > 2 && GeneratorInList.Contains("PARAM_PHOS") && (i == 4 || i == 6) ){
                      if (fDebugLevel > 2 ) cout << "accepted PHOS header "<< endl;
                      if (fDebugLevel > 0 ) cout << "cond 1: "<< fnHeaders << endl;
                      fnHeaders++;
                      continue;

                    }
                    continue;
                  }
                } else {
                  if (fDebugLevel > 0 ) cout << "cond 2: " << fnHeaders << endl;
                  fnHeaders++;
                  continue;
                }
              }
              continue;
            }
            if(GeneratorName.CompareTo(GeneratorInList) == 0 ){
              if (fDebugLevel > 0 ) cout << "cond 3: "<< fnHeaders << endl;
              fnHeaders++;
              continue;
            }
          }
        }
        firstindexA       = firstindexA + gh->NProduced();
      }
    }
    if (fDebugLevel > 0 ) cout << "number of headers: " <<fnHeaders << endl;

    fNotRejectedStart       = new Int_t[fnHeaders];
    fNotRejectedEnd         = new Int_t[fnHeaders];
    fGeneratorNames         = new TString[fnHeaders];

    if(rejection == 1 || rejection == 3){
      if (fPeriodEnum==kLHC20g10){
        TString lMinBiasInjectorName("Hijing_0");
        auto hasDesiredName = [&](Int_t thePos){
          return thePos<genHeaders->GetEntries() && genHeaders->At(thePos)->GetName()==lMinBiasInjectorName;
        };

        // find out position of min bias injector. Try with 7 first
        Int_t lPos = 7;
        Bool_t lFound = hasDesiredName(lPos);

        if (!lFound){
          // search for it
          if (fDebugLevel>0) cout << "Didn't find MB header at initial guess position. Looping to find it.\n";
          for (lPos=0; lPos!=genHeaders->GetEntries() && !(lFound=hasDesiredName(lPos)); ++lPos){}
        }
        if (!lFound){
          AliFatal("Could not find min bias injector header");
        }
        if (fDebugLevel>0) cout << "Found MB header at position " << lPos << endl;
        Int_t lStartIndex = 0;
        for (Int_t i=0; i<lPos; ++i){lStartIndex += ((AliGenEventHeader*)genHeaders->At(i))->NProduced();}
        fNotRejectedStart[0]    = lStartIndex;
        fNotRejectedEnd[0]      = lStartIndex + ((AliGenEventHeader*)genHeaders->At(lPos))->NProduced() - 1;
        fGeneratorNames[0]      = lMinBiasInjectorName;
      }
      else {
        // for all other productions MB header is at position 0 (thats at least how the current code reads)
        fNotRejectedStart[0]    = 0;
        fNotRejectedEnd[0]      = ((AliGenEventHeader*)genHeaders->At(0))->NProduced()-1;
        fGeneratorNames[0]      = ((AliGenEventHeader*)genHeaders->At(0))->GetName();
      }
      if (fDebugLevel > 0 ) cout << 0 << "\t" <<fGeneratorNames[0] << "\t" << fNotRejectedStart[0] << "\t" <<fNotRejectedEnd[0] << endl;
      return;
    }

    Int_t firstindex        = 0;
    Int_t lastindex         =  -1;
    Int_t number            = 0;

    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName = gh->GetName();
      lastindex             = lastindex + gh->NProduced();
      if (fDebugLevel > 0 ) cout << i << "\t" << GeneratorName.Data() << endl;
      for(Int_t j = 0; j<HeaderList->GetEntries();j++){
        TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
        if(GeneratorInList.Contains(GeneratorName) ){
          if (GeneratorInList.Contains("PARAM") || GeneratorInList.CompareTo("BOX") == 0 ){
            if(fMCEvent){
              if (fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c ){
                if (fMCEvent->Particle(firstindex)->GetPdgCode() == fAddedSignalPDGCode ) {
                  if (fDebugLevel > 0 ) cout << "produced " << gh->NProduced() << " with box generator" << endl;
                  if (gh->NProduced() > 10 && fMCEvent->Particle(firstindex+10)->GetPdgCode() == fAddedSignalPDGCode && GeneratorInList.CompareTo("BOX") == 0){
                    if (fDebugLevel > 0 ) cout << "one of them was a pi0 or eta" <<  endl;
                    fNotRejectedStart[number] = firstindex;
                    fNotRejectedEnd[number] = lastindex;
                    fGeneratorNames[number] = GeneratorName;
                    number++;
                    if (fDebugLevel > 0 ) cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
                    continue;
                  } else if (gh->NProduced() == 3 && GeneratorInList.Contains("PARAM_EMC") && (i == 3 || i == 5) ){
                    fNotRejectedStart[number] = firstindex;
                    fNotRejectedEnd[number] = lastindex;
                    fGeneratorNames[number] = GeneratorName;
                    number++;
                    if (fDebugLevel > 0 ) cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
                    continue;
                  } else if (gh->NProduced() > 2 && GeneratorInList.Contains("PARAM_PHOS") && (i == 4 || i == 6) ){
                    fNotRejectedStart[number] = firstindex;
                    fNotRejectedEnd[number] = lastindex;
                    fGeneratorNames[number] = GeneratorName;
                    number++;
                    if (fDebugLevel > 0 ) cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
                    continue;
                  }
                }
                continue;
              } else {
                fNotRejectedStart[number] = firstindex;
                fNotRejectedEnd[number] = lastindex;
                fGeneratorNames[number] = GeneratorName;
                number++;
                continue;
              }
            }
            if ( fMCEventAOD){
              AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fMCEventAOD->At(firstindex));
              if (fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c ){
                if (  aodMCParticle->GetPdgCode() == fAddedSignalPDGCode ){
                  if (gh->NProduced() > 10 && GeneratorInList.CompareTo("BOX") == 0) {
                    AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(fMCEventAOD->At(firstindex+10));
                    if ( aodMCParticle2->GetPdgCode() == fAddedSignalPDGCode ){
                      fNotRejectedEnd[number] = lastindex;
                      fNotRejectedStart[number] = firstindex;
                      fGeneratorNames[number] = GeneratorName;
                      number++;
                    }
                    continue;
                  } else if (gh->NProduced() == 3 && GeneratorInList.Contains("PARAM_EMC") && (i == 3 || i == 5) ){
                    fNotRejectedStart[number] = firstindex;
                    fNotRejectedEnd[number] = lastindex;
                    fGeneratorNames[number] = GeneratorName;
                    number++;
                    if (fDebugLevel > 0 ) cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
                    continue;
                  } else if (gh->NProduced() > 2 && GeneratorInList.Contains("PARAM_PHOS") && (i == 4 || i == 6) ){
                    fNotRejectedStart[number] = firstindex;
                    fNotRejectedEnd[number] = lastindex;
                    fGeneratorNames[number] = GeneratorName;
                    number++;
                    if (fDebugLevel > 0 ) cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
                    continue;
                  }
                  continue;
                }
              } else {
                fNotRejectedStart[number] = firstindex;
                fNotRejectedEnd[number] = lastindex;
                fGeneratorNames[number] = GeneratorName;
                number++;
                continue;
              }
            }
            continue;
          } else if(GeneratorName.CompareTo(GeneratorInList) == 0 ){
            fNotRejectedStart[number] = firstindex;
            fNotRejectedEnd[number] = lastindex;
            fGeneratorNames[number] = GeneratorName;
            if (fDebugLevel > 0 )  cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
            number++;
            continue;
          }

        }
      }
      firstindex           = firstindex + gh->NProduced();
    }
    if (fDebugLevel > 0 ) {
      for (Int_t i = 0; i < number; i++){
        cout << i << "\t" <<fGeneratorNames[i] << "\t" << fNotRejectedStart[i] << "\t" <<fNotRejectedEnd[i] << endl;
      }
    }
  } else { // No Cocktail Header Found
    fNotRejectedStart         = new Int_t[1];
    fNotRejectedEnd         = new Int_t[1];

    fnHeaders             = 1;
    fNotRejectedStart[0]       = 0;
    fNotRejectedEnd[0]         = static_cast<AliMCEvent*>(event)->GetNumberOfPrimaries()-1;
    if (rejection > 1){
      fNotRejectedStart[0]     = -1;
      fNotRejectedEnd[0]       = -1;
    }

    fGeneratorNames         = new TString[1];
    fGeneratorNames[0]         = "NoCocktailGeneratorFound";
//     SetRejectExtraSignalsCut(0);
  }

}

/* checks for a photon candidate if its tracks come from an accepted injector.
 * return value kFALSE: photon gets discarded completely
 *              kTRUE:  photon will get added to fGammaCandidates
 *
 * (only for return value kTRUE):
 * theIsFromSelectedHeader=kTRUE:  photon will get processed fully (gamma histos will get filled)
 * theIsFromSelectedHeader=kFALSE: will only get added to fGammaCandidates
 *
 * note: IsParticleFromBGEvent()>0 : particle comes from an injector on AliConvEventCuts::fHeaderList
 *                              ==2: injector is first on that list */
//________________________________________________________________________
Bool_t AliConvEventCuts::PhotonPassesAddedParticlesCriterion(AliMCEvent             *theMCEvent,
                                                             AliVEvent              *theInputEvent,
                                                             AliAODConversionPhoton &thePhoton,
                                                             Bool_t                 &theIsFromSelectedHeader)
{
  theIsFromSelectedHeader = kTRUE;
  if (!fRejectExtraSignals) return kTRUE;

  auto bothFromFirstHeader = [](Int_t thePos, Int_t theNeg){ return (thePos+theNeg)==4; };
  Int_t lIsPosFromMBHeader = IsParticleFromBGEvent(thePhoton.GetMCLabelPositive(), theMCEvent, theInputEvent);

  if (fRejectExtraSignals==3){
    Int_t lIsNegFromMBHeader = IsParticleFromBGEvent(thePhoton.GetMCLabelNegative(), theMCEvent, theInputEvent);
    theIsFromSelectedHeader = bothFromFirstHeader(lIsNegFromMBHeader, lIsPosFromMBHeader);
  }
  else{ // 1,2,4
    if (!lIsPosFromMBHeader) return kFALSE;
    Int_t lIsNegFromMBHeader = IsParticleFromBGEvent(thePhoton.GetMCLabelNegative(), theMCEvent, theInputEvent);
    if (!lIsNegFromMBHeader) return kFALSE;

    if (fRejectExtraSignals!=2) {
      theIsFromSelectedHeader = bothFromFirstHeader(lIsNegFromMBHeader, lIsPosFromMBHeader);
    }
  }
  return kTRUE;
}

//_________________________________________________________________________
Int_t AliConvEventCuts::IsParticleFromBGEvent(Int_t index, AliMCEvent *mcEvent, AliVEvent *InputEvent, Int_t debug ){

  //   if (debug > 2 ) cout << index << endl;
  if(index < 0) return 0; // No Particle

  Int_t accepted = 0;
  if(!InputEvent || InputEvent->IsA()==AliESDEvent::Class()){
    if(!mcEvent) return 0; // no mcEvent available, return 0
    if(index >= mcEvent->GetNumberOfPrimaries()){ // initial particle is secondary particle
      if( ((TParticle*)mcEvent->Particle(index))->GetMother(0) < 0) return 0; // material particle, return 0
      return IsParticleFromBGEvent(((TParticle*)mcEvent->Particle(index))->GetMother(0),mcEvent,InputEvent, debug);
    }
    for(Int_t i = 0;i<fnHeaders;i++){
      //       if (debug > 2 ) cout << "header " << fGeneratorNames[i].Data() << ":"<< fNotRejectedStart[i] << "\t" << fNotRejectedEnd[i] << endl;
      if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
        if (debug > 1 ) cout << "accepted:" << index << "\t header " << fGeneratorNames[i].Data()  << ": "<< fNotRejectedStart[i] << "\t" << fNotRejectedEnd[i] << endl;
        accepted = 1;
        if(i == 0) accepted = 2; // MB Header
      }
    }
    if (debug > 1 && !accepted) cout << "rejected:" << index << endl;
  }
  else if(InputEvent->IsA()==AliAODEvent::Class()){
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(InputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (fAODMCTrackArray){
      AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
      if(!aodMCParticle) return 0; // no particle

      if(!aodMCParticle->IsPrimary()){
        if( aodMCParticle->GetMother() < 0) return 0;// material particle, return 0
        return IsParticleFromBGEvent(aodMCParticle->GetMother(),mcEvent,InputEvent, debug);
      }
      index = TMath::Abs(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index))->GetLabel());

      for(Int_t i = 0;i<fnHeaders;i++){
        if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
          accepted = 1;
          if(i == 0) accepted = 2; // MB Header
        }
      }
    }
  }

  return accepted;
}

// returns name of header that particle belongs to
// by looping over all headers available
// DOES NOT YET WORK FOR EMBEDDED IN AOD
//_________________________________________________________________________
TString AliConvEventCuts::GetParticleHeaderName(Int_t index, AliMCEvent *mcEvent, AliVEvent *InputEvent, Int_t debug ){

  //   if (debug > 2 ) cout << index << endl;
  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound           = kFALSE;
  TString headername = "";
  if(index < 0) return headername; // No Particle

  // Get ALL headers
  if(mcEvent){
     cHeader           =  dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
     if(cHeader) headerFound = kTRUE;
  } else{
    AliFatal("No MC event found");
  }

  if(headerFound){
    TList* genHeaders = 0x0;
    if(cHeader) genHeaders    = cHeader->GetHeaders();
    AliGenEventHeader* gh     = 0;
    if(!InputEvent || InputEvent->IsA()==AliESDEvent::Class()){
      if(!mcEvent) return headername; // no mcEvent available, return 0
      if(index >= mcEvent->GetNumberOfPrimaries()){ // initial particle is secondary particle
        if( ((TParticle*)mcEvent->Particle(index))->GetMother(0) < 0) return headername; // material particle, return 0
        return GetParticleHeaderName(((TParticle*)mcEvent->Particle(index))->GetMother(0),mcEvent,InputEvent, debug);
      }

      // Loop over gen headers
      Int_t firstindex        = 0;
      Int_t lastindex         =  -1;
      for(Int_t i = 0; i<genHeaders->GetEntries();i++){
        gh = (AliGenEventHeader*)genHeaders->At(i);
        TString GeneratorName = gh->GetName();
        lastindex             = lastindex + gh->NProduced();
        if(index >= firstindex && index <= lastindex){
          // cout << "accepted:" << index << "\t header " << GeneratorName.Data() << endl;
          headername = GeneratorName;
        }
        firstindex           = firstindex + gh->NProduced();
      }
    }
    else if(InputEvent->IsA()==AliAODEvent::Class()){
      if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(InputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (fAODMCTrackArray){
        AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
        if(!aodMCParticle) return headername; // no particle
        index = TMath::Abs(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index))->GetLabel());

        // Loop over gen headers
        Int_t firstindex        = 0;
        Int_t lastindex         =  -1;
        for(Int_t i = 0; i<genHeaders->GetEntries();i++){
          gh = (AliGenEventHeader*)genHeaders->At(i);
          TString GeneratorName = gh->GetName();
          lastindex             = lastindex + gh->NProduced();
          if(index >= firstindex && index <= lastindex){
            headername = GeneratorName;
          }
          firstindex           = firstindex + gh->NProduced();
        }
      }
    }
  } else {
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    headername     = eventHeader->ClassName();
  }
  return headername;
}

//_________________________________________________________________________
Int_t AliConvEventCuts::IsEventAcceptedByCut(AliConvEventCuts *ReaderCuts, AliVEvent *event, AliMCEvent *mcEvent, Int_t isHeavyIon, Bool_t isEMCALAnalysis){

  Bool_t isMC = kFALSE;
  if (mcEvent){isMC = kTRUE;}

  if ( !IsTriggerSelected(event, isMC) )
    return 3;

  if( !(IsCentralitySelected(event,mcEvent)))
    return 1; // Check Centrality --> Not Accepted => eventQuality = 1

  Bool_t hasV0And = ReaderCuts->HasV0AND();
  Bool_t isSDDFired = ReaderCuts->IsSDDFired();

  if( ( (IsSpecialTrigger() == 0 && IsSpecialSubTrigger() == 1) || (IsSpecialTrigger() == 1 && IsSpecialSubTrigger() == 1) ) && !isSDDFired && !mcEvent)
  //if V0OR with SDD requested or V0AND with SDD request but the SDD has not fired
  return 7; // V0 with SDD requested but no fired

  if( ( (IsSpecialTrigger() == 1 && IsSpecialSubTrigger() == 0) || (IsSpecialTrigger() == 1 && IsSpecialSubTrigger() == 1) ) && !hasV0And)
  //if V0AND (only) or V0AND with SDD requested but V0AND requested but no fired
  return 8; // V0AND requested but no fired


  if( (IsSpecialTrigger() == 2 || IsSpecialTrigger() == 3) && !isSDDFired && !mcEvent)
    return 7; // With SDD requested but no fired

  if( (IsSpecialTrigger() == 1 || IsSpecialTrigger() == 3) && !hasV0And)
    return 8; // V0AND requested but no fired

  // Special EMCAL checks due to hardware issues in LHC11a or LHC12x or LHC16rs
  if (isEMCALAnalysis || IsSpecialTrigger() == 5 || IsSpecialTrigger() == 8 || IsSpecialTrigger() == 9 ){
    Int_t runnumber = event->GetRunNumber();
    if ((runnumber>=144871) && (runnumber<=146860)) {

      AliVCaloCells *cells   = event->GetEMCALCells();
      const Short_t nCells   = cells->GetNumberOfCells();

      if (event->IsA()==AliESDEvent::Class()) AliAnalysisManager::GetAnalysisManager()->LoadBranch("EMCALCells.");

      AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      if (!fInputHandler) return 3;

      // count cells above threshold
      Int_t nCellCount[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
      for(Int_t iCell=0; iCell<nCells; ++iCell) {
        Short_t cellId = cells->GetCellNumber(iCell);
        Double_t cellE = cells->GetCellAmplitude(cellId);
        Int_t sm       = cellId / (24*48);
        if (cellE>0.1) ++nCellCount[sm];
      }

      Bool_t fIsLedEvent = kFALSE;
      if (nCellCount[4] > 100) {
        fIsLedEvent = kTRUE;
      } else {
        if ((runnumber>=146858) && (runnumber<=146860)) {
          if ((fInputHandler->IsEventSelected() & AliVEvent::kMB) && (nCellCount[3]>=21))
            fIsLedEvent = kTRUE;
          else if ((fInputHandler->IsEventSelected() & AliVEvent::kEMC1) && (nCellCount[3]>=35))
            fIsLedEvent = kTRUE;
        }
      }
      if (fIsLedEvent) {
        return 9;
      }
    }
    Bool_t fRejectEMCalLEDevents = kTRUE;
    if (fRejectEMCalLEDevents && (fPeriodEnum == kLHC12 || fPeriodEnum == kLHC16NomB || fPeriodEnum == kLHC17NomB || fPeriodEnum == kLHC18NomB)) {
      AliVCaloCells *cells   = event->GetEMCALCells();
      const Short_t nCells   = cells->GetNumberOfCells();

      if(!fGeomEMCAL) fGeomEMCAL = AliEMCALGeometry::GetInstance();
      if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}

      // count cells above threshold
      Int_t nStripsLED = 0;
      Int_t nCellCountInStrip[480] = {0};

      Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0, row=0, column=0;
      for(Int_t iCell=0; iCell<nCells; ++iCell) {
        // Get SM number and relative row/column for SM
        fGeomEMCAL->GetCellIndex(cells->GetCellNumber(iCell), nSupMod,nModule,nIphi,nIeta);
        fGeomEMCAL->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, row,column);

        Short_t cellId = cells->GetCellNumber(iCell);
        Double_t cellE = cells->GetCellAmplitude(cellId);
        Int_t strip       = nSupMod*24 + column/2;
        if (cellE>0.1) nCellCountInStrip[strip]++;
      }

      Bool_t fIsLedEvent = kFALSE;
      for(Int_t istrip=0; istrip<480; istrip++) {
        if(nCellCountInStrip[istrip]>40) nStripsLED++;
      }
      if(nStripsLED>=3) fIsLedEvent = kTRUE;
      if (fIsLedEvent) {
        return 9;
      }
    }
  }

  // SPD clusters vs tracklets to check for pileup/background
  Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
  Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
  if(hSPDClusterTrackletBackgroundBefore) hSPDClusterTrackletBackgroundBefore->Fill(nTracklets, (nClustersLayer0 + nClustersLayer1));


  Double_t distZMax     = 0;
  if(event->IsA()==AliESDEvent::Class()){
    Int_t nPileVert = ((AliESDEvent*)event)->GetNumberOfPileupVerticesSPD();
    if (hNPileupVertices) hNPileupVertices->Fill(nPileVert);
    if (nPileVert > 0){
      for(Int_t i=0; i<nPileVert;i++){
        const AliESDVertex* pv= ((AliESDEvent*)event)->GetPileupVertexSPD(i);
        Int_t nc2             = pv->GetNContributors();
        if(nc2>=3){
          Double_t z1 = ((AliESDEvent*)event)->GetPrimaryVertexSPD()->GetZ();
          Double_t z2 = pv->GetZ();
          Double_t distZ  = z2-z1;
          if (TMath::Abs(distZMax) <  TMath::Abs(distZ) ){
            distZMax      = distZ;
          }
        }
      }
      if (hPileupVertexToPrimZ) hPileupVertexToPrimZ->Fill(distZMax);
    }
  }
  if(GetPastFutureLowBC()!=0 && GetPastFutureHighBC()!=0 ){
    if(IsOutOfBunchPileupPastFuture(event))
      return 12;
  }

  if( isHeavyIon != 2 && GetIsFromPileupSPD()){
    if(event->IsPileupFromSPD(3,0.8,3.,2.,5.) ){
      if (hPileupVertexToPrimZSPDPileup) hPileupVertexToPrimZSPDPileup->Fill(distZMax);
      return 6; // Check Pileup --> Not Accepted => eventQuality = 6
    }
    if (fUtils->IsSPDClusterVsTrackletBG(event)){
      if (hPileupVertexToPrimZTrackletvsHits) hPileupVertexToPrimZTrackletvsHits->Fill(distZMax);
      return 11; // Check Pileup --> Not Accepted => eventQuality = 11
    }
  }
  if(isHeavyIon == 2 && GetIsFromPileupSPD()){
    if(fUtils->IsPileUpEvent(event) ){
      if (hPileupVertexToPrimZSPDPileup) hPileupVertexToPrimZSPDPileup->Fill(distZMax);
      return 6; // Check Pileup --> Not Accepted => eventQuality = 6
    }
    if (fUtils->IsSPDClusterVsTrackletBG(event)){
      if (hPileupVertexToPrimZTrackletvsHits) hPileupVertexToPrimZTrackletvsHits->Fill(distZMax);
      return 11; // Check Pileup --> Not Accepted => eventQuality = 11
    }
  }

  if(fRemovePileUp){
    if (   (GetDoPileUpRejectV0MTPCout() && IsPileUpV0MTPCout(event))
        || (fRemovePileUpSDDSSDTPC && IsPileUpSDDSSDTPC(event))) {
      return 13;
     }
  }

  if(fUseSphericity > 0){
    Double_t eventSphericity  = -1;
    Int_t nPrimTracks         = GetV0Reader()->GetNumberOfPrimaryTracks();
    Double_t InAcceptance     = GetV0Reader()->IsSphericityAxisInEMCalAcceptance();
    if(fUseSphericityTrue){
      eventSphericity         = GetV0Reader()->GetSphericityTrue();
    } else if(!fUseSphericityTrue){
      eventSphericity         = GetV0Reader()->GetSphericity();
    }
    if (eventSphericity == -1) return 14;

    if (fUseSphericity == 2){
      if (nPrimTracks > 20) return 14;
    } else if (fUseSphericity == 3){
      if (nPrimTracks <= 20) return 14;
    } else if (fUseSphericity == 4){
      if (!InAcceptance) return 14;
    } else if (fUseSphericity == 5){
      if (InAcceptance) return 14;
    }

    if ((fCentralityMin < fCentralityMax) ){
      if ( eventSphericity < (Double_t)fCentralityMin/10 || eventSphericity > (Double_t)fCentralityMax/10)
        return 14;
    }
  }

  if(hCentrality)hCentrality->Fill(GetCentrality(event));

  if(hVertexZ)hVertexZ->Fill(event->GetPrimaryVertex()->GetZ());
//  if(hCentralityVsNumberOfPrimaryTracks)
//    hCentralityVsNumberOfPrimaryTracks->Fill(GetCentrality(event),
//                        ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()
//                        ->GetTask(fV0ReaderName.Data()))->GetNumberOfPrimaryTracks());

  if(fIsHeavyIon == 1){
    AliEventplane *EventPlane = event->GetEventplane();
    fEventPlaneAngle = EventPlane->GetEventplane("V0",event,2);
    if(hEventPlaneAngle)hEventPlaneAngle->Fill(TMath::Abs(fEventPlaneAngle));
  }
  if(hSPDClusterTrackletBackground) hSPDClusterTrackletBackground->Fill(nTracklets, (nClustersLayer0 + nClustersLayer1));

  return 0;
}


//_________________________________________________________________________
Float_t AliConvEventCuts::GetWeightForCentralityFlattening(AliVEvent *event){

  AliInfo("Inside the GetWeightForCentralityFlattening function");
  Double_t centrality = 0.;
  //obtain centrality for ESD or AOD
  if(!event || event->IsA()==AliESDEvent::Class()){
    AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
    if(esdEvent){
      AliCentrality *fESDCentrality=(AliCentrality*)esdEvent->GetCentrality();
      if(fDetectorCentrality==0 && fIsHeavyIon==1){
          centrality = fESDCentrality->GetCentralityPercentile("V0M"); // default for PbPb
      }
    }
  } else if(event->IsA()==AliAODEvent::Class()){
    AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);
    if(aodEvent){
      if(aodEvent->GetHeader()){
        centrality = ((AliVAODHeader*)aodEvent->GetHeader())->GetCentrality();
      }
    }
  }

  //Get the maximum vlaue from the reference distribution and interpolated value
  Float_t GetValueForWeight = 1.;
  Float_t maximum = 1.;
  Double_t weightCentrality = 1.;
  Bool_t CorrCentrLoop = kFALSE;

  //depending on the value of the flag, flattening in different cent. range
  if ( fDoCentralityFlat == 1 && (centrality >= 0. && centrality <= 20.) ){
    GetValueForWeight = hCentralityNotFlat->Interpolate(centrality);
    maximum = hCentralityNotFlat->GetMaximum();
    CorrCentrLoop = kTRUE;
  } else if ( fDoCentralityFlat == 8 ){
    GetValueForWeight = hCentralityNotFlat->Interpolate(centrality);
    maximum = hCentralityNotFlat->GetMaximum();
    CorrCentrLoop = kTRUE;
  } else {
    CorrCentrLoop = kFALSE;
  }

  if (CorrCentrLoop && GetValueForWeight != 0. && maximum !=0. && isfinite(GetValueForWeight) && isfinite(maximum) ){
      weightCentrality = maximum/GetValueForWeight;
      if (!isfinite(GetValueForWeight)) weightCentrality = 1.;
      if (!isfinite(weightCentrality)) weightCentrality = 1.;
  }

  return weightCentrality;
}

//_________________________________________________________________________
Float_t AliConvEventCuts::GetWeightForMultiplicity(Int_t mult){

  Double_t weightMult         = 1.;

  if (hReweightMultData == NULL || hReweightMultMC == NULL ) return weightMult;

  // get mult values for weights
  Float_t valueMultData       = -1.;
  Float_t valueMultMC         = -1.;
  valueMultData               = hReweightMultData->Interpolate(mult);
  valueMultMC                 = hReweightMultMC->Interpolate(mult);

  // calculate relative error for data and MC
  Float_t valueMC   = 0;
  Float_t valueData = 0;
  Float_t errorMC   = 0;
  Float_t errorData = 0;
  valueMC   = hReweightMultMC->GetBinContent(hReweightMultMC->FindBin(mult));
  valueData = hReweightMultData->GetBinContent(hReweightMultData->FindBin(mult));
  errorMC   = hReweightMultMC->GetBinError(hReweightMultMC->FindBin(mult));
  errorData = hReweightMultData->GetBinError(hReweightMultData->FindBin(mult));
  Float_t relativeErrorMC   = 1;
  Float_t relativeErrorData = 1;
  if(valueMC!=0)   relativeErrorMC   = errorMC / valueMC;
  if(valueData!=0) relativeErrorData = errorData / valueData;

  Double_t errorTolerance = 0.2;
  Double_t cutOff = -1.0;  // if > 0 : if rel error too large => weight with 0 instead of 1 above this value

  //  For these periods allow larger statistical error in the MC to apply the multiplicity weight
  if ( fPeriodEnum == kLHC16NomB || fPeriodEnum == kLHC16P1Pyt8 || fPeriodEnum == kLHC16P1PHO ||
       fPeriodEnum == kLHC17pq  ||  fPeriodEnum == kLHC17P1PHO  || fPeriodEnum == kLHC17l3b  || fPeriodEnum == kLHC18j2  || fPeriodEnum == kLHC17l4b){
    errorTolerance = 0.6;
  }

  if ( fPeriodEnum == kLHC18e1  || fPeriodEnum == kLHC16h4  ||
       fPeriodEnum == kLHC18e1a || fPeriodEnum == kLHC18e1b || fPeriodEnum == kLHC18e1c ||
       fPeriodEnum == kLHC16i1a || fPeriodEnum == kLHC16i1b || fPeriodEnum == kLHC16i1c ||
       fPeriodEnum == kLHC16i2a || fPeriodEnum == kLHC16i2b || fPeriodEnum == kLHC16i2c ||
       fPeriodEnum == kLHC16i3a || fPeriodEnum == kLHC16i3b || fPeriodEnum == kLHC16i3c ) {
    errorTolerance = 0.6;
    cutOff = 2800;  // MC distribution still significant while data error already too large
      }

  if (relativeErrorData < errorTolerance && relativeErrorMC < errorTolerance ){
    if (isfinite(valueMultData) && isfinite(valueMultMC) ){
      weightMult               = valueMultData/valueMultMC;
    }
  } else if(cutOff>0){
    if(mult>cutOff) {
      weightMult = 0;
    }
  }

  return weightMult;
}



//_________________________________________________________________________
Float_t AliConvEventCuts::GetWeightForMeson(Int_t index, AliMCEvent *mcEvent, AliVEvent *event){
  if(index < 0) return 0; // No Particle

  // check if MC production should be weighted. If it is with added particles check that particle is not rejected
  Int_t kCaseGen = 0;
  if (fPeriodEnum == kLHC13d2   || fPeriodEnum == kLHC13d2b       ||                                                                      // LHC10h MCs
      fPeriodEnum == kLHC14a1a  || fPeriodEnum == kLHC14a1b       || fPeriodEnum == kLHC14a1c   ||                                        // LHC11h MCs
      fPeriodEnum == kLHC13e7   || fPeriodEnum == kLHC13b2_efix   || fPeriodEnum == kLHC14b2    ||  fPeriodEnum == kLHC18j5  ||           // LHC13bc MCs
      fPeriodEnum == kLHC14e2b  ||                                                                                                        // LHC12[a-i] pass 1 MCs
      fPeriodEnum == kLHC12f1a  || fPeriodEnum == kLHC12f1b       || fPeriodEnum == kLHC12i3    ||                                        // LHC11a MCs
      fPeriodEnum == kLHC16h4   || fPeriodEnum == kLHC19h3        || fPeriodEnum == kLHC20g10 )                                           // LHC15o, LHC18qr pass1, LHC18qr pass3  MCs
    kCaseGen = 1;  // added particles MC
  if( fPeriodEnum == kLHC18e1 || fPeriodEnum == kLHC18e1a || fPeriodEnum == kLHC18e1b || fPeriodEnum == kLHC18e1c || fPeriodEnum == kLHC16i1a || fPeriodEnum == kLHC16i1b || fPeriodEnum == kLHC16i1c || fPeriodEnum == kLHC16i2a || fPeriodEnum == kLHC16i2b || fPeriodEnum == kLHC16i2c || fPeriodEnum == kLHC16i3a || fPeriodEnum == kLHC16i3b || fPeriodEnum == kLHC16i3c || // LHC15o MCs
      fPeriodEnum == kLHC12P2JJ || fPeriodEnum == kLHC16h3 || fPeriodEnum == kLHC18b8 || fPeriodEnum == kLHC16rP1JJ || fPeriodEnum == kLHC16sP1JJ  || fPeriodEnum == kLHC16rsGJ || fPeriodEnum == kLHC17g8a ||
      fPeriodEnum == kLHC18f3 ||  //LHC16qt MCs
      fPeriodEnum == kLHC17l3b || fPeriodEnum == kLHC18j2 || //LHC17pq MCs
      fPeriodEnum == kLHC18l8a || fPeriodEnum == kLHC18l8b ||  fPeriodEnum == kLHC18l8c || // LHC18qr MC
      fPeriodEnum == kLHC19h2a || fPeriodEnum == kLHC19h2b ||  fPeriodEnum == kLHC19h2c || // LHC18qr MC
      fPeriodEnum == kLHC20e3a || fPeriodEnum == kLHC20e3b ||  fPeriodEnum == kLHC20e3c )  // LHC18qr MC pass3
    kCaseGen = 2;  // regular MC



  if (kCaseGen == 0) return 1.;
  if(kCaseGen==1 && !IsParticleFromBGEvent(index, mcEvent, event)) return 1.;

  // get pT and pdg code
  Double_t mesonPt = 0;
  //Double_t mesonMass = 0;
  Int_t PDGCode = 0;
  if(!event || event->IsA()==AliESDEvent::Class()){
    mesonPt = ((TParticle*)mcEvent->Particle(index))->Pt();
    //mesonMass = ((TParticle*)mcEvent->Particle(index))->GetCalcMass();
    PDGCode = ((TParticle*)mcEvent->Particle(index))->GetPdgCode();
  } else if(event->IsA()==AliAODEvent::Class()){
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (fAODMCTrackArray){
      AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
      mesonPt = aodMCParticle->Pt();
      //mesonMass = aodMCParticle->GetCalcMass();
      PDGCode = aodMCParticle->GetPdgCode();
    } else {
      return 1;
    }
  }

  // get MC value
  Float_t functionResultMC = 1.;
  if ( PDGCode ==  111 && fDoReweightHistoMCPi0 && hReweightMCHistPi0!= 0x0){
    functionResultMC = hReweightMCHistPi0->Interpolate(mesonPt);
  }
  if ( PDGCode ==  221 && fDoReweightHistoMCEta && hReweightMCHistEta!= 0x0){
    functionResultMC = hReweightMCHistEta->Interpolate(mesonPt);
  }
  if ( PDGCode ==  310 && fDoReweightHistoMCK0s && hReweightMCHistK0s!= 0x0){
    functionResultMC = hReweightMCHistK0s->Interpolate(mesonPt);
  }

  // get data value
  Float_t functionResultData = 1;
  if ( PDGCode ==  111 && fDoReweightHistoMCPi0 && fFitDataPi0!= 0x0){
    functionResultData = fFitDataPi0->Eval(mesonPt);
  }
  if ( PDGCode ==  221 && fDoReweightHistoMCEta && fFitDataEta!= 0x0){
    functionResultData = fFitDataEta->Eval(mesonPt);
  }
  if ( PDGCode ==  310 && fDoReweightHistoMCK0s && fFitDataK0s!= 0x0){
    functionResultData = fFitDataK0s->Eval(mesonPt);
  }

  // calculate weight from data and MC
  Double_t weight = 1;
  if (PDGCode ==  111 || PDGCode ==  221){
    if (functionResultData != 0. && functionResultMC != 0. && isfinite(functionResultData) && isfinite(functionResultMC)){
      weight = functionResultData/functionResultMC;
      if ( kCaseGen == 3){   // never true ?
        if (PDGCode ==  111){
      if (!(fDoReweightHistoMCPi0 && hReweightMCHistPi0!= 0x0 && PDGCode ==  111)){
        weight = 1.;
      }
        }
        if (PDGCode ==  221){
      if (!(fDoReweightHistoMCEta && hReweightMCHistEta!= 0x0 && PDGCode ==  221)){
        weight = 1.;
      }
        }
      }
      if (!isfinite(functionResultData)) weight = 1.;
      if (!isfinite(weight)) weight = 1.;
    }
  } else if (PDGCode ==  310 && functionResultMC != 0 && isfinite(functionResultMC)){
    weight = functionResultMC;
  }

  return weight;
}


//_________________________________________________________________________
Float_t AliConvEventCuts::GetWeightForGamma(Int_t index, Double_t gammaPTrec, AliMCEvent *mcEvent, AliVEvent *event){
  // Gamma pT weighting for the material budget weights

  if(index < 0) return 0; // No Particle

  // check if MC production should be weighted. If it is with added particles check that particle is not rejected
  Int_t kCaseGen = 0;
  if ( fPeriodEnum == kLHC16NomB || fPeriodEnum == kLHC16P1Pyt8 || fPeriodEnum == kLHC16P1PHO ||
    fPeriodEnum == kLHC17pq  ||  fPeriodEnum == kLHC17P1PHO  || fPeriodEnum == kLHC17l3b  )
    kCaseGen = 2;  // regular MC


  if (kCaseGen == 0) return 1.;
  if (kCaseGen == 1 && !IsParticleFromBGEvent(index, mcEvent, event)) return 1.;

  // get pdg code
  Int_t PDGCode = 0;

  if(!event || event->IsA()==AliESDEvent::Class()){
    PDGCode = ((TParticle*)mcEvent->Particle(index))->GetPdgCode();

  } else if(event->IsA()==AliAODEvent::Class()){
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (fAODMCTrackArray){
      AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
      PDGCode = aodMCParticle->GetPdgCode();
    } else {
      return 1;
    }
  }

  // get MC value
  Float_t functionResultMC = 1.;
  if ( PDGCode == 22 && fDoReweightHistoMCGamma && hReweightMCHistGamma!= 0x0){
    //    cout << "gammaPt,PDGCode::"<< gammaPt<< "  " << PDGCode << endl;
    // AM 20.06.29 replaced Interpolate by BinContent
    // Use reconstructed pT as this is the momentum used to calculate the gamma pT weights;

    // functionResultMC = hReweightMCHistGamma->Interpolate(gammaPTrec);
    functionResultMC = hReweightMCHistGamma->GetBinContent(hReweightMCHistGamma->GetXaxis()->FindBin(gammaPTrec+0.001));
    // cout<< "functionResultMC::"<< functionResultMC  << " myBinMC::"<< hReweightMCHistGamma->GetXaxis()->FindBin(gammaPTrec+0.001)<< endl;
  }


  // get data value
  Float_t functionResultData = 1.;
  if ( PDGCode ==  22 && fDoReweightHistoMCGamma && hReweightDataHistGamma!= 0x0){
    // AM 20.06.29 replaced Interpolate by BinContent
    // Use reconstructed pT as this is the momenta used to calculate the gamma pT weights;

    // functionResultData = hReweightDataHistGamma->Interpolate(gammaPTrec);
    functionResultData = hReweightDataHistGamma->GetBinContent(hReweightDataHistGamma->GetXaxis()->FindBin(gammaPTrec+0.001));
    // cout<< "functionResultData::"<< functionResultData << " myBinData::"<< hReweightDataHistGamma->GetXaxis()->FindBin(gammaPTrec+0.001)<<endl;
  }


  // calculate weight from data and MC
  Double_t weight = 1.;
  if (PDGCode ==  22 ){
    if (functionResultData != 0. && functionResultMC != 0. && isfinite(functionResultData) && isfinite(functionResultMC)){
      weight = functionResultData/functionResultMC;
      //      cout << "weight::"<< weight << "  "   <<  endl;
      // if ( kCaseGen == 3){   // never true ?
      //   if (PDGCode ==  22){
      // 	  if (!(fDoReweightHistoMCGamma && hReweightMCHistGamma!= 0x0 && PDGCode ==  22)){
      // 	    weight = 1.;
      // 	  }
      //   }
      // }
      if (!isfinite(functionResultData)) weight = 1.;
      if (!isfinite(weight)) weight = 1.;
    }
  }

  return weight;
}


///________________________________________________________________________
void AliConvEventCuts::GetCorrectEtaShiftFromPeriod(){

  if( fPeriodEnum == kLHC13bc ||                                      // mainly minimum bias
      fPeriodEnum == kLHC13de ||                                      // mainly triggered
      fPeriodEnum == kLHC13b4_fix || fPeriodEnum == kLHC13b4_plus ||  // MC Pythia 6 (Jet-Jet), anchor LHC13b-e
      fPeriodEnum == kLHC13b2_efix ||                                 //MC DPMJET, anchr LHC13b+c
      fPeriodEnum == kLHC13e7 ||                                      //MC HIJING, anchr LHC13b+c
      fPeriodEnum == kLHC14b2 ||                                      //MC HIJING, anchr LHC13b+c
      fPeriodEnum == kLHC18j5 ||                                      //MC HIJING, anchr LHC13b+c
      fPeriodEnum == kLHC19a4                                         //MC JJ, anchr LHC13b-f
    ){
      printf(" Gamma Conversion Cuts %s :: pPb Run doing Eta Shift of %f \n\n",(GetCutNumber()).Data(),-0.465);
      SetEtaShift(-0.465);
  } else if(  fPeriodEnum == kLHC13f ) {
      printf(" Gamma Conversion Cuts %s :: Pbp Run doing Eta Shift of %f \n\n",(GetCutNumber()).Data(),0.465);
      SetEtaShift(+0.465);
    }
  else printf(" Gamma Conversion Cuts %s :: Automatic Eta Shift requested but Period is not known -> No Shift \n\n",(GetCutNumber()).Data());
}

//________________________________________________________________________
AliEMCALTriggerPatchInfo* AliConvEventCuts::GetMainTriggerPatch()
{
  //get main trigger match; if not known yet, look for it and cache

  if (fMainTriggerPatchEMCAL)
    return fMainTriggerPatchEMCAL;

  if (!fTriggerPatchInfo) {
    AliError(Form("%s: fTriggerPatchInfo not available",GetName()));
    return 0;
  }

  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntries();

  //extract main trigger patch
  AliEMCALTriggerPatchInfo *patch;
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
    patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    if (patch->IsMainTrigger()) {
      fMainTriggerPatchEMCAL = patch;
      break;
    }
  }

  return fMainTriggerPatchEMCAL;
}


//________________________________________________________________________
void AliConvEventCuts::InitializeEMCALTrigger(AliVEvent *event)
{
  //   cout << "entered EMCAL trigger initialization" << endl;

  // Init the analysis.
  if (fCaloTriggersName.IsNull()){
    if (event->IsA()==AliESDEvent::Class()){
      fCaloTriggersName = "EMCALTrigger";
    } else {
      fCaloTriggersName = "emcalTrigger";
    }
  }

  if (!fCaloTriggersName.IsNull() && !fCaloTriggers) {
    fCaloTriggers =  dynamic_cast<AliVCaloTrigger*>(event->FindListObject(fCaloTriggersName));
    if (!fCaloTriggers) {
      AliError(Form("%s: Could not retrieve calo triggers %s!", GetName(), fCaloTriggersName.Data()));
    return;
    }
  }

  if (fCaloTriggerPatchInfoName.IsNull()){
    if (event->IsA()==AliESDEvent::Class()){
      fCaloTriggerPatchInfoName = "EmcalTriggers";
    } else {
      fCaloTriggerPatchInfoName = "EmcalTriggers";
    }
  }

  if (!fCaloTriggerPatchInfoName.IsNull() && !fTriggerPatchInfo) {
    fTriggerPatchInfo = GetArrayFromEvent(event, fCaloTriggerPatchInfoName.Data(), "AliEMCALTriggerPatchInfo");
    if (!fTriggerPatchInfo) {
      AliError(Form("%s: Could not retrieve calo trigger patch info %s!", GetName(), fCaloTriggerPatchInfoName.Data()));
    return;
    }

  }

  fEMCALTrigInitialized = kTRUE;
}

//________________________________________________________________________
ULong_t AliConvEventCuts::GetTriggerList(){
  if (!fTriggerPatchInfo)
  return 0;
  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntries();

  //loop over patches to define trigger type of event
  Int_t nG1 = 0;
  Int_t nG2 = 0;
  Int_t nJ1 = 0;
  Int_t nJ2 = 0;
  Int_t nL0 = 0;
  AliEMCALTriggerPatchInfo *patch;
  // if (nPatch> 0) {cout << "NEW Triggers in this event*********************************" << endl;}
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
    patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    // cout << "Patch energy: "<<patch->GetPatchE() << "\t ADC counts: " << patch->GetADCAmp() << endl;
    // cout << "Phi: " << patch->GetPhiMin() << " - " << patch->GetPhiMax() << " delta phi: " <<TMath::Abs(patch->GetPhiMin()-patch->GetPhiMax())<< endl;
    // cout << "Eta: " << patch->GetEtaMin() << " - " << patch->GetEtaMax() << " delta eta: " <<TMath::Abs(patch->GetEtaMin()-patch->GetEtaMax())<< endl;
    if (patch->IsGammaHigh()){
      // cout << "fired L1GA high" << endl;
      nG1++;
    }
    if (patch->IsGammaLow()){
      // cout << "fired L1GA low" << endl;
      nG2++;
    }
    if (patch->IsJetHigh()){
      // cout << "fired L1JE high" << endl;
      nJ1++;
    }
    if (patch->IsJetLow()){
      // cout << "fired L1JE low" << endl;
      nJ2++;
    }
    if (patch->IsLevel0()){
      // cout << "fired L0" << endl;
      nL0++;
    }
    // cout << patch->GetPatchE()   << "\t" << patch->GetADCAmp()  << "\t" << patch->IsGammaHigh() << "\t" << patch->IsGammaLow()
        //  << "\t" << patch->IsJetHigh()  << "\t" << patch->IsJetLow()  << "\t" << patch->IsLevel0()
      //  << "\t" << patch->GetPhiMin()  << "\t" << patch->GetPhiMax()  << "\t" << TMath::Abs(patch->GetPhiMin()-patch->GetPhiMax())
      //  << "\t" << patch->GetEtaMin()  << "\t" << patch->GetEtaMax()  << "\t" << TMath::Abs(patch->GetEtaMin()-patch->GetEtaMax()) << endl;
  }

  if (nPatch > 0){
    AliDebug(2, "Patch summary: ");
    AliDebug(2, Form("Number of patches: %d", nPatch));
    AliDebug(2, Form("Level0: [%d]" ,nL0));
    AliDebug(2, Form("Jet:    low[%d], high[%d]" ,nJ2, nJ1));
    AliDebug(2, Form("Gamma:  low[%d], high[%d]" ,nG2, nG1));
  }

  // if (nPatch > 0){
  //   cout <<     Form("Number of patches: %d", nPatch) << endl;
  //   cout <<     Form("Level0: [%d]" ,nL0) << endl;
  //   cout <<     Form("Jet:    low[%d], high[%d]" ,nJ2, nJ1) << endl;
  //   cout <<     Form("Gamma:  low[%d], high[%d]" ,nG2, nG1) << endl;
  // }

  ULong_t triggers(0);
  if (nG1>0)
    SETBIT(triggers, kG1);
  if (nG2>0)
    SETBIT(triggers, kG2);
  if (nJ1>0)
    SETBIT(triggers, kJ1);
  if (nJ2>0)
    SETBIT(triggers, kJ2);
  if (nL0>0)
    SETBIT(triggers, kL0);
  return triggers;
}

//________________________________________________________________________
Bool_t AliConvEventCuts::HasTriggerType(TriggerTypeEMCAL t){
  // Check if event has a given trigger type
  if(t == kND){
    return fTriggersEMCAL == 0;
  }
  return TESTBIT(fTriggersEMCAL, int(t));
}


//________________________________________________________________________
TClonesArray *AliConvEventCuts::GetArrayFromEvent(AliVEvent* event, const char *name, const char *clname)
{
  // Get array from event.

  TClonesArray *arr = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    arr = dynamic_cast<TClonesArray*>(event->FindListObject(sname));
    if (!arr) {
    AliWarning(Form("%s: Could not retrieve array with name %s!", GetName(), name));
    return 0;
    }
  } else {
    return 0;
  }

  if (!clname)
    return arr;

  TString objname(arr->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom(clname)) {
    AliWarning(Form("%s: Objects of type %s in %s are not inherited from %s!",
            GetName(), cls.GetName(), name, clname));
    return 0;
  }
  return arr;
}

//_________________________________________________________________________
Bool_t AliConvEventCuts::IsConversionPrimaryESD( AliMCEvent *mcEvent, Long_t eventpos, Double_t prodVtxX, Double_t prodVtxY, Double_t prodVtxZ){

  if (eventpos < 0) return kFALSE;
  TParticle* particle = (TParticle *)mcEvent->Particle(eventpos);
  if (!particle) return kFALSE;
  if (TMath::Abs(particle->GetPdgCode()) == 11 ){
    if (particle->GetMother(0) != -1){
      TParticle* particleMother = (TParticle *)mcEvent->Particle(particle->GetMother(0));
      if (particleMother){
        if (TMath::Abs(particleMother->GetPdgCode()) == 22)
          particle = particleMother;
      }
    }
  }
  if (particle->GetMother(0) != -1){
    Double_t deltaX = particle->Vx() - prodVtxX;
    Double_t deltaY = particle->Vy() - prodVtxY;
    Double_t deltaZ = particle->Vz() - prodVtxZ;

    //Double_t realRadius2D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY);
    Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);


    Bool_t dalitzCand = kFALSE;

    TParticle* firstmother = (TParticle *)mcEvent->Particle(particle->GetMother(0));
    if (!firstmother) return kFALSE;
    Int_t pdgCodeFirstMother     = firstmother->GetPdgCode();
    Bool_t intDecay = kFALSE;
    if ( pdgCodeFirstMother == 111 || pdgCodeFirstMother == 221 ) intDecay = kTRUE;
    if ( intDecay && TMath::Abs(particle->GetPdgCode()) == 11 ){
      dalitzCand = kTRUE;
      // cout << "dalitz candidate found" << endl;
    }

    Long_t source = particle->GetMother(0);
    Bool_t foundExcludedPart = kFALSE;
    Bool_t foundShower = kFALSE;
    Int_t pdgCodeMotherPrev = 0;
    Int_t pdgCodeMotherPPrevMother = 0;
    Int_t depth = 0;
    if (dalitzCand || realRadius3D < fSecProdBoundary ){
      // if (particle->GetPdgCode() == 22){
      //   cout << endl << endl << "new particle: " << eventpos <<endl;
      //   cout << particle->GetPdgCode() << "\t" << particle->R() << "\t" << realRadius2D << "\t" << realRadius3D << endl;
      // }
      while (depth < 20){
        TParticle* mother   = (TParticle *)mcEvent->Particle(source);
        source         = mother->GetMother(0);
        // if (particle->GetPdgCode() == 22)cout << "eventposition: "<< source << endl;
        Int_t pdgCodeMother     = mother->GetPdgCode();
        // if (particle->GetPdgCode() == 22)cout << "Previous mothers: " << pdgCodeMother << "\t"<< pdgCodeMotherPrev<< "\t" << pdgCodeMotherPPrevMother << endl;
        if (pdgCodeMother == pdgCodeMotherPrev && pdgCodeMother == pdgCodeMotherPPrevMother) depth = 20;
        if (TMath::Abs(pdgCodeMother) == 11 && TMath::Abs(pdgCodeMotherPrev) == 22 && TMath::Abs(pdgCodeMotherPPrevMother) == 11 ){
          foundShower = kTRUE;
          depth =20;
        }
        if (TMath::Abs(pdgCodeMother) == 22 && TMath::Abs(pdgCodeMotherPrev) == 11 && TMath::Abs(pdgCodeMotherPPrevMother) == 22 ){
          foundShower = kTRUE;
          depth =20;
        }

        // particles to be excluded:
        // K0s     - 310
        // K0l     - 130
        // K+/-    - 321
        // Lambda  - 3122
        // Sigma0  - 3212
        // Sigma+/-  - 3222, 3112
        // Cascades  - 3322, 3312
        if (TMath::Abs(pdgCodeMother) == 310   || TMath::Abs(pdgCodeMother) == 130   || TMath::Abs(pdgCodeMother) == 321  ||
          TMath::Abs(pdgCodeMother) == 3122   || TMath::Abs(pdgCodeMother) == 3212   || TMath::Abs(pdgCodeMother) == 3222 ||
          TMath::Abs(pdgCodeMother) == 3112   || TMath::Abs(pdgCodeMother) == 3322   || TMath::Abs(pdgCodeMother) == 3312
        ) {
          foundExcludedPart = kTRUE;
        }
        // if (particle->GetPdgCode() == 22)cout << mother->GetPdgCode() << "\t" <<  source << "\t" << foundExcludedPart<< endl;
        pdgCodeMotherPPrevMother = pdgCodeMotherPrev;
        pdgCodeMotherPrev = pdgCodeMother;
        if (source == -1) depth = 20;

        // if (particle->GetPdgCode() == 22)cout << depth << endl;
        depth++;
      }
    }
    if (foundExcludedPart){
      // if (particle->GetPdgCode() == 22)cout << "This is definitely a secondary, manually excluded" << endl;
      return kFALSE;
    } else if (dalitzCand && realRadius3D < fSecProdBoundary ){
      // if (particle->GetPdgCode() == 22)cout << "This was a decay via a virtual photon" << endl;
      return kTRUE;
    } else if (foundShower){
      // if (particle->GetPdgCode() == 22)cout << "This is a shower" << endl;
      return kFALSE;
    } else if (realRadius3D >= fSecProdBoundary){
      // cout << "This is a secondary, to large production radius" << endl;
      return kFALSE;
    }
  }

  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliConvEventCuts::IsConversionPrimaryAOD(AliVEvent *event, AliAODMCParticle* AODMCParticle,  Double_t prodVtxX, Double_t prodVtxY, Double_t prodVtxZ){

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return kFALSE;
  AliAODMCParticle* currentParticle = AODMCParticle;
  if (TMath::Abs(currentParticle->GetPdgCode()) == 11 ){
    if (currentParticle->GetMother() != -1){
      AliAODMCParticle* particleMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(currentParticle->GetMother()));
      if (particleMother){
        if (TMath::Abs(particleMother->GetPdgCode()) == 22)
          currentParticle = particleMother;
      }
    }
  }
  if (currentParticle->GetMother() > -1){
    Double_t deltaX = currentParticle->Xv() - prodVtxX;
    Double_t deltaY = currentParticle->Yv() - prodVtxY;
    Double_t deltaZ = currentParticle->Zv() - prodVtxZ;

    //Double_t realRadius2D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY);
    Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);

    Bool_t dalitzCand = kFALSE;

    AliAODMCParticle* firstmother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(currentParticle->GetMother()));
    if (!firstmother) return kFALSE;
    Int_t pdgCodeFirstMother = firstmother->GetPdgCode();
    Bool_t intDecay = kFALSE;
    if ( pdgCodeFirstMother == 111 || pdgCodeFirstMother == 221 ) intDecay = kTRUE;
    if ( intDecay && TMath::Abs(currentParticle->GetPdgCode()) == 11 ){
      dalitzCand = kTRUE;
      // cout << "dalitz candidate found" << endl;
    }

    Long_t source = currentParticle->GetMother();
    Bool_t foundExcludedPart = kFALSE;
    Bool_t foundShower = kFALSE;
    Int_t pdgCodeMotherPrev = 0;
    Int_t pdgCodeMotherPPrevMother = 0;
    Int_t depth = 0;
    if (dalitzCand || realRadius3D < fSecProdBoundary ){
      // if (currentParticle->GetPdgCode() == 22){
      //   cout << endl << endl << "new particle: " << eventpos <<endl;
      //   cout << currentParticle->GetPdgCode() << "\t" << currentParticle->R() << "\t" << realRadius2D << "\t" << realRadius3D << endl;
      // }
      while (depth < 20){
        AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(source));
        source = mother->GetMother();
        // if (currentParticle->GetPdgCode() == 22)cout << "eventposition: "<< source << endl;
        Int_t pdgCodeMother     = mother->GetPdgCode();
        // if (currentParticle->GetPdgCode() == 22)cout << "Previous mothers: " << pdgCodeMother << "\t"<< pdgCodeMotherPrev<< "\t" << pdgCodeMotherPPrevMother << endl;
        if (pdgCodeMother == pdgCodeMotherPrev && pdgCodeMother == pdgCodeMotherPPrevMother) depth = 20;
        if (TMath::Abs(pdgCodeMother) == 11 && TMath::Abs(pdgCodeMotherPrev) == 22 && TMath::Abs(pdgCodeMotherPPrevMother) == 11 ){
          foundShower = kTRUE;
          depth =20;
        }
        if (TMath::Abs(pdgCodeMother) == 22 && TMath::Abs(pdgCodeMotherPrev) == 11 && TMath::Abs(pdgCodeMotherPPrevMother) == 22 ){
          foundShower = kTRUE;
          depth =20;
        }

        // particles to be excluded:
        // K0s      - 310
        // K0l      - 130
        // K+/-     - 321
        // Lambda   - 3122
        // Sigma0   - 3212
        // Sigma+/- - 3222, 3112
        // Cascades - 3322, 3312
        if (TMath::Abs(pdgCodeMother) == 310   || TMath::Abs(pdgCodeMother) == 130   || TMath::Abs(pdgCodeMother) == 321  ||
          TMath::Abs(pdgCodeMother) == 3122   || TMath::Abs(pdgCodeMother) == 3212   || TMath::Abs(pdgCodeMother) == 3222 ||
          TMath::Abs(pdgCodeMother) == 3112   || TMath::Abs(pdgCodeMother) == 3322   || TMath::Abs(pdgCodeMother) == 3312)
        {
          foundExcludedPart = kTRUE;
        }
    //  if (currentParticle->GetPdgCode() == 22)cout << mother->GetPdgCode() << "\t" <<  source << "\t" << foundExcludedPart<< endl;
        pdgCodeMotherPPrevMother = pdgCodeMotherPrev;
        pdgCodeMotherPrev = pdgCodeMother;
        if (source == -1) depth = 20;

    //  if (currentParticle->GetPdgCode() == 22)cout << depth << endl;
        depth++;
      }
    }
    if (foundExcludedPart){
    //  if (currentParticle->GetPdgCode() == 22)cout << "This is definitely a secondary, manually excluded" << endl;
      return kFALSE;
    } else if (dalitzCand && realRadius3D < fSecProdBoundary ){
    //  if (currentParticle->GetPdgCode() == 22)cout << "This was a decay via a virtual photon" << endl;
      return kTRUE;
    } else if (foundShower){
    //  if (currentParticle->GetPdgCode() == 22)cout << "This is a shower" << endl;
      return kFALSE;
    } else if (realRadius3D >= fSecProdBoundary){
    //  cout << "This is a secondary, too large production radius" << endl;
      return kFALSE;
    }
  }

  return kTRUE;
}


//________________________________________________________________________
Int_t AliConvEventCuts::SecondaryClassificationPhoton( TParticle *particle, AliMCEvent* mcEvent, Bool_t isConversion ){
  if (particle != NULL && mcEvent != NULL){
    Int_t pdgSecondary      = 0;
    if (!isConversion){
      //Bool_t hasMother        = kFALSE;
      //Bool_t hasGrandMother   = kFALSE;
      Long_t motherID         = particle->GetMother(0);
      Long_t grandMotherID    = -1;
      // is the photon a direct photons, without a mother?
      if (motherID > -1){
        //hasMother             = kTRUE;
        grandMotherID         = mcEvent->Particle(motherID)->GetMother(0);
        // is the meson a primary?
        if (grandMotherID > -1){
          //hasGrandMother      = kTRUE;
          pdgSecondary        = mcEvent->Particle(grandMotherID)->GetPdgCode();
        }
      }
    } else {
      //Bool_t hasMother            = kFALSE;
      //Bool_t hasGrandMother       = kFALSE;
      //Bool_t hasGreatGrandMother  = kFALSE;
      Long_t motherID             = particle->GetMother(0);
      Long_t grandMotherID        = -1;
      Long_t greatGrandMotherID   = -1;
      // is the electron a direct electron, without a mother?
      if (motherID > -1){
        //hasMother                 = kTRUE;
        grandMotherID             = mcEvent->Particle(motherID)->GetMother(0);
        // is the photon a direct photons, without a mother?
        if (grandMotherID > -1){
          //hasGrandMother          = kTRUE;
          greatGrandMotherID      = mcEvent->Particle(grandMotherID)->GetMother(0);
          // is the meson a primary?
          if (greatGrandMotherID > -1){
            //hasGreatGrandMother   = kTRUE;
            pdgSecondary          = mcEvent->Particle(greatGrandMotherID)->GetPdgCode();
          }
        }
      }
    }
    // is the secondary photon from a lambda
    if (TMath::Abs(pdgSecondary) == 3122 )
      return 3;
    // is the secondary photon from a K0s
    else if ( TMath::Abs(pdgSecondary) == 310 )
      return 2;
    // is the secondary photon from a K0l
    else if ( TMath::Abs(pdgSecondary) == 130 )
      return 5;
    // is the secondary photon from a eta
    else if ( TMath::Abs(pdgSecondary) == 221 )
      return 4;
    // is the secondary photon from something else
    else if ( TMath::Abs(pdgSecondary) != 0 )
      return 1;

  }

  return 0;
}

//________________________________________________________________________
Int_t AliConvEventCuts::SecondaryClassificationPhotonAOD( AliAODMCParticle *particle, TClonesArray *aodmcArray, Bool_t isConversion ){
  if (particle != NULL && aodmcArray != NULL){
    Int_t pdgSecondary      = 0;
    if (!isConversion){
      //Bool_t hasMother        = kFALSE;
      //Bool_t hasGrandMother   = kFALSE;
      Long_t motherID         = particle->GetMother();
      Long_t grandMotherID    = -1;
      // is the photon a direct photons, without a mother?
      if (motherID > -1){
        //hasMother             = kTRUE;
        grandMotherID         = ((AliAODMCParticle*)aodmcArray->At(motherID))->GetMother();
        // is the meson a primary?
        if (grandMotherID > -1){
          //hasGrandMother      = kTRUE;
          pdgSecondary        = ((AliAODMCParticle*)aodmcArray->At(grandMotherID))->GetPdgCode();
        }
      }
    } else {
      //Bool_t hasMother            = kFALSE;
      //Bool_t hasGrandMother       = kFALSE;
      //Bool_t hasGreatGrandMother  = kFALSE;
      Long_t motherID             = particle->GetMother();
      Long_t grandMotherID        = -1;
      Long_t greatGrandMotherID   = -1;
      // is the electron a direct electron, without a mother?
      if (motherID > -1){
        //hasMother                 = kTRUE;
        grandMotherID             = ((AliAODMCParticle*)aodmcArray->At(motherID))->GetMother();
        // is the photon a direct photons, without a mother?
        if (grandMotherID > -1){
          //hasGrandMother          = kTRUE;
          greatGrandMotherID      = ((AliAODMCParticle*)aodmcArray->At(grandMotherID))->GetMother();
          // is the meson a primary?
          if (greatGrandMotherID > -1){
            //hasGreatGrandMother   = kTRUE;
            pdgSecondary          = ((AliAODMCParticle*)aodmcArray->At(greatGrandMotherID))->GetPdgCode();
          }
        }
      }
    }
    // is the secondary photon from a lambda
    if (TMath::Abs(pdgSecondary) == 3122 )
      return 3;
    // is the secondary photon from a K0s
    else if ( TMath::Abs(pdgSecondary) == 310 )
      return 2;
    // is the secondary photon from a K0l
    else if ( TMath::Abs(pdgSecondary) == 130 )
      return 5;
    // is the secondary photon from a eta
    else if ( TMath::Abs(pdgSecondary) == 221 )
      return 4;
    // is the secondary photon from something else
    else if ( TMath::Abs(pdgSecondary) != 0 )
      return 1;

  }

  return 0;
}

void AliConvEventCuts::SetPeriodEnum (TString periodName){

  if (periodName.CompareTo("") == 0){
    periodName = GetV0Reader()->GetPeriodName();
  }
  if (periodName.CompareTo("") == 0) {
    fPeriodEnum = kNoPeriod;
    fEnergyEnum = kUnset;
    AliError("No correct period could be set, periodName string empty");
    return;
  }

  // Data
  if (periodName.CompareTo("LHC10b") == 0 || periodName.CompareTo("LHC10c") == 0 || periodName.CompareTo("LHC10d") == 0 || periodName.CompareTo("LHC10e") == 0 ||
      periodName.CompareTo("LHC10f") == 0 || periodName.CompareTo("LHC10g") == 0 || periodName.CompareTo("LHC10bg") == 0
  ){
    fPeriodEnum = kLHC10bg;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10h") == 0) {
    fPeriodEnum = kLHC10h;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC11a") == 0) {
    fPeriodEnum = kLHC11a;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC11b") == 0) {
    fPeriodEnum = kLHC11b;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC11c") == 0 || periodName.CompareTo("LHC11d") == 0 || periodName.CompareTo("LHC11e") == 0 || periodName.CompareTo("LHC11f") == 0 ||
            periodName.CompareTo("LHC11g") == 0
  ) {
    fPeriodEnum = kLHC11cg;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC11h") == 0) {
    fPeriodEnum = kLHC11h;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a") == 0 || periodName.CompareTo("LHC12b") == 0 || periodName.CompareTo("LHC12c") == 0 || periodName.CompareTo("LHC12d") == 0 ||
            periodName.CompareTo("LHC12e") == 0 || periodName.CompareTo("LHC12f") == 0 || periodName.CompareTo("LHC12g") == 0 || periodName.CompareTo("LHC12h") == 0 ||
            periodName.CompareTo("LHC12i") == 0 || periodName.CompareTo("LHC12ai") == 0
  ) {
    fPeriodEnum = kLHC12;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC13b") == 0 || periodName.CompareTo("LHC13c") == 0 || periodName.CompareTo("LHC13bc") == 0){
    fPeriodEnum = kLHC13bc;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13d") == 0 || periodName.CompareTo("LHC13e") == 0 || periodName.CompareTo("LHC13de") == 0){
    fPeriodEnum = kLHC13de;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13f") == 0 ){
    fPeriodEnum = kLHC13f;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13g") == 0 ){
    fPeriodEnum = kLHC13g;
    fEnergyEnum = k2760GeV;
  } else if ( periodName.CompareTo("LHC15f") == 0 || periodName.CompareTo("LHC15g") == 0 || periodName.CompareTo("LHC15h") == 0 || periodName.CompareTo("LHC15i") == 0 ||
              periodName.CompareTo("LHC15j") == 0 || periodName.CompareTo("LHC15k") == 0 || periodName.CompareTo("LHC15l") == 0 || periodName.CompareTo("LHC15m") == 0 ||
              periodName.CompareTo("LHC15fm") == 0
  ) {
    fPeriodEnum = kLHC15fm;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC15n") == 0 ){
    fPeriodEnum = kLHC15n;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC15o") == 0 ){
    fPeriodEnum = kLHC15o;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC16dp") == 0 || periodName.CompareTo("LHC16d") == 0 || periodName.CompareTo("LHC16e") == 0 || periodName.CompareTo("LHC16g") == 0 ||
              periodName.CompareTo("LHC16h") == 0 || periodName.CompareTo("LHC16i") == 0 || periodName.CompareTo("LHC16j") == 0 || periodName.CompareTo("LHC16k") == 0 || periodName.CompareTo("LHC16l") == 0 || periodName.CompareTo("LHC16m") == 0 || periodName.CompareTo("LHC16n") == 0 || periodName.CompareTo("LHC16o") == 0 ||
              periodName.CompareTo("LHC16p") == 0){
    fPeriodEnum = kLHC16NomB;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16f") == 0 ){
    fPeriodEnum = kLHC16LowB;
    fEnergyEnum = k13TeVLowB;
  } else if (periodName.CompareTo("LHC16qt") == 0  || periodName.CompareTo("LHC16q") == 0 ||  periodName.CompareTo("LHC16t") == 0 ){
    fPeriodEnum = kLHC16qt;
    fEnergyEnum = kpPb5TeVR2;
  } else if (periodName.CompareTo("LHC16r") == 0 ){
    fPeriodEnum = kLHC16r;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC16s") == 0 ){
    fPeriodEnum = kLHC16s;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17cr") == 0 ||  periodName.CompareTo("LHC17c") == 0 || periodName.CompareTo("LHC17d") == 0 || periodName.CompareTo("LHC17e") == 0 ||
             periodName.CompareTo("LHC17f") == 0 || periodName.CompareTo("LHC17h") == 0 || periodName.CompareTo("LHC17i") == 0 || periodName.CompareTo("LHC17j") == 0 ||
             periodName.CompareTo("LHC17k") == 0 || periodName.CompareTo("LHC17l") == 0 || periodName.CompareTo("LHC17m") == 0 || periodName.CompareTo("LHC17o") == 0 ||
             periodName.CompareTo("LHC17r") == 0 ){
    fPeriodEnum = kLHC17NomB;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17g") == 0 ){
    fPeriodEnum = kLHC17LowB;
    fEnergyEnum = k13TeVLowB;
  } else if ( periodName.CompareTo("LHC17n") == 0 ){
    fPeriodEnum = kLHC17n;
    fEnergyEnum = kXeXe5440GeV;
  } else if ( periodName.Contains("LHC17pq")  || periodName.Contains("LHC17p") || periodName.Contains("LHC17q")){
    fPeriodEnum = kLHC17pq;
    fEnergyEnum = k5TeV;
  } else if ( periodName.CompareTo("LHC18b") == 0 || periodName.CompareTo("LHC18d") == 0 || periodName.CompareTo("LHC18e") == 0 || periodName.CompareTo("LHC18f") == 0 ||
              periodName.CompareTo("LHC18g") == 0 || periodName.CompareTo("LHC18h") == 0 || periodName.CompareTo("LHC18i") == 0 || periodName.CompareTo("LHC18j") == 0 ){
    fPeriodEnum = kLHC18NomB;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC18c") == 0 ){
    fPeriodEnum = kLHC18LowB;
    fEnergyEnum = k13TeVLowB;
  } else if ( periodName.CompareTo("LHC18q") == 0 || periodName.CompareTo("LHC18r") == 0 || periodName.CompareTo("LHC18qr") == 0){
    fPeriodEnum = kLHC18qr;
    fEnergyEnum = kPbPb5TeV;

  // LHC10x anchored MCs
  } else if (periodName.CompareTo("LHC10d1") == 0){
    fPeriodEnum = kLHC10d1;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10d2") == 0){
    fPeriodEnum = kLHC10d2;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10d4a") == 0){
    fPeriodEnum = kLHC10d4a;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10d4") == 0){
    fPeriodEnum = kLHC10d4;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10e12") == 0){
    fPeriodEnum = kLHC10e12;
    fEnergyEnum = k900GeV;
  } else if (periodName.CompareTo("LHC10e13") == 0){
    fPeriodEnum = kLHC10e13;
    fEnergyEnum = k900GeV;
  } else if (periodName.CompareTo("LHC10e20") == 0){
    fPeriodEnum = kLHC10e20;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10e21") == 0){
    fPeriodEnum = kLHC10e21;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10f6a") == 0){
    fPeriodEnum = kLHC10f6a;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC10f6") == 0){
    fPeriodEnum = kLHC10f6;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC14b7") == 0){
    fPeriodEnum = kLHC14b7;
    fEnergyEnum = k7TeV;
  } else if (periodName.Contains("LHC14j4")){
    fPeriodEnum = kLHC14j4;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC13d2") == 0){
    fPeriodEnum = kLHC13d2;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC13d2b") == 0){
    fPeriodEnum = kLHC13d2b;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a11a") == 0){
    fPeriodEnum = kLHC12a11a;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a11b") == 0){
    fPeriodEnum = kLHC12a11b;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a11c") == 0){
    fPeriodEnum = kLHC12a11c;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a11d") == 0){
    fPeriodEnum = kLHC12a11d;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a11e") == 0){
    fPeriodEnum = kLHC12a11e;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC12a11f") == 0){
    fPeriodEnum = kLHC12a11f;
    fEnergyEnum = kPbPb2760GeV;
  // LHC11x anchored MCs
  } else if (periodName.CompareTo("LHC12a15c") == 0){
    fPeriodEnum = kLHC12a15c;
    fEnergyEnum = k2760GeV;
  } else if (periodName.Contains("LHC12f1a") ){
    fPeriodEnum = kLHC12f1a;
    fEnergyEnum = k2760GeV;
  } else if (periodName.Contains("LHC12f1b") ){
    fPeriodEnum = kLHC12f1b;
    fEnergyEnum = k2760GeV;
  } else if (periodName.Contains("LHC12i3") ){
    fPeriodEnum = kLHC12i3;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15g1a") == 0){
    fPeriodEnum = kLHC15g1a;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15g1b") == 0){
    fPeriodEnum = kLHC15g1b;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC13e4") == 0){
    fPeriodEnum = kLHC13e4;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC13e5") == 0){
    fPeriodEnum = kLHC13e5;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC14k1a") == 0){
    fPeriodEnum = kLHC14k1a;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC14k1b") == 0){
    fPeriodEnum = kLHC14k1b;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC12a15f") == 0){
    fPeriodEnum = kLHC12a15f;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC12a15g") == 0){
    fPeriodEnum = kLHC12a15g;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC12f2a") == 0){
    fPeriodEnum = kLHC12f2a;
    fEnergyEnum = k7TeV;
  } else if (periodName.CompareTo("LHC14a1a") == 0){
    fPeriodEnum = kLHC14a1a;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC14a1b") == 0){
    fPeriodEnum = kLHC14a1b;
    fEnergyEnum = kPbPb2760GeV;
  } else if (periodName.CompareTo("LHC14a1c") == 0){
    fPeriodEnum = kLHC14a1c;
    fEnergyEnum = kPbPb2760GeV;
  // LHC12x anchored MCs
  } else if (periodName.CompareTo("LHC14e2b") == 0){
    fPeriodEnum = kLHC14e2b;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC12P2Pyt8")==0 || periodName.Contains("LHC15h1")){
    fPeriodEnum = kLHC15h1;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC12P2Pho")==0 || periodName.Contains("LHC15h2")){
    fPeriodEnum = kLHC15h2;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC12P2JJ") == 0 || periodName.CompareTo("LHC16c2") == 0 || periodName.CompareTo("LHC16c2_plus") == 0){
    fPeriodEnum = kLHC12P2JJ;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC17g5a1") == 0){
    fPeriodEnum = kLHC17g5a1;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC17g5a2") == 0){
    fPeriodEnum = kLHC17g5a2;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC17g5b") == 0){
    fPeriodEnum = kLHC17g5b;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC17g5c") == 0){
    fPeriodEnum = kLHC17g5c;
    fEnergyEnum = k8TeV;
  // LHC13x anchored MCs
  } else if (periodName.Contains("LHC13b2_efix")){
    fPeriodEnum = kLHC13b2_efix;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13e7") == 0){
    fPeriodEnum = kLHC13e7;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC14b2") == 0){
    fPeriodEnum = kLHC14b2;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC18j5") == 0){
    fPeriodEnum = kLHC18j5;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13b4_fix") == 0){
    fPeriodEnum = kLHC13b4_fix;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13b4_plus") == 0){
    fPeriodEnum = kLHC13b4_plus;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.Contains("LHC19a4")){
    fPeriodEnum = kLHC19a4;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16c3a") == 0 || periodName.CompareTo("LHC16c3a2") == 0){
    fPeriodEnum = kLHC16c3a;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16c3b") == 0 || periodName.CompareTo("LHC16c3b2") == 0){
    fPeriodEnum = kLHC16c3b;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16c3c") == 0 || periodName.CompareTo("LHC16c3c2") == 0){
    fPeriodEnum = kLHC16c3c;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17g6a1") == 0 ){
    fPeriodEnum = kLHC17g6a1;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17g6a2") == 0 ){
    fPeriodEnum = kLHC17g6a2;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17g6a3") == 0 ){
    fPeriodEnum = kLHC17g6a3;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC15g2") == 0){
    fPeriodEnum = kLHC15g2;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15a3a") == 0){
    fPeriodEnum = kLHC15a3a;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15a3a_plus") == 0){
    fPeriodEnum = kLHC15a3a_plus;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15a3b") == 0){
    fPeriodEnum = kLHC15a3b;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15d3a") == 0){
    fPeriodEnum = kLHC15d3a;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15d3b") == 0){
    fPeriodEnum = kLHC15d3b;
    fEnergyEnum = k2760GeV;
  // LHC15x anchored MCs
  } else if (periodName.CompareTo("LHC15g3a3") == 0){
    fPeriodEnum = kLHC15g3a3;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC15g3a") == 0){
    fPeriodEnum = kLHC15g3a;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC15g3c2") == 0){
    fPeriodEnum = kLHC15g3c2;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC15g3c3") == 0){
    fPeriodEnum = kLHC15g3c3;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC15g3") == 0){
    fPeriodEnum = kLHC15g3;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16a2a") == 0){
    fPeriodEnum = kLHC16a2a;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16a2b") == 0){
    fPeriodEnum = kLHC16a2b;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16a2c") == 0){
    fPeriodEnum = kLHC16a2c;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16d3") == 0){
    fPeriodEnum = kLHC15P2EPos;
    fEnergyEnum = k13TeV; //LHC15f
  } else if (periodName.CompareTo("LHC17i4") == 0){
    fPeriodEnum = kLHC15P2Pyt8;
    fEnergyEnum = k13TeV; //LHC15i
  } else if (periodName.CompareTo("LHC17g7") == 0){
    fPeriodEnum = kLHC15P2Pyt8;
    fEnergyEnum = k13TeV; //LHC15h
  } else if (periodName.CompareTo("LHC15l1a2") == 0){
    fPeriodEnum = kLHC15l1a2;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC15l1b2") == 0){
    fPeriodEnum = kLHC15l1b2;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16h3") == 0){
    fPeriodEnum = kLHC16h3;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16h8a") == 0){
    fPeriodEnum = kLHC16h8a;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16h8b") == 0){
    fPeriodEnum = kLHC16h8b;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16k5a") == 0){
    fPeriodEnum = kLHC16k5a;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16k5b") == 0){
    fPeriodEnum = kLHC16k5b;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16k3a") == 0){
    fPeriodEnum = kLHC16k3a;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC16k3a2") == 0){
    fPeriodEnum = kLHC16k3a2;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC17e2") == 0){
    fPeriodEnum = kLHC17e2;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC18j3") == 0){
    fPeriodEnum = kLHC18j3;
    fEnergyEnum = k5TeV;
  } else if (periodName.Contains("LHC15k1a1")){
    fPeriodEnum = kLHC15k1a1;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC15k1a2")){
    fPeriodEnum = kLHC15k1a2;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC15k1a3")){
    fPeriodEnum = kLHC15k1a3;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16j7")){
    fPeriodEnum = kLHC16j7;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g2")){
    fPeriodEnum = kLHC16g2;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g3")){
    fPeriodEnum = kLHC16g3;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16h4")){
    fPeriodEnum = kLHC16h4;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i1a")){
    fPeriodEnum = kLHC16i1a;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i1b")){
    fPeriodEnum = kLHC16i1b;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i1c")){
    fPeriodEnum = kLHC16i1c;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i2a")){
    fPeriodEnum = kLHC16i2a;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i2b")){
    fPeriodEnum = kLHC16i2b;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i2c")){
    fPeriodEnum = kLHC16i2c;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i3a")){
    fPeriodEnum = kLHC16i3a;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i3b")){
    fPeriodEnum = kLHC16i3b;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16i3c")){
    fPeriodEnum = kLHC16i3c;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16h2a")){
    fPeriodEnum = kLHC16h2a;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16h2b")){
    fPeriodEnum = kLHC16h2b;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16h2c")){
    fPeriodEnum = kLHC16h2c;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.CompareTo("LHC16k3b") == 0){
    fPeriodEnum = kLHC16k3b;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.CompareTo("LHC16k3b2") == 0){
    fPeriodEnum = kLHC16k3b2;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC15k5a") == 0 || periodName.CompareTo("LHC15k5a2") == 0){
    fPeriodEnum = kLHC15k5a;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC15k5b") == 0 || periodName.CompareTo("LHC15k5b2") == 0){
    fPeriodEnum = kLHC15k5b;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC15k5c") == 0 || periodName.CompareTo("LHC15k5c2") == 0){
    fPeriodEnum = kLHC15k5c;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC18b11a") == 0){
    fPeriodEnum = kLHC18b11a;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18b11b") == 0){
    fPeriodEnum = kLHC18b11b;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18b11c") == 0){
    fPeriodEnum = kLHC18b11c;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18e1") == 0){
    fPeriodEnum = kLHC18e1;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18e1a") == 0){
    fPeriodEnum = kLHC18e1a;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18e1b") == 0){
    fPeriodEnum = kLHC18e1b;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18e1c") == 0){
    fPeriodEnum = kLHC18e1c;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18l8a") == 0 ){
    fPeriodEnum = kLHC18l8a;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18l8b") == 0 ){
    fPeriodEnum = kLHC18l8b;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC18l8c") == 0 ){
    fPeriodEnum = kLHC18l8c;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC19h2a") == 0 ){
    fPeriodEnum = kLHC19h2a;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC19h2b") == 0 ){
    fPeriodEnum = kLHC19h2b;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC19h2c") == 0 ){
    fPeriodEnum = kLHC19h2c;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC19h3") == 0 ){
    fPeriodEnum = kLHC19h3;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC20e3a") == 0 ){
    fPeriodEnum = kLHC20e3a;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC20e3b") == 0 ){
    fPeriodEnum = kLHC20e3b;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC20e3c") == 0 ){
    fPeriodEnum = kLHC20e3c;
    fEnergyEnum = kPbPb5TeV;
  } else if ( periodName.CompareTo("LHC20g10") == 0 ){
    fPeriodEnum = kLHC20g10;
    fEnergyEnum = kPbPb5TeV;



  // LHC16x anchored MCs
  // 13TeV LHC16* anchors full field Pythia 8 MB
  } else if ( periodName.CompareTo("LHC16P1Pyt8") == 0 ||
    periodName.CompareTo("LHC17f6") == 0 ||  periodName.CompareTo("LHC17f6_extra") == 0 ||
    periodName.CompareTo("LHC17d17") == 0 || periodName.CompareTo("LHC17d17_extra") == 0 ||
    periodName.CompareTo("LHC17f5") == 0 ||  periodName.CompareTo("LHC17f5_extra") == 0 ||
    periodName.CompareTo("LHC17d3") == 0 ||  periodName.CompareTo("LHC17d3_extra") == 0 ||
    periodName.CompareTo("LHC17e5") == 0 ||  periodName.CompareTo("LHC17e5_extra") == 0 ||
    periodName.CompareTo("LHC17d20a1") == 0 || periodName.CompareTo("LHC17d20a1_extra") == 0 ||
    periodName.CompareTo("LHC17d20a2") == 0 || periodName.CompareTo("LHC17d20a2_extra") == 0 ||
    periodName.CompareTo("LHC17d16") == 0 ||  periodName.CompareTo("LHC17d16_extra") == 0 ||
    periodName.CompareTo("LHC17d18") == 0 ||  periodName.CompareTo("LHC17d18_extra") == 0  ||
    periodName.CompareTo("LHC17f9") == 0 ||  periodName.CompareTo("LHC17f9_extra") == 0 ||
    periodName.CompareTo("LHC17f9_test") == 0 ||
    periodName.CompareTo("LHC18f1") == 0 || periodName.CompareTo("LHC18d8") == 0 ||
    periodName.CompareTo("LHC19g2") == 0 ||  periodName.CompareTo("LHC19g2a") == 0 ||
    periodName.CompareTo("LHC20f2") == 0    // pass2 test of ITS/TPC geometry changes
  ){
    fPeriodEnum = kLHC16P1Pyt8;
    fEnergyEnum = k13TeV;
  // 13TeV LHC16* anchors low field Pythia 8 MB
  } else if ( periodName.CompareTo("LHC16P1Pyt8LowB") == 0 || periodName.CompareTo("LHC17d1") == 0 ){
    fPeriodEnum = kLHC16P1Pyt8LowB;
    fEnergyEnum = k13TeVLowB;
  // 13TeV LHC16* anchors full field EPOS MB
  } else if ( periodName.CompareTo("LHC16P1EPOS") == 0 || periodName.CompareTo("LHC17d20b1") == 0 || periodName.CompareTo("LHC17d20b2") == 0 ){
      fPeriodEnum = kLHC16P1EPOS;
      fEnergyEnum = k13TeV;
  // 13TeV LHC16d anchors full field Phojet MB
  } else if ( periodName.CompareTo("LHC16P1PHO") == 0 || periodName.CompareTo("LHC18d6a") == 0 || periodName.CompareTo("LHC18d6a2") == 0  ){
      fPeriodEnum = kLHC16P1PHO;
      fEnergyEnum = k13TeV;
  // 13TeV LHC16* anchors full field JJ Pythia 8 MB
  } else if ( periodName.CompareTo("LHC16P1JJ") == 0 || periodName.CompareTo("LHC17f8a") == 0 || periodName.CompareTo("LHC17f8c") == 0 || periodName.CompareTo("LHC17f8d") == 0 ||
              periodName.CompareTo("LHC17f8e") == 0 || periodName.CompareTo("LHC17f8f") == 0 || periodName.CompareTo("LHC17f8g") == 0 || periodName.CompareTo("LHC17f8h") == 0 ||
              periodName.CompareTo("LHC17f8i") == 0 || periodName.CompareTo("LHC17f8j") == 0 || periodName.CompareTo("LHC17f8k") == 0){
    fPeriodEnum = kLHC16P1JJ;
    fEnergyEnum = k13TeV;
  // 13TeV LHC16* anchors low field JJ Pythia 8 MB
  } else if ( periodName.CompareTo("LHC16P1JJLowB") == 0 || periodName.CompareTo("LHC17f8b") == 0  ){
    fPeriodEnum = kLHC16P1JJLowB;
    fEnergyEnum = k13TeVLowB;

  // 13TeV HF-MC anchors LHC16i,j,o,p
  } else if ( periodName.CompareTo("LHC17h8c") == 0){
    fPeriodEnum = kLHC17h8c;
    fEnergyEnum = k13TeV;
  // 13TeV HF-MC anchors LHC16d,e,g,h,j,o,p
  } else if ( periodName.CompareTo("LHC17h8b") == 0){
    fPeriodEnum = kLHC17h8b;
    fEnergyEnum = k13TeV;
  // 13TeV HF-MC anchors LHC16d,e,g,h,j,o,p
  } else if ( periodName.CompareTo("LHC17h8a") == 0){
    fPeriodEnum = kLHC17h8a;
    fEnergyEnum = k13TeV;

  // 13TeV HF-MC anchors LHC16k
  } else if ( periodName.CompareTo("LHC17c3b1") == 0){
    fPeriodEnum = kLHC17c3b1;
    fEnergyEnum = k13TeV;
  // 13TeV HF-MC anchors LHC16k
  } else if ( periodName.CompareTo("LHC17c3a1") == 0){
    fPeriodEnum = kLHC17c3a1;
    fEnergyEnum = k13TeV;
  // 13TeV HF-MC anchors LHC16l
  } else if ( periodName.CompareTo("LHC17c3b2") == 0){
    fPeriodEnum = kLHC17c3b2;
    fEnergyEnum = k13TeV;
  // 13TeV HF-MC anchors LHC16l
  } else if ( periodName.CompareTo("LHC17c3a2") == 0){
    fPeriodEnum = kLHC17c3a2;
    fEnergyEnum = k13TeV;

  // 13TeV GJ-MC anchors LHC16i,j,k,l,o,p
  } else if (periodName.CompareTo("LHC17i3a1") == 0){
    fPeriodEnum = kLHC17i3a1;
    fEnergyEnum = k13TeV;

    // 13TeV JJ-MC anchors LHC16i,j,k,l,o,p  (with decay photon in EMCal acc.)
  } else if (periodName.CompareTo("LHC17i3b1") == 0){
      fPeriodEnum = kLHC17i3b1;
      fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17i3c1") == 0){
      fPeriodEnum = kLHC17i3c1;
      fEnergyEnum = k13TeV;

    // 13TeV JJ-MC anchors LHC16i,j,k,l,o,p  (with decay photon in DCal/PHOS acc.)
  } else if (periodName.CompareTo("LHC17i3b2") == 0){
      fPeriodEnum = kLHC17i3b2;
      fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17i3c2") == 0){
      fPeriodEnum = kLHC17i3c2;
      fEnergyEnum = k13TeV;

    // 13TeV JJ-MC anchors LHC16i,j,k,l,o,p new production (with decay photon in EMCal acc.)
  } else if (periodName.CompareTo("LHC20b1b1") == 0){
      fPeriodEnum = kLHC20b1b1;
      fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC20b1c1") == 0){
      fPeriodEnum = kLHC20b1c1;
      fEnergyEnum = k13TeV;

    // 13TeV JJ-MC anchors LHC16i,j,k,l,o,p new production (with decay photon in DCal/PHOS acc.)
  } else if (periodName.CompareTo("LHC20b1b2") == 0){
      fPeriodEnum = kLHC20b1b2;
      fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC20b1c2") == 0){
      fPeriodEnum = kLHC20b1c2;
      fEnergyEnum = k13TeV;

  // LHC16qt anchored MCs
  } else if ( periodName.CompareTo("LHC17f2a") == 0          || periodName.CompareTo("LHC17f2a_fast") == 0 ||
              periodName.CompareTo("LHC17f2a_cent") == 0     || periodName.CompareTo("LHC17f2a_cent_woSDD") == 0 ||
              periodName.CompareTo("LHC17f2a_fast_fix") == 0 || periodName.CompareTo("LHC17f2a_cent_fix") == 0 ||
              periodName.CompareTo("LHC17f2a_cent_woSDD_fix") == 0){
    fPeriodEnum = kLHC17f2a;
    fEnergyEnum = kpPb5TeVR2;
  } else if (periodName.CompareTo("LHC17f2b") == 0          || periodName.CompareTo("LHC17f2b_fast") == 0 ||
              periodName.CompareTo("LHC17f2b_cent") == 0     || periodName.CompareTo("LHC17f2b_cent_woSDD") == 0 ||
              periodName.CompareTo("LHC17f2b_cent_woSDD_fix") == 0){
    fPeriodEnum = kLHC17f2b;
    fEnergyEnum = kpPb5TeVR2;
  } else if (periodName.Contains("LHC18f3b") || periodName.Contains("LHC18f3c") ){
    fPeriodEnum = kLHC18f3bc;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.Contains("LHC18f3") ){
    fPeriodEnum = kLHC18f3;
    fEnergyEnum = kpPb5TeVR2;
  } else if ( periodName.CompareTo("LHC17g8a") == 0 || periodName.CompareTo("LHC17g8a_fast") == 0 ||
              periodName.CompareTo("LHC17g8a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17g8a;
    fEnergyEnum = kpPb5TeVR2;
  } else if (periodName.CompareTo("LHC17d2a") == 0 || periodName.CompareTo("LHC17d2a_fast") == 0 || periodName.CompareTo("LHC17d2a_cent") == 0 ){
    fPeriodEnum = kLHC17d2a;
    fEnergyEnum = kpPb5TeVR2;
  } else if (periodName.CompareTo("LHC17d2b") == 0 || periodName.CompareTo("LHC17d2b_fast") == 0 || periodName.CompareTo("LHC17d2b_cent") == 0 ){
    fPeriodEnum = kLHC17d2b;
    fEnergyEnum = kpPb5TeVR2;
    // LHC16r anchored MCs
  } else if (periodName.CompareTo("LHC17a3a") == 0      || periodName.CompareTo("LHC17a3a_fast") == 0 ||
             periodName.CompareTo("LHC17a3a_cent") == 0 || periodName.CompareTo("LHC17a3a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a3a;
    fEnergyEnum = kpPb8TeV;
  } else if ( periodName.CompareTo("LHC17a3b") == 0      || periodName.CompareTo("LHC17a3b_fast") == 0 ||
              periodName.CompareTo("LHC17a3b_cent") == 0 || periodName.CompareTo("LHC17a3b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a3b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.Contains("LHC17f3")){
    fPeriodEnum = kLHC17f3;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.Contains("LHC17f4")){
    fPeriodEnum = kLHC17f4;
    fEnergyEnum = kpPb8TeV;
  // LHC16s anchored MCs
  } else if (periodName.CompareTo("LHC17a4a") == 0      || periodName.CompareTo("LHC17a4a_fast") == 0 ||
              periodName.CompareTo("LHC17a4a_cent") == 0 || periodName.CompareTo("LHC17a4a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a4a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4b") == 0){
    fPeriodEnum = kLHC17a4b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g6b2a") == 0){
    fPeriodEnum = kLHC17g6b2a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g6b2b") == 0){
    fPeriodEnum = kLHC17g6b2b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g6b3a") == 0){
    fPeriodEnum = kLHC17g6b3a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g6b3b") == 0){
    fPeriodEnum = kLHC17g6b3b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g8b") == 0 || periodName.CompareTo("LHC18b9b") == 0){
    fPeriodEnum = kLHC16rP1JJ;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g8c") == 0 || periodName.CompareTo("LHC18b9c") == 0){
    fPeriodEnum = kLHC16sP1JJ;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17g6b1a") == 0){
    fPeriodEnum = kLHC16rsGJ;
    fEnergyEnum = kpPb8TeV;

  //pp 13 TeV anchored to LHC17
  } else if ( periodName.CompareTo("LHC17P1Pyt8NomB") == 0 ||
              periodName.CompareTo("LHC17k4") ==0 || periodName.CompareTo("LHC17h11") ==0  || periodName.CompareTo("LHC17h1") == 0 ||
              periodName.CompareTo("LHC17l5") == 0 || periodName.CompareTo("LHC18c13") == 0|| periodName.CompareTo("LHC18a8") == 0 ||
              periodName.CompareTo("LHC18a9") == 0 || periodName.CompareTo("LHC18a1") == 0 || periodName.CompareTo("LHC18c12") == 0
            ){
    fPeriodEnum = kLHC17P1Pyt8NomB;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17P1PHONomB13TeV") == 0 ||
              periodName.CompareTo("LHC17h7b") ==0){
    fPeriodEnum = kLHC17P1PHONomB13TeV;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17P1Pyt6NomB") == 0 ||
              periodName.CompareTo("LHC17h7a") ==0 ){
    fPeriodEnum = kLHC17P1Pyt6NomB;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17P1Pyt8LowB") == 0 ||
              periodName.CompareTo("LHC17h3") ==0 ){
    fPeriodEnum = kLHC17P1Pyt8LowB;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17P1JJ") == 0 ||
              periodName.CompareTo("LHC18f5") == 0 ){
    fPeriodEnum = kLHC17P1JJ;
    fEnergyEnum = k13TeV;
  // special JJ MC with decay photon in EMCal acc
  } else if ( periodName.CompareTo("LHC18l6b1") == 0){
    fPeriodEnum = kLHC18l6b1;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC18l6c1") == 0){
    fPeriodEnum = kLHC18l6c1;
    fEnergyEnum = k13TeV;
  // special JJ MC with decay photon in DCal/PHOS acc
  } else if ( periodName.CompareTo("LHC18l6b2") == 0){
    fPeriodEnum = kLHC18l6b2;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC18l6c2") == 0){
    fPeriodEnum = kLHC18l6c2;
    fEnergyEnum = k13TeV;

  // pp 13 TeV LHC17 special MC - strangeness enhanced
  } else if ( periodName.CompareTo("LHC17j5a") ==0 ){
      fPeriodEnum = kLHC17j5a;
      fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17j5b") ==0 ){
      fPeriodEnum = kLHC17j5b;
      fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC17j5c") ==0 ){
      fPeriodEnum = kLHC17j5c;
      fEnergyEnum = k13TeV;

  // MC for Xe-Xe
  } else if (periodName.CompareTo("LHC17j7") == 0){         // HIJING
    fPeriodEnum = kLHC17XeXeHi;
    fEnergyEnum = kXeXe5440GeV;
  } else if (periodName.CompareTo("LHC18d2") == 0 || periodName.CompareTo("LHC18d2_1") == 0 || periodName.CompareTo("LHC18d2_2") == 0 || periodName.CompareTo("LHC18d2_3") == 0){    // HIJING
    fPeriodEnum = kLHC17XeXeHi;
    fEnergyEnum = kXeXe5440GeV;

  // LHC17pq anchored MCs
  } else if ( periodName.CompareTo("LHC17l3b") == 0 || periodName.CompareTo("LHC17l3b_fast") == 0 || periodName.CompareTo("LHC17l3b_cent") == 0 ||
              periodName.CompareTo("LHC17l3b_cent_woSDD") == 0 ||
              periodName.Contains("LHC18d6c")  ){
    fPeriodEnum = kLHC17l3b;
    fEnergyEnum = k5TeV;
  } else if ( periodName.CompareTo("LHC18j2") == 0 || periodName.CompareTo("LHC18j2_fast") == 0 || periodName.CompareTo("LHC18j2_cent") == 0 ||
              periodName.CompareTo("LHC18j2_cent_woSDD") == 0){
    fPeriodEnum = kLHC18j2;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC17l4b") == 0 || periodName.CompareTo("LHC17l4b_fast") == 0 || periodName.CompareTo("LHC17l4b_cent") == 0 ||
            periodName.CompareTo("LHC17l4b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17l4b;
    fEnergyEnum = k5TeV;
  } else if ( periodName.CompareTo("LHC18b8") == 0 || periodName.CompareTo("LHC18b8_fast") == 0 || periodName.CompareTo("LHC18b8_cent") == 0 ||
              periodName.CompareTo("LHC18b8_cent_woSDD") == 0){
    fPeriodEnum = kLHC18b8;
    fEnergyEnum = k5TeV;
  } else if ( periodName.Contains("LHC18b10")){
    fPeriodEnum = kLHC18b10;
    fEnergyEnum = k5TeV;
  } else if ( periodName.Contains("LHC18l2") || periodName.Contains("LHC18g7") ){
    fPeriodEnum = kLHC18l2;
    fEnergyEnum = k5TeV;

  // LHC17p Low Intensity MC using Phojet
  } else if ( periodName.CompareTo("LHC17P1PHO") == 0 ||  periodName.Contains("LHC18d6b")  ){
    fPeriodEnum = kLHC17P1PHO;
    fEnergyEnum = k5TeV;

  //pp 13 TeV anchored to LHC18
  } else if ( periodName.CompareTo("LHC18P1Pyt8NomB") == 0 ||
              periodName.CompareTo("LHC18g4") ==0 || periodName.CompareTo("LHC18g5") ==0  || periodName.CompareTo("LHC18g6") == 0 ||
          periodName.CompareTo("LHC18h2") ==0 || periodName.CompareTo("LHC18h4") ==0  ||
          periodName.CompareTo("LHC18j1") ==0 || periodName.CompareTo("LHC18j4") == 0 ||
          periodName.CompareTo("LHC18k1") ==0 || periodName.CompareTo("LHC18k2") == 0 || periodName.CompareTo("LHC18k3") == 0
  ){
    fPeriodEnum = kLHC18P1Pyt8NomB;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC18P1Pyt8LowB") == 0 ||
              periodName.CompareTo("LHC18h1") ==0 ){
    fPeriodEnum = kLHC18P1Pyt8LowB;
    fEnergyEnum = k13TeV;

  } else if ( periodName.Contains("LHC19d3") || periodName.Contains("LHC18P1JJ") ){
    fPeriodEnum = kLHC18P1JJ;
    fEnergyEnum = k13TeV;

  // special JJ MC with decay photon in EMCal acc
  } else if ( periodName.CompareTo("LHC19i3b1") == 0){
    fPeriodEnum = kLHC19i3b1;
    fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC19i3c1") == 0){
    fPeriodEnum = kLHC19i3c1;
    fEnergyEnum = k13TeV;
  // special JJ MC with decay photon in DCal/PHOS acc
  } else if ( periodName.CompareTo("LHC19i3b2") == 0){
      fPeriodEnum = kLHC19i3b2;
      fEnergyEnum = k13TeV;
  } else if ( periodName.CompareTo("LHC19i3c2") == 0){
    fPeriodEnum = kLHC19i3c2;
    fEnergyEnum = k13TeV;


  // MC upgrade
  } else if (periodName.Contains("LHC13d19")){
    fPeriodEnum = kLHC13d19;
    fEnergyEnum = kPbPb5TeV;
  // fall back
  } else {
    AliWarning("No correct period could be set");
    fPeriodEnum = kUnknownPeriod;
    fEnergyEnum = kUnset;
  }
  return;
}
