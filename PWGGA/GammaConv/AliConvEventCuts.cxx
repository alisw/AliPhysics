/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Authors: Friederike Bock, Daniel Muehlheim                              *
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
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
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
#include "AliEMCALTriggerPatchInfo.h"

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
  fDoLightOutput(kFALSE),
  fEventQuality(-1),
  fIsHeavyIon(0),
  fDetectorCentrality(0),
  fModCentralityClass(0),
  fEnableVertexCut(kTRUE),
  fMaxVertexZ(10),
  fCentralityMin(0),
  fCentralityMax(0),
  fMultiplicityMethod(0),
  fSpecialTrigger(0),
  fSpecialSubTrigger(0),
  fRemovePileUp(kFALSE),
  fPastFutureRejectionLow(0),
  fPastFutureRejectionHigh(0),
  fDoPileUpRejectV0MTPCout(0),
  fFPileUpRejectV0MTPCout(0),
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
  fCutString(NULL),
  fCutStringRead(""),
  fUtils(NULL),
  fEtaShift(0.0),
  fDoEtaShift(kFALSE),
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
  fAddedSignalPDGCode(0),
  fPreSelCut(kFALSE),
  fTriggerSelectedManually(kFALSE),
  fSpecialTriggerName(""),
  fSpecialSubTriggerName(""),
  fNSpecialSubTriggerOptions(0),
  hSPDClusterTrackletBackgroundBefore(NULL),
  hSPDClusterTrackletBackground(NULL),
  fV0ReaderName(""),
  fCaloTriggers(NULL),
  fTriggerPatchInfo(NULL),
  fMainTriggerPatchEMCAL(NULL),
  fCaloTriggersName(""),
  fCaloTriggerPatchInfoName(""),
  fTriggersEMCAL(0),
  fTriggersEMCALSelected(-1),
  fEMCALTrigInitialized(kFALSE),
  fSecProdBoundary(1.0),
  fMaxPtJetMC(0),
  fMaxFacPtHard(2.5),
  fMaxFacPtHardSingleParticle(1.5),
  fMimicTrigger(kFALSE),
  fRejectTriggerOverlap(kFALSE),
  fDoMultiplicityWeighting(kFALSE),
  fPathReweightingMult(""),
  fNameHistoReweightingMultData(""),
  fNameHistoReweightingMultMC(""), 
  hReweightMultData(NULL),
  hReweightMultMC(NULL),
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
  fPastFutureRejectionLow(ref.fPastFutureRejectionLow),
  fPastFutureRejectionHigh(ref.fPastFutureRejectionHigh),
  fDoPileUpRejectV0MTPCout(ref.fDoPileUpRejectV0MTPCout),
  fFPileUpRejectV0MTPCout(ref.fFPileUpRejectV0MTPCout),
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
  fCutString(NULL),
  fCutStringRead(""),
  fUtils(NULL),
  fEtaShift(ref.fEtaShift),
  fDoEtaShift(ref.fDoEtaShift),
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
  fAddedSignalPDGCode(ref.fAddedSignalPDGCode),
  fPreSelCut(ref.fPreSelCut),
  fTriggerSelectedManually(ref.fTriggerSelectedManually),
  fSpecialTriggerName(ref.fSpecialTriggerName),
  fSpecialSubTriggerName(ref.fSpecialSubTriggerName),
  fNSpecialSubTriggerOptions(ref.fNSpecialSubTriggerOptions),
  hSPDClusterTrackletBackgroundBefore(NULL),
  hSPDClusterTrackletBackground(NULL),
  fV0ReaderName(ref.fV0ReaderName),
  fCaloTriggers(NULL),
  fTriggerPatchInfo(NULL),
  fMainTriggerPatchEMCAL(NULL),
  fCaloTriggersName(ref.fCaloTriggersName),
  fCaloTriggerPatchInfoName(ref.fCaloTriggerPatchInfoName),
  fTriggersEMCAL(ref.fTriggersEMCAL),
  fTriggersEMCALSelected(ref.fTriggersEMCALSelected),
  fEMCALTrigInitialized(kFALSE),
  fSecProdBoundary(ref.fSecProdBoundary),
  fMaxPtJetMC(ref.fMaxPtJetMC),
  fMaxFacPtHard(ref.fMaxFacPtHard),
  fMaxFacPtHardSingleParticle(ref.fMaxFacPtHardSingleParticle),
  fMimicTrigger(ref.fMimicTrigger),
  fRejectTriggerOverlap(ref.fRejectTriggerOverlap),
  fDoMultiplicityWeighting(ref.fDoMultiplicityWeighting),
  fPathReweightingMult(ref.fPathReweightingMult),
  fNameHistoReweightingMultData(ref.fNameHistoReweightingMultData),
  fNameHistoReweightingMultMC(ref.fNameHistoReweightingMultMC), 
  hReweightMultData(ref.hReweightMultData),
  hReweightMultMC(ref.hReweightMultMC),
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

}

//________________________________________________________________________
void AliConvEventCuts::InitCutHistograms(TString name, Bool_t preCut){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
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

  if (hReweightMultData){
    hReweightMultData->SetName(Form("hReweightMultData_%s",GetCutNumber().Data()));
    fHistograms->Add(hReweightMultData);
  }
  if (hReweightMultMC){
    hReweightMultMC->SetName(Form("hReweightMultMC_%s",GetCutNumber().Data()));
    fHistograms->Add(hReweightMultMC);
  }
  
  if(!fDoLightOutput){
    hSPDClusterTrackletBackgroundBefore = new TH2F(Form("SPD tracklets vs SPD clusters %s before Pileup Cut",GetCutNumber().Data()),"SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
    fHistograms->Add(hSPDClusterTrackletBackgroundBefore);

    hSPDClusterTrackletBackground = new TH2F(Form("SPD tracklets vs SPD clusters %s",GetCutNumber().Data()),"SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
    fHistograms->Add(hSPDClusterTrackletBackground);
  }

  if(fIsHeavyIon > 0){
    hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",400,0,100);
    fHistograms->Add(hCentrality);
  }
    
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

    hTriggerClass= new TH1F(Form("OfflineTrigger %s",GetCutNumber().Data()),"OfflineTrigger",36,-0.5,35.5);
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
    hTriggerClass->GetXaxis()->SetBinLabel(35,"failed Physics Selection");
    hTriggerClass->GetXaxis()->SetBinLabel(36,"mimickedTrigger");
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
    
    if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9){
      hTriggerClassesCorrelated= new TH1F(Form("TriggerCorrelations %s",GetCutNumber().Data()),"Triggers Correlated with EMCal triggers",10,-0.5,9.5);
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
    if (!mcHandler->TreeTR() ) {
      fEventQuality = 2;
      return kFALSE;
    }
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
    if(fUtils->IsFirstEventInChunk(event)){
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      fEventQuality = 6;
      return kFALSE;
    }
    if(fRemovePileUp){
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
  } else if(fRemovePileUp){
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
    if (hReweightMultDatatemp == NULL) AliError(Form("%s was not contained in %s", fNameHistoReweightingMultData.Data(),fPathReweightingMult.Data() ));
    hReweightMultData = new TH1D(*hReweightMultDatatemp);
    if (hReweightMultData) AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingMultData.Data(),fPathReweightingMult.Data() ));
    else AliWarning(Form("%s not found in %s", fNameHistoReweightingMultData.Data() ,fPathReweightingMult.Data()));
    hReweightMultData->SetDirectory(0);
  }
  if (fNameHistoReweightingMultMC.CompareTo("") != 0 && (fDoMultiplicityWeighting > 0)){
    cout << "I have to find: " <<  fNameHistoReweightingMultMC.Data() << endl;
    TH1D *hReweightMultMCtemp = (TH1D*)w->Get(fNameHistoReweightingMultMC.Data());
    if (hReweightMultMCtemp == NULL) AliError(Form("%s was not contained in %s", fNameHistoReweightingMultMC.Data(),fPathReweightingMult.Data() ));
    hReweightMultMC = new TH1D(*hReweightMultMCtemp);
    if (hReweightMultData) AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingMultMC.Data(),fPathReweightingMult.Data() ));
    else AliWarning(Form("%s not found in %s", fNameHistoReweightingMultMC.Data() ,fPathReweightingMult.Data()));
    hReweightMultMC->SetDirectory(0);
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
    hReweightMCHistPi0 = new TH1D(*hReweightMCHistPi0temp);
    if (hReweightMCHistPi0) AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingPi0.Data(),fPathTrFReweighting.Data() ));
    else AliWarning(Form("%s not found in %s", fNameHistoReweightingPi0.Data() ,fPathTrFReweighting.Data()));
    hReweightMCHistPi0->SetDirectory(0);
  }
  if (fNameFitDataPi0.CompareTo("") != 0 && fDoReweightHistoMCPi0 ){
    cout << "I have to find: " <<  fNameFitDataPi0.Data() << endl;
    TF1 *fFitDataPi0temp = (TF1*)f->Get(fNameFitDataPi0.Data());
    fFitDataPi0 = new TF1(*fFitDataPi0temp);
    if (fFitDataPi0) AliInfo(Form("%s has been loaded from %s", fNameFitDataPi0.Data(),fPathTrFReweighting.Data() ));
    else AliWarning(Form("%s not found in %s",fPathTrFReweighting.Data(), fNameFitDataPi0.Data() ));
  }

  if (fNameHistoReweightingEta.CompareTo("") != 0 && fDoReweightHistoMCEta){
    cout << "I have to find: " <<  fNameHistoReweightingEta.Data() << endl;
    TH1D *hReweightMCHistEtatemp = (TH1D*)f->Get(fNameHistoReweightingEta.Data());
    hReweightMCHistEta = new TH1D(*hReweightMCHistEtatemp);
    if (hReweightMCHistEta) AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingEta.Data(),fPathTrFReweighting.Data() ));
    else AliWarning(Form("%s not found in %s", fNameHistoReweightingEta.Data(),fPathTrFReweighting.Data() ));
    hReweightMCHistEta->SetDirectory(0);
  }

  if (fNameFitDataEta.CompareTo("") != 0 && fDoReweightHistoMCEta){
    cout << "I have to find: " <<  fNameFitDataEta.Data() << endl;
    TF1 *fFitDataEtatemp = (TF1*)f->Get(fNameFitDataEta.Data());
    fFitDataEta = new TF1(*fFitDataEtatemp);
    if (fFitDataEta) AliInfo(Form("%s has been loaded from %s", fNameFitDataEta.Data(),fPathTrFReweighting.Data() ));
    else AliWarning(Form("%s not found in %s", fNameFitDataEta.Data(),fPathTrFReweighting.Data() ));

  }
  if (fNameHistoReweightingK0s.CompareTo("") != 0 && fDoReweightHistoMCK0s){
    cout << "I have to find: " <<  fNameHistoReweightingK0s.Data() << endl;
    TH1D *hReweightMCHistK0stemp = (TH1D*)f->Get(fNameHistoReweightingK0s.Data());
    hReweightMCHistK0s = new TH1D(*hReweightMCHistK0stemp);
    if (hReweightMCHistK0s) AliInfo(Form("%s has been loaded from %s", fNameHistoReweightingK0s.Data(),fPathTrFReweighting.Data() ));
    else AliWarning(Form("%s not found in %s", fNameHistoReweightingK0s.Data(),fPathTrFReweighting.Data() ));
    hReweightMCHistK0s->SetDirectory(0);
  }

  if (fNameFitDataK0s.CompareTo("") != 0 && fDoReweightHistoMCK0s){
    cout << "I have to find: " <<  fNameFitDataK0s.Data() << endl; 
    TF1 *fFitDataK0stemp = (TF1*)f->Get(fNameFitDataK0s.Data());
    fFitDataK0s = new TF1(*fFitDataK0stemp);
    if (fFitDataK0s) AliInfo(Form("%s has been loaded from %s", fNameFitDataK0s.Data(),fPathTrFReweighting.Data() ));
    else AliWarning(Form("%s not found in %s", fNameFitDataK0s.Data(),fPathTrFReweighting.Data() ));
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
  ///Set individual cut ID
//   "HeavyIon",                     //0
//   "CentralityMin",                //1
//   "CentralityMax",                //2
//   "SelectSpecialTrigger",         //3
//   "SelectSpecialSubTriggerClass", //4
//   "RemovePileUp",                 //5
//   "RejectExtraSignals",           //6
//   "VertexCut",                    //7

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
  if (fIsHeavyIon == 0) {
    printf("Running in pp mode \n");
    if (fSpecialTrigger == 0){
      if (fSpecialSubTrigger == 0){
        printf("\t only events triggered by V0OR will be analysed \n");
      } else if (fSpecialSubTrigger == 1){
        printf("\t only events where SDD was present will be analysed \n");
      }
    } else if (fSpecialTrigger == 1){
      if (fSpecialSubTrigger == 0){
        printf("\t only events triggered by V0AND will be analysed \n");
      } else if(fSpecialSubTrigger == 1){
        printf("\t only events where SDD was present will be analysed and triggered by VOAND\n");
      }
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps");
    } else if (fSpecialTrigger > 1){ 
      printf("\t only events triggered by %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data());
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
    } 
    if (fModCentralityClass == 0){
      printf("\t %d - %d \n", fCentralityMin*10, fCentralityMax*10);
    } else if ( fModCentralityClass == 1){ 
      printf("\t %d - %d \n", fCentralityMin*5, fCentralityMax*5);
    } else if ( fModCentralityClass == 2){ 
      printf("\t %d - %d \n", fCentralityMin*5+45, fCentralityMax*5+45);
    } else if (fModCentralityClass == 3){
      printf("\t %d - %d, with Track mult in MC as data \n", fCentralityMin*10, fCentralityMax*10);
    } else if ( fModCentralityClass == 4){ 
      printf("\t %d - %d, with Track mult in MC as data \n", fCentralityMin*5, fCentralityMax*5);
    } else if ( fModCentralityClass == 5){ 
      printf("\t %d - %d, with Track mult in MC as data \n", fCentralityMin*5+45, fCentralityMax*5+45);
    }
    if (fSpecialTrigger == 0){
      printf("\t only events triggered by kMB, kCentral, kSemiCentral will be analysed \n");
    } else if (fSpecialTrigger > 1){
      printf("\t only events triggered by %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data());
      printf("\n\t        SpecialTrigger is:  %s\n", fSpecialTriggerName.Data());
      printf("\t        SpecialSubTrigger is: %s\n", fSpecialSubTriggerName.Data());
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps\n\n");
    }
  } else if (fIsHeavyIon == 2){
    printf("Running in pPb mode \n");
    if (fDetectorCentrality == 0){
      printf("\t centrality selection based on V0A \n");
    } else if (fDetectorCentrality == 1){
      printf("\t centrality selection based on Cl1 \n");
    } 
    if (fModCentralityClass == 0){
      printf("\t %d - %d \n", fCentralityMin*10, fCentralityMax*10);
    }
    if (fSpecialTrigger == 0){
      printf("\t only events triggered by kINT7 will be analysed \n");
    } else if (fSpecialTrigger > 1){ 
      printf("\t only events triggered by %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data());
      if (fRejectTriggerOverlap) printf("\t        reject trigger overlaps\n\n");
    }
  }
  if (fEnableVertexCut) printf("\t Vertex cut with |Z_{vtx}| <%2.2f \n",fMaxVertexZ);
    else printf("\t No vertex cut \n");

  if (fRemovePileUp ==1 ) {
     printf("\t Doing pile up removal  \n");
     if (fDoPileUpRejectV0MTPCout ==1 ){
       printf("\t Doing extra pile up removal V0M vs TPCout  \n");
     }
  }

  printf("MC event cuts: \n");
  if (fRejectExtraSignals == 0) printf("\t no rejection was applied \n");
    else if (fRejectExtraSignals == 1) printf("\t only MB header will be inspected \n");
    else if (fRejectExtraSignals > 1) printf("\t special header have been selected \n");
  printf("\t maximum factor between jet and pt hard = %2.2f \n", fMaxFacPtHard);
}

///________________________________________________________________________
Bool_t AliConvEventCuts::SetIsHeavyIon(Int_t isHeavyIon)
{   // Set Cut
  switch(isHeavyIon){
  case 0:
    fIsHeavyIon=0;
    break;
  case 1:
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    break;
  case 2:
    fIsHeavyIon=1;
    fDetectorCentrality=1;
    break;
  case 3: //allows to select centrality 0-45% in steps of 5% for V0 Multiplicity
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=1;
    break;
  case 4: //allows to select centrality 45-90% in steps of 5% for V0 Multiplicity
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=2;
    break;
  case 5: //strict cut on v0 tracks for MC
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=3;
    break;
  case 6: //allows to select centrality 0-45% in steps of 5% for track mult
    //strict cut on v0 tracks for MC
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=4;
    break;
  case 7: //allows to select centrality 45-90% in steps of 5% for V0 Multiplicity
    //strict cut on v0 tracks for MC
    fIsHeavyIon=1;
    fDetectorCentrality=0;
    fModCentralityClass=5;
    break;
  case 8:
    fIsHeavyIon=2;
    fDetectorCentrality=0;
    break;
  case 9:
    fIsHeavyIon=2;
    fDetectorCentrality=1;
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
  if(minCentrality<0||minCentrality>9){
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
  if(maxCentrality<0||maxCentrality>9){
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
//   case 2:
//     fSpecialTrigger=2; // 
//     break;
  case 3:       
    fSpecialTrigger=3; //specific centrality trigger selection
    fSpecialTriggerName="AliVEvent::kCentral/kSemiCentral/kMB";
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
  default:
    AliError("Warning: Special Trigger Not known");
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
//       AliInfo("Info: Nothing to be done");
      break;
    case 3: //V0OR with SDD requested (will only work with LHC11a dataset)
      fSpecialSubTrigger=1; 
//       cout << "V0OR with SDD requested" << endl;
      break;
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    }
  } else if (fSpecialTrigger == 1){ //AND with different detectors
    switch(selectSpecialSubTriggerClass){
    case 0:  //with VZERO general implementation of V0AND (periods LHC11c onwards)
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
//       AliInfo("Info: Nothing to be done");
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
//       AliInfo("Info: Nothing to be done");
      break;
    case 1: // kCentral - no vertex restriction
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHN";
//       cout << "kCentralOpen" << endl;
      break;
    case 2: // kCentral - T00 +- 10 cm
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CCENT";
//       cout << "kCentralVertex" << endl;
      break;
    case 3: // kCentral - both
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVHN|CCENT|CSEMI|CVLN";
//       cout << "kCentral both" << endl;
      break;
    case 4: // kSemiCentral - no vertex restriction
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CVLN";
//       cout << "kSemiCentralOpen" << endl;
      break;
    case 5: // kSemiCentral - T00 +- 10 cm
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSEMI";
//       cout << "kSemiCentralVertex" << endl;
      break;
    case 6: // kSemiCentral - both
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CSEMI%CVLN";
//       cout << "kSemiCentral both" << endl;
      break;
    case 7: // kMB
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPBI1_|CPBI1-";
//       cout << "kMB 1" << endl;
      break;
    case 8: // kMB
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPBI2_|CPBI2-";
//       cout << "kMB 2" << endl;
      break;
    case 9: // kMB
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPBI2_@CPBI2-@CPBI2_@CPBI2-";
//       cout << "kMB both" << endl;
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
//       AliInfo("Info: Nothing to be done");
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
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    } 
  } else if (fSpecialTrigger == 5){ // Subdivision of kEMC trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0; 
      fSpecialSubTriggerName="";
//       AliInfo("Info: Nothing to be done");
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
    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    } 
  }else if (fSpecialTrigger == 6){ // Subdivision of kPHI trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0; 
      fSpecialSubTriggerName="";
//       AliInfo("Info: Nothing to be done");
      break;
    case 1: // CEMC1 - V0OR and EMCAL fired
      fOfflineTriggerMask=AliVEvent::kPHI1;
      fSpecialTriggerName="AliVEvent::kPHI1";
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPHI1";
      break;
    case 2: // CEMC7 - V0AND and EMCAL fired 
      fSpecialSubTrigger=1; 
      fOfflineTriggerMask=AliVEvent::kPHI7;
      fSpecialTriggerName="AliVEvent::kPHI7";
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="CPHI7";
      break;
    case 3: // CEMC8  - T0OR and EMCAL fired
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
  } else if (fSpecialTrigger == 7){ // Subdivision of kHighMult trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0; 
      fSpecialSubTriggerName="";
//       AliInfo("Info: Nothing to be done");
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

    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    } 
  }else if (fSpecialTrigger == 8){ // Subdivision of kEMCEGA trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0; 
      fSpecialSubTriggerName="";
//       AliInfo("Info: Nothing to be done");
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
    case 10: // 8DG1 - CINT8 DG1
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DG1";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG1);
      break;
    case 11: // 7DG2 - CINT7 DG2
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="7DG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;
    case 12: // 8DG2 - CINT8 DG2
      fSpecialSubTrigger=1; 
      fNSpecialSubTriggerOptions=1;
      fSpecialSubTriggerName="8DG2";
      fTriggersEMCALSelected= 0;
      SETBIT(fTriggersEMCALSelected, kG2);
      break;

    default:
      AliError("Warning: Special Subtrigger Class Not known");
      return 0;
    } 
  } else if (fSpecialTrigger == 9){ // Subdivision of kEMCEGA trigger classes
    switch(selectSpecialSubTriggerClass){
    case 0: // all together
      fSpecialSubTrigger=0; 
      fSpecialSubTriggerName="";
//       AliInfo("Info: Nothing to be done");
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

    default:
      AliError("Warning: Special Subtrigger Class Not known");
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
    break;
  case 2:
    fRemovePileUp           = kTRUE;
    fPastFutureRejectionLow =-89;
    fPastFutureRejectionHigh= 89;
    break;
  case 3:
    fRemovePileUp           = kTRUE;
    fPastFutureRejectionLow = -4;
    fPastFutureRejectionHigh=  7;
    break;
  case 4:
    fRemovePileUp           = kTRUE;
    fPastFutureRejectionLow = -10;
    fPastFutureRejectionHigh=  13;
    break;
  case 5:
    fRemovePileUp           = kTRUE;
    fPastFutureRejectionLow = -40;
    fPastFutureRejectionHigh=  43;
    break;
  case 6:
    fRemovePileUp           = kTRUE;
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
  default:
    AliError("RemovePileUpCut not defined");
    return kFALSE;
  }
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
  default:
    AliError(Form("Vertex Cut not defined %d",vertexCut));
    return kFALSE;
  }
  return kTRUE;
}

//-------------------------------------------------------------
Bool_t AliConvEventCuts::GetUseNewMultiplicityFramework(){ 
  if (fPeriodEnum == kLHC15o ||                                                                                            // PbPb 5TeV
      fPeriodEnum == kLHC15k1a1 || fPeriodEnum == kLHC15k1a2 || fPeriodEnum == kLHC15k1a3  || fPeriodEnum == kLHC16j7 ||   // MC PbPb 5TeV LowIR
      fPeriodEnum == kLHC16h4 ||                                                                                           // MC PbPb 5TeV added signals
      fPeriodEnum == kLHC16g1 || fPeriodEnum == kLHC16g1a || fPeriodEnum == kLHC16g1b || fPeriodEnum == kLHC16g1c ||       // MC PbPb 5TeV general purpose
      fPeriodEnum == kLHC16g2 || fPeriodEnum == kLHC16g3 ||
      fPeriodEnum == kLHC16h2a || fPeriodEnum ==  kLHC16h2b || fPeriodEnum ==  kLHC16h2c ||                                // MC PbPb 5TeV jet-jet
      fPeriodEnum == kLHC15fm ||                                                                                           // pp 13TeV
      fPeriodEnum == kLHC15g3a3 || fPeriodEnum == kLHC15g3c3 ||                                                            // MC pp 13TeV
      fPeriodEnum == kLHC16q || fPeriodEnum == kLHC16t ||                                                                  // pPb 5TeV LHC16qt
      fPeriodEnum == kLHC17a2a || fPeriodEnum == kLHC17a2a_fast || fPeriodEnum == kLHC17a2a_cent || fPeriodEnum == kLHC17a2a_cent_woSDD || // MC pPb 5TeV LHC16qt
      fPeriodEnum == kLHC17a2b || fPeriodEnum == kLHC17a2b_fast || fPeriodEnum == kLHC17a2b_cent || fPeriodEnum == kLHC17a2b_cent_woSDD    // MC pPb 5TeV LHC16qt
      ){
      return kTRUE;
  } else {
     return kFALSE;
  } 
}

//-------------------------------------------------------------
Float_t AliConvEventCuts::GetCentrality(AliVEvent *event)
{   // Get Event Centrality

  AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
  if(esdEvent){
    if(GetUseNewMultiplicityFramework()){
      AliMultSelection *MultSelection = (AliMultSelection*)event->FindListObject("MultSelection");
      if(fDetectorCentrality==0){
	                               return MultSelection->GetMultiplicityPercentile("V0M",kTRUE); // default for pPb
      }else if(fDetectorCentrality==1) return MultSelection->GetMultiplicityPercentile("CL1",kTRUE);
    }else{
      AliCentrality *fESDCentrality = (AliCentrality*)esdEvent->GetCentrality();
      if(fDetectorCentrality==0){
        if(fIsHeavyIon==2)             return fESDCentrality->GetCentralityPercentile("V0A"); // default for pPb
        else                           return fESDCentrality->GetCentralityPercentile("V0M"); // default
      }else if(fDetectorCentrality==1) return fESDCentrality->GetCentralityPercentile("CL1");
    }
  }

  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);
  if(aodEvent){
    if(GetUseNewMultiplicityFramework()){
      AliMultSelection *MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
      if(fDetectorCentrality==0) return MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
      else if(fDetectorCentrality==1) return MultSelection->GetMultiplicityPercentile("CL1",kTRUE);
    }else{
      if(aodEvent->GetHeader()){return ((AliVAODHeader*)aodEvent->GetHeader())->GetCentrality();}
    }
  }

  return -1;
}

//_____________________________________________________________________________________
Bool_t AliConvEventCuts::IsCentralitySelected(AliVEvent *event, AliMCEvent *mcEvent)
{   // Centrality Selection
  if(!fIsHeavyIon){
    if ((fCentralityMin == 0 && fCentralityMax == 0) || (fCentralityMin > fCentralityMax) ){
      return kTRUE;
    } else {
      Int_t primaryTracksPP[9] = { 0,   2,   5,    10,   15, 
                                  30,  50,  100,  1000 
                                  };
      Int_t nprimaryTracks = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetNumberOfPrimaryTracks();
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
  if(centrality<0)return kFALSE;

  Int_t centralityC=0;
  if (fModCentralityClass == 0){
    centralityC= Int_t(centrality/10);
    if(centralityC >= fCentralityMin && centralityC < fCentralityMax)
      return kTRUE;
    else return kFALSE;
  }
  else if (fModCentralityClass ==1){
    centralityC= Int_t(centrality);
    if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
      return kTRUE;
    } else return kFALSE;
  }
  else if (fModCentralityClass ==2){
    centralityC= Int_t(centrality);
    if(centralityC >= ((fCentralityMin*5)+45) && centralityC < ((fCentralityMax*5)+45))
      return kTRUE;
    else return kFALSE;
  }

  Int_t nprimaryTracks = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetNumberOfPrimaryTracks();
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
  Int_t PrimaryTracks5a[11][2] =
    {
      {9999,9999}, // 0 ///1550 changed to 9999 on 9 Dec 
      {1485,1168}, // 5
      {1210, 928}, // 10
      { 995, 795}, // 15
      { 817, 658}, // 20
      { 666, 538}, // 25
      { 536, 435}, // 30
      { 428, 350}, // 35
      { 337, 276}, // 40
      { 260, 214},  // 45
      { 0, 162}// 50 only max accessible
    };
  Int_t PrimaryTracksLHC11h5a[11][2] =
    {
      {9999,9999}, // 0 ///1550 changed to 9999 on 9 Dec
      {1166,1168}, // 5
      { 953, 928}, // 10
      { 805, 795}, // 15
      { 655, 658}, // 20
      { 535, 538}, // 25
      { 435, 435}, // 30
      { 349, 350}, // 35
      { 275, 276}, // 40
      { 214, 214},  // 45
      { 165, 162}// 50 only max accessible
    };
  Int_t PrimaryTracks5b[11][2] =
    {
      { 260, 214}, // 45
      { 197, 162}, // 50
      { 147, 125}, // 55
      { 106, 100}, // 60
      {  75,  63}, // 65
      {  51,  44}, // 70
      {  34,  29}, // 75
      {  21,  18}, // 80
      {  13,  11}, // 85
      {   0,   0},  // 90
      {   0,   0}// 100 only max accessible
    };
  Int_t PrimaryTracksLHC11h5b[11][2] =
    {
      { 214, 214}, // 45
      { 165, 162}, // 50
      { 127, 125}, // 55
      {  93, 100}, // 60
      {  64,  63}, // 65
      {  44,  44}, // 70
      {  30,  29}, // 75
      {  18,  18}, // 80
      {  11,  11}, // 85
      {   0,   0},  // 90
      {   0,   0}// 100 only max accessible
    };
  Int_t column = 0;
  if(event->IsA()==AliESDEvent::Class()) column = 0;
  if(event->IsA()==AliAODEvent::Class()) column = 1;

  if (fModCentralityClass == 3){
    if(mcEvent){
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h10[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h10[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      } else {
        if(nprimaryTracks > PrimaryTracks10[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks10[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      }
    }
    else{
      centralityC= Int_t(centrality/10);
      if(centralityC >= fCentralityMin && centralityC < fCentralityMax)
        return kTRUE;
      else return kFALSE;
    }
  }
  else if (fModCentralityClass ==4){
    if(mcEvent){
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h5a[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h5a[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      } else {
        if(nprimaryTracks > PrimaryTracks5a[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks5a[fCentralityMin][column])
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
  else if (fModCentralityClass ==5){
    if(mcEvent){
      if(fPeriodEnum == kLHC14a1a || fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c){
        if(nprimaryTracks > PrimaryTracksLHC11h5b[fCentralityMax][column] && nprimaryTracks <= PrimaryTracksLHC11h5b[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      } else {
        if(nprimaryTracks > PrimaryTracks5b[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks5b[fCentralityMin][column])
          return kTRUE;
        else return kFALSE;
      }
    }
    else{
      centralityC= Int_t(centrality);
      if(centralityC >= ((fCentralityMin*5)+45) && centralityC < ((fCentralityMax*5)+45))
        return kTRUE;
      else return kFALSE;
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

Bool_t AliConvEventCuts::IsPileUpV0MTPCout(AliVEvent *event)
{
  Bool_t isPileUpV0MTPCout=0;

  Double_t multV0M;
  Double_t valFunc;  
  if (fIsHeavyIon==2){
      multV0M =  event->GetVZEROData()->GetMTotV0A();
  }else{
      multV0M = event->GetVZEROData()->GetMTotV0A() + event->GetVZEROData()->GetMTotV0C() ;
  }
 
  if ( fFPileUpRejectV0MTPCout != 0x0 ){
  valFunc= fFPileUpRejectV0MTPCout->Eval(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetNumberOfTPCoutTracks());
    if (multV0M < valFunc  ) isPileUpV0MTPCout=1;
  }

  return isPileUpV0MTPCout;

}
//________________________________________________________________________
Int_t AliConvEventCuts::GetNumberOfContributorsVtx(AliVEvent *event){
  // returns number of contributors to the vertex

  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(event);
  if(fESDEvent){
    if (fESDEvent->GetPrimaryVertex() != NULL){
      if(fESDEvent->GetPrimaryVertex()->GetNContributors()>0) {
  //     cout << "accepted global" << fESDEvent->GetEventNumberInFile() << " with NCont: " << fESDEvent->GetPrimaryVertex()->GetNContributors() << endl;
        return fESDEvent->GetPrimaryVertex()->GetNContributors();
      }
    }

    if(fESDEvent->GetPrimaryVertexSPD() !=NULL){
      if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
  //     cout << "accepted SPD" << fESDEvent->GetEventNumberInFile() << " with NCont: " << fESDEvent->GetPrimaryVertexSPD()->GetNContributors() << endl;
        return fESDEvent->GetPrimaryVertexSPD()->GetNContributors();
      }else {
        AliWarning(Form("Number of contributors from bad vertex type:: %s",fESDEvent->GetPrimaryVertex()->GetName()));
  //            cout << "rejected " << fESDEvent->GetEventNumberInFile() << endl;
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
Bool_t AliConvEventCuts::IsJetJetMCEventAccepted(AliMCEvent *mcEvent, Double_t& weight){
  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;
  weight                               = -1;
  fMaxPtJetMC                          = 0;
  
  if (  fPeriodEnum != kLHC17f8a &&  fPeriodEnum != kLHC17f8b  && fPeriodEnum != kLHC17f8c &&       // LHC16X Jet Jet MC's
        fPeriodEnum != kLHC17f8d &&  fPeriodEnum != kLHC17f8e &&
        fPeriodEnum != kLHC16h3  &&                                                                 // LHC15n Jet Jet MC's
        fPeriodEnum != kLHC15a3a && fPeriodEnum != kLHC15a3a_plus && fPeriodEnum != kLHC15a3b &&    // LHC13g Jet Jet MC's
        fPeriodEnum != kLHC15g1a && fPeriodEnum != kLHC15g1b &&                                     // LHC11a Jet Jet MC's
        fPeriodEnum != kLHC13b4_fix && fPeriodEnum != kLHC13b4_plus &&                              // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC16c3a && fPeriodEnum != kLHC16c3b && fPeriodEnum != kLHC16c3c &&         // LHC13 pPb Jet Jet MC's        
        fPeriodEnum != kLHC16c2 && fPeriodEnum != kLHC16c2_plus                                     // LHC12 JetJet MC
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
      if (GeneratorName.CompareTo("AliGenPythiaEventHeader") == 0){
        Bool_t eventAccepted = kTRUE;
        TParticle * jet = 0;
        Int_t nTriggerJets = dynamic_cast<AliGenPythiaEventHeader*>(gh)->NTriggerJets();
        Float_t ptHard = dynamic_cast<AliGenPythiaEventHeader*>(gh)->GetPtHard();
        Float_t tmpjet[]={0,0,0,0};
        for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
          dynamic_cast<AliGenPythiaEventHeader*>(gh)->TriggerJet(ijet, tmpjet);
          jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
          //Compare jet pT and pt Hard
          if(jet->Pt() > fMaxFacPtHard * ptHard){
            eventAccepted= kFALSE;
          }
          if (jet->Pt() > fMaxPtJetMC) fMaxPtJetMC = jet->Pt(); 
        }
        if (jet) delete jet;
        if (mcEvent){
          for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
            TParticle* particle = (TParticle *)mcEvent->Particle(i);
            if (!particle) continue;
            if (TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221){
              if (particle->Pt() > fMaxFacPtHardSingleParticle*ptHard){
                eventAccepted= kFALSE;
              }
            }

          }
        }
        
        if ( fPeriodEnum == kLHC16h3 ){
          Double_t ptHardBinRanges[21] = {  5,  7,  9, 12, 16,
                                           21, 28, 36, 45, 57,
                                           70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
          Double_t weightsBins[20]     = {  0.957702, 0.41837, 0.406279, 0.266936, 0.135179,
                                         6.4687e-02, 2.27254e-02, 8.30769e-03, 3.56008e-03, 1.22934e-03,
                                         4.91352e-04, 1.77601e-04, 8.79608e-05, 4.13652e-05, 2.02997e-05,
                                         1.03682e-06, 5.64732e-06, 2.96158e-06, 1.5999e-06, 2.08374e-06};

          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];


        } else if ( fPeriodEnum == kLHC15a3b || fPeriodEnum == kLHC15g1b ){
          Double_t ptHardBinRanges[13] = {  5,  7,  9, 12, 16,
                                           21, 28, 36, 45, 57,
                                           70, 85, 1000};
          Double_t weightsBins[12]     = {  7.858393e-03, 4.718691e-03, 4.077575e-03, 2.814527e-03, 1.669625e-03,
                                            1.007535e-03, 4.536554e-04, 2.111041e-04, 1.094840e-04, 4.404973e-05,
                                            1.933238e-05, 1.562895e-05};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 12) weight = weightsBins[bin];

        } else if ( fPeriodEnum == kLHC15g1a ){
          Double_t ptHardBinRanges[20]    = { 5,      11,     21,     36,      57,
                                             84,     117,    152,    191, 1000000,
                                              5,       7,      9,     12,      16,
                                             21,      28,     36,     45,      57 };
    //                     Double_t weightsBins[19]        = { 4.407782 , 4.946649e-01, 3.890474e-02, 3.826300e-03, 4.429376e-04,
    //                                                         6.306745e-05, 1.031527e-05, 2.267429e-06, 7.552074e-07, 0,
    //                                                         2.4635e+00, 1.1483e+00, 6.5069e-01, 2.7130e-01,  8.1947e-02, 
    //                                                         3.1536e-02, 9.3139e-03, 2.9779e-03, 1.1252e-03};
                    // LHC15g1a                                    
          Double_t weightsBins[19]        = { 4.43629 , 0.49523, 0.0394921, 0.00383174, 0.000446559,
                                              6.37374e-05, 1.03134e-05, 2.27012e-06, 7.59281e-07, 0,
                                              2.62906, 1.12884, 0.656873, 0.262822,  0.0876732,
                                              0.0307759, 0.0087083, 0.0027664, 0.00106203};

          Int_t bin = 0;
          Int_t binFromFile = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPtHardFromFile();
          if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 19) weight = weightsBins[bin];

        } else if ( fPeriodEnum == kLHC15a3a || fPeriodEnum == kLHC15a3a_plus ) {
          Double_t ptHardBinRanges[20]    = { 5,      11,     21,     36,      57,
                                             84,     117,    152,    191, 1000000,
                                              5,       7,      9,     12,      16,
                                             21,      28,     36,     45,      57 };
          // LHC15a3a
          Double_t weightsBins[19]        = { 4.43897 , 0.495766, 0.039486, 0.00383011, 0.000447104,
                                           6.37277e-05, 1.03166e-05, 2.26971e-06, 7.59023e-07, 0,
                                               2.63331, 1.12815, 0.657034, 0.262756,  0.0877227,
                                             0.0307638, 0.00870635, 0.00276658, 0.00106229};
          Int_t bin = 0;
          Int_t binFromFile = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPtHardFromFile();
          if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 19) weight = weightsBins[bin];

        } else if ( fPeriodEnum == kLHC16c2 || fPeriodEnum == kLHC16c2_plus ){
            Double_t ptHardBinRanges[21] = {  5,  7,  9, 12, 16,
                                             21, 28, 36, 45, 57,
                                             70, 85, 99, 115, 132,
                                             150, 169, 190, 212, 235,
                                             1000000};
            Double_t weightsBins[20]     = {  28.3084, 8.43277, 4.07753, 1.54359, 0.543318,
                                              0.208394, 0.0652349, 0.0186904, 0.00834528, 0.00301414,
                                              0.00125939, 0.000474403, 0.000244052, 0.00011924, 6.09838e-05,
                                              3.24148e-05, 1.84314e-05, 1.00926e-05, 5.68632e-06, 8.38092e-06};
            Int_t bin = 0;
            while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
            if (bin < 20) weight = weightsBins[bin];

        } else if ( fPeriodEnum == kLHC16c3a ){
            Double_t ptHardBinRanges[6] = {  7, 9, 12, 16, 21, 1000};
            Double_t weightsBins[5]     = {  0.00672445, 0.00799158, 0.00678934, 0.00463908, 0.00600068}; //preliminary estimates
            Int_t bin = 0;
            while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
            if (bin < 5) weight = weightsBins[bin];

        } else if (fPeriodEnum == kLHC16c3b ){
            Double_t ptHardBinRanges[7] = {  14, 19, 26, 35, 48, 66, 1000};
            Double_t weightsBins[6]     = {  0.00608281, 0.00393646, 0.00200138, 0.000986267, 0.000389051, 0.0001863}; //preliminary estimates
            Int_t bin = 0;
            while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
            if (bin < 6) weight = weightsBins[bin];

        } else if (fPeriodEnum == kLHC16c3c ){
            Double_t ptHardBinRanges[8] = {  0, 5, 11, 21, 36, 57, 84, 1000};
            Double_t weightsBins[7]     = {  0.00151999, 0.000100346, 1.27688e-05, 1.82388e-06, 3.08506e-07, 6.00308e-08, 1.88414e-08}; //preliminary estimates
            Int_t bin = 0;
            while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
            if (bin < 7) weight = weightsBins[bin];

        } else if ( fPeriodEnum == kLHC13b4_fix || fPeriodEnum == kLHC13b4_plus ){
          Double_t ptHardBinRanges[11]   = {  5,     11,   21,   36,   57, 
                                             84,    117,   152,  191,    234,
                                             1000};
          Double_t weightsBins[10]     = {  2.24185e-6 , 2.48463e-7, 2.23171e-8, 2.43667e-9, 3.29934e-10,
                                            5.34592e-11, 1.00937e-11, 2.6493e-12, 8.53912e-13, 5.43077e-13};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 10) weight = weightsBins[bin];

        } else {
          weight = 1;
        }

        if (weight == -1) return kFALSE;
        else return eventAccepted;

      }
    }
  } else {    
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    TString eventHeaderName     = eventHeader->ClassName();
    if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0){
      Bool_t eventAccepted = kTRUE;
      TParticle * jet =  0;
      Int_t nTriggerJets =  dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->NTriggerJets();
      Float_t ptHard = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->GetPtHard();
      Float_t tmpjet[]={0,0,0,0};
      for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
        dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->TriggerJet(ijet, tmpjet);
        jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
        //Compare jet pT and pt Hard
        if(jet->Pt() > fMaxFacPtHard * ptHard){
          eventAccepted= kFALSE;
        }
        if (jet->Pt() > fMaxPtJetMC) fMaxPtJetMC = jet->Pt(); 
      }
      if (mcEvent){
        for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
          TParticle* particle = (TParticle *)mcEvent->Particle(i);
          if (!particle) continue;
          if (TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221){
            if (particle->Pt() > fMaxFacPtHardSingleParticle*ptHard){
              eventAccepted= kFALSE;
            }
          }
          
        }
      }
      
      if ( fPeriodEnum == kLHC16h3 ){
          Double_t ptHardBinRanges[21] = {  5,  7,  9, 12, 16,
                                           21, 28, 36, 45, 57,
                                           70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
          Double_t weightsBins[20]     = {  0.957702, 0.41837, 0.406279, 0.266936, 0.135179,
                                         6.4687e-02, 2.27254e-02, 8.30769e-03, 3.56008e-03, 1.22934e-03,
                                         4.91352e-04, 1.77601e-04, 8.79608e-05, 4.13652e-05, 2.02997e-05,
                                         1.03682e-06, 5.64732e-06, 2.96158e-06, 1.5999e-06, 2.08374e-06};

          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC15a3b || fPeriodEnum == kLHC15g1b ){
        Double_t ptHardBinRanges[13]   = {  5,   7,   9,   12, 16, 
                          21,  28, 36, 45, 57, 
                          70, 85, 1000};
        Double_t weightsBins[12]     = {  7.858393e-03, 4.718691e-03, 4.077575e-03, 2.814527e-03, 1.669625e-03,
                          1.007535e-03, 4.536554e-04, 2.111041e-04, 1.094840e-04, 4.404973e-05,
                          1.933238e-05, 1.562895e-05};
        Int_t bin = 0;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 12) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC15g1a ){
        Double_t ptHardBinRanges[20]    = { 5,      11,     21,     36,     57,
                                           84,     117,    152,    191,    1000000,
                                            5,      7,      9,      12,     16,
                                           21,     28,     36,     45,     57 };
//                     Double_t weightsBins[19]        = { 4.407782 , 4.946649e-01, 3.890474e-02, 3.826300e-03, 4.429376e-04,
//                                                         6.306745e-05, 1.031527e-05, 2.267429e-06, 7.552074e-07, 0,
//                                                         2.4635e+00, 1.1483e+00, 6.5069e-01, 2.7130e-01,  8.1947e-02, 
//                                                         3.1536e-02, 9.3139e-03, 2.9779e-03, 1.1252e-03};
                // LHC15g1a                                    
        Double_t weightsBins[19]        = { 4.43629 , 0.49523, 0.0394921, 0.00383174, 0.000446559,
                                            6.37374e-05, 1.03134e-05, 2.27012e-06, 7.59281e-07, 0,
                                            2.62906, 1.12884, 0.656873, 0.262822,  0.0876732,
                                            0.0307759, 0.0087083, 0.0027664, 0.00106203};

        Int_t bin = 0;
        Int_t binFromFile = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPtHardFromFile();
        if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 19) weight = weightsBins[bin];

      } else if (fPeriodEnum == kLHC15a3a || fPeriodEnum == kLHC15a3a_plus ) {
        Double_t ptHardBinRanges[20]    = { 5,      11,     21,     36,     57,
                                           84,     117,    152,    191,    1000000,
                                            5,      7,      9,      12,     16,
                                           21,     28,     36,     45,     57 };
                // LHC15a3a                                    
        Double_t weightsBins[19]        = { 4.43897 , 0.495766, 0.039486, 0.00383011, 0.000447104,
                                            6.37277e-05, 1.03166e-05, 2.26971e-06, 7.59023e-07, 0,
                                            2.63331, 1.12815, 0.657034, 0.262756,  0.0877227,
                                            0.0307638, 0.00870635, 0.00276658, 0.00106229};
        Int_t bin = 0;
        Int_t binFromFile = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPtHardFromFile();
        if (binFromFile != -1 && binFromFile >9 && ptHard < 57) bin = 9;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 19) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC16c2 || fPeriodEnum == kLHC16c2_plus ){
          Double_t ptHardBinRanges[21] = {  5,  7,  9, 12, 16,
                                           21, 28, 36, 45, 57,
                                           70, 85, 99, 115, 132,
                                           150, 169, 190, 212, 235,
                                           1000000};
          Double_t weightsBins[20]     = {  28.3084, 8.43277, 4.07753, 1.54359, 0.543318,
                                            0.208394, 0.0652349, 0.0186904, 0.00834528, 0.00301414,
                                            0.00125939, 0.000474403, 0.000244052, 0.00011924, 6.09838e-05,
                                            3.24148e-05, 1.84314e-05, 1.00926e-05, 5.68632e-06, 8.38092e-06};
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC16c3a ){
          Double_t ptHardBinRanges[6] = {  7, 9, 12, 16, 21, 1000};
          Double_t weightsBins[5]     = {  0.00672445, 0.00799158, 0.00678934, 0.00463908, 0.00600068}; //preliminary estimates
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 5) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC16c3b ){
          Double_t ptHardBinRanges[7] = {  14, 19, 26, 35, 48, 66, 1000};
          Double_t weightsBins[6]     = {  0.00608281, 0.00393646, 0.00200138, 0.000986267, 0.000389051, 0.0001863}; //preliminary estimates
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 6) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC16c3c ){
          Double_t ptHardBinRanges[8] = {  0, 5, 11, 21, 36, 57, 84, 1000};
          Double_t weightsBins[7]     = {  0.00151999, 0.000100346, 1.27688e-05, 1.82388e-06, 3.08506e-07, 6.00308e-08, 1.88414e-08}; //preliminary estimates
          Int_t bin = 0;
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 7) weight = weightsBins[bin];

      } else if ( fPeriodEnum == kLHC13b4_plus ||  fPeriodEnum == kLHC13b4_fix ){
        Double_t ptHardBinRanges[11]   = {  5,     11,   21,   36,   57, 
                                           84,    117,  152,  191,  234,
                                          1000};
        Double_t weightsBins[10]     = {  2.24185e-6 , 2.48463e-7, 2.23171e-8, 2.43667e-9, 3.29934e-10,
                                          5.34592e-11, 1.00937e-11, 2.6493e-12, 8.53912e-13, 5.43077e-13};
        Int_t bin = 0;
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 10) weight = weightsBins[bin];
      } else {
        weight = 1;
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
void AliConvEventCuts::GetXSectionAndNTrials(AliMCEvent *mcEvent, Float_t &XSection, Float_t &NTrials){

  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;

  if (  fPeriodEnum != kLHC17f8a && fPeriodEnum != kLHC17f8b && fPeriodEnum != kLHC17f8c &&         // LHC16X Jet Jet MC's
        fPeriodEnum != kLHC17f8d && fPeriodEnum != kLHC17f8e &&
        fPeriodEnum != kLHC16h3 &&                                                                  // LHC15n Jet Jet MC's
        fPeriodEnum != kLHC15a3a && fPeriodEnum != kLHC15a3a_plus && fPeriodEnum != kLHC15a3b &&    // LHC13g Jet Jet MC's
        fPeriodEnum != kLHC15g1a && fPeriodEnum != kLHC15g1b &&                                     // LHC11a Jet Jet MC's
        fPeriodEnum != kLHC13b4_fix && fPeriodEnum != kLHC13b4_plus &&                              // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC16c3a && fPeriodEnum != kLHC16c3b && fPeriodEnum != kLHC16c3c &&         // LHC13 pPb Jet Jet MC's        
        fPeriodEnum != kLHC16c2 && fPeriodEnum != kLHC16c2_plus                                     // LHC12 JetJet MC
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
      if (GeneratorName.CompareTo("AliGenPythiaEventHeader") == 0){
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
      if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0){
        AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader);
        NTrials = gPythia->Trials();
        XSection = gPythia->GetXsection();
        return;
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
Float_t AliConvEventCuts::GetPtHard(AliMCEvent *mcEvent){
  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;
  
  if (  fPeriodEnum != kLHC17f8a && fPeriodEnum != kLHC17f8b && fPeriodEnum != kLHC17f8c &&         // LHC16X Jet Jet MC's
        fPeriodEnum != kLHC17f8d && fPeriodEnum != kLHC17f8e &&
        fPeriodEnum != kLHC16h3 &&                                                                  // LHC15n Jet Jet MC's
        fPeriodEnum != kLHC15a3a && fPeriodEnum != kLHC15a3a_plus && fPeriodEnum != kLHC15a3b &&    // LHC13g Jet Jet MC's
        fPeriodEnum != kLHC15g1a && fPeriodEnum != kLHC15g1b &&                                     // LHC11a Jet Jet MC's
        fPeriodEnum != kLHC13b4_fix && fPeriodEnum != kLHC13b4_plus &&                              // LHC13 pPb Jet Jet MC's
        fPeriodEnum != kLHC16c3a && fPeriodEnum != kLHC16c3b && fPeriodEnum != kLHC16c3c &&         // LHC13 pPb Jet Jet MC's        
        fPeriodEnum != kLHC16c2 && fPeriodEnum != kLHC16c2_plus                                     // LHC12 JetJet MC
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
      if (GeneratorName.CompareTo("AliGenPythiaEventHeader") == 0){
        return dynamic_cast<AliGenPythiaEventHeader*>(gh)->GetPtHard();
      } 
    }
  } else {    
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    if(eventHeader){
      TString eventHeaderName     = eventHeader->ClassName();
      if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0){
        return dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->GetPtHard();
      }
    }
  }
  
  return -1;
}


//________________________________________________________________________
Bool_t AliConvEventCuts::MimicTrigger(AliVEvent *event, Bool_t isMC ){
  // abort if mimicing not enabled

  if (!fMimicTrigger) return kTRUE;
  
  Int_t runRangesEMCalL0 [35]   = { 144871, 145288, 146375, 146382,  // LHC11a
                                    146502, 148522,         // LHC11a
                                    150209, 153056, 153911, 153915, // LHC11b,c,d
                                    158135, 158136, 158178, 158182, 160683,
                                    160764, 161139, 161256, 161379, 161457,
                                    161525, 161556, 161558, 161609, 161630,
                                    161724, // LHC11d,e
                                    173731, 177144, 177147, 177653, 177724, 178327,
                                    195180,              // LHC13b-f  
                                    197469, 197692            // LHC13g
  };
  
  Double_t thresholdEMCalL0[34] = { 2.11, 3.43, 1.71, 2.05,   // LHC11a 7 TeV
                                    3.43,           // LHC11a  2.76TeV
                                    1.94, 3.39, 4.01, 5.25, 5.5,     // LHC11b, LHC11c, LHC11d
                                    2.05, 5.50, 2.05, 5.50, 2.05, 1.71, 5.50, 1.71, 5.50, 1.71, 5.50, 1.71, 5.50, 1.71, 5.50, 1.71,
                                    2.01, 1.75, 1.52, 2.01, 1.52, 1.85,
                                    3.2,
                                    /*2.01*/1.8 
  };
  Double_t spreadEMCalL0[34]    = { 0., 0., 0, 0,   // LHC11a 7TeV
                                    /*0.7*/0.65,           // LHC11a 2.76TeV    
                                    0., 0., 0., 0., 0.,     // LHC11b, LHC11c, LHC11d
                                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                    0., 0., 0., 0., 0.2, 0.2,/*0.,0.,*/
                                    0.1,
                                    /*0.1*/0.12 
  };

  Int_t runRangesEMCalL1[4]     = { 179796,             // LHC12c-i
                                    195180,              // LHC13b-f  
                                    197469, 197692            // LHC13g
  };
  
  Double_t thresholdEMCalL1[3]  = { 9.5/*8.398*/, 11.5, /*6.*/5.5};
  Double_t spreadEMCalL1[3]     = { 1.0/*0.*/, 0.5, /*0.4*/0.6};
  
  Int_t runRangesEMCalL1G2[3]   = { 195180,              // LHC13b-f  
                                    197469, 197692            // LHC13g
  };
  
  Double_t thresholdEMCalL1G2[2]  = { 7.2, /*3.9*/3.75};
  Double_t spreadEMCalL1G2[2]     = { 0.3, /*0.2*/0.25};
  
  Int_t runnumber = event->GetRunNumber();
  
  if (fSpecialTrigger == 5 ){
    if (runnumber < runRangesEMCalL0[0]) return kTRUE;
    Int_t binRun = 0;
    while (!(runnumber >= runRangesEMCalL0[binRun] && runnumber < runRangesEMCalL0[binRun+1] ) && binRun < 34 ){
//       cout << runnumber << "\t" << binRun << "\t" << runRangesEMCalL0[binRun] << "\t" << runRangesEMCalL0[binRun+1] << endl;
      binRun++;
    }
    if (binRun==34) return kFALSE;
    Double_t threshold = thresholdEMCalL0[binRun];
    
    if (isMC && spreadEMCalL0[binRun] != 0.){
      TF1* triggerSmearing =  new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)",0,15);
      triggerSmearing->SetParameter(0, 1/(spreadEMCalL0[binRun]*TMath::Sqrt(TMath::Pi()*2)));
      triggerSmearing->SetParameter(1, thresholdEMCalL0[binRun]);
      triggerSmearing->SetParameter(2, spreadEMCalL0[binRun]);
      threshold = triggerSmearing->GetRandom();
      delete triggerSmearing;
    }
    
//     cout << runnumber << "\t"<< binRun << "\t"<< threshold << endl;
    
    Int_t nclus = 0;
    nclus = event->GetNumberOfCaloClusters();
  
    if(nclus == 0)  return kFALSE;
    
    // Loop over EMCal clusters
    Bool_t eventIsAccepted = kFALSE;
    for(Int_t i = 0; i < nclus; i++){  
      AliVCluster* clus = NULL;
      clus = event->GetCaloCluster(i);
      if (!clus) continue;
      if (!clus->IsEMCAL()) continue;
      if (clus->GetM02()<0.1) continue;
      if (clus->GetNCells()<2) continue;
      if (clus->E() > threshold ){
//         cout << "found L0" << endl;
        eventIsAccepted = kTRUE;
      }
    }
    return eventIsAccepted;
    
  } else if (fSpecialTrigger == 6 ) {

    return kTRUE;
  } else if (fSpecialTrigger == 8 ) {
    if (fSpecialSubTriggerName.CompareTo("7EGA")==0 || fSpecialSubTriggerName.CompareTo("8EGA")==0 || fSpecialSubTriggerName.CompareTo("7EG1")==0 ||fSpecialSubTriggerName.CompareTo("8EG1")==0 ){
      if (runnumber < runRangesEMCalL1[0]) return kTRUE;
      Int_t binRun = 0;
      while (!(runnumber >= runRangesEMCalL1[binRun] && runnumber < runRangesEMCalL1[binRun+1] ) && binRun < 3 ){
  //       cout << runnumber << "\t" << binRun << "\t" << runRangesEMCalL0[binRun] << "\t" << runRangesEMCalL0[binRun+1] << endl;
        binRun++;
      }
      if (binRun==3) return kFALSE;
      Double_t threshold = thresholdEMCalL1[binRun];

      if (isMC && spreadEMCalL1[binRun] != 0.){
        TF1* triggerSmearing =  new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)",0,15);
        triggerSmearing->SetParameter(0, 1/(spreadEMCalL1[binRun]*TMath::Sqrt(TMath::Pi()*2)));
        triggerSmearing->SetParameter(1, thresholdEMCalL1[binRun]);
        triggerSmearing->SetParameter(2, spreadEMCalL1[binRun]);
        threshold = triggerSmearing->GetRandom();
        delete triggerSmearing;
      }
      
//       cout << runnumber << "\t"<< binRun << "\t L1 \t"<< threshold << endl;
      
      Int_t nclus = 0;
      nclus = event->GetNumberOfCaloClusters();
    
      if(nclus == 0)  return kFALSE;
      
      // Loop over EMCal clusters
      Bool_t eventIsAccepted = kFALSE;
      for(Int_t i = 0; i < nclus; i++){  
        AliVCluster* clus = NULL;
        clus = event->GetCaloCluster(i);
        if (!clus) continue;
        if (!clus->IsEMCAL()) continue;
        if (clus->GetM02()<0.1) continue;
        if (clus->GetNCells()<2) continue;
        if (clus->E() > threshold ){
//           cout << "found L1G1" << endl;
          eventIsAccepted = kTRUE;
        }
      }
      return eventIsAccepted;
    } else if ( fSpecialSubTriggerName.CompareTo("7EG2")==0 ||fSpecialSubTriggerName.CompareTo("8EG2")==0 ){  
      if (runnumber < runRangesEMCalL1G2[0]) return kTRUE;
      Int_t binRun = 0;
      while (!(runnumber >= runRangesEMCalL1G2[binRun] && runnumber < runRangesEMCalL1G2[binRun+1] ) && binRun < 2 ){
  //       cout << runnumber << "\t" << binRun << "\t" << runRangesEMCalL0[binRun] << "\t" << runRangesEMCalL0[binRun+1] << endl;
        binRun++;
      }
      if (binRun==2) return kFALSE;
      Double_t threshold = thresholdEMCalL1G2[binRun];
      if (isMC && spreadEMCalL1G2[binRun] != 0.){
        TF1* triggerSmearing =  new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)",0,15);
        triggerSmearing->SetParameter(0, 1/(spreadEMCalL1G2[binRun]*TMath::Sqrt(TMath::Pi()*2)));
        triggerSmearing->SetParameter(1, thresholdEMCalL1G2[binRun]);
        triggerSmearing->SetParameter(2, spreadEMCalL1G2[binRun]);
        threshold = triggerSmearing->GetRandom();
        delete triggerSmearing;
      }
//       cout << runnumber << "\t"<< binRun << "\t L2 \t"<< threshold << endl;
      
      Int_t nclus = 0;
      nclus = event->GetNumberOfCaloClusters();
    
      if(nclus == 0)  return kFALSE;
      
      // Loop over EMCal clusters
      Bool_t eventIsAccepted = kFALSE;
      for(Int_t i = 0; i < nclus; i++){  
        AliVCluster* clus = NULL;
        clus = event->GetCaloCluster(i);
        if (!clus) continue;
        if (!clus->IsEMCAL()) continue;
        if (clus->GetM02()<0.1) continue;
        if (clus->GetNCells()<2) continue;
        if (clus->E() > threshold ){
//           cout << "found L1G2" << endl;
          eventIsAccepted = kTRUE;
        }
      }
      return eventIsAccepted;
    }
    return kTRUE;
  } else if (fSpecialTrigger == 9 ) {
    return kTRUE;
  } else {
    return kTRUE;
  } 
  
  return kTRUE;
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
            fOfflineTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;    
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
  
    if (fOfflineTriggerMask){
      isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected(); 
      if (isSelected && !fPreSelCut){
//         cout << firedTrigClass.Data() << endl;
//         cout << "Special trigger: "<< fSpecialTrigger << " initialized " << fEMCALTrigInitialized << endl;
//         if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9){ // EMCAL triggers
//           if (!fEMCALTrigInitialized ) InitializeEMCALTrigger(event);
//           fTriggersEMCAL= GetTriggerList();
//         }
        if (fSpecialSubTrigger>0 && !isMC){
          if (!firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 0;
          if (fRejectTriggerOverlap){            
            // trigger rejection EMC1,7,8
            if (fSpecialTrigger == 5 && fSpecialSubTriggerName.CompareTo("CEMC7") == 0){
              if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
            } else if (fSpecialTrigger == 5 && fSpecialSubTriggerName.CompareTo("CEMC1") == 0){
              if (fInputHandler->IsEventSelected() & AliVEvent::kMB) isSelected = 0;
            } else if (fSpecialTrigger == 5 && fSpecialSubTriggerName.CompareTo("CEMC8") == 0){
              if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
            }
            // trigger rejection EGA
            if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("7EGA") == 0){
              if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
              if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
            } else if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("8EGA") == 0){
              if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
              if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
            }
            // trigger rejection EG1 & EG2
            if (fPeriodEnum == kLHC13g){
              // EG1 is the trigger with the highest threshold
              if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("7EG1") == 0){
//                 cout << firedTrigClass.Data() << endl;
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
//                 cout << "INT7? " << isSelected << endl;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
//                 cout << "CEM7? " << isSelected << endl;
                if (firedTrigClass.Contains("7EG2"))  isSelected = 0;
//                 cout << "7EG2? " << isSelected << endl;
              } else if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("8EG1") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                if (firedTrigClass.Contains("8EG2"))  isSelected = 0;
              } else   if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("7EG2") == 0){
//                 cout << firedTrigClass.Data() << endl;
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
//                 cout << "INT7? " << isSelected << endl;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
//                 cout << "CEM7? " << isSelected << endl;
              } else   if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("8EG2") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              }
            } else {
              // EG2 is the trigger with the highest threshold
              if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("7EG2") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                if (firedTrigClass.Contains("7EG1"))  isSelected = 0;
              } else if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("8EG2") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
                if (firedTrigClass.Contains("8EG1"))  isSelected = 0;
              } else   if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("7EG1") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              } else   if (fSpecialTrigger == 8 && fSpecialSubTriggerName.CompareTo("8EG1") == 0){
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) isSelected = 0;
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) isSelected = 0;
              }
            }
          }
          if (isSelected != 0 ){
//             cout << "I am here" << " :" << fSpecialSubTriggerName.Data() <<endl;
            if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9 ){
              if (hTriggerClassesCorrelated){
                if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClassesCorrelated->Fill(0);
                if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)hTriggerClassesCorrelated->Fill(1);
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)hTriggerClassesCorrelated->Fill(2);
                if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7)hTriggerClassesCorrelated->Fill(3);
                if (firedTrigClass.Contains("7EJE") || firedTrigClass.Contains("8EJE")) hTriggerClassesCorrelated->Fill(4);
                if (firedTrigClass.Contains("7EJ1") || firedTrigClass.Contains("8EJ1")) hTriggerClassesCorrelated->Fill(5);
                if (firedTrigClass.Contains("7EJ2") || firedTrigClass.Contains("8EJ2")) hTriggerClassesCorrelated->Fill(6);
                if (firedTrigClass.Contains("7EGA") || firedTrigClass.Contains("8EGA")) hTriggerClassesCorrelated->Fill(7);
                if (firedTrigClass.Contains("7EG1") || firedTrigClass.Contains("8EG1")) hTriggerClassesCorrelated->Fill(8);
                if (firedTrigClass.Contains("7EG2") || firedTrigClass.Contains("8EG2")) hTriggerClassesCorrelated->Fill(9);
              }
            }
          }
          
        } else if (isMC){
          if (fSpecialTrigger == 5 || fSpecialTrigger == 8 || fSpecialTrigger == 9){ // EMCAL triggers
//             isSelected = 0;
//             if (fTriggersEMCAL > 0)cout << "Special Trigger " << fSpecialTrigger << " triggers: " << fTriggersEMCAL << "    selected triggers: " << fTriggersEMCALSelected << " run number: " <<event->GetRunNumber()<<endl;
//             if (fTriggersEMCAL&fTriggersEMCALSelected){
//               cout << "accepted ++++++++++++++++++++" << endl;
              isSelected = 1;
//             }
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
  }
  fIsSDDFired = !(fInputHandler->IsEventSelected() & AliVEvent::kFastOnly);

  Bool_t mimickedTrigger = kTRUE;
  if (fMimicTrigger) mimickedTrigger = MimicTrigger(event, isMC);
//   cout << "mimicked decision \t" << mimickedTrigger << "expect decision? "<< fMimicTrigger<< endl;
  
  // Fill Histogram
  if(hTriggerClass){
    if (fIsSDDFired) hTriggerClass->Fill(33);
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
      if (fInputHandler->IsEventSelected() & AliVEvent::kAnyINT)hTriggerClass->Fill(30);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAny)hTriggerClass->Fill(31);
      if (!fInputHandler->IsEventSelected()) hTriggerClass->Fill(34);
    }
    if (mimickedTrigger && fMimicTrigger) hTriggerClass->Fill(35);
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
      fMCEvent = dynamic_cast<AliMCEvent*>(event);
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
    if(rejection == 2){ // TList of Headers Names
      for(Int_t i = 0; i<genHeaders->GetEntries();i++){
        gh                    = (AliGenEventHeader*)genHeaders->At(i);
        TString GeneratorName = gh->GetName();
        lastindexA            = lastindexA + gh->NProduced();
        if (fDebugLevel > 0 ) cout << i << "\t" << GeneratorName.Data() << endl;
        for(Int_t j = 0; j<HeaderList->GetEntries();j++){
          TString GeneratorInList   = ((TObjString*)HeaderList->At(j))->GetString();
          if (fDebugLevel > 0 )  cout << GeneratorInList.Data() << endl;
          if(GeneratorName.CompareTo(GeneratorInList) == 0){
            if (fDebugLevel > 0 ) cout << "accepted" << endl;
            if (GeneratorInList.CompareTo("PARAM") == 0 || GeneratorInList.CompareTo("BOX") == 0 ){
              if(fMCEvent){
                if (fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c ){
                  if (fMCEvent->Particle(firstindexA)->GetPdgCode() == fAddedSignalPDGCode ) {
                    if (gh->NProduced() > 10 && fMCEvent->Particle(firstindexA+10)->GetPdgCode() == fAddedSignalPDGCode ){
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
                    if (gh->NProduced() > 10){
                      AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(fMCEventAOD->At(firstindexA+10));
                      if (  aodMCParticle2->GetPdgCode() == fAddedSignalPDGCode ){
                        if (fDebugLevel > 0 ) cout << "cond 1: " << fnHeaders << endl;
                        fnHeaders++;
                        continue;
                      } 
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
            if (fDebugLevel > 0 ) cout << "cond 3: "<< fnHeaders << endl;
            fnHeaders++;
            continue;
          }
        }
        firstindexA       = firstindexA + gh->NProduced();
      }
    }
    if (fDebugLevel > 0 ) cout << "number of headers: " <<fnHeaders << endl;
    
    fNotRejectedStart         = new Int_t[fnHeaders];
    fNotRejectedEnd         = new Int_t[fnHeaders];
    fGeneratorNames         = new TString[fnHeaders];

    if(rejection == 1 || rejection == 3){
      fNotRejectedStart[0]     = 0;
      fNotRejectedEnd[0]       = ((AliGenEventHeader*)genHeaders->At(0))->NProduced()-1;
      fGeneratorNames[0]       = ((AliGenEventHeader*)genHeaders->At(0))->GetName();
      if (fDebugLevel > 0 ) cout << 0 << "\t" <<fGeneratorNames[0] << "\t" << fNotRejectedStart[0] << "\t" <<fNotRejectedEnd[0] << endl;
      return;
    }

    Int_t firstindex         = 0;
    Int_t lastindex         =  -1;
    Int_t number           = 0;
    
    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName     = gh->GetName();
      lastindex           = lastindex + gh->NProduced();
      for(Int_t j = 0; j<HeaderList->GetEntries();j++){
        TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
        if (fDebugLevel > 0 ) cout << i << "\t" << GeneratorName.Data() << endl;
        if(GeneratorName.CompareTo(GeneratorInList) == 0){
          if (GeneratorInList.CompareTo("PARAM") == 0 || GeneratorInList.CompareTo("BOX") == 0 ){
            if(fMCEvent){
              if (fPeriodEnum == kLHC14a1b || fPeriodEnum == kLHC14a1c ){
                if (fMCEvent->Particle(firstindex)->GetPdgCode() == fAddedSignalPDGCode ) {
                  if (fDebugLevel > 0 ) cout << "produced " << gh->NProduced() << " with box generator" << endl;
                  if (gh->NProduced() > 10 && fMCEvent->Particle(firstindex+10)->GetPdgCode() == fAddedSignalPDGCode){
                    if (fDebugLevel > 0 ) cout << "one of them was a pi0 or eta" <<  endl;
                    fNotRejectedStart[number] = firstindex;
                    fNotRejectedEnd[number] = lastindex;
                    fGeneratorNames[number] = GeneratorName;
                    number++;
                    if (fDebugLevel > 0 ) cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
                    continue;
                  }
                }
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
                  if (gh->NProduced() > 10) {
                    AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(fMCEventAOD->At(firstindex+10));
                    if ( aodMCParticle2->GetPdgCode() == fAddedSignalPDGCode ){
                      fNotRejectedEnd[number] = lastindex;
                      fNotRejectedStart[number] = firstindex;
                      fGeneratorNames[number] = GeneratorName;
                      number++;
                    } 
                    continue;
                  }
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
          } else {
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

//_________________________________________________________________________
Int_t AliConvEventCuts::IsParticleFromBGEvent(Int_t index, AliMCEvent *mcEvent, AliVEvent *InputEvent){

  if (fDebugLevel > 2 ) cout << index << endl;
  if(index < 0) return 0; // No Particle

  Int_t accepted = 0;
  if(!InputEvent || InputEvent->IsA()==AliESDEvent::Class()){
    if(!mcEvent) return 0; // no mcEvent available, return 0
    if(index >= mcEvent->GetNumberOfPrimaries()){ // initial particle is secondary particle
      if( ((TParticle*)mcEvent->Particle(index))->GetMother(0) < 0) return 0; // material particle, return 0
      return IsParticleFromBGEvent(((TParticle*)mcEvent->Particle(index))->GetMother(0),mcEvent,InputEvent);
    }
    for(Int_t i = 0;i<fnHeaders;i++){
      if (fDebugLevel > 2 ) cout << "header " << i << ":"<< fNotRejectedStart[i] << "\t" << fNotRejectedEnd[i] << endl;
      if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
        if (fDebugLevel > 1 ) cout << "accepted:" << index << "\t header " << i << ": "<< fNotRejectedStart[i] << "\t" << fNotRejectedEnd[i] << endl;
        accepted = 1;
        if(i == 0) accepted = 2; // MB Header
      }
    }
  }
  else if(InputEvent->IsA()==AliAODEvent::Class()){
    TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(InputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray){
      AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(index));
      if(!aodMCParticle) return 0; // no particle
      if(!aodMCParticle->IsPrimary()){
        if( aodMCParticle->GetMother() < 0) return 0;// material particle, return 0
        return IsParticleFromBGEvent(aodMCParticle->GetMother(),mcEvent,InputEvent);
      }
      index = TMath::Abs(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(index))->GetLabel());
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

  // Special EMCAL checks due to hardware issues in LHC11a  
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

  if( isHeavyIon != 2 && GetIsFromPileup()){
    if(event->IsPileupFromSPD(3,0.8,3.,2.,5.) ){
      if (hPileupVertexToPrimZSPDPileup) hPileupVertexToPrimZSPDPileup->Fill(distZMax);
      return 6; // Check Pileup --> Not Accepted => eventQuality = 6
    }
    if (fUtils->IsSPDClusterVsTrackletBG(event)){
      if (hPileupVertexToPrimZTrackletvsHits) hPileupVertexToPrimZTrackletvsHits->Fill(distZMax);
      return 11; // Check Pileup --> Not Accepted => eventQuality = 11
    }
  }
  if(isHeavyIon == 2 && GetIsFromPileup()){
    if(fUtils->IsPileUpEvent(event) ){
      if (hPileupVertexToPrimZSPDPileup) hPileupVertexToPrimZSPDPileup->Fill(distZMax);
      return 6; // Check Pileup --> Not Accepted => eventQuality = 6
    }
    if (fUtils->IsSPDClusterVsTrackletBG(event)){
      if (hPileupVertexToPrimZTrackletvsHits) hPileupVertexToPrimZTrackletvsHits->Fill(distZMax);
      return 11; // Check Pileup --> Not Accepted => eventQuality = 11
    }
  }

  if(GetIsFromPileup() && GetDoPileUpRejectV0MTPCout() ){
     if( IsPileUpV0MTPCout(event) ){
       return 13;
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
  if ( fDoCentralityFlat == 1 && (centrality >= 0. && centrality <= 10.) ){
    GetValueForWeight = hCentralityNotFlat->Interpolate(centrality);
    maximum = hCentralityNotFlat->GetMaximum();
    CorrCentrLoop = kTRUE;
  } else if ( fDoCentralityFlat == 2 && (centrality >=10. && centrality <= 20.) ){
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
  
  Float_t valueMultData       = -1.;
  Float_t valueMultMC         = -1.;
  
  if (hReweightMultData == NULL || hReweightMultMC == NULL ) return weightMult;
  
  valueMultData               = hReweightMultData->Interpolate(mult);
  valueMultMC                 = hReweightMultMC->Interpolate(mult);
  
  Float_t relativeErrorMC     = hReweightMultMC->GetBinError(hReweightMultMC->FindBin(mult))/hReweightMultMC->GetBinContent(hReweightMultMC->FindBin(mult));
  Float_t relativeErrorData   = hReweightMultData->GetBinError(hReweightMultData->FindBin(mult))/hReweightMultData->GetBinContent(hReweightMultData->FindBin(mult));
  
  if (relativeErrorData < 0.2 && relativeErrorMC < 0.2 ){
     if (isfinite(valueMultData) && isfinite(valueMultMC) ){
        weightMult               = valueMultData/valueMultMC;
     } 
  }
  
  return weightMult;
}



//_________________________________________________________________________
Float_t AliConvEventCuts::GetWeightForMeson(Int_t index, AliMCEvent *mcEvent, AliVEvent *event){
  if (!(  fPeriodEnum == kLHC13d2   || fPeriodEnum == kLHC13d2b       ||                                            // LHC10h MCs
          fPeriodEnum == kLHC14a1a  || fPeriodEnum == kLHC14a1b       || fPeriodEnum == kLHC14a1c   ||             // LHC11h MCs
          fPeriodEnum == kLHC13e7   || fPeriodEnum == kLHC13b2_efix   || fPeriodEnum == kLHC14b2      ||             // LHC13bc MCs
          fPeriodEnum == kLHC14e2a  || fPeriodEnum == kLHC14e2b       || fPeriodEnum == kLHC14e2c    ||             // LHC12[a-i] pass 1 MCs
          fPeriodEnum == kLHC12f1a  || fPeriodEnum == kLHC12f1b       || fPeriodEnum == kLHC12i3                    // LHC11a MCs
     ) ) return 1.;
  Int_t kCaseGen = 0;

  if(index < 0) return 0; // No Particle
    
  if (IsParticleFromBGEvent(index, mcEvent, event)){
    if (fPeriodEnum == kLHC13d2 || fPeriodEnum == kLHC13d2b || fPeriodEnum == kLHC13e7 || fPeriodEnum == kLHC13b2_efix || fPeriodEnum == kLHC14a1a || fPeriodEnum ==  kLHC14a1b || fPeriodEnum ==  kLHC14a1c       ||
      fPeriodEnum == kLHC14b2 || fPeriodEnum == kLHC14e2a || fPeriodEnum == kLHC14e2b || fPeriodEnum == kLHC14e2c || fPeriodEnum == kLHC12f1a || fPeriodEnum == kLHC12f1b || fPeriodEnum == kLHC12i3){
      kCaseGen = 1;
    }
  }
  if (kCaseGen == 0) return 1;

  Double_t mesonPt = 0;
  Double_t mesonMass = 0;
  Int_t PDGCode = 0;
  if(!event || event->IsA()==AliESDEvent::Class()){
    mesonPt = ((TParticle*)mcEvent->Particle(index))->Pt();
    mesonMass = ((TParticle*)mcEvent->Particle(index))->GetCalcMass();
    PDGCode = ((TParticle*)mcEvent->Particle(index))->GetPdgCode();
  } else if(event->IsA()==AliAODEvent::Class()){
    TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray){
      AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(index));
      mesonPt = aodMCParticle->Pt();
      mesonMass = aodMCParticle->GetCalcMass();
      PDGCode = aodMCParticle->GetPdgCode();
    } else {
      return 1;
    }
  }

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

  Double_t weight = 1;
  if (PDGCode ==  111 || PDGCode ==  221){
    if (functionResultData != 0. && functionResultMC != 0. && isfinite(functionResultData) && isfinite(functionResultMC)){
      weight = functionResultData/functionResultMC;
      if ( kCaseGen == 3){
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


///________________________________________________________________________
void AliConvEventCuts::GetCorrectEtaShiftFromPeriod(){

  if( fPeriodEnum == kLHC13bc ||                                      // mainly minimum bias
      fPeriodEnum == kLHC13de ||                                      // mainly triggered
      fPeriodEnum == kLHC13b4_fix || fPeriodEnum == kLHC13b4_plus ||  // MC Pythia 6 (Jet-Jet), anchor LHC13b-e
      fPeriodEnum == kLHC13b2_efix ||                                 //MC DPMJET, anchr LHC13b+c
      fPeriodEnum == kLHC13e7 ||                                      //MC HIJING, anchr LHC13b+c
      fPeriodEnum == kLHC14b2                                         //MC HIJING, anchr LHC13b+c
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
//   if (nPatch> 0) {cout << "NEW Triggers in this event*********************************" << endl;}
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
    patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
//     cout << "Patch energy: "<<patch->GetPatchE() << "\t ADC counts: " << patch->GetADCAmp() << endl;
//     cout << "Phi: " << patch->GetPhiMin() << " - " << patch->GetPhiMax() << " delta phi: " <<TMath::Abs(patch->GetPhiMin()-patch->GetPhiMax())<< endl;
//     cout << "Eta: " << patch->GetEtaMin() << " - " << patch->GetEtaMax() << " delta eta: " <<TMath::Abs(patch->GetEtaMin()-patch->GetEtaMax())<< endl;
    if (patch->IsGammaHigh()){
//       cout << "fired L1GA high" << endl;
      nG1++;
    }
    if (patch->IsGammaLow()){
//       cout << "fired L1GA low" << endl;
      nG2++;
    }
    if (patch->IsJetHigh()){
//       cout << "fired L1JE high" << endl;
      nJ1++;
    }
    if (patch->IsJetLow()){
//       cout << "fired L1JE low" << endl;
      nJ2++;
    }
    if (patch->IsLevel0()){
//       cout << "fired L0" << endl;
      nL0++;
    }
//     cout << patch->GetPatchE()   << "\t" << patch->GetADCAmp()  << "\t" << patch->IsGammaHigh() << "\t" << patch->IsGammaLow()  
//          << "\t" << patch->IsJetHigh()  << "\t" << patch->IsJetLow()  << "\t" << patch->IsLevel0() 
//        << "\t" << patch->GetPhiMin()  << "\t" << patch->GetPhiMax()  << "\t" << TMath::Abs(patch->GetPhiMin()-patch->GetPhiMax())
//        << "\t" << patch->GetEtaMin()  << "\t" << patch->GetEtaMax()  << "\t" << TMath::Abs(patch->GetEtaMin()-patch->GetEtaMax()) << endl;
  }

  if (nPatch > 0){
    AliDebug(2, "Patch summary: ");
    AliDebug(2, Form("Number of patches: %d", nPatch));
    AliDebug(2, Form("Level0: [%d]" ,nL0));
    AliDebug(2, Form("Jet:    low[%d], high[%d]" ,nJ2, nJ1));
    AliDebug(2, Form("Gamma:  low[%d], high[%d]" ,nG2, nG1));
  }
    
//   if (nPatch > 0){
//     cout <<     Form("Number of patches: %d", nPatch) << endl;
//     cout <<     Form("Level0: [%d]" ,nL0) << endl;
//     cout <<     Form("Jet:    low[%d], high[%d]" ,nJ2, nJ1) << endl;
//     cout <<     Form("Gamma:  low[%d], high[%d]" ,nG2, nG1) << endl;
//   }
    
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
//       cout << "dalitz candidate found" << endl;
    }
  
    Long_t source = particle->GetMother(0);
    Bool_t foundExcludedPart = kFALSE;
    Bool_t foundShower = kFALSE;
    Int_t pdgCodeMotherPrev = 0;
    Int_t pdgCodeMotherPPrevMother = 0;
    Int_t depth = 0;
    if (dalitzCand || realRadius3D < fSecProdBoundary ){
//       if (particle->GetPdgCode() == 22){
//         cout << endl << endl << "new particle: " << eventpos <<endl;
//         cout << particle->GetPdgCode() << "\t" << particle->R() << "\t" << realRadius2D << "\t" << realRadius3D << endl;
//       }
      while (depth < 20){
        TParticle* mother   = (TParticle *)mcEvent->Particle(source);
        source         = mother->GetMother(0); 
//         if (particle->GetPdgCode() == 22)cout << "eventposition: "<< source << endl;
        Int_t pdgCodeMother     = mother->GetPdgCode();
//         if (particle->GetPdgCode() == 22)cout << "Previous mothers: " << pdgCodeMother << "\t"<< pdgCodeMotherPrev<< "\t" << pdgCodeMotherPPrevMother << endl;
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
//         if (particle->GetPdgCode() == 22)cout << mother->GetPdgCode() << "\t" <<  source << "\t" << foundExcludedPart<< endl;
        pdgCodeMotherPPrevMother = pdgCodeMotherPrev;
        pdgCodeMotherPrev = pdgCodeMother;
        if (source == -1) depth = 20;
        
//         if (particle->GetPdgCode() == 22)cout << depth << endl;
        depth++;
      }
    }
    if (foundExcludedPart){
//       if (particle->GetPdgCode() == 22)cout << "This is definitely a secondary, manually excluded" << endl;
      return kFALSE;
    } else if (dalitzCand && realRadius3D < fSecProdBoundary ){
//       if (particle->GetPdgCode() == 22)cout << "This was a decay via a virtual photon" << endl;
      return kTRUE;
    } else if (foundShower){
//       if (particle->GetPdgCode() == 22)cout << "This is a shower" << endl;
      return kFALSE;
    } else if (realRadius3D >= fSecProdBoundary){
//       cout << "This is a secondary, to large production radius" << endl;
      return kFALSE;
    }
  }

  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliConvEventCuts::IsConversionPrimaryAOD(AliVEvent *event, AliAODMCParticle* AODMCParticle,  Double_t prodVtxX, Double_t prodVtxY, Double_t prodVtxZ){

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return kFALSE;
  AliAODMCParticle* currentParticle = AODMCParticle;
  if (TMath::Abs(currentParticle->GetPdgCode()) == 11 ){
    if (currentParticle->GetMother() != -1){
      AliAODMCParticle* particleMother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(currentParticle->GetMother()));
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

    AliAODMCParticle* firstmother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(currentParticle->GetMother()));
    if (!firstmother) return kFALSE;
    Int_t pdgCodeFirstMother = firstmother->GetPdgCode();
    Bool_t intDecay = kFALSE;
    if ( pdgCodeFirstMother == 111 || pdgCodeFirstMother == 221 ) intDecay = kTRUE;
    if ( intDecay && TMath::Abs(currentParticle->GetPdgCode()) == 11 ){
      dalitzCand = kTRUE;
//       cout << "dalitz candidate found" << endl;
    }

    Long_t source = currentParticle->GetMother();
    Bool_t foundExcludedPart = kFALSE;
    Bool_t foundShower = kFALSE;
    Int_t pdgCodeMotherPrev = 0;
    Int_t pdgCodeMotherPPrevMother = 0;
    Int_t depth = 0;
    if (dalitzCand || realRadius3D < fSecProdBoundary ){
//       if (currentParticle->GetPdgCode() == 22){
//         cout << endl << endl << "new particle: " << eventpos <<endl;
//         cout << currentParticle->GetPdgCode() << "\t" << currentParticle->R() << "\t" << realRadius2D << "\t" << realRadius3D << endl;
//       }
      while (depth < 20){
        AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(source));
        source = mother->GetMother();
//         if (currentParticle->GetPdgCode() == 22)cout << "eventposition: "<< source << endl;
        Int_t pdgCodeMother     = mother->GetPdgCode();
//         if (currentParticle->GetPdgCode() == 22)cout << "Previous mothers: " << pdgCodeMother << "\t"<< pdgCodeMotherPrev<< "\t" << pdgCodeMotherPPrevMother << endl;
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
//      if (currentParticle->GetPdgCode() == 22)cout << mother->GetPdgCode() << "\t" <<  source << "\t" << foundExcludedPart<< endl;
        pdgCodeMotherPPrevMother = pdgCodeMotherPrev;
        pdgCodeMotherPrev = pdgCodeMother;
        if (source == -1) depth = 20;

//      if (currentParticle->GetPdgCode() == 22)cout << depth << endl;
        depth++;
      }
    }
    if (foundExcludedPart){
//      if (currentParticle->GetPdgCode() == 22)cout << "This is definitely a secondary, manually excluded" << endl;
      return kFALSE;
    } else if (dalitzCand && realRadius3D < fSecProdBoundary ){
//      if (currentParticle->GetPdgCode() == 22)cout << "This was a decay via a virtual photon" << endl;
      return kTRUE;
    } else if (foundShower){
//      if (currentParticle->GetPdgCode() == 22)cout << "This is a shower" << endl;
      return kFALSE;
    } else if (realRadius3D >= fSecProdBoundary){
//      cout << "This is a secondary, too large production radius" << endl;
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
      Bool_t hasMother        = kFALSE;
      Bool_t hasGrandMother   = kFALSE;
      Long_t motherID         = particle->GetMother(0);
      Long_t grandMotherID    = -1;
      // is the photon a direct photons, without a mother?
      if (motherID > -1){
        hasMother             = kTRUE;
        grandMotherID         = mcEvent->Particle(motherID)->GetMother(0);
        // is the meson a primary?
        if (grandMotherID > -1){
          hasGrandMother      = kTRUE;
          pdgSecondary        = mcEvent->Particle(grandMotherID)->GetPdgCode();
        }
      }
    } else {
      Bool_t hasMother            = kFALSE;
      Bool_t hasGrandMother       = kFALSE;
      Bool_t hasGreatGrandMother  = kFALSE;
      Long_t motherID             = particle->GetMother(0);
      Long_t grandMotherID        = -1;
      Long_t greatGrandMotherID   = -1;
      // is the electron a direct electron, without a mother?
      if (motherID > -1){
        hasMother                 = kTRUE;
        grandMotherID             = mcEvent->Particle(motherID)->GetMother(0);
        // is the photon a direct photons, without a mother?
        if (grandMotherID > -1){
          hasGrandMother          = kTRUE;
          greatGrandMotherID      = mcEvent->Particle(grandMotherID)->GetMother(0);
          // is the meson a primary?
          if (greatGrandMotherID > -1){
            hasGreatGrandMother   = kTRUE;
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
      Bool_t hasMother        = kFALSE;
      Bool_t hasGrandMother   = kFALSE;
      Long_t motherID         = particle->GetMother();
      Long_t grandMotherID    = -1;
      // is the photon a direct photons, without a mother?
      if (motherID > -1){
        hasMother             = kTRUE;
        grandMotherID         = ((AliAODMCParticle*)aodmcArray->At(motherID))->GetMother();
        // is the meson a primary?
        if (grandMotherID > -1){
          hasGrandMother      = kTRUE;
          pdgSecondary        = ((AliAODMCParticle*)aodmcArray->At(grandMotherID))->GetPdgCode();
        }
      }
    } else {
      Bool_t hasMother            = kFALSE;
      Bool_t hasGrandMother       = kFALSE;
      Bool_t hasGreatGrandMother  = kFALSE;
      Long_t motherID             = particle->GetMother();
      Long_t grandMotherID        = -1;
      Long_t greatGrandMotherID   = -1;
      // is the electron a direct electron, without a mother?
      if (motherID > -1){
        hasMother                 = kTRUE;
        grandMotherID             = ((AliAODMCParticle*)aodmcArray->At(motherID))->GetMother();
        // is the photon a direct photons, without a mother?
        if (grandMotherID > -1){
          hasGrandMother          = kTRUE;
          greatGrandMotherID      = ((AliAODMCParticle*)aodmcArray->At(grandMotherID))->GetMother();
          // is the meson a primary?
          if (greatGrandMotherID > -1){
            hasGreatGrandMother   = kTRUE;
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
    periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPeriodName();
  } 
  
  if (periodName.CompareTo("") == 0) {
    fPeriodEnum = kNoPeriod;
    fEnergyEnum = kUnset;
    AliError("No correct period could be set, periodName string empty");
    return;
  }
  
  // Data
  if (periodName.CompareTo("LHC10b") == 0 || periodName.CompareTo("LHC10c") == 0 || periodName.CompareTo("LHC10d") == 0 || periodName.CompareTo("LHC10e") == 0 ||
      periodName.CompareTo("LHC10f") == 0 || periodName.CompareTo("LHC10g") == 0 
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
            periodName.CompareTo("LHC12i") == 0 
  ) {
    fPeriodEnum = kLHC12;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC13b") == 0 || periodName.CompareTo("LHC13c") == 0 ){
    fPeriodEnum = kLHC13bc;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13d") == 0 || periodName.CompareTo("LHC13e") == 0 ){
    fPeriodEnum = kLHC13de;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13f") == 0 ){
    fPeriodEnum = kLHC13f;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13g") == 0 ){
    fPeriodEnum = kLHC13g;
    fEnergyEnum = k2760GeV;
  } else if (periodName.CompareTo("LHC15f") == 0 || periodName.CompareTo("LHC15g") == 0 || periodName.CompareTo("LHC15h") == 0 || periodName.CompareTo("LHC15i") == 0 ||
            periodName.CompareTo("LHC15j") == 0 || periodName.CompareTo("LHC15k") == 0 || periodName.CompareTo("LHC15l") == 0 || periodName.CompareTo("LHC15m") == 0 
  ) {
    fPeriodEnum = kLHC15fm;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC15n") == 0 ){
    fPeriodEnum = kLHC15n;
    fEnergyEnum = k5TeV;
  } else if (periodName.CompareTo("LHC15o") == 0 ){
    fPeriodEnum = kLHC15o;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.CompareTo("LHC16k") == 0 || periodName.CompareTo("LHC16l") == 0 ){
    fPeriodEnum = kLHC16kl;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16d") == 0 ){
    fPeriodEnum = kLHC16d;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16e") == 0 ){
    fPeriodEnum = kLHC16e;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16f") == 0 ){
    fPeriodEnum = kLHC16f;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16g") == 0 ){
    fPeriodEnum = kLHC16g;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16h") == 0 ){
    fPeriodEnum = kLHC16h;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16i") == 0 ){
    fPeriodEnum = kLHC16i;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16j") == 0 ){
    fPeriodEnum = kLHC16j;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16o") == 0 ){
    fPeriodEnum = kLHC16o;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16p") == 0 ){
    fPeriodEnum = kLHC16p;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16q") == 0 ){
    fPeriodEnum = kLHC16q;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16r") == 0 ){
    fPeriodEnum = kLHC16r;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC16s") == 0 ){
    fPeriodEnum = kLHC16s;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC16t") == 0 ){
    fPeriodEnum = kLHC16t;
    fEnergyEnum = kpPb5TeV;
    
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
  } else if (periodName.CompareTo("LHC14e2a") == 0){
    fPeriodEnum = kLHC14e2a;   
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC14e2b") == 0){
    fPeriodEnum = kLHC14e2b;   
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC14e2c") == 0){
    fPeriodEnum = kLHC14e2c;
    fEnergyEnum = k8TeV;
  } else if (periodName.Contains("LHC15h1")){
    fPeriodEnum = kLHC15h1;
    fEnergyEnum = k8TeV;
  } else if (periodName.Contains("LHC15h2")){
    fPeriodEnum = kLHC15h2;
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC16c2") == 0){
    fPeriodEnum = kLHC16c2;   
    fEnergyEnum = k8TeV;
  } else if (periodName.CompareTo("LHC16c2_plus") == 0){
    fPeriodEnum = kLHC16c2_plus;
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
  } else if (periodName.CompareTo("LHC13b4_fix") == 0){
    fPeriodEnum = kLHC13b4_fix;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC13b4_plus") == 0){
    fPeriodEnum = kLHC13b4_plus;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16c3a") == 0){
    fPeriodEnum = kLHC16c3a;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16c3b") == 0){
    fPeriodEnum = kLHC16c3b;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC16c3c") == 0){
    fPeriodEnum = kLHC16c3c;
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
  } else if (periodName.Contains("LHC16h4")){
    fPeriodEnum = kLHC16h4;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g1")){
    fPeriodEnum = kLHC16g1;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g1a")){
    fPeriodEnum = kLHC16g1a;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g1b")){
    fPeriodEnum = kLHC16g1b;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g1c")){
    fPeriodEnum = kLHC16g1c;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g2")){
    fPeriodEnum = kLHC16g2;
    fEnergyEnum = kPbPb5TeV;
  } else if (periodName.Contains("LHC16g3")){
    fPeriodEnum = kLHC16g3;
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
  // LHC16x anchored MCs
  } else if (periodName.CompareTo("LHC16j2a1") == 0){
    fPeriodEnum = kLHC16j2a1;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16j2b1") == 0){
    fPeriodEnum = kLHC16j2b1;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16j2a2") == 0){
    fPeriodEnum = kLHC16j2a2;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC16j2b2") == 0){
    fPeriodEnum = kLHC16j2b2;
    fEnergyEnum = k13TeV;

  } else if (periodName.CompareTo("LHC17f6") == 0){
    fPeriodEnum = kLHC17f6;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f9") == 0){
    fPeriodEnum = kLHC17f9;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d1") == 0){
    fPeriodEnum = kLHC17d1;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d17") == 0){
    fPeriodEnum = kLHC17d17;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f5") == 0){
    fPeriodEnum = kLHC17f5;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d3") == 0){
    fPeriodEnum = kLHC17d3;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17e5") == 0){
    fPeriodEnum = kLHC17e5;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d20a1") == 0){
    fPeriodEnum = kLHC17d20a1;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d20a1_extra") == 0){
    fPeriodEnum = kLHC17d20a1_extra;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d20a2") == 0){
    fPeriodEnum = kLHC17d20a2;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d20a2_extra") == 0){
    fPeriodEnum = kLHC17d20a2_extra;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d16") == 0){
    fPeriodEnum = kLHC17d16;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17d18") == 0){
    fPeriodEnum = kLHC17d18;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f8a") == 0){
    fPeriodEnum = kLHC17f8a;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f8b") == 0){
    fPeriodEnum = kLHC17f8b;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f8c") == 0){
    fPeriodEnum = kLHC17f8c;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f8d") == 0){
    fPeriodEnum = kLHC17f8d;
    fEnergyEnum = k13TeV;
  } else if (periodName.CompareTo("LHC17f8e") == 0){
    fPeriodEnum = kLHC17f8e;
    fEnergyEnum = k13TeV;
  // LHC16qt anchored MCs
  } else if (periodName.CompareTo("LHC17a2a") == 0){
    fPeriodEnum = kLHC17a2a;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2a_fast") == 0){
    fPeriodEnum = kLHC17a2a_fast;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2a_cent") == 0){
    fPeriodEnum = kLHC17a2a_cent;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a2a_cent_woSDD;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2b") == 0){
    fPeriodEnum = kLHC17a2b;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2b_fast") == 0){
    fPeriodEnum = kLHC17a2b_fast;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2b_cent") == 0){
    fPeriodEnum = kLHC17a2b_cent;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17a2b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a2b_cent_woSDD;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2a") == 0){
    fPeriodEnum = kLHC17f2a;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2a_fast") == 0){
    fPeriodEnum = kLHC17f2a_fast;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2a_cent") == 0){
    fPeriodEnum = kLHC17f2a_cent;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17f2a_cent_woSDD;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2b") == 0){
    fPeriodEnum = kLHC17f2b;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2b_fast") == 0){
    fPeriodEnum = kLHC17f2b_fast;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2b_cent") == 0){
    fPeriodEnum = kLHC17f2b_cent;
    fEnergyEnum = kpPb5TeV;
  } else if (periodName.CompareTo("LHC17f2b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17f2b_cent_woSDD;
    fEnergyEnum = kpPb5TeV;
  // LHC16r anchored MCs
  } else if (periodName.CompareTo("LHC17a3a") == 0){
    fPeriodEnum = kLHC17a3a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3a_fast") == 0){
    fPeriodEnum = kLHC17a3a_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3a_cent") == 0){
    fPeriodEnum = kLHC17a3a_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a3a_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3b") == 0){
    fPeriodEnum = kLHC17a3b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3b_fast") == 0){
    fPeriodEnum = kLHC17a3b_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3b_cent") == 0){
    fPeriodEnum = kLHC17a3b_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a3b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a3b_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3a") == 0){
    fPeriodEnum = kLHC17f3a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3a_fast") == 0){
    fPeriodEnum = kLHC17f3a_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3a_cent") == 0){
    fPeriodEnum = kLHC17f3a_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17f3a_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3b") == 0){
    fPeriodEnum = kLHC17f3b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3b_fast") == 0){
    fPeriodEnum = kLHC17f3b_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3b_cent") == 0){
    fPeriodEnum = kLHC17f3b_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f3b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17f3b_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  // LHC16s anchored MCs
  } else if (periodName.CompareTo("LHC17a4a") == 0){
    fPeriodEnum = kLHC17a4a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4a_fast") == 0){
    fPeriodEnum = kLHC17a4a_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4a_cent") == 0){
    fPeriodEnum = kLHC17a4a_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a4a_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4b") == 0){
    fPeriodEnum = kLHC17a4b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4b_fast") == 0){
    fPeriodEnum = kLHC17a4b_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4b_cent") == 0){
    fPeriodEnum = kLHC17a4b_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17a4b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17a4b_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4a") == 0){
    fPeriodEnum = kLHC17f4a;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4a_fast") == 0){
    fPeriodEnum = kLHC17f4a_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4a_cent") == 0){
    fPeriodEnum = kLHC17f4a_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4a_cent_woSDD") == 0){
    fPeriodEnum = kLHC17f4a_cent_woSDD;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4b") == 0){
    fPeriodEnum = kLHC17f4b;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4b_fast") == 0){
    fPeriodEnum = kLHC17f4b_fast;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4b_cent") == 0){
    fPeriodEnum = kLHC17f4b_cent;
    fEnergyEnum = kpPb8TeV;
  } else if (periodName.CompareTo("LHC17f4b_cent_woSDD") == 0){
    fPeriodEnum = kLHC17f4b_cent_woSDD;
    fEnergyEnum = kpPb8TeV;


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
