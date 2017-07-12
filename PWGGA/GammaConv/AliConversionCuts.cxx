/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *					                                  *
 * Authors: Svein Lindal, Daniel Lohner 		                  *
 * Version 1.0                        					  *
 *                           						  *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is    	  *
 * provided "as is" without express or implied warranty.       		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
// WARNING: this class is no longer supported,
// please use AliConversionPhotonCuts and AliConvEventCuts
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
#include "AliConversionCuts.h"

#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
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
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliTRDTriggerAnalysis.h"

class iostream;

using namespace std;

ClassImp(AliConversionCuts)

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
// WARNING: this class is no longer supported,
// please use AliConversionPhotonCuts and AliConvEventCuts
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

const char* AliConversionCuts::fgkCutNames[AliConversionCuts::kNCuts] = {
   "HeavyIon",//0
   "CentralityMin",//1
   "CentralityMax",//2
   "SelectSpecialTrigger",//3
   "SelectSpecialSubTriggerClass",//4
   "RemovePileUp",//5
   "RejectExtraSignals",//6
   "V0FinderType",//7
   "EtaCut",//8
   "MinRCut",//9
   "SinglePtCut",//10
   "ClsTPCCut", //11
   "ededxSigmaCut",//12
   "pidedxSigmaCut",//13
   "piMomdedxSigmaCut",//14
   "piMaxMomdedxSigmaCut",//15
   "LowPRejectionSigmaCut",//16
   "TOFelectronPID",//17
   "QtMaxCut",//18
   "Chi2GammaCut", //19
   "PsiPair", //20
   "DoPhotonAsymmetryCut",//21
   "CosinePointingAngle", //22
   "SharedElectronCuts", //23
   "RejectToCloseV0s", //24
   "DcaRPrimVtx", //25
   "DcaZPrimVtx", //26
   "EvetPlane" //27
};


//________________________________________________________________________
AliConversionCuts::AliConversionCuts(const char *name,const char *title) :
   AliAnalysisCuts(name,title),
   fHistograms(NULL),
   fHeaderList(NULL),
   fPIDResponse(NULL),
   fEventQuality(-1),
   fMaxR(200),
   fMinR(0),
   fEtaCut(0.9),
   fEtaCutMin(-0.1),
   fPtCut(0.02),
   fSinglePtCut(0),
   fMaxZ(1000),
   fMinClsTPC(0.),
   fMinClsTPCToF(0.),
   fLineCutZRSlope(0.),
   fLineCutZValue(0),
   fLineCutZRSlopeMin(0.),
   fLineCutZValueMin(0),
   fChi2CutConversion(1000),
   fPIDProbabilityCutNegativeParticle(0),
   fPIDProbabilityCutPositiveParticle(0),
   fDodEdxSigmaCut(kTRUE),
   fDoTOFsigmaCut(kFALSE),
   fPIDTRDEfficiency(1),
   fDoTRDPID(kFALSE),
   fPIDnSigmaAboveElectronLine(100),
   fPIDnSigmaBelowElectronLine(-100),
   fTofPIDnSigmaAboveElectronLine(100),
   fTofPIDnSigmaBelowElectronLine(-100),
   fPIDnSigmaAbovePionLine(0),
   fPIDnSigmaAbovePionLineHighPt(-100),
   fPIDMinPnSigmaAbovePionLine(0),
   fPIDMaxPnSigmaAbovePionLine(0),
   fDoKaonRejectionLowP(kFALSE),
   fDoProtonRejectionLowP(kFALSE),
   fDoPionRejectionLowP(kFALSE),
   fPIDnSigmaAtLowPAroundKaonLine(0),
   fPIDnSigmaAtLowPAroundProtonLine(0),
   fPIDnSigmaAtLowPAroundPionLine(0),
   fPIDMinPKaonRejectionLowP(1.5),
   fPIDMinPProtonRejectionLowP(2),
   fPIDMinPPionRejectionLowP(0),
   fDoQtGammaSelection(kTRUE),
   fDo2DQt(kFALSE),
   fQtMax(100),
   fXVertexCut(0.),
   fYVertexCut(0.),
   fZVertexCut(0.),
   fNSigmaMass(0.),
   fUseEtaMinCut(kFALSE),
   fUseOnFlyV0Finder(kTRUE),
   fDoPhotonAsymmetryCut(kTRUE),
   fMinPPhotonAsymmetryCut(100.),
   fMinPhotonAsymmetry(0.),
   fIsHeavyIon(0),
   fDetectorCentrality(0),
   fModCentralityClass(0),
   fMaxVertexZ(10),
   fCentralityMin(0),
   fCentralityMax(0),
   fUseCorrectedTPCClsInfo(kFALSE),
   fUseTOFpid(kFALSE),
   fMultiplicityMethod(0),
   fSpecialTrigger(0),
   fSpecialSubTrigger(0),
   fRemovePileUp(kFALSE),
   fOpeningAngle(0.005),
   fPsiPairCut(10000),
   fDo2DPsiPairChi2(kFALSE),
   fCosPAngleCut(10000),
   fDoToCloseV0sCut(kFALSE),
   fRejectExtraSignals(0),
   fminV0Dist(200.),
   fDoSharedElecCut(kFALSE),
   fDoPhotonQualitySelectionCut(kFALSE),
   fPhotonQualityCut(0),
   fOfflineTriggerMask(0),
   fHasV0AND(kTRUE),
   fIsSDDFired(kTRUE),
   fRandom(0),
   fElectronArraySize(500),
   fElectronLabelArray(NULL),
   fDCAZPrimVtxCut(1000),
   fDCARPrimVtxCut(1000),
   fInPlaneOutOfPlane(0),
   fConversionPointXArray(0.0),
   fConversionPointYArray(0.0),
   fConversionPointZArray(0.0),
   fnHeaders(0),
   fNotRejectedStart(NULL),
   fNotRejectedEnd(NULL),
   fGeneratorNames(NULL),
   fCutString(NULL),
   fCutStringRead(""),
   fUtils(NULL),
   fEtaShift(0.0),
   fDoEtaShift(kFALSE),
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
   hEtaDistV0s(NULL),
   hEtaDistV0sAfterdEdxCuts(NULL),
   hdEdxCuts(NULL),
   hTPCdEdxbefore(NULL),
   hTPCdEdxafter(NULL),
   hTPCdEdxSigbefore(NULL),
   hTPCdEdxSigafter(NULL),
   hTOFbefore(NULL),
   hTOFSigbefore(NULL),
   hTOFSigafter(NULL),
   hPsiPairDeltaPhiafter(NULL),
   hTrackCuts(NULL),
   hPhotonCuts(NULL),
   hInvMassbefore(NULL),
   hArmenterosbefore(NULL),
   hInvMassafter(NULL),
   hArmenterosafter(NULL),
   hAcceptanceCuts(NULL),
   hCutIndex(NULL),
   hV0EventCuts(NULL),
   hCentrality(NULL),
   hCentralityVsNumberOfPrimaryTracks(NULL),
   hVertexZ(NULL),
   hEventPlanePhi(NULL),
   hTriggerClass(NULL),
   hTriggerClassSelected(NULL),
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
   fNSpecialSubTriggerOptions(0)

{
   InitPIDResponse();
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
   fCutString=new TObjString((GetCutNumber()).Data());

   fElectronLabelArray = new Int_t[fElectronArraySize];
   fUtils = new AliAnalysisUtils();
   //if you do not want to apply the cut on the distance between the SPD and TRK vertex:
   //fUtils->SetCutOnZVertexSPD(kFALSE);


}

//________________________________________________________________________
AliConversionCuts::AliConversionCuts(const AliConversionCuts &ref) :
   AliAnalysisCuts(ref),
   fHistograms(NULL),
   fHeaderList(ref.fHeaderList),
   fPIDResponse(NULL),
   fEventQuality(ref.fEventQuality),
   fMaxR(ref.fMaxR),
   fMinR(ref.fMinR),
   fEtaCut(ref.fEtaCut),
   fEtaCutMin(ref.fEtaCutMin),
   fPtCut(ref.fPtCut),
   fSinglePtCut(ref.fSinglePtCut),
   fMaxZ(ref.fMaxZ),
   fMinClsTPC(ref.fMinClsTPC),
   fMinClsTPCToF(ref.fMinClsTPCToF),
   fLineCutZRSlope(ref.fLineCutZRSlope),
   fLineCutZValue(ref.fLineCutZValue),
   fLineCutZRSlopeMin(ref.fLineCutZRSlopeMin),
   fLineCutZValueMin(ref.fLineCutZValueMin),
   fChi2CutConversion(ref.fChi2CutConversion),
   fPIDProbabilityCutNegativeParticle(ref.fPIDProbabilityCutNegativeParticle),
   fPIDProbabilityCutPositiveParticle(ref.fPIDProbabilityCutPositiveParticle),
   fDodEdxSigmaCut(ref. fDodEdxSigmaCut),
   fDoTOFsigmaCut(ref.fDoTOFsigmaCut),
   fPIDTRDEfficiency(ref.fPIDTRDEfficiency),
   fDoTRDPID(ref.fDoTRDPID),
   fPIDnSigmaAboveElectronLine(ref.fPIDnSigmaAboveElectronLine),
   fPIDnSigmaBelowElectronLine(ref.fPIDnSigmaBelowElectronLine),
   fTofPIDnSigmaAboveElectronLine(ref.fTofPIDnSigmaAboveElectronLine),
   fTofPIDnSigmaBelowElectronLine(ref.fTofPIDnSigmaBelowElectronLine),
   fPIDnSigmaAbovePionLine(ref.fPIDnSigmaAbovePionLine),
   fPIDnSigmaAbovePionLineHighPt(ref.fPIDnSigmaAbovePionLineHighPt),
   fPIDMinPnSigmaAbovePionLine(ref.fPIDMinPnSigmaAbovePionLine),
   fPIDMaxPnSigmaAbovePionLine(ref.fPIDMaxPnSigmaAbovePionLine),
   fDoKaonRejectionLowP(ref.fDoKaonRejectionLowP),
   fDoProtonRejectionLowP(ref.fDoProtonRejectionLowP),
   fDoPionRejectionLowP(ref.fDoPionRejectionLowP),
   fPIDnSigmaAtLowPAroundKaonLine(ref.fPIDnSigmaAtLowPAroundKaonLine),
   fPIDnSigmaAtLowPAroundProtonLine(ref.fPIDnSigmaAtLowPAroundProtonLine),
   fPIDnSigmaAtLowPAroundPionLine(ref.fPIDnSigmaAtLowPAroundPionLine),
   fPIDMinPKaonRejectionLowP(ref.fPIDMinPKaonRejectionLowP),
   fPIDMinPProtonRejectionLowP(ref.fPIDMinPProtonRejectionLowP),
   fPIDMinPPionRejectionLowP(ref.fPIDMinPPionRejectionLowP),
   fDoQtGammaSelection(ref.fDoQtGammaSelection),
   fDo2DQt(ref.fDo2DQt),
   fQtMax(ref.fQtMax),
   fXVertexCut(ref.fXVertexCut),
   fYVertexCut(ref.fYVertexCut),
   fZVertexCut(ref.fZVertexCut),
   fNSigmaMass(ref.fNSigmaMass),
   fUseEtaMinCut(ref.fUseEtaMinCut),
   fUseOnFlyV0Finder(ref.fUseOnFlyV0Finder),
   fDoPhotonAsymmetryCut(ref.fDoPhotonAsymmetryCut),
   fMinPPhotonAsymmetryCut(ref.fMinPPhotonAsymmetryCut),
   fMinPhotonAsymmetry(ref.fMinPhotonAsymmetry),
   fIsHeavyIon(ref.fIsHeavyIon),
   fDetectorCentrality(ref.fDetectorCentrality),
   fModCentralityClass(ref.fModCentralityClass),
   fMaxVertexZ(ref.fMaxVertexZ),
   fCentralityMin(ref.fCentralityMin),
   fCentralityMax(ref.fCentralityMax),
   fUseCorrectedTPCClsInfo(ref.fUseCorrectedTPCClsInfo),
   fUseTOFpid(ref.fUseTOFpid),
   fMultiplicityMethod(ref.fMultiplicityMethod),
   fSpecialTrigger(ref.fSpecialTrigger),
   fSpecialSubTrigger(ref.fSpecialSubTrigger),
   fRemovePileUp(ref.fRemovePileUp),
   fOpeningAngle(ref.fOpeningAngle),
   fPsiPairCut(ref.fPsiPairCut),
   fDo2DPsiPairChi2(ref.fDo2DPsiPairChi2),
   fCosPAngleCut(ref.fCosPAngleCut),
   fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
   fRejectExtraSignals(ref.fRejectExtraSignals),
   fminV0Dist(ref.fminV0Dist),
   fDoSharedElecCut(ref.fDoSharedElecCut),
   fDoPhotonQualitySelectionCut(ref.fDoPhotonQualitySelectionCut),
   fPhotonQualityCut(ref.fPhotonQualityCut),
   fOfflineTriggerMask(ref.fOfflineTriggerMask),
   fHasV0AND(ref.fHasV0AND),
   fIsSDDFired(ref.fIsSDDFired),
   fRandom(ref.fRandom),
   fElectronArraySize(ref.fElectronArraySize),
   fElectronLabelArray(NULL),
   fDCAZPrimVtxCut(ref.fDCAZPrimVtxCut),
   fDCARPrimVtxCut(ref.fDCAZPrimVtxCut),
   fInPlaneOutOfPlane(ref.fInPlaneOutOfPlane),
   fConversionPointXArray(ref.fConversionPointXArray),
   fConversionPointYArray(ref.fConversionPointYArray),
   fConversionPointZArray(ref.fConversionPointZArray),
   fnHeaders(ref.fnHeaders),
   fNotRejectedStart(NULL),
   fNotRejectedEnd(NULL),
   fGeneratorNames(ref.fGeneratorNames),
   fCutString(NULL),
   fCutStringRead(""),
   fUtils(NULL),
   fEtaShift(ref.fEtaShift),
   fDoEtaShift(ref.fDoEtaShift),
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
   hEtaDistV0s(NULL),
   hEtaDistV0sAfterdEdxCuts(NULL),
   hdEdxCuts(NULL),
   hTPCdEdxbefore(NULL),
   hTPCdEdxafter(NULL),
   hTPCdEdxSigbefore(NULL),
   hTPCdEdxSigafter(NULL),
   hTOFbefore(NULL),
   hTOFSigbefore(NULL),
   hTOFSigafter(NULL),
   hPsiPairDeltaPhiafter(NULL),
   hTrackCuts(NULL),
   hPhotonCuts(NULL),
   hInvMassbefore(NULL),
   hArmenterosbefore(NULL),
   hInvMassafter(NULL),
   hArmenterosafter(NULL),
   hAcceptanceCuts(NULL),
   hCutIndex(NULL),
   hV0EventCuts(NULL),
   hCentrality(NULL),
   hCentralityVsNumberOfPrimaryTracks(NULL),
   hVertexZ(NULL),
   hEventPlanePhi(NULL),
   hTriggerClass(NULL),
   hTriggerClassSelected(NULL),
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
   fNSpecialSubTriggerOptions(ref.fNSpecialSubTriggerOptions)
{
   // Copy Constructor
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
   fCutString=new TObjString((GetCutNumber()).Data());
   fElectronLabelArray = new Int_t[fElectronArraySize];
   fUtils = new AliAnalysisUtils();
   // dont copy histograms (if you like histograms, call InitCutHistograms())

}


//________________________________________________________________________
AliConversionCuts::~AliConversionCuts() {
   // Destructor
   //Deleting fHistograms leads to seg fault it it's added to output collection of a task
   // if(fHistograms)
   //    delete fHistograms;
   // fHistograms = NULL;
   if(fCutString != NULL){
      delete fCutString;
      fCutString = NULL;
   }
   if(fElectronLabelArray){
      delete fElectronLabelArray;
      fElectronLabelArray = NULL;
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
void AliConversionCuts::InitCutHistograms(TString name, Bool_t preCut){

   // Initialize Cut Histograms for QA (only initialized and filled if function is called)
   TH1::AddDirectory(kFALSE);

   if(fHistograms != NULL){
      delete fHistograms;
      fHistograms=NULL;
   }
   if(fHistograms==NULL){
      fHistograms=new TList();
      fHistograms->SetOwner(kTRUE);
      if(name=="")fHistograms->SetName(Form("ConvCuts_%s",GetCutNumber().Data()));
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
//    if (fFitDataPi0){
//       fFitDataPi0->SetName("DataFitForWeightingPi0");
//       fHistograms->Add(fFitDataPi0);
//    }
//    if (fFitDataEta){
//       fFitDataEta->SetName("DataFitForWeightingEta");
//       fHistograms->Add(fFitDataEta);
//    }
//    if (fFitDataK0s){
//       fFitDataK0s->SetName("DataFitForWeightingK0s");
//       fHistograms->Add(fFitDataK0s);
//    }
   // IsPhotonSelected
   hCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",11,-0.5,10.5);
   hCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
   hCutIndex->GetXaxis()->SetBinLabel(kOnFly+1,"onfly");
   hCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
   hCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"dEdx");
   hCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
   hCutIndex->GetXaxis()->SetBinLabel(kConvPointFail+1,"ConvPoint fail");
   hCutIndex->GetXaxis()->SetBinLabel(kPhotonCuts+1,"PhotonCuts");
   hCutIndex->GetXaxis()->SetBinLabel(kEventPlane+1,"EventPlane");
   hCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
   fHistograms->Add(hCutIndex);

   // Track Cuts
   hTrackCuts=new TH1F(Form("TrackCuts %s",GetCutNumber().Data()),"TrackCuts",9,-0.5,8.5);
   hTrackCuts->GetXaxis()->SetBinLabel(1,"in");
   hTrackCuts->GetXaxis()->SetBinLabel(2,"likesign");
   hTrackCuts->GetXaxis()->SetBinLabel(3,"ntpccl");
   hTrackCuts->GetXaxis()->SetBinLabel(4,"acceptance");
   hTrackCuts->GetXaxis()->SetBinLabel(5,"singlept");
   hTrackCuts->GetXaxis()->SetBinLabel(6,"TPCrefit");
   hTrackCuts->GetXaxis()->SetBinLabel(7,"kink");
   hTrackCuts->GetXaxis()->SetBinLabel(8,"out");
   fHistograms->Add(hTrackCuts);

   // Photon Cuts
   hPhotonCuts=new TH1F(Form("PhotonCuts %s",GetCutNumber().Data()),"PhotonCuts",15,-0.5,14.5);
   hPhotonCuts->GetXaxis()->SetBinLabel(1,"in");
   hPhotonCuts->GetXaxis()->SetBinLabel(2,"qtcut");
   hPhotonCuts->GetXaxis()->SetBinLabel(3,"chi2");
   hPhotonCuts->GetXaxis()->SetBinLabel(4,"acceptance");
   hPhotonCuts->GetXaxis()->SetBinLabel(5,"asymmetry");
   hPhotonCuts->GetXaxis()->SetBinLabel(6,"pidprob");
   hPhotonCuts->GetXaxis()->SetBinLabel(7,"cortpcclinfo");
   hPhotonCuts->GetXaxis()->SetBinLabel(8,"PsiPair");
   hPhotonCuts->GetXaxis()->SetBinLabel(9,"CosPAngle");
   hPhotonCuts->GetXaxis()->SetBinLabel(10,"DCA R");
   hPhotonCuts->GetXaxis()->SetBinLabel(11,"DCA Z");
   hPhotonCuts->GetXaxis()->SetBinLabel(12,"Photon Quality");
   hPhotonCuts->GetXaxis()->SetBinLabel(13,"out");
   fHistograms->Add(hPhotonCuts);

   if(preCut){
      hInvMassbefore=new TH1F(Form("InvMass_before %s",GetCutNumber().Data()),"InvMass_before",1000,0,0.3);
      fHistograms->Add(hInvMassbefore);
      hArmenterosbefore=new TH2F(Form("Armenteros_before %s",GetCutNumber().Data()),"Armenteros_before",200,-1,1,1000,0,1.);
      fHistograms->Add(hArmenterosbefore);
      hEtaDistV0s = new TH1F(Form("Eta_before %s",GetCutNumber().Data()),"Eta_before",2000,-2,2);
      fHistograms->Add(hEtaDistV0s);

   }
   hInvMassafter=new TH1F(Form("InvMass_after %s",GetCutNumber().Data()),"InvMass_after",1000,0,0.3);
   fHistograms->Add(hInvMassafter);
   hArmenterosafter=new TH2F(Form("Armenteros_after %s",GetCutNumber().Data()),"Armenteros_after",200,-1,1,250,0,0.25);
   fHistograms->Add(hArmenterosafter);

   hAcceptanceCuts=new TH1F(Form("PhotonAcceptanceCuts %s",GetCutNumber().Data()),"PhotonAcceptanceCuts",10,-0.5,9.5);
   hAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(2,"maxR");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(3,"minR");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(4,"line");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(5,"maxZ");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(6,"eta");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(7,"minpt");
   hAcceptanceCuts->GetXaxis()->SetBinLabel(8,"out");
   fHistograms->Add(hAcceptanceCuts);

   // dEdx Cuts
   hdEdxCuts=new TH1F(Form("dEdxCuts %s",GetCutNumber().Data()),"dEdxCuts",10,-0.5,9.5);
   hdEdxCuts->GetXaxis()->SetBinLabel(1,"in");
   hdEdxCuts->GetXaxis()->SetBinLabel(2,"TPCelectron");
   hdEdxCuts->GetXaxis()->SetBinLabel(3,"TPCpion");
   hdEdxCuts->GetXaxis()->SetBinLabel(4,"TPCpionhighp");
   hdEdxCuts->GetXaxis()->SetBinLabel(5,"TPCkaonlowprej");
   hdEdxCuts->GetXaxis()->SetBinLabel(6,"TPCprotonlowprej");
   hdEdxCuts->GetXaxis()->SetBinLabel(7,"TPCpionlowprej");
   hdEdxCuts->GetXaxis()->SetBinLabel(8,"TOFelectron");
   hdEdxCuts->GetXaxis()->SetBinLabel(9,"TRDelectron");
   hdEdxCuts->GetXaxis()->SetBinLabel(10,"out");
   fHistograms->Add(hdEdxCuts);

   TAxis *AxisBeforedEdx = NULL;
   TAxis *AxisBeforedEdxSig = NULL;
   TAxis *AxisBeforeTOF = NULL;
   TAxis *AxisBeforeTOFSig = NULL;
   if(preCut){
      hTPCdEdxbefore=new TH2F(Form("Gamma_dEdx_before %s",GetCutNumber().Data()),"dEdx Gamma before" ,150,0.03,20,800,0,200);
      fHistograms->Add(hTPCdEdxbefore);
      AxisBeforedEdx = hTPCdEdxbefore->GetXaxis();
      hTPCdEdxSigbefore=new TH2F(Form("Gamma_dEdxSig_before %s",GetCutNumber().Data()),"dEdx Sigma Gamma before" ,150,0.03,20,400,-10,10);
      fHistograms->Add(hTPCdEdxSigbefore);
      AxisBeforedEdxSig = hTPCdEdxSigbefore->GetXaxis();

      hTOFbefore=new TH2F(Form("Gamma_TOF_before %s",GetCutNumber().Data()),"TOF Gamma before" ,150,0.03,20,11000,-1000,10000);
      fHistograms->Add(hTOFbefore);
      AxisBeforeTOF = hTOFbefore->GetXaxis();
      hTOFSigbefore=new TH2F(Form("Gamma_TOFSig_before %s",GetCutNumber().Data()),"TOF Sigma Gamma before" ,150,0.03,20,400,-6,10);
      fHistograms->Add(hTOFSigbefore);
      AxisBeforeTOFSig = hTOFSigbefore->GetXaxis();

   }
   hTPCdEdxSigafter=new TH2F(Form("Gamma_dEdxSig_after %s",GetCutNumber().Data()),"dEdx Sigma Gamma after" ,150,0.03,20,400, -10,10);
   fHistograms->Add(hTPCdEdxSigafter);

   hTPCdEdxafter=new TH2F(Form("Gamma_dEdx_after %s",GetCutNumber().Data()),"dEdx Gamma after" ,150,0.03,20,800,0,200);
   fHistograms->Add(hTPCdEdxafter);

   hTOFSigafter=new TH2F(Form("Gamma_TOFSig_after %s",GetCutNumber().Data()),"TOF Sigma Gamma after" ,150,0.03,20,400,-6,10);
   fHistograms->Add(hTOFSigafter);

   hEtaDistV0sAfterdEdxCuts = new TH1F(Form("Eta_afterdEdx %s",GetCutNumber().Data()),"Eta_afterdEdx",2000,-2,2);
   fHistograms->Add(hEtaDistV0sAfterdEdxCuts);

   hPsiPairDeltaPhiafter=new TH2F(Form("Gamma_PsiPairDeltaPhi_after %s",GetCutNumber().Data()),"Psi Pair vs Delta Phi Gamma after" ,200,-2,2,200,-2,2);
   fHistograms->Add(hPsiPairDeltaPhiafter);

   TAxis *AxisAfter = hTPCdEdxSigafter->GetXaxis();
   Int_t bins = AxisAfter->GetNbins();
   Double_t from = AxisAfter->GetXmin();
   Double_t to = AxisAfter->GetXmax();
   Double_t *newBins = new Double_t[bins+1];
   newBins[0] = from;
   Double_t factor = TMath::Power(to/from, 1./bins);
   for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
   AxisAfter->Set(bins, newBins);
   AxisAfter = hTOFSigafter->GetXaxis();
   AxisAfter->Set(bins, newBins);
   AxisAfter = hTPCdEdxafter->GetXaxis();
   AxisAfter->Set(bins, newBins);
   if(preCut){
      AxisBeforedEdx->Set(bins, newBins);
      AxisBeforeTOF->Set(bins, newBins);
      AxisBeforedEdxSig->Set(bins, newBins);
      AxisBeforeTOFSig->Set(bins, newBins);
   }
   delete [] newBins;

   hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",100,0,100);
   fHistograms->Add(hCentrality);
   hCentralityVsNumberOfPrimaryTracks=new TH2F(Form("Centrality vs Primary Tracks %s",GetCutNumber().Data()),"Centrality vs Primary Tracks ",100,0,100,4000,0,4000);
   fHistograms->Add(hCentralityVsNumberOfPrimaryTracks);

   // Event Cuts and Info
   if(preCut){
      hV0EventCuts=new TH1F(Form("ESD_EventCuts %s",GetCutNumber().Data()),"Event Cuts",7,-0.5,6.5);
      hV0EventCuts->GetXaxis()->SetBinLabel(1,"in");
      hV0EventCuts->GetXaxis()->SetBinLabel(2,"OfflineTrigger");
      hV0EventCuts->GetXaxis()->SetBinLabel(3,"nvtxcontr");
      hV0EventCuts->GetXaxis()->SetBinLabel(4,"VertexZ");
      hV0EventCuts->GetXaxis()->SetBinLabel(5,"pileup");
      hV0EventCuts->GetXaxis()->SetBinLabel(6,"centrsel");
      hV0EventCuts->GetXaxis()->SetBinLabel(7,"out");
      fHistograms->Add(hV0EventCuts);

      hVertexZ=new TH1F(Form("VertexZ %s",GetCutNumber().Data()),"VertexZ",1000,-50,50);
      fHistograms->Add(hVertexZ);

      hTriggerClass= new TH1F(Form("OfflineTrigger %s",GetCutNumber().Data()),"OfflineTrigger",35,-0.5,34.5);
      hTriggerClass->GetXaxis()->SetBinLabel( 1,"kMB");
      hTriggerClass->GetXaxis()->SetBinLabel( 2,"kINT7");
      hTriggerClass->GetXaxis()->SetBinLabel( 3,"kMUON");
      hTriggerClass->GetXaxis()->SetBinLabel( 4,"kHighMult");
      hTriggerClass->GetXaxis()->SetBinLabel( 5,"kKEMC1");
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
      fHistograms->Add(hTriggerClass);
   }
   if(!preCut){
      hTriggerClassSelected= new TH1F(Form("OfflineTriggerSelected %s",GetCutNumber().Data()),"OfflineTriggerSelected",34,-0.5,33.5);
      hTriggerClassSelected->GetXaxis()->SetBinLabel( 1,"kMB");
      hTriggerClassSelected->GetXaxis()->SetBinLabel( 2,"kINT7");
      hTriggerClassSelected->GetXaxis()->SetBinLabel( 3,"kMUON");
      hTriggerClassSelected->GetXaxis()->SetBinLabel( 4,"kHighMult");
      hTriggerClassSelected->GetXaxis()->SetBinLabel( 5,"kKEMC1");
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
      fHistograms->Add(hTriggerClassSelected);
      
      hEventPlanePhi=new TH1F(Form("EventPlaneMinusPhotonAngle %s",GetCutNumber().Data()),"EventPlaneMinusPhotonAngle",360,-TMath::Pi(),TMath::Pi());
      fHistograms->Add(hEventPlanePhi);

      
   }
   TH1::AddDirectory(kTRUE);
}

//________________________________________________________________________
Bool_t AliConversionCuts::InitPIDResponse(){
   // Set Pointer to AliPIDResponse

   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   if(man) {
      AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
      fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
      if(fPIDResponse)return kTRUE;

   }


   return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::EventIsSelected(AliVEvent *fInputEvent, AliMCEvent *fMCEvent){
   // Process Event Selection

   Int_t cutindex=0;
   if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
   cutindex++;

   // Check for MC event
   if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
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
   }

   // Event Trigger
//    cout << "before event trigger" << endl;
   if(!IsTriggerSelected(fInputEvent)){
      if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
      fEventQuality = 3;
      return kFALSE;
   }
   cutindex++;

   if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliTriggerAnalysis fTriggerAnalysis;// = new AliTriggerAnalysis;
      fHasV0AND = fTriggerAnalysis.IsOfflineTriggerFired((AliESDEvent*)fInputEvent, AliTriggerAnalysis::kV0AND);
      if(fHasV0AND&&hTriggerClass)hTriggerClass->Fill(32);
   }
//   cout << "event number " << ((AliESDEvent*)fInputEvent)->GetEventNumberInFile() << " entered"<< endl;


   // Number of Contributors Cut
   if(GetNumberOfContributorsVtx(fInputEvent)<=0) {
      if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
      fEventQuality = 5;
      return kFALSE;
   }
   cutindex++;

   // Z Vertex Position Cut
   if(!VertexZCut(fInputEvent)){
      if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
      fEventQuality = 4;
      return kFALSE;
   }
   cutindex++;

   // Pile Up Rejection

   if(fRemovePileUp){
      if(fInputEvent->IsPileupFromSPD(3,0.8,3.,2.,5.)){
         if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
         fEventQuality = 6;
         return kFALSE;
      }
   }
   cutindex++;

   // Centrality Selection
   if(!IsCentralitySelected(fInputEvent,fMCEvent)){
      if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
      fEventQuality = 1;
      return kFALSE;
   }
   cutindex++;

   // Fill Event Histograms
   if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
   if(hVertexZ)hVertexZ->Fill(fInputEvent->GetPrimaryVertex()->GetZ());
   if(hCentrality)hCentrality->Fill(GetCentrality(fInputEvent));
   if(hCentralityVsNumberOfPrimaryTracks)
      hCentralityVsNumberOfPrimaryTracks->Fill(GetCentrality(fInputEvent),
                                               ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()
                                                ->GetTask("V0ReaderV1"))->GetNumberOfPrimaryTracks());
   fEventQuality = 0;
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PhotonIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent,Bool_t checkForConvertedGamma){
   // MonteCarlo Photon Selection

   if(!mcEvent)return kFALSE;

   if (particle->GetPdgCode() == 22){


     if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
         return kFALSE;
      if(fEtaCutMin>-0.1){
         if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
            return kFALSE;
      }

      if(particle->GetMother(0) >-1 && mcEvent->Particle(particle->GetMother(0))->GetPdgCode() == 22){
         return kFALSE; // no photon as mothers!
      }

      if(particle->GetMother(0) >= mcEvent->GetNumberOfPrimaries()){
         return kFALSE; // the gamma has a mother, and it is not a primary particle
      }

      if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

      // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
      TParticle* ePos = NULL;
      TParticle* eNeg = NULL;

      if(particle->GetNDaughters() >= 2){
         for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
            TParticle *tmpDaughter = mcEvent->Particle(daughterIndex);
            if(tmpDaughter->GetUniqueID() == 5){
               if(tmpDaughter->GetPdgCode() == 11){
                  eNeg = tmpDaughter;
               } else if(tmpDaughter->GetPdgCode() == -11){
                  ePos = tmpDaughter;
               }
            }
         }
      }

      if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
         return kFALSE;
      }

      if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
         return kFALSE; // no reconstruction below the Pt cut
      }

      if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
          eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
         return kFALSE;

      if(fEtaCutMin > -0.1){
         if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
             (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
            return kFALSE;
      }

      if(ePos->R()>fMaxR){
         return kFALSE; // cuts on distance from collision point
      }

      if(TMath::Abs(ePos->Vz()) > fMaxZ){
         return kFALSE;  // outside material
      }
      if(TMath::Abs(eNeg->Vz()) > fMaxZ){
         return kFALSE;  // outside material
      }

      if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
         return kFALSE;  // line cut to exclude regions where we do not reconstruct
      } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
         return kFALSE;
      }

      if( eNeg->R() <= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
         return kFALSE; // line cut to exclude regions where we do not reconstruct
      } else if ( fEtaCutMin != -0.1 &&   eNeg->R() >= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
         return kFALSE;
      }

      return kTRUE;
      //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
   }
   return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::PhotonIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray,Bool_t checkForConvertedGamma){
   // MonteCarlo Photon Selection

   if(!aodmcArray)return kFALSE;

   if (particle->GetPdgCode() == 22){
      if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
         return kFALSE;
      if(fEtaCutMin>-0.1){
         if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
            return kFALSE;
      }

      if(particle->GetMother() > -1){
         if((static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
            return kFALSE; // no photon as mothers!
		 }
		 if(!(static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother()))->IsPrimary())){
            return kFALSE; // the gamma has a mother, and it is not a primary particle
         }
      }

      if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

      // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
      AliAODMCParticle* ePos = NULL;
      AliAODMCParticle* eNeg = NULL;

      if(particle->GetNDaughters() >= 2){
         for(Int_t daughterIndex=particle->GetDaughter(0);daughterIndex<=particle->GetDaughter(1);daughterIndex++){
            AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(aodmcArray->At(daughterIndex));
            if(!tmpDaughter) continue;
            if(((tmpDaughter->GetMCProcessCode())) == 5){    // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX
               if(tmpDaughter->GetPdgCode() == 11){
                  eNeg = tmpDaughter;
               } else if(tmpDaughter->GetPdgCode() == -11){
                  ePos = tmpDaughter;
               }
            }
         }
      }

      if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
         return kFALSE;
      }

      if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
         return kFALSE; // no reconstruction below the Pt cut
      }

      if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
          eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
         return kFALSE;

      if(fEtaCutMin > -0.1){
         if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
             (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
            return kFALSE;
      }

      Double_t rPos = sqrt( (ePos->Xv()*ePos->Xv()) + (ePos->Yv()*ePos->Yv()) );
      Double_t rNeg = sqrt( (eNeg->Xv()*eNeg->Xv()) + (eNeg->Yv()*eNeg->Yv()) );

      if(rPos>fMaxR){
         return kFALSE; // cuts on distance from collision point
      }
      if(TMath::Abs(ePos->Zv()) > fMaxZ){
         return kFALSE;  // outside material
      }
      if(TMath::Abs(eNeg->Zv()) > fMaxZ){
         return kFALSE;  // outside material
      }

      if( rPos <= ((TMath::Abs(ePos->Zv()) * fLineCutZRSlope) - fLineCutZValue)){
         return kFALSE;  // line cut to exclude regions where we do not reconstruct
      } else if ( fEtaCutMin != -0.1 &&   rPos >= ((TMath::Abs(ePos->Zv()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
         return kFALSE;
      }

      if( rNeg <= ((TMath::Abs(eNeg->Zv()) * fLineCutZRSlope) - fLineCutZValue)){
         return kFALSE; // line cut to exclude regions where we do not reconstruct
      } else if ( fEtaCutMin != -0.1 &&   rNeg >= ((TMath::Abs(eNeg->Zv()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
         return kFALSE;
      }

      return kTRUE;
      //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
   }
   return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event)
{   // Specific Photon Cuts

   Int_t cutIndex = 0;
   if(hPhotonCuts)hPhotonCuts->Fill(cutIndex);
   cutIndex++;

   // Fill Histos before Cuts
   if(hInvMassbefore)hInvMassbefore->Fill(photon->GetMass());
   if(hArmenterosbefore)hArmenterosbefore->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());

   // Gamma selection based on QT from Armenteros
   if(fDoQtGammaSelection == kTRUE){
      if(!ArmenterosQtCut(photon)){
         if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //1
         return kFALSE;
      }
   }
   cutIndex++; //2

   // Chi Cut
   if(photon->GetChi2perNDF() > fChi2CutConversion || photon->GetChi2perNDF() <=0){
      {
         if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //2
         return kFALSE;
      }
   }
   cutIndex++;//3

   // Reconstruction Acceptance Cuts
   if(!AcceptanceCuts(photon)){
      if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //3
      return kFALSE;
   }

   cutIndex++; //4
   // Asymmetry Cut
   if(fDoPhotonAsymmetryCut == kTRUE){
      if(!AsymmetryCut(photon,event)){
         if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //4
         return kFALSE;
      }
   }

   //Check the pid probability
   cutIndex++; //5
   if(!PIDProbabilityCut(photon, event)) {
      if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //5
      return kFALSE;
   }

   cutIndex++; //6
   if(!CorrectedTPCClusterCut(photon, event)) {
      if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //6
      return kFALSE;
   }

   Double_t magField = event->GetMagneticField();
   if( magField  < 0.0 ){
      magField =  1.0;
   } else {
      magField =  -1.0;
   }
   
   AliVTrack * electronCandidate = GetTrack(event,photon->GetTrackLabelNegative() );
   AliVTrack * positronCandidate = GetTrack(event,photon->GetTrackLabelPositive() );
   Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

   cutIndex++; //7
   if(!PsiPairCut(photon)) {
      if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //7
      return kFALSE;
   }

   cutIndex++; //8
   if(!CosinePAngleCut(photon, event)) {
      if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //8
      return kFALSE;
   }

   AliAODConversionPhoton* photonAOD = dynamic_cast<AliAODConversionPhoton*>(photon);
   if (photonAOD){
      photonAOD->CalculateDistanceOfClossetApproachToPrimVtx(event->GetPrimaryVertex());

      cutIndex++; //9
      if(photonAOD->GetDCArToPrimVtx() > fDCARPrimVtxCut) { //DCA R cut of photon to primary vertex
         if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //9
         return kFALSE;
      }

      cutIndex++; //10
      if(TMath::Abs(photonAOD->GetDCAzToPrimVtx()) > fDCAZPrimVtxCut) { //DCA Z cut of photon to primary vertex
         if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //10
         return kFALSE;
      }
   } else {
      cutIndex++; //9
      cutIndex++; //10
   }
   cutIndex++; //11

   if (photonAOD){
	   UChar_t photonQuality = 0;
	    AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(event);
		if(aodEvent) {
			photonQuality = DeterminePhotonQualityAOD(photonAOD, event);
		} else {
			photonQuality = photonAOD->GetPhotonQuality();
		}	
		if (fDoPhotonQualitySelectionCut && photonQuality != fPhotonQualityCut){
			if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //11
			return kFALSE;
		}	
   } 
   cutIndex++; //12
   if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //11

   // Histos after Cuts
   if(hInvMassafter)hInvMassafter->Fill(photon->GetMass());
   if(hArmenterosafter)hArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());
   if(hPsiPairDeltaPhiafter)hPsiPairDeltaPhiafter->Fill(deltaPhi,photon->GetPsiPair());
   return kTRUE;

}

///________________________________________________________________________
Bool_t AliConversionCuts::CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event)
{   //Cut on corrected TPC Cluster Info

   AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
   AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

   if(!negTrack||!posTrack)return kFALSE;

   Double_t negclsToF=0;

   if (!fUseCorrectedTPCClsInfo ){
      if(negTrack->GetTPCNclsF()!=0){
         negclsToF = (Double_t)negTrack->GetNcls(1)/(Double_t)negTrack->GetTPCNclsF();}// Ncluster/Nfindablecluster
   }
   else {
      negclsToF = negTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
   }

   Double_t posclsToF = 0.;
   if (!fUseCorrectedTPCClsInfo ){
      if(posTrack->GetTPCNclsF()!=0){
         posclsToF = (Double_t)posTrack->GetNcls(1)/(Double_t)posTrack->GetTPCNclsF();
      }
   }else{
      posclsToF = posTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
   }

   if( negclsToF < fMinClsTPCToF || posclsToF < fMinClsTPCToF ){
      return kFALSE;
   }

   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PhotonIsSelected(AliConversionPhotonBase *photon, AliVEvent * event)
{
   //Selection of Reconstructed Photons

   FillPhotonCutIndex(kPhotonIn);

   if(event->IsA()==AliESDEvent::Class()) {
     if(!SelectV0Finder( ( ((AliESDEvent*)event)->GetV0(photon->GetV0Index()))->GetOnFlyStatus() ) ){
         FillPhotonCutIndex(kOnFly);
         return kFALSE;
      }
   }
   // else if(event->IsA()==AliAODEvent::Class()) {
   //    if(!SelectV0Finder( ( ((AliAODEvent*)event)->GetV0(photon->GetV0Index())) ) ){
   //       FillPhotonCutIndex(kOnFly);
   //       return kFALSE;
   //    }
   // }

   // Get Tracks
   AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
   AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

   if(!negTrack || !posTrack) {
      FillPhotonCutIndex(kNoTracks);
      return kFALSE;
   }
   photon->DeterminePhotonQuality(negTrack,posTrack);
   // Track Cuts
   if(!TracksAreSelected(negTrack, posTrack)){
      FillPhotonCutIndex(kTrackCuts);
      return kFALSE;
   }
   if (hEtaDistV0s)hEtaDistV0s->Fill(photon->GetPhotonEta());
   // dEdx Cuts
   if(!dEdxCuts(negTrack) || !dEdxCuts(posTrack)) {
      FillPhotonCutIndex(kdEdxCuts);
      return kFALSE;
   }
   if (hEtaDistV0sAfterdEdxCuts)hEtaDistV0sAfterdEdxCuts->Fill(photon->GetPhotonEta());
   // Photon Cuts
   if(!PhotonCuts(photon,event)){
      FillPhotonCutIndex(kPhotonCuts);
      return kFALSE;
   }

   // Photon passed cuts
   FillPhotonCutIndex(kPhotonOut);
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::ArmenterosQtCut(AliConversionPhotonBase *photon)
{   // Armenteros Qt Cut
   if(fDo2DQt){
      if ( !(TMath::Power(photon->GetArmenterosAlpha()/0.95,2)+TMath::Power(photon->GetArmenterosQt()/fQtMax,2) < 1) ){
         return kFALSE;
      }
   } else {
      if(photon->GetArmenterosQt()>fQtMax){
         return kFALSE;
      }
   }
   return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionCuts::AcceptanceCuts(AliConversionPhotonBase *photon) {
   // Exclude certain areas for photon reconstruction

   Int_t cutIndex=0;
   if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
   cutIndex++;

   if(photon->GetConversionRadius()>fMaxR){ // cuts on distance from collision point
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(photon->GetConversionRadius()<fMinR){ // cuts on distance from collision point
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(photon->GetConversionRadius() <= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlope)-fLineCutZValue)){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   else if (fUseEtaMinCut &&  photon->GetConversionRadius() >= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(TMath::Abs(photon->GetConversionZ()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;


   if( photon->GetPhotonEta() > (fEtaCut)    || photon->GetPhotonEta() < (-fEtaCut) ){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   if(fEtaCutMin>-0.1){
      if( photon->GetPhotonEta() < (fEtaCutMin) && photon->GetPhotonEta() > (-fEtaCutMin) ){
         if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
         return kFALSE;
      }
   }
   cutIndex++;

   if(photon->GetPhotonPt()<fPtCut){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);

   return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionCuts::SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex) {
   // Track Cuts which require AOD/ESD specific implementation

   if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   AliAODVertex * NegVtxType=negTrack->GetProdVertex();
   AliAODVertex * PosVtxType=posTrack->GetProdVertex();
   if( (NegVtxType->GetType())==AliAODVertex::kKink || (PosVtxType->GetType())==AliAODVertex::kKink) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   return kTRUE;

}


///________________________________________________________________________
Bool_t AliConversionCuts::SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex) {
   // Track Cuts which require AOD/ESD specific implementation

   if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(negTrack->GetKinkIndex(0) > 0  || posTrack->GetKinkIndex(0) > 0 ) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   return kTRUE;
}



///________________________________________________________________________
Bool_t AliConversionCuts::TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack) {
   // Track Selection for Photon Reconstruction

   Int_t cutIndex=0;
   if(hTrackCuts)hTrackCuts->Fill(cutIndex);
   cutIndex++;

   // avoid like sign
   if(negTrack->Charge() == posTrack->Charge()) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   // Number of TPC Clusters


   if( negTrack->GetNcls(1) < fMinClsTPC || posTrack->GetNcls(1) < fMinClsTPC ) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   // Acceptance
   if( posTrack->Eta() > (fEtaCut) || posTrack->Eta() < (-fEtaCut) ||
       negTrack->Eta() > (fEtaCut) || negTrack->Eta() < (-fEtaCut) ){
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   if(fEtaCutMin>-0.1){
      if( (posTrack->Eta() < (fEtaCutMin) && posTrack->Eta() > (-fEtaCutMin)) ||
          (negTrack->Eta() < (fEtaCutMin) && negTrack->Eta() > (-fEtaCutMin)) ){
         if(hTrackCuts)hTrackCuts->Fill(cutIndex);
         return kFALSE;
      }
   }
   cutIndex++;

   // Single Pt Cut
   if( negTrack->Pt()< fSinglePtCut || posTrack->Pt()< fSinglePtCut){
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   // AOD ESD specific cuts
   Bool_t passCuts = kTRUE;

   if(negTrack->IsA()==AliAODTrack::Class()) {
      passCuts = passCuts * SpecificTrackCuts(static_cast<AliAODTrack*>(negTrack), static_cast<AliAODTrack*>(posTrack),cutIndex);
   } else {
      passCuts = passCuts * SpecificTrackCuts(static_cast<AliESDtrack*>(negTrack), static_cast<AliESDtrack*>(posTrack),cutIndex);
   }

   if(!passCuts){
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(hTrackCuts)hTrackCuts->Fill(cutIndex);

   return kTRUE;

}

///________________________________________________________________________
Bool_t AliConversionCuts::dEdxCuts(AliVTrack *fCurrentTrack){
   // Electron Identification Cuts for Photon reconstruction
   if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
   if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error

   Int_t cutIndex=0;
   if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
   if(hTPCdEdxSigbefore)hTPCdEdxSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
   if(hTPCdEdxbefore)hTPCdEdxbefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
   cutIndex++;

   if(fDodEdxSigmaCut == kTRUE){
      // TPC Electron Line
      if( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
          fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine){

         if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
         return kFALSE;
      }
      cutIndex++;

      // TPC Pion Line
      if( fCurrentTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentTrack->P()<fPIDMaxPnSigmaAbovePionLine ){
         if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
            fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
            fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
      cutIndex++;

      // High Pt Pion rej
      if( fCurrentTrack->P()>fPIDMaxPnSigmaAbovePionLine ){
         if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
            fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine &&
            fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
      cutIndex++;
   }
   else{cutIndex+=3;}

   if(fDoKaonRejectionLowP == kTRUE){
      if(fCurrentTrack->P()<fPIDMinPKaonRejectionLowP ){
         if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
   }
   cutIndex++;
   if(fDoProtonRejectionLowP == kTRUE){
      if( fCurrentTrack->P()<fPIDMinPProtonRejectionLowP ){
         if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
   }
   cutIndex++;

   if(fDoPionRejectionLowP == kTRUE){
      if( fCurrentTrack->P()<fPIDMinPPionRejectionLowP ){
         if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
   }
   cutIndex++;


   // cout<<"Start"<<endl;
   // AliPIDResponse::EDetPidStatus status=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,fCurrentTrack);

   // if( ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) && ( (status & AliVTrack::kTIME) == AliVTrack::kTIME ))
   //    {cout<<"TOF DA"<<endl;}
   // if(status == AliPIDResponse::kDetPidOk){
   //    Float_t probMis = fPIDResponse->GetTOFMismatchProbability(fCurrentTrack);
   //    cout<<"--> "<<probMis<<endl;
   //    if(probMis > 0.01){

   //    }
   // }

   if((fCurrentTrack->GetStatus() & AliESDtrack::kTOFpid) && !(fCurrentTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
      if(hTOFbefore){
         Double_t t0 = fPIDResponse->GetTOFResponse().GetStartTime(fCurrentTrack->P());
         Double_t  times[AliPID::kSPECIESC];
         fCurrentTrack->GetIntegratedTimes(times,AliPID::kSPECIESC);
         Double_t TOFsignal = fCurrentTrack->GetTOFsignal();
         Double_t dT = TOFsignal - t0 - times[0];
         hTOFbefore->Fill(fCurrentTrack->P(),dT);
      }
      if(hTOFSigbefore) hTOFSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
      if(fUseTOFpid){
         if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)>fTofPIDnSigmaAboveElectronLine ||
            fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)<fTofPIDnSigmaBelowElectronLine ){
            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
      if(hTOFSigafter)hTOFSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
   }
   cutIndex++;
   // Apply TRD PID
   if(fDoTRDPID){
      if(!fPIDResponse->IdentifiedAsElectronTRD(fCurrentTrack,fPIDTRDEfficiency)){
         if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
         return kFALSE;
      }
   }
   cutIndex++;

   if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
   if(hTPCdEdxSigafter)hTPCdEdxSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
   if(hTPCdEdxafter)hTPCdEdxafter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
   
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::AsymmetryCut(AliConversionPhotonBase * photon,AliVEvent *event) {
   // Cut on Energy Assymetry

   for(Int_t ii=0;ii<2;ii++){

      AliVTrack *track=GetTrack(event,photon->GetTrackLabel(ii));

      if( track->P() > fMinPPhotonAsymmetryCut ){
         Double_t trackNegAsy=0;
         if (photon->GetPhotonP()!=0.){
            trackNegAsy= track->P()/photon->GetPhotonP();
         }

         if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
            return kFALSE;
         }
      }

   }
   return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliConversionCuts::GetTrack(AliVEvent * event, Int_t label){
   //Returns pointer to the track with given ESD label
   //(Important for AOD implementation, since Track array in AOD data is different
   //from ESD array, but ESD tracklabels are stored in AOD Tracks)

   AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);
   if(esdEvent) {
      if(label > event->GetNumberOfTracks() ) return NULL;
      AliESDtrack * track = esdEvent->GetTrack(label);
      return track;

   } else {
      AliVTrack * track = 0x0;
      if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->AreAODsRelabeled()){
         track = dynamic_cast<AliVTrack*>(event->GetTrack(label));
         return track;
      }
      else{
         for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
            track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
            if(track){
               if(track->GetID() == label) {
                  return track;
               }
            }
         }
      }
   }
   //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
   return NULL;
}

///________________________________________________________________________
AliESDtrack *AliConversionCuts::GetESDTrack(AliESDEvent * event, Int_t label){
   //Returns pointer to the track with given ESD label
   //(Important for AOD implementation, since Track array in AOD data is different
   //from ESD array, but ESD tracklabels are stored in AOD Tracks)

   if(event) {
      if(label > event->GetNumberOfTracks() ) return NULL;
      AliESDtrack * track = event->GetTrack(label);
      return track;
   }
   //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
   return NULL;
}



///________________________________________________________________________
Bool_t AliConversionCuts::PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event){
   // Cut on Electron Probability for Photon Reconstruction

   AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);

   if(esdEvent){

      Bool_t iResult=kFALSE;

      Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
      Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

      AliESDtrack* negTrack   = esdEvent->GetTrack(photon->GetTrackLabelNegative());
      AliESDtrack* posTrack   = esdEvent->GetTrack(photon->GetTrackLabelPositive());

      if(negProbArray && posProbArray){

         negTrack->GetTPCpid(negProbArray);
         posTrack->GetTPCpid(posProbArray);

         if(negProbArray[AliPID::kElectron]>=fPIDProbabilityCutNegativeParticle && posProbArray[AliPID::kElectron]>=fPIDProbabilityCutPositiveParticle){
            iResult=kTRUE;
         }
      }

      delete [] posProbArray;
      delete [] negProbArray;
      return iResult;

   } else {
      ///Not possible for AODs
      return kTRUE;
   }




}


///________________________________________________________________________
Bool_t AliConversionCuts::AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg){
   // MC Acceptance Cuts
   //(Certain areas were excluded for photon reconstruction)

   if(particle->R()>fMaxR){
      return kFALSE;}

   if(ePos->R()>fMaxR){
      return kFALSE;
   }

   if(ePos->R()<fMinR){
      return kFALSE;
   }

   if( ePos->R() <= ((TMath::Abs(ePos->Vz())*fLineCutZRSlope)-fLineCutZValue)){
      return kFALSE;
   }
   else if (fUseEtaMinCut &&  ePos->R() >= ((TMath::Abs(ePos->Vz())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
      return kFALSE;
   }

   if(TMath::Abs(eNeg->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      return kFALSE;
   }

   if(eNeg->Vz()!=ePos->Vz()||eNeg->R()!=ePos->R()){
      return kFALSE;
   }

   if(TMath::Abs(ePos->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      return kFALSE;
   }


   if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) ){
      return kFALSE;
   }
   if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ){
      return kFALSE;
   }
   if( eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) ){
      return kFALSE;
   }
   if(fEtaCutMin>-0.1){
      if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) ){
         return kFALSE;
      }
      if( ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin) ){
         return kFALSE;
      }
      if( eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin) ){
         return kFALSE;
      }
   }

   if( ePos->Pt()< fSinglePtCut ||  eNeg->Pt()< fSinglePtCut){
      return kFALSE;
   }

   if(particle->Pt()<fPtCut){
      return kFALSE;
   }

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::UpdateCutString() {
   ///Update the cut string (if it has been created yet)

   if(fCutString && fCutString->GetString().Length() == kNCuts) {
      fCutString->SetString(GetCutNumber());
   } else {
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
void AliConversionCuts::LoadReweightingHistosMCFromFile() {

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
Bool_t AliConversionCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
   fCutStringRead = Form("%s",analysisCutSelection.Data());

   // Initialize Cuts from a given Cut string
   if(fDoReweightHistoMCPi0 || fDoReweightHistoMCEta || fDoReweightHistoMCK0s) {
      AliInfo("Weighting was enabled");
      LoadReweightingHistosMCFromFile();
   }

   AliInfo(Form("Set Photoncut Number: %s",analysisCutSelection.Data()));
   if(analysisCutSelection.Length()!=kNCuts) {
      AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
      return kFALSE;
   }
   if(!analysisCutSelection.IsAlnum()){
      AliError("Cut selection is not alphanumeric");
      return kFALSE;
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
Bool_t AliConversionCuts::SetCut(cutIds cutID, const Int_t value) {
   ///Set individual cut ID

   switch (cutID) {

   case kv0FinderType:
      if( SetV0Finder(value)) {
         fCuts[kv0FinderType] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kededxSigmaCut:
      if( SetTPCdEdxCutElectronLine(value)) {
         fCuts[kededxSigmaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kpidedxSigmaCut:
      if( SetTPCdEdxCutPionLine(value)) {
         fCuts[kpidedxSigmaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kpiMomdedxSigmaCut:
      if( SetMinMomPiondEdxCut(value)) {
         fCuts[kpiMomdedxSigmaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kchi2GammaCut:
      if( SetChi2GammaCut(value)) {
         fCuts[kchi2GammaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case ksinglePtCut:
      if( SetSinglePtCut(value)) {
         fCuts[ksinglePtCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kclsTPCCut:
      if( SetTPCClusterCut(value)) {
         fCuts[kclsTPCCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case ketaCut:
      if( SetEtaCut(value)) {
         fCuts[ketaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kLowPRejectionSigmaCut:
      if( SetLowPRejectionCuts(value)) {
         fCuts[kLowPRejectionSigmaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kQtMaxCut:
      if( SetQtMaxCut(value)) {
         fCuts[kQtMaxCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kpiMaxMomdedxSigmaCut:
      if( SetMaxMomPiondEdxCut(value)) {
         fCuts[kpiMaxMomdedxSigmaCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kRCut:
      if( SetRCut(value)) {
         fCuts[kRCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kremovePileUp:
      if( SetRemovePileUp(value)) {
         fCuts[kremovePileUp] = value;
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

   case kTOFelectronPID:
      if( SetTOFElectronPIDCut(value)) {
         fCuts[kTOFelectronPID] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kdoPhotonAsymmetryCut:
      if( SetPhotonAsymmetryCut(value)) {
         fCuts[kdoPhotonAsymmetryCut] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kPsiPair:
      if( SetPsiPairCut(value)) {
         fCuts[kPsiPair] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kCosPAngle:
      if( SetCosPAngleCut(value)) {
         fCuts[kCosPAngle] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;


   case kElecShare:
      if( SetSharedElectronCut(value)) {
         fCuts[kElecShare] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kToCloseV0s:
      if( SetToCloseV0sCut(value)) {
         fCuts[kToCloseV0s] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kExtraSignals:
      if( SetRejectExtraSignalsCut(value)) {
         fCuts[kExtraSignals] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kDcaRPrimVtx:
      if( SetDCARPhotonPrimVtxCut(value)) {
         fCuts[kDcaRPrimVtx] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kDcaZPrimVtx:
      if( SetDCAZPhotonPrimVtxCut(value)) {
         fCuts[kDcaZPrimVtx] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kInPlaneOutOfPlane:
   if( SetInPlaneOutOfPlane(value)) {
      fCuts[kInPlaneOutOfPlane] = value;
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
void AliConversionCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

void AliConversionCuts::PrintCutsWithValues() {
   // Print out current Cut Selection with value
	printf("\nConversion cutnumber \n");
	for(Int_t ic = 0; ic < kNCuts; ic++) {
		printf("%d",fCuts[ic]);
	}
	printf("\n\n");

	
   if (fIsHeavyIon == 0) {
      printf("Running in pp mode \n");
      if (fSpecialTrigger == 0){
        printf("\t only events triggered by V0OR will be analysed \n");
      } else if (fSpecialTrigger == 1){
        printf("\t only events triggered by V0AND will be analysed \n");
      } else if (fSpecialTrigger == 2){
         printf("\t only events where SDD was present will be analysed \n");
      } else if (fSpecialTrigger == 3){
         printf("\t only events where SDD was present will be analysed and triggered by VOAND\n");
      } else if (fSpecialTrigger > 3){ 
         printf("\t only events triggered by %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data());
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
      } else if (fSpecialTrigger > 4){   
         printf("\t only events triggered by %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data());
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
      } else if (fSpecialTrigger > 4){   
         printf("\t only events triggered by %s %s\n", fSpecialTriggerName.Data(), fSpecialSubTriggerName.Data());
      }
   }
   printf("MC event cuts: \n");
   if (fRejectExtraSignals == 0) printf("\t no rejection was applied \n");
	else if (fRejectExtraSignals == 1) printf("\t only MB header will be inspected \n");
	else if (fRejectExtraSignals > 1) printf("\t special header have been selected \n");
	
   printf("Electron cuts & Secondary Track Cuts - only track from secondaries enter analysis: \n");
   printf("\t no like sign pairs from V0s \n");
   if (!fUseCorrectedTPCClsInfo) printf("\t # TPC clusters > %3.2f \n", fMinClsTPC);
   if (fEtaCutMin > -0.1) printf("\t %3.2f < eta_{e} < %3.2f\n", fEtaCutMin, fEtaCut );
     else printf("\t eta_{e} < %3.2f\n", fEtaCut );
   printf("\t p_{T,e} > %3.2f\n", fSinglePtCut );
   printf("\t TPC refit \n");
   printf("\t no kinks \n");
   printf("\t accept: %3.2f < n sigma_{e,TPC} < %3.2f\n", fPIDnSigmaBelowElectronLine, fPIDnSigmaAboveElectronLine );
   printf("\t reject: %3.2f < p_{e,T} < %3.2f, n sigma_{pi,TPC} < %3.2f\n", fPIDMinPnSigmaAbovePionLine, fPIDMaxPnSigmaAbovePionLine, fPIDnSigmaAbovePionLine );
   printf("\t reject: p_{e,T} > %3.2f, n sigma_{pi,TPC} < %3.2f\n", fPIDMaxPnSigmaAbovePionLine, fPIDnSigmaAbovePionLineHighPt );
   if (fDoPionRejectionLowP) printf("\t reject: p_{e,T} < %3.2f, -%3.2f < n sigma_{pi,TPC} < %3.2f\n", fPIDMinPPionRejectionLowP, fPIDnSigmaAtLowPAroundPionLine, fPIDnSigmaAtLowPAroundPionLine );
   if (fDoKaonRejectionLowP) printf("\t reject: -%3.2f < n sigma_{K,TPC} < %3.2f\n", fPIDnSigmaAtLowPAroundKaonLine, fPIDnSigmaAtLowPAroundKaonLine );
   if (fDoProtonRejectionLowP) printf("\t reject: -%3.2f < n sigma_{K,TPC} < %3.2f\n", fPIDnSigmaAtLowPAroundProtonLine, fPIDnSigmaAtLowPAroundProtonLine );
   if (fUseTOFpid) printf("\t accept: %3.2f < n sigma_{e,TOF} < %3.2f\n", fTofPIDnSigmaBelowElectronLine, fTofPIDnSigmaAboveElectronLine);
   
   printf("Photon cuts: \n");
   if (fUseOnFlyV0Finder) printf("\t using Onfly V0 finder \n");
   else printf("\t using Offline V0 finder \n");
   if (fDo2DQt){
	  printf("\t 2 dimensional q_{T} cut applied with maximum of %3.2f \n", fQtMax );
   } else {
	  printf("\t 1 dimensional q_{T} cut applied with maximum of %3.2f \n", fQtMax );
   }
   if (fDo2DPsiPairChi2){
	  printf("\t 2 dimensional triangle chi^{2} and psi_{pair} cut applied with maximum of chi^{2} = %3.2f and |psi_{pair}| = %3.2f \n", fChi2CutConversion, fPsiPairCut ); 
   } else {
      printf("\t chi^{2} max cut chi^{2} < %3.2f \n", fChi2CutConversion ); 
	  printf("\t psi_{pair} max cut |psi_{pair}| < %3.2f \n", fPsiPairCut ); 
   }	      
   printf("\t %3.2f < R_{conv} < %3.2f\n", fMinR, fMaxR );
   printf("\t Z_{conv} < %3.2f\n", fMaxZ );
   if (fEtaCutMin > -0.1) printf("\t %3.2f < eta_{conv} < %3.2f\n", fEtaCutMin, fEtaCut );
     else printf("\t eta_{conv} < %3.2f\n", fEtaCut );
   if (fDoPhotonAsymmetryCut) printf("\t for p_{T,track} > %3.2f,  A_{gamma} < %3.2f \n", fMinPPhotonAsymmetryCut, fMinPhotonAsymmetry  );
   if (fUseCorrectedTPCClsInfo) printf("\t #cluster TPC/ #findable clusters TPC (corrected for radius) > %3.2f\n", fMinClsTPCToF );
   printf("\t p_{T,gamma} > %3.2f\n", fPtCut );	 
   printf("\t cos(Theta_{point}) > %3.2f \n", fCosPAngleCut );
   printf("\t dca_{R} < %3.2f \n", fDCARPrimVtxCut );
   printf("\t dca_{Z} < %3.2f \n", fDCAZPrimVtxCut );
   if (fDoPhotonQualitySelectionCut) printf("\t selection based on photon quality with quality %d \n", fPhotonQualityCut );
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetIsHeavyIon(Int_t isHeavyIon)
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
Bool_t AliConversionCuts::SetCentralityMin(Int_t minCentrality)
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
Bool_t AliConversionCuts::SetCentralityMax(Int_t maxCentrality)
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
Bool_t AliConversionCuts::SetSelectSpecialTrigger(Int_t selectSpecialTrigger)
{// Set Cut

   switch(selectSpecialTrigger){
   case 0:
      fSpecialTrigger=0; // dont care
      break;
   case 1:
      fSpecialTrigger=1; // V0AND
      break;
   case 2:
      fSpecialTrigger=2; // with SDD requested
      break;
   case 3:
      fSpecialTrigger=3; // V0AND plus with SDD requested
      break;
   // allows to run MB & 6 other different trigger classes in parallel with the same photon cut
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
	  fSpecialTriggerName="AliVEvent::kEMCEGA";
      break;
    case 9:
      fSpecialTrigger=9; // trigger alias kEMCEJE
      fOfflineTriggerMask=AliVEvent::kEMCEJE;
      fTriggerSelectedManually = kTRUE;
	  fSpecialTriggerName="AliVEvent::kEMCEJE";
      break;
   default:
      AliError("Warning: Special Trigger Not known");
      return 0;
   }
   return 1;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetSelectSubTriggerClass(Int_t selectSpecialSubTriggerClass)
{// Set Cut

	if (fSpecialTrigger == 1){ //V0AND with different detectors
		switch(selectSpecialSubTriggerClass){
		case 0: //with VZERO
			fSpecialTrigger=1;
			fSpecialSubTrigger=0; 
// 			AliInfo("Info: Nothing to be done");
			break;
		case 1: //with TZERO
			fSpecialTrigger=0;
			fSpecialSubTrigger=0; 
			fOfflineTriggerMask=AliVEvent::kINT8;
			fTriggerSelectedManually = kTRUE;
			fSpecialTriggerName="AliVEvent::kINT8";
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
// 			AliInfo("Info: Nothing to be done");
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
// 			AliInfo("Info: Nothing to be done");
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
		default:
			AliError("Warning: Special Subtrigger Class Not known");
			return 0;
		}		   
	}  else if (fSpecialTrigger == 6){ // Subdivision of kPHI trigger classes
		switch(selectSpecialSubTriggerClass){
		case 0: // all together
			fSpecialSubTrigger=0; 
			fSpecialSubTriggerName="";
// 			AliInfo("Info: Nothing to be done");
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
// 			AliInfo("Info: Nothing to be done");
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
		default:
			AliError("Warning: Special Subtrigger Class Not known");
			return 0;
		}		   
	}  else if (fSpecialTrigger == 8){ // Subdivision of kEMCEGA trigger classes
		switch(selectSpecialSubTriggerClass){
		case 0: // all together
			fSpecialSubTrigger=0; 
			fSpecialSubTriggerName="";
// 			AliInfo("Info: Nothing to be done");
			break;
		case 1: // 7EGA - CINT7 EGA
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="7EGA";
			break;
		case 2: // 8EGA - CINT8 EGA
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="8EGA";
			break;
		case 3: // 7EG1 - CINT7 EG1
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="7EG1";
			break;
		case 4: // 8EG1 - CINT8 EG1
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="8EG1";
			break;
		case 5: // 7EG2 - CINT7 EG2
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="7EG2";
			break;
		case 6: // 8EG2 - CINT8 EG2
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="8EG2";
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
// 			AliInfo("Info: Nothing to be done");
			break;
		case 1: // 7EJE - CINT7 EJE
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="7EJE";
			break;
		case 2: // 8EJE - CINT8 EJE
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="8EJE";
			break;
		case 3: // 7EJ1 - CINT7 EJ1
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="7EJ1";
			break;
		case 4: // 8EJ1 - CINT8 EJ1
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="8EJ1";
			break;
		case 5: // 7EJ2 - CINT7 EJ2
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="7EJ2";
			break;
		case 6: // 8EJ2 - CINT8 EJ2
			fSpecialSubTrigger=1; 
			fNSpecialSubTriggerOptions=1;
			fSpecialSubTriggerName="8EJ2";
			break;
		default:
			AliError("Warning: Special Subtrigger Class Not known");
			return 0;
		}		   
	}
	return 1;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetMultiplicityMethod(Int_t multiplicityMethod)
{
   // Set Cut
   fMultiplicityMethod=multiplicityMethod;

   // 0 Photon Multiplicity
   // 1 TPC Track multiplicity
   // 2 V0 Mult
   // 3 SPD Mult

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetRemovePileUp(Int_t removePileUp)
{// Set Cut
   switch(removePileUp){
   case 0:
      fRemovePileUp=kFALSE;
      break;
   case 1:
      fRemovePileUp=kTRUE;
      break;
   default:
      AliError("RemovePileUpCut not defined");
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetRejectExtraSignalsCut(Int_t extraSignal) {

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
///________________________________________________________________________
Bool_t AliConversionCuts::SetV0Finder(Int_t v0FinderType)
{   // Set Cut
   switch (v0FinderType){
   case 0:  // on fly V0 finder
      cout << "have chosen onfly V0" << endl;
      fUseOnFlyV0Finder=kTRUE;
      break;
   case 1:  // offline V0 finder
      cout << "have chosen offline V0" << endl;
      fUseOnFlyV0Finder=kFALSE;
      break;
   default:
      AliError(Form(" v0FinderType not defined %d",v0FinderType));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetEtaCut(Int_t etaCut)
{   // Set Cut

   //Set Standard LineCutZValues
   fLineCutZValueMin = -2;
   fLineCutZValue = 7.;

   switch(etaCut){
   case 0: // 0.9
      fEtaCut     = 0.9;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 1:  // 0.6  // changed from 1.2 to 0.6 on 2013.06.10
      fEtaCut     = 0.6;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 2:  // 1.4
      fEtaCut     = 1.4;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 3: // 0.65
      fEtaCut     = 0.65;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 4: // 0.75
      fEtaCut     = 0.75;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 5: // 0.5
      fEtaCut     = 0.5;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 6: // 5.
      fEtaCut     = 5.;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 7:
      if (fIsHeavyIon==1){
         fEtaCut     = 0.7;
         fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
         fEtaCutMin     = -0.1;
         fLineCutZRSlopeMin = 0.;
         break;
      } else {   
         fEtaCut     = 0.3;
         fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
         fEtaCutMin     = -0.1;
         fLineCutZRSlopeMin = 0.;
         break;
      }
   // case 8: // 0.1 - 0.8
   //    fEtaCut     = 0.9;
   //    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
   //    fEtaCutMin     = 0.1;
   //    fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
   //    break;
   case 8: // 0.4
      fEtaCut     = 0.4;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 9: // 10
      fEtaCut     = 10;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   default:
      AliError(Form(" EtaCut not defined %d",etaCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetRCut(Int_t RCut){
   // Set Cut
   switch(RCut){
   case 0:
      fMinR=0;
      fMaxR = 180.;
      break;
   case 1:
      fMinR=2.8;
      fMaxR = 180.;
      break;
   case 2:
      fMinR=5.;
      fMaxR = 180.;
      break;
   case 3:
      fMaxR = 70.;
      fMinR = 10.;
      break;
   case 4:
      fMaxR = 70.;
      fMinR = 5.;
      break;
   case 5:
      fMaxR = 180.;
      fMinR = 10.;
      break;
   case 6:
      fMaxR = 180.;
      fMinR = 20.;
      break;
   case 7:
      fMaxR = 180.;
      fMinR = 35.; //old 26.
      break;
   case 8:
      fMaxR = 180.;
      fMinR = 12.5;
      break;
   case 9:
      fMaxR = 180.;
      fMinR = 7.5;
      break;

   default:
      AliError("RCut not defined");
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetSinglePtCut(Int_t singlePtCut)
{   // Set Cut
   switch(singlePtCut){
   case 0: // 0.050 GeV
      fSinglePtCut = 0.050;
      break;
   case 1:  // 0.100 GeV
      fSinglePtCut = 0.100;
      break;
   case 2:  // 0.150 GeV
      fSinglePtCut = 0.150;
      break;
   case 3:  // 0.200 GeV
      fSinglePtCut = 0.200;
      break;
   case 4:  // 0.075 GeV
      fSinglePtCut = 0.075;
      break;
   case 5:  // 0.125 GeV
      fSinglePtCut = 0.125;
      break;
   case 6:  // 0.04 GeV
      fSinglePtCut = 0.040;
      break;
   case 7:  // 0.0 GeV
      fSinglePtCut = 0.0;
      break;
   default:
      AliError(Form("singlePtCut not defined %d",singlePtCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetTPCClusterCut(Int_t clsTPCCut)
{   // Set Cut
   switch(clsTPCCut){
   case 0: // 0
      fMinClsTPC= 0.;
      break;
   case 1:  // 60
      fMinClsTPC= 60.;
      break;
   case 2:  // 80
      fMinClsTPC= 80.;
      break;
   case 3:  // 100
      fMinClsTPC= 100.;
      break;
   case 4:  // 95% of findable clusters
      fMinClsTPCToF= 0.95;
      fUseCorrectedTPCClsInfo=1;
      break;
   case 5:  // 0% of findable clusters
      fMinClsTPCToF= 0.0;
      fUseCorrectedTPCClsInfo=1;
      break;
   case 6:  // 70% of findable clusters
      fMinClsTPCToF= 0.7;
      fUseCorrectedTPCClsInfo=1;
      break;
   case 7:  // 0% of findable clusters
      fMinClsTPCToF= 0.35;
      fUseCorrectedTPCClsInfo=0;
      break;
   case 8:
      fMinClsTPCToF= 0.35;
      fUseCorrectedTPCClsInfo=1;
      break;
   case 9:
      fMinClsTPCToF= 0.6;
      fUseCorrectedTPCClsInfo=1;
      break;
   default:
      AliError(Form("Warning: clsTPCCut not defined %d",clsTPCCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut)
{   // Set Cut
   switch(ededxSigmaCut){
   case 0: // -10,10
      fPIDnSigmaBelowElectronLine=-10;
      fPIDnSigmaAboveElectronLine=10;
      break;
   case 1: // -5,5
      fPIDnSigmaBelowElectronLine=-5;
      fPIDnSigmaAboveElectronLine=5;
      break;
   case 2: // -3,5
      fPIDnSigmaBelowElectronLine=-3;
      fPIDnSigmaAboveElectronLine=5;
      break;
   case 3: // -4,5
      fPIDnSigmaBelowElectronLine=-4;
      fPIDnSigmaAboveElectronLine=5;
      break;
   case 4: // -6,7
      fPIDnSigmaBelowElectronLine=-6;
      fPIDnSigmaAboveElectronLine=7;
      break;
   case 5: // -4,4
      fPIDnSigmaBelowElectronLine=-4;
      fPIDnSigmaAboveElectronLine=4;
      break;
   case 6: // -2.5,4
      fPIDnSigmaBelowElectronLine=-2.5;
      fPIDnSigmaAboveElectronLine=4;
      break;
   case 7: // -2,3.5
      fPIDnSigmaBelowElectronLine=-2;
      fPIDnSigmaAboveElectronLine=3.5;
      break;
   default:
      AliError("TPCdEdxCutElectronLine not defined");
      return kFALSE;

   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut)
{   // Set Cut

   switch(pidedxSigmaCut){
   case 0:  // -10
      fPIDnSigmaAbovePionLine=-10;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 1:   // 0
      fPIDnSigmaAbovePionLine=0;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 2:  // 1
      fPIDnSigmaAbovePionLine=1;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 3:  // 1
      fPIDnSigmaAbovePionLine=2.5;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 4:  // 1
      fPIDnSigmaAbovePionLine=0.5;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 5:  // 1
      fPIDnSigmaAbovePionLine=2.;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 6:  // 1
      fPIDnSigmaAbovePionLine=2.;
      fPIDnSigmaAbovePionLineHighPt=0.5;
      break;
   case 7:  // 1
      fPIDnSigmaAbovePionLine=3.5;
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   case 8:  // 1
      fPIDnSigmaAbovePionLine=2.;
      fPIDnSigmaAbovePionLineHighPt=1.;
      break;
   case 9:
      fPIDnSigmaAbovePionLine=3.0; // We need a bit less tight cut on dE/dx
      fPIDnSigmaAbovePionLineHighPt=-10;
      break;
   default:
      AliError(Form("Warning: pidedxSigmaCut not defined %d",pidedxSigmaCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetMinMomPiondEdxCut(Int_t piMomdedxSigmaCut)
{   // Set Cut
   switch(piMomdedxSigmaCut){
   case 0:  // 0.5 GeV
      fPIDMinPnSigmaAbovePionLine=0.5;
      break;
   case 1:  // 1. GeV
      fPIDMinPnSigmaAbovePionLine=1.;
      break;
   case 2:  // 1.5 GeV
      fPIDMinPnSigmaAbovePionLine=1.5;
      break;
   case 3:  // 20.0 GeV
      fPIDMinPnSigmaAbovePionLine=20.;
      break;
   case 4:  // 50.0 GeV
      fPIDMinPnSigmaAbovePionLine=50.;
      break;
   case 5:  // 0.3 GeV
      fPIDMinPnSigmaAbovePionLine=0.3;
      break;
   case 6:  // 0.25 GeV
      fPIDMinPnSigmaAbovePionLine=0.25;
      break;
   case 7:  // 0.4 GeV
      fPIDMinPnSigmaAbovePionLine=0.4;
      break;
   case 8:  // 0.2 GeV
      fPIDMinPnSigmaAbovePionLine=0.2;
      break;
   default:
      AliError(Form("piMomdedxSigmaCut not defined %d",piMomdedxSigmaCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut)
{   // Set Cut
   switch(piMaxMomdedxSigmaCut){
   case 0:  // 100. GeV
      fPIDMaxPnSigmaAbovePionLine=100.;
      break;
   case 1:  // 5. GeV
      fPIDMaxPnSigmaAbovePionLine=5.;
      break;
   case 2:  // 4. GeV
      fPIDMaxPnSigmaAbovePionLine=4.;
      break;
   case 3:  // 3.5 GeV
      fPIDMaxPnSigmaAbovePionLine=3.5;
      break;
   case 4:  // 3. GeV
      fPIDMaxPnSigmaAbovePionLine=3.;
      break;
   case 5:  // 7. GeV
      fPIDMaxPnSigmaAbovePionLine=7.;
      break;
   default:
      AliError(Form("piMaxMomdedxSigmaCut not defined %d",piMaxMomdedxSigmaCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut)
{   // Set Cut
   switch(LowPRejectionSigmaCut){
   case 0:  //
      fPIDnSigmaAtLowPAroundKaonLine=0;
      fPIDnSigmaAtLowPAroundProtonLine=0;
      fPIDnSigmaAtLowPAroundPionLine=0;
      fDoKaonRejectionLowP = kFALSE;
      fDoProtonRejectionLowP = kFALSE;
      fDoPionRejectionLowP = kFALSE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 1:  //
      fPIDnSigmaAtLowPAroundKaonLine=0.5;
      fPIDnSigmaAtLowPAroundProtonLine=0.5;
      fPIDnSigmaAtLowPAroundPionLine=0.5;
      fDoKaonRejectionLowP = kTRUE;
      fDoProtonRejectionLowP = kTRUE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 2:  //
      fPIDnSigmaAtLowPAroundKaonLine=1;
      fPIDnSigmaAtLowPAroundProtonLine=1;
      fPIDnSigmaAtLowPAroundPionLine=1;
      fDoKaonRejectionLowP = kTRUE;
      fDoProtonRejectionLowP = kTRUE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 3:  //
      fPIDnSigmaAtLowPAroundKaonLine=2.;
      fPIDnSigmaAtLowPAroundProtonLine=2.;
      fPIDnSigmaAtLowPAroundPionLine=2.;
      fDoKaonRejectionLowP = kTRUE;
      fDoProtonRejectionLowP = kTRUE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 4:  //
      fPIDnSigmaAtLowPAroundKaonLine=0.;
      fPIDnSigmaAtLowPAroundProtonLine=0.;
      fPIDnSigmaAtLowPAroundPionLine=1;
      fDoKaonRejectionLowP = kFALSE;
      fDoProtonRejectionLowP = kFALSE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 5:  //
      fPIDnSigmaAtLowPAroundKaonLine=0.;
      fPIDnSigmaAtLowPAroundProtonLine=0.;
      fPIDnSigmaAtLowPAroundPionLine=1.5;
      fDoKaonRejectionLowP = kFALSE;
      fDoProtonRejectionLowP = kFALSE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 6:  //
      fPIDnSigmaAtLowPAroundKaonLine=0.;
      fPIDnSigmaAtLowPAroundProtonLine=0.;
      fPIDnSigmaAtLowPAroundPionLine=2.;
      fDoKaonRejectionLowP = kFALSE;
      fDoProtonRejectionLowP = kFALSE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   case 7:  //
      fPIDnSigmaAtLowPAroundKaonLine=0.;
      fPIDnSigmaAtLowPAroundProtonLine=0.;
      fPIDnSigmaAtLowPAroundPionLine=0.5;
      fDoKaonRejectionLowP = kFALSE;
      fDoProtonRejectionLowP = kFALSE;
      fDoPionRejectionLowP = kTRUE;
      fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
      break;
   default:
      AliError(Form("LowPRejectionSigmaCut not defined %d",LowPRejectionSigmaCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetTOFElectronPIDCut(Int_t TOFelectronPID){
   // Set Cut
   switch(TOFelectronPID){
   case 0: // no cut
      fUseTOFpid = kFALSE;
      fTofPIDnSigmaBelowElectronLine=-100;
      fTofPIDnSigmaAboveElectronLine=100;
      break;
   case 1: // -7,7
      fUseTOFpid = kTRUE;
      fTofPIDnSigmaBelowElectronLine=-7;
      fTofPIDnSigmaAboveElectronLine=7;
      break;
   case 2: // -5,5
      fUseTOFpid = kTRUE;
      fTofPIDnSigmaBelowElectronLine=-5;
      fTofPIDnSigmaAboveElectronLine=5;
      break;
   case 3: // -3,5
      fUseTOFpid = kTRUE;
      fTofPIDnSigmaBelowElectronLine=-3;
      fTofPIDnSigmaAboveElectronLine=5;
      break;
   case 4: // -2,3
      fUseTOFpid = kTRUE;
      fTofPIDnSigmaBelowElectronLine=-2;
      fTofPIDnSigmaAboveElectronLine=3;
      break;
   case 5: // -3,3
      fUseTOFpid = kTRUE;
      fTofPIDnSigmaBelowElectronLine=-3;
      fTofPIDnSigmaAboveElectronLine=3;
      break;
   default:
      AliError(Form("TOFElectronCut not defined %d",TOFelectronPID));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetQtMaxCut(Int_t QtMaxCut)
{   // Set Cut
   switch(QtMaxCut){
   case 0: //
      fQtMax=1.;
      fDoQtGammaSelection=kFALSE;
      fDo2DQt=kFALSE;
      break;
   case 1:
      fQtMax=0.1;
      fDo2DQt=kFALSE;
      break;
   case 2:
      fQtMax=0.07;
      fDo2DQt=kFALSE;
      break;
   case 3:
      fQtMax=0.05;
      fDo2DQt=kFALSE;
      break;
   case 4:
      fQtMax=0.03;
      fDo2DQt=kFALSE;
      break;
   case 5:
      fQtMax=0.02;
      fDo2DQt=kFALSE;
      break;
   case 6:
      fQtMax=0.02;
      fDo2DQt=kTRUE;
      break;
   case 7:
      fQtMax=0.15;
      fDo2DQt=kFALSE;
      break;
   case 8:
      fQtMax=0.05;
      fDo2DQt=kTRUE;
      break;   
   case 9:
      fQtMax=0.03;
      fDo2DQt=kTRUE;
      break;      
   default:
      AliError(Form("Warning: QtMaxCut not defined %d",QtMaxCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetChi2GammaCut(Int_t chi2GammaCut)
{   // Set Cut

   switch(chi2GammaCut){
   case 0: // 100
      fChi2CutConversion = 100.;
      break;
   case 1:  // 50
      fChi2CutConversion = 50.;
      break;
   case 2:  // 30
      fChi2CutConversion = 30.;
      break;
   case 3:
      fChi2CutConversion = 200.;
      break;
   case 4:
      fChi2CutConversion = 500.;
      break;
   case 5:
      fChi2CutConversion = 100000.;
      break;
   case 6:
      fChi2CutConversion = 5.;
      break;
   case 7:
      fChi2CutConversion = 10.;
      break;
   case 8:
      fChi2CutConversion = 20.;
      break;
   case 9:
      fChi2CutConversion = 15.;
      break;
   default:
      AliError(Form("Warning: Chi2GammaCut not defined %d",chi2GammaCut));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetPsiPairCut(Int_t psiCut) {

   switch(psiCut) {
   case 0:
      fPsiPairCut = 10000; //
      break;
   case 1:
      fPsiPairCut = 0.1; //
      break;
   case 2:
      fPsiPairCut = 0.05; // Standard
      break;
   case 3:
      fPsiPairCut = 0.035; //
      break;
   case 4:
      fPsiPairCut = 0.2; //
      break;   
   case 5:
      fPsiPairCut = 0.1; //
      fDo2DPsiPairChi2 = kTRUE;
      break;
   case 6:
      fPsiPairCut = 0.05; //
      fDo2DPsiPairChi2 = kTRUE;
      break;
   case 7:
      fPsiPairCut = 0.035; //
      fDo2DPsiPairChi2 = kTRUE;
      break;
   case 8:
      fPsiPairCut = 0.2; //
      fDo2DPsiPairChi2 = kTRUE; //
      break;
   case 9:
      fPsiPairCut = 0.5; //
      break;
   default:
      AliError(Form("PsiPairCut not defined %d",psiCut));
      return kFALSE;
   }

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut){
   // Set Cut
   switch(doPhotonAsymmetryCut){
   case 0:
      fDoPhotonAsymmetryCut=0;
      fMinPPhotonAsymmetryCut=100.;
      fMinPhotonAsymmetry=0.;
      break;
   case 1:
      fDoPhotonAsymmetryCut=1;
      fMinPPhotonAsymmetryCut=3.5;
      fMinPhotonAsymmetry=0.04;
      break;
   case 2:
      fDoPhotonAsymmetryCut=1;
      fMinPPhotonAsymmetryCut=3.5;
      fMinPhotonAsymmetry=0.06;
      break;
   case 3:
      fDoPhotonAsymmetryCut=1;
      fMinPPhotonAsymmetryCut=0.0;
      fMinPhotonAsymmetry=0.05;
      break; 
   default:
      AliError(Form("PhotonAsymmetryCut not defined %d",doPhotonAsymmetryCut));
      return kFALSE;
   }
   fCuts[kdoPhotonAsymmetryCut]=doPhotonAsymmetryCut;
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetCosPAngleCut(Int_t cosCut) {

   switch(cosCut){
   case 0:
      fCosPAngleCut = -1; 
      break;
   case 1:
      fCosPAngleCut = 0; 
      break;
   case 2:
      fCosPAngleCut = 0.5; 
      break;
   case 3:
      fCosPAngleCut = 0.75; 
      break;
   case 4:
      fCosPAngleCut = 0.85; 
      break;
   case 5:
      fCosPAngleCut = 0.88; 
      break;
   case 6:
      fCosPAngleCut = 0.9;
      break;
   case 7:
      fCosPAngleCut = 0.95;
      break;
   default:
      AliError(Form("Cosine Pointing Angle cut not defined %d",cosCut));
      return kFALSE;
   }

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetSharedElectronCut(Int_t sharedElec) {

	switch(sharedElec){
	case 0:
		fDoSharedElecCut = kFALSE;
		fDoPhotonQualitySelectionCut = kFALSE;
		fPhotonQualityCut = 0;
		break;
	case 1:
		fDoSharedElecCut = kTRUE;
		fDoPhotonQualitySelectionCut = kFALSE;
		fPhotonQualityCut = 0;
		break;
	case 2:
		fDoSharedElecCut = kFALSE;
		fDoPhotonQualitySelectionCut = kTRUE;
		fPhotonQualityCut = 1;
		break;
	case 3:
		fDoSharedElecCut = kFALSE;
		fDoPhotonQualitySelectionCut = kTRUE;	  
		fPhotonQualityCut = 2;
		break;
	case 4:
		fDoSharedElecCut = kFALSE;
		fDoPhotonQualitySelectionCut = kTRUE;	  
		fPhotonQualityCut = 3;
		break;
	default:
		AliError(Form("Shared Electron Cut not defined %d",sharedElec));	
		return kFALSE;
   }

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetToCloseV0sCut(Int_t toClose) {

   switch(toClose){
   case 0:
      fDoToCloseV0sCut = kFALSE;
      fminV0Dist = 250;
      break;
   case 1:
      fDoToCloseV0sCut = kTRUE;
      fminV0Dist = 1;
      break;
   case 2:
      fDoToCloseV0sCut = kTRUE;
      fminV0Dist = 2;
      break;
   case 3:
      fDoToCloseV0sCut = kTRUE;
      fminV0Dist = 3;
      break;
   default:
      AliError(Form("Shared Electron Cut not defined %d",toClose));
      return kFALSE;
   }
   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::SetTRDElectronCut(Int_t TRDElectronCut)
{   // Set Cut
   switch(TRDElectronCut){
   case 0:
      fDoTRDPID=kFALSE;
      break;
   case 1:
      fDoTRDPID=kTRUE;
      fPIDTRDEfficiency=0.1;
      break;
   case 8:
      fDoTRDPID=kTRUE;
      fPIDTRDEfficiency=0.8;
      break;
   case 9:
      fDoTRDPID=kTRUE;
      fPIDTRDEfficiency=0.9;
      break;
   default:
      AliError(Form("TRDElectronCut not defined %d",TRDElectronCut));
      return kFALSE;
   }

   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetDCAZPhotonPrimVtxCut(Int_t DCAZPhotonPrimVtx){
   // Set Cut
   switch(DCAZPhotonPrimVtx){
   case 0:  //
      fDCAZPrimVtxCut   = 1000;
      break;
   case 1:  //
      fDCAZPrimVtxCut   = 10;
      break;
   case 2:  //
      fDCAZPrimVtxCut   = 5;
      break;
   case 3:  //
      fDCAZPrimVtxCut   = 4;
      break;
   case 4:  //
      fDCAZPrimVtxCut   = 3;
      break;
   case 5:  //
      fDCAZPrimVtxCut   = 2.5;
      break;
   case 6:  //
      fDCAZPrimVtxCut   = 2;
      break;
   case 7:  //
      fDCAZPrimVtxCut   = 1.5;
      break;
   case 8:  //
      fDCAZPrimVtxCut   = 1;
      break;
   case 9:  //
      fDCAZPrimVtxCut   = 0.5;
      break;
   default:
      cout<<"Warning: DCAZPhotonPrimVtx not defined "<<DCAZPhotonPrimVtx<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetDCARPhotonPrimVtxCut(Int_t DCARPhotonPrimVtx){
   // Set Cut
   switch(DCARPhotonPrimVtx){
   case 0:  //
      fDCARPrimVtxCut   = 1000;
      break;
   case 1:  //
      fDCARPrimVtxCut   = 10;
      break;
   case 2:  //
      fDCARPrimVtxCut   = 5;
      break;
   case 3:  //
      fDCARPrimVtxCut   = 4;
      break;
   case 4:  //
      fDCARPrimVtxCut   = 3;
      break;
   case 5:  //
      fDCARPrimVtxCut   = 2.5;
      break;
   case 6:  //
      fDCARPrimVtxCut   = 2;
      break;
   case 7:  //
      fDCARPrimVtxCut   = 1.5;
      break;
   case 8:  //
      fDCARPrimVtxCut   = 1;
      break;
   case 9:  //
      fDCARPrimVtxCut   = 0.5;
      break;
   default:
      cout<<"Warning: DCARPhotonPrimVtx not defined "<<DCARPhotonPrimVtx<<endl;
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::SetInPlaneOutOfPlane(Int_t inOutPlane){
   // Set Cut
   switch(inOutPlane){
   case 0:  //
      fInPlaneOutOfPlane = 0; // No Event Plane
      break;
   case 1:  //
      fInPlaneOutOfPlane = 1; // In-Plane
      break;
   case 2:  //
      fInPlaneOutOfPlane = 2; // Out-Of-Plane
      break;
   default:
      cout<<"Warning: In-Plane or Out-Of-Plane not defined "<<inOutPlane<<endl;
      return kFALSE;
   }
   return kTRUE;
}

//-------------------------------------------------------------
Bool_t AliConversionCuts::GetUseNewMultiplicityFramework(TString period){ 

  if((period.CompareTo("LHC15o")==0)||(period.CompareTo("LHC15k1_plus3")==0)||(period.CompareTo("LHC15k1a1")==0)
     ||(period.CompareTo("LHC15k1b1")==0)||(period.CompareTo("LHC15k1a2")==0)||(period.CompareTo("LHC15k1a3")==0)
     ||(period.CompareTo("LHC15k1b2")==0)||(period.CompareTo("LHC15k1b3")==0)||(period.CompareTo("LHC15k1_plus")==0)
     ||(period.CompareTo("LHC15k1")==0)||(period.CompareTo("LHC15k1_plus2")==0)||(period.CompareTo("LHC15k1_plus21")==0)
     ||(period.CompareTo("LHC15k1_plus31")==0)||(period.CompareTo("LHC15k1_plus32")==0)){ 
     return kTRUE;
  }
  else{
     return kFALSE;
  } 
}

//-------------------------------------------------------------
Double_t AliConversionCuts::GetCentrality(AliVEvent *event)
{   // Get Event Centrality

   AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
   if(esdEvent){
     TString periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->GetPeriodName();
     if(GetUseNewMultiplicityFramework(periodName)){
       AliMultSelection *MultSelection = (AliMultSelection*)event->FindListObject("MultSelection");
       AliCentrality *fESDCentrality = (AliCentrality*)esdEvent->GetCentrality();
       if(fDetectorCentrality==0){
         if(fIsHeavyIon==2)             return fESDCentrality->GetCentralityPercentile("V0A"); // default for pPb
         else                           return MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
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
     if(aodEvent->GetHeader()){return ((AliVAODHeader*)aodEvent->GetHeader())->GetCentrality();}
   }

   return -1;
}
//-------------------------------------------------------------
Bool_t AliConversionCuts::IsCentralitySelected(AliVEvent *event, AliMCEvent *fMCEvent)
{   // Centrality Selection
   if(!fIsHeavyIon)return kTRUE;

   if(fCentralityMin == fCentralityMax ) return kTRUE;//0-100%
   else if(fCentralityMax==0) fCentralityMax=10; //CentralityRange = fCentralityMin-100%

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

   Int_t nprimaryTracks = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->GetNumberOfPrimaryTracks();
	Int_t PrimaryTracks10[11][2] =
		{
			{9999,9999}, //  0
			{1210, 928}, // 10
			{ 817, 658}, // 20
			{ 536, 435}, // 30
			{ 337, 276}, // 40
			{ 197, 162}, // 50
			{ 106, 100}, // 60
			{  51,  44}, // 70
			{  21,  18}, // 80
			{   0,   0},  // 90
			{   0,   0}  // 100 // only max accessible
		};
	Int_t PrimaryTracks5a[11][2] =
		{
			{9999,9999}, // 0
			{1485,1168}, // 5
			{1210, 928}, // 10
			{ 995, 795}, // 15
			{ 817, 658}, // 20
			{ 666, 538}, // 25
			{ 536, 435}, // 30
			{ 428, 350}, // 35
			{ 337, 276}, // 40
			{ 260, 214},  // 45
			{ 0, 162}  // 50 only max accessible
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
			{   0,   0}  // 100 only max accessible
		};
   Int_t column = 0;
   if(event->IsA()==AliESDEvent::Class()) column = 0;
   if(event->IsA()==AliAODEvent::Class()) column = 1;

   if (fModCentralityClass == 3){
      if(fMCEvent){
         if(nprimaryTracks > PrimaryTracks10[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks10[fCentralityMin][column])
            return kTRUE;
         else return kFALSE;
      }
      else{
         centralityC= Int_t(centrality/10);
         if(centralityC >= fCentralityMin && centralityC < fCentralityMax)
            return kTRUE;
         else return kFALSE;
      }
   }
   else if (fModCentralityClass ==4){
      if(fMCEvent){
         if(nprimaryTracks > PrimaryTracks5a[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks5a[fCentralityMin][column])
            return kTRUE;
         else return kFALSE;
      }
      else{
         centralityC= Int_t(centrality);
         if(centralityC >= fCentralityMin*5 && centralityC < fCentralityMax*5){
            return kTRUE;
         } else return kFALSE;
      }
   }
   else if (fModCentralityClass ==5){
      if(fMCEvent){
         if(nprimaryTracks > PrimaryTracks5b[fCentralityMax][column] && nprimaryTracks <= PrimaryTracks5b[fCentralityMin][column])
            return kTRUE;
         else return kFALSE;
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
///________________________________________________________________________
Bool_t AliConversionCuts::VertexZCut(AliVEvent *event){
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

   TString periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()
                                                ->GetTask("V0ReaderV1"))->GetPeriodName();
   if (periodName.CompareTo("LHC11h")==0){
     if (TMath::Abs(fVertexZ-fVertexZSPD) > 0.1) return kFALSE;
   }						
   if (fIsHeavyIon == 2){
     if(fUtils->IsFirstEventInChunk(event)) return kFALSE;
     if(!fUtils->IsVertexSelected2013pA(event)) return kFALSE;
     if(fUtils->IsPileUpEvent(event)) return kFALSE;
   }

   return kTRUE;
}
///________________________________________________________________________

Int_t AliConversionCuts::GetNumberOfContributorsVtx(AliVEvent *event){
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
         }  else {
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

///________________________________________________________________________

Bool_t AliConversionCuts::IsTriggerSelected(AliVEvent *fInputEvent)
{

	AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
// 	AliTRDTriggerAnalysis *trdSelection=  new AliTRDTriggerAnalysis();
// 	trdSelection->CalcTriggers(fInputEvent);
	
	UInt_t isSelected = AliVEvent::kAny;
	TString periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->GetPeriodName();
	//    cout << 	periodName.Data() << endl;
	
	if (fInputHandler==NULL) return kFALSE;
	if( fInputHandler->GetEventSelection() || fInputEvent->IsA()==AliAODEvent::Class()) {
	
		TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();  
		if (!fTriggerSelectedManually){
			if (fPreSelCut) fOfflineTriggerMask = AliVEvent::kAny;
			else {
				if (fIsHeavyIon == 1) fOfflineTriggerMask = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
				else if (fIsHeavyIon == 2) fOfflineTriggerMask = AliVEvent::kINT7;
				else if (periodName.CompareTo("LHC11c") == 0 || periodName.CompareTo("LHC11d") == 0 || periodName.CompareTo("LHC11e") == 0 || periodName.CompareTo("LHC11f") == 0 || periodName.CompareTo("LHC11g") == 0  || periodName.CompareTo("LHC12a") == 0 || periodName.CompareTo("LHC12b") == 0 || periodName.CompareTo("LHC12c") == 0 || periodName.CompareTo("LHC12d") == 0 || periodName.CompareTo("LHC12f") == 0  || periodName.CompareTo("LHC12g") == 0  || periodName.CompareTo("LHC12h") == 0  || periodName.CompareTo("LHC12i") == 0  ||periodName.CompareTo("LHC13g") == 0 ) {
					fOfflineTriggerMask = AliVEvent::kINT7;      
	// 				cout << "will take kINT7 as trigger mask" << endl; 
				}	
				else fOfflineTriggerMask = AliVEvent::kMB;
			}
		}
		// Get the actual offline trigger mask for the event and AND it with the
		// requested mask. If no mask requested select by default the event.
	//       if (fPreSelCut) cout << "Trigger selected from outside: "<< fTriggerSelectedManually <<"\t Offline Trigger mask for Precut: " << fOfflineTriggerMask << endl;
	//       else cout << "Trigger selected from outside: "<< fTriggerSelectedManually <<"\t Offline Trigger mask: " << fOfflineTriggerMask << endl;

		if (fOfflineTriggerMask){
			isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();		 
			if (isSelected && !fPreSelCut){
// 				if (fSpecialTriggerName.Contains("kTRD")){ // make special selection for TRD trigger
// 					Bool_t bTRDHasFiredConfirmed= 0; // bool whether the TRD has triggered and has been read out due to that trigger & should have triggered
// 					Bool_t bTRDClassContainedInTriggerList= 1; //check whether the trigger list contains the requested trigger
// 					if (fSpecialSubTrigger>0){
// 						if (fSpecialSubTriggerName.Contains("HSE")){
// 							bTRDHasFiredConfirmed = trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHSE);
// 						} else if (fSpecialSubTriggerName.Contains("HJT")){
// 							bTRDHasFiredConfirmed = trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHJT);
// 						} else if (fSpecialSubTriggerName.Contains("HEE")){
// 							bTRDHasFiredConfirmed = trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHEE);
// 						} else if (fSpecialSubTriggerName.Contains("HQU")){	
// 							bTRDHasFiredConfirmed = trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHQU);
// 						}
// 						if (!firedTrigClass.Contains(fSpecialSubTriggerName.Data())) bTRDClassContainedInTriggerList = 0;
// 					} else {
// 						bTRDHasFiredConfirmed = trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHQU)  || trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHSE) || trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHJT) || trdSelection->HasTriggeredConfirmed(AliTRDTriggerAnalysis::kHEE);
// 					}	
// 					if (!bTRDHasFiredConfirmed || !bTRDClassContainedInTriggerList)  isSelected = 0;									
// 				} else { // more general condition for all other sub-triggers
// 					if (fSpecialTriggerName.Contains("kEMCE"))cout << fSpecialTriggerName.Data() << "\t" <<fSpecialSubTriggerName.Data()<< endl;
					if (fSpecialSubTrigger>0){
						if (!firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 0;
					}	 
// 					if (fSpecialTriggerName.Contains("kEMCE"))cout <<firedTrigClass << endl;
// 				}
			}				
		} 	 
	}
	fIsSDDFired = !(fInputHandler->IsEventSelected() & AliVEvent::kFastOnly);

	// Fill Histogram
	if(hTriggerClass){
		if (fIsSDDFired) hTriggerClass->Fill(33);
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

	if(hTriggerClassSelected && isSelected){
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

	if(!isSelected)return kFALSE;

	return kTRUE;

}

///________________________________________________________________________
Int_t AliConversionCuts::GetFirstTPCRow(Double_t radius){
   // Get first TPC row
   Int_t firstTPCRow = 0;
   Double_t radiusI = 84.8;
   Double_t radiusO = 134.6;
   Double_t radiusOB = 198.;
   Double_t rSizeI = 0.75;
   Double_t rSizeO = 1.;
   Double_t rSizeOB = 1.5;
   Int_t nClsI = 63;
   Int_t nClsIO = 127;

   if(radius <= radiusI){
      return firstTPCRow;
   }
   if(radius>radiusI && radius<=radiusO){
      firstTPCRow = (Int_t)((radius-radiusI)/rSizeI);
   }
   if(radius>radiusO && radius<=radiusOB){
      firstTPCRow = (Int_t)(nClsI+(radius-radiusO)/rSizeO);
   }

   if(radius>radiusOB){
      firstTPCRow =(Int_t)(nClsIO+(radius-radiusOB)/rSizeOB);
   }

   return firstTPCRow;
}

Bool_t AliConversionCuts::CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const {
   ///Check if passes cosine of pointing angle cut
   if(GetCosineOfPointingAngle(photon, event) < fCosPAngleCut){
      return kFALSE;
   }
   return kTRUE;
}

Double_t AliConversionCuts::GetCosineOfPointingAngle( const AliConversionPhotonBase * photon, AliVEvent * event) const{
   // calculates the pointing angle of the recalculated V0

   Double_t momV0[3] = {0,0,0};
   if(event->IsA()==AliESDEvent::Class()){
      AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
      if(!esdEvent) return -999;
      AliESDv0 *v0 = esdEvent->GetV0(photon->GetV0Index());
      if(!v0) return -999;
      v0->GetPxPyPz(momV0[0],momV0[1],momV0[2]);
   }
   if(event->IsA()==AliAODEvent::Class()){
      momV0[0] = photon->GetPx();
      momV0[1] = photon->GetPy();
      momV0[2] = photon->GetPz();
   }

   //Double_t momV0[3] = { photon->GetPx(), photon->GetPy(), photon->GetPz() }; //momentum of the V0
   Double_t PosV0[3] = { photon->GetConversionX() - event->GetPrimaryVertex()->GetX(),
                         photon->GetConversionY() - event->GetPrimaryVertex()->GetY(),
                         photon->GetConversionZ() - event->GetPrimaryVertex()->GetZ() }; //Recalculated V0 Position vector

   Double_t momV02 = momV0[0]*momV0[0] + momV0[1]*momV0[1] + momV0[2]*momV0[2];
   Double_t PosV02 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1] + PosV0[2]*PosV0[2];


   Double_t cosinePointingAngle = -999;
   if(momV02*PosV02 > 0.0)
      cosinePointingAngle = (PosV0[0]*momV0[0] +  PosV0[1]*momV0[1] + PosV0[2]*momV0[2] ) / TMath::Sqrt(momV02 * PosV02);

   return cosinePointingAngle;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PsiPairCut(const AliConversionPhotonBase * photon) const {

   if (fDo2DPsiPairChi2){
      if (TMath::Abs(photon->GetPsiPair()) < -fPsiPairCut/fChi2CutConversion*photon->GetChi2perNDF() + fPsiPairCut ){
         return kTRUE;
      } else {
         return kFALSE;
      }    
   } else {
      if(TMath::Abs(photon->GetPsiPair()) > fPsiPairCut){
         return kFALSE;}
      else{return kTRUE;}
   } 
}

///________________________________________________________________________
TString AliConversionCuts::GetCutNumber(){
   // returns TString with current cut number
   return fCutStringRead;
}

///________________________________________________________________________
void AliConversionCuts::FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0){

   Int_t posLabel = photon->GetTrackLabelPositive();
   Int_t negLabel = photon->GetTrackLabelNegative();

   fElectronLabelArray[nV0*2] = posLabel;
   fElectronLabelArray[(nV0*2)+1] = negLabel;
}
///________________________________________________________________________
Bool_t AliConversionCuts::RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s){

   Int_t posLabel = photon->GetTrackLabelPositive();
   Int_t negLabel = photon->GetTrackLabelNegative();

   for(Int_t i = 0; i<nV0s*2;i++){
      if(i==nV0*2)     continue;
      if(i==(nV0*2)+1) continue;
      if(fElectronLabelArray[i] == posLabel){
         return kFALSE;}
      if(fElectronLabelArray[i] == negLabel){
         return kFALSE;}
   }

   return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionCuts::RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0){


   Double_t posX = photon->GetConversionX();
   Double_t posY = photon->GetConversionY();
   Double_t posZ = photon->GetConversionZ();

   for(Int_t i = 0;i<photons->GetEntries();i++){
      if(nV0 == i) continue;
      AliAODConversionPhoton *photonComp = (AliAODConversionPhoton*) photons->At(i);
      Double_t posCompX = photonComp->GetConversionX();
      Double_t posCompY = photonComp->GetConversionY();
      Double_t posCompZ = photonComp->GetConversionZ();

      Double_t dist = pow((posX - posCompX),2)+pow((posY - posCompY),2)+pow((posZ - posCompZ),2);

      if(dist < fminV0Dist*fminV0Dist){
         if(photon->GetChi2perNDF() < photonComp->GetChi2perNDF()) return kTRUE;
         else {
            return kFALSE;}
      }

   }
   return kTRUE;
}
///________________________________________________________________________
void AliConversionCuts::GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *MCEvent){

	TString periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->GetPeriodName();


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

	if(rejection == 0) return; // No Rejection

	AliGenCocktailEventHeader *cHeader = 0x0;
	AliAODMCHeader *cHeaderAOD = 0x0;
	Bool_t headerFound = kFALSE;
    AliMCEvent *mcEvent = 0x0;
    TClonesArray *mcEventAOD = 0x0;
	if(MCEvent->IsA()==AliMCEvent::Class()){
		if (dynamic_cast<AliMCEvent*>(MCEvent)){
			cHeader = dynamic_cast<AliGenCocktailEventHeader*>(dynamic_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
			if(cHeader) headerFound = kTRUE;
		}	
	}
	if(MCEvent->IsA()==AliAODEvent::Class()){ // MCEvent is a AODEvent in case of AOD
		cHeaderAOD = dynamic_cast<AliAODMCHeader*>(MCEvent->FindListObject(AliAODMCHeader::StdBranchName()));
        mcEventAOD = dynamic_cast<TClonesArray*>(MCEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		
		
		if(cHeaderAOD) headerFound = kTRUE;
	}

	if(headerFound){
		TList *genHeaders = 0x0;
		if(cHeader) genHeaders = cHeader->GetHeaders();
		if(cHeaderAOD){
			genHeaders = cHeaderAOD->GetCocktailHeaders();
			if(genHeaders->GetEntries()==1){
				SetRejectExtraSignalsCut(0);
				return;
			}
		}
		AliGenEventHeader* gh = 0;
		fnHeaders = 0;
		Int_t firstindexA = 0;
		Int_t lastindexA =  -1;
		if(rejection == 1 || rejection == 3) fnHeaders = 1; // MinBiasHeader
		if(rejection == 2){ // TList of Headers Names
			for(Int_t i = 0; i<genHeaders->GetEntries();i++){
				gh = (AliGenEventHeader*)genHeaders->At(i);
				TString GeneratorName = gh->GetName();
				lastindexA = lastindexA + gh->NProduced();
// 				cout << i << "\t" << GeneratorName.Data() << endl;
				for(Int_t j = 0; j<HeaderList->GetEntries();j++){
					TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
					if(GeneratorName.CompareTo(GeneratorInList) == 0){
						if (GeneratorInList.CompareTo("PARAM") == 0 || GeneratorInList.CompareTo("BOX") == 0 ){
                            if(mcEvent){
                                if (mcEvent->Particle(firstindexA)->GetPdgCode() == fAddedSignalPDGCode ) {
									if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
                                        if (gh->NProduced() > 10 && mcEvent->Particle(firstindexA+10)->GetPdgCode() == fAddedSignalPDGCode ){
// 											cout << "cond 1: "<< fnHeaders << endl;
											fnHeaders++;
											continue;
										}	
										continue;
									} else {
// 										cout << "cond 2: " << fnHeaders << endl;
										fnHeaders++;
										continue;
									}	
								}
							}   
                            if ( mcEventAOD){
                                AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(mcEventAOD->At(firstindexA));
								if (  aodMCParticle->GetPdgCode() == fAddedSignalPDGCode ){
									if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
										if (gh->NProduced() > 10){
                                            AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(mcEventAOD->At(firstindexA+10));
											if (  aodMCParticle2->GetPdgCode() == fAddedSignalPDGCode ){
// 												cout << "cond 1: " << fnHeaders << endl;
												fnHeaders++;
												continue;
											} 
										}	
										continue;
									} else {
// 										cout << "cond 2: " << fnHeaders << endl;
										fnHeaders++;
										continue;
									}	
								}   
							}
							continue;
						}
// 						cout << "cond 3: "<< fnHeaders << endl;
						fnHeaders++;
						continue;
					}
				}
				firstindexA = firstindexA + gh->NProduced();
			}
		}
// 		cout << "number of headers: " <<fnHeaders << endl;
		
		fNotRejectedStart = new Int_t[fnHeaders];
		fNotRejectedEnd = new Int_t[fnHeaders];
		fGeneratorNames = new TString[fnHeaders];

		if(rejection == 1 || rejection == 3){
			fNotRejectedStart[0] = 0;
			fNotRejectedEnd[0] = ((AliGenEventHeader*)genHeaders->At(0))->NProduced()-1;
			fGeneratorNames[0] = ((AliGenEventHeader*)genHeaders->At(0))->GetName();
			return;
		}

		Int_t firstindex = 0;
		Int_t lastindex =  -1;
		Int_t number = 0;
		
		for(Int_t i = 0; i<genHeaders->GetEntries();i++){
			gh = (AliGenEventHeader*)genHeaders->At(i);
			TString GeneratorName = gh->GetName();
			lastindex = lastindex + gh->NProduced();
			for(Int_t j = 0; j<HeaderList->GetEntries();j++){
				TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
// 				cout << i << "\t" << GeneratorName.Data() << endl;
				if(GeneratorName.CompareTo(GeneratorInList) == 0){
					if (GeneratorInList.CompareTo("PARAM") == 0 || GeneratorInList.CompareTo("BOX") == 0 ){
                        if(mcEvent){
                            if (mcEvent->Particle(firstindex)->GetPdgCode() == fAddedSignalPDGCode ) {
								if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
// 									cout << "produced " << gh->NProduced() << " with box generator" << endl;
                                    if (gh->NProduced() > 10 && mcEvent->Particle(firstindex+10)->GetPdgCode() == fAddedSignalPDGCode){
// 										cout << "one of them was a pi0 or eta" <<  endl;
										fNotRejectedStart[number] = firstindex;
										fNotRejectedEnd[number] = lastindex;
										fGeneratorNames[number] = GeneratorName;
										number++;
// 										cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
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
						}   
                        if ( mcEventAOD){
                            AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(mcEventAOD->At(firstindex));
							if (  aodMCParticle->GetPdgCode() == fAddedSignalPDGCode ){
								if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
									if (gh->NProduced() > 10) {
                                        AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(mcEventAOD->At(firstindex+10));
										if ( aodMCParticle2->GetPdgCode() == fAddedSignalPDGCode ){
											fNotRejectedEnd[number] = lastindex;
											fNotRejectedStart[number] = firstindex;
											fGeneratorNames[number] = GeneratorName;
											number++;
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
						}
						continue;
					} else {
						fNotRejectedStart[number] = firstindex;
						fNotRejectedEnd[number] = lastindex;
						fGeneratorNames[number] = GeneratorName;
// 						cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
						number++;
						continue;
					}
					
				}
			}
			firstindex = firstindex + gh->NProduced();
		}
// 		for (Int_t i = 0; i < number; i++){
// 			cout << i << "\t" <<fGeneratorNames[i] << "\t" << fNotRejectedStart[i] << "\t" <<fNotRejectedEnd[i] << endl;
// 		}	
	
	} else { // No Cocktail Header Found
		fNotRejectedStart = new Int_t[1];
		fNotRejectedEnd = new Int_t[1];

		fnHeaders = 1;
		fNotRejectedStart[0] = 0;
        fNotRejectedEnd[0] = mcEvent->GetNumberOfPrimaries()-1;
		fGeneratorNames = new TString[1];
		fGeneratorNames[0] = "NoCocktailGeneratorFound";

		AliGenPythiaEventHeader *mcHeaderPythia = dynamic_cast<AliGenPythiaEventHeader*>(static_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
		if (mcHeaderPythia) fGeneratorNames[0] = "NoCocktailGeneratorFound_Pythia";
		AliGenDPMjetEventHeader *mcHeaderPhojet = dynamic_cast<AliGenDPMjetEventHeader*>(static_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
		if (mcHeaderPhojet) fGeneratorNames[0] = "NoCocktailGeneratorFound_Phojet";
		AliGenHijingEventHeader *mcHeaderHijing = dynamic_cast<AliGenHijingEventHeader*>(static_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
		if (mcHeaderHijing) fGeneratorNames[0] = "NoCocktailGeneratorFound_Hijing";

		SetRejectExtraSignalsCut(0);
	}
	
}

//_________________________________________________________________________
Int_t AliConversionCuts::IsParticleFromBGEvent(Int_t index, AliMCEvent *mcEvent, AliVEvent *InputEvent){

	// Not Accepted == kFALSE == 0
	//     Accepted ==  kTRUE == 1
	//  FirstHeader ==  kTRUE == 3
	if(index < 0) return 0; // No Particle

// 	if (index == 100){
// 		cout << "possible headers" << endl;
// 		for(Int_t i = 0;i<fnHeaders;i++){
// 			cout << i << "\t" <<fGeneratorNames[i] << "\t" << fNotRejectedStart[i] << "\t" <<fNotRejectedEnd[i] << endl;
// 		}
// 	}	
	Int_t accepted = 0;
	if(!InputEvent || InputEvent->IsA()==AliESDEvent::Class()){
        if( index >= mcEvent->GetNumberOfPrimaries()){ // Secondary Particle
            if( ((TParticle*)mcEvent->Particle(index))->GetMother(0) < 0) return 1; // Secondary Particle without Mother??
            return IsParticleFromBGEvent(((TParticle*)mcEvent->Particle(index))->GetMother(0),mcEvent,InputEvent);
		}
		for(Int_t i = 0;i<fnHeaders;i++){
			if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
				accepted = 1;
				if(i == 0) accepted = 2; // MB Header
			}
		}
	}
	else if(InputEvent->IsA()==AliAODEvent::Class()){
		TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(InputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		if (AODMCTrackArray){
			AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(index));
			if(!aodMCParticle) return 1; // Photon Without a Mother ? --> Accepted
			if(!aodMCParticle->IsPrimary()){
				if( aodMCParticle->GetMother() < 0) return 1;// Secondary Particle without Mother??
                return IsParticleFromBGEvent(aodMCParticle->GetMother(),mcEvent,InputEvent);
			}
			index = abs(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(index))->GetLabel());
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
Int_t AliConversionCuts::IsEventAcceptedByConversionCut(AliConversionCuts *ReaderCuts, AliVEvent *InputEvent, AliMCEvent *MCEvent, Int_t isHeavyIon){

   if ( !IsTriggerSelected(InputEvent) )
      return 3;

   if(isHeavyIon != 0 && !(IsCentralitySelected(InputEvent,MCEvent)))
      return 1; // Check Centrality --> Not Accepted => eventQuality = 1
      
      
   if(isHeavyIon == 0 && GetIsFromPileup()){
      if(InputEvent->IsPileupFromSPD(3,0.8,3.,2.,5.)){

         return 6; // Check Pileup --> Not Accepted => eventQuality = 6
      }
   }

   Bool_t hasV0And = ReaderCuts->HasV0AND();
   Bool_t isSDDFired = ReaderCuts->IsSDDFired();
   if( (IsSpecialTrigger() == 2 || IsSpecialTrigger() == 3) && !isSDDFired && !MCEvent)
      return 7; // With SDD requested but no fired

   if( (IsSpecialTrigger() == 1 || IsSpecialTrigger() == 3) && !hasV0And)
      return 8; // V0AND requested but no fired

   if(hCentrality)hCentrality->Fill(GetCentrality(InputEvent));
   if(hCentralityVsNumberOfPrimaryTracks)
      hCentralityVsNumberOfPrimaryTracks->Fill(GetCentrality(InputEvent),
                                             ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()
                                                ->GetTask("V0ReaderV1"))->GetNumberOfPrimaryTracks());     

   return 0;
}

//_________________________________________________________________________
Float_t AliConversionCuts::GetWeightForMeson(TString period, Int_t index, AliMCEvent *mcEvent, AliVEvent *InputEvent){
	if (!(period.CompareTo("LHC12f1a") == 0 || period.CompareTo("LHC12f1b") == 0  || period.CompareTo("LHC12i3") == 0 || period.CompareTo("LHC11a10a") == 0 || period.CompareTo("LHC11a10b") == 0 || period.CompareTo("LHC11a10b_bis") == 0 || period.CompareTo("LHC11a10a_bis") == 0 || period.CompareTo("LHC11a10b_plus") == 0 || period.Contains("LHC13d2")|| period.Contains("LHC14a1") || 
	period.CompareTo("LHC13e7") == 0 || period.Contains("LHC13b2_efix") || period.CompareTo("LHC14b2") == 0 )) return 1.;

	Int_t kCaseGen = 0;
	for (Int_t i = 0; i < fnHeaders; i++){
		if (index >= fNotRejectedStart[i] && index < fNotRejectedEnd[i]+1){
			if (fGeneratorNames[i].CompareTo("Pythia") == 0){
				kCaseGen = 1;
			} else if (fGeneratorNames[i].CompareTo("DPMJET") == 0){
				kCaseGen = 2;
			} else if (fGeneratorNames[i].CompareTo("HIJING") == 0 ||
						fGeneratorNames[i].CompareTo("Hijing") == 0 ||
						fGeneratorNames[i].Contains("hijing")){
				kCaseGen = 3;
			} else if (fGeneratorNames[i].CompareTo("BOX") == 0){
				kCaseGen = 4;
			} else if (fGeneratorNames[i].CompareTo("PARAM") == 0){
				kCaseGen = 5;
			} else if (fGeneratorNames[i].CompareTo("NoCocktailGeneratorFound") == 0){
				kCaseGen = 6;
			} else if (fGeneratorNames[i].CompareTo("NoCocktailGeneratorFound_Pythia") == 0){
				kCaseGen = 1;
			} else if (fGeneratorNames[i].CompareTo("NoCocktailGeneratorFound_Phojet") == 0){
				kCaseGen = 2;
			} else if (fGeneratorNames[i].CompareTo("NoCocktailGeneratorFound_Hijing") == 0){
				kCaseGen = 3;
			}
			if (period.Contains("LHC13d2") || period.CompareTo("LHC13e7") == 0 || period.Contains("LHC13b2_efix")  || period.Contains("LHC14a1") || period.CompareTo("LHC14b2") == 0 ){
				kCaseGen = 3;
			}
		}
	}
	if (kCaseGen == 0) return 1;


	Double_t mesonPt = 0;
	Double_t mesonMass = 0;
	Int_t PDGCode = 0;
	if(!InputEvent || InputEvent->IsA()==AliESDEvent::Class()){
        mesonPt = ((TParticle*)mcEvent->Particle(index))->Pt();
        mesonMass = ((TParticle*)mcEvent->Particle(index))->GetCalcMass();
        PDGCode = ((TParticle*)mcEvent->Particle(index))->GetPdgCode();
	}
	else if(InputEvent->IsA()==AliAODEvent::Class()){
		TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(InputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
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
	if (kCaseGen == 1){ // Pythia 6
		Float_t dNdyMC = 2.1462;
		Float_t nMC = 7.06055;
		Float_t tMC = 0.12533;
		if ( PDGCode ==  111){
			dNdyMC = 2.1462;
			nMC = 7.06055;
			tMC = 0.12533;
		} else if ( PDGCode ==  221){
			dNdyMC = 0.2357;
			nMC = 5.9105;
			tMC = 0.1525;
		}
		functionResultMC = dNdyMC / ( 2 * TMath::Pi())*(nMC-1.)*(nMC-2.) / (nMC*tMC*(nMC*tMC+mesonMass*(nMC-2.)))  * TMath::Power(1.+(TMath::Sqrt(mesonPt*mesonPt+mesonMass*mesonMass)-mesonMass)/(nMC*tMC), -nMC);
	} else if (kCaseGen == 2){ // Phojet
		Float_t dNdyMC = 2.35978;
		Float_t nMC = 6.81795;
		Float_t tMC = 0.11492;
		if ( PDGCode ==  111){
			dNdyMC = 2.35978;
			nMC = 6.81795;
			tMC = 0.11492;
		} else if ( PDGCode ==  221){
			dNdyMC = 0.3690;
			nMC = 5.55809;
			tMC = 0.13387;
		}
		functionResultMC = dNdyMC / ( 2 * TMath::Pi())*(nMC-1.)*(nMC-2.) / (nMC*tMC*(nMC*tMC+mesonMass*(nMC-2.)))  * TMath::Power(1.+(TMath::Sqrt(mesonPt*mesonPt+mesonMass*mesonMass)-mesonMass)/(nMC*tMC), -nMC);
	} else if (kCaseGen == 4){ // BOX generators pp
	//       functionResultMC = 1./sqrt(1.-mesonMass*mesonMass/((mesonMass*mesonMass+mesonPt*mesonPt)*cosh(mesonY)*cosh(mesonY)));
		Float_t a = 0.23437;
		Float_t b = 5.6661;
		Float_t c = -1430.5863;
		Float_t d = -0.6966624;
		Float_t e = 252.3742;
		if ( PDGCode ==  111){
			a = 0.23437;
			b = 5.6661;
			c = -1430.5863;
			d = -0.6966624;
			e = 252.3742;
		} else if ( PDGCode ==  221){
			a = 0.10399;
			b = 4.35311;
			c = -12.17723;
			d = -0.01172;
			e =1.85140;
		}
		functionResultMC = a*TMath::Power(mesonPt,-1.*(b+c/(TMath::Power(mesonPt,d)+e)))*1./mesonPt *1./1.6 *1./(2.* TMath::Pi());
	//       cout << functionResultMC << endl;
	} else if (kCaseGen == 3 ){ // HIJING
		if ( PDGCode ==  111 && fDoReweightHistoMCPi0 && hReweightMCHistPi0!= 0x0){
			functionResultMC = hReweightMCHistPi0->Interpolate(mesonPt);
		}
		if ( PDGCode ==  221 && fDoReweightHistoMCEta && hReweightMCHistEta!= 0x0){
			functionResultMC = hReweightMCHistEta->Interpolate(mesonPt);
		}
		if ( PDGCode ==  310 && fDoReweightHistoMCK0s && hReweightMCHistK0s!= 0x0){
			functionResultMC = hReweightMCHistK0s->Interpolate(mesonPt);
		}
	}

	Float_t functionResultData = 1;
	if (kCaseGen == 1 || kCaseGen == 2 || kCaseGen == 4 ){
		Float_t dNdyData = 2.2328;
		Float_t nData = 7.1473;
		Float_t tData = 0.1346;
		if ( PDGCode ==  111){
			dNdyData = 2.2328;
			nData = 7.1473;
			tData = 0.1346;
		} else if ( PDGCode ==  221){
			dNdyData = 0.38992; //be careful this fit is not optimal, eta in data still has problems
			nData = 5.72778;
			tData = 0.13835;
		}
		functionResultData = dNdyData / ( 2 * TMath::Pi())*(nData-1.)*(nData-2.) / (nData*tData*(nData*tData+mesonMass*(nData-2.)))  * TMath::Power(1.+(TMath::Sqrt(mesonPt*mesonPt+mesonMass*mesonMass)-mesonMass)/(nData*tData), -nData);
	//       cout << functionResultData << endl;
	} else {
		if ( PDGCode ==  111 && fDoReweightHistoMCPi0 && fFitDataPi0!= 0x0){
			functionResultData = fFitDataPi0->Eval(mesonPt);
		}
		if ( PDGCode ==  221 && fDoReweightHistoMCEta && fFitDataEta!= 0x0){
			functionResultData = fFitDataEta->Eval(mesonPt);
		}
		if ( PDGCode ==  310 && fDoReweightHistoMCK0s && fFitDataK0s!= 0x0){
			functionResultData = fFitDataK0s->Eval(mesonPt);
		}

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

	//    if (fModCentralityClass == 0 && fCentralityMin == 4 && fCentralityMax == 6 && PDGCode ==  111){
	//        cout << period.Data() << "\t" << kCaseGen << "\t" <<fModCentralityClass<< "\t" <<fCentralityMin<< "\t" <<fCentralityMax << "\t" << mesonPt << "\t" <<mesonMass<< "\t"<<functionResultData << "\t"<< functionResultMC << "\t" << weight <<endl;
	//    }
	return weight;
}
///________________________________________________________________________
AliConversionCuts* AliConversionCuts::GetStandardCuts2010PbPb(){
    //Create and return standard 2010 PbPb cuts
    AliConversionCuts *cuts=new AliConversionCuts("StandardCuts2010PbPb","StandardCuts2010PbPb");
    if(!cuts->InitializeCutsFromCutString("100000204209297002322000000")){
   cout<<"Warning: Initialization of Standardcuts2010PbPb failed"<<endl;}
    return cuts;
}

///________________________________________________________________________
AliConversionCuts* AliConversionCuts::GetStandardCuts2010pp(){
    //Create and return standard 2010 PbPb cuts
    AliConversionCuts *cuts=new AliConversionCuts("StandardCuts2010pp","StandardCuts2010pp");
    if(!cuts->InitializeCutsFromCutString("000001100209366300380000000")){
   cout<<"Warning: Initialization of Standardcuts2010pp failed"<<endl;}
    return cuts;
}
///________________________________________________________________________
void AliConversionCuts::GetCorrectEtaShiftFromPeriod(TString periodName){

   if(periodName.CompareTo("LHC12g") == 0 || //pilot run 2012
      periodName.CompareTo("LHC13b") == 0 || //mainly minimum bias
      periodName.CompareTo("LHC13c") == 0 || //mainly minimum bias
      periodName.CompareTo("LHC13d") == 0 || //mainly triggered
      periodName.CompareTo("LHC13e") == 0 || //mainly triggered
      periodName.CompareTo("LHC13c3") == 0 || //MC Starlight, anchor LHC13d+e
      periodName.CompareTo("LHC13c2") == 0 || //MC Starlight, coherent J/Psi, UPC muon anchor LHC13d+e
      periodName.CompareTo("LHC13b4") == 0 || //MC Pythia 6 (Jet-Jet), anchor LHC13b
      periodName.CompareTo("LHC13b2_fix_1") == 0 || //MC DPMJET, anchr LHC13b+c
      periodName.CompareTo("LHC13b2_efix_p1") == 0 || //MC DPMJET, anchr LHC13b+c
      periodName.CompareTo("LHC13b2_efix_p2") == 0 || //MC DPMJET, anchr LHC13b+c
      periodName.CompareTo("LHC13b2_efix_p3") == 0 || //MC DPMJET, anchr LHC13b+c
      periodName.CompareTo("LHC13b2_efix_p4") == 0 || //MC DPMJET, anchr LHC13b+c
      periodName.CompareTo("LHC13e7") == 0 || //MC DPMJET, anchr LHC13b+c
      periodName.CompareTo("LHC13b3") == 0 || //MC HIJING, weighted to number of events per run, anchor LHC13b
      periodName.CompareTo("LHC13b2") == 0 ||  // MC DPMJET, wrong energy, anchor LHC13b
      periodName.CompareTo("LHC13b2_plus") == 0 || // MC DPMJET, weighted to number event per run, anchor LHC13b
      periodName.CompareTo("LHC13c1_bis") == 0 || // MC AMPT fast generation, pT hardbin, anchor ?
      periodName.CompareTo("LHC13c1") == 0 || // MC AMPT fast generation, anchor ?
      periodName.CompareTo("LHC13b1") == 0 || // MC DPMJET, fragments, with fixed label 0, anchor LHC12g
      periodName.CompareTo("LHC12g4b_fix") == 0 || // MC DPMJET, with fixed label 0, anchor LHC12g
      periodName.CompareTo("LHC12g1_fix") == 0 || // MC ?, with fixed label 0, anchor LHC12g
      periodName.CompareTo("LHC12g4c") == 0 || // MC DPMJET, shifted vertex runs, anchor LHC12g
      periodName.CompareTo("LHC12h6") == 0 || // MC muon cocktail, anchor LHC12g
      periodName.CompareTo("LHC12g4b") == 0 || // MC DPMJET 3rd iteration, anchor LHC12g
      periodName.CompareTo("LHC12g4a") == 0 || // MC DPMJET improved, anchor LHC12g
      periodName.CompareTo("LHC12g4") == 0 || // MC DPMJET, anchor LHC12g
      periodName.CompareTo("LHC12g5") == 0 || // MC PHOJET, anchor LHC12g
      periodName.CompareTo("LHC12g2") == 0 || // MC Starlight background, anchor LHC12g
      periodName.CompareTo("LHC12g1") == 0 ) // MC ?, anchor LHC12g
      {
         printf(" Gamma Conversion Cuts %s :: pPb Run doing Eta Shift of %f \n\n",(GetCutNumber()).Data(),-0.465);
         SetEtaShift(-0.465);
      }
   else if(periodName.CompareTo("LHC13f") == 0 ||
           periodName.CompareTo("LHC13c6b") == 0 ||// MC Jpsi -> mumu, anchor LHC13f
           periodName.CompareTo("LHC13c5") == 0 || //MC Starlight, gamma gamma UPC muon, anchor LHC13f
           periodName.CompareTo("LHC13c4") == 0 )//MC Starlight, coherent JPsi, UPC muon, anchor LHC13f
      {
         printf(" Gamma Conversion Cuts %s :: Pbp Run doing Eta Shift of %f \n\n",(GetCutNumber()).Data(),0.465);
         SetEtaShift(+0.465);
      }
   else printf(" Gamma Conversion Cuts %s :: Automatic Eta Shift requested but Period is not known -> No Shift \n\n",(GetCutNumber()).Data());
}

///________________________________________________________________________
Bool_t AliConversionCuts::InPlaneOutOfPlaneCut(Double_t photonPhi, Double_t eventPlaneAngle, Bool_t fill){
   
   //GetPhotonPhi() 0-2 Pi  //eventPlaneAngle -1pi-1pi
   eventPlaneAngle=eventPlaneAngle+TMath::Pi();
   Double_t gammaToEPAngle = eventPlaneAngle-photonPhi;
   if(gammaToEPAngle < 0) gammaToEPAngle=gammaToEPAngle+2*TMath::Pi();
   gammaToEPAngle = gammaToEPAngle-TMath::Pi(); // angle from -pi +pi

   if(!fInPlaneOutOfPlane){
      if(fill&&hEventPlanePhi)hEventPlanePhi->Fill(gammaToEPAngle);
      return kTRUE;
   }
   else if(fInPlaneOutOfPlane == 1){
      if(TMath::Abs(gammaToEPAngle)<=0.25*TMath::Pi() || TMath::Abs(gammaToEPAngle)>=0.75*TMath::Pi()){
         if(fill&&hEventPlanePhi)hEventPlanePhi->Fill(gammaToEPAngle);
         return kTRUE;
      }
      else return kFALSE;
   }
   else if(fInPlaneOutOfPlane == 2){
      if(TMath::Abs(gammaToEPAngle)>0.25*TMath::Pi() && TMath::Abs(gammaToEPAngle)<0.75*TMath::Pi()){
         if(fill&&hEventPlanePhi)hEventPlanePhi->Fill(gammaToEPAngle);
         return kTRUE;
      }
      else return kFALSE;
   }
   return kFALSE;

}

///________________________________________________________________________
UChar_t AliConversionCuts::DeterminePhotonQualityAOD(AliAODConversionPhoton* photon, AliVEvent* eventDummy){

   AliAODTrack * negTrack = static_cast<AliAODTrack*>(GetTrack(eventDummy, photon->GetTrackLabelNegative()));
   AliAODTrack * posTrack = static_cast<AliAODTrack*>(GetTrack(eventDummy, photon->GetTrackLabelPositive()));

   if(!negTrack || !posTrack) {
        return 0;
   }
   if(negTrack->Charge() == posTrack->Charge()){
        return 0;
   }   
   Int_t nClusterITSneg = negTrack->GetITSNcls();
   Int_t nClusterITSpos = posTrack->GetITSNcls();
//    cout << nClusterITSneg << "\t" << nClusterITSpos <<endl;
   
   if (nClusterITSneg > 1 && nClusterITSpos > 1){
      return 3;
   } else if (nClusterITSneg > 1 || nClusterITSpos > 1){
      return 2;
   } else {
      return 1;
   }
   return 0;
   
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
// WARNING: this class is no longer supported,
// please use AliConversionPhotonCuts and AliConvEventCuts
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
