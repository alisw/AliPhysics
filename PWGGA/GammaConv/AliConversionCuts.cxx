/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *				       					  *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

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
#include "AliStack.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "AliLog.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"

class iostream;

using namespace std;

ClassImp(AliConversionCuts)


const char* AliConversionCuts::fgkCutNames[AliConversionCuts::kNCuts] = {
   "HeavyIon",
   "CentralityMin",
   "CentralityMax",
   "SelectV0AND",
   "MultiplicityMethod",
   "RemovePileUp",
   "RejectExtraSignals",
   "V0FinderType",
   "EtaCut",
   "MinRCut",
   "SinglePtCut",
   "ClsTPCCut",
   "ededxSigmaCut",
   "pidedxSigmaCut",
   "piMomdedxSigmaCut",
   "piMaxMomdedxSigmaCut",
   "LowPRejectionSigmaCut",
   "TOFelectronPID",
   "QtMaxCut",
   "Chi2GammaCut",
   "PsiPair",
   "DoPhotonAsymmetryCut",
   "CosinePointingAngle",
   "SharedElectronCuts",
   "RejectToCloseV0s"
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
   fDoHighPtQtGammaSelection(kFALSE),
   fQtMax(100),
   fHighPtQtMax(0.),
   fPtBorderForQt(0),
   fXVertexCut(0.),
   fYVertexCut(0.),
   fZVertexCut(0.),
   fNSigmaMass(0.),
   fUseEtaMinCut(kFALSE),
   fUseOnFlyV0Finder(kTRUE),
   fDoPhotonAsymmetryCut(kTRUE),
   fMinPPhotonAsymmetryCut(100.),
   fMinPhotonAsymmetry(0.),
   fIsHeavyIon(kFALSE),
   fDetectorCentrality(0),
   fModCentralityClass(0),
   fMaxVertexZ(10),
   fCentralityMin(0),
   fCentralityMax(0),
   fUseCorrectedTPCClsInfo(kFALSE),
   fUseTOFpid(kFALSE),
   fMultiplicityMethod(0),
   fSpecialTrigger(0),
   fRemovePileUp(kFALSE),
   fOpeningAngle(0.005),
   fPsiPairCut(10000),
   fCosPAngleCut(10000),
   fDoToCloseV0sCut(kFALSE),
   fRejectExtraSignals(0),
   fminV0Dist(200.),
   fDoSharedElecCut(kFALSE),
   fOfflineTriggerMask(0),
   fHasV0AND(kTRUE),
   fIsSDDFired(kTRUE),
   fRandom(0),
   fElectronArraySize(500),
   fElectronLabelArray(NULL),
   fConversionPointXArray(0.0),
   fConversionPointYArray(0.0),
   fConversionPointZArray(0.0),
   fnHeaders(0),
   fNotRejectedStart(NULL),
   fNotRejectedEnd(NULL),
   fGeneratorNames(NULL),
   fCutString(NULL),
   hdEdxCuts(NULL),
   hTPCdEdxbefore(NULL),
   hTPCdEdxafter(NULL),
   hTPCdEdxSigbefore(NULL),
   hTPCdEdxSigafter(NULL),
   hTOFbefore(NULL),
   hTOFSigbefore(NULL),
   hTOFSigafter(NULL),
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
   hVertexZ(NULL),
   hTriggerClass(NULL),
   hTriggerClassSelected(NULL)
{
   InitPIDResponse();
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
   fCutString=new TObjString((GetCutNumber()).Data());

   fElectronLabelArray = new Int_t[fElectronArraySize];

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
   fDoHighPtQtGammaSelection(ref.fDoHighPtQtGammaSelection),
   fQtMax(ref.fQtMax),
   fHighPtQtMax(ref.fHighPtQtMax),
   fPtBorderForQt(ref.fPtBorderForQt),
   fXVertexCut(ref.fXVertexCut),
   fYVertexCut(ref.fYVertexCut),
   fZVertexCut(ref.fZVertexCut),
   fNSigmaMass(ref.fNSigmaMass),
   fUseEtaMinCut(ref.fUseEtaMinCut),
   fUseOnFlyV0Finder(ref.fUseOnFlyV0Finder),
   fDoPhotonAsymmetryCut(fDoPhotonAsymmetryCut),
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
   fRemovePileUp(ref.fRemovePileUp),
   fOpeningAngle(ref.fOpeningAngle),
   fPsiPairCut(ref.fPsiPairCut),
   fCosPAngleCut(ref.fCosPAngleCut),
   fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
   fRejectExtraSignals(ref.fRejectExtraSignals),
   fminV0Dist(ref.fminV0Dist),
   fDoSharedElecCut(ref.fDoSharedElecCut),
   fOfflineTriggerMask(ref.fOfflineTriggerMask),
   fHasV0AND(ref.fHasV0AND),
   fIsSDDFired(ref.fIsSDDFired),
   fRandom(ref.fRandom),
   fElectronArraySize(ref.fElectronArraySize),
   fElectronLabelArray(NULL),
   fConversionPointXArray(ref.fConversionPointXArray),
   fConversionPointYArray(ref.fConversionPointYArray),
   fConversionPointZArray(ref.fConversionPointZArray),
   fnHeaders(ref.fnHeaders),
   fNotRejectedStart(NULL),
   fNotRejectedEnd(NULL),
   fGeneratorNames(ref.fGeneratorNames),
   fCutString(NULL),
   hdEdxCuts(NULL),
   hTPCdEdxbefore(NULL),
   hTPCdEdxafter(NULL),
   hTPCdEdxSigbefore(NULL),
   hTPCdEdxSigafter(NULL),
   hTOFbefore(NULL),
   hTOFSigbefore(NULL),
   hTOFSigafter(NULL),
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
   hVertexZ(NULL),
   hTriggerClass(NULL),
   hTriggerClassSelected(NULL)
{
   // Copy Constructor
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
   fCutString=new TObjString((GetCutNumber()).Data());
   fElectronLabelArray = new Int_t[fElectronArraySize];
   // dont copy histograms (if you like histograms, call InitCutHistograms())
}


//________________________________________________________________________
AliConversionCuts::~AliConversionCuts() {
   // Destructor
   //Deleting fHistograms leads to seg fault it it's added to output collection of a task
   // if(fHistograms)
   // 	delete fHistograms;
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
}

//________________________________________________________________________
void AliConversionCuts::InitCutHistograms(TString name, Bool_t preCut){

   // Initialize Cut Histograms for QA (only initialized and filled if function is called)

   if(fHistograms != NULL){
      delete fHistograms;
      fHistograms=NULL;
   }
   if(fHistograms==NULL){
      fHistograms=new TList();
      if(name=="")fHistograms->SetName(Form("ConvCuts_%s",GetCutNumber().Data()));
      else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
   }

   // IsPhotonSelected
   hCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",10,-0.5,9.5);
   hCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
   hCutIndex->GetXaxis()->SetBinLabel(kOnFly+1,"onfly");
   hCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
   hCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"dEdx");
   hCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
   hCutIndex->GetXaxis()->SetBinLabel(kConvPointFail+1,"ConvPoint fail");
   hCutIndex->GetXaxis()->SetBinLabel(kPhotonCuts+1,"PhotonCuts");
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
   hPhotonCuts=new TH1F(Form("PhotonCuts %s",GetCutNumber().Data()),"PhotonCuts",12,-0.5,11.5);
   hPhotonCuts->GetXaxis()->SetBinLabel(1,"in");
   hPhotonCuts->GetXaxis()->SetBinLabel(2,"qtcut");
   hPhotonCuts->GetXaxis()->SetBinLabel(3,"chi2");
   hPhotonCuts->GetXaxis()->SetBinLabel(4,"acceptance");
   hPhotonCuts->GetXaxis()->SetBinLabel(5,"asymmetry");
   hPhotonCuts->GetXaxis()->SetBinLabel(6,"pidprob");
   hPhotonCuts->GetXaxis()->SetBinLabel(7,"cortpcclinfo");
   hPhotonCuts->GetXaxis()->SetBinLabel(8,"PsiPair");
   hPhotonCuts->GetXaxis()->SetBinLabel(9,"CosPAngle");
   hPhotonCuts->GetXaxis()->SetBinLabel(10,"out");
   fHistograms->Add(hPhotonCuts);

   if(preCut){
      hInvMassbefore=new TH1F(Form("InvMass_before %s",GetCutNumber().Data()),"InvMass_before",1000,0,0.3);
      fHistograms->Add(hInvMassbefore);
      hArmenterosbefore=new TH2F(Form("Armenteros_before %s",GetCutNumber().Data()),"Armenteros_before",200,-1,1,1000,0,1.);
      fHistograms->Add(hArmenterosbefore);
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

      hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",100,0,100);
      fHistograms->Add(hCentrality);
      hVertexZ=new TH1F(Form("VertexZ %s",GetCutNumber().Data()),"VertexZ",1000,-50,50);
      fHistograms->Add(hVertexZ);

      hTriggerClass= new TH1F(Form("OfflineTrigger %s",GetCutNumber().Data()),"OfflineTrigger",34,-0.5,33.5);
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
      hTriggerClass->GetXaxis()->SetBinLabel(29,"kFastOnly");
      hTriggerClass->GetXaxis()->SetBinLabel(30,"kAnyINT");
      hTriggerClass->GetXaxis()->SetBinLabel(31,"kAny");
      hTriggerClass->GetXaxis()->SetBinLabel(32,"V0AND");
      hTriggerClass->GetXaxis()->SetBinLabel(33,"NOT kFastOnly");
      hTriggerClass->GetXaxis()->SetBinLabel(34,"failed Physics Selection");
      fHistograms->Add(hTriggerClass);

      hTriggerClassSelected= new TH1F(Form("OfflineTriggerSelected %s",GetCutNumber().Data()),"OfflineTriggerSelected",33,-0.5,32.5);
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
      hTriggerClassSelected->GetXaxis()->SetBinLabel(29,"kFastOnly");
      hTriggerClassSelected->GetXaxis()->SetBinLabel(30,"kAnyINT");
      hTriggerClassSelected->GetXaxis()->SetBinLabel(31,"kAny");
      hTriggerClassSelected->GetXaxis()->SetBinLabel(32,"V0AND");
      hTriggerClassSelected->GetXaxis()->SetBinLabel(33,"NOT kFastOnly");
      fHistograms->Add(hTriggerClassSelected);
   }
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
Bool_t AliConversionCuts::EventIsSelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent){
   // Process Event Selection

   Int_t cutindex=0;
   if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
   cutindex++;

   // Check for MC event
   if(fMCEvent){
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
   if(!IsTriggerSelected()){
      if(hV0EventCuts)hV0EventCuts->Fill(cutindex);
      fEventQuality = 3;
      return kFALSE;
   }
   cutindex++;

   if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliTriggerAnalysis fTriggerAnalysis;// = new AliTriggerAnalysis;
      fHasV0AND = fTriggerAnalysis.IsOfflineTriggerFired((AliESDEvent*)fInputEvent, AliTriggerAnalysis::kV0AND);
      if(fHasV0AND&&hTriggerClass)hTriggerClass->Fill(31);
   }

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

   fEventQuality = 0;
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionCuts::PhotonIsSelectedMC(TParticle *particle,AliStack *fMCStack,Bool_t checkForConvertedGamma){
   // MonteCarlo Photon Selection

   if(!fMCStack)return kFALSE;

   if (particle->GetPdgCode() == 22){

      if(particle->R() > fMaxR) return kFALSE;
      if(abs(particle->Eta())> fEtaCut || abs(particle->Eta())< fEtaCutMin) return kFALSE;

      if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
         return kFALSE; // no photon as mothers!
      }

      if(particle->GetMother(0) >= fMCStack->GetNprimary()){
         return kFALSE; // the gamma has a mother, and it is not a primary particle
      }

      if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

      // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
      TParticle* ePos = NULL;
      TParticle* eNeg = NULL;

      if(particle->GetNDaughters() >= 2){
         for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
            TParticle *tmpDaughter = fMCStack->Particle(daughterIndex);
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

      if( abs(ePos->Eta())> fEtaCut || abs(ePos->Eta())< fEtaCutMin ||
          abs(eNeg->Eta())> fEtaCut || abs(eNeg->Eta())< fEtaCutMin ) {
         return kFALSE;
      }

      if(ePos->R()>fMaxR){
         return kFALSE; // cuts on distance from collision point
      }

      if(abs(ePos->Vz()) > fMaxZ){
         return kFALSE;	 // outside material
      }
      if(abs(eNeg->Vz()) > fMaxZ){
         return kFALSE;	 // outside material
      }

      if( ePos->R() <= ((abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
         return kFALSE;  // line cut to exclude regions where we do not reconstruct
      } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((abs(ePos->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
         return kFALSE;
      }

      if( eNeg->R() <= ((abs(eNeg->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
         return kFALSE; // line cut to exclude regions where we do not reconstruct
      } else if ( fEtaCutMin != -0.1 &&   eNeg->R() >= ((abs(eNeg->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
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

   cutIndex++; //9
   if(hPhotonCuts)hPhotonCuts->Fill(cutIndex); //9

   // Histos after Cuts
   if(hInvMassafter)hInvMassafter->Fill(photon->GetMass());
   if(hArmenterosafter)hArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());


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

   // Get Tracks
   AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
   AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

   if(!negTrack || !posTrack) {
      FillPhotonCutIndex(kNoTracks);
      return kFALSE;
   }

   // dEdx Cuts
   if(!dEdxCuts(negTrack) || !dEdxCuts(posTrack)) {
      FillPhotonCutIndex(kdEdxCuts);
      return kFALSE;
   }

   // Track Cuts
   if(!TracksAreSelected(negTrack, posTrack)){
      FillPhotonCutIndex(kTrackCuts);
      return kFALSE;
   }

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

   if(fDoHighPtQtGammaSelection){
      if(photon->GetPhotonPt() < fPtBorderForQt){
         if(photon->GetArmenterosQt()>fQtMax){
            return kFALSE;
         }
      } else {
         if(photon->GetArmenterosQt()>fHighPtQtMax){
            return kFALSE;
         }
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

   if(photon->GetConversionRadius() <= ((abs(photon->GetConversionZ())*fLineCutZRSlope)-fLineCutZValue)){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   else if (fUseEtaMinCut &&  photon->GetConversionRadius() >= ((abs(photon->GetConversionZ())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   if(abs(photon->GetConversionZ()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;


   if(abs(photon->GetPhotonEta())> fEtaCut || abs(photon->GetPhotonEta())< fEtaCutMin){
      if(hAcceptanceCuts)hAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
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
   if((NegVtxType->GetType())==AliAODVertex::kKink  || (PosVtxType->GetType())==AliAODVertex::kKink) {
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

   if(abs(negTrack->Eta()) > fEtaCut || abs(negTrack->Eta()) < fEtaCutMin ||
      abs(posTrack->Eta())> fEtaCut || abs(posTrack->Eta())< fEtaCutMin) {
      if(hTrackCuts)hTrackCuts->Fill(cutIndex);
      return kFALSE;
   }
   cutIndex++;

   // Single Pt Cut
   if( negTrack->Pt()< fSinglePtCut ||	posTrack->Pt()< fSinglePtCut){
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
         if( abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
   }
   cutIndex++;

   if(fDoProtonRejectionLowP == kTRUE){
      if( fCurrentTrack->P()<fPIDMinPProtonRejectionLowP ){
         if( abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){

            if(hdEdxCuts)hdEdxCuts->Fill(cutIndex);
            return kFALSE;
         }
      }
   }
   cutIndex++;

   if(fDoPionRejectionLowP == kTRUE){
      if( fCurrentTrack->P()<fPIDMinPPionRejectionLowP ){
         if( abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){

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
         Double_t times[5];
         fCurrentTrack->GetIntegratedTimes(times);
         Double_t TOFsignal =	fCurrentTrack->GetTOFsignal();
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
      for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
         AliVTrack * track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));

         if(track) {
            if(track->GetID() == label) {
               return track;
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

      AliESDtrack* negTrack	= esdEvent->GetTrack(photon->GetTrackLabelNegative());
      AliESDtrack* posTrack	= esdEvent->GetTrack(photon->GetTrackLabelPositive());

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

   if( ePos->R() <= ((abs(ePos->Vz())*fLineCutZRSlope)-fLineCutZValue)){
      return kFALSE;
   }
   else if (fUseEtaMinCut &&  ePos->R() >= ((abs(ePos->Vz())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
      return kFALSE;
   }

   if(abs(eNeg->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      return kFALSE;
   }

   if(eNeg->Vz()!=ePos->Vz()||eNeg->R()!=ePos->R()){
      return kFALSE;
   }

   if(abs(ePos->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
      return kFALSE;
   }

   if(abs(particle->Eta())> fEtaCut || abs(particle->Eta())< fEtaCutMin){
      return kFALSE;
   }

   if(abs(ePos->Eta())> fEtaCut || abs(ePos->Eta())< fEtaCutMin){
      return kFALSE;
   }

   if(abs(eNeg->Eta())> fEtaCut || abs(eNeg->Eta())< fEtaCutMin){
      return kFALSE;
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
Bool_t AliConversionCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
   // Initialize Cuts from a given Cut string

   AliInfo(Form("Set Photoncut Number: %s",analysisCutSelection.Data()));
   if(analysisCutSelection.Length()!=kNCuts) {
      AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
      return kFALSE;
   }
   if(!analysisCutSelection.IsDigit()){
      AliError("Cut selection contains characters");
      return kFALSE;
   }

   const char *cutSelection = analysisCutSelection.Data();
#define ASSIGNARRAY(i)	fCuts[i] = cutSelection[i] - '0'
   for(Int_t ii=0;ii<kNCuts;ii++){
      ASSIGNARRAY(ii);
   }

   // Set Individual Cuts
   for(Int_t ii=0;ii<kNCuts;ii++){
      if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
   }

   //PrintCuts();

   // Set StandardTriggers

   if(fIsHeavyIon)SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   else SelectCollisionCandidates(AliVEvent::kMB);

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

   case kselectV0AND:
      if( SetSelectSpecialTrigger(value)) {
         fCuts[kselectV0AND] = value;
         UpdateCutString();
         return kTRUE;
      } else return kFALSE;

   case kmultiplicityMethod:
      if( SetMultiplicityMethod(value)) {
         fCuts[kmultiplicityMethod] = value;
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
   case 6: //allows to select centrality 0-45% in steps of 5% for V0 Multiplicity
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
Int_t AliConversionCuts::SetSelectSpecialTrigger(Int_t selectSpecialTrigger)
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

   default:
      AliError("Warning: Special Trigger Not known");
      return kFALSE;
   }
   return kTRUE;
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
      fUseOnFlyV0Finder=kTRUE;
      break;
   case 1:  // offline V0 finder
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
      fEtaCut		= 0.9;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 1:	// 1.2
      fEtaCut		= 1.2;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 2:	// 1.4
      fEtaCut		= 1.4;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 3: // 0.8
      fEtaCut		= 0.8;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 4: // 0.75
      fEtaCut		= 0.75;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 5: // 0.9 - 1.4
      fEtaCut		= 1.4;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= 0.9;
      fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
      break;
   case 6: // 5.
      fEtaCut		= 5.;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
   case 7: // 0.1 - 0.8
      fEtaCut		= 0.8;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= 0.1;
      fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
      break;
   case 8: // 0.1 - 0.8
      fEtaCut		= 0.9;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= 0.1;
      fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
      break;
   case 9: // 10
      fEtaCut		= 10;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin		= -0.1;
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
      // High purity cuts for PbPb (remove first layers of material)
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
      fMinR = 26.;
      break;
   case 8:
      fMaxR = 180.;
      fMinR = 35.;
      break;
   case 9:
      fMaxR = 35.;
      fMinR = 5.;
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
   case 1:  // 70
      fMinClsTPC= 70.;
      break;
   case 2:  // 80
      fMinClsTPC= 80.;
      break;
   case 3:  // 100
      fMinClsTPC= 100.;
      break;
   case 4:  // 60% of findable clusters
      fMinClsTPCToF= 0.6;
      fUseCorrectedTPCClsInfo=0;
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
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
      break;
   case 1:
      fQtMax=0.1;
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
      break;
   case 2:
      fQtMax=0.07;
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
      break;
   case 3:
      fQtMax=0.05;
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
      break;
   case 4:
      fQtMax=0.03;
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
      break;
   case 5:
      fQtMax=0.02;
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
      break;
   case 6:
      fQtMax=0.02;
      fDoHighPtQtGammaSelection=kTRUE;
      fHighPtQtMax=0.06;
      fPtBorderForQt=2.5;
      break;
   case 7:
      fQtMax=0.15;
      fDoHighPtQtGammaSelection=kFALSE;
      fHighPtQtMax=100.;
      fPtBorderForQt=100.;
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
      fPsiPairCut = 0.15; //
      break;
   case 5:
      fPsiPairCut = 0.2; //
      break;
   case 6:
      fPsiPairCut = 0.03; //
      break;
   case 7:
      fPsiPairCut = 0.025; //
      break;
   case 8:
      fPsiPairCut = 0.01; //
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
      fCosPAngleCut = TMath::Pi(); // -1
      break;
   case 1:
      fCosPAngleCut = 0.1; // 0.99500
      break;
   case 2:
      fCosPAngleCut = 0.05; // 0.99875
      break;
   case 3:
      fCosPAngleCut = 0.025; // 0.99969
      break;
   case 4:
      fCosPAngleCut = 0.01; // 0.99995
      break;
   case 5:
      fCosPAngleCut = 0.2; // 0.98007
      break;
   case 6:
      fCosPAngleCut = 0.5; // 0.87758
      break;
   case 7:
      fCosPAngleCut = 0.075; // 0.73169
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
      break;
   case 1:
      fDoSharedElecCut = kTRUE;
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
//-------------------------------------------------------------
Double_t AliConversionCuts::GetCentrality(AliVEvent *event)
{   // Get Event Centrality

   AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
   if(esdEvent){
      AliCentrality *fESDCentrality=(AliCentrality*)esdEvent->GetCentrality();

      if(fDetectorCentrality==0){
         return fESDCentrality->GetCentralityPercentile("V0M"); // default
      }
      if(fDetectorCentrality==1){
         return fESDCentrality->GetCentralityPercentile("CL1");
      }
   }

   AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(event);
   if(aodEvent){
      if(aodEvent->GetHeader()){return aodEvent->GetHeader()->GetCentrality();}
   }

   return -1;
}
//-------------------------------------------------------------
Bool_t AliConversionCuts::IsCentralitySelected(AliVEvent *event, AliVEvent *fMCEvent)
{   // Centrality Selection
   if(!fIsHeavyIon)return kTRUE;

   if(fCentralityMin == 0 && fCentralityMax == 0) return kTRUE;//0-100%
   if(fCentralityMin >= fCentralityMax) return kTRUE;//0-100%

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

   // Use strict V0 amplitude cut for MC centrality
   Float_t nv0amplitude = event->GetVZEROData()->GetMTotV0A()+event->GetVZEROData()->GetMTotV0C();
   Float_t V0Amplitude10[10] = {9999999.0,13670,9345,6209,3944,2352,1272,611,255, 83};
   //                                    0    10   20   30   40   50   60  70  80  90%
   Float_t V0Amplitude5a[10] = {9999999.0,16612,13670,11290,9345,7650,6209,4984,3944,3074};
   //                                    0     5    10    15   20   25   30   35   40   45%
   Float_t V0Amplitude5b[10] = {3074,2352,1725,1272,899,611,402,255,152,83};
   //                             45   50   55   60  65  70  75  80  85 90%

   if (fModCentralityClass == 3){
      if(fMCEvent){
         if(nv0amplitude > V0Amplitude10[fCentralityMax] && nv0amplitude <= V0Amplitude10[fCentralityMin])
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
         if(nv0amplitude > V0Amplitude5a[fCentralityMax] && nv0amplitude <= V0Amplitude5a[fCentralityMin])
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
         if(nv0amplitude > V0Amplitude5b[fCentralityMax] && nv0amplitude <= V0Amplitude5b[fCentralityMin])
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

   if(abs(fVertexZ)>fMaxVertexZ)return kFALSE;
   return kTRUE;
}
///________________________________________________________________________

Int_t AliConversionCuts::GetNumberOfContributorsVtx(AliVEvent *event){
   // returns number of contributors to the vertex

   AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(event);
   if(fESDEvent){
      if (fESDEvent->GetPrimaryVertex() != NULL){
         if(fESDEvent->GetPrimaryVertex()->GetNContributors()>0) {
            return fESDEvent->GetPrimaryVertex()->GetNContributors();
         }
      }

      if(fESDEvent->GetPrimaryVertexSPD() !=NULL){
         if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
            return fESDEvent->GetPrimaryVertexSPD()->GetNContributors();
         }  else {
            AliWarning(Form("Number of contributors from bad vertex type:: %s",fESDEvent->GetPrimaryVertex()->GetName()));
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

   return 0;
}

///________________________________________________________________________

Bool_t AliConversionCuts::IsTriggerSelected()
{

   AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

   UInt_t isSelected = AliVEvent::kAny;
   if( fInputHandler && fInputHandler->GetEventSelection()) {
      // Get the actual offline trigger mask for the event and AND it with the
      // requested mask. If no mask requested select by default the event.
      if (fOfflineTriggerMask)
         isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
   }
   fIsSDDFired = !(fInputHandler->IsEventSelected() & AliVEvent::kFastOnly);

   // Fill Histogram
   if(hTriggerClass){
      if (fIsSDDFired) hTriggerClass->Fill(32);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClass->Fill(0);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)hTriggerClass->Fill(1);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUON)hTriggerClass->Fill(2);
      if (fInputHandler->IsEventSelected() & AliVEvent::kHighMult)hTriggerClass->Fill(3);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)hTriggerClass->Fill(4);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCINT5)hTriggerClass->Fill(5);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCMUS5)hTriggerClass->Fill(6);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSPB)hTriggerClass->Fill(6);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSH7)hTriggerClass->Fill(7);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSHPB)hTriggerClass->Fill(7);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUL7)hTriggerClass->Fill(8);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikePB)hTriggerClass->Fill(8);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUU7)hTriggerClass->Fill(9);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikePB)hTriggerClass->Fill(9);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7)hTriggerClass->Fill(10);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC8)hTriggerClass->Fill(10);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUS7)hTriggerClass->Fill(11);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI1)hTriggerClass->Fill(12);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI7)hTriggerClass->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI8)hTriggerClass->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHOSPb)hTriggerClass->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)hTriggerClass->Fill(14);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA)hTriggerClass->Fill(15);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)hTriggerClass->Fill(16);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)hTriggerClass->Fill(17);
      if (fInputHandler->IsEventSelected() & AliVEvent::kDG5)hTriggerClass->Fill(18);
      if (fInputHandler->IsEventSelected() & AliVEvent::kZED)hTriggerClass->Fill(19);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSPI7)hTriggerClass->Fill(20);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSPI)hTriggerClass->Fill(20);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT8)hTriggerClass->Fill(21);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt8)hTriggerClass->Fill(22);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt8)hTriggerClass->Fill(23);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt8)hTriggerClass->Fill(24);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8)hTriggerClass->Fill(25);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt0)hTriggerClass->Fill(26);
      if (fInputHandler->IsEventSelected() & AliVEvent::kUserDefined)hTriggerClass->Fill(27);
      if (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly)hTriggerClass->Fill(28);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAnyINT)hTriggerClass->Fill(29);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAny)hTriggerClass->Fill(30);
      if (!fInputHandler->IsEventSelected()) hTriggerClass->Fill(33);
   }

   if(hTriggerClassSelected && isSelected){
      if (!fIsSDDFired) hTriggerClassSelected->Fill(32);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMB)hTriggerClassSelected->Fill(0);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)hTriggerClassSelected->Fill(1);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUON)hTriggerClassSelected->Fill(2);
      if (fInputHandler->IsEventSelected() & AliVEvent::kHighMult)hTriggerClassSelected->Fill(3);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1)hTriggerClassSelected->Fill(4);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCINT5)hTriggerClassSelected->Fill(5);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCMUS5)hTriggerClassSelected->Fill(6);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSPB)hTriggerClassSelected->Fill(6);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSH7)hTriggerClassSelected->Fill(7);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUSHPB)hTriggerClassSelected->Fill(7);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUL7)hTriggerClassSelected->Fill(8);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikePB)hTriggerClassSelected->Fill(8);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUU7)hTriggerClassSelected->Fill(9);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikePB)hTriggerClassSelected->Fill(9);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7)hTriggerClassSelected->Fill(10);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMC8)hTriggerClassSelected->Fill(10);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMUS7)hTriggerClassSelected->Fill(11);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI1)hTriggerClassSelected->Fill(12);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI7)hTriggerClassSelected->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHI8)hTriggerClassSelected->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kPHOSPb)hTriggerClassSelected->Fill(13);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)hTriggerClassSelected->Fill(14);
      if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA)hTriggerClassSelected->Fill(15);
      if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)hTriggerClassSelected->Fill(16);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)hTriggerClassSelected->Fill(17);
      if (fInputHandler->IsEventSelected() & AliVEvent::kDG5)hTriggerClassSelected->Fill(18);
      if (fInputHandler->IsEventSelected() & AliVEvent::kZED)hTriggerClassSelected->Fill(19);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSPI7)hTriggerClassSelected->Fill(20);
      if (fInputHandler->IsEventSelected() & AliVEvent::kSPI)hTriggerClassSelected->Fill(20);
      if (fInputHandler->IsEventSelected() & AliVEvent::kINT8)hTriggerClassSelected->Fill(21);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt8)hTriggerClassSelected->Fill(22);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt8)hTriggerClassSelected->Fill(23);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt8)hTriggerClassSelected->Fill(24);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8)hTriggerClassSelected->Fill(25);
      if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt0)hTriggerClassSelected->Fill(26);
      if (fInputHandler->IsEventSelected() & AliVEvent::kUserDefined)hTriggerClassSelected->Fill(27);
      if (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly)hTriggerClassSelected->Fill(28);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAnyINT)hTriggerClassSelected->Fill(29);
      if (fInputHandler->IsEventSelected() & AliVEvent::kAny)hTriggerClassSelected->Fill(30);
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
   if(GetCosineOfPointingAngle(photon, event) < (TMath::Cos(fCosPAngleCut))){
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

   if(photon->GetPsiPair() > fPsiPairCut){
      return kFALSE;}
   else{return kTRUE;}

}

///________________________________________________________________________
TString AliConversionCuts::GetCutNumber(){
   // returns TString with current cut number
   TString a(kNCuts);
   for(Int_t ii=0;ii<kNCuts;ii++){
      a.Append(Form("%d",fCuts[ii]));
   }
   return a;
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
void AliConversionCuts::GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliMCEvent *MCEvent){

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
   AliGenCocktailEventHeader *cHeader = dynamic_cast<AliGenCocktailEventHeader*>(MCEvent->GenEventHeader());
   if(cHeader){
      TList *genHeaders = cHeader->GetHeaders();
      AliGenEventHeader* gh = 0;
      fnHeaders = 0;
      if(rejection == 1 || rejection == 3) fnHeaders = 1; // MinBiasHeader
      if(rejection == 2){ // TList of Headers Names
         for(Int_t i = 0; i<genHeaders->GetEntries();i++){
            gh = (AliGenEventHeader*)genHeaders->At(i);
            TString GeneratorName = gh->GetName();
            for(Int_t j = 0; j<HeaderList->GetEntries();j++){
               TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
               if(GeneratorName.CompareTo(GeneratorInList) == 0){
                  fnHeaders++;
                  continue;
               }
            }
         }
      }

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
      Int_t nummer = 0;
      for(Int_t i = 0; i<genHeaders->GetEntries();i++){
         gh = (AliGenEventHeader*)genHeaders->At(i);
         TString GeneratorName = gh->GetName();
         lastindex = lastindex + gh->NProduced();
         for(Int_t j = 0; j<HeaderList->GetEntries();j++){
            TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
            if(GeneratorName.CompareTo(GeneratorInList) == 0){
               fNotRejectedStart[nummer] = firstindex;
               fNotRejectedEnd[nummer] = lastindex;
               fGeneratorNames[nummer] = GeneratorName;
//                cout << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << endl;
               nummer++;
               continue;
            }
         }
         firstindex = firstindex + gh->NProduced();
      }
   } else { // No Cocktail Header Found
      fNotRejectedStart = new Int_t[1];
      fNotRejectedEnd = new Int_t[1];

      fnHeaders = 1;
      fNotRejectedStart[0] = 0;
      fNotRejectedEnd[0] = MCEvent->Stack()->GetNprimary()-1;
//       if(rejection == 2){
      fGeneratorNames = new TString[1];
      fGeneratorNames[0] = "NoCocktailGeneratorFound";
//       }
            
      AliGenPythiaEventHeader *mcHeaderPythia = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent->GenEventHeader());      
      if (mcHeaderPythia) fGeneratorNames[0] = "NoCocktailGeneratorFound_Pythia";
      AliGenDPMjetEventHeader *mcHeaderPhojet = dynamic_cast<AliGenDPMjetEventHeader*>(MCEvent->GenEventHeader());      
      if (mcHeaderPhojet) fGeneratorNames[0] = "NoCocktailGeneratorFound_Phojet";
      AliGenHijingEventHeader *mcHeaderHijing = dynamic_cast<AliGenHijingEventHeader*>(MCEvent->GenEventHeader());      
      if (mcHeaderHijing) fGeneratorNames[0] = "NoCocktailGeneratorFound_Hijing";

      SetRejectExtraSignalsCut(0);
   }

}
//_________________________________________________________________________
Int_t AliConversionCuts::IsParticleFromBGEvent(Int_t index, AliStack *MCStack){

   // Not Accepted == kFALSE == 0
   //     Accepted ==  kTRUE == 1
   //  FirstHeader ==  kTRUE == 3

   if(index < 0) return 0; // No Particle

   Int_t accepted = 0;
   if( index >= MCStack->GetNprimary()){ // Secondary Particle
      if( ((TParticle*)MCStack->Particle(index))->GetMother(0) < 0) return 1; // Secondary Particle without Mother??
      return IsParticleFromBGEvent(((TParticle*)MCStack->Particle(index))->GetMother(0),MCStack);
   }
   for(Int_t i = 0;i<fnHeaders;i++){
      if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
         accepted = 1;
         if(i == 0) accepted = 2; // MB Header
      }
   }

   return accepted;
}
//_________________________________________________________________________
Int_t AliConversionCuts::IsEventAcceptedByConversionCut(AliConversionCuts *ReaderCuts, AliVEvent *InputEvent, AliMCEvent *MCEvent, Bool_t isHeavyIon){

   if(isHeavyIon && !(IsCentralitySelected(InputEvent,MCEvent)))
      return 1; // Check Centrality --> Not Accepted => eventQuality = 1

   if(!isHeavyIon && GetIsFromPileup()){
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


   return 0;
}

Float_t AliConversionCuts::GetWeightForMeson(TString period, Int_t index, AliStack *MCStack){
   if (!(period.CompareTo("LHC12f1a") == 0 || period.CompareTo("LHC12f1b") == 0  || period.CompareTo("LHC12i3") == 0 )) return 1.;
   
   Int_t kCaseGen = 0;
   for (Int_t i = 0; i < fnHeaders; i++){
      if (index >= fNotRejectedStart[i] && index < fNotRejectedEnd[i]+1){
//          cout << fGeneratorNames[i].Data() << endl;
         if (fGeneratorNames[i].CompareTo("Pythia") == 0){
            kCaseGen = 1;
         } else if (fGeneratorNames[i].CompareTo("DPMJET") == 0){
            kCaseGen = 2;
         } else if (fGeneratorNames[i].CompareTo("HIJING") == 0){
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
//          cout << "resulting kCaseGen :" << kCaseGen << endl;
      }
   }
   
   Double_t mesonPt = ((TParticle*)MCStack->Particle(index))->Pt();
   
//   Double_t mesonY = 10.;
//   if(((TParticle*)MCStack->Particle(index))->Energy() - ((TParticle*)MCStack->Particle(index))->Pz() == 0 || ((TParticle*)MCStack->Particle(index))->Energy() + ((TParticle*)MCStack->Particle(index))->Pz() == 0){
//      mesonY=10.;
//   } else{
//      mesonY = 0.5*(TMath::Log((((TParticle*)MCStack->Particle(index))->Energy()+((TParticle*)MCStack->Particle(index))->Pz()) / (((TParticle*)MCStack->Particle(index))->Energy()-((TParticle*)MCStack->Particle(index))->Pz())));
//   }
   Double_t mesonMass = ((TParticle*)MCStack->Particle(index))->GetCalcMass();
   Float_t functionResult = 1.;
   if (kCaseGen == 1){
      Float_t dNdyMC = 2.1462;
      Float_t nMC = 7.06055;
      Float_t tMC = 0.12533;
      if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  111){
         dNdyMC = 2.1462;
         nMC = 7.06055;
         tMC = 0.12533;
      } else if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  221){
         dNdyMC = 0.2357;
         nMC = 5.9105;
         tMC = 0.1525;
      }
      functionResult = dNdyMC / ( 2 * TMath::Pi())*(nMC-1.)*(nMC-2.) / (nMC*tMC*(nMC*tMC+mesonMass*(nMC-2.)))  * pow(1.+(sqrt(mesonPt*mesonPt+mesonMass*mesonMass)-mesonMass)/(nMC*tMC), -nMC);
   } else if (kCaseGen == 2){
      Float_t dNdyMC = 2.35978;
      Float_t nMC = 6.81795;
      Float_t tMC = 0.11492;
      if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  111){
         dNdyMC = 2.35978;
         nMC = 6.81795;
         tMC = 0.11492;
      } else if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  221){
         dNdyMC = 0.3690;
         nMC = 5.55809;
         tMC = 0.13387;
      }
      functionResult = dNdyMC / ( 2 * TMath::Pi())*(nMC-1.)*(nMC-2.) / (nMC*tMC*(nMC*tMC+mesonMass*(nMC-2.)))  * pow(1.+(sqrt(mesonPt*mesonPt+mesonMass*mesonMass)-mesonMass)/(nMC*tMC), -nMC);
   } else if (kCaseGen == 4){
//       functionResult = 1./sqrt(1.-mesonMass*mesonMass/((mesonMass*mesonMass+mesonPt*mesonPt)*cosh(mesonY)*cosh(mesonY)));
      Float_t a = 0.23437;
      Float_t b = 5.6661;
      Float_t c = -1430.5863;
      Float_t d = -0.6966624;
      Float_t e = 252.3742;
      if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  111){
         a = 0.23437;
         b = 5.6661;
         c = -1430.5863;
         d = -0.6966624;
         e = 252.3742;
      } else if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  221){
         a = 0.10399;
         b = 4.35311;
         c = -12.17723;
         d = -0.01172;
         e =1.85140;
      }
      functionResult = a*TMath::Power(mesonPt,-1.*(b+c/(TMath::Power(mesonPt,d)+e)))*1./mesonPt *1./1.6 *1./(2.* TMath::Pi());
   }
   
   Float_t dNdyData = 2.2328;
   Float_t nData = 7.1473;
   Float_t tData = 0.1346;
   if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  111){
         dNdyData = 2.2328;
         nData = 7.1473;
         tData = 0.1346;
      } else if ( ((TParticle*)MCStack->Particle(index))->GetPdgCode() ==  221){
         dNdyData = 0.38992; //be careful this fit is not optimal, eta in data still has problems
         nData = 5.72778;
         tData = 0.13835;
      }
   Float_t tsallisData = dNdyData / ( 2 * TMath::Pi())*(nData-1.)*(nData-2.) / (nData*tData*(nData*tData+mesonMass*(nData-2.)))  * pow(1.+(sqrt(mesonPt*mesonPt+mesonMass*mesonMass)-mesonMass)/(nData*tData), -nData);
//    cout << "tsallisData/functionResult: " << tsallisData/functionResult << endl;
   return tsallisData/functionResult;
//       return 
}

///________________________________________________________________________
AliConversionCuts* AliConversionCuts::GetStandardCuts2010PbPb(){
    //Create and return standard 2010 PbPb cuts
    AliConversionCuts *cuts=new AliConversionCuts("StandardCuts2010PbPb","StandardCuts2010PbPb");
    if(!cuts->InitializeCutsFromCutString("1000002042092970023220000")){
	cout<<"Warning: Initialization of Standardcuts2010PbPb failed"<<endl;}
    return cuts;
}

///________________________________________________________________________
AliConversionCuts* AliConversionCuts::GetStandardCuts2010pp(){
    //Create and return standard 2010 PbPb cuts
    AliConversionCuts *cuts=new AliConversionCuts("StandardCuts2010pp","StandardCuts2010pp");
    if(!cuts->InitializeCutsFromCutString("0000011002093663003800000")){
	cout<<"Warning: Initialization of Standardcuts2010pp failed"<<endl;}
    return cuts;
}

