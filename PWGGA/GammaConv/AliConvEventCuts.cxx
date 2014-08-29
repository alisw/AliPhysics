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
#include "AliStack.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
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

class iostream;

using namespace std;

ClassImp(AliConvEventCuts)


const char* AliConvEventCuts::fgkCutNames[AliConvEventCuts::kNCuts] = {
   "HeavyIon",						//0
   "CentralityMin",					//1
   "CentralityMax",					//2
   "SelectSpecialTrigger",			//3
   "SelectSpecialSubTriggerClass",	//4
   "RemovePileUp",					//5
   "RejectExtraSignals",			//6
};


//________________________________________________________________________
AliConvEventCuts::AliConvEventCuts(const char *name,const char *title) :
   AliAnalysisCuts(name,title),
   fHistograms(NULL),
   fHeaderList(NULL),
   fEventQuality(-1),
   fIsHeavyIon(0),
   fDetectorCentrality(0),
   fModCentralityClass(0),
   fMaxVertexZ(10),
   fCentralityMin(0),
   fCentralityMax(0),
   fMultiplicityMethod(0),
   fSpecialTrigger(0),
   fSpecialSubTrigger(0),
   fRemovePileUp(kFALSE),
   fRejectExtraSignals(0),
   fOfflineTriggerMask(0),
   fHasV0AND(kTRUE),
   fIsSDDFired(kTRUE),
   fRandom(0),
   fnHeaders(0),
   fNotRejectedStart(NULL),
   fNotRejectedEnd(NULL),
   fGeneratorNames(NULL),
   fCutString(NULL),
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
   fHistoEventCuts(NULL),
   hCentrality(NULL),
   hCentralityVsNumberOfPrimaryTracks(NULL),
   hVertexZ(NULL),
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
   fEventQuality(ref.fEventQuality),
   fIsHeavyIon(ref.fIsHeavyIon),
   fDetectorCentrality(ref.fDetectorCentrality),
   fModCentralityClass(ref.fModCentralityClass),
   fMaxVertexZ(ref.fMaxVertexZ),
   fCentralityMin(ref.fCentralityMin),
   fCentralityMax(ref.fCentralityMax),
   fMultiplicityMethod(ref.fMultiplicityMethod),
   fSpecialTrigger(ref.fSpecialTrigger),
   fSpecialSubTrigger(ref.fSpecialSubTrigger),
   fRemovePileUp(ref.fRemovePileUp),
   fRejectExtraSignals(ref.fRejectExtraSignals),
   fOfflineTriggerMask(ref.fOfflineTriggerMask),
   fHasV0AND(ref.fHasV0AND),
   fIsSDDFired(ref.fIsSDDFired),
   fRandom(ref.fRandom),
   fnHeaders(ref.fnHeaders),
   fNotRejectedStart(NULL),
   fNotRejectedEnd(NULL),
   fGeneratorNames(ref.fGeneratorNames),
   fCutString(NULL),
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
   fHistoEventCuts(NULL),
   hCentrality(NULL),
   hCentralityVsNumberOfPrimaryTracks(NULL),
   hVertexZ(NULL),
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

   hCentrality=new TH1F(Form("Centrality %s",GetCutNumber().Data()),"Centrality",100,0,100);
   fHistograms->Add(hCentrality);
   hCentralityVsNumberOfPrimaryTracks=new TH2F(Form("Centrality vs Primary Tracks %s",GetCutNumber().Data()),"Centrality vs Primary Tracks ",100,0,100,4000,0,4000);
   fHistograms->Add(hCentralityVsNumberOfPrimaryTracks);

   // Event Cuts and Info
   if(preCut){
      fHistoEventCuts=new TH1F(Form("ESD_EventCuts %s",GetCutNumber().Data()),"Event Cuts",7,-0.5,6.5);
      fHistoEventCuts->GetXaxis()->SetBinLabel(1,"in");
      fHistoEventCuts->GetXaxis()->SetBinLabel(2,"OfflineTrigger");
      fHistoEventCuts->GetXaxis()->SetBinLabel(3,"nvtxcontr");
      fHistoEventCuts->GetXaxis()->SetBinLabel(4,"VertexZ");
      fHistoEventCuts->GetXaxis()->SetBinLabel(5,"pileup");
      fHistoEventCuts->GetXaxis()->SetBinLabel(6,"centrsel");
      fHistoEventCuts->GetXaxis()->SetBinLabel(7,"out");
      fHistograms->Add(fHistoEventCuts);

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
      
   }
   TH1::AddDirectory(kTRUE);
}

///________________________________________________________________________
Bool_t AliConvEventCuts::EventIsSelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent){
   // Process Event Selection

   Int_t cutindex=0;
   if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
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
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
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
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      fEventQuality = 5;
      return kFALSE;
   }
   cutindex++;

   // Z Vertex Position Cut
   if(!VertexZCut(fInputEvent)){
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      fEventQuality = 4;
      return kFALSE;
   }
   cutindex++;

   // Pile Up Rejection

   if(fRemovePileUp){
      if(fInputEvent->IsPileupFromSPD(3,0.8,3.,2.,5.)){
         if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
         fEventQuality = 6;
         return kFALSE;
      }
   }
   cutindex++;

   // Centrality Selection
   if(!IsCentralitySelected(fInputEvent,fMCEvent)){
      if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
      fEventQuality = 1;
      return kFALSE;
   }
   cutindex++;

   // Fill Event Histograms
   if(fHistoEventCuts)fHistoEventCuts->Fill(cutindex);
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
	// Initialize Cuts from a given Cut string
	if(fDoReweightHistoMCPi0 || fDoReweightHistoMCEta || fDoReweightHistoMCK0s) {
		AliInfo("Weighting was enabled");
		LoadReweightingHistosMCFromFile();
	}

	AliInfo(Form("Set Event Cut Number: %s",analysisCutSelection.Data()));
	if(analysisCutSelection.Length()!=kNCuts) {
		AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
		return kFALSE;
	}
	if(!analysisCutSelection.IsDigit()){
		AliError("Cut selection contains characters");
		return kFALSE;
	}

	const char *cutSelection = analysisCutSelection.Data();
	#define ASSIGNARRAY(i)  fCuts[i] = cutSelection[i] - '0'
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

	switch (cutID) {
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
Bool_t AliConvEventCuts::SetSelectSubTriggerClass(Int_t selectSpecialSubTriggerClass)
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

///________________________________________________________________________
Bool_t AliConvEventCuts::SetRemovePileUp(Int_t removePileUp)
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

//-------------------------------------------------------------
Double_t AliConvEventCuts::GetCentrality(AliVEvent *event)
{   // Get Event Centrality

	AliESDEvent *esdEvent=dynamic_cast<AliESDEvent*>(event);
	if(esdEvent){
		AliCentrality *fESDCentrality=(AliCentrality*)esdEvent->GetCentrality();
		if(fDetectorCentrality==0){
			if (fIsHeavyIon==2){
				return fESDCentrality->GetCentralityPercentile("V0A"); // default for pPb
			} else{
				return fESDCentrality->GetCentralityPercentile("V0M"); // default
			}
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

//_____________________________________________________________________________________
Bool_t AliConvEventCuts::IsCentralitySelected(AliVEvent *event, AliVEvent *fMCEvent)
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
	Int_t PrimaryTracks10[10][2] =
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
			{   0,   0}  // 90
		};
	Int_t PrimaryTracks5a[10][2] =
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
			{ 260, 214}  // 45
		};
	Int_t PrimaryTracks5b[10][2] =
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
			{   0,   0}  // 90
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
	
	if(abs(fVertexZ)>fMaxVertexZ)return kFALSE;

	TString periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()
													->GetTask("V0ReaderV1"))->GetPeriodName();
	if (periodName.CompareTo("LHC11h")==0){
		if (abs(fVertexZ-fVertexZSPD) > 0.1) return kFALSE;
	}						
	if (fIsHeavyIon == 2){
		if(fUtils->IsFirstEventInChunk(event)) return kFALSE;
		if(!fUtils->IsVertexSelected2013pA(event)) return kFALSE;
		if(fUtils->IsPileUpEvent(event)) return kFALSE;
	}

	return kTRUE;
}

///________________________________________________________________________
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
Bool_t AliConvEventCuts::IsTriggerSelected(AliVEvent *fInputEvent)
{

	AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	
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
				if (fSpecialSubTrigger>0){
					if (!firedTrigClass.Contains(fSpecialSubTriggerName.Data())) isSelected = 0;
				}	 
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
TString AliConvEventCuts::GetCutNumber(){
   // returns TString with current cut number
   TString a(kNCuts);
   for(Int_t ii=0;ii<kNCuts;ii++){
      a.Append(Form("%d",fCuts[ii]));
   }
   return a;
}

///________________________________________________________________________
void AliConvEventCuts::GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *MCEvent){

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
	AliStack *fMCStack = 0x0;
	TClonesArray *fMCStackAOD = 0x0;
	if(MCEvent->IsA()==AliMCEvent::Class()){
		cHeader = dynamic_cast<AliGenCocktailEventHeader*>(dynamic_cast<AliMCEvent*>(MCEvent)->GenEventHeader());
		if(cHeader) headerFound = kTRUE;
		if(dynamic_cast<AliMCEvent*>(MCEvent))fMCStack = dynamic_cast<AliStack*>(dynamic_cast<AliMCEvent*>(MCEvent)->Stack());
	}
	if(MCEvent->IsA()==AliAODEvent::Class()){ // MCEvent is a AODEvent in case of AOD
		cHeaderAOD = dynamic_cast<AliAODMCHeader*>(MCEvent->FindListObject(AliAODMCHeader::StdBranchName()));
		fMCStackAOD = dynamic_cast<TClonesArray*>(MCEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		
		
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
// 					cout << GeneratorInList.Data() << endl;
					if(GeneratorName.CompareTo(GeneratorInList) == 0){
// 						cout << "accepted" << endl;
						if (GeneratorInList.CompareTo("PARAM") == 0 || GeneratorInList.CompareTo("BOX") == 0 ){
							if(fMCStack){
								if (fMCStack->Particle(firstindexA)->GetPdgCode() == fAddedSignalPDGCode ) {
									if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
										if (gh->NProduced() > 10 && fMCStack->Particle(firstindexA+10)->GetPdgCode() == fAddedSignalPDGCode ){
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
							if ( fMCStackAOD){
								AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fMCStackAOD->At(firstindexA));
								if (  aodMCParticle->GetPdgCode() == fAddedSignalPDGCode ){
									if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
										if (gh->NProduced() > 10){
											AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(fMCStackAOD->At(firstindexA+10));
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
						if(fMCStack){
							if (fMCStack->Particle(firstindex)->GetPdgCode() == fAddedSignalPDGCode ) {
								if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
// 									cout << "produced " << gh->NProduced() << " with box generator" << endl;
									if (gh->NProduced() > 10 && fMCStack->Particle(firstindex+10)->GetPdgCode() == fAddedSignalPDGCode){
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
						if ( fMCStackAOD){
							AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fMCStackAOD->At(firstindex));
							if (  aodMCParticle->GetPdgCode() == fAddedSignalPDGCode ){
								if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0 ){
									if (gh->NProduced() > 10) {
										AliAODMCParticle *aodMCParticle2 = static_cast<AliAODMCParticle*>(fMCStackAOD->At(firstindex+10));
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
		fNotRejectedEnd[0] = static_cast<AliMCEvent*>(MCEvent)->Stack()->GetNprimary()-1;
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
Int_t AliConvEventCuts::IsParticleFromBGEvent(Int_t index, AliStack *MCStack, AliVEvent *InputEvent){

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
		if( index >= MCStack->GetNprimary()){ // Secondary Particle
			if( ((TParticle*)MCStack->Particle(index))->GetMother(0) < 0) return 1; // Secondary Particle without Mother??
			return IsParticleFromBGEvent(((TParticle*)MCStack->Particle(index))->GetMother(0),MCStack,InputEvent);
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
				return IsParticleFromBGEvent(aodMCParticle->GetMother(),MCStack,InputEvent);
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
Int_t AliConvEventCuts::IsEventAcceptedByCut(AliConvEventCuts *ReaderCuts, AliVEvent *InputEvent, AliMCEvent *MCEvent, Int_t isHeavyIon){

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
Float_t AliConvEventCuts::GetWeightForMeson(TString period, Int_t index, AliStack *MCStack, AliVEvent *InputEvent){
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
		mesonPt = ((TParticle*)MCStack->Particle(index))->Pt();
		mesonMass = ((TParticle*)MCStack->Particle(index))->GetCalcMass();
		PDGCode = ((TParticle*)MCStack->Particle(index))->GetPdgCode();
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
void AliConvEventCuts::GetCorrectEtaShiftFromPeriod(TString periodName){

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

