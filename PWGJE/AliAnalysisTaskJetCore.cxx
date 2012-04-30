
// ******************************************
// This task computes several jet observables like 
// the fraction of energy in inner and outer coronnas,
// jet-track correlations,triggered jet shapes and 
// correlation strength distribution of particles inside jets.    
// Author: lcunquei@cern.ch
// *******************************************


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


#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliAnalysisTaskFastEmbedding.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"

#include "AliAnalysisTaskJetCore.h"

ClassImp(AliAnalysisTaskJetCore)

AliAnalysisTaskJetCore::AliAnalysisTaskJetCore() :
AliAnalysisTaskSE(),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fBackgroundBranch(""),
fNonStdFile(""),
fIsPbPb(kTRUE),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.),
fVtxZMax(10.),
fEvtClassMin(0),
fEvtClassMax(4),
fFilterMask(0),
fRadioFrac(0.2),
fMinDist(0.1),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fAngStructCloseTracks(0),
fCheckMethods(0),
fDoEventMixing(0), 
fJetEtaMin(-.5),
fJetEtaMax(.5),
fNevents(0x0),
fTindex(0x0),
fTrigBufferIndex(0x0),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fHistEvtSelection(0x0),
fhnDeltaR(0x0),
fhnMixedEvents(0x0),
fh2JetCoreMethod1C10(0x0),
fh2JetCoreMethod2C10(0x0),
fh2JetCoreMethod1C20(0x0),
fh2JetCoreMethod2C20(0x0),
fh2JetCoreMethod1C30(0x0),
fh2JetCoreMethod2C30(0x0),
fh2JetCoreMethod1C60(0x0),
fh2JetCoreMethod2C60(0x0),
fh2AngStructpt1C10(0x0),
fh2AngStructpt2C10(0x0),
fh2AngStructpt3C10(0x0),
fh2AngStructpt4C10(0x0),
fh2AngStructpt1C20(0x0),
fh2AngStructpt2C20(0x0),
fh2AngStructpt3C20(0x0),
fh2AngStructpt4C20(0x0),    
fh2AngStructpt1C30(0x0),
fh2AngStructpt2C30(0x0),
fh2AngStructpt3C30(0x0),
fh2AngStructpt4C30(0x0),   
fh2AngStructpt1C60(0x0),
fh2AngStructpt2C60(0x0),
fh2AngStructpt3C60(0x0),
fh2AngStructpt4C60(0x0),
fh2JetsumHT3R2a(0x0),
fh2JetsumHT3R2ap(0x0),
fh2JetsumHT3R4a(0x0),
fh2JetsumHT3R4ap(0x0),
fh2JetsumHT3R6a(0x0),
fh2JetsumHT3R6ap(0x0),
fh2JetsumHT3R8a(0x0),
fh2JetsumHT3R8ap(0x0),
fh2JetsumHT3R10a(0x0),
fh2JetsumHT3R10ap(0x0),
fh2JetsumHT3R2aa(0x0),
fh2JetsumHT3R2aap(0x0),
fh2JetsumHT3R4aa(0x0),
fh2JetsumHT3R4aap(0x0),
fh2JetsumHT3R6aa(0x0),
fh2JetsumHT3R6aap(0x0),
fh2JetsumHT3R8aa(0x0),
fh2JetsumHT3R8aap(0x0),
fh2JetsumHT3R10aa(0x0),
fh2JetsumHT3R10aap(0x0),
fh2JetsumHT3R2aaa(0x0),
fh2JetsumHT3R2aaap(0x0),
fh2JetsumHT3R4aaa(0x0),
fh2JetsumHT3R4aaap(0x0),
fh2JetsumHT3R6aaa(0x0),
fh2JetsumHT3R6aaap(0x0),
fh2JetsumHT3R8aaa(0x0),
fh2JetsumHT3R8aaap(0x0),
fh2JetsumHT3R10aaa(0x0),
fh2JetsumHT3R10aaap(0x0),
fh2JetsumHT3R2b(0x0),
fh2JetsumHT3R2bp(0x0),
fh2JetsumHT3R4b(0x0),
fh2JetsumHT3R4bp(0x0),
fh2JetsumHT3R6b(0x0),
fh2JetsumHT3R6bp(0x0),
fh2JetsumHT3R8b(0x0),
fh2JetsumHT3R8bp(0x0),
fh2JetsumHT3R10b(0x0),
fh2JetsumHT3R10bp(0x0),
fh2JetsumHT3R2bb(0x0),
fh2JetsumHT3R2bbp(0x0),
fh2JetsumHT3R4bb(0x0),
fh2JetsumHT3R4bbp(0x0),
fh2JetsumHT3R6bb(0x0),
fh2JetsumHT3R6bbp(0x0),
fh2JetsumHT3R8bb(0x0),
fh2JetsumHT3R8bbp(0x0),
fh2JetsumHT3R10bb(0x0),
fh2JetsumHT3R10bbp(0x0),
fh2JetsumHT3R2bbb(0x0),
fh2JetsumHT3R2bbbp(0x0),
fh2JetsumHT3R4bbb(0x0),
fh2JetsumHT3R4bbbp(0x0),
fh2JetsumHT3R6bbb(0x0),
fh2JetsumHT3R6bbbp(0x0),
fh2JetsumHT3R8bbb(0x0),
fh2JetsumHT3R8bbbp(0x0),
fh2JetsumHT3R10bbb(0x0),
fh2JetsumHT3R10bbbp(0x0),
fh3spectriggeredC10(0x0),
fh3spectriggeredC20(0x0),
fh3spectriggeredC3060(0x0),
fh3specbiased(0x0),
fh3spectot(0x0),
fh3spectotb(0x0)
 
{
   // default Constructor


 // Trigger buffer.
   for(Int_t i=0; i<10; i++) {
		for(Int_t j=0; j<7; j++) {
			fTrigBuffer[i][j]=0;
		}				
	}	





   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;
}

AliAnalysisTaskJetCore::AliAnalysisTaskJetCore(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fBackgroundBranch(""),
fNonStdFile(""),
fIsPbPb(kTRUE),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.),
fVtxZMax(10.),
fEvtClassMin(0),
fEvtClassMax(4),
fFilterMask(0),
fRadioFrac(0.2),
fMinDist(0.1),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fAngStructCloseTracks(0),
fCheckMethods(0),
fDoEventMixing(0),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fNevents(0x0),
fTindex(0x0),
fTrigBufferIndex(0x0),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fHistEvtSelection(0x0),
fhnDeltaR(0x0),
fhnMixedEvents(0x0),
fh2JetCoreMethod1C10(0x0),
fh2JetCoreMethod2C10(0x0),
fh2JetCoreMethod1C20(0x0),
fh2JetCoreMethod2C20(0x0),
fh2JetCoreMethod1C30(0x0),
fh2JetCoreMethod2C30(0x0),
fh2JetCoreMethod1C60(0x0),
fh2JetCoreMethod2C60(0x0),
fh2AngStructpt1C10(0x0),
fh2AngStructpt2C10(0x0),
fh2AngStructpt3C10(0x0),
fh2AngStructpt4C10(0x0),
fh2AngStructpt1C20(0x0),
fh2AngStructpt2C20(0x0),
fh2AngStructpt3C20(0x0),
fh2AngStructpt4C20(0x0),    
fh2AngStructpt1C30(0x0),
fh2AngStructpt2C30(0x0),
fh2AngStructpt3C30(0x0),
fh2AngStructpt4C30(0x0),   
fh2AngStructpt1C60(0x0),
fh2AngStructpt2C60(0x0),
fh2AngStructpt3C60(0x0),
fh2AngStructpt4C60(0x0),    
fh2JetsumHT3R2a(0x0),
fh2JetsumHT3R2ap(0x0),
fh2JetsumHT3R4a(0x0),
fh2JetsumHT3R4ap(0x0),
fh2JetsumHT3R6a(0x0),
fh2JetsumHT3R6ap(0x0),
fh2JetsumHT3R8a(0x0),
fh2JetsumHT3R8ap(0x0),
fh2JetsumHT3R10a(0x0),
fh2JetsumHT3R10ap(0x0),
fh2JetsumHT3R2aa(0x0),
fh2JetsumHT3R2aap(0x0),
fh2JetsumHT3R4aa(0x0),
fh2JetsumHT3R4aap(0x0),
fh2JetsumHT3R6aa(0x0),
fh2JetsumHT3R6aap(0x0),
fh2JetsumHT3R8aa(0x0),
fh2JetsumHT3R8aap(0x0),
fh2JetsumHT3R10aa(0x0),
fh2JetsumHT3R10aap(0x0),
fh2JetsumHT3R2aaa(0x0),
fh2JetsumHT3R2aaap(0x0),
fh2JetsumHT3R4aaa(0x0),
fh2JetsumHT3R4aaap(0x0),
fh2JetsumHT3R6aaa(0x0),
fh2JetsumHT3R6aaap(0x0),
fh2JetsumHT3R8aaa(0x0),
fh2JetsumHT3R8aaap(0x0),
fh2JetsumHT3R10aaa(0x0),
fh2JetsumHT3R10aaap(0x0),
fh2JetsumHT3R2b(0x0),
fh2JetsumHT3R2bp(0x0),
fh2JetsumHT3R4b(0x0),
fh2JetsumHT3R4bp(0x0),
fh2JetsumHT3R6b(0x0),
fh2JetsumHT3R6bp(0x0),
fh2JetsumHT3R8b(0x0),
fh2JetsumHT3R8bp(0x0),
fh2JetsumHT3R10b(0x0),
fh2JetsumHT3R10bp(0x0),
fh2JetsumHT3R2bb(0x0),
fh2JetsumHT3R2bbp(0x0),
fh2JetsumHT3R4bb(0x0),
fh2JetsumHT3R4bbp(0x0),
fh2JetsumHT3R6bb(0x0),
fh2JetsumHT3R6bbp(0x0),
fh2JetsumHT3R8bb(0x0),
fh2JetsumHT3R8bbp(0x0),
fh2JetsumHT3R10bb(0x0),
fh2JetsumHT3R10bbp(0x0),
fh2JetsumHT3R2bbb(0x0),
fh2JetsumHT3R2bbbp(0x0),
fh2JetsumHT3R4bbb(0x0),
fh2JetsumHT3R4bbbp(0x0),
fh2JetsumHT3R6bbb(0x0),
fh2JetsumHT3R6bbbp(0x0),
fh2JetsumHT3R8bbb(0x0),
fh2JetsumHT3R8bbbp(0x0),
fh2JetsumHT3R10bbb(0x0),
fh2JetsumHT3R10bbbp(0x0),
fh3spectriggeredC10(0x0),
fh3spectriggeredC20(0x0),
fh3spectriggeredC3060(0x0),
fh3specbiased(0x0),
fh3spectot(0x0),
fh3spectotb(0x0)
 {
   // Constructor


    for(Int_t i=0; i<10; i++) {
       for(Int_t j=0; j<7; j++) {
	    fTrigBuffer[i][j]=0;
		}				
    }	



   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;

   DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetCore::~AliAnalysisTaskJetCore()
{
   delete fListJets[0];
   delete fListJets[1];
}

void AliAnalysisTaskJetCore::SetBranchNames(const TString &branch1, const TString &branch2)
{
   fJetBranchName[0] = branch1;
   fJetBranchName[1] = branch2;
}

void AliAnalysisTaskJetCore::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskJetCore::UserCreateOutputObjects()
{
   // Create histograms
   // Called once
   OpenFile(1);
   if(!fOutputList) fOutputList = new TList;
   fOutputList->SetOwner(kTRUE);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);


   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");

     UInt_t entries = 0; // bit coded, see GetDimParams() below 
     entries = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 |1<<7; 
     fhnDeltaR = NewTHnSparseF("fhnDeltaR", entries);

     //change binning in pTtrack
     Double_t *xPt3 = new Double_t[10];
     xPt3[0] = 0.;
     for(int i = 1; i<=9;i++){
      if(xPt3[i-1]<1)xPt3[i] = xPt3[i-1] + 0.2; // 1 - 5
      else if(xPt3[i-1]<10)xPt3[i] = xPt3[i-1] + 3; // 5 - 12
      else xPt3[i] = xPt3[i-1] + 150.; // 18 
     }
    fhnDeltaR->SetBinEdges(2,xPt3);
    delete [] xPt3;

    //change binning in HTI
     Double_t *xPt4 = new Double_t[14];
     xPt4[0] = 0.;
     for(int i = 1; i<=13;i++){
      if(xPt4[i-1]<10)xPt4[i] = xPt4[i-1] + 1; // 1 - 10
      else if(xPt4[i-1]<20)xPt4[i] = xPt4[i-1] + 5; // 10 - 12
      else xPt4[i] = xPt4[i-1] + 30.; // 13
     }
    fhnDeltaR->SetBinEdges(6,xPt4);
    delete [] xPt4;

    




   
     if(fDoEventMixing){    
     UInt_t cifras = 0; // bit coded, see GetDimParams() below 
     cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 |1<<7; 
     fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", cifras);}

    if(fCheckMethods){

    fh2JetCoreMethod1C10 = new TH2F("JetCoreMethod1C10","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C10 = new TH2F("JetCoreMethod2C10","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod1C20 = new TH2F("JetCoreMethod1C20","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C20 = new TH2F("JetCoreMethod2C20","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod1C30 = new TH2F("JetCoreMethod1C30","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C30 = new TH2F("JetCoreMethod2C30","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod1C60 = new TH2F("JetCoreMethod1C60","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C60 = new TH2F("JetCoreMethod2C60","",150, 0., 150.,100, 0., 1.5);}

   
    if(fAngStructCloseTracks>0){
    fh2AngStructpt1C10 = new TH2F("Ang struct pt1 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C10 = new TH2F("Ang struct pt2 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C10 = new TH2F("Ang struct pt3 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C10 = new TH2F("Ang struct pt4 C10","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt1C20 = new TH2F("Ang struct pt1 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C20 = new TH2F("Ang struct pt2 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C20 = new TH2F("Ang struct pt3 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C20 = new TH2F("Ang struct pt4 C20","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt1C30 = new TH2F("Ang struct pt1 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C30 = new TH2F("Ang struct pt2 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C30 = new TH2F("Ang struct pt3 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C30 = new TH2F("Ang struct pt4 C30","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt1C60 = new TH2F("Ang struct pt1 C60","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt2C60 = new TH2F("Ang struct pt2 C60","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt3C60 = new TH2F("Ang struct pt3 C60","",15,0.,1.5,150,0.,10.);
    fh2AngStructpt4C60 = new TH2F("Ang struct pt4 C60","",15,0.,1.5,150,0.,10.); }

    fh2JetsumHT3R2a = new TH2F("Pt sum R02 HT0 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R2ap = new TH2F("Pt sum R02 HT0 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4a = new TH2F("Pt sum R04 HT0 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4ap = new TH2F("Pt sum R04 HT0 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6a = new TH2F("Pt sum R06 HT0 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6ap = new TH2F("Pt sum R06 HT0 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8a = new TH2F("Pt sum R08 HT0 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8ap = new TH2F("Pt sum R08 HT0 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10a = new TH2F("Pt sum R10 HT0 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10ap = new TH2F("Pt sum R10 HT0 TT10 p","",20,0.,200.,100,0.,10.);  

    fh2JetsumHT3R2aa = new TH2F("Pt sum R02 HT0 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R2aap = new TH2F("Pt sum R02 HT0 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4aa = new TH2F("Pt sum R04 HT0 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4aap = new TH2F("Pt sum R04 HT0 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6aa = new TH2F("Pt sum R06 HT0 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6aap = new TH2F("Pt sum R06 HT0 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8aa = new TH2F("Pt sum R08 HT0 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8aap = new TH2F("Pt sum R08 HT0 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10aa = new TH2F("Pt sum R10 HT0 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10aap = new TH2F("Pt sum R10 HT0 TT20 p","",20,0.,200.,100,0.,10.);  

    fh2JetsumHT3R2aaa = new TH2F("Pt sum R02 HT0 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R2aaap = new TH2F("Pt sum R02 HT0 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4aaa = new TH2F("Pt sum R04 HT0 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4aaap = new TH2F("Pt sum R04 HT0 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6aaa = new TH2F("Pt sum R06 HT0 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6aaap = new TH2F("Pt sum R06 HT0 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8aaa = new TH2F("Pt sum R08 HT0 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8aaap = new TH2F("Pt sum R08 HT0 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10aaa = new TH2F("Pt sum R10 HT0 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10aaap = new TH2F("Pt sum R10 HT0 TT0 p","",20,0.,200.,100,0.,10.);  

    fh2JetsumHT3R2b = new TH2F("Pt sum R02 HT6 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R2bp = new TH2F("Pt sum R02 HT6 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4b = new TH2F("Pt sum R04 HT6 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4bp = new TH2F("Pt sum R04 HT6 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6b = new TH2F("Pt sum R06 HT6 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6bp = new TH2F("Pt sum R06 HT6 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8b = new TH2F("Pt sum R08 HT6 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8bp = new TH2F("Pt sum R08 HT6 TT10 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10b = new TH2F("Pt sum R10 HT6 TT10","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10bp = new TH2F("Pt sum R10 HT6 TT10 p","",20,0.,200.,100,0.,10.);  

    fh2JetsumHT3R2bb = new TH2F("Pt sum R02 HT6 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R2bbp = new TH2F("Pt sum R02 HT6 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4bb = new TH2F("Pt sum R04 HT6 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4bbp = new TH2F("Pt sum R04 HT6 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6bb = new TH2F("Pt sum R06 HT6 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6bbp = new TH2F("Pt sum R06 HT6 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8bb = new TH2F("Pt sum R08 HT6 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8bbp = new TH2F("Pt sum R08 HT6 TT20 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10bb = new TH2F("Pt sum R10 HT6 TT20","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10bbp = new TH2F("Pt sum R10 HT6 TT20 p","",20,0.,200.,100,0.,10.);  

    fh2JetsumHT3R2bbb = new TH2F("Pt sum R02 HT6 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R2bbbp = new TH2F("Pt sum R02 HT6 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4bbb = new TH2F("Pt sum R04 HT6 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R4bbbp = new TH2F("Pt sum R04 HT6 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6bbb = new TH2F("Pt sum R06 HT6 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R6bbbp = new TH2F("Pt sum R06 HT6 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8bbb = new TH2F("Pt sum R08 HT6 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R8bbbp = new TH2F("Pt sum R08 HT6 TT0 p","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10bbb = new TH2F("Pt sum R10 HT6 TT0","",20,0.,200.,100,0.,10.);
    fh2JetsumHT3R10bbbp = new TH2F("Pt sum R10 HT6 TT0 p","",20,0.,200.,100,0.,10.);  

   



     fh3spectriggeredC10 = new TH3F("Triggered spectrumC10","",10,0.,1.,100,-200.,200.,50,0.,50.);
     fh3spectriggeredC20 = new TH3F("Triggered spectrumC20","",10,0.,1.,100,-200.,200.,50,0.,50.);
     fh3spectriggeredC3060 = new TH3F("Triggered spectrumC3060","",10,0.,1.,100,-200.,200.,50,0.,50.);

    fh3specbiased = new TH3F("Biased spectrum","",10,0,100,50,0.,200.,50,0.,50.);
    fh3spectot = new TH3F("Total spectrum 0-10","",50,0.,200.,50,0.,50.,50,0.,50.);
    fh3spectotb = new TH3F("Total spectrum 30-60","",50,0.,200.,50,0.,50.,50,0.,50.);    
   fOutputList->Add(fHistEvtSelection);

   fOutputList->Add(fhnDeltaR);
   
   fOutputList->Add(fhnMixedEvents);
         
     
  
      if(fCheckMethods){

      fOutputList->Add(fh2JetCoreMethod1C10);
      fOutputList->Add(fh2JetCoreMethod2C10);
      fOutputList->Add(fh2JetCoreMethod1C20);
      fOutputList->Add(fh2JetCoreMethod2C20);
      fOutputList->Add(fh2JetCoreMethod1C30);
      fOutputList->Add(fh2JetCoreMethod2C30);
      fOutputList->Add(fh2JetCoreMethod1C60);
      fOutputList->Add(fh2JetCoreMethod2C60);}
      
      
     


        if(fAngStructCloseTracks>0){
       fOutputList->Add(fh2AngStructpt1C10);
       fOutputList->Add(fh2AngStructpt2C10);
       fOutputList->Add(fh2AngStructpt3C10);
       fOutputList->Add(fh2AngStructpt4C10); 
       fOutputList->Add(fh2AngStructpt1C20);
       fOutputList->Add(fh2AngStructpt2C20);
       fOutputList->Add(fh2AngStructpt3C20);
       fOutputList->Add(fh2AngStructpt4C20); 
       fOutputList->Add(fh2AngStructpt1C30);
       fOutputList->Add(fh2AngStructpt2C30);
       fOutputList->Add(fh2AngStructpt3C30);
       fOutputList->Add(fh2AngStructpt4C30);
       fOutputList->Add(fh2AngStructpt1C60);
       fOutputList->Add(fh2AngStructpt2C60);
       fOutputList->Add(fh2AngStructpt3C60);
       fOutputList->Add(fh2AngStructpt4C60);}  


	fOutputList->Add(fh2JetsumHT3R2a);
        fOutputList->Add(fh2JetsumHT3R2ap);
       	fOutputList->Add(fh2JetsumHT3R4a);
        fOutputList->Add(fh2JetsumHT3R4ap);
        fOutputList->Add(fh2JetsumHT3R6a);
        fOutputList->Add(fh2JetsumHT3R6ap);
       	fOutputList->Add(fh2JetsumHT3R8a);
        fOutputList->Add(fh2JetsumHT3R8ap);
        fOutputList->Add(fh2JetsumHT3R10a);
        fOutputList->Add(fh2JetsumHT3R10ap);
	fOutputList->Add(fh2JetsumHT3R2aa);
        fOutputList->Add(fh2JetsumHT3R2aap);
       	fOutputList->Add(fh2JetsumHT3R4aa);
        fOutputList->Add(fh2JetsumHT3R4aap);
        fOutputList->Add(fh2JetsumHT3R6aa);
        fOutputList->Add(fh2JetsumHT3R6aap);
       	fOutputList->Add(fh2JetsumHT3R8aa);
        fOutputList->Add(fh2JetsumHT3R8aap);
        fOutputList->Add(fh2JetsumHT3R10aa);
        fOutputList->Add(fh2JetsumHT3R10aap);
        fOutputList->Add(fh2JetsumHT3R2aaa);
        fOutputList->Add(fh2JetsumHT3R2aaap);
       	fOutputList->Add(fh2JetsumHT3R4aaa);
        fOutputList->Add(fh2JetsumHT3R4aaap);
        fOutputList->Add(fh2JetsumHT3R6aaa);
        fOutputList->Add(fh2JetsumHT3R6aaap);
       	fOutputList->Add(fh2JetsumHT3R8aaa);
        fOutputList->Add(fh2JetsumHT3R8aaap);
        fOutputList->Add(fh2JetsumHT3R10aaa);
        fOutputList->Add(fh2JetsumHT3R10aaap);

        fOutputList->Add(fh2JetsumHT3R2b);
        fOutputList->Add(fh2JetsumHT3R2bp);
       	fOutputList->Add(fh2JetsumHT3R4b);
        fOutputList->Add(fh2JetsumHT3R4bp);
        fOutputList->Add(fh2JetsumHT3R6b);
        fOutputList->Add(fh2JetsumHT3R6bp);
       	fOutputList->Add(fh2JetsumHT3R8b);
        fOutputList->Add(fh2JetsumHT3R8bp);
        fOutputList->Add(fh2JetsumHT3R10b);
        fOutputList->Add(fh2JetsumHT3R10bp);
	fOutputList->Add(fh2JetsumHT3R2bb);
        fOutputList->Add(fh2JetsumHT3R2bbp);
       	fOutputList->Add(fh2JetsumHT3R4bb);
        fOutputList->Add(fh2JetsumHT3R4bbp);
        fOutputList->Add(fh2JetsumHT3R6bb);
        fOutputList->Add(fh2JetsumHT3R6bbp);
       	fOutputList->Add(fh2JetsumHT3R8bb);
        fOutputList->Add(fh2JetsumHT3R8bbp);
        fOutputList->Add(fh2JetsumHT3R10bb);
        fOutputList->Add(fh2JetsumHT3R10bbp);
        fOutputList->Add(fh2JetsumHT3R2bbb);
        fOutputList->Add(fh2JetsumHT3R2bbbp);
       	fOutputList->Add(fh2JetsumHT3R4bbb);
        fOutputList->Add(fh2JetsumHT3R4bbbp);
        fOutputList->Add(fh2JetsumHT3R6bbb);
        fOutputList->Add(fh2JetsumHT3R6bbbp);
       	fOutputList->Add(fh2JetsumHT3R8bbb);
        fOutputList->Add(fh2JetsumHT3R8bbbp);
        fOutputList->Add(fh2JetsumHT3R10bbb);
        fOutputList->Add(fh2JetsumHT3R10bbbp);

 

       fOutputList->Add(fh3spectriggeredC10);
       fOutputList->Add(fh3spectriggeredC20); 
       fOutputList->Add(fh3spectriggeredC3060);   

       fOutputList->Add(fh3specbiased);
       fOutputList->Add(fh3spectot);
       fOutputList->Add(fh3spectotb); 
   // =========== Switch on Sumw2 for all histos ===========
   for (Int_t i=0; i<fOutputList->GetEntries(); ++i) {
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if (h1){
         h1->Sumw2();
         continue;
      }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
      if (hn){
         hn->Sumw2();
      }	  
   }
   TH1::AddDirectory(oldStatus);

   PostData(1, fOutputList);
}

void AliAnalysisTaskJetCore::UserExec(Option_t *)
{
   

   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
      return;
   }

   fESD=dynamic_cast<AliESDEvent*>(InputEvent());
   if (!fESD) {
      AliError("ESD not available");
      fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
   } 
      fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());

       static AliAODEvent* aod = 0;
       // take all other information from the aod we take the tracks from
       if(!aod){
       if(!fESD)aod = fAODIn;
       else aod = fAODOut;}

   
 
    if(fNonStdFile.Length()!=0){
    // case that we have an AOD extension we need can fetch the jets from the extended output
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
    if(!fAODExtension){
      if(fDebug>1)Printf("AODExtension found for %s",fNonStdFile.Data());
    }}
    




   // -- event selection --
   fHistEvtSelection->Fill(1); // number of events before event selection

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
   ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   cout<<inputHandler->IsEventSelected()<<" "<<fOfflineTrgMask<<endl;
   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      fHistEvtSelection->Fill(2);
      PostData(1, fOutputList);
      return;
   }

   // vertex selection
   if(!aod){
     if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
     fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
   }
   AliAODVertex* primVtx = aod->GetPrimaryVertex();

   if(!primVtx){
     if(fDebug) Printf("%s:%d No primVtx",(char*)__FILE__,__LINE__);
     fHistEvtSelection->Fill(3);
     PostData(1, fOutputList);
     return;
   }

   Int_t nTracksPrim = primVtx->GetNContributors();
   if ((nTracksPrim < fMinContribVtx) ||
         (primVtx->GetZ() < fVtxZMin) ||
         (primVtx->GetZ() > fVtxZMax) ){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }

   // event class selection (from jet helper task)
   Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   if(fDebug) Printf("Event class %d", eventClass);
   if (eventClass < fEvtClassMin || eventClass > fEvtClassMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fOutputList);
      return;
   }

   // centrality selection
   AliCentrality *cent = 0x0;
   Double_t centValue = 0.; 
   if(fESD) {cent = fESD->GetCentrality();
     if(cent) centValue = cent->GetCentralityPercentile("V0M");}
   else     centValue=aod->GetHeader()->GetCentrality();
   
   if(fDebug) printf("centrality: %f\n", centValue);
      if (centValue < fCentMin || centValue > fCentMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fOutputList);
      return;
    }


   fHistEvtSelection->Fill(0); 
   // accepted events  
   // -- end event selection --
  
   // get background
   AliAODJetEventBackground* externalBackground = 0;
   if(fAODOut&&!externalBackground&&fBackgroundBranch.Length()){
      externalBackground =  (AliAODJetEventBackground*)(fAODOut->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   if(fAODExtension&&!externalBackground&&fBackgroundBranch.Length()){
     externalBackground =  (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }

    if(fAODIn&&!externalBackground&&fBackgroundBranch.Length()){
      externalBackground =  (AliAODJetEventBackground*)(fAODIn->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
    }
   
   Float_t rho = 0;
   if(externalBackground)rho = externalBackground->GetBackground(0);


   // fetch jets
   TClonesArray *aodJets[2];
   aodJets[0]=0;
   if(fAODOut&&!aodJets[0]){
   aodJets[0] = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName[0].Data())); 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName[1].Data()));  }
   if(fAODExtension && !aodJets[0]){ 
   aodJets[0] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[0].Data())); 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[1].Data()));  }
     if(fAODIn&&!aodJets[0]){
   aodJets[0] = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName[0].Data())); 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName[1].Data()));  } 


   //Double_t ptsub[aodJets[0]->GetEntriesFast()];
   //Int_t inord[aodJets[0]->GetEntriesFast()];
   //for(Int_t n=0;n<aodJets[0]->GetEntriesFast();n++){
   //  ptsub[n]=0;
   //  inord[n]=0;}   

   TList ParticleList;
   Int_t nT = GetListOfTracks(&ParticleList);
     for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
      fListJets[iJetType]->Clear();
      if (!aodJets[iJetType]) continue;

      if(fDebug) Printf("%s: %d jets",fJetBranchName[iJetType].Data(),aodJets[iJetType]->GetEntriesFast());
      
   
      for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
         if (jet) fListJets[iJetType]->Add(jet);
	 // if(iJetType==0){
	 // ptsub[iJet]=jet->Pt()-rho*jet->EffectiveAreaCharged();}
      }}
   
   Double_t etabig=0;
   Double_t ptbig=0;
   Double_t areabig=0;
   Double_t phibig=0.;
   Double_t etasmall=0;
   Double_t ptsmall=0;
   Double_t areasmall=0;
   //Double_t distr=0.;
   Double_t phismall=0.;
         
  
   // Double_t up1[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t up2[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t up3[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t up4[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t down1[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t down2[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t down3[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // Double_t down4[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Int_t iCount=0; 
   Int_t trigJet=-1;
   Int_t trigBBTrack=-1;
   Int_t trigInTrack=-1;
     
   for(Int_t i=0; i<fListJets[0]->GetEntries(); ++i){
           AliAODJet* jetbig = (AliAODJet*)(fListJets[0]->At(i));
           etabig  = jetbig->Eta();
           phibig  = jetbig->Phi();
           ptbig   = jetbig->Pt();
           if(ptbig==0) continue; 
           areabig = jetbig->EffectiveAreaCharged();
           Double_t ptcorr=ptbig-rho*areabig;
       
      	   if((etabig<fJetEtaMin)||(etabig>fJetEtaMax)) continue;
                   Double_t dismin=100.;
                   Double_t ptmax=-10.; 
                   Int_t index1=-1;
                   Int_t index2=-1;
                  
           Int_t point=GetHardestTrackBackToJet(jetbig);    
	   AliVParticle *partback = (AliVParticle*)ParticleList.At(point);                            
           if(!partback) continue; 
	   if(centValue<10.)  fh3spectriggeredC10->Fill(jetbig->EffectiveAreaCharged(),ptcorr,partback->Pt());
           if(centValue<20.)  fh3spectriggeredC20->Fill(jetbig->EffectiveAreaCharged(),ptcorr,partback->Pt());
           if(centValue>30. && centValue<60.)  fh3spectriggeredC3060->Fill(jetbig->EffectiveAreaCharged(),ptcorr,partback->Pt());

                   if(ptcorr<=0) continue;
		   //if(partback->Pt()<6.) continue;
                       AliAODTrack* leadtrack; 
                       Int_t ippt=0;
                       Double_t ppt=-10;   
		       TRefArray *genTrackList = jetbig->GetRefTracks();
                       Int_t nTracksGenJet = genTrackList->GetEntriesFast();
                       AliAODTrack* genTrack;
                       for(Int_t ir=0; ir<nTracksGenJet; ++ir){
                       genTrack = (AliAODTrack*)(genTrackList->At(ir));
		       if(genTrack->Pt()>ppt){ppt=genTrack->Pt();
		       ippt=ir;}}
                        leadtrack=(AliAODTrack*)(genTrackList->At(ippt));
                        if(!leadtrack) continue;
                        fh3specbiased->Fill(centValue,ptcorr,leadtrack->Pt());
                        if(centValue<10)fh3spectot->Fill(ptcorr,leadtrack->Pt(),partback->Pt());
                        if(centValue>30. && centValue<60.)fh3spectotb->Fill(ptcorr,leadtrack->Pt(),partback->Pt()); 
			//store one trigger info                   
                        if((partback->Pt()>10.)&&(iCount==0)){                        
			  trigJet=i;
                          trigBBTrack=point;
                          trigInTrack=ippt;
                          iCount=iCount+1;} 


		  if(fCheckMethods){
                  for(Int_t j=0; j<fListJets[1]->GetEntries(); ++j){
                  AliAODJet* jetsmall = (AliAODJet*)(fListJets[1]->At(j));
                  etasmall  = jetsmall->Eta();
                  phismall = jetsmall->Phi();
                  ptsmall   = jetsmall->Pt();
                  areasmall = jetsmall->EffectiveAreaCharged();
                  Double_t tmpDeltaR=(phismall-phibig)*(phismall-phibig)+(etasmall-etabig)*(etasmall-etabig);
		  tmpDeltaR=TMath::Sqrt(tmpDeltaR);
		     //Fraction in the jet core  
                    if((ptsmall>ptmax)&&(tmpDeltaR<=fRadioFrac)){ptmax=ptsmall;  
		    index2=j;}  
                    if(tmpDeltaR<=dismin){ dismin=tmpDeltaR;
		      index1=j;}} //en of loop over R=0.2 jets
                //method1:most concentric jet=core 
		if(dismin<fMinDist){ AliAODJet* jetmethod1 = (AliAODJet*)(fListJets[1]->At(index1));       
      		  if(centValue<10) fh2JetCoreMethod1C10->Fill(ptcorr,jetmethod1->Pt()/ptbig);
		  if((centValue>20)&&(centValue<40)) fh2JetCoreMethod1C20->Fill(ptcorr,jetmethod1->Pt()/ptbig);
		  if((centValue>30)&&(centValue<60)) fh2JetCoreMethod1C30->Fill(ptcorr,jetmethod1->Pt()/ptbig);
		  if(centValue>60) fh2JetCoreMethod1C60->Fill(ptcorr,jetmethod1->Pt()/ptbig); }
                //method2:hardest contained jet=core   
		if(index2!=-1){ 
                  AliAODJet* jetmethod2 = (AliAODJet*)(fListJets[1]->At(index2));
                  if(centValue<10) fh2JetCoreMethod2C10->Fill(ptcorr,jetmethod2->Pt()/ptbig);
                  if((centValue>20)&&(centValue<40)) fh2JetCoreMethod2C20->Fill(ptcorr,jetmethod2->Pt()/ptbig); 
		  if((centValue>30)&&(centValue<60)) fh2JetCoreMethod2C30->Fill(ptcorr,jetmethod2->Pt()/ptbig);
		  if(centValue>60) fh2JetCoreMethod2C60->Fill(ptcorr,jetmethod2->Pt()/ptbig); }}  
		  Double_t sumpt2a=0.;
                  Double_t sumpt2b=0.;
                  Double_t sumpt4a=0.;
                  Double_t sumpt4b=0.;
                  Double_t sumpt6a=0.;
                  Double_t sumpt6b=0.;
                  Double_t sumpt8a=0.;
                  Double_t sumpt8b=0.;
                  Double_t sumpt10a=0.;
                  Double_t sumpt10b=0.;
                  Double_t sumpt2aa=0.;
                  Double_t sumpt2bb=0.;
                  Double_t sumpt4aa=0.;
                  Double_t sumpt4bb=0.;
                  Double_t sumpt6aa=0.;
                  Double_t sumpt6bb=0.;
                  Double_t sumpt8aa=0.;
                  Double_t sumpt8bb=0.;
                  Double_t sumpt10aa=0.;
                  Double_t sumpt10bb=0.;
                  Double_t sumpt2aaa=0.;
                  Double_t sumpt2bbb=0.;
                  Double_t sumpt4aaa=0.;
                  Double_t sumpt4bbb=0.;
                  Double_t sumpt6aaa=0.;
                  Double_t sumpt6bbb=0.;
                  Double_t sumpt8aaa=0.;
                  Double_t sumpt8bbb=0.;
                  Double_t sumpt10aaa=0.;
                  Double_t sumpt10bbb=0.; 
                  Double_t sumpt2ap=0.;
                  Double_t sumpt2bp=0.;
                  Double_t sumpt4ap=0.;
                  Double_t sumpt4bp=0.;
                  Double_t sumpt6ap=0.;
                  Double_t sumpt6bp=0.;
                  Double_t sumpt8ap=0.;
                  Double_t sumpt8bp=0.;
                  Double_t sumpt10ap=0.;
                  Double_t sumpt10bp=0.;
                  Double_t sumpt2aap=0.;
                  Double_t sumpt2bbp=0.;
                  Double_t sumpt4aap=0.;
                  Double_t sumpt4bbp=0.;
                  Double_t sumpt6aap=0.;
                  Double_t sumpt6bbp=0.;
                  Double_t sumpt8aap=0.;
                  Double_t sumpt8bbp=0.;
                  Double_t sumpt10aap=0.;
                  Double_t sumpt10bbp=0.;
                  Double_t sumpt2aaap=0.;
                  Double_t sumpt2bbbp=0.;
                  Double_t sumpt4aaap=0.;
                  Double_t sumpt4bbbp=0.;
                  Double_t sumpt6aaap=0.;
                  Double_t sumpt6bbbp=0.;
                  Double_t sumpt8aaap=0.;
                  Double_t sumpt8bbbp=0.;
                  Double_t sumpt10aaap=0.;
                  Double_t sumpt10bbbp=0.; 



                  

          for(int it = 0;it<nT;++it){
	  AliVParticle *part = (AliVParticle*)ParticleList.At(it);
       	  Double_t deltaR = jetbig->DeltaR(part);
          Double_t deltaEta = etabig-part->Eta();
          if(partback->Pt()>10.){

          if(centValue<10){
	    //for one more centrality
               if(deltaR<0.2){if(leadtrack->Pt()>0.)sumpt2a=sumpt2a+part->Pt();
		              if(leadtrack->Pt()>6) sumpt2b=sumpt2b+part->Pt();}
	       if(deltaR>=0.2 && deltaR<0.4){if(leadtrack->Pt()>0.)sumpt4a=sumpt4a+part->Pt();
		                               if(leadtrack->Pt()>6.)sumpt4b=sumpt4b+part->Pt();}

               if(deltaR>=0.4 && deltaR<0.6){if(leadtrack->Pt()>0.)sumpt6a=sumpt6a+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt6b=sumpt6b+part->Pt();} 
               if(deltaR>=0.6 && deltaR<0.8){if(leadtrack->Pt()>0.)sumpt8a=sumpt8a+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt8b=sumpt8b+part->Pt();}                               
               if(deltaR>=0.8 && deltaR<1.2){if(leadtrack->Pt()>0.)sumpt10a=sumpt10a+part->Pt();
		 if(leadtrack->Pt()>6.) sumpt10b=sumpt10b+part->Pt();}}

               if(centValue>30. && centValue<60.){
	    //for one more centrality
               if(deltaR<0.2){if(leadtrack->Pt()>0.)sumpt2ap=sumpt2ap+part->Pt();
		              if(leadtrack->Pt()>6) sumpt2bp=sumpt2bp+part->Pt();}
	       if(deltaR>=0.2 && deltaR<0.4){if(leadtrack->Pt()>0.)sumpt4ap=sumpt4ap+part->Pt();
		                               if(leadtrack->Pt()>6.)sumpt4bp=sumpt4bp+part->Pt();}

               if(deltaR>=0.4 && deltaR<0.6){if(leadtrack->Pt()>0.)sumpt6ap=sumpt6ap+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt6bp=sumpt6bp+part->Pt();} 
               if(deltaR>=0.6 && deltaR<0.8){if(leadtrack->Pt()>0.)sumpt8ap=sumpt8ap+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt8bp=sumpt8bp+part->Pt();}                               
               if(deltaR>=0.8 && deltaR<1.2){if(leadtrack->Pt()>0.)sumpt10ap=sumpt10ap+part->Pt();
		 if(leadtrack->Pt()>6.) sumpt10bp=sumpt10bp+part->Pt();}}}
 

                if(partback->Pt()>20.){

              if(centValue<10){
	    //for one more centrality
               if(deltaR<0.2){if(leadtrack->Pt()>0.)sumpt2aa=sumpt2aa+part->Pt();
		              if(leadtrack->Pt()>6.) sumpt2bb=sumpt2bb+part->Pt();}
	       if(deltaR>=0.2 && deltaR<0.4){if(leadtrack->Pt()>0.)sumpt4aa=sumpt4aa+part->Pt();
		                               if(leadtrack->Pt()>6.)sumpt4bb=sumpt4bb+part->Pt();}

               if(deltaR>=0.4 && deltaR<0.6){if(leadtrack->Pt()>0.)sumpt6aa=sumpt6aa+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt6bb=sumpt6bb+part->Pt();} 
               if(deltaR>=0.6 && deltaR<0.8){if(leadtrack->Pt()>0.)sumpt8aa=sumpt8aa+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt8bb=sumpt8bb+part->Pt();}                               
               if(deltaR>=0.8 && deltaR<1.2){if(leadtrack->Pt()>0.)sumpt10aa=sumpt10aa+part->Pt();
		 if(leadtrack->Pt()>6.) sumpt10bb=sumpt10bb+part->Pt();}}

               if(centValue>30. && centValue<60.){
	    //for one more centrality
               if(deltaR<0.2){if(leadtrack->Pt()>0.)sumpt2aap=sumpt2aap+part->Pt();
		              if(leadtrack->Pt()>6) sumpt2bbp=sumpt2bbp+part->Pt();}
	       if(deltaR>=0.2 && deltaR<0.4){if(leadtrack->Pt()>0.)sumpt4aap=sumpt4aap+part->Pt();
		                               if(leadtrack->Pt()>6.)sumpt4bbp=sumpt4bbp+part->Pt();}

               if(deltaR>=0.4 && deltaR<0.6){if(leadtrack->Pt()>0.)sumpt6aap=sumpt6aap+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt6bbp=sumpt6bbp+part->Pt();} 
               if(deltaR>=0.6 && deltaR<0.8){if(leadtrack->Pt()>0.)sumpt8aap=sumpt8aap+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt8bbp=sumpt8bbp+part->Pt();}                               
               if(deltaR>=0.8 && deltaR<1.2){if(leadtrack->Pt()>0.)sumpt10aap=sumpt10aap+part->Pt();
		 if(leadtrack->Pt()>6.) sumpt10bbp=sumpt10bbp+part->Pt();}}}
 



             if(partback->Pt()<1.){

              if(centValue<10){
	    //for one more centrality
               if(deltaR<0.2){if(leadtrack->Pt()>0.)sumpt2aaa=sumpt2aaa+part->Pt();
		              if(leadtrack->Pt()>6.) sumpt2bbb=sumpt2bbb+part->Pt();}
	       if(deltaR>=0.2 && deltaR<0.4){if(leadtrack->Pt()>0.)sumpt4aaa=sumpt4aaa+part->Pt();
		                               if(leadtrack->Pt()>6.)sumpt4bbb=sumpt4bbb+part->Pt();}

               if(deltaR>=0.4 && deltaR<0.6){if(leadtrack->Pt()>0.)sumpt6aaa=sumpt6aaa+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt6bbb=sumpt6bbb+part->Pt();} 
               if(deltaR>=0.6 && deltaR<0.8){if(leadtrack->Pt()>0.)sumpt8aaa=sumpt8aaa+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt8bbb=sumpt8bbb+part->Pt();}                               
               if(deltaR>=0.8 && deltaR<1.2){if(leadtrack->Pt()>0.)sumpt10aaa=sumpt10aaa+part->Pt();
		 if(leadtrack->Pt()>6.) sumpt10bbb=sumpt10bbb+part->Pt();}}

               if(centValue>30. && centValue<60.){
	    //for one more centrality
               if(deltaR<0.2){if(leadtrack->Pt()>0.)sumpt2aaap=sumpt2aaap+part->Pt();
		              if(leadtrack->Pt()>6) sumpt2bbbp=sumpt2bbbp+part->Pt();}
	       if(deltaR>=0.2 && deltaR<0.4){if(leadtrack->Pt()>0.)sumpt4aaap=sumpt4aaap+part->Pt();
		                               if(leadtrack->Pt()>6.)sumpt4bbbp=sumpt4bbbp+part->Pt();}

               if(deltaR>=0.4 && deltaR<0.6){if(leadtrack->Pt()>0.)sumpt6aaap=sumpt6aaap+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt6bbbp=sumpt6bbbp+part->Pt();} 
               if(deltaR>=0.6 && deltaR<0.8){if(leadtrack->Pt()>0.)sumpt8aaap=sumpt8aaap+part->Pt();
		                             if(leadtrack->Pt()>6.) sumpt8bbbp=sumpt8bbbp+part->Pt();}                               
               if(deltaR>=0.8 && deltaR<1.2){if(leadtrack->Pt()>0.)sumpt10aaap=sumpt10aaap+part->Pt();
		 if(leadtrack->Pt()>6.) sumpt10bbbp=sumpt10bbbp+part->Pt();}}}
 

 



          Double_t deltaPhi=phibig-part->Phi();
          if(deltaPhi<-0.5*TMath::Pi()) deltaPhi+=2.*TMath::Pi();
          if(deltaPhi>3./2.*TMath::Pi()) deltaPhi-=2.*TMath::Pi();
       	  Double_t jetEntries[8] = {centValue,ptcorr,part->Pt(),deltaR,deltaEta,deltaPhi,leadtrack->Pt(),partback->Pt()};                     fhnDeltaR->Fill(jetEntries);
          }
          //end of track loop
            Double_t rhoin2=rho*TMath::Pi()*0.2*0.2; 
            Double_t rhoin4=rho*TMath::Pi()*(0.4*0.4-0.2*0.2);
            Double_t rhoin6=rho*TMath::Pi()*(0.6*0.6-0.4*0.4);
            Double_t rhoin8=rho*TMath::Pi()*(0.8*0.8-0.6*0.6);
            Double_t rhoin10=rho*TMath::Pi()*(1.2*1.2-0.8*0.8); 
           
              
          if(rho!=0){

          if(partback->Pt()>10.){
	    if(centValue<10.){
	    if(leadtrack->Pt()>0.){
          fh2JetsumHT3R2a->Fill(ptcorr,sumpt2a/rhoin2);
          fh2JetsumHT3R4a->Fill(ptcorr,sumpt4a/rhoin4);
          fh2JetsumHT3R6a->Fill(ptcorr,sumpt6a/rhoin6);
          fh2JetsumHT3R8a->Fill(ptcorr,sumpt8a/rhoin8);
          fh2JetsumHT3R10a->Fill(ptcorr,sumpt10a/rhoin10);}
	    if(leadtrack->Pt()>6.){
          fh2JetsumHT3R2b->Fill(ptcorr,sumpt2b/rhoin2);
          fh2JetsumHT3R4b->Fill(ptcorr,sumpt4b/rhoin4);
          fh2JetsumHT3R6b->Fill(ptcorr,sumpt6b/rhoin6);
          fh2JetsumHT3R8b->Fill(ptcorr,sumpt8b/rhoin8);
          fh2JetsumHT3R10b->Fill(ptcorr,sumpt10b/rhoin10);}}

          	    if(centValue>30 && centValue<60.){
	     if(leadtrack->Pt()>0.){
          fh2JetsumHT3R2ap->Fill(ptcorr,sumpt2ap/rhoin2);
          fh2JetsumHT3R4ap->Fill(ptcorr,sumpt4ap/rhoin4);
          fh2JetsumHT3R6ap->Fill(ptcorr,sumpt6ap/rhoin6);
          fh2JetsumHT3R8ap->Fill(ptcorr,sumpt8ap/rhoin8);
          fh2JetsumHT3R10ap->Fill(ptcorr,sumpt10ap/rhoin10);}
	    if(leadtrack->Pt()>6.){
          fh2JetsumHT3R2bp->Fill(ptcorr,sumpt2bp/rhoin2);
          fh2JetsumHT3R4bp->Fill(ptcorr,sumpt4bp/rhoin4);
          fh2JetsumHT3R6bp->Fill(ptcorr,sumpt6bp/rhoin6);
          fh2JetsumHT3R8bp->Fill(ptcorr,sumpt8bp/rhoin8);
          fh2JetsumHT3R10bp->Fill(ptcorr,sumpt10bp/rhoin10);}}}





               if(partback->Pt()>20.){
	    if(centValue<10.){
	      if(leadtrack->Pt()>0.){
          fh2JetsumHT3R2aa->Fill(ptcorr,sumpt2aa/rhoin2);
          fh2JetsumHT3R4aa->Fill(ptcorr,sumpt4aa/rhoin4);
          fh2JetsumHT3R6aa->Fill(ptcorr,sumpt6aa/rhoin6);
          fh2JetsumHT3R8aa->Fill(ptcorr,sumpt8aa/rhoin8);
          fh2JetsumHT3R10aa->Fill(ptcorr,sumpt10aa/rhoin10);}
	    if(leadtrack->Pt()>6.){
          fh2JetsumHT3R2bb->Fill(ptcorr,sumpt2bb/rhoin2);
          fh2JetsumHT3R4bb->Fill(ptcorr,sumpt4bb/rhoin4);
          fh2JetsumHT3R6bb->Fill(ptcorr,sumpt6bb/rhoin6);
          fh2JetsumHT3R8bb->Fill(ptcorr,sumpt8bb/rhoin8);
          fh2JetsumHT3R10bb->Fill(ptcorr,sumpt10bb/rhoin10);}}

          	    if(centValue>30 && centValue<60.){
	     if(leadtrack->Pt()>0.){
          fh2JetsumHT3R2aap->Fill(ptcorr,sumpt2aap/rhoin2);
          fh2JetsumHT3R4aap->Fill(ptcorr,sumpt4aap/rhoin4);
          fh2JetsumHT3R6aap->Fill(ptcorr,sumpt6aap/rhoin6);
          fh2JetsumHT3R8aap->Fill(ptcorr,sumpt8aap/rhoin8);
          fh2JetsumHT3R10aap->Fill(ptcorr,sumpt10aap/rhoin10);}
	    if(leadtrack->Pt()>6.){
          fh2JetsumHT3R2bbp->Fill(ptcorr,sumpt2bbp/rhoin2);
          fh2JetsumHT3R4bbp->Fill(ptcorr,sumpt4bbp/rhoin4);
          fh2JetsumHT3R6bbp->Fill(ptcorr,sumpt6bbp/rhoin6);
          fh2JetsumHT3R8bbp->Fill(ptcorr,sumpt8bbp/rhoin8);
          fh2JetsumHT3R10bbp->Fill(ptcorr,sumpt10bbp/rhoin10);}}}



                        if(partback->Pt()<1.){
	    if(centValue<10.){
	      if(leadtrack->Pt()>0.){
          fh2JetsumHT3R2aaa->Fill(ptcorr,sumpt2aaa/rhoin2);
          fh2JetsumHT3R4aaa->Fill(ptcorr,sumpt4aaa/rhoin4);
          fh2JetsumHT3R6aaa->Fill(ptcorr,sumpt6aaa/rhoin6);
          fh2JetsumHT3R8aaa->Fill(ptcorr,sumpt8aaa/rhoin8);
          fh2JetsumHT3R10aaa->Fill(ptcorr,sumpt10aaa/rhoin10);}
	    if(leadtrack->Pt()>6.){
          fh2JetsumHT3R2bbb->Fill(ptcorr,sumpt2bbb/rhoin2);
          fh2JetsumHT3R4bbb->Fill(ptcorr,sumpt4bbb/rhoin4);
          fh2JetsumHT3R6bbb->Fill(ptcorr,sumpt6bbb/rhoin6);
          fh2JetsumHT3R8bbb->Fill(ptcorr,sumpt8bbb/rhoin8);
          fh2JetsumHT3R10bbb->Fill(ptcorr,sumpt10bbb/rhoin10);}}

          	    if(centValue>30 && centValue<60.){
	     if(leadtrack->Pt()>0.){
          fh2JetsumHT3R2aaap->Fill(ptcorr,sumpt2aaap/rhoin2);
          fh2JetsumHT3R4aaap->Fill(ptcorr,sumpt4aaap/rhoin4);
          fh2JetsumHT3R6aaap->Fill(ptcorr,sumpt6aaap/rhoin6);
          fh2JetsumHT3R8aaap->Fill(ptcorr,sumpt8aaap/rhoin8);
          fh2JetsumHT3R10aaap->Fill(ptcorr,sumpt10aaap/rhoin10);}
	    if(leadtrack->Pt()>6.){
          fh2JetsumHT3R2bbbp->Fill(ptcorr,sumpt2bbbp/rhoin2);
          fh2JetsumHT3R4bbbp->Fill(ptcorr,sumpt4bbbp/rhoin4);
          fh2JetsumHT3R6bbbp->Fill(ptcorr,sumpt6bbbp/rhoin6);
          fh2JetsumHT3R8bbbp->Fill(ptcorr,sumpt8bbbp/rhoin8);
          fh2JetsumHT3R10bbbp->Fill(ptcorr,sumpt10bbbp/rhoin10);}}}

	  }




       



   }


          //end of jet loop




          if(fDoEventMixing){
            //check before if the trigger exists
            // fTrigBuffer[i][0] = zvtx
            // fTrigBuffer[i][1] = phi
            // fTrigBuffer[i][2] = eta
            // fTrigBuffer[i][3] = pt_jet
            // fTrigBuffer[i][4] = pt_trig
            // fTrigBuffer[i][5]= pt_track_in
            // fTrigBuffer[i][6]= centrality
	    if(fTindex==11) fTindex=0;
            if(fTrigBuffer[fTindex][3]>0){
	    if (TMath::Abs(fTrigBuffer[fTindex][0]-primVtx->GetZ()<2.)){
	    if (TMath::Abs(fTrigBuffer[fTindex][6]-centValue<10)){  
	      
                        for(int it = 0;it<nT;++it){
	                AliVParticle *part = (AliVParticle*)ParticleList.At(it);         
                        Double_t DPhi = fTrigBuffer[fTindex][1] - part->Phi();
                        Double_t DEta = fTrigBuffer[fTindex][2] - part->Eta();
                        Double_t DR=TMath::Sqrt(DPhi*DPhi+DEta*DEta);
                        if(DPhi<-0.5*TMath::Pi()) DPhi+=2.*TMath::Pi();
                        if(DPhi>3./2.*TMath::Pi()) DPhi-=2.*TMath::Pi();
                        Double_t triggerEntries[8] = {centValue,fTrigBuffer[fTindex][3],part->Pt(),DR,DEta,DPhi,fTrigBuffer[fTindex][4],fTrigBuffer[fTindex][5]};                      
                        fhnMixedEvents->Fill(triggerEntries);
                        }
                        fNevents=fNevents+1;  
                        if(fNevents==9) {fTindex=fTindex+1;
			fNevents=0;} 
	    }}}
        

               // Copy the triggers from the current event into the buffer.
               //again, only if the trigger exists:
	        if(trigJet>-1){
                AliAODJet* jetT = (AliAODJet*)(fListJets[0]->At(trigJet));                 
                AliVParticle *partL = (AliVParticle*)ParticleList.At(trigInTrack);
                AliVParticle *partT = (AliVParticle*)ParticleList.At(trigBBTrack);         
                fTrigBuffer[fTrigBufferIndex][0] = primVtx->GetZ();
                fTrigBuffer[fTrigBufferIndex][1] = jetT->Phi();
                fTrigBuffer[fTrigBufferIndex][2] = jetT->Eta();
                fTrigBuffer[fTrigBufferIndex][3] = jetT->Pt()-rho*jetT->EffectiveAreaCharged();
                fTrigBuffer[fTrigBufferIndex][4] = partT->Pt();
                fTrigBuffer[fTrigBufferIndex][5] = partL->Pt();
                fTrigBuffer[fTrigBufferIndex][6] = centValue;
                fTrigBufferIndex++;
		if(fTrigBufferIndex==9) fTrigBufferIndex=0;
		}
	  }
	  




     //////////////////ANGULAR STRUCTURE//////////////////////////////////////

     //tracks up to R=0.8 distant from the jet axis
   //   if(fAngStructCloseTracks==1){
   //    TList CloseTrackList;
   //    Int_t nn=GetListOfTracksCloseToJet(&CloseTrackList,jetbig);
   //    Double_t difR=0.04;
   //    for(Int_t l=0;l<15;l++){
   //    Double_t rr=l*0.1+0.1;
   //     for(int it = 0;it<nn;++it){
   //         AliVParticle *part1 = (AliVParticle*)CloseTrackList.At(it);
   //         for(int itu=it+1;itu<CloseTrackList.GetEntries();itu++){      
   //         AliVParticle *part2 = (AliVParticle*)CloseTrackList.At(itu);  
   //         Double_t ptm=part1->Pt();
   //         Double_t ptn=part2->Pt();	
   //         Double_t Rnm = (part1->Eta()-part2->Eta())*(part1->Eta()-part2->Eta())+(part1->Phi()-part2->Phi())*(part1->Phi()-part2->Phi());
   //                    Rnm=TMath::Sqrt(Rnm);
   //                    Double_t deltag=(1./(TMath::Sqrt(2*TMath::Pi())*difR))*TMath::Exp(-1.*(rr-Rnm)*(rr-Rnm)/(2.*difR*difR));      
   //                    Double_t stepf=0.5*(1.+TMath::Erf((rr-Rnm)/(TMath::Sqrt(2.)*difR)));
   //                    if((ptcorr<85.) && (ptcorr>=70.)){up1[l]=up1[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			                                down1[l]=down1[l]+ptm*ptn*Rnm*Rnm*stepf;}
   //                    if((ptcorr<100.) && (ptcorr>=85.)){up2[l]=up2[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			                                down2[l]=down2[l]+ptm*ptn*Rnm*Rnm*stepf;}  
   //                    if((ptcorr<120.) && (ptcorr>=100.)){up3[l]=up3[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			                                 down3[l]=down3[l]+ptm*ptn*Rnm*Rnm*stepf;}
   //                    if((ptcorr<140.) && (ptcorr>=120.)){up4[l]=up4[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			down4[l]=down4[l]+ptm*ptn*Rnm*Rnm*stepf;}}}}
   //   }
    
   //   //only jet constituents
   //    if(fAngStructCloseTracks==2){

   //    Double_t difR=0.04;
   //    for(Int_t l=0;l<15;l++){
   //    Double_t rr=l*0.1+0.1;

    
   //    AliAODTrack* part1;
   //    AliAODTrack* part2;

   //        TRefArray *genTrackListb = jetbig->GetRefTracks();
   //        Int_t nTracksGenJetb = genTrackListb->GetEntriesFast();
          
             

   //        for(Int_t it=0; it<nTracksGenJetb; ++it){
   //           part1 = (AliAODTrack*)(genTrackListb->At(it));
   //         for(Int_t itu=0; itu<nTracksGenJetb; ++itu){
   //           part2 = (AliAODTrack*)(genTrackListb->At(itu));
   //         Double_t ptm=part1->Pt();
   //         Double_t ptn=part2->Pt();	
   //         Double_t Rnm = (part1->Eta()-part2->Eta())*(part1->Eta()-part2->Eta())+(part1->Phi()-part2->Phi())*(part1->Phi()-part2->Phi());
   //                    Rnm=TMath::Sqrt(Rnm);
   //                    Double_t deltag=(1./(TMath::Sqrt(2*TMath::Pi())*difR))*TMath::Exp(-1.*(rr-Rnm)*(rr-Rnm)/(2.*difR*difR));
   //                    Double_t stepf=0.5*(1.+TMath::Erf((rr-Rnm)/(TMath::Sqrt(2.)*difR))); 
   //                    if((ptcorr<85.) && (ptcorr>=70.)){up1[l]=up1[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			                                down1[l]=down1[l]+ptm*ptn*Rnm*Rnm*stepf;}
   //                    if((ptcorr<100.) && (ptcorr>=85.)){up2[l]=up2[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			                                down2[l]=down2[l]+ptm*ptn*Rnm*Rnm*stepf;}  
   //                    if((ptcorr<120.) && (ptcorr>=100.)){up3[l]=up3[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			                                 down3[l]=down3[l]+ptm*ptn*Rnm*Rnm*stepf;}
   //                    if((ptcorr<140.) && (ptcorr>=120.)){up4[l]=up4[l]+ptm*ptn*Rnm*Rnm*deltag;
   // 			down4[l]=down4[l]+ptm*ptn*Rnm*Rnm*stepf;}}}}}
   // }
     // //end loop over R=0.4 jets	
     // if(fAngStructCloseTracks>0){
     // for(Int_t l=0;l<15;l++){
     // Double_t rr=l*0.1+0.1;
     //    if(down1[l]!=0){  
     // 	if(centValue<10.)fh2AngStructpt1C10->Fill(rr,rr*up1[l]/down1[l]);
     //    if(centValue>20. && centValue<40.) fh2AngStructpt1C20->Fill(rr,rr*up1[l]/down1[l]);
     //    if(centValue>30. && centValue<60.) fh2AngStructpt1C30->Fill(rr,rr*up1[l]/down1[l]);
     //    if(centValue>60.) fh2AngStructpt1C60->Fill(rr,rr*up1[l]/down1[l]);}
     //    if(down2[l]!=0){  
     // 	if(centValue<10.) fh2AngStructpt2C10->Fill(rr,rr*up2[l]/down2[l]);
     //    if(centValue>20. && centValue<40.) fh2AngStructpt2C20->Fill(rr,rr*up2[l]/down2[l]);
     //    if(centValue>30. && centValue<60.) fh2AngStructpt2C30->Fill(rr,rr*up2[l]/down2[l]);
     //    if(centValue>60.) fh2AngStructpt2C60->Fill(rr,rr*up2[l]/down2[l]);}
     //    if(down3[l]!=0){  
     // 	if(centValue<10.) fh2AngStructpt3C10->Fill(rr,rr*up3[l]/down3[l]);
     //    if(centValue>20. && centValue<40.) fh2AngStructpt3C20->Fill(rr,rr*up3[l]/down3[l]);
     //    if(centValue>30. && centValue<60.) fh2AngStructpt3C30->Fill(rr,rr*up3[l]/down3[l]);
     //    if(centValue>60.) fh2AngStructpt3C60->Fill(rr,rr*up3[l]/down3[l]);}
     //    if(down4[l]!=0){  
     // 	if(centValue<10.) fh2AngStructpt4C10->Fill(rr,rr*up4[l]/down4[l]);
     //    if(centValue>20. && centValue<40.) fh2AngStructpt4C20->Fill(rr,rr*up4[l]/down4[l]);
     //    if(centValue>30. && centValue<60.) fh2AngStructpt4C30->Fill(rr,rr*up4[l]/down4[l]);
     //    if(centValue>60.) fh2AngStructpt4C60->Fill(rr,rr*up4[l]/down4[l]);}}}

    





   PostData(1, fOutputList);
}

void AliAnalysisTaskJetCore::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if (!GetOutputData(1))
   return;
}



  


Int_t  AliAnalysisTaskJetCore::GetListOfTracks(TList *list){

     Int_t iCount = 0;
     AliAODEvent *aod = 0;
     if(!fESD)aod = fAODIn;
     else aod = fAODOut;   
    
    for(int it = 0;it < aod->GetNumberOfTracks();++it){
      AliAODTrack *tr = aod->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      list->Add(tr);
      //cout<<fAOD->GetNumberOfTracks()<<" "<<tr->Pt()<<endl;
      iCount++;
    }
  
   
  return iCount;
 
}

   Int_t  AliAnalysisTaskJetCore::GetHardestTrackBackToJet(AliAODJet *jetbig){
 
    AliAODEvent *aod = 0;
    if(!fESD)aod = fAODIn;
    else aod = fAODOut;     
    Int_t index=-1;
    Double_t ptmax=-10;
    Double_t dphi=0;
    Double_t dif=0;
    Int_t iCount=0;
    for(int it = 0;it < aod->GetNumberOfTracks();++it){
      AliAODTrack *tr = aod->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      iCount=iCount+1;
      dphi=RelativePhi(tr->Phi(),jetbig->Phi());  
      if(TMath::Abs(dphi)<TMath::Pi()-0.2) continue;
      if(tr->Pt()>ptmax){ ptmax=tr->Pt();
      index=iCount-1;
      dif=dphi;  }}
  
      return index;

   }









 Int_t  AliAnalysisTaskJetCore::GetListOfTracksCloseToJet(TList *list,AliAODJet *jetbig){

    Int_t iCount = 0;
      AliAODEvent *aod = 0;
     if(!fESD)aod = fAODIn;
     else aod = fAODOut;   
  
    for(int it = 0;it < aod->GetNumberOfTracks();++it){
      AliAODTrack *tr = aod->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      Double_t disR=jetbig->DeltaR(tr);
      if(disR>0.8)  continue;
      list->Add(tr);
      //cout<<fAOD->GetNumberOfTracks()<<" "<<tr->Pt()<<endl;
      iCount++;
    }
  
   list->Sort();
   return iCount;

}











Int_t AliAnalysisTaskJetCore::GetNInputTracks()
{

   Int_t nInputTracks = 0;
     AliAODEvent *aod = 0;
     if(!fESD)aod = fAODIn;
     else aod = fAODOut;   
   TString jbname(fJetBranchName[1]);
   //needs complete event, use jets without background subtraction
   for(Int_t i=1; i<=3; ++i){
      if(jbname.Contains(Form("B%d",i))) jbname.ReplaceAll(Form("B%d",i),"B0");
   }
   // use only HI event
   if(jbname.Contains("AODextraonly")) jbname.ReplaceAll("AODextraonly","AOD");
   if(jbname.Contains("AODextra")) jbname.ReplaceAll("AODextra","AOD");

   if(fDebug) Printf("Multiplicity from jet branch %s", jbname.Data());
   TClonesArray *tmpAODjets = dynamic_cast<TClonesArray*>(aod->FindListObject(jbname.Data()));
   if(!tmpAODjets){
      Printf("Jet branch %s not found", jbname.Data());
      Printf("AliAnalysisTaskJetCore::GetNInputTracks FAILED");
      return -1;
   }
   
   for (Int_t iJet=0; iJet<tmpAODjets->GetEntriesFast(); iJet++){
      AliAODJet *jet = dynamic_cast<AliAODJet*>((*tmpAODjets)[iJet]);
      if(!jet) continue;
      TRefArray *trackList = jet->GetRefTracks();
      Int_t nTracks = trackList->GetEntriesFast();
      nInputTracks += nTracks;
      if(fDebug) Printf("#jet%d: %d tracks", iJet, nTracks);
   }
   if(fDebug) Printf("---> input tracks: %d", nInputTracks);

   return nInputTracks;  
}



Double_t AliAnalysisTaskJetCore::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}



THnSparse* AliAnalysisTaskJetCore::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()

   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){
      
         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }
      
      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}

void AliAnalysisTaskJetCore::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse

   const Double_t pi = TMath::Pi();
   
   switch(iEntry){
      
   case 0:
      label = "V0 centrality (%)";
     
         nbins = 10;
         xmin = 0.;
         xmax = 100.;
         break;
      
      
   case 1:
      label = "corrected jet pt";
         nbins = 20;
         xmin = 0.;
         xmax = 200.;
          break;
      
      
   case 2:
      label = "track pT";
     
         nbins = 9;
         xmin = 0.;
         xmax = 150;
         break;
      
      
    case 3:
      label = "deltaR";
      nbins = 15;
      xmin = 0.;
      xmax = 1.5;
      break;



   case 4:
      label = "deltaEta";
      nbins = 8;
      xmin = -1.6;
      xmax = 1.6;
      break;


  case 5:
      label = "deltaPhi";
      nbins = 90;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;   
   
      
        
    case 6:
      label = "leading track";
      nbins = 13;
      xmin = 0;
      xmax = 50;
      break;
           
     case 7:
    
      label = "trigger track";
      nbins =10;
      xmin = 0;
      xmax = 50;
      break;


   
  




   }

}

