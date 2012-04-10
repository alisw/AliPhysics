
// ******************************************
// This task searches for events with large dijet imbalance
// and then looks to the jet structure of the b-t-b jets.    
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
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

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

#include "AliAnalysisTaskAj.h"

ClassImp(AliAnalysisTaskAj)

AliAnalysisTaskAj::AliAnalysisTaskAj() :
AliAnalysisTaskSE(),
fESD(0x0),
fAOD(0x0),
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
fh2Pt1Pt2C10(0x0),
fh2Pt1Pt2C40(0x0), 
fh3LocalCoordinates(0x0), 
fh2Sum2pt20(0x0),
fh2Sum4pt20(0x0),
fh2Sum6pt20(0x0),
fh2Sum8pt20(0x0),
fh2Sum12pt20(0x0),
fh2Sum2lpt20(0x0),
fh2Sum4lpt20(0x0),
fh2Sum6lpt20(0x0),
fh2Sum8lpt20(0x0),
fh2Sum12lpt20(0x0),
fh2Sum2pt40(0x0),
fh2Sum4pt40(0x0),
fh2Sum6pt40(0x0),
fh2Sum8pt40(0x0),
fh2Sum12pt40(0x0),
fh2Sum2lpt40(0x0),
fh2Sum4lpt40(0x0),
fh2Sum6lpt40(0x0),
fh2Sum8lpt40(0x0),
fh2Sum12lpt40(0x0),
fhnDeltaR(0x0)
 {
   // default Constructor


   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;
}

AliAnalysisTaskAj::AliAnalysisTaskAj(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAOD(0x0),
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
fh2Pt1Pt2C10(0x0),
fh2Pt1Pt2C40(0x0), 
fh3LocalCoordinates(0x0), 
fh2Sum2pt20(0x0),
fh2Sum4pt20(0x0),
fh2Sum6pt20(0x0),
fh2Sum8pt20(0x0),
fh2Sum12pt20(0x0),
fh2Sum2lpt20(0x0),
fh2Sum4lpt20(0x0),
fh2Sum6lpt20(0x0),
fh2Sum8lpt20(0x0),
fh2Sum12lpt20(0x0),
fh2Sum2pt40(0x0),
fh2Sum4pt40(0x0),
fh2Sum6pt40(0x0),
fh2Sum8pt40(0x0),
fh2Sum12pt40(0x0),
fh2Sum2lpt40(0x0),
fh2Sum4lpt40(0x0),
fh2Sum6lpt40(0x0),
fh2Sum8lpt40(0x0),
fh2Sum12lpt40(0x0),
fhnDeltaR(0x0)
 {
   // Constructor


   




   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;

   DefineOutput(1, TList::Class());
}

AliAnalysisTaskAj::~AliAnalysisTaskAj()
{
   delete fListJets[0];
   delete fListJets[1];
}

void AliAnalysisTaskAj::SetBranchNames(const TString &branch1, const TString &branch2)
{
   fJetBranchName[0] = branch1;
   fJetBranchName[1] = branch2;
}

void AliAnalysisTaskAj::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskAj::UserCreateOutputObjects()
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
     entries = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 |1<<7|1<<8; 
     fhnDeltaR = NewTHnSparseF("fhnDeltaR", entries);

    
     //change binning in centrality
     Double_t *xPt0 = new Double_t[4];
     xPt0[0] = 0.;
     for(int i = 1; i<=3;i++){
      if(xPt0[i-1]<10)xPt0[i] = xPt0[i-1] + 10; // 1 - 5
      else if(xPt0[i-1]<40)xPt0[i] = xPt0[i-1] + 30; // 5 - 12
      else xPt0[i] = xPt0[i-1] + 60.; // 18 
     }
     fhnDeltaR->SetBinEdges(0,xPt0);
     delete [] xPt0;



     //change binning in jTtrack
     Double_t *xPt3 = new Double_t[3];
     xPt3[0] = 0.;
     for(int i = 1; i<=2;i++){
      if(xPt3[i-1]<1)xPt3[i] = xPt3[i-1] + 1.; // 1 - 5
      else xPt3[i] = xPt3[i-1] + 149.; // 18 
     }
     fhnDeltaR->SetBinEdges(3,xPt3);
     delete [] xPt3;



     //change binning in pTtrack
     Double_t *xPt4 = new Double_t[5];
     xPt4[0] = 0.;
     for(int i = 1; i<=4;i++){
      if(xPt4[i-1]<0.4)xPt4[i] = xPt4[i-1] + 0.4; // 1 - 5
      else if(xPt4[i-1]<3)xPt4[i] = xPt4[i-1] + 2.6; // 5 - 12
      else if(xPt4[i-1]<10)xPt4[i] = xPt4[i-1] + 7.4; // 5 - 12
      else xPt4[i] = xPt4[i-1] + 150.; // 18 
     }
    fhnDeltaR->SetBinEdges(4,xPt4);
    delete [] xPt4;






    //change binning in HTI
     Double_t *xPt5 = new Double_t[4];
     xPt5[0] = 0.;
     for(int i = 1; i<=3;i++){
      if(xPt5[i-1]<10)xPt5[i] = xPt5[i-1] + 5; // 10 - 12
      else xPt5[i] = xPt5[i-1] + 40.; // 13
     }
    fhnDeltaR->SetBinEdges(8,xPt5);
    delete [] xPt5;

    
   
    fh2Pt1Pt2C10 = new TH2F("Dijet spectra central","",20,0.,200.,20,0.,200.);
    fh2Pt1Pt2C40 = new TH2F("Dijet spectra peripheral","",20,0.,200.,20,0.,200.);
    fh3LocalCoordinates = new TH3F("Local coordinates","",10,-2,2,10,-2,2,10,0,100);
    fh2Sum2pt20 = new TH2F("pL R<0.2 pt20","",10,0.,1.,100,0.,100.);
    fh2Sum4pt20 = new TH2F("pL R<0.4 pt20","",10,0.,1.,100,0.,100.);
    fh2Sum6pt20 = new TH2F("pL R<0.6 pt20","",10,0.,1.,100,0.,100.);
    fh2Sum8pt20 = new TH2F("pL R<0.8 pt20","",10,0.,1.,100,0.,100.);
    fh2Sum12pt20 = new TH2F("pL R<1.2 pt20","",10,0.,1.,100,0.,100.);
    fh2Sum2lpt20 = new TH2F("pL R<0.2 low pt pt20","",10,0.,1.,100,0,100);
    fh2Sum4lpt20 = new TH2F("pL R<0.4 low pt pt20","",10,0.,1.,100,0,100);
    fh2Sum6lpt20 = new TH2F("pL R<0.6 low pt pt20","",10,0.,1.,100,0,100);
    fh2Sum8lpt20 = new TH2F("pL R<0.8 low pt pt20","",10,0.,1.,100,0,100);
    fh2Sum12lpt20 = new TH2F("pL R<1.2 low pt pt20","",10,0.,1.,100,0,100);
    fh2Sum2pt40 = new TH2F("pL R<0.2 pt40","",10,0.,1.,100,0.,100.);
    fh2Sum4pt40 = new TH2F("pL R<0.4 pt40","",10,0.,1.,100,0.,100.);
    fh2Sum6pt40 = new TH2F("pL R<0.6 pt40","",10,0.,1.,100,0.,100.);
    fh2Sum8pt40 = new TH2F("pL R<0.8 pt40","",10,0.,1.,100,0.,100.);
    fh2Sum12pt40 = new TH2F("pL R<1.2 pt40","",10,0.,1.,100,0.,100.);
    fh2Sum2lpt40 = new TH2F("pL R<0.2 low pt pt40","",10,0.,1.,100,0,100);
    fh2Sum4lpt40 = new TH2F("pL R<0.4 low pt pt40","",10,0.,1.,100,0,100);
    fh2Sum6lpt40 = new TH2F("pL R<0.6 low pt pt40","",10,0.,1.,100,0,100);
    fh2Sum8lpt40 = new TH2F("pL R<0.8 low pt pt40","",10,0.,1.,100,0,100);
    fh2Sum12lpt40 = new TH2F("pL R<1.2 low pt pt40","",10,0.,1.,100,0,100);


   fOutputList->Add(fHistEvtSelection);
   fOutputList->Add(fh2Pt1Pt2C10);
   fOutputList->Add(fh2Pt1Pt2C40);
   fOutputList->Add(fh3LocalCoordinates);

   fOutputList->Add(fh2Sum2pt20);
   fOutputList->Add(fh2Sum4pt20);
   fOutputList->Add(fh2Sum6pt20);      
   fOutputList->Add(fh2Sum8pt20);
   fOutputList->Add(fh2Sum12pt20);  
   fOutputList->Add(fh2Sum2lpt20);
   fOutputList->Add(fh2Sum4lpt20);
   fOutputList->Add(fh2Sum6lpt20);      
   fOutputList->Add(fh2Sum8lpt20);
   fOutputList->Add(fh2Sum12lpt20);  

   fOutputList->Add(fh2Sum2pt40);
   fOutputList->Add(fh2Sum4pt40);
   fOutputList->Add(fh2Sum6pt40);      
   fOutputList->Add(fh2Sum8pt40);
   fOutputList->Add(fh2Sum12pt40);  
   fOutputList->Add(fh2Sum2lpt40);
   fOutputList->Add(fh2Sum4lpt40);
   fOutputList->Add(fh2Sum6lpt40);      
   fOutputList->Add(fh2Sum8lpt40);
   fOutputList->Add(fh2Sum12lpt40);  

   fOutputList->Add(fhnDeltaR);
        
      
     
   // fOutputList->Add(fh3specbiased);

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

void AliAnalysisTaskAj::UserExec(Option_t *)
{
   

   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
      return;
   }

   fESD=dynamic_cast<AliESDEvent*>(InputEvent());
   if (!fESD) {
      AliError("ESD not available");
      fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
   } else {
      fAOD = dynamic_cast<AliAODEvent*>(AODEvent());
   }
 
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
   if(!fAOD){
     if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
     fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
   }
   AliAODVertex* primVtx = fAOD->GetPrimaryVertex();

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
   else     centValue=fAOD->GetHeader()->GetCentrality();
   
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
   if(fAOD&&!externalBackground&&fBackgroundBranch.Length()){
      externalBackground =  (AliAODJetEventBackground*)(fAOD->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   if(fAODExtension&&!externalBackground&&fBackgroundBranch.Length()){
     externalBackground =  (AliAODJetEventBackground*)(fAODExtension->GetAOD()->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   
   Float_t rho = 0;
   if(externalBackground)rho = externalBackground->GetBackground(0);


   // fetch jets
   TClonesArray *aodJets[2];
   aodJets[0]=0;
   if(fAOD&&!aodJets[0]){
   aodJets[0] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[0].Data())); 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[1].Data()));  }
   if(fAODExtension && !aodJets[0]){ 
   aodJets[0] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[0].Data())); 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName[1].Data()));  }

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
   

   Double_t eta1=0;
   Int_t selec=-1;
   Double_t ptmax=-10; 
   Double_t areaj=0;
   Double_t phij=0;
   Double_t etaj=0;
   Double_t ptj=0;
   Double_t ptcorrj=0;  
   for(Int_t i=0; i<fListJets[0]->GetEntries(); ++i){
           AliAODJet* jetj = (AliAODJet*)(fListJets[0]->At(i));
           etaj  = jetj->Eta();
           phij  = jetj->Phi();
           ptj = jetj->Pt();
           if(ptj==0) continue; 
           areaj = jetj->EffectiveAreaCharged();
           ptcorrj=ptj-rho*areaj;
           if(ptcorrj<=0) continue;
      	   if((etaj<fJetEtaMin)||(eta1>fJetEtaMax)) continue;
	   if(ptcorrj>ptmax){ptmax=ptcorrj;
	                     selec=i;}}
   ///hardest jet selected
      if(selec<0){PostData(1, fOutputList);
			  return;} 
    AliAODJet* jet1 = (AliAODJet*)(fListJets[0]->At(selec));
    //What is the hardest constituent track?
                       AliAODTrack* leadtrack1; 
                       Int_t ippt=0;
                       Double_t ppt=-10;   
		       TRefArray *genTrackList = jet1->GetRefTracks();
                       Int_t nTracksGenJet = genTrackList->GetEntriesFast();
                       AliAODTrack* genTrack;
                       for(Int_t ir=0; ir<nTracksGenJet; ++ir){
                       genTrack = (AliAODTrack*)(genTrackList->At(ir));
		       if(genTrack->Pt()>ppt){ppt=genTrack->Pt();
			 ippt=ir;}}
                        leadtrack1=(AliAODTrack*)(genTrackList->At(ippt));
			//If it doesn't exist or if it is greater that 100 GeV, discard.
                        if(!leadtrack1)  {PostData(1, fOutputList);
			  return;}
                        if(leadtrack1->Pt()>=100.){ PostData(1, fOutputList);
			  return;}             

	   //Look to the back-to-back jet 
	   Int_t btb=-1;
           for(Int_t j=1;j<fListJets[0]->GetEntries();j++){
	   if(j==selec) continue;
           AliAODJet* jetb = (AliAODJet*)(fListJets[0]->At(j));
           etaj  = jetb->Eta();
           phij  = jetb->Phi();
           ptj   = jetb->Pt();
           if(ptj<=0) continue; 
           areaj = jetb->EffectiveAreaCharged();
           ptcorrj=ptj-rho*areaj;
           if(ptcorrj<=0) continue;
      	   if((etaj<fJetEtaMin)||(etaj>fJetEtaMax)) continue;
           Double_t dphij=RelativePhi(jetb->Phi(),jet1->Phi());  
           if(TMath::Abs(dphij)>TMath::Pi()-0.2) { btb=j;
	     break;}}

            AliAODJet* jet2 = (AliAODJet*)(fListJets[0]->At(btb));
           //the back-to-back jet is also identified
           
           if(btb<0){PostData(1, fOutputList);
			  return;} 

            Double_t ptcorr1=jet1->Pt()-rho*jet1->EffectiveAreaCharged();
            Double_t ptcorr2=jet2->Pt()-rho*jet2->EffectiveAreaCharged();
          
	    if(centValue<10.) fh2Pt1Pt2C10->Fill(ptcorr1,ptcorr2);
            if(centValue>40.) fh2Pt1Pt2C40->Fill(ptcorr1,ptcorr2);
 
              
           Double_t px2=jet2->Px();
           Double_t py2=jet2->Py();
           Double_t pz2=jet2->Pz();  
           Double_t phi2=jet2->Phi();
           Double_t eta2=jet2->Eta();
	   //Once we have have a dijet event,look to the structure of the back-to-back jet:

      TVector3  ppJ1(px2, py2, pz2);
      TVector3  ppJ3(- px2 * pz2, - py2 * pz2, px2 * px2 + py2 * py2);
      ppJ3.SetMag(1.);
      TVector3  ppJ2(-py2, px2, 0);
      ppJ2.SetMag(1.);
      Float_t mxx    = 0.;
      Float_t myy    = 0.;
      Float_t mxy    = 0.;
      Int_t   nc     = 0;
      Float_t sump2  = 0.;
      Float_t ptMax  = 0.;
      Float_t etaMax = 0.;
      Float_t phiMax = 0.;
      Int_t   iMax   = -1;
      //1st loop over all tracks      
         for(int it = 0;it<nT;++it){
         AliVParticle *track = (AliVParticle*)ParticleList.At(it);
         TVector3 pp(track->Px(), track->Py(), track->Pz());
      	  Float_t phi = track->Phi();
	  Float_t eta = track->Eta();
	  Float_t pt  = track->Pt();
	  Float_t jT  = pp.Perp(ppJ1);
          Float_t deta = eta - eta2;
	  Float_t dphi = phi - phi2;
	  if (dphi >   TMath::Pi()) dphi =   2. * TMath::Pi() - dphi;
	  if (dphi < - TMath::Pi()) dphi = - 2. * TMath::Pi() - dphi;
	  Float_t r = TMath::Sqrt(dphi * dphi + deta * deta);
	  /////////To compute the TM axis, we use particles with large jT
         
              //cout<<"Jt spectrum............"<<ptbig<<" "<<jT<<endl;
	  /////////We compute the TM with large jT particles only
	     if((jT >1.)&&(r<1.2)){
	    //longitudinal and perpendicular component of the track pT in the
            //local frame
	      TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
	      TVector3 pPerp = pp - pLong;
	      //projection onto the two perpendicular vectors defined above
	      Float_t ppjX = pPerp.Dot(ppJ2);
	      Float_t ppjY = pPerp.Dot(ppJ3);
	      Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
              //components of the sphericity matrix
	      mxx += (ppjX * ppjX / ppjT);
	      myy += (ppjY * ppjY / ppjT);
	      mxy += (ppjX * ppjY / ppjT);
	      nc++;
	      sump2 += ppjT;
	      // max pt
	      if (pt > ptMax) {
		  ptMax  = pt;
		  iMax   = it;
		  etaMax = deta;
		  phiMax = dphi;
	      }
		} // R < 0.4
	 } // 1st Track Loop
	    
  // At this point we have mxx, myy, mxy
      if (nc == 0) return;      
// Shericity Matrix	
      const Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};	
      TMatrixDSym m0(2,ele);
// Find eigenvectors
      TMatrixDSymEigen m(m0);
      TVectorD eval(2);
      TMatrixD evecm = m.GetEigenVectors();
      eval  = m.GetEigenValues();
// Largest eigenvector
      Int_t jev = 0;
      if (eval[0] < eval[1]) jev = 1;
      TVectorD evec0(2);
// Principle axis
      evec0 = TMatrixDColumn(evecm, jev);
      TVector2 evec(evec0[0], evec0[1]); 
// Principle axis from leading partice
      Float_t phiM = TMath::ATan2(phiMax, etaMax);
      TVector2 evecM(TMath::Cos(phiM), TMath::Sin(phiM)); 
      Float_t phistM = evecM.DeltaPhi(evec);
      if (TMath::Abs(phistM) > TMath::Pi()/2.) evec*=(-1.);
	   
      //////we have now the direction 
      /////along which the sum of the projections of the particle
      ///momentum is higher.
      Double_t sum2lpt20=0;
      Double_t sum4lpt20=0;
      Double_t sum6lpt20=0;
      Double_t sum8lpt20=0;
      Double_t sum12lpt20=0;
      Double_t sum2pt20=0;
      Double_t sum4pt20=0;
      Double_t sum6pt20=0;
      Double_t sum8pt20=0;
      Double_t sum12pt20=0;

       Double_t sum2lpt40=0;
      Double_t sum4lpt40=0;
      Double_t sum6lpt40=0;
      Double_t sum8lpt40=0;
      Double_t sum12lpt40=0;
      Double_t sum2pt40=0;
      Double_t sum4pt40=0;
      Double_t sum6pt40=0;
      Double_t sum8pt40=0;
      Double_t sum12pt40=0;

	  for (Int_t ip = 0; ip < nT; ip++) {
	  AliVParticle *track = (AliVParticle*)ParticleList.At(ip);
	  TVector3 pp(track->Px(), track->Py(), track->Pz());
	  Float_t phi = track->Phi();
	  Float_t eta = track->Eta();
	  Float_t pt  = track->Pt();
	  Float_t jT  = pp.Perp(ppJ1);
          Double_t deta = eta - eta2;
	  Double_t deltaPhi = phi - phi2;
          Double_t r = TMath::Sqrt(deltaPhi * deltaPhi + deta * deta);
          
            if(ptcorr2>20. && ptcorr2<40.){ 
	    if(pt<0.4){
                 if(r<0.2) sum2lpt20=sum2lpt20-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.4) sum4lpt20=sum4lpt20-1.*pt*TMath::Cos(phi-jet1->Phi());   
                 if(r<0.6) sum6lpt20=sum6lpt20-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.8) sum8lpt20=sum8lpt20-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<1.2) sum12lpt20=sum12lpt20-1.*pt*TMath::Cos(phi-jet1->Phi());}
            
                 if(r<0.2) sum2pt20=sum2pt20-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.4) sum4pt20=sum4pt20-1.*pt*TMath::Cos(phi-jet1->Phi());   
                 if(r<0.6) sum6pt20=sum6pt20-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.8) sum8pt20=sum8pt20-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<1.2) sum12pt20=sum12pt20-1.*pt*TMath::Cos(phi-jet1->Phi());}


                if(ptcorr2>40. && ptcorr2<60.){ 
	    if(pt<0.4){
                 if(r<0.2) sum2lpt40=sum2lpt40-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.4) sum4lpt40=sum4lpt40-1.*pt*TMath::Cos(phi-jet1->Phi());   
                 if(r<0.6) sum6lpt40=sum6lpt40-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.8) sum8lpt40=sum8lpt40-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<1.2) sum12lpt40=sum12lpt40-1.*pt*TMath::Cos(phi-jet1->Phi());}
            
                 if(r<0.2) sum2pt40=sum2pt40-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.4) sum4pt40=sum4pt40-1.*pt*TMath::Cos(phi-jet1->Phi());   
                 if(r<0.6) sum6pt40=sum6pt40-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<0.8) sum8pt40=sum8pt40-1.*pt*TMath::Cos(phi-jet1->Phi());
                 if(r<1.2) sum12pt40=sum12pt40-1.*pt*TMath::Cos(phi-jet1->Phi());}
	 


          
          TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
	  TVector3 pPerp = pp - pLong;
	  Float_t ppjX = pPerp.Dot(ppJ2);
	  Float_t ppjY = pPerp.Dot(ppJ3);
       	  TVector2 vr(ppjX, ppjY) ;
          //and this is the angle between the particle and the TM axis. 
	  //	  Float_t phistr = evec.DeltaPhi(vr);

          Double_t phistr=vr.Phi()-evec.Phi();

	  if(centValue<10.) fh3LocalCoordinates->Fill(ppjX,ppjY,ptcorr2); 
          Double_t deltaEta = eta2-track->Eta();
          if(phistr<-0.5*TMath::Pi()) phistr+=2.*TMath::Pi();
          if(phistr>3./2.*TMath::Pi()) phistr-=2.*TMath::Pi();

          if(deltaPhi<-0.5*TMath::Pi()) deltaPhi+=2.*TMath::Pi();
          if(deltaPhi>3./2.*TMath::Pi()) deltaPhi-=2.*TMath::Pi();
       	  Double_t jetEntries[9] = {centValue,ptcorr1,ptcorr2,jT,pt,deltaEta,deltaPhi,phistr,ptMax}; 
          fhnDeltaR->Fill(jetEntries);
	  }
          if(centValue<10.){
	    if(ptcorr2>20.){
	      if(ptcorr1>80.){
          Double_t aj=(ptcorr1-ptcorr2)/(ptcorr1+ptcorr2); 
	  
	  fh2Sum2pt20->Fill(aj,sum2pt20);
	  fh2Sum4pt20->Fill(aj,sum4pt20);
          fh2Sum6pt20->Fill(aj,sum6pt20);
          fh2Sum8pt20->Fill(aj,sum8pt20);
          fh2Sum12pt20->Fill(aj,sum12pt20);  
          fh2Sum2lpt20->Fill(aj,sum2lpt20);
          fh2Sum4lpt20->Fill(aj,sum4lpt20);
          fh2Sum6lpt20->Fill(aj,sum6lpt20);
          fh2Sum8lpt20->Fill(aj,sum8lpt20);
          fh2Sum12lpt20->Fill(aj,sum12lpt20);
	  
          fh2Sum2pt40->Fill(aj,sum2pt40);
          fh2Sum4pt40->Fill(aj,sum4pt40);
          fh2Sum6pt40->Fill(aj,sum6pt40);
          fh2Sum8pt40->Fill(aj,sum8pt40);
          fh2Sum12pt40->Fill(aj,sum12pt40);  
          fh2Sum2lpt40->Fill(aj,sum2lpt40);
          fh2Sum4lpt40->Fill(aj,sum4lpt40);
          fh2Sum6lpt40->Fill(aj,sum6lpt40);
          fh2Sum8lpt40->Fill(aj,sum8lpt40);
          fh2Sum12lpt40->Fill(aj,sum12lpt40);
	      }}}







   PostData(1, fOutputList);
}

void AliAnalysisTaskAj::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if (!GetOutputData(1))
   return;
}











Int_t  AliAnalysisTaskAj::GetListOfTracks(TList *list){

    Int_t iCount = 0;
 
    
    for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      list->Add(tr);
      //cout<<fAOD->GetNumberOfTracks()<<" "<<tr->Pt()<<endl;
      iCount++;
    }
  
   
  return iCount;
 
}

   Int_t  AliAnalysisTaskAj::GetHardestTrackBackToJet(AliAODJet *jetbig){

   
    Int_t index=-1;
    Double_t ptmax=-10;
    Double_t dphi=0;
    Double_t dif=0;
    Int_t iCount=0;
    for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
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









 Int_t  AliAnalysisTaskAj::GetListOfTracksCloseToJet(TList *list,AliAODJet *jetbig){

    Int_t iCount = 0;
 
  
    for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
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











Int_t AliAnalysisTaskAj::GetNInputTracks()
{

   Int_t nInputTracks = 0;

   TString jbname(fJetBranchName[1]);
   //needs complete event, use jets without background subtraction
   for(Int_t i=1; i<=3; ++i){
      if(jbname.Contains(Form("B%d",i))) jbname.ReplaceAll(Form("B%d",i),"B0");
   }
   // use only HI event
   if(jbname.Contains("AODextraonly")) jbname.ReplaceAll("AODextraonly","AOD");
   if(jbname.Contains("AODextra")) jbname.ReplaceAll("AODextra","AOD");

   if(fDebug) Printf("Multiplicity from jet branch %s", jbname.Data());
   TClonesArray *tmpAODjets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(jbname.Data()));
   if(!tmpAODjets){
      Printf("Jet branch %s not found", jbname.Data());
      Printf("AliAnalysisTaskAj::GetNInputTracks FAILED");
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



Double_t AliAnalysisTaskAj::RelativePhi(Double_t mphi,Double_t vphi){

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}



THnSparse* AliAnalysisTaskAj::NewTHnSparseF(const char* name, UInt_t entries)
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

void AliAnalysisTaskAj::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse

   const Double_t pi = TMath::Pi();
   
   switch(iEntry){
      
   case 0:
      label = "V0 centrality (%)";
     
         nbins = 3;
         xmin = 0.;
         xmax = 100.;
         break;
      
      
   case 1:
      label = "corrected jet pt1";
         nbins = 10;
         xmin = 0.;
         xmax = 200.;
          break;
      
      case 2:
      label = "corrected jet pt2";
         nbins = 10;
         xmin = 0.;
         xmax = 200.;
          break;
      

      
   case 3:
      label = "track jT";
     
         nbins = 2;
         xmin = 0.;
         xmax = 150.;
         break;

   case 4:
      label = "track pT";
     
         nbins = 4;
         xmin = 0.;
         xmax = 150.;
         break;
      
      
      case 5:
      label = "deltaEta";
      nbins = 8;
      xmin = -1.6;
      xmax = 1.6;
      break;


      case 6:
      label = "deltaPhi";
      nbins = 60;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;   
   
      case 7:
      label = "deltaPhiTM";
      nbins = 60;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;   

     
        
      case 8:
      label = "leading track";
      nbins = 3;
      xmin = 0;
      xmax = 50;
      break;
           



   
  




   }

}

