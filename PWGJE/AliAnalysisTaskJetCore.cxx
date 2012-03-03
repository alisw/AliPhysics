
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
fJetEtaMin(-.5),
fJetEtaMax(.5),
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
fhnSumBkg(0x0),  
fhnJetCoreMethod3(0x0),
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
fh3spectriggered(0x0),
fh3specbiased(0x0),
fh3specleadsublead(0x0)
 
{
   // default Constructor

   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;
}

AliAnalysisTaskJetCore::AliAnalysisTaskJetCore(const char *name) :
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
fJetEtaMin(-.5),
fJetEtaMax(.5),
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
fhnSumBkg(0x0),  
fhnJetCoreMethod3(0x0),
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
fh3spectriggered(0x0),
fh3specbiased(0x0),
fh3specleadsublead(0x0)
 {
   // Constructor

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
    
     entries = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 |1<<7|1<<8; 
     fhnDeltaR = NewTHnSparseF("fhnDeltaR", entries);
     // fhnSumBkg = NewTHnSparseF("fhnDeltaR", entries);

    if(fCheckMethods){
      entries = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5  ; 
     fhnJetCoreMethod3 = NewTHnSparseF("fhnEvent", entries);

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
    fh3spectriggered = new TH3F("Triggered spectrum","",10,0,100,50,0.,200,50,0.,50.);
    fh3specbiased = new TH3F("Biased spectrum","",10,0,100,50,0.,200.,50,0.,50.);
    fh3specleadsublead = new TH3F("Leading/subleading spectrum","",10,0,100,50,0.,200.,3,0,3);

   fOutputList->Add(fHistEvtSelection);

   fOutputList->Add(fhnDeltaR);
   
   //fOutputList->Add(fhnSumBkg);
         
     
  
      if(fCheckMethods){
      fOutputList->Add(fhnJetCoreMethod3);
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
       fOutputList->Add(fh3spectriggered);
       fOutputList->Add(fh3specbiased);
       fOutputList->Add(fh3specleadsublead);

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
     PostData(1, fOutputList);
     if(fDebug) Printf(" No AOD found ... ");
     return;
   }
   AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
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
   if(fESD) cent = fESD->GetCentrality();
   if(cent) centValue = cent->GetCentralityPercentile("V0M");
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

   Double_t ptsub[aodJets[0]->GetEntriesFast()];
   Int_t inord[aodJets[0]->GetEntriesFast()];
   for(Int_t n=0;n<aodJets[0]->GetEntriesFast();n++){
     ptsub[n]=0;
     inord[n]=0;}   

   TList ParticleList;
   Int_t nT = GetListOfTracks(&ParticleList);
     for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
      fListJets[iJetType]->Clear();
      if (!aodJets[iJetType]) continue;
      if(fDebug) Printf("%s: %d jets",fJetBranchName[iJetType].Data(),aodJets[iJetType]->GetEntriesFast());
      
   
      for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
	 if(!jet)continue;
          fListJets[iJetType]->Add(jet);
         if(iJetType==0){
	   ptsub[iJet]=jet->Pt()-rho*jet->EffectiveAreaCharged();
	 }}}
   
   Double_t etabig=0;
   Double_t ptbig=0;
   Double_t areabig=0;
   Double_t phibig=0.;
   Double_t etasmall=0;
   Double_t ptsmall=0;
   Double_t areasmall=0;
   //Double_t distr=0.;
   Double_t phismall=0.;
   Int_t indexlead=-1;
   Int_t indexsublead=-1;
   Int_t indexstop=-1;
  
   if(fListJets[0]->GetEntries()>0) TMath::Sort(fListJets[0]->GetEntries(),ptsub,inord);

   for(Int_t jj=0;jj<fListJets[0]->GetEntries();jj++){
   AliAODJet* jetlead = (AliAODJet*)(fListJets[0]->At(inord[jj]));
   if(jetlead->Pt()-rho*jetlead->EffectiveAreaCharged()<=0) continue;
   if((jetlead->Eta()<fJetEtaMin)||(jetlead->Eta()>fJetEtaMax)) continue;
   indexlead=inord[jj];
   indexstop=jj;
   break;}
   if((indexstop>-1)&&(indexstop+1<fListJets[0]->GetEntries()-1)){
   for(Int_t k=indexstop+1;k<fListJets[0]->GetEntries();k++){
   AliAODJet* jetsublead = (AliAODJet*)(fListJets[0]->At(inord[k]));
   if(jetsublead->Pt()-rho*jetsublead->EffectiveAreaCharged()<=0) continue;
   if((jetsublead->Eta()<fJetEtaMin)||(jetsublead->Eta()>fJetEtaMax)) continue;
   indexsublead=inord[k];
   break;}}
         
  
   Double_t up1[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t up2[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t up3[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t up4[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down1[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down2[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down3[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t down4[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    
  

   for(Int_t i=0; i<fListJets[0]->GetEntries(); ++i){
           AliAODJet* jetbig = (AliAODJet*)(fListJets[0]->At(i));
           etabig  = jetbig->Eta();
           phibig  = jetbig->Phi();
           ptbig   = jetbig->Pt();
           if(ptbig==0) continue; 
           areabig = jetbig->EffectiveAreaCharged();
           Double_t ptcorr=ptbig-rho*areabig;
           if(ptcorr<=0) continue;
      	   if((etabig<fJetEtaMin)||(etabig>fJetEtaMax)) continue;
                   Double_t dismin=100.;
                   Double_t ptmax=-10.; 
                   Int_t index1=-1;
                   Int_t index2=-1;
                   //Double_t fracin=0.;
                  
	       Int_t point=GetHardestTrackBackToJet(jetbig);    
	       AliVParticle *partback = (AliVParticle*)ParticleList.At(point);                            
                   if(!partback) continue; 
        	   fh3spectriggered->Fill(centValue,ptcorr,partback->Pt());
                    Int_t  flaglead=0;
                       if(i==indexlead) flaglead=1;
                       if(i==indexsublead) flaglead=2;

		       fh3specleadsublead->Fill(centValue,ptcorr,flaglead);

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
	               //Float_t etr=genTrack->Eta();
                       //Float_t phir=genTrack->Phi();
                       //distr=(etr-etabig)*(etr-etabig)+(phir-phibig)*(phir-phibig);
                       //distr=TMath::Sqrt(distr);
	               //if(distr<=fRadioFrac){ fracin=fracin+genTrack->Pt();}}
                        leadtrack=(AliAODTrack*)(genTrackList->At(ippt));
                        fh3specbiased->Fill(centValue,ptcorr,leadtrack->Pt());
                       //fhnJetCoreMethod3->Fill(centValue,ptcorr,fracin/ptbig,partback->Pt(),flaglead);
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

		   
  
          for(int it = 0;it<nT;++it){
	  AliVParticle *part = (AliVParticle*)ParticleList.At(it);
     	  Double_t deltaR = jetbig->DeltaR(part);
          Double_t deltaEta = etabig-part->Eta();
          Double_t deltaPhi=phibig-part->Phi();
          if(deltaPhi<-0.5*TMath::Pi()) deltaPhi+=2.*TMath::Pi();
          if(deltaPhi>3./2.*TMath::Pi()) deltaPhi-=2.*TMath::Pi();

       	  Double_t jetEntries[9] = {centValue,ptcorr,part->Pt(),deltaR,deltaEta,deltaPhi,flaglead,leadtrack->Pt(),partback->Pt()};            fhnDeltaR->Fill(jetEntries);
	  }  
     //end of track loop
   
     //fhnSumBkg->Fill(centValue,ptcorr,bkg/jetbig->Pt(),partback->Pt(),flaglead);       
   
   
     //////////////////ANGULAR STRUCTURE//////////////////////////////////////

     //tracks up to R=0.8 distant from the jet axis
     if(fAngStructCloseTracks==1){
      TList CloseTrackList;
      Int_t nn=GetListOfTracksCloseToJet(&CloseTrackList,jetbig);
      Double_t difR=0.04;
      for(Int_t l=0;l<15;l++){
      Double_t rr=l*0.1+0.1;
       for(int it = 0;it<nn;++it){
           AliVParticle *part1 = (AliVParticle*)CloseTrackList.At(it);
           for(int itu=it+1;itu<CloseTrackList.GetEntries();itu++){      
           AliVParticle *part2 = (AliVParticle*)CloseTrackList.At(itu);  
           Double_t ptm=part1->Pt();
           Double_t ptn=part2->Pt();	
           Double_t Rnm = (part1->Eta()-part2->Eta())*(part1->Eta()-part2->Eta())+(part1->Phi()-part2->Phi())*(part1->Phi()-part2->Phi());
                      Rnm=TMath::Sqrt(Rnm);
                      Double_t deltag=(1./(TMath::Sqrt(2*TMath::Pi())*difR))*TMath::Exp(-1.*(rr-Rnm)*(rr-Rnm)/(2.*difR*difR));      
                      Double_t stepf=0.5*(1.+TMath::Erf((rr-Rnm)/(TMath::Sqrt(2.)*difR)));
                      if((ptcorr<85.) && (ptcorr>=70.)){up1[l]=up1[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down1[l]=down1[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<100.) && (ptcorr>=85.)){up2[l]=up2[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down2[l]=down2[l]+ptm*ptn*Rnm*Rnm*stepf;}  
                      if((ptcorr<120.) && (ptcorr>=100.)){up3[l]=up3[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                 down3[l]=down3[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<140.) && (ptcorr>=120.)){up4[l]=up4[l]+ptm*ptn*Rnm*Rnm*deltag;
			down4[l]=down4[l]+ptm*ptn*Rnm*Rnm*stepf;}}}}
     }
    
     //only jet constituents
      if(fAngStructCloseTracks==2){

      Double_t difR=0.04;
      for(Int_t l=0;l<15;l++){
      Double_t rr=l*0.1+0.1;

    
      AliAODTrack* part1;
      AliAODTrack* part2;

          TRefArray *genTrackListb = jetbig->GetRefTracks();
          Int_t nTracksGenJetb = genTrackListb->GetEntriesFast();
          
             

          for(Int_t it=0; it<nTracksGenJetb; ++it){
             part1 = (AliAODTrack*)(genTrackListb->At(it));
           for(Int_t itu=0; itu<nTracksGenJetb; ++itu){
             part2 = (AliAODTrack*)(genTrackListb->At(itu));
           Double_t ptm=part1->Pt();
           Double_t ptn=part2->Pt();	
           Double_t Rnm = (part1->Eta()-part2->Eta())*(part1->Eta()-part2->Eta())+(part1->Phi()-part2->Phi())*(part1->Phi()-part2->Phi());
                      Rnm=TMath::Sqrt(Rnm);
                      Double_t deltag=(1./(TMath::Sqrt(2*TMath::Pi())*difR))*TMath::Exp(-1.*(rr-Rnm)*(rr-Rnm)/(2.*difR*difR));
                      Double_t stepf=0.5*(1.+TMath::Erf((rr-Rnm)/(TMath::Sqrt(2.)*difR))); 
                      if((ptcorr<85.) && (ptcorr>=70.)){up1[l]=up1[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down1[l]=down1[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<100.) && (ptcorr>=85.)){up2[l]=up2[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                down2[l]=down2[l]+ptm*ptn*Rnm*Rnm*stepf;}  
                      if((ptcorr<120.) && (ptcorr>=100.)){up3[l]=up3[l]+ptm*ptn*Rnm*Rnm*deltag;
			                                 down3[l]=down3[l]+ptm*ptn*Rnm*Rnm*stepf;}
                      if((ptcorr<140.) && (ptcorr>=120.)){up4[l]=up4[l]+ptm*ptn*Rnm*Rnm*deltag;
			down4[l]=down4[l]+ptm*ptn*Rnm*Rnm*stepf;}}}}}


   












   }


     //end loop over R=0.4 jets	
     if(fAngStructCloseTracks>0){
     for(Int_t l=0;l<15;l++){
     Double_t rr=l*0.1+0.1;
        if(down1[l]!=0){  
	if(centValue<10.)fh2AngStructpt1C10->Fill(rr,rr*up1[l]/down1[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt1C20->Fill(rr,rr*up1[l]/down1[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt1C30->Fill(rr,rr*up1[l]/down1[l]);
        if(centValue>60.) fh2AngStructpt1C60->Fill(rr,rr*up1[l]/down1[l]);}
        if(down2[l]!=0){  
	if(centValue<10.) fh2AngStructpt2C10->Fill(rr,rr*up2[l]/down2[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt2C20->Fill(rr,rr*up2[l]/down2[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt2C30->Fill(rr,rr*up2[l]/down2[l]);
        if(centValue>60.) fh2AngStructpt2C60->Fill(rr,rr*up2[l]/down2[l]);}
        if(down3[l]!=0){  
	if(centValue<10.) fh2AngStructpt3C10->Fill(rr,rr*up3[l]/down3[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt3C20->Fill(rr,rr*up3[l]/down3[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt3C30->Fill(rr,rr*up3[l]/down3[l]);
        if(centValue>60.) fh2AngStructpt3C60->Fill(rr,rr*up3[l]/down3[l]);}
        if(down4[l]!=0){  
	if(centValue<10.) fh2AngStructpt4C10->Fill(rr,rr*up4[l]/down4[l]);
        if(centValue>20. && centValue<40.) fh2AngStructpt4C20->Fill(rr,rr*up4[l]/down4[l]);
        if(centValue>30. && centValue<60.) fh2AngStructpt4C30->Fill(rr,rr*up4[l]/down4[l]);
        if(centValue>60.) fh2AngStructpt4C60->Fill(rr,rr*up4[l]/down4[l]);}}}

    





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
 
    
    for(int it = 0;it < fAOD->GetNumberOfTracks();++it){
      AliAODTrack *tr = fAOD->GetTrack(it);
      if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
      if(TMath::Abs(tr->Eta())>0.9)continue;
      if(tr->Pt()<0.15)continue;
      list->Add(tr);
      //cout<<fAOD->GetNumberOfTracks()<<" "<<tr->Pt()<<endl;
      iCount++;
    }
  
   list->Sort();
  return iCount;
 
}

   Int_t  AliAnalysisTaskJetCore::GetHardestTrackBackToJet(AliAODJet *jetbig){

   
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
	index=iCount;
        dif=dphi;  }}
  
      return index;

   }









 Int_t  AliAnalysisTaskJetCore::GetListOfTracksCloseToJet(TList *list,AliAODJet *jetbig){

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











Int_t AliAnalysisTaskJetCore::GetNInputTracks()
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
         nbins = 50;
         xmin = 0.;
         xmax = 200.;
          break;
      
      
   case 2:
      label = "track pT";
     
      nbins = 50;
         xmin = 0.;
         xmax = 50;
         break;
      
      
   case 3:
      label = "deltaR";
      nbins = 15;
      xmin = 0.;
      xmax = 1.5;
      break;
      
   case 4:
      label = "deltaEta";
      nbins = 30;
      xmin = -1.5;
      xmax = 1.5;
      break;


  case 5:
      label = "deltaPhi";
      nbins = 90;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;   
   
      
   case 6:
      label="flagleadname";
      nbins=3;
      xmin=0;
      xmax=3;
      break;

    case 7:
    
      label = "leading track";
      nbins = 50;
      xmin = 0;
      xmax = 50;
      break;
           
     case 8:
    
      label = "trigger track";
      nbins =50;
      xmin = 0;
      xmax = 50;
      break;
      
        
   }

}

