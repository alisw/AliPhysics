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
#include "AliAODJet.h"

#include "AliAnalysisTaskJetCore.h"

ClassImp(AliAnalysisTaskJetCore)

AliAnalysisTaskJetCore::AliAnalysisTaskJetCore() :
AliAnalysisTaskSE(),
fESD(0x0),
fAOD(0x0),
fBackgroundBranch(""),
fIsPbPb(kTRUE),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-8.),
fVtxZMax(8.),
fEvtClassMin(0),
fEvtClassMax(4),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fRadioFrac(0.2),
fMinDist(0.1),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fHistEvtSelection(0x0),
fHistJetSelection(0x0),
fh2JetSelection(0x0),
fh2JetCoreMethod1C10(0x0),
fh2JetCoreMethod2C10(0x0),
fh2JetCoreMethod3C10(0x0),
fh2JetCoreMethod1C30(0x0),
fh2JetCoreMethod2C30(0x0),
fh2JetCoreMethod3C30(0x0),
fh2JetCoreMethod1C60(0x0),
fh2JetCoreMethod2C60(0x0),
fh2JetCoreMethod3C60(0x0)

 
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
fBackgroundBranch(""),
fIsPbPb(kTRUE),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-8.),
fVtxZMax(8.),
fEvtClassMin(0),
fEvtClassMax(4),
fCentMin(0.),
fCentMax(100.),
fNInputTracksMin(0),
fNInputTracksMax(-1),
fJetEtaMin(-.5),
fJetEtaMax(.5),
fJetPtMin(20.),
fJetTriggerExcludeMask(AliAODJet::kHighTrackPtTriggered),
fJetPtFractionMin(0.5),
fNMatchJets(4),
fMatchMaxDist(0.8),
fKeepJets(kFALSE),
fRadioFrac(0.2),
fMinDist(0.1),
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fHistEvtSelection(0x0),
fHistJetSelection(0x0),
fh2JetSelection(0x0),
fh2JetCoreMethod1C10(0x0),
fh2JetCoreMethod2C10(0x0),
fh2JetCoreMethod3C10(0x0),
fh2JetCoreMethod1C30(0x0),
fh2JetCoreMethod2C30(0x0),
fh2JetCoreMethod3C30(0x0),
fh2JetCoreMethod1C60(0x0),
fh2JetCoreMethod2C60(0x0),
fh2JetCoreMethod3C60(0x0)


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

   fHistJetSelection = new TH1I("fHistJetSelection", "jet selection", 8, -0.5, 7.5);
   fHistJetSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistJetSelection->GetXaxis()->SetBinLabel(2,"probes IN");
   fHistJetSelection->GetXaxis()->SetBinLabel(3,"no matching jet");
   fHistJetSelection->GetXaxis()->SetBinLabel(4,"not in list");
   fHistJetSelection->GetXaxis()->SetBinLabel(5,"fraction cut");
   fHistJetSelection->GetXaxis()->SetBinLabel(6,"acceptance cut");
   fHistJetSelection->GetXaxis()->SetBinLabel(7,"p_{T} cut");
   fHistJetSelection->GetXaxis()->SetBinLabel(8,"trigger exclude mask");

   fh2JetSelection = new TH2F("fh2JetSelection", "jet selection", 8, -0.5, 7.5,100,0.,200.);
   fh2JetSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fh2JetSelection->GetXaxis()->SetBinLabel(2,"probes IN");
   fh2JetSelection->GetXaxis()->SetBinLabel(3,"no matching jet");
   fh2JetSelection->GetXaxis()->SetBinLabel(4,"not in list");
   fh2JetSelection->GetXaxis()->SetBinLabel(5,"fraction cut");
   fh2JetSelection->GetXaxis()->SetBinLabel(6,"acceptance cut");
   fh2JetSelection->GetXaxis()->SetBinLabel(7,"p_{T} cut");
   fh2JetSelection->GetXaxis()->SetBinLabel(8,"trigger exclude mask");


   //   UInt_t entries = 0; // bit coded, see GetDimParams() below
   //   UInt_t opt = 0;  // bit coded, default (0) or high resolution (1)

    fh2JetCoreMethod1C10 = new TH2F("JetCoreMethod1C10","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C10 = new TH2F("JetCoreMethod2C10","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C10 = new TH2F("JetCoreMethod3C10","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod1C30 = new TH2F("JetCoreMethod1C30","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C30 = new TH2F("JetCoreMethod2C30","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C30 = new TH2F("JetCoreMethod3C30","",150, 0., 150.,100, 0., 1.5); 
    fh2JetCoreMethod1C60 = new TH2F("JetCoreMethod1C60","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod2C60 = new TH2F("JetCoreMethod2C60","",150, 0., 150.,100, 0., 1.5);
    fh2JetCoreMethod3C60 = new TH2F("JetCoreMethod3C60","",150, 0., 150.,100, 0., 1.5); 


   fOutputList->Add(fHistEvtSelection);
   fOutputList->Add(fHistJetSelection);
   fOutputList->Add(fh2JetSelection);
  
  
      fOutputList->Add(fh2JetCoreMethod1C10);
      fOutputList->Add(fh2JetCoreMethod2C10);
      fOutputList->Add(fh2JetCoreMethod3C10);
      fOutputList->Add(fh2JetCoreMethod1C30);
      fOutputList->Add(fh2JetCoreMethod2C30);
      fOutputList->Add(fh2JetCoreMethod3C30);
      fOutputList->Add(fh2JetCoreMethod1C60);
      fOutputList->Add(fh2JetCoreMethod2C60);
      fOutputList->Add(fh2JetCoreMethod3C60);
     
     

   // =========== Switch on Sumw2 for all histos ===========
   for (Int_t i=0; i<fOutputList->GetEntries(); ++i) {
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if (h1){
         h1->Sumw2();
         continue;
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
   if (!fAOD) {
      AliError("AOD not available");
      return;
   }

   // -- event selection --
   fHistEvtSelection->Fill(1); // number of events before event selection

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
   ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      fHistEvtSelection->Fill(2);
      PostData(1, fOutputList);
      return;
   }

   // vertex selection
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
   Float_t centValue = 0.; 
   if(fESD) cent = fESD->GetCentrality();
   if(cent) centValue = cent->GetCentralityPercentile("V0M");
   if(fDebug) printf("centrality: %f\n", centValue);
   if (centValue < fCentMin || centValue > fCentMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fOutputList);
      return;
   }


   // multiplicity due to input tracks
   Int_t nInputTracks = GetNInputTracks();

   if (nInputTracks < fNInputTracksMin || (fNInputTracksMax > -1 && nInputTracks > fNInputTracksMax)){
      fHistEvtSelection->Fill(5);
      PostData(1, fOutputList);
      return;
   }

   
   fHistEvtSelection->Fill(0); // accepted events  
   // -- end event selection --


   // get background
   AliAODJetEventBackground* externalBackground = 0;
   if(!externalBackground&&fBackgroundBranch.Length()){
      externalBackground =  (AliAODJetEventBackground*)(fAOD->FindListObject(fBackgroundBranch.Data()));
      if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   Float_t rho = 0;
   if(externalBackground)rho = externalBackground->GetBackground(0);


   // fetch jets
   TClonesArray *aodJets[2];
   aodJets[0] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[0].Data())); // 
   aodJets[1] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[1].Data())); // 



   for (Int_t iJetType = 0; iJetType < 2; iJetType++) {
      fListJets[iJetType]->Clear();
      if (!aodJets[iJetType]) continue;

      if(fDebug) Printf("%s: %d jets",fJetBranchName[iJetType].Data(),aodJets[iJetType]->GetEntriesFast());
      
      for (Int_t iJet = 0; iJet < aodJets[iJetType]->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets[iJetType])[iJet]);
         if (jet) fListJets[iJetType]->Add(jet);
      }
      fListJets[iJetType]->Sort();
   }
   
   Double_t etabig=0;
   Double_t ptbig=0;
   Double_t areabig=0;
   Double_t phibig=0.;
   Double_t etasmall=0;
   Double_t ptsmall=0;
   Double_t areasmall=0;
   Double_t distr=0.;
   Double_t phismall=0.;
   for(Int_t i=0; i<fListJets[0]->GetEntries(); ++i){
           AliAODJet* jetbig = (AliAODJet*)(fListJets[0]->At(i));
           etabig  = jetbig->Eta();
           phibig  = jetbig->Phi();
           ptbig   = jetbig->Pt();
           if(ptbig==0) continue; 
           areabig = jetbig->EffectiveAreaCharged();
      	   if((etabig<fJetEtaMin)||(etabig>fJetEtaMax)) continue;
                   Double_t dismin=100.;
                   Double_t ptmax=-10.; 
                   Int_t index1=-1;
                   Int_t index2=-1;
                   Double_t fracin=0.; 
	   //compute sum of the pt of the tracks in a concentric cone
           TRefArray *genTrackList = jetbig->GetRefTracks();
           Int_t nTracksGenJet = genTrackList->GetEntriesFast();
           AliAODTrack* genTrack;
             for(Int_t ir=0; ir<nTracksGenJet; ++ir){
             genTrack = (AliAODTrack*)(genTrackList->At(ir));
	     Float_t etr=genTrack->Eta();
             Float_t phir=genTrack->Phi();
             distr=(etr-etabig)*(etr-etabig)+(phir-phibig)*(phir-phibig);
             distr=TMath::Sqrt(distr);
	     if(distr<=fRadioFrac){ fracin=fracin+genTrack->Pt();}}    
	     if(centValue<10) fh2JetCoreMethod3C10->Fill(ptbig-rho*areabig,fracin/ptbig);
             if((centValue>30)&&(centValue<60)) fh2JetCoreMethod3C30->Fill(ptbig-rho*areabig,fracin/ptbig);
             if(centValue>60) fh2JetCoreMethod3C60->Fill(ptbig-rho*areabig,fracin/ptbig); 
             //cout<<"method3  "<<fracin/ptbig<<" "<<ptbig-rho*areabig<<endl;
	     /////////////////////////////////////////////////////////////

                for(Int_t j=0; j<fListJets[1]->GetEntries(); ++j){
                  AliAODJet* jetsmall = (AliAODJet*)(fListJets[1]->At(j));
                  etasmall  = jetsmall->Eta();
                  phismall = jetsmall->Phi();
                  ptsmall   = jetsmall->Pt();
                  areasmall = jetsmall->EffectiveAreaCharged();

                  if((etasmall<fJetEtaMin)||(etasmall>fJetEtaMax)) continue;
                  
                   Float_t tmpDeltaR=(phismall-phibig)*(phismall-phibig)+(etasmall-etabig)*(etasmall-etabig);
		   tmpDeltaR=TMath::Sqrt(tmpDeltaR);

		  if((ptsmall>ptmax)&&(tmpDeltaR<=fRadioFrac)){ptmax=ptsmall;  
		    index2=j;}  
  
                  if(tmpDeltaR<=dismin){ dismin=tmpDeltaR;
		    index1=j;}}

                //method1:most concentric jet=core 
		if(dismin<fMinDist){ AliAODJet* jetmethod1 = (AliAODJet*)(fListJets[1]->At(index1));       
	
		  if(centValue<10) fh2JetCoreMethod1C10->Fill(ptbig-rho*areabig,jetmethod1->Pt()/ptbig);
		  if((centValue>30)&&(centValue<60)) fh2JetCoreMethod1C30->Fill(ptbig-rho*areabig,jetmethod1->Pt()/ptbig);
		  if(centValue>60) fh2JetCoreMethod1C60->Fill(ptbig-rho*areabig,jetmethod1->Pt()/ptbig); 

                  //cout<<"method1  "<<jetmethod1->Pt()/ptbig<<" "<<ptbig<<endl;   
		}
                //method2:hardest contained jet=core   
		if(index2!=-1){ 
                  AliAODJet* jetmethod2 = (AliAODJet*)(fListJets[1]->At(index2));
                  if(centValue<10) fh2JetCoreMethod2C10->Fill(ptbig-rho*areabig,jetmethod2->Pt()/ptbig);
		  if((centValue>30)&&(centValue<60)) fh2JetCoreMethod2C30->Fill(ptbig-rho*areabig,jetmethod2->Pt()/ptbig);
		  if(centValue>60) fh2JetCoreMethod2C60->Fill(ptbig-rho*areabig,jetmethod2->Pt()/ptbig); 

		  //cout<<"method2  "<<jetmethod2->Pt()/ptbig<<" "<<ptbig<<endl;
                               }  

   
      
	
   }
     
      


   PostData(1, fOutputList);
}

void AliAnalysisTaskJetCore::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if (!GetOutputData(1))
   return;
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





