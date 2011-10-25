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

#include "AliAnalysisTaskJetResponseV2.h"

ClassImp(AliAnalysisTaskJetResponseV2)

AliAnalysisTaskJetResponseV2::AliAnalysisTaskJetResponseV2() :
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
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fbJetsMismatch1(kTRUE),
fbJetsMismatch2(kTRUE),
fbJetsRp(kTRUE),
fbJetsDeltaPt(kTRUE),
fbJetsEta(kTRUE),
fbJetsPhi(kTRUE),
fbJetsArea(kTRUE),
fbJetsBeforeCut1(kTRUE),
fbJetsBeforeCut2(kTRUE),
fHistEvtSelection(0x0),
fHistJetSelection(0x0),
fh2JetSelection(0x0),
fhnEvent(0x0),
fhnJetsMismatch1(0x0),
fhnJetsMismatch2(0x0),
fhnJetsRp(0x0),
fhnJetsDeltaPt(0x0),
fhnJetsEta(0x0),
fhnJetsPhi(0x0),
fhnJetsArea(0x0),
fhnJetsBeforeCut1(0x0),
fhnJetsBeforeCut2(0x0)
{
   // default Constructor

   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;
}

AliAnalysisTaskJetResponseV2::AliAnalysisTaskJetResponseV2(const char *name) :
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
fkNbranches(2),
fkEvtClasses(12),
fOutputList(0x0),
fbEvent(kTRUE),
fbJetsMismatch1(kTRUE),
fbJetsMismatch2(kTRUE),
fbJetsRp(kTRUE),
fbJetsDeltaPt(kTRUE),
fbJetsEta(kTRUE),
fbJetsPhi(kTRUE),
fbJetsArea(kTRUE),
fbJetsBeforeCut1(kTRUE),
fbJetsBeforeCut2(kTRUE),
fHistEvtSelection(0x0),
fHistJetSelection(0x0),
fh2JetSelection(0x0),
fhnEvent(0x0),
fhnJetsMismatch1(0x0),
fhnJetsMismatch2(0x0),
fhnJetsRp(0x0),
fhnJetsDeltaPt(0x0),
fhnJetsEta(0x0),
fhnJetsPhi(0x0),
fhnJetsArea(0x0),
fhnJetsBeforeCut1(0x0),
fhnJetsBeforeCut2(0x0)
{
   // Constructor

   fJetBranchName[0] = "";
   fJetBranchName[1] = "";

   fListJets[0] = new TList;
   fListJets[1] = new TList;

   DefineOutput(1, TList::Class());
}

AliAnalysisTaskJetResponseV2::~AliAnalysisTaskJetResponseV2()
{
   delete fListJets[0];
   delete fListJets[1];
}

void AliAnalysisTaskJetResponseV2::SetBranchNames(const TString &branch1, const TString &branch2)
{
   fJetBranchName[0] = branch1;
   fJetBranchName[1] = branch2;
}

void AliAnalysisTaskJetResponseV2::Init()
{

   // check for jet branches
   if(!strlen(fJetBranchName[0].Data()) || !strlen(fJetBranchName[1].Data())){
      AliError("Jet branch name not set.");
   }

}

void AliAnalysisTaskJetResponseV2::UserCreateOutputObjects()
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


   UInt_t entries = 0; // bit coded, see GetDimParams() below
   UInt_t opt = 0;  // bit coded, default (0) or high resolution (1)

   if(fbEvent){
      entries = 1<<0 | 1<<1 | 1<<2 | 1<<26;  // cent : nInpTrks : rp psi : pT hard bin
      opt = 1<<0 | 1<<1; // centrality and nInpTrks in high resolution
      fhnEvent = NewTHnSparseF("fhnEvent", entries, opt);
   }
   
   if(fbJetsMismatch1){ // real mismatch (no related rec jet found)
      // cent : nInpTrks : rp bins : probe pt : probe eta : probe phi : probe area : pT hard bin
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<6 | 1<<8 | 1<<10 | 1<<12 | 1<<26;
      opt =  1<<6 | 1<<8 | 1<<10;
      fhnJetsMismatch1 = NewTHnSparseF("fhnJetsMismatch1", entries, opt);
   }

   if(fbJetsMismatch2){  // acceptance + fraction cut
      // cent : nInpTrks : rp bins : jetPt(2x) : jetEta(2x) : deltaEta : deltaR : fraction : pT hard bin
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<15 | 1<<17 | 1<<19 | 1<<26;
      opt = 1<<6 | 1<<7 | 1<<8 | 1<<9;
      fhnJetsMismatch2 = NewTHnSparseF("fhnJetsMismatch2", entries, opt);

   }
   
   
   if(fbJetsRp){
      // cent : nInpTrks : rp bins : rp wrt jet(2x) : probe pT : probe area: deltaPt : rp phi : rho : correction for RP : local rho
      /*
   entries = 1<<0 | 1<<1 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<12 | 1<<14 | 1<<2 | 1<<22 | 1<<23 | 1<<24 | 1<<25;
   opt = 1<<4 | 1<<5;*/
      // cent : nInpTrks : rp bins : rp wrt jet(2x) : probe pT : deltaPt : rp phi : rho : correction for RP | pT hard bin
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<14 | 1<<2 | 1<<22 | 1<<23 | 1<<26;
      opt = 1<<4 | 1<<5;
      fhnJetsRp = NewTHnSparseF("fhnJetsRp", entries, opt);
   }

   // cent : nInpTrks : rp bins:  deltaPt : jetPt(2x) : deltaArea : pT hard bin (hr delta pt)
   if(fbJetsDeltaPt){
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<14 | 1<<6 | 1<<7 | 1<<18 | 1<<26;
      opt = 1<<1 | 1<<14 | 1<<6 | 1<<7 | 1<<18;
      fhnJetsDeltaPt = NewTHnSparseF("fhnJetsDeltaPt", entries, opt);
   }

   // cent : nInpTrks : rp bins : deltaPt : jetPt(2x) : deltaR : deltaEta : jetEta(2x) : pT hard bin (hr for eta)
   if(fbJetsEta){
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<14 | 1<<6 | 1<<7 | 1<<17 | 1<<15 | 1<<8 | 1<<9 | 1<<26;
      opt = 1<<15 | 1<<8 | 1<<9;
      fhnJetsEta = NewTHnSparseF("fhnJetsEta", entries, opt);
   }

   // cent : nInpTrks : rp bins : jetPt(2x) : jetPhi(2x) : deltaPt : deltaPhi : pT hard bin
   if(fbJetsPhi){
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<6 | 1<<7 | 1<<10 | 1<<11 | 1<<14 | 1<<16 | 1<<26;
      opt = 1<<10 | 1<<11;
      fhnJetsPhi = NewTHnSparseF("fhnJetsPhi", entries, opt);
   }

   // cent : nInpTrks : rp bins : deltaArea : jetArea(2x) : deltaR : fraction : distance next rec jet : pT next jet : deltaPt : jetPt(2x) : pT hard bin (hr for area)
   if(fbJetsArea){
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<18 | 1<<12 | 1<<13 | 1<<17 | 1<<19 | 1<<20 | 1<<21 | 1<<14 | 1<<6 | 1<<7 | 1<<26;
      opt = 1<<18 | 1<<12 | 1<<13;
      fhnJetsArea = NewTHnSparseF("fhnJetsArea", entries, opt);
   }


   //before cut

   // cent : nInpTrks : rp bins : fraction : jetPt(2x) : jetEta(2x) : jetPhi(2x) (low resolution) (with fraction, eta, phi, pt cuts possible)
   if(fbJetsBeforeCut1){
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<19 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11;
      opt = 0;
      fhnJetsBeforeCut1 = NewTHnSparseF("fhnJetsBeforeCut1", entries, opt);
   }

   // cent : nInpTrks : rp bins : deltaPt : jetPt(2x) : deltaR : deltaEta : jetEta(2x) (low resolution)
   if(fbJetsBeforeCut2){
      entries = 1<<0 | 1<<1 | 1<<3 | 1<<14 | 1<<6 | 1<<7 | 1<<17 | 1<<15 | 1<<8 | 1<<9;
      opt = 0;
      fhnJetsBeforeCut2 = NewTHnSparseF("fhnJetsBeforeCut2", entries, opt);
   }

   fOutputList->Add(fHistEvtSelection);
   fOutputList->Add(fHistJetSelection);
   fOutputList->Add(fh2JetSelection);
   if(fbEvent)           fOutputList->Add(fhnEvent);
   if(fbJetsMismatch1)   fOutputList->Add(fhnJetsMismatch1);
   if(fbJetsMismatch2)   fOutputList->Add(fhnJetsMismatch2);
   if(fbJetsRp)          fOutputList->Add(fhnJetsRp);
   if(fbJetsDeltaPt)     fOutputList->Add(fhnJetsDeltaPt);
   if(fbJetsEta)         fOutputList->Add(fhnJetsEta);
   if(fbJetsPhi)         fOutputList->Add(fhnJetsPhi);
   if(fbJetsArea)        fOutputList->Add(fhnJetsArea);
   if(fbJetsBeforeCut1)  fOutputList->Add(fhnJetsBeforeCut1);
   if(fbJetsBeforeCut2)  fOutputList->Add(fhnJetsBeforeCut2);

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

void AliAnalysisTaskJetResponseV2::UserExec(Option_t *)
{
   // load events, apply event cuts, then compare leading jets

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

   // pt hard
   Double_t pthard = AliAnalysisTaskFastEmbedding::GetPtHard();
   Int_t pthardbin = GetPtHardBin(pthard);

   // reaction plane
   Double_t rp = AliAnalysisHelperJetTasks::ReactionPlane(kFALSE);
   if(fbEvent){
      Double_t eventEntries[3] = { (Double_t)centValue, (Double_t)nInputTracks, rp };
      fhnEvent->Fill(eventEntries);
   }


   // get background
   AliAODJetEventBackground* externalBackground = 0;
   if(!externalBackground&&fBackgroundBranch.Length()){
      externalBackground =  (AliAODJetEventBackground*)(fAOD->FindListObject(fBackgroundBranch.Data()));
      //if(!externalBackground)Printf("%s:%d Background branch not found %s",(char*)__FILE__,__LINE__,fBackgroundBranch.Data());;
   }
   Float_t rho = 0;
   if(externalBackground)rho = externalBackground->GetBackground(0);


   // fetch jets
   TClonesArray *aodJets[2];
   aodJets[0] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[0].Data())); // in general: embedded jet
   aodJets[1] = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fJetBranchName[1].Data())); // in general: embedded jet + UE



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
   
   
   // jet matching 
   static TArrayI aMatchIndex(fListJets[0]->GetEntries());
   static TArrayF aPtFraction(fListJets[0]->GetEntries());
   if(aMatchIndex.GetSize()<fListJets[0]->GetEntries()) aMatchIndex.Set(fListJets[0]->GetEntries());
   if(aPtFraction.GetSize()<fListJets[0]->GetEntries()) aPtFraction.Set(fListJets[0]->GetEntries());

   // stores matched jets in 'aMatchIndex' and fraction of pT in 'aPtFraction'
   AliAnalysisHelperJetTasks::GetJetMatching(fListJets[0], TMath::Min((Int_t)fNMatchJets,(Int_t)fListJets[0]->GetEntries()),
   fListJets[1], TMath::Min((Int_t)fNMatchJets,(Int_t)fListJets[1]->GetEntries()),
   aMatchIndex, aPtFraction, fDebug, fMatchMaxDist, fIsPbPb?1:2);
   
   // loop over matched jets
   Int_t ir = -1; // index of matched reconstruced jet
   Float_t fraction = -1.;
   AliAODJet *jet[2]  = { 0x0, 0x0 };
   Float_t jetEta[2]  = { -990., -990. };
   Float_t jetPhi[2]  = { -990., -990. };
   Float_t jetPt[2]   = { -990., -990. };
   Float_t jetArea[2] = { -990., -990. };
   Float_t rpJet[2]   = { -990., -990. };
   Int_t rpBin = -990;
   
   for(Int_t ig=0; ig<fListJets[0]->GetEntries(); ++ig){
      ir = aMatchIndex[ig];

      //fetch jets
      jet[0] = (AliAODJet*)(fListJets[0]->At(ig));
      if(ir>=0) jet[1] = (AliAODJet*)(fListJets[1]->At(ir));
      else      jet[1] = 0x0;

      for(Int_t i=0; i<fkNbranches; ++i){
         if(!jet[i]){
            jetEta[i]  = -990;
            jetPhi[i]  = -990.;
            jetPt[i]   = -990.;
            jetArea[i] = -990.;
            rpJet[i]   = -990.;
            if(i==1) rpBin = -990;
         } 
         else {
            jetEta[i]  = jet[i]->Eta();
            jetPhi[i]  = jet[i]->Phi();
            jetPt[i]   = GetPt(jet[i], i);
            jetArea[i] = jet[i]->EffectiveAreaCharged();
            rpJet[i]   = TVector2::Phi_mpi_pi(rp-jetPhi[i]);
            if(i==1) rpBin = AliAnalysisHelperJetTasks::GetPhiBin(TVector2::Phi_mpi_pi(rp-jetPhi[i]), 3);
         }
      }
      fraction = aPtFraction[ig];

      // jet statistics
      fHistJetSelection->Fill(1); // all probe jets
      if(jet[0]) fh2JetSelection->Fill(1.,jetPt[0]); // all probe jets

      if(ir<0){
         fHistJetSelection->Fill(2);
         if(jet[0]) fh2JetSelection->Fill(2.,jetPt[0]);
         
         if(fbJetsMismatch1){
            if(!jet[0]) continue;
            Double_t jetEntriesMismatch1[7] = {
               (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
               (Double_t)jetPt[0], (Double_t)jetEta[0], (Double_t)jetPhi[0], (Double_t)jetArea[0]
            };		 
            fhnJetsMismatch1->Fill(jetEntriesMismatch1);
         }
         
         continue;
      }
      
      if(!jet[0] || !jet[1]){
         fHistJetSelection->Fill(3);
         if(jet[0]) fh2JetSelection->Fill(3.,jetPt[0]);
         continue;
      }
      
      // look for distance to next rec jet
      Float_t distNextJet = -0.01; // no neighbor
      Float_t ptNextJet = -1.;
      for(Int_t r=0; r<fListJets[1]->GetEntries(); ++r){
         if(r==ir) continue;
         Float_t tmpDeltaR = jet[1]->DeltaR((AliAODJet*)fListJets[1]->At(r));
         if(distNextJet<0. || distNextJet>tmpDeltaR){
            distNextJet = tmpDeltaR;
            if(fKeepJets) ptNextJet   = ((AliAODJet*)fListJets[1]->At(r))->GetPtSubtracted(0);
            else          ptNextJet   = ((AliAODJet*)fListJets[1]->At(r))->Pt();
         }
      }
      


      //Float_t localRho = jetArea[1]>0. ? (jetPt[1]+rho*jetArea[1] - jetPt[0]) / jetArea[1] : 0.;
      //Float_t relRho = rho>0. ? localRho / rho : 0.;

      
      // calculate parameters of associated jets
      /* from Leticia, kT clusterizer
   Float_t par0[4] = { 0.00409,       0.01229,      0.05,         0.26 };
   Float_t par1[4] = { -2.97035e-03, -2.03182e-03, -1.25702e-03, -9.95107e-04 };
   Float_t par2[4] = { 1.02865e-01,   1.49039e-01,  1.53910e-01,  1.51109e-01 };
   */
      // own, from embedded tracks
      Float_t par0[4] = { 0.02841,       0.05039,      0.09092,      0.24089     };
      Float_t par1[4] = { -4.26725e-04, -1.15273e-03, -1.56827e-03, -3.08003e-03 };
      Float_t par2[4] = { 4.95415e-02,   9.79538e-02,  1.32814e-01,  1.71743e-01 };
      
      Float_t rpCorr = 0.;         
      
      if(eventClass>0&&eventClass<4){
         rpCorr = par0[eventClass-1] + par1[eventClass-1] * 2*TMath::Cos(rpJet[1]) + par2[eventClass-1] * 2*TMath::Cos(2*rpJet[1]);
         rpCorr *= rho * jetArea[1];
      }

      Float_t deltaPt    = jetPt[1]-jetPt[0];
      Float_t deltaEta   = jetEta[1]-jetEta[0];
      Float_t deltaPhi   = TVector2::Phi_mpi_pi(jetPhi[1]-jetPhi[0]);
      Float_t deltaR     = TMath::Sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
      Float_t deltaArea  = jetArea[1]-jetArea[0];
      
      
      // fill thnsparse before acceptance cut
      if(fbJetsBeforeCut1){
         Double_t jetBeforeCutEntries1[10] = { 
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
            (Double_t)jetPt[0], (Double_t)jetPt[1], (Double_t)jetEta[0], (Double_t)jetEta[1], 
            (Double_t)jetPhi[0], (Double_t)jetPhi[1], (Double_t)fraction
         };
         fhnJetsBeforeCut1->Fill(jetBeforeCutEntries1);
      }
      
      if(fbJetsBeforeCut2){
         Double_t jetBeforeCutEntries2[10] = {
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin, 
            (Double_t)jetPt[0], (Double_t)jetPt[1],
            (Double_t)jetEta[0], (Double_t)jetEta[1], 
            (Double_t)deltaPt, (Double_t)deltaEta, (Double_t)deltaR
         };
         fhnJetsBeforeCut2->Fill(jetBeforeCutEntries2);
      }
      
      // jet selection
      Bool_t jetAccepted = kTRUE;
      // minimum fraction required
      if(fraction<fJetPtFractionMin){
         fHistJetSelection->Fill(4);
         fh2JetSelection->Fill(4.,jetPt[0]);
         jetAccepted = kFALSE;
      }
      
      if(jetAccepted){
         // jet acceptance + minimum pT check
         if(jetEta[0]>fJetEtaMax || jetEta[0]<fJetEtaMin ||
               jetEta[1]>fJetEtaMax || jetEta[1]<fJetEtaMin){
            
            if(fDebug){
               Printf("Jet not in eta acceptance.");
               Printf("[0]: jet %d eta %.2f", ig, jetEta[0]);
               Printf("[1]: jet %d eta %.2f", ir, jetEta[1]);
            }
            fHistJetSelection->Fill(5);
            fh2JetSelection->Fill(5.,jetPt[0]);
            jetAccepted = kFALSE;
         }
      }
      if(jetAccepted){
         if(jetPt[1] < fJetPtMin){
            if(fDebug) Printf("Jet %d (pT %.1f GeV/c) has less than required pT.", ir, jetPt[1]);
            fHistJetSelection->Fill(6);
            fh2JetSelection->Fill(6.,jetPt[0]);
            jetAccepted = kFALSE;
         }
      }

      if(jetAccepted){
         if(jet[1]->Trigger()&fJetTriggerExcludeMask){
            fHistJetSelection->Fill(7);
            fh2JetSelection->Fill(7.,jetPt[0]);
            jetAccepted = kFALSE;
         }
         
      }
      
      if(!jetAccepted){
         if(fbJetsMismatch2){
            Double_t jetEntriesMismatch2[10] = {
               (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
               (Double_t)jetPt[0], (Double_t)jetPt[1],
               (Double_t)jetEta[0], (Double_t)jetEta[1],
               (Double_t)deltaEta, (Double_t)deltaR,
               (Double_t)fraction
            };
            fhnJetsMismatch2->Fill(jetEntriesMismatch2);     
         }
         continue;
      }
      
      // all accepted jets
      fHistJetSelection->Fill(0);
      fh2JetSelection->Fill(0.,jetPt[0]);
      
      // fill thnsparse
      if(fbJetsRp){
         Double_t jetEntriesRp[11] = {
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t) rp,
            (Double_t)rpBin, (Double_t)rpJet[0], (Double_t)rpJet[1], 
            (Double_t)jetPt[0], (Double_t)deltaPt, (Double_t)rho, (Double_t)rpCorr,
            (Double_t)pthardbin
         };
         fhnJetsRp->Fill(jetEntriesRp);
      }
      
      if(fbJetsDeltaPt){
         Double_t jetEntriesDeltaPt[8] = {
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
            (Double_t)jetPt[0], (Double_t)jetPt[1], (Double_t)deltaPt, (Double_t)deltaArea,
            (Double_t)pthardbin
         };		 
         fhnJetsDeltaPt->Fill(jetEntriesDeltaPt);
      }
      
      if(fbJetsEta){
         Double_t jetEntriesEta[11] = {
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
            (Double_t)jetPt[0], (Double_t)jetPt[1], (Double_t)jetEta[0], (Double_t)jetEta[1], 
            (Double_t)deltaPt, (Double_t)deltaEta, (Double_t)deltaR,
            (Double_t)pthardbin
         };				 
         fhnJetsEta->Fill(jetEntriesEta);
      }
      
      if(fbJetsPhi){
         Double_t jetEntriesPhi[10] = {
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
            (Double_t)jetPt[0], (Double_t)jetPt[1], (Double_t)jetPhi[0], (Double_t)jetPhi[1],
            (Double_t)deltaPt, (Double_t)deltaPhi,
            (Double_t)pthardbin
         };
         fhnJetsPhi->Fill(jetEntriesPhi);
      }
      
      if(fbJetsArea){
         Double_t jetEntriesArea[14] = {
            (Double_t)centValue, (Double_t)nInputTracks, (Double_t)rpBin,
            (Double_t)jetPt[0], (Double_t)jetPt[1], (Double_t)jetArea[0], (Double_t)jetArea[1], 
            (Double_t)deltaPt, (Double_t)deltaR, (Double_t)deltaArea, 
            (Double_t)fraction, (Double_t)distNextJet, (Double_t)ptNextJet,
            (Double_t)pthardbin
         };				 
         fhnJetsArea->Fill(jetEntriesArea);
      }
      
   }

   PostData(1, fOutputList);
}

void AliAnalysisTaskJetResponseV2::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if (!GetOutputData(1))
   return;
}

Int_t AliAnalysisTaskJetResponseV2::GetNInputTracks()
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
      Printf("AliAnalysisTaskJetResponseV2::GetNInputTracks FAILED");
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

THnSparse* AliAnalysisTaskJetResponseV2::NewTHnSparseF(const char* name, UInt_t entries, UInt_t opt)
{
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
         Bool_t highres = opt&(1<<i);
         TString label("");
         GetDimParams(i, highres, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }
      
      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}

void AliAnalysisTaskJetResponseV2::GetDimParams(Int_t iEntry, Bool_t hr, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{

   const Double_t pi = TMath::Pi();
   
   switch(iEntry){
      
   case 0:
      label = "V0 centrality (%)";
      if(hr){
         nbins = 100;
         xmin = 0.;
         xmax = 100.;
      } else {
         nbins = 8;
         xmin = 0.;
         xmax = 80.;
      }
      break;
      
      
   case 1:
      label = "nb. of input tracks";
      if(fIsPbPb){
         if(hr){
            nbins = 400;
            xmin = 0.;
            xmax = 4000.;
         } else {
            nbins = 40;
            xmin = 0.;
            xmax = 4000.;
         }
      } else {
         nbins = 40;
         xmin = 0.;
         xmax = 400.;
      }
      break;
      
      
   case 2:
      label = "event plane #psi";
      if(hr){
         nbins = 30;
         xmin = 0.;
         xmax = pi;
      } else {
         nbins = 30;
         xmin = 0.;
         xmax = pi;
      }
      break;
      
      
   case 3:
      label = "event plane bin";
      nbins = 3;
      xmin = -.5;
      xmax = 2.5;
      break;
      
      
   case 4:
   case 5:
      if(iEntry==4)label = "#Delta#phi(RP-jet) (probe)";
      if(iEntry==5)label = "#Delta#phi(RP-jet) (rec)";
      nbins = 48;
      xmin = -pi;
      xmax =  pi;
      break;
      
      
   case 6:
   case 7:
      if(iEntry==6)label = "probe p_{T} (GeV/c)";
      if(iEntry==7)label = "rec p_{T} (GeV/c)";
      if(hr){
         nbins = 300;
         xmin = -50.;
         xmax = 250.;
      } else {
         nbins = 50;
         xmin = 0.;
         xmax = 250.;
      }
      break;
      
      
   case 8:
   case 9:
      if(iEntry==8)label = "probe #eta";
      if(iEntry==9)label = "rec #eta";
      if(hr){
         nbins = 56;
         xmin = -.7;
         xmax =  .7;
      } else {
         nbins = 28;
         xmin = -.7;
         xmax =  .7;
      }
      break;
      
      
   case 10:
   case 11:
      if(iEntry==10)label = "probe #phi";
      if(iEntry==11)label = "rec #phi";
      if(hr){
         nbins = 90;  // modulo 18 (sectors)
         xmin = 0.;
         xmax = 2*pi;
      } else {
         nbins = 25;
         xmin = 0.;
         xmax = 2*pi;
      }
      break;
      
      
   case 12:
   case 13:
      if(iEntry==12)label = "probe area";
      if(iEntry==13)label = "rec area";
      if(hr){
         nbins = 100;
         xmin = 0.;
         xmax = 1.;
      } else {
         nbins = 25;
         xmin = 0.;
         xmax = 1.;
      }
      break;
      
   case 14:
      label = "#Delta p_{T}";
      if(hr){
         nbins = 241;
         xmin = -120.5;
         xmax =  120.5;
      } else {
         nbins = 101;
         xmin = -101.;
         xmax =  101.;
      }
      break;
      
   case 15:
      label = "#Delta#eta";
      if(hr){
         nbins = 51;
         xmin = -1.02;
         xmax =  1.02;
      } else {
         nbins = 51;
         xmin = -1.02;
         xmax =  1.02;
      }
      break;
      
      
   case 16:
      label = "#Delta#phi";
      if(hr){
         nbins = 45;
         xmin = -pi;
         xmax =  pi;
      } else {
         nbins = 45;
         xmin = -pi;
         xmax =  pi;
      }
      break;
      
      
   case 17:
      label = "#DeltaR";
      if(hr){
         nbins = 50;
         xmin = 0.;
         xmax = 1.;
      } else {
         nbins = 50;
         xmin = 0.;
         xmax = 1.;
      }
      break;
      
      
   case 18:
      label = "#Deltaarea";
      if(hr){
         nbins = 81;
         xmin = -.81;
         xmax =  .81;
      } else {
         nbins = 33;
         xmin = -.825;
         xmax =  .825;
      }
      break;
      
      
   case 19:
      label = "fraction";
      if(hr){
         nbins = 52;
         xmin = 0.;
         xmax = 1.04;
      } else {
         nbins = 52;
         xmin = 0.;
         xmax = 1.04;
      }
      break;
      
      
   case 20:
      label = "distance to closest rec jet";
      if(hr){
         nbins = 51;
         xmin = -0.02;
         xmax =  1.;
      } else {
         nbins = 51;
         xmin = -0.02;
         xmax = 1.;
      }
      break;
      
      
   case 21:
      label = "p_{T} of closest rec jet";
      nbins = 100;
      xmin =  0.;
      xmax =  200.;
      break;
      
   case 22:
      label = "#rho";
      nbins = 125;
      xmin  = 0.;
      xmax  = 250.;
      break;
      
   case 23:
      label = "abs. correction of #rho for RP";
      nbins =  51;
      xmin  = -51.;
      xmax  =  51.;
      break;
      
   case 24:
      label = "local #rho";
      nbins = 125;
      xmin  = 0.;
      xmax  = 250.;
      break;
      
   case 25:
      label = "local #rho / #rho";
      nbins = 500;
      xmin  = 0.;
      xmax  = 5.;
      break;
      
   case 26:
      label = "p_{T,hard} bin";
      nbins = 10;
      xmin  = -.5;
      xmax  = 9.5;
      break;
      
   }

}

//____________________________________________________________________________
Int_t AliAnalysisTaskJetResponseV2::GetPtHardBin(Double_t ptHard){

   const Int_t nBins = 10;
   Double_t binLimits[nBins] = { 5., 11., 21., 36., 57., 84., 117., 156., 200., 249. }; // lower limits
   
   Int_t bin = -1;
   while(bin<nBins-1 && binLimits[bin+1]<ptHard){
      bin++;
   }
   
   return bin;
}

//____________________________________________________________________________
Double_t AliAnalysisTaskJetResponseV2::GetPt(AliAODJet *j, Int_t mode=0){

   Double_t pt = 0.;

   if(fKeepJets && mode==1){  // background subtracted pt, can be negative
      pt = j->GetPtSubtracted(0);
   }
   else{
      pt = j->Pt();  // if negative, pt=0.01
   }

   return pt;
}
