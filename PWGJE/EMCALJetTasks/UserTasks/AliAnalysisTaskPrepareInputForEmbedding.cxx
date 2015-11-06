#include <TTree.h>
#include <TLorentzVector.h>

#include <AliLog.h>
#include <AliEmcalJet.h>
#include "AliJetContainer.h"
#include "AliRhoParameter.h"
#include "AliAnalysisTaskPrepareInputForEmbedding.h"


ClassImp(AliAnalysisTaskPrepareInputForEmbedding)

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding() : AliAnalysisTaskEmcalJet("AliAnalysisTaskPrepareInputForEmbedding", kTRUE),
fContainer(0),
fMinFractionShared(0.5),
fLeadingJetOnly(0),
fTreeJets(0),
fJetDet(0),
fJetPart(0),
fJetDetL(0),
fJetPartL(0),
fNumberOfJets(0),
fhFractionSharedpT(0)
{
   /// default constructor
   //fJetGenSub = new TLorentzVector();
   //fJetConSub = new TLorentzVector();
   //fJetPart1   = new TLorentzVector();
   //fJetPart2   = new TLorentzVector();
   SetMakeGeneralHistograms(kTRUE);
   
   //DefineOutput(2, TTree::Class());
}

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
fContainer(0),
fMinFractionShared(0.5),
fLeadingJetOnly(0),
fTreeJets(0),
fJetDet(0),
fJetPart(0),
fJetDetL(0),
fJetPartL(0),
fNumberOfJets(0),
fhFractionSharedpT(0)

{
   /// standard constructor
   
   SetMakeGeneralHistograms(kTRUE);
   AliInfo("standard constructor");
}
//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::UserCreateOutputObjects(){
   /// Create output
   AliAnalysisTaskEmcal::UserCreateOutputObjects();
   AliInfo("Sei qui?");
   fJetDet = new TLorentzVector();
   fJetPart  = new TLorentzVector();

   fJetDetL = new TLorentzVector();
   fJetPartL  = new TLorentzVector();

   fTreeJets = new TTree("fTreeJet", "fTreeJet");
   if(fLeadingJetOnly){
      fTreeJets->Branch("fJetDetL", &fJetDetL);
      fTreeJets->Branch("fJetPartL",&fJetPartL);
   } else {
      fTreeJets->Branch("fJetDet", &fJetDet);
      fTreeJets->Branch("fJetPart",&fJetPart);
   }
   fOutput->Add(fTreeJets);

   //check
   const Int_t nBinsPt          = 40;
   const Int_t nBinsFraction    = 101;
   
   const Double_t minPt       = -50.;
   const Double_t maxPt       = 150.;
   const Double_t minFraction =  -0.005;
   const Double_t maxFraction =  1.005;
   
   fNumberOfJets = new TH1F("fNumberOfJets", "Number of Jets", 6, -0.5, 5.5);
   fOutput->Add(fNumberOfJets);
   fhFractionSharedpT = new TH2F("fhFractionSharedpT", "Reco/particle shared #it{p}_{T} fraction;#it{p}_{T} (GeV/c); Fraction shared #it{p}_{T}", nBinsPt, minPt, maxPt, nBinsFraction, minFraction, maxFraction);
   fOutput->Add(fhFractionSharedpT);
   
}


//________________________________________________________________________________________________
Bool_t AliAnalysisTaskPrepareInputForEmbedding::Run(){
   /// Run code before FillHistograms()
      
   
   return kTRUE;
}
//________________________________________________________________________________________________
Bool_t AliAnalysisTaskPrepareInputForEmbedding::FillHistograms(){
   /// Fill histograms defined
   /// In this case it's a TTree
   
   // Get container with area-based subtraction
   AliEmcalJet* jet = NULL;
   AliJetContainer *jetCont = GetJetContainer(fContainer);
   if(!jetCont) {
      AliError(Form("Container position %d (area based) not found", fContainer));
      return kFALSE;
   }
   
   //TVector for leading jet
   fJetDetL ->SetPtEtaPhiM(0,0,0,0);
   fJetPartL->SetPtEtaPhiM(0,0,0,0);
   
   if(!fLeadingJetOnly) {
      fTreeJets->SetBranchAddress("fJetDet", &fJetDet);
      fTreeJets->SetBranchAddress("fJetPart",&fJetPart);
   } else {
      fTreeJets->SetBranchAddress("fJetDetL", &fJetDetL);
      fTreeJets->SetBranchAddress("fJetPartL",&fJetPartL);
   }
      
   while((jet = jetCont->GetNextAcceptJet())) {
      fJetDet  ->SetPtEtaPhiM(0,0,0,0);
      fJetPart ->SetPtEtaPhiM(0,0,0,0);
      
      fNumberOfJets->Fill(0);
      
      AliEmcalJet *jetP = jet->ClosestJet();
      if(!jetP) continue;
      fNumberOfJets->Fill(1);
      Double_t fraction = jetCont->GetFractionSharedPt(jet);
      //fill the TLorentsVectors with the jet 4-vectors
      //Printf("MC jet %p, Fraction %.4f, pT %f, eta %f, phi %f",jetP, fraction, jetA->Pt(), jetA->Eta(), jetA->Phi());
      fhFractionSharedpT->Fill(jet->Pt(), fraction);
      //if(fMinFractionShared>0. && fraction>fMinFractionShared) {
      	 fNumberOfJets->Fill(2);
      	 fJetDet ->SetPtEtaPhiM(jet->Pt(), jet->Eta(), jet->Phi(), jet->M());
      	 fJetPart->SetPtEtaPhiM(jetP->Pt(), jetP->Eta(), jetP->Phi(), jetP->M());
      	
      	 if(fJetDet->Pt() > fJetDetL->Pt()){
      	    fJetDetL->SetPtEtaPhiM(fJetDet->Pt(), fJetDet->Eta(), fJetDet->Phi(), fJetDet->M());
      	    fJetPartL->SetPtEtaPhiM(fJetPart->Pt(), fJetPart->Eta(), fJetPart->Phi(), fJetPart->M());
      	 
      	 }
      	 if(!fLeadingJetOnly) fTreeJets->Fill();
      //}
      
      
   }
   
   if(fLeadingJetOnly) {
      if(fJetPartL->Pt() > 0){ // fill only when there is a leading jet
      	 fTreeJets->Fill();
      }
   }

   return kTRUE;
}


//________________________________________________________________________________________________

AliAnalysisTaskPrepareInputForEmbedding::~AliAnalysisTaskPrepareInputForEmbedding() {
   /// Destructor.
   
   delete fJetDet  ;
   delete fJetPart ;

   delete fJetDetL ;
   delete fJetPartL;
}

//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::Terminate(Option_t *) 
{
   /// Called once at the end of the analysis.
}
