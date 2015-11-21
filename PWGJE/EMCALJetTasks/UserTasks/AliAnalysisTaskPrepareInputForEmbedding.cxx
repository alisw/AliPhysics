#include <TTree.h>
#include <TLorentzVector.h>

#include <AliLog.h>
#include <AliEmcalJet.h>
#include "AliJetContainer.h"
#include "AliAnalysisTaskPrepareInputForEmbedding.h"


ClassImp(AliAnalysisTaskPrepareInputForEmbedding)

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding() : AliAnalysisTaskEmcalJet("AliAnalysisTaskPrepareInputForEmbedding", kTRUE),
fContainer(0),
fMinFractionShared(-1),
fLeadingJetOnly(0),
fTreeJets(0),
fJetDet(0),
fJetPart(0),
fJetDetL(0),
fJetPartL(0),
fNumberOfJets(0),
fhFractionSharedpT(0),
fNAccJets(0)
{
   /// default constructor
   
   SetMakeGeneralHistograms(kTRUE);
   
   //DefineOutput(2, TTree::Class());
}

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
fContainer(0),
fMinFractionShared(-1),
fLeadingJetOnly(0),
fTreeJets(0),
fJetDet(0),
fJetPart(0),
fJetDetL(0),
fJetPartL(0),
fNumberOfJets(0),
fhFractionSharedpT(0),
fNAccJets(0)

{
   /// standard constructor
   
   SetMakeGeneralHistograms(kTRUE);
   AliInfo("standard constructor");
}
//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::UserCreateOutputObjects(){
   /// Create output
   AliAnalysisTaskEmcal::UserCreateOutputObjects();

   fJetDet   = new TLorentzVector();
   fJetPart  = new TLorentzVector();
   fJetDetL  = new TLorentzVector();
   fJetPartL = new TLorentzVector();
   
   fTreeJets = new TTree("fTreeJet", "fTreeJet");
   //Important! 'dot' after the name needed -- see TTree doc! otherwise the objects in the branches have all the same name and the output cannot be properly retrieved (only with the TBrowser)
   if(fLeadingJetOnly){
      fTreeJets->Branch("fJetDetL.", fJetDetL);
      fTreeJets->Branch("fJetPartL.",fJetPartL);
   } else {
      fTreeJets->Branch("fJetDet.", fJetDet);
      fTreeJets->Branch("fJetPart.",fJetPart);
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
   
   fNAccJets = new TH1F("fNAccJets","fNAccJets;N/ev",11,-0.5, 9.5);
   fOutput->Add(fNAccJets);
   
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
      
   Int_t count = 0;
   jetCont->ResetCurrentID();
   
   while((jet = jetCont->GetNextJet())) {
      Bool_t acc = jetCont->AcceptJet(jet);
      if(!acc) {
      	 continue;
      }
      fJetDet  ->SetPtEtaPhiM(0,0,0,0);
      fJetPart ->SetPtEtaPhiM(0,0,0,0);
      
      count++;
      fNumberOfJets->Fill(0);
      
      AliEmcalJet *jetP = jet->ClosestJet();
      if(!jetP) continue;
      fNumberOfJets->Fill(1);
      Double_t fraction = jetCont->GetFractionSharedPt(jet);
      //fill the TLorentsVectors with the jet 4-vectors
      //Printf("MC jet %p, Fraction %.4f, pT %f, eta %f, phi %f, m  %f",jetP, fraction, jet->Pt(), jet->Eta(), jet->Phi(), jet->M());
      fhFractionSharedpT->Fill(jet->Pt(), fraction);
      if(fMinFractionShared<0. || fraction>fMinFractionShared) {
      	 fNumberOfJets->Fill(2);
      	 fJetDet ->SetPtEtaPhiM(jet->Pt(), jet->Eta(), jet->Phi(), jet->M());
      	 fJetPart->SetPtEtaPhiM(jetP->Pt(), jetP->Eta(), jetP->Phi(), jetP->M());
      	
      	 if(fJetDet->Pt() > fJetDetL->Pt()){
      	    fJetDetL->SetPtEtaPhiM(fJetDet->Pt(), fJetDet->Eta(), fJetDet->Phi(), fJetDet->M());
      	    fJetPartL->SetPtEtaPhiM(fJetPart->Pt(), fJetPart->Eta(), fJetPart->Phi(), fJetPart->M());
      	 
      	 }
      	 if(!fLeadingJetOnly) fTreeJets->Fill();
      	 
      }
      
   }
   fNAccJets->Fill(count);
   
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
   
 }

//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::Terminate(Option_t *) 
{
   /// Called once at the end of the analysis.

}
