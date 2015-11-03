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
fContainerArea(0),
fContainerConst(1),
fMinFractionShared(0.5),
fLeadingJetOnly(0),
fTreeJetsA(0),
fTreeJetsC(0),
fJetGenSub(0),
fJetConSub(0),
fJetPart1(0),
fJetPart2(0),
fJetGenSubL(0),
fJetConSubL(0),
fJetPart1L(0),
fJetPart2L(0),
fNumberOfJets(0)
{
   /// default constructor
   //fJetGenSub = new TLorentzVector();
   //fJetConSub = new TLorentzVector();
   //fJetPart1   = new TLorentzVector();
   //fJetPart2   = new TLorentzVector();
   SetMakeGeneralHistograms(kTRUE);
   DefineOutput(2, TTree::Class());
}

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
fContainerArea(0),
fContainerConst(1),
fMinFractionShared(0.5),
fLeadingJetOnly(0),
fTreeJetsA(0),
fTreeJetsC(0),
fJetGenSub(0),
fJetConSub(0),
fJetPart1(0),
fJetPart2(0),
fJetGenSubL(0),
fJetConSubL(0),
fJetPart1L(0),
fJetPart2L(0),

fNumberOfJets(0)

{
   /// standard constructor
   
   SetMakeGeneralHistograms(kTRUE);
   
}
//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::UserCreateOutputObjects(){
   /// Create output
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
   fJetGenSub = new TLorentzVector();
   fJetConSub = new TLorentzVector();
   fJetPart1  = new TLorentzVector();
   fJetPart2  = new TLorentzVector();

   fJetGenSubL = new TLorentzVector();
   fJetConSubL = new TLorentzVector();
   fJetPart1L  = new TLorentzVector();
   fJetPart2L  = new TLorentzVector();

   fTreeJetsA = new TTree("fTreeJetA", "fTreeJetA");
   if(fLeadingJetOnly){
      fTreeJetsA->Branch("fJetGenSubL",&fJetGenSubL);
      fTreeJetsA->Branch("fJetPart1L"  ,&fJetPart1L);
   } else {
      fTreeJetsA->Branch("fJetGenSub",&fJetGenSub);
      fTreeJetsA->Branch("fJetPart1"  ,&fJetPart1);
   }
   fOutput->Add(fTreeJetsA);
   
   fTreeJetsC = new TTree("fTreeJetC", "fTreeJetC");
   if(fLeadingJetOnly){
      fTreeJetsC->Branch("fJetConSubL",&fJetConSubL);
      fTreeJetsC->Branch("fJetPart2L"  ,&fJetPart2L);
   } else {
      fTreeJetsC->Branch("fJetConSub", &fJetConSub);
      fTreeJetsC->Branch("fJetPart2"  ,&fJetPart2);
   }
   fOutput->Add(fTreeJetsC);
   
   //check
   const Int_t nBinsPt          = 40;
   const Int_t nBinsFraction    = 101;
   
   const Double_t minPt       = -50.;
   const Double_t maxPt       = 150.;
   const Double_t minFraction =  -0.005;
   const Double_t maxFraction =  1.005;
   
   fNumberOfJets = new TH1F("fNumberOfJets", "Number of Jets", 6, -0.5, 5.5);
   fOutput->Add(fNumberOfJets);
   fhFractionSharedpTA = new TH2F("fhFractionSharedpTA", "Area Based;#it{p}_{T} (GeV/c); Fraction shared #it{p}_{T}", nBinsPt, minPt, maxPt, nBinsFraction, minFraction, maxFraction);
   fOutput->Add(fhFractionSharedpTA);
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
   AliEmcalJet* jetA = NULL;
   AliJetContainer *jetContA = GetJetContainer(fContainerArea);
   if(!jetContA) {
      AliError(Form("Container position %d (area based) not found", fContainerArea));
      return kFALSE;
   } else {
      AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(jetContA->GetRhoMassName()));
      Double_t rhomassval = 0;
      if (!rhomParam) {
      	 AliError(Form("%s: Could not retrieve rho_m %s!", GetName(), jetContA->GetRhoMassName().Data()));
      	 
      } else rhomassval = rhomParam->GetVal();
   }
   
   //TVector for leading jet
   fJetGenSubL->SetPtEtaPhiM(0,0,0,0);
   fJetConSubL->SetPtEtaPhiM(0,0,0,0);
   fJetPart1L ->SetPtEtaPhiM(0,0,0,0);
   fJetPart2L ->SetPtEtaPhiM(0,0,0,0);
   
   if(!fLeadingJetOnly) {
      fTreeJetsA->SetBranchAddress("fJetGenSub",&fJetGenSub);
      fTreeJetsA->SetBranchAddress("fJetPart1"  ,&fJetPart1);
   } else {
      fTreeJetsA->SetBranchAddress("fJetGenSubL",&fJetGenSubL);
      fTreeJetsA->SetBranchAddress("fJetPart1L"  ,&fJetPart1L);
   }
   // Get container with constituent-based subtraction
   AliEmcalJet* jetC = NULL;
   AliJetContainer *jetContC = GetJetContainer(fContainerConst);
   if(!jetContC) {
      AliError(Form("Container position %d (constituent based) not found", fContainerConst));
      return kFALSE;
   }
   
      
   while((jetA = jetContA->GetNextAcceptJet())) {
      fJetGenSub->SetPtEtaPhiM(0,0,0,0);
      fJetPart1  ->SetPtEtaPhiM(0,0,0,0);
      
      fNumberOfJets->Fill(0);
      
      AliEmcalJet *jetP = jetA->ClosestJet();
      if(!jetP) continue;
      fNumberOfJets->Fill(1);
      Double_t fraction = jetContA->GetFractionSharedPt(jetA);
      //fill the TLorentsVectors with the jet 4-vectors
      Printf("MC jet %p, Fraction %.4f, pT %f, eta %f, phi %f",jetP, fraction, jetA->Pt(), jetA->Eta(), jetA->Phi());
      fhFractionSharedpTA->Fill(jetA->Pt(), fraction);
      //if(fMinFractionShared>0. && fraction>fMinFractionShared) {
      	 fNumberOfJets->Fill(2);
      	 fJetGenSub->SetPtEtaPhiM(jetA->Pt()-GetRhoVal(fContainerArea)*jetA->Area(), jetA->Eta(), jetA->Phi(), jetA->GetSecondOrderSubtracted()); // apply the subtraction to the pT and get the subtracted mass (second order here)
      	 fJetPart1 ->SetPtEtaPhiM(jetP->Pt(), jetP->Eta(), jetP->Phi(), jetP->M());
      	 Printf("Area %.3f > %.3f ?", fJetGenSub->Pt(), fJetGenSubL->Pt());
      	 if(fJetGenSub->Pt() > fJetGenSubL->Pt()){
      	    Printf("YES!");
      	    fJetGenSubL->SetPtEtaPhiM(fJetGenSub->Pt(), fJetGenSub->Eta(), fJetGenSub->Phi(), fJetGenSub->M());
      	    fJetPart1L->SetPtEtaPhiM(fJetPart1->Pt(), fJetPart1->Eta(), fJetPart1->Phi(), fJetPart1->M());
      	 
      	 }
      	 if(!fLeadingJetOnly) fTreeJetsA->Fill();
      //}
      
      
   }
   
   if(fLeadingJetOnly) {
      if(fJetPart1L->Pt() > 0){ // fill only when there is a leading jet
      	 fTreeJetsA->Fill();
      }
   }
   
   if(!fLeadingJetOnly) {
      fTreeJetsC->SetBranchAddress("fJetConSub", &fJetConSub);
      fTreeJetsC->SetBranchAddress("fJetPart2"  ,&fJetPart2);
   } else {
      fTreeJetsC->SetBranchAddress("fJetConSubL",&fJetConSubL);
      fTreeJetsC->SetBranchAddress("fJetPart2L"  ,&fJetPart2L);
   }
   while((jetC = jetContC->GetNextAcceptJet())) {
      
      fJetConSub->SetPtEtaPhiM(0,0,0,0);
      fJetPart2  ->SetPtEtaPhiM(0,0,0,0);
      fNumberOfJets->Fill(3);
      
      AliEmcalJet *jetP = jetC->ClosestJet();
      if(!jetP) continue;
      fNumberOfJets->Fill(4);
      Double_t fraction = jetContC->GetFractionSharedPt(jetC);
      //fill the TLorentsVectors with the jet 4-vectors
      
      //if(fMinFractionShared>0. && fraction>fMinFractionShared) { //to be added eventually
      	 fNumberOfJets->Fill(5);
      	 fJetConSub->SetPtEtaPhiM(jetC->Pt(), jetC->Eta(), jetC->Phi(), jetC->M()); //remeber: in the constituent method the jet is already subtracted, so normal pT and M are taken
      	 fJetPart2 ->SetPtEtaPhiM(jetP->Pt(), jetP->Eta(), jetP->Phi(), jetP->M());
      	 Printf("Const %.3f > %.3f ?", fJetConSub->Pt(), fJetConSubL->Pt());
      	 if(fJetConSub->Pt() > fJetConSubL->Pt()){
      	    Printf("YES!");
      	    fJetConSubL->SetPtEtaPhiM(fJetConSub->Pt(), fJetConSub->Eta(), fJetConSub->Phi(), fJetConSub->M());
      	    fJetPart2L->SetPtEtaPhiM(fJetPart2->Pt(), fJetPart2->Eta(), fJetPart2->Phi(), fJetPart2->M());
      	 
      	 }
      	 
      	 if(!fLeadingJetOnly) fTreeJetsC->Fill();
      //}
      
      
   }
   if(fLeadingJetOnly) {
      if(fJetPart2L->Pt() > 0){ // fill only when there is a leading jet
      	 Printf("Fill with %.3f, %.3f, %.3f, %.3f ", fJetConSubL->Pt(), fJetConSubL->Eta(), fJetConSubL->Phi(), fJetConSubL->M());
      	 fTreeJetsC->Fill();
      }
   }
   
   //delete fJetGenSub; fJetGenSub = 0;
   //delete fJetConSub; fJetConSub = 0;
   //delete fJetPart; fJetPart = 0;

   return kTRUE;
}


//________________________________________________________________________________________________

AliAnalysisTaskPrepareInputForEmbedding::~AliAnalysisTaskPrepareInputForEmbedding() {
   /// Empty destructor.
}

//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::Terminate(Option_t *) 
{
   /// Called once at the end of the analysis.
}
