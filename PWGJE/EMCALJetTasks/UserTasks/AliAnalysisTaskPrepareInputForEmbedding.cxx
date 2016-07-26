#include <TTree.h>
#include <TLorentzVector.h>

#include <AliLog.h>
#include <AliEmcalJet.h>
#include <AliAnalysisManager.h>
#include <TFile.h>
#include <THnSparse.h>
#include "AliJetContainer.h"
#include "AliAnalysisTaskPrepareInputForEmbedding.h"


ClassImp(AliAnalysisTaskPrepareInputForEmbedding)

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding() : AliAnalysisTaskEmcalJet("AliAnalysisTaskPrepareInputForEmbedding", kTRUE),
fContainer(0),
fMinFractionShared(-1),
fLeadingJetOnly(0),
fHardCoreTag(0),
fTreeJets(0),
fJetDet(0),
fJetPart(0),
fJetDetL(0),
fJetPartL(0),
fNumberOfJets(0),
fhFractionSharedpT(0),
fNAccJets(0),
fhResponse(0)
{
   /// default constructor
   
   SetMakeGeneralHistograms(kTRUE);
   
   DefineOutput(2, TTree::Class());
}

//________________________________________________________________________________________________
AliAnalysisTaskPrepareInputForEmbedding::AliAnalysisTaskPrepareInputForEmbedding(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
fContainer(0),
fMinFractionShared(-1),
fLeadingJetOnly(0),
fHardCoreTag(0),
fTreeJets(0),
fJetDet(0),
fJetPart(0),
fJetDetL(0),
fJetPartL(0),
fNumberOfJets(0),
fhFractionSharedpT(0),
fNAccJets(0),
fhResponse(0)

{
   /// standard constructor
   
   SetMakeGeneralHistograms(kTRUE);
   AliInfo("standard constructor");
   DefineOutput(2, TTree::Class());
}
//________________________________________________________________________________________________

void AliAnalysisTaskPrepareInputForEmbedding::UserCreateOutputObjects(){
   /// Create output
   AliAnalysisTaskEmcal::UserCreateOutputObjects();

   fJetDet   = new TLorentzVector();
   fJetPart  = new TLorentzVector();
   fJetDetL  = new TLorentzVector();
   fJetPartL = new TLorentzVector();
   
   fTreeJets = new TTree(Form("fTreeJet%s%s", fLeadingJetOnly ? "Lj" : "", fHardCoreTag ? "HC" : ""), "fTreeJet");
   //Important! 'dot' after the name needed -- see TTree doc! otherwise the objects in the branches have all the same name and the output cannot be properly retrieved (only with the TBrowser)
   if(fLeadingJetOnly){
      fTreeJets->Branch("fJetDetL.", fJetDetL);
      fTreeJets->Branch("fJetPartL.",fJetPartL);
   } else {
      fTreeJets->Branch("fJetDet.", fJetDet);
      fTreeJets->Branch("fJetPart.",fJetPart);
   }
   //fOutput->Add(fTreeJets);
   PostData(2, fTreeJets);
   
   //check
   const Int_t nBinsPt          = 40;
   const Int_t nBinsFraction    = 101;
   const Int_t nBinsfineM       = 200;
   const Int_t nBinsfinePt      = 200;
   
   const Double_t minPt       = 0.;
   const Double_t maxPt       = 200.;
   const Double_t minM        = 0.;
   const Double_t maxM        = 50.;
   const Double_t minFraction =  -0.005;
   const Double_t maxFraction =  1.005;
   
   fNumberOfJets = new TH1F("fNumberOfJets", "Number of Jets", 6, -0.5, 5.5);
   fOutput->Add(fNumberOfJets);
   
   fhFractionSharedpT = new TH2F("fhFractionSharedpT", "Reco/particle shared #it{p}_{T} fraction;#it{p}_{T} (GeV/c); Fraction shared #it{p}_{T}", nBinsPt, minPt, maxPt, nBinsFraction, minFraction, maxFraction);
   fOutput->Add(fhFractionSharedpT);
   
   fNAccJets = new TH1F("fNAccJets","fNAccJets;N/ev",11,-0.5, 9.5);
   fOutput->Add(fNAccJets);
   
   const Int_t dim = 5;
   Int_t nBins[dim] = {nBinsfineM, nBinsfineM, nBinsfinePt, nBinsfinePt, 40};
   Double_t xmin[dim]  = {minM, minM, minPt, minPt, 0.};
   Double_t xmax[dim]  = {maxM, maxM, maxPt, maxPt, 2.};
   TString hsptitle = "Mass-pT response; #it{M}_{part}; #it{M}_{det}; #it{p}_{T,part}; #it{p}_{T,det}; #it{N}_{const}^{det}/#it{N}_{const}^{part}";
   fhResponse = new THnSparseF("hResponse", hsptitle.Data(), dim, nBins, xmin, xmax);
   fOutput->Add(fhResponse);
   
}


//________________________________________________________________________________________________
Bool_t AliAnalysisTaskPrepareInputForEmbedding::Run(){
   /// Run code before FillHistograms()
   
   //TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
   //if (!tree) {
   //   AliError(Form("%s : No current tree!",GetName()));
   //   return kFALSE;
   //}
   //TFile *curfile = tree->GetCurrentFile();
   //if (!curfile) {
   //   AliError(Form("%s : No current file!",GetName()));
   //   return kFALSE;
   //}
   //PythiaInfoFromFile(curfile->GetName(), fXsec, fNtrials, fPthardB);
   
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
      //jet hard core tagging (if requested)
      if(fHardCoreTag && jet->GetTagStatus()<1 && !jet->GetTaggedJet()) continue;

      fJetDet  ->SetPtEtaPhiM(0,0,0,0);
      fJetPart ->SetPtEtaPhiM(0,0,0,0);
      
      count++;
      fNumberOfJets->Fill(0.);
      
      AliEmcalJet *jetP = jet->ClosestJet();
      if(!jetP) continue;
      fNumberOfJets->Fill(1.);
      Double_t fraction = jetCont->GetFractionSharedPt(jet);
      //fill the TLorentsVectors with the jet 4-vectors
      //Printf("MC jet %p, Fraction %.4f, pT %f, eta %f, phi %f, m  %f",jetP, fraction, jet->Pt(), jet->Eta(), jet->Phi(), jet->M());
      fhFractionSharedpT->Fill(jet->Pt(), fraction);
      if(fMinFractionShared<0. || fraction>fMinFractionShared) {
      	 fNumberOfJets->Fill(2.);
      	 fJetDet ->SetPtEtaPhiM(jet->Pt(), jet->Eta(), jet->Phi(), jet->M());
      	 fJetPart->SetPtEtaPhiM(jetP->Pt(), jetP->Eta(), jetP->Phi(), jetP->M());
      	
      	 if(fJetDet->Pt() > fJetDetL->Pt()){
      	    fJetDetL->SetPtEtaPhiM(fJetDet->Pt(), fJetDet->Eta(), fJetDet->Phi(), fJetDet->M());
      	    fJetPartL->SetPtEtaPhiM(fJetPart->Pt(), fJetPart->Eta(), fJetPart->Phi(), fJetPart->M());
      	 
      	 }
      	 if(!fLeadingJetOnly) fTreeJets->Fill();
      	 
      	 Double_t response[5] = {jetP->M(), jet->M(), jetP->Pt(), jet->Pt(), (Double_t)jet->GetNumberOfConstituents()/(Double_t)jetP->GetNumberOfConstituents()};
      	 fhResponse->Fill(response);
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
