#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliLog.h"

#include "AliAnalysisTaskEfficiencyBF.h"

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
//
// Authors: Panos Christakoglou, Michael Weber
// 
// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEfficiencyBF)

//________________________________________________________________________
AliAnalysisTaskEfficiencyBF::AliAnalysisTaskEfficiencyBF(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fQAList(0), fOutputList(0), 
  fHistEventStats(0), fHistCentrality(0), fHistNMult(0), 
  fHistGeneratedEtaPtPhiPlus(0), fHistFindableEtaPtPhiPlus(0), 
  fHistReconstructedEtaPtPhiPlus(0), fHistSurvivedEtaPtPhiPlus(0),
  fHistGeneratedEtaPtPhiMinus(0), fHistFindableEtaPtPhiMinus(0), 
  fHistReconstructedEtaPtPhiMinus(0), fHistSurvivedEtaPtPhiMinus(0),
  fHistGeneratedEtaPtPlusControl(0), fHistFindableEtaPtPlusControl(0), 
  fHistReconstructedEtaPtPlusControl(0), fHistSurvivedEtaPtPlusControl(0),
  fHistGeneratedEtaPtMinusControl(0), fHistFindableEtaPtMinusControl(0), 
  fHistReconstructedEtaPtMinusControl(0), fHistSurvivedEtaPtMinusControl(0),
  fHistGeneratedEtaPtPlusPlus(0), fHistFindableEtaPtPlusPlus(0), 
  fHistReconstructedEtaPtPlusPlus(0), fHistSurvivedEtaPtPlusPlus(0),
  fHistGeneratedEtaPtMinusMinus(0), fHistFindableEtaPtMinusMinus(0), 
  fHistReconstructedEtaPtMinusMinus(0), fHistSurvivedEtaPtMinusMinus(0),
  fHistGeneratedEtaPtPlusMinus(0), fHistFindableEtaPtPlusMinus(0), 
  fHistReconstructedEtaPtPlusMinus(0), fHistSurvivedEtaPtPlusMinus(0),
  fHistGeneratedPhiEtaPlusPlus(0), fHistFindablePhiEtaPlusPlus(0), 
  fHistReconstructedPhiEtaPlusPlus(0), fHistSurvivedPhiEtaPlusPlus(0),
  fHistGeneratedPhiEtaMinusMinus(0), fHistFindablePhiEtaMinusMinus(0), 
  fHistReconstructedPhiEtaMinusMinus(0), fHistSurvivedPhiEtaMinusMinus(0),
  fHistGeneratedPhiEtaPlusMinus(0), fHistFindablePhiEtaPlusMinus(0), 
  fHistReconstructedPhiEtaPlusMinus(0), fHistSurvivedPhiEtaPlusMinus(0),
  fESDtrackCuts(0), fAnalysisMode(0), 
  fCentralityEstimator("V0M"), fCentralityPercentileMin(0.0), fCentralityPercentileMax(5.0), 
  fVxMax(3.0), fVyMax(3.0), fVzMax(10.), 
  fMinNumberOfTPCClusters(80), fMaxChi2PerTPCCluster(4.0), fMaxDCAxy(3.0), fMaxDCAz(3.0),
    fMinPt(0.3), fMaxPt(1.5), fMinEta(-0.8), fMaxEta(0.8),fEtaRangeMin(0.0), fEtaRangeMax(1.6), fPtRangeMin(0.3), fPtRangeMax(1.5), fPhiRangeMin(0.0),fPhiRangeMax(360.),fdPhiRangeMax(180.), fEtaBin(100),fdEtaBin(64),fPtBin(49),fPhiBin(100),fdPhiBin(90)   {  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fQAList = new TList();
  fQAList->SetName("QAList");
  fQAList->SetOwner();

  fOutputList = new TList();
  fOutputList->SetName("OutputList");
  fOutputList->SetOwner();

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fQAList->Add(fHistEventStats);

  //ESD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			     20,0.5,20.5);
  fQAList->Add(fHistCentrality);
  
  //multiplicity (good MC tracks)
  TString histName;
  histName = "fHistNMult";
  fHistNMult = new TH1F(histName.Data(), 
			";N_{mult.}",
			200,0,20000);
  fQAList->Add(fHistNMult);
  
  //eta vs pt for MC positives
  fHistGeneratedEtaPtPhiPlus = new TH3D("fHistGeneratedEtaPtPhiPlus",
				     "Generated positive primaries;#eta;p_{T} (GeV/c);#phi",
					fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtPhiPlus);
  fHistFindableEtaPtPhiPlus = new TH3D("fHistFindableEtaPtPhiPlus",
				     "Findable positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistFindableEtaPtPhiPlus);
  fHistReconstructedEtaPtPhiPlus = new TH3D("fHistReconstructedEtaPtPhiPlus",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtPhiPlus);
  fHistSurvivedEtaPtPhiPlus = new TH3D("fHistSurvivedEtaPtPhiPlus",
				     "Survived positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtPhiPlus);

  //eta vs pt for MC negatives
  fHistGeneratedEtaPtPhiMinus = new TH3D("fHistGeneratedEtaPtPhiMinus",
				     "Generated positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtPhiMinus);
  fHistFindableEtaPtPhiMinus = new TH3D("fHistFindableEtaPtPhiMinus",
				     "Findable positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistFindableEtaPtPhiMinus);
  fHistReconstructedEtaPtPhiMinus = new TH3D("fHistReconstructedEtaPtPhiMinus",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtPhiMinus);
  fHistSurvivedEtaPtPhiMinus = new TH3D("fHistSurvivedEtaPtPhiMinus",
				     "Survived positive primaries;#eta;p_{T} (GeV/c);#phi",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax,fPhiBin,fPhiRangeMin,fPhiRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtPhiMinus);

  //eta vs pt for MC positives (control)
  fHistGeneratedEtaPtPlusControl = new TH2F("fHistGeneratedEtaPtPlusControl",
				     "Generated positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtPlusControl);
  fHistFindableEtaPtPlusControl = new TH2F("fHistFindableEtaPtPlusControl",
				     "Findable positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistFindableEtaPtPlusControl);
  fHistReconstructedEtaPtPlusControl = new TH2F("fHistReconstructedEtaPtPlusControl",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtPlusControl);
  fHistSurvivedEtaPtPlusControl = new TH2F("fHistSurvivedEtaPtPlusControl",
				     "Survived positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtPlusControl);

  //eta vs pt for MC negatives (control)
  fHistGeneratedEtaPtMinusControl = new TH2F("fHistGeneratedEtaPtMinusControl",
				     "Generated positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtMinusControl);
  fHistFindableEtaPtMinusControl = new TH2F("fHistFindableEtaPtMinusControl",
				     "Findable positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistFindableEtaPtMinusControl);
  fHistReconstructedEtaPtMinusControl = new TH2F("fHistReconstructedEtaPtMinusControl",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtMinusControl);
  fHistSurvivedEtaPtMinusControl = new TH2F("fHistSurvivedEtaPtMinusControl",
				     "Survived positive primaries;#eta;p_{T} (GeV/c)",
				     fEtaBin,fMinEta,fMaxEta,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtMinusControl);

  //eta vs pt for MC ++
  fHistGeneratedEtaPtPlusPlus = new TH2F("fHistGeneratedEtaPtPlusPlus",
				     "Generated ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtPlusPlus);
  fHistFindableEtaPtPlusPlus = new TH2F("fHistFindableEtaPtPlusPlus",
				     "Findable ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistFindableEtaPtPlusPlus);
  fHistReconstructedEtaPtPlusPlus = new TH2F("fHistReconstructedEtaPtPlusPlus",
				     "Reconstructed ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtPlusPlus);
  fHistSurvivedEtaPtPlusPlus = new TH2F("fHistSurvivedEtaPtPlusPlus",
				     "Survived ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtPlusPlus);

  //eta vs pt for MC --
  fHistGeneratedEtaPtMinusMinus = new TH2F("fHistGeneratedEtaPtMinusMinus",
				     "Generated -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtMinusMinus);
  fHistFindableEtaPtMinusMinus = new TH2F("fHistFindableEtaPtMinusMinus",
				     "Findable -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistFindableEtaPtMinusMinus);
  fHistReconstructedEtaPtMinusMinus = new TH2F("fHistReconstructedEtaPtMinusMinus",
				     "Reconstructed -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtMinusMinus);
  fHistSurvivedEtaPtMinusMinus = new TH2F("fHistSurvivedEtaPtMinusMinus",
				     "Survived -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtMinusMinus);

  //eta vs pt for MC +-
  fHistGeneratedEtaPtPlusMinus = new TH2F("fHistGeneratedEtaPtPlusMinus",
				     "Generated +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistGeneratedEtaPtPlusMinus);
  fHistFindableEtaPtPlusMinus = new TH2F("fHistFindableEtaPtPlusMinus",
				     "Findable +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistFindableEtaPtPlusMinus);
  fHistReconstructedEtaPtPlusMinus = new TH2F("fHistReconstructedEtaPtPlusMinus",
				     "Reconstructed +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistReconstructedEtaPtPlusMinus);
  fHistSurvivedEtaPtPlusMinus = new TH2F("fHistSurvivedEtaPtPlusMinus",
				     "Survived +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     fdEtaBin,fEtaRangeMin,fEtaRangeMax,fPtBin,fPtRangeMin,fPtRangeMax);
  fOutputList->Add(fHistSurvivedEtaPtPlusMinus);

  //=============================//
   //phi vs eta for MC ++
  fHistGeneratedPhiEtaPlusPlus = new TH2F("fHistGeneratedPhiEtaPlusPlus",
				     "Generated ++ primaries;#Delta#phi",
					  fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistGeneratedPhiEtaPlusPlus);
  fHistFindablePhiEtaPlusPlus = new TH2F("fHistFindablePhiEtaPlusPlus",
				     "Findable ++ primaries;#Delta#phi;#Delta#eta",
					 fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistFindablePhiEtaPlusPlus);
  fHistReconstructedPhiEtaPlusPlus = new TH2F("fHistReconstructedPhiEtaPlusPlus",
				     "Reconstructed ++ primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistReconstructedPhiEtaPlusPlus);
  fHistSurvivedPhiEtaPlusPlus = new TH2F("fHistSurvivedPhiEtaPlusPlus",
				     "Survived ++ primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistSurvivedPhiEtaPlusPlus);

  //phi vs eta for MC --
  fHistGeneratedPhiEtaMinusMinus = new TH2F("fHistGeneratedPhiEtaMinusMinus",
				     "Generated -- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistGeneratedPhiEtaMinusMinus);
  fHistFindablePhiEtaMinusMinus = new TH2F("fHistFindablePhiEtaMinusMinus",
				     "Findable -- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistFindablePhiEtaMinusMinus);
  fHistReconstructedPhiEtaMinusMinus = new TH2F("fHistReconstructedPhiEtaMinusMinus",
				     "Reconstructed -- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistReconstructedPhiEtaMinusMinus);
  fHistSurvivedPhiEtaMinusMinus = new TH2F("fHistSurvivedPhiEtaMinusMinus",
				     "Survived -- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistSurvivedPhiEtaMinusMinus);

  //phi vs eta for MC +-
  fHistGeneratedPhiEtaPlusMinus = new TH2F("fHistGeneratedPhiEtaPlusMinus",
				     "Generated +- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistGeneratedPhiEtaPlusMinus);
  fHistFindablePhiEtaPlusMinus = new TH2F("fHistFindablePhiEtaPlusMinus",
				     "Findable +- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistFindablePhiEtaPlusMinus);
  fHistReconstructedPhiEtaPlusMinus = new TH2F("fHistReconstructedPhiEtaPlusMinus",
				     "Reconstructed +- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistReconstructedPhiEtaPlusMinus);
  fHistSurvivedPhiEtaPlusMinus = new TH2F("fHistSurvivedPhiEtaPlusMinus",
				     "Survived +- primaries;#Delta#phi;#Delta#eta",
				     fdPhiBin,fPhiRangeMin,fdPhiRangeMax,fdEtaBin,fEtaRangeMin,fEtaRangeMax);
  fOutputList->Add(fHistSurvivedPhiEtaPlusMinus);
  //=============================//

  fQAList->Print();
  fOutputList->Print();

  PostData(1, fQAList);
  PostData(2, fOutputList);

  TH1::AddDirectory(oldStatus);

}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  // Post output data.
  //ESD analysis
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    AliError("ERROR: Could not retrieve MC event");
    return;
  }
  AliStack* stack = mcEvent->Stack();
  if (!stack) {
    AliError("ERROR: Could not retrieve MC stack");
    return;
  }

  // arrays for 2 particle histograms
  Int_t nMCLabelCounter         = 0;
  const Int_t maxMCLabelCounter = 20000;

  Double_t eta[maxMCLabelCounter];
  Double_t pt[maxMCLabelCounter];
  Double_t phi[maxMCLabelCounter];
  Int_t level[maxMCLabelCounter];
  Int_t charge[maxMCLabelCounter];

  //AliInfo(Form("%d %d",mcEvent->GetNumberOfTracks(),fESD->GetNumberOfTracks()));
  fHistEventStats->Fill(1); //all events
    
  //Centrality stuff
  AliCentrality *centrality = fESD->GetCentrality();
  Int_t nCentrality = 0;
  nCentrality = (Int_t)(centrality->GetCentralityPercentile(fCentralityEstimator.Data())/10.);

  //Printf("Centrality: %lf",centrality->GetCentralityPercentile(fCentralityEstimator.Data()));

  if(centrality->IsEventInCentralityClass(fCentralityPercentileMin,
					  fCentralityPercentileMax,
					  fCentralityEstimator.Data())) {
    fHistEventStats->Fill(2); //triggered + centrality
    fHistCentrality->Fill(nCentrality+1);

    //Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);
  
    if(fAnalysisMode.CompareTo("TPC") == 0 ) {
      const AliESDVertex *vertex = fESD->GetPrimaryVertexTPC();
      if(vertex) {
	if(vertex->GetNContributors() > 0) {
	  if(vertex->GetZRes() != 0) {
	    fHistEventStats->Fill(3); //events with a proper vertex
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4); //analyzed events
		  
		  Int_t nMCParticles = mcEvent->GetNumberOfTracks();
		  TArrayI labelMCArray(nMCParticles);
		  
		  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
		    AliMCParticle *mcTrack = (AliMCParticle*) mcEvent->GetTrack(iTracks);
		    if (!mcTrack) {
		      AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
		      continue;
		    }
		    
		    //exclude particles generated out of the acceptance
		    Double_t vz = mcTrack->Zv();
		    if (TMath::Abs(vz) > 50.) continue;
		   
		    //acceptance
		    if(TMath::Abs(mcTrack->Eta()) > fMaxEta) 
		      continue;
		    if((mcTrack->Pt() > fMaxPt)||(mcTrack->Pt() < fMinPt)) 
		      continue;
		    if((mcTrack->Phi() > fPhiRangeMax)||(mcTrack->Phi() < fPhiRangeMin)) 
		      continue;
		    
		    TParticle* particle = mcTrack->Particle();
		    if(!particle) continue;
		    if(!stack->IsPhysicalPrimary(iTracks)) continue;

		    if(iTracks <= stack->GetNprimary()) {		      
		      Short_t gMCCharge = mcTrack->Charge();
		      Double_t phiRad = particle->Phi();
		      Double_t phiDeg = phiRad*TMath::RadToDeg();	
			
		      if(gMCCharge > 0)
			fHistGeneratedEtaPtPhiPlus->Fill(particle->Eta(),
							 particle->Pt(),
							 phiDeg);
		      else if(gMCCharge < 0)
			fHistGeneratedEtaPtPhiMinus->Fill(particle->Eta(),
							  particle->Pt(),
							  phiDeg);

		      
		      // findable tracks --> DOES NOT WORK????
		      // Loop over Track References
		      Bool_t labelTPC = kTRUE;
		      AliTrackReference* trackRef = 0;

		      for (Int_t iTrackRef = 0; iTrackRef < mcTrack->GetNumberOfTrackReferences(); iTrackRef++) {
			trackRef = mcTrack->GetTrackReference(iTrackRef);
			if(trackRef) {
			  Int_t detectorId = trackRef->DetectorId();
			  if (detectorId == AliTrackReference::kTPC) {
			    labelTPC = kTRUE;
			    break;
			  }
			}
		      }//loop over track references

		      if(labelTPC) {
			labelMCArray.AddAt(iTracks,nMCLabelCounter);

			if(nMCLabelCounter >= maxMCLabelCounter){
			  AliWarning(Form("MC Label Counter > Limit (%d) --> stop loop here",maxMCLabelCounter));
			  break;
			}

			//fill the arrays for 2 particle analysis
			eta[nMCLabelCounter]    = particle->Eta();
			pt[nMCLabelCounter]     = particle->Pt();
			phi[nMCLabelCounter]     = particle->Phi()*TMath::RadToDeg();
			charge[nMCLabelCounter] = gMCCharge;
		
			// findable = generated in this case!

			level[nMCLabelCounter]  = 1;
			nMCLabelCounter += 1;
				

			if(gMCCharge > 0)
			  fHistFindableEtaPtPhiPlus->Fill(particle->Eta(),
							  particle->Pt(),
							  phiDeg);
			else if(gMCCharge < 0)
			  fHistFindableEtaPtPhiMinus->Fill(particle->Eta(),
							   particle->Pt(),
							   phiDeg);
		      }
		    }//primaries
		  }//loop over MC particles
		
		  fHistNMult->Fill(nMCLabelCounter);
		  
		  // not used so far
		  //Float_t dcaXY = 0.0, dcaZ = 0.0;

		  //ESD track loop
		  Int_t nGoodTracks = fESD->GetNumberOfTracks();
		  
		  TArrayI labelArray(nGoodTracks);
		  Int_t labelCounter = 0;
		  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
		    AliESDtrack* track = fESD->GetTrack(iTracks);
		    //AliESDtrack* track = fESDtrackCuts->GetTPCOnlyTrack(fESD,iTracks);
		    if(!track) continue;

		    AliESDtrack *tpcOnlyTrack = new AliESDtrack();
		    
		    if (!track->FillTPCOnlyTrack(*tpcOnlyTrack)) {
		      delete tpcOnlyTrack;
		      continue;
		    }

		    Int_t label = TMath::Abs(track->GetTPCLabel());
		    if(IsLabelUsed(labelArray,label)) continue;
		    labelArray.AddAt(label,labelCounter);
		    labelCounter += 1;

		    Int_t mcGoods = nMCLabelCounter;
		    for (Int_t k = 0; k < mcGoods; k++) {
		      Int_t mcLabel = labelMCArray.At(k);
		    			      
		      if (mcLabel != TMath::Abs(label)) continue;
		      if(mcLabel != label) continue;
		      if(label > stack->GetNtrack()) continue;
		      
		      TParticle *particle = stack->Particle(label);
		      if(!particle) continue;
		      
		      //acceptance
		      if(TMath::Abs(particle->Eta()) > fMaxEta) 
			continue;
		      if((particle->Pt() > fMaxPt)||(particle->Pt() <  fMinPt)) 
			continue;
		      if((particle->Phi() > fPhiRangeMax)||(particle->Phi() < fPhiRangeMin)) 
			continue;

		      if(!stack->IsPhysicalPrimary(label)) continue;
		      
		      if(label <= stack->GetNprimary()) {
			
			// reconstructed
			level[k]  = 2;
			
			Short_t gCharge = track->Charge();
			Double_t phiRad = particle->Phi();
			Double_t phiDeg = phiRad*TMath::RadToDeg();

			if(gCharge > 0)
			  fHistReconstructedEtaPtPhiPlus->Fill(particle->Eta(),
							       particle->Pt(),
							       phiDeg);
			else if(gCharge < 0)
			  fHistReconstructedEtaPtPhiMinus->Fill(particle->Eta(),
								particle->Pt(),
								phiDeg);
			
			// track cuts + analysis kinematic cuts
			if(fESDtrackCuts->AcceptTrack(track) && TMath::Abs(track->Eta()) < fMaxEta && track->Pt() > fMinPt && track->Pt() < fMaxPt ){

			  // survived
			  level[k]  = 3;

			  if(gCharge > 0)
			    fHistSurvivedEtaPtPhiPlus->Fill(particle->Eta(),
							    particle->Pt(),
							    phiDeg);
			  else if(gCharge < 0)
			    fHistSurvivedEtaPtPhiMinus->Fill(particle->Eta(),
							     particle->Pt(),
							     phiDeg);
			  
			}//track cuts
		      }//primary particles
		    }//findable track loop
		  }//ESD track loop

		  labelMCArray.Reset();
		  labelArray.Reset();
		  
		  
		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//Vz resolution
	}//number of contributors
      }//valid vertex
    }//TPC analysis mode
  }//centrality  
      


  // Here comes the 2 particle analysis
  // loop over all good MC particles
  for (Int_t i = 0; i < nMCLabelCounter ; i++) {
    
    // control 1D histograms (charge might be different?)
    if(charge[i] > 0){
      if(level[i] > 0) fHistGeneratedEtaPtPlusControl->Fill(eta[i],pt[i]);
      if(level[i] > 1) fHistReconstructedEtaPtPlusControl->Fill(eta[i],pt[i]);
      if(level[i] > 2) fHistSurvivedEtaPtPlusControl->Fill(eta[i],pt[i]);
    }
    else if(charge[i] < 0){
      if(level[i] > 0) fHistGeneratedEtaPtMinusControl->Fill(eta[i],pt[i]);
      if(level[i] > 1) fHistReconstructedEtaPtMinusControl->Fill(eta[i],pt[i]);
      if(level[i] > 2) fHistSurvivedEtaPtMinusControl->Fill(eta[i],pt[i]);
    }
    
    
    for (Int_t j = i+1; j < nMCLabelCounter ; j++) {
      
      if(charge[i] > 0 && charge[j] > 0 ){
	if(level[i] > 0 && level[j] > 0) {  
	  fHistGeneratedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistGeneratedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
	if(level[i] > 1 && level[j] > 1) {
	  fHistReconstructedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistReconstructedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
	if(level[i] > 2 && level[j] > 2) {
	  fHistSurvivedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistSurvivedPhiEtaPlusPlus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
      }
      
      else if(charge[i] < 0 && charge[j] < 0 ){
	if(level[i] > 0 && level[j] > 0) {
	  fHistGeneratedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);	    
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistGeneratedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));  	    
	}
	if(level[i] > 1 && level[j] > 1) {
	  fHistReconstructedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);   
	  fHistReconstructedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));	   
	}
	if(level[i] > 2 && level[j] > 2) {
	  fHistSurvivedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistSurvivedPhiEtaMinusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}
      }
      
      else if((charge[i] > 0 && charge[j] < 0)||(charge[i] < 0 && charge[j] > 0)){
	if(level[i] > 0 && level[j] > 0) {
	  fHistGeneratedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistGeneratedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));	
	}
	if(level[i] > 1 && level[j] > 1) {
	  fHistReconstructedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  fHistReconstructedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j])); 
	}
	if(level[i] > 2 && level[j] > 2) {
	  fHistSurvivedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	  if (TMath::Abs(phi[i]-phi[j]) < 180)
	  fHistSurvivedPhiEtaPlusMinus->Fill(TMath::Abs(phi[i]-phi[j]),TMath::Abs(eta[i]-eta[j]));
	}	      
      }
    }
  }
  
  //Checking Entries for generated and survived
  /*TH1F* GeneratedEtaPt = (TH1F*)fHistGeneratedEtaPtPlusMinus->ProjectionX("GeneratedEtaPt",5,15);
    Int_t xGeneratedPt = fHistGeneratedEtaPtPlusMinus->GetNbinsX();
    for (Int_t h=1; h < xGeneratedPt+1; h++){
    Double_t binEntriesGenerated = GeneratedEtaPt->GetBinContent(h);
    Printf("binEntriesGenerated: %lf - xGeneratedPt: %d",binEntriesGenerated,h);  
    }
    TH1F* GeneratedPhiEta = (TH1F*)fHistGeneratedPhiEtaPlusMinus->ProjectionY("GeneratedPhiEta",5,15);
    Int_t yGeneratedPhi = fHistGeneratedPhiEtaPlusMinus->GetNbinsY();
    for (Int_t h=1; h < yGeneratedPhi+1; h++){
    Double_t binEntriesGenerated = GeneratedPhiEta->GetBinContent(h);
    Printf("binEntriesGenerated: %lf - yGeneratedPhi: %d",binEntriesGenerated,h);  
    }*/
  
  /*TH1F* SurvivedEtaPt = (TH1F*)fHistSurvivedEtaPtPlusMinus->ProjectionX("SurvivedEtaPt",5,15);
    Int_t xSurvivedPt = fHistSurvivedEtaPtPlusMinus->GetNbinsX();
    for (Int_t h=1; h < xSurvivedPt+1; h++){
    Double_t binEntriesSurvived = SurvivedEtaPt->GetBinContent(h);
    Printf("binEntriesSurvived: %lf - xSurvivedPt: %d",binEntriesSurvived,h);
    }
    TH1F* SurvivedPhiEta = (TH1F*)fHistSurvivedPhiEtaPlusMinus->ProjectionY("SurvivedPhiEta",5,15);
    Int_t ySurvivedPhi = fHistSurvivedPhiEtaPlusMinus->GetNbinsY();
    for (Int_t h=1; h < ySurvivedPhi+1; h++){
    Double_t binEntriesSurvived = SurvivedPhiEta->GetBinContent(h);
    Printf("binEntriesSurvived: %lf - ySurvivedPhi: %d",binEntriesSurvived,h);
    }*/  
}
  //________________________________________________________________________
void AliAnalysisTaskEfficiencyBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//____________________________________________________________________//
Bool_t AliAnalysisTaskEfficiencyBF::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
