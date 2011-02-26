//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD & MC analysis
//  - reconstruction and MonteCarlo output
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEt.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <iostream>
#include "AliAnalysisEtCuts.h"
#include "AliESDtrackCuts.h"
#include "AliVEvent.h"
#include "Rtypes.h"
#include "TString.h"

using namespace std;
ClassImp(AliAnalysisEt);

/* Auxiliary Histogram variables */
Int_t AliAnalysisEt::fgnumOfEtaBins = 16;
Float_t AliAnalysisEt::fgEtaAxis[17]={-0.78, -0.7, -0.58, -0.46, -0.34, -0.22, -0.12, -0.06, -0.0, 0.06, 0.12, 0.22, 0.34, 0.46, 0.58, 0.7, 0.78};

Int_t AliAnalysisEt::fgNumOfEBins = 78;
Float_t AliAnalysisEt::fgEAxis[79]={0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
									1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,11.,
									12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,
									32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.};


AliAnalysisEt::AliAnalysisEt() : AliAnalysisEtCommon()
			       ,fTotEt(0)
			       ,fTotEtAcc(0)
			       ,fTotNeutralEt(0)
			       ,fTotNeutralEtAcc(0)
			       ,fTotChargedEt(0)
			       ,fTotChargedEtAcc(0)
			       ,fMultiplicity(0)
			       ,fChargedMultiplicity(0)
			       ,fNeutralMultiplicity(0)
			       ,fBaryonEt(0)
			       ,fAntiBaryonEt(0)
			       ,fMesonEt(0)
			       ,fProtonEt(0)
			       ,fPionEt(0)
			       ,fChargedKaonEt(0)
			       ,fMuonEt(0)
			       ,fElectronEt(0)
			       ,fNeutronEt(0)
			       ,fAntiNeutronEt(0)
			       ,fGammaEt(0)
			       ,fProtonEtAcc(0)
			       ,fPionEtAcc(0)
			       ,fChargedKaonEtAcc(0)
			       ,fMuonEtAcc(0)
			       ,fElectronEtAcc(0)
			       ,fEnergyDeposited(0)
			       ,fEnergyTPC(0)
			       ,fCharge(0)
			       ,fParticlePid(0)
			       ,fPidProb(0)
			       ,fTrackPassedCut(kFALSE)
			       ,fEtaCut(0)
			       ,fEtaCutAcc(0)
			       ,fPhiCutAccMin(0)
			       ,fPhiCutAccMax(0)
			       ,fDetectorRadius(0)
			       ,fClusterEnergyCut(0) 
			       ,fSingleCellEnergyCut(0)
			       ,fHistEt(0)
			       ,fHistChargedEt(0)
			       ,fHistNeutralEt(0)
			       ,fHistEtAcc(0)
			       ,fHistChargedEtAcc(0)
			       ,fHistNeutralEtAcc(0)
			       ,fHistMult(0)
			       ,fHistChargedMult(0)
			       ,fHistNeutralMult(0)
			       ,fHistPhivsPtPos(0)
			       ,fHistPhivsPtNeg(0)
			       ,fHistBaryonEt(0)
			       ,fHistAntiBaryonEt(0)
			       ,fHistMesonEt(0)
			       ,fHistProtonEt(0)
			       ,fHistPionEt(0)
			       ,fHistChargedKaonEt(0)
			       ,fHistMuonEt(0)
			       ,fHistElectronEt(0)
			       ,fHistNeutronEt(0)
			       ,fHistAntiNeutronEt(0)
			       ,fHistGammaEt(0)
			       ,fHistProtonEtAcc(0)
			       ,fHistPionEtAcc(0)
			       ,fHistChargedKaonEtAcc(0)
			       ,fHistMuonEtAcc(0)
			       ,fHistElectronEtAcc(0)
			       ,fHistTMDeltaR(0)
			       ,fTree(0)
			       ,fTreeDeposit(0)
			       ,fCentrality(0)
{}

AliAnalysisEt::~AliAnalysisEt()
{//Destructor
  if(fTreeDeposit){
    fTreeDeposit->Clear();
    delete fTreeDeposit; // optional TTree
  }
  if(fTree){
    fTree->Clear();
    delete fTree; // optional TTree
  }
  delete fHistEt; //Et spectrum
  delete fHistChargedEt; //Charged Et spectrum 
  delete fHistNeutralEt; //Neutral Et spectrum
  delete fHistEtAcc; //Et in acceptance
  delete fHistChargedEtAcc; //Charged Et in acceptance
  delete fHistNeutralEtAcc; //Et in acceptance
  delete fHistMult; //Multiplicity
  delete fHistChargedMult; //Charged multiplicity
  delete fHistNeutralMult; //Neutral multiplicity
  delete fHistPhivsPtPos; //phi vs pT plot for positive tracks
  delete fHistPhivsPtNeg; //phi vs pT plot for negative tracks
  delete fHistBaryonEt; /** Et of identified baryons */    
  delete fHistAntiBaryonEt; /** Et of identified anti-baryons */
  delete fHistMesonEt; /** Et of identified mesons */
  delete fHistProtonEt; /** Et of identified protons */
  delete fHistPionEt; /** Et of identified protons */
  delete fHistChargedKaonEt; /** Et of identified charged kaons */
  delete fHistMuonEt; /** Et of identified muons */
  delete fHistElectronEt; /** Et of identified electrons */
  delete fHistNeutronEt; /** Et of neutrons (MC only for now) */
  delete fHistAntiNeutronEt; /** Et of anti-neutrons (MC only for now) */
  delete fHistGammaEt; /** Et of gammas (MC only for now) */
  delete fHistProtonEtAcc; /** Et of identified protons in calorimeter acceptance */    
  delete fHistPionEtAcc; /** Et of identified protons in calorimeter acceptance */    
  delete fHistChargedKaonEtAcc; /** Et of identified charged kaons in calorimeter acceptance */    
  delete fHistMuonEtAcc; /** Et of identified muons in calorimeter acceptance */
  delete fHistElectronEtAcc; /** Et of identified electrons in calorimeter acceptance */
  delete fHistTMDeltaR; /* Track matching plots; Rec only for now */
}

void AliAnalysisEt::FillOutputList(TList *list)
{ // histograms to be added to output
    list->Add(fHistEt);
    list->Add(fHistChargedEt);
    list->Add(fHistNeutralEt);

    list->Add(fHistEtAcc);
    list->Add(fHistChargedEtAcc);
    list->Add(fHistNeutralEtAcc);

    list->Add(fHistMult);
    list->Add(fHistChargedMult);
    list->Add(fHistNeutralMult);

    list->Add(fHistPhivsPtPos);
    list->Add(fHistPhivsPtNeg);

    list->Add(fHistBaryonEt);
    list->Add(fHistAntiBaryonEt);
    list->Add(fHistMesonEt);

    list->Add(fHistProtonEt);
    list->Add(fHistPionEt);
    list->Add(fHistChargedKaonEt);
    list->Add(fHistMuonEt);
    list->Add(fHistElectronEt);
    
    list->Add(fHistNeutronEt);
    list->Add(fHistAntiNeutronEt);
    list->Add(fHistGammaEt);
    
    list->Add(fHistProtonEtAcc);
    list->Add(fHistPionEtAcc);
    list->Add(fHistChargedKaonEtAcc);
    list->Add(fHistMuonEtAcc);
    list->Add(fHistElectronEtAcc);

    list->Add(fHistTMDeltaR);
	
    if (fCuts) {
      if (fCuts->GetHistMakeTree()) {
	list->Add(fTree);
      }
    }
    list->Add(fTreeDeposit);

}

void AliAnalysisEt::Init()
{// clear variables, set up cuts and PDG info
  AliAnalysisEtCommon::Init();
  ResetEventValues();
}

void AliAnalysisEt::CreateHistograms()
{ // create histograms..
  // histogram binning for E_T, p_T and Multiplicity: defaults for p+p
  Int_t nbinsEt = 1000;
  Double_t minEt = 0.0001;
  Double_t maxEt = 100;
  Int_t nbinsPt = 200;
  Double_t minPt = 0;
  Double_t maxPt = 20;
  Int_t nbinsMult = 200;
  Double_t minMult = -0.5; // offset -0.5 to have integer bins centered around 0
  Double_t maxMult = nbinsMult + minMult; // 1 bin per integer value

  // see if we should change histogram limits etc, and possibly create a tree
  if (fCuts) {
    //if (fCuts->GetHistMakeTree()) {
      CreateTrees();
    //}

    nbinsMult = fCuts->GetHistNbinsMult();
    minMult = fCuts->GetHistMinMult();
    maxMult = fCuts->GetHistMaxMult();

    nbinsEt = fCuts->GetHistNbinsTotEt();
    minEt = fCuts->GetHistMinTotEt();
    maxEt = fCuts->GetHistMaxTotEt();

    nbinsPt = fCuts->GetHistNbinsParticlePt();
    minPt = fCuts->GetHistMinParticlePt();
    maxPt = fCuts->GetHistMaxParticlePt();
  }

    TString histname = "fHistEt" + fHistogramNameSuffix;
    fHistEt = new TH1F(histname.Data(), "Total E_{T} Distribution", nbinsEt, minEt, maxEt);
    fHistEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistChargedEt" + fHistogramNameSuffix;
    fHistChargedEt = new TH1F(histname.Data(), "Total Charged E_{T} Distribution", nbinsEt, minEt, maxEt);
    fHistChargedEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistChargedEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistNeutralEt" + fHistogramNameSuffix;
    fHistNeutralEt = new TH1F(histname.Data(), "Total Neutral E_{T} Distribution", nbinsEt, minEt, maxEt);
    fHistNeutralEt->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistNeutralEt->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistEtAcc" + fHistogramNameSuffix;
    fHistEtAcc = new TH1F(histname.Data(), "Total E_{T} Distribution in Acceptance", nbinsEt, minEt, maxEt);
    fHistEtAcc->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistEtAcc->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistChargedEtAcc" + fHistogramNameSuffix;
    fHistChargedEtAcc = new TH1F(histname.Data(), "Total Charged E_{T} Distribution in Acceptance", nbinsEt, minEt, maxEt);
    fHistChargedEtAcc->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistChargedEtAcc->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");

    histname = "fHistNeutralEtAcc" + fHistogramNameSuffix;
    fHistNeutralEtAcc = new TH1F(histname.Data(), "Total Neutral E_{T} Distribution in Acceptance", nbinsEt, minEt, maxEt);
    fHistNeutralEtAcc->GetXaxis()->SetTitle("E_{T} (GeV/c^{2})");
    fHistNeutralEtAcc->GetYaxis()->SetTitle("dN/dE_{T} (c^{2}/GeV)");
    std::cout << histname << std::endl;

	histname = "fHistMult" + fHistogramNameSuffix;
    fHistMult = new TH1F(histname.Data(), "Total Multiplicity", nbinsMult, minMult, maxMult);
    fHistMult->GetXaxis()->SetTitle("N");
    fHistMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistChargedMult" + fHistogramNameSuffix;
    fHistChargedMult = new TH1F(histname.Data(), "Charged Multiplicity", nbinsMult, minMult, maxMult);
    fHistChargedMult->GetXaxis()->SetTitle("N");
    fHistChargedMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistNeutralMult" + fHistogramNameSuffix;
    fHistNeutralMult = new TH1F(histname.Data(), "Neutral Multiplicity", nbinsMult, minMult, maxMult);
    fHistNeutralMult->GetXaxis()->SetTitle("N");
    fHistNeutralMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistPhivsPtPos" + fHistogramNameSuffix;
    fHistPhivsPtPos = new TH2F(histname.Data(), "Phi vs pT of positively charged tracks hitting the calorimeter", 	200, 0, 2*TMath::Pi(), nbinsPt, minPt, maxPt);

    histname = "fHistPhivsPtNeg" + fHistogramNameSuffix;
    fHistPhivsPtNeg = new TH2F(histname.Data(), "Phi vs pT of negatively charged tracks hitting the calorimeter", 	200, 0, 2*TMath::Pi(), nbinsPt, minPt, maxPt);

    histname = "fHistBaryonEt" + fHistogramNameSuffix;
    fHistBaryonEt = new TH1F(histname.Data(), "E_{T} for baryons",  nbinsEt, minEt, maxEt);

    histname = "fHistAntiBaryonEt" + fHistogramNameSuffix;
    fHistAntiBaryonEt = new TH1F(histname.Data(), "E_{T} for anti baryons",  nbinsEt, minEt, maxEt);

    histname = "fHistMesonEt" + fHistogramNameSuffix;
    fHistMesonEt = new TH1F(histname.Data(), "E_{T} for mesons",  nbinsEt, minEt, maxEt);

    histname = "fHistProtonEt" + fHistogramNameSuffix;
    fHistProtonEt = new TH1F(histname.Data(), "E_{T} for (anti-)protons", nbinsEt, minEt, maxEt);

    histname = "fHistPionEt" + fHistogramNameSuffix;
    fHistPionEt = new TH1F(histname.Data(), "E_{T} for #pi^+/#pi^-", nbinsEt, minEt, maxEt);

    histname = "fHistKaonEt" + fHistogramNameSuffix;
    fHistChargedKaonEt = new TH1F(histname.Data(), "E_{T} for charged kaons", nbinsEt, minEt, maxEt);

    histname = "fHistMuonEt" + fHistogramNameSuffix;
    fHistMuonEt = new TH1F(histname.Data(), "E_{T} for muons", nbinsEt, minEt, maxEt);

    histname = "fHistElectronEt" + fHistogramNameSuffix;
    fHistElectronEt = new TH1F(histname.Data(), "E_{T} for electrons/positrons", nbinsEt, minEt, maxEt);

    histname = "fHistNeutronEt" + fHistogramNameSuffix;
    fHistNeutronEt = new TH1F(histname.Data(), "E_{T} for neutrons", nbinsEt, minEt, maxEt);

    histname = "fHistAntiNeutronEt" + fHistogramNameSuffix;
    fHistAntiNeutronEt = new TH1F(histname.Data(), "E_{T} for anti-neutrons", nbinsEt, minEt, maxEt);

    histname = "fHistGammaEt" + fHistogramNameSuffix;
    fHistGammaEt = new TH1F(histname.Data(), "E_{T} for gammas", nbinsEt, minEt, maxEt);

    histname = "fHistProtonEtAcc" + fHistogramNameSuffix;
    fHistProtonEtAcc = new TH1F(histname.Data(), "E_{T} for (anti-)protons in calorimeter acceptance", nbinsEt, minEt, maxEt);

    histname = "fHistPionEtAcc" + fHistogramNameSuffix;
    fHistPionEtAcc = new TH1F(histname.Data(), "E_{T} for #pi^+/#pi^- in calorimeter acceptance", nbinsEt, minEt, maxEt);

    histname = "fHistKaonEtAcc" + fHistogramNameSuffix;
    fHistChargedKaonEtAcc = new TH1F(histname.Data(), "E_{T} for charged kaons in calorimeter acceptance", nbinsEt, minEt, maxEt);

    histname = "fHistMuonEtAcc" + fHistogramNameSuffix;
    fHistMuonEtAcc = new TH1F(histname.Data(), "E_{T} for muons in calorimeter acceptance", nbinsEt, minEt, maxEt);

    histname = "fHistElectronEtAcc" + fHistogramNameSuffix;
    fHistElectronEtAcc = new TH1F(histname.Data(), "E_{T} for electrons/positrons in calorimeter acceptance", nbinsEt, minEt, maxEt);

    //
    histname = "fHistTMDeltaR" + fHistogramNameSuffix;
    fHistTMDeltaR = new TH1F(histname.Data(), "#Delta R for calorimeter clusters", 200, 0, 50);
	
}

TH2F* AliAnalysisEt::CreateEtaEHisto2D(TString name, TString title, TString ztitle)
{     //creates a 2-d histogram in eta and E and adds it to the list of histograms to be saved
	TString *histoname   = new TString();
	TString *histotitle   = new TString();
	TString *zaxistitle   = new TString();
	
	histoname->Append(name);
	histotitle->Append(title);
	zaxistitle->Append(ztitle);
	
	TH2F *histo = new TH2F(histoname->Data(),histotitle->Data(),fgNumOfEBins, fgEAxis, fgnumOfEtaBins, fgEtaAxis);
	histo->SetYTitle("#eta");
	histo->SetXTitle("E (GeV)");
	histo->SetZTitle(zaxistitle->Data());
	histo->Sumw2();
	//fhistoList->Add(histo);
	delete histoname;
	delete histotitle;
    delete zaxistitle;
	
	return histo; 
}


void AliAnalysisEt::CreateTrees()
{ // create tree..
  TString treename = "fTree" + fHistogramNameSuffix;
  if(fCuts->GetHistMakeTree())
  {
  
    fTree = new TTree(treename, treename);
    fTree->Branch("fTotEt",&fTotEt,"fTotEt/D");
    fTree->Branch("fTotEtAcc",&fTotEtAcc,"fTotEtAcc/D");
    fTree->Branch("fTotNeutralEt",&fTotNeutralEt,"fTotNeutralEt/D");
    fTree->Branch("fTotNeutralEtAcc",&fTotNeutralEtAcc,"fTotNeutralEtAcc/D");
    fTree->Branch("fTotChargedEt",&fTotChargedEt,"fTotChargedEt/D");
    fTree->Branch("fTotChargedEtAcc",&fTotChargedEtAcc,"fTotChargedEtAcc/D");
    fTree->Branch("fMultiplicity",&fMultiplicity,"fMultiplicity/I");
    fTree->Branch("fChargedMultiplicity",&fChargedMultiplicity,"fChargedMultiplicity/I");
    fTree->Branch("fNeutralMultiplicity",&fNeutralMultiplicity,"fNeutralMultiplicity/I");
    fTree->Branch("fBaryonEt",&fBaryonEt,"fBaryonEt/D");
    fTree->Branch("fAntiBaryonEt",&fAntiBaryonEt,"fAntiBaryonEt/D");
    fTree->Branch("fMesonEt",&fMesonEt,"fMesonEt/D");
    fTree->Branch("fProtonEt",&fProtonEt,"fProtonEt/D");
    fTree->Branch("fChargedKaonEt",&fChargedKaonEt,"fChargedKaonEt/D");
    fTree->Branch("fMuonEt",&fMuonEt,"fMuonEt/D");
    fTree->Branch("fElectronEt",&fElectronEt,"fElectronEt/D");
    fTree->Branch("fProtonEtAcc",&fProtonEtAcc,"fProtonEtAcc/D");
    fTree->Branch("fChargedKaonEtAcc",&fChargedKaonEtAcc,"fChargedKaonEtAcc/D");
    fTree->Branch("fMuonEtAcc",&fMuonEtAcc,"fMuonEtAcc/D");
    fTree->Branch("fElectronEtAcc",&fElectronEtAcc,"fElectronEtAcc/D");
  }
  
  if(fCuts->GetHistMakeTreeDeposit())
  {
    treename = "fTreeDeposit" + fHistogramNameSuffix;
    fTreeDeposit = new TTree(treename, treename);
  
    fTreeDeposit->Branch("fEnergyDeposited", &fEnergyDeposited, "fEnergyDeposited/F");
    fTreeDeposit->Branch("fEnergyTPC", &fEnergyTPC, "fEnergyTPC/F");
    fTreeDeposit->Branch("fCharge", &fCharge, "fCharge/S");
    fTreeDeposit->Branch("fParticlePid", &fParticlePid, "fParticlePid/S");
    fTreeDeposit->Branch("fPidProb", &fPidProb, "fPidProb/F");
    fTreeDeposit->Branch("fTrackPassedCut", &fTrackPassedCut, "fTrackPassedCut/B");
 
  }

  return;
}
void AliAnalysisEt::FillHistograms()
{ // fill histograms..
    fHistEt->Fill(fTotEt);
    fHistChargedEt->Fill(fTotChargedEt);
    fHistNeutralEt->Fill(fTotNeutralEt);

    fHistEtAcc->Fill(fTotEtAcc);
    fHistChargedEtAcc->Fill(fTotChargedEtAcc);
    fHistNeutralEtAcc->Fill(fTotNeutralEtAcc);

    fHistMult->Fill(fMultiplicity);
    fHistChargedMult->Fill(fChargedMultiplicity);
    fHistNeutralMult->Fill(fNeutralMultiplicity);

    fHistBaryonEt->Fill(fBaryonEt);
    fHistAntiBaryonEt->Fill(fAntiBaryonEt);
    fHistMesonEt->Fill(fMesonEt);

    fHistProtonEt->Fill(fProtonEt);
    fHistPionEt->Fill(fPionEt);
    fHistChargedKaonEt->Fill(fChargedKaonEt);
    fHistMuonEt->Fill(fMuonEt);
    fHistElectronEt->Fill(fElectronEt);
    fHistNeutronEt->Fill(fNeutronEt);
    fHistAntiNeutronEt->Fill(fAntiNeutronEt);
    fHistGammaEt->Fill(fGammaEt);
    
    fHistProtonEtAcc->Fill(fProtonEtAcc);
    fHistPionEtAcc->Fill(fPionEtAcc);
    fHistChargedKaonEtAcc->Fill(fChargedKaonEtAcc);
    fHistMuonEtAcc->Fill(fMuonEtAcc);
    fHistElectronEtAcc->Fill(fElectronEtAcc);

    if (fCuts) {
      if (fCuts->GetHistMakeTree()) {
	fTree->Fill();
      }
    }

}

Int_t AliAnalysisEt::AnalyseEvent(AliVEvent *event)
{ //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
  AliAnalysisEtCommon::AnalyseEvent(event);
  ResetEventValues();
  return 0;
}

void AliAnalysisEt::ResetEventValues()
{ // clear
  AliAnalysisEtCommon::ResetEventValues();
  fTotEt = 0;
  fTotEtAcc = 0;
  fTotNeutralEt = 0;
  fTotNeutralEtAcc = 0;
  fTotChargedEt  = 0;
  fTotChargedEtAcc = 0;
  fMultiplicity = 0;
  fChargedMultiplicity = 0;
  fNeutralMultiplicity = 0;
  fBaryonEt = 0;
  fAntiBaryonEt = 0;
  fMesonEt = 0;
  fProtonEt = 0;
  fPionEt = 0;
  fChargedKaonEt = 0;
  fMuonEt = 0;
  fElectronEt = 0;
  fNeutronEt = 0;
  fAntiNeutronEt = 0;
  fGammaEt = 0;
  fProtonEtAcc = 0;
  fPionEtAcc = 0;
  fChargedKaonEtAcc = 0;
  fMuonEtAcc = 0;
  fElectronEtAcc = 0;
  return;
}

