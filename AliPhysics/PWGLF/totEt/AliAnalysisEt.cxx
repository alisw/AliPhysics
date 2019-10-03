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
#include "TH1I.h" 
#include "TTree.h"
#include <iostream>
#include "AliAnalysisEtCuts.h"
#include "AliESDtrackCuts.h"
#include "AliESDCaloCluster.h"
#include "AliVEvent.h"
#include "Rtypes.h"
#include "TString.h"
#include "AliCentrality.h"
#include "AliAnalysisEtSelector.h"
#include "AliAnalysisEtTrackMatchCorrections.h"
#include "AliAnalysisEtRecEffCorrection.h"
#include "TFile.h"
#include "TVector3.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h" 
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

using namespace std;
ClassImp(AliAnalysisEt);

/* Auxiliary Histogram variables */
Int_t AliAnalysisEt::fgnumOfEtaBins = 16;
Float_t AliAnalysisEt::fgEtaAxis[17]={-0.78, -0.7, -0.58, -0.46, -0.34, -0.22, -0.12, -0.06, -0.0, 0.06, 0.12, 0.22, 0.34, 0.46, 0.58, 0.7, 0.78};

Int_t AliAnalysisEt::fgNumOfPtBins = 111;
Float_t AliAnalysisEt::fgPtAxis[117]=
   {0.0,0.01,0.02,0.03,0.04, 0.05, 0.06,0.07,0.08,0.09, 0.10,0.11, .12,0.13, .14,0.15, .16,0.17, .18,0.19,
	0.2, .22, .24, .26, .28, 0.30, 0.32, .34, .36, .38, 0.40, .42, .44, .46, .48,
	0.5, .52, .54, .56, .58, 0.60, 0.62, .64, .66, .68, 0.70, .72, .74, .76, .78,
	.80, .82, .84, .86, .88, 0.90, 0.92, .94, .96, .98, 1.00,1.05, 1.1,1.15, 1.2,
	1.25, 1.3,1.35,1.40,1.45, 1.50, 1.55, 1.6,1.65, 1.7, 1.75, 1.8,1.85, 1.9,1.95,
	2.0, 2.2, 2.4, 2.6, 2.8, 3.00, 3.20, 3.4, 3.6, 3.8, 4.00, 4.2, 4.4, 4.6, 4.8,
	5.0, 5.5, 6.0, 6.5, 7.0, 7.50, 8.00, 8.5, 9.0, 9.5, 10.0,12.0,14.0,16.0,18.0,
	20.0,25.0,30.0,35.0,40.0, 45.0, 50.0}; 

Int_t AliAnalysisEt::fgNumOfEBins = 78;
Float_t AliAnalysisEt::fgEAxis[79]={0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
									1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,11.,
									12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,
									32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.};

Int_t AliAnalysisEt::fgNumOfRBins = 47;
Float_t AliAnalysisEt::fgRAxis[48]={-2.,-1.,0.,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,0.01,
									0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.020,0.022,0.024,0.026,0.028,0.03,0.032,0.034,0.036,0.038,0.04,0.042,0.044,0.046,0.048,0.05};


AliAnalysisEt::AliAnalysisEt() : AliAnalysisEtCommon()
			       ,fTmCorrections(0)
			       ,fReCorrections(0)
			       ,fEventSummaryTree(0)
			       ,fAcceptedTree(0)
			       ,fDepositTree(0)
			       ,fTotEt(0)
			       ,fTotEtAcc(0)
			       ,fTotNeutralEt(0)
			       ,fTotNeutralEtAcc(0)
			       ,fTotChargedEt(0)
			       ,fTotChargedEtAcc(0)
			       ,fMultiplicity(0)
			       ,fChargedMultiplicity(0)
			       ,fNeutralMultiplicity(0)
			       ,fProtonEt(0)
			       ,fAntiProtonEt(0)
			       ,fNeutronEt(0)
			       ,fAntiNeutronEt(0)
			       ,fPi0Et(0)
			       ,fPiPlusEt(0)
			       ,fPiMinusEt(0)
			       ,fKPlusEt(0)
			       ,fKMinusEt(0)
			       ,fK0sEt(0)
			       ,fK0lEt(0)
			       ,fMuMinusEt(0)
			       ,fMuPlusEt(0)
			       ,fEMinusEt(0)
			       ,fEPlusEt(0)
			       ,fGammaEt(0)
			       ,fProtonRemovedEt(0)
			       ,fAntiProtonRemovedEt(0)
			       ,fNeutronRemovedEt(0)
			       ,fAntiNeutronRemovedEt(0)
			       ,fPi0RemovedEt(0)
			       ,fPiPlusRemovedEt(0)
			       ,fPiMinusRemovedEt(0)
			       ,fKPlusRemovedEt(0)
			       ,fKMinusRemovedEt(0)
			       ,fK0sRemovedEt(0)
			       ,fK0lRemovedEt(0)
			       ,fMuMinusRemovedEt(0)
			       ,fMuPlusRemovedEt(0)
			       ,fEMinusRemovedEt(0)
			       ,fEPlusRemovedEt(0)
			       ,fGammaRemovedEt(0)
			       ,fProtonMult(0)
			       ,fAntiProtonMult(0)
			       ,fNeutronMult(0)
			       ,fAntiNeutronMult(0)
			       ,fPi0Mult(0)
			       ,fPiPlusMult(0)
			       ,fPiMinusMult(0)
			       ,fKPlusMult(0)
			       ,fKMinusMult(0)
			       ,fK0sMult(0)
			       ,fK0lMult(0)
			       ,fMuMinusMult(0)
			       ,fMuPlusMult(0)
			       ,fEMinusMult(0)
			       ,fEPlusMult(0)
			       ,fGammaMult(0)
			       ,fProtonRemovedMult(0)
			       ,fAntiProtonRemovedMult(0)
			       ,fNeutronRemovedMult(0)
			       ,fAntiNeutronRemovedMult(0)
			       ,fPi0RemovedMult(0)
			       ,fPiPlusRemovedMult(0)
			       ,fPiMinusRemovedMult(0)
			       ,fKPlusRemovedMult(0)
			       ,fKMinusRemovedMult(0)
			       ,fK0sRemovedMult(0)
			       ,fK0lRemovedMult(0)
			       ,fMuMinusRemovedMult(0)
			       ,fMuPlusRemovedMult(0)
			       ,fEMinusRemovedMult(0)
			       ,fEPlusRemovedMult(0)
			       ,fGammaRemovedMult(0)
			       ,fEnergyDeposited(0)
			       ,fMomentumTPC(0)
			       ,fCharge(0)
			       ,fParticlePid(0)
			       ,fPidProb(0)
			       ,fTrackPassedCut(kFALSE)
			       ,fCentClass(0)
			       ,fDetectorRadius(0)
			       ,fSingleCellEnergyCut(0)
			       ,fChargedEnergyRemoved(0)
			       ,fNeutralEnergyRemoved(0)
			       ,fGammaEnergyAdded(0)
			       ,fHistEt(0)
			       ,fHistNeutralMult(0)
			       ,fHistPhivsPtPos(0)
			       ,fHistPhivsPtNeg(0)
			       ,fCentrality(0)
			       ,fMakeSparse(kFALSE)
			       ,fCutFlow(0)
			       ,fSelector(0)
			       ,fPIDResponse(0)
			       ,fsub(1.0)
			       ,fsubmeanhade(1.0)
	       
{}

AliAnalysisEt::~AliAnalysisEt()
{//Destructor
  delete fTmCorrections;
  delete fReCorrections;
  if(fDepositTree){
    fDepositTree->Clear();
    delete fDepositTree; // optional TTree
  }
  if(fEventSummaryTree)
  {
    fEventSummaryTree->Clear();
    delete fEventSummaryTree;
  }
  if(fAcceptedTree)
  {
    fAcceptedTree->Clear();
    delete fAcceptedTree;
  }
  delete fHistEt; //Et spectrum
  delete fHistNeutralMult; //Neutral multiplicity
  delete fHistPhivsPtPos; //phi vs pT plot for positive tracks
  delete fHistPhivsPtNeg; //phi vs pT Moplot for negative tracks
  //delete fCentrality;//this code does not actually own AliCentrality so we don't have to worry about deleting it...  we just borrow it...
  delete fCutFlow;
  delete fSelector;
  delete fPIDResponse;
}

void AliAnalysisEt::FillOutputList(TList *list)
{ // histograms to be added to output
    list->Add(fHistEt);
    list->Add(fHistNeutralMult);

    list->Add(fHistPhivsPtPos);
    list->Add(fHistPhivsPtNeg);

    if (fCuts) {
      if (fCuts->GetHistMakeTree()) {
	//list->Add(fTree);
	list->Add(fEventSummaryTree);
      }
      if (fCuts->GetHistMakeTreeDeposit()) {
	list->Add(fDepositTree);
      }
    }
    
    list->Add(fCutFlow);

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    if (!man) {
      AliFatal("Analysis manager needed");
      return;
    }

    AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler) {
      AliFatal("Input handler needed");
      return;
    }

    //pid response object
    fPIDResponse=inputHandler->GetPIDResponse();
    if (!fPIDResponse) AliError("PIDResponse object was not created");

}

void AliAnalysisEt::Init()
{// clear variables, set up cuts and PDG info
  AliAnalysisEtCommon::Init();
  if(ReadCorrections("calocorrections.root") != 0)
  {
    // Shouldn't do this, why oh why are exceptions not allowed?
    exit(-1);
  }
  ResetEventValues();

}

void AliAnalysisEt::CreateHistograms()
{ // create histograms..
  // histogram binning for E_T, p_T and Multiplicity: defaults for p+p
  Int_t nbinsEt = 10000;
  Double_t minEt = 0.0;
  Double_t maxEt = 1000;
  Int_t nbinsPt = 200;
  Double_t minPt = 0;
  Double_t maxPt = 20;
  Int_t nbinsMult = 200;
  Double_t minMult = -0.5; // offset -0.5 to have integer bins centered around 0
  Double_t maxMult = nbinsMult + minMult; // 1 bin per integer value

  // see if we should change histogram limits etc, and possibly create a tree
  if (fCuts) {
    if (fCuts->GetHistMakeTree()) {
      CreateTrees();
    }

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

    histname = "fHistNeutralMult" + fHistogramNameSuffix;
    fHistNeutralMult = new TH1F(histname.Data(), "Neutral Multiplicity", nbinsMult, minMult, maxMult);
    fHistNeutralMult->GetXaxis()->SetTitle("N");
    fHistNeutralMult->GetYaxis()->SetTitle("Multiplicity");

    histname = "fHistPhivsPtPos" + fHistogramNameSuffix;
    fHistPhivsPtPos = new TH2F(histname.Data(), "Phi vs pT of positively charged tracks hitting the calorimeter", 	200, 0, 2*TMath::Pi(), nbinsPt, minPt, maxPt);

    histname = "fHistPhivsPtNeg" + fHistogramNameSuffix;
    fHistPhivsPtNeg = new TH2F(histname.Data(), "Phi vs pT of negatively charged tracks hitting the calorimeter", 	200, 0, 2*TMath::Pi(), nbinsPt, minPt, maxPt);

    histname = "fCutFlow" + fHistogramNameSuffix;
    fCutFlow = new TH1I(histname, histname, 20, -0.5, 19.5);
}

TH2F* AliAnalysisEt::CreateEtaEHisto2D(TString name, TString title, TString ztitle)
{     //creates a 2-d histogram in eta and E and adds it to the list of histograms to be saved
	TString histoname   = name + fHistogramNameSuffix;
	
	TH2F *histo = new TH2F(histoname.Data(),title.Data(),fgNumOfEBins, fgEAxis, fgnumOfEtaBins, fgEtaAxis);
	histo->SetYTitle("#eta");
	histo->SetXTitle("E (GeV)");
	histo->SetZTitle(ztitle.Data());
	histo->Sumw2();
	
	return histo; 
}


TH2F* AliAnalysisEt::CreateEtaPtHisto2D(TString name, TString title, TString ztitle)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
	TString histoname   = name + fHistogramNameSuffix;
	
	TH2F *histo = new TH2F(histoname.Data(),title.Data(),fgNumOfPtBins, fgPtAxis, fgnumOfEtaBins, fgEtaAxis);
	histo->SetYTitle("#eta");
	histo->SetXTitle("p_{T}");
	histo->SetZTitle(ztitle.Data());
	histo->Sumw2();
	
	return histo; 
}

TH2F* AliAnalysisEt::CreateEtaEtHisto2D(TString name, TString title, TString ztitle)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
	TString histoname   = name + fHistogramNameSuffix;
	
	TH2F *histo = new TH2F(histoname.Data(),title.Data(),fgNumOfEBins, fgEAxis, fgnumOfEtaBins, fgEtaAxis);
	histo->SetYTitle("#eta");
	histo->SetXTitle("E_{T}");
	histo->SetZTitle(ztitle.Data());
	histo->Sumw2();
	
	return histo; 
}

TH2F* AliAnalysisEt::CreateResEHisto2D(TString name, TString title, TString ztitle)
{     //creates a 2-d histogram in eta and E and adds it to the list of histograms to be saved
	TString histoname   = name + fHistogramNameSuffix;
	
	TH2F *histo = new TH2F(histoname.Data(),title.Data(),fgNumOfEBins, fgEAxis, fgNumOfRBins, fgRAxis);
	histo->SetYTitle("R");
	histo->SetXTitle("E (GeV)");
	histo->SetZTitle(ztitle.Data());
	histo->Sumw2();
	
	return histo; 
}


TH2F* AliAnalysisEt::CreateResPtHisto2D(TString name, TString title, TString ztitle)
{     //creates a 2-d histogram in eta and phi and adds it to the list of histograms to be saved
	TString histoname   = name + fHistogramNameSuffix;
	
	TH2F *histo = new TH2F(histoname.Data(),title.Data(),fgNumOfPtBins, fgPtAxis, fgNumOfRBins, fgRAxis);
	histo->SetYTitle("R");
	histo->SetXTitle("p_{T}");
	histo->SetZTitle(ztitle.Data());
	histo->Sumw2();
	
	return histo; 
}

THnSparseF* AliAnalysisEt::CreateClusterHistoSparse(TString name, TString title)
{     //creates a 2D sparse histogram
	TString histoname   = name + fHistogramNameSuffix;
	
	Int_t nBins[4] = {200,200,240,20};
	Double_t min[4] = {0.,-1.,70.,0.5};
	Double_t max[4] = {50.,1.,190,20.5};
	
	THnSparseF *histo = new THnSparseF(histoname.Data(),title.Data(),4,nBins,min, max);
    histo->GetAxis(0)->SetTitle("E");
    histo->GetAxis(1)->SetTitle("#eta_{cluster}");
    histo->GetAxis(2)->SetTitle("#phi_{cluster}");
    histo->GetAxis(3)->SetTitle("n_{cells}");
	histo->Sumw2();
	
	return histo; 
}

THnSparseF* AliAnalysisEt::CreateNeutralPartHistoSparse(TString name, TString title)
{     //creates a sparse neutral particle histogram
	TString histoname   = name + fHistogramNameSuffix;
	
	Int_t nBins[6] = {20,200,200,200,240,20};
	Double_t min[6] = {-1.,0.,0.,-1.,70.,0.5};
	Double_t max[6] = {1.,50.,50.,1.,190,20.5};
	
	THnSparseF *histo = new THnSparseF(histoname.Data(),title.Data(),6,nBins,min, max);
    histo->GetAxis(0)->SetTitle("#eta");
    histo->GetAxis(1)->SetTitle("p_{T}");
    histo->GetAxis(2)->SetTitle("E");
    histo->GetAxis(3)->SetTitle("#eta_{cluster}");
    histo->GetAxis(4)->SetTitle("#phi_{cluster}");
    histo->GetAxis(5)->SetTitle("n_{cells}");
	histo->Sumw2();
	
	return histo; 
}

THnSparseF* AliAnalysisEt::CreateChargedPartHistoSparse(TString name, TString title)
{     //creates a sparse charged particle histogram
	TString histoname   = name + fHistogramNameSuffix;
	
	Int_t nBins[7] = {20,200,200,200,240,20,100};
	Double_t min[7] = {-1.,0.,0.,-1.,70.,0.5,0.};
	Double_t max[7] = {1.,50.,50.,1.,190,20.5,0.1};
	
	THnSparseF *histo = new THnSparseF(histoname.Data(),title.Data(),7,nBins,min, max);
    histo->GetAxis(0)->SetTitle("#eta");
    histo->GetAxis(1)->SetTitle("p_{T}");
    histo->GetAxis(2)->SetTitle("E");
    histo->GetAxis(3)->SetTitle("#eta_{cluster}");
    histo->GetAxis(4)->SetTitle("#phi_{cluster}");
    histo->GetAxis(5)->SetTitle("n_{cells}");
    histo->GetAxis(6)->SetTitle("R_{match}");
	histo->Sumw2();
	
	return histo; 
}


void AliAnalysisEt::CreateTrees()
{ // create tree..
  TString treename = "fEventSummaryTree" + fHistogramNameSuffix;
  if(fCuts->GetHistMakeTree())
  {
  
    fEventSummaryTree = new TTree(treename, treename);
    fEventSummaryTree->Branch("fTotEt",&fTotEt,"fTotEt/D");
    fEventSummaryTree->Branch("fTotEtAcc",&fTotEtAcc,"fTotEtAcc/D");
    fEventSummaryTree->Branch("fTotNeutralEt",&fTotNeutralEt,"fTotNeutralEt/D");
    fEventSummaryTree->Branch("fTotNeutralEtAcc",&fTotNeutralEtAcc,"fTotNeutralEtAcc/D");
    fEventSummaryTree->Branch("fTotChargedEt",&fTotChargedEt,"fTotChargedEt/D");
    fEventSummaryTree->Branch("fTotChargedEtAcc",&fTotChargedEtAcc,"fTotChargedEtAcc/D");
    fEventSummaryTree->Branch("fMultiplicity",&fMultiplicity,"fMultiplicity/I");
    fEventSummaryTree->Branch("fChargedMultiplicity",&fChargedMultiplicity,"fChargedMultiplicity/I");
    fEventSummaryTree->Branch("fNeutralMultiplicity",&fNeutralMultiplicity,"fNeutralMultiplicity/I");
    fEventSummaryTree->Branch("fCentClass",&fCentClass,"fCentClass/I");
    fEventSummaryTree->Branch("fChargedEnergyRemoved", &fChargedEnergyRemoved, "fChargedEnergyRemoved/D");
    fEventSummaryTree->Branch("fNeutralEnergyRemoved", &fNeutralEnergyRemoved, "fNeutralEnergyRemoved/D");
    fEventSummaryTree->Branch("fGammaEnergyAdded", &fGammaEnergyAdded, "fGammaEnergyAdded/D");
    

    fEventSummaryTree->Branch("fProtonEt",&fProtonEt,"fProtonEt/D");
    fEventSummaryTree->Branch("fAntiProtonEt",&fAntiProtonEt,"fAntiProtonEt/D");

    fEventSummaryTree->Branch("fNeutronEt",&fNeutronEt,"fNeutronEt/D");
    fEventSummaryTree->Branch("fAntiNeutronEt",&fAntiNeutronEt,"fAntiNeutronEt/D");

    fEventSummaryTree->Branch("fPi0Et",&fPi0Et,"fPi0Et/D");
    fEventSummaryTree->Branch("fPiPlusEt",&fPiPlusEt,"fPiPlusEt/D");
    fEventSummaryTree->Branch("fPiMinusEt",&fPiMinusEt,"fPiMinusEt/D");

    fEventSummaryTree->Branch("fKPlusEt",&fKPlusEt,"fKPlusEt/D");
    fEventSummaryTree->Branch("fKMinusEt",&fKMinusEt,"fKMinusEt/D");
    fEventSummaryTree->Branch("fK0sEt",&fK0sEt,"fK0sEt/D");
    fEventSummaryTree->Branch("fK0lEt",&fK0lEt,"fK0lEt/D");

    fEventSummaryTree->Branch("fMuMinusEt",&fMuMinusEt,"fMuMinusEt/D");
    fEventSummaryTree->Branch("fMuPlusEt",&fMuPlusEt,"fMuPlusEt/D");
    
    fEventSummaryTree->Branch("fEMinusEt",&fEMinusEt,"fEMinusEt/D");
    fEventSummaryTree->Branch("fEPlusEt",&fEPlusEt,"fEPlusEt/D");
    
    fEventSummaryTree->Branch("fGammaEt", &fGammaEt, "fGammaEt/D");
    
    fEventSummaryTree->Branch("fProtonRemovedEt",&fProtonRemovedEt,"fProtonRemovedEt/D");
    fEventSummaryTree->Branch("fAntiProtonRemovedEt",&fAntiProtonRemovedEt,"fAntiProtonRemovedEt/D");

    fEventSummaryTree->Branch("fNeutronRemovedEt",&fNeutronRemovedEt,"fNeutronRemovedEt/D");
    fEventSummaryTree->Branch("fAntiNeutronRemovedEt",&fAntiNeutronRemovedEt,"fAntiNeutronRemovedEt/D");

    fEventSummaryTree->Branch("fPi0RemovedEt",&fPi0RemovedEt,"fPi0RemovedEt/D");
    fEventSummaryTree->Branch("fPiPlusRemovedEt",&fPiPlusRemovedEt,"fPiPlusRemovedEt/D");
    fEventSummaryTree->Branch("fPiMinusRemovedEt",&fPiMinusRemovedEt,"fPiMinusRemovedEt/D");

    fEventSummaryTree->Branch("fKPlusRemovedEt",&fKPlusRemovedEt,"fKPlusRemovedEt/D");
    fEventSummaryTree->Branch("fKMinusRemovedEt",&fKMinusEt,"fKMinusRemovedEt/D");
    fEventSummaryTree->Branch("fK0sRemovedEt",&fK0sEt,"fK0sRemovedEt/D");
    fEventSummaryTree->Branch("fK0lRemovedEt",&fK0lRemovedEt,"fK0lRemovedEt/D");

    fEventSummaryTree->Branch("fMuMinusRemovedEt",&fMuMinusRemovedEt,"fMuMinusRemovedEt/D");
    fEventSummaryTree->Branch("fMuPlusRemovedEt",&fMuPlusRemovedEt,"fMuPlusRemovedEt/D");
    
    fEventSummaryTree->Branch("fEMinusRemovedEt",&fEMinusRemovedEt,"fEMinusRemovedEt/D");
    fEventSummaryTree->Branch("fEPlusRemovedEt",&fEPlusRemovedEt,"fEPlusRemovedEt/D");
    
    fEventSummaryTree->Branch("fGammaRemovedEt", &fGammaRemovedEt, "fGammaEtRemoved/D");

    fEventSummaryTree->Branch("fProtonMult",&fProtonMult,"fProtonMult/D");
    fEventSummaryTree->Branch("fAntiProtonMult",&fAntiProtonMult,"fAntiProtonMult/D");

    fEventSummaryTree->Branch("fNeutronMult",&fNeutronMult,"fNeutronMult/D");
    fEventSummaryTree->Branch("fAntiNeutronMult",&fAntiNeutronMult,"fAntiNeutronMult/D");

    fEventSummaryTree->Branch("fPi0Mult",&fPi0Mult,"fPi0Mult/D");
    fEventSummaryTree->Branch("fPiPlusMult",&fPiPlusMult,"fPiPlusMult/D");
    fEventSummaryTree->Branch("fPiMinusMult",&fPiMinusMult,"fPiMinusMult/D");

    fEventSummaryTree->Branch("fKPlusMult",&fKPlusMult,"fKPlusMult/D");
    fEventSummaryTree->Branch("fKMinusMult",&fKMinusMult,"fKMinusMult/D");
    fEventSummaryTree->Branch("fK0sMult",&fK0sMult,"fK0sMult/D");
    fEventSummaryTree->Branch("fK0lMult",&fK0lMult,"fK0lMult/D");

    fEventSummaryTree->Branch("fMuMinusMult",&fMuMinusMult,"fMuMinusMult/D");
    fEventSummaryTree->Branch("fMuPlusMult",&fMuPlusMult,"fMuPlusMult/D");
    
    fEventSummaryTree->Branch("fEMinusMult",&fEMinusMult,"fEMinusMult/D");
    fEventSummaryTree->Branch("fEPlusMult",&fEPlusMult,"fEPlusMult/D");
    
    fEventSummaryTree->Branch("fGammaMult", &fGammaMult, "fGammaMult/D");
    
    fEventSummaryTree->Branch("fProtonRemovedMult",&fProtonRemovedMult,"fProtonRemovedMult/D");
    fEventSummaryTree->Branch("fAntiProtonRemovedMult",&fAntiProtonRemovedMult,"fAntiProtonRemovedMult/D");

    fEventSummaryTree->Branch("fNeutronRemovedMult",&fNeutronRemovedMult,"fNeutronRemovedMult/D");
    fEventSummaryTree->Branch("fAntiNeutronRemovedMult",&fAntiNeutronRemovedMult,"fAntiNeutronRemovedMult/D");

    fEventSummaryTree->Branch("fPi0RemovedMult",&fPi0RemovedMult,"fPi0RemovedMult/D");
    fEventSummaryTree->Branch("fPiPlusRemovedMult",&fPiPlusRemovedMult,"fPiPlusRemovedMult/D");
    fEventSummaryTree->Branch("fPiMinusRemovedMult",&fPiMinusRemovedMult,"fPiMinusRemovedMult/D");

    fEventSummaryTree->Branch("fKPlusRemovedMult",&fKPlusRemovedMult,"fKPlusRemovedMult/D");
    fEventSummaryTree->Branch("fKMinusRemovedMult",&fKMinusMult,"fKMinusRemovedMult/D");
    fEventSummaryTree->Branch("fK0sRemovedMult",&fK0sMult,"fK0sRemovedMult/D");
    fEventSummaryTree->Branch("fK0lRemovedMult",&fK0lRemovedMult,"fK0lRemovedMult/D");

    fEventSummaryTree->Branch("fMuMinusRemovedMult",&fMuMinusRemovedMult,"fMuMinusRemovedMult/D");
    fEventSummaryTree->Branch("fMuPlusRemovedMult",&fMuPlusRemovedMult,"fMuPlusRemovedMult/D");
    
    fEventSummaryTree->Branch("fEMinusRemovedMult",&fEMinusRemovedMult,"fEMinusRemovedMult/D");
    fEventSummaryTree->Branch("fEPlusRemovedMult",&fEPlusRemovedMult,"fEPlusRemovedMult/D");
    
    fEventSummaryTree->Branch("fGammaRemovedMult", &fGammaRemovedMult, "fGammaMultRemoved/D");

    
    
  }
  
  if(fCuts->GetHistMakeTreeDeposit())
  {
    treename = "fTreeDeposit" + fHistogramNameSuffix;
    fDepositTree = new TTree(treename, treename);
  
    fDepositTree->Branch("fEnergyDeposited", &fEnergyDeposited, "fEnergyDeposited/F");
    fDepositTree->Branch("fMomentumTPC", &fMomentumTPC, "fMomentumTPC/F");
    fDepositTree->Branch("fCharge", &fCharge, "fCharge/S");
    fDepositTree->Branch("fParticlePid", &fParticlePid, "fParticlePid/S");
    fDepositTree->Branch("fPidProb", &fPidProb, "fPidProb/F");
    fDepositTree->Branch("fTrackPassedCut", &fTrackPassedCut, "fTrackPassedCut/B");
  }

  return;
}
void AliAnalysisEt::FillHistograms()
{ // fill histograms..
    fHistEt->Fill(fTotEt);

    fHistNeutralMult->Fill(fNeutralMultiplicity);

    if (fCuts) {
      if (fCuts->GetHistMakeTree()) {
	fEventSummaryTree->Fill();
      }
    }
    if(fCuts)
    {
      if(fCuts->GetHistMakeTreeDeposit())
      {
	fDepositTree->Fill();
      }
    }
      
}

Int_t AliAnalysisEt::AnalyseEvent(AliVEvent *event)
{ //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
  AliAnalysisEtCommon::AnalyseEvent(event);
  //fSelector->SetEvent(event);
  ResetEventValues();
  return 0;
}

void AliAnalysisEt::ResetEventValues()
{ // clear
  AliAnalysisEtCommon::ResetEventValues();
  fTotEt = 0;
  fTotEtAcc = 0;
  fTotNeutralEt = 0;
  fTotChargedEt = 0;
  fTotNeutralEtAcc = 0;
  fTotChargedEtAcc  = 0;
  fMultiplicity = 0;
  fChargedMultiplicity = 0;
  fNeutralMultiplicity = 0;
  
  fProtonEt = 0;
  fAntiProtonEt = 0;
  
  fNeutronEt = 0;
  fAntiNeutronEt = 0;
  
  fPi0Et = 0;
  fPiPlusEt = 0;
  fPiMinusEt = 0;
  
  fKPlusEt = 0;
  fKMinusEt = 0;
  fK0sEt = 0;
  fK0lEt = 0;
  
  fMuMinusEt = 0;
  fMuPlusEt = 0;
  
  fEMinusEt = 0;
  fEPlusEt = 0;
  
  fGammaEt = 0;
  
  fProtonRemovedEt = 0;
  fAntiProtonRemovedEt = 0;
  
  fNeutronRemovedEt = 0;
  fAntiNeutronRemovedEt = 0;
  
  fPi0RemovedEt = 0;
  fPiPlusRemovedEt = 0;
  fPiMinusRemovedEt = 0;
  
  fKPlusRemovedEt = 0;
  fKMinusRemovedEt = 0;
  fK0sRemovedEt = 0;
  fK0lRemovedEt = 0;
  
  fMuMinusRemovedEt = 0;
  fMuPlusRemovedEt = 0;
  
  fEMinusRemovedEt = 0;
  fEPlusRemovedEt = 0;
  
  fGammaRemovedEt = 0;
  
  
  fProtonMult = 0;
  fAntiProtonMult = 0;
  
  fNeutronMult = 0;
  fAntiNeutronMult = 0;
  
  fPi0Mult = 0;
  fPiPlusMult = 0;
  fPiMinusMult = 0;
  
  fKPlusMult = 0;
  fKMinusMult = 0;
  fK0sMult = 0;
  fK0lMult = 0;
  
  fMuMinusMult = 0;
  fMuPlusMult = 0;
  
  fEMinusMult = 0;
  fEPlusMult = 0;
  
  fGammaMult = 0;
  
  fProtonRemovedMult = 0;
  fAntiProtonRemovedMult = 0;
  
  fNeutronRemovedMult = 0;
  fAntiNeutronRemovedMult = 0;
  
  fPi0RemovedMult = 0;
  fPiPlusRemovedMult = 0;
  fPiMinusRemovedMult = 0;
  
  fKPlusRemovedMult = 0;
  fKMinusRemovedMult = 0;
  fK0sRemovedMult = 0;
  fK0lRemovedMult = 0;
  
  fMuMinusRemovedMult = 0;
  fMuPlusRemovedMult = 0;
  
  fEMinusRemovedMult = 0;
  fEPlusRemovedMult = 0;
  
  fGammaRemovedMult = 0;
  
  return;
}

Int_t AliAnalysisEt::ReadCorrections(TString filename)
{
  TFile *f = TFile::Open(filename, "READ");
  if(f)
  {
    TString det = "Phos";
    if(fHistogramNameSuffix.Contains("Emcal"))
    {
      det = "Emcal";
    }
    cout<<"Histo name suffix "<<fHistogramNameSuffix<<endl;
    TString name = "TmCorrections" + det;
    std::cout << name << std::endl;
    fTmCorrections = dynamic_cast<AliAnalysisEtTrackMatchCorrections*>(f->Get(name));
    if(!fTmCorrections)
    {
      cout<<"No corrections with name "<<name<<endl;
      Printf("Could not load TM corrections");
      return -1;
    }
    name = "ReCorrections" + det;
    fReCorrections = dynamic_cast<AliAnalysisEtRecEffCorrection*>(f->Get(name));
    if(!fReCorrections)
    {
      Printf("Could not load rec eff corrections");
      return -1;
    }
    return 0;
  }
 return -1; 
}

Double_t AliAnalysisEt::CorrectForReconstructionEfficiency(const AliESDCaloCluster& cluster, Int_t cent)
{
  Float_t pos[3];
  cluster.GetPosition(pos);
  TVector3 cp(pos);
  Double_t corrEnergy = fReCorrections->CorrectedEnergy(cluster.E(), cent);
  
  //std::cout << "Original energy: " << cluster.E() << ", corrected energy: " << corrEnergy << std::endl;
  return TMath::Sin(cp.Theta())*corrEnergy;
}


Double_t AliAnalysisEt::CorrectForReconstructionEfficiency(const AliESDCaloCluster& cluster, Float_t eReco, Int_t cent)
{
  Float_t pos[3];
  cluster.GetPosition(pos);
  TVector3 cp(pos);
  Double_t corrEnergy = fReCorrections->CorrectedEnergy(eReco, cent);
  
  //std::cout << "Original energy: " << cluster.E() << ", corrected energy: " << corrEnergy << std::endl;
  return TMath::Sin(cp.Theta())*corrEnergy;
}

