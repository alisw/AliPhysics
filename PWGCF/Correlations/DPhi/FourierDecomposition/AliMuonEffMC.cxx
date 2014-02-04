// MUON track QA referring AliMuonEffMC.cxx
// Author : Saehanseul Oh

#include "AliMuonEffMC.h"

#include <TList.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>

#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliAODMCParticle.h"

using std::cout;
using std::endl;

ClassImp(AliMuonEffMC)

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC() :
  AliAnalysisTaskSE(), fESD(0), fAOD(0), fMC(0), fStack(0), fCentrality(99), fZVertex(99), 
  fOutputList(0x0),fHEventStat(0), fHXsec(0), fHTrials(0), fHEvt(0x0), fIsMc(kTRUE), fIsPythia(kFALSE), fMDProcess(kFALSE), fPlotMode(0),
  fCentralityEstimator("V0M"), fNEtaBins(100), fNpTBins(50), fNCentBins(1), fNZvtxBins(1), fNPhiBins(12), 
  fHFPM(0x0), fHPP(0)
{
  // Constructor
  //DefineInput(0, TChain::Class());
  //DefineOutput(1, TList::Class());
  for(Int_t i=0; i<2; i++)
  {
    fHDetRecMu[i] = NULL;
    fHDetRecMuFPM[i] = NULL;
    fHDetRecMuPP[i] = NULL;
    fHMuFPM[i] = NULL;
    fHMuPP[i] = NULL;
  }
  for(Int_t i=0; i<2; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      fHMuMotherRecPt[i][j] = NULL;
      fHMuMotherRecPhi[i][j] = NULL;
      fHMuMotherRecEta[i][j] = NULL;
       for(Int_t k=0; k<3; k++)
      {
	fHMuMohterPtDifRec[i][j][k] = NULL;
	fHMuMohterPhiDifRec[i][j][k] = NULL;
	fHMuMohterEtaDifRec[i][j][k] = NULL;
      }
    }
  }
  for(Int_t i=0; i<3; i++)
  {
    fHZvRv[i] = NULL;
    fHXvYv[i] = NULL;
  }
}

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC(const char *name) :
  AliAnalysisTaskSE(name), fESD(0), fAOD(0), fMC(0), fStack(0), fCentrality(99), fZVertex(99), 
  fOutputList(0x0),fHEventStat(0), fHXsec(0), fHTrials(0), fHEvt(0x0), fIsMc(kTRUE), fIsPythia(kFALSE), fMDProcess(kFALSE), fPlotMode(0),
  fCentralityEstimator("V0M"), fNEtaBins(100), fNpTBins(50), fNCentBins(1), fNZvtxBins(1), fNPhiBins(12), 
  fHFPM(0x0), fHPP(0)
{
  // Constructor
  for(Int_t i=0; i<2; i++)
  {
    fHDetRecMu[i] = NULL;
    fHDetRecMuFPM[i] = NULL;
    fHDetRecMuPP[i] = NULL;
    fHMuFPM[i] = NULL;
    fHMuPP[i] = NULL;
  }
  for(Int_t i=0; i<2; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      fHMuMotherRecPt[i][j] = NULL;
      fHMuMotherRecPhi[i][j] = NULL;
      fHMuMotherRecEta[i][j] = NULL;
      for(Int_t k=0; k<3; k++)
      {
	fHMuMohterPtDifRec[i][j][k] = NULL;
	fHMuMohterPhiDifRec[i][j][k] = NULL;
	fHMuMohterEtaDifRec[i][j][k] = NULL;
      }
    }
  }
  for(Int_t i=0; i<3; i++)
  {
    fHZvRv[i] = NULL;
    fHXvYv[i] = NULL;
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliMuonEffMC::~AliMuonEffMC()
{
  //Destructor
  if(fOutputList) delete fOutputList;
}

//________________________________________________________________________
void AliMuonEffMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (per slave on PROOF!)
  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fHEventStat = new TH1D("fHEventStat","Event statistics for analysis",18,0,18);
  fHEventStat->GetXaxis()->SetBinLabel(1,"Event");
  fHEventStat->GetXaxis()->SetBinLabel(2,"SelectedEvent");
  fHEventStat->GetXaxis()->SetBinLabel(3,"File");
  fHEventStat->GetXaxis()->SetBinLabel(4,"fSPHighpt");  //!Global Trigger Single plus High p_T
  fHEventStat->GetXaxis()->SetBinLabel(5,"fSPAllpt");   //!Global Trigger Single plus All p_T
  fHEventStat->GetXaxis()->SetBinLabel(6,"fSMLowpt");   //!Global Trigger Single minus Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(7,"fSMHighpt");  //!Global Trigger Single minus High p_T
  fHEventStat->GetXaxis()->SetBinLabel(8,"fSMAllpt");   //!Global Trigger Single minus All p_T
  fHEventStat->GetXaxis()->SetBinLabel(9,"fSULowpt");   //!Global Trigger Single undefined Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(10,"fSUHighpt"); //!Global Trigger Single undefined High p_T
  fHEventStat->GetXaxis()->SetBinLabel(11,"fSUAllpt");  //!Global Trigger Single undefined All p_T
  fHEventStat->GetXaxis()->SetBinLabel(12,"fUSLowpt");  //!Global Trigger UnlikeSign pair Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(13,"fUSHighpt"); //!Global Trigger UnlikeSign pair High p_T
  fHEventStat->GetXaxis()->SetBinLabel(14,"fUSAllpt");  //!Global Trigger UnlikeSign pair All p_T
  fHEventStat->GetXaxis()->SetBinLabel(15,"fLSLowpt");  //!Global Trigger LikeSign pair pair Low p_T
  fHEventStat->GetXaxis()->SetBinLabel(16,"fLSHighpt"); //!Global Trigger LikeSign pair pair High p_T
  fHEventStat->GetXaxis()->SetBinLabel(17,"fLSAllpt");  //!Global Trigger LikeSign pair pair All p_T
  fHEventStat->GetXaxis()->SetBinLabel(18,"fSPLowpt");   //!Global Trigger Single plus Low p_T
  fOutputList->Add(fHEventStat);

  if(fIsPythia)
  {
    fHXsec = new TH1F("fHXsec", "Cross section from pyxsec.root", 1, 0, 1);
    fHXsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
    fOutputList->Add(fHXsec);
    
    fHTrials = new TH1F("fHTrials", "Number of Trials", 1, 0, 1);
    fHTrials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    fOutputList->Add(fHTrials);  
  }

  fHEvt = new TH2F("fHEvt", "Event-level variables; Zvtx; Cent", 30, -15, 15, 103, -2, 101);
  fOutputList->Add(fHEvt);

  // Define THn's
  Int_t iTrackBinFPM[7];
  Double_t* trackBinsFPM[7];
  const char* trackAxisTitleFPM[7];

  Int_t iTrackBinPP[7];
  Double_t* trackBinsPP[7];
  const char* trackAxisTitlePP[7];

  Int_t iTrackBinMu[7];
  Double_t* trackBinsMu[7];
  const char* trackAxisTitleMu[7];

  Int_t iTrackBinMuFPM[6];
  Double_t* trackBinsMuFPM[6];
  const char* trackAxisTitleMuFPM[6];

  Int_t iTrackBinMuPP[6];
  Double_t* trackBinsMuPP[6];
  const char* trackAxisTitleMuPP[6];

  // eta
  Double_t etaBins[fNEtaBins+1];
  for(Int_t i=0; i<=fNEtaBins; i++) { etaBins[i] = (Double_t)(-5.0 + 10.0/fNEtaBins*i); }
  iTrackBinFPM[0] = fNEtaBins;
  trackBinsFPM[0] = etaBins;
  trackAxisTitleFPM[0] = "#eta";

  iTrackBinPP[0] = fNEtaBins;
  trackBinsPP[0] = etaBins;
  trackAxisTitlePP[0] = "#eta";

  iTrackBinMu[0] = fNEtaBins;
  trackBinsMu[0] = etaBins;
  trackAxisTitleMu[0] = "#eta";

  iTrackBinMuFPM[0] = fNEtaBins;
  trackBinsMuFPM[0] = etaBins;
  trackAxisTitleMuFPM[0] = "#eta";

  iTrackBinMuPP[0] = fNEtaBins;
  trackBinsMuPP[0] = etaBins;
  trackAxisTitleMuPP[0] = "#eta";

  // p_T
  Double_t pTBins[fNpTBins+1];
  for(Int_t i=0; i<=fNpTBins; i++) { pTBins[i] = (Double_t)(5.0/fNpTBins * i); }
  iTrackBinFPM[1] = fNpTBins;
  trackBinsFPM[1] = pTBins;
  trackAxisTitleFPM[1] = "p_{T} (GeV/c)";

  iTrackBinPP[1] = fNpTBins;
  trackBinsPP[1] = pTBins;
  trackAxisTitlePP[1] = "p_{T} (GeV/c)";

  iTrackBinMu[1] = fNpTBins;
  trackBinsMu[1] = pTBins;
  trackAxisTitleMu[1] = "p_{T} (GeV/c)";

  iTrackBinMuFPM[1] = fNpTBins;
  trackBinsMuFPM[1] = pTBins;
  trackAxisTitleMuFPM[1] = "p_{T} (GeV/c)";

  iTrackBinMuPP[1] = fNpTBins;
  trackBinsMuPP[1] = pTBins;
  trackAxisTitleMuPP[1] = "p_{T} (GeV/c)";

  // centrality
  Double_t CentBins[fNCentBins+1];
  for (Int_t i=0; i<=fNCentBins; i++) { CentBins[i] = (Double_t)(100.0/fNCentBins * i); }
  iTrackBinFPM[2] = fNCentBins;
  trackBinsFPM[2] = CentBins;
  trackAxisTitleFPM[2] = "Cent";

  iTrackBinPP[2] = fNCentBins;
  trackBinsPP[2] = CentBins;
  trackAxisTitlePP[2] = "Cent";

  iTrackBinMu[2] = fNCentBins;
  trackBinsMu[2] = CentBins;
  trackAxisTitleMu[2] = "Cent";

  // Z-vertex
  Double_t ZvtxBins[fNZvtxBins+1];
  for(Int_t i=0; i<=fNZvtxBins; i++) { ZvtxBins[i] = (Double_t)(-10.0 + 20.0/fNZvtxBins * i); }
  iTrackBinFPM[3] = fNZvtxBins;
  trackBinsFPM[3] = ZvtxBins;
  trackAxisTitleFPM[3] = "Zvtx";

  iTrackBinPP[3] = fNZvtxBins;
  trackBinsPP[3] = ZvtxBins;
  trackAxisTitlePP[3] = "Zvtx";

  iTrackBinMu[3] = fNZvtxBins;
  trackBinsMu[3] = ZvtxBins;
  trackAxisTitleMu[3] = "Zvtx";

  // phi
  Double_t phiBins[fNPhiBins+1];
  for(Int_t i=0; i<=fNPhiBins; i++) { phiBins[i] = (Double_t)(TMath::TwoPi()/fNPhiBins * i); }
  iTrackBinFPM[4] = fNPhiBins;
  trackBinsFPM[4] = phiBins;
  trackAxisTitleFPM[4] = "#phi";

  iTrackBinPP[4] = fNPhiBins;
  trackBinsPP[4] = phiBins;
  trackAxisTitlePP[4] = "#phi";

  iTrackBinMu[4] = fNPhiBins;
  trackBinsMu[4] = phiBins;
  trackAxisTitleMu[4] = "#phi";

  iTrackBinMuFPM[2] = fNPhiBins;
  trackBinsMuFPM[2] = phiBins;
  trackAxisTitleMuFPM[2] = "#phi";

  iTrackBinMuPP[2] = fNPhiBins;
  trackBinsMuPP[2] = phiBins;
  trackAxisTitleMuPP[2] = "#phi";

  // charge
  Double_t chargeBins[4] = {-10.0, -0.5, 0.5, 10.0};
  iTrackBinFPM[5] = 3;
  trackBinsFPM[5] = chargeBins;
  trackAxisTitleFPM[5] = "charge";

  iTrackBinPP[5] = 3;
  trackBinsPP[5] = chargeBins;
  trackAxisTitlePP[5] = "charge";

  iTrackBinMu[5] = 3;
  trackBinsMu[5] = chargeBins;
  trackAxisTitleMu[5] = "charge";

  iTrackBinMuFPM[3] = 3;
  trackBinsMuFPM[3] = chargeBins;
  trackAxisTitleMuFPM[3] = "charge";

  iTrackBinMuPP[3] = 3;
  trackBinsMuPP[3] = chargeBins;
  trackAxisTitleMuPP[3] = "charge";

  // Muon type
  Double_t MuSpeciesBins[4] = {0.0, 1.0, 2.0, 3.0};
  iTrackBinMu[6] = 3;
  trackBinsMu[6] = MuSpeciesBins;
  trackAxisTitleMu[6] = "MUON type";

  iTrackBinMuFPM[4] = 3;
  trackBinsMuFPM[4] = MuSpeciesBins;
  trackAxisTitleMuFPM[4] = "MUON type";

  iTrackBinMuPP[4] = 3;
  trackBinsMuPP[4] = MuSpeciesBins;
  trackAxisTitleMuPP[4] = "MUON type";

  // FPM species
  Double_t FPMSpecies[8] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  iTrackBinFPM[6] = 7;
  trackBinsFPM[6] = FPMSpecies;
  trackAxisTitleFPM[6] = "FPM species";

  iTrackBinMuFPM[5] = 7;
  trackBinsMuFPM[5] = FPMSpecies;
  trackAxisTitleMuFPM[5] = "FPM species";

  // First Physical Primary  species
  Double_t PPMSpecies[8] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  iTrackBinPP[6] = 7;
  trackBinsPP[6] = PPMSpecies;
  trackAxisTitlePP[6] = "PPM species";

  iTrackBinMuPP[5] = 7;
  trackBinsMuPP[5] = PPMSpecies;
  trackAxisTitleMuPP[5] = "PPM species";
 
  const char* MuonType[3] = {"Prim","Sec","Had"};
  const char *MuPt[3] = {"0005","0520","2040"};
  const char *cutlabel[2] = {"cut1", "cut2"};

  if(fIsMc)
  {
    // THn for tracking efficiency
    if(fPlotMode==0)
    {
      fHFPM = new THnF("fHFPM", "", 7, iTrackBinFPM, 0, 0);
      for (Int_t i=0; i<7; i++)
      {
	fHFPM->SetBinEdges(i, trackBinsFPM[i]);
	fHFPM->GetAxis(i)->SetTitle(trackAxisTitleFPM[i]);
      }
      fOutputList->Add(fHFPM); 
      
      fHPP = new THnF("fHPP", "", 7, iTrackBinPP, 0, 0);
      for (Int_t i=0; i<7; i++)
      {
	fHPP->SetBinEdges(i, trackBinsPP[i]);
	fHPP->GetAxis(i)->SetTitle(trackAxisTitlePP[i]);
      }
      fOutputList->Add(fHPP); 
    }

    for(Int_t j=0; j<2; j++)
    {
      if(fPlotMode==1)
      {
	fHDetRecMu[j] = new THnF(Form("fHDetRecMu_%s",cutlabel[j]),"", 7, iTrackBinMu, 0, 0);
	for (Int_t i=0; i<7; i++)
	{
	  fHDetRecMu[j]->SetBinEdges(i, trackBinsMu[i]);
	  fHDetRecMu[j]->GetAxis(i)->SetTitle(trackAxisTitleMu[i]);
	}
	fOutputList->Add(fHDetRecMu[j]); 
      }
      if(fPlotMode==2)
      {	
	fHDetRecMuFPM[j] = new THnF(Form("fHDetRecMuFPM_%s",cutlabel[j]),"", 6, iTrackBinMuFPM, 0, 0);
	for (Int_t i=0; i<6; i++)
	{
	  fHDetRecMuFPM[j]->SetBinEdges(i, trackBinsMuFPM[i]);
	  fHDetRecMuFPM[j]->GetAxis(i)->SetTitle(trackAxisTitleMuFPM[i]);
	}
	fOutputList->Add(fHDetRecMuFPM[j]); 

	fHDetRecMuPP[j] = new THnF(Form("fHDetRecMuPP_%s",cutlabel[j]),"", 6, iTrackBinMuPP, 0, 0);
	for (Int_t i=0; i<6; i++)
	{
	  fHDetRecMuPP[j]->SetBinEdges(i, trackBinsMuPP[i]);
	  fHDetRecMuPP[j]->GetAxis(i)->SetTitle(trackAxisTitleMuPP[i]);
	}
	fOutputList->Add(fHDetRecMuPP[j]); 
      }
      if(fPlotMode==3)
      {
	fHMuFPM[j] = new THnF(Form("fHMuFPM_%s",cutlabel[j]),"", 6, iTrackBinMuFPM, 0, 0);
	for (Int_t i=0; i<6; i++)
	{
	  fHMuFPM[j]->SetBinEdges(i, trackBinsMuFPM[i]);
	  fHMuFPM[j]->GetAxis(i)->SetTitle(trackAxisTitleMuFPM[i]);
	}
	fOutputList->Add(fHMuFPM[j]); 

	fHMuPP[j] = new THnF(Form("fHMuPP_%s",cutlabel[j]),"", 6, iTrackBinMuPP, 0, 0);
	for (Int_t i=0; i<6; i++)
	{
	  fHMuPP[j]->SetBinEdges(i, trackBinsMuPP[i]);
	  fHMuPP[j]->GetAxis(i)->SetTitle(trackAxisTitleMuPP[i]);
	}
	fOutputList->Add(fHMuPP[j]); 
      }
    }

    if(fMDProcess)
    {
      for(Int_t i=0; i<2; i++)
      {
	for(Int_t j=0; j<3; j++)
	{
	  fHMuMotherRecPt[i][j] = new TH2F(Form("fHMuMotherRecPt_%s_%s",cutlabel[i], MuonType[j]),";p_{T,muon}^{rec} (GeV/c);p_{T,mother}^{Truth} (GeV/c);",500, 0, 50, 500, 0, 50);
	  fOutputList->Add(fHMuMotherRecPt[i][j]);
	  fHMuMotherRecPhi[i][j] = new TH2F(Form("fHMuMotherRecPhi_%s_%s",cutlabel[i], MuonType[j]),";#phi_{rec};mother #phi;",100, 0, TMath::TwoPi(), 100, 0, TMath::TwoPi());
	  fOutputList->Add(fHMuMotherRecPhi[i][j]);
	  fHMuMotherRecEta[i][j] = new TH2F(Form("fHMuMotherRecEta_%s_%s",cutlabel[i], MuonType[j]),";#eta_{rec};mother #eta;",100, -5., -1., 100, -5., -1.);
	  fOutputList->Add(fHMuMotherRecEta[i][j]);

	  for(Int_t k=0; k<3; k++)
	  {
	    fHMuMohterPtDifRec[i][j][k] = new TH1F(Form("fHMuMohterPtDifRec_%s_%s_%s",cutlabel[i], MuonType[j], MuPt[k]),";#Delta#phi",200, -10.0, 10.0);
	    fOutputList->Add(fHMuMohterPtDifRec[i][j][k]);

	    fHMuMohterPhiDifRec[i][j][k] = new TH1F(Form("fHMuMohterPhiDifRec_%s_%s_%s",cutlabel[i], MuonType[j], MuPt[k]),";#Delta#phi",100, -1.0*TMath::Pi(), TMath::Pi());
	    fOutputList->Add(fHMuMohterPhiDifRec[i][j][k]);

	    fHMuMohterEtaDifRec[i][j][k] = new TH1F(Form("fHMuMohterEtaDifRec_%s_%s_%s",cutlabel[i], MuonType[j], MuPt[k]),";#Delta#eta",100, -5.0, 5.0);
	    fOutputList->Add(fHMuMohterEtaDifRec[i][j][k]);
	  }
	}
      }
      for(Int_t i = 0; i<3; i++)
      {
	fHZvRv[i] = new TH2F(Form("fHZvRv_%s",MuonType[i]), "", 300, -500, 100, 200, 0, 800);
	fOutputList->Add(fHZvRv[i]);
      	fHXvYv[i] = new TH2F(Form("fHXvYv_%s",MuonType[i]), "", 200, -500, 500, 200, -500, 500);
	fOutputList->Add(fHXvYv[i]);
      }	
    }
  }
  else
  {
    for(Int_t j=0; j<2; j++)
    {
      fHDetRecMu[j] = new THnF(Form("fHDetRecMu_%s",cutlabel[j]),"", 7, iTrackBinMu, 0, 0);
      for (Int_t i=0; i<7; i++)
      {
	fHDetRecMu[j]->SetBinEdges(i, trackBinsMu[i]);
	fHDetRecMu[j]->GetAxis(i)->SetTitle(trackAxisTitleMu[i]);
      }
      fOutputList->Add(fHDetRecMu[j]); 
    }
  }
  PostData(1, fOutputList);
}

//________________________________________________________________________
Bool_t AliMuonEffMC::Notify()
{
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  if(fIsPythia)
  {
    fHEventStat->Fill(2.5);
    
    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    Float_t xsection = 0;
    Float_t ftrials  = 1;
    
    if(tree){
      TFile *curfile = tree->GetCurrentFile();
      if (!curfile) {
	Error("Notify","No current file");
	return kFALSE;
      }

      AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
      fHXsec->Fill("<#sigma>",xsection);
      fHTrials->Fill("#sum{ntrials}",ftrials);
    }  
  }
  return kTRUE;
}
//________________________________________________________________________
void AliMuonEffMC::UserExec(Option_t *)
{
  // Main loop, Called for each event
  Int_t ntrks = 0; // number of tracks in an event

  if(((TString)InputEvent()->IsA()->GetName())=="AliAODEvent")
  {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) { AliError("AOD event not found. Nothing done!"); return; }
    ntrks = fAOD->GetNTracks();
  }
  else
  {
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) { AliError("ESD event not found. Nothing done!"); return; }
    ntrks = fESD->GetNumberOfMuonTracks();
  }
       
  fHEventStat->Fill(0.5);
  if(fIsMc)
  {
    fMC = MCEvent();
    if (!fMC) { AliError("MC event not avaliable."); return; }
    fStack = fMC->Stack();
  }

  // Centrality, vertex, other event variables...
  if(fAOD)
  {
    const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
    fZVertex = vertex->GetZ();
    if(fAOD->GetCentrality())  fCentrality = fAOD->GetCentrality()->GetCentralityPercentile(fCentralityEstimator);
  }
  else if(fESD)
  {
    const AliESDVertex* vertex = fESD->GetPrimaryVertex();
    fZVertex = vertex->GetZ();
    if(fESD->GetCentrality()) fCentrality = fESD->GetCentrality()->GetCentralityPercentile(fCentralityEstimator);
  }  

  if ((fESD && !VertexOk(fESD)) || (fAOD && !VertexOk(fAOD))) { 
    //AliInfo(Form("Event REJECTED. z = %.1f", fZVertex));
    return; 
  }
  if (fCentrality > 100. || fCentrality < -1.5) { 
    //AliInfo(Form("Event REJECTED. fCentrality = %.1f", fCentrality));
    return; 
  }
 
  if(fCentrality < 0) fCentrality = 0.5; //ad hoc centrality for pp
  // Fill Event histogram
  fHEvt->Fill(fZVertex, fCentrality);
  fHEventStat->Fill(1.5);
 
  ULong64_t trigword = 0;
  if(fAOD) trigword=fAOD->GetTriggerMask();
  else if(fESD) trigword=fESD->GetTriggerMask();
 
  if (trigword & 0x01) fHEventStat->Fill(17.5);
  if (trigword & 0x02) fHEventStat->Fill(3.5);
  if (trigword & 0x04) fHEventStat->Fill(4.5);
  if (trigword & 0x08) fHEventStat->Fill(5.5);      
  if (trigword & 0x010) fHEventStat->Fill(6.5);
  if (trigword & 0x020) fHEventStat->Fill(7.5);
  if (trigword & 0x040) fHEventStat->Fill(8.5);
  if (trigword & 0x080) fHEventStat->Fill(9.5);
  if (trigword & 0x100) fHEventStat->Fill(10.5);
  if (trigword & 0x200) fHEventStat->Fill(11.5);
  if (trigword & 0x400) fHEventStat->Fill(12.5);
  if (trigword & 0x800) fHEventStat->Fill(13.5);
  if (trigword & 0x1000) fHEventStat->Fill(14.5);
  if (trigword & 0x2000) fHEventStat->Fill(15.5);
  if (trigword & 0x4000) fHEventStat->Fill(16.5);

  if(fIsMc && fPlotMode==0)
  {
    Double_t PPEta = 0.0;
    Double_t PPPt = 0.0;
    Double_t PPPhi = 0.0;
    Double_t PPCharge = 0.0;
    Double_t PPSpecies = 0.0;
    
    Double_t FPMEta = 0.0;
    Double_t FPMPt = 0.0;
    Double_t FPMPhi = 0.0;
    Double_t FPMCharge = 0.0;
    Double_t FPMSpecies = 0.0;

    // generated level loop
    for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++)
    {
      if(fESD)
      {
	AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(ipart);
	if(!fMC->IsPhysicalPrimary(ipart) || McParticle->Charge()==0) continue;

	PPEta = McParticle->Eta();
	PPPt = McParticle->Pt();
	PPPhi = McParticle->Phi();
	PPCharge = McParticle->Charge();
	PPSpecies = GetSpecies(McParticle->PdgCode());

	AliMCParticle *FPMParticle  = (AliMCParticle*)fMC->GetTrack(GetFirstPrimaryMother(ipart));
	if(PPSpecies==1.5 || PPSpecies==2.5 || PPSpecies==4.5)
	{
	  if(ipart < fStack->GetNprimary())
	  {
	    FPMEta = PPEta;
	    FPMPt = PPPt;
	    FPMPhi = PPPhi;
	    FPMCharge = PPCharge;
	    FPMSpecies = PPSpecies;
	  }
	  else
	  {
	    FPMEta = FPMParticle->Eta();
	    FPMPt = FPMParticle->Pt();
	    FPMPhi = FPMParticle->Phi();
	    FPMCharge = FPMParticle->Charge();
	    FPMSpecies = GetSpecies(FPMParticle->PdgCode());
	  }
	}
	else
	{
	  FPMEta = FPMParticle->Eta();
	  FPMPt = FPMParticle->Pt();
	  FPMPhi = FPMParticle->Phi();
	  FPMCharge = FPMParticle->Charge();
	  FPMSpecies = GetSpecies(FPMParticle->PdgCode());
	}
      }
      if(fAOD)
      {
	AliAODMCParticle *AodMcParticle  = (AliAODMCParticle*)fMC->GetTrack(ipart);
	if(!fMC->IsPhysicalPrimary(ipart) || AodMcParticle->Charge()==0) continue;

	PPEta = AodMcParticle->Eta();
	PPPt = AodMcParticle->Pt();
	PPPhi = AodMcParticle->Phi();
	PPCharge = AodMcParticle->Charge();
	PPSpecies = GetSpecies(AodMcParticle->PdgCode());

	AliAODMCParticle *AODFPMParticle  = (AliAODMCParticle*)fMC->GetTrack(GetFirstPrimaryMother(ipart));
	if(PPSpecies==1.5 || PPSpecies==2.5)
	{
	  if(ipart < fStack->GetNprimary())
	  {
	    FPMEta = PPEta;
	    FPMPt = PPPt;
	    FPMPhi = PPPhi;
	    FPMCharge = PPCharge;
	    FPMSpecies = PPSpecies;
	  }
	  else
	  {
	    FPMEta = AODFPMParticle->Eta();
	    FPMPt = AODFPMParticle->Pt();
	    FPMPhi = AODFPMParticle->Phi();
	    FPMCharge = AODFPMParticle->Charge();
	    FPMSpecies = GetSpecies(AODFPMParticle->PdgCode());
	  }
	}
	else
	{
	  FPMEta = AODFPMParticle->Eta();
	  FPMPt = AODFPMParticle->Pt();
	  FPMPhi = AODFPMParticle->Phi();
	  FPMCharge = AODFPMParticle->Charge();
	  FPMSpecies = GetSpecies(AODFPMParticle->PdgCode());
	}
      }
      Double_t fillArrayPP[7] = { PPEta, PPPt, fCentrality, fZVertex, PPPhi, PPCharge, PPSpecies };
      fHPP->Fill(fillArrayPP);
      Double_t fillArrayFPM[7] = { FPMEta, FPMPt, fCentrality, fZVertex, FPMPhi, FPMCharge, FPMSpecies };
      fHFPM->Fill(fillArrayFPM);
    } // end of generated level loop
  }
  
  // reconstructed level loop
  for (Int_t iTrack = 0; iTrack<ntrks; iTrack++)
  {
    Int_t label = 0;
    Int_t CutType = 0;
    Double_t MuonType = 0.0;
    Double_t FPMSpecies = 0.0;
    Double_t PPSpecies = 0.0;

    Double_t RecEta = 0.0;
    Double_t RecPt = 0.0;
    Double_t RecPhi = 0.0;
    Double_t RecCharge = 0.0;

    Double_t MuFPMEta = 0.0;
    Double_t MuFPMPt = 0.0;
    Double_t MuFPMPhi = 0.0;
    Double_t MuFPMCharge = 0.0;

    Double_t MuPPEta = 0.0;
    Double_t MuPPPt = 0.0;
    Double_t MuPPPhi = 0.0;
    Double_t MuPPCharge = 0.0;

    Double_t motherXv = 0.0;
    Double_t motherYv = 0.0;
    Double_t motherZv = 0.0;
    
    if(fESD)
    {
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
      if(muonTrack)
      {
	CutType = GetMUONCutType(*muonTrack);
	if(CutType > 1) continue;
	RecEta = muonTrack->Eta();
        RecPt = muonTrack->Pt();
	RecPhi = muonTrack->Phi();
	RecCharge = muonTrack->Charge();
        label =  TMath::Abs(muonTrack->GetLabel());
	if (label>=fMC->GetNumberOfTracks()) {
          AliError(Form("Label %d larger than number of particles on stack %d\n",label,fMC->GetNumberOfTracks()));
          continue;
	}
      }
    }
    else if(fAOD)
    {
      AliAODTrack* muonTrack = (AliAODTrack*)fAOD->GetTrack(iTrack);
      if(muonTrack)
      {
        if(!(muonTrack->IsMuonTrack())) continue;
	CutType =  GetMUONCutType(*muonTrack);
	if(CutType > 1) continue;
	RecEta = muonTrack->Eta();
        RecPt = muonTrack->Pt();
	RecPhi = muonTrack->Phi();
	RecCharge = muonTrack->Charge();
        label =  TMath::Abs(muonTrack->GetLabel());
        if (label>=fMC->GetNumberOfTracks()) {
          AliError(Form("Label %d larger than number of particles on stack %d\n",label,fMC->GetNumberOfTracks()));
          continue;
        }
      }
    }
    
    if(fIsMc)
    {
      AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(label);
      MuonType = GetMuonTrackType(*McParticle);
      
      if(GetFirstPrimaryMother(label) < 0) continue;
      AliMCParticle *MuFPMParticle  = (AliMCParticle*)fMC->GetTrack(GetFirstPrimaryMother(label));
      MuFPMEta = MuFPMParticle->Eta();
      MuFPMPt = MuFPMParticle->Pt();
      MuFPMPhi = MuFPMParticle->Phi();
      MuFPMCharge = MuFPMParticle->Charge();
      FPMSpecies = GetSpecies(MuFPMParticle->PdgCode());
      
      if(GetFirstPPMother(label) < 0) continue;
      AliMCParticle *MuPPParticle  = (AliMCParticle*)fMC->GetTrack(GetFirstPPMother(label));
      MuPPEta = MuPPParticle->Eta();
      MuPPPt = MuPPParticle->Pt();
      MuPPPhi = MuPPParticle->Phi();
      MuPPCharge = MuPPParticle->Charge();
      PPSpecies = GetSpecies(MuPPParticle->PdgCode());
      
      if(MuFPMParticle->GetFirstDaughter() > 0)
      {
	AliMCParticle *DaughtParticle  = (AliMCParticle*)fMC->GetTrack(MuFPMParticle->GetFirstDaughter());
	motherXv = DaughtParticle->Xv();
	motherYv = DaughtParticle->Yv();
	motherZv = DaughtParticle->Zv();
      }
  
      if(fPlotMode==1)
      {
	Double_t fillArrayDetRecMu[7] = { RecEta, RecPt, fCentrality, fZVertex, RecPhi, RecCharge, MuonType };
	fHDetRecMu[CutType]->Fill(fillArrayDetRecMu);
      }
      if(fPlotMode==2)
      {	
	Double_t fillArrayDetRecMuFPM[6] = { RecEta, RecPt, RecPhi, RecCharge, MuonType, FPMSpecies };
	fHDetRecMuFPM[CutType]->Fill(fillArrayDetRecMuFPM);
      
	Double_t fillArrayDetRecMuPP[6] = { RecEta, RecPt, RecPhi, RecCharge, MuonType, PPSpecies };
	fHDetRecMuPP[CutType]->Fill(fillArrayDetRecMuPP);
      }
      if(fPlotMode==3)
      {
	Double_t fillArrayMuFPM[6] = { MuFPMEta, MuFPMPt, MuFPMPhi, MuFPMCharge, MuonType, FPMSpecies };
	fHMuFPM[CutType]->Fill(fillArrayMuFPM);
	
	Double_t fillArrayMuPP[6] = { MuPPEta, MuPPPt, MuPPPhi, MuPPCharge, MuonType, PPSpecies };
	fHMuPP[CutType]->Fill(fillArrayMuPP);
      }

      // mother-daughter kinematic relation
      if(fMDProcess)
      {
	if(MuFPMParticle->GetFirstDaughter() > 0)
	{
	  if(CutType==0 || CutType==1)
	  {
	    fHZvRv[(Int_t)(MuonType - 0.5)]->Fill(motherZv, TMath::Sqrt(motherXv*motherXv + motherYv*motherYv));
	    fHXvYv[(Int_t)(MuonType - 0.5)]->Fill(motherXv, motherYv);
	  }
	}
	MDProcess((Int_t)(MuonType - 0.5), (Int_t)(CutType - 0.5), RecPt, RecPhi, RecEta, MuFPMPt, MuFPMPhi, MuFPMEta);
      }
    }// end of MC process
  }// end of reconstructed loop
  PostData(1, fOutputList);
  return;
}

//________________________________________________________________________
void AliMuonEffMC::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
}
//________________________________________________________________________
Bool_t AliMuonEffMC::VertexOk(TObject* obj) const
{
  // Modified from AliAnalyseLeadingTrackUE::VertexSelection()
 
  Int_t nContributors  = 0;
  Double_t zVertex     = 999;
  TString name("");
 
  if (obj->InheritsFrom("AliESDEvent")) {
    AliESDEvent* esdevt = (AliESDEvent*) obj;
    const AliESDVertex* vtx = esdevt->GetPrimaryVertex();
    if (!vtx)
      return 0;
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
  else if (obj->InheritsFrom("AliAODEvent")) {
    AliAODEvent* aodevt = (AliAODEvent*) obj;
    if (aodevt->GetNumberOfVertices() < 1)
      return 0;
    const AliAODVertex* vtx = aodevt->GetPrimaryVertex();
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
 
  // Reject if TPC-only vertex
  if (name.CompareTo("TPCVertex")==0)
    return kFALSE;
 
  // Check # contributors and range...
  if( nContributors < 1 || TMath::Abs(zVertex) > 10 ) {
    return kFALSE;
  }
 
  return kTRUE;
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetMUONCutType(AliESDMuonTrack &track)
{
  Int_t cutNum = 4;
  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  Double_t eta = track.Eta();

  if(track.ContainTrackerData()) cutNum = 3;
  if(track.ContainTrackerData() && eta > -4. && -2.5 > eta) cutNum = 2;
  if(track.ContainTrackerData() && eta > -4. && -2.5 > eta && thetaTrackAbsEnd > 2. &&  10. > thetaTrackAbsEnd) cutNum =1;
  if(track.ContainTrackerData() && eta > -4. && -2.5 > eta && thetaTrackAbsEnd > 2. &&  10. > thetaTrackAbsEnd && track.GetMatchTrigger() > 0) cutNum = 0;
  return cutNum;
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetMUONCutType(AliAODTrack &track)
{
  Int_t cutNum = 4;
  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  Double_t eta = track.Eta();

  if(track.IsMuonTrack()) cutNum = 3;
  if(track.IsMuonTrack() && eta > -4. && -2.5 > eta) cutNum = 2;
  if(track.IsMuonTrack() && eta > -4. && -2.5 > eta && thetaTrackAbsEnd > 2. &&  10. > thetaTrackAbsEnd) cutNum = 1;
  if(track.IsMuonTrack() && eta > -4. && -2.5 > eta && thetaTrackAbsEnd > 2. &&  10. > thetaTrackAbsEnd && track.GetMatchTrigger() > 0) cutNum = 0;
   
  return cutNum;
}

//________________________________________________________________________
Double_t AliMuonEffMC::GetMuonTrackType(AliMCParticle &track)
{
  if(track.GetMother() < fStack->GetNprimary() && track.PdgCode() == 13) return 0.5; 
  else if(track.GetMother() >= fStack->GetNprimary() && track.PdgCode() == 13) return 1.5;
  else return 2.5;
}

//________________________________________________________________________
Double_t AliMuonEffMC::GetSpecies(Int_t PdgCode)
{
  Int_t code = TMath::Abs(PdgCode);
  if(code==13) return 0.5;
  else if(code==211) return 1.5;
  else if(code==321) return 2.5;
  else if(code==411 || code==413 || code==421 || code==423 || code==431 || code==433 || code==10413 || code==10411 || code==10423 || code==10421 || code==10433 || code==10431 || code==20413 || code==415 || code==20423 || code==425 || code==20433 || code==435) return 3.5;
  else if(code==2212) return 4.5;
  else if(code==11) return 5.5;
  else return 6.5;
}

//________________________________________________________________________
void AliMuonEffMC::MDProcess(Int_t isprimary, Int_t cutNum, Double_t trackpt, Double_t trackphi, Double_t tracketa, Double_t motherpt, Double_t motherphi, Double_t mothereta)
{
  Int_t recptbin = -1;

  if((0. <= trackpt) && (trackpt < 0.5)) recptbin = 0;
  else if((0.5 <= trackpt) && (trackpt < 2.0)) recptbin = 1;
  else recptbin = 2;

  fHMuMotherRecPt[cutNum][isprimary]->Fill(trackpt, motherpt);
  fHMuMotherRecPhi[cutNum][isprimary]->Fill(trackphi, motherphi);
  fHMuMotherRecEta[cutNum][isprimary]->Fill(tracketa, mothereta);

  fHMuMohterPtDifRec[cutNum][isprimary][recptbin]->Fill(motherpt-trackpt);
  fHMuMohterPhiDifRec[cutNum][isprimary][recptbin]->Fill(deltaphi(motherphi-trackphi));
  fHMuMohterEtaDifRec[cutNum][isprimary][recptbin]->Fill(mothereta-tracketa);
}

//________________________________________________________________________
Double_t AliMuonEffMC::deltaphi(Double_t phi)
{
  if(phi < -1.0*TMath::Pi()) { return (phi + TMath::TwoPi()); }
  else if(phi > TMath::Pi()) { return (phi - TMath::TwoPi()); }
  else { return phi; }
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetFirstPrimaryMother(Int_t muonlabel)
{
  if(fAOD) return 1;
  else if(fESD)
  {
    AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(muonlabel);
    if(McParticle->GetMother()<fStack->GetNprimary()) return McParticle->GetMother();
    else
    {
      Int_t motherlabel = McParticle->GetMother();
      while(motherlabel > -1)
      {
	AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
	if(MotherParticle->GetMother()<fStack->GetNprimary()) break;
	else motherlabel = MotherParticle->GetMother();
      }
      AliMCParticle *FirstSecondaryMotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
      return FirstSecondaryMotherParticle->GetMother();
    }
  }
  else return -1;
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetFirstPPMother(Int_t muonlabel)
{
  if(fAOD) return 1;
  else if(fESD)
  {
    AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(muonlabel);
    if(fMC->IsPhysicalPrimary(muonlabel)) return muonlabel;
    else
    {
      Int_t motherlabel = McParticle->GetMother();
      while(motherlabel > -1)
      {
	AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
	if(fMC->IsPhysicalPrimary(motherlabel)) break;
	else motherlabel = MotherParticle->GetMother();
      }
      return motherlabel;
    }
  }
  else return -1;
}
