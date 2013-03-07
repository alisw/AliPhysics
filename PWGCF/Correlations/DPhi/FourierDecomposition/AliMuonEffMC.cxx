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
  AliAnalysisTaskSE(), fESD(0), fAOD(0), fMC(0), fStack(0), fCentrality(99), fZVertex(99), fOutputList(0x0),      
  fHEventStat(0), fHXsec(0), fHTrials(0), fHEvt(0x0), fIsMc(kTRUE), fIsPythia(kFALSE), fMDProcess(kFALSE), fFeynmanX(kFALSE), fScatFX(kFALSE), fZvProcess(kTRUE), 
  fCentralityEstimator("V0M"), fNEtaBins(15), fNpTBins(50), fNCentBins(1), fNZvtxBins(1), fNPhiBins(12), fNPBins(150),
  fHMuonParGen(0x0), fHMuonParGenSec(0x0), fHCutSpecies(0x0), fHMuonParGenP(0x0), fHMuonDetGenP(0x0), fHMuonDetRecP(0x0), 
  fHFXu(0), fHFXantiu(0), fHFXd(0), fHFXantid(0), fHFXg(0), fHFXetc(0), 
  fHaFx(0), fHbFx(0), fHabFxRatio(0), fHabDeltaFx(0), fHabRelDeltaFx(0)
{
  // Constructor
  //DefineInput(0, TChain::Class());
  //DefineOutput(1, TList::Class());
  for(Int_t i=0; i<5; i++)
  {
    fHMuonDetGen[i] = NULL;
    fHMuonDetRec[i] = NULL;
    fHHadDetRec[i] = NULL;
    fHSecDetRec[i] = NULL;
  }
  for(Int_t i=0; i<4; i++)
  {
    fHMuMotherGenPt[i] = NULL;
    fHMuMotherRecPt[i] = NULL;
    fHMuMotherGenPhi[i] = NULL;
    fHMuMotherRecPhi[i] = NULL;
    fHMuMotherGenEta[i] = NULL;
    fHMuMotherRecEta[i] = NULL;
    fHMuDCA[i] = NULL;
    fHMuZv[i] = NULL;
    fHMuRelZv[i] = NULL;
    for(Int_t j=0; j<3; j++)
    {
      fHMuMohterPhiDifGen[i][j] = NULL;
      fHMuMohterPhiDifRec[i][j] = NULL;
      fHMuMohterEtaDifGen[i][j] = NULL;
      fHMuMohterEtaDifRec[i][j] = NULL;
    }
  }
  for(Int_t i=0; i<2; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      fHFXmuonP[i][j] = NULL;
      fHFXmuonM[i][j] = NULL;
    }
  }
  for(Int_t i=0; i<3; i++)
  {
    fHabFxMu[i] = NULL;
    fHabFxRatioMu[i] = NULL;
    fHabDeltaFxMu[i] = NULL;
    fHabRelDeltaFxMu[i] = NULL;
  }
}

//________________________________________________________________________
AliMuonEffMC::AliMuonEffMC(const char *name) :
  AliAnalysisTaskSE(name), fESD(0), fAOD(0), fMC(0), fStack(0), fCentrality(99), fZVertex(99), fOutputList(0x0),      
  fHEventStat(0), fHXsec(0), fHTrials(0), fHEvt(0x0),  fIsMc(kTRUE), fIsPythia(kFALSE), fMDProcess(kFALSE), fFeynmanX(kFALSE), fScatFX(kFALSE), fZvProcess(kTRUE), 
  fCentralityEstimator("V0M"), fNEtaBins(15), fNpTBins(50), fNCentBins(1), fNZvtxBins(1), fNPhiBins(12), fNPBins(150),
  fHMuonParGen(0x0), fHMuonParGenSec(0x0), fHCutSpecies(0x0), fHMuonParGenP(0x0), fHMuonDetGenP(0x0), fHMuonDetRecP(0x0), 
  fHFXu(0), fHFXantiu(0), fHFXd(0),  fHFXantid(0), fHFXg(0), fHFXetc(0),
  fHaFx(0), fHbFx(0), fHabFxRatio(0), fHabDeltaFx(0), fHabRelDeltaFx(0)
{
  // Constructor
  for(Int_t i=0; i<5; i++)
  {
    fHMuonDetGen[i] = NULL;
    fHMuonDetRec[i] = NULL;
    fHHadDetRec[i] = NULL;
    fHSecDetRec[i] = NULL;
  }
  for(Int_t i=0; i<4; i++)
  {
    fHMuMotherGenPt[i] = NULL;
    fHMuMotherRecPt[i] = NULL;
    fHMuMotherGenPhi[i] = NULL;
    fHMuMotherRecPhi[i] = NULL;
    fHMuMotherGenEta[i] = NULL;
    fHMuMotherRecEta[i] = NULL;
    fHMuDCA[i] = NULL;
    fHMuZv[i] = NULL;
    fHMuRelZv[i] = NULL;
    for(Int_t j=0; j<3; j++)
    {
      fHMuMohterPhiDifGen[i][j] = NULL;
      fHMuMohterPhiDifRec[i][j] = NULL;
      fHMuMohterEtaDifGen[i][j] = NULL;
      fHMuMohterEtaDifRec[i][j] = NULL;
    }
  }
  for(Int_t i=0; i<2; i++)
  {
    for(Int_t j=0; j<3; j++)
    {
      fHFXmuonP[i][j] = NULL;
      fHFXmuonM[i][j] = NULL;
    }
  }
  for(Int_t i=0; i<3; i++)
  {
    fHabFxMu[i] = NULL;
    fHabFxRatioMu[i] = NULL;
    fHabDeltaFxMu[i] = NULL;
    fHabRelDeltaFxMu[i] = NULL;
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

  Int_t iTrackBin[6];
  Double_t* trackBins[6];
  const char* trackAxisTitle[6];

  Int_t iTrackBinP[6];
  Double_t* trackBinsP[6];
  const char* trackAxisTitleP[6];

  // eta
  Double_t etaBins[fNEtaBins+1];
  for(Int_t i=0; i<=fNEtaBins; i++) { etaBins[i] = (Double_t)(-4.0 + 1.5/fNEtaBins*i); }
  iTrackBin[0] = fNEtaBins;
  iTrackBinP[0] = fNEtaBins;
  trackBins[0] = etaBins;
  trackBinsP[0] = etaBins;
  trackAxisTitle[0] = "#eta";
  trackAxisTitleP[0] = "#eta";

  // p_T
  Double_t pTBins[fNpTBins+1];
  for(Int_t i=0; i<=fNpTBins; i++) { pTBins[i] = (Double_t)(5.0/fNpTBins * i); }
  iTrackBin[1] = fNpTBins;
  trackBins[1] = pTBins;
  trackAxisTitle[1] = "p_{T} (GeV/c)";

  // P
  Double_t PBins[fNPBins+1];
  for(Int_t i=0; i<=fNPBins; i++) { PBins[i] = (Double_t)(150.0/fNPBins * i); }
  iTrackBinP[1] = fNPBins;
  trackBinsP[1] = PBins;
  trackAxisTitleP[1] = "P (GeV/c)";

  // centrality
  Double_t CentBins[fNCentBins+1];
  for (Int_t i=0; i<=fNCentBins; i++) { CentBins[i] = (Double_t)(100.0/fNCentBins * i); }
  iTrackBin[2] = fNCentBins;
  iTrackBinP[2] = fNCentBins;
  trackBins[2] = CentBins;
  trackBinsP[2] = CentBins;
  trackAxisTitle[2] = "Cent";
  trackAxisTitleP[2] = "Cent";

  // Z-vertex
  Double_t ZvtxBins[fNZvtxBins+1];
  for(Int_t i=0; i<=fNZvtxBins; i++) { ZvtxBins[i] = (Double_t)(-10.0 + 20.0/fNZvtxBins * i); }
  iTrackBin[3] = fNZvtxBins;
  iTrackBinP[3] = fNZvtxBins;
  trackBins[3] = ZvtxBins;
  trackBinsP[3] = ZvtxBins;
  trackAxisTitle[3] = "Zvtx";
  trackAxisTitleP[3] = "Zvtx";

  // phi
  Double_t phiBins[fNPhiBins+1];
  for(Int_t i=0; i<=fNPhiBins; i++) { phiBins[i] = (Double_t)(TMath::TwoPi()/fNPhiBins * i); }
  iTrackBin[4] = fNPhiBins;
  iTrackBinP[4] = fNPhiBins;
  trackBins[4] = phiBins;
  trackBinsP[4] = phiBins;
  trackAxisTitle[4] = "#phi";
  trackAxisTitleP[4] = "#phi";

  // charge
  Double_t chargeBins[3] = {-10.0, 0, 10.0};
  iTrackBin[5] = 2;
  iTrackBinP[5] = 2;
  trackBins[5] = chargeBins;
  trackBinsP[5] = chargeBins;
  trackAxisTitle[5] = "charge";
  trackAxisTitleP[5] = "charge";

  const char *cutlabel[5] = {"NoCut", "Muon", "Trigger", "eta", "ThetaAbs"};
  if(fIsMc)
  {
    // THn for tracking efficiency
    fHMuonParGen = new THnF("fHMuonParGen", "", 6, iTrackBin, 0, 0);
    for (Int_t i=0; i<6; i++)
    {
      fHMuonParGen->SetBinEdges(i, trackBins[i]);
      fHMuonParGen->GetAxis(i)->SetTitle(trackAxisTitle[i]);
    }
    fHMuonParGen->Sumw2();
    fOutputList->Add(fHMuonParGen);
    
    fHMuonParGenSec = (THnF*) fHMuonParGen->Clone("fHMuonParGenSec");
    fHMuonParGenSec->Sumw2();
    fOutputList->Add(fHMuonParGenSec);
    
    for(Int_t i=0; i<5; i++)
    {
      fHMuonDetGen[i] = (THnF*)fHMuonParGen->Clone(Form("fHMuonDetGen_%s",cutlabel[i]));
      fHMuonDetGen[i]->Sumw2();
      fOutputList->Add(fHMuonDetGen[i]);
      
      fHMuonDetRec[i] = (THnF*) fHMuonParGen->Clone(Form("fHMuonDetRec_%s",cutlabel[i]));
      fHMuonDetRec[i]->Sumw2();
      fOutputList->Add(fHMuonDetRec[i]);
      
      fHHadDetRec[i] = (THnF*) fHMuonParGen->Clone(Form("fHHadDetRec_%s",cutlabel[i]));
      fHHadDetRec[i]->Sumw2();
      fOutputList->Add(fHHadDetRec[i]);
      
      fHSecDetRec[i] = (THnF*) fHMuonParGen->Clone(Form("fHSecDetRec_%s",cutlabel[i]));
      fHSecDetRec[i]->Sumw2();
      fOutputList->Add(fHSecDetRec[i]);
    }
    
    fHCutSpecies = new TH2F("fHCutSpecies", "", 4, 0, 4, 5, 0, 5);
    fOutputList->Add(fHCutSpecies);
    
    fHMuonParGenP = new THnF("fHMuonParGenP", "", 6, iTrackBinP, 0, 0);
    for (Int_t i=0; i<6; i++)
    {
      fHMuonParGenP->SetBinEdges(i, trackBinsP[i]);
      fHMuonParGenP->GetAxis(i)->SetTitle(trackAxisTitleP[i]);
    }
    fHMuonParGenP->Sumw2();
    fOutputList->Add(fHMuonParGenP);
    
    fHMuonDetGenP = (THnF*) fHMuonParGenP->Clone("fHMuonDetGenP");
    fHMuonDetGenP->Sumw2();
    fOutputList->Add(fHMuonDetGenP);
    
    fHMuonDetRecP = (THnF*) fHMuonParGenP->Clone("fHMuonDetRecP");
    fHMuonDetRecP->Sumw2();
    fOutputList->Add(fHMuonDetRecP);
    
    const char* MotherSpecies[4] = {"Pion","Kaon","D", "Etc"};
    const char *MuPt[3] = {"0510","1020","2040"};
    if(fZvProcess)
    {
      for(Int_t i=0; i<4; i++)
      {
	fHMuZv[i] = new TH1F(Form("fHMuZv_%s",MotherSpecies[i]), ";Z_{V}(cm)", 600, -500, 100);
	fOutputList->Add(fHMuZv[i]);
	fHMuRelZv[i] = new TH1F(Form("fHMuRelZv_%s",MotherSpecies[i]), ";|Z_{V,muon}-Zvtx|(cm)", 500, 0, 500);
	fOutputList->Add(fHMuRelZv[i]);
      }
    }
    if(fMDProcess)
    {
      for(Int_t i=0; i<4; i++)
      {
	fHMuMotherGenPt[i] = new TH2F(Form("fHMuMotherGenPt_%s",MotherSpecies[i]),";p_{T,muon}^{gen} (GeV/c);p_{T,mother}^{Truth} (GeV/c);",500, 0, 50, 500, 0, 50);
	fOutputList->Add(fHMuMotherGenPt[i]);
	fHMuMotherRecPt[i] = new TH2F(Form("fHMuMotherRecPt_%s",MotherSpecies[i]),";p_{T,muon}^{rec} (GeV/c);p_{T,mother}^{Truth} (GeV/c);",500, 0, 50, 500, 0, 50);
	fOutputList->Add(fHMuMotherRecPt[i]);
	fHMuMotherGenPhi[i] = new TH2F(Form("fHMuMotherGenPhi_%s",MotherSpecies[i]),";#phi_{gen};mother #phi;",100, 0, TMath::TwoPi(), 100, 0, TMath::TwoPi());
	fOutputList->Add(fHMuMotherGenPhi[i]);
	fHMuMotherRecPhi[i] = new TH2F(Form("fHMuMotherRecPhi_%s",MotherSpecies[i]),";#phi_{rec};mother #phi;",100, 0, TMath::TwoPi(), 100, 0, TMath::TwoPi());
	fOutputList->Add(fHMuMotherRecPhi[i]);
	fHMuMotherGenEta[i] = new TH2F(Form("fHMuMotherGenEta_%s",MotherSpecies[i]),";#eta_{gen};mother #eta;",100, -5., -1., 100, -5., -1.);
	fOutputList->Add(fHMuMotherGenEta[i]);
	fHMuMotherRecEta[i] = new TH2F(Form("fHMuMotherRecEta_%s",MotherSpecies[i]),";#eta_{rec};mother #eta;",100, -5., -1., 100, -5., -1.);
	fOutputList->Add(fHMuMotherRecEta[i]);
	fHMuDCA[i] =  new TH1F(Form("fHMuDCA_%s",MotherSpecies[i]), ";DCA", 100, 0, 50);
	fOutputList->Add(fHMuDCA[i]);
	
	for(Int_t j=0; j<3; j++)
	{
	  fHMuMohterPhiDifGen[i][j] = new TH1F(Form("fHMuMohterPhiDifGen_%s_%s",MotherSpecies[i], MuPt[j]),";#Delta#phi",100, -1.0*TMath::Pi(), TMath::Pi());
	  fOutputList->Add(fHMuMohterPhiDifGen[i][j]);
	  fHMuMohterPhiDifRec[i][j] = new TH1F(Form("fHMuMohterPhiDifRec_%s_%s",MotherSpecies[i], MuPt[j]),";#Delta#phi",100, -1.0*TMath::Pi(), TMath::Pi());
	  fOutputList->Add(fHMuMohterPhiDifRec[i][j]);
	  fHMuMohterEtaDifGen[i][j] = new TH1F(Form("fHMuMohterEtaDifGen_%s_%s",MotherSpecies[i], MuPt[j]),";#Delta#eta",100, -5.0, 5.0);
	  fOutputList->Add(fHMuMohterEtaDifGen[i][j]);
	  fHMuMohterEtaDifRec[i][j] = new TH1F(Form("fHMuMohterEtaDifRec_%s_%s",MotherSpecies[i], MuPt[j]),";#Delta#eta",100, -5.0, 5.0);
	  fOutputList->Add(fHMuMohterEtaDifRec[i][j]);
	}
      }
    }
    if(fFeynmanX) 
    {
      fHFXu = new TH1F("fHFXu",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHFXu);
      
      fHFXantiu = new TH1F("fHFXantiu",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHFXantiu);
      
      fHFXd = new TH1F("fHFXd",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHFXd);
      
      fHFXantid = new TH1F("fHFXantid",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHFXantid);
      
      fHFXg = new TH1F("fHFXg",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHFXg);
      
      fHFXetc = new TH1F("fHFXetc",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHFXetc);
      
      const char* muonptbin[2] = {"02","2"};
      const char *muonetabin[3] = {"4035","3530","3025"};
      for(Int_t i=0; i<2; i++)
      {
	for(Int_t j=0; j<3; j++)
	{
	  fHFXmuonP[i][j] = new TH1F(Form("fHFXmuonP_%s_%s",muonptbin[i], muonetabin[j]),";x_{F}",1000, 0.0, 1.0);
	  fOutputList->Add(fHFXmuonP[i][j]);
	  fHFXmuonM[i][j] = new TH1F(Form("fHFXmuonM_%s_%s",muonptbin[i], muonetabin[j]),";x_{F}",1000, 0.0, 1.0);
	  fOutputList->Add(fHFXmuonM[i][j]);
	}
      }
    }
    if(fScatFX)
    {
      fHaFx = new TH1F("fHaFx",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHaFx);
      fHbFx = new TH1F("fHbFx",";x_{F}",1000, 0.0, 1.0);
      fOutputList->Add(fHbFx);
      fHabFxRatio = new TH1F("fHabFxRatio",";x_{F,a}/x_{F,b}",1000, 0.0, 5.0);
      fOutputList->Add(fHabFxRatio);
      fHabDeltaFx = new TH1F("fHabDeltaFx",";#Deltax_{F}",1000, -2.0, 2.0);
      fOutputList->Add(fHabDeltaFx);
      fHabRelDeltaFx = new TH1F("fHabRelDeltaFx",";(x_{a}-x_{b})/(x_{a}+x{b})",1000, -1.0, 1.0);
      fOutputList->Add(fHabRelDeltaFx);
      
      const char* MuonPt[3] = {"Mu01","Mu12", "Mu2"};
      for(Int_t i=0; i<3; i++)
      {
	Int_t NabBin = 22;
	Double_t abBins[22+1] = { 0.0, 0.000001, 0.000002, 0.000004, 0.00001, 0.00002, 0.00004, 0.0001, 0.0002, 0.0004, 0.001, 0.002, 0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1};
	
	fHabFxMu[i] = new TH2F(Form("fHabFxMu_%s",MuonPt[i]),";x_{F,a};x_{F,b}",NabBin, abBins, NabBin, abBins);
	fOutputList->Add(fHabFxMu[i]);
	
	fHabFxRatioMu[i] = new TH1F(Form("fHabFxRatioMu_%s",MuonPt[i]),";x_{F,a}/x_{F,b}",1000, 0.0, 5.0);
	fOutputList->Add(fHabFxRatioMu[i]);
	
	fHabDeltaFxMu[i] = new TH1F(Form("fHabDeltaFxMu_%s",MuonPt[i]),";#Deltax_{F}",1000, -2.0, 2.0);
	fOutputList->Add(fHabDeltaFxMu[i]);
	
	fHabRelDeltaFxMu[i] = new TH1F(Form("fHabRelDeltaFxMu_%s",MuonPt[i]),";(x_{a}-x_{b})/(x_{a}+x{b})",1000, -1.0, 1.0);
	fOutputList->Add(fHabRelDeltaFxMu[i]);
      }
    }
  }
  else
  {
    fHMuonDetRec[0] = new THnF(Form("fHMuonDetRec_%s", cutlabel[0]), "", 6, iTrackBin, 0, 0);
    for (Int_t i=0; i<6; i++)
    {
      fHMuonDetRec[0]->SetBinEdges(i, trackBins[i]);
      fHMuonDetRec[0]->GetAxis(i)->SetTitle(trackAxisTitle[i]);
    }
    fHMuonDetRec[0]->Sumw2();
    fOutputList->Add(fHMuonDetRec[0]);

    for(Int_t i=1; i<5; i++)
    {
      
      fHMuonDetRec[i] = (THnF*) fHMuonDetRec[0]->Clone(Form("fHMuonDetRec_%s",cutlabel[i]));
      fHMuonDetRec[i]->Sumw2();
      fOutputList->Add(fHMuonDetRec[i]);
    }

   fHMuonDetRecP = new THnF("fHMuonDetRecP", "", 6, iTrackBinP, 0, 0);
    for (Int_t i=0; i<6; i++)
    {
      fHMuonDetRecP->SetBinEdges(i, trackBinsP[i]);
      fHMuonDetRecP->GetAxis(i)->SetTitle(trackAxisTitleP[i]);
    }
    fHMuonDetRecP->Sumw2();
    fOutputList->Add(fHMuonDetRecP);
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

  if ((fESD && !VertexOk(fESD)) || (fAOD && !VertexOk(fAOD))) { //AliInfo(Form("Event REJECTED. z = %.1f", fZVertex));
    return; }
  if (fCentrality > 100. || fCentrality < -1.5) { //AliInfo(Form("Event REJECTED. fCentrality = %.1f", fCentrality));
    return; }
 
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

  if(fIsMc)
  {
    // generated level loop
    for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++)
    {
      if(fAOD)
      {
	AliAODMCParticle *AodMcParticle  = (AliAODMCParticle*)fMC->GetTrack(ipart);
	if(AodMcParticle->Eta() < -4.0 || AodMcParticle->Eta() > -2.5) continue;
	
	if(TMath::Abs(AodMcParticle->PdgCode())==13)
	{
	  Double_t fillArrayParGen[6] = { AodMcParticle->Eta(), AodMcParticle->Pt(), fCentrality, fZVertex, AodMcParticle->Phi(), (Double_t)AodMcParticle->Charge() };
	  if(AodMcParticle->IsPrimary()) fHMuonParGen->Fill(fillArrayParGen);
	  else fHMuonParGenSec->Fill(fillArrayParGen);
	  
	  Double_t fillArrayParGenP[6] = { AodMcParticle->Eta(), AodMcParticle->P(), fCentrality, fZVertex, AodMcParticle->Phi(), (Double_t)AodMcParticle->Charge() };
	  if(AodMcParticle->IsPrimary()) fHMuonParGenP->Fill(fillArrayParGenP);
	}
      }
      else if(fESD)
      {
	AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(ipart);
	if(McParticle->Eta() < -4.0 || McParticle->Eta() > -2.5) continue;
	
	if(TMath::Abs(McParticle->PdgCode())==13)
	{
	  Double_t fillArrayParGen[6] = { McParticle->Eta(), McParticle->Pt(), fCentrality, fZVertex, McParticle->Phi(), (Double_t)McParticle->Charge() };
	  if(McParticle->GetMother()< fStack->GetNprimary()) fHMuonParGen->Fill(fillArrayParGen);
	  else fHMuonParGenSec->Fill(fillArrayParGen);
	  
	  Double_t fillArrayParGenP[6] = { McParticle->Eta(), McParticle->P(), fCentrality, fZVertex, McParticle->Phi(), (Double_t)McParticle->Charge() };
	  if(McParticle->GetMother()< fStack->GetNprimary()) fHMuonParGenP->Fill(fillArrayParGenP);
	}
      }
    }
  }
  
  // reconstructed level loop
  for (Int_t iTrack = 0; iTrack<ntrks; iTrack++)
  {
    Int_t label = 0;
    Double_t trackpt = 0;
    Double_t tracketa = 0;
    Double_t trackphi = 0;
    Double_t trackcharge = 0;
    Double_t trackp = 0;
    Double_t dcavalue = 0;

    Int_t cutNum = 0;   

    Double_t mcpt = 0;
    Double_t mceta = 0;
    Double_t mcphi = 0;
    Double_t mcp = 0;
    Double_t mccharge = 0;
    Double_t mcZv = 0;

    Int_t motherlabel =0;
    Int_t motherpdg = 0;
    Double_t motherpt = 0;
    Double_t mothereta = 0;
    Double_t motherphi = 0;
   
    if(fAOD)
    {
      AliAODTrack* muonTrack = (AliAODTrack*)fAOD->GetTrack(iTrack);
      if(muonTrack)
      {
        if(!(muonTrack->IsMuonTrack())) continue;
	cutNum = IsGoodMUONtrack(*muonTrack);
        trackpt = muonTrack->Pt();
        tracketa = muonTrack->Eta();
        trackphi = muonTrack->Phi();
	trackcharge = muonTrack->Charge();
	trackp = muonTrack->P();
        label =  TMath::Abs(muonTrack->GetLabel());
      }
    }
    else if(fESD)
    {
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
      if(muonTrack)
      {
	cutNum = IsGoodMUONtrack(*muonTrack);
        trackpt = muonTrack->Pt();
        tracketa = muonTrack->Eta();
        trackphi = muonTrack->Phi();
	trackcharge = muonTrack->Charge();
	trackp = muonTrack->P();
        label =  TMath::Abs(muonTrack->GetLabel());
        dcavalue = muonTrack->GetDCA();
      }
    }
    Double_t fillArrayDetRec[6] = { tracketa, trackpt, fCentrality, fZVertex, trackphi, trackcharge }; 
    Double_t fillArrayDetRecP[6] = { tracketa, trackp, fCentrality, fZVertex, trackphi, trackcharge };
    for(Int_t icut = 0; icut <= cutNum; icut++) { fHMuonDetRec[icut]->Fill(fillArrayDetRec); }
    if(cutNum==4) fHMuonDetRecP->Fill(fillArrayDetRecP);

    if(fIsMc)
    {
      if(fAOD)
      {
	AliAODMCParticle *aodMcParticle  = (AliAODMCParticle*)fMC->GetTrack(label);
	if(TMath::Abs(aodMcParticle->PdgCode()) != 13) { for(Int_t icut=0; icut<=cutNum; icut++) fHHadDetRec[icut]->Fill(fillArrayDetRec); }
	else if(!aodMcParticle->IsPrimary()) { for(Int_t icut=0; icut<=cutNum; icut++) fHSecDetRec[icut]->Fill(fillArrayDetRec); }
	
	mcpt = aodMcParticle->Pt();
	mceta = aodMcParticle->Eta();
	mcphi = aodMcParticle->Phi();
	mccharge = aodMcParticle->Charge();
	mcp = aodMcParticle->P();
	mcZv = aodMcParticle->Zv();
	motherlabel = aodMcParticle->GetMother();
	if(motherlabel > 0)
	{
	  AliAODMCParticle *aodMotherParticle  = (AliAODMCParticle*)fMC->GetTrack(motherlabel);
	  motherpdg = TMath::Abs(aodMotherParticle->PdgCode());
	  motherpt = aodMotherParticle->Pt();
	  mothereta = aodMotherParticle->Eta();
	  motherphi = aodMotherParticle->Phi();
	}      
      }
      else if(fESD)
      {
	AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(label);
	if(TMath::Abs(McParticle->PdgCode())!=13) { for(Int_t icut=0; icut<=cutNum; icut++) fHHadDetRec[icut]->Fill(fillArrayDetRec); }
	if(McParticle->GetMother()>=fStack->GetNprimary()) { for(Int_t icut=0; icut<=cutNum; icut++) fHSecDetRec[icut]->Fill(fillArrayDetRec);}
	
	mcpt = McParticle->Pt();
	mceta = McParticle->Eta();
	mcphi = McParticle->Phi();
	mccharge = McParticle->Charge();
	mcp = McParticle->P();
	mcZv = McParticle->Zv();
	motherlabel = McParticle->GetMother();
	if(motherlabel > 0)
	{
	  AliMCParticle *MotherParticle  = (AliMCParticle*)fMC->GetTrack(motherlabel);
	  motherpdg = TMath::Abs(MotherParticle->PdgCode());
	  motherpt = MotherParticle->Pt();
	  mothereta = MotherParticle->Eta();
	  motherphi = MotherParticle->Phi();
	}
      }

      Double_t fillArrayDetGen[6] = { mceta, mcpt, fCentrality, fZVertex, mcphi, mccharge };
      for(Int_t icut=0; icut <= cutNum; icut++) { fHMuonDetGen[icut]->Fill(fillArrayDetGen); }
      
      Double_t fillArrayDetGenP[6] = { mceta, mcp, fCentrality, fZVertex, mcphi, mccharge };
      if(cutNum==4) fHMuonDetGenP->Fill(fillArrayDetGenP);
    
      Int_t motherbin = GetMotherBin(motherpdg);
      for(Int_t icut=0; icut<=cutNum; icut++) fHCutSpecies->Fill((Double_t)motherbin+0.5, (Double_t)icut+0.5);

      // muon Z-Vertex process
      if(fZvProcess && cutNum==4)
      {
	fHMuZv[motherbin]->Fill(mcZv);
	fHMuRelZv[motherbin]->Fill(TMath::Abs(mcZv-fZVertex));
      }
      // mother-daughter kinematic relation
      if(fMDProcess && motherlabel>0) MDProcess(motherpdg, mcpt, mcphi, mceta, trackpt, trackphi, tracketa, motherpt, motherphi, mothereta, dcavalue);
    }
  }

  if(fIsMc)
  {
    if(fFeynmanX) FeynmanX();
    if(fScatFX) ScatFX();
  }

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
Int_t AliMuonEffMC::IsGoodMUONtrack(AliESDMuonTrack &track)
{
  // Applying track cuts for MUON tracks
  Int_t cutnum = 0;
  if(track.ContainTrackerData()) cutnum = 1;
  else return cutnum;

  if(track.GetMatchTrigger() > 0) cutnum = 2;
  else return cutnum;

  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  Double_t eta = track.Eta();

  // Eta cut
  if(eta > -4. && eta < -2.5) cutnum = 3;
  else return cutnum;

  // Theta cut at absorber end
  if(thetaTrackAbsEnd > 2. && thetaTrackAbsEnd < 10.) cutnum = 4;
  else return cutnum;
 
  return cutnum;
}

//________________________________________________________________________
Int_t AliMuonEffMC::IsGoodMUONtrack(AliAODTrack &track)
{
  Int_t cutnum = 0;
  if(track.IsMuonTrack()) cutnum = 1;
  else return cutnum;

  if (track.GetMatchTrigger() > 0) cutnum = 2;
  else return cutnum;

  Double_t thetaTrackAbsEnd = TMath::ATan(track.GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
  Double_t eta = track.Eta();

  // Eta cut
  if(eta > -4. && eta < -2.5) cutnum = 3;
  else return cutnum;

  // Theta cut at absorber end
  if(thetaTrackAbsEnd > 2. && thetaTrackAbsEnd < 10.) cutnum = 4;
  else return cutnum;

  return cutnum;
}

//________________________________________________________________________
void AliMuonEffMC::MDProcess(Int_t motherpdg, Double_t mcpt, Double_t mcphi, Double_t mceta, Double_t trackpt, Double_t trackphi, Double_t tracketa, Double_t motherpt, Double_t motherphi, Double_t mothereta, Double_t dcavalue)
{
  Int_t motherbin = -1;
  Int_t genptbin = -1;
  Int_t recptbin = -1;

  if(motherpdg==211) motherbin = 0;
  else if(motherpdg==321) motherbin = 1;
  else if(motherpdg==411 || motherpdg==413 || motherpdg==421 || motherpdg==423 || motherpdg==431 || motherpdg==433 || motherpdg==10413 || motherpdg==10411 || motherpdg==10423 || motherpdg==10421 || motherpdg==10433 || motherpdg==10431 || motherpdg==20413 || motherpdg==415 || motherpdg==20423 || motherpdg==425 || motherpdg==20433 || motherpdg==435) motherbin = 2;
  else motherbin = 3;
 
  if((0.5 <= mcpt) && (mcpt < 1.0)) genptbin = 0;
  else if((1.0 <= mcpt) && (mcpt < 2.0)) genptbin = 1;
  else genptbin = 2;

  if((0.5 <= trackpt) && (trackpt < 1.0)) recptbin = 0;
  else if((1.0 <= trackpt) && (trackpt < 2.0)) recptbin = 1;
  else recptbin = 2;

  fHMuMotherGenPt[motherbin]->Fill(mcpt, motherpt);
  fHMuMotherRecPt[motherbin]->Fill(trackpt, motherpt);
  fHMuMotherGenPhi[motherbin]->Fill(mcphi, motherphi);
  fHMuMotherRecPhi[motherbin]->Fill(trackphi, motherphi);
  fHMuMotherGenEta[motherbin]->Fill(mceta, mothereta);
  fHMuMotherRecEta[motherbin]->Fill(tracketa, mothereta);
  
  // generated
  fHMuMohterPhiDifGen[motherbin][genptbin]->Fill(deltaphi(motherphi-mcphi));
  fHMuMohterEtaDifGen[motherbin][genptbin]->Fill(mothereta-mceta);

  // reconstructed
  fHMuMohterPhiDifRec[motherbin][recptbin]->Fill(deltaphi(motherphi-trackphi));
  fHMuMohterEtaDifRec[motherbin][recptbin]->Fill(mothereta-tracketa);
  
  // DCA
  if(fESD) fHMuDCA[motherbin]->Fill(dcavalue);
}

//________________________________________________________________________
Double_t AliMuonEffMC::deltaphi(Double_t phi)
{
  if(phi < -1.0*TMath::Pi()) { return (phi + TMath::TwoPi()); }
  else if(phi > TMath::Pi()) { return (phi - TMath::TwoPi()); }
  else { return phi; }
}

//________________________________________________________________________
Int_t AliMuonEffMC::GetMotherBin(Int_t motherpdg)
{
  Int_t motherbin = -1;

  if(motherpdg==211) motherbin = 0;
  else if(motherpdg==321) motherbin = 1;
  else if(motherpdg==411 || motherpdg==413 || motherpdg==421 || motherpdg==423 || motherpdg==431 || motherpdg==433 || motherpdg==10413 || motherpdg==10411 || motherpdg==10423 || motherpdg==10421 || motherpdg==10433 || motherpdg==10431 || motherpdg==20413 || motherpdg==415 || motherpdg==20423 || motherpdg==425 || motherpdg==20433 || motherpdg==435) motherbin = 2;
  else motherbin = 3;

  return motherbin;
}

//________________________________________________________________________
void AliMuonEffMC::FeynmanX()
{  
  for (Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++)
  {
    TParticle *particle = fStack->Particle(ipart);
    if(particle->GetFirstMother()==0 || particle->GetFirstMother()==1)
    {
      Int_t pdgcode = particle->GetPdgCode();
      if(pdgcode==2) fHFXu->Fill(TMath::Abs(particle->Pz())/1380.0);
      else if(pdgcode==-2) fHFXantiu->Fill(TMath::Abs(particle->Pz())/1380.0);
      else if(pdgcode==1) fHFXd->Fill(TMath::Abs(particle->Pz())/1380.0);
      else if(pdgcode==-1) fHFXantid->Fill(TMath::Abs(particle->Pz())/1380.0);
      else if(TMath::Abs(pdgcode==21)) fHFXg->Fill(TMath::Abs(particle->Pz())/1380.0);
      else fHFXetc->Fill(TMath::Abs(particle->Pz())/1380.0);
    }
  }

  if(fESD)
  {
    for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfMuonTracks(); iTrack++)
    {
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
      if(muonTrack)
      {
        if(!IsGoodMUONtrack(*muonTrack)) continue;
	Int_t label = muonTrack->GetLabel();
        if(TMath::Abs(label) > fMC->GetNumberOfTracks()) continue;
       
        TParticle *particle = fStack->Particle(TMath::Abs(label));
        Int_t motherlabel = particle->GetFirstMother();
        if(!(motherlabel==-1) && !(motherlabel==0) && !(motherlabel==1))
        {
          while(1)
          {
            particle = fStack->Particle(motherlabel);
            motherlabel = particle->GetFirstMother();
            if(motherlabel==-1 || motherlabel==0 || motherlabel==1) break;
          }
        }
        if(motherlabel==-1) continue;

	Int_t ptbin = -1;
	Int_t etabin = -1;

	if(muonTrack->Pt() < 2.0) ptbin=0;
	else ptbin=1;

	if(muonTrack->Eta() <= -3.5) etabin=0;
	else if(-3.5<muonTrack->Eta() && muonTrack->Eta()<=-3.0) etabin=1;
	else etabin=2;

        if(motherlabel==0) fHFXmuonP[ptbin][etabin]->Fill(TMath::Abs(particle->Pz())/1380.0);
        else if(motherlabel==1) fHFXmuonM[ptbin][etabin]->Fill(TMath::Abs(particle->Pz())/1380.0);
      }
    }
  }
  else if(fAOD)
  {
    for (Int_t iTrack = 0; iTrack<fAOD->GetNTracks(); iTrack++)
    {
      AliAODTrack* muonTrack = (AliAODTrack*)fAOD->GetTrack(iTrack);
      if(muonTrack)
      {
        if(!(IsGoodMUONtrack(*muonTrack)) || !(muonTrack->IsMuonTrack())) continue;
        Int_t label = muonTrack->GetLabel();
        if(TMath::Abs(label) > fMC->GetNumberOfTracks()) continue;
       
        TParticle *particle = fStack->Particle(TMath::Abs(label));
        Int_t motherlabel = particle->GetFirstMother();
        if(!(motherlabel==-1) && !(motherlabel==0) && !(motherlabel==1))
        {
          while(1)
          {
            particle = fStack->Particle(motherlabel);
            motherlabel = particle->GetFirstMother();
            if(motherlabel==-1 || motherlabel==0 || motherlabel==1) break;
          }
        }
        if(motherlabel==-1) continue;

	Int_t ptbin = -1;
	Int_t etabin = -1;
	if(muonTrack->Pt() < 2.0) ptbin=0;
	else ptbin=1;

	if(muonTrack->Eta() <= -3.5) etabin=0;
	else if(-3.5<muonTrack->Eta() && muonTrack->Eta()<=-3.0) etabin=1;
	else etabin=2;

	if(motherlabel==0) fHFXmuonP[ptbin][etabin]->Fill(TMath::Abs(particle->Pz())/1380.0);
        else if(motherlabel==1) fHFXmuonM[ptbin][etabin]->Fill(TMath::Abs(particle->Pz())/1380.0);
      }
    }
  }
}
//________________________________________________________________________
void AliMuonEffMC::ScatFX()
{
  TParticle *parta = fStack->Particle(6);
  TParticle *partb = fStack->Particle(7);
  if(parta->GetFirstMother()!=-1 || parta->GetFirstMother()!=-1) return;
    
  Double_t xa = TMath::Abs(0.5*(parta->Energy() + partb->Energy() + parta->Pz() + partb->Pz())/1380.0);
  Double_t xb = TMath::Abs(0.5*(parta->Energy() + partb->Energy() - parta->Pz() - partb->Pz())/1380.0);

  fHaFx->Fill(xa);
  fHbFx->Fill(xb);
  fHabFxRatio->Fill(xa/xb);
  fHabDeltaFx->Fill(xa-xb);
  fHabRelDeltaFx->Fill((xa-xb)/(xa+xb));

  if(fESD)
  {
    for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfMuonTracks(); iTrack++)
    {
      Int_t MuonPtBin = 0;
      AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack);
      if(muonTrack)
      {
	if(!IsGoodMUONtrack(*muonTrack)) continue;

	if(muonTrack->Pt() < 1) MuonPtBin = 0;
	else if(muonTrack->Pt() < 2) MuonPtBin = 1;
	else MuonPtBin = 2;

	fHabFxMu[MuonPtBin]->Fill(xa, xb);
	fHabFxRatioMu[MuonPtBin]->Fill(xa/xb);
	fHabDeltaFxMu[MuonPtBin]->Fill(xa-xb);
	fHabRelDeltaFxMu[MuonPtBin]->Fill((xa-xb)/(xa+xb));
      }
    }
  }
}
