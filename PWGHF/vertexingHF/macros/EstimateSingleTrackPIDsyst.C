#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TFractionFitter.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "AliAnalysisTaskSEHFSystPID.h"

#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief macro for estimation of PID Systematic uncertainty of the single track (pions/kaons/protons)  //
// \main function: EstimateSingleTrackPIDsyst                                                           //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                                            //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

//_____________________________________________________
//GLOBAL VARIABLES 
//NOT TO EDIT
enum partPID {kElectron,kMuon,kPion,kKaon,kProton,kAll};
const int nPDGcodes = 6;
const int pdgcodes[nPDGcodes]={11,13,211,321,2212,-100};
const int pdgcolors[nPDGcodes]={kOrange+7,kGray,kRed,kBlue,kGreen+2,kBlack};
const int pdgfillcolors[nPDGcodes]={kOrange+7,kGray,kRed,kBlue,kGreen+2,kWhite};
const TString pdgnames[nPDGcodes]={"Electron","Muon","Pion","Kaon","Proton","All"};
//TO EDIT
const double ptlims[] = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};
const TString infileNameData = "LHC17pq.root"; 
const TString indirNameData = "PWGHF_D2H_SystNsigmaPID"; 
const TString inlistNameData = "coutputPIDhistos_ppMB_kINT7"; 
const TString infileNameMC = "LHC17l3b.root"; 
const TString indirNameMC = "PWGHF_D2H_SystNsigmaPID"; 
const TString inlistNameMC = "coutputPIDhistos_ppMB_kINT7"; 

//_____________________________________________________
//METHOD PROTOTYPES
int EstimateSingleTrackPIDsyst(int maxEntries=1.e8);
int GetHistoParticleIndex(short pdgcode);
int FindPtbin(float pt, const double ptlims[], int nPtbins);
void ComputeEfficiency(double num, double den, double &eff, double &effunc);
void GetTOFFractionsFromData(int whichpart, int iPt, TH1F* hFractionMC[nPDGcodes-1], TH1F* hFractionData[nPDGcodes-1], TH1F* hNsigmaMC[nPDGcodes], TH1F* hNsigmaData, TFractionFitter *&fNsigmaFitter, vector<int> &templUsed);
double PDFnsigmaTPCtot(double* nsigma, double* pars);
void PlotQAhistos(TList* listMC, TList* listData);
void DivideCanvas(TCanvas* c, int nPtbins);
void SetStyle();
void SetTH1Style(TH1F* histo, int markerstyle, int markercolor, float markersize, int linewidth, int linecolor, int fillcolor, float labelsize=-1, float titlesize=-1);

//_____________________________________________________
//METHOD IMPLEMENTATIONS
int EstimateSingleTrackPIDsyst(int maxEntries) {

  SetStyle();
  
  //********************************************************************************************************************************************//
  //get pt bins
  const int nPtbins = sizeof(ptlims)/sizeof(ptlims[0])-1;
  
  cout << "\n\n********************\n" << endl;
  cout << nPtbins << " p_T bins (GeV/c):" << endl;
  for(int iPt=0; iPt<nPtbins; iPt++) {
    cout << ptlims[iPt]<<"-"<< ptlims[iPt+1] <<endl;
  }
  
  //********************************************************************************************************************************************//
  //define histos

  //MC truth
  TH2F *hNsigmaTPCPionVsPtMCTrue, *hNsigmaTPCKaonVsPtMCTrue, *hNsigmaTPCProtonVsPtMCTrue;
  TH2F *hNsigmaTOFPionVsPtMCTrue, *hNsigmaTOFKaonVsPtMCTrue, *hNsigmaTOFProtonVsPtMCTrue;

  TH1F *hNsigmaTPCPionMCTrue[nPtbins], *hNsigmaTPCKaonMCTrue[nPtbins], *hNsigmaTPCProtonMCTrue[nPtbins];
  TH1F *hNsigmaTOFPionMCTrue[nPtbins], *hNsigmaTOFKaonMCTrue[nPtbins], *hNsigmaTOFProtonMCTrue[nPtbins];

  //MC tagged
  TH1F *hNsigmaTPCPionMCV0tag[nPtbins][nPDGcodes], *hNsigmaTPCKaonMCKinktag[nPtbins][nPDGcodes], *hNsigmaTPCKaonMCTOFtag[nPtbins][nPDGcodes], *hNsigmaTPCProtonMCV0tag[nPtbins][nPDGcodes];
  TH1F *hNsigmaTOFPionMCV0tag[nPtbins][nPDGcodes], *hNsigmaTOFKaonMCKinktag[nPtbins][nPDGcodes], *hNsigmaTOFKaonMCTPCtag[nPtbins][nPDGcodes], *hNsigmaTOFProtonMCV0tag[nPtbins][nPDGcodes];
  
  //data tagged
  TH1F *hNsigmaTPCPionDataV0tag[nPtbins], *hNsigmaTPCKaonDataKinktag[nPtbins], *hNsigmaTPCKaonDataTOFtag[nPtbins], *hNsigmaTPCProtonDataV0tag[nPtbins];
  TH1F *hNsigmaTOFPionDataV0tag[nPtbins], *hNsigmaTOFKaonDataKinktag[nPtbins], *hNsigmaTOFKaonDataTPCtag[nPtbins], *hNsigmaTOFProtonDataV0tag[nPtbins];

  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=0; iPart<nPDGcodes; iPart++) {
      hNsigmaTPCPionMCV0tag[iPt][iPart] = new TH1F(Form("hNsigmaTPCPionMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(#pi);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      hNsigmaTPCKaonMCKinktag[iPt][iPart] = new TH1F(Form("hNsigmaTPCKaonMCKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      hNsigmaTPCKaonMCTOFtag[iPt][iPart] = new TH1F(Form("hNsigmaTPCKaonMCTOFtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      hNsigmaTPCProtonMCV0tag[iPt][iPart] = new TH1F(Form("hNsigmaTPCProtonMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      
      hNsigmaTOFPionMCV0tag[iPt][iPart] = new TH1F(Form("hNsigmaTOFPionMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(#pi);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      hNsigmaTOFKaonMCKinktag[iPt][iPart] = new TH1F(Form("hNsigmaTOFKaonMCKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      hNsigmaTOFKaonMCTPCtag[iPt][iPart] = new TH1F(Form("hNsigmaTOFKaonMCTPCtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      hNsigmaTOFProtonMCV0tag[iPt][iPart] = new TH1F(Form("hNsigmaTOFProtonMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(p);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
      
      SetTH1Style(hNsigmaTPCPionMCV0tag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTPCKaonMCKinktag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTPCKaonMCTOFtag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTPCProtonMCV0tag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFPionMCV0tag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFKaonMCKinktag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFKaonMCTPCtag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFProtonMCV0tag[iPt][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
    }
    
    hNsigmaTPCPionDataV0tag[iPt] = new TH1F(Form("hNsigmaTPCPionDataV0tag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(#pi);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    hNsigmaTPCKaonDataKinktag[iPt] = new TH1F(Form("hNsigmaTPCKaonDataKinktag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    hNsigmaTPCKaonDataTOFtag[iPt] = new TH1F(Form("hNsigmaTPCKaonDataTOFtag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    hNsigmaTPCProtonDataV0tag[iPt] = new TH1F(Form("hNsigmaTPCProtonDataV0tag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    
    hNsigmaTOFPionDataV0tag[iPt] = new TH1F(Form("hNsigmaTOFPionDataV0tag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(#pi);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    hNsigmaTOFKaonDataKinktag[iPt] = new TH1F(Form("hNsigmaTOFKaonDataKinktag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    hNsigmaTOFKaonDataTPCtag[iPt] = new TH1F(Form("hNsigmaTOFKaonDataTPCtag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    hNsigmaTOFProtonDataV0tag[iPt] = new TH1F(Form("hNsigmaTOFProtonDataV0tag_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(p);Normalised entries",ptlims[iPt],ptlims[iPt+1]),1000,-50,50);
    
    SetTH1Style(hNsigmaTPCPionDataV0tag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCKaonDataKinktag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCKaonDataTOFtag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCProtonDataV0tag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFPionDataV0tag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataKinktag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataTPCtag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFProtonDataV0tag[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
  }

  //********************************************************************************************************************************************//
  //load MC inputs

  TFile* infileMC = TFile::Open(infileNameMC.Data());
  if(!infileMC || !infileMC->IsOpen()) return 2;
  TDirectoryFile* indirMC = (TDirectoryFile*)infileMC->Get(indirNameMC.Data());
  if(!indirMC) return 3;
  TList* listMC = (TList*)indirMC->Get(inlistNameMC.Data());
  hNsigmaTPCPionVsPtMCTrue   = (TH2F*)listMC->FindObject("fHistNsigmaTPCvsPt_Pion");
  hNsigmaTPCKaonVsPtMCTrue   = (TH2F*)listMC->FindObject("fHistNsigmaTPCvsPt_Kaon");
  hNsigmaTPCProtonVsPtMCTrue = (TH2F*)listMC->FindObject("fHistNsigmaTPCvsPt_Proton");
  hNsigmaTOFPionVsPtMCTrue   = (TH2F*)listMC->FindObject("fHistNsigmaTOFvsPt_Pion");
  hNsigmaTOFKaonVsPtMCTrue   = (TH2F*)listMC->FindObject("fHistNsigmaTOFvsPt_Kaon");
  hNsigmaTOFProtonVsPtMCTrue = (TH2F*)listMC->FindObject("fHistNsigmaTOFvsPt_Proton");
  TTree* treePIDMC = (TTree*)indirMC->Get("fPIDtree");
  if(!treePIDMC) return 4;
  treePIDMC->SetName("treePIDMC");

  for(int iPt=0; iPt<nPtbins; iPt++) {
    int ptbinmin = hNsigmaTPCPionVsPtMCTrue->GetXaxis()->FindBin(ptlims[iPt]*1.0001);
    int ptbinmax = hNsigmaTPCPionVsPtMCTrue->GetXaxis()->FindBin(ptlims[iPt+1]*0.9999);

    hNsigmaTPCPionMCTrue[iPt] = (TH1F*)hNsigmaTPCPionVsPtMCTrue->ProjectionY(Form("hNsigmaTPCPionMCTrue_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),ptbinmin,ptbinmax);
    hNsigmaTPCKaonMCTrue[iPt] = (TH1F*)hNsigmaTPCKaonVsPtMCTrue->ProjectionY(Form("hNsigmaTPCKaonMCTrue_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),ptbinmin,ptbinmax);
    hNsigmaTPCProtonMCTrue[iPt] = (TH1F*)hNsigmaTPCProtonVsPtMCTrue->ProjectionY(Form("hNsigmaTPCProtonMCTrue_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),ptbinmin,ptbinmax);

    hNsigmaTOFPionMCTrue[iPt] = (TH1F*)hNsigmaTOFPionVsPtMCTrue->ProjectionY(Form("hNsigmaTOFPionMCTrue_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),ptbinmin,ptbinmax);
    hNsigmaTOFKaonMCTrue[iPt] = (TH1F*)hNsigmaTOFKaonVsPtMCTrue->ProjectionY(Form("hNsigmaTOFKaonMCTrue_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),ptbinmin,ptbinmax);
    hNsigmaTOFProtonMCTrue[iPt] = (TH1F*)hNsigmaTOFProtonVsPtMCTrue->ProjectionY(Form("hNsigmaTOFProtonMCTrue_p%0.f_%0.f",ptlims[iPt]*10,ptlims[iPt+1]*10),ptbinmin,ptbinmax);

    SetTH1Style(hNsigmaTPCPionMCTrue[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCKaonMCTrue[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCProtonMCTrue[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFPionMCTrue[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonMCTrue[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFProtonMCTrue[iPt],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
  }
  
  short pdgMC             = -1;
  unsigned char tagMC     = 0;
  short n_sigma_TPC_pi_MC = -999;
  short n_sigma_TPC_K_MC  = -999;
  short n_sigma_TPC_p_MC  = -999;
  short n_sigma_TOF_pi_MC = -999;
  short n_sigma_TOF_K_MC  = -999;
  short n_sigma_TOF_p_MC  = -999;
  unsigned short pT_MC    = 0;

  treePIDMC->SetBranchAddress("PDGcode",&pdgMC);
  treePIDMC->SetBranchAddress("tag",&tagMC);
  treePIDMC->SetBranchAddress("n_sigma_TPC_pi",&n_sigma_TPC_pi_MC);
  treePIDMC->SetBranchAddress("n_sigma_TPC_K",&n_sigma_TPC_K_MC);
  treePIDMC->SetBranchAddress("n_sigma_TPC_p",&n_sigma_TPC_p_MC);
  treePIDMC->SetBranchAddress("n_sigma_TOF_pi",&n_sigma_TOF_pi_MC);
  treePIDMC->SetBranchAddress("n_sigma_TOF_K",&n_sigma_TOF_K_MC);
  treePIDMC->SetBranchAddress("n_sigma_TOF_p",&n_sigma_TOF_p_MC);
  treePIDMC->SetBranchAddress("pT",&pT_MC);

  cout << "\n********** Loop on MC tracks **********\n" << endl;

  for(int iEntry=0; iEntry<treePIDMC->GetEntriesFast(); iEntry++) {
    
    if(iEntry>maxEntries) break;

    if(iEntry%1000000==0 || iEntry==treePIDMC->GetEntriesFast()-1) cout << Form("MC Track %010d",iEntry) << endl;

    treePIDMC->GetEntry(iEntry);
    int iPt = FindPtbin(static_cast<float>(pT_MC)/1000,ptlims,nPtbins);
    if(iPt<0 || iPt>=static_cast<int>(nPtbins)) continue;

    if((tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromK0s || tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromL)) {
      hNsigmaTPCPionMCV0tag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TPC_pi_MC)/100);
      hNsigmaTOFPionMCV0tag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TOF_pi_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks) {
      hNsigmaTPCKaonMCKinktag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
      hNsigmaTOFKaonMCKinktag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF) {
      hNsigmaTPCKaonMCTOFtag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC) {
      hNsigmaTOFKaonMCTPCtag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsProtonFromL) {
      hNsigmaTPCProtonMCV0tag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TPC_p_MC)/100);
      hNsigmaTOFProtonMCV0tag[iPt][kAll]->Fill(static_cast<float>(n_sigma_TOF_p_MC)/100);
    }

    int iHisto = GetHistoParticleIndex(pdgMC);
    if(iHisto>=kElectron && iHisto<kAll) {
      if((tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromK0s || tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromL)) {
        hNsigmaTPCPionMCV0tag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TPC_pi_MC)/100);
        hNsigmaTOFPionMCV0tag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TOF_pi_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks) {
        hNsigmaTPCKaonMCKinktag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
        hNsigmaTOFKaonMCKinktag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF) {
        hNsigmaTPCKaonMCTOFtag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC) {
        hNsigmaTOFKaonMCTPCtag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsProtonFromL) {
        hNsigmaTPCProtonMCV0tag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TPC_p_MC)/100);
        hNsigmaTOFProtonMCV0tag[iPt][iHisto]->Fill(static_cast<float>(n_sigma_TOF_p_MC)/100);
      }
    }
  }
  
  //********************************************************************************************************************************************//
  //load data inputs

  TFile* infileData = TFile::Open(infileNameData.Data());
  if(!infileData || !infileData->IsOpen()) return 5;
  TDirectoryFile* indirData = (TDirectoryFile*)infileData->Get(indirNameData.Data());
  if(!indirData) return 6;
  TList* listData = (TList*)indirData->Get(inlistNameData.Data());
  TTree* treePIDData = (TTree*)indirData->Get("fPIDtree");
  if(!treePIDData) return 7;
  treePIDData->SetName("treePIDData");

  unsigned char tagData     = 0;
  short n_sigma_TPC_pi_Data = -999;
  short n_sigma_TPC_K_Data  = -999;
  short n_sigma_TPC_p_Data  = -999;
  short n_sigma_TOF_pi_Data = -999;
  short n_sigma_TOF_K_Data  = -999;
  short n_sigma_TOF_p_Data  = -999;
  unsigned short pT_Data    = 0;

  treePIDData->SetBranchAddress("tag",&tagData);
  treePIDData->SetBranchAddress("n_sigma_TPC_pi",&n_sigma_TPC_pi_Data);
  treePIDData->SetBranchAddress("n_sigma_TPC_K",&n_sigma_TPC_K_Data);
  treePIDData->SetBranchAddress("n_sigma_TPC_p",&n_sigma_TPC_p_Data);
  treePIDData->SetBranchAddress("n_sigma_TOF_pi",&n_sigma_TOF_pi_Data);
  treePIDData->SetBranchAddress("n_sigma_TOF_K",&n_sigma_TOF_K_Data);
  treePIDData->SetBranchAddress("n_sigma_TOF_p",&n_sigma_TOF_p_Data);
  treePIDData->SetBranchAddress("pT",&pT_Data);

  cout << "\n********** Loop on data tracks **********\n" << endl;

  TH2F* hNsigmaTPCKaonTPCtagged = new TH2F("hNsigmaTPCKaonTPCtagged",";#it{p}_{T} (GeV/#it{c});N_{#sigma}(K)",100,0.,10.,1000,-50,50.);
  hNsigmaTPCKaonTPCtagged->SetMarkerColor(kBlue);
  TH2F* hNsigmaTPCKaon = new TH2F("hNsigmaTPCKaon",";#it{p}_{T} (GeV/#it{c});N_{#sigma}^{TPC}(K)",100,0.,10.,1000,-50,50.);
  hNsigmaTPCKaon->SetMarkerColor(kBlack);
  TH2F* hNsigmaTOFKaonTOFtagged = new TH2F("hNsigmaTOFKaonTOFtagged",";#it{p}_{T} (GeV/#it{c});N_{#sigma}(K)",100,0.,10.,1000,-50,50.);
  hNsigmaTOFKaonTOFtagged->SetMarkerColor(kBlue);
  TH2F* hNsigmaTOFKaon = new TH2F("hNsigmaTOFKaon",";#it{p}_{T} (GeV/#it{c});N_{#sigma}^{TOF}(K)",100,0.,10.,1000,-50,50.);
  hNsigmaTOFKaon->SetMarkerColor(kBlack);
  
  for(int iEntry=0; iEntry<treePIDData->GetEntriesFast(); iEntry++) {

    if(iEntry>maxEntries) break;

    if(iEntry%1000000==0 || iEntry==treePIDData->GetEntriesFast()-1) cout << Form("Data Track %010d",iEntry) << endl;

    treePIDData->GetEntry(iEntry);
    int iPt = FindPtbin(static_cast<float>(pT_Data)/1000,ptlims,nPtbins);
    if(iPt<0 || iPt>=static_cast<int>(nPtbins)) continue;
    hNsigmaTPCKaon->Fill(static_cast<float>(pT_Data)/1000,static_cast<float>(n_sigma_TPC_K_Data)/100);
    hNsigmaTOFKaon->Fill(static_cast<float>(pT_Data)/1000,static_cast<float>(n_sigma_TOF_K_Data)/100);
    if((tagData&AliAnalysisTaskSEHFSystPID::kIsPionFromK0s || tagData&AliAnalysisTaskSEHFSystPID::kIsPionFromL)) {
      hNsigmaTPCPionDataV0tag[iPt]->Fill(static_cast<float>(n_sigma_TPC_pi_Data)/100);
      hNsigmaTOFPionDataV0tag[iPt]->Fill(static_cast<float>(n_sigma_TOF_pi_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks) {
      hNsigmaTPCKaonDataKinktag[iPt]->Fill(static_cast<float>(n_sigma_TPC_K_Data)/100);
      hNsigmaTOFKaonDataKinktag[iPt]->Fill(static_cast<float>(n_sigma_TOF_K_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF) {
      hNsigmaTPCKaonDataTOFtag[iPt]->Fill(static_cast<float>(n_sigma_TPC_K_Data)/100);
      hNsigmaTOFKaonTOFtagged->Fill(static_cast<float>(pT_Data)/1000,static_cast<float>(n_sigma_TOF_K_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC) {
      hNsigmaTOFKaonDataTPCtag[iPt]->Fill(static_cast<float>(n_sigma_TOF_K_Data)/100);
      hNsigmaTPCKaonTPCtagged->Fill(static_cast<float>(pT_Data)/1000,static_cast<float>(n_sigma_TPC_K_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsProtonFromL) {
      hNsigmaTPCProtonDataV0tag[iPt]->Fill(static_cast<float>(n_sigma_TPC_p_Data)/100);
      hNsigmaTOFProtonDataV0tag[iPt]->Fill(static_cast<float>(n_sigma_TOF_p_Data)/100);
    }
  }
  cout << "\n********************\n" << endl;
  
  //********************************************************************************************************************************************//
  //check QA plots
  PlotQAhistos(listMC,listData);
  
  TLatex* latAll = new TLatex();
  latAll->SetNDC();
  latAll->SetTextFont(42);
  latAll->SetTextColor(kBlack);
  latAll->SetTextSize(0.045);
  TLatex* latTag = new TLatex();
  latTag->SetNDC();
  latTag->SetTextFont(42);
  latTag->SetTextColor(kBlue);
  latTag->SetTextSize(0.045);

  TCanvas* cProva = new TCanvas("cProva","",1920,1080);
  cProva->Divide(2,1);
  cProva->cd(1);
  hNsigmaTPCKaon->DrawCopy("");
  hNsigmaTPCKaonTPCtagged->DrawCopy("same");
  latAll->DrawLatex(0.65,0.25,"All tags");
  latTag->DrawLatex(0.65,0.2,"TPC tag");
  cProva->cd(2);
  hNsigmaTOFKaon->DrawCopy("");
  hNsigmaTOFKaonTOFtagged->DrawCopy("same");
  latAll->DrawLatex(0.65,0.25,"All tags");
  latTag->DrawLatex(0.65,0.2,"TOF tag");
  
  cProva->SaveAs("TPC_TOF_tag.pdf");
  cProva->SaveAs("TPC_TOF_tag.png");

  //********************************************************************************************************************************************//
  //compute Fraction and contamination in MC
  double intNsigmaTPCPionMCV0tag[nPtbins][nPDGcodes];
  double intNsigmaTPCKaonMCKinktag[nPtbins][nPDGcodes];
  double intNsigmaTPCKaonMCTOFtag[nPtbins][nPDGcodes];
  double intNsigmaTPCProtonMCV0tag[nPtbins][nPDGcodes];
  double intNsigmaTOFPionMCV0tag[nPtbins][nPDGcodes];
  double intNsigmaTOFKaonMCKinktag[nPtbins][nPDGcodes];
  double intNsigmaTOFKaonMCTPCtag[nPtbins][nPDGcodes];
  double intNsigmaTOFProtonMCV0tag[nPtbins][nPDGcodes];

  double intNsigmaTPCPionDataV0tag[nPtbins];
  double intNsigmaTPCKaonDataKinktag[nPtbins];
  double intNsigmaTPCKaonDataTOFtag[nPtbins];
  double intNsigmaTPCProtonDataV0tag[nPtbins];
  double intNsigmaTOFPionDataV0tag[nPtbins];
  double intNsigmaTOFKaonDataKinktag[nPtbins];
  double intNsigmaTOFKaonDataTPCtag[nPtbins];
  double intNsigmaTOFProtonDataV0tag[nPtbins];
  
  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      intNsigmaTPCPionMCV0tag[iPt][iPart]   = hNsigmaTPCPionMCV0tag[iPt][iPart]->Integral() / hNsigmaTPCPionMCV0tag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTPCKaonMCKinktag[iPt][iPart] = hNsigmaTPCKaonMCKinktag[iPt][iPart]->Integral() / hNsigmaTPCKaonMCKinktag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTPCKaonMCTOFtag[iPt][iPart]  = hNsigmaTPCKaonMCTOFtag[iPt][iPart]->Integral() / hNsigmaTPCKaonMCTOFtag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTPCProtonMCV0tag[iPt][iPart] = hNsigmaTPCProtonMCV0tag[iPt][iPart]->Integral() / hNsigmaTPCProtonMCV0tag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTOFPionMCV0tag[iPt][iPart]   = hNsigmaTOFPionMCV0tag[iPt][iPart]->Integral() / hNsigmaTOFPionMCV0tag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTOFKaonMCKinktag[iPt][iPart] = hNsigmaTOFKaonMCKinktag[iPt][iPart]->Integral() / hNsigmaTOFKaonMCKinktag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTOFKaonMCTPCtag[iPt][iPart]  = hNsigmaTOFKaonMCTPCtag[iPt][iPart]->Integral() / hNsigmaTOFKaonMCTPCtag[iPt][iPart]->GetBinWidth(1);
      intNsigmaTOFProtonMCV0tag[iPt][iPart] = hNsigmaTOFProtonMCV0tag[iPt][iPart]->Integral() / hNsigmaTOFProtonMCV0tag[iPt][iPart]->GetBinWidth(1);
    }
    
    intNsigmaTPCPionDataV0tag[iPt] = hNsigmaTPCPionDataV0tag[iPt]->Integral() / hNsigmaTPCPionDataV0tag[iPt]->GetBinWidth(1);
    intNsigmaTPCKaonDataKinktag[iPt] = hNsigmaTPCKaonDataKinktag[iPt]->Integral() / hNsigmaTPCKaonDataKinktag[iPt]->GetBinWidth(1);
    intNsigmaTPCKaonDataTOFtag[iPt] = hNsigmaTPCKaonDataTOFtag[iPt]->Integral() / hNsigmaTPCKaonDataTOFtag[iPt]->GetBinWidth(1);
    intNsigmaTPCProtonDataV0tag[iPt] = hNsigmaTPCProtonDataV0tag[iPt]->Integral() / hNsigmaTPCProtonDataV0tag[iPt]->GetBinWidth(1);
    intNsigmaTOFPionDataV0tag[iPt] = hNsigmaTOFPionDataV0tag[iPt]->Integral() / hNsigmaTOFPionDataV0tag[iPt]->GetBinWidth(1);
    intNsigmaTOFKaonDataKinktag[iPt] = hNsigmaTOFKaonDataKinktag[iPt]->Integral() / hNsigmaTOFKaonDataKinktag[iPt]->GetBinWidth(1);
    intNsigmaTOFKaonDataTPCtag[iPt] = hNsigmaTOFKaonDataTPCtag[iPt]->Integral() / hNsigmaTOFKaonDataTPCtag[iPt]->GetBinWidth(1);
    intNsigmaTOFProtonDataV0tag[iPt] = hNsigmaTOFProtonDataV0tag[iPt]->Integral() / hNsigmaTOFProtonDataV0tag[iPt]->GetBinWidth(1);
  }
  
  TH1F *hFractionTPCPionMCV0tag[nPDGcodes-1], *hFractionTPCKaonMCKinktag[nPDGcodes-1], *hFractionTPCKaonMCTOFtag[nPDGcodes-1], *hFractionTPCProtonMCV0tag[nPDGcodes-1];
  TH1F *hFractionTOFPionMCV0tag[nPDGcodes-1], *hFractionTOFKaonMCKinktag[nPDGcodes-1], *hFractionTOFKaonMCTPCtag[nPDGcodes-1], *hFractionTOFProtonMCV0tag[nPDGcodes-1];

  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    hFractionTPCPionMCV0tag[iPart] = new TH1F(Form("hFractionTPCPionMCV0tag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTOFPionMCV0tag[iPart] = new TH1F(Form("hFractionTOFPionMCV0tag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTPCKaonMCKinktag[iPart] = new TH1F(Form("hFractionTPCKaonMCKinktag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTPCKaonMCTOFtag[iPart] = new TH1F(Form("hFractionTPCKaonMCTOFtag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTOFKaonMCKinktag[iPart] = new TH1F(Form("hFractionTOFKaonMCKinktag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTOFKaonMCTPCtag[iPart] = new TH1F(Form("hFractionTOFKaonMCTPCtag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTPCProtonMCV0tag[iPart] = new TH1F(Form("hFractionTPCProtonMCV0tag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);
    hFractionTOFProtonMCV0tag[iPart] = new TH1F(Form("hFractionTOFProtonMCV0tag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Contamination",nPtbins,ptlims);

    SetTH1Style(hFractionTPCPionMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFPionMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTPCKaonMCKinktag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTPCKaonMCTOFtag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFKaonMCKinktag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFKaonMCTPCtag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTPCProtonMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFProtonMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
  }
  
  for(int iPt=0; iPt<nPtbins; iPt++) {
    double eff=-1, unc=-1;
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      ComputeEfficiency(intNsigmaTPCPionMCV0tag[iPt][iPart],intNsigmaTPCPionMCV0tag[iPt][kAll],eff,unc);
      hFractionTPCPionMCV0tag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTPCPionMCV0tag[iPart]->SetBinError(iPt+1,unc);

      ComputeEfficiency(intNsigmaTOFPionMCV0tag[iPt][iPart],intNsigmaTOFPionMCV0tag[iPt][kAll],eff,unc);
      hFractionTOFPionMCV0tag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTOFPionMCV0tag[iPart]->SetBinError(iPt+1,unc);

      ComputeEfficiency(intNsigmaTPCKaonMCKinktag[iPt][iPart],intNsigmaTPCKaonMCKinktag[iPt][kAll],eff,unc);
      hFractionTPCKaonMCKinktag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTPCKaonMCKinktag[iPart]->SetBinError(iPt+1,unc);

      ComputeEfficiency(intNsigmaTPCKaonMCTOFtag[iPt][iPart],intNsigmaTPCKaonMCTOFtag[iPt][kAll],eff,unc);
      hFractionTPCKaonMCTOFtag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTPCKaonMCTOFtag[iPart]->SetBinError(iPt+1,unc);

      ComputeEfficiency(intNsigmaTOFKaonMCKinktag[iPt][iPart],intNsigmaTOFKaonMCKinktag[iPt][kAll],eff,unc);
      hFractionTOFKaonMCKinktag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTOFKaonMCKinktag[iPart]->SetBinError(iPt+1,unc);

      ComputeEfficiency(intNsigmaTOFKaonMCTPCtag[iPt][iPart],intNsigmaTOFKaonMCTPCtag[iPt][kAll],eff,unc);
      hFractionTOFKaonMCTPCtag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTOFKaonMCTPCtag[iPart]->SetBinError(iPt+1,unc);
    
      ComputeEfficiency(intNsigmaTPCProtonMCV0tag[iPt][iPart],intNsigmaTPCProtonMCV0tag[iPt][kAll],eff,unc);
      hFractionTPCProtonMCV0tag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTPCProtonMCV0tag[iPart]->SetBinError(iPt+1,unc);
        
      ComputeEfficiency(intNsigmaTOFProtonMCV0tag[iPt][iPart],intNsigmaTOFProtonMCV0tag[iPt][kAll],eff,unc);
      hFractionTOFProtonMCV0tag[iPart]->SetBinContent(iPt+1,eff);
      hFractionTOFProtonMCV0tag[iPart]->SetBinError(iPt+1,unc);
    }
  }
  
  TLegend* legFracMC = new TLegend(0.2,0.4,0.45,0.7);
  legFracMC->SetTextSize(0.05);
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    legFracMC->AddEntry(hFractionTOFProtonMCV0tag[iPart],pdgnames[iPart].Data(),"lpe");
  }

  TCanvas* cFractionTPCMC = new TCanvas("cFractionTPCMC","cFractionTPCMC",1920,1080);
  cFractionTPCMC->Divide(2,2);
  cFractionTPCMC->cd(1)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TPC #pi from K_{s}^{0};#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTPCMC->cd(1)->SetLogy();
  cFractionTPCMC->cd(1)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCPionMCV0tag[iPart]->Draw("same");
  legFracMC->Draw("same");
  
  cFractionTPCMC->cd(2)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TPC K from kinks;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTPCMC->cd(2)->SetLogy();
  cFractionTPCMC->cd(2)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCKaonMCKinktag[iPart]->Draw("same");
  
  cFractionTPCMC->cd(3)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TPC K TOF tagged;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTPCMC->cd(3)->SetLogy();
  cFractionTPCMC->cd(3)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCKaonMCTOFtag[iPart]->Draw("same");

  cFractionTPCMC->cd(4)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TPC p from #Lambda;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTPCMC->cd(4)->SetLogy();
  cFractionTPCMC->cd(4)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCProtonMCV0tag[iPart]->Draw("same");
  
  TCanvas* cFractionTOFMC = new TCanvas("cFractionTOFMC","cFractionTOFMC",1920,1080);
  cFractionTOFMC->Divide(2,2);
  cFractionTOFMC->cd(1)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TOF #pi from K_{s}^{0};#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTOFMC->cd(1)->SetLogy();
  cFractionTOFMC->cd(1)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFPionMCV0tag[iPart]->Draw("same");
  legFracMC->Draw("same");
  
  cFractionTOFMC->cd(2)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TOF K from kinks;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTOFMC->cd(2)->SetLogy();
  cFractionTOFMC->cd(2)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFKaonMCKinktag[iPart]->Draw("same");
  
  cFractionTOFMC->cd(3)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TOF K TPC tag;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTOFMC->cd(3)->SetLogy();
  cFractionTOFMC->cd(3)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFKaonMCTPCtag[iPart]->Draw("same");
  
  cFractionTOFMC->cd(4)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TOF p from #Lambda;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cFractionTOFMC->cd(4)->SetLogy();
  cFractionTOFMC->cd(45)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFProtonMCV0tag[iPart]->Draw("same");

  //********************************************************************************************************************************************//
  //compute fractions in data (only TOF)
  
  TFractionFitter *fNsigmaTOFKaonTPCtagFitter[nPtbins];
  TFractionFitter *fNsigmaTOFProtonV0tagFitter[nPtbins];
  TH1F *hFractionTOFKaonDataTPCtag[nPDGcodes-1];
  TH1F *hFractionTOFProtonDataV0tag[nPDGcodes-1];
  TH1F *hNsigmaTOFKaonDataTPCtag_Fit[nPtbins][nPDGcodes];
  TH1F *hNsigmaTOFProtonDataV0tag_Fit[nPtbins][nPDGcodes];
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    hFractionTOFKaonDataTPCtag[iPart] = new TH1F(Form("hFractionTOFKaonDataTPCtag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Fraction",nPtbins,ptlims);
    hFractionTOFProtonDataV0tag[iPart] = new TH1F(Form("hFractionTOFProtonDataV0tag_%s",pdgnames[iPart].Data()),";#it{p}_{T} (GeV/#it{c});Fraction",nPtbins,ptlims);
    SetTH1Style(hFractionTOFKaonDataTPCtag[iPart],kOpenSquare,pdgcolors[iPart]+1,1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFProtonDataV0tag[iPart],kOpenSquare,pdgcolors[iPart]+1,1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
  }

  TCanvas* cFitResultTOFKaonFromTPCtag = new TCanvas("cFitResultTOFKaonFromTPCtag","cFitResultTOFKaonFromTPCtag",1920,1080);
  DivideCanvas(cFitResultTOFKaonFromTPCtag,nPtbins);
  TCanvas* cFitResultTOFProtonFromV0tag = new TCanvas("cFitResultTOFProtonFromV0tag","cFitResultTOFProtonFromV0tag",1920,1080);
  DivideCanvas(cFitResultTOFProtonFromV0tag,nPtbins);
  TCanvas* cTOFFractionData = new TCanvas("cTOFFractionData","cTOFFractionData",1920,1080);
  cTOFFractionData->Divide(3,1);

  TLegend* legTOFFitter = new TLegend(0.14,0.62,0.52,0.86);
  legTOFFitter->SetTextSize(0.04);
  legTOFFitter->AddEntry(hNsigmaTOFKaonDataTPCtag[0],"Data","p");
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    legTOFFitter->AddEntry(hNsigmaTOFKaonMCTPCtag[0][iPart],Form("Templ %s",pdgnames[iPart].Data()),"l");
  }

  for(int iPt=0; iPt<nPtbins; iPt++) {

    vector<int> templUsedTOFKaonMCTPCtag;
    GetTOFFractionsFromData(kKaon,iPt,hFractionTOFKaonMCTPCtag,hFractionTOFKaonDataTPCtag,hNsigmaTOFKaonMCTPCtag[iPt],hNsigmaTOFKaonDataTPCtag[iPt],fNsigmaTOFKaonTPCtagFitter[iPt],templUsedTOFKaonMCTPCtag);

    if(templUsedTOFKaonMCTPCtag.size()>1) {
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll] = (TH1F*)fNsigmaTOFKaonTPCtagFitter[iPt]->GetPlot();
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->SetLineColor(kRed);
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->SetFillStyle(0);
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->SetFillColor(kWhite);
      
      for(unsigned int iTempl=0; iTempl<templUsedTOFKaonMCTPCtag.size(); iTempl++) {
        double frac, err;
        fNsigmaTOFKaonTPCtagFitter[iPt]->GetResult(iTempl, frac, err);
        hNsigmaTOFKaonDataTPCtag_Fit[iPt][templUsedTOFKaonMCTPCtag[iTempl]] = (TH1F*)fNsigmaTOFKaonTPCtagFitter[iPt]->GetMCPrediction(iTempl);
        hNsigmaTOFKaonDataTPCtag_Fit[iPt][templUsedTOFKaonMCTPCtag[iTempl]]->Scale(frac/hNsigmaTOFKaonDataTPCtag_Fit[iPt][templUsedTOFKaonMCTPCtag[iTempl]]->Integral()*hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->Integral());
        hNsigmaTOFKaonDataTPCtag_Fit[iPt][templUsedTOFKaonMCTPCtag[iTempl]]->SetFillColor(kWhite);
        hNsigmaTOFKaonDataTPCtag_Fit[iPt][templUsedTOFKaonMCTPCtag[iTempl]]->SetFillStyle(0);
      }
    }
    else {
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll] = (TH1F*)hNsigmaTOFKaonDataTPCtag[iPt]->Clone();
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->SetLineColor(kRed);
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->SetFillStyle(0);
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kAll]->SetFillColor(kWhite);
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kKaon] = (TH1F*)hNsigmaTOFKaonDataTPCtag[iPt]->Clone();
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kKaon]->SetFillColor(kWhite);
      hNsigmaTOFKaonDataTPCtag_Fit[iPt][kKaon]->SetFillStyle(0);
    }
    cFitResultTOFKaonFromTPCtag->cd(iPt+1)->SetLogy();
    hNsigmaTOFKaonDataTPCtag[iPt]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
    vector<int>::iterator it = find(templUsedTOFKaonMCTPCtag.begin(),templUsedTOFKaonMCTPCtag.end(),iPart);
      if(it!=templUsedTOFKaonMCTPCtag.end()) {
        hNsigmaTOFKaonDataTPCtag_Fit[iPt][iPart]->DrawCopy("hist same");
      }
      else {
        hNsigmaTOFKaonDataTPCtag_Fit[iPt][iPart]=nullptr;
      }
    }
    legTOFFitter->Draw("same");

    vector<int> templUsedTOFProtonMCV0tag;
    GetTOFFractionsFromData(kProton,iPt,hFractionTOFProtonMCV0tag,hFractionTOFProtonDataV0tag,hNsigmaTOFProtonMCV0tag[iPt],hNsigmaTOFProtonDataV0tag[iPt],fNsigmaTOFProtonV0tagFitter[iPt],templUsedTOFProtonMCV0tag);

    if(templUsedTOFProtonMCV0tag.size()>1) {
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll] = (TH1F*)fNsigmaTOFProtonV0tagFitter[iPt]->GetPlot();
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->SetLineColor(kRed);
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->SetFillStyle(0);
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->SetFillColor(kWhite);
      
      for(unsigned int iTempl=0; iTempl<templUsedTOFProtonMCV0tag.size(); iTempl++) {
        double frac, err;
        fNsigmaTOFProtonV0tagFitter[iPt]->GetResult(iTempl, frac, err);
        hNsigmaTOFProtonDataV0tag_Fit[iPt][templUsedTOFProtonMCV0tag[iTempl]] = (TH1F*)fNsigmaTOFProtonV0tagFitter[iPt]->GetMCPrediction(iTempl);
        hNsigmaTOFProtonDataV0tag_Fit[iPt][templUsedTOFProtonMCV0tag[iTempl]]->Scale(frac/hNsigmaTOFProtonDataV0tag_Fit[iPt][templUsedTOFProtonMCV0tag[iTempl]]->Integral()*hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->Integral());
        hNsigmaTOFProtonDataV0tag_Fit[iPt][templUsedTOFProtonMCV0tag[iTempl]]->SetFillColor(kWhite);
        hNsigmaTOFProtonDataV0tag_Fit[iPt][templUsedTOFProtonMCV0tag[iTempl]]->SetFillStyle(0);
      }
    }
    else {
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll] = (TH1F*)hNsigmaTOFProtonDataV0tag[iPt]->Clone();
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->SetLineColor(kRed);
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->SetFillStyle(0);
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kAll]->SetFillColor(kWhite);
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kProton] = (TH1F*)hNsigmaTOFProtonDataV0tag[iPt]->Clone();
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kProton]->SetFillColor(kWhite);
      hNsigmaTOFProtonDataV0tag_Fit[iPt][kProton]->SetFillStyle(0);
    }
    cFitResultTOFProtonFromV0tag->cd(iPt+1)->SetLogy();
    hNsigmaTOFProtonDataV0tag[iPt]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
    vector<int>::iterator it = find(templUsedTOFProtonMCV0tag.begin(),templUsedTOFProtonMCV0tag.end(),iPart);
      if(it!=templUsedTOFProtonMCV0tag.end()) {
        hNsigmaTOFProtonDataV0tag_Fit[iPt][iPart]->DrawCopy("hist same");
      }
      else {
        hNsigmaTOFProtonDataV0tag_Fit[iPt][iPart]=nullptr;
      }
    }
    legTOFFitter->Draw("same");
  }

  cTOFFractionData->cd(1)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TOF K TPC tag;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cTOFFractionData->cd(1)->SetLogy();
  cTOFFractionData->cd(1)->SetLogx();
  for(int iPart=0; iPart<=kProton; iPart++) {
    hFractionTOFKaonDataTPCtag[iPart]->DrawCopy("same");
  }
  cTOFFractionData->cd(2)->DrawFrame(ptlims[0],1.e-5,ptlims[nPtbins],10.,"TOF p from #Lambda;#it{p}_{T} (GeV/#it{c});Purity / Contamination");
  cTOFFractionData->cd(2)->SetLogy();
  cTOFFractionData->cd(2)->SetLogx();
  for(int iPart=0; iPart<=kProton; iPart++) {
    hFractionTOFProtonDataV0tag[iPart]->DrawCopy("same");
  }
  TLegend* legFracData = new TLegend(0.3,0.3,0.8,0.8);
  legFracData->SetTextSize(0.05);
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    legFracData->AddEntry(hFractionTOFProtonDataV0tag[iPart],pdgnames[iPart].Data(),"lpe");
  }
  cTOFFractionData->cd(3);
  legFracData->Draw("same");

  //********************************************************************************************************************************************//
  //normalise histograms
  
  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      hNsigmaTPCPionMCV0tag[iPt][iPart]->Scale(1./hNsigmaTPCPionMCV0tag[iPt][kAll]->Integral());
      hNsigmaTPCKaonMCKinktag[iPt][iPart]->Scale(1./hNsigmaTPCKaonMCKinktag[iPt][kAll]->Integral());
      hNsigmaTPCKaonMCTOFtag[iPt][iPart]->Scale(1./hNsigmaTPCKaonMCTOFtag[iPt][kAll]->Integral());
      hNsigmaTPCProtonMCV0tag[iPt][iPart]->Scale(1./hNsigmaTPCProtonMCV0tag[iPt][kAll]->Integral());
      hNsigmaTOFPionMCV0tag[iPt][iPart]->Scale(1./hNsigmaTOFPionMCV0tag[iPt][kAll]->Integral());
      hNsigmaTOFKaonMCKinktag[iPt][iPart]->Scale(1./hNsigmaTOFKaonMCKinktag[iPt][kAll]->Integral());
      hNsigmaTOFKaonMCTPCtag[iPt][iPart]->Scale(1./hNsigmaTOFKaonMCTPCtag[iPt][kAll]->Integral());
      hNsigmaTOFProtonMCV0tag[iPt][iPart]->Scale(1./hNsigmaTOFProtonMCV0tag[iPt][kAll]->Integral());

      if(iPart==kAll || hFractionTOFKaonDataTPCtag[iPart]->GetBinContent(iPt+1)>0) hNsigmaTOFKaonDataTPCtag_Fit[iPt][iPart]->Scale(1./hNsigmaTOFKaonDataTPCtag[iPt]->Integral());
      if(iPart==kAll || hFractionTOFProtonDataV0tag[iPart]->GetBinContent(iPt+1)>0) hNsigmaTOFProtonDataV0tag_Fit[iPt][iPart]->Scale(1./hNsigmaTOFProtonDataV0tag[iPt]->Integral());
    }

    hNsigmaTPCPionDataV0tag[iPt]->Scale(1./hNsigmaTPCPionDataV0tag[iPt]->Integral());
    hNsigmaTPCKaonDataKinktag[iPt]->Scale(1./hNsigmaTPCKaonDataKinktag[iPt]->Integral());
    hNsigmaTPCKaonDataTOFtag[iPt]->Scale(1./hNsigmaTPCKaonDataTOFtag[iPt]->Integral());
    hNsigmaTPCProtonDataV0tag[iPt]->Scale(1./hNsigmaTPCProtonDataV0tag[iPt]->Integral());
    hNsigmaTOFPionDataV0tag[iPt]->Scale(1./hNsigmaTOFPionDataV0tag[iPt]->Integral());
    hNsigmaTOFKaonDataKinktag[iPt]->Scale(1./hNsigmaTOFKaonDataKinktag[iPt]->Integral());
    hNsigmaTOFKaonDataTPCtag[iPt]->Scale(1./hNsigmaTOFKaonDataTPCtag[iPt]->Integral());
    hNsigmaTOFProtonDataV0tag[iPt]->Scale(1./hNsigmaTOFProtonDataV0tag[iPt]->Integral());
  }
  
  //********************************************************************************************************************************************//
  //draw and fit distributions

  TCanvas* cPionMCV0tagTPC = new TCanvas("cPionMCV0tagTPC","cPionMCV0tagTPC",1920,1080);
  DivideCanvas(cPionMCV0tagTPC,nPtbins);
  TCanvas* cKaonMCKinkstagTPC = new TCanvas("cKaonMCKinkstagTPC","cKaonMCKinkstagTPC",1920,1080);
  DivideCanvas(cKaonMCKinkstagTPC,nPtbins);
  TCanvas* cKaonMCTOFtagTPC = new TCanvas("cKaonMCTOFtagTPC","cKaonMCTOFtagTPC",1920,1080);
  DivideCanvas(cKaonMCTOFtagTPC,nPtbins);
  TCanvas* cProtonMCV0tagTPC = new TCanvas("cProtonMCV0tagTPC","cProtonMCV0tagTPC",1920,1080);
  DivideCanvas(cProtonMCV0tagTPC,nPtbins);
  TCanvas* cPionMCV0tagTOF = new TCanvas("cPionMCV0tagTOF","cPionMCV0tagTOF",1920,1080);
  DivideCanvas(cPionMCV0tagTOF,nPtbins);
  TCanvas* cKaonMCKinkstagTOF = new TCanvas("cKaonMCKinkstagTOF","cKaonMCKinkstagTOF",1920,1080);
  DivideCanvas(cKaonMCKinkstagTOF,nPtbins);
  TCanvas* cKaonMCTPCtagTOF = new TCanvas("cKaonMCTPCtagTOF","cKaonMCTPCtagTOF",1920,1080);
  DivideCanvas(cKaonMCTPCtagTOF,nPtbins);
  TCanvas* cProtonMCV0tagTOF = new TCanvas("cProtonMCV0tagTOF","cProtonMCV0tagTOF",1920,1080);
  DivideCanvas(cProtonMCV0tagTOF,nPtbins);

  TCanvas* cPionDataV0tagTPC = new TCanvas("cPionDataV0tagTPC","cPionDataV0tagTPC",1920,1080);
  DivideCanvas(cPionDataV0tagTPC,nPtbins);
  TCanvas* cKaonDataKinkstagTPC = new TCanvas("cKaonDataKinkstagTPC","cKaonDataKinkstagTPC",1920,1080);
  DivideCanvas(cKaonDataKinkstagTPC,nPtbins);
  TCanvas* cKaonDataTOFtagTPC = new TCanvas("cKaonDataTOFtagTPC","cKaonDataTOFtagTPC",1920,1080);
  DivideCanvas(cKaonDataTOFtagTPC,nPtbins);
  TCanvas* cProtonDataV0tagTPC = new TCanvas("cProtonDataV0tagTPC","cProtonDataV0tagTPC",1920,1080);
  DivideCanvas(cProtonDataV0tagTPC,nPtbins);
  TCanvas* cPionDataV0tagTOF = new TCanvas("cPionDataV0tagTOF","cPionDataV0tagTOF",1920,1080);
  DivideCanvas(cPionDataV0tagTOF,nPtbins);
  TCanvas* cKaonDataKinkstagTOF = new TCanvas("cKaonDataKinkstagTOF","cKaonDataKinkstagTOF",1920,1080);
  DivideCanvas(cKaonDataKinkstagTOF,nPtbins);
  TCanvas* cKaonDataTPCtagTOF = new TCanvas("cKaonDataTPCtagTOF","cKaonDataTPCtagTOF",1920,1080);
  DivideCanvas(cKaonDataTPCtagTOF,nPtbins);
  TCanvas* cProtonDataV0tagTOF = new TCanvas("cProtonDataV0tagTOF","cProtonDataV0tagTOF",1920,1080);
  DivideCanvas(cProtonDataV0tagTOF,nPtbins);

  //fit TPC nsigma
  TF1 *fNsigmaTPCPionMCV0tag[nPtbins][nPDGcodes-1], *fNsigmaTPCKaonMCKinktag[nPtbins][nPDGcodes-1], *fNsigmaTPCKaonMCTOFtag[nPtbins][nPDGcodes-1], *fNsigmaTPCProtonMCV0tag[nPtbins][nPDGcodes-1];

  TF1* fNsigmaTPCPionDataV0tag[nPtbins][nPDGcodes], *fNsigmaTPCKaonDataKinktag[nPtbins][nPDGcodes], *fNsigmaTPCKaonDataTOFtag[nPtbins][nPDGcodes], *fNsigmaTPCProtonDataV0tag[nPtbins][nPDGcodes];

  //subtract bkg template for TOF nsigma
  TH1F *hNsigmaTOFPionDataV0tag_sub[nPtbins], *hNsigmaTOFKaonDataKinktag_sub[nPtbins], *hNsigmaTOFKaonDataTPCtag_sub[nPtbins], *hNsigmaTOFProtonDataV0tag_sub[nPtbins];

  TLegend* legMCdist = new TLegend(0.14,0.62,0.44,0.86);
  legMCdist->SetTextSize(0.04);
  legMCdist->AddEntry(hNsigmaTPCPionMCV0tag[0][kAll],"MC All","p");
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    legMCdist->AddEntry(hNsigmaTPCPionMCV0tag[0][iPart],Form("MC %s",pdgnames[iPart].Data()),"f");
  }

  TLegend* legTPCFitter = new TLegend(0.14,0.62,0.48,0.86);
  legTPCFitter->SetTextSize(0.04);
  legTPCFitter->AddEntry(hNsigmaTPCPionDataV0tag[0],"Data","p");

  TLegend* legTOFPionDataV0tag = new TLegend(0.14,0.7,0.48,0.8);
  legTOFPionDataV0tag->SetTextSize(0.05);

  TLegend* legTOFKaonDataKinkstag = new TLegend(0.14,0.7,0.48,0.8);
  legTOFKaonDataKinkstag->SetTextSize(0.05);

  TLegend* legTOFKaonDataTPCtag = new TLegend(0.14,0.7,0.48,0.8);
  legTOFKaonDataTPCtag->SetTextSize(0.05);

  TLegend* legTOFProtonDataV0tag = new TLegend(0.14,0.7,0.48,0.8);
  legTOFProtonDataV0tag->SetTextSize(0.05);

  for(int iPt=0; iPt<nPtbins; iPt++) {

    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      fNsigmaTPCPionMCV0tag[iPt][iPart] = new TF1(Form("fNsigmaTPCPionMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonMCKinktag[iPt][iPart] = new TF1(Form("fNsigmaTPCKaonMCKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonMCTOFtag[iPt][iPart] = new TF1(Form("fNsigmaTPCKaonMCTOFtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
      fNsigmaTPCProtonMCV0tag[iPt][iPart] = new TF1(Form("fNsigmaTPCProtonMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
    }
    
    cPionMCV0tagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCPionMCV0tag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCPionMCV0tag[iPt][iPart]->DrawCopy("hist same");
      if(iPart==kPion || hFractionTPCPionMCV0tag[iPart]->GetBinContent(iPt+1)>1.e-5) hNsigmaTPCPionMCV0tag[iPt][iPart]->Fit(fNsigmaTPCPionMCV0tag[iPt][iPart],"L0R");
      else fNsigmaTPCPionMCV0tag[iPt][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCPionMCV0tag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCKinkstagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCKaonMCKinktag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCKaonMCKinktag[iPt][iPart]->DrawCopy("hist same");
      if(iPart==kKaon || hFractionTPCKaonMCKinktag[iPart]->GetBinContent(iPt+1)>1.e-5) hNsigmaTPCKaonMCKinktag[iPt][iPart]->Fit(fNsigmaTPCKaonMCKinktag[iPt][iPart],"L0R");
      else fNsigmaTPCKaonMCKinktag[iPt][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCKaonMCKinktag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCTOFtagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCKaonMCTOFtag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCKaonMCTOFtag[iPt][iPart]->DrawCopy("hist same");
      if(iPart==kKaon || hFractionTPCKaonMCTOFtag[iPart]->GetBinContent(iPt+1)>1.e-5) hNsigmaTPCKaonMCTOFtag[iPt][iPart]->Fit(fNsigmaTPCKaonMCTOFtag[iPt][iPart],"L0R");
      else fNsigmaTPCKaonMCTOFtag[iPt][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCKaonMCTOFtag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cProtonMCV0tagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCProtonMCV0tag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCProtonMCV0tag[iPt][iPart]->DrawCopy("hist same");
      if(iPart==kProton || hFractionTPCProtonMCV0tag[iPart]->GetBinContent(iPt+1)>1.e-5) hNsigmaTPCProtonMCV0tag[iPt][iPart]->Fit(fNsigmaTPCProtonMCV0tag[iPt][iPart],"L0R");
      else fNsigmaTPCProtonMCV0tag[iPt][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCProtonMCV0tag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    
    cPionMCV0tagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFPionMCV0tag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFPionMCV0tag[iPt][iPart]->DrawCopy("hist same");
    hNsigmaTOFPionMCV0tag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCKinkstagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFKaonMCKinktag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFKaonMCKinktag[iPt][iPart]->DrawCopy("hist same");
    hNsigmaTOFKaonMCKinktag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCTPCtagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFKaonMCTPCtag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFKaonMCTPCtag[iPt][iPart]->DrawCopy("hist same");
    hNsigmaTOFKaonMCTPCtag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cProtonMCV0tagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFProtonMCV0tag[iPt][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFProtonMCV0tag[iPt][iPart]->DrawCopy("hist same");
    hNsigmaTOFProtonMCV0tag[iPt][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    
    fNsigmaTPCPionDataV0tag[iPt][kAll] = new TF1(Form("fNsigmaTPCPionDataV0tag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);
    fNsigmaTPCKaonDataKinktag[iPt][kAll] = new TF1(Form("fNsigmaTPCKaonDataKinktag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);
    fNsigmaTPCKaonDataTOFtag[iPt][kAll] = new TF1(Form("fNsigmaTPCKaonDataTOFtag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);
    fNsigmaTPCProtonDataV0tag[iPt][kAll] = new TF1(Form("fNsigmaTPCProtonDataV0tag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);

    for(int iPart=kElectron; iPart<=kProton; iPart++) {
     
      double toll = 0.5;
      for(int iPar=0; iPar<3; iPar++) {
        fNsigmaTPCPionDataV0tag[iPt][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCPionMCV0tag[iPt][iPart]->GetParameter(iPar));
        fNsigmaTPCKaonDataKinktag[iPt][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCKaonMCKinktag[iPt][iPart]->GetParameter(iPar));
        fNsigmaTPCKaonDataTOFtag[iPt][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCKaonMCTOFtag[iPt][iPart]->GetParameter(iPar));
        fNsigmaTPCProtonDataV0tag[iPt][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCProtonMCV0tag[iPt][iPart]->GetParameter(iPar));
        
        if(iPart*3+iPar!=kPion*3+1) {
          fNsigmaTPCPionDataV0tag[iPt][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCPionMCV0tag[iPt][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCPionMCV0tag[iPt][iPart]->GetParameter(iPar))*toll,fNsigmaTPCPionMCV0tag[iPt][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCPionMCV0tag[iPt][iPart]->GetParameter(iPar))*toll);
        }
        if(iPart*3+iPar!=kKaon*3+1) {
        fNsigmaTPCKaonDataKinktag[iPt][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCKaonMCKinktag[iPt][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCKaonMCKinktag[iPt][iPart]->GetParameter(iPar))*toll,fNsigmaTPCKaonMCKinktag[iPt][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCKaonMCKinktag[iPt][iPart]->GetParameter(iPar))*toll);
        
        fNsigmaTPCKaonDataTOFtag[iPt][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCKaonMCTOFtag[iPt][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCKaonMCTOFtag[iPt][iPart]->GetParameter(iPar))*toll,fNsigmaTPCKaonMCTOFtag[iPt][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCKaonMCTOFtag[iPt][iPart]->GetParameter(iPar))*toll);
        }
        if(iPart*3+iPar!=kProton*3+1) {
        fNsigmaTPCProtonDataV0tag[iPt][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCProtonMCV0tag[iPt][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCProtonMCV0tag[iPt][iPart]->GetParameter(iPar))*toll,fNsigmaTPCProtonMCV0tag[iPt][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCProtonMCV0tag[iPt][iPart]->GetParameter(iPar))*toll);
        }
      }
    }

    hNsigmaTPCPionDataV0tag[iPt]->Fit(fNsigmaTPCPionDataV0tag[iPt][kAll],"L0R");
    hNsigmaTPCKaonDataKinktag[iPt]->Fit(fNsigmaTPCKaonDataKinktag[iPt][kAll],"L0R");
    hNsigmaTPCKaonDataTOFtag[iPt]->Fit(fNsigmaTPCKaonDataTOFtag[iPt][kAll],"L0R");
    hNsigmaTPCProtonDataV0tag[iPt]->Fit(fNsigmaTPCProtonDataV0tag[iPt][kAll],"L0R");

    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      fNsigmaTPCPionDataV0tag[iPt][iPart] = new TF1(Form("fNsigmaTPCPionDataV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonDataKinktag[iPt][iPart] = new TF1(Form("fNsigmaTPCKaonDataKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonDataTOFtag[iPt][iPart] = new TF1(Form("fNsigmaTPCKaonDataTOFtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);
      fNsigmaTPCProtonDataV0tag[iPt][iPart] = new TF1(Form("fNsigmaTPCProtonDataV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),ptlims[iPt]*10,ptlims[iPt+1]*10),"gaus",-50,50);

      for(int iPar=0; iPar<3; iPar++) {
        fNsigmaTPCPionDataV0tag[iPt][iPart]->SetParameter(iPar,fNsigmaTPCPionDataV0tag[iPt][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCKaonDataKinktag[iPt][iPart]->SetParameter(iPar,fNsigmaTPCKaonDataKinktag[iPt][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCKaonDataTOFtag[iPt][iPart]->SetParameter(iPar,fNsigmaTPCKaonDataTOFtag[iPt][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCProtonDataV0tag[iPt][iPart]->SetParameter(iPar,fNsigmaTPCProtonDataV0tag[iPt][kAll]->GetParameter(iPart*3+iPar));
      }
      fNsigmaTPCPionDataV0tag[iPt][iPart]->SetLineColor(pdgcolors[iPart]);
      fNsigmaTPCKaonDataKinktag[iPt][iPart]->SetLineColor(pdgcolors[iPart]);
      fNsigmaTPCKaonDataTOFtag[iPt][iPart]->SetLineColor(pdgcolors[iPart]);
      fNsigmaTPCProtonDataV0tag[iPt][iPart]->SetLineColor(pdgcolors[iPart]);
      if(iPt==0) legTPCFitter->AddEntry(fNsigmaTPCProtonDataV0tag[iPt][iPart],Form("Func %s",pdgnames[iPart].Data()),"l");
    }

    cPionDataV0tagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCPionDataV0tag[iPt]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCPionDataV0tag[iPt][iPart]->Draw("same");
    legTPCFitter->Draw("same");
    cKaonDataKinkstagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCKaonDataKinktag[iPt]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCKaonDataKinktag[iPt][iPart]->Draw("same");
    legTPCFitter->Draw("same");
    cKaonDataTOFtagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCKaonDataTOFtag[iPt]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCKaonDataTOFtag[iPt][iPart]->Draw("same");
    legTPCFitter->Draw("same");
    cProtonDataV0tagTPC->cd(iPt+1)->SetLogy();
    hNsigmaTPCProtonDataV0tag[iPt]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCProtonDataV0tag[iPt][iPart]->Draw("same");
    legTPCFitter->Draw("same");

    //subtract bkg template for TOF nsigma
    TString name = hNsigmaTOFPionDataV0tag[iPt]->GetName();
    name.ReplaceAll("hNsigmaTOFPionDataV0tag","hNsigmaTOFPionDataV0tag_sub");
    hNsigmaTOFPionDataV0tag_sub[iPt] = (TH1F*)hNsigmaTOFPionDataV0tag[iPt]->Clone(name.Data());
    
    name = hNsigmaTOFKaonDataKinktag[iPt]->GetName();
    name.ReplaceAll("hNsigmaTOFKaonDataKinktag","hNsigmaTOFKaonDataKinktag_sub");
    hNsigmaTOFKaonDataKinktag_sub[iPt] = (TH1F*)hNsigmaTOFKaonDataKinktag[iPt]->Clone(name.Data());

    name = hNsigmaTOFKaonDataTPCtag[iPt]->GetName();
    name.ReplaceAll("hNsigmaTOFKaonDataTPCtag","hNsigmaTOFKaonDataTPCtag_sub");
    hNsigmaTOFKaonDataTPCtag_sub[iPt] = (TH1F*)hNsigmaTOFKaonDataTPCtag[iPt]->Clone(name.Data());

    name = hNsigmaTOFProtonDataV0tag[iPt]->GetName();
    name.ReplaceAll("hNsigmaTOFProtonDataV0tag","hNsigmaTOFProtonDataV0tag_sub");
    hNsigmaTOFProtonDataV0tag_sub[iPt] = (TH1F*)hNsigmaTOFProtonDataV0tag[iPt]->Clone(name.Data());
    
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      if(iPart!=kPion) hNsigmaTOFPionDataV0tag_sub[iPt]->Add(hNsigmaTOFPionMCV0tag[iPt][iPart],-1.);
      if(iPart!=kKaon) {
        hNsigmaTOFKaonDataKinktag_sub[iPt]->Add(hNsigmaTOFKaonMCKinktag[iPt][iPart],-1.);
        if(hFractionTOFKaonDataTPCtag[iPart]->GetBinContent(iPt+1)>0) {
          hNsigmaTOFKaonDataTPCtag_sub[iPt]->Add(hNsigmaTOFKaonDataTPCtag_Fit[iPt][iPart],-1.);
        }
        else {
          hNsigmaTOFKaonDataTPCtag_sub[iPt]->Add(hNsigmaTOFKaonMCTPCtag[iPt][iPart],-1.);
        }
      }
      if(iPart!=kProton) {
        if(hFractionTOFProtonDataV0tag[iPart]->GetBinContent(iPt+1)>0) {
          hNsigmaTOFProtonDataV0tag_sub[iPt]->Add(hNsigmaTOFProtonDataV0tag_Fit[iPt][iPart],-1.);
       }
       else {
          hNsigmaTOFProtonDataV0tag_sub[iPt]->Add(hNsigmaTOFProtonMCV0tag[iPt][iPart],-1.);
        }
      }
    }
    
    SetTH1Style(hNsigmaTOFPionDataV0tag_sub[iPt],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataKinktag_sub[iPt],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataTPCtag_sub[iPt],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    SetTH1Style(hNsigmaTOFProtonDataV0tag_sub[iPt],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    if(iPt==0) {
      legTOFPionDataV0tag->AddEntry(hNsigmaTOFPionDataV0tag_sub[iPt],"Data","p");
      legTOFPionDataV0tag->AddEntry(hNsigmaTOFPionMCV0tag[iPt][kPion],"MC Pion","p");
      legTOFKaonDataKinkstag->AddEntry(hNsigmaTOFKaonDataKinktag_sub[iPt],"Data","p");
      legTOFKaonDataKinkstag->AddEntry(hNsigmaTOFKaonMCKinktag[iPt][kKaon],"MC Kaon","p");
      legTOFKaonDataTPCtag->AddEntry(hNsigmaTOFKaonDataTPCtag_sub[iPt],"Data","p");
      legTOFKaonDataTPCtag->AddEntry(hNsigmaTOFKaonMCTPCtag[iPt][kKaon],"MC Kaon","p");      
      legTOFProtonDataV0tag->AddEntry(hNsigmaTOFProtonDataV0tag_sub[iPt],"Data","p");
      legTOFProtonDataV0tag->AddEntry(hNsigmaTOFKaonMCTPCtag[iPt][kProton],"MC Proton","p");      
    }

    cPionDataV0tagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFPionDataV0tag_sub[iPt]->DrawCopy("E");
    hNsigmaTOFPionMCV0tag[iPt][kPion]->DrawCopy("hist same");
    legTOFPionDataV0tag->Draw("same");
    cKaonDataKinkstagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFKaonDataKinktag_sub[iPt]->DrawCopy("E");
    hNsigmaTOFKaonMCKinktag[iPt][kKaon]->DrawCopy("hist same");
    legTOFKaonDataKinkstag->Draw("same");
    cKaonDataTPCtagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFKaonDataTPCtag_sub[iPt]->DrawCopy("E");
    hNsigmaTOFKaonMCTPCtag[iPt][kKaon]->DrawCopy("hist same");
    legTOFKaonDataTPCtag->Draw("same");
    cProtonDataV0tagTOF->cd(iPt+1)->SetLogy();
    hNsigmaTOFProtonDataV0tag_sub[iPt]->DrawCopy("E");
    hNsigmaTOFProtonMCV0tag[iPt][kProton]->DrawCopy("hist same");
    legTOFProtonDataV0tag->Draw("same");
  }
  
  //********************************************************************************************************************************************//
  //compute efficiencies
 
  const int nEff = 3;
  int nSigma[nEff] = {3,2,1};
  int markersMC[nEff] = {kFullCircle,kFullSquare,kFullDiamond};
  int markersData[nEff] = {kOpenCircle,kOpenSquare,kOpenDiamond};
  TH1F *hEffPionTPCMCtrue[nEff], *hEffPionTOFMCtrue[nEff], *hEffKaonTPCMCtrue[nEff], *hEffKaonTOFMCtrue[nEff],*hEffProtonTPCMCtrue[nEff], *hEffProtonTOFMCtrue[nEff];
  TH1F *hEffPionTPCDataV0tag[nEff], *hEffPionTOFDataV0tag[nEff], *hEffKaonTPCDataKinktag[nEff], *hEffKaonTOFDataKinktag[nEff], *hEffKaonTPCDataTOFtag[nEff], *hEffKaonTOFDataTPCtag[nEff],*hEffProtonTPCDataV0tag[nEff], *hEffProtonTOFDataV0tag[nEff];
  TH1F *hRatioEffPionTPCDataV0tag[nEff], *hRatioEffPionTOFDataV0tag[nEff], *hRatioEffKaonTPCDataKinktag[nEff], *hRatioEffKaonTOFDataKinktag[nEff], *hRatioEffKaonTPCDataTOFtag[nEff], *hRatioEffKaonTOFDataTPCtag[nEff],*hRatioEffProtonTPCDataV0tag[nEff], *hRatioEffProtonTOFDataV0tag[nEff];
  TLegend* legEffPion = new TLegend(0.18,0.15,0.925,0.38); 
  legEffPion->SetTextSize(0.045);
  legEffPion->SetNColumns(2);
  legEffPion->AddEntry("","MC","");
  legEffPion->AddEntry("","V0 tag","");
  TLegend* legEffKaonTPC = new TLegend(0.18,0.15,0.925,0.38); 
  legEffKaonTPC->SetTextSize(0.045);
  legEffKaonTPC->SetNColumns(3);
  legEffKaonTPC->AddEntry("","MC","");
  legEffKaonTPC->AddEntry("","TOF tag","");
  legEffKaonTPC->AddEntry("","kink tag","");
  TLegend* legEffKaonTOF = new TLegend(0.18,0.15,0.925,0.38); 
  legEffKaonTOF->SetTextSize(0.045);
  legEffKaonTOF->SetNColumns(2);
  legEffKaonTOF->AddEntry("","MC","");
  legEffKaonTOF->AddEntry("","TPC tag","");
  TLegend* legEffProton = new TLegend(0.18,0.15,0.925,0.38); 
  legEffProton->SetTextSize(0.045);
  legEffProton->SetNColumns(2);
  legEffProton->AddEntry("","MC","");
  legEffProton->AddEntry("","V0 tag","");

  for(int iEff=0; iEff<nEff; iEff++) {
    hEffPionTPCMCtrue[iEff] = new TH1F(Form("hEffPionTPCMCtrue_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC #pi efficiency",nPtbins,ptlims);
    hEffPionTOFMCtrue[iEff] = new TH1F(Form("hEffPionTOFMCtrue_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF #pi efficiency",nPtbins,ptlims);

    hEffKaonTPCMCtrue[iEff] = new TH1F(Form("hEffKaonTPCMCtrue_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC K efficiency",nPtbins,ptlims);
    hEffKaonTOFMCtrue[iEff] = new TH1F(Form("hEffKaonTOFMCtrue_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF K efficiency",nPtbins,ptlims);
    
    hEffProtonTPCMCtrue[iEff] = new TH1F(Form("hEffProtonTPCMCtrue_%dsigma;",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC p efficiency",nPtbins,ptlims);
    hEffProtonTOFMCtrue[iEff] = new TH1F(Form("hEffProtonTOFMCtrue_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF p efficiency",nPtbins,ptlims);

    hEffPionTPCDataV0tag[iEff] = new TH1F(Form("hEffPionTPCDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC #pi efficiency",nPtbins,ptlims);
    hEffPionTOFDataV0tag[iEff] = new TH1F(Form("hEffPionTOFDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF #pi efficiency",nPtbins,ptlims);
    
    hEffKaonTPCDataKinktag[iEff] = new TH1F(Form("hEffKaonTPCDataKinktag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC K efficiency",nPtbins,ptlims);
    hEffKaonTOFDataKinktag[iEff] = new TH1F(Form("hEffKaonTOFDataKinktag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF K efficiency",nPtbins,ptlims);

    hEffKaonTPCDataTOFtag[iEff] = new TH1F(Form("hEffKaonTPCDataTOFtag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC K efficiency",nPtbins,ptlims);
    hEffKaonTOFDataTPCtag[iEff] = new TH1F(Form("hEffKaonTOFDataTPCtag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF K efficiency",nPtbins,ptlims);

    hEffProtonTPCDataV0tag[iEff] = new TH1F(Form("hEffProtonTPCDataV0tag_%dsigma;",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC p efficiency",nPtbins,ptlims);
    hEffProtonTOFDataV0tag[iEff] = new TH1F(Form("hEffProtonTOFDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF p efficiency",nPtbins,ptlims);

    hRatioEffPionTPCDataV0tag[iEff] = new TH1F(Form("hRatioEffPionTPCDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC #pi efficiency",nPtbins,ptlims);
    hRatioEffPionTOFDataV0tag[iEff] = new TH1F(Form("hRatioEffPionTOFDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF #pi efficiency",nPtbins,ptlims);

    hRatioEffKaonTPCDataKinktag[iEff] = new TH1F(Form("hRatioEffKaonTPCDataKinktag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC K efficiency",nPtbins,ptlims);
    hRatioEffKaonTOFDataKinktag[iEff] = new TH1F(Form("hRatioEffKaonTOFDataKinktag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF K efficiency",nPtbins,ptlims);

    hRatioEffKaonTPCDataTOFtag[iEff] = new TH1F(Form("hRatioEffKaonTPCDataTOFtag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC K efficiency",nPtbins,ptlims);
    hRatioEffKaonTOFDataTPCtag[iEff] = new TH1F(Form("hRatioEffKaonTOFDataTPCtag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF K efficiency",nPtbins,ptlims);

    hRatioEffProtonTPCDataV0tag[iEff] = new TH1F(Form("hRatioEffProtonTPCDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TPC #pi efficiency",nPtbins,ptlims);
    hRatioEffProtonTOFDataV0tag[iEff] = new TH1F(Form("hRatioEffProtonTOFDataV0tag_%dsigma",nSigma[iEff]),";#it{p}_{T} (GeV/#it{c});TOF #pi efficiency",nPtbins,ptlims);

    int nsigmabinlow = hNsigmaTOFPionMCTrue[0]->GetXaxis()->FindBin(-nSigma[iEff]+1.e-5);
    int nsigmabinhigh = hNsigmaTOFPionMCTrue[0]->GetXaxis()->FindBin(nSigma[iEff]-1.e-5);

    legEffPion->AddEntry(hEffPionTPCMCtrue[iEff],Form("|N_{#sigma}(#pi)| < %d",nSigma[iEff]),"p");
    legEffPion->AddEntry(hEffPionTPCDataV0tag[iEff],Form("|N_{#sigma}(#pi)| < %d",nSigma[iEff]),"p");
    legEffKaonTPC->AddEntry(hEffKaonTPCMCtrue[iEff],Form("|N_{#sigma}(K)| < %d",nSigma[iEff]),"p");
    legEffKaonTPC->AddEntry(hEffKaonTPCDataTOFtag[iEff],Form("|N_{#sigma}(K)| < %d",nSigma[iEff]),"p");
    legEffKaonTPC->AddEntry(hEffKaonTPCDataKinktag[iEff],Form("|N_{#sigma}(K)| < %d",nSigma[iEff]),"p");
    legEffKaonTOF->AddEntry(hEffKaonTPCMCtrue[iEff],Form("|N_{#sigma}(K)| < %d",nSigma[iEff]),"p");
    legEffKaonTOF->AddEntry(hEffKaonTPCDataTOFtag[iEff],Form("|N_{#sigma}(K)| < %d",nSigma[iEff]),"p");
    legEffProton->AddEntry(hEffProtonTPCMCtrue[iEff],Form("|N_{#sigma}(p)| < %d",nSigma[iEff]),"p");
    legEffProton->AddEntry(hEffProtonTOFDataV0tag[iEff],Form("|N_{#sigma}(p)| < %d",nSigma[iEff]),"p");

    for(int iPt=0; iPt<nPtbins; iPt++) {
      double eff=-1, unc=-1;
      ComputeEfficiency(hNsigmaTPCPionMCTrue[iPt]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTPCPionMCTrue[iPt]->Integral(),eff,unc);
      hEffPionTPCMCtrue[iEff]->SetBinContent(iPt+1,eff);
      hEffPionTPCMCtrue[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(hNsigmaTOFPionMCTrue[iPt]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTOFPionMCTrue[iPt]->Integral(),eff,unc);
      hEffPionTOFMCtrue[iEff]->SetBinContent(iPt+1,eff);
      hEffPionTOFMCtrue[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(hNsigmaTPCKaonMCTrue[iPt]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTPCKaonMCTrue[iPt]->Integral(),eff,unc);
      hEffKaonTPCMCtrue[iEff]->SetBinContent(iPt+1,eff);
      hEffKaonTPCMCtrue[iEff]->SetBinError(iPt+1,unc);
      
      ComputeEfficiency(hNsigmaTOFKaonMCTrue[iPt]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTOFKaonMCTrue[iPt]->Integral(),eff,unc);
      hEffKaonTOFMCtrue[iEff]->SetBinContent(iPt+1,eff);
      hEffKaonTOFMCtrue[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(hNsigmaTPCProtonMCTrue[iPt]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTPCProtonMCTrue[iPt]->Integral(),eff,unc);
      hEffProtonTPCMCtrue[iEff]->SetBinContent(iPt+1,eff);
      hEffProtonTPCMCtrue[iEff]->SetBinError(iPt+1,unc);
      
      ComputeEfficiency(hNsigmaTOFProtonMCTrue[iPt]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTOFProtonMCTrue[iPt]->Integral(),eff,unc);
      hEffProtonTOFMCtrue[iEff]->SetBinContent(iPt+1,eff);
      hEffProtonTOFMCtrue[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(fNsigmaTPCPionDataV0tag[iPt][kPion]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCPionDataV0tag[iPt],fNsigmaTPCPionDataV0tag[iPt][kPion]->Integral(-50,50)*intNsigmaTPCPionDataV0tag[iPt],eff,unc);
      hEffPionTPCDataV0tag[iEff]->SetBinContent(iPt+1,eff);
      hEffPionTPCDataV0tag[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(fNsigmaTPCKaonDataKinktag[iPt][kKaon]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCKaonDataKinktag[iPt],fNsigmaTPCKaonDataKinktag[iPt][kKaon]->Integral(-50,50)*intNsigmaTPCKaonDataKinktag[iPt],eff,unc);
      hEffKaonTPCDataKinktag[iEff]->SetBinContent(iPt+1,eff);
      hEffKaonTPCDataKinktag[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(fNsigmaTPCKaonDataTOFtag[iPt][kKaon]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCKaonDataTOFtag[iPt],fNsigmaTPCKaonDataTOFtag[iPt][kKaon]->Integral(-50,50)*intNsigmaTPCKaonDataTOFtag[iPt],eff,unc);
      hEffKaonTPCDataTOFtag[iEff]->SetBinContent(iPt+1,eff);
      hEffKaonTPCDataTOFtag[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(fNsigmaTPCProtonDataV0tag[iPt][kProton]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCProtonDataV0tag[iPt],fNsigmaTPCProtonDataV0tag[iPt][kProton]->Integral(-50,50)*intNsigmaTPCProtonDataV0tag[iPt],eff,unc);
      hEffProtonTPCDataV0tag[iEff]->SetBinContent(iPt+1,eff);
      hEffProtonTPCDataV0tag[iEff]->SetBinError(iPt+1,unc);
      
      ComputeEfficiency(hNsigmaTOFPionDataV0tag_sub[iPt]->Integral(nsigmabinlow,nsigmabinhigh)*intNsigmaTOFPionDataV0tag[iPt],hNsigmaTOFPionDataV0tag_sub[iPt]->Integral()*intNsigmaTOFPionDataV0tag[iPt],eff,unc);
      hEffPionTOFDataV0tag[iEff]->SetBinContent(iPt+1,eff);
      hEffPionTOFDataV0tag[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(hNsigmaTOFKaonDataTPCtag_sub[iPt]->Integral(nsigmabinlow,nsigmabinhigh)*intNsigmaTOFKaonDataTPCtag[iPt],hNsigmaTOFKaonDataTPCtag_sub[iPt]->Integral()*intNsigmaTOFKaonDataTPCtag[iPt],eff,unc);
      hEffKaonTOFDataTPCtag[iEff]->SetBinContent(iPt+1,eff);
      hEffKaonTOFDataTPCtag[iEff]->SetBinError(iPt+1,unc);

      ComputeEfficiency(hNsigmaTOFProtonDataV0tag_sub[iPt]->Integral(nsigmabinlow,nsigmabinhigh)*intNsigmaTOFProtonDataV0tag[iPt],hNsigmaTOFProtonDataV0tag_sub[iPt]->Integral()*intNsigmaTOFProtonDataV0tag[iPt],eff,unc);
      hEffProtonTOFDataV0tag[iEff]->SetBinContent(iPt+1,eff);
      hEffProtonTOFDataV0tag[iEff]->SetBinError(iPt+1,unc);
    }

    hRatioEffPionTPCDataV0tag[iEff]->Divide(hEffPionTPCDataV0tag[iEff],hEffPionTPCMCtrue[iEff],1.,1.,"");
    hRatioEffPionTOFDataV0tag[iEff]->Divide(hEffPionTOFDataV0tag[iEff],hEffPionTOFMCtrue[iEff],1.,1.,"");
    hRatioEffKaonTPCDataKinktag[iEff]->Divide(hEffKaonTPCDataKinktag[iEff],hEffKaonTPCMCtrue[iEff],1.,1.,"");
    hRatioEffKaonTPCDataTOFtag[iEff]->Divide(hEffKaonTPCDataTOFtag[iEff],hEffKaonTPCMCtrue[iEff],1.,1.,"");
    hRatioEffKaonTOFDataTPCtag[iEff]->Divide(hEffKaonTOFDataTPCtag[iEff],hEffKaonTOFMCtrue[iEff],1.,1.,"");
    hRatioEffProtonTPCDataV0tag[iEff]->Divide(hEffProtonTPCDataV0tag[iEff],hEffProtonTPCMCtrue[iEff],1.,1.,"");
    hRatioEffProtonTOFDataV0tag[iEff]->Divide(hEffProtonTOFDataV0tag[iEff],hEffProtonTOFMCtrue[iEff],1.,1.,"");

    SetTH1Style(hEffPionTPCMCtrue[iEff],markersMC[iEff],pdgcolors[kPion],1.,2,pdgcolors[kPion],kWhite,0.045,0.055);
    SetTH1Style(hEffPionTOFMCtrue[iEff],markersMC[iEff],pdgcolors[kPion],1.,2,pdgcolors[kPion],kWhite,0.045,0.055);
    SetTH1Style(hEffKaonTPCMCtrue[iEff],markersMC[iEff],pdgcolors[kKaon],1.,2,pdgcolors[kKaon],kWhite,0.045,0.055);
    SetTH1Style(hEffKaonTOFMCtrue[iEff],markersMC[iEff],pdgcolors[kKaon],1.,2,pdgcolors[kKaon],kWhite,0.045,0.055);
    SetTH1Style(hEffProtonTPCMCtrue[iEff],markersMC[iEff],pdgcolors[kProton],1.,2,pdgcolors[kProton],kWhite,0.045,0.055);
    SetTH1Style(hEffProtonTOFMCtrue[iEff],markersMC[iEff],pdgcolors[kProton],1.,2,pdgcolors[kProton],kWhite,0.045,0.055);

    SetTH1Style(hEffPionTPCDataV0tag[iEff],markersData[iEff],pdgcolors[kPion]+1,1.,2,pdgcolors[kPion]+1,kWhite,0.045,0.055);
    SetTH1Style(hEffPionTOFDataV0tag[iEff],markersData[iEff],pdgcolors[kPion]+1,1.,2,pdgcolors[kPion]+1,kWhite,0.045,0.055);
    SetTH1Style(hEffKaonTPCDataKinktag[iEff],markersData[iEff],pdgcolors[kKaon]+1,1.,2,pdgcolors[kKaon]+1,kWhite,0.045,0.055);
    SetTH1Style(hEffKaonTPCDataTOFtag[iEff],markersData[iEff],pdgcolors[kKaon]+3,1.,2,pdgcolors[kKaon]+3,kWhite,0.045,0.055);
    SetTH1Style(hEffKaonTOFDataTPCtag[iEff],markersData[iEff],pdgcolors[kKaon]+3,1.,2,pdgcolors[kKaon]+3,kWhite,0.045,0.055);
    SetTH1Style(hEffProtonTPCDataV0tag[iEff],markersData[iEff],pdgcolors[kProton]+1,1.,2,pdgcolors[kProton]+1,kWhite,0.045,0.055);
    SetTH1Style(hEffProtonTOFDataV0tag[iEff],markersData[iEff],pdgcolors[kProton]+1,1.,2,pdgcolors[kProton]+1,kWhite,0.045,0.055);

    SetTH1Style(hRatioEffPionTPCDataV0tag[iEff],markersData[iEff],pdgcolors[kPion]+1,1.,2,pdgcolors[kPion]+1,kWhite,0.045,0.055);
    SetTH1Style(hRatioEffPionTOFDataV0tag[iEff],markersData[iEff],pdgcolors[kPion]+1,1.,2,pdgcolors[kPion]+1,kWhite,0.045,0.055);
    SetTH1Style(hRatioEffKaonTPCDataKinktag[iEff],markersData[iEff],pdgcolors[kKaon]+1,1.,2,pdgcolors[kKaon]+1,kWhite,0.045,0.055);
    SetTH1Style(hRatioEffKaonTPCDataTOFtag[iEff],markersData[iEff],pdgcolors[kKaon]+3,1.,2,pdgcolors[kKaon]+3,kWhite,0.045,0.055);
    SetTH1Style(hRatioEffKaonTOFDataTPCtag[iEff],markersData[iEff],pdgcolors[kKaon]+3,1.,2,pdgcolors[kKaon]+3,kWhite,0.045,0.055);
    SetTH1Style(hRatioEffProtonTPCDataV0tag[iEff],markersData[iEff],pdgcolors[kProton]+1,1.,2,pdgcolors[kProton]+1,kWhite,0.045,0.055);
    SetTH1Style(hRatioEffProtonTOFDataV0tag[iEff],markersData[iEff],pdgcolors[kProton]+1,1.,2,pdgcolors[kProton]+1,kWhite,0.045,0.055);
  }

  TCanvas* cEffPion = new TCanvas("cEffPion","cEffPion",800,800);
  TCanvas* cEffKaon = new TCanvas("cEffKaon","cEffKaon",800,800);
  TCanvas* cEffProton = new TCanvas("cEffProton","cEffProton",800,800);
  cEffPion->Divide(2,2);
  cEffPion->cd(1)->DrawFrame(ptlims[0],0.,ptlims[nPtbins],1.,";#it{p}_{T} (GeV/#it{c});pion TPC PID efficiency");
  cEffPion->cd(2)->DrawFrame(ptlims[0],0.,ptlims[nPtbins],1.,";#it{p}_{T} (GeV/#it{c});pion TOF PID efficiency");
  cEffPion->cd(3)->DrawFrame(ptlims[0],0.85,ptlims[nPtbins],1.15,";#it{p}_{T} (GeV/#it{c}); pion data / MC TPC PID efficiency");
  cEffPion->cd(4)->DrawFrame(ptlims[0],0.85,ptlims[nPtbins],1.15,";#it{p}_{T} (GeV/#it{c});pion data / MC TOF PID efficiency");
  cEffKaon->Divide(2,2);
  cEffKaon->cd(1)->DrawFrame(ptlims[0],0.,ptlims[nPtbins],1.,";#it{p}_{T} (GeV/#it{c});kaon TPC PID efficiency");
  cEffKaon->cd(2)->DrawFrame(ptlims[0],0.,ptlims[nPtbins],1.,";#it{p}_{T} (GeV/#it{c});kaon TOF PID efficiency");
  cEffKaon->cd(3)->DrawFrame(ptlims[0],0.85,ptlims[nPtbins],1.15,";#it{p}_{T} (GeV/#it{c});kaon data / MC TPC PID efficiency");
  cEffKaon->cd(4)->DrawFrame(ptlims[0],0.85,ptlims[nPtbins],1.15,";#it{p}_{T} (GeV/#it{c});kaon data / MC TOF PID efficiency");
  cEffProton->Divide(2,2);
  cEffProton->cd(1)->DrawFrame(ptlims[0],0.,ptlims[nPtbins],1.,";#it{p}_{T} (GeV/#it{c});proton TPC PID efficiency");
  cEffProton->cd(2)->DrawFrame(ptlims[0],0.,ptlims[nPtbins],1.,";#it{p}_{T} (GeV/#it{c});proton TOF PID efficiency");
  cEffProton->cd(3)->DrawFrame(ptlims[0],0.85,ptlims[nPtbins],1.15,";#it{p}_{T} (GeV/#it{c});proton data / MC TPC PID efficiency");
  cEffProton->cd(4)->DrawFrame(ptlims[0],0.85,ptlims[nPtbins],1.15,";#it{p}_{T} (GeV/#it{c});proton data / MC TOF PID efficiency");

  for(int iEff=0; iEff<nEff; iEff++) {
    cEffPion->cd(1)->SetLogx();
    cEffPion->cd(1)->SetTopMargin(0.035);
    cEffPion->cd(1)->SetBottomMargin(0.12);
    cEffPion->cd(1)->SetRightMargin(0.035);
    hEffPionTPCMCtrue[iEff]->Draw("same");
    hEffPionTPCDataV0tag[iEff]->Draw("same");
    legEffPion->Draw("same");

    cEffPion->cd(2)->SetLogx();
    cEffPion->cd(2)->SetTopMargin(0.035);
    cEffPion->cd(2)->SetBottomMargin(0.12);
    cEffPion->cd(2)->SetRightMargin(0.035);
    hEffPionTOFMCtrue[iEff]->Draw("same");
    hEffPionTOFDataV0tag[iEff]->Draw("same");
    legEffPion->Draw("same");

    cEffPion->cd(3)->SetLogx();
    cEffPion->cd(3)->SetTopMargin(0.035);
    cEffPion->cd(3)->SetBottomMargin(0.12);
    cEffPion->cd(3)->SetRightMargin(0.035);
    hRatioEffPionTPCDataV0tag[iEff]->Draw("same");

    cEffPion->cd(4)->SetLogx();
    cEffPion->cd(4)->SetTopMargin(0.035);
    cEffPion->cd(4)->SetBottomMargin(0.12);
    cEffPion->cd(4)->SetRightMargin(0.035);
    hRatioEffPionTOFDataV0tag[iEff]->Draw("same");

    cEffKaon->cd(1)->SetLogx();
    cEffKaon->cd(1)->SetTopMargin(0.035);
    cEffKaon->cd(1)->SetBottomMargin(0.12);
    cEffKaon->cd(1)->SetRightMargin(0.035);
    hEffKaonTPCMCtrue[iEff]->Draw("same");
    hEffKaonTPCDataKinktag[iEff]->Draw("same");
    hEffKaonTPCDataTOFtag[iEff]->Draw("same");
    legEffKaonTPC->Draw("same");

    cEffKaon->cd(2)->SetLogx();
    cEffKaon->cd(2)->SetTopMargin(0.035);
    cEffKaon->cd(2)->SetBottomMargin(0.12);
    cEffKaon->cd(2)->SetRightMargin(0.035);
    hEffKaonTOFMCtrue[iEff]->Draw("same");
    hEffKaonTOFDataTPCtag[iEff]->Draw("same");
    legEffKaonTOF->Draw("same");

    cEffKaon->cd(3)->SetLogx();
    cEffKaon->cd(3)->SetTopMargin(0.035);
    cEffKaon->cd(3)->SetBottomMargin(0.12);
    cEffKaon->cd(3)->SetRightMargin(0.035);
    hRatioEffKaonTPCDataKinktag[iEff]->Draw("same");
    hRatioEffKaonTPCDataTOFtag[iEff]->Draw("same");

    cEffKaon->cd(4)->SetLogx();
    cEffKaon->cd(4)->SetTopMargin(0.035);
    cEffKaon->cd(4)->SetBottomMargin(0.12);
    cEffKaon->cd(4)->SetRightMargin(0.035);
    hRatioEffKaonTOFDataTPCtag[iEff]->Draw("same");

    cEffProton->cd(1)->SetLogx();
    cEffProton->cd(1)->SetTopMargin(0.035);
    cEffProton->cd(1)->SetBottomMargin(0.12);
    cEffProton->cd(1)->SetRightMargin(0.035);
    hEffProtonTPCMCtrue[iEff]->Draw("same");
    hEffProtonTPCDataV0tag[iEff]->Draw("same");
    legEffProton->Draw("same");

    cEffProton->cd(2)->SetLogx();
    cEffProton->cd(2)->SetTopMargin(0.035);
    cEffProton->cd(2)->SetBottomMargin(0.12);
    cEffProton->cd(2)->SetRightMargin(0.035);
    hEffProtonTOFMCtrue[iEff]->Draw("same");
    hEffProtonTOFDataV0tag[iEff]->Draw("same");
    legEffProton->Draw("same");

    cEffProton->cd(3)->SetLogx();
    cEffProton->cd(3)->SetTopMargin(0.035);
    cEffProton->cd(3)->SetBottomMargin(0.12);
    cEffProton->cd(3)->SetRightMargin(0.035);
    hRatioEffProtonTPCDataV0tag[iEff]->Draw("same");

    cEffProton->cd(4)->SetLogx();
    cEffProton->cd(4)->SetTopMargin(0.035);
    cEffProton->cd(4)->SetBottomMargin(0.12);
    cEffProton->cd(4)->SetRightMargin(0.035);
    hRatioEffProtonTOFDataV0tag[iEff]->Draw("same");
  }

  //********************************************************************************************************************************************//
  //output files
  TFile outfile("PIDSystSingleTrack.root","RECREATE");
  for(int iEff=0; iEff<nEff; iEff++) {
    hEffPionTPCMCtrue[iEff]->Write();
    hEffPionTOFMCtrue[iEff]->Write();
    hEffKaonTPCMCtrue[iEff]->Write();
    hEffKaonTOFMCtrue[iEff]->Write();
    hEffProtonTPCMCtrue[iEff]->Write();
    hEffProtonTOFMCtrue[iEff]->Write();
    hEffPionTPCDataV0tag[iEff]->Write();
    hEffPionTOFDataV0tag[iEff]->Write();
    hEffKaonTPCDataKinktag[iEff]->Write();
    hEffKaonTPCDataTOFtag[iEff]->Write();
    hEffKaonTOFDataTPCtag[iEff]->Write();
    hEffProtonTPCDataV0tag[iEff]->Write();
    hEffProtonTOFDataV0tag[iEff]->Write();
    hRatioEffPionTPCDataV0tag[iEff]->Write();
    hRatioEffPionTOFDataV0tag[iEff]->Write();
    hRatioEffKaonTPCDataKinktag[iEff]->Write();
    hRatioEffKaonTPCDataTOFtag[iEff]->Write();
    hRatioEffKaonTOFDataTPCtag[iEff]->Write();
    hRatioEffProtonTPCDataV0tag[iEff]->Write();
    hRatioEffProtonTOFDataV0tag[iEff]->Write();
  }
  outfile.Close();

  TFile outfiledist("NsigmaPIDdist.root","RECREATE");
  for(int iPt=0; iPt<nPtbins; iPt++) {
    hNsigmaTPCPionMCTrue[iPt]->Write();
    hNsigmaTPCKaonMCTrue[iPt]->Write();
    hNsigmaTPCProtonMCTrue[iPt]->Write();
    hNsigmaTOFPionMCTrue[iPt]->Write();
    hNsigmaTOFKaonMCTrue[iPt]->Write();
    hNsigmaTOFProtonMCTrue[iPt]->Write();
  }
  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      hNsigmaTPCPionMCV0tag[iPt][iPart]->Write();
      hNsigmaTPCKaonMCKinktag[iPt][iPart]->Write();
      hNsigmaTPCKaonMCTOFtag[iPt][iPart]->Write();
      hNsigmaTPCProtonMCV0tag[iPt][iPart]->Write();
      hNsigmaTOFPionMCV0tag[iPt][iPart]->Write();
      hNsigmaTOFKaonMCKinktag[iPt][iPart]->Write();
      hNsigmaTOFKaonMCTPCtag[iPt][iPart]->Write();
      hNsigmaTOFProtonMCV0tag[iPt][iPart]->Write();

    }
  }
  for(int iPt=0; iPt<nPtbins; iPt++) {
    hNsigmaTPCPionDataV0tag[iPt]->Write();
    hNsigmaTPCKaonDataKinktag[iPt]->Write();
    hNsigmaTPCKaonDataTOFtag[iPt]->Write();
    hNsigmaTPCProtonDataV0tag[iPt]->Write();
    hNsigmaTOFPionDataV0tag[iPt]->Write();
    hNsigmaTOFKaonDataKinktag[iPt]->Write();
    hNsigmaTOFKaonDataTPCtag[iPt]->Write();
    hNsigmaTOFProtonDataV0tag[iPt]->Write();
  }
  for(int iPt=0; iPt<nPtbins; iPt++) {
    hNsigmaTOFPionDataV0tag_sub[iPt]->Write();
    hNsigmaTOFKaonDataKinktag_sub[iPt]->Write();
    hNsigmaTOFKaonDataTPCtag_sub[iPt]->Write();
    hNsigmaTOFProtonDataV0tag_sub[iPt]->Write();
  }
  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      fNsigmaTPCPionMCV0tag[iPt][iPart]->Write();
      fNsigmaTPCKaonMCKinktag[iPt][iPart]->Write();
      fNsigmaTPCKaonMCTOFtag[iPt][iPart]->Write();
      fNsigmaTPCProtonMCV0tag[iPt][iPart]->Write();
    }
  }
  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      fNsigmaTPCPionDataV0tag[iPt][iPart]->Write();
      fNsigmaTPCKaonDataKinktag[iPt][iPart]->Write();
      fNsigmaTPCKaonDataTOFtag[iPt][iPart]->Write();
      fNsigmaTPCProtonDataV0tag[iPt][iPart]->Write();
    }
  }
  for(int iPt=0; iPt<nPtbins; iPt++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      if(hNsigmaTOFKaonDataTPCtag_Fit[iPt][iPart]) hNsigmaTOFKaonDataTPCtag_Fit[iPt][iPart]->Write();
      if(hNsigmaTOFProtonDataV0tag_Fit[iPt][iPart]) hNsigmaTOFProtonDataV0tag_Fit[iPt][iPart]->Write();
    }
  }
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    hFractionTPCPionMCV0tag[iPart]->Write();
    hFractionTPCKaonMCKinktag[iPart]->Write();
    hFractionTPCKaonMCTOFtag[iPart]->Write();
    hFractionTPCProtonMCV0tag[iPart]->Write();
    hFractionTOFPionMCV0tag[iPart]->Write();
    hFractionTOFKaonMCKinktag[iPart]->Write();
    hFractionTOFKaonMCTPCtag[iPart]->Write();
    hFractionTOFProtonMCV0tag[iPart]->Write();
  }
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    hFractionTOFKaonDataTPCtag[iPart]->Write();
    hFractionTOFProtonDataV0tag[iPart]->Write();
  }
  outfiledist.Close();

  
  cFractionTPCMC->SaveAs("ParticleFractionsTPCMC.pdf");
  cFractionTOFMC->SaveAs("ParticleFractionsTOFMC.pdf");
  cTOFFractionData->SaveAs("ParticleFractionsTOFdata.pdf");

  cFitResultTOFKaonFromTPCtag->SaveAs("TFractionFitterResultTOFKaonFromTPCtag.pdf");
  cFitResultTOFProtonFromV0tag->SaveAs("TFractionFitterResultTOFProtonFromV0tag.pdf");

  cPionMCV0tagTPC->SaveAs("PionMCV0tagTPC.pdf");
  cKaonMCKinkstagTPC->SaveAs("KaonMCKinkstagTPC.pdf");
  cKaonMCTOFtagTPC->SaveAs("KaonMCTOFtagTPC.pdf");
  cProtonMCV0tagTPC->SaveAs("ProtonMCV0tagTPC.pdf");
  cPionMCV0tagTOF->SaveAs("PionMCV0tagTOF.pdf");
  cKaonMCKinkstagTOF->SaveAs("KaonMCKinkstagTOF.pdf");
  cKaonMCTPCtagTOF->SaveAs("KaonMCTPCtagTOF.pdf");
  cProtonMCV0tagTOF->SaveAs("ProtonMCV0tagTOF.pdf");

  cPionDataV0tagTPC->SaveAs("PionDataV0tagTPC.pdf");
  cKaonDataKinkstagTPC->SaveAs("KaonDataKinkstagTPC.pdf");
  cKaonDataTOFtagTPC->SaveAs("KaonDataTOFtagTPC.pdf");
  cProtonDataV0tagTPC->SaveAs("ProtonDataV0tagTPC.pdf");
  cPionDataV0tagTOF->SaveAs("PionDataV0tagTOF.pdf");
  cKaonDataKinkstagTOF->SaveAs("KaonDataKinkstagTOF.pdf");
  cKaonDataTPCtagTOF->SaveAs("KaonDataTPCtagTOF.pdf");
  cProtonDataV0tagTOF->SaveAs("ProtonDataV0tagTOF.pdf");

  cEffPion->SaveAs("PionPIDefficiency_data_MC.pdf");
  cEffKaon->SaveAs("KaonPIDefficiency_data_MC.pdf");
  cEffProton->SaveAs("ProtonPIDefficiency_data_MC.pdf");

  return 0;
}

//_____________________________________________________
int GetHistoParticleIndex(short pdgcode) {
  switch(pdgcode) {
    case 11:
      return kElectron;
    break;
    case 13:
      return kMuon;
    break;
    case 211:
      return kPion;
    break;
    case 321:
      return kKaon;
    break;
    case 2212:
      return kProton;
    break;
  }
  return -1;
}

//_____________________________________________________
int FindPtbin(float pt, const double ptlims[], int nPtbins) {

  int ptbin=-1;
  for(int iPt=0; iPt<nPtbins; iPt++) {
    if(pt>ptlims[iPt] && pt<ptlims[iPt+1]) {
      ptbin=iPt;
      break;
    }
  }

  return ptbin;
}

//______________________________________________________
void ComputeEfficiency(double num, double den, double &eff, double &effunc) {
  
  TH1F htmpnum("htmpnum","",1,0,1);
  TH1F htmpden("htmpden","",1,0,1);
  TH1F htmpratio("htmpratio","",1,0,1);

  htmpnum.SetBinContent(1,num);
  htmpden.SetBinContent(1,den);
  htmpnum.SetBinError(1,TMath::Sqrt(num));
  htmpden.SetBinError(1,TMath::Sqrt(den));
  htmpratio.Divide(&htmpnum,&htmpden,1.,1.,"B");
  eff=htmpratio.GetBinContent(1);
  effunc=htmpratio.GetBinError(1);
}

//______________________________________________________
void GetTOFFractionsFromData(int whichpart, int iPt, TH1F* hFractionMC[nPDGcodes-1], TH1F* hFractionData[nPDGcodes-1], TH1F* hNsigmaMC[nPDGcodes], TH1F* hNsigmaData, TFractionFitter *&fNsigmaFitter, vector<int> &templUsed) {

  TObjArray* oNsigmaMC = new TObjArray(0);
    
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      if(iPart==whichpart || hFractionMC[iPart]->GetBinContent(iPt+1)>1.e-5) {
        TH1F* hMC = (TH1F*)hNsigmaMC[iPart]->Clone(Form("hMC_%d_%s",iPt,pdgnames[iPart].Data()));
        oNsigmaMC->AddLast(hMC);
        templUsed.push_back(iPart);
      }
    }
    
    if(templUsed.size()==1) {
      for(int iPart=kElectron; iPart<=kProton; iPart++) {
        if(iPart==whichpart) {
          hFractionData[iPart]->SetBinContent(iPt+1,1.);
        }
       else {
          hFractionData[iPart]->SetBinContent(iPt+1,0.);
        }
      }
      return;
    }

    if(oNsigmaMC->GetEntries()>1) {
      
      fNsigmaFitter = new TFractionFitter(hNsigmaData,oNsigmaMC);
      for(int iEntry=0; iEntry<oNsigmaMC->GetEntries(); iEntry++) {
        fNsigmaFitter->Constrain(iEntry,hFractionMC[templUsed[iEntry]]->GetBinContent(iPt+1)*0.5,hFractionMC[templUsed[iEntry]]->GetBinContent(iPt+1)*2);
      }
      fNsigmaFitter->Fit();
      
      for(int iPart=kElectron; iPart<=kProton; iPart++) {

        vector<int>::iterator it = find(templUsed.begin(), templUsed.end(), iPart);
        if(it!=templUsed.end()) {
          double frac, err;
          int iEntry = static_cast<int>(distance(templUsed.begin(),it));

          fNsigmaFitter->GetResult(iEntry, frac, err);
          hFractionData[iPart]->SetBinContent(iPt+1,frac);
          hFractionData[iPart]->SetBinError(iPt+1,err);
        }
        else {
          hFractionData[iPart]->SetBinContent(iPt+1,0);
          hFractionData[iPart]->SetBinError(iPt+1,0);
        }
      }
    }
    else {
      for(int iPart=kElectron; iPart<=kProton; iPart++) {
        if(iPart!=whichpart) {
          hFractionData[iPart]->SetBinContent(iPt+1,1);
          hFractionData[iPart]->SetBinError(iPt+1,1);
        }
        else {
          hFractionData[iPart]->SetBinContent(iPt+1,0);
          hFractionData[iPart]->SetBinError(iPt+1,0);
        }
      }
    }
  
}

//______________________________________________________
double PDFnsigmaTPCtot(double* nsigma, double* pars) {
  
  double partpdf[nPDGcodes-1];
  double totalpdf = 0;
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    partpdf[iPart] = pars[3*iPart]*TMath::Gaus(nsigma[0],pars[3*iPart+1],pars[3*iPart+2]);
    totalpdf += partpdf[iPart];
  }
  
  return totalpdf;
}

//_____________________________________________________
void PlotQAhistos(TList* listMC, TList* listData) {
  
  //V0 QA histos
  TString ArmenterosName[5] = {"All","K0s","Lambda","AntiLambda","Gamma"};
  TString legV0name[4] = {"K_{s}^{0} #rightarrow #pi^{+}#pi^{-}","#Lambda #rightarrow p#pi^{-}","#bar{#Lambda} #rightarrow #bar{p}#pi^{+}","#gamma #rightarrow e^{+}e^{-}"};
  int ArmenterosColor[5] = {kBlack,kOrange+7,kRed,kBlue,kGreen+3};
  TH2F* hArmenterosMC[5];
  TH2F* hArmenterosData[5];

  //kinks QA histos
  TString KinksName[5] = {"QtVsMassKinks","PDaughterVsMotherKink","OpeningAngleVsPMotherKink","dEdxVsPMotherKink","NTPCclsVsRadius"};
  TH2F* hKinksMC[5];
  TH2F* hKinksData[5];

  TCanvas* cArmentero = new TCanvas("cArmenteros","cArmenteros",1920,1080);
  TCanvas* cKinks = new TCanvas("cKinks","cKinks",1920,1080);
  cArmentero->Divide(2,1);
  cKinks->Divide(5,2);
  TLatex* latV0[5];

  for(int iHisto=0; iHisto<5; iHisto++) {
    hArmenterosMC[iHisto] = (TH2F*)listMC->FindObject(Form("fHistArmenteroPlot%s",ArmenterosName[iHisto].Data()));
    hArmenterosData[iHisto] = (TH2F*)listData->FindObject(Form("fHistArmenteroPlot%s",ArmenterosName[iHisto].Data()));
    hArmenterosMC[iHisto]->SetDirectory(0);
    hArmenterosData[iHisto]->SetDirectory(0);
    hArmenterosData[iHisto]->SetMarkerColor(ArmenterosColor[iHisto]);
    hArmenterosMC[iHisto]->SetMarkerColor(ArmenterosColor[iHisto]);
    hArmenterosData[iHisto]->GetYaxis()->SetTitleOffset(1.4);
    hArmenterosMC[iHisto]->GetYaxis()->SetTitleOffset(1.4);
    latV0[iHisto] = new TLatex();
    latV0[iHisto]->SetNDC();
    latV0[iHisto]->SetTextFont(42);
    latV0[iHisto]->SetTextColor(ArmenterosColor[iHisto]);
    latV0[iHisto]->SetTextSize(0.045);
    cArmentero->cd(1);
    hArmenterosMC[iHisto]->Draw("same");
    cArmentero->cd(2);
    hArmenterosData[iHisto]->Draw("same");
    if(iHisto!=0) {
      cArmentero->cd(1);
      latV0[iHisto]->DrawLatex(0.6,0.8-0.05*iHisto,legV0name[iHisto-1].Data());
      cArmentero->cd(2);
      latV0[iHisto]->DrawLatex(0.6,0.8-0.05*iHisto,legV0name[iHisto-1].Data());
    }

    hKinksMC[iHisto] = (TH2F*)listMC->FindObject(Form("fHist%s",KinksName[iHisto].Data()));
    hKinksData[iHisto] = (TH2F*)listData->FindObject(Form("fHist%s",KinksName[iHisto].Data()));
    hKinksMC[iHisto]->SetTitle("MC");
    hKinksData[iHisto]->SetTitle("Data");

    cKinks->cd(1+iHisto)->SetLogz();
    hKinksMC[iHisto]->Draw("colz");
    cKinks->cd(6+iHisto)->SetLogz();
    hKinksData[iHisto]->Draw("colz");
  }

  TH1F* hKinksMassMC = (TH1F*)hKinksMC[0]->ProjectionX("hKinksMassMC");
  hKinksMassMC->Scale(1./hKinksMassMC->Integral());
  hKinksMassMC->SetLineColor(kBlue);
  hKinksMassMC->SetLineWidth(2);
  hKinksMassMC->SetMarkerColor(kBlue);
  hKinksMassMC->SetMarkerStyle(kFullSquare);
  hKinksMassMC->SetTitle("");
  hKinksMassMC->GetYaxis()->SetTitle("Normalised entries");
  hKinksMassMC->GetYaxis()->SetTitleSize(0.045);
  hKinksMassMC->GetYaxis()->SetLabelSize(0.04);
  hKinksMassMC->GetXaxis()->SetTitleSize(0.045);
  hKinksMassMC->GetXaxis()->SetLabelSize(0.04);
  
  TH1F* hKinksMassData = (TH1F*)hKinksData[0]->ProjectionX("hKinksMassData");
  hKinksMassData->Scale(1./hKinksMassData->Integral());
  hKinksMassData->SetLineColor(kRed);
  hKinksMassData->SetLineWidth(2);
  hKinksMassData->SetMarkerColor(kRed);
  hKinksMassData->SetMarkerStyle(kFullDiamond);
  hKinksMassData->SetTitle("");
  hKinksMassData->GetYaxis()->SetTitle("Normalised entries");
  hKinksMassData->GetYaxis()->SetTitleSize(0.045);
  hKinksMassData->GetYaxis()->SetLabelSize(0.04);
  hKinksMassData->GetXaxis()->SetTitleSize(0.045);
  hKinksMassData->GetXaxis()->SetLabelSize(0.04);

  TLegend* leg = new TLegend(0.2,0.7,0.4,0.8);
  leg->SetTextSize(0.045);
  leg->AddEntry(hKinksMassMC,"MC","lpe");
  leg->AddEntry(hKinksMassData,"Data","lpe");
  TCanvas* cKinksMass = new TCanvas("cKinksMass","",800,800);
  cKinksMass->SetTopMargin(0.1);
  hKinksMassMC->Draw();
  hKinksMassData->Draw("same");
  leg->Draw("same");

  cArmentero->SaveAs("V0QAplots.pdf");
  cKinks->SaveAs("KinksQAplots.pdf");
  cKinksMass->SaveAs("KinksMass.pdf");
  cArmentero->SaveAs("V0QAplots.png");
  cKinks->SaveAs("KinksQAplots.png");
}

//_____________________________________________________
void DivideCanvas(TCanvas* c, int nPtbins) {
  if(nPtbins<2)
    c->cd();
  else if(nPtbins==2 || nPtbins==3)
    c->Divide(nPtbins,1);
  else if(nPtbins==4 || nPtbins==6 || nPtbins==8)
    c->Divide(nPtbins/2,2);
  else if(nPtbins==5 || nPtbins==7)
    c->Divide((nPtbins+1)/2,2);
  else if(nPtbins==9 || nPtbins==12 || nPtbins==15)
    c->Divide(nPtbins/3,3);
  else if(nPtbins==10 || nPtbins==11)
    c->Divide(4,3);
  else if(nPtbins==13 || nPtbins==14)
    c->Divide(5,3);
  else if(nPtbins>15 && nPtbins<=20 && nPtbins%4==0)
    c->Divide(nPtbins/4,4);
  else if(nPtbins>15 && nPtbins<=20 && nPtbins%4!=0)
    c->Divide(5,4);
  else if(nPtbins==21)
    c->Divide(7,3);
  else if(nPtbins>21 && nPtbins<=25)
    c->Divide(5,5);
  else if(nPtbins>25 && nPtbins%2==0)
    c->Divide(nPtbins/2,2);
  else
    c->Divide((nPtbins+1)/2,2);
}

//_____________________________________________________
void SetStyle() {
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetPadTopMargin(0.14);
  gStyle->SetPadRightMargin(0.035);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetTitleSize(0.05,"xyzt");
  gStyle->SetLabelSize(0.045,"xyz");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  
  TGaxis::SetMaxDigits(3);
}

//_____________________________________________________
void SetTH1Style(TH1F* histo, int markerstyle, int markercolor, float markersize, int linewidth, int linecolor, int fillcolor, float labelsize, float titlesize) {

  histo->SetMarkerStyle(markerstyle);
  histo->SetMarkerSize(markersize);
  histo->SetMarkerColor(markercolor);
  histo->SetLineWidth(linewidth);
  histo->SetLineColor(linecolor);
  histo->SetFillColorAlpha(fillcolor,0.25);
  
  if(labelsize>0) {
    histo->GetXaxis()->SetLabelSize(labelsize);
    histo->GetYaxis()->SetLabelSize(labelsize);
  }
  if(titlesize>0) {
    histo->GetXaxis()->SetTitleSize(labelsize);
    histo->GetYaxis()->SetTitleSize(labelsize);
    histo->SetTitleSize(labelsize);
  }
}
