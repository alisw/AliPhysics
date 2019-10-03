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
#include <TProfile.h>
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
enum vars {kPt,kP};
enum partPID {kElectron,kMuon,kPion,kKaon,kProton,kAll};
const int nPDGcodes = 6;
const int pdgcodes[nPDGcodes]={11,13,211,321,2212,-100};
const int pdgcolors[nPDGcodes]={kOrange+7,kGray,kRed,kBlue,kGreen+2,kBlack};
const int pdgfillcolors[nPDGcodes]={kOrange+7,kGray,kRed,kBlue,kGreen+2,kWhite};
const TString pdgnames[nPDGcodes]={"Electron","Muon","Pion","Kaon","Proton","All"};
//TO EDIT
const double binlims[] = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};
const TString infileNameData = "/AnalysisResults_data.root"; 
const TString indirNameData = "PWGHF_D2H_SystNsigmaPID"; 
const TString inlistNameData = "coutputPIDhistos"; 
const TString infileNameMC = "AnalysisResults_MC.root"; 
const TString indirNameMC = "PWGHF_D2H_SystNsigmaPID"; 
const TString inlistNameMC = "coutputPIDhistos"; 
const TString outputdirName = "outputs/";
//applied only to data
double absetamin = 0.;
double absetamax = 1.;

//_____________________________________________________
//METHOD PROTOTYPES
int EstimateSingleTrackPIDsyst(int maxEntries=1e9, int var4proj=kP);
int GetHistoParticleIndex(short pdgcode);
int FindPtbin(float pt, const double binlims[], int nBins);
void ComputeEfficiency(double num, double den, double &eff, double &effunc);
void GetTOFFractionsFromData(int whichpart, int iBin, TH1F* hFractionMC[nPDGcodes-1], TH1F* hFractionData[nPDGcodes-1], TH1F* hNsigmaMC[nPDGcodes], TH1F* hNsigmaData, TFractionFitter *&fNsigmaFitter, vector<int> &templUsed);
double PDFnsigmaTPCtot(double* nsigma, double* pars);
void PlotQAhistos(TList* listMC, TList* listData);
void DivideCanvas(TCanvas* c, int nBins);
void SetStyle();
void SetTH1Style(TH1F* histo, int markerstyle, int markercolor, float markersize, int linewidth, int linecolor, int fillcolor, float labelsize=-1, float titlesize=-1);

//_____________________________________________________
//METHOD IMPLEMENTATIONS
int EstimateSingleTrackPIDsyst(int maxEntries, int var4proj) {

  SetStyle();
  
  TString varname = "";
  if(var4proj==kPt)
    varname = "#it{p}_{T}";
  else if(var4proj==kP)
    varname = "#it{p}";
  else {
    cerr << "You can chose to project the Nsigma distributions vs. pT or p only. Exit" << endl;
    return -1;
  }

  //********************************************************************************************************************************************//
  //get pt bins
  const int nBins = sizeof(binlims)/sizeof(binlims[0])-1;

  cout << "\n\n********************\n" << endl;
  cout << nBins << " bins (GeV/c):" << endl;
  for(int iBin=0; iBin<nBins; iBin++) {
    cout << binlims[iBin]<<"-"<< binlims[iBin+1] <<endl;
  }
  
  //********************************************************************************************************************************************//
  //define histos

  //MC truth
  TH2F *hNsigmaTPCPionVsPtMCTrue, *hNsigmaTPCKaonVsPtMCTrue, *hNsigmaTPCProtonVsPtMCTrue;
  TH2F *hNsigmaTOFPionVsPtMCTrue, *hNsigmaTOFKaonVsPtMCTrue, *hNsigmaTOFProtonVsPtMCTrue;

  TH1F *hNsigmaTPCPionMCTrue[nBins], *hNsigmaTPCKaonMCTrue[nBins], *hNsigmaTPCProtonMCTrue[nBins];
  TH1F *hNsigmaTOFPionMCTrue[nBins], *hNsigmaTOFKaonMCTrue[nBins], *hNsigmaTOFProtonMCTrue[nBins];

  //MC tagged
  TH1F *hNsigmaTPCPionMCV0tag[nBins][nPDGcodes], *hNsigmaTPCKaonMCKinktag[nBins][nPDGcodes], *hNsigmaTPCKaonMCTOFtag[nBins][nPDGcodes], *hNsigmaTPCProtonMCV0tag[nBins][nPDGcodes];
  TH1F *hNsigmaTOFPionMCV0tag[nBins][nPDGcodes], *hNsigmaTOFKaonMCKinktag[nBins][nPDGcodes], *hNsigmaTOFKaonMCTPCtag[nBins][nPDGcodes], *hNsigmaTOFProtonMCV0tag[nBins][nPDGcodes];
  
  //data tagged
  TH1F *hNsigmaTPCPionDataV0tag[nBins], *hNsigmaTPCKaonDataKinktag[nBins], *hNsigmaTPCKaonDataTOFtag[nBins], *hNsigmaTPCProtonDataV0tag[nBins];
  TH1F *hNsigmaTOFPionDataV0tag[nBins], *hNsigmaTOFKaonDataKinktag[nBins], *hNsigmaTOFKaonDataTPCtag[nBins], *hNsigmaTOFProtonDataV0tag[nBins];

  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=0; iPart<nPDGcodes; iPart++) {
      hNsigmaTPCPionMCV0tag[iBin][iPart] = new TH1F(Form("hNsigmaTPCPionMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(#pi);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      hNsigmaTPCKaonMCKinktag[iBin][iPart] = new TH1F(Form("hNsigmaTPCKaonMCKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      hNsigmaTPCKaonMCTOFtag[iBin][iPart] = new TH1F(Form("hNsigmaTPCKaonMCTOFtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      hNsigmaTPCProtonMCV0tag[iBin][iPart] = new TH1F(Form("hNsigmaTPCProtonMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      
      hNsigmaTOFPionMCV0tag[iBin][iPart] = new TH1F(Form("hNsigmaTOFPionMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(#pi);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      hNsigmaTOFKaonMCKinktag[iBin][iPart] = new TH1F(Form("hNsigmaTOFKaonMCKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      hNsigmaTOFKaonMCTPCtag[iBin][iPart] = new TH1F(Form("hNsigmaTOFKaonMCTPCtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      hNsigmaTOFProtonMCV0tag[iBin][iPart] = new TH1F(Form("hNsigmaTOFProtonMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(p);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
      
      SetTH1Style(hNsigmaTPCPionMCV0tag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTPCKaonMCKinktag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTPCKaonMCTOFtag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTPCProtonMCV0tag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFPionMCV0tag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFKaonMCKinktag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFKaonMCTPCtag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
      SetTH1Style(hNsigmaTOFProtonMCV0tag[iBin][iPart],kFullCircle,pdgcolors[iPart],0.6,2,pdgcolors[iPart],pdgfillcolors[iPart],0.055,0.06);
    }
    
    hNsigmaTPCPionDataV0tag[iBin] = new TH1F(Form("hNsigmaTPCPionDataV0tag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(#pi);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    hNsigmaTPCKaonDataKinktag[iBin] = new TH1F(Form("hNsigmaTPCKaonDataKinktag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    hNsigmaTPCKaonDataTOFtag[iBin] = new TH1F(Form("hNsigmaTPCKaonDataTOFtag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    hNsigmaTPCProtonDataV0tag[iBin] = new TH1F(Form("hNsigmaTPCProtonDataV0tag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TPC}(p);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    
    hNsigmaTOFPionDataV0tag[iBin] = new TH1F(Form("hNsigmaTOFPionDataV0tag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(#pi);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    hNsigmaTOFKaonDataKinktag[iBin] = new TH1F(Form("hNsigmaTOFKaonDataKinktag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    hNsigmaTOFKaonDataTPCtag[iBin] = new TH1F(Form("hNsigmaTOFKaonDataTPCtag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(K);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    hNsigmaTOFProtonDataV0tag[iBin] = new TH1F(Form("hNsigmaTOFProtonDataV0tag_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),Form("%0.2f < %s < %0.2f GeV/#it{c};N_{#sigma}^{TOF}(p);Normalised entries",binlims[iBin],varname.Data(),binlims[iBin+1]),1000,-50,50);
    
    SetTH1Style(hNsigmaTPCPionDataV0tag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCKaonDataKinktag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCKaonDataTOFtag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCProtonDataV0tag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFPionDataV0tag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataKinktag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataTPCtag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFProtonDataV0tag[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
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

  for(int iBin=0; iBin<nBins; iBin++) {
    int ptbinmin = hNsigmaTPCPionVsPtMCTrue->GetXaxis()->FindBin(binlims[iBin]*1.0001);
    int ptbinmax = hNsigmaTPCPionVsPtMCTrue->GetXaxis()->FindBin(binlims[iBin+1]*0.9999);

    hNsigmaTPCPionMCTrue[iBin] = (TH1F*)hNsigmaTPCPionVsPtMCTrue->ProjectionY(Form("hNsigmaTPCPionMCTrue_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),ptbinmin,ptbinmax);
    hNsigmaTPCKaonMCTrue[iBin] = (TH1F*)hNsigmaTPCKaonVsPtMCTrue->ProjectionY(Form("hNsigmaTPCKaonMCTrue_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),ptbinmin,ptbinmax);
    hNsigmaTPCProtonMCTrue[iBin] = (TH1F*)hNsigmaTPCProtonVsPtMCTrue->ProjectionY(Form("hNsigmaTPCProtonMCTrue_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),ptbinmin,ptbinmax);

    hNsigmaTOFPionMCTrue[iBin] = (TH1F*)hNsigmaTOFPionVsPtMCTrue->ProjectionY(Form("hNsigmaTOFPionMCTrue_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),ptbinmin,ptbinmax);
    hNsigmaTOFKaonMCTrue[iBin] = (TH1F*)hNsigmaTOFKaonVsPtMCTrue->ProjectionY(Form("hNsigmaTOFKaonMCTrue_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),ptbinmin,ptbinmax);
    hNsigmaTOFProtonMCTrue[iBin] = (TH1F*)hNsigmaTOFProtonVsPtMCTrue->ProjectionY(Form("hNsigmaTOFProtonMCTrue_p%0.f_%0.f",binlims[iBin]*10,binlims[iBin+1]*10),ptbinmin,ptbinmax);

    SetTH1Style(hNsigmaTPCPionMCTrue[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCKaonMCTrue[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTPCProtonMCTrue[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFPionMCTrue[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonMCTrue[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
    SetTH1Style(hNsigmaTOFProtonMCTrue[iBin],kFullCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],kWhite,0.055,0.06);
  }
  
  short pdgMC             = -1;
  unsigned short tagMC    = 0;
  short n_sigma_TPC_pi_MC = -999;
  short n_sigma_TPC_K_MC  = -999;
  short n_sigma_TPC_p_MC  = -999;
  short n_sigma_TOF_pi_MC = -999;
  short n_sigma_TOF_K_MC  = -999;
  short n_sigma_TOF_p_MC  = -999;
  unsigned short pT_MC    = 0;
  unsigned short pTPC_MC  = 0;
  unsigned short pTOF_MC  = 0;

  treePIDMC->SetBranchAddress("PDGcode",&pdgMC);
  treePIDMC->SetBranchAddress("tag",&tagMC);
  treePIDMC->SetBranchAddress("n_sigma_TPC_pi",&n_sigma_TPC_pi_MC);
  treePIDMC->SetBranchAddress("n_sigma_TPC_K",&n_sigma_TPC_K_MC);
  treePIDMC->SetBranchAddress("n_sigma_TPC_p",&n_sigma_TPC_p_MC);
  treePIDMC->SetBranchAddress("n_sigma_TOF_pi",&n_sigma_TOF_pi_MC);
  treePIDMC->SetBranchAddress("n_sigma_TOF_K",&n_sigma_TOF_K_MC);
  treePIDMC->SetBranchAddress("n_sigma_TOF_p",&n_sigma_TOF_p_MC);
  treePIDMC->SetBranchAddress("pT",&pT_MC);
  treePIDMC->SetBranchAddress("pTPC",&pTPC_MC);
  treePIDMC->SetBranchAddress("pTOF",&pTOF_MC);

  cout << "\n********** Loop on MC tracks **********\n" << endl;

  for(int iEntry=0; iEntry<treePIDMC->GetEntriesFast(); iEntry++) {
    
    if(iEntry>maxEntries) break;

    if(iEntry%1000000==0 || iEntry==treePIDMC->GetEntriesFast()-1) cout << Form("MC Track %010d",iEntry) << endl;

    treePIDMC->GetEntry(iEntry);
    int iBin = -1;
    if(var4proj==kPt) 
      iBin = FindPtbin(static_cast<float>(pT_MC)/1000,binlims,nBins);
    else if(var4proj==kP)
      iBin = FindPtbin(static_cast<float>(pTPC_MC)/1000,binlims,nBins);
    if(iBin<0 || iBin>=static_cast<int>(nBins)) continue;

    if((tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromK0s || tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromL)) {
      hNsigmaTPCPionMCV0tag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TPC_pi_MC)/100);
      hNsigmaTOFPionMCV0tag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TOF_pi_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks) {
      hNsigmaTPCKaonMCKinktag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
      hNsigmaTOFKaonMCKinktag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF) {
      hNsigmaTPCKaonMCTOFtag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC) {
      hNsigmaTOFKaonMCTPCtag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
    }
    if(tagMC&AliAnalysisTaskSEHFSystPID::kIsProtonFromL) {
      hNsigmaTPCProtonMCV0tag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TPC_p_MC)/100);
      hNsigmaTOFProtonMCV0tag[iBin][kAll]->Fill(static_cast<float>(n_sigma_TOF_p_MC)/100);
    }

    int iHisto = GetHistoParticleIndex(pdgMC);
    if(iHisto>=kElectron && iHisto<kAll) {
      if((tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromK0s || tagMC&AliAnalysisTaskSEHFSystPID::kIsPionFromL)) {
        hNsigmaTPCPionMCV0tag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TPC_pi_MC)/100);
        hNsigmaTOFPionMCV0tag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TOF_pi_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks) {
        hNsigmaTPCKaonMCKinktag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
        hNsigmaTOFKaonMCKinktag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF) {
        hNsigmaTPCKaonMCTOFtag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TPC_K_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC) {
        hNsigmaTOFKaonMCTPCtag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TOF_K_MC)/100);
      }
      if(tagMC&AliAnalysisTaskSEHFSystPID::kIsProtonFromL) {
        hNsigmaTPCProtonMCV0tag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TPC_p_MC)/100);
        hNsigmaTOFProtonMCV0tag[iBin][iHisto]->Fill(static_cast<float>(n_sigma_TOF_p_MC)/100);
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

  unsigned short tagData    = 0;
  short n_sigma_TPC_pi_Data = -999;
  short n_sigma_TPC_K_Data  = -999;
  short n_sigma_TPC_p_Data  = -999;
  short n_sigma_TOF_pi_Data = -999;
  short n_sigma_TOF_K_Data  = -999;
  short n_sigma_TOF_p_Data  = -999;
  unsigned short pT_Data    = 0;
  unsigned short pTPC_Data  = 0;
  unsigned short pTOF_Data  = 0;
  short eta_Data            = 9;

  treePIDData->SetBranchAddress("tag",&tagData);
  treePIDData->SetBranchAddress("n_sigma_TPC_pi",&n_sigma_TPC_pi_Data);
  treePIDData->SetBranchAddress("n_sigma_TPC_K",&n_sigma_TPC_K_Data);
  treePIDData->SetBranchAddress("n_sigma_TPC_p",&n_sigma_TPC_p_Data);
  treePIDData->SetBranchAddress("n_sigma_TOF_pi",&n_sigma_TOF_pi_Data);
  treePIDData->SetBranchAddress("n_sigma_TOF_K",&n_sigma_TOF_K_Data);
  treePIDData->SetBranchAddress("n_sigma_TOF_p",&n_sigma_TOF_p_Data);
  treePIDData->SetBranchAddress("pT",&pT_Data);
  treePIDData->SetBranchAddress("pTPC",&pTPC_Data);
  treePIDData->SetBranchAddress("pTOF",&pTOF_Data);
  treePIDData->SetBranchAddress("eta",&eta_Data);

  cout << "\n********** Loop on data tracks **********\n" << endl;

  TH2F* hNsigmaTPCKaonTPCtagged = new TH2F("hNsigmaTPCKaonTPCtagged",Form(";%s (GeV/#it{c});N_{#sigma}(K)",varname.Data()),100,0.,10.,1000,-50,50);
  hNsigmaTPCKaonTPCtagged->SetMarkerColor(kBlue);
  TH2F* hNsigmaTPCKaon = new TH2F("hNsigmaTPCKaon",Form(";%s (GeV/#it{c});N_{#sigma}^{TPC}(K)",varname.Data()),100,0.,10.,1000,-50,50);
  hNsigmaTPCKaon->SetMarkerColor(kBlack);
  TH2F* hNsigmaTOFKaonTOFtagged = new TH2F("hNsigmaTOFKaonTOFtagged",Form(";%s (GeV/#it{c});N_{#sigma}(K)",varname.Data()),100,0.,10.,1000,-50,50);
  hNsigmaTOFKaonTOFtagged->SetMarkerColor(kBlue);
  TH2F* hNsigmaTOFKaon = new TH2F("hNsigmaTOFKaon",Form(";%s (GeV/#it{c});N_{#sigma}^{TOF}(K)",varname.Data()),100,0.,10.,1000,-50,50);
  hNsigmaTOFKaon->SetMarkerColor(kBlack);
  
  TH2F* hNsigmaTPCPionDataV0tagVsEta = new TH2F("hNsigmaTPCPionDataV0tagVsEta",Form("#pi V0 tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TPC}(#pi)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
  TH2F* hNsigmaTPCKaonDataKinktagVsEta = new TH2F("hNsigmaTPCKaonDataKinktagVsEta",Form("K kinks V0 tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TPC}(K)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
  TH2F* hNsigmaTPCKaonDataTOFtagVsEta = new TH2F("hNsigmaTPCKaonDataTOFtagVsEta",Form("K TOF tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TPC}(K)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
  TH2F* hNsigmaTPCProtonDataV0tagVsEta = new TH2F("hNsigmaTPCProtonDataV0tagVsEta",Form("p V0 tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TPC}(p)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
    
  TH2F* hNsigmaTOFPionDataV0tagVsEta = new TH2F("hNsigmaTOFPionDataV0tagVsEta",Form("#pi V0 tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TOF}(#pi)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
  TH2F* hNsigmaTOFKaonDataKinktagVsEta = new TH2F("hNsigmaTOFKaonDataKinktagVsEta",Form("K kinks tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TOF}(K)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
  TH2F* hNsigmaTOFKaonDataTPCtagVsEta = new TH2F("hNsigmaTOFKaonDataTPCtagVsEta",Form("K TPC tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TOF}(K)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);
  TH2F* hNsigmaTOFProtonDataV0tagVsEta = new TH2F("hNsigmaTOFProtonDataV0tagVsEta",Form("p V0 tag %0.2f < %s < %0.2f GeV/#it{c};#eta;N_{#sigma}^{TOF}(p)",binlims[0],varname.Data(),binlims[nBins]),100,-1,1,100,-5,5);

  for(int iEntry=0; iEntry<treePIDData->GetEntriesFast(); iEntry++) {

    if(iEntry>maxEntries) break;

    if(iEntry%1000000==0 || iEntry==treePIDData->GetEntriesFast()-1) cout << Form("Data Track %010d",iEntry) << endl;

    treePIDData->GetEntry(iEntry);
    float etaf = static_cast<float>(eta_Data)/1000;
    if(TMath::Abs(etaf)<absetamin || TMath::Abs(etaf)>absetamax) continue;
    int iBin = -1;
    float xvar = -1.;
    if(var4proj==kPt) {
      xvar = static_cast<float>(pT_Data)/1000;
      iBin = FindPtbin(xvar,binlims,nBins);
    }
    else if(var4proj==kP) {
      xvar = static_cast<float>(pTPC_Data)/1000;
      iBin = FindPtbin(xvar,binlims,nBins);
    }
    if(iBin<0 || iBin>=static_cast<int>(nBins)) continue;

    hNsigmaTPCKaon->Fill(xvar,static_cast<float>(n_sigma_TPC_K_Data)/100);
    hNsigmaTOFKaon->Fill(xvar,static_cast<float>(n_sigma_TOF_K_Data)/100);
    if((tagData&AliAnalysisTaskSEHFSystPID::kIsPionFromK0s || tagData&AliAnalysisTaskSEHFSystPID::kIsPionFromL)) {
      hNsigmaTPCPionDataV0tag[iBin]->Fill(static_cast<float>(n_sigma_TPC_pi_Data)/100);
      hNsigmaTOFPionDataV0tag[iBin]->Fill(static_cast<float>(n_sigma_TOF_pi_Data)/100);
      hNsigmaTPCPionDataV0tagVsEta->Fill(etaf,static_cast<float>(n_sigma_TPC_pi_Data)/100);
      hNsigmaTOFPionDataV0tagVsEta->Fill(etaf,static_cast<float>(n_sigma_TOF_pi_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsKaonFromKinks) {
      hNsigmaTPCKaonDataKinktag[iBin]->Fill(static_cast<float>(n_sigma_TPC_K_Data)/100);
      hNsigmaTOFKaonDataKinktag[iBin]->Fill(static_cast<float>(n_sigma_TOF_K_Data)/100);
      hNsigmaTPCKaonDataKinktagVsEta->Fill(etaf,static_cast<float>(n_sigma_TPC_K_Data)/100);
      hNsigmaTOFKaonDataKinktagVsEta->Fill(etaf,static_cast<float>(n_sigma_TOF_K_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsKaonFromTOF) {
      hNsigmaTPCKaonDataTOFtag[iBin]->Fill(static_cast<float>(n_sigma_TPC_K_Data)/100);
      hNsigmaTOFKaonTOFtagged->Fill(xvar,static_cast<float>(n_sigma_TOF_K_Data)/100);
      hNsigmaTPCKaonDataTOFtagVsEta->Fill(etaf,static_cast<float>(n_sigma_TPC_K_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsKaonFromTPC) {
      hNsigmaTOFKaonDataTPCtag[iBin]->Fill(static_cast<float>(n_sigma_TOF_K_Data)/100);
      hNsigmaTPCKaonTPCtagged->Fill(xvar,static_cast<float>(n_sigma_TPC_K_Data)/100);
      hNsigmaTOFKaonDataTPCtagVsEta->Fill(etaf,static_cast<float>(n_sigma_TOF_K_Data)/100);
    }
    if(tagData&AliAnalysisTaskSEHFSystPID::kIsProtonFromL) {
      hNsigmaTPCProtonDataV0tag[iBin]->Fill(static_cast<float>(n_sigma_TPC_p_Data)/100);
      hNsigmaTOFProtonDataV0tag[iBin]->Fill(static_cast<float>(n_sigma_TOF_p_Data)/100);
      hNsigmaTPCProtonDataV0tagVsEta->Fill(etaf,static_cast<float>(n_sigma_TPC_p_Data)/100);
      hNsigmaTOFProtonDataV0tagVsEta->Fill(etaf,static_cast<float>(n_sigma_TOF_p_Data)/100);
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

  TCanvas* cNsigmaTag = new TCanvas("cNsigmaTag","",1920,1080);
  cNsigmaTag->Divide(2,1);
  cNsigmaTag->cd(1);
  hNsigmaTPCKaon->SetMarkerStyle(kFullCircle);
  hNsigmaTPCKaonTPCtagged->SetMarkerStyle(kFullCircle);
  hNsigmaTPCKaon->SetMarkerSize(0.5);
  hNsigmaTPCKaonTPCtagged->SetMarkerSize(0.5);
  hNsigmaTPCKaon->DrawCopy("");
  hNsigmaTPCKaonTPCtagged->DrawCopy("same");
  latAll->DrawLatex(0.65,0.25,"All tags");
  latTag->DrawLatex(0.65,0.2,"TPC tag");
  cNsigmaTag->cd(2);
  hNsigmaTOFKaon->SetMarkerStyle(kFullCircle);
  hNsigmaTOFKaonTOFtagged->SetMarkerStyle(kFullCircle);
  hNsigmaTOFKaon->SetMarkerSize(0.5);
  hNsigmaTOFKaonTOFtagged->SetMarkerSize(0.5);
  hNsigmaTOFKaon->DrawCopy("");
  hNsigmaTOFKaonTOFtagged->DrawCopy("same");
  latAll->DrawLatex(0.65,0.25,"All tags");
  latTag->DrawLatex(0.65,0.2,"TOF tag");
  
  cNsigmaTag->SaveAs(Form("%s/TPC_TOF_tag.pdf",outputdirName.Data()));
  cNsigmaTag->SaveAs(Form("%s/TPC_TOF_tag.png",outputdirName.Data()));

  TProfile* hProfileNsigmaTPCPionDataV0tagVsEta = hNsigmaTPCPionDataV0tagVsEta->ProfileX("hProfileNsigmaTPCPionDataV0tagVsEta");
  hProfileNsigmaTPCPionDataV0tagVsEta->SetLineColor(kRed);
  hProfileNsigmaTPCPionDataV0tagVsEta->SetLineWidth(2);
  TProfile* hProfileNsigmaTPCProtonDataV0tagVsEta = hNsigmaTPCProtonDataV0tagVsEta->ProfileX("hProfileNsigmaTPCProtonDataV0tagVsEta");
  hProfileNsigmaTPCProtonDataV0tagVsEta->SetLineColor(kRed);
  hProfileNsigmaTPCProtonDataV0tagVsEta->SetLineWidth(2);
  TProfile* hProfileNsigmaTPCKaonDataTOFtagVsEta = hNsigmaTPCKaonDataTOFtagVsEta->ProfileX("hProfileNsigmaTPCKaonDataTOFtagVsEta");
  hProfileNsigmaTPCKaonDataTOFtagVsEta->SetLineColor(kRed);
  hProfileNsigmaTPCKaonDataTOFtagVsEta->SetLineWidth(2);
  TProfile* hProfileNsigmaTPCKaonDataKinktagVsEta = hNsigmaTPCKaonDataKinktagVsEta->ProfileX("hProfileNsigmaTPCKaonDataKinktagVsEta");
  hProfileNsigmaTPCKaonDataKinktagVsEta->SetLineColor(kRed);
  hProfileNsigmaTPCKaonDataKinktagVsEta->SetLineWidth(2);
  
  TProfile* hProfileNsigmaTOFPionDataV0tagVsEta = hNsigmaTOFPionDataV0tagVsEta->ProfileX("hProfileNsigmaTOFPionDataV0tagVsEta");
  hProfileNsigmaTOFPionDataV0tagVsEta->SetLineColor(kRed);
  hProfileNsigmaTOFPionDataV0tagVsEta->SetLineWidth(2);
  TProfile* hProfileNsigmaTOFProtonDataV0tagVsEta = hNsigmaTOFProtonDataV0tagVsEta->ProfileX("hProfileNsigmaTOFProtonDataV0tagVsEta");
  hProfileNsigmaTOFProtonDataV0tagVsEta->SetLineColor(kRed);
  hProfileNsigmaTOFProtonDataV0tagVsEta->SetLineWidth(2);
  TProfile* hProfileNsigmaTOFKaonDataTPCtagVsEta = hNsigmaTOFKaonDataTPCtagVsEta->ProfileX("hProfileNsigmaTOFKaonDataTPCtagVsEta");
  hProfileNsigmaTOFKaonDataTPCtagVsEta->SetLineColor(kRed);
  hProfileNsigmaTOFKaonDataTPCtagVsEta->SetLineWidth(2);
  TProfile* hProfileNsigmaTOFKaonDataKinktagVsEta = hNsigmaTOFKaonDataKinktagVsEta->ProfileX("hProfileNsigmaTOFKaonDataKinktagVsEta");
  hProfileNsigmaTOFKaonDataKinktagVsEta->SetLineColor(kRed);
  hProfileNsigmaTOFKaonDataKinktagVsEta->SetLineWidth(2);

  TCanvas* cNsigmaTPCVsEta = new TCanvas("cNsigmaTPCVsEta","",1920,1080);
  cNsigmaTPCVsEta->Divide(2,2);
  cNsigmaTPCVsEta->cd(1)->SetLogz();
  hNsigmaTPCPionDataV0tagVsEta->Draw("colz");
  hProfileNsigmaTPCPionDataV0tagVsEta->Draw("same");
  cNsigmaTPCVsEta->cd(2)->SetLogz();
  hNsigmaTPCProtonDataV0tagVsEta->Draw("colz");
  hProfileNsigmaTPCProtonDataV0tagVsEta->Draw("same");
  cNsigmaTPCVsEta->cd(3)->SetLogz();
  hNsigmaTPCKaonDataTOFtagVsEta->Draw("colz");
  hProfileNsigmaTPCKaonDataTOFtagVsEta->Draw("same");
  cNsigmaTPCVsEta->cd(4)->SetLogz();
  hNsigmaTPCKaonDataKinktagVsEta->Draw("colz");
  hProfileNsigmaTPCKaonDataKinktagVsEta->Draw("same");

  TCanvas* cNsigmaTOFVsEta = new TCanvas("cNsigmaTOFVsEta","",1920,1080);
  cNsigmaTOFVsEta->Divide(2,2);
  cNsigmaTOFVsEta->cd(1)->SetLogz();
  hNsigmaTOFPionDataV0tagVsEta->Draw("colz");
  hProfileNsigmaTOFPionDataV0tagVsEta->Draw("same");
  cNsigmaTOFVsEta->cd(2)->SetLogz();
  hNsigmaTOFProtonDataV0tagVsEta->Draw("colz");
  hProfileNsigmaTOFProtonDataV0tagVsEta->Draw("same");
  cNsigmaTOFVsEta->cd(3)->SetLogz();
  hNsigmaTOFKaonDataTPCtagVsEta->Draw("colz");
  hProfileNsigmaTOFKaonDataTPCtagVsEta->Draw("same");
  cNsigmaTOFVsEta->cd(4)->SetLogz();
  hNsigmaTOFKaonDataKinktagVsEta->Draw("colz");
  hProfileNsigmaTOFKaonDataKinktagVsEta->Draw("same");

  cNsigmaTPCVsEta->SaveAs(Form("%s/NsigmaTPCvsEta.pdf",outputdirName.Data()));
  cNsigmaTOFVsEta->SaveAs(Form("%s/NsigmaTOFvsEta.pdf",outputdirName.Data()));

  //********************************************************************************************************************************************//
  //compute Fraction and contamination in MC
  double intNsigmaTPCPionMCV0tag[nBins][nPDGcodes];
  double intNsigmaTPCKaonMCKinktag[nBins][nPDGcodes];
  double intNsigmaTPCKaonMCTOFtag[nBins][nPDGcodes];
  double intNsigmaTPCProtonMCV0tag[nBins][nPDGcodes];
  double intNsigmaTOFPionMCV0tag[nBins][nPDGcodes];
  double intNsigmaTOFKaonMCKinktag[nBins][nPDGcodes];
  double intNsigmaTOFKaonMCTPCtag[nBins][nPDGcodes];
  double intNsigmaTOFProtonMCV0tag[nBins][nPDGcodes];

  double intNsigmaTPCPionDataV0tag[nBins];
  double intNsigmaTPCKaonDataKinktag[nBins];
  double intNsigmaTPCKaonDataTOFtag[nBins];
  double intNsigmaTPCProtonDataV0tag[nBins];
  double intNsigmaTOFPionDataV0tag[nBins];
  double intNsigmaTOFKaonDataKinktag[nBins];
  double intNsigmaTOFKaonDataTPCtag[nBins];
  double intNsigmaTOFProtonDataV0tag[nBins];
  
  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      intNsigmaTPCPionMCV0tag[iBin][iPart]   = hNsigmaTPCPionMCV0tag[iBin][iPart]->Integral() / hNsigmaTPCPionMCV0tag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTPCKaonMCKinktag[iBin][iPart] = hNsigmaTPCKaonMCKinktag[iBin][iPart]->Integral() / hNsigmaTPCKaonMCKinktag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTPCKaonMCTOFtag[iBin][iPart]  = hNsigmaTPCKaonMCTOFtag[iBin][iPart]->Integral() / hNsigmaTPCKaonMCTOFtag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTPCProtonMCV0tag[iBin][iPart] = hNsigmaTPCProtonMCV0tag[iBin][iPart]->Integral() / hNsigmaTPCProtonMCV0tag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTOFPionMCV0tag[iBin][iPart]   = hNsigmaTOFPionMCV0tag[iBin][iPart]->Integral() / hNsigmaTOFPionMCV0tag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTOFKaonMCKinktag[iBin][iPart] = hNsigmaTOFKaonMCKinktag[iBin][iPart]->Integral() / hNsigmaTOFKaonMCKinktag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTOFKaonMCTPCtag[iBin][iPart]  = hNsigmaTOFKaonMCTPCtag[iBin][iPart]->Integral() / hNsigmaTOFKaonMCTPCtag[iBin][iPart]->GetBinWidth(1);
      intNsigmaTOFProtonMCV0tag[iBin][iPart] = hNsigmaTOFProtonMCV0tag[iBin][iPart]->Integral() / hNsigmaTOFProtonMCV0tag[iBin][iPart]->GetBinWidth(1);
    }
    
    intNsigmaTPCPionDataV0tag[iBin] = hNsigmaTPCPionDataV0tag[iBin]->Integral() / hNsigmaTPCPionDataV0tag[iBin]->GetBinWidth(1);
    intNsigmaTPCKaonDataKinktag[iBin] = hNsigmaTPCKaonDataKinktag[iBin]->Integral() / hNsigmaTPCKaonDataKinktag[iBin]->GetBinWidth(1);
    intNsigmaTPCKaonDataTOFtag[iBin] = hNsigmaTPCKaonDataTOFtag[iBin]->Integral() / hNsigmaTPCKaonDataTOFtag[iBin]->GetBinWidth(1);
    intNsigmaTPCProtonDataV0tag[iBin] = hNsigmaTPCProtonDataV0tag[iBin]->Integral() / hNsigmaTPCProtonDataV0tag[iBin]->GetBinWidth(1);
    intNsigmaTOFPionDataV0tag[iBin] = hNsigmaTOFPionDataV0tag[iBin]->Integral() / hNsigmaTOFPionDataV0tag[iBin]->GetBinWidth(1);
    intNsigmaTOFKaonDataKinktag[iBin] = hNsigmaTOFKaonDataKinktag[iBin]->Integral() / hNsigmaTOFKaonDataKinktag[iBin]->GetBinWidth(1);
    intNsigmaTOFKaonDataTPCtag[iBin] = hNsigmaTOFKaonDataTPCtag[iBin]->Integral() / hNsigmaTOFKaonDataTPCtag[iBin]->GetBinWidth(1);
    intNsigmaTOFProtonDataV0tag[iBin] = hNsigmaTOFProtonDataV0tag[iBin]->Integral() / hNsigmaTOFProtonDataV0tag[iBin]->GetBinWidth(1);
  }
  
  TH1F *hFractionTPCPionMCV0tag[nPDGcodes-1], *hFractionTPCKaonMCKinktag[nPDGcodes-1], *hFractionTPCKaonMCTOFtag[nPDGcodes-1], *hFractionTPCProtonMCV0tag[nPDGcodes-1];
  TH1F *hFractionTOFPionMCV0tag[nPDGcodes-1], *hFractionTOFKaonMCKinktag[nPDGcodes-1], *hFractionTOFKaonMCTPCtag[nPDGcodes-1], *hFractionTOFProtonMCV0tag[nPDGcodes-1];

  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    hFractionTPCPionMCV0tag[iPart] = new TH1F(Form("hFractionTPCPionMCV0tag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTOFPionMCV0tag[iPart] = new TH1F(Form("hFractionTOFPionMCV0tag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTPCKaonMCKinktag[iPart] = new TH1F(Form("hFractionTPCKaonMCKinktag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTPCKaonMCTOFtag[iPart] = new TH1F(Form("hFractionTPCKaonMCTOFtag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTOFKaonMCKinktag[iPart] = new TH1F(Form("hFractionTOFKaonMCKinktag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTOFKaonMCTPCtag[iPart] = new TH1F(Form("hFractionTOFKaonMCTPCtag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTPCProtonMCV0tag[iPart] = new TH1F(Form("hFractionTPCProtonMCV0tag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);
    hFractionTOFProtonMCV0tag[iPart] = new TH1F(Form("hFractionTOFProtonMCV0tag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Contamination",varname.Data()),nBins,binlims);

    SetTH1Style(hFractionTPCPionMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFPionMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTPCKaonMCKinktag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTPCKaonMCTOFtag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFKaonMCKinktag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFKaonMCTPCtag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTPCProtonMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFProtonMCV0tag[iPart],kFullCircle,pdgcolors[iPart],1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
  }
  
  for(int iBin=0; iBin<nBins; iBin++) {
    double eff=-1, unc=-1;
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      ComputeEfficiency(intNsigmaTPCPionMCV0tag[iBin][iPart],intNsigmaTPCPionMCV0tag[iBin][kAll],eff,unc);
      hFractionTPCPionMCV0tag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTPCPionMCV0tag[iPart]->SetBinError(iBin+1,unc);

      ComputeEfficiency(intNsigmaTOFPionMCV0tag[iBin][iPart],intNsigmaTOFPionMCV0tag[iBin][kAll],eff,unc);
      hFractionTOFPionMCV0tag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTOFPionMCV0tag[iPart]->SetBinError(iBin+1,unc);

      ComputeEfficiency(intNsigmaTPCKaonMCKinktag[iBin][iPart],intNsigmaTPCKaonMCKinktag[iBin][kAll],eff,unc);
      hFractionTPCKaonMCKinktag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTPCKaonMCKinktag[iPart]->SetBinError(iBin+1,unc);

      ComputeEfficiency(intNsigmaTPCKaonMCTOFtag[iBin][iPart],intNsigmaTPCKaonMCTOFtag[iBin][kAll],eff,unc);
      hFractionTPCKaonMCTOFtag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTPCKaonMCTOFtag[iPart]->SetBinError(iBin+1,unc);

      ComputeEfficiency(intNsigmaTOFKaonMCKinktag[iBin][iPart],intNsigmaTOFKaonMCKinktag[iBin][kAll],eff,unc);
      hFractionTOFKaonMCKinktag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTOFKaonMCKinktag[iPart]->SetBinError(iBin+1,unc);

      ComputeEfficiency(intNsigmaTOFKaonMCTPCtag[iBin][iPart],intNsigmaTOFKaonMCTPCtag[iBin][kAll],eff,unc);
      hFractionTOFKaonMCTPCtag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTOFKaonMCTPCtag[iPart]->SetBinError(iBin+1,unc);
    
      ComputeEfficiency(intNsigmaTPCProtonMCV0tag[iBin][iPart],intNsigmaTPCProtonMCV0tag[iBin][kAll],eff,unc);
      hFractionTPCProtonMCV0tag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTPCProtonMCV0tag[iPart]->SetBinError(iBin+1,unc);
        
      ComputeEfficiency(intNsigmaTOFProtonMCV0tag[iBin][iPart],intNsigmaTOFProtonMCV0tag[iBin][kAll],eff,unc);
      hFractionTOFProtonMCV0tag[iPart]->SetBinContent(iBin+1,eff);
      hFractionTOFProtonMCV0tag[iPart]->SetBinError(iBin+1,unc);
    }
  }
  
  TLegend* legFracMC = new TLegend(0.2,0.4,0.45,0.7);
  legFracMC->SetTextSize(0.05);
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    legFracMC->AddEntry(hFractionTOFProtonMCV0tag[iPart],pdgnames[iPart].Data(),"lpe");
  }

  TCanvas* cFractionTPCMC = new TCanvas("cFractionTPCMC","cFractionTPCMC",1920,1080);
  cFractionTPCMC->Divide(2,2);
  cFractionTPCMC->cd(1)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TPC #pi from K_{s}^{0};%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTPCMC->cd(1)->SetLogy();
  cFractionTPCMC->cd(1)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCPionMCV0tag[iPart]->Draw("same");
  legFracMC->Draw("same");
  
  cFractionTPCMC->cd(2)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TPC K from kinks;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTPCMC->cd(2)->SetLogy();
  cFractionTPCMC->cd(2)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCKaonMCKinktag[iPart]->Draw("same");
  
  cFractionTPCMC->cd(3)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TPC K TOF tagged;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTPCMC->cd(3)->SetLogy();
  cFractionTPCMC->cd(3)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCKaonMCTOFtag[iPart]->Draw("same");

  cFractionTPCMC->cd(4)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TPC p from #Lambda;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTPCMC->cd(4)->SetLogy();
  cFractionTPCMC->cd(4)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTPCProtonMCV0tag[iPart]->Draw("same");
  
  TCanvas* cFractionTOFMC = new TCanvas("cFractionTOFMC","cFractionTOFMC",1920,1080);
  cFractionTOFMC->Divide(2,2);
  cFractionTOFMC->cd(1)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TOF #pi from K_{s}^{0};%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTOFMC->cd(1)->SetLogy();
  cFractionTOFMC->cd(1)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFPionMCV0tag[iPart]->Draw("same");
  legFracMC->Draw("same");
  
  cFractionTOFMC->cd(2)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TOF K from kinks;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTOFMC->cd(2)->SetLogy();
  cFractionTOFMC->cd(2)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFKaonMCKinktag[iPart]->Draw("same");
  
  cFractionTOFMC->cd(3)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TOF K TPC tag;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTOFMC->cd(3)->SetLogy();
  cFractionTOFMC->cd(3)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFKaonMCTPCtag[iPart]->Draw("same");
  
  cFractionTOFMC->cd(4)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TOF p from #Lambda;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cFractionTOFMC->cd(4)->SetLogy();
  cFractionTOFMC->cd(45)->SetLogx();
  for(int iPart=kElectron; iPart<=kProton; iPart++) hFractionTOFProtonMCV0tag[iPart]->Draw("same");

  //********************************************************************************************************************************************//
  //compute fractions in data (only TOF)
  
  TFractionFitter *fNsigmaTOFKaonTPCtagFitter[nBins];
  TFractionFitter *fNsigmaTOFProtonV0tagFitter[nBins];
  TH1F *hFractionTOFKaonDataTPCtag[nPDGcodes-1];
  TH1F *hFractionTOFProtonDataV0tag[nPDGcodes-1];
  TH1F *hNsigmaTOFKaonDataTPCtag_Fit[nBins][nPDGcodes];
  TH1F *hNsigmaTOFProtonDataV0tag_Fit[nBins][nPDGcodes];
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    hFractionTOFKaonDataTPCtag[iPart] = new TH1F(Form("hFractionTOFKaonDataTPCtag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Fraction",varname.Data()),nBins,binlims);
    hFractionTOFProtonDataV0tag[iPart] = new TH1F(Form("hFractionTOFProtonDataV0tag_%s",pdgnames[iPart].Data()),Form(";%s (GeV/#it{c});Fraction",varname.Data()),nBins,binlims);
    SetTH1Style(hFractionTOFKaonDataTPCtag[iPart],kOpenSquare,pdgcolors[iPart]+1,1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
    SetTH1Style(hFractionTOFProtonDataV0tag[iPart],kOpenSquare,pdgcolors[iPart]+1,1.,2,pdgcolors[iPart],kWhite,0.045,0.05);
  }

  TCanvas* cFitResultTOFKaonFromTPCtag = new TCanvas("cFitResultTOFKaonFromTPCtag","cFitResultTOFKaonFromTPCtag",1920,1080);
  DivideCanvas(cFitResultTOFKaonFromTPCtag,nBins);
  TCanvas* cFitResultTOFProtonFromV0tag = new TCanvas("cFitResultTOFProtonFromV0tag","cFitResultTOFProtonFromV0tag",1920,1080);
  DivideCanvas(cFitResultTOFProtonFromV0tag,nBins);
  TCanvas* cTOFFractionData = new TCanvas("cTOFFractionData","cTOFFractionData",1920,1080);
  cTOFFractionData->Divide(3,1);

  TLegend* legTOFFitter = new TLegend(0.14,0.62,0.52,0.86);
  legTOFFitter->SetTextSize(0.04);
  legTOFFitter->AddEntry(hNsigmaTOFKaonDataTPCtag[0],"Data","p");
  for(int iPart=kElectron; iPart<=kProton; iPart++) {
    legTOFFitter->AddEntry(hNsigmaTOFKaonMCTPCtag[0][iPart],Form("Templ %s",pdgnames[iPart].Data()),"l");
  }

  for(int iBin=0; iBin<nBins; iBin++) {

    vector<int> templUsedTOFKaonMCTPCtag;
    GetTOFFractionsFromData(kKaon,iBin,hFractionTOFKaonMCTPCtag,hFractionTOFKaonDataTPCtag,hNsigmaTOFKaonMCTPCtag[iBin],hNsigmaTOFKaonDataTPCtag[iBin],fNsigmaTOFKaonTPCtagFitter[iBin],templUsedTOFKaonMCTPCtag);

    if(templUsedTOFKaonMCTPCtag.size()>1) {
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll] = (TH1F*)fNsigmaTOFKaonTPCtagFitter[iBin]->GetPlot();
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->SetLineColor(kRed);
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->SetFillStyle(0);
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->SetFillColor(kWhite);
      
      for(unsigned int iTempl=0; iTempl<templUsedTOFKaonMCTPCtag.size(); iTempl++) {
        double frac, err;
        fNsigmaTOFKaonTPCtagFitter[iBin]->GetResult(iTempl, frac, err);
        hNsigmaTOFKaonDataTPCtag_Fit[iBin][templUsedTOFKaonMCTPCtag[iTempl]] = (TH1F*)fNsigmaTOFKaonTPCtagFitter[iBin]->GetMCPrediction(iTempl);
        hNsigmaTOFKaonDataTPCtag_Fit[iBin][templUsedTOFKaonMCTPCtag[iTempl]]->Scale(frac/hNsigmaTOFKaonDataTPCtag_Fit[iBin][templUsedTOFKaonMCTPCtag[iTempl]]->Integral()*hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->Integral());
        hNsigmaTOFKaonDataTPCtag_Fit[iBin][templUsedTOFKaonMCTPCtag[iTempl]]->SetFillColor(kWhite);
        hNsigmaTOFKaonDataTPCtag_Fit[iBin][templUsedTOFKaonMCTPCtag[iTempl]]->SetFillStyle(0);
      }
    }
    else {
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll] = (TH1F*)hNsigmaTOFKaonDataTPCtag[iBin]->Clone();
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->SetLineColor(kRed);
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->SetFillStyle(0);
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kAll]->SetFillColor(kWhite);
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kKaon] = (TH1F*)hNsigmaTOFKaonDataTPCtag[iBin]->Clone();
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kKaon]->SetFillColor(kWhite);
      hNsigmaTOFKaonDataTPCtag_Fit[iBin][kKaon]->SetFillStyle(0);
    }
    cFitResultTOFKaonFromTPCtag->cd(iBin+1)->SetLogy();
    hNsigmaTOFKaonDataTPCtag[iBin]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
    vector<int>::iterator it = find(templUsedTOFKaonMCTPCtag.begin(),templUsedTOFKaonMCTPCtag.end(),iPart);
      if(it!=templUsedTOFKaonMCTPCtag.end()) {
        hNsigmaTOFKaonDataTPCtag_Fit[iBin][iPart]->DrawCopy("hist same");
      }
      else {
        hNsigmaTOFKaonDataTPCtag_Fit[iBin][iPart]=nullptr;
      }
    }
    legTOFFitter->Draw("same");

    vector<int> templUsedTOFProtonMCV0tag;
    GetTOFFractionsFromData(kProton,iBin,hFractionTOFProtonMCV0tag,hFractionTOFProtonDataV0tag,hNsigmaTOFProtonMCV0tag[iBin],hNsigmaTOFProtonDataV0tag[iBin],fNsigmaTOFProtonV0tagFitter[iBin],templUsedTOFProtonMCV0tag);

    if(templUsedTOFProtonMCV0tag.size()>1) {
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll] = (TH1F*)fNsigmaTOFProtonV0tagFitter[iBin]->GetPlot();
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->SetLineColor(kRed);
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->SetFillStyle(0);
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->SetFillColor(kWhite);
      
      for(unsigned int iTempl=0; iTempl<templUsedTOFProtonMCV0tag.size(); iTempl++) {
        double frac, err;
        fNsigmaTOFProtonV0tagFitter[iBin]->GetResult(iTempl, frac, err);
        hNsigmaTOFProtonDataV0tag_Fit[iBin][templUsedTOFProtonMCV0tag[iTempl]] = (TH1F*)fNsigmaTOFProtonV0tagFitter[iBin]->GetMCPrediction(iTempl);
        hNsigmaTOFProtonDataV0tag_Fit[iBin][templUsedTOFProtonMCV0tag[iTempl]]->Scale(frac/hNsigmaTOFProtonDataV0tag_Fit[iBin][templUsedTOFProtonMCV0tag[iTempl]]->Integral()*hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->Integral());
        hNsigmaTOFProtonDataV0tag_Fit[iBin][templUsedTOFProtonMCV0tag[iTempl]]->SetFillColor(kWhite);
        hNsigmaTOFProtonDataV0tag_Fit[iBin][templUsedTOFProtonMCV0tag[iTempl]]->SetFillStyle(0);
      }
    }
    else {
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll] = (TH1F*)hNsigmaTOFProtonDataV0tag[iBin]->Clone();
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->SetLineColor(kRed);
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->SetFillStyle(0);
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kAll]->SetFillColor(kWhite);
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kProton] = (TH1F*)hNsigmaTOFProtonDataV0tag[iBin]->Clone();
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kProton]->SetFillColor(kWhite);
      hNsigmaTOFProtonDataV0tag_Fit[iBin][kProton]->SetFillStyle(0);
    }
    cFitResultTOFProtonFromV0tag->cd(iBin+1)->SetLogy();
    hNsigmaTOFProtonDataV0tag[iBin]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
    vector<int>::iterator it = find(templUsedTOFProtonMCV0tag.begin(),templUsedTOFProtonMCV0tag.end(),iPart);
      if(it!=templUsedTOFProtonMCV0tag.end()) {
        hNsigmaTOFProtonDataV0tag_Fit[iBin][iPart]->DrawCopy("hist same");
      }
      else {
        hNsigmaTOFProtonDataV0tag_Fit[iBin][iPart]=nullptr;
      }
    }
    legTOFFitter->Draw("same");
  }

  cTOFFractionData->cd(1)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TOF K TPC tag;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
  cTOFFractionData->cd(1)->SetLogy();
  cTOFFractionData->cd(1)->SetLogx();
  for(int iPart=0; iPart<=kProton; iPart++) {
    hFractionTOFKaonDataTPCtag[iPart]->DrawCopy("same");
  }
  cTOFFractionData->cd(2)->DrawFrame(binlims[0],1.e-5,binlims[nBins],10.,Form("TOF p from #Lambda;%s (GeV/#it{c});Purity / Contamination",varname.Data()));
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
  
  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      hNsigmaTPCPionMCV0tag[iBin][iPart]->Scale(1./hNsigmaTPCPionMCV0tag[iBin][kAll]->Integral());
      hNsigmaTPCKaonMCKinktag[iBin][iPart]->Scale(1./hNsigmaTPCKaonMCKinktag[iBin][kAll]->Integral());
      hNsigmaTPCKaonMCTOFtag[iBin][iPart]->Scale(1./hNsigmaTPCKaonMCTOFtag[iBin][kAll]->Integral());
      hNsigmaTPCProtonMCV0tag[iBin][iPart]->Scale(1./hNsigmaTPCProtonMCV0tag[iBin][kAll]->Integral());
      hNsigmaTOFPionMCV0tag[iBin][iPart]->Scale(1./hNsigmaTOFPionMCV0tag[iBin][kAll]->Integral());
      hNsigmaTOFKaonMCKinktag[iBin][iPart]->Scale(1./hNsigmaTOFKaonMCKinktag[iBin][kAll]->Integral());
      hNsigmaTOFKaonMCTPCtag[iBin][iPart]->Scale(1./hNsigmaTOFKaonMCTPCtag[iBin][kAll]->Integral());
      hNsigmaTOFProtonMCV0tag[iBin][iPart]->Scale(1./hNsigmaTOFProtonMCV0tag[iBin][kAll]->Integral());

      if(iPart==kAll || hFractionTOFKaonDataTPCtag[iPart]->GetBinContent(iBin+1)>0) hNsigmaTOFKaonDataTPCtag_Fit[iBin][iPart]->Scale(1./hNsigmaTOFKaonDataTPCtag[iBin]->Integral());
      if(iPart==kAll || hFractionTOFProtonDataV0tag[iPart]->GetBinContent(iBin+1)>0) hNsigmaTOFProtonDataV0tag_Fit[iBin][iPart]->Scale(1./hNsigmaTOFProtonDataV0tag[iBin]->Integral());
    }

    hNsigmaTPCPionDataV0tag[iBin]->Scale(1./hNsigmaTPCPionDataV0tag[iBin]->Integral());
    hNsigmaTPCKaonDataKinktag[iBin]->Scale(1./hNsigmaTPCKaonDataKinktag[iBin]->Integral());
    hNsigmaTPCKaonDataTOFtag[iBin]->Scale(1./hNsigmaTPCKaonDataTOFtag[iBin]->Integral());
    hNsigmaTPCProtonDataV0tag[iBin]->Scale(1./hNsigmaTPCProtonDataV0tag[iBin]->Integral());
    hNsigmaTOFPionDataV0tag[iBin]->Scale(1./hNsigmaTOFPionDataV0tag[iBin]->Integral());
    hNsigmaTOFKaonDataKinktag[iBin]->Scale(1./hNsigmaTOFKaonDataKinktag[iBin]->Integral());
    hNsigmaTOFKaonDataTPCtag[iBin]->Scale(1./hNsigmaTOFKaonDataTPCtag[iBin]->Integral());
    hNsigmaTOFProtonDataV0tag[iBin]->Scale(1./hNsigmaTOFProtonDataV0tag[iBin]->Integral());
  }
  
  //********************************************************************************************************************************************//
  //draw and fit distributions

  TCanvas* cPionMCV0tagTPC = new TCanvas("cPionMCV0tagTPC","cPionMCV0tagTPC",1920,1080);
  DivideCanvas(cPionMCV0tagTPC,nBins);
  TCanvas* cKaonMCKinkstagTPC = new TCanvas("cKaonMCKinkstagTPC","cKaonMCKinkstagTPC",1920,1080);
  DivideCanvas(cKaonMCKinkstagTPC,nBins);
  TCanvas* cKaonMCTOFtagTPC = new TCanvas("cKaonMCTOFtagTPC","cKaonMCTOFtagTPC",1920,1080);
  DivideCanvas(cKaonMCTOFtagTPC,nBins);
  TCanvas* cProtonMCV0tagTPC = new TCanvas("cProtonMCV0tagTPC","cProtonMCV0tagTPC",1920,1080);
  DivideCanvas(cProtonMCV0tagTPC,nBins);
  TCanvas* cPionMCV0tagTOF = new TCanvas("cPionMCV0tagTOF","cPionMCV0tagTOF",1920,1080);
  DivideCanvas(cPionMCV0tagTOF,nBins);
  TCanvas* cKaonMCKinkstagTOF = new TCanvas("cKaonMCKinkstagTOF","cKaonMCKinkstagTOF",1920,1080);
  DivideCanvas(cKaonMCKinkstagTOF,nBins);
  TCanvas* cKaonMCTPCtagTOF = new TCanvas("cKaonMCTPCtagTOF","cKaonMCTPCtagTOF",1920,1080);
  DivideCanvas(cKaonMCTPCtagTOF,nBins);
  TCanvas* cProtonMCV0tagTOF = new TCanvas("cProtonMCV0tagTOF","cProtonMCV0tagTOF",1920,1080);
  DivideCanvas(cProtonMCV0tagTOF,nBins);

  TCanvas* cPionDataV0tagTPC = new TCanvas("cPionDataV0tagTPC","cPionDataV0tagTPC",1920,1080);
  DivideCanvas(cPionDataV0tagTPC,nBins);
  TCanvas* cKaonDataKinkstagTPC = new TCanvas("cKaonDataKinkstagTPC","cKaonDataKinkstagTPC",1920,1080);
  DivideCanvas(cKaonDataKinkstagTPC,nBins);
  TCanvas* cKaonDataTOFtagTPC = new TCanvas("cKaonDataTOFtagTPC","cKaonDataTOFtagTPC",1920,1080);
  DivideCanvas(cKaonDataTOFtagTPC,nBins);
  TCanvas* cProtonDataV0tagTPC = new TCanvas("cProtonDataV0tagTPC","cProtonDataV0tagTPC",1920,1080);
  DivideCanvas(cProtonDataV0tagTPC,nBins);
  TCanvas* cPionDataV0tagTOF = new TCanvas("cPionDataV0tagTOF","cPionDataV0tagTOF",1920,1080);
  DivideCanvas(cPionDataV0tagTOF,nBins);
  TCanvas* cKaonDataKinkstagTOF = new TCanvas("cKaonDataKinkstagTOF","cKaonDataKinkstagTOF",1920,1080);
  DivideCanvas(cKaonDataKinkstagTOF,nBins);
  TCanvas* cKaonDataTPCtagTOF = new TCanvas("cKaonDataTPCtagTOF","cKaonDataTPCtagTOF",1920,1080);
  DivideCanvas(cKaonDataTPCtagTOF,nBins);
  TCanvas* cProtonDataV0tagTOF = new TCanvas("cProtonDataV0tagTOF","cProtonDataV0tagTOF",1920,1080);
  DivideCanvas(cProtonDataV0tagTOF,nBins);

  //fit TPC nsigma
  TF1 *fNsigmaTPCPionMCV0tag[nBins][nPDGcodes-1], *fNsigmaTPCKaonMCKinktag[nBins][nPDGcodes-1], *fNsigmaTPCKaonMCTOFtag[nBins][nPDGcodes-1], *fNsigmaTPCProtonMCV0tag[nBins][nPDGcodes-1];
  TF1* fNsigmaTPCPionDataV0tag[nBins][nPDGcodes], *fNsigmaTPCKaonDataKinktag[nBins][nPDGcodes], *fNsigmaTPCKaonDataTOFtag[nBins][nPDGcodes], *fNsigmaTPCProtonDataV0tag[nBins][nPDGcodes];

  //subtract bkg template for TOF nsigma
  TH1F *hNsigmaTOFPionDataV0tag_sub[nBins], *hNsigmaTOFKaonDataKinktag_sub[nBins], *hNsigmaTOFKaonDataTPCtag_sub[nBins], *hNsigmaTOFProtonDataV0tag_sub[nBins];

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

  TString fitopt="ML0R";
  for(int iBin=0; iBin<nBins; iBin++) {

    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      fNsigmaTPCPionMCV0tag[iBin][iPart] = new TF1(Form("fNsigmaTPCPionMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonMCKinktag[iBin][iPart] = new TF1(Form("fNsigmaTPCKaonMCKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonMCTOFtag[iBin][iPart] = new TF1(Form("fNsigmaTPCKaonMCTOFtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
      fNsigmaTPCProtonMCV0tag[iBin][iPart] = new TF1(Form("fNsigmaTPCProtonMCV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
    }
    
    cPionMCV0tagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCPionMCV0tag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCPionMCV0tag[iBin][iPart]->DrawCopy("hist same");
      if(iPart==kPion || hFractionTPCPionMCV0tag[iPart]->GetBinContent(iBin+1)>1.e-5) hNsigmaTPCPionMCV0tag[iBin][iPart]->Fit(fNsigmaTPCPionMCV0tag[iBin][iPart],fitopt.Data());
      else fNsigmaTPCPionMCV0tag[iBin][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCPionMCV0tag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCKinkstagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCKaonMCKinktag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCKaonMCKinktag[iBin][iPart]->DrawCopy("hist same");
      if(iPart==kKaon || hFractionTPCKaonMCKinktag[iPart]->GetBinContent(iBin+1)>1.e-5) hNsigmaTPCKaonMCKinktag[iBin][iPart]->Fit(fNsigmaTPCKaonMCKinktag[iBin][iPart],fitopt.Data());
      else fNsigmaTPCKaonMCKinktag[iBin][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCKaonMCKinktag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCTOFtagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCKaonMCTOFtag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCKaonMCTOFtag[iBin][iPart]->DrawCopy("hist same");
      if(iPart==kKaon || hFractionTPCKaonMCTOFtag[iPart]->GetBinContent(iBin+1)>1.e-5) hNsigmaTPCKaonMCTOFtag[iBin][iPart]->Fit(fNsigmaTPCKaonMCTOFtag[iBin][iPart],fitopt.Data());
      else fNsigmaTPCKaonMCTOFtag[iBin][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCKaonMCTOFtag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cProtonMCV0tagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCProtonMCV0tag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      hNsigmaTPCProtonMCV0tag[iBin][iPart]->DrawCopy("hist same");
      if(iPart==kProton || hFractionTPCProtonMCV0tag[iPart]->GetBinContent(iBin+1)>1.e-5) hNsigmaTPCProtonMCV0tag[iBin][iPart]->Fit(fNsigmaTPCProtonMCV0tag[iBin][iPart],fitopt.Data());
      else fNsigmaTPCProtonMCV0tag[iBin][iPart]->SetParameter(0,0);
    }
    hNsigmaTPCProtonMCV0tag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    
    cPionMCV0tagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFPionMCV0tag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFPionMCV0tag[iBin][iPart]->DrawCopy("hist same");
    hNsigmaTOFPionMCV0tag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCKinkstagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFKaonMCKinktag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFKaonMCKinktag[iBin][iPart]->DrawCopy("hist same");
    hNsigmaTOFKaonMCKinktag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cKaonMCTPCtagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFKaonMCTPCtag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFKaonMCTPCtag[iBin][iPart]->DrawCopy("hist same");
    hNsigmaTOFKaonMCTPCtag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    cProtonMCV0tagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFProtonMCV0tag[iBin][nPDGcodes-1]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) hNsigmaTOFProtonMCV0tag[iBin][iPart]->DrawCopy("hist same");
    hNsigmaTOFProtonMCV0tag[iBin][nPDGcodes-1]->DrawCopy("Esame");
    legMCdist->Draw("same");
    
    fNsigmaTPCPionDataV0tag[iBin][kAll] = new TF1(Form("fNsigmaTPCPionDataV0tag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),binlims[iBin]*10,binlims[iBin+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);
    fNsigmaTPCKaonDataKinktag[iBin][kAll] = new TF1(Form("fNsigmaTPCKaonDataKinktag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),binlims[iBin]*10,binlims[iBin+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);
    fNsigmaTPCKaonDataTOFtag[iBin][kAll] = new TF1(Form("fNsigmaTPCKaonDataTOFtag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),binlims[iBin]*10,binlims[iBin+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);
    fNsigmaTPCProtonDataV0tag[iBin][kAll] = new TF1(Form("fNsigmaTPCProtonDataV0tag_%s_p%0.f_%0.f",pdgnames[kAll].Data(),binlims[iBin]*10,binlims[iBin+1]*10),PDFnsigmaTPCtot,-50,50,(nPDGcodes-1)*3);

    for(int iPart=kElectron; iPart<=kProton; iPart++) {
     
      double toll = 100.;
      for(int iPar=0; iPar<3; iPar++) {
        if(iPar==0)
          toll=5.;
        else
          toll=0.5;
        fNsigmaTPCPionDataV0tag[iBin][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCPionMCV0tag[iBin][iPart]->GetParameter(iPar));
        fNsigmaTPCKaonDataKinktag[iBin][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCKaonMCKinktag[iBin][iPart]->GetParameter(iPar));
        fNsigmaTPCKaonDataTOFtag[iBin][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCKaonMCTOFtag[iBin][iPart]->GetParameter(iPar));
        fNsigmaTPCProtonDataV0tag[iBin][kAll]->SetParameter(iPart*3+iPar,fNsigmaTPCProtonMCV0tag[iBin][iPart]->GetParameter(iPar));
        
        if(iPart*3+iPar!=kPion*3+1) {
          fNsigmaTPCPionDataV0tag[iBin][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCPionMCV0tag[iBin][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCPionMCV0tag[iBin][iPart]->GetParameter(iPar))*toll,fNsigmaTPCPionMCV0tag[iBin][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCPionMCV0tag[iBin][iPart]->GetParameter(iPar))*toll);
        }
        if(iPart*3+iPar!=kKaon*3+1) {
        fNsigmaTPCKaonDataKinktag[iBin][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCKaonMCKinktag[iBin][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCKaonMCKinktag[iBin][iPart]->GetParameter(iPar))*toll,fNsigmaTPCKaonMCKinktag[iBin][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCKaonMCKinktag[iBin][iPart]->GetParameter(iPar))*toll);
        
        fNsigmaTPCKaonDataTOFtag[iBin][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCKaonMCTOFtag[iBin][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCKaonMCTOFtag[iBin][iPart]->GetParameter(iPar))*toll,fNsigmaTPCKaonMCTOFtag[iBin][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCKaonMCTOFtag[iBin][iPart]->GetParameter(iPar))*toll);
        }
        if(iPart*3+iPar!=kProton*3+1) {
        fNsigmaTPCProtonDataV0tag[iBin][kAll]->SetParLimits(iPart*3+iPar,fNsigmaTPCProtonMCV0tag[iBin][iPart]->GetParameter(iPar)-TMath::Abs(fNsigmaTPCProtonMCV0tag[iBin][iPart]->GetParameter(iPar))*toll,fNsigmaTPCProtonMCV0tag[iBin][iPart]->GetParameter(iPar)+TMath::Abs(fNsigmaTPCProtonMCV0tag[iBin][iPart]->GetParameter(iPar))*toll);
        }
      }
    }

    hNsigmaTPCPionDataV0tag[iBin]->Fit(fNsigmaTPCPionDataV0tag[iBin][kAll],fitopt.Data());
    hNsigmaTPCKaonDataKinktag[iBin]->Fit(fNsigmaTPCKaonDataKinktag[iBin][kAll],fitopt.Data());
    hNsigmaTPCKaonDataTOFtag[iBin]->Fit(fNsigmaTPCKaonDataTOFtag[iBin][kAll],fitopt.Data());
    hNsigmaTPCProtonDataV0tag[iBin]->Fit(fNsigmaTPCProtonDataV0tag[iBin][kAll],fitopt.Data());

    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      fNsigmaTPCPionDataV0tag[iBin][iPart] = new TF1(Form("fNsigmaTPCPionDataV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonDataKinktag[iBin][iPart] = new TF1(Form("fNsigmaTPCKaonDataKinktag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
      fNsigmaTPCKaonDataTOFtag[iBin][iPart] = new TF1(Form("fNsigmaTPCKaonDataTOFtag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);
      fNsigmaTPCProtonDataV0tag[iBin][iPart] = new TF1(Form("fNsigmaTPCProtonDataV0tag_%s_p%0.f_%0.f",pdgnames[iPart].Data(),binlims[iBin]*10,binlims[iBin+1]*10),"gaus",-50,50);

      for(int iPar=0; iPar<3; iPar++) {
        fNsigmaTPCPionDataV0tag[iBin][iPart]->SetParameter(iPar,fNsigmaTPCPionDataV0tag[iBin][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCKaonDataKinktag[iBin][iPart]->SetParameter(iPar,fNsigmaTPCKaonDataKinktag[iBin][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCKaonDataTOFtag[iBin][iPart]->SetParameter(iPar,fNsigmaTPCKaonDataTOFtag[iBin][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCProtonDataV0tag[iBin][iPart]->SetParameter(iPar,fNsigmaTPCProtonDataV0tag[iBin][kAll]->GetParameter(iPart*3+iPar));
        fNsigmaTPCPionDataV0tag[iBin][iPart]->SetParError(iPar,fNsigmaTPCPionDataV0tag[iBin][kAll]->GetParError(iPart*3+iPar));
        fNsigmaTPCKaonDataKinktag[iBin][iPart]->SetParError(iPar,fNsigmaTPCKaonDataKinktag[iBin][kAll]->GetParError(iPart*3+iPar));
        fNsigmaTPCKaonDataTOFtag[iBin][iPart]->SetParError(iPar,fNsigmaTPCKaonDataTOFtag[iBin][kAll]->GetParError(iPart*3+iPar));
        fNsigmaTPCProtonDataV0tag[iBin][iPart]->SetParError(iPar,fNsigmaTPCProtonDataV0tag[iBin][kAll]->GetParError(iPart*3+iPar));
      }
      fNsigmaTPCPionDataV0tag[iBin][iPart]->SetLineColor(pdgcolors[iPart]);
      fNsigmaTPCKaonDataKinktag[iBin][iPart]->SetLineColor(pdgcolors[iPart]);
      fNsigmaTPCKaonDataTOFtag[iBin][iPart]->SetLineColor(pdgcolors[iPart]);
      fNsigmaTPCProtonDataV0tag[iBin][iPart]->SetLineColor(pdgcolors[iPart]);
      if(iBin==0) legTPCFitter->AddEntry(fNsigmaTPCProtonDataV0tag[iBin][iPart],Form("Func %s",pdgnames[iPart].Data()),"l");
    }

    cPionDataV0tagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCPionDataV0tag[iBin]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCPionDataV0tag[iBin][iPart]->Draw("same");
    legTPCFitter->Draw("same");
    cKaonDataKinkstagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCKaonDataKinktag[iBin]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCKaonDataKinktag[iBin][iPart]->Draw("same");
    legTPCFitter->Draw("same");
    cKaonDataTOFtagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCKaonDataTOFtag[iBin]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCKaonDataTOFtag[iBin][iPart]->Draw("same");
    legTPCFitter->Draw("same");
    cProtonDataV0tagTPC->cd(iBin+1)->SetLogy();
    hNsigmaTPCProtonDataV0tag[iBin]->DrawCopy("E");
    for(int iPart=kElectron; iPart<=kProton; iPart++) fNsigmaTPCProtonDataV0tag[iBin][iPart]->Draw("same");
    legTPCFitter->Draw("same");

    //subtract bkg template for TOF nsigma
    TString name = hNsigmaTOFPionDataV0tag[iBin]->GetName();
    name.ReplaceAll("hNsigmaTOFPionDataV0tag","hNsigmaTOFPionDataV0tag_sub");
    hNsigmaTOFPionDataV0tag_sub[iBin] = (TH1F*)hNsigmaTOFPionDataV0tag[iBin]->Clone(name.Data());
    
    name = hNsigmaTOFKaonDataKinktag[iBin]->GetName();
    name.ReplaceAll("hNsigmaTOFKaonDataKinktag","hNsigmaTOFKaonDataKinktag_sub");
    hNsigmaTOFKaonDataKinktag_sub[iBin] = (TH1F*)hNsigmaTOFKaonDataKinktag[iBin]->Clone(name.Data());

    name = hNsigmaTOFKaonDataTPCtag[iBin]->GetName();
    name.ReplaceAll("hNsigmaTOFKaonDataTPCtag","hNsigmaTOFKaonDataTPCtag_sub");
    hNsigmaTOFKaonDataTPCtag_sub[iBin] = (TH1F*)hNsigmaTOFKaonDataTPCtag[iBin]->Clone(name.Data());

    name = hNsigmaTOFProtonDataV0tag[iBin]->GetName();
    name.ReplaceAll("hNsigmaTOFProtonDataV0tag","hNsigmaTOFProtonDataV0tag_sub");
    hNsigmaTOFProtonDataV0tag_sub[iBin] = (TH1F*)hNsigmaTOFProtonDataV0tag[iBin]->Clone(name.Data());
    
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      if(iPart!=kPion) hNsigmaTOFPionDataV0tag_sub[iBin]->Add(hNsigmaTOFPionMCV0tag[iBin][iPart],-1.);
      if(iPart!=kKaon) {
        hNsigmaTOFKaonDataKinktag_sub[iBin]->Add(hNsigmaTOFKaonMCKinktag[iBin][iPart],-1.);
        if(hFractionTOFKaonDataTPCtag[iPart]->GetBinContent(iBin+1)>0) {
          hNsigmaTOFKaonDataTPCtag_sub[iBin]->Add(hNsigmaTOFKaonDataTPCtag_Fit[iBin][iPart],-1.);
        }
        else {
          hNsigmaTOFKaonDataTPCtag_sub[iBin]->Add(hNsigmaTOFKaonMCTPCtag[iBin][iPart],-1.);
        }
      }
      if(iPart!=kProton) {
        if(hFractionTOFProtonDataV0tag[iPart]->GetBinContent(iBin+1)>0) {
          hNsigmaTOFProtonDataV0tag_sub[iBin]->Add(hNsigmaTOFProtonDataV0tag_Fit[iBin][iPart],-1.);
       }
       else {
          hNsigmaTOFProtonDataV0tag_sub[iBin]->Add(hNsigmaTOFProtonMCV0tag[iBin][iPart],-1.);
        }
      }
    }
    
    SetTH1Style(hNsigmaTOFPionDataV0tag_sub[iBin],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataKinktag_sub[iBin],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    SetTH1Style(hNsigmaTOFKaonDataTPCtag_sub[iBin],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    SetTH1Style(hNsigmaTOFProtonDataV0tag_sub[iBin],kOpenCircle,pdgcolors[kAll],0.6,2,pdgcolors[kAll],pdgfillcolors[kAll],0.055,0.06);
    if(iBin==0) {
      legTOFPionDataV0tag->AddEntry(hNsigmaTOFPionDataV0tag_sub[iBin],"Data","p");
      legTOFPionDataV0tag->AddEntry(hNsigmaTOFPionMCV0tag[iBin][kPion],"MC Pion","p");
      legTOFKaonDataKinkstag->AddEntry(hNsigmaTOFKaonDataKinktag_sub[iBin],"Data","p");
      legTOFKaonDataKinkstag->AddEntry(hNsigmaTOFKaonMCKinktag[iBin][kKaon],"MC Kaon","p");
      legTOFKaonDataTPCtag->AddEntry(hNsigmaTOFKaonDataTPCtag_sub[iBin],"Data","p");
      legTOFKaonDataTPCtag->AddEntry(hNsigmaTOFKaonMCTPCtag[iBin][kKaon],"MC Kaon","p");      
      legTOFProtonDataV0tag->AddEntry(hNsigmaTOFProtonDataV0tag_sub[iBin],"Data","p");
      legTOFProtonDataV0tag->AddEntry(hNsigmaTOFKaonMCTPCtag[iBin][kProton],"MC Proton","p");      
    }

    cPionDataV0tagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFPionDataV0tag_sub[iBin]->DrawCopy("E");
    hNsigmaTOFPionMCV0tag[iBin][kPion]->DrawCopy("hist same");
    legTOFPionDataV0tag->Draw("same");
    cKaonDataKinkstagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFKaonDataKinktag_sub[iBin]->DrawCopy("E");
    hNsigmaTOFKaonMCKinktag[iBin][kKaon]->DrawCopy("hist same");
    legTOFKaonDataKinkstag->Draw("same");
    cKaonDataTPCtagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFKaonDataTPCtag_sub[iBin]->DrawCopy("E");
    hNsigmaTOFKaonMCTPCtag[iBin][kKaon]->DrawCopy("hist same");
    legTOFKaonDataTPCtag->Draw("same");
    cProtonDataV0tagTOF->cd(iBin+1)->SetLogy();
    hNsigmaTOFProtonDataV0tag_sub[iBin]->DrawCopy("E");
    hNsigmaTOFProtonMCV0tag[iBin][kProton]->DrawCopy("hist same");
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
    hEffPionTPCMCtrue[iEff] = new TH1F(Form("hEffPionTPCMCtrue_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC #pi efficiency",varname.Data()),nBins,binlims);
    hEffPionTOFMCtrue[iEff] = new TH1F(Form("hEffPionTOFMCtrue_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF #pi efficiency",varname.Data()),nBins,binlims);

    hEffKaonTPCMCtrue[iEff] = new TH1F(Form("hEffKaonTPCMCtrue_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC K efficiency",varname.Data()),nBins,binlims);
    hEffKaonTOFMCtrue[iEff] = new TH1F(Form("hEffKaonTOFMCtrue_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF K efficiency",varname.Data()),nBins,binlims);
    
    hEffProtonTPCMCtrue[iEff] = new TH1F(Form("hEffProtonTPCMCtrue_%dsigma;",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC p efficiency",varname.Data()),nBins,binlims);
    hEffProtonTOFMCtrue[iEff] = new TH1F(Form("hEffProtonTOFMCtrue_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF p efficiency",varname.Data()),nBins,binlims);

    hEffPionTPCDataV0tag[iEff] = new TH1F(Form("hEffPionTPCDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC #pi efficiency",varname.Data()),nBins,binlims);
    hEffPionTOFDataV0tag[iEff] = new TH1F(Form("hEffPionTOFDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF #pi efficiency",varname.Data()),nBins,binlims);
    
    hEffKaonTPCDataKinktag[iEff] = new TH1F(Form("hEffKaonTPCDataKinktag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC K efficiency",varname.Data()),nBins,binlims);
    hEffKaonTOFDataKinktag[iEff] = new TH1F(Form("hEffKaonTOFDataKinktag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF K efficiency",varname.Data()),nBins,binlims);

    hEffKaonTPCDataTOFtag[iEff] = new TH1F(Form("hEffKaonTPCDataTOFtag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC K efficiency",varname.Data()),nBins,binlims);
    hEffKaonTOFDataTPCtag[iEff] = new TH1F(Form("hEffKaonTOFDataTPCtag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF K efficiency",varname.Data()),nBins,binlims);

    hEffProtonTPCDataV0tag[iEff] = new TH1F(Form("hEffProtonTPCDataV0tag_%dsigma;",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC p efficiency",varname.Data()),nBins,binlims);
    hEffProtonTOFDataV0tag[iEff] = new TH1F(Form("hEffProtonTOFDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF p efficiency",varname.Data()),nBins,binlims);

    hRatioEffPionTPCDataV0tag[iEff] = new TH1F(Form("hRatioEffPionTPCDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC #pi efficiency",varname.Data()),nBins,binlims);
    hRatioEffPionTOFDataV0tag[iEff] = new TH1F(Form("hRatioEffPionTOFDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF #pi efficiency",varname.Data()),nBins,binlims);

    hRatioEffKaonTPCDataKinktag[iEff] = new TH1F(Form("hRatioEffKaonTPCDataKinktag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC K efficiency",varname.Data()),nBins,binlims);
    hRatioEffKaonTOFDataKinktag[iEff] = new TH1F(Form("hRatioEffKaonTOFDataKinktag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF K efficiency",varname.Data()),nBins,binlims);

    hRatioEffKaonTPCDataTOFtag[iEff] = new TH1F(Form("hRatioEffKaonTPCDataTOFtag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC K efficiency",varname.Data()),nBins,binlims);
    hRatioEffKaonTOFDataTPCtag[iEff] = new TH1F(Form("hRatioEffKaonTOFDataTPCtag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF K efficiency",varname.Data()),nBins,binlims);

    hRatioEffProtonTPCDataV0tag[iEff] = new TH1F(Form("hRatioEffProtonTPCDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TPC #pi efficiency",varname.Data()),nBins,binlims);
    hRatioEffProtonTOFDataV0tag[iEff] = new TH1F(Form("hRatioEffProtonTOFDataV0tag_%dsigma",nSigma[iEff]),Form(";%s (GeV/#it{c});TOF #pi efficiency",varname.Data()),nBins,binlims);

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

    for(int iBin=0; iBin<nBins; iBin++) {
      double eff=-1, unc=-1;
      ComputeEfficiency(hNsigmaTPCPionMCTrue[iBin]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTPCPionMCTrue[iBin]->Integral(),eff,unc);
      hEffPionTPCMCtrue[iEff]->SetBinContent(iBin+1,eff);
      hEffPionTPCMCtrue[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(hNsigmaTOFPionMCTrue[iBin]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTOFPionMCTrue[iBin]->Integral(),eff,unc);
      hEffPionTOFMCtrue[iEff]->SetBinContent(iBin+1,eff);
      hEffPionTOFMCtrue[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(hNsigmaTPCKaonMCTrue[iBin]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTPCKaonMCTrue[iBin]->Integral(),eff,unc);
      hEffKaonTPCMCtrue[iEff]->SetBinContent(iBin+1,eff);
      hEffKaonTPCMCtrue[iEff]->SetBinError(iBin+1,unc);
      
      ComputeEfficiency(hNsigmaTOFKaonMCTrue[iBin]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTOFKaonMCTrue[iBin]->Integral(),eff,unc);
      hEffKaonTOFMCtrue[iEff]->SetBinContent(iBin+1,eff);
      hEffKaonTOFMCtrue[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(hNsigmaTPCProtonMCTrue[iBin]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTPCProtonMCTrue[iBin]->Integral(),eff,unc);
      hEffProtonTPCMCtrue[iEff]->SetBinContent(iBin+1,eff);
      hEffProtonTPCMCtrue[iEff]->SetBinError(iBin+1,unc);
      
      ComputeEfficiency(hNsigmaTOFProtonMCTrue[iBin]->Integral(nsigmabinlow,nsigmabinhigh),hNsigmaTOFProtonMCTrue[iBin]->Integral(),eff,unc);
      hEffProtonTOFMCtrue[iEff]->SetBinContent(iBin+1,eff);
      hEffProtonTOFMCtrue[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(fNsigmaTPCPionDataV0tag[iBin][kPion]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCPionDataV0tag[iBin],fNsigmaTPCPionDataV0tag[iBin][kPion]->Integral(-50,50)*intNsigmaTPCPionDataV0tag[iBin],eff,unc);
      hEffPionTPCDataV0tag[iEff]->SetBinContent(iBin+1,eff);
      hEffPionTPCDataV0tag[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(fNsigmaTPCKaonDataKinktag[iBin][kKaon]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCKaonDataKinktag[iBin],fNsigmaTPCKaonDataKinktag[iBin][kKaon]->Integral(-50,50)*intNsigmaTPCKaonDataKinktag[iBin],eff,unc);
      hEffKaonTPCDataKinktag[iEff]->SetBinContent(iBin+1,eff);
      hEffKaonTPCDataKinktag[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(fNsigmaTPCKaonDataTOFtag[iBin][kKaon]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCKaonDataTOFtag[iBin],fNsigmaTPCKaonDataTOFtag[iBin][kKaon]->Integral(-50,50)*intNsigmaTPCKaonDataTOFtag[iBin],eff,unc);
      hEffKaonTPCDataTOFtag[iEff]->SetBinContent(iBin+1,eff);
      hEffKaonTPCDataTOFtag[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(fNsigmaTPCProtonDataV0tag[iBin][kProton]->Integral(-nSigma[iEff],nSigma[iEff])*intNsigmaTPCProtonDataV0tag[iBin],fNsigmaTPCProtonDataV0tag[iBin][kProton]->Integral(-50,50)*intNsigmaTPCProtonDataV0tag[iBin],eff,unc);
      hEffProtonTPCDataV0tag[iEff]->SetBinContent(iBin+1,eff);
      hEffProtonTPCDataV0tag[iEff]->SetBinError(iBin+1,unc);
      
      ComputeEfficiency(hNsigmaTOFPionDataV0tag_sub[iBin]->Integral(nsigmabinlow,nsigmabinhigh)*intNsigmaTOFPionDataV0tag[iBin],hNsigmaTOFPionDataV0tag_sub[iBin]->Integral()*intNsigmaTOFPionDataV0tag[iBin],eff,unc);
      hEffPionTOFDataV0tag[iEff]->SetBinContent(iBin+1,eff);
      hEffPionTOFDataV0tag[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(hNsigmaTOFKaonDataTPCtag_sub[iBin]->Integral(nsigmabinlow,nsigmabinhigh)*intNsigmaTOFKaonDataTPCtag[iBin],hNsigmaTOFKaonDataTPCtag_sub[iBin]->Integral()*intNsigmaTOFKaonDataTPCtag[iBin],eff,unc);
      hEffKaonTOFDataTPCtag[iEff]->SetBinContent(iBin+1,eff);
      hEffKaonTOFDataTPCtag[iEff]->SetBinError(iBin+1,unc);

      ComputeEfficiency(hNsigmaTOFProtonDataV0tag_sub[iBin]->Integral(nsigmabinlow,nsigmabinhigh)*intNsigmaTOFProtonDataV0tag[iBin],hNsigmaTOFProtonDataV0tag_sub[iBin]->Integral()*intNsigmaTOFProtonDataV0tag[iBin],eff,unc);
      hEffProtonTOFDataV0tag[iEff]->SetBinContent(iBin+1,eff);
      hEffProtonTOFDataV0tag[iEff]->SetBinError(iBin+1,unc);
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
  cEffPion->cd(1)->DrawFrame(binlims[0],0.,binlims[nBins],1.,Form(";%s (GeV/#it{c});pion TPC PID efficiency",varname.Data()));
  cEffPion->cd(2)->DrawFrame(binlims[0],0.,binlims[nBins],1.,Form(";%s (GeV/#it{c});pion TOF PID efficiency",varname.Data()));
  cEffPion->cd(3)->DrawFrame(binlims[0],0.5,binlims[nBins],1.15,Form(";%s (GeV/#it{c}); pion data / MC TPC PID efficiency",varname.Data()));
  cEffPion->cd(4)->DrawFrame(binlims[0],0.5,binlims[nBins],1.15,Form(";%s (GeV/#it{c});pion data / MC TOF PID efficiency",varname.Data()));
  cEffKaon->Divide(2,2);
  cEffKaon->cd(1)->DrawFrame(binlims[0],0.,binlims[nBins],1.,Form(";%s (GeV/#it{c});kaon TPC PID efficiency",varname.Data()));
  cEffKaon->cd(2)->DrawFrame(binlims[0],0.,binlims[nBins],1.,Form(";%s (GeV/#it{c});kaon TOF PID efficiency",varname.Data()));
  cEffKaon->cd(3)->DrawFrame(binlims[0],0.5,binlims[nBins],1.15,Form(";%s (GeV/#it{c});kaon data / MC TPC PID efficiency",varname.Data()));
  cEffKaon->cd(4)->DrawFrame(binlims[0],0.5,binlims[nBins],1.15,Form(";%s (GeV/#it{c});kaon data / MC TOF PID efficiency",varname.Data()));
  cEffProton->Divide(2,2);
  cEffProton->cd(1)->DrawFrame(binlims[0],0.,binlims[nBins],1.,Form(";%s (GeV/#it{c});proton TPC PID efficiency",varname.Data()));
  cEffProton->cd(2)->DrawFrame(binlims[0],0.,binlims[nBins],1.,Form(";%s (GeV/#it{c});proton TOF PID efficiency",varname.Data()));
  cEffProton->cd(3)->DrawFrame(binlims[0],0.5,binlims[nBins],1.15,Form(";%s (GeV/#it{c});proton data / MC TPC PID efficiency",varname.Data()));
  cEffProton->cd(4)->DrawFrame(binlims[0],0.5,binlims[nBins],1.15,Form(";%s (GeV/#it{c});proton data / MC TOF PID efficiency",varname.Data()));

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

  //mean and width of Nsigma distributions
  TH1F* hMeanPionTPCMCV0tag = new TH1F("hMeanPionTPCMCV0tag",Form(";%s (GeV/#it{c});Mean(#pi)",varname.Data()),nBins,binlims);
  TH1F* hSigmaPionTPCMCV0tag = new TH1F("hSigmaPionTPCMCV0tag",Form(";%s (GeV/#it{c});Sigma(#pi)",varname.Data()),nBins,binlims);
  TH1F* hMeanPionTPCDataV0tag = new TH1F("hMeanPionTPCDataV0tag",Form(";%s (GeV/#it{c});Mean(#pi)",varname.Data()),nBins,binlims);
  TH1F* hSigmaPionTPCDataV0tag = new TH1F("hSigmaPionTPCDataV0tag",Form(";%s (GeV/#it{c});Sigma(#pi)",varname.Data()),nBins,binlims);
  SetTH1Style(hMeanPionTPCMCV0tag,markersMC[0],pdgcolors[kPion],1.,2,pdgcolors[kPion],kWhite,0.045,0.055);
  SetTH1Style(hSigmaPionTPCMCV0tag,markersMC[0],pdgcolors[kPion],1.,2,pdgcolors[kPion],kWhite,0.045,0.055);
  SetTH1Style(hMeanPionTPCDataV0tag,markersData[0],pdgcolors[kPion]+1,1.,2,pdgcolors[kPion]+1,kWhite,0.045,0.055);
  SetTH1Style(hSigmaPionTPCDataV0tag,markersData[0],pdgcolors[kPion]+1,1.,2,pdgcolors[kPion]+1,kWhite,0.045,0.055);

  TH1F* hMeanKaonTPCMCTOFtag = new TH1F("hMeanKaonTPCMCTOFtag",Form(";%s (GeV/#it{c});Mean(#pi)",varname.Data()),nBins,binlims);
  TH1F* hSigmaKaonTPCMCTOFtag = new TH1F("hSigmaKaonTPCMCTOFtag",Form(";%s (GeV/#it{c});Sigma(#pi)",varname.Data()),nBins,binlims);
  TH1F* hMeanKaonTPCDataTOFtag = new TH1F("hMeanKaonTPCDataTOFtag",Form(";%s (GeV/#it{c});Mean(#pi)",varname.Data()),nBins,binlims);
  TH1F* hSigmaKaonTPCDataTOFtag = new TH1F("hSigmaKaonTPCDataTOFtag",Form(";%s (GeV/#it{c});Sigma(#pi)",varname.Data()),nBins,binlims);
  SetTH1Style(hMeanKaonTPCMCTOFtag,markersMC[0],pdgcolors[kKaon],1.,2,pdgcolors[kKaon],kWhite,0.045,0.055);
  SetTH1Style(hSigmaKaonTPCMCTOFtag,markersMC[0],pdgcolors[kKaon],1.,2,pdgcolors[kKaon],kWhite,0.045,0.055);
  SetTH1Style(hMeanKaonTPCDataTOFtag,markersData[0],pdgcolors[kKaon]+1,1.,2,pdgcolors[kKaon]+1,kWhite,0.045,0.055);
  SetTH1Style(hSigmaKaonTPCDataTOFtag,markersData[0],pdgcolors[kKaon]+1,1.,2,pdgcolors[kKaon]+1,kWhite,0.045,0.055);

  TH1F* hMeanProtonTPCMCV0tag = new TH1F("hMeanProtonTPCMCV0tag",Form(";%s (GeV/#it{c});Mean(#pi)",varname.Data()),nBins,binlims);
  TH1F* hSigmaProtonTPCMCV0tag = new TH1F("hSigmaProtonTPCMCV0tag",Form(";%s (GeV/#it{c});Sigma(#pi)",varname.Data()),nBins,binlims);
  TH1F* hMeanProtonTPCDataV0tag = new TH1F("hMeanProtonTPCDataV0tag",Form(";%s (GeV/#it{c});Mean(#pi)",varname.Data()),nBins,binlims);
  TH1F* hSigmaProtonTPCDataV0tag = new TH1F("hSigmaProtonTPCDataV0tag",Form(";%s (GeV/#it{c});Sigma(#pi)",varname.Data()),nBins,binlims);
  SetTH1Style(hMeanProtonTPCMCV0tag,markersMC[0],pdgcolors[kProton],1.,2,pdgcolors[kProton],kWhite,0.045,0.055);
  SetTH1Style(hSigmaProtonTPCMCV0tag,markersMC[0],pdgcolors[kProton],1.,2,pdgcolors[kProton],kWhite,0.045,0.055);
  SetTH1Style(hMeanProtonTPCDataV0tag,markersData[0],pdgcolors[kProton]+1,1.,2,pdgcolors[kProton]+1,kWhite,0.045,0.055);
  SetTH1Style(hSigmaProtonTPCDataV0tag,markersData[0],pdgcolors[kProton]+1,1.,2,pdgcolors[kProton]+1,kWhite,0.045,0.055);

  TLegend* legPionPars = new TLegend(0.8,0.7,0.99,0.89);
  legPionPars->SetTextSize(0.045);
  legPionPars->AddEntry(hMeanPionTPCMCV0tag,"MC","p");
  legPionPars->AddEntry(hMeanPionTPCDataV0tag,"V0 tag","p");
  TLegend* legKaonPars = new TLegend(0.8,0.7,0.99,0.89);
  legKaonPars->SetTextSize(0.045);
  legKaonPars->AddEntry(hMeanKaonTPCMCTOFtag,"MC","p");
  legKaonPars->AddEntry(hMeanKaonTPCDataTOFtag,"TOF tag","p");
  TLegend* legProtonPars = new TLegend(0.8,0.7,0.99,0.89);
  legProtonPars->SetTextSize(0.045);
  legProtonPars->AddEntry(hMeanProtonTPCMCV0tag,"MC","p");
  legProtonPars->AddEntry(hMeanProtonTPCDataV0tag,"V0 tag","p");

  for(int iBin=0; iBin<nBins; iBin++) {
    hMeanPionTPCMCV0tag->SetBinContent(iBin+1,fNsigmaTPCPionMCV0tag[iBin][kPion]->GetParameter(1));
    hMeanPionTPCMCV0tag->SetBinError(iBin+1,1.e-20);
    hSigmaPionTPCMCV0tag->SetBinContent(iBin+1,fNsigmaTPCPionMCV0tag[iBin][kPion]->GetParameter(2));
    hSigmaPionTPCMCV0tag->SetBinError(iBin+1,1.e-20);
    hMeanPionTPCDataV0tag->SetBinContent(iBin+1,fNsigmaTPCPionDataV0tag[iBin][kPion]->GetParameter(1));
    hMeanPionTPCDataV0tag->SetBinError(iBin+1,1.e-20);
    hSigmaPionTPCDataV0tag->SetBinContent(iBin+1,fNsigmaTPCPionDataV0tag[iBin][kPion]->GetParameter(2));
    hSigmaPionTPCDataV0tag->SetBinError(iBin+1,1.e-20);

    hMeanKaonTPCMCTOFtag->SetBinContent(iBin+1,fNsigmaTPCKaonMCTOFtag[iBin][kKaon]->GetParameter(1));
    hMeanKaonTPCMCTOFtag->SetBinError(iBin+1,1.e-20);
    hSigmaKaonTPCMCTOFtag->SetBinContent(iBin+1,fNsigmaTPCKaonMCTOFtag[iBin][kKaon]->GetParameter(2));
    hSigmaKaonTPCMCTOFtag->SetBinError(iBin+1,1.e-20);
    hMeanKaonTPCDataTOFtag->SetBinContent(iBin+1,fNsigmaTPCKaonDataTOFtag[iBin][kKaon]->GetParameter(1));
    hMeanKaonTPCDataTOFtag->SetBinError(iBin+1,1.e-20);
    hSigmaKaonTPCDataTOFtag->SetBinContent(iBin+1,fNsigmaTPCKaonDataTOFtag[iBin][kKaon]->GetParameter(2));
    hSigmaKaonTPCDataTOFtag->SetBinError(iBin+1,1.e-20);

    hMeanProtonTPCMCV0tag->SetBinContent(iBin+1,fNsigmaTPCProtonMCV0tag[iBin][kProton]->GetParameter(1));
    hMeanProtonTPCMCV0tag->SetBinError(iBin+1,1.e-20);
    hSigmaProtonTPCMCV0tag->SetBinContent(iBin+1,fNsigmaTPCProtonMCV0tag[iBin][kProton]->GetParameter(2));
    hSigmaProtonTPCMCV0tag->SetBinError(iBin+1,1.e-20);
    hMeanProtonTPCDataV0tag->SetBinContent(iBin+1,fNsigmaTPCProtonDataV0tag[iBin][kProton]->GetParameter(1));
    hMeanProtonTPCDataV0tag->SetBinError(iBin+1,1.e-20);
    hSigmaProtonTPCDataV0tag->SetBinContent(iBin+1,fNsigmaTPCProtonDataV0tag[iBin][kProton]->GetParameter(2));
    hSigmaProtonTPCDataV0tag->SetBinError(iBin+1,1.e-20);
  }

  TCanvas* cMeanSigma = new TCanvas("cMeanSigma","",1920,1080);
  cMeanSigma->Divide(3,2);
  cMeanSigma->cd(1)->DrawFrame(binlims[0],-3.,binlims[nBins],3.,Form(";%s (GeV/#it{c});mean #it{N}_{#sigma}(#pi)",varname.Data()));
  hMeanPionTPCMCV0tag->Draw("same");
  hMeanPionTPCDataV0tag->Draw("same");
  legPionPars->Draw();
  cMeanSigma->cd(2)->DrawFrame(binlims[0],-3.,binlims[nBins],3.,Form(";%s (GeV/#it{c});mean #it{N}_{#sigma}(K)",varname.Data()));
  hMeanKaonTPCMCTOFtag->Draw("same");
  hMeanKaonTPCDataTOFtag->Draw("same");
  legKaonPars->Draw();
  cMeanSigma->cd(3)->DrawFrame(binlims[0],-3.,binlims[nBins],3.,Form(";%s (GeV/#it{c});mean #it{N}_{#sigma}(p)",varname.Data()));
  hMeanProtonTPCMCV0tag->Draw("same");
  hMeanProtonTPCDataV0tag->Draw("same");
  legProtonPars->Draw();
  cMeanSigma->cd(4)->DrawFrame(binlims[0],0.,binlims[nBins],3.,Form(";%s (GeV/#it{c});width #it{N}_{#sigma}(#pi)",varname.Data()));
  hSigmaPionTPCMCV0tag->Draw("same");
  hSigmaPionTPCDataV0tag->Draw("same");
  legPionPars->Draw();
  cMeanSigma->cd(5)->DrawFrame(binlims[0],0.,binlims[nBins],3.,Form(";%s (GeV/#it{c});width #it{N}_{#sigma}(K)",varname.Data()));
  hSigmaKaonTPCMCTOFtag->Draw("same");
  hSigmaKaonTPCDataTOFtag->Draw("same");
  legKaonPars->Draw();
  cMeanSigma->cd(6)->DrawFrame(binlims[0],0.,binlims[nBins],3.,Form(";%s (GeV/#it{c});width #it{N}_{#sigma}(p)",varname.Data()));
  hSigmaProtonTPCMCV0tag->Draw("same");
  hSigmaProtonTPCDataV0tag->Draw("same");
  legProtonPars->Draw();

  //********************************************************************************************************************************************//
  //output files
  TFile outfile(Form("%s/PIDSystSingleTrack.root",outputdirName.Data()),"RECREATE");
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

  TFile outfiledist(Form("%s/NsigmaPIDdist.root",outputdirName.Data()),"RECREATE");
  hMeanPionTPCMCV0tag->Write();
  hMeanPionTPCDataV0tag->Write();
  hMeanKaonTPCMCTOFtag->Write();
  hMeanKaonTPCDataTOFtag->Write();
  hMeanProtonTPCMCV0tag->Write();
  hMeanProtonTPCDataV0tag->Write();
  hSigmaPionTPCMCV0tag->Write();
  hSigmaPionTPCDataV0tag->Write();
  hSigmaKaonTPCMCTOFtag->Write();
  hSigmaKaonTPCDataTOFtag->Write();
  hSigmaProtonTPCMCV0tag->Write();
  hSigmaProtonTPCDataV0tag->Write();

  for(int iBin=0; iBin<nBins; iBin++) {
    hNsigmaTPCPionMCTrue[iBin]->Write();
    hNsigmaTPCKaonMCTrue[iBin]->Write();
    hNsigmaTPCProtonMCTrue[iBin]->Write();
    hNsigmaTOFPionMCTrue[iBin]->Write();
    hNsigmaTOFKaonMCTrue[iBin]->Write();
    hNsigmaTOFProtonMCTrue[iBin]->Write();
  }
  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      hNsigmaTPCPionMCV0tag[iBin][iPart]->Write();
      hNsigmaTPCKaonMCKinktag[iBin][iPart]->Write();
      hNsigmaTPCKaonMCTOFtag[iBin][iPart]->Write();
      hNsigmaTPCProtonMCV0tag[iBin][iPart]->Write();
      hNsigmaTOFPionMCV0tag[iBin][iPart]->Write();
      hNsigmaTOFKaonMCKinktag[iBin][iPart]->Write();
      hNsigmaTOFKaonMCTPCtag[iBin][iPart]->Write();
      hNsigmaTOFProtonMCV0tag[iBin][iPart]->Write();

    }
  }
  for(int iBin=0; iBin<nBins; iBin++) {
    hNsigmaTPCPionDataV0tag[iBin]->Write();
    hNsigmaTPCKaonDataKinktag[iBin]->Write();
    hNsigmaTPCKaonDataTOFtag[iBin]->Write();
    hNsigmaTPCProtonDataV0tag[iBin]->Write();
    hNsigmaTOFPionDataV0tag[iBin]->Write();
    hNsigmaTOFKaonDataKinktag[iBin]->Write();
    hNsigmaTOFKaonDataTPCtag[iBin]->Write();
    hNsigmaTOFProtonDataV0tag[iBin]->Write();
  }
  for(int iBin=0; iBin<nBins; iBin++) {
    hNsigmaTOFPionDataV0tag_sub[iBin]->Write();
    hNsigmaTOFKaonDataKinktag_sub[iBin]->Write();
    hNsigmaTOFKaonDataTPCtag_sub[iBin]->Write();
    hNsigmaTOFProtonDataV0tag_sub[iBin]->Write();
  }
  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      fNsigmaTPCPionMCV0tag[iBin][iPart]->Write();
      fNsigmaTPCKaonMCKinktag[iBin][iPart]->Write();
      fNsigmaTPCKaonMCTOFtag[iBin][iPart]->Write();
      fNsigmaTPCProtonMCV0tag[iBin][iPart]->Write();
    }
  }
  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      fNsigmaTPCPionDataV0tag[iBin][iPart]->Write();
      fNsigmaTPCKaonDataKinktag[iBin][iPart]->Write();
      fNsigmaTPCKaonDataTOFtag[iBin][iPart]->Write();
      fNsigmaTPCProtonDataV0tag[iBin][iPart]->Write();
    }
  }
  for(int iBin=0; iBin<nBins; iBin++) {
    for(int iPart=kElectron; iPart<=kAll; iPart++) {
      if(hNsigmaTOFKaonDataTPCtag_Fit[iBin][iPart]) hNsigmaTOFKaonDataTPCtag_Fit[iBin][iPart]->Write();
      if(hNsigmaTOFProtonDataV0tag_Fit[iBin][iPart]) hNsigmaTOFProtonDataV0tag_Fit[iBin][iPart]->Write();
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

  cFractionTPCMC->SaveAs(Form("%s/ParticleFractionsTPCMC.pdf",outputdirName.Data()));
  cFractionTOFMC->SaveAs(Form("%s/ParticleFractionsTOFMC.pdf",outputdirName.Data()));
  cTOFFractionData->SaveAs(Form("%s/ParticleFractionsTOFdata.pdf",outputdirName.Data()));

  cFitResultTOFKaonFromTPCtag->SaveAs(Form("%s/TFractionFitterResultTOFKaonFromTPCtag.pdf",outputdirName.Data()));
  cFitResultTOFProtonFromV0tag->SaveAs(Form("%s/TFractionFitterResultTOFProtonFromV0tag.pdf",outputdirName.Data()));

  cPionMCV0tagTPC->SaveAs(Form("%s/PionMCV0tagTPC.pdf",outputdirName.Data()));
  cKaonMCKinkstagTPC->SaveAs(Form("%s/KaonMCKinkstagTPC.pdf",outputdirName.Data()));
  cKaonMCTOFtagTPC->SaveAs(Form("%s/KaonMCTOFtagTPC.pdf",outputdirName.Data()));
  cProtonMCV0tagTPC->SaveAs(Form("%s/ProtonMCV0tagTPC.pdf",outputdirName.Data()));
  cPionMCV0tagTOF->SaveAs(Form("%s/PionMCV0tagTOF.pdf",outputdirName.Data()));
  cKaonMCKinkstagTOF->SaveAs(Form("%s/KaonMCKinkstagTOF.pdf",outputdirName.Data()));
  cKaonMCTPCtagTOF->SaveAs(Form("%s/KaonMCTPCtagTOF.pdf",outputdirName.Data()));
  cProtonMCV0tagTOF->SaveAs(Form("%s/ProtonMCV0tagTOF.pdf",outputdirName.Data()));

  cPionDataV0tagTPC->SaveAs(Form("%s/PionDataV0tagTPC.pdf",outputdirName.Data()));
  cKaonDataKinkstagTPC->SaveAs(Form("%s/KaonDataKinkstagTPC.pdf",outputdirName.Data()));
  cKaonDataTOFtagTPC->SaveAs(Form("%s/KaonDataTOFtagTPC.pdf",outputdirName.Data()));
  cProtonDataV0tagTPC->SaveAs(Form("%s/ProtonDataV0tagTPC.pdf",outputdirName.Data()));
  cPionDataV0tagTOF->SaveAs(Form("%s/PionDataV0tagTOF.pdf",outputdirName.Data()));
  cKaonDataKinkstagTOF->SaveAs(Form("%s/KaonDataKinkstagTOF.pdf",outputdirName.Data()));
  cKaonDataTPCtagTOF->SaveAs(Form("%s/KaonDataTPCtagTOF.pdf",outputdirName.Data()));
  cProtonDataV0tagTOF->SaveAs(Form("%s/ProtonDataV0tagTOF.pdf",outputdirName.Data()));

  cEffPion->SaveAs(Form("%s/PionPIDefficiency_data_MC.pdf",outputdirName.Data()));
  cEffKaon->SaveAs(Form("%s/KaonPIDefficiency_data_MC.pdf",outputdirName.Data()));
  cEffProton->SaveAs(Form("%s/ProtonPIDefficiency_data_MC.pdf",outputdirName.Data()));

  cMeanSigma->SaveAs(Form("%s/NsigmaTPCDistrPars_data_MC.pdf",outputdirName.Data()));
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
int FindPtbin(float pt, const double binlims[], int nBins) {

  int ptbin=-1;
  for(int iBin=0; iBin<nBins; iBin++) {
    if(pt>binlims[iBin] && pt<binlims[iBin+1]) {
      ptbin=iBin;
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
void GetTOFFractionsFromData(int whichpart, int iBin, TH1F* hFractionMC[nPDGcodes-1], TH1F* hFractionData[nPDGcodes-1], TH1F* hNsigmaMC[nPDGcodes], TH1F* hNsigmaData, TFractionFitter *&fNsigmaFitter, vector<int> &templUsed) {

  TObjArray* oNsigmaMC = new TObjArray(0);
    
    for(int iPart=kElectron; iPart<=kProton; iPart++) {
      if(iPart==whichpart || hFractionMC[iPart]->GetBinContent(iBin+1)>1.e-5) {
        TH1F* hMC = (TH1F*)hNsigmaMC[iPart]->Clone(Form("hMC_%d_%s",iBin,pdgnames[iPart].Data()));
        oNsigmaMC->AddLast(hMC);
        templUsed.push_back(iPart);
      }
    }
    
    if(templUsed.size()==1) {
      for(int iPart=kElectron; iPart<=kProton; iPart++) {
        if(iPart==whichpart) {
          hFractionData[iPart]->SetBinContent(iBin+1,1.);
        }
       else {
          hFractionData[iPart]->SetBinContent(iBin+1,0.);
        }
      }
      return;
    }

    if(oNsigmaMC->GetEntries()>1) {
      
      fNsigmaFitter = new TFractionFitter(hNsigmaData,oNsigmaMC);
      for(int iEntry=0; iEntry<oNsigmaMC->GetEntries(); iEntry++) {
        fNsigmaFitter->Constrain(iEntry,hFractionMC[templUsed[iEntry]]->GetBinContent(iBin+1)*0.5,hFractionMC[templUsed[iEntry]]->GetBinContent(iBin+1)*2);
      }
      fNsigmaFitter->Fit();
      
      for(int iPart=kElectron; iPart<=kProton; iPart++) {

        vector<int>::iterator it = find(templUsed.begin(), templUsed.end(), iPart);
        if(it!=templUsed.end()) {
          double frac, err;
          int iEntry = static_cast<int>(distance(templUsed.begin(),it));

          fNsigmaFitter->GetResult(iEntry, frac, err);
          hFractionData[iPart]->SetBinContent(iBin+1,frac);
          hFractionData[iPart]->SetBinError(iBin+1,err);
        }
        else {
          hFractionData[iPart]->SetBinContent(iBin+1,0);
          hFractionData[iPart]->SetBinError(iBin+1,0);
        }
      }
    }
    else {
      for(int iPart=kElectron; iPart<=kProton; iPart++) {
        if(iPart!=whichpart) {
          hFractionData[iPart]->SetBinContent(iBin+1,1);
          hFractionData[iPart]->SetBinError(iBin+1,1);
        }
        else {
          hFractionData[iPart]->SetBinContent(iBin+1,0);
          hFractionData[iPart]->SetBinError(iBin+1,0);
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
    hArmenterosData[iHisto]->SetMarkerStyle(kFullCircle);
    hArmenterosMC[iHisto]->SetMarkerStyle(kFullCircle);
    hArmenterosData[iHisto]->SetMarkerSize(0.5);
    hArmenterosMC[iHisto]->SetMarkerSize(0.5);
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

  cArmentero->SaveAs(Form("%s/V0QAplots.pdf",outputdirName.Data()));
  cKinks->SaveAs(Form("%s/KinksQAplots.pdf",outputdirName.Data()));
  cKinksMass->SaveAs(Form("%s/KinksMass.pdf",outputdirName.Data()));
  cArmentero->SaveAs(Form("%s/V0QAplots.png",outputdirName.Data()));
  cKinks->SaveAs(Form("%s/KinksQAplots.png",outputdirName.Data()));
}

//_____________________________________________________
void DivideCanvas(TCanvas* c, int nBins) {
  if(nBins<2)
    c->cd();
  else if(nBins==2 || nBins==3)
    c->Divide(nBins,1);
  else if(nBins==4 || nBins==6 || nBins==8)
    c->Divide(nBins/2,2);
  else if(nBins==5 || nBins==7)
    c->Divide((nBins+1)/2,2);
  else if(nBins==9 || nBins==12 || nBins==15)
    c->Divide(nBins/3,3);
  else if(nBins==10 || nBins==11)
    c->Divide(4,3);
  else if(nBins==13 || nBins==14)
    c->Divide(5,3);
  else if(nBins>15 && nBins<=20 && nBins%4==0)
    c->Divide(nBins/4,4);
  else if(nBins>15 && nBins<=20 && nBins%4!=0)
    c->Divide(5,4);
  else if(nBins==21)
    c->Divide(7,3);
  else if(nBins>21 && nBins<=25)
    c->Divide(5,5);
  else if(nBins>25 && nBins%2==0)
    c->Divide(nBins/2,2);
  else
    c->Divide((nBins+1)/2,2);
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
