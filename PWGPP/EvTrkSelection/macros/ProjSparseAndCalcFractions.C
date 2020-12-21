#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TObjArray.h"

#include <Riostream.h>
#include <fstream>
#include <string>
using std::cout;
using std::endl;


/////////////////////////////////////////////////////////////////////////////////////
//  This scirpt produces output txt files containg values for primary particles    //
//  fraction in data.                                                              //
//  How to use:                                                                    //
//  .L ProjSparseAndCalcFractions.C+                                               //
//  1) projectSparse() *** to be called once: it produces thnsparse projections    //
//  and obtains DCA distributions vs desired variables (after desired cuts)        //
//  Then saves histograms in a root file.                                          //
//  arg1: root file with thnsparse on data                                         //
//  arg2: root file with thnsparse on MC                                           //
//  arg3: set the suffix containg your production/period name                      //
//  arg4: fWhichVar == 1 vs pt, 2 vs phi, 3 vs eta                                 //
//  arg5: fUseMCPt: true = particle MC pt used, false = track global pt used       //
//  2) calcFractions() to extract fraction through a TFractionFitter.              //
//  You can launch this function multiple times in order to adjust fit             //
//  parametrs and obtains primary fractions in all your analysis bins.             //
//  arg1: fWhichVar == 1 vs pt, 2 vs phi, 3 vs eta                                 //
//  arg2: fUseMCPt: true = particle MC pt used, false = track global pt used       //
//  arg3: set the suffix containg your production/period name                      //
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Set all the variables here /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

TString suffix="_PiK"; //suffix to match with that of thnsparse you want to project
Bool_t fUseTOFbc = kFALSE;

const int nPhi = 18; // nbins for vs phi analysis
const int nEta = 16; // nbins for vs eta analysis
const int npt  = 7;  // nbins for vs pt analysis

double etaRange = 0.8;  //set eta range
double minPhi   = 0.;   //set min phi
double maxPhi   = 6.28; //set max phi

double ptlims[npt+1]   = {0.7,1.,2.,3.,4.,6.,8.,15.};
double xlowsec[nPhi]   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double xhighprim[nPhi] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

//Usually, it is enough to adjust these lines below: these are the constraints for the
//parameters of the TFractionFitter, namely primary fraction lower limit and
//secondary fraction upper limit.
//***** if vs pt analysis
double xlowprimVsPt[npt]     = {0.9,0.9,0.9,0.9,0.9,0.9,0.9};
double xhighsecVsPt[npt]     = {0.2,0.2,0.2,0.2,0.1,0.1,0.1};
double xhighsecMatVsPt[npt]  = {0.10,0.1,0.1,0.05,0.05,0.02,0.02};

//***** if vs phi
double xlowprimVsPhi[nPhi]  = {0.6,0.9,0.9,0.9,0.9,0.9,0.5,0.7,0.9,0.8,0.8,0.9,0.9,0.9,0.9,0.9,0.5,0.7};
double xhighsecVsPhi[nPhi]  = {0.3,0.1,0.5,0.5,0.5,0.5,0.5,0.3,0.3,0.3,0.3,0.2,0.3,0.3,0.3,0.3,0.5,0.1};
double xhighsecMatVsPhi[nPhi]  = {0.3,0.1,0.5,0.5,0.5,0.5,0.5,0.3,0.3,0.3,0.3,0.2,0.3,0.3,0.3,0.3,0.5,0.1};

//***** if vs eta
double xlowprimVsEta[nEta]  = {0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95};
double xhighsecVsEta[nEta]  = {0.3,0.2,0.3,0.2,0.3,0.3,0.2,0.3,0.3,0.3,0.2,0.2,0.3,0.2,0.3,0.3};
double xhighsecMatVsEta[nEta]  = {0.3,0.2,0.3,0.2,0.3,0.3,0.2,0.3,0.3,0.3,0.2,0.2,0.3,0.2,0.3,0.3};

/////////////Don't touch the 3 lines below (no of axis of thnsparses)/////////////
const Int_t nVarData   = 6;  //nVars for data
const Int_t nVarMC     = 10;  //nVars for MC
const Int_t nVarMCTPC  = 10;  //nVars for MC, TPC only

//do not touch unless you want fraction in a smaller DCAxy range
double fitRange[npt] = {2.4,2.4,2.4,2.4,2.4,2.4,2.4};

//////////////////////////////////////////////////////////////////////////////////

Bool_t fUseMCPt;    //true = particle MC pt used, false = track global pt used
Int_t  fWhichVar;   //1 == vs pt,2 == vs phi,3 == vs eta

TH1D* projectMC(THnSparse *hSparse, Double_t partType, Int_t ipt);
TH1D* projectData(THnSparse *hSparse, Int_t);
void calcFractions(Int_t , Int_t, TString, TString);

//___________________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
void projectSparse(TString filenameData = "15o/AnalysisResults_Data_megarun.root",
                   TString filenameMC   = "15o/AnalysisResults_MC_megarun.root",
                   TString period = "_highIR",
                   TString dir = "15o",
                   Int_t vsWhichVar = 3,
                   Int_t fPt = 1) {
    
    vector<double> vectorPrim;
    vector<double> vectorSec;
    vectorPrim.clear();
    vectorSec.clear();
    
    TString whichPt;
    fUseMCPt   = kFALSE;
    Int_t npx, npy;
    npx = 5;
    npy = 2;
    
    if(fPt==1) {
        fUseMCPt = kFALSE;
        whichPt = "globPt";
    }
    else if(fPt==2) {
        fUseMCPt = kTRUE;
        whichPt  = "MCPt";
    }
    
    TFile *_file0 = TFile::Open(filenameMC.Data());
    TFile *_file1 = TFile::Open(filenameData.Data());
    if(!_file0 || !_file1) {
        Printf("Input files for MC or data not found!");
        return;
    }
    TList *list0  = (TList*)_file0->Get(Form("trackingUncert%s",suffix.Data()));
    TList *list1  = (TList*)_file1->Get(Form("trackingUncert%s",suffix.Data()));
    if(!list0 || !list1) {
        Printf("Lists not found!");
        return;
    }
    
    THnSparse* hSparseMC    = (THnSparse*)list0->FindObject("fHistMC");
    THnSparse* hSparseMCTPC = (THnSparse*)list0->FindObject("fHistMCTPConly");
    THnSparse* hSparseData  = (THnSparse*)list1->FindObject("fHistData");
    if(!hSparseMC || !hSparseMCTPC || !hSparseData) {
        Printf("Sparses not found!");
        return;
    }
    TH1F *hData;
    TH1F *hMC[3], *hMCTPC[3];
    TH1F *hPrimMC, *hSecMC, *hSecMatMC, *hTotMC, *hPrimMCTPConly, *hSecMCTPConly, *hSecMatMCTPConly, *hTotMCTPConly;
    
    fWhichVar = vsWhichVar;
    TString whichVar;
    Int_t nBinTot;
    
    if(vsWhichVar == 1) {
        printf("**** Analysis versus Pt ****\n");
        
        whichVar = "Pt";
        nBinTot = npt;
        
        hPrimMC           = new TH1F("hPrimMC","",npt,ptlims);
        hSecMC            = new TH1F("hSecMC","",npt,ptlims);
        hSecMatMC         = new TH1F("hSecMatMC","",npt,ptlims);
        hTotMC            = new TH1F("hTotMC","",npt,ptlims);
        hPrimMCTPConly    = new TH1F("hPrimMCTPConly","",npt,ptlims);
        hSecMCTPConly     = new TH1F("hSecMCTPConly","",npt,ptlims);
        hSecMatMCTPConly  = new TH1F("hSecMatMCTPConly","",npt,ptlims);
        hTotMCTPConly     = new TH1F("hTotMCTPConly","",npt,ptlims);
        
    }
    else if(vsWhichVar == 2) {
        printf("**** Analysis versus Phi ****\n");
        
        nBinTot  = nPhi;
        whichVar = "Phi";
        
        hPrimMC           = new TH1F("hPrimMC","",nPhi,minPhi,maxPhi);
        hSecMC            = new TH1F("hSecMC","",nPhi,minPhi,maxPhi);
        hSecMatMC         = new TH1F("hSecMatMC","",nPhi,minPhi,maxPhi);
        hTotMC            = new TH1F("hTotMC","",nPhi,minPhi,maxPhi);
        hPrimMCTPConly    = new TH1F("hPrimMCTPConly","",nPhi,minPhi,maxPhi);
        hSecMCTPConly     = new TH1F("hSecMCTPConly","",nPhi,minPhi,maxPhi);
        hSecMatMCTPConly  = new TH1F("hSecMatMCTPConly","",nPhi,minPhi,maxPhi);
        hTotMCTPConly     = new TH1F("hTotMCTPConly","",nPhi,minPhi,maxPhi);
        
        npx = 5;
        npy = 4;
    }
    else if(vsWhichVar == 3) {
        printf("**** Analysis versus Eta ****\n");
        
        nBinTot  = nEta;
        whichVar = "Eta";
        
        hPrimMC           = new TH1F("hPrimMC","",nEta,-etaRange,etaRange);
        hSecMC            = new TH1F("hSecMC","",nEta,-etaRange,etaRange);
        hSecMatMC         = new TH1F("hSecMatMC","",nEta,-etaRange,etaRange);
        hTotMC            = new TH1F("hTotMC","",nEta,-etaRange,etaRange);
        hPrimMCTPConly    = new TH1F("hPrimMCTPConly","",nEta,-etaRange,etaRange);
        hSecMCTPConly     = new TH1F("hSecMCTPConly","",nEta,-etaRange,etaRange);
        hSecMatMCTPConly  = new TH1F("hSecMatMCTPConly","",nEta,-etaRange,etaRange);
        hTotMCTPConly     = new TH1F("hTotMCTPConly","",nEta,-etaRange,etaRange);
        
        npx = 4;
        npy = 4;
    }
    else {
        Printf("Invalid setting for analysis variable!");
        return;
    }
    
    TFile *fout = new TFile(Form("%s/Proj%s_%s_%s_%s.root",dir.Data(),suffix.Data(),whichPt.Data(),period.Data(),whichVar.Data()),"recreate");
    TCanvas *c1 = new TCanvas("c1","",1300,700);
    c1->Divide(npx,npy);
    

    Double_t FracTPC[4];
    Double_t Frac[4];
    for(int ibin = 0; ibin < nBinTot; ibin++) {
 
        Printf("\n*** Bin %d ***",ibin);
        Frac[3]    = 0.;
        FracTPC[3] = 0.;
        for(int iSpecie = 0; iSpecie < 3; iSpecie++) {
            hMC[iSpecie]    = (TH1F*)projectMC(hSparseMC,(double)iSpecie,ibin); //prim
            hMCTPC[iSpecie] = (TH1F*)projectMC(hSparseMCTPC,(double)iSpecie,ibin); //prim
            
            Frac[iSpecie]    = hMC[iSpecie]->GetEntries();
            FracTPC[iSpecie] = hMCTPC[iSpecie]->GetEntries();
       
            Frac[3]    += Frac[iSpecie];
            FracTPC[3] += FracTPC[iSpecie];
        }
        
        hData = (TH1F*)projectData(hSparseData,ibin);
        
        hPrimMC->SetBinContent(ibin+1,Frac[0]);
        hSecMC->SetBinContent(ibin+1,Frac[1]);
        hSecMatMC->SetBinContent(ibin+1,Frac[2]);
        hTotMC->SetBinContent(ibin+1,Frac[3]);
        hPrimMCTPConly->SetBinContent(ibin+1,FracTPC[0]);
        hSecMCTPConly->SetBinContent(ibin+1,FracTPC[1]);
        hSecMatMCTPConly->SetBinContent(ibin+1,FracTPC[2]);
        hTotMCTPConly->SetBinContent(ibin+1,FracTPC[3]);
        
        hPrimMC->SetBinError(ibin+1,TMath::Sqrt(Frac[0]));
        hSecMC->SetBinError(ibin+1,TMath::Sqrt(Frac[1]));
        hSecMatMC->SetBinError(ibin+1,TMath::Sqrt(Frac[2]));
        hTotMC->SetBinError(ibin+1,TMath::Sqrt(Frac[3]));
        hPrimMCTPConly->SetBinError(ibin+1,TMath::Sqrt(FracTPC[0]));
        hSecMCTPConly->SetBinError(ibin+1,TMath::Sqrt(FracTPC[1]));
        hSecMatMCTPConly->SetBinError(ibin+1,TMath::Sqrt(FracTPC[2]));
        hTotMCTPConly->SetBinError(ibin+1,TMath::Sqrt(FracTPC[3]));
        
        
        hMC[0]->Write(Form("hMC%s%d_Prim_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        hMC[1]->Write(Form("hMC%s%d_Sec_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        hMC[2]->Write(Form("hMC%s%d_Mat_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        hData->Write(Form("hData%s%d_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        
        c1->cd(ibin+1);
        gPad->SetLogy();
        hMC[0]->Add(hMC[1]);
        hMC[0]->Add(hMC[2]);
        hMC[0]->SetLineColor(kRed);
        hData->SetLineColor(kBlue);
        hMC[0]->Scale(hData->Integral()/hMC[0]->Integral());
        hMC[0]->DrawCopy();
        hData->DrawCopy("same");
    }
    fout->Close();
    
    
    
    c1->SaveAs(Form("%s/DCA_%s_ESDTrOnly_%s_Vs%s.eps",dir.Data(),suffix.Data(),period.Data(),whichVar.Data()));
    
    
    hPrimMC->Sumw2();
    hSecMC->Sumw2();
    hSecMatMC->Sumw2();
    hTotMC->Sumw2();
    hPrimMCTPConly->Sumw2();
    hSecMCTPConly->Sumw2();
    hSecMatMCTPConly->Sumw2();
    hTotMCTPConly->Sumw2();
    
    hPrimMC->Divide(hPrimMC,hTotMC,1,1,"B");
    hSecMC->Divide(hSecMC,hTotMC,1,1,"B");
    hSecMatMC->Divide(hSecMatMC,hTotMC,1,1,"B");
    hPrimMCTPConly->Divide(hPrimMCTPConly,hTotMCTPConly,1,1,"B");
    hSecMCTPConly->Divide(hSecMCTPConly,hTotMCTPConly,1,1,"B");
    hSecMatMCTPConly->Divide(hSecMatMCTPConly,hTotMCTPConly,1,1,"B");
    
    TCanvas *c = new TCanvas("c","",800,700);
    gPad->SetLogy();
    
    hPrimMC->SetLineColor(kBlue);
    hSecMC->SetLineColor(kRed);
    hPrimMCTPConly->SetLineColor(kBlue);
    hSecMCTPConly->SetLineColor(kRed);
    hPrimMCTPConly->SetLineStyle(kDashed);
    hSecMCTPConly->SetLineStyle(kDashed);
    hSecMatMCTPConly->SetLineStyle(kDashed);
    hSecMatMCTPConly->SetLineColor(kMagenta);
    hSecMatMC->SetLineColor(kMagenta);
    hPrimMC->GetYaxis()->SetRangeUser(0.0005,1.2);
    hPrimMC->GetYaxis()->SetTitle("MC fractions");
    if(vsWhichVar==1) hPrimMC->GetXaxis()->SetTitle("p_{T}");
    else if(vsWhichVar==2) hPrimMC->GetXaxis()->SetTitle("phi");
    else if(vsWhichVar==3) hPrimMC->GetXaxis()->SetTitle("eta");
    hPrimMC->SetStats(0);
    hPrimMC->DrawCopy();
    hPrimMCTPConly->DrawCopy("same");
    hSecMC->DrawCopy("same");
    hSecMCTPConly->DrawCopy("same");
    hSecMatMCTPConly->DrawCopy("same");
    hSecMatMC->DrawCopy("same");
    
    TLegend *t = new TLegend(0.65,0.65,0.88,0.85);
    t->AddEntry(hPrimMC,"prim.");
    t->AddEntry(hSecMC,"sec.");
    t->AddEntry(hSecMatMC,"sec. mat");
    t->AddEntry(hPrimMCTPConly,"prim. TPC only");
    t->AddEntry(hSecMCTPConly,"sec. TPC only");
    t->AddEntry(hSecMatMCTPConly,"sec. mat. TPC only");
    t->Draw("same");
    
    c->SaveAs(Form("%s/MCfractions_ESDTrOnly_Vs%s%s.eps",dir.Data(),whichVar.Data(),suffix.Data()));
    
    TCanvas *cl = new TCanvas("cl","",700,600);
    
    hPrimMCTPConly->Sumw2();
    hPrimMC->Sumw2();
    hPrimMCTPConly->Divide(hPrimMCTPConly,hPrimMC,1,1,"B");
    hPrimMCTPConly->DrawCopy();
    
    for(int i = 1; i <= hPrimMCTPConly->GetNbinsX(); i++) printf("%f\n",hPrimMCTPConly->GetBinContent(i));
    
    
    ofstream myfile (Form("%s/corrTPC_Vs%s%s.txt",dir.Data(),whichVar.Data(),period.Data()));
    if (myfile.is_open())
    {
        for(int i = 1; i <= hPrimMCTPConly->GetNbinsX(); i++) myfile << Form("%f\n",hPrimMCTPConly->GetBinContent(i));
        myfile.close();
    }
    
    ofstream myfile2 (Form("%s/corrTPCErr_Vs%s%s.txt",dir.Data(),whichVar.Data(),period.Data()));
    if (myfile2.is_open())
    {
        for(int i = 1; i <= hPrimMCTPConly->GetNbinsX(); i++) myfile2 << Form("%f\n",hPrimMCTPConly->GetBinError(i));
        myfile2.close();
    }
    
    
    
}
//__________________________________________________________________
void calcFractions(Int_t vsWhichVar, Int_t fPt, TString period, TString dir) {

    vector<double> vectorPrim;
    vector<double> vectorPrimErr;
    vector<double> vectorSec;
    vector<double> vectorSecErr;
    vector<double> vectorPrimInt;

    vector<double> vectorSecPhInt;
    vector<double> vectorSecMatInt;
    vector<double> vectorSecPh;
    vector<double> vectorSecMat;

    vectorPrimInt.clear();
    vectorSecPhInt.clear();
    vectorSecMatInt.clear();
    vectorPrim.clear();
    vectorPrimErr.clear();
    vectorSec.clear();
    vectorSecPh.clear();
    vectorSecMat.clear();
    vectorSecErr.clear();
    
    TString whichPt;
    fUseMCPt   = kFALSE;
    
    if(fPt==1) {
        fUseMCPt = kFALSE;
        whichPt = "globPt";
    }
    else if(fPt==2) {
        fUseMCPt = kTRUE;
        whichPt  = "MCPt";
    }
    TH1F *hPrimFit;
    TH1F *hPrimMC;
    TH1F *hSecFit;
    TH1F *hSecPhFit;
    TH1F *hSecMatFit;
    TH1F *hSecMC;
    TH1F *hSecPhMC;
    TH1F *hSecMatMC;
    TH1F *hRedChi2;
    TH1F *hChi2;
    TH1F *hNDF;
    
    TH1F *hPrimLim;
    TH1F *hSecPhLim;
    TH1F *hSecMatLim;
    
    TString whichVar;
    Int_t nBinTot;
    Int_t npx, npy;
    if(vsWhichVar == 1) {
        hPrimFit = new TH1F("hPrimFit","",npt,ptlims);
        hPrimMC  = new TH1F("hPrimMC","",npt,ptlims);
        hSecFit  = new TH1F("hSecFit","",npt,ptlims);
        hSecPhFit  = new TH1F("hSecPhFit","",npt,ptlims);
        hSecMatFit = new TH1F("hSecMatFit","",npt,ptlims);
        hSecMC     = new TH1F("hSecMC","",npt,ptlims);
        hSecPhMC   = new TH1F("hSecPhMC","",npt,ptlims);
        hSecMatMC  = new TH1F("hSecMatMC","",npt,ptlims);
        hRedChi2 = new TH1F("hRedChi2","",npt,ptlims);
        hChi2    = new TH1F("hChi2","",npt,ptlims);
        hNDF     = new TH1F("hNDF","",npt,ptlims);
        
        hPrimLim = new TH1F("hPrimLim","",npt,ptlims);
        hSecMatLim = new TH1F("hSecMatLim","",npt,ptlims);
        hSecPhLim  = new TH1F("hSecPhLim","",npt,ptlims);

        whichVar = "Pt";
        nBinTot  = npt;
        npx = 3;
        npy = 3;
        printf("**** Analysis versus Pt ****\n");
        
    }
    else if(vsWhichVar == 2) {
        hPrimFit = new TH1F("hPrimFit","",nPhi,0.,6.28);
        hPrimMC  = new TH1F("hPrimMC","",nPhi,0.,6.28);
        hSecFit  = new TH1F("hSecFit","",nPhi,0.,6.28);
        hSecPhFit  = new TH1F("hSecPhFit","",nPhi,0.,6.28);
        hSecMatFit = new TH1F("hSecMatFit","",nPhi,0.,6.28);
        hSecMC     = new TH1F("hSecMC","",nPhi,0.,6.28);
        hSecPhMC   = new TH1F("hSecPhMC","",nPhi,0.,6.28);
        hSecMatMC  = new TH1F("hSecMatMC","",nPhi,0.,6.28);
        hRedChi2 = new TH1F("hRedChi2","",nPhi,0.,6.28);
        hChi2    = new TH1F("hChi2","",nPhi,0.,6.28);
        hNDF     = new TH1F("hNDF","",nPhi,0.,6.28);
        
        nBinTot  = nPhi;
        whichVar = "Phi";
        npx = 5;
        npy = 4;
        printf("**** Analysis versus Phi ****\n");
    }
    else if(vsWhichVar == 3) {
        hPrimFit = new TH1F("hPrimFit","",nEta,-0.8,0.8);
        hPrimMC  = new TH1F("hPrimMC","",nEta,-0.8,0.8);
        hSecFit  = new TH1F("hSecFit","",nEta,-0.8,0.8);
        hSecPhFit  = new TH1F("hSecPhFit","",nEta,-0.8,0.8);
        hSecMatFit = new TH1F("hSecMatFit","",nEta,-0.8,0.8);
        hSecMC     = new TH1F("hSecMC","",nEta,-0.8,0.8);
        hSecPhMC   = new TH1F("hSecPhMC","",nEta,-0.8,0.8);
        hSecMatMC  = new TH1F("hSecMatMC","",nEta,-0.8,0.8);
        hRedChi2 = new TH1F("hRedChi2","",nEta,-0.8,0.8);
        hChi2    = new TH1F("hChi2","",nEta,-0.8,0.8);
        hNDF     = new TH1F("hNDF","",nEta,-0.8,0.8);
        
        nBinTot  = nEta;
        whichVar = "Eta";
        npx = 4;
        npy = 4;
        printf("**** Analysis versus Eta ****\n");
    }
    
    
    TCanvas *cc = new TCanvas("cc","",1300,800);
    cc->Divide(npx,npy);
    TCanvas *c5 = new TCanvas("c5","",1300,800);
    c5->Divide(npx,npy);
    TCanvas *c6 = new TCanvas("c6","",1300,800);
    c6->Divide(npx,npy);
    TCanvas *cPulls = new TCanvas("cPulls","",1300,800);
    cPulls->Divide(npx,npy);
    
    TCanvas *cTFResults = new TCanvas("cTFResults","",1300,800);
    cTFResults->Divide(npx,npy);

    
    TFile *_file0 = TFile::Open(Form("%s/Proj%s_%s_%s_%s.root",dir.Data(),suffix.Data(),whichPt.Data(),period.Data(),whichVar.Data()));
    
    TH1F *hData;
    TH1F *hMC[3];
    double DCAcut = 2.4;
    
    TLegend *legRes = new TLegend(0.6,0.55,0.88,0.88);

    for(int ibin = 0; ibin < nBinTot; ibin++) {
        printf("Bin %d/%d\n",ibin,nBinTot);
        
        hData = (TH1F*)_file0->Get(Form("hData%s%d_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        hMC[0] = (TH1F*)_file0->Get(Form("hMC%s%d_Prim_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        hMC[1] = (TH1F*)_file0->Get(Form("hMC%s%d_Sec_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        hMC[2] = (TH1F*)_file0->Get(Form("hMC%s%d_Mat_%s_%s",whichVar.Data(),ibin,suffix.Data(),whichPt.Data()));
        
        hMC[0]->Sumw2();
        hMC[1]->Sumw2();
        hMC[2]->Sumw2();

        if(hData && hMC[0] && hMC[1] && hMC[2]) {
            TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
            mc->Add(hMC[0]);
            mc->Add(hMC[1]);
            mc->Add(hMC[2]);
            
            ///////////// data-MC comparison ////////////////
            TH1F *hdata  = (TH1F*)hData->Clone(Form("hdata%d",ibin));
            TH1F *hmc    = (TH1F*)hMC[0]->Clone(Form("hmc%d",ibin));
            TH1F *hmcprim    = (TH1F*)hMC[0]->Clone(Form("hmc%dprim",ibin));
            TH1F *hmcsec     = (TH1F*)hMC[1]->Clone(Form("hmc%dsec",ibin));
            TH1F *hmcmat     = (TH1F*)hMC[2]->Clone(Form("hmc%dmat",ibin));
            
            hmc->Add(hMC[1]);
            hmc->Add(hMC[2]);
            
            
            hdata->Sumw2();
            hmc->Sumw2();
            
//            hmc->Scale(hdata->Integral()/hmc->Integral());
//            hmcprim->Scale(hdata->Integral()/hmc->Integral());
//            hmcsec->Scale(hdata->Integral()/hmc->Integral());
//            hmcmat->Scale(hdata->Integral()/hmc->Integral());
            
            hdata->Divide(hdata,hmc,1,1);
            
            c6->cd(ibin+1);
            if(vsWhichVar==1)hdata->SetTitle(Form("data/MC, %.0f < p_{t} < %.0f GeV/c",ptlims[ibin],ptlims[ibin+1]));
            else hdata->SetTitle(Form("data/MC, %s Bin %d",whichVar.Data(),ibin));
            hdata->SetTitleSize(0.06);
            hdata->SetLineColor(kBlue);
            hdata->GetXaxis()->SetLabelSize(0.05);
            hdata->GetXaxis()->SetTitleSize(0.05);
            hdata->GetXaxis()->SetRangeUser(-2,2.);
            hdata->GetYaxis()->SetLabelSize(0.07);
            hdata->DrawCopy();
            TLegend *leg2 = new TLegend(0.1,0.1,0.45,0.2);
            leg2->AddEntry(hdata,"data/MC");
            leg2->Draw("same");
            ///////////////////////////
            cc->cd(ibin+1);
            gPad->SetLogy();
            hData->GetXaxis()->SetLabelSize(0.05);
            hData->GetXaxis()->SetTitleSize(0.05);
            hData->GetXaxis()->SetRangeUser(-2,2.);
            hData->GetYaxis()->SetLabelSize(0.07);
            if(vsWhichVar==1)hData->SetTitle(Form("%.0f < p_{t} < %.0f GeV/c",ptlims[ibin],ptlims[ibin+1]));
            else hData->SetTitle(Form("%s Bin %d",whichVar.Data(),ibin));
            hData->DrawCopy();
            // hmc->Scale(hData->Integral());
            hmc->SetLineColor(kOrange+1);
            hmcprim->SetLineColor(kMagenta+1);
            hmcsec->SetLineColor(kRed+1);
            hmcmat->SetLineColor(kGreen+2);
            hmc->DrawCopy("same");
            hmcprim->DrawCopy("same");
            hmcsec->DrawCopy("same");
            hmcmat->DrawCopy("same");
            TLegend *leg = new TLegend(0.1,0.65,0.45,0.88);
            leg->AddEntry(hData,Form("%s",period.Data()));
            leg->AddEntry(hmc,"MC, all");
            leg->AddEntry(hmcprim,"MC, prim.");
            leg->AddEntry(hmcsec,"MC, sec.");
            leg->AddEntry(hmcmat,"MC, mat.");
            leg->Draw();
            /////////////////////////
            
            TFractionFitter* fit = new TFractionFitter(hData, mc);       // initialise
            if(vsWhichVar==1) {
                fit->Constrain(1,xlowprimVsPt[ibin],xhighprim[ibin]);             // constrain fraction 1 to be between 0 and 1
                fit->Constrain(2,xlowsec[ibin],xhighsecVsPt[ibin]);               // constrain fraction 1 to be between 0 and 1
                fit->Constrain(3,xlowsec[ibin],xhighsecMatVsPt[ibin]);               // constrain fraction 1 to be between 0 and 1
            }
            else if(vsWhichVar==2) {
                fit->Constrain(1,xlowprimVsPhi[ibin],xhighprim[ibin]);             // constrain fraction 1 to be between 0 and 1
                fit->Constrain(2,xlowsec[ibin],xhighsecVsPhi[ibin]);               // constrain fraction 1 to be between 0 and 1
                fit->Constrain(3,xlowsec[ibin],xhighsecMatVsPhi[ibin]);               // constrain fraction 1 to be between 0 and 1
            }
            else if(vsWhichVar==3) {
                fit->Constrain(1,xlowprimVsEta[ibin],xhighprim[ibin]);             // constrain fraction 1 to be between 0 and 1
                fit->Constrain(2,xlowsec[ibin],xhighsecVsEta[ibin]);               // constrain fraction 1 to be between 0 and 1
                fit->Constrain(3,xlowsec[ibin],xhighsecMatVsEta[ibin]);               // constrain fraction 1 to be between 0 and 1
            }
            fit->SetRangeX(hData->FindBin(-fitRange[ibin]),hData->FindBin(fitRange[ibin]));
      
            Int_t status = fit->Fit();                                  // perform the fit
            std::cout << "fit status: " << status << std::endl;
            
            double prim, secMat = 0, sec, secAll = 0,totres,totdata, errMat, errSec;
            if (status == 0) {                                          // check on fit status
                TH1F* result = (TH1F*) fit->GetPlot();
                
                TH1F* PrimMCPred   = (TH1F*)fit->GetMCPrediction(0);
                TH1F* SecMCPred    = (TH1F*)fit->GetMCPrediction(1);
                TH1F* SecMatMCPred = (TH1F*)fit->GetMCPrediction(2);
                PrimMCPred->SetLineColor(kRed);
                SecMCPred->SetLineColor(kGreen+1);
                SecMatMCPred->SetLineColor(kBlue-1);
                
                Double_t value,error;
                
                fit->GetResult(0,prim,error);
                vectorPrim.push_back(prim);
                vectorPrimErr.push_back(error);
                PrimMCPred->Scale(hData->GetSumOfWeights()*prim/PrimMCPred->GetSumOfWeights());

                fit->GetResult(1,sec,errSec);
                vectorSecPh.push_back(sec);
                SecMCPred->Scale(hData->GetSumOfWeights()*sec/SecMCPred->GetSumOfWeights());
         
                fit->GetResult(2,secMat,errMat);
                vectorSecMat.push_back(secMat);
                SecMatMCPred->Scale(hData->GetSumOfWeights()*secMat/SecMatMCPred->GetSumOfWeights());
                
                secAll = sec+secMat;

                vectorSec.push_back(secAll);
                vectorSecErr.push_back(TMath::Sqrt(errSec*errSec+errMat*errMat));
                
                vectorPrimInt.push_back(PrimMCPred->GetSumOfWeights());
                vectorSecPhInt.push_back(SecMCPred->GetSumOfWeights());
                vectorSecMatInt.push_back(SecMatMCPred->GetSumOfWeights());

                
//                Double_t errPrim, errSecAll;
//                prim = PrimMCPred->IntegralAndError(PrimMCPred->FindBin(-DCAcut),PrimMCPred->FindBin(DCAcut),errPrim);
//                secMat = SecMatMCPred->Integral(PrimMCPred->FindBin(-DCAcut),PrimMCPred->FindBin(DCAcut));
//                SecMCPred->Add(SecMatMCPred);
//                secAll  = SecMCPred->IntegralAndError(PrimMCPred->FindBin(-DCAcut),PrimMCPred->FindBin(DCAcut),errSecAll);
                
                cTFResults->cd(ibin+1);
                gPad->SetLogy();
                result->SetLineColor(kOrange);
                PrimMCPred->SetLineColor(kMagenta);
                SecMCPred->SetLineColor(kRed);
                SecMatMCPred->SetLineColor(kGreen+1);
                PrimMCPred->SetTitle(Form("pt %d",ibin));
                PrimMCPred->SetStats(0);
                PrimMCPred->GetXaxis()->SetRangeUser(-2.4,2.4);
                PrimMCPred->DrawCopy("e");
                SecMCPred->DrawCopy("esame");
                SecMatMCPred->DrawCopy("esame");
                hData->DrawCopy("same");
                if(ibin==0) {
                    legRes->AddEntry(hData,"Data");
                    legRes->AddEntry(PrimMCPred,"Prim.");
                    legRes->AddEntry(SecMCPred,"Sec.");
                    legRes->AddEntry(SecMatMCPred,"Mat.");
                    legRes->AddEntry(result,"MC inclusive");
                    legRes->SetTextSize(0.045);
                    
                }
                legRes->Draw("same");
                result->DrawCopy("esame");

                
                hRedChi2->SetBinContent(ibin+1,fit->GetChisquare()/fit->GetNDF());
                hChi2->SetBinContent(ibin+1,fit->GetChisquare());
                hNDF->SetBinContent(ibin+1,fit->GetNDF());
//                vectorPrim.push_back(prim/(prim+secAll));
//                vectorSec.push_back(secAll/(prim+secAll));
//                
//                hPrimFit->SetBinContent(ibin+1,prim/(prim+secAll));
//                hPrimFit->SetBinError(ibin+1,vectorPrimErr[ibin]);
//                hSecFit->SetBinContent(ibin+1,secAll/(prim+secAll));
//                hSecFit->SetBinError(ibin+1,0.00001);
                
                
                hPrimFit->SetBinContent(ibin+1,prim/(prim+sec+secMat));
                hPrimFit->SetBinError(ibin+1,vectorPrimErr[ibin]);
                
                hSecFit->SetBinContent(ibin+1,secAll/(prim+sec+secMat));
                hSecFit->SetBinError(ibin+1,vectorSecErr[ibin]);
                
                hSecPhFit->SetBinContent(ibin+1,sec/(prim+sec+secMat));
                hSecPhFit->SetBinError(ibin+1,errSec);
                
                hSecMatFit->SetBinContent(ibin+1,secMat/(prim+sec+secMat));
                hSecMatFit->SetBinError(ibin+1,errMat);
                
                if(vsWhichVar==1) {
                    hPrimLim->SetBinContent(ibin+1,xlowprimVsPt[ibin]);
                    hSecPhLim->SetBinContent(ibin+1,xhighsecVsPt[ibin]);
                    hSecMatLim->SetBinContent(ibin+1,xhighsecMatVsPt[ibin]);
                }
                if(vsWhichVar==2) {
                    hPrimLim->SetBinContent(ibin+1,xlowprimVsPhi[ibin]);
                    hSecPhLim->SetBinContent(ibin+1,xhighsecVsPhi[ibin]);
                    hSecMatLim->SetBinContent(ibin+1,xhighsecMatVsPhi[ibin]);
                }
                if(vsWhichVar==3) {
                    hPrimLim->SetBinContent(ibin+1,xlowprimVsPhi[ibin]);
                    hSecPhLim->SetBinContent(ibin+1,xhighsecVsPhi[ibin]);
                    hSecMatLim->SetBinContent(ibin+1,xhighsecMatVsPhi[ibin]);
                }
                hPrimLim->SetBinError(ibin+1,0.00000001);
                hSecPhLim->SetBinError(ibin+1,0.00000001);
                hSecMatLim->SetBinError(ibin+1,0.00000001);

                
                c5->cd(ibin+1);
                TH1F *hDCAdata = (TH1F*)hData->Clone(Form("hDCAdata%d",ibin));
                hDCAdata->Sumw2();
                result->Sumw2();
                hDCAdata->Divide(hDCAdata,result,1,1);
                hDCAdata->SetTitle("data / fit");
                hDCAdata->GetXaxis()->SetLabelSize(0.04);
                hDCAdata->GetYaxis()->SetLabelSize(0.07);
                hDCAdata->GetYaxis()->SetLabelSize(0.07);
                if(vsWhichVar==1)hDCAdata->SetTitle(Form("%.0f < p_{t} < %.0f GeV/c",ptlims[ibin],ptlims[ibin+1]));
                else hDCAdata->SetTitle(Form("%s Bin %d",whichVar.Data(),ibin));
                
                if(vsWhichVar==1)hDCAdata->SetTitle(Form("data/fit, %.0f < p_{t} < %.0f GeV/c",ptlims[ibin],ptlims[ibin+1]));
                else hDCAdata->SetTitle(Form("data/fit, %s Bin %d",whichVar.Data(),ibin));
                hDCAdata->SetTitleSize(0.06);
                hDCAdata->GetXaxis()->SetTitleSize(0.05);
                
                hDCAdata->DrawCopy();
                
                cPulls->cd(ibin+1);
                TH1F *hPulls = (TH1F*)hData->Clone(Form("hPullsPt%d",ibin));
                hPulls->Sumw2();
                for(int ii=1; ii<= hPulls->GetNbinsX(); ii++) {
                    double yData = hPulls->GetBinContent(ii);
                    double eData = hPulls->GetBinError(ii);
                    double yRes  = result->GetBinContent(ii);
                    double eRes  = result->GetBinError(ii);
                    if(yData>0. && yRes>0.) {
                        //                        Printf("ydata - yres = %f - %f, errData %f, errres %f",yData,yRes,eData,eRes);
                        hPulls->SetBinContent(ii,(yData-yRes)/TMath::Sqrt(eData*eData+eRes*eRes));
                        hPulls->SetBinError(ii,0.00000001);
                    }
                }
                hPulls->SetMarkerColor(kBlue+2);
                hPulls->SetMarkerStyle(20);
                hPulls->DrawCopy();

            }
            else {
                vectorPrim.push_back(-1.);
                vectorPrimErr.push_back(-1.);
                vectorSec.push_back(-1.);
            }
            
            
            prim   = hMC[0]->Integral();
            sec    = hMC[1]->Integral();
            secMat = hMC[2]->Integral();
            
            hPrimMC->SetBinContent(ibin+1,prim/(prim+sec+secMat));
            hPrimMC->SetBinError(ibin+1,0.0001);
            hSecMC->SetBinContent(ibin+1,(secMat+sec)/(prim+sec+secMat));
            hSecMC->SetBinError(ibin+1,0.0001);
            hSecPhMC->SetBinContent(ibin+1,sec/(prim+sec+secMat));
            hSecPhMC->SetBinError(ibin+1,0.0001);
            hSecMatMC->SetBinContent(ibin+1,secMat/(prim+sec+secMat));
            hSecMatMC->SetBinError(ibin+1,0.0001);

        } else {
            printf("no histos!\n");
            return;
        }
    }
    
    
    cTFResults->SaveAs("Fitresults_TFF.eps");
  
    TCanvas *cFractionsPrim = new TCanvas("cFractionsPrim","Fractions from TFF",700,800);
    gPad->SetLogy();
    hPrimFit->GetYaxis()->SetRangeUser(0.7,1.);

    hPrimFit->SetLineColor(kBlue);
    hSecPhFit->SetLineColor(kRed);
    hSecMatFit->SetLineColor(kMagenta);
    hPrimFit->SetMarkerColor(kBlue);
    hSecPhFit->SetMarkerColor(kRed);
    hSecMatFit->SetMarkerColor(kMagenta);
    hPrimFit->SetMarkerStyle(20);
    hSecPhFit->SetMarkerStyle(20);
    hSecMatFit->SetMarkerStyle(20);
    
    
    hPrimMC->SetLineColor(kBlue);
    hSecPhMC->SetLineColor(kRed);
    hSecMatMC->SetLineColor(kMagenta);
    hPrimMC->SetMarkerColor(kBlue);
    hSecPhMC->SetMarkerColor(kRed);
    hSecMatMC->SetMarkerColor(kMagenta);
    hPrimMC->SetMarkerStyle(3);
    hSecPhMC->SetMarkerStyle(3);
    hSecMatMC->SetMarkerStyle(3);
    
    
    hPrimLim->SetLineColor(kBlue);
    hSecPhLim->SetLineColor(kRed);
    hSecMatLim->SetLineColor(kMagenta);
    
    hPrimMC->SetLineStyle(kDashed);
    hSecMC->SetLineStyle(kDashed);
    hSecPhMC->SetLineStyle(kDashed);
    hSecMatMC->SetLineStyle(kDashed);
    hPrimLim->SetLineStyle(kDashed);
    hSecPhLim->SetLineStyle(kDashed);
    hSecMatLim->SetLineStyle(kDashed);
    hPrimLim->SetMarkerColor(kBlue);
    hSecPhLim->SetMarkerColor(kRed);
    hSecMatLim->SetMarkerColor(kMagenta);
    
    hPrimLim->SetMarkerStyle(4);
    hSecPhLim->SetMarkerStyle(4);
    hSecMatLim->SetMarkerStyle(4);
    
    hPrimFit->SetTitle("Fraction resulting from TFF");
    hPrimFit->SetStats(0);
    hPrimFit->DrawCopy();
    hPrimLim->DrawCopy("psame");
    hPrimMC->DrawCopy("psame");
    
    //    cFractions->SaveAs("ClosureTestMC.eps");
    cFractionsPrim->SaveAs("FractionsPrimFromFit.eps");
    
    
    TH1F *hcpPrim = (TH1F*)hPrimFit->Clone("hPrimCl");
    TH1F *hcpPrimMC = (TH1F*)hPrimMC->Clone("hPrimClMC");
    TLegend *legFrac = new TLegend(0.6,0.3,0.88,0.55);
    legFrac->AddEntry(hcpPrim,"prim.","pl");
    legFrac->AddEntry(hcpPrimMC,"MC prim.","pl");
    legFrac->AddEntry(hPrimLim,"Input fit prim.","pl");
    
    legFrac->Draw();
    
    TCanvas *cFractionsSec = new TCanvas("cFractionsSec","Fractions from TFF",700,800);
    gPad->SetLogy();
    hSecPhFit->GetYaxis()->SetRangeUser(0.00005,0.4);
    
    
    hSecPhFit->SetTitle("Fraction resulting from TFF");
    hSecPhFit->SetStats(0);
    hSecPhFit->DrawCopy();
    hSecMatFit->DrawCopy("same");
    hSecPhLim->DrawCopy("psame");
    hSecMatLim->DrawCopy("psame");
    
    hSecPhMC->DrawCopy("psame");
    hSecMatMC->DrawCopy("psame");
    
    
    TLegend *legFrac2 = new TLegend(0.6,0.12,0.88,0.4);
    legFrac2->AddEntry(hSecPhFit,"sec.","pl");
    legFrac2->AddEntry(hSecMatFit,"sec. mat.","pl");
    legFrac2->AddEntry(hSecPhMC,"MC sec.","pl");
    legFrac2->AddEntry(hSecMatMC,"MC sec. mat.","pl");
    legFrac2->AddEntry(hSecPhLim,"Input fit sec.","pl");
    legFrac2->AddEntry(hSecMatLim,"Input fit sec. mat.","pl");
    
    legFrac2->Draw();
    
    cFractionsSec->SaveAs("FractionsSecFromFit.eps");
    
    
    
    hPrimFit->SetLineColor(kRed);
    hPrimFit->SetMarkerColor(kRed);
    hPrimMC->SetLineColor(kRed);
    hPrimMC->SetMarkerColor(kRed);
    hPrimFit->SetMarkerStyle(20);
    hPrimFit->SetMarkerSize(0.4);
    hPrimMC->SetMarkerStyle(4);
    hPrimMC->SetMarkerSize(0.4);
    
    hSecFit->SetLineColor(kBlue);
    hSecFit->SetMarkerColor(kBlue);
    hSecMC->SetLineColor(kBlue);
    hSecMC->SetMarkerColor(kBlue);
    hSecFit->SetMarkerStyle(20);
    hSecFit->SetMarkerSize(0.4);
    hSecMC->SetMarkerStyle(4);
    hSecMC->SetMarkerSize(0.4);
    hPrimFit->GetYaxis()->SetRangeUser(0.,1.);
    hPrimFit->GetYaxis()->SetTitle("Fractions");
    if(vsWhichVar==1) hPrimFit->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    else if(vsWhichVar==2) hPrimFit->GetXaxis()->SetTitle("phi");
    else if(vsWhichVar==3) hPrimFit->GetXaxis()->SetTitle("eta");
    hPrimFit->GetYaxis()->SetRangeUser(-0.1,1.6);
    hPrimFit->SetStats(0);
    
    TLegend *l = new TLegend(0.15,0.7,0.45,0.88);
    l->AddEntry(hPrimFit,"TFracFitt primaries");
    l->AddEntry(hPrimMC,"MC primaries");
    l->AddEntry(hSecFit,"TFracFitt sec.");
    l->AddEntry(hSecMC,"MC sec.");
    
    TCanvas *c1 = new TCanvas("c1","",1200,800);
    c1->Divide(2,2);
    c1->cd(1);
    hPrimFit->DrawCopy("EP");
    hPrimMC->DrawCopy("EPsame");
    hSecFit->DrawCopy("EPsame");
    hSecMC->DrawCopy("EPsame");
    l->Draw("same");
    c1->cd(2);
    hRedChi2->SetStats(0);
    hRedChi2->SetTitle("#chi^{2}/ndf");
    hRedChi2->DrawCopy("EP");
    c1->cd(3);
    hPrimMC->Sumw2();
    hPrimFit->Sumw2();
    hSecFit->Sumw2();
    hSecMC->Sumw2();
    hPrimFit->Divide(hPrimMC);
    hSecFit->Divide(hSecMC);
    hPrimFit->GetYaxis()->SetTitle("fit/MC fractions");
    hPrimFit->DrawCopy();
    hSecFit->DrawCopy("same");
    c1->cd(4);
    hNDF->SetMarkerStyle(20);
    hNDF->GetYaxis()->SetRangeUser(hNDF->GetBinContent(1)-50,hNDF->GetBinContent(1)+50);
    hNDF->SetStats(0);
    hNDF->SetTitle("ndf");
    hNDF->DrawCopy("EP");
    
    c1->SaveAs(Form("%s/fractions_%s_ESDTrOnly_%s_Vs%s.eps",dir.Data(),suffix.Data(),period.Data(),whichVar.Data()));
    
    
    Printf("*********************");
    Printf("Primary fractions:");
    for(int i=0; i<(int)vectorPrim.size();i++) {
        Printf("%f +- %f, --> re-normalized: %f\n",vectorPrim[i],vectorPrimErr[i],vectorPrim[i]/(vectorPrim[i]+vectorSecPh[i]+vectorSecMat[i]));
    }
    Printf("*********************");
    Printf("Secondary fractions:");
    for(int i=0; i<(int)vectorSecPh.size();i++) {
        Printf("%f --> re-normalized: %f\n",vectorSecPh[i],vectorSecPh[i]/(vectorPrim[i]+vectorSecPh[i]+vectorSecMat[i]));
    }
    Printf("*********************");
    Printf("Material fractions:");
    for(int i=0; i<(int)vectorSecMat.size();i++) {
        Printf("%f --> re-normalized: %f\n",vectorSecMat[i],vectorSecMat[i]/(vectorPrim[i]+vectorSecPh[i]+vectorSecMat[i]));
    }
    ofstream myfile (Form("%s/fractions_%s_%s_Vs%s.txt",dir.Data(),suffix.Data(),period.Data(),whichVar.Data()));
    if (myfile.is_open())
    {
        for(int i = 0; i <(int)vectorPrim.size(); i++) myfile << Form("%f\n",vectorPrim[i]/(vectorPrim[i]+vectorSecPh[i]+vectorSecMat[i]));
        myfile.close();
    }
    
    ofstream myfile2 (Form("%s/fractionErrs_%s_%s_Vs%s.txt",dir.Data(),suffix.Data(),period.Data(),whichVar.Data()));
    if (myfile2.is_open())
    {
       for(int i = 0; i <(int)vectorPrimErr.size(); i++) myfile2 << Form("%f\n",vectorPrimErr[i]);
        myfile2.close();
    }
    
    ofstream myfile3 (Form("%s/Secfractions_%s_%s_Vs%s.txt",dir.Data(),suffix.Data(),period.Data(),whichVar.Data()));
    if (myfile3.is_open()) {
        for(int i = 0; i <(int)vectorSec.size(); i++){
            myfile3 << Form("%f\n",vectorSec[i]/(vectorPrim[i]+vectorSecPh[i]+vectorSecMat[i]));
        }
        myfile3.close();
    }
    
    ofstream myfile4 (Form("%s/SecfractionErrs_%s_%s_Vs%s.txt",dir.Data(),suffix.Data(),period.Data(),whichVar.Data()));
    if (myfile4.is_open()) {
        for(int i = 0; i <(int)vectorSecErr.size(); i++){
            myfile4 << Form("%f\n",vectorSecErr[i]);
        }
        myfile4.close();
    }

    
    delete hData;
    for(int i=0; i<3; i++) delete hMC[i];
}


//__________________________________________________________________
TH1D* projectMC(THnSparse *hSparse, Double_t partType, Int_t ibin) {
    
    Printf("\n*** MC, bin %d ***",ibin);
    
    Int_t nbins;
    TAxis* ax[nVarMC];
    for (Int_t ii=0; ii<nVarMC; ii++) {
        ax[ii] = (TAxis*)hSparse->GetAxis(ii);
        nbins = ax[ii]->GetNbins();
        ax[ii]->SetRange(0,nbins+1);
    }
    
    
    Int_t first=0;
    Int_t last=0;
    Double_t binWidth;
    
    Int_t iax;
    
    //cut on pt
    if(fUseMCPt)  iax = 3;
    else iax = 2;
    binWidth = ax[iax]->GetBinWidth(1);
    if(fWhichVar==1) {  //analysis vs pt
        first = ax[iax]->FindBin(ptlims[ibin]+binWidth/2.);
        last  = ax[iax]->FindBin(ptlims[ibin+1]-binWidth/2.);
    }
    else {
        first = ax[iax]->FindBin(ptlims[0]+binWidth/2.);
        last  = ax[iax]->FindBin(ptlims[npt]-binWidth/2.);
    }
    ax[iax]->SetRange(first,last);
    Printf("Pt : %g, %g",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    
    //cut on phi
    iax = 4;
    binWidth = ax[iax]->GetBinWidth(1);
    if(fWhichVar==2) {  //analysis vs phi
        first = ibin+1;
        last  = ibin+1;
    }
    else {
        first = 1;
        last  = nPhi;
    }
    ax[iax]->SetRange(first,last);
    Printf("Phi : %g, %g",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    
    //cut on eta
    iax = 5;
    binWidth = ax[iax]->GetBinWidth(1);
    if(fWhichVar==3) {  //analysis vs eta
        first = ibin+1;
        last  = ibin+1;
    }
    else {
        first = ax[iax]->FindBin(-etaRange);
        last  = ax[iax]->FindBin(etaRange-binWidth/2.);
    }
    ax[iax]->SetRange(first,last);
    Printf("Eta : %g, %g",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    
    
    //prim, sec, mat
    iax = 6;
    first = ax[iax]->FindBin(partType);
    last  = ax[iax]->FindBin(partType);
    ax[iax]->SetRange(first,last);
    Printf("Part.type : %g, %g",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    

    //TOF bc
    if(fUseTOFbc) {
        iax = 9;
        binWidth = ax[iax]->GetBinWidth(1);
        first = ax[iax]->FindBin(-0.5 + binWidth/2.);
        last  = ax[iax]->FindBin(0.5  - binWidth/2.);
        ax[iax]->SetRange(first,last);
        Printf("TOF bc : %g, %g",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    }

    TH1D* hproj = (TH1D*)hSparse->Projection(0);
    return hproj;
}

//__________________________________________________________________
TH1D* projectData(THnSparse *hSparse, Int_t ibin) {
    
    Printf("\n *** data, bin %d ***",ibin);
    Int_t nbins;
    TAxis* ax[nVarData];
    for (Int_t ii=0; ii<nVarData; ii++) {
        ax[ii] = (TAxis*)hSparse->GetAxis(ii);
        nbins = ax[ii]->GetNbins();
        ax[ii]->SetRange(0,nbins+1);
    }
    
    
    Int_t first=0;
    Int_t last=0;
    Double_t binWidth;
    
    Int_t iax;

    //cut on pt
    iax = 2;
    binWidth = ax[iax]->GetBinWidth(1);
    if(fWhichVar==1) {  //analysis vs pt
        first = ax[iax]->FindBin(ptlims[ibin]+binWidth/2.);
        last  = ax[iax]->FindBin(ptlims[ibin+1]-binWidth/2.);
    }
    else {
        first = ax[iax]->FindBin(ptlims[0]+binWidth/2.);
        last  = ax[iax]->FindBin(ptlims[npt]-binWidth/2.);
    }
    ax[iax]->SetRange(first,last);
    printf("Pt : %g, %g \n",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    
    //cut on phi
    iax = 3;
    binWidth = ax[iax]->GetBinWidth(1);
    if(fWhichVar==2) {  //analysis vs phi
        first = ibin+1;
        last  = ibin+1;
    }
    else {
        first = 1;
        last  = nPhi;
    }
    ax[iax]->SetRange(first,last);
    printf("Phi : %g, %g \n",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    
    //cut on eta
    iax = 4;
    binWidth = ax[iax]->GetBinWidth(1);
    if(fWhichVar==3) {  //analysis vs eta
        first = ibin+1;
        last  = ibin+1;
    }
    else {
        first = ax[iax]->FindBin(-0.8);
        last  = ax[iax]->FindBin(0.8-binWidth/2.);
    }
    ax[iax]->SetRange(first,last);
    printf("Eta : %g, %g \n",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    
    //request on TOF bc==0
    if(fUseTOFbc) {
        iax = 5;
        binWidth = ax[iax]->GetBinWidth(1);
        first = ax[iax]->FindBin(-0.5 + binWidth/2.);
        last  = ax[iax]->FindBin(0.5  - binWidth/2.);
        ax[iax]->SetRange(first,last);
        Printf("TOF bc : %g, %g",ax[iax]->GetBinLowEdge(first),ax[iax]->GetBinUpEdge(last));
    }

    
    TH1D* hproj = (TH1D*)hSparse->Projection(0);
    return hproj;
    
}

