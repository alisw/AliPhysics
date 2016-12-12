#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THn.h"
#include "TList.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include <Riostream.h>
#include <fstream>
#include <string>
#endif

using std::cout;
using std::endl;

enum particles {kPion,kPionKaon,kKaon,kElectron,kProton,kAll}; //the user can add his own preferred combination
enum var4analysis {kPt, kPhi, kEta}; //preferred variable in the analysis (vs pt, vs phi or vs eta)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Here set file names, variable, particle species and types and eventual suffix      ////////////////
///////////// To launch the script: .x CalculateMathcingEfficiency.C+                            ////////////////
///////////// The script produces output txt files with matching efficiency for data, inclusive  ////////////////
///////////// MC, MC primary and secondary particles.                                            ////////////////
////////////////////////////////////////////////////////////////////////////////x/////////////////////////////////

TString period = "_lowIR"; // suffix to append in the ouput txt files
void CalculateMathcingEfficiency(TString inFileNameData = "15o/AnalysisResults_Data_lowIR.root",//or mc inclusive or mc primaries
                                 TString inFileNameMC = "15o/AnalysisResults_MC_lowIR.root", // mc secondaries
                                 Int_t particleType = kPionKaon,
                                 Int_t variable = kPt);

TH1D *GetITSTPCMatchingHisto(TString inFileNameData = "15o/AnalysisResults_Data_lowIR.root",
                             Int_t particleType = kPionKaon,
                             Int_t species = -1,
                             Bool_t isMatched = kTRUE,
                             Float_t pt=10.,
                             Float_t pt1=11.,
                             Int_t variable = kPt,
                             TString listname="");

//_______________________________________________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CalculateMathcingEfficiency(TString inFileNameData,
                                 TString inFileNameMC,
                                 Int_t particleType,
                                 Int_t variable){
    TString suffix;
    if(particleType==kPion)          suffix = "_Pi";
    else if(particleType==kPionKaon) suffix = "_PiK";
    else if(particleType==kKaon)     suffix = "_K";
    else if(particleType==kElectron) suffix = "_e";
    else if(particleType==kProton)   suffix = "_p";
    else if(particleType==kAll)      suffix = "_All";
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    //
    // make a single plot
    //
    
    const int nptbins = 16;
    TCanvas *canvEff = new TCanvas("canvEff","matching efficiency",1000,800);
    //    canvEff->Divide(2,2);
    canvEff->Divide(4,4);
    TCanvas *c[nptbins];
    
    Float_t coeffData[nptbins];
    Float_t coeffMCall[nptbins];
    Float_t coeffPrim[nptbins];
    Float_t coeffSec[nptbins];
    Float_t errData[nptbins];
    Float_t errMC[nptbins];
    Float_t errPrim[nptbins];
    Float_t errSec[nptbins];
    
    Float_t Pt[nptbins]    = {0.1,1.2,1.95,2.9,4.,5,6.2,7.,8.,8.6,9.6,10.6,11.8,13.09,14.55,16.2};
    Float_t array[nptbins] = {0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};
    
    int species[4] = {-1,99,1,0};//data, mcall(prim+sec), prim, sec
    int colors[4]  = {1,2,4,3};
    
    double err[nptbins][4];
    double coeff[nptbins][4];
    
    TH1D *hITSTPCtracks[4];
    TH1D *hTPCtracks[4];
    
    TF1 *fitFunc[4];
    
    TString names[4] = {"Data","MCall","MCprim","MCsec"};
    
    for(Int_t pt = 0; pt < nptbins; pt++) {
        
        TString listname  = Form("trackingUncert%s",suffix.Data());//data directory name  or MC primaries or MC inclusive
        
        for(int is = 0; is < 4; is++) {
            //ITSTPC
            if(is==0) hITSTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameData, particleType, species[is], kTRUE, Pt[pt], Pt[pt+1], variable,listname);
            else hITSTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameMC, particleType, species[is], kTRUE, Pt[pt], Pt[pt+1], variable,listname);
            if(!hITSTPCtracks[is])return;
            hITSTPCtracks[is]->SetNameTitle(Form("hITSTPCtracks%sPt%d",names[is].Data(),pt), Form("ITS-TPC, %s: pt %d",names[is].Data(),pt));
            hITSTPCtracks[is]->Sumw2();
            
            //TPC only
            if(is==0) hTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameData, particleType,species[is], kFALSE, Pt[pt],Pt[pt+1],variable,listname);
            else hTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameMC, particleType, species[is], kFALSE, Pt[pt], Pt[pt+1], variable,listname);
            if(!hTPCtracks[is])return;
            hTPCtracks[is]->SetNameTitle(Form("hTPCtracks%sPt%d",names[is].Data(),pt), Form("TPC, %s: pt %d",names[is].Data(),pt));
            hTPCtracks[is]->Sumw2();
        }
        c[pt] = new TCanvas(Form("c%d",pt));
        c[pt]->cd(1);
        for(int is = 0; is < 4; is++) {
            hITSTPCtracks[is]->Divide(hITSTPCtracks[is],hTPCtracks[is],1,1,"B");
            hITSTPCtracks[is]->GetYaxis()->SetRangeUser(0.,1.);
            hITSTPCtracks[is]->SetLineColor(colors[is]);
            hITSTPCtracks[is]->SetMarkerColor(colors[is]);
            
            err[pt][is] = hITSTPCtracks[is]->GetBinError(1);
            if(is == 0)hITSTPCtracks[is]->Draw();
            else hITSTPCtracks[is]->Draw("same");
            
            if(fitFunc[is]) {
                fitFunc[is] = 0x0;
                delete fitFunc[is];
            }
            fitFunc[is] = new TF1(Form("fitFunc%d",is), "pol0",array[pt], array[pt+1]);
            hITSTPCtracks[is]->Fit(fitFunc[is], "0 R W");
            coeff[pt][is]  = TMath::Abs(fitFunc[is]->GetParameter(0));
        }
        
        canvEff->cd(pt+1);
        hITSTPCtracks[2]->DrawCopy();
        hITSTPCtracks[3]->DrawCopy("sames");
        
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
        
        
        ofstream f1(Form("MatchEff_data%s.txt",period.Data()));
        ofstream f2(Form("MatchEff_data_err%s.txt",period.Data()));
        if (f1.is_open())
            for(Int_t i=0;i<nptbins;i++){
                f1 << Form("%f\n",coeff[i][0]);
                f2 << Form("%f\n",err[i][0]);
            }
        f1.close();
        f2.close();
        
        ofstream f3(Form("MatchEff_MCall%s.txt",period.Data()));
        ofstream f4(Form("MatchEff_MCall_err%s.txt",period.Data()));
        if (f3.is_open())
            for(Int_t i=0;i<nptbins;i++){
                f3 << Form("%f\n",coeff[i][1]);
                f4 << Form("%f\n",err[i][1]);
            }
        f3.close();
        f4.close();
        
        ofstream f5(Form("MatchEff_MCprim%s.txt",period.Data()));
        ofstream f6(Form("MatchEff_MCprim_err%s.txt",period.Data()));
        if (f5.is_open())
            for(Int_t i=0;i<nptbins;i++){
                f5 << Form("%f\n",coeff[i][2]);
                f6 << Form("%f\n",err[i][2]);
            }
        f5.close();
        f6.close();
        
        ofstream f7(Form("MatchEff_MCsec%s.txt",period.Data()));
        ofstream f8(Form("MatchEff_MCsec_err%s.txt",period.Data()));
        if (f7.is_open())
            for(Int_t i=0;i<nptbins;i++){
                f7 << Form("%f\n",coeff[i][3]);
                f8 << Form("%f\n",err[i][3]);
            }
        f7.close();
        f8.close();
    }
    return;
    
    
}

//______________________________________________________________________________
TH1D *GetITSTPCMatchingHisto(TString inFileNameData,
                             Int_t particleType,
                             Int_t species,
                             Bool_t isMatched,
                             Float_t pt,
                             Float_t pt1,
                             Int_t variable,
                             TString listname){
    
    //
    // ITS-TPC matching histograms as a funct. of pT, eta, phi, for each species
    //
    //axis 0 isMatched: -0.5 0.5, 1.5 (central value 0, 1)
    //axis 1 pT 0 - 20  0.1 --> 20 (50 bins)
    //axis 2 eta -1, 1 20 bins
    //axis 3 phi 0, 6.28 18 0, 2*TMath::Pi()
    //axis 4 pid 0 electrons, 1 pions,2 kaons 3 protont   4 all  -0.5, 5.5 6 bins
    //axis 5 : -1 all, 0 secondaries, 1 primaries, 3 bins -0.5, 2.5
    
    TFile * inFileData = TFile::Open(inFileNameData);
    TList * l = (TList * ) inFileData->Get(listname);
    
    //trackingUncert
    THnF * histITSTPC = (THnF *) l->FindObject("histTpcItsMatch");
    
    // select particleType
    double xmin,xmax;
    if(particleType==kElectron) {xmin=-0.5; xmax=0.5;}
    else if(particleType==kPion)     {xmin=0.5;  xmax=1.5;}
    else if(particleType==kKaon)     {xmin=1.5;  xmax=2.5;}
    else if(particleType==kProton)   {xmin=2.5;  xmax=3.5;}
    else if(particleType==kAll)      {xmin=3.5;  xmax=4.5;}
    else if(particleType==kPionKaon) {xmin=0.5;  xmax=2.5;}
    else {Printf("particleType not yet impelemented! Please provide it first!"); return 0x0;}
    
    histITSTPC->GetAxis(4)->SetRangeUser(xmin,xmax);
    histITSTPC->GetAxis(5)->SetRangeUser(species-0.5,species+0.5);
 
    
    //ptbins
    histITSTPC->GetAxis(1)->SetRangeUser(pt,pt1);
    if(isMatched)
        histITSTPC->GetAxis(0)->SetRangeUser(1,1);
    else
        histITSTPC->GetAxis(0)->SetRangeUser(0,0);
    
    int projectionAxis=1;
    if(variable==kPt) projectionAxis = 1;
    else if(variable==kPhi) projectionAxis = 2;
    else if(variable==kEta) projectionAxis = 3;

    TH1D * hProj  = histITSTPC->Projection(projectionAxis);
    hProj->SetDirectory(0);
    hProj->SetNameTitle(Form("h%d",projectionAxis),Form("h%d",projectionAxis));
    //    hProj->Draw();
    //
    delete l;
    inFileData->Close();
    return hProj;
}
