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
///////////// To launch the script: .x CalculateMatchingEfficiency.C+                            ////////////////
///////////// The script produces output txt files with matching efficiency for data, inclusive  ////////////////
///////////// MC, MC primary and secondary particles.                                            ////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CalculateMatchingEfficiency(TString inFileNameData = "15o/AnalysisResults_Data_pidfix_Golden.root",//or mc inclusive or mc primaries
                                 TString inFileNameMC = "15o/AnalysisResults_MC_pidfix_Golden.root", // mc secondaries
                                 TString dir = "15o",
                                 TString period = "_highIR_pidfix",         // suffix to append in the ouput txt files
                                 Int_t particleType = kPionKaon,
                                 Int_t variable = kPhi);

TH1D *GetITSTPCMatchingHisto(TString inFileNameData = "15o/AnalysisResults_Data_megarun.root",
                             Int_t particleType = kPionKaon,
                             Int_t species = -1,
                             Bool_t isMatched = kTRUE,
                             Float_t pt=10.,
                             Float_t pt1=11.,
                             Int_t variable = kPt,
                             TString listname="");


//Set here the variables
const int nptbins = 16;  // nbins for vs pt analysis
const int netabins = 16; // nbins for vs eta analysis
const int nphibins = 18; // nbins for vs phi analysis

Float_t minEta = -0.8;
Float_t maxEta =  0.8;
Float_t minPhi = 0.;
Float_t maxPhi = 2*TMath::Pi();
Float_t DCAxy  = 2.4;
Float_t DCAz   = 3.2;
Bool_t fIsSelectTOFBC=kFALSE;

//Don't touch these. they come from log. axis in the sparse
Float_t Pt[nptbins]    = {0.1,1.2,1.95,2.9,4.,5,6.2,7.,8.,8.6,9.6,10.6,11.8,13.09,14.55,16.2};
Float_t Ptarray[nptbins] = {0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};

//_______________________________________________________________________________________________________________
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CalculateMatchingEfficiency(TString inFileNameData,
                                 TString inFileNameMC,
                                 TString dir,
                                 TString period,
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
    Int_t nBinTot;
    Double_t lowLim, upLim;
    TString suffVar = "";
    if(variable == kPt)       {nBinTot = nptbins; suffVar = "vsPt";}
    else if(variable == kPhi) {nBinTot = nphibins; suffVar = "vsPhi";}
    else if(variable == kEta) {nBinTot = netabins; suffVar = "vsEta";}
    else {Printf("Please define the variable under study!"); return;}
    
    TCanvas *canvEff = new TCanvas("canvEff","matching efficiency",1000,800);
    canvEff->Divide(4,4);
    TCanvas *c[nBinTot];
    
    Float_t coeffData[nBinTot];
    Float_t coeffMCall[nBinTot];
    Float_t coeffPrim[nBinTot];
    Float_t coeffSec[nBinTot];
    Float_t errData[nBinTot];
    Float_t errMC[nBinTot];
    Float_t errPrim[nBinTot];
    Float_t errSec[nBinTot];
    
    
    int species[4] = {-1,2,1,0};//data, mcall(prim+sec), prim, sec
    int colors[4]  = {1,2,4,3};
    
    double err[nBinTot][4];
    double coeff[nBinTot][4];
    
    TH1D *hITSTPCtracks[4];
    TH1D *hTPCtracks[4];
    
    TF1 *fitFunc[4];
    
    TString names[4] = {"Data","MCall","MCprim","MCsec"};
    
    for(Int_t ibin = 0; ibin < nBinTot; ibin++) {
        
        TString listname  = Form("trackingUncert%s",suffix.Data());//data directory name  or MC primaries or MC inclusive
        if(variable == kPt)       {lowLim = Pt[ibin]; upLim = Pt[ibin+1];}
        else if(variable == kPhi) {lowLim = minPhi + ibin*(maxPhi-minPhi)/nphibins; upLim = minPhi + (ibin+1)*(maxPhi-minPhi)/nphibins;}
        else if(variable == kEta) {lowLim = minEta + ibin*(maxEta-minEta)/netabins; upLim = minEta + (ibin+1)*(maxEta-minEta)/netabins;}
        
        for(int is = 0; is < 4; is++) {
            //ITSTPC
            if(is==0) hITSTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameData, particleType, species[is], kTRUE, lowLim, upLim, variable,listname);
            else hITSTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameMC, particleType, species[is], kTRUE, lowLim, upLim, variable,listname);
            if(!hITSTPCtracks[is])return;
            hITSTPCtracks[is]->SetNameTitle(Form("hITSTPCtracks%s%s%d",names[is].Data(),suffVar.Data(),ibin), Form("ITS-TPC, %s: %s %d",names[is].Data(),suffVar.Data(),ibin));
            hITSTPCtracks[is]->Sumw2();
            
            //TPC only
            if(is==0) hTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameData, particleType,species[is], kFALSE, lowLim, upLim,variable,listname);
            else hTPCtracks[is] = GetITSTPCMatchingHisto(inFileNameMC, particleType, species[is], kFALSE, lowLim, upLim, variable,listname);
            if(!hTPCtracks[is])return;
            hTPCtracks[is]->SetNameTitle(Form("hTPCtracks%s%s%d",names[is].Data(),suffVar.Data(),ibin), Form("TPC, %s: %s %d",names[is].Data(),suffVar.Data(),ibin));
            hTPCtracks[is]->Sumw2();
        }
        c[ibin] = new TCanvas(Form("c%d",ibin));
        c[ibin]->cd(1);
        for(int is = 0; is < 4; is++) {
            hITSTPCtracks[is]->Divide(hITSTPCtracks[is],hTPCtracks[is],1,1,"B");
            hITSTPCtracks[is]->GetYaxis()->SetRangeUser(0.,1.);
            hITSTPCtracks[is]->SetLineColor(colors[is]);
            hITSTPCtracks[is]->SetMarkerColor(colors[is]);
            
            //            if(is == 0)hITSTPCtracks[is]->Draw();
            //            else hITSTPCtracks[is]->Draw("same");
            //
            //            to be adapted if you have log axis
            //            if(variable == kPt) {
            //                if(fitFunc[is]) {
            //                    fitFunc[is] = 0x0;
            //                    delete fitFunc[is];
            //                }
            //                fitFunc[is] = new TF1(Form("fitFunc%d",is), "pol0",Ptarray[ibin], Ptarray[ibin+1]);
            //                hITSTPCtracks[is]->Fit(fitFunc[is], "0 R W");
            //                coeff[ibin][is]  = TMath::Abs(fitFunc[is]->GetParameter(0));
            //            } else {coeff[ibin][is]  = hITSTPCtracks[is]->GetBinContent(1);Printf("ME: %f",coeff[ibin][is]);}
            //
            //            err[ibin][is] = hITSTPCtracks[is]->GetBinError(ibin+1);
            //        }
            //
            //        canvEff->cd(ibin+1);
            //        hITSTPCtracks[2]->DrawCopy();
            //        hITSTPCtracks[3]->DrawCopy("sames");
            
            double averME  = 0;
            double eaverME = 0;
            for(int ib = 1; ib <= hITSTPCtracks[is]->GetNbinsX(); ib++) {
                averME  += hITSTPCtracks[is]->GetBinContent(ib);
                eaverME += hITSTPCtracks[is]->GetBinError(ib);
            }
            averME /= hITSTPCtracks[is]->GetNbinsX();
            eaverME /= hITSTPCtracks[is]->GetNbinsX();
            
            coeff[ibin][is]  = averME;
            err[ibin][is] = eaverME;
        }

        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
        
    }
    ofstream f1(Form("%s/MatchEff_data%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    ofstream f2(Form("%s/MatchEff_data_err%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    if (f1.is_open())
        for(Int_t i=0;i<nBinTot;i++){
            f1 << Form("%f\n",coeff[i][0]);
            f2 << Form("%f\n",err[i][0]);
        }
    f1.close();
    f2.close();
    
    ofstream f3(Form("%s/MatchEff_MCall%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    ofstream f4(Form("%s/MatchEff_MCall_err%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    if (f3.is_open())
        for(Int_t i=0;i<nBinTot;i++){
            f3 << Form("%f\n",coeff[i][1]);
            f4 << Form("%f\n",err[i][1]);
        }
    f3.close();
    f4.close();
    
    ofstream f5(Form("%s/MatchEff_MCprim%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    ofstream f6(Form("%s/MatchEff_MCprim_err%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    if (f5.is_open())
        for(Int_t i=0;i<nBinTot;i++){
            f5 << Form("%f\n",coeff[i][2]);
            f6 << Form("%f\n",err[i][2]);
        }
    f5.close();
    f6.close();
    
    ofstream f7(Form("%s/MatchEff_MCsec%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    ofstream f8(Form("%s/MatchEff_MCsec_err%s_%s.txt",dir.Data(),period.Data(),suffVar.Data()));
    if (f7.is_open())
        for(Int_t i=0;i<nBinTot;i++){
            f7 << Form("%f\n",coeff[i][3]);
            f8 << Form("%f\n",err[i][3]);
        }
    f7.close();
    f8.close();
    
    return;
    
    
}

//______________________________________________________________________________
TH1D *GetITSTPCMatchingHisto(TString inFileNameData,
                             Int_t particleType,
                             Int_t species,
                             Bool_t isMatched,
                             Float_t lowLim,
                             Float_t upLim,
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
    //axis 6 : 0: select on TOFBC, 1: no TOFbc
    //axis 7 : DCAxy
    //axis 8 : DCAz
    
    TFile * inFileData = TFile::Open(inFileNameData);
    TList * l = (TList * ) inFileData->Get(listname);
    
    //trackingUncert
    THnF * histITSTPC = (THnF *) l->FindObject("histTpcItsMatch");
    
    //TPC-ITS or TPC only
    if(isMatched) histITSTPC->GetAxis(0)->SetRangeUser(1,1);
    else histITSTPC->GetAxis(0)->SetRangeUser(0,0);
    
    double binw;
    //ptbins
    binw    = histITSTPC->GetAxis(1)->GetBinWidth(1)/2.;
    if(variable == kPt) {
        lowLim += binw;
        upLim  -= binw;
        histITSTPC->GetAxis(1)->SetRangeUser(lowLim,upLim);
    }
    else histITSTPC->GetAxis(1)->SetRangeUser(Ptarray[0]+binw,Ptarray[nptbins-1]-binw);
    
    //etabins
    binw    = histITSTPC->GetAxis(2)->GetBinWidth(1)/2.;
    if(variable == kEta) {
        lowLim += binw;
        upLim  -= binw;
        histITSTPC->GetAxis(2)->SetRangeUser(lowLim,upLim);
    }
    else histITSTPC->GetAxis(2)->SetRangeUser(minEta+binw,maxEta-binw);
    
    //phibins
    binw    = histITSTPC->GetAxis(3)->GetBinWidth(1)/2.;
    if(variable == kPhi) {
        lowLim += binw;
        upLim  -= binw;
        histITSTPC->GetAxis(3)->SetRangeUser(lowLim,upLim);
    }
    else histITSTPC->GetAxis(3)->SetRangeUser(minPhi+binw,maxPhi-binw);

    // select particleType
    double xmin,xmax;
    if(particleType==kElectron) {xmin=-0.5; xmax=0.5;}
    else if(particleType==kPion)     {xmin=0.5;  xmax=1.5;}
    else if(particleType==kKaon)     {xmin=1.5;  xmax=2.5;}
    else if(particleType==kProton)   {xmin=2.5;  xmax=3.5;}
    else if(particleType==kAll)      {xmin=3.5;  xmax=4.5;}
    else if(particleType==kPionKaon) {xmin=0.5;  xmax=2.5;}
    else {Printf("particleType not yet impelemented! Please provide it first!"); return 0x0;}
    binw    = histITSTPC->GetAxis(4)->GetBinWidth(1)/2.;
    histITSTPC->GetAxis(4)->SetRangeUser(xmin+binw,xmax-binw);
    
    //specie: primary, secondary, all
    binw    = histITSTPC->GetAxis(5)->GetBinWidth(1)/2.;
    if(species!=2)histITSTPC->GetAxis(5)->SetRangeUser(species-0.5+binw,species+0.5-binw);
    else histITSTPC->GetAxis(5)->SetRangeUser(-0.5+binw,1.5-binw);

    //TOF bc
    binw    = histITSTPC->GetAxis(6)->GetBinWidth(1)/2.;
    if(fIsSelectTOFBC)histITSTPC->GetAxis(6)->SetRangeUser(-0.5+binw,0.5-binw);

    //DCAxy
    binw    = histITSTPC->GetAxis(7)->GetBinWidth(1)/2.;
    histITSTPC->GetAxis(7)->SetRangeUser(-DCAxy+binw,DCAxy-binw);

    //DCAz
    binw    = histITSTPC->GetAxis(8)->GetBinWidth(1)/2.;
    histITSTPC->GetAxis(8)->SetRangeUser(-DCAz+binw,DCAz-binw);

    
    int projectionAxis=1;
    if(variable==kPt) projectionAxis = 1;
    else if(variable==kEta) projectionAxis = 2;
    else if(variable==kPhi) projectionAxis = 3;
    
    TH1D * hProj  = histITSTPC->Projection(projectionAxis);
    hProj->SetDirectory(0);
    hProj->SetNameTitle(Form("h%d",projectionAxis),Form("h%d",projectionAxis));

    
    delete l;
    inFileData->Close();
    return hProj;
}
