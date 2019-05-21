#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include "TH1D.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;
#endif

Double_t calcdelta(TH1D* data, TH1D* pri, TH1D* dcy, TH1D* mat, Double_t wdcy, Double_t wmat, Double_t exponent, Int_t method=0) {
    
    Double_t deltac = 0;
    for (int i=1; i<=data->GetNbinsX(); i++) {
//         if (TMath::Abs((pri->GetBinCenter(i) < 0.25))) continue;
        if (pri->GetBinContent(i) == 0) continue;
        if (dcy->GetBinContent(i) == 0) continue;
        if (mat->GetBinContent(i) == 0) continue;
        if (data->GetBinContent(i) == 0) continue;     
        Double_t ratio = (pri->GetBinContent(i)*(1.0-wdcy-wmat) + dcy->GetBinContent(i)*wdcy + mat->GetBinContent(i)*wmat) / data->GetBinContent(i);
        Double_t diff = pri->GetBinContent(i)*(1.0-wdcy-wmat) + dcy->GetBinContent(i)*wdcy + mat->GetBinContent(i)*wmat - data->GetBinContent(i);
        Double_t err = TMath::Sqrt(TMath::Power(pri->GetBinError(i)*(1.0-wdcy-wmat),2)+TMath::Power(dcy->GetBinError(i)*(wdcy),2)+TMath::Power(mat->GetBinError(i)*(wmat),2)+TMath::Power(data->GetBinError(i),2));
 
        if (method == 0) { deltac += TMath::Power(TMath::Abs(diff)/(err),exponent); } // default method ~chi2        
        if (method == 1) { deltac += TMath::Power(TMath::Abs(ratio-1.0),exponent); }   // only the difference to unity, no weight!
    }
    return deltac;
}

Double_t dcafits_3t_ptbins(Double_t ptmin, Double_t ptmax, Int_t ngroup = 1, Double_t exponent = 1.0, Double_t multmin = 2, Double_t multmax = 1000,  Int_t method=0) {
    
Double_t accuracy = 1e-4;

// Double_t exponent = 1.0;    
    
const char*   mcfilename = "/home/mknichel/charged-particle-spectra/train_output/LF_pp_MC/1111_LHC17l3b_cent/AnalysisResults.root";
const char* datafilename = "/home/mknichel/charged-particle-spectra/train_output/LF_pp/1247_LHC17p_pass1_CENT_wSDD/AnalysisResults.root";


// Double_t ptmin = 1.0;
// Double_t ptmax = 1.2;

Double_t binpt1 = 22;
Double_t binpt2 = 81;


// Double_t multmin = 2;
// Double_t multmax = 1000;

Double_t binmult1 = 1;
Double_t binmult2 = 10;


TFile * fmc = TFile::Open(mcfilename,"READ");
TFile * fdata = TFile::Open(datafilename,"READ");

THnSparseD* hmc   = (THnSparseD*)((TList*)((TDirectoryFile*)fmc->Get("TaskDCArStudy"))->Get("TaskDCArStudy"))->FindObject("fHistDCA");
THnSparseD* hdata = (THnSparseD*)((TList*)((TDirectoryFile*)fdata->Get("TaskDCArStudy"))->Get("TaskDCArStudy"))->FindObject("fHistDCA");

if ((ptmin>=0) && (ptmax > ptmin)) {
    binpt1 = hmc->GetAxis(1)->FindBin(ptmin);
    binpt2 = hmc->GetAxis(1)->FindBin(ptmax);
}

if ((multmin>=0) && (multmax > multmin)) {
    binmult1 = hmc->GetAxis(2)->FindBin(multmin);
    binmult2 = hmc->GetAxis(2)->FindBin(multmax);
}

hmc->GetAxis(1)->SetRange(binpt1,binpt2);
hmc->GetAxis(2)->SetRange(binmult1,binmult2);
hmc->GetAxis(3)->SetRangeUser(0,0);
TH1D* pri = (TH1D*) hmc->Projection(0)->Rebin(ngroup,"pri");
hmc->GetAxis(3)->SetRangeUser(1,2);
TH1D* sec = (TH1D*) hmc->Projection(0)->Rebin(ngroup,"sec");
hmc->GetAxis(3)->SetRangeUser(1,1);
TH1D* dcy = (TH1D*) hmc->Projection(0)->Rebin(ngroup,"dcy");
hmc->GetAxis(3)->SetRangeUser(2,2);
TH1D* mat = (TH1D*) hmc->Projection(0)->Rebin(ngroup,"mat");

hdata->GetAxis(1)->SetRange(binpt1,binpt2);
hdata->GetAxis(2)->SetRange(binmult1,binmult2);
hdata->GetAxis(3)->SetRange();
TH1D* data = (TH1D*) hdata->Projection(0)->Rebin(ngroup,"data");

pri->Sumw2();
sec->Sumw2();
dcy->Sumw2();
mat->Sumw2();
data->Sumw2();


Double_t fracmc = sec->Integral()/(pri->Integral()+sec->Integral());

Double_t fracmmc = mat->Integral()/(pri->Integral()+sec->Integral());
Double_t fracdmc = dcy->Integral()/(pri->Integral()+sec->Integral());


pri->Scale(1./pri->Integral());
sec->Scale(1./sec->Integral());
dcy->Scale(1./dcy->Integral());
mat->Scale(1./mat->Integral());
data->Scale(1./data->Integral());



pri->SetLineColor(kGreen);
sec->SetLineColor(kMagenta);
dcy->SetLineColor(kRed);
mat->SetLineColor(kBlue);
data->SetLineColor(kBlack);

//

Double_t wm = 0;
Double_t wd = 0;
Double_t wmmin = 0;
Double_t wdmin = 0;
Double_t deltamin = 1e20;

Double_t step=1./16.;

Double_t wm1 = 0;
Double_t wm2 = 1;

Double_t wd1 = 0;
Double_t wd2 = 1;

while (step>=accuracy) {
    for (wm = wm1; wm <= wm2; wm += step) {
        for (wd = wd1; wd <= wd2; wd += step) {
            Double_t d = calcdelta(data,pri,dcy,mat,wd,wm,exponent,method);
            if (d<deltamin) { deltamin=d; wdmin=wd; wmmin=wm; }
        }
    }
    wd1 = wdmin-5*step;
    if (wd1<0) wd1 = 0;
    wd2 = wdmin+5*step;
    if (wd2>1) wd2 = 1;
    wm1 = wmmin-5*step;
    if (wm1<0) wm1 = 0;
    wm2 = wmmin+5*step;
    if (wm2>1) wm2 = 1;
    step = step/2;    
    
}

Double_t dcaptmin = 0.0182+0.0350/TMath::Power(ptmin,1.01);
Double_t dcaptmax = 0.0182+0.0350/TMath::Power(ptmax,1.01);


cout<<"w DCY           = "<<wdmin<<endl;
cout<<"w MAT           = "<<wmmin<<endl;
cout<<"w               = "<<wdmin+wmmin<<endl;
cout<<"deltamin        = "<<deltamin<<endl;
cout<<"fracmc          = "<<fracmc<<endl;
cout<<"scalefactor DCY = "<<wdmin/fracdmc<<endl;
cout<<"scalefactor MAT = "<<wmmin/fracmmc<<endl;
cout<<"scalefactor     = "<<(wdmin+wmmin)/fracmc<<endl;

// cout<<"dca cut ptmin   = "<<-dcaptmin<<"..."<<dcaptmin<<endl;
// cout<<"dca cut ptmax   = "<<-dcaptmax<<"..."<<dcaptmax<<endl;



TH1D* fit = (TH1D*) pri->Clone("fit");
fit->Scale(1.-wdmin-wmmin);
fit->Add(dcy,wdmin);
fit->Add(mat,wmmin);


TH1D* fitratio = (TH1D*) fit->Clone("fitratio");

fitratio->Divide(data);

fit->SetLineColor(kOrange);
fit->SetMarkerColor(kOrange);
fitratio->SetLineColor(kBlack);
fitratio->SetMarkerColor(kBlack);
sec->SetLineColor(kMagenta);
sec->SetMarkerColor(kMagenta);
mat->SetLineColor(kBlue);
mat->SetMarkerColor(kBlue);
dcy->SetLineColor(kRed);
dcy->SetMarkerColor(kRed);
pri->SetLineColor(kGreen);
pri->SetMarkerColor(kGreen);
data->SetLineColor(kBlack);
data->SetMarkerColor(kBlack);


data->SetTitle("DATA");
fit->SetTitle("FIT (3 templates)");
pri->SetTitle("MC primaries");
sec->SetTitle("MC secondaries");
mat->SetTitle("MC secondaries from material");
dcy->SetTitle("MC secondaries from decay");
fitratio->SetTitle("FIT/DATA");

gROOT->SetStyle("Modern");
gStyle->SetOptStat(kFALSE);
gStyle->SetOptTitle(kTRUE);
TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
c1->Divide(0,2,0.01,0.00);

c1->cd(1);
gPad->SetLogy();
gPad->SetTicks(1,1);
pri->GetYaxis()->SetRangeUser(2e-7,9.9);
pri->GetYaxis()->SetTitle("");
pri->Draw("");
sec->Draw(" SAME");
mat->Draw(" SAME");
dcy->Draw(" SAME");
fit->Draw(" SAME");
data->Draw(" SAME");

TLegend* l1 = new TLegend(0.45,0.04,0.66,0.28);
l1->AddEntry(data);
l1->AddEntry(fit);
l1->AddEntry(pri);
l1->AddEntry(sec);
l1->AddEntry(mat);
l1->AddEntry(dcy);
l1->Draw();



c1->cd(2);
gPad->SetTicks(1,1);
fitratio->GetYaxis()->SetRangeUser(0,2.1);
fitratio->Draw();
fitratio->GetYaxis()->SetTitle("FIT/DATA");
fitratio->GetXaxis()->SetTitle("DCA_{xy} (cm)");
TF1* one = new TF1("one","1",-1,1);
one->SetLineColor(kBlack);
one->SetLineWidth(1);
one->Draw("SAME");

/*
TLegend* l2 = new TLegend(0.38,0.04,0.67,0.30);
l2->AddEntry(data);
l2->AddEntry(fit);
l2->AddEntry(pri);
l2->AddEntry(sec);
l2->AddEntry(mat);
l2->AddEntry(dcy);
l2->Draw();
*/

Double_t ssf = (wdmin+wmmin)/fracmc;

TString s1 = "";
s1 += multmin;
s1 += " < Nacc < ";
s1 += multmax;
s1 += ",   ";
s1 += Form("%g",ptmin);
s1 += " < pT < ";
s1 += Form("%g",ptmax);
s1 += ",   equalweights=";
s1 += method;

TString s2 = "resulting scaling factor = ";
s2 += ssf;

pri->SetTitle(s1.Data());
fitratio->SetTitle(s2.Data());

TString filename= "plots/dcafit_";
filename += multmin;
filename += "nacc";
filename += multmax;
filename += "_";
filename += Form("%g",ptmin);
filename += "pt";
filename += Form("%g",ptmax);
filename += "_method";
filename += method;

filename.ReplaceAll(".","");

cout<<filename.Data()<<endl;
cout<<filename.Data()<<endl;


c1->SaveAs(filename + ".root");
c1->SaveAs(filename + ".png");

return (wdmin+wmmin)/fracmc;

}
