//________________________________________________________________________________________________//
// Brief: Macro for correct observed vn for prompt fraction and application of syst unc           //
// Main Function: ComputePromptVnWithSyst                                                         //
// Inputs:                                                                                        //
// - inFileNameVn -> output of CharmHadronVnAnalysis.C                                            //
// - inFileNamefPrompt -> output of ComputeDmesonYield.C                                          //
// Author: Fabrizio Grosa, fabrizio.grosa@cern.ch                                                 //
//________________________________________________________________________________________________//

#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>

#endif

using namespace std;

//________________________________________________________________________________________________//
//systematic uncertainties
const double absSystUncFit[] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02};
const double relResolUnc = 0.01;

//________________________________________________________________________________________________//
//method prototypes
void ComputePromptVnWithSyst(TString outFileName, TString inFileNameVn, TString graphNameVn="gvnSimFit", TString inFileNamefPrompt = "");
TGraphAsymmErrors* CorrectVnForFeedDown(TGraphAsymmErrors *&gVnPromptSyst, TGraphAsymmErrors* gVnObsStat, TH1F* hfPromptCent, TH1F* hfPromptMin, TH1F* hfPromptMax);
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, bool dofill=false, float markersize = 1.5, int linewidth = 2);
void SetStyle();

//________________________________________________________________________________________________//
//main function to get prompt vn with all syst unc
void ComputePromptVnWithSyst(TString outFileName, TString inFileNameVn, TString graphNameVn, TString inFileNamefPrompt) {

    SetStyle();

    const int nPtBins = sizeof(absSystUncFit) / sizeof(absSystUncFit[0]);

    TFile* inFileVn = TFile::Open(inFileNameVn.Data());
    if(!inFileVn || !inFileVn->IsOpen())
        return;
    TGraphAsymmErrors* gVnObsStat = static_cast<TGraphAsymmErrors*>(inFileVn->Get(graphNameVn.Data()));
    if(!gVnObsStat) {
        cerr << Form("ERROR: graph %s does not exist! Exit",graphNameVn.Data()) << endl;
        return;
    }
    gVnObsStat->SetName("gVnObsStat");
    SetGraphStyle(gVnObsStat,kRed,kOpenCircle);
    if(nPtBins!=gVnObsStat->GetN()) {
        cerr << "ERROR: different number of pT bins in systematic-uncertainty array and input vn graph! Exit" << endl;
        return;
    }

    //Apply prompt correction, if fprompt file is passed
    TH1F *hfPromptCent = NULL, *hfPromptMin = NULL, *hfPromptMax = NULL;
    TGraphAsymmErrors *gVnPromptStat = NULL, *gVnPromptSystFeedDown = NULL, *gfPrompt = NULL;
    TFile* inFilefPrompt = NULL;
    if(inFileNamefPrompt!="") {
        inFilefPrompt = TFile::Open(inFileNamefPrompt.Data());
        if(!inFilefPrompt || !inFilefPrompt->IsOpen())
            return;
        hfPromptCent = static_cast<TH1F*>(inFilefPrompt->Get("hfPromptCent"));
        hfPromptMin = static_cast<TH1F*>(inFilefPrompt->Get("hfPromptMinNb"));
        hfPromptMax = static_cast<TH1F*>(inFilefPrompt->Get("hfPromptMaxNb"));
        hfPromptCent->SetDirectory(0);
        hfPromptMin->SetDirectory(0);
        hfPromptMax->SetDirectory(0);
        if(gVnObsStat->GetN()!=hfPromptCent->GetNbinsX()) {
            cerr << "ERROR: different number of pT bins in input vn graph and fprompt histo! Exit" << endl;
            return;
        }
        gVnPromptStat = CorrectVnForFeedDown(gVnPromptSystFeedDown, gVnObsStat,hfPromptCent,hfPromptMin,hfPromptMax);
        gVnPromptStat->SetName("gVnPromptStat");
        gVnPromptSystFeedDown->SetName("gVnPromptSystFeedDown");

        gfPrompt = new TGraphAsymmErrors(0);
        gfPrompt->SetName("gfPrompt");
        SetGraphStyle(gfPrompt,kBlack,kFullCircle);
        for(int iPt=0; iPt<nPtBins; iPt++) {
            gfPrompt->SetPoint(iPt,hfPromptCent->GetBinCenter(iPt+1),hfPromptCent->GetBinContent(iPt+1));
            gfPrompt->SetPointError(iPt,hfPromptCent->GetBinWidth(iPt+1)/2,hfPromptCent->GetBinWidth(iPt+1)/2,hfPromptCent->GetBinContent(iPt+1)-hfPromptMin->GetBinContent(iPt+1),hfPromptMax->GetBinContent(iPt+1)-hfPromptCent->GetBinContent(iPt+1));
        }
    }
    else {
        gVnPromptStat = static_cast<TGraphAsymmErrors*>(gVnObsStat->Clone("gVnPromptStat"));
        gVnPromptSystFeedDown = static_cast<TGraphAsymmErrors*>(gVnObsStat->Clone("gVnPromptSystFeedDown"));
        for(int iPt=0; iPt<nPtBins; iPt++)
            gVnPromptSystFeedDown->SetPointError(iPt,0.3,0.3,0.,0.);
    }
    SetGraphStyle(gVnPromptStat,kBlack,kFullCircle);
    SetGraphStyle(gVnPromptSystFeedDown,kGray+2,kFullCircle,true);

    //compute fit syst and resol syst on vn prompt
    TGraphAsymmErrors *gVnPromptSystFit = static_cast<TGraphAsymmErrors*>(gVnPromptStat->Clone("gVnPromptSystFit"));
    TGraphAsymmErrors *gVnPromptSystResol = static_cast<TGraphAsymmErrors*>(gVnPromptStat->Clone("gVnPromptSystResol"));
    TGraphAsymmErrors *gVnPromptSystData = static_cast<TGraphAsymmErrors*>(gVnPromptStat->Clone("gVnPromptSystData"));
    TGraphAsymmErrors *gVnPromptSystTot = static_cast<TGraphAsymmErrors*>(gVnPromptStat->Clone("gVnPromptSystTot"));

    double ptmin = -1., ptmax = -1.;
    for(int iPt=0; iPt<nPtBins; iPt++) {
        double pt=-1., vnobs=-1., vnprompt=-1.;
        gVnPromptStat->GetPoint(iPt,pt,vnprompt);
        gVnObsStat->GetPoint(iPt,pt,vnobs);
        double scalefactor = vnprompt/vnobs;
        double fitsyst = absSystUncFit[iPt]*scalefactor;
        double resolsyst = relResolUnc*vnprompt;
        double feeddownsystlow = gVnPromptSystFeedDown->GetErrorYlow(iPt);
        double feeddownsysthigh = gVnPromptSystFeedDown->GetErrorYhigh(iPt);
        double datasyst = TMath::Sqrt(fitsyst*fitsyst+resolsyst*resolsyst);
        double totsystlow = TMath::Sqrt(datasyst*datasyst+feeddownsystlow*feeddownsystlow);
        double totsysthigh = TMath::Sqrt(datasyst*datasyst+feeddownsysthigh*feeddownsysthigh);
        gVnPromptSystFit->SetPointError(iPt,0.5,0.5,fitsyst,fitsyst);
        gVnPromptSystResol->SetPointError(iPt,0.5,0.5,resolsyst,resolsyst);
        gVnPromptSystData->SetPointError(iPt,0.5,0.5,datasyst,datasyst);
        gVnPromptSystTot->SetPointError(iPt,0.5,0.5,totsystlow,totsysthigh);
        if(iPt==0)
            ptmin = pt-gVnObsStat->GetErrorXlow(iPt);
        else if(iPt==nPtBins-1)
            ptmax = pt+gVnObsStat->GetErrorXhigh(iPt);
    }

    //plots
    TLine* lineatzero = new TLine(ptmin, 0., ptmax, 0.);
    lineatzero->SetLineWidth(2);
    lineatzero->SetLineColor(kBlack);
    lineatzero->SetLineStyle(9);

    TLegend* legPrompt = new TLegend(0.5,0.75,0.85,0.85);
    legPrompt->SetFillStyle(0);
    legPrompt->SetTextSize(0.04);
    legPrompt->SetBorderSize(0);
    legPrompt->AddEntry(gVnObsStat,"Observed #it{v}_{2}","p");
    legPrompt->AddEntry(gVnPromptStat,"Prompt #it{v}_{2}","p");

    TLegend* leg = new TLegend(0.5,0.7,0.85,0.85);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->AddEntry(gVnPromptStat,"Prompt #it{v}_{2}","p");
    leg->AddEntry(gVnPromptSystData,"Data syst. unc.","f");
    leg->AddEntry(gVnPromptSystFeedDown,"B feed-down syst. unc.","f");

    TCanvas* cVnPrompt = new TCanvas("cVnPrompt","",800,800);
    cVnPrompt->DrawFrame(ptmin,-0.2,ptmax,0.4,";#it{p}_{T} (GeV/#it{c});#it{v}_{2}");
    lineatzero->Draw("same");
    gVnObsStat->Draw("P");
    gVnPromptStat->Draw("P");
    legPrompt->Draw();

    TCanvas* cVn = new TCanvas("cVn","",800,800);
    cVn->DrawFrame(ptmin,-0.2,ptmax,0.4,";#it{p}_{T} (GeV/#it{c});#it{v}_{2}");
    lineatzero->Draw("same");
    gVnPromptSystFeedDown->Draw("2");
    gVnPromptSystData->Draw("2");
    gVnPromptStat->Draw("P");
    leg->Draw();

    //output files
    TFile outFile(outFileName.Data(),"recreate");
    gVnObsStat->Write();
    gVnPromptStat->Write();
    gVnPromptSystFit->Write();
    gVnPromptSystResol->Write();
    gVnPromptSystData->Write();
    gVnPromptSystFeedDown->Write();
    gVnPromptSystTot->Write();
    if(gfPrompt) gfPrompt->Write();
    outFile.Close();
}

//________________________________________________________________________________________________//
//method to compute prompt vn and associated syst unc
TGraphAsymmErrors* CorrectVnForFeedDown(TGraphAsymmErrors *&gVnPromptSyst, TGraphAsymmErrors* gVnObsStat, TH1F* hfPromptCent, TH1F* hfPromptMin, TH1F* hfPromptMax) {

    TGraphAsymmErrors* gVnPromptStat = new TGraphAsymmErrors(0);
    gVnPromptSyst = new TGraphAsymmErrors(0);

    for(int iPt=0; iPt<gVnObsStat->GetN(); iPt++) {
        double pt=-1., vnobs=-1.;
        gVnObsStat->GetPoint(iPt,pt,vnobs);
        double vnobsstatunc = gVnObsStat->GetErrorYlow(iPt);
        double ptunclow = gVnObsStat->GetErrorXlow(iPt);
        double ptunchigh = gVnObsStat->GetErrorXhigh(iPt);

        double fpromptcent = hfPromptCent->GetBinContent(iPt+1);
        double fpromptmin = hfPromptMin->GetBinContent(iPt+1);
        double fpromptmax = hfPromptMax->GetBinContent(iPt+1);

        double vnprompt = 2*vnobs/(1+fpromptcent);
        double vnpromptstatunc = 2*vnobsstatunc/(1+fpromptcent);
        double vnpromptmin = 2*TMath::Sqrt(3)*vnobs/((TMath::Sqrt(3)+1)+fpromptmax*(TMath::Sqrt(3)-1));
        double vnpromptmax = 2*TMath::Sqrt(3)*vnobs/((TMath::Sqrt(3)-1)+fpromptmin*(TMath::Sqrt(3)+1));

        gVnPromptStat->SetPoint(iPt,pt,vnprompt);
        gVnPromptStat->SetPointError(iPt,ptunclow,ptunchigh,vnpromptstatunc,vnpromptstatunc);
        gVnPromptSyst->SetPoint(iPt,pt,vnprompt);
        gVnPromptSyst->SetPointError(iPt,0.3,0.3,vnprompt-vnpromptmin,vnpromptmax-vnprompt);
    }

    return gVnPromptStat;
}

//___________________________________________________________________________________//
//method to set plots style
void SetGraphStyle(TGraphAsymmErrors* graph, int color, int markerstyle, bool dofill, float markersize, int linewidth) {
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    if(dofill) {
        graph->SetFillColor(color);
        graph->SetFillStyle(1000);
    }
    else {
        graph->SetFillStyle(0);
    }
    graph->SetMarkerStyle(markerstyle);
    graph->SetMarkerSize(markersize);
    graph->SetLineWidth(linewidth);
}

//________________________________________________________________________________________________//
//method to set plots style
void SetStyle() {
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetTitleSize(0.045,"xy");
    gStyle->SetLabelSize(0.04,"xy");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
}
