#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TStyle.h>
#include "AdditionalFunctions.h"
#include "PlotXi1530.C"
#include "SystematicHelper.cxx"
#include "YieldMean.C"

TString path = "./";
TString workdirectory = "data/";
TString finalfile;
bool isINEL = false;
// Default Canvas
Double_t w = 1920;
Double_t h = 1080;
double zero = 1e-25;
vector<int> vcolors;
vector<double> ptbin = {0.0, 0.8, 1.2, 1.6, 2.0, 2.4,
                        3.2, 4.0, 4.8, 5.6, 8.8, 15};
vector<vector<double>> multibincheck = {{0, 100}, {0, 10},  {10, 30},
                                        {30, 50}, {50, 70}, {70, 100}, {0, 0.01}, {0.01, 0.05}, {0.05,0.1}};
vector<int> Markerset = {20, 21, 33, 34, 29, 24, 25, 26, 27, 28,
                         30, 3,  5,  42, 43, 46, 47, 48, 49, 50,
                         51, 20, 21, 33, 34, 29, 24, 25, 26};
vector<int> Markerset_size = {1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
enum {
    kFitExpPt,
    kFitLevi,
    fFitExpMt,
    kFitBoltzmann,
    kFitBlastWave,
    kFitBoseEinstein,
    kFitFermiDirac
};
vector<TString> fitfuctions = {"fLevi", "fMtExp",        "fBoltzmann",
                               "fBGBW", "fBoseEinstein", "fFermiDirac"};
vector<TString> exhistograms = {"extremehard_max", "extremehard_min",
                                "extremehigh", "extremelow"};
// functions
TCanvas* plotHistsAndRatio(vector<TH1D*> numeratorHistograms,
                           TH1D* denominatorHist,
                           TString title = "",
                           TString xTitle = "",
                           TString yTitle = "",
                           TString options = "");
TCanvas* plotHistsAndRatio(TH1D* numeratorHist,
                           TH1D* denominatorHist,
                           TString title = "",
                           TString xTitle = "",
                           TString yTitle = "",
                           TString options = "");
vector<double> GetTrigEfficiency(double multi_start, double multi_end);
vector<double> GetVertexEfficiency(double multi_start, double multi_end);
TH1D* GetSigLoss(double multi_start, double multi_end);
TH1D* GetMCEfficiency(double multi_start, double multi_end);
TH1D* GetMCEfficiencyfromName(TString inputfilename);
TH1D* GetSpectra(double multi_start, double multi_end);
TH1D* GetSpectrasys(double multi_start, double multi_end);
TH1D* GetSpectrastat(double multi_start, double multi_end);
TH1D* GetSysError(double multi_start, double multi_end);
TF1* GetFittedFunction(double multi_start, double multi_end, int fitfunction);
TH1* GetExtraSpectra(double multi_start, double multi_end, int histooption);
TString GetBinInfo(double multis, double multie);
double myLevyPtFunc(Double_t* x, Double_t* par);
void PlotMiscXi1530() {
    TLatex* t_big = new TLatex();
    t_big->SetNDC();
    t_big->SetTextSize(0.07);

    // for memo, small
    TLatex* t = new TLatex();
    t->SetNDC();
    t->SetTextSize(0.04);

    TLatex* t0 = new TLatex();
    t0->SetNDC();
    t0->SetTextSize(0.030);

    vector<double> multibinsloop = {0, 10, 30, 50, 70, 100};
     // vector<double> multibinsloop = {0, 10, 30, 50, 100};

    auto hDef = GetMCEfficiency(0, 100);
    vector<TH1D*> buffer;
    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        auto temp_multi = GetMCEfficiency(multibinsloop[imultibin],
                                          multibinsloop[imultibin + 1]);
        buffer.push_back(temp_multi);
    }

    TCanvas* ratio_mult = plotHistsAndRatio(buffer, hDef, "", "p_{T}(GeV/c)",
                                            "Acceptance x Efficiency x B.R.",
                                            "NLOG Efficiency FINAL fine");

    ratio_mult->cd(1);
    auto legend_effi_ratio = new TLegend(.18, .01, .5, .20);
    legend_effi_ratio->SetNColumns(2);
    legend_effi_ratio->SetBorderSize(0);
    legend_effi_ratio->SetFillStyle(0);
    legend_effi_ratio->AddEntry(hDef, "MB(0-100)", "PL");

    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        buffer.at(imultibin)->SetFillColor(vcolors[imultibin]);
        buffer.at(imultibin)->SetLineColor(vcolors[imultibin]);
        buffer.at(imultibin)->SetMarkerColor(vcolors[imultibin]);
        buffer.at(imultibin)->SetMarkerStyle(Markerset[imultibin]);
        buffer.at(imultibin)->SetMarkerSize(Markerset_size[imultibin]);
    }

    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        legend_effi_ratio->AddEntry(
            buffer.at(imultibin),
            Form("%.0f - %.0f", multibinsloop.at(imultibin),
                 multibinsloop.at(imultibin + 1)),
            "PL");
    }
    legend_effi_ratio->Draw();

    ratio_mult->SaveAs(Form("%s1_Multi_0.00-100.00_Default1/MCEffi_Multi.pdf",
                            workdirectory.Data()));
    // Temp
    TCanvas* cMultiVert = new TCanvas("cMultiVert", "cMultiVert", 540, 960);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cMultiVert->SetTickx();
    cMultiVert->SetLogy();
    cMultiVert->SetTicky();
    cMultiVert->SetTopMargin(0.05);
    cMultiVert->SetLeftMargin(0.15);
    // cSigbkg->SetBottomMargin(0.01);
    cMultiVert->SetRightMargin(0.01);
    cMultiVert->SetFillStyle(0);
    cMultiVert->Draw();

    hDef->Draw();
    cMultiVert->SaveAs(Form("%s1_Multi_0.00-100.00_Default1/MCEffi_Vert.pdf",
                            workdirectory.Data()));


    //---------------------------------------------
    // Recon.effi. trend
    vector<vector<double>> multibin = {
        {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    vector<double> dNdetaAxis; //[multibin]
    vector<double> dNdetaAxis_e; //[multibin]
    vector<double> dNdetaAxis_e_half; //[multibin]
    vector<double> zeroerror; // [multibin]
    vector<vector<double>> Efficiency; // [ptbin][multibin]
    vector<vector<double>> Efficiency_state; //[ptbin][multibin]
    vector<TString> fMulti_text;

    int nptbins = hDef->GetNbinsX();
    int pTColorPallet = GetSerialColors(nptbins);
    cout << "# of pt bins:" << nptbins << endl;

    // 1. Get dNdeta(x axis) vectors
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        cout << "mutlibin" << imultibin << ", " << multibin[imultibin][0] << " - " << multibin[imultibin][1] << endl;

        vector<double> dndetav_temp  = GetdNdetawithError(multibin[imultibin][0],multibin[imultibin][1]);
        dNdetaAxis.push_back(dndetav_temp[0]);
        dNdetaAxis_e.push_back(dndetav_temp[1]);
        dNdetaAxis_e_half.push_back(dndetav_temp[1]/2);

        // text
        TString text_temp = Form("#bf{%.0f-%.0f}", multibin[imultibin][0],multibin[imultibin][1]);
        fMulti_text.push_back(text_temp);
    }
    dNdetaAxis.push_back(35.3);
    dNdetaAxis_e.push_back(1e-3);
    dNdetaAxis_e.push_back(5e-4);
    // 2. Get Efficiency vectors
    for (int i = 0; i < nptbins-1; i++) {
        if(i == 0) continue;
        Efficiency.push_back({});
        Efficiency_state.push_back({});
        for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
            auto temp_multi = GetMCEfficiency(multibin[imultibin][0],multibin[imultibin][1]);
            cout << "ptbin " << i-1 << ", mutlibin" << imultibin << ", " << multibin[imultibin][0] << " - " << multibin[imultibin][1] <<  ", Efficiency " << temp_multi->GetBinContent(i+1) << " +-" << temp_multi->GetBinError(i+1) << endl;
            Efficiency[i-1].push_back(temp_multi->GetBinContent(i+1));
            Efficiency_state[i-1].push_back(temp_multi->GetBinError(i+1));
        }
    }
    // 3. Make vecor of effi and legend
    vector<TGraphErrors*> gMulti_effi;
    vector<TF1*> fMulti_effi;
    vector<TF1*> fMulti_effi_MB;

    auto legend_effi_mult = new TLegend(.15, .15, .6, .4);
    legend_effi_mult->SetNColumns(2);
    legend_effi_mult->SetBorderSize(0);
    legend_effi_mult->SetFillStyle(0);

    // 4. Make graph and add into vectors
    for (int i = 0; i < nptbins-2; i++) {
        // graph
        TGraphErrors* g_temp = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &Efficiency[i][0], &dNdetaAxis_e[0], &Efficiency_state[i][0]);
        g_temp->SetMaximum(4e-1);
        g_temp->SetMinimum(4e-4);
        g_temp->GetXaxis()->SetLimits(0,45);
        g_temp->SetLineColor(pTColorPallet+i);
        g_temp->SetMarkerColor(pTColorPallet+i);
        g_temp->SetMarkerStyle(20);
        g_temp->SetMarkerSize(1.5);
        g_temp->GetYaxis()->SetTitle(
        "Rec.Efficiency of #Xi(1530)^{0}");
        g_temp->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
        gMulti_effi.push_back(g_temp);

        // fit
        TF1* ftemp = new TF1("linear fit", "[0]*x+[1]", 0, 45);
        g_temp->Fit(ftemp, "RQN");
        ftemp->SetLineWidth(1);
        ftemp->SetLineStyle(2);
        ftemp->SetLineColor(pTColorPallet+i);
        fMulti_effi.push_back(ftemp);
        legend_effi_mult->AddEntry(gMulti_effi[i], Form("#it{p}_{T} %.2f-%.2f, Slope: %.2f#times10^{-4}",ptbin[i],ptbin[i+1],ftemp->GetParameter(0)*1e4), "PL");

        double tempvalue =ftemp->GetParameter(0)*35+ftemp->GetParameter(1);

        cout << "ptbin: " << i << " MB effi: " << hDef->GetBinContent(i+2) << ", [70-100] value: " << Efficiency[i][multibin.size()-1] << " (" << Form("%.2f%%",(hDef->GetBinContent(i+2)-Efficiency[i][multibin.size()-1])/Efficiency[i][multibin.size()-1]*100) << ")" <<", [0, 0.01](35.3) value:" << tempvalue << " (" << Form("%.2f%%",(hDef->GetBinContent(i+2)-tempvalue)/hDef->GetBinContent(i+2)*100) << ")" << endl;
        Efficiency[i].push_back(tempvalue);
        Efficiency_state[i].push_back(1e-10);

    }

    TCanvas* cMultiEffi = new TCanvas("cMultiEffi", "cMultiEffi", w, h);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cMultiEffi->SetTickx();
    cMultiEffi->SetLogy();
    cMultiEffi->SetTicky();
    cMultiEffi->SetTopMargin(0.05);
    cMultiEffi->SetLeftMargin(0.10);
    // cSigbkg->SetBottomMargin(0.01);
    cMultiEffi->SetRightMargin(0.01);
    cMultiEffi->SetFillStyle(0);
    cMultiEffi->Draw();

    // 5. Draw togther
    gMulti_effi[0]->Draw("AP");
    fMulti_effi[0]->Draw("same");

    for (int i = 1; i < nptbins-2; i++) {
        gMulti_effi[i]->Draw("P");
        fMulti_effi[i]->Draw("same");
    }
    legend_effi_mult->Draw();

    // multiplicity text
    t0->DrawLatex(0.12, 0.88,  fMulti_text[multibin.size()-1].Data());
    t0->DrawLatex(0.17, 0.88,  fMulti_text[multibin.size()-2].Data());
    t0->DrawLatex(0.225, 0.88,  fMulti_text[multibin.size()-3].Data());
    t0->DrawLatex(0.31, 0.88,  fMulti_text[multibin.size()-4].Data());
    t0->DrawLatex(0.455, 0.88,  fMulti_text[multibin.size()-5].Data());

    // HM line text
    TLine* lHM001 = new TLine(35.3, 4e-4, 35.3, 4e-1);
    lHM001->SetLineStyle(2);
    lHM001->Draw();
    TString text_HM = "0-0.01 bin";
    t0->DrawLatex(0.805, 0.88,  text_HM.Data());

    cMultiEffi->SaveAs(Form("%s1_Multi_0.00-100.00_Default1/MCEffi_Multi_trend.pdf",
                            workdirectory.Data()));

    cMultiEffi->SetLogy(false);
    // Additional
    vector<TGraphErrors*> gMulti_effi_HM;
    
    for (int i = 1; i < nptbins-1; i++) {
        // graph
        TGraphErrors* g_temp = new TGraphErrors(multibin.size()+1, &dNdetaAxis[0], &Efficiency[i-1][0], &dNdetaAxis_e[0], &Efficiency_state[i-1][0]);
        g_temp->GetXaxis()->SetLimits(0,45);
        g_temp->SetLineColor(pTColorPallet+i);
        g_temp->SetMarkerColor(pTColorPallet+i);
        g_temp->SetMarkerStyle(20);
        g_temp->SetMarkerSize(1.5);
        g_temp->GetYaxis()->SetTitle(
        "Rec.Efficiency of #Xi(1530)^{0}");
        g_temp->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
        gMulti_effi_HM.push_back(g_temp);

        TF1* templine =
            new TF1("templine", Form("%f",hDef->GetBinContent(i+1)), -1, 45);

        g_temp->Draw("AP");
        templine->Draw("same");
        fMulti_effi[i-1]->Draw("same");

        cMultiEffi->SaveAs(Form("%s1_Multi_0.00-100.00_Default1/MCEffi_Multi_trend_bin%i.pdf",
                            workdirectory.Data(),i));
    }

    //---------------------------------------------

    TFile* inputfile_MB = new TFile(
        Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default1.root",
             workdirectory.Data()));
    auto hXiSpectrum = (TH1D*)inputfile_MB->Get("hXiSpectrum");
    hXiSpectrum->SetMaximum(5e-3);
    hXiSpectrum->SetMinimum(1e-6);
    hXiSpectrum->GetXaxis()->SetRangeUser(0.5, 9.0);
    hXiSpectrum->SetFillColor(kRed);
    hXiSpectrum->SetLineColor(kRed);
    hXiSpectrum->SetMarkerColor(kRed);

    // 7TeV Results
    // Old data
    // from HEPDATA https://www.hepdata.net/record/ins1300380
    vector<double> CorrectedYeild_7TeV = {1e-10,   0.00138,  0.00109, 0.00078,
                                          0.00046, 0.000226, 6.6e-05, 2.34e-05,
                                          8.1e-06, 1e-10,    1e-10};
    vector<double> CorrectedYeild_syserr_7TeV = {
        1e-10,    0.00017,  6.0e-05,  5.0e-05, 2.2e-05, 1.15e-05,
        3.73e-06, 1.96e-06, 4.54e-07, 1e-10,   1e-10};
    vector<double> CorrectedYeild_staterr_7TeV = {
        1e-10,   4.0e-05, 3.0e-05, 2.0e-05, 1.5e-05, 3.38e-06,
        2.08e-6, 8.11e-7, 9.09e-7, 1e-10,   1e-10};

    auto hSpectra_7TeV_syserr = MakeHistfromArray(
        "7TeV Spectra with systematic error", CorrectedYeild_7TeV,
        ptbin, CorrectedYeild_syserr_7TeV);
    auto hSpectra_7TeV_staterr = MakeHistfromArray(
        "7TeV Spectra with statistical error", CorrectedYeild_7TeV, ptbin,
        CorrectedYeild_staterr_7TeV);

    //---------------------------------------------
    vector<TH1D*> buffer_multi;
    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        auto temp_spectra_multi =
            GetSpectra(multibinsloop[imultibin], multibinsloop[imultibin + 1]);
        buffer_multi.push_back(temp_spectra_multi);
    }


    bool mceffi_period = false;
    if (mceffi_period) {
        auto tempreceffi_default = (TH1D*)GetMCEfficiencyfromName(Form(
            "%s/AnalysisResults_Extracted_1_Multi_0.00-100.00_Default1.root",
            path.Data()));

        auto tempreceffi_LHC18c6a_1 = (TH1D*)GetMCEfficiencyfromName(Form(
            "%sAnalysisResults_Xi1530MCQA_1_Multi_0.00-100.00_LHC18c6a_part1.root",
            path.Data()));
        auto tempreceffi_LHC18c6a_2 = (TH1D*)GetMCEfficiencyfromName(Form(
            "%sAnalysisResults_Xi1530MCQA_1_Multi_0.00-100.00_LHC18c6a_part2.root",
            path.Data()));
        auto tempreceffi_LHC18c6b_1 = (TH1D*)GetMCEfficiencyfromName(Form(
            "%sAnalysisResults_Xi1530MCQA_1_Multi_0.00-100.00_LHC18c6b_part1.root",
            path.Data()));
        auto tempreceffi_LHC18c6b_2 = (TH1D*)GetMCEfficiencyfromName(Form(
            "%sAnalysisResults_Xi1530MCQA_1_Multi_0.00-100.00_LHC18c6b_part2.root",
            path.Data()));
        auto tempreceffi_LHC18c6c = (TH1D*)GetMCEfficiencyfromName(Form(
            "%sAnalysisResults_Xi1530MCQA_1_Multi_0.00-100.00_LHC18c6c_novertexer.root",
            path.Data()));
        auto tempreceffi_LHC18c6b4 = (TH1D*)GetMCEfficiencyfromName(Form(
            "%sAnalysisResults_Xi1530MCQA_1_Multi_0.00-100.00_LHC18c6b4_novertexer.root",
            path.Data()));
        
        tempreceffi_default->SetFillColor(kBlack);
        tempreceffi_default->SetLineColor(kBlack);
        tempreceffi_default->SetMarkerColor(kBlack);
        tempreceffi_default->SetMarkerStyle(Markerset[0]);
        tempreceffi_default->SetMarkerSize(Markerset_size[0]);

        vector<TH1D*> periodeffi;
        periodeffi.push_back(tempreceffi_LHC18c6a_1);
        periodeffi.push_back(tempreceffi_LHC18c6a_2);
        periodeffi.push_back(tempreceffi_LHC18c6b_1);
        periodeffi.push_back(tempreceffi_LHC18c6b_2);
        periodeffi.push_back(tempreceffi_LHC18c6c);
        periodeffi.push_back(tempreceffi_LHC18c6b4);

        TCanvas* ratio_period = plotHistsAndRatio(
            periodeffi, tempreceffi_default, "", "p_{T}(GeV/c)",
            "Acceptance x Efficiency x B.R.", "NLOG Efficiency FINAL");

        auto legend_effi_period = new TLegend(.18, .015, .8, .20);
        // legend_effi_vertexer->SetNColumns(2);
        legend_effi_period->SetBorderSize(0);
        legend_effi_period->SetFillStyle(0);
        legend_effi_period->AddEntry(tempreceffi_default,
                                       "Default(18c6b_1)", "PL");
        int runs = 1;
        vector<TString> testnames = {"18c6a_1","18c6a_2","18c6b_1", "18c6b_2", "18c6c", "18c6b4"};
        for (auto const& hTemp_vertex : periodeffi) {
            hTemp_vertex->SetFillColor(vcolors[runs-1]);
            hTemp_vertex->SetLineColor(vcolors[runs-1]);
            hTemp_vertex->SetMarkerColor(vcolors[runs-1]);
            hTemp_vertex->SetMarkerStyle(Markerset[runs-1]);
            hTemp_vertex->SetMarkerSize(Markerset_size[runs-1]);
            legend_effi_period->AddEntry(hTemp_vertex,
                                       Form("Vertexer with %s", testnames[runs-1].Data()), "PL");
            runs++;
        }
        legend_effi_period->Draw();


        ratio_period->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/MCEffi_period.pdf",
                 workdirectory.Data()));

    }

    bool mceffi_vertexer = false;
    if (mceffi_vertexer) {
        auto tempreceffi_novertexer = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "NoVertexer1.root",
                 path.Data()));
        auto tempreceffi_vertexer_test1 = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Allon1.root",
                 path.Data()));

        auto tempreceffi_vertexer_test2 = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Extraoff1.root",
                 path.Data()));
        auto tempreceffi_vertexer_test3 = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Preseloff1.root",
                 path.Data()));
        auto tempreceffi_vertexer_test4 = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Alloff1.root",
                 path.Data()));

        vector<TH1D*> vertexertests;
        vertexertests.push_back(tempreceffi_vertexer_test1);
        vertexertests.push_back(tempreceffi_vertexer_test2);
        vertexertests.push_back(tempreceffi_vertexer_test3);
        vertexertests.push_back(tempreceffi_vertexer_test4);
        
        //vertexertests.push_back(tempreceffi_vertexer_test3);
        tempreceffi_novertexer->SetFillColor(kBlack);
        tempreceffi_novertexer->SetLineColor(kBlack);
        tempreceffi_novertexer->SetMarkerColor(kBlack);
        tempreceffi_novertexer->SetMarkerStyle(Markerset[0]);
        tempreceffi_novertexer->SetMarkerSize(Markerset_size[0]);

        TCanvas* ratio_vertexer = plotHistsAndRatio(
            vertexertests, tempreceffi_novertexer, "", "p_{T}(GeV/c)",
            "Acceptance x Efficiency x B.R.", "NLOG Efficiency FINAL");

        ratio_vertexer->cd(1);
        auto legend_effi_vertexer = new TLegend(.18, .015, .8, .20);
        // legend_effi_vertexer->SetNColumns(2);
        legend_effi_vertexer->SetBorderSize(0);
        legend_effi_vertexer->SetFillStyle(0);
        legend_effi_vertexer->AddEntry(tempreceffi_novertexer,
                                       "Without Vertexer", "PL");
        int runs = 1;
        vector<TString> testnames = {
            "All on(Default)", "SetExtraCleanup(kFALSE)",
            "SetPreselectDedx(kFALSE)", "All Off"};
        for (auto const& hTemp_vertex : vertexertests) {
            hTemp_vertex->SetFillColor(vcolors[runs-1]);
            hTemp_vertex->SetLineColor(vcolors[runs-1]);
            hTemp_vertex->SetMarkerColor(vcolors[runs-1]);
            hTemp_vertex->SetMarkerStyle(Markerset[runs-1]);
            hTemp_vertex->SetMarkerSize(Markerset_size[runs-1]);
            legend_effi_vertexer->AddEntry(hTemp_vertex,
                                       Form("Vertexer with %s", testnames[runs-1].Data()), "PL");
            runs++;
        }
        legend_effi_vertexer->Draw();

        ratio_vertexer->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/MCEffi_vertexer.pdf",
                 workdirectory.Data()));

        /// RAW
        TFile* novertexer = new TFile(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "NoVertexer1.root",
                 path.Data()));

        TFile* vertexer_test1 = new TFile(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Allon1.root",
                 path.Data()));
        TFile* vertexer_test2 = new TFile(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Extraoff1.root",
                 path.Data()));
        TFile* vertexer_test3 = new TFile(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Preseloff1.root",
                 path.Data()));
        TFile* vertexer_test4 = new TFile(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default_"
                 "Alloff1.root",
                 path.Data()));

        TH1D* hRawnovertexer = (TH1D*)novertexer->Get("hDataRawYield");
        TH1* hNofEventnovertexer = (TH1*)novertexer->Get("hNumberofEvent");
        hRawnovertexer->Scale(1/(hNofEventnovertexer->GetBinContent(9)));

        TH1D* hRaw_Test1 = (TH1D*)vertexer_test1->Get("hDataRawYield");
        TH1* hNofEventTest1 = (TH1*)vertexer_test1->Get("hNumberofEvent");
        hRaw_Test1->Scale(1 / (hNofEventTest1->GetBinContent(9)));

        TH1D* hRaw_Test2 = (TH1D*)vertexer_test2->Get("hDataRawYield");
        TH1* hNofEventTest2 = (TH1*)vertexer_test2->Get("hNumberofEvent");
        hRaw_Test2->Scale(1 / (hNofEventTest2->GetBinContent(9)));

        TH1D* hRaw_Test3 = (TH1D*)vertexer_test3->Get("hDataRawYield");
        TH1* hNofEventTest3 = (TH1*)vertexer_test3->Get("hNumberofEvent");
        hRaw_Test3->Scale(1 / (hNofEventTest3->GetBinContent(9)));

        TH1D* hRaw_Test4 = (TH1D*)vertexer_test4->Get("hDataRawYield");
        TH1* hNofEventTest4 = (TH1*)vertexer_test4->Get("hNumberofEvent");
        hRaw_Test4->Scale(1 / (hNofEventTest4->GetBinContent(9)));

        vector<TH1D*> vertexertests_raw;
        vertexertests_raw.push_back(hRaw_Test1);
        vertexertests_raw.push_back(hRaw_Test2);
        vertexertests_raw.push_back(hRaw_Test3);
        vertexertests_raw.push_back(hRaw_Test4);

        hRawnovertexer->SetFillColor(kBlack);
        hRawnovertexer->SetLineColor(kBlack);
        hRawnovertexer->SetMarkerColor(kBlack);
        hRawnovertexer->SetMarkerStyle(Markerset[0]);
        hRawnovertexer->SetMarkerSize(Markerset_size[0]);

        TCanvas* ratio_vertexer_raw = plotHistsAndRatio(
            vertexertests_raw, hRawnovertexer, "", "p_{T}(GeV/c)",
            "Raw Yield", "RAW FINAL NLOG");

        ratio_vertexer_raw->cd(1);
        auto legend_effi_vertexer_raw = new TLegend(.5, .75, .99, .95);
        // legend_effi_vertexer->SetNColumns(2);
        legend_effi_vertexer_raw->SetBorderSize(0);
        legend_effi_vertexer_raw->SetFillStyle(0);
        legend_effi_vertexer_raw->AddEntry(hRawnovertexer,
                                       "Raw Yield with vAN-novertexer", "PL");
        runs = 1;
        for (auto const& hTemp_vertex : vertexertests_raw) {
            hTemp_vertex->SetFillColor(vcolors[runs-1]);
            hTemp_vertex->SetLineColor(vcolors[runs-1]);
            hTemp_vertex->SetMarkerColor(vcolors[runs-1]);
            hTemp_vertex->SetMarkerStyle(Markerset[runs-1]);
            hTemp_vertex->SetMarkerSize(Markerset_size[runs-1]);
            legend_effi_vertexer_raw->AddEntry(
                hTemp_vertex,
                Form("Vertexer with %s", testnames[runs - 1].Data()), "PL");
            runs++;
        }
        
        legend_effi_vertexer_raw->Draw();


        ratio_vertexer_raw->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/Raw_vertexer.pdf",
                 workdirectory.Data()));

        // Corrected spectra
        TH1D* hCornovertexer = (TH1D*)novertexer->Get("hXiSpectrum");

        TH1D* hCor_Test1 = (TH1D*)vertexer_test1->Get("hXiSpectrum");
        TH1D* hCor_Test2 = (TH1D*)vertexer_test2->Get("hXiSpectrum");
        TH1D* hCor_Test3 = (TH1D*)vertexer_test3->Get("hXiSpectrum");
        TH1D* hCor_Test4 = (TH1D*)vertexer_test4->Get("hXiSpectrum");

        vector<TH1D*> vertexertests_cor;
        vertexertests_cor.push_back(hCor_Test1);
        vertexertests_cor.push_back(hCor_Test2);
        vertexertests_cor.push_back(hCor_Test3);
        vertexertests_cor.push_back(hCor_Test4);

        hCornovertexer->SetFillColor(kBlack);
        hCornovertexer->SetLineColor(kBlack);
        hCornovertexer->SetMarkerColor(kBlack);
        hCornovertexer->SetMarkerStyle(Markerset[0]);
        hCornovertexer->SetMarkerSize(Markerset_size[0]);

        TCanvas* ratio_vertexer_Cor =
            plotHistsAndRatio(vertexertests_cor, hCornovertexer, "",
                              "p_{T}(GeV/c)", "Raw Yield", "RAW NLOG FINAL");

        ratio_vertexer_Cor->cd(1);
        auto legend_effi_vertexer_cor = new TLegend(.5, .75, .99, .95);
        // legend_effi_vertexer->SetNColumns(2);
        legend_effi_vertexer_cor->SetBorderSize(0);
        legend_effi_vertexer_cor->SetFillStyle(0);
        legend_effi_vertexer_cor->AddEntry(
            hCornovertexer, "Corrected Yield with vAN-novertexer", "PL");
        runs = 1;
        for (auto const& hTemp_vertex : vertexertests_cor) {
            hTemp_vertex->SetFillColor(vcolors[runs - 1]);
            hTemp_vertex->SetLineColor(vcolors[runs - 1]);
            hTemp_vertex->SetMarkerColor(vcolors[runs - 1]);
            hTemp_vertex->SetMarkerStyle(Markerset[runs - 1]);
            hTemp_vertex->SetMarkerSize(Markerset_size[runs - 1]);
            legend_effi_vertexer_cor->AddEntry(
                hTemp_vertex,
                Form("Vertexer with %s", testnames[runs - 1].Data()), "PL");
            runs++;
        }

        legend_effi_vertexer_cor->Draw();

        ratio_vertexer_Cor->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/Cor_vertexer.pdf",
                 workdirectory.Data()));
    }
    bool mceffi_fixcut = false;
    if (mceffi_fixcut) {
        auto tempreceffi_old = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default1"
                 "_old.root",
                 workdirectory.Data()));
        auto tempreceffi_fix = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_Default"
                 "1.root",
                 workdirectory.Data()));
        tempreceffi_old->SetFillColor(kBlack);
        tempreceffi_old->SetLineColor(kBlack);
        tempreceffi_old->SetMarkerColor(kBlack);
        tempreceffi_old->SetMarkerStyle(Markerset[0]);
        tempreceffi_old->SetMarkerSize(Markerset_size[0]);

        TCanvas* ratio_fixcut = plotHistsAndRatio(
            tempreceffi_fix, tempreceffi_old, "", "p_{T}(GeV/c)",
            "Acceptance x Efficiency x B.R.", "NLOG Efficiency FINAL");

        ratio_fixcut->cd(1);
        auto legend_effi_fixcut = new TLegend(.18, .05, .8, .15);
        // legend_effi_vertexer->SetNColumns(2);
        legend_effi_fixcut->SetBorderSize(0);
        legend_effi_fixcut->SetFillStyle(0);
        legend_effi_fixcut->AddEntry(tempreceffi_old,
                                       "Old (wrong #Lambda CPA)", "PL");
        tempreceffi_fix->SetFillColor(vcolors[0]);
        tempreceffi_fix->SetLineColor(vcolors[0]);
        tempreceffi_fix->SetMarkerColor(vcolors[0]);
        tempreceffi_fix->SetMarkerStyle(Markerset[0]);
        tempreceffi_fix->SetMarkerSize(Markerset_size[0]);
        legend_effi_fixcut->AddEntry(tempreceffi_fix,
                                   "New (fixed #Lambda CPA)", "PL");
        legend_effi_fixcut->Draw();

        ratio_fixcut->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/MCEffi_fixcut.pdf",
                 workdirectory.Data()));

    }

    bool PreliminaryCopare = true;
    if (PreliminaryCopare) {
        std::cout << "test" <<std::endl;
        auto tempSpectra_old = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_studyoff"
                 ".root",
                 path.Data()));
        auto tempSpectra_fix = (TH1D*)GetMCEfficiencyfromName(
            Form("%sAnalysisResults_Extracted_1_Multi_0.00-100.00_studyon"
                 ".root",
                 path.Data()));
        tempSpectra_old->SetFillColor(kBlack);
        tempSpectra_old->SetLineColor(kBlack);
        tempSpectra_old->SetMarkerColor(kBlack);
        tempSpectra_old->SetMarkerStyle(Markerset[0]);
        tempSpectra_old->SetMarkerSize(Markerset_size[0]);

        TCanvas* ratio_pre = plotHistsAndRatio(
            tempSpectra_fix, tempSpectra_old, "", "p_{T}(GeV/c)",
            "Acceptance x Efficiency x B.R.", "NLOG Efficiency FINAL");

        ratio_pre->cd(1);
        auto legend_fixcut = new TLegend(.18, .05, .8, .15);
        // legend_effi_vertexer->SetNColumns(2);
        legend_fixcut->SetBorderSize(0);
        legend_fixcut->SetFillStyle(0);
        legend_fixcut->AddEntry(tempSpectra_old,
                                       "Preliminary result", "PL");
        tempSpectra_fix->SetFillColor(vcolors[0]);
        tempSpectra_fix->SetLineColor(vcolors[0]);
        tempSpectra_fix->SetMarkerColor(vcolors[0]);
        tempSpectra_fix->SetMarkerStyle(Markerset[0]);
        tempSpectra_fix->SetMarkerSize(Markerset_size[0]);
        legend_fixcut->AddEntry(tempSpectra_fix,
                                   "Updated result (work on progress)", "PL");
        legend_fixcut->Draw();

        ratio_pre->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/Spectra_precompare.pdf",
                 workdirectory.Data()));

    }

    // Vertex Effi
    TCanvas* cVertexEff = new TCanvas("cVertexEff", "cVertexEff", w, h);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cVertexEff->SetTickx();
    cVertexEff->SetTicky();
    cVertexEff->SetTopMargin(0.05);
    cVertexEff->SetLeftMargin(0.10);
    // cSigbkg->SetBottomMargin(0.01);
    cVertexEff->SetRightMargin(0.01);
    cVertexEff->SetFillStyle(0);
    cVertexEff->Draw();

    //vector<double> tempVertexEff = {0.995591, 0.989234, 0.96235, 0.892998, 0.860893};
    TH1D* vertex_eff_Cent =
        new TH1D("vertex_eff_Cent", "vertex_eff_Cent",
                 multibinsloop.size() - 1, &multibinsloop[0]);
    for (auto i = 0; i < multibinsloop.size() - 1; i++) {
        auto temptrigeff = GetVertexEfficiency(multibinsloop[i], multibinsloop[i + 1]);
        vertex_eff_Cent->SetBinContent(i + 1, temptrigeff[0]);
        vertex_eff_Cent->SetBinError(i + 1, temptrigeff[1]);
    }
    vertex_eff_Cent->SetMinimum(0.80);
    vertex_eff_Cent->SetMaximum(1.05);
    vertex_eff_Cent->GetYaxis()->SetTitle(
        "Vertex Loss Correction(N_{evt;trig}/N_{evt;GoodVertex})");
    vertex_eff_Cent->GetXaxis()->SetTitle("Multiplicity Percentile(V0M) [%]");
    cVertexEff->cd();
    vertex_eff_Cent->Draw("HIST text");
    double vtxeffi0100 =
        GetVertexEfficiency(0, 100)[0];
    t->DrawLatex(0.7, 0.87, Form("#bf{MB(0-100): %.3f}", vtxeffi0100));

    cVertexEff->SaveAs(
        Form("%s1_Multi_0.00-100.00_Default1/hMulti_VetexEffi.pdf",
             workdirectory.Data()));

    // Trigger Effi
    TCanvas* cTriggerEff = new TCanvas("cTriggerEff", "cTriggerEff", w, h);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cTriggerEff->SetTickx();
    cTriggerEff->SetTicky();
    cTriggerEff->SetTopMargin(0.05);
    cTriggerEff->SetLeftMargin(0.10);
    // cSigbkg->SetBottomMargin(0.01);
    cTriggerEff->SetRightMargin(0.01);
    cTriggerEff->SetFillStyle(0);
    cTriggerEff->Draw();

    vector<double> temptrigeff;
    TH1D* trigger_eff_Cent =
        new TH1D("trigger_eff_Cent", "trigger_eff_Cent",
                 multibinsloop.size() - 1, &multibinsloop[0]);
    for (auto i = 0; i < multibinsloop.size() - 1; i++) {
        temptrigeff = GetTrigEfficiency(multibinsloop[i], multibinsloop[i + 1]);
        trigger_eff_Cent->SetBinContent(i + 1, temptrigeff[0]);
        trigger_eff_Cent->SetBinError(i + 1, temptrigeff[1]);
    }
    trigger_eff_Cent->SetMinimum(0.6);
    trigger_eff_Cent->SetMaximum(1.2);
    trigger_eff_Cent->GetYaxis()->SetTitle(
        "Trigger Efficiency(INT7|INEL>0 / INEL>0)");
    trigger_eff_Cent->GetXaxis()->SetTitle("Multiplicity Percentile(V0M) [%]");
    temptrigeff = GetTrigEfficiency(0, 100);
    cTriggerEff->cd();
    trigger_eff_Cent->Draw("HIST text");
    t->DrawLatex(0.7, 0.87, Form("#bf{MB(0-100): %.3f}", temptrigeff[0]));

    cTriggerEff->SaveAs(
        Form("%s1_Multi_0.00-100.00_Default1/hMulti_TrigEffi.pdf",
             workdirectory.Data()));
    // Signal losss

    auto hDefSigloss = GetSigLoss(0, 100);
    vector<TH1D*> buffer_sigloss;
    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        auto temp_multi =
            GetSigLoss(multibinsloop[imultibin], multibinsloop[imultibin + 1]);
        buffer_sigloss.push_back(temp_multi);
    }
    TCanvas* ratio_sigloss_mult =
        plotHistsAndRatio(buffer_sigloss, hDefSigloss, "", "p_{T}(GeV/c)",
                          "Signal loss", "NLOG sigloss FINAL CORR");

    ratio_sigloss_mult->cd(1);
    auto legend_sigloss_ratio = new TLegend(.18, .01, .5, .20);
    legend_sigloss_ratio->SetNColumns(2);
    legend_sigloss_ratio->SetBorderSize(0);
    legend_sigloss_ratio->SetFillStyle(0);
    legend_sigloss_ratio->AddEntry(hDef, "MB(0-100)", "PL");

    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        buffer_sigloss.at(imultibin)->SetFillColor(vcolors[imultibin]);
        buffer_sigloss.at(imultibin)->SetLineColor(vcolors[imultibin]);
        buffer_sigloss.at(imultibin)->SetMarkerColor(vcolors[imultibin]);
        buffer_sigloss.at(imultibin)->SetMarkerStyle(Markerset[imultibin]);
        buffer_sigloss.at(imultibin)->SetMarkerSize(Markerset_size[imultibin]);
    }

    for (int imultibin = 0; imultibin < multibinsloop.size() - 1; imultibin++) {
        legend_sigloss_ratio->AddEntry(
            buffer_sigloss.at(imultibin),
            Form("%.0f - %.0f", multibinsloop.at(imultibin),
                 multibinsloop.at(imultibin + 1)),
            "PL");
    }
    legend_sigloss_ratio->Draw();

    ratio_sigloss_mult->SaveAs(
        Form("%s1_Multi_0.00-100.00_Default1/Sigloss_Multi.pdf",
             workdirectory.Data()));

    // Get SysError
    TCanvas* cSysError = new TCanvas("cSysError", "cSysError", w, h);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cSysError->SetTickx();
    cSysError->SetTicky();
    cSysError->SetTopMargin(0.05);
    cSysError->SetLeftMargin(0.10);
    // cSigbkg->SetBottomMargin(0.01);
    cSysError->SetRightMargin(0.01);
    cSysError->SetFillStyle(0);
    cSysError->Draw();

    auto legend_sysError = new TLegend(.6, .6, .9, .9);
    // legend_sysError->SetNColumns(2);
    legend_sysError->SetBorderSize(0);
    legend_sysError->SetFillStyle(0);


    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    Int_t FIh = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 multibincheck.size());
    vcolors = {};
    for (int color = 0; color < multibincheck.size(); color++) {
        // histColors.push_back(FIh+color);
        // vcolors.push_back(FIh+color);
        vcolors.insert(vcolors.begin(), FIh + color);
    }

    vector<TH1D*> SysErrorHists;
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        SysErrorHists.push_back(GetSysError(multibincheck[imultibin][0],
                                            multibincheck[imultibin][1]));
        SysErrorHists.at(imultibin)->SetMaximum(0.4);
        SysErrorHists.at(imultibin)->SetFillColor(0);
        SysErrorHists.at(imultibin)->SetLineColor(vcolors[imultibin]);
        SysErrorHists.at(imultibin)->SetMarkerColor(vcolors[imultibin]);
        SysErrorHists.at(imultibin)->SetMarkerStyle(Markerset[imultibin]);
        SysErrorHists.at(imultibin)->SetMarkerSize(Markerset_size[imultibin]);

        legend_sysError->AddEntry(
            SysErrorHists.at(imultibin),
            Form("systematic Error (%.2f-%.2f)", multibincheck[imultibin][0],
                 multibincheck[imultibin][1]),
            "PL");
    }

    cSysError->cd();
    SysErrorHists[0]->SetMaximum(0.5);
    SysErrorHists[0]->Draw("HIST");
    for (int imultibin = 1; imultibin < multibinsloop.size(); imultibin++) {
        SysErrorHists[imultibin]->Draw("HIST same");
    }
    legend_sysError->Draw();
    cSysError->SaveAs(Form("%s1_Multi_0.00-100.00_Default1/SysErrors_mult.pdf",
                           workdirectory.Data()));

}
vector<double> GetTrigEfficiency(double multi_start, double multi_end) {
    vector<double> multibin = {multi_start, multi_end};

    TString inputOptions = "Default";
    int cutbin = 1;
    TString finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
             workdirectory.Data(), cutbin, multibin[0], multibin[1],
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* htrigeffi = (TH1D*)inputfile->Get("hTriggerEffi");
    vector<double> trigeffi;
    trigeffi.push_back(htrigeffi->GetBinContent(1));
    trigeffi.push_back(htrigeffi->GetBinError(1));
    inputfile->Close();
    return trigeffi;
}
vector<double> GetVertexEfficiency(double multi_start, double multi_end) {
    vector<double> multibin = {multi_start, multi_end};

    TString inputOptions = "Default";
    int cutbin = 1;
    TString finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
             workdirectory.Data(), cutbin, multibin[0], multibin[1],
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* htrigeffi = (TH1D*)inputfile->Get("hVertexEffi");
    vector<double> vtxeffi;
    vtxeffi.push_back(htrigeffi->GetBinContent(1));
    vtxeffi.push_back(htrigeffi->GetBinError(1));
    inputfile->Close();
    return vtxeffi;
}
TH1D* GetSigLoss(double multi_start, double multi_end) {
    vector<double> multibin = {multi_start, multi_end};

    TString inputOptions = "Default";
    int cutbin = 1;
    TString finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
             workdirectory.Data(), cutbin, multibin[0], multibin[1],
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* base = (TH1D*)inputfile->Get("hMCSigLoss");
    gROOT->cd();
    TH1D* hReturn = (TH1D*)base->Clone();
    inputfile->Close();
    return hReturn;
}
TH1D* GetMCEfficiency(double multi_start, double multi_end) {
    vector<double> multibin = {multi_start, multi_end};

    TString inputOptions = "MCcheck";
    int cutbin = 1;
    TString finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
             workdirectory.Data(), cutbin, multibin[0], multibin[1],
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* base = (TH1D*)inputfile->Get("hMCReconEffi");
    gROOT->cd();
    TH1D* hReturn = (TH1D*)base->Clone();
    inputfile->Close();
    return hReturn;
}
TH1D* GetMCEfficiencyfromName(TString inputfilename) {
    TFile* inputfile = new TFile(inputfilename.Data());
    TH1D* base = (TH1D*)inputfile->Get("hMCReconEffi");
    cout << "name: " << inputfilename.Data() << endl;
    gROOT->cd();
    TH1D* hReturn = (TH1D*)base->Clone();
    inputfile->Close();
    return hReturn;
}

TH1D* GetSpectra(double multi_start, double multi_end) {
    vector<double> multibin = {multi_start, multi_end};

    TString inputOptions = "Default";
    int cutbin = 1;
    TString finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
             workdirectory.Data(), cutbin, multibin[0], multibin[1],
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* base = (TH1D*)inputfile->Get("hXiSpectrum");
    gROOT->cd();
    TH1D* hReturn = (TH1D*)base->Clone();
    inputfile->Close();
    return hReturn;
}
TH1D* GetSpectrastat(double multi_start, double multi_end){
    TString bininfo = GetBinInfo(multi_start,multi_end);
    finalfile = Form("./AnalysisResults_Xi1530_systematic%s.root", bininfo.Data());
    TFile* inputfile = new TFile(finalfile.Data());
    TString text;
    if(isINEL)
        text = "hSpectra_INEL_stat";
    else
        text = Form("hSpectra_%.2f_%.2f_stat",multi_start,multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(text);
    gROOT->cd();
    TH1D* hReturn = (TH1D*)hr->Clone();
    inputfile->Close();
    return hReturn;
}
TH1D* GetSysError(double multi_start, double multi_end) {
    TH1D* hr = (TH1D*)GetSpectrasys(multi_start,multi_end);
    vector<double> systerror;
    systerror.push_back(0);
    for (int i = 0; i < hr->GetNbinsX(); i++) {
        systerror.push_back(hr->GetBinError(i + 1) / hr->GetBinContent(i + 1));
    }
    auto htempsys =
        (TH1D*)MakeHistfromArray(Form("SysError_%.2f-%.2f", multi_start, multi_end),
                          systerror, ptbin);
    return htempsys;
}

TF1* GetFittedFunction(double multi_start, double multi_end, int fitfunction) {
    TFile* inputfile = new TFile(Form("%sXi1530YieldMean.root", path.Data()));
    TCanvas* temp = (TCanvas*)inputfile->Get(
        Form("c_h_%.2f-%.2f_stat_tot_%s_fit_0.80_8.80", multi_start, multi_end,
             fitfuctions[fitfunction].Data()));
    TF1* Ffunc =
        (TF1*)temp->GetPrimitive(Form("%s", fitfuctions[fitfunction].Data()));
    gROOT->cd();
    TF1* Fr = (TF1*)Ffunc->Clone();
    inputfile->Close();
    return Fr;
}
TH1* GetExtraSpectra(double multi_start, double multi_end, int histooption) {
    cout << "intput: " << multi_start << "-" << multi_end
         << ", option: " << exhistograms[histooption] << endl;
    TFile* inputfile = new TFile(Form("%sXi1530YieldMean.root", path.Data()));
    TCanvas* temp = (TCanvas*)inputfile->Get(
        Form("c_h_%.2f-%.2f_sys_%s_fLevi_fit_0.80_8.80", multi_start, multi_end,
             exhistograms[histooption].Data()));
    TH1* hEx = (TH1*)temp->GetPrimitive(Form("h_%.2f-%.2f_sys_%s", multi_start,
                                             multi_end,
                                             exhistograms[histooption].Data()));
    gROOT->cd();
    TH1* hr = (TH1*)hEx->Clone();
    inputfile->Close();
    return hr;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TCanvas* plotHistsAndRatio(TH1D* numeratorHist,
                           TH1D* denominatorHist,
                           TString title = "",
                           TString xTitle = "",
                           TString yTitle = "",
                           TString options = "") {
    vector<TH1D*> numeratorHistograms;
    numeratorHistograms.push_back(numeratorHist);
    TCanvas* c1 = plotHistsAndRatio(numeratorHistograms, denominatorHist, title,
                                    xTitle, yTitle, options);
    return c1;
}
TCanvas* plotHistsAndRatio(vector<TH1D*> numeratorHistograms,
                           TH1D* denominatorHist,
                           TString title = "",
                           TString xTitle = "",
                           TString yTitle = "",
                           TString options = "") {
    int numberOfNumeratorHists = numeratorHistograms.size();
    if (numberOfNumeratorHists > 3) {
        cout << "Too many histograms for numerator (currently only supports up "
                "to 3)"
             << endl;
        exit;
    }
    if (!denominatorHist) {
        cout << "denominatorHist provided does not exist" << endl;
        exit;
    }

    //*************************************************
    // Variables
    bool topPlotLogY = 1;  // 0 = no log; 1= log
    if (options.Contains("sigloss"))
        topPlotLogY = 0;
    bool bottomPlotLogY = 1;                 // 0 = no log; 1= log
    TString yTitle2 = "Ratio to MB(0-100)";  // bottom plot y axis title
    if (options.Contains("SPECTRA"))
        yTitle2 = "Ratio to INEL>0";
    if (options.Contains("DIFF"))
        yTitle2 = "Ratio to 7 TeV Results";
    if (options.Contains("HARDSOFT"))
        yTitle2 = "Ratio to Default";
    if (options.Contains("RAW")){
        topPlotLogY = 0;
        yTitle2 = "Ratio to novertexer";
    }
    if (options.Contains("NLOG"))
        bottomPlotLogY = 0;

    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    Int_t FIh = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 numberOfNumeratorHists);

    vector<int> histColors = {};
    vcolors = {};
    for (int color = 0; color < numberOfNumeratorHists; color++) {
        // histColors.push_back(FIh+color);
        // vcolors.push_back(FIh+color);
        histColors.insert(histColors.begin(), FIh + color);
        vcolors.insert(vcolors.begin(), FIh + color);
    }

    int histDenominatorColor = kBlue;

    float defaultRatioYmin = 4e-2;
    float defaultRatioYmax = 2e1;
    if (options.Contains("Efficiency")) {
        defaultRatioYmin = 0.7;
        defaultRatioYmax = 1.4;
    }
    if (options.Contains("RAW")) {
        defaultRatioYmin = 0.4;
        defaultRatioYmax = 1.6;
    }
    if (options.Contains("SPECTRA")) {
        defaultRatioYmin = 0.08;
        defaultRatioYmax = 8;
    }
    if (options.Contains("7T")) {
        defaultRatioYmin = 0.9;
        defaultRatioYmax = 3.0;
    }
    if (options.Contains("CUTS")) {
        defaultRatioYmin = 0.8;
        defaultRatioYmax = 1.2;
    }
    if (options.Contains("sigloss")) {
        defaultRatioYmin = 0.9;
        defaultRatioYmax = 1.25;
    }
    if (options.Contains("fine")) {
        defaultRatioYmin = 0.9;
        defaultRatioYmax = 1.1;
    }
    if (options.Contains("HARDSOFT")) {
        defaultRatioYmin = 0.5;
        defaultRatioYmax = 1.5;
    }

    // END of Variables
    //*************************************************
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TString name = yTitle + options;
    TCanvas* c1 = new TCanvas(name.Data(), name.Data(), 0, 0, 850, 1150);
    c1->Range(0, 0, 1, 1);

    vector<TH1D*> hists;
    for (int i = 0; i < numberOfNumeratorHists; i++) {
        hists.push_back((TH1D*)numeratorHistograms[i]->Clone());
    }
    TH1D* denominatorHistogram = (TH1D*)denominatorHist->Clone();

    // Create ratio histograms
    vector<TH1D*> hist_over_denomHist;
    for (int i = 0; i < numberOfNumeratorHists; i++) {
        hist_over_denomHist.push_back((TH1D*)numeratorHistograms[i]->Clone());
        hist_over_denomHist[i]->GetTitle();
        if (options.Contains("CORR")){
            hist_over_denomHist[i]->Divide(hist_over_denomHist[i],denominatorHistogram,1.0,1.0,"B");
        }
        else
            hist_over_denomHist[i]->Divide(denominatorHistogram);
    }

    //*************************************************
    // Bottom plot
    TPad* c1_1 = new TPad("c1_1", "newpad", 0.01, 0.01, 1, 0.32);
    c1_1->Draw();
    c1_1->cd();
    c1_1->SetTickx(1);
    c1_1->SetTicky(1);
    c1_1->SetTopMargin(0.01);
    c1_1->SetLeftMargin(0.15);
    c1_1->SetBottomMargin(0.2);
    c1_1->SetRightMargin(0.01);
    c1_1->SetFillStyle(0);
    c1_1->SetLogy(bottomPlotLogY);

    if (options.Contains("SYS"))
        hist_over_denomHist[0]->Draw("E2");
    else
        hist_over_denomHist[0]->Draw("E");
    hist_over_denomHist[0]->SetLineWidth(1);
    hist_over_denomHist[0]->SetLineColor(histColors[0]);
    hist_over_denomHist[0]->SetMarkerColor(histColors[0]);
    hist_over_denomHist[0]->SetMinimum(defaultRatioYmin);
    hist_over_denomHist[0]->SetMaximum(defaultRatioYmax);
    hist_over_denomHist[0]->GetYaxis()->SetNdivisions(5);
    hist_over_denomHist[0]->SetTitle(";" + xTitle + ";" + yTitle2);
    hist_over_denomHist[0]->GetXaxis()->SetTitleSize(0.08);
    hist_over_denomHist[0]->GetXaxis()->SetLabelSize(0.08);
    hist_over_denomHist[0]->GetYaxis()->SetLabelSize(0.08);
    hist_over_denomHist[0]->GetYaxis()->SetTitleSize(0.08);
    hist_over_denomHist[0]->GetYaxis()->SetTitleOffset(0.6);
    if (options.Contains("FINAL"))
        hist_over_denomHist[0]->GetXaxis()->SetRangeUser(0.8, 7.4);

    for (int i = 1; i < numberOfNumeratorHists; i++) {
        hist_over_denomHist[i]->SetLineWidth(1);
        hist_over_denomHist[i]->SetLineColor(histColors[i]);
        hist_over_denomHist[i]->SetFillColor(histColors[i]);
        hist_over_denomHist[i]->SetMarkerColor(histColors[i]);
        hist_over_denomHist[i]->SetMarkerStyle(Markerset[i]);
        hist_over_denomHist[i]->SetMarkerSize(Markerset_size[i]);
        if (options.Contains("SYS"))
            hist_over_denomHist[i]->Draw("E2 same");
        else
            hist_over_denomHist[i]->Draw("E same");
    }
    TF1* line1 = new TF1("line1", "1", -1, 15);
    line1->SetLineColor(1);
    line1->SetLineWidth(1);
    line1->SetLineStyle(2);
    line1->Draw("same");
    // End bottom plot
    //*************************************************

    //*************************************************
    // Top Plot
    c1->cd();
    TPad* c1_2 = new TPad("c1_2", "newpad", 0.01, 0.32, 1, 1);
    c1_2->SetLogy(true);
    c1_2->SetTickx(1);
    c1_2->SetTicky(1);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetTopMargin(0.01);
    c1_2->SetLeftMargin(0.15);
    c1_2->SetBottomMargin(0.01);
    c1_2->SetRightMargin(0.01);
    c1_1->SetFillStyle(0);

    denominatorHistogram->SetLineWidth(1);
    if (options.Contains("SPECTRA"))
        histDenominatorColor = 0;
    if (options.Contains("HARDSOFT"))
        histDenominatorColor = 1;
    //denominatorHistogram->SetLineColor(histDenominatorColor);
    //denominatorHistogram->SetMarkerColor(histDenominatorColor);
    denominatorHistogram->SetLineColor(kBlack);
    denominatorHistogram->SetMarkerColor(kBlack);
    if (options.Contains("SYS"))
        denominatorHistogram->Draw("E2");
    else
        denominatorHistogram->Draw("E");
    denominatorHistogram->SetLabelSize(0.0);
    denominatorHistogram->GetXaxis()->SetTitleSize(0.00);
    denominatorHistogram->GetYaxis()->SetLabelSize(0.04);
    denominatorHistogram->GetYaxis()->SetTitleSize(0.05);
    denominatorHistogram->GetYaxis()->SetTitleOffset(1.2);
    denominatorHistogram->SetTitle(title + ";;" + yTitle);
    if (options.Contains("FINAL"))
        denominatorHistogram->GetXaxis()->SetRangeUser(0.8, 7.4);
    denominatorHistogram->SetMaximum(2e-1);
    if (options.Contains("SPECTRA"))
        denominatorHistogram->SetMaximum(5e-1);
    if (options.Contains("Efficiency"))
        denominatorHistogram->SetMaximum(5e-1);
    if (options.Contains("7T"))
        denominatorHistogram->SetMaximum(5e-3);
    if (options.Contains("LOW"))
        denominatorHistogram->SetMinimum(1e-6);
    if (options.Contains("FINAL"))
        denominatorHistogram->SetMinimum(5e-8);
    if (options.Contains("Efficiency"))
        denominatorHistogram->SetMinimum(5e-4);

    if (options.Contains("CUTS")) {
        denominatorHistogram->SetMaximum(5e-3);
    }
    if (options.Contains("sigloss")) {
        denominatorHistogram->SetMaximum(1.2);
        denominatorHistogram->SetMinimum(0.9);
    }

    if (options.Contains("RAW")) {
        denominatorHistogram->SetMaximum(hists[0]->GetMaximum()*1.2);
        denominatorHistogram->SetMinimum(hists[0]->GetMinimum()*0.7);
    }
    for (int i = 0; i < numberOfNumeratorHists; i++) {
        hists[i]->SetLineWidth(1);
        if (options.Contains("POW"))
            hists[i]->Scale(pow(2, numberOfNumeratorHists - i));
        hists[i]->SetLineColor(histColors[i]);
        hists[i]->SetFillColor(histColors[i]);
        hists[i]->SetMarkerColor(histColors[i]);
        hists[i]->SetMarkerStyle(Markerset[i]);
        hists[i]->SetMarkerSize(Markerset_size[i]);
        if (options.Contains("SYS"))
            hists[i]->Draw("2E same");
        else
            hists[i]->Draw("E same");
    }

    c1_2->SetLogy(topPlotLogY);
    // End bottom plot
    //*************************************************

    return c1;
}
double myLevyPtFunc(Double_t* x, Double_t* par) {
    Double_t lMass = 0;
    lMass = 1.5318;  // Xi* mass

    Double_t ldNdy = par[0];          // dN/dy
    Double_t l2pi = 2 * TMath::Pi();  // 2pi (cancels in return statement)
    Double_t lTemp = par[1];          // Temperature
    Double_t lPower = par[2];         // power=n

    Double_t lBigCoef =
        ((lPower - 1) * (lPower - 2)) /
        (l2pi * lPower * lTemp * (lPower * lTemp + lMass * (lPower - 2)));
    Double_t lInPower = 1 + (TMath::Sqrt(x[0] * x[0] + lMass * lMass) - lMass) /
                                (lPower * lTemp);

    return l2pi * ldNdy * x[0] * lBigCoef *
           TMath::Power(lInPower, (-1) * lPower);
}
TH1D* GetSpectrasys(double multi_start, double multi_end){
    TString bininfo = GetBinInfo(multi_start,multi_end);
    finalfile = Form("./AnalysisResults_Xi1530_systematic%s.root", bininfo.Data());
    TFile* inputfile = new TFile(finalfile.Data());
    TString text;
    if(isINEL)
        text = "hSpectra_INEL_sys";
    else
        text = Form("hSpectra_%.2f_%.2f_sys", multi_start, multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(text);
    gROOT->cd();
    TH1D* hreturn = (TH1D*)hr->Clone();
    inputfile->Close();

    return hreturn;
}
TString GetBinInfo(double multis, double multie){
    TString bininfo = "_";
    if ((multis == 0) && (multie == 0)){ 
        //INEL case
        bininfo += "INEL";
        isINEL = true;
    }
    else if (multie < 0.5) {
        bininfo += Form("%.2f", multis);
        bininfo += "-";
        bininfo += Form("%.2f", multie);
    }
    else {
        bininfo += Form("%.0f", multis);
        bininfo += "-";
        bininfo += Form("%.0f", multie);   
    }
    return bininfo;
}