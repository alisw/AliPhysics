#include "DrawingHelper.C"

// Public variables
int canvasCount = 1;
vector<int> MaterialColors
    = {kOrange+10, kOrange-3, kTeal+2, kAzure+2, kViolet+4};
//TString PhysicsFile = "AnalysisResults_Xi1530_PhysicsResult_0.00-0.00_0-100_0-10_10-30_30-50_50-70_70-100_bins.root";
TString PhysicsFile = "AnalysisResults_Xi1530_PhysicsResult_0-10_10-30_30-50_50-70_70-100_bins.root";
//General Text
TLatex* t0 = new TLatex();
TLatex* t0R = new TLatex();
TLatex* t = new TLatex();
TLatex* tR = new TLatex();
TLatex* t2 = new TLatex();
TLatex* t2R = new TLatex();
TLatex* t3 = new TLatex();
TLatex* t3R = new TLatex();
TLatex* t4 = new TLatex();
TLatex* t4R = new TLatex();

// Common macro
TCanvas* GetCanvas(TString name = "Canvas", double w = 960, double h = 720);

// Draw functions
void DrawSignalBackground();
void DrawFittedSignal();
void DrawRecEffi();
void DrawSpectraMB();
void DrawSpectraMBRatio();
void DrawSpectraMulti();
void DrawSpectraMultiPi();
void DrawdNdy();
void DrawMeanpT();
void DrawRatioToPi();
void DrawRatioToXi();

void PlotXi1530ApprovalFigure(){
    // Common latex
    // for memo, super small
    t0->SetNDC();
    t0->SetTextSize(0.030);
    t0R->SetNDC();
    t0R->SetTextSize(0.030);
    t0R->SetTextAlign(33);
    // for memo, small
    t->SetNDC();
    t->SetTextSize(0.035);
    tR->SetNDC();
    tR->SetTextSize(0.035);
    tR->SetTextAlign(33);
    // for warning, big
    t2->SetNDC();
    t2->SetTextSize(0.04);
    t2R->SetNDC();
    t2R->SetTextSize(0.04);
    t2R->SetTextAlign(33);
    // for small pad, huge
    t3->SetNDC();
    t3->SetTextSize(0.06);
    t3R->SetNDC();
    t3R->SetTextSize(0.06);
    t3R->SetTextAlign(33);
    // 
    t4->SetNDC();
    t4->SetTextSize(0.08);
    t4R->SetNDC();
    t4R->SetTextSize(0.08);
    t4R->SetTextAlign(33);

    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);

    //DrawRatioToXi();
    //gSystem->Exit(1);
    
    //DrawSpectraMultiPi
    DrawSpectraMultiPi();
    
    // 1. Signal to Background figure
    DrawSignalBackground();

    // 2. Fitting figure
    DrawFittedSignal();
    //gSystem->Exit(1);
    // 3. Rec.Effi
    DrawRecEffi();
    
    // 4. pT Spectra INEL
    DrawSpectraMB();

    // 4.1. Ratio to 7 TeV Result
    DrawSpectraMBRatio();

    // 5. pT Spectra with multiplicity percentile
    DrawSpectraMulti();

    // 6. pT-integrated yield dN/dy
    DrawdNdy();
    
    // 7. mean pT
    DrawMeanpT();
    
    // 8. Ratio to Pi
    DrawRatioToPi();

    // 9. Ratio to Xi
    DrawRatioToXi();
    
        
    //DrawSpectraMBRatio();
    
}
void DrawSignalBackground(){
    vector<double> fDrawRange = {1.48, 1.79};


    TCanvas* cSigbkg = GetCanvas("cSigbkg");
    cSigbkg->SetLeftMargin(0.10);
    cSigbkg->Draw();
    TFile* fSigbkg = new TFile(
        "data/AnalysisResults_Extracted_1_Multi_10.00-30.00_Default1.root");

    auto hsignal = (TH1*)fSigbkg->Get("hSignalOnly_3");   // signal
    auto hbkg = (TH1*)fSigbkg->Get("hBkgOnly_3");         // Bkg
    auto hbkgnorm = (TH1*)fSigbkg->Get("hBkgNorm_3");     // Bkg norm

    hsignal->SetMarkerStyle(20);
    hsignal->SetMarkerSize(1);
    hsignal->SetMarkerColor(kBlack);
    hsignal->SetLineColor(kBlack);
    hsignal->SetMaximum(hsignal->GetMaximum()*1.1);
    hsignal->GetYaxis()->SetTitleOffset(1.2);
    hsignal->GetXaxis()->SetTitleOffset(1.2);
    hsignal->GetYaxis()->SetTitleSize(0.04);
    hsignal->GetXaxis()->SetTitleSize(0.04);
    hsignal->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);
    hsignal->GetYaxis()->SetTitle("Counts / (4.0 MeV/#it{c} ^{2})");
    hsignal->GetXaxis()->SetTitle("#it{M}_{#Xi#pi} (GeV/#it{c}^{2})");

    hbkg->SetMarkerStyle(4);
    hbkg->SetMarkerColor(MaterialColors[0]);
    hbkg->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);
    hbkgnorm->SetFillColor(MaterialColors[0]);
    hbkgnorm->SetFillStyle(3005);
    hbkgnorm->SetLineColor(0);
    hbkgnorm->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);

    cSigbkg->cd();
    hsignal->Draw("PZ");
    hbkgnorm->Draw("BAR same");
    hbkg->Draw("same");

    t2R->DrawLatex(0.72, 0.91, "ALICE Preliminary");
    t2R->DrawLatex(0.93, 0.91, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    t2R->DrawLatex(0.93, 0.85, "#bf{V0M Multiplicity Event Classes 10 - 30% }");
    t2R->DrawLatex(0.93, 0.79,
                  "#bf{1.6 < #it{p}_{T} < 2.0 GeV/#it{c}, |#it{y}| < 0.5}");

    t2R->DrawLatex(0.93, 0.73,
                  "#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} #rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}");
    //auto lSigBkg = new TLegend(0.13, 0.74, 0.45, 0.89);
    auto lSigBkg = new TLegend(0.56, 0.50, 0.93, 0.66);
    lSigBkg->SetFillStyle(0);
    lSigBkg->AddEntry(hsignal, "Data (stat. uncert.)", "PLE");
    lSigBkg->AddEntry(hbkg, "Mixed-event background", "PLE");
    lSigBkg->AddEntry(hbkgnorm, "Normalisation Region", "F");
    lSigBkg->Draw();

    SaveCanvas(cSigbkg,"figure_sigbkg_preliminary","figs/Approval/");
    // 1.1 0-100 case
    TFile* fSigbkg0100 = new TFile(
        "data/AnalysisResults_Extracted_1_Multi_0.00-100.00_Default1.root");

    auto hsignal0100 = (TH1*)fSigbkg0100->Get("hSignalOnly_3");  // signal
    auto hbkg0100 = (TH1*)fSigbkg0100->Get("hBkgOnly_3");        // Bkg
    auto hbkgnorm0100 = (TH1*)fSigbkg0100->Get("hBkgNorm_3");    // Bkg norm

    hsignal0100->SetMarkerStyle(20);
    hsignal0100->SetMarkerSize(1);
    hsignal0100->SetMarkerColor(kBlack);
    hsignal0100->SetLineColor(kBlack);
    hsignal0100->SetMaximum(hsignal0100->GetMaximum()*1.1);
    hsignal0100->GetYaxis()->SetTitleOffset(1.2);
    hsignal0100->GetXaxis()->SetTitleOffset(1.2);
    hsignal0100->GetYaxis()->SetTitleSize(0.04);
    hsignal0100->GetXaxis()->SetTitleSize(0.04);
    hsignal0100->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);
    hsignal0100->GetYaxis()->SetTitle("Counts / (4.0 MeV/#it{c} ^{2})");
    hsignal0100->GetXaxis()->SetTitle("#it{M}_{#Xi#pi} (GeV/#it{c}^{2})");

    hbkg0100->SetMarkerStyle(4);
    hbkg0100->SetMarkerColor(MaterialColors[0]);
    hbkg0100->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);
    hbkgnorm0100->SetFillColor(MaterialColors[0]);
    hbkgnorm0100->SetFillStyle(3005);
    hbkgnorm0100->SetLineColor(0);
    hbkgnorm0100->SetLineStyle(0);
    //hbkgnorm0100->SetMarkerSyle(0);
    hbkgnorm0100->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);

    cSigbkg->cd();
    hsignal0100->Draw("PZ");
    hbkgnorm0100->Draw("BAR same");
    hbkg0100->Draw("same");

    t2R->DrawLatex(0.62, 0.91, "ALICE Preliminary");
    t2R->DrawLatex(0.93, 0.91, "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV}");
    //t2R->DrawLatex(0.93, 0.85, "#bf{V0M Multiplicity Event Classes 0 - 100% }");
    t2R->DrawLatex(0.93, 0.85,
                  "#bf{1.6 < #it{p}_{T} < 2.0 GeV/#it{c}, |#it{y}| < 0.5}");

    t2R->DrawLatex(0.93, 0.78,
                  "#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} #rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}");
    //t3R->DrawLatex(0.93, 0.72,"#bf{#Xi(1530)^{0} + cc #rightarrow #Xi^{#mp} + #pi^{#pm}}");

    auto lSigBkg0100 = new TLegend(0.56, 0.52, 0.93, 0.68);
    lSigBkg0100->SetFillStyle(0);
    lSigBkg0100->AddEntry(hsignal0100, "Data (stat. uncert.)", "PLE");
    lSigBkg0100->AddEntry(hbkg0100, "Mixed-event background", "PLE");
    lSigBkg0100->AddEntry(hbkgnorm0100, "Normalisation Region", "F");
    lSigBkg0100->Draw();

    SaveCanvas(cSigbkg, "figure_sigbkgMB_preliminary", "figs/Approval/");
}
void DrawFittedSignal(){
    vector<double> fDrawRange = {1.484, 1.76};

    TCanvas* cSignalFit = GetCanvas("cSignalFit");
    cSignalFit->SetLeftMargin(0.10);
    cSignalFit->Draw();
    TFile* fSignalFit = new TFile(
        "data/AnalysisResults_Extracted_1_Multi_10.00-30.00_Default1.root");

    auto hfitsignal = (TH1*)fSignalFit->Get("hSignalBkgSubtraction_3");  // signal
    auto fFitSum = (TF1*)fSignalFit->Get("fDataFitResult_3");         // Fit
    auto fFitPeak = (TF1*)fSignalFit->Get("fDataFitOnlyResult_3");         // Fit
    auto fFitbkg = (TF1*)fSignalFit->Get("fDataFitBkgResult_3");      // Bkg

    hfitsignal->SetMarkerStyle(20);
    hfitsignal->SetMarkerSize(1);
    hfitsignal->SetMarkerColor(kBlack);
    hfitsignal->SetLineColor(kBlack);
    hfitsignal->SetMaximum(hfitsignal->GetMaximum()*1.3);
    hfitsignal->GetYaxis()->SetTitleOffset(1.2);
    hfitsignal->GetXaxis()->SetTitleOffset(1.2);
    hfitsignal->GetYaxis()->SetTitleSize(0.04);
    hfitsignal->GetXaxis()->SetTitleSize(0.04);
    hfitsignal->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);
    hfitsignal->GetYaxis()->SetTitle("Counts / (4.0 MeV/#it{c} ^{2})");
    hfitsignal->GetXaxis()->SetTitle("#it{M}_{#Xi#pi} (GeV/#it{c}^{2})");

    fFitbkg->SetLineColor(kBlack);
    fFitbkg->SetLineWidth(2);
    fFitbkg->SetRange(1.484, 1.6);
    fFitSum->SetLineColor(MaterialColors[0]);
    fFitSum->SetLineWidth(2);
    fFitSum->SetNpx(500);

    cSignalFit->cd();
    hfitsignal->Draw("PZ");
    fFitbkg->Draw("same");
    fFitSum->Draw("same");
    //fFitPeak->Draw("same");

    t2R->DrawLatex(0.72, 0.91, "ALICE Preliminary");
    t2R->DrawLatex(0.93, 0.91, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    t2R->DrawLatex(0.93, 0.85, "#bf{V0M Multiplicity Event Classes 10 - 30% }");
    t2R->DrawLatex(0.93, 0.79,
                  "#bf{1.6 < #it{p}_{T} < 2.0 GeV/#it{c}, |#it{y}| < 0.5}");
    t2R->DrawLatex(0.93, 0.70,
                  "#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} #rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}");
    //t3R->DrawLatex(0.93, 0.72,"#bf{#Xi(1530)^{0} #rightarrow #Xi^{#mp} + #pi^{#pm}}");

    auto lSignalFit = new TLegend(0.46, 0.3, 0.97, 0.48);
    lSignalFit->SetFillStyle(0);
    lSignalFit->AddEntry(hfitsignal, "Data (stat. uncert.)", "PLE");
    lSignalFit->AddEntry(fFitSum, Form("Voigtian peak + background"), "L");
    lSignalFit->AddEntry(fFitbkg, Form("Residual background"), "L");
    lSignalFit->Draw();

    SaveCanvas(cSignalFit,"figure_finalfit_preliminary","figs/Approval/");

    // 0-100 case:
    TFile* fSignalFit0100 = new TFile(
        "data/AnalysisResults_Extracted_1_Multi_0.00-100.00_Default1.root");

    auto hfitsignal0100 =
        (TH1*)fSignalFit0100->Get("hSignalBkgSubtraction_3");     // signal
    auto fFitSum0100 = (TF1*)fSignalFit0100->Get("fDataFitResult_3");  // Fit
    auto fFitbkg0100 = (TF1*)fSignalFit0100->Get("fDataFitBkgResult_3");  // Bkg

    hfitsignal0100->SetMarkerStyle(20);
    hfitsignal0100->SetMarkerSize(1);
    hfitsignal0100->SetMaximum(hfitsignal0100->GetMaximum()*1.3);
    hfitsignal0100->SetMarkerColor(kBlack);
    hfitsignal0100->SetLineColor(kBlack);
    hfitsignal0100->GetYaxis()->SetTitleOffset(1.2);
    hfitsignal0100->GetXaxis()->SetTitleOffset(1.2);
    hfitsignal0100->GetYaxis()->SetTitleSize(0.04);
    hfitsignal0100->GetXaxis()->SetTitleSize(0.04);
    hfitsignal0100->GetXaxis()->SetRangeUser(fDrawRange[0], fDrawRange[1]);
    hfitsignal0100->GetYaxis()->SetTitle("Counts / (4.0 MeV/#it{c} ^{2})");
    hfitsignal0100->GetXaxis()->SetTitle(
        "#it{M}_{#Xi#pi} (GeV/#it{c}^{2})");

    fFitbkg0100->SetLineColor(kBlack);
    fFitbkg0100->SetLineWidth(2);
    fFitbkg0100->SetRange(1.484, 1.6);
    fFitSum0100->SetLineColor(MaterialColors[0]);
    fFitSum0100->SetLineWidth(2);
    fFitSum0100->SetNpx(1000);

    cSignalFit->cd();
    hfitsignal0100->Draw("PZ");
    fFitbkg0100->Draw("same");
    fFitSum0100->Draw("same");

        t2R->DrawLatex(0.62, 0.91, "ALICE Preliminary");
    t2R->DrawLatex(0.93, 0.91, "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV}");
    //t2R->DrawLatex(0.93, 0.85, "#bf{V0M Multiplicity Event Classes 0 - 100% }");
    t2R->DrawLatex(0.93, 0.85,
                  "#bf{1.6 < #it{p}_{T} < 2.0 GeV/#it{c}, |#it{y}| < 0.5}");

    t2R->DrawLatex(0.93, 0.76,
                  "#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} #rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}");
    //t3R->DrawLatex(0.93, 0.72,"#bf{#Xi(1530)^{0} #rightarrow #Xi^{#mp} + #pi^{#pm}}");

    auto lSignalFit0100 = new TLegend(0.46, 0.3, 0.97, 0.48);
    lSignalFit0100->SetFillStyle(0);
    lSignalFit0100->AddEntry(hfitsignal0100, "Data (stat. uncert.)", "PLE");
    lSignalFit0100->AddEntry(fFitSum0100, Form("Voigtian peak + background"),
                             "L");
    lSignalFit0100->AddEntry(fFitbkg0100, Form("Residual background"), "L");
    lSignalFit0100->Draw();

    SaveCanvas(cSignalFit, "figure_finalfitMB_preliminary", "figs/Approval/");

}
void DrawRecEffi(){
    TCanvas* cRecEffi = GetCanvas("cRecEffi");
    cRecEffi->SetLeftMargin(0.15);
    cRecEffi->SetBottomMargin(0.13);
    cRecEffi->SetLogy();
    cRecEffi->Draw();
    TFile* fResults = new TFile("data/AnalysisResults_Extracted_1_Multi_0.00-100.00_Default1.root");

    auto hSpectraSys =
        (TH1*)fResults->Get("hMCReconEffi");  // sys
    hSpectraSys->SetMaximum(1);
    hSpectraSys->SetMinimum(1e-3);
    hSpectraSys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSpectraSys->GetYaxis()->SetTitle("Acceptance #times Efficiency #times B.R.");
    hSpectraSys->GetXaxis()->SetRangeUser(0.1, 8.8);
    hSpectraSys->GetYaxis()->SetTitleOffset(1.2);
    hSpectraSys->GetYaxis()->SetLabelSize(0.05);
    hSpectraSys->GetXaxis()->SetLabelSize(0.05);
    hSpectraSys->GetYaxis()->SetTitleSize(0.06);
    hSpectraSys->GetXaxis()->SetTitleSize(0.06);
    hSpectraSys->GetXaxis()->SetTitleOffset(0.95);
    hSpectraSys->SetLineWidth(2);
    hSpectraSys->SetLineColor(kBlack);
    hSpectraSys->SetMarkerColor(kBlack);
    hSpectraSys->SetMarkerStyle(20);
    hSpectraSys->SetMarkerSize(1.5);

    cRecEffi->cd();
    hSpectraSys->Draw("E");

    t3R->DrawLatex(0.93, 0.42, "ALICE Preliminary");
    t2R->DrawLatex(0.93, 0.34, "#bf{pp #sqrt{#it{s}} = 13 TeV, INEL>0, |#it{y}| < 0.5}");
    t3R->DrawLatex(0.93, 0.28,
                  "#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} #rightarrow "
                  "#Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}");

    //t2->DrawLatex(0.5, 0.15, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");
    
    SaveCanvas(cRecEffi,"figure_RecEffi_preliminary","figs/Approval/");
}
void DrawSpectraMB(){
    TCanvas* cSpectra = GetCanvas("cSpectra",850,900);
    cSpectra->SetLeftMargin(0.18);
    cSpectra->SetLogy();
    cSpectra->Draw();
    TFile* fSpectra = new TFile("AnalysisResults_Xi1530_YieldMean_INEL.root");

    auto hSpectraSys =
        (TH1*)fSpectra->Get("INEL_SYS_corrected");  // sys
    auto hSpectraStat =
        (TH1*)fSpectra->Get("INEL_stat_corrected");  // signal
    auto hSpectraLeviFit =
        (TF1*)fSpectra->Get("INEL_kFitLevi_0.80-4.80");  // signal

    hSpectraSys->GetXaxis()->SetRangeUser(0, 10.0);

    TGraphErrors* gSpectraSys = new TGraphErrors(hSpectraSys);
    TGraphErrors* gSpectraStat = new TGraphErrors(hSpectraStat);

    for (int j = 0; j < gSpectraStat->GetN(); j++) {
        gSpectraStat->SetPointError(j, 0, gSpectraStat->GetErrorY(j));
    }
    gSpectraSys->GetYaxis()->SetTitle(
        "1/#it{N}_{event}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    gSpectraSys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gSpectraSys->GetXaxis()->SetTitleSize(0.04);
    gSpectraSys->GetXaxis()->SetLabelSize(0.04);
    gSpectraSys->GetXaxis()->SetTitleOffset(1.2);
    gSpectraSys->GetYaxis()->SetTitleSize(0.04);
    gSpectraSys->GetYaxis()->SetLabelSize(0.04);
    gSpectraSys->GetYaxis()->SetTitleOffset(1.8);
    //gSpectraSys->GetXaxis()->SetTitleOffset(0.85);

    gSpectraSys->GetYaxis()->SetRangeUser(1e-7, 1e-2);
    gSpectraSys->SetLineColor(kBlack);
    gSpectraSys->SetMarkerColor(kBlack);
    gSpectraSys->SetMarkerSize(1);
    gSpectraSys->SetMarkerStyle(0);
    gSpectraSys->SetFillColor(0);
    gSpectraSys->SetLineWidth(1);
    gSpectraStat->SetLineColor(kBlack);
    gSpectraStat->SetMarkerColor(kBlack);
    gSpectraStat->SetMarkerStyle(20);
    gSpectraStat->SetMarkerSize(1);
    gSpectraStat->SetFillColor(kBlack);
    gSpectraStat->SetLineWidth(1);

    hSpectraLeviFit->SetLineColor(MaterialColors[0]);
    hSpectraLeviFit->SetLineWidth(2);
    hSpectraLeviFit->SetLineStyle(2);
    hSpectraLeviFit->SetRange(0,10);

    gSpectraSys->Draw("A5");
    gSpectraStat->Draw("P");
    hSpectraLeviFit->Draw("same");

    t2R->DrawLatex(0.95, 0.92, "ALICE Preliminary");
    t2R->DrawLatex(0.95, 0.87, "#bf{pp (INEL) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    //t0R->DrawLatex(0.95, 0.81,"#bf{V0M Multiplicity Event Classes 0 - 100%}");
    t0R->DrawLatex(0.95, 0.81, "#bf{Uncertainties: stat.(bar), syst.(box)}");

    t4R->DrawLatex(0.93, 0.72, "#Xi(1530)^{0}");
    // t2->DrawLatex(0.64, 0.69,"#bf{#Xi(1530)^{0} #rightarrow #Xi^{#mp} + #pi^{#pm}}");
    //t2->DrawLatex(0.15, 0.15, "#bf{Uncertainties: stat.(bar), syst.(box)}");

    gSpectraSys->SetMarkerSize(1);
    gSpectraSys->SetMarkerStyle(20);

    auto lSpectra = new TLegend(0.22, 0.18, 0.6, 0.33);
    lSpectra->SetFillStyle(0);
    lSpectra->AddEntry(gSpectraSys, "#frac{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0}}{2}",
                       "FP");
    lSpectra->AddEntry(hSpectraLeviFit, "L#acute{e}vy-Tsallis function", "L");
    lSpectra->Draw();

    SaveCanvas(cSpectra, "figure_Spectra_preliminary", "figs/Approval/");
}
void DrawSpectraMBRatio(){
    TCanvas* cSpectraRatio = GetCanvas("cSpectraRatio", 960, 450);
    cSpectraRatio->SetLeftMargin(0.10);
    cSpectraRatio->Draw();

    TFile* fSpectra = new TFile("AnalysisResults_Xi1530_YieldMean_INEL.root");
    auto hSpectraSys =
        (TH1*)fSpectra->Get("INEL_SYS_corrected");  // sys
    auto hSpectraStat =
        (TH1*)fSpectra->Get("INEL_stat_corrected");  // signal
    auto hSpectraSys_cor =
        (TH1*)fSpectra->Get("INEL_SYS_uncor_corrected");  // only pT dependent sys error.
    
    TFile* fSpectra7TeV = new TFile("AnalysisResults_Xi1530_YieldMean_7TeV.root");
    auto hSpectra7TeV_stat = (TH1*)fSpectra7TeV->Get("0.00-100.00_stat_corrected_7TeV");
    auto hSpectra7TeV_sys = (TH1*)fSpectra7TeV->Get("0.00-100.00_SYS_full_corrected_7TeV");
    auto hSpectra7TeV_sys_cor = (TH1*)fSpectra7TeV->Get("0.00-100.00_SYS_corrected_7TeV");

    const TArrayD* ptbin7TeV_array = hSpectra7TeV_stat->GetXaxis()->GetXbins();
    vector<double> ptbin;
    for (int i = 0; i < ptbin7TeV_array->GetSize(); i++)
        ptbin.push_back(ptbin7TeV_array->GetAt(i));

    TH1D* hRatio7TeV_stat = new TH1D("hRatio7TeV_stat", "hRatio7TeV_stat", ptbin7TeV_array->GetSize(), &ptbin[0]);
    TH1D* hRatio7TeV_sys = new TH1D("hRatio7TeV_sys", "hRatio7TeV_sys", ptbin7TeV_array->GetSize(), &ptbin[0]);
    TH1D* hRatio7TeV_sys_cor = new TH1D("hRatio7TeV_sys", "hRatio7TeV_sys", ptbin7TeV_array->GetSize(), &ptbin[0]);
    for (int i = 0; i < ptbin.size(); i++){
        double tempratio = hSpectraSys->GetBinContent(i+1)/hSpectra7TeV_stat->GetBinContent(i+1);
        hRatio7TeV_stat->SetBinContent(i + 1, tempratio);
        hRatio7TeV_sys->SetBinContent(i + 1, tempratio);
        hRatio7TeV_sys_cor->SetBinContent(i + 1, tempratio);

        double temperror_sys = tempratio
                                *sqrt( 
                                    pow(hSpectraSys->GetBinError(i+1)/hSpectraSys->GetBinContent(i+1),2)
                                    +pow(hSpectra7TeV_sys->GetBinError(i+1)/hSpectra7TeV_sys->GetBinContent(i+1),2)
                                    );
        double temperror_sys_cor = tempratio
                                *sqrt( 
                                    pow(hSpectraSys_cor->GetBinError(i+1)/hSpectraSys_cor->GetBinContent(i+1),2)
                                    +pow(hSpectra7TeV_sys_cor->GetBinError(i+1)/hSpectra7TeV_sys_cor->GetBinContent(i+1),2)
                                    );
        double temperror_stat = tempratio
                                *sqrt( 
                                    pow(hSpectraStat->GetBinError(i+1)/hSpectraStat->GetBinContent(i+1),2)
                                    +pow(hSpectra7TeV_stat->GetBinError(i+1)/hSpectra7TeV_stat->GetBinContent(i+1),2)
                                    );
        hRatio7TeV_stat->SetBinError(i + 1, temperror_stat);
        hRatio7TeV_sys->SetBinError(i + 1, temperror_sys);
        hRatio7TeV_sys_cor->SetBinError(i + 1, temperror_sys_cor);
    }
    hRatio7TeV_sys->GetXaxis()->SetRangeUser(0.5, 6);
    TGraphErrors* gRatio7TeV_stat = new TGraphErrors(hRatio7TeV_stat);
    TGraphErrors* gRatio7TeV_sys = new TGraphErrors(hRatio7TeV_sys);
    TGraphErrors* gRatio7TeV_sys_cor = new TGraphErrors(hRatio7TeV_sys_cor);
    for (int j = 0; j < gRatio7TeV_stat->GetN(); j++) {
        gRatio7TeV_stat->SetPointError(j, 0, gRatio7TeV_stat->GetErrorY(j));
    }
    gRatio7TeV_sys_cor->GetYaxis()->SetTitle(
        "Yield ratio to 7 TeV");
    gRatio7TeV_sys_cor->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gRatio7TeV_sys_cor->GetXaxis()->SetLimits(0.5,6);
    gRatio7TeV_sys_cor->GetYaxis()->SetRangeUser(0.7,3.5);
    gRatio7TeV_sys_cor->GetXaxis()->SetTitleSize(0.05);
    gRatio7TeV_sys_cor->GetXaxis()->SetLabelSize(0.05);
    gRatio7TeV_sys_cor->GetYaxis()->SetTitleSize(0.05);
    gRatio7TeV_sys_cor->GetYaxis()->SetLabelSize(0.05);
    gRatio7TeV_sys_cor->SetFillColorAlpha(MaterialColors[0], 0.3);
    //gRatio7TeV_sys_cor->SetFillStyle(3001);
    gRatio7TeV_sys_cor->SetLineColor(0);

    gRatio7TeV_sys->GetYaxis()->SetTitle(
        "Yield ratio to 7 TeV");
    gRatio7TeV_sys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gRatio7TeV_sys->GetXaxis()->SetLimits(0.5,6);
    gRatio7TeV_sys->GetYaxis()->SetRangeUser(0.7,2.7);
    gRatio7TeV_sys->GetXaxis()->SetTitleSize(0.05);
    gRatio7TeV_sys->GetXaxis()->SetLabelSize(0.05);
    gRatio7TeV_sys->GetYaxis()->SetTitleSize(0.05);
    gRatio7TeV_sys->GetYaxis()->SetLabelSize(0.05);

    gRatio7TeV_sys->SetMarkerColor(MaterialColors[0]);
    gRatio7TeV_sys->SetMarkerStyle(20);
    gRatio7TeV_sys->SetMarkerSize(0.5);
    //gRatio7TeV_sys->SetFillColor(0);
    gRatio7TeV_sys->SetLineColor(MaterialColors[0]);
    gRatio7TeV_sys->SetLineWidth(1);
    gRatio7TeV_sys->SetFillColorAlpha(MaterialColors[0], 0.0);
    gRatio7TeV_stat->SetMarkerColor(MaterialColors[0]);
    gRatio7TeV_stat->SetMarkerStyle(21);
    gRatio7TeV_stat->SetMarkerSize(1);
    gRatio7TeV_stat->SetFillColor(MaterialColors[0]);
    gRatio7TeV_stat->SetLineColor(MaterialColors[0]);
    gRatio7TeV_stat->SetLineWidth(1);
    
    cSpectraRatio->cd();
    gRatio7TeV_sys_cor->Draw("A5");
    gRatio7TeV_sys->Draw("5");
    gRatio7TeV_stat->Draw("P");

    TF1* tOneline =
        new TF1("tOneline", "1.0", -1, 100);
    tOneline->SetLineColor(kBlack);
    tOneline->SetLineWidth(1);
    tOneline->SetLineStyle(2);
    tOneline->Draw("same");

    auto lMBRatio = new TLegend(0.13, 0.5, 0.5, 0.7);
    lMBRatio->SetFillStyle(0);
    lMBRatio->AddEntry(gRatio7TeV_stat, "stat.", "LE");
    lMBRatio->AddEntry(gRatio7TeV_sys, "sys.","F");
    lMBRatio->AddEntry(gRatio7TeV_sys_cor, "uncorr. sys.", "F");
    lMBRatio->Draw();

    t3->DrawLatex(0.13, 0.86, "ALICE Preliminary");
    t3->DrawLatex(0.36, 0.86,
                  "#bf{pp (INEL) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    //t3->DrawLatex(0.13, 0.79, "#bf{Uncertainties: stat.(bar), syst.(box)}");

    t4->DrawLatex(0.13, 0.74, "#Xi(1530)^{0} + #bar{#Xi}(1530)^{0}");

    SaveCanvas(cSpectraRatio, "figure_SpectraRatio_preliminary", "figs/Approval/");
}
void DrawSpectraMulti(){
    vector<vector<double>> multibincheck = {
        {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    vector<int> Markerset = {20, 21, 33, 34, 29, 24, 25, 26, 27, 28,
                         30, 3,  5,  42, 43, 46, 47, 48, 49, 50,
                         51, 20, 21, 33, 34, 29, 24, 25, 26};
    vector<int> Markerset_size = {1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    vector<TFile*> buffer_file;
    vector<TH1*> buffer_multi_final_sys;
    vector<TH1*> buffer_multi_final_stat;
    vector<TF1*> buffer_multi_final_func;
    TFile* fSpectra = new TFile("AnalysisResults_Xi1530_YieldMean_0-100.root");

    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        TFile* fSpectraMulti = new TFile(Form("AnalysisResults_Xi1530_YieldMean_%.0f-%.0f.root",
                            multibincheck[imultibin][0],
                            multibincheck[imultibin][1]));
        buffer_file.push_back(fSpectraMulti);
        auto htemp_sys =
            (TH1*)buffer_file[imultibin]->Get(Form("%.2f-%.2f_SYS_corrected",
                            multibincheck[imultibin][0],
                            multibincheck[imultibin][1]));  // sys
        auto htemp_stat =
            (TH1*)buffer_file[imultibin]->Get(Form("%.2f-%.2f_stat_corrected",
                            multibincheck[imultibin][0],
                            multibincheck[imultibin][1]));  // sys
        auto htempLeviFit =
            (TF1*)buffer_file[imultibin]->Get(Form("%.2f-%.2f_kFitLevi_0.80-8.80",
                            multibincheck[imultibin][0],
                            multibincheck[imultibin][1]));  // sys
        buffer_multi_final_sys.push_back(htemp_sys);
        buffer_multi_final_stat.push_back(htemp_stat);
        buffer_multi_final_func.push_back(htempLeviFit);

    }

    //int sysColorPallet = GetSerialColors(buffer_multi_final_sys.size());
    

    TCanvas* c2 = new TCanvas("Final", "Final", 0, 0, 850, 1150);
    c2->Range(0, 0, 1, 1);

    auto denominatorHistogram_sys = 
            (TH1*)fSpectra->Get(Form("%.2f-%.2f_SYS_corrected", 0.0, 100.0));  // sys;
    auto denominatorHistogram_stat = 
            (TH1*)fSpectra->Get(Form("%.2f-%.2f_stat_corrected", 0.0, 100.0));  // sys;

    // Create ratio histograms
    vector<TH1*> hist_over_denomHist_sys;
    vector<TH1*> hist_over_denomHist_stat;
    vector<int> vcolors;
    for (int i = 0; i < buffer_multi_final_sys.size(); i++) {
        //vcolors.push_back(sysColorPallet + buffer_multi_final_sys.size() - i - 1);
        vcolors.push_back(MaterialColors[i]);
        hist_over_denomHist_sys.push_back((TH1*)buffer_multi_final_sys[i]->Clone());
        hist_over_denomHist_sys[i]->GetTitle();
        hist_over_denomHist_sys[i]->Divide(denominatorHistogram_sys);

        hist_over_denomHist_stat.push_back((TH1*)buffer_multi_final_stat[i]->Clone());
        hist_over_denomHist_stat[i]->GetTitle();
        hist_over_denomHist_stat[i]->Divide(denominatorHistogram_stat);
    }

    // Bottom plot
    c2->cd();
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
    c1_1->SetLogy(true);

    vector<TGraphAsymmErrors*> hdenomgsys;
    vector<TGraphAsymmErrors*> hdenomgstat;
    TGraphAsymmErrors* htmps1 =
        new TGraphAsymmErrors(hist_over_denomHist_sys[0]);
    TGraphAsymmErrors* htmpstat1 =
        new TGraphAsymmErrors(hist_over_denomHist_stat[0]);
    hdenomgsys.push_back(htmps1);
    hdenomgstat.push_back(htmpstat1);
    hdenomgsys[0]->SetLineColor(vcolors[0]);
    hdenomgsys[0]->SetMarkerColor(vcolors[0]);
    hdenomgsys[0]->SetMarkerSize(0);
    hdenomgsys[0]->SetMarkerStyle(0);
    hdenomgsys[0]->SetFillColor(0);
    hdenomgsys[0]->SetLineWidth(1);
    hdenomgsys[0]->GetYaxis()->SetNdivisions(5);
    hdenomgsys[0]->SetTitle(";#it{p}_{T} (GeV/#it{c});Raito to INEL>0");
    hdenomgsys[0]->GetXaxis()->SetRangeUser(0, 10);
    hdenomgsys[0]->GetXaxis()->SetTitleSize(0.1);
    hdenomgsys[0]->GetXaxis()->SetLabelSize(0.1);
    hdenomgsys[0]->GetXaxis()->SetTitleOffset(0.95);
    hdenomgsys[0]->GetYaxis()->SetLabelSize(0.1);
    hdenomgsys[0]->GetYaxis()->SetTitleSize(0.1);
    hdenomgsys[0]->GetYaxis()->SetTitleOffset(0.6);
    hdenomgsys[0]->SetMinimum(2e-2);
    hdenomgsys[0]->SetMaximum(7);
    for (int j = 0; j < hdenomgsys[0]->GetN(); j++) {
        hdenomgstat[0]->SetPointError(
            j, 0, 0,
            hdenomgstat[0]->GetErrorYhigh(j),
            hdenomgstat[0]->GetErrorYlow(j));
    }
    hdenomgstat[0]->SetLineColor(vcolors[0]);
    hdenomgstat[0]->SetMarkerColor(vcolors[0]);
    hdenomgstat[0]->SetMarkerStyle(20);
    hdenomgstat[0]->SetMarkerSize(0.5);
    hdenomgstat[0]->SetFillColor(vcolors[0]);
    hdenomgstat[0]->SetLineWidth(1);
    hdenomgstat[0]->SetMarkerSize(Markerset_size[0]);
    hdenomgstat[0]->SetMarkerStyle(Markerset[0]);
    hdenomgstat[0]->GetXaxis()->SetRangeUser(0, 10);
    hdenomgsys[0]->Draw("A5");
    hdenomgstat[0]->Draw("P");


    TF1* line1 = new TF1("line1", "1", -1, 15);
    line1->SetLineColor(1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");

    // hist_over_denomHist_stat[0]->SetMinimum(8e-2);
    // hist_over_denomHist_sys[0]->Draw("E2");
    // hist_over_denomHist_stat[0]->Draw("EP");
    hist_over_denomHist_sys[0]->SetLineWidth(1);
    hist_over_denomHist_sys[0]->SetLineColor(vcolors[0]);
    hist_over_denomHist_sys[0]->SetMarkerColor(vcolors[0]);

    // vector<TGraphAsymmErrors*> hdenomgsys;
    // vector<TGraphAsymmErrors*> hdenomgstat;
    for (int i = 1; i < buffer_multi_final_sys.size(); i++) {
        TGraphAsymmErrors* htmps =
            new TGraphAsymmErrors(hist_over_denomHist_sys[i]);
        TGraphAsymmErrors* htmpstat =
            new TGraphAsymmErrors(hist_over_denomHist_stat[i]);
        hdenomgsys.push_back(htmps);
        hdenomgstat.push_back(htmpstat);

        for (int j = 0; j < hdenomgstat[i]->GetN(); j++) {
            hdenomgstat[i]->SetPointError(j, 0, 0,
                                          hdenomgstat[i]->GetErrorYhigh(j),
                                          hdenomgstat[i]->GetErrorYlow(j));
            hdenomgsys[i]->SetPointError(j, hdenomgsys[i]->GetErrorXhigh(j),
                                         hdenomgsys[i]->GetErrorXlow(j),
                                         hdenomgsys[i]->GetErrorYhigh(j),
                                         hdenomgsys[i]->GetErrorYlow(j));
        }

        hdenomgsys[i]->SetLineColor(vcolors[i]);
        hdenomgsys[i]->SetMarkerColor(vcolors[i]);
        hdenomgsys[i]->SetMarkerSize(0);
        hdenomgsys[i]->SetMarkerStyle(0);
        hdenomgsys[i]->SetFillColor(0);
        hdenomgsys[i]->SetLineWidth(1);
        hdenomgstat[i]->SetLineColor(vcolors[i]);
        hdenomgstat[i]->SetMarkerColor(vcolors[i]);
        hdenomgstat[i]->SetMarkerStyle(20);
        hdenomgstat[i]->SetMarkerSize(0.5);
        hdenomgstat[i]->SetFillColor(vcolors[i]);
        hdenomgstat[i]->SetLineWidth(1);
        hdenomgstat[i]->SetMarkerSize(Markerset_size[i]);
        hdenomgstat[i]->SetMarkerStyle(Markerset[i]);

        hdenomgsys[i]->Draw("5");
        hdenomgstat[i]->Draw("P");
    }

    //*************************************************
    // Top Plot
    c2->cd();
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

    vector<TGraphAsymmErrors*> hgsys;
    vector<TGraphAsymmErrors*> hgstat;
    TGraphAsymmErrors* htmps2 =
        new TGraphAsymmErrors(buffer_multi_final_sys[0]);
    TGraphAsymmErrors* htmpstat2 =
        new TGraphAsymmErrors(buffer_multi_final_stat[0]);
    hgsys.push_back(htmps2);
    hgstat.push_back(htmpstat2);
    for (int j = 0; j < hgstat[0]->GetN(); j++) {
        hgstat[0]->GetY()[j] *= pow(2, buffer_multi_final_sys.size());
        hgstat[0]->SetPointError(
            j, 0, 0,
            hgstat[0]->GetErrorYhigh(j) * pow(2, buffer_multi_final_sys.size()),
            hgstat[0]->GetErrorYlow(j) * pow(2, buffer_multi_final_sys.size()));
        hgsys[0]->GetY()[j] *= pow(2, buffer_multi_final_sys.size());
        hgsys[0]->SetPointError(
            j, hgsys[0]->GetErrorXhigh(j), hgsys[0]->GetErrorXlow(j),
            hgsys[0]->GetErrorYhigh(j) * pow(2, buffer_multi_final_sys.size()),
            hgsys[0]->GetErrorYlow(j) * pow(2, buffer_multi_final_sys.size()));
    }
    hgsys[0]->SetLineColor(vcolors[0]);
    hgsys[0]->SetMarkerColor(vcolors[0]);
    hgsys[0]->SetMarkerSize(0);
    hgsys[0]->SetMarkerStyle(0);
    hgsys[0]->SetFillColor(0);
    hgsys[0]->SetLineWidth(1);
    hgsys[0]->GetXaxis()->SetTitleSize(0.00);
    hgsys[0]->GetYaxis()->SetLabelSize(0.04);
    hgsys[0]->GetYaxis()->SetTitleSize(0.05);
    hgsys[0]->GetYaxis()->SetTitleOffset(1.2);
    hgsys[0]->GetYaxis()->SetTitle(
        "1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hgsys[0]->GetXaxis()->SetRangeUser(0, 10);
    hgsys[0]->SetMaximum(5e-1);
    hgsys[0]->SetMinimum(5e-9);
    hgstat[0]->SetLineColor(vcolors[0]);
    hgstat[0]->SetMarkerColor(vcolors[0]);
    hgstat[0]->SetMarkerStyle(20);
    hgstat[0]->SetMarkerSize(0.5);
    hgstat[0]->SetFillColor(vcolors[0]);
    hgstat[0]->SetLineWidth(1);
    hgstat[0]->SetMarkerSize(Markerset_size[0]);
    hgstat[0]->SetMarkerStyle(Markerset[0]);

    hgsys[0]->Draw("A5");
    hgstat[0]->Draw("P");

    for (int i = 1; i < buffer_multi_final_sys.size(); i++) {
        TGraphAsymmErrors* htmps =
            new TGraphAsymmErrors(buffer_multi_final_sys[i]);
        TGraphAsymmErrors* htmpstat =
            new TGraphAsymmErrors(buffer_multi_final_stat[i]);
        hgsys.push_back(htmps);
        hgstat.push_back(htmpstat);
        for (int j = 0; j < hgstat[i]->GetN(); j++) {
            hgstat[i]->GetY()[j] *= pow(2, buffer_multi_final_sys.size() - i);
            hgstat[i]->SetPointError(
                j, 0, 0,
                hgstat[i]->GetErrorYhigh(j) *
                    pow(2, buffer_multi_final_sys.size() - i),
                hgstat[i]->GetErrorYlow(j) *
                    pow(2, buffer_multi_final_sys.size() - i));
            hgsys[i]->GetY()[j] *= pow(2, buffer_multi_final_sys.size() - i);
            hgsys[i]->SetPointError(
                j, hgsys[i]->GetErrorXhigh(j), hgsys[i]->GetErrorXlow(j),
                hgsys[i]->GetErrorYhigh(j) *
                    pow(2, buffer_multi_final_sys.size() - i),
                hgsys[i]->GetErrorYlow(j) *
                    pow(2, buffer_multi_final_sys.size() - i));
        }
        hgsys[i]->SetLineColor(vcolors[i]);
        hgsys[i]->SetMarkerColor(vcolors[i]);
        hgsys[i]->SetMarkerSize(0);
        hgsys[i]->SetMarkerStyle(0);
        hgsys[i]->SetFillColor(0);
        hgsys[i]->SetLineWidth(1);
        hgsys[i]->GetXaxis()->SetTitleSize(0.00);
        hgsys[i]->GetYaxis()->SetLabelSize(0.04);
        hgsys[i]->GetYaxis()->SetTitleSize(0.05);
        hgsys[i]->GetYaxis()->SetTitleOffset(1.2);
        hgsys[i]->GetYaxis()->SetTitle(
            "1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
        hgsys[i]->GetXaxis()->SetRangeUser(0, 10);
        hgsys[i]->SetMaximum(5e-1);
        hgsys[i]->SetMinimum(5e-9);
        hgstat[i]->SetLineColor(vcolors[i]);
        hgstat[i]->SetMarkerColor(vcolors[i]);
        hgstat[i]->SetMarkerStyle(20);
        hgstat[i]->SetMarkerSize(0.5);
        hgstat[i]->SetFillColor(vcolors[i]);
        hgstat[i]->SetLineWidth(1);
        hgstat[i]->SetMarkerSize(Markerset_size[i]);
        hgstat[i]->SetMarkerStyle(Markerset[i]);

        hgsys[i]->Draw("5");
        hgstat[i]->Draw("P");
    }
    // t->DrawLatex(0.42, 0.92, "ALICE Preliminary, #bf{pp #sqrt{#it{s}} = 13
    // TeV}");

    t2R->DrawLatex(0.71, 0.95, "ALICE Preliminary");
    t2R->DrawLatex(0.95, 0.95, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    t2R->DrawLatex(0.95, 0.89, "#bf{INEL>0, |#it{y}| < 0.5}");

    t4R->DrawLatex(0.95, 0.82, "#Xi(1530)^{0}");


    //t->DrawLatex(0.7, 0.87, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    //t->DrawLatex(0.82, 0.86, "#bf{|#it{y}| < 0.5}");
    //t4->DrawLatex(0.72, 0.78, "#Xi(1530)^{0}");
    t->DrawLatex(0.19, 0.28, "#bf{V0M Multiplicity Percentile(%)}");
    t->DrawLatex(0.19, 0.05, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    auto legend_hSpectra_ratio_final_spectra = new TLegend(.18, .1, .58, .27);
    legend_hSpectra_ratio_final_spectra->SetNColumns(2);
    legend_hSpectra_ratio_final_spectra->SetBorderSize(0);
    legend_hSpectra_ratio_final_spectra->SetFillStyle(0);

    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        buffer_multi_final_sys[imultibin]->SetLineColor(
            vcolors[imultibin]);
        buffer_multi_final_sys[imultibin]->SetMarkerColor(
            vcolors[imultibin]);
        buffer_multi_final_sys[imultibin]->SetMarkerStyle(
            Markerset[imultibin]);
        buffer_multi_final_sys[imultibin]->SetMarkerSize(
            Markerset_size[imultibin]);

        legend_hSpectra_ratio_final_spectra->AddEntry(
            buffer_multi_final_sys[imultibin],
            Form("%.0f - %.0f (x2^{%d})", multibincheck[imultibin][0],
                 multibincheck[imultibin][1],
                 int(buffer_multi_final_sys.size() - imultibin)),
            "PF");
    }
    legend_hSpectra_ratio_final_spectra->Draw();

    c1_2->SetLogy(true);

    SaveCanvas(c2, "figure_Spectra_Multi_preliminary", "figs/Approval/");

}
void DrawSpectraMultiPi(){
    vector<double> multibincheck = {0,1,5,10,15,20,30,40,50,70,100};
    vector<double> multibincheck_draw = {0,10,30,50,70,100};
    vector<int> Markerset = {20, 21, 33, 34, 29, 24, 25, 26, 27, 28,
                         30, 3,  5,  42, 43, 46, 47, 48, 49, 50,
                         51, 20, 21, 33, 34, 29, 24, 25, 26};
    vector<int> Markerset_size = {1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    vector<TFile*> buffer_file;
    vector<TH1*> buffer_multi_final_sys_temp;
    vector<TH1*> buffer_multi_final_stat_temp;
    TFile* fSpectra = new TFile("Final_combined_spectra_TPCTOFTOFonlyrTPCKinksITSsa_pp13TeV.root");

    for (int imultibin = 0; imultibin < multibincheck.size()-1; imultibin++) {
        TString temps = "hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pi_";
        if (multibincheck[imultibin] < 10)
            temps += Form("%1.f",multibincheck[imultibin]);
        else
            temps += Form("%02.f",multibincheck[imultibin]);
        
        if (multibincheck[imultibin+1] < 10)
            temps += Form("to%1.f",multibincheck[imultibin+1]);
        else
            temps += Form("to%02.f",multibincheck[imultibin+1]);
        auto htemp_sys =
            (TH1*)fSpectra->Get(Form("%s_syst",temps.Data()));  // sys
        auto htemp_stat =
            (TH1*)fSpectra->Get(Form("%s_stat",temps.Data()));  // sys
        
        buffer_multi_final_sys_temp.push_back(htemp_sys);
        buffer_multi_final_stat_temp.push_back(htemp_stat);

    }

    vector<TH1*> buffer_multi_final_sys;
    vector<TH1*> buffer_multi_final_stat;

    // 0-10
    TH1* h010sys = (TH1*)buffer_multi_final_sys_temp[0]->Clone();
    TH1* h010stat = (TH1*)buffer_multi_final_stat_temp[0]->Clone();
    
    TH1* temp15sys = (TH1*)buffer_multi_final_sys_temp[1]->Clone();
    TH1* temp15stat = (TH1*)buffer_multi_final_sys_temp[1]->Clone();
    TH1* temp510sys = (TH1*)buffer_multi_final_sys_temp[2]->Clone();
    TH1* temp510stat = (TH1*)buffer_multi_final_sys_temp[2]->Clone();

    temp15sys->Scale(4);
    temp15stat->Scale(4);
    temp510sys->Scale(5);
    temp510stat->Scale(5);

    h010sys->Add(temp15sys);
    h010sys->Add(temp510sys);
    h010stat->Add(temp15stat);
    h010stat->Add(temp510stat);

    h010stat->Scale(0.1);
    h010sys->Scale(0.1);
    buffer_multi_final_sys.push_back(h010sys);
    buffer_multi_final_stat.push_back(h010stat);

    //10-30
    TH1* h1030sys = (TH1*)buffer_multi_final_sys_temp[3]->Clone();
    TH1* h1030stat = (TH1*)buffer_multi_final_stat_temp[3]->Clone();
    
    TH1* temp1520sys = (TH1*)buffer_multi_final_sys_temp[4]->Clone();
    TH1* temp1520stat = (TH1*)buffer_multi_final_sys_temp[4]->Clone();
    TH1* temp2030sys = (TH1*)buffer_multi_final_sys_temp[5]->Clone();
    TH1* temp2030stat = (TH1*)buffer_multi_final_sys_temp[5]->Clone();

    h1030sys->Scale(5);
    h1030stat->Scale(5);
    temp1520sys->Scale(5);
    temp1520stat->Scale(5);
    temp2030sys->Scale(10);
    temp2030stat->Scale(10);

    h1030sys->Add(temp1520sys);
    h1030sys->Add(temp2030sys);
    h1030stat->Add(temp1520stat);
    h1030stat->Add(temp2030stat);

    h1030sys->Scale(0.05);
    h1030stat->Scale(0.05);
    buffer_multi_final_sys.push_back(h1030sys);
    buffer_multi_final_stat.push_back(h1030stat);
    
    //30-50
    TH1* h3050sys = (TH1*)buffer_multi_final_sys_temp[6]->Clone();
    TH1* h3050stat = (TH1*)buffer_multi_final_stat_temp[6]->Clone();
    
    TH1* temp4050sys = (TH1*)buffer_multi_final_sys_temp[7]->Clone();
    TH1* temp4050stat = (TH1*)buffer_multi_final_sys_temp[7]->Clone();

    h3050sys->Add(temp4050sys);
    h3050stat->Add(temp4050stat);

    h3050sys->Scale(0.5);
    h3050stat->Scale(0.5);
    buffer_multi_final_sys.push_back(h3050sys);
    buffer_multi_final_stat.push_back(h3050stat);

    //

    for (int imultibin = 8; imultibin < multibincheck.size()-1; imultibin++) {
        cout << imultibin << endl;
        buffer_multi_final_sys.push_back(buffer_multi_final_sys_temp[imultibin]);
        buffer_multi_final_stat.push_back(buffer_multi_final_stat_temp[imultibin]);
    }

    int sysColorPallet = GetSerialColors(buffer_multi_final_sys.size());
    

    TCanvas* c2 = new TCanvas("Final", "Final", 0, 0, 850, 1150);
    c2->Range(0, 0, 1, 1);

    auto denominatorHistogram_sys = 
            (TH1*)fSpectra->Get(Form("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pi_%1.fto%02.f_syst", 0.0, 100.0));  // sys;
    auto denominatorHistogram_stat = 
            (TH1*)fSpectra->Get(Form("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pi_%1.fto%02.f_stat", 0.0, 100.0));  // sys;
    
    // Create ratio histograms
    vector<TH1*> hist_over_denomHist_sys;
    vector<TH1*> hist_over_denomHist_stat;
    vector<int> vcolors;
    for (int i = 0; i < buffer_multi_final_sys.size(); i++) {
        cout << i << endl;
        vcolors.push_back(sysColorPallet + buffer_multi_final_sys.size() - i - 1);
        //vcolors.push_back(MaterialColors[i]);
        hist_over_denomHist_sys.push_back((TH1*)buffer_multi_final_sys[i]->Clone());
        hist_over_denomHist_sys[i]->GetTitle();
        hist_over_denomHist_sys[i]->Divide(denominatorHistogram_sys);

        hist_over_denomHist_stat.push_back((TH1*)buffer_multi_final_stat[i]->Clone());
        hist_over_denomHist_stat[i]->GetTitle();
        hist_over_denomHist_stat[i]->Divide(denominatorHistogram_stat);
    }
    cout << 1 << endl;
    
    // Bottom plot
    c2->cd();
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
    c1_1->SetLogy(true);

    vector<TGraphAsymmErrors*> hdenomgsys;
    vector<TGraphAsymmErrors*> hdenomgstat;
    TGraphAsymmErrors* htmps1 =
        new TGraphAsymmErrors(hist_over_denomHist_sys[0]);
    TGraphAsymmErrors* htmpstat1 =
        new TGraphAsymmErrors(hist_over_denomHist_stat[0]);
    hdenomgsys.push_back(htmps1);
    hdenomgstat.push_back(htmpstat1);
    hdenomgsys[0]->SetLineColor(vcolors[0]);
    hdenomgsys[0]->SetMarkerColor(vcolors[0]);
    hdenomgsys[0]->SetMarkerSize(0);
    hdenomgsys[0]->SetMarkerStyle(0);
    hdenomgsys[0]->SetFillColor(0);
    hdenomgsys[0]->SetLineWidth(1);
    hdenomgsys[0]->GetYaxis()->SetNdivisions(5);
    hdenomgsys[0]->SetTitle(";#it{p}_{T} (GeV/#it{c});Raito to INEL>0");
    hdenomgsys[0]->GetXaxis()->SetRangeUser(0, 9.5);
    hdenomgsys[0]->GetXaxis()->SetTitleSize(0.1);
    hdenomgsys[0]->GetXaxis()->SetLabelSize(0.1);
    hdenomgsys[0]->GetXaxis()->SetTitleOffset(0.95);
    hdenomgsys[0]->GetYaxis()->SetLabelSize(0.1);
    hdenomgsys[0]->GetYaxis()->SetTitleSize(0.1);
    hdenomgsys[0]->GetYaxis()->SetTitleOffset(0.6);
    hdenomgsys[0]->SetMinimum(6e-3);
    hdenomgsys[0]->SetMaximum(10);
    for (int j = 0; j < hdenomgsys[0]->GetN(); j++) {
        hdenomgstat[0]->SetPointError(
            j, 0, 0,
            hdenomgstat[0]->GetErrorYhigh(j),
            hdenomgstat[0]->GetErrorYlow(j));
    }
    hdenomgstat[0]->SetLineColor(vcolors[0]);
    hdenomgstat[0]->SetMarkerColor(vcolors[0]);
    hdenomgstat[0]->SetMarkerStyle(20);
    hdenomgstat[0]->SetMarkerSize(0.5);
    hdenomgstat[0]->SetFillColor(vcolors[0]);
    hdenomgstat[0]->SetLineWidth(1);
    hdenomgstat[0]->SetMarkerSize(Markerset_size[0]);
    hdenomgstat[0]->SetMarkerStyle(Markerset[0]);
    hdenomgstat[0]->GetXaxis()->SetRangeUser(0, 10);
    hdenomgsys[0]->Draw("A5");
    hdenomgstat[0]->Draw("P");


    TF1* line1 = new TF1("line1", "1", -1, 15);
    line1->SetLineColor(1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");

    // hist_over_denomHist_stat[0]->SetMinimum(8e-2);
    // hist_over_denomHist_sys[0]->Draw("E2");
    // hist_over_denomHist_stat[0]->Draw("EP");
    hist_over_denomHist_sys[0]->SetLineWidth(1);
    hist_over_denomHist_sys[0]->SetLineColor(vcolors[0]);
    hist_over_denomHist_sys[0]->SetMarkerColor(vcolors[0]);

    // vector<TGraphAsymmErrors*> hdenomgsys;
    // vector<TGraphAsymmErrors*> hdenomgstat;
    for (int i = 1; i < buffer_multi_final_sys.size(); i++) {
        TGraphAsymmErrors* htmps =
            new TGraphAsymmErrors(hist_over_denomHist_sys[i]);
        TGraphAsymmErrors* htmpstat =
            new TGraphAsymmErrors(hist_over_denomHist_stat[i]);
        hdenomgsys.push_back(htmps);
        hdenomgstat.push_back(htmpstat);

        for (int j = 0; j < hdenomgstat[i]->GetN(); j++) {
            hdenomgstat[i]->SetPointError(j, 0, 0,
                                          hdenomgstat[i]->GetErrorYhigh(j),
                                          hdenomgstat[i]->GetErrorYlow(j));
            hdenomgsys[i]->SetPointError(j, hdenomgsys[i]->GetErrorXhigh(j),
                                         hdenomgsys[i]->GetErrorXlow(j),
                                         hdenomgsys[i]->GetErrorYhigh(j),
                                         hdenomgsys[i]->GetErrorYlow(j));
        }

        hdenomgsys[i]->SetLineColor(vcolors[i]);
        hdenomgsys[i]->SetMarkerColor(vcolors[i]);
        hdenomgsys[i]->SetMarkerSize(0);
        hdenomgsys[i]->SetMarkerStyle(0);
        hdenomgsys[i]->SetFillColor(0);
        hdenomgsys[i]->SetLineWidth(1);
        hdenomgstat[i]->SetLineColor(vcolors[i]);
        hdenomgstat[i]->SetMarkerColor(vcolors[i]);
        hdenomgstat[i]->SetMarkerStyle(20);
        hdenomgstat[i]->SetMarkerSize(0.5);
        hdenomgstat[i]->SetFillColor(vcolors[i]);
        hdenomgstat[i]->SetLineWidth(1);
        hdenomgstat[i]->SetMarkerSize(Markerset_size[i]);
        hdenomgstat[i]->SetMarkerStyle(Markerset[i]);

        hdenomgsys[i]->Draw("5");
        hdenomgstat[i]->Draw("P");
    }

    //*************************************************
    // Top Plot
    c2->cd();
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

    vector<TGraphAsymmErrors*> hgsys;
    vector<TGraphAsymmErrors*> hgstat;
    TGraphAsymmErrors* htmps2 =
        new TGraphAsymmErrors(buffer_multi_final_sys[0]);
    TGraphAsymmErrors* htmpstat2 =
        new TGraphAsymmErrors(buffer_multi_final_stat[0]);
    hgsys.push_back(htmps2);
    hgstat.push_back(htmpstat2);
    for (int j = 0; j < hgstat[0]->GetN(); j++) {
        hgstat[0]->GetY()[j] *= pow(2, buffer_multi_final_sys.size());
        hgstat[0]->SetPointError(
            j, 0, 0,
            hgstat[0]->GetErrorYhigh(j) * pow(2, buffer_multi_final_sys.size()),
            hgstat[0]->GetErrorYlow(j) * pow(2, buffer_multi_final_sys.size()));
        hgsys[0]->GetY()[j] *= pow(2, buffer_multi_final_sys.size());
        hgsys[0]->SetPointError(
            j, hgsys[0]->GetErrorXhigh(j), hgsys[0]->GetErrorXlow(j),
            hgsys[0]->GetErrorYhigh(j) * pow(2, buffer_multi_final_sys.size()),
            hgsys[0]->GetErrorYlow(j) * pow(2, buffer_multi_final_sys.size()));
    }
    hgsys[0]->SetLineColor(vcolors[0]);
    hgsys[0]->SetMarkerColor(vcolors[0]);
    hgsys[0]->SetMarkerSize(0);
    hgsys[0]->SetMarkerStyle(0);
    hgsys[0]->SetFillColor(0);
    hgsys[0]->SetLineWidth(1);
    hgsys[0]->GetXaxis()->SetTitleSize(0.00);
    hgsys[0]->GetYaxis()->SetLabelSize(0.04);
    hgsys[0]->GetYaxis()->SetTitleSize(0.05);
    hgsys[0]->GetYaxis()->SetTitleOffset(1.2);
    hgsys[0]->GetYaxis()->SetTitle(
        "1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
    hgsys[0]->GetXaxis()->SetRangeUser(0, 9.5);
    hgsys[0]->SetMaximum(5e+3);
    hgsys[0]->SetMinimum(5e-6);
    hgstat[0]->SetLineColor(vcolors[0]);
    hgstat[0]->SetMarkerColor(vcolors[0]);
    hgstat[0]->SetMarkerStyle(20);
    hgstat[0]->SetMarkerSize(0.5);
    hgstat[0]->SetFillColor(vcolors[0]);
    hgstat[0]->SetLineWidth(1);
    hgstat[0]->SetMarkerSize(Markerset_size[0]);
    hgstat[0]->SetMarkerStyle(Markerset[0]);

    hgsys[0]->Draw("A5");
    hgstat[0]->Draw("P");

    for (int i = 1; i < buffer_multi_final_sys.size(); i++) {
        TGraphAsymmErrors* htmps =
            new TGraphAsymmErrors(buffer_multi_final_sys[i]);
        TGraphAsymmErrors* htmpstat =
            new TGraphAsymmErrors(buffer_multi_final_stat[i]);
        hgsys.push_back(htmps);
        hgstat.push_back(htmpstat);
        for (int j = 0; j < hgstat[i]->GetN(); j++) {
            hgstat[i]->GetY()[j] *= pow(2, buffer_multi_final_sys.size() - i);
            hgstat[i]->SetPointError(
                j, 0, 0,
                hgstat[i]->GetErrorYhigh(j) *
                    pow(2, buffer_multi_final_sys.size() - i),
                hgstat[i]->GetErrorYlow(j) *
                    pow(2, buffer_multi_final_sys.size() - i));
            hgsys[i]->GetY()[j] *= pow(2, buffer_multi_final_sys.size() - i);
            hgsys[i]->SetPointError(
                j, hgsys[i]->GetErrorXhigh(j), hgsys[i]->GetErrorXlow(j),
                hgsys[i]->GetErrorYhigh(j) *
                    pow(2, buffer_multi_final_sys.size() - i),
                hgsys[i]->GetErrorYlow(j) *
                    pow(2, buffer_multi_final_sys.size() - i));
        }
        hgsys[i]->SetLineColor(vcolors[i]);
        hgsys[i]->SetMarkerColor(vcolors[i]);
        hgsys[i]->SetMarkerSize(0);
        hgsys[i]->SetMarkerStyle(0);
        hgsys[i]->SetFillColor(0);
        hgsys[i]->SetLineWidth(1);
        hgsys[i]->GetXaxis()->SetTitleSize(0.00);
        hgsys[i]->GetYaxis()->SetLabelSize(0.04);
        hgsys[i]->GetYaxis()->SetTitleSize(0.05);
        hgsys[i]->GetYaxis()->SetTitleOffset(1.2);
        hgsys[i]->GetYaxis()->SetTitle(
            "1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
        hgsys[i]->GetXaxis()->SetRangeUser(0, 10);
        hgsys[i]->SetMaximum(5e-1);
        hgsys[i]->SetMinimum(5e-9);
        hgstat[i]->SetLineColor(vcolors[i]);
        hgstat[i]->SetMarkerColor(vcolors[i]);
        hgstat[i]->SetMarkerStyle(20);
        hgstat[i]->SetMarkerSize(0.5);
        hgstat[i]->SetFillColor(vcolors[i]);
        hgstat[i]->SetLineWidth(1);
        hgstat[i]->SetMarkerSize(Markerset_size[i]);
        hgstat[i]->SetMarkerStyle(Markerset[i]);

        hgsys[i]->Draw("5");
        hgstat[i]->Draw("P");
    }
    // t->DrawLatex(0.42, 0.92, "ALICE Preliminary, #bf{pp #sqrt{#it{s}} = 13
    // TeV}");

    t2R->DrawLatex(0.71, 0.95, "ALICE Preliminary");
    t2R->DrawLatex(0.95, 0.95, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    t2R->DrawLatex(0.95, 0.89, "#bf{INEL>0, |#it{y}| < 0.5}");


    //t->DrawLatex(0.7, 0.87, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    //t->DrawLatex(0.82, 0.86, "#bf{|#it{y}| < 0.5}");
    //t4->DrawLatex(0.72, 0.78, "#Xi(1530)^{0}");
    t->DrawLatex(0.19, 0.28, "#bf{V0M Multiplicity Percentile(%)}");
    t->DrawLatex(0.19, 0.05, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    auto legend_hSpectra_ratio_final_spectra = new TLegend(.18, .1, .58, .27);
    legend_hSpectra_ratio_final_spectra->SetNColumns(2);
    legend_hSpectra_ratio_final_spectra->SetBorderSize(0);
    legend_hSpectra_ratio_final_spectra->SetFillStyle(0);

    for (int imultibin = 0; imultibin < buffer_multi_final_sys.size(); imultibin++) {
        buffer_multi_final_sys[imultibin]->SetLineColor(
            vcolors[imultibin]);
        buffer_multi_final_sys[imultibin]->SetMarkerColor(
            vcolors[imultibin]);
        buffer_multi_final_sys[imultibin]->SetMarkerStyle(
            Markerset[imultibin]);
        buffer_multi_final_sys[imultibin]->SetMarkerSize(
            Markerset_size[imultibin]);

        legend_hSpectra_ratio_final_spectra->AddEntry(
            buffer_multi_final_sys[imultibin],
            Form("%.0f - %.0f (x2^{%d})", multibincheck_draw[imultibin],
                 multibincheck_draw[imultibin+1],
                 int(buffer_multi_final_sys.size() - imultibin)),
            "PF");
    }
    legend_hSpectra_ratio_final_spectra->Draw();

    c1_2->SetLogy(true);

    SaveCanvas(c2, "figure_Spectra_Multi_Pi", "figs/Approval/");
}
void DrawdNdy(){
    double ppYmax = 18e-3;
    double ppYmin = 0;
    double pPbYmax = 35e-3;
    double pPbYmin = 0;
    double PbPbYmax = 3;
    double PbPbYmin = 5e-4;
    bool draw7TeV = false;

    TCanvas* cdNdy = GetCanvas("cdNdy", 960, 720);
    cdNdy->SetLeftMargin(0.10);
    cdNdy->Draw();

    TFile* fdNdy = new TFile(PhysicsFile.Data());
    auto fdNdy_stat = (TGraphErrors*)fdNdy->Get("gYield_stat");  // signal
    auto fdNdy_sys = (TGraphErrors*)fdNdy->Get("gYield_syse");         // Fit
    auto fdNdy_sys_cor = (TGraphErrors*)fdNdy->Get("gYield_sys_cor");         // Fit
    auto fdNdy_stat_7Tev = (TGraphErrors*)fdNdy->Get("gYield7TeV_stat");      // Bkg
    auto fdNdy_sys_7Tev = (TGraphErrors*)fdNdy->Get("gYield7TeV_syse");      // Bkg

    fdNdy_stat_7Tev->SetPointError(
                0, 0, fdNdy_stat_7Tev->GetErrorY(0));    
    /*
    for (int j = 0; j < fdNdy_stat->GetN(); j++) {
        fdNdy_stat->GetY()[j] *= 2;
        fdNdy_sys->GetY()[j] *= 2;
        fdNdy_sys_cor->GetY()[j] *= 2;
    }
    */
    fdNdy_sys_cor->SetMinimum(0.0);
    fdNdy_sys_cor->SetMarkerColor(MaterialColors[0]);
    fdNdy_sys_cor->SetFillColorAlpha(MaterialColors[0], 0.3);
    //fdNdy_sys_cor->SetFillStyle(3001);
    fdNdy_sys_cor->SetLineColor(0);
    fdNdy_sys_cor->GetYaxis()->SetLabelSize(0.05);
    fdNdy_sys_cor->GetYaxis()->SetTitleOffset(0.7);
    fdNdy_sys_cor->GetYaxis()->SetTitleSize(0.05);

    fdNdy_sys_cor->GetXaxis()->SetLabelSize(0.05);
    //fdNdy_sys_cor->GetXaxis()->SetTitleOffset(0.95);
    fdNdy_sys_cor->GetXaxis()->SetTitleSize(0.05);
    //fdNdy_sys->SetMaximum(0.030);
    fdNdy_sys_cor->SetMaximum(ppYmax);
    fdNdy_sys_cor->SetMarkerStyle(20);
    fdNdy_sys_cor->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    fdNdy_sys_cor->GetYaxis()->SetTitle("d#it{N}/d#it{y}");

    //fdNdy_stat->SetMarkerStyle(20);
    //fdNdy_stat->SetMarkerSize(1);
    fdNdy_sys->SetMarkerColor(MaterialColors[0]);
    fdNdy_sys->SetLineColor(MaterialColors[0]);
    fdNdy_sys->SetFillColorAlpha(MaterialColors[0], 0.0);
    //fdNdy_sys->SetMaximum(0.030);
    fdNdy_sys->SetMaximum(ppYmax);
    fdNdy_sys->SetMarkerStyle(20);
    fdNdy_sys->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    fdNdy_sys->GetYaxis()->SetTitle("d#it{N}/d#it{y}");



    fdNdy_stat->SetMarkerColor(MaterialColors[0]);
    fdNdy_stat->SetLineColor(MaterialColors[0]);

    fdNdy_stat_7Tev->SetMarkerColor(kBlack);
    fdNdy_stat_7Tev->SetLineColor(kBlack);

    fdNdy_sys_7Tev->SetMarkerColor(kBlack);
    fdNdy_sys_7Tev->SetLineColor(kBlack);
    fdNdy_sys_7Tev->SetMarkerStyle(20);
    fdNdy_sys_7Tev->SetFillColor(0);

    //0-100 bin
    fdNdy_sys_7Tev->SetMinimum(0.0);
    fdNdy_sys_7Tev->SetMarkerColor(MaterialColors[0]);
    fdNdy_sys_7Tev->SetFillColorAlpha(MaterialColors[0], 0.3);
    //fdNdy_sys_7Tev->SetFillStyle(3001);
    fdNdy_sys_7Tev->SetLineColor(0);
    fdNdy_sys_7Tev->GetYaxis()->SetLabelSize(0.05);
    fdNdy_sys_7Tev->GetYaxis()->SetTitleOffset(0.6);
    fdNdy_sys_7Tev->GetYaxis()->SetTitleSize(0.06);

    fdNdy_sys_7Tev->GetXaxis()->SetLabelSize(0.05);
    //fdNdy_sys_7Tev->GetXaxis()->SetTitleOffset(0.95);
    fdNdy_sys_7Tev->GetXaxis()->SetTitleSize(0.05);
    //fdNdy_sys->SetMaximum(0.030);
    fdNdy_sys_7Tev->SetMarkerStyle(20);
    fdNdy_sys_7Tev->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    fdNdy_sys_7Tev->GetYaxis()->SetTitle("d#it{N}/d#it{y}");

    //fdNdy_stat->SetMarkerStyle(20);
    //fdNdy_stat->SetMarkerSize(1);
    fdNdy_sys_7Tev->SetMarkerColor(kBlack);
    fdNdy_sys_7Tev->SetLineColor(kBlack);
    fdNdy_sys_7Tev->SetFillColorAlpha(kBlack, 0.0);
    //fdNdy_sys->SetMaximum(0.030);
    fdNdy_sys_7Tev->SetMaximum(ppYmax);
    fdNdy_sys_7Tev->SetMarkerStyle(20);
    fdNdy_sys_7Tev->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    fdNdy_sys_7Tev->GetYaxis()->SetTitle("d#it{N}/d#it{y}");

    TFile* fdNdy0100 = new TFile("AnalysisResults_Xi1530_YieldMean_0-100.root");
    auto fdNdy_YieldMean0100 = (TH1D*)fdNdy0100->Get("hYield_0.00-100.00");  // signal
    vector<double> yield0100;
    yield0100.push_back(fdNdy_YieldMean0100->GetBinContent(1));
    vector<double> yield0100state;
    yield0100state.push_back(fdNdy_YieldMean0100->GetBinContent(2));
    vector<double> yield0100syse;
    double ftemp_yieldsyse = 0.0;
    ftemp_yieldsyse += pow(fdNdy_YieldMean0100->GetBinContent(3),2);
    ftemp_yieldsyse += pow(fdNdy_YieldMean0100->GetBinContent(4),2);
    yield0100syse.push_back(sqrt(ftemp_yieldsyse)/2);
    vector<double> yield0100syse_full;
    ftemp_yieldsyse += pow(0.0538*yield0100[0],2);
    yield0100syse_full.push_back(sqrt(ftemp_yieldsyse)/2);

    vector<double> x0100 = {6.94};
    vector<double> x0100e = {0.10};
    vector<double> x0100zero = {0};

    TGraphErrors* ge_stat_0100 = new TGraphErrors(yield0100.size(), &x0100[0], &yield0100[0], &x0100zero[0], &yield0100state[0]);
    TGraphErrors* ge_sys_0100 = new TGraphErrors(yield0100.size(), &x0100[0], &yield0100[0], &x0100e[0], &yield0100syse[0]);
    TGraphErrors* ge_sys_full_0100 = new TGraphErrors(yield0100.size(), &x0100[0], &yield0100[0], &x0100e[0], &yield0100syse_full[0]);
    ge_sys_0100->SetLineColor(0);
    ge_sys_0100->SetMarkerColor(MaterialColors[0]);
    ge_sys_0100->SetFillColorAlpha(MaterialColors[0], 0.3);
    ge_sys_0100->SetMarkerStyle(20);
    ge_sys_0100->SetMarkerSize(0.5);

    ge_stat_0100->SetMinimum(0.0);
    ge_stat_0100->SetMarkerStyle(20);
    ge_stat_0100->SetMarkerSize(0.5);
    ge_stat_0100->SetMarkerColor(MaterialColors[0]);
    ge_stat_0100->SetLineColor(MaterialColors[0]);
    ge_sys_full_0100->SetMarkerColor(MaterialColors[0]);
    ge_sys_full_0100->SetMarkerStyle(20);
    //ge_sys_full_0100->SetMarkerSize(0.5);

    ge_sys_full_0100->SetLineColor(MaterialColors[0]);
    ge_sys_full_0100->SetMinimum(0.0);
    ge_sys_full_0100->SetMaximum(0.01);
    ge_sys_full_0100->SetFillColorAlpha(MaterialColors[3], 0.0);
    ge_sys_full_0100->GetYaxis()->SetTitleOffset(0.7);
    ge_sys_full_0100->GetXaxis()->SetLimits(0,30);

    
    //

    cdNdy->cd();
    fdNdy_sys_cor->Draw("a5");
    fdNdy_sys->Draw("5");
    fdNdy_stat->Draw("P");
    if(draw7TeV){
        fdNdy_sys_7Tev->Draw("5");
        fdNdy_stat_7Tev->Draw("P");
    }  
    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t->DrawLatex(0.5, 0.17, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* ldNdy;
    if(draw7TeV) 
        ldNdy = new TLegend(0.13, 0.56, 0.5, 0.73);
    else
        ldNdy = new TLegend(0.13, 0.64, 0.5, 0.73); 
    ldNdy->SetFillStyle(0);
    ldNdy->AddEntry(fdNdy_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) ldNdy->AddEntry(fdNdy_sys_7Tev, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) ldNdy->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    ldNdy->Draw();
    
    SaveCanvas(cdNdy,"figure_dNdy_preliminary","figs/Approval/");


    // Further figure
    TCanvas* cdNdy2 = GetCanvas("cdNdy2", 960, 720);
    cdNdy2->SetLeftMargin(0.12);
    cdNdy2->Draw();
    // pPb Result
    // https://www.hepdata.net/record/ins1510878
    vector<double> xpPb502 = {35.6, 23.2, 16.1, 7.1};
    vector<double> xpPb502e = {0.8, 0.5, 0.4, 0.2};
    vector<double> xpPb502zero = {0, 0, 0, 0};
     
    // HEP data Xi(1530) pPb5.02TeV
    vector<double> ypPb502 = {2.73e-2, 1.77e-2, 1.07e-2, 3.6e-3};
    vector<double> ypPb502e = {0.6e-3, 0.5e-3, 0.3e-3, 0.1e-3};
    vector<double> ypPb502se = {0.0028, 0.0024, 0.0016, 0.0005};
    
    TGraphErrors* ge_stat_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502zero[0], &ypPb502e[0]);
    TGraphErrors* ge_sys_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502e[0], &ypPb502se[0]);
    
    /*
    for (int j = 0; j < ge_stat_pPb502TeV->GetN(); j++) {
        ge_stat_pPb502TeV->GetY()[j] /= 1.5;
        ge_sys_pPb502TeV->GetY()[j] /= 1.5;
    }
    */
    
    ge_sys_pPb502TeV->SetLineColor(MaterialColors[3]);
    ge_sys_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    ge_sys_pPb502TeV->SetMarkerStyle(20);
    //ge_sys_pPb502TeV->SetMarkerSize(0.5);

    ge_stat_pPb502TeV->SetLineColor(MaterialColors[3]);
    ge_stat_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    ge_stat_pPb502TeV->SetMarkerStyle(20);
    ge_stat_pPb502TeV->SetMarkerSize(0.5);
    
    auto fdNdy_sys_cor2 = (TGraphErrors*)fdNdy_sys_cor->Clone();
    fdNdy_sys_cor2->SetMaximum(pPbYmax);
    fdNdy_sys_cor2->SetMinimum(pPbYmin);
    fdNdy_sys_cor2->GetYaxis()->SetTitleOffset(0.9);
    fdNdy_sys_cor2->GetXaxis()->SetLimits(1,46);

    //fdNdy_sys_7Tev->GetYaxis()->SetTitle("d#it{N}/d#it{y}}");
    fdNdy_sys_7Tev->GetYaxis()->SetTitleOffset(1.0);
    fdNdy_sys_7Tev->SetMaximum(pPbYmax);
    fdNdy_sys_7Tev->SetMinimum(pPbYmin);
    fdNdy_sys_7Tev->GetXaxis()->SetLimits(0,5e1);
    
    cdNdy2->cd();

    fdNdy_sys_cor2->Draw("A5");
    fdNdy_sys->Draw("5");
    fdNdy_stat->Draw("P");
    ge_sys_pPb502TeV->Draw("5");
    ge_stat_pPb502TeV->Draw("P");

    if(draw7TeV){
        fdNdy_sys_7Tev->Draw("5");
        fdNdy_stat_7Tev->Draw("P");
    }  

    //ge_sys_0100->Draw("5");
    //ge_sys_full_0100->Draw("5");
    //ge_stat_0100->Draw("P");

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* ldNdy2;
    if(draw7TeV)
        ldNdy2 = new TLegend(0.55, 0.20, 0.95, 0.45);
    else
        ldNdy2 = new TLegend(0.55, 0.20, 0.95, 0.35);
    ldNdy2->SetFillStyle(0);
    //ldNdy2->AddEntry(fdNdy_sys, "#Xi(1530)^{0} pp 13 TeV (INEL>0)", "PEF");
    ldNdy2->AddEntry(fdNdy_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) ldNdy2->AddEntry(fdNdy_sys_7Tev, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) ldNdy2->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    ldNdy2->AddEntry(ge_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    ldNdy2->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    ldNdy2->Draw();
    //cdNdy2->SetLogx();
    //cdNdy2->SetLogy();
    SaveCanvas(cdNdy2,"figure_dNdy_withpPb502_preliminary","figs/Approval/");

    // PbPb Result
    TCanvas* cdNdyPbPb = GetCanvas("cdNdy2", 960, 720);
    cdNdyPbPb->SetLeftMargin(0.10);
    cdNdyPbPb->SetBottomMargin(0.15);
    cdNdyPbPb->Draw();
    cdNdyPbPb->cd();
    // charged particle multiplicity
    // https://arxiv.org/pdf/1012.1657.pdf
    vector<double> xPbPb276 = {1447.5, 680.3, 130.25};
    vector<double> xPbPb276e = {54.5, 25, 5.25};
    vector<double> xPbPb276zero = {0, 0, 0};
     
    // Analysis note PbPb 2.76
    vector<double> yPbPb276 = {686.4e-3, 338.7e-3, 82e-3};
    vector<double> yPbPb276e = {18.2e-3, 5.9e-3, 2.5e-3};
    vector<double> yPbPb276se = {228.5e-3, 80.4e-3, 21.4e-3};
    
    TGraphErrors* ge_stat_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276zero[0], &yPbPb276e[0]);
    TGraphErrors* ge_sys_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276e[0], &yPbPb276se[0]);
    
    fdNdy_sys_cor2->GetYaxis()->SetTitleOffset(0.6);
    fdNdy_sys_cor2->GetXaxis()->SetTitleOffset(1.35);
    fdNdy_sys_cor2->GetXaxis()->SetLimits(2,5e3);
    fdNdy_sys_cor2->SetMaximum(PbPbYmax);
    fdNdy_sys_cor2->SetMinimum(PbPbYmin);

    fdNdy_sys_cor2->SetMaximum(PbPbYmax);
    fdNdy_sys_cor2->SetMinimum(PbPbYmin);

    fdNdy_sys_cor2->Draw("A5");
    fdNdy_sys->Draw("5");
    fdNdy_stat->Draw("P");
    ge_sys_pPb502TeV->Draw("5");
    ge_stat_pPb502TeV->Draw("P");

    //ge_sys_0100->Draw("5");
    //ge_sys_full_0100->Draw("5");
    //ge_stat_0100->Draw("P");
    
    if(draw7TeV){
        fdNdy_sys_7Tev->Draw("5");
        fdNdy_stat_7Tev->Draw("P");
    }  
    

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");
    //t2->DrawLatex(0.5, 0.15, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");
    ge_sys_PbPb276TeV->Draw("5");
    ge_stat_PbPb276TeV->Draw("P");
    ge_sys_PbPb276TeV->SetMarkerStyle(20);
    //ge_sys_PbPb276TeV->SetMarkerSize(0.5);
    ge_stat_PbPb276TeV->SetMarkerStyle(20);
    ge_stat_PbPb276TeV->SetMarkerSize(0.5);
    ge_sys_PbPb276TeV->SetLineColor(MaterialColors[1]);
    ge_sys_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    ge_stat_PbPb276TeV->SetLineColor(MaterialColors[1]);
    ge_stat_PbPb276TeV->SetMarkerColor(MaterialColors[1]);

    fdNdy_sys_7Tev->GetYaxis()->SetTitleOffset(0.6);
    fdNdy_sys_7Tev->GetXaxis()->SetTitleOffset(1.35);
    fdNdy_sys_7Tev->GetXaxis()->SetLimits(2,5e3);
    //fdNdy_sys_7Tev->SetMaximum(PbPbYmax);
    //fdNdy_sys_7Tev->SetMinimum(PbPbYmin);
    cdNdyPbPb->SetLogx();
    cdNdyPbPb->SetLogy();

    TLegend* ldNdyPbPb;
    if(draw7TeV)
        ldNdyPbPb = new TLegend(0.6, 0.18, 0.95, 0.5);
    else
        ldNdyPbPb = new TLegend(0.6, 0.18, 0.95, 0.4);
    ldNdyPbPb->AddEntry(fdNdy_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) ldNdyPbPb->AddEntry(fdNdy_sys_7Tev, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) ldNdyPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    ldNdyPbPb->AddEntry(ge_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    ldNdyPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    ldNdyPbPb->AddEntry(ge_sys_PbPb276TeV, "Pb-Pb 2.76TeV (INEL)", "PEF");
    ldNdyPbPb->AddEntry((TObject*)0, "ALICE prelimiary", "");
    ldNdyPbPb->Draw();
    ldNdyPbPb->SetFillStyle(0);

    SaveCanvas(cdNdyPbPb,"figure_dNdy_withpPb502_preliminary_PbPb","figs/Approval/");
    
    // INEL case
    TFile* fINEL = new TFile("AnalysisResults_Xi1530_YieldMean_INEL.root");
    TH1D* hINEL = (TH1D*)fINEL->Get("hYield_INEL");

    // x pp 13TeV inel https://doi.org/10.1016/j.physletb.2015.12.030
    vector<double> xINEL = {5.31};
    vector<double> xINELe_zero = {0};
    vector<double> xINELe = {0.18};
    vector<double> yieldINEL = {};
    vector<double> yieldINEL_state = {};
    vector<double> yieldINEL_syshe = {};
    vector<double> yieldINEL_sysle = {};
    vector<double> yieldINEL_sys_cor = {};

    // INEL yield
    yieldINEL.push_back(hINEL->GetBinContent(1));
    yieldINEL_state.push_back(hINEL->GetBinContent(0));

    double tempfinalsyserror = 0;
    tempfinalsyserror += pow(hINEL->GetBinContent(2),2);
    tempfinalsyserror += pow(hINEL->GetBinContent(3),2);

    double averageerror_correl = sqrt(tempfinalsyserror)/2;
    double averageerror = sqrt(pow(averageerror_correl,2) + pow(0.0538*hINEL->GetBinContent(1),2));
    double averageerror_uncorr = 0.0538*hINEL->GetBinContent(1);

    yieldINEL_syshe.push_back(averageerror);
    yieldINEL_sysle.push_back(averageerror);
    yieldINEL_sys_cor.push_back(averageerror_uncorr);

    TGraphErrors* grdndyINEL_stat = new TGraphErrors(1, &xINEL[0], &yieldINEL[0], &xINELe_zero[0], &yieldINEL_state[0]);
    TGraphErrors* grdndyINEL_sys = new TGraphErrors(1, &xINEL[0], &yieldINEL[0], &xINELe[0], &yieldINEL_sysle[0]);
    TGraphErrors* grdndyINEL_sys_cor = new TGraphErrors(1, &xINEL[0], &yieldINEL[0], &xINELe[0], &yieldINEL_sys_cor[0]);
    grdndyINEL_sys->SetTitle("");
    grdndyINEL_sys->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    grdndyINEL_sys->GetYaxis()->SetTitle("d#it{N}_{#Xi^{*0}}/dy");
    grdndyINEL_sys->GetXaxis()->SetLimits(1,10);
    //ge_stat->GetYaxis()->SetRangeUser(0,0.007);
    grdndyINEL_sys->SetMinimum(0);
    grdndyINEL_sys->SetMaximum(0.005);
    grdndyINEL_stat->SetLineColor(MaterialColors[3]);
    grdndyINEL_stat->SetMarkerColor(MaterialColors[3]);
    grdndyINEL_stat->SetMarkerStyle(20);
    grdndyINEL_stat->SetMarkerSize(0.5);
    grdndyINEL_stat->SetFillColor(MaterialColors[3]);

    grdndyINEL_sys->SetLineColor(MaterialColors[3]);
    grdndyINEL_sys->SetMarkerColor(MaterialColors[3]);
    grdndyINEL_sys->SetFillColor(0);

    grdndyINEL_sys_cor->SetMarkerColor(MaterialColors[3]);
    grdndyINEL_sys_cor->SetFillColorAlpha(MaterialColors[3], 0.3);
    grdndyINEL_sys_cor->SetLineColor(0);
    grdndyINEL_sys_cor->SetMarkerStyle(20);

    TCanvas* cdNdy3 = GetCanvas("cdNdy3", 960, 720);
    cdNdy3->SetLeftMargin(0.10);
    cdNdy3->Draw();

    fdNdy_sys_7Tev->SetMaximum(0.006);
    fdNdy_sys_7Tev->GetXaxis()->SetLimits(4,9);

    cdNdy3->cd();
    
    fdNdy_sys_cor->Draw("a5");
    fdNdy_sys->Draw("5");
    fdNdy_stat->Draw("P");
    grdndyINEL_sys->Draw("5");
    grdndyINEL_sys_cor->Draw("5");
    fdNdy_sys_7Tev->Draw("5");
    fdNdy_stat_7Tev->Draw("P");
    grdndyINEL_stat->Draw("P");
    
    
    
    /*
    fdNdy_sys_7Tev->Draw("A5");
    fdNdy_stat_7Tev->Draw("P");
    ge_sys_0100->Draw("5");
    ge_sys_full_0100->Draw("5");
    ge_stat_0100->Draw("P");
    */
    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    auto ldNdy3 = new TLegend(0.6, 0.18, 0.95, 0.4);
    ldNdy3->SetFillStyle(0);
    //ldNdy3->AddEntry(fdNdy_sys, "#Xi(1530)^{0} pp 13 TeV", "PEF");
    ldNdy3->AddEntry(fdNdy_sys_cor, "#Xi(1530)^{0} pp (INEL>0) 13 TeV", "PEF");
    ldNdy3->AddEntry(grdndyINEL_sys_cor, "#Xi(1530)^{0} pp (inel) 13 TeV", "PEF");
    ldNdy3->AddEntry(fdNdy_sys_7Tev, "#Xi(1530)^{0} pp (inel) 7 TeV", "PEF");
    ldNdy3->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    ldNdy3->Draw();

    SaveCanvas(cdNdy3,"figure_dNdy_0100_preliminary","figs/Approval/");
    
}
void DrawMeanpT(){
    double meanpTYmax = 2.2;
    double meanpTYmin = 0.6;
    bool draw7TeV = false;
    
    TCanvas* cMeanPt = GetCanvas("cMeanPt", 960, 720);
    cMeanPt->SetLeftMargin(0.10);
    cMeanPt->Draw();

    TFile* fMeanPt = new TFile(PhysicsFile.Data());
    auto fMeanPt_stat = (TGraphErrors*)fMeanPt->Get("gMeanpT_stat");
    auto fMeanPt_sys = (TGraphErrors*)fMeanPt->Get("gMeanpT_syse");
    auto fMeanPt_sys_cor = (TGraphErrors*)fMeanPt->Get("gMeanpT_sys_cor");
    
    auto fMeanPt_stat_7TeV = (TGraphErrors*)fMeanPt->Get("gMeanpT7TeV_stat");
    auto fMeanPt_sys_7TeV = (TGraphErrors*)fMeanPt->Get("gMeanpT7TeV_syse");

    fMeanPt_sys->SetMarkerColor(MaterialColors[0]);
    fMeanPt_sys->SetLineColor(MaterialColors[0]);
    fMeanPt_sys->SetFillColorAlpha(MaterialColors[0], 0.0);

    fMeanPt_sys_cor->SetMarkerColor(MaterialColors[0]);
    fMeanPt_sys_cor->SetLineColor(0);
    fMeanPt_sys_cor->SetFillColorAlpha(MaterialColors[0], 0.3);
    fMeanPt_sys_cor->SetFillStyle(1000);
    fMeanPt_sys_cor->GetYaxis()->SetLabelSize(0.05);
    fMeanPt_sys_cor->GetYaxis()->SetTitleOffset(0.9);
    fMeanPt_sys_cor->GetYaxis()->SetTitleSize(0.05);
    fMeanPt_sys_cor->GetXaxis()->SetLabelSize(0.05);
    fMeanPt_sys_cor->GetXaxis()->SetTitleSize(0.05);
    fMeanPt_sys_cor->SetMaximum(meanpTYmax);
    fMeanPt_sys_cor->SetMinimum(meanpTYmin);
    fMeanPt_sys->SetMarkerStyle(20);
    fMeanPt_sys_cor->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    fMeanPt_sys_cor->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");

    fMeanPt_stat->SetMarkerColor(MaterialColors[0]);
    fMeanPt_stat->SetLineColor(MaterialColors[0]);
    fMeanPt_stat->SetMarkerStyle(20);
    fMeanPt_stat->SetMarkerSize(0.5);

    fMeanPt_stat_7TeV->SetMarkerColor(kBlack);
    fMeanPt_stat_7TeV->SetLineColor(kBlack);

    fMeanPt_sys_7TeV->SetMarkerColor(kBlack);
    fMeanPt_sys_7TeV->SetLineColor(kBlack);
    fMeanPt_sys_7TeV->SetMarkerStyle(20);
    fMeanPt_sys_7TeV->SetFillColor(0);

    fMeanPt_sys_7TeV->GetXaxis()->SetLimits(1,21);


    fMeanPt_sys_7TeV->GetYaxis()->SetLabelSize(0.05);
    fMeanPt_sys_7TeV->GetYaxis()->SetTitleOffset(0.9);
    fMeanPt_sys_7TeV->GetYaxis()->SetTitleSize(0.05);
    fMeanPt_sys_7TeV->GetXaxis()->SetLabelSize(0.05);
    fMeanPt_sys_7TeV->GetXaxis()->SetTitleSize(0.05);
    fMeanPt_sys_7TeV->SetMaximum(2.2);
    fMeanPt_sys_7TeV->SetMinimum(0.7);
    fMeanPt_sys_7TeV->SetMarkerStyle(20);
    fMeanPt_sys_7TeV->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    fMeanPt_sys_7TeV->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/c)");

    cMeanPt->cd();
    fMeanPt_sys_cor->Draw("A5");
    fMeanPt_sys->Draw("5");
    fMeanPt_stat->Draw("P");
    if(draw7TeV){
        fMeanPt_sys_7TeV->Draw("5");
        fMeanPt_stat_7TeV->Draw("P");
    }

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* lMeanPt;
    if(draw7TeV) 
        lMeanPt = new TLegend(0.6, 0.26, 0.95, 0.42);
    else 
        lMeanPt = new TLegend(0.6, 0.26, 0.95, 0.36);
    lMeanPt->SetFillStyle(0);
    lMeanPt->AddEntry(fMeanPt_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lMeanPt->AddEntry(fMeanPt_sys_7TeV, "pp 7 TeV", "PEF");
    if(draw7TeV) lMeanPt->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lMeanPt->Draw();

    SaveCanvas(cMeanPt,"figure_MeanPt_preliminary","figs/Approval/");

    // 0 - 100%
    TCanvas* cMeanPt2 = GetCanvas("cMeanPt2", 960, 720);
    cMeanPt2->SetLeftMargin(0.10);
    cMeanPt2->Draw();

    TFile* fdNdy0100 = new TFile("AnalysisResults_Xi1530_YieldMean_0-100.root");
    auto fdNdy_YieldMean0100 = (TH1D*)fdNdy0100->Get("hMeanPt_0.00-100.00");  // signal
    vector<double> x0100 = {6.94};
    vector<double> x0100e = {0.10};
    vector<double> x0100zero = {0};

    vector<double> mean0100;
    mean0100.push_back(fdNdy_YieldMean0100->GetBinContent(1));
    vector<double> mean0100state;
    mean0100state.push_back(fdNdy_YieldMean0100->GetBinContent(2));
    vector<double> mean0100syse;
    double ftemp_mean0100syse = 0.0;
    ftemp_mean0100syse += pow(fdNdy_YieldMean0100->GetBinContent(3),2);
    ftemp_mean0100syse += pow(fdNdy_YieldMean0100->GetBinContent(4),2);
    mean0100syse.push_back(sqrt(ftemp_mean0100syse)/2);
    vector<double> ftemp_mean0100syse_full;
    ftemp_mean0100syse += pow(0.0538*mean0100[0],2);
    ftemp_mean0100syse_full.push_back(sqrt(ftemp_mean0100syse)/2);

    TGraphErrors* gpt_stat_0100 = new TGraphErrors(mean0100.size(), &x0100[0], &mean0100[0], &x0100zero[0], &mean0100state[0]);
    TGraphErrors* gpt_sys_0100 = new TGraphErrors(mean0100.size(), &x0100[0], &mean0100[0], &x0100e[0], &mean0100syse[0]);
    TGraphErrors* gpt_sys_full_0100 = new TGraphErrors(mean0100.size(), &x0100[0], &mean0100[0], &x0100e[0], &ftemp_mean0100syse_full[0]);

    gpt_sys_0100->SetLineColor(0);
    gpt_sys_0100->SetMarkerColor(MaterialColors[0]);
    gpt_sys_0100->SetFillColor(kRed-9);
    gpt_sys_0100->SetMarkerStyle(20);
    gpt_sys_0100->SetMarkerSize(1);

    gpt_stat_0100->SetLineColor(MaterialColors[0]);
    gpt_stat_0100->SetMarkerStyle(20);
    gpt_stat_0100->SetMarkerSize(1);
    gpt_stat_0100->SetMarkerColor(MaterialColors[0]);
    gpt_sys_full_0100->SetMarkerColor(MaterialColors[0]);
    gpt_sys_full_0100->SetMarkerStyle(20);

    gpt_sys_full_0100->SetLineColor(MaterialColors[0]);
    gpt_sys_full_0100->SetFillColor(0);

    cMeanPt2->cd();
    //fMeanPt_sys_cor->Draw("a5");
    //fMeanPt_sys->Draw("5");
    //fMeanPt_stat->Draw("P");
    fMeanPt_sys_7TeV->Draw("A5");
    fMeanPt_stat_7TeV->Draw("P");;
    gpt_sys_0100->Draw("5");
    gpt_sys_full_0100->Draw("5");
    gpt_stat_0100->Draw("P");

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    auto lMeanPt2 = new TLegend(0.6, 0.2, 0.95, 0.42);
    lMeanPt2->SetFillStyle(0);
    lMeanPt2->AddEntry(fMeanPt_sys, "pp 13 TeV (INEL>0)", "PEF");
    lMeanPt2->AddEntry(gpt_sys_full_0100, "pp 13 TeV (INEL)", "PEF");
    lMeanPt2->AddEntry(fMeanPt_sys_7TeV, "pp 7 TeV", "PEF");
    lMeanPt2->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lMeanPt2->Draw();

    SaveCanvas(cMeanPt2,"figure_MeanPt_0100_preliminary","figs/Approval/");
    
    // pPb 5.02 TeV
    TCanvas* cMeanPt3 = GetCanvas("cMeanPt3", 960, 720);
    cMeanPt3->SetLeftMargin(0.10);
    cMeanPt3->Draw();
    // https://www.hepdata.net/record/ins1510878
    vector<double> xpPb502 = {35.6, 23.2, 16.1, 7.1};
    vector<double> xpPb502e = {0.8, 0.5, 0.4, 0.2};
    vector<double> xpPb502zero = {0, 0, 0, 0};
    
    vector<double> ptpPb502 = {1.626, 1.482, 1.459, 1.377};
    vector<double> ptpPb502e = {0.016, 0.02, 0.025, 0.023};
    vector<double> ptpPb502se = {0.068, 0.1, 0.114, 0.089};

    TGraphErrors* gpt_stat_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ptpPb502[0], &xpPb502zero[0], &ptpPb502e[0]);
    TGraphErrors* gpt_sys_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ptpPb502[0], &xpPb502e[0], &ptpPb502se[0]);

    gpt_sys_pPb502TeV->SetLineColor(MaterialColors[3]);
    gpt_sys_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    //gpt_sys_pPb502TeV->SetFillColorAlpha(MaterialColors[3], 0.3);
    gpt_stat_pPb502TeV->SetMarkerStyle(20);
    gpt_stat_pPb502TeV->SetMarkerSize(0.5);

    gpt_stat_pPb502TeV->SetLineColor(MaterialColors[3]);
    gpt_stat_pPb502TeV->SetMarkerColor(MaterialColors[3]);

    auto fMeanPt_sys_cor2 = (TGraphErrors*)fMeanPt_sys_cor->Clone();
    //fMeanPt_sys_cor2->SetMaximum(0.040);
    //fMeanPt_sys_cor2->GetYaxis()->SetTitleOffset(0.7);
    fMeanPt_sys_cor2->GetXaxis()->SetLimits(0,46);
    fMeanPt_sys_cor2->SetMinimum(meanpTYmin);

    cMeanPt3->cd();
    fMeanPt_sys_cor2->Draw("a5");
    gpt_sys_pPb502TeV->Draw("5");
    gpt_stat_pPb502TeV->Draw("P");
    fMeanPt_sys->Draw("5");
    fMeanPt_stat->Draw("P");
    if(draw7TeV) fMeanPt_sys_7TeV->Draw("A5");
    if(draw7TeV) fMeanPt_stat_7TeV->Draw("P");;
    //gpt_sys_0100->Draw("5");
    //gpt_sys_full_0100->Draw("5");
    //gpt_stat_0100->Draw("P");

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* lMeanPt3;
    if(draw7TeV) 
        lMeanPt3 = new TLegend(0.6, 0.26, 0.95, 0.48);
    else 
        lMeanPt3 = new TLegend(0.6, 0.26, 0.95, 0.42);
    lMeanPt3->SetFillStyle(0);
    lMeanPt3->AddEntry(fMeanPt_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lMeanPt3->AddEntry(fMeanPt_sys_7TeV, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) lMeanPt3->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lMeanPt3->AddEntry(gpt_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    lMeanPt3->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    lMeanPt3->Draw();

    SaveCanvas(cMeanPt3,"figure_MeanPt_withpPb502_preliminary","figs/Approval/");

    // PbPb Result
    TCanvas* cMeanPtPbPb = GetCanvas("cdNdy2", 960, 720);
    cMeanPtPbPb->SetLeftMargin(0.10);
    cMeanPtPbPb->SetBottomMargin(0.15);
    cMeanPtPbPb->SetLogx();
    cMeanPtPbPb->Draw();
    cMeanPtPbPb->cd();
    // charged particle multiplicity
    // https://arxiv.org/pdf/1012.1657.pdf
    vector<double> xPbPb276 = {1447.5, 680.3, 130.25};
    vector<double> xPbPb276e = {54.5, 25, 5.25};
    vector<double> xPbPb276zero = {0, 0, 0};
     
    // Analysis note PbPb 2.76
    vector<double> yPbPb276 = {1.748, 1.778, 1.515};
    vector<double> yPbPb276e = {0.015, 0.017, 0.016};
    vector<double> yPbPb276se = {0.250, 0.197, 0.156};
    
    TGraphErrors* gpt_stat_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276zero[0], &yPbPb276e[0]);
    TGraphErrors* gpt_sys_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276e[0], &yPbPb276se[0]);

    gpt_sys_PbPb276TeV->SetMarkerStyle(20);
    //ge_sys_PbPb276TeV->SetMarkerSize(0.5);
    gpt_stat_PbPb276TeV->SetMarkerStyle(20);
    gpt_stat_PbPb276TeV->SetMarkerSize(0.5);
    gpt_sys_PbPb276TeV->SetLineColor(MaterialColors[1]);
    gpt_sys_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    gpt_stat_PbPb276TeV->SetLineColor(MaterialColors[1]);
    gpt_stat_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    fMeanPt_sys_cor2->GetXaxis()->SetTitleOffset(1.3);
    fMeanPt_sys_cor2->GetXaxis()->SetLimits(1,5e3);
    fMeanPt_sys_cor2->SetMaximum(2.4);

    fMeanPt_sys_cor2->Draw("A5");
    fMeanPt_sys->Draw("5");
    fMeanPt_stat->Draw("P");
    gpt_sys_pPb502TeV->Draw("5");
    gpt_stat_pPb502TeV->Draw("P");
    if(draw7TeV) fMeanPt_sys_7TeV->Draw("5");
    if(draw7TeV) fMeanPt_stat_7TeV->Draw("P");;
    gpt_sys_PbPb276TeV->Draw("5");
    gpt_stat_PbPb276TeV->Draw("P");

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp #sqrt{#it{s}} = 13 TeV, INEL>0, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    auto lMeanPtPbPb = new TLegend(0.6, 0.20, 0.95, 0.5);
    lMeanPtPbPb->SetFillStyle(0);
    lMeanPtPbPb->AddEntry(fMeanPt_sys, "pp 13 TeV (INEL>0)", "PEF");
    lMeanPtPbPb->AddEntry(fMeanPt_sys_7TeV, "pp 7 TeV (INEL)", "PEF");
    lMeanPtPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lMeanPtPbPb->AddEntry(gpt_sys_pPb502TeV, "p-Pb 5.02TeV", "PEF");
    lMeanPtPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    lMeanPtPbPb->AddEntry(gpt_sys_PbPb276TeV, "Pb-Pb 2.76TeV", "PEF");
    lMeanPtPbPb->AddEntry((TObject*)0, "ALICE prelimiary", "");
    lMeanPtPbPb->Draw();

    SaveCanvas(cMeanPtPbPb,"figure_MeanPt_withpPb502_PbPb_preliminary","figs/Approval/");
}
void DrawRatioToPi(){
    double ratioPiYmax = 2.5e-3;
    double ratioPiYmin = 0.4e-3;
    bool draw7TeV = false;
    TCanvas* cPiRatio = GetCanvas("cPiRatio", 960, 720);
    cPiRatio->SetLeftMargin(0.10);
    cPiRatio->Draw();
    
    TFile* fPiRatio = new TFile(PhysicsFile.Data());
    
    auto gRatioPi_stat = (TGraphErrors*)fPiRatio->Get("gRatioPi_stat");  // signal
    auto gRatioPi_sys_cor = (TGraphErrors*)fPiRatio->Get("gRatioPi_sys");  // signal
    auto gRatioPi_sys = (TGraphErrors*)fPiRatio->Get("gRatioPi_sys_full");  // signal
    auto gRatioPi_7TeV_stat = (TGraphErrors*)fPiRatio->Get("gRatioPi_7TeV_stat");  // signal
    auto gRatioPi_7TeV_sys = (TGraphErrors*)fPiRatio->Get("gRatioPi_7TeV_sys");  // signal

    gRatioPi_sys->SetMarkerColor(MaterialColors[0]);
    gRatioPi_sys->SetLineColor(MaterialColors[0]);
    gRatioPi_sys->SetFillColorAlpha(MaterialColors[0], 0.0);

    gRatioPi_sys_cor->SetMarkerColor(MaterialColors[0]);
    gRatioPi_sys_cor->SetLineColor(0);
    gRatioPi_sys_cor->SetFillColorAlpha(MaterialColors[0], 0.3);
    gRatioPi_sys_cor->GetYaxis()->SetLabelSize(0.05);
    gRatioPi_sys_cor->GetYaxis()->SetTitleOffset(0.9);
    gRatioPi_sys_cor->GetYaxis()->SetTitleSize(0.05);
    gRatioPi_sys_cor->GetXaxis()->SetLabelSize(0.05);
    gRatioPi_sys_cor->GetXaxis()->SetTitleSize(0.05);
    // fMeanPt_sys->SetMaximum(0.008);
    gRatioPi_sys_cor->SetMaximum(ratioPiYmax);
    gRatioPi_sys_cor->SetMinimum(ratioPiYmin);
    gRatioPi_sys->SetMarkerStyle(20);
    gRatioPi_sys_cor->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    gRatioPi_sys_cor->GetYaxis()->SetTitle("Ratio to #pi");

    gRatioPi_stat->SetMarkerColor(MaterialColors[0]);
    gRatioPi_stat->SetLineColor(MaterialColors[0]);

    gRatioPi_7TeV_sys->SetMarkerColor(MaterialColors[0]);
    gRatioPi_7TeV_sys->SetLineColor(0);
    gRatioPi_7TeV_sys->SetFillColorAlpha(MaterialColors[0], 0.3);
    gRatioPi_7TeV_sys->GetYaxis()->SetLabelSize(0.05);
    gRatioPi_7TeV_sys->GetYaxis()->SetTitleOffset(0.7);
    gRatioPi_7TeV_sys->GetYaxis()->SetTitleSize(0.05);
    gRatioPi_7TeV_sys->GetXaxis()->SetLabelSize(0.05);
    gRatioPi_7TeV_sys->GetXaxis()->SetTitleSize(0.05);
    // fMeanPt_sys->SetMaximum(0.008);
    gRatioPi_7TeV_sys->SetMaximum(2.0e-3);
    gRatioPi_7TeV_sys->SetMinimum(0.0);
    gRatioPi_7TeV_sys->SetMarkerStyle(20);
    gRatioPi_7TeV_sys->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    gRatioPi_7TeV_sys->GetYaxis()->SetTitle("Ratio to #pi");

    gRatioPi_7TeV_stat->SetMarkerColor(kBlack);
    gRatioPi_7TeV_stat->SetLineColor(kBlack);

    gRatioPi_7TeV_sys->SetMarkerColor(kBlack);
    gRatioPi_7TeV_sys->SetLineColor(kBlack);
    gRatioPi_7TeV_sys->SetMarkerStyle(20);
    gRatioPi_7TeV_sys->SetFillColor(0);

    cPiRatio->cd();
    gRatioPi_sys_cor->Draw("a5");
    gRatioPi_sys->Draw("5");
    gRatioPi_stat->Draw("P");
    if(draw7TeV) gRatioPi_7TeV_sys->Draw("5");
    if(draw7TeV) gRatioPi_7TeV_stat->Draw("P");

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* lPiRatio;
    if(draw7TeV)
        lPiRatio = new TLegend(0.55, 0.15, 0.90, 0.31);
    else
        lPiRatio = new TLegend(0.55, 0.2, 0.90, 0.25);
    lPiRatio->SetFillStyle(0);
    lPiRatio->AddEntry(gRatioPi_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lPiRatio->AddEntry(gRatioPi_7TeV_sys, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) lPiRatio->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lPiRatio->Draw();

    SaveCanvas(cPiRatio,"figure_RatioToPi_preliminary","figs/Approval/");


    // Further figure
    TCanvas* cPiRatio2 = GetCanvas("cPiRatio2", 960, 720);
    cPiRatio2->SetLeftMargin(0.10);
    cPiRatio2->Draw();
    // pPb Result
    // https://www.hepdata.net/record/ins1510878
    vector<double> xpPb502 = {35.6, 23.2, 16.1, 7.1};
    vector<double> xpPb502e = {0.8, 0.5, 0.4, 0.2};
    vector<double> xpPb502zero = {0, 0, 0, 0};
    
    // HEP data Xi(1530) pPb5.02TeV
    vector<double> ypPb502 = { 0.00168, 0.00163, 0.0014,0.00105};
    vector<double> ypPb502e = { 4.0e-05, 5.0e-05, 4.0e-05,3.0e-05};
    vector<double> ypPb502se = { 0.00019, 0.00023, 0.00022,0.00015};
    vector<double> ypPb502se_uncor = { 0.00009, 0.00012, 0.00011,0.00008};
    vector<double> ypPb502se_full;
    for(int i = 0; i < ypPb502.size(); i++){
        double temp = ypPb502[i]*sqrt( pow(ypPb502se[i]/ypPb502[i],2) + pow(ypPb502se_uncor[i]/ypPb502[i],2));
        ypPb502se_full.push_back(temp);
    }
    
    TGraphErrors* gRatioPi_stat_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502zero[0], &ypPb502e[0]);
    TGraphErrors* gRatioPi_sys_cor_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502e[0], &ypPb502se[0]);
    TGraphErrors* gRatioPi_sys_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502e[0], &ypPb502se_full[0]);
    
    //MB
    TFile* fdNdy0100 = new TFile("AnalysisResults_Xi1530_YieldMean_0-100.root");
    auto fdNdy_YieldMean0100 = (TH1D*)fdNdy0100->Get("hYield_0.00-100.00");  // signal
    vector<double> yield0100;
    yield0100.push_back(fdNdy_YieldMean0100->GetBinContent(1));
    vector<double> yield0100state;
    yield0100state.push_back(fdNdy_YieldMean0100->GetBinContent(2));
    vector<double> yield0100syse;
    double ftemp_yieldsyse = 0.0;
    ftemp_yieldsyse += pow(fdNdy_YieldMean0100->GetBinContent(3),2);
    ftemp_yieldsyse += pow(fdNdy_YieldMean0100->GetBinContent(4),2);
    yield0100syse.push_back(sqrt(ftemp_yieldsyse)/2);
    vector<double> yield0100syse_full;
    ftemp_yieldsyse += pow(0.0538*yield0100[0],2);
    yield0100syse_full.push_back(sqrt(ftemp_yieldsyse)/2);

    vector<double> xppMB = {6.94};
    vector<double> xppMBe = {0.10};
    vector<double> xppMBzero = {0};

    double piMB = 6.36830;
    double piMBsys = 0.45657;

    vector<double> yRatio;
    vector<double> yRatioe;
    vector<double> yRatiosyse;
    yRatio.push_back(yield0100[0]/(piMB/2));
    yRatioe.push_back(yRatio[0]*(yield0100state[0]/yield0100[0]));
    yRatiosyse.push_back(yRatio[0]*sqrt(pow(piMBsys/piMB,2)+pow(yield0100syse_full[0]/yield0100[0],2)));
    cout << yRatio[0] << endl;
    TGraphErrors* ge_stat_MB = new TGraphErrors(yRatio.size(), &xppMB[0], &yRatio[0], &xppMBzero[0], &yRatioe[0]);
    TGraphErrors* ge_sys_MB = new TGraphErrors(yRatio.size(), &xppMB[0], &yRatio[0], &xppMBe[0], &yRatiosyse[0]);
    ge_stat_MB->SetLineColor(MaterialColors[0]);
    ge_stat_MB->SetMarkerColor(MaterialColors[0]);
    ge_sys_MB->SetLineColor(MaterialColors[0]);
    ge_sys_MB->SetMarkerColor(MaterialColors[0]);
    ge_sys_MB->SetFillColorAlpha(MaterialColors[3], 0.3);
    ge_stat_MB->SetMarkerStyle(20);
    ge_stat_MB->SetMarkerSize(0.5);

    //



    gRatioPi_sys_cor_pPb502TeV->SetLineColor(0);
    gRatioPi_sys_cor_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    gRatioPi_sys_cor_pPb502TeV->SetFillColorAlpha(MaterialColors[3], 0.3);
    gRatioPi_sys_cor_pPb502TeV->SetMarkerStyle(20);
    gRatioPi_sys_cor_pPb502TeV->SetMarkerSize(0.5);

    gRatioPi_sys_pPb502TeV->SetLineColor(MaterialColors[3]);
    gRatioPi_sys_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    gRatioPi_sys_pPb502TeV->SetFillColorAlpha(MaterialColors[3], 0.0);
    gRatioPi_sys_pPb502TeV->SetMarkerStyle(20);
    //gRatioPi_sys_pPb502TeV->SetMarkerSize(0.5);

    gRatioPi_stat_pPb502TeV->SetLineColor(MaterialColors[3]);
    gRatioPi_stat_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    gRatioPi_stat_pPb502TeV->SetMarkerStyle(20);
    gRatioPi_stat_pPb502TeV->SetMarkerSize(0.5);
    
    auto gRatioPi_sys_cor2 = (TGraphErrors*)gRatioPi_sys_cor->Clone();
    //gRatioPi_sys_cor2->GetXaxis()->SetLimits(1,5e3);
    gRatioPi_sys_cor2->GetXaxis()->SetLimits(1,46);

    cPiRatio2->cd();

    gRatioPi_sys_cor2->Draw("A5");
    gRatioPi_sys->Draw("5");
    gRatioPi_stat->Draw("P");

    if(draw7TeV) gRatioPi_7TeV_sys->Draw("5");
    if(draw7TeV) gRatioPi_7TeV_stat->Draw("P");

    gRatioPi_sys_cor_pPb502TeV->Draw("5");
    gRatioPi_sys_pPb502TeV->Draw("5");
    gRatioPi_stat_pPb502TeV->Draw("P");



    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* lPiRatio2;
    if(draw7TeV)
        lPiRatio2 = new TLegend(0.55, 0.15, 0.90, 0.45);
    else
        lPiRatio2 = new TLegend(0.55, 0.15, 0.90, 0.35);
    lPiRatio2->SetFillStyle(0);
    lPiRatio2->AddEntry(gRatioPi_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lPiRatio2->AddEntry(gRatioPi_7TeV_sys, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) lPiRatio2->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lPiRatio2->AddEntry(gRatioPi_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    lPiRatio2->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    lPiRatio2->Draw();

    SaveCanvas(cPiRatio2,"figure_RatioToPi_withpPb502_preliminary","figs/Approval/");

    // PbPb Result
    TCanvas* cPiRatioPbPb = GetCanvas("cdNdy2", 960, 720);
    cPiRatioPbPb->SetLeftMargin(0.10);
    cPiRatioPbPb->SetBottomMargin(0.15);
    cPiRatioPbPb->Draw();
    cPiRatioPbPb->cd();
    // charged particle multiplicity
    // https://arxiv.org/pdf/1012.1657.pdf
    vector<double> xPbPb276 = {1447.5, 680.3, 130.25};
    vector<double> xPbPb276e = {54.5, 25, 5.25};
    vector<double> xPbPb276zero = {0, 0, 0};
     
    // Analysis note PbPb 2.76
    vector<double> yPbPb276 = {1.03e-3, 1.06e-3, 1.36e-3};
    vector<double> yPbPb276e = {0.02678e-3, 0.018e-3, 0.0408e-3};
    vector<double> yPbPb276se = {0.35e-3, 0.26e-3, 0.36e-3};
    
    TGraphErrors* ge_stat_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276zero[0], &yPbPb276e[0]);
    TGraphErrors* ge_sys_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276e[0], &yPbPb276se[0]);

    ge_sys_PbPb276TeV->SetMarkerStyle(20);
    ge_stat_PbPb276TeV->SetMarkerStyle(20);
    ge_stat_PbPb276TeV->SetMarkerSize(0.5);
    ge_sys_PbPb276TeV->SetLineColor(MaterialColors[1]);
    ge_sys_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    ge_stat_PbPb276TeV->SetLineColor(MaterialColors[1]);
    ge_stat_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    gRatioPi_sys_cor2->GetXaxis()->SetLimits(1,5e3);
    gRatioPi_sys_cor2->SetMaximum(3e-3);
    gRatioPi_sys_cor2->SetMinimum(0.5e-3);

    gRatioPi_sys_cor2->GetYaxis()->SetTitleOffset(0.9);
    gRatioPi_sys_cor2->GetXaxis()->SetTitleOffset(1.3);

    gRatioPi_sys_cor2->Draw("A5");
    gRatioPi_sys->Draw("5");
    gRatioPi_stat->Draw("P");

    if(draw7TeV) gRatioPi_7TeV_sys->Draw("A5");
    if(draw7TeV) gRatioPi_7TeV_stat->Draw("P");

    gRatioPi_sys_cor_pPb502TeV->Draw("5");
    gRatioPi_sys_pPb502TeV->Draw("5");
    gRatioPi_stat_pPb502TeV->Draw("P");
    
    //ge_sys_MB->Draw("5");
    //ge_stat_MB->Draw("P");
    ge_sys_PbPb276TeV->Draw("5");
    ge_stat_PbPb276TeV->Draw("P");

    cPiRatioPbPb->SetLogx();
    //cPiRatioPbPb->SetLogy();
    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");


    TLegend* lPiRatioPbPb;
    if(draw7TeV)
        lPiRatioPbPb = new TLegend(0.6, 0.57, 0.95, 0.92);
    else
        lPiRatioPbPb = new TLegend(0.6, 0.67, 0.95, 0.92);
    lPiRatioPbPb->SetFillStyle(0);
    lPiRatioPbPb->AddEntry(gRatioPi_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lPiRatioPbPb->AddEntry(gRatioPi_7TeV_sys, "pp 7 TeV (INEL>0)", "PEF");
    if(draw7TeV) lPiRatioPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lPiRatioPbPb->AddEntry(gRatioPi_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    lPiRatioPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    lPiRatioPbPb->AddEntry(ge_sys_PbPb276TeV, "Pb-Pb 2.76TeV (INEL)", "PEF");
    lPiRatioPbPb->AddEntry((TObject*)0, "ALICE prelimiary", "");
    lPiRatioPbPb->Draw();

    SaveCanvas(cPiRatioPbPb,"figure_RatioToPi_withpPb502_PbPb_preliminary","figs/Approval/");
    
}
void DrawRatioToXi(){
    double ratioXiYmax = 0.55;
    double ratioXiYmin = 0.1;
    bool draw7TeV = false;

    TCanvas* cXiRatio = GetCanvas("cXiRatio", 960, 720);
    cXiRatio->SetLeftMargin(0.10);
    cXiRatio->Draw();
    
    TFile* fXiRatio = new TFile(PhysicsFile.Data());
    
    auto gRatioToXi_stat = (TGraphErrors*)fXiRatio->Get("gRatioToXi_stat");  // signal
    auto gRatioToXi_sys_cor = (TGraphErrors*)fXiRatio->Get("gRatioToXi_sys_cor");  // signal
    auto gRatioToXi_sys = (TGraphErrors*)fXiRatio->Get("gRatioToXi_sys");  // signal
    auto gRatioToXi_7TeV_stat = (TGraphErrors*)fXiRatio->Get("gRatioToXi_7TeV_stat");  // signal
    auto gRatioToXi_7TeV_sys = (TGraphErrors*)fXiRatio->Get("gRatioToXi_7TeV_sys");  // signal

    //
    //MB
    TFile* fdNdy0100 = new TFile("AnalysisResults_Xi1530_YieldMean_0-100.root");
    auto fdNdy_YieldMean0100 = (TH1D*)fdNdy0100->Get("hYield_0.00-100.00");  // signal
    vector<double> yield0100;
    yield0100.push_back(fdNdy_YieldMean0100->GetBinContent(1));
    vector<double> yield0100state;
    yield0100state.push_back(fdNdy_YieldMean0100->GetBinContent(2));
    vector<double> yield0100syse;
    double ftemp_yieldsyse = 0.0;
    ftemp_yieldsyse += pow(fdNdy_YieldMean0100->GetBinContent(3),2);
    ftemp_yieldsyse += pow(fdNdy_YieldMean0100->GetBinContent(4),2);
    yield0100syse.push_back(sqrt(ftemp_yieldsyse)/2);
    vector<double> yield0100syse_full;
    ftemp_yieldsyse += pow(0.0538*yield0100[0],2);
    yield0100syse_full.push_back(sqrt(ftemp_yieldsyse)/2);

    vector<double> xppMB = {6.94};
    vector<double> xppMBe = {0.10};
    vector<double> xppMBzero = {0};

    double XiMB = 0.0133295;
    double XiMBe = 0.0001160;
    double XiMBsys = 0.0005296;

    vector<double> yRatio;
    vector<double> yRatioe;
    vector<double> yRatiosyse;
    yRatio.push_back(yield0100[0]/(XiMB));
    yRatioe.push_back(yRatio[0]*sqrt(pow(XiMBe/XiMB,2)+pow(yield0100state[0]/yield0100[0],2)));
    yRatiosyse.push_back(yRatio[0]*sqrt(pow(XiMBsys/XiMB,2)+pow(yield0100syse_full[0]/yield0100[0],2)));
    cout << yRatio[0] << endl;
    TGraphErrors* ge_stat_MB = new TGraphErrors(yRatio.size(), &xppMB[0], &yRatio[0], &xppMBzero[0], &yRatioe[0]);
    TGraphErrors* ge_sys_MB = new TGraphErrors(yRatio.size(), &xppMB[0], &yRatio[0], &xppMBe[0], &yRatiosyse[0]);
    ge_stat_MB->SetLineColor(MaterialColors[0]);
    ge_stat_MB->SetMarkerColor(MaterialColors[0]);
    ge_sys_MB->SetLineColor(MaterialColors[0]);
    ge_sys_MB->SetMarkerColor(MaterialColors[0]);
    ge_sys_MB->SetFillColorAlpha(MaterialColors[0], 0.3);
    ge_stat_MB->SetMarkerStyle(20);
    ge_stat_MB->SetMarkerSize(0.5);

    //
    //  


    gRatioToXi_sys->SetMarkerColor(MaterialColors[0]);
    gRatioToXi_sys->SetLineColor(MaterialColors[0]);
    gRatioToXi_sys->SetFillColorAlpha(MaterialColors[0], 0.0);

    gRatioToXi_sys_cor->SetMarkerColor(MaterialColors[0]);
    gRatioToXi_sys_cor->SetLineColor(0);
    gRatioToXi_sys_cor->SetFillColorAlpha(MaterialColors[0], 0.3);
    gRatioToXi_sys_cor->SetFillStyle(1000);
    gRatioToXi_sys_cor->GetYaxis()->SetLabelSize(0.04);
    gRatioToXi_sys_cor->GetYaxis()->SetTitleOffset(1.0);
    gRatioToXi_sys_cor->GetYaxis()->SetTitleSize(0.05);
    gRatioToXi_sys_cor->GetXaxis()->SetTitleOffset(1.0);
    gRatioToXi_sys_cor->GetXaxis()->SetTitleSize(0.05);
    // fMeanPt_sys->SetMaximum(0.008);
    gRatioToXi_sys_cor->SetMaximum(ratioXiYmax);
    gRatioToXi_sys_cor->SetMinimum(ratioXiYmin);
    gRatioToXi_sys->SetMarkerStyle(20);
    gRatioToXi_sys_cor->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    gRatioToXi_sys_cor->GetYaxis()->SetTitle("Ratio to #Xi");

    gRatioToXi_stat->SetMarkerColor(MaterialColors[0]);
    gRatioToXi_stat->SetLineColor(MaterialColors[0]);
    gRatioToXi_stat->SetMarkerStyle(20);
    gRatioToXi_stat->SetMarkerSize(0.5);

    gRatioToXi_7TeV_sys->GetYaxis()->SetLabelSize(0.05);
    gRatioToXi_7TeV_sys->GetYaxis()->SetTitleOffset(1.2);
    gRatioToXi_7TeV_sys->GetYaxis()->SetTitleSize(0.05);
    gRatioToXi_7TeV_sys->GetXaxis()->SetLabelSize(0.05);
    gRatioToXi_7TeV_sys->GetXaxis()->SetTitleOffset(1.3);
    gRatioToXi_7TeV_sys->GetXaxis()->SetTitleSize(0.07);
    // fMeanPt_sys->SetMaximum(0.008);
    gRatioToXi_7TeV_sys->SetMarkerStyle(20);
    gRatioToXi_7TeV_sys->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{y}|<0.5}");
    gRatioToXi_7TeV_sys->GetYaxis()->SetTitle("Ratio to #Xi");
    gRatioToXi_7TeV_stat->SetMarkerColor(kBlack);
    gRatioToXi_7TeV_stat->SetLineColor(kBlack);

    gRatioToXi_7TeV_sys->SetMarkerColor(kBlack);
    gRatioToXi_7TeV_sys->SetLineColor(kBlack);
    gRatioToXi_7TeV_sys->SetMarkerStyle(20);
    gRatioToXi_7TeV_sys->SetFillColor(0);

    gRatioToXi_sys_cor->GetXaxis()->SetLimits(0,25);

    cXiRatio->cd();
    gRatioToXi_sys_cor->Draw("A5");
    gRatioToXi_sys->Draw("5");
    gRatioToXi_stat->Draw("P");
    if(draw7TeV) gRatioToXi_7TeV_sys->Draw("5");
    if(draw7TeV) gRatioToXi_7TeV_stat->Draw("P");
    
    //ge_sys_MB->Draw("5");
    //ge_stat_MB->Draw("P");

    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");

    TLegend* lXiRatio;
    if(draw7TeV)
        lXiRatio = new TLegend(0.55, 0.15, 0.90, 0.31);
    else
        lXiRatio = new TLegend(0.55, 0.20, 0.90, 0.25);
    lXiRatio->SetFillStyle(0);
    lXiRatio->AddEntry(gRatioToXi_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lXiRatio->AddEntry(gRatioToXi_7TeV_sys, "pp 7 TeV", "PEF");
    if(draw7TeV) lXiRatio->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lXiRatio->Draw();

    SaveCanvas(cXiRatio,"figure_RatioToXi_preliminary","figs/Approval/");


    // Further figure
    TCanvas* cXiRatio2 = GetCanvas("cXiRatio2", 960, 720);
    cXiRatio2->SetLeftMargin(0.10);
    cXiRatio2->Draw();
    // pPb Result
    // https://www.hepdata.net/record/ins1510878
    vector<double> xpPb502 = {35.6, 23.2, 16.1, 7.1};
    vector<double> xpPb502e = {0.8, 0.5, 0.4, 0.2};
    vector<double> xpPb502zero = {0, 0, 0, 0};
    
    // HEP data Xi(1530) pPb5.02TeV
    vector<double> ypPb502 = { 0.301, 0.317, 0.288, 0.259};
    vector<double> ypPb502e = { 0.007, 0.009, 0.009,0.008};
    vector<double> ypPb502se = { 0.038, 0.048, 0.044, 0.042};
    vector<double> ypPb502se_uncor = { 0.019, 0.024, 0.022,0.021};
    vector<double> ypPb502se_full;
    for(int i = 0; i < ypPb502.size(); i++){
        double temp = ypPb502[i]*sqrt( pow(ypPb502se[i]/ypPb502[i],2) + pow(ypPb502se_uncor[i]/ypPb502[i],2));
        ypPb502se_full.push_back(temp);
    }
    
    TGraphErrors* gRatioToXi_stat_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502zero[0], &ypPb502e[0]);
    TGraphErrors* gRatioToXi_sys_cor_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502e[0], &ypPb502se[0]);
    TGraphErrors* gRatioToXi_sys_pPb502TeV = new TGraphErrors(xpPb502.size(), &xpPb502[0], &ypPb502[0], &xpPb502e[0], &ypPb502se_full[0]);

    gRatioToXi_sys_cor_pPb502TeV->SetLineColor(0);
    gRatioToXi_sys_cor_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    gRatioToXi_sys_cor_pPb502TeV->SetFillColorAlpha(MaterialColors[3], 0.3);
    gRatioToXi_sys_cor_pPb502TeV->SetMarkerStyle(20);
    gRatioToXi_sys_cor_pPb502TeV->SetMarkerSize(0.5);

    gRatioToXi_sys_pPb502TeV->SetLineColor(MaterialColors[3]);
    gRatioToXi_sys_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    gRatioToXi_sys_pPb502TeV->SetFillColorAlpha(MaterialColors[3], 0.0);
    gRatioToXi_sys_pPb502TeV->SetMarkerStyle(20);
    //gRatioToXi_sys_pPb502TeV->SetMarkerSize(0.5);

    gRatioToXi_stat_pPb502TeV->SetLineColor(MaterialColors[3]);
    gRatioToXi_stat_pPb502TeV->SetMarkerColor(MaterialColors[3]);
    gRatioToXi_stat_pPb502TeV->SetMarkerStyle(20);
    gRatioToXi_stat_pPb502TeV->SetMarkerSize(0.5);
    
    auto gRatioToXi_sys_cor2 = (TGraphErrors*)gRatioToXi_sys_cor->Clone();
    gRatioToXi_sys_cor2->GetXaxis()->SetLimits(0,46);
    //gRatioToXi_sys_cor2->SetMaximum(0.7);
    //cXiRatio2->SetLogx();
    cXiRatio2->cd();

    gRatioToXi_sys_cor2->Draw("A5");
    gRatioToXi_sys->Draw("5");
    gRatioToXi_stat->Draw("P");

    gRatioToXi_sys_cor_pPb502TeV->Draw("5");
    gRatioToXi_sys_pPb502TeV->Draw("5");
    gRatioToXi_stat_pPb502TeV->Draw("P");
    


    if(draw7TeV) gRatioToXi_7TeV_sys->Draw("5");
    if(draw7TeV) gRatioToXi_7TeV_stat->Draw("P");

    //ge_sys_MB->Draw("5");
    //ge_stat_MB->Draw("P");

    
    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");
    TLegend* lXiRatio2;
    if(draw7TeV)
        lXiRatio2 = new TLegend(0.6, 0.15, 0.95, 0.37);
    else
        lXiRatio2 = new TLegend(0.6, 0.15, 0.95, 0.27);
    lXiRatio2->SetFillStyle(0);
    lXiRatio2->AddEntry(gRatioToXi_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lXiRatio2->AddEntry(gRatioToXi_7TeV_sys, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) lXiRatio2->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lXiRatio2->AddEntry(gRatioToXi_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    lXiRatio2->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    lXiRatio2->Draw();

    SaveCanvas(cXiRatio2,"figure_RatioToXi_withpPb502_preliminary","figs/Approval/");

    // PbPb Result
    TCanvas* cXiRatioPbPb = GetCanvas("cdNdy2", 960, 720);
    cXiRatioPbPb->SetLeftMargin(0.1);
    cXiRatioPbPb->Draw();
    cXiRatioPbPb->cd();
    // charged particle multiplicity
    // https://arxiv.org/pdf/1012.1657.pdf
    vector<double> xPbPb276 = {1447.5, 680.3, 130.25};
    vector<double> xPbPb276e = {54.5, 25, 5.25};
    vector<double> xPbPb276zero = {0, 0, 0};
     
    // Analysis note PbPb 2.76
    vector<double> yPbPb276 = {0.206, 0.166, 0.25};
    vector<double> yPbPb276e = {0.006, 0.004, 0.009};
    vector<double> yPbPb276se = {0.070, 0.041, 0.074};
    
    TGraphErrors* ge_stat_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276zero[0], &yPbPb276e[0]);
    TGraphErrors* ge_sys_PbPb276TeV = new TGraphErrors(xPbPb276.size(), &xPbPb276[0], &yPbPb276[0], &xPbPb276e[0], &yPbPb276se[0]);

    ge_sys_PbPb276TeV->SetMarkerStyle(20);
    ge_stat_PbPb276TeV->SetMarkerStyle(20);
    ge_stat_PbPb276TeV->SetMarkerSize(0.5);
    ge_sys_PbPb276TeV->SetLineColor(MaterialColors[1]);
    ge_sys_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    ge_stat_PbPb276TeV->SetLineColor(MaterialColors[1]);
    ge_stat_PbPb276TeV->SetMarkerColor(MaterialColors[1]);
    gRatioToXi_sys_cor2->GetXaxis()->SetLimits(1,5e3);
    //gRatioPi_7TeV_sys->SetMaximum(4e-3);
    //gRatioPi_7TeV_sys->SetMinimum(0);
    gRatioToXi_sys_cor2->GetYaxis()->SetLabelSize(0.04);
    gRatioToXi_sys_cor2->GetYaxis()->SetTitleOffset(1.0);
    gRatioToXi_sys_cor2->GetYaxis()->SetTitleSize(0.05);


    gRatioToXi_sys_cor2->Draw("A5");
    gRatioToXi_sys->Draw("5");
    gRatioToXi_stat->Draw("P");
    gRatioToXi_sys_cor_pPb502TeV->Draw("5");
    gRatioToXi_sys_pPb502TeV->Draw("5");
    gRatioToXi_stat_pPb502TeV->Draw("P");
    //ge_sys_MB->Draw("5");
    //ge_stat_MB->Draw("P");
    ge_sys_PbPb276TeV->Draw("5");
    ge_stat_PbPb276TeV->Draw("P");  
    if(draw7TeV) gRatioToXi_7TeV_sys->Draw("5");
    if(draw7TeV) gRatioToXi_7TeV_stat->Draw("P");
    cXiRatioPbPb->SetLogx();
    
    t2->DrawLatex(0.15, 0.88, "ALICE Preliminary");
    t->DrawLatex(0.15, 0.82,
                 "#bf{pp (INEL>0) #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5}");
    t2->DrawLatex(0.15, 0.76,
                  "#bf{#Xi(1530)^{0} + #bar{#Xi}(1530)^{0} "
                  "#rightarrow #Xi^{-}#pi^{+} + #bar{#Xi}^{+}#pi^{-}}");

    t0->DrawLatex(0.15, 0.71, "#bf{Uncertainties: stat.(bars), syst.(boxes)}");
    TLegend* lXiRatioPbPb;
    if(draw7TeV)
        lXiRatioPbPb = new TLegend(0.6, 0.57, 0.95, 0.92);
    else
        lXiRatioPbPb = new TLegend(0.6, 0.67, 0.95, 0.92);
    lXiRatioPbPb->SetFillStyle(0);
    lXiRatioPbPb->AddEntry(gRatioToXi_sys, "pp 13 TeV (INEL>0)", "PEF");
    if(draw7TeV) lXiRatioPbPb->AddEntry(gRatioToXi_7TeV_sys, "pp 7 TeV (INEL)", "PEF");
    if(draw7TeV) lXiRatioPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 75 (2015) 1", "");
    lXiRatioPbPb->AddEntry(gRatioToXi_sys_pPb502TeV, "p-Pb 5.02TeV (NSD)", "PEF");
    lXiRatioPbPb->AddEntry((TObject*)0, "Eur. Phys. J. C 77 (2017) 389", "");
    lXiRatioPbPb->AddEntry(ge_sys_PbPb276TeV, "Pb-Pb 2.76TeV (INEL)", "PEF");
    lXiRatioPbPb->AddEntry((TObject*)0, "ALICE prelimiary", "");
    lXiRatioPbPb->Draw();

    SaveCanvas(cXiRatioPbPb,"figure_RatioToXi_withpPb502_PbPb_preliminary","figs/Approval/");
}
TCanvas* GetCanvas(TString cname, double w, double h){
    TCanvas* c = new TCanvas(Form("%s%d",cname.Data(),canvasCount), Form("%s%d",cname.Data(),canvasCount), w, h);
    c->SetTickx();
    c->SetTicky();
    c->SetTopMargin(0.05);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.01);
    c->SetFillStyle(0);
    canvasCount++;
    return c;
}