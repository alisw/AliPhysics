#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TString.h>
#include <AliPWGFunc.h>
#include "DrawingHelper.C"
// Save?
bool savefile = true;
bool corey = false;
TString savetype = "pdf";
TString savepath = "./data";

//functions
double myLevyPt(Double_t* x, Double_t* par);
void PlotXi1530INEL(TString finputfile, char const* options = "SAVE") {
    // Default Canvas
    Double_t w = 1920;
    Double_t h = 1080;

    //Default draw range
    //vector<double> DrawRange = {1.484, 1.592};
    vector<double> DrawRange = {1.484, 1.8};
    string ftemp = finputfile.Data();

    savepath = finputfile(finputfile.Index("INEL_") - 2,
                   finputfile.Index("root") - finputfile.Index("INEL_") + 1);

    double multi_start =0;
    double multi_end = 100;
    cout << "Multi bin: " << multi_start << "-" << multi_end << endl;
    
    TFile* inputfile = new TFile(finputfile.Data());
    TString Options = options;
    if (Options.Contains("SAVE"))
        savefile = true;
    if (Options.Contains("PNG"))
        savetype = "png";
    if (finputfile.Contains("Corey"))
        corey = true;
    //=========Basic QA Plots=============
    TCanvas* cQA = new TCanvas("cQA", "", w, h);
    cQA->Draw();
    cQA->SetTickx();
    cQA->SetLogy(false);
    gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(0);
    cQA->cd();
    // N of events
    auto hNumberofEvent = (TH1D*)inputfile->Get("hNumberofEvent");
    hNumberofEvent->Draw("HIST text");
    SaveCanvas(cQA, "hNumberofEvent", Form("%s/QA/", savepath.Data()),
               savetype);
    
    // Trigger Effi
    auto TriggerEffi = (TH1D*)inputfile->Get("hTriggerEffi");
    TriggerEffi->SetTitle(Form("Int7 Trigger Efficiency in %.2f-%.2f bin",
                               multi_start, multi_end));
    TriggerEffi->GetXaxis()->SetTitle("Charged particle multiplicty percentile(%)");
    TriggerEffi->Draw("E TEXT");
    SaveCanvas(cQA, "hTriggerEffi", Form("%s/QA/", savepath.Data()),
               savetype);

    //=========Data Plots=============
    vector<double> ptbin = {0, 0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6, 8.8, 15};
    if (corey)
        ptbin = {0, 20};
    double* ptbin_array = &ptbin[0]; // for histogram axis
    Int_t colors[] = {633, 810, 800, 830, 840, 840, 870, 864, 890, 617, 619};

    // for memo, small
    TLatex* t = new TLatex();
    t->SetNDC();
    t->SetTextSize(0.035);
    // for warning, big
    TLatex* t2 = new TLatex();
    t2->SetNDC();
    t2->SetTextSize(0.04);
    // for small pad, huge
    TLatex* t3 = new TLatex();
    t3->SetNDC();
    t3->SetTextSize(0.09);

    // Signal to Background
    TCanvas* cSigBkg = new TCanvas("cSigBkg", "cSigBkg", w, h);
    cSigBkg->Draw();
    cSigBkg->Divide(4, 2, 0.0001, 0.0001);
    cSigBkg->SetTickx();
    cSigBkg->SetLogy(false);
    cSigBkg->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);

    for (int bin = 0; bin < ptbin.size()-1; bin++) {  // ptbin = 10
        if ((!corey) && ((bin == 0) || (bin == 10)))
                continue;
        cout << "ptBin: " << bin<< endl;
        cSigBkg->cd(bin);
        gPad->SetLeftMargin(0.10);
        TGaxis::SetMaxDigits(3);

        auto hSig = (TH1D*)inputfile->Get(Form("hSignalOnly_%i", bin));
        auto hBkg = (TH1D*)inputfile->Get(Form("hBkgOnly_%i", bin));
        auto hBkgNorm = (TH1D*)inputfile->Get(Form("hBkgNorm_%i", bin));
        hBkgNorm->SetFillColorAlpha(kRed, 0.1);
        hSig->GetYaxis()->SetTitleOffset(1.2);
        hSig->GetYaxis()->SetTitleSize(0.04);
        hSig->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hBkg->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hBkgNorm->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hSig->Draw("PZ");
        hBkg->Draw("PZsame");
        hBkgNorm->Draw("BAR same");

        // Legend
        auto legend = new TLegend(0.13, 0.74, 0.5, 0.89);
        legend->SetFillStyle(0);
        legend->AddEntry(hSig, "data", "LE");
        legend->AddEntry(hBkg, "Mixed Bkg", "LE");
        legend->AddEntry(hBkgNorm, "Normalization Region", "F");
        legend->Draw();
        t2->DrawLatex(0.2, 0.92, "#bf{#Xi(1530)^{0} #rightarrow #Xi + #pi}");
        //t2->DrawLatex(0.75, 0.92, "#bf{pp 13 TeV}");
        t2->DrawLatex(0.55, 0.92, Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", ptbin.at(bin), ptbin.at(bin + 1)));
        if (bin == 0 || bin == 10) {
            t2->DrawLatex(0.45, 0.14, "#color[2]{For Feasibility Check Only}");
        }
        if ((bin == 3) || (bin == 1))
            SavePad((TPad*)gPad, Form("hSigBkg_%d", bin),
                    Form("%s/", savepath.Data()), savetype);
    }
    SaveCanvas(cSigBkg, "hSigBkg", Form("%s/", savepath.Data()), savetype);

    // Fitted result
    TCanvas* cFit = new TCanvas("cFit", "cFit", w, h);
    cFit->Draw();
    cFit->Divide(4, 2, 0.00001, 0.00001);
    cFit->SetTickx();
    cFit->SetLogy(false);
    cFit->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    
    vector<double> fitmean;
    vector<double> fitmean_err;
    vector<double> fitsigma;
    vector<double> fitsigma_err;
    for (int bin = 0; bin < ptbin.size()-1; bin++) {  // ptbin = 10
        if ((!corey) && ((bin == 0) || (bin == 10))) {
            fitmean.push_back(1e-25);
            fitsigma.push_back(1e-25);
            fitmean_err.push_back(1e-25);
            fitsigma_err.push_back(1e-25);
            continue;
        }
        cFit->cd(bin);
        gPad->SetLeftMargin(0.1);
        TGaxis::SetMaxDigits(3);

        auto Sigfit = (TH1D*)inputfile->Get(Form("hSignalBkgSubtraction_%i", bin));
        auto fitresult = (TF1*)inputfile->Get(Form("fDataFitResult_%i", bin));
        auto fitbkg = (TF1*)inputfile->Get(Form("fDataFitBkgResult_%i", bin));
        Sigfit->GetYaxis()->SetTitleOffset(1.2);
        Sigfit->GetYaxis()->SetTitleSize(0.04);
        Sigfit->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        fitresult->SetLineWidth(1);
        fitresult->SetLineColor(628);
        fitresult->SetNpx(500);
        fitbkg->SetLineWidth(1);
        fitbkg->SetRange(1.484,1.6);
        fitbkg->SetLineColor(819);
        fitbkg->SetNpx(500);
        Sigfit->Draw("PZ");
        fitbkg->Draw("same");
        fitresult->Draw("same");

        // Legend
        auto legend_fit = new TLegend(0.42, 0.59, 0.7, 0.72);
        legend_fit->SetFillStyle(0);
        legend_fit->AddEntry(Sigfit, "data", "LE");
        legend_fit->AddEntry(fitresult, Form("Voigt+Pol(2)"), "L");
        legend_fit->AddEntry(fitbkg, Form("Pol(2) Bkg"), "L");
        legend_fit->Draw();
        t2->DrawLatex(0.2, 0.92, "#bf{#Xi(1530)^{0} #rightarrow #Xi + #pi}");
        //t2->DrawLatex(0.75, 0.915, "#bf{pp 13 TeV}");
        t2->DrawLatex(0.55, 0.92,
                     Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", ptbin.at(bin),
                          ptbin.at(bin + 1)));
        if (bin == 0 || bin == ptbin.size() - 2) {
            t2->DrawLatex(0.45, 0.14, "#color[2]{For Feasibility Check Only}");
        }
        int bkgfactor = 4;
        if(finputfile.Contains("BkgFit"))
            bkgfactor = 4;
        if(finputfile.Contains("Corey"))
            bkgfactor = 2;
        t2->DrawLatex(0.42, 0.855,
                      Form("#bf{#chi^{2}/NDF: %.1f/%.1f }",
                           fitresult->GetChisquare(), (double)fitresult->GetNDF()));
        t->DrawLatex(
            0.42, 0.81,
            Form("#bf{Mean: %.4f #pm %.2f(#times10^{-4})}", fitresult->GetParameter(bkgfactor),
                 fitresult->GetParError(bkgfactor)*1e4));
        t->DrawLatex(0.42, 0.77,
                     Form("#bf{Sigma(#times10^{3}): %.2f #pm %.2f}",
                          fitresult->GetParameter(bkgfactor + 1)*1e3,
                          fitresult->GetParError(bkgfactor+1)*1e3));
        t->DrawLatex(0.42, 0.73,
                     Form("#bf{Gamma(#times10^{3}): %.2f #pm %.2f}",
                          fitresult->GetParameter(bkgfactor + 2)*1e3,
                          fitresult->GetParError(bkgfactor+2)*1e3));
        if ((bin == 3) || (bin == 1))
            SavePad((TPad*)gPad, Form("hFit_%d", bin),
                    Form("%s/", savepath.Data()), savetype);
        fitmean.push_back(fitresult->GetParameter(bkgfactor));
        fitmean_err.push_back(fitresult->GetParError(bkgfactor));
        fitsigma.push_back(fitresult->GetParameter(bkgfactor+1));
        fitsigma_err.push_back(fitresult->GetParError(bkgfactor + 1));
    }
    SaveCanvas(cFit, "hFit", Form("%s/", savepath.Data()), savetype);

    // Fitted MC
    TCanvas* cMCFit = new TCanvas("cMCFit", "cMCFit", w, h);
    cMCFit->Draw();
    cMCFit->Divide(4, 2, 0.0001, 0.0001);
    cMCFit->SetTickx();
    cMCFit->SetLogy(false);
    cMCFit->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);

    vector<double> fitmeanMC;
    vector<double> fitsigmaMC;
    vector<double> fitmeanMC_err;
    vector<double> fitsigmaMC_err;
    for (int bin = 0; bin < ptbin.size()-1; bin++) {  // ptbin = 10
        if ((!corey) && ((bin == 0) || (bin == 10))) {
            fitmeanMC.push_back(1e-25);
            fitsigmaMC.push_back(1e-25);
            fitmeanMC_err.push_back(1e-25);
            fitsigmaMC_err.push_back(1e-25);
            continue;
        }
        cFit->cd(bin);
        TGaxis::SetMaxDigits(3);

        auto hReco = (TH1D*)inputfile->Get(Form("hMCRecon_%i", bin));
        auto fitresult = (TF1*)inputfile->Get(Form("fMCFit_%i", bin));
        hReco->GetXaxis()->SetRangeUser(DrawRange[0], DrawRange[1]);
        hReco->Draw("PZ");
        fitresult->Draw("same");

        // Legend
        auto legend_MCfit = new TLegend(0.13, 0.76, 0.5, 0.89);
        legend_MCfit->SetFillStyle(0);
        legend_MCfit->AddEntry(hReco, "MC recon", "LE");
        legend_MCfit->AddEntry(fitresult, Form("Voigt only"), "L");
        legend_MCfit->Draw();
        t2->DrawLatex(0.2, 0.92, "#bf{#Xi(1530)^{0} #rightarrow #Xi + #pi}");
        //t2->DrawLatex(0.75, 0.915, "#bf{pp 13 TeV}");
        t2->DrawLatex(0.55, 0.92,
                      Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", ptbin.at(bin),
                           ptbin.at(bin + 1)));
        t->DrawLatex(0.14, 0.71,
                     Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", ptbin.at(bin),
                          ptbin.at(bin + 1)));
        if (bin == 0 || bin == ptbin.size() - 2) {
            t2->DrawLatex(0.45, 0.14, "#color[2]{For Feasibility Check Only}");
        }
        t->DrawLatex(0.5, 0.81,
                     Form("#bf{Mean: %.4f}", fitresult->GetParameter(2)));
        t->DrawLatex(0.5, 0.77,
                     Form("#bf{Sigma: %.4f}", fitresult->GetParameter(3)));
        t->DrawLatex(0.5, 0.73,
                     Form("#bf{Gamma: %.4f}", fitresult->GetParameter(4)));
        fitmeanMC.push_back(fitresult->GetParameter(2));
        fitsigmaMC.push_back(fitresult->GetParameter(3));
        fitmeanMC_err.push_back(fitresult->GetParError(2));
        fitsigmaMC_err.push_back(fitresult->GetParError(3));
    }
    SaveCanvas(cFit, "hMCFit", Form("%s/", savepath.Data()), savetype);

    //=========Result Plots=============
    TCanvas* clogy = new TCanvas("clogy", "", w, h);
    clogy->Draw();
    clogy->SetTickx();
    clogy->SetLogy(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas* cNology = new TCanvas("cNology", "", w, h);
    cNology->Draw();
    cNology->SetTickx();
    cNology->SetLogy(false);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Fit results
    TH1D* hFitMean_data =
        new TH1D("Fit Mean", "", ptbin.size() - 1, ptbin_array);
    TH1D* hFitMean_MC = new TH1D("Fit Mean MC", "", ptbin.size() - 1, ptbin_array);
    TH1D* hFitSigma_data =
        new TH1D("Fit Sigma", "", ptbin.size() - 1, ptbin_array);
    TH1D* hFitSigma_MC =
        new TH1D("Fit Sigma MC", "", ptbin.size() - 1, ptbin_array);
    for (int i = 0; i < fitmean.size(); i++) {
        hFitMean_data->SetBinContent(i + 1, fitmean.at(i));
        hFitMean_data->SetBinError(i + 1, fitmean_err.at(i));

        hFitMean_MC->SetBinContent(i + 1, fitmeanMC.at(i));
        hFitMean_MC->SetBinError(i + 1, fitmeanMC_err.at(i));

        hFitSigma_data->SetBinContent(i + 1, fitsigma.at(i));
        hFitSigma_data->SetBinError(i + 1, fitsigma_err.at(i));

        hFitSigma_MC->SetBinContent(i + 1, fitsigmaMC.at(i));
        hFitSigma_MC->SetBinError(i + 1, fitsigmaMC_err.at(i));
    }
    hFitMean_MC->SetMarkerColor(628);
    hFitMean_MC->SetLineColor(628);
    hFitSigma_MC->SetMarkerColor(628);
    hFitSigma_MC->SetLineColor(628);

    cNology->cd();
    hFitMean_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFitMean_data->GetYaxis()->SetTitle("Fit mean");
    hFitMean_data->SetMinimum(1.53);
    hFitMean_data->SetMaximum(1.536);
    hFitMean_data->GetXaxis()->SetRangeUser(0.8, 7.4);
    hFitMean_MC->GetXaxis()->SetRangeUser(0.8, 7.4);
    hFitMean_data->Draw("E");
    hFitMean_MC->Draw("E SAME");
    TF1* pdg_massline =
        new TF1("pdg_massline", "1.53178", -1, ptbin.back() + 1);
    pdg_massline->SetLineWidth(1);
    pdg_massline->SetLineStyle(2);
    pdg_massline->Draw("same");

    auto legend_Mean = new TLegend(0.55, 0.70, 0.95, 0.85);
    legend_Mean->SetFillStyle(0);
    legend_Mean->AddEntry(hFitMean_data, "Data", "LE");
    legend_Mean->AddEntry(hFitMean_MC, "MC", "LE");
    legend_Mean->AddEntry(pdg_massline, "PDG Mass(1.53178)", "L");
    legend_Mean->Draw();

    SaveCanvas(cNology, "FitMean", Form("%s/", savepath.Data()), savetype);

    hFitSigma_data->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hFitSigma_data->GetYaxis()->SetTitle("Fit sigma");
    hFitSigma_data->SetMinimum(0);
    hFitSigma_data->SetMaximum(0.01);
    hFitSigma_data->GetXaxis()->SetRangeUser(0.8, 7.4);
    hFitSigma_MC->GetXaxis()->SetRangeUser(0.8, 7.4);
    hFitSigma_data->Draw("E");
    hFitSigma_MC->Draw("E SAME");

    auto legend_sigma = new TLegend(0.55, 0.75, 0.95, 0.85);
    legend_sigma->SetFillStyle(0);
    legend_sigma->AddEntry(hFitSigma_data, "Data", "LE");
    legend_sigma->AddEntry(hFitSigma_MC, "MC", "LE");
    legend_sigma->Draw();
    SaveCanvas(cNology, "FitSigma", Form("%s/", savepath.Data()), savetype);

    // MC Recon effi.
    clogy->cd();
    auto hMCReconEffi = (TH1D*)inputfile->Get("hMCReconEffi");

    vector<double> MCEfficiency_7TeV = {1e-10,   0.006046, 0.02355, 0.04314,
                                  0.07004, 0.1002,   0.1412,  0.1579,
                                  0.1569,  1e-10,    1e-10};
    TH1D* MCeffi_7TeV =
        new TH1D("7TeV MC Efficiency", "", ptbin.size() - 1, ptbin_array);
    for (int i = 0; i < MCEfficiency_7TeV.size(); i++) {
        MCeffi_7TeV->SetBinContent(i + 1, MCEfficiency_7TeV.at(i));
        MCeffi_7TeV->SetBinError(i + 1, 1e-20);
    }

    MCeffi_7TeV->SetMarkerColor(2);
    MCeffi_7TeV->SetLineColor(2);
    MCeffi_7TeV->SetMinimum(4e-4);
    MCeffi_7TeV->SetMaximum(4e-1);
    MCeffi_7TeV->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    MCeffi_7TeV->GetYaxis()->SetTitle("Acceptance x Efficiency x BR");
    MCeffi_7TeV->GetXaxis()->SetRangeUser(0.8, 7.4);
    MCeffi_7TeV->Draw("E");
    hMCReconEffi->Draw("SAME E");

    auto legend_MCEffi = new TLegend(0.55, 0.25, 0.85, 0.40);
    legend_MCEffi->SetFillStyle(0);
    legend_MCEffi->AddEntry(hMCReconEffi, Form("13 TeV %.2f-%.2f", multi_start, multi_end), "LE");
    legend_MCEffi->AddEntry(MCeffi_7TeV, "7 TeV MB", "LE");
    legend_MCEffi->Draw();
    t2->DrawLatex(0.12, 0.92, "#bf{#Xi(1530)^{0} #rightarrow #Xi + #pi}");
    SaveCanvas(clogy, "MCEfficiency", Form("%s/", savepath.Data()), savetype);

    // Singal loss
    cNology->cd();
    auto hMCSigLoss = (TH1D*)inputfile->Get("hMCSigLoss");
    hMCSigLoss->SetMaximum(1.1);
    hMCSigLoss->SetMinimum(0.9);
    hMCSigLoss->Draw("E TEXT");
    SaveCanvas(cNology, "SignalLoss", Form("%s/", savepath.Data()), savetype);

    // Raw yield spectra
    clogy->cd();
    auto hDataRawYield = (TH1D*)inputfile->Get("hDataRawYield");
    hDataRawYield->GetXaxis()->SetRangeUser(0.8, 7.4);
    hDataRawYield->Draw("HIST TEXT");
    SaveCanvas(clogy, "hRawSpectra", Form("%s/", savepath.Data()), savetype);

    // Corrected spectra
    auto hXiSpectrum = (TH1D*)inputfile->Get("hXiSpectrum");
    hXiSpectrum->Draw("E");
    hXiSpectrum->Draw("E");
    SaveCanvas(clogy, "hXiSpectrum", Form("%s/", savepath.Data()), savetype);

    AliPWGFunc * fm = new AliPWGFunc;
    fm->SetVarType(AliPWGFunc::VarType_t(AliPWGFunc::kdNdpt));
    TF1* func = 0;
    
    func = fm->GetLevi (1.5318, 0.4, 750,3);
    func->SetParLimits(1,0.0001,20000);
    /*
    func = fm->GetBGBW(1.5318 ,0.6,0.3, 1, 1e5);// beta, T, n, norm 
    func->SetParLimits(1, 0.1, 0.99);
    func->SetParLimits(2, 0.01, 1);
    func->SetParLimits(3, 0.01, 2);
    */
    /*
    TF1* myLevy = new TF1("myLevy", myLevyPt, 0, 8.5, 3);
    myLevy->SetParName(0, "dN/dy");
    myLevy->SetParName(1, "C");
    myLevy->SetParName(2, "n");
    myLevy->SetParameter(0, .008);
    myLevy->SetParameter(1, .3);
    myLevy->SetParameter(2, 15);
    myLevy->SetParLimits(0, .001, .2);
    myLevy->SetParLimits(1, .1, 1);
    myLevy->SetParLimits(2, 1, 500);
    */
    hXiSpectrum->Fit(func, "IME", "", 0.8, 8.8);
    /*
    TF1* Levy_full = new TF1("myLevy", myLevyPt, 0, 8.5, 3);
    Levy_full->SetParameter(0, func->GetParameter(0));
    Levy_full->SetParameter(1, func->GetParameter(1));
    Levy_full->SetParameter(2, func->GetParameter(2));
    */
    TF1* Levy_full = fm->GetLevi (1.5318, 0.4, 750,3);
    //TF1* Levy_full = fm->GetBGBW(1.5318 ,0.6,0.3, 1, 1e5);// beta, T, n, norm 
    Levy_full->SetParameter(0, func->GetParameter(0));
    Levy_full->SetParameter(1, func->GetParameter(1));
    Levy_full->SetParameter(2, func->GetParameter(2));
    Levy_full->SetRange(0,10);
    hXiSpectrum->Draw("E");
    Levy_full->Draw("same");

    auto legend_Spetra_fit = new TLegend(0.57, 0.70, 0.85, 0.83);
    legend_Spetra_fit->SetFillStyle(0);
    legend_Spetra_fit->AddEntry(
        hXiSpectrum,
        Form("#Xi(1530)^{0} yield 13TeV (%.2f-%.2f)", multi_start, multi_end),
        "PEL");
    legend_Spetra_fit->AddEntry(Levy_full, "Levy fit", "L");
    legend_Spetra_fit->Draw();
    t2->DrawLatex(0.18, 0.92, "#bf{#Xi(1530)^{0} #rightarrow #Xi + #pi}");
    t2->DrawLatex(0.75, 0.92, "#bf{pp 13 TeV}");
    t->DrawLatex(
        0.59, 0.67,
        Form("#bf{dN/dy: %.3f #pm %.3f (x10^{-3})}",
             func->GetParameter(0) * 1e3, func->GetParError(0) * 1e3));
    t->DrawLatex(
        0.59, 0.62,
        Form("#bf{n: %.2f #pm %.2f}", func->GetParameter(1),
             func->GetParError(1)));
    t->DrawLatex(0.59, 0.57,
                 Form("#bf{T: %.2f #pm %.2f}", func->GetParameter(2),
                      func->GetParError(2)));
    SaveCanvas(clogy, "hSpectra_Fit", Form("%s/", savepath.Data()), savetype);

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

    hSpectra_7TeV_syserr->SetLineColor(kRed);
    hSpectra_7TeV_syserr->Draw("same");

    SaveCanvas(clogy, "hSpectra_Fit_7TeV", Form("%s/", savepath.Data()), savetype);

    cout << "===========================" << endl;
    cout << "\\begin{table}[]" << endl;
    cout << "\\centering" << endl;
    cout << "\\begin{tabular}{|l|l|}" << endl;
    cout << "\\hline" << endl;
    cout << "\\ensuremath{p_{\\mathrm{T}}} & Raw yield \\\\ \\hline"
         << endl;
    for (int bin = 1; bin < hDataRawYield->GetNbinsX(); bin++) {
        cout << "$" << ptbin[bin-1] << " - " << ptbin[bin] << "$ & $"
             << hDataRawYield->GetBinContent(bin) << " \\pm "
             << hDataRawYield->GetBinError(bin) << "$ \\\\" << endl;
    }
    cout << "\\hline" << endl;
    cout << "\\end{tabular}" << endl;
    cout << "\\caption{Rawyield in " << multi_start << " - " << multi_end << " bin}" << endl;
    cout << "\\label{tab:rawyield" << multi_start << multi_end << "}" << endl;
    cout << "\\end{table}" << endl;
}
double myLevyPt(Double_t* x, Double_t* par) {
    double lMass = 1.5318;  // Xi mass

    Double_t ldNdy = par[0];          // dN/dy
    Double_t l2pi = 2 * TMath::Pi();  // 2pi
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