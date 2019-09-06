#include <AliPWGHistoTools.h>
#include "PlotXi1530.C"

// Default Canvas
Double_t w = 1920;
Double_t h = 1080;
Int_t colors[] = {633, 810, 800, 830, 840, 840, 870, 864, 890, 617, 619};
vector<int> vcolors;
vector<int> Markerset={24, 21, 33, 34, 29, 2, 25, 26, 27};
vector<int> Markerset_size={1, 1, 2, 2, 2, 1, 1, 1, 1, 1};
const double inf = 1e20;
const double zero = 1e-20;
bool save = false;
TString path = "./data/";

bool fUseCutSysfrom0100 = true;
vector<double> TopolCutSysSum;
vector<double> PIDCutSysSum;

vector<double> ptbin = {0.0, 0.8, 1.2, 1.6, 2.0, 2.4,
                        3.2, 4.0, 4.8, 5.6, 8.8, 15};
vector<double> multibin = {0, 100};
vector<double> multibin0100 = {0, 100};

TLatex* t_big = new TLatex();
// for memo, small
TLatex* t = new TLatex();
vector<double> baseStatError;

//for sysplot
TCanvas* cPlotOut = new TCanvas("cPlotOut", "cPlotOut", w, h);
auto legeondPlotOut = new TLegend(0.13, 0.60, 0.5, 0.88);
int linestyle = 0;
bool GlobalFirst = true;
vector<double> row;
vector<vector<double>> ptbinStatErrorAvg(ptbin.size(), row);

TH1D* CheckBarlowStatistic(TString fvariationname = "NormRange",
                           vector<double> fmultibin = {0, 10, 30, 50, 70, 100},
                           int options_barlow = 1,
                           int fitvar = 1);
TH1D* Xi1530PlotSystematics_multi(double multi_start = 0,
                                  double multi_end = 100,
                                  bool correl_error_skip = false);
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
vector<TH1D*> Xi1530SysFitVarCheck(TString nameFitVar = "BinCountR");
TH1D* Xi1530SysFitVar(TString nameFitVar = "FitVar",
                      double multi_start = 0,
                      double multi_end = 100);
vector<TH1D*> Xi1530SysCheck(TString fdeafultinputfile,
                             TString fvariationinputfile,
                             char const* options = "Var");
void Xi1530PlotSystematics() {
    vector<double> multibinsloop = {0, 10, 30, 50, 70, 100};
    //vector<double> multibinsloop = {70, 100};
    vector<vector<double>> multibincheck = {
        //{0, 100}};
        //{0, 100}, {0, 0.1}, {0, 10},  {10, 30}, {30, 50}, {50, 70}, {70, 100}};
        {0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    //{0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 100}};
    TFile* fout_systematic = new TFile("./AnalysisResults_Xi1530_systematic.root", "RECREATE");
    t_big->SetNDC();
    t_big->SetTextSize(0.07);
    t->SetNDC();
    t->SetTextSize(0.05);
    TCanvas* cCanvas = new TCanvas("cCanvas", "cCanvas", w, h);
    //gStyle->SetOptTitle(0);
    cCanvas->SetTickx();
    cCanvas->Draw();
    cCanvas->cd();
    cCanvas->SetLogy(true);
    cout << "START OF DEFAULT CONFIGURE =========================" << endl;
    auto hDef = Xi1530PlotSystematics_multi(0, 100, false);
    fout_systematic->cd();
    //hDef->Write("hSpectra_0.00_100.00_sys");
    vector<TH1D*> buffer;
    cout << "END OF DEFAULT CONFIGURE =========================" << endl;

    vector<double> ptbin_final = {0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6, 8.8};
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        auto temp_multi = Xi1530PlotSystematics_multi(
            multibincheck[imultibin][0], multibincheck[imultibin][1], false);
        vector<double> temp_y;
        vector<double> temp_e;
        for (int bin = 0; bin <= ptbin.size(); bin++) {
            // barlow calculation
            if( (bin == 0) || (bin == 10) || (bin == 11) ) continue;
            temp_y.push_back(temp_multi->GetBinContent(bin+1));
            temp_e.push_back(temp_multi->GetBinError(bin + 1));
        }
        auto htemp = (TH1D*)MakeHistfromArray(Form("h_%.2f-%.2f_sys", multibincheck[imultibin][0],
                               multibincheck[imultibin][1]),
                                       temp_y, ptbin_final, temp_e);
        buffer.push_back(htemp);
        fout_systematic->cd();
        htemp->Write(Form("hSpectra_%.2f_%.2f_sys", multibincheck[imultibin][0],
                          multibincheck[imultibin][1]));
    }
    vector<TH1D*> buffer_multi;
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        auto temp_spectra_multi = GetSpectrafromName(Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s%i.root",
                 path.Data(), 1, multibincheck[imultibin][0], multibincheck[imultibin][1],
                 "Default", 1));
        vector<double> temp_y;
        vector<double> temp_e;
        for (int bin = 0; bin <= ptbin.size(); bin++) {
            // barlow calculation
            if ((bin == 0) || (bin == 10) || (bin == 11))
                continue;
            temp_y.push_back(temp_spectra_multi->GetBinContent(bin + 1));
            temp_e.push_back(temp_spectra_multi->GetBinError(bin + 1));
        }
        auto htemp = (TH1D*)MakeHistfromArray(Form("h_%.2f-%.2f_stat", multibincheck[imultibin][0],
                               multibincheck[imultibin][1]),
                                       temp_y, ptbin_final, temp_e);
        buffer_multi.push_back(htemp);
        fout_systematic->cd();
        htemp->Write(Form("hSpectra_%.2f_%.2f_stat",
                          multibincheck[imultibin][0],
                          multibincheck[imultibin][1]));
    }
    /*
    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    Int_t FIh = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 multibincheck.size());
    cCanvas->cd();
    for (int imultibin = 0; imultibin < buffer.size() - 1; imultibin++) {
        buffer[imultibin]->SetLineColor(FIh + imultibin);
        //buffer[imultibin]->Scale(pow(2, multibinsloop.size() - imultibin));
        buffer[imultibin]->SetMaximum(1);
        buffer[imultibin]->SetMinimum(1e-7);
        buffer[imultibin]->Draw("E2 same");
    }
    TCanvas* ratio_mult = plotHistsAndRatio(
        buffer, buffer[0], "", "p_{T}(GeV/c)",
        "1/N_{event}d^{2}N/(dydp_{T}) (GeV/c)^{-1}", "POW FINAL");
    ratio_mult->cd(1);
    t->DrawLatex(0.68, 0.92, "#bf{pp #sqrt{#it{s}} = 13 TeV}");
    t->DrawLatex(0.80, 0.86, "#bf{|#it{y}| < 0.5}");
    t_big->DrawLatex(0.70, 0.78, "#Xi(1530)^{0}");
    t->DrawLatex(0.17, 0.27, "#bf{V0M Multiplicity}");

    auto legend_hSpectra_ratio = new TLegend(.16, .03, .58, .27);
    legend_hSpectra_ratio->SetNColumns(2);
    legend_hSpectra_ratio->SetBorderSize(0);
    legend_hSpectra_ratio->SetFillStyle(0);

    for (int imultibin = 0; imultibin < buffer.size(); imultibin++) {
        buffer.at(imultibin)->SetFillColor(vcolors[imultibin]);
        buffer.at(imultibin)->SetLineColor(vcolors[imultibin]);
        buffer.at(imultibin)->SetMarkerColor(vcolors[imultibin]);
        buffer.at(imultibin)->SetMarkerStyle(Markerset[imultibin]);
        buffer.at(imultibin)->SetMarkerSize(Markerset_size[imultibin]);
        legend_hSpectra_ratio->AddEntry(
            buffer.at(imultibin),
            Form("%.2f - %.2f (x2^{%d})", multibincheck[imultibin+1][0],
                 multibincheck[imultibin+1][1], int(buffer.size() - imultibin)),
            "PL");
    }
    legend_hSpectra_ratio->Draw();

    cCanvas->SaveAs("./test.pdf");
    ratio_mult->SaveAs("./test2.pdf");
    
    */
    // Make a gaussian with 

    

    TCanvas* cCanvas_barlow =
        new TCanvas("cCanvas_barlow", "cCanvas_barlow", 960, 720);
    cCanvas_barlow->SetTickx();
    cCanvas_barlow->SetTicky();
    cCanvas_barlow->SetTopMargin(0.05);
    cCanvas_barlow->SetLeftMargin(0.1);
    cCanvas_barlow->SetBottomMargin(0.1);
    cCanvas_barlow->SetRightMargin(0.01);
    cCanvas_barlow->SetFillStyle(0);
    
    gStyle->SetOptStat(11111);
    TCanvas* cCanvas_barlow_KSCheck =
        new TCanvas("cCanvas_barlow_KSCheck", "cCanvas_barlow_KSCheck", w, h);
    cCanvas_barlow_KSCheck->SetTickx();

    cCanvas_barlow->Draw();
    cCanvas_barlow->cd();
    double fullintegral;
    double I1integral;
    double I2integral;
    // -----------------------------------------------------------------------
    vector<TString> TopologicalCutSystematic_bins = {"TPCNsigmaXi1530PionLoose",
                                          "TPCNsigmaXi1530PionTight",
                                          "TPCNsigmaXiLoose",
                                          "TPCNsigmaXiTight",
                                          "Xi1530PionZVertexLoose",
                                          "Xi1530PionZVertexTight",
                                          "DCADistLambdaDaughtersLoose",
                                          "DCADistLambdaDaughtersTight",
                                          "DCADistXiDaughtersLoose",
                                          "DCADistXiDaughtersTight",
                                          "DCADistLambdaPVLoose",
                                          "DCADistLambdaPVTight",
                                          "V0CosineOfPointingAngleLoose",
                                          "V0CosineOfPointingAngleTight",
                                          "CascadeCosineOfPointingAngleLoose",
                                          "CascadeCosineOfPointingAngleTight",
                                          "XiMassWindowLoose",
                                          "XiMassWindowTight"};
    for (int topological_cutvar = 0; topological_cutvar < TopologicalCutSystematic_bins.size(); topological_cutvar++) {
        auto hbarlowcheck =
            (TH1D*)CheckBarlowStatistic(TopologicalCutSystematic_bins[topological_cutvar].Data(), multibinsloop, topological_cutvar+2);
        // K-S Check
        int nEntries = hbarlowcheck->GetEntries();
        auto hKSresults =
            new TH1F("hKSresult", "Result of KS test prob", 1000, 0, 1);
        for (Int_t irnd = 0; irnd < 1000; irnd++) {
            auto hRndGaus = new TH1F(Form("hRndGaus%d", irnd),
                                     "Random Normaldistribution", 200, -10, 10);
            hRndGaus->FillRandom("gaus", nEntries);
            hKSresults->Fill(hbarlowcheck->KolmogorovTest(hRndGaus));
        }
        cCanvas_barlow_KSCheck->cd();
        hKSresults->Draw();
        t->DrawLatex(0.65, 0.52,
                     Form("KS test: %.3f +- %.3f", hKSresults->GetMean(),
                          hKSresults->GetStdDev()));
        /*
        cCanvas_barlow_KSCheck->SaveAs(
            Form("figs/ConsistencyCheck/check_%s_KS.pdf",
                 TopologicalCutSystematic_bins[topological_cutvar].Data()));
        */
        //
        cCanvas_barlow->cd();
        hbarlowcheck->Draw();
        hbarlowcheck->GetYaxis()->SetTitle(TopologicalCutSystematic_bins[topological_cutvar].Data());
        hbarlowcheck->GetXaxis()->SetTitle(
            "Barlow Criteria(#Delta/#sigma_{cc})");
        //Consistency check
        fullintegral =
            hbarlowcheck->Integral(hbarlowcheck->GetXaxis()->FindBin(-99),
                                   hbarlowcheck->GetXaxis()->FindBin(99));
        I1integral =
            hbarlowcheck->Integral(hbarlowcheck->GetXaxis()->FindBin(-1),
                                   hbarlowcheck->GetXaxis()->FindBin(1));
        I2integral =
            hbarlowcheck->Integral(hbarlowcheck->GetXaxis()->FindBin(-2),
                                   hbarlowcheck->GetXaxis()->FindBin(2));
        auto I1value = I1integral / fullintegral;
        auto I2value = I2integral / fullintegral;
        TString I1text;
        TString I2text;
        int check = 0;
        if(I1value > 0.55){
            I1text = Form("#bf{I1= %.3f / Pass!}", I1value);
            check++;
        }
        else
            I1text = Form("#bf{I1= %.3f / Fail!}", I1value);

        if (I2value > 0.75){
            I2text = Form("#bf{I2= %.3f / Pass!}", I2value);
            check++;
        }
        else
            I2text = Form("#bf{I2= %.3f / Fail!}", I2value);

        TString Meantext;
        TString Stdtext;
        if (abs(hbarlowcheck->GetMean()) < 0.8){
            Meantext = Form("#bf{Mean: %.3f / Pass!}", hbarlowcheck->GetMean());
            check++;
        }
        else
            Meantext = Form("#bf{Mean: %.3f / Fail!}", hbarlowcheck->GetMean());

        if (hbarlowcheck->GetStdDev() < 1.3){
            Stdtext =
                Form("#bf{StdDev: %.3f / Pass!}", hbarlowcheck->GetStdDev());
            check++;;
        }
        else
            Stdtext = Form("#bf{StdDev: %.3f / Fail!}", hbarlowcheck->GetStdDev());

        t->DrawLatex(0.15, 0.87, Meantext);
        t->DrawLatex(0.15, 0.82, Stdtext);
        t->DrawLatex(0.15, 0.77, I1text);
        t->DrawLatex(0.15, 0.72, I2text);
        if (check > 2)
            t->DrawLatex(0.15, 0.67, "Consistency Check pass!");

        // Lines
        TLine* rightI1 = new TLine(1, 0, 1, hbarlowcheck->GetMaximum());
        TLine* leftI1 = new TLine(-1, 0, -1, hbarlowcheck->GetMaximum());
        TLine* rightI2 = new TLine(2, 0, 2, hbarlowcheck->GetMaximum());
        TLine* leftI2 = new TLine(-2, 0, -2, hbarlowcheck->GetMaximum());

        rightI1->SetLineColor(2);
        leftI1->SetLineColor(2);
        rightI2->SetLineColor(3);
        leftI2->SetLineColor(3);

        rightI1->Draw();
        leftI1->Draw();
        rightI2->Draw();
        leftI2->Draw();
        //
        cCanvas_barlow->SaveAs(
            Form("figs/ConsistencyCheck/check_%s.pdf", TopologicalCutSystematic_bins[topological_cutvar].Data()));
    }
    
    
    //--------------------------------------------------
    vector<TString> FitSystematic_bins = {"BkgFit", "LikeSignBkg", "FitVarLm",
                                          "BinCountLm", "NormVarLm"};
    for (int topological_cutvar = 0; topological_cutvar < FitSystematic_bins.size(); topological_cutvar++) {
        TH1D* hbarlowcheck;
        if (FitSystematic_bins[topological_cutvar].Contains("Var") ||
            FitSystematic_bins[topological_cutvar].Contains("Count"))
            hbarlowcheck = (TH1D*)CheckBarlowStatistic(
                FitSystematic_bins[topological_cutvar].Data(), multibinsloop,
                1, 3);
        else hbarlowcheck = (TH1D*)CheckBarlowStatistic(
            FitSystematic_bins[topological_cutvar].Data(), multibinsloop, 1, 1);
        // K-S Check
        int nEntries = hbarlowcheck->GetEntries();
        auto hKSresults =
            new TH1F("hKSresult", "Result of KS test prob", 1000, 0, 1);
        for (Int_t irnd = 0; irnd < 1000; irnd++) {
            auto hRndGaus = new TH1F(Form("hRndGaus%d", irnd),
                                     "Random Normaldistribution", 200, -10, 10);
            hRndGaus->FillRandom("gaus", nEntries);
            hKSresults->Fill(hbarlowcheck->KolmogorovTest(hRndGaus));
        }
        cCanvas_barlow_KSCheck->cd();
        hKSresults->Draw();
        t->DrawLatex(0.65, 0.52,
                     Form("KS test: %.3f +- %.3f", hKSresults->GetMean(),
                          hKSresults->GetStdDev()));
        /*
        cCanvas_barlow_KSCheck->SaveAs(
            Form("figs/ConsistencyCheck/check_%s_KS.pdf",
                 FitSystematic_bins[topological_cutvar].Data()));
        */
        //

        cCanvas_barlow->cd();
        hbarlowcheck->Draw();
        hbarlowcheck->GetYaxis()->SetTitle(
            FitSystematic_bins[topological_cutvar].Data());
        hbarlowcheck->GetXaxis()->SetTitle(
            "Barlow Criteria(#Delta/#sigma_{cc})");

        // Consistency check
        fullintegral =
            hbarlowcheck->Integral(hbarlowcheck->GetXaxis()->FindBin(-99),
                                   hbarlowcheck->GetXaxis()->FindBin(99));
        I1integral =
            hbarlowcheck->Integral(hbarlowcheck->GetXaxis()->FindBin(-1),
                                   hbarlowcheck->GetXaxis()->FindBin(1));
        I2integral =
            hbarlowcheck->Integral(hbarlowcheck->GetXaxis()->FindBin(-2),
                                   hbarlowcheck->GetXaxis()->FindBin(2));
        auto I1value = I1integral / fullintegral;
        auto I2value = I2integral / fullintegral;
        TString I1text;
        TString I2text;
        int check = 0;
        if (I1value > 0.55) {
            I1text = Form("#bf{I1= %.3f / Pass!}", I1value);
            check++;
        } else
            I1text = Form("#bf{I1= %.3f / Fail!}", I1value);

        if (I2value > 0.75) {
            I2text = Form("#bf{I2= %.3f / Pass!}", I2value);
            check++;
        } else
            I2text = Form("#bf{I2= %.3f / Fail!}", I2value);

        TString Meantext;
        TString Stdtext;
        if (abs(hbarlowcheck->GetMean()) < 0.8) {
            Meantext = Form("#bf{Mean: %.3f / Pass!}", hbarlowcheck->GetMean());
            check++;
        } else
            Meantext = Form("#bf{Mean: %.3f / Fail!}", hbarlowcheck->GetMean());

        if (hbarlowcheck->GetStdDev() < 1.3) {
            Stdtext =
                Form("#bf{StdDev: %.3f / Pass!}", hbarlowcheck->GetStdDev());
            check++;
            ;
        } else
            Stdtext =
                Form("#bf{StdDev: %.3f / Fail!}", hbarlowcheck->GetStdDev());

        t->DrawLatex(0.15, 0.87, Meantext);
        t->DrawLatex(0.15, 0.82, Stdtext);
        t->DrawLatex(0.15, 0.77, I1text);
        t->DrawLatex(0.15, 0.72, I2text);
        if (check > 2)
            t->DrawLatex(0.15, 0.67, "Consistency Check pass!");

        // Lines
        TLine* rightI1 = new TLine(1, 0, 1, hbarlowcheck->GetMaximum());
        TLine* leftI1 = new TLine(-1, 0, -1, hbarlowcheck->GetMaximum());
        TLine* rightI2 = new TLine(2, 0, 2, hbarlowcheck->GetMaximum());
        TLine* leftI2 = new TLine(-2, 0, -2, hbarlowcheck->GetMaximum());

        rightI1->SetLineColor(2);
        leftI1->SetLineColor(2);
        rightI2->SetLineColor(3);
        leftI2->SetLineColor(3);

        rightI1->Draw();
        leftI1->Draw();
        rightI2->Draw();
        leftI2->Draw();
        //
        cCanvas_barlow->SaveAs(
            Form("figs/ConsistencyCheck/check_%s.pdf", FitSystematic_bins[topological_cutvar].Data()));
    }
    

    fout_systematic->Close();
}
TH1D* CheckBarlowStatistic(TString fvariationname,
                           vector<double> fmultibin,
                           int fitsys,
                           int fitvar) {
    // will make a (error between default and varation) / (differency of
    // stat_Error) like a barlow check, but want to see how it will be
    // distributed in all pT/multiplicity bins.
    
    TCanvas* cBarlowCanvas =
        new TCanvas("cBarlowCanvas", "cBarlowCanvas", w, h);
    cBarlowCanvas->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // cTestCanvas->SetLogy(true);
    cBarlowCanvas->Draw();

    
    TString inputOptions_barlow = "Default";
    int cutbin_barlow = 1;

    if(fitsys > 1)
        cutbin_barlow = fitsys;
    // init StatError Array
    double difference_between_two_spectra = 0.;
    double baseStatError_barlow = 0.;
    double VariationStatError_barlow = 0.;
    double error_barlow = 0.;
    double deltasigma = 0.;
    TString fvarfile;
    TString fdeafultinputfile_barlow;
    TH1D* houtputhist =
        new TH1D(Form("%s_variations", fvariationname.Data()),
                 Form("%s_variations", fvariationname.Data()), 200, -10, 10);
    TFile* fDefault_barlow;
    TFile* tempfile;

    for (int fmulti_barlow = 0; fmulti_barlow < fmultibin.size()-1;
         fmulti_barlow++) {
        fdeafultinputfile_barlow =
            Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
                 path.Data(), 1, fmultibin[fmulti_barlow],
                 fmultibin[fmulti_barlow + 1], inputOptions_barlow.Data());
        fDefault_barlow = new TFile(fdeafultinputfile_barlow.Data());
        auto base_barlow = (TH1D*)fDefault_barlow->Get("hXiSpectrum");

        fvarfile =
            Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s%i.root",
                 path.Data(), cutbin_barlow, fmultibin[fmulti_barlow],
                 fmultibin[fmulti_barlow + 1], fvariationname.Data(), fitvar);
        tempfile = new TFile(fvarfile.Data());
        auto htemp_barlow = (TH1D*)tempfile->Get("hXiSpectrum");
        for (int bin = 0; bin <= ptbin.size(); bin++) {
            // barlow calculation
            if( (bin == 0) || (bin == 10) || (bin == 11) || (bin == 12) ) continue;
            

            //difference_between_two_spectra = temp[2]->GetBinContent(bin + 1);
            difference_between_two_spectra =
                htemp_barlow->GetBinContent(bin + 1) -
                base_barlow->GetBinContent(bin + 1);
            baseStatError_barlow = base_barlow->GetBinError(bin + 1);
            VariationStatError_barlow = htemp_barlow->GetBinError(bin + 1);

            deltasigma = sqrt(abs(pow(VariationStatError_barlow, 2) -
                                  pow(baseStatError_barlow, 2)));
            houtputhist->Fill(difference_between_two_spectra / deltasigma);

            cout << fvariationname.Data() << " - "
                 << "Multi: " << fmultibin[fmulti_barlow] << "-"
                 << fmultibin[fmulti_barlow + 1] << ", bin " << bin
                 << ", default vaule: " << base_barlow->GetBinContent(bin + 1)
                 << " +- " << baseStatError_barlow << " (stat err)"
                 << ", var value: " << htemp_barlow->GetBinContent(bin + 1)
                 << " +- " << VariationStatError_barlow << " (stat err)"
                 << ", value: " << difference_between_two_spectra / deltasigma
                 << endl;
        }
    }
    return houtputhist;
}

TH1D* Xi1530PlotSystematics_multi(double multi_start,
                                  double multi_end,
                                  bool correl_error_skip) {
    // vector<double> multibins = {0, 10, 30, 50, 70, 100};
    cout << "Input multi start: " << multi_start << ", multi end: " << multi_end
          << endl;
    multibin = {multi_start, multi_end};
    vector<double> multibin_local = {multi_start, multi_end};
    //vector<TH1D*> returnHist;
    vector<TH1D*> totalsystematic; // Systematics, will be used!
    
    // smoothing?
    bool smoothing = true;
    double smoothingCriteron = 1.3;

    //
    TFile* fout = new TFile(
        Form("Systematics_%2.f-%2.f.root", multi_start, multi_end), "RECREATE");
    TCanvas* cPlotCanvas = new TCanvas("cPlotCanvas", "cPlotCanvas", w, h);
    cPlotCanvas->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // cTestCanvas->SetLogy(true);
    cPlotCanvas->Draw();

    // Import base -> (TH1D*)base
    TString inputOptions = "Default";
    int cutbin = 1;
    TString finputfile = Form(
        "%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root", path.Data(),
        cutbin, multibin_local[0], multibin_local[1], inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* base = (TH1D*)inputfile->Get("hXiSpectrum");
    TH1D* base_for_error = (TH1D*)inputfile->Get("hXiSpectrum");

    // inputfile->Close();
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    // init Error Array
    vector<double> ErrorArray;
    for (int bin = 0; bin <= ptbin.size(); bin++)
        ErrorArray.push_back(0);

    baseStatError = {};
    cout << "File: " << finputfile << endl;
    for (int bin = 0; bin < ptbin.size(); bin++){
        if(bin < 1) baseStatError.push_back(0);
        else if(bin > 9) baseStatError.push_back(0);
        else{
            double tempe = base_for_error->GetBinError(bin + 1);
            double tempy = base_for_error->GetBinContent(bin + 1);
            double tempr = tempe/tempy;
            cout << "Tempe: " << tempe << ", tempy: " << tempy << ", tempe/tempy: " << tempr << endl;
            baseStatError.push_back(tempr);
        }
    }
    auto hStatError = (TH1D*)MakeHistfromArray("hStatError", baseStatError, ptbin);
    hStatError->SetFillColorAlpha(kBlack, 0.1);
    hStatError->SetLineColorAlpha(kBlack, 0.1);
    hStatError->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    cPlotCanvas->cd();
    hStatError->Draw("BAR");
    hStatError->SetMaximum(0.3);
    //cPlotCanvas->SaveAs(Form("./Check_%2.f-%2.f.pdf", multibin[0], multibin[1]));
    // Fit variations
    //vector<TString> FitSystematic_bins = {"FitRange", "NormRange", "NormRight","NormBoth", "BinCount"};
    //vector<TString> FitSystematic_bins = {"NormRange", "NormRight", "BinCount"};
    vector<TString> FitSystematic_bins = {"FitVar", "NormVar", "BinCount"};
    //vector<TString> FitSystematic_bins = {"FitVar"};  // "NormVar", "BinCount" has been disabled.
    vector<TString> FitSystematic_bins2 = {
        "LikeSignBkg",
        "BkgFit"};  // LikeSignBkg, BkgFit "LikeSignBkg", "BkgFit"
    vector<TH1D*> hFitsys_fraction;

    // Maximum Error
    vector<double> sigExMaxError;
    for (int bin = 0; bin < ptbin.size(); bin++) {
        sigExMaxError.push_back(0);
    }

    auto hFitSysSum = (TH1D*)MakeHistfromArray("FitSys", ErrorArray, ptbin);
    TString fvarfile;
    const int colorbins = FitSystematic_bins.size() + FitSystematic_bins2.size();
    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    Int_t FIf = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 colorbins);

    auto FitSyslegend = new TLegend(0.6, 0.6, 0.9, 0.88);
    FitSyslegend->SetNColumns(2);
    FitSyslegend->SetBorderSize(0);
    FitSyslegend->SetFillStyle(0);
    FitSyslegend->SetFillStyle(0);
    for (int fitvar2 = 0; fitvar2 < FitSystematic_bins2.size(); fitvar2++) {
        fvarfile =
            Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
                 path.Data(), cutbin, multibin[0], multibin[1],
                 FitSystematic_bins2[fitvar2].Data());
        vector<TH1D*> temp = Xi1530SysCheck(
            finputfile, fvarfile, FitSystematic_bins2[fitvar2].Data());
        cPlotCanvas->cd();
        hFitsys_fraction.push_back(temp[0]);
        hFitsys_fraction[fitvar2]->SetLineColor(FIf + fitvar2);
        hFitsys_fraction[fitvar2]->SetLineColor(FIf + fitvar2);
        hFitsys_fraction[fitvar2]->SetMarkerColor(FIf + fitvar2);
        hFitsys_fraction[fitvar2]->SetMaximum(0.3);
        hFitsys_fraction[fitvar2]->Draw("same");
        FitSyslegend->AddEntry(hFitsys_fraction[fitvar2],
                               FitSystematic_bins2[fitvar2], "L");
        for (int bin = 0; bin < ptbin.size(); bin++) {
            if (hFitsys_fraction[fitvar2]->GetBinContent(bin + 1) >
                sigExMaxError[bin])
                sigExMaxError[bin] =
                    hFitsys_fraction[fitvar2]->GetBinContent(bin + 1);
            cout << "Check bin " << ptbin[bin] << " - " << ptbin[bin + 1]
                 << " fitvar " << FitSystematic_bins2[fitvar2].Data() << " : "
                 << hFitsys_fraction[fitvar2]->GetBinContent(bin + 1) << endl;
            hFitSysSum->SetBinContent(
                bin + 1,
                hFitSysSum->GetBinContent(bin + 1) +
                    pow(hFitsys_fraction[fitvar2]->GetBinContent(bin + 1), 2));
        }
    }
    int fitSysVectorSize = hFitsys_fraction.size();
    for (int fitvar = 0; fitvar < FitSystematic_bins.size(); fitvar++) {
        auto fittemp = (TH1D*)Xi1530SysFitVar(FitSystematic_bins[fitvar].Data(),
                                              multibin[0], multibin[1]);
        hFitsys_fraction.push_back(fittemp);

        hFitsys_fraction[fitvar + fitSysVectorSize]->SetLineColor(
            FIf + fitvar + fitSysVectorSize);
        hFitsys_fraction[fitvar + fitSysVectorSize]->SetMarkerColor(
            FIf + fitvar + fitSysVectorSize);
        hFitsys_fraction[fitvar + fitSysVectorSize]->SetMaximum(0.3);

        cPlotCanvas->cd();
        hFitsys_fraction[fitvar + fitSysVectorSize]->Draw("same");
        FitSyslegend->AddEntry(hFitsys_fraction[fitvar + fitSysVectorSize],
                               FitSystematic_bins[fitvar].Data(), "L");
        for (int bin = 0; bin < ptbin.size(); bin++) {
            if (hFitsys_fraction[fitvar + fitSysVectorSize]->GetBinContent(
                    bin + 1) > sigExMaxError[bin])
                sigExMaxError[bin] =
                    hFitsys_fraction[fitvar + fitSysVectorSize]->GetBinContent(
                        bin + 1);
            cout << "Check bin " << bin + 1 << " fitvar "
                 << FitSystematic_bins[fitvar].Data() << " : "
                 << hFitsys_fraction[fitvar + fitSysVectorSize]->GetBinContent(
                        bin + 1)
                 << endl;
            hFitSysSum->SetBinContent(
                bin + 1, hFitSysSum->GetBinContent(bin + 1) +
                             pow(hFitsys_fraction[fitvar + fitSysVectorSize]
                                     ->GetBinContent(bin + 1),
                                 2));
        }
    }
    /*
    // Quadrature sum
    for (int bin = 0; bin < ptbin.size(); bin++) {
        hFitSysSum->SetBinContent(bin+1, sqrt(hFitSysSum->GetBinContent(bin+1)));
        cout << "Check bin(final) " << bin+1 << " : "
             << hFitSysSum->GetBinContent(bin+1) << endl;
    }
    */
    for (int bin = 0; bin < ptbin.size(); bin++) {
        hFitSysSum->SetBinContent(bin + 1,
                                  sigExMaxError[bin]);
        cout << "Check bin(final) " << bin + 1 << " : "
             << hFitSysSum->GetBinContent(bin + 1) << endl;
    }
    if (smoothing) {
        for (int bin = 1; bin < ptbin.size()-1; bin++) {
            if( (ptbin[bin] < 1.2) || (ptbin[bin] > 5.5) ) continue;
            double checkValue = hFitSysSum->GetBinContent(bin + 1);
            double checkValue_before = hFitSysSum->GetBinContent(bin);
            double checkValue_after = hFitSysSum->GetBinContent(bin + 2);
            double smoothedValue = hFitSysSum->GetBinContent(bin + 1);
            double correctionfactor = 3;

            // divided to "2" since one of neighbor is zero.
            double averageOfNeighbor =
                (checkValue + checkValue_before + checkValue_after) /
                correctionfactor;
            if (checkValue > smoothingCriteron * averageOfNeighbor){
                // only if the value is significantly larger than neighbor.
                smoothedValue = averageOfNeighbor;
            }
            if (smoothingCriteron * checkValue < averageOfNeighbor) {
                // only if the value is significantly lower than neighbor.
                smoothedValue = averageOfNeighbor;
                cout << "corrected!" << endl;
            }
            hFitSysSum->SetBinContent(bin + 1, smoothedValue);
        }
    }

    totalsystematic.push_back(hFitSysSum);

    fout->cd();
    cPlotCanvas->cd();
    hFitSysSum->SetLineColor(kBlack);
    hFitSysSum->SetLineStyle(7);
    hFitSysSum->Draw("same");
    hFitSysSum->Write("FitSysSum");
    FitSyslegend->AddEntry(hFitSysSum, "FitSys Sum", "L");
    FitSyslegend->AddEntry(hStatError, "Stat Error", "F");
    FitSyslegend->Draw("same");
    cPlotCanvas->Write(Form("FitVar(%2.f-%2.f)", multi_start, multi_end));
    cPlotCanvas->SaveAs(
        Form("./FitSysFraction(%2.f-%2.f).pdf", multi_start, multi_end));

    // Topological Cut variations
    TCanvas* cPlotCanvas_TopologyCut =
        new TCanvas("cPlotCanvas_TopologyCut", "cPlotCanvas_TopologyCut", w, h);
    cPlotCanvas_TopologyCut->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // cTestCanvas->SetLogy(true);
    cPlotCanvas_TopologyCut->Draw();
    cPlotCanvas_TopologyCut->cd();
    hStatError->Draw("BAR");
    hStatError->SetMaximum(0.3);
    vector<TString> TopologicalCutSystematic_bins = {
        "TPCNsigmaXi1530PionLoose",
        "TPCNsigmaXi1530PionTight",
        "TPCNsigmaXiLoose",
        "TPCNsigmaXiTight",
        "Xi1530PionZVertexLoose",
        "Xi1530PionZVertexTight",
        "DCADistLambdaDaughtersLoose",
        "DCADistLambdaDaughtersTight",
        "DCADistXiDaughtersLoose",
        "DCADistXiDaughtersTight",
        "DCADistLambdaPVLoose",
        "DCADistLambdaPVTight",
        "V0CosineOfPointingAngleLoose",
        "V0CosineOfPointingAngleTight",
        "CascadeCosineOfPointingAngleLoose",
        "CascadeCosineOfPointingAngleTight",
        "XiMassWindowLoose",
        "XiMassWindowTight"};
    vector<TString> TopologicalCutSystematic_notuse = {
        "Xi1530PionZVertexLoose",
        "DCADistLambdaDaughtersLoose",
        "DCADistXiDaughtersLoose",
        "DCADistLambdaPVLoose",
        "V0CosineOfPointingAngleTight",
        "CascadeCosineOfPointingAngleTight",
        "XiMassWindowLoose"};
    vector<TH1D*> hTopologicalCutSys_fraction;
    auto TopologicalCutSyslegend = new TLegend(0.18, 0.6, 0.8, 0.88);
    TopologicalCutSyslegend->SetNColumns(2);
    TopologicalCutSyslegend->SetBorderSize(0);
    TopologicalCutSyslegend->SetFillStyle(0);
    TopologicalCutSyslegend->SetFillStyle(0);

    TH1D* hMultiplicityIndependentTopologicalCutSysSum;
    if ((fUseCutSysfrom0100) && (multi_end - multi_start > 99)){
        cout << "Check Sys Error for Cut sys " << endl;
        hMultiplicityIndependentTopologicalCutSysSum = (TH1D*)MakeHistfromArray("CutSys", ErrorArray, ptbin);
        Int_t FI_topocut = TColor::CreateGradientColorTable(
            NRGBs, stops, red, green, blue, TopologicalCutSystematic_bins.size()-4);
        int arraycount = 0;
        for (int topological_cutvar = 0; topological_cutvar < TopologicalCutSystematic_bins.size(); topological_cutvar++) {
            if (topological_cutvar< 4) 
                continue;
            if (std::find(TopologicalCutSystematic_notuse.begin(),
                          TopologicalCutSystematic_notuse.end(),
                          TopologicalCutSystematic_bins[topological_cutvar]) !=
                TopologicalCutSystematic_notuse.end())
                continue;
            // if (topological_cutvar == 14) continue;
            fvarfile =
                Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
                     path.Data(), topological_cutvar + 2, multibin[0], multibin[1],
                     TopologicalCutSystematic_bins[topological_cutvar].Data());
            vector<TH1D*> temp = Xi1530SysCheck(
                finputfile, fvarfile, TopologicalCutSystematic_bins[topological_cutvar].Data());
            cPlotCanvas_TopologyCut->cd();
            hTopologicalCutSys_fraction.push_back(temp[0]);
            hTopologicalCutSys_fraction[arraycount]->SetLineColor(FI_topocut + arraycount);
            hTopologicalCutSys_fraction[arraycount]->SetMarkerColor(FI_topocut + arraycount);
            hTopologicalCutSys_fraction[arraycount]->Draw("same");
            TopologicalCutSyslegend->AddEntry(hTopologicalCutSys_fraction[arraycount],
                                   TopologicalCutSystematic_bins[topological_cutvar], "L");
            for (int bin = 0; bin <= ptbin.size(); bin++) {
                hMultiplicityIndependentTopologicalCutSysSum->SetBinContent(
                    bin+1,
                    hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(bin+1) +
                        pow(hTopologicalCutSys_fraction[arraycount]->GetBinContent(bin+1),
                            2));
            }
            arraycount++;
        }
        TopologicalCutSyslegend->AddEntry(
            hMultiplicityIndependentTopologicalCutSysSum,
            Form("Cut Systematic sum (%2.f-%2.f)", multi_start, multi_end),
            "L");
        for (int bin = 0; bin < ptbin.size(); bin++) {
            hMultiplicityIndependentTopologicalCutSysSum->SetBinContent(bin+1,
                                      sqrt(hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(bin+1)));
            cout << "TEST, bin" << bin+1 << ", " << hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(bin+1)
                 << endl;
            TopolCutSysSum.push_back(hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(bin+1));
        }
        cout << "End of Check Sys Error for Topological Cut sys " << endl;
    }
    else{
        hMultiplicityIndependentTopologicalCutSysSum = (TH1D*)MakeHistfromArray("TopologicalCutSys", TopolCutSysSum, ptbin);
    }
    fout->cd();
    if (!(multi_end - multi_start > 99))
        TopologicalCutSyslegend->AddEntry(hMultiplicityIndependentTopologicalCutSysSum,
                               Form("Topological Cut Systematic sum (%2.f-%2.f)",
                                    multibin0100[0], multibin0100[1]),
                               "L");
    cPlotCanvas_TopologyCut->cd();
    TopologicalCutSyslegend->AddEntry(hStatError, "Stat Error", "F");
    if (smoothing) {
        for (int bin = 1; bin < ptbin.size() - 1; bin++) {
            if ((ptbin[bin] < 0.8) || (ptbin[bin] > 8.8))
                continue;
            double checkValue =
                hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(
                    bin + 1);
            double checkValue_before =
                hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(
                    bin);
            double checkValue_after =
                hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(
                    bin + 2);
            double smoothedValue =
                hMultiplicityIndependentTopologicalCutSysSum->GetBinContent(
                    bin + 1);
            double correctionfactor = 3;

            double averageOfNeighbor =
                (checkValue + checkValue_before + checkValue_after) /
                correctionfactor;

            if (checkValue > smoothingCriteron * averageOfNeighbor)
                smoothedValue = averageOfNeighbor;
            // only if the value is significantly larger than neighbor.
            hMultiplicityIndependentTopologicalCutSysSum->SetBinContent(
                bin + 1, smoothedValue);
        }
    }
    hMultiplicityIndependentTopologicalCutSysSum->SetLineColor(kBlack);
    hMultiplicityIndependentTopologicalCutSysSum->SetLineStyle(7);
    hMultiplicityIndependentTopologicalCutSysSum->Draw("same");
    hMultiplicityIndependentTopologicalCutSysSum->Write("TopologicalCutSysSum");

    TopologicalCutSyslegend->Draw();
    cPlotCanvas_TopologyCut->Write(Form("TopologicalCutVar(%2.f-%2.f)", multi_start, multi_end));
    cPlotCanvas_TopologyCut->SaveAs(
        Form("./TopologicalCutSysFraction(%2.f-%2.f).pdf", multi_start, multi_end));

    totalsystematic.push_back(hMultiplicityIndependentTopologicalCutSysSum);

    // PID Systematics
    TCanvas* cPlotCanvas_PIDCut =
        new TCanvas("cPlotCanvas_PIDCut", "cPlotCanvas_PIDCut", w, h);
    cPlotCanvas_PIDCut->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // cTestCanvas->SetLogy(true);
    cPlotCanvas_PIDCut->Draw();
    cPlotCanvas_PIDCut->cd();
    hStatError->Draw("BAR");
    hStatError->SetMaximum(0.3);
    vector<TString> PIDCutSystematic_bins = {"TPCNsigmaXi1530PionLoose",
                                          "TPCNsigmaXi1530PionTight",
                                          "TPCNsigmaXiLoose",
                                          "TPCNsigmaXiTight"};
    vector<TString> PIDCutSystematic_notuse = {
        "TPCNsigmaXiLoose",
        "TPCNsigmaXi1530PionLoose"};
    vector<TH1D*> hPIDCutSys_fraction;
    auto PIDCutSyslegend = new TLegend(0.4, 0.6, 0.9, 0.88);
    PIDCutSyslegend->SetNColumns(2);
    PIDCutSyslegend->SetBorderSize(0);
    PIDCutSyslegend->SetFillStyle(0);
    PIDCutSyslegend->SetFillStyle(0);
    TH1D* hMultiplicityIndependentPIDCutSysSum;

    for (int bin = 0; bin < ptbin.size(); bin++) {
        sigExMaxError[bin] = 0;
    }
    if ((fUseCutSysfrom0100) && (multi_end - multi_start > 99)){
        cout << "Check Sys Error for Cut sys " << endl;
        hMultiplicityIndependentPIDCutSysSum = (TH1D*)MakeHistfromArray("PIDCutSys", ErrorArray, ptbin);
        Int_t FI_PIDcut = TColor::CreateGradientColorTable(
            NRGBs, stops, red, green, blue, PIDCutSystematic_bins.size());
        int arraycount = 0;
        for (int PID_cutvar = 0; PID_cutvar < PIDCutSystematic_bins.size(); PID_cutvar++) {
            // if (PID_cutvar == 2) continue;
            if (std::find(PIDCutSystematic_notuse.begin(),
                          PIDCutSystematic_notuse.end(),
                          PIDCutSystematic_bins[PID_cutvar]) !=
                PIDCutSystematic_notuse.end())
                continue;
            fvarfile =
                Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
                     path.Data(), PID_cutvar + 2, multibin[0], multibin[1],
                     PIDCutSystematic_bins[PID_cutvar].Data());
            vector<TH1D*> temp = Xi1530SysCheck(
                finputfile, fvarfile, PIDCutSystematic_bins[PID_cutvar].Data());
            cPlotCanvas_PIDCut->cd();
            hPIDCutSys_fraction.push_back(temp[0]);
            hPIDCutSys_fraction[arraycount]->SetLineColor(FI_PIDcut + arraycount);
            hPIDCutSys_fraction[arraycount]->SetMarkerColor(FI_PIDcut + arraycount);
            hPIDCutSys_fraction[arraycount]->Draw("same");
            PIDCutSyslegend->AddEntry(hPIDCutSys_fraction[arraycount],
                                   PIDCutSystematic_bins[PID_cutvar], "L");
            for (int bin = 0; bin <= ptbin.size(); bin++) {
                if (hPIDCutSys_fraction[arraycount]->GetBinContent(bin + 1) >
                    sigExMaxError[bin])
                    sigExMaxError[bin] =
                        hPIDCutSys_fraction[arraycount]->GetBinContent(bin + 1);
                hMultiplicityIndependentPIDCutSysSum->SetBinContent(
                    bin+1,
                    hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin+1) +
                        pow(hPIDCutSys_fraction[arraycount]->GetBinContent(bin+1),
                            2));
            }
            arraycount++;
        }
        PIDCutSyslegend->AddEntry(
            hMultiplicityIndependentPIDCutSysSum,
            Form("PID Cut Systematic sum (%2.f-%2.f)", multi_start, multi_end),
            "L");
        /*
    // Quadrature sum
    for (int bin = 0; bin < ptbin.size(); bin++) {
            hMultiplicityIndependentPIDCutSysSum->SetBinContent(bin+1,
                                      sqrt(hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin+1)));
            cout << "TEST, bin" << bin+1 << ", " <<
    hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin+1)
                 << endl;
            PIDCutSysSum.push_back(hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin+1));
        }
    */
        for (int bin = 0; bin < ptbin.size(); bin++) {
            hMultiplicityIndependentPIDCutSysSum->SetBinContent(
                bin + 1, sigExMaxError[bin]);
            cout << "Check bin(final) " << bin + 1 << " : "
                 << hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin + 1)
                 << endl;
            PIDCutSysSum.push_back(
                hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin + 1));
        }
        
        cout << "End of Check Sys Error for PID Cut sys " << endl;
    }
    else{
        hMultiplicityIndependentPIDCutSysSum = (TH1D*)MakeHistfromArray("PIDCutSys", PIDCutSysSum, ptbin);
    }
    fout->cd();
    if (!(multi_end - multi_start > 99))
        PIDCutSyslegend->AddEntry(hMultiplicityIndependentPIDCutSysSum,
                               Form("PID Cut Systematic sum (%2.f-%2.f)",
                                    multibin0100[0], multibin0100[1]),
                               "L");
    cPlotCanvas_PIDCut->cd();
    if (smoothing) {
        for (int bin = 1; bin < ptbin.size() - 1; bin++) {
            if ((ptbin[bin] < 0.8) || (ptbin[bin] > 8.8))
                continue;
            double checkValue =
                hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin + 1);
            double checkValue_before =
                hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin);
            double checkValue_after =
                hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin + 2);
            double smoothedValue =
                hMultiplicityIndependentPIDCutSysSum->GetBinContent(bin + 1);
            double correctionfactor = 3;
            

            double averageOfNeighbor =
                (checkValue + checkValue_before + checkValue_after) /
                correctionfactor;

            if (checkValue > smoothingCriteron * averageOfNeighbor)
                smoothedValue = averageOfNeighbor;
            // only if the value is significantly larger than neighbor.
            hMultiplicityIndependentPIDCutSysSum->SetBinContent(bin + 1,
                                                                smoothedValue);
        }
    }
    PIDCutSyslegend->AddEntry(hStatError, "Stat Error", "F");
    hMultiplicityIndependentPIDCutSysSum->SetLineColor(kBlack);
    hMultiplicityIndependentPIDCutSysSum->SetLineStyle(7);
    hMultiplicityIndependentPIDCutSysSum->Draw("same");
    hMultiplicityIndependentPIDCutSysSum->Write("PIDCutSysSum");

    PIDCutSyslegend->Draw();
    cPlotCanvas_PIDCut->Write(Form("PIDCutVar(%2.f-%2.f)", multi_start, multi_end));
    cPlotCanvas_PIDCut->SaveAs(
        Form("./PIDCutSysFraction(%2.f-%2.f).pdf", multi_start, multi_end));

    totalsystematic.push_back(hMultiplicityIndependentPIDCutSysSum);

    // Other Systematics
    // 1. Material Budget (Xi+-)
    // Study of strangeness production as a function of the charged particle
    // multiplicity in pp collisions at \sqrt{s} = 13 TeV
    // https://alice-notes.web.cern.ch/node/478
    // 4%, pT independent, due to the lack of knowledge of the material budget.
    vector<double> materialbudget;
    for (int bin = 0; bin <= ptbin.size(); bin++)
        materialbudget.push_back(0.04);
    auto hSysMaterialBudget = (TH1D*)MakeHistfromArray("hSysMaterialBudget", materialbudget, ptbin);
    // 2. Mult. indipendent efficiencies
    // Study of strangeness production as a function of the charged particle
    // multiplicity in pp collisions at \sqrt{s} = 13 TeV
    // https://alice-notes.web.cern.ch/node/478
    // 2%, pT independent, due to the computation of the efficiencies in the
    // integrated case over multi-plicity
    vector<double> multiindepeffi;
    for (int bin = 0; bin <= ptbin.size(); bin++)
        multiindepeffi.push_back(0.02);
    auto hSysMultiIndepEffi = (TH1D*)MakeHistfromArray("hSysMultiIndepEffi", multiindepeffi, ptbin);
    // 3. Tracking efficiencies
    // https://twiki.cern.ch/twiki/bin/view/ALICE/TrackingEfficiencyCharged
    // 3%, pT independent
    vector<double> trackingeffi;
    for (int bin = 0; bin <= ptbin.size(); bin++)
        trackingeffi.push_back(0.03);
    auto hSysTrackingEffi = (TH1D*)MakeHistfromArray(
        "hSysTrackingEffi", trackingeffi, ptbin);
    //------------------------------------------------------------------------------
    // Total Uncertainty
    TCanvas* cPlotCanvas_TotalUncertainty =
        new TCanvas("cPlotCanvas_TotalUncertainty", "cPlotCanvas_TotalUncertainty", w, h);
    cPlotCanvas_TotalUncertainty->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cPlotCanvas_TotalUncertainty->Draw();
    cPlotCanvas_TotalUncertainty->cd();
    hStatError->Draw("BAR");
    hStatError->SetMaximum(0.3);
    hStatError->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    auto TotalSyslegend = new TLegend(0.6, 0.6, 0.9, 0.88);
    TotalSyslegend->SetNColumns(1);
    TotalSyslegend->SetBorderSize(0);
    TotalSyslegend->SetFillStyle(0);

    Int_t FI_topocut = TColor::CreateGradientColorTable(
            NRGBs, stops, red, green, blue, 6);
    // 1. Fit Sys sum
    hFitSysSum->SetMaximum(0.3);
    hFitSysSum->SetLineColor(FI_topocut);
    hFitSysSum->SetLineStyle(1);
    hFitSysSum->SetBinContent(1,0);
    hFitSysSum->SetBinContent(11,0);
    hFitSysSum->SetAxisRange(0, 8.8, "X");
    hFitSysSum->Draw("same");
    TotalSyslegend->AddEntry(hFitSysSum,"Fit Variation(Signal Extraction)","L");
    // 2. Topological Sys sum
    hMultiplicityIndependentTopologicalCutSysSum->SetLineColor(FI_topocut + 1);
    hMultiplicityIndependentTopologicalCutSysSum->SetLineStyle(1);
    hMultiplicityIndependentTopologicalCutSysSum->SetBinContent(1,0);
    hMultiplicityIndependentTopologicalCutSysSum->SetBinContent(11,0);
    hMultiplicityIndependentTopologicalCutSysSum->Draw("same");
    TotalSyslegend->AddEntry(hMultiplicityIndependentTopologicalCutSysSum,"Topological Cut variation","L");
    // 3. PID Sys sum
    hMultiplicityIndependentPIDCutSysSum->SetLineColor(FI_topocut + 2);
    hMultiplicityIndependentPIDCutSysSum->SetLineStyle(1);
    hMultiplicityIndependentPIDCutSysSum->SetBinContent(1,0);
    hMultiplicityIndependentPIDCutSysSum->SetBinContent(11,0);
    hMultiplicityIndependentPIDCutSysSum->Draw("same");
    TotalSyslegend->AddEntry(hMultiplicityIndependentPIDCutSysSum,"PID Cut variation","L");
    // 4. Material Budget
    hSysMaterialBudget->SetLineColor(FI_topocut + 3);
    hSysMaterialBudget->SetLineStyle(1);
    hSysMaterialBudget->SetBinContent(1,0);
    hSysMaterialBudget->SetBinContent(11,0);
    hSysMaterialBudget->Draw("same");
    TotalSyslegend->AddEntry(hSysMaterialBudget,"Material Budget","L");
    // 5. Mult. indipendent efficiencies
    hSysMultiIndepEffi->SetLineColor(FI_topocut + 4);
    hSysMultiIndepEffi->SetLineStyle(2);
    hSysMultiIndepEffi->SetBinContent(1,0);
    hSysMultiIndepEffi->SetBinContent(11,0);
    hSysMultiIndepEffi->Draw("same");
    TotalSyslegend->AddEntry(hSysMultiIndepEffi,"Mult. indipendent efficiencies","L");
    // 5. Mult. indipendent efficiencies
    hSysTrackingEffi->SetLineColor(FI_topocut + 5);
    hSysTrackingEffi->SetLineStyle(2);
    hSysTrackingEffi->SetBinContent(1, 0);
    hSysTrackingEffi->SetBinContent(11, 0);
    hSysTrackingEffi->Draw("same");
    TotalSyslegend->AddEntry(hSysTrackingEffi, "ITS-TPC tracking efficiency",
                             "L");

    if(!correl_error_skip){
        totalsystematic.push_back(hSysMaterialBudget);
        totalsystematic.push_back(hSysMultiIndepEffi);
        totalsystematic.push_back(hSysTrackingEffi);
    }

    auto base_state = (TH1D*)base->Clone();
    for (int bin = 0; bin <= ptbin.size(); bin++) {
        base->SetBinError(bin + 1, 0);
    }

    vector<double> totalSysError;
    for (int bin = 0; bin < ptbin.size(); bin++) {
        //double TotError = pow(base->GetBinError(bin + 1), 2);  // Stat
        double TotError = 0.;  // Non Stat
        cout << "Spectrum ptbin:" << bin << "  "
             << base->GetBinContent(bin + 1);

        for (int syserrorbin = 0; syserrorbin < totalsystematic.size();
             syserrorbin++) {
            TotError +=
                pow(totalsystematic[syserrorbin]->GetBinContent(bin + 1) *
                        base->GetBinContent(bin + 1),
                    2);
        }
        base->SetBinError(bin + 1, sqrt(TotError));
        cout << "  +- " << base->GetBinError(bin + 1) << "(sys)" << endl;
    }
    for (int bin = 0; bin < ptbin.size(); bin++) {
        double temp = 0.;
        for (int syserrorbin = 0; syserrorbin < totalsystematic.size();
             syserrorbin++) {
            temp +=
                pow(totalsystematic[syserrorbin]->GetBinContent(bin + 1), 2);
        }
        totalSysError.push_back(sqrt(temp));
    }
    auto hSysErrorSum = (TH1D*)MakeHistfromArray("hSysErrorSum", totalSysError, ptbin);
    hSysErrorSum->SetLineColor(kBlack);
    hSysErrorSum->SetLineStyle(2);
    TotalSyslegend->AddEntry(hSysErrorSum,"Total Error","L");
    TotalSyslegend->AddEntry(hStatError, "Stat Error", "F");
    hSysErrorSum->Draw("same");
    TotalSyslegend->Draw("same");
    cPlotCanvas_TotalUncertainty->SaveAs(Form("./TotalSysError(%2.f-%2.f).pdf", multi_start, multi_end));

    TCanvas* cResultCanvas =
        new TCanvas("cResultCanvas", "cResultCanvas", w, h);
    cResultCanvas->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cResultCanvas->Draw();
    cResultCanvas->cd();
    cResultCanvas->SetLogy(true);
    fout->cd();

    base->Draw("E2");  // SYS ERROR
    base->SetMarkerSize(0.);
    base->SetFillStyle(0);
    base->SetMinimum(1e-7);
    base->SetMaximum(5e-2);
    base->SetFillColor(kBlue);
    base->SetLineColor(kBlue);
    base->GetYaxis()->SetTitle("1/N_{event}d^{2}N/(dydp_{T}) (GeV/c)^{-1}");
    base->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    //base->GetXaxis()->SetRangeUser(0.8, 7.4);
    base->Write(Form("hSpectra_SYSe(%2.f-%2.f)", multi_start, multi_end));

    base_state->Draw("E1 SAME");  // stat ERROR
    base_state->SetMarkerSize(0.);
    // base_state->SetFillStyle(3003);
    base_state->SetMarkerColor(kBlue);
    base_state->SetLineColor(kBlue);
    base_state->SetMinimum(1e-7);
    base_state->SetMaximum(5e-2);
    base_state->GetYaxis()->SetTitle(
        "1/N_{event}d^{2}N/(dydp_{T}) (GeV/c)^{-1}");
    base_state->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    //base_state->GetXaxis()->SetRangeUser(0.8, 7.4);
    base_state->Write(
        Form("hSpectra_STATe(%2.f-%2.f)", multi_start, multi_end));

    //-------------------------------------------------------------------------------

    TCanvas* cResultCanvas_gr =
        new TCanvas("cResultCanvas_gr", "cResultCanvas_gr", w, h);
    cResultCanvas_gr->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cResultCanvas_gr->Draw();
    cResultCanvas_gr->cd();
    cResultCanvas_gr->SetLogy(true);

    base->SetBinContent(1, 0);
    base->SetBinContent(ptbin.size() - 1, 0);
    base_state->SetBinContent(1, 0);
    base_state->SetBinContent(ptbin.size() - 1, 0);
    // auto htest =
    // (TGraphErrors*)AliPWGHistoTools::GetGraphFromHisto((TH1D*)base); auto
    // htest2
    // =(TGraphErrors*)AliPWGHistoTools::GetGraphFromHisto((TH1D*)base_state);

    TGraphAsymmErrors* htest = new TGraphAsymmErrors(base);
    htest->SetMinimum(1e-7);
    htest->SetMaximum(5e-2);
    htest->GetYaxis()->SetTitle("1/N_{event}d^{2}N/(dydp_{T}) (GeV/c)^{-1}");
    htest->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    htest->GetXaxis()->SetRangeUser(0.5, 10);
    TGraphAsymmErrors* htest2 = new TGraphAsymmErrors(base_state);

    htest->SetMarkerColor(kBlack);
    htest->SetFillColor(0);
    htest->SetLineColor(kBlack);
    htest2->SetMarkerColor(kBlack);
    htest2->SetFillColor(0);
    htest2->SetLineColor(kBlack);

    htest->Draw("A5");
    htest2->Draw("P same");

    // t->DrawLatex(0.42, 0.92, "ALICE Preliminary, #bf{pp #sqrt{#it{s}} = 13
    // TeV}");
    t->DrawLatex(0.6, 0.82, "#bf{ALICE, pp, #sqrt{#it{s}} = 13 TeV}");
    t->DrawLatex(0.66, 0.75, "#bf{|#it{y}| < 0.5, INEL > 0}");
    t->DrawLatex(0.13, 0.15, "#bf{Uncertainties: stat.(bars), syst.(boxex)}");
    t_big->DrawLatex(0.58, 0.66, "1/2(#Xi(1530)^{0} + cc)");

    htest->Write(Form("hSpectra_sys(%2.f-%2.f)", multi_start, multi_end));
    htest2->Write(Form("hSpectra_stat(%2.f-%2.f)", multi_start, multi_end));
    cResultCanvas_gr->SaveAs(
        Form("hSpectra(%2.f-%2.f).pdf", multi_start, multi_end));

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
        "7TeV Spectra with systematic error", CorrectedYeild_7TeV, ptbin, CorrectedYeild_syserr_7TeV);
    auto hSpectra_7TeV_staterr = MakeHistfromArray(
        "7TeV Spectra with statistical error", CorrectedYeild_7TeV, ptbin, CorrectedYeild_staterr_7TeV);

    TGraphAsymmErrors* h7TeV_sys = new TGraphAsymmErrors(hSpectra_7TeV_syserr);
    h7TeV_sys->SetMinimum(1e-7);
    h7TeV_sys->SetMaximum(5e-2);
    //h7TeV_sys->GetXaxis()->SetRangeUser(0.8, 7.4);
    TGraphAsymmErrors* h7TeV_stat =
        new TGraphAsymmErrors(hSpectra_7TeV_staterr);

    h7TeV_sys->SetMarkerColor(kBlue);
    h7TeV_sys->SetFillColorAlpha(kBlue, 0.00001);
    h7TeV_sys->SetLineColor(kBlue);
    h7TeV_stat->SetMarkerColor(kBlue);
    h7TeV_stat->SetFillColorAlpha(kBlue, 0.00001);
    h7TeV_stat->SetLineColor(kBlue);

    h7TeV_sys->Draw("5 same");
    h7TeV_stat->Draw("P same");
    h7TeV_sys->Write("h7TeV_sys");
    h7TeV_stat->Write("h7TeV_stat");
    cResultCanvas_gr->SaveAs(
        Form("hSpectra_with7TeV(%2.f-%2.f).pdf", multi_start, multi_end));

    //returnHist.push_back(base);
    //returnHist.push_back(base_state);

    fout->Close();
    //return returnHist;
    return base;
}
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

vector<TH1D*> Xi1530SysCheck(TString fdeafultinputfile,
                             TString fvariationinputfile,
                             char const* options) {
    // Return TH1D* Histogram vector
    //   0. ErrorFraction (1 - variation/origin) spectra
    //   1. Ratio pT spectra between origin and variation
    //   2. Error spectra
    //   3. StatError of variation spectra

    bool fastbarlowcheck = true; // true: will check barlow cut every systematics/pT bin.

    cout << "input1: " << fdeafultinputfile.Data() << endl;
    cout << "input2: " << fvariationinputfile.Data() << endl;
    // Output file
    TFile* fROOTout;
    vector<TH1D*> hReturn;
    if (save)
        fROOTout = new TFile(Form("./Var_%s.root", options), "RECREATE");

    TCanvas* cTestCanvas = new TCanvas("cTestCanvas", "cTestCanvas", w, h);
    cTestCanvas->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // cTestCanvas->SetLogy(true);
    cTestCanvas->Draw();

    // inputs
    TFile* fDefault = new TFile(fdeafultinputfile.Data());
    TFile* fVariation = new TFile(fvariationinputfile.Data());
    auto base = (TH1D*)fDefault->Get("hXiSpectrum");
    auto variation = (TH1D*)fVariation->Get("hXiSpectrum");

    if (save) {
        fROOTout->cd();
        base->Write("Default");
        variation->Write(Form("%s", options));
    }

    // Let's calcuate!!
    base->SetMarkerSize(1.);
    int nbins = base->GetNbinsX();

    vector<double> StatError;  // from base
    vector<double> StatError_var;  // from varation
    vector<double> VarError;
    vector<double> VarErrorCheck;
    vector<double> VarErrorFraction;

    // init StatError Array
    for (int bin = 0; bin < nbins; bin++){
        StatError.push_back(base->GetBinError(bin + 1));
        StatError_var.push_back(variation->GetBinError(bin + 1));
    }
        
    for (int bin = 1; bin < nbins+1; bin++) {  // pt bins
        double error =  // abs(Variation - Default) -> "Error"
            abs(variation->GetBinContent(bin) - base->GetBinContent(bin));
        double errorFraction =  // abs(Variation - Default)/Deafult - >"Error
                                // fraction"
            error / variation->GetBinContent(bin);
        double deltasigma =  // sqrt(abs(variation_stat.error^2 -
                             // default_stat.error^2)) -> "delta sigma"
            sqrt(abs(pow(variation->GetBinError(bin), 2) -
                     pow(base->GetBinError(bin), 2)));
        double deltasigmafraction = deltasigma / variation->GetBinContent(bin);
        VarErrorCheck.push_back(variation->GetBinContent(bin) -
                                base->GetBinContent(bin));
        if (fastbarlowcheck) {  // Barlow check in dhevan case!
            cout << "FAST Barlow check in bin " << ptbin[bin-1] << " - "
                 << ptbin[bin] << endl;
            cout << "Variation value: " << variation->GetBinContent(bin)
                 << ", Default value: " << base->GetBinContent(bin)
                 << ", Variation error: " << variation->GetBinError(bin)
                 << ", Default error: " << base->GetBinError(bin)
                 << endl;
            if (errorFraction / deltasigmafraction > 1) {
                cout << "-> PASS! errorfraction: " << errorFraction
                     << ", deltasigmafraction:" << deltasigmafraction
                     << ", check: " << errorFraction / deltasigmafraction
                     << endl;
                VarError.push_back(error);
                VarErrorFraction.push_back(errorFraction);
            } else {
                cout << "-> FAIL! errorfraction: " << errorFraction
                     << ", deltasigmafraction:" << deltasigmafraction
                     << ", check: " << errorFraction / deltasigmafraction
                     << endl;
                VarError.push_back(zero);
                VarErrorFraction.push_back(zero);
            }
        }
        else{ // or, just check them later! stat error will be returned as a third member of array.
            VarError.push_back(error);
            VarErrorFraction.push_back(errorFraction);
        }
    }

    // Drawing QA Plots
    auto hSystematics =
        (TH1D*)MakeHistfromArray(Form("Var_%s", options), VarErrorFraction, ptbin);
    hSystematics->SetFillStyle(0);
    hSystematics->SetMinimum(0);
    hSystematics->SetMaximum(0.3);
    hSystematics->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hSystematics->GetYaxis()->SetTitle("Uncertainty(fraction)");
    hSystematics->Draw();

    auto hVarStatError =
        (TH1D*)MakeHistfromArray(Form("VarStat_%s", options), StatError_var, ptbin);
    auto hVarError =
        (TH1D*)MakeHistfromArray(Form("VarStat_%s", options), VarErrorCheck, ptbin);

    auto ratio = (TH1D*)variation->Clone();
    ratio->Divide(base);
    ratio->SetMarkerSize(0.);
    ratio->SetBinContent(1, zero);
    ratio->SetBinError(1, zero);
    ratio->SetBinContent(nbins, zero);
    ratio->SetBinError(nbins, zero);
    ratio->SetMinimum(0.5);
    ratio->SetMaximum(1.5);
    double totalerrorsum = 0.;
    for (int bin = 1; bin < nbins; bin++) {  // pt bins
        totalerrorsum += ratio->GetBinContent(bin);
    }
    cout << "TOTAL Error Sum from " << options << " : " << totalerrorsum
         << endl;
    if (save) {
        hSystematics->Write(Form("Var_%s", options));
        ratio->Write(Form("ratio_%s", options));
        fROOTout->Close();
    }
    hReturn.push_back(hSystematics);
    hReturn.push_back(ratio);
    hReturn.push_back(hVarError);
    hReturn.push_back(hVarStatError);
    return hReturn;
}
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
    bool topPlotLogY = 1;                 // 0 = no log; 1= log
    bool bottomPlotLogY = 1;              // 0 = no log; 1= log
    TString yTitle2 = "Ratio to no W.D.Vertexer";  // bottom plot y axis title
    if (options.Contains("DIFF"))
        yTitle2 = "Ratio to 7 TeV Results";
    if (options.Contains("FINAL"))
        yTitle2 = "Ratio to INEL>0";
    if (options.Contains("NLOG"))
        bottomPlotLogY = 0;

    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    Int_t FIh = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 numberOfNumeratorHists);
    
    vector<int> histColors;
    vcolors = {};
    for (int color = 0; color < numberOfNumeratorHists; color++) {
        // histColors.push_back(FIh+color);
        // vcolors.push_front(FIh+color);
        histColors.insert(histColors.begin(), FIh + color);
        vcolors.insert(vcolors.begin(), FIh + color);
    }
    // histColors.push_back(kBlue); // change colors as you like
    // histColors.push_back(kRed);
    // histColors.push_back(kGreen - 1);

    int histDenominatorColor = kBlue;

    float defaultRatioYmin = 4e-2;
    float defaultRatioYmax = 2e1;
    if (options.Contains("Efficiency")) {
        defaultRatioYmin = 0.7;
        defaultRatioYmax = 1.5;
    }
    if (options.Contains("7T")) {
        defaultRatioYmin = 0.9;
        defaultRatioYmax = 3.0;
    }
    // END of Variables
    //*************************************************

    TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 850, 1150);
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
    c1_1->SetBottomMargin(0.2);
    c1_1->SetRightMargin(0.01);
    c1_1->SetFillStyle(0);
    c1_1->SetLogy(bottomPlotLogY);

    hist_over_denomHist[0]->Draw("E2");
    hist_over_denomHist[0]->SetLineWidth(1);
    hist_over_denomHist[0]->SetLineColor(histColors[0]);
    hist_over_denomHist[0]->SetFillColor(histColors[0]);
    hist_over_denomHist[0]->SetMarkerColor(histColors[0]);
    hist_over_denomHist[0]->SetMarkerStyle(Markerset[0]);
    hist_over_denomHist[0]->SetMarkerSize(Markerset_size[0]);
    hist_over_denomHist[0]->SetMinimum(defaultRatioYmin);
    hist_over_denomHist[0]->SetMaximum(defaultRatioYmax);
    hist_over_denomHist[0]->GetYaxis()->SetNdivisions(5);
    hist_over_denomHist[0]->SetTitle(";" + xTitle + ";" + yTitle2);
    hist_over_denomHist[0]->GetXaxis()->SetTitleSize(0.07);
    hist_over_denomHist[0]->GetXaxis()->SetLabelSize(0.07);
    hist_over_denomHist[0]->GetYaxis()->SetLabelSize(0.07);
    hist_over_denomHist[0]->GetYaxis()->SetTitleSize(0.07);
    hist_over_denomHist[0]->GetYaxis()->SetTitleOffset(0.6);
    if (options.Contains("FINAL"))
        hist_over_denomHist[0]->GetXaxis()->SetRangeUser(0.5, 8.8);
    for (int i = 1; i < numberOfNumeratorHists; i++) {
        hist_over_denomHist[i]->SetLineWidth(1);
        hist_over_denomHist[i]->SetLineColor(histColors[i]);
        hist_over_denomHist[i]->SetFillColor(histColors[i]);
        hist_over_denomHist[i]->SetMarkerColor(histColors[i]);
        hist_over_denomHist[i]->SetMarkerStyle(Markerset[i]);
        hist_over_denomHist[i]->SetMarkerSize(Markerset_size[i]);
        hist_over_denomHist[i]->Draw("E2 same");
    }
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
    c1_2->SetBottomMargin(0.01);
    c1_2->SetRightMargin(0.01);
    c1_1->SetFillStyle(0);
    
    denominatorHistogram->SetLineWidth(1);
    if (options.Contains("FINAL"))
        histDenominatorColor = 0;
    denominatorHistogram->SetLineColor(histDenominatorColor);
    denominatorHistogram->SetMarkerColor(histDenominatorColor);
    denominatorHistogram->Draw("E2");
    denominatorHistogram->Draw("SAME");
    denominatorHistogram->SetLabelSize(0.0);
    denominatorHistogram->GetXaxis()->SetTitleSize(0.00);
    denominatorHistogram->GetYaxis()->SetLabelSize(0.03);
    denominatorHistogram->GetYaxis()->SetTitleSize(0.04);
    denominatorHistogram->GetYaxis()->SetTitleOffset(1.2);
    denominatorHistogram->SetTitle(title + ";;" + yTitle);
    //if (options.Contains("FINAL")) denominatorHistogram->GetXaxis()->SetRangeUser(0.5, 8.8);
    denominatorHistogram->SetMaximum(2e-1);
    if (options.Contains("FINAL"))
        denominatorHistogram->SetMaximum(2);
    if (options.Contains("LOW"))
        denominatorHistogram->SetMinimum(5e-8);

    for (int i = 0; i < numberOfNumeratorHists; i++) {
        hists[i]->SetLineWidth(1);
        if (options.Contains("POW"))
            hists[i]->Scale(pow(2, numberOfNumeratorHists - i));
        hists[i]->SetLineColor(histColors[i]);
        hists[i]->SetFillColor(histColors[i]);
        hists[i]->SetMarkerColor(histColors[i]);
        hists[i]->SetMarkerStyle(Markerset[i]);
        hists[i]->SetMarkerSize(Markerset_size[i]);
        hists[i]->Draw("same");
        hists[i]->Draw("E2 same");
    }

    c1_2->SetLogy(topPlotLogY);
    // End bottom plot
    //*************************************************

    return c1;
}
TH1D* Xi1530SysFitVar(TString nameFitVar = "FitVar",
                     double multi_start = 0,
                     double multi_end = 100) {
    vector<TString> FitSysVar_bins;
    if (nameFitVar.Contains("BinCount")) {
        FitSysVar_bins = {"BinCountL", "BinCountR", "BinCountBoth"};
    }
    if (nameFitVar.Contains("FitVar")) {
        FitSysVar_bins = {"FitVarL", "FitVarR", "FitVarBoth"};
    }
    if (nameFitVar.Contains("NormVar")) {
        FitSysVar_bins = {"NormVar", "NormVarL", "NormVarR"};
    }
    multibin[0] = multi_start;
    multibin[1] = multi_end;

    legeondPlotOut->SetNColumns(2);
    legeondPlotOut->SetBorderSize(0);
    legeondPlotOut->SetFillStyle(0);
    TCanvas* cPlotbin = new TCanvas("cPlotbin", "cPlotbin", w, h);
    cPlotbin->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cPlotbin->Draw();
    vector<TH1D*> fBinByBinError;
    //initialize
    for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
        ptbinStatErrorAvg[bin].clear();
    }

    linestyle = 0;
    GlobalFirst = true;
    bool firstloopr = true;
    for (int i = 0; i < FitSysVar_bins.size(); i++) {
        cout << "Check for " << FitSysVar_bins[i] << endl;
        linestyle++;
        if (firstloopr) {
            vector<TH1D*> temp = Xi1530SysFitVarCheck(FitSysVar_bins[i]);
            for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
                fBinByBinError.push_back(temp[bin]);
            }
            firstloopr = false;
        } else {
            vector<TH1D*> temp = Xi1530SysFitVarCheck(FitSysVar_bins[i]);
            for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
                auto temp2 = (TH1D*)temp.at(bin)->Clone();
                fBinByBinError[bin]->Add(temp2);
            }
        }
        cout << "Check size ===: " << ptbinStatErrorAvg[0].size() << endl;
    }
    
    
    // for memo, small
    TLatex* t = new TLatex();
    t->SetNDC();
    t->SetTextSize(0.04);

    // Default configure histotgram for barlow check
    TString fDefault =
        Form("%sAnalysisResults_Extracted_1_Multi_%.2f-%.2f_Default1.root",
             path.Data(), multibin[0], multibin[1]);
    TFile* Defaultfile = new TFile(fDefault.Data());
    TH1D* hbase = (TH1D*)Defaultfile->Get("hXiSpectrum");

    vector<double> ErrorArray2;
    for (int bin = 0; bin <= ptbin.size(); bin++)
        ErrorArray2.push_back(0);
    auto temphist =
        (TH1D*)MakeHistfromArray(Form("%s", nameFitVar.Data()), ErrorArray2,ptbin);

    // Average Statistical Error from each test
    vector<double> fstatErrorAvg;
    for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
        double tempstaterror = 0.;
        for (int chek = 0; chek < ptbinStatErrorAvg[bin].size(); chek++) {  // for eacbin
            tempstaterror += ptbinStatErrorAvg[bin][chek];
        }
        fstatErrorAvg.push_back(tempstaterror/ptbinStatErrorAvg[bin].size());
        //cout << "Error in test[" << nameFitVar.Data() <<"], bin(" << bin << "): "  << fstatErrorAvg[bin] << " from total " << ptbinStatErrorAvg[bin].size() << "entry, default error: " << hbase->GetBinError(bin + 1) << endl;
    }
    vector<double> stderror;
    for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
        cPlotbin->cd();
        fBinByBinError[bin]->Draw();
        double Delta = abs(fBinByBinError[bin]->GetMean() + fBinByBinError[bin]->GetStdDev()) * hbase->GetBinContent(bin + 1);
        double realsigma = fstatErrorAvg[bin];
        double sigmadiffsum = sqrt(abs(pow(realsigma, 2) - pow(hbase->GetBinError(bin + 1), 2)));
        double barlowfactor = Delta / sigmadiffsum;
        cout << "Mean: " << fBinByBinError[bin]->GetMean() << endl;
        cout << "Std: " << fBinByBinError[bin]->GetStdDev() << endl;
        cout << "Delta(Mean+Std)): " << Delta << endl;
        cout << "StatE(default): " << hbase->GetBinError(bin + 1) << endl;
        cout << "Error(Delta, deafult - variation): " << Delta << endl;
        cout << "Average sigma: " << realsigma << endl;
        cout << "Error sum(subtract): " << sigmadiffsum << endl;
        cout << "Barlow factor: " << barlowfactor << endl;
        t->DrawLatex(0.68, 0.77,
                     Form("#bf{Mean: %.3f}", fBinByBinError[bin]->GetMean()));
        t->DrawLatex(
            0.68, 0.72,
            Form("#bf{StdDev: %.3f}", fBinByBinError[bin]->GetStdDev()));
        t->DrawLatex(0.68, 0.67, Form("#bf{#Delta: %.3f (x10^{6})}", Delta*1e6));
        t->DrawLatex(0.68, 0.62, Form("#bf{#sigma: %.3f (x10^{6})}", sigmadiffsum*1e6));
        t->DrawLatex(0.68, 0.57,
                     Form("#bf{#Delta/#sigma: %.3f}", barlowfactor));
        cPlotbin->SaveAs(Form("./%s/FitVarbin%i(%2.f-%2.f).pdf",
                              nameFitVar.Data(), bin, multibin[0],
                              multibin[1]));
        /*
        if(barlowfactor > 1){
            temphist->SetBinContent(bin + 1,
                                    abs(fBinByBinError[bin]->GetMean())
                                   +fBinByBinError[bin]->GetStdDev());
        }
        else{
            temphist->SetBinContent(bin + 1, 0);
        }*/
        
        temphist->SetBinContent(bin + 1, fBinByBinError[bin]->GetStdDev());
        double temperror = fBinByBinError[bin]->GetStdDev() +
                           abs(fBinByBinError[bin]->GetMean());
        stderror.push_back(temperror);
    }
    auto hstderror =
        (TH1D*)MakeHistfromArray("StdError", ErrorArray2, ptbin);
    for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
        hstderror->SetBinContent(bin + 1, 1);
        hstderror->SetBinError(bin + 1, stderror[bin]);
    }
    hstderror->SetFillColorAlpha(kBlack, 0.2);
    hstderror->SetLineColorAlpha(kBlack, 0.2);
    legeondPlotOut->AddEntry(hstderror, "Systematic Error(Stdev)", "F");
    cPlotOut->cd();
    hstderror->Draw("E3 Same");
    legeondPlotOut->Draw();
    cPlotOut->SaveAs(Form("./%s/FullVar(%.2f-%.2f).pdf", nameFitVar.Data(),
                          multibin[0], multibin[1]));
    legeondPlotOut->Clear();
    cPlotOut->Clear();
    return temphist;
}

vector<TH1D*> Xi1530SysFitVarCheck(TString nameFitVar) {
    bool plotit = false;
    bool savecheck = true;

    bool furthercheck = false;
    int nBins = 8;  // 8 default
    int variationrange_m = 4;
    int variationrange_p = 4;
    double Ymax = 1.1;
    double Ymin = 0.9;

    TString outputfolder = "out";

    if (nameFitVar.Contains("NormVar")) {
        outputfolder = "NormVar";
        nBins = 6;
        Ymax = 1.02;
        Ymin = 0.98;
        variationrange_m = 4;
        variationrange_p = 4;
    }
    if (nameFitVar.Contains("BinCount")) {
        outputfolder = "BinCount";
        Ymax = 1.1;
        Ymin = 0.9;
    }
    if (nameFitVar.Contains("BinCountL")) {
        nBins = 12;
        variationrange_m = 10;
        variationrange_p = 3;
    }
    if (nameFitVar.Contains("BinCountR")) {
        nBins = 12;
        variationrange_m = 3;
        variationrange_p = 10;
    }
    if (nameFitVar.Contains("BinCountBoth")) {
        nBins = 12;
        variationrange_m = 3;
        variationrange_p = 10;
    }

    if (nameFitVar.Contains("FitVar")) {
        outputfolder = "FitVar";
        nBins = 6;
        Ymax = 1.2;
        Ymin = 0.8;
        variationrange_m = 4;
        variationrange_p = 4;
    }

    // vector<TH1D*> returnHist;

    TString nameFitVar_p = nameFitVar + "p";
    TString nameFitVar_m = nameFitVar + "m";
    TCanvas* cPlot = new TCanvas("cPlot", "cPlot", w, h);
    cPlot->SetTickx();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cPlot->Draw();
    TF1* oneline = new TF1("oneline", "1", 0.5, 7.4);
    oneline->SetLineStyle(2);
    oneline->SetLineColor(kBlack);

    // Import base -> (TH1D*)base
    TString inputOptions = "Default";
    int cutbin = 1;
    TString finputfile = Form(
        "%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root", path.Data(),
        cutbin, multibin[0], multibin[1], inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1D* base = (TH1D*)inputfile->Get("hXiSpectrum");
    TH1D* base_bak = (TH1D*)base->Clone();
    // inputfile->Close();
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    // init Error Array
    vector<double> ErrorArray_fitvar;
    for (int bin = 0; bin <= ptbin.size(); bin++)
        ErrorArray_fitvar.push_back(0);

    const Int_t NRGBs = 5;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    Int_t FIf =
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, nBins);
    vector<TH1D*> hFitvarsys_fraction;

    TString ffitvarfile_p;
    TString ffitvarfile_pf;
    TString ffitvarfile_m;
    TString ffitvarfile_zero;

    vector<TH1D*> totalhistogram;

    auto FitVarSyslegend = new TLegend(0.13, 0.74, 0.5, 0.88);
    FitVarSyslegend->SetFillStyle(0);
    for (int variation = 1; variation < variationrange_m; variation++) {
        int vari_temp = variationrange_m - variation;
        cout << vari_temp << "in " << variationrange_m << endl;
        ffitvarfile_m =
            Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s%i.root",
                 path.Data(), cutbin, multibin[0], multibin[1],
                 nameFitVar_m.Data(), vari_temp);
        if (plotit)
            PlotXi1530(ffitvarfile_m);

        totalhistogram.push_back((TH1D*)GetSpectrafromName(ffitvarfile_m));

        vector<TH1D*> temp_m =
            Xi1530SysCheck(finputfile, ffitvarfile_m, nameFitVar_m.Data());
        temp_m[1]->SetLineColor(FIf + variation);
        temp_m[1]->SetMarkerColor(FIf + variation);
        temp_m[1]->SetLineStyle(linestyle);
        hFitvarsys_fraction.push_back(temp_m[1]);
        FitVarSyslegend->AddEntry(hFitvarsys_fraction[variation - 1],
                                  Form("%s_%i", nameFitVar.Data(), -vari_temp),
                                  "L");
        legeondPlotOut->AddEntry(hFitvarsys_fraction[variation - 1],
                                 Form("%s_%i", nameFitVar.Data(), -vari_temp),
                                 "L");
        hFitvarsys_fraction[variation - 1]->SetMaximum(Ymax);
        hFitvarsys_fraction[variation - 1]->SetMinimum(Ymin);
        for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
            ptbinStatErrorAvg[bin].push_back(temp_m[3]->GetBinContent(bin+1));
        }
    }
    int vectorsize = hFitvarsys_fraction.size();
    for (int variation = 1; variation < variationrange_p; variation++) {
        cout << variationrange_p + 3 << endl;
        cout << variationrange_p << "in " << variationrange_p << endl;
        ffitvarfile_p =
            Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s%i.root",
                 path.Data(), cutbin, multibin[0], multibin[1],
                 nameFitVar_p.Data(), variation);
        if (plotit)
            PlotXi1530(ffitvarfile_p);

        totalhistogram.push_back((TH1D*)GetSpectrafromName(ffitvarfile_p));

        vector<TH1D*> temp_p =
            Xi1530SysCheck(finputfile, ffitvarfile_p, nameFitVar_p.Data());
        temp_p[1]->SetLineColor(FIf + variation + vectorsize - 1);
        temp_p[1]->SetMarkerColor(FIf + variation + vectorsize - 1);
        temp_p[1]->SetLineStyle(linestyle);
        hFitvarsys_fraction.push_back(temp_p[1]);
        FitVarSyslegend->AddEntry(
            hFitvarsys_fraction[variation + vectorsize - 1],
            Form("%s_+%i", nameFitVar.Data(), variation), "L");
        legeondPlotOut->AddEntry(
            hFitvarsys_fraction[variation + vectorsize - 1],
            Form("%s_+%i", nameFitVar.Data(), variation), "L");
        hFitvarsys_fraction[variation + vectorsize - 1]->SetMaximum(Ymax);
        hFitvarsys_fraction[variation + vectorsize - 1]->SetMinimum(Ymin);

        for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
            ptbinStatErrorAvg[bin].push_back(temp_p[3]->GetBinContent(bin+1));
        }
    }
    if (furthercheck) {
        vectorsize = hFitvarsys_fraction.size();
        for (int variation = 1; variation < 5; variation++) {
            ffitvarfile_pf =
                Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s%i.root",
                     path.Data(), cutbin, multibin[0], multibin[1],
                     nameFitVar_p.Data(), variation + 6);
            if (plotit)
                PlotXi1530(ffitvarfile_pf);

            totalhistogram.push_back((TH1D*)GetSpectrafromName(ffitvarfile_pf));

            vector<TH1D*> temp_pf =
                Xi1530SysCheck(finputfile, ffitvarfile_pf, nameFitVar_p.Data());
            temp_pf[1]->SetLineColor(FIf + variation + vectorsize - 1);
            temp_pf[1]->SetMarkerColor(FIf + variation + vectorsize - 1);
            temp_pf[1]->SetLineStyle(linestyle);
            hFitvarsys_fraction.push_back(temp_pf[1]);
            FitVarSyslegend->AddEntry(
                hFitvarsys_fraction[variation + vectorsize - 1],
                Form("%s_+%i", nameFitVar.Data(), variation + 6), "L");
            legeondPlotOut->AddEntry(
                hFitvarsys_fraction[variation + vectorsize - 1],
                Form("%s_+%i", nameFitVar.Data(), variation + 6), "L");
            hFitvarsys_fraction[variation + vectorsize - 1]->SetMaximum(1.25);
            hFitvarsys_fraction[variation + vectorsize - 1]->SetMinimum(0.7);
            for (int bin = 0; bin < ptbin.size(); bin++) {  // for eacbin
                ptbinStatErrorAvg[bin].push_back(temp_pf[3]->GetBinContent(bin+1));
            }
        }
    }
    cPlot->cd();
    bool first = true;
    for (int variation = 0; variation < hFitvarsys_fraction.size();
         variation++) {
        cPlot->cd();
        hFitvarsys_fraction[variation]->SetAxisRange(0.8, 7.4, "X");
        hFitvarsys_fraction[variation]->GetYaxis()->SetTitle(
            "Uncertainty(Fraction)");
        if (first) {
            hFitvarsys_fraction[variation]->Draw("HIST L");
            oneline->Draw("same");
            first = false;
        } else
            hFitvarsys_fraction[variation]->Draw("same HIST L");
        cPlotOut->cd();
        if (GlobalFirst){
            hFitvarsys_fraction[variation]->Draw("HIST L");
            GlobalFirst = false;
        }
        else
            hFitvarsys_fraction[variation]->Draw("same HIST L");
    }
    cPlot->cd();
    FitVarSyslegend->Draw();
    if (savecheck)
        cPlot->SaveAs(Form("./%s/FitVarSysFraction_(%2.f-%2.f).pdf",
                           nameFitVar.Data(), multibin[0], multibin[1]));

    vector<TH1D*> fReturnHist;
    for (int bin = 1; bin <= ptbin.size(); bin++) {  // for eacbin
        TH1D* temphist;
        for (int variation = 0;
             variation < hFitvarsys_fraction.size();  // check it's tendency
             variation++) {
            auto tempdata =
                hFitvarsys_fraction[variation]->GetBinContent(bin) - 1;
            if (variation == 0)
                temphist =
                    new TH1D(Form("%i%i", variation, bin), "", 200, -0.5, 0.5);
            temphist->Fill(tempdata);
        }
        temphist->Draw();
        fReturnHist.push_back(temphist);
        if (savecheck)
            cPlot->SaveAs(Form("./%s/FitVarbin_%i_%s_(%2.f-%2.f).pdf",
                               outputfolder.Data(), bin, nameFitVar.Data(),
                               multibin[0], multibin[1]));
    }
    /*
    TCanvas* ratio_firvar = plotHistsAndRatio(
        totalhistogram, base_bak, "", "p_{T}(GeV/c)",
        "1/N_{event}d^{2}N/(dydp_{T}) (GeV/c)^{-1}", "SPECTRA NLOG FINAL CUTS");
    */
    TCanvas* ratio_firvar = new TCanvas("cCanvas", "cCanvas", w, h);
    // gStyle->SetOptTitle(0);
    ratio_firvar->SetTickx();
    ratio_firvar->Draw();
    ratio_firvar->cd();
    ratio_firvar->SetLogy(true);
    totalhistogram[0]->Draw("E");
    for (int i = 1; i < totalhistogram.size(); i++) {
        totalhistogram[i]->Draw("E same");
    }
    if (savecheck)
        //ratio_firvar->SaveAs(Form("./hRatio%s_(%2.f-%2.f).pdf",nameFitVar.Data(), multibin[0], multibin[1]));
    return fReturnHist;
}