//#include "DrawingHelper.C" // included in PlotXi1530.C
#include "PlotXi1530.C"
#include "SystematicHelper.cxx"

TString path = "./data/";

vector<vector<double>> multibincheck = {
    // {0, 100}}; // for test
    // {0, 100}, {0, 0.1}, {0, 10},  {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    {0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    //{0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 100}};

TH1* Xi1530PlotSystematics_multi(double multi_start = 0,
                                 double multi_end = 100,
                                 bool correl_error_skip = false,
                                 bool DrawQAplots = true);
TH1* ErrorSmoothing(TH1* hinput,
                    double smoothingCriteron = 1.5,
                    int binstart = 2,
                    int binend = 8);
TH1* CleanZeroBinhisto(TH1* hInput);
void Xi1530Systematics(bool QAPlots = false) {
    
    // output file
    TString bininfo = "";
    for (auto const& bin : multibincheck)
        bininfo += Form("%.0f", bin[0]);
    bininfo += Form("%.0f", multibincheck.back()[1]);
    TFile* fout_systematic =
        new TFile(Form("./AnalysisResults_Xi1530_systematic%s.root", bininfo.Data()), "RECREATE");
    
    //
    TString statfile;
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        statfile =
            Form("%sAnalysisResults_Extracted_1_Multi_%.2f-%.2f_Default1.root",
                 path.Data(), multibincheck[imultibin][0],
                 multibincheck[imultibin][1]);

        auto temp_multi = Xi1530PlotSystematics_multi(
            multibincheck[imultibin][0], multibincheck[imultibin][1], false,
            QAPlots);
        auto temp_multi_nocorerror = Xi1530PlotSystematics_multi(
            multibincheck[imultibin][0], multibincheck[imultibin][1], true,
            QAPlots);
        auto hBase = (TH1*)GetSpectrafromName(statfile);

        auto temp_multi_clean = (TH1*)CleanZeroBinhisto(temp_multi);
        auto temp_multi_nocorerror_clean =
            (TH1*)CleanZeroBinhisto(temp_multi_nocorerror);
        auto hBase_clean = (TH1*)CleanZeroBinhisto(hBase);
        fout_systematic->cd();
        hBase_clean->Write(Form("hSpectra_%.2f_%.2f_stat",
                                multibincheck[imultibin][0],
                                multibincheck[imultibin][1]));
        temp_multi_clean->Write(Form("hSpectra_%.2f_%.2f_sys",
                                     multibincheck[imultibin][0],
                                     multibincheck[imultibin][1]));
        temp_multi_nocorerror_clean->Write(
            Form("hSpectra_%.2f_%.2f_sys_noCorrelation",
                 multibincheck[imultibin][0], multibincheck[imultibin][1]));
    }

    fout_systematic->Close();
}

TH1* Xi1530PlotSystematics_multi(double multi_start,
                                 double multi_end,
                                 bool correl_error_skip,
                                 bool DrawQAplots) {
    // cout << "Input multi start: " << multi_start << ", multi end: " <<
    // multi_end << endl;
    double fYmax = 0.5;
    bool skipCorrelatedError = correl_error_skip;
    bool smoothing = true;

    double multi_start_topol = 0.; // for topol, use 0-100 case
    double multi_end_topol = 100.; //

    vector<TH1*> totalsystematic;
    vector<TString> totalsystematicname;

    // Deafult Spectra
    // Given multiplicity
    TString defaultfile = Form("%sAnalysisResults_Extracted_1_Multi_%.2f-%.2f_Default1.root",
             path.Data(), multi_start, multi_end);
    TH1* hBase = GetSpectrafromName(defaultfile);
    // 0-100 for topological/PID variation
    TFile* inputfile0100 = new TFile(
        Form("%sAnalysisResults_Extracted_1_Multi_%.2f-%.2f_Default1.root",
             path.Data(), multi_start_topol, multi_end_topol));
    TH1* hBase0100 = (TH1*)inputfile0100->Get("hXiSpectrum");  //

    // pT bin
    const TArrayD* ptbinarray = hBase->GetXaxis()->GetXbins();
    vector<double> ptbin;
    for (int i = 0; i < ptbinarray->GetSize(); i++)
        ptbin.push_back(ptbinarray->GetAt(i));
    TString fvarfile;

    // Signal Extraction
    // -----------------------------------------------------------
    vector<TH1*> SigExtErrorHistos;
    // Type1 variation: No detailed study needed. ----------------
    vector<TH1*> hSigExtSys_variations_type1;  // spectra
    vector<TString> SigExtSys_type1bins = {"LikeSignBkg", "BkgFit"};

    for (int i = 0; i < SigExtSys_type1bins.size(); i++) {
        fvarfile = Form(
            "%sAnalysisResults_Extracted_1_Multi_%.2f-%.2f_%s1.root",
            path.Data(), multi_start, multi_end, SigExtSys_type1bins[i].Data());
        auto temp = GetSpectrafromName(fvarfile);
        hSigExtSys_variations_type1.push_back(temp);
    }
    auto Systematic_SigExt = SystematicHelper(hBase);
    Systematic_SigExt.AddVariationtHistograms(hSigExtSys_variations_type1);
    Systematic_SigExt.SetBarlowCheck();
    // Systematic_SigExt.SetVerbose(); // for debuging
    Systematic_SigExt.SetDefaultYRange({0.5, 1.5});
    Systematic_SigExt.SetYaxisTitle("Ratio");
    Systematic_SigExt.InitAbsDiffRatioColors(false);

    auto hSigExtErrors = Systematic_SigExt.GetAbsDiffRatio();
    // auto hSigExtErrorMax = Systematic_SigExt.GetMaxAbsDiffRatio();
    // auto hSigExtErrorSum = Systematic_SigExt.GetSumAbsDiffRatio();
    for (auto const& hInput : hSigExtErrors) {
        SigExtErrorHistos.push_back(hInput);
    }

    // Type2 variation: detailed fit variation study. -------------
    //
    vector<TString> SigExtSys_type2bins = {"FitVar", "NormVar"};
    // "BinCount" variation is now excluded due to it doesn't affect so much.
    // and can't pass the barlow check.
    vector<TString> SigExtSys_type2subbins = {"Lm",    "Lp",    "Rm", "Rp",
                                             "Bothm", "Bothp", "p",  "m"};
    vector<vector<TH1*>> hSigExtSys_variations_type2(
        SigExtSys_type2bins.size(), {});

    int maximumVariation = 4; // from 5, variation can't describe spectra well.
    for (int i = 0; i < SigExtSys_type2bins.size(); i++) {
        for (int subbins = 0; subbins < SigExtSys_type2subbins.size();
             subbins++) {
            for (int var = 1; var < maximumVariation; var++) {
                fvarfile = Form(
                    "%sAnalysisResults_Extracted_1_Multi_%.2f-%.2f_%s%s%d.root",
                    path.Data(), multi_start, multi_end,
                    SigExtSys_type2bins[i].Data(),
                    SigExtSys_type2subbins[subbins].Data(), var);
                // Example)
                // AnalysisResults_Extracted_1_Multi_0.00-100.00_BinCountBothp4.root
                auto temp = GetSpectrafromName(fvarfile);
                if(!temp)
                    continue; // Anyway, we'll use available error only.
                hSigExtSys_variations_type2[i].push_back(temp);
            }
        }
    }

    vector<SystematicHelper> Systematic_SigExt_type2;
    vector<TH1*> hSigExtError_type2_Stdev;  // use stdev error
    for (int i = 0; i < SigExtSys_type2bins.size(); i++) {
        // The loop can make more flexible code structure.

        Systematic_SigExt_type2.push_back(
            SystematicHelper(hBase,hSigExtSys_variations_type2[i]));
        Systematic_SigExt_type2[i].SetBarlowCheck();
        // Systematic_SigExt_type2[i].SetVerbose(); // for debuging
        Systematic_SigExt_type2[i].SetDefaultYRange({0.5, 1.5});
        Systematic_SigExt_type2[i].SetYaxisTitle("Ratio");
        Systematic_SigExt_type2[i].InitAbsDiffRatioColors(false);
        hSigExtError_type2_Stdev.push_back(
            Systematic_SigExt_type2[i].GetStdevRatio());
    }
    for (auto const& hInput : hSigExtError_type2_Stdev)
        SigExtErrorHistos.push_back(hInput);
    
    // Error sum
    vector<double> SigExtErrorSum;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++) {
        SigExtErrorSum.push_back(0);
        for (auto const& hInput : SigExtErrorHistos) {
            if (hInput->GetBinContent(bin+1) > SigExtErrorSum[bin])
                SigExtErrorSum[bin] = hInput->GetBinContent(bin+1);
        }
    }
    auto hSigExtTotalError =
        (TH1*)MakeHistfromArray("SigExtErrorSum", SigExtErrorSum, ptbin);
    // -----------------------------------------------------------

    // Topological Cut variations
    // -----------------------------------------------------------
    vector<TH1*> hTopolsys_variation;
    vector<TString> TopologicalCutSystematic_bins = {
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
        // Reason: one source, one error.
        "Xi1530PionZVertexLoose",
        "DCADistLambdaDaughtersLoose",
        "DCADistXiDaughtersLoose",
        "DCADistLambdaPVLoose",
        "V0CosineOfPointingAngleTight",
        "CascadeCosineOfPointingAngleTight",
        "XiMassWindowLoose"};
    vector<TString> TopologicalCutSystematic_used;
    std::set_difference(TopologicalCutSystematic_bins.begin(),
                        TopologicalCutSystematic_bins.end(),
                        TopologicalCutSystematic_notuse.begin(),
                        TopologicalCutSystematic_notuse.end(),
                        std::inserter(TopologicalCutSystematic_used,
                                      TopologicalCutSystematic_used.begin()));
    // TString fvarfile; 
    for (int topological_cutvar = 0;
         topological_cutvar < TopologicalCutSystematic_bins.size();
         topological_cutvar++) {
        if (std::find(TopologicalCutSystematic_notuse.begin(),
                      TopologicalCutSystematic_notuse.end(),
                      TopologicalCutSystematic_bins[topological_cutvar]) !=
            TopologicalCutSystematic_notuse.end())
            continue;  // Not using sources.
        fvarfile =
            Form("%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
                 path.Data(), topological_cutvar + 6, multi_start_topol,
                 multi_end_topol,
                 TopologicalCutSystematic_bins[topological_cutvar].Data());
        auto temp = GetSpectrafromName(fvarfile);
        hTopolsys_variation.push_back(temp);
    }
    auto Systematic_TopolCuts = SystematicHelper(hBase0100);
    Systematic_TopolCuts.AddVariationtHistograms(hTopolsys_variation);
    Systematic_TopolCuts.SetBarlowCheck();
    // Systematic_TopolCuts.SetVerbose(); // for debuging
    Systematic_TopolCuts.SetDefaultYRange({0.8, 1.2});
    Systematic_TopolCuts.SetYaxisTitle("Ratio");
    Systematic_TopolCuts.InitAbsDiffRatioColors(false);

    auto hTopolCutErrors = Systematic_TopolCuts.GetAbsDiffRatio();
    // auto hTopolCutErrorMax = Systematic_TopolCuts.GetMaxAbsDiffRatio();
    auto hTopolCutErrorSum =
        Systematic_TopolCuts.GetSumAbsDiffRatio();  // use sum
    // -----------------------------------------------------------

    // PID Cut variations
    // -----------------------------------------------------------
    vector<TH1*> hPIDCutsys_variation;
    vector<TString> PIDCutSystematic_bins = {
        "TPCNsigmaXi1530PionLoose", "TPCNsigmaXi1530PionTight",
        "TPCNsigmaXiLoose", "TPCNsigmaXiTight"};
    vector<TString> PIDCutSystematic_notuse = {"TPCNsigmaXiLoose",
                                               "TPCNsigmaXi1530PionLoose"};
    vector<TString> PIDCutSystematic_used;
    std::set_difference(
        PIDCutSystematic_bins.begin(), PIDCutSystematic_bins.end(),
        PIDCutSystematic_notuse.begin(), PIDCutSystematic_notuse.end(),
        std::inserter(PIDCutSystematic_used, PIDCutSystematic_used.begin()));
    // TString fvarfile;
    for (int topological_cutvar = 0;
         topological_cutvar < PIDCutSystematic_bins.size();
         topological_cutvar++) {
        if (std::find(PIDCutSystematic_notuse.begin(),
                      PIDCutSystematic_notuse.end(),
                      PIDCutSystematic_bins[topological_cutvar]) !=
            PIDCutSystematic_notuse.end())
            continue;  // Not using sources.
        fvarfile = Form(
            "%sAnalysisResults_Extracted_%i_Multi_%.2f-%.2f_%s1.root",
            path.Data(), topological_cutvar + 2, multi_start_topol,
            multi_end_topol, PIDCutSystematic_bins[topological_cutvar].Data());
        auto temp = GetSpectrafromName(fvarfile);
        hPIDCutsys_variation.push_back(temp);
    }
    auto Systematic_PIDCuts = SystematicHelper(hBase0100);
    Systematic_PIDCuts.AddVariationtHistograms(hPIDCutsys_variation);
    Systematic_PIDCuts.SetBarlowCheck();
    // Systematic_PIDCuts.SetVerbose(); // for debuging
    Systematic_PIDCuts.SetDefaultYRange({0.8, 1.2});
    Systematic_PIDCuts.SetYaxisTitle("Ratio");
    Systematic_PIDCuts.InitAbsDiffRatioColors(false);

    auto hPIDCutErrors = Systematic_PIDCuts.GetAbsDiffRatio();
    auto hPIDCutErrorMax = Systematic_PIDCuts.GetMaxAbsDiffRatio();
    // -----------------------------------------------------------

    // Other Systematics
    // -----------------------------------------------------------
    //vector<double> ptbin = {0,   0.8, 1.2, 1.6, 2.0, 2.4,3.2, 4.0, 4.8, 5.6, 8.8, 15};
    // 1. Material Budget (Xi+-)
    // Study of strangeness production as a function of the charged particle
    // multiplicity in pp collisions at \sqrt{s} = 13 TeV
    // https://alice-notes.web.cern.ch/node/478
    // 4%, pT independent, due to the lack of knowledge of the material budget.
    vector<double> materialbudget;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
        materialbudget.push_back(0.04);
    auto hSysMaterialBudget =
        (TH1*)MakeHistfromArray("hSysMaterialBudget", materialbudget, ptbin);
    // 2. Mult. indipendent efficiencies
    // 2%, pT independent, due to the computation of the efficiencies in the
    // integrated case over multi-plicity
    vector<double> multiindepeffi;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
        multiindepeffi.push_back(0.02);
    auto hSysMultiIndepEffi =
        (TH1*)MakeHistfromArray("hSysMultiIndepEffi", multiindepeffi, ptbin);
    // 3. Tracking efficiencies
    // https://twiki.cern.ch/twiki/bin/view/ALICE/TrackingEfficiencyCharged
    // 3%, pT independent
    vector<double> trackingeffi;
    for (int bin = 0; bin < hBase->GetNbinsX(); bin++)
        trackingeffi.push_back(0.03);
    auto hSysTrackingEffi =
        (TH1*)MakeHistfromArray("hSysTrackingEffi", trackingeffi, ptbin);
    // -----------------------------------------------------------

    // Total systematic uncertainty
    // -----------------------------------------------------------
    totalsystematic.push_back(hSigExtTotalError);
    totalsystematicname.push_back("Signal Extraction");
    totalsystematic.push_back(hTopolCutErrorSum);
    totalsystematicname.push_back("Topological Cut variation");
    totalsystematic.push_back(hPIDCutErrorMax);
    totalsystematicname.push_back("TPC PID Cut variation");
    if (!skipCorrelatedError) {
        totalsystematic.push_back(hSysMaterialBudget);
        totalsystematicname.push_back("Material budget");
        totalsystematic.push_back(hSysMultiIndepEffi);
        totalsystematicname.push_back("Multiplicity indep. efficiency");
        totalsystematic.push_back(hSysTrackingEffi);
        totalsystematicname.push_back("Tracking Efficiency");
    }

    vector<TH1*> totalsystematic_beforesmoothing;
    if (smoothing) {
        for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++) {
            // keep old syserror
            totalsystematic_beforesmoothing.push_back(totalsystematic[sysbin]);
            // smoothed error
            totalsystematic[sysbin] = ErrorSmoothing(totalsystematic[sysbin],1.3);
        }
    }

    vector<double> totalSysError;  // Total sys error sum array
    for (int bin = 0; bin < totalsystematic[0]->GetNbinsX(); bin++) {
        double TotError = 0.;
        for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++) {
            TotError += pow(totalsystematic[sysbin]->GetBinContent(bin + 1), 2);
        }
        totalSysError.push_back(sqrt(TotError));
    }

    auto hSysSpectra = (TH1*)hBase->Clone();  // Spectrum with systematic error
    for (int bin = 0; bin < hSysSpectra->GetNbinsX(); bin++)
        hSysSpectra->SetBinError(
            bin + 1, totalSysError[bin] * hSysSpectra->GetBinContent(bin + 1));

    // return hSysSpectra; // It's okay to finish here.
    // -----------------------------------------------------------

    // QA Plots
    // -----------------------------------------------------------
    if (DrawQAplots) {
        auto hSysTotalSum =  // systematic error fraction
            (TH1*)MakeHistfromArray(
                Form("hSysTotalSum_%.2f-%.2f", multi_start, multi_end),
                totalSysError, ptbin);
        hSysTotalSum->SetLineColor(kBlack);
        hSysTotalSum->SetLineStyle(2);
        hSysTotalSum->SetLineWidth(2);

        vector<double> baseStatError = GetBinErrorFrations(hBase);
        auto hStatError =  // statistical error fraction
            (TH1*)MakeHistfromArray(
                Form("hSysTotalSum_%.2f-%.2f", multi_start, multi_end),
                baseStatError, ptbin);
        hStatError->SetFillColorAlpha(kBlack, 0.2);
        hStatError->SetLineColorAlpha(kBlack, 0.2);
        hStatError->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hStatError->SetMaximum(0.4);
        hStatError->SetMinimum(0);
        hStatError->GetXaxis()->SetRangeUser(0.8, 8.8);

        // Initialize the total systematic errors
        // color, axis range
        int sysColorPallet = GetSerialColors(totalsystematic.size());
        for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++) {
            totalsystematic[sysbin]->SetMarkerColor(
                sysColorPallet + totalsystematic.size() - sysbin-1);
            totalsystematic[sysbin]->SetLineColor(
                sysColorPallet + totalsystematic.size() - sysbin -1);
            totalsystematic[sysbin]->GetXaxis()->SetRangeUser(0.8, 8.8);
            totalsystematic[sysbin]->SetTitle("");
            totalsystematic[sysbin]->SetLineWidth(2);
        }

        TCanvas* cQA = new TCanvas(
            Form("cFitResults_%.2f-%.2f", multi_start, multi_end),
            Form("cFitResults_%.2f-%.2f", multi_start, multi_end), 960, 720);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        cQA->Draw();
        cQA->cd();
        // cTest->SetLogy();
        auto legend_TotalSys = new TLegend(.6, .63, .9, .88);
        legend_TotalSys->SetBorderSize(0);
        legend_TotalSys->SetFillStyle(0);
        legend_TotalSys->AddEntry(hStatError, "Statistical Error", "F");
        hStatError->Draw("BAR");
        for (int i = 0; i < totalsystematic.size(); i++) {
            legend_TotalSys->AddEntry(totalsystematic[i],
                                      totalsystematicname[i], "L");
            totalsystematic[i]->Draw("hist same");
        }
        legend_TotalSys->AddEntry(hSysTotalSum, "Total systematic uncertainty", "L");
        hSysTotalSum->Draw("hist same");
        legend_TotalSys->Draw();
        SaveCanvas(cQA,
                   Form("cQA_TotalError_%.2f-%.2f", multi_start, multi_end));

        // Signal Extraction QA histo
        auto legend_SigExt = new TLegend(.6, .73, .9, .88);
        legend_SigExt->SetBorderSize(0);
        legend_SigExt->SetFillStyle(0);
        vector<TString> SigExtSys_totalbins;
        for (auto const& hInput : SigExtSys_type1bins)
            SigExtSys_totalbins.push_back(hInput);
        for (auto const& hInput : SigExtSys_type2bins)
            SigExtSys_totalbins.push_back(hInput);
        int sysColorPallet_SigExt = GetSerialColors(SigExtSys_totalbins.size());
        for (int sysbin = 0; sysbin < SigExtErrorHistos.size(); sysbin++) {
            legend_SigExt->AddEntry(SigExtErrorHistos[sysbin],
                                 SigExtSys_totalbins[sysbin], "L");
            SigExtErrorHistos[sysbin]->SetMarkerColor(
                sysColorPallet_SigExt + SigExtErrorHistos.size() - sysbin - 1);
            SigExtErrorHistos[sysbin]->SetLineColor(
                sysColorPallet_SigExt + SigExtErrorHistos.size() - sysbin - 1);
            SigExtErrorHistos[sysbin]->GetXaxis()->SetRangeUser(0.8, 8.8);
            SigExtErrorHistos[sysbin]->SetTitle("");
            SigExtErrorHistos[sysbin]->SetLineWidth(2);
        }
        SigExtErrorHistos[0]->Draw("hist");
        for (int bin = 1; bin < SigExtErrorHistos.size(); bin++) {
            SigExtErrorHistos[bin]->Draw("hist same");
        }
        hSigExtTotalError->SetLineColor(kBlack);
        hSigExtTotalError->SetLineStyle(2);
        hSigExtTotalError->SetLineWidth(2);
        hSigExtTotalError->Draw("hist same");
        legend_SigExt->AddEntry(hSigExtTotalError, "Sum of signal extraction uncertainty",
                                "L");
        legend_SigExt->Draw();
        SaveCanvas(cQA,
                   Form("cQA_SigExt_%.2f-%.2f", multi_start, multi_end),"figs/SignalExtraction/");
        //-------------------------

        // Topological Cut variation QA histo
        auto legend_Topolcut = new TLegend(.6, .73, .9, .88);
        legend_Topolcut->SetBorderSize(0);
        legend_Topolcut->SetFillStyle(0);
        int sysColorPallet_TopolCuts = GetSerialColors(hTopolCutErrors.size());
        for (int sysbin = 0; sysbin < hTopolCutErrors.size(); sysbin++) {
            legend_Topolcut->AddEntry(hTopolCutErrors[sysbin],
                                      TopologicalCutSystematic_used[sysbin],
                                      "L");
            hTopolCutErrors[sysbin]->SetMarkerColor(
                sysColorPallet_TopolCuts + hTopolCutErrors.size() - sysbin - 1);
            hTopolCutErrors[sysbin]->SetLineColor(
                sysColorPallet_TopolCuts + hTopolCutErrors.size() - sysbin - 1);
            hTopolCutErrors[sysbin]->GetXaxis()->SetRangeUser(0.8, 8.8);
            hTopolCutErrors[sysbin]->SetTitle("");
            hTopolCutErrors[sysbin]->SetLineWidth(2);
        }
        hTopolCutErrors[0]->Draw("hist");
        for (int bin = 1; bin < hTopolCutErrors.size(); bin++) {
            hTopolCutErrors[bin]->Draw("hist same");
        }
        hTopolCutErrorSum->SetLineColor(kBlack);
        hTopolCutErrorSum->SetLineStyle(2);
        hTopolCutErrorSum->SetLineWidth(2);
        hTopolCutErrorSum->Draw("hist same");
        legend_Topolcut->AddEntry(hTopolCutErrorSum,
                                  "Sum of topological cut uncertainty", "L");
        legend_Topolcut->Draw();
        SaveCanvas(cQA, Form("cQA_Topolcut_%.2f-%.2f", multi_start, multi_end),
                   "figs/TopologicalCutVar/");
        //-------------------------

        // TPC PID Cut variation QA histo
        auto legend_TPCPIDcut = new TLegend(.6, .73, .9, .88);
        legend_TPCPIDcut->SetBorderSize(0);
        legend_TPCPIDcut->SetFillStyle(0);
        int sysColorPallet_PIDCuts = GetSerialColors(hPIDCutErrors.size());
        for (int sysbin = 0; sysbin < hPIDCutErrors.size(); sysbin++) {
            legend_TPCPIDcut->AddEntry(hPIDCutErrors[sysbin],
                                       PIDCutSystematic_used[sysbin],
                                       "L");
            hPIDCutErrors[sysbin]->SetMarkerColor(
                sysColorPallet_PIDCuts + hPIDCutErrors.size() - sysbin - 1);
            hPIDCutErrors[sysbin]->SetLineColor(
                sysColorPallet_PIDCuts + hPIDCutErrors.size() - sysbin - 1);
            hPIDCutErrors[sysbin]->GetXaxis()->SetRangeUser(0.8, 8.8);
            hPIDCutErrors[sysbin]->SetTitle("");
            hPIDCutErrors[sysbin]->SetLineWidth(2);
        }
        hPIDCutErrors[0]->Draw("hist");
        for (int bin = 1; bin < hPIDCutErrors.size(); bin++) {
            hPIDCutErrors[bin]->Draw("hist same");
        }
        hPIDCutErrorMax->SetLineColor(kBlack);
        hPIDCutErrorMax->SetLineStyle(2);
        hPIDCutErrorMax->SetLineWidth(2);
        hPIDCutErrorMax->Draw("hist same");
        legend_TPCPIDcut->AddEntry(hPIDCutErrorMax,
                                   "Sum of TPC PID cut uncertainty", "L");
        legend_TPCPIDcut->Draw();
        SaveCanvas(cQA, Form("cQA_TPCPIDcut_%.2f-%.2f", multi_start, multi_end),
                   "figs/TPCPIDCutVar/");
        //-------------------------

        // Smoothing Check
        if (smoothing) {
            vector<double> totalSysError_noSmoothing;
            for (int bin = 0; bin < totalsystematic[0]->GetNbinsX(); bin++) {
                double TotError = 0.;
                for (int sysbin = 0; sysbin < totalsystematic.size();
                     sysbin++) {
                    TotError += pow(
                        totalsystematic_beforesmoothing[sysbin]->GetBinContent(
                            bin + 1),
                        2);
                }
                totalSysError_noSmoothing.push_back(sqrt(TotError));
            }

            // Initialize the total systematic errors
            // color, axis range
            for (int sysbin = 0; sysbin < totalsystematic.size(); sysbin++) {
                totalsystematic_beforesmoothing[sysbin]->SetMarkerColor(
                    sysColorPallet + totalsystematic.size() - sysbin - 1);
                totalsystematic_beforesmoothing[sysbin]->SetLineColor(
                    sysColorPallet + totalsystematic.size() - sysbin - 1);
                totalsystematic_beforesmoothing[sysbin]
                    ->GetXaxis()
                    ->SetRangeUser(0.8, 8.8);
                totalsystematic_beforesmoothing[sysbin]->SetTitle("");
                totalsystematic_beforesmoothing[sysbin]->SetLineWidth(2);
            }

            auto hSysTotalSum_noSmoothing =  // systematic error fraction
                (TH1*)MakeHistfromArray(
                    Form("hSysTotalSum_%.2f-%.2f", multi_start, multi_end),
                    totalSysError_noSmoothing, ptbin);
            hSysTotalSum_noSmoothing->SetLineColor(kRed);
            hSysTotalSum_noSmoothing->SetLineStyle(2);
            hSysTotalSum_noSmoothing->SetLineWidth(2);

            hStatError->Draw("BAR");
            for (int i = 0; i < totalsystematic_beforesmoothing.size(); i++) {
                totalsystematic_beforesmoothing[i]->Draw("hist same");
            }
            hSysTotalSum_noSmoothing->Draw("hist same");
            legend_TotalSys->Draw();
            SaveCanvas(cQA, Form("cQA_TotalError_%.2f-%.2f_noSmoothing",
                                 multi_start, multi_end),"figs/NoSmoothing/");
        }

        // QA Histgram of the each variations.
        Systematic_SigExt.SetVariationNames(SigExtSys_type1bins);
        Systematic_TopolCuts.SetVariationNames(TopologicalCutSystematic_used);
        Systematic_PIDCuts.SetVariationNames(PIDCutSystematic_used);
        Systematic_SigExt.InitVariationColors();
        Systematic_TopolCuts.InitVariationColors();
        Systematic_PIDCuts.InitVariationColors();
        SaveCanvas(
            Systematic_SigExt.GetQAPlot("AbsDiffRatio"),
            Form("cQA_SignalExtraction_%.2f-%.2f", multi_start, multi_end),
            "figs/detailedQA/");
        SaveCanvas(Systematic_TopolCuts.GetQAPlot("AbsDiffRatio"),
                   Form("cQA_ToplogicalCuts_%.2f-%.2f", multi_start, multi_end),
                   "figs/detailedQA/");
        SaveCanvas(Systematic_PIDCuts.GetQAPlot("AbsDiffRatio"),
                   Form("cQA_PIDCuts_%.2f-%.2f", multi_start,
                        multi_end), "figs/detailedQA/");
    }
    return hSysSpectra;
}
TH1* ErrorSmoothing(TH1* hinput,
                    double smoothingCriteron,  // 1.5
                    int binstart,              // 2
                    int binend) {              // 8
    auto htemp = (TH1*)hinput->Clone();
    for (int bin = 1; bin < htemp->GetNbinsX(); bin++) {
        if (bin < binstart)
            continue;  // like 0 - 0.8 GeV/c bin
        if (binend < bin)
            continue;  // like 8.8 - 15.0 GeV/c bin
        double Error1 = htemp->GetBinContent(bin + 1);   // current bin
        double Error0 = htemp->GetBinContent(bin);       // post bin
        double Error2 = htemp->GetBinContent(bin + 2);   // next bin
        double smoothedValue = htemp->GetBinContent(bin + 1);

        double averageOfNeighbor = (Error1 + Error0 + Error2) / 3;

        if (Error1 > smoothingCriteron * averageOfNeighbor)
            // only if the value is 'significantly' larger than neighbor.
            smoothedValue = averageOfNeighbor;
        htemp->SetBinContent(bin + 1, smoothedValue);
    }

    return htemp;
}
TH1* CleanZeroBinhisto(TH1* hInput){
    vector<double> temp_y;
    vector<double> temp_e;
    vector<double> ptbin_final = {0.8, 1.2, 1.6, 2.0, 2.4,
                                  3.2, 4.0, 4.8, 5.6, 8.8};
    for (int bin = 0; bin < hInput->GetNbinsX(); bin++) {
        if ((bin == 0) || (bin == 10) || (bin == 11))
            continue;
        temp_y.push_back(hInput->GetBinContent(bin + 1));
        temp_e.push_back(hInput->GetBinError(bin + 1));
    }
    auto tempReturn = MakeHistfromArray("", temp_y, ptbin_final, temp_e);
    return tempReturn;
}