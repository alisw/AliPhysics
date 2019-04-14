#include <AliFigure.h>
#include <AliPWGHistoTools.h>
#include "AdditionalFunctions.h"
#include "FitParticle.C"

// Configure
vector<vector<double>> multibin = {{0, 100}, {0, 10},  {10, 30},
                                   {30, 50}, {50, 70}, {70, 100}};
                                   //{30, 50}, {50, 100}};
//{0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}, {30, 100}};
//{0, 100}};
vector<double> spectrafitrange = {0.8, 8.8};
vector<double> integralrange = {0, 100};
TString finalfile = "AnalysisResults_Xi1530_systematic.root";
TString workdirectory = "/Users/blim/alidock/Postprocessing/data/";
bool staterror = true;

const char* outPut = "Xi1530YieldMean.root";

void Xi1530ExtractYieldMean();
void Xi1530ExtractYieldMean(TString inputfile);
vector<TH1D*> GetSysSpectra(vector<vector<double>> input_multibin);
vector<TH1D*> GetStatSpectra(vector<vector<double>> input_multibin);
TH1D* GetSpectrasys(double multi_start, double multi_end);
TH1D* GetSpectrastat(double multi_start, double multi_end);
TH1* GetSpectraTotalError(TH1* hstat, TH1* hsys);
TH1* ReturnExtremeHighHisto(TH1* hin);
TH1* ReturnExtremeLowHisto(TH1* hin);
TH1* ReturnExtremeHisto(TH1* hin, Float_t sign);
TH1* ReturnRandom(TH1* hin);
TH1* ReturnCoherentRandom(TH1* hin);
TH1D* MakeHistfromArray(char const* name,
                        vector<double> dArray,
                        vector<double> eArray,
                        vector<double> ptbin,
                        const char* foption = "");
vector<double> GetdNdetawithError(double multi_start, double multi_end);
vector<double> GetPidNdetawithError(double multi_start, double multi_end);
void Xi1530ExtractYieldMean() {
    TFile* foutput = new TFile(outPut);
    foutput->Close();
    delete foutput;

    vector<TH1D*> hspectra_sys = GetSysSpectra(multibin);
    vector<TH1D*> hspectra_stat = GetStatSpectra(multibin);
    vector<TH1D*> hspectra_tot;       // sys+stat
    vector<TH1D*> hspectra_hextreme;  // extreme hard case
    vector<TH1D*> hspectra_lextreme;  // extreme soft case
    vector<TH1D*> hspectra_bextreme;  // extreme high case
    vector<TH1D*> hspectra_sextreme;  // extreme low case

    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        // total error
        hspectra_tot.push_back((TH1D*)GetSpectraTotalError(
            hspectra_stat[imultibin], hspectra_sys[imultibin]));
        // extremely soft case
        hspectra_lextreme.push_back(
            (TH1D*)ReturnExtremeHisto(hspectra_sys[imultibin], -1));
        // extremely hard case
        hspectra_hextreme.push_back(
            (TH1D*)ReturnExtremeHisto(hspectra_sys[imultibin], +1));
        // extremely high case
        hspectra_bextreme.push_back(
            (TH1D*)ReturnExtremeHighHisto(hspectra_sys[imultibin]));
        // extremely low case
        hspectra_sextreme.push_back(
            (TH1D*)ReturnExtremeLowHisto(hspectra_sys[imultibin]));
    }

    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        cout << "Bin: " << multibin[imultibin][0] << " - "
             << multibin[imultibin][1] << endl;
        FitParticle(hspectra_sys[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);
        FitParticle(hspectra_stat[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);
        FitParticle(hspectra_tot[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);
        FitParticle(hspectra_lextreme[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);
        FitParticle(hspectra_hextreme[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);
        FitParticle(hspectra_bextreme[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);
        FitParticle(hspectra_sextreme[imultibin], "Xi*0", spectrafitrange[0],
                    spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                    outPut, 0, integralrange[0], integralrange[1]);

    }
}
void Xi1530ExtractYieldMean(TString inputfile) {
    TFile* finput = new TFile(inputfile.Data(), "open");

    // ch <dn/deta>
    vector<double> dNdetaAxis;
    vector<double> dNdetaAxis_e;
    vector<double> dNdetaAxis_e_half;
    vector<double> zeroerror;

    // pion <dn/dy>
    vector<double> dNdeta_pi;
    vector<double> dNdeta_pi_e;

    vector<double> yield;
    vector<double> yield_state;
    vector<double> yield_syshe;
    vector<double> yield_sysle;
    vector<double> meanpt;
    vector<double> meanpt_state;
    vector<double> meanpt_syshe;
    vector<double> meanpt_sysle;
    vector<double> extraYields;

    vector<TH1D*> hspectra_sys = GetSysSpectra(multibin);
    vector<TH1D*> hspectra_stat = GetStatSpectra(multibin);
    vector<TH1D*> hspectra_tot;       // sys+stat

    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        TH1* hcentral = (TH1*)finput->Get(
            Form("hFit_h_%.2f-%.2f_stat_tot_fLevi_fit_%.2f_%.2f",
                 multibin[imultibin][0], multibin[imultibin][1],
                 spectrafitrange[0], spectrafitrange[1]));
        TH1* hhard = (TH1*)finput->Get(
            Form("hFit_h_%.2f-%.2f_sys_extremehard_max_fLevi_fit_%.2f_%.2f",
                 multibin[imultibin][0], multibin[imultibin][1],
                 spectrafitrange[0], spectrafitrange[1]));
        TH1* hsoft = (TH1*)finput->Get(
            Form("hFit_h_%.2f-%.2f_sys_extremehard_min_fLevi_fit_%.2f_%.2f",
                 multibin[imultibin][0], multibin[imultibin][1],
                 spectrafitrange[0], spectrafitrange[1]));
        TH1* hhigh = (TH1*)finput->Get(
            Form("hFit_h_%.2f-%.2f_sys_extremehigh_fLevi_fit_%.2f_%.2f",
                 multibin[imultibin][0], multibin[imultibin][1],
                 spectrafitrange[0], spectrafitrange[1]));
        TH1* hlow = (TH1*)finput->Get(
            Form("hFit_h_%.2f-%.2f_sys_extremelow_fLevi_fit_%.2f_%.2f",
                 multibin[imultibin][0], multibin[imultibin][1],
                 spectrafitrange[0], spectrafitrange[1]));

        yield.push_back(hcentral->GetBinContent(1));
        yield_state.push_back(hcentral->GetBinError(1));
        yield_syshe.push_back(
            abs(hhigh->GetBinContent(1) - hcentral->GetBinContent(1)));
        yield_sysle.push_back(
            abs(hlow->GetBinContent(1) - hcentral->GetBinContent(1)));

        meanpt.push_back(hcentral->GetBinContent(2));
        meanpt_state.push_back(hcentral->GetBinError(2));
        meanpt_syshe.push_back(
            abs(hhard->GetBinContent(2) - hcentral->GetBinContent(2)));
        meanpt_sysle.push_back(
            abs(hsoft->GetBinContent(2) - hcentral->GetBinContent(2)));

        vector<double> temp =
            GetdNdetawithError(multibin[imultibin][0], multibin[imultibin][1]);
        dNdetaAxis.push_back(temp[0]);
        dNdetaAxis_e.push_back(temp[1]);
        dNdetaAxis_e_half.push_back(temp[1] / 2);

        vector<double> temp2 = GetPidNdetawithError(multibin[imultibin][0],
                                                    multibin[imultibin][1]);
        dNdeta_pi.push_back(temp2[0]);
        dNdeta_pi_e.push_back(temp2[1]);

        zeroerror.push_back(0);

        // Get the statistical error from spectra.
        // from AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C
        if(staterror){
        hspectra_tot.push_back((TH1D*)GetSpectraTotalError(
        hspectra_stat[imultibin], hspectra_sys[imultibin]));

        double integral = yield[imultibin];
        double mean = meanpt[imultibin];
        vector<double> temp_vector;
        // random generation with integration (coarse) 
        TH1* hIntegral_tmp =
            new TH1F("hIntegral_tmp", "", 1000, 0.75 * integral, 1.25 * integral);
        TH1* hMean_tmp = new TH1F("hMean_tmp", "", 1000, 0.75 * mean, 1.25 * mean);
        for (Int_t irnd = 0; irnd < 100; irnd++) {
            TH1* hrnd = ReturnRandom(hspectra_stat[imultibin]);
            temp_vector = FitParticleReturn(hrnd, "Xi*0", spectrafitrange[0],
                spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                outPut, 0, integralrange[0], integralrange[1]);
            hIntegral_tmp->Fill(temp_vector[0]); // yield
            hMean_tmp->Fill(temp_vector[2]);     // mean pT
            delete hrnd;
        }
        // random generation with integration (fine) 
        TH1* hIntegral = new TH1F("hIntegral", "", 100,
                 hIntegral_tmp->GetMean() - 10. * hIntegral_tmp->GetRMS(),
                 hIntegral_tmp->GetMean() + 10. * hIntegral_tmp->GetRMS());
        TH1* hMean = new TH1F("hMean", "", 100,
                          hMean_tmp->GetMean() - 10. * hMean_tmp->GetRMS(),
                          hMean_tmp->GetMean() + 10. * hMean_tmp->GetRMS());
        for (Int_t irnd = 0; irnd < 1000; irnd++) {
            TH1* hrnd = ReturnRandom(hspectra_stat[imultibin]);
            temp_vector = FitParticleReturn(hrnd, "Xi*0", spectrafitrange[0],
                spectrafitrange[1], -1, kFitLevi, AliPWGFunc::kdNdpt,
                outPut, 0, integralrange[0], integralrange[1]);
            hIntegral->Fill(temp_vector[0]); // yield
            hMean->Fill(temp_vector[2]);     // mean pT
            delete hrnd;
        }
        TF1* gaus = (TF1*)gROOT->GetFunction("gaus");

        hIntegral->Fit(gaus, "q");
        hMean->Fit(gaus, "q");

        yield_state.push_back(integral * gaus->GetParameter(2) /
                   gaus->GetParameter(1));
        meanpt_state.push_back(mean * gaus->GetParameter(2) /
               gaus->GetParameter(1));
        }
    }
    

    TGraphErrors* ge_stat =
        new TGraphErrors(multibin.size(), &dNdetaAxis[0], &yield[0],
                         &zeroerror[0], &yield_state[0]);
    TGraphAsymmErrors* ge_sys = new TGraphAsymmErrors(
        multibin.size(), &dNdetaAxis[0], &yield[0], &dNdetaAxis_e[0],
        &dNdetaAxis_e[0], &yield_sysle[0], &yield_syshe[0]);
    ge_stat->SetTitle("");
    ge_stat->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    ge_stat->GetYaxis()->SetTitle("d#it{N}/dy");
    ge_stat->GetXaxis()->SetRangeUser(0, 30);
    ge_stat->GetYaxis()->SetRangeUser(0, 0.007);
    ge_stat->SetLineColor(2);
    ge_stat->SetMarkerColor(2);
    ge_stat->SetMarkerStyle(20);
    ge_stat->SetMarkerSize(0.5);
    ge_stat->SetFillColor(2);

    TGaxis::SetMaxDigits(3);
    ge_sys->SetTitle("");
    ge_sys->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    ge_sys->GetYaxis()->SetTitle("d#it{N}/dy");
    ge_sys->GetXaxis()->SetRangeUser(0, 30);
    ge_sys->GetYaxis()->SetRangeUser(0, 0.007);
    // ge_sys->SetFillStyle(3001);
    ge_sys->SetLineColor(2);
    ge_sys->SetMarkerColor(2);
    ge_sys->SetFillColor(0);
    ge_sys->GetXaxis()->SetLabelSize(0.05);
    ge_sys->GetXaxis()->SetTitleSize(0.05);
    ge_sys->GetYaxis()->SetLabelSize(0.05);
    ge_sys->GetYaxis()->SetTitleSize(0.05);
    ge_sys->GetYaxis()->SetTitleOffset(1.0);

    vector<double> x7 = {6.01};
    vector<double> x7e = {0.10};

    // HEP data
    vector<double> y7 = {2.56e-3};
    vector<double> y7e = {0.07e-3};
    vector<double> y7l = {0.37e-3};
    vector<double> y7h = {0.40e-3};

    vector<double> pt7 = {1.31};
    vector<double> pt7e = {0.02};
    vector<double> pt7l = {0.09};
    vector<double> pt7h = {0.09};
    /*
    // Refit
    vector<double> y7 = {0.00245381};
    vector<double> y7e = {2.53791e-05};
    vector<double> y7l = {0.000199095};
    vector<double> y7h = {0.000845923};
    */
    /*
    ge_stat->SetPoint(multibin.size()+1, 6.01, 2.56e-3);
    ge_stat->SetPointError(multibin.size()+1, 0.1, 0.07e-3);

    ge_sys->SetPoint(multibin.size()+1, 6.01, 2.56e-3);
    ge_sys->SetPointError(multibin.size()+1,0.1,0.1,0.37e-3,0.40e-3);
    */

    TGraphErrors* ge_stat_7TeV =
        new TGraphErrors(1, &x7[0], &y7[0], &x7e[0], &y7e[0]);
    ge_stat_7TeV->SetLineColor(1);
    ge_stat_7TeV->SetMarkerColor(1);
    ge_stat_7TeV->SetMarkerStyle(20);
    ge_stat_7TeV->SetMarkerSize(0.5);
    ge_stat_7TeV->SetFillColor(1);
    TGraphAsymmErrors* ge_sys_7TeV = new TGraphAsymmErrors(
        1, &x7[0], &y7[0], &x7e[0], &x7e[0], &y7l[0], &y7h[0]);
    ge_sys_7TeV->SetLineColor(1);
    ge_sys_7TeV->SetMarkerColor(1);
    ge_sys_7TeV->SetFillColor(0);

    // 7TeV Results
    // Old data
    // from HEPDATA https://www.hepdata.net/record/ins1300380

    vector<double> ptbin2 = {0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6};
    vector<double> CorrectedYeild_7TeV2 = {0.00138,  0.00109,  0.00078,
                                           0.00046,  0.000226, 6.6e-05,
                                           2.34e-05, 8.1e-06};
    vector<double> CorrectedYeild_syserr_7TeV = {0.00017,  6.0e-05,  5.0e-05,
                                                 2.2e-05,  1.15e-05, 3.73e-06,
                                                 1.96e-06, 4.54e-07};
    vector<double> CorrectedYeild_staterr_7TeV = {4.0e-05, 3.0e-05,  2.0e-05,
                                                  1.5e-05, 3.38e-06, 2.08e-6,
                                                  8.11e-7, 9.09e-7};

    auto hSpectra_7TeV_syserr = MakeHistfromArray(
        "7TeV Spectra with systematic error", CorrectedYeild_7TeV2,
        CorrectedYeild_syserr_7TeV, ptbin2);
    auto hSpectra_7TeV_staterr = MakeHistfromArray(
        "7TeV Spectra with statistical error", CorrectedYeild_7TeV2,
        CorrectedYeild_staterr_7TeV, ptbin2);

    /*
    AliFigure *fig = new AliFigure("fig", "AliFigure", 800, 600);
    fig->SetStatus(AliFigure::kWorkInProgress);
    fig->SetLogoPos(AliFigure::kN);
    fig->SetCollSystem("pp #sqrt{s} = 13 TeV");
    //fig->SetDataSample("LHC15fi_16deghijklop_17cefgijklmor");
    fig->SetDataSample("");

    fig->cd();
    ge_sys->Draw("a5");
    ge_stat->Draw("P");
    ge_sys_7TeV->Draw("5");
    ge_stat_7TeV->Draw("P");


    gPad->Update();
    gPad->SaveAs("plot.pdf");
    */
    TCanvas* cCanvas = new TCanvas("cCanvas", "cCanvas", 1200, 720);
    // gStyle->SetOptTitle(0);
    cCanvas->SetTickx();
    cCanvas->Draw();
    cCanvas->cd();
    // cCanvas->SetLogy(true);
    ge_sys->Draw("a5");
    ge_stat->Draw("P");
    ge_sys_7TeV->Draw("5");
    ge_stat_7TeV->Draw("P");

    TText* fStatusPad = new TText(0.5, 0.15, "work in progress");
    fStatusPad->SetNDC();
    fStatusPad->SetTextSize(25);
    fStatusPad->SetTextFont(43);
    fStatusPad->SetTextAlign(22);
    fStatusPad->SetTextColor(kBlack);
    fStatusPad->Draw();

    auto legendyield = new TLegend(.6, .2, .8, .30);
    legendyield->SetBorderSize(0);
    legendyield->SetFillStyle(0);
    legendyield->AddEntry(ge_sys, "pp 13 TeV", "F");
    legendyield->AddEntry(ge_sys_7TeV, "pp 7 TeV", "F");
    legendyield->Draw();
    cCanvas->SaveAs(
        Form("%s1_Multi_0.00-100.00_Default1/yield.pdf", workdirectory.Data()));

    // Mean pt
    TGraphErrors* gpt_stat =
        new TGraphErrors(multibin.size(), &dNdetaAxis[0], &meanpt[0],
                         &zeroerror[0], &meanpt_state[0]);
    TGraphAsymmErrors* gpt_sys = new TGraphAsymmErrors(
        multibin.size(), &dNdetaAxis[0], &meanpt[0], &dNdetaAxis_e[0],
        &dNdetaAxis_e[0], &meanpt_sysle[0], &meanpt_syshe[0]);
    gpt_stat->SetTitle("");
    gpt_stat->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    gpt_stat->GetYaxis()->SetTitle("< p_{T} > (GeV/c)");
    gpt_stat->GetXaxis()->SetRangeUser(0, 30);
    gpt_stat->GetYaxis()->SetRangeUser(1, 1.7);
    gpt_stat->SetLineColor(2);
    gpt_stat->SetMarkerColor(2);
    gpt_stat->SetMarkerStyle(20);
    gpt_stat->SetMarkerSize(0.5);
    gpt_stat->SetFillColor(2);

    gpt_sys->SetTitle("");
    gpt_sys->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    gpt_sys->GetYaxis()->SetTitle("< p_{T} > (GeV/c)");
    gpt_sys->GetXaxis()->SetRangeUser(0, 30);
    gpt_sys->GetYaxis()->SetRangeUser(1, 1.7);
    // ge_sys->SetFillStyle(3001);
    gpt_sys->SetLineColor(2);
    gpt_sys->SetMarkerColor(2);
    gpt_sys->SetFillColor(0);
    gpt_sys->GetXaxis()->SetLabelSize(0.05);
    gpt_sys->GetXaxis()->SetTitleSize(0.05);
    gpt_sys->GetYaxis()->SetLabelSize(0.05);
    gpt_sys->GetYaxis()->SetTitleSize(0.05);
    gpt_sys->GetYaxis()->SetTitleOffset(1.0);

    TGraphErrors* gpt_stat_7TeV =
        new TGraphErrors(1, &x7[0], &pt7[0], &x7e[0], &pt7e[0]);
    gpt_stat_7TeV->SetLineColor(1);
    gpt_stat_7TeV->SetMarkerColor(1);
    gpt_stat_7TeV->SetMarkerStyle(20);
    gpt_stat_7TeV->SetMarkerSize(0.5);
    gpt_stat_7TeV->SetFillColor(1);
    TGraphAsymmErrors* gpt_sys_7TeV = new TGraphAsymmErrors(
        1, &x7[0], &pt7[0], &x7e[0], &x7e[0], &pt7l[0], &pt7h[0]);
    gpt_sys_7TeV->SetLineColor(1);
    gpt_sys_7TeV->SetMarkerColor(1);
    gpt_sys_7TeV->SetFillColor(0);

    TCanvas* cCanvas2 = new TCanvas("cCanvas2", "cCanvas2", 1200, 720);
    // gStyle->SetOptTitle(0);
    cCanvas2->SetTickx();
    cCanvas2->Draw();
    cCanvas2->cd();
    // cCanvas->SetLogy(true);
    gpt_sys->Draw("a5");
    gpt_stat->Draw("P");
    gpt_sys_7TeV->Draw("5");
    gpt_stat_7TeV->Draw("P");

    fStatusPad->Draw();

    legendyield->Draw();
    cCanvas2->SaveAs(Form("%s1_Multi_0.00-100.00_Default1/meanpt.pdf",
                          workdirectory.Data()));

    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        cout << "Multi bin: " << multibin[imultibin][0] << " - "
             << multibin[imultibin][1] << ", Yield: " << yield[imultibin]
             << " +- " << yield_state[imultibin] << " + "
             << yield_syshe[imultibin] << " - " << yield_sysle[imultibin]
             << ", mean pT: " << meanpt[imultibin] << " +- "
             << meanpt_state[imultibin] << " + " << meanpt_syshe[imultibin]
             << " - " << meanpt_sysle[imultibin] << endl;
    }
}
vector<TH1D*> GetSysSpectra(vector<vector<double>> input_multibin) {
    vector<TH1D*> buffer;
    for (int imultibin = 0; imultibin < input_multibin.size(); imultibin++) {
        auto temp_spectra_multi = GetSpectrasys(input_multibin[imultibin][0],
                                                input_multibin[imultibin][1]);
        buffer.push_back(temp_spectra_multi);
    }
    return buffer;
}
vector<TH1D*> GetStatSpectra(vector<vector<double>> input_multibin) {
    vector<TH1D*> buffer;
    for (int imultibin = 0; imultibin < input_multibin.size(); imultibin++) {
        auto temp_spectra_multi = GetSpectrastat(input_multibin[imultibin][0],
                                                 input_multibin[imultibin][1]);
        buffer.push_back(temp_spectra_multi);
    }
    return buffer;
}
TH1D* GetSpectrasys(double multi_start, double multi_end) {
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_sys", multi_start, multi_end));

    return hr;
}
TH1D* GetSpectrastat(double multi_start, double multi_end) {
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_stat", multi_start, multi_end));

    return hr;
}
TH1* GetSpectraTotalError(TH1* hstat, TH1* hsys) {
    TH1* htot = (TH1*)hstat->Clone(Form("%s_tot", hstat->GetName()));
    for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
        htot->SetBinError(ibin + 1,
                          TMath::Sqrt(pow(hsys->GetBinError(ibin + 1), 2) +
                                      pow(hstat->GetBinError(ibin + 1), 2)));
    }
    return htot;
}
vector<double> GetdNdetawithError(double multi_start, double multi_end) {
    // Return dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}

    vector<double> returnarray;

    //--dNdeta histogram
    // Ref: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/
    //      ReferenceMult#Multiplicity_dependent_pp_at_AN2
    // LHC16k data.
    // Error was asym error, so choosed bigger error for sym error.

    vector<double> dNchdeta_multibin = {0,  1,  5,  10, 15, 20,
                                        30, 40, 50, 70, 100};
    vector<double> dNchdeta = {0,     26.02, 20.02, 16.17, 13.77, 12.04,
                               10.02, 7.95,  6.32,  4.50,  2.55};
    vector<double> dNchdeta_e = {0,    0.35, 0.27, 0.22, 0.19, 0.17,
                                 0.14, 0.11, 0.09, 0.07, 0.04};

    // input must be in the multiplicity range
    if (std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(),
                  multi_start) == end(dNchdeta_multibin))
        return {99, 99};
    if (std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(),
                  multi_end) == end(dNchdeta_multibin))
        return {99, 99};

    // special cases
    if ((multi_start == 0) && (multi_end == 0.01)) {
        returnarray = {35.37, 0.92};
    } else if ((multi_start == 0.01) && (multi_end == 0.1)) {
        returnarray = {30.89, 0.57};
    } else if ((multi_start == 0) && (multi_end == 5)) {
        returnarray = {21.20, 0.28};
    } else if ((multi_start == 0) && (multi_end == 100)) {
        returnarray = {6.94, 0.10};
    } else if ((multi_start == 0.1) && (multi_end == 0.5)) {
        returnarray = {26.96, 0.37};
    } else if ((multi_start == 0.5) && (multi_end == 1)) {
        returnarray = {24.23, 0.36};
    }
    // Common case
    else {
        // Value
        vector<double>::iterator itr_left = find(
            dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_start);
        vector<double>::iterator itr_right =
            find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_end);
        int left = distance(dNchdeta_multibin.begin(), itr_left);
        int right = distance(dNchdeta_multibin.begin(), itr_right);

        int gap = right - left;

        double result = 0.;
        for (int i = 1; i < gap + 1; i++)
            result += dNchdeta[i + left] * (dNchdeta_multibin[i + left] -
                                            dNchdeta_multibin[i + left - 1]);

        result /= (multi_end - multi_start);
        returnarray.push_back(result);

        // Error
        double error = 0.;
        for (int i = 1; i < gap + 1; i++)
            error += pow(dNchdeta_e[i + left], 2);

        error = sqrt(error);
        returnarray.push_back(error);
    }

    return returnarray;
}
vector<double> GetPidNdetawithError(double multi_start, double multi_end) {
    // Return pion's dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}

    vector<double> returnarray;

    //--dNdeta histogram of pion+ + pion-
    // Ref: given by Anders, need to update ref.
    // Error was stat+sys error, so choosed bigger error(sym) error.

    vector<double> dNchdeta_multibin = {0,  1,  5,  10, 15, 20,
                                        30, 40, 50, 70, 100};
    vector<double> dNchdeta = {0,
                               24.605455025,
                               19.016207266,
                               15.477779961,
                               13.241832514,
                               11.612550017,
                               9.742647922,
                               7.779575692,
                               6.241633459,
                               4.530678113,
                               2.713659699};
    vector<double> dNchdeta_e = {0,           1.121689195, 0.856354796,
                                 0.685848455, 0.582627504, 0.498773083,
                                 0.415988997, 0.327417792, 0.261067034,
                                 0.186885663, 0.110000678};

    // input must be in the multiplicity range
    if (std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(),
                  multi_start) == end(dNchdeta_multibin))
        return {99, 99};
    if (std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(),
                  multi_end) == end(dNchdeta_multibin))
        return {99, 99};

    // special cases
    // Common case
    // Value
    vector<double>::iterator itr_left =
        find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_start);
    vector<double>::iterator itr_right =
        find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_end);
    int left = distance(dNchdeta_multibin.begin(), itr_left);
    int right = distance(dNchdeta_multibin.begin(), itr_right);

    int gap = right - left;

    double result = 0.;
    for (int i = 1; i < gap + 1; i++)
        result += dNchdeta[i + left] * (dNchdeta_multibin[i + left] -
                                        dNchdeta_multibin[i + left - 1]);

    result /= (multi_end - multi_start);
    returnarray.push_back(result);

    // Error
    double error = 0.;
    for (int i = 1; i < gap + 1; i++)
        error += pow(dNchdeta_e[i + left], 2);

    error = sqrt(error);
    returnarray.push_back(error);

    return returnarray;
}
TH1* ReturnExtremeHighHisto(TH1* hin) {
    // from AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C
    // input: histogram with sys error
    // output: maximum high error histogram.
    TH1* hout = (TH1*)hin->Clone(Form("%s_extremehigh", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        Double_t val = hin->GetBinContent(ibin + 1);
        Double_t err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, val + err);
    }
    return hout;
}
TH1* ReturnExtremeLowHisto(TH1* hin) {
    // from AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C
    // input: histogram with sys error
    // output: minimum low error histogram.
    TH1* hout = (TH1*)hin->Clone(Form("%s_extremelow", hin->GetName()));
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        Double_t val = hin->GetBinContent(ibin + 1);
        Double_t err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, val - err);
    }
    return hout;
}
TH1* ReturnExtremeHisto(TH1* hin, Float_t sign) {
    // from AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C
    // input: histogram with sys eror, sign(+1 or -1)
    // output maximum hard(+1)/soft(-1)error histogram.
    Double_t ptlow, pthigh;
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        ptlow = hin->GetBinLowEdge(ibin + 1);
        break;
    }
    for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        pthigh = hin->GetBinLowEdge(ibin + 2);
        break;
    }

    Double_t mean = hin->GetMean();
    Double_t maxdiff = 0.;
    TH1* hmax = NULL;
    for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {
        Double_t ptnode = hin->GetBinCenter(inode + 1);
        TH1* hout = (TH1*)hin->Clone(Form("%s_extremehard", hin->GetName()));

        for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
            if (hin->GetBinError(ibin + 1) <= 0.)
                continue;
            Double_t val = hin->GetBinContent(ibin + 1);
            Double_t err = hin->GetBinError(ibin + 1);
            Double_t cen = hin->GetBinCenter(ibin + 1);
            if (cen < ptnode)
                err *= -1. + (cen - ptlow) / (ptnode - ptlow);
            else
                err *= (cen - ptnode) / (pthigh - ptnode);

            hout->SetBinContent(ibin + 1, val + sign * err);
        }

        Double_t diff = TMath::Abs(mean - hout->GetMean());
        if (diff > maxdiff) {
            //      printf("found max at %f\n", ptnode);
            if (hmax)
                delete hmax;
            const char* nameout;
            if (sign < 0)
                nameout = "min";
            else
                nameout = "max";
            hmax = (TH1*)hout->Clone(
                Form("%s_extremehard_%s", hin->GetName(), nameout));
            maxdiff = diff;
        }
        delete hout;
    }
    return hmax;
}
TH1* ReturnRandom(TH1* hin) {
    // from AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C
    // input: histogram with stat error
    // output: histogram with randomly extracted value from
    //         N(contents, error)
    TH1* hout = (TH1*)hin->Clone("hout");
    hout->Reset();
    Double_t cont, err;
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        cont = hin->GetBinContent(ibin + 1);
        err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, gRandom->Gaus(cont, err));
        hout->SetBinError(ibin + 1, err);
    }
    return hout;
}

TH1* ReturnCoherentRandom(TH1* hin) {
    // from AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C
    // input: histogram with stat error
    // output: histogram with value + coherent randome error
    //         from N(0, 1) * Error
    if (!hin)
        return 0x0;
    TH1* hout = (TH1*)hin->Clone("hout");
    hout->Reset();
    Double_t cont, err, cohe;
    cohe = gRandom->Gaus(0., 1.);
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.)
            continue;
        cont = hin->GetBinContent(ibin + 1);
        err = hin->GetBinError(ibin + 1);
        hout->SetBinContent(ibin + 1, cont + cohe * err);
        hout->SetBinError(ibin + 1, err);
    }
    return hout;
}
TH1D* MakeHistfromArray(char const* name,
                        vector<double> dArray,
                        vector<double> eArray,
                        vector<double> ptbin,
                        const char* foption) {
    TString option = foption;
    double* ptbin_array = &ptbin[0];

    TH1D* htemp = new TH1D(Form("%s", name), "", ptbin.size() - 1, ptbin_array);
    for (int i = 0; i < dArray.size(); i++) {
        if (dArray.at(i) < 1e-9)
            cout << "skip zero bin" << endl;
        else
            htemp->SetBinContent(i + 1, dArray.at(i));
        if (!option.Contains("NOERROR"))
            htemp->SetBinError(i + 1, eArray.at(i));
    }
    return htemp;
}