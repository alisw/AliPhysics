#include <AliPWGHistoTools.h>
#include <AliFigure.h>
#include "YieldMean.C"
#include "AdditionalFunctions.h"

vector<TH1D*> GetSysSpectra(vector<vector<double>> multibin);
vector<TH1D*> GetSysSpectraNocor(vector<vector<double>> multibin);
vector<TH1D*> GetStatSpectra(vector<vector<double>> multibin);
TH1D* GetSpectrasys(double multi_start, double multi_end);
TH1D* GetSpectrasysNocor(double multi_start, double multi_end);
TH1D* GetSpectrastat(double multi_start, double multi_end);
TH1D *MakeHistfromArray(char const *name, vector<double> dArray, vector<double> eArray, vector<double> ptbin, const char* foption = "");
vector<double> GetdNdetawithError(double multi_start, double multi_end);
vector<double> GetPidNdetawithError(double multi_start, double multi_end);
TString finalfile = "AnalysisResults_Xi1530_systematic0010305070100.root";
TString workdirectory = "/Users/blim/alidock/Postprocessing/data/";
void DrawXi1530PhysicsPlots(){
    vector<vector<double>> multibin = {
        {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
        //{0, 10}, {10, 30}, {30, 50}, {50, 100}};
    //{0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}, {30, 100}};
    //{0, 100}};
    TString bininfo = "";
    for (auto const& bin : multibin)
        bininfo += Form("%.0f", bin[0]);
    bininfo += Form("%.0f", multibin.back()[1]);

    Int_t maxtrial = 1000;
    TString fitoption = "0q"; //default "0q"
    
    TH1* hout = new TH1D("hout", "", 9, 0, 9);
    TFile* resultfile =
        TFile::Open(Form("AnalysisResults_Xi1530_YieldMean_%s_try%d%s.root",
                         bininfo.Data(), maxtrial, fitoption.Data()),
                    "RECREATE");

    bool recalculate7TeV = false;

    vector<TH1D*> hspectra_sys = GetSysSpectra(multibin);
    vector<TH1D*> hspectra_sys_nocor = GetSysSpectraNocor(multibin);
    vector<TH1D*> hspectra_stat  = GetStatSpectra(multibin);
    // ch <dn/deta>
    vector<double> dNdetaAxis;
    vector<double> dNdetaAxis_e;
    vector<double> dNdetaAxis_e_half;
    vector<double> zeroerror;

    //pion <dn/dy>
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

    TF1* myLevy;
    // for memo, small
    TLatex* tm = new TLatex();
    tm->SetNDC();
    tm->SetTextSize(0.03);

    TCanvas* cFitResults = new TCanvas("cFitResults", "cFitResults", 1280, 720);
    gStyle->SetOptStat(0);
    cFitResults->Draw();
    cFitResults->SetLogy();

    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        cout << "Bin: " << multibin[imultibin][0] << " - "
             << multibin[imultibin][1] << endl;
        // testing!
        if(imultibin > 2){
            hspectra_sys_nocor[imultibin]->SetBinError(
                1, 5*hspectra_sys_nocor[imultibin]->GetBinError(1));
            hspectra_stat[imultibin]->SetBinError(
                1, 5*hspectra_stat[imultibin]->GetBinError(1));
        }
        myLevy = LevyTsallis("Levy", 1.5318, 15, 0.4, 1.5);
        auto hfinal =
            YieldMean(hspectra_stat[imultibin], hspectra_sys_nocor[imultibin],
                      myLevy, 0.0001, 10, 0.01, 0.1, fitoption.Data(), "log.root",
                      0.8, 8.8, maxtrial);
        cFitResults->cd();
        hspectra_sys_nocor[imultibin]->SetMarkerStyle(20);
        hspectra_sys_nocor[imultibin]->GetXaxis()->SetRangeUser(0, 10);
        hspectra_sys_nocor[imultibin]->SetMaximum(5e-1);
        hspectra_sys_nocor[imultibin]->SetMinimum(5e-9);
        hspectra_sys_nocor[imultibin]->Draw("E");
        myLevy->SetRange(0,10);
        myLevy->SetLineColor(kRed);
        myLevy->SetLineWidth(2);
        myLevy->Draw("same");
        tm->DrawLatex(0.62, 0.87,
                      Form("#bf{#chi^{2}/NDF: %.1f/%d}", myLevy->GetChisquare(),
                           myLevy->GetNDF()));
        tm->DrawLatex(
            0.75, 0.87,
            Form("#bf{0-0.8 bin ratio:%.2f}",
                 (myLevy->Integral(0, 0.8) / myLevy->Integral(0, 10))));
        cFitResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Spectrafit_%.2f-%.2f_%s_maxtrial_%d_yield.pdf",
                 workdirectory.Data(), multibin[imultibin][0],
                 multibin[imultibin][1], fitoption.Data(), maxtrial));
        hspectra_sys[imultibin]->GetXaxis()->SetRangeUser(0, 3.2);
        hspectra_sys[imultibin]->SetMaximum(
            2 * hspectra_sys[imultibin]->GetBinContent(1));
        hspectra_sys[imultibin]->SetMinimum(
            0.9 * hspectra_sys[imultibin]->GetBinContent(6));
        hspectra_sys[imultibin]->Draw("E");
        myLevy->SetRange(0, 3.2);
        myLevy->Draw("same");
        tm->DrawLatex(0.62, 0.87,
                      Form("#bf{#chi^{2}/NDF: %.1f/%d}", myLevy->GetChisquare(),
                           myLevy->GetNDF()));
        tm->DrawLatex(
            0.75, 0.87,
            Form("#bf{0-0.8 bin ratio:%.2f}",
                 (myLevy->Integral(0, 0.8) / myLevy->Integral(0, 10))));
        cFitResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Spectrafit_%.2f-%.2f_%s_maxtrial_%d_yield_zoom.pdf",
                 workdirectory.Data(), multibin[imultibin][0],
                 multibin[imultibin][1], fitoption.Data(), maxtrial));
        /*
        myLevy = LevyTsallis("Levy", 1.5318, 15, 0.4, 1.0);
        auto hfinal2 =
            YieldMean(hspectra_stat[imultibin], hspectra_sys_nocor[imultibin],
                      myLevy, 0.0, 10, 0.01, 0.1, fitoption.Data(), "log.root",
                      0.8, 8.8, maxtrial);
        cFitResults->cd();
        hspectra_sys_nocor[imultibin]->GetXaxis()->SetRangeUser(0, 10);
        hspectra_sys_nocor[imultibin]->SetMaximum(5e-1);
        hspectra_sys_nocor[imultibin]->SetMinimum(5e-9);
        hspectra_sys_nocor[imultibin]->Draw("E");
        myLevy->SetRange(0, 10);
        myLevy->SetLineColor(kRed);
        myLevy->SetLineWidth(2);
        myLevy->Draw("same");
        tm->DrawLatex(0.62, 0.87,
                      Form("#bf{#chi^{2}/NDF: %.1f/%d}", myLevy->GetChisquare(),
                           myLevy->GetNDF()));
        tm->DrawLatex(
            0.75, 0.87,
            Form("#bf{0-0.8 bin ratio:%.2f}",
                 (myLevy->Integral(0, 0.8) / myLevy->Integral(0, 10))));
        cFitResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Spectrafit_%.2f-%.2f_%s_maxtrial_%d_meanpt.pdf",
                 workdirectory.Data(), multibin[imultibin][0],
                 multibin[imultibin][1], fitoption.Data(), maxtrial));
        */

        yield.push_back(hfinal->GetBinContent(1));
        yield_state.push_back(hfinal->GetBinContent(2));
        yield_syshe.push_back(hfinal->GetBinContent(3));
        yield_sysle.push_back(hfinal->GetBinContent(4));
        meanpt.push_back(hfinal->GetBinContent(5));
        meanpt_state.push_back(hfinal->GetBinContent(6));
        meanpt_syshe.push_back(hfinal->GetBinContent(7));
        meanpt_sysle.push_back(hfinal->GetBinContent(8));
        extraYields.push_back(hfinal->GetBinContent(9));
        /*
        meanpt.push_back(hfinal2->GetBinContent(5));
        meanpt_state.push_back(hfinal2->GetBinContent(6));
        meanpt_syshe.push_back(hfinal2->GetBinContent(7));
        meanpt_sysle.push_back(hfinal2->GetBinContent(8));
        extraYields.push_back(hfinal->GetBinContent(9));
        */
        vector<double> temp = GetdNdetawithError(multibin[imultibin][0],multibin[imultibin][1]);
        dNdetaAxis.push_back(temp[0]);
        dNdetaAxis_e.push_back(temp[1]);
        dNdetaAxis_e_half.push_back(temp[1]/2);

        vector<double> temp2 = GetPidNdetawithError(multibin[imultibin][0],multibin[imultibin][1]);
        dNdeta_pi.push_back(temp2[0]);
        dNdeta_pi_e.push_back(temp2[1]);

        zeroerror.push_back(0);

        // save to outputfile
        resultfile->cd();
        hout->SetBinContent(kYield, yield[imultibin]);
        hout->SetBinContent(kYieldStat, yield_state[imultibin]);
        hout->SetBinContent(kYieldSysHi, yield_syshe[imultibin]);
        hout->SetBinContent(kYieldSysLo, yield_sysle[imultibin]);
        hout->SetBinContent(kMean, meanpt[imultibin]);
        hout->SetBinContent(kMeanStat, meanpt_state[imultibin]);
        hout->SetBinContent(kMeanSysHi, meanpt_syshe[imultibin]);
        hout->SetBinContent(kMeanSysLo, meanpt_sysle[imultibin]);

        hout->Write(
            Form("YieldMean_%.2f-%.2f",
                 multibin[imultibin][0],
                 multibin[imultibin][1]));
    }
    TGraphErrors* ge_stat = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &yield[0], &zeroerror[0], &yield_state[0]);
    TGraphAsymmErrors* ge_sys = new TGraphAsymmErrors(multibin.size(), &dNdetaAxis[0], &yield[0], &dNdetaAxis_e[0], &dNdetaAxis_e[0], &yield_sysle[0], &yield_syshe[0]);
    ge_stat->SetTitle("");
    ge_stat->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    ge_stat->GetYaxis()->SetTitle("< d#it{N}_{#Xi^{*0}}/dy >");
    ge_stat->GetXaxis()->SetRangeUser(0,30);
    //ge_stat->GetYaxis()->SetRangeUser(0,0.007);
    ge_stat->SetMinimum(0);
    ge_stat->SetMaximum(0.008);
    ge_stat->SetLineColor(2);
    ge_stat->SetMarkerColor(2);
    ge_stat->SetMarkerStyle(20);
    ge_stat->SetMarkerSize(0.5);
    ge_stat->SetFillColor(2);
    
    resultfile->cd();
    ge_stat->Write("gYield_stat");

    TGaxis::SetMaxDigits(3);
    ge_sys->SetTitle("");
    ge_sys->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    ge_sys->GetYaxis()->SetTitle("< d#it{N}_{#Xi^{*0}}/dy >");
    ge_sys->GetXaxis()->SetRangeUser(0,30);
    ge_sys->GetYaxis()->SetRangeUser(0,0.007);
    ge_sys->SetMinimum(0);
    ge_sys->SetMaximum(0.008);
    //ge_sys->SetFillStyle(3001);
    ge_sys->SetLineColor(2);
    ge_sys->SetMarkerColor(2);
    ge_sys->SetFillColor(0);
    ge_sys->GetXaxis()->SetLabelSize(0.05);
    ge_sys->GetXaxis()->SetTitleSize(0.05);
    ge_sys->GetYaxis()->SetLabelSize(0.05);
    ge_sys->GetYaxis()->SetTitleSize(0.05);
    ge_sys->GetYaxis()->SetTitleOffset(1.0);

    resultfile->cd();
    ge_sys->Write("gYield_syse");

    // 7TeV Results
    // Old data
    // from HEPDATA https://www.hepdata.net/record/ins1300380

    vector<double> ptbin2 = {0.8, 1.2, 1.6, 2.0, 2.4,
                        3.2, 4.0, 4.8, 5.6};
    vector<double> CorrectedYeild_7TeV2 = {0.00138,  0.00109, 0.00078,
                                          0.00046, 0.000226, 6.6e-05, 2.34e-05,
                                          8.1e-06};
    vector<double> CorrectedYeild_syserr_7TeV = {
        0.00017,  6.0e-05,  5.0e-05, 2.2e-05, 1.15e-05,
        3.73e-06, 1.96e-06, 4.54e-07};
    vector<double> CorrectedYeild_staterr_7TeV = {
        4.0e-05, 3.0e-05, 2.0e-05, 1.5e-05, 3.38e-06,
        2.08e-6, 8.11e-7, 9.09e-7};

    auto hSpectra_7TeV_syserr = MakeHistfromArray(
        "7TeV Spectra with systematic error", CorrectedYeild_7TeV2,
        CorrectedYeild_syserr_7TeV, ptbin2);
    auto hSpectra_7TeV_staterr = MakeHistfromArray(
        "7TeV Spectra with statistical error", CorrectedYeild_7TeV2,
        CorrectedYeild_staterr_7TeV, ptbin2);

    auto hfinal = YieldMean(hSpectra_7TeV_staterr, hSpectra_7TeV_syserr, myLevy, 0.0, 10, 0.01, 0.1);
    cout << "7TeV Yield: " << hfinal->GetBinContent(1) << " +- " << hfinal->GetBinContent(2) << " + " << hfinal->GetBinContent(3) << " - " << hfinal->GetBinContent(4) << endl;

    vector<double> x7 = {6.01};
    vector<double> x7e = {0.10};
     
    // HEP data
    vector<double> y7 = {2.56e-3};
    vector<double> y7e = {0.07e-3};
    vector<double> y7l = {0.37e-3};
    vector<double> y7h = {0.40e-3};

    // This macro
    vector<double> y7new = {hfinal->GetBinContent(1)};
    vector<double> y7enew = {hfinal->GetBinContent(2)};
    vector<double> y7lnew = {hfinal->GetBinContent(4)};
    vector<double> y7hnew = {hfinal->GetBinContent(3)};    

    vector<double> pt7 = {1.31};
    vector<double> pt7e = {0.02};
    vector<double> pt7l = {0.09};
    vector<double> pt7h = {0.09};

    //this macro

    vector<double> pt7new = {hfinal->GetBinContent(5)};
    vector<double> pt7enew = {hfinal->GetBinContent(6)};
    vector<double> pt7lnew = {hfinal->GetBinContent(8)};
    vector<double> pt7hnew = {hfinal->GetBinContent(7)};    

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

    //TGraphErrors* ge_stat_7TeV = new TGraphErrors(1, &x7[0], &y7[0], &x7e[0], &y7e[0]);
    TGraphErrors* ge_stat_7TeV = new TGraphErrors(1, &x7[0], &y7new[0], &x7e[0], &y7enew[0]);
    ge_stat_7TeV->SetLineColor(4);
    ge_stat_7TeV->SetMarkerColor(4);
    ge_stat_7TeV->SetMarkerStyle(20);
    ge_stat_7TeV->SetMarkerSize(0.5);
    ge_stat_7TeV->SetFillColor(4);

    TGraphErrors* ge_stat_7TeV_HEP = new TGraphErrors(1, &x7[0], &y7[0], &x7e[0], &y7e[0]);
    ge_stat_7TeV_HEP->SetLineColor(1);
    ge_stat_7TeV_HEP->SetMarkerColor(1);
    ge_stat_7TeV_HEP->SetMarkerStyle(20);
    ge_stat_7TeV_HEP->SetMarkerSize(0.5);
    ge_stat_7TeV_HEP->SetFillColor(1);

    resultfile->cd();
    ge_stat_7TeV_HEP->Write("gYield7TeV_stat");

    TGraphAsymmErrors* ge_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &y7new[0], &x7e[0], &x7e[0], &y7lnew[0], &y7hnew[0]);
    //TGraphAsymmErrors* ge_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &y7[0], &x7e[0], &x7e[0], &y7l[0], &y7h[0]);
    ge_sys_7TeV->SetLineColor(4);
    ge_sys_7TeV->SetMarkerColor(4);
    ge_sys_7TeV->SetFillColor(0);


    TGraphAsymmErrors* ge_sys_7TeV_HEP = new TGraphAsymmErrors(1, &x7[0], &y7[0], &x7e[0], &x7e[0], &y7l[0], &y7h[0]);
    //TGraphAsymmErrors* ge_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &y7[0], &x7e[0], &x7e[0], &y7l[0], &y7h[0]);
    ge_sys_7TeV_HEP->SetLineColor(1);
    ge_sys_7TeV_HEP->SetMarkerColor(1);
    ge_sys_7TeV_HEP->SetFillColor(0);

    resultfile->cd();
    ge_stat_7TeV_HEP->Write("gYield7TeV_syse");

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
    cCanvas->SetTickx();
    cCanvas->Draw();
    cCanvas->cd();
    //cCanvas->SetLogy(true);
    ge_sys->Draw("a5");
    ge_stat->Draw("P");
    ge_sys_7TeV_HEP->Draw("5");
    ge_stat_7TeV_HEP->Draw("P");
    if(recalculate7TeV){
        ge_sys_7TeV->Draw("5");
        ge_stat_7TeV->Draw("P");
    }

    TText* fStatusPad  = new TText(0.5, 0.15, "work in progress");
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
    if (recalculate7TeV)
        legendyield->AddEntry(ge_sys_7TeV, "pp 7 TeV[Re-calculated]", "F");
    legendyield->AddEntry(ge_sys_7TeV_HEP, "pp 7 TeV[Paper]", "F");
    legendyield->Draw();
    cCanvas->SaveAs(
        Form("%s1_Multi_0.00-100.00_Default1/yield_%s_maxtrial_%d.pdf",
             workdirectory.Data(), fitoption.Data(), maxtrial));
    resultfile->cd();
    cCanvas->Write("cYield");

    // Mean pt
    TGraphErrors* gpt_stat = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &meanpt[0], &zeroerror[0], &meanpt_state[0]);
    TGraphAsymmErrors* gpt_sys = new TGraphAsymmErrors(multibin.size(), &dNdetaAxis[0], &meanpt[0], &dNdetaAxis_e[0], &dNdetaAxis_e[0], &meanpt_sysle[0], &meanpt_syshe[0]);
    gpt_stat->SetTitle("");
    gpt_stat->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    gpt_stat->GetYaxis()->SetTitle("< #it{p}_{T} > (GeV/c)");
    gpt_stat->GetXaxis()->SetRangeUser(0,30);
    gpt_stat->GetYaxis()->SetRangeUser(1,1.7);
    gpt_stat->SetLineColor(2);
    gpt_stat->SetMarkerColor(2);
    gpt_stat->SetMarkerStyle(20);
    gpt_stat->SetMarkerSize(0.5);
    gpt_stat->SetFillColor(2);

    resultfile->cd();
    gpt_stat->Write("gMeanpT_stat");

    gpt_sys->SetTitle("");
    gpt_sys->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    gpt_sys->GetYaxis()->SetTitle("< #it{p}_{T} > (GeV/c)");
    gpt_sys->GetXaxis()->SetRangeUser(0,30);
    gpt_sys->GetYaxis()->SetRangeUser(1,1.7);
    //ge_sys->SetFillStyle(3001);
    gpt_sys->SetLineColor(2);
    gpt_sys->SetMarkerColor(2);
    gpt_sys->SetFillColor(0);
    gpt_sys->GetXaxis()->SetLabelSize(0.05);
    gpt_sys->GetXaxis()->SetTitleSize(0.05);
    gpt_sys->GetYaxis()->SetLabelSize(0.05);
    gpt_sys->GetYaxis()->SetTitleSize(0.05);
    gpt_sys->GetYaxis()->SetTitleOffset(1.0);

    resultfile->cd();
    gpt_sys->Write("gMeanpT_syse");

    TGraphErrors* gpt_stat_7TeV = new TGraphErrors(1, &x7[0], &pt7new[0], &x7e[0], &pt7enew[0]);
    gpt_stat_7TeV->SetLineColor(4);
    gpt_stat_7TeV->SetMarkerColor(4);
    gpt_stat_7TeV->SetMarkerStyle(20);
    gpt_stat_7TeV->SetMarkerSize(0.5);
    gpt_stat_7TeV->SetFillColor(4);

    TGraphErrors* gpt_stat_7TeV_HEP = new TGraphErrors(1, &x7[0], &pt7[0], &x7e[0], &pt7e[0]);
    gpt_stat_7TeV_HEP->SetLineColor(1);
    gpt_stat_7TeV_HEP->SetMarkerColor(1);
    gpt_stat_7TeV_HEP->SetMarkerStyle(20);
    gpt_stat_7TeV_HEP->SetMarkerSize(0.5);
    gpt_stat_7TeV_HEP->SetFillColor(1);

    resultfile->cd();
    gpt_stat_7TeV_HEP->Write("gMeanpT7TeV_stat");

    TGraphAsymmErrors* gpt_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &pt7new[0], &x7e[0], &x7e[0], &pt7lnew[0], &pt7hnew[0]);
    gpt_sys_7TeV->SetLineColor(4);
    gpt_sys_7TeV->SetMarkerColor(4);
    gpt_sys_7TeV->SetFillColor(0);

    TGraphAsymmErrors* gpt_sys_7TeV_HEP = new TGraphAsymmErrors(1, &x7[0], &pt7[0], &x7e[0], &x7e[0], &pt7l[0], &pt7h[0]);
    gpt_sys_7TeV_HEP->SetLineColor(1);
    gpt_sys_7TeV_HEP->SetMarkerColor(1);
    gpt_sys_7TeV_HEP->SetFillColor(0);

    resultfile->cd();
    gpt_sys_7TeV_HEP->Write("gMeanpT7TeV_syse");

    TCanvas* cCanvas2 = new TCanvas("cCanvas2", "cCanvas2", 1200, 720);
    //gStyle->SetOptTitle(0);
    cCanvas2->SetTickx();
    cCanvas2->Draw();
    cCanvas2->cd();
    //cCanvas->SetLogy(true);
    gpt_sys->Draw("a5");
    gpt_stat->Draw("P");
    gpt_sys_7TeV_HEP->Draw("5");
    gpt_stat_7TeV_HEP->Draw("P");
    if(recalculate7TeV){
        gpt_sys_7TeV->Draw("5");
        gpt_stat_7TeV->Draw("P");
    }


    fStatusPad->Draw();

    legendyield->Draw();
    cCanvas2->SaveAs(
        Form("%s1_Multi_0.00-100.00_Default1/meanpt_%s_maxtrial_%d.pdf",
             workdirectory.Data(), fitoption.Data(), maxtrial));
    resultfile->cd();
    cCanvas2->Write("cMeanpT");
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        cout << "Multi bin: " << multibin[imultibin][0] << " - " << multibin[imultibin][1] << ", Yield: " << yield[imultibin] << " +- " << yield_state[imultibin] << " + " << yield_syshe[imultibin] << " - " << yield_sysle[imultibin] << ", mean pT: " << meanpt[imultibin] << " +- " << meanpt_state[imultibin] << " + " << meanpt_syshe[imultibin] << " - " << meanpt_sysle[imultibin]  << "| Extra yields: " << extraYields[imultibin] << endl; 
    }

    resultfile->Close();
}
vector<TH1D*> GetSysSpectra(vector<vector<double>> multibin){
    vector<TH1D*> buffer = {};
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        auto temp_spectra_multi2 = GetSpectrasys(multibin[imultibin][0],
                                    multibin[imultibin][1]);
        buffer.push_back(temp_spectra_multi2);
    }
    return buffer;
}
vector<TH1D*> GetSysSpectraNocor(vector<vector<double>> multibin) {
    vector<TH1D*> buffer = {};
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        auto temp_spectra_multi =
            GetSpectrasysNocor(multibin[imultibin][0], multibin[imultibin][1]);
        buffer.push_back(temp_spectra_multi);
    }
    return buffer;
}
vector<TH1D*> GetStatSpectra(vector<vector<double>> multibin){
    vector<TH1D*> buffer = {};
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        auto temp_spectra_multi = GetSpectrastat(multibin[imultibin][0],
                                    multibin[imultibin][1]);
        buffer.push_back(temp_spectra_multi);
    }
    return buffer;
}
TH1D* GetSpectrasys(double multi_start, double multi_end){
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_sys", multi_start, multi_end));

    return hr;
}
TH1D* GetSpectrasysNocor(double multi_start, double multi_end) {
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_sys_noCorrelation", multi_start, multi_end));

    return hr;
}
TH1D* GetSpectrastat(double multi_start, double multi_end){
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(Form("hSpectra_%.2f_%.2f_stat",
                               multi_start,
                               multi_end));

    return hr;
}
vector<double> GetdNdetawithError(double multi_start, double multi_end){
    // Return dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}

    vector<double> returnarray;

    //--dNdeta histogram
    // Ref: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/
    //      ReferenceMult#Multiplicity_dependent_pp_at_AN2
    // LHC16k data.
    // Error was asym error, so choosed bigger error for sym error.

    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta = 
    {0, 26.02, 20.02, 16.17, 13.77, 12.04, 10.02, 7.95, 6.32, 4.50, 2.55};
    vector<double> dNchdeta_e = 
    {0,  0.35,  0.27,  0.22,  0.19,  0.17,  0.14, 0.11, 0.09, 0.07, 0.04};

    // input must be in the multiplicity range
    if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_start) == end(dNchdeta_multibin))
        return {99,99};
    if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_end) == end(dNchdeta_multibin))
        return {99,99};

    // special cases
    if((multi_start == 0) && (multi_end == 0.01)){
        returnarray = {35.37, 0.92};
    }
    else if((multi_start == 0.01) && (multi_end == 0.1)){
        returnarray = {30.89, 0.57};
    }
    else if((multi_start == 0) && (multi_end == 5)){
        returnarray = {21.20, 0.28};
    }
    else if((multi_start == 0) && (multi_end == 100)){
        returnarray = {6.94, 0.10};
    }
    else if((multi_start == 0.1) && (multi_end == 0.5)){
        returnarray = {26.96, 0.37};
    }
    else if((multi_start == 0.5) && (multi_end == 1)){
        returnarray = {24.23, 0.36};
    }
    // Common case
    else{
        // Value
        vector<double>::iterator itr_left = find(dNchdeta_multibin.begin(),
                                dNchdeta_multibin.end(),
                                multi_start);
        vector<double>::iterator itr_right = find(dNchdeta_multibin.begin(),
                                dNchdeta_multibin.end(),
                                multi_end);
        int left = distance(dNchdeta_multibin.begin(), itr_left);
        int right = distance(dNchdeta_multibin.begin(), itr_right);

        int gap = right - left;

        double result = 0.;
        for(int i = 1; i < gap+1; i++)
            result += dNchdeta[i+left]*(dNchdeta_multibin[i+left] - dNchdeta_multibin[i+left-1]);
            
        result /= (multi_end - multi_start);
        returnarray.push_back(result);

        // Error
        double error = 0.;
        for(int i = 1; i < gap+1; i++)
            error += pow( dNchdeta_e[i+left], 2); 
            
        error = sqrt(error);
        returnarray.push_back(error);
    }

    return returnarray;
}
vector<double> GetPidNdetawithError(double multi_start, double multi_end){
    // Return pion's dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}

    vector<double> returnarray;

    //--dNdeta histogram of pion+ + pion-
    // Ref: given by Anders, need to update ref.
    // Error was stat+sys error, so choosed bigger error(sym) error.

    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta = 
    {0, 24.605455025, 19.016207266, 15.477779961, 13.241832514, 11.612550017, 9.742647922, 7.779575692,  6.241633459, 4.530678113, 2.713659699};
    vector<double> dNchdeta_e = 
    {0,  1.121689195,  0.856354796,  0.685848455,  0.582627504,   0.498773083,  0.415988997, 0.327417792, 0.261067034, 0.186885663, 0.110000678};

    // input must be in the multiplicity range
    if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_start) == end(dNchdeta_multibin))
        return {99,99};
    if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_end) == end(dNchdeta_multibin))
        return {99,99};

    // special cases
    // Common case
    // Value
    vector<double>::iterator itr_left = find(dNchdeta_multibin.begin(),
                            dNchdeta_multibin.end(),
                            multi_start);
    vector<double>::iterator itr_right = find(dNchdeta_multibin.begin(),
                            dNchdeta_multibin.end(),
                            multi_end);
    int left = distance(dNchdeta_multibin.begin(), itr_left);
    int right = distance(dNchdeta_multibin.begin(), itr_right);

    int gap = right - left;

    double result = 0.;
    for(int i = 1; i < gap+1; i++)
        result += dNchdeta[i+left]*(dNchdeta_multibin[i+left] - dNchdeta_multibin[i+left-1]);
        
    result /= (multi_end - multi_start);
    returnarray.push_back(result);

    // Error
    double error = 0.;
    for(int i = 1; i < gap+1; i++)
        error += pow( dNchdeta_e[i+left], 2); 
        
    error = sqrt(error);
    returnarray.push_back(error);
    

    return returnarray;
}
TH1* ReturnExtremeHisto(TH1* hin, Float_t sign) {
    // from YieldMean.C
    // input: histogram with sys+stat eror, sign(+1 or -1)
    // output maximum error histogram.
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
            hmax = (TH1*)hout->Clone("hmax");
            maxdiff = diff;
        }
        delete hout;
    }
    return hmax;
}
TH1D *MakeHistfromArray(char const *name, vector<double> dArray, vector<double> eArray, vector<double> ptbin, const char* foption)
{
    TString option = foption;
    double *ptbin_array = &ptbin[0];

    TH1D *htemp = new TH1D(Form("%s", name), "", ptbin.size() - 1, ptbin_array);
    for (int i = 0; i < dArray.size(); i++)
    {
        if(dArray.at(i) < 1e-9)
            cout << "skip zero bin" << endl;
        else
            htemp->SetBinContent(i + 1, dArray.at(i));
        if(!option.Contains("NOERROR"))htemp->SetBinError(i + 1, eArray.at(i));
    }
    return htemp;
}