#include <AliPWGHistoTools.h>
#include <AliFigure.h>
#include <AliPWGFunc.h>
#include "YieldMean.C"
#include "AdditionalFunctions.h"
#include "DrawingHelper.C"

vector<TH1D*> GetSysSpectra(vector<vector<double>> multibin);
vector<TH1D*> GetSysSpectraNocor(vector<vector<double>> multibin);
vector<TH1D*> GetStatSpectra(vector<vector<double>> multibin);
TH1D* GetSpectrasys(double multi_start, double multi_end);
TH1D* GetSpectrasysNocor(double multi_start, double multi_end);
TH1D* GetSpectrastat(double multi_start, double multi_end);
vector<double> GetPidNdetawithError(double multi_start, double multi_end, int errortype = 1);
vector<double> GetXidNdetawithError(double multi_start, double multi_end, int errortype = 1);
vector<double> GetYieldError(double multi_start, double multi_end);
vector<double> GetMeanPtError(double multi_start, double multi_end);
TString GetFinalFileName(double multi_start, double multi_end);
TString finalfile;
//TString finalfile = "AnalysisResults_Xi1530_systematic0020507050100.root";
TString workdirectory = "/Users/blim/alidock/Postprocessing/data/";
enum {kFitExpPt=1, kFitLevi, fFitExpMt, kFitBoltzmann, kFitBlastWave, kFitBoseEinstein, kFitFermiDirac};
vector<TString> functions = {"", "kFitExpPt", "kFitLevi", "fFitExpMt", "kFitBoltzmann", "kFitBlastWave", "kFitBoseEinstein", "kFitFermiDirac"};
Int_t maxtrial = 10000;
TString fitoption = "0qEI"; //default "0q"
double mass = 1.5318;
void DrawXi1530PhysicsPlots(int fitFunc = 1){
    TCanvas* cResults = new TCanvas("cResults", "cResults", 960, 720);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cResults->SetTickx();
    cResults->SetTicky();
    //cResults->SetLogy();
    cResults->SetTopMargin(0.05);
    cResults->SetLeftMargin(0.10);
    //cSigbkg->SetBottomMargin(0.01);
    cResults->SetRightMargin(0.01);
    cResults->SetFillStyle(0);
    cResults->Draw();

    vector<vector<double>> multibin = {
        //{0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
        // {0, 20}, {20, 50}, {50, 70}, {70, 100}};
        // {0, 10}, {10, 30}, {30, 50}, {50, 100}};
        {0, 0.01}, {0.01, 0.05}, {0.05, 0.1}, {0, 10},
        {10, 30},  {30, 50},     {50, 70},    {70, 100}};
    //{0, 0.01}, {0.01, 0.1}, {0, 0}, {0, 100}, {0, 10}, {10, 30}, {30, 50},
    //{50, 70}, {70, 100}}; {0, 100}};


    
    /*
    TFile* fReweighting =
        TFile::Open("AnalysisResults_Xi1530_efficiencyReweighing.root","READ");
    
    vector<TH1D*> hdata_reweight;
    hdata_reweight.push_back((TH1D*)fReweighting->Get("_i0")); 
    hdata_reweight.push_back((TH1D*)fReweighting->Get("_i1")); 
    hdata_reweight.push_back((TH1D*)fReweighting->Get("_i2")); 
    hdata_reweight.push_back((TH1D*)fReweighting->Get("_i3"));

    vector<TH1D*> hinput_reweight;
    hinput_reweight.push_back((TH1D*)fReweighting->Get("inputSys001001Type006006Cent001001_i0")); 
    hinput_reweight.push_back((TH1D*)fReweighting->Get("inputSys001001Type006006Cent001001_i1")); 
    hinput_reweight.push_back((TH1D*)fReweighting->Get("inputSys001001Type006006Cent001001_i2")); 
    hinput_reweight.push_back((TH1D*)fReweighting->Get("inputSys001001Type006006Cent001001_i3")); 


    vector<TH1D*> hrecon_reweight;
    hrecon_reweight.push_back((TH1D*)fReweighting->Get("reconSys001001Type004004Cent001001_i0")); 
    hrecon_reweight.push_back((TH1D*)fReweighting->Get("reconSys001001Type004004Cent001001_i1")); 
    hrecon_reweight.push_back((TH1D*)fReweighting->Get("reconSys001001Type004004Cent001001_i2")); 
    hrecon_reweight.push_back((TH1D*)fReweighting->Get("reconSys001001Type004004Cent001001_i3")); 

    vector<TH1D*> heffi_reweight;
    for (int bin = 0; bin < hinput_reweight.size(); bin++){
        auto temp = (TH1D*)fReweighting->Get(Form("inputSys001001Type006006Cent001001_rebin_i%d",bin));
        auto temp2 = (TH1D*)fReweighting->Get(Form("reconSys001001Type004004Cent001001_rebin_i%d",bin));
        temp2->Divide(temp);
        heffi_reweight.push_back(temp2);
    }


    vector<TF1*> fLevy_fit;
    fLevy_fit.push_back((TF1*)fReweighting->Get("Levy_i1")); 
    fLevy_fit.push_back((TF1*)fReweighting->Get("Levy_i2")); 

    vector<TH1D*> hCorrectionfactor_rewight;
    hCorrectionfactor_rewight.push_back((TH1D*)fReweighting->Get("_correction_i1")); 
    hCorrectionfactor_rewight.push_back((TH1D*)fReweighting->Get("_correction_i2")); 
    hCorrectionfactor_rewight.push_back((TH1D*)fReweighting->Get("_correction_i3")); 
    
    cResults->cd();
    TH1D* hinput_temp = (TH1D*)hinput_reweight[0]->Clone();
    hinput_temp->SetMarkerStyle(20);
    hinput_temp->SetMarkerColor(kOrange + 10);
    hinput_temp->Scale(1e-7);
    hinput_temp->GetXaxis()->SetRangeUser(0, 10);
    hinput_temp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hinput_temp->SetMinimum(1e-6);
    hinput_temp->Draw();

    hdata_reweight[0]->SetMarkerStyle(20);
    hdata_reweight[0]->SetLineColor(kBlack);
    hdata_reweight[0]->SetMarkerColor(kBlack);
    hdata_reweight[0]->SetFillColor(0);
    hdata_reweight[0]->Draw("same");
    

    fLevy_fit[0]->SetLineColor(kBlack);
    fLevy_fit[0]->Draw("same");
    
    hrecon_reweight[0]->Scale(1e-7);
    hrecon_reweight[0]->SetMarkerStyle(24);
    hrecon_reweight[0]->SetMarkerColor(kOrange + 10);
    hrecon_reweight[0]->Draw("same");

    auto legendweight1 = new TLegend(.6, .35, .85, .6);
    legendweight1->SetBorderSize(0);
    legendweight1->SetFillStyle(0);
    legendweight1->AddEntry(hinput_temp, "MC Generated(x10^{-8})", "P");
    legendweight1->AddEntry(hrecon_reweight[0], "MC Reconstructed(x10^{-8})", "P");
    legendweight1->AddEntry(hdata_reweight[0], "Real data", "P");
    legendweight1->AddEntry(fLevy_fit[0], "Real data fit(Levy)", "L");
    legendweight1->Draw();

    cResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Reweight_Allinputs.pdf",
                 workdirectory.Data()));
    
    // Second
    auto legendweight2 = new TLegend(.6, .35, .85, .6);
    legendweight2->SetBorderSize(0);
    legendweight2->SetFillStyle(0);
    legendweight2->AddEntry(hinput_reweight[0], "Unweightend", "P");

    int sysColorPallet_reweight = GetSerialColors(hinput_reweight.size());
    int iloop = 1;
    for (auto const& hInput : hinput_reweight){
        if(iloop == 1) hInput->Scale(1e-7);
        hInput->SetMarkerStyle(iloop+19);
        hInput->SetMinimum(1e-8);
        hInput->GetXaxis()->SetRangeUser(0, 10);
        hInput->SetMarkerColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);

        if(iloop > 1)legendweight2->AddEntry(hInput, Form("Iteration %d", iloop-1), "P");

        iloop++;
    }
    hinput_reweight[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hinput_reweight[0]->Draw();

    for (auto const& hInput : hinput_reweight)
        hInput->Draw("same");
    
    legendweight2->Draw();

    cResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Reweight_inputiterate.pdf",
                 workdirectory.Data()));

    // Third
    heffi_reweight[0]->SetMarkerStyle(5);
    heffi_reweight[0]->SetMarkerColor(kOrange + 10);
    heffi_reweight[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    heffi_reweight[0]->Draw();
    auto legendweight3 = new TLegend(.6, .35, .85, .6);
    legendweight3->SetBorderSize(0);
    legendweight3->SetFillStyle(0);
    legendweight3->AddEntry(heffi_reweight[0], "Unweightend", "P");
    iloop = 1;
    for (auto const& hInput : heffi_reweight){
        hInput->SetMarkerStyle(iloop+19);
        hInput->SetMarkerColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);
        hInput->SetLineColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);
        hInput->Draw("same");
        if(iloop > 1) legendweight3->AddEntry(hInput, Form("Iteration %d", iloop-1), "P");
        iloop++;
    }
    legendweight3->Draw();


    cResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Reweight_effi.pdf",
                 workdirectory.Data()));
    // Fourth
    hdata_reweight[0]->SetMarkerStyle(5);
    hdata_reweight[0]->SetMarkerColor(kOrange + 10);
    hdata_reweight[0]->GetXaxis()->SetRangeUser(0, 4.8);
    hdata_reweight[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hdata_reweight[0]->Draw();
    auto legendweight4 = new TLegend(.6, .6, .85, .85);
    legendweight4->SetBorderSize(0);
    legendweight4->SetFillStyle(0);
    legendweight4->AddEntry(hdata_reweight[0], "Unweightend", "P");
    iloop = 1;
    for (auto const& hInput : hdata_reweight){
        hInput->SetMarkerStyle(iloop+19);
        hInput->SetMarkerColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);
        hInput->SetLineColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);
        hInput->GetXaxis()->SetRangeUser(0, 8.8);
        hInput->Draw("same");
        if(iloop > 1) legendweight4->AddEntry(hInput, Form("Iteration %d", iloop-1), "P");
        iloop++;
    }
    legendweight4->Draw();
        

    cResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Reweight_datachange.pdf",
                 workdirectory.Data()));

    // Fifth
    cResults->SetLogy(false);
    //hCorrectionfactor_rewight[0]->SetMarkerStyle(5);
    hCorrectionfactor_rewight[0]->SetMarkerColor(kOrange + 10);
    hCorrectionfactor_rewight[0]->SetLineColor(kOrange + 10);
    hCorrectionfactor_rewight[0]->GetXaxis()->SetRangeUser(0, 8.8);
    hCorrectionfactor_rewight[0]->SetMaximum(1.1);
    hCorrectionfactor_rewight[0]->SetMinimum(0.9);
    hCorrectionfactor_rewight[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    hCorrectionfactor_rewight[0]->Draw("E");
    auto legendweight5 = new TLegend(.6, .6, .85, .85);
    legendweight5->SetBorderSize(0);
    legendweight5->SetFillStyle(0);
    legendweight5->AddEntry(hCorrectionfactor_rewight[0], "Iteration 1", "PL");
    iloop = 1;
    for (auto const& hInput : hCorrectionfactor_rewight){
        //hInput->SetMarkerStyle(iloop+19);
        hInput->SetMarkerColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);
        hInput->SetLineColor(sysColorPallet_reweight + hinput_reweight.size() - iloop);
        hInput->GetXaxis()->SetRangeUser(0, 8.8);
        hInput->Draw("E same");
        if(iloop > 1) legendweight5->AddEntry(hInput, Form("Iteration %d", iloop), "PL");
        iloop++;
    }
    TF1* line1 = new TF1("line1", "1", -1, 15);
    line1->SetLineColor(1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");
    legendweight5->Draw();
        

    cResults->SaveAs(
            Form("%s1_Multi_0.00-100.00_Default1/"
                 "Reweight_factor.pdf",
                 workdirectory.Data()));
    gSystem->Exit(1);
    */
    

    // Results
    vector<double> dNdetaAxis;
    vector<double> dNdetaAxis_e;
    vector<double> dNdetaAxis_e_half;
    vector<double> zeroerror;

    //pion <dn/dy>
    vector<double> dNdeta_pi;
    vector<double> dNdeta_pi_e;
    vector<double> dNdeta_pi_sys;

    //Xi <dn/dy>
    vector<double> dNdeta_Xi;
    vector<double> dNdeta_Xi_state; // total sys
    vector<double> dNdeta_Xi_sys; // total sys
    vector<double> dNdeta_Xi_sys_uncor; // total sys

    vector<double> yield;
    vector<double> yield_state;
    vector<double> yield_syshe;
    vector<double> yield_sysle;
    vector<double> yield_sys_cor;


    vector<double> meanpt;
    vector<double> meanpt_state;
    vector<double> meanpt_syshe;
    vector<double> meanpt_sysle;
    vector<double> meanpt_sys_cor;
    vector<double> extraYields;

    // output file
    TString bininfo = "";
    for (auto const& bin : multibin){
        if(bin[1] < 0.5){
            bininfo += Form("%.2f", bin[0]);
            bininfo += "-";
            bininfo += Form("%.2f", bin[1]);
        }
        else {
            bininfo += Form("%.0f", bin[0]);
            bininfo += "-";
            bininfo += Form("%.0f", bin[1]);
        }
        bininfo += "_";
    }
    bininfo += "bins";
    
    TH1* hout = new TH1D("hout", "", 9, 0, 9);
    TFile* resultfile =
        TFile::Open(Form("AnalysisResults_Xi1530_PhysicsResult_%s.root",
                         bininfo.Data()),
                        "RECREATE");

    bool recalculate7TeV = false;
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        auto tempdndeta = GetdNdetawithError(multibin[imultibin][0],multibin[imultibin][1]);
        auto tempdndetapi = GetPidNdetawithError(multibin[imultibin][0],multibin[imultibin][1],1);
        auto tempdndetapi_sys = GetPidNdetawithError(multibin[imultibin][0],multibin[imultibin][1],2);
        auto tempdndetaXi = GetXidNdetawithError(multibin[imultibin][0],multibin[imultibin][1],1);
        auto tempdndetaXi_sys = GetXidNdetawithError(multibin[imultibin][0],multibin[imultibin][1],2);
        auto tempdndetaXi_uncor = GetXidNdetawithError(multibin[imultibin][0],multibin[imultibin][1],3);

        dNdetaAxis.push_back(tempdndeta[0]);
        dNdetaAxis_e.push_back(tempdndeta[1]);
        dNdetaAxis_e_half.push_back(tempdndeta[1]);
        zeroerror.push_back(0);

        dNdeta_pi.push_back(tempdndetapi[0]);
        dNdeta_pi_e.push_back(tempdndetapi[1]);
        dNdeta_pi_sys.push_back(tempdndetapi_sys[1]);

        dNdeta_Xi.push_back(tempdndetaXi[0]);
        dNdeta_Xi_state.push_back(tempdndetaXi[1]);
        dNdeta_Xi_sys.push_back(tempdndetaXi_sys[1]);
        dNdeta_Xi_sys_uncor.push_back(tempdndetaXi_uncor[1]);


        auto tempyield = GetYieldError(multibin[imultibin][0],multibin[imultibin][1]);
        yield.push_back(tempyield[0]);
        yield_state.push_back(tempyield[1]);

        double tempfinalsyserror = 0;
        tempfinalsyserror += pow(tempyield[2],2);
        tempfinalsyserror += pow((tempyield[3]+tempyield[4])/2,2); // avg(high+low)

        double averageerror_correl = sqrt(tempfinalsyserror);
        double averageerror_uncorr = 0.0538*tempyield[0];
        tempfinalsyserror += pow(averageerror_uncorr,2);
        double averageerror = sqrt(tempfinalsyserror);
        cout << "yield check: " << tempyield[0] << ", " << tempyield[1] << "("
             << tempyield[1] / tempyield[0] * 100 << "%), " << tempyield[2]
             << "(" << tempyield[2] / tempyield[0] * 100 << "%), "
             << tempyield[3] << "(" << tempyield[3] / tempyield[0] * 100
             << "%), " << tempyield[4] << "("
             << tempyield[4] / tempyield[0] * 100 << "%)" << endl;

        yield_syshe.push_back(averageerror);
        yield_sysle.push_back(averageerror);
        yield_sys_cor.push_back(averageerror_uncorr);

        auto tempmeanpt = GetMeanPtError(multibin[imultibin][0],multibin[imultibin][1]);
        meanpt.push_back(tempmeanpt[0]);
        meanpt_state.push_back(tempmeanpt[1]);

        double mtempfinalsyserror = 0;
        mtempfinalsyserror += pow(tempmeanpt[2],2);
        mtempfinalsyserror += pow((tempmeanpt[3]+tempmeanpt[4])/2,2); // avg(high+low)

        double averagpteerror_correl = sqrt(mtempfinalsyserror);
        double averagpteerror_uncorr = 0.0538*tempmeanpt[0];
        mtempfinalsyserror += pow(averagpteerror_uncorr,2);
        double averagpteerror = sqrt(mtempfinalsyserror);

        cout << "mpT check: " << tempmeanpt[0] << ", " << tempmeanpt[1] << "("
             << tempmeanpt[1] / tempmeanpt[0] * 100 << "%), " << tempmeanpt[2]
             << "(" << tempmeanpt[2] / tempmeanpt[0] * 100 << "%), "
             << tempmeanpt[3] << "(" << tempmeanpt[3] / tempmeanpt[0] * 100
             << "%), " << tempmeanpt[4] << "("
             << tempmeanpt[4] / tempmeanpt[0] * 100 << "%)" << endl;

        meanpt_syshe.push_back(averagpteerror);
        meanpt_sysle.push_back(averagpteerror);
        meanpt_sys_cor.push_back(averagpteerror_uncorr);

        cout << Form("%.2f", multibin[imultibin][0]) << " - "
             << Form("%.2f", multibin[imultibin][1]) << " ("
             << dNdetaAxis[imultibin] << ") "
             << "| " << yield[imultibin]*1000 << " +- " << yield_state[imultibin]*1000
             << " (" << 100 * (yield_state[imultibin] / yield[imultibin])
             << "%%)"
             << " + " << yield_syshe[imultibin]*1000 << " ("
             << 100 * (yield_syshe[imultibin] / yield[imultibin]) << "%%)"
             << " - " << yield_sysle[imultibin]*1000 << " ("
             << 100 * (yield_sysle[imultibin] / yield[imultibin]) << "%%)"
             << " +- " << averageerror_uncorr*1000 << " ("
             << 100 * (averageerror_uncorr / yield[imultibin]) << "%%) "
             << (tempyield[5]/tempyield[0])*100 << "%"

             << "| " << meanpt[imultibin] << " +- " << meanpt_state[imultibin]
             << " (" << 100 * (meanpt_state[imultibin] / meanpt[imultibin])
             << "%%)"
             << " + " << meanpt_syshe[imultibin] << " ("
             << 100 * (meanpt_syshe[imultibin] / meanpt[imultibin]) << "%%)"
             << " - " << meanpt_sysle[imultibin] << " ("
             << 100 * (meanpt_sysle[imultibin] / meanpt[imultibin]) << "%%)"
             << " +- " << averagpteerror_uncorr << " ("
             << 100 * (averagpteerror_uncorr / meanpt[imultibin]) << "%%)" << endl;
    }


    TGraphErrors* ge_stat = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &yield[0], &zeroerror[0], &yield_state[0]);
    TGraphErrors* ge_sys = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &yield[0], &dNdetaAxis_e[0], &yield_sysle[0]); // yield_sysle = yield_syshe
    TGraphErrors* ge_sys_cor = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &yield[0], &dNdetaAxis_e[0], &yield_sys_cor[0]); // yield_sysle = yield_syshe
    ge_stat->SetTitle("");
    ge_stat->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    ge_stat->GetYaxis()->SetTitle("d#it{N}_{#Xi^{*0}}/dy");
    ge_stat->GetXaxis()->SetLimits(0,45);
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
    ge_sys->GetYaxis()->SetTitle("d#it{N}_{#Xi^{*0}}/dy");
    ge_sys->GetXaxis()->SetLimits(0,45);
    //ge_sys->GetYaxis()->SetRangeUser(0,0.007);
    ge_sys->SetMinimum(1e-4);
    ge_sys->SetMaximum(4e-2);
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

    ge_sys_cor->Write("gYield_sys_cor");

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
        ptbin2, CorrectedYeild_syserr_7TeV);
    auto hSpectra_7TeV_staterr = MakeHistfromArray(
        "7TeV Spectra with statistical error", CorrectedYeild_7TeV2,
        ptbin2, CorrectedYeild_staterr_7TeV);
    /*
    AliPWGFunc * fm = new AliPWGFunc;
    fm->SetVarType(AliPWGFunc::VarType_t(AliPWGFunc::kdNdpt));
    TF1* func = fm->GetLevi (mass, 0.4, 750,3);
    func->SetParLimits(1,0.0001,20000);
    auto hfinal = YieldMean(hSpectra_7TeV_staterr, hSpectra_7TeV_syserr, func, 0.0, 10, 0.01, 0.1);
    cout << "7TeV Yield: " << hfinal->GetBinContent(1) << " +- " << hfinal->GetBinContent(2) << " + " << hfinal->GetBinContent(3) << " - " << hfinal->GetBinContent(4) << endl;
    */
    vector<double> x7 = {4.6};
    vector<double> x7e = {0.34};
    vector<double> x7z = {0.0};
    vector<double> x7l = {0.17};
    vector<double> x7h = {0.34};

    // HEP data
    vector<double> y7 = {2.56e-3};
    vector<double> y7e = {0.07e-3};
    vector<double> y7l = {0.37e-3};
    vector<double> y7h = {0.40e-3};
    vector<double> y7s_full = {0.385e-3}; // avg. of l,h

    vector<double> pt7 = {1.31};
    vector<double> pt7e = {0.02};
    vector<double> pt7l = {0.09};
    vector<double> pt7h = {0.09};
    vector<double> pt7s_full = {0.09};// avg. of l,h
    /*
    // This macro
    vector<double> y7new = {hfinal->GetBinContent(1)};
    vector<double> y7enew = {hfinal->GetBinContent(2)};
    vector<double> y7lnew = {hfinal->GetBinContent(4)};
    vector<double> y7hnew = {hfinal->GetBinContent(3)};    
    

    //this macro

    vector<double> pt7new = {hfinal->GetBinContent(5)};
    vector<double> pt7enew = {hfinal->GetBinContent(6)};
    vector<double> pt7lnew = {hfinal->GetBinContent(8)};
    vector<double> pt7hnew = {hfinal->GetBinContent(7)};    
    */
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
    /*
    TGraphErrors* ge_stat_7TeV = new TGraphErrors(1, &x7[0], &y7new[0], &x7e[0], &y7enew[0]);
    ge_stat_7TeV->SetLineColor(4);
    ge_stat_7TeV->SetMarkerColor(4);
    ge_stat_7TeV->SetMarkerStyle(20);
    ge_stat_7TeV->SetMarkerSize(0.5);
    ge_stat_7TeV->GetXaxis()->SetLimits(0,45);
    ge_stat_7TeV->SetFillColor(4);
    */
    //TGraphErrors* ge_stat_7TeV_HEP = new TGraphErrors(1, &x7[0], &y7[0], &x7e[0], &y7e[0]);
    TGraphAsymmErrors* ge_stat_7TeV_HEP = new TGraphAsymmErrors(
        1, &x7[0], &y7[0], &x7z[0], &x7z[0], &y7e[0], &y7e[0]);
    ge_stat_7TeV_HEP->SetLineColor(1);
    ge_stat_7TeV_HEP->SetMarkerColor(1);
    ge_stat_7TeV_HEP->SetMarkerStyle(20);
    ge_stat_7TeV_HEP->SetMarkerSize(0.5);
    ge_stat_7TeV_HEP->GetXaxis()->SetLimits(0,45);
    ge_stat_7TeV_HEP->SetFillColor(1);

    resultfile->cd();
    ge_stat_7TeV_HEP->Write("gYield7TeV_stat");

    /*
    TGraphAsymmErrors* ge_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &y7new[0], &x7e[0], &x7e[0], &y7lnew[0], &y7hnew[0]);
    //TGraphAsymmErrors* ge_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &y7[0], &x7e[0], &x7e[0], &y7l[0], &y7h[0]);
    ge_sys_7TeV->SetLineColor(4);
    ge_sys_7TeV->SetMarkerColor(4);
    ge_sys_7TeV->GetXaxis()->SetLimits(0,45);
    ge_sys_7TeV->SetFillColor(0);
    */

    //TGraphErrors* ge_sys_7TeV_HEP = new TGraphErrors(1, &x7[0], &y7[0], &x7e[0], &y7s_full[0]); // y7l = y7h
    TGraphAsymmErrors* ge_sys_7TeV_HEP = new TGraphAsymmErrors(
        1, &x7[0], &y7[0], &x7l[0], &x7h[0], &y7l[0], &y7h[0]);
    ge_sys_7TeV_HEP->SetLineColor(1);
    ge_sys_7TeV_HEP->SetMarkerColor(1);
    ge_sys_7TeV_HEP->SetFillColor(0);
    ge_sys_7TeV_HEP->GetXaxis()->SetLimits(0,45);

    resultfile->cd();
    ge_sys_7TeV_HEP->Write("gYield7TeV_syse");

    //TGraphErrors* ge_sys_cor_7TeV_HEP = new TGraphErrors(1, &x7[0], &y7[0], &x7e[0], &y7s_full[0]); // y7l = y7h
    TGraphAsymmErrors* ge_sys_cor_7TeV_HEP = new TGraphAsymmErrors(
        1, &x7[0], &y7[0], &x7l[0], &x7h[0], &y7l[0], &y7h[0]);
    ge_sys_7TeV_HEP->SetLineColor(1);
    ge_sys_7TeV_HEP->SetMarkerColor(1);
    ge_sys_7TeV_HEP->SetFillColor(0);
    ge_sys_7TeV_HEP->GetXaxis()->SetLimits(0,45);

    resultfile->cd();
    ge_sys_7TeV_HEP->Write("gYield7TeV_syse");

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
    ge_sys->Draw("a5");
    ge_stat->Draw("P");
    //ge_sys_7TeV_HEP->Draw("5");
    //ge_stat_7TeV_HEP->Draw("P");
    if(recalculate7TeV){
        //ge_sys_7TeV->Draw("5");
        //ge_stat_7TeV->Draw("P");
    }
    ge_sys->GetXaxis()->SetLimits(0,45);

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
    //if (recalculate7TeV)
        //legendyield->AddEntry(ge_sys_7TeV, "pp 7 TeV[Re-calculated]", "F");
    legendyield->AddEntry(ge_sys_7TeV_HEP, "pp 7 TeV[Paper]", "F");
    legendyield->Draw();
    SaveCanvas(cResults,Form("totalYields_%s",bininfo.Data()),"figs/");
    resultfile->cd();
    cResults->Write("cYield");

    //gSystem->Exit(1);
    // Mean pt
    TGraphErrors* gpt_stat = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &meanpt[0], &zeroerror[0], &meanpt_state[0]);
    TGraphErrors* gpt_sys = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &meanpt[0], &dNdetaAxis_e[0], &meanpt_syshe[0]); //meanpt_sysle = meanpt_syshe
    TGraphErrors* gpt_sys_cor = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &meanpt[0], &dNdetaAxis_e[0], &meanpt_sys_cor[0]);
    gpt_stat->SetTitle("");
    gpt_stat->GetXaxis()->SetTitle("< d#it{N}_{ch}/d#eta >");
    gpt_stat->GetYaxis()->SetTitle("< #it{p}_{T} > (GeV/c)");
    gpt_stat->GetXaxis()->SetRangeUser(0,30);
    //gpt_stat->GetYaxis()->SetRangeUser(1,1.7);
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
    gpt_sys->GetXaxis()->SetRangeUser(0,45);
    gpt_sys->SetMaximum(2.2);
    gpt_sys->SetMinimum(0.8);
    //gpt_sys->GetYaxis()->SetRangeUser(1,1.7);
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

    gpt_sys_cor->Write("gMeanpT_sys_cor");
    /*
    TGraphErrors* gpt_stat_7TeV = new TGraphErrors(1, &x7[0], &pt7new[0], &x7e[0], &pt7enew[0]);
    gpt_stat_7TeV->SetLineColor(4);
    gpt_stat_7TeV->SetMarkerColor(4);
    gpt_stat_7TeV->SetMarkerStyle(20);
    gpt_stat_7TeV->SetMarkerSize(0.5);
    gpt_stat_7TeV->SetFillColor(4);
    gpt_stat_7TeV->GetXaxis()->SetRangeUser(0,30);
    */
    TGraphAsymmErrors* gpt_stat_7TeV_HEP = new TGraphAsymmErrors(
        1, &x7[0], &pt7[0], &x7z[0], &x7z[0], &pt7e[0], &pt7e[0]);
    //TGraphErrors* gpt_stat_7TeV_HEP = new TGraphErrors(1, &x7[0], &pt7[0], &x7e[0], &pt7e[0]);
    gpt_stat_7TeV_HEP->SetLineColor(1);
    gpt_stat_7TeV_HEP->SetMarkerColor(1);
    gpt_stat_7TeV_HEP->SetMarkerStyle(20);
    gpt_stat_7TeV_HEP->SetMarkerSize(0.5);
    gpt_stat_7TeV_HEP->SetFillColor(1);
    gpt_stat_7TeV_HEP->GetXaxis()->SetRangeUser(0,30);

    resultfile->cd();
    gpt_stat_7TeV_HEP->Write("gMeanpT7TeV_stat");
    /*
    TGraphAsymmErrors* gpt_sys_7TeV = new TGraphAsymmErrors(1, &x7[0], &pt7new[0], &x7e[0], &x7e[0], &pt7lnew[0], &pt7hnew[0]);
    gpt_sys_7TeV->SetLineColor(4);
    gpt_sys_7TeV->SetMarkerColor(4);
    gpt_sys_7TeV->SetFillColor(0);
    gpt_sys_7TeV->GetXaxis()->SetRangeUser(0,30);
    */
    TGraphAsymmErrors* gpt_sys_7TeV_HEP = new TGraphAsymmErrors(
        1, &x7[0], &pt7[0], &x7l[0], &x7h[0], &pt7l[0], &pt7h[0]);
    //TGraphErrors* gpt_sys_7TeV_HEP = new TGraphErrors(1, &x7[0], &pt7[0], &x7e[0], &pt7s_full[0]); // pt7l = pt7h
    gpt_sys_7TeV_HEP->SetLineColor(1);
    gpt_sys_7TeV_HEP->SetMarkerColor(1);
    gpt_sys_7TeV_HEP->SetFillColor(0);
    gpt_sys_7TeV_HEP->GetXaxis()->SetRangeUser(0,30);

    resultfile->cd();
    gpt_sys_7TeV_HEP->Write("gMeanpT7TeV_syse");

    cResults->cd();
    cResults->SetLogy(false);
    gpt_sys->Draw("a5");
    gpt_stat->Draw("P");
    //gpt_sys_7TeV_HEP->Draw("5");
    //gpt_stat_7TeV_HEP->Draw("P");
    if(recalculate7TeV){
        //gpt_sys_7TeV->Draw("5");
        //gpt_stat_7TeV->Draw("P");
    }
    //gpt_sys->GetXaxis()->SetLimits(0,45);

    fStatusPad->Draw();

    legendyield->Draw();
    SaveCanvas(cResults,Form("totalMeanpT_%s",bininfo.Data()),"figs/");
    resultfile->cd();
    cResults->Write("cMeanpT");


    // Particle ratio    
    vector<double> RatioToXi;
    vector<double> RatioToXi_e;
    vector<double> RatioToXi_sys;
    vector<double> RatioToXi_sys_cor;

    vector<double> RatioToXi_7TeV;
    vector<double> RatioToXi_7TeV_e;
    vector<double> RatioToXi_7TeV_sys;
    
    vector<double> RatioToPi;
    vector<double> RatioToPi_e;
    vector<double> RatioToPi_sys;
    vector<double> RatioToPi_sys_cor;

    vector<double> RatioToPi_7TeV;
    vector<double> RatioToPi_7TeV_e;
    vector<double> RatioToPi_7TeV_sys;

    //7TeV 
    // pion https://arxiv.org/pdf/1504.00024.pdf
    double dndy_pi_7TeV = 2.245;
    double dndy_pi_7TeV_sys = 0.1; // stat.e negligible
    // Xi https://arxiv.org/abs/1204.0282
    double dndy_Xi_7TeV = 7.9e-3; // (8.0 + 7.8) /2
    double dndy_Xi_7TeV_e = 0.1e-3;
    double dndy_Xi_7TeV_sys = 0.6e-3; // +0.7, -0.5 -> avg. 0.6

    double tempRatioToPi_7TeV = y7[0]/dndy_pi_7TeV;
    RatioToPi_7TeV.push_back(tempRatioToPi_7TeV); 
    RatioToPi_7TeV_e.push_back( tempRatioToPi_7TeV*(y7e[0]/y7[0]));
    RatioToPi_7TeV_sys.push_back( tempRatioToPi_7TeV*sqrt( pow(0.385e-3/y7[0],2) + pow(dndy_pi_7TeV_sys/dndy_pi_7TeV,2)) );
    
    RatioToXi_7TeV.push_back(y7[0]/dndy_Xi_7TeV); 
    RatioToXi_7TeV_e.push_back( (y7[0]/dndy_Xi_7TeV)*sqrt( pow(y7e[0]/y7[0],2) + pow(dndy_Xi_7TeV_e/dndy_Xi_7TeV,2)));
    RatioToXi_7TeV_sys.push_back( (y7[0]/dndy_Xi_7TeV)*sqrt( pow(0.385e-3/y7[0],2) + pow(dndy_Xi_7TeV_sys/dndy_Xi_7TeV,2)) );


    for (int j = 0; j < ge_sys->GetN(); j++) {
        double tempRatioToXi = yield[j]/dNdeta_Xi[j];
        double tempRatioToXi_e = tempRatioToXi
                                 *sqrt( pow(yield_state[j]/yield[j],2)
                                 + pow(dNdeta_Xi_state[j]/dNdeta_Xi[j],2) );
        double tempRatioToXi_sys = tempRatioToXi
                                 *sqrt( pow(yield_sysle[j]/yield[j],2)
                                 + pow(dNdeta_Xi_sys[j]/dNdeta_Xi[j],2) );
        double tempRatioToXi_sys_cor = tempRatioToXi
                                 *sqrt( pow(yield_sys_cor[j]/yield[j],2)
                                 + pow(dNdeta_Xi_sys_uncor[j]/dNdeta_Xi[j],2) );
        
        double tempRatioToPi = yield[j]/dNdeta_pi[j];
        double tempRatioToPi_e = tempRatioToPi
                                 *sqrt( pow(yield_state[j]/yield[j],2)
                                 + pow(dNdeta_pi_e[j]/dNdeta_pi[j],2) );
        double tempRatioToPi_sys = tempRatioToPi
                                 *sqrt( pow(yield_sysle[j]/yield[j],2)
                                 + pow(dNdeta_pi_sys[j]/dNdeta_pi[j],2) );
        double tempRatioToPi_sys_cor = tempRatioToPi
                                 *sqrt( pow(yield_sys_cor[j]/yield[j],2)
                                 + pow(dNdeta_pi_sys[j]/dNdeta_pi[j],2) );
        
        RatioToXi.push_back(tempRatioToXi);
        RatioToXi_e.push_back(tempRatioToXi_e);
        RatioToXi_sys.push_back(tempRatioToXi_sys);
        RatioToXi_sys_cor.push_back(tempRatioToXi_sys_cor);

        RatioToPi.push_back(tempRatioToPi);
        RatioToPi_e.push_back(tempRatioToPi_e);
        RatioToPi_sys.push_back(tempRatioToPi_sys);
        RatioToPi_sys_cor.push_back(tempRatioToPi_sys_cor);

        cout << "Xi" << endl;
        cout << "$"<< Form("%2.f", multibin[j][0]) << " - "
             << Form("%2.f", multibin[j][1]) << " ("
             << Form("%.2f",dNdetaAxis[j]) << ")$ & $"
             << Form("%.2f",yield[j]*1000)
             << "$ & $" << Form("%.3f",dNdeta_Xi[j]*1000)
             << "$ & $" << Form("%.3f",tempRatioToXi)
             << " \\pm " << Form("%.3f",tempRatioToXi_e)
             << " \\pm " << Form("%.3f",tempRatioToXi_sys)
             << " \\pm " << Form("%.3f",tempRatioToXi_sys_cor)
             << "$ \\\\" << endl;
        cout << "Pi" << endl;
        cout << "$"<< Form("%2.f", multibin[j][0]) << " - "
             << Form("%2.f", multibin[j][1]) << " ("
             << Form("%.2f",dNdetaAxis[j]) << ")$ & $"
             << Form("%.2f",yield[j]*1000)
             << "$ & $" << Form("%.2f",dNdeta_pi[j])
             << "$ & $" << Form("%.2f",tempRatioToPi*1000)
             << " \\pm " << Form("%.2f",tempRatioToPi_sys*1000)
             << " \\pm " << Form("%.2f",tempRatioToPi_sys_cor*1000)
             << "$ \\\\" << endl;

    }

    TGraphErrors* grRatioXistat = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &RatioToXi[0], &zeroerror[0], &RatioToXi_e[0]);
    TGraphErrors* grRatioXisys = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &RatioToXi[0], &dNdetaAxis_e[0], &RatioToXi_sys[0]);
    TGraphErrors* grRatioXisys_uncor = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &RatioToXi[0], &dNdetaAxis_e[0], &RatioToXi_sys_cor[0]);

    TGraphAsymmErrors* grRatioXistat_7TeV =
        new TGraphAsymmErrors(1, &x7[0], &RatioToXi_7TeV[0], &x7z[0], &x7z[0],
                              &RatioToXi_7TeV_e[0], &RatioToXi_7TeV_e[0]);
    TGraphAsymmErrors* grRatioXisys_7TeV =
        new TGraphAsymmErrors(1, &x7[0], &RatioToXi_7TeV[0], &x7l[0], &x7h[0],
                              &RatioToXi_7TeV_sys[0], &RatioToXi_7TeV_sys[0]);
    //TGraphErrors* grRatioXistat_7TeV = new TGraphErrors(multibin.size(), &x7[0], &RatioToXi_7TeV[0], &x7e[0], &RatioToXi_7TeV_e[0]);
    //TGraphErrors* grRatioXisys_7TeV = new TGraphErrors(multibin.size(), &x7[0], &RatioToXi_7TeV[0], &x7e[0], &RatioToXi_7TeV_sys[0]);

    grRatioXisys_uncor->SetMarkerColor(kOrange + 10);
    grRatioXisys_uncor->SetFillColorAlpha(kOrange + 10, 0.3);
    //grRatioXisys_uncor->SetFillStyle(3001);
    grRatioXisys_uncor->SetLineColor(0);
    
    grRatioXisys->SetMaximum(0.5);
    grRatioXisys->SetMinimum(0.1);
    grRatioXisys->SetFillColorAlpha(kOrange + 10, 0.0);
    grRatioXisys->SetMarkerColor(kOrange + 10);
    grRatioXisys->SetLineColor(kOrange + 10);

    grRatioXistat->SetMarkerColor(kOrange + 10);
    grRatioXistat->SetLineColor(kOrange + 10);

    grRatioXisys_7TeV->SetMarkerColor(kBlack);
    grRatioXisys_7TeV->SetLineColor(kBlack);

    grRatioXistat_7TeV->SetMarkerColor(kBlack);
    grRatioXistat_7TeV->SetLineColor(kBlack);

    cResults->cd();

    grRatioXisys->GetXaxis()->SetLimits(0,45);

    grRatioXisys->Draw("a5");
    grRatioXisys_uncor->Draw("5");
    grRatioXistat->Draw("P");
    //grRatioXisys_7TeV->Draw("5");
    //grRatioXistat_7TeV->Draw("P");

    SaveCanvas(cResults,"RatioToXi","figs/");
    resultfile->cd();
    grRatioXistat->Write("gRatioToXi_stat");
    grRatioXisys->Write("gRatioToXi_sys");
    grRatioXisys_uncor->Write("gRatioToXi_sys_cor");
    grRatioXistat_7TeV->Write("gRatioToXi_7TeV_stat");
    grRatioXisys_7TeV->Write("gRatioToXi_7TeV_sys");


    TGraphErrors* grRatioPistat = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &RatioToPi[0], &zeroerror[0], &RatioToPi_e[0]);
    TGraphErrors* grRatioPisys = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &RatioToPi[0], &dNdetaAxis_e[0], &RatioToPi_sys_cor[0]);
    TGraphErrors* grRatioPisys_full = new TGraphErrors(multibin.size(), &dNdetaAxis[0], &RatioToPi[0], &dNdetaAxis_e[0], &RatioToPi_sys[0]);

    TGraphErrors* grRatioPistat_7TeV = new TGraphErrors(multibin.size(), &x7[0], &RatioToPi_7TeV[0], &x7e[0], &RatioToPi_7TeV_e[0]);
    TGraphErrors* grRatioPisys_7TeV = new TGraphErrors(multibin.size(), &x7[0], &RatioToPi_7TeV[0], &x7e[0], &RatioToPi_7TeV_sys[0]);
    
    //grRatioPisys->SetMaximum(2e-3);
    grRatioPisys->GetXaxis()->SetLimits(0,45);
    grRatioPisys->SetMaximum(2e-3);
    grRatioPisys->SetMinimum(5e-4);
    //grRatioPisys->GetXaxis()->SetLimits(1,2e3);
    //grRatioPisys->SetFillColorAlpha(kBlack, 0.0);
    grRatioPisys->SetMarkerColor(kOrange + 10);
    grRatioPisys->SetLineColor(kOrange + 10);

    grRatioPistat->SetMarkerColor(kOrange + 10);
    grRatioPistat->SetLineColor(kOrange + 10);

    grRatioPisys_7TeV->SetMarkerColor(kBlack);
    grRatioPisys_7TeV->SetLineColor(kBlack);

    grRatioPistat_7TeV->SetMarkerColor(kBlack);
    grRatioPistat_7TeV->SetLineColor(kBlack);


    cResults->cd();
    //cResults->SetLogx();
    //cResults->SetLogy();
    grRatioPisys->Draw("a5");
    grRatioPistat->Draw("P");
    //grRatioPisys_7TeV->Draw("5");
    //grRatioPistat_7TeV->Draw("P");

    SaveCanvas(cResults,"RatioToPi","figs/");
    resultfile->cd();
    grRatioPistat->Write("gRatioPi_stat");
    grRatioPisys->Write("gRatioPi_sys");
    grRatioPisys_full->Write("gRatioPi_sys_full");
    grRatioPistat_7TeV->Write("gRatioPi_7TeV_stat");
    grRatioPisys_7TeV->Write("gRatioPi_7TeV_sys");

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
    finalfile = GetFinalFileName(multi_start,multi_end);
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_sys", multi_start, multi_end));
    gROOT->cd();
    TH1D* hReturn = (TH1D*)hr->Clone();
    inputfile->Close();
    return hReturn;
}
TH1D* GetSpectrasysNocor(double multi_start, double multi_end) {
    finalfile = GetFinalFileName(multi_start,multi_end);
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(
        Form("hSpectra_%.2f_%.2f_sys_noCorrelation", multi_start, multi_end));
    gROOT->cd();
    TH1D* hReturn = (TH1D*)hr->Clone();
    inputfile->Close();
    return hReturn;
}
TH1D* GetSpectrastat(double multi_start, double multi_end){
    finalfile = GetFinalFileName(multi_start,multi_end);
    TFile* inputfile = new TFile(finalfile.Data());
    TH1D* hr = (TH1D*)inputfile->Get(Form("hSpectra_%.2f_%.2f_stat",
                               multi_start,
                               multi_end));

    gROOT->cd();
    TH1D* hReturn = (TH1D*)hr->Clone();
    inputfile->Close();
    return hReturn;
}
vector<double> GetYieldError(double multi_start, double multi_end){
    vector<double> temp;
    bool isINEL = false;
    TString bininfo;
    if ((multi_start == 0) && (multi_end == 0)){ 
        //INEL case
        bininfo = "INEL";
        isINEL = true;
    } 
    else if(multi_end < 0.5)
        bininfo = Form("%.2f-%.2f",multi_start,multi_end);
    else
        bininfo = Form("%.0f-%.0f",multi_start,multi_end);
    TFile* inputfile = new TFile(Form("AnalysisResults_Xi1530_YieldMean_%s.root",bininfo.Data()));

    TString histname;
    if (isINEL)
        histname = "hYield_INEL";
    else
        histname = Form("hYield_%.2f-%.2f", multi_start, multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(histname);
    for(int i = 0; i < hr->GetNbinsX(); i++)
        temp.push_back(hr->GetBinContent(i+1));
    inputfile->Close();
    return temp;
}
vector<double> GetMeanPtError(double multi_start, double multi_end){
    vector<double> temp;
    bool isINEL = false;
    TString bininfo;
    if ((multi_start == 0) && (multi_end == 0)){ 
        //INEL case
        bininfo = "INEL";
        isINEL = true;
    } 
    else if(multi_end < 0.5)
        bininfo = Form("%.2f-%.2f",multi_start,multi_end);
    else
        bininfo = Form("%.0f-%.0f",multi_start,multi_end);
    TFile* inputfile = new TFile(Form("AnalysisResults_Xi1530_YieldMean_%s.root",bininfo.Data()));

    TString histname;
    if (isINEL)
        histname = "hMeanPt_INEL";
    else
        histname = Form("hMeanPt_%.2f-%.2f", multi_start, multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(histname);
    for(int i = 0; i < hr->GetNbinsX(); i++)
        temp.push_back(hr->GetBinContent(i+1));
    inputfile->Close();
    return temp;
}
vector<double> GetPidNdetawithError(double multi_start, double multi_end, int errortype){
    // Return pion's dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}
    // errortype:
    //   1: stat err.
    //   2: total sys err.
    TFile* fPi =
        TFile::Open("OutputYields_pp13mult.root","READ");
    TGraphAsymmErrors* pisys;
    if(errortype == 2)
        pisys = (TGraphAsymmErrors*)fPi->Get("pi_Syst");
    else
        pisys = (TGraphAsymmErrors*)fPi->Get("pi_Stat");
    
    vector<double> returnarray;
    int lenpisys = pisys->GetN();

    //--dNdeta histogram of pion+ + pion-
    // Ref: given by Anders, need to update ref.
    // Error was stat+sys error, so choosed bigger error(sym) error.

    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta = {0};
    vector<double> dNchdeta_e = {0};
    for(int i = 0; i < lenpisys+1; i++){
        dNchdeta.push_back(pisys->GetY()[i]);
        dNchdeta_e.push_back(pisys->GetErrorYhigh(i));
        
        cout << "bin :" << i << " " << dNchdeta[i] << " +- " << dNchdeta_e[i] << endl;
    }
    
    /*
    vector<double> dNchdeta = 
    {0, 24.605455025, 19.016207266, 15.477779961, 13.241832514, 11.612550017, 9.742647922, 7.779575692,  6.241633459, 4.530678113, 2.713659699};
    vector<double> dNchdeta_e;
    if(errortype == 2)
        dNchdeta_e = {0,  1.121689195,  0.856354796,  0.685848455,  0.582627504,   0.498773083,  0.415988997, 0.327417792, 0.261067034, 0.186885663, 0.110000678};
    else
        dNchdeta_e = {0,  0.009062988,  0.004189384,  0.003506052,  0.003237967,   0.003090847,  0.002062999, 0.001795797, 0.001551303, 0.000883705, 0.000526020};
    */
    // input must be in the multiplicity range
    if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_start) == end(dNchdeta_multibin))
        return {9999,9999};
    if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_end) == end(dNchdeta_multibin))
        return {9999,9999};

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
    returnarray.push_back(result/2);

    // Error
    double error = 0.;
    for(int i = 1; i < gap+1; i++)
        error += pow( dNchdeta_e[i+left], 2); 
        
    error = sqrt(error);
    returnarray.push_back(error/2);
    
    cout << "Pi yield" << endl;
    cout << "$"<< multi_start << " - " << multi_end << "$ & $" << result << " \\pm " << error/2 << "$ \\\\" << endl;

    return returnarray;
}
vector<double> GetXidNdetawithError(double multi_start, double multi_end, int errortype){
    // Return 13 TeV Xi's dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}
    // errortype:
    //   1: stat err.
    //   2: total sys err.
    //   3: un correlated sys err.

    TFile* fXi13TeV = TFile::Open("YieldMean-Xi-V0M.root", "READ");
    TGraphErrors* grXistat = (TGraphErrors*)fXi13TeV->Get("grYieldXiStat");
    TGraphErrors* grXisys = (TGraphErrors*)fXi13TeV->Get("grYieldXiSystTot");
    TGraphErrors* grXisys_uncor = (TGraphErrors*)fXi13TeV->Get("grYieldXiSystUncorr");

    vector<double> returnarray;

    //--dNdeta histogram of pion+ + pion-
    // Ref: given by Fiorella Fionda, https://alice-notes.web.cern.ch/node/478

    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta;
    vector<double> dNchdeta_e;

    dNchdeta.push_back(0);
    dNchdeta_e.push_back(0);
    for (int j = 0; j < grXistat->GetN(); j++) {
        dNchdeta.push_back(grXistat->GetY()[j]);
        if(errortype == 3)
            dNchdeta_e.push_back(grXisys_uncor->GetErrorY(j));
        else if(errortype == 2)
            dNchdeta_e.push_back(grXisys->GetErrorY(j));
        else
            dNchdeta_e.push_back(grXistat->GetErrorY(j));
    }
    /*
    // input must be in the multiplicity range
    if((multi_start == 0) && (multi_end == 0.01)){
        if(errortype == 3)
            return {0.166411/2, 0.00861725/2};
        if(errortype == 2)
            return {0.166411/2, 0.0136379/2};
        return {0.166411/2, 0.00273244/2};
    }
    else if((multi_start == 0.01) && (multi_end == 0.05)){
        if(errortype == 3)
            return {0.1485/2, 0.00768977/2};
        if(errortype == 2)
            return {0.1485/2, 0.0121701/2};
        return {0.1485/2, 0.00243834/2};
    } else if ((multi_start == 0.05) && (multi_end == 0.1)) {
        if(errortype == 3)
            return {0.138608/2, 0.00717754/2};
        if(errortype == 2)
            return {0.138608/2, 0.0113594/2};
        return {0.138608/2, 0.00227592/2};
    }
    */
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
    returnarray.push_back(result/2);

    // Error
    double error = 0.;
    for(int i = 1; i < gap+1; i++)
        error += pow( dNchdeta_e[i+left], 2); 
        
    error = sqrt(error);
    returnarray.push_back(error/2);
    
    cout << "Xi yield" << endl;
    cout << "$"<< multi_start << " - " << multi_end << "$ & $" << result << " \\pm " << error/2 << "$ \\\\" << endl;
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
TString GetFinalFileName(double multi_start, double multi_end){
    TString finalfile_text;
    TString bininfo = "_";
    if ((multi_start == 0) && (multi_end == 0)){ 
        //INEL case
        bininfo += "INEL";
    }
    else if (multi_end < 0.5) {
        bininfo += Form("%.2f", multi_start);
        bininfo += "-";
        bininfo += Form("%.2f", multi_end);
    }
    else {
        bininfo += Form("%.0f", multi_start);
        bininfo += "-";
        bininfo += Form("%.0f", multi_end);   
    }
    finalfile_text = Form("./AnalysisResults_Xi1530_systematic%s.root", bininfo.Data());

    return finalfile_text;
}