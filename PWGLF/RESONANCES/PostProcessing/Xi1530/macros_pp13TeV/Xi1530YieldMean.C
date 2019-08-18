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
vector<double> GetdNdetawithError(double multi_start, double multi_end);
vector<double> GetPidNdetawithError(double multi_start, double multi_end);
TString finalfile = "AnalysisResults_Xi1530_systematic001030507000.root";
TString workdirectory = "/Users/blim/alidock/Postprocessing/data/";
enum {kFitExpPt=1, kFitLevi, fFitExpMt, kFitBoltzmann, kFitBlastWave, kFitBoseEinstein, kFitFermiDirac};
vector<TString> functions = {"", "kFitExpPt", "kFitLevi", "fFitExpMt", "kFitBoltzmann", "kFitBlastWave", "kFitBoseEinstein", "kFitFermiDirac"};

vector<vector<vector<double>>> // { {first pT bin variation}, {last pTbin variation} } // for each function
            functionFitRangeVar = { { {}, {} },            // null
                                    { {0.8}, {1.6, 2.0, 2.4} }, // ExpPt, only front part.
                                    { {0.8}, {2.4,3.2,4.0,4.8,5.6,8.8} }, // Levy, all region.
                                    { {0.8}, {1.6, 2.0, 2.4} }, // ExpMt, only front part.
                                    { {0.8}, {1.6, 2.0, 2.4} }, // Boltzmann, only front part.
                                    { {0.8}, {2.4, 3.2,4.0,4.8,5.6,8.8} }, // BGBW, all region.
                                    { {0.8}, {1.6, 2.0} }, // BoseEinstein, only front part.
                                    { {0.8}, {1.6, 2.0, 2.4} } // FermiDirac, only front part.
                                    };
Int_t maxtrial = 10000;
double mass = 1.5318;
TString fitoption = "0q"; //default "0q"

void Xi1530YieldMean(double multis = 0, double multie = 100){

    // Reweighting efficiency
    TFile* fReweighting =
        TFile::Open("AnalysisResults_Xi1530_efficiencyReweighing.root","READ");
    vector<TH1D*> hCorrectionfactor_rewight;
    hCorrectionfactor_rewight.push_back((TH1D*)fReweighting->Get("_correction_i1")); 
    hCorrectionfactor_rewight.push_back((TH1D*)fReweighting->Get("_correction_i2")); 
    hCorrectionfactor_rewight.push_back((TH1D*)fReweighting->Get("_correction_i3")); 

    TCanvas* cResults = new TCanvas("cResults", "cResults", 960, 720);
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    cResults->SetTickx();
    cResults->SetTicky();
    cResults->SetTopMargin(0.05);
    cResults->SetLeftMargin(0.10);
    //cSigbkg->SetBottomMargin(0.01);
    cResults->SetRightMargin(0.01);
    cResults->SetFillStyle(0);
    cResults->Draw();
    cResults->SetLogy();

    vector<vector<double>> multibin = {
        // {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
        // {0, 10}, {10, 30}, {30, 50}, {50, 100}};
        // {0, 100}, {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}, {30, 100}};
        // {0, 100}};
         {multis, multie}};
    TString bininfo = "";
    for (auto const& bin : multibin)
        bininfo += Form("%.0f", bin[0]);
    bininfo += Form("%.0f", multibin.back()[1]);

    bool recalculate7TeV = false;
    TString filename;
    if(recalculate7TeV)
        filename = "AnalysisResults_Xi1530_YieldMean_7TeV.root";
    else
        filename = Form("AnalysisResults_Xi1530_YieldMean_%s.root", bininfo.Data());

    TFile* resultfile = TFile::Open(filename.Data(),"RECREATE");


    vector<TH1D*> hspectra_sys;
    vector<TH1D*> hspectra_sys_nocor;
    vector<TH1D*> hspectra_stat;
    if(recalculate7TeV){
        // 7TeV Results
        // Old data
        // from HEPDATA https://www.hepdata.net/record/ins1300380
        vector<double> ptbin7TeV = {0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6};
        double* ptbin7TeV_array = &ptbin7TeV[0];
        vector<double> CorrectedYeild_7TeV = {0.00138,  0.00109, 0.00078,
                                            0.00046, 0.000226, 6.6e-05, 2.34e-05,
                                            8.1e-06};
        vector<double> CorrectedYeild_syserr_7TeV = {
            0.00017,  6.0e-05,  5.0e-05, 2.2e-05, 1.15e-05,
            3.73e-06, 1.96e-06, 4.54e-07};
        vector<double> CorrectedYeild_staterr_7TeV = {
            4.0e-05, 3.0e-05, 2.0e-05, 1.5e-05, 3.38e-06,
            2.08e-6, 8.11e-7, 9.09e-7};
        vector<double> CorrectedYeild_total_syserr_7TeV; // add sys+un.cor.sys, 
        // value was +9, +6 -> try use 7.5%(average)
        for (int value = 0; value < CorrectedYeild_syserr_7TeV.size(); value++) {
            double tempe = 0.0;
            tempe += pow(CorrectedYeild_syserr_7TeV[value]/CorrectedYeild_7TeV[value], 2);
            tempe += pow(0.075, 2); // 7.5%
            CorrectedYeild_total_syserr_7TeV.push_back(CorrectedYeild_7TeV[value]*sqrt(tempe));
        }

        

        TH1D* hRatio7TeV_stat = new TH1D("hRatio7TeV_stat", "hRatio7TeV_stat", ptbin7TeV.size() - 1, ptbin7TeV_array);
        TH1D* hRatio7TeV_sys = new TH1D("hRatio7TeV_sys", "hRatio7TeV_sys", ptbin7TeV.size() - 1, ptbin7TeV_array);
        TH1D* hRatio7TeV_sys_full = new TH1D("hRatio7TeV_sys_full", "hRatio7TeV_sys_full", ptbin7TeV.size() - 1, ptbin7TeV_array);
        for (int i = 0; i < CorrectedYeild_7TeV.size(); i++){
            hRatio7TeV_stat->SetBinContent(i + 1, CorrectedYeild_7TeV[i]);
            hRatio7TeV_sys->SetBinContent(i + 1, CorrectedYeild_7TeV[i]);
            hRatio7TeV_sys_full->SetBinContent(i + 1, CorrectedYeild_7TeV[i]);

            hRatio7TeV_stat->SetBinError(i + 1, CorrectedYeild_staterr_7TeV[i]);
            hRatio7TeV_sys->SetBinError(i + 1, CorrectedYeild_syserr_7TeV[i]);
            hRatio7TeV_sys_full->SetBinError(i + 1, CorrectedYeild_total_syserr_7TeV[i]);
        }

        hspectra_stat.push_back(hRatio7TeV_stat);
        hspectra_sys_nocor.push_back(hRatio7TeV_sys);
        hspectra_sys.push_back(hRatio7TeV_sys_full);

        resultfile->cd();
        hRatio7TeV_sys->Write("0.00-100.00_SYS_corrected_7TeV");
        hRatio7TeV_sys_full->Write("0.00-100.00_SYS_full_corrected_7TeV");
        hRatio7TeV_stat->Write("0.00-100.00_stat_corrected_7TeV");
    }
    else{
        hspectra_sys = GetSysSpectra(multibin);
        hspectra_sys_nocor = GetSysSpectraNocor(multibin);
        hspectra_stat  = GetStatSpectra(multibin);
        for (int multbin = 0; multbin < hspectra_sys.size(); multbin++){
            for (int bin = 0; bin < hspectra_sys[0]->GetNbinsX(); bin++) {
                double weightcorrectionfactor = hCorrectionfactor_rewight.back()->GetBinContent(bin+1);
                hspectra_sys[multbin]->SetBinContent(bin+1, hspectra_sys[multbin]->GetBinContent(bin+1)/weightcorrectionfactor);
                hspectra_sys_nocor[multbin]->SetBinContent(bin+1, hspectra_sys_nocor[multbin]->GetBinContent(bin+1)/weightcorrectionfactor);
                hspectra_stat[multbin]->SetBinContent(bin+1, hspectra_stat[multbin]->GetBinContent(bin+1)/weightcorrectionfactor);

                hspectra_sys[multbin]->GetXaxis()->SetTitle("#it{P}_{T} (GeV/#it{c})");
                hspectra_sys[multbin]->GetYaxis()->SetTitle("1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
                hspectra_sys_nocor[multbin]->GetXaxis()->SetTitle("#it{P}_{T} (GeV/#it{c})");
                hspectra_sys_nocor[multbin]->GetYaxis()->SetTitle("1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
                hspectra_stat[multbin]->GetXaxis()->SetTitle("#it{P}_{T} (GeV/#it{c})");
                hspectra_stat[multbin]->GetYaxis()->SetTitle("1/N_{event}d^{2}N/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");
            }
            resultfile->cd();
            hspectra_sys[multbin]->Write(
                    Form("%.2f-%.2f_SYS_corrected",
                        multibin[multbin][0],
                        multibin[multbin][1]
                        ));
            hspectra_stat[multbin]->Write(
                    Form("%.2f-%.2f_stat_corrected",
                        multibin[multbin][0],
                        multibin[multbin][1]
                        ));
            hspectra_sys_nocor[multbin]->Write(
                    Form("%.2f-%.2f_SYS_uncor_corrected",
                        multibin[multbin][0],
                        multibin[multbin][1]
                        ));
        }
    }

    // ch <dn/deta>
    vector<double> dNdetaAxis;
    vector<double> dNdetaAxis_e;
    vector<double> dNdetaAxis_e_half;
    vector<double> zeroerror;

    //pion <dn/dy>
    vector<double> dNdeta_pi;
    vector<double> dNdeta_pi_e;
    
    vector<vector<double>> yield(multibin.size());
    vector<double> yield_state;

    vector<vector<double>> meanpt(multibin.size());
    vector<double> meanpt_state;
    vector<double> extraYields;

    AliPWGFunc * fm = new AliPWGFunc;
    fm->SetVarType(AliPWGFunc::VarType_t(AliPWGFunc::kdNdpt));
    TF1* func = 0;
    int slopePar=0;
    vector<vector<TF1*>> funcs(multibin.size());
    vector<TString> fitnames;
    //vector<double> fitrange = {0.8,8.8};
    // Stat error and central value?
    TString nameinput;
    int flinestyle = 1;
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        for (int fitvar = 1; fitvar < functions.size(); fitvar++){
            //if(fitvar > 2) continue;
            flinestyle = 1;
            for (auto const& fitstart : functionFitRangeVar[fitvar][0]) { // fit bin start
                for (auto const& fitend : functionFitRangeVar[fitvar][1]) { // fit bin end
                    nameinput = Form("%s %.2f - %.2f GeV/#it{c}", functions[fitvar].Data(), fitstart, fitend);
                    fitnames.push_back(nameinput);
                    if (fitvar == kFitLevi)   {
                        func = fm->GetLevi (mass, 0.4, 750,3);
                        func->SetParLimits(1,0.0001,20000);
                        slopePar = 2;
                        //fitrange = {0.8,8.8};
                    }
                    if (fitvar == kFitExpPt)  {
                        func = fm->GetPTExp(0.2, 20);
                        slopePar = 1;
                        //fitrange = {0.8,1.6};
                    }
                    if (fitvar == fFitExpMt)  {
                        func = fm->GetMTExp(mass,0.2, 20);
                        slopePar = 1;
                    }
                    if (fitvar == kFitBoltzmann)  {
                        func = fm->GetBoltzmann(mass, 0.2, 20);
                        slopePar = 1;
                        //fitrange = {0.8,1.6};
                    }
                    if (fitvar == kFitBlastWave)  {
                        func = fm->GetBGBW(mass,0.6,0.3, 1, 1e5);// beta, T, n, norm 
                        func->SetParLimits(1, 0.1, 0.99);
                        func->SetParLimits(2, 0.01, 1);
                        func->SetParLimits(3, 0.01, 2);
                        slopePar = 2;
                        //fitrange = {0.8,8.8};
                    }
                    if (fitvar == kFitBoseEinstein)  {
                        func = fm->GetBoseEinstein(mass,0.3,20);
                        slopePar = 1;
                        //fitrange = {0.8,1.6};
                    }
                    if (fitvar == kFitFermiDirac)  {
                        func = fm->GetFermiDirac(mass,0.3,20);
                        slopePar = 1;
                        //fitrange = {0.8,1.6};
                    }
                    auto hfinal =
                            YieldMean(hspectra_stat[imultibin], hspectra_sys_nocor[imultibin],
                                    func, 0.05, 10, 0.01, 0.1, fitoption.Data(), Form("log_%s.root",bininfo.Data()),fitstart, fitend);
                    yield[imultibin].push_back(hfinal->GetBinContent(1));
                    meanpt[imultibin].push_back(hfinal->GetBinContent(5));
                    func->SetLineWidth(1);
                    func->SetLineStyle(flinestyle);
                    func->SetRange(0,fitend);
                    resultfile->cd();
                    func->Write(
                        Form("%.2f-%.2f_%s_%.2f-%.2f",
                        multibin[imultibin][0],
                        multibin[imultibin][1],
                        functions[fitvar].Data(),
                        fitstart,fitend
                        ));
                    funcs[imultibin].push_back(func);
                    flinestyle++;
                }
            }
        }
    }

    int sysColorPallet = GetSerialColors(functions.size()-1);
    // for memo, small
    TLatex* tm = new TLatex();
    tm->SetNDC();
    tm->SetTextSize(0.03);
    auto legendhead = new TLegend(.65, .65, .95, .9);
    legendhead->SetNColumns(2);
    legendhead->SetBorderSize(0);
    legendhead->SetFillStyle(0);
    vector<double> fitsyserrory;
    vector<double> fitsyserrorm;
    //QA
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        cout << "Mutli: " << multibin[imultibin][0] << " - " << multibin[imultibin][0] << endl;
        legendhead->Clear();
        cResults->cd();
        hspectra_sys_nocor[imultibin]->GetXaxis()->SetRangeUser(0, 3.2);
        hspectra_sys_nocor[imultibin]->SetMaximum(
            5 * hspectra_sys_nocor[imultibin]->GetBinContent(2));
        hspectra_sys_nocor[imultibin]->SetMinimum(
            0.5 * hspectra_sys_nocor[imultibin]->GetBinContent(6));
        hspectra_sys_nocor[imultibin]->Draw("E");
        legendhead->AddEntry(hspectra_sys_nocor[imultibin], "Data", "PLE");

        // Result?
        auto ymax = max_element(std::begin(yield[imultibin]), std::end(yield[imultibin]));
        auto ymin = min_element(std::begin(yield[imultibin]), std::end(yield[imultibin]));
        auto mmax = max_element(std::begin(meanpt[imultibin]), std::end(meanpt[imultibin]));
        auto mmin = min_element(std::begin(meanpt[imultibin]), std::end(meanpt[imultibin]));

        TH1* systematicyields =
            new TH1D("systematic_yields", "", 100, *ymin*0.75 , *ymax*1.25);
        systematicyields->GetXaxis()->SetTitle("<d#it{N}/d#it{y}>");
        systematicyields->GetYaxis()->SetTitle("Counts");
        TH1* systematicmeanpTs =
            new TH1D("systematic_meanpTs", "", 100, *mmin*0.75 , *mmax*1.25);
        systematicmeanpTs->GetXaxis()->SetTitle("<#it{p}_{T}>");
        systematicmeanpTs->GetYaxis()->SetTitle("Counts");

        double centralyield;
        double centralmeanpT;
        int loopcheck = 0;
        //for (int fitvar = 0; fitvar < yield[0].size(); fitvar++){ // old loop, loop for the contents
        for (int fitvar = 1; fitvar < functions.size(); fitvar++){
            cout << "fit funcion: " << functions[fitvar].Data() << endl;
            for (auto const& fitstart : functionFitRangeVar[fitvar][0]) { // fit bin start
                for (auto const& fitend : functionFitRangeVar[fitvar][1]) { // fit bin end

                    funcs[imultibin][loopcheck]->SetLineColor(sysColorPallet + functions.size() - fitvar);
                    funcs[imultibin][loopcheck]->Draw("same");
                    legendhead->AddEntry(funcs[imultibin][loopcheck], fitnames[loopcheck], "L");
                    
                    cout << "yields: " << yield[imultibin][loopcheck] << " meanpT" << meanpt[imultibin][loopcheck]  <<  ", Chi^2: " << funcs[imultibin][loopcheck]->GetChisquare() << ", NDF: " << funcs[imultibin][loopcheck]->GetNDF() << ", Chi^2/NDF: " << funcs[imultibin][loopcheck]->GetChisquare()/funcs[imultibin][loopcheck]->GetNDF() << endl;

                    if( (fitvar == kFitLevi)
                        && (fitstart == 0.8)
                        && (fitend == 8.8)
                    ) {
                        centralyield = yield[imultibin][loopcheck];
                        centralmeanpT = meanpt[imultibin][loopcheck];
                    }

                    if(fitvar != kFitBoseEinstein) systematicyields->Fill(yield[imultibin][loopcheck]);
                    if(fitvar != kFitBoseEinstein) systematicmeanpTs->Fill(meanpt[imultibin][loopcheck]);
                    loopcheck++;

                }
            }
        }
        legendhead->Draw();
        SaveCanvas(cResults,
                   Form("Spectrafit_%.2f-%.2f_zoom", multibin[imultibin][0],
                    multibin[imultibin][1]),"figs/LowpTExtrapolate/");
        resultfile->cd();
        cResults->Write(
                Form("Spectrafit_%.2f-%.2f_zoom",
                    multibin[imultibin][0],
                    multibin[imultibin][1]));
        systematicyields->Write(
                Form("Yieldsys_%.2f-%.2f",
                    multibin[imultibin][0],
                    multibin[imultibin][1]));
        systematicmeanpTs->Write(
                Form("MeanPtsys_%.2f-%.2f",
                    multibin[imultibin][0],
                    multibin[imultibin][1]));
        cResults->SetLogy(false);
        systematicyields->Draw();
        TLine* centralarrow1 = new TLine(centralyield, 0, centralyield, 
                                        systematicyields->GetMaximum());
        centralarrow1->SetLineColor(kRed);
        centralarrow1->SetLineStyle(2);
        centralarrow1->Draw();
        SaveCanvas(cResults,
                   Form("Yields_systematic_%.2f-%.2f", multibin[imultibin][0],
                    multibin[imultibin][1]),"figs/LowpTExtrapolate/");
        cout << "mean: " << systematicyields->GetMean() << ", stdev: " << systematicyields->GetRMS() 
             << ", percentage: " << systematicyields->GetRMS()/systematicyields->GetMean()  << endl;
        systematicmeanpTs->Draw();
        TLine* centralarrow2 = new TLine(centralmeanpT, 0, centralmeanpT, 
                                        systematicyields->GetMaximum());
        centralarrow2->SetLineColor(kRed);
        centralarrow2->SetLineStyle(2);
        centralarrow2->Draw();
        SaveCanvas(cResults,
                   Form("MeanPts_systematic_%.2f-%.2f", multibin[imultibin][0],
                    multibin[imultibin][1]),"figs/LowpTExtrapolate/");
        fitsyserrory.push_back(systematicyields->GetRMS());
        fitsyserrorm.push_back(systematicmeanpTs->GetRMS());
    }

    // Systematic error
    func = fm->GetLevi (mass, 0.4, 750, 3);
    func->SetParLimits(1,0.0001,20000);
    slopePar = 2;
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        auto hResult = YieldMean(hspectra_stat[imultibin], hspectra_sys_nocor[imultibin],
                            func, 0.05, 10, 0.01, 0.1, fitoption.Data(), Form("log_%s.root",bininfo.Data()),
                            0.8, 8.8);
        auto houty = new TH1D(Form("hYield_%.2f-%.2f",multibin[imultibin][0],
                    multibin[imultibin][1]), "", 6, 0, 6);
        houty->SetBinContent(1, hResult->GetBinContent(1)); // central value
        houty->SetBinContent(2, hResult->GetBinContent(2)); // stat.err.
        houty->SetBinContent(3, fitsyserrory[imultibin]); // fit syst. err.
        houty->SetBinContent(4, hResult->GetBinContent(3)); // extreme high
        houty->SetBinContent(5, hResult->GetBinContent(4)); // extreme low
        houty->SetBinContent(6, hResult->GetBinContent(9)); // extrapolated yield
        resultfile->cd();
        houty->Write(Form("hYield_%.2f-%.2f",multibin[imultibin][0],
                    multibin[imultibin][1]));
        
        auto houtm = new TH1D(Form("hMeanPt_%.2f-%.2f",multibin[imultibin][0],
                    multibin[imultibin][1]), "", 5, 0, 5);
        houtm->SetBinContent(1, hResult->GetBinContent(5)); // central value
        houtm->SetBinContent(2, hResult->GetBinContent(6)); // stat.err.
        houtm->SetBinContent(3, fitsyserrorm[imultibin]); // fit syst. err.
        houtm->SetBinContent(4, hResult->GetBinContent(7)); // extreme high
        houtm->SetBinContent(5, hResult->GetBinContent(8)); // extreme low
        resultfile->cd();
        houtm->Write(Form("hMeanPt_%.2f-%.2f",multibin[imultibin][0],
                    multibin[imultibin][1]));
    }
    resultfile->Close();
    gSystem->Exit(1);
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