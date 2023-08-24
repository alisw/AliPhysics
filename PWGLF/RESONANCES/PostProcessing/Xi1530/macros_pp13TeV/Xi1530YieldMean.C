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
vector<double> GetPidNdetawithError(double multi_start, double multi_end);
TString finalfile;
TString multistring;
bool isINEL = false;
bool skipcustomfit = false;
//TString finalfile = "AnalysisResults_Xi1530_systematic0-100_0-10_10-30_30-50_50-70_70-100_0.00-0.01_0.01-0.05_0.05-0.10_bins.root";
TString workdirectory = "/Users/blim/alidock/Postprocessing/data/";
enum {kFitExpPt=1, kFitLevi, fFitExpMt, kFitBoltzmann, kFitBlastWave, kFitBoseEinstein, kFitFermiDirac};
vector<TString> functions = {"", "kFitExpPt", "kFitLevi", "fFitExpMt", "kFitBoltzmann", "kFitBlastWave", "kFitBoseEinstein", "kFitFermiDirac"};

vector<vector<vector<double>>> // { {first pT bin variation}, {last pTbin variation} } // for each function
            functionFitRangeVar = { { {}, {} },            // null
                                    { {0.8}, {2.0, 2.4, 3.2} }, // ExpPt, only front part.
                                    { {0.8}, {2.4, 3.2, 4.0, 4.8, 5.6, 8.8} }, // Levy, all region.
                                    { {0.8}, {2.0, 2.4, 3.2} }, // ExpMt, only front part.
                                    { {0.8}, {2.0, 2.4, 3.2} }, // Boltzmann, only front part.
                                    { {0.8}, {2.4, 3.2, 4.0, 4.8, 5.6, 8.8} }, // BGBW, all region.
                                    { {0.8}, {2.0, 2.4, 3.2} }, // BoseEinstein, only front part.
                                    { {0.8}, {2.0, 2.4, 3.2} } // FermiDirac, only front part.
                                    };
vector<vector<double>> centralvalue = {
    {0.8, 3.2},  // 0: INEL
    {0.8, 3.2},  // 1: 0-10
    {0.8, 4.0},  // 2: 10-30
    {0.8, 5.6},  // 3: 30-50
    {0.8, 2.4},  // 4: 50-70
    {0.8, 4.0},  // 5: 70-100
    {0.8, 4.8},  // 6: 0-0.01 (HM)
    {0.8, 3.2},   // 7: 0.01-0.05 (HM)
    {0.8, 3.2}   // 8: 0.05-0.1 (HM)
};
Int_t maxtrial = 10000;
double mass = 1.5318;
TString fitoption = "0qi"; //default "0q"

void Xi1530YieldMean(double multis = 0, double multie = 100){
    int centralvaluebin = 0;
    TString bininfo = "_";
    if ((multis == 0) && (multie == 0)){ 
        //INEL case
        bininfo += "INEL";
        isINEL = true;
        multistring = "INEL";
    }
    else if (multie < 0.5) {
        bininfo += Form("%.2f", multis);
        bininfo += "-";
        bininfo += Form("%.2f", multie);
        multistring = Form("%.2f-%.2f", multis, multie);
    }
    else {
        bininfo += Form("%.0f", multis);
        bininfo += "-";
        bininfo += Form("%.0f", multie);   
        multistring = Form("%.2f-%.2f", multis, multie);
    }
    double binsum = multis + multie;
    if (binsum < 0.011) // 0-0.01
        centralvaluebin = 6;
    else if (binsum < 0.061) // 0.01-0.05
        centralvaluebin = 7;
    else if (binsum < 0.151)  // 0.05-0.1
        centralvaluebin = 8;
    else if (binsum < 11) // 0-10
        centralvaluebin = 1;
    else if (binsum < 41) // 10-30
        centralvaluebin = 2;
    else if (binsum < 81) // 30-50
        centralvaluebin = 3;
    else if (binsum < 101) // 0-100
        centralvaluebin = 0;
    else if (binsum < 121) // 50-70
        centralvaluebin = 4;
    else if (binsum < 171) // 70-100
        centralvaluebin = 5;
    else
        centralvaluebin = 0;
    finalfile = Form("./AnalysisResults_Xi1530_systematic%s.root", bininfo.Data());

    // Multibin Customizing
    vector<vector<double>> functionFitCustom = {
        {0},         // null
        {0},         // ExpPt, only front part.
        {0},         // Levy, all region.
        {0},         // ExpMt, only front part.
        {0},         // Boltzmann, only front part.
        {0},         // BGBW, all region.
        {0},         // BoseEinstein, only front part.
        {0}          // FermiDirac, only front part.
    };
    if (skipcustomfit)
        cout << "skip custom" << endl;
    else if (binsum < 0.011)  // 0-0.01
        cout << "skip custom" << endl;
    else if (binsum < 0.061)  // 0.01-0.05
        functionFitCustom = {
            {0},    // null
            {0},    // ExpPt, only front part.
            {2.4},  // Levy, all region.
            {0},    // ExpMt, only front part.
            {0},    // Boltzmann, only front part.
            {0},    // BGBW, all region.
            {0},    // BoseEinstein, only front part.
            {0}     // FermiDirac, only front part.
        };
    else if (binsum < 0.151)  // 0.05-0.1
        functionFitCustom = {
            {0},         // null
            {0},  // ExpPt, only front part.
            {2.4},  // Levy, all region.
            {0},         // ExpMt, only front part.
            {0},         // Boltzmann, only front part.
            {0},         // BGBW, all region.
            {0},         // BoseEinstein, only front part.
            {0}          // FermiDirac, only front part.
        };
    else if (binsum < 11)  // 0-10
        functionFitCustom = {
            {0},         // null
            {2.4, 3.2},       // ExpPt, only front part.
            {2.4, 5.6},  // Levy, all region.
            {0},         // ExpMt, only front part.
            {0},         // Boltzmann, only front part.
            {0},         // BGBW, all region.
            {0},         // BoseEinstein, only front part.
            {0}          // FermiDirac, only front part.
        };
    else if (binsum < 41)  // 10-30
        functionFitCustom = {
            {0},         // null
            {2.0, 2.4, 3.2},  // ExpPt, only front part.
            {0},         // Levy, all region.
            {0},         // ExpMt, only front part.
            {0},         // Boltzmann, only front part.
            {0},         // BGBW, all region.
            {0},         // BoseEinstein, only front part.
            {0}          // FermiDirac, only front part.
        };
    else if (binsum < 81)  // 30-50
        functionFitCustom = {
            {0},              // null
            {2.0, 2.4, 3.2},  // ExpPt, only front part.
            {0},              // Levy, all region.
            {0},              // ExpMt, only front part.
            {0},              // Boltzmann, only front part.
            {0},              // BGBW, all region.
            {0},              // BoseEinstein, only front part.
            {0}               // FermiDirac, only front part.
        };
    else if (binsum < 101)  // 0-100
        functionFitCustom = {
            {0},              // null
            {2.0, 2.4, 3.2},  // ExpPt, only front part.
            {0},              // Levy, all region.
            {0},              // ExpMt, only front part.
            {0},              // Boltzmann, only front part.
            {0},              // BGBW, all region.
            {0},              // BoseEinstein, only front part.
            {0}               // FermiDirac, only front part.
        };
    else if (binsum < 121)  // 50-70
        functionFitCustom = {
            {0},              // null
            {2.0, 2.4, 3.2},  // ExpPt, only front part.
            {0},              // Levy, all region.
            {0},              // ExpMt, only front part.
            {0},              // Boltzmann, only front part.
            {0},              // BGBW, all region.
            {0},              // BoseEinstein, only front part.
            {0}               // FermiDirac, only front part.
        };
    else if (binsum < 171)  // 70-100
        functionFitCustom = {
            {0},              // null
            {2.0, 2.4, 3.2},  // ExpPt, only front part.
            {0},              // Levy, all region.
            {0},              // ExpMt, only front part.
            {0},              // Boltzmann, only front part.
            {0},              // BGBW, all region.
            {0},              // BoseEinstein, only front part.
            {0}               // FermiDirac, only front part.
        };
    else
        cout << "skip custom" << endl;

    for (int fitvar = 1; fitvar < functions.size(); fitvar++) {
        for (int binvar = 0; binvar < functionFitCustom[fitvar].size();
             binvar++) {
            functionFitRangeVar[fitvar][1].erase(
                std::remove(functionFitRangeVar[fitvar][1].begin(),
                            functionFitRangeVar[fitvar][1].end(),
                            functionFitCustom[fitvar][binvar]),
                functionFitRangeVar[fitvar][1].end());
            cout << "Fit function(" << functions[fitvar] << " bin "
                 << functionFitCustom[fitvar][binvar] << "removed!" << endl;
        }
    }

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

    bool recalculate7TeV = false;
    TString filename;
    if(recalculate7TeV)
        filename = "AnalysisResults_Xi1530_YieldMean_7TeV.root";
    else
        filename = Form("AnalysisResults_Xi1530_YieldMean%s.root", bininfo.Data());

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
                    Form("%s_SYS_corrected",
                        multistring.Data()
                        ));
            hspectra_stat[multbin]->Write(
                    Form("%s_stat_corrected",
                        multistring.Data()
                        ));
            hspectra_sys_nocor[multbin]->Write(
                    Form("%s_SYS_uncor_corrected",
                        multistring.Data()
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
    double fitchi2, fitNDF;
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
                                    func, 0.05, 10, 0.01, 0.1, fitoption.Data(), Form("log%s.root",bininfo.Data()),fitstart, fitend);
                    yield[imultibin].push_back(hfinal->GetBinContent(1));
                    meanpt[imultibin].push_back(hfinal->GetBinContent(5));
                    func->SetLineWidth(1);
                    func->SetLineStyle(flinestyle);
                    func->SetRange(0,fitend);
                    resultfile->cd();
                    func->Write(
                        Form("%s_%s_%.2f-%.2f",
                        multistring.Data(),
                        functions[fitvar].Data(),
                        fitstart,fitend
                        ));
                    fitchi2 = func->GetChisquare();
                    fitNDF = func->GetNDF();
                    cout << "Fit Result!!: " << fitchi2 / fitNDF
                         << ", fitchi2: " << fitchi2 << ", fitNDF: " << fitNDF
                         << endl;
                    funcs[imultibin].push_back(func);
                    flinestyle++;
                }
            }
        }
    }

    int sysColorPallet = GetSerialColors(functions.size());
    // for memo, small
    TLatex* tm = new TLatex();
    tm->SetNDC();
    tm->SetTextSize(0.03);
    auto legendhead = new TLegend(.4, .65, .95, .9);
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
            new TH1D("systematic_yields", "", 100, *ymin*0.75 , *ymax*1.4);
        systematicyields->GetXaxis()->SetTitle("<d#it{N}/d#it{y}>");
        systematicyields->GetYaxis()->SetTitle("Counts");
        TH1* systematicmeanpTs =
            new TH1D("systematic_meanpTs", "", 100, *mmin*0.75 , *mmax*1.4);
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
                    fitnames[loopcheck] = Form(
                        "%s(dydn %.2f, mpt: %.2f)", fitnames[loopcheck].Data(),
                        1e2 * yield[imultibin][loopcheck],
                        meanpt[imultibin][loopcheck]);
                    funcs[imultibin][loopcheck]->SetLineColor(sysColorPallet + functions.size() - fitvar);
                    funcs[imultibin][loopcheck]->Draw("same");
                    legendhead->AddEntry(funcs[imultibin][loopcheck], fitnames[loopcheck], "L");
                    
                    cout << "yields: " << yield[imultibin][loopcheck] << " meanpT" << meanpt[imultibin][loopcheck]  <<  ", Chi^2: " << funcs[imultibin][loopcheck]->GetChisquare() << ", NDF: " << funcs[imultibin][loopcheck]->GetNDF() << ", Chi^2/NDF: " << funcs[imultibin][loopcheck]->GetChisquare()/funcs[imultibin][loopcheck]->GetNDF() << endl;

                    if( (fitvar == kFitLevi)
                        && (fitstart == centralvalue[centralvaluebin][0])
                        && (fitend == centralvalue[centralvaluebin][1])
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
                   Form("Spectrafit_%s_zoom", multistring.Data()),"figs/LowpTExtrapolate/");
        resultfile->cd();
        cResults->Write(
                Form("Spectrafit_%s_zoom",
                    multistring.Data()));
        systematicyields->Write(
                Form("Yieldsys_%s",
                    multistring.Data()));
        systematicmeanpTs->Write(
                Form("MeanPtsys_%s",
                    multistring.Data()));
        cResults->SetLogy(false);
        systematicyields->Draw();
        TLine* centralarrow1 = new TLine(centralyield, 0, centralyield, 
                                        systematicyields->GetMaximum());
        centralarrow1->SetLineColor(kRed);
        centralarrow1->SetLineStyle(2);
        centralarrow1->Draw();
        SaveCanvas(cResults,
                   Form("Yields_systematic_%s", multistring.Data()),"figs/LowpTExtrapolate/");
        cout << "mean: " << systematicyields->GetMean() << ", stdev: " << systematicyields->GetRMS() 
             << ", percentage: " << systematicyields->GetRMS()/systematicyields->GetMean()  << endl;
        systematicmeanpTs->Draw();
        TLine* centralarrow2 = new TLine(centralmeanpT, 0, centralmeanpT, 
                                        systematicyields->GetMaximum());
        centralarrow2->SetLineColor(kRed);
        centralarrow2->SetLineStyle(2);
        centralarrow2->Draw();
        SaveCanvas(cResults,
                   Form("MeanPts_systematic_%s", multistring.Data()),"figs/LowpTExtrapolate/");
        fitsyserrory.push_back(systematicyields->GetRMS());
        fitsyserrorm.push_back(systematicmeanpTs->GetRMS());
    }

    // Systematic error
    func = fm->GetLevi (mass, 0.4, 750, 3);
    func->SetParLimits(1,0.0001,20000);
    slopePar = 2;
    for (int imultibin = 0; imultibin < multibin.size(); imultibin++) {
        auto hResult = YieldMean(hspectra_stat[imultibin], hspectra_sys_nocor[imultibin],
                            func, 0.05, 10, 0.01, 0.1, fitoption.Data(), Form("log%s.root",bininfo.Data()),
                            centralvalue[centralvaluebin][0], centralvalue[centralvaluebin][1]);
        auto houty = new TH1D(Form("hYield_%s",multistring.Data()), "", 6, 0, 6);
        houty->SetBinContent(1, hResult->GetBinContent(1)); // central value
        houty->SetBinContent(2, hResult->GetBinContent(2)); // stat.err.
        houty->SetBinContent(3, fitsyserrory[imultibin]); // fit syst. err.
        houty->SetBinContent(4, hResult->GetBinContent(3)); // extreme high
        houty->SetBinContent(5, hResult->GetBinContent(4)); // extreme low
        houty->SetBinContent(6, hResult->GetBinContent(9)); // extrapolated yield
        resultfile->cd();
        houty->Write(Form("hYield_%s",multistring.Data()));
        cout << "dNdy: Central value: " << hResult->GetBinContent(1) << ", stat.e: " << hResult->GetBinContent(2)
             << ", fit sys err.: " << fitsyserrory[imultibin] << ", extreme high: " << hResult->GetBinContent(3)
             << ", extreme low: " << hResult->GetBinContent(4) << endl;
        
        auto houtm = new TH1D(Form("hMeanPt_%s",multistring.Data()), "", 5, 0, 5);
        houtm->SetBinContent(1, hResult->GetBinContent(5)); // central value
        houtm->SetBinContent(2, hResult->GetBinContent(6)); // stat.err.
        houtm->SetBinContent(3, fitsyserrorm[imultibin]); // fit syst. err.
        houtm->SetBinContent(4, hResult->GetBinContent(7)); // extreme high
        houtm->SetBinContent(5, hResult->GetBinContent(8)); // extreme low
        resultfile->cd();
        houtm->Write(Form("hMeanPt_%s",multistring.Data()));

        cout << "Mean pT: Central value: " << hResult->GetBinContent(5) << ", stat.e: " << hResult->GetBinContent(6)
             << ", fit sys err.: " << fitsyserrorm[imultibin] << ", extreme high: " << hResult->GetBinContent(7)
             << ", extreme low: " << hResult->GetBinContent(8) << endl;
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
    TString text;
    if(isINEL)
        text = "hSpectra_INEL_sys";
    else
        text = Form("hSpectra_%.2f_%.2f_sys", multi_start, multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(text);

    return hr;
}
TH1D* GetSpectrasysNocor(double multi_start, double multi_end) {
    TFile* inputfile = new TFile(finalfile.Data());
    TString text;
    if(isINEL)
        text = "hSpectra_INEL_sys_noCorrelation";
    else
        text = Form("hSpectra_%.2f_%.2f_sys_noCorrelation", multi_start, multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(text);

    return hr;
}
TH1D* GetSpectrastat(double multi_start, double multi_end){
    TFile* inputfile = new TFile(finalfile.Data());
    TString text;
    if(isINEL)
        text = "hSpectra_INEL_stat";
    else
        text = Form("hSpectra_%.2f_%.2f_stat",multi_start,multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(text);

    return hr;
}

vector<double> GetPidNdetawithError(double multi_start, double multi_end){
    // Return pion's dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}
    TFile* fPi =
        TFile::Open("OutputYields_pp13mult.root","READ");
    auto pisys = (TGraphAsymmErrors*)fPi->Get("pi_Syst");
    vector<double> returnarray;
    int lenpisys = pisys->GetN();

    //--dNdeta histogram of pion+ + pion-
    // Ref: given by Anders, need to update ref.
    // Error was stat+sys error, so choosed bigger error(sym) error.

    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta = {0};
    vector<double> dNchdeta_e = {0};
    for(int i = 0; i < lenpisys; i++){
        dNchdeta.push_back(pisys->GetY()[lenpisys-i]);
        dNchdeta_e.push_back(pisys->GetErrorYhigh(lenpisys-i));
    }

    /*
    vector<double> dNchdeta = 
    {0, 24.605455025, 19.016207266, 15.477779961, 13.241832514, 11.612550017, 9.742647922, 7.779575692,  6.241633459, 4.530678113, 2.713659699};
    vector<double> dNchdeta_e = 
    {0,  1.121689195,  0.856354796,  0.685848455,  0.582627504,   0.498773083,  0.415988997, 0.327417792, 0.261067034, 0.186885663, 0.110000678};
    */
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