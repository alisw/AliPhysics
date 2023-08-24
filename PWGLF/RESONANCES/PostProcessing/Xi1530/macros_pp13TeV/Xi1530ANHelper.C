// This will make a several tables for AN

TString workdirectory = "data/";

TH1* GethYeilds(double multi_start,
                double multi_end,
                TString inputOptions = "Default");
TH1* GetEfficiency(double multi_start,
                double multi_end,
                TString inputOptions = "MCCheck");
TH1* GetfNorm(TString inputOptions = "LHC17_GenMC_novertexer", int option = 1);
TH1* GetfVtx(TString inputOptions = "LHC17_GenMC_novertexer", int option = 1);
TH1* GetResults(double multi_start,
                double multi_end,
                TString inputOptions = "Yield"); //Yield, MeanPt
void PrintRawYields(vector<vector<double>> multibincheck);
void PrintdNdympT(vector<vector<double>> multibincheck);
void Printfnorm();
void Printfvtx();
void Xi1530ANHelper(){
    vector<vector<double>> multibincheck = {
        {0, 10},  {10, 30}, {30, 50},  {50, 70},  {70, 100},
        {0, 100}, {0, 0},   {0, 0.01}, {0.01, 0.05}, {0.05, 0.1}};

    //PrintRawYields(multibincheck);

    PrintdNdympT(multibincheck);

    //Printfnorm();
    //Printfvtx();

}
TH1* GethYeilds(double multi_start,
                double multi_end,
                TString inputOptions) {
    vector<double> multibin = {multi_start, multi_end};
    TString bininfo;
    bool isINEL = false;
    if ((multi_start == 0) && (multi_end == 0)){ 
        //INEL case
        bininfo += "INEL";
        isINEL=true;
    }
    else {
        bininfo += Form("%.2f", multi_start);
        bininfo += "-";
        bininfo += Form("%.2f", multi_end);
    }
    int cutbin = 1;
    TString finputfile; 
    if(isINEL)
        finputfile =
        Form("%sAnalysisResults_Extracted_%i_%s_%s1.root",
             workdirectory.Data(), cutbin, bininfo.Data(),
             inputOptions.Data());
    else
        finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%s_%s1.root",
             workdirectory.Data(), cutbin, bininfo.Data(),
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1* base = (TH1*)inputfile->Get("hDataRawYield");

    return base;
}
TH1* GetResults(double multi_start,
                double multi_end,
                TString inputOptions) {
    vector<double> multibin = {multi_start, multi_end};
    TString bininfo;
    bool isINEL = false;
    if ((multi_start == 0) && (multi_end == 0)){ 
        //INEL case
        bininfo += "INEL";
        isINEL=true;
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
    int cutbin = 1;
    TString finputfile =
        Form("AnalysisResults_Xi1530_YieldMean_%s.root",bininfo.Data());
    TFile* inputfile = new TFile(finputfile.Data());

    TString histname;
    if (isINEL)
        histname = Form("h%s_INEL", inputOptions.Data());
    else
        histname = Form("h%s_%.2f-%.2f", inputOptions.Data(), multi_start, multi_end);
    TH1D* hr = (TH1D*)inputfile->Get(histname);

    return hr;
}
TH1* GetfNorm(TString inputOptions, int option) {
    TString finputfile =
        Form("AnalysisResults_Xi1530MCQA_%s.root",inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    if(option == 1){
        TH1* base = (TH1*)inputfile->Get("hTrigEffi");
        return base;
    }
    if(option == 2){
        TH1* base = (TH1*)inputfile->Get("hTrigEffiINEL");
        return base;
    }
    else{
        TH1* base = (TH1*)inputfile->Get("hTrigEffiMB");
        return base;
    }   
}
TH1* GetfVtx(TString inputOptions, int option) {
    TString finputfile =
        Form("AnalysisResults_Xi1530MCQA_%s.root",inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    if(option == 1){
        TH1* base = (TH1*)inputfile->Get("hVertexEffi");
        return base;
    }
    if(option == 2){
        TH1* base = (TH1*)inputfile->Get("hVertexEffiINEL");
        return base;
    }
    else{
        TH1* base = (TH1*)inputfile->Get("hVertexEffiMB");
        return base;
    } 
}
TH1* GetEfficiency(double multi_start,
                double multi_end,
                TString inputOptions) {
    vector<double> multibin = {multi_start, multi_end};
    TString bininfo;
    bool isINEL = false;
    if ((multi_start == 0) && (multi_end == 0)){ 
        //INEL case
        bininfo += "INEL";
        isINEL=true;
    }
    else {
        bininfo += Form("%.2f", multi_start);
        bininfo += "-";
        bininfo += Form("%.2f", multi_end);
    }
    TString finputfile; 
    int cutbin = 1;
    if(isINEL)
        finputfile =
        Form("%sAnalysisResults_Extracted_%i_%s_%s1.root",
             workdirectory.Data(), cutbin, bininfo.Data(),
             inputOptions.Data());
    else
        finputfile =
        Form("%sAnalysisResults_Extracted_%i_Multi_%s_%s1.root",
             workdirectory.Data(), cutbin, bininfo.Data(),
             inputOptions.Data());
    TFile* inputfile = new TFile(finputfile.Data());
    TH1* base = (TH1*)inputfile->Get("hMCReconEffi");

    return base;
}
void PrintRawYields(vector<vector<double>> multibincheck){
    vector<TH1*> hyields;
    vector<TString> checklist = {"Default", "LikeSignBkg", "BkgFit", "BinCountLm"};

    auto hBase = GethYeilds(multibincheck[0][0],
                            multibincheck[0][1], checklist[0].Data());
    const TArrayD* ptbinarray = hBase->GetXaxis()->GetXbins();
    vector<double> ptbin;
    for (int i = 0; i < ptbinarray->GetSize(); i++)
        ptbin.push_back(ptbinarray->GetAt(i));
    TString fvarfile;
    
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        hyields.clear();
        for (auto const& check : checklist) {
            auto htemp = GethYeilds(multibincheck[imultibin][0],
                                    multibincheck[imultibin][1], check.Data());
            hyields.push_back(htemp);
        }
        cout << "\\begin{table}[]" << endl;
        cout << "\\centering" << endl;
        cout << "\\begin{tabular}{|l|l|l|l|l|}" << endl;
        cout << "\\hline" << endl;
        cout << "\\ensuremath{p_{\\mathrm{T}}} & Default(Event Mixing)($\\times \\num{e3}$) & Like-Sign Bkg($\\times \\num{e3}$) & BkgFit($\\times \\num{e3}$) & BinCount(Lm-1)($\\times \\num{e3}$) \\\\ \\hline"
             << endl;
        for (int bin = 2; bin < hyields[0]->GetNbinsX(); bin++) {
            cout << "$" << Form("%.1f", ptbin[bin - 1]) << " - "
                 << Form("%.1f", ptbin[bin]) << "$ & $"
                 << Form("%.2f", hyields[0]->GetBinContent(bin)*1e-3) << " \\pm "
                 << Form("%.2f", hyields[0]->GetBinError(bin)*1e-3) << "$ & $"
                 << Form("%.2f", hyields[1]->GetBinContent(bin)*1e-3) << " \\pm "
                 << Form("%.2f", hyields[1]->GetBinError(bin)*1e-3) << "$ & $"
                 << Form("%.2f", hyields[2]->GetBinContent(bin)*1e-3) << " \\pm "
                 << Form("%.2f", hyields[2]->GetBinError(bin)*1e-3) << "$ & $"
                 << Form("%.2f", hyields[3]->GetBinContent(bin)*1e-3) << " \\pm "
                 << Form("%.2f", hyields[3]->GetBinError(bin)*1e-3) << "$ \\\\" << endl;
        }
        cout << "\\hline" << endl;
        cout << "\\end{tabular}" << endl;
        cout << "\\caption{Rawyield in " << multibincheck[imultibin][0] << " - "
             << multibincheck[imultibin][1] << " bin}" << endl;
        cout << "\\label{tab:rawyield" << multibincheck[imultibin][0]
             << multibincheck[imultibin][1] << "}" << endl;
        cout << "\\end{table}" << endl;
    }

}
void PrintdNdympT(vector<vector<double>> multibincheck){
    vector<TString> checklist_result = {"Yield", "MeanPt"};
    vector<TH1*> hResults;

    cout << "\\begin{table}[]" << endl;
        cout << "\\centering" << endl;
        cout << "\\begin{tabular}{|l|l|l|l|}" << endl;
        cout << "\\hline" << endl;
        cout << "Multiplicity & ${\\rm d}N/{\\rm d}y$ ($\\times \\num{e-3}$) & $\\mpt$ & ratio [\\%] \\\\ \\hline"
             << endl;
    for (int imultibin = 0; imultibin < multibincheck.size(); imultibin++) {
        hResults.clear();
        for (auto const& check : checklist_result) {
            auto htemp = GetResults(multibincheck[imultibin][0],
                                    multibincheck[imultibin][1], check.Data());
            hResults.push_back(htemp);
        }
        double tempfinalsyserror = 0;
        tempfinalsyserror += pow(hResults[0]->GetBinContent(3),2);
        tempfinalsyserror += pow(hResults[0]->GetBinContent(4),2);

        double averageerror_correl = sqrt(tempfinalsyserror)/2;
        double averageerror = sqrt(pow(averageerror_correl,2) + pow(0.0538*hResults[0]->GetBinContent(1),2));
        double averageerror_uncorr = 0.0538*hResults[0]->GetBinContent(1);

        double mtempfinalsyserror = 0;
        mtempfinalsyserror += pow(hResults[1]->GetBinContent(3),2);
        mtempfinalsyserror += pow(hResults[1]->GetBinContent(4),2);

        double averagpteerror_correl = sqrt(mtempfinalsyserror)/2;
        double averagpteerror = sqrt(pow(averagpteerror_correl,2) + pow(0.0538*hResults[1]->GetBinContent(1),2));
        double averagpteerror_uncorr = 0.0538*hResults[1]->GetBinContent(1);

        if(imultibin == multibincheck.size()-1) cout << "\\hline \\hline" << endl;
        cout << "$" << multibincheck[imultibin][0] << "$ - $"
             << multibincheck[imultibin][1] << "$ & $" 
             << Form("%.3f", hResults[0]->GetBinContent(1)*1e3)
             << " \\pm \\num{" << Form("%.3f", hResults[0]->GetBinContent(2)*1e3)
             << "} \\pm \\num{" << Form("%.3f", averageerror*1e3 )
             << "} \\pm \\num{" << Form("%.3f", averageerror_uncorr*1e3 )
             << "}$ & $"
             << Form("%.3f", hResults[1]->GetBinContent(1))
             << " \\pm \\num{" << Form("%.3f", hResults[1]->GetBinContent(2))
             << "} \\pm \\num{" << Form("%.3f", averagpteerror )
             << "} \\pm \\num{" << Form("%.3f", averagpteerror_uncorr )
             << "}$ & " << 100*hResults[0]->GetBinContent(6)/hResults[0]->GetBinContent(1) <<  " \\\\" << endl;
    }
        cout << "\\hline" << endl;
        cout << "\\end{tabular}" << endl;
        cout << "\\caption{ ${\\rm d}N/{\\rm d}y$ and $\\mpt$ with each multiplicity bins and ratio of the extrapolated yield (value $\\pm$ stat.e $\\pm$ sys.e(total) $\\pm$ sys.e(uncorr) )}" << endl;
        cout << "\\label{tab:physicsreult}" << endl;
        cout << "\\end{table}" << endl;
}
void Printfnorm(){
    vector<TString> checklist_result = {"LHC16_GenMC_AOD310", "LHC17_GenMC_AOD309", "LHC18_GenMC_AOD308", "LHC15_LHC15g3a3_AOD311", "LHC15_LHC15g3c3_AOD311"};
    vector<vector<double>> multibin = {
        {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    vector<TH1*> hResults;
    vector<TH1*> hResultsMB;
    vector<TH1*> hResultsINEL;
    
    cout << "\\begin{table}[]" << endl;
        cout << "\\centering" << endl;
        cout << "\\begin{tabular}{|l|l|l|l|l|l|l|}" << endl;
        cout << "\\hline" << endl;
        cout << "Multiplicity & $f_{\\rm{norm}}$ $\\pm$ stat.e (LHC16) & (LHC17) & (LHC18) & (LHC15, PYTHIA 8) & (LHC15, PYTHIA 6) \\\\ \\hline"
             << endl;
    for (int checkbin = 0; checkbin < checklist_result.size(); checkbin++) {
        auto htemp = GetfNorm(checklist_result[checkbin].Data(),1);
        hResults.push_back(htemp);
        auto htempMB = GetfNorm(checklist_result[checkbin].Data(),3);
        hResultsMB.push_back(htempMB);
        auto htempINEL = GetfNorm(checklist_result[checkbin].Data(),2);
        hResultsINEL.push_back(htempINEL);
    }
    for (int j = 0; j < hResults[0]->GetNbinsX(); j++) {
        cout << "$" << multibin[j][0] << "$ - $"
             << multibin[j][1] << "$ & $" 
             << Form("%.3f", hResults[0]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[0]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[1]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[1]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[2]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[2]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[3]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[4]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[4]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[4]->GetBinError(j+1)*1e3)
             << "}$ \\\\" << endl;
    }
        cout << "\\hline \\hline" << endl;

        //MB
        cout << "$" << 0 << "$ - $"
             << 100 << "$ & $" 
             << Form("%.3f", hResultsMB[0]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[0]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[1]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[1]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[2]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[2]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[3]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[3]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[4]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[4]->GetBinError(1)*1e3)
             << "}$ \\\\ \\hline" << endl;

        cout << "\\end{tabular}" << endl;
        cout << "\\caption{ $f_{\\rm{norm}}$ on each multiplicity bins with different MC data, the errors are muliplied to $\\num{e3}$ }" << endl;
        cout << "\\label{tab:fnormINELg0}" << endl;
        cout << "\\end{table}" << endl;
}
void Printfvtx(){
    vector<TString> checklist_result = {"LHC16_GenMC_AOD310", "LHC17_GenMC_AOD309", "LHC18_GenMC_AOD308", "LHC15_LHC15g3a3_AOD311", "LHC15_LHC15g3c3_AOD311"};
    vector<vector<double>> multibin = {
        {0, 10}, {10, 30}, {30, 50}, {50, 70}, {70, 100}};
    vector<TH1*> hResults;
    vector<TH1*> hResultsMB;
    vector<TH1*> hResultsINEL;
    
    cout << "\\begin{table}[]" << endl;
        cout << "\\centering" << endl;
        cout << "\\begin{tabular}{|l|l|l|l|l|l|l|}" << endl;
        cout << "\\hline" << endl;
        cout << "Multiplicity & $f_{\\rm{vtx}}$ $\\pm$ stat.e (LHC16) & (LHC17) & (LHC18) & (LHC15, PYTHIA 8) & (LHC15, PYTHIA 6) \\\\ \\hline"
             << endl;
    for (int checkbin = 0; checkbin < checklist_result.size(); checkbin++) {
        auto htemp = GetfVtx(checklist_result[checkbin].Data(),1);
        hResults.push_back(htemp);
        auto htempMB = GetfVtx(checklist_result[checkbin].Data(),3);
        hResultsMB.push_back(htempMB);
        auto htempINEL = GetfVtx(checklist_result[checkbin].Data(),2);
        hResultsINEL.push_back(htempINEL);
    }
    for (int j = 0; j < hResults[0]->GetNbinsX(); j++) {
        cout << "$" << multibin[j][0] << "$ - $"
             << multibin[j][1] << "$ & $" 
             << Form("%.3f", hResults[0]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[0]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[1]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[1]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[2]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[2]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[3]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[3]->GetBinError(j+1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResults[4]->GetBinContent(j+1)) << " \\pm \\num{" << Form("%.2f", hResults[4]->GetBinError(j+1)*1e3)
             << "}$ \\\\" << endl;
    }
        cout << "\\hline \\hline" << endl;

        //MB
        cout << "$" << 0 << "$ - $"
             << 100 << "$ & $" 
             << Form("%.3f", hResultsMB[0]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[0]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[1]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[1]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[2]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[2]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[3]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[3]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsMB[4]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsMB[4]->GetBinError(1)*1e3)
             << "}$ \\\\ \\hline" << endl;
        //MB
        cout << "inel & $" 
             << Form("%.3f", hResultsINEL[0]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsINEL[0]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsINEL[1]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsINEL[1]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsINEL[2]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsINEL[2]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsINEL[3]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsINEL[3]->GetBinError(1)*1e3)
             << "}$ & $" 
             << Form("%.3f", hResultsINEL[4]->GetBinContent(1)) << " \\pm \\num{" << Form("%.2f", hResultsINEL[4]->GetBinError(1)*1e3)
             << "}$ \\\\ \\hline" << endl;

        cout << "\\end{tabular}" << endl;
        cout << "\\caption{ $f_{\\rm{vtx}}$ on each multiplicity bins with different MC data, the errors are muliplied to $\\num{e3}$ }" << endl;
        cout << "\\label{tab:fvtxINELg0}" << endl;
        cout << "\\end{table}" << endl;
}