#include "DrawingHelper.C"
int canvasCount = 1;
vector<int> Markerset={24, 21, 33, 34, 29, 2, 25, 26, 27};
TCanvas* GetCanvas(TString name = "Canvas", double w = 960, double h = 720);
bool isHM = false;
vector<TH1*> LoadAllMultQAHist(vector<TString> inputnames);
vector<TH1*> LoadAllVzQAHist(vector<TString> inputnames);
vector<TH1*> LoadAllV0MQAHist(vector<TString> inputnames);
vector<TH1*> LoadAllXiRawQAHist(vector<TString> inputnames);
vector<TH1*> LoadAllXiQAMeanHist(vector<TString> inputnames);
vector<TH1*> LoadAllXiQASigmaHist(vector<TString> inputnames);
TH1* LoadMultQAHist(TString inputnames);
TH1* LoadVzQAHist(TString inputnames);
TH1* LoadV0MQAHist(TString inputnames);
TH1* LoadXiRawQAHist(TString inputnames);
TH1* LoadXiQAMeanHist(TString inputnames);
TH1* LoadXiQASigmaHist(TString inputnames);

void DrawQAPlots(){
    TFile* resultfile_new =
        TFile::Open("AnalysisResults_Xi1530QA_MB.root","RECREATE");
    resultfile_new->Close();
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);

    vector<TString> periods15 = 
        {
         "pp_AOD439/LHC15f", 
         "pp_AOD439/LHC15i"
        };  
    vector<TString> periods16 = 
        {
         "pp_AOD441/LHC16d",
         "pp_AOD441/LHC16e",
         "pp_AOD441/LHC16g_p",
         "pp_AOD441/LHC16h",
         "pp_AOD441/LHC16i",
         "pp_AOD441/LHC16j",
         "pp_AOD441/LHC16k",
         "pp_AOD441/LHC16l",
         "pp_AOD441/LHC16o",
         "pp_AOD441/LHC16p"
        };
    vector<TString> periods17 = 
        {
         "pp_AOD441/LHC17c",
         "pp_AOD441/LHC17e",
         "pp_AOD441/LHC17f_p",
         "pp_AOD441/LHC17g",
         "pp_AOD441/LHC17h",
         "pp_AOD441/LHC17i",
         "pp_AOD441/LHC17j_p",
         "pp_AOD441/LHC17k",
         "pp_AOD441/LHC17l",
         "pp_AOD441/LHC17m",
         "pp_AOD441/LHC17o",
         "pp_AOD441/LHC17r"
        };
    vector<TString> periods18 = 
        {
         "pp_AOD441/LHC18b",
         "pp_AOD441/LHC18d",
         "pp_AOD441/LHC18e",
         "pp_AOD441/LHC18f",
         "pp_AOD441/LHC18g",
         "pp_AOD441/LHC18h",
         "pp_AOD441/LHC18i",
         "pp_AOD441/LHC18k",
         "pp_AOD441/LHC18l",
         "pp_AOD441/LHC18m_p",
         "pp_AOD441/LHC18n",
         "pp_AOD441/LHC18p"
        };
    // All List
    vector<vector<TString>> allperiods;
    //allperiods.push_back(periods15);
    //allperiods.push_back(periods16);
    //allperiods.push_back(periods17);
    allperiods.push_back(periods18);
    isHM = false;

    // Load QA Histogram
    vector<vector<TH1*>> AllMultQA;
    vector<vector<TH1*>> AllVzQA;
    vector<vector<TH1*>> AllXiRawQA;
    for (auto const& hList : allperiods){
        AllMultQA.push_back(LoadAllMultQAHist(hList));
        //AllVzQA.push_back(LoadAllVzQAHist(hList));
        //AllXiRawQA.push_back(LoadAllXiQASigmaHist(hList));
    }
    gSystem->Exit(1);
    
    //AllXiRawQA.push_back(LoadAllXiRawQAHist(periods15));

    //vector<TString> periodname = {"LHC15_FitSigma", "LHC16_FitSigma", "LHC17_FitSigma", "LHC18_FitSigma"};
    vector<TString> periodname = {"LHC15"};
    int periodn = 0;
    
    // Multiplicity
    
    for (auto const& hList : AllMultQA){
        TCanvas* cMultQA = GetCanvas("cMultQA", 960, 720);
        cMultQA->SetLeftMargin(0.10);
        cMultQA->Draw();
        cMultQA->SetLogy();

        int MultiColorPallet = GetSerialColors(hList.size());
        auto lMultQA = new TLegend(0.13, 0.13, 0.9, 0.7);
        lMultQA->SetFillStyle(0);
        lMultQA->SetNColumns(3);

        int colorp = 0;
        for (auto const& hInput : hList) {
            auto tempyield = hInput->GetEntries();
            hInput->GetXaxis()->SetTitle("Multiplicity Percentile (%)");
            hInput->GetYaxis()->SetTitle("# of Event(normalized)");
            hInput->SetMarkerStyle(Markerset[colorp]);
            hInput->SetMarkerSize(0.5);
            hInput->SetLineColor(MultiColorPallet + colorp);
            hInput->SetMarkerColor(MultiColorPallet + colorp);
            if(!isHM) hInput->Rebin(10);
            if(!isHM) 
                hInput->GetXaxis()->SetRangeUser(0,100);
            else
                hInput->GetXaxis()->SetRangeUser(0,0.1);

            if(tempyield > 100){
                double centralyield = 0;
                if(isHM){
                    centralyield = hInput->Integral(hInput->GetXaxis()->FindBin(0.0),
                                    hInput->GetXaxis()->FindBin(0.01));
                    hInput->Scale(1/centralyield);
                }
                else
                    hInput->Scale(1/tempyield);                
                lMultQA->AddEntry(hInput, Form("%s_Yield(%.1fM event)", allperiods[periodn][colorp].Data(), tempyield*1e-6), "PEL");
                cout << Form("%s_Yield(%.1fM event)", allperiods[periodn][colorp].Data(), tempyield*1e-6) << endl;
                if(isHM){
                    hInput->SetMaximum(8e-1);
                    hInput->SetMinimum(1e-4);
                }
                else{
                    hInput->SetMaximum(8e-2);
                    hInput->SetMinimum(1e-5);   
                }
                if(colorp < 1) hInput->Draw();
                else hInput->Draw("same");
            }
            else{
                lMultQA->AddEntry(hInput, Form("%s_Yield(No HM trigger)", allperiods[periodn][colorp].Data()), "PEL");
            }
            colorp++;
            }
        
        
        lMultQA->Draw();
        if(isHM)
            SaveCanvas(cMultQA,Form("cQA_Multi_%s_HM",periodname[periodn].Data()),"figs/");
        else
            SaveCanvas(cMultQA,Form("cQA_Multi_%s",periodname[periodn].Data()),"figs/");
        periodn++;
    }

    /*
    // Vz
    periodn = 0;
    for (auto const& hList : AllVzQA){
        TCanvas* cVzQA = GetCanvas("cVzQA", 960, 720);
        cVzQA->SetLeftMargin(0.10);
        cVzQA->Draw();

        int MultiColorPallet = GetSerialColors(hList.size());
        auto lVzQA = new TLegend(0.13, 0.13, 0.9, 0.7);
        lVzQA->SetFillStyle(0);
        lVzQA->SetNColumns(3);

        int colorp = 0;
        for (auto const& hInput : hList) {
            auto tempyield = hInput->GetEntries();
            hInput->Scale(1/tempyield);
            hInput->SetLineColor(MultiColorPallet + colorp);
            hInput->SetMarkerColor(MultiColorPallet + colorp);
            //hInput->Rebin(10);
            //hInput->SetMinimum(7e-3);
            //hInput->GetXaxis()->SetRangeUser(0,100);
            hInput->SetMaximum(12e-3);

            lVzQA->AddEntry(hInput, Form("%s_Yield(%.1fM event)", allperiods[periodn][colorp].Data(), tempyield*1e-6), "EL");
            cout << Form("%s_Yield(%.1fM event)", allperiods[periodn][colorp].Data(), tempyield*1e-6) << endl;

            if(colorp < 1) hInput->Draw();
            else hInput->Draw("same");
            colorp++;
        }
        
        lVzQA->Draw();
        SaveCanvas(cVzQA,Form("cQA_Vz_%s",periodname[periodn].Data()),"figs/");
        periodn++;
    }
    */
    /*
    // Xi Raw
    periodn = 0;
    for (auto const& hList : AllXiRawQA){
        TCanvas* cVzQA = GetCanvas("cXiRawQA", 960, 720);
        cVzQA->SetLeftMargin(0.10);
        cVzQA->Draw();
        //cVzQA->SetLogy();

        int MultiColorPallet = GetSerialColors(hList.size());
        auto lVzQA = new TLegend(0.3, 0.13, 0.7, 0.5);
        lVzQA->SetFillStyle(0);
        lVzQA->SetNColumns(3);

        int colorp = 0;
        for (auto const& hInput : hList) {
            hInput->SetLineColor(MultiColorPallet + colorp);
            hInput->SetMarkerColor(MultiColorPallet + colorp);
            //hInput->Rebin(10);
            //hInput->SetMinimum(7e-3);
            //hInput->GetXaxis()->SetRangeUser(0,100);
            hInput->SetMaximum(0.01);
            hInput->SetMinimum(0.001);

            lVzQA->AddEntry(hInput, Form("%s", allperiods[periodn][colorp].Data()), "EL");
            cout << Form("%s", allperiods[periodn][colorp].Data()) << endl;

            if(colorp < 1) hInput->Draw();
            else hInput->Draw("same");
            colorp++;
        }
        
        lVzQA->Draw();
        SaveCanvas(cVzQA,Form("cQA_XiRaw_%s",periodname[periodn].Data()),"figs/");
        periodn++;
    }
    */


}
vector<TH1*> LoadAllMultQAHist(vector<TString> inputnames){
    vector<TH1*> vReturn;
    for(int i=0;i<inputnames.size();i++)
        vReturn.push_back(LoadMultQAHist(inputnames[i]));

    return vReturn;

}
vector<TH1*> LoadAllVzQAHist(vector<TString> inputnames){
    vector<TH1*> vReturn;
    for(int i=0;i<inputnames.size();i++)
        vReturn.push_back(LoadVzQAHist(inputnames[i]));

    return vReturn;

}
vector<TH1*> LoadAllV0MQAHist(vector<TString> inputnames){
    vector<TH1*> vReturn;
    for(int i=0;i<inputnames.size();i++)
        vReturn.push_back(LoadV0MQAHist(inputnames[i]));

    return vReturn;

}
vector<TH1*> LoadAllXiRawQAHist(vector<TString> inputnames){
    vector<TH1*> vReturn;
    for(int i=0;i<inputnames.size();i++)
        vReturn.push_back(LoadXiRawQAHist(inputnames[i](7,13)));

    return vReturn;

}
vector<TH1*> LoadAllXiQAMeanHist(vector<TString> inputnames){
    vector<TH1*> vReturn;
    for(int i=0;i<inputnames.size();i++)
        vReturn.push_back(LoadXiQAMeanHist(inputnames[i](7,13)));

    return vReturn;

}
vector<TH1*> LoadAllXiQASigmaHist(vector<TString> inputnames){
    vector<TH1*> vReturn;
    for(int i=0;i<inputnames.size();i++)
        vReturn.push_back(LoadXiQASigmaHist(inputnames[i](7,13)));

    return vReturn;

}
TH1* LoadMultQAHist(TString inputnames){
    bool save = true;
    TFile* resultfile =
        TFile::Open("AnalysisResults_Xi1530QA_MB.root","UPDATE");
    TH1* vReturn;
    
    
    TFile* tempfile = new TFile(Form(
    "traindata/%s.root",inputnames.Data()));
    cout << Form(
    "traindata/%s.root",inputnames.Data()) << endl;
    TList* templist;
    if(isHM)
        templist = (TList*)tempfile->Get("Xi1530HM");
    else
        templist = (TList*)tempfile->Get("Xi1530MB");
    auto temphist = (TH1*)templist->FindObject("hMult_QA_onlyMult");
    temphist->SetTitle(Form("hMultQA_%s",inputnames.Data()));
    vReturn = (TH1*)temphist->Clone();  // signal
    resultfile->cd();
    vReturn->Write(Form("hMultQA_%s",inputnames.Data()));
    resultfile->Close();
    
    if(save){
        TCanvas* cSave = GetCanvas();
        cSave->Draw();
        vReturn->GetXaxis()->SetTitle("Multiplicity Percentile (%)");
        vReturn->GetYaxis()->SetTitle("# of Event");
        if(!isHM)
            vReturn->Rebin(10);
        vReturn->Draw();
        TString savename;
        if(isHM)
            savename = Form("cQA_Multi_%s_HM",inputnames(10,13).Data());
        else
            savename = Form("cQA_Multi_%s",inputnames(10,13).Data());
        SaveCanvas(cSave,savename,"figs/MultiQA/");
    }

    return vReturn;
}
TH1* LoadVzQAHist(TString inputnames){
    TFile* resultfile =
        TFile::Open("AnalysisResults_Xi1530QA_MB.root","UPDATE");
    TH1* vReturn;
    
    
    TFile* tempfile = new TFile(Form(
    "traindata/%s.root",inputnames.Data()));
    cout << Form(
    "traindata/%s.root",inputnames.Data()) << endl;
    auto templist = (TList*)tempfile->Get("RsnOut_Sigma1385MB");
    auto temphist = (TH2*)templist->FindObject("hVzVsCent");
    auto temphist_prj = (TH1*)temphist->ProjectionY();
    temphist_prj->SetTitle(Form("hVzVsCent_%s",inputnames.Data()));
    vReturn = (TH1*)temphist_prj->Clone();  // signal
    resultfile->cd();
    vReturn->Write(Form(
    "hVzVsCent_%s",inputnames.Data()));
    resultfile->Close();

    return vReturn;
}
TH1* LoadV0MQAHist(TString inputnames){
    TFile* resultfile =
        TFile::Open("AnalysisResults_Xi1530QA_MB.root","UPDATE");
    TH1* vReturn;
    
    
    TFile* tempfile = new TFile(Form(
    "traindata/%s.root",inputnames.Data()));
    cout << Form(
    "traindata/%s.root",inputnames.Data()) << endl;
    

    auto templist = (TList*)tempfile->Get("RsnOut_Sigma1385MB");
    auto temphist = (TH2*)templist->FindObject("hVzVsCent");
    auto temphist_prj = (TH1*)temphist->ProjectionY();
    temphist_prj->SetTitle(Form("hVzVsCent_%s",inputnames.Data()));
    vReturn = (TH1*)temphist_prj->Clone();  // signal
    resultfile->cd();
    vReturn->Write(Form(
    "hVzVsCent_%s",inputnames.Data()));
    resultfile->Close();

    return vReturn;
}
TH1* LoadXiRawQAHist(TString inputnames){
    
    TFile* tempfile = new TFile(Form(
    "AnalysisResults_Extracted_Xi_Raw_%s.root",inputnames.Data()));

    auto temphist = (TH1*)tempfile->Get("hDataRawYield_norm");
    return temphist;
}
TH1* LoadXiQAMeanHist(TString inputnames){
    
    TFile* tempfile = new TFile(Form(
    "AnalysisResults_Extracted_Xi_Raw_%s.root",inputnames.Data()));

    auto temphist = (TH1*)tempfile->Get("hFitMean");
    return temphist;
}
TH1* LoadXiQASigmaHist(TString inputnames){
    
    TFile* tempfile = new TFile(Form(
    "AnalysisResults_Extracted_Xi_Raw_%s.root",inputnames.Data()));

    auto temphist = (TH1*)tempfile->Get("hFitSigma");
    return temphist;
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