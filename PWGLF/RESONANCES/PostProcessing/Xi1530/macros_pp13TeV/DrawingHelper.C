void SaveCanvas(TCanvas* c,
                char const* name = "temp",
                TString path = "figs/",
                char const* type = "pdf") {
    // Save canvas with path. if path is not there, make a folder.
    if (gSystem->Exec(Form("ls %s", path.Data())) == 256)
        gSystem->Exec(Form("mkdir -p %s", path.Data()));
    c->SaveAs(Form("%s%s.%s", path.Data(), name, type));
}
void SavePad(TPad* p,
             char const* name = "temp",
             TString path = "figs/",
             char const* type = "pdf") {
    // Save Pad with path.
    TCanvas* ctemp = new TCanvas();
    TPad* clone = (TPad*)p->DrawClone();
    clone->SetPad(0, 0, 1, 1);
    SaveCanvas(ctemp, name, path, type);
}
TH1* MakeHistfromArray(char const* name,
                        vector<double> dArray,
                        vector<double> ptbin,
                        vector<double> eArray = {}) {
    double* ptbin_array = &ptbin[0];
    bool IsError = eArray.size() > 1;
    TH1* htemp = new TH1D(Form("%s", name), "", ptbin.size() - 1, ptbin_array);
    for (int i = 0; i < dArray.size(); i++) {
        htemp->SetBinContent(i + 1, dArray[i]);
        if (IsError)
            htemp->SetBinError(i + 1, eArray[i]);
    }
    return htemp;
}
int GetSerialColors(int colorbins = 1){
    const int NRGBs = 5;
    double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    double red[NRGBs] = {0.00, 0.00, 0.87, 0.9 * 1.00, 0.51};
    double green[NRGBs] = {0.00, 0.81, 0.9 * 1.00, 0.20, 0.00};
    double blue[NRGBs] = {0.51, 0.9 * 1.00, 0.12, 0.00, 0.00};
    int FIf = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                                 colorbins);
    return FIf;
}
vector<double> GetBinError(TH1* hinput) {
    // return vector array of bin errors
    vector<double> templist;
    for (int bin = 0; bin < hinput->GetNbinsX(); bin++) {
        double tempe = hinput->GetBinError(bin + 1);
        templist.push_back(tempe);
    }
    return templist;
}
vector<double> GetBinContents(TH1* hinput) {
    // return vector array of bin contents
    vector<double> templist;
    for (int bin = 0; bin < hinput->GetNbinsX(); bin++) {
        double tempv = hinput->GetBinContent(bin + 1);
        templist.push_back(tempv);
    }
    return templist;
}
vector<double> GetBinErrorFrations(TH1* hinput) {
    // return vector array of bin error/content
    vector<double> templist;
    vector<double> tempy = GetBinContents(hinput);
    vector<double> tempe = GetBinError(hinput);
    for (int bin = 0; bin < tempy.size(); bin++) {
        double tempr = tempe[bin] / tempy[bin];
        templist.push_back(tempr);
    }
    return templist;
}
TH1* GetSpectrafromName(TString inputfilename) {
    TFile* inputfile = new TFile(inputfilename.Data());
    TH1* base = (TH1*)inputfile->Get("hXiSpectrum");

    return base;
}