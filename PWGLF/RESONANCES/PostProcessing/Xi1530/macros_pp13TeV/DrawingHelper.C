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
    auto base = (TH1*)inputfile->Get("hXiSpectrum");
    gROOT->cd();
    TH1D* hReturn = (TH1D*)base->Clone();
    inputfile->Close();
    return hReturn;
}
vector<double> GetdNdetawithError(double multi_start, double multi_end){
    // Return dN/deta with give Multiplicity bin.
    // it works with only dedicated multiplicit bins(see below)
    // return {value, err}

    vector<double> returnarray;

    //--dNdeta histogram
    // Ref: https://arxiv.org/pdf/1908.01861.pdf
    // LHC16k data.
    // Error was asym error, so choosed bigger error for sym error.

    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta = 
    {0, 26.32, 19.51, 15.45, 13.14, 11.63, 9.50, 7.68, 6.35, 4.36, 2.67};
    vector<double> dNchdeta_e = 
    {0,  0.40,  0.29,  0.23,  0.20,  0.17,  0.14, 0.11, 0.10, 0.06, 0.04};
    

    // input must be in the multiplicity range
    if( multi_end > 0.2){
        if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_start) == end(dNchdeta_multibin))
            return {99,99};
        if(std::find(dNchdeta_multibin.begin(), dNchdeta_multibin.end(), multi_end) == end(dNchdeta_multibin))
            return {99,99};
    }

    // special cases
    if((multi_start == 0) && (multi_end == 0.01)){
        returnarray = {35.92, 0.78};
    }
    else if((multi_start == 0.01) && (multi_end == 0.05)){
        returnarray = {32.19, 0.56};
    } else if ((multi_start == 0.05) && (multi_end == 0.1)) {
        returnarray = {30.13, 0.49};
    } else if ((multi_start == 0.01) && (multi_end == 0.1)) {
        returnarray = {30.89, 0.57};
    } else if ((multi_start == 0) && (multi_end == 5)) {
        returnarray = {21.20, 0.28};
    }
    else if((multi_start == 0) && (multi_end == 0)){
        returnarray = {5.31, 0.18}; // INEL // https://doi.org/10.1016/j.physletb.2015.12.030
    }
    else if((multi_start == 0) && (multi_end == 100)){
        returnarray = {6.94, 0.10}; // INEL>0
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