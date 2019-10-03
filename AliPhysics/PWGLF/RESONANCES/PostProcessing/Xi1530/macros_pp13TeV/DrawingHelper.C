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
    return base;
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
    

    // 2nd Reference (LHC15 study)
    // Ref: https://aliceinfo.cern.ch/Notes/node/510
    /*
    vector<double> dNchdeta_multibin = 
    {0,     1,     5,    10,    15,    20,    30,   40,   50,   70, 100};
    vector<double> dNchdeta = 
    {0, 26.18, 20.16, 16.40, 14.00, 12.28, 10.31, 8.24, 6.62, 4.77, 2.76};
    vector<double> dNchdeta_e = 
    {0,  0.55,  0.41,  0.31,  0.29,  0.25,  0.21, 0.17, 0.13, 0.09, 0.05};
    */

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