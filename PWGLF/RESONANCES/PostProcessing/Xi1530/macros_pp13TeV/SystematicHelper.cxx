// Systematic Helper Class - Bong-Hwi Lim(2019)
//
// This class will take:
//   - 1 Default histogram (B)
//   - 1 or more variation histograms (A)
//
// This calls will make:
//   - Ratio between default histogram and variations
//     - Ratio(A/B)
//     - Differance(|A-B|)
//     - Diffrerance Ratio((A-B)/B)
//     - Abas Differance Ratio (|A-B|/B)
//   - Standard deviation of Diffrerance Ratio
//   - QA Canvas
//
// How to use:
/*
 #include "SystematicHelper.cxx" // don't forget to include it.

 auto Systematic_SigExt = SystematicHelper(hDefault); // add Default histogram
 Systematic_SigExt.AddVariationtHistograms(hSigExtSys_variations_type1);
                                                 // add variations.
                                                 // it can be a TH1* or
                                                 // vector<TH1*>
 Systematic_SigExt.SetBarlowCheck();             // turn on the barlow check
                                                 // for all pT(X) bins
 // Systematic_SigExt.SetVerbose();              // for debuging
 Systematic_SigExt.SetDefaultYRange({0.5, 1.5}); // Adjust output histogram
 Systematic_SigExt.SetYaxisTitle("Ratio");       //
 Systematic_SigExt.InitAbsDiffRatioColors(false);

 auto hSigExtErrors = Systematic_SigExt.GetAbsDiffRatio();
                                                 //will return vector<TH1*>
 Systematic_SigExt.SetVariationNames(SigExtSys_type1bins);
                                                 // Input vector<TString> names
 Systematic_SigExt.InitVariationColors();        // Initialize the color of
                                                 //outputs
 auto Systematic_SigExt.GetQAPlot("AbsDiffRatio"); // output canvas
*/
class SystematicHelper {
   public:
    // Constructors
    SystematicHelper(TH1* h) { SetDefaultHistogram(h); }
    SystematicHelper(TH1* h, TH1* hinput) {
        SetDefaultHistogram(h);
        AddVariationtHistograms(hinput);
    }
    SystematicHelper(TH1* h, vector<TH1*> hList) {
        SetDefaultHistogram(h);
        AddVariationtHistograms(hList);
    }

    //// Setters
    // default histogram
    void SetDefaultHistogram(TH1* h) { hDefault = h; };
    void SetDefaultHistogram(TObject* o) { hDefault = (TH1*)o; };

    // Variation histograms
    void AddVariationtHistograms(TH1* h) { hVariations.push_back((TH1*)h); };
    void AddVariationtHistograms(TObject* o) {
        hVariations.push_back((TH1*)o);
    };
    void AddVariationtHistograms(vector<TH1*> hList) {
        for (auto const& hInput : hList)
            hVariations.push_back((TH1*)hInput);
    }
    
    // Default Setter
    void SetVariationNames(vector<TString> inputnamelist){
        hVariationNames = inputnamelist;};
    void SetBarlowCheck(bool setbarlow = true) { fSetBarlow = setbarlow; };
    void SetVerbose(bool setverbose = true) { fSetVerbose = setverbose; };

    // Default axis range
    void SetDefaultYRange(vector<double> iRange) {
        SetDefaultYmax(iRange[1]);
        SetDefaultYmin(iRange[0]);
    };
    void SetDefaultYmax(double imax) { fYmax = imax; };
    void SetDefaultYmin(double imin) { fYmin = imin; };
    void SetDefaultXRange(vector<double> iRange) {
        SetDefaultXmax(iRange[1]);
        SetDefaultXmin(iRange[0]);
    };
    void SetDefaultXmax(double imax) { fXmax = imax; };
    void SetDefaultXmin(double imin) { fXmin = imin; };

    // Default axis title
    void SetYaxisTitle(TString iTitle) { fDefaultYtitle = iTitle; };
    void SetXaxisTitle(TString iTitle) { fDefaultXtitle = iTitle; };

    // Default colors
    void InitColors() { InitColors(hVariations.size()); }
    void InitColors(int nbins) {
        FIh =
            TColor::CreateGradientColorTable(5, stops, red, green, blue, nbins);
    };

    // Functions
    vector<TH1*> GetRatio() {
        if (hRatios.size() < 1) {
            CalculateRatio();
            return hRatios;
        } else
            return hRatios;
    }
    vector<TH1*> GetDiffRatio() {
        if (hDiffRatios.size() < 1) {
            CalculateDiffRatio();
            return hDiffRatios;
        } else
            return hDiffRatios;
    }
    vector<TH1*> GetDifference() {
        if (hDifferences.size() < 1) {
            CalculateDifference();
            return hDifferences;
        } else
            return hDifferences;
    }
    vector<TH1*> GetAbsDiffRatio() {
        if (hAbsDiffRatios.size() < 1) {
            CalculateAbsDiffRatio();
            return hAbsDiffRatios;
        } else
            return hAbsDiffRatios;
    }
    TH1* GetStdevRatio(bool AddMean = false) {
        auto temp = (TH1*)hDefault->Clone();
        if (hRatios.size() < 1) {
            if (fSetVerbose)
                cout << "Get hRatios first";
            CalculateRatio();
        }
        for (int bin = 1; bin < temp->GetNbinsX(); bin++)
            temp->SetBinContent(bin, 0);

        for (int bin = 1; bin < temp->GetNbinsX(); bin++) {
            TH1* tempStdev = new TH1F(Form("hMean_%d",bin), "", 200, -1, 1);
            for (auto const& hInput : hRatios)
                tempStdev->Fill(hInput->GetBinContent(bin)-1);
            if (fSetVerbose)
                cout << "bin: " << bin << ", mean: " << tempStdev->GetMean()
                     << ", Std: " << tempStdev->GetStdDev() << endl;
            double tempcontent = tempStdev->GetStdDev();
            if(AddMean) 
                tempcontent += TMath::Abs(tempStdev->GetMean());
            temp->SetBinContent(bin, tempcontent);
        }
        temp->SetMaximum(fYmax - 1);
        temp->SetMinimum(0);
        temp->GetXaxis()->SetRangeUser(fXmin, fXmax);
        return temp;
    }
    TH1* GetMaxAbsDiffRatio() {
        auto temp = (TH1*) hDefault->Clone();
        if (hAbsDiffRatios.size() < 1) {
            if (fSetVerbose)
                cout << "Get AbsDiffRatio first";
            CalculateAbsDiffRatio();
        }
        temp->SetMaximum(1);
        for (int bin = 1; bin < temp->GetNbinsX(); bin++)
            temp->SetBinContent(bin, 0);
        
        for (auto const& hInput : hAbsDiffRatios) {
            for (int bin = 1; bin < temp->GetNbinsX(); bin++){
                if (hInput->GetBinContent(bin) > temp->GetBinContent(bin))
                    temp->SetBinContent(bin, hInput->GetBinContent(bin));
                if (fSetVerbose)
                    cout << "bin: " << bin
                         << ", content: " << temp->GetBinContent(bin) << endl;
            }
        }
        temp->SetMaximum(fYmax - 1);
        temp->SetMinimum(0);
        temp->GetXaxis()->SetRangeUser(fXmin, fXmax);
        return temp;
    }
    TH1* GetSumAbsDiffRatio() {
        auto temp = (TH1*)hDefault->Clone();
        if (hAbsDiffRatios.size() < 1) {
            if (fSetVerbose)
                cout << "Get AbsDiffRatio first";
            CalculateAbsDiffRatio();
        }
        temp->SetMaximum(1);
        for (int bin = 1; bin < temp->GetNbinsX(); bin++)
            temp->SetBinContent(bin, 0);

        for (auto const& hInput : hAbsDiffRatios) {
            for (int bin = 1; bin < temp->GetNbinsX(); bin++) {
                temp->SetBinContent(
                    bin, temp->GetBinContent(bin) + TMath::Power(hInput->GetBinContent(bin),2));
            }
        }
        for (int bin = 1; bin < temp->GetNbinsX(); bin++)
            temp->SetBinContent(bin, TMath::Sqrt(temp->GetBinContent(bin)));
        temp->SetMaximum(fYmax-1);
        temp->SetMinimum(0);
        temp->GetXaxis()->SetRangeUser(fXmin, fXmax);
        return temp;
    }
    void CalculateRatio() {
        if (!hDefault) {
            if (fSetVerbose)
                cout << "No base histogram" << endl;
            return;
        }
        if (hVariations.size() < 1) {
            if (fSetVerbose)
                cout << "No variation histograms" << endl;
            return;
        }
        if (hRatios.size() > 1) {
            if (fSetVerbose)
                cout << "clear old values" << endl;
            hRatios.clear();
        }
        for (auto const& hInput : hVariations) {
            TH1* temp = (TH1*)hInput->Clone();
            TH1* hDivided = (TH1*)hDefault->Clone();
            hDivided->Divide(temp);
            hDivided->SetMaximum(fYmax);
            hDivided->SetMinimum(fYmin);
            hDivided->GetXaxis()->SetRangeUser(fXmin, fXmax);
            hDivided->SetTitle("");
            hDivided->SetLineWidth(2);
            hDivided->GetXaxis()->SetTitle(fDefaultXtitle.Data());
            hDivided->GetYaxis()->SetTitle(fDefaultYtitle.Data());
            hRatios.push_back(hDivided);
        }
    }
    void CalculateDiffRatio() {
        if (!hDefault) {
            if (fSetVerbose)
                cout << "No base histogram" << endl;
            return;
        }
        if (hVariations.size() < 1) {
            if (fSetVerbose)
                cout << "No variation histograms" << endl;
            return;
        }
        if (hDiffRatios.size() > 1) {
            if (fSetVerbose)
                cout << "clear old values" << endl;
            hDiffRatios.clear();
        }
        for (auto const& hInput : hVariations) {
            TH1* temp = (TH1*)hInput->Clone();
            TH1* hDivided = (TH1*)hDefault->Clone();
            hDivided->Divide(temp);
            hDivided->SetMaximum(fYmax - 1);
            hDivided->SetMinimum(fYmin - 1);
            hDivided->GetXaxis()->SetRangeUser(fXmin, fXmax);
            hDivided->SetTitle("");
            hDivided->SetLineWidth(2);
            hDivided->GetXaxis()->SetTitle(fDefaultXtitle.Data());
            hDivided->GetYaxis()->SetTitle(fDefaultYtitle.Data());
            for (int bin = 1; bin < hDivided->GetNbinsX(); bin++) {
                hDivided->SetBinContent(bin, hDivided->GetBinContent(bin) - 1);
            }
            hDiffRatios.push_back(hDivided);
        }
    }
    void CalculateAbsDiffRatio() {
        if (!hDefault) {
            if (fSetVerbose)
                cout << "No base histogram" << endl;
            return;
        }
        if (hVariations.size() < 1) {
            if (fSetVerbose)
                cout << "No variation histograms" << endl;
            return;
        }
        if (hAbsDiffRatios.size() > 1) {
            if (fSetVerbose)
                cout << "clear old values";
            hAbsDiffRatios.clear();
        }
        for (auto const& hInput : hVariations) {
            TH1* temp = (TH1*)hInput->Clone();
            TH1* hDivided = (TH1*)hDefault->Clone();
            hDivided->Divide(temp);
            hDivided->SetMaximum(fYmax - 1);
            hDivided->SetMinimum(0);
            hDivided->GetXaxis()->SetRangeUser(fXmin, fXmax);
            hDivided->SetTitle("");
            hDivided->SetLineWidth(2);
            hDivided->GetXaxis()->SetTitle(fDefaultXtitle.Data());
            hDivided->GetYaxis()->SetTitle(fDefaultYtitle.Data());
            for (int bin = 1; bin < hDivided->GetNbinsX(); bin++) {
                if (fSetBarlow) {
                    double error =  // abs(Variation - Default) -> "Error"
                        abs(temp->GetBinContent(bin) -
                            hDefault->GetBinContent(bin));
                    double errorFraction =  // abs(Variation - Default)/Deafult
                                            // - >"Error fraction"
                        error / temp->GetBinContent(bin);
                    double
                        deltasigma =  // sqrt(abs(variation_stat.error^2 -
                                      // default_stat.error^2)) -> "delta sigma"
                        sqrt(abs(pow(temp->GetBinError(bin), 2) -
                                 pow(hDefault->GetBinError(bin), 2)));
                    double deltasigmafraction =
                        deltasigma / temp->GetBinContent(bin);
                    if (fSetVerbose)
                        cout
                            << "Variation value: " << temp->GetBinContent(bin)
                            << ", Default value: "
                            << hDefault->GetBinContent(bin)
                            << ", Variation error: " << temp->GetBinError(bin)
                            << ", Default error: " << hDefault->GetBinError(bin)
                            << endl;
                    if (errorFraction / deltasigmafraction > 1) {
                        if (fSetVerbose)
                            cout
                                << "-> PASS! errorfraction: " << errorFraction
                                << ", deltasigmafraction:" << deltasigmafraction
                                << ", check: "
                                << errorFraction / deltasigmafraction << endl;
                        hDivided->SetBinContent(
                            bin, TMath::Abs(hDivided->GetBinContent(bin) - 1));
                        }
                        else
                        {
                            hDivided->SetBinContent(bin,0);
                        }
                        
                    } else
                        hDivided->SetBinContent(
                            bin, TMath::Abs(hDivided->GetBinContent(bin) - 1));
                }
                hAbsDiffRatios.push_back(hDivided);
            }
        }
        void CalculateDifference() {
            if (!hDefault) {
                if (fSetVerbose)
                    cout << "No base histogram" << endl;
                return;
            }
            if (hVariations.size() < 1) {
                if (fSetVerbose)
                    cout << "No variation histograms" << endl;
                return;
            }
            if (hDifferences.size() > 1) {
                if (fSetVerbose)
                    cout << "clear old values" << endl;
                hDifferences.clear();
            }
            for (auto const& hInput : hVariations) {
                TH1* temp = (TH1*)hInput->Clone();
                TH1* hSubtracted = (TH1*)hDefault->Clone();
                hSubtracted->Add(temp, -1);
                hSubtracted->GetXaxis()->SetRangeUser(fXmin, fXmax);
                hSubtracted->SetTitle("");
                hSubtracted->SetLineWidth(2);
                hSubtracted->GetXaxis()->SetTitle(fDefaultXtitle.Data());
                hDifferences.push_back(hSubtracted);
            }
        }

        // Colors
        void InitAllListColors(bool fillColorChange = false) {
            InitRatioColors(fillColorChange);
            InitDiffRatioColors(fillColorChange);
            InitAbsDiffRatioColors(fillColorChange);
            InitDifferenceColors(fillColorChange);
            InitVariationColors(fillColorChange);
        }
        void InitRatioColors(bool fillColorChange = false) {
            if (hRatios.size() < 1)
                CalculateRatio();
            InitListColors(hRatios, fillColorChange);
        }
        void InitDiffRatioColors(bool fillColorChange = false) {
            if (hDiffRatios.size() < 1)
                CalculateDiffRatio();
            InitListColors(hDiffRatios, fillColorChange);
        }
        void InitAbsDiffRatioColors(bool fillColorChange = false) {
            if (hAbsDiffRatios.size() < 1)
                CalculateAbsDiffRatio();
            InitListColors(hAbsDiffRatios, fillColorChange);
        }
        void InitDifferenceColors(bool fillColorChange = false) {
            if (hDifferences.size() < 1)
                CalculateDifference();
            InitListColors(hDifferences, fillColorChange);
        }
        void InitVariationColors(bool fillColorChange = false) {
            if (hVariations.size() < 1) {
                if (fSetVerbose)
                    cout << "No variation histograms" << endl;
                return;
            }
            InitListColors(hVariations, fillColorChange);
        }
        void InitListColors(vector<TH1*> hList, bool fillColorChange = false) {
            int totalbin = hList.size();
            if (totalbin < 1) {
                if (fSetVerbose)
                    cout << "No entries" << endl;
                return;
            }
            InitColors(totalbin);
            
            int i = 0;
            for (auto const& hInput : hList) {
                hInput->SetMarkerColor(FIh + totalbin - i);
                hInput->SetLineColor(FIh + totalbin -  i);
                if (fillColorChange)
                    hInput->SetFillColor(FIh + totalbin - i);
                i++;
            }
        }

        void PurgeVariations() { hVariations.clear(); }

        TCanvas* GetQAPlot(TString name = "RatioQA") {
            gStyle->SetOptStat(0);
            gStyle->SetOptTitle(0);
            TCanvas* c1 =
                new TCanvas(name.Data(), name.Data(), 0, 0, 850, 1150);
            c1->Range(0, 0, 1, 1);

            //*************************************************
            // Bottom plot
            TPad* c1_1 = new TPad("c1_1", "newpad", 0.01, 0.01, 1, 0.32);
            c1_1->Draw();
            c1_1->cd();
            c1_1->SetTickx(1);
            c1_1->SetTicky(1);
            c1_1->SetTopMargin(0.01);
            c1_1->SetLeftMargin(0.15);
            c1_1->SetBottomMargin(0.2);
            c1_1->SetRightMargin(0.01);
            c1_1->SetFillStyle(0);
            // c1_1->SetLogy(true);
            if (name.Contains("AbsDiffRatio")) {
                hAbsDiffRatios[0]->Draw("HIST");
                hAbsDiffRatios[0]->SetLineWidth(2);
                hAbsDiffRatios[0]->GetXaxis()->SetTitleSize(0.08);
                hAbsDiffRatios[0]->GetXaxis()->SetLabelSize(0.08);
                hAbsDiffRatios[0]->GetYaxis()->SetLabelSize(0.08);
                hAbsDiffRatios[0]->GetYaxis()->SetTitleSize(0.08);
                hAbsDiffRatios[0]->GetYaxis()->SetTitleOffset(0.6);
                for (int i = 1; i < hAbsDiffRatios.size(); i++)
                    hAbsDiffRatios[i]->Draw("HIST same");
            } else if (name.Contains("DiffRatio")) {
                hDiffRatios[0]->Draw("E");
                hDiffRatios[0]->SetLineWidth(2);
                hDiffRatios[0]->GetXaxis()->SetTitleSize(0.08);
                hDiffRatios[0]->GetXaxis()->SetLabelSize(0.08);
                hDiffRatios[0]->GetYaxis()->SetLabelSize(0.08);
                hDiffRatios[0]->GetYaxis()->SetTitleSize(0.08);
                hDiffRatios[0]->GetYaxis()->SetTitleOffset(0.6);
                for (int i = 1; i < hDiffRatios.size(); i++)
                    hDiffRatios[i]->Draw("E same");
                TF1* line1 = new TF1("line1", "0", -1, 15);
                line1->SetLineColor(1);
                line1->SetLineWidth(2);
                line1->SetLineStyle(2);
                line1->Draw("same");
            } else if (name.Contains("Ratio")) {
                hRatios[0]->Draw("E");
                hRatios[0]->SetLineWidth(2);
                hRatios[0]->GetXaxis()->SetTitleSize(0.08);
                hRatios[0]->GetXaxis()->SetLabelSize(0.08);
                hRatios[0]->GetYaxis()->SetLabelSize(0.08);
                hRatios[0]->GetYaxis()->SetTitleSize(0.08);
                hRatios[0]->GetYaxis()->SetTitleOffset(0.6);
                for (int i = 1; i < hRatios.size(); i++)
                    hRatios[i]->Draw("E same");

                TF1* line1 = new TF1("line1", "1", -1, 15);
                line1->SetLineColor(1);
                line1->SetLineWidth(2);
                line1->SetLineStyle(2);
                line1->Draw("same");
            } else if (name.Contains("Difference")) {
                hDifferences[0]->Draw("E");
                hDifferences[0]->SetLineWidth(2);
                hDifferences[0]->GetXaxis()->SetTitleSize(0.08);
                hDifferences[0]->GetXaxis()->SetLabelSize(0.08);
                hDifferences[0]->GetYaxis()->SetLabelSize(0.08);
                hDifferences[0]->GetYaxis()->SetTitleSize(0.08);
                hDifferences[0]->GetYaxis()->SetTitleOffset(0.6);
                for (int i = 1; i < hDifferences.size(); i++)
                    hDifferences[i]->Draw("E same");
            }
            // End bottom plot
            //*************************************************

            //*************************************************
            // Top Plot
            c1->cd();
            TPad* c1_2 = new TPad("c1_2", "newpad", 0.01, 0.32, 1, 1);
            c1_2->SetLogy(true);
            c1_2->SetTickx(1);
            c1_2->SetTicky(1);
            c1_2->Draw();
            c1_2->cd();
            c1_2->SetTopMargin(0.01);
            c1_2->SetLeftMargin(0.15);
            c1_2->SetBottomMargin(0.01);
            c1_2->SetRightMargin(0.01);
            c1_1->SetFillStyle(0);

            hDefault->Draw("E");
            hDefault->SetLabelSize(0.0);
            hDefault->GetXaxis()->SetRangeUser(fXmin, fXmax);
            hDefault->GetXaxis()->SetTitleSize(0.00);
            hDefault->GetYaxis()->SetLabelSize(0.04);
            hDefault->GetYaxis()->SetTitleSize(0.05);
            hDefault->GetYaxis()->SetTitleOffset(1.2);
            for (int i = 0; i < hVariations.size(); i++)
                hVariations[i]->Draw("E same");
            if(hVariationNames.size() > 0){
                auto legendtemp = new TLegend(.18, .03, .58, .27);
                legendtemp->SetBorderSize(0);
                legendtemp->SetFillStyle(0);
                int check = 0;
                for (auto const& hName : hVariationNames) {
                    legendtemp->AddEntry(hVariations[check], hName, "L");
                    check++;
                    if (check + 1  > hVariations.size()) break;
                }
                legendtemp->Draw();
            }
            return c1;
        }

       private:
        TH1* hDefault;
        vector<TH1*> hVariations;

        vector<TH1*> hRatios;
        vector<TH1*> hDiffRatios;
        vector<TH1*> hAbsDiffRatios;
        vector<TH1*> hDifferences;

        vector<TString> hVariationNames;

        bool fSetBarlow = false;
        bool fSetVerbose = false;
        double fXmin = 0.8;
        double fXmax = 8.8;
        double fYmin = 0.5;
        double fYmax = 1.5;
        TString fDefaultXtitle = "#it{p}_{T} (GeV/#it{c})";
        TString fDefaultYtitle = "Ratio to default";

        // Auto Color
        Double_t stops[5] = {0.00, 0.25, 0.5, 0.75, 1.00};
        Double_t red[5] = {0.00, 0.00, 0.87, 1.00, 0.51};
        Double_t green[5] = {0.00, 0.81, 1.00, 0.20, 0.00};
        Double_t blue[5] = {0.51, 1.00, 0.12, 0.00, 0.00};
        Int_t FIh =
            TColor::CreateGradientColorTable(5, stops, red, green, blue, 10);
    };