////////////////////////////////////////////////////////////
//
// Macro for determining which raw V0M value corresponds
// to a certain high-multiplicity-percentile of a min-bias
// data sample for anchoring high-multiplicity-triggered
// data
//
////////////////////////////////////////////////////////////

// global pointer to calib histogram
TH1D *hCalib = 0x0;

void DetermineAnchorsPP(const Char_t* inputDir, TString lPeriodName = "LHC18f", Int_t runNo, const Char_t* chunkName = "", Bool_t automaticMode=kFALSE) {
   //
   // In automatic mode the function does not request any standard input and creates anchor points if the fit was ok
   // using the anchor value determined automatically.
   // One must check the QA plots for individual run to make sure the automatic values are fine and eventually run this 
   // function again in manual mode.
   //
   
    Bool_t lUseDefaultAnchorPercentile = kFALSE;
    Double_t  lDefaultAnchorPercentile = 0.10;
    Double_t  lMinimumAnchorPercentile = 0.05;

    // open minimum bias OADB file
    TString lOADBfile = Form("OADB-%s-MB.root", lPeriodName.Data());

    cout << "Opening minimum bias info ... " << endl;
    TFile *foadb = new TFile( lOADBfile.Data(), "READ" );
    AliOADBContainer *lOADBcontainer = (AliOADBContainer*)foadb->Get("MultSel");

    // set percentile boundaries for the estimator histos
    // (based on what is implemented in the calibration)
    Double_t lDesiredBoundaries[1000];
    Long_t   lNDesiredBoundaries=0;
    lDesiredBoundaries[0] = 0.0;
    //From High To Low Multiplicity
    for( Int_t ib = 1; ib < 101; ib++) { // 100 bins  ] 0.0 , 0.1 ]
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.01;
    }
    for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 1.0 , 10. ]
        lNDesiredBoundaries++;
        lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 0.1;
    }
    for( Int_t ib = 1; ib < 91; ib++) { // 90 bins ] 10.0 , 100.0 ]
        lNDesiredBoundaries++;
        lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] + 1.0;
    }
    
    FILE *fap = 0x0;

    // auxiliary objects
    TLegend *legEstimator = 0x0;
    //
    TLine *anchorLine = new TLine();
    anchorLine->SetLineStyle(2);
    //
    TLatex *latex = new TLatex();
    latex->SetTextFont(42);
    latex->SetTextSize(0.025);

    // constant function for the scaling factor determination
    TF1 *fpol0 = new TF1("fpol0", "[0]", 0.005, lMinimumAnchorPercentile);
    fpol0->SetLineStyle(3);
    fpol0->SetLineWidth(1);
    fpol0->SetLineColor(kBlack);
    TF1 *fpol0_hi = (TF1*)fpol0->Clone("fpol0_hi");
    TF1 *fpol0_lo = (TF1*)fpol0->Clone("fpol0_lo");
    //
    Int_t npar = 3;
    TF1 *fturnon = new TF1("fturnon", func_turnon, 0., 1., npar);
    fturnon->SetParameters(1., 0.1, -1.);
    fturnon->SetParLimits(1, lMinimumAnchorPercentile, 1.0);
    fturnon->SetParLimits(2, -1.e15, 0.);
    fturnon->SetLineColor(1);

    // open input AnalysisResults.root file for the VHM sample
    TString fileIdentifier = Form("%d", runNo);
    if(chunkName[0]!='\0') fileIdentifier = chunkName;
    TFile *fin = TFile::Open(Form("%s/AnalysisResults_%s.root", inputDir, fileIdentifier.Data()), "READ");
    TTree *treeEvent = (TTree*)fin->Get("MultSelection/fTreeEvent");

    cout << "   - run number....................: " << runNo << endl;

    // define estimator histo for this run
    TH1D* hEstimator = new TH1D(Form("hEstimator_%d", runNo), "", lNDesiredBoundaries, lDesiredBoundaries);
    hEstimator->Sumw2();
    hEstimator->GetXaxis()->SetTitle("V0M Percentile");
    hEstimator->GetYaxis()->SetTitle("Counts");
    hEstimator->SetStats(0);
    hEstimator->SetLineColor(kRed);

    // get corresponding calibration histogram from OADB
    AliOADBMultSelection* lOADB = (AliOADBMultSelection*)lOADBcontainer->GetObject( runNo, "Default" );
    if( (Int_t)lOADBcontainer->GetIndexForRun( runNo )<0 ) {
      cout << "   ---> Warning: no calibration histo found for this run - skipping..." << endl;
      return;
    }

    // set the pointer to the calib histo for this run
    hCalib = (TH1D*)lOADB->GetCalibHisto( "hCalib_V0M" );;
    //
    Double_t nall = treeEvent->Draw(Form("get_percentile(fAmplitude_V0A+fAmplitude_V0C)>>hEstimator_%d", runNo), 
                                        Form("fRunNumber==%d && fEvSel_Triggered && fEvSel_IsNotPileupInMultBins && fEvSel_PassesTrackletVsCluster && fEvSel_INELgtZERO && fEvSel_HasNoInconsistentVertices && TMath::Abs(fEvSel_VtxZ)<=10.0 && isSelectedHM(fEvSel_TriggerMask)", runNo),
                                        "goff"); 

    hEstimator->Scale(1., "width");
    Double_t nevents = (Double_t)hEstimator->GetEntries();
    cout << "   - number of events (selected)...: " << nevents << endl;

    // draw histogram
    TCanvas *cEstimator = new TCanvas(Form("cEstimator_%d", runNo), "Estimator Distribution", 10, 10, 1000, 750);
    cEstimator->SetRightMargin(0.05);
    cEstimator->SetTopMargin(0.11);

    hEstimator->GetXaxis()->SetRangeUser(0., 0.2);
    hEstimator->Draw("hist e0");
    latex->SetNDC();
    latex->SetTextSize(0.06);
    latex->DrawLatex(0.1, 0.93, Form("Run: %d", runNo));

    // first, fit a pol0 in the flat region (usually up to 0.05)
    hEstimator->Fit(fpol0, "RQ0");
    Double_t flat_top = fpol0->GetParameter(0);

    // get standard deviantion of bin contents in the flat region
    Double_t flat_top_stdev = 0.;
    for(Int_t ibin=1; ibin<=hEstimator->FindBin(lMinimumAnchorPercentile); ++ibin) {
      Double_t content = hEstimator->GetBinContent(ibin);
      Double_t   width = hEstimator->GetBinWidth(ibin);
      flat_top_stdev += TMath::Power((content-flat_top), 2.)*width;
    }
    flat_top_stdev = TMath::Sqrt(flat_top_stdev/lMinimumAnchorPercentile) / 2.;
    fpol0_hi->SetParameter(0, flat_top+flat_top_stdev);
    fpol0_lo->SetParameter(0, flat_top-flat_top_stdev);
        
    // now, fix the constant parameter in the turnon function
    fturnon->SetParameters(1., 0.1, -1.);
    fturnon->FixParameter(0, flat_top);
        
    // get the maximum range to perform the fit
    Double_t range_max = (hEstimator->GetBinLowEdge(hEstimator->FindLastBinAbove())) / 1.8;
    fturnon->SetRange(0.005, (range_max>0.1) ? range_max : 0.1);

    // get anchor percentile
    Double_t anchor_percentile = -1.;
    TString fitstatus = "";
    if(nevents>0) {
      TFitResultPtr fitr = hEstimator->Fit(fturnon, "RQM");
      fturnon->Draw("lsame");
      fitstatus = gMinuit->fCstatu;
    }
    cEstimator->Flush();
    cEstimator->Update();
    cout << "   - fit status....................: " << fitstatus << endl;
    if( !fitstatus.Contains("OK") ) {
      if(gROOT->IsBatch()) {
         cout << "   ---> Warning: fit failed! -- skipping this run..." << endl;
         if(!automaticMode) {
            fap = fopen(Form("temp/anchors/Anchor_%s_%d_VHM.txt", lPeriodName.Data(), runNo), "w");
            fprintf(fap, "%d %d %.2lf %lf\n", runNo, runNo, -1., -1.);
         }
         return;
      }
      
      if(!automaticMode) {
         cout << "   - Please, provide an anchor percentile to continue: " << endl;
         cout << "     (entering a negative value will skip this run)" << endl;
         cout << "     >>>> anchor percentile: "; 
         cin >> anchor_percentile;
         if(anchor_percentile<0.) {
            cout << "   ---> Warning: percentile provided is negative -- skipping this run..." << endl;
            fap = fopen(Form("temp/anchors/Anchor_%s_%d_VHM.txt", lPeriodName.Data(), runNo), "w");
            fprintf(fap, "%d %d %.2lf %lf\n", runNo, runNo, -1., -1.);
            return;
         }
      }
      else return;      // in automatic mode we do not create an anchor file
    }
    else {
      if(lUseDefaultAnchorPercentile) anchor_percentile = lDefaultAnchorPercentile;
      else {
         anchor_percentile = fturnon->GetParameter(1);
         if( anchor_percentile >= lDefaultAnchorPercentile ) anchor_percentile = lDefaultAnchorPercentile;
         else {
            Double_t this_percentile = lDefaultAnchorPercentile;
            while(this_percentile>=lMinimumAnchorPercentile) {
               Double_t diff = flat_top - fturnon->Eval(this_percentile);
               if( diff < (flat_top_stdev) ) { 
                  anchor_percentile = this_percentile; 
                  break; 
               }
               else this_percentile -= 0.01; 
            }
         }
       }
    }
    //fturnon->DrawCopy("lsame");
    fpol0_hi->SetRange(0.005, range_max);
    fpol0_hi->DrawCopy("lsame");
    fpol0_lo->SetRange(0.005, range_max);
    fpol0_lo->DrawCopy("lsame");
    cout << "   - anchor percentile found.......: " << anchor_percentile << endl;
    Double_t y_max = fturnon->Eval(anchor_percentile);
    anchorLine->DrawLine(anchor_percentile, 0., anchor_percentile, y_max);

    // find corresponding anchor point
    Double_t anchor_point = -1.; 
    for(Int_t i=0; i<hCalib->GetNbinsX(); ++i) {
       if(hCalib->GetBinContent( i+1 ) < anchor_percentile ) {
          anchor_point = hCalib->GetBinLowEdge( i+1 );
          break;
       }
    }
    cout << "   - corresponding anchor point....: " << anchor_point << endl;

    // please, confirm
    cEstimator->Update();
    if(!automaticMode) {
      Bool_t isOk = kFALSE;
      cout << "   ---> Is this Ok? (0->no; 1->yes): ";
      cin >> isOk;
      if(!isOk) {
         cout << "   - Please, provide an anchor percentile to continue: " << endl;
         cout << "     (entering a negative value will skip this run)" << endl;
         cout << "     >>>> anchor percentile: "; 
         cin >> anchor_percentile;
         anchorLine->SetLineStyle(1);
         anchorLine->DrawLine(anchor_percentile, 0., anchor_percentile, y_max);
         cEstimator->Update();
      }

      if(anchor_percentile<lMinimumAnchorPercentile) {
         cout << "   ---> Warning: anchor percentile too low -- skipping this run..." << endl;
         fap = fopen(Form("temp/anchors/Anchor_%s_%d_VHM.txt", lPeriodName.Data(), runNo), "w");
         fprintf(fap, "%d %d %.2lf %lf\n", runNo, runNo, -1., -1.);
         return;
      }
    }
        
    // find corresponding anchor point
    Double_t anchor_point = -1.; 
    for(Int_t i=0; i<hCalib->GetNbinsX(); ++i) {
       if(hCalib->GetBinContent( i+1 ) < anchor_percentile ) {
         anchor_point = hCalib->GetBinLowEdge( i+1 );
         break;
       }
    }
    cout << "   - confirmed anchor point........: " << anchor_point << endl;

    latex->SetTextAlign(31);
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.88, 0.85, Form("%.3g events", nevents));
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.88, 0.96, "anchor percentile found: ");
    latex->DrawLatex(0.88, 0.92, "corresponding anchor point: ");
    latex->SetTextAlign(11);
    latex->DrawLatex(0.88, 0.96,  Form("%.2lf%%", anchor_percentile));
    latex->DrawLatex(0.88, 0.92, Form("%.1lf", anchor_point));

    // write anchors for this run to a txt file
    fap = fopen(Form("temp/anchors/Anchor_%s_%d_VHM.txt", lPeriodName.Data(), runNo), "w");
    fprintf(fap, "%d %d %.2lf %lf\n", runNo, runNo, anchor_percentile, anchor_point);

    // save output
    cEstimator->Update();
    cEstimator->SaveAs(Form("temp/anchors/Anchor_%s_%d_VHM.png", lPeriodName.Data(), runNo));
        
    fclose(fap);
}

// function to get V0M percentile
Double_t get_percentile(Double_t amplitude) 
{
    
    Double_t percentile = hCalib->GetBinContent( hCalib->FindBin( amplitude ) );

    return percentile;

}


// function to fit high-mult. trigger turn-on curve
Double_t func_turnon(Double_t *x, Double_t *par)
{

    Double_t f = par[0];
    if(x[0]>par[1]) f = par[0] + par[2]*TMath::Power(x[0]-par[1], 2.);

    return f;

}

Bool_t isSelectedHM(AliVEvent::EOfflineTriggerTypes trgMask) {
   
   Bool_t isSel = kFALSE;
   isSel = trgMask & AliVEvent::kHighMultV0;
   
   return isSel;
}
