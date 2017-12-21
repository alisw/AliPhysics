/*

 >> Macro to calculate the systematics from the fitting procedure
    based on the AliHFCorrFitSystematics class
    
 Developed by S.Bjelogrlic(Utrecht) 
 latest version 23/04/2015
 
 Instructions for use
 For the Setting functions, see instructions
 
 various enums defined in the class
 
 AliHFCorrFitSystematics::SetCombineSystematicsMode(AliHFCorrFitSystematics::kMax); // Set mode to compute the final systematics, kMax (use the maximum variation), kRMS (use the RMS), kSumQuadr (use sum in quadrature)
 
 
 AliHFCorrFitSystematics::AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kTRUE); 
 // modes to use for estimating the systeamtic due to the baseline:
 -kTransverse (use transverse region) - specifiy also the range
 -kMinVar (fits the histogram obtained from the minimum variation due to delta-phi dependendt systematics)
 -kMaxVar (fits the histogram obtained from the maximum variation due to delta-phi dependendt systematics)
 -kLowestPoint (sets the baseline to the lowest point)
 -kBinCount (uses bin counting - yields only)
 -k2PointsAtPiHalf (two points at pi/2)
 -k4PointsAtPiHalf (four points at pi/2)
 After set up
 



*/


TString inputpath = "";
Bool_t useReflectedPlots = kTRUE;
TString outputpath="";
TString strAverage;
Bool_t v2systOn=kFALSE;
Bool_t plotv2Separated=kFALSE;
Double_t v2D=0.,v2had=0.;
Bool_t saveAwaySide=kFALSE;
void SetSaveAwaySidePlots(Bool_t saveAS){
  saveAwaySide=saveAS;
}

void SetV2values(Bool_t usev2=kTRUE,Double_t v2hadrons=0.1,Double_t v2Dmes=0.1, Bool_t plotseparate=kFALSE){
  v2systOn=usev2;
  plotv2Separated=plotseparate;
  v2D=v2Dmes;
  v2had=v2hadrons;
}
void SetInputPath(TString str){
  inputpath=str;
}
void SetOutputPath(TString str){
  outputpath=str;
}
void SetAverageString(TString str){
  strAverage=str;
}

//=========================================================================================================================================================================================================================

void Systematics_pPb_lowpthad(Bool_t useReflected){
    
    TString path = inputpath;
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    AliHFCorrFitSystematics * systfitter = new AliHFCorrFitSystematics();
    
    // set up the sysfitter
    gROOT->ProcessLine(Form(".!mkdir -p %s/FitOutputs_pPb",outputpath.Data()));
    systfitter->CreateOutputFile(Form("%s/FitOutputs_pPb/FitSystematics_pPb_%sAverage_0.3_1.0.root",outputpath.Data(),strAverage.Data())); // set the name of the output file
    systfitter->SetUseCorrelatedSystematics(kTRUE); // kTRUE (will use the Delta phi independent systematics (from the correlation plots) to sum in quadrature for the final systematics)
    systfitter->SetUseCorrelatedSystematicsForWidths(kFALSE); //Delta phi independent systematics for the final systematics, ktrue will use them also for the widths, kfalse not (in my opinion, we should use kfalse)
    systfitter->SetCombineSystematicsMode(AliHFCorrFitSystematics::kRMS); // set method to compute the systematics due to the fixing of the baseline
    systfitter->SetIspPb(kTRUE); // kTrue for pPb, kFALSE for pp
    if(v2systOn){
      systfitter->Setv2ForSystematics(v2had, v2D); // Set v2 values for the systematics (v2had, v2 D)
      systfitter->SetPlotV2SystSeparately(plotv2Separated); // set true if you want the v2 box to be plotted separately
    }

 //   systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    if(useReflected) systfitter->SetFitRange(0,TMath::Pi()); // set the fit range
    if(!useReflected) systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    systfitter->SetReferenceBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi()); // set the default baseline estimation range
    
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi()); //add systematic method (ktrue if is default to compare to, kfalse if is used to estimated the variations)
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.25*TMath::Pi(),0.375*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.375*TMath::Pi(),0.5*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.5*TMath::Pi(),0.625*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.375*TMath::Pi(),0.625*TMath::Pi());
   // systfitter->AddSystematicMode(AliHFCorrFitSystematics::kLowestPoint,kFALSE);
   // systfitter->AddSystematicMode(AliHFCorrFitSystematics::k2PointsAtPiHalf,kFALSE);
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kMinVar,kFALSE);
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kMaxVar,kFALSE);
    //systfitter->AddSystematicMode(AliHFCorrFitSystematics::kBinCount,kFALSE);
    
    //systfitter->CheckBaselineRanges();
    //return;
    systfitter->SetAssociatedTrackPtMin(0.3); // set the min pt assoc (only for plotting purposes)
    systfitter->SetAssociatedTrackPtMax(1); // set the max pt assoc (only for plotting purposes)
    
/*    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
   // ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt3to5assocPt0.3to1.0.root");
    ishistoadded = systfitter->AddHistoToFit(Form("%sCanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt3to5assocPt0.3to1.0.root",path.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.12); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.13); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(3); systfitter->SetDptUpEdge(5); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    */
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
    //ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt0.3to1.0.root");
    ishistoadded = systfitter->AddHistoToFit(Form("%s/CanvaAndVariedHisto%sAverageDzeroDstarDplus_pPb_Pt5to8assocPt0.3to1.0.root",path.Data(),strAverage.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.09); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.12); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(5); systfitter->SetDptUpEdge(8); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
  //  ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt8to16assocPt0.3to1.0.root");
    ishistoadded = systfitter->AddHistoToFit(Form("%s/CanvaAndVariedHisto%sAverageDzeroDstarDplus_pPb_Pt8to16assocPt0.3to1.0.root",path.Data(),strAverage.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.10); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.13); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(8); systfitter->SetDptUpEdge(16); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    
    
    

    
    Bool_t setupfitter = systfitter->SetUpSystematicsFitter(); // set up the fitter (initializes variables)
    if(!setupfitter) {cout << "cannot set up fitter " << endl;
    return;
    }
    Bool_t fitted = systfitter->RunFits(); // main function, does all the calculations
    if(!fitted) {cout << "something went wrong with the fits - sorry " << endl;
        return;
    }
 

    if(saveAwaySide)systfitter->DrawFinalPlots(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE);
    else systfitter->DrawFinalPlots(); // draws the final trend plots


    systfitter->PrintAllSystematicsOnShell(); // print all the systematics on shell
   // systfitter->DrawFinalPlots(kTRUE,kFALSE,kFALSE,kFALSE,kFALSE);
    gROOT->ProcessLine(Form(".!mkdir -p %s/Trends_pPb",outputpath.Data()));
    systfitter->SetOutputDirectory(Form("%s/Trends_pPb",outputpath.Data()));
    systfitter->SaveCanvasesDotC();
    systfitter->SaveCanvasesDotRoot();
    systfitter->SaveCanvasesDotPng();
    systfitter->SaveCanvasesDotEps();
    systfitter->SaveCanvasesDotPdf();
    
  
    if(saveAwaySide)systfitter->DrawFinalCorrelationPlot(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE); // draws the final correlation plots
	else systfitter->DrawFinalCorrelationPlot(); // draws the final correlation plots
   // systfitter->CheckHisto(1);
}



//=========================================================================================================================================================================================================================
//=========================================================================================================================================================================================================================
void Systematics_pPb_highpthad(Bool_t useReflected){
    
    TString path = inputpath;    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    AliHFCorrFitSystematics * systfitter = new AliHFCorrFitSystematics();
    
    // set up the sysfitter
    gROOT->ProcessLine(Form(".!mkdir -p %s/FitOutputs_pPb",outputpath.Data()));
    systfitter->CreateOutputFile(Form("%s/FitOutputs_pPb/FitSystematics_pPb_%sAverage_1.0_99.0.root",outputpath.Data(),strAverage.Data())); // set the name of the output file
    systfitter->SetUseCorrelatedSystematics(kTRUE); // kTRUE (will use the Delta phi independent systematics (from the correlation plots) to sum in quadrature for the final systematics)
    systfitter->SetUseCorrelatedSystematicsForWidths(kFALSE); //Delta phi independent systematics for the final systematics, ktrue will use them also for the widths, kfalse not (in my opinion, we should use kfalse)
    systfitter->SetCombineSystematicsMode(AliHFCorrFitSystematics::kRMS); // set method to compute the systematics due to the fixing of the baseline
    systfitter->SetIspPb(kTRUE); // kTrue for pPb, kFALSE for pp
    if(v2systOn){
      systfitter->Setv2ForSystematics(v2had, v2D); // Set v2 values for the systematics (v2had, v2 D)
      systfitter->SetPlotV2SystSeparately(plotv2Separated); // set true if you want the v2 box to be plotted separately
    }
    //   systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    if(useReflected) systfitter->SetFitRange(0,TMath::Pi()); // set the fit range
    if(!useReflected) systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    systfitter->SetReferenceBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi()); // set the default baseline estimation range
    
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi()); //add systematic method (ktrue if is default to compare to, kfalse if is used to estimated the variations)
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.25*TMath::Pi(),0.375*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.375*TMath::Pi(),0.5*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.5*TMath::Pi(),0.625*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.375*TMath::Pi(),0.625*TMath::Pi());
 //   systfitter->AddSystematicMode(AliHFCorrFitSystematics::kLowestPoint,kFALSE);
 //   systfitter->AddSystematicMode(AliHFCorrFitSystematics::k2PointsAtPiHalf,kFALSE);
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kMinVar,kFALSE);
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kMaxVar,kFALSE);
    //systfitter->AddSystematicMode(AliHFCorrFitSystematics::kBinCount,kFALSE);
    
    //systfitter->CheckBaselineRanges();
    //return;
    systfitter->SetAssociatedTrackPtMin(1); // set the min pt assoc (only for plotting purposes)
    systfitter->SetAssociatedTrackPtMax(99); // set the max pt assoc (only for plotting purposes)
    
 /*   Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
    //ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt3to5assocPt1.0to99.0.root");
    ishistoadded = systfitter->AddHistoToFit(Form("%sCanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt3to5assocPt1.0to99.0.root",path.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.10); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.13); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(3); systfitter->SetDptUpEdge(5); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    */
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
    //ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0to99.0.root");
    ishistoadded = systfitter->AddHistoToFit(Form("%s/CanvaAndVariedHisto%sAverageDzeroDstarDplus_pPb_Pt5to8assocPt1.0to99.0.root",path.Data(),strAverage.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.1); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.11); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(5); systfitter->SetDptUpEdge(8); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
   // ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt8to16assocPt1.0to99.0.root");
    ishistoadded = systfitter->AddHistoToFit(Form("%s/CanvaAndVariedHisto%sAverageDzeroDstarDplus_pPb_Pt8to16assocPt1.0to99.0.root",path.Data(),strAverage.Data()));//
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.1); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.11); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(8); systfitter->SetDptUpEdge(16); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    
    
    
    
    
    Bool_t setupfitter = systfitter->SetUpSystematicsFitter(); // set up the fitter (initializes variables)
    if(!setupfitter) {cout << "cannot set up fitter " << endl;
        return;
    }
    Bool_t fitted = systfitter->RunFits(); // main function, does all the calculations
    if(!fitted) {cout << "something went wrong with the fits - sorry " << endl;
        return;
    }
    
    systfitter->PrintAllSystematicsOnShell(); // print all the systematics on shell
   
    if(saveAwaySide)systfitter->DrawFinalPlots(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE);
    else systfitter->DrawFinalPlots(); // draws the final trend plots
    
    // systfitter->DrawFinalPlots(kTRUE,kFALSE,kFALSE,kFALSE,kFALSE);
    gROOT->ProcessLine(Form(".!mkdir -p %s/Trends_pPb",outputpath.Data()));
    systfitter->SetOutputDirectory(Form("%s/Trends_pPb",outputpath.Data())); 
    systfitter->SaveCanvasesDotC();
    systfitter->SaveCanvasesDotRoot();
    systfitter->SaveCanvasesDotPng();
    systfitter->SaveCanvasesDotEps();
    systfitter->SaveCanvasesDotPdf();
    
    if(saveAwaySide)systfitter->DrawFinalCorrelationPlot(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE); // draws the final correlation plots
    else systfitter->DrawFinalCorrelationPlot(); // draws the final correlation plots

    // systfitter->CheckHisto(1);
}


//=========================================================================================================================================================================================================================
//=========================================================================================================================================================================================================================


void Systematics_pPb_integratedpthad(Bool_t useReflected){
    
    TString path = inputpath;
    
    cout << "Opening files from " << path << endl;
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    AliHFCorrFitSystematics * systfitter = new AliHFCorrFitSystematics();
    cout << "crash 1 " << endl;
     // set up the sysfitter
    gROOT->ProcessLine(Form(".!mkdir -p %s/FitOutputs_pPb",outputpath.Data()));
    systfitter->CreateOutputFile(Form("%s/FitOutputs_pPb/FitSystematics_pPb_%sAverage_0.3_99.0.root",outputpath.Data(),strAverage.Data())); // set the name of the output file
    systfitter->SetUseCorrelatedSystematics(kTRUE); // kTRUE (will use the Delta phi independent systematics (from the correlation plots) to sum in quadrature for the final systematics)
    systfitter->SetUseCorrelatedSystematicsForWidths(kFALSE); //Delta phi independent systematics for the final systematics, ktrue will use them also for the widths, kfalse not (in my opinion, we should use kfalse)
    systfitter->SetCombineSystematicsMode(AliHFCorrFitSystematics::kRMS); // set method to compute the systematics due to the fixing of the baseline
    systfitter->SetIspPb(kTRUE); // kTrue for pPb, kFALSE for pp
    if(v2systOn){
      systfitter->Setv2ForSystematics(v2had, v2D); // Set v2 values for the systematics (v2had, v2 D)
      systfitter->SetPlotV2SystSeparately(plotv2Separated); // set true if you want the v2 box to be plotted separately
    }
    //   systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    if(useReflected) systfitter->SetFitRange(0,TMath::Pi()); // set the fit range
    if(!useReflected) systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    cout << "crash 2 " << endl;
    systfitter->SetReferenceBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi()); // set the default baseline estimation range
    
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi()); //add systematic method (ktrue if is default to compare to, kfalse if is used to estimated the variations)
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.25*TMath::Pi(),0.375*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.375*TMath::Pi(),0.5*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.5*TMath::Pi(),0.625*TMath::Pi());
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.375*TMath::Pi(),0.625*TMath::Pi());
    //systfitter->AddSystematicMode(AliHFCorrFitSystematics::kLowestPoint,kFALSE);
    //systfitter->AddSystematicMode(AliHFCorrFitSystematics::k2PointsAtPiHalf,kFALSE);
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kMinVar,kFALSE);
    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kMaxVar,kFALSE);
    //systfitter->AddSystematicMode(AliHFCorrFitSystematics::kBinCount,kFALSE);
    cout << "crash 3 " << endl;
    //systfitter->CheckBaselineRanges();
    //return;
    systfitter->SetAssociatedTrackPtMin(0.3); // set the min pt assoc (only for plotting purposes)
    systfitter->SetAssociatedTrackPtMax(99); // set the max pt assoc (only for plotting purposes)
    
 /*   Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
   // ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt3to5assocPt0.3to99.0.root");
      ishistoadded = systfitter->AddHistoToFit(Form("%sCanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt3to5assocPt0.3to99.0.root",path.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.11); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.13); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(3); systfitter->SetDptUpEdge(5); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    cout << "crash 4 " << endl;*/
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
  //  ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt5to8assocPt0.3to99.0.root");
      ishistoadded = systfitter->AddHistoToFit(Form("%s/CanvaAndVariedHisto%sAverageDzeroDstarDplus_pPb_Pt5to8assocPt0.3to99.0.root",path.Data(),strAverage.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.09); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.14); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(5); systfitter->SetDptUpEdge(8); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
   // ishistoadded = systfitter->AddHistoToFit("~/Desktop/FinalPlots/ReflectedAverages/CanvaAndVariedHistoWeightedAverageDzeroDstarDplus_pPb_Pt8to16assocPt0.3to99.0.root");
      ishistoadded = systfitter->AddHistoToFit(Form("%s/CanvaAndVariedHisto%sAverageDzeroDstarDplus_pPb_Pt8to16assocPt0.3to99.0.root",path.Data(),strAverage.Data()));
    systfitter->SetMinCorrelatedSystematicsDeltaPhi(0.1); // set min dphi scale syst uncertainty for the related histo
    systfitter->SetMaxCorrelatedSystematicsDeltaPhi(0.13); // set max dphi scale syst uncertainty for the related histo
    systfitter->SetDptLowEdge(8); systfitter->SetDptUpEdge(16); // set the pt range for the related histo
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    
    
    
    
    
    Bool_t setupfitter = systfitter->SetUpSystematicsFitter(); // set up the fitter (initializes variables)
    if(!setupfitter) {cout << "cannot set up fitter " << endl;
        return;
    }
    Bool_t fitted = systfitter->RunFits(); // main function, does all the calculations
    if(!fitted) {cout << "something went wrong with the fits - sorry " << endl;
        return;
    }
    
    systfitter->PrintAllSystematicsOnShell(); // print all the systematics on shell
    
    if(saveAwaySide)systfitter->DrawFinalPlots(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE);
    else systfitter->DrawFinalPlots(); // draws the final trend plots
    
    // systfitter->DrawFinalPlots(kTRUE,kFALSE,kFALSE,kFALSE,kFALSE);
    gROOT->ProcessLine(Form(".!mkdir -p %s/Trends_pPb",outputpath.Data()));
    systfitter->SetOutputDirectory(Form("%s/Trends_pPb",outputpath.Data()));     
    systfitter->SaveCanvasesDotC();
    systfitter->SaveCanvasesDotRoot();
    systfitter->SaveCanvasesDotPng();
    systfitter->SaveCanvasesDotEps();
    systfitter->SaveCanvasesDotPdf();
    
    if(saveAwaySide)  systfitter->DrawFinalCorrelationPlot(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE); // draws the final correlation plots
    else systfitter->DrawFinalCorrelationPlot(); // draws the final correlation plots
    
    // systfitter->CheckHisto(1);
}





//////// TEMPLATE CASE ///////////////
void SystematicsMC_pPb_highpthad(Bool_t useReflected,TString strTemplRootName){
    
  
    TString path = inputpath;
    TString mcCase;
    if(strTemplRootName.Contains("Perugia2010"))mcCase="Perugia2010";
    else if(strTemplRootName.Contains("Perugia2011"))mcCase="Perugia2011";
    else if(strTemplRootName.Contains("Perugia0"))mcCase="Perugia0";
    else if(strTemplRootName.Contains("Herwig"))mcCase="Herwig";
    else if(strTemplRootName.Contains("PYTHIA8"))mcCase="PYTHIA8";
    else if(strTemplRootName.Contains("POWHEG")){
      mcCase="POWHEG";
      strTemplRootName.Remove(0,3);
    }
    else if(strTemplRootName.Contains("EPOS3"))mcCase="EPOS3";
    else {
      Printf("MC case not foreseen");
      return;
    }
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    AliHFCorrFitSystematics * systfitter = new AliHFCorrFitSystematics();
    
    // set up the sysfitter
    gROOT->ProcessLine(Form(".!mkdir -p %s/FitOutputs_pPb",outputpath.Data()));
    systfitter->SetIsFittingTemplate(kTRUE); 
    systfitter->CreateOutputFile(Form("%s/FitOutputs_pPb/FitSystematics_pPb_%s_1.0_99.0.root",outputpath.Data(),mcCase.Data())); // set the name of the output file
    systfitter->SetUseCorrelatedSystematics(kFALSE); // kTRUE (will use the Delta phi independent systematics (from the correlation plots) to sum in quadrature for the final systematics)
    systfitter->SetUseCorrelatedSystematicsForWidths(kFALSE); //Delta phi independent systematics for the final systematics, ktrue will use them also for the widths, kfalse not (in my opinion, we should use kfalse)
    systfitter->SetCombineSystematicsMode(AliHFCorrFitSystematics::kRMS); // set method to compute the systematics due to the fixing of the baseline
    systfitter->SetIspPb(kFALSE); // kTrue for pPb, kFALSE for pp
    //    systfitter->Setv2ForSystematics(0.1, 0.1); // Set v2 values for the systematics (v2had, v2 D)
    systfitter->SetPlotV2SystSeparately(kFALSE); // set true if you want the v2 box to be plotted separately

    //   systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    if(useReflected) systfitter->SetFitRange(0,TMath::Pi()); // set the fit range
    if(!useReflected) systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    systfitter->SetReferenceBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi()); // set the default baseline estimation range
    
    if(useReflected){
      if(strTemplRootName.Contains("EPOS3")){
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),1);//,AliHFCorrFitter::kConstThreeGausPeriodicity); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2);//,AliHFCorrFitter::kConstThreeGausPeriodicity);
      }
      else{
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),1); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2); 
      }
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),0);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),1,AliHFCorrFitter::kConstThreeGausPeriodicity);  
    }
    else {
      if(strTemplRootName.Contains("EPOS3")){
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),2);//,AliHFCorrFitter::kConstThreeGausPeriodicity); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),4);//,AliHFCorrFitter::kConstThreeGausPeriodicity);
      }
      else{
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),2); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),4); 
      }
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),0);  
      //systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2,AliHFCorrFitter::kConstThreeGausPeriodicity);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2,AliHFCorrFitter::kConstThreeGausPeriodicity);  
    }
    //    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi());
    
    //systfitter->CheckBaselineRanges();
    //return;
    systfitter->SetAssociatedTrackPtMin(1); // set the min pt assoc (only for plotting purposes)
    systfitter->SetAssociatedTrackPtMax(99); // set the max pt assoc (only for plotting purposes)
    
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s3To5_ptAssall1.0to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s3To5_ptAssall1.0to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(3); systfitter->SetDptUpEdge(5); // set the pt range for the related histo


    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s5To8_ptAssall1.0to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s5To8_ptAssall1.0to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(5); systfitter->SetDptUpEdge(8); // set the pt range for the related histo


    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s8To16_ptAssall1.0to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s8To16_ptAssall1.0to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(8); systfitter->SetDptUpEdge(16); // set the pt range for the related histo

        
    
    Bool_t setupfitter = systfitter->SetUpSystematicsFitter(); // set up the fitter (initializes variables)
    if(!setupfitter) {cout << "cannot set up fitter " << endl;
        return;
    }
    Bool_t fitted = systfitter->RunFits(); // main function, does all the calculations
    if(!fitted) {cout << "something went wrong with the fits - sorry " << endl;
        return;
    }
    
    systfitter->PrintAllSystematicsOnShell(); // print all the systematics on shell
   
    if(saveAwaySide)systfitter->DrawFinalPlots(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE);
    else systfitter->DrawFinalPlots(); // draws the final trend plots
    
    // systfitter->DrawFinalPlots(kTRUE,kFALSE,kFALSE,kFALSE,kFALSE);
    gROOT->ProcessLine(Form(".!mkdir -p %s/Trends_pPb/%s",outputpath.Data(),mcCase.Data()));
    systfitter->SetOutputDirectory(Form("%s/Trends_pPb/%s/",outputpath.Data(),mcCase.Data()));

    systfitter->SaveCanvasesDotC();
    systfitter->SaveCanvasesDotRoot();
    systfitter->SaveCanvasesDotPng();
    systfitter->SaveCanvasesDotEps();
    systfitter->SaveCanvasesDotPdf();
    
    if(saveAwaySide)systfitter->DrawFinalCorrelationPlot(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE); // draws the final correlation plots
	else systfitter->DrawFinalCorrelationPlot(); // draws the final correlation plots

    // systfitter->CheckHisto(1);
}




void SystematicsMC_pPb_lowpthad(Bool_t useReflected,TString strTemplRootName){
    
  
    TString path = inputpath;
    TString mcCase;
    if(strTemplRootName.Contains("Perugia2010"))mcCase="Perugia2010";
    else if(strTemplRootName.Contains("Perugia2011"))mcCase="Perugia2011";
    else if(strTemplRootName.Contains("Perugia0"))mcCase="Perugia0";
    else if(strTemplRootName.Contains("Herwig"))mcCase="Herwig";
    else if(strTemplRootName.Contains("PYTHIA8"))mcCase="PYTHIA8";
    else if(strTemplRootName.Contains("POWHEG")){
      mcCase="POWHEG";
      strTemplRootName.Remove(0,3);
    }
    else if(strTemplRootName.Contains("EPOS3"))mcCase="EPOS3";
    else {
      Printf("MC case not foreseen");
      return;
    }
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    AliHFCorrFitSystematics * systfitter = new AliHFCorrFitSystematics();
    
    // set up the sysfitter
    gROOT->ProcessLine(Form(".!mkdir -p %s/FitOutputs_pPb",outputpath.Data()));
    systfitter->SetIsFittingTemplate(kTRUE); 
    systfitter->CreateOutputFile(Form("%s/FitOutputs_pPb/FitSystematics_pPb_%s_0.3_1.0.root",outputpath.Data(),mcCase.Data())); // set the name of the output file
    systfitter->SetUseCorrelatedSystematics(kFALSE); // kTRUE (will use the Delta phi independent systematics (from the correlation plots) to sum in quadrature for the final systematics)
    systfitter->SetUseCorrelatedSystematicsForWidths(kFALSE); //Delta phi independent systematics for the final systematics, ktrue will use them also for the widths, kfalse not (in my opinion, we should use kfalse)
    systfitter->SetCombineSystematicsMode(AliHFCorrFitSystematics::kRMS); // set method to compute the systematics due to the fixing of the baseline
    systfitter->SetIspPb(kFALSE); // kTrue for pPb, kFALSE for pp
    //    systfitter->Setv2ForSystematics(0.1, 0.1); // Set v2 values for the systematics (v2had, v2 D)
    systfitter->SetPlotV2SystSeparately(kFALSE); // set true if you want the v2 box to be plotted separately
    //   systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    if(useReflected) systfitter->SetFitRange(0,TMath::Pi()); // set the fit range
    if(!useReflected) systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    systfitter->SetReferenceBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi()); // set the default baseline estimation range
    
    if(useReflected){
      if(strTemplRootName.Contains("EPOS3")){
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),1);//,AliHFCorrFitter::kConstThreeGausPeriodicity); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2);//,AliHFCorrFitter::kConstThreeGausPeriodicity);
      }
      else{
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),1); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2); 
      }
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),0);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),1,AliHFCorrFitter::kConstThreeGausPeriodicity);  
    }
    else {
      if(strTemplRootName.Contains("EPOS3")){
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),2);//,AliHFCorrFitter::kConstThreeGausPeriodicity); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),4);//,AliHFCorrFitter::kConstThreeGausPeriodicity);
      }
      else{
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),2); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),4); 
      }
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),0);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2,AliHFCorrFitter::kConstThreeGausPeriodicity);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2,AliHFCorrFitter::kConstThreeGausPeriodicity);  
    }
    //    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi());
    
    //systfitter->CheckBaselineRanges();
    //return;
    systfitter->SetAssociatedTrackPtMin(0.3); // set the min pt assoc (only for plotting purposes)
    systfitter->SetAssociatedTrackPtMax(1.); // set the max pt assoc (only for plotting purposes)
    
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s3To5_ptAssall0.3to1.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s3To5_ptAssall0.3to1.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(3); systfitter->SetDptUpEdge(5); // set the pt range for the related histo


    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s5To8_ptAssall0.3to1.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s5To8_ptAssall0.3to1.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(5); systfitter->SetDptUpEdge(8); // set the pt range for the related histo

    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s8To16_ptAssall0.3to1.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s8To16_ptAssall0.3to1.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(8); systfitter->SetDptUpEdge(16); // set the pt range for the related histo
        
    
    Bool_t setupfitter = systfitter->SetUpSystematicsFitter(); // set up the fitter (initializes variables)
    if(!setupfitter) {cout << "cannot set up fitter " << endl;
        return;
    }
    Bool_t fitted = systfitter->RunFits(); // main function, does all the calculations
    if(!fitted) {cout << "something went wrong with the fits - sorry " << endl;
        return;
    }
    
    systfitter->PrintAllSystematicsOnShell(); // print all the systematics on shell
   
    if(saveAwaySide)systfitter->DrawFinalPlots(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE);
    else systfitter->DrawFinalPlots(); // draws the final trend plots
    
    // systfitter->DrawFinalPlots(kTRUE,kFALSE,kFALSE,kFALSE,kFALSE);
    gROOT->ProcessLine(Form(".!mkdir -p %s/Trends_pPb/%s",outputpath.Data(),mcCase.Data()));
    systfitter->SetOutputDirectory(Form("%s/Trends_pPb/%s/",outputpath.Data(),mcCase.Data()));

    systfitter->SaveCanvasesDotC();
    systfitter->SaveCanvasesDotRoot();
    systfitter->SaveCanvasesDotPng();
    systfitter->SaveCanvasesDotEps();
    systfitter->SaveCanvasesDotPdf();
    
    if(saveAwaySide)systfitter->DrawFinalCorrelationPlot(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE); // draws the final correlation plots
    else systfitter->DrawFinalCorrelationPlot(); // draws the final correlation plots

    // systfitter->CheckHisto(1);
}




void SystematicsMC_pPb_integrpthad(Bool_t useReflected,TString strTemplRootName){
    
  
    TString path = inputpath;
    TString mcCase;
    if(strTemplRootName.Contains("Perugia2010"))mcCase="Perugia2010";
    else if(strTemplRootName.Contains("Perugia2011"))mcCase="Perugia2011";
    else if(strTemplRootName.Contains("Perugia0"))mcCase="Perugia0";
    else if(strTemplRootName.Contains("Herwig"))mcCase="Herwig";
    else if(strTemplRootName.Contains("PYTHIA8"))mcCase="PYTHIA8";
    else if(strTemplRootName.Contains("POWHEG")){
      mcCase="POWHEG";
      strTemplRootName.Remove(0,3);
    }
    else if(strTemplRootName.Contains("EPOS3"))mcCase="EPOS3";
    else {
      Printf("MC case not foreseen");
      return;
    }
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    AliHFCorrFitSystematics * systfitter = new AliHFCorrFitSystematics();
    
    // set up the sysfitter
    gROOT->ProcessLine(Form(".!mkdir -p %s/FitOutputs_pPb",outputpath.Data()));
    systfitter->SetIsFittingTemplate(kTRUE); 
    systfitter->CreateOutputFile(Form("%s/FitOutputs_pPb/FitSystematics_pPb_%s_0.3_99.0.root",outputpath.Data(),mcCase.Data())); // set the name of the output file
    systfitter->SetUseCorrelatedSystematics(kFALSE); // kTRUE (will use the Delta phi independent systematics (from the correlation plots) to sum in quadrature for the final systematics)
    systfitter->SetUseCorrelatedSystematicsForWidths(kFALSE); //Delta phi independent systematics for the final systematics, ktrue will use them also for the widths, kfalse not (in my opinion, we should use kfalse)
    systfitter->SetCombineSystematicsMode(AliHFCorrFitSystematics::kRMS); // set method to compute the systematics due to the fixing of the baseline
    systfitter->SetIspPb(kFALSE); // kTrue for pPb, kFALSE for pp
    //    systfitter->Setv2ForSystematics(0.1, 0.1); // Set v2 values for the systematics (v2had, v2 D)
    systfitter->SetPlotV2SystSeparately(kFALSE); // set true if you want the v2 box to be plotted separately
    //   systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    if(useReflected) systfitter->SetFitRange(0,TMath::Pi()); // set the fit range
    if(!useReflected) systfitter->SetFitRange(-0.5*TMath::Pi(),1.5*TMath::Pi()); // set the fit range
    
    systfitter->SetReferenceBaselineEstimationRange(0.25*TMath::Pi(),0.5*TMath::Pi()); // set the default baseline estimation range
    
    if(useReflected){
      if(strTemplRootName.Contains("EPOS3")){
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),1);//,AliHFCorrFitter::kConstThreeGausPeriodicity); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2);//,AliHFCorrFitter::kConstThreeGausPeriodicity);
      }
      else{
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),1); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2); 
      }
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),0);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),1,AliHFCorrFitter::kConstThreeGausPeriodicity);  
    }
    else {
      if(strTemplRootName.Contains("EPOS3")){
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),2);//,AliHFCorrFitter::kConstThreeGausPeriodicity); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),4);//,AliHFCorrFitter::kConstThreeGausPeriodicity);
      }
      else{
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kTRUE,0.25*TMath::Pi(),0.5*TMath::Pi(),2); 
	systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),4); 
      }
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),0);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kNLowest,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2,AliHFCorrFitter::kConstThreeGausPeriodicity);  
      //      systfitter->AddSystematicMode(AliHFCorrFitSystematics::kFree,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi(),2,AliHFCorrFitter::kConstThreeGausPeriodicity);  
    }
    //    systfitter->AddSystematicMode(AliHFCorrFitSystematics::kTransverse,kFALSE,0.25*TMath::Pi(),0.5*TMath::Pi());
    
    //systfitter->CheckBaselineRanges();
    //return;
    systfitter->SetAssociatedTrackPtMin(0.3); // set the min pt assoc (only for plotting purposes)
    systfitter->SetAssociatedTrackPtMax(99.); // set the max pt assoc (only for plotting purposes)
    
    Bool_t ishistoadded;
    // add histograms to fit (different pt bins of the
    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s3To5_ptAssall0.3to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s3To5_ptAssall0.3to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(3); systfitter->SetDptUpEdge(5); // set the pt range for the related histo

    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s5To8_ptAssall0.3to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s5To8_ptAssall0.3to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(5); systfitter->SetDptUpEdge(8); // set the pt range for the related histo

    ishistoadded = systfitter->AddHistoToFit(Form("%s/%s8To16_ptAssall0.3to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"","hCorrDeltaPhi","","");
    if(!ishistoadded){
      ishistoadded = systfitter->AddHistoToFit(Form("%s/%s8To16_ptAssall0.3to99.0_DeltaEta10.root",path.Data(),strTemplRootName.Data()),"cDeltaPhi","hCorrDeltaPhi","","");
    }
    if(!ishistoadded){cout << "cannot add histogram, exiting... " << endl; return;}
    systfitter->SetDptLowEdge(8); systfitter->SetDptUpEdge(16); // set the pt range for the related histo
        
    
    Bool_t setupfitter = systfitter->SetUpSystematicsFitter(); // set up the fitter (initializes variables)
    if(!setupfitter) {cout << "cannot set up fitter " << endl;
        return;
    }
    Bool_t fitted = systfitter->RunFits(); // main function, does all the calculations
    if(!fitted) {cout << "something went wrong with the fits - sorry " << endl;
        return;
    }
    
    systfitter->PrintAllSystematicsOnShell(); // print all the systematics on shell
   
    if(saveAwaySide)systfitter->DrawFinalPlots(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE);
    else systfitter->DrawFinalPlots(); // draws the final trend plots
    
    // systfitter->DrawFinalPlots(kTRUE,kFALSE,kFALSE,kFALSE,kFALSE);
    gROOT->ProcessLine(Form(".!mkdir -p %s/Trends_pPb/%s",outputpath.Data(),mcCase.Data()));
    systfitter->SetOutputDirectory(Form("%s/Trends_pPb/%s/",outputpath.Data(),mcCase.Data()));

    systfitter->SaveCanvasesDotC();
    systfitter->SaveCanvasesDotRoot();
    systfitter->SaveCanvasesDotPng();
    systfitter->SaveCanvasesDotEps();
    systfitter->SaveCanvasesDotPdf();
    
    if(saveAwaySide)systfitter->DrawFinalCorrelationPlot(kTRUE,kTRUE,kTRUE,kTRUE,kTRUE); // draws the final correlation plots
    else systfitter->DrawFinalCorrelationPlot(); // draws the final correlation plots

    // systfitter->CheckHisto(1);
}
