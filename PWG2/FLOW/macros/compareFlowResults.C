//type of analysis can be: ESD, AOD, MC, ESDMC0, ESDMC1
//const TString type = "ESD"; 

enum libModes {mLocal,mLocalSource};
//mLocal: Analyze data on your computer using aliroot
//mLocalSource: Analyze data on your computer using root + source files

//void compareFlowResults(TString type="",Int_t mode=mLocalSource)
void compareFlowResults(TString type="",Int_t mode=mLocal)
{ 

  // load needed libraries:                       
  LoadPlotLibraries(mode);

  Bool_t plotLegendIntFlow = kTRUE; // plot legend with average multiplicity and number of events for each method in all plots for integrated flow

  //==================================================================================
  //             set here which plots will be shown by default
  //==================================================================================
  Bool_t plotIntFlow = kTRUE;                  // integrated flow (no-name) // to be improved
  Bool_t plotIntFlowRelativeToMC = kTRUE;      // plot |v{MC}-v{method}/v{MC}| for integrated flow (no-name) // to be improved
  // RP = particles used to determine the reaction plane
  Bool_t plotIntFlowRP = kTRUE;             // integrated flow RP
  Bool_t plotIntFlowRelativeToMCRP = kTRUE;   // plot |v{MC}-v{method}/v{MC}| for integrated flow (RP) // to be improved
  Bool_t plotDiffFlowPtRP = kTRUE;             // differential flow (Pt,RP)
  Bool_t plotDiffFlowEtaRP = kTRUE;             // differential flow (Eta,RP)
  Bool_t plotDiffFlowPtRelativeToMCRP = kTRUE;  // plot |v{MC}-v{method}/v{MC}| as a function of pt for RPs   
  // POI = particle of interest
  Bool_t plotIntFlowPOI = kTRUE;              // integrated flow POI
  Bool_t plotIntFlowRelativeToMCPOI = kTRUE;   // plot |v{MC}-v{method}/v{MC}| for integrated flow (POI) // to be improved  
  Bool_t plotDiffFlowPtPOI = kTRUE;           // differential flow (Pt,POI)
  Bool_t plotDiffFlowEtaPOI = kTRUE;         // differential flow (Eta,POI)
  //==================================================================================
  
  
  //==================================================================================
  // set here which methods will be plotted by default for differential flow (Pt,RP):
  Bool_t plotMCPtRP       = kFALSE;
  Bool_t plotSPPtRP       = kTRUE;
  Bool_t plotGFC2PtRP     = kTRUE;
  Bool_t plotGFC4PtRP     = kTRUE;
  Bool_t plotGFC6PtRP     = kTRUE;
  Bool_t plotGFC8PtRP     = kTRUE;
  Bool_t plotQC2PtRP      = kTRUE;
  Bool_t plotQC4PtRP      = kTRUE;
  Bool_t plotQC6PtRP      = kFALSE; // not calculated yet
  Bool_t plotQC8PtRP      = kFALSE; // not calculated yet
  Bool_t plotLYZ2SUMPtRP  = kTRUE;
  Bool_t plotLYZ2PRODPtRP = kTRUE;
  Bool_t plotLYZEPPtRP    = kTRUE; 
  
  // set here which methods will be plotted by default for differential flow (Eta,RP):
  Bool_t plotMCEtaRP       = kFALSE;
  Bool_t plotSPEtaRP       = kTRUE;
  Bool_t plotGFC2EtaRP     = kTRUE;
  Bool_t plotGFC4EtaRP     = kTRUE;
  Bool_t plotGFC6EtaRP     = kTRUE;
  Bool_t plotGFC8EtaRP     = kTRUE;
  Bool_t plotQC2EtaRP      = kTRUE;
  Bool_t plotQC4EtaRP      = kTRUE;
  Bool_t plotQC6EtaRP      = kFALSE; // not calculated yet
  Bool_t plotQC8EtaRP      = kFALSE; // not calculated yet
  Bool_t plotLYZ2SUMEtaRP  = kTRUE;
  Bool_t plotLYZ2PRODEtaRP = kTRUE;
  Bool_t plotLYZEPEtaRP    = kTRUE;
  
  // set here which methods will be plotted by default for |v{MC}-v{method}/v{MC}| as a function of pt for RPs 
  Bool_t plotSPRelativeToMCRP       = kTRUE;
  Bool_t plotGFC2RelativeToMCRP     = kTRUE;
  Bool_t plotGFC4RelativeToMCRP     = kTRUE;
  Bool_t plotGFC6RelativeToMCRP     = kTRUE;
  Bool_t plotGFC8RelativeToMCRP     = kTRUE;
  Bool_t plotQC2RelativeToMCRP      = kTRUE;
  Bool_t plotQC4RelativeToMCRP      = kTRUE;
  Bool_t plotQC6RelativeToMCRP      = kFALSE; // not calculated yet
  Bool_t plotQC8RelativeToMCRP      = kFALSE; // not calculated yet
  Bool_t plotLYZ2SUMRelativeToMCRP  = kTRUE;
  Bool_t plotLYZ2PRODRelativeToMCRP = kTRUE;
  Bool_t plotLYZEPRelativeToMCRP    = kTRUE;  
  
  // set here which methods will be plotted by default for differential flow (Pt,POI):
  Bool_t plotMCPtPOI       = kFALSE;
  Bool_t plotSPPtPOI       = kTRUE;
  Bool_t plotGFC2PtPOI     = kTRUE;
  Bool_t plotGFC4PtPOI     = kTRUE;
  Bool_t plotGFC6PtPOI     = kTRUE;
  Bool_t plotGFC8PtPOI     = kTRUE;
  Bool_t plotQC2PtPOI      = kTRUE;
  Bool_t plotQC4PtPOI      = kTRUE;
  Bool_t plotQC6PtPOI      = kFALSE; // not calculated yet
  Bool_t plotQC8PtPOI      = kFALSE; // not calculated yet
  Bool_t plotLYZ2SUMPtPOI  = kTRUE;
  Bool_t plotLYZ2PRODPtPOI = kTRUE;
  Bool_t plotLYZEPPtPOI    = kTRUE; 
  
  // set here which methods will be plotted by default for differential flow (Eta,POI):
  Bool_t plotMCEtaPOI       = kFALSE;
  Bool_t plotSPEtaPOI       = kTRUE;
  Bool_t plotGFC2EtaPOI     = kTRUE;
  Bool_t plotGFC4EtaPOI     = kTRUE;
  Bool_t plotGFC6EtaPOI     = kTRUE;
  Bool_t plotGFC8EtaPOI     = kTRUE;
  Bool_t plotQC2EtaPOI      = kTRUE;
  Bool_t plotQC4EtaPOI      = kTRUE;
  Bool_t plotQC6EtaPOI      = kFALSE; // not calculated yet
  Bool_t plotQC8EtaPOI      = kFALSE; // not calculated yet
  Bool_t plotLYZ2SUMEtaPOI  = kTRUE;
  Bool_t plotLYZ2PRODEtaPOI = kTRUE;
  Bool_t plotLYZEPEtaPOI    = kTRUE;
  //==================================================================================
 
  
  //==================================================================================  
  // cosmetics: marker style (see TMarker) and color (see TAttFill) for each method:
  // MC:
  Int_t markerStyleMC = 20; // full circle
  Int_t markerColorMC = kRed;
  // SP:
  Int_t markerStyleSP = 3; // star
  Int_t markerColorSP = kViolet-6;
  // GFC{2}
  Int_t markerStyleGFC2 = 21; // full square
  Int_t markerColorGFC2 = kAzure-7;
  // GFC{4}
  Int_t markerStyleGFC4 = 20; // full circle
  Int_t markerColorGFC4 = kAzure+3;
  // GFC{6}
  Int_t markerStyleGFC6 = 25; // open square
  Int_t markerColorGFC6 = kAzure-7;
  // GFC{8}
  Int_t markerStyleGFC8 = 24; // open circle
  Int_t markerColorGFC8 = kAzure+3;
  // QC{2}
  Int_t markerStyleQC2 = 21; // full square
  Int_t markerColorQC2 = kOrange-7;
  // QC{4}
  Int_t markerStyleQC4 = 20; // full circle
  Int_t markerColorQC4 = kOrange+3;
  // QC{6}
  Int_t markerStyleQC6 = 25; // open square
  Int_t markerColorQC6 = kOrange-7;
  // QC{8}
  Int_t markerStyleQC8 = 24; // open circle
  Int_t markerColorQC8 = kOrange+3;
  // LYZ2SUM
  Int_t markerStyleLYZ2SUM = 22; // full triangle
  Int_t markerColorLYZ2SUM = kYellow+3;
  // LYZ2PROD
  Int_t markerStyleLYZ2PROD = 22; // full triangle
  Int_t markerColorLYZ2PROD = kGreen+3;
  // LYZEP
  Int_t markerStyleLYZEP = 26; // open triangle
  Int_t markerColorLYZEP = kYellow+3; 
  //==================================================================================


  //==================================================================================  
  // set here which result goes in which bin in the plot for integrated flow (no-name) 
  // MC:
  Int_t binMC = 1; 
  // SP:
  Int_t binSP = 2;
  // GFC{2}
  Int_t binGFC2 = 3; 
  // GFC{4}
  Int_t binGFC4 = 5; 
  // GFC{6}
  Int_t binGFC6 = 7; 
  // GFC{8}
  Int_t binGFC8 = 9; 
  // QC{2}
  Int_t binQC2 = 4; 
  // QC{4}
  Int_t binQC4 = 6; 
  // QC{6}
  Int_t binQC6 = 8; 
  // QC{8}
  Int_t binQC8 = 10; 
  // FQD 
  Int_t binFQD = 11; 
  // LYZ1SUM
  Int_t binLYZ1SUM = 12; 
  // LYZ1PROD
  Int_t binLYZ1PROD = 13; 
  // LYZEP
  Int_t binLYZEP = 14; 
  //==================================================================================


  //==================================================================================  
  // set here which result goes in which bin in the plot for integrated flow (RP) 
  // MC:
  Int_t binMCRP = 1; 
  // SP:
  Int_t binSPRP = 2;
  // GFC{2}
  Int_t binGFC2RP = 3; 
  // GFC{4}
  Int_t binGFC4RP = 5; 
  // GFC{6}
  Int_t binGFC6RP = 7; 
  // GFC{8}
  Int_t binGFC8RP = 9; 
  // QC{2}
  Int_t binQC2RP = 4; 
  // QC{4}
  Int_t binQC4RP = 6; 
  // QC{6}
  Int_t binQC6RP = 8; 
  // QC{8}
  Int_t binQC8RP = 10; 
  // FQD 
  Int_t binFQDRP = 11; 
  // LYZ2SUM
  Int_t binLYZ2SUMRP = 12; 
  // LYZ2PROD
  Int_t binLYZ2PRODRP = 13; 
  // LYZEP
  Int_t binLYZEPRP = 14; 
  //==================================================================================


  //==================================================================================  
  // set here which result goes in which bin in the plot for integrated flow (POI) 
  // MC:
  Int_t binMCPOI = 1; 
  // SP:
  Int_t binSPPOI = 2;
  // GFC{2}
  Int_t binGFC2POI = 3; 
  // GFC{4}
  Int_t binGFC4POI = 5; 
  // GFC{6}
  Int_t binGFC6POI = 7; 
  // GFC{8}
  Int_t binGFC8POI = 9; 
  // QC{2}
  Int_t binQC2POI = 4; 
  // QC{4}
  Int_t binQC4POI = 6; 
  // QC{6}
  Int_t binQC6POI = 8; 
  // QC{8}
  Int_t binQC8POI = 10; 
  // FQD 
  Int_t binFQDPOI = 11; 
  // LYZ2SUM
  Int_t binLYZ2SUMPOI = 12;  
  // LYZ2PROD
  Int_t binLYZ2PRODPOI = 13;    
  // LYZEP
  Int_t binLYZEPPOI = 14; 
  //==================================================================================
 
                                        
  //==================================================================================
  //                         accessing output files
  //==================================================================================
  TString outputFileName = "AnalysisResults.root"; // final output file name holding final results for large statistics sample
  // access the merged, large statistics file obtained with macro mergeOutput.C:
  TString pwd(gSystem->pwd());
  pwd+="/";
  pwd+=outputFileName.Data();
  TFile *outputFile = NULL;
  if(gSystem->AccessPathName(pwd.Data(),kFileExists))
  {
   cout<<"WARNING: You do not have an output file:"<<endl;
   cout<<"         "<<pwd.Data()<<endl;
   exit(0);
  } else 
    {
     outputFile = TFile::Open(pwd.Data(),"READ");
    }
  
  //open the output files for each method:
  TString inputFileNameMCEP = "outputMCEPanalysis";
  TFile* fileMCEP = (TFile*)outputFile->FindObjectAny((inputFileNameMCEP.Append(type)).Data());
   
  TString inputFileNameSP = "outputSPanalysis";
  TFile* fileSP = (TFile*)outputFile->FindObjectAny((inputFileNameSP.Append(type)).Data());
    
  TString inputFileNameLYZ1SUM = "outputLYZ1SUManalysis";
  TFile* fileLYZ1SUM = (TFile*)outputFile->FindObjectAny((inputFileNameLYZ1SUM.Append(type)).Data());
  
  TString inputFileNameLYZ2SUM = "outputLYZ2SUManalysis";
  TFile* fileLYZ2SUM = (TFile*)outputFile->FindObjectAny((inputFileNameLYZ2SUM.Append(type)).Data());
    
  TString inputFileNameLYZ1PROD = "outputLYZ1PRODanalysis";
  TFile* fileLYZ1PROD = (TFile*)outputFile->FindObjectAny((inputFileNameLYZ1PROD.Append(type)).Data());
  
  TString inputFileNameLYZ2PROD = "outputLYZ2PRODanalysis";
  TFile* fileLYZ2PROD = (TFile*)outputFile->FindObjectAny((inputFileNameLYZ2PROD.Append(type)).Data());
  
  TString inputFileNameLYZEP = "outputLYZEPanalysis";
  TFile* fileLYZEP = (TFile*)outputFile->FindObjectAny((inputFileNameLYZEP.Append(type)).Data());
  
  TString inputFileNameFQD = "outputFQDanalysis";
  TFile* fileFQD = (TFile*)outputFile->FindObjectAny((inputFileNameFQD.Append(type)).Data());
  
  TString inputFileNameGFC = "outputGFCanalysis";
  TFile* fileGFC = (TFile*)outputFile->FindObjectAny((inputFileNameGFC.Append(type)).Data());
  
  TString inputFileNameQC = "outputQCanalysis";
  TFile* fileQC = (TFile*)outputFile->FindObjectAny((inputFileNameQC.Append(type)).Data());
  //==================================================================================
 
 
 
 
  //==================================================================================
  //                                 cosmetics
  //==================================================================================
  //removing the title and stat. box from all histograms:
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  // plot for 'no-name' integrated flow:
  // choosing the style and color of mesh for MC error bands
  Int_t meshStyle = 1001;
  Int_t meshColor = kGray;
  
  // marker style and color  
  Int_t markerStyle = 21;
  Int_t markerColor = kBlack;
  
  // plot for RP's integrated flow:
  // choosing the style and color of mesh for MC error bands
  Int_t meshStyleRP = 1001;
  Int_t meshColorRP = kRed-10;
  
  // marker style and color  
  Int_t markerStyleRP = 21;
  Int_t markerColorRP = kRed-3;
  
  // plot for POI's integrated flow:
  // choosing the style and color of mesh for MC error bands
  Int_t meshStylePOI = 1001;
  Int_t meshColorPOI = kBlue-10;
  
  // marker style and color  
  Int_t markerStylePOI = 21;
  Int_t markerColorPOI = kBlue-3;
  
  // choosing the style and color of mesh for MC error bands in the plots for diff. flow
  // plot for differential flow (Pt,RP):
  Int_t meshStyleDiffFlowPtRP = 1001;
  Int_t meshColorDiffFlowPtRP = kRed-10;
  
  // plot for differential flow (Eta,RP):
  Int_t meshStyleDiffFlowEtaRP = 1001;
  Int_t meshColorDiffFlowEtaRP = kRed-10;
  
  // plot for differential flow (Pt,POI):
  Int_t meshStyleDiffFlowPtPOI = 1001;
  Int_t meshColorDiffFlowPtPOI = kRed-10;

  // plot for differential flow (Eta,POI):
  Int_t meshStyleDiffFlowEtaPOI = 1001;
  Int_t meshColorDiffFlowEtaPOI = kRed-10;
  //==================================================================================
  
  
  
  
  //==================================================================================
  //                             INTEGRATED FLOW (no-name, RP and POI)
  //==================================================================================
  // the number of different methods:
  const Int_t nMethods=14;
  
  // booking the histogram for the integrated flow results from all methods:
  TH1D* intFlowAll = new TH1D("intFlowAll","Integrated Flow",nMethods,0,nMethods);      
  if(!plotLegendIntFlow) intFlowAll->SetLabelSize(0.044,"X");
  // intFlowAll->SetLabelSize(0.036,"Y");
  intFlowAll->SetMarkerStyle(markerStyle);
  intFlowAll->SetMarkerColor(markerColor);
  (intFlowAll->GetXaxis())->SetBinLabel(binMC,"v_{2}{MC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binSP,"v_{2}{SP}");
  (intFlowAll->GetXaxis())->SetBinLabel(binGFC2,"v_{2}{2,GFC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binQC2,"v_{2}{2,QC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binGFC4,"v_{2}{4,GFC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binQC4,"v_{2}{4,QC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binGFC6,"v_{2}{6,GFC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binQC6,"v_{2}{6,QC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binGFC8,"v_{2}{8,GFC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binQC8,"v_{2}{8,QC}");
  (intFlowAll->GetXaxis())->SetBinLabel(binFQD,"v_{2}{FQD}");
  (intFlowAll->GetXaxis())->SetBinLabel(binLYZ1SUM,"v_{2}{LYZ,sum}");
  (intFlowAll->GetXaxis())->SetBinLabel(binLYZ1PROD,"v_{2}{LYZ,prod}");
  (intFlowAll->GetXaxis())->SetBinLabel(binLYZEP,"v_{2}{LYZEP}");
  
  // booking the graph to store flow values and errors from all methods:    
  Double_t x[nMethods] = {0.};
  for(Int_t i=0;i<nMethods;i++)
  {
   x[i]=i+0.5;
  }
  Double_t xError[nMethods] = {0.};
  Double_t flowValue[nMethods] = {0.};
  Double_t flowError[nMethods] = {0.};
  Double_t flowValueRP[nMethods] = {0.};
  Double_t flowErrorRP[nMethods] = {0.};
  Double_t flowValuePOI[nMethods] = {0.};
  Double_t flowErrorPOI[nMethods] = {0.};
  
  // accessing the results for integrated flow for each method:
  // MCEP = Monte Carlo Event Plane
  TList *pListMCEP = NULL;
  AliFlowCommonHist *mcepCommonHist = NULL;
  AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
  if(fileMCEP) {
    fileMCEP->GetObject("cobjMCEP",pListMCEP); 
    if(pListMCEP) {
      mcepCommonHist = dynamic_cast<AliFlowCommonHist*> (pListMCEP->FindObject("AliFlowCommonHistMCEP"));
      mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListMCEP->FindObject("AliFlowCommonHistResultsMCEP"));
      if(mcepCommonHistRes) {
	flowValue[binMC-1] = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
	flowError[binMC-1] = (mcepCommonHistRes->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binMCRP-1] = (mcepCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binMCRP-1] = (mcepCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binMCPOI-1] = (mcepCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binMCPOI-1] = (mcepCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
  
  // SP = Scalar Product
  TList *pListSP = NULL;
  AliFlowCommonHist *spCommonHist = NULL;
  AliFlowCommonHistResults *spCommonHistRes = NULL; 
  if(fileSP) {
    fileSP->GetObject("cobjSP",pListSP); 
    if(pListSP) {
      spCommonHist = dynamic_cast<AliFlowCommonHist*> (pListSP->FindObject("AliFlowCommonHistSP"));
      spCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListSP->FindObject("AliFlowCommonHistResultsSP"));
      if(spCommonHistRes) {
	flowValue[binSP-1] = (spCommonHistRes->GetHistIntFlow())->GetBinContent(1);
	flowError[binSP-1] = (spCommonHistRes->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binSPRP-1] = (spCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binSPRP-1] = (spCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binSPPOI-1] = (spCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binSPPOI-1] = (spCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
  
  // LYZ1SUM = Lee-Yang Zeros (1st run, sum) is used to get only 'no-name' integrated flow
  TList *pListLYZ1SUM = NULL;
  AliFlowCommonHist *lyz1sumCommonHist = NULL;
  AliFlowCommonHistResults *lyz1sumCommonHistRes = NULL; 
  if(fileLYZ1SUM) {
    fileLYZ1SUM->GetObject("cobjLYZ1SUM",pListLYZ1SUM); 
    if(pListLYZ1SUM) {
      lyz1sumCommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZ1SUM->FindObject("AliFlowCommonHistLYZ1SUM"));
      lyz1sumCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ1SUM->FindObject("AliFlowCommonHistResultsLYZ1SUM"));
      if(lyz1sumCommonHistRes) {
	flowValue[binLYZ1SUM-1] = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinContent(1);
	flowError[binLYZ1SUM-1] = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinError(1);
      }
    }
  }
  
  // LYZ2SUM = Lee-Yang Zeros (2nd run, sum) is used to get RP's and POI's integrated flow
  TList *pListLYZ2SUM = NULL;
  AliFlowCommonHist *lyz2sumCommonHist = NULL;
  AliFlowCommonHistResults *lyz2sumCommonHistRes = NULL; 
  if(fileLYZ2SUM) {
    fileLYZ2SUM->GetObject("cobjLYZ2SUM",pListLYZ2SUM); 
    if(pListLYZ2SUM) {
      lyz2sumCommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZ2SUM->FindObject("AliFlowCommonHistLYZ2SUM"));
      lyz2sumCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ2SUM->FindObject("AliFlowCommonHistResultsLYZ2SUM"));
      if(lyz2sumCommonHistRes) {
	flowValueRP[binLYZ2SUMRP-1] = (lyz2sumCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binLYZ2SUMRP-1] = (lyz2sumCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binLYZ2SUMPOI-1] = (lyz2sumCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binLYZ2SUMPOI-1] = (lyz2sumCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
 
  // LYZ1PROD = Lee-Yang Zeros (1st run, product) is used to get only 'no-name' integrated flow
  TList *pListLYZ1PROD = NULL;
  AliFlowCommonHist *lyz1prodCommonHist = NULL;
  AliFlowCommonHistResults *lyz1prodCommonHistRes = NULL; 
  if(fileLYZ1PROD) {
    fileLYZ1PROD->GetObject("cobjLYZ1PROD",pListLYZ1PROD); 
    if(pListLYZ1PROD) {
      lyz1prodCommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZ1PROD->FindObject("AliFlowCommonHistLYZ1PROD"));
      lyz1prodCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ1PROD->FindObject("AliFlowCommonHistResultsLYZ1PROD"));
      if(lyz1prodCommonHistRes) {
	flowValue[binLYZ1PROD-1] = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinContent(1);
	flowError[binLYZ1PROD-1] = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinError(1);
      }
    }
  }
  
  // LYZ2PROD = Lee-Yang Zeros (2nd run, product) is used to get RP's and POI's integrated flow
  TList *pListLYZ2PROD = NULL;
  AliFlowCommonHist *lyz2prodCommonHist = NULL;
  AliFlowCommonHistResults *lyz2prodCommonHistRes = NULL; 
  if(fileLYZ2PROD) {
    fileLYZ2PROD->GetObject("cobjLYZ2PROD",pListLYZ2PROD); 
    if(pListLYZ2PROD) {
      lyz2prodCommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZ2PROD->FindObject("AliFlowCommonHistLYZ2PROD"));
      lyz2prodCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ2PROD->FindObject("AliFlowCommonHistResultsLYZ2PROD"));
      if(lyz2prodCommonHistRes) {
	flowValueRP[binLYZ2PRODRP-1] = (lyz2prodCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binLYZ2PRODRP-1] = (lyz2prodCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binLYZ2PRODPOI-1] = (lyz2prodCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binLYZ2PRODPOI-1] = (lyz2prodCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
  
  // LYZEP = Lee-Yang Zeros Event Plane
  TList *pListLYZEP = NULL;
  AliFlowCommonHist *lyzepCommonHist = NULL;
  AliFlowCommonHistResults *lyzepCommonHistRes = NULL; 
  if(fileLYZEP) {
    fileLYZEP->GetObject("cobjLYZEP",pListLYZEP); 
    if(pListLYZEP) {
      lyzepCommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZEP->FindObject("AliFlowCommonHistLYZEP"));
      lyzepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZEP->FindObject("AliFlowCommonHistResultsLYZEP"));
      if(lyzepCommonHistRes) {
	flowValue[binLYZEP-1] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
	flowError[binLYZEP-1] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binLYZEPRP-1] = (lyzepCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binLYZEPRP-1] = (lyzepCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binLYZEPPOI-1] = (lyzepCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binLYZEPPOI-1] = (lyzepCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
 
  // FQD = Fitting q-distribution
  TList *pListFQD = NULL;
  AliFlowCommonHist *fqdCommonHist = NULL;
  AliFlowCommonHistResults *fqdCommonHistRes = NULL; 
  if(fileFQD) {
    fileFQD->GetObject("cobjFQD",pListFQD); 
    if(pListFQD) {
      fqdCommonHist = dynamic_cast<AliFlowCommonHist*> (pListFQD->FindObject("AliFlowCommonHistFQD"));
      fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListFQD->FindObject("AliFlowCommonHistResultsFQD"));
      if(fqdCommonHistRes) {
	flowValue[binFQD-1] = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);
	flowError[binFQD-1] = (fqdCommonHistRes->GetHistIntFlow())->GetBinError(1);
	//flowValueRP[binFQDRP-1] = (fqdCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
	//flowErrorRP[binFQDRP-1] = (fqdCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
	//flowValuePOI[binFQDPOI-1] = (fqdCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
	//flowErrorPOI[binFQDPOI-1] = (fqdCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
 
  // GFC = Generating Function Cumulants
  TList *pListGFC = NULL;
  AliFlowCommonHist *gfcCommonHist = NULL;
  AliFlowCommonHistResults *gfcCommonHistRes2 = NULL; 
  AliFlowCommonHistResults *gfcCommonHistRes4 = NULL; 
  AliFlowCommonHistResults *gfcCommonHistRes6 = NULL; 
  AliFlowCommonHistResults *gfcCommonHistRes8 = NULL; 
  if(fileGFC) {
    fileGFC->GetObject("cobjGFC",pListGFC);
    if(pListGFC) {
      gfcCommonHist = dynamic_cast<AliFlowCommonHist*> (pListGFC->FindObject("AliFlowCommonHistGFC"));
      gfcCommonHistRes2 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC"));
      if(gfcCommonHistRes2) {
	flowValue[binGFC2-1] = (gfcCommonHistRes2->GetHistIntFlow())->GetBinContent(1);
	flowError[binGFC2-1] = (gfcCommonHistRes2->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binGFC2RP-1] = (gfcCommonHistRes2->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binGFC2RP-1] = (gfcCommonHistRes2->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binGFC2POI-1] = (gfcCommonHistRes2->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binGFC2POI-1] = (gfcCommonHistRes2->GetHistIntFlowPOI())->GetBinError(1);
      }
      gfcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults4thOrderGFC"));
      if(gfcCommonHistRes4) {
	flowValue[binGFC4-1] = (gfcCommonHistRes4->GetHistIntFlow())->GetBinContent(1);
	flowError[binGFC4-1] = (gfcCommonHistRes4->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binGFC4RP-1] = (gfcCommonHistRes4->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binGFC4RP-1] = (gfcCommonHistRes4->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binGFC4POI-1] = (gfcCommonHistRes4->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binGFC4POI-1] = (gfcCommonHistRes4->GetHistIntFlowPOI())->GetBinError(1);
      }
      gfcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults6thOrderGFC"));
      if(gfcCommonHistRes6) {
	flowValue[binGFC6-1] = (gfcCommonHistRes6->GetHistIntFlow())->GetBinContent(1);
	flowError[binGFC6-1] = (gfcCommonHistRes6->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binGFC6RP-1] = (gfcCommonHistRes6->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binGFC6RP-1] = (gfcCommonHistRes6->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binGFC6POI-1] = (gfcCommonHistRes6->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binGFC6POI-1] = (gfcCommonHistRes6->GetHistIntFlowPOI())->GetBinError(1);
      }
      gfcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults8thOrderGFC"));
      if(gfcCommonHistRes8)  {
	flowValue[binGFC8-1] = (gfcCommonHistRes8->GetHistIntFlow())->GetBinContent(1);
	flowError[binGFC8-1] = (gfcCommonHistRes8->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binGFC8RP-1] = (gfcCommonHistRes8->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binGFC8RP-1] = (gfcCommonHistRes8->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binGFC8POI-1] = (gfcCommonHistRes8->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binGFC8POI-1] = (gfcCommonHistRes8->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }
 
  //QC = Q-cumulants
  TList *pListQC = NULL;
  AliFlowCommonHist *qcCommonHist2 = NULL; 
  AliFlowCommonHist *qcCommonHist4 = NULL; 
  AliFlowCommonHist *qcCommonHist6 = NULL; 
  AliFlowCommonHist *qcCommonHist8 = NULL; 
  AliFlowCommonHistResults *qcCommonHistRes2 = NULL; 
  AliFlowCommonHistResults *qcCommonHistRes4 = NULL; 
  AliFlowCommonHistResults *qcCommonHistRes6 = NULL; 
  AliFlowCommonHistResults *qcCommonHistRes8 = NULL; 
  
  if(fileQC) {
    fileQC->GetObject("cobjQC",pListQC);
    if(pListQC) {
      qcCommonHist2 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist2ndOrderQC"));
      qcCommonHistRes2 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults2ndOrderQC"));
      if(qcCommonHistRes2)  {
	flowValue[binQC2-1] = (qcCommonHistRes2->GetHistIntFlow())->GetBinContent(1);
	flowError[binQC2-1] = (qcCommonHistRes2->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binQC2RP-1] = (qcCommonHistRes2->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binQC2RP-1] = (qcCommonHistRes2->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binQC2POI-1] = (qcCommonHistRes2->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binQC2POI-1] = (qcCommonHistRes2->GetHistIntFlowPOI())->GetBinError(1);
      }
      qcCommonHist4 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist4thOrderQC"));
      qcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
      if(qcCommonHistRes4) {
	flowValue[binQC4-1] = (qcCommonHistRes4->GetHistIntFlow())->GetBinContent(1);
	flowError[binQC4-1] = (qcCommonHistRes4->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binQC4RP-1] = (qcCommonHistRes4->GetHistIntFlowRP())->GetBinContent(1);
	flowErrorRP[binQC4RP-1] = (qcCommonHistRes4->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binQC4POI-1] = (qcCommonHistRes4->GetHistIntFlowPOI())->GetBinContent(1);
	flowErrorPOI[binQC4POI-1] = (qcCommonHistRes4->GetHistIntFlowPOI())->GetBinError(1);
      }
      qcCommonHist6 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist6thOrderQC"));
      qcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
      if(qcCommonHistRes6) {
	flowValue[binQC6-1] = (qcCommonHistRes6->GetHistIntFlow())->GetBinContent(1);
	flowError[binQC6-1] = (qcCommonHistRes6->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binQC6RP-1] = (qcCommonHistRes6->GetHistIntFlowRP())->GetBinContent(1);
	//flowErrorRP[binQC6RP-1] = (qcCommonHistRes6->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binQC6POI-1] = (qcCommonHistRes6->GetHistIntFlowPOI())->GetBinContent(1);
	//flowErrorPOI[binQC6POI-1] = (qcCommonHistRes6->GetHistIntFlowPOI())->GetBinError(1);
      }
      qcCommonHist8 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist8thOrderQC"));
      qcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults8thOrderQC"));
      if(qcCommonHistRes8)  {
	flowValue[binQC8-1] = (qcCommonHistRes8->GetHistIntFlow())->GetBinContent(1);
	flowError[binQC8-1] = (qcCommonHistRes8->GetHistIntFlow())->GetBinError(1);
	flowValueRP[binQC8RP-1] = (qcCommonHistRes8->GetHistIntFlowRP())->GetBinContent(1);
	//flowErrorRP[binQC8RP-1] = (qcCommonHistRes8->GetHistIntFlowRP())->GetBinError(1);
	flowValuePOI[binQC8POI-1] = (qcCommonHistRes8->GetHistIntFlowPOI())->GetBinContent(1);
	//flowErrorPOI[binQC8POI-1] = (qcCommonHistRes8->GetHistIntFlowPOI())->GetBinError(1);
      }
    }
  }        
  
  // ranges on y-axis for 'no-name' plot:
  Double_t dMax=flowValue[binMC-1]+flowError[binMC-1];
  Double_t dMin=flowValue[binMC-1]-flowError[binMC-1];
  
  for(Int_t i=1;i<nMethods;i++) {
    if(!(flowValue[i]==0. && flowError[i]==0.)) {
      if(dMax<flowValue[i]+flowError[i]) dMax=flowValue[i]+flowError[i];
      if(dMin>flowValue[i]-flowError[i]) dMin=flowValue[i]-flowError[i];
    } 
  }  
  
  // ranges on y-axis for RP plot:
  Double_t dMaxRP=flowValueRP[binMCRP-1]+flowErrorRP[binMCRP-1];
  Double_t dMinRP=flowValueRP[binMCRP-1]-flowErrorRP[binMCRP-1];
  
  for(Int_t i=1;i<nMethods;i++) {
    if(!(flowValueRP[i]==0. && flowErrorRP[i]==0.)) {
      if(dMaxRP<flowValueRP[i]+flowErrorRP[i]) dMaxRP=flowValueRP[i]+flowErrorRP[i];
      if(dMinRP>flowValueRP[i]-flowErrorRP[i]) dMinRP=flowValueRP[i]-flowErrorRP[i];
    } 
  }  

  // ranges on y-axis for POI plot:
  Double_t dMaxPOI=flowValuePOI[binMCPOI-1]+flowErrorPOI[binMCPOI-1];
  Double_t dMinPOI=flowValuePOI[binMCPOI-1]-flowErrorPOI[binMCPOI-1];
  
  for(Int_t i=1;i<nMethods;i++) {
    if(!(flowValuePOI[i]==0. && flowErrorPOI[i]==0.)) {
      if(dMaxPOI<flowValuePOI[i]+flowErrorPOI[i]) dMaxPOI=flowValuePOI[i]+flowErrorPOI[i];
      if(dMinPOI>flowValuePOI[i]-flowErrorPOI[i]) dMinPOI=flowValuePOI[i]-flowErrorPOI[i];
    } 
  }  
 
  // no-name:
  TGraph* flowResults = new TGraphErrors(nMethods, x, flowValue, xError, flowError);
  flowResults->SetMarkerStyle(markerStyle);
  flowResults->SetMarkerColor(markerColor);
  
  // RP:
  TGraph* flowResultsRP = new TGraphErrors(nMethods, x, flowValueRP, xError, flowErrorRP); 
  flowResultsRP->SetMarkerStyle(markerStyleRP);
  flowResultsRP->SetMarkerColor(markerColorRP);
  
  // POI:
  TGraph* flowResultsPOI = new TGraphErrors(nMethods, x, flowValuePOI, xError, flowErrorPOI);
  flowResultsPOI->SetMarkerStyle(markerStylePOI);
  flowResultsPOI->SetMarkerColor(markerColorPOI);
  
  // for plot |v{MC}-v{method}/v{MC}|
  // no-name, RP and POI
  Double_t flowRelativeToMC[nMethods] = {0.};
  Double_t flowRelativeToMCError[nMethods] = {0.}; 
  Double_t flowRelativeToMCRP[nMethods] = {0.};
  Double_t flowRelativeToMCErrorRP[nMethods] = {0.}; 
  Double_t flowRelativeToMCPOI[nMethods] = {0.};
  Double_t flowRelativeToMCErrorPOI[nMethods] = {0.}; 
  flowRelativeToMC[0] = 0.; // MC relative to itself (to be improved)
  flowRelativeToMCRP[0] = 0.; // MC relative to itself (to be improved)
  flowRelativeToMCPOI[0] = 0.; // MC relative to itself (to be improved)
  for(Int_t i=1;i<nMethods;i++) 
  {
   if(flowValue[0] != 0)
   {
    if(flowValue[i] != 0) flowRelativeToMC[i] = (flowValue[i]-flowValue[0])/(flowValue[0]);
    if(flowError[i] != 0) flowRelativeToMCError[i] = flowError[i]/flowValue[0];
   }
   if(flowValueRP[0] != 0)
   {
    if(flowValueRP[i] != 0) flowRelativeToMCRP[i] = (flowValueRP[i]-flowValueRP[0])/(flowValueRP[0]);
    if(flowErrorRP[i] != 0) flowRelativeToMCErrorRP[i] = flowErrorRP[i]/flowValueRP[0];
   }
   if(flowValuePOI[0] != 0)
   {
    if(flowValuePOI[i] != 0) flowRelativeToMCPOI[i] = (flowValuePOI[i]-flowValuePOI[0])/(flowValuePOI[0]);
    if(flowErrorPOI[i] != 0) flowRelativeToMCErrorPOI[i] = flowErrorPOI[i]/flowValuePOI[0];
   }
  } // for(Int_t i=1;i<nMethods;i++) 
  
  // integrated flow (no-name) relative to MC:
  TGraph* flowResultsRelativeToMC = new TGraphErrors(nMethods, x, flowRelativeToMC, xError, flowRelativeToMCError);
  flowResultsRelativeToMC->SetMarkerStyle(markerStyle);
  flowResultsRelativeToMC->SetMarkerColor(markerColor);
  
  // integrated flow (RP) relative to MC:
  TGraph* flowResultsRelativeToMCRP = new TGraphErrors(nMethods, x, flowRelativeToMCRP, xError, flowRelativeToMCErrorRP);
  flowResultsRelativeToMCRP->SetMarkerStyle(markerStyleRP);
  flowResultsRelativeToMCRP->SetMarkerColor(markerColorRP);
 
  // integrated flow (POI) relative to MC:
  TGraph* flowResultsRelativeToMCPOI = new TGraphErrors(nMethods, x, flowRelativeToMCPOI, xError, flowRelativeToMCErrorPOI);
  flowResultsRelativeToMCPOI->SetMarkerStyle(markerStylePOI);
  flowResultsRelativeToMCPOI->SetMarkerColor(markerColorPOI);
  //-----------------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------------
  // cosmetics: mesh for MC error bands (integrated flow)
  TGraph* pMesh = NULL;
  TGraph* pMeshRP = NULL;
  TGraph* pMeshPOI = NULL;
  
  if(intFlowAll && mcepCommonHistRes) {
    Int_t nPts       = nMethods;
    Double_t valueMC = flowValue[binMC-1];
    Double_t errorMC = flowError[binMC-1]; 
    Double_t valueMCRP = flowValueRP[binMCRP-1];
    Double_t errorMCRP = flowErrorRP[binMCRP-1]; 
    Double_t valueMCPOI = flowValuePOI[binMCPOI-1];
    Double_t errorMCPOI = flowErrorPOI[binMCPOI-1]; 
    
    pMesh = new TGraph(nPts);
    pMeshRP = new TGraph(nPts);
    pMeshPOI = new TGraph(nPts);
    
    // no-name:
    pMesh->SetPoint(1,0,valueMC+errorMC);
    pMesh->SetPoint(2,nPts+1,valueMC+errorMC);
    pMesh->SetPoint(3,nPts+1,valueMC-errorMC);
    pMesh->SetPoint(4,0,valueMC-errorMC);
    pMesh->SetPoint(5,0,valueMC+errorMC);    
    pMesh->SetFillStyle(meshStyle);
    pMesh->SetFillColor(meshColor);
     
    // RP:
    pMeshRP->SetPoint(1,0,valueMCRP+errorMCRP);
    pMeshRP->SetPoint(2,nPts+1,valueMCRP+errorMCRP);
    pMeshRP->SetPoint(3,nPts+1,valueMCRP-errorMCRP);
    pMeshRP->SetPoint(4,0,valueMCRP-errorMCRP);
    pMeshRP->SetPoint(5,0,valueMCRP+errorMCRP);   
    pMeshRP->SetFillStyle(meshStyleRP);
    pMeshRP->SetFillColor(meshColorRP);
    
    // POI:
    pMeshPOI->SetPoint(1,0,valueMCPOI+errorMCPOI);
    pMeshPOI->SetPoint(2,nPts+1,valueMCPOI+errorMCPOI);
    pMeshPOI->SetPoint(3,nPts+1,valueMCPOI-errorMCPOI);
    pMeshPOI->SetPoint(4,0,valueMCPOI-errorMCPOI);
    pMeshPOI->SetPoint(5,0,valueMCPOI+errorMCPOI);    
    pMeshPOI->SetFillStyle(meshStylePOI);
    pMeshPOI->SetFillColor(meshColorPOI);     
  }
  //---------------------------------------------------------------------------------- 
 
  
  //----------------------------------------------------------------------------------
  //cosmetics: text (integrated flow) 
  //default text:
  TPaveText *textDefault = new TPaveText(0.05,0.77,0.95,0.90,"NDC");
  textDefault->SetTextFont(72);
  textDefault->SetTextSize(0.08);
  
  TString *entryDefaultAvM = new TString("Average Multiplicity");
  TString *entryDefaultAnd = new TString("and"); 
  TString *entryDefaultNumOfEvts = new TString("Number of Events:");
  
  textDefault->AddText(entryDefaultAvM->Data());
  textDefault->AddText(entryDefaultAnd->Data());
  textDefault->AddText(entryDefaultNumOfEvts->Data());
  
  // results (no-name):
  TPaveText *textResults = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
  textResults->SetTextFont(72);
  textResults->SetTextSize(0.06);
  
  // results (RP):
  TPaveText *textResultsRP = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
  textResultsRP->SetTextFont(72);
  textResultsRP->SetTextSize(0.06);
  
  // results (POI):
  TPaveText *textResultsPOI = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
  textResultsPOI->SetTextFont(72);
  textResultsPOI->SetTextSize(0.06);
  
  // no-name:              
  TString *entryMC       = new TString("MC ........ ");
  TString *entrySP       = new TString("SP ........ "); 
  TString *entryGFC      = new TString("GFC ....... "); 
  TString *entryQC2      = new TString("QC{2} ..... ");
  TString *entryQC4      = new TString("QC{4} ..... ");
  TString *entryQC6      = new TString("QC{6} ..... ");
  TString *entryQC8      = new TString("QC{8} ..... ");
  TString *entryFQD      = new TString("FQD ....... "); 
  TString *entryLYZ1SUM  = new TString("LYZ{sum} .. "); 
  TString *entryLYZ1PROD = new TString("LYZ{prod} . "); 
  TString *entryLYZEP    = new TString("LYZEP ..... "); 
  
  // RP: 
  TString *entryMCRP       = new TString("MC ........ ");
  TString *entrySPRP       = new TString("SP ........ ");
  TString *entryGFCRP      = new TString("GFC ....... "); 
  TString *entryQC2RP      = new TString("QC{2} ..... ");
  TString *entryQC4RP      = new TString("QC{4} ..... ");
  TString *entryQC6RP      = new TString("QC{6} ..... ");
  TString *entryQC8RP      = new TString("QC{8} ..... ");
  TString *entryFQDRP      = new TString("FQD ....... "); 
  TString *entryLYZ1SUMRP  = new TString("LYZ{sum} .. "); 
  TString *entryLYZ1PRODRP = new TString("LYZ{prod} . ");
  TString *entryLYZEPRP    = new TString("LYZEP ..... "); 
  
  // POI: 
  TString *entryMCPOI       = new TString("MC ........ ");
  TString *entrySPPOI       = new TString("SP ........ ");
  TString *entryGFCPOI      = new TString("GFC ....... "); 
  TString *entryQC2POI      = new TString("QC{2} ..... ");
  TString *entryQC4POI      = new TString("QC{4} ..... ");
  TString *entryQC6POI      = new TString("QC{6} ..... ");
  TString *entryQC8POI      = new TString("QC{8} ..... ");
  TString *entryFQDPOI      = new TString("FQD ....... "); 
  TString *entryLYZ1SUMPOI  = new TString("LYZ{sum} .. "); 
  TString *entryLYZ1PRODPOI = new TString("LYZ{prod} . "); 
  TString *entryLYZEPPOI    = new TString("LYZEP ..... "); 
  
  // no-name:
  Double_t avMultMC=0.;
  Long_t nEvtsMC=0;
  Double_t avMultSP=0.;
  Long_t nEvtsSP=0;
  Double_t avMultGFC=0.;
  Long_t nEvtsGFC=0;
  Double_t avMultQC2=0., avMultQC4=0., avMultQC6=0., avMultQC8=0.;
  Long_t nEvtsQC2=0, nEvtsQC4=0, nEvtsQC6=0, nEvtsQC8=0;
  Double_t avMultFQD=0.;
  Long_t nEvtsFQD=0;
  Double_t avMultLYZ1SUM=0.;
  Long_t nEvtsLYZ1SUM=0;
  Double_t avMultLYZ1PROD=0.;
  Long_t nEvtsLYZ1PROD=0;
  Double_t avMultLYZEP=0.;
  Long_t nEvtsLYZEP=0;
  
  // RP:
  Double_t avMultMCRP=0.;
  Long_t nEvtsMCRP=0;
  Double_t avMultSPRP=0.;
  Long_t nEvtsSPRP=0;
  Double_t avMultGFCRP=0.;
  Long_t nEvtsGFCRP=0;
  Double_t avMultQC2RP=0., avMultQC4RP=0., avMultQC6RP=0., avMultQC8RP=0.;
  Long_t nEvtsQC2RP=0, nEvtsQC4RP=0, nEvtsQC6RP=0, nEvtsQC8RP=0;
  Double_t avMultFQDRP=0.;
  Long_t nEvtsFQDRP=0;
  Double_t avMultLYZ1SUMRP=0.;
  Long_t nEvtsLYZ1SUMRP=0;
  Double_t avMultLYZ1PRODRP=0.;
  Long_t nEvtsLYZ1PRODRP=0;
  Double_t avMultLYZEPRP=0.;
  Long_t nEvtsLYZEPRP=0;
  
  // POI:
  Double_t avMultMCPOI=0.;
  Long_t nEvtsMCPOI=0;
  Double_t avMultSPPOI=0.;
  Long_t nEvtsSPPOI=0;
  Double_t avMultGFCPOI=0.;
  Long_t nEvtsGFCPOI=0;
  Double_t avMultQC2POI=0., avMultQC4POI=0., avMultQC6POI=0., avMultQC8POI=0.;
  Long_t nEvtsQC2POI=0, nEvtsQC4POI=0, nEvtsQC6POI=0, nEvtsQC8POI=0;
  Double_t avMultFQDPOI=0.;
  Long_t nEvtsFQDPOI=0;
  Double_t avMultLYZ1SUMPOI=0.;
  Long_t nEvtsLYZ1SUMPOI=0;
  Double_t avMultLYZ1PRODPOI=0.;
  Long_t nEvtsLYZ1PRODPOI=0;
  Double_t avMultLYZEPPOI=0.;
  Long_t nEvtsLYZEPPOI=0;
  
  // MC:  
  if(mcepCommonHist) {
    avMultMC = (mcepCommonHist->GetHistMultRP())->GetMean();
    nEvtsMC  = (mcepCommonHist->GetHistMultRP())->GetEntries();
    avMultMCRP = (mcepCommonHist->GetHistMultRP())->GetMean();
    nEvtsMCRP  = (mcepCommonHist->GetHistMultRP())->GetEntries();
    avMultMCPOI = (mcepCommonHist->GetHistMultPOI())->GetMean();
    nEvtsMCPOI  = (mcepCommonHist->GetHistMultPOI())->GetEntries();
  }
  
  if(entryMC) {   
    entryMC->Append("M = ");
    (*entryMC)+=(Long_t)avMultMC;
    entryMC->Append(", N = ");
    (*entryMC)+=(Long_t)nEvtsMC;
  }
 
  if(entryMCRP) {   
    entryMCRP->Append("M = ");
    (*entryMCRP)+=(Long_t)avMultMCRP;
    entryMCRP->Append(", N = ");
    (*entryMCRP)+=(Long_t)nEvtsMCRP;
  }
  
  if(entryMCPOI) {   
   entryMCPOI->Append("M = ");
   (*entryMCPOI)+=(Long_t)avMultMCPOI;
   entryMCPOI->Append(", N = ");
   (*entryMCPOI)+=(Long_t)nEvtsMCPOI;
  }
  
  // SP:  
  if(spCommonHist) {
    avMultSP = (spCommonHist->GetHistMultRP())->GetMean();
    nEvtsSP  = (spCommonHist->GetHistMultRP())->GetEntries();
    avMultSPRP = (spCommonHist->GetHistMultRP())->GetMean();
    nEvtsSPRP  = (spCommonHist->GetHistMultRP())->GetEntries();
    avMultSPPOI = (spCommonHist->GetHistMultPOI())->GetMean();
    nEvtsSPPOI  = (spCommonHist->GetHistMultPOI())->GetEntries();
  }
  
  if(entrySP) {   
    entrySP->Append("M = ");
    (*entrySP)+=(Long_t)avMultSP;
    entrySP->Append(", N = ");
    (*entrySP)+=(Long_t)nEvtsSP;
  }
 
  if(entrySPRP) {   
    entrySPRP->Append("M = ");
    (*entrySPRP)+=(Long_t)avMultSPRP;
    entrySPRP->Append(", N = ");
    (*entrySPRP)+=(Long_t)nEvtsSPRP;
  }
  
  if(entrySPPOI) {   
   entrySPPOI->Append("M = ");
   (*entrySPPOI)+=(Long_t)avMultSPPOI;
   entrySPPOI->Append(", N = ");
   (*entrySPPOI)+=(Long_t)nEvtsSPPOI;
  }
  
 // GFC:
 if(gfcCommonHist) {
   avMultGFC = (gfcCommonHist->GetHistMultRP())->GetMean();
   nEvtsGFC  = (gfcCommonHist->GetHistMultRP())->GetEntries();
   avMultGFCRP = (gfcCommonHist->GetHistMultRP())->GetMean();
   nEvtsGFCRP  = (gfcCommonHist->GetHistMultRP())->GetEntries();
   avMultGFCPOI = (gfcCommonHist->GetHistMultPOI())->GetMean();
   nEvtsGFCPOI  = (gfcCommonHist->GetHistMultPOI())->GetEntries();
 }
 
 if(entryGFC) { 
   entryGFC->Append("M = ");
   (*entryGFC)+=(Long_t)avMultGFC;
   entryGFC->Append(", N = ");
   (*entryGFC)+=(Long_t)nEvtsGFC;
 }
 
 if(entryGFCRP) { 
   entryGFCRP->Append("M = ");
   (*entryGFCRP)+=(Long_t)avMultGFCRP;
   entryGFCRP->Append(", N = ");
   (*entryGFCRP)+=(Long_t)nEvtsGFCRP;
 }
 if(entryGFCPOI) { 
   entryGFCPOI->Append("M = ");
   (*entryGFCPOI)+=(Long_t)avMultGFCPOI;
   entryGFCPOI->Append(", N = ");
   (*entryGFCPOI)+=(Long_t)nEvtsGFCPOI;
 }
 
 // QC:
 if(qcCommonHist2) {
   avMultQC2 = (qcCommonHist2->GetHistMultRP())->GetMean();
   nEvtsQC2  = (qcCommonHist2->GetHistMultRP())->GetEntries();
   avMultQC2RP = (qcCommonHist2->GetHistMultRP())->GetMean();
   nEvtsQC2RP  = (qcCommonHist2->GetHistMultRP())->GetEntries();
   avMultQC2POI = (qcCommonHist2->GetHistMultPOI())->GetMean();
   nEvtsQC2POI  = (qcCommonHist2->GetHistMultPOI())->GetEntries();
 }
 
 if(entryQC2)
 { 
  entryQC2->Append("M = ");
  (*entryQC2)+=(Long_t)avMultQC2;
  entryQC2->Append(", N = ");
  (*entryQC2)+=(Long_t)nEvtsQC2;
 }
 
 if(entryQC2RP)
 { 
  entryQC2RP->Append("M = ");
  (*entryQC2RP)+=(Long_t)avMultQC2RP;
  entryQC2RP->Append(", N = ");
  (*entryQC2RP)+=(Long_t)nEvtsQC2RP;
 } 
 
 if(entryQC2POI)
 { 
  entryQC2POI->Append("M = ");
  (*entryQC2POI)+=(Long_t)avMultQC2POI;
  entryQC2POI->Append(", N = ");
  (*entryQC2POI)+=(Long_t)nEvtsQC2POI;
 } 

 if(qcCommonHist4)
 {
  avMultQC4 = (qcCommonHist4->GetHistMultRP())->GetMean();
  nEvtsQC4  = (qcCommonHist4->GetHistMultRP())->GetEntries();
  avMultQC4RP = (qcCommonHist4->GetHistMultRP())->GetMean();
  nEvtsQC4RP  = (qcCommonHist4->GetHistMultRP())->GetEntries();
  avMultQC4POI = (qcCommonHist4->GetHistMultPOI())->GetMean();
  nEvtsQC4POI  = (qcCommonHist4->GetHistMultPOI())->GetEntries();
 }
 
 if(entryQC4)
 {
  entryQC4->Append("M = ");
  (*entryQC4)+=(Long_t)avMultQC4;
  entryQC4->Append(", N = ");
  (*entryQC4)+=(Long_t)nEvtsQC4;
 }
 
 if(entryQC4RP)
 {
  entryQC4RP->Append("M = ");
  (*entryQC4RP)+=(Long_t)avMultQC4RP;
  entryQC4RP->Append(", N = ");
  (*entryQC4RP)+=(Long_t)nEvtsQC4RP;
 }
 
 if(entryQC4POI)
 {
  entryQC4POI->Append("M = ");
  (*entryQC4POI)+=(Long_t)avMultQC4POI;
  entryQC4POI->Append(", N = ");
  (*entryQC4POI)+=(Long_t)nEvtsQC4POI;
 }
   
 if(qcCommonHist6)
 {
  avMultQC6 = (qcCommonHist6->GetHistMultRP())->GetMean();
  nEvtsQC6  = (qcCommonHist6->GetHistMultRP())->GetEntries();
  avMultQC6RP = (qcCommonHist6->GetHistMultRP())->GetMean();
  nEvtsQC6RP  = (qcCommonHist6->GetHistMultRP())->GetEntries();
  avMultQC6POI = (qcCommonHist6->GetHistMultPOI())->GetMean();
  nEvtsQC6POI  = (qcCommonHist6->GetHistMultPOI())->GetEntries();
 }
 
 if(entryQC6)
 {  
  entryQC6->Append("M = ");
  (*entryQC6)+=(Long_t)avMultQC6;
  entryQC6->Append(", N = ");
  (*entryQC6)+=(Long_t)nEvtsQC6;
 }
 
 if(entryQC6RP)
 {  
  entryQC6RP->Append("M = ");
  (*entryQC6RP)+=(Long_t)avMultQC6RP;
  entryQC6RP->Append(", N = ");
  (*entryQC6RP)+=(Long_t)nEvtsQC6RP;
 }
 
 if(entryQC6POI)
 {  
  entryQC6POI->Append("M = ");
  (*entryQC6POI)+=(Long_t)avMultQC6POI;
  entryQC6POI->Append(", N = ");
  (*entryQC6POI)+=(Long_t)nEvtsQC6POI;
 }
   
 if(qcCommonHist8)
 {
  avMultQC8 = (qcCommonHist8->GetHistMultRP())->GetMean();
  nEvtsQC8  = (qcCommonHist8->GetHistMultRP())->GetEntries();
  avMultQC8RP = (qcCommonHist8->GetHistMultRP())->GetMean();
  nEvtsQC8RP  = (qcCommonHist8->GetHistMultRP())->GetEntries();
  avMultQC8POI = (qcCommonHist8->GetHistMultPOI())->GetMean();
  nEvtsQC8POI  = (qcCommonHist8->GetHistMultPOI())->GetEntries();    
 }
  
 if(entryQC8)
 {
  entryQC8->Append("M = ");
  (*entryQC8)+=(Long_t)avMultQC8;
  entryQC8->Append(", N = ");
  (*entryQC8)+=(Long_t)nEvtsQC8;
 }
 
 if(entryQC8RP)
 {
  entryQC8RP->Append("M = ");
  (*entryQC8RP)+=(Long_t)avMultQC8RP;
  entryQC8RP->Append(", N = ");
  (*entryQC8RP)+=(Long_t)nEvtsQC8RP;
 }
 
 if(entryQC8POI)
 {
  entryQC8POI->Append("M = ");
  (*entryQC8POI)+=(Long_t)avMultQC8POI;
  entryQC8POI->Append(", N = ");
  (*entryQC8POI)+=(Long_t)nEvtsQC8POI;
 }
  
 // FQD:
 if(fqdCommonHist)
 {
  avMultFQD = (fqdCommonHist->GetHistMultRP())->GetMean();
  nEvtsFQD  = (fqdCommonHist->GetHistMultRP())->GetEntries();
  avMultFQDRP = (fqdCommonHist->GetHistMultRP())->GetMean();
  nEvtsFQDRP  = (fqdCommonHist->GetHistMultRP())->GetEntries();
  avMultFQDPOI = (fqdCommonHist->GetHistMultPOI())->GetMean();
  nEvtsFQDPOI  = (fqdCommonHist->GetHistMultPOI())->GetEntries();
 } 
 
 if(entryFQD)
 {
  entryFQD->Append("M = ");
  (*entryFQD)+=(Long_t)avMultFQD;
  entryFQD->Append(", N = ");
  (*entryFQD)+=(Long_t)nEvtsFQD;
 }
 
 if(entryFQDRP)
 {
  entryFQDRP->Append("M = ");
  (*entryFQDRP)+=(Long_t)avMultFQDRP;
  entryFQDRP->Append(", N = ");
  (*entryFQDRP)+=(Long_t)nEvtsFQDRP;
 }
 
 if(entryFQDPOI)
 {
  entryFQDPOI->Append("M = ");
  (*entryFQDPOI)+=(Long_t)avMultFQDPOI;
  entryFQDPOI->Append(", N = ");
  (*entryFQDPOI)+=(Long_t)nEvtsFQDPOI;
 }  
  
 // LYZ1SUM:
 if(lyz1sumCommonHist)
 {
  avMultLYZ1SUM = (lyz1sumCommonHist->GetHistMultRP())->GetMean();
  nEvtsLYZ1SUM  = (lyz1sumCommonHist->GetHistMultRP())->GetEntries();
  avMultLYZ1SUMRP = (lyz1sumCommonHist->GetHistMultRP())->GetMean();
  nEvtsLYZ1SUMRP  = (lyz1sumCommonHist->GetHistMultRP())->GetEntries();
  avMultLYZ1SUMPOI = (lyz1sumCommonHist->GetHistMultPOI())->GetMean();
  nEvtsLYZ1SUMPOI  = (lyz1sumCommonHist->GetHistMultPOI())->GetEntries();
 }
 
 if(entryLYZ1SUM) 
 {
  entryLYZ1SUM->Append("M = ");
  (*entryLYZ1SUM)+=(Long_t)avMultLYZ1SUM;
  entryLYZ1SUM->Append(", N = ");
  (*entryLYZ1SUM)+=(Long_t)nEvtsLYZ1SUM;
 }
 
 if(entryLYZ1SUMRP) 
 {
  entryLYZ1SUMRP->Append("M = ");
  (*entryLYZ1SUMRP)+=(Long_t)avMultLYZ1SUMRP;
  entryLYZ1SUMRP->Append(", N = ");
  (*entryLYZ1SUMRP)+=(Long_t)nEvtsLYZ1SUMRP;
 }
 
 if(entryLYZ1SUMPOI) 
 {
  entryLYZ1SUMPOI->Append("M = ");
  (*entryLYZ1SUMPOI)+=(Long_t)avMultLYZ1SUMPOI;
  entryLYZ1SUMPOI->Append(", N = ");
  (*entryLYZ1SUMPOI)+=(Long_t)nEvtsLYZ1SUMPOI;
 }
 
 // LYZ1PROD:
 if(lyz1prodCommonHist)
 {
  avMultLYZ1PROD = (lyz1prodCommonHist->GetHistMultRP())->GetMean();
  nEvtsLYZ1PROD  = (lyz1prodCommonHist->GetHistMultRP())->GetEntries();
  avMultLYZ1PRODRP = (lyz1prodCommonHist->GetHistMultRP())->GetMean();
  nEvtsLYZ1PRODRP  = (lyz1prodCommonHist->GetHistMultRP())->GetEntries();
  avMultLYZ1PRODPOI = (lyz1prodCommonHist->GetHistMultPOI())->GetMean();
  nEvtsLYZ1PRODPOI  = (lyz1prodCommonHist->GetHistMultPOI())->GetEntries();
 }
  
 if(entryLYZ1PROD) 
 {
  entryLYZ1PROD->Append("M = ");
  (*entryLYZ1PROD)+=(Long_t)avMultLYZ1PROD;
  entryLYZ1PROD->Append(", N = ");
  (*entryLYZ1PROD)+=(Long_t)nEvtsLYZ1PROD;
 }
 
 if(entryLYZ1PRODRP) 
 {
  entryLYZ1PRODRP->Append("M = ");
  (*entryLYZ1PRODRP)+=(Long_t)avMultLYZ1PRODRP;
  entryLYZ1PRODRP->Append(", N = ");
  (*entryLYZ1PRODRP)+=(Long_t)nEvtsLYZ1PRODRP;
 }
 
 if(entryLYZ1PRODPOI) 
 {
  entryLYZ1PRODPOI->Append("M = ");
  (*entryLYZ1PRODPOI)+=(Long_t)avMultLYZ1PRODPOI;
  entryLYZ1PRODPOI->Append(", N = ");
  (*entryLYZ1PRODPOI)+=(Long_t)nEvtsLYZ1PRODPOI;
 } 
 
 // LYZEP:
 if(lyzepCommonHist)
 {
  avMultLYZEP = (lyzepCommonHist->GetHistMultRP())->GetMean();
  nEvtsLYZEP  = (lyzepCommonHist->GetHistMultRP())->GetEntries();
  avMultLYZEPRP = (lyzepCommonHist->GetHistMultRP())->GetMean();
  nEvtsLYZEPRP  = (lyzepCommonHist->GetHistMultRP())->GetEntries();
  avMultLYZEPPOI = (lyzepCommonHist->GetHistMultPOI())->GetMean();
  nEvtsLYZEPPOI  = (lyzepCommonHist->GetHistMultPOI())->GetEntries();    
 }
 
 
 if(entryLYZEP)
 {
  entryLYZEP->Append("M = ");
  (*entryLYZEP)+=(Long_t)avMultLYZEP;
  entryLYZEP->Append(", N = ");
  (*entryLYZEP)+=(Long_t)nEvtsLYZEP;
 }
 
 if(entryLYZEPRP)
 {
  entryLYZEPRP->Append("M = ");
  (*entryLYZEPRP)+=(Long_t)avMultLYZEPRP;
  entryLYZEPRP->Append(", N = ");
  (*entryLYZEPRP)+=(Long_t)nEvtsLYZEPRP;
 }
 
 if(entryLYZEPPOI)
 {
  entryLYZEPPOI->Append("M = ");
  (*entryLYZEPPOI)+=(Long_t)avMultLYZEPPOI;
  entryLYZEPPOI->Append(", N = ");
  (*entryLYZEPPOI)+=(Long_t)nEvtsLYZEPPOI;
 }

 // no-name:
 if(textResults)
 {
  textResults->AddText(entryMC->Data());
  textResults->AddText(entrySP->Data());
  textResults->AddText(entryGFC->Data());
  textResults->AddText(entryQC2->Data());
  textResults->AddText(entryQC4->Data());
  textResults->AddText(entryQC6->Data());
  textResults->AddText(entryQC8->Data());
  textResults->AddText(entryFQD->Data());
  textResults->AddText(entryLYZ1SUM->Data());
  textResults->AddText(entryLYZ1PROD->Data());
  textResults->AddText(entryLYZEP->Data());
 }
 
 // RP:
 if(textResultsRP)
 {
  textResultsRP->AddText(entryMCRP->Data());
  textResultsRP->AddText(entrySPRP->Data());
  textResultsRP->AddText(entryGFCRP->Data());
  textResultsRP->AddText(entryQC2RP->Data());
  textResultsRP->AddText(entryQC4RP->Data());
  textResultsRP->AddText(entryQC6RP->Data());
  textResultsRP->AddText(entryQC8RP->Data());
  textResultsRP->AddText(entryFQDRP->Data());
  textResultsRP->AddText(entryLYZ1SUMRP->Data());
  textResultsRP->AddText(entryLYZ1PRODRP->Data());
  textResultsRP->AddText(entryLYZEPRP->Data());
 }
 
 // POI:
 if(textResultsPOI)
 {
  textResultsPOI->AddText(entryMCPOI->Data());
  textResultsPOI->AddText(entrySPPOI->Data());
  textResultsPOI->AddText(entryGFCPOI->Data());
  textResultsPOI->AddText(entryQC2POI->Data());
  textResultsPOI->AddText(entryQC4POI->Data());
  textResultsPOI->AddText(entryQC6POI->Data());
  textResultsPOI->AddText(entryQC8POI->Data());
  textResultsPOI->AddText(entryFQDPOI->Data());
  textResultsPOI->AddText(entryLYZ1SUMPOI->Data());
  textResultsPOI->AddText(entryLYZ1PRODPOI->Data());
  textResultsPOI->AddText(entryLYZEPPOI->Data());
 }
 //----------------------------------------------------------------------------------
 
 
 //----------------------------------------------------------------------------------
 // final drawing for integrated flow (no-name):
 if(plotIntFlow)
 {
  TCanvas* intFlowAllCanvas = new TCanvas("Integrated Flow","Integrated Flow",1000,600);
 
  if(plotLegendIntFlow) 
  {
   intFlowAllCanvas->Divide(2,1);
   // 1st pad is for plot:
   (intFlowAllCanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  } 
    
  if(intFlowAll)
  {
   if(dMin>0. && dMax>0.)
   {
    (intFlowAll->GetYaxis())->SetRangeUser(0.9744*dMin,1.0144*dMax);
   } else if(dMin<0. && dMax>0.)
     {
      if(!(-1.*dMin<4.*dMax))
      {  
       (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMin,1.0144*dMax);
      } else {(intFlowAll->GetYaxis())->SetRangeUser(1.1266*dMin,1.0144*dMax);}  
     } else if(dMin<0. && dMax<0.)
       {
        (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMin,0.9866*dMax);      
       }
   intFlowAll->Draw("E1");
  }                    
                                                    
  if(pMesh) pMesh->Draw("LFSAME");
   
  if(flowResults) flowResults->Draw("PSAME");

  // 2nd pad is for legend:
  if(plotLegendIntFlow)
  {
   (intFlowAllCanvas->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
   if(textDefault) textDefault->Draw();
 
   if(textResults) textResults->Draw();
  }
  
 }// end of if(plotIntFlow) 
 //----------------------------------------------------------------------------------
 
 
 //----------------------------------------------------------------------------------
 // final drawing for integrated flow relative to MC (no-name):
 if(plotIntFlowRelativeToMC)
 {
  TCanvas* intFlowAllRelativeToMCCanvas = new TCanvas("Integrated Flow Relative To MC","Integrated Flow Relative To MC",1000,600);
 
  intFlowAllRelativeToMCCanvas->Divide(2,1);
 
  // 1st pad is for plot:
  (intFlowAllRelativeToMCCanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  
  TH1D *intFlowAllRelativeToMC = new TH1D(*intFlowAll);
  (intFlowAllRelativeToMC->GetYaxis())->SetRangeUser(-1,1);
  (intFlowAllRelativeToMC->GetYaxis())->SetTitle("(v_{n}\{method\} - v_{n}\{MC\})/v_{n}\{MC\}"); 
  intFlowAllRelativeToMC->Draw("E1");            
  
  if(flowResultsRelativeToMC) flowResultsRelativeToMC->Draw("PSAME");

  // 2nd pad is for legend:
  (intFlowAllRelativeToMCCanvas->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
  if(textDefault) textDefault->Draw();
 
  if(textResults) textResults->Draw();
 
 }// end of if(plotIntFlowRelativeToMC) 
 //----------------------------------------------------------------------------------
 
 
 //----------------------------------------------------------------------------------
 //final drawing for integrated flow of RP (i.e. of particles used to determine the reaction plane):
 if(plotIntFlowRP)
 {
  TCanvas* intFlowAllCanvasRP = new TCanvas("Integrated Flow RP","Integrated Flow RP",1000,600);
 
  if(plotLegendIntFlow)
  {
   intFlowAllCanvasRP->Divide(2,1);
 
   //1st pad is for plot:
   (intFlowAllCanvasRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  }
 
  TH1D *intFlowAllRP = new TH1D(*intFlowAll);
  intFlowAllRP->SetMarkerStyle(markerStyleRP);
  intFlowAllRP->SetMarkerColor(markerColorRP);
  (intFlowAllRP->GetXaxis())->SetBinLabel(binMCRP,"v_{2}{MC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binSPRP,"v_{2}{SP}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binGFC2RP,"v_{2}{2,GFC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binQC2RP,"v_{2}{2,QC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binGFC4RP,"v_{2}{4,GFC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binQC4RP,"v_{2}{4,QC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binGFC6RP,"v_{2}{6,GFC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binQC6RP,"v_{2}{6,QC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binGFC8RP,"v_{2}{8,GFC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binQC8RP,"v_{2}{8,QC}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binFQDRP,"v_{2}{FQD}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binLYZ2SUMRP,"v_{2}{LYZ,sum}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binLYZ2PRODRP,"v_{2}{LYZ,prod}");
  (intFlowAllRP->GetXaxis())->SetBinLabel(binLYZEPRP,"v_{2}{LYZEP}");

  if(intFlowAllRP)
  {
   if(dMinRP>0. && dMaxRP>0.)
   {
    (intFlowAllRP->GetYaxis())->SetRangeUser(0.9744*dMinRP,1.0144*dMaxRP);
   } else if(dMinRP<0. && dMaxRP>0.)
     {
      if(!(-1.*dMinRP<4.*dMaxRP))
      {  
       (intFlowAllRP->GetYaxis())->SetRangeUser(1.0266*dMinRP,1.0144*dMaxRP);
      } else {(intFlowAllRP->GetYaxis())->SetRangeUser(1.1266*dMinRP,1.0144*dMaxRP);}  
     } else if(dMinRP<0. && dMaxRP<0.)
       {
        (intFlowAllRP->GetYaxis())->SetRangeUser(1.0266*dMinRP,0.9866*dMaxRP);      
       } 
   intFlowAllRP->Draw("E1");
  }
                                                                                                                                                                                                                                                                                   
  if(pMeshRP) pMeshRP->Draw("LFSAME");
   
  if(flowResultsRP) flowResultsRP->Draw("PSAME");

  if(plotLegendIntFlow)
  {
   //2nd pad is for legend:
   (intFlowAllCanvasRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);
  
   if(textDefault) textDefault->Draw();
 
   if(textResultsRP) textResultsRP->Draw();
  }
   
 }//end of if(plotIntFlowRP} 
 //----------------------------------------------------------------------------------
 
 
 //----------------------------------------------------------------------------------
 // final drawing for integrated flow relative to MC (RP):
 if(plotIntFlowRelativeToMCRP)
 {
  TCanvas* intFlowAllRelativeToMCRPCanvas = new TCanvas("Integrated Flow (RP) Relative To MC","Integrated Flow (RP) Relative To MC",1000,600);
 
  intFlowAllRelativeToMCRPCanvas->Divide(2,1);
 
  // 1st pad is for plot:
  (intFlowAllRelativeToMCRPCanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  
  TH1D *intFlowAllRelativeToMCRP = new TH1D(*intFlowAll);
  (intFlowAllRelativeToMCRP->GetYaxis())->SetRangeUser(-1,1); 
  (intFlowAllRelativeToMCRP->GetYaxis())->SetTitle("(v_{n}\{method\} - v_{n}\{MC\})/v_{n}\{MC\}");
  intFlowAllRelativeToMCRP->Draw("E1");            
  
  if(flowResultsRelativeToMCRP) flowResultsRelativeToMCRP->Draw("PSAME");

  // 2nd pad is for legend:
  (intFlowAllRelativeToMCRPCanvas->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
  if(textDefault) textDefault->Draw();
 
  if(textResultsRP) textResultsRP->Draw();
 
 }// end of if(plotIntFlowRelativeToMCRP) 
 //----------------------------------------------------------------------------------
  
 
 //----------------------------------------------------------------------------------
 //final drawing for integrated flow of POI (i.e. of particles of interest):
 if(plotIntFlowPOI)
 {
  TCanvas* intFlowAllCanvasPOI = new TCanvas("Integrated Flow POI","Integrated Flow POI",1000,600);
  
  if(plotLegendIntFlow)
  {
   intFlowAllCanvasPOI->Divide(2,1);
 
   //1st pad is for plot:
   (intFlowAllCanvasPOI->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  }
  
  TH1D *intFlowAllPOI = new TH1D(*intFlowAll);
  intFlowAllPOI->SetMarkerStyle(markerStylePOI);
  intFlowAllPOI->SetMarkerColor(markerColorPOI);
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binMCPOI,"v_{2}{MC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binSPPOI,"v_{2}{SP}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binGFC2POI,"v_{2}{2,GFC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binQC2POI,"v_{2}{2,QC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binGFC4POI,"v_{2}{4,GFC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binQC4POI,"v_{2}{4,QC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binGFC6POI,"v_{2}{6,GFC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binQC6POI,"v_{2}{6,QC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binGFC8POI,"v_{2}{8,GFC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binQC8POI,"v_{2}{8,QC}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binFQDPOI,"v_{2}{FQD}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binLYZ2SUMPOI,"v_{2}{LYZ,sum}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binLYZ2PRODPOI,"v_{2}{LYZ,prod}");
  (intFlowAllPOI->GetXaxis())->SetBinLabel(binLYZEPPOI,"v_{2}{LYZEP}");
  
  if(intFlowAllPOI)
  {
   if(dMinPOI>0. && dMaxPOI>0.)
   {
    (intFlowAllPOI->GetYaxis())->SetRangeUser(0.9744*dMinPOI,1.0144*dMaxPOI);
   } else if(dMinPOI<0. && dMaxPOI>0.)
     {
      if(!(-1.*dMinPOI<4.*dMaxPOI))
      {  
       (intFlowAllPOI->GetYaxis())->SetRangeUser(1.0266*dMinPOI,1.0144*dMaxPOI);
      } else {(intFlowAllPOI->GetYaxis())->SetRangeUser(1.1266*dMinPOI,1.0144*dMaxPOI);}  
     } else if(dMinPOI<0. && dMaxPOI<0.)
       {
        (intFlowAllPOI->GetYaxis())->SetRangeUser(1.0266*dMinPOI,0.9866*dMaxPOI);      
       }
   intFlowAllPOI->Draw("E1");
  }
                            
  if(pMeshPOI) pMeshPOI->Draw("LFSAME");
   
  if(flowResultsPOI) flowResultsPOI->Draw("PSAME");
 
  if(plotLegendIntFlow)
  {
   //2nd pad is for legend:
   (intFlowAllCanvasPOI->cd(2))->SetPad(0.75,0.0,1.0,1.0);
  
   if(textDefault) textDefault->Draw();

   if(textResultsPOI) textResultsPOI->Draw();
  } 
   
 }// end of if(plotIntFlowPOI) 
 //----------------------------------------------------------------------------------      
 
    
 //----------------------------------------------------------------------------------
 // final drawing for integrated flow relative to MC (POI):
 if(plotIntFlowRelativeToMCPOI)
 {
  TCanvas* intFlowAllRelativeToMCPOICanvas = new TCanvas("Integrated Flow (POI) Relative To MC","Integrated Flow (POI) Relative To MC",1000,600);
 
  intFlowAllRelativeToMCPOICanvas->Divide(2,1);
 
  // 1st pad is for plot:
  (intFlowAllRelativeToMCPOICanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  
  TH1D *intFlowAllRelativeToMCPOI = new TH1D(*intFlowAll);
  (intFlowAllRelativeToMCPOI->GetYaxis())->SetRangeUser(-1,1); 
  (intFlowAllRelativeToMCPOI->GetYaxis())->SetTitle("(v_{n}\{method\} - v_{n}\{MC\})/v_{n}\{MC\}");
  intFlowAllRelativeToMCPOI->Draw("E1");            
  
  if(flowResultsRelativeToMCPOI) flowResultsRelativeToMCPOI->Draw("PSAME");

  // 2nd pad is for legend:
  (intFlowAllRelativeToMCPOICanvas->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
  if(textDefault) textDefault->Draw();
 
  if(textResultsPOI) textResultsPOI->Draw();
 
 }// end of if(plotIntFlowRelativeToMCPOI) 
 //----------------------------------------------------------------------------------
         
 //==================================================================================   




 //==================================================================================
 //                            DIFFERENTIAL FLOW
 //==================================================================================
 Int_t iNbinsPt  = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
 Double_t dPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();
 Double_t dPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
 
 Int_t iNbinsEta  = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
 Double_t dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();
 Double_t dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();
 
 //----------------------------------------------------------------------------------
 //cosmetics: the style histogram for differential flow (pt):
 TH1D *styleHistPt = new TH1D("styleHistPt","styleHistPt",iNbinsPt,dPtMin,dPtMax);
 styleHistPt->SetTitle("Differential Flow");
 styleHistPt->SetXTitle("p_{t} [GeV]");
 styleHistPt->SetYTitle("v_{n}");
 
 //cosmetics: the style histogram for differential flow (eta):
 TH1D *styleHistEta = new TH1D("styleHistEta","styleHistEta",iNbinsEta,dEtaMin,dEtaMax);
 styleHistEta->SetTitle("Differential Flow");
 styleHistEta->SetXTitle("#eta");
 styleHistEta->SetYTitle("v_{n}");
 //----------------------------------------------------------------------------------
 
 

 //----------------------------------------------------------------------------------
 //RP:
 //cosmetics: Monte Carlo error bands for differential flow (Pt)
 TGraph* pMeshDiffFlowPtRP = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nBinsDiffFlowPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetNbinsX();
  Double_t binWidthPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinWidth(1);//assuming that all bins have the same width

  //counting the non-empty bins: 
  Int_t nNonEmptyBinsDiffFlowPtRP=0;
  for(Int_t i=1;i<nBinsDiffFlowPtRP+1;i++)
  {
   if(!(mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(i)==0.0))
   {
    nNonEmptyBinsDiffFlowPtRP++;
   }
  }    
       
  pMeshDiffFlowPtRP = new TGraph(2*nNonEmptyBinsDiffFlowPtRP+1);
  
  Double_t valueMCPtRP=0.,errorMCPtRP=0.;
  Int_t countDiffFlowPtRP=1;
  Double_t xFirstDiffFlowPtRP=0.,yUpFirstDiffFlowPtRP=0.;//needed to close up the mesh
  for(Int_t i=1;i<nBinsDiffFlowPtRP+1;i++)
  {
   //setting up the upper limit of the mesh:
   valueMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(i);
   errorMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(i);   
   if(!(errorMCPtRP==0.0))
   {    
    pMeshDiffFlowPtRP->SetPoint(countDiffFlowPtRP++,(i-0.5)*binWidthPtRP+dPtMin,valueMCPtRP+errorMCPtRP);
    if(xFirstDiffFlowPtRP==0.)
    {
     xFirstDiffFlowPtRP=(i-0.5)*binWidthPtRP+dPtMin;
     yUpFirstDiffFlowPtRP=valueMCPtRP+errorMCPtRP;
    }
   } 
  }   
  for(Int_t i=nBinsDiffFlowPtRP+1;i<2*nBinsDiffFlowPtRP+1;i++)
  {
   //setting up the lower limit of the mesh:
   valueMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(2*nBinsDiffFlowPtRP+1-i);
   errorMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(2*nBinsDiffFlowPtRP+1-i); 
   if(!(errorMCPtRP==0.0))
   {      
    pMeshDiffFlowPtRP->SetPoint(countDiffFlowPtRP++,(2*nBinsDiffFlowPtRP-i+0.5)*binWidthPtRP+dPtMin,valueMCPtRP-errorMCPtRP);
   }  
  }
  //closing the mesh area:
  pMeshDiffFlowPtRP->SetPoint(2*nNonEmptyBinsDiffFlowPtRP+1,xFirstDiffFlowPtRP,yUpFirstDiffFlowPtRP);   
  
  //setting the mesh style and color:               
  pMeshDiffFlowPtRP->SetFillStyle(meshStyleDiffFlowPtRP);
  pMeshDiffFlowPtRP->SetFillColor(meshColorDiffFlowPtRP);
 }

 //cosmetics: Monte Carlo error bands for differential flow (Eta)
 TGraph* pMeshDiffFlowEtaRP = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nBinsDiffFlowEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetNbinsX();
  Double_t binWidthEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinWidth(1);//assuming that all bins have the same width

  //counting the non-empty bins: 
  Int_t nNonEmptyBinsDiffFlowEtaRP=0;
  for(Int_t i=1;i<nBinsDiffFlowEtaRP+1;i++)
  {
   if(!(mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinError(i)==0.0))
   {
    nNonEmptyBinsDiffFlowEtaRP++;
   }
  }    
       
  pMeshDiffFlowEtaRP = new TGraph(2*nNonEmptyBinsDiffFlowEtaRP+1);
  
  Double_t valueMCEtaRP=0.,errorMCEtaRP=0.;
  Int_t countDiffFlowEtaRP=1;
  Double_t xFirstDiffFlowEtaRP=0.,yUpFirstDiffFlowEtaRP=0.;//needed to close up the mesh
  for(Int_t i=1;i<nBinsDiffFlowEtaRP+1;i++)
  {
   //setting up the upper limit of the mesh:
   valueMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinContent(i);
   errorMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinError(i);   
   if(!(errorMCEtaRP==0.0))
   {    
    pMeshDiffFlowEtaRP->SetPoint(countDiffFlowEtaRP++,(i-0.5)*binWidthEtaRP+dEtaMin,valueMCEtaRP+errorMCEtaRP);
    if(xFirstDiffFlowEtaRP==0.)
    {
     xFirstDiffFlowEtaRP=(i-0.5)*binWidthEtaRP+dEtaMin;
     yUpFirstDiffFlowEtaRP=valueMCEtaRP+errorMCEtaRP;
    }
   } 
  }   
  for(Int_t i=nBinsDiffFlowEtaRP+1;i<2*nBinsDiffFlowEtaRP+1;i++)
  {
   //setting up the lower limit of the mesh:
   valueMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinContent(2*nBinsDiffFlowEtaRP+1-i);
   errorMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinError(2*nBinsDiffFlowEtaRP+1-i); 
   if(!(errorMCEtaRP==0.0))
   {      
    pMeshDiffFlowEtaRP->SetPoint(countDiffFlowEtaRP++,(2*nBinsDiffFlowEtaRP-i+0.5)*binWidthEtaRP+dEtaMin,valueMCEtaRP-errorMCEtaRP);
   }  
  }
  //closing the mesh area:
  pMeshDiffFlowEtaRP->SetPoint(2*nNonEmptyBinsDiffFlowEtaRP+1,xFirstDiffFlowEtaRP,yUpFirstDiffFlowEtaRP);   
  
  //setting the mesh style and color:               
  pMeshDiffFlowEtaRP->SetFillStyle(meshStyleDiffFlowEtaRP);
  pMeshDiffFlowEtaRP->SetFillColor(meshColorDiffFlowEtaRP);
 } 
 //----------------------------------------------------------------------------------




 //----------------------------------------------------------------------------------
 //POI:
 //cosmetics: Monte Carlo error bands for differential flow (Pt)
 TGraph* pMeshDiffFlowPtPOI = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nBinsDiffFlowPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetNbinsX();
  Double_t binWidthPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinWidth(1);//assuming that all bins have the same width

  //counting the non-empty bins: 
  Int_t nNonEmptyBinsDiffFlowPtPOI=0;
  for(Int_t i=1;i<nBinsDiffFlowPtPOI+1;i++)
  {
   if(!(mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinError(i)==0.0))
   {
    nNonEmptyBinsDiffFlowPtPOI++;
   }
  }    
       
  pMeshDiffFlowPtPOI = new TGraph(2*nNonEmptyBinsDiffFlowPtPOI+1);
  
  Double_t valueMCPtPOI=0.,errorMCPtPOI=0.;
  Int_t countDiffFlowPtPOI=1;
  Double_t xFirstDiffFlowPtPOI=0.,yUpFirstDiffFlowPtPOI=0.;//needed to close up the mesh
  for(Int_t i=1;i<nBinsDiffFlowPtPOI+1;i++)
  {
   //setting up the upper limit of the mesh:
   valueMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinContent(i);
   errorMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinError(i);   
   if(!(errorMCPtPOI==0.0))
   {    
    pMeshDiffFlowPtPOI->SetPoint(countDiffFlowPtPOI++,(i-0.5)*binWidthPtPOI+dPtMin,valueMCPtPOI+errorMCPtPOI);
    if(xFirstDiffFlowPtPOI==0.)
    {
     xFirstDiffFlowPtPOI=(i-0.5)*binWidthPtPOI+dPtMin;
     yUpFirstDiffFlowPtPOI=valueMCPtPOI+errorMCPtPOI;
    }
   } 
  }   
  for(Int_t i=nBinsDiffFlowPtPOI+1;i<2*nBinsDiffFlowPtPOI+1;i++)
  {
   //setting up the lower limit of the mesh:
   valueMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinContent(2*nBinsDiffFlowPtPOI+1-i);
   errorMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinError(2*nBinsDiffFlowPtPOI+1-i); 
   if(!(errorMCPtPOI==0.0))
   {      
    pMeshDiffFlowPtPOI->SetPoint(countDiffFlowPtPOI++,(2*nBinsDiffFlowPtPOI-i+0.5)*binWidthPtPOI+dPtMin,valueMCPtPOI-errorMCPtPOI);
   }  
  }
  //closing the mesh area:
  pMeshDiffFlowPtPOI->SetPoint(2*nNonEmptyBinsDiffFlowPtPOI+1,xFirstDiffFlowPtPOI,yUpFirstDiffFlowPtPOI);   
  
  //setting the mesh style and color:               
  pMeshDiffFlowPtPOI->SetFillStyle(meshStyleDiffFlowPtPOI);
  pMeshDiffFlowPtPOI->SetFillColor(meshColorDiffFlowPtPOI);
 }

 //cosmetics: Monte Carlo error bands for differential flow (Eta)
 TGraph* pMeshDiffFlowEtaPOI = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nBinsDiffFlowEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetNbinsX();
  Double_t binWidthEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinWidth(1);//assuming that all bins have the same width

  //counting the non-empty bins: 
  Int_t nNonEmptyBinsDiffFlowEtaPOI=0;
  for(Int_t i=1;i<nBinsDiffFlowEtaPOI+1;i++)
  {
   if(!(mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinError(i)==0.0))
   {
    nNonEmptyBinsDiffFlowEtaPOI++;
   }
  }    
       
  pMeshDiffFlowEtaPOI = new TGraph(2*nNonEmptyBinsDiffFlowEtaPOI+1);
  
  Double_t valueMCEtaPOI=0.,errorMCEtaPOI=0.;
  Int_t countDiffFlowEtaPOI=1;
  Double_t xFirstDiffFlowEtaPOI=0.,yUpFirstDiffFlowEtaPOI=0.;//needed to close up the mesh
  for(Int_t i=1;i<nBinsDiffFlowEtaPOI+1;i++)
  {
   //setting up the upper limit of the mesh:
   valueMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinContent(i);
   errorMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinError(i);   
   if(!(errorMCEtaPOI==0.0))
   {    
    pMeshDiffFlowEtaPOI->SetPoint(countDiffFlowEtaPOI++,(i-0.5)*binWidthEtaPOI+dEtaMin,valueMCEtaPOI+errorMCEtaPOI);
    if(xFirstDiffFlowEtaPOI==0.)
    {
     xFirstDiffFlowEtaPOI=(i-0.5)*binWidthEtaPOI+dEtaMin;
     yUpFirstDiffFlowEtaPOI=valueMCEtaPOI+errorMCEtaPOI;
    }
   } 
  }   
  for(Int_t i=nBinsDiffFlowEtaPOI+1;i<2*nBinsDiffFlowEtaPOI+1;i++)
  {
   //setting up the lower limit of the mesh:
   valueMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinContent(2*nBinsDiffFlowEtaPOI+1-i);
   errorMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinError(2*nBinsDiffFlowEtaPOI+1-i); 
   if(!(errorMCEtaPOI==0.0))
   {      
    pMeshDiffFlowEtaPOI->SetPoint(countDiffFlowEtaPOI++,(2*nBinsDiffFlowEtaPOI-i+0.5)*binWidthEtaPOI+dEtaMin,valueMCEtaPOI-errorMCEtaPOI);
   }  
  }
  //closing the mesh area:
  pMeshDiffFlowEtaPOI->SetPoint(2*nNonEmptyBinsDiffFlowEtaPOI+1,xFirstDiffFlowEtaPOI,yUpFirstDiffFlowEtaPOI);   
  
  //setting the mesh style and color:               
  pMeshDiffFlowEtaPOI->SetFillStyle(meshStyleDiffFlowEtaPOI);
  pMeshDiffFlowEtaPOI->SetFillColor(meshColorDiffFlowEtaPOI);
 }
 //----------------------------------------------------------------------------------
   
 //MCEP = Monte Carlo Event Plane
 Double_t avMultDiffFlowMCRP=0.;
 Double_t nEvtsDiffFlowMCRP=0;
 Double_t avMultDiffFlowMCPOI=0.;
 Double_t nEvtsDiffFlowMCPOI=0;
 if(fileMCEP)
 {
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorMC);
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleMC);
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorMC);
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleMC);
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorMC);
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleMC);
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorMC);
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleMC);
  } 
  if(mcepCommonHist)
  {
   avMultDiffFlowMCRP = (mcepCommonHist->GetHistMultRP())->GetMean();
   nEvtsDiffFlowMCRP  = (mcepCommonHist->GetHistMultRP())->GetEntries();
   avMultDiffFlowMCPOI = (mcepCommonHist->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowMCPOI  = (mcepCommonHist->GetHistMultPOI())->GetEntries();      
  } 
 } 
 
 //SP = Scalar Product
 Double_t avMultDiffFlowSPRP=0.;
 Double_t nEvtsDiffFlowSPRP=0;
 Double_t avMultDiffFlowSPPOI=0.;
 Double_t nEvtsDiffFlowSPPOI=0;
 if(fileSP)
 {
  if(spCommonHistRes)
  {
   (spCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorSP);
   (spCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleSP);
   (spCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorSP);
   (spCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleSP);
   (spCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorSP);
   (spCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleSP);
   (spCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorSP);
   (spCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleSP);
  } 
  if(spCommonHist)
  {
   avMultDiffFlowSPRP = (spCommonHist->GetHistMultRP())->GetMean();
   nEvtsDiffFlowSPRP = (spCommonHist->GetHistMultRP())->GetEntries();
   avMultDiffFlowSPPOI = (spCommonHist->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowSPPOI = (spCommonHist->GetHistMultPOI())->GetEntries();      
  } 
 } 

 //GFC = Generating Function Cumulants
 Double_t avMultDiffFlowGFC=0.;//to be removed
 Double_t nEvtsDiffFlowGFC=0.;//to be removed
 Double_t avMultDiffFlowGFCRP=0.;
 Double_t nEvtsDiffFlowGFCRP=0.;
 Double_t avMultDiffFlowGFCPOI=0.;
 Double_t nEvtsDiffFlowGFCPOI=0.; 
 if(fileGFC)
 {
  if(gfcCommonHistRes2)
  {
   (gfcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorGFC2);
   (gfcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleGFC2);
  }
  if(gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorGFC4);
   (gfcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleGFC4);         
  }
  if(gfcCommonHistRes6)
  { 
   (gfcCommonHistRes6->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorGFC6);
   (gfcCommonHistRes6->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleGFC6);
  }
  if(gfcCommonHistRes8)
  { 
   (gfcCommonHistRes8->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorGFC8);
   (gfcCommonHistRes8->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleGFC8);
  }
  if(gfcCommonHist)
  {
   avMultDiffFlowGFCRP = (gfcCommonHist->GetHistMultRP())->GetMean();   
   nEvtsDiffFlowGFCRP  = (gfcCommonHist->GetHistMultRP())->GetEntries();
   avMultDiffFlowGFCPOI = (gfcCommonHist->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowGFCPOI  = (gfcCommonHist->GetHistMultPOI())->GetEntries();   
  } 
 }
  
 //QC = Q-cumulants
 Double_t avMultDiffFlowQC2RP=0.;
 Double_t nEvtsDiffFlowQC2RP=0.;
 Double_t avMultDiffFlowQC2POI=0.;
 Double_t nEvtsDiffFlowQC2POI=0.;
 Double_t avMultDiffFlowQC4RP=0.;
 Double_t nEvtsDiffFlowQC4RP=0.;
 Double_t avMultDiffFlowQC4POI=0.;
 Double_t nEvtsDiffFlowQC4POI=0.;
 Double_t avMultDiffFlowQC6RP=0.;
 Double_t nEvtsDiffFlowQC6RP=0.;
 Double_t avMultDiffFlowQC6POI=0.;
 Double_t nEvtsDiffFlowQC6POI=0.;
 Double_t avMultDiffFlowQC8RP=0.;
 Double_t nEvtsDiffFlowQC8RP=0.;
 Double_t avMultDiffFlowQC8POI=0.;
 Double_t nEvtsDiffFlowQC8POI=0.;

 if(fileQC)
 {
  //QC{2}
  if(qcCommonHistRes2)
  {
   (qcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorQC2);
   (qcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleQC2);
   (qcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorQC2);
   (qcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleQC2);
   (qcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorQC2);
   (qcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleQC2);
   (qcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorQC2);
   (qcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleQC2);
  }
  if(qcCommonHist2)
  {
   avMultDiffFlowQC2RP = (qcCommonHist2->GetHistMultRP())->GetMean();
   nEvtsDiffFlowQC2RP  = (qcCommonHist2->GetHistMultRP())->GetEntries();
   avMultDiffFlowQC2POI = (qcCommonHist2->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowQC2POI  = (qcCommonHist2->GetHistMultPOI())->GetEntries();
  }
  //QC{4}
  if(qcCommonHistRes4)
  {
   (qcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorQC4);
   (qcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleQC4);
   (qcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorQC4);
   (qcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleQC4);
   (qcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorQC4);
   (qcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleQC4);
   (qcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorQC4);
   (qcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleQC4);
  }
  if(qcCommonHist4)
  {
   avMultDiffFlowQC4RP = (qcCommonHist4->GetHistMultRP())->GetMean();
   nEvtsDiffFlowQC4RP  = (qcCommonHist4->GetHistMultRP())->GetEntries();
   avMultDiffFlowQC4POI = (qcCommonHist4->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowQC4POI  = (qcCommonHist4->GetHistMultPOI())->GetEntries();
  }
  //QC{6}
  if(qcCommonHistRes6)
  {
   (qcCommonHistRes6->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorQC6);
   (qcCommonHistRes6->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleQC6);
   (qcCommonHistRes6->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorQC6);
   (qcCommonHistRes6->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleQC6);
   (qcCommonHistRes6->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorQC6);
   (qcCommonHistRes6->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleQC6);
   (qcCommonHistRes6->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorQC6);
   (qcCommonHistRes6->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleQC6);
  }
  if(qcCommonHist6)
  {
   avMultDiffFlowQC6RP = (qcCommonHist6->GetHistMultRP())->GetMean();
   nEvtsDiffFlowQC6RP  = (qcCommonHist6->GetHistMultRP())->GetEntries();
   avMultDiffFlowQC6POI = (qcCommonHist6->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowQC6POI  = (qcCommonHist6->GetHistMultPOI())->GetEntries();
  }
  //QC{8}
  if(qcCommonHistRes8)
  {
   (qcCommonHistRes8->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorQC8);
   (qcCommonHistRes8->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleQC8);
   (qcCommonHistRes8->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorQC8);
   (qcCommonHistRes8->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleQC8);
   (qcCommonHistRes8->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorQC8);
   (qcCommonHistRes8->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleQC8);
   (qcCommonHistRes8->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorQC8);
   (qcCommonHistRes8->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleQC8);
  }
  if(qcCommonHist8)
  {
   avMultDiffFlowQC8RP = (qcCommonHist8->GetHistMultRP())->GetMean();
   nEvtsDiffFlowQC8RP  = (qcCommonHist8->GetHistMultRP())->GetEntries();
   avMultDiffFlowQC8POI = (qcCommonHist8->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowQC8POI  = (qcCommonHist8->GetHistMultPOI())->GetEntries();
  }
 } 

 //LYZ2SUM = Lee-Yang Zeros (2nd run, sum)
 Double_t avMultDiffFlowLYZ2SUMRP=0.;
 Double_t nEvtsDiffFlowLYZ2SUMRP=0;
 Double_t avMultDiffFlowLYZ2SUMPOI=0.;
 Double_t nEvtsDiffFlowLYZ2SUMPOI=0;
 if(fileLYZ2SUM)
 {
  if(lyz2sumCommonHistRes)
  {
   (lyz2sumCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorLYZ2SUM);
   (lyz2sumCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleLYZ2SUM);
  } 
  if(lyz2sumCommonHist)
  {
   avMultDiffFlowLYZ2SUMRP = (lyz2sumCommonHist->GetHistMultRP())->GetMean();
   nEvtsDiffFlowLYZ2SUMRP  = (lyz2sumCommonHist->GetHistMultRP())->GetEntries();
   avMultDiffFlowLYZ2SUMPOI = (lyz2sumCommonHist->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowLYZ2SUMPOI  = (lyz2sumCommonHist->GetHistMultPOI())->GetEntries();
  } 
 } 
 
 //LYZ2PROD = Lee-Yang Zeros (2nd run, product)
 Double_t avMultDiffFlowLYZ2PRODRP=0.;
 Double_t nEvtsDiffFlowLYZ2PRODRP=0;
 Double_t avMultDiffFlowLYZ2PRODPOI=0.;
 Double_t nEvtsDiffFlowLYZ2PRODPOI=0;
 if(fileLYZ2PROD)
 {
  if(lyz2prodCommonHistRes)
  {
   (lyz2prodCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorLYZ2PROD);
   (lyz2prodCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleLYZ2PROD);
  } 
  if(lyz2prodCommonHist)
  {
   avMultDiffFlowLYZ2PRODRP = (lyz2prodCommonHist->GetHistMultRP())->GetMean();
   nEvtsDiffFlowLYZ2PRODRP  = (lyz2prodCommonHist->GetHistMultRP())->GetEntries();
   avMultDiffFlowLYZ2PRODPOI = (lyz2prodCommonHist->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowLYZ2PRODPOI  = (lyz2prodCommonHist->GetHistMultPOI())->GetEntries();
  } 
 } 

 //LYZEP = Lee-Yang Zeros Event Plane
 Double_t avMultDiffFlowLYZEPRP=0.;
 Double_t nEvtsDiffFlowLYZEPRP=0;
 Double_t avMultDiffFlowLYZEPPOI=0.;
 Double_t nEvtsDiffFlowLYZEPPOI=0;
 if(fileLYZEP)
 {
  if(lyzepCommonHistRes)
  {
   (lyzepCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(markerColorLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(markerStyleLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(markerColorLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(markerStyleLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(markerColorLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(markerStyleLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(markerColorLYZEP);
   (lyzepCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(markerStyleLYZEP);
  } 
  if(lyzepCommonHist)
  {
   avMultDiffFlowLYZEPRP = (lyzepCommonHist->GetHistMultRP())->GetMean();
   nEvtsDiffFlowLYZEPRP  = (lyzepCommonHist->GetHistMultRP())->GetEntries();
   avMultDiffFlowLYZEPPOI = (lyzepCommonHist->GetHistMultPOI())->GetMean();
   nEvtsDiffFlowLYZEPPOI  = (lyzepCommonHist->GetHistMultPOI())->GetEntries();
  } 
 } 


 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Pt, RP):
 if(plotDiffFlowPtRP)
 {
  TCanvas* diffFlowPtAllCanvasRP = new TCanvas("Differential Flow (Pt) of RP","Differential Flow (Pt) of RP ",1000,600);
  
  diffFlowPtAllCanvasRP->Divide(2,1);
 
  //1st pad is for plot:
  (diffFlowPtAllCanvasRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
  if(styleHistPt)
  {
   (styleHistPt->GetYaxis())->SetRangeUser(-0.3,1.0);
   styleHistPt->Draw();
  }
  if(pMeshDiffFlowPtRP)
  {
   pMeshDiffFlowPtRP->Draw("LFSAME");
  }
 
  //MC 
  if(plotMCPtRP && mcepCommonHistRes)
  { 
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  //SP 
  if(plotSPPtRP && spCommonHistRes)
  { 
   (spCommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  //GFC
  if(plotGFC2PtRP && gfcCommonHistRes2)Pt
  { 
   (gfcCommonHistRes2->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
  } 
  if(plotGFC4PtRP && gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
  } 
  if(plotGFC6PtRP && gfcCommonHistRes6)
  { 
   (gfcCommonHistRes6->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
  } 
  if(plotGFC8PtRP && gfcCommonHistRes8)
  { 
   (gfcCommonHistRes8->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
  }    
  //QC
  if(plotQC2PtRP && qcCommonHistRes2)
  { 
   (qcCommonHistRes2->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  if(plotQC4PtRP && qcCommonHistRes4)
  { 
   (qcCommonHistRes4->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  if(plotQC6PtRP && qcCommonHistRes6)
  { 
   //(qcCommonHistRes6->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  if(plotQC8PtRP && qcCommonHistRes8)
  { 
   //(qcCommonHistRes8->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  //LYZ2SUM
  if(plotLYZ2SUMPtRP && lyz2sumCommonHistRes)
  { 
   (lyz2sumCommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  //LYZ2PROD
  if(plotLYZ2PRODPtRP && lyz2prodCommonHistRes)
  { 
   (lyz2prodCommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }
  //LYZEP
  if(plotLYZEPPtRP && lyzepCommonHistRes)
  { 
   (lyzepCommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
  }

  //2nd pad is for legend:
  (diffFlowPtAllCanvasRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);  
    
  TLegend* legendDiffFlowPtRP = new TLegend(0.02,0.12,0.97,0.70);
  legendDiffFlowPtRP->SetTextFont(72);
  legendDiffFlowPtRP->SetTextSize(0.06);
 
  //legend's entries:Pt
  TString *entryDiffMCPtRP       = new TString("MC ........ ");
  TString *entryDiffSPPtRP       = new TString("SP ........ ");
  TString *entryDiffGFC2PtRP     = new TString("GFC{2} .... ");
  TString *entryDiffGFC4PtRP     = new TString("GFC{4} .... ");
  TString *entryDiffGFC6PtRP     = new TString("GFC{6} .... ");
  TString *entryDiffGFC8PtRP     = new TString("GFC{8} .... "); 
  TString *entryDiffQC2PtRP      = new TString("QC{2} ..... ");
  TString *entryDiffQC4PtRP      = new TString("QC{4} ..... ");
  TString *entryDiffQC6PtRP      = new TString("QC{6} ..... ");
  TString *entryDiffQC8PtRP      = new TString("QC{8} ..... ");
  TString *entryDiffLYZ2SUMPtRP  = new TString("LYZ{sum} .. ");
  TString *entryDiffLYZ2PRODPtRP = new TString("LYZ{prod} . ");
  TString *entryDiffLYZEPPtRP    = new TString("LYZEP ..... ");
  
  //MC
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetFillStyle(meshStyleDiffFlowPtRP);
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetFillColor(meshColorDiffFlowPtRP);
   entryDiffMCPtRP->Append("M = ");
   (*entryDiffMCPtRP)+=(Long_t)avMultDiffFlowMCRP;
   entryDiffMCPtRP->Append(", N = ");
   (*entryDiffMCPtRP)+=(Long_t)nEvtsDiffFlowMCRP; 
   legendDiffFlowPtRP->AddEntry(mcepCommonHistRes->GetHistDiffFlowPtRP(),entryDiffMCPtRP->Data(),"f");
  }
  
  //SP
  if(plotSPPtRP && spCommonHistRes)
  {
   entryDiffSPPtRP->Append("M = ");
   (*entryDiffSPPtRP)+=(Long_t)avMultDiffFlowSPRP;
   entryDiffSPPtRP->Append(", N = ");
   (*entryDiffSPPtRP)+=(Long_t)nEvtsDiffFlowSPRP; 
   legendDiffFlowPtRP->AddEntry(spCommonHistRes->GetHistDiffFlowPtRP(),entryDiffSPPtRP->Data(),"p");
  }

  //GFC
  if(plotGFC2PtRP && gfcCommonHistRes2)
  {
   entryDiffGFC2PtRP->Append("M = ");
   (*entryDiffGFC2PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC2PtRP->Append(", N = ");
   (*entryDiffGFC2PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowPtRP->AddEntry(gfcCommonHistRes2->GetHistDiffFlowPtRP(),entryDiffGFC2PtRP->Data(),"p");
  }
  if(plotGFC4PtRP && gfcCommonHistRes4)
  {
   entryDiffGFC4PtRP->Append("M = ");
   (*entryDiffGFC4PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC4PtRP->Append(", N = ");
   (*entryDiffGFC4PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowPtRP->AddEntry(gfcCommonHistRes4->GetHistDiffFlowPtRP(),entryDiffGFC4PtRP->Data(),"p");
  }
  if(plotGFC6PtRP && gfcCommonHistRes6)
  {
   entryDiffGFC6PtRP->Append("M = ");
   (*entryDiffGFC6PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC6PtRP->Append(", N = ");
   (*entryDiffGFC6PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowPtRP->AddEntry(gfcCommonHistRes6->GetHistDiffFlowPtRP(),entryDiffGFC6PtRP->Data(),"p");
  } 
  if(plotGFC8PtRP && gfcCommonHistRes8)
  {
   entryDiffGFC8PtRP->Append("M = ");
   (*entryDiffGFC8PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC8PtRP->Append(", N = ");
   (*entryDiffGFC8PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowPtRP->AddEntry(gfcCommonHistRes8->GetHistDiffFlowPtRP(),entryDiffGFC8PtRP->Data(),"p");
  }  
  
  //QC
  if(plotQC2PtRP && qcCommonHistRes2)
  {
   entryDiffQC2PtRP->Append("M = ");
   (*entryDiffQC2PtRP)+=(Long_t)avMultDiffFlowQC2RP;
   entryDiffQC2PtRP->Append(", N = ");
   (*entryDiffQC2PtRP)+=(Long_t)nEvtsDiffFlowQC2RP; 
   legendDiffFlowPtRP->AddEntry(qcCommonHistRes2->GetHistDiffFlowPtRP(),entryDiffQC2PtRP->Data(),"p");
  }
  if(plotQC4PtRP && qcCommonHistRes4)
  {
   entryDiffQC4PtRP->Append("M = ");
   (*entryDiffQC4PtRP)+=(Long_t)avMultDiffFlowQC4RP;
   entryDiffQC4PtRP->Append(", N = ");
   (*entryDiffQC4PtRP)+=(Long_t)nEvtsDiffFlowQC4RP; 
   legendDiffFlowPtRP->AddEntry(qcCommonHistRes4->GetHistDiffFlowPtRP(),entryDiffQC4PtRP->Data(),"p");
  }
  if(plotQC6PtRP && qcCommonHistRes6)
  {
   entryDiffQC6PtRP->Append("M = ");
   (*entryDiffQC6PtRP)+=(Long_t)avMultDiffFlowQC6RP;
   entryDiffQC6PtRP->Append(", N = ");
   (*entryDiffQC6PtRP)+=(Long_t)nEvtsDiffFlowQC6RP; 
   legendDiffFlowPtRP->AddEntry(qcCommonHistRes6->GetHistDiffFlowPtRP(),entryDiffQC6PtRP->Data(),"p");
  }
  if(plotQC8PtRP && qcCommonHistRes8)
  {
   entryDiffQC8PtRP->Append("M = ");
   (*entryDiffQC8PtRP)+=(Long_t)avMultDiffFlowQC8RP;
   entryDiffQC8PtRP->Append(", N = ");
   (*entryDiffQC8PtRP)+=(Long_t)nEvtsDiffFlowQC8RP; 
   legendDiffFlowPtRP->AddEntry(qcCommonHistRes8->GetHistDiffFlowPtRP(),entryDiffQC8PtRP->Data(),"p");
  }
  
  //LYZ2SUM
  if(plotLYZ2SUMPtRP && lyz2sumCommonHistRes)
  {
   entryDiffLYZ2SUMPtRP->Append("M = ");
   (*entryDiffLYZ2SUMPtRP)+=(Long_t)avMultDiffFlowLYZ2SUMRP;
   entryDiffLYZ2SUMPtRP->Append(", N = ");
   (*entryDiffLYZ2SUMPtRP)+=(Long_t)nEvtsDiffFlowLYZ2SUMRP; 
   legendDiffFlowPtRP->AddEntry(lyz2sumCommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZ2SUMPtRP->Data(),"p");
  }
  
  //LYZ2PROD
  if(plotLYZ2PRODPtRP && lyz2prodCommonHistRes)
  {
   entryDiffLYZ2PRODPtRP->Append("M = ");
   (*entryDiffLYZ2PRODPtRP)+=(Long_t)avMultDiffFlowLYZ2PRODRP;
   entryDiffLYZ2PRODPtRP->Append(", N = ");
   (*entryDiffLYZ2PRODPtRP)+=(Long_t)nEvtsDiffFlowLYZ2PRODRP; 
   legendDiffFlowPtRP->AddEntry(lyz2prodCommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZ2PRODPtRP->Data(),"p");
  }
  
  //LYZEP
  if(plotLYZEPPtRP && lyzepCommonHistRes)
  {
   entryDiffLYZEPPtRP->Append("M = ");
   (*entryDiffLYZEPPtRP)+=(Long_t)avMultDiffFlowLYZEPRP;
   entryDiffLYZEPPtRP->Append(", N = ");
   (*entryDiffLYZEPPtRP)+=(Long_t)nEvtsDiffFlowLYZEPRP; 
   legendDiffFlowPtRP->AddEntry(lyzepCommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZEPPtRP->Data(),"p");
  }

  //drawing finally the legend in the 2nd pad:         
  if(textDefault) textDefault->Draw();
  
  if(legendDiffFlowPtRP)
  {
   legendDiffFlowPtRP->SetMargin(0.15);
   legendDiffFlowPtRP->Draw();
  }
 }// end of if(plotDiffFlowPtRP)
 //----------------------------------------------------------------------------------
 
 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Eta, RP):
 if(plotDiffFlowEtaRP)
 {
  TCanvas* diffFlowEtaAllCanvasRP = new TCanvas("Differential Flow (Eta) of RP","Differential Flow (Eta) of RP ",1000,600);
 
  diffFlowEtaAllCanvasRP->Divide(2,1);
 
  //1st pad is for plot:
  (diffFlowEtaAllCanvasRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  
  if(styleHistEta)
  {
   (styleHistEta->GetYaxis())->SetRangeUser(-0.3,1.0);
   styleHistEta->Draw();
  }
  if(pMeshDiffFlowEtaRP)
  {
   pMeshDiffFlowEtaRP->Draw("LFSAME");
  }
 
  //MC 
  if(plotMCEtaRP && mcepCommonHistRes)
  { 
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  //SP 
  if(plotSPEtaRP && spCommonHistRes)
  { 
   (spCommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  //GFC
  if(plotGFC2EtaRP && gfcCommonHistRes2)
  { 
   (gfcCommonHistRes2->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
  } 
  if(plotGFC4EtaRP && gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
  } 
  if(plotGFC6EtaRP && gfcCommonHistRes6)
  { 
   (gfcCommonHistRes6->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
  } 
  if(plotGFC8EtaRP && gfcCommonHistRes8)
  { 
   (gfcCommonHistRes8->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
  }    
  //QC
  if(plotQC2EtaRP && qcCommonHistRes2)
  { 
   (qcCommonHistRes2->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  if(plotQC4EtaRP && qcCommonHistRes4)
  { 
   (qcCommonHistRes4->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  if(plotQC6EtaRP && qcCommonHistRes6)
  { 
   //(qcCommonHistRes6->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  if(plotQC8EtaRP && qcCommonHistRes8)
  { 
   //(qcCommonHistRes8->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  //LYZ2SUM
  if(plotLYZ2SUMEtaRP && lyz2sumCommonHistRes)
  { 
   (lyz2sumCommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  //LYZ2PROD
  if(plotLYZ2PRODEtaRP && lyz2prodCommonHistRes)
  { 
   (lyz2prodCommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
  //LYZEP
  if(plotLYZEPEtaRP && lyzepCommonHistRes)
  { 
   (lyzepCommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
  }
 
  //2nd pad is for legend:
  (diffFlowEtaAllCanvasRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);

  TLegend* legendDiffFlowEtaRP = new TLegend(0.02,0.12,0.97,0.70);
  legendDiffFlowEtaRP->SetTextFont(72);
  legendDiffFlowEtaRP->SetTextSize(0.06);
  
  //legend's entries:
  TString *entryDiffMCEtaRP       = new TString("MC ........ ");
  TString *entryDiffSPEtaRP       = new TString("SP ........ ");
  TString *entryDiffGFC2EtaRP     = new TString("GFC{2} .... ");
  TString *entryDiffGFC4EtaRP     = new TString("GFC{4} .... ");
  TString *entryDiffGFC6EtaRP     = new TString("GFC{6} .... ");
  TString *entryDiffGFC8EtaRP     = new TString("GFC{8} .... "); 
  TString *entryDiffQC2EtaRP      = new TString("QC{2} ..... ");
  TString *entryDiffQC4EtaRP      = new TString("QC{4} ..... ");
  TString *entryDiffQC6EtaRP      = new TString("QC{6} ..... ");
  TString *entryDiffQC8EtaRP      = new TString("QC{8} ..... ");
  TString *entryDiffLYZ2SUMEtaRP  = new TString("LYZ{sum} .. ");
  TString *entryDiffLYZ2PRODEtaRP = new TString("LYZ{prod} . ");
  TString *entryDiffLYZEPEtaRP    = new TString("LYZEP ..... ");
 
  //MC
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->SetFillStyle(meshStyleDiffFlowEtaRP);
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->SetFillColor(meshColorDiffFlowEtaRP);
   entryDiffMCEtaRP->Append("M = ");
   (*entryDiffMCEtaRP)+=(Long_t)avMultDiffFlowMCRP;
   entryDiffMCEtaRP->Append(", N = ");
   (*entryDiffMCEtaRP)+=(Long_t)nEvtsDiffFlowMCRP; 
   legendDiffFlowEtaRP->AddEntry(mcepCommonHistRes->GetHistDiffFlowEtaRP(),entryDiffMCEtaRP->Data(),"f");
  }
  
  //SP
  if(plotSPEtaRP && spCommonHistRes)
  {
   entryDiffSPEtaRP->Append("M = ");
   (*entryDiffSPEtaRP)+=(Long_t)avMultDiffFlowSPRP;
   entryDiffSPEtaRP->Append(", N = ");
   (*entryDiffSPEtaRP)+=(Long_t)nEvtsDiffFlowSPRP; 
   legendDiffFlowEtaRP->AddEntry(spCommonHistRes->GetHistDiffFlowEtaRP(),entryDiffSPEtaRP->Data(),"p");
  }
 
  //GFC
  if(plotGFC2EtaRP && gfcCommonHistRes2)
  {
   entryDiffGFC2EtaRP->Append("M = ");
   (*entryDiffGFC2EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC2EtaRP->Append(", N = ");
   (*entryDiffGFC2EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes2->GetHistDiffFlowEtaRP(),entryDiffGFC2EtaRP->Data(),"p");
  }
  if(plotGFC4EtaRP && gfcCommonHistRes4)
  {
   entryDiffGFC4EtaRP->Append("M = ");
   (*entryDiffGFC4EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC4EtaRP->Append(", N = ");
   (*entryDiffGFC4EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes4->GetHistDiffFlowEtaRP(),entryDiffGFC4EtaRP->Data(),"p");
  }
  if(plotGFC6EtaRP && gfcCommonHistRes6)
  {
   entryDiffGFC6EtaRP->Append("M = ");
   (*entryDiffGFC6EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC6EtaRP->Append(", N = ");
   (*entryDiffGFC6EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes6->GetHistDiffFlowEtaRP(),entryDiffGFC6EtaRP->Data(),"p");
  } 
  if(plotGFC8EtaRP && gfcCommonHistRes8)
  {
   entryDiffGFC8EtaRP->Append("M = ");
   (*entryDiffGFC8EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC8EtaRP->Append(", N = ");
   (*entryDiffGFC8EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes8->GetHistDiffFlowEtaRP(),entryDiffGFC8EtaRP->Data(),"p");
  }  
  
  //QC
  if(plotQC2EtaRP && qcCommonHistRes2)
  {
   entryDiffQC2EtaRP->Append("M = ");
   (*entryDiffQC2EtaRP)+=(Long_t)avMultDiffFlowQC2RP;
   entryDiffQC2EtaRP->Append(", N = ");
   (*entryDiffQC2EtaRP)+=(Long_t)nEvtsDiffFlowQC2RP; 
   legendDiffFlowEtaRP->AddEntry(qcCommonHistRes2->GetHistDiffFlowEtaRP(),entryDiffQC2EtaRP->Data(),"p");
  }
  if(plotQC4EtaRP && qcCommonHistRes4)
  {
   entryDiffQC4EtaRP->Append("M = ");
   (*entryDiffQC4EtaRP)+=(Long_t)avMultDiffFlowQC4RP;
   entryDiffQC4EtaRP->Append(", N = ");
   (*entryDiffQC4EtaRP)+=(Long_t)nEvtsDiffFlowQC4RP; 
   legendDiffFlowEtaRP->AddEntry(qcCommonHistRes4->GetHistDiffFlowEtaRP(),entryDiffQC4EtaRP->Data(),"p");
  }
  if(plotQC6EtaRP && qcCommonHistRes6)
  {
   entryDiffQC6EtaRP->Append("M = ");
   (*entryDiffQC6EtaRP)+=(Long_t)avMultDiffFlowQC6RP;
   entryDiffQC6EtaRP->Append(", N = ");
   (*entryDiffQC6EtaRP)+=(Long_t)nEvtsDiffFlowQC6RP; 
   legendDiffFlowEtaRP->AddEntry(qcCommonHistRes6->GetHistDiffFlowEtaRP(),entryDiffQC6EtaRP->Data(),"p");
  }
  if(plotQC8EtaRP && qcCommonHistRes8)
  {
   entryDiffQC8EtaRP->Append("M = ");
   (*entryDiffQC8EtaRP)+=(Long_t)avMultDiffFlowQC8RP;
   entryDiffQC8EtaRP->Append(", N = ");
   (*entryDiffQC8EtaRP)+=(Long_t)nEvtsDiffFlowQC8RP; 
   legendDiffFlowEtaRP->AddEntry(qcCommonHistRes8->GetHistDiffFlowEtaRP(),entryDiffQC8EtaRP->Data(),"p");
  }
 
  //LYZ2SUM
  if(plotLYZ2SUMEtaRP && lyz2sumCommonHistRes)
  {
   entryDiffLYZ2SUMEtaRP->Append("M = ");
   (*entryDiffLYZ2SUMEtaRP)+=(Long_t)avMultDiffFlowLYZ2SUMRP;
   entryDiffLYZ2SUMEtaRP->Append(", N = ");
   (*entryDiffLYZ2SUMEtaRP)+=(Long_t)nEvtsDiffFlowLYZ2SUMRP; 
   legendDiffFlowEtaRP->AddEntry(lyz2sumCommonHistRes->GetHistDiffFlowEtaRP(),entryDiffLYZ2SUMEtaRP->Data(),"p");
  }
  
  //LYZ2PROD
  if(plotLYZ2PRODEtaRP && lyz2prodCommonHistRes)
  {
   entryDiffLYZ2PRODEtaRP->Append("M = ");
   (*entryDiffLYZ2PRODEtaRP)+=(Long_t)avMultDiffFlowLYZ2PRODRP;
   entryDiffLYZ2PRODEtaRP->Append(", N = ");
   (*entryDiffLYZ2PRODEtaRP)+=(Long_t)nEvtsDiffFlowLYZ2PRODRP; 
   legendDiffFlowEtaRP->AddEntry(lyz2prodCommonHistRes->GetHistDiffFlowEtaRP(),entryDiffLYZ2PRODEtaRP->Data(),"p");
  }
  
  //LYZEP
  if(plotLYZEPEtaRP && lyzepCommonHistRes)
  {
   entryDiffLYZEPEtaRP->Append("M = ");
   (*entryDiffLYZEPEtaRP)+=(Long_t)avMultDiffFlowLYZEPRP;
   entryDiffLYZEPEtaRP->Append(", N = ");
   (*entryDiffLYZEPEtaRP)+=(Long_t)nEvtsDiffFlowLYZEPRP; 
   legendDiffFlowEtaRP->AddEntry(lyzepCommonHistRes->GetHistDiffFlowEtaRP(),entryDiffLYZEPEtaRP->Data(),"p");
  }

  //drawing finally the legend in the 2nd pad: 
  if(textDefault) textDefault->Draw();    
  
  if(legendDiffFlowEtaRP)
  {
   legendDiffFlowEtaRP->SetMargin(0.15);
   legendDiffFlowEtaRP->Draw();
  }
 }// end of if(plotDiffFlowEtaRP)
 //----------------------------------------------------------------------------------
 
 //----------------------------------------------------------------------------------
 // final drawing for plot |(v{method}-v{MC})/v{MC}| as a function of pt for RPs  
 if(plotDiffFlowPtRelativeToMCRP)
 {  
  TCanvas* diffFlowPtRelativeToMCRP = new TCanvas("Differential Flow (Pt) of RPs relative to MC","Differential Flow (Pt) of RPs relative to MC",1000,600);
  
  diffFlowPtRelativeToMCRP->Divide(2,1);
 
  //1st pad is for plot:
  (diffFlowPtRelativeToMCRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
  if(styleHistPt)
  {
   TH1D *styleHistPtReleativeToMC = new TH1D(*styleHistPt);
   (styleHistPtReleativeToMC->GetYaxis())->SetRangeUser(-4.0,4.0);
   (styleHistPtReleativeToMC->GetYaxis())->SetTitle("(v_{n}\{method\} - v_{n}\{MC\})/v_{n}\{MC\}");
   styleHistPtReleativeToMC->Draw();
  }
  
  TH1D *spDiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *gfc2DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *gfc4DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *gfc6DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *gfc8DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *qc2DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *qc4DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *qc6DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *qc8DiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *lyz2sumDiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *lyz2prodDiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
  TH1D *lyzepDiffFlowPtRelativeToMCRP = new TH1D("","",iNbinsPt,dPtMin,dPtMax);
    
  if(mcepCommonHistRes && mcepCommonHistRes->GetHistDiffFlowPtRP())
  {
   for(Int_t p=1;p<iNbinsPt+1;p++)
   {
    Double_t dvnMC = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(p);
    
    // SP:
    Double_t dvnSP = 0.;
    Double_t dvnErrorSP = 0.;
    if(spCommonHistRes && spCommonHistRes->GetHistDiffFlowPtRP())
    {
     dvnSP = (spCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorSP = (spCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnSP!=0.) 
    {
     spDiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnSP-dvnMC)/dvnMC);
     spDiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorSP/dvnMC);
    }
    
    // GFC{2}:
    Double_t dvnGFC2 = 0.;
    Double_t dvnErrorGFC2 = 0.;
    if(gfcCommonHistRes2 && gfcCommonHistRes2->GetHistDiffFlowPtRP())
    {
     dvnGFC2 = (gfcCommonHistRes2->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorGFC2 = (gfcCommonHistRes2->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnGFC2!=0.) 
    {
     gfc2DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnGFC2-dvnMC)/dvnMC);
     gfc2DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorGFC2/dvnMC);
    }
    
    // GFC{4}:
    Double_t dvnGFC4 = 0.;
    Double_t dvnErrorGFC4 = 0.;
    if(gfcCommonHistRes4 && gfcCommonHistRes4->GetHistDiffFlowPtRP())
    {
     dvnGFC4 = (gfcCommonHistRes4->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorGFC4 = (gfcCommonHistRes4->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnGFC4!=0.) 
    {
     gfc4DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnGFC4-dvnMC)/dvnMC);
     gfc4DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorGFC4/dvnMC);
    }
    
    // GFC{6}:
    Double_t dvnGFC6 = 0.;
    Double_t dvnErrorGFC6 = 0.;
    if(gfcCommonHistRes6 && gfcCommonHistRes6->GetHistDiffFlowPtRP())
    {
     dvnGFC6 = (gfcCommonHistRes6->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorGFC6 = (gfcCommonHistRes6->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnGFC6!=0.) 
    {
     gfc6DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnGFC6-dvnMC)/dvnMC);
     gfc6DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorGFC6/dvnMC);
    }
    
    // GFC{8}:
    Double_t dvnGFC8 = 0.;
    Double_t dvnErrorGFC8 = 0.;
    if(gfcCommonHistRes8 && gfcCommonHistRes8->GetHistDiffFlowPtRP())
    {
     dvnGFC8 = (gfcCommonHistRes8->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorGFC8 = (gfcCommonHistRes8->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnGFC8!=0.) 
    {
     gfc8DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnGFC8-dvnMC)/dvnMC);
     gfc8DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorGFC8/dvnMC);
    }
    
    // QC{2}:
    Double_t dvnQC2 = 0.;
    Double_t dvnErrorQC2 = 0.;
    if(qcCommonHistRes2 && qcCommonHistRes2->GetHistDiffFlowPtRP())
    {
     dvnQC2 = (qcCommonHistRes2->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorQC2 = (qcCommonHistRes2->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnQC2!=0.) 
    {
     qc2DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnQC2-dvnMC)/dvnMC);
     qc2DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorQC2/dvnMC);
    }
    
    // QC{4}:
    Double_t dvnQC4 = 0.;
    Double_t dvnErrorQC4 = 0.;
    if(qcCommonHistRes4 && qcCommonHistRes4->GetHistDiffFlowPtRP())
    {
     dvnQC4 = (qcCommonHistRes4->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorQC4 = (qcCommonHistRes4->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnQC4!=0.) 
    {
     qc4DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnQC4-dvnMC)/dvnMC);
     qc4DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorQC4/dvnMC);
    }
    
    // QC{6}:
    Double_t dvnQC6 = 0.;
    Double_t dvnErrorQC6 = 0.;
    if(qcCommonHistRes6 && qcCommonHistRes6->GetHistDiffFlowPtRP())
    {
     dvnQC6 = (qcCommonHistRes6->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorQC6 = (qcCommonHistRes6->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnQC6!=0.) 
    {
     qc6DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnQC6-dvnMC)/dvnMC);
     qc6DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorQC6/dvnMC);
    }
    
    // QC{8}:
    Double_t dvnQC8 = 0.;
    Double_t dvnErrorQC8 = 0.;
    if(qcCommonHistRes8 && qcCommonHistRes8->GetHistDiffFlowPtRP())
    {
     dvnQC8 = (qcCommonHistRes8->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorQC8 = (qcCommonHistRes8->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnQC8!=0.) 
    {
     qc8DiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnQC8-dvnMC)/dvnMC);
     qc8DiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorQC8/dvnMC);
    }
    
    // LYZ2SUM:
    Double_t dvnLYZ2SUM = 0.;
    Double_t dvnErrorLYZ2SUM = 0.;
    if(lyz2sumCommonHistRes && lyz2sumCommonHistRes->GetHistDiffFlowPtRP())
    {
     dvnLYZ2SUM = (lyz2sumCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorLYZ2SUM = (lyz2sumCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnLYZ2SUM!=0.) 
    {
     lyz2sumDiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnLYZ2SUM-dvnMC)/dvnMC);
     lyz2sumDiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorLYZ2SUM/dvnMC);
    }
    
    // LYZ2PROD:
    Double_t dvnLYZ2PROD = 0.;
    Double_t dvnErrorLYZ2PROD = 0.;
    if(lyz2prodCommonHistRes && lyz2prodCommonHistRes->GetHistDiffFlowPtRP())
    {
     dvnLYZ2PROD = (lyz2prodCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorLYZ2PROD = (lyz2prodCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnLYZ2PROD!=0.) 
    {
     lyz2prodDiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnLYZ2PROD-dvnMC)/dvnMC);
     lyz2prodDiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorLYZ2PROD/dvnMC);
    }
      
    // LYZEP:
    Double_t dvnLYZEP = 0.;
    Double_t dvnErrorLYZEP = 0.;
    if(lyzepCommonHistRes && lyzepCommonHistRes->GetHistDiffFlowPtRP())
    {
     dvnLYZEP = (lyzepCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(p);
     dvnErrorLYZEP = (lyzepCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(p);
    }
    if(dvnMC!=0. && dvnLYZEP!=0.) 
    {
     lyzepDiffFlowPtRelativeToMCRP->SetBinContent(p,(dvnLYZEP-dvnMC)/dvnMC);
     lyzepDiffFlowPtRelativeToMCRP->SetBinError(p,dvnErrorLYZEP/dvnMC);
    }  
 
   } // end of for(Int_t p=1;p<iNbinsPt+1;p++)
  } // end of if(mcepCommonHistRes->GetHistDiffFlowPtRP())
  
  
  // final drawings: 
  spDiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorSP);
  spDiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleSP);
  if(plotSPRelativeToMCRP && spCommonHistRes) spDiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  gfc2DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorGFC2);
  gfc2DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleGFC2);
  if(plotGFC2RelativeToMCRP && gfcCommonHistRes2) gfc2DiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  gfc4DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorGFC4);
  gfc4DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleGFC4);
  if(plotGFC4RelativeToMCRP && gfcCommonHistRes4) gfc4DiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  gfc6DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorGFC6);
  gfc6DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleGFC6);
  if(plotGFC6RelativeToMCRP && gfcCommonHistRes6) gfc6DiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  gfc8DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorGFC8);
  gfc8DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleGFC8);
  if(plotGFC8RelativeToMCRP && gfcCommonHistRes8) gfc8DiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  qc2DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorQC2);
  qc2DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleQC2);
  if(plotQC2RelativeToMCRP && qcCommonHistRes2) qc2DiffFlowPtRelativeToMCRP->Draw("PSAME");
  
  qc4DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorQC4);
  qc4DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleQC4);
  if(plotQC4RelativeToMCRP && qcCommonHistRes4) qc4DiffFlowPtRelativeToMCRP->Draw("PSAME");
  
  qc6DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorQC6);
  qc6DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleQC6);
  if(plotQC6RelativeToMCRP && qcCommonHistRes6) qc6DiffFlowPtRelativeToMCRP->Draw("PSAME");
  
  qc8DiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorQC8);
  qc8DiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleQC8);
  if(plotQC8RelativeToMCRP && qcCommonHistRes8) qc8DiffFlowPtRelativeToMCRP->Draw("PSAME");
  
  lyz2sumDiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorLYZ2SUM);
  lyz2sumDiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleLYZ2SUM);
  if(plotLYZ2SUMRelativeToMCRP && lyz2sumCommonHistRes) lyz2sumDiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  lyz2prodDiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorLYZ2SUM);
  lyz2prodDiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleLYZ2SUM);
  if(plotLYZ2SUMRelativeToMCRP && lyz2prodCommonHistRes) lyz2prodDiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  lyzepDiffFlowPtRelativeToMCRP->SetMarkerColor(markerColorLYZEP);
  lyzepDiffFlowPtRelativeToMCRP->SetMarkerStyle(markerStyleLYZEP);
  if(plotLYZEPRelativeToMCRP && lyzepCommonHistRes) lyzepDiffFlowPtRelativeToMCRP->Draw("E1PSAME");
  
  //2nd pad is for legend:
  (diffFlowPtRelativeToMCRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);
  
  TLegend* legendDiffFlowPtRP = new TLegend(0.02,0.12,0.97,0.70);
  legendDiffFlowPtRP->SetTextFont(72);
  legendDiffFlowPtRP->SetTextSize(0.06);
 
  //legend's entries:Pt 
  TString *entryDiffMCPtRP       = new TString("MC ........ ");
  TString *entryDiffSPPtRP       = new TString("SP ........ ");
  TString *entryDiffGFC2PtRP     = new TString("GFC{2} .... ");
  TString *entryDiffGFC4PtRP     = new TString("GFC{4} .... ");
  TString *entryDiffGFC6PtRP     = new TString("GFC{6} .... ");
  TString *entryDiffGFC8PtRP     = new TString("GFC{8} .... "); 
  TString *entryDiffQC2PtRP      = new TString("QC{2} ..... ");
  TString *entryDiffQC4PtRP      = new TString("QC{4} ..... ");
  TString *entryDiffQC6PtRP      = new TString("QC{6} ..... ");
  TString *entryDiffQC8PtRP      = new TString("QC{8} ..... ");
  TString *entryDiffLYZ2SUMPtRP  = new TString("LYZ{sum} .. ");
  TString *entryDiffLYZ2PRODPtRP = new TString("LYZ{prod} . ");
  TString *entryDiffLYZEPPtRP    = new TString("LYZEP ..... ");
  
  //MC
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetFillStyle(meshStyleDiffFlowPtRP);
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetFillColor(meshColorDiffFlowPtRP);
   entryDiffMCPtRP->Append("M = ");
   (*entryDiffMCPtRP)+=(Long_t)avMultDiffFlowMCRP;
   entryDiffMCPtRP->Append(", N = ");
   (*entryDiffMCPtRP)+=(Long_t)nEvtsDiffFlowMCRP; 
   //legendDiffFlowPtRP->AddEntry(mcepCommonHistRes->GetHistDiffFlowPtRP(),entryDiffMCPtRP->Data(),"f");
  }
  
  //SP
  if(plotSPPtRP && spCommonHistRes)
  {
   entryDiffSPPtRP->Append("M = ");
   (*entryDiffSPPtRP)+=(Long_t)avMultDiffFlowSPRP;
   entryDiffSPPtRP->Append(", N = ");
   (*entryDiffSPPtRP)+=(Long_t)nEvtsDiffFlowSPRP; 
   if(plotSPRelativeToMCRP) legendDiffFlowPtRP->AddEntry(spCommonHistRes->GetHistDiffFlowPtRP(),entryDiffSPPtRP->Data(),"p");
  }

  //GFC
  if(plotGFC2PtRP && gfcCommonHistRes2)
  {
   entryDiffGFC2PtRP->Append("M = ");
   (*entryDiffGFC2PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC2PtRP->Append(", N = ");
   (*entryDiffGFC2PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   if(plotGFC2RelativeToMCRP) legendDiffFlowPtRP->AddEntry(gfcCommonHistRes2->GetHistDiffFlowPtRP(),entryDiffGFC2PtRP->Data(),"p");
  }
  if(plotGFC4PtRP && gfcCommonHistRes4)
  {
   entryDiffGFC4PtRP->Append("M = ");
   (*entryDiffGFC4PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC4PtRP->Append(", N = ");
   (*entryDiffGFC4PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   if(plotGFC4RelativeToMCRP) legendDiffFlowPtRP->AddEntry(gfcCommonHistRes4->GetHistDiffFlowPtRP(),entryDiffGFC4PtRP->Data(),"p");
  }
  if(plotGFC6PtRP && gfcCommonHistRes6)
  {
   entryDiffGFC6PtRP->Append("M = ");
   (*entryDiffGFC6PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC6PtRP->Append(", N = ");
   (*entryDiffGFC6PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   if(plotGFC6RelativeToMCRP) legendDiffFlowPtRP->AddEntry(gfcCommonHistRes6->GetHistDiffFlowPtRP(),entryDiffGFC6PtRP->Data(),"p");
  } 
  if(plotGFC8PtRP && gfcCommonHistRes8)
  {
   entryDiffGFC8PtRP->Append("M = ");
   (*entryDiffGFC8PtRP)+=(Long_t)avMultDiffFlowGFCRP;
   entryDiffGFC8PtRP->Append(", N = ");
   (*entryDiffGFC8PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
   if(plotGFC8RelativeToMCRP) legendDiffFlowPtRP->AddEntry(gfcCommonHistRes8->GetHistDiffFlowPtRP(),entryDiffGFC8PtRP->Data(),"p");
  }  
  
  //QC
  if(plotQC2PtRP && qcCommonHistRes2)
  {
   entryDiffQC2PtRP->Append("M = ");
   (*entryDiffQC2PtRP)+=(Long_t)avMultDiffFlowQC2RP;
   entryDiffQC2PtRP->Append(", N = ");
   (*entryDiffQC2PtRP)+=(Long_t)nEvtsDiffFlowQC2RP; 
   if(plotQC2RelativeToMCRP) legendDiffFlowPtRP->AddEntry(qcCommonHistRes2->GetHistDiffFlowPtRP(),entryDiffQC2PtRP->Data(),"p");
  }
  if(plotQC4PtRP && qcCommonHistRes4)
  {
   entryDiffQC4PtRP->Append("M = ");
   (*entryDiffQC4PtRP)+=(Long_t)avMultDiffFlowQC4RP;
   entryDiffQC4PtRP->Append(", N = ");
   (*entryDiffQC4PtRP)+=(Long_t)nEvtsDiffFlowQC4RP; 
   if(plotQC4RelativeToMCRP) legendDiffFlowPtRP->AddEntry(qcCommonHistRes4->GetHistDiffFlowPtRP(),entryDiffQC4PtRP->Data(),"p");
  }
  if(plotQC6PtRP && qcCommonHistRes6)
  {
   entryDiffQC6PtRP->Append("M = ");
   (*entryDiffQC6PtRP)+=(Long_t)avMultDiffFlowQC6RP;
   entryDiffQC6PtRP->Append(", N = ");
   (*entryDiffQC6PtRP)+=(Long_t)nEvtsDiffFlowQC6RP; 
   if(plotQC6RelativeToMCRP) legendDiffFlowPtRP->AddEntry(qcCommonHistRes6->GetHistDiffFlowPtRP(),entryDiffQC6PtRP->Data(),"p");
  }
  if(plotQC8PtRP && qcCommonHistRes8)
  {
   entryDiffQC8PtRP->Append("M = ");
   (*entryDiffQC8PtRP)+=(Long_t)avMultDiffFlowQC8RP;
   entryDiffQC8PtRP->Append(", N = ");
   (*entryDiffQC8PtRP)+=(Long_t)nEvtsDiffFlowQC8RP; 
   if(plotQC8RelativeToMCRP) legendDiffFlowPtRP->AddEntry(qcCommonHistRes8->GetHistDiffFlowPtRP(),entryDiffQC8PtRP->Data(),"p");
  }
  
  //LYZ2SUM
  if(plotLYZ2SUMPtRP && lyz2sumCommonHistRes)
  {
   entryDiffLYZ2SUMPtRP->Append("M = ");
   (*entryDiffLYZ2SUMPtRP)+=(Long_t)avMultDiffFlowLYZ2SUMRP;
   entryDiffLYZ2SUMPtRP->Append(", N = ");
   (*entryDiffLYZ2SUMPtRP)+=(Long_t)nEvtsDiffFlowLYZ2SUMRP; 
   if(plotLYZ2SUMRelativeToMCRP) legendDiffFlowPtRP->AddEntry(lyz2sumCommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZ2SUMPtRP->Data(),"p");
  }
  
  //LYZ2PROD
  if(plotLYZ2PRODPtRP && lyz2prodCommonHistRes)
  {
   entryDiffLYZ2PRODPtRP->Append("M = ");
   (*entryDiffLYZ2PRODPtRP)+=(Long_t)avMultDiffFlowLYZ2PRODRP;
   entryDiffLYZ2PRODPtRP->Append(", N = ");
   (*entryDiffLYZ2PRODPtRP)+=(Long_t)nEvtsDiffFlowLYZ2PRODRP; 
   if(plotLYZ2PRODRelativeToMCRP) legendDiffFlowPtRP->AddEntry(lyz2prodCommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZ2PRODPtRP->Data(),"p");
  }
  
  //LYZEP
  if(plotLYZEPPtRP && lyzepCommonHistRes)
  {
   entryDiffLYZEPPtRP->Append("M = ");
   (*entryDiffLYZEPPtRP)+=(Long_t)avMultDiffFlowLYZEPRP;
   entryDiffLYZEPPtRP->Append(", N = ");
   (*entryDiffLYZEPPtRP)+=(Long_t)nEvtsDiffFlowLYZEPRP; 
   if(plotLYZEPRelativeToMCRP) legendDiffFlowPtRP->AddEntry(lyzepCommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZEPPtRP->Data(),"p");
  }

  //drawing finally the legend in the 2nd pad:         
  if(textDefault) textDefault->Draw();
  
  if(legendDiffFlowPtRP)
  {
   legendDiffFlowPtRP->SetMargin(0.15);
   legendDiffFlowPtRP->Draw();
  } 
 } 
 //----------------------------------------------------------------------------------
 
 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Pt, POI):
 if(plotDiffFlowPtPOI)
 {
  TCanvas* diffFlowPtAllCanvasPOI = new TCanvas("Differential Flow (Pt) of POI","Differential Flow (Pt) of POI ",1000,600);
 
  diffFlowPtAllCanvasPOI->Divide(2,1);
 
  //1st pad is for plot:
  (diffFlowPtAllCanvasPOI->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  
  if(styleHistPt)
  {
   (styleHistPt->GetYaxis())->SetRangeUser(-0.3,1.0);
   styleHistPt->Draw();
  }
  if(pMeshDiffFlowPtPOI)
  {
   pMeshDiffFlowPtPOI->Draw("LFSAME");
  }
 
  //MC 
  if(plotMCPtPOI && mcepCommonHistRes)
  { 
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  //SP 
  if(plotSPPtPOI && spCommonHistRes)
  { 
   (spCommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  //GFC
  if(plotGFC2PtPOI && gfcCommonHistRes2)
  { 
   (gfcCommonHistRes2->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
  } 
  if(plotGFC4PtPOI && gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
  } 
  if(plotGFC6PtPOI && gfcCommonHistRes6)
  { 
   (gfcCommonHistRes6->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
  } 
  if(plotGFC8PtPOI && gfcCommonHistRes8)
  { 
   (gfcCommonHistRes8->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
  }    
  //QC
  if(plotQC2PtPOI && qcCommonHistRes2)
  { 
   (qcCommonHistRes2->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  if(plotQC4PtPOI && qcCommonHistRes4)
  { 
   (qcCommonHistRes4->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  if(plotQC6PtPOI && qcCommonHistRes6)
  { 
   //(qcCommonHistRes6->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  if(plotQC8PtPOI && qcCommonHistRes8)
  { 
   //(qcCommonHistRes8->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  //LYZ2SUM
  if(plotLYZ2SUMPtPOI && lyz2sumCommonHistRes)
  { 
   (lyz2sumCommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  //LYZ2PROD
  if(plotLYZ2PRODPtPOI && lyz2prodCommonHistRes)
  { 
   (lyz2prodCommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
  //LYZEP
  if(plotLYZEPPtPOI && lyzepCommonHistRes)
  { 
   (lyzepCommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
  }
 
  //2nd pad is for legend:
  (diffFlowPtAllCanvasPOI->cd(2))->SetPad(0.75,0.0,1.0,1.0);

  TLegend* legendDiffFlowPtPOI = new TLegend(0.02,0.12,0.97,0.70);
  legendDiffFlowPtPOI->SetTextFont(72);
  legendDiffFlowPtPOI->SetTextSize(0.06);
 
  //legend's entries:
  TString *entryDiffMCPtPOI       = new TString("MC ........ ");
  TString *entryDiffSPPtPOI       = new TString("SP ........ ");
  TString *entryDiffGFC2PtPOI     = new TString("GFC{2} .... ");
  TString *entryDiffGFC4PtPOI     = new TString("GFC{4} .... ");
  TString *entryDiffGFC6PtPOI     = new TString("GFC{6} .... ");
  TString *entryDiffGFC8PtPOI     = new TString("GFC{8} .... "); 
  TString *entryDiffQC2PtPOI      = new TString("QC{2} ..... ");
  TString *entryDiffQC4PtPOI      = new TString("QC{4} ..... ");
  TString *entryDiffQC6PtPOI      = new TString("QC{6} ..... ");
  TString *entryDiffQC8PtPOI      = new TString("QC{8} ..... ");
  TString *entryDiffLYZ2SUMPtPOI  = new TString("LYZ{sum} .. ");
  TString *entryDiffLYZ2PRODPtPOI = new TString("LYZ{prod} . ");
  TString *entryDiffLYZEPPtPOI    = new TString("LYZEP ..... "); 
 
  //MC
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->SetFillStyle(meshStyleDiffFlowPtPOI);
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->SetFillColor(meshColorDiffFlowPtPOI);
   entryDiffMCPtPOI->Append("M = ");
   (*entryDiffMCPtPOI)+=(Long_t)avMultDiffFlowMCPOI;
   entryDiffMCPtPOI->Append(", N = ");
   (*entryDiffMCPtPOI)+=(Long_t)nEvtsDiffFlowMCPOI; 
   legendDiffFlowPtPOI->AddEntry(mcepCommonHistRes->GetHistDiffFlowPtPOI(),entryDiffMCPtPOI->Data(),"f");
  }
  
  //SP
  if(plotSPPtPOI && spCommonHistRes)
  {
   entryDiffSPPtPOI->Append("M = ");
   (*entryDiffSPPtPOI)+=(Long_t)avMultDiffFlowSPPOI;
   entryDiffSPPtPOI->Append(", N = ");
   (*entryDiffSPPtPOI)+=(Long_t)nEvtsDiffFlowSPPOI; 
   legendDiffFlowPtPOI->AddEntry(spCommonHistRes->GetHistDiffFlowPtPOI(),entryDiffSPPtPOI->Data(),"p");
  }

  //GFC
  if(plotGFC2PtPOI && gfcCommonHistRes2)
  {
   entryDiffGFC2PtPOI->Append("M = ");
   (*entryDiffGFC2PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC2PtPOI->Append(", N = ");
   (*entryDiffGFC2PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes2->GetHistDiffFlowPtPOI(),entryDiffGFC2PtPOI->Data(),"p");
  }
  if(plotGFC4PtPOI && gfcCommonHistRes4)
  {
   entryDiffGFC4PtPOI->Append("M = ");
   (*entryDiffGFC4PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC4PtPOI->Append(", N = ");
   (*entryDiffGFC4PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes4->GetHistDiffFlowPtPOI(),entryDiffGFC4PtPOI->Data(),"p");
  }
  if(plotGFC6PtPOI && gfcCommonHistRes6)
  {
   entryDiffGFC6PtPOI->Append("M = ");
   (*entryDiffGFC6PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC6PtPOI->Append(", N = ");
   (*entryDiffGFC6PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes6->GetHistDiffFlowPtPOI(),entryDiffGFC6PtPOI->Data(),"p");
  } 
  if(plotGFC8PtPOI && gfcCommonHistRes8)
  {
   entryDiffGFC8PtPOI->Append("M = ");
   (*entryDiffGFC8PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC8PtPOI->Append(", N = ");
   (*entryDiffGFC8PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes8->GetHistDiffFlowPtPOI(),entryDiffGFC8PtPOI->Data(),"p");
  }  
  
  //QC
  if(plotQC2PtPOI && qcCommonHistRes2)
  {
   entryDiffQC2PtPOI->Append("M = ");
   (*entryDiffQC2PtPOI)+=(Long_t)avMultDiffFlowQC2POI;
   entryDiffQC2PtPOI->Append(", N = ");
   (*entryDiffQC2PtPOI)+=(Long_t)nEvtsDiffFlowQC2POI; 
   legendDiffFlowPtPOI->AddEntry(qcCommonHistRes2->GetHistDiffFlowPtPOI(),entryDiffQC2PtPOI->Data(),"p");
  }
  if(plotQC4PtPOI && qcCommonHistRes4)
  {
   entryDiffQC4PtPOI->Append("M = ");
   (*entryDiffQC4PtPOI)+=(Long_t)avMultDiffFlowQC4POI;
   entryDiffQC4PtPOI->Append(", N = ");
   (*entryDiffQC4PtPOI)+=(Long_t)nEvtsDiffFlowQC4POI; 
   legendDiffFlowPtPOI->AddEntry(qcCommonHistRes4->GetHistDiffFlowPtPOI(),entryDiffQC4PtPOI->Data(),"p");
  }
  if(plotQC6PtPOI && qcCommonHistRes6)
  {
   entryDiffQC6PtPOI->Append("M = ");
   (*entryDiffQC6PtPOI)+=(Long_t)avMultDiffFlowQC6POI;
   entryDiffQC6PtPOI->Append(", N = ");
   (*entryDiffQC6PtPOI)+=(Long_t)nEvtsDiffFlowQC6POI; 
   legendDiffFlowPtPOI->AddEntry(qcCommonHistRes6->GetHistDiffFlowPtPOI(),entryDiffQC6PtPOI->Data(),"p");
  }
  if(plotQC8PtPOI && qcCommonHistRes8)
  {
   entryDiffQC8PtPOI->Append("M = ");
   (*entryDiffQC8PtPOI)+=(Long_t)avMultDiffFlowQC8POI;
   entryDiffQC8PtPOI->Append(", N = ");
   (*entryDiffQC8PtPOI)+=(Long_t)nEvtsDiffFlowQC8POI; 
   legendDiffFlowPtPOI->AddEntry(qcCommonHistRes8->GetHistDiffFlowPtPOI(),entryDiffQC8PtPOI->Data(),"p");
  }
 
  //LYZ2SUM
  if(plotLYZ2SUMPtPOI && lyz2sumCommonHistRes)
  {
   entryDiffLYZ2SUMPtPOI->Append("M = ");
   (*entryDiffLYZ2SUMPtPOI)+=(Long_t)avMultDiffFlowLYZ2SUMPOI;
   entryDiffLYZ2SUMPtPOI->Append(", N = ");
   (*entryDiffLYZ2SUMPtPOI)+=(Long_t)nEvtsDiffFlowLYZ2SUMPOI; 
   legendDiffFlowPtPOI->AddEntry(lyz2sumCommonHistRes->GetHistDiffFlowPtPOI(),entryDiffLYZ2SUMPtPOI->Data(),"p");
  }
  
  //LYZ2PROD
  if(plotLYZ2PRODPtPOI && lyz2prodCommonHistRes)
  {
   entryDiffLYZ2PRODPtPOI->Append("M = ");
   (*entryDiffLYZ2PRODPtPOI)+=(Long_t)avMultDiffFlowLYZ2PRODPOI;
   entryDiffLYZ2PRODPtPOI->Append(", N = ");
   (*entryDiffLYZ2PRODPtPOI)+=(Long_t)nEvtsDiffFlowLYZ2PRODPOI; 
   legendDiffFlowPtPOI->AddEntry(lyz2prodCommonHistRes->GetHistDiffFlowPtPOI(),entryDiffLYZ2PRODPtPOI->Data(),"p");
  }
  
  //LYZEP
  if(plotLYZEPPtPOI && lyzepCommonHistRes)
  {
   entryDiffLYZEPPtPOI->Append("M = ");
   (*entryDiffLYZEPPtPOI)+=(Long_t)avMultDiffFlowLYZEPPOI;
   entryDiffLYZEPPtPOI->Append(", N = ");
   (*entryDiffLYZEPPtPOI)+=(Long_t)nEvtsDiffFlowLYZEPPOI; 
   legendDiffFlowPtPOI->AddEntry(lyzepCommonHistRes->GetHistDiffFlowPtPOI(),entryDiffLYZEPPtPOI->Data(),"p");
  }

  //drawing finally the legend in the 2nd pad: 
  if(textDefault) textDefault->Draw();
          
  if(legendDiffFlowPtPOI)
  {
   legendDiffFlowPtPOI->SetMargin(0.15);
   legendDiffFlowPtPOI->Draw();
  }
 }//end of if(plotDiffFlowPtPOI)
 //----------------------------------------------------------------------------------


 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Eta, POI):
 if(plotDiffFlowEtaPOI)
 {
  TCanvas* diffFlowEtaAllCanvasPOI = new TCanvas("Differential Flow (Eta) of POI","Differential Flow (Eta) of POI ",1000,600);
 
  diffFlowEtaAllCanvasPOI->Divide(2,1);
  
  //1st pad is for plot:
  (diffFlowEtaAllCanvasPOI->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
  if(styleHistEta)
  {
   (styleHistEta->GetYaxis())->SetRangeUser(-0.3,1.0);
   styleHistEta->Draw();
  }
  if(pMeshDiffFlowEtaPOI)
  {
   pMeshDiffFlowEtaPOI->Draw("LFSAME");
  }
 
  //MC 
  if(plotMCEtaPOI && mcepCommonHistRes)
  { 
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  //SP 
  if(plotSPEtaPOI && spCommonHistRes)
  { 
   (spCommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  //GFC
  if(plotGFC2EtaPOI && gfcCommonHistRes2)
  { 
   (gfcCommonHistRes2->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
  } 
  if(plotGFC4EtaPOI && gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
  } 
  if(plotGFC6EtaPOI && gfcCommonHistRes6)
  { 
   (gfcCommonHistRes6->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
  } 
  if(plotGFC8EtaPOI && gfcCommonHistRes8)
  { 
   (gfcCommonHistRes8->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
  }    
  //QC
  if(plotQC2EtaPOI && qcCommonHistRes2)
  { 
   (qcCommonHistRes2->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  if(plotQC4EtaPOI && qcCommonHistRes4)
  { 
   (qcCommonHistRes4->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  if(plotQC6EtaPOI && qcCommonHistRes6)
  { 
   //(qcCommonHistRes6->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  if(plotQC8EtaPOI && qcCommonHistRes8)
  { 
   //(qcCommonHistRes8->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  //LYZ2SUM
  if(plotLYZ2SUMEtaPOI && lyz2sumCommonHistRes)
  { 
   (lyz2sumCommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  //LYZ2PROD
  if(plotLYZ2PRODEtaPOI && lyz2prodCommonHistRes)
  { 
   (lyz2prodCommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
  //LYZEP
  if(plotLYZEPEtaPOI && lyzepCommonHistRes)
  { 
   (lyzepCommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
  }
 
  //2nd pad is for legend:
  (diffFlowEtaAllCanvasPOI->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
  TLegend* legendDiffFlowEtaPOI = new TLegend(0.02,0.12,0.97,0.70);
  legendDiffFlowEtaPOI->SetTextFont(72);
  legendDiffFlowEtaPOI->SetTextSize(0.06);
 
  //legend's entries:
  TString *entryDiffMCEtaPOI       = new TString("MC ........ ");
  TString *entryDiffSPEtaPOI       = new TString("SP ........ ");
  TString *entryDiffGFC2EtaPOI     = new TString("GFC{2} .... ");
  TString *entryDiffGFC4EtaPOI     = new TString("GFC{4} .... ");
  TString *entryDiffGFC6EtaPOI     = new TString("GFC{6} .... ");
  TString *entryDiffGFC8EtaPOI     = new TString("GFC{8} .... "); 
  TString *entryDiffQC2EtaPOI      = new TString("QC{2} ..... ");
  TString *entryDiffQC4EtaPOI      = new TString("QC{4} ..... ");
  TString *entryDiffQC6EtaPOI      = new TString("QC{6} ..... ");
  TString *entryDiffQC8EtaPOI      = new TString("QC{8} ..... ");
  TString *entryDiffLYZ2SUMEtaPOI  = new TString("LYZ{sum} .. ");
  TString *entryDiffLYZ2PRODEtaPOI = new TString("LYZ{prod} . ");
  TString *entryDiffLYZEPEtaPOI    = new TString("LYZEP ..... ");
 
  //MC
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->SetFillStyle(meshStyleDiffFlowEtaPOI);
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->SetFillColor(meshColorDiffFlowEtaPOI);
   entryDiffMCEtaPOI->Append("M = ");
   (*entryDiffMCEtaPOI)+=(Long_t)avMultDiffFlowMCPOI;
   entryDiffMCEtaPOI->Append(", N = ");
   (*entryDiffMCEtaPOI)+=(Long_t)nEvtsDiffFlowMCPOI; 
   legendDiffFlowEtaPOI->AddEntry(mcepCommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffMCEtaPOI->Data(),"f");
  }
  
  //SP
  if(plotSPEtaPOI && spCommonHistRes)
  {
   entryDiffSPEtaPOI->Append("M = ");
   (*entryDiffSPEtaPOI)+=(Long_t)avMultDiffFlowSPPOI;
   entryDiffSPEtaPOI->Append(", N = ");
   (*entryDiffSPEtaPOI)+=(Long_t)nEvtsDiffFlowSPPOI; 
   legendDiffFlowEtaPOI->AddEntry(spCommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffSPEtaPOI->Data(),"p");
  }

  //GFC
  if(plotGFC2EtaPOI && gfcCommonHistRes2)
  {
   entryDiffGFC2EtaPOI->Append("M = ");
   (*entryDiffGFC2EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC2EtaPOI->Append(", N = ");
   (*entryDiffGFC2EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes2->GetHistDiffFlowEtaPOI(),entryDiffGFC2EtaPOI->Data(),"p");
  }
  if(plotGFC4EtaPOI && gfcCommonHistRes4)
  {
   entryDiffGFC4EtaPOI->Append("M = ");
   (*entryDiffGFC4EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC4EtaPOI->Append(", N = ");
   (*entryDiffGFC4EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes4->GetHistDiffFlowEtaPOI(),entryDiffGFC4EtaPOI->Data(),"p");
  }
  if(plotGFC6EtaPOI && gfcCommonHistRes6)
  {
   entryDiffGFC6EtaPOI->Append("M = ");
   (*entryDiffGFC6EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC6EtaPOI->Append(", N = ");
   (*entryDiffGFC6EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes6->GetHistDiffFlowEtaPOI(),entryDiffGFC6EtaPOI->Data(),"p");
  } 
  if(plotGFC8EtaPOI && gfcCommonHistRes8)
  {
   entryDiffGFC8EtaPOI->Append("M = ");
   (*entryDiffGFC8EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
   entryDiffGFC8EtaPOI->Append(", N = ");
   (*entryDiffGFC8EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
   legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes8->GetHistDiffFlowEtaPOI(),entryDiffGFC8EtaPOI->Data(),"p");
  }  
 
  //QC
  if(plotQC2EtaPOI && qcCommonHistRes2)
  {
   entryDiffQC2EtaPOI->Append("M = ");
   (*entryDiffQC2EtaPOI)+=(Long_t)avMultDiffFlowQC2POI;
   entryDiffQC2EtaPOI->Append(", N = ");
   (*entryDiffQC2EtaPOI)+=(Long_t)nEvtsDiffFlowQC2POI; 
   legendDiffFlowEtaPOI->AddEntry(qcCommonHistRes2->GetHistDiffFlowEtaPOI(),entryDiffQC2EtaPOI->Data(),"p");
  }
  if(plotQC4EtaPOI && qcCommonHistRes4)
  {
   entryDiffQC4EtaPOI->Append("M = ");
   (*entryDiffQC4EtaPOI)+=(Long_t)avMultDiffFlowQC4POI;
   entryDiffQC4EtaPOI->Append(", N = ");
   (*entryDiffQC4EtaPOI)+=(Long_t)nEvtsDiffFlowQC4POI; 
   legendDiffFlowEtaPOI->AddEntry(qcCommonHistRes4->GetHistDiffFlowEtaPOI(),entryDiffQC4EtaPOI->Data(),"p");
  }
  if(plotQC6EtaPOI && qcCommonHistRes6)
  {
   entryDiffQC6EtaPOI->Append("M = ");
   (*entryDiffQC6EtaPOI)+=(Long_t)avMultDiffFlowQC6POI;
   entryDiffQC6EtaPOI->Append(", N = ");
   (*entryDiffQC6EtaPOI)+=(Long_t)nEvtsDiffFlowQC6POI; 
   legendDiffFlowEtaPOI->AddEntry(qcCommonHistRes6->GetHistDiffFlowEtaPOI(),entryDiffQC6EtaPOI->Data(),"p");
  }
  if(plotQC8EtaPOI && qcCommonHistRes8)
  {
   entryDiffQC8EtaPOI->Append("M = ");
   (*entryDiffQC8EtaPOI)+=(Long_t)avMultDiffFlowQC8POI;
   entryDiffQC8EtaPOI->Append(", N = ");
   (*entryDiffQC8EtaPOI)+=(Long_t)nEvtsDiffFlowQC8POI; 
   legendDiffFlowEtaPOI->AddEntry(qcCommonHistRes8->GetHistDiffFlowEtaPOI(),entryDiffQC8EtaPOI->Data(),"p");
  }
 
  //LYZ2SUM
  if(plotLYZ2SUMEtaPOI && lyz2sumCommonHistRes)
  {
   entryDiffLYZ2SUMEtaPOI->Append("M = ");
   (*entryDiffLYZ2SUMEtaPOI)+=(Long_t)avMultDiffFlowLYZ2SUMPOI;
   entryDiffLYZ2SUMEtaPOI->Append(", N = ");
   (*entryDiffLYZ2SUMEtaPOI)+=(Long_t)nEvtsDiffFlowLYZ2SUMPOI; 
   legendDiffFlowEtaPOI->AddEntry(lyz2sumCommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffLYZ2SUMEtaPOI->Data(),"p");
  }
  
  //LYZ2PROD
  if(plotLYZ2PRODEtaPOI && lyz2prodCommonHistRes)
  {
   entryDiffLYZ2PRODEtaPOI->Append("M = ");
   (*entryDiffLYZ2PRODEtaPOI)+=(Long_t)avMultDiffFlowLYZ2PRODPOI;
   entryDiffLYZ2PRODEtaPOI->Append(", N = ");
   (*entryDiffLYZ2PRODEtaPOI)+=(Long_t)nEvtsDiffFlowLYZ2PRODPOI; 
   legendDiffFlowEtaPOI->AddEntry(lyz2prodCommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffLYZ2PRODEtaPOI->Data(),"p");
  }
  
  //LYZEP
  if(plotLYZEPEtaPOI && lyzepCommonHistRes)
  {
   entryDiffLYZEPEtaPOI->Append("M = ");
   (*entryDiffLYZEPEtaPOI)+=(Long_t)avMultDiffFlowLYZEPPOI;
   entryDiffLYZEPEtaPOI->Append(", N = ");
   (*entryDiffLYZEPEtaPOI)+=(Long_t)nEvtsDiffFlowLYZEPPOI; 
   legendDiffFlowEtaPOI->AddEntry(lyzepCommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffLYZEPEtaPOI->Data(),"p");
  }

  //drawing finally the legend in the 2nd pad:   
  if(textDefault) textDefault->Draw();
      
  if(legendDiffFlowEtaPOI)
  {
   legendDiffFlowEtaPOI->SetMargin(0.15);
   legendDiffFlowEtaPOI->Draw();
  }
 }//end of if(plotDiffFlowEtaPOI)
 //----------------------------------------------------------------------------------


 //=====================================================================================
 
 
}

void LoadPlotLibraries(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------

  //==================================================================================  
  //load needed libraries:
  gSystem->AddIncludePath("-I$ROOTSYS/include");
  //gSystem->Load("libTree");

  // for AliRoot
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG2flowCommon");
  cerr<<"libPWG2flowCommon loaded ..."<<endl;
  
  }
  
  else if (mode==mLocalSource) {
 
    // In root inline compile
  
    // Constants  
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCumuConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEvent.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
    
    // Cuts
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    
    // Output histosgrams
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
       
    cout << "finished loading macros!" << endl;  
    
  }  
  
}


