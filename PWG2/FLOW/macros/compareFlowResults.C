// add comment
//
// RF = Reference Flow
// RP = Reference Particles
// POI = Particles Of Interest

// Set here which plots will be shown by default:
// Results:
Bool_t plotReferenceFlow = kTRUE; // reference flow 
Bool_t plotIntFlowPOI = kTRUE; // integrated flow of POIs
Bool_t plotDiffFlowPtPOI = kTRUE; // differential flow v(pt) for POIs
Bool_t plotDiffFlowEtaPOI = kTRUE; // differential flow v(eta) for POIs
Bool_t plotIntFlowRP = kTRUE; // integrated flow of RPs
Bool_t plotDiffFlowPtRP = kTRUE; // differential flow v(pt) for RPs
Bool_t plotDiffFlowEtaRP = kTRUE; // differential flow v(eta) for RPs
// Results relative to MC |(v{MC}-v{method})/v{MC}|:
Bool_t plotReferenceFlowRelativeToMC = kTRUE; // plot |(v{MC}-v{method})/v{MC}| for reference flow
Bool_t plotIntFlowRelativeToMCPOI = kTRUE; // plot |(v{MC}-v{method})/v{MC}| for integrated flow of POIs   
Bool_t plotDiffFlowPtRelativeToMCPOI = kTRUE; // plot |(v{MC}-v{method})/v{MC}| as a function of pt for POIs 
Bool_t plotDiffFlowEtaRelativeToMCPOI = kTRUE; // plot |(v{MC}-v{method})/v{MC}| as a function of eta for POIs
Bool_t plotIntFlowRelativeToMCRP = kTRUE; // plot |(v{MC}-v{method})/v{MC}| for integrated flow of RPs
Bool_t plotDiffFlowPtRelativeToMCRP = kTRUE; // plot |(v{MC}-v{method})/v{MC}| as a function of pt for RPs   
Bool_t plotDiffFlowEtaRelativeToMCRP = kTRUE; // plot |(v{MC}-v{method})/v{MC}| as a function of eta for RPs   
// Set here if the legends will be shown on the plots:
Bool_t showLegend = kTRUE; 
Bool_t showLegendDiffFlow = kTRUE;
// Some quick settings:
Bool_t showOnlyReferenceFlow = kFALSE;
Bool_t showResultsRelativeToMC = kTRUE;
Bool_t showOnlyPlotsForPOIs = kFALSE;
Bool_t showOnlyPlotsForRPs = kFALSE;

const Int_t nMethods = 12;
TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
TList *list[nMethods] = {NULL}; // lists holding histograms for each flow analysis method

enum libModes{mLocal,mLocalSource};

//void newCompare(TString analysisType="",Int_t analysisMode=mLocalSource)
void compareFlowResults(TString analysisType="",Int_t analysisMode=mLocal)
{
 // 1. analysisType: "ESD", "AOD", "ESDMC0", "ESDMC1"; for Monte Carlo and 'on the fly' use simply "";
 // 2. analysisMode: if analysisMode = mLocal -> analyze data on your computer using aliroot
 //                  if analysisMode = mLocalSource -> analyze data on your computer using root + source files 
  
 // Load needed libraries:
 LoadLibrariesCFR(analysisMode); 

 // Access common output file:
 TString outputFileName = "AnalysisResults.root"; // name of the common output file
 TFile *outputFile = AccessOutputFile(outputFileName);
 
 // Access from common output file the TDirectoryFile's for each flow analysis method
 // and from them the lists holding histograms with final results:
 GetListsWithHistograms(outputFile,analysisType);
 
 // Global settings which will affect all plots:
 GlobalSettings();
  
 // Calling the functions to produce the final plots:
 if(plotReferenceFlow) PlotReferenceFlow();
 if(!showOnlyReferenceFlow)
 { 
  if(!showOnlyPlotsForRPs)
  {
   if(plotIntFlowPOI) PlotIntFlowPOI(); 
   if(plotDiffFlowPtPOI) PlotDiffFlowPtPOI();
   if(plotDiffFlowEtaPOI) PlotDiffFlowEtaPOI();
  }
  if(!showOnlyPlotsForPOIs)
  {
   if(plotIntFlowRP) PlotIntFlowRP(); 
   if(plotDiffFlowPtRP) PlotDiffFlowPtRP();
   if(plotDiffFlowEtaRP) PlotDiffFlowEtaRP();  
  }
  if(showResultsRelativeToMC)
  {
   if(plotReferenceFlowRelativeToMC) PlotReferenceFlowRelativeToMC();
   if(!showOnlyPlotsForRPs)
   {
    if(plotIntFlowRelativeToMCPOI) PlotIntFlowRelativeToMCPOI();
    if(plotDiffFlowPtRelativeToMCPOI) PlotDiffFlowPtRelativeToMCPOI();
    if(plotDiffFlowEtaRelativeToMCPOI) PlotDiffFlowEtaRelativeToMCPOI();
   }
   if(!showOnlyPlotsForPOIs)
   {
    if(plotIntFlowRelativeToMCRP) PlotIntFlowRelativeToMCRP();
    if(plotDiffFlowPtRelativeToMCRP) PlotDiffFlowPtRelativeToMCRP();
    if(plotDiffFlowEtaRelativeToMCRP) PlotDiffFlowEtaRelativeToMCRP();
   }
  }
 }

} // end of void newCompare(TString analysisType="",Int_t analysisMode=mLocal)

// ===========================================================================================

void PlotReferenceFlow()
{
 // Make a plot which compares the results for reference flow.
  
 // Settings for methods:
 const Int_t nMethods = 13;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","6,QC","8,GFC","8,QC","FQD","LYZ1SUM","LYZ1PROD"};
 Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21,21,21,21,21}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kGray;
   
 Plot(nMethods,method,methodMarkerStyle,methodMarkerColor,
      methodUsedToMakeErrorMesh,meshStyle,meshColor,"RF");
 
} // end of void PlotReferenceFlow()

// ===========================================================================================

void PlotIntFlowPOI()
{
 // Make a plot which compares the results for reference flow.
  
 // Settings for methods:
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21,21}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kRed-10;
   
 Plot(nMethods,method,methodMarkerStyle,methodMarkerColor,
      methodUsedToMakeErrorMesh,meshStyle,meshColor,"POI");
  
} // end of void PlotIntFlowPOI()

// ===========================================================================================

void PlotIntFlowRP()
{
 // Make a plot which compares the results for reference flow.
  
 // Settings for methods:
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21,21}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kBlue-10;
   
 Plot(nMethods,method,methodMarkerStyle,methodMarkerColor,
      methodUsedToMakeErrorMesh,meshStyle,meshColor,"RP");
  
} // end of void PlotIntFlowRP()

// ===========================================================================================

void PlotDiffFlowPtPOI()
{
 // Make a plot which compares the results for differential flow of POIs vs pt.

 // Settings for methods:
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {20,3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kRed,kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kRed-10;
 
 PlotDiffFlow(nMethods,method,methodMarkerStyle,methodMarkerColor,
              methodUsedToMakeErrorMesh,meshStyle,meshColor,"Pt","POI");
 
} // end of void PlotDiffFlowPtPOI()

// ===========================================================================================

void PlotDiffFlowEtaPOI()
{
 // Make a plot which compares the results for differential flow of POIs vs eta.

 // Settings for methods:
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {20,3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kRed,kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kRed-10;
 
 PlotDiffFlow(nMethods,method,methodMarkerStyle,methodMarkerColor,
              methodUsedToMakeErrorMesh,meshStyle,meshColor,"Eta","POI");
 
} // end of void PlotDiffFlowEtaPOI()

// ===========================================================================================

void PlotDiffFlowPtRP()
{
 // Make a plot which compares the results for differential flow of RPs vs pt.

 // Settings for methods:
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {20,3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kRed,kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kBlue-10;
 
 PlotDiffFlow(nMethods,method,methodMarkerStyle,methodMarkerColor,
              methodUsedToMakeErrorMesh,meshStyle,meshColor,"Pt","RP");
 
} // end of void PlotDiffFlowPtRP()

// ===========================================================================================

void PlotDiffFlowEtaRP()
{
 // Make a plot which compares the results for differential flow of RPs vs eta.

 // Settings for methods:
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {20,3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kRed,kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 // Settings for error mesh:
 TString methodUsedToMakeErrorMesh = "MCEP";
 Int_t meshStyle = 1001;
 Int_t meshColor = kBlue-10;
 
 PlotDiffFlow(nMethods,method,methodMarkerStyle,methodMarkerColor,
              methodUsedToMakeErrorMesh,meshStyle,meshColor,"Eta","RP");
 
} // end of void PlotDiffFlowEtaRP()

// ===========================================================================================

void PlotReferenceFlowRelativeToMC()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for reference flow.
  
 // Settings for methods:
 const Int_t nMethods = 12;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","6,QC","8,GFC","8,QC","FQD","LYZ1SUM","LYZ1PROD"};
 Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21,21,21,21}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack};
   
 PlotRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"RF");

} // end of void PlotReferenceFlowRelativeToMC()

// ===========================================================================================

void PlotIntFlowRelativeToMCPOI()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for integrated flow of POIs.
  
 // Settings for methods:
 const Int_t nMethods = 9;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3,kRed-3};
   
 PlotRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"POI");

} // end of void PlotIntFlowRelativeToMCPOI()

// ===========================================================================================

void PlotIntFlowRelativeToMCRP()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for integrated flow of RPs.
  
 // Settings for methods:
 const Int_t nMethods = 9;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3,kBlue-3};
   
 PlotRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"RP");

} // end of void PlotIntFlowRelativeToMCRP()

// ===========================================================================================

void PlotDiffFlowPtRelativeToMCPOI()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for differential flow of POIs vs pt.

 // Settings for methods:
 const Int_t nMethods = 9;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 
 PlotDiffFlowRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"Pt","POI");

} // end of void PlotDiffFlowPtRelativeToMCPOI()

// ===========================================================================================
  
void PlotDiffFlowEtaRelativeToMCPOI()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for differential flow of POIs vs eta.

 // Settings for methods:
 const Int_t nMethods = 9;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 
 PlotDiffFlowRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"Eta","POI");

} // end of void PlotDiffFlowEtaRelativeToMCPOI()

// ===========================================================================================
  
void PlotDiffFlowPtRelativeToMCRP()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for differential flow of RPs vs pt.

 // Settings for methods:
 const Int_t nMethods = 9;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 
 PlotDiffFlowRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"Pt","RP");

} // end of void PlotDiffFlowPtRelativeToMCRP()

// ===========================================================================================
  
void PlotDiffFlowEtaRelativeToMCRP()
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for differential flow of RPs vs eta.

 // Settings for methods:
 const Int_t nMethods = 9;
 TString method[nMethods] = {"SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","8,GFC","LYZ2SUM","LYZ2PROD"};
 Int_t methodMarkerStyle[nMethods] = {3,21,21,20,20,25,24,22,22}; // see available marker styles in TAttMarker
 Int_t methodMarkerColor[nMethods] = {kViolet-6,kAzure-7,kOrange-7,kAzure+3,kOrange+3,kAzure-7,kAzure+3,kYellow+3,kGreen+3};
 
 PlotDiffFlowRelativeToMC(nMethods,method,methodMarkerStyle,methodMarkerColor,"Eta","RP");

} // end of void PlotDiffFlowPtRelativeToMCRP()

// ===========================================================================================
        
TGraph* GetErrorMesh(Int_t nPts, Double_t result, Double_t error, Int_t meshStyle, Int_t meshColor)
{
 // Make an error mesh from the specified method.
 
 TGraph *g = new TGraph(nPts);
 g->SetPoint(1,0,result+error);
 g->SetPoint(2,nPts+1,result+error);
 g->SetPoint(3,nPts+1,result-error);
 g->SetPoint(4,0,result-error);
 g->SetPoint(5,0,result-error);    
 g->SetFillStyle(meshStyle);
 g->SetFillColor(meshColor);

 return g;

} // end of GetErrorMesh(Int_t nPoints, Double_t result,Double_t error)

// ===========================================================================================

TGraphErrors* GetGraphErrors(Double_t x, Double_t result, Double_t error, Int_t markerStyle, Int_t markerColor)
{
 // From the result and error for each method make TGraphErrors.
 
 TGraphErrors *ge = NULL;
 
 ge = new TGraphErrors(1);
 ge->SetPoint(0,x,result);
 ge->SetPointError(0,0,error);
 ge->SetMarkerStyle(markerStyle);
 ge->SetMarkerColor(markerColor);

 return ge;

} // end of TGraphErrors* GetGraphErrors(Double_t x, Double_t result, Double_t error, Int_t markerStyle, Int_t markerColor)

// ===========================================================================================
 
void PlotDiffFlow(Int_t nMethods, TString *method, Int_t *methodMarkerStyle, Int_t *methodMarkerColor,
                  TString methodUsedToMakeErrorMesh, Int_t meshStyle, Int_t meshColor, TString ptEta, TString rpPoi)
{
 // Make plot for differential flow.

 TCanvas *c = NULL;
 Int_t sizeX = 1000; // canvas size in pixels along x
 Int_t sizeY = 600; // canvas size in pixels along y
 TString title = Form("Differential Flow vs %s (%s)",ptEta.Data(),rpPoi.Data()); 
 if(!showLegendDiffFlow) sizeX = 0.75*sizeX;
 c = new TCanvas(title.Data(),title.Data(),sizeX,sizeY);
 if(showLegendDiffFlow)
 {
  c->Divide(2,1);
  c->cd(1)->SetPad(0.0,0.0,0.75,1.0);
 } 
 // Style histogram:
 StyleHistDiffFlow(ptEta.Data(),rpPoi.Data())->Draw();
 // Error mesh:  
 TGraph *errorMesh = GetErrorMeshDiffFlow(methodUsedToMakeErrorMesh.Data(),rpPoi.Data(),ptEta.Data());
 if(errorMesh) 
 {
  errorMesh->SetFillStyle(meshStyle);
  errorMesh->SetFillColor(meshColor);
  errorMesh->Draw("lfsame");
 }
 // Results of methods:
 for(Int_t b=0;b<nMethods;b++)
 {
  if(method[b]==methodUsedToMakeErrorMesh) continue;
  TH1D *hist = GetResultHistogram(method[b].Data(),rpPoi.Data(),ptEta.Data());
  if(hist)
  {
   hist->SetMarkerStyle(methodMarkerStyle[b]);
   hist->SetMarkerColor(methodMarkerColor[b]);
   hist->Draw("e1psame");
  }
 } 
 if(showLegendDiffFlow)
 {
  c->cd(2)->SetPad(0.73,0.0,0.97,1.0);
  DefaultTextInLegend()->Draw();
  LegendDiffFlow(nMethods,method,methodMarkerStyle,methodMarkerColor,
                 methodUsedToMakeErrorMesh,meshStyle,meshColor,ptEta,rpPoi)->Draw();
 }
} // end of void PlotDiffFlow(...)

// ===========================================================================================

void PlotDiffFlowRelativeToMC(Int_t nMethods, TString *method, Int_t *methodMarkerStyle,
                              Int_t *methodMarkerColor, TString ptEta, TString rpPoi)
{
 // Make plot for differential flow.
 
 TString title = Form("Differential Flow vs %s (%s) relative to MCEP",ptEta.Data(),rpPoi.Data()); 
 // MCEP:
 TH1D *mcep = GetResultHistogram("MCEP",rpPoi.Data(),ptEta.Data()); // MCEP result and error:
 if(!mcep)
 {
  cout<<"WARNING: MCEP histogram not available in making the plot for "<<title.Data()<<" !!!!"<<endl;   
  return;    
 } 
 TCanvas *c = NULL;
 Int_t sizeX = 1000; // canvas size in pixels along x
 Int_t sizeY = 600; // canvas size in pixels along y
 if(!showLegendDiffFlow) sizeX = 0.75*sizeX;
 c = new TCanvas(title.Data(),title.Data(),sizeX,sizeY);
 if(showLegendDiffFlow)
 {
  c->Divide(2,1);
  c->cd(1)->SetPad(0.0,0.0,0.75,1.0);
 } 
 // Style histogram:
 Int_t n = 2; // to be improved (accessed from common control histograms)
 TH1D *styleHist = StyleHistDiffFlow(ptEta.Data(),rpPoi.Data());
 styleHist->GetYaxis()->SetTitle(Form("(v_{%d}\{MCEP\}-v_{%d}\{method\})/v_{%d}\{MCEP\}",n,n,n));
 styleHist->SetTitle(Form("Differential Flow #font[72]{vs} %s (%s) relative to MCEP",ptEta.Data(),rpPoi.Data())); 
 styleHist->SetMinimum(-10.); // to be improved
 styleHist->SetMaximum(10.); // to be improved
 styleHist->SetTitle(title.Data());
 styleHist->Draw();
 // Methods:
 for(Int_t nm=0;nm<nMethods;nm++)
 {
  TH1D *hist = NULL;
  if(GetResultHistogram(method[nm].Data(),rpPoi.Data(),ptEta.Data()))
  {
   hist = (TH1D*)(GetResultHistogram(method[nm].Data(),rpPoi.Data(),ptEta.Data())->Clone());
  }
  if(hist)
  {
   Int_t nBins = hist->GetNbinsX(); 
   for(Int_t b=1;b<=nBins;b++)
   { 
    Double_t mcepResult = mcep->GetBinContent(b);
    Double_t mcepError = mcep->GetBinError(b);
    if((TMath::Abs(mcepResult) > 1.e-44) && (TMath::Abs(mcepError) > 1.e-44)) 
    {
     Double_t result = hist->GetBinContent(b);
     Double_t error = hist->GetBinError(b); 
     //if((TMath::Abs(result)>1.e-44) && TMath::Abs(error)>1.e-44)) 
     if(TMath::Abs(result)>1.e-44)
     {
      error = pow(pow(error/mcepResult,2.)+pow(result*mcepError/pow(mcepResult,2.),2.),0.5); // Do not switch with the next line!
      result = (mcepResult-result)/mcepResult;
      hist->SetBinContent(b,result);
      hist->SetBinError(b,error);
     } else // end of if((TMath::Abs(result)>1.e-44) && TMath::Abs(error)>1.e-44))
       {
        hist->SetBinContent(b,0.);
        hist->SetBinError(b,0.);          
       }
    } else // end of if(TMath::Abs(mcepResult) > 1.e-44 && TMath::Abs(mcepError) > 1.e-44)  
      {
       hist->SetBinContent(b,0.);
       hist->SetBinError(b,0.);     
      }
   } // end of for(b=1;b<=nBins;b++) 
   hist->SetMarkerStyle(methodMarkerStyle[nm]);
   hist->SetMarkerColor(methodMarkerColor[nm]);
   hist->Draw("e1psame");
  } // end of if(hist)
 } // end of for(Int_t nm=0;nm<nMethods;nm++)  
 if(showLegendDiffFlow)
 {
  c->cd(2)->SetPad(0.73,0.0,0.97,1.0);
  DefaultTextInLegend()->Draw();
  LegendDiffFlow(nMethods,method,methodMarkerStyle,methodMarkerColor,"",-1,-1,ptEta,rpPoi)->Draw();
 }
} // end of void PlotDiffFlowRelativeToMC(...)

// ===========================================================================================

TGraph* GetErrorMeshDiffFlow(TString methodUsedToMakeErrorMesh, TString rpPoi, TString ptEta)
{
 // Error mesh for differential flow.
 
 TH1D *hist = GetResultHistogram(methodUsedToMakeErrorMesh.Data(),rpPoi.Data(),ptEta.Data());
 
 Double_t dMin = 0.;
 if(ptEta == "Pt")
 { 
  dMin = AliFlowCommonConstants::GetMaster()->GetPtMin();
 } else if(ptEta == "Eta")
   {
    dMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();
   }  
 TGraph* errorMesh = NULL;
 if(hist)
 {
  Int_t nBins = hist->GetNbinsX();
  Double_t binWidth = hist->GetBinWidth(1); // assuming that all bins have the same width
  // Counting the non-empty bins: 
  Int_t nNonEmptyBins=0;
  for(Int_t i=1;i<nBins+1;i++)
  {
   if(!(hist)->GetBinError(i)==0.0))
   {
    nNonEmptyBins++;
   }
  }   
  errorMesh = new TGraph(2*nNonEmptyBins+1); 
  Double_t value=0.,error=0.;
  Int_t count=1;
  Double_t xFirst=0.,yUpFirst=0.; // needed to close up the mesh
  for(Int_t i=1;i<nBins+1;i++)
  {
   // Setting up the upper limit of the mesh:
   value = hist->GetBinContent(i);
   error = hist->GetBinError(i);   
   if(!(error==0.0))
   {    
    errorMesh->SetPoint(count++,(i-0.5)*binWidth+dMin,value+error);
    if(xFirst==0.)
    {
     xFirst=(i-0.5)*binWidth+dMin;
     yUpFirst=value+error;
    }
   } 
  }   
  for(Int_t i=nBins+1;i<2*nBins+1;i++)
  {
   // Setting up the lower limit of the mesh:
   value = hist->GetBinContent(2*nBins+1-i);
   error = hist->GetBinError(2*nBins+1-i); 
   if(!(error==0.0))
   {      
    errorMesh->SetPoint(count++,(2*nBins-i+0.5)*binWidth+dMin,value-error);
   }  
  }
  // Closing the mesh area:
  errorMesh->SetPoint(2*nNonEmptyBins+1,xFirst,yUpFirst);   
 } // end if(hist)
 
 return errorMesh;
 
} // end of TGraph* GetErrorMeshDiffFlow(TString methodUsedToMakeErrorMesh, TString ptEta, TString rpPoi)

// ===========================================================================================

TH1D* StyleHistDiffFlow(TString ptEta, TString rpPoi)
{
 // Style histogram for differential flow.

 Int_t n = 2; // to be improved - access from common control histogram
 TH1D *styleHistDiffFlow = NULL;
 if(ptEta == "Pt")
 {
  Int_t iNbinsPt  = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Double_t dPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();
  Double_t dPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
  styleHistDiffFlow = new TH1D("","styleHistDiffFlow",iNbinsPt,dPtMin,dPtMax);
  styleHistDiffFlow->SetTitle(Form("Differential Flow #font[72]{vs} p_{t} (%s)",rpPoi.Data()));
  styleHistDiffFlow->SetXTitle("p_{t} [GeV]");
  styleHistDiffFlow->SetYTitle(Form("v_{%d}",n));
 } 
 else if(ptEta == "Eta")
 {
  Int_t iNbinsEta  = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  Double_t dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();
  Double_t dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();
  styleHistDiffFlow = new TH1D("","",iNbinsEta,dEtaMin,dEtaMax);
  styleHistDiffFlow->SetTitle(Form("Differential Flow #font[72]{vs} #eta (%s)",rpPoi.Data()));
  styleHistDiffFlow->SetXTitle("#eta");
  styleHistDiffFlow->SetYTitle(Form("v_{%d}",n));
 }
 if(styleHistDiffFlow)
 {
  styleHistDiffFlow->SetMinimum(-0.25); // to be improved - implement algorithm for this 
  styleHistDiffFlow->SetMaximum(1.); // to be improved - implement algorithm for this 
  //styleHistDiffFlow->GetYaxis()->SetLabelSize(0.05);
  //styleHistDiffFlow->GetYaxis()->SetTitleSize(0.06);
  //styleHistDiffFlow->GetYaxis()->SetTitleOffset(0.55);
  //styleHistDiffFlow->GetXaxis()->SetLabelSize(0.05);
  //styleHistDiffFlow->GetXaxis()->SetTitleSize(0.06);
  //styleHistDiffFlow->GetXaxis()->SetTitleOffset(0.6);
  //styleHistDiffFlow->GetXaxis()->SetLabelOffset(0.02);
 }
 
 return styleHistDiffFlow;

} // end of TH1D* StyleHistDiffFlow(TString ptEta, TString rpPoi)

// ===========================================================================================

void Plot(const Int_t nMethods,TString *method,Int_t *methodMarkerStyle,Int_t *methodMarkerColor,
         TString methodUsedToMakeErrorMesh,Int_t meshStyle,Int_t meshColor,TString rfRpPoi)
{
 // Make a plot for reference and integrated flow.
 
 TString title = "";
 if(rfRpPoi == "RF")
 {
  title = "Reference Flow";
 } else if(rfRpPoi == "POI") 
   {
    title = "Integrated Flow (POI)";
   } else if(rfRpPoi == "RP") 
     {
      title = "Integrated Flow (RP)";
     }

 Double_t x = 0.; // determines position of the marker on x axis
 Double_t results[nMethods] = {0.};
 Double_t errors[nMethods] = {0.};
 TGraphErrors *ge[nMethods] = {NULL};
 TGraph *errorMesh = NULL;
 for(Int_t b=0;b<nMethods;b++)
 {
  x = b+0.5; 
  TH1D *hist = NULL;
  hist = GetResultHistogram(method[b].Data(),rfRpPoi.Data()); 
  if(hist)
  {
   results[b] = hist->GetBinContent(1);
   errors[b] = hist->GetBinError(1);
   if(TMath::Abs(results[b])>1.e-44) 
   {
    ge[b] = GetGraphErrors(x,results[b],errors[b],methodMarkerStyle[b],methodMarkerColor[b]);
   }
   if(strcmp(method[b].Data(),methodUsedToMakeErrorMesh.Data()) == 0)
   {
    errorMesh = GetErrorMesh(nMethods+1,results[b],errors[b],meshStyle,meshColor);
   }
  } else
    {
     //cout<<"WARNING: For a method "<<method[b].Data()<<" couldn't get the histogram with result"<<endl;
     //cout<<"         for "<<title.Data()<<" !!!! "<<endl;
    } 
 } // end of for(Int_t b=0;b<nMethods;b++)

 // Final drawing:
 TCanvas *c = NULL;
 // Settings for canvas:
 Int_t sizeX = 1000; // canvas size in pixels along x
 Int_t sizeY = 600; // canvas size in pixels along y
 if(!showLegend) sizeX = 0.75*sizeX;
 c = new TCanvas(title.Data(),title.Data(),sizeX,sizeY);
 if(showLegend)
 {
  c->Divide(2,1);
  c->cd(1)->SetPad(0.0,0.0,0.75,1.0);
 } 
 StyleHist(title,nMethods,method,results,errors)->Draw();
 if(errorMesh)errorMesh->Draw("lfsame");
 for(Int_t b=0;b<nMethods;b++)
 {
  if(ge[b])ge[b]->Draw("psame");
 } 
 if(showLegend)
 {
  c->cd(2)->SetPad(0.73,0.0,0.97,1.0);
  DefaultTextInLegend()->Draw();
  Legend(nMethods,method,rfRpPoi)->Draw();
 }
 
 return;

} // end of Plot(...)

// ===========================================================================================
    
void PlotRelativeToMC(const Int_t nMethods, TString *method, Int_t *methodMarkerStyle,
                      Int_t *methodMarkerColor, TString rfRpPoi)
{
 // Make a plot |(v{MC}-v{method})/v{MC}| for reference and integrated flow.

 TString title = "";
 if(rfRpPoi == "RF")
 {
  title = "Reference Flow relative to MCEP";
 } else if(rfRpPoi == "POI") 
   {
    title = "Integrated Flow (POI) relative to MCEP";
   } else if(rfRpPoi == "RP") 
     {
      title = "Integrated Flow (RP) relative to MCEP";
     }

 Double_t x = 0.; // determines position of the marker on x axis
 Double_t results[nMethods] = {0.};
 Double_t errors[nMethods] = {0.};
 // MCEP result and error:
 TH1D *mcep = GetResultHistogram("MCEP","RF");
 Double_t mcepResult = 0.;
 Double_t mcepError = 0.;
 if(mcep)
 {
  mcepResult = mcep->GetBinContent(1);
  mcepError = mcep->GetBinError(1);
 } else
   {
    cout<<"WARNING: MCEP histogram not available in making the plot for "<<title.Data()<<" !!!!"<<endl;   
    return;    
   }
 if(TMath::Abs(mcepResult) < 1.e-44 || TMath::Abs(mcepError) < 1.e-44) 
 {
  cout<<"WARNING: Result or error for v{MCEP} is zero in making the plot for "<<title.Data()<<" !!!!"<<endl;   
  return;
 } 
 TGraphErrors *ge[nMethods] = {NULL};
 for(Int_t b=0;b<nMethods;b++)
 {
  x = b+0.5; 
  TH1D *hist = NULL;
  hist = GetResultHistogram(method[b].Data(),rfRpPoi.Data()); 
  if(hist)
  {
   results[b] = hist->GetBinContent(1);
   errors[b] = hist->GetBinError(1);
   if(TMath::Abs(results[b])>1.e-44) 
   {
    errors[b] = pow(pow(errors[b]/mcepResult,2.)+pow(results[b]*mcepError/pow(mcepResult,2.),2.),0.5); // Do not switch with the next line!
    results[b] = (mcepResult-results[b])/mcepResult;
    ge[b] = GetGraphErrors(x,results[b],errors[b],methodMarkerStyle[b],methodMarkerColor[b]);
   }
  } else
    {
     //cout<<"WARNING: For a method "<<method[b].Data()<<" couldn't get the histogram with result"<<endl;
     //cout<<"         for "<<title.Data()<<" !!!! "<<endl;
    } 
 } // end of for(Int_t b=0;b<nMethods;b++)

 // Final drawing:
 TCanvas *c = NULL;
 // Settings for canvas:
 Int_t sizeX = 1000; // canvas size in pixels along x
 Int_t sizeY = 600; // canvas size in pixels along y
 if(!showLegend) sizeX = 0.75*sizeX;
 c = new TCanvas(title.Data(),title.Data(),sizeX,sizeY);
 if(showLegend)
 {
  c->Divide(2,1);
  c->cd(1)->SetPad(0.0,0.0,0.75,1.0);
 } 
 // Style histogram:
 Int_t n = 2; // to be improved (accessed from common control histograms)
 TH1D *styleHist = StyleHist(title,nMethods,method,results,errors);
 styleHist->GetYaxis()->SetTitleOffset(1.25);
 styleHist->GetYaxis()->SetTitleSize(0.03);
 styleHist->GetYaxis()->SetLabelSize(0.03);
 styleHist->GetYaxis()->SetTitle(Form("(v_{%d}\{MCEP\}-v_{%d}\{method\})/v_{%d}\{MCEP\}",n,n,n));
 styleHist->Draw();
 // Methods:
 for(Int_t b=0;b<nMethods;b++)
 {
  if(ge[b])ge[b]->Draw("psame");
 } 
 if(showLegend)
 {
  c->cd(2)->SetPad(0.73,0.0,0.97,1.0);
  DefaultTextInLegend()->Draw();
  Legend(nMethods,method,rfRpPoi)->Draw();
 }
 
 return;

} // end of void PlotRelativeToMC(...) 

// ===========================================================================================

TPaveText* DefaultTextInLegend()
{
 // Determine the default text in legend.
 
 TPaveText *textDefault = new TPaveText(0.05,0.77,0.95,0.90,"NDC");
 textDefault->SetTextFont(72);
 textDefault->SetTextSize(0.08);
 textDefault->AddText("Average Multiplicity");
 textDefault->AddText("and");
 textDefault->AddText("Number of Events");
 textDefault->SetFillStyle(0); // white instead of default grey

 return textDefault;
 
} // end of TPaveText* DefaultTextInLegend()
 
// =========================================================================================== 

TPaveText* Legend(Int_t nMethods, TString *method, TString rfRpPoi)
{
 // Make a legend.
 
 TPaveText *legend = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
 legend->SetTextFont(72);
 legend->SetTextSize(0.06);
 legend->SetFillStyle(0); // white instead of default grey
 const Int_t nLegendEntries = 11;
 TString legendEntries[nLegendEntries] = {"MCEP ...... ","SP ........ ","GFC ....... ",
                                          "QC{2} ..... ","QC{4} ..... ","QC{6} ..... ",
                                          "QC{8} ..... ","FQD ....... ","LYZ{sum} .. ",
                                          "LYZ{prod} . ","LYZEP ..... "}; 
 TString temp = "";
 Int_t gfcCounter = 0; // represent "2,GFC", "4,GFC", "6,GFC" and "8,GFC" with the same entry GFC
 for(Int_t b=0;b<nMethods;b++)
 {
  for(Int_t le=0;le<nLegendEntries;le++)
  {
   if(legendEntries[le].Contains(method[b].Data())) // this is either MCEP, FQD, SP or LYZEP
   { 
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rfRpPoi);
    legend->AddText(temp.Data());
   }
   else if(method[b].Contains("GFC") && gfcCounter == 0 && legendEntries[le].Contains("GFC")) // GFC
   {
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rfRpPoi);
    legend->AddText(temp.Data()); 
    gfcCounter++;  
   } 
   else if(method[b].Contains("QC") && legendEntries[le].Contains("QC")) // QC
   {
    for(Int_t o=1;o<=4;o++) // QC order
    {
     if(method[b].Contains(Form("%d",2*o)) && legendEntries[le].Contains(Form("%d",2*o)))
     {
      temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rfRpPoi);
      legend->AddText(temp.Data());
     }
    } // end of for(Int o=1;o<=4;o++) // QC order
   }
   else if((method[b].Contains("LYZ1SUM")||method[b].Contains("LYZ2SUM")) && legendEntries[le].Contains("LYZ{sum}"))
   {
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rfRpPoi);
    legend->AddText(temp.Data());
   }
   else if((method[b].Contains("LYZ1PROD")||method[b].Contains("LYZ2PROD")) && legendEntries[le].Contains("LYZ{prod}"))
   {
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rfRpPoi);
    legend->AddText(temp.Data());
   }
  } // end of for(Int_t le=0;le<nLegendEntries;le++)
 } // end of for(Int_t b=0;b<nMethods;b++)
  
 return legend;

} // end of TPaveText* Legend(Int_t nMethods, TString *method, TString rfRpPoi)
// =========================================================================================== 

TLegend* LegendDiffFlow(Int_t nMethods, TString *method, Int_t *methodMarkerStyle, Int_t *methodMarkerColor,
                        TString methodUsedToMakeErrorMesh, Int_t meshStyle, Int_t meshColor, TString ptEta, TString rpPoi)
{
 // Make a legend for differential flow.
 
 TLegend *legend = new TLegend(0.0,0.12,0.99,0.70);
 legend->SetMargin(0.15);
 legend->SetTextFont(72);
 legend->SetTextSize(0.06);
 legend->SetFillStyle(0); // white instead of default grey
 const Int_t nLegendEntries = 14;
 TString legendEntries[nLegendEntries] = {"MCEP ...... ","SP ........ ","GFC{2} .... ","QC{2} ..... ",
                                          "GFC{4} .... ","QC{4} ..... ","GFC{6} .... ","QC{6} ..... ",
                                          "GFC{8} .... ","QC{8} ..... ","FQD ....... ","LYZ{sum} .. ",
                                          "LYZ{prod} . ","LYZEP ..... "};                                                                                     
 
 TH1D *hist = NULL;
 TString temp = "";
 for(Int_t b=0;b<nMethods;b++)
 {
  for(Int_t le=0;le<nLegendEntries;le++)
  {
   if(legendEntries[le].Contains(method[b].Data())) // this is either MCEP, FQD, SP or LYZEP
   { 
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rpPoi);
    hist = GetResultHistogram(method[b].Data(),rpPoi.Data(),ptEta.Data());
    if(hist)
    {
     hist->SetMarkerStyle(methodMarkerStyle[b]);
     hist->SetMarkerColor(methodMarkerColor[b]);  
     if(methodUsedToMakeErrorMesh == method[b].Data())
     {
      hist->SetFillStyle(meshStyle);
      hist->SetFillColor(meshColor);
      legend->AddEntry(hist,temp.Data(),"f");    
     } else
       {
        legend->AddEntry(hist,temp.Data(),"p");    
       }
    } // end of if(hist)   
   }       
   else if(method[b].Contains("GFC") && legendEntries[le].Contains("GFC")) // GFC
   {
    for(Int_t o=1;o<=4;o++) // GFC order
    {
     if(method[b].Contains(Form("%d",2*o)) && legendEntries[le].Contains(Form("%d",2*o)))
     { 
      temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rpPoi);
      hist = GetResultHistogram(method[b].Data(),rpPoi.Data(),ptEta.Data());
      if(hist)
      {
       hist->SetMarkerStyle(methodMarkerStyle[b]);
       hist->SetMarkerColor(methodMarkerColor[b]); 
       if(methodUsedToMakeErrorMesh == method[b].Data())
       {
        hist->SetFillStyle(meshStyle);
        hist->SetFillColor(meshColor);
        legend->AddEntry(hist,temp.Data(),"f");    
       } else
         {
          legend->AddEntry(hist,temp.Data(),"p");    
         }
      } // end of if(hist)           
     }
    }
   } 
   else if(method[b].Contains("QC") && legendEntries[le].Contains("QC")) // QC
   {
    for(Int_t o=1;o<=4;o++) // QC order
    {
     if(method[b].Contains(Form("%d",2*o)) && legendEntries[le].Contains(Form("%d",2*o)))
     {
      temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rpPoi);
      hist = GetResultHistogram(method[b].Data(),rpPoi.Data(),ptEta.Data());
      if(hist)
      {
       hist->SetMarkerStyle(methodMarkerStyle[b]);
       hist->SetMarkerColor(methodMarkerColor[b]);   
       if(methodUsedToMakeErrorMesh == method[b].Data())
       {
        hist->SetFillStyle(meshStyle);
        hist->SetFillColor(meshColor);
        legend->AddEntry(hist,temp.Data(),"f");    
       } else
         {
          legend->AddEntry(hist,temp.Data(),"p");    
         }
      } // end of if(hist)             
     }
    } // end of for(Int o=1;o<=4;o++) // QC order
   }
   else if((method[b].Contains("LYZ1SUM")||method[b].Contains("LYZ2SUM")) && legendEntries[le].Contains("LYZ{sum}"))
   {
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rpPoi);
    hist = GetResultHistogram(method[b].Data(),rpPoi.Data(),ptEta.Data());
    if(hist)
    {
     hist->SetMarkerStyle(methodMarkerStyle[b]);
     hist->SetMarkerColor(methodMarkerColor[b]);  
     if(methodUsedToMakeErrorMesh == method[b].Data())
     {
      hist->SetFillStyle(meshStyle);
      hist->SetFillColor(meshColor);
      legend->AddEntry(hist,temp.Data(),"f");    
     } else
       {
        legend->AddEntry(hist,temp.Data(),"p");    
       }         
    } // end of if(hist)
   }
   else if((method[b].Contains("LYZ1PROD")||method[b].Contains("LYZ2PROD")) && legendEntries[le].Contains("LYZ{prod}"))
   {
    temp = legendEntries[le]+GetAvMultiplicityAndNoOfEvents(method[b].Data(),rpPoi);
    hist = GetResultHistogram(method[b].Data(),rpPoi.Data(),ptEta.Data());
    if(hist)
    {
     hist->SetMarkerStyle(methodMarkerStyle[b]);
     hist->SetMarkerColor(methodMarkerColor[b]);  
     if(methodUsedToMakeErrorMesh == method[b].Data())
     {
      hist->SetFillStyle(meshStyle);
      hist->SetFillColor(meshColor);
      legend->AddEntry(hist,temp.Data(),"f");    
     } else
       {
        legend->AddEntry(hist,temp.Data(),"p");    
       }      
    } // end of if(hist)      
   }
  } // end of for(Int_t le=0;le<nLegendEntries;le++)
 } // end of for(Int_t b=0;b<nMethods;b++)
   
 return legend;

} // end of TLegend* Legend(...)

// =========================================================================================== 

TString GetAvMultiplicityAndNoOfEvents(TString method, TString rfRpPoi)
{
 // Get average multiplicity and number of events for specified method and return it as "M = <AvM>, N = <N>".
 
 TString MN = ""; // string to hold "M = <AvM>, N = <N>"
 Long_t N = 0; // number of events 
 Double_t M = 0.; // average multiplicity
 
 TH1F *hist = GetControlHistogram(method,rfRpPoi); 
 if(hist)
 {
  N = hist->GetEntries();
  M = hist->GetMean();
  MN.Append("M = ");
  MN+=(Long_t)M;
  MN.Append(", N = ");
  MN+=N;
 } else 
   {
    MN.Append("n/a");
   }
 
 return MN;
 
} // end of TString GetAvMultiplicityAndNoOfEvents(TString *method, TString rfRpPoi)

// ===========================================================================================

TH1F* GetControlHistogram(TString method, TString rfRpPoi)
{
 // Get the control histogram for specified method holding the average multiplicity and number of events.
 
 AliFlowCommonHist *commonHist = NULL; 
 TString methodName = method.Data(); 
 Int_t cumulantOrder = 0;
 if(method.Contains("GFC"))
 {
  TString methodNameTemp1 = method;   
  TString methodNameTemp2 = method;   
  methodName = methodNameTemp1.Remove(0,2);
  cumulantOrder = methodNameTemp2.Remove(1,4).Atoi();      
 } else if(method.Contains("QC"))
   {
    TString methodNameTemp1 = method;   
    TString methodNameTemp2 = method;   
    methodName = methodNameTemp1.Remove(0,2);
    cumulantOrder = methodNameTemp2.Remove(1,3).Atoi();      
   }

 // Get for specified methodName (and cumulantOrder, if needed) the common control histogram:
 // (to be improved - this can certainly be implemented better, but some redesign of the flow code is first in order)
 for(Int_t l=0;l<nMethods;l++)
 {
  TString temp = "";
  if(list[l]) 
  {
   temp = TString(list[l]->GetName());
  }
  if(temp.Contains(methodName.Data()))
  {
   // Access from the common list the needed common result histogram:
   if(!(methodName.Contains("QC")))
   {
    commonHist = dynamic_cast<AliFlowCommonHist*> list[l]->FindObject(Form("AliFlowCommonHist%s",methodName.Data()));
   } 
   else if(methodName=="QC")
   {
    if(cumulantOrder==2)
    {
     commonHist = dynamic_cast<AliFlowCommonHist*> list[l]->FindObject("AliFlowCommonHist2ndOrderQC");    
    } 
    else if(cumulantOrder==4)
    {
     commonHist = dynamic_cast<AliFlowCommonHist*> list[l]->FindObject("AliFlowCommonHist4thOrderQC");    
    }
    else if(cumulantOrder==6)
    {
     commonHist = dynamic_cast<AliFlowCommonHist*> list[l]->FindObject("AliFlowCommonHist6thOrderQC");    
    }
    else if(cumulantOrder==8)
    {
     commonHist = dynamic_cast<AliFlowCommonHist*> list[l]->FindObject("AliFlowCommonHist8thOrderQC");    
    }
    else
    {
     cout<<"WARNING: You have specified cumulant order to be "<<cumulantOrder<<" !!!!"<<endl;
     cout<<"         There are no results for this cumulant order. "<<endl;
    }
   } // end of else if(methodName=="QC")
  } // end of if(temp.Contains(methodName.Data()))
 } // end of for(Int_t l=0;l<nMethods;l++)
 
 if(!commonHist)
 {
  //cout<<"WARNING: For a method "<<method.Data()<<" couldn't access the common hist result !!!!"<<endl;
  //cout<<"         File absent? Or perhaps a typo in the method's name?"<<endl;
  if(methodName.Contains("QC"))
  {
   //cout<<"         Or impossible cumulant order? Otherwise.... :'( "<<endl;
  }
 }

 // Finally, from commonHist access the required histogram: 
 TH1F *hist = NULL;
 if(commonHist) 
 { 
  if(rfRpPoi == "RF" || rfRpPoi == "RP")
  {
   hist = commonHist->GetHistMultRP();
  } else if(rfRpPoi == "POI")
    {
     hist = commonHist->GetHistMultPOI(); 
    } else 
      {
       cout<<"WARNING: Impossible TString rfRpPoi in GetControlHistogram() !!!!"<<endl;
       exit(0);
      } 
 }

 return hist;

} // end of TH1D* GetControlHistogram(TString method, TString rfRpPoi)

// =========================================================================================== 

TH1D* StyleHist(TString title, Int_t nMethods, TString *method, Double_t *results, Double_t *errors)
{
 // Make style histogram.
 TH1D *styleHist = new TH1D("",title.Data(),nMethods,0,nMethods);
 Int_t n=2; // to be improved (access this from common control histogram)
 Double_t styleHistMin = 44.;
 Double_t styleHistMax = -44.;
 for(Int_t b=0;b<nMethods;b++)
 {
  // Form bin labels from method's names:
  styleHist->GetXaxis()->SetBinLabel(b+1,Form("v_{%d}\{%s\}",n,method[b].Data()));
  if(method[b].Contains("LYZ1SUM") || method[b].Contains("LYZ2SUM"))
  {
   styleHist->GetXaxis()->SetBinLabel(b+1,Form("v_{%d}\{LYZ,sum\}",n)); 
  }
  if(method[b].Contains("LYZ1PROD") || method[b].Contains("LYZ2PROD"))
  {
   styleHist->GetXaxis()->SetBinLabel(b+1,Form("v_{%d}\{LYZ,prod\}",n)); 
  }
  // Establish min and max values for style histograms:
  if(TMath::Abs(results[b])>1.e-44)
  {
   if(styleHistMin > results[b]-errors[b]) styleHistMin = results[b]-errors[b]; // min value
   if(styleHistMax < results[b]+errors[b]) styleHistMax = results[b]+errors[b]; // max value    
  }
 }
 
 styleHist->SetMinimum(0.99*styleHistMin);
 styleHist->SetMaximum(1.01*styleHistMax);
 /*
 styleHist->GetYaxis()->SetNdivisions(5,5,0);
 styleHist->GetYaxis()->SetLabelSize(0.05);
 styleHist->GetXaxis()->SetLabelSize(0.07);
 styleHist->GetXaxis()->SetLabelOffset(0.01);
 */
 
 return styleHist;

} // end of TH1D* StyleHist(TString *title, Int_t nMethods, TString *method, Double_t *results, Double_t *errors)

// ===========================================================================================

TFile* AccessOutputFile(TString outputFileName)
{
 // Access the common output file.
    
 TFile *outputFile = NULL;
 if(!(gSystem->AccessPathName(Form("%s%s%s",gSystem->pwd(),"/",outputFileName.Data()),kFileExists)))
 {
  outputFile = TFile::Open(outputFileName.Data(),"READ");
 } else
   {
    cout<<endl;
    cout<<"WARNING: Couldn't find the file "<<outputFileName.Data()<<" in "<<endl;
    cout<<"         directory "<<gSystem->pwd()<<" !!!!"<<endl;
    cout<<endl;
    exit(0);
   }
  
 return outputFile;
 
} // end of TFile* AccessOutputFile(TString outputFileName)
  
// ===========================================================================================
 
void GetListsWithHistograms(TFile *outputFile, TString analysisType)
{
 // Access from common output file the TDirectoryFile's for each flow analysis method
 // and from them the lists holding histograms with final results:
 
 TString fileName[nMethods]; 
 TDirectoryFile *dirFile[nMethods] = {NULL}; 
 TString listName[nMethods]; 
 Int_t fileCounter = 0;
 for(Int_t i=0;i<nMethods;i++)
 {
  // Form a file name for each method:
  fileName[i]+="output";
  fileName[i]+=method[i].Data();
  fileName[i]+="analysis";
  fileName[i]+=analysisType.Data();
  // Access this file:
  dirFile[i] = (TDirectoryFile*)outputFile->FindObjectAny(fileName[i].Data());
  // Form a list name for each method:
  listName[i]+="cobj";
  listName[i]+=method[i].Data();
  // Access this list:
  if(dirFile[i])
  {
   if(!(dirFile[i]->GetNkeys() == 0))
   {
    dirFile[i]->GetObject(listName[i].Data(),list[i]);
   } else
     {
      TString temp = listName[i].Data();
      cout<<"WARNING: Couldn't access a list holding histograms for "<<temp.Remove(0,4).Data()<<" method !!!!"<<endl;     
     } 
   fileCounter++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileName[i].Data()<<".root !!!!"<<endl;
    }
 } // end of for(Int_t i=0;i<nMethods;i++) 
 
 // If no files were found most probably the 'TString analysisType' was specified wrongly:
 if(fileCounter == 0)
 {
  cout<<endl; 
  cout<<"Did you specify 'TString analysisType' correctly? Can be \"\", \"ESD\", \"AOD\", etc."<<endl;
  cout<<endl; 
  exit(0);
 }
    
} // end of void GetListsWithHistograms(TFile *outputFile, TString analysisType)
 
// ===========================================================================================

TH1D* GetResultHistogram(TString method, TString rfRpPoi, TString ptEta="")
{
 // Get the specified histogram holding result.
 
 AliFlowCommonHistResults *commonHistRes = NULL; 
 TString methodName = method.Data(); 
 Int_t cumulantOrder = 0;
 if(method.Contains("GFC"))
 {
  TString methodNameTemp1 = method;   
  TString methodNameTemp2 = method;   
  methodName = methodNameTemp1.Remove(0,2);
  cumulantOrder = methodNameTemp2.Remove(1,4).Atoi();      
 } else if(method.Contains("QC"))
   {
    TString methodNameTemp1 = method;   
    TString methodNameTemp2 = method;   
    methodName = methodNameTemp1.Remove(0,2);
    cumulantOrder = methodNameTemp2.Remove(1,3).Atoi();      
   }
 // Get for specified methodName (and cumulantOrder, if needed) the common hist result:
 // (to be improved - this can certainly be implemented better, but some redesign of the flow code is first in order)
 for(Int_t l=0;l<nMethods;l++)
 {
  TString temp = "";
  if(list[l]) 
  {
   temp = TString(list[l]->GetName());
  }
  if(temp.Contains(methodName.Data()))
  {
   // Access from the common list the needed common result histogram:
   if(!(methodName.Contains("GFC") || methodName.Contains("QC")))
   {
    commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject(Form("AliFlowCommonHistResults%s",methodName.Data()));
   } 
   else if(methodName=="GFC")
   {
    if(cumulantOrder==2)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults2ndOrderGFC");    
    } 
    else if(cumulantOrder==4)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults4thOrderGFC");    
    }
    else if(cumulantOrder==6)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults6thOrderGFC");    
    }
    else if(cumulantOrder==8)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults8thOrderGFC");    
    } 
    else
    {
     cout<<"WARNING: You have specified cumulant order to be "<<cumulantOrder<<" !!!!"<<endl;
     cout<<"         That is really not funny.... "<<endl;
    }
   } // end of else if(methodName=="GFC")
   else if(methodName=="QC")
   {
    if(cumulantOrder==2)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults2ndOrderQC");    
    } 
    else if(cumulantOrder==4)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults4thOrderQC");    
    }
    else if(cumulantOrder==6)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults6thOrderQC");    
    }
    else if(cumulantOrder==8)
    {
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> list[l]->FindObject("AliFlowCommonHistResults8thOrderQC");    
    }
    else
    {
     cout<<"WARNING: You have specified cumulant order to be "<<cumulantOrder<<" !!!!"<<endl;
     cout<<"         There are no results for this cumulant order. "<<endl;
    }
   } // end of else if(methodName=="QC")
  } // end of if(temp.Contains(methodName.Data()))
 } // end of for(Int_t l=0;l<nMethods;l++)
 
 if(!commonHistRes)
 {
  //cout<<"WARNING: For a method "<<method.Data()<<" couldn't access the common hist result !!!!"<<endl;
  //cout<<"         File absent? Or perhaps a typo in the method's name?"<<endl;
  if(methodName.Contains("GFC") || methodName.Contains("QC"))
  {
   //cout<<"         Or impossible cumulant order? Otherwise.... :'( "<<endl;
  }
 }
 

 // Finally, from commonHistRes access the required histogram: 
 TH1D *hist = NULL;
 if(commonHistRes) 
 { 
  if(rfRpPoi == "RF")
  {
   hist = commonHistRes->GetHistIntFlow();
  } else if(rfRpPoi == "RP")
    {
     if(ptEta == "")
     {
      hist = commonHistRes->GetHistIntFlowRP(); 
     } else if(ptEta == "Pt")
       {
        hist = commonHistRes->GetHistDiffFlowPtRP();      
       } else if(ptEta == "Eta")
         {
          hist = commonHistRes->GetHistDiffFlowEtaRP();
         }
    } else if(rfRpPoi == "POI")
      {
       if(ptEta == "")
       {
        hist = commonHistRes->GetHistIntFlowPOI(); 
       } else if(ptEta == "Pt") 
         {
          hist = commonHistRes->GetHistDiffFlowPtPOI();
         } else if(ptEta == "Eta")
           {
            hist = commonHistRes->GetHistDiffFlowEtaPOI();
           }  
      }
 }

 return hist;

} // end of TH1D* GetResultHistogram(TString method, TString rfRpPoi, TString ptEta="")

// ===========================================================================================

void GlobalSettings()
{
 // Settings which will affect all plots.
 
 gROOT->SetStyle("Plain"); // default color is white instead of gray
 gStyle->SetOptStat(0); // remove stat. box from all histos
 
} // end of void GlobalSettings()

// ===========================================================================================

void LoadLibrariesCFR(const libModes analysisMode) {
  
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
  if (analysisMode==mLocal) {
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
    //cerr<<"libPWG2flowCommon loaded ..."<<endl;
    
  }
  
  else if (analysisMode==mLocalSource) {
    
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
  
} // end of void LoadLibrariesCFR(const libModes analysisMode) 

