// ... 

enum libModes {mLocal,mLocalSource};

const Int_t nFilesMax = 10; // number of files to be accessed to estimate spread

Bool_t showPlotForReferenceFlow = kTRUE;
Bool_t showPlotForIntegratedFlowRP = kTRUE;
Bool_t showPlotForIntegratedFlowPOI = kTRUE;
Bool_t showErrorOnMergedResult = kTRUE; 

const Int_t nMethods = 13;
TString method[nMethods] = {"MCEP","SP","2,GFC","2,QC","4,GFC","4,QC","6,GFC","6,QC","8,GFC","8,QC","FQD","LYZ1SUM","LYZ1PROD"};
Int_t methodMarkerStyle[nMethods] = {21,21,21,21,21,21,21,21,21,21,21,21,21};
Int_t methodMarkerColor[nMethods] = {kGray+1,kViolet-8,kBlue-9,kRed-7,kBlue-9,kRed-7,kBlue-9,kRed-7,kBlue-9,kRed-7,kOrange-8,kYellow-5,kYellow-2};
Int_t methodMeshColor[nMethods] = {kGray,kViolet-9,kBlue-10,kRed-10,kBlue-10,kRed-10,kBlue-10,kRed-10,kBlue-10,kRed-10,kOrange-9,kYellow-8,kYellow-6};

/*
const Int_t nMethods = 4;
TString method[nMethods] = {"2,QC","4,QC","6,QC","8,QC"};
Int_t methodMarkerStyle[nMethods] = {21,21,21,21};
Int_t methodMarkerColor[nMethods] = {kRed-7,kRed-7,kRed-7,kRed-7};
Int_t methodMeshColor[nMethods] = {kRed-10,kRed-10,kRed-10,kRed-10};
*/

void showSpread(TString type="", Int_t mode=mLocal)
{
 // type: type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //       (if type="" output files are from MC simulation (default))
 // mode: if mode = mLocal: analyze data on your computer using aliroot
 //       if mode = mLocalSource: analyze data on your computer using root + source files 

 // Cross-check if the user's settings make sense:
 CrossCheckUserSettings();
 
 // Load needed libraries:                       
 LoadLibrariesSS(mode);  
  
 // Output file name:
 TString outputFileName = "AnalysisResults.root"; 
 
 // Labels for reference flow, integrated flow of RPs and of POIs:
 TString label[3] = {"","RP","POI"};
  
 // Standard magic:
 TString *baseDirPath = new TString(gSystem->pwd());
 TSystemDirectory *baseDir = new TSystemDirectory(".",baseDirPath->Data());          
 TList *listOfFilesInBaseDir = baseDir->GetListOfFiles();
 TStopwatch timer;
 timer.Start();
 // listOfFilesInBaseDir->Print();
 Int_t nFiles = listOfFilesInBaseDir->GetEntries();
 Int_t fileCounter = 0;
 Double_t result[nMethods][3][nFilesMax] = {{{0.}}}; // [3 = "", "RP" or "POI"]
 Double_t error[nMethods][3][nFilesMax] = {{{0.}}}; // [3 = "", "RP" or "POI"]
 Double_t resultMinMax[nMethods][3][2] = {{{0.}}}; // [3 = "", "RP" or "POI"], [2 = min value, max value]
 for(Int_t m=0;m<nMethods;m++)
 {
  for(Int_t l=0;l<3;l++)
  {
   resultMinMax[m][l][0] = 44.;
   resultMinMax[m][l][1] = -44.;
  }
 } 
 Double_t styleHistMinMax[3][2] = {{0.}}; // [3 = "", "RP" or "POI"], [2 = min value, max value]
 for(Int_t l=0;l<3;l++)
 {
  styleHistMinMax[l][0] = 44.;
  styleHistMinMax[l][1] = -44.;
 }  
 cout<<endl;
 for(Int_t iFile=0;iFile<nFiles;iFile++)
 {
  TSystemFile *currentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
  // Consider only subdirectories: 
  if(!currentFile || 
     !currentFile->IsDirectory() || 
     strcmp(currentFile->GetName(), ".") == 0 || 
    strcmp(currentFile->GetName(), "..") == 0) continue; 
  // Accessing the output file "AnalysisResults.root" in current subdirectory: 
  TString currentSubDirName = baseDirPath->Data();
  (currentSubDirName+="/")+=currentFile->GetName();
  currentSubDirName+="/";
  TString fileName = currentSubDirName; 
  fileName+=outputFileName.Data();
  if(!(gSystem->AccessPathName(fileName.Data(),kFileExists)))
  {
   TFile *file = NULL; 
   file = TFile::Open(fileName.Data(),"READ");
   for(Int_t m=0;m<nMethods;m++)
   {
    // Access from common output file the output file for particular method:
    TDirectoryFile *methodFile = NULL;
    methodFile = GetMethodFile(file,method[m],type);  
    for(Int_t l=0;l<3;l++)
    {
     TH1D *histResult = NULL;
     histResult = GetHistogramWithResult(method[m],label[l],methodFile);
     if(histResult)
     {
      // Access the results:
      result[m][l][fileCounter] = histResult->GetBinContent(1); 
      // Access the errors:
      error[m][l][fileCounter] = histResult->GetBinError(1);
      if(TMath::Abs(result[m][l][fileCounter])>pow(10.,-6.)) // take into account only if != 0 (to be improved - special care for < 0 is required)
      {
       // Establish min and max values for results:
       if(resultMinMax[m][l][0] > result[m][l][fileCounter]) resultMinMax[m][l][0] = result[m][l][fileCounter]; // min value
       if(resultMinMax[m][l][1] < result[m][l][fileCounter]) resultMinMax[m][l][1] = result[m][l][fileCounter]; // max value
       // Establish min and max values for style histograms:
       if(styleHistMinMax[l][0] > result[m][l][fileCounter]) styleHistMinMax[l][0] = result[m][l][fileCounter]; // min value
       if(styleHistMinMax[l][1] < result[m][l][fileCounter]) styleHistMinMax[l][1] = result[m][l][fileCounter]; // max value    
      }
     }
    } // end of for(Int_t l=0;l<3;l++)
   } // end of for(Int_t m=0;m<nMethods;m++)
   if(file) file->Close();
   fileCounter++;
   //if(fileCounter%10==0)
   //{
    cout<<Form("Accessed %d files \"AnalysisResults.root\" so far....",fileCounter)<<"\r"<<flush;
   //}   
  } 
  if(fileCounter == nFilesMax) break;
 } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)   
  
 cout<<Form("Accessed %d files \"AnalysisResults.root\" in total to estimate spread. ",fileCounter)<<endl;
 cout<<endl;
 const Int_t nFilesFinal = fileCounter;

 // Make for each method graph holding results:
 TGraph *methodGraph[nMethods][3] = {{NULL}}; // [3 = "", "RP" or "POI"] 
 Double_t x[nMethods][nFilesFinal] = {{0.}};
 for(Int_t m=0;m<nMethods;m++)
 {
  for(Int_t f=0;f<nFilesFinal;f++)
  {
   x[m][f]=m+0.5;
  } 
 }
 for(Int_t m=0;m<nMethods;m++)
 {
  for(Int_t l=0;l<3;l++)
  {
   methodGraph[m][l] = new TGraph(nFilesFinal,x[m],result[m][l]);
   methodGraph[m][l]->SetMarkerStyle(methodMarkerStyle[m]);
   methodGraph[m][l]->SetMarkerColor(methodMarkerColor[m]);
  } // end of for(Int_t l=0;l<3;l++)
 } // for(Int_t m=0;m<nMethods;m++)
 
 // Make for each method coloured mesh out of min and max values:
 Double_t meshWidth = 0.25;
 TGraph *methodMesh[nMethods][3] = {{NULL}}; // [3 = "", "RP" or "POI"]
 for(Int_t m=0;m<nMethods;m++)
 {
  for(Int_t l=0;l<3;l++)
  {
   if(resultMinMax[m][l][0]<44. && resultMinMax[m][l][1]>-44.)
   {
    methodMesh[m][l] = new TGraph(5);
    methodMesh[m][l]->SetPoint(0,(m+1-0.5)-meshWidth,resultMinMax[m][l][0]);
    methodMesh[m][l]->SetPoint(1,(m+1-0.5)+meshWidth,resultMinMax[m][l][0]);
    methodMesh[m][l]->SetPoint(2,(m+1-0.5)+meshWidth,resultMinMax[m][l][1]);
    methodMesh[m][l]->SetPoint(3,(m+1-0.5)-meshWidth,resultMinMax[m][l][1]);
    methodMesh[m][l]->SetPoint(4,(m+1-0.5)-meshWidth,resultMinMax[m][l][0]);    
    methodMesh[m][l]->SetFillStyle(1001);
    methodMesh[m][l]->SetFillColor(methodMeshColor[m]);
   } 
  }
 }
  
 // Access for each method the results from the merged, large statistics file:
 Double_t resultMerged[nMethods][3] = {{0.}}; // [3 = "", "RP" or "POI"]
 TString mergedFileName = Form("%s%s%s",gSystem->pwd(),"/",outputFileName.Data());
 TFile *mergedFile = NULL; 
 mergedFile = TFile::Open(mergedFileName.Data(),"READ"); 
 for(Int_t m=0;m<nMethods;m++)
 {
  TDirectoryFile *methodFile = NULL;
  if(!(gSystem->AccessPathName(fileName.Data(),kFileExists)))
  {
   if(mergedFile) methodFile = GetMethodFile(mergedFile,method[m],type);  
  } else
    {
     cout<<"WARNING: Couldn't find the merged, large statistics file "<<endl;
     cout<<"         "<<fileName.Data()<<endl;
     cout<<"         in directory "<<gSystem->pwd()<<" !!!!"<<endl;     
     cout<<"         Use macros mergeOuput.C and redoFinish.C to get it."<<endl;
     cout<<endl;
     break;
    }
  for(Int_t l=0;l<3;l++)
  {
   TH1D *histResult = NULL;
   if(methodFile)
   { 
    histResult = GetHistogramWithResult(method[m],label[l],methodFile);
   }
   if(histResult)
   {
    // Access the results from the merged, large statistics file:
    resultMerged[m][l] = histResult->GetBinContent(1);
   } // end of for(Int_t l=0;l<3;l++)
  }
 } // end of for(Int_t m=0;m<nMethods;m++)  
 if(mergedFile) mergedFile->Close();
        
 // Make for each method graph holding results from the merged, large statistics file
 // and the errors from the randomly chosen small statistics file:
 TGraphErrors *methodMergedGraph[nMethods][3] = {{NULL}}; // [3 = "", "RP" or "POI"] 
 Double_t xMerged[nMethods] = {0.};
 for(Int_t m=0;m<nMethods;m++)
 {
  xMerged[m]=m+0.5; 
 }
 // Select randomly small statistics file:
 gRandom->SetSeed((UInt_t) (4400*timer.RealTime()/fileCounter));
 Int_t randomFile = (Int_t)gRandom->Uniform(0,fileCounter); 
 for(Int_t m=0;m<nMethods;m++)
 {
  for(Int_t l=0;l<3;l++)
  {
   methodMergedGraph[m][l] = new TGraphErrors(1);
   methodMergedGraph[m][l]->SetPoint(0,xMerged[m],resultMerged[m][l]);
   if(showErrorOnMergedResult) methodMergedGraph[m][l]->SetPointError(0,0.,error[m][l][randomFile]);
   methodMergedGraph[m][l]->SetMarkerStyle(25);
   methodMergedGraph[m][l]->SetMarkerColor(kBlack);
  } // end of for(Int_t l=0;l<3;l++)
 } // for(Int_t m=0;m<nMethods;m++)
 
 // Final drawing:
 gROOT->SetStyle("Plain"); // removing default gray color and setting white instead
 gStyle->SetOptStat(0); // removing statistics box from all histograms
 Bool_t showPlot[3] = {showPlotForReferenceFlow,showPlotForIntegratedFlowRP,showPlotForIntegratedFlowPOI};
 TString title[3] = {"Reference Flow","Integrated Flow (RP)","Integrated Flow (POI)"}
 TCanvas *canvas[3] = {NULL}; // [3 = "", "RP" or "POI"]
 for(Int_t l=0;l<3;l++)
 {
  if(showPlot[l])
  {
   canvas[l] = new TCanvas(Form("%s",title[l].Data()),Form("%s",title[l].Data()));
   TH1D *styleHist = StyleHist(title[l]);
   styleHist->SetMinimum(0.99*styleHistMinMax[l][0]);
   styleHist->SetMaximum(1.01*styleHistMinMax[l][1]);
   styleHist->Draw();
   for(Int_t m=0;m<nMethods;m++)
   {
    if(methodMesh[m][l]) methodMesh[m][l]->Draw("lfsame");
    if(methodGraph[m][l]) methodGraph[m][l]->Draw("psame");  
    if(TMath::Abs(*(methodMergedGraph[m][l]->GetY()))>pow(10.,-6.)) // draw only if not == 0.
    {
     methodMergedGraph[m][l]->Draw("psame");
    } 
   } // end of for(Int_t m=0;m<nMethods;m++)
  } // end of if(showPlot[l])
 } // end of for(Int_t l=0;l<3;l++)
 
 timer.Stop();
 timer.Print(); 
 cout<<endl;
 
} // end of void showSpread(TString type="", Int_t mode=mLocal)

// =============================================================================================

TH1D *StyleHist(TString histTitle)
{
 // Make the style histogram.
 Int_t n = 2; // harmonic (to be improved - access this from common control hist)
 TH1D *styleHist = new TH1D(Form("%s",histTitle.Data()),Form("%s",histTitle.Data()),nMethods,0,nMethods);
 styleHist->GetYaxis()->SetTickLength(0.01);
 for(Int_t m=0;m<nMethods;m++)
 {
  styleHist->GetXaxis()->SetBinLabel(m+1,Form("v_{%d}{%s}",n,method[m].Data()));
  if(method[m]=="LYZ1SUM" || method[m]=="LYZ2SUM")
  {
   styleHist->GetXaxis()->SetBinLabel(m+1,Form("v_{%d}{%s}",n,"LYZ,sum"));
  } 
  else if(method[m]=="LYZ1PROD" || method[m]=="LYZ2PROD")
  {
   styleHist->GetXaxis()->SetBinLabel(m+1,Form("v_{%d}{%s}",n,"LYZ,prod"));
  } 
 } // end of for(Int_t m=0;m<nMethods;m++)
 
 return styleHist;
}

// =============================================================================================

TDirectoryFile* GetMethodFile(TFile *commonFile, TString method, TString type)
{
 // Form a file name for each method:
 TString methodFileName = "output";
 if(method.Contains("GFC"))
 {
  methodFileName+="GFC";
 } else if(method.Contains("QC"))
   {
    methodFileName+="QC";
   } else
     {
      methodFileName+=method.Data();
     } 
 methodFileName+="analysis";
 methodFileName+=type.Data();
 TDirectoryFile *methodFile = NULL;
 if(commonFile)
 {
  methodFile = (TDirectoryFile*)commonFile->FindObjectAny(methodFileName.Data());
 } 
 
 return methodFile;

} // end of TDirectoryFile* AccessMethodFile(TString commonFile, TString method, TString type)

// =============================================================================================

TH1D* GetHistogramWithResult(TString method, TString label, TDirectoryFile *methodFile)
{
 // Access first the common list holding all output histograms:
 TList *methodList = NULL;
 if(method.Contains("GFC"))
 {
  methodFile->GetObject("cobjGFC",methodList);
 } else if(method.Contains("QC"))
   {
    methodFile->GetObject("cobjQC",methodList);
   } else
     {
      methodFile->GetObject(Form("cobj%s",method.Data()),methodList);
     } 
 // Access from the common list the needed histogram:
 if(methodList)
 {
  AliFlowCommonHistResults *commonHistRes = NULL; 
  if(!(method.Contains("GFC") || method.Contains("QC")))
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject(Form("AliFlowCommonHistResults%s",method.Data()));
  } 
  else if(method=="2,GFC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults2ndOrderGFC");    
  } 
  else if(method=="4,GFC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults4thOrderGFC");    
  } 
  else if(method=="6,GFC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults6thOrderGFC");    
  } 
  else if(method=="8,GFC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults8thOrderGFC");    
  }    
  else if(method=="2,QC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults2ndOrderQC");    
  } 
  else if(method=="4,QC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults4thOrderQC");    
  } 
  else if(method=="6,QC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults6thOrderQC");    
  } 
  else if(method=="8,QC")
  {
   commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject("AliFlowCommonHistResults8thOrderQC");    
  }  
 } // end of if(methodList)
 
 // Access histogram with results for reference flow or integrated flow of RPs or POIs:
 TH1D *hist = NULL;
 if(label=="")
 {
  hist = commonHistRes->GetHistIntFlow();
 } else if(label=="RP")
   {
    hist = commonHistRes->GetHistIntFlowRP(); 
   } else if(label=="POI")
     {
      hist = commonHistRes->GetHistIntFlowPOI(); 
     }

 if(hist) return hist;

} // end of TH1D *GetHistogramWithResult(TString method, TString rf_rp_poi, TString file);

// =============================================================================================

void CrossCheckUserSettings() 
{
 // Check in this method if the user settings make sense.
 
 if(nFilesMax<=0)
 {
  cout<<endl;
  cout<<"WARNING: nFilesMax must be a positive integer (not too large, though) !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(nFilesMax>44)
 { 
  cout<<endl;
  cout<<"WARNING: You may want to set nFilesMax to the smaller value."<<endl;
  cout<<"         Otherwise you might wait forever to see the plots."<<endl;
 } 
 
} // end of void CrossCheckUserSettings()

// =============================================================================================

void LoadLibrariesSS(const libModes mode) {
  
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
  //cerr<<"libPWG2flowCommon loaded ..."<<endl;
  
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
    
    // Output histograms
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
    
    // Functions needed for various methods
    gROOT->LoadMacro("AliFlowCommon/AliCumulantsFunctions.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZEventPlane.cxx+");
    
    // Flow Analysis code for various methods
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMCEventPlane.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithScalarProduct.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLYZEventPlane.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithLeeYangZeros.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithCumulants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithQCumulants.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithFittingQDistribution.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithMixedHarmonics.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowAnalysisWithNestedLoops.cxx+");          
    
    cout << "finished loading macros!" << endl;  
    
  } // end of else if (mode==mLocalSource) 
  
} // end of void LoadLibrariesSS(const libModes mode)
