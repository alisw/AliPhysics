//==========================================================================//
// Macro to refit the q-distribution in FQD method. By using this macro the // 
// fitting parameters for q-distribution can be tuned and correspondingly   //
// the estimate for reference flow harmonic obtained from FQD method can be //
// significantly improved.                                                  // 
//==========================================================================//

// Fitting parameters for q-distribution:
Double_t vStart = 0.05; // starting value for v fit
Double_t vMin = 0.0; // lower bound for v fit 
Double_t vMax = 0.25; // upper bound for v fit  
Double_t sigma2Start = 0.75; // starting value for sigma^2 fit 
Double_t sigma2Min = 0.5; // lower bound for sigma^2 fit (according to theorists must be >= 0.5)   
Double_t sigma2Max = 2.5; // upper bound for sigma^2 fit
Double_t treshold = 5.; // first and last bin taken for fitting are determined as the first and last bin with more than 'treshold' number of entries
Bool_t finalResultIsFromSigma2Fitted = kTRUE; // final saved result is obtained with sigma^2 fitted (kTRUE) or sigma^2 fixed (kFALSE)
Bool_t printOnTheScreen = kTRUE; // print or not the final results on the screen
Bool_t plotResults = kTRUE; // plot on the screen q-distribution and resulting fitting functions (for non-fixed sigma^2 and fixed sigma^2)

enum libModes {mLocal,mLocalSource};

TFile *commonOutputFile = NULL; // common output file "AnalysisResults.root"
TDirectoryFile *dirFileFQD = NULL; // FQD's TDirectoryFile in "AnalysisResults.root"
TList *outputListFQD = NULL; // output list holding all FQD objects

void fqd(TString analysisType="", Int_t analysisMode=mLocal)
{
 // 1. analysisType: "ESD", "AOD", "ESDMC0", "ESDMC1"; for Monte Carlo and 'on the fly' use simply "";
 // 2. analysisMode: if analysisMode = mLocal -> analyze data on your computer using aliroot
 //                  if analysisMode = mLocalSource -> analyze data on your computer using root + source files 
  
 // Cross-check if the user's settings make sense:
 CrossCheckSettings();
 // Load needed libraries:                       
 LoadLibrariesFQD(analysisMode); 
 // Get output list which holds all FQD objects:
 GetOutputList(analysisType); 
 // Redo fit of q-distribution:
 RedoFit();
 // Plot q-distribution and resulting fitting functions:
 if(plotResults) Plot(); 
 // Save the new results in the common output file: 
 Save();
  
} // end of void fqd(TString type="ESD", Int_t mode=mLocal)   

// =========================================================================================== 

void RedoFit()
{
 // Redo the fit of q-distribution.

 AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
 fqd->GetOutputHistograms(outputListFQD);
 // Set new fitting parameters:
 fqd->SetvStart(vStart);  
 fqd->SetvMin(vMin);
 fqd->SetvMax(vMax);
 fqd->SetSigma2Start(sigma2Start);  
 fqd->SetSigma2Min(sigma2Min);
 fqd->SetSigma2Max(sigma2Max);
 fqd->SetTreshold(treshold);
 fqd->SetFinalResultIsFromSigma2Fitted(finalResultIsFromSigma2Fitted);
 fqd->SetPrintOnTheScreen(printOnTheScreen); 
 // Save new fitting parameters:  
 fqd->StoreFittingParameters();
 // Redo fit:
 fqd->Finish(kTRUE);
 
} // end of void RedoFit(TList *outputList)
 
// =========================================================================================== 
 
void GetOutputList(TString analysisType)
{
 // Get output list which holds all FQD objects.
  
 // Access common output file:
 TString outputFileName = "AnalysisResults.root"; // name of the common output file
 commonOutputFile = AccessOutputFile(outputFileName);
 
 // Access from common output file the TDirectoryFile for FQD method
 // and from it the output list holding all objects:
 GetListWithHistograms(commonOutputFile,analysisType,"FQD");
   
} // end of TList* GetOutputList(TString analysisType)
 
// ===========================================================================================

void Plot()
{
 // Plot q-distribution and resulting fitting functions.
 
 gROOT->SetStyle("Plain"); // default color is white instead of gray
 gStyle->SetOptStat(0); // remove stat. box from all histos
 fLegend = new TLegend(0.6,0.55,0.85,0.7); 
 // q-distribution:
 TH1D *qDistribution = dynamic_cast<TH1D*>(outputListFQD->FindObject("fqDistribution"));
 Cosmetics(qDistribution);
 fLegend->AddEntry(qDistribution,"q-distribution","f");
 qDistribution->Draw();
 // resulting fitting functions:
 TF1 *fittingFun[2] = {NULL};
 TString sigmaFlag[2] = {"#sigma^{2} not fitted","#sigma^{2} fitted"};
 Double_t lineColors[2] = {kBlue,kRed};
 for(Int_t s2F=0;s2F<2;s2F++)
 {
  fittingFun[s2F] = dynamic_cast<TF1*>(outputListFQD->FindObject(Form("fFittingFunction, %s",sigmaFlag[s2F].Data())));
  fittingFun[s2F]->SetLineColor(lineColors[s2F]);
  fittingFun[s2F]->Draw("SAME");
  fLegend->AddEntry(fittingFun[s2F]->GetName(),sigmaFlag[s2F].Data(),"l");
 }
 fLegend->Draw("SAME"); 
  
} // end of Plot()

// ===========================================================================================

void Save()
{
 // Save the new results of fitting for FQD method in the common output file.
 
 dirFileFQD->Add(outputListFQD,kTRUE);
 dirFileFQD->Write(dirFileFQD->GetName(),TObject::kSingleKey+TObject::kOverwrite);
 //delete commonOutputFile;

} // end of void Save()

// ===========================================================================================

void Cosmetics(TH1D *hist)
{
 // Set cosmetics for the q-distribution.
 
 Int_t firstNonEmptyBin = hist->FindFirstBinAbove(0);
 Double_t lowRange = hist->GetBinLowEdge(firstNonEmptyBin);
 Int_t lastNonEmptyBin = hist->FindLastBinAbove(0);
 Double_t upperRange = hist->GetBinLowEdge(lastNonEmptyBin+10);
 hist->GetXaxis()->SetRangeUser(lowRange,upperRange);  
 hist->SetFillColor(16);  
 hist->SetTitle("Fitted q-distribution");
 
} // end of void Cosmetics(TH1D *hist)

// ===========================================================================================

TFile* AccessOutputFile(TString outputFileName)
{
 // Access the common output file.
    
 TFile *outputFile = NULL;
 if(!(gSystem->AccessPathName(Form("%s%s%s",gSystem->pwd(),"/",outputFileName.Data()),kFileExists)))
 {
  outputFile = TFile::Open(outputFileName.Data(),"UPDATE");
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
 
void GetListWithHistograms(TFile *outputFile, TString analysisType, TString methodName)
{
 // Access from common output file the TDirectoryFile for FQD method and from it
 // the output list holding all objects:

 TString fileName = ""; 
 TString listName = ""; 
 // Form a file name: 
 fileName+="output";
 fileName+=methodName.Data();
 fileName+="analysis";
 fileName+=analysisType.Data();
 // Access this file:
 dirFileFQD = (TDirectoryFile*)outputFile->FindObjectAny(fileName.Data());
 // Form a list name:
 listName+="cobj";
 listName+=methodName.Data();
 // Access this list:
 if(dirFileFQD)
 {
  if(!(dirFileFQD->GetNkeys() == 0))
  {
   dirFileFQD->GetObject(listName.Data(),outputListFQD);
  } else
    {
     TString temp = listName.Data();
     cout<<"WARNING: Couldn't access a list holding histograms for "<<temp.Remove(0,4).Data()<<" method !!!!"<<endl;     
    } 
 } else 
   {
    cout<<"WARNING: Couldn't find a file "<<fileName.Data()<<".root !!!!"<<endl;
   }
     
 if(!outputListFQD) 
 {
  cout<<endl;
  cout<<"WARNING: Couldn't access the output list "<<listName.Data()<<" !!!!"<<endl;
  cout<<endl;
  exit(0);
 }      
     
} // end of void GetListWithHistograms(TFile *outputFile, TString analysisType, TString methodName)
 
//===========================================================================================

void CrossCheckSettings() 
{
 // Check in this method if the settings make sense.
 
 if(sigma2Min<0.5)
 {
  cout<<endl; 
  cout<<"WARNING: According to theorists sigma2Min >= 0.5 !!!!"<<endl;
  cout<<endl; 
  exit(0);
 }
 
} // end of void CrossCheckSettings()

//===========================================================================================

void LoadLibrariesFQD(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");
  
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
  gSystem->Load("libTree.so");

  // for AliRoot
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libPWG2flowCommon.so");
  cerr<<"libPWG2flowCommon.so loaded ..."<<endl;
  
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
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
    
    // Cuts
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    
    // Output histosgrams
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
       
    cout << "finished loading macros!" << endl;  
    
  } // end of else if (mode==mLocalSource) 
  
} // end of void LoadLibrariesFQD(const libModes mode)
