enum libModes {mLocal,mLocalSource};

void fqd(TString type="", Int_t mode=mLocal)
{
 // macro to refit q-distribution
 
 Bool_t plotResults = kTRUE; // plot on the screen q-distribution and resulting fitting functions (for non-fixed sigma^2 and fixed sigma^2)
 Double_t treshold = 5; // first and last bin taken for fitting are determined as the first and last bin with more than 'treshold' number of entries
 Double_t vStart = 0.05; // starting value for v fit
 Double_t vMin = 0.0; // lower bound for v fit 
 Double_t vMax = 0.10; // upper bound for v fit  
 Double_t Sigma2Start = 0.75; // starting value for sigma^2 fit 
 Double_t Sigma2Min = 0.5; // lower bound for sigma^2 fit (according to theorists must be >= 0.5)   
 Double_t Sigma2Max = 2.5; // upper bound for sigma^2 fit
 
 if(Sigma2Min < 0.5)
 {
  cout<<"WARNING: According to theorists Sigma2Min >= 0.5 ..."<<endl;
  exit(0);
 }
 
 // load needed libraries:                       
 LoadLibrariesFQD(mode); 
 
 // access the path of current diretory:
 TString pwd(gSystem->pwd());
 pwd+="/";
 
 // FQD:
 TString outputFileNameFQD("outputFQDanalysis");
 (outputFileNameFQD+=(type.Data()))+=(".root");
 TFile *outputFileFQD = NULL;
 TList *outputListFQD = NULL;
 TString pwdFQD = pwd.Data();
 pwdFQD += outputFileNameFQD;
 if(gSystem->AccessPathName(pwdFQD.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have an output file "<<pwdFQD.Data()<<" !!!!"<<endl;
 } else 
   {
    outputFileFQD = TFile::Open(pwdFQD.Data(),"READ");
    if(outputFileFQD) 
    {
     outputFileFQD->GetObject("cobjFQD",outputListFQD);
     outputFileFQD->Close();
    } 
   }
 if(outputListFQD)
 {
  AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
  fqd->GetOutputHistograms(outputListFQD);
  // set fitting parameters:
  fqd->SetPlotResults(plotResults);
  fqd->SetTreshold(treshold);
  fqd->SetvStart(vStart);  
  fqd->SetvMin(vMin);
  fqd->SetvMax(vMax);
  fqd->SetSigma2Start(Sigma2Start);  
  fqd->SetSigma2Min(Sigma2Min);
  fqd->SetSigma2Max(Sigma2Max);
  // store fitting parameters:  
  fqd->StoreFittingParameters();
  // redo fit:
  fqd->Finish(kTRUE);
  
  // save the new results of fitting in the FQD output file: 
  TString newOutputFileNameFQD("outputFQDanalysis");
  (newOutputFileNameFQD+=(type.Data()))+=(".root");
  TString pwdFinalFQD=pwd.Data();
  pwdFinalFQD+=newOutputFileNameFQD;
  TFile *newOutputFQD = new TFile(pwdFinalFQD.Data(),"RECREATE");
  outputListFQD->SetName("cobjFQD");
  outputListFQD->Write(outputListFQD->GetName(),TObject::kSingleKey);
  newOutputFQD->Close();
  delete fqd;
  delete newOutputFQD;
 } else 
   {
    cout<<"WARNING: outputListFQD is NULL !!!!"<<endl;
   }              
 
} // end of void fqd(TString type="", Int_t mode=mLocal)

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
