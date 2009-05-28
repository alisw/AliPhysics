enum libModes {mLocal,mLocalSource};

void redoFinish(TString type="", Int_t mode=mLocal)
{
 // type:  type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //        (if type="" output files are from MC simulation (default))
 // mode:  if mode=mLocal analyze data on your computer using aliroot
 //        if mode=mLocalSource analyze data on your computer using root + source files
 
 // load needed libraries:                       
 LoadLibrariesRDF(mode);  
 
 // access the path of current diretory:
 TString pwd(gSystem->pwd());
 pwd+="/";
 
 // access the merged output files and redo the flow analysis on them:
 // MCEP:
 TString mergedFileNameMCEP("mergedMCEPanalysis");
 (mergedFileNameMCEP+=(type.Data()))+=(".root");
 TFile *mergedFileMCEP = NULL;
 TList *mergedListMCEP = NULL;
 TString pwdMCEP=pwd.Data();
 pwdMCEP+=mergedFileNameMCEP;
 if(gSystem->AccessPathName(pwdMCEP.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdMCEP.Data()<<endl;
 } else 
   {
    mergedFileMCEP = TFile::Open(pwdMCEP.Data(),"READ");
    if(mergedFileMCEP) 
    {
     mergedFileMCEP->GetObject("cobjMCEP",mergedListMCEP);
     mergedFileMCEP->Close();
    } 
   }
 if(mergedListMCEP)
 {
  AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
  mcep->GetOutputHistograms(mergedListMCEP);
  mcep->Finish();
  // save the final results for MCEP in final output file: 
  TString finalOutputFileNameMCEP("outputMCEPanalysis");
  (finalOutputFileNameMCEP+=(type.Data()))+=(".root");
  TString pwdFinalMCEP=pwd.Data();
  pwdFinalMCEP+=finalOutputFileNameMCEP;
  TFile *finalOutputMCEP = new TFile(pwdFinalMCEP.Data(),"NEW");
  mergedListMCEP->SetName("cobjMCEP");
  mergedListMCEP->Write(mergedListMCEP->GetName(),TObject::kSingleKey);
  finalOutputMCEP->Close();
 } else 
   {
    cout<<"WARNING: mergedListMCEP is NULL !!!!"<<endl;
   }
   
 // SP:
 TString mergedFileNameSP("mergedSPanalysis");
 (mergedFileNameSP+=(type.Data()))+=(".root");
 TFile *mergedFileSP = NULL;
 TList *mergedListSP = NULL;
 TString pwdSP=pwd.Data();
 pwdSP+=mergedFileNameSP;
 if(gSystem->AccessPathName(pwdSP.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdSP.Data()<<endl;
 } else 
   {
    mergedFileSP = TFile::Open(pwdSP.Data(),"READ");
    if(mergedFileSP) 
    {
     mergedFileSP->GetObject("cobjSP",mergedListSP);
     mergedFileSP->Close();
    } 
   }
 if(mergedListSP)
 {
  AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
  sp->GetOutputHistograms(mergedListSP);
  sp->Finish();
  // save the final results for SP in final output file: 
  TString finalOutputFileNameSP("outputSPanalysis");
  (finalOutputFileNameSP+=(type.Data()))+=(".root");
  TString pwdFinalSP=pwd.Data();
  pwdFinalSP+=finalOutputFileNameSP;
  TFile *finalOutputSP = new TFile(pwdFinalSP.Data(),"NEW");
  mergedListSP->SetName("cobjSP");
  mergedListSP->Write(mergedListSP->GetName(),TObject::kSingleKey);
  finalOutputSP->Close();
 } else 
   {
    cout<<"WARNING: mergedListSP is NULL !!!!"<<endl;
   }     
   
 // GFC:
 TString mergedFileNameGFC("mergedGFCanalysis");
 (mergedFileNameGFC+=(type.Data()))+=(".root");
 TFile *mergedFileGFC = NULL;
 TList *mergedListGFC = NULL;
 TString pwdGFC=pwd.Data();
 pwdGFC+=mergedFileNameGFC;
 if(gSystem->AccessPathName(pwdGFC.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdGFC.Data()<<endl;
 } else 
   {
    mergedFileGFC = TFile::Open(pwdGFC.Data(),"READ");
    if(mergedFileGFC) 
    {
     mergedFileGFC->GetObject("cobjGFC",mergedListGFC);
     mergedFileGFC->Close();
    } 
   }
 if(mergedListGFC)
 {
  AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
  gfc->GetOutputHistograms(mergedListGFC);
  gfc->Finish();
  // save the final results for GFC in final output file: 
  TString finalOutputFileNameGFC("outputGFCanalysis");
  (finalOutputFileNameGFC+=(type.Data()))+=(".root");
  TString pwdFinalGFC=pwd.Data();
  pwdFinalGFC+=finalOutputFileNameGFC;
  TFile *finalOutputGFC = new TFile(pwdFinalGFC.Data(),"NEW");
  mergedListGFC->SetName("cobjGFC");
  mergedListGFC->Write(mergedListGFC->GetName(),TObject::kSingleKey);
  finalOutputGFC->Close();
 } else 
   {
    cout<<"WARNING: mergedListGFC is NULL !!!!"<<endl;
   }      
   
 // QC:
 TString mergedFileNameQC("mergedQCanalysis");
 (mergedFileNameQC+=(type.Data()))+=(".root");
 TFile *mergedFileQC = NULL;
 TList *mergedListQC = NULL;
 TString pwdQC=pwd.Data();
 pwdQC+=mergedFileNameQC;
 if(gSystem->AccessPathName(pwdQC.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdQC.Data()<<endl;
 } else 
   {
    mergedFileQC = TFile::Open(pwdQC.Data(),"READ");
    if(mergedFileQC) 
    {
     mergedFileQC->GetObject("cobjQC",mergedListQC);
     mergedFileQC->Close();
    } 
   }
 if(mergedListQC)
 {
  AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
  qc->GetOutputHistograms(mergedListQC);
  qc->Finish();
  // save the final results for QC in final output file: 
  TString finalOutputFileNameQC("outputQCanalysis");
  (finalOutputFileNameQC+=(type.Data()))+=(".root");
  TString pwdFinalQC=pwd.Data();
  pwdFinalQC+=finalOutputFileNameQC;
  TFile *finalOutputQC = new TFile(pwdFinalQC.Data(),"NEW");
  mergedListQC->SetName("cobjQC");
  mergedListQC->Write(mergedListQC->GetName(),TObject::kSingleKey);
  finalOutputQC->Close();
 } else 
   {
    cout<<"WARNING: mergedListQC is NULL !!!!"<<endl;
   }  
   
 // FQD:
 TString mergedFileNameFQD("mergedFQDanalysis");
 (mergedFileNameFQD+=(type.Data()))+=(".root");
 TFile *mergedFileFQD = NULL;
 TList *mergedListFQD = NULL;
 TString pwdFQD=pwd.Data();
 pwdFQD+=mergedFileNameFQD;
 if(gSystem->AccessPathName(pwdFQD.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdFQD.Data()<<endl;
 } else 
   {
    mergedFileFQD = TFile::Open(pwdFQD.Data(),"READ");
    if(mergedFileFQD) 
    {
     mergedFileFQD->GetObject("cobjFQD",mergedListFQD);
     mergedFileFQD->Close();
    } 
   }
 if(mergedListFQD)
 {
  AliFittingQDistribution* fqd = new AliFittingQDistribution();
  fqd->GetOutputHistograms(mergedListFQD);
  fqd->Finish();
  // save the final results for FQD in final output file: 
  TString finalOutputFileNameFQD("outputFQDanalysis");
  (finalOutputFileNameFQD+=(type.Data()))+=(".root");
  TString pwdFinalFQD=pwd.Data();
  pwdFinalFQD+=finalOutputFileNameFQD;
  TFile *finalOutputFQD = new TFile(pwdFinalFQD.Data(),"NEW");
  mergedListFQD->SetName("cobjFQD");
  mergedListFQD->Write(mergedListFQD->GetName(),TObject::kSingleKey);
  finalOutputFQD->Close();
 } else 
   {
    cout<<"WARNING: mergedListFQD is NULL !!!!"<<endl;
   }                               

 // LYZ1:
 TString mergedFileNameLYZ1("mergedLYZ1analysis");
 (mergedFileNameLYZ1+=(type.Data()))+=(".root");
 TFile *mergedFileLYZ1 = NULL;
 TList *mergedListLYZ1 = NULL;
 TString pwdLYZ1=pwd.Data();
 pwdLYZ1+=mergedFileNameLYZ1;
 if(gSystem->AccessPathName(pwdLYZ1.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdLYZ1.Data()<<endl;
 } else 
   {
    mergedFileLYZ1 = TFile::Open(pwdLYZ1.Data(),"READ");
    if(mergedFileLYZ1) 
    {
     mergedFileLYZ1->GetObject("cobjLYZ1",mergedListLYZ1);
     mergedFileLYZ1->Close();
    } 
   }
 if(mergedListLYZ1)
 {
  AliFlowAnalysisWithLeeYangZeros* lyz1 = new AliFlowAnalysisWithLeeYangZeros();
  lyz1->SetFirstRun(kTRUE);   
  lyz1->SetUseSum(kTRUE);       
  lyz1->GetOutputHistograms(mergedListLYZ1);
  lyz1->Finish();
  // save the final results for LYZ1 in final output file: 
  TString finalOutputFileNameLYZ1("outputLYZ1analysis");
  (finalOutputFileNameLYZ1+=(type.Data()))+=(".root");
  TString pwdFinalLYZ1=pwd.Data();
  pwdFinalLYZ1+=finalOutputFileNameLYZ1;
  TFile *finalOutputLYZ1 = new TFile(pwdFinalLYZ1.Data(),"NEW");
  mergedListLYZ1->SetName("cobjLYZ1");
  mergedListLYZ1->Write(mergedListLYZ1->GetName(),TObject::kSingleKey);
  finalOutputLYZ1->Close();
 } else 
   {
    cout<<"WARNING: mergedListLYZ1 is NULL !!!!"<<endl;
   }                               

 
} // end of void reCallFinish(TString type="", Int_t mode=mLocal)
 
void LoadLibrariesRDF(const libModes mode) {
  
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
  
} // end of void LoadLibrariesRDF(const libModes mode)

 
