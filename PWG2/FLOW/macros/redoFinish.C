// Macro redoFinish.C is used after macro mergeOutput.C has been used for merging. 
// Before using this macro first read the explanation at the beginning of macro mergeOutput.C. 

enum libModes {mLocal,mLocalSource};

void redoFinish(TString type="", Int_t mode=mLocal)
{
 // type: type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //       (if type="" output files are from MC simulation (default))
 // mode: if mode = mLocal: analyze data on your computer using aliroot
 //       if mode = mLocalSource: analyze data on your computer using root + source files

 TString mergedFileName = "mergedAnalysisResults.root"; // hardwired name of merged, large statistics file obtained with macro mergeOutput.C 
 TString outputFileName = "AnalysisResults.root"; // final output file name holding final results for large statistics sample
 
 const Int_t nMethods = 10;
 TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP"};
 
 // load needed libraries:                       
 LoadLibrariesRF(mode);  
  
 // access the merged, large statistics file obtained with macro mergeOutput.C:
 TString pwd(gSystem->pwd());
 pwd+="/";
 pwd+=mergedFileName.Data();
 TFile *mergedFile = NULL;
 if(gSystem->AccessPathName(pwd.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file:"<<endl;
  cout<<"         "<<pwd.Data()<<endl;
  exit(0);
 } else 
   {
    mergedFile = TFile::Open(pwd.Data(),"READ");
   }

 // access from mergedFile the merged files for each method and from them the lists holding histograms:
 TString fileName[nMethods]; 
 TDirectoryFile *dirFile[nMethods] = {NULL}; 
 TString listName[nMethods]; 
 TList *list[nMethods] = {NULL};
 for(Int_t i=0;i<nMethods;i++)
 {
  // form a file name for each method:
  fileName[i]+="output";
  fileName[i]+=method[i].Data();
  fileName[i]+="analysis";
  fileName[i]+=type.Data();
  // access this file:
  dirFile[i] = (TDirectoryFile*)mergedFile->FindObjectAny(fileName[i].Data());
  // form a list name for each method:
  listName[i]+="cobj";
  listName[i]+=method[i].Data();
  // access this list and close the file:
  if(dirFile[i])
  {
   dirFile[i]->GetObject(listName[i].Data(),list[i]);
   dirFile[i]->Close();
  }
 } 
 
 // close the mergedFile:
 mergedFile->Close();
 
 // create a new file which will hold the final results of all methods:
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");

 // create a new file for each method wich will hold list with final results:
 TDirectoryFile *dirFileFinal[nMethods] = {NULL};
 for(Int_t i=0;i<nMethods;i++)
 {
  if(dirFile[i]) dirFileFinal[i] = new TDirectoryFile(fileName[i].Data(),fileName[i].Data());
 } 

 // redo finish for each method (REMARK: this implementation can be dramatically improved!):
 // MCEP:
 for(Int_t i=0;i<nMethods;i++)
 {
  if(list[i] && strcmp(list[i]->GetName(),"cobjMCEP")==0)
  {
   AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
   mcep->GetOutputHistograms(list[i]);
   mcep->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  } 
  // SP:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjSP")==0)
  {
   AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
   sp->GetOutputHistograms(list[i]);
   sp->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  } 
  // GFC:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjGFC")==0)
  {
   AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
   gfc->GetOutputHistograms(list[i]);
   gfc->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  } 
  // QC:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjQC")==0)
  {
   AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
   qc->GetOutputHistograms(list[i]);
   qc->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  } 
  // FQD:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjFQD")==0)
  {
   AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
   fqd->GetOutputHistograms(list[i]);
   fqd->Finish(kTRUE);
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  }
  // LYZ1SUM:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ1SUM")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
   lyz1sum->SetFirstRun(kTRUE);   
   lyz1sum->SetUseSum(kTRUE);       
   lyz1sum->GetOutputHistograms(list[i]);
   lyz1sum->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  } 
  // LYZ2SUM:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ2SUM")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
   lyz2sum->SetFirstRun(kFALSE);   
   lyz2sum->SetUseSum(kTRUE);       
   lyz2sum->GetOutputHistograms(list[i]);
   lyz2sum->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  }
  // LYZ1PROD:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ1PROD")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
   lyz1prod->SetFirstRun(kTRUE);   
   lyz1prod->SetUseSum(kFALSE);       
   lyz1prod->GetOutputHistograms(list[i]);
   lyz1prod->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  }    
  // LYZ2PROD:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ2PROD")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
   lyz2prod->SetFirstRun(kFALSE);   
   lyz2prod->SetUseSum(kFALSE);       
   lyz2prod->GetOutputHistograms(list[i]);
   lyz2prod->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  }
  // LYZEP:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZEP")==0)
  {
   AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
   lyzep->GetOutputHistograms(list[i]);
   lyzep->Finish();
   dirFileFinal[i]->Add(list[i]);
   dirFileFinal[i]->Write(dirFileFinal[i]->GetName(),TObject::kSingleKey);
  }                
 } // end of for(Int_t i=0;i<nMethods;i++)

 // close the final output file:
 outputFile->Close();
 delete outputFile;
     
} // end of void reCallFinish(Int_t mode=mLocal)
 
void LoadLibrariesRF(const libModes mode) {
  
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
  
} // end of void LoadLibrariesRF(const libModes mode)

 
