// Macro redoFinish.C is used after macro mergeOutput.C has been used for merging. 
// Before using this macro first read the explanation at the beginning of macro mergeOutput.C. 

enum libModes {mLocal,mLocalSource};

void redoFinish(TString type="ESD", Int_t mode=mLocal)
{
 // type: type of analysis can be ESD, AOD, MC, ESDMCkineESD, ESDMCkineMC
 //       (if type="" output files are from MC simulation (default))
 // mode: if mode = mLocal: analyze data on your computer using aliroot
 //       if mode = mLocalSource: analyze data on your computer using root + source files

 // Name of merged, large statistics file obtained with macro mergeOutput.C: 
 TString mergedFileName = "mergedAnalysisResults.root";
 // Final output file name holding final results for large statistics sample:
 TString outputFileName = "AnalysisResults.root"; 
 
 const Int_t nMethods = 12;
 TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
 
 // Load needed libraries:                       
 LoadLibrariesRF(mode);  
   
 // Accessing the merged, large statistics file obtained with macro mergeOutput.C.
 // On this file the flow analysis will be redone, its content modified to store
 // the correct final results and eventually renamed into TString outputFileName.
 TString pwd(gSystem->pwd());
 pwd+="/";
 pwd+=mergedFileName.Data();
 TFile *mergedFile = NULL;
 if(gSystem->AccessPathName(pwd.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file:"<<endl;
  cout<<"         "<<pwd.Data()<<endl;
  cout<<endl;
  cout<<"In order to get that file use macro mergeOutput.C first."<<endl;
  cout<<endl;
  exit(0);
 } else 
   {
    // Create temporarily copy of "mergedAnalysisResults.root":
    TSystemFile *fileTemp = new TSystemFile(mergedFileName.Data(),".");
    fileTemp->Copy("mergedAnalysisResultsTemp.root");
    delete fileTemp;
    // Access merged file:
    mergedFile = TFile::Open(pwd.Data(),"UPDATE");
   }
   
 // Access from mergedFile the merged files for each method and from them the lists holding histograms:
 TString fileName[nMethods]; 
 TDirectoryFile *dirFile[nMethods] = {NULL}; 
 TString listName[nMethods]; 
 TList *list[nMethods] = {NULL};
 for(Int_t i=0;i<nMethods;i++)
 {
  // Form a file name for each method:
  fileName[i]+="output";
  fileName[i]+=method[i].Data();
  fileName[i]+="analysis";
  fileName[i]+=type.Data();
  // Access this file:
  dirFile[i] = (TDirectoryFile*)mergedFile->FindObjectAny(fileName[i].Data());
  // Form a list name for each method:
  listName[i]+="cobj";
  listName[i]+=method[i].Data();
  // Access this list:
  if(dirFile[i])
  {
   dirFile[i]->GetObject(listName[i].Data(),list[i]);
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileName[i].Data()<<".root !!!!"<<endl;
    }
 } // End of for(Int_t i=0;i<nMethods;i++)

 // Redo finish for each method (REMARK: this implementation can be dramatically improved!):
 // MCEP:
 for(Int_t i=0;i<nMethods;i++)
 {
  if(list[i] && strcmp(list[i]->GetName(),"cobjMCEP")==0)
  {
   AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
   mcep->GetOutputHistograms(list[i]);
   mcep->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // SP:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjSP")==0)
  {
   AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
   sp->GetOutputHistograms(list[i]);
   sp->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // GFC:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjGFC")==0)
  {
   AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
   gfc->GetOutputHistograms(list[i]);
   gfc->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // QC:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjQC")==0)
  {
   AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
   qc->GetOutputHistograms(list[i]);
   qc->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // FQD:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjFQD")==0)
  {
   AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
   fqd->GetOutputHistograms(list[i]);
   fqd->Finish(kTRUE);
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }
  // LYZ1SUM:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ1SUM")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
   lyz1sum->SetFirstRun(kTRUE);   
   lyz1sum->SetUseSum(kTRUE);       
   lyz1sum->GetOutputHistograms(list[i]);
   lyz1sum->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // LYZ2SUM:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ2SUM")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
   lyz2sum->SetFirstRun(kFALSE);   
   lyz2sum->SetUseSum(kTRUE);       
   lyz2sum->GetOutputHistograms(list[i]);
   lyz2sum->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }
  // LYZ1PROD:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ1PROD")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
   lyz1prod->SetFirstRun(kTRUE);   
   lyz1prod->SetUseSum(kFALSE);       
   lyz1prod->GetOutputHistograms(list[i]);
   lyz1prod->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }    
  // LYZ2PROD:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZ2PROD")==0)
  {
   AliFlowAnalysisWithLeeYangZeros* lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
   lyz2prod->SetFirstRun(kFALSE);   
   lyz2prod->SetUseSum(kFALSE);       
   lyz2prod->GetOutputHistograms(list[i]);
   lyz2prod->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }
  // LYZEP:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjLYZEP")==0)
  {
   AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
   lyzep->GetOutputHistograms(list[i]);
   lyzep->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // MH:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjMH")==0)
  {
   AliFlowAnalysisWithMixedHarmonics* mh = new AliFlowAnalysisWithMixedHarmonics();
   mh->GetOutputHistograms(list[i]);
   mh->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // NL:
  else if(list[i] && strcmp(list[i]->GetName(),"cobjNL")==0)
  {
   AliFlowAnalysisWithNestedLoops* nl = new AliFlowAnalysisWithNestedLoops();
   nl->GetOutputHistograms(list[i]);
   nl->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }                
 } // End of for(Int_t i=0;i<nMethods;i++)

 // Close the final output file:
 mergedFile->Close();
 delete mergedFile;
 
 // Giving the final names:
 TSystemFile *outputFileFinal = new TSystemFile(mergedFileName.Data(),".");
 outputFileFinal->Rename(outputFileName.Data());
 delete outputFileFinal; 
 TSystemFile *mergedFileFinal = new TSystemFile("mergedAnalysisResultsTemp.root",".");
 mergedFileFinal->Rename(mergedFileName.Data());
 delete mergedFileFinal;
     
} // End of void redoFinish(Int_t mode=mLocal)
 
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
  
} // end of void LoadLibrariesRF(const libModes mode)

 
