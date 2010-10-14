// Macro redoFinish.C is typically used after the merging macros (mergeOutput.C or 
// mergeOutputOnGrid.C) have been used to produce the merged, large statistics
// file of flow analysis. Results stored in merged file are WRONG because after 
// merging the results from small statistics files are trivially summed up in all
// histograms. This is taken into account and corrected for with macro redoFinish.C.
// Another typical use of the macro redoFinish.C is to repeat the call to Finish()
// in all classes, but with different values of some settings which might modify
// the final results (Example: redo the Finish() and apply correction for detector
// effects in QC code because by default this correction is switched off).

// Name of the merged, large statistics file obtained with the merging macros: 
TString mergedFileName = "mergedAnalysisResults.root";
// Final output file name holding correct final results for large statistics sample:
TString outputFileName = "AnalysisResults.root"; 
// Methods for which the flow analysis will be redone:
const Int_t nMethods = 12;
TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
// Methods for which some settings can be modified before redoing Finish():
// QC:
Bool_t bApplyCorrectionForNUA = kFALSE; // apply correction for non-uniform acceptance
Bool_t bApplyCorrectionForNUAVsM = kFALSE; // apply correction for non-uniform acceptance in each multiplicity bin independently
Bool_t bPropagateErrorAlsoFromNIT = kFALSE; // propagate error also from non-isotropic terms
Bool_t bFinalRefFlowResultIsRebinnedInM = kFALSE; // store in CRH for reference flow the result obtained after rebinning in multiplicity

enum libModes {mLocal,mLocalSource};

void redoFinish(TString type="", Int_t mode=mLocal)
{
 // type: type of analysis can be ESD, AOD, MC, ....
 //       (if type="" output files are from simulation 'on the fly')
 // mode: if mode = mLocal: analyze data on your computer using aliroot
 //       if mode = mLocalSource: analyze data on your computer using root + source files

 // Load needed libraries:                       
 LoadLibrariesRF(mode);  

 // Accessing <mergedFileName>:
 TString mergedFileFullPathName(gSystem->pwd());
 mergedFileFullPathName+="/";
 mergedFileFullPathName+=mergedFileName.Data();
 TFile *mergedFile = NULL;
 if(gSystem->AccessPathName(mergedFileFullPathName.Data(),kFileExists))
 {
  cout<<endl;
  cout<<" WARNING: Couldn't find a file: "<<mergedFileName.Data()<<endl;
  cout<<"          in directory "<<gSystem->pwd()<<" !!!!"<<endl;
  cout<<endl;
  exit(0);
 } else 
   {
    // Create temporarily copy of <mergedFileName> if neccessary:
    if(!(mergedFileName == outputFileName))
    {
     TSystemFile *fileTemp = new TSystemFile(mergedFileFullPathName.Data(),".");
     fileTemp->Copy("mergedAnalysisResultsTemp.root");
     delete fileTemp;
    }
    // Access <mergedFileName>:
    mergedFile = TFile::Open(mergedFileFullPathName.Data(),"UPDATE");
   }
   
 // Access from <mergedFileName> the merged TDirectoryFile's for each method and from them the lists holding histograms:
 TString fileName[nMethods]; 
 TDirectoryFile *dirFile[nMethods] = {NULL}; 
 TString listName[nMethods]; 
 TList *list[nMethods] = {NULL};
 Int_t failureCounter = 0;
 cout<<endl;
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
  if(dirFile[i])
  {
   TList* listTemp = dirFile[i]->GetListOfKeys();
   if(listTemp && listTemp->GetEntries() == 1)
   {
    listName[i] = listTemp->At(0)->GetName(); // to be improved - implemented better
    dirFile[i]->GetObject(listName[i].Data(),list[i]);
   } else
     {
      cout<<" WARNING: Accessing TList from TDirectoryFile failed for method "<<method[i].Data()<<" !!!!"<<endl;
      cout<<"          Did you actually specify "<<method[i].Data()<<" in the analysis?"<<endl;
      cout<<endl;
     }
  } else 
    {
     cout<<" WARNING: Couldn't find a TDirectoryFile "<<fileName[i].Data()<<".root !!!!"<<endl;
     failureCounter++;
    }   
 } // End of for(Int_t i=0;i<nMethods;i++)
 
 if(failureCounter == nMethods)
 {       
  cout<<endl;
  cout<<" Did you specify 'TString type' correctly? Can be \"ESD\",\"AOD\",\"MC\",\"\", .... "<<endl;
  cout<<endl;  
 } else if(failureCounter < nMethods && failureCounter > 0)
   {
    cout<<endl; // cosmetics
   }

 // Redo Finish() for each method (to be improved - implemented better):
 // MCEP:
 for(Int_t i=0;i<nMethods;i++)
 {
  if(list[i] && TString(list[i]->GetName()).Contains("cobjMCEP"))
  {
   AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
   mcep->GetOutputHistograms(list[i]);
   mcep->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // SP:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjSP"))
  {
   AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
   sp->GetOutputHistograms(list[i]);
   sp->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // GFC:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjGFC"))
  {
   AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants();
   gfc->GetOutputHistograms(list[i]);
   gfc->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // QC:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjQC"))
  {
   AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants();
   qc->GetOutputHistograms(list[i]);
   if(bApplyCorrectionForNUA){qc->GetIntFlowFlags()->SetBinContent(3,1);} 
   if(bApplyCorrectionForNUAVsM){qc->GetIntFlowFlags()->SetBinContent(8,1);} 
   if(bPropagateErrorAlsoFromNIT){qc->GetIntFlowFlags()->SetBinContent(9,1);}
   if(bFinalRefFlowResultIsRebinnedInM){qc->GetIntFlowFlags()->SetBinContent(11,0);}
   qc->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // FQD:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjFQD"))
  {
   AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
   fqd->GetOutputHistograms(list[i]);
   fqd->Finish(kTRUE);
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }
  // LYZ1SUM:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjLYZ1SUM"))
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
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjLYZ2SUM"))
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
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjLYZ1PROD"))
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
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjLYZ2PROD"))
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
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjLYZEP"))
  {
   AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
   lyzep->GetOutputHistograms(list[i]);
   lyzep->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // MH:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjMH"))
  {
   AliFlowAnalysisWithMixedHarmonics* mh = new AliFlowAnalysisWithMixedHarmonics();
   mh->GetOutputHistograms(list[i]);
   mh->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  } 
  // NL:
  else if(list[i] && TString(list[i]->GetName()).Contains("cobjNL"))
  {
   AliFlowAnalysisWithNestedLoops* nl = new AliFlowAnalysisWithNestedLoops();
   nl->GetOutputHistograms(list[i]);
   nl->Finish();
   dirFile[i]->Add(list[i],kTRUE);
   dirFile[i]->Write(dirFile[i]->GetName(),TObject::kSingleKey+TObject::kOverwrite);
  }                
 } // end of for(Int_t i=0;i<nMethods;i++)

 // Close the final output file:
 delete mergedFile;
 
 // Giving the final names if neccessary:
 if(!(mergedFileName == outputFileName))
 {
  TSystemFile *outputFileFinal = new TSystemFile(mergedFileName.Data(),".");
  outputFileFinal->Rename(outputFileName.Data());
  delete outputFileFinal;
  TSystemFile *mergedFileFinal = new TSystemFile("mergedAnalysisResultsTemp.root",".");
  mergedFileFinal->Rename(mergedFileName.Data());
  delete mergedFileFinal;
 } // end of if(!(mergedFileName == outputFileName))
    
 cout<<endl;
    
} // end of void redoFinish(Int_t mode=mLocal)
 
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
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
        
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

 
