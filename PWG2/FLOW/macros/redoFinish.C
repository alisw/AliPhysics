enum libModes {mLocal,mLocalSource};

void redoFinish(TString type="", Int_t mode=mLocal)
{
 // type:  type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //        (if type="" output files are from MC simulation (default))
 // mode:  if mode=mLocal analyze data on your computer using aliroot
 //        if mode=mLocalSource analyze data on your computer using root + source files
 
 Bool_t redoFinishAlsoInSubsets = kFALSE;
 
 // load needed libraries:                       
 LoadLibrariesRDF(mode);  
 
 // access the path of current diretory:
 TString pwd(gSystem->pwd());
 pwd+="/";
 
 cout<<endl;
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
  delete mcep;
  delete finalOutputMCEP;
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
  delete sp;
  delete finalOutputSP;
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
  delete gfc;
  delete finalOutputGFC;
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
  delete qc;
  delete finalOutputQC;
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
  AliFlowAnalysisWithFittingQDistribution* fqd = new AliFlowAnalysisWithFittingQDistribution();
  fqd->GetOutputHistograms(mergedListFQD);
  fqd->Finish(kTRUE);
  // save the final results for FQD in final output file: 
  TString finalOutputFileNameFQD("outputFQDanalysis");
  (finalOutputFileNameFQD+=(type.Data()))+=(".root");
  TString pwdFinalFQD=pwd.Data();
  pwdFinalFQD+=finalOutputFileNameFQD;
  TFile *finalOutputFQD = new TFile(pwdFinalFQD.Data(),"NEW");
  mergedListFQD->SetName("cobjFQD");
  mergedListFQD->Write(mergedListFQD->GetName(),TObject::kSingleKey);
  finalOutputFQD->Close();
  delete fqd;
  delete finalOutputFQD;
 } else 
   {
    cout<<"WARNING: mergedListFQD is NULL !!!!"<<endl;
   }                               

 // LYZ1SUM:
 TString mergedFileNameLYZ1SUM("mergedLYZ1SUManalysis");
 (mergedFileNameLYZ1SUM+=(type.Data()))+=(".root");
 TFile *mergedFileLYZ1SUM = NULL;
 TList *mergedListLYZ1SUM = NULL;
 TString pwdLYZ1SUM=pwd.Data();
 pwdLYZ1SUM+=mergedFileNameLYZ1SUM;
 if(gSystem->AccessPathName(pwdLYZ1SUM.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdLYZ1SUM.Data()<<endl;
 } else 
   {
    mergedFileLYZ1SUM = TFile::Open(pwdLYZ1SUM.Data(),"READ");
    if(mergedFileLYZ1SUM) 
    {
     mergedFileLYZ1SUM->GetObject("cobjLYZ1SUM",mergedListLYZ1SUM);
     mergedFileLYZ1SUM->Close();
    } 
   }
 if(mergedListLYZ1SUM)
 {
  AliFlowAnalysisWithLeeYangZeros* lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
  lyz1sum->SetFirstRun(kTRUE);   
  lyz1sum->SetUseSum(kTRUE);       
  lyz1sum->GetOutputHistograms(mergedListLYZ1SUM);
  lyz1sum->Finish();
  // save the final results for LYZ1SUM in final output file: 
  TString finalOutputFileNameLYZ1SUM("outputLYZ1SUManalysis");
  (finalOutputFileNameLYZ1SUM+=(type.Data()))+=(".root");
  TString pwdFinalLYZ1SUM=pwd.Data();
  pwdFinalLYZ1SUM+=finalOutputFileNameLYZ1SUM;
  TFile *finalOutputLYZ1SUM = new TFile(pwdFinalLYZ1SUM.Data(),"NEW");
  mergedListLYZ1SUM->SetName("cobjLYZ1SUM");
  mergedListLYZ1SUM->Write(mergedListLYZ1SUM->GetName(),TObject::kSingleKey);
  finalOutputLYZ1SUM->Close();
  delete lyz1sum;
  delete finalOutputLYZ1SUM;
 } else 
   {
    cout<<"WARNING: mergedListLYZ1SUM is NULL !!!!"<<endl;
   }  
   
 // LYZ2SUM:
 TString mergedFileNameLYZ2SUM("mergedLYZ2SUManalysis");
 (mergedFileNameLYZ2SUM+=(type.Data()))+=(".root");
 TFile *mergedFileLYZ2SUM = NULL;
 TList *mergedListLYZ2SUM = NULL;
 TString pwdLYZ2SUM=pwd.Data();
 pwdLYZ2SUM+=mergedFileNameLYZ2SUM;
 if(gSystem->AccessPathName(pwdLYZ2SUM.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdLYZ2SUM.Data()<<endl;
 } else 
   {
    mergedFileLYZ2SUM = TFile::Open(pwdLYZ2SUM.Data(),"READ");
    if(mergedFileLYZ2SUM) 
    {
     mergedFileLYZ2SUM->GetObject("cobjLYZ2SUM",mergedListLYZ2SUM);
     mergedFileLYZ2SUM->Close();
    } 
   }
 if(mergedListLYZ2SUM)
 {
  AliFlowAnalysisWithLeeYangZeros* lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
  lyz2sum->SetFirstRun(kFALSE);   
  lyz2sum->SetUseSum(kTRUE);       
  lyz2sum->GetOutputHistograms(mergedListLYZ2SUM);
  lyz2sum->Finish();
  // save the final results for LYZ2SUM in final output file: 
  TString finalOutputFileNameLYZ2SUM("outputLYZ2SUManalysis");
  (finalOutputFileNameLYZ2SUM+=(type.Data()))+=(".root");
  TString pwdFinalLYZ2SUM=pwd.Data();
  pwdFinalLYZ2SUM+=finalOutputFileNameLYZ2SUM;
  TFile *finalOutputLYZ2SUM = new TFile(pwdFinalLYZ2SUM.Data(),"NEW");
  mergedListLYZ2SUM->SetName("cobjLYZ2SUM");
  mergedListLYZ2SUM->Write(mergedListLYZ2SUM->GetName(),TObject::kSingleKey);
  finalOutputLYZ2SUM->Close();
  delete lyz2sum;
  delete finalOutputLYZ2SUM ;
 } else 
   {
    cout<<"WARNING: mergedListLYZ2SUM is NULL !!!!"<<endl;
   }                                                                                           
 
 // LYZ1PROD:
 TString mergedFileNameLYZ1PROD("mergedLYZ1PRODanalysis");
 (mergedFileNameLYZ1PROD+=(type.Data()))+=(".root");
 TFile *mergedFileLYZ1PROD = NULL;
 TList *mergedListLYZ1PROD = NULL;
 TString pwdLYZ1PROD=pwd.Data();
 pwdLYZ1PROD+=mergedFileNameLYZ1PROD;
 if(gSystem->AccessPathName(pwdLYZ1PROD.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdLYZ1PROD.Data()<<endl;
 } else 
   {
    mergedFileLYZ1PROD = TFile::Open(pwdLYZ1PROD.Data(),"READ");
    if(mergedFileLYZ1PROD) 
    {
     mergedFileLYZ1PROD->GetObject("cobjLYZ1PROD",mergedListLYZ1PROD);
     mergedFileLYZ1PROD->Close();
    } 
   }
 if(mergedListLYZ1PROD)
 {
  AliFlowAnalysisWithLeeYangZeros* lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
  lyz1prod->SetFirstRun(kTRUE);   
  lyz1prod->SetUseSum(kFALSE);       
  lyz1prod->GetOutputHistograms(mergedListLYZ1PROD);
  lyz1prod->Finish();
  // save the final results for LYZ1PROD in final output file: 
  TString finalOutputFileNameLYZ1PROD("outputLYZ1PRODanalysis");
  (finalOutputFileNameLYZ1PROD+=(type.Data()))+=(".root");
  TString pwdFinalLYZ1PROD=pwd.Data();
  pwdFinalLYZ1PROD+=finalOutputFileNameLYZ1PROD;
  TFile *finalOutputLYZ1PROD = new TFile(pwdFinalLYZ1PROD.Data(),"NEW");
  mergedListLYZ1PROD->SetName("cobjLYZ1PROD");
  mergedListLYZ1PROD->Write(mergedListLYZ1PROD->GetName(),TObject::kSingleKey);
  finalOutputLYZ1PROD->Close();
  delete lyz1prod;
  delete finalOutputLYZ1PROD;
 } else 
   {
    cout<<"WARNING: mergedListLYZ1PROD is NULL !!!!"<<endl;
   }  
   
 // LYZ2PROD:
 TString mergedFileNameLYZ2PROD("mergedLYZ2PRODanalysis");
 (mergedFileNameLYZ2PROD+=(type.Data()))+=(".root");
 TFile *mergedFileLYZ2PROD = NULL;
 TList *mergedListLYZ2PROD = NULL;
 TString pwdLYZ2PROD=pwd.Data();
 pwdLYZ2PROD+=mergedFileNameLYZ2PROD;
 if(gSystem->AccessPathName(pwdLYZ2PROD.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdLYZ2PROD.Data()<<endl;
 } else 
   {
    mergedFileLYZ2PROD = TFile::Open(pwdLYZ2PROD.Data(),"READ");
    if(mergedFileLYZ2PROD) 
    {
     mergedFileLYZ2PROD->GetObject("cobjLYZ2PROD",mergedListLYZ2PROD);
     mergedFileLYZ2PROD->Close();
    } 
   }
 if(mergedListLYZ2PROD)
 {
  AliFlowAnalysisWithLeeYangZeros* lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
  lyz2prod->SetFirstRun(kFALSE);   
  lyz2prod->SetUseSum(kFALSE);       
  lyz2prod->GetOutputHistograms(mergedListLYZ2PROD);
  lyz2prod->Finish();
  // save the final results for LYZ2PROD in final output file: 
  TString finalOutputFileNameLYZ2PROD("outputLYZ2PRODanalysis");
  (finalOutputFileNameLYZ2PROD+=(type.Data()))+=(".root");
  TString pwdFinalLYZ2PROD=pwd.Data();
  pwdFinalLYZ2PROD+=finalOutputFileNameLYZ2PROD;
  TFile *finalOutputLYZ2PROD = new TFile(pwdFinalLYZ2PROD.Data(),"NEW");
  mergedListLYZ2PROD->SetName("cobjLYZ2PROD");
  mergedListLYZ2PROD->Write(mergedListLYZ2PROD->GetName(),TObject::kSingleKey);
  finalOutputLYZ2PROD->Close();
  delete lyz2prod;
  delete finalOutputLYZ2PROD;
 } else 
   {
    cout<<"WARNING: mergedListLYZ2PROD is NULL !!!!"<<endl;
   } 
   
 // LYZEP:
 TString mergedFileNameLYZEP("mergedLYZEPanalysis");
 (mergedFileNameLYZEP+=(type.Data()))+=(".root");
 TFile *mergedFileLYZEP = NULL;
 TList *mergedListLYZEP = NULL;
 TString pwdLYZEP=pwd.Data();
 pwdLYZEP+=mergedFileNameLYZEP;
 if(gSystem->AccessPathName(pwdLYZEP.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<pwdLYZEP.Data()<<endl;
 } else 
   {
    mergedFileLYZEP = TFile::Open(pwdLYZEP.Data(),"READ");
    if(mergedFileLYZEP) 
    {
     mergedFileLYZEP->GetObject("cobjLYZEP",mergedListLYZEP);
     mergedFileLYZEP->Close();
    } 
   }
 if(mergedListLYZEP)
 {
  AliFlowAnalysisWithLYZEventPlane* lyzep = new AliFlowAnalysisWithLYZEventPlane();
  lyzep->GetOutputHistograms(mergedListLYZEP);
  lyzep->Finish();
  // save the final results for LYZEP in final output file: 
  TString finalOutputFileNameLYZEP("outputLYZEPanalysis");
  (finalOutputFileNameLYZEP+=(type.Data()))+=(".root");
  TString pwdFinalLYZEP=pwd.Data();
  pwdFinalLYZEP+=finalOutputFileNameLYZEP;
  TFile *finalOutputLYZEP = new TFile(pwdFinalLYZEP.Data(),"NEW");
  mergedListLYZEP->SetName("cobjLYZEP");
  mergedListLYZEP->Write(mergedListLYZEP->GetName(),TObject::kSingleKey);
  finalOutputLYZEP->Close();
  delete lyzep;
  delete finalOutputLYZEP;
 } else 
   {
    cout<<"WARNING: mergedListLYZEP is NULL !!!!"<<endl;
   }                                                                                            
 
 // redoing Finish() also in subsets:
 if(redoFinishAlsoInSubsets)
 {
  // standard magic:
  TString *baseDirPath = new TString(gSystem->pwd());
  TSystemDirectory* baseDir = new TSystemDirectory(".",baseDirPath->Data());          
  TList* listOfFilesInBaseDir = baseDir->GetListOfFiles();
  // listOfFilesInBaseDir->Print();
  Int_t nFiles = listOfFilesInBaseDir->GetEntries();
  Int_t nSubsets = 0;
  gSystem->cd(baseDirPath->Data());
  // count subsets:
  for(Int_t iFile=0;iFile<nFiles;iFile++)
  {
   TSystemFile* presentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
   TString presentFileName = baseDirPath->Data();
   (presentFileName+="/")+=presentFile->GetName(); 
   if(presentFileName.Contains("subset"))
   {
    nSubsets++;
   } 
  } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)
  
  // redo finish in each subset:
  gSystem->cd(baseDirPath->Data());
  for(Int_t i=0;i<nSubsets;i++)
  {
   cout<<endl;
   // access the merged output files and redo the flow analysis on them:
   // MCEP:
   TString *mergedFileNameInSubsetMCEP = new TString("subset"); 
   (*mergedFileNameInSubsetMCEP)+=(i+1);
   (*mergedFileNameInSubsetMCEP)+="/";
   (*mergedFileNameInSubsetMCEP)+="mergedMCEPanalysis";
   (*mergedFileNameInSubsetMCEP)+=type.Data();
   (*mergedFileNameInSubsetMCEP)+=".root";
   
   TFile *mergedFileInSubsetMCEP = NULL;
   TList *mergedListInSubsetMCEP = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetMCEP->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetMCEP->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetMCEP = TFile::Open(mergedFileNameInSubsetMCEP->Data(),"READ");
      delete mergedFileNameInSubsetMCEP;
      if(mergedFileinSubsetMCEP) 
      {
       mergedFileinSubsetMCEP->GetObject("cobjMCEP",mergedListInSubsetMCEP);
       mergedFileinSubsetMCEP->Close();
      } 
     }
     if(mergedListInSubsetMCEP)
     {
      AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
      mcep->GetOutputHistograms(mergedListInSubsetMCEP);
      mcep->Finish();
      // save the final results for MCEP in final output file: 
      TString *finalOutputFileNameinSubsetMCEP = new TString("subset");
      (*finalOutputFileNameinSubsetMCEP)+=(i+1);  
      (*finalOutputFileNameinSubsetMCEP)+="/";
      (*finalOutputFileNameinSubsetMCEP)+="outputMCEPanalysis";
      (*finalOutputFileNameinSubsetMCEP)+=type.Data();
      (*finalOutputFileNameinSubsetMCEP)+=".root";
      TFile *finalOutputInSubsetMCEP = new TFile(finalOutputFileNameinSubsetMCEP->Data(),"NEW");
      mergedListInSubsetMCEP->SetName("cobjMCEP");
      mergedListInSubsetMCEP->Write(mergedListInSubsetMCEP->GetName(),TObject::kSingleKey);
      finalOutputInSubsetMCEP->Close();
      delete finalOutputFileNameinSubsetMCEP;
      delete finalOutputInSubsetMCEP;
      delete mcep;
     } else 
       {
        cout<<"WARNING: mergedListMCEP is NULL !!!!"<<endl;
       }
           
   // SP:
   TString *mergedFileNameInSubsetSP = new TString("subset"); 
   (*mergedFileNameInSubsetSP)+=(i+1);
   (*mergedFileNameInSubsetSP)+="/";
   (*mergedFileNameInSubsetSP)+="mergedSPanalysis";
   (*mergedFileNameInSubsetSP)+=type.Data();
   (*mergedFileNameInSubsetSP)+=".root";
   
   TFile *mergedFileInSubsetSP = NULL;
   TList *mergedListInSubsetSP = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetSP->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetSP->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetSP = TFile::Open(mergedFileNameInSubsetSP->Data(),"READ");
      delete mergedFileNameInSubsetSP;
      if(mergedFileinSubsetSP) 
      {
       mergedFileinSubsetSP->GetObject("cobjSP",mergedListInSubsetSP);
       mergedFileinSubsetSP->Close();
      } 
     }
     if(mergedListInSubsetSP)
     {
      AliFlowAnalysisWithScalarProduct *sp = new AliFlowAnalysisWithScalarProduct();
      sp->GetOutputHistograms(mergedListInSubsetSP);
      sp->Finish();
      // save the final results for SP in final output file: 
      TString *finalOutputFileNameinSubsetSP = new TString("subset");
      (*finalOutputFileNameinSubsetSP)+=(i+1);  
      (*finalOutputFileNameinSubsetSP)+="/";
      (*finalOutputFileNameinSubsetSP)+="outputSPanalysis";
      (*finalOutputFileNameinSubsetSP)+=type.Data();
      (*finalOutputFileNameinSubsetSP)+=".root";
      TFile *finalOutputInSubsetSP = new TFile(finalOutputFileNameinSubsetSP->Data(),"NEW");
      mergedListInSubsetSP->SetName("cobjSP");
      mergedListInSubsetSP->Write(mergedListInSubsetSP->GetName(),TObject::kSingleKey);
      finalOutputInSubsetSP->Close();
      delete finalOutputFileNameinSubsetSP;
      delete finalOutputInSubsetSP;
      delete sp;
     } else 
       {
        cout<<"WARNING: mergedListSP is NULL !!!!"<<endl;
       }    
      
   // GFC:
   TString *mergedFileNameInSubsetGFC = new TString("subset"); 
   (*mergedFileNameInSubsetGFC)+=(i+1);
   (*mergedFileNameInSubsetGFC)+="/";
   (*mergedFileNameInSubsetGFC)+="mergedGFCanalysis";
   (*mergedFileNameInSubsetGFC)+=type.Data();
   (*mergedFileNameInSubsetGFC)+=".root";
   
   TFile *mergedFileInSubsetGFC = NULL;
   TList *mergedListInSubsetGFC = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetGFC->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetGFC->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetGFC = TFile::Open(mergedFileNameInSubsetGFC->Data(),"READ");
      delete mergedFileNameInSubsetGFC;
      if(mergedFileinSubsetGFC) 
      {
       mergedFileinSubsetGFC->GetObject("cobjGFC",mergedListInSubsetGFC);
       mergedFileinSubsetGFC->Close();
      } 
     }
     if(mergedListInSubsetGFC)
     {
      AliFlowAnalysisWithCumulants *gfc = new AliFlowAnalysisWithCumulants();
      gfc->GetOutputHistograms(mergedListInSubsetGFC);
      gfc->Finish();
      // save the final results for GFC in final output file: 
      TString *finalOutputFileNameinSubsetGFC = new TString("subset");
      (*finalOutputFileNameinSubsetGFC)+=(i+1);  
      (*finalOutputFileNameinSubsetGFC)+="/";
      (*finalOutputFileNameinSubsetGFC)+="outputGFCanalysis";
      (*finalOutputFileNameinSubsetGFC)+=type.Data();
      (*finalOutputFileNameinSubsetGFC)+=".root";
      TFile *finalOutputInSubsetGFC = new TFile(finalOutputFileNameinSubsetGFC->Data(),"NEW");
      mergedListInSubsetGFC->SetName("cobjGFC");
      mergedListInSubsetGFC->Write(mergedListInSubsetGFC->GetName(),TObject::kSingleKey);
      finalOutputInSubsetGFC->Close();
      delete finalOutputFileNameinSubsetGFC;
      delete finalOutputInSubsetGFC;
      delete gfc;
     } else 
       {
        cout<<"WARNING: mergedListGFC is NULL !!!!"<<endl;
       }        
   
   // QC:
   TString *mergedFileNameInSubsetQC = new TString("subset"); 
   (*mergedFileNameInSubsetQC)+=(i+1);
   (*mergedFileNameInSubsetQC)+="/";
   (*mergedFileNameInSubsetQC)+="mergedQCanalysis";
   (*mergedFileNameInSubsetQC)+=type.Data();
   (*mergedFileNameInSubsetQC)+=".root";
   
   TFile *mergedFileInSubsetQC = NULL;
   TList *mergedListInSubsetQC = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetQC->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetQC->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetQC = TFile::Open(mergedFileNameInSubsetQC->Data(),"READ");
      delete mergedFileNameInSubsetQC;
      if(mergedFileinSubsetQC) 
      {
       mergedFileinSubsetQC->GetObject("cobjQC",mergedListInSubsetQC);
       mergedFileinSubsetQC->Close();
      } 
     }
     if(mergedListInSubsetQC)
     {
      AliFlowAnalysisWithQCumulants *qc = new AliFlowAnalysisWithQCumulants();
      qc->GetOutputHistograms(mergedListInSubsetQC);
      qc->Finish();
      // save the final results for QC in final output file: 
      TString *finalOutputFileNameinSubsetQC = new TString("subset");
      (*finalOutputFileNameinSubsetQC)+=(i+1);  
      (*finalOutputFileNameinSubsetQC)+="/";
      (*finalOutputFileNameinSubsetQC)+="outputQCanalysis";
      (*finalOutputFileNameinSubsetQC)+=type.Data();
      (*finalOutputFileNameinSubsetQC)+=".root";
      TFile *finalOutputInSubsetQC = new TFile(finalOutputFileNameinSubsetQC->Data(),"NEW");
      mergedListInSubsetQC->SetName("cobjQC");
      mergedListInSubsetQC->Write(mergedListInSubsetQC->GetName(),TObject::kSingleKey);
      finalOutputInSubsetQC->Close();
      delete finalOutputFileNameinSubsetQC;
      delete finalOutputInSubsetQC;
      delete qc;
     } else 
       {
        cout<<"WARNING: mergedListQC is NULL !!!!"<<endl;
       }        
       
   // FQD:
   TString *mergedFileNameInSubsetFQD = new TString("subset"); 
   (*mergedFileNameInSubsetFQD)+=(i+1);
   (*mergedFileNameInSubsetFQD)+="/";
   (*mergedFileNameInSubsetFQD)+="mergedFQDanalysis";
   (*mergedFileNameInSubsetFQD)+=type.Data();
   (*mergedFileNameInSubsetFQD)+=".root";
   
   TFile *mergedFileInSubsetFQD = NULL;
   TList *mergedListInSubsetFQD = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetFQD->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetFQD->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetFQD = TFile::Open(mergedFileNameInSubsetFQD->Data(),"READ");
      if(mergedFileinSubsetFQD) 
      {
       mergedFileinSubsetFQD->GetObject("cobjFQD",mergedListInSubsetFQD);
       mergedFileinSubsetFQD->Close();
      } 
     }
     if(mergedListInSubsetFQD)
     {
      AliFlowAnalysisWithFittingQDistribution *fqd = new AliFlowAnalysisWithFittingQDistribution();
      fqd->GetOutputHistograms(mergedListInSubsetFQD);
      cout<<endl;
      cout<<"!!!! WARNING WARNING WARNING WARNING !!!!"<<endl;
      cout<<endl;
      cout<<"     If Minuit crashed here, remove the file "<<mergedFileNameInSubsetFQD->Data()<<endl;
      cout<<"     and lunch redoFinish.C again. Good luck! "<<endl; 
      cout<<endl;
      cout<<"!!!! WARNING WARNING WARNING WARNING !!!!"<<endl;
      cout<<endl;
      delete mergedFileNameInSubsetFQD;
      fqd->Finish(kTRUE);
      // save the final results for FQD in final output file: 
      TString *finalOutputFileNameinSubsetFQD = new TString("subset");
      (*finalOutputFileNameinSubsetFQD)+=(i+1);  
      (*finalOutputFileNameinSubsetFQD)+="/";
      (*finalOutputFileNameinSubsetFQD)+="outputFQDanalysis";
      (*finalOutputFileNameinSubsetFQD)+=type.Data();
      (*finalOutputFileNameinSubsetFQD)+=".root";
      TFile *finalOutputInSubsetFQD = new TFile(finalOutputFileNameinSubsetFQD->Data(),"NEW");
      mergedListInSubsetFQD->SetName("cobjFQD");
      mergedListInSubsetFQD->Write(mergedListInSubsetFQD->GetName(),TObject::kSingleKey);
      finalOutputInSubsetFQD->Close();
      delete finalOutputFileNameinSubsetFQD;
      delete finalOutputInSubsetFQD;
      delete fqd;
     } else 
       {
        cout<<"WARNING: mergedListFQD is NULL !!!!"<<endl;
       }            
      
   // LYZ1SUM:
   TString *mergedFileNameInSubsetLYZ1SUM = new TString("subset"); 
   (*mergedFileNameInSubsetLYZ1SUM)+=(i+1);
   (*mergedFileNameInSubsetLYZ1SUM)+="/";
   (*mergedFileNameInSubsetLYZ1SUM)+="mergedLYZ1SUManalysis";
   (*mergedFileNameInSubsetLYZ1SUM)+=type.Data();
   (*mergedFileNameInSubsetLYZ1SUM)+=".root";
   
   TFile *mergedFileInSubsetLYZ1SUM = NULL;
   TList *mergedListInSubsetLYZ1SUM = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetLYZ1SUM->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetLYZ1SUM->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetLYZ1SUM = TFile::Open(mergedFileNameInSubsetLYZ1SUM->Data(),"READ");
      delete mergedFileNameInSubsetLYZ1SUM;
      if(mergedFileinSubsetLYZ1SUM) 
      {
       mergedFileinSubsetLYZ1SUM->GetObject("cobjLYZ1SUM",mergedListInSubsetLYZ1SUM);
       mergedFileinSubsetLYZ1SUM->Close();
      } 
     }
     if(mergedListInSubsetLYZ1SUM)
     {
      AliFlowAnalysisWithLeeYangZeros *lyz1sum = new AliFlowAnalysisWithLeeYangZeros();
      lyz1sum->SetFirstRun(kTRUE);   
      lyz1sum->SetUseSum(kTRUE);     
      lyz1sum->GetOutputHistograms(mergedListInSubsetLYZ1SUM);
      lyz1sum->Finish();
      // save the final results for LYZ1SUM in final output file: 
      TString *finalOutputFileNameinSubsetLYZ1SUM = new TString("subset");
      (*finalOutputFileNameinSubsetLYZ1SUM)+=(i+1);  
      (*finalOutputFileNameinSubsetLYZ1SUM)+="/";
      (*finalOutputFileNameinSubsetLYZ1SUM)+="outputLYZ1SUManalysis";
      (*finalOutputFileNameinSubsetLYZ1SUM)+=type.Data();
      (*finalOutputFileNameinSubsetLYZ1SUM)+=".root";
      TFile *finalOutputInSubsetLYZ1SUM = new TFile(finalOutputFileNameinSubsetLYZ1SUM->Data(),"NEW");
      mergedListInSubsetLYZ1SUM->SetName("cobjLYZ1SUM");
      mergedListInSubsetLYZ1SUM->Write(mergedListInSubsetLYZ1SUM->GetName(),TObject::kSingleKey);
      finalOutputInSubsetLYZ1SUM->Close();
      delete finalOutputFileNameinSubsetLYZ1SUM;
      delete finalOutputInSubsetLYZ1SUM;
      delete lyz1sum;
     } else 
       {
        cout<<"WARNING: mergedListLYZ1SUM is NULL !!!!"<<endl;
       }                    
                                                     
   // LYZ1PROD:
   TString *mergedFileNameInSubsetLYZ1PROD = new TString("subset"); 
   (*mergedFileNameInSubsetLYZ1PROD)+=(i+1);
   (*mergedFileNameInSubsetLYZ1PROD)+="/";
   (*mergedFileNameInSubsetLYZ1PROD)+="mergedLYZ1PRODanalysis";
   (*mergedFileNameInSubsetLYZ1PROD)+=type.Data();
   (*mergedFileNameInSubsetLYZ1PROD)+=".root";
   
   TFile *mergedFileInSubsetLYZ1PROD = NULL;
   TList *mergedListInSubsetLYZ1PROD = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetLYZ1PROD->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetLYZ1PROD->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetLYZ1PROD = TFile::Open(mergedFileNameInSubsetLYZ1PROD->Data(),"READ");
      delete mergedFileNameInSubsetLYZ1PROD;
      if(mergedFileinSubsetLYZ1PROD) 
      {
       mergedFileinSubsetLYZ1PROD->GetObject("cobjLYZ1PROD",mergedListInSubsetLYZ1PROD);
       mergedFileinSubsetLYZ1PROD->Close();
      } 
     }
     if(mergedListInSubsetLYZ1PROD)
     {
      AliFlowAnalysisWithLeeYangZeros *lyz1prod = new AliFlowAnalysisWithLeeYangZeros();
      lyz1prod->SetFirstRun(kTRUE);   
      lyz1prod->SetUseSum(kFALSE);     
      lyz1prod->GetOutputHistograms(mergedListInSubsetLYZ1PROD);
      lyz1prod->Finish();
      // save the final results for LYZ1PROD in final output file: 
      TString *finalOutputFileNameinSubsetLYZ1PROD = new TString("subset");
      (*finalOutputFileNameinSubsetLYZ1PROD)+=(i+1);  
      (*finalOutputFileNameinSubsetLYZ1PROD)+="/";
      (*finalOutputFileNameinSubsetLYZ1PROD)+="outputLYZ1PRODanalysis";
      (*finalOutputFileNameinSubsetLYZ1PROD)+=type.Data();
      (*finalOutputFileNameinSubsetLYZ1PROD)+=".root";
      TFile *finalOutputInSubsetLYZ1PROD = new TFile(finalOutputFileNameinSubsetLYZ1PROD->Data(),"NEW");
      mergedListInSubsetLYZ1PROD->SetName("cobjLYZ1PROD");
      mergedListInSubsetLYZ1PROD->Write(mergedListInSubsetLYZ1PROD->GetName(),TObject::kSingleKey);
      finalOutputInSubsetLYZ1PROD->Close();
      delete finalOutputFileNameinSubsetLYZ1PROD;
      delete finalOutputInSubsetLYZ1PROD;
      delete lyz1prod;
     } else 
       {
        cout<<"WARNING: mergedListLYZ1PROD is NULL !!!!"<<endl;
       }                                                                                                                            
                                                                                                                                                               
   // LYZ2SUM:
   TString *mergedFileNameInSubsetLYZ2SUM = new TString("subset"); 
   (*mergedFileNameInSubsetLYZ2SUM)+=(i+1);
   (*mergedFileNameInSubsetLYZ2SUM)+="/";
   (*mergedFileNameInSubsetLYZ2SUM)+="mergedLYZ2SUManalysis";
   (*mergedFileNameInSubsetLYZ2SUM)+=type.Data();
   (*mergedFileNameInSubsetLYZ2SUM)+=".root";
   
   TFile *mergedFileInSubsetLYZ2SUM = NULL;
   TList *mergedListInSubsetLYZ2SUM = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetLYZ2SUM->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetLYZ2SUM->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetLYZ2SUM = TFile::Open(mergedFileNameInSubsetLYZ2SUM->Data(),"READ");
      delete mergedFileNameInSubsetLYZ2SUM;
      if(mergedFileinSubsetLYZ2SUM) 
      {
       mergedFileinSubsetLYZ2SUM->GetObject("cobjLYZ2SUM",mergedListInSubsetLYZ2SUM);
       mergedFileinSubsetLYZ2SUM->Close();
      } 
     }
     if(mergedListInSubsetLYZ2SUM)
     {
      AliFlowAnalysisWithLeeYangZeros *lyz2sum = new AliFlowAnalysisWithLeeYangZeros();
      lyz2sum->SetFirstRun(kFALSE);   
      lyz2sum->SetUseSum(kTRUE);     
      lyz2sum->GetOutputHistograms(mergedListInSubsetLYZ2SUM);
      lyz2sum->Finish();
      // save the final results for LYZ2SUM in final output file: 
      TString *finalOutputFileNameinSubsetLYZ2SUM = new TString("subset");
      (*finalOutputFileNameinSubsetLYZ2SUM)+=(i+1);  
      (*finalOutputFileNameinSubsetLYZ2SUM)+="/";
      (*finalOutputFileNameinSubsetLYZ2SUM)+="outputLYZ2SUManalysis";
      (*finalOutputFileNameinSubsetLYZ2SUM)+=type.Data();
      (*finalOutputFileNameinSubsetLYZ2SUM)+=".root";
      TFile *finalOutputInSubsetLYZ2SUM = new TFile(finalOutputFileNameinSubsetLYZ2SUM->Data(),"NEW");
      mergedListInSubsetLYZ2SUM->SetName("cobjLYZ2SUM");
      mergedListInSubsetLYZ2SUM->Write(mergedListInSubsetLYZ2SUM->GetName(),TObject::kSingleKey);
      finalOutputInSubsetLYZ2SUM->Close();
      delete finalOutputFileNameinSubsetLYZ2SUM;
      delete finalOutputInSubsetLYZ2SUM;
      delete lyz2sum;
     } else 
       {
        cout<<"WARNING: mergedListLYZ2SUM is NULL !!!!"<<endl;
       }                    
                                                     
   // LYZ2PROD:
   TString *mergedFileNameInSubsetLYZ2PROD = new TString("subset"); 
   (*mergedFileNameInSubsetLYZ2PROD)+=(i+1);
   (*mergedFileNameInSubsetLYZ2PROD)+="/";
   (*mergedFileNameInSubsetLYZ2PROD)+="mergedLYZ2PRODanalysis";
   (*mergedFileNameInSubsetLYZ2PROD)+=type.Data();
   (*mergedFileNameInSubsetLYZ2PROD)+=".root";
   
   TFile *mergedFileInSubsetLYZ2PROD = NULL;
   TList *mergedListInSubsetLYZ2PROD = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetLYZ2PROD->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetLYZ2PROD->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetLYZ2PROD = TFile::Open(mergedFileNameInSubsetLYZ2PROD->Data(),"READ");
      delete mergedFileNameInSubsetLYZ2PROD;
      if(mergedFileinSubsetLYZ2PROD) 
      {
       mergedFileinSubsetLYZ2PROD->GetObject("cobjLYZ2PROD",mergedListInSubsetLYZ2PROD);
       mergedFileinSubsetLYZ2PROD->Close();
      } 
     }
     if(mergedListInSubsetLYZ2PROD)
     {
      AliFlowAnalysisWithLeeYangZeros *lyz2prod = new AliFlowAnalysisWithLeeYangZeros();
      lyz2prod->SetFirstRun(kFALSE);   
      lyz2prod->SetUseSum(kFALSE);     
      lyz2prod->GetOutputHistograms(mergedListInSubsetLYZ2PROD);
      lyz2prod->Finish();
      // save the final results for LYZ2PROD in final output file: 
      TString *finalOutputFileNameinSubsetLYZ2PROD = new TString("subset");
      (*finalOutputFileNameinSubsetLYZ2PROD)+=(i+1);  
      (*finalOutputFileNameinSubsetLYZ2PROD)+="/";
      (*finalOutputFileNameinSubsetLYZ2PROD)+="outputLYZ2PRODanalysis";
      (*finalOutputFileNameinSubsetLYZ2PROD)+=type.Data();
      (*finalOutputFileNameinSubsetLYZ2PROD)+=".root";
      TFile *finalOutputInSubsetLYZ2PROD = new TFile(finalOutputFileNameinSubsetLYZ2PROD->Data(),"NEW");
      mergedListInSubsetLYZ2PROD->SetName("cobjLYZ2PROD");
      mergedListInSubsetLYZ2PROD->Write(mergedListInSubsetLYZ2PROD->GetName(),TObject::kSingleKey);
      finalOutputInSubsetLYZ2PROD->Close();
      delete finalOutputFileNameinSubsetLYZ2PROD;
      delete finalOutputInSubsetLYZ2PROD;
      delete lyz2prod;
     } else 
       {
        cout<<"WARNING: mergedListLYZ2PROD is NULL !!!!"<<endl;
       }                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
   // LYZEP:
   TString *mergedFileNameInSubsetLYZEP = new TString("subset"); 
   (*mergedFileNameInSubsetLYZEP)+=(i+1);
   (*mergedFileNameInSubsetLYZEP)+="/";
   (*mergedFileNameInSubsetLYZEP)+="mergedLYZEPanalysis";
   (*mergedFileNameInSubsetLYZEP)+=type.Data();
   (*mergedFileNameInSubsetLYZEP)+=".root";
   
   TFile *mergedFileInSubsetLYZEP = NULL;
   TList *mergedListInSubsetLYZEP = NULL;
   if(gSystem->AccessPathName(mergedFileNameInSubsetLYZEP->Data(),kFileExists))  
   {
    cout<<"WARNING: You do not have a merged output file "<<mergedFileNameInSubsetLYZEP->Data()<<endl;
    cout<<"         Couldn't redo Finish() on this file in subdirectory subset"<<i+1<<endl;        
   } else 
     {
      mergedFileinSubsetLYZEP = TFile::Open(mergedFileNameInSubsetLYZEP->Data(),"READ");
      delete mergedFileNameInSubsetLYZEP;
      if(mergedFileinSubsetLYZEP) 
      {
       mergedFileinSubsetLYZEP->GetObject("cobjLYZEP",mergedListInSubsetLYZEP);
       mergedFileinSubsetLYZEP->Close();
      } 
     }
     if(mergedListInSubsetLYZEP)
     {
      AliFlowAnalysisWithLYZEventPlane *lyzep = new AliFlowAnalysisWithLYZEventPlane();
      lyzep->GetOutputHistograms(mergedListInSubsetLYZEP);
      lyzep->Finish(kTRUE);
      // save the final results for LYZEP in final output file: 
      TString *finalOutputFileNameinSubsetLYZEP = new TString("subset");
      (*finalOutputFileNameinSubsetLYZEP)+=(i+1);  
      (*finalOutputFileNameinSubsetLYZEP)+="/";
      (*finalOutputFileNameinSubsetLYZEP)+="outputLYZEPanalysis";
      (*finalOutputFileNameinSubsetLYZEP)+=type.Data();
      (*finalOutputFileNameinSubsetLYZEP)+=".root";
      TFile *finalOutputInSubsetLYZEP = new TFile(finalOutputFileNameinSubsetLYZEP->Data(),"NEW");
      mergedListInSubsetLYZEP->SetName("cobjLYZEP");
      mergedListInSubsetLYZEP->Write(mergedListInSubsetLYZEP->GetName(),TObject::kSingleKey);
      finalOutputInSubsetLYZEP->Close();
      delete finalOutputFileNameinSubsetLYZEP;
      delete finalOutputInSubsetLYZEP;
      delete lyzep;
     } else 
       {
        cout<<"WARNING: mergedListLYZEP is NULL !!!!"<<endl;
       }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
   } // end of for(Int_t i=0,i<nSubsets;i++)
  
  delete baseDirPath;
  delete baseDir; 
  
 } // end of if(redoFinishAlsoInSubsets)

 
} // end of void reCallFinish(TString type="", Int_t mode=mLocal)
 
void LoadLibrariesRDF(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
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
  gSystem->Load("libTree");

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

 
