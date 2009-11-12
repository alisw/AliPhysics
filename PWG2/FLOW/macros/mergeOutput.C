enum libModes {mLocal,mLocalSource};

void mergeOutput(TString type="", Int_t mode=mLocal)
{
 // type:  type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //        (if type="" than output files are from MC simulation (default))
 // mode:  if mode=mLocal than analyze data on your computer using aliroot
 //        if mode=mLocalSource than analyze data on your computer using root + source files

 // settings for merging:
 Bool_t devideOutputFilesInEqualSubsets = kFALSE; // if kFALSE: All output files will be merged in a single output file
                                                 //            (this is deafult for merging on Grid).
                                                 // if kTRUE: All output files will be merged in a single output file, but ALSO
                                                 //           all output files will be devided in N equal subsets (nSubsets) and for 
                                                 //           each subset merging will be done separately. Merged output files for subsets
                                                 //           are stored in subdirectories subset1, subset2, ..., subsetN.  
                                                 //           This feature is used to estimate the spread. 
                                                                                                         
 const Int_t nSubsets = 10; //number od subsets used to estimate spread
 const Int_t cycle = 10; // merging is done in cycles and this is the cycle period 

 // default names of subdirectories in which merged output files for subsets are stored:
 TString *dirNameSubset = NULL;
 if(devideOutputFilesInEqualSubsets)
 {
  dirNameSubset = new TString("subset");  
 }
 
 if(cycle>100 && !devideOutputFilesInEqualSubsets)
 {
  cout<<"WARNING: Cycle is too big !!!!"<<endl; 
  cout<<"         Set const Int_t cycle to smaller value in the macro."<<endl;
  exit(0);
 }
 
 if(cycle>10 && devideOutputFilesInEqualSubsets)
 {
  cout<<"WARNING: Cycle is too big !!!!"<<endl; 
  cout<<"         Set const Int_t cycle to smaller value in the macro."<<endl;
  exit(0);
 }

 if(cycle<=0)
 {
  cout<<"WARNING: Cycle must be a positive integer !!!!"<<endl;
  cout<<"         Set const Int_t cycle to a positive integer in the macro."<<endl;
  exit(0);
 }
 
 if(devideOutputFilesInEqualSubsets && nSubsets<=0)
 {
  cout<<"WARNING: Number od subsets must be a positive integer !!!!"<<endl;
  cout<<"         Set const Int_t nSubsets to a positive integer in the macro."<<endl;
  exit(0);
 }
 
 // load needed libraries:                       
 LoadLibrariesMO(mode);  

 // standard magic:
 TString *baseDirPath = new TString(gSystem->pwd());
 TSystemDirectory* baseDir = new TSystemDirectory(".",baseDirPath->Data());          
 TList* listOfFilesInBaseDir = baseDir->GetListOfFiles();
 // listOfFilesInBaseDir->Print();
 Int_t nFiles = listOfFilesInBaseDir->GetEntries();
 
 // checking if subdirectories with default names for subsets already exist: 
 gSystem->cd(baseDirPath->Data()); 
 if(devideOutputFilesInEqualSubsets)
 {
  for(Int_t iFile=0;iFile<nFiles;iFile++)
  {
   TSystemFile* presentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
   TString presentFileName = baseDirPath->Data();
   (presentFileName+="/")+=presentFile->GetName(); 
   if(presentFileName.Contains("subset"))
   {
    cout<<"WARNING: You already have directory "<<presentFileName.Data()<<" !!!!"<<endl;
    cout<<"         Remove this directory and start macro again."<<endl;
    exit(0);
   } 
  } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)
 } // end of if(devideOutputFilesInEqualSubsets)

 // create subdirectories for storing merged results for subsets:
 gSystem->cd(baseDirPath->Data());
 if(devideOutputFilesInEqualSubsets)
 {
  for(Int_t i=0;i<nSubsets;i++)
  {
   TString currentDirForSpreadName(dirNameSubset->Data());
   Bool_t bDirectoryMade = gSystem->mkdir((currentDirForSpreadName += (i+1)).Data());
   if(bDirectoryMade)
   {
    cout<<"WARNING: The directory "<<currentDirForSpreadName.Data()<<" could not be created (illegal path name) !!!!"<<endl;
    exit(0);
   }
  } 
 }
 
 // loop over directories and add output *.root files to TFileMerger:
 gSystem->cd(baseDirPath->Data());
 
 // MCEP:
 Int_t fileCounterMCEP = 0;
 Int_t fileCounterLastSubsetMCEP = 0;
 TFileMerger *fileMergerMCEP = new TFileMerger();
 TFileMerger fileMergerForSubsetsMCEP[nSubsets];
 
 // SP:
 Int_t fileCounterSP = 0;
 Int_t fileCounterLastSubsetSP = 0;
 TFileMerger *fileMergerSP = new TFileMerger();
 TFileMerger fileMergerForSubsetsSP[nSubsets];
 
 // GFC:
 Int_t fileCounterGFC = 0;
 Int_t fileCounterLastSubsetGFC = 0;
 TFileMerger *fileMergerGFC = new TFileMerger();
 TFileMerger fileMergerForSubsetsGFC[nSubsets];
 
 // QC:
 Int_t fileCounterQC = 0;
 Int_t fileCounterLastSubsetQC = 0;
 TFileMerger *fileMergerQC = new TFileMerger();
 TFileMerger fileMergerForSubsetsQC[nSubsets];
 
 // FQD:
 Int_t fileCounterFQD = 0;
 Int_t fileCounterLastSubsetFQD = 0;
 TFileMerger *fileMergerFQD = new TFileMerger();
 TFileMerger fileMergerForSubsetsFQD[nSubsets];
 
 // LYZ1SUM:
 Int_t fileCounterLYZ1SUM = 0;
 Int_t fileCounterLastSubsetLYZ1SUM = 0;
 TFileMerger *fileMergerLYZ1SUM = new TFileMerger();
 TFileMerger fileMergerForSubsetsLYZ1SUM[nSubsets];
 
 // LYZ1PROD:
 Int_t fileCounterLYZ1PROD = 0;
 Int_t fileCounterLastSubsetLYZ1PROD = 0;
 TFileMerger *fileMergerLYZ1PROD = new TFileMerger();
 TFileMerger fileMergerForSubsetsLYZ1PROD[nSubsets];
 
 // LYZEP:
 Int_t fileCounterLYZEP = 0;
 Int_t fileCounterLastSubsetLYZEP = 0;
 TFileMerger *fileMergerLYZEP = new TFileMerger();
 TFileMerger fileMergerForSubsetsLYZEP[nSubsets];
 
 for(Int_t iFile=0;iFile<nFiles;iFile++)
 {
  TSystemFile* presentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
  
  if(!presentFile || !presentFile->IsDirectory() || strcmp(presentFile->GetName(), ".") == 0 || 
     strcmp(presentFile->GetName(), "..") == 0) continue; 
  
  TString presentDirName = baseDirPath->Data();
  (presentDirName+="/")+=presentFile->GetName();
  presentDirName+="/";
  
  // accessing the output *.root file in this directory:
  // MCEP:
  TString fileNameMCEP = presentDirName; 
  ((fileNameMCEP+="outputMCEPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameMCEP.Data(),kFileExists)))
  {
   fileCounterMCEP++;
   fileMergerMCEP->AddFile(fileNameMCEP.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsMCEP[(fileCounterMCEP-1) % nSubsets].AddFile(fileNameMCEP.Data());
    if(fileCounterMCEP % nSubsets == 0)
    {
     fileCounterLastSubsetMCEP++;
    }
   } 
   
   // global merging:
   if(fileCounterMCEP % cycle == 0)
   {
    TString *mergedFileForPreviousCycleMCEP = new TString("mergedMCEPanalysis"); 
    (*mergedFileForPreviousCycleMCEP)+=type.Data();
    (*mergedFileForPreviousCycleMCEP)+=(fileCounterMCEP/cycle - 1);
    (*mergedFileForPreviousCycleMCEP)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleMCEP->Data(),kFileExists)))
    {
     fileMergerMCEP->AddFile(mergedFileForPreviousCycleMCEP->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleMCEP->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleMCEP = new TString("mergedMCEPanalysis"); 
    (*mergedFileForCurrentCycleMCEP)+=type.Data();
    (*mergedFileForCurrentCycleMCEP)+=(fileCounterMCEP/cycle);
    (*mergedFileForCurrentCycleMCEP)+=".root";
    
    fileMergerMCEP->OutputFile(mergedFileForCurrentCycleMCEP->Data());
    fileMergerMCEP->Merge();
    fileMergerMCEP->Reset();
    
    delete mergedFileForPreviousCycleMCEP;
    delete mergedFileForCurrentCycleMCEP;
   } // end of if(fileCounterMCEP % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetMCEP > 0 && fileCounterLastSubsetMCEP % cycle == 0 && fileCounterMCEP % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetMCEP[nSubsets];
     TString mergedFileForCurrentCycleForSubsetMCEP[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetMCEP[i]="mergedMCEPanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetMCEP[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetMCEP[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetMCEP[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetMCEP[i]+=(fileCounterLastSubsetMCEP/cycle - 1);      
      mergedFileForPreviousCycleForSubsetMCEP[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetMCEP[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsMCEP[i].AddFile(mergedFileForPreviousCycleForSubsetMCEP[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetMCEP[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetMCEP[i]="mergedMCEPanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetMCEP[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetMCEP[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetMCEP[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetMCEP[i]+=(fileCounterLastSubsetMCEP/cycle);
      mergedFileForCurrentCycleForSubsetMCEP[i]+=".root";
      
      fileMergerForSubsetsMCEP[i].OutputFile(mergedFileForCurrentCycleForSubsetMCEP[i].Data());
      fileMergerForSubsetsMCEP[i].Merge();
      fileMergerForSubsetsMCEP[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetMCEP % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameMCEP.Data(),kFileExists))) 
  
  // SP:
  TString fileNameSP = presentDirName; 
  ((fileNameSP+="outputSPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameSP.Data(),kFileExists)))
  {
   fileCounterSP++;
   fileMergerSP->AddFile(fileNameSP.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsSP[(fileCounterSP-1) % nSubsets].AddFile(fileNameSP.Data());
    if(fileCounterSP % nSubsets == 0)
    {
     fileCounterLastSubsetSP++;
    }
   } 
   
   // global merging:
   if(fileCounterSP % cycle == 0)
   {
    TString *mergedFileForPreviousCycleSP = new TString("mergedSPanalysis"); 
    (*mergedFileForPreviousCycleSP)+=type.Data();
    (*mergedFileForPreviousCycleSP)+=(fileCounterSP/cycle - 1);
    (*mergedFileForPreviousCycleSP)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleSP->Data(),kFileExists)))
    {
     fileMergerSP->AddFile(mergedFileForPreviousCycleSP->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleSP->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleSP = new TString("mergedSPanalysis"); 
    (*mergedFileForCurrentCycleSP)+=type.Data();
    (*mergedFileForCurrentCycleSP)+=(fileCounterSP/cycle);
    (*mergedFileForCurrentCycleSP)+=".root";
    
    fileMergerSP->OutputFile(mergedFileForCurrentCycleSP->Data());
    fileMergerSP->Merge();
    fileMergerSP->Reset();
    
    delete mergedFileForPreviousCycleSP;
    delete mergedFileForCurrentCycleSP;
   } // end of if(fileCounterSP % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetSP > 0 && fileCounterLastSubsetSP % cycle == 0  && fileCounterSP % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetSP[nSubsets];
     TString mergedFileForCurrentCycleForSubsetSP[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetSP[i]="mergedSPanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetSP[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetSP[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetSP[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetSP[i]+=(fileCounterLastSubsetSP/cycle - 1);      
      mergedFileForPreviousCycleForSubsetSP[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetSP[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsSP[i].AddFile(mergedFileForPreviousCycleForSubsetSP[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetSP[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetSP[i]="mergedSPanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetSP[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetSP[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetSP[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetSP[i]+=(fileCounterLastSubsetSP/cycle);
      mergedFileForCurrentCycleForSubsetSP[i]+=".root";
      
      fileMergerForSubsetsSP[i].OutputFile(mergedFileForCurrentCycleForSubsetSP[i].Data());
      fileMergerForSubsetsSP[i].Merge();
      fileMergerForSubsetsSP[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetSP % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameSP.Data(),kFileExists)))
  
  // GFC:
  TString fileNameGFC = presentDirName; 
  ((fileNameGFC+="outputGFCanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameGFC.Data(),kFileExists)))
  {
   fileCounterGFC++;
   fileMergerGFC->AddFile(fileNameGFC.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsGFC[(fileCounterGFC-1) % nSubsets].AddFile(fileNameGFC.Data());
    if(fileCounterGFC % nSubsets == 0)
    {
     fileCounterLastSubsetGFC++;
    }
   } 
   
   // global merging:
   if(fileCounterGFC % cycle == 0)
   {
    TString *mergedFileForPreviousCycleGFC = new TString("mergedGFCanalysis"); 
    (*mergedFileForPreviousCycleGFC)+=type.Data();
    (*mergedFileForPreviousCycleGFC)+=(fileCounterGFC/cycle - 1);
    (*mergedFileForPreviousCycleGFC)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleGFC->Data(),kFileExists)))
    {
     fileMergerGFC->AddFile(mergedFileForPreviousCycleGFC->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleGFC->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleGFC = new TString("mergedGFCanalysis"); 
    (*mergedFileForCurrentCycleGFC)+=type.Data();
    (*mergedFileForCurrentCycleGFC)+=(fileCounterGFC/cycle);
    (*mergedFileForCurrentCycleGFC)+=".root";
    
    fileMergerGFC->OutputFile(mergedFileForCurrentCycleGFC->Data());
    fileMergerGFC->Merge();
    fileMergerGFC->Reset();
    
    delete mergedFileForPreviousCycleGFC;
    delete mergedFileForCurrentCycleGFC;
   } // end of if(fileCounterGFC % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetGFC > 0 && fileCounterLastSubsetGFC % cycle == 0  && fileCounterGFC % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetGFC[nSubsets];
     TString mergedFileForCurrentCycleForSubsetGFC[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetGFC[i]="mergedGFCanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetGFC[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetGFC[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetGFC[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetGFC[i]+=(fileCounterLastSubsetGFC/cycle - 1);      
      mergedFileForPreviousCycleForSubsetGFC[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetGFC[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsGFC[i].AddFile(mergedFileForPreviousCycleForSubsetGFC[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetGFC[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetGFC[i]="mergedGFCanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetGFC[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetGFC[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetGFC[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetGFC[i]+=(fileCounterLastSubsetGFC/cycle);
      mergedFileForCurrentCycleForSubsetGFC[i]+=".root";
      
      fileMergerForSubsetsGFC[i].OutputFile(mergedFileForCurrentCycleForSubsetGFC[i].Data());
      fileMergerForSubsetsGFC[i].Merge();
      fileMergerForSubsetsGFC[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetGFC % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameGFC.Data(),kFileExists)))
   
  // QC:
  TString fileNameQC = presentDirName; 
  ((fileNameQC+="outputQCanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameQC.Data(),kFileExists)))
  {
   fileCounterQC++;
   fileMergerQC->AddFile(fileNameQC.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsQC[(fileCounterQC-1) % nSubsets].AddFile(fileNameQC.Data());
    if(fileCounterQC % nSubsets == 0)
    {
     fileCounterLastSubsetQC++;
    }
   } 
   
   // global merging:
   if(fileCounterQC % cycle == 0)
   {
    TString *mergedFileForPreviousCycleQC = new TString("mergedQCanalysis"); 
    (*mergedFileForPreviousCycleQC)+=type.Data();
    (*mergedFileForPreviousCycleQC)+=(fileCounterQC/cycle - 1);
    (*mergedFileForPreviousCycleQC)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleQC->Data(),kFileExists)))
    {
     fileMergerQC->AddFile(mergedFileForPreviousCycleQC->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleQC->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleQC = new TString("mergedQCanalysis"); 
    (*mergedFileForCurrentCycleQC)+=type.Data();
    (*mergedFileForCurrentCycleQC)+=(fileCounterQC/cycle);
    (*mergedFileForCurrentCycleQC)+=".root";
    
    fileMergerQC->OutputFile(mergedFileForCurrentCycleQC->Data());
    fileMergerQC->Merge();
    fileMergerQC->Reset();
    
    delete mergedFileForPreviousCycleQC;
    delete mergedFileForCurrentCycleQC;
   } // end of if(fileCounterQC % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetQC > 0 && fileCounterLastSubsetQC % cycle == 0  && fileCounterQC % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetQC[nSubsets];
     TString mergedFileForCurrentCycleForSubsetQC[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetQC[i]="mergedQCanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetQC[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetQC[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetQC[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetQC[i]+=(fileCounterLastSubsetQC/cycle - 1);      
      mergedFileForPreviousCycleForSubsetQC[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetQC[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsQC[i].AddFile(mergedFileForPreviousCycleForSubsetQC[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetQC[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetQC[i]="mergedQCanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetQC[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetQC[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetQC[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetQC[i]+=(fileCounterLastSubsetQC/cycle);
      mergedFileForCurrentCycleForSubsetQC[i]+=".root";
      
      fileMergerForSubsetsQC[i].OutputFile(mergedFileForCurrentCycleForSubsetQC[i].Data());
      fileMergerForSubsetsQC[i].Merge();
      fileMergerForSubsetsQC[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetQC % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameQC.Data(),kFileExists)))
  
  // FQD:
  TString fileNameFQD = presentDirName; 
  ((fileNameFQD+="outputFQDanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameFQD.Data(),kFileExists)))
  {
   fileCounterFQD++;
   fileMergerFQD->AddFile(fileNameFQD.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsFQD[(fileCounterFQD-1) % nSubsets].AddFile(fileNameFQD.Data());
    if(fileCounterFQD % nSubsets == 0)
    {
     fileCounterLastSubsetFQD++;
    }
   } 
   
   // global merging:
   if(fileCounterFQD % cycle == 0)
   {
    TString *mergedFileForPreviousCycleFQD = new TString("mergedFQDanalysis"); 
    (*mergedFileForPreviousCycleFQD)+=type.Data();
    (*mergedFileForPreviousCycleFQD)+=(fileCounterFQD/cycle - 1);
    (*mergedFileForPreviousCycleFQD)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleFQD->Data(),kFileExists)))
    {
     fileMergerFQD->AddFile(mergedFileForPreviousCycleFQD->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleFQD->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleFQD = new TString("mergedFQDanalysis"); 
    (*mergedFileForCurrentCycleFQD)+=type.Data();
    (*mergedFileForCurrentCycleFQD)+=(fileCounterFQD/cycle);
    (*mergedFileForCurrentCycleFQD)+=".root";
    
    fileMergerFQD->OutputFile(mergedFileForCurrentCycleFQD->Data());
    fileMergerFQD->Merge();
    fileMergerFQD->Reset();
    
    delete mergedFileForPreviousCycleFQD;
    delete mergedFileForCurrentCycleFQD;
   } // end of if(fileCounterFQD % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetFQD > 0 && fileCounterLastSubsetFQD % cycle == 0  && fileCounterFQD % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetFQD[nSubsets];
     TString mergedFileForCurrentCycleForSubsetFQD[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetFQD[i]="mergedFQDanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetFQD[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetFQD[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetFQD[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetFQD[i]+=(fileCounterLastSubsetFQD/cycle - 1);      
      mergedFileForPreviousCycleForSubsetFQD[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetFQD[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsFQD[i].AddFile(mergedFileForPreviousCycleForSubsetFQD[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetFQD[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetFQD[i]="mergedFQDanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetFQD[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetFQD[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetFQD[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetFQD[i]+=(fileCounterLastSubsetFQD/cycle);
      mergedFileForCurrentCycleForSubsetFQD[i]+=".root";
      
      fileMergerForSubsetsFQD[i].OutputFile(mergedFileForCurrentCycleForSubsetFQD[i].Data());
      fileMergerForSubsetsFQD[i].Merge();
      fileMergerForSubsetsFQD[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetFQD % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameFQD.Data(),kFileExists)))
  
  // LYZ1SUM:
  TString fileNameLYZ1SUM = presentDirName; 
  ((fileNameLYZ1SUM+="outputLYZ1SUManalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZ1SUM.Data(),kFileExists)))
  {
   fileCounterLYZ1SUM++;
   fileMergerLYZ1SUM->AddFile(fileNameLYZ1SUM.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsLYZ1SUM[(fileCounterLYZ1SUM-1) % nSubsets].AddFile(fileNameLYZ1SUM.Data());
    if(fileCounterLYZ1SUM % nSubsets == 0)
    {
     fileCounterLastSubsetLYZ1SUM++;
    }
   } 
   
   // global merging:
   if(fileCounterLYZ1SUM % cycle == 0)
   {
    TString *mergedFileForPreviousCycleLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
    (*mergedFileForPreviousCycleLYZ1SUM)+=type.Data();
    (*mergedFileForPreviousCycleLYZ1SUM)+=(fileCounterLYZ1SUM/cycle - 1);
    (*mergedFileForPreviousCycleLYZ1SUM)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleLYZ1SUM->Data(),kFileExists)))
    {
     fileMergerLYZ1SUM->AddFile(mergedFileForPreviousCycleLYZ1SUM->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleLYZ1SUM->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
    (*mergedFileForCurrentCycleLYZ1SUM)+=type.Data();
    (*mergedFileForCurrentCycleLYZ1SUM)+=(fileCounterLYZ1SUM/cycle);
    (*mergedFileForCurrentCycleLYZ1SUM)+=".root";
    
    fileMergerLYZ1SUM->OutputFile(mergedFileForCurrentCycleLYZ1SUM->Data());
    fileMergerLYZ1SUM->Merge();
    fileMergerLYZ1SUM->Reset();
    
    delete mergedFileForPreviousCycleLYZ1SUM;
    delete mergedFileForCurrentCycleLYZ1SUM;
   } // end of if(fileCounterLYZ1SUM % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetLYZ1SUM > 0 && fileCounterLastSubsetLYZ1SUM % cycle == 0  && fileCounterLYZ1SUM % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetLYZ1SUM[nSubsets];
     TString mergedFileForCurrentCycleForSubsetLYZ1SUM[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]="mergedLYZ1SUManalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=(fileCounterLastSubsetLYZ1SUM/cycle - 1);      
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetLYZ1SUM[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsLYZ1SUM[i].AddFile(mergedFileForPreviousCycleForSubsetLYZ1SUM[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZ1SUM[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetLYZ1SUM[i]="mergedLYZ1SUManalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetLYZ1SUM[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetLYZ1SUM[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetLYZ1SUM[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetLYZ1SUM[i]+=(fileCounterLastSubsetLYZ1SUM/cycle);
      mergedFileForCurrentCycleForSubsetLYZ1SUM[i]+=".root";
      
      fileMergerForSubsetsLYZ1SUM[i].OutputFile(mergedFileForCurrentCycleForSubsetLYZ1SUM[i].Data());
      fileMergerForSubsetsLYZ1SUM[i].Merge();
      fileMergerForSubsetsLYZ1SUM[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetLYZ1SUM % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameLYZ1SUM.Data(),kFileExists)))
  
  // LYZ1PROD:
  TString fileNameLYZ1PROD = presentDirName; 
  ((fileNameLYZ1PROD+="outputLYZ1PRODanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZ1PROD.Data(),kFileExists)))
  {
   fileCounterLYZ1PROD++;
   fileMergerLYZ1PROD->AddFile(fileNameLYZ1PROD.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsLYZ1PROD[(fileCounterLYZ1PROD-1) % nSubsets].AddFile(fileNameLYZ1PROD.Data());
    if(fileCounterLYZ1PROD % nSubsets == 0)
    {
     fileCounterLastSubsetLYZ1PROD++;
    }
   } 
   
   // global merging:
   if(fileCounterLYZ1PROD % cycle == 0)
   {
    TString *mergedFileForPreviousCycleLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
    (*mergedFileForPreviousCycleLYZ1PROD)+=type.Data();
    (*mergedFileForPreviousCycleLYZ1PROD)+=(fileCounterLYZ1PROD/cycle - 1);
    (*mergedFileForPreviousCycleLYZ1PROD)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleLYZ1PROD->Data(),kFileExists)))
    {
     fileMergerLYZ1PROD->AddFile(mergedFileForPreviousCycleLYZ1PROD->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleLYZ1PROD->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
    (*mergedFileForCurrentCycleLYZ1PROD)+=type.Data();
    (*mergedFileForCurrentCycleLYZ1PROD)+=(fileCounterLYZ1PROD/cycle);
    (*mergedFileForCurrentCycleLYZ1PROD)+=".root";
    
    fileMergerLYZ1PROD->OutputFile(mergedFileForCurrentCycleLYZ1PROD->Data());
    fileMergerLYZ1PROD->Merge();
    fileMergerLYZ1PROD->Reset();
    
    delete mergedFileForPreviousCycleLYZ1PROD;
    delete mergedFileForCurrentCycleLYZ1PROD;
   } // end of if(fileCounterLYZ1PROD % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetLYZ1PROD > 0 && fileCounterLastSubsetLYZ1PROD % cycle == 0  && fileCounterLYZ1PROD % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetLYZ1PROD[nSubsets];
     TString mergedFileForCurrentCycleForSubsetLYZ1PROD[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]="mergedLYZ1PRODanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=(fileCounterLastSubsetLYZ1PROD/cycle - 1);      
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetLYZ1PROD[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsLYZ1PROD[i].AddFile(mergedFileForPreviousCycleForSubsetLYZ1PROD[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZ1PROD[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetLYZ1PROD[i]="mergedLYZ1PRODanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetLYZ1PROD[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetLYZ1PROD[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetLYZ1PROD[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetLYZ1PROD[i]+=(fileCounterLastSubsetLYZ1PROD/cycle);
      mergedFileForCurrentCycleForSubsetLYZ1PROD[i]+=".root";
      
      fileMergerForSubsetsLYZ1PROD[i].OutputFile(mergedFileForCurrentCycleForSubsetLYZ1PROD[i].Data());
      fileMergerForSubsetsLYZ1PROD[i].Merge();
      fileMergerForSubsetsLYZ1PROD[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetLYZ1PROD % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameLYZ1PROD.Data(),kFileExists)))
   
  // LYZEP:
  TString fileNameLYZEP = presentDirName; 
  ((fileNameLYZEP+="outputLYZEPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZEP.Data(),kFileExists)))
  {
   fileCounterLYZEP++;
   fileMergerLYZEP->AddFile(fileNameLYZEP.Data());
   
   if(devideOutputFilesInEqualSubsets)
   {
    fileMergerForSubsetsLYZEP[(fileCounterLYZEP-1) % nSubsets].AddFile(fileNameLYZEP.Data());
    if(fileCounterLYZEP % nSubsets == 0)
    {
     fileCounterLastSubsetLYZEP++;
    }
   } 
   
   // global merging:
   if(fileCounterLYZEP % cycle == 0)
   {
    TString *mergedFileForPreviousCycleLYZEP = new TString("mergedLYZEPanalysis"); 
    (*mergedFileForPreviousCycleLYZEP)+=type.Data();
    (*mergedFileForPreviousCycleLYZEP)+=(fileCounterLYZEP/cycle - 1);
    (*mergedFileForPreviousCycleLYZEP)+=".root";
    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleLYZEP->Data(),kFileExists)))
    {
     fileMergerLYZEP->AddFile(mergedFileForPreviousCycleLYZEP->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleLYZEP->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    TString *mergedFileForCurrentCycleLYZEP = new TString("mergedLYZEPanalysis"); 
    (*mergedFileForCurrentCycleLYZEP)+=type.Data();
    (*mergedFileForCurrentCycleLYZEP)+=(fileCounterLYZEP/cycle);
    (*mergedFileForCurrentCycleLYZEP)+=".root";
    
    fileMergerLYZEP->OutputFile(mergedFileForCurrentCycleLYZEP->Data());
    fileMergerLYZEP->Merge();
    fileMergerLYZEP->Reset();
    
    delete mergedFileForPreviousCycleLYZEP;
    delete mergedFileForCurrentCycleLYZEP;
   } // end of if(fileCounterLYZEP % cycle == 0) 
    
   // merging for subsets:
   if(devideOutputFilesInEqualSubsets)
   {
    if(fileCounterLastSubsetLYZEP > 0 && fileCounterLastSubsetLYZEP % cycle == 0  && fileCounterLYZEP % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetLYZEP[nSubsets];
     TString mergedFileForCurrentCycleForSubsetLYZEP[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetLYZEP[i]="mergedLYZEPanalysisForSubsetNo";
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetLYZEP[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=(fileCounterLastSubsetLYZEP/cycle - 1);      
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=".root";
      
      if(!(gSystem->AccessPathName(mergedFileForPreviousCycleForSubsetLYZEP[i].Data(),kFileExists)))
      {
       fileMergerForSubsetsLYZEP[i].AddFile(mergedFileForPreviousCycleForSubsetLYZEP[i].Data());
       // delete merged output from previous cycle:
       TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZEP[i].Data(),gSystem->pwd());
       file->Delete();
       delete file;
      }
      
      mergedFileForCurrentCycleForSubsetLYZEP[i]="mergedLYZEPanalysisForSubsetNo";
      mergedFileForCurrentCycleForSubsetLYZEP[i]+=(i+1);
      mergedFileForCurrentCycleForSubsetLYZEP[i]+=type.Data();
      mergedFileForCurrentCycleForSubsetLYZEP[i]+="CycleNo";
      mergedFileForCurrentCycleForSubsetLYZEP[i]+=(fileCounterLastSubsetLYZEP/cycle);
      mergedFileForCurrentCycleForSubsetLYZEP[i]+=".root";
      
      fileMergerForSubsetsLYZEP[i].OutputFile(mergedFileForCurrentCycleForSubsetLYZEP[i].Data());
      fileMergerForSubsetsLYZEP[i].Merge();
      fileMergerForSubsetsLYZEP[i].Reset();
       
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } // end of if(fileCounterLastSubsetLYZEP % cycle == 0)
   } // end of if(devideOutputFilesInEqualSubsets)
  } // end of if(!(gSystem->AccessPathName(fileNameLYZEP.Data(),kFileExists)))
  
 } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)
 
 
 //=================================================================================================
 
 
 // global merging at the end of the day:
 gSystem->cd(baseDirPath->Data());
 // MCEP:
 if(fileCounterMCEP < cycle)
 {
  TString *mergedFileFinalMCEP = new TString("mergedMCEPanalysis"); 
  (*mergedFileFinalMCEP)+=type.Data();
  (*mergedFileFinalMCEP)+=".root";
  fileMergerMCEP->OutputFile(mergedFileFinalMCEP->Data());
  fileMergerMCEP->Merge();
  delete mergedFileFinalMCEP;
 } else if (fileCounterMCEP % cycle == 0)
   {
    TString *mergedFileForPreviousCycleMCEP = new TString("mergedMCEPanalysis"); 
    (*mergedFileForPreviousCycleMCEP)+=type.Data();
    (*mergedFileForPreviousCycleMCEP)+=(fileCounterMCEP/cycle);
    (*mergedFileForPreviousCycleMCEP)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleMCEP->Data(),kFileExists)))
    {
     TString *mergedFileFinalMCEP = new TString("mergedMCEPanalysis"); 
     (*mergedFileFinalMCEP)+=type.Data();
     (*mergedFileFinalMCEP)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleMCEP->Data(),mergedFileFinalMCEP->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleMCEP = new TString("mergedMCEPanalysis"); 
      (*mergedFileForPreviousCycleMCEP)+=type.Data();
      (*mergedFileForPreviousCycleMCEP)+=((Int_t)fileCounterMCEP/cycle);
      (*mergedFileForPreviousCycleMCEP)+=".root";
      
      fileMergerMCEP->AddFile(mergedFileForPreviousCycleMCEP->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleMCEP->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalMCEP = new TString("mergedMCEPanalysis"); 
      (*mergedFileFinalMCEP)+=type.Data(); 
      (*mergedFileFinalMCEP)+=".root";
      fileMergerMCEP->OutputFile(mergedFileFinalMCEP->Data());
      fileMergerMCEP->Merge();
           
      delete mergedFileForPreviousCycleMCEP;
      delete mergedFileFinalMCEP;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetMCEP < cycle)
   {
    TString mergedFileFinalForSubsetMCEP[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetMCEP[i]="subset";
     mergedFileFinalForSubsetMCEP[i]+=(i+1);
     mergedFileFinalForSubsetMCEP[i]+="/";
     mergedFileFinalForSubsetMCEP[i]+="mergedMCEPanalysis";
     mergedFileFinalForSubsetMCEP[i]+=type.Data();
     mergedFileFinalForSubsetMCEP[i]+=".root";
       
     fileMergerForSubsetsMCEP[i].OutputFile(mergedFileFinalForSubsetMCEP[i].Data());
     fileMergerForSubsetsMCEP[i].Merge();
     fileMergerForSubsetsMCEP[i].Reset();  
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetMCEP % cycle == 0 && fileCounterMCEP % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetMCEP[nSubsets];
     TString mergedFileFinalForSubsetMCEP[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetMCEP[i]="mergedMCEPanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetMCEP[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetMCEP[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetMCEP[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetMCEP[i]+=(fileCounterLastSubsetMCEP/cycle);
      mergedFileForPreviousCycleForSubsetMCEP[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetMCEP[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedMCEPanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetMCEP[nSubsets];
       TString mergedFileFinalForSubsetMCEP[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetMCEP[i]="mergedMCEPanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetMCEP[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetMCEP[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetMCEP[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetMCEP[i]+=((Int_t)fileCounterLastSubsetMCEP/cycle);
        mergedFileForPreviousCycleForSubsetMCEP[i]+=".root";
      
        fileMergerForSubsetsMCEP[i].AddFile(mergedFileForPreviousCycleForSubsetMCEP[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetMCEP[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetMCEP[i]="subset";
        mergedFileFinalForSubsetMCEP[i]+=(i+1);
        mergedFileFinalForSubsetMCEP[i]+="/";
        mergedFileFinalForSubsetMCEP[i]+="mergedMCEPanalysis";
        mergedFileFinalForSubsetMCEP[i]+=type.Data();
        mergedFileFinalForSubsetMCEP[i]+=".root";
        
        fileMergerForSubsetsMCEP[i].OutputFile(mergedFileFinalForSubsetMCEP[i].Data());
        
        fileMergerForSubsetsMCEP[i].Merge();
        fileMergerForSubsetsMCEP[i].Reset();                                      
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets) 
   
 // SP:
 if(fileCounterSP < cycle)
 {
  TString *mergedFileFinalSP = new TString("mergedSPanalysis"); 
  (*mergedFileFinalSP)+=type.Data();
  (*mergedFileFinalSP)+=".root";
  fileMergerSP->OutputFile(mergedFileFinalSP->Data());
  fileMergerSP->Merge();
  delete mergedFileFinalSP;
 } else if (fileCounterSP % cycle == 0)
   {
    TString *mergedFileForPreviousCycleSP = new TString("mergedSPanalysis"); 
    (*mergedFileForPreviousCycleSP)+=type.Data();
    (*mergedFileForPreviousCycleSP)+=(fileCounterSP/cycle);
    (*mergedFileForPreviousCycleSP)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleSP->Data(),kFileExists)))
    {
     TString *mergedFileFinalSP = new TString("mergedSPanalysis"); 
     (*mergedFileFinalSP)+=type.Data();
     (*mergedFileFinalSP)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleSP->Data(),mergedFileFinalSP->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleSP = new TString("mergedSPanalysis"); 
      (*mergedFileForPreviousCycleSP)+=type.Data();
      (*mergedFileForPreviousCycleSP)+=((Int_t)fileCounterSP/cycle);
      (*mergedFileForPreviousCycleSP)+=".root";
      
      fileMergerSP->AddFile(mergedFileForPreviousCycleSP->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleSP->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalSP = new TString("mergedSPanalysis"); 
      (*mergedFileFinalSP)+=type.Data(); 
      (*mergedFileFinalSP)+=".root";
      fileMergerSP->OutputFile(mergedFileFinalSP->Data());
      fileMergerSP->Merge();
           
      delete mergedFileForPreviousCycleSP;
      delete mergedFileFinalSP;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetSP < cycle)
   {
    TString mergedFileFinalForSubsetSP[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetSP[i]="subset";
     mergedFileFinalForSubsetSP[i]+=(i+1);
     mergedFileFinalForSubsetSP[i]+="/";
     mergedFileFinalForSubsetSP[i]+="mergedSPanalysis";
     mergedFileFinalForSubsetSP[i]+=type.Data();
     mergedFileFinalForSubsetSP[i]+=".root";
       
     fileMergerForSubsetsSP[i].OutputFile(mergedFileFinalForSubsetSP[i].Data());
     fileMergerForSubsetsSP[i].Merge();
     fileMergerForSubsetsSP[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetSP % cycle == 0  && fileCounterSP % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetSP[nSubsets];
     TString mergedFileFinalForSubsetSP[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetSP[i]="mergedSPanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetSP[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetSP[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetSP[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetSP[i]+=(fileCounterLastSubsetSP/cycle);
      mergedFileForPreviousCycleForSubsetSP[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetSP[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedSPanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetSP[nSubsets];
       TString mergedFileFinalForSubsetSP[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetSP[i]="mergedSPanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetSP[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetSP[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetSP[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetSP[i]+=((Int_t)fileCounterLastSubsetSP/cycle);
        mergedFileForPreviousCycleForSubsetSP[i]+=".root";
      
        fileMergerForSubsetsSP[i].AddFile(mergedFileForPreviousCycleForSubsetSP[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetSP[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetSP[i]="subset";
        mergedFileFinalForSubsetSP[i]+=(i+1);
        mergedFileFinalForSubsetSP[i]+="/";
        mergedFileFinalForSubsetSP[i]+="mergedSPanalysis";
        mergedFileFinalForSubsetSP[i]+=type.Data();
        mergedFileFinalForSubsetSP[i]+=".root";
        
        fileMergerForSubsetsSP[i].OutputFile(mergedFileFinalForSubsetSP[i].Data());
        fileMergerForSubsetsSP[i].Merge();
        fileMergerForSubsetsSP[i].Reset();      
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets)   
        
 // GFC:
 if(fileCounterGFC < cycle)
 {
  TString *mergedFileFinalGFC = new TString("mergedGFCanalysis"); 
  (*mergedFileFinalGFC)+=type.Data();
  (*mergedFileFinalGFC)+=".root";
  fileMergerGFC->OutputFile(mergedFileFinalGFC->Data());
  fileMergerGFC->Merge();
  delete mergedFileFinalGFC;
 } else if (fileCounterGFC % cycle == 0)
   {
    TString *mergedFileForPreviousCycleGFC = new TString("mergedGFCanalysis"); 
    (*mergedFileForPreviousCycleGFC)+=type.Data();
    (*mergedFileForPreviousCycleGFC)+=(fileCounterGFC/cycle);
    (*mergedFileForPreviousCycleGFC)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleGFC->Data(),kFileExists)))
    {
     TString *mergedFileFinalGFC = new TString("mergedGFCanalysis"); 
     (*mergedFileFinalGFC)+=type.Data();
     (*mergedFileFinalGFC)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleGFC->Data(),mergedFileFinalGFC->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleGFC = new TString("mergedGFCanalysis"); 
      (*mergedFileForPreviousCycleGFC)+=type.Data();
      (*mergedFileForPreviousCycleGFC)+=((Int_t)fileCounterGFC/cycle);
      (*mergedFileForPreviousCycleGFC)+=".root";
      
      fileMergerGFC->AddFile(mergedFileForPreviousCycleGFC->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleGFC->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalGFC = new TString("mergedGFCanalysis"); 
      (*mergedFileFinalGFC)+=type.Data(); 
      (*mergedFileFinalGFC)+=".root";
      fileMergerGFC->OutputFile(mergedFileFinalGFC->Data());
      fileMergerGFC->Merge();
           
      delete mergedFileForPreviousCycleGFC;
      delete mergedFileFinalGFC;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetGFC < cycle)
   {
    TString mergedFileFinalForSubsetGFC[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetGFC[i]="subset";
     mergedFileFinalForSubsetGFC[i]+=(i+1);
     mergedFileFinalForSubsetGFC[i]+="/";
     mergedFileFinalForSubsetGFC[i]+="mergedGFCanalysis";
     mergedFileFinalForSubsetGFC[i]+=type.Data();
     mergedFileFinalForSubsetGFC[i]+=".root";
       
     fileMergerForSubsetsGFC[i].OutputFile(mergedFileFinalForSubsetGFC[i].Data());
     fileMergerForSubsetsGFC[i].Merge();
     fileMergerForSubsetsGFC[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetGFC % cycle == 0  && fileCounterGFC % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetGFC[nSubsets];
     TString mergedFileFinalForSubsetGFC[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetGFC[i]="mergedGFCanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetGFC[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetGFC[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetGFC[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetGFC[i]+=(fileCounterLastSubsetGFC/cycle);
      mergedFileForPreviousCycleForSubsetGFC[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetGFC[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedGFCanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetGFC[nSubsets];
       TString mergedFileFinalForSubsetGFC[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetGFC[i]="mergedGFCanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetGFC[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetGFC[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetGFC[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetGFC[i]+=((Int_t)fileCounterLastSubsetGFC/cycle);
        mergedFileForPreviousCycleForSubsetGFC[i]+=".root";
      
        fileMergerForSubsetsGFC[i].AddFile(mergedFileForPreviousCycleForSubsetGFC[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetGFC[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetGFC[i]="subset";
        mergedFileFinalForSubsetGFC[i]+=(i+1);
        mergedFileFinalForSubsetGFC[i]+="/";
        mergedFileFinalForSubsetGFC[i]+="mergedGFCanalysis";
        mergedFileFinalForSubsetGFC[i]+=type.Data();
        mergedFileFinalForSubsetGFC[i]+=".root";
        
        fileMergerForSubsetsGFC[i].OutputFile(mergedFileFinalForSubsetGFC[i].Data());
        fileMergerForSubsetsGFC[i].Merge();
        fileMergerForSubsetsGFC[i].Reset();         
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets)   
                             
 // QC:
 if(fileCounterQC < cycle)
 {
  TString *mergedFileFinalQC = new TString("mergedQCanalysis"); 
  (*mergedFileFinalQC)+=type.Data();
  (*mergedFileFinalQC)+=".root";
  fileMergerQC->OutputFile(mergedFileFinalQC->Data());
  fileMergerQC->Merge();
  delete mergedFileFinalQC;
 } else if (fileCounterQC % cycle == 0)
   {
    TString *mergedFileForPreviousCycleQC = new TString("mergedQCanalysis"); 
    (*mergedFileForPreviousCycleQC)+=type.Data();
    (*mergedFileForPreviousCycleQC)+=(fileCounterQC/cycle);
    (*mergedFileForPreviousCycleQC)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleQC->Data(),kFileExists)))
    {
     TString *mergedFileFinalQC = new TString("mergedQCanalysis"); 
     (*mergedFileFinalQC)+=type.Data();
     (*mergedFileFinalQC)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleQC->Data(),mergedFileFinalQC->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleQC = new TString("mergedQCanalysis"); 
      (*mergedFileForPreviousCycleQC)+=type.Data();
      (*mergedFileForPreviousCycleQC)+=((Int_t)fileCounterQC/cycle);
      (*mergedFileForPreviousCycleQC)+=".root";
      
      fileMergerQC->AddFile(mergedFileForPreviousCycleQC->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleQC->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalQC = new TString("mergedQCanalysis"); 
      (*mergedFileFinalQC)+=type.Data(); 
      (*mergedFileFinalQC)+=".root";
      fileMergerQC->OutputFile(mergedFileFinalQC->Data());
      fileMergerQC->Merge();
           
      delete mergedFileForPreviousCycleQC;
      delete mergedFileFinalQC;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetQC < cycle)
   {
    TString mergedFileFinalForSubsetQC[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetQC[i]="subset";
     mergedFileFinalForSubsetQC[i]+=(i+1);
     mergedFileFinalForSubsetQC[i]+="/";
     mergedFileFinalForSubsetQC[i]+="mergedQCanalysis";
     mergedFileFinalForSubsetQC[i]+=type.Data();
     mergedFileFinalForSubsetQC[i]+=".root";
       
     fileMergerForSubsetsQC[i].OutputFile(mergedFileFinalForSubsetQC[i].Data());
     fileMergerForSubsetsQC[i].Merge();
     fileMergerForSubsetsQC[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetQC % cycle == 0  && fileCounterQC % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetQC[nSubsets];
     TString mergedFileFinalForSubsetQC[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetQC[i]="mergedQCanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetQC[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetQC[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetQC[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetQC[i]+=(fileCounterLastSubsetQC/cycle);
      mergedFileForPreviousCycleForSubsetQC[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetQC[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedQCanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetQC[nSubsets];
       TString mergedFileFinalForSubsetQC[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetQC[i]="mergedQCanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetQC[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetQC[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetQC[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetQC[i]+=((Int_t)fileCounterLastSubsetQC/cycle);
        mergedFileForPreviousCycleForSubsetQC[i]+=".root";
      
        fileMergerForSubsetsQC[i].AddFile(mergedFileForPreviousCycleForSubsetQC[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetQC[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetQC[i]="subset";
        mergedFileFinalForSubsetQC[i]+=(i+1);
        mergedFileFinalForSubsetQC[i]+="/";
        mergedFileFinalForSubsetQC[i]+="mergedQCanalysis";
        mergedFileFinalForSubsetQC[i]+=type.Data();
        mergedFileFinalForSubsetQC[i]+=".root";
        
        fileMergerForSubsetsQC[i].OutputFile(mergedFileFinalForSubsetQC[i].Data());
        fileMergerForSubsetsQC[i].Merge();
        fileMergerForSubsetsQC[i].Reset();        
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets)                                  
                                                                                                                                                                         
 // FQD:
 if(fileCounterFQD < cycle)
 {
  TString *mergedFileFinalFQD = new TString("mergedFQDanalysis"); 
  (*mergedFileFinalFQD)+=type.Data();
  (*mergedFileFinalFQD)+=".root";
  fileMergerFQD->OutputFile(mergedFileFinalFQD->Data());
  fileMergerFQD->Merge();
  delete mergedFileFinalFQD;
 } else if (fileCounterFQD % cycle == 0)
   {
    TString *mergedFileForPreviousCycleFQD = new TString("mergedFQDanalysis"); 
    (*mergedFileForPreviousCycleFQD)+=type.Data();
    (*mergedFileForPreviousCycleFQD)+=(fileCounterFQD/cycle);
    (*mergedFileForPreviousCycleFQD)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleFQD->Data(),kFileExists)))
    {
     TString *mergedFileFinalFQD = new TString("mergedFQDanalysis"); 
     (*mergedFileFinalFQD)+=type.Data();
     (*mergedFileFinalFQD)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleFQD->Data(),mergedFileFinalFQD->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleFQD = new TString("mergedFQDanalysis"); 
      (*mergedFileForPreviousCycleFQD)+=type.Data();
      (*mergedFileForPreviousCycleFQD)+=((Int_t)fileCounterFQD/cycle);
      (*mergedFileForPreviousCycleFQD)+=".root";
      
      fileMergerFQD->AddFile(mergedFileForPreviousCycleFQD->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleFQD->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalFQD = new TString("mergedFQDanalysis"); 
      (*mergedFileFinalFQD)+=type.Data(); 
      (*mergedFileFinalFQD)+=".root";
      fileMergerFQD->OutputFile(mergedFileFinalFQD->Data());
      fileMergerFQD->Merge();
           
      delete mergedFileForPreviousCycleFQD;
      delete mergedFileFinalFQD;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetFQD < cycle)
   {
    TString mergedFileFinalForSubsetFQD[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetFQD[i]="subset";
     mergedFileFinalForSubsetFQD[i]+=(i+1);
     mergedFileFinalForSubsetFQD[i]+="/";
     mergedFileFinalForSubsetFQD[i]+="mergedFQDanalysis";
     mergedFileFinalForSubsetFQD[i]+=type.Data();
     mergedFileFinalForSubsetFQD[i]+=".root";
       
     fileMergerForSubsetsFQD[i].OutputFile(mergedFileFinalForSubsetFQD[i].Data());
     fileMergerForSubsetsFQD[i].Merge();
     fileMergerForSubsetsFQD[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetFQD % cycle == 0  && fileCounterFQD % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetFQD[nSubsets];
     TString mergedFileFinalForSubsetFQD[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetFQD[i]="mergedFQDanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetFQD[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetFQD[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetFQD[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetFQD[i]+=(fileCounterLastSubsetFQD/cycle);
      mergedFileForPreviousCycleForSubsetFQD[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetFQD[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedFQDanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetFQD[nSubsets];
       TString mergedFileFinalForSubsetFQD[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetFQD[i]="mergedFQDanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetFQD[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetFQD[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetFQD[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetFQD[i]+=((Int_t)fileCounterLastSubsetFQD/cycle);
        mergedFileForPreviousCycleForSubsetFQD[i]+=".root";
      
        fileMergerForSubsetsFQD[i].AddFile(mergedFileForPreviousCycleForSubsetFQD[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetFQD[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetFQD[i]="subset";
        mergedFileFinalForSubsetFQD[i]+=(i+1);
        mergedFileFinalForSubsetFQD[i]+="/";
        mergedFileFinalForSubsetFQD[i]+="mergedFQDanalysis";
        mergedFileFinalForSubsetFQD[i]+=type.Data();
        mergedFileFinalForSubsetFQD[i]+=".root";

        fileMergerForSubsetsFQD[i].OutputFile(mergedFileFinalForSubsetFQD[i].Data());
        fileMergerForSubsetsFQD[i].Merge();
        fileMergerForSubsetsFQD[i].Reset();                      
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets) 
                                                                                                     
 // LYZ1SUM:
 if(fileCounterLYZ1SUM < cycle)
 {
  TString *mergedFileFinalLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
  (*mergedFileFinalLYZ1SUM)+=type.Data();
  (*mergedFileFinalLYZ1SUM)+=".root";
  fileMergerLYZ1SUM->OutputFile(mergedFileFinalLYZ1SUM->Data());
  fileMergerLYZ1SUM->Merge();
  delete mergedFileFinalLYZ1SUM;
 } else if (fileCounterLYZ1SUM % cycle == 0)
   {
    TString *mergedFileForPreviousCycleLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
    (*mergedFileForPreviousCycleLYZ1SUM)+=type.Data();
    (*mergedFileForPreviousCycleLYZ1SUM)+=(fileCounterLYZ1SUM/cycle);
    (*mergedFileForPreviousCycleLYZ1SUM)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleLYZ1SUM->Data(),kFileExists)))
    {
     TString *mergedFileFinalLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
     (*mergedFileFinalLYZ1SUM)+=type.Data();
     (*mergedFileFinalLYZ1SUM)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleLYZ1SUM->Data(),mergedFileFinalLYZ1SUM->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
      (*mergedFileForPreviousCycleLYZ1SUM)+=type.Data();
      (*mergedFileForPreviousCycleLYZ1SUM)+=((Int_t)fileCounterLYZ1SUM/cycle);
      (*mergedFileForPreviousCycleLYZ1SUM)+=".root";
      
      fileMergerLYZ1SUM->AddFile(mergedFileForPreviousCycleLYZ1SUM->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleLYZ1SUM->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalLYZ1SUM = new TString("mergedLYZ1SUManalysis"); 
      (*mergedFileFinalLYZ1SUM)+=type.Data(); 
      (*mergedFileFinalLYZ1SUM)+=".root";
      fileMergerLYZ1SUM->OutputFile(mergedFileFinalLYZ1SUM->Data());
      fileMergerLYZ1SUM->Merge();
           
      delete mergedFileForPreviousCycleLYZ1SUM;
      delete mergedFileFinalLYZ1SUM;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetLYZ1SUM < cycle)
   {
    TString mergedFileFinalForSubsetLYZ1SUM[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetLYZ1SUM[i]="subset";
     mergedFileFinalForSubsetLYZ1SUM[i]+=(i+1);
     mergedFileFinalForSubsetLYZ1SUM[i]+="/";
     mergedFileFinalForSubsetLYZ1SUM[i]+="mergedLYZ1SUManalysis";
     mergedFileFinalForSubsetLYZ1SUM[i]+=type.Data();
     mergedFileFinalForSubsetLYZ1SUM[i]+=".root";
       
     fileMergerForSubsetsLYZ1SUM[i].OutputFile(mergedFileFinalForSubsetLYZ1SUM[i].Data());
     fileMergerForSubsetsLYZ1SUM[i].Merge();
     fileMergerForSubsetsLYZ1SUM[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetLYZ1SUM % cycle == 0  && fileCounterLYZ1SUM % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetLYZ1SUM[nSubsets];
     TString mergedFileFinalForSubsetLYZ1SUM[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]="mergedLYZ1SUManalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=(fileCounterLastSubsetLYZ1SUM/cycle);
      mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZ1SUM[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedLYZ1SUManalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetLYZ1SUM[nSubsets];
       TString mergedFileFinalForSubsetLYZ1SUM[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetLYZ1SUM[i]="mergedLYZ1SUManalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=((Int_t)fileCounterLastSubsetLYZ1SUM/cycle);
        mergedFileForPreviousCycleForSubsetLYZ1SUM[i]+=".root";
      
        fileMergerForSubsetsLYZ1SUM[i].AddFile(mergedFileForPreviousCycleForSubsetLYZ1SUM[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZ1SUM[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetLYZ1SUM[i]="subset";
        mergedFileFinalForSubsetLYZ1SUM[i]+=(i+1);
        mergedFileFinalForSubsetLYZ1SUM[i]+="/";
        mergedFileFinalForSubsetLYZ1SUM[i]+="mergedLYZ1SUManalysis";
        mergedFileFinalForSubsetLYZ1SUM[i]+=type.Data();
        mergedFileFinalForSubsetLYZ1SUM[i]+=".root";
        
        fileMergerForSubsetsLYZ1SUM[i].OutputFile(mergedFileFinalForSubsetLYZ1SUM[i].Data());
        fileMergerForSubsetsLYZ1SUM[i].Merge();
        fileMergerForSubsetsLYZ1SUM[i].Reset();       
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets)  
   
 // LYZ1PROD:
 if(fileCounterLYZ1PROD < cycle)
 {
  TString *mergedFileFinalLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
  (*mergedFileFinalLYZ1PROD)+=type.Data();
  (*mergedFileFinalLYZ1PROD)+=".root";
  fileMergerLYZ1PROD->OutputFile(mergedFileFinalLYZ1PROD->Data());
  fileMergerLYZ1PROD->Merge();
  delete mergedFileFinalLYZ1PROD;
 } else if (fileCounterLYZ1PROD % cycle == 0)
   {
    TString *mergedFileForPreviousCycleLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
    (*mergedFileForPreviousCycleLYZ1PROD)+=type.Data();
    (*mergedFileForPreviousCycleLYZ1PROD)+=(fileCounterLYZ1PROD/cycle);
    (*mergedFileForPreviousCycleLYZ1PROD)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleLYZ1PROD->Data(),kFileExists)))
    {
     TString *mergedFileFinalLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
     (*mergedFileFinalLYZ1PROD)+=type.Data();
     (*mergedFileFinalLYZ1PROD)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleLYZ1PROD->Data(),mergedFileFinalLYZ1PROD->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
      (*mergedFileForPreviousCycleLYZ1PROD)+=type.Data();
      (*mergedFileForPreviousCycleLYZ1PROD)+=((Int_t)fileCounterLYZ1PROD/cycle);
      (*mergedFileForPreviousCycleLYZ1PROD)+=".root";
      
      fileMergerLYZ1PROD->AddFile(mergedFileForPreviousCycleLYZ1PROD->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleLYZ1PROD->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalLYZ1PROD = new TString("mergedLYZ1PRODanalysis"); 
      (*mergedFileFinalLYZ1PROD)+=type.Data(); 
      (*mergedFileFinalLYZ1PROD)+=".root";
      fileMergerLYZ1PROD->OutputFile(mergedFileFinalLYZ1PROD->Data());
      fileMergerLYZ1PROD->Merge();
           
      delete mergedFileForPreviousCycleLYZ1PROD;
      delete mergedFileFinalLYZ1PROD;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetLYZ1PROD < cycle)
   {
    TString mergedFileFinalForSubsetLYZ1PROD[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetLYZ1PROD[i]="subset";
     mergedFileFinalForSubsetLYZ1PROD[i]+=(i+1);
     mergedFileFinalForSubsetLYZ1PROD[i]+="/";
     mergedFileFinalForSubsetLYZ1PROD[i]+="mergedLYZ1PRODanalysis";
     mergedFileFinalForSubsetLYZ1PROD[i]+=type.Data();
     mergedFileFinalForSubsetLYZ1PROD[i]+=".root";
       
     fileMergerForSubsetsLYZ1PROD[i].OutputFile(mergedFileFinalForSubsetLYZ1PROD[i].Data());
     fileMergerForSubsetsLYZ1PROD[i].Merge();
     fileMergerForSubsetsLYZ1PROD[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetLYZ1PROD % cycle == 0  && fileCounterLYZ1PROD % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetLYZ1PROD[nSubsets];
     TString mergedFileFinalForSubsetLYZ1PROD[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]="mergedLYZ1PRODanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=(fileCounterLastSubsetLYZ1PROD/cycle);
      mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZ1PROD[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedLYZ1PRODanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetLYZ1PROD[nSubsets];
       TString mergedFileFinalForSubsetLYZ1PROD[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetLYZ1PROD[i]="mergedLYZ1PRODanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=((Int_t)fileCounterLastSubsetLYZ1PROD/cycle);
        mergedFileForPreviousCycleForSubsetLYZ1PROD[i]+=".root";
      
        fileMergerForSubsetsLYZ1PROD[i].AddFile(mergedFileForPreviousCycleForSubsetLYZ1PROD[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZ1PROD[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetLYZ1PROD[i]="subset";
        mergedFileFinalForSubsetLYZ1PROD[i]+=(i+1);
        mergedFileFinalForSubsetLYZ1PROD[i]+="/";
        mergedFileFinalForSubsetLYZ1PROD[i]+="mergedLYZ1PRODanalysis";
        mergedFileFinalForSubsetLYZ1PROD[i]+=type.Data();
        mergedFileFinalForSubsetLYZ1PROD[i]+=".root";
        
        fileMergerForSubsetsLYZ1PROD[i].OutputFile(mergedFileFinalForSubsetLYZ1PROD[i].Data());
        fileMergerForSubsetsLYZ1PROD[i].Merge();
        fileMergerForSubsetsLYZ1PROD[i].Reset();      
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets) 
 
 // LYZEP:
 if(fileCounterLYZEP < cycle)
 {
  TString *mergedFileFinalLYZEP = new TString("mergedLYZEPanalysis"); 
  (*mergedFileFinalLYZEP)+=type.Data();
  (*mergedFileFinalLYZEP)+=".root";
  fileMergerLYZEP->OutputFile(mergedFileFinalLYZEP->Data());
  fileMergerLYZEP->Merge();
  delete mergedFileFinalLYZEP;
 } else if (fileCounterLYZEP % cycle == 0)
   {
    TString *mergedFileForPreviousCycleLYZEP = new TString("mergedLYZEPanalysis"); 
    (*mergedFileForPreviousCycleLYZEP)+=type.Data();
    (*mergedFileForPreviousCycleLYZEP)+=(fileCounterLYZEP/cycle);
    (*mergedFileForPreviousCycleLYZEP)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycleLYZEP->Data(),kFileExists)))
    {
     TString *mergedFileFinalLYZEP = new TString("mergedLYZEPanalysis"); 
     (*mergedFileFinalLYZEP)+=type.Data();
     (*mergedFileFinalLYZEP)+=".root";
     gSystem->Rename(mergedFileForPreviousCycleLYZEP->Data(),mergedFileFinalLYZEP->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycleLYZEP = new TString("mergedLYZEPanalysis"); 
      (*mergedFileForPreviousCycleLYZEP)+=type.Data();
      (*mergedFileForPreviousCycleLYZEP)+=((Int_t)fileCounterLYZEP/cycle);
      (*mergedFileForPreviousCycleLYZEP)+=".root";
      
      fileMergerLYZEP->AddFile(mergedFileForPreviousCycleLYZEP->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleLYZEP->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinalLYZEP = new TString("mergedLYZEPanalysis"); 
      (*mergedFileFinalLYZEP)+=type.Data(); 
      (*mergedFileFinalLYZEP)+=".root";
      fileMergerLYZEP->OutputFile(mergedFileFinalLYZEP->Data());
      fileMergerLYZEP->Merge();
           
      delete mergedFileForPreviousCycleLYZEP;
      delete mergedFileFinalLYZEP;
     }
 
  // merging for subsets at the end of the day:
  gSystem->cd(baseDirPath->Data());
  if(devideOutputFilesInEqualSubsets) 
  {
   if(fileCounterLastSubsetLYZEP < cycle)
   {
    TString mergedFileFinalForSubsetLYZEP[nSubsets];
    for(Int_t i=0;i<nSubsets;i++)
    {
     mergedFileFinalForSubsetLYZEP[i]="subset";
     mergedFileFinalForSubsetLYZEP[i]+=(i+1);
     mergedFileFinalForSubsetLYZEP[i]+="/";
     mergedFileFinalForSubsetLYZEP[i]+="mergedLYZEPanalysis";
     mergedFileFinalForSubsetLYZEP[i]+=type.Data();
     mergedFileFinalForSubsetLYZEP[i]+=".root";
       
     fileMergerForSubsetsLYZEP[i].OutputFile(mergedFileFinalForSubsetLYZEP[i].Data());
     fileMergerForSubsetsLYZEP[i].Merge();
     fileMergerForSubsetsLYZEP[i].Reset();
    } // end of for(Int_t i=0;i<nSubsets;i++) 
   } else if (fileCounterLastSubsetLYZEP % cycle == 0  && fileCounterLYZEP % cycle == 0)
    {
     TString mergedFileForPreviousCycleForSubsetLYZEP[nSubsets];
     TString mergedFileFinalForSubsetLYZEP[nSubsets];
     for(Int_t i=0;i<nSubsets;i++)
     {
      mergedFileForPreviousCycleForSubsetLYZEP[i]="mergedLYZEPanalysisForSubsetNo"; 
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=(i+1);
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=type.Data();
      mergedFileForPreviousCycleForSubsetLYZEP[i]+="CycleNo";
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=(fileCounterLastSubsetLYZEP/cycle);
      mergedFileForPreviousCycleForSubsetLYZEP[i]+=".root";
      // move and rename this file to subdirectory subset*:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZEP[i].Data(),gSystem->pwd());
      TString *subdir = new TString("subset");
      (*subdir)+=(i+1);
      (*subdir)+="/";
      (*subdir)+="mergedLYZEPanalysis";
      (*subdir)+=type.Data();
      (*subdir)+=".root";
      file->Move(subdir->Data());
      delete file;
      delete subdir;    
     } // end of for(Int_t i=0;i<nSubsets;i++)
    } else
      {
       TString mergedFileForPreviousCycleForSubsetLYZEP[nSubsets];
       TString mergedFileFinalForSubsetLYZEP[nSubsets];
       for(Int_t i=0;i<nSubsets;i++)
       {
        mergedFileForPreviousCycleForSubsetLYZEP[i]="mergedLYZEPanalysisForSubsetNo"; 
        mergedFileForPreviousCycleForSubsetLYZEP[i]+=(i+1);
        mergedFileForPreviousCycleForSubsetLYZEP[i]+=type.Data();
        mergedFileForPreviousCycleForSubsetLYZEP[i]+="CycleNo";
        mergedFileForPreviousCycleForSubsetLYZEP[i]+=((Int_t)fileCounterLastSubsetLYZEP/cycle);
        mergedFileForPreviousCycleForSubsetLYZEP[i]+=".root";
      
        fileMergerForSubsetsLYZEP[i].AddFile(mergedFileForPreviousCycleForSubsetLYZEP[i].Data());
    
        // delete merged output from previous cycle:
        TSystemFile *file = new TSystemFile(mergedFileForPreviousCycleForSubsetLYZEP[i].Data(),gSystem->pwd());
        file->Delete();
        delete file;
    
        mergedFileFinalForSubsetLYZEP[i]="subset";
        mergedFileFinalForSubsetLYZEP[i]+=(i+1);
        mergedFileFinalForSubsetLYZEP[i]+="/";
        mergedFileFinalForSubsetLYZEP[i]+="mergedLYZEPanalysis";
        mergedFileFinalForSubsetLYZEP[i]+=type.Data();
        mergedFileFinalForSubsetLYZEP[i]+=".root";
        
        fileMergerForSubsetsLYZEP[i].OutputFile(mergedFileFinalForSubsetLYZEP[i].Data());
        fileMergerForSubsetsLYZEP[i].Merge();
        fileMergerForSubsetsLYZEP[i].Reset();          
       } // end of for(Int_t i=0;i<nSubsets;i++)
      } // end of last else  
  } // end of if(devideOutputFilesInEqualSubsets) 
 
 delete dirNameSubset;
 delete baseDirPath;
 delete baseDir;
 
} // end of void mergeOutput(TString type="", Int_t mode=mLocal)

void LoadLibrariesMO(const libModes mode) {
  
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
  
} // end of void LoadLibrariesMO(const libModes mode)
