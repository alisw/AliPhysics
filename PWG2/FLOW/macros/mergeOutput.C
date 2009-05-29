enum libModes {mLocal,mLocalSource};

void mergeOutput(TString type="", Int_t nRuns=-1, Int_t mode=mLocal)
{
 // type:  type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //        (if type="" output files are from MC simulation (default))
 // nRuns: specify here how many output .root files will be merged 
 //        (if nRuns = -1 all of them will be merged)
 // mode:  if mode=mLocal analyze data on your computer using aliroot
 //        if mode=mLocalSource analyze data on your computer using root + source files
 
 Int_t cycle = 50; // if you cannot have >= 50 open files in memory, simply decrease this number 
 
 // load needed libraries:                       
 LoadSpreadLibrariesMO(mode);  
 
 // access the path of current directory:
 TString pwd(gSystem->pwd());
 pwd+="/";
    
 // file mergers for the output file of each method separately: 
 // MCEP:                       
 TFileMerger *mcepFileMerger = new TFileMerger();
 TString mergedFileNameMCEP("mergedMCEPanalysis");
 (mergedFileNameMCEP+=(type.Data()))+=(".root");
 TString pwdMCEP=pwd.Data();
 if(!(gSystem->AccessPathName((pwdMCEP+=mergedFileNameMCEP).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for MCEP !!!!"<<endl;
  break;
 }
 //mcepFileMerger->OutputFile(mergedFileNameMCEP);
 
 // SP:                       
 TFileMerger *spFileMerger = new TFileMerger();
 TString mergedFileNameSP("mergedSPanalysis");
 (mergedFileNameSP+=(type.Data()))+=(".root");
 TString pwdSP=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdSP+=mergedFileNameSP).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for SP !!!!"<<endl;
  break;
 }
 spFileMerger->OutputFile(mergedFileNameSP);
 
 // GFC:                       
 TFileMerger *gfcFileMerger = new TFileMerger();
 TString mergedFileNameGFC("mergedGFCanalysis");
 (mergedFileNameGFC+=(type.Data()))+=(".root");
 TString pwdGFC=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdGFC+=mergedFileNameGFC).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for GFC !!!!"<<endl;
  break;
 }
 gfcFileMerger->OutputFile(mergedFileNameGFC);
 
 // QC:                       
 TFileMerger *qcFileMerger = new TFileMerger();
 TString mergedFileNameQC("mergedQCanalysis");
 (mergedFileNameQC+=(type.Data()))+=(".root");
 TString pwdQC=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdQC+=mergedFileNameQC).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for QC !!!!"<<endl;
  break;
 }
 qcFileMerger->OutputFile(mergedFileNameQC);
 
 // FQD:                       
 TFileMerger *fqdFileMerger = new TFileMerger();
 TString mergedFileNameFQD("mergedFQDanalysis");
 (mergedFileNameFQD+=(type.Data()))+=(".root");
 TString pwdFQD=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdFQD+=mergedFileNameFQD).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for FQD !!!!"<<endl;
  break;
 }
 fqdFileMerger->OutputFile(mergedFileNameFQD);
 
 // LYZ1:                       
 TFileMerger *lyz1FileMerger = new TFileMerger();
 TString mergedFileNameLYZ1("mergedLYZ1analysis");
 (mergedFileNameLYZ1+=(type.Data()))+=(".root");
 TString pwdLYZ1=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdLYZ1+=mergedFileNameLYZ1).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for LYZ1 !!!!"<<endl;
  break;
 }
 lyz1FileMerger->OutputFile(mergedFileNameLYZ1);
 
 // LYZ2:                       
 TFileMerger *lyz2FileMerger = new TFileMerger();
 TString mergedFileNameLYZ2("mergedLYZ2analysis");
 (mergedFileNameLYZ2+=(type.Data()))+=(".root");
 TString pwdLYZ2=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdLYZ2+=mergedFileNameLYZ2).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for LYZ2 !!!!"<<endl;
  break;
 }
 lyz2FileMerger->OutputFile(mergedFileNameLYZ2);
 
 // LYZEP:                       
 TFileMerger *lyzepFileMerger = new TFileMerger();
 TString mergedFileNameLYZEP("mergedLYZEPanalysis");
 (mergedFileNameLYZEP+=(type.Data()))+=(".root");
 TString pwdLYZEP=pwd.Data(); 
 if(!(gSystem->AccessPathName((pwdLYZEP+=mergedFileNameLYZEP).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for LYZEP !!!!"<<endl;
  break;
 }
 lyzepFileMerger->OutputFile(mergedFileNameLYZEP);
 
 // standard magic:
 TString execDir(gSystem->pwd());  
 TSystemDirectory* baseDir = new TSystemDirectory(".",execDir.Data());          
 TList* dirList = baseDir->GetListOfFiles();
 Int_t nDirs = dirList->GetEntries();
 gSystem->cd(execDir);          

 Int_t counter = 0;
 
 Int_t counterMCEP = 0;
 Int_t maxCycleMCEP = 0;
 Int_t counterSP = 0;
 Int_t maxCycleSP = 0;
 Int_t counterGFC = 0;
 Int_t maxCycleGFC = 0;
 Int_t counterQC = 0;
 Int_t maxCycleQC = 0;
 Int_t counterFQD = 0;
 Int_t maxCycleFQD = 0;
 Int_t counterLYZ1 = 0;
 Int_t maxCycleLYZ1 = 0;
  
 if(cycle>nDirs) // to be improved (temporary workaroud)
 {
  cout<<"WARNING: Cycle "<<cycle<<" is too big. Decrease it's value in the declaration Int_t cycle in the macro !!!!"<<endl;
  break;
 }
   
 for(Int_t iDir=0;iDir<nDirs;++iDir)
 {
  TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir);
  if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || 
     strcmp(presentDir->GetName(), "..") == 0) continue; 
   
  if(nRuns != -1)
  {               
   if (counter >= nRuns) break;       
  } 
                       
  TString presentDirName(gSystem->pwd()); 
  presentDirName += "/";
  presentDirName += presentDir->GetName();
  presentDirName += "/";
   
  // accessing the output .root files from independent analysis for each method:
  // MCEP:     
  TString fileNameMCEP = presentDirName;   
  ((fileNameMCEP+="outputMCEPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameMCEP.Data(),kFileExists)))
  {
   mcepFileMerger->AddFile(fileNameMCEP.Data());
   counterMCEP++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameMCEP.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  if(counterMCEP % cycle == 0)
  {
   maxCycleMCEP = counterMCEP/cycle; 
   TString mergedFileNameCurrentCycleMCEP("mergedMCEPanalysis");   
   ((mergedFileNameCurrentCycleMCEP+=(type.Data()))+=counterMCEP/cycle)+=(".root");    
   mcepFileMerger->OutputFile(mergedFileNameCurrentCycleMCEP);
   
   TString mergedFileNamePreviousCycleMCEP("mergedMCEPanalysis");   
   ((mergedFileNamePreviousCycleMCEP+=(type.Data()))+=(counterMCEP/cycle-1))+=(".root");
   if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleMCEP.Data(),kFileExists)))
   {
    mcepFileMerger->AddFile(mergedFileNamePreviousCycleMCEP.Data());
   }
   
   Bool_t mcepMerged = kFALSE;
   if(mcepFileMerger)
   {  
    mcepMerged = mcepFileMerger->Merge();
    mcepFileMerger->Reset();
   }
   
   if(mcepMerged)
   {   
    TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleMCEP.Data(),gSystem->pwd());
    if(previousCycle) previousCycle->Delete();
    delete previousCycle;
   }
  } // end of if(counterMCEP % cycle == 0)

  // SP:     
  TString fileNameSP = presentDirName;   
  ((fileNameSP+="outputSPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameSP.Data(),kFileExists)))
  {
   spFileMerger->AddFile(fileNameSP.Data());
   counterSP++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameSP.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  if(counterSP % cycle == 0)
  {
   maxCycleSP = counterSP/cycle; 
   TString mergedFileNameCurrentCycleSP("mergedSPanalysis");   
   ((mergedFileNameCurrentCycleSP+=(type.Data()))+=counterSP/cycle)+=(".root");    
   spFileMerger->OutputFile(mergedFileNameCurrentCycleSP);
   
   TString mergedFileNamePreviousCycleSP("mergedSPanalysis");   
   ((mergedFileNamePreviousCycleSP+=(type.Data()))+=(counterSP/cycle-1))+=(".root");
   if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleSP.Data(),kFileExists)))
   {
    spFileMerger->AddFile(mergedFileNamePreviousCycleSP.Data());
   }
   
   Bool_t spMerged = kFALSE;
   if(spFileMerger)
   {  
    spMerged = spFileMerger->Merge();
    spFileMerger->Reset();
   }
   
   if(spMerged)
   {   
    TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleSP.Data(),gSystem->pwd());
    if(previousCycle) previousCycle->Delete();
    delete previousCycle;
   }
  } // end of if(counterSP % cycle == 0)
  
  // GFC:     
  TString fileNameGFC = presentDirName;   
  ((fileNameGFC+="outputGFCanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameGFC.Data(),kFileExists)))
  {
   gfcFileMerger->AddFile(fileNameGFC.Data());
   counterGFC++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameGFC.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  if(counterGFC % cycle == 0)
  {
   maxCycleGFC = counterGFC/cycle; 
   TString mergedFileNameCurrentCycleGFC("mergedGFCanalysis");   
   ((mergedFileNameCurrentCycleGFC+=(type.Data()))+=counterGFC/cycle)+=(".root");    
   gfcFileMerger->OutputFile(mergedFileNameCurrentCycleGFC);
   
   TString mergedFileNamePreviousCycleGFC("mergedGFCanalysis");   
   ((mergedFileNamePreviousCycleGFC+=(type.Data()))+=(counterGFC/cycle-1))+=(".root");
   if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleGFC.Data(),kFileExists)))
   {
    gfcFileMerger->AddFile(mergedFileNamePreviousCycleGFC.Data());
   }
   
   Bool_t gfcMerged = kFALSE;
   if(gfcFileMerger)
   {  
    gfcMerged = gfcFileMerger->Merge();
    gfcFileMerger->Reset();
   }
   
   if(gfcMerged)
   {   
    TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleGFC.Data(),gSystem->pwd());
    if(previousCycle) previousCycle->Delete();
    delete previousCycle;
   }
  } // end of if(counterGFC % cycle == 0)
  
  // QC:     
  TString fileNameQC = presentDirName;   
  ((fileNameQC+="outputQCanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameQC.Data(),kFileExists)))
  {
   qcFileMerger->AddFile(fileNameQC.Data());
   counterQC++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameQC.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  if(counterQC % cycle == 0)
  {
   maxCycleQC = counterQC/cycle; 
   TString mergedFileNameCurrentCycleQC("mergedQCanalysis");   
   ((mergedFileNameCurrentCycleQC+=(type.Data()))+=counterQC/cycle)+=(".root");    
   qcFileMerger->OutputFile(mergedFileNameCurrentCycleQC);
   
   TString mergedFileNamePreviousCycleQC("mergedQCanalysis");   
   ((mergedFileNamePreviousCycleQC+=(type.Data()))+=(counterQC/cycle-1))+=(".root");
   if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleQC.Data(),kFileExists)))
   {
    qcFileMerger->AddFile(mergedFileNamePreviousCycleQC.Data());
   }
   
   Bool_t qcMerged = kFALSE;
   if(qcFileMerger)
   {  
    qcMerged = qcFileMerger->Merge();
    qcFileMerger->Reset();
   }
   
   if(qcMerged)
   {   
    TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleQC.Data(),gSystem->pwd());
    if(previousCycle) previousCycle->Delete();
    delete previousCycle;
   }
  } // end of if(counterQC % cycle == 0)
  
  // FQD:     
  TString fileNameFQD = presentDirName;   
  ((fileNameFQD+="outputFQDanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameFQD.Data(),kFileExists)))
  {
   fqdFileMerger->AddFile(fileNameFQD.Data());
   counterFQD++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameFQD.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  if(counterFQD % cycle == 0)
  {
   maxCycleFQD = counterFQD/cycle; 
   TString mergedFileNameCurrentCycleFQD("mergedFQDanalysis");   
   ((mergedFileNameCurrentCycleFQD+=(type.Data()))+=counterFQD/cycle)+=(".root");    
   fqdFileMerger->OutputFile(mergedFileNameCurrentCycleFQD);
   
   TString mergedFileNamePreviousCycleFQD("mergedFQDanalysis");   
   ((mergedFileNamePreviousCycleFQD+=(type.Data()))+=(counterFQD/cycle-1))+=(".root");
   if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleFQD.Data(),kFileExists)))
   {
    fqdFileMerger->AddFile(mergedFileNamePreviousCycleFQD.Data());
   }
   
   Bool_t fqdMerged = kFALSE;
   if(fqdFileMerger)
   {  
    fqdMerged = fqdFileMerger->Merge();
    fqdFileMerger->Reset();
   }
   
   if(fqdMerged)
   {   
    TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleFQD.Data(),gSystem->pwd());
    if(previousCycle) previousCycle->Delete();
    delete previousCycle;
   }
  } // end of if(counterFQD % cycle == 0)
  
  // LYZ1:     
  TString fileNameLYZ1 = presentDirName;   
  ((fileNameLYZ1+="outputLYZ1analysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZ1.Data(),kFileExists)))
  {
   lyz1FileMerger->AddFile(fileNameLYZ1.Data());
   counterLYZ1++;
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameLYZ1.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  if(counterLYZ1 % cycle == 0)
  {
   maxCycleLYZ1 = counterLYZ1/cycle; 
   TString mergedFileNameCurrentCycleLYZ1("mergedLYZ1analysis");   
   ((mergedFileNameCurrentCycleLYZ1+=(type.Data()))+=counterLYZ1/cycle)+=(".root");    
   lyz1FileMerger->OutputFile(mergedFileNameCurrentCycleLYZ1);
   
   TString mergedFileNamePreviousCycleLYZ1("mergedLYZ1analysis");   
   ((mergedFileNamePreviousCycleLYZ1+=(type.Data()))+=(counterLYZ1/cycle-1))+=(".root");
   if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleLYZ1.Data(),kFileExists)))
   {
    lyz1FileMerger->AddFile(mergedFileNamePreviousCycleLYZ1.Data());
   }
   
   Bool_t lyz1Merged = kFALSE;
   if(lyz1FileMerger)
   {  
    lyz1Merged = lyz1FileMerger->Merge();
    lyz1FileMerger->Reset();
   }
   
   if(lyz1Merged)
   {   
    TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleLYZ1.Data(),gSystem->pwd());
    if(previousCycle) previousCycle->Delete();
    delete previousCycle;
   }
  } // end of if(counterLYZ1 % cycle == 0)

  counter++;
  
 } // end of for(Int_t iDir=0;iDir<nDirs;++iDir)
 
 // last cycle MCEP:
 if(counterMCEP % cycle != 0)
 { 
  TString mergedFileNameCurrentCycleMCEP("mergedMCEPanalysis");   
  ((mergedFileNameCurrentCycleMCEP+=(type.Data()))+=(maxCycleMCEP+1))+=(".root");    
  mcepFileMerger->OutputFile(mergedFileNameCurrentCycleMCEP);
   
  TString mergedFileNamePreviousCycleMCEP("mergedMCEPanalysis");   
  ((mergedFileNamePreviousCycleMCEP+=(type.Data()))+=(maxCycleMCEP))+=(".root");
  if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleMCEP.Data(),kFileExists)))
  {
   mcepFileMerger->AddFile(mergedFileNamePreviousCycleMCEP.Data());
  }
   
  Bool_t mcepMerged = kFALSE;
  if(mcepFileMerger)
  {  
   mcepMerged = mcepFileMerger->Merge();
   mcepFileMerger->Reset();
  }
   
  if(mcepMerged)
  {   
   TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleMCEP.Data(),gSystem->pwd());
   if(previousCycle) previousCycle->Delete();
   delete previousCycle;
  }
  maxCycleMCEP++;
 } // end of if(counterMCEP % cycle != 0)
  
 // renaming the final merged output for MCEP: 
 TString finalMergedFileNameMCEP("mergedMCEPanalysis");   
 ((finalMergedFileNameMCEP+=(type.Data()))+=maxCycleMCEP)+=(".root"); 
 TSystemFile *finalMergedFileMCEP = new TSystemFile(finalMergedFileNameMCEP.Data(),gSystem->pwd());
 TString defaultMergedFileNameMCEP("mergedMCEPanalysis");   
 (defaultMergedFileNameMCEP+=(type.Data()))+=(".root"); 
 if(finalMergedFileMCEP) finalMergedFileMCEP->Rename(defaultMergedFileNameMCEP.Data());
 
 // last cycle SP:
 if(counterSP % cycle != 0)
 { 
  TString mergedFileNameCurrentCycleSP("mergedSPanalysis");   
  ((mergedFileNameCurrentCycleSP+=(type.Data()))+=(maxCycleSP+1))+=(".root");    
  spFileMerger->OutputFile(mergedFileNameCurrentCycleSP);
   
  TString mergedFileNamePreviousCycleSP("mergedSPanalysis");   
  ((mergedFileNamePreviousCycleSP+=(type.Data()))+=(maxCycleSP))+=(".root");
  if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleSP.Data(),kFileExists)))
  {
   spFileMerger->AddFile(mergedFileNamePreviousCycleSP.Data());
  }
   
  Bool_t spMerged = kFALSE;
  if(spFileMerger)
  {  
   spMerged = spFileMerger->Merge();
   spFileMerger->Reset();
  }
   
  if(spMerged)
  {   
   TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleSP.Data(),gSystem->pwd());
   if(previousCycle) previousCycle->Delete();
   delete previousCycle;
  }
  maxCycleSP++;
 } // end of if(counterSP % cycle != 0)
 
 // renaming the final merged output for SP: 
 TString finalMergedFileNameSP("mergedSPanalysis");   
 ((finalMergedFileNameSP+=(type.Data()))+=maxCycleSP)+=(".root"); 
 TSystemFile *finalMergedFileSP = new TSystemFile(finalMergedFileNameSP.Data(),gSystem->pwd());
 TString defaultMergedFileNameSP("mergedSPanalysis");   
 (defaultMergedFileNameSP+=(type.Data()))+=(".root"); 
 if(finalMergedFileSP) finalMergedFileSP->Rename(defaultMergedFileNameSP.Data());
 
 // last cycle GFC:
 if(counterGFC % cycle != 0)
 { 
  TString mergedFileNameCurrentCycleGFC("mergedGFCanalysis");   
  ((mergedFileNameCurrentCycleGFC+=(type.Data()))+=(maxCycleGFC+1))+=(".root");    
  gfcFileMerger->OutputFile(mergedFileNameCurrentCycleGFC);
   
  TString mergedFileNamePreviousCycleGFC("mergedGFCanalysis");   
  ((mergedFileNamePreviousCycleGFC+=(type.Data()))+=(maxCycleGFC))+=(".root");
  if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleGFC.Data(),kFileExists)))
  {
   gfcFileMerger->AddFile(mergedFileNamePreviousCycleGFC.Data());
  }
   
  Bool_t gfcMerged = kFALSE;
  if(gfcFileMerger)
  {  
   gfcMerged = gfcFileMerger->Merge();
   gfcFileMerger->Reset();
  }
   
  if(gfcMerged)
  {   
   TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleGFC.Data(),gSystem->pwd());
   if(previousCycle) previousCycle->Delete();
   delete previousCycle;
  }
  maxCycleGFC++;
 } // end of if(counterGFC % cycle != 0)
 
 // renaming the final merged output for GFC: 
 TString finalMergedFileNameGFC("mergedGFCanalysis");   
 ((finalMergedFileNameGFC+=(type.Data()))+=maxCycleGFC)+=(".root"); 
 TSystemFile *finalMergedFileGFC = new TSystemFile(finalMergedFileNameGFC.Data(),gSystem->pwd());
 TString defaultMergedFileNameGFC("mergedGFCanalysis");   
 (defaultMergedFileNameGFC+=(type.Data()))+=(".root"); 
 if(finalMergedFileGFC) finalMergedFileGFC->Rename(defaultMergedFileNameGFC.Data());
 
 // last cycle QC:
 if(counterQC % cycle != 0)
 { 
  TString mergedFileNameCurrentCycleQC("mergedQCanalysis");   
  ((mergedFileNameCurrentCycleQC+=(type.Data()))+=(maxCycleQC+1))+=(".root");    
  qcFileMerger->OutputFile(mergedFileNameCurrentCycleQC);
   
  TString mergedFileNamePreviousCycleQC("mergedQCanalysis");   
  ((mergedFileNamePreviousCycleQC+=(type.Data()))+=(maxCycleQC))+=(".root");
  if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleQC.Data(),kFileExists)))
  {
   qcFileMerger->AddFile(mergedFileNamePreviousCycleQC.Data());
  }
   
  Bool_t qcMerged = kFALSE;
  if(qcFileMerger)
  {  
   qcMerged = qcFileMerger->Merge();
   qcFileMerger->Reset();
  }
   
  if(qcMerged)
  {   
   TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleQC.Data(),gSystem->pwd());
   if(previousCycle) previousCycle->Delete();
   delete previousCycle;
  }
  maxCycleQC++;
 } // end of if(counterQC % cycle != 0)
 
 // renaming the final merged output for QC: 
 TString finalMergedFileNameQC("mergedQCanalysis");   
 ((finalMergedFileNameQC+=(type.Data()))+=maxCycleQC)+=(".root"); 
 TSystemFile *finalMergedFileQC = new TSystemFile(finalMergedFileNameQC.Data(),gSystem->pwd());
 TString defaultMergedFileNameQC("mergedQCanalysis");   
 (defaultMergedFileNameQC+=(type.Data()))+=(".root"); 
 if(finalMergedFileQC) finalMergedFileQC->Rename(defaultMergedFileNameQC.Data());
 
 // last cycle FQD:
 if(counterFQD % cycle != 0)
 { 
  TString mergedFileNameCurrentCycleFQD("mergedFQDanalysis");   
  ((mergedFileNameCurrentCycleFQD+=(type.Data()))+=(maxCycleFQD+1))+=(".root");    
  fqdFileMerger->OutputFile(mergedFileNameCurrentCycleFQD);
   
  TString mergedFileNamePreviousCycleFQD("mergedFQDanalysis");   
  ((mergedFileNamePreviousCycleFQD+=(type.Data()))+=(maxCycleFQD))+=(".root");
  if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleFQD.Data(),kFileExists)))
  {
   fqdFileMerger->AddFile(mergedFileNamePreviousCycleFQD.Data());
  }
   
  Bool_t fqdMerged = kFALSE;
  if(fqdFileMerger)
  {  
   fqdMerged = fqdFileMerger->Merge();
   fqdFileMerger->Reset();
  }
   
  if(fqdMerged)
  {   
   TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleFQD.Data(),gSystem->pwd());
   if(previousCycle) previousCycle->Delete();
   delete previousCycle;
  }
  maxCycleFQD++;
 } // end of if(counterFQD % cycle != 0)
 
 // renaming the final merged output for FQD: 
 TString finalMergedFileNameFQD("mergedFQDanalysis");   
 ((finalMergedFileNameFQD+=(type.Data()))+=maxCycleFQD)+=(".root"); 
 TSystemFile *finalMergedFileFQD = new TSystemFile(finalMergedFileNameFQD.Data(),gSystem->pwd());
 TString defaultMergedFileNameFQD("mergedFQDanalysis");   
 (defaultMergedFileNameFQD+=(type.Data()))+=(".root"); 
 if(finalMergedFileFQD) finalMergedFileFQD->Rename(defaultMergedFileNameFQD.Data());
 
 // last cycle LYZ1:
 if(counterLYZ1 % cycle != 0)
 { 
  TString mergedFileNameCurrentCycleLYZ1("mergedLYZ1analysis");   
  ((mergedFileNameCurrentCycleLYZ1+=(type.Data()))+=(maxCycleLYZ1+1))+=(".root");    
  lyz1FileMerger->OutputFile(mergedFileNameCurrentCycleLYZ1);
   
  TString mergedFileNamePreviousCycleLYZ1("mergedLYZ1analysis");   
  ((mergedFileNamePreviousCycleLYZ1+=(type.Data()))+=(maxCycleLYZ1))+=(".root");
  if(!(gSystem->AccessPathName(mergedFileNamePreviousCycleLYZ1.Data(),kFileExists)))
  {
   lyz1FileMerger->AddFile(mergedFileNamePreviousCycleLYZ1.Data());
  }
   
  Bool_t lyz1Merged = kFALSE;
  if(lyz1FileMerger)
  {  
   lyz1Merged = lyz1FileMerger->Merge();
   lyz1FileMerger->Reset();
  }
   
  if(lyz1Merged)
  {   
   TSystemFile *previousCycle = new TSystemFile(mergedFileNamePreviousCycleLYZ1.Data(),gSystem->pwd());
   if(previousCycle) previousCycle->Delete();
   delete previousCycle;
  }
  maxCycleLYZ1++;
 } // end of if(counterLYZ1 % cycle != 0)
 
 // renaming the final merged output for LYZ1: 
 TString finalMergedFileNameLYZ1("mergedLYZ1analysis");   
 ((finalMergedFileNameLYZ1+=(type.Data()))+=maxCycleLYZ1)+=(".root"); 
 TSystemFile *finalMergedFileLYZ1 = new TSystemFile(finalMergedFileNameLYZ1.Data(),gSystem->pwd());
 TString defaultMergedFileNameLYZ1("mergedLYZ1analysis");   
 (defaultMergedFileNameLYZ1+=(type.Data()))+=(".root"); 
 if(finalMergedFileLYZ1) finalMergedFileLYZ1->Rename(defaultMergedFileNameLYZ1.Data());
 
} // end of void mergeOutput(TString type="", const Int_t nRuns=-1, Int_t mode=mLocal)

void LoadSpreadLibrariesMO(const libModes mode) {
  
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
  //cerr<<"libPWG2flowCommon.so loaded ..."<<endl;
  
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
  
} // end of void LoadSpreadLibrariesMO(const libModes mode)


