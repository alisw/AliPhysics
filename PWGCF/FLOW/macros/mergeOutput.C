// Macro mergeOutput.C is used to merge output files "AnalysisResults.root" locally. 
// (To merge output files "AnalysisResults.root" on Grid use macro mergeOutputOnGrid.C.)
// Its usage relies on the following directory structure and naming conventions:
//  1.) In directory <dir> there are subdirectories <subdir1>, <subdir2>, ...;   
//  2.) In subdirectory <subdir*> there is a file named "AnalysisResults.root";
//  3.) Copy macro mergeOutput.C in directory <dir> and launch it. It will automatically 
//      access from each subdirectory a file "AnalysisResults.root" and by making use of 
//      TFileMerger utility it will produce a merged file "mergedAnalysisResults.root" 
//      and save it in directory <dir>;
//  4.) IMPORTANT: Macro mergeOutput.C must be used in a combination with macro redoFinish.C
//      in order to get the analysis done on merged, large statistics sample. In particular,
//      after you got the file "mergedAnalysisResults.root", copy macro redoFinish.C in 
//      directory <dir> and launch it. This macro will access "mergedAnalysisResults.root"
//      and produce the output file "AnalysisResults.root" in directory <dir>. This file will 
//      hold the final results for merged, large statistics sample. 
//  5.) REMARK 1: To see plots for some of results use macro compareFlowResults.C. This macro
//      accesses file "AnalysisResults.root" and produces couple of predefined example plots.        
//  6.) REMARK 2: To use this macro for merging of other root files (having other name than
//      the default "AnalysisResults.root") modify the declaration of TString outputFileName.

// Name of the output files to be merged (sitting in <subdir1>, <subdir2>, ...):
TString outputFileName = "AnalysisResults.root";
// Name of the merged, large statistics file (to be saved in <dir>):
TString mergedFileName = "mergedAnalysisResults.root";
// Optionally set maximum number of files to be merged:
Int_t maxNoOfFiles = -1; // -1 = ignore
// For a large number of output files merging is done in cycles and this is the cycle period: 
const Int_t cycle = 50;

enum libModes {mLocal,mLocalSource};

void mergeOutput(Int_t mode=mLocal)
{
 // Mode: if mode = mLocal: analyze data on your computer using aliroot
 //       if mode = mLocalSource: analyze data on your computer using root + source files 
 
 // Cross-check user settings before starting:
 CrossCheckUserSettings(); 
 
 // Load needed flow libraries:
 LoadLibrariesMO(mode);  
 
 // Standard magic:
 TString *baseDirPath = new TString(gSystem->pwd());
 TSystemDirectory *baseDir = new TSystemDirectory(".",baseDirPath->Data());          
 TList *listOfFilesInBaseDir = baseDir->GetListOfFiles();
 TStopwatch timer;
 timer.Start();
 //listOfFilesInBaseDir->Print();
 
 // Count number of files to be merged:
 CountNumberOfFilesToBeMerged(baseDirPath,listOfFilesInBaseDir);
 
 // Loop over all files and from each subdirectory add file <outputFileName> to TFileMerger:
 Int_t nFiles = listOfFilesInBaseDir->GetEntries();
 Int_t fileCounter = 0;
 TFileMerger *fileMerger = new TFileMerger(); 
 cout<<" .... merging in progress, silence please ...."<<"\r"<<flush;
 for(Int_t iFile=0;iFile<nFiles;iFile++)
 {
  TSystemFile *currentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
  // Consider only subdirectories: 
  if(!currentFile || 
     !currentFile->IsDirectory() || 
     strcmp(currentFile->GetName(), ".") == 0 || 
     strcmp(currentFile->GetName(), "..") == 0){continue;} 
  // Accessing the output file <outputFileName> in current subdirectory: 
  TString currentSubDirName = baseDirPath->Data();
  (currentSubDirName+="/")+=currentFile->GetName();
  currentSubDirName+="/";
  TString fileName = currentSubDirName; 
  fileName+=outputFileName.Data();
  if(!(gSystem->AccessPathName(fileName.Data(),kFileExists)))
  {
   if(gROOT->GetVersionDate()>=20100404) // to be improved - removed eventually (this is temporary protection)
   { 
    Bool_t success = fileMerger->AddFile(fileName.Data(),kFALSE); // kFALSE switches off cp printout (not supported in older Root)
    if(success){fileCounter++;}
   } else 
     {
      Bool_t success = fileMerger->AddFile(fileName.Data());
      if(success){fileCounter++;}
     }    
   // Merging in cycles:
   if(fileCounter % cycle == 0)
   {
    // If the merged output from previous cycle exists add it to TFileMerger:
    TString *mergedFileForPreviousCycle = new TString("mergedCycle"); 
    (*mergedFileForPreviousCycle)+=(fileCounter/cycle - 1);
    (*mergedFileForPreviousCycle)+=".root";
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycle->Data(),kFileExists)))
    {
     if(gROOT->GetVersionDate()>=20100404) // to be improved - removed eventually (this is temporary protection)
     {  
      fileMerger->AddFile(mergedFileForPreviousCycle->Data(),kFALSE); // kFALSE switches off cp printout (not supported in older Root)
     } else
       {
        fileMerger->AddFile(mergedFileForPreviousCycle->Data());
       }
     // Delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycle->Data(),baseDirPath->Data());
     file->Delete();
     delete file;
    }    
    // Create merged output for the current cycle:
    TString *mergedFileForCurrentCycle = new TString("mergedCycle"); 
    (*mergedFileForCurrentCycle)+=(fileCounter/cycle);
    (*mergedFileForCurrentCycle)+=".root";    
    fileMerger->OutputFile(mergedFileForCurrentCycle->Data());
    fileMerger->Merge();
    cout<<" .... merged "<<fileCounter<<" files so far, marching on ....                       "<<"\r"<<flush;
    fileMerger->Reset();    
    delete mergedFileForPreviousCycle;
    delete mergedFileForCurrentCycle;
   } // end of if(fileCounter % cycle == 0) 
  } // end of if(!(gSystem->AccessPathName(fileName.Data(),kFileExists))) 
  if(fileCounter==maxNoOfFiles){break;}
 } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)
 
 //=================================================================================================
 
 // Final merging at the end of the day (3 distinct cases):
 gSystem->cd(baseDirPath->Data());
 if(fileCounter < cycle)
 {
  fileMerger->OutputFile(mergedFileName.Data());
  fileMerger->Merge();
 } else if(fileCounter % cycle == 0)
   {
    TString *mergedFileForPreviousCycle = new TString("mergedCycle"); 
    (*mergedFileForPreviousCycle)+=(fileCounter/cycle);
    (*mergedFileForPreviousCycle)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycle->Data(),kFileExists)))
    {
     gSystem->Rename(mergedFileForPreviousCycle->Data(),mergedFileName.Data());
    }     
   } else
     {
      TString *mergedFileForPreviousCycle = new TString("mergedCycle"); 
      (*mergedFileForPreviousCycle)+=((Int_t)fileCounter/cycle);
      (*mergedFileForPreviousCycle)+=".root";      
      if(gROOT->GetVersionDate()>=20100404) // to be improved - removed eventually (this is temporary protection)
      {  
       fileMerger->AddFile(mergedFileForPreviousCycle->Data(),kFALSE);
      } else
        {
         fileMerger->AddFile(mergedFileForPreviousCycle->Data());      
        }
      // Delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycle->Data(),baseDirPath->Data());
      file->Delete();
      delete file;
      fileMerger->OutputFile(mergedFileName.Data());
      fileMerger->Merge();
      delete mergedFileForPreviousCycle;
     }
 delete fileMerger;
 delete baseDirPath;
 delete baseDir;

 if(!(gSystem->AccessPathName(mergedFileName.Data(),kFileExists)) && fileCounter > 0)
 {
  cout<<" Merging went successfully: "<<fileCounter<<" files \""<<outputFileName.Data()<<"\" were"<<endl;
  cout<<" merged into the newly created file \""<<mergedFileName.Data()<<"\"."<<endl;
  cout<<endl;
  cout<<" Launch now macro redoFinish.C to get the correct final results."<<endl;
 } else
   {
    cout<<" WARNING: Merging failed miserably !!!!"<<endl;
   }   
 cout<<endl;
 timer.Stop();
 timer.Print(); 
 cout<<endl;

} // End of void mergeOutput(Int_t mode=mLocal)

// ===================================================================================

void CrossCheckUserSettings()
{
 // Cross-check user settings before starting:
 
 if(cycle>100)
 {
  cout<<endl;
  cout<<" WARNING: Cycle is too big !!!!"<<endl; 
  cout<<"          Set \"const Int_t cycle\" to smaller value in the macro."<<endl;
  cout<<endl;
  exit(0);
 }

} // void CrossCheckUserSettings() 
 
// ===================================================================================

void CountNumberOfFilesToBeMerged(TString *baseDirPath,TList *listOfFilesInBaseDir)
{
 // Count number of files to be merged and print info on the screen.
 
 if(!(baseDirPath && listOfFilesInBaseDir))
 {
  cout<<endl;
  cout<<" WARNING: baseDirPath||listOfFilesInBaseDir is NULL in CountNumberOfFilesToBeMerged(...) !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 Int_t nFiles = listOfFilesInBaseDir->GetEntries();
 Int_t fileCounter = 0;
 for(Int_t iFile=0;iFile<nFiles;iFile++)
 {
  TSystemFile *currentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
  // Consider only subdirectories: 
  if(!currentFile || 
     !currentFile->IsDirectory() || 
     strcmp(currentFile->GetName(), ".") == 0 || 
     strcmp(currentFile->GetName(), "..") == 0) continue; 
  // Accessing the output file <outputFileName> in current subdirectory: 
  TString currentSubDirName = baseDirPath->Data();
  (currentSubDirName+="/")+=currentFile->GetName();
  currentSubDirName+="/";
  TString fileName = currentSubDirName; 
  fileName+=outputFileName.Data();
  if(!(gSystem->AccessPathName(fileName.Data(),kFileExists)))
  {
   fileCounter++;
  }
 } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)

 // Print info on the screen:
 if(fileCounter>0)
 { 
  cout<<endl;
  cout<<" In subdirectores of directory "<<baseDirPath->Data()<<endl;
  cout<<" there are "<<fileCounter<<" files \""<<outputFileName.Data()<<"\" to be merged."<<endl;
  cout<<" Let the merging begins!"<<endl;
  cout<<endl;
 }
 else if(fileCounter==0)
 {
  cout<<endl;
  cout<<" TFileMerger wasn't lucky :'( "<<endl;
  cout<<" Couldn't find a single file \""<<outputFileName.Data()<<"\" to merge"<<endl;
  cout<<" in subdirectories of directory "<<baseDirPath->Data()<<" !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 
 return;
 
} // void CountNumberOfFilesToBeMerged(TString *dirName,TList *listOfAllFiles)

// ===================================================================================

void LoadLibrariesMO(const libModes mode) {
  
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
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWGflowBase");
  //cerr<<"libPWGflowBase loaded ..."<<endl;
  
  }
  
  else if (mode==mLocalSource) {
 
    // In root inline compile
  
    // Constants  
    gROOT->LoadMacro("BaseAliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("BaseAliFlowLYZConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("BaseAliFlowVector.cxx+"); 
    gROOT->LoadMacro("BaseAliFlowTrackSimple.cxx+");    
    gROOT->LoadMacro("BaseAliFlowTrackSimpleCuts.cxx+");    
    gROOT->LoadMacro("BaseAliFlowEventSimple.cxx+");
    
    // Output histosgrams
    gROOT->LoadMacro("BaseAliFlowCommonHist.cxx+");
    gROOT->LoadMacro("BaseAliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("BaseAliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("BaseAliFlowLYZHist2.cxx+");
       
    cout << "finished loading macros!" << endl;  
    
  } // end of else if (mode==mLocalSource) 
  
} // end of void LoadLibrariesMO(const libModes mode)

