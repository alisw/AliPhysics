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
//  5.) REMARK: To see plots for some of the results use macro compareFlowResults.C. This macro
//      accesses file "AnalysisResults.root" and produces couple of predefined example plots.        

void mergeOutput()
{
 // Name of the output files to be merged:
 TString outputFileName = "AnalysisResults.root";
 // Name of the merged, large statistics file:
 TString mergedFileName = "mergedAnalysisResults.root";
 // Load needed flow libraries:
 gSystem->AddIncludePath("-I$ROOTSYS/include");
 gSystem->AddIncludePath("-I$ALICE_ROOT/include");
 gSystem->Load("libPWG2flowCommon");
 cerr<<"Library \"libPWG2flowCommon\" loaded ...."<<endl;
 // For a large number of output files merging is done in cycles
 // and this is the cycle period: 
 const Int_t cycle = 500;
 if(cycle>500)
 {
  cout<<"WARNING: Cycle is too big !!!!"<<endl; 
  cout<<"         Set \"const Int_t cycle\" to smaller value in the macro."<<endl;
  exit(0);
 }
 // Standard magic:
 TString *baseDirPath = new TString(gSystem->pwd());
 TSystemDirectory *baseDir = new TSystemDirectory(".",baseDirPath->Data());          
 TList *listOfFilesInBaseDir = baseDir->GetListOfFiles();
 TStopwatch timer;
 timer.Start();
 // listOfFilesInBaseDir->Print();
 Int_t nFiles = listOfFilesInBaseDir->GetEntries();
 // loop over all files and from each subdirectory add file AnalysisResults.root to TFileMerger:
 Int_t fileCounter = 0;
 TFileMerger *fileMerger = new TFileMerger(); 
 for(Int_t iFile=0;iFile<nFiles;iFile++)
 {
  TSystemFile *currentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
  // Consider only subdirectories: 
  if(!currentFile || 
     !currentFile->IsDirectory() || 
     strcmp(currentFile->GetName(), ".") == 0 || 
     strcmp(currentFile->GetName(), "..") == 0) continue; 
  // Accessing the output file "AnalysisResults.root" in current subdirectory: 
  TString currentSubDirName = baseDirPath->Data();
  (currentSubDirName+="/")+=currentFile->GetName();
  currentSubDirName+="/";
  TString fileName = currentSubDirName; 
  fileName+=outputFileName.Data();
  if(!(gSystem->AccessPathName(fileName.Data(),kFileExists)))
  {
   fileCounter++;
   fileMerger->AddFile(fileName.Data());
   // Merging in cycles:
   if(fileCounter % cycle == 0)
   {
    // If the merged output from previous cycle exists add it to TFileMerger:
    TString *mergedFileForPreviousCycle = new TString("mergedCycle"); 
    (*mergedFileForPreviousCycle)+=(fileCounter/cycle - 1);
    (*mergedFileForPreviousCycle)+=".root";
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycle->Data(),kFileExists)))
    {
     fileMerger->AddFile(mergedFileForPreviousCycle->Data());
     // Delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycle->Data(),baseDirPath->Data());
     file->Delete();
     delete file;
    }    
    // Create merged output for current cycle:
    TString *mergedFileForCurrentCycle = new TString("mergedCycle"); 
    (*mergedFileForCurrentCycle)+=(fileCounter/cycle);
    (*mergedFileForCurrentCycle)+=".root";    
    fileMerger->OutputFile(mergedFileForCurrentCycle->Data());
    fileMerger->Merge();
    fileMerger->Reset();    
    delete mergedFileForPreviousCycle;
    delete mergedFileForCurrentCycle;
   } // end of if(fileCounter % cycle == 0) 
  } // end of if(!(gSystem->AccessPathName(fileName.Data(),kFileExists))) 
 } // end of for(Int_t iFile=0;iFile<nFiles;iFile++)
 
 
 //=================================================================================================
 
 
 // Final merging at the end of the day (3 distinct cases):
 if(fileCounter==0)
 {
  cout<<endl;
  cout<<"Merger wasn't lucky: Couldn't find a single file "<<outputFileName.Data()<<" to merge"<<endl;
  cout<<"in subdirectories of directory "<<baseDirPath->Data()<<endl;
  cout<<endl;
  exit(0);
 }
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
      fileMerger->AddFile(mergedFileForPreviousCycle->Data());
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
 
 cout<<endl;
 timer.Stop();
 timer.Print(); 
 cout<<endl;
 if(!(gSystem->AccessPathName(mergedFileName.Data(),kFileExists)))
 {
  cout<<"Merging went successfully: "<<fileCounter<<" files "<<outputFileName.Data()<<" were"<<endl;
  cout<<"merged into the newly created file "<<mergedFileName.Data()<<"."<<endl;
  cout<<endl;
  cout<<"Launch now macro redoFinish.C to get the correct final results."<<endl;
 } else
   {
    cout<<"WARNING: Merging failed !!!!"<<endl;
   } 
 cout<<endl;
} // End of void mergeOutput()
