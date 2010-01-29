// Macro mergeOutputOnGrid.C is used to merge output files "AnalysisResults.root" 
// from the AliEn file catalog. It produces the merged, larged statistics file
// "mergedAnalysisResults.root". To get the final flow estimates for the large
// statistics file use macro redoFinish.C (see comments in macro mergeOutput.C
// for more details). 

void mergeOutputOnGrid(const char* gridPath = "/alice/cern.ch/user/a/abilandz/flowAnalysisOnGrid/output*") 
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
 // Connect to the GRID:
 TGrid::Connect("alien://");  
 // Query the file catalog and get a TGridResult:
 TString queryPattern = outputFileName;
 TGridResult *result = gGrid->Query(gridPath,queryPattern.Data());
 Int_t nEntries = result->GetEntries();
 Printf("\nFound %d files %s in %s ....\n",nEntries,outputFileName.Data(),gridPath);   
 // Create a TFileMerger:
 TFileMerger *fileMerger = new TFileMerger(); 
 // For a large number of output files merging is done in cycles
 // and this is the cycle period: 
 const Int_t cycle = 500;
 if(cycle>500)
 {
  cout<<"WARNING: Cycle is too big !!!!"<<endl; 
  cout<<"         Set \"const Int_t cycle\" to smaller value in the macro."<<endl;
  exit(0);
 }
 Int_t fileCounter = 0;
 TString *baseDirPath = new TString(gSystem->pwd());
 // Loop over the TGridResult's entries and add the files to the TFileMerger:
 TString alienUrl;
 for(Int_t i=0;i<nEntries;i++) 
 {
  alienUrl = result->GetKey(i,"turl");
  //Printf("%s",alienUrl.Data());
  fileMerger->AddFile(alienUrl.Data());
  fileCounter++;
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
 } // end of for(Int_t i=0;i<nEntries;i++) 

 
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
}
