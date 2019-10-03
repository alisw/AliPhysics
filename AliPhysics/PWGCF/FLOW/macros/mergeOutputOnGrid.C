// Macro mergeOutputOnGrid.C is used to merge output files "AnalysisResults.root" 
// from the AliEn file catalog. Remarks:
//  1.) It uses your local AliRoot version - try to have locally the same AliRoot
//      version as the one you have specified for the analysis on Grid, otherwise
//      strange things might happen ....  
//  2.) It produces the merged, large statistics file "mergedAnalysisResults.root" 
//      in the directory <dir> from which the macro was launched. Results stored 
//      in merged file are WRONG because after merging the results from small 
//      statistics files are trivially summed up in all histograms. To get the 
//      final correct flow estimates for the large statistics file copy in <dir> 
//      and execute macro redoFinish.C (see comments in redoFinish.C). 

// Name of the output files to be merged:
TString outputFileName = "AnalysisResults.root";
// Name of the merged, large statistics file:
TString mergedFileName = "mergedAnalysisResults.root";
// Optionally set maximum number of files to be merged:
Int_t maxNoOfFiles = -1; // -1 = ignore
// For a large number of output files merging is done in cycles and this is the cycle period: 
const Int_t cycle = 25;

void mergeOutputOnGrid(const char* gridPath = "/alice/cern.ch/user/a/abilandz/sim/LHC10e16/default/output/*") 
{
 // Cross-check user settings before starting:
 CrossCheckUserSettings(); 
 // Time: 
 TStopwatch timer;
 timer.Start();
 // Load needed flow libraries:
 gSystem->AddIncludePath("-I$ROOTSYS/include");
 gSystem->AddIncludePath("-I$ALICE_ROOT/include");
 gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
 gSystem->Load("libPWGflowBase");
 //cerr<<"Library \"libPWGflowBase\" loaded ...."<<endl;
 // Connect to the GRID:
 TGrid::Connect("alien://");  
 // Query the file catalog and get a TGridResult:
 TString queryPattern = outputFileName;
 TGridResult *result = gGrid->Query(gridPath,queryPattern.Data());
 Int_t nEntries = result->GetEntries();
 Printf("\n Found %d files \"%s\" in subdirectories of directory \n %s \n Let the merging begins! \n",nEntries,outputFileName.Data(),gridPath);   
 gSystem->Sleep(4400); // in ms
 // Create a TFileMerger:
 TFileMerger *fileMerger = new TFileMerger(); 
 Int_t fileCounter = 0;
 TString *baseDirPath = new TString(gSystem->pwd());
 // Loop over the TGridResult's entries and add the files to the TFileMerger:
 TString alienUrl;
 for(Int_t i=0;i<nEntries;i++) 
 {
  alienUrl = result->GetKey(i,"turl");
  //Printf("%s",alienUrl.Data());  
  if(gROOT->GetVersionDate()>=20100404) // to be improved - removed eventually (this is temporary protection)
  { 
   Bool_t success = fileMerger->AddFile(alienUrl.Data(),kFALSE); // kFALSE switches off cp printout (not supported in older Root)
   if(success){fileCounter++;}
  } else 
    {
     Bool_t success = fileMerger->AddFile(alienUrl.Data());
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
   // Create merged output for current cycle:
   TString *mergedFileForCurrentCycle = new TString("mergedCycle"); 
   (*mergedFileForCurrentCycle)+=(fileCounter/cycle);
   (*mergedFileForCurrentCycle)+=".root";    
   fileMerger->OutputFile(mergedFileForCurrentCycle->Data());
   fileMerger->Merge();
   Int_t percentage = 0;
   if(maxNoOfFiles > 0)
   {
    percentage = (Int_t)100.*fileCounter/maxNoOfFiles;
   } else
     {
      percentage = (Int_t)100.*fileCounter/nEntries;   
     }   
   Int_t max = 0;
   if(maxNoOfFiles > 0)
   {
    max = maxNoOfFiles;
   } else
     {
      max = nEntries;
     }     
   if(max-fileCounter > 0)
   {  
    cout<<endl;
    cout<<" .... merged "<<percentage<<"% requested files so far, marching on ....                       "<<endl;   
    cout<<endl;
   }  
   fileMerger->Reset();    
   delete mergedFileForPreviousCycle;
   delete mergedFileForCurrentCycle;
  } // end of if(fileCounter % cycle == 0) 
  if(fileCounter==maxNoOfFiles){break;}
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
 
 if(!(gSystem->AccessPathName(mergedFileName.Data(),kFileExists)) && fileCounter > 0)
 {
  cout<<endl;
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

} // end of void mergeOutputOnGrid(const char* gridPath = "...")

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
