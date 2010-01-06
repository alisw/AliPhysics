// Macro mergeOutput.C is used to merge output files. Its usage relies
// on the following directory structure and naming conventions:
//  1.) In directory <dir> there are subdirectories <subdir1>, <subdir2>, ...;   
//  2.) In each subdirectory <subdir*> there is a file named "AnalysisResults.root";
//  3.) Copy macro mergeOutput.C in directory <dir> and launch it. It will 
//      automatically access from each subdirectory a file "AnalysisResults.root"
//      and by making use of TFileMerger utility it will produce a merged, large
//      statisics file "mergedAnalysisResults.root" and save it in directory <dir>;
//  4.) IMPORTANT: Macro mergeOutput.C must be used in a combination with macro redoFinish.C
//      in order to get the analysis done on merged, large statistics sample. In particular,
//      after you got the file "mergedAnalysisResults.root", copy macro redoFinish.C in 
//      directory <dir> and launch it. This macro will access "mergedAnalysisResults.root"
//      and produce the output file "AnalysisResults.root" in directory <dir>. This file will 
//      hold the final results for merged, large statistics sample. 
//  5.) REMARK: To see plots for some of the results use macro compareFlowResults.C. This macro
//      accesses file "AnalysisResults.root" and produces couple of predefined example plots.        

enum libModes {mLocal,mLocalSource};

void mergeOutput(Int_t mode=mLocal)
{
 // mode: if mode = mLocal: analyze data on your computer using aliroot
 //       if mode = mLocalSource: analyze data on your computer using root + source files

 TString outputFileName = "AnalysisResults.root"; // hardwired name of the output files to be merged
 
 const Int_t cycle = 10; // merging is done in cycles and this is the cycle period 

 if(cycle>100)
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

 // load needed libraries:                       
 LoadLibrariesMO(mode);  

 // standard magic:
 TString *baseDirPath = new TString(gSystem->pwd());
 TSystemDirectory *baseDir = new TSystemDirectory(".",baseDirPath->Data());          
 TList *listOfFilesInBaseDir = baseDir->GetListOfFiles();
 // listOfFilesInBaseDir->Print();
 Int_t nFiles = listOfFilesInBaseDir->GetEntries();
 
 // loop over all files and add files AnalysisResults.root to TFileMerger:
 gSystem->cd(baseDirPath->Data());
 Int_t fileCounter = 0;
 TFileMerger *fileMerger = new TFileMerger();
 
 for(Int_t iFile=0;iFile<nFiles;iFile++)
 {
  TSystemFile *presentFile = (TSystemFile*)listOfFilesInBaseDir->At(iFile);
 
  // consider only subdirectories: 
  if(!presentFile || !presentFile->IsDirectory() || strcmp(presentFile->GetName(), ".") == 0 || 
     strcmp(presentFile->GetName(), "..") == 0) continue; 
  
  TString presentDirName = baseDirPath->Data();
  (presentDirName+="/")+=presentFile->GetName();
  presentDirName+="/";
 
  // accessing the output AnalysisResults.root file in current subdirectory:
  TString fileName = presentDirName; 
  fileName+=outputFileName.Data();
  
  if(!(gSystem->AccessPathName(fileName.Data(),kFileExists)))
  {
   fileCounter++;
   fileMerger->AddFile(fileName.Data());
 
   // merging in cycles:
   if(fileCounter % cycle == 0)
   {
    TString *mergedFileForPreviousCycle = new TString("mergedAnalysisResults"); 
    (*mergedFileForPreviousCycle)+=(fileCounter/cycle - 1);
    (*mergedFileForPreviousCycle)+=".root";
    
    // if the merged output from previous cycle exists add it to TFileMerger:
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycle->Data(),kFileExists)))
    {
     fileMerger->AddFile(mergedFileForPreviousCycle->Data());
     // delete merged output from previous cycle:
     TSystemFile *file = new TSystemFile(mergedFileForPreviousCycle->Data(),gSystem->pwd());
     file->Delete();
     delete file;
    }
    
    // create merged output for current cycle:
    TString *mergedFileForCurrentCycle = new TString("mergedAnalysisResults"); 
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
 
 
 // final merging at the end of the day:
 gSystem->cd(baseDirPath->Data());

 if(fileCounter < cycle)
 {
  TString *mergedFileFinal = new TString("mergedAnalysisResults"); 
  (*mergedFileFinal)+=".root";
  fileMerger->OutputFile(mergedFileFinal->Data());
  fileMerger->Merge();
  delete mergedFileFinal;
 } else if (fileCounter % cycle == 0)
   {
    TString *mergedFileForPreviousCycle = new TString("mergedAnalysisResults"); 
    (*mergedFileForPreviousCycle)+=(fileCounter/cycle);
    (*mergedFileForPreviousCycle)+=".root";    
    if(!(gSystem->AccessPathName(mergedFileForPreviousCycle->Data(),kFileExists)))
    {
     TString *mergedFileFinal = new TString("mergedAnalysisResults"); 
     (*mergedFileFinal)+=".root";
     gSystem->Rename(mergedFileForPreviousCycle->Data(),mergedFileFinal->Data());
    }    
   } else
     {
      TString *mergedFileForPreviousCycle = new TString("mergedAnalysisResults"); 
      (*mergedFileForPreviousCycle)+=((Int_t)fileCounter/cycle);
      (*mergedFileForPreviousCycle)+=".root";
      
      fileMerger->AddFile(mergedFileForPreviousCycle->Data());
      
      // delete merged output from previous cycle:
      TSystemFile *file = new TSystemFile(mergedFileForPreviousCycle->Data(),gSystem->pwd());
      file->Delete();
      delete file;
      
      TString *mergedFileFinal = new TString("mergedAnalysisResults"); 
      (*mergedFileFinal)+=".root";
      fileMerger->OutputFile(mergedFileFinal->Data());
      fileMerger->Merge();
           
      delete mergedFileForPreviousCycle;
      delete mergedFileFinal;
     }
 
 delete fileMerger;
 delete baseDirPath;
 delete baseDir;

} // end of void mergeOutput(Int_t mode=mLocal)

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
