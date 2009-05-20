enum libModes {mLocal,mLocalSource};
//mLocal: Analyze data on your computer using aliroot
//mLocalSource: Analyze data on your computer using root + source files

void mergeOutput(const Int_t nRuns=-1, TString type="",Int_t mode=mLocal)
{ 
 // load needed libraries:                       
 LoadSpreadLibraries(mode);  
  
 // file mergers for the output file of each method separately: 
 // MCEP:                       
 TFileMerger *mcepFileMerger = new TFileMerger();
 TString mergedFileNameMCEP("outputMCEPanalysis.root");
 mcepFileMerger->OutputFile(mergedFileNameMCEP);
 
 // SP:                       
 TFileMerger *spFileMerger = new TFileMerger();
 TString mergedFileNameSP("outputSPanalysis.root");
 spFileMerger->OutputFile(mergedFileNameSP);

 // GFC:                       
 TFileMerger *gfcFileMerger = new TFileMerger();
 TString mergedFileNameGFC("outputGFCanalysis.root");
 gfcFileMerger->OutputFile(mergedFileNameGFC);
 
 // QC:                       
 TFileMerger *qcFileMerger = new TFileMerger();
 TString mergedFileNameQC("outputQCanalysis.root");
 qcFileMerger->OutputFile(mergedFileNameQC);
 
 // FQD:                       
 TFileMerger *fqdFileMerger = new TFileMerger();
 TString mergedFileNameFQD("outputFQDanalysis.root");
 fqdFileMerger->OutputFile(mergedFileNameFQD);
 
 // LYZ1:                       
 TFileMerger *lyz1FileMerger = new TFileMerger();
 TString mergedFileNameLYZ1("outputLYZ1analysis.root");
 lyz1FileMerger->OutputFile(mergedFileNameLYZ1);
 
 // LYZ2:                       
 TFileMerger *lyz2FileMerger = new TFileMerger();
 TString mergedFileNameLYZ2("outputLYZ2analysis.root");
 lyz2FileMerger->OutputFile(mergedFileNameLYZ2);
 
 // LYZEP:                       
 TFileMerger *lyzepFileMerger = new TFileMerger();
 TString mergedFileNameLYZEP("outputLYZEPanalysis.root");
 lyzepFileMerger->OutputFile(mergedFileNameLYZEP);
 
 // standard magic:
 TString execDir(gSystem->pwd());  
 TSystemDirectory* baseDir = new TSystemDirectory(".",execDir.Data());          
 TList* dirList = baseDir->GetListOfFiles();
 Int_t nDirs = dirList->GetEntries();
 gSystem->cd(execDir);        

 Int_t counter = 0;
  
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
  fileNameMCEP+="outputMCEPanalysis.root";
  mcepFileMerger->AddFile(fileNameMCEP.Data());
  // SP:     
  TString fileNameSP = presentDirName;   
  fileNameSP+="outputSPanalysis.root";
  spFileMerger->AddFile(fileNameSP.Data());
  // GFC:     
  TString fileNameGFC = presentDirName;   
  fileNameGFC+="outputGFCanalysis.root";
  gfcFileMerger->AddFile(fileNameGFC.Data());
  // QC:     
  TString fileNameQC = presentDirName;   
  fileNameQC+="outputQCanalysis.root";
  qcFileMerger->AddFile(fileNameQC.Data());
  // FQD:     
  TString fileNameFQD = presentDirName;   
  fileNameFQD+="outputFQDanalysis.root";
  fqdFileMerger->AddFile(fileNameFQD.Data());
  // LYZ1:     
  TString fileNameLYZ1 = presentDirName;   
  fileNameLYZ1+="outputLYZ1analysis.root";
  lyz1FileMerger->AddFile(fileNameLYZ1.Data());
  // LYZ2:     
  TString fileNameLYZ2 = presentDirName;   
  fileNameLYZ2+="outputLYZ2analysis.root";
  lyz2FileMerger->AddFile(fileNameLYZ2.Data());
  // LYZEP:     
  TString fileNameLYZEP = presentDirName;   
  fileNameLYZEP+="outputLYZEPanalysis.root";
  lyzepFileMerger->AddFile(fileNameLYZEP.Data());
 
  counter++;
  
 } // end of for(Int_t iDir=0;iDir<nDirs;++iDir)
 
 // merge everything:
 mcepFileMerger->Merge();
 spFileMerger->Merge();
 gfcFileMerger->Merge();
 qcFileMerger->Merge();
 fqdFileMerger->Merge();
 lyz1FileMerger->Merge();
 lyz2FileMerger->Merge();
 lyzepFileMerger->Merge();
}


void LoadSpreadLibraries(const libModes mode) {
  
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
  cerr<<"libPWG2flowCommon.so loaded ..."<<endl;
  
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
    
  }  
  
}


