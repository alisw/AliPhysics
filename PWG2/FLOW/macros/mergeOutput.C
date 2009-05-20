enum libModes {mLocal,mLocalSource};
//mLocal: Analyze data on your computer using aliroot
//mLocalSource: Analyze data on your computer using root + source files

void mergeOutput(const Int_t nRuns=-1, TString type="",Int_t mode=mLocal)
{ 
 // load needed libraries:                       
 LoadSpreadLibraries(mode);  

 TString pwd(gSystem->pwd());
         
 // file mergers for the output file of each method separately: 
 // MCEP:                       
 TFileMerger *mcepFileMerger = new TFileMerger();
 TString mergedFileNameMCEP("outputMCEPanalysis.root");
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameMCEP).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for MCEP !!!!"<<endl;
  break;
 }
 mcepFileMerger->OutputFile(mergedFileNameMCEP);
 
 // SP:                       
 TFileMerger *spFileMerger = new TFileMerger();
 TString mergedFileNameSP("outputSPanalysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameSP).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for SP !!!!"<<endl;
  break;
 }
 spFileMerger->OutputFile(mergedFileNameSP);

 // GFC:                       
 TFileMerger *gfcFileMerger = new TFileMerger();
 TString mergedFileNameGFC("outputGFCanalysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameGFC).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for GFC !!!!"<<endl;
  break;
 }
 gfcFileMerger->OutputFile(mergedFileNameGFC);
 
 // QC:                       
 TFileMerger *qcFileMerger = new TFileMerger();
 TString mergedFileNameQC("outputQCanalysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameQC).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for QC !!!!"<<endl;
  break;
 } 
 qcFileMerger->OutputFile(mergedFileNameQC);
 
 // FQD:                       
 TFileMerger *fqdFileMerger = new TFileMerger();
 TString mergedFileNameFQD("outputFQDanalysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameFQD).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for FQD !!!!"<<endl;
  break;
 } 
 fqdFileMerger->OutputFile(mergedFileNameFQD);
 
 // LYZ1:                       
 TFileMerger *lyz1FileMerger = new TFileMerger();
 TString mergedFileNameLYZ1("outputLYZ1analysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameLYZ1).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for LYZ1 !!!!"<<endl;
  break;
 } 
 lyz1FileMerger->OutputFile(mergedFileNameLYZ1);
 
 // LYZ2:                       
 TFileMerger *lyz2FileMerger = new TFileMerger();
 TString mergedFileNameLYZ2("outputLYZ2analysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameLYZ2).Data(),kFileExists)))
 {
  cout<<"WARNING: You already have a merged output for LYZ2 !!!!"<<endl;
  break;
 }
 lyz2FileMerger->OutputFile(mergedFileNameLYZ2);
 
 // LYZEP:                       
 TFileMerger *lyzepFileMerger = new TFileMerger();
 TString mergedFileNameLYZEP("outputLYZEPanalysis.root");
 pwd=gSystem->pwd();
 if(!(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameLYZEP).Data(),kFileExists)))
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
  if(!(gSystem->AccessPathName(fileNameMCEP.Data(),kFileExists)))
  {
   mcepFileMerger->AddFile(fileNameMCEP.Data());
  }  
  // SP:     
  TString fileNameSP = presentDirName;   
  fileNameSP+="outputSPanalysis.root";
  if(!(gSystem->AccessPathName(fileNameSP.Data(),kFileExists)))
  {
   spFileMerger->AddFile(fileNameSP.Data());
  } 
  // GFC:     
  TString fileNameGFC = presentDirName;   
  fileNameGFC+="outputGFCanalysis.root";
  if(!(gSystem->AccessPathName(fileNameGFC.Data(),kFileExists)))
  { 
   gfcFileMerger->AddFile(fileNameGFC.Data());
  } 
  // QC:     
  TString fileNameQC = presentDirName;   
  fileNameQC+="outputQCanalysis.root";
  if(!(gSystem->AccessPathName(fileNameQC.Data(),kFileExists)))
  { 
   qcFileMerger->AddFile(fileNameQC.Data());
  } 
  // FQD:     
  TString fileNameFQD = presentDirName;   
  fileNameFQD+="outputFQDanalysis.root";
  if(!(gSystem->AccessPathName(fileNameFQD.Data(),kFileExists)))
  { 
   fqdFileMerger->AddFile(fileNameFQD.Data());
  } 
  // LYZ1:     
  TString fileNameLYZ1 = presentDirName;   
  fileNameLYZ1+="outputLYZ1analysis.root";
  if(!(gSystem->AccessPathName(fileNameLYZ1.Data(),kFileExists)))
  { 
   lyz1FileMerger->AddFile(fileNameLYZ1.Data());
  } 
  // LYZ2:     
  TString fileNameLYZ2 = presentDirName;   
  fileNameLYZ2+="outputLYZ2analysis.root";
  if(!(gSystem->AccessPathName(fileNameLYZ2.Data(),kFileExists)))
  {  
   lyz2FileMerger->AddFile(fileNameLYZ2.Data());
  } 
  // LYZEP:     
  TString fileNameLYZEP = presentDirName;   
  fileNameLYZEP+="outputLYZEPanalysis.root";
  if(!(gSystem->AccessPathName(fileNameLYZEP.Data(),kFileExists)))
  {  
   lyzepFileMerger->AddFile(fileNameLYZEP.Data());
  }
   
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


