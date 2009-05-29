enum libModes {mLocal,mLocalSource};

void mergeOutput(TString type="", Int_t nRuns=-1, Int_t mode=mLocal)
{
 // type:  type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //        (if type="" output files are from MC simulation (default))
 // nRuns: specify here how many output .root files will be merged 
 //        (if nRuns = -1 all of them will be merged)
 // mode:  if mode=mLocal analyze data on your computer using aliroot
 //        if mode=mLocalSource analyze data on your computer using root + source files
 
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
 mcepFileMerger->OutputFile(mergedFileNameMCEP);
 
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
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameMCEP.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  // SP:     
  TString fileNameSP = presentDirName;   
  ((fileNameSP+="outputSPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameSP.Data(),kFileExists)))
  {
   spFileMerger->AddFile(fileNameSP.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameSP.Data()<<". Merging will continue without this file."<<endl;
    }  
    
  // GFC:     
  TString fileNameGFC = presentDirName;   
  ((fileNameGFC+="outputGFCanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameGFC.Data(),kFileExists)))
  {
   gfcFileMerger->AddFile(fileNameGFC.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameGFC.Data()<<". Merging will continue without this file."<<endl;
    }    
        
  // QC:     
  TString fileNameQC = presentDirName;   
  ((fileNameQC+="outputQCanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameQC.Data(),kFileExists)))
  {
   qcFileMerger->AddFile(fileNameQC.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameQC.Data()<<". Merging will continue without this file."<<endl;
    }   
    
  // FQD:     
  TString fileNameFQD = presentDirName;   
  ((fileNameFQD+="outputFQDanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameFQD.Data(),kFileExists)))
  {
   fqdFileMerger->AddFile(fileNameFQD.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameFQD.Data()<<". Merging will continue without this file."<<endl;
    }   
  
  // LYZ1:     
  TString fileNameLYZ1 = presentDirName;   
  ((fileNameLYZ1+="outputLYZ1analysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZ1.Data(),kFileExists)))
  {
   lyz1FileMerger->AddFile(fileNameLYZ1.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameLYZ1.Data()<<". Merging will continue without this file."<<endl;
    }   
    
  // LYZ2:     
  TString fileNameLYZ2 = presentDirName;   
  ((fileNameLYZ2+="outputLYZ2analysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZ2.Data(),kFileExists)))
  {
   lyz2FileMerger->AddFile(fileNameLYZ2.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameLYZ2.Data()<<". Merging will continue without this file."<<endl;
    }     
    
  // LYZEP:     
  TString fileNameLYZEP = presentDirName;   
  ((fileNameLYZEP+="outputLYZEPanalysis")+=type.Data())+=".root";
  if(!(gSystem->AccessPathName(fileNameLYZEP.Data(),kFileExists)))
  {
   lyzepFileMerger->AddFile(fileNameLYZEP.Data());
  } else 
    {
     cout<<"WARNING: Couldn't find a file "<<fileNameLYZEP.Data()<<". Merging will continue without this file."<<endl;
    }       
   
  counter++;
  
 } // end of for(Int_t iDir=0;iDir<nDirs;++iDir)
 
 // merge everything:
 if(mcepFileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge MCEP files ----"<<endl;
  mcepFileMerger->Merge();
 } 
 if(spFileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge SP files ----"<<endl;
  spFileMerger->Merge();
 }
 if(gfcFileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge GFC files ----"<<endl;
  gfcFileMerger->Merge();
 }
 if(qcFileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge QC files ----"<<endl;
  qcFileMerger->Merge();
 }
 if(fqdFileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge FQD files ----"<<endl;
  fqdFileMerger->Merge();
 }
 if(lyz1FileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge LYZ1 files ----"<<endl;
  lyz1FileMerger->Merge();
 }
 if(lyz2FileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge LYZ2 files ----"<<endl;
  lyz2FileMerger->Merge();
 }
 if(lyzepFileMerger)
 {
  cout<<endl;
  cout<<" ---- Starting to merge LYZEP files ----"<<endl;
  lyzepFileMerger->Merge();
 }
 
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
    
  } // end of else if (mode==mLocalSource) 
  
} // end of void LoadSpreadLibrariesMO(const libModes mode)


