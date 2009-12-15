void AliHMPIDAnalysisTask()
{
  
  
  //
  // Macro to run the AliHMPIDAnalysisTask
  //
  //______ Input settings
  Int_t nEvents = 10000000;           //on Proof
  Int_t nOffset = 0;               //on Proof
  Int_t nFile2Xml = 10;         //on Grid - how many ESD files to chain in xml
  TString xmlName="run104157.xml";//new.xml";       //name of the xml collection created
  // char *indataset = "/home/lmolnar/CERN/alice/Dec01_rv5-25-04_gcc/test/cosmic/100923/09000100923018.10";
//   char *indataset = "/COMMON/COMMON/LHC09a4_run8101X";
  // char *indataset = "/PWG4/kleinb/LHC09a3_90023_checked";///HMPID/lmolnar/run100923";///PWG4/kleinb/LHC09a3_90023_checked";
   char *indataset = "/home/lmolnar/CERN/alice/Dec01_rv5-25-04_gcc/test/cosmic/100923/09000100923018.10/AliESDs.root";
 // char *indataset = "/media/data/Data/09000100923018.10/AliESDs.root";
 // char *indataset = "/media/data/Data/alice/sim/LHC09a1/70213/001/AliESDs.root";
  
  //______ Select Local, Proof or Intercative Grid Analysis (but just only one :-))
  Bool_t bLOCAL        = kFALSE;
  Bool_t bPROOF        = kFALSE;
  Bool_t bGRIDINT      = kTRUE;

  AliLog::SetGlobalDebugLevel(2);


  //______ Init data set chain
  TString dataset(indataset);
  TChain *chain = new TChain("esdTree");

  //______ Define settings for PROOF
  const char* proofNode = "lmolnar@alicecaf";
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");  
 //______ Create analysis manager
   AliAnalysisManager *mgr  = new AliAnalysisManager("HMPID Analysis Train", "HMPID Analysis Train");
   //______ Create input handler, default ESD
   AliESDInputHandler *esdHandler = new AliESDInputHandler();
   mgr->SetInputEventHandler(esdHandler);
   mgr->SetDebugLevel(0);
   AliLog::SetGlobalLogLevel(0);
   //______ Create default input container
   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *hmpoutput= mgr->CreateContainer("hmpoutput", TList::Class(),AliAnalysisManager::kOutputContainer,"HmpidOutput.root");
  
   if(bLOCAL) 
  {
   
   Printf("========> Running LOCAL HMPID Analysis <=========");
   gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");  
   //chain = CreateESDChain(dataset.Data());//,1,0,kFALSE,kFALSE,0);
    chain->Add(indataset);  
   gSystem->SetIncludePath("-I. -I$HOME/HMPID -I$ALICE_ROOT/include -I$ROOTSYS/include");
   gROOT->LoadMacro(Form("%s/HMPID/AliHMPIDAnalysisTask.cxx++g",gSystem->Getenv("ALICE_ROOT")));//"AliHMPIDAnalysisTask.cxx++");
   AliHMPIDAnalysisTask *hmpTask = new AliHMPIDAnalysisTask("HMPIDAnalysisTask");
   mgr->AddTask(hmpTask);  
  }
  else if(bPROOF) 
  {
    Printf("========> Running PROOF HMPID Analysis <=========");
    TProof::Open(proofNode);  
    /*
    gProof->UploadPackage("STEERBase.par");
    gProof->EnablePackage("STEERBase");	   
    gProof->UploadPackage("ESD.par");	   
    gProof->EnablePackage("ESD");	
    gProof->UploadPackage("ANALYSIS.par"); 
    gProof->EnablePackage("ANALYSIS");	   
    gProof->UploadPackage("ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");
    */
    gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-17-Release/AF-v4-17");
    gProof->EnablePackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-17-Release/AF-v4-17");   
    gProof->Load(Form("%s/HMPID/AliHMPIDAnalysisTask.cxx++g",gSystem->Getenv("ALICE_ROOT")));
    AliHMPIDAnalysisTask *hmpTask = new AliHMPIDAnalysisTask("HMPIDAnalysisTask");
    mgr->AddTask(hmpTask);    
  }
  else if(bGRIDINT) 
  {
    Printf("========> Running INTERCATIVE GRID HMPID Analysis <=========");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");  
    //______ Add HMPID task 
    gSystem->SetIncludePath("-I. -I$HOME/HMPID -I$ALICE_ROOT/include -I$ROOTSYS/include");
    gROOT->LoadMacro("$ALICE_ROOT/HMPID/AliHMPIDAnalysisTask.cxx++");
  
       
    TGrid::Connect("alien://");
    chain = CreateChainFromCollection(xmlName.Data(),"esdTree",nFile2Xml); 
    //gROOT->LoadMacro(Form("%s/HMPID/AliHMPIDAnalysisTask.cxx++g",gSystem->Getenv("ALICE_ROOT")));//"AliHMPIDAnalysisTask.cxx++");
    //gROOT->LoadMacro("$ALICE_ROOT/HMPID/AliHMPIDAnalysisTask.cxx++");
    AliHMPIDAnalysisTask *hmpTask = new AliHMPIDAnalysisTask("HMPIDAnalysisTask");
    mgr->AddTask(hmpTask);  
     
  }
  
   
  if(bLOCAL) 
  {
   //______ Add HMPID task 
   gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");  
   

  }
  if ( bPROOF ) 
  {
    
  } 
  
   mgr->ConnectInput(hmpTask,0,cinput);
   mgr->ConnectOutput(hmpTask,0,hmpoutput);
 
   
      
    if (mgr->InitAnalysis()) {
       mgr->PrintStatus();
       if(bLOCAL || bGRIDINT) mgr->StartAnalysis("local",chain);
       else if(bPROOF) mgr->StartAnalysis("proof",dataset.Data(), nEvents,nOffset);
     }
     
}   
//__________________________________________________________________________________
TChain *CreateChainFromCollection(const char* xmlfile, const char *treeName="esdTree",Int_t nFiles = 0)
{
// Create a chain from an alien collection.                                                                          
   TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);

   if (!myCollection) {
      ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
     return NULL ;
   }

  TChain* chain = new TChain(treeName);
  myCollection->Reset() ;
  Int_t iCount = 0;
  while ( myCollection->Next() ){
    if(nFiles!=0)iCount++;
    if(iCount > nFiles)break;
    chain->Add(myCollection->GetTURL("")) ;
    Printf("Adding %s",myCollection->GetTURL(""));
  }
  chain->ls();
  return chain;
}    
//__________________________________________________________________________________  
void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries                                                         
  //For testing, if par file already decompressed and modified                                        
  //classes then do not decompress.                                                                   

  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
  TString parpar(Form("%s.par", pararchivename)) ;
  
  if (!gSystem->AccessPathName(pararchivename) ) {
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }

  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);

  // check for BUILD.sh and execute                                                                   
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");

    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return -1;
    }
  }
  // check for SETUP.C and execute                                                                    
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }

  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}
  
