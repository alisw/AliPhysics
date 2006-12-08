Bool_t gIsAnalysisLoaded = kFALSE ; 

//______________________________________________________________________
Bool_t LoadLib( const char* pararchivename) 
{
  // Loads the AliRoot required libraries from a tar file 

  Bool_t rv = kTRUE ; 
 
  char cdir[1024] ; 
  sprintf(cdir, "%s", gSystem->WorkingDirectory() ) ; 
  
  // Setup par File
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*** Building PAR archive  %s  ***\n", pararchivename);

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	AliError(Form("Cannot Build the PAR Archive %s! - Abort!", pararchivename) );

        return kFALSE ;
      }
    }

    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*** Setup PAR archive  %s     ***\n", pararchivename);
      gROOT->Macro("PROOF-INF/SETUP.C");
    }    
  }

  if ( strstr(pararchivename, "ESD") ) {
    gSystem->Load("libVMC.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libRAliEn.so") ;
    gSystem->Load("libProof.so") ;
  }

  printf("*** %s library loaded *** %s **\n", pararchivename);

  gSystem->ChangeDirectory(cdir);

  gIsAnalysisLoaded = kTRUE ; 
  return rv ;  ; 
}

//______________________________________________________________________
void ana() 
{  
  if (! gIsAnalysisLoaded ) {
    LoadLib("ESD") ; 
    LoadLib("ANALYSIS_NEW") ; 
    printf("Include path = %s\n", gSystem->GetIncludePath()) ; 
    LoadLib("AnalysisCheck") ; 
  }
  
  // create the analysis goodies object
  AliAnalysisGoodies * ag = new AliAnalysisGoodies() ; 

  // definition of analysis tasks
 
  const Int_t knumberOfTasks = 8 ; 
  AliAnalysisTask * taskList[knumberOfTasks] ; 
  TClass * taskInputList[knumberOfTasks]  ; 
  TClass * taskOutputList[knumberOfTasks] ; 

  taskList[0]       = new AliPHOSQATask("PHOS") ;
  taskInputList[0]  = TChain::Class() ; 
  taskOutputList[0] = TObjArray::Class() ; 

  taskList[1]       = new AliEMCALQATask("EMCal") ;
  taskInputList[1]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[1] = TObjArray::Class() ; 

  taskList[2]       = new AliPMDQATask("PMD") ;
  taskInputList[2]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[2] = TObjArray::Class() ; 

  taskList[3]       = new AliAnalysisTaskPt("Pt") ;
  taskInputList[3]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[3] = TObjArray::Class() ; 
  
  taskList[4]       = new AliHMPIDQATask("HMPID") ;
  taskInputList[4]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[4] = TObjArray::Class() ; 

  taskList[5]       = new AliT0QATask("T0") ;
  taskInputList[5]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[5] = TObjArray::Class() ; 

  taskList[6]       = new AliMUONQATask("MUON") ;
  taskInputList[6]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[6] = TObjArray::Class() ; 
  
  taskList[7]       = new AliFMDQATask("FMD") ;
  taskInputList[7]  = taskInputList[0] ; // only one top input container allowed 
  taskOutputList[7] = TObjArray::Class() ; 

  ag->SetTasks(knumberOfTasks, taskList, taskInputList, taskOutputList) ; 

  // get the data to analyze

  // definition of Tag cuts 
  const char * runCuts = 0x0 ; 
  const char * evtCuts = "fEventTag.fNPHOSClustersMin == 1 && fEventTag.fNEMCALClustersMin == 1" ; 

  
  TString input = gSystem->Getenv("ANA_INPUT") ; 
  if ( input != "") {
    char argument[1024] ;  
    if ( input.Contains("tag?") ) {
      //create the ESD collection from the tag collection 
      input.ReplaceAll("tag?", "") ; 
      const char * collESD = "esdCollection.xml" ;
      ag->MakeEsdCollectionFromTagCollection(runCuts, evtCuts, input.Data(), collESD) ;
      sprintf(argument, "esd?%s", collESD) ; 
    } else if ( input.Contains("esd?") ) 
      argument = input.Data() ; 
    ag->Process(argument) ;

  } else {

    TChain* analysisChain = new TChain("esdTree") ;
    //   input = "alien:///alice/cern.ch/user/a/aliprod/prod2006_2/output_pp/105/411/AliESDs.root" ; 
    //   analysisChain->AddFile(input);
    input = "AliESDs.root" ; 
    analysisChain->AddFile(input);
    ag->Process(analysisChain) ; 
  }
}

//______________________________________________________________________
void Merge(const char * xml, const char * sub, const char * out) 
{
  if (! gIsAnalysisLoaded ) 
    LoadLib("ESD") ; 
  
  AliAnalysisGoodies * ag = new AliAnalysisGoodies() ; 
  ag->Merge(xml, sub, out) ;
}

//______________________________________________________________________
void test(const char * fcollection1) 
{
  AliXMLCollection collection1(fcollection1);
 TChain* analysisChain = new TChain("esdTree");
 
 collection1.Reset();
 while (collection1.Next()) {
   cout<<"Adding "<<collection1.GetTURL()<<endl;
   analysisChain->Add(collection1.GetTURL());
 }
 
 return ;
}
