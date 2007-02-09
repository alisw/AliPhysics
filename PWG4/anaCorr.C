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
  return rv ;  
}

//______________________________________________________________________
void anaCorr() 
{  

  if (! gIsAnalysisLoaded ) {
    cout<<"Load ESD"<<endl;  
    LoadLib("ESD") ;   
    cout<<"Load ANALYSIS"<<endl;  
    LoadLib("ANALYSIS") ; 
    printf("Include path = %s\n", gSystem->GetIncludePath()) ; 
    LoadLib("PWG4Gamma") ; 
  }
  //AliLog::SetGlobalDebugLevel(3);
  
  // create the analysis goodies object
  AliAnalysisGoodies * ag = new AliAnalysisGoodies() ; 

  // definition of analysis tasks
  const Int_t knumberOfTasks = 4 ; 
  AliAnalysisTask * taskList[knumberOfTasks] ; 
  TClass * taskInputList[knumberOfTasks]  ; 
  TClass * taskOutputList[knumberOfTasks] ; 

//   taskList[0]       = new AliAnaGammaHadron("GammaTest") ;
//   taskInputList[0]  = TChain::Class() ; 
//   taskOutputList[0] = TObjArray::Class() ; 

//   AliAnaGammaDirect * gd = new AliAnaGammaDirect("GammaDirect") ;
//   gd->SetMinGammaPt(1.);
//   taskList[0]       = gd ;
  taskList[0]          = new AliAnaGammaDirect("GammaDirect") ;
  taskInputList[0]  = TChain::Class() ; 
  taskOutputList[0] = TObjArray::Class() ; 

  taskList[1]       = new AliAnaGammaIsolCut("GammaIC") ;
  taskInputList[1]  = TChain::Class() ; 
  taskOutputList[1] = TObjArray::Class() ; 

  taskList[2]       = new AliAnaGammaHadron("GammaHadron") ;
  taskInputList[2]  = TChain::Class() ; 
  taskOutputList[2] = TObjArray::Class() ; 

 taskList[3]       = new AliAnaGammaJet("GammaJet") ;
  taskInputList[3]  = TChain::Class() ; 
  taskOutputList[3] = TObjArray::Class() ; 

  ag->SetTasks(knumberOfTasks, taskList, taskInputList, taskOutputList) ; 

  // get the data to analyze

  // definition of Tag cuts 
  const char * runCuts = 0x0 ; 
  const char * evtCuts = "fEventTag.fNPHOSClustersMin == 1 && fEventTag.fNEMCALClustersMin == 1" ; 

  
  TString input = gSystem->Getenv("ANA_INPUT") ; 
  cout<<"INPUT 0"<<input<<endl;
  if ( input != "") {
    char argument[1024] ;  
    if ( input.Contains("tag?") ) {
      //create the ESD collection from the tag collection 
      input.ReplaceAll("tag?", "") ; 
      const char * collESD = "esdCollection.xml" ;
      ag->MakeEsdCollectionFromTagCollection(runCuts, evtCuts, input.Data(), collESD) ;
      sprintf(argument, "esd?%s", collESD) ; 
    } else if ( input.Contains("esd?") ) 
      //argument = input.Data() ; 
      sprintf(argument, "%s", input.Data()) ;
    printf("................ %s", input.Data()) ;
    ag->Process(argument) ;

  } else {

    TChain* analysisChain = new TChain("esdTree") ;
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
