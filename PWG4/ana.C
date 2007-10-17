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
   char parpar[80] ; 
   sprintf(parpar, "%s.par", pararchivename) ;
   if ( gSystem->AccessPathName(parpar) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    char processline[1024];
    sprintf(processline, ".! make %s", parpar) ; 
    cout << processline << endl ; 
    gROOT->ProcessLine(processline) ;
    gSystem->ChangeDirectory(cdir) ; 
    sprintf(processline, ".! mv /tmp/%s .", parpar) ;
    //gROOT->ProcessLine(processline) ; 	
    sprintf(processline,".! tar xvzf %s",parpar);
    //gROOT->ProcessLine(processline);
   }
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

  if ( strstr(pararchivename, "AnalysisCheck") ) {
    gSystem->Load("libSpectrum.so");
  }
  
  printf("lib%s done\n", pararchivename);

  gSystem->ChangeDirectory(cdir);

  gIsAnalysisLoaded = kTRUE ; 
  return rv ; 
}

//________________________________ana______________________________________

void ana(const Int_t kEvent=10)  
{ 

  //AliLog::SetGlobalDebugLevel(1);
  gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");
  gSystem->Load("libANALYSIS.so");
 
  if (! gIsAnalysisLoaded ) {
    //    LoadLib("ESD") ; 
    //   LoadLib("AOD") ;
    //   LoadLib("ANALYSIS") ;
    LoadLib("PWG4Gamma") ; 
  }
  
  // create the analysis goodies object
  AliAnalysisGoodies * ag = new AliAnalysisGoodies() ; 
  
  AliAnalysisTaskGamma * task = new AliAnalysisTaskGamma ("Gamma");
  ag->ConnectInput(task, TChain::Class(), 0) ; 
  ag->ConnectOuput(task, TTree::Class(), 0, "AOD") ;  
  ag->ConnectOuput(task, TList::Class(), 1) ;  
   
  // get the data to analyze

  // definition of Tag cuts 
  const char * runCuts = 0x0 ; 
  const char * evtCuts = 0x0 ; 
  const char * lhcCuts = 0x0 ; 
  const char * detCuts = 0x0 ; 
  
  //"fEventTag.fNPHOSClustersMin == 1 && fEventTag.fNEMCALClustersMin == 1" ; 
  
  TString input = gSystem->Getenv("ANA_INPUT") ; 
  cout<<"ANA_INPUT >>> "<<input<<endl;
  if ( input != "") {
    cout<<"ANALYSIS 1"<<endl;

    char argument[1024] ;  

    if ( input.Contains("tag?") ) {
      //create the ESD collection from the tag collection 
      input.ReplaceAll("tag?", "") ; 

      const char * collESD = "esdCollection.xml" ;
      ag->MakeEsdCollectionFromTagCollection(runCuts, lhcCuts, detCuts, evtCuts, input.Data(), collESD) ;
      sprintf(argument, "esd?%s", collESD) ; 
      cout<<"TAG collection >> "<<argument<<endl; 
    } 

    else if ( input.Contains("esd?") ){
      sprintf(argument, "%s", input.Data()) ;
      cout<<"ESD collection >> "<<argument<<endl; 
    } 
    
    ag->Process(argument) ;
    
  } 
  else {
    cout<<"ANALYSIS 2"<<endl;
    TChain* analysisChain = new TChain("esdTree") ;
    //   input = "alien:///alice/cern.ch/user/a/aliprod/prod2006_2/output_pp/105/411/AliESDs.root" ; 
    //   analysisChain->AddFile(input);
    input = "AliESDs.root" ; 
    const char * kInDir = gSystem->Getenv("OUTDIR") ; 
    if ( kInDir ) {
      if ( ! gSystem->cd(kInDir) ) {
	printf("%s does not exist\n", kInDir) ;
	return ;
      }     
      Int_t event, skipped=0 ; 
      char file[120] ;
      for (event = 0 ; event < kEvent ; event++) {
        sprintf(file, "%s/%d/AliESDs.root", kInDir,event) ; 
	TFile * fESD = 0 ; 
	if ( fESD = TFile::Open(file)) 
	  if ( fESD->Get("esdTree") ) { 
            printf("++++ Adding %s\n", file) ;
            analysisChain->AddFile(file);
	  }
	  else { 
            printf("---- Skipping %s\n", file) ;
            skipped++ ;
	  }
      }
      printf("number of entries # %lld, skipped %d\n", analysisChain->GetEntries(), skipped*100) ; 	
    }
    else  
      analysisChain->AddFile(input);
    
    ag->Process(analysisChain) ; 
  }
  return ;
}

//______________________________________________________________________
void Merge(const char * xml, const char * sub, const char * out) 
{
  if (! gIsAnalysisLoaded ) 
    LoadLib("ESD") ; 
  
  AliAnalysisGoodies * ag = new AliAnalysisGoodies() ; 
  ag->Merge(xml, sub, out) ;
}

