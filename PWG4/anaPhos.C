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
    gROOT->ProcessLine(processline) ; 	
    sprintf(processline,".! tar xvzf %s",parpar);
    gROOT->ProcessLine(processline);
   }
   gSystem->ChangeDirectory(pararchivename);
   
    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*** Building PAR archive  %s  ***\n", pararchivename);

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	printf("Cannot Build the PAR Archive %s! - Abort!", pararchivename) );
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

//______________________________________________________________________
void anaPhos(const Int_t kEvent=10)  
{ 
  if (! gIsAnalysisLoaded ) {
    LoadLib("ESD") ; 
    LoadLib("AOD") ;
    LoadLib("ANALYSIS") ; 
    LoadLib("PWG4Gamma") ; 
  }
  
  // create the analysis goodies object
  AliAnalysisGoodies * ag = new AliAnalysisGoodies() ; 

  // definition of analysis tasks
 
  // first task 
  AliAnaGammaPhos * phostask = new AliAnaGammaPhos("GammaPhos") ;
  ag->ConnectInput(phostask, TChain::Class(), 0) ; 
  ag->ConnectOuput(phostask, TTree::Class(), 0, "AOD") ;  
  AliAnalysisDataContainer * outGammaPhos = ag->ConnectOuput(phostask, TList::Class(), 1) ;

  AliAnaScale * scale = new AliAnaScale("ScaledGammaPhos") ; 
  ag->ConnectInput(scale, outGammaPhos, 0) ; 
  ag->ConnectOuput(scale, TList::Class(), 0) ; 

 
  // get the data to analyze

  // definition of Tag cuts 
  const char * runCuts = 0x0 ; 
  const char * evtCuts = 0x0 ; 
  const char * lhcCuts = 0x0 ; 
  const char * detCuts = 0x0 ; 
  
//"fEventTag.fNPHOSClustersMin == 1 && fEventTag.fNEMCALClustersMin == 1" ; 

  
  TString input = gSystem->Getenv("ANA_INPUT") ; 
  if ( input != "") {
    char argument[1024] ;  
    if ( input.Contains("tag?") ) {
      //create the ESD collection from the tag collection 
      input.ReplaceAll("tag?", "") ; 
      const char * collESD = "esdCollection.xml" ;
      ag->MakeEsdCollectionFromTagCollection(runCuts, lhcCuts, detCuts, evtCuts, input.Data(), collESD) ;
      sprintf(argument, "esd?%s", collESD) ; 
    } else if ( input.Contains("esd?") ) 
      sprintf(argument, "%s", input.Data()) ; 
    ag->Process(argument) ;

  } else {

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
      Double_t xsection = 0., ntrials = 0. ; 
      for (event = 0 ; event < kEvent ; event++) {
        sprintf(file, "%s/%d/AliESDs.root", kInDir,event) ; 
	TFile * fESD = 0 ; 
	if ( fESD = TFile::Open(file)) 
	  if ( fESD->Get("esdTree") ) { 
            printf("++++ Adding %s\n", file) ;
            analysisChain->AddFile(file);
	    Double_t * rv = ReadXsection(kInDir, event) ; 
	    xsection = xsection + rv[0] ;  
	    ntrials  = ntrials  + rv[1] ; 
	    cout << xsection << " " << ntrials << endl ; 
           
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
    
    scale->Set(xsection/ntrials) ; 
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

Double_t * ReadXsection(const char * inDir, const Int_t event)
{
  // Read the PYTHIA statistics from the file pyxsec.root created by
  // the function WriteXsection():
  // integrated cross section (xsection) and
  // the  number of Pyevent() calls (ntrials)
  // and calculate the weight per one event xsection/ntrials
  // The spectrum calculated by a user should be
  // multiplied by this weight, something like this:
  // TH1F *userSpectrum ... // book and fill the spectrum
  // userSpectrum->Scale(weight)
  //
  // Yuri Kharlov 19 June 2007

  Double_t xsection;
  UInt_t    ntrials;

  char cfile[80] ; 
  sprintf(cfile, "%s/%d/pyxsec.root", inDir, event) ; 

  TFile *file = new TFile(cfile,"readonly");
  if ( ! file ) {
    AliFatal(Form("could not open %s", cfile)) ; 
    exit(1) ;
  }
  TTree *tree = file->Get("Xsection");
  tree->SetBranchAddress("xsection",&xsection);
  tree->SetBranchAddress("ntrials",&ntrials);
  tree->GetEntry(0);
  cout << "Cross section = "<<xsection<<" mb, N trials = "<<ntrials<<endl;
  Double_t rv[2] ; 
  rv[0] = xsection ; 
  rv[1] = ntrials ; 
  return rv ;
}
