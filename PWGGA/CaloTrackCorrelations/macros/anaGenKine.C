/// \file anaGenKine.C
/// \ingroup CaloTrackCorrMacros
/// \brief Example of execution macro for analysis at generation level
///
/// Example of simple analysis of a single file
/// accessing directly the MC kinematics information
/// nothing else.
/// Execute the analysis in class AliAnaGeneratorKine
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//___________________________
/// Main execution method.
//________________________
void anaGenKine(Int_t mode)
{
  // Main
  
  //--------------------------------------------------------------------
  // Load analysis libraries
  // Look at the method below, 
  // ------------------------------------------------------------------
  
  LoadLibraries() ;
  //gSystem->ListLibraries();
  
  TChain * chain   = new TChain("TE") ;
  TChain * chainxs = new TChain("Xsection") ;
  
  chain  ->AddFile("galice.root");
  chainxs->AddFile("pyxsec.root");

  
  Double_t scale = -1;
  if(chainxs && chainxs->GetEntries() > 0)
  {
    Int_t nfiles = chainxs->GetEntries();
    
    //Get the cross section
    Double_t xsection = 0; 
    Float_t ntrials   = 0;
    
    GetAverageXsection(chainxs, xsection, ntrials);
    
    Int_t    nEventsPerFile = chain->GetEntries() / nfiles;
    
    Double_t trials = ntrials / nEventsPerFile ;
    
    scale = xsection/trials;
    
    printf("Get Cross section : nfiles  %d, nevents %d, nevents per file %d \n",nfiles, chain->GetEntries(),nEventsPerFile);     
    printf("                    ntrials %d, trials %2.2f, xs %2.2e, scale factor %2.2e\n", ntrials,trials,xsection,scale);
    
  } 
  
  printf("*********************************************\n");
  printf("number of entries # %lld, skipped %d\n", chain->GetEntries()) ; 	
  printf("*********************************************\n");
  
  if(!chain)
  { 
    printf("STOP, no chain available\n"); 
    return;
  }
  
  AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
  
  //--------------------------------------
  // Make the analysis manager
  //-------------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
  
  // MC handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
  mgr->SetMCtruthEventHandler(mcHandler);
  mgr->SetInputEventHandler(NULL);

  //mgr->SetDebugLevel(1); // For debugging, do not uncomment if you want no messages.
  
  TString outputFile = AliAnalysisManager::GetCommonFileName(); 
   
  gROOT->LoadMacro("AddTaskGenKine.C");  
  
  AliAnalysisTaskCaloTrackCorrelation *ana   = AddTaskGenKine(outputFile, scale);
  
  //-----------------------
  // Run the analysis
  //-----------------------    
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  mgr->StartAnalysis("local",chain);
  
  cout <<" Analysis ended sucessfully "<< endl ;
  
}

//--------------------------------------
/// Load the needed libraries most of them already loaded by aliroot
//--------------------------------------
void  LoadLibraries()
{
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libMatrix");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit"); // Root + libraries to if reclusterization is done
  
  gSystem->Load("libSTEERBase");
  gSystem->Load("libGui"); // Root + libraries to if reclusterization is done
  gSystem->Load("libCDB"); // Root + libraries to if reclusterization is done
  gSystem->Load("libESD"); // Root + libraries to if reclusterization is done
  gSystem->Load("libAOD");
  gSystem->Load("libRAWDatabase"); // Root + libraries to if reclusterization is done
  gSystem->Load("libProof"); 
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEER"); // Root + libraries to if reclusterization is done
  
  gSystem->Load("libRAWDatarec"); // Root + libraries to if reclusterization is done
  gSystem->Load("libRAWDatasim"); // Root + libraries to if reclusterization is done
  gSystem->Load("libVZERObase");  // Root + libraries to if reclusterization is done
  gSystem->Load("libVZEROrec");   // Root + libraries to if reclusterization is done
  
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libPWGGACaloTrackCorrelations");
  //SetupPar("PWGCaloTrackCorrBase");
  //SetupPar("PWGGACaloTrackCorrelations");
  
  // needed for plugin?
  gSystem->AddIncludePath("-I$ALICE_ROOT");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS");
  gSystem->AddIncludePath("-I./");
}

//_________________________________
/// Load par files, create analysis libraries
/// For testing, if par file already decompressed and modified
/// classes then do not decompress.
//_________________________________
void SetupPar(char* pararchivename)
{
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  
  if ( gSystem->AccessPathName(pararchivename) )
  {
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

//______________________________________________________________________________
/// Read the PYTHIA statistics from the file pyxsec.root created by
/// the function WriteXsection():
/// integrated cross section (xsection) and
/// the  number of Pyevent() calls (ntrials)
/// and calculate the weight per one event xsection/ntrials
/// The spectrum calculated by a user should be
/// multiplied by this weight, something like this:
/// TH1F *userSpectrum ... // book and fill the spectrum
/// userSpectrum->Scale(weight)
///
/// Yuri Kharlov 19 June 2007
/// Gustavo Conesa 15 April 2008
/// Add recovery of xs from pyxsec_hists.root file 15/jan/2015
//______________________________________________________________________________
void GetAverageXsection(TTree * tree, Double_t & xs, Float_t & ntr)
{
  Double_t xsection = 0;
  UInt_t    ntrials = 0;
  xs = 0;
  ntr = 0;
  
  Int_t nfiles =  tree->GetEntries()  ;
  if (tree && nfiles > 0)
  {
    tree->SetBranchAddress("xsection",&xsection);
    tree->SetBranchAddress("ntrials" ,&ntrials );
    for(Int_t i = 0; i < nfiles; i++)
    {
      tree->GetEntry(i);
      xs  += xsection ;
      ntr += ntrials ;
      cout << "xsection " <<xsection<<" ntrials "<<ntrials<<endl; 
    }
    
    xs =   xs /  nfiles;
    ntr =  ntr / nfiles;
    cout << "-----------------------------------------------------------------"<<endl;
    cout << "Average of "<< nfiles<<" files: xsection " <<xs<<" ntrials "<<ntr<<endl; 
    cout << "-----------------------------------------------------------------"<<endl;
  } 
  else cout << " >>>> Empty tree !!!! <<<<< "<<endl;
}



