TFile* AccessFile(TString inFile="galice.root", TString acctype="R");

void AliITSHits2FastRecPoints (Int_t evNumber1=0,Int_t evNumber2=0, TString inFile = "galice.root", TString outFile="galice.root", Int_t nsignal=25, Int_t size=-1) 
{
  /////////////////////////////////////////////////////////////////////////
  //   
  //   This macro creates fast recpoints, optionally on a separate file
  //   
  /////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs

  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file;
  if(outFile.Data() == inFile.Data()){
    file = AccessFile(inFile,"U");
  }
  else {
    file = AccessFile(inFile);
  }

  TFile *file2 = 0;   // possible output file for TreeR
  if(!(outFile.Data() == inFile.Data())){
      // open output file and create TreeR on it
    file2 = gAlice->InitTreeFile("R",outFile);
  }

  AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
  if (!ITS) return;

  // Set the simulation model

  for (Int_t i=0;i<3;i++) {
    ITS->SetSimulationModel(i,new AliITSsimulationFastPoints());
  }
   

  //
  // Event Loop
  //

  Int_t nbgr_ev=0;
  TStopwatch timer;

  cout << "Creating fast reconstructed points from hits for the ITS..." << endl;

  for (int ev=evNumber1; ev<= evNumber2; ev++) {
    cout << "...working on event "<< ev << " ..." << endl;
    Int_t nparticles = gAlice->GetEvent(ev);
    cout << "event         " <<ev<<endl;
    cout << "nparticles  " <<nparticles<<endl;
    gAlice->SetEvent(ev);
    if(!gAlice->TreeR() && file2 == 0) gAlice-> MakeTree("R");
    if(!gAlice->TreeR() && file2 != 0) gAlice->MakeTree("R",file2);
    ITS->MakeBranch("RF");
    if (ev < evNumber1) continue;
    if (nparticles <= 0) return;

    Int_t bgr_ev=Int_t(ev/nsignal);
    //printf("bgr_ev %d\n",bgr_ev);
    timer.Start();
    ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
    timer.Stop(); timer.Print();
  } // event loop 

  delete gAlice; gAlice=0;
  file->Close();
}


//-------------------------------------------------------------------
TFile * AccessFile(TString FileName, TString acctype){

  // Function used to open the input file and fetch the AliRun object

  TFile *retfil = 0;
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
  if (file) {file->Close(); delete file; file = 0;}
  if(acctype.Contains("U")){
    file = new TFile(FileName,"update");
  }
  if(acctype.Contains("N") && !file){
    file = new TFile(FileName,"recreate");
  }
  if(!file) file = new TFile(FileName);   // default readonly
  if (!file->IsOpen()) {
	cerr<<"Can't open "<<FileName<<" !" << endl;
	return retfil;
  } 

  // Get AliRun object from file or return if not on file
  if (gAlice) {delete gAlice; gAlice = 0;}
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
	cerr << "AliRun object not found on file"<< endl;
	return retfil;
  } 
  return file;
}
