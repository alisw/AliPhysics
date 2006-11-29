void VertexMaker(Int_t evNumber=1) 
{
  
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  char filename[100];
  sprintf(filename,"galice.root");
  AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  if (rl == 0x0)
   {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  
  AliT0* T0  = (AliT0 *)gAlice->GetDetector("T0");
  /*  
  rl->LoadHeader();
  
  rl->LoadKinematics("READ");
  Int_t retval;
  AliLoader* lstart = rl->GetLoader("T0Loader");
  lstart->LoadDigits("READ");
  //  lstart->Dump();
  //  lstart->LoadRecPoints("UPDATE");
  */  
  vertex= new AliT0vertex();
  TBranch *bRec=0;
  vertex->Dump();
  
 // Event ------------------------- LOOP  
  //     Int_t iNevents=rl->GetNumberOfEvents();
  // cout<<"  nevents   "<<iNevents<<endl;

  // for (Int_t j=0; j<iNevents; j++){
  //  lstart->LoadDigits("READ");
    vertex->Reconstruct();
    //  }

} // end of macro




