//	------------------------------------------------------------
//		    	Macro for ZDC digitization 
//	Macro arguments: totnev = num. of events to be processed
//	mode = 0 -> Only digitization (without merging) is performed
//	mode = 1 -> Signal for spectators is added to the generated
//		background and then the complete event is digitized
//	------------------------------------------------------------
void ZDCH2D2f(Int_t totnev=1, const char *filedig="ZDCdigits.root", Int_t mode=1) 
{
  delete gAlice;
  gAlice=0;
  
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
      
  // Connect the Root Galice file containing Geometry, Kine, Hits and Digits
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("galice.root");
  } else {
	printf("\n galice.root found in file list");
  }

  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) {
      printf("\n create new gAlice object");
      gAlice = new AliRun("gAlice","Alice test program");
	}
  }
  // File where SDigits and Digits Tree are written
  TFile *file2 = gAlice->InitTreeFile("SD", filedig);

  AliZDC *ZDC  = (AliZDC*) gAlice->GetModule("ZDC");
  if (ZDC && mode) {
    //printf("\n	Initializing AliZDCMerger\n");
    merger = new AliZDCMerger();
    merger->SetMode(mode);
    merger->SetBackgroundFileName("galice.root");
    ZDC->SetMerger(merger);
  }

  //
  //   Loop over events              
  //
  for(Int_t iev=0; iev<totnev; iev++) {
    gAlice->GetEvent(iev);           
    if(mode){
      merger->SetBackgroundEventNum(iev);
    }
    gAlice->MakeTree("S",file2);
    gAlice->MakeTree("D",file2);
    ZDC->MakeBranch("S");
    ZDC->Hits2SDigits();
    ZDC->MakeBranch("D");
    ZDC->SDigits2Digits();
  }   // event loop 
}






