//	------------------------------------------------------------
//		    	Macro for ZDC reconstruction 
//	------------------------------------------------------------
void ZDCD2R2f(Int_t totnev=1, const char *filedig="ZDCdigits.root", 
		const char *filerec="ZDCreco.root") 
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
  gAlice->SetTreeDFileName(filedig);
  // File where SDigits and Digits Tree are written
  TFile *file2 = gAlice->InitTreeFile("R",filerec);

  AliZDC *ZDC  = (AliZDC*) gAlice->GetModule("ZDC");
  merger = new AliZDCMerger();
  merger->SetMode(1);
  merger->SetBackgroundFileName("galice.root");
  ZDC->SetMerger(merger);
  printf("\n AliZDCMerger set \n");

  //
  //   Loop over events              
  //
  for(Int_t iev=0; iev<totnev; iev++) {
    gAlice->GetEvent(iev);           
    merger->SetBackgroundEventNum(iev);
    printf("	iev = %d -> Background Event Num setted\n",iev);
    gAlice->MakeTree("R",file2);
    printf("		TreeR made\n");
    ZDC->MakeBranch("R");
    printf("		ZDC branch in TreeR made\n");
    ZDC->Digits2Reco();
  }   // event loop 
}






