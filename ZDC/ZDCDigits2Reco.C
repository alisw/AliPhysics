//		Macro to perform ZDC reconstruction
void ZDCDigits2Reco(Int_t nev=1) 
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
	file = new TFile("galice.root","UPDATE");
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
   
   AliZDC *ZDC  = (AliZDC*) gAlice->GetModule("ZDC");
   AliZDCMerger *merger = new AliZDCMerger();
   merger->SetMode(1);
   merger->SetBackgroundFileName("galice.root");
   ZDC->SetMerger(merger);
   
//   Loop over events to be reconstructed          
   for(Int_t iev=0; iev<nev; iev++) {
       merger->SetBackgroundEventNum(iev);
       gAlice->Digits2Reco("ZDC");
   }   // event loop 
}
