Int_t diaglevel=2;         // 1->Hits, 2->Spectra, 3->Statistics 


void RICHpadtestC (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
 
    if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    }
    else {
      //delete gAlice;
      gAlice = 0;
    }

    gAlice=0;
    
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (!file) file = new TFile("galice.root","UPDATE");
    
// Get AliRun object from file or create it if not on file
    
    if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    else {
      delete gAlice;
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }


    // Get pointers to RICH detector and Hits containers
       
    AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
    
    RICH->DiagnosticsFE(evNumber1,evNumber2);
    
    file->Write();

}
