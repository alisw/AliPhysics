void DigitMaker (Int_t evNumber=1) 
{
  /////////////////////////////////////////////////////////////////////////
  //   This macro is a small example of a ROOT macro
  //   illustrating how to read the output of GALICE
  //   and fill some histograms.
  //   
  //     Root > .L anal.C   //this loads the macro in memory
  //     Root > anal();     //by default process first event   
  //     Root > anal(2);    //process third event
  //Begin_Html
  /*
    <img src="gif/anal.gif">
  */
  //End_Html
  /////////////////////////////////////////////////////////////////////////
  
  
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file =  (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) file = new TFile("galice.root","UPDATE");
  
  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }

  TParticle *particle;
  AliT0hit  *startHit;

  Int_t buffersize=256;
  Int_t split=1;

  digits= new AliT0digit();
  TBranch *bDig=0;
  printf("Branch\n");

  
  AliT0 *T0  = (AliT0*) gAlice->GetDetector("T0");
  
 // Event ------------------------- LOOP  
  
  for (j=0; j<evNumber; j++){

   
    T0->Hit2digit(j); 

  }// event loop

  file->Write();
  file->Close();
 

}//endmacro













