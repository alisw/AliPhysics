void VertexMaker(Int_t evNumber=1) 
{
  
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


  vertex= new AliT0vertex();
  TBranch *bRec=0;

  
 // Event ------------------------- LOOP  
  for (j=0; j<evNumber; j++){
    vertex->Reconstruct(j);
  }
  file->Write();
  file->Close();

} // end of macro




