void TPCHits2Digits()
{
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gSystem->Load("libGeant3Dummy.so");      // a dummy version of Geant3
      gSystem->Load("PHOS/libPHOSdummy.so");   // the standard Alice classes 
      gSystem->Load("libgalice.so");           // the standard Alice classes 
   } else {
      delete gAlice;
      gAlice=0;
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close();
   file = new TFile("galice.root","UPDATE");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

   gAlice->GetEvent(0);

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
   TPC->Hits2Digits();
   printf("entries = %d\n",gAlice->TreeD()->GetEntries());

   gAlice->TreeD()->Write("TreeD1");
}
