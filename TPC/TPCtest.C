void TPCtest()
{
// Connect the Root Galice file containing Geometry, Kine and Hits
TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root","UPDATE");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

gAlice->GetEvent(0);
AliDetector *TPC = gAlice->GetDetector("TPC");
Int_t en;

((AliTPC*)TPC)->Hits2Digits();

TClonesArray *Digits = TPC->Digits();

en = gAlice->TreeD()->GetEntries();
printf("entries = %d\n",en);

gAlice->TreeD()->Write();

}
