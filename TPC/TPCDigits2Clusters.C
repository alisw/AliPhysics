void TPCDigits2Clusters() {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gSystem->Load("libGeant3Dummy.so");      // a dummy version of Geant3
      gSystem->Load("PHOS/libPHOSdummy.so");   // the standard Alice classes 
      gSystem->Load("libgalice.so");           // the standard Alice classes 
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close();
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

   gAlice->GetEvent(0);

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
   TPC->Digits2Clusters();
   TClonesArray *c=TPC->Clusters();
   int n=c->GetEntriesFast();
   cout<<"Number of clusters "<<n<<endl;

   TPolyMarker3D *pm=new TPolyMarker3D(n);
   for (int i=0; i<n; i++) {
       AliTPCcluster *cl=(AliTPCcluster *)c->UncheckedAt(i);
       pm->SetPoint(i,cl->fX,cl->fY,cl->fZ);
   }

   c1=new TCanvas("c1", "Cluster display",0,0,575,750);
   TView *v=new TView(1);
   v->SetRange(-430,-560,-430,430,560,1710);

   c1->Clear();
   c1->SetFillColor(1);
   c1->SetTheta(90.);
   c1->SetPhi(0.);

   pm->SetMarkerSize(1);
   pm->SetMarkerColor(2);
   pm->SetMarkerStyle(1);
   pm->Draw();

   gAlice->GetGeometry()->Draw("same");
}
