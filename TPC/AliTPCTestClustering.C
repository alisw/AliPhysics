void AliTPCTestClustering() {
   const char *pname="Param1";
   const char *tname="TreeD0_Param1";

// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice=0;
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

   gAlice->GetEvent(0);

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
   int ver=TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

   AliTPCD *dig=(AliTPCD*)file->Get(pname);
   if (dig!=0) TPC->SetDigParam(dig);
   else cerr<<"Warning: default TPC parameters will be used !\n";

   switch (ver) {
   case 1:
      cerr<<"Making clusters...\n";
      TPC->Hits2Clusters();
      break;
   case 2:
      cerr<<"Looking for clusters...\n";
      TPC->Digits2Clusters();
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      return;
   }
   TClonesArray *c=TPC->Clusters();
   int n=c->GetEntriesFast();
   cerr<<"Number of clusters "<<n<<"                            \n";

   AliTPCParam *par=&TPC->GetDigParam()->GetParam();
   Float_t x[3];
   TPolyMarker3D *pm=new TPolyMarker3D(n);
   for (int i=0; i<n; i++) {
       AliTPCcluster *cl=(AliTPCcluster *)c->UncheckedAt(i);
       cl->GetXYZ(x,par);
       Double_t xx=x[0], yy=x[1], zz=x[2];
       pm->SetPoint(i,xx,yy,zz);
   }

   c1=new TCanvas("c1", "Cluster display",0,0,700,730);
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
