Int_t AliTPCDisplayClusters() {
   cerr<<"Displaying clusters...\n";

   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

   TFile *cf=TFile::Open("AliTPCclusters.root");
   if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n"; return 3;}

   AliTPCParam *dig=(AliTPCParam *)cf->Get("75x40_100x60");
   if (!dig) {cerr<<"TPC parameters have not been found !\n"; return 2;}

   TCanvas *c1=new TCanvas("cdisplay", "Cluster display",0,0,700,730);
   TView *v=new TView(1);
   v->SetRange(-430,-560,-430,430,560,1710);
   c1->Clear();
   c1->SetFillColor(1);
   c1->SetTheta(90.);
   c1->SetPhi(0.);

   AliTPCClustersArray *ca=new AliTPCClustersArray;
   ca->Setup(dig);
   ca->SetClusterType("AliTPCcluster");
   ca->ConnectTree("Segment Tree");
   Int_t nrows=Int_t(ca->GetTree()->GetEntries());
   for (Int_t n=0; n<nrows; n++) {
       AliSegmentID *s=ca->LoadEntry(n);
       Int_t sec,row;
       dig->AdjustSectorRow(s->GetID(),sec,row);
       AliTPCClustersRow &clrow = *ca->GetRow(sec,row);
       Int_t ncl=clrow.GetArray()->GetEntriesFast();
       TPolyMarker3D *pm=new TPolyMarker3D(ncl);
       while (ncl--) {
           AliTPCcluster *cl=(AliTPCcluster*)clrow[ncl];
           Double_t x=dig->GetPadRowRadii(sec,row), y=cl->GetY(), z=cl->GetZ();
           Float_t cs, sn, tmp;
           dig->AdjustCosSin(sec,cs,sn);
           tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp;
           pm->SetPoint(ncl,x,y,z);
       }
       ca->ClearRow(sec,row);
       pm->SetMarkerSize(1); pm->SetMarkerColor(2); pm->SetMarkerStyle(1);
       pm->Draw();
   }
   delete ca;
   cf->Close();

   TGeometry *geom=(TGeometry*)file->Get("AliceGeom");
   geom->Draw("same");
   c1->Modified(); c1->Update(); 

   file->Close();
   return 0;
}
