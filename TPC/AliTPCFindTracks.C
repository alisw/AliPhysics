Int_t AliTPCFindTracks() {
   cerr<<"Looking for tracks...\n";

   TFile *out=TFile::Open("AliTPCtracks.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliTPCtracks.root !\n"; return 1;}

   TFile *in=TFile::Open("AliTPCclusters.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliTPCclusters.root !\n"; return 2;}

   AliTPCv2 TPC;
   AliTPCParam *digp= (AliTPCParam*)in->Get("75x40_100x60");
   if (!digp) {cerr<<"TPC parameters have not been found !\n"; return 3;}
   TPC.SetParam(digp);

   TStopwatch timer;
   TPC.Clusters2Tracks(out);
   timer.Stop(); timer.Print();

   in->Close();
   out->Close();

   return 0;
}
