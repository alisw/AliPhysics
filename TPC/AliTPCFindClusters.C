Int_t AliTPCFindClusters() {
   TFile *out=TFile::Open("AliTPCclusters.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliTPCclusters.root !\n"; return 1;}
   TFile *in=TFile::Open("galice.root");
   if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 2;}

   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     return 3;
   }

   TDirectory *cwd = gDirectory;

   gAlice->GetEvent(0);

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC"); 
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

   AliTPCParam *dig=(AliTPCParam *)in->Get("75x40_100x60");
   if (!dig) {cerr<<"TPC parameters have not been found !\n"; return 4;}

   TStopwatch timer;

   switch (ver) {
   case 1:
      cerr<<"Making clusters...\n";
      {
       AliTPCv1 &tpc=*((AliTPCv1*)TPC);
       tpc.SetParam(dig); timer.Start(); cwd->cd(); tpc.Hits2Clusters(out); 
      }
      break;
   case 2:
      cerr<<"Looking for clusters...\n";
      {
       delete gAlice; gAlice=0;
       AliTPCv2 tpc; 
       tpc.SetParam(dig); timer.Start(); cwd->cd(); tpc.Digits2Clusters(out); 
      }
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      return 5;
   }

   timer.Stop(); timer.Print();

   delete gAlice; gAlice=0;

   out->Close();
   in->Close();

   return 0;
}
