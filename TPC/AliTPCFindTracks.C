#ifndef __CINT__
  #include <iostream.h>
  #include "AliTPCtracker.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliTPCFindTracks(Int_t eventn=1) { 
   cerr<<"Looking for tracks...\n";

   TFile *out=TFile::Open("AliTPCtracks.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliTPCtracks.root !\n"; return 1;}

   TFile *in=TFile::Open("AliTPCclusters.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliTPCclusters.root !\n"; return 2;}

   AliTPCParam *par=(AliTPCParam*)in->Get("75x40_100x60");
   if (!par) {cerr<<"Can't get TPC parameters !\n"; return 3;}
 
   TStopwatch timer;

   for (Int_t i=0;i<eventn;i++){
     printf("Processing event %d\n",i);
     AliTPCtracker *tracker = new AliTPCtracker(par,i);
     Int_t rc=tracker->Clusters2Tracks(0,out);
     delete tracker;
   }
   timer.Stop(); timer.Print();
 
   in->Close();
   out->Close();

   return 1;
}
