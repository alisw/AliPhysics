/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#ifndef __CINT__
  #include <iostream.h>
  #include "AliTPCParam.h"
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

   AliTPCParam *par=(AliTPCParam*)in->Get("75x40_100x60_150x60");
   if (!par) {cerr<<"Can't get TPC parameters !\n"; return 3;}
 
   TStopwatch timer;

   Int_t rc=0;
   for (Int_t i=0;i<eventn;i++){
     printf("Processing event %d\n",i);
     AliTPCtracker *tracker = new AliTPCtracker(par,i);
     //Double_t xyz[]={0.,0.,0.}; tracker->SetVertex(xyz); //primary vertex
     rc=tracker->Clusters2Tracks(0,out);
     delete tracker;
   }
   timer.Stop(); timer.Print();
 
   delete par; //Thanks to Mariana Bondila

   in->Close();
   out->Close();

   return rc;
}
