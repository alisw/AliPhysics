/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#ifndef __CINT__
  #include <Riostream.h>
  #include "AliTPCtracker.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliTPCPropagateBack() {
   cerr<<"Propagating tracks back through the TPC...\n";

   TFile *in=TFile::Open("AliTPCtracks.root");
   if (!in->IsOpen()) {
      cerr<<"Can't open AliTPCtracks.root !\n"; 
      return 1;
   }
   TFile *out=TFile::Open("AliTPCBackTracks.root","new");
   if (!out->IsOpen()) {
      cerr<<"Delete old AliTPCBackTracks.root !\n"; return 2;
   }
   TFile *file=TFile::Open("AliTPCclusters.root");
   if (!file->IsOpen()) {
      cerr<<"Can't open AliTPCclusters.root !\n";return 3;
   }
   AliTPCParam *param=(AliTPCParam*)file->Get("75x40_100x60_150x60");
   if (!param) {cerr<<"Can't get TPC parameters !\n"; return 4;}

   TStopwatch timer;
   AliTPCtracker *tracker=new AliTPCtracker(param);
   Int_t rc=tracker->PropagateBack(in,out);
   delete tracker;
   timer.Stop(); timer.Print();

   file->Close();
   in->Close();
   out->Close();

   return rc;
}
