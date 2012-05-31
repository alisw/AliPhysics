#ifndef __CINT__
  #include <Riostream.h>
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliITSrefitV2(Int_t nev=1) {
   cerr<<"Propagating tracks inward through the ITS...\n";

   TFile *in=TFile::Open("AliTPCrefited.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliTPCrefited.root !\n"; return 1;}

   TFile *out=TFile::Open("AliITStracksV2.root","update");
   if (!out->IsOpen()) {
      cerr<<"Can't open AliITStracksV2.root !\n"; return 2;
   }
   TFile *file=TFile::Open("AliITSclustersV2.root");
   if (!file->IsOpen()) {
      cerr<<"Can't open AliITSclustersV2.root !\n";return 3;
   }
   AliITSgeom *geom=(AliITSgeom*)file->Get("AliITSgeom");

   Int_t rc=0;
   TStopwatch timer;
   AliITStrackerV2 tracker(geom);
   for (Int_t i=0; i<nev; i++) {
     cerr<<"Processing event number : "<<i<<endl;
     tracker.SetEventNumber(i);
     rc=tracker.RefitInward(in,out);
   }
   timer.Stop(); timer.Print();

   delete geom;

   file->Close();
   in->Close();
   out->Close();

   return rc;
}
