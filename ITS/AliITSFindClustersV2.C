/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#ifndef __CINT__
  #include <iostream.h>

  #include "AliRun.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSclustererV2.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliITSFindClustersV2(Char_t SlowOrFast='s',Int_t eventn=1) {

   cerr<<"Looking for clusters...\n";

   TFile *in=TFile::Open("galice.root");
   if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 2;}
   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
      cerr<<"Can't find gAlice !\n";
      return 2;
   }
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't find the ITS !\n"; return 3; }

   AliITSgeom *geom=ITS->GetITSgeom();

   TFile *out=TFile::Open("AliITSclustersV2.root","new");
   if (!out->IsOpen()) {
      cerr<<"Delete old AliITSclustersV2.root !\n"; return 1;}
   geom->Write();

   TStopwatch timer;
   AliITSclustererV2 clusterer(geom);
   for (Int_t i=0; i<eventn; i++) {
       cerr<<"Processing event number: "<<i<<endl;
       gAlice->GetEvent(i);
       //ITS->MakeTreeC(); //To make the V1 cluster finders happy
       clusterer.SetEvent(i);
       if (SlowOrFast=='s') clusterer.Digits2Clusters(in,out);
       else                 clusterer.Hits2Clusters(in,out);
   }
   timer.Stop(); timer.Print();

   delete gAlice; gAlice=0;
   out->Close();
   in->Close();

   return 0;
}



