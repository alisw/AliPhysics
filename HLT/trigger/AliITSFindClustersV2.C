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

Int_t AliITSFindClustersV2(Char_t SlowOrFast='f',Int_t eventn=1,Char_t* path="./") {

   cerr<<"Looking for clusters...\n";

   Char_t fname[1024];
   sprintf(fname,"%s/galice.root",path);
   TFile *in=TFile::Open(fname);
   if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 2;}
   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
      cerr<<"Can't find gAlice !\n";
      return 2;
   }
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't find the ITS !\n"; return 3; }

   AliITSgeom *geom=ITS->GetITSgeom();

   Char_t fnameCluster[1024];
   sprintf(fnameCluster,"%s/AliITSclustersV2.root",path);
   TFile *out=TFile::Open(fnameCluster,"recreate");
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



