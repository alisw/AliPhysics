#ifndef __CINT__
  #include <iostream.h>
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliITSFindTracksV2() {
   cerr<<"Looking for tracks...\n";

   TFile *out=TFile::Open("AliITStracksV2.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliITStracksV2.root !\n"; return 1;}

   TFile *in=TFile::Open("AliTPCtracks.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return 2;}

   TFile *file=TFile::Open("AliITSclustersV2.root");
   if (!file->IsOpen()) {cerr<<"Can't open AliITSclustersV2.root !\n";return 3;}

   AliITSgeom *geom=(AliITSgeom*)file->Get("AliITSgeom");

   TStopwatch timer;
   AliITStrackerV2 tracker(geom);

   //Double_t xyz[]={0.,0.,0.}; tracker.SetVertex(xyz);  //primary vertex
   //Int_t flag[]={1};                                   //some default flags
   //flag[0]= 0; tracker.SetupFirstPass(flag);           //no constraint
   //flag[0]=-1; tracker.SetupSecondPass(flag);          //skip second pass

   Int_t rc=tracker.Clusters2Tracks(in,out);
   timer.Stop(); timer.Print();

   delete geom; //Thanks to Mariana Bondila

   file->Close();
   in->Close();
   out->Close();

   return rc;
}
