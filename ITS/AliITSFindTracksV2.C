#ifndef __CINT__
  #include "Riostream.h"
  #include "TFile.h"
  #include "TStopwatch.h"

  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
#endif

Int_t AliITSFindTracksV2(Int_t nev=1) {  //number of events to process
   cerr<<"Looking for tracks...\n";

   TFile *out=TFile::Open("AliITStracksV2.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliITStracksV2.root !\n"; return 1;}

   TFile *in=TFile::Open("AliTPCtracks.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return 2;}

   TFile *file=TFile::Open("AliITSclustersV2.root");
   if (!file->IsOpen()) {cerr<<"Can't open AliITSclustersV2.root !\n";return 3;}

   AliITSgeom *geom=(AliITSgeom*)file->Get("AliITSgeom");
   if (!geom) {cerr<<"Can't get AliITSgeom !\n"; return 4;}

   Int_t rc=0;
   TStopwatch timer;
   AliITStrackerV2 tracker(geom);
   for (Int_t i=0; i<nev; i++) {
     cerr<<"Processing event number : "<<i<<endl;
     tracker.SetEventNumber(i);
     //Double_t xyz[]={0.,0.,0.}, ers[]={0.,0.,0.01};//main vertex with errors
     //tracker.SetVertex(xyz,ers);
     //Int_t flag[]={1};                                 //some default flags
     //flag[0]= 0; tracker.SetupFirstPass(flag);         //no constraint
     //flag[0]=-1; tracker.SetupSecondPass(flag);        //skip second pass
     //tracker.SetLastLayerToTrackTo(2);            //track down to the layer 2
     //Int_t mask[6]={1,1,0,0,0,0};                 //not to skip pixels !
     //tracker.SetLayersNotToSkip(mask);            //
     rc=tracker.Clusters2Tracks(in,out);
   }
   timer.Stop(); timer.Print();

   delete geom; //Thanks to Mariana Bondila

   file->Close();
   in->Close();
   out->Close();

   return rc;
}
