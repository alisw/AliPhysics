//********************************************************************
//     Example of the reconstruction that generates the ESD
// Input files: 
//   a) AliTPCclusters.root containing the TPC clusters
//      (the AliTPCFindClusters.C macro can be used to generate it)
//   b) AliITSclustersV2.root containing the ITS clusters
//      (the AliITSFindClustersV2.C macro can be used to generate it)
// Ouput file:
//      AliESDs.root containing the ESD events 
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#ifndef __CINT__
  #include <Riostream.h>
  #include "TFile.h"
  #include "TStopwatch.h"

  #include "AliESD.h"
  #include "AliTPCParam.h"
  #include "AliTPCtracker.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
#endif

Int_t AliESDtest(Int_t nev=1) { 
   //File with the TPC clusters
   TFile *tpccf=TFile::Open("AliTPCclusters.root");
   if (!tpccf->IsOpen()) {
      cerr<<"Can't open AliTPCclusters.root !\n"; 
      return 2;
   }
   AliTPCParam *par=(AliTPCParam*)tpccf->Get("75x40_100x60_150x60");
   if (!par) {cerr<<"Can't get TPC parameters !\n"; return 3;}

   //An instance of the TPC tracker
   AliTPCtracker tpcTracker(par);


   //File with the ITS clusters
   TFile *itscf=TFile::Open("AliITSclustersV2.root");
   if (!itscf->IsOpen()) {
      cerr<<"Can't open AliITSclustersV2.root !\n"; 
      return 4;
   }
   AliITSgeom *geom=(AliITSgeom*)itscf->Get("AliITSgeom");
   if (!geom) {cerr<<"Can't get AliITSgeom !\n"; return 5;}

   //An instance of the ITS tracker
   AliITStrackerV2 itsTracker(geom);

   TFile *ef=TFile::Open("AliESDs.root","new");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TStopwatch timer;
   Int_t rc=0;
   //The loop over events
   for (Int_t i=0; i<nev; i++) {
     cerr<<"\n\nProcessing event number : "<<i<<endl;
     AliESD *event=new AliESD(); 
 
     tpcTracker.SetEventNumber(i); 
     tpcTracker.LoadClusters(tpccf);

     itsTracker.SetEventNumber(i); 
     itsTracker.LoadClusters(itscf);

     rc+=tpcTracker.Clusters2Tracks(event);

     rc+=itsTracker.Clusters2Tracks(event);

     rc+=itsTracker.PropagateBack(event); 
     itsTracker.UnloadClusters();
     
     rc+=tpcTracker.PropagateBack(event);
     tpcTracker.UnloadClusters();

     if (rc==0) {
        Char_t ename[100]; 
        sprintf(ename,"%d",i);
        if (!event->Write(ename)) rc++;
     } 
     if (rc) {
        cerr<<"Something bad happened...\n";
     }
     delete event;
   }
   timer.Stop(); timer.Print();

   delete geom;
   itscf->Close();
   delete par;
   tpccf->Close();
   ef->Close();

   return rc;
}
