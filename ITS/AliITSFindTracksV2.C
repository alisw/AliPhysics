/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "Riostream.h"

  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"
  #include "AliITSLoader.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
#endif

extern AliRun *gAlice;

Int_t AliITSFindTracksV2(Int_t nev=5) {  //number of events to process
   cerr<<"Looking for tracks...\n";
   
   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }
 
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"AliITSFindTracks.C : Can not open session RL=NULL"<< endl;
      return 3;
   }
     
   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"AliITSFindTracksV2.C : LoadgAlice returned error"<<endl;
      delete rl;
       return 3;
   }
   retval = rl->LoadHeader();
   if (retval) {
      cerr<<"AliITSFindTracksV2.C : LoadHeader returned error"<<endl;
      delete rl;
      return 3;
   }
   gAlice=rl->GetAliRun();
       
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliITSFindTracksV2.C : Can not get ITS loader"<<endl;
      return 4;
   }

   AliTPCLoader* tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
      cerr<<"AliITSFindTracksV2.C : can not get TPC loader"<<endl;
      return 5;
   }

   AliITS *dITS = (AliITS*)gAlice->GetDetector("ITS");
   if (!dITS) {
      cerr<<"AliITSFindClusters.C : Can not find the ITS detector !"<<endl;
      return 6;
   }
   AliITSgeom *geom = dITS->GetITSgeom();

   AliITStrackerV2 tracker(geom);

   tpcl->LoadTracks("read"); 
   itsl->LoadTracks("recreate");
   itsl->LoadRecPoints("read");

   TStopwatch timer; 
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   Int_t rc=0;
   for (Int_t i=0; i<nev; i++) {
       cerr<<"Processing event number: "<<i<<endl;
       rl->GetEvent(i);

       TTree *cTree=itsl->TreeR();
       if (!cTree) {
	  cerr<<"AliITSFindTracksV2.C : Can't get the clusters tree !"<<endl;
          return 4;
       }
       TTree *tpcTree=tpcl->TreeT();
       if (!tpcTree) {
	  cerr<<"AliITSFindTracksV2.C : Can't get the TPC track tree !"<<endl;
          return 4;
       }
       TTree *itsTree=itsl->TreeT();
       if (!itsTree) {
          itsl->MakeTree("T");
          itsTree=itsl->TreeT();
       }

       tracker.LoadClusters(cTree);
       rc=tracker.Clusters2Tracks(tpcTree,itsTree);
       tracker.UnloadClusters();

       itsl->WriteTracks("OVERWRITE");
   }
   timer.Stop(); timer.Print();

   delete rl;

   return rc;
}
