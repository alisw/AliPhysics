/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "Riostream.h"
  #include "TKey.h"
  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliMagF.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"
  #include "AliITSLoader.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliESD.h"
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
     
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliITSFindTracksV2.C : Can not get ITS loader"<<endl;
      return 4;
   }

   if (rl->LoadgAlice()) {
      cerr<<"AliITSFindTracksV2.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 3;
   }

   AliKalmanTrack::SetConvConst(
      1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField()
   );
       
   AliITS *dITS = (AliITS*)rl->GetAliRun()->GetDetector("ITS");
   if (!dITS) {
      cerr<<"AliITSFindClusters.C : Can not find the ITS detector !"<<endl;
      return 6;
   }
   AliITSgeom *geom = dITS->GetITSgeom();

   AliITStrackerV2 tracker(geom);

   itsl->LoadRecPoints("read");

   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   Int_t rc=0;

   TFile *itsf=TFile::Open("AliESDits.root","RECREATE");
   if ((!itsf)||(!itsf->IsOpen())) {
      cerr<<"Can't AliESDits.root !\n"; return 1;
   }
   TFile *tpcf=TFile::Open("AliESDtpc.root");
   if ((!tpcf)||(!tpcf->IsOpen())) {
      cerr<<"Can't AliESDtpc.root !\n"; return 1;
   }
   TKey *key=0;
   TIter next(tpcf->GetListOfKeys());
   TStopwatch timer; 
   for (Int_t i=0; i<nev; i++) {
       tpcf->cd();
       if ((key=(TKey*)next())==0) break;
       cerr<<"Processing event number: "<<i<<endl;
       AliESD *event=(AliESD*)key->ReadObj();

       rl->GetEvent(i);

       TTree *cTree=itsl->TreeR();
       if (!cTree) {
	  cerr<<"AliITSFindTracksV2.C : Can't get the clusters tree !"<<endl;
          return 4;
       }

       tracker.LoadClusters(cTree);
       rc=tracker.Clusters2Tracks(event);
       tracker.UnloadClusters();

       if (rc==0) {
          Char_t ename[100]; 
          sprintf(ename,"%d",i);
          itsf->cd();
          if (!event->Write(ename)) rc++;
       } 
       if (rc) {
          cerr<<"Something bad happened...\n";
       }
       delete event;
   }
   timer.Stop(); timer.Print();

   tpcf->Close();
   itsf->Close();

   delete rl;

   return rc;
}
