/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <iostream.h>
  #include "AliTPCParam.h"
  #include "AliTPCtracker.h"
  #include "AliRun.h"
  #include "AliMagF.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

extern AliRun *gAlice;

Int_t AliTPCFindTracks(Int_t nev=5) {

   cerr<<"Looking for tracks...\n";

   if (gAlice) {
     delete gAlice->GetRunLoader();
     delete gAlice;
     gAlice = 0x0;
   }

   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }

   AliTPCLoader *tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
      cerr<<"Can not get TPC Loader"<<endl;
      return 1;
   }
   
   if (rl->LoadgAlice()) {
      cerr<<"Error occured while loading gAlice"<<endl;
      return 1;
   }
   AliKalmanTrack::SetConvConst(
      1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField()
   );
   rl->CdGAFile();
   AliTPCParam *dig=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
   if (!dig) { 
        cerr<<"TPC parameters have not been found !\n";
        return 1;
   }

   rl->UnloadgAlice();   

   tpcl->LoadRecPoints("read");
   tpcl->LoadTracks("recreate");

   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
    
   TStopwatch timer;
   Int_t rc=0;
   AliTPCtracker tracker(dig);
   for (Int_t i=0;i<nev;i++){
     printf("Processing event %d\n",i);
     rl->GetEvent(i);

     TTree *in=tpcl->TreeR();
     if (!in) {
        cerr<<"Can't get clusters tree !\n";
        return 4;
     }

     TTree *out=tpcl->TreeT();
     if (!out) {
        tpcl->MakeTree("T");
        out=tpcl->TreeT();
     }
         
     rc=tracker.Clusters2Tracks(in,out);

     tpcl->WriteTracks("OVERWRITE");
   }

   timer.Stop(); timer.Print();
 
   delete dig; //Thanks to Mariana Bondila

   delete rl;

   return rc;
}
