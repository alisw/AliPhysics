/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#ifndef __CINT__
  #include <iostream.h>
  #include "AliTPCParam.h"
  #include "AliTPCtracker.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliTPCFindTracks(Int_t N=-1) {

   cerr<<"Looking for tracks...\n";

   if (gAlice)
    {
     delete gAlice->GetRunLoader();
     delete gAlice;
     gAlice = 0x0;
    }
    
   rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0)
    {
      cerr<<"Can not open session"<<endl;
      return 1;
    }
   tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0)
    {
      cerr<<"Can not get TPC Loader"<<endl;
      return 1;
    }
   
   if (rl->LoadgAlice())
    {
      cerr<<"Error occured while l"<<endl;
      return 1;
    }
   AliKalmanTrack::SetConvConst(1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField());
   rl->CdGAFile();
   AliTPCParam *dig=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
   if (!dig) 
    {
     dig=(AliTPCParam *)gDirectory->Get("75x40_100x60");
     if (!param) 
      {
        cerr<<"TPC parameters have not been found !\n";
        return 1;
      }
     else
      {
        cout<<"TPC 75x40_100x60 geometry found"<<endl;
      }
    }
   else
    {
      cout<<"TPC 75x40_100x60_150x60  geometry found"<<endl;
    }

   rl->UnloadgAlice()   
   
   tpcl->LoadRecPoints("read");
   tpcl->LoadTracks("recreate");

   Int_t eventn;
   if (N<=0) 
    {
     eventn = rl->GetNumberOfEvents();
     rl->UnloadHeader();
    }
   else
    eventn = N;
    
   TStopwatch timer;
   Int_t rc=0;
   for (Int_t i=0;i<eventn;i++){
     printf("Processing event %d\n",i);
     AliTPCtracker *tracker = new AliTPCtracker(dig,i);
     rc=tracker->Clusters2Tracks();
     delete tracker;
   }
   timer.Stop(); timer.Print();
 
   delete dig; //Thanks to Mariana Bondila
   delete rl;
   return rc;
}
