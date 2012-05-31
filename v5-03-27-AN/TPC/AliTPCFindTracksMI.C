/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "AliTPCParam.h"
  #include "AliTPCtracker.h"
  #include "AliTPCtrackerMI.h"
  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"
  #include "AliESD.h"
  #include "TFile.h"
  #include "TStopwatch.h"
#endif

extern AliRun *gAlice;


Int_t AliTPCFindTracksMI(Int_t N=-1) {

   cerr<<"Looking for tracks...\n";

   if (gAlice)
    {
     delete AliRunLoader::Instance();
     delete gAlice;
     gAlice = 0x0;
    }
    
   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0)
    {
      cerr<<"Can not open session"<<endl;
      return 1;
    }
   AliTPCLoader *tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
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
   AliKalmanTrack::SetFieldMap(rl->GetAliRun()->Field());
   
   rl->CdGAFile();
   
   AliTPCParam *param=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
   if (!param) 
    {
     param=(AliTPCParam *)gDirectory->Get("75x40_100x60");
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
   for (Int_t i=0;i<eventn;i++)
    { 
      rl->GetEvent(i);
      TTree * input = tpcl->TreeR();
      if (input == 0x0)
       {
         tpcl->LoadRecPoints("read");
         input = tpcl->TreeR();
         if (input == 0x0)
          {
            cerr << "Problems with input tree (TreeR) for event " << i <<endl;
            continue;
          }
       }
      TTree * output = tpcl->TreeT();
      if (output == 0x0)
       {
         tpcl->MakeTree("T");
         output = tpcl->TreeT();
         if (output == 0x0)
          {
            cerr << "Problems with output tree (TreeT) for event " << i <<endl;
            continue;
          }
       }

      printf("Processing event %d\n",i);
      AliTPCtrackerMI *tracker = new AliTPCtrackerMI(param);
      tracker->SetIO();
      tracker->LoadClusters();
      rc=tracker->Clusters2Tracks();
      tracker->WriteTracks(output);
      tracker->UnloadClusters();
      tpcl->WriteTracks("OVERWRITE");
      //output->GetDirectory()->cd();
      //output->Write();
      delete tracker;
    }
   timer.Stop(); timer.Print();
   rl->UnloadgAlice();
   
   delete param; 
   delete rl;
   return rc;
}
