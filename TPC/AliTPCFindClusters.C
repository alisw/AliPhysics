/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#ifndef __CINT__
  #include <Riostream.h>
  #include "AliRun.h"
  #include "AliTPCv1.h"
  #include "AliTPCv2.h"
  #include "AliTPCParam.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliTPCFindClusters(Int_t N=-1) 
 {
   
   if (gAlice)
    {
      delete gAlice->GetRunLoader();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }
    
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0)
    {
      cerr<<"Can not open session"<<endl;
      return 1;
    }
   
   if (rl->LoadgAlice())
    {
      cerr<<"Error occured while l"<<endl;
      return 1;
    }
   AliKalmanTrack::SetConvConst(1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField());
   
   tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0)
    {
      cerr<<"Can not get TPC Loader"<<endl;
      return 1;
    }

   gAlice=rl->GetAliRun();
   if (!gAlice) {
      cerr<<"Can't get gAlice !\n";
      return 1;
   }

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC"); 
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

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

   Int_t n;
   if (N<=0) n = rl->GetNumberOfEvents();
   else n=N;
   
   TStopwatch timer;

   switch (ver) {
   case 1:
      cerr<<"Making clusters...\n";
      {
       AliTPCv1 &tpc=*((AliTPCv1*)TPC);
       tpc.SetParam(dig); timer.Start(); cwd->cd(); 
       tpc.SetLoader(tpcl);
       tpcl->LoadHits("read");
       tpcl->LoadRecPoints("recreate");
       for(Int_t i=0;i<n;i++){
         printf("Processing event %d\n",i);
         rl->GetEvent(i);
         tpc.Hits2Clusters(out,i);
       } 
      }
      break;
   case 2:
      cerr<<"Looking for clusters...\n";
      {
	// delete gAlice; gAlice=0;
       AliTPCv2 * tpc = new AliTPCv2();
       tpc->SetLoader(tpcl);
       tpcl->LoadDigits("read");
       tpcl->LoadRecPoints("recreate");

       tpc->SetParam(dig); timer.Start();
       for (Int_t i=0;i<n;i++)
        {
         printf("Processing event %d\n",i);
         tpc->Digits2Clusters(i);
         //AliTPCclusterer::Digits2Clusters(dig, out, i);
       }
       delete tpc;
      }
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      delete rl;
      return 5;
   }

   timer.Stop(); timer.Print();

   delete rl;
   return 0;
}
