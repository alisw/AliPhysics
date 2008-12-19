/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"
  #include "AliTPCv1.h"
  #include "AliTPCParam.h"
  #include "AliTPCclusterer.h"

  #include "TTree.h"
  #include "TStopwatch.h"
#endif

extern AliRun *gAlice;

Int_t AliTPCFindClusters(Int_t nev=5) {

   if (gAlice) {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
   }
    
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   
   if (rl->LoadgAlice()) {
      cerr<<"Error occured while l"<<endl;
      return 1;
   }
   
   AliTPCLoader *tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
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
   if (!dig) {
        cerr<<"TPC parameters have not been found !\n";
        return 1;
   }

   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   
   tpcl->LoadRecPoints("recreate");
   if (ver==1) tpcl->LoadHits("read");
   else tpcl->LoadDigits("read");

   TStopwatch timer;

   if (ver==1) {
      cerr<<"Making clusters...\n";
      AliTPCv1 &tpc=*((AliTPCv1*)TPC);
      tpc.SetParam(dig);
      tpc.SetLoader(tpcl);
      rl->LoadKinematics();
      timer.Start();
      for(Int_t i=0;i<nev;i++) {
         printf("Processing event %d\n",i);
         rl->GetEvent(i);
         tpc.Hits2Clusters(i);
      }
   } else if (ver==2) {
      cerr<<"Looking for clusters...\n";
      AliTPCclusterer *dummy=new AliTPCclusterer(dig), &clusterer=*dummy; 
      timer.Start();
      for (Int_t i=0;i<nev;i++) {
         printf("Processing event %d\n",i);
         rl->GetEvent(i);

         TTree *out=tpcl->TreeR();
         if (!out) {
            tpcl->MakeTree("R");
            out=tpcl->TreeR();
         }
         TTree *in=tpcl->TreeD();
         if (!in) {
            cerr<<"Can't get digits tree !\n";
            return 4;
         }

         clusterer.Digits2Clusters(in,out);
         
         tpcl->WriteRecPoints("OVERWRITE");
      }
      delete dummy;
      delete dig;
   } else {
      cerr<<"Invalid TPC version !\n";
      delete rl;
      return 5;
   }

   timer.Stop(); timer.Print();

   delete rl;

   return 0;
}
