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

Int_t AliTPCFindClusters(Int_t n=1) {
   TFile *out=TFile::Open("AliTPCclusters.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliTPCclusters.root !\n"; return 1;}
   TFile *in=TFile::Open("rfio:galice.root");
   if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 2;}

   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     return 3;
   }

   TDirectory *cwd = gDirectory;

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC"); 
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

   AliTPCParam *dig=(AliTPCParam *)in->Get("75x40_100x60_150x60");
   if (!dig) {cerr<<"TPC parameters have not been found !\n"; return 4;}

   TStopwatch timer;

   switch (ver) {
   case 1:
      cerr<<"Making clusters...\n";
      {
       AliTPCv1 &tpc=*((AliTPCv1*)TPC);
       tpc.SetParam(dig); timer.Start(); cwd->cd(); 
       for(Int_t i=0;i<n;i++){
         printf("Processing event %d\n",i);
         gAlice->GetEvent(i);
         tpc.Hits2Clusters(out,i);
       } 
      }
      break;
   case 2:
      cerr<<"Looking for clusters...\n";
      {
	// delete gAlice; gAlice=0;
       AliTPCv2 tpc; 
       tpc.SetParam(dig); timer.Start(); cwd->cd();  
       for (Int_t i=0;i<n;i++){
	 printf("Processing event %d\n",i);
         tpc.Digits2Clusters(out,i);
	 //	 AliTPCclusterer::Digits2Clusters(dig, out, i);
       }
      }
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      return 5;
   }

   timer.Stop(); timer.Print();

   delete gAlice; gAlice=0;

   out->Close();

   in->Close();

   return 0;
}
