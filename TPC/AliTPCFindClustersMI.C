/****************************************************************************
 *           Origin: M.Ivanov marian.ivanov@cern.ch                         *
 ****************************************************************************/

/*

  macro to create array of clusters from TPC digits
  input files - galice.root 
                digits.root - file with digits - usualy use link to galice.root
		            - in splitted mode - neccesary to create link to proper file
			    
   output file - AliTPCclusters.root
               - to be used by AliTPCTrackFinderMI.C

  Warning - if cluster file AliTPCclusters.root already exist - macro exit and don't produce anything
	       
 
*/


#ifndef __CINT__
#include <iostream.h>
#include "AliRun.h"
#include "AliTPCv1.h"
#include "AliTPCv2.h"
#include "AliTPCParam.h"
#include "AliTPCclustererMI.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#endif

Int_t AliTPCFindClustersMI(Int_t n=1) {
   TFile *out=TFile::Open("AliTPCclusters.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliTPCclusters.root !\n"; return 1;}
   TFile *in=TFile::Open("rfio:galice.root");
   if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 2;}

   TFile *ind=TFile::Open("rfio:digits.root");
   if (!ind->IsOpen()) {cerr<<"Can't open digits file !\n"; return 2;}


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
	 AliTPCclustererMI clusterer;
	 char dname[100];
	 char cname[100];
	 sprintf(dname,"TreeD_75x40_100x60_150x60_%d",i);
	 sprintf(cname,"TreeC_TPC_%d",i);
	 TTree * input = (TTree*)ind->Get(dname);
	 out->cd();
	 TTree * output = new TTree(cname,cname); 
	 
	 printf("Processing event %d\n",i); 
	 clusterer.SetInput(input);
	 clusterer.SetOutput(output);
	 clusterer.Digits2Clusters(dig, i);
         //tpc.Digits2Clusters(out,i);
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
