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
#include "AliTPCclusterer.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#endif

Int_t AliTPCFindClustersMI(Int_t n=1) {
   
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   
   AliTPCLoader *tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
      cerr<<"Can not get TPC Loader"<<endl;
      return 1;
   }

   if (tpcl->LoadDigits()) {
      cerr<<"Error occured while loading digits"<<endl;
      return 1;
   }

   if (tpcl->LoadRecPoints("recreate")) {
      cerr<<"Error occured while loading digits"<<endl;
      return 1;
   }
   
   if (rl->LoadgAlice()) {
      cerr<<"Error occured while l"<<endl;
      return 1;
   }
   
   gAlice=rl->GetAliRun();
   if (!gAlice) {
      cerr<<"Can't get gAlice !\n";
      return 1;
   }

   TDirectory *cwd = gDirectory;

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC"); 
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";
   
   rl->CdGAFile();
   
   AliTPCParam *dig=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
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
       
       n = rl->GetNumberOfEvents();
       for (Int_t i=0;i<n;i++)
        { 
          rl->GetEvent(i);
          AliTPCclusterer clusterer(dig);
          
          TTree * input = tpcl->TreeD();
          if (input == 0x0)
           {
             cerr << "Can not get TreeD for event " << i <<endl;
             continue;
           }
          
          TTree * output = tpcl->TreeR();
          if (output == 0x0)
           {
             tpcl->MakeTree("R");
             output = tpcl->TreeR();
             if (output == 0x0)
              {
                cerr << "Problems with output tree (TreeR) for event " << i <<endl;
                continue;
              }
           }

          printf("Processing event %d\n",i); 
          clusterer.SetInput(input);
          clusterer.SetOutput(output);
          clusterer.Digits2Clusters();
          
          tpcl->WriteRecPoints("OVERWRITE");
       }
     }
     break;
   default:
     cerr<<"Invalid TPC version !\n";
     return 5;
   }
   
   timer.Stop(); timer.Print();
   
   delete rl;//cleans everything

   return 0;
}
