/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliITSLoader.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSclustererV2.h"

  #include "TTree.h"
  #include "TStopwatch.h"
#endif

extern AliRun *gAlice;

Int_t AliITSFindClustersV2(Int_t nev=5, Char_t SlowOrFast='s') {

   cerr<<"Looking for clusters...\n";

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }
 
   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"AliITSFindClustersV2.C : Can not open session RL=NULL"<< endl;
      return 1;
   }
     
   AliITSLoader *itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliITSFindClustersV2.C : can not get ITS loader"<< endl;
      return 2;
   }

   rl->LoadKinematics();

   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"AliITSFindClustersV2.C : LoadgAlice returned error"<< endl;
      delete rl;
      return 3;
   }

   gAlice=rl->GetAliRun();
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't find the ITS !\n"; delete rl; return 3; }
   AliITSgeom *geom=ITS->GetITSgeom();

   itsl->LoadRecPoints("recreate");
   if (SlowOrFast=='s') itsl->LoadDigits("read");
   else itsl->LoadHits("read");
   
   if(SlowOrFast=='s'){
     AliITSclustererV2 clusterer(geom);

     TStopwatch timer;
     if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
     for (Int_t i=0; i<nev; i++) {
       cerr<<"Processing event number: "<<i<<endl;
       rl->GetEvent(i);
       
       TTree *out=itsl->TreeR();
       if (!out) {
	 itsl->MakeTree("R");
	 out=itsl->TreeR();
       }
       
       TTree *in=itsl->TreeD();
       if (!in) {
	 cerr<<"Can't get digits tree !\n";
	 return 4;
       }
       clusterer.Digits2Clusters(in,out);       
       itsl->WriteRecPoints("OVERWRITE");
       timer.Stop(); timer.Print();
     }

   } else{
     
     for(Int_t i=0;i<3;i++){
       ITS->SetSimulationModel(i,new AliITSsimulationFastPoints());
     }

     TStopwatch timer;
     for (Int_t i=0; i<nev; i++) {
       rl->GetEvent(i);
       if(itsl->TreeR()== 0x0) itsl->MakeTree("R");
       TTree* in = (TTree*)itsl->TreeH();
       TTree* out= (TTree*)itsl->TreeR();
       timer.Start();
       ITS->Hits2Clusters(in,out);
       timer.Stop(); timer.Print();
       itsl->WriteRecPoints("OVERWRITE");
     }
   }

   delete rl;

   return 0;
}



