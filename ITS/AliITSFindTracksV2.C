#ifndef __CINT__
  #include "Riostream.h"
  #include "TFile.h"
  #include "TStopwatch.h"

  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
#endif

Int_t AliITSFindTracksV2(Int_t nev=1) {  //number of events to process
   cerr<<"Looking for tracks...\n";
   
   if (gAlice) 
    {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
    }
 
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    if (rl == 0x0)
     {
      cerr<<"AliITSHits2DigitsDefault.C : Can not open session RL=NULL"
           << endl;
       return 3;
     }
     
    Int_t retval = rl->LoadgAlice();
    if (retval)
     {
       ::Error("AliITSHits2DigitsDefault.C","LoadgAlice returned error");
       delete rl;
       return 3;
     }
    retval = rl->LoadHeader();
    if (retval)
     {
      ::Error("AliITSHits2DigitsDefault.C","LoadHeader returned error");
      delete rl;
      return 3;
     }
    gAlice=rl->GetAliRun();
   
    
    AliITSLoader* itsloader = (AliITSLoader*)rl->GetLoader("ITSLoader");
    if (itsloader == 0x0)
     {
      ::Error("AliITSHits2DigitsDefault.C","can not get ITS loader");
      return 4;
     }

    AliLoader* tpcloader = rl->GetLoader("TPCLoader");
    if (tpcloader == 0x0)
     {
      cerr<<"AliITSHits2DigitsDefault.C : can not get TPC loader"
           << endl;
     }

   rl->GetEvent(0);

   itsloader->LoadTracks("recreate");
   tpcloader->LoadTracks("read"); 
   itsloader->LoadRawClusters("read");

   AliITS* dITS = (AliITS*)gAlice->GetDetector("ITS");
   if(!dITS)
    {
      ::Error("AliITSHits2DigitsDefault.C","Can not find ITS detector.");
      return 6;
    } // end if !fITS

   AliITSgeom *geom = dITS->GetITSgeom();
   if(geom == 0x0)
    {
      ::Error("AliITSHits2DigitsDefault.C","Can not get geometry from ITS detector.");
      return 6;
    } // end if !GetITSgeom()


   TStopwatch timer;
 
   for (Int_t i = 0;i < rl->GetNumberOfEvents(); i++)
    {
      AliITStrackerV2* tracker = new AliITStrackerV2(geom,i);
      Int_t rc=tracker->Clusters2Tracks();
      if (rc) 
       {
         ::Error("AliITSHits2DigitsDefault.C",
                 "AliITStrackerV2::Clusters2Tracks returned errror for event %d",i);
         delete tracker;
         break;
       }
    }
   timer.Stop(); timer.Print();
   delete tracker;
   delete rl;
   return rc;
}
