//********************************************************************
//     Example of the reconstruction that generates the ESD
// Input files: 
//   a) AliTPCclusters.root containing the TPC clusters
//      (the AliTPCFindClusters.C macro can be used to generate it)
//   b) AliITSclustersV2.root containing the ITS clusters
//      (the AliITSFindClustersV2.C macro can be used to generate it)
// Ouput file:
//      AliESDs.root containing the ESD events 
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#ifndef __CINT__
  #include <Riostream.h>
  #include "TFile.h"
  #include "TStopwatch.h"

  #include "AliESD.h"
  #include "AliESDpid.h"
  #include "AliTPCpidESD.h"
  #include "AliTPCParam.h"
  #include "AliTPCtracker.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliTRDtracker.h"
  #include "AliTRDPartID.h"
  #include "AliITSpidESD.h"
#endif

Int_t AliESDtest(Int_t nev=1, 
		 const char* fileNameITSClusters = "its.clusters.root",
		 const char* fileNameTPCClusters = "tpc.clusters.root",
		 const char* fileNameTRDClusters = "trd.clusters.root") { 

   //File with the TPC clusters
   TFile *tpccf=TFile::Open(fileNameTPCClusters);
   if (!tpccf->IsOpen()) {
      cerr<<"Can't open "<<fileNameTPCClusters<<" !\n"; 
      return 2;
   }
   AliTPCParam *par=(AliTPCParam*)tpccf->Get("75x40_100x60_150x60");
   if (!par) {cerr<<"Can't get TPC parameters !\n"; return 3;}

   //An instance of the TPC tracker
   AliTPCtracker tpcTracker(par);

   //An instance of the TPC PID maker
   Double_t parTPC[]={47.,0.1,3.};
   AliTPCpidESD tpcPID(parTPC);

   //File with the ITS clusters
   TFile *itscf=TFile::Open(fileNameITSClusters);
   if (!itscf->IsOpen()) {
      cerr<<"Can't open "<<fileNameITSClusters<<".root !\n"; 
      return 4;
   }
   AliITSgeom *geom=(AliITSgeom*)itscf->Get("AliITSgeom");
   if (!geom) {cerr<<"Can't get AliITSgeom !\n"; return 5;}

   //An instance of the ITS tracker
   AliITStrackerV2 itsTracker(geom);
   
   //An instance of the ITS PID maker
   Double_t parITS[]={34.,0.12,3.};
   AliITSpidESD itsPID(parITS);

   //File with the TRD clusters
   TFile *trdcf=TFile::Open(fileNameTRDClusters);
   if (!trdcf->IsOpen()) {
      cerr<<"Can't open "<<fileNameTRDClusters<<".root !\n"; 
      return 6;
   }

   //An instance of the TRD tracker
   AliTRDtracker trdTracker(trdcf);

   //An instance of the TRD PID maker
   AliTRDPartID* trdPID = AliTRDPartID::GetFromFile();
   if (!trdPID) return 7;

   TFile *ef=TFile::Open("AliESDs.root","RECREATE");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TStopwatch timer;
   Int_t rc=0;
   //The loop over events
   for (Int_t i=0; i<nev; i++) {
     cerr<<"\n\nProcessing event number : "<<i<<endl;
     AliESD *event=new AliESD(); 
 
     tpcTracker.SetEventNumber(i); 
     tpcTracker.LoadClusters(tpccf);

     itsTracker.SetEventNumber(i); 
     itsTracker.LoadClusters(itscf);

     rc+=tpcTracker.Clusters2Tracks(event);

     rc+=itsTracker.Clusters2Tracks(event);

     rc+=itsTracker.PropagateBack(event); 
     itsTracker.UnloadClusters();

     itsPID.MakePID(event);
     
     rc+=tpcTracker.PropagateBack(event);
     tpcTracker.UnloadClusters();

     tpcPID.MakePID(event);

     trdTracker.SetEventNumber(i);
     trdcf->cd();
     trdTracker.LoadClusters();

     rc+=trdTracker.PropagateBack(event);
     trdTracker.UnloadClusters();

     for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
       AliESDtrack* track = event->GetTrack(iTrack);
       trdPID->MakePID(track);
     }

    //Here is the combined PID
     AliESDpid::MakePID(event);

     if (rc==0) {
        Char_t ename[100]; 
        sprintf(ename,"%d",i);
	ef->cd();
        if (!event->Write(ename)) rc++;
     } 
     if (rc) {
        cerr<<"Something bad happened...\n";
     }
     delete event;
   }
   timer.Stop(); timer.Print();

   trdcf->Close();
   delete geom;
   itscf->Close();
   delete par;
   tpccf->Close();
   ef->Close();

   return rc;
}
