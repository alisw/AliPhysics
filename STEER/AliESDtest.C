//********************************************************************
//     Example of the reconstruction that generates the ESD
// Input files: 
//   a) file containing the ITS clusters
//      (the AliITSFindClustersV2.C macro can be used to generate it)
//   b) file containing the TPC clusters
//      (the AliTPCFindClusters.C macro can be used to generate it)
//   c) file containing the TRD clusters
//      (the AliTRDdigits2cluster.C macro can be used to generate it)
//   d) file containing the TOF digits
//      (the AliTOFSDigits2Digits.C macro can be used to generate it)
// Ouput file:
//      AliESDs.root containing the ESD events 
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "TFile.h"
  #include "TSystem.h"
  #include "TStopwatch.h"
  #include "TArrayF.h"

  #include "AliMagF.h"
  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliHeader.h"
  #include "AliGenEventHeader.h"

  #include "AliESD.h"
  #include "AliESDpid.h"

  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliV0vertexer.h"
  #include "AliCascadeVertexer.h"
  #include "AliITSpidESD.h"
  #include "AliITSLoader.h"

  #include "AliTPCParam.h"
  #include "AliTPCtracker.h"
  #include "AliTPCpidESD.h"
  #include "AliTPCLoader.h"

  #include "AliTRDtracker.h"
  #include "AliTRDPartID.h"

  #include "AliTOFpidESD.h"
  #include "AliTOF.h"
  #include "AliTOFGeometry.h"
#endif

extern AliRun *gAlice;
extern TFile *gFile;

Int_t AliESDtest(Int_t nev=1,Int_t run=0) {

/**** Initialization of the NewIO *******/

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }

   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"AliESDtest.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 1;
   }
   retval = rl->LoadHeader();
   if (retval) {
      cerr<<"AliESDtest.C : LoadHeader returned error"<<endl;
      delete rl;
      return 2;
   }
   gAlice=rl->GetAliRun();
       

   AliKalmanTrack::SetConvConst(
      1000/0.299792458/gAlice->Field()->SolenoidField()
   );


/**** The ITS corner ********************/

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliESDtest.C : Can not get the ITS loader"<<endl;
      return 3;
   }
   itsl->LoadRecPoints("read");

   AliITS *dITS = (AliITS*)gAlice->GetDetector("ITS");
   if (!dITS) {
      cerr<<"AliESDtest.C : Can not find the ITS detector !"<<endl;
      return 4;
   }
   AliITSgeom *geom = dITS->GetITSgeom();

   //An instance of the ITS tracker
   AliITStrackerV2 itsTracker(geom);
   
   //An instance of the ITS PID maker
   Double_t parITS[]={35.5,0.11,10.};
   AliITSpidESD itsPID(parITS);

   //An instance of the V0 finder
   Double_t cuts[]={33,  // max. allowed chi2
                    0.16,// min. allowed negative daughter's impact parameter 
                    0.05,// min. allowed positive daughter's impact parameter 
                    0.080,// max. allowed DCA between the daughter tracks
                    0.998,// max. allowed cosine of V0's pointing angle
                    0.9,  // min. radius of the fiducial volume
                    2.9   // max. radius of the fiducial volume
                   };
   AliV0vertexer vtxer(cuts);

   Double_t cts[]={33.,    // max. allowed chi2
                    0.05,   // min. allowed V0 impact parameter 
                    0.008,  // window around the Lambda mass 
                    0.035,  // min. allowed bachelor's impact parameter 
                    0.10,   // max. allowed DCA between a V0 and a track
                    0.9985, //max. allowed cosine of the cascade pointing angle
                    0.9,    // min. radius of the fiducial volume
                    2.9     // max. radius of the fiducial volume
                   };
   AliCascadeVertexer cvtxer=AliCascadeVertexer(cts);

/**** The TPC corner ********************/

   AliTPCLoader* tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
      cerr<<"AliESDtest.C : can not get the TPC loader"<<endl;
      return 5;
   }
   tpcl->LoadRecPoints("read");

   rl->CdGAFile();
   AliTPCParam *par=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
   if (!par) { 
      cerr<<"TPC parameters have not been found !\n";
      return 6;
   }
   
   //An instance of the TPC tracker
   AliTPCtracker tpcTracker(par);

   //An instance of the TPC PID maker
   Double_t parTPC[]={45.0,0.08,10.};
   AliTPCpidESD tpcPID(parTPC);


/**** The TRD corner ********************/

   AliLoader* trdl = rl->GetLoader("TRDLoader");
   if (trdl == 0x0) {
      cerr<<"AliESDtest.C : can not get the TRD loader"<<endl;
      return 5;
   }
   trdl->LoadRecPoints("read");

   //An instance of the TRD tracker
   rl->CdGAFile();
   AliTRDtracker trdTracker(gFile);  //galice.root file

/*
   //An instance of the TRD PID maker
   TFile* pidFile = TFile::Open("pid.root");
   if (!pidFile->IsOpen()) {
     cerr << "Can't get pid.root !\n";
     return 7;
   }
   AliTRDPartID* trdPID = (AliTRDPartID*) pidFile->Get("AliTRDPartID");
   if (!trdPID) {
     cerr << "Can't get PID object !\n";
     return 8;
   }
*/


/**** The TOF corner ********************/
   AliTOF *dTOF = (AliTOF*)gAlice->GetDetector("TOF");
   if (!dTOF) {
      cerr<<"AliESDtest.C : Can not find the TOF detector !"<<endl;
      return 4;
   }
   AliTOFGeometry *tofGeo = dTOF->GetGeometry();
   if (!tofGeo) {
      cerr<<"AliESDtest.C : Can not find the TOF geometry !"<<endl;
      return 4;
   }

   AliLoader* tofl = rl->GetLoader("TOFLoader");
   if (tofl == 0x0) {
      cerr<<"AliESDtest.C : can not get the TOF loader"<<endl;
      return 5;
   }
   tofl->LoadDigits("read");

   //Instance of the TOF PID maker
   Double_t parTOF[]={130.,5.};
   AliTOFpidESD tofPID(parTOF);


   //rl->UnloadgAlice();


   TFile *ef=TFile::Open("AliESDs.root","RECREATE");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TStopwatch timer;
   Int_t rc=0;
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   //The loop over events
   for (Int_t i=0; i<nev; i++) {
     cerr<<"\n\nProcessing event number : "<<i<<endl;
     AliESD *event=new AliESD(); 
     event->SetRunNumber(run);
     event->SetEventNumber(i);

     rl->GetEvent(i);
 
//***** Primary vertex reconstruction (MC vertex position, for the moment)
     TArrayF v(3);     
     rl->GetHeader()->GenEventHeader()->PrimaryVertex(v);
     Double_t vtx[3]={v[0],v[1],v[2]};
     Double_t cvtx[6]={
       0.005,
       0.000, 0.005,
       0.000, 0.000, 0.010
     };
     event->SetVertex(vtx,cvtx);
     cvtx[1]=cvtx[0]; cvtx[2]=cvtx[5]; //trackers use only the diag.elements

//***** Initial path towards the primary vertex
     tpcTracker.SetVertex(vtx,cvtx);
     TTree *tpcTree=tpcl->TreeR();
     if (!tpcTree) {
        cerr<<"Can't get the TPC cluster tree !\n";
        return 4;
     }     
     tpcTracker.LoadClusters(tpcTree);
     rc+=tpcTracker.Clusters2Tracks(event);
     tpcPID.MakePID(event);                 // preliminary PID
     AliESDpid::MakePID(event);             // for the ITS tracker

     itsTracker.SetVertex(vtx,cvtx);
     TTree *itsTree=itsl->TreeR();
     if (!itsTree) {
        cerr<<"Can't get the ITS cluster tree !\n";
        return 4;
     }     
     itsTracker.LoadClusters(itsTree);
     rc+=itsTracker.Clusters2Tracks(event);


//***** Back propagation towards the outer barrel detectors
     rc+=itsTracker.PropagateBack(event); 
     
     rc+=tpcTracker.PropagateBack(event);

     TTree *trdTree=trdl->TreeR();
     if (!trdTree) {
        cerr<<"Can't get the TRD cluster tree !\n";
        return 4;
     } 
     trdTracker.LoadClusters(trdTree);
     rc+=trdTracker.PropagateBack(event);

/*
     for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
       AliESDtrack* track = event->GetTrack(iTrack);
       trdPID->MakePID(track);
     }
*/

     TTree *tofTree=tofl->TreeD();
     if (!tofTree) {
        cerr<<"Can't get the TOF cluster tree !\n";
        return 4;
     } 
     tofPID.LoadClusters(tofTree,tofGeo);
     tofPID.MakePID(event);
     tofPID.UnloadClusters();


//***** Now the final refit at the primary vertex...
     rc+=trdTracker.RefitInward(event);
     trdTracker.UnloadClusters();

     rc+=tpcTracker.RefitInward(event);
     tpcTracker.UnloadClusters();
     tpcPID.MakePID(event);

     rc+=itsTracker.RefitInward(event); 
     itsTracker.UnloadClusters();
     itsPID.MakePID(event);


//***** Here is the combined PID
     AliESDpid::MakePID(event);


//***** Hyperon reconstruction 
     vtxer.SetVertex(vtx);
     rc+=vtxer.Tracks2V0vertices(event);            // V0 finding
     rc+=cvtxer.V0sTracks2CascadeVertices(event);   // cascade finding



//***** Some final manipulations with this event 
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

   //pidFile->Close();

   delete par;

   ef->Close();

   delete rl;

   return rc;
}
