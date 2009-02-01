/****************************************************************************
 * This macro is used to create a DataBase for the TPC tracking             *
 * parameterization.                                                        * 
 * 1) the function CreateAllGeantTracks gives all tracks at the 1st hit of  *
 *    the TPC                                                               *
 * 2) the function TrackCompare compares them with track found by the       *
 *    Kalman filter for the same event and computes efficiency and          *
 *    resolution on the track parameters for the Kalman filter.             *
 * 3) the function BuildDataBase calls many functions of AliTPCtrackerParam:*
 *    - merge results from TrackCompare for many events and compute         *
 *      average efficiency.                                                 *
 *    - analyze the pulls of the covariance matrix                          *
 *    - analyze the dE/dx                                                   *
 *    - regularize the covariance matrix as a function of the momentum      *
 *    - write all the informations and the trees with regularized cov.      *
 *      matrices in the DataBase file.                                      *
 *                                                                          *
 *  Origin: A.Dainese, Padova, andrea.dainese@pd,infn.it                    * 
 *                                                                          *
 ****************************************************************************/

#ifndef __CINT__
#include "Riostream.h"
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TObject.h>
#include "alles.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliMagF.h"
#include "AliModule.h"
#include "AliArrayI.h"
#include "AliDigits.h"
#include "AliITS.h"
#include "AliTPC.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliITSclusterV2.h"
#include "AliITSsimulationFastPoints.h"
#include "AliITStrackerV2.h"
#include "AliKalmanTrack.h"
#include "AliTPCtrackerParam.h"
#include "AliTracker.h"
#include "AliESD.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#endif

Int_t TPCParamTracks(const Int_t coll=1,Int_t firstEvent=0,Int_t lastEvent=0);
void CreateAllGeantTracks(const Int_t coll=1,Int_t nev=1);
void TrackCompare(const Int_t coll,const Double_t Bfield,Int_t n);
void BuildDataBase(const Int_t coll,const Double_t Bfield);

//_____________________________________________________________________________
Int_t TPCParamTracks(const Int_t coll,Int_t firstEvent,Int_t lastEvent) {
//
// Ordinary TPC tracking parameterization
//

   /**** Initialization of the NewIO *******/

   if (gAlice) {
      delete AliRunLoader::GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }

   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return;
   }
   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"LoadgAlice returned error"<<endl;
      delete rl;
      return;
   }
   retval = rl->LoadHeader();
   if (retval) {
      cerr<<"LoadHeader returned error"<<endl;
      delete rl;
      return;
   }
   gAlice=rl->GetAliRun();
       

   TDatabasePDG *DataBase = TDatabasePDG::Instance();

   // Get field from galice.root
   AliMagF *fiel = TGeoGlobalMagField::Instance()->GetField();
   Double_t fieval=TMath::Abs((Double_t)fiel->SolenoidField()/10.);

   /**** The TPC corner ********************/

   AliTPCtrackerParam tpcTrackerPar(coll,fieval);
   tpcTrackerPar.Init();

   /***** The TREE is born *****/
   
   TTree *esdTree=new TTree("esdTree","Tree with ESD objects");
   AliESD *event=0;
   esdTree->Branch("ESD","AliESD",&event);
   
   if(firstEvent>rl->GetNumberOfEvents()) firstEvent=rl->GetNumberOfEvents()-1;
   if(lastEvent>rl->GetNumberOfEvents())  lastEvent=rl->GetNumberOfEvents()-1;
   cout<<" Number of events: "<<1+lastEvent-firstEvent<<endl;
   
   //<----------------------------------The Loop over events begins
   TStopwatch timer;
   Int_t trc;
   for(Int_t i=firstEvent; i<=lastEvent; i++) { 
     
     cerr<<" Processing event number : "<<i<<endl;
     AliESD *event = new AliESD(); 
     event->SetRunNumber(gAlice->GetRunNumber());
     event->SetEventNumber(i);
     event->SetMagneticField(gAlice->Field()->SolenoidField());
     rl->GetEvent(i);

     if ( (trc=tpcTrackerPar.BuildTPCtracks(event)) ) {
       printf("exiting tracker with code %d in event %d\n",trc,i);
       esdTree->Fill(); delete event;
       continue;
     }

     esdTree->Fill();
     delete event;

   }//<-----------------------------------The Loop over events ends here
   timer.Stop(); timer.Print();

   //        The AliESDs.root is born
   TFile *ef = TFile::Open("AliESDs.root","RECREATE"); 
   if (!ef || !ef->IsOpen()) {cerr<<"Can't open AliESDs.root !\n"; return;}

   esdTree->Write(); //Write the TREE and close everything
   delete esdTree;
   ef->Close();

   delete rl;

   return;
}
//_____________________________________________________________________________
void CreateAllGeantTracks(const Int_t coll,Int_t nev) {
//
// Get all tracks at TPC 1st hit w/o selection and smearing
//
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="CreateAllGeantTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  TFile *outfile=TFile::Open(outname,"recreate");
  TFile *infile =TFile::Open(galice);

  AliTPCtrackerParam tracker(coll,Bfield,n);
  tracker.AllGeantTracks(); // this is to switch-off selection and smearing
  tracker.BuildTPCtracks(infile,outfile);

  delete gAlice; gAlice=0;

  infile->Close();
  outfile->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//_____________________________________________________________________________
void TrackCompare(const Int_t coll,const Double_t Bfield,Int_t n) {
//
// Compare Kalman tracks with tracks at TPC 1st hit
//
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="TrackCompare";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  AliTPCtrackerParam tracker(coll,Bfield,n);
  tracker.CompareTPCtracks();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//_____________________________________________________________________________
void BuildDataBase(const Int_t coll,const Double_t Bfield) {
//
//
//
  cerr<<"\n*******************************************************************\n";

  AliTPCtrackerParam tracker(coll,Bfield);

  // Merge files with cov. matrix and compute average efficiencies
  cerr<<"\n   --- Merging Events ---\n\n";
  tracker.MergeEvents(1,5);
 
  // Compute the pulls for pions, kaons, electrons
  cerr<<"\n   --- Analyzing Pulls ---\n\n";
  tracker.AnalyzePulls("pulls.root");

  // Draw pulls and efficiencies  
  tracker.DrawPulls("CovMatrixDB_PbPb6000_B0.4T.root",211,0);
  tracker.DrawEffs("CovMatrixDB_PbPb6000_B0.4T.root",13);

  // Regularize the covariance matrix
  tracker.RegularizeCovMatrix("regPi.root",211);

  // Analyze the dE/dx
  tracker.AnalyzedEdx("dEdxPi.root",211);


  // Put everything together and create the DB file
  tracker.MakeDataBase();

  return;
}







