#include "AliAnaTPCTrackCalib.h"
//
// This class is ment as an example of how to get access to the 
// clusters associated with reconstructed ESD tracks 
//
//

// ROOT includes
#include <TChain.h>
#include <TMath.h>

// ALIROOT includes
#include <AliTPCclusterMI.h>
#include <AliTPCcalibTracksCuts.h>
#include <AliTPCClusterParam.h>

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliAnaTPCTrackCalib)
  
//________________________________________________________________________
AliAnaTPCTrackCalib::AliAnaTPCTrackCalib() : 
  AliAnaTPCTrackBase(), fNtracks(0), fNClusters(0), fCalibTracks(0)
{
  //
  // Default constructor (should not be used)
  //
}

//________________________________________________________________________
AliAnaTPCTrackCalib::AliAnaTPCTrackCalib(const char *name) : 
  AliAnaTPCTrackBase(name), fNtracks(0), fNClusters(0), fCalibTracks(0)
{
  //
  // Normal constructor
  //
  printf("Normal constructor called with name %s \n", name);  
  // Input slot #1 works with a AliTPCCalibTracksCuts
  DefineInput(1, AliTPCcalibTracksCuts::Class());
  // Input slot #2 works with a AliTPCClusterParam
  DefineInput(2, AliTPCClusterParam::Class());
  
  // Output slot #0 writes into a TList
  // DefineOutput(0, TList::Class());
 }

//________________________________________________________________________
void AliAnaTPCTrackCalib::CreateOutputObjects() 
{
   //
   // Connect the output objects
   //
   if(fDebug>0)
    cout << "AliAnaTPCTrackCalib::CreateOutputObjects()" << endl;

   AliAnaTPCTrackBase::CreateOutputObjects();
  
   // user code to go here
   if (GetInputData(0)) printf("Input slo 0, Class_Name: %s\n", GetInputData(0)->Class_Name());
   
   fNtracks       = new TH1I("ntracks","Number of tracks", 100, 0, 400);
   fListOfHists->Add(fNtracks);
   fNClusters     = new TH1I("ncluster","Number of clusters",100, 0, 200);
   fListOfHists->Add(fNClusters);
   
   AliTPCcalibTracksCuts *cuts = (AliTPCcalibTracksCuts*)GetInputData(1);
   AliTPCClusterParam *clusterParam = (AliTPCClusterParam*)GetInputData(2);
   if (!cuts) Error("CreateOutputObjects", "No CUTS found in input slot 1");
   else {
      printf("\nCuts found :-) \n");
      cuts->Print();
   }
   if (!clusterParam) Error("CreateOutputObjects", "No CLUSTERPARAM found in input slot 2");
   
/*   if ( !fCalibTracks ) {
      OpenFile(0, "RECREATE");*/
      
   fCalibTracks = new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts);
   fListOfHists->Add(fCalibTracks);
}


//________________________________________________________________________
Int_t AliAnaTPCTrackCalib::FillTrackHistograms(Int_t nTracks, AliESDtrack* track, AliESDfriendTrack* friendTrack, AliTPCseed* seed) {
  //
  // This is the main method which rejects noise tracks and fills 
  // the histograms
  //

  if(!nTracks || !track || !friendTrack) {
    if (fDebug > 1) AliWarning("WARNING: missing track information in AliAnaTPCTrackCalib");
    return 0;
  }
  if(seed==0) {
    if (fDebug > 1) AliWarning("WARNING: Missing seed in AliAnaTPCTrackCalib");
    return 0;
  }
   
   // calibration components to go here:  
  
   if (seed) {
      fNClusters->Fill(seed->GetNumberOfClusters());
      fCalibTracks->Process(seed);   // analysis is done in fCalibTracks
   }
 
  return 1;
}      

