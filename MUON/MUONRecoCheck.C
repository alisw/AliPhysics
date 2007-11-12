/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// ROOT includes
#include "TClonesArray.h"
#include "TH1.h"
#include "TFile.h"
#include <TGeoManager.h>

// STEER includes
#include "AliRun.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliStack.h"
#include "AliMagFMaps.h"
#include "AliTracker.h"

// MUON includes
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"

#include "AliMUONVTrackStore.h"

Int_t TrackCheck( Bool_t *compTrack);

void MUONRecoCheck (Int_t nEvent = 1, char* geoFilename = "geometry.root", 
                    char * pathSim="./generated/", char * filename="galice.root"){
  
  // Utility macro to check the muon reconstruction. Reconstructed tracks are compared
  // to reference tracks. The reference tracks are built from AliTrackReference for the
  // hit in chamber (0..9) and from kinematics (TreeK) for the vertex parameters.     
  
  Bool_t *compTrack;
  Bool_t compTrackOK[10];
  Int_t nClusterOk = 0;
  Int_t testTrack = 0;	
  Int_t iTrack = 0;
  AliMUONTrack* trackOK(0x0);
  Int_t trackID = 0;
  Double_t sigma2Cut = 16;  // 4 sigmas cut, sigma2Cut = 4*4
  AliMUONTrackParam *trackParam;
  Double_t x1,y1,z1,pX1,pY1,pZ1,p1;
  Double_t x2,y2,z2,pX2,pY2,pZ2,p2;
  
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("MUONRecoCheck.root", "RECREATE");
  
  TH1F *hReconstructible = new TH1F("hReconstructible"," Nb of reconstructible tracks ",15,-0.5,14.5);
  TH1F *hReco = new TH1F("hReco"," Nb of reconstructed tracks / evt",15,-0.5,14.5);
  TH1F *hNClusterComp = new TH1F("hNClusterComp"," Nb of compatible clusters / track ",15,-0.5,14.5);
  TH1F *hTestTrack = new TH1F("hTestTrack"," Reconstruction requirement / track",15,-0.5,14.5);
  TH1F *hTrackRefID = new TH1F("hTrackRefID"," track reference ID ",100,-0.5,99.5);
  
  TH1F *hResMomVertex = new TH1F("hResMomVertex"," delta P vertex (GeV/c)",100,-10.,10);
  TH1F *hResMomFirstCluster = new TH1F("hResMomFirstCluster"," delta P first cluster (GeV/c)",100,-10.,10);
  
  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if (!gGeoManager) {
    TGeoManager::Import(geoFilename);
    if (!gGeoManager) {
      Error("MUONmass_ESD", "getting geometry from file %s failed", geoFilename);
      return;
    }
  }
  
  // set  mag field 
  // waiting for mag field in CDB 
  printf("Loading field map...\n");
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());
    
  AliMUONRecoCheck rc(filename,pathSim);
  
  Int_t nevents = rc.NumberOfEvents();
  
  if (nevents < nEvent) nEvent = nevents;
  
  Int_t ievent;
  Int_t nReconstructibleTracks = 0;
  Int_t nReconstructedTracks = 0;
  Int_t nReconstructibleTracksCheck = 0;
  
  for (ievent=0; ievent<nEvent; ievent++)
  {
    if (!(ievent%10)) printf(" **** event # %d  \n",ievent);
    
    AliMUONVTrackStore* trackStore = rc.ReconstructedTracks(ievent);
    AliMUONVTrackStore* trackRefStore = rc.ReconstructibleTracks(ievent);
    
    hReconstructible->Fill(trackRefStore->GetSize());
    hReco->Fill(trackStore->GetSize());
    
    nReconstructibleTracks += trackRefStore->GetSize();
    nReconstructedTracks += trackStore->GetSize();
    
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef;
    
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) )
    {
      testTrack = 0;
      trackOK = 0x0;
      for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) 
      {
        compTrackOK[ch] = kFALSE;
      }
      
      TIter next2(trackStore->CreateIterator());
      AliMUONTrack* trackReco;
      
      while ( ( trackReco = static_cast<AliMUONTrack*>(next2()) ) )
      {
        // check if trackRef is compatible with trackReco
        compTrack = trackRef->CompatibleTrack(trackReco,sigma2Cut);
        
        iTrack = TrackCheck(compTrack);
        
        if (iTrack > testTrack) 
        {
          nClusterOk = 0;
          for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) 
          {
            if (compTrack[ch]) nClusterOk++;
            compTrackOK[ch] = compTrack[ch];
          }
          testTrack = iTrack;
          trackOK = trackReco;
        }
      }
      
      hTestTrack->Fill(testTrack);
      trackID = trackRef->GetTrackID();
      hTrackRefID->Fill(trackID);
      
      if (testTrack == 4) {     // tracking requirements verified, track is found
        nReconstructibleTracksCheck++;
        hNClusterComp->Fill(nClusterOk);
        trackParam = trackRef->GetTrackParamAtVertex();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        
        // 	printf(" Ref. track at vertex: x,y,z: %f %f %f px,py,pz,p: %f %f %f %f \n",x1,y1,z1,pX1,pY1,pZ1,p1);
        trackReco = trackOK;
        trackParam = new AliMUONTrackParam(*((AliMUONTrackParam*)(trackReco->GetTrackParamAtCluster()->First())));
        AliMUONTrackExtrap::ExtrapToVertex(trackParam,x1,y1,z1);
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        delete trackParam;
        // 	printf(" Reconst. track at vertex: x,y,z: %f %f %f px,py,pz: %f %f %f %f \n",x2,y2,z2,pX2,pY2,pZ2,p2);
        
        hResMomVertex->Fill(p2-p1);
        
        trackParam = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->First();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        // 	printf(" Ref. track at 1st cluster: x,y,z: %f %f %f px,py,pz: %f %f %f \n",x1,y1,z1,pX1,pY1,pZ1);
        trackParam = (AliMUONTrackParam*) trackOK->GetTrackParamAtCluster()->First();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        // 	printf(" Reconst. track at 1st cluster: x,y,z: %f %f %f px,py,pz: %f %f %f \n",x2,y2,z2,pX2,pY2,pZ2);
        
        hResMomFirstCluster->Fill(p2-p1);
	       
      }
    } // end loop track ref.

  } // end loop on event  
  
  delete field;
  
  printf(" nb of reconstructible tracks: %d \n", nReconstructibleTracks);
  printf(" nb of reconstructed tracks: %d \n", nReconstructedTracks);
  printf(" nb of reconstructible tracks which are reconstructed: %d \n", nReconstructibleTracksCheck);
  
  histoFile->Write();
  histoFile->Close();
}


Int_t TrackCheck( Bool_t *compTrack)
{
  // Apply reconstruction requirements
  // Return number of validated conditions 
  // If all the tests are verified then TrackCheck = 4 (good track)
  Int_t iTrack = 0;
  Int_t nCompClustersInLastStations = 0;
  
  // apply reconstruction requirements
  if (compTrack[0] || compTrack[1]) iTrack++; // at least one compatible cluster in st. 0
  if (compTrack[2] || compTrack[3]) iTrack++; // at least one compatible cluster in st. 1
  if (compTrack[4] || compTrack[5]) iTrack++; // at least one compatible cluster in st. 2
  for (Int_t ch = 6; ch < AliMUONConstants::NTrackingCh(); ch++) {
    if (compTrack[ch]) nCompClustersInLastStations++; 
  }
  if (nCompClustersInLastStations > 2) iTrack++; // at least 3 compatible clusters in st. 3 & 4
  
  return iTrack;
  
}

