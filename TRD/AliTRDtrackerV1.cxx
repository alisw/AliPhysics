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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Track finder                                                             //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdio.h>
#include <string.h>

#include <TBranch.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLinearFitter.h>
#include <TObjArray.h> 
#include <TROOT.h>
#include <TTree.h>  
#include <TClonesArray.h>
#include <TRandom.h>
#include <TTreeStream.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAlignObj.h"
#include "AliRieman.h"
#include "AliTrackPointArray.h"

#include "AliTRDtracker.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliTRDseed.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDReconstructor.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliTRDtrackerFitter.h"
#include "AliTRDstackLayer.h"
#include "AliTRDrecoParam.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "Cal/AliTRDCalDet.h"

#define DEBUG

ClassImp(AliTRDtrackerV1)
Double_t AliTRDtrackerV1::fgTopologicQA[kNConfigs] = {
		0.1112, 0.1112, 0.1112, 0.0786, 0.0786,
		0.0786, 0.0786, 0.0579, 0.0579, 0.0474,
		0.0474, 0.0408, 0.0335, 0.0335, 0.0335
};

//____________________________________________________________________
AliTRDtrackerV1::AliTRDtrackerV1(AliTRDrecoParam *p) 
  :AliTRDtracker()
  ,fSieveSeeding(0)
  ,fTracklets(0x0)
  ,fRecoParam(p)
  ,fFitter(0x0)
{
  //
  // Default constructor. Nothing is initialized.
  //

}

//____________________________________________________________________
AliTRDtrackerV1::AliTRDtrackerV1(const TFile *in, AliTRDrecoParam *p) 
  :AliTRDtracker(in)
  ,fSieveSeeding(0)
  ,fTracklets(0x0)
  ,fRecoParam(p)
  ,fFitter(0x0)
{
  //
  // Standard constructor.
  // Setting of the geometry file, debug output (if enabled)
  // and initilize fitter helper.
  //

	fFitter = new AliTRDtrackerFitter();

#ifdef DEBUG
	fFitter->SetDebugStream(fDebugStreamer);
#endif

}
  
//____________________________________________________________________
AliTRDtrackerV1::~AliTRDtrackerV1()
{ 
  //
  // Destructor
  //

	if(fFitter) delete fFitter;
	if(fRecoParam) delete fRecoParam;
	if(fTracklets) {fTracklets->Delete(); delete fTracklets;}
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::Clusters2Tracks(AliESDEvent *esd)
{
  //
  // Steering stand alone tracking for full TRD detector
  //
  // Parameters :
  //   esd     : The ESD event. On output it contains 
  //             the ESD tracks found in TRD.
  //
  // Output :
  //   Number of tracks found in the TRD detector.
  // 
  // Detailed description
  // 1. Launch individual SM trackers. 
  //    See AliTRDtrackerV1::Clusters2TracksSM() for details.
  //

	if(!fRecoParam){
		AliError("Reconstruction configuration not initialized. Call first AliTRDtrackerV1::SetRecoParam().");
		return 0;
	}
	
	//AliInfo("Start Track Finder ...");
	Int_t ntracks = 0;
	for(int ism=0; ism<AliTRDtracker::kTrackingSectors; ism++){
			//AliInfo(Form("Processing supermodule %i ...", ism));
			ntracks += Clusters2TracksSM(fTrSec[ism], esd);
	}
	AliInfo(Form("Found %d TRD tracks.", ntracks));
	return ntracks;
}


//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::GetTrackPoint(Int_t /*index*/, AliTrackPoint &/*p*/) const
{
	//AliInfo(Form("Asking for tracklet %d", index));
	
	if(index<0) return kFALSE;
	//AliTRDseedV1 *tracklet = (AliTRDseedV1*)fTracklets->UncheckedAt(index);
	// etc
	return kTRUE;
}


//_____________________________________________________________________________
Int_t AliTRDtrackerV1::PropagateBack(AliESDEvent *event) 
{
  //
  // Gets seeds from ESD event. The seeds are AliTPCtrack's found and
  // backpropagated by the TPC tracker. Each seed is first propagated 
  // to the TRD, and then its prolongation is searched in the TRD.
  // If sufficiently long continuation of the track is found in the TRD
  // the track is updated, otherwise it's stored as originaly defined 
  // by the TPC tracker.   
  //  

	Int_t   found    = 0;     // number of tracks found
	Float_t foundMin = 20.0;
	
	AliTRDseed::SetNTimeBins(fTimeBinsPerPlane);
	Int_t    nSeed   = event->GetNumberOfTracks();
	if(!nSeed){
		// run stand alone tracking
		if (AliTRDReconstructor::SeedingOn()) Clusters2Tracks(event);
		return 0;
	}
	
	Float_t *quality = new Float_t[nSeed];
	Int_t   *index   = new Int_t[nSeed];
	for (Int_t iSeed = 0; iSeed < nSeed; iSeed++) {
		AliESDtrack *seed = event->GetTrack(iSeed);
		Double_t covariance[15];
		seed->GetExternalCovariance(covariance);
		quality[iSeed] = covariance[0] + covariance[2];
	}
	// Sort tracks according to covariance of local Y and Z
	TMath::Sort(nSeed,quality,index,kFALSE);
	
	// Backpropagate all seeds
	for (Int_t iSeed = 0; iSeed < nSeed; iSeed++) {
	
		// Get the seeds in sorted sequence
		AliESDtrack *seed = event->GetTrack(index[iSeed]);
	
		// Check the seed status
		ULong_t status = seed->GetStatus();
		if ((status & AliESDtrack::kTPCout) == 0) continue;
		if ((status & AliESDtrack::kTRDout) != 0) continue;
	
		// Do the back prolongation
		Int_t   lbl         = seed->GetLabel();
		AliTRDtrackV1 *track  = new AliTRDtrackV1(*seed);
		//track->Print();
		track->SetSeedLabel(lbl);
		seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup); // Make backup
		fNseeds++;
		Float_t p4          = track->GetC();
		Int_t   expectedClr = FollowBackProlongation(*track);
		//AliInfo(Form("\nTRACK %d Clusters %d [%d] in chi2 %f", index[iSeed], expectedClr, track->GetNumberOfClusters(), track->GetChi2()));
		//track->Print();

		//Double_t cov[15];
		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 0", index[iSeed], cov[0], cov[2]));

		if ((TMath::Abs(track->GetC() - p4) / TMath::Abs(p4) < 0.2) ||
				(track->Pt() > 0.8)) {
			//
			// Make backup for back propagation
			//
			Int_t foundClr = track->GetNumberOfClusters();
			if (foundClr >= foundMin) {
				//AliInfo(Form("Making backup track ncls [%d]...", foundClr));
				track->CookdEdx();
				track->CookdEdxTimBin(seed->GetID()); // A.Bercuci 25.07.07
				CookLabel(track,1 - fgkLabelFraction);
				if (track->GetBackupTrack()) UseClusters(track->GetBackupTrack());
				
				
		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 0 test", index[iSeed], cov[0], cov[2]));

				// Sign only gold tracks
				if (track->GetChi2() / track->GetNumberOfClusters() < 4) {
					if ((seed->GetKinkIndex(0)      ==   0) &&
							(track->Pt()                <  1.5)) UseClusters(track);
				}
				Bool_t isGold = kFALSE;
	
				// Full gold track
				if (track->GetChi2() / track->GetNumberOfClusters() < 5) {
					if (track->GetBackupTrack()) seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);

					isGold = kTRUE;
				}
		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 00", index[iSeed], cov[0], cov[2]));
	
				// Almost gold track
				if ((!isGold)  && (track->GetNCross() == 0) &&
						(track->GetChi2() / track->GetNumberOfClusters()  < 7)) {
					//seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
					if (track->GetBackupTrack()) seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);
					
					isGold = kTRUE;
				}
		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 01", index[iSeed], cov[0], cov[2]));
				
				if ((!isGold) && (track->GetBackupTrack())) {
					if ((track->GetBackupTrack()->GetNumberOfClusters() > foundMin) && ((track->GetBackupTrack()->GetChi2()/(track->GetBackupTrack()->GetNumberOfClusters()+1)) < 7)) {
						seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);
						isGold = kTRUE;
					}
				}
		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 02", index[iSeed], cov[0], cov[2]));
	
				//if ((track->StatusForTOF() > 0) && (track->GetNCross() == 0) && (Float_t(track->GetNumberOfClusters()) / Float_t(track->GetNExpected())  > 0.4)) {
					//seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
				//}
			}
		}
		/**/
	
		/**/
		// Debug part of tracking
/*		TTreeSRedirector &cstream = *fDebugStreamer;
		Int_t eventNrInFile = event->GetEventNumberInFile(); // This is most likely NOT the event number you'd like to use. It has nothing to do with the 'real' event number.
		if (AliTRDReconstructor::StreamLevel() > 0) {
			if (track->GetBackupTrack()) {
				cstream << "Tracks"
				<< "EventNrInFile="  << eventNrInFile
				<< "ESD.="     << seed
				<< "trd.="     << track
				<< "trdback.=" << track->GetBackupTrack()
				<< "\n";
			}
			else {
				cstream << "Tracks"
				<< "EventNrInFile="  << eventNrInFile
				<< "ESD.="     << seed
				<< "trd.="     << track
				<< "trdback.=" << track
				<< "\n";
			}
		}*/
		/**/
	
		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 1", index[iSeed], cov[0], cov[2]));

		// Propagation to the TOF (I.Belikov)
		if (track->GetStop() == kFALSE) {
			//AliInfo("Track not stopped in TRD ...");
			Double_t xtof  = 371.0;
			Double_t xTOF0 = 370.0;
		
			Double_t c2    = track->GetSnp() + track->GetC() * (xtof - track->GetX());
			if (TMath::Abs(c2) >= 0.99) {
				delete track;
				continue;
			}
			
			PropagateToX(*track,xTOF0,fgkMaxStep);
	
			// Energy losses taken to the account - check one more time
			c2 = track->GetSnp() + track->GetC() * (xtof - track->GetX());
			if (TMath::Abs(c2) >= 0.99) {
				delete track;
				continue;
			}
			
			//if (!PropagateToX(*track,xTOF0,fgkMaxStep)) {
			//	fHBackfit->Fill(7);
			//delete track;
			//	continue;
			//}
	
			Double_t ymax = xtof * TMath::Tan(0.5 * AliTRDgeometry::GetAlpha());
			Double_t y;
			track->GetYAt(xtof,GetBz(),y);
			if (y >  ymax) {
				if (!track->Rotate( AliTRDgeometry::GetAlpha())) {
					delete track;
					continue;
				}
			}else if (y < -ymax) {
				if (!track->Rotate(-AliTRDgeometry::GetAlpha())) {
					delete track;
					continue;
				}
			}
					
			if (track->PropagateTo(xtof)) {
				//AliInfo("set kTRDout");
				seed->UpdateTrackParams(track,AliESDtrack::kTRDout);
	
				for (Int_t i = 0; i < AliESDtrack::kNPlane; i++) {
					for (Int_t j = 0; j < AliESDtrack::kNSlice; j++) {
						seed->SetTRDsignals(track->GetPIDsignals(i,j),i,j);
					}
					seed->SetTRDTimBin(track->GetPIDTimBin(i),i);
				}
				//seed->SetTRDtrack(new AliTRDtrack(*track));
				if (track->GetNumberOfClusters() > foundMin) found++;
			}
		} else {
			//AliInfo("Track stopped in TRD ...");
			
			if ((track->GetNumberOfClusters() >              15) &&
					(track->GetNumberOfClusters() > 0.5*expectedClr)) {
				seed->UpdateTrackParams(track,AliESDtrack::kTRDout);
	
				//seed->SetStatus(AliESDtrack::kTRDStop);
				for (Int_t i = 0; i < AliESDtrack::kNPlane; i++) {
					for (Int_t j = 0; j <AliESDtrack::kNSlice; j++) {
						seed->SetTRDsignals(track->GetPIDsignals(i,j),i,j);
					}
					seed->SetTRDTimBin(track->GetPIDTimBin(i),i);
				}
				//seed->SetTRDtrack(new AliTRDtrack(*track));
				found++;
			}
		}

		//if (((t->GetStatus()&AliESDtrack::kTRDout)!=0 )
	
		seed->SetTRDQuality(track->StatusForTOF());
		seed->SetTRDBudget(track->GetBudget(0));
		
// 		if ((seed->GetStatus()&AliESDtrack::kTRDin)!=0 ) printf("TRDin ");	
// 		if ((seed->GetStatus()&AliESDtrack::kTRDbackup)!=0 ) printf("TRDbackup ");	
// 		if ((seed->GetStatus()&AliESDtrack::kTRDStop)!=0 ) printf("TRDstop ");	
// 		if ((seed->GetStatus()&AliESDtrack::kTRDout)!=0 ) printf("TRDout ");	
// 		printf("\n");
		delete track;

		//seed->GetExternalCovariance(cov);
		//AliInfo(Form("track %d cov[%f %f] 2", index[iSeed], cov[0], cov[2]));
	}
	

	AliInfo(Form("Number of seeds: %d",fNseeds));
	AliInfo(Form("Number of back propagated TRD tracks: %d",found));
		
	//fSeeds->Clear(); 
	fNseeds = 0;
	
	delete [] index;
	delete [] quality;
	
  return 0;
}


//____________________________________________________________________
Int_t AliTRDtrackerV1::RefitInward(AliESDEvent *event)
{
  //
  // Refits tracks within the TRD. The ESD event is expected to contain seeds 
  // at the outer part of the TRD. 
  // The tracks are propagated to the innermost time bin 
  // of the TRD and the ESD event is updated
  // Origin: Thomas KUHR (Thomas.Kuhr@cern.ch)
  //

  Int_t   nseed    = 0; // contor for loaded seeds
  Int_t   found    = 0; // contor for updated TRD tracks
  AliTRDtrackV1 track;
  for (Int_t itrack = 0; itrack < event->GetNumberOfTracks(); itrack++) {
    AliESDtrack *seed = event->GetTrack(itrack);
		new(&track) AliTRDtrackV1(*seed);

    if (track.GetX() < 270.0) {
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDbackup);
      //AliInfo(Form("Remove for X = %7.3f [270.]\n", track.GetX()));
			continue;
    }

    ULong_t status = seed->GetStatus();
    if((status & AliESDtrack::kTRDout) == 0) continue;
    if((status & AliESDtrack::kTRDin)  != 0) continue;
    nseed++; 

    track.ResetCovariance(50.0);

		// do the propagation and processing
    FollowProlongation(track);
    // computes PID for track
    track.CookPID();
    // update calibration references using this track
		//track.Calibrate();

		// Prolongate to TPC
    Double_t xTPC = 250.0;
    if (PropagateToX(track, xTPC, fgkMaxStep)) { //  -with update
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDrefit);
      track.UpdateESDtrack(seed);
    	// Add TRD track to ESDfriendTrack
			if (AliTRDReconstructor::StreamLevel() > 0 /*&& quality TODO*/){ 
				AliTRDtrackV1 *calibTrack = new AliTRDtrackV1(track);
				calibTrack->SetOwner();
				seed->AddCalibObject(calibTrack);
			}
			found++;
    } else {  // - without update
      AliTRDtrackV1 tt(*seed);
      if (PropagateToX(tt, xTPC, fgkMaxStep)) seed->UpdateTrackParams(&tt, AliESDtrack::kTRDrefit);
    }
  }
  AliInfo(Form("Number of loaded seeds: %d",nseed));
  AliInfo(Form("Number of found tracks from loaded seeds: %d",found));
  
	return 0;
}


//____________________________________________________________________
Int_t AliTRDtrackerV1::FollowProlongation(AliTRDtrackV1 &t)
{
// Extrapolates the TRD track in the TPC direction.
//
// Parameters
//   t : the TRD track which has to be extrapolated
// 
// Output
//   number of clusters attached to the track
//
// Detailed description
//
// Starting from current radial position of track <t> this function
// extrapolates the track through the 6 TRD layers. The following steps
// are being performed for each plane:
// 1. prepare track:
//   a. get plane limits in the local x direction
//   b. check crossing sectors 
//   c. check track inclination
// 2. search tracklet in the tracker list (see GetTracklet() for details)
// 3. evaluate material budget using the geo manager
// 4. propagate and update track using the tracklet information.
//
// Debug level 2
//
  
	//AliInfo("");
	Int_t    nClustersExpected = 0;
	Int_t lastplane = 5; //GetLastPlane(&t);
	for (Int_t iplane = lastplane; iplane >= 0; iplane--) {
		//AliInfo(Form("plane %d", iplane));
    Int_t row1 = GetGlobalTimeBin(0, iplane, 0); // to be modified to the true time bin in the geometrical acceptance
		//AliInfo(Form("row1 %d", row1));

    // Propagate track close to the plane if neccessary
		AliTRDpropagationLayer *layer = fTrSec[0]->GetLayer(row1);
		Double_t currentx  = layer->GetX();
    if (currentx < (-fgkMaxStep + t.GetX())) 
      if (!PropagateToX(t, currentx+fgkMaxStep, fgkMaxStep)) break;

    if (!AdjustSector(&t)) break;
     
		Int_t row0    = GetGlobalTimeBin(0,iplane,GetTimeBinsPerPlane()-1);
		//AliInfo(Form("row0 %d", row0));

    // Start global position
    Double_t xyz0[3];
    t.GetXYZ(xyz0);

		// End global position
    Double_t x = fTrSec[0]->GetLayer(row0)->GetX(), y, z;
    if (!t.GetProlongation(x,y,z)) break;    
    Double_t xyz1[3];
    xyz1[0] =  x * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha());
    xyz1[1] = +x * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;
		
    // Get material budget
    Double_t param[7];
    AliTracker::MeanMaterialBudget(xyz0,xyz1,param);
    Double_t xrho= param[0]*param[4];
    Double_t xx0 = param[1]; // Get mean propagation parameters

    // Propagate and update
    //Int_t sector = t.GetSector();
    Int_t   index   = 0;
		//AliInfo(Form("sector %d", sector));
    AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
	  //AliInfo(Form("tracklet %p @ %d", tracklet, index));
		if(!tracklet) continue;
		//AliInfo(Form("reco %p", tracklet->GetRecoParam()));
		t.SetTracklet(tracklet, iplane, index);
		
		t.PropagateTo(tracklet->GetX0(), xx0, xrho); // not correct
	  if (!AdjustSector(&t)) break;
	  
    Double_t maxChi2 = t.GetPredictedChi2(tracklet);
	  if (maxChi2 < 1e+10 && t.Update(tracklet, maxChi2)){ 
	  	nClustersExpected += tracklet->GetN();
  	}
  }

#ifdef DEBUG
	if(AliTRDReconstructor::StreamLevel() > 1){
		Int_t index;
		for(int iplane=0; iplane<6; iplane++){
			AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
			if(!tracklet) continue;
			t.SetTracklet(tracklet, iplane, index);
		}

		TTreeSRedirector &cstreamer = *fDebugStreamer;
		cstreamer << "FollowProlongation"
			<< "ncl="      << nClustersExpected
			<< "track.="   << &t
			<< "\n";
	}
#endif

  return nClustersExpected;

}

//_____________________________________________________________________________
Int_t AliTRDtrackerV1::FollowBackProlongation(AliTRDtrackV1 &t)
{
// Extrapolates the TRD track in the TOF direction.
//
// Parameters
//   t : the TRD track which has to be extrapolated
// 
// Output
//   number of clusters attached to the track
//
// Detailed description
//
// Starting from current radial position of track <t> this function
// extrapolates the track through the 6 TRD layers. The following steps
// are being performed for each plane:
// 1. prepare track:
//   a. get plane limits in the local x direction
//   b. check crossing sectors 
//   c. check track inclination
// 2. build tracklet (see AliTRDseed::AttachClusters() for details)
// 3. evaluate material budget using the geo manager
// 4. propagate and update track using the tracklet information.
//
// Debug level 2
//

	Int_t nClustersExpected = 0;

  // Loop through the TRD planes
  for (Int_t iplane = 0; iplane < AliTRDgeometry::Nplan(); iplane++) {
		//AliInfo(Form("Processing plane %d ...", iplane));
		// Get the global time bin for the first local time bin of the given plane
    Int_t row0 = GetGlobalTimeBin(0, iplane, fTimeBinsPerPlane-1);

		// Retrive first propagation layer in the chamber
		AliTRDpropagationLayer *layer = fTrSec[0]->GetLayer(row0);

		// Get the X coordinates of the propagation layer for the first time bin
    Double_t currentx = layer->GetX(); // what if X is not defined ???
    if (currentx < t.GetX()) continue;
    
		// Get the global time bin for the last local time bin of the given plane
		Int_t row1 = GetGlobalTimeBin(0, iplane, 0);
    
		// Propagate closer to the current chamber if neccessary 
    if (currentx > (fgkMaxStep + t.GetX()) && !PropagateToX(t, currentx-fgkMaxStep, fgkMaxStep)) break;
   
		// Rotate track to adjacent sector if neccessary
    if (!AdjustSector(&t)) break;
		Int_t sector =  Int_t(TMath::Abs(t.GetAlpha()/AliTRDgeometry::GetAlpha()));
		if(t.GetAlpha() < 0) sector = AliTRDgeometry::Nsect() - sector-1;
		
		//AliInfo(Form("sector %d [%f]", sector, t.GetAlpha()));

		// Check whether azimuthal angle is getting too large
    if (TMath::Abs(t.GetSnp()) > fgkMaxSnp) break;
   
    //Calculate global entry and exit positions of the track in chamber (only track prolongation)
    Double_t xyz0[3]; // entry point 
		t.GetXYZ(xyz0);   
		//printf("Entry global x[%7.3f] y[%7.3f] z[%7.3f]\n", xyz0[0], xyz0[1], xyz0[2]);
				
		// Get local Y and Z at the X-position of the end of the chamber
		Double_t x0 = fTrSec[sector]->GetLayer(row1)->GetX(), y, z;
		if (!t.GetProlongation(x0, y, z)) break;
		//printf("Exit  local x[%7.3f] y[%7.3f] z[%7.3f]\n", x0, y, z);

		Double_t xyz1[3]; // exit point
		xyz1[0] =  x0 * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha()); 
    xyz1[1] = +x0 * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;

		//printf("Exit  global x[%7.3f] y[%7.3f] z[%7.3f]\n", xyz1[0], xyz1[1], xyz1[2]);
		// Find tracklet along the path inside the chamber
    AliTRDseedV1 tracklet(*t.GetTracklet(iplane));
		// if the track is not already build (e.g. stand alone tracker) we build it now.
		if(!tracklet.GetN()){ // a better check has to be implemented TODO!!!!!!!
		
			//AliInfo(Form("Building tracklet for plane %d ...", iplane));
			// check if we are inside detection volume
			Int_t ichmb = fGeom->GetChamber(xyz0[2], iplane); 
			if(ichmb<0) ichmb = fGeom->GetChamber(xyz1[2], iplane);
			if(ichmb<0){
				// here we should decide what to do with the track. The space between the pads in 2 chambers is 4cm+. Is it making sense to continue building the tracklet here TODO????
				AliWarning(Form("Track prolongated in the interspace between TRD detectors in plane %d. Skip plane. To be fixed !", iplane));
				continue;
			}
	
			// temporary until the functionalities of AliTRDpropagationLayer and AliTRDstackLayer are merged TODO
			AliTRDpadPlane *pp = fGeom->GetPadPlane(iplane, ichmb);
			Int_t nrows = pp->GetNrows();
			Double_t stacklength = pp->GetRow0ROC() - pp->GetRowEndROC();/*(nrows - 2) * pp->GetLengthIPad()	+ 2 * pp->GetLengthOPad() + (nrows - 1) * pp->GetRowSpacing();*/
			Double_t z0  = fGeom->GetRow0(iplane, ichmb, 0);

			Int_t nClustersChmb = 0;
			AliTRDstackLayer stackLayer[35];
			for(int itb=0; itb<fTimeBinsPerPlane; itb++){
				const AliTRDpropagationLayer ksmLayer(*(fTrSec[sector]->GetLayer(row1 - itb)));
				stackLayer[itb] = ksmLayer;
#ifdef DEBUG
				stackLayer[itb].SetDebugStream(fDebugStreamer);
#endif			
				stackLayer[itb].SetRange(z0 - stacklength, stacklength);
				stackLayer[itb].SetSector(sector);
				stackLayer[itb].SetStackNr(ichmb);
				stackLayer[itb].SetNRows(nrows);
				stackLayer[itb].SetRecoParam(fRecoParam);
				stackLayer[itb].BuildIndices();
				nClustersChmb += stackLayer[itb].GetNClusters();
			}
			//AliInfo(Form("Detector p[%d] c[%d]. Building tracklet from %d clusters ... ", iplane, ichmb, nClustersChmb));

			tracklet.SetRecoParam(fRecoParam);
			tracklet.SetTilt(TMath::Tan(-TMath::DegToRad()*pp->GetTiltingAngle()));
			tracklet.SetPadLength(pp->GetLengthIPad());
			tracklet.SetPlane(iplane);
			Int_t tbRange   = fTimeBinsPerPlane; //Int_t(AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght() * AliTRDCommonParam::Instance()->GetSamplingFrequency()/AliTRDcalibDB::Instance()->GetVdriftDet()->GetValue(det));
			//printf("%d hl[%f] pl[%f] tb[%d]\n", il, hL[il], padlength[il], tbRange[il]);
			tracklet.SetNTimeBinsRange(tbRange);
			tracklet.SetX0(x0);
			tracklet.Init(&t);
			if(!tracklet.AttachClustersIter(stackLayer, 1000.)) continue;

			//if(!tracklet.AttachClusters(stackLayer, kTRUE)) continue;
			//if(!tracklet.Fit()) continue;
		}
		Int_t ncl = tracklet.GetN();
		//AliInfo(Form("N clusters %d", ncl));
		
		// Discard tracklet if bad quality.
		//Check if this condition is not already checked during building of the tracklet
    if(ncl < fTimeBinsPerPlane * fRecoParam->GetFindableClusters()){
			//AliInfo(Form("Discard tracklet for %d nclusters", ncl));
			continue;
		}
		
		// load tracklet to the tracker and the track
		Int_t index = SetTracklet(&tracklet);
		t.SetTracklet(&tracklet, iplane, index);
		
		// Calculate the mean material budget along the path inside the chamber
    Double_t param[7];
		AliTracker::MeanMaterialBudget(xyz0, xyz1, param);	
    // The mean propagation parameters
    Double_t xrho = param[0]*param[4]; // density*length
    Double_t xx0  = param[1]; // radiation length
		
		// Propagate and update track
		t.PropagateTo(tracklet.GetX0(), xx0, xrho);
	  if (!AdjustSector(&t)) break;
		Double_t maxChi2 = t.GetPredictedChi2(&tracklet);
		if (maxChi2<1e+10 && t.Update(&tracklet, maxChi2)){ 
			nClustersExpected += ncl;
		}
		// Reset material budget if 2 consecutive gold
		if(iplane>0 && ncl + t.GetTracklet(iplane-1)->GetN() > 20) t.SetBudget(2, 0.);

		// Make backup of the track until is gold
		// TO DO update quality check of the track.
		// consider comparison with fTimeBinsRange
		Float_t ratio0 = ncl / Float_t(fTimeBinsPerPlane);
		//Float_t ratio1 = Float_t(t.GetNumberOfClusters()+1) / Float_t(t.GetNExpected()+1);	
    //printf("tracklet.GetChi2() %f     [< 18.0]\n", tracklet.GetChi2()); 
		//printf("ratio0    %f              [>   0.8]\n", ratio0);
		//printf("ratio1     %f             [>   0.6]\n", ratio1); 
		//printf("ratio0+ratio1 %f          [>   1.5]\n", ratio0+ratio1); 
		//printf("t.GetNCross()  %d         [==    0]\n", t.GetNCross()); 
		//printf("TMath::Abs(t.GetSnp()) %f [<  0.85]\n", TMath::Abs(t.GetSnp()));
		//printf("t.GetNumberOfClusters() %d [>    20]\n", t.GetNumberOfClusters());
    
		if (//(tracklet.GetChi2()      <  18.0) && TO DO check with FindClusters and move it to AliTRDseed::Update 
        (ratio0                  >   0.8) && 
        //(ratio1                  >   0.6) && 
        //(ratio0+ratio1           >   1.5) && 
        (t.GetNCross()           ==    0) && 
        (TMath::Abs(t.GetSnp())  <  0.85) &&
        (t.GetNumberOfClusters() >    20)) t.MakeBackupTrack();
		
	} // end planes loop

#ifdef DEBUG
	if(AliTRDReconstructor::StreamLevel() > 1){
		TTreeSRedirector &cstreamer = *fDebugStreamer;
		cstreamer << "FollowBackProlongation"
			<< "ncl="      << nClustersExpected
			<< "track.="   << &t
			<< "\n";
	}
#endif
	
	return nClustersExpected;
}

//____________________________________________________________________
void AliTRDtrackerV1::UnloadClusters() 
{ 
  //
  // Clears the arrays of clusters and tracks. Resets sectors and timebins 
  //

  Int_t i;
  Int_t nentr;

  nentr = fClusters->GetEntriesFast();
	//AliInfo(Form("clearing %d clusters", nentr));
  for (i = 0; i < nentr; i++) {
    delete fClusters->RemoveAt(i);
  }
  fNclusters = 0;

  nentr = fTracklets->GetEntriesFast();
	//AliInfo(Form("clearing %d tracklets", nentr));
  for (i = 0; i < nentr; i++) {
    delete fTracklets->RemoveAt(i);
  }

  nentr = fSeeds->GetEntriesFast();
	//AliInfo(Form("clearing %d seeds", nentr));
  for (i = 0; i < nentr; i++) {
    delete fSeeds->RemoveAt(i);
  }

  nentr = fTracks->GetEntriesFast();
  //AliInfo(Form("clearing %d tracks", nentr));
	for (i = 0; i < nentr; i++) {
    delete fTracks->RemoveAt(i);
  }

  Int_t nsec = AliTRDgeometry::kNsect;
  for (i = 0; i < nsec; i++) {    
    for(Int_t pl = 0; pl < fTrSec[i]->GetNumberOfLayers(); pl++) {
      fTrSec[i]->GetLayer(pl)->Clear();
    }
  }

}

//____________________________________________________________________
AliTRDseedV1* AliTRDtrackerV1::GetTracklet(AliTRDtrackV1 *track, Int_t p, Int_t &idx)
{
// Find tracklet for TRD track <track>
// Parameters
// - track
// - sector
// - plane
// - index
// Output
// tracklet
// index
// Detailed description
//
	idx = track->GetTrackletIndex(p);
	//AliInfo(Form("looking for tracklet in plane %d idx %d [%d]", p, idx, track->GetTrackletIndex(p)));
	AliTRDseedV1 *tracklet = idx<0 ? 0x0 : (AliTRDseedV1*)fTracklets->UncheckedAt(idx);
	//AliInfo(Form("found 0x%x @ %d", tracklet, idx));

//   Int_t *index = track->GetTrackletIndexes();
//   for (UInt_t i = 0; i < 6; i++) AliInfo(Form("index[%d] = %d", i, index[i]));
// 
// 	for (UInt_t i = 0; i < 6/*kMaxTimeBinIndex*/; i++) {
//     if (index[i] < 0) continue;
// 
// 		tracklet = (AliTRDseedV1*)fTracklets->UncheckedAt(index[i]);
// 		if(!tracklet) break;
// 
// 		if(tracklet->GetPlane() != p) continue;
// 
// 		idx = index[i];
// 	}

	return tracklet;
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::SetTracklet(AliTRDseedV1 *tracklet)
{
// Add this tracklet to the list of tracklets stored in the tracker
//
// Parameters
//   - tracklet : pointer to the tracklet to be added to the list
//
// Output
//   - the index of the new tracklet in the tracker tracklets list
//
// Detailed description
// Build the tracklets list if it is not yet created (late initialization)
// and adds the new tracklet to the list.
//
	if(!fTracklets){
		fTracklets = new TClonesArray("AliTRDseedV1", AliTRDgeometry::Nsect()*kMaxTracksStack);
		fTracklets->SetOwner(kTRUE);
	}
	Int_t nentries = fTracklets->GetEntriesFast();
	AliTRDseedV1 *t = new ((*fTracklets)[nentries]) AliTRDseedV1(*tracklet);
	//AliInfo(Form("0x%x @ %d", t, nentries));
	return nentries;
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::Clusters2TracksSM(AliTRDtracker::AliTRDtrackingSector *sector
                                       , AliESDEvent *esd)
{
  //
  // Steer tracking for one SM.
  //
  // Parameters :
  //   sector  : Array of (SM) propagation layers containing clusters
  //   esd     : The current ESD event. On output it contains the also
  //             the ESD (TRD) tracks found in this SM. 
  //
  // Output :
  //   Number of tracks found in this TRD supermodule.
  // 
  // Detailed description
  //
  // 1. Unpack AliTRDpropagationLayers objects for each stack.
  // 2. Launch stack tracking. 
  //    See AliTRDtrackerV1::Clusters2TracksStack() for details.
  // 3. Pack results in the ESD event.
  //
	
	AliTRDpadPlane *pp = 0x0;
	
	// allocate space for esd tracks in this SM
	TClonesArray esdTrackList("AliESDtrack", 2*kMaxTracksStack);
	esdTrackList.SetOwner();
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	const Int_t kFindable = Int_t(fRecoParam->GetFindableClusters()*6.*nTimeBins);
	
	Int_t ntracks = 0;
	Int_t nClStack = 0;
	for(int istack = 0; istack<AliTRDpropagationLayer::kZones; istack++){
		AliTRDstackLayer stackLayer[kNPlanes*kNTimeBins];
		
		nClStack = 0;
		//AliInfo(Form("Processing stack %i ...",istack));
		//AliInfo("Building stack propagation layers ...");
		for(int ilayer=0; ilayer<kNPlanes*nTimeBins; ilayer++){
			pp = fGeom->GetPadPlane((Int_t)(ilayer/nTimeBins), istack);
			Double_t stacklength = (pp->GetNrows() - 2) * pp->GetLengthIPad() 
                                             + 2 * pp->GetLengthOPad() + 2 * pp->GetLengthRim();
			//Debug
			Double_t z0  = fGeom->GetRow0((Int_t)(ilayer/nTimeBins),istack,0);
			const AliTRDpropagationLayer ksmLayer(*(sector->GetLayer(ilayer)));
			stackLayer[ilayer] = ksmLayer;
#ifdef DEBUG
			stackLayer[ilayer].SetDebugStream(fDebugStreamer);
#endif			
			stackLayer[ilayer].SetRange(z0 - stacklength, stacklength);
			stackLayer[ilayer].SetSector(sector->GetSector());
			stackLayer[ilayer].SetStackNr(istack);
			stackLayer[ilayer].SetNRows(pp->GetNrows());
			stackLayer[ilayer].SetRecoParam(fRecoParam);
			stackLayer[ilayer].BuildIndices();
			nClStack += stackLayer[ilayer].GetNClusters();
		}
		//AliInfo(Form("Finish building stack propagation layers. nClusters %d.", nClStack));
		if(nClStack < kFindable) continue;
		ntracks += Clusters2TracksStack(&stackLayer[0], &esdTrackList);
	}
	//AliInfo(Form("Found %d tracks in SM", ntracks));
	
	for(int itrack=0; itrack<ntracks; itrack++) 
          esd->AddTrack((AliESDtrack*)esdTrackList[itrack]);

	return ntracks;
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::Clusters2TracksStack(AliTRDstackLayer *layer
                                          , TClonesArray *esdTrackList)
{
  //
  // Make tracks in one TRD stack.
  //
  // Parameters :
  //   layer  : Array of stack propagation layers containing clusters
  //   esdTrackList  : Array of ESD tracks found by the stand alone tracker. 
  //                   On exit the tracks found in this stack are appended.
  //
  // Output :
  //   Number of tracks found in this stack.
  // 
  // Detailed description
  //
  // 1. Find the 3 most useful seeding chambers. See BuildSeedingConfigs() for details.
  // 2. Steer AliTRDtrackerV1::MakeSeeds() for 3 seeding layer configurations. 
  //    See AliTRDtrackerV1::MakeSeeds() for more details.
  // 3. Arrange track candidates in decreasing order of their quality
  // 4. Classify tracks in 5 categories according to:
  //    a) number of layers crossed
  //    b) track quality 
  // 5. Sign clusters by tracks in decreasing order of track quality
  // 6. Build AliTRDtrack out of seeding tracklets
  // 7. Cook MC label
  // 8. Build ESD track and register it to the output list
  //

	AliTRDseedV1 sseed[kMaxTracksStack*6]; // to be initialized
	Int_t pars[4]; // MakeSeeds parameters

	//Double_t alpha = AliTRDgeometry::GetAlpha();
	//Double_t shift = .5 * alpha;
	Int_t configs[kNConfigs];
	
	// Build initial seeding configurations
	Double_t quality = BuildSeedingConfigs(layer, configs);
#ifdef DEBUG
		if(AliTRDReconstructor::StreamLevel() > 1) 
                  AliInfo(Form("Plane config %d %d %d Quality %f"
                              , configs[0], configs[1], configs[2], quality));
#endif

	// Initialize contors
	Int_t ntracks,      // number of TRD track candidates
	      ntracks1,     // number of registered TRD tracks/iter
	      ntracks2 = 0; // number of all registered TRD tracks in stack
	fSieveSeeding = 0;
	do{
		// Loop over seeding configurations
		ntracks = 0; ntracks1 = 0;
		for (Int_t iconf = 0; iconf<3; iconf++) {
			pars[0] = configs[iconf];
			pars[1] = layer->GetStackNr();
			pars[2] = ntracks;
			ntracks = MakeSeeds(layer, &sseed[6*ntracks], pars);
			if(ntracks == kMaxTracksStack) break;
		}
#ifdef DEBUG
		if(AliTRDReconstructor::StreamLevel() > 1) AliInfo(Form("Candidate TRD tracks %d in stack %d iteration %d.", ntracks, pars[1], fSieveSeeding));
#endif		
		if(!ntracks) break;
		
		// Sort the seeds according to their quality
		Int_t sort[kMaxTracksStack];
		TMath::Sort(ntracks, fTrackQuality, sort, kTRUE);
	
		// Initialize number of tracks so far and logic switches
		Int_t ntracks0 = esdTrackList->GetEntriesFast();
		Bool_t signedTrack[kMaxTracksStack];
		Bool_t fakeTrack[kMaxTracksStack];
		for (Int_t i=0; i<ntracks; i++){
			signedTrack[i] = kFALSE;
			fakeTrack[i] = kFALSE;
		}
		//AliInfo("Selecting track candidates ...");
		
		// Sieve clusters in decreasing order of track quality
		Double_t trackParams[7];
// 		AliTRDseedV1 *lseed = 0x0;
		Int_t jSieve = 0, candidates;
		do{
			//AliInfo(Form("\t\tITER = %i ", jSieve));

			// Check track candidates
			candidates = 0;
			for (Int_t itrack = 0; itrack < ntracks; itrack++) {
				Int_t trackIndex = sort[itrack];
				if (signedTrack[trackIndex] || fakeTrack[trackIndex]) continue;
	
				
				// Calculate track parameters from tracklets seeds
				Int_t labelsall[1000];
				Int_t nlabelsall = 0;
				Int_t naccepted  = 0;
				Int_t ncl        = 0;
				Int_t nused      = 0;
				Int_t nlayers    = 0;
				Int_t findable   = 0;
				for (Int_t jLayer = 0; jLayer < kNPlanes; jLayer++) {
					Int_t jseed = kNPlanes*trackIndex+jLayer;
					if (TMath::Abs(sseed[jseed].GetYref(0) / sseed[jseed].GetX0()) < 0.15) 
                                          findable++;
	
					if(!sseed[jseed].IsOK()) continue;
					sseed[jseed].UpdateUsed();
					ncl   += sseed[jseed].GetN2();
					nused += sseed[jseed].GetNUsed();
					nlayers++;
	
					// Cooking label
					for (Int_t itime = 0; itime < fTimeBinsPerPlane; itime++) {
						if(!sseed[jseed].IsUsable(itime)) continue;
						naccepted++;
						Int_t tindex = 0, ilab = 0;
						while(ilab<3 && (tindex = sseed[jseed].GetClusters(itime)->GetLabel(ilab)) >= 0){
							labelsall[nlabelsall++] = tindex;
							ilab++;
						}
					}
				}
				// Filter duplicated tracks
				if (nused > 30){
					printf("Skip %d nused %d\n", trackIndex, nused);
					fakeTrack[trackIndex] = kTRUE;
					continue;
				}
				if (Float_t(nused)/ncl >= .25){
					printf("Skip %d nused/ncl >= .25\n", trackIndex);
					fakeTrack[trackIndex] = kTRUE;
				 	continue;
				}
				
				// Classify tracks
				Bool_t skip = kFALSE;
				switch(jSieve){
				case 0:
					if(nlayers < 6) {skip = kTRUE; break;}
					if(TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -5.){skip = kTRUE; break;}
					break;
	
				case 1:
					if(nlayers < findable){skip = kTRUE; break;}
					if(TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -4.){skip = kTRUE; break;}
					break;
	
				case 2:
					if ((nlayers == findable) || (nlayers == 6)) { skip = kTRUE; break;}
					if (TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -6.0){skip = kTRUE; break;}
					break;
	
				case 3:
					if (TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -5.){skip = kTRUE; break;}
					break;
	
				case 4:
					if (nlayers == 3){skip = kTRUE; break;}
					//if (TMath::Log(1.E-9+fTrackQuality[trackIndex]) - nused/(nlayers-3.0) < -15.0){skip = kTRUE; break;}
					break;
				}
				if(skip){
					candidates++;
					printf("REJECTED : %d [%d] nlayers %d trackQuality = %e nused %d\n", itrack, trackIndex, nlayers, fTrackQuality[trackIndex], nused);
					continue;
				}
				signedTrack[trackIndex] = kTRUE;
						

				// Build track label - what happens if measured data ???
				Int_t labels[1000];
				Int_t outlab[1000];
				Int_t nlab = 0;
				for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
					Int_t jseed = kNPlanes*trackIndex+iLayer;
					if(!sseed[jseed].IsOK()) continue;
					for(int ilab=0; ilab<2; ilab++){
						if(sseed[jseed].GetLabels(ilab) < 0) continue;
						labels[nlab] = sseed[jseed].GetLabels(ilab);
						nlab++;
					}
				}
				Freq(nlab,labels,outlab,kFALSE);
				Int_t   label     = outlab[0];
				Int_t   frequency = outlab[1];
				Freq(nlabelsall,labelsall,outlab,kFALSE);
				Int_t   label1    = outlab[0];
				Int_t   label2    = outlab[2];
				Float_t fakeratio = (naccepted - outlab[1]) / Float_t(naccepted);
	
				
				// Sign clusters
				AliTRDcluster *cl = 0x0; Int_t clusterIndex = -1;
				for (Int_t jLayer = 0; jLayer < 6; jLayer++) {
					Int_t jseed = kNPlanes*trackIndex+jLayer;
					if(!sseed[jseed].IsOK()) continue;
					if(TMath::Abs(sseed[jseed].GetYfit(1) - sseed[jseed].GetYfit(1)) >= .2) continue; // check this condition with Marian
					sseed[jseed].UseClusters();
					if(!cl){
						Int_t ic = 0;
						while(!(cl = sseed[jseed].GetClusters(ic))) ic++;
						clusterIndex =  sseed[jseed].GetIndexes(ic);
					}
				}
				if(!cl) continue;

				
				// Build track parameters
				AliTRDseedV1 *lseed =&sseed[trackIndex*6];
				Int_t idx = 0;
				while(idx<3 && !lseed->IsOK()) {
					idx++;
					lseed++;
				}
				Double_t cR = lseed->GetC();
				trackParams[1] = lseed->GetYref(0);
				trackParams[2] = lseed->GetZref(0);
				trackParams[3] = lseed->GetX0() * cR - TMath::Sin(TMath::ATan(lseed->GetYref(1)));
				trackParams[4] = lseed->GetZref(1) / TMath::Sqrt(1. + lseed->GetYref(1) * lseed->GetYref(1));
				trackParams[5] = cR;
				trackParams[0] = lseed->GetX0();
				trackParams[6] = layer[0].GetSector();/* *alpha+shift;	// Supermodule*/

#ifdef DEBUG
				if(AliTRDReconstructor::StreamLevel() > 1) printf("Track %d [%d] nlayers %d trackQuality = %e nused %d, yref = %3.3f\n", itrack, trackIndex, nlayers, fTrackQuality[trackIndex], nused, trackParams[1]);
				
				if(AliTRDReconstructor::StreamLevel() >= 1){
					Int_t sector = layer[0].GetSector();
					Int_t nclusters = 0;
					AliTRDseedV1 *dseed[6];
					for(int is=0; is<6; is++){
						dseed[is] = new AliTRDseedV1(sseed[trackIndex*6+is]);
						dseed[is]->SetOwner();
						nclusters += sseed[is].GetN2();
						//for(int ic=0; ic<30; ic++) if(sseed[trackIndex*6+is].GetClusters(ic)) printf("l[%d] tb[%d] cptr[%p]\n", is, ic, sseed[trackIndex*6+is].GetClusters(ic));
					}
					//Int_t eventNrInFile = esd->GetEventNumberInFile();
					//AliInfo(Form("Number of clusters %d.", nclusters));
					TTreeSRedirector &cstreamer = *fDebugStreamer;
					cstreamer << "Clusters2TracksStack"
						<< "Iter="      << fSieveSeeding
						<< "Like="      << fTrackQuality[trackIndex]
						<< "S0.="       << dseed[0]
						<< "S1.="       << dseed[1]
						<< "S2.="       << dseed[2]
						<< "S3.="       << dseed[3]
						<< "S4.="       << dseed[4]
						<< "S5.="       << dseed[5]
						<< "p0=" << trackParams[0]
						<< "p1=" << trackParams[1]
						<< "p2=" << trackParams[2]
						<< "p3=" << trackParams[3]
						<< "p4=" << trackParams[4]
						<< "p5=" << trackParams[5]
						<< "p6=" << trackParams[6]
						<< "Sector="    << sector
						<< "Stack="     << pars[1]
						<< "Label="     << label
						<< "Label1="    << label1
						<< "Label2="    << label2
						<< "FakeRatio=" << fakeratio
						<< "Freq="      << frequency
						<< "Ncl="       << ncl
						<< "NLayers="   << nlayers
						<< "Findable="  << findable
						<< "NUsed="     << nused
						<< "\n";
					//???for(int is=0; is<6; is++) delete dseed[is];
				}
#endif
			
				AliTRDtrackV1 *track = AliTRDtrackerV1::MakeTrack(&sseed[trackIndex*kNPlanes], trackParams);
				if(!track){
					AliWarning("Fail to build a TRD Track.");
					continue;
				}
				AliInfo("End of MakeTrack()");
				AliESDtrack esdTrack;
				esdTrack.UpdateTrackParams(track, AliESDtrack::kTRDout);
				esdTrack.SetLabel(track->GetLabel());
				new ((*esdTrackList)[ntracks0++]) AliESDtrack(esdTrack);
				ntracks1++;
			}

			jSieve++;
		} while(jSieve<5 && candidates); // end track candidates sieve
		if(!ntracks1) break;

		// increment counters
		ntracks2 += ntracks1;
		fSieveSeeding++;

		// Rebuild plane configurations and indices taking only unused clusters into account
		quality = BuildSeedingConfigs(layer, configs);
		//if(quality < fRecoParam->GetPlaneQualityThreshold()) break;
		
		for(Int_t il = 0; il < kNPlanes * fTimeBinsPerPlane; il++) layer[il].BuildIndices(fSieveSeeding);

#ifdef DEBUG
				if(AliTRDReconstructor::StreamLevel() > 1) AliInfo(Form("Sieve level %d Plane config %d %d %d Quality %f", fSieveSeeding, configs[0], configs[1], configs[2], quality));
#endif
	} while(fSieveSeeding<10); // end stack clusters sieve
	


	//AliInfo(Form("Registered TRD tracks %d in stack %d.", ntracks2, pars[1]));

	return ntracks2;
}

//___________________________________________________________________
Double_t AliTRDtrackerV1::BuildSeedingConfigs(AliTRDstackLayer *layers
                                            , Int_t *configs)
{
  //
  // Assign probabilities to chambers according to their
  // capability of producing seeds.
  // 
  // Parameters :
  //
  //   layers : Array of stack propagation layers for all 6 chambers in one stack
  //   configs : On exit array of configuration indexes (see GetSeedingConfig()
  // for details) in the decreasing order of their seeding probabilities. 
  //
  // Output :
  //
  //  Return top configuration quality 
  //
  // Detailed description:
  //
  // To each chamber seeding configuration (see GetSeedingConfig() for
  // the list of all configurations) one defines 2 quality factors:
  //  - an apriori topological quality (see GetSeedingConfig() for details) and
  //  - a data quality based on the uniformity of the distribution of
  //    clusters over the x range (time bins population). See CookChamberQA() for details.
  // The overall chamber quality is given by the product of this 2 contributions.
  // 

	Double_t chamberQA[kNPlanes];
	for(int iplane=0; iplane<kNPlanes; iplane++){
		chamberQA[iplane] = MakeSeedingPlanes(&layers[iplane*fTimeBinsPerPlane]);
		//printf("chamberQA[%d] = %f\n", iplane, chamberQA[iplane]);
	}

	Double_t tconfig[kNConfigs];
	Int_t planes[4];
	for(int iconf=0; iconf<kNConfigs; iconf++){
		GetSeedingConfig(iconf, planes);
		tconfig[iconf] = fgTopologicQA[iconf];
		for(int iplane=0; iplane<4; iplane++) tconfig[iconf] *= chamberQA[planes[iplane]]; 
	}
	
	TMath::Sort(kNConfigs, tconfig, configs, kTRUE);
	return tconfig[configs[0]];
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::MakeSeeds(AliTRDstackLayer *layers
                               , AliTRDseedV1 *sseed
                               , Int_t *ipar)
{
  //
  // Make tracklet seeds in the TRD stack.
  //
  // Parameters :
  //   layers : Array of stack propagation layers containing clusters
  //   sseed  : Array of empty tracklet seeds. On exit they are filled.
  //   ipar   : Control parameters:
  //       ipar[0] -> seeding chambers configuration
  //       ipar[1] -> stack index
  //       ipar[2] -> number of track candidates found so far
  //
  // Output :
  //   Number of tracks candidates found.
  // 
  // Detailed description
  //
  // The following steps are performed:
  // 1. Select seeding layers from seeding chambers
  // 2. Select seeding clusters from the seeding AliTRDpropagationLayerStack.
  //   The clusters are taken from layer 3, layer 0, layer 1 and layer 2, in
  //   this order. The parameters controling the range of accepted clusters in
  //   layer 0, 1, and 2 are defined in AliTRDstackLayer::BuildCond().
  // 3. Helix fit of the cluster set. (see AliTRDtrackerFitter::FitRieman(AliTRDcluster**))
  // 4. Initialize seeding tracklets in the seeding chambers.
  // 5. Filter 0.
  //   Chi2 in the Y direction less than threshold ... (1./(3. - sLayer))
  //   Chi2 in the Z direction less than threshold ... (1./(3. - sLayer))
  // 6. Attach clusters to seeding tracklets and find linear approximation of
  //   the tracklet (see AliTRDseedV1::AttachClustersIter()). The number of used
  //   clusters used by current seeds should not exceed ... (25).
  // 7. Filter 1.
  //   All 4 seeding tracklets should be correctly constructed (see
  //   AliTRDseedV1::AttachClustersIter())
  // 8. Helix fit of the seeding tracklets
  // 9. Filter 2.
  //   Likelihood calculation of the fit. (See AliTRDtrackerV1::CookLikelihood() for details)
  // 10. Extrapolation of the helix fit to the other 2 chambers:
  //    a) Initialization of extrapolation tracklet with fit parameters
  //    b) Helix fit of tracklets
  //    c) Attach clusters and linear interpolation to extrapolated tracklets
  //    d) Helix fit of tracklets
  // 11. Improve seeding tracklets quality by reassigning clusters.
  //      See AliTRDtrackerV1::ImproveSeedQuality() for details.
  // 12. Helix fit of all 6 seeding tracklets and chi2 calculation
  // 13. Hyperplane fit and track quality calculation. See AliTRDtrackerFitter::FitHyperplane() for details.
  // 14. Cooking labels for tracklets. Should be done only for MC
  // 15. Register seeds.
  //

	AliTRDcluster *c[4] = {0x0, 0x0, 0x0, 0x0}; // initilize seeding clusters
	AliTRDseedV1 *cseed = &sseed[0]; // initialize tracklets for first track
	Int_t ncl, mcl; // working variable for looping over clusters
	Int_t index[AliTRDstackLayer::kMaxClustersLayer], jndex[AliTRDstackLayer::kMaxClustersLayer];
	// chi2 storage
	// chi2[0] = tracklet chi2 on the Z direction
	// chi2[1] = tracklet chi2 on the R direction
	Double_t chi2[4];


	// this should be data member of AliTRDtrack
	Double_t seedQuality[kMaxTracksStack];
	
	// unpack control parameters
	Int_t config  = ipar[0];
	Int_t istack  = ipar[1];
	Int_t ntracks = ipar[2];
	Int_t planes[kNSeedPlanes]; GetSeedingConfig(config, planes);	
#ifdef DEBUG
		if(AliTRDReconstructor::StreamLevel() > 1) AliInfo(Form("Making seeds Stack[%d] Config[%d] Tracks[%d]...", istack, config, ntracks));
#endif
	
	// Init chambers geometry
	Int_t det, tbRange[6]; // time bins inside the detector geometry
	Double_t hL[kNPlanes];       // Tilting angle
	Float_t padlength[kNPlanes]; // pad lenghts
	AliTRDpadPlane *pp;
	for(int il=0; il<kNPlanes; il++){
		pp = fGeom->GetPadPlane(il, istack); // istack has to be imported
		hL[il]        = TMath::Tan(-TMath::DegToRad()*pp->GetTiltingAngle());
		padlength[il] = pp->GetLengthIPad();
		det           = il; // to be fixed !!!!!
		tbRange[il]   = fTimeBinsPerPlane; //Int_t(AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght() * AliTRDCommonParam::Instance()->GetSamplingFrequency()/AliTRDcalibDB::Instance()->GetVdriftDet()->GetValue(det));
		//printf("%d hl[%f] pl[%f] tb[%d]\n", il, hL[il], padlength[il], tbRange[il]);
	}

	Double_t cond0[4], cond1[4], cond2[4];
	// make seeding layers (to be moved in Clusters2TracksStack)
	AliTRDstackLayer *layer[] = {0x0, 0x0, 0x0, 0x0};
	for(int isl=0; isl<kNSeedPlanes; isl++) layer[isl] = MakeSeedingLayer(&layers[planes[isl] * fTimeBinsPerPlane], planes[isl]);


	// Start finding seeds
	Int_t icl = 0;
	while((c[3] = (*layer[3])[icl++])){
		if(!c[3]) continue;
		layer[0]->BuildCond(c[3], cond0, 0);
		layer[0]->GetClusters(cond0, index, ncl);
		Int_t jcl = 0;
		while(jcl<ncl) {
			c[0] = (*layer[0])[index[jcl++]];
			if(!c[0]) continue;
			Double_t dx    = c[3]->GetX() - c[0]->GetX();
			Double_t theta = (c[3]->GetZ() - c[0]->GetZ())/dx;
			Double_t phi   = (c[3]->GetY() - c[0]->GetY())/dx;
			layer[1]->BuildCond(c[0], cond1, 1, theta, phi);
			layer[1]->GetClusters(cond1, jndex, mcl);

			Int_t kcl = 0;
			while(kcl<mcl) {
				c[1] = (*layer[1])[jndex[kcl++]];
				if(!c[1]) continue;
				layer[2]->BuildCond(c[1], cond2, 2, theta, phi);
				c[2] = layer[2]->GetNearestCluster(cond2);
				if(!c[2]) continue;
				
				//AliInfo("Seeding clusters found. Building seeds ...");
				//for(Int_t i = 0; i < kNSeedPlanes; i++) printf("%i. coordinates: x = %3.3f, y = %3.3f, z = %3.3f\n", i, c[i]->GetX(), c[i]->GetY(), c[i]->GetZ());
				for (Int_t il = 0; il < 6; il++) cseed[il].Reset();

				fFitter->Reset();

				fFitter->FitRieman(c, kNSeedPlanes);

				chi2[0] = 0.; chi2[1] = 0.;
				AliTRDseedV1 *tseed = 0x0;
				for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
					Int_t jLayer = planes[iLayer];
					tseed = &cseed[jLayer];
					tseed->SetRecoParam(fRecoParam);
					tseed->SetPlane(jLayer);
					tseed->SetTilt(hL[jLayer]);
					tseed->SetPadLength(padlength[jLayer]);
					tseed->SetNTimeBinsRange(tbRange[jLayer]);
					tseed->SetX0(layer[iLayer]->GetX());//layers[jLayer*fTimeBinsPerPlane].GetX());

					tseed->Init(fFitter->GetRiemanFitter());
					// temporary until new AttachClusters()
					tseed->SetX0(layers[(jLayer+1)*fTimeBinsPerPlane-1].GetX());
					chi2[0] += tseed->GetChi2Z(c[iLayer]->GetZ());
					chi2[1] += tseed->GetChi2Y(c[iLayer]->GetY());
				}

				Bool_t isFake = kFALSE;
				if (c[0]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
				if (c[1]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
				if (c[2]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
#ifdef DEBUG
				if(AliTRDReconstructor::StreamLevel() >= 2){
					Float_t yref[4], ycluster[4];
					for(int il=0; il<4; il++){
						tseed = &cseed[planes[il]];
						yref[il] = tseed->GetYref(0);
						ycluster[il] = c[il]->GetY();
					}
					Float_t threshold = .5;//1./(3. - sLayer);
					Int_t ll = c[3]->GetLabel(0);
					TTreeSRedirector &cs0 = *fDebugStreamer;
							cs0 << "MakeSeeds0"
							<<"isFake=" << isFake
							<<"label=" << ll
							<<"threshold=" << threshold
							<<"chi2=" << chi2[1]
							<<"yref0="<<yref[0]
							<<"yref1="<<yref[1]
							<<"yref2="<<yref[2]
							<<"yref3="<<yref[3]
							<<"ycluster0="<<ycluster[0]
							<<"ycluster1="<<ycluster[1]
							<<"ycluster2="<<ycluster[2]
							<<"ycluster3="<<ycluster[3]
							<<"\n";
				}
#endif

				if(chi2[0] > fRecoParam->GetChi2Z()/*7./(3. - sLayer)*//*iter*/){
					//AliInfo(Form("Failed chi2 filter on chi2Z [%f].", chi2[0]));
					continue;
				}
				if(chi2[1] > fRecoParam->GetChi2Y()/*1./(3. - sLayer)*//*iter*/){
					//AliInfo(Form("Failed chi2 filter on chi2Y [%f].", chi2[1]));
					continue;
				}
				//AliInfo("Passed chi2 filter.");

#ifdef DEBUG
				if(AliTRDReconstructor::StreamLevel() >= 2){
					Float_t minmax[2] = { -100.0,  100.0 };
					for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
						Float_t max = c[iLayer]->GetZ() + cseed[planes[iLayer]].GetPadLength() * 0.5 + 1.0 - cseed[planes[iLayer]].GetZref(0);
						if (max < minmax[1]) minmax[1] = max;
						Float_t min = c[iLayer]->GetZ()-cseed[planes[iLayer]].GetPadLength() * 0.5 - 1.0 - cseed[planes[iLayer]].GetZref(0);
						if (min > minmax[0]) minmax[0] = min;
					}
					Double_t xpos[4];
					for(Int_t l = 0; l < kNSeedPlanes; l++) xpos[l] = layer[l]->GetX();
					TTreeSRedirector &cstreamer = *fDebugStreamer;
							cstreamer << "MakeSeeds1"
						<< "isFake=" << isFake
						<< "config="   << config
						<< "Cl0.="   << c[0]
						<< "Cl1.="   << c[1]
						<< "Cl2.="   << c[2]
						<< "Cl3.="   << c[3]
						<< "X0="     << xpos[0] //layer[sLayer]->GetX()
						<< "X1="     << xpos[1] //layer[sLayer + 1]->GetX()
						<< "X2="     << xpos[2] //layer[sLayer + 2]->GetX()
						<< "X3="     << xpos[3] //layer[sLayer + 3]->GetX()
						<< "Y2exp="  << cond2[0]
						<< "Z2exp="  << cond2[1]
						<< "Chi2R="  << chi2[0]
						<< "Chi2Z="  << chi2[1]
						<< "Seed0.=" << &cseed[planes[0]]
						<< "Seed1.=" << &cseed[planes[1]]
						<< "Seed2.=" << &cseed[planes[2]]
						<< "Seed3.=" << &cseed[planes[3]]
						<< "Zmin="   << minmax[0]
						<< "Zmax="   << minmax[1]
						<< "\n" ;
				}		
#endif
				// try attaching clusters to tracklets
				Int_t nUsedCl = 0;
				Int_t nlayers = 0;
				for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
					Int_t jLayer = planes[iLayer];
					if(!cseed[jLayer].AttachClustersIter(&layers[jLayer*fTimeBinsPerPlane], 5., kFALSE, c[iLayer])) continue;
					nUsedCl += cseed[jLayer].GetNUsed();
					if(nUsedCl > 25) break;
					nlayers++;
				}
				if(nlayers < kNSeedPlanes){ 
					//AliInfo("Failed updating all seeds.");
					continue;
				}
				// fit tracklets and cook likelihood
				chi2[0] = 0.; chi2[1] = 0.;
				fFitter->FitRieman(&cseed[0], &planes[0]);
				AliRieman *rim = fFitter->GetRiemanFitter();
				for(int iLayer=0; iLayer<4; iLayer++){
					cseed[planes[iLayer]].Init(rim);
					chi2[0] += (Float_t)cseed[planes[iLayer]].GetChi2Z();
					chi2[1] += cseed[planes[iLayer]].GetChi2Y();
				}
				Double_t chi2r = chi2[1], chi2z = chi2[0];
				Double_t like = CookLikelihood(&cseed[0], planes, chi2); // to be checked
				if (TMath::Log(1.E-9 + like) < fRecoParam->GetTrackLikelihood()){
					//AliInfo(Form("Failed likelihood %f[%e].", TMath::Log(1.E-9 + like), like));
					continue;
				}
				//AliInfo(Form("Passed likelihood %f[%e].", TMath::Log(1.E-9 + like), like));


				// book preliminary results
				seedQuality[ntracks] = like;
				fSeedLayer[ntracks]  = config;/*sLayer;*/

				// attach clusters to the extrapolation seeds
				Int_t lextrap[2];
				GetExtrapolationConfig(config, lextrap);
				Int_t nusedf   = 0; // debug value
				for(int iLayer=0; iLayer<2; iLayer++){
					Int_t jLayer = lextrap[iLayer];
					
					// prepare extrapolated seed
					cseed[jLayer].Reset();
					cseed[jLayer].SetRecoParam(fRecoParam);
					cseed[jLayer].SetPlane(jLayer);
					cseed[jLayer].SetTilt(hL[jLayer]);
					cseed[jLayer].SetX0(layers[(jLayer +1) * fTimeBinsPerPlane-1].GetX());
					cseed[jLayer].SetPadLength(padlength[jLayer]);
					cseed[jLayer].SetNTimeBinsRange(tbRange[jLayer]);
					cseed[jLayer].Init(rim);
// 					AliTRDcluster *cd = FindSeedingCluster(&layers[jLayer*fTimeBinsPerPlane], &cseed[jLayer]);
// 					if(cd == 0x0) continue;

					// fit extrapolated seed
					AliTRDseedV1::FitRiemanTilt(cseed, kTRUE);
					if ((jLayer == 0) && !(cseed[1].IsOK())) continue;
					if ((jLayer == 5) && !(cseed[4].IsOK())) continue;
					AliTRDseedV1 tseed = cseed[jLayer];
					if(!tseed.AttachClustersIter(&layers[jLayer*fTimeBinsPerPlane], 1000.)) continue;
					cseed[jLayer] = tseed;
					nusedf += cseed[jLayer].GetNUsed(); // debug value
					AliTRDseedV1::FitRiemanTilt(cseed, kTRUE);
				}
				//AliInfo("Extrapolation done.");

				ImproveSeedQuality(layers, cseed);
				//AliInfo("Improve seed quality done.");

				nlayers   = 0;
				Int_t nclusters = 0;
				Int_t findable  = 0;
				for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
					if (TMath::Abs(cseed[iLayer].GetYref(0) / cseed[iLayer].GetX0()) < 0.15) findable++;
					if (!cseed[iLayer].IsOK()) continue;
					nclusters += cseed[iLayer].GetN2();
					nlayers++;
				}
				if (nlayers < 3){ 
					//AliInfo("Failed quality check on seeds.");
					continue;
				}

				// fit full track and cook likelihoods
				fFitter->FitRieman(&cseed[0]);
				Double_t chi2ZF = 0., chi2RF = 0.;
				for(int ilayer=0; ilayer<6; ilayer++){
					cseed[ilayer].Init(fFitter->GetRiemanFitter());
					if (!cseed[ilayer].IsOK()) continue;
					//tchi2 = cseed[ilayer].GetChi2Z();
					//printf("layer %d chi2 %e\n", ilayer, tchi2);
					chi2ZF += cseed[ilayer].GetChi2Z();
					chi2RF += cseed[ilayer].GetChi2Y();
				}
				chi2ZF /= TMath::Max((nlayers - 3.), 1.);
				chi2RF /= TMath::Max((nlayers - 3.), 1.);

				// do the final track fitting
				fFitter->SetLayers(nlayers);
#ifdef DEBUG
				fFitter->SetDebugStream(fDebugStreamer);
#endif
				fTrackQuality[ntracks] = fFitter->FitHyperplane(&cseed[0], chi2ZF, GetZ());
				Double_t param[3];
				Double_t chi2[2];
				fFitter->GetHyperplaneFitResults(param);
				fFitter->GetHyperplaneFitChi2(chi2);
				//AliInfo("Hyperplane fit done\n");

				// finalize tracklets
				Int_t labels[12];
				Int_t outlab[24];
				Int_t nlab = 0;
				for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
					if (!cseed[iLayer].IsOK()) continue;

					if (cseed[iLayer].GetLabels(0) >= 0) {
						labels[nlab] = cseed[iLayer].GetLabels(0);
						nlab++;
					}

					if (cseed[iLayer].GetLabels(1) >= 0) {
						labels[nlab] = cseed[iLayer].GetLabels(1);
						nlab++;
					}
				}
				Freq(nlab,labels,outlab,kFALSE);
				Int_t label     = outlab[0];
				Int_t frequency = outlab[1];
				for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
					cseed[iLayer].SetFreq(frequency);
					cseed[iLayer].SetC(param[1]/*cR*/);
					cseed[iLayer].SetCC(param[0]/*cC*/);
					cseed[iLayer].SetChi2(chi2[0]);
					cseed[iLayer].SetChi2Z(chi2ZF);
				}
	    
#ifdef DEBUG
				if(AliTRDReconstructor::StreamLevel() >= 2){
					Double_t curv = (fFitter->GetRiemanFitter())->GetC();
					TTreeSRedirector &cstreamer = *fDebugStreamer;
					cstreamer << "MakeSeeds2"
						<< "C="       << curv
						<< "Chi2R="   << chi2r
						<< "Chi2Z="   << chi2z
						<< "Chi2TR="  << chi2[0]
						<< "Chi2TC="  << chi2[1]
						<< "Chi2RF="  << chi2RF
						<< "Chi2ZF="  << chi2ZF
						<< "Ncl="     << nclusters
						<< "Nlayers=" << nlayers
						<< "NUsedS="  << nUsedCl
						<< "NUsed="   << nusedf
						<< "Findable" << findable
						<< "Like="    << like
						<< "S0.="     << &cseed[0]
						<< "S1.="     << &cseed[1]
						<< "S2.="     << &cseed[2]
						<< "S3.="     << &cseed[3]
						<< "S4.="     << &cseed[4]
						<< "S5.="     << &cseed[5]
						<< "Label="   << label
						<< "Freq="    << frequency
						<< "\n";
				}
#endif
				
				ntracks++;
				if(ntracks == kMaxTracksStack){
					AliWarning(Form("Number of seeds reached maximum allowed (%d) in stack.", kMaxTracksStack));
					return ntracks;
				}
				cseed += 6;
			}
		}
	}
	for(int isl=0; isl<4; isl++) delete layer[isl];
	
	return ntracks;
}

//_____________________________________________________________________________
AliTRDtrackV1* AliTRDtrackerV1::MakeTrack(AliTRDseedV1 *seeds, Double_t *params)
{
  //
  // Build a TRD track out of tracklet candidates
  //
  // Parameters :
  //   seeds  : array of tracklets
  //   params : track parameters (see MakeSeeds() function body for a detailed description)
  //
  // Output :
  //   The TRD track.
  //
  // Detailed description
  //
  // To be discussed with Marian !!
  //

  Double_t alpha = AliTRDgeometry::GetAlpha();
  Double_t shift = AliTRDgeometry::GetAlpha()/2.0;
  Double_t c[15];

  c[ 0] = 0.2;
  c[ 1] = 0.0; c[ 2] = 2.0;
  c[ 3] = 0.0; c[ 4] = 0.0; c[ 5] = 0.02;
  c[ 6] = 0.0; c[ 7] = 0.0; c[ 8] = 0.0;  c[ 9] = 0.1;
  c[10] = 0.0; c[11] = 0.0; c[12] = 0.0;  c[13] = 0.0; c[14] = params[5]*params[5]*0.01;

  AliTRDtrackV1 *track = new AliTRDtrackV1(seeds, &params[1], c, params[0], params[6]*alpha+shift);
	track->PropagateTo(params[0]-5.0);
  track->ResetCovariance(1);
  Int_t nc = FollowBackProlongation(*track);
	AliInfo(Form("N clusters for track %d", nc));
	if (nc < 30) {
    delete track;
    track = 0x0;
  } else {
//     track->CookdEdx();
//     track->CookdEdxTimBin(-1);
//     CookLabel(track, 0.9);
  }

  return track;
}

//____________________________________________________________________
void AliTRDtrackerV1::CookLabel(AliKalmanTrack */*pt*/, Float_t /*wrong*/) const
{
	// to be implemented, preferably at the level of TRD tracklet. !!!!!!!
}

//____________________________________________________________________
void AliTRDtrackerV1::ImproveSeedQuality(AliTRDstackLayer *layers
                                       , AliTRDseedV1 *cseed)
{
  //
  // Sort tracklets according to "quality" and try to "improve" the first 4 worst
  //
  // Parameters :
  //  layers : Array of propagation layers for a stack/supermodule
  //  cseed  : Array of 6 seeding tracklets which has to be improved
  // 
  // Output :
  //   cssed : Improved seeds
  // 
  // Detailed description
  //
  // Iterative procedure in which new clusters are searched for each
  // tracklet seed such that the seed quality (see AliTRDseed::GetQuality())
  // can be maximized. If some optimization is found the old seeds are replaced.
  //
	
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	
	// make a local working copy
	AliTRDseedV1 bseed[6];
	for (Int_t jLayer = 0; jLayer < 6; jLayer++) bseed[jLayer] = cseed[jLayer];


	Float_t lastquality = 10000.0;
	Float_t lastchi2    = 10000.0;
	Float_t chi2        =  1000.0;

	for (Int_t iter = 0; iter < 4; iter++) {
		Float_t sumquality = 0.0;
		Float_t squality[6];
		Int_t   sortindexes[6];

		for (Int_t jLayer = 0; jLayer < 6; jLayer++) {
			squality[jLayer]  = bseed[jLayer].IsOK() ? bseed[jLayer].GetQuality(kTRUE) : -1.;
			sumquality +=squality[jLayer];
		}
		if ((sumquality >= lastquality) || (chi2       >     lastchi2)) break;

		
		lastquality = sumquality;
		lastchi2    = chi2;
		if (iter > 0) for (Int_t jLayer = 0; jLayer < 6; jLayer++) cseed[jLayer] = bseed[jLayer];

		
		TMath::Sort(6, squality, sortindexes, kFALSE);
		for (Int_t jLayer = 5; jLayer > 1; jLayer--) {
			Int_t bLayer = sortindexes[jLayer];
			bseed[bLayer].AttachClustersIter(&layers[bLayer*nTimeBins], squality[bLayer], kTRUE);
		}

		chi2 = AliTRDseedV1::FitRiemanTilt(bseed,kTRUE);
	} // Loop: iter
}

//____________________________________________________________________
Double_t  AliTRDtrackerV1::MakeSeedingPlanes(AliTRDstackLayer *layers)
{
  //
  // Calculate plane quality for seeding.
  // 
  //
  // Parameters :
  //   layers : Array of propagation layers for this plane.
  //
  // Output :
  //   plane quality factor for seeding
  // 
  // Detailed description
  //
  // The quality of the plane for seeding is higher if:
  //  1. the average timebin population is closer to an integer number
  //  2. the distribution of clusters/timebin is closer to a uniform distribution.
  //    - the slope of the first derivative of a parabolic fit is small or
  //    - the slope of a linear fit is small
  //

	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();

// 	Double_t x;
// 	TLinearFitter fitter(1, "pol1");
// 	fitter.ClearPoints();
	Int_t ncl = 0;
	Int_t nused = 0;
	Int_t nClLayer;
	for(int itb=0; itb<nTimeBins; itb++){
		//x = layer[itb].GetX();
		//printf("x[%d] = %f nCls %d\n", itb, x, layer[itb].GetNClusters());
		//if(!layer[itb].GetNClusters()) continue;
		//fitter.AddPoint(&x, layer[itb].GetNClusters(), 1.);
		nClLayer = layers[itb].GetNClusters();
		ncl += nClLayer;
		for(Int_t incl = 0; incl < nClLayer; incl++)
			if((layers[itb].GetCluster(incl))->IsUsed()) nused++;
	}
	
	// calculate the deviation of the mean number of clusters from the
	// closest integer values
	Float_t nclMed = float(ncl-nused)/nTimeBins;
	Int_t ncli = Int_t(nclMed);
	Float_t nclDev = TMath::Abs(nclMed - TMath::Max(ncli, 1));
	nclDev -= (nclDev>.5) && ncli ? 0. : 1.; 
	/*Double_t quality = */ return TMath::Exp(2.*nclDev);

// 	// get slope of the derivative
// 	if(!fitter.Eval()) return quality;
// 	fitter.PrintResults(3);
// 	Double_t a = fitter.GetParameter(1);
// 
// 	printf("ncl_dev(%f)  a(%f)\n", ncl_dev, a);
// 	return quality*TMath::Exp(-a);
}

//____________________________________________________________________
Double_t AliTRDtrackerV1::CookLikelihood(AliTRDseedV1 *cseed
                                       , Int_t planes[4]
                                       , Double_t *chi2)
{
  //
  // Calculate the probability of this track candidate.
  //
  // Parameters :
  //   cseeds : array of candidate tracklets
  //   planes : array of seeding planes (see seeding configuration)
  //   chi2   : chi2 values (on the Z and Y direction) from the rieman fit of the track.
  //
  // Output :
  //   likelihood value
  // 
  // Detailed description
  //
  // The track quality is estimated based on the following 4 criteria:
  //  1. precision of the rieman fit on the Y direction (likea)
  //  2. chi2 on the Y direction (likechi2y)
  //  3. chi2 on the Z direction (likechi2z)
  //  4. number of attached clusters compared to a reference value 
  //     (see AliTRDrecoParam::fkFindable) (likeN)
  //
  // The distributions for each type of probabilities are given below as of
  // (date). They have to be checked to assure consistency of estimation.
  //
 
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	// ratio of the total number of clusters/track which are expected to be found by the tracker.
	Float_t fgFindable = fRecoParam->GetFindableClusters();

	
	Int_t nclusters = 0;
	Double_t sumda = 0.;
	for(UChar_t ilayer = 0; ilayer < 4; ilayer++){
		Int_t jlayer = planes[ilayer];
		nclusters += cseed[jlayer].GetN2();
		sumda += TMath::Abs(cseed[jlayer].GetYfitR(1) - cseed[jlayer].GetYref(1));
	}
	Double_t likea     = TMath::Exp(-sumda*10.6);
	Double_t likechi2y  = 0.0000000001;
	if (chi2[1] < 0.5) likechi2y += TMath::Exp(-TMath::Sqrt(chi2[1]) * 7.73);
	Double_t likechi2z = TMath::Exp(-chi2[0] * 0.088) / TMath::Exp(-chi2[0] * 0.019);
	Int_t enc = Int_t(fgFindable*4.*nTimeBins); 	// Expected Number Of Clusters, normally 72
	Double_t likeN     = TMath::Exp(-(enc - nclusters) * 0.19);
	
	Double_t like      = likea * likechi2y * likechi2z * likeN;

#ifdef DEBUG
	//AliInfo(Form("sumda(%f) chi2[0](%f) chi2[1](%f) likea(%f) likechi2y(%f) likechi2z(%f) nclusters(%d) likeN(%f)", sumda, chi2[0], chi2[1], likea, likechi2y, likechi2z, nclusters, likeN));
	if(AliTRDReconstructor::StreamLevel() >= 2){
		TTreeSRedirector &cstreamer = *fDebugStreamer;
		cstreamer << "CookLikelihood"
			<< "sumda="     << sumda
			<< "chi0="      << chi2[0]
			<< "chi1="      << chi2[1]
			<< "likea="     << likea
			<< "likechi2y=" << likechi2y
			<< "likechi2z=" << likechi2z
			<< "nclusters=" << nclusters
			<< "likeN="     << likeN
			<< "like="      << like
			<< "\n";
 	}
#endif

	return like;
}

//___________________________________________________________________
void AliTRDtrackerV1::GetMeanCLStack(AliTRDstackLayer *layers
                                   , Int_t *planes
                                   , Double_t *params)
{
  //
  // Determines the Mean number of clusters per layer.
  // Needed to determine good Seeding Layers
  //
  // Parameters:
  //    - Array of AliTRDstackLayers
  //    - Container for the params
  //
  // Detailed description
  //
  // Two Iterations:
  // In the first Iteration the mean is calculted using all layers.
  // After this, all layers outside the 1-sigma-region are rejected.
  // Then the mean value and the standard-deviation are calculted a second
  // time in order to select all layers in the 1-sigma-region as good-candidates.
  //

	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	
	Float_t mean = 0, stdev = 0;
	Double_t ncl[kNTimeBins*kNSeedPlanes], mcl[kNTimeBins*kNSeedPlanes];
	Int_t position = 0;
	memset(ncl, 0, sizeof(Int_t)*kNTimeBins*kNSeedPlanes);
	memset(mcl, 0, sizeof(Int_t)*kNTimeBins*kNSeedPlanes);
	Int_t nused = 0;
	for(Int_t ipl = 0; ipl < kNSeedPlanes; ipl++){
		for(Int_t ils = 0; ils < nTimeBins; ils++){
			position = planes[ipl]*nTimeBins + ils;
			ncl[ipl * nTimeBins + ils] = layers[position].GetNClusters();
			nused = 0;
			for(Int_t icl = 0; icl < ncl[ipl * nTimeBins + ils]; icl++)
				if((layers[position].GetCluster(icl))->IsUsed()) nused++;
			ncl[ipl * nTimeBins + ils] -= nused;
		}
	}
	// Declaration of quartils:
	//Double_t qvals[3] = {0.0, 0.0, 0.0};
	//Double_t qprop[3] = {0.16667, 0.5, 0.83333};
	// Iterations
	Int_t counter;
	Double_t *array;
	Int_t *limit;
	Int_t nLayers = nTimeBins * kNSeedPlanes;
	for(Int_t iter = 0; iter < 2; iter++){
		array = (iter == 0) ? &ncl[0] : &mcl[0];
		limit = (iter == 0) ? &nLayers : &counter;
		counter = 0;
		if(iter == 1){
			for(Int_t i = 0; i < nTimeBins *kNSeedPlanes; i++){
				if((ncl[i] >  mean + stdev) || (ncl[i] <  mean - stdev)) continue; // Outside 1-sigma region
// 				if((ncl[i] >  qvals[2]) || (ncl[i] <  qvals[0])) continue; // Outside 1-sigma region
				if(ncl[i] == 0) continue;	                                         // 0-Layers also rejected
				mcl[counter] = ncl[i];
				counter++;
			}
		}
		if(*limit == 0) break;
		printf("Limit = %d\n", *limit);
		//using quartils instead of mean and RMS 
// 		TMath::Quantiles(*limit,3,array,qvals,qprop,kFALSE);
 		mean = TMath::Median(*limit, array, 0x0);
 		stdev  = TMath::RMS(*limit, array);
	}
// 	printf("Quantiles: 0.16667 = %3.3f, 0.5 = %3.3f, 0.83333 = %3.3f\n", qvals[0],qvals[1],qvals[2]);
// 	memcpy(params,qvals,sizeof(Double_t)*3);
	params[1] = (Double_t)TMath::Nint(mean);
	params[0] = (Double_t)TMath::Nint(mean - stdev);
	params[2] = (Double_t)TMath::Nint(mean + stdev);

}

//___________________________________________________________________
Int_t AliTRDtrackerV1::GetSeedingLayers(AliTRDstackLayer *layers
                                      , Double_t *params)
{
  //
  // Algorithm to find optimal seeding layer
  // Layers inside one sigma region (given by Quantiles) are sorted
  // according to their difference.
  // All layers outside are sorted according t
  //
  // Parameters:
  //     - Array of AliTRDstackLayers (in the current plane !!!)
  //     - Container for the Indices of the seeding Layer candidates
  //
  // Output:
  //     - Number of Layers inside the 1-sigma-region
  //
  // The optimal seeding layer should contain the mean number of
  // custers in the layers in one chamber.
  //

	//printf("Params: %3.3f, %3.3f, %3.3f\n", params[0], params[1], params[2]);
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	const Int_t kMaxClustersLayer = AliTRDstackLayer::kMaxClustersLayer;
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	Int_t ncl[kNTimeBins], indices[kNTimeBins], bins[kMaxClustersLayer];
	memset(ncl, 0, sizeof(Int_t)*kNTimeBins);
	memset(indices, 0, sizeof(Int_t)*kNTimeBins);
	memset(bins, 0, sizeof(Int_t)*kMaxClustersLayer);
	Int_t nused = 0;
	for(Int_t ils = 0; ils < nTimeBins; ils++){
		ncl[ils] = layers[ils].GetNClusters();
		nused = 0;
		for(Int_t icl = 0; icl < ncl[ils]; icl++)
			if((layers[ils].GetCluster(icl))->IsUsed()) nused++;
		ncl[ils] -= nused;
	}
	
	Float_t mean = params[1];
	for(Int_t ils = 0; ils < nTimeBins; ils++){
		memmove(indices + bins[ncl[ils]+1] + 1, indices + bins[ncl[ils]+1], sizeof(Int_t)*(nTimeBins - ils));
		indices[bins[ncl[ils]+1]] = ils;
		for(Int_t i = ncl[ils]+1; i < kMaxClustersLayer; i++)
			bins[i]++;
	}
	
	//for(Int_t i = 0; i < nTimeBins; i++) printf("Bin %d = %d\n", i, bins[i]);
	Int_t sbin = -1;
	Int_t nElements;
	Int_t position = 0;
	TRandom *r = new TRandom();
	Int_t iter = 0;
	while(1){
		while(sbin < (Int_t)params[0] || sbin > (Int_t)params[2]){
			// Randomly selecting one bin
			sbin = (Int_t)r->Poisson(mean);
		}
		printf("Bin = %d\n",sbin);
		//Randomly selecting one Layer in the bin
		nElements = bins[sbin + 1] - bins[sbin];
		printf("nElements = %d\n", nElements);
		if(iter == 5){
			position = (Int_t)(gRandom->Rndm()*(nTimeBins-1));
			break;
		}
		else if(nElements==0){
			iter++;
			continue;
		}
		position = (Int_t)(gRandom->Rndm()*(nElements-1)) + bins[sbin];
		break;
	}
	delete r;
	return indices[position];
}

//____________________________________________________________________
AliTRDcluster *AliTRDtrackerV1::FindSeedingCluster(AliTRDstackLayer *layers
                                                 , AliTRDseedV1/*AliRieman*/ *reference)
{
  //
  // Finds a seeding Cluster for the extrapolation chamber.
  //
  // The seeding cluster should be as close as possible to the assumed
  // track which is represented by a Rieman fit.
  // Therefore the selecting criterion is the minimum distance between
  // the best fitting cluster and the Reference which is derived from
  // the AliTRDseed. Because all layers are assumed to be equally good
  // a linear search is performed.
  //
  // Imput parameters: - layers: array of AliTRDstackLayers (in one chamber!!!)
  //                   - sfit: the reference
  //
  // Output:           - the best seeding cluster
  //

	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	
	// distances as squared distances
	Int_t index = 0;
	Float_t ypos = 0.0, zpos = 0.0, distance = 0.0, nearestDistance =100000.0; 
	ypos = reference->GetYref(0);
	zpos = reference->GetZref(0);
	AliTRDcluster *currentBest = 0x0, *temp = 0x0;
	for(Int_t ils = 0; ils < nTimeBins; ils++){
		// Reference positions
// 		ypos = reference->GetYat(layers[ils].GetX());
// 		zpos = reference->GetZat(layers[ils].GetX());
		index = layers[ils].SearchNearestCluster(ypos, zpos, fRecoParam->GetRoad2y(), fRecoParam->GetRoad2z());
		if(index == -1) continue;
		temp = layers[ils].GetCluster(index);
		if(!temp) continue;
		distance = (temp->GetY() - ypos) * (temp->GetY() - ypos) + (temp->GetZ() - zpos) * (temp->GetZ() - zpos);
		if(distance < nearestDistance){
			nearestDistance = distance;
			currentBest = temp;
		}
	}
	return currentBest;
}

//____________________________________________________________________
AliTRDstackLayer *AliTRDtrackerV1::MakeSeedingLayer(AliTRDstackLayer *layers
                                                  , Int_t plane)
{
  //
  // Creates a seeding layer
  //
	
	// constants
	const Int_t kMaxRows = 16;
	const Int_t kMaxCols = 144;
	const Int_t kMaxPads = 2304;
	
	// Get the calculation
	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	
	// Get the geometrical data of the chamber
	AliTRDpadPlane *pp = fGeom->GetPadPlane(plane, layers[0].GetStackNr());
	Int_t nCols = pp->GetNcols();
	Float_t ymin = TMath::Min(pp->GetCol0(), pp->GetColEnd());
	Float_t ymax = TMath::Max(pp->GetCol0(), pp->GetColEnd());
	Float_t zmin = TMath::Min(pp->GetRow0(), pp->GetRowEnd());
	Float_t zmax = TMath::Max(pp->GetRow0(), pp->GetRowEnd());
	Int_t nRows = pp->GetNrows();
	Float_t binlength = (ymax - ymin)/nCols; 
	//AliInfo(Form("ymin(%f) ymax(%f) zmin(%f) zmax(%f) nRows(%d) binlength(%f)", ymin, ymax, zmin, zmax, nRows, binlength));
	
	// Fill the histogram
	Int_t arrpos;
	Float_t ypos;
	Int_t irow, nClusters;	
	Int_t *histogram[kMaxRows];											// 2D-Histogram
	Int_t hvals[kMaxPads];	memset(hvals, 0, sizeof(Int_t)*kMaxPads);	
	Float_t *sigmas[kMaxRows];
	Float_t svals[kMaxPads];	memset(svals, 0, sizeof(Float_t)*kMaxPads);	
	AliTRDcluster *c = 0x0;
	for(Int_t irs = 0; irs < kMaxRows; irs++){
		histogram[irs] = &hvals[irs*kMaxCols];
		sigmas[irs] = &svals[irs*kMaxCols];
	}
	for(Int_t iTime = 0; iTime < nTimeBins; iTime++){
		nClusters = layers[iTime].GetNClusters();
		for(Int_t incl = 0; incl < nClusters; incl++){
			c = layers[iTime].GetCluster(incl);	
			ypos = c->GetY();
			if(ypos > ymax && ypos < ymin) continue;
			irow = pp->GetPadRowNumber(c->GetZ());				// Zbin
			if(irow < 0)continue;
			arrpos = static_cast<Int_t>((ypos - ymin)/binlength);
			if(ypos == ymax) arrpos = nCols - 1;
			histogram[irow][arrpos]++;
			sigmas[irow][arrpos] += c->GetSigmaZ2();
		}
	}
	
// Now I have everything in the histogram, do the selection
// 	printf("Starting the analysis\n");
	//Int_t nPads = nCols * nRows;
	// This is what we are interested in: The center of gravity of the best candidates
	Float_t cogyvals[kMaxPads]; memset(cogyvals, 0, sizeof(Float_t)*kMaxPads);
	Float_t cogzvals[kMaxPads]; memset(cogzvals, 0, sizeof(Float_t)*kMaxPads);
	Float_t *cogy[kMaxRows];
	Float_t *cogz[kMaxRows];
	// Lookup-Table storing coordinates according ti the bins
	Float_t yLengths[kMaxCols];
	Float_t zLengths[kMaxRows];
	for(Int_t icnt = 0; icnt < nCols; icnt++){
 		yLengths[icnt] = pp->GetColPos(nCols - 1 - icnt) + binlength/2;
	}
	for(Int_t icnt = 0; icnt < nRows; icnt++){
		zLengths[icnt] = pp->GetRowPos(icnt) - pp->GetRowSize(icnt)/2;
	}

	// A bitfield is used to mask the pads as usable
	Short_t mask[kMaxCols]; memset(mask, 0 ,sizeof(Short_t) * kMaxCols);//bool mvals[kMaxPads];
	for(UChar_t icount = 0; icount < nRows; icount++){
		cogy[icount] = &cogyvals[icount*kMaxCols];
		cogz[icount] = &cogzvals[icount*kMaxCols];
	}
	// In this array the array position of the best candidates will be stored
	Int_t cand[kMaxTracksStack];
	Float_t sigcands[kMaxTracksStack];
 	
	// helper variables
	Int_t indices[kMaxPads]; memset(indices, 0, sizeof(Int_t)*kMaxPads);
	Int_t nCandidates = 0;
 	Float_t norm, cogv;
	// histogram filled -> Select best bins
	TMath::Sort(kMaxPads, hvals, indices);			// bins storing a 0 should not matter
	// Set Threshold
	Int_t maximum = hvals[indices[0]];	// best
	Int_t threshold = static_cast<UChar_t>(maximum * fRecoParam->GetFindableClusters());
	Int_t col, row, lower, lower1, upper, upper1;
	for(Int_t ib = 0; ib < kMaxPads; ib++){
		if(nCandidates >= kMaxTracksStack){
			AliWarning(Form("Number of seed candidates %d exceeded maximum allowed per stack %d", nCandidates, kMaxTracksStack));
			break;
		}
		// Positions
		row = indices[ib]/nCols;
		col = indices[ib]%nCols;
		// here will be the threshold condition:
		if((mask[col] & (1 << row)) != 0) continue;		// Pad is masked: continue
		if(histogram[row][col] < TMath::Max(threshold, 1)){	// of course at least one cluster is needed
			break;			// number of clusters below threshold: break;
		} 
		// passing: Mark the neighbors
		lower  = TMath::Max(col - 1, 0); upper  = TMath::Min(col + 2, nCols);
		lower1 = TMath::Max(row - 1, 0); upper1 = TMath::Min(row + 2, nCols);
		for(Int_t ic = lower; ic < upper; ++ic)
			for(Int_t ir = lower1; ir < upper1; ++ir){
				if(ic == col && ir == row) continue;
				mask[ic] |= (1 << ir);
			}
		// Storing the position in an array
		// testing for neigboring
		cogv = 0;
		norm = 0;
		lower = TMath::Max(col - 1,0);
		upper = TMath::Min(col + 2, nCols);
		for(Int_t inb = lower; inb < upper; ++inb){
			cogv += yLengths[inb] * histogram[row][inb];
			norm += histogram[row][inb];
		}
		cogy[row][col] = cogv / norm;
		cogv = 0; norm = 0;
		lower = TMath::Max(row - 1, 0);
		upper = TMath::Min(row + 2, nRows);
		for(Int_t inb = lower; inb < upper; ++inb){
			cogv += zLengths[inb] * histogram[inb][col];
			norm += histogram[inb][col];
		}
		cogz[row][col] = cogv /  norm;
		// passed the filter
		cand[nCandidates] = row*kMaxCols + col;	// store the position of a passig candidate into an Array
		sigcands[nCandidates] = sigmas[row][col] / histogram[row][col]; // never be a floating point exeption
		// Analysis output
		nCandidates++;
	}
	AliTRDstackLayer *fakeLayer = new AliTRDstackLayer(layers[0].GetZ0(), layers[0].GetDZ0(), layers[0].GetStackNr());
	fakeLayer->SetX((TMath::Abs(layers[nTimeBins-1].GetX() + layers[0].GetX()))/2);
	fakeLayer->SetSector(layers[0].GetSector());
	AliTRDcluster **fakeClusters = 0x0;
	UInt_t *fakeIndices = 0x0;
	if(nCandidates){
		fakeClusters = new AliTRDcluster*[nCandidates];
	 	fakeIndices = new UInt_t[nCandidates];
		UInt_t fakeIndex = 0;
		for(Int_t ican = 0; ican < nCandidates; ican++){
			fakeClusters[ican] = new AliTRDcluster();
			fakeClusters[ican]->SetX(fakeLayer->GetX());
			fakeClusters[ican]->SetY(cogyvals[cand[ican]]);
			fakeClusters[ican]->SetZ(cogzvals[cand[ican]]);
			fakeClusters[ican]->SetSigmaZ2(sigcands[ican]);
			fakeIndices[ican] = fakeIndex++;// fantasy number
		}
	}
	fakeLayer->SetRecoParam(fRecoParam);
	fakeLayer->SetClustersArray(fakeClusters, nCandidates);
	fakeLayer->SetIndexArray(fakeIndices);
	fakeLayer->SetNRows(nRows);
	fakeLayer->BuildIndices();
	//fakeLayer->PrintClusters();
	
#ifdef DEBUG
	if(AliTRDReconstructor::StreamLevel() >= 3){
		TMatrixD hist(nRows, nCols);
		for(Int_t i = 0; i < nRows; i++)
			for(Int_t j = 0; j < nCols; j++)
				hist(i,j) = histogram[i][j];
		TTreeSRedirector &cstreamer = *fDebugStreamer;
		cstreamer << "MakeSeedingLayer"
			<< "Iteration="  << fSieveSeeding
			<< "plane="      << plane
			<< "ymin="       << ymin
			<< "ymax="       << ymax
			<< "zmin="       << zmin
			<< "zmax="       << zmax
			<< "L.="         << fakeLayer
			<< "Histogram.=" << &hist
			<< "\n";
	}
#endif
	return fakeLayer;
}

//____________________________________________________________________
void AliTRDtrackerV1::GetSeedingConfig(Int_t iconfig, Int_t planes[4])
{
  //
  // Map seeding configurations to detector planes.
  //
  // Parameters :
  //   iconfig : configuration index
  //   planes  : member planes of this configuration. On input empty.
  //
  // Output :
  //   planes : contains the planes which are defining the configuration
  // 
  // Detailed description
  //
  // Here is the list of seeding planes configurations together with
  // their topological classification:
  //
  //  0 - 5432 TQ 0
  //  1 - 4321 TQ 0
  //  2 - 3210 TQ 0
  //  3 - 5321 TQ 1
  //  4 - 4210 TQ 1
  //  5 - 5431 TQ 1
  //  6 - 4320 TQ 1
  //  7 - 5430 TQ 2
  //  8 - 5210 TQ 2
  //  9 - 5421 TQ 3
  // 10 - 4310 TQ 3
  // 11 - 5410 TQ 4
  // 12 - 5420 TQ 5
  // 13 - 5320 TQ 5
  // 14 - 5310 TQ 5
  //
  // The topologic quality is modeled as follows:
  // 1. The general model is define by the equation:
  //  p(conf) = exp(-conf/2)
  // 2. According to the topologic classification, configurations from the same
  //    class are assigned the agerage value over the model values.
  // 3. Quality values are normalized.
  // 
  // The topologic quality distribution as function of configuration is given below:
  //Begin_Html
  // <img src="gif/topologicQA.gif">
  //End_Html
  //

	switch(iconfig){
	case 0: // 5432 TQ 0
		planes[0] = 2;
		planes[1] = 3;
		planes[2] = 4;
		planes[3] = 5;
		break;
	case 1: // 4321 TQ 0
		planes[0] = 1;
		planes[1] = 2;
		planes[2] = 3;
		planes[3] = 4;
		break;
	case 2: // 3210 TQ 0
		planes[0] = 0;
		planes[1] = 1;
		planes[2] = 2;
		planes[3] = 3;
		break;
	case 3: // 5321 TQ 1
		planes[0] = 1;
		planes[1] = 2;
		planes[2] = 3;
		planes[3] = 5;
		break;
	case 4: // 4210 TQ 1
		planes[0] = 0;
		planes[1] = 1;
		planes[2] = 2;
		planes[3] = 4;
		break;
	case 5: // 5431 TQ 1
		planes[0] = 1;
		planes[1] = 3;
		planes[2] = 4;
		planes[3] = 5;
		break;
	case 6: // 4320 TQ 1
		planes[0] = 0;
		planes[1] = 2;
		planes[2] = 3;
		planes[3] = 4;
		break;
	case 7: // 5430 TQ 2
		planes[0] = 0;
		planes[1] = 3;
		planes[2] = 4;
		planes[3] = 5;
		break;
	case 8: // 5210 TQ 2
		planes[0] = 0;
		planes[1] = 1;
		planes[2] = 2;
		planes[3] = 5;
		break;
	case 9: // 5421 TQ 3
		planes[0] = 1;
		planes[1] = 2;
		planes[2] = 4;
		planes[3] = 5;
		break;
	case 10: // 4310 TQ 3
		planes[0] = 0;
		planes[1] = 1;
		planes[2] = 3;
		planes[3] = 4;
		break;
	case 11: // 5410 TQ 4
		planes[0] = 0;
		planes[1] = 1;
		planes[2] = 4;
		planes[3] = 5;
		break;
	case 12: // 5420 TQ 5
		planes[0] = 0;
		planes[1] = 2;
		planes[2] = 4;
		planes[3] = 5;
		break;
	case 13: // 5320 TQ 5
		planes[0] = 0;
		planes[1] = 2;
		planes[2] = 3;
		planes[3] = 5;
		break;
	case 14: // 5310 TQ 5
		planes[0] = 0;
		planes[1] = 1;
		planes[2] = 3;
		planes[3] = 5;
		break;
	}
}

//____________________________________________________________________
void AliTRDtrackerV1::GetExtrapolationConfig(Int_t iconfig, Int_t planes[2])
{
  //
  // Returns the extrapolation planes for a seeding configuration.
  //
  // Parameters :
  //   iconfig : configuration index
  //   planes  : planes which are not in this configuration. On input empty.
  //
  // Output :
  //   planes : contains the planes which are not in the configuration
  // 
  // Detailed description
  //

	switch(iconfig){
	case 0: // 5432 TQ 0
		planes[0] = 1;
		planes[1] = 0;
		break;
	case 1: // 4321 TQ 0
		planes[0] = 5;
		planes[1] = 0;
		break;
	case 2: // 3210 TQ 0
		planes[0] = 4;
		planes[1] = 5;
		break;
	case 3: // 5321 TQ 1
		planes[0] = 4;
		planes[1] = 0;
		break;
	case 4: // 4210 TQ 1
		planes[0] = 5;
		planes[1] = 3;
		break;
	case 5: // 5431 TQ 1
		planes[0] = 2;
		planes[1] = 0;
		break;
	case 6: // 4320 TQ 1
		planes[0] = 5;
		planes[1] = 1;
		break;
	case 7: // 5430 TQ 2
		planes[0] = 2;
		planes[1] = 1;
		break;
	case 8: // 5210 TQ 2
		planes[0] = 4;
		planes[1] = 3;
		break;
	case 9: // 5421 TQ 3
		planes[0] = 3;
		planes[1] = 0;
		break;
	case 10: // 4310 TQ 3
		planes[0] = 5;
		planes[1] = 2;
		break;
	case 11: // 5410 TQ 4
		planes[0] = 3;
		planes[1] = 2;
		break;
	case 12: // 5420 TQ 5
		planes[0] = 3;
		planes[1] = 1;
		break;
	case 13: // 5320 TQ 5
		planes[0] = 4;
		planes[1] = 1;
		break;
	case 14: // 5310 TQ 5
		planes[0] = 4;
		planes[1] = 2;
		break;
	}
}
