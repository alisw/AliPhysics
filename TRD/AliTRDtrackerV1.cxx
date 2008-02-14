
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

#include "AliTRDtrackerV1.h"
#include "AliTRDtrackingChamber.h"
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
#include "AliTRDchamberTimeBin.h"
#include "AliTRDrecoParam.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "Cal/AliTRDCalDet.h"


ClassImp(AliTRDtrackerV1)


const  Float_t  AliTRDtrackerV1::fgkMinClustersInTrack =  0.5;  //
const  Float_t  AliTRDtrackerV1::fgkLabelFraction      =  0.8;  //
const  Double_t AliTRDtrackerV1::fgkMaxChi2            = 12.0;  //
const  Double_t AliTRDtrackerV1::fgkMaxSnp             =  0.95; // Maximum local sine of the azimuthal angle
const  Double_t AliTRDtrackerV1::fgkMaxStep            =  2.0;  // Maximal step size in propagation 
Double_t AliTRDtrackerV1::fgTopologicQA[kNConfigs] = {
		0.1112, 0.1112, 0.1112, 0.0786, 0.0786,
		0.0786, 0.0786, 0.0579, 0.0579, 0.0474,
		0.0474, 0.0408, 0.0335, 0.0335, 0.0335
};
TTreeSRedirector *AliTRDtrackerV1::fgDebugStreamer = 0x0;
AliRieman* AliTRDtrackerV1::fgRieman = 0x0;
TLinearFitter* AliTRDtrackerV1::fgTiltedRieman = 0x0;
TLinearFitter* AliTRDtrackerV1::fgTiltedRiemanConstrained = 0x0;

//____________________________________________________________________
AliTRDtrackerV1::AliTRDtrackerV1() 
  :AliTracker()
  ,fGeom(new AliTRDgeometry())
  ,fClusters(0x0)
  ,fTracklets(0x0)
  ,fTracks(0x0)
  ,fTimeBinsPerPlane(0)
  ,fSieveSeeding(0)
{
  //
  // Default constructor.
  // 
  if (!AliTRDcalibDB::Instance()) {
    AliFatal("Could not get calibration object");
  }
  fTimeBinsPerPlane = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();

	for (Int_t isector = 0; isector < AliTRDgeometry::kNsect; isector++) new(&fTrSec[isector]) AliTRDtrackingSector(fGeom, isector, fTimeBinsPerPlane);
  
  if(AliTRDReconstructor::StreamLevel() > 1){
		TDirectory *savedir = gDirectory; 
 		fgDebugStreamer    = new TTreeSRedirector("TRD.TrackerDebug.root");
  	savedir->cd();
	}
}

//____________________________________________________________________
AliTRDtrackerV1::~AliTRDtrackerV1()
{ 
  //
  // Destructor
  //
	
	if(fgDebugStreamer) delete fgDebugStreamer;
	if(fgRieman) delete fgRieman;
	if(fgTiltedRieman) delete fgTiltedRieman;
	if(fgTiltedRiemanConstrained) delete fgTiltedRiemanConstrained;
	if(fTracks) {fTracks->Delete(); delete fTracks;}
	if(fTracklets) {fTracklets->Delete(); delete fTracklets;}
	if(fClusters) {fClusters->Delete(); delete fClusters;}
	if(fGeom) delete fGeom;
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

	if(!AliTRDReconstructor::RecoParam()){
		AliError("Reconstruction configuration not initialized. Call first AliTRDReconstructor::SetRecoParam().");
		return 0;
	}
	
	//AliInfo("Start Track Finder ...");
	Int_t ntracks = 0;
	for(int ism=0; ism<AliTRDgeometry::kNsect; ism++){
//	for(int ism=1; ism<2; ism++){
			//AliInfo(Form("Processing supermodule %i ...", ism));
			ntracks += Clusters2TracksSM(ism, esd);
	}
  AliInfo(Form("Number of found tracks : %d", ntracks));
	return ntracks;
}


//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::GetTrackPoint(Int_t index, AliTrackPoint &p) const
{
	//AliInfo(Form("Asking for tracklet %d", index));
	
	if(index<0) return kFALSE;
	AliTRDseedV1 *tracklet = 0x0; 
	if(!(tracklet = (AliTRDseedV1*)fTracklets->UncheckedAt(index))) return kFALSE;
	
	// get detector for this tracklet
	AliTRDcluster *cl = 0x0;
	Int_t ic = 0; do; while(!(cl = tracklet->GetClusters(ic++)));	 
	Int_t  idet     = cl->GetDetector();
		
	Double_t local[3];
	local[0] = tracklet->GetX0(); 
	local[1] = tracklet->GetYfit(0);
	local[2] = tracklet->GetZfit(0);
	Double_t global[3];
	fGeom->RotateBack(idet, local, global);
	p.SetXYZ(global[0],global[1],global[2]);
	
	
	// setting volume id
	AliGeomManager::ELayerID iLayer = AliGeomManager::kTRD1;
	switch (fGeom->GetPlane(idet)) {
	case 0:
		iLayer = AliGeomManager::kTRD1;
		break;
	case 1:
		iLayer = AliGeomManager::kTRD2;
		break;
	case 2:
		iLayer = AliGeomManager::kTRD3;
		break;
	case 3:
		iLayer = AliGeomManager::kTRD4;
		break;
	case 4:
		iLayer = AliGeomManager::kTRD5;
		break;
	case 5:
		iLayer = AliGeomManager::kTRD6;
		break;
	};
	Int_t    modId = fGeom->GetSector(idet) * fGeom->Ncham() + fGeom->GetChamber(idet);
	UShort_t volid = AliGeomManager::LayerToVolUID(iLayer, modId);
	p.SetVolumeID(volid);
		
	return kTRUE;
}

//____________________________________________________________________
TLinearFitter* AliTRDtrackerV1::GetTiltedRiemanFitter()
{
	if(!fgTiltedRieman) fgTiltedRieman = new TLinearFitter(4, "hyp4");
	return fgTiltedRieman;
}

//____________________________________________________________________
TLinearFitter* AliTRDtrackerV1::GetTiltedRiemanFitterConstraint()
{
	if(!fgTiltedRiemanConstrained) fgTiltedRiemanConstrained = new TLinearFitter(2, "hyp2");
	return fgTiltedRiemanConstrained;
}
	
//____________________________________________________________________	
AliRieman* AliTRDtrackerV1::GetRiemanFitter()
{
	if(!fgRieman) fgRieman = new AliRieman(AliTRDtrackingChamber::kNTimeBins * AliTRDgeometry::kNplan);
	return fgRieman;
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
/*		TTreeSRedirector &cstream = *fgDebugStreamer;
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
		delete track;
	}
	

	AliInfo(Form("Number of seeds: %d", nSeed));
	AliInfo(Form("Number of back propagated TRD tracks: %d", found));
			
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
  
  // Calibration monitor
  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if (!calibra) AliInfo("Could not get Calibra instance\n");
  
  
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
    Bool_t kUPDATE = kFALSE;
		Double_t xTPC = 250.0;
    if(FollowProlongation(track)){
			// computes PID for track
			track.CookPID();
			// update calibration references using this track
			if(calibra->GetHisto2d()) calibra->UpdateHistogramsV1(&track);
	
			// Prolongate to TPC
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
				kUPDATE = kTRUE;
			}
		}	 
		
		// Prolongate to TPC without update
		if(!kUPDATE) {
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
    Int_t   index   = 0;
    AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
		if(!tracklet) continue;
		if(!tracklet->IsOK()) AliWarning("tracklet not OK");
		
		t.SetTracklet(tracklet, iplane, index);
		
		Double_t x  = tracklet->GetX0();
    if (x < (t.GetX()-fgkMaxStep) && !PropagateToX(t, x+fgkMaxStep, fgkMaxStep)) break;
    if (!AdjustSector(&t)) break;
     
    // Start global position
    Double_t xyz0[3];
    t.GetXYZ(xyz0);

		// End global position
    Double_t alpha = t.GetAlpha(), y, z;
    if (!t.GetProlongation(x,y,z)) break;    
    Double_t xyz1[3];
    xyz1[0] =  x * TMath::Cos(alpha) - y * TMath::Sin(alpha);
    xyz1[1] =  x * TMath::Sin(alpha) + y * TMath::Cos(alpha);
    xyz1[2] =  z;
				
    // Get material budget
    Double_t param[7];
    AliTracker::MeanMaterialBudget(xyz0, xyz1, param);
    Double_t xrho= param[0]*param[4];
    Double_t xx0 = param[1]; // Get mean propagation parameters

    // Propagate and update		
		t.PropagateTo(x, xx0, xrho);
	  if (!AdjustSector(&t)) break;
	  
    Double_t maxChi2 = t.GetPredictedChi2(tracklet);
	  if (maxChi2 < 1e+10 && t.Update(tracklet, maxChi2)){ 
	  	nClustersExpected += tracklet->GetN();
  	}
  }

	if(AliTRDReconstructor::StreamLevel() > 1){
		Int_t index;
		for(int iplane=0; iplane<6; iplane++){
			AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
			if(!tracklet) continue;
			t.SetTracklet(tracklet, iplane, index);
		}

		TTreeSRedirector &cstreamer = *fgDebugStreamer;
		cstreamer << "FollowProlongation"
			<< "ncl="      << nClustersExpected
			<< "track.="   << &t
			<< "\n";
	}

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
  Double_t clength = AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick();
  AliTRDtrackingChamber *chamber = 0x0;
  
  // Loop through the TRD planes
  for (Int_t iplane = 0; iplane < AliTRDgeometry::Nplan(); iplane++) {
		// BUILD TRACKLET IF NOT ALREADY BUILT
		Double_t x = 0., y, z, alpha;
    AliTRDseedV1 tracklet(*t.GetTracklet(iplane));
		if(!tracklet.IsOK()){
  		alpha = t.GetAlpha();
  		Int_t sector = Int_t(alpha/AliTRDgeometry::GetAlpha() + (alpha>0. ? 0 : AliTRDgeometry::kNsect));

  		if(!fTrSec[sector].GetNChambers()) continue;
  		
  		if((x = fTrSec[sector].GetX(iplane)) < 1.) continue;
		
			if (!t.GetProlongation(x, y, z)) break;
			Int_t stack = fGeom->GetChamber(z, iplane);
			Int_t nCandidates = stack >= 0 ? 1 : 2;
			z -= stack >= 0 ? 0. : 4.; 
			
			for(int icham=0; icham<nCandidates; icham++, z+=8){
				if((stack = fGeom->GetChamber(z, iplane)) < 0) continue;
			
				if(!(chamber = fTrSec[sector].GetChamber(stack, iplane))) continue;
			
				if(chamber->GetNClusters() < fTimeBinsPerPlane*AliTRDReconstructor::AliTRDReconstructor::RecoParam()->GetFindableClusters()) continue;
			
				x = chamber->GetX();
			
				AliTRDpadPlane *pp = fGeom->GetPadPlane(iplane, stack);
				tracklet.SetTilt(TMath::Tan(-TMath::DegToRad()*pp->GetTiltingAngle()));
				tracklet.SetPadLength(pp->GetLengthIPad());
				tracklet.SetPlane(iplane);
				tracklet.SetX0(x);
				tracklet.Init(&t);
				if(!tracklet.AttachClustersIter(chamber, 1000.)) continue;
				tracklet.Init(&t);
	    	
	    	if(tracklet.GetN() < fTimeBinsPerPlane * AliTRDReconstructor::RecoParam()->GetFindableClusters()) continue;
			
				break;
			}
		}
    if(!tracklet.IsOK()){
			if(x < 1.) continue; //temporary
			if(!PropagateToX(t, x-fgkMaxStep, fgkMaxStep)) break;
			if(!AdjustSector(&t)) break;
			if(TMath::Abs(t.GetSnp()) > fgkMaxSnp) break;
    	continue;
    }
    
		// Propagate closer to the current chamber if neccessary 
    x -= clength;
    if (x > (fgkMaxStep + t.GetX()) && !PropagateToX(t, x-fgkMaxStep, fgkMaxStep)) break;
    if (!AdjustSector(&t)) break;
    if (TMath::Abs(t.GetSnp()) > fgkMaxSnp) break;
		
		// load tracklet to the tracker and the track
		Int_t index = SetTracklet(&tracklet);
		t.SetTracklet(&tracklet, iplane, index);
   
   
		// Calculate the mean material budget along the path inside the chamber
    //Calculate global entry and exit positions of the track in chamber (only track prolongation)
    Double_t xyz0[3]; // entry point 
		t.GetXYZ(xyz0);
		alpha = t.GetAlpha();
		x = tracklet.GetX0();
		if (!t.GetProlongation(x, y, z)) break;
		Double_t xyz1[3]; // exit point
		xyz1[0] =  x * TMath::Cos(alpha) - y * TMath::Sin(alpha); 
    xyz1[1] = +x * TMath::Sin(alpha) + y * TMath::Cos(alpha);
    xyz1[2] =  z;
    Double_t param[7];
		AliTracker::MeanMaterialBudget(xyz0, xyz1, param);	
    // The mean propagation parameters
    Double_t xrho = param[0]*param[4]; // density*length
    Double_t xx0  = param[1]; // radiation length
		
		// Propagate and update track
		t.PropagateTo(x, xx0, xrho);
	  if (!AdjustSector(&t)) break;
		Double_t maxChi2 = t.GetPredictedChi2(&tracklet);
		if (maxChi2<1e+10 && t.Update(&tracklet, maxChi2)){ 
			nClustersExpected += tracklet.GetN();
		}
		// Reset material budget if 2 consecutive gold
		if(iplane>0 && tracklet.GetN() + t.GetTracklet(iplane-1)->GetN() > 20) t.SetBudget(2, 0.);

		// Make backup of the track until is gold
		// TO DO update quality check of the track.
		// consider comparison with fTimeBinsRange
		Float_t ratio0 = tracklet.GetN() / Float_t(fTimeBinsPerPlane);
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

	if(AliTRDReconstructor::StreamLevel() > 1){
		TTreeSRedirector &cstreamer = *fgDebugStreamer;
		cstreamer << "FollowBackProlongation"
			<< "ncl="      << nClustersExpected
			<< "track.="   << &t
			<< "\n";
	}
	
	return nClustersExpected;
}

//_________________________________________________________________________
Float_t AliTRDtrackerV1::FitRieman(AliTRDseedV1 *tracklets, Double_t *chi2, Int_t *planes){
//
// Fits a Riemann-circle to the given points without tilting pad correction.
// The fit is performed using an instance of the class AliRieman (equations 
// and transformations see documentation of this class)
// Afterwards all the tracklets are Updated
//
// Parameters: - Array of tracklets (AliTRDseedV1)
//             - Storage for the chi2 values (beginning with direction z)  
//             - Seeding configuration
// Output:     - The curvature
//
  AliRieman *fitter = AliTRDtrackerV1::GetRiemanFitter();
	fitter->Reset();
  Int_t allplanes[] = {0, 1, 2, 3, 4, 5};
  Int_t *ppl = &allplanes[0];
	Int_t maxLayers = 6;
  if(planes){
    maxLayers = 4;
    ppl = planes;
  }
  for(Int_t il = 0; il < maxLayers; il++){
		if(!tracklets[ppl[il]].IsOK()) continue;
    fitter->AddPoint(tracklets[ppl[il]].GetX0(), tracklets[ppl[il]].GetYfitR(0), tracklets[ppl[il]].GetZProb(),1,10);
  }
  fitter->Update();
  // Set the reference position of the fit and calculate the chi2 values
  memset(chi2, 0, sizeof(Double_t) * 2);
	for(Int_t il = 0; il < maxLayers; il++){
		// Reference positions
		tracklets[ppl[il]].Init(fitter);
		
		// chi2
		if((!tracklets[ppl[il]].IsOK()) && (!planes)) continue;
		chi2[0] += tracklets[ppl[il]].GetChi2Z();
		chi2[1] += tracklets[ppl[il]].GetChi2Y();
	}
	return fitter->GetC();
}

//_________________________________________________________________________
void AliTRDtrackerV1::FitRieman(AliTRDcluster **seedcl, Double_t chi2[2])
{
//
// Performs a Riemann helix fit using the seedclusters as spacepoints
// Afterwards the chi2 values are calculated and the seeds are updated
//
// Parameters: - The four seedclusters
//             - The tracklet array (AliTRDseedV1)
//             - The seeding configuration
//             - Chi2 array
//
// debug level 2
//
	AliRieman *fitter = AliTRDtrackerV1::GetRiemanFitter();
	fitter->Reset();
	for(Int_t i = 0; i < 4; i++)
		fitter->AddPoint(seedcl[i]->GetX(), seedcl[i]->GetY(), seedcl[i]->GetZ(), 1, 10);
	fitter->Update();
	
	
	// Update the seed and calculated the chi2 value
	chi2[0] = 0; chi2[1] = 0;
	for(Int_t ipl = 0; ipl < kNSeedPlanes; ipl++){
		// chi2
		chi2[0] += (seedcl[ipl]->GetZ() - fitter->GetZat(seedcl[ipl]->GetX())) * (seedcl[ipl]->GetZ() - fitter->GetZat(seedcl[ipl]->GetX()));
		chi2[1] += (seedcl[ipl]->GetY() - fitter->GetYat(seedcl[ipl]->GetX())) * (seedcl[ipl]->GetY() - fitter->GetYat(seedcl[ipl]->GetX()));
	}	
}


//_________________________________________________________________________
Float_t AliTRDtrackerV1::FitTiltedRiemanConstraint(AliTRDseedV1 *tracklets, Double_t zVertex)
{
//
// Fits a helix to the clusters. Pad tilting is considered. As constraint it is 
// assumed that the vertex position is set to 0.
// This method is very usefull for high-pt particles
// Basis for the fit: (x - x0)^2 + (y - y0)^2 - R^2 = 0
//      x0, y0: Center of the circle
// Measured y-position: ymeas = y - tan(phiT)(zc - zt)
//      zc: center of the pad row
// Equation which has to be fitted (after transformation):
// a + b * u + e * v + 2*(ymeas + tan(phiT)(z - zVertex))*t = 0
// Transformation:
// t = 1/(x^2 + y^2)
// u = 2 * x * t
// v = 2 * x * tan(phiT) * t
// Parameters in the equation: 
//    a = -1/y0, b = x0/y0, e = dz/dx
//
// The Curvature is calculated by the following equation:
//               - curv = a/Sqrt(b^2 + 1) = 1/R
// Parameters:   - the 6 tracklets
//               - the Vertex constraint
// Output:       - the Chi2 value of the track
//
// debug level 5
//

	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();
	
	TLinearFitter *fitter = GetTiltedRiemanFitterConstraint();
	fitter->StoreData(kTRUE);
	fitter->ClearPoints();

	Float_t x, y, z, w, t, error, tilt;
	Double_t uvt[2];
	Int_t nPoints = 0;
  for(Int_t ipl = 0; ipl < AliTRDgeometry::kNplan; ipl++){
		if(!tracklets[ipl].IsOK()) continue;
		for(Int_t itb = 0; itb < nTimeBins; itb++){
			if(!tracklets[ipl].IsUsable(itb)) continue;
			x = tracklets[ipl].GetX(itb) + tracklets[ipl].GetX0();
			y = tracklets[ipl].GetY(itb);
			z = tracklets[ipl].GetZ(itb);
			tilt = tracklets[ipl].GetTilt();
			// Transformation
			t = 1/(x * x + y * y);
			uvt[0] = 2 * x* t;      
			uvt[1] =  2.0 * tilt * x * t;
			w = 2.0 * (y + tilt * (z - zVertex)) * t;
			error = 2 * 0.2 * t;
			fitter->AddPoint(uvt, w, error);
			nPoints++;
		}
	}
	fitter->Eval();

	// Calculate curvature
	Double_t a = fitter->GetParameter(0);
	Double_t b = fitter->GetParameter(0);
	Double_t curvature = a/TMath::Sqrt(b*b + 1);

	Float_t chi2track = fitter->GetChisquare()/Double_t(nPoints);
	for(Int_t ip = 0; ip < AliTRDtrackerV1::kNPlanes; ip++)
		tracklets[ip].SetCC(curvature);

	if(AliTRDReconstructor::StreamLevel() >= 5){
		//Linear Model on z-direction
	  Double_t xref = (tracklets[2].GetX0() + tracklets[3].GetX0())/2;		// Relative to the middle of the stack
		Double_t slope = fitter->GetParameter(2);
		Double_t zref = slope * xref;
		Float_t chi2Z = CalculateChi2Z(tracklets, zref, slope);
		TTreeSRedirector &treeStreamer = *fgDebugStreamer;
		treeStreamer << "FitTiltedRiemanConstraint"
			<< "Curvature="	<< curvature
			<< "Chi2Track="	<< chi2track
			<< "Chi2Z="			<< chi2Z
			<< "zref="			<< zref
			<< "\n";
	}
	return chi2track;
}

//_________________________________________________________________________
Float_t AliTRDtrackerV1::FitTiltedRieman(AliTRDseedV1 *tracklets, Bool_t sigError)
{
//
// Performs a Riemann fit taking tilting pad correction into account
// The equation of a Riemann circle, where the y position is substituted by the 
// measured y-position taking pad tilting into account, has to be transformed
// into a 4-dimensional hyperplane equation
// Riemann circle: (x-x0)^2 + (y-y0)^2 -R^2 = 0
// Measured y-Position: ymeas = y - tan(phiT)(zc - zt)
//          zc: center of the pad row
//          zt: z-position of the track
// The z-position of the track is assumed to be linear dependent on the x-position
// Transformed equation: a + b * u + c * t + d * v  + e * w - 2 * (ymeas + tan(phiT) * zc) * t = 0
// Transformation:       u = 2 * x * t
//                       v = 2 * tan(phiT) * t
//                       w = 2 * tan(phiT) * (x - xref) * t
//                       t = 1 / (x^2 + ymeas^2)
// Parameters:           a = -1/y0
//                       b = x0/y0
//                       c = (R^2 -x0^2 - y0^2)/y0
//                       d = offset
//                       e = dz/dx
// If the offset respectively the slope in z-position is impossible, the parameters are fixed using 
// results from the simple riemann fit. Afterwards the fit is redone.
// The curvature is calculated according to the formula:
//                       curv = a/(1 + b^2 + c*a) = 1/R
//
// Paramters:   - Array of tracklets (connected to the track candidate)
//              - Flag selecting the error definition
// Output:      - Chi2 value of the track
//
// debug level 5
//

	AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
	Int_t nTimeBins = cal->GetNumberOfTimeBins();

	TLinearFitter *fitter = GetTiltedRiemanFitter();
  fitter->StoreData(kTRUE);
	fitter->ClearPoints();
	
	// Calculate the reference position:
	Int_t nDistances = 0;
	Float_t meanDistance = 0.;
	Int_t startIndex = 5;
	for(Int_t il =5; il > 0; il--){
		if(tracklets[il].IsOK() && tracklets[il -1].IsOK()){
			meanDistance += tracklets[il].GetX0() - tracklets[il -1].GetX0();
			nDistances++;
		}
		if(tracklets[il].IsOK()) startIndex = il;
	}
	meanDistance /= nDistances;
	if(tracklets[0].IsOK()) startIndex = 0;
	Double_t xref = tracklets[startIndex].GetX0() + (2.5 - startIndex) * meanDistance - 0.5 * (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());

	Float_t x, y, z, t, tilt, xdelta, rhs, error;
	Float_t dzMean = 0;	Int_t dzcounter = 0;	// A reference z and a reference slope is used if the fitresults in z-direction are not acceptable
	Double_t uvt[4];
	Int_t nPoints = 0;
	for(Int_t ipl = 0; ipl < AliTRDgeometry::kNplan; ipl++){
		if(!tracklets[ipl].IsOK()) continue;
		dzMean += tracklets[ipl].GetZfitR(1);
		dzcounter++;
		for(Int_t itb = 0; itb < nTimeBins; itb++){
			if (!tracklets[ipl].IsUsable(itb)) continue;
			x = tracklets[ipl].GetX(itb) + tracklets[ipl].GetX0();
			y = tracklets[ipl].GetY(itb);
			z = tracklets[ipl].GetZ(itb);
			tilt = tracklets[ipl].GetTilt();
			xdelta = x - xref;
			// Transformation
			t = 1/(x*x + y*y);
			uvt[0] = 2.0 * x * t;
			uvt[1] = t;
			uvt[2] = 2.0 * tilt * t;
			uvt[3] = 2.0 * tilt * xdelta * t;
			rhs = 2.0 * (y + tilt*z) * t;
			// error definition changes for the different calls
			error = 2.0 * t;
			error *= sigError ? tracklets[ipl].GetSigmaY() : 0.2;
			fitter->AddPoint(uvt, rhs, error);
			nPoints++;
		}
	}
	
	fitter->Eval();

	Double_t offset = fitter->GetParameter(3);
	Double_t slope  = fitter->GetParameter(4);

	// Linear fitter  - not possible to make boundaries
	// Do not accept non possible z and dzdx combinations
	Bool_t acceptablez = kTRUE;
	Double_t zref = 0.0;
	for (Int_t iLayer = 0; iLayer < AliTRDgeometry::kNplan; iLayer++) {
		if(!tracklets[iLayer].IsOK()) continue;
		zref = offset + slope * (tracklets[iLayer].GetX0() - xref);
		if (TMath::Abs(tracklets[iLayer].GetZProb() - zref) > tracklets[iLayer].GetPadLength() * 0.5 + 1.0) 
			acceptablez = kFALSE;
	}
	if (!acceptablez) {
		dzMean /= dzcounter;
		Double_t zmf  = tracklets[startIndex].GetZfitR(0) + dzMean * (xref - tracklets[startIndex].GetX0()); // Z-Position of the track at the middle of a stack assuming a linear dependence on x (approximation)
		fgTiltedRieman->FixParameter(3, zmf);
		fgTiltedRieman->FixParameter(4, dzMean);
		fitter->Eval();
		fitter->ReleaseParameter(3);
		fitter->ReleaseParameter(4);
		offset = fitter->GetParameter(3);
		slope = fitter->GetParameter(4);
	}

	// Calculate Curvarture
	Double_t a     =  fitter->GetParameter(0);
	Double_t b     =  fitter->GetParameter(1);
	Double_t c     =  fitter->GetParameter(2);
	Double_t curvature =  1.0 + b*b - c*a;
	Double_t dca  =  0.0;															// Distance to closest approach
	if (curvature > 0.0) {
		dca = -c / (TMath::Sqrt(1.0 + b*b - c*a) + TMath::Sqrt(1.0 + b*b));
		curvature  =  a / TMath::Sqrt(curvature);
	}

	Double_t chi2track = fitter->GetChisquare()/Double_t(nPoints);

	// Update the tracklets
	Double_t dy, dz;
	for(Int_t iLayer = 0; iLayer < AliTRDtrackerV1::kNPlanes; iLayer++) {

		x  = tracklets[iLayer].GetX0();
		y  = 0;
		z  = 0;
		dy = 0;
		dz = 0;

		// y:     R^2 = (x - x0)^2 + (y - y0)^2
		//     =>   y = y0 +/- Sqrt(R^2 - (x - x0)^2)
		//          R = Sqrt() = 1/Curvature
		//     =>   y = y0 +/- Sqrt(1/Curvature^2 - (x - x0)^2)  
		Double_t res = (x * a + b);								// = (x - x0)/y0
		res *= res;
		res  = 1.0 - c * a + b * b - res;					// = (R^2 - (x - x0)^2)/y0^2
		if (res >= 0) {
			res = TMath::Sqrt(res);
			y    = (1.0 - res) / a;
		}

		// dy:      R^2 = (x - x0)^2 + (y - y0)^2
		//     =>     y = +/- Sqrt(R^2 - (x - x0)^2) + y0
		//     => dy/dx = (x - x0)/Sqrt(R^2 - (x - x0)^2) 
		// Curvature: cr = 1/R = a/Sqrt(1 + b^2 - c*a)
		//     => dy/dx =  (x - x0)/(1/(cr^2) - (x - x0)^2) 
		Double_t x0 = -b / a;
		if (-c * a + b * b + 1 > 0) {
			if (1.0/(curvature * curvature) - (x - x0) * (x - x0) > 0.0) {
				Double_t yderiv = (x - x0) / TMath::Sqrt(1.0/(curvature * curvature) - (x - x0) * (x - x0));
				if (a < 0) yderiv *= -1.0;
				dy = yderiv;
			}
		}
		z  = offset + slope * (x - xref);
		dz = slope;
		tracklets[iLayer].SetYref(0, y);
		tracklets[iLayer].SetYref(1, dy);
		tracklets[iLayer].SetZref(0, z);
		tracklets[iLayer].SetZref(1, dz);
		tracklets[iLayer].SetC(curvature);
		tracklets[iLayer].SetChi2(chi2track);
  }


	if(AliTRDReconstructor::StreamLevel() >= 5){
		Double_t chi2Z = CalculateChi2Z(tracklets, offset, slope);
		TTreeSRedirector &treeStreamer = *fgDebugStreamer;
		treeStreamer << "FitTiltedRieman"
			<< "error="         << sigError
			<< "Curvature="     << curvature
			<< "Chi2track="     << chi2track
			<< "Chi2Z="         << chi2Z
			<< "D="             << c
			<< "DCA="           << dca
			<< "Offset="        << offset
			<< "Slope="         << slope
			<< "\n";
	}

	return chi2track;
}

//_________________________________________________________________________
Float_t AliTRDtrackerV1::CalculateChi2Z(AliTRDseedV1 *tracklets, Double_t offset, Double_t slope)
{
//
// Calculates the chi2-value of the track in z-Direction including tilting pad correction.
// A linear dependence on the x-value serves as a model.
// The parameters are related to the tilted Riemann fit.
// Parameters: - Array of tracklets (AliTRDseedV1) related to the track candidate
//             - the offset for the reference x
//             - the slope
// Output:     - The Chi2 value of the track in z-Direction
//
	Double_t xref = .5 * (tracklets[2].GetX0() + tracklets[3].GetX0());
	Float_t chi2Z = 0, nLayers = 0;
	for (Int_t iLayer = 0; iLayer < AliTRDgeometry::kNplan; iLayer++) {
		if(!tracklets[iLayer].IsOK()) continue;
		Double_t z = offset + slope * (tracklets[iLayer].GetX0() - xref);
		chi2Z += TMath::Abs(tracklets[iLayer].GetMeanz() - z);
		nLayers++;
	}
	chi2Z /= TMath::Max((nLayers - 3.0),1.0);
	return chi2Z;
}



//_____________________________________________________________________________
Int_t AliTRDtrackerV1::PropagateToX(AliTRDtrackV1 &t, Double_t xToGo, Double_t maxStep)
{
  //
  // Starting from current X-position of track <t> this function
  // extrapolates the track up to radial position <xToGo>. 
  // Returns 1 if track reaches the plane, and 0 otherwise 
  //

  const Double_t kEpsilon = 0.00001;

  // Current track X-position
  Double_t xpos = t.GetX();

  // Direction: inward or outward
  Double_t dir  = (xpos < xToGo) ? 1.0 : -1.0;

  while (((xToGo - xpos) * dir) > kEpsilon) {

    Double_t xyz0[3];
    Double_t xyz1[3];
    Double_t param[7];
    Double_t x;
    Double_t y;
    Double_t z;

    // The next step size
    Double_t step = dir * TMath::Min(TMath::Abs(xToGo-xpos),maxStep);

    // Get the global position of the starting point
    t.GetXYZ(xyz0);

    // X-position after next step
    x = xpos + step;

    // Get local Y and Z at the X-position of the next step
    if (!t.GetProlongation(x,y,z)) {
      return 0; // No prolongation possible
    }

    // The global position of the end point of this prolongation step
    xyz1[0] =  x * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha()); 
    xyz1[1] = +x * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;

    // Calculate the mean material budget between start and
    // end point of this prolongation step
    AliTracker::MeanMaterialBudget(xyz0, xyz1, param);

    // Propagate the track to the X-position after the next step
    if (!t.PropagateTo(x,param[1],param[0]*param[4])) {
      return 0;
    }

    // Rotate the track if necessary
    AdjustSector(&t);

    // New track X-position
    xpos = t.GetX();

  }

  return 1;

}


//_____________________________________________________________________________
Int_t AliTRDtrackerV1::ReadClusters(TClonesArray* &array, TTree *clusterTree) const
{
  //
  // Reads AliTRDclusters from the file. 
  // The names of the cluster tree and branches 
  // should match the ones used in AliTRDclusterizer::WriteClusters()
  //

  Int_t nsize = Int_t(clusterTree->GetTotBytes() / (sizeof(AliTRDcluster))); 
  TObjArray *clusterArray = new TObjArray(nsize+1000); 
  
  TBranch *branch = clusterTree->GetBranch("TRDcluster");
  if (!branch) {
    AliError("Can't get the branch !");
    return 1;
  }
  branch->SetAddress(&clusterArray); 
  
  if(!fClusters){ 
  	array = new TClonesArray("AliTRDcluster", nsize);
  	array->SetOwner(kTRUE);
  }
  
  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t) clusterTree->GetEntries();
  Int_t nbytes     = 0;
  Int_t ncl        = 0;
  AliTRDcluster *c = 0x0;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);  
    
    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      if(!(c = (AliTRDcluster *) clusterArray->UncheckedAt(iCluster))) continue;
      new((*fClusters)[ncl++]) AliTRDcluster(*c);
      clusterArray->RemoveAt(iCluster); 
    }

  }
  delete clusterArray;

  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDtrackerV1::LoadClusters(TTree *cTree)
{
  //
  // Fills clusters into TRD tracking_sectors 
  // Note that the numbering scheme for the TRD tracking_sectors 
  // differs from that of TRD sectors
  //

	
  if (ReadClusters(fClusters, cTree)) {
    AliError("Problem with reading the clusters !");
    return 1;
  }
  Int_t ncl  = fClusters->GetEntriesFast(), nin = 0;
  Int_t icl = ncl;
  while (icl--) {
    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(icl);
		if(c->IsInChamber()) nin++;
    Int_t detector       = c->GetDetector();
    Int_t sector         = fGeom->GetSector(detector);
    Int_t stack          = fGeom->GetChamber(detector);
    Int_t plane          = fGeom->GetPlane(detector);
		
		fTrSec[sector].GetChamber(stack, plane, kTRUE)->InsertCluster(c, icl);
  }
  AliInfo(Form("Clusters %d in %6.2f %%", ncl, 100.*float(nin)/ncl));
	
	for(int isector =0; isector<AliTRDgeometry::kNsect; isector++){ 
		if(!fTrSec[isector].GetNChambers()) continue;
		fTrSec[isector].Init();
  }
  
  return 0;
}


//____________________________________________________________________
void AliTRDtrackerV1::UnloadClusters() 
{ 
  //
  // Clears the arrays of clusters and tracks. Resets sectors and timebins 
  //

	if(fTracks) fTracks->Delete(); 
  if(fTracklets) fTracklets->Delete();
  if(fClusters) fClusters->Delete();

  for (int i = 0; i < AliTRDgeometry::kNsect; i++) fTrSec[i].Clear();

}

//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::AdjustSector(AliTRDtrackV1 *track) 
{
  //
  // Rotates the track when necessary
  //

  Double_t alpha = AliTRDgeometry::GetAlpha(); 
  Double_t y     = track->GetY();
  Double_t ymax  = track->GetX()*TMath::Tan(0.5*alpha);

  if      (y >  ymax) {
    if (!track->Rotate( alpha)) {
      return kFALSE;
    }
  } 
  else if (y < -ymax) {
    if (!track->Rotate(-alpha)) {
      return kFALSE;   
    }
  } 

  return kTRUE;

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
	AliTRDseedV1 *tracklet = idx<0 ? 0x0 : (AliTRDseedV1*)fTracklets->UncheckedAt(idx);

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
	new ((*fTracklets)[nentries]) AliTRDseedV1(*tracklet);
	return nentries;
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::Clusters2TracksSM(Int_t sector, AliESDEvent *esd)
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
	
	// allocate space for esd tracks in this SM
	TClonesArray esdTrackList("AliESDtrack", 2*kMaxTracksStack);
	esdTrackList.SetOwner();
	
	Int_t nTracks   = 0;
	Int_t nChambers = 0;
	AliTRDtrackingChamber **stack = 0x0, *chamber = 0x0;
	for(int istack = 0; istack<AliTRDgeometry::kNcham; istack++){
		if(!(stack = fTrSec[sector].GetStack(istack))) continue;
		nChambers = 0;
		for(int iplane=0; iplane<AliTRDgeometry::kNplan; iplane++){
			if(!(chamber = stack[iplane])) continue;
			if(chamber->GetNClusters() < fTimeBinsPerPlane * AliTRDReconstructor::RecoParam()->GetFindableClusters()) continue;
			nChambers++;
			//AliInfo(Form("sector %d stack %d plane %d clusters %d", sector, istack, iplane, chamber->GetNClusters()));
		}
		if(nChambers < 4) continue;
		//AliInfo(Form("Doing stack %d", istack));
		nTracks += Clusters2TracksStack(stack, &esdTrackList);
	}
	//AliInfo(Form("Found %d tracks in SM %d [%d]\n", nTracks, sector, esd->GetNumberOfTracks()));
	
	for(int itrack=0; itrack<nTracks; itrack++) 
          esd->AddTrack((AliESDtrack*)esdTrackList[itrack]);

	return nTracks;
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::Clusters2TracksStack(AliTRDtrackingChamber **stack, TClonesArray *esdTrackList)
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

	AliTRDtrackingChamber *chamber = 0x0;
	AliTRDseedV1 sseed[kMaxTracksStack*6]; // to be initialized
	Int_t pars[4]; // MakeSeeds parameters

	//Double_t alpha = AliTRDgeometry::GetAlpha();
	//Double_t shift = .5 * alpha;
	Int_t configs[kNConfigs];
	
	// Build initial seeding configurations
	Double_t quality = BuildSeedingConfigs(stack, configs);
	if(AliTRDReconstructor::StreamLevel() > 1){
  	AliInfo(Form("Plane config %d %d %d Quality %f"
    , configs[0], configs[1], configs[2], quality));
	}
	
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
			pars[1] = ntracks;
			ntracks = MakeSeeds(stack, &sseed[6*ntracks], pars);
			if(ntracks == kMaxTracksStack) break;
		}
		if(AliTRDReconstructor::StreamLevel() > 1) AliInfo(Form("Candidate TRD tracks %d in iteration %d.", ntracks, fSieveSeeding));
		
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
					if(!sseed[jseed].IsOK()) continue;
					if (TMath::Abs(sseed[jseed].GetYref(0) / sseed[jseed].GetX0()) < 0.15) findable++;
	
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
					//printf("Skip %d nused %d\n", trackIndex, nused);
					fakeTrack[trackIndex] = kTRUE;
					continue;
				}
				if (Float_t(nused)/ncl >= .25){
					//printf("Skip %d nused/ncl >= .25\n", trackIndex);
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
					//printf("REJECTED : %d [%d] nlayers %d trackQuality = %e nused %d\n", itrack, trackIndex, nlayers, fTrackQuality[trackIndex], nused);
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
				Int_t ich = 0; while(!(chamber = stack[ich])) ich++;
				trackParams[6] = fGeom->GetSector(chamber->GetDetector());/* *alpha+shift;	// Supermodule*/

				if(AliTRDReconstructor::StreamLevel() > 1){
					AliInfo(Form("Track %d [%d] nlayers %d trackQuality = %e nused %d, yref = %3.3f", itrack, trackIndex, nlayers, fTrackQuality[trackIndex], nused, trackParams[1]));
					
					Int_t nclusters = 0;
					AliTRDseedV1 *dseed[6];
					for(int is=0; is<6; is++){
						dseed[is] = new AliTRDseedV1(sseed[trackIndex*6+is]);
						dseed[is]->SetOwner();
						nclusters += sseed[is].GetN2();
					}
					//Int_t eventNrInFile = esd->GetEventNumberInFile();
					//AliInfo(Form("Number of clusters %d.", nclusters));
					TTreeSRedirector &cstreamer = *fgDebugStreamer;
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
				}
			
				AliTRDtrackV1 *track = MakeTrack(&sseed[trackIndex*kNPlanes], trackParams);
				if(!track){
					//AliWarning("Fail to build a TRD Track.");
					continue;
				}
				//AliInfo("End of MakeTrack()");
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
		quality = BuildSeedingConfigs(stack, configs);
		if(quality < 1.E-7) break; //AliTRDReconstructor::RecoParam()->GetPlaneQualityThreshold()) break;
		
		for(Int_t ip = 0; ip < kNPlanes; ip++){ 
			if(!(chamber = stack[ip])) continue;
			chamber->Build(fGeom);//Indices(fSieveSeeding);
		}

		if(AliTRDReconstructor::StreamLevel() > 1){ 
			AliInfo(Form("Sieve level %d Plane config %d %d %d Quality %f", fSieveSeeding, configs[0], configs[1], configs[2], quality));
		}
	} while(fSieveSeeding<10); // end stack clusters sieve
	


	//AliInfo(Form("Registered TRD tracks %d in stack %d.", ntracks2, pars[1]));

	return ntracks2;
}

//___________________________________________________________________
Double_t AliTRDtrackerV1::BuildSeedingConfigs(AliTRDtrackingChamber **stack, Int_t *configs)
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

	Double_t chamberQ[kNPlanes];
	AliTRDtrackingChamber *chamber = 0x0;
	for(int iplane=0; iplane<kNPlanes; iplane++){
		if(!(chamber = stack[iplane])) continue;
		chamberQ[iplane] = (chamber = stack[iplane]) ?  chamber->GetQuality(fTimeBinsPerPlane) : 0.;
	}

	Double_t tconfig[kNConfigs];
	Int_t planes[4];
	for(int iconf=0; iconf<kNConfigs; iconf++){
		GetSeedingConfig(iconf, planes);
		tconfig[iconf] = fgTopologicQA[iconf];
		for(int iplane=0; iplane<4; iplane++) tconfig[iconf] *= chamberQ[planes[iplane]]; 
	}
	
	TMath::Sort(kNConfigs, tconfig, configs, kTRUE);
// 	AliInfo(Form("q[%d] = %f", configs[0], tconfig[configs[0]]));
// 	AliInfo(Form("q[%d] = %f", configs[1], tconfig[configs[1]]));
// 	AliInfo(Form("q[%d] = %f", configs[2], tconfig[configs[2]]));
	
	return tconfig[configs[0]];
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::MakeSeeds(AliTRDtrackingChamber **stack, AliTRDseedV1 *sseed, Int_t *ipar)
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
  //   layer 0, 1, and 2 are defined in AliTRDchamberTimeBin::BuildCond().
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

	AliTRDtrackingChamber *chamber = 0x0;
	AliTRDcluster *c[4] = {0x0, 0x0, 0x0, 0x0}; // initilize seeding clusters
	AliTRDseedV1 *cseed = &sseed[0]; // initialize tracklets for first track
	Int_t ncl, mcl; // working variable for looping over clusters
	Int_t index[AliTRDchamberTimeBin::kMaxClustersLayer], jndex[AliTRDchamberTimeBin::kMaxClustersLayer];
	// chi2 storage
	// chi2[0] = tracklet chi2 on the Z direction
	// chi2[1] = tracklet chi2 on the R direction
	Double_t chi2[4];


	// this should be data member of AliTRDtrack
	Double_t seedQuality[kMaxTracksStack];
	
	// unpack control parameters
	Int_t config  = ipar[0];
	Int_t ntracks = ipar[1];
	Int_t planes[kNSeedPlanes]; GetSeedingConfig(config, planes);	
	
	// Init chambers geometry
	Int_t ic = 0; while(!(chamber = stack[ic])) ic++;
	Int_t istack = fGeom->GetChamber(chamber->GetDetector());
	Double_t hL[kNPlanes];       // Tilting angle
	Float_t padlength[kNPlanes]; // pad lenghts
	AliTRDpadPlane *pp = 0x0;
	for(int iplane=0; iplane<kNPlanes; iplane++){
		pp                = fGeom->GetPadPlane(iplane, istack);
		hL[iplane]        = TMath::Tan(-TMath::DegToRad()*pp->GetTiltingAngle());
		padlength[iplane] = pp->GetLengthIPad();
	}
	
	if(AliTRDReconstructor::StreamLevel() > 1){
		AliInfo(Form("Making seeds Stack[%d] Config[%d] Tracks[%d]...", istack, config, ntracks));
	}

	Int_t nlayers = 0;
	AliTRDchamberTimeBin *layer[] = {0x0, 0x0, 0x0, 0x0};
	for(int isl=0; isl<kNSeedPlanes; isl++){ 
		if(!(chamber = stack[planes[isl]])) continue;
		if(!(layer[isl] = chamber->GetSeedingLayer(fGeom))) continue;
		nlayers++;
		//AliInfo(Form("seeding plane %d clusters %d", planes[isl], Int_t(*layer[isl])));
	}
	if(nlayers < 4) return 0;
	
	
	// Start finding seeds
	Double_t cond0[4], cond1[4], cond2[4];
	Int_t icl = 0;
	while((c[3] = (*layer[3])[icl++])){
		if(!c[3]) continue;
		layer[0]->BuildCond(c[3], cond0, 0);
		layer[0]->GetClusters(cond0, index, ncl);
		//printf("Found c[3] candidates 0 %d\n", ncl);
		Int_t jcl = 0;
		while(jcl<ncl) {
			c[0] = (*layer[0])[index[jcl++]];
			if(!c[0]) continue;
			Double_t dx    = c[3]->GetX() - c[0]->GetX();
			Double_t theta = (c[3]->GetZ() - c[0]->GetZ())/dx;
			Double_t phi   = (c[3]->GetY() - c[0]->GetY())/dx;
			layer[1]->BuildCond(c[0], cond1, 1, theta, phi);
			layer[1]->GetClusters(cond1, jndex, mcl);
			//printf("Found c[0] candidates 1 %d\n", mcl);

			Int_t kcl = 0;
			while(kcl<mcl) {
				c[1] = (*layer[1])[jndex[kcl++]];
				if(!c[1]) continue;
				layer[2]->BuildCond(c[1], cond2, 2, theta, phi);
				c[2] = layer[2]->GetNearestCluster(cond2);
				//printf("Found c[1] candidate 2 %p\n", c[2]);
				if(!c[2]) continue;
				
// 				AliInfo("Seeding clusters found. Building seeds ...");
// 				for(Int_t i = 0; i < kNSeedPlanes; i++) printf("%i. coordinates: x = %6.3f, y = %6.3f, z = %6.3f\n", i, c[i]->GetX(), c[i]->GetY(), c[i]->GetZ());
				
				for (Int_t il = 0; il < 6; il++) cseed[il].Reset();

				FitRieman(c, chi2);

				AliTRDseedV1 *tseed = 0x0;
				for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
					Int_t jLayer = planes[iLayer];
					tseed = &cseed[jLayer];
					tseed->SetPlane(jLayer);
					tseed->SetTilt(hL[jLayer]);
					tseed->SetPadLength(padlength[jLayer]);
					tseed->SetX0(stack[jLayer]->GetX());
					tseed->Init(GetRiemanFitter());
				}

				Bool_t isFake = kFALSE;
				if(AliTRDReconstructor::StreamLevel() >= 2){
					if (c[0]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
					if (c[1]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
					if (c[2]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
					
					Float_t yref[4];
					for(int il=0; il<4; il++) yref[il] = cseed[planes[il]].GetYref(0);
					Int_t ll = c[3]->GetLabel(0);
					TTreeSRedirector &cs0 = *fgDebugStreamer;
							cs0 << "MakeSeeds0"
							<<"isFake=" << isFake
							<<"label=" << ll
							<<"chi2z=" << chi2[0]
							<<"chi2y=" << chi2[1]
							<<"yref0=" << yref[0]
							<<"yref1=" << yref[1]
							<<"yref2=" << yref[2]
							<<"yref3=" << yref[3]
							<<"c0.="   << c[0]
							<<"c1.="   << c[1]
							<<"c2.="   << c[2]
							<<"c3.="   << c[3]
							<<"\n";
				}

				if(chi2[0] > AliTRDReconstructor::RecoParam()->GetChi2Z()/*7./(3. - sLayer)*//*iter*/){
					//AliInfo(Form("Failed chi2 filter on chi2Z [%f].", chi2[0]));
					continue;
				}
				if(chi2[1] > AliTRDReconstructor::RecoParam()->GetChi2Y()/*1./(3. - sLayer)*//*iter*/){
					//AliInfo(Form("Failed chi2 filter on chi2Y [%f].", chi2[1]));
					continue;
				}
				//AliInfo("Passed chi2 filter.");

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
					TTreeSRedirector &cstreamer = *fgDebugStreamer;
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
				
				// try attaching clusters to tracklets
				Int_t nUsedCl = 0;
				Int_t nlayers = 0;
				for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
					Int_t jLayer = planes[iLayer];
					if(!cseed[jLayer].AttachClustersIter(stack[jLayer], 5., kFALSE, c[iLayer])) continue;
					nUsedCl += cseed[jLayer].GetNUsed();
					if(nUsedCl > 25) break;
					nlayers++;
				}
				if(nlayers < kNSeedPlanes){ 
					//AliInfo(Form("Failed updating all seeds %d [%d].", nlayers, kNSeedPlanes));
					continue;
				}
				// fit tracklets and cook likelihood
				FitRieman(&cseed[0], chi2, &planes[0]);
				Double_t like = CookLikelihood(&cseed[0], planes, chi2); // to be checked
				if (TMath::Log(1.E-9 + like) < AliTRDReconstructor::RecoParam()->GetTrackLikelihood()){
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
					if(!(chamber = stack[jLayer])) continue;
						
					// prepare extrapolated seed
					cseed[jLayer].Reset();
					cseed[jLayer].SetPlane(jLayer);
					cseed[jLayer].SetTilt(hL[jLayer]);
					cseed[jLayer].SetX0(chamber->GetX());
					cseed[jLayer].SetPadLength(padlength[jLayer]);

					// fit extrapolated seed
					FitTiltedRieman(cseed, kTRUE);
					if ((jLayer == 0) && !(cseed[1].IsOK())) continue;
					if ((jLayer == 5) && !(cseed[4].IsOK())) continue;
					AliTRDseedV1 tseed = cseed[jLayer];
					if(!tseed.AttachClustersIter(chamber, 1000.)) continue;
					cseed[jLayer] = tseed;
					nusedf += cseed[jLayer].GetNUsed(); // debug value
				}
				FitTiltedRieman(cseed, kTRUE);
				//AliInfo("Extrapolation done.");

				if(ImproveSeedQuality(stack, cseed) < 4) continue;
				//AliInfo("Improve seed quality done.");

				// fit full track and cook likelihoods
				Double_t curv = FitRieman(&cseed[0], chi2);
				Double_t chi2ZF = chi2[0] / TMath::Max((nlayers - 3.), 1.);
				Double_t chi2RF = chi2[1] / TMath::Max((nlayers - 3.), 1.);

				// do the final track fitting (Once with vertex constraint and once without vertex constraint)
				Double_t chi2Vals[3];
				chi2Vals[0] = FitTiltedRieman(&cseed[0], kFALSE);
				chi2Vals[1] = FitTiltedRiemanConstraint(&cseed[0], GetZ());
				chi2Vals[2] = chi2ZF;
				fTrackQuality[ntracks] = CalculateTrackLikelihood(&cseed[0], &chi2Vals[0]);
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
					cseed[iLayer].SetChi2Z(chi2ZF);
				}
	    
				if(AliTRDReconstructor::StreamLevel() >= 2){
					TTreeSRedirector &cstreamer = *fgDebugStreamer;
					cstreamer << "MakeSeeds2"
						<< "C="       << curv
						<< "Chi2TR="  << chi2[0]
						<< "Chi2TC="  << chi2[1]
						<< "Chi2RF="  << chi2RF
						<< "Chi2ZF="  << chi2ZF
						<< "Nlayers=" << nlayers
						<< "NUsedS="  << nUsedCl
						<< "NUsed="   << nusedf
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
				
				ntracks++;
				if(ntracks == kMaxTracksStack){
					AliWarning(Form("Number of seeds reached maximum allowed (%d) in stack.", kMaxTracksStack));
					for(int isl=0; isl<4; isl++) delete layer[isl];
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
	//AliInfo(Form("N clusters for track %d", nc));
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
Int_t AliTRDtrackerV1::ImproveSeedQuality(AliTRDtrackingChamber **stack, AliTRDseedV1 *cseed)
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
	
	// make a local working copy
	AliTRDtrackingChamber *chamber = 0x0;
	AliTRDseedV1 bseed[6];
	Int_t nLayers = 0;
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
			sumquality += squality[jLayer];
		}
		if ((sumquality >= lastquality) || (chi2       >     lastchi2)) break;

		nLayers = 0;
		lastquality = sumquality;
		lastchi2    = chi2;
		if (iter > 0) for (Int_t jLayer = 0; jLayer < 6; jLayer++) cseed[jLayer] = bseed[jLayer];

		TMath::Sort(6, squality, sortindexes, kFALSE);
		for (Int_t jLayer = 5; jLayer > 1; jLayer--) {
			Int_t bLayer = sortindexes[jLayer];
			if(!(chamber = stack[bLayer])) continue;
			bseed[bLayer].AttachClustersIter(chamber, squality[bLayer], kTRUE);
			if(bseed[bLayer].IsOK()) nLayers++;
		}

		chi2 = FitTiltedRieman(bseed, kTRUE);
	} // Loop: iter
	
	// we are sure that at least 2 tracklets are OK !
	return nLayers+2;
}

//_________________________________________________________________________
Double_t AliTRDtrackerV1::CalculateTrackLikelihood(AliTRDseedV1 *tracklets, Double_t *chi2){
//
// Calculates the Track Likelihood value. This parameter serves as main quality criterion for 
// the track selection
// The likelihood value containes:
//    - The chi2 values from the both fitters and the chi2 values in z-direction from a linear fit
//    - The Sum of the Parameter  |slope_ref - slope_fit|/Sigma of the tracklets
// For all Parameters an exponential dependency is used
//
// Parameters: - Array of tracklets (AliTRDseedV1) related to the track candidate
//             - Array of chi2 values: 
//                 * Non-Constrained Tilted Riemann fit
//                 * Vertex-Constrained Tilted Riemann fit
//                 * z-Direction from Linear fit
// Output:     - The calculated track likelihood
//
// debug level 2
//

	Double_t sumdaf = 0, nLayers = 0;
	for (Int_t iLayer = 0; iLayer < kNPlanes; iLayer++) {
		if(!tracklets[iLayer].IsOK()) continue;
		sumdaf += TMath::Abs((tracklets[iLayer].GetYfit(1) - tracklets[iLayer].GetYref(1))/ tracklets[iLayer].GetSigmaY2());
		nLayers++;
	}
	sumdaf /= Float_t (nLayers - 2.0);
	
	Double_t likeChi2Z  = TMath::Exp(-chi2[2] * 0.14);			// Chi2Z 
	Double_t likeChi2TC = TMath::Exp(-chi2[1] * 0.677);			// Constrained Tilted Riemann
	Double_t likeChi2TR = TMath::Exp(-chi2[0] * 0.78);			// Non-constrained Tilted Riemann
	Double_t likeAF     = TMath::Exp(-sumdaf * 3.23);
	Double_t trackLikelihood     = likeChi2Z * likeChi2TR * likeAF;

	if(AliTRDReconstructor::StreamLevel() >= 2){
		TTreeSRedirector &cstreamer = *fgDebugStreamer;
		cstreamer << "CalculateTrackLikelihood0"
			<< "LikeChi2Z=" 	<< likeChi2Z
			<< "LikeChi2TR="	<< likeChi2TR
			<< "LikeChi2TC="	<< likeChi2TC
			<< "LikeAF="			<< likeAF
			<< "TrackLikelihood=" << trackLikelihood
			<< "\n";
	}

	return trackLikelihood;
}

//____________________________________________________________________
Double_t AliTRDtrackerV1::CookLikelihood(AliTRDseedV1 *cseed, Int_t planes[4]
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
	Float_t fgFindable = AliTRDReconstructor::RecoParam()->GetFindableClusters();

	
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

	//AliInfo(Form("sumda(%f) chi2[0](%f) chi2[1](%f) likea(%f) likechi2y(%f) likechi2z(%f) nclusters(%d) likeN(%f)", sumda, chi2[0], chi2[1], likea, likechi2y, likechi2z, nclusters, likeN));
	if(AliTRDReconstructor::StreamLevel() >= 2){
		TTreeSRedirector &cstreamer = *fgDebugStreamer;
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

	return like;
}


//___________________________________________________________________
void AliTRDtrackerV1::GetMeanCLStack(AliTRDtrackingChamber *chamber, Int_t *planes, Double_t *params)
{
  //
  // Determines the Mean number of clusters per layer.
  // Needed to determine good Seeding Layers
  //
  // Parameters:
  //    - Array of AliTRDchamberTimeBins
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

	Float_t mean = 0, stdev = 0;
	Double_t ncl[kNTimeBins*kNSeedPlanes], mcl[kNTimeBins*kNSeedPlanes];
	Int_t position = 0;
	memset(ncl, 0, sizeof(Int_t)*kNTimeBins*kNSeedPlanes);
	memset(mcl, 0, sizeof(Int_t)*kNTimeBins*kNSeedPlanes);
	Int_t nused = 0;
	AliTRDchamberTimeBin *layers = chamber->GetTB(0);
	for(Int_t ipl = 0; ipl < kNSeedPlanes; ipl++){
		for(Int_t ils = 0; ils < fTimeBinsPerPlane; ils++){
			position = planes[ipl]*fTimeBinsPerPlane + ils;
			ncl[ipl * fTimeBinsPerPlane + ils] = layers[position].GetNClusters();
			nused = 0;
			for(Int_t icl = 0; icl < ncl[ipl * fTimeBinsPerPlane + ils]; icl++)
				if((layers[position].GetCluster(icl))->IsUsed()) nused++;
			ncl[ipl * fTimeBinsPerPlane + ils] -= nused;
		}
	}
	// Declaration of quartils:
	//Double_t qvals[3] = {0.0, 0.0, 0.0};
	//Double_t qprop[3] = {0.16667, 0.5, 0.83333};
	// Iterations
	Int_t counter;
	Double_t *array;
	Int_t *limit;
	Int_t nLayers = fTimeBinsPerPlane * kNSeedPlanes;
	for(Int_t iter = 0; iter < 2; iter++){
		array = (iter == 0) ? &ncl[0] : &mcl[0];
		limit = (iter == 0) ? &nLayers : &counter;
		counter = 0;
		if(iter == 1){
			for(Int_t i = 0; i < fTimeBinsPerPlane *kNSeedPlanes; i++){
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
Int_t AliTRDtrackerV1::GetSeedingLayers(AliTRDtrackingChamber *chamber, Double_t *params)
{
  //
  // Algorithm to find optimal seeding layer
  // Layers inside one sigma region (given by Quantiles) are sorted
  // according to their difference.
  // All layers outside are sorted according t
  //
  // Parameters:
  //     - Array of AliTRDchamberTimeBins (in the current plane !!!)
  //     - Container for the Indices of the seeding Layer candidates
  //
  // Output:
  //     - Number of Layers inside the 1-sigma-region
  //
  // The optimal seeding layer should contain the mean number of
  // custers in the layers in one chamber.
  //

	//printf("Params: %3.3f, %3.3f, %3.3f\n", params[0], params[1], params[2]);
	const Int_t kMaxClustersLayer = AliTRDchamberTimeBin::kMaxClustersLayer;
	Int_t ncl[kNTimeBins], indices[kNTimeBins], bins[kMaxClustersLayer];
	memset(ncl, 0, sizeof(Int_t)*kNTimeBins);
	memset(indices, 0, sizeof(Int_t)*kNTimeBins);
	memset(bins, 0, sizeof(Int_t)*kMaxClustersLayer);
	
	AliTRDchamberTimeBin *layers = chamber->GetTB(0);
	Int_t nused = 0;
	for(Int_t ils = 0; ils < fTimeBinsPerPlane; ils++){
		ncl[ils] = layers[ils].GetNClusters();
		nused = 0;
		for(Int_t icl = 0; icl < ncl[ils]; icl++)
			if((layers[ils].GetCluster(icl))->IsUsed()) nused++;
		ncl[ils] -= nused;
	}
	
	Float_t mean = params[1];
	for(Int_t ils = 0; ils < fTimeBinsPerPlane; ils++){
		memmove(indices + bins[ncl[ils]+1] + 1, indices + bins[ncl[ils]+1], sizeof(Int_t)*(fTimeBinsPerPlane - ils));
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
			position = (Int_t)(gRandom->Rndm()*(fTimeBinsPerPlane-1));
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
AliTRDcluster *AliTRDtrackerV1::FindSeedingCluster(AliTRDtrackingChamber *chamber, AliTRDseedV1* reference) const
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
  // Imput parameters: - layers: array of AliTRDchamberTimeBins (in one chamber!!!)
  //                   - sfit: the reference
  //
  // Output:           - the best seeding cluster
  //

	
	// distances as squared distances
	Int_t index = 0;
	Float_t ypos = 0.0, zpos = 0.0, distance = 0.0, nearestDistance =100000.0; 
	ypos = reference->GetYref(0);
	zpos = reference->GetZref(0);
	AliTRDcluster *currentBest = 0x0, *temp = 0x0;
	AliTRDchamberTimeBin *layers = chamber->GetTB(0);
	for(Int_t ils = 0; ils < fTimeBinsPerPlane; ils++){
		// Reference positions
// 		ypos = reference->GetYat(layers[ils].GetX());
// 		zpos = reference->GetZat(layers[ils].GetX());
		index = layers[ils].SearchNearestCluster(ypos, zpos, AliTRDReconstructor::RecoParam()->GetRoad2y(), AliTRDReconstructor::RecoParam()->GetRoad2z());
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

//____________________________________________________________________
AliCluster* AliTRDtrackerV1::GetCluster(Int_t idx) const
{
	Int_t ncls = fClusters->GetEntriesFast();
	return idx >= 0 || idx < ncls ? (AliCluster*)fClusters->UncheckedAt(idx) : 0x0;
}


//_____________________________________________________________________________
Int_t AliTRDtrackerV1::Freq(Int_t n, const Int_t *inlist
                        , Int_t *outlist, Bool_t down)
{    
  //
  // Sort eleements according occurancy 
  // The size of output array has is 2*n 
  //

  if (n <= 0) {
    return 0;
  }

  Int_t *sindexS = new Int_t[n];   // Temporary array for sorting
  Int_t *sindexF = new Int_t[2*n];   
  for (Int_t i = 0; i < n; i++) {
    sindexF[i] = 0;
  }

  TMath::Sort(n,inlist,sindexS,down); 
 
  Int_t last     = inlist[sindexS[0]];
  Int_t val      = last;
  sindexF[0]     = 1;
  sindexF[0+n]   = last;
  Int_t countPos = 0;

  // Find frequency
  for (Int_t i = 1; i < n; i++) {
    val = inlist[sindexS[i]];
    if (last == val) {
      sindexF[countPos]++;
    }
    else {      
      countPos++;
      sindexF[countPos+n] = val;
      sindexF[countPos]++;
      last                = val;
    }
  }
  if (last == val) {
    countPos++;
  }

  // Sort according frequency
  TMath::Sort(countPos,sindexF,sindexS,kTRUE);

  for (Int_t i = 0; i < countPos; i++) {
    outlist[2*i  ] = sindexF[sindexS[i]+n];
    outlist[2*i+1] = sindexF[sindexS[i]];
  }

  delete [] sindexS;
  delete [] sindexF;
  
  return countPos;

}
