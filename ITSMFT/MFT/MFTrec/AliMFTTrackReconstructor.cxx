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

#include "TClonesArray.h"
#include "TMath.h"

#include "AliCodeTimer.h"
#include "AliLog.h"

#include "AliMFTTrackReconstructor.h"
#include "AliMFTTrackParam.h"
#include "AliMFTTrack.h"
#include "AliMFTTrackExtrap.h"
#include "AliMFTCATrack.h"
#include "AliMFTCACell.h"
#include "AliMFTConstants.h"

/// \cond CLASSIMP
ClassImp(AliMFTTrackReconstructor); // Class implementation in ROOT context
/// \endcond


//=============================================================================================

AliMFTTrackReconstructor::AliMFTTrackReconstructor():TObject()
{
	/// Default constructor
	
	// set the magnetic field for track extrapolations
	AliMFTTrackExtrap::SetField();

}


//=============================================================================================


AliMFTTrackReconstructor::~AliMFTTrackReconstructor() {
	
}
//__________________________________________________________________________
Bool_t AliMFTTrackReconstructor::TraceTrack(AliMFTTrack *currentTrack ){
	
	Bool_t extrapStatus = kTRUE;
	Double_t addChi2TrackAtCluster;

	AliMFTCATrack * currentCATrack = currentTrack->GetCATrack();
	Int_t nCells = currentCATrack->GetNcells();
	
	if (nCells < 2) return kFALSE; // Skip tracks wih less than 2 cells
	
	// Initiate the seed  starting from the last cells (the more downstream one)
	AliMFTCACell * caCell = currentCATrack->GetCell(0);
	caCell->PrintCell("");
	
	// Evaluate the track sign and Pt
	currentCATrack->EvalSignedPt();
	
	Double_t *caHit1, *caHit2;
	caHit1 = caCell->GetHit1();
	caHit2 = caCell->GetHit2();
		
	Double_t dX = caHit1[0] - caHit2[0];
	Double_t dY = caHit1[1] - caHit2[1];
	Double_t dZ = caHit1[2] - caHit2[2];
	Double_t dr = TMath::Sqrt(dX*dX+dY*dY);
	Double_t slopeX_Z = dX / dZ;
	Double_t slopeY_Z = dY / dZ;
	Double_t slopeY_R = dY / dr;
	Double_t slope2 = slopeX_Z*slopeX_Z + slopeY_Z*slopeY_Z;

	// Inverse  momentum
		Double_t inversePt = 1./currentCATrack->GetPt(); // Signed Pt estimation

	// Set track parameters at second cluster
	AliMFTTrackParam* trackParamAtLastCluster = (AliMFTTrackParam*) currentTrack->GetTrackParamAtCluster()->First();

	trackParamAtLastCluster->SetX(caHit2[0]);
	trackParamAtLastCluster->SetY(caHit2[1]);
	trackParamAtLastCluster->SetZ(caHit2[2]);
	trackParamAtLastCluster->SetSlopeX(slopeX_Z);
	trackParamAtLastCluster->SetSlopeY(slopeY_Z);
	trackParamAtLastCluster->SetInverseTransverseMomentum(inversePt);
	
	Double_t inverseMomentum;
	
	// Compute and set track parameters covariances at first cluster
	TMatrixD paramCov(5,5);
	paramCov.Zero();
	Double_t errX2 = AliMFTConstants::kXPixelPitch * AliMFTConstants::kXPixelPitch / 12.;
	Double_t errY2 = AliMFTConstants::kYPixelPitch * AliMFTConstants::kYPixelPitch / 12.;


	paramCov(0,0) = errX2;
	paramCov(1,1) = errY2;
	paramCov(2,2) = ( 1001. * errX2 )/ dZ / dZ;// Weight 1 for the last cluster and 1000 for the first one. the first cluster error will be taken into account at the first step of the Kalman filter
	paramCov(2,0) = -errX2 / dZ;
	paramCov(0,2) = paramCov(2,0);
	
	paramCov(3,3) = ( 1001. * errY2 )/ dZ / dZ;// Weight 1 for the last cluster and 1000 for the first one. the first cluster error will be taken into account at the first step of the Kalman filter
	paramCov(3,1) = -errY2 / dZ;
	paramCov(1,3) = paramCov(3,1);

	
	
	//take 100% error on inverse momentum
	
	paramCov(4,4) =   inversePt*inversePt
		*(0.1*0.1+ (slopeX_Z*slopeX_Z *errX2 *1001. +	slopeY_Z*slopeY_Z *errY2 *1001.)/dZ/dZ / slope2 / slope2);
	paramCov(4,0) =  inversePt / dZ / slope2 * errX2 * slopeX_Z;
	paramCov(4,1) =  inversePt / dZ / slope2 * errY2 * slopeY_Z;
	paramCov(0,4) = paramCov(4,0);
	paramCov(1,4) = paramCov(4,1);
	paramCov(4,2) = - inversePt  / slope2 * slopeX_Z * paramCov(2,2);
	paramCov(4,3) = - inversePt  / slope2 * slopeY_Z * paramCov(3,3);
	paramCov(2,4) = paramCov(4,2);
	paramCov(3,4) = paramCov(4,3);

	trackParamAtLastCluster->SetCovariances(paramCov);
	AliInfo("Starting Covariance Matrix");
	paramCov.Print();
	// Reset the track chi2
	trackParamAtLastCluster->SetTrackChi2(0.);
	trackParamAtLastCluster->Print("FULL");
	

	// Follow the track going upstream
	AliMFTTrackParam * startingTrackParam = NULL;
	startingTrackParam = trackParamAtLastCluster;
	AliMFTTrackParam * trackParamAtCluster = (AliMFTTrackParam*) currentTrack->GetTrackParamAtCluster()->After(startingTrackParam);
	currentTrack->SetMCLabel(currentCATrack->GetMCindex());
	
	for (Int_t iCell = 0  ; iCell < nCells; iCell++) {
		caCell = currentCATrack->GetCell(iCell);
		caHit1 = caCell->GetHit1();
		caHit2 = caCell->GetHit2();

		caCell->PrintCell("MC");

		// reset track parameters and their covariances
		trackParamAtCluster->SetParameters(startingTrackParam->GetParameters());
		trackParamAtCluster->SetZ(startingTrackParam->GetZ());
		trackParamAtCluster->SetCovariances(startingTrackParam->GetCovariances());

		// add MCS effect
		extrapStatus = AddMCSEffect(caCell, trackParamAtCluster);
		
		// extrapolation to the plane of the cluster attached to the current trackParamAtCluster (update the propagator)
		if (!AliMFTTrackExtrap::ExtrapToZCov(trackParamAtCluster, caHit1[2],
																					kFALSE)) extrapStatus = kFALSE;
		
		// Compute new track parameters using kalman filter
		addChi2TrackAtCluster = RunKalmanFilter(*trackParamAtCluster);

		// Update the track chi2
		currentTrack->SetChi2(currentTrack->GetChi2() + addChi2TrackAtCluster);

		trackParamAtCluster->SetTrackChi2(currentTrack->GetChi2());
		trackParamAtCluster->SetLocalChi2(addChi2TrackAtCluster);
		AliInfo("Param at  cluster");
		trackParamAtCluster->Print("FULL");

		
		// prepare next step
		startingTrackParam = trackParamAtCluster;
		trackParamAtCluster = (AliMFTTrackParam*) (currentTrack->GetTrackParamAtCluster()->After(startingTrackParam));
		
		if( (iCell == nCells-1) && trackParamAtCluster )
			AliWarning("Reaching last cell but still some AliMFTTrackParameter objects in the array ...");

	}
	
	currentTrack->SetChi2(((AliMFTTrackParam*) currentTrack->GetTrackParamAtCluster()->Last())->GetTrackChi2()/(nCells+1));
	currentTrack->SetPhi(((AliMFTTrackParam*) currentTrack->GetTrackParamAtCluster()->Last())->GetPhi());
	currentTrack->SetTheta(((AliMFTTrackParam*) currentTrack->GetTrackParamAtCluster()->Last())->GetTheta());
	currentTrack->SetP(((AliMFTTrackParam*) currentTrack->GetTrackParamAtCluster()->Last())->P());
	AliInfo(Form("Input Pt %f ",currentCATrack->GetPt()));

	currentTrack->SetPt(currentCATrack->GetPt());
	currentTrack->Print();
	return extrapStatus;
}
//__________________________________________________________________________
Bool_t AliMFTTrackReconstructor::AddMCSEffect(AliMFTCACell *currentCell,  AliMFTTrackParam *trackParam)
{
	Bool_t returnStatus = kTRUE;
	Int_t *planeIDs = currentCell->GetLayers();
	Int_t startingPlaneID = planeIDs[1];
	Int_t startingDiskID = startingPlaneID/2;
	AliInfo(Form("Layer ID = %d ; Length = %d ",startingPlaneID,currentCell->GetLength()));
	
	if(currentCell->GetLength()==1){
		// No MCS effect between front face and back face of 2 different disks (it is the air between two disks, neglect it for now)
		if(startingPlaneID%2==1){
			AliInfo(Form("Add MCS of Disk %d ",startingDiskID));

			AliMFTTrackExtrap::AddMCSEffect(trackParam,-1.4,AliMFTConstants::DiskThicknessInX0(startingDiskID));
		}
		return kTRUE;
	}
	

	// Applying MCS taking into account missing planes if any
	Int_t currentPlaneID = startingPlaneID;
	AliInfo(Form("Layer IDs %d - %d ",planeIDs[0], planeIDs[1]));
	while (currentPlaneID != planeIDs[0]) {
		// extrapolation to the missing chamber (update the propagator)
		// add MCS effect if second hit is in the back plane face of a disk
		if(currentPlaneID%2==1){
			AliInfo(Form("Add MCS of Disk %d ",currentPlaneID/2));
			AliMFTTrackExtrap::AddMCSEffect(trackParam,-1.4,AliMFTConstants::DiskThicknessInX0(currentPlaneID/2));

		}
		AliInfo(Form("Extrapolate to Plane ID %d ",currentPlaneID-1));
		if (!AliMFTTrackExtrap::ExtrapToZCov(trackParam, AliMFTConstants::DefaultPlaneZ(currentPlaneID-1),																					 kFALSE)) returnStatus = kFALSE;
		AliInfo(Form("After extrapolation to Z = %f",trackParam->GetZ()));
		trackParam->Print("FULL");

//		// add MCS effect
//		AliInfo(Form("Add MCS of Disk %d ",currentPlaneID/2));
//		AliMFTTrackExtrap::AddMCSEffect(trackParam,-1.4,AliMFTConstants::DiskThicknessInX0(startingDiskID));
		
		currentPlaneID--;
	}

	return returnStatus;

}
//__________________________________________________________________________
void AliMFTTrackReconstructor::EventReconstruct(TClonesArray *fMFTTracks )
{
	/// To reconstruct one event
	AliDebug(1,"");
	AliCodeTimerAuto("",0);

	// Stop tracking if no track candidate found
	
	Int_t nTracks = fMFTTracks->GetEntriesFast();
	if (nTracks == 0) return;
	

	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
		AliInfo(Form("Treating Track %d",iTrack));
		TraceTrack((AliMFTTrack *) fMFTTracks->UncheckedAt(iTrack));
		((AliMFTTrack *) fMFTTracks->UncheckedAt(iTrack))->Print("");
	}
	
	

	AliInfo("I am HERE ....");
	// Reset array of tracks
//	ResetTracks();
//	
//	// Look for candidates from clusters in stations(1..) 4 and 5 (abort in case of failure)
//	if (!MakeTrackCandidates(clusterStore)) return;
//	
//	// Look for extra candidates from clusters in stations(1..) 4 and 5 (abort in case of failure)
//	if (GetRecoParam()->MakeMoreTrackCandidates()) {
//		if (!MakeMoreTrackCandidates(clusterStore)) return;
//	}
//	
//
//	// Follow tracks in stations(1..) 3, 2 and 1 (abort in case of failure)
//	if (!FollowTracks(clusterStore)) return;
//	
//	// Complement the reconstructed tracks
//	if (GetRecoParam()->ComplementTracks()) {
//		if (ComplementTracks(clusterStore)) RemoveIdenticalTracks();
//	}
//	
//	// Improve the reconstructed tracks
//	if (GetRecoParam()->ImproveTracks()) ImproveTracks();
//	
//	// Remove connected tracks
//	RemoveConnectedTracks(3, 4, kFALSE);
//	RemoveConnectedTracks(2, 2, kFALSE);
//	if (GetRecoParam()->RemoveConnectedTracksInSt12()) RemoveConnectedTracks(0, 1, kFALSE);
//	
//	// Fill AliMUONTrack data members
//	Finalize();
//	if (!GetRecoParam()->RemoveConnectedTracksInSt12()) TagConnectedTracks(0, 1, kTRUE);
//	
//	// Make sure there is no bad track left
//	RemoveBadTracks();
//	
//	// Refit the reconstructed tracks with a different resolution for mono-cathod clusters
//	if (GetRecoParam()->DiscardMonoCathodClusters()) DiscardMonoCathodClusters();
//	
//	// Add tracks to MUON data container
//	for (Int_t i=0; i<fNRecTracks; ++i)
//	{
//		AliMUONTrack * track = (AliMUONTrack*) fRecTracksPtr->At(i);
//		track->SetUniqueID(i+1);
//		trackStore.Add(*track);
//	}
	
}

//__________________________________________________________________________
Double_t AliMFTTrackReconstructor::RunKalmanFilter(AliMFTTrackParam &trackParamAtCluster)
{
	/// Compute new track parameters and their covariances including new cluster using kalman filter
	/// return the additional track chi2
	AliInfo("Enter RunKalmanFilter");
	
	// Get actual track parameters (p)
	TMatrixD param(trackParamAtCluster.GetParameters());
	
	// Get new cluster parameters (m)
	TMatrixD clusterParam(5,1);
	clusterParam.Zero();
	clusterParam(0,0) = trackParamAtCluster.GetClusterX();
	clusterParam(1,0) = trackParamAtCluster.GetClusterY();
	
	AliInfo(Form("Cluster X,Y : %f   %f ",clusterParam(0,0), clusterParam(1,0)));
	

	
	// Compute the actual parameter weight (W)
	TMatrixD paramWeight(trackParamAtCluster.GetCovariances());
	AliInfo("actual parameter covariance matrix");
	paramWeight.Print();
	AliInfo(Form("paramWeight.Determinant  = %e ",paramWeight.Determinant()));

	if (paramWeight.Determinant() != 0) {
		paramWeight.Invert();
	} else {
		AliWarning(" paramWeight Determinant = 0");
		return 1.e6;
	}
	AliInfo("actual parameter weight");
	paramWeight.Print();

	// Compute the new cluster weight (U)
	TMatrixD clusterWeight(5,5);
	clusterWeight.Zero();
	Double_t errX2 = AliMFTConstants::kXPixelPitch * AliMFTConstants::kXPixelPitch / 12.;
	Double_t errY2 = AliMFTConstants::kYPixelPitch * AliMFTConstants::kYPixelPitch / 12.;
	clusterWeight(0,0) = 1./errX2;
	clusterWeight(1,1) = 1./errY2;
	AliInfo("new cluster weight");
	clusterWeight.Print();
	// Compute the new parameters covariance matrix ( (W+U)^-1 )
	TMatrixD newParamCov(paramWeight,TMatrixD::kPlus,clusterWeight);
	AliInfo("new parameters weight");
	newParamCov.Print();

	if (newParamCov.Determinant() != 0) {
		newParamCov.Invert();
	} else {
		AliWarning(" newParamCov Determinant = 0");
		return 1.e6;
	}

	// Save the new parameters covariance matrix
	AliInfo("new parameters covariance matrix");
	newParamCov.Print();

	trackParamAtCluster.SetCovariances(newParamCov);
	
	// Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
	TMatrixD tmp(clusterParam,TMatrixD::kMinus,param);
	TMatrixD tmp2(clusterWeight,TMatrixD::kMult,tmp); // U(m-p)
	TMatrixD newParam(newParamCov,TMatrixD::kMult,tmp2); // ((W+U)^-1)U(m-p)
	AliInfo("delta parameters");
	newParam.Print();
	AliInfo("old parameters");
	param.Print();

	newParam += param; // ((W+U)^-1)U(m-p) + p
	
	// Save the new parameters
	trackParamAtCluster.SetParameters(newParam);
	
	// Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
	tmp = newParam; // p'
	tmp -= param; // (p'-p)
	TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp); // W(p'-p)
	TMatrixD addChi2Track(tmp,TMatrixD::kTransposeMult,tmp3); // ((p'-p)^-1)W(p'-p)
	tmp = newParam; // p'
	tmp -= clusterParam; // (p'-m)
	TMatrixD tmp4(clusterWeight,TMatrixD::kMult,tmp); // U(p'-m)
	addChi2Track += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
	
	return addChi2Track(0,0);
	
}

