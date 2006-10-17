//========================================================================
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
//                                                                        
// Author: The ALICE Off-line Project.                                    
// Contributors are mentioned in the code where appropriate.              
//                                                                        
// Permission to use, copy, modify and distribute this software and its   
// documentation strictly for non-commercial purposes is hereby granted   
// without fee, provided that the above copyright notice appears in all   
// copies and that both the copyright notice and this permission notice   
// appear in the supporting documentation. The authors make no claims     
// about the suitability of this software for any purpose. It is          
// provided "as is" without express or implied warranty.                  
//======================================================================== 
//                       
//                       Class AliEMCALTracker 
//                      -----------------------
// Implementation of the track matching method between barrel tracks and
// EMCAL clusters.
// Besides algorithm implementation, some cuts are required to be set
// in order to define, for each track, an acceptance window where clusters
// are searched to find best match (if any).
// The class accepts as input an ESD container, and works directly on it,
// simply setting, for each of its tracks, the fEMCALindex flag, for each
// track which is matched to a cluster.
// In order to use method, one must launch PropagateBack().
//
// ------------------------------------------------------------------------
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//=========================================================================

#include <Riostream.h>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TList.h>
#include <TString.h>
#include <TVector3.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliKalmanTrack.h"
#include "AliEMCALRecPoint.h"
#include "AliRunLoader.h"
#include "AliEMCALTrack.h"
#include "AliEMCALLoader.h"
#include "AliEMCALTracker.h"

ClassImp(AliEMCALTracker)
//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALTracker() 
  : AliTracker(),
    fNPropSteps(0),
    fTrackCorrMode(kTrackCorrNone),
    fCutX(50.0),
    fCutY(50.0),
    fCutZ(50.0),
    fCutAlphaMin(-200.0),
    fCutAlphaMax(200.0),
    fCutAngle(100.0),
    fMaxDist(100.0),
    fRho(1.0),
    fX0(1.0),
    fTracks(0),
    fClusters(0),
    fMatches(0)
{
	//
	// Default constructor.
	// Initializes al simple data members to default values,
	// and all collections to NULL.
	// Output file name is set to a default value.
	//
}
//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALTracker(const AliEMCALTracker& copy) 
  : AliTracker(),
    fNPropSteps(copy.fNPropSteps),
    fTrackCorrMode(copy.fTrackCorrMode),
    fCutX(copy.fCutX),
    fCutY(copy.fCutY),
    fCutZ(copy.fCutZ),
    fCutAlphaMin(copy.fCutAlphaMin),
    fCutAlphaMax(copy.fCutAlphaMax),
    fCutAngle(copy.fCutAngle),
    fMaxDist(copy.fMaxDist),
    fRho(copy.fRho),
    fX0(copy.fX0),
    fTracks((TObjArray*)copy.fTracks->Clone()),
    fClusters((TObjArray*)copy.fClusters->Clone()),
    fMatches((TList*)copy.fMatches->Clone())
{
	//
	// Copy constructor
	// Besides copying all parameters, duplicates all collections.
	//
}
//
//------------------------------------------------------------------------------
//
AliEMCALTracker& AliEMCALTracker::operator=(const AliEMCALTracker& copy)
{
	//
	// Assignment operator.
	// Besides copying all parameters, duplicates all collections.	
	//
	
	fCutX = copy.fCutX;
	fCutY = copy.fCutY;
	fCutZ = copy.fCutZ;
	fCutAlphaMin = copy.fCutAlphaMin;
	fCutAlphaMax = copy.fCutAlphaMax;
	fCutAngle = copy.fCutAngle;
	fMaxDist = copy.fMaxDist;
	
	fTracks = (TObjArray*)copy.fTracks->Clone();
	fClusters = (TObjArray*)copy.fClusters->Clone();
	fMatches = (TList*)copy.fMatches->Clone();
	
	return (*this);
}
//
//------------------------------------------------------------------------------
//
void AliEMCALTracker::Clear(Option_t* option)
{
	//
	// Clearing method
	// Clears all specified arrays and the containers themselves.
	//

	TString opt(option);
	Bool_t clearTracks = opt.Contains("TRACKS");
	Bool_t clearClusters = opt.Contains("CLUSTERS");
	Bool_t clearMatches = opt.Contains("MATCHES");
	if (opt.Contains("ALL")) {
		clearTracks = kTRUE;
		clearClusters = kTRUE;
		clearMatches = kTRUE;
	}
	
	if (fTracks != 0x0 && clearTracks) {
		if (!fTracks->IsEmpty()) fTracks->Delete();
		delete fTracks;
		fTracks = 0;
	}
	if (fClusters != 0x0 && clearClusters) {
		if (!fClusters->IsEmpty()) fClusters->Delete();
		delete fClusters;
		fClusters = 0;
	}
	if (fMatches != 0x0 && clearMatches) {
		if (!fMatches->IsEmpty()) fMatches->Delete();
		delete fMatches;
		fMatches = 0;
	}
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::LoadClusters(TTree *cTree) 
{
	//
	// Load EMCAL clusters in the form of AliEMCALRecPoint,
	// from simulation temporary files.
	// (When included in reconstruction chain, this method is used automatically)
	//
	
	Clear("CLUSTERS");

	TBranch *branch = cTree->GetBranch("EMCALECARP");
	if (!branch) {
		AliError("can't get the branch with the EMCAL clusters");
		return 1;
	}
	
	TClonesArray dummy("AliEMCALRecPoint", 10000);
	TClonesArray *clusters = &dummy;
	branch->SetAddress(&clusters);
	Int_t nClusters = (Int_t)clusters->GetEntries();
	
	cTree->GetEvent(0);
	fClusters = new TObjArray(0);
	for (Int_t i = 0; i < nClusters; i++) {
		AliEMCALRecPoint *cluster = (AliEMCALRecPoint*)clusters->At(i);
		if (!cluster) continue;
		if (cluster->GetClusterType() != AliESDCaloCluster::kClusterv1) continue;
		AliEMCALMatchCluster *matchCluster = new AliEMCALMatchCluster(i, cluster);
		fClusters->AddLast(matchCluster);
	}
	if (fClusters->IsEmpty()) {
		AliError("No clusters collected");
		return 1;
	}
	
	AliInfo(Form("Collected %d clusters", fClusters->GetEntries()));

	return 0;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::LoadClusters(AliESD *esd) 
{
	//
	// Load EMCAL clusters in the form of AliESDCaloClusters,
	// from an AliESD object.
	//

	// make sure that tracks/clusters collections are empty
	Clear("CLUSTERS");
	
	Int_t start = esd->GetFirstEMCALCluster();
	Int_t nClustersEMC = esd->GetNumberOfEMCALClusters();
	Int_t end = start + nClustersEMC;
	
	fClusters = new TObjArray(0);
		
	Int_t i;
	for (i = start; i < end; i++) {
		AliESDCaloCluster *cluster = esd->GetCaloCluster(i);
		if (!cluster) continue;
		if (cluster->GetClusterType() != AliESDCaloCluster::kClusterv1) continue;
		AliEMCALMatchCluster *matchCluster = new AliEMCALMatchCluster(i, cluster);
		fClusters->AddLast(matchCluster);
	}
	if (fClusters->IsEmpty()) {
		AliError("No clusters collected");
		return 1;
	}
	
	AliInfo(Form("Collected %d clusters", fClusters->GetEntries()));

	return 0;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::LoadTracks(AliESD *esd)
{
	//
	// Load ESD tracks.
	//
	
	Clear("TRACKS");
	
	Int_t nTracks = esd->GetNumberOfTracks();
	fTracks = new TObjArray(0);
	
	Int_t i, j;
	Bool_t isKink;
	Double_t alpha; 
	for (i = 0; i < nTracks; i++) {
		AliESDtrack *esdTrack = esd->GetTrack(i);
		// set by default the value corresponding to "no match"
		esdTrack->SetEMCALcluster(-99999);
//		if (esdTrack->GetLabel() < 0) continue;
//		if (!(esdTrack->GetStatus() & AliESDtrack::kTOFout)) continue;
		isKink = kFALSE;
		for (j = 0; j < 3; j++) {
			if (esdTrack->GetKinkIndex(j) > 0) isKink = kTRUE;
		}
		if (isKink) continue;
		AliEMCALTrack *track = new AliEMCALTrack(*esdTrack);
		// check alpha and reject the tracks which fall outside EMCAL acceptance
		alpha = track->GetAlpha() * TMath::RadToDeg();
		if (alpha >  -155.0 && alpha < 67.0) {
			delete track;
			continue;
		}
		track->SetSeedIndex(i);
		track->SetSeedLabel(esdTrack->GetLabel());
		fTracks->AddLast(track);
	}
	if (fTracks->IsEmpty()) {
		AliError("No tracks collected");
		return 1;
	}
	
	AliInfo(Form("Collected %d tracks", fTracks->GetEntries()));

	return 0;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::PropagateBack(AliESD* esd)
{
	//
	// Main operation method.
	// Gets external AliESD containing tracks to be matched.
	// After executing match finding, stores in the same ESD object all infos
	// and releases the object for further reconstruction steps.
	//
	
	if (!esd) {
		AliError("NULL ESD passed");
		return 1;
	}
	
	// step 1: 
	// if cluster array is empty, cluster are collected
	// from the passed ESD, and work is done with ESDCaloClusters
	Int_t okLoadClusters;
	if (!fClusters || (fClusters && fClusters->IsEmpty())) {
		AliInfo("Cluster array is empty. Loading clusters...");
		okLoadClusters = LoadClusters(esd);
		if (okLoadClusters) return 2;
	}
	
	// step 2:
	// collect ESD tracks
	Int_t okLoadTracks = LoadTracks(esd);
	if (okLoadTracks) return 3;
	
	// step 3:
	// each track is propagated to the "R" position of each cluster.
	// The closest cluster is assigned as match.
	// IF no clusters lie within the maximum allowed distance, no matches are assigned.
	Int_t nMatches = CreateMatches();
	if (!nMatches) {
		AliInfo("No good matches found.");
		return 4;
	}
	else {
		AliInfo(Form("Found %d matches", nMatches));
	}
	
	// step 4:
	// when more than 1 track share the same matched cluster, only the closest one is kept.
	Int_t nRemoved = SolveCompetitions();
	AliInfo(Form("Removed %d duplicate matches", nRemoved));
	if (nRemoved >= nMatches) {
		AliError("Removed ALL matches! Check the algorithm or data. Nothing to save");
		return 5;
	}
	
	// step 5:
	// save obtained information setting the 'fEMCALindex' field of AliESDtrack object
	Int_t nSaved = 0, trackID;
	TListIter iter(fMatches);
	AliEMCALMatch *match = 0;
	while ( (match = (AliEMCALMatch*)iter.Next()) ) {
		if (!match->CanBeSaved()) continue;
		AliEMCALTrack *track = (AliEMCALTrack*)fTracks->At(match->GetIndexT());
		AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(match->GetIndexC());
		trackID = track->GetSeedIndex();
		AliESDtrack *esdTrack = esd->GetTrack(trackID);
		if (!esdTrack) continue;
		if (esdTrack->GetLabel() == cluster->Label()) {
			esdTrack->SetEMCALcluster(cluster->Index());
		}
		else {
			esdTrack->SetEMCALcluster(-cluster->Index());
		}
		nSaved++;
	}
	/*
	AliEMCALTrack *track = 0;
	TObjArrayIter tracks(fTracks);
	while ( (track = (AliEMCALTrack*)tracks.Next()) ) {
		trackID = track->GetSeedIndex();
		clusterID = track->GetMatchedClusterIndex();
		AliESDtrack *esdTrack = esd->GetTrack(trackID);
		if (!esdTrack) continue;
		if (clusterID < 0) {
			esdTrack->SetEMCALcluster(-99999);
		}
		else {
			AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(clusterID);
			if (!cluster) continue;
			if (esdTrack->GetLabel() == cluster->Label()) {
				nGood++;
				esdTrack->SetEMCALcluster(cluster->Index());
			}
			else {
				esdTrack->SetEMCALcluster(-cluster->Index());
			}
			nSaved++;
		}
	}
	*/
	AliInfo(Form("Saved %d matches", nSaved));

	return 0;
}
//
//------------------------------------------------------------------------------
//
void AliEMCALTracker::SetTrackCorrectionMode(Option_t *option)
{
	//
	// Set track correction mode
	// gest the choice in string format and converts into 
	// internal enum
	//
	
	TString opt(option);
	opt.ToUpper();
	
	if (!opt.CompareTo("NONE")) {
		fTrackCorrMode = kTrackCorrNone;
	}
	else if (!opt.CompareTo("MMB")) {
		fTrackCorrMode = kTrackCorrMMB;
	}
	else if (!opt.CompareTo("FIXED")) {
		fTrackCorrMode = kTrackCorrFixed;
	}
	else {
		cerr << "E-AliEMCALTracker::SetTrackCorrectionMode '" << option << "': Unrecognized option" << endl;
	}
}

//
//------------------------------------------------------------------------------
//
Double_t AliEMCALTracker::AngleDiff(Double_t angle1, Double_t angle2)
{
	// 
	// [PRIVATE]
	// Given two angles in radiants, it converts them in the range 0-2pi
	// then computes their true difference, i.e. if the difference a1-a2
	// results to be larger than 180 degrees, it returns 360 - diff.
	//
	
	if (angle1 < 0.0) angle1 += TMath::TwoPi();
	if (angle1 > TMath::TwoPi()) angle1 -= TMath::TwoPi();
	if (angle2 < 0.0) angle2 += TMath::TwoPi();
	if (angle2 > TMath::TwoPi()) angle2 -= TMath::TwoPi();
	
	Double_t diff = TMath::Abs(angle1 - angle2);
	if (diff > TMath::Pi()) diff = TMath::TwoPi() - diff;
	
	if (angle2 > angle1) diff = -diff;
	
	return diff;
}
//
//------------------------------------------------------------------------------
//
Double_t AliEMCALTracker::CheckPair
(AliEMCALTrack *track, AliEMCALMatchCluster *cl)
{
	//
	// Given a track and a cluster,
	// propagates the first to the radius of the second.
	// Then, checks the propagation point against all cuts.
	// If at least a cut is not passed, a valuer equal to 
	// twice the maximum allowed distance is passed (so the value returned
	// will not be taken into account when creating matches)
	//
	
	// TEMP
	Bool_t isTrue = kFALSE;
//	if (tr->GetSeedLabel() == cl->Label()) {
//		isTrue = kTRUE;
//		cout << "TRUE MATCH!!!" << endl;
//	}
	
	// copy track into temporary variable
	AliEMCALTrack *tr = new AliEMCALTrack(*track);
	
	Double_t distance = 2.0 * fMaxDist;
	
	// check against cut on difference 'alpha - phi'
	Double_t phi = TMath::ATan2(cl->Y(), cl->X());
	phi = AngleDiff(phi, tr->GetAlpha());
	if (phi < fCutAlphaMin || phi > fCutAlphaMax) return distance;
	
	// try to propagate to cluster radius
	// (return the 'distance' value if it fails)
	Double_t pos[3], &x = pos[0], &y = pos[1], &z = pos[2];
	Double_t x0, rho;
	tr->GetXYZ(pos);
	Double_t rt = TMath::Sqrt(x*x + y*y);
	Double_t rc = TMath::Sqrt(cl->X()*cl->X() + cl->Y()*cl->Y());
	
	if (fTrackCorrMode == kTrackCorrMMB) {
		Double_t pos1[3], pos2[3], param[6];
		pos1[0] = x;
		pos1[1] = y;
		pos1[2] = z;
		pos2[0] = cl->X();
		pos2[1] = cl->Y();
		pos2[2] = cl->Z();
		AliKalmanTrack::MeanMaterialBudget(pos1, pos2, param);
		rho = param[0];
		x0 = param[1];
	}
	else if (fTrackCorrMode == kTrackCorrFixed) {
		rho = fRho;
		x0 = fX0;
	}
	else {
		rho = 0.0;
		x0 = 0.0;
	}
	if (fNPropSteps) {
		Int_t i;
		Double_t r;
		cout.setf(ios::fixed);
		cout.precision(5);
		if (isTrue) cout << "Init : " << rt << ' ' << x << ' ' << y << ' ' << z << endl;
		for (i = 0; i < fNPropSteps; i++) {
			r = rt + (rc - rt) * ((Double_t)(i+1)/(Double_t)fNPropSteps);
			if (!tr->PropagateTo(r, x0, rho)) return distance;
			tr->GetXYZ(pos);
			if (isTrue) cout << "Step : " << r << ' ' << x << ' ' << y << ' ' << z << endl;
		}
		if (isTrue) cout << "Clstr: " << rc << ' ' << cl->X() << ' ' << cl->Y() << ' ' << cl->Z() << endl;
	}
	else {
		// when no steps are used, no correction makes sense
		if (!tr->PropagateTo(rc, 0.0, 0.0)) return distance;
	}
	
	// get global propagation of track at end of propagation
	tr->GetXYZ(pos);
	
	// check angle cut
	TVector3 vc(cl->X(), cl->Y(), cl->Z());
	TVector3 vt(x, y, z);
	Double_t angle = TMath::Abs(vc.Angle(vt)) * TMath::RadToDeg();
	if (angle > fCutAngle) {
		//cout << "angle" << endl;
		return distance;
	}
		
	// compute differences wr to each coordinate
	x -= cl->X();
	if (x > fCutX) {
		//cout << "cut X" << endl;
		return distance;
	}
	y -= cl->Y();
	if (y > fCutY) {
		//cout << "cut Y" << endl;
		return distance;
	}
	z -= cl->Z();
	if (z > fCutZ) {
		//cout << "cut Z" << endl;
		return distance;
	}
	
	// compute true distance
	distance = TMath::Sqrt(x*x + y*y + z*z);
	//Double_t temp = CheckPairV2(tr, cl);
	//if (temp < distance) return temp; else 
	
	// delete temporary object
	delete tr;
	
	return distance;
}
//
//------------------------------------------------------------------------------
//
Double_t AliEMCALTracker::CheckPairV2
(AliEMCALTrack *tr, AliEMCALMatchCluster *cl)
{
	//
	// Given a track and a cluster,
	// propagates the first to the radius of the second.
	// Then, checks the propagation point against all cuts.
	// If at least a cut is not passed, a valuer equal to 
	// twice the maximum allowed distance is passed (so the value returned
	// will not be taken into account when creating matches)
	//
	
	// TEMP
//	Bool_t isTrue = kFALSE;
//	if (tr->GetSeedLabel() == cl->Label()) {
//		isTrue = kTRUE;
//		cout << "TRUE MATCH!!!" << endl;
//	}
	
	Double_t distance = 2.0 * fMaxDist;
	
	Double_t x0, rho;
	if (fTrackCorrMode == kTrackCorrMMB) {
		Double_t pos1[3], pos2[3], param[6];
		tr->GetXYZ(pos1);
//		pos1[0] = x;
//		pos1[1] = y;
//		pos1[2] = z;
		pos2[0] = cl->X();
		pos2[1] = cl->Y();
		pos2[2] = cl->Z();
		AliKalmanTrack::MeanMaterialBudget(pos1, pos2, param);
		rho = param[0];
		x0 = param[1];
	}
	else if (fTrackCorrMode == kTrackCorrFixed) {
		rho = fRho;
		x0 = fX0;
	}
	else {
		rho = 0.0;
		x0 = 0.0;
	}
	
	// check against cut on difference 'alpha - phi'
	Double_t phi = TMath::ATan2(cl->Y(), cl->X());
	phi = AngleDiff(phi, tr->GetAlpha());
	if (phi < fCutAlphaMin || phi > fCutAlphaMax) return distance;
	
	// get cluster position and put them into a vector
	TVector3 vc(cl->X(), cl->Y(), cl->Z());
	// rotate the vector in order to put all clusters on a plane intersecting 
	// vertically the X axis; the angle depends on the sector
	Double_t clusterRot, clusterPhi = vc.Phi() * TMath::RadToDeg();
	if (clusterPhi < 0.0) clusterPhi += 360.0;
	if (clusterPhi < 100.0) {
		clusterRot = -90.0;
	}
	else if (clusterPhi < 120.0) {
		clusterRot = -110.0;
	}
	else if (clusterPhi < 140.0) {
		clusterRot = -130.0;
	}
	else if (clusterPhi < 160.0) {
		clusterRot = -150.0;
	}
	else if (clusterPhi < 180.0) {
		clusterRot = -170.0;
	}
	else {
		clusterRot = -190.0;
	}
	vc.RotateZ(clusterRot * TMath::DegToRad());
	// generate a track from the ESD track selected
	AliEMCALTrack *track = new AliEMCALTrack(*tr);
	// compute the 'phi' coordinate of the intersection point to 
	// the EMCAL surface
	Double_t x = vc.X();
	Double_t y;
	track->GetYAt(vc.X(), track->GetBz(), y);
	Double_t tmp = x*TMath::Cos(track->GetAlpha()) - y*TMath::Sin(track->GetAlpha());
	y = x*TMath::Sin(track->GetAlpha()) + y*TMath::Cos(track->GetAlpha());
	x = tmp;
	Double_t trackPhi = TMath::ATan2(y, x) * TMath::RadToDeg();
	// compute phi difference
	Double_t dphi = trackPhi - clusterPhi;
	if (TMath::Abs(dphi) > 180.0) {
		dphi = 360.0 - TMath::Abs(dphi);
		if (clusterPhi > trackPhi) dphi = -dphi;
	}
	// propagate track to the X position of rotated cluster
	// and get the vector of X, Y, Z in the local ref. frame of the track
	track->PropagateTo(vc.X(), x0, rho);
	TVector3 vt(track->GetX(), track->GetY(), track->GetZ());
	vt.RotateZ((clusterPhi - trackPhi) * TMath::DegToRad());
	TVector3 vdiff = vt-vc;
		
	// compute differences wr to each coordinate
	if (vdiff.X() > fCutX) return distance;
	if (vdiff.Y() > fCutY) return distance;
	if (vdiff.Z() > fCutZ) return distance;
	
	// compute true distance
	distance = vdiff.Mag();
	return distance;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::CreateMatches()
{
	//
	// Creation of matches between tracks and clusters.
	// For each ESD track collected by ReadESD(), an AliEMCALTrack is made.
	// If it finds a cluster close enough to its propagation to EMCAL,
	// which passes all cuts, its index is stored.
	// If many clusters are found which satisfy the criteria described above, 
	// only the closest one is stored.
	// At this level, it is possible that two tracks share the same cluster.
	//
	
	// if matches collection is already present, it is deleted
	if (fMatches) {
		fMatches->Delete();
		delete fMatches;
	}
	fMatches = new TList;
	
	// initialize counters and indexes
	Int_t count = 0;
	Int_t ic, nClusters = (Int_t)fClusters->GetEntries();
	Int_t it, nTracks = fTracks->GetEntries();
	
	// external loop on clusters, internal loop on tracks
	Double_t dist;
	for (ic = 0; ic < nClusters; ic++) {
		cout << "\rMatching cluster " << ic+1 << " of " << nClusters << flush;
		AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(ic);
		for (it = 0; it < nTracks; it++) {
			AliEMCALTrack *track = (AliEMCALTrack*)fTracks->At(it);
			dist = CheckPair(track, cluster);
			if (dist <= fMaxDist) {
				AliEMCALMatch *candidate = new AliEMCALMatch;
				candidate->SetIndexT(it);
				candidate->SetIndexC(ic);
				candidate->SetDistance(dist);
				fMatches->Add(candidate);
				count++;
			}
		}
	}
	cout << endl;		
	
	/*
	// loop on clusters and tracks
	Int_t icBest;
	Double_t dist, distBest;
	for (it = 0; it < nTracks; it++) {
		AliEMCALTrack *track = (AliEMCALTrack*)fTracks->At(it);
		if (!track) continue;
		icBest = -1;
		distBest = fMaxDist;
		for (ic = 0; ic < nClusters; ic++) {
			AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(ic);
			if (!cluster) continue;
			dist = CheckPair(track, cluster);
			if (dist < distBest) {
				distBest = dist;
				icBest = ic;
			}
		}
		if (icBest >= 0) {
			track->SetMatchedClusterIndex(icBest);
			track->SetMatchedClusterDist(distBest);
			count++;
		}
		else {
			track->SetMatchedClusterIndex(-1);
		}
	}
	*/
	
	return count;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::SolveCompetitions()
{
	//
	// Match selector.
	// The match list is sorted from the best to the worst match, w.r. to the 
	// distance between track prolongation and cluster position.
	// Based on this criterion, starting from the first (best) match, a flag
	// is set to both the involved track and cluster, and all matches containing
	// an already used track or cluster are removed, leaving only the best match
	// for each cluster.
	//
	
	// sort matches with respect to track-cluster distance
	fMatches->Sort(kSortAscending);
	
	// keep track of eliminated matches
	Int_t count = 0;
	
	// initialize flags to check repetitions
	Int_t ic, nClusters = (Int_t)fClusters->GetEntries();
	Int_t it, nTracks = fTracks->GetEntries();
	Bool_t *usedC = new Bool_t[nClusters];
	Bool_t *usedT = new Bool_t[nTracks];
	for (ic = 0; ic < nClusters; ic++) usedC[ic] = kFALSE;
	for (it = 0; it < nTracks; it++) usedT[it] = kFALSE;
	
	// loop on matches
	TListIter iter(fMatches);
	AliEMCALMatch *match = 0;
	while ( (match = (AliEMCALMatch*)iter.Next()) ) {
		ic = match->GetIndexC();
		it = match->GetIndexT();
		if (!usedT[it] && !usedC[ic]) {
			usedT[it] = kTRUE;
			usedC[ic] = kTRUE;
			match->CanBeSaved() = kTRUE;
		}
		else {
			count++;
		}
	}
	
	/*
	Int_t it1, it2, nTracks = (Int_t)fTracks->GetEntries();
	AliEMCALTrack *track1 = 0, *track2 = 0;
	for (it1 = 0; it1 < nTracks; it1++) {
		track1 = (AliEMCALTrack*)fTracks->At(it1);
		if (!track1) continue;
		if (track1->GetMatchedClusterIndex() < 0) continue;
		for (it2 = it1+1; it2 < nTracks; it2++) {
			track2 = (AliEMCALTrack*)fTracks->At(it2);
			if (!track2) continue;
			if (track2->GetMatchedClusterIndex() < 0) continue;
			if (track1->GetMatchedClusterIndex() != track2->GetMatchedClusterIndex()) continue;
			count++;
			if (track1->GetMatchedClusterDist() < track2->GetMatchedClusterDist()) {
				track2->SetMatchedClusterIndex(-1);
			}
			else if (track2->GetMatchedClusterDist() < track1->GetMatchedClusterDist()) {
				track1->SetMatchedClusterIndex(-1);
			}
		}
	}
	*/
	
	return count;
}
//
//------------------------------------------------------------------------------
//
void AliEMCALTracker::UnloadClusters() 
{
	//
	// Free memory from clusters
	//
	
  	Clear("CLUSTERS");
}
//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALMatchCluster::AliEMCALMatchCluster(Int_t index, AliEMCALRecPoint *recPoint)
  : fIndex(index),
    fLabel(recPoint->GetPrimaryIndex()),
    fX(0.),
    fY(0.),
    fZ(0.)
{
	//
	// Translates an AliEMCALRecPoint object into the internal format.
	// Index of passed cluster in its native array must be specified.
	//
	TVector3 clpos;
	recPoint->GetGlobalPosition(clpos);
	
	fX = clpos.X();
	fY = clpos.Y();
	fZ = clpos.Z();
}
//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALMatchCluster::AliEMCALMatchCluster(Int_t index, AliESDCaloCluster *caloCluster)
  : fIndex(index),
    fLabel(caloCluster->GetPrimaryIndex()),
    fX(0.),
    fY(0.),
    fZ(0.)
{
	//
	// Translates an AliESDCaloCluster object into the internal format.
	// Index of passed cluster in its native array must be specified.
	//
	Float_t clpos[3];
	caloCluster->GetGlobalPosition(clpos);
	
	fX = (Double_t)clpos[0];
	fY = (Double_t)clpos[1];
	fZ = (Double_t)clpos[2];
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::AliEMCALMatch::Compare(const TObject *obj) const 
{
	//
	// Tracks compared wrt their distance from matched point
	//
	
	AliEMCALTracker::AliEMCALMatch *that = (AliEMCALTracker::AliEMCALMatch*)obj;
	
	Double_t thisDist = fDistance;
	Double_t thatDist = that->GetDistance();
	
	if (thisDist > thatDist) return 1;
	else if (thisDist < thatDist) return -1;
	return 0;
}

AliEMCALTracker::AliEMCALMatch::AliEMCALMatch() 
  : TObject(),
    fCanBeSaved(kFALSE), 
    fIndexC(0), 
    fIndexT(0), 
    fDistance(0.)
{
  //default constructor

}

AliEMCALTracker::AliEMCALMatch::AliEMCALMatch(const AliEMCALMatch& copy)
  : TObject(),
    fCanBeSaved(copy.fCanBeSaved),
    fIndexC(copy.fIndexC),
    fIndexT(copy.fIndexT),
    fDistance(copy.fDistance)
{
  //copy ctor
}
