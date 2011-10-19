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
// Revised by Rongrong 2010-05-31 (rongrong.ma@cern.ch)
//=========================================================================

#include <Riostream.h>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TString.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include <TGeoMatrix.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALRecPoint.h"
#include "AliRunLoader.h"
#include "AliEMCALTrack.h"
#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALReconstructor.h"
#include "AliEMCALRecParam.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALReconstructor.h"

#include "AliEMCALTracker.h"

ClassImp(AliEMCALTracker)

//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALTracker() 
: AliTracker(),
  fCutPt(0),
  fCutNITS(0),
  fCutNTPC(50),
  fStep(10),
  fTrackCorrMode(kTrackCorrMMB),
  fClusterWindow(50),
  fCutEta(0.025),
  fCutPhi(0.05),
  fTracks(0),
  fClusters(0),
  fGeom(0)
{
  //
  // Default constructor.
  // Initializes all simple data members to default values,
   // and all collections to NULL.
  // Output file name is set to a default value.
  //
  InitParameters();
}
//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALTracker(const AliEMCALTracker& copy) 
  : AliTracker(),
    fCutPt(copy.fCutPt),
    fCutNITS(copy.fCutNITS),
    fCutNTPC(copy.fCutNTPC),
    fStep(copy.fStep),
    fTrackCorrMode(copy.fTrackCorrMode),
    fClusterWindow(copy.fClusterWindow),
    fCutEta(copy.fCutEta),
    fCutPhi(copy.fCutPhi),
    fTracks((TObjArray*)copy.fTracks->Clone()),
    fClusters((TObjArray*)copy.fClusters->Clone()),
    fGeom(copy.fGeom)
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

  fCutPt  = copy.fCutPt;
  fClusterWindow = copy.fClusterWindow;
  fCutEta = copy.fCutEta;
  fCutPhi = copy.fCutPhi;	
  fStep = copy.fStep;
  fTrackCorrMode = copy.fTrackCorrMode;

  fCutNITS = copy.fCutNITS;
  fCutNTPC = copy.fCutNTPC;
  
  fTracks = (TObjArray*)copy.fTracks->Clone();
  fClusters = (TObjArray*)copy.fClusters->Clone();
  fGeom = copy.fGeom;
  
  return (*this);
}
//
//------------------------------------------------------------------------------
//
void AliEMCALTracker::InitParameters()
{
  //
  // Retrieve initialization parameters
  //
	
  // Check if the instance of AliEMCALRecParam exists, 
  const AliEMCALRecParam* recParam = AliEMCALReconstructor::GetRecParam();

  if(!recParam){
    AliFatal("Reconstruction parameters for EMCAL not set!");
  }
  else{
 
    fCutEta  =  recParam->GetMthCutEta();
    fCutPhi  =  recParam->GetMthCutPhi();
    fStep    =  recParam->GetExtrapolateStep();
    fCutPt   =  recParam->GetTrkCutPt();
    fCutNITS =  recParam->GetTrkCutNITS();
    fCutNTPC =  recParam->GetTrkCutNTPC();
  }
	
}

//
//------------------------------------------------------------------------------
//
void AliEMCALTracker::Clear(Option_t* option)
{
	//
	// Clearing method
        // Deletes all objects in arrays and the arrays themselves
	//

	TString opt(option);
	Bool_t clearTracks = opt.Contains("TRACKS");
	Bool_t clearClusters = opt.Contains("CLUSTERS");
	if (opt.Contains("ALL")) {
		clearTracks = kTRUE;
		clearClusters = kTRUE;
	}
	
	//fTracks is a collection of esdTrack
	//When clearing this array, the linked objects should not be deleted
	if (fTracks != 0x0 && clearTracks) {
	   fTracks->Clear();
	   delete fTracks;
           fTracks = 0;
	}
	if (fClusters != 0x0 && clearClusters) {
	   fClusters->Delete();
	   delete fClusters;
	   fClusters = 0;
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

	cTree->SetBranchStatus("*",0); //disable all branches
	cTree->SetBranchStatus("EMCALECARP",1); //Enable only the branch we need

	TBranch *branch = cTree->GetBranch("EMCALECARP");
	if (!branch) {
		AliError("Can't get the branch with the EMCAL clusters");
		return 1;
	}
	
	TClonesArray *clusters = new TClonesArray("AliEMCALRecPoint", 1000);
	branch->SetAddress(&clusters);
	
	//cTree->GetEvent(0);
	branch->GetEntry(0);
	Int_t nClusters = (Int_t)clusters->GetEntries();
	if(fClusters) fClusters->Delete();
	else fClusters = new TObjArray(0);
	for (Int_t i = 0; i < nClusters; i++) {
		AliEMCALRecPoint *cluster = (AliEMCALRecPoint*)clusters->At(i);
		if (!cluster) continue;
		AliEMCALMatchCluster *matchCluster = new AliEMCALMatchCluster(i, cluster);
		fClusters->AddLast(matchCluster);
	}

	branch->SetAddress(0);
        clusters->Delete();
        delete clusters;

	AliInfo(Form("Collected %d RecPoints from Tree", fClusters->GetEntries()));

	return 0;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::LoadClusters(AliESDEvent *esd) 
{
  //
  // Load EMCAL clusters in the form of AliESDCaloClusters,
  // from an AliESD object.
  //
  
  // make sure that tracks/clusters collections are empty
  Clear("CLUSTERS");
  fClusters = new TObjArray(0);
  
  Int_t nClusters = esd->GetNumberOfCaloClusters();       		
  for (Int_t i=0; i<nClusters; i++) 
    {
      AliESDCaloCluster *cluster = esd->GetCaloCluster(i);
      if (!cluster || !cluster->IsEMCAL()) continue ; 
      AliEMCALMatchCluster *matchCluster = new AliEMCALMatchCluster(i, cluster);
      fClusters->AddLast(matchCluster);
    }
  
  AliInfo(Form("Collected %d clusters from ESD", fClusters->GetEntries()));
  return 0;
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::LoadTracks(AliESDEvent *esd)
{
  //
  // Load ESD tracks.
  //
	
  Clear("TRACKS");
  fTracks = new TObjArray(0);
	
  Int_t nTracks = esd->GetNumberOfTracks();
  //Bool_t isKink=kFALSE;
  for (Int_t i = 0; i < nTracks; i++) 
    {
      AliESDtrack *esdTrack = esd->GetTrack(i);
      // set by default the value corresponding to "no match"
      esdTrack->SetEMCALcluster(kUnmatched);
      esdTrack->ResetStatus(AliESDtrack::kEMCALmatch);

      //Select good quaulity tracks
      if(esdTrack->Pt()<fCutPt) continue;
      if(esdTrack->GetNcls(1)<fCutNTPC)continue;

      //Loose geometric cut
      Double_t phi = esdTrack->Phi()*TMath::RadToDeg();
      if(TMath::Abs(esdTrack->Eta())>0.8 || phi <= 20 || phi >= 240 ) continue;

      fTracks->AddLast(esdTrack);
    }

      AliInfo(Form("Collected %d tracks", fTracks->GetEntries()));
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
  
  if (!opt.CompareTo("NONE")) 
    {
      fTrackCorrMode = kTrackCorrNone;
    }
  else if (!opt.CompareTo("MMB")) 
    {
      fTrackCorrMode = kTrackCorrMMB;
    }
  else 
    {
      cerr << "E-AliEMCALTracker::SetTrackCorrectionMode '" << option << "': Unrecognized option" << endl;
    }
}
//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::PropagateBack(AliESDEvent* esd)
{
	//
	// Main operation method.
	// Gets external AliESD containing tracks to be matched.
	// After executing match finding, stores in the same ESD object all infos
	// and releases the object for further reconstruction steps.
	//
        //
        // Note: should always return 0=OK, because otherwise all tracking
        // is aborted for this event
  
	if (!esd) {
		AliError("NULL ESD passed");
		return 1;
	}
	
	// step 1: collect clusters
	Int_t okLoadClusters, nClusters;
	if (!fClusters || (fClusters && fClusters->IsEmpty())) {
		okLoadClusters = LoadClusters(esd);
	}
	nClusters = fClusters->GetEntries();
		
	// step 2: collect ESD tracks
	Int_t nTracks, okLoadTracks;
	okLoadTracks = LoadTracks(esd);
	nTracks = fTracks->GetEntries();
	
	// step 3: for each track, find the closest cluster as matched within residual cuts
	Int_t index=-1;
	for (Int_t it = 0; it < nTracks; it++) 
	  {
	    AliESDtrack *track = (AliESDtrack*)fTracks->At(it);
	    index = FindMatchedCluster(track);
	    if (index>-1) 
	      {
		AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(index);
		track->SetEMCALcluster(cluster->Index());
		track->SetStatus(AliESDtrack::kEMCALmatch);
	      }
	  }

	return 0;
}

//
//------------------------------------------------------------------------------
//
Int_t AliEMCALTracker::FindMatchedCluster(AliESDtrack *track)
{         
  //
  // For each track, extrapolate it to all the clusters
  // Find the closest one as matched if the residuals (dEta, dPhi) satisfy the cuts
  //

  Double_t maxEta=fCutEta;
  Double_t maxPhi=fCutPhi;
  Int_t index = -1;
  
  // If the esdFriend is available, use the TPCOuter point as the starting point of extrapolation
  // Otherwise use the TPCInner point
  AliExternalTrackParam *trkParam = 0;
  const AliESDfriendTrack*  friendTrack = track->GetFriendTrack();
  if(friendTrack && friendTrack->GetTPCOut())
    trkParam = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
  else
    trkParam = const_cast<AliExternalTrackParam*>(track->GetInnerParam());
  if(!trkParam) return index;


  AliExternalTrackParam trkParamTmp(*trkParam);
  Double_t eta, phi;
  if(!ExtrapolateTrackToEMCalSurface(&trkParamTmp, 430., track->GetMass(), fStep, eta, phi)) return index;
  if(TMath::Abs(eta)>0.75 || (phi) < 70*TMath::DegToRad() || (phi) > 190*TMath::DegToRad()) return index;

  //Perform extrapolation
  Double_t trkPos[3];
  trkParamTmp.GetXYZ(trkPos);
  Int_t nclusters = fClusters->GetEntries();
  for(Int_t ic=0; ic<nclusters; ic++)
    {
      AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(ic);
      Float_t clsPos[3] = {cluster->X(),cluster->Y(),cluster->Z()};
      Double_t dR = TMath::Sqrt(TMath::Power(trkPos[0]-clsPos[0],2)+TMath::Power(trkPos[1]-clsPos[1],2)+TMath::Power(trkPos[2]-clsPos[2],2));
      if(dR > fClusterWindow) continue;
      
      AliExternalTrackParam trkParTmp(trkParamTmp);

      Double_t tmpEta, tmpPhi;
      if(!ExtrapolateTrackToPosition(&trkParTmp, clsPos,track->GetMass(), fStep, tmpEta, tmpPhi)) continue;
      if(TMath::Abs(tmpPhi)<TMath::Abs(maxPhi) && TMath::Abs(tmpEta)<TMath::Abs(maxEta))
        {
          maxPhi=tmpPhi;
          maxEta=tmpEta;
          index=ic;
        }
      }
  return index;
}


//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTracker::ExtrapolateTrackToEMCalSurface(AliExternalTrackParam *trkParam, Double_t emcalR, Double_t mass, Double_t step, Double_t &eta, Double_t &phi)
{
  eta = -999, phi = -999;
  if(!trkParam) return kFALSE;
  if(!AliTrackerBase::PropagateTrackToBxByBz(trkParam, emcalR, mass, step, kTRUE, 0.8, -1)) return kFALSE;
  Double_t trkPos[3] = {0.,0.,0.};
  if(!trkParam->GetXYZ(trkPos)) return kFALSE;
  TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
  eta = trkPosVec.Eta();
  phi = trkPosVec.Phi();
  if(phi<0)
    phi += 2*TMath::Pi();

  return kTRUE;
}


//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTracker::ExtrapolateTrackToPosition(AliExternalTrackParam *trkParam, Float_t *clsPos, Double_t mass, Double_t step, Double_t &tmpEta, Double_t &tmpPhi)
{
  //
  //Return the residual by extrapolating a track param to a global position
  //
  tmpEta = -999;
  tmpPhi = -999;
  if(!trkParam) return kFALSE;
  Double_t trkPos[3] = {0.,0.,0.};
  TVector3 vec(clsPos[0],clsPos[1],clsPos[2]);
  Double_t alpha =  ((int)(vec.Phi()*TMath::RadToDeg()/20)+0.5)*20*TMath::DegToRad();
  vec.RotateZ(-alpha); //Rotate the cluster to the local extrapolation coordinate system
  if(!AliTrackerBase::PropagateTrackToBxByBz(trkParam, vec.X(), mass, step,kTRUE, 0.8, -1)) return kFALSE;
  if(!trkParam->GetXYZ(trkPos)) return kFALSE; //Get the extrapolated global position

  TVector3 clsPosVec(clsPos[0],clsPos[1],clsPos[2]);
  TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);

  // track cluster matching
  tmpPhi = clsPosVec.DeltaPhi(trkPosVec);    // tmpPhi is between -pi and pi
  tmpEta = clsPosVec.Eta()-trkPosVec.Eta();

  return kTRUE;
}


//
//------------------------------------------------------------------------------
//
Bool_t AliEMCALTracker::ExtrapolateTrackToCluster(AliExternalTrackParam *trkParam, AliVCluster *cluster, Double_t mass, Double_t step, Double_t &tmpEta, Double_t &tmpPhi)
{
  //
  //Return the residual by extrapolating a track param to a cluster
  //
  tmpEta = -999;
  tmpPhi = -999;
  if(!cluster) return kFALSE;

  Float_t clsPos[3] = {0.,0.,0.};
  cluster->GetPosition(clsPos);

  return ExtrapolateTrackToPosition(trkParam, clsPos, mass, step, tmpEta, tmpPhi);
}


//
//------------------------------------------------------------------------------
//
void AliEMCALTracker::UnloadClusters() 
{
	//
	// Free memory from all arrays
        // This method is called after the local tracking step
        // so we can safely delete everything 
	//
	
  	Clear();
}

//
//------------------------------------------------------------------------------
//
AliEMCALTracker::AliEMCALMatchCluster::AliEMCALMatchCluster(Int_t index, AliEMCALRecPoint *recPoint)
  : fIndex(index),
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
    fX(0.),
    fY(0.),
    fZ(0.)
{
	//
	// Translates an AliESDCaloCluster object into the internal format.
	// Index of passed cluster in its native array must be specified.
	//
	Float_t clpos[3]= {0., 0., 0.};
	caloCluster->GetPosition(clpos);
	
	fX = (Double_t)clpos[0];
	fY = (Double_t)clpos[1];
	fZ = (Double_t)clpos[2];
}
