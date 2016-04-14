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
#include "AliEMCALRecoUtils.h"

#include "AliEMCALTracker.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTracker) ;
/// \endcond

///
/// Default constructor.
/// Initializes all simple data members to default values,
/// and all collections to NULL.
/// Output file name is set to a default value.
///
//------------------------------------------------------------------------------
AliEMCALTracker::AliEMCALTracker() : 
  AliTracker(),
  fCutPt(0),
  fCutNITS(0),
  fCutNTPC(50),
  fStep(20),
  fTrackCorrMode(kTrackCorrMMB),
  fEMCalSurfaceDistance(440),
  fClusterWindow(50),
  fCutEta(0.025),
  fCutPhi(0.05),
  fITSTrackSA(kFALSE),
  fTrackInITS(kFALSE),
  fTracks(0),
  fClusters(0),
  fGeom(0)
{
  InitParameters();
}

///
/// Copy constructor.
/// Besides copying all parameters, duplicates all collections.
///
//------------------------------------------------------------------------------
AliEMCALTracker::AliEMCALTracker(const AliEMCALTracker& copy) : 
  AliTracker(),
  fCutPt(copy.fCutPt),
  fCutNITS(copy.fCutNITS),
  fCutNTPC(copy.fCutNTPC),
  fStep(copy.fStep),
  fTrackCorrMode(copy.fTrackCorrMode),
  fEMCalSurfaceDistance(copy.fEMCalSurfaceDistance),
  fClusterWindow(copy.fClusterWindow),
  fCutEta(copy.fCutEta),
  fCutPhi(copy.fCutPhi),
  fITSTrackSA(copy.fITSTrackSA), 
  fTrackInITS(copy.fTrackInITS),  
  fTracks((TObjArray*)copy.fTracks->Clone()),
  fClusters((TObjArray*)copy.fClusters->Clone()),
  fGeom(copy.fGeom)
{
}

///
/// Assignment operator; use copy ctor.
///
//------------------------------------------------------------------------------
AliEMCALTracker& AliEMCALTracker::operator=(const AliEMCALTracker& source)
{ 
  if (&source == this) return *this;

  new (this) AliEMCALTracker(source);
  return *this;
}

///
/// Retrieve initialization parameters.
///
//------------------------------------------------------------------------------
void AliEMCALTracker::InitParameters()
{
  // Check if the instance of AliEMCALRecParam exists, 
  const AliEMCALRecParam* recParam = AliEMCALReconstructor::GetRecParam();

  if (!recParam) 
  {
    AliFatal("Reconstruction parameters for EMCAL not set!");
    return; // not needed, but for coverity.
  }

  fCutEta     = recParam->GetMthCutEta();
  fCutPhi     = recParam->GetMthCutPhi();
  fStep       = recParam->GetExtrapolateStep();
  fCutPt      = recParam->GetTrkCutPt();   
  fCutNITS    = recParam->GetTrkCutNITS(); // Not used?
  fCutNTPC    = recParam->GetTrkCutNTPC();
  fTrackInITS = recParam->GetTrkInITS();
}

///
/// Clearing method.
/// Deletes all objects in arrays and the arrays themselves.
///
//------------------------------------------------------------------------------
void AliEMCALTracker::Clear(Option_t* option)
{
  TString opt(option);
  Bool_t clearTracks = opt.Contains("TRACKS");
  Bool_t clearClusters = opt.Contains("CLUSTERS");
  if (opt.Contains("ALL")) 
  {
    clearTracks = kTRUE;
    clearClusters = kTRUE;
  }
	
  // fTracks is a collection of esdTrack
  // When clearing this array, the linked objects should not be deleted
  if (fTracks != 0x0 && clearTracks) 
  {
    fTracks->Clear();
    delete fTracks;
    fTracks = 0;
  }

  if (fClusters != 0x0 && clearClusters) 
  {
    fClusters->Delete();
    delete fClusters;
    fClusters = 0;
  }
}

///
/// Load EMCAL clusters in the form of AliEMCALRecPoint,
/// from simulation temporary files.
/// (When included in reconstruction chain, this method is used automatically).
///
//------------------------------------------------------------------------------
Int_t AliEMCALTracker::LoadClusters(TTree *cTree) 
{	
  Clear("CLUSTERS");
  
  cTree->SetBranchStatus("*",0); //disable all branches
  cTree->SetBranchStatus("EMCALECARP",1); //Enable only the branch we need

  TBranch *branch = cTree->GetBranch("EMCALECARP");
  if (!branch) 
  {
    AliError("Can't get the branch with the EMCAL clusters");
    return 1;
  }
	
  TClonesArray *clusters = new TClonesArray("AliEMCALRecPoint", 1000);
  branch->SetAddress(&clusters);
	
  //cTree->GetEvent(0);
  branch->GetEntry(0);
  
  Int_t nClusters = (Int_t)clusters->GetEntries();
  if (fClusters) fClusters->Delete();
  else           fClusters = new TObjArray(0);
  
  for (Int_t i = 0; i < nClusters; i++) 
  {
    AliEMCALRecPoint *cluster = (AliEMCALRecPoint*)clusters->At(i);
    if (!cluster) continue;
    
    AliEMCALMatchCluster *matchCluster = new AliEMCALMatchCluster(i, cluster);
    fClusters->AddLast(matchCluster);
  }

  branch->SetAddress(0);
  clusters->Delete();
  delete clusters;

  AliDebug(1,Form("Collected %d RecPoints from Tree", fClusters->GetEntries()));
  
  return 0;
}

///
/// Load EMCAL clusters in the form of AliESDCaloClusters,
/// from an AliESD object.
//
//------------------------------------------------------------------------------
Int_t AliEMCALTracker::LoadClusters(AliESDEvent *esd) 
{ 
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
  
  AliDebug(1,Form("Collected %d clusters from ESD", fClusters->GetEntries()));
  return 0;
}

///
/// Load ESD tracks.
/// Filter them depending on loose track and acceptance cuts.
///
/// If TPC is not present, try matching with ITS tracks only.
/// Set the flag fITSTrackSA.
///
//------------------------------------------------------------------------------
Int_t AliEMCALTracker::LoadTracks(AliESDEvent *esd)
{  
  // In case of runs without TPC but ITS on, 
  // do matching with different cuts
  UInt_t mask1 = esd->GetESDRun()->GetDetectorsInDAQ();
  UInt_t mask2 = esd->GetESDRun()->GetDetectorsInReco();
  Bool_t desc1 = (mask1 >> 3) & 0x1;
  Bool_t desc2 = (mask2 >> 3) & 0x1;
  
  if (desc1==0 || desc2==0) 
  {
    //     AliError(Form("TPC not in DAQ/RECO: %u (%u)/%u (%u)",
    //                   mask1, esd->GetESDRun()->GetDetectorsInReco(),
    //                   mask2, esd->GetESDRun()->GetDetectorsInDAQ()));
    fITSTrackSA = kTRUE;
  }  
    
  Clear("TRACKS");
  fTracks = new TObjArray(0);
  
  Int_t nTracks = esd->GetNumberOfTracks();
  //printf("N tracks %d\n",nTracks);
  
  for (Int_t i = 0; i < nTracks; i++) 
  {
    AliESDtrack *esdTrack = esd->GetTrack(i);
  
    // Set by default the value corresponding to "no match"
    esdTrack->SetEMCALcluster(kUnmatched);
    esdTrack->ResetStatus(AliESDtrack::kEMCALmatch);
    
    //
    // Select tracks depending on status bit 
    // and number of clusters.
    //
    if ( esdTrack->Pt() < fCutPt ) continue;

    // When there are TPC tracks, select them depeding on the fit bit
    // Take TPCout and if requested ITSout (default)
    if ( !fITSTrackSA )
    {
      Int_t nClustersITS = esdTrack->GetITSclusters(0);
      Int_t nClustersTPC = esdTrack->GetTPCclusters(0);

      ULong_t status = esdTrack->GetStatus();

      //Bool_t  tpcIn  = (status&AliESDtrack::kTPCin )==AliESDtrack::kTPCin ;
      //Bool_t  itsIn  = (status&AliESDtrack::kITSin )==AliESDtrack::kITSin ;
      Bool_t  tpcOut = (status&AliESDtrack::kTPCout)==AliESDtrack::kTPCout;
      Bool_t  itsOut = (status&AliESDtrack::kITSout)==AliESDtrack::kITSout;

      // Check ITSout bit, if requested
      if ( fTrackInITS )
      {
	 // if(nClustersITS > 0) // fCutNITS?
         //    printf("ITS clusters %d, Cut %2.0f, kITSin? %d, kITSout? %d - select %d\n", 
	 //           nClustersITS,fCutNITS,itsIn,itsOut,fTrackInITS);

	if ( !itsOut || nClustersITS <= 0 ) continue; // fCutNITS
      }

      // if(nClustersTPC > 0) 
      //  	 printf("TPC clusters %d, Cut %2.0f, kTPCin? %d; kTPCout? %d\n", 
      //  	        nClustersTPC,fCutNTPC,tpcIn,tpcOut);

      if ( !tpcOut || nClustersTPC < fCutNTPC ) continue;        

    } // Do not match with ITS SA track

    //printf("\t accepted!\n");
    
    //
    // Loose geometric cut
    //
    if ( TMath::Abs(esdTrack->Eta()) > 0.9 ) continue;

    // Save some time and memory in case of no DCal present
    
    if(fGeom->GetNumberOfSuperModules() < 13)
    {
      Double_t phi = esdTrack->Phi()*TMath::RadToDeg();
      if ( phi <= 10 || phi >= 250 ) continue;
    }
        
    fTracks->AddLast(esdTrack);
  }
  
  AliDebug(1,Form("Collected %d tracks", fTracks->GetEntries()));
  return 0;
}

///
/// Set track correction mode
/// gets the choice in string format and converts into 
/// internal enum.
///
//------------------------------------------------------------------------------
void AliEMCALTracker::SetTrackCorrectionMode(Option_t *option)
{
  TString opt(option);
  opt.ToUpper();
  
  if      ( !opt.CompareTo("NONE") ) 
    fTrackCorrMode = kTrackCorrNone;
  else if ( !opt.CompareTo("MMB") ) 
    fTrackCorrMode = kTrackCorrMMB;
  else 
    AliError("Unrecognized option");
}

///
/// Main operation method.
/// Gets external AliESD containing tracks to be matched.
/// After executing match finding, stores in the same ESD object all infos
/// and releases the object for further reconstruction steps.
///
/// Note: should always return 0=OK, because otherwise all tracking
/// is aborted for this event.
//------------------------------------------------------------------------------
Int_t AliEMCALTracker::PropagateBack(AliESDEvent* esd)
{ 
  if (!esd)
  {
    AliError("NULL ESD passed");
    return 1;
  }
	
  if( !fGeom ) fGeom = AliEMCALGeometry::GetInstance();
  
  // step 1: collect clusters
  Int_t okLoadClusters = 0;  
  if (!fClusters || (fClusters && fClusters->IsEmpty()))
    okLoadClusters = LoadClusters(esd);
  
  Int_t nClusters = fClusters->GetEntries();
		
  // step 2: collect ESD tracks
  Int_t nTracks, okLoadTracks;
  okLoadTracks = LoadTracks(esd);
  nTracks = fTracks->GetEntries();
  
  AliDebug(5,Form("Propagate back %d tracks ok %d, for %d clusters ok %d",
                  nTracks,okLoadTracks,nClusters,okLoadClusters));
  
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

///
/// For each track, extrapolate it to all the clusters.
/// Find the closest one as matched if the residuals (dEta, dPhi) satisfy the cuts.
///
//------------------------------------------------------------------------------
Int_t AliEMCALTracker::FindMatchedCluster(AliESDtrack *track)
{         
  Float_t maxEta=fCutEta;
  Float_t maxPhi=fCutPhi;
  Int_t index = -1;
  
  // If the esdFriend is available, use the TPCOuter point as the starting point of extrapolation
  // Otherwise use the TPCInner point
  AliExternalTrackParam *trkParam = 0;
  
  if (!fITSTrackSA)
  {
    const AliESDfriendTrack*  friendTrack = track->GetFriendTrack();
  
    if (friendTrack && friendTrack->GetTPCOut())
      trkParam = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
    else if (track->GetInnerParam())
      trkParam = const_cast<AliExternalTrackParam*>(track->GetInnerParam());
  }
  else
    trkParam = new AliExternalTrackParam(*track);
  
  if (!trkParam) return index;
  
  AliExternalTrackParam trkParamTmp(*trkParam);
  Float_t eta, phi, pt;
  if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&trkParamTmp, fEMCalSurfaceDistance, track->GetMass(kTRUE), fStep, eta, phi, pt))
  {
    if (fITSTrackSA) delete trkParam;
    return index;
  }
  
  track->SetTrackPhiEtaPtOnEMCal(phi,eta,pt);

  if ( TMath::Abs(eta) > 0.75 ) 
  {
    if (fITSTrackSA) delete trkParam;
    return index;
  }
  
  // Save some time and memory in case of no DCal present
  if ( fGeom->GetNumberOfSuperModules() < 13 && 
      ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())) 
  {
    if (fITSTrackSA) delete trkParam;
    return index;
  }

  //Perform extrapolation
  Double_t trkPos[3];
  trkParamTmp.GetXYZ(trkPos);
  Int_t nclusters = fClusters->GetEntries();
  for (Int_t ic=0; ic<nclusters; ic++)
  {
    AliEMCALMatchCluster *cluster = (AliEMCALMatchCluster*)fClusters->At(ic);
    
    Float_t clsPos[3] = {static_cast<Float_t>(cluster->X()),
                         static_cast<Float_t>(cluster->Y()),
                         static_cast<Float_t>(cluster->Z())};
    
    Double_t dR = TMath::Sqrt(TMath::Power(trkPos[0]-clsPos[0],2)+TMath::Power(trkPos[1]-clsPos[1],2)+TMath::Power(trkPos[2]-clsPos[2],2));
    //printf("\n dR=%f,wind=%f\n",dR,fClusterWindow); //MARCEL
    
    if (dR > fClusterWindow) continue;
      
    AliExternalTrackParam trkParTmp(trkParamTmp);

    Float_t tmpEta, tmpPhi;
    if (!AliEMCALRecoUtils::ExtrapolateTrackToPosition(&trkParTmp, clsPos,track->GetMass(kTRUE), 5, tmpEta, tmpPhi)) continue;
    
    if (TMath::Abs(tmpPhi)<TMath::Abs(maxPhi) && TMath::Abs(tmpEta)<TMath::Abs(maxEta))
    {
      maxPhi=tmpPhi;
      maxEta=tmpEta;
      index=ic;
    }
  }

  if (fITSTrackSA) delete trkParam;
  
  return index;
}

///
/// Free memory from all arrays.
/// This method is called after the local tracking step.
/// so we can safely delete everything.
///
//------------------------------------------------------------------------------
void AliEMCALTracker::UnloadClusters() 
{	
  Clear();
}

///
/// Translates an AliEMCALRecPoint object into the internal format.
/// Index of passed cluster in its native array must be specified.
///
//------------------------------------------------------------------------------
AliEMCALTracker::AliEMCALMatchCluster::AliEMCALMatchCluster(Int_t index, AliEMCALRecPoint *recPoint) : 
  fIndex(index),
  fX(0.),
  fY(0.),
  fZ(0.)
{
  TVector3 clpos;
  recPoint->GetGlobalPosition(clpos);
  
  fX = clpos.X();
  fY = clpos.Y();
  fZ = clpos.Z();
}

///
/// Translates an AliESDCaloCluster object into the internal format.
/// Index of passed cluster in its native array must be specified.
///
//------------------------------------------------------------------------------
AliEMCALTracker::AliEMCALMatchCluster::AliEMCALMatchCluster(Int_t index, AliESDCaloCluster *caloCluster) : 
  fIndex(index),
  fX(0.),
  fY(0.),
  fZ(0.)
{
  Float_t clpos[3]= {0., 0., 0.};
  caloCluster->GetPosition(clpos);
	
  fX = (Double_t)clpos[0];
  fY = (Double_t)clpos[1];
  fZ = (Double_t)clpos[2];
}
