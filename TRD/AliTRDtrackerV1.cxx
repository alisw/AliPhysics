
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

// #include <Riostream.h>
// #include <stdio.h>
// #include <string.h>

#include <TBranch.h>
#include <TDirectory.h>
#include <TLinearFitter.h>
#include <TTree.h>  
#include <TClonesArray.h>
#include <TTreeStream.h>

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliGeomManager.h"
#include "AliRieman.h"
#include "AliTrackPointArray.h"

#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcalibDB.h"
#include "AliTRDReconstructor.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliTRDrecoParam.h"

#include "AliTRDcluster.h" 
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDtrackerDebug.h"
#include "AliTRDtrackingChamber.h"
#include "AliTRDchamberTimeBin.h"



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
Int_t AliTRDtrackerV1::fgNTimeBins = 0;
TTreeSRedirector *AliTRDtrackerV1::fgDebugStreamer = 0x0;
AliRieman* AliTRDtrackerV1::fgRieman = 0x0;
TLinearFitter* AliTRDtrackerV1::fgTiltedRieman = 0x0;
TLinearFitter* AliTRDtrackerV1::fgTiltedRiemanConstrained = 0x0;

//____________________________________________________________________
AliTRDtrackerV1::AliTRDtrackerV1(AliTRDReconstructor *rec) 
  :AliTracker()
  ,fReconstructor(0x0)
  ,fGeom(new AliTRDgeometry())
  ,fClusters(0x0)
  ,fTracklets(0x0)
  ,fTracks(0x0)
  ,fSieveSeeding(0)
{
  //
  // Default constructor.
  // 
  AliTRDcalibDB *trd = 0x0;
  if (!(trd = AliTRDcalibDB::Instance())) {
    AliFatal("Could not get calibration object");
  }

  if(!fgNTimeBins) fgNTimeBins = trd->GetNumberOfTimeBins();

  for (Int_t isector = 0; isector < AliTRDgeometry::kNsector; isector++) new(&fTrSec[isector]) AliTRDtrackingSector(fGeom, isector);
  
  for(Int_t isl =0; isl<kNSeedPlanes; isl++) fSeedTB[isl] = 0x0;

  // Initialize debug stream
  if(rec) SetReconstructor(rec);
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
  for(Int_t isl =0; isl<kNSeedPlanes; isl++) if(fSeedTB[isl]) delete fSeedTB[isl];
  if(fTracks) {fTracks->Delete(); delete fTracks;}
  if(fTracklets) {fTracklets->Delete(); delete fTracklets;}
  if(fClusters) {
    fClusters->Delete(); delete fClusters;
  }
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

  if(!fReconstructor->GetRecoParam() ){
    AliError("Reconstruction configuration not initialized. Call first AliTRDReconstructor::SetRecoParam().");
    return 0;
  }
  
  //AliInfo("Start Track Finder ...");
  Int_t ntracks = 0;
  for(int ism=0; ism<AliTRDgeometry::kNsector; ism++){
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
  
  // reset position of the point before using it
  p.SetXYZ(0., 0., 0.);
  AliTRDseedV1 *tracklet = GetTracklet(index); 
  if (!tracklet) return kFALSE;

  // get detector for this tracklet
  Int_t  idet     = tracklet->GetDetector();
    
  Double_t local[3];
  local[0] = tracklet->GetX0(); 
  local[1] = tracklet->GetYfit(0);
  local[2] = tracklet->GetZfit(0);
  Double_t global[3];
  fGeom->RotateBack(idet, local, global);
  p.SetXYZ(global[0],global[1],global[2]);
  
  
  // setting volume id
  AliGeomManager::ELayerID iLayer = AliGeomManager::kTRD1;
  switch (fGeom->GetLayer(idet)) {
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
  Int_t    modId = fGeom->GetSector(idet) * fGeom->Nstack() + fGeom->GetStack(idet);
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
  if(!fgRieman) fgRieman = new AliRieman(AliTRDtrackingChamber::kNTimeBins * AliTRDgeometry::kNlayer);
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

  // Calibration monitor
  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if (!calibra) AliInfo("Could not get Calibra instance\n");
  
  Int_t   found    = 0;     // number of tracks found
  Float_t foundMin = 20.0;
  
  Float_t *quality = 0x0;
  Int_t   *index   = 0x0;
  Int_t    nSeed   = event->GetNumberOfTracks();
  if(nSeed){  
    quality = new Float_t[nSeed];
    index   = new Int_t[nSeed];
    for (Int_t iSeed = 0; iSeed < nSeed; iSeed++) {
      AliESDtrack *seed = event->GetTrack(iSeed);
      Double_t covariance[15];
      seed->GetExternalCovariance(covariance);
      quality[iSeed] = covariance[0] + covariance[2];
    }
    // Sort tracks according to covariance of local Y and Z
    TMath::Sort(nSeed,quality,index,kFALSE);
  }
  
  // Backpropagate all seeds
  Int_t   expectedClr;
  AliTRDtrackV1 track;
  for (Int_t iSeed = 0; iSeed < nSeed; iSeed++) {
  
    // Get the seeds in sorted sequence
    AliESDtrack *seed = event->GetTrack(index[iSeed]);
  
    // Check the seed status
    ULong_t status = seed->GetStatus();
    if ((status & AliESDtrack::kTPCout) == 0) continue;
    if ((status & AliESDtrack::kTRDout) != 0) continue;
  
    // Do the back prolongation
    new(&track) AliTRDtrackV1(*seed);
    track.SetReconstructor(fReconstructor);

    //Int_t   lbl         = seed->GetLabel();
    //track.SetSeedLabel(lbl);

    // Make backup and mark entrance in the TRD
    seed->UpdateTrackParams(&track, AliESDtrack::kTRDin);
    seed->UpdateTrackParams(&track, AliESDtrack::kTRDbackup);
    Float_t p4          = track.GetC();
    expectedClr = FollowBackProlongation(track);

    if (expectedClr<0) continue; // Back prolongation failed

    if(expectedClr){
      found++;  
      // computes PID for track
      track.CookPID();
      // update calibration references using this track
      if(calibra->GetHisto2d()) calibra->UpdateHistogramsV1(&track);
      // save calibration object
      if (fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 0 /*&& quality TODO*/){ 
        AliTRDtrackV1 *calibTrack = new AliTRDtrackV1(track);
        calibTrack->SetOwner();
        seed->AddCalibObject(calibTrack);
      }
      //update ESD track
      if ((track.GetNumberOfClusters() > 15) && (track.GetNumberOfClusters() > 0.5*expectedClr)) {
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDout);
        track.UpdateESDtrack(seed);
      }
    }

    if ((TMath::Abs(track.GetC() - p4) / TMath::Abs(p4) < 0.2) ||(track.Pt() > 0.8)) {
      //
      // Make backup for back propagation
      //
      Int_t foundClr = track.GetNumberOfClusters();
      if (foundClr >= foundMin) {
        //AliInfo(Form("Making backup track ncls [%d]...", foundClr));
        //track.CookdEdx();
        //track.CookdEdxTimBin(seed->GetID());
        track.CookLabel(1. - fgkLabelFraction);
        if(track.GetBackupTrack()) UseClusters(track.GetBackupTrack());

        // Sign only gold tracks
        if (track.GetChi2() / track.GetNumberOfClusters() < 4) {
          if ((seed->GetKinkIndex(0)      ==   0) && (track.Pt() <  1.5)){
            //UseClusters(&track);
          }
        }
        Bool_t isGold = kFALSE;
  
        // Full gold track
        if (track.GetChi2() / track.GetNumberOfClusters() < 5) {
          if (track.GetBackupTrack()) seed->UpdateTrackParams(track.GetBackupTrack(),AliESDtrack::kTRDbackup);

          isGold = kTRUE;
        }
  
        // Almost gold track
        if ((!isGold)  && (track.GetNCross() == 0) &&	(track.GetChi2() / track.GetNumberOfClusters()  < 7)) {
          //seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
          if (track.GetBackupTrack()) seed->UpdateTrackParams(track.GetBackupTrack(),AliESDtrack::kTRDbackup);
  
          isGold = kTRUE;
        }
        
        if ((!isGold) && (track.GetBackupTrack())) {
          if ((track.GetBackupTrack()->GetNumberOfClusters() > foundMin) && ((track.GetBackupTrack()->GetChi2()/(track.GetBackupTrack()->GetNumberOfClusters()+1)) < 7)) {
            seed->UpdateTrackParams(track.GetBackupTrack(),AliESDtrack::kTRDbackup);
            isGold = kTRUE;
          }
        }
  
        //if ((track->StatusForTOF() > 0) && (track->GetNCross() == 0) && (Float_t(track->GetNumberOfClusters()) / Float_t(track->GetNExpected())  > 0.4)) {
        //seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
        //}
      }
    }
    
    // Propagation to the TOF (I.Belikov)
    if (track.IsStopped() == kFALSE) {
      Double_t xtof  = 371.0;
      Double_t xTOF0 = 370.0;
    
      Double_t c2    = track.GetSnp() + track.GetC() * (xtof - track.GetX());
      if (TMath::Abs(c2) >= 0.99) continue;
      
      if (!PropagateToX(track, xTOF0, fgkMaxStep)) continue;
  
      // Energy losses taken to the account - check one more time
      c2 = track.GetSnp() + track.GetC() * (xtof - track.GetX());
      if (TMath::Abs(c2) >= 0.99) continue;
      
      //if (!PropagateToX(*track,xTOF0,fgkMaxStep)) {
      //	fHBackfit->Fill(7);
      //delete track;
      //	continue;
      //}
  
      Double_t ymax = xtof * TMath::Tan(0.5 * AliTRDgeometry::GetAlpha());
      Double_t y;
      track.GetYAt(xtof,GetBz(),y);
      if (y >  ymax) {
        if (!track.Rotate( AliTRDgeometry::GetAlpha())) continue;	
      }else if (y < -ymax) {
        if (!track.Rotate(-AliTRDgeometry::GetAlpha())) continue;
      }
          
      if (track.PropagateTo(xtof)) {
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDout);
        track.UpdateESDtrack(seed);
      }
    } else {			
      if ((track.GetNumberOfClusters() > 15) && (track.GetNumberOfClusters() > 0.5*expectedClr)) {
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDout);
  
        track.UpdateESDtrack(seed);
      }
    }
  
    seed->SetTRDQuality(track.StatusForTOF());
    seed->SetTRDBudget(track.GetBudget(0));
  }
  if(index) delete [] index;
  if(quality) delete [] quality;
  

  AliInfo(Form("Number of seeds: %d", nSeed));
  AliInfo(Form("Number of back propagated TRD tracks: %d", found));
      
  // run stand alone tracking
  if (fReconstructor->IsSeeding()) Clusters2Tracks(event);
  
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
      continue;
    }

    ULong_t status = seed->GetStatus();
    // reject tracks which failed propagation in the TRD
    if((status & AliESDtrack::kTRDout) == 0) continue;

    // reject tracks which are produced by the TRD stand alone track finder.
    if((status & AliESDtrack::kTRDin)  == 0) continue;
    nseed++; 

    track.ResetCovariance(50.0);

    // do the propagation and processing
    Bool_t kUPDATE = kFALSE;
    Double_t xTPC = 250.0;
    if(FollowProlongation(track)){	
      // Prolongate to TPC
      if (PropagateToX(track, xTPC, fgkMaxStep)) { //  -with update
  seed->UpdateTrackParams(&track, AliESDtrack::kTRDrefit);
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
  
  Int_t    nClustersExpected = 0;
  Int_t lastplane = 5; //GetLastPlane(&t);
  for (Int_t iplane = lastplane; iplane >= 0; iplane--) {
    Int_t   index   = 0;
    AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
    if(!tracklet) continue;
    if(!tracklet->IsOK()) AliWarning("tracklet not OK");
    
    Double_t x  = tracklet->GetX0();
    // reject tracklets which are not considered for inward refit
    if(x > t.GetX()+fgkMaxStep) continue;

    // append tracklet to track
    t.SetTracklet(tracklet, index);
    
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
        
    Double_t length = TMath::Sqrt(
      (xyz0[0]-xyz1[0])*(xyz0[0]-xyz1[0]) +
      (xyz0[1]-xyz1[1])*(xyz0[1]-xyz1[1]) +
      (xyz0[2]-xyz1[2])*(xyz0[2]-xyz1[2])
    );
    if(length>0.){
      // Get material budget
      Double_t param[7];
      if(AliTracker::MeanMaterialBudget(xyz0, xyz1, param)<=0.) break;
      Double_t xrho= param[0]*param[4];
      Double_t xx0 = param[1]; // Get mean propagation parameters
  
      // Propagate and update		
      t.PropagateTo(x, xx0, xrho);
      if (!AdjustSector(&t)) break;
    }
    
    Double_t maxChi2 = t.GetPredictedChi2(tracklet);
    if (maxChi2 < 1e+10 && t.Update(tracklet, maxChi2)){ 
      nClustersExpected += tracklet->GetN();
    }
  }

  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){
    Int_t index;
    for(int iplane=0; iplane<AliTRDgeometry::kNlayer; iplane++){
      AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
      if(!tracklet) continue;
      t.SetTracklet(tracklet, index);
    }

    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    TTreeSRedirector &cstreamer = *fgDebugStreamer;
    cstreamer << "FollowProlongation"
        << "EventNumber="	<< eventNumber
        << "ncl="					<< nClustersExpected
        //<< "track.="			<< &t
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
  
  AliTRDseedV1 tracklet, *ptrTracklet = 0x0;
  // in case of stand alone tracking we store all the pointers to the tracklets in a temporary array
  AliTRDseedV1 *tracklets[kNPlanes];
  memset(tracklets, 0, sizeof(AliTRDseedV1 *) * kNPlanes);
  for(Int_t ip = 0; ip < kNPlanes; ip++){
    tracklets[ip] = t.GetTracklet(ip);
    t.UnsetTracklet(ip);
  } 

  // Loop through the TRD layers
  for (Int_t ilayer = 0; ilayer < AliTRDgeometry::Nlayer(); ilayer++) {
    // BUILD TRACKLET IF NOT ALREADY BUILT
    Double_t x = 0., y, z, alpha;
    ptrTracklet  = tracklets[ilayer];
    if(!ptrTracklet){
      ptrTracklet = new(&tracklet) AliTRDseedV1(ilayer);
      ptrTracklet->SetReconstructor(fReconstructor);
      alpha = t.GetAlpha();
      Int_t sector = Int_t(alpha/AliTRDgeometry::GetAlpha() + (alpha>0. ? 0 : AliTRDgeometry::kNsector));

      if(!fTrSec[sector].GetNChambers()) continue;
      
      if((x = fTrSec[sector].GetX(ilayer)) < 1.) continue;
    
      if (!t.GetProlongation(x, y, z)) return -1/*nClustersExpected*/;
      Int_t stack = fGeom->GetStack(z, ilayer);
      Int_t nCandidates = stack >= 0 ? 1 : 2;
      z -= stack >= 0 ? 0. : 4.; 
      
      for(int icham=0; icham<nCandidates; icham++, z+=8){
        if((stack = fGeom->GetStack(z, ilayer)) < 0) continue;
      
        if(!(chamber = fTrSec[sector].GetChamber(stack, ilayer))) continue;
      
        if(chamber->GetNClusters() < fgNTimeBins*fReconstructor->GetRecoParam() ->GetFindableClusters()) continue;
      
        x = chamber->GetX();
      
        AliTRDpadPlane *pp = fGeom->GetPadPlane(ilayer, stack);
        tracklet.SetTilt(TMath::Tan(TMath::DegToRad()*pp->GetTiltingAngle()));
        tracklet.SetPadLength(pp->GetLengthIPad());
        tracklet.SetDetector(chamber->GetDetector());
        tracklet.SetX0(x);
        if(!tracklet.Init(&t)){
          t.SetStopped(kTRUE);
          return nClustersExpected;
        }
        if(!tracklet.AttachClustersIter(chamber, 1000./*, kTRUE*/)) continue;
        tracklet.Init(&t);
        
        if(tracklet.GetN() < fgNTimeBins*fReconstructor->GetRecoParam() ->GetFindableClusters()) continue;
      
        break;
      }
      //ptrTracklet->UseClusters();
    }
    if(!ptrTracklet->IsOK()){
      if(x < 1.) continue; //temporary
      if(!PropagateToX(t, x-fgkMaxStep, fgkMaxStep)) return -1/*nClustersExpected*/;
      if(!AdjustSector(&t)) return -1/*nClustersExpected*/;
      if(TMath::Abs(t.GetSnp()) > fgkMaxSnp) return -1/*nClustersExpected*/;
      continue;
    }
    
    // Propagate closer to the current chamber if neccessary 
    x -= clength;
    if (x > (fgkMaxStep + t.GetX()) && !PropagateToX(t, x-fgkMaxStep, fgkMaxStep)) return -1/*nClustersExpected*/;
    if (!AdjustSector(&t)) return -1/*nClustersExpected*/;
    if (TMath::Abs(t.GetSnp()) > fgkMaxSnp) return -1/*nClustersExpected*/;
    
    // load tracklet to the tracker and the track
    ptrTracklet = SetTracklet(ptrTracklet);
    t.SetTracklet(ptrTracklet, fTracklets->GetEntriesFast()-1);
  
  
    // Calculate the mean material budget along the path inside the chamber
    //Calculate global entry and exit positions of the track in chamber (only track prolongation)
    Double_t xyz0[3]; // entry point 
    t.GetXYZ(xyz0);
    alpha = t.GetAlpha();
    x = ptrTracklet->GetX0();
    if (!t.GetProlongation(x, y, z)) return -1/*nClustersExpected*/;
    Double_t xyz1[3]; // exit point
    xyz1[0] =  x * TMath::Cos(alpha) - y * TMath::Sin(alpha); 
    xyz1[1] = +x * TMath::Sin(alpha) + y * TMath::Cos(alpha);
    xyz1[2] =  z;
    Double_t param[7];
    if(AliTracker::MeanMaterialBudget(xyz0, xyz1, param)<=0.) return -1;	
    // The mean propagation parameters
    Double_t xrho = param[0]*param[4]; // density*length
    Double_t xx0  = param[1]; // radiation length
    
    // Propagate and update track
    if (!t.PropagateTo(x, xx0, xrho)) return -1/*nClustersExpected*/;
    if (!AdjustSector(&t)) return -1/*nClustersExpected*/;
    Double_t maxChi2 = t.GetPredictedChi2(ptrTracklet);
    if (!t.Update(ptrTracklet, maxChi2)) return -1/*nClustersExpected*/;
    if (maxChi2<1e+10) { 
      nClustersExpected += ptrTracklet->GetN();
      //t.SetTracklet(&tracklet, index);
    }
    // Reset material budget if 2 consecutive gold
    if(ilayer>0 && t.GetTracklet(ilayer-1) && ptrTracklet->GetN() + t.GetTracklet(ilayer-1)->GetN() > 20) t.SetBudget(2, 0.);

    // Make backup of the track until is gold
    // TO DO update quality check of the track.
    // consider comparison with fTimeBinsRange
    Float_t ratio0 = ptrTracklet->GetN() / Float_t(fgNTimeBins);
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
    
  } // end layers loop

  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){
    TTreeSRedirector &cstreamer = *fgDebugStreamer;
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    //AliTRDtrackV1 *debugTrack = new AliTRDtrackV1(t);
    //debugTrack->SetOwner();
    cstreamer << "FollowBackProlongation"
        << "EventNumber="			<< eventNumber
        << "ncl="							<< nClustersExpected
        //<< "track.="					<< debugTrack
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
    chi2[0] += tracklets[ppl[il]].GetChi2Y();
    chi2[1] += tracklets[ppl[il]].GetChi2Z();
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

  TLinearFitter *fitter = GetTiltedRiemanFitterConstraint();
  fitter->StoreData(kTRUE);
  fitter->ClearPoints();
  AliTRDcluster *cl = 0x0;
  
  Float_t x, y, z, w, t, error, tilt;
  Double_t uvt[2];
  Int_t nPoints = 0;
  for(Int_t ilr = 0; ilr < AliTRDgeometry::kNlayer; ilr++){
    if(!tracklets[ilr].IsOK()) continue;
    for(Int_t itb = 0; itb < fgNTimeBins; itb++){
      if(!tracklets[ilr].IsUsable(itb)) continue;
      cl = tracklets[ilr].GetClusters(itb);
      x = cl->GetX();
      y = cl->GetY();
      z = cl->GetZ();
      tilt = tracklets[ilr].GetTilt();
      // Transformation
      t = 1./(x * x + y * y);
      uvt[0] = 2. * x * t;
      uvt[1] = 2. * x * t * tilt ;
      w = 2. * (y + tilt * (z - zVertex)) * t;
      error = 2. * 0.2 * t;
      fitter->AddPoint(uvt, w, error);
      nPoints++;
    }
  }
  fitter->Eval();

  // Calculate curvature
  Double_t a = fitter->GetParameter(0);
  Double_t b = fitter->GetParameter(1);
  Double_t curvature = a/TMath::Sqrt(b*b + 1);

  Float_t chi2track = fitter->GetChisquare()/Double_t(nPoints);
  for(Int_t ip = 0; ip < AliTRDtrackerV1::kNPlanes; ip++)
    tracklets[ip].SetCC(curvature);

/*  if(fReconstructor->GetStreamLevel() >= 5){
    //Linear Model on z-direction
    Double_t xref = CalculateReferenceX(tracklets);		// Relative to the middle of the stack
    Double_t slope = fitter->GetParameter(2);
    Double_t zref = slope * xref;
    Float_t chi2Z = CalculateChi2Z(tracklets, zref, slope, xref);
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    TTreeSRedirector &treeStreamer = *fgDebugStreamer;
    treeStreamer << "FitTiltedRiemanConstraint"
    << "EventNumber=" 		<< eventNumber
    << "CandidateNumber="	<< candidateNumber
    << "Curvature="				<< curvature
    << "Chi2Track="				<< chi2track
    << "Chi2Z="						<< chi2Z
    << "zref="						<< zref
    << "\n";
  }*/
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
  // Output:      - Chi2 values of the track (in Parameter list)
  //
  TLinearFitter *fitter = GetTiltedRiemanFitter();
  fitter->StoreData(kTRUE);
  fitter->ClearPoints();
  AliTRDLeastSquare zfitter;
  AliTRDcluster *cl = 0x0;

  Double_t xref = CalculateReferenceX(tracklets);
  Double_t x, y, z, t, tilt, dx, w, we;
  Double_t uvt[4];
  Int_t nPoints = 0;
  // Containers for Least-square fitter
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    for(Int_t itb = 0; itb < fgNTimeBins; itb++){
      if(!(cl = tracklets[ipl].GetClusters(itb))) continue;
      if (!tracklets[ipl].IsUsable(itb)) continue;
      x = cl->GetX();
      y = cl->GetY();
      z = cl->GetZ();
      tilt = tracklets[ipl].GetTilt();
      dx = x - xref;
      // Transformation
      t = 1./(x*x + y*y);
      uvt[0] = 2. * x * t;
      uvt[1] = t;
      uvt[2] = 2. * tilt * t;
      uvt[3] = 2. * tilt * dx * t;
      w = 2. * (y + tilt*z) * t;
      // error definition changes for the different calls
      we = 2. * t;
      we *= sigError ? tracklets[ipl].GetSigmaY() : 0.2;
      fitter->AddPoint(uvt, w, we);
      zfitter.AddPoint(&x, z, static_cast<Double_t>(TMath::Sqrt(cl->GetSigmaZ2())));
      nPoints++;
    }
  }
  fitter->Eval();
  zfitter.Eval();

  Double_t offset = fitter->GetParameter(3);
  Double_t slope  = fitter->GetParameter(4);

  // Linear fitter  - not possible to make boundaries
  // Do not accept non possible z and dzdx combinations
  Bool_t acceptablez = kTRUE;
  Double_t zref = 0.0;
  for (Int_t iLayer = 0; iLayer < kNPlanes; iLayer++) {
    if(!tracklets[iLayer].IsOK()) continue;
    zref = offset + slope * (tracklets[iLayer].GetX0() - xref);
    if (TMath::Abs(tracklets[iLayer].GetZProb() - zref) > tracklets[iLayer].GetPadLength() * 0.5 + 1.0) 
      acceptablez = kFALSE;
  }
  if (!acceptablez) {
    Double_t dzmf	= zfitter.GetFunctionParameter(1);
    Double_t zmf	= zfitter.GetFunctionValue(&xref);
    fgTiltedRieman->FixParameter(3, zmf);
    fgTiltedRieman->FixParameter(4, dzmf);
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
  if (curvature > 0.0) 
    curvature  =  a / TMath::Sqrt(curvature);

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
  
/*  if(fReconstructor->GetStreamLevel() >=5){
    TTreeSRedirector &cstreamer = *fgDebugStreamer;
    Int_t eventNumber			= AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber	= AliTRDtrackerDebug::GetCandidateNumber();
    Double_t chi2z = CalculateChi2Z(tracklets, offset, slope, xref);
    cstreamer << "FitTiltedRieman0"
        << "EventNumber="			<< eventNumber
        << "CandidateNumber="	<< candidateNumber
        << "xref="						<< xref
        << "Chi2Z="						<< chi2z
        << "\n";
  }*/
  return chi2track;
}


//____________________________________________________________________
Double_t AliTRDtrackerV1::FitLine(const AliTRDtrackV1 *track, AliTRDseedV1 *tracklets, Bool_t err, Int_t np, AliTrackPoint *points)
{
  AliTRDLeastSquare yfitter, zfitter;
  AliTRDcluster *cl = 0x0;

  AliTRDseedV1 work[kNPlanes], *tracklet = 0x0;
  if(!tracklets){
    for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
      if(!(tracklet = track->GetTracklet(ipl))) continue;
      if(!tracklet->IsOK()) continue;
      new(&work[ipl]) AliTRDseedV1(*tracklet);
    }
    tracklets = &work[0];
  }

  Double_t xref = CalculateReferenceX(tracklets);
  Double_t x, y, z, dx, ye, yr, tilt;
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    for(Int_t itb = 0; itb < fgNTimeBins; itb++){
      if(!(cl = tracklets[ipl].GetClusters(itb))) continue;
      if (!tracklets[ipl].IsUsable(itb)) continue;
      x = cl->GetX();
      z = cl->GetZ();
      dx = x - xref;
      zfitter.AddPoint(&dx, z, static_cast<Double_t>(TMath::Sqrt(cl->GetSigmaZ2())));
    }
  }
  zfitter.Eval();
  Double_t z0    = zfitter.GetFunctionParameter(0);
  Double_t dzdx  = zfitter.GetFunctionParameter(1);
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    for(Int_t itb = 0; itb < fgNTimeBins; itb++){
      if(!(cl = tracklets[ipl].GetClusters(itb))) continue;
      if (!tracklets[ipl].IsUsable(itb)) continue;
      x = cl->GetX();
      y = cl->GetY();
      z = cl->GetZ();
      tilt = tracklets[ipl].GetTilt();
      dx = x - xref;
      yr = y + tilt*(z - z0 - dzdx*dx); 
      // error definition changes for the different calls
      ye = tilt*TMath::Sqrt(cl->GetSigmaZ2());
      ye += err ? tracklets[ipl].GetSigmaY() : 0.2;
      yfitter.AddPoint(&dx, yr, ye);
    }
  }
  yfitter.Eval();
  Double_t y0   = yfitter.GetFunctionParameter(0);
  Double_t dydx = yfitter.GetFunctionParameter(1);
  Double_t chi2 = 0.;//yfitter.GetChisquare()/Double_t(nPoints);

  //update track points array
  if(np && points){
    Float_t xyz[3];
    for(int ip=0; ip<np; ip++){
      points[ip].GetXYZ(xyz);
      xyz[1] = y0 + dydx * (xyz[0] - xref);
      xyz[2] = z0 + dzdx * (xyz[0] - xref);
      points[ip].SetXYZ(xyz);
    }
  }
  return chi2;
}


//_________________________________________________________________________
Double_t AliTRDtrackerV1::FitRiemanTilt(const AliTRDtrackV1 *track, AliTRDseedV1 *tracklets, Bool_t sigError, Int_t np, AliTrackPoint *points)
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
  // Output:      - Chi2 values of the track (in Parameter list)
  //
  TLinearFitter *fitter = GetTiltedRiemanFitter();
  fitter->StoreData(kTRUE);
  fitter->ClearPoints();
  AliTRDLeastSquare zfitter;
  AliTRDcluster *cl = 0x0;

  AliTRDseedV1 work[kNPlanes], *tracklet = 0x0;
  if(!tracklets){
    for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
      if(!(tracklet = track->GetTracklet(ipl))) continue;
      if(!tracklet->IsOK()) continue;
      new(&work[ipl]) AliTRDseedV1(*tracklet);
    }
    tracklets = &work[0];
  }

  Double_t xref = CalculateReferenceX(tracklets);
  Double_t x, y, z, t, tilt, dx, w, we;
  Double_t uvt[4];
  Int_t nPoints = 0;
  // Containers for Least-square fitter
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    for(Int_t itb = 0; itb < fgNTimeBins; itb++){
      if(!(cl = tracklets[ipl].GetClusters(itb))) continue;
      if (!tracklets[ipl].IsUsable(itb)) continue;
      x = cl->GetX();
      y = cl->GetY();
      z = cl->GetZ();
      tilt = tracklets[ipl].GetTilt();
      dx = x - xref;
      // Transformation
      t = 1./(x*x + y*y);
      uvt[0] = 2. * x * t;
      uvt[1] = t;
      uvt[2] = 2. * tilt * t;
      uvt[3] = 2. * tilt * dx * t;
      w = 2. * (y + tilt*z) * t;
      // error definition changes for the different calls
      we = 2. * t;
      we *= sigError ? tracklets[ipl].GetSigmaY() : 0.2;
      fitter->AddPoint(uvt, w, we);
      zfitter.AddPoint(&x, z, static_cast<Double_t>(TMath::Sqrt(cl->GetSigmaZ2())));
      nPoints++;
    }
  }
  if(fitter->Eval()) return 1.E10;

  Double_t z0    = fitter->GetParameter(3);
  Double_t dzdx  = fitter->GetParameter(4);


  // Linear fitter  - not possible to make boundaries
  // Do not accept non possible z and dzdx combinations
  Bool_t accept = kTRUE;
  Double_t zref = 0.0;
  for (Int_t iLayer = 0; iLayer < kNPlanes; iLayer++) {
    if(!tracklets[iLayer].IsOK()) continue;
    zref = z0 + dzdx * (tracklets[iLayer].GetX0() - xref);
    if (TMath::Abs(tracklets[iLayer].GetZProb() - zref) > tracklets[iLayer].GetPadLength() * 0.5 + 1.0) 
      accept = kFALSE;
  }
  if (!accept) {
    zfitter.Eval();
    Double_t dzmf	= zfitter.GetFunctionParameter(1);
    Double_t zmf	= zfitter.GetFunctionValue(&xref);
    fitter->FixParameter(3, zmf);
    fitter->FixParameter(4, dzmf);
    fitter->Eval();
    fitter->ReleaseParameter(3);
    fitter->ReleaseParameter(4);
    z0   = fitter->GetParameter(3); // = zmf ?
    dzdx = fitter->GetParameter(4); // = dzmf ?
  }

  // Calculate Curvature
  Double_t a    =  fitter->GetParameter(0);
  Double_t b    =  fitter->GetParameter(1);
  Double_t c    =  fitter->GetParameter(2);
  Double_t y0   = 1. / a;
  Double_t x0   = -b * y0;
  Double_t tmp  = y0*y0 + x0*x0 - c*y0;
  if(tmp<=0.) return 1.E10;
  Double_t R    = TMath::Sqrt(tmp);
  Double_t C    =  1.0 + b*b - c*a;
  if (C > 0.0) C  =  a / TMath::Sqrt(C);

  // Calculate chi2 of the fit 
  Double_t chi2 = fitter->GetChisquare()/Double_t(nPoints);

  // Update the tracklets
  if(!track){
    for(Int_t ip = 0; ip < kNPlanes; ip++) {
      x = tracklets[ip].GetX0();
      tmp = R*R-(x-x0)*(x-x0);  
      if(tmp <= 0.) continue;
      tmp = TMath::Sqrt(tmp);  

      // y:     R^2 = (x - x0)^2 + (y - y0)^2
      //     =>   y = y0 +/- Sqrt(R^2 - (x - x0)^2)
      tracklets[ip].SetYref(0, y0 - (y0>0.?1.:-1)*tmp);
      //     => dy/dx = (x - x0)/Sqrt(R^2 - (x - x0)^2) 
      tracklets[ip].SetYref(1, (x - x0) / tmp);
      tracklets[ip].SetZref(0, z0 + dzdx * (x - xref));
      tracklets[ip].SetZref(1, dzdx);
      tracklets[ip].SetC(C);
      tracklets[ip].SetChi2(chi2);
    }
  }
  //update track points array
  if(np && points){
    Float_t xyz[3];
    for(int ip=0; ip<np; ip++){
      points[ip].GetXYZ(xyz);
      xyz[1] = TMath::Abs(xyz[0] - x0) > R ? 100. : y0 - (y0>0.?1.:-1.)*TMath::Sqrt(R*R-(xyz[0]-x0)*(xyz[0]-x0));
      xyz[2] = z0 + dzdx * (xyz[0] - xref);
      points[ip].SetXYZ(xyz);
    }
  }
  
  return chi2;
}


//____________________________________________________________________
Double_t AliTRDtrackerV1::FitKalman(AliTRDtrackV1 *track, AliTRDseedV1 *tracklets, Bool_t up, Int_t np, AliTrackPoint *points)
{
//   Kalman filter implementation for the TRD.
//   It returns the positions of the fit in the array "points"
// 
//   Author : A.Bercuci@gsi.de

  // printf("Start track @ x[%f]\n", track->GetX());
	
  //prepare marker points along the track
  Int_t ip = np ? 0 : 1;
  while(ip<np){
    if((up?-1:1) * (track->GetX() - points[ip].GetX()) > 0.) break;
    //printf("AliTRDtrackerV1::FitKalman() : Skip track marker x[%d] = %7.3f. Before track start ( %7.3f ).\n", ip, points[ip].GetX(), track->GetX());
    ip++;
  }
  //if(points) printf("First marker point @ x[%d] = %f\n", ip, points[ip].GetX());


  AliTRDseedV1 tracklet, *ptrTracklet = 0x0;

  //Loop through the TRD planes
  for (Int_t jplane = 0; jplane < kNPlanes; jplane++) {
    // GET TRACKLET OR BUILT IT		
    Int_t iplane = up ? jplane : kNPlanes - 1 - jplane;
    if(tracklets){ 
      if(!(ptrTracklet = &tracklets[iplane])) continue;
    }else{
      if(!(ptrTracklet  = track->GetTracklet(iplane))){ 
      /*AliTRDtrackerV1 *tracker = 0x0;
        if(!(tracker = dynamic_cast<AliTRDtrackerV1*>( AliTRDReconstructor::Tracker()))) continue;
        ptrTracklet = new(&tracklet) AliTRDseedV1(iplane);
        if(!tracker->MakeTracklet(ptrTracklet, track)) */
        continue;
      }
    }
    if(!ptrTracklet->IsOK()) continue;

    Double_t x = ptrTracklet->GetX0();

    while(ip < np){
      //don't do anything if next marker is after next update point.
      if((up?-1:1) * (points[ip].GetX() - x) - fgkMaxStep < 0) break;
      if(((up?-1:1) * (points[ip].GetX() - track->GetX()) < 0) && !PropagateToX(*track, points[ip].GetX(), fgkMaxStep)) return -1.;
      
      Double_t xyz[3]; // should also get the covariance
      track->GetXYZ(xyz);
      track->Global2LocalPosition(xyz, track->GetAlpha());
      points[ip].SetXYZ(xyz[0], xyz[1], xyz[2]);
      ip++;
    }
    // printf("plane[%d] tracklet[%p] x[%f]\n", iplane, ptrTracklet, x);

    // Propagate closer to the next update point 
    if(((up?-1:1) * (x - track->GetX()) + fgkMaxStep < 0) && !PropagateToX(*track, x + (up?-1:1)*fgkMaxStep, fgkMaxStep)) return -1.;

    if(!AdjustSector(track)) return -1;
    if(TMath::Abs(track->GetSnp()) > fgkMaxSnp) return -1;
    
    //load tracklet to the tracker and the track
/*    Int_t index;
    if((index = FindTracklet(ptrTracklet)) < 0){
      ptrTracklet = SetTracklet(&tracklet);
      index = fTracklets->GetEntriesFast()-1;
    }
    track->SetTracklet(ptrTracklet, index);*/


    // register tracklet to track with tracklet creation !!
    // PropagateBack : loaded tracklet to the tracker and update index 
    // RefitInward : update index 
    // MakeTrack   : loaded tracklet to the tracker and update index 
    if(!tracklets) track->SetTracklet(ptrTracklet, -1);
    
  
    //Calculate the mean material budget along the path inside the chamber
    Double_t xyz0[3]; track->GetXYZ(xyz0);
    Double_t alpha = track->GetAlpha();
    Double_t xyz1[3], y, z;
    if(!track->GetProlongation(x, y, z)) return -1;
    xyz1[0] =  x * TMath::Cos(alpha) - y * TMath::Sin(alpha); 
    xyz1[1] = +x * TMath::Sin(alpha) + y * TMath::Cos(alpha);
    xyz1[2] =  z;
    if((xyz0[0] - xyz1[9] < 1e-3) && (xyz0[0] - xyz1[9] < 1e-3)) continue; // check wheter we are at the same global x position
    Double_t param[7];
    if(AliTracker::MeanMaterialBudget(xyz0, xyz1, param) <=0.) break;	
    Double_t xrho = param[0]*param[4]; // density*length
    Double_t xx0  = param[1]; // radiation length
    
    //Propagate the track
    track->PropagateTo(x, xx0, xrho);
    if (!AdjustSector(track)) break;
  
    //Update track
    Double_t chi2 = track->GetPredictedChi2(ptrTracklet);
    if(chi2<1e+10) track->Update(ptrTracklet, chi2);
    if(!up) continue;

		//Reset material budget if 2 consecutive gold
		if(iplane>0 && track->GetTracklet(iplane-1) && ptrTracklet->GetN() + track->GetTracklet(iplane-1)->GetN() > 20) track->SetBudget(2, 0.);
	} // end planes loop

  // extrapolation
  while(ip < np){
    if(((up?-1:1) * (points[ip].GetX() - track->GetX()) < 0) && !PropagateToX(*track, points[ip].GetX(), fgkMaxStep)) return -1.;
    
    Double_t xyz[3]; // should also get the covariance
    track->GetXYZ(xyz); 
    track->Global2LocalPosition(xyz, track->GetAlpha());
    points[ip].SetXYZ(xyz[0], xyz[1], xyz[2]);
    ip++;
  }

	return track->GetChi2();
}

//_________________________________________________________________________
Float_t AliTRDtrackerV1::CalculateChi2Z(AliTRDseedV1 *tracklets, Double_t offset, Double_t slope, Double_t xref)
{
  //
  // Calculates the chi2-value of the track in z-Direction including tilting pad correction.
  // A linear dependence on the x-value serves as a model.
  // The parameters are related to the tilted Riemann fit.
  // Parameters: - Array of tracklets (AliTRDseedV1) related to the track candidate
  //             - the offset for the reference x
  //             - the slope
  //             - the reference x position
  // Output:     - The Chi2 value of the track in z-Direction
  //
  Float_t chi2Z = 0, nLayers = 0;
  for (Int_t iLayer = 0; iLayer < AliTRDgeometry::kNlayer; iLayer++) {
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
    if(AliTracker::MeanMaterialBudget(xyz0, xyz1, param)<=0.) return 0;

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
    Float_t nclusters =  fReconstructor->GetRecoParam()->GetNClusters();
    if(fReconstructor->IsHLT()) nclusters /= AliTRDgeometry::kNsector;
    array = new TClonesArray("AliTRDcluster", Int_t(nclusters));
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
      c->SetInChamber();
      new((*fClusters)[ncl++]) AliTRDcluster(*c);
      delete (clusterArray->RemoveAt(iCluster)); 
    }

  }
  delete clusterArray;

  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDtrackerV1::LoadClusters(TTree *cTree)
{
  //
  // Fills clusters into TRD tracking sectors
  //
  
  if(!fReconstructor->IsWritingClusters()){ 
    fClusters = AliTRDReconstructor::GetClusters();
  } else {
    if (ReadClusters(fClusters, cTree)) {
      AliError("Problem with reading the clusters !");
      return 1;
    }
  }
  SetClustersOwner();

  if(!fClusters || !fClusters->GetEntriesFast()){ 
    AliInfo("No TRD clusters");
    return 1;
  }

  //Int_t nin = 
  BuildTrackingContainers();  

  //Int_t ncl  = fClusters->GetEntriesFast();
  //AliInfo(Form("Clusters %d [%6.2f %% in the active volume]", ncl, 100.*float(nin)/ncl));

  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDtrackerV1::LoadClusters(TClonesArray *clusters)
{
  //
  // Fills clusters into TRD tracking sectors
  // Function for use in the HLT
  
  if(!clusters || !clusters->GetEntriesFast()){ 
    AliInfo("No TRD clusters");
    return 1;
  }

  fClusters = clusters;
  SetClustersOwner();

  //Int_t nin = 
  BuildTrackingContainers();  

  //Int_t ncl  = fClusters->GetEntriesFast();
  //AliInfo(Form("Clusters %d [%6.2f %% in the active volume]", ncl, 100.*float(nin)/ncl));

  return 0;
}


//____________________________________________________________________
Int_t AliTRDtrackerV1::BuildTrackingContainers()
{
// Building tracking containers for clusters

  Int_t nin =0, icl = fClusters->GetEntriesFast();
  while (icl--) {
    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(icl);
    if(c->IsInChamber()) nin++;
    Int_t detector       = c->GetDetector();
    Int_t sector         = fGeom->GetSector(detector);
    Int_t stack          = fGeom->GetStack(detector);
    Int_t layer          = fGeom->GetLayer(detector);
    
    fTrSec[sector].GetChamber(stack, layer, kTRUE)->InsertCluster(c, icl);
  }

  const AliTRDCalDet *cal = AliTRDcalibDB::Instance()->GetT0Det();
  for(int isector =0; isector<AliTRDgeometry::kNsector; isector++){ 
    if(!fTrSec[isector].GetNChambers()) continue;
    fTrSec[isector].Init(fReconstructor, cal);
  }

  return nin;
}



//____________________________________________________________________
void AliTRDtrackerV1::UnloadClusters() 
{ 
  //
  // Clears the arrays of clusters and tracks. Resets sectors and timebins 
  //

  if(fTracks) fTracks->Delete(); 
  if(fTracklets) fTracklets->Delete();
  if(fClusters){ 
    if(IsClustersOwner()) fClusters->Delete();
    
    // save clusters array in the reconstructor for further use.
    if(!fReconstructor->IsWritingClusters()){
      AliTRDReconstructor::SetClusters(fClusters);
      SetClustersOwner(kFALSE);
    } else AliTRDReconstructor::SetClusters(0x0);
  }

  for (int i = 0; i < AliTRDgeometry::kNsector; i++) fTrSec[i].Clear();

  // Increment the Event Number
  AliTRDtrackerDebug::SetEventNumber(AliTRDtrackerDebug::GetEventNumber()  + 1);
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
  AliTRDseedV1 *tracklet = (idx==0xffff) ? 0x0 : (AliTRDseedV1*)fTracklets->UncheckedAt(idx);

  return tracklet;
}

//____________________________________________________________________
AliTRDseedV1* AliTRDtrackerV1::SetTracklet(AliTRDseedV1 *tracklet)
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
    fTracklets = new TClonesArray("AliTRDseedV1", AliTRDgeometry::Nsector()*kMaxTracksStack);
    fTracklets->SetOwner(kTRUE);
  }
  Int_t nentries = fTracklets->GetEntriesFast();
  return new ((*fTracklets)[nentries]) AliTRDseedV1(*tracklet);
}

//____________________________________________________________________
AliTRDtrackV1* AliTRDtrackerV1::SetTrack(AliTRDtrackV1 *track)
{
  // Add this track to the list of tracks stored in the tracker
  //
  // Parameters
  //   - track : pointer to the track to be added to the list
  //
  // Output
  //   - the pointer added
  //
  // Detailed description
  // Build the tracks list if it is not yet created (late initialization)
  // and adds the new track to the list.
  //
  if(!fTracks){
    fTracks = new TClonesArray("AliTRDtrackV1", AliTRDgeometry::Nsector()*kMaxTracksStack);
    fTracks->SetOwner(kTRUE);
  }
  Int_t nentries = fTracks->GetEntriesFast();
  return new ((*fTracks)[nentries]) AliTRDtrackV1(*track);
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
  for(int istack = 0; istack<AliTRDgeometry::kNstack; istack++){
    if(!(stack = fTrSec[sector].GetStack(istack))) continue;
    nChambers = 0;
    for(int ilayer=0; ilayer<AliTRDgeometry::kNlayer; ilayer++){
      if(!(chamber = stack[ilayer])) continue;
      if(chamber->GetNClusters() < fgNTimeBins * fReconstructor->GetRecoParam() ->GetFindableClusters()) continue;
      nChambers++;
      //AliInfo(Form("sector %d stack %d layer %d clusters %d", sector, istack, ilayer, chamber->GetNClusters()));
    }
    if(nChambers < 4) continue;
    //AliInfo(Form("Doing stack %d", istack));
    nTracks += Clusters2TracksStack(stack, &esdTrackList);
  }
  //AliInfo(Form("Found %d tracks in SM %d [%d]\n", nTracks, sector, esd->GetNumberOfTracks()));
  
  for(int itrack=0; itrack<nTracks; itrack++)
    esd->AddTrack((AliESDtrack*)esdTrackList[itrack]);

  // Reset Track and Candidate Number
  AliTRDtrackerDebug::SetCandidateNumber(0);
  AliTRDtrackerDebug::SetTrackNumber(0);
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

  const AliTRDCalDet *cal = AliTRDcalibDB::Instance()->GetT0Det();
  AliTRDtrackingChamber *chamber = 0x0;
  AliTRDseedV1 sseed[kMaxTracksStack*6]; // to be initialized
  Int_t pars[4]; // MakeSeeds parameters

  //Double_t alpha = AliTRDgeometry::GetAlpha();
  //Double_t shift = .5 * alpha;
  Int_t configs[kNConfigs];
  
  // Build initial seeding configurations
  Double_t quality = BuildSeedingConfigs(stack, configs);
  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){
    AliInfo(Form("Plane config %d %d %d Quality %f"
    , configs[0], configs[1], configs[2], quality));
  }

  
  // Initialize contors
  Int_t ntracks,      // number of TRD track candidates
    ntracks1,     // number of registered TRD tracks/iter
    ntracks2 = 0; // number of all registered TRD tracks in stack
  fSieveSeeding = 0;

  // Get stack index
  Int_t ic = 0; AliTRDtrackingChamber **cIter = &stack[0];
  while(ic<kNPlanes && !(*cIter)){ic++; cIter++;}
  if(!(*cIter)) return ntracks2;
  Int_t istack = fGeom->GetStack((*cIter)->GetDetector());

  do{
    // Loop over seeding configurations
    ntracks = 0; ntracks1 = 0;
    for (Int_t iconf = 0; iconf<3; iconf++) {
      pars[0] = configs[iconf];
      pars[1] = ntracks;
      pars[2] = istack;
      ntracks = MakeSeeds(stack, &sseed[6*ntracks], pars);
      if(ntracks == kMaxTracksStack) break;
    }
    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1) AliInfo(Form("Candidate TRD tracks %d in iteration %d.", ntracks, fSieveSeeding));
    
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
            
        
  // Sign clusters
  AliTRDcluster *cl = 0x0; Int_t clusterIndex = -1;
  for (Int_t jLayer = 0; jLayer < kNPlanes; jLayer++) {
    Int_t jseed = kNPlanes*trackIndex+jLayer;
    if(!sseed[jseed].IsOK()) continue;
    if(TMath::Abs(sseed[jseed].GetYfit(1) - sseed[jseed].GetYfit(1)) >= .2) continue; // check this condition with Marian
    sseed[jseed].UseClusters();
    if(!cl){
      ic = 0;
      while(!(cl = sseed[jseed].GetClusters(ic))) ic++;
      clusterIndex =  sseed[jseed].GetIndexes(ic);
    }
  }
  if(!cl) continue;

        
  // Build track parameters
  AliTRDseedV1 *lseed =&sseed[trackIndex*6];
/*  Int_t idx = 0;
  while(idx<3 && !lseed->IsOK()) {
    idx++;
    lseed++;
  }*/
  Double_t x = lseed->GetX0();// - 3.5;
  trackParams[0] = x; //NEW AB
  trackParams[1] = lseed->GetYref(0); // lseed->GetYat(x);  
  trackParams[2] = lseed->GetZref(0); // lseed->GetZat(x); 
  trackParams[3] = TMath::Sin(TMath::ATan(lseed->GetYref(1)));
  trackParams[4] = lseed->GetZref(1) / TMath::Sqrt(1. + lseed->GetYref(1) * lseed->GetYref(1));
  trackParams[5] = lseed->GetC();
  Int_t ich = 0; while(!(chamber = stack[ich])) ich++;
  trackParams[6] = fGeom->GetSector(chamber->GetDetector());/* *alpha+shift;	// Supermodule*/

  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){
    AliInfo(Form("Track %d [%d] nlayers %d trackQuality = %e nused %d, yref = %3.3f", itrack, trackIndex, nlayers, fTrackQuality[trackIndex], nused, trackParams[1]));
          
    Int_t nclusters = 0;
    AliTRDseedV1 *dseed[6];

    // Build track label - what happens if measured data ???
    Int_t labels[1000];
    Int_t outlab[1000];
    Int_t nlab = 0;

    Int_t labelsall[1000];
    Int_t nlabelsall = 0;
    Int_t naccepted  = 0;

    for (Int_t iLayer = 0; iLayer < kNPlanes; iLayer++) {
      Int_t jseed = kNPlanes*trackIndex+iLayer;
      dseed[iLayer] = new AliTRDseedV1(sseed[jseed]);
      dseed[iLayer]->SetOwner();
      nclusters += sseed[jseed].GetN2();
      if(!sseed[jseed].IsOK()) continue;
      for(int ilab=0; ilab<2; ilab++){
        if(sseed[jseed].GetLabels(ilab) < 0) continue;
        labels[nlab] = sseed[jseed].GetLabels(ilab);
        nlab++;
      }

      // Cooking label
      for (Int_t itime = 0; itime < fgNTimeBins; itime++) {
        if(!sseed[jseed].IsUsable(itime)) continue;
        naccepted++;
        Int_t tindex = 0, ilab = 0;
        while(ilab<3 && (tindex = sseed[jseed].GetClusters(itime)->GetLabel(ilab)) >= 0){
          labelsall[nlabelsall++] = tindex;
          ilab++;
        }
      }
    }
    Freq(nlab,labels,outlab,kFALSE);
    Int_t   label     = outlab[0];
    Int_t   frequency = outlab[1];
    Freq(nlabelsall,labelsall,outlab,kFALSE);
    Int_t   label1    = outlab[0];
    Int_t   label2    = outlab[2];
    Float_t fakeratio = (naccepted - outlab[1]) / Float_t(naccepted);

    //Int_t eventNrInFile = esd->GetEventNumberInFile();
    //AliInfo(Form("Number of clusters %d.", nclusters));
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t trackNumber = AliTRDtrackerDebug::GetTrackNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    TTreeSRedirector &cstreamer = *fgDebugStreamer;
    cstreamer << "Clusters2TracksStack"
        << "EventNumber="		<< eventNumber
        << "TrackNumber="		<< trackNumber
        << "CandidateNumber="	<< candidateNumber
        << "Iter="				<< fSieveSeeding
        << "Like="				<< fTrackQuality[trackIndex]
        << "S0.="				<< dseed[0]
        << "S1.="				<< dseed[1]
        << "S2.="				<< dseed[2]
        << "S3.="				<< dseed[3]
        << "S4.="				<< dseed[4]
        << "S5.="				<< dseed[5]
        << "p0="				<< trackParams[0]
        << "p1="				<< trackParams[1]
        << "p2="				<< trackParams[2]
        << "p3="				<< trackParams[3]
        << "p4="				<< trackParams[4]
        << "p5="				<< trackParams[5]
        << "p6="				<< trackParams[6]
        << "Label="				<< label
        << "Label1="			<< label1
        << "Label2="			<< label2
        << "FakeRatio="			<< fakeratio
        << "Freq="				<< frequency
        << "Ncl="				<< ncl
        << "NLayers="			<< nlayers
        << "Findable="			<< findable
        << "NUsed="				<< nused
        << "\n";
  }
      
  AliTRDtrackV1 *track = MakeTrack(&sseed[trackIndex*kNPlanes], trackParams);
  if(!track){
    AliWarning("Fail to build a TRD Track.");
    continue;
  }

  //AliInfo("End of MakeTrack()");
  AliESDtrack *esdTrack = new ((*esdTrackList)[ntracks0++]) AliESDtrack();
  esdTrack->UpdateTrackParams(track, AliESDtrack::kTRDout);
  esdTrack->SetLabel(track->GetLabel());
  track->UpdateESDtrack(esdTrack);
  // write ESD-friends if neccessary
  if (fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 0){
    AliTRDtrackV1 *calibTrack = new AliTRDtrackV1(*track);
    calibTrack->SetOwner();
    esdTrack->AddCalibObject(calibTrack);
  }
  ntracks1++;
  AliTRDtrackerDebug::SetTrackNumber(AliTRDtrackerDebug::GetTrackNumber() + 1);
      }

      jSieve++;
    } while(jSieve<5 && candidates); // end track candidates sieve
    if(!ntracks1) break;

    // increment counters
    ntracks2 += ntracks1;

    if(fReconstructor->IsHLT()) break;
    fSieveSeeding++;

    // Rebuild plane configurations and indices taking only unused clusters into account
    quality = BuildSeedingConfigs(stack, configs);
    if(quality < 1.E-7) break; //fReconstructor->GetRecoParam() ->GetPlaneQualityThreshold()) break;
    
    for(Int_t ip = 0; ip < kNPlanes; ip++){ 
      if(!(chamber = stack[ip])) continue;
      chamber->Build(fGeom, cal);//Indices(fSieveSeeding);
    }

    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){ 
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
    chamberQ[iplane] = (chamber = stack[iplane]) ?  chamber->GetQuality() : 0.;
  }

  Double_t tconfig[kNConfigs];
  Int_t planes[4];
  for(int iconf=0; iconf<kNConfigs; iconf++){
    GetSeedingConfig(iconf, planes);
    tconfig[iconf] = fgTopologicQA[iconf];
    for(int iplane=0; iplane<4; iplane++) tconfig[iconf] *= chamberQ[planes[iplane]]; 
  }
  
  TMath::Sort((Int_t)kNConfigs, tconfig, configs, kTRUE);
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
  AliTRDcluster *c[kNSeedPlanes] = {0x0, 0x0, 0x0, 0x0}; // initilize seeding clusters
  AliTRDseedV1 *cseed = &sseed[0]; // initialize tracklets for first track
  Int_t ncl, mcl; // working variable for looping over clusters
  Int_t index[AliTRDchamberTimeBin::kMaxClustersLayer], jndex[AliTRDchamberTimeBin::kMaxClustersLayer];
  // chi2 storage
  // chi2[0] = tracklet chi2 on the Z direction
  // chi2[1] = tracklet chi2 on the R direction
  Double_t chi2[4];

	// Default positions for the anode wire in all 6 Layers in case of a stack with missing clusters
	// Positions taken using cosmic data taken with SM3 after rebuild
  Double_t x_def[kNPlanes] = {300.2, 312.8, 325.4, 338.0, 350.6, 363.2};

  // this should be data member of AliTRDtrack
  Double_t seedQuality[kMaxTracksStack];
  
  // unpack control parameters
  Int_t config  = ipar[0];
  Int_t ntracks = ipar[1];
  Int_t istack  = ipar[2];
  Int_t planes[kNSeedPlanes]; GetSeedingConfig(config, planes);	
  Int_t planesExt[kNPlanes-kNSeedPlanes];         GetExtrapolationConfig(config, planesExt);


  // Init chambers geometry
  Double_t hL[kNPlanes];       // Tilting angle
  Float_t padlength[kNPlanes]; // pad lenghts
  AliTRDpadPlane *pp = 0x0;
  for(int iplane=0; iplane<kNPlanes; iplane++){
    pp                = fGeom->GetPadPlane(iplane, istack);
    hL[iplane]        = TMath::Tan(TMath::DegToRad()*pp->GetTiltingAngle());
    padlength[iplane] = pp->GetLengthIPad();
  }
  
  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){
    AliInfo(Form("Making seeds Stack[%d] Config[%d] Tracks[%d]...", istack, config, ntracks));
  }

  // Build seeding layers
  ResetSeedTB();
  Int_t nlayers = 0;
  for(int isl=0; isl<kNSeedPlanes; isl++){ 
    if(!(chamber = stack[planes[isl]])) continue;
    if(!chamber->GetSeedingLayer(fSeedTB[isl], fGeom, fReconstructor)) continue;
    nlayers++;
  }
  if(nlayers < 4) return ntracks;
  
  
  // Start finding seeds
  Double_t cond0[4], cond1[4], cond2[4];
  Int_t icl = 0;
  while((c[3] = (*fSeedTB[3])[icl++])){
    if(!c[3]) continue;
    fSeedTB[0]->BuildCond(c[3], cond0, 0);
    fSeedTB[0]->GetClusters(cond0, index, ncl);
    //printf("Found c[3] candidates 0 %d\n", ncl);
    Int_t jcl = 0;
    while(jcl<ncl) {
      c[0] = (*fSeedTB[0])[index[jcl++]];
      if(!c[0]) continue;
      Double_t dx    = c[3]->GetX() - c[0]->GetX();
      Double_t theta = (c[3]->GetZ() - c[0]->GetZ())/dx;
      Double_t phi   = (c[3]->GetY() - c[0]->GetY())/dx;
      fSeedTB[1]->BuildCond(c[0], cond1, 1, theta, phi);
      fSeedTB[1]->GetClusters(cond1, jndex, mcl);
      //printf("Found c[0] candidates 1 %d\n", mcl);

      Int_t kcl = 0;
      while(kcl<mcl) {
        c[1] = (*fSeedTB[1])[jndex[kcl++]];
        if(!c[1]) continue;
        fSeedTB[2]->BuildCond(c[1], cond2, 2, theta, phi);
        c[2] = fSeedTB[2]->GetNearestCluster(cond2);
        //printf("Found c[1] candidate 2 %p\n", c[2]);
        if(!c[2]) continue;
              
        // 				AliInfo("Seeding clusters found. Building seeds ...");
        // 				for(Int_t i = 0; i < kNSeedPlanes; i++) printf("%i. coordinates: x = %6.3f, y = %6.3f, z = %6.3f\n", i, c[i]->GetX(), c[i]->GetY(), c[i]->GetZ());
              
        for (Int_t il = 0; il < kNPlanes; il++) cseed[il].Reset();
      
        FitRieman(c, chi2);
      
        AliTRDseedV1 *tseed = &cseed[0];
        AliTRDtrackingChamber **cIter = &stack[0];
        for(int iLayer=0; iLayer<kNPlanes; iLayer++, tseed++, cIter++){
          tseed->SetDetector((*cIter) ? (*cIter)->GetDetector() : -1);
          tseed->SetTilt(hL[iLayer]);
          tseed->SetPadLength(padlength[iLayer]);
          tseed->SetReconstructor(fReconstructor);
          tseed->SetX0((*cIter) ? (*cIter)->GetX() : x_def[iLayer]);
          tseed->Init(GetRiemanFitter());
        }
      
        Bool_t isFake = kFALSE;
        if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) >= 2){
          if (c[0]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
          if (c[1]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
          if (c[2]->GetLabel(0) != c[3]->GetLabel(0)) isFake = kTRUE;
      
          Double_t xpos[4];
          for(Int_t l = 0; l < kNSeedPlanes; l++) xpos[l] = fSeedTB[l]->GetX();
          Float_t yref[4];
          for(int il=0; il<4; il++) yref[il] = cseed[planes[il]].GetYref(0);
          Int_t ll = c[3]->GetLabel(0);
          Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
          Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
          AliRieman *rim = GetRiemanFitter();
          TTreeSRedirector &cs0 = *fgDebugStreamer;
          cs0 << "MakeSeeds0"
              <<"EventNumber="		<< eventNumber
              <<"CandidateNumber="	<< candidateNumber
              <<"isFake="				<< isFake
              <<"config="				<< config
              <<"label="				<< ll
              <<"chi2z="				<< chi2[0]
              <<"chi2y="				<< chi2[1]
              <<"Y2exp="				<< cond2[0]	
              <<"Z2exp="				<< cond2[1]
              <<"X0="					<< xpos[0] //layer[sLayer]->GetX()
              <<"X1="					<< xpos[1] //layer[sLayer + 1]->GetX()
              <<"X2="					<< xpos[2] //layer[sLayer + 2]->GetX()
              <<"X3="					<< xpos[3] //layer[sLayer + 3]->GetX()
              <<"yref0="				<< yref[0]
              <<"yref1="				<< yref[1]
              <<"yref2="				<< yref[2]
              <<"yref3="				<< yref[3]
              <<"c0.="				<< c[0]
              <<"c1.="				<< c[1]
              <<"c2.="				<< c[2]
              <<"c3.="				<< c[3]
              <<"Seed0.="				<< &cseed[planes[0]]
              <<"Seed1.="				<< &cseed[planes[1]]
              <<"Seed2.="				<< &cseed[planes[2]]
              <<"Seed3.="				<< &cseed[planes[3]]
              <<"RiemanFitter.="		<< rim
              <<"\n";
        }
        if(chi2[0] > fReconstructor->GetRecoParam() ->GetChi2Z()/*7./(3. - sLayer)*//*iter*/){
//          //AliInfo(Form("Failed chi2 filter on chi2Z [%f].", chi2[0]));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        if(chi2[1] > fReconstructor->GetRecoParam() ->GetChi2Y()/*1./(3. - sLayer)*//*iter*/){
//          //AliInfo(Form("Failed chi2 filter on chi2Y [%f].", chi2[1]));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        //AliInfo("Passed chi2 filter.");
      
        // try attaching clusters to tracklets
        Int_t nUsedCl = 0;
        Int_t mlayers = 0;
        for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
          Int_t jLayer = planes[iLayer];
          if(!cseed[jLayer].AttachClustersIter(stack[jLayer], 5., kFALSE, c[iLayer])) continue;
          nUsedCl += cseed[jLayer].GetNUsed();
          if(nUsedCl > 25) break;
          mlayers++;
        }

        if(mlayers < kNSeedPlanes){ 
          //AliInfo(Form("Failed updating all seeds %d [%d].", mlayers, kNSeedPlanes));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }

        // temporary exit door for the HLT
        if(fReconstructor->IsHLT()){ 
          // attach clusters to extrapolation chambers
          for(int iLayer=0; iLayer<kNPlanes-kNSeedPlanes; iLayer++){
            Int_t jLayer = planesExt[iLayer];
            if(!(chamber = stack[jLayer])) continue;
            cseed[jLayer].AttachClustersIter(chamber, 1000.);
          }
          fTrackQuality[ntracks] = 1.; // dummy value
          ntracks++;
          if(ntracks == kMaxTracksStack) return ntracks;
          cseed += 6; 
          continue;
        }


        // fit tracklets and cook likelihood
        FitTiltedRieman(&cseed[0], kTRUE);// Update Seeds and calculate Likelihood
        Double_t like = CookLikelihood(&cseed[0], planes); // to be checked
      
        if (TMath::Log(1.E-9 + like) < fReconstructor->GetRecoParam() ->GetTrackLikelihood()){
          //AliInfo(Form("Failed likelihood %f[%e].", TMath::Log(1.E-9 + like), like));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        //AliInfo(Form("Passed likelihood %f[%e].", TMath::Log(1.E-9 + like), like));
      
        // book preliminary results
        seedQuality[ntracks] = like;
        fSeedLayer[ntracks]  = config;/*sLayer;*/
      
        // attach clusters to the extrapolation seeds
        Int_t nusedf   = 0; // debug value
        for(int iLayer=0; iLayer<kNPlanes-kNSeedPlanes; iLayer++){
          Int_t jLayer = planesExt[iLayer];
          if(!(chamber = stack[jLayer])) continue;
      
          // fit extrapolated seed
          if ((jLayer == 0) && !(cseed[1].IsOK())) continue;
          if ((jLayer == 5) && !(cseed[4].IsOK())) continue;
          AliTRDseedV1 pseed = cseed[jLayer];
          if(!pseed.AttachClustersIter(chamber, 1000.)) continue;
          cseed[jLayer] = pseed;
          nusedf += cseed[jLayer].GetNUsed(); // debug value
          FitTiltedRieman(cseed,  kTRUE);
        }
      
        // AliInfo("Extrapolation done.");
        // Debug Stream containing all the 6 tracklets
        if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) >= 2){
          TTreeSRedirector &cstreamer = *fgDebugStreamer;
          TLinearFitter *tiltedRieman = GetTiltedRiemanFitter();
          Int_t eventNumber 		= AliTRDtrackerDebug::GetEventNumber();
          Int_t candidateNumber	= AliTRDtrackerDebug::GetCandidateNumber();
          cstreamer << "MakeSeeds1"
              << "EventNumber="		<< eventNumber
              << "CandidateNumber="	<< candidateNumber
              << "S0.=" 				<< &cseed[0]
              << "S1.=" 				<< &cseed[1]
              << "S2.=" 				<< &cseed[2]
              << "S3.=" 				<< &cseed[3]
              << "S4.=" 				<< &cseed[4]
              << "S5.=" 				<< &cseed[5]
              << "FitterT.="			<< tiltedRieman
              << "\n";
        }
              
        if(fReconstructor->GetRecoParam()->HasImproveTracklets() && ImproveSeedQuality(stack, cseed) < 4){
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        //AliInfo("Improve seed quality done.");
      
        // fit full track and cook likelihoods
        // 				Double_t curv = FitRieman(&cseed[0], chi2);
        // 				Double_t chi2ZF = chi2[0] / TMath::Max((mlayers - 3.), 1.);
        // 				Double_t chi2RF = chi2[1] / TMath::Max((mlayers - 3.), 1.);
      
        // do the final track fitting (Once with vertex constraint and once without vertex constraint)
        Double_t chi2Vals[3];
        chi2Vals[0] = FitTiltedRieman(&cseed[0], kFALSE);
        if(fReconstructor->GetRecoParam()->IsVertexConstrained())
          chi2Vals[1] = FitTiltedRiemanConstraint(&cseed[0], GetZ()); // Do Vertex Constrained fit if desired
        else
          chi2Vals[1] = 1.;
        chi2Vals[2] = GetChi2Z(&cseed[0]) / TMath::Max((mlayers - 3.), 1.);
        // Chi2 definitions in testing stage
        //chi2Vals[2] = GetChi2ZTest(&cseed[0]);
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
          cseed[iLayer].SetChi2Z(chi2[1]);
        }
            
        if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) >= 2){
          TTreeSRedirector &cstreamer = *fgDebugStreamer;
          Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
          Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
          TLinearFitter *fitterTC = GetTiltedRiemanFitterConstraint();
          TLinearFitter *fitterT = GetTiltedRiemanFitter();
          Int_t ncls = 0; 
          for(Int_t iseed = 0; iseed < kNPlanes; iseed++){
          	ncls += cseed[iseed].IsOK() ? cseed[iseed].GetN2() : 0;
          }
          cstreamer << "MakeSeeds2"
              << "EventNumber=" 		<< eventNumber
              << "CandidateNumber="	<< candidateNumber
              << "Chi2TR="			<< chi2Vals[0]
              << "Chi2TC="			<< chi2Vals[1]
              << "Nlayers="			<< mlayers
              << "NClusters="   << ncls
              << "NUsedS="			<< nUsedCl
              << "NUsed="				<< nusedf
              << "Like="				<< like
              << "S0.="				<< &cseed[0]
              << "S1.="				<< &cseed[1]
              << "S2.="				<< &cseed[2]
              << "S3.="				<< &cseed[3]
              << "S4.="				<< &cseed[4]
              << "S5.="				<< &cseed[5]
              << "Label="				<< label
              << "Freq="				<< frequency
              << "FitterT.="			<< fitterT
              << "FitterTC.="			<< fitterTC
              << "\n";
        }
              
        ntracks++;
        AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
        if(ntracks == kMaxTracksStack){
          AliWarning(Form("Number of seeds reached maximum allowed (%d) in stack.", kMaxTracksStack));
          return ntracks;
        }
        cseed += 6;
      }
    }
  }
  
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

  AliTRDtrackV1 track(seeds, &params[1], c, params[0], params[6]*alpha+shift);
  track.PropagateTo(params[0]-5.0);
  if(fReconstructor->IsHLT()){ 
    AliTRDseedV1 *ptrTracklet = 0x0;
    for(Int_t ip=0; ip<kNPlanes; ip++){
      track.UnsetTracklet(ip);
      ptrTracklet = SetTracklet(&seeds[ip]);
      track.SetTracklet(ptrTracklet, fTracklets->GetEntriesFast()-1);
    }
    AliTRDtrackV1 *ptrTrack = SetTrack(&track);
    ptrTrack->SetReconstructor(fReconstructor);
    return ptrTrack;
  }

  track.ResetCovariance(1);
  Int_t nc = TMath::Abs(FollowBackProlongation(track));
  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 5){
    Int_t eventNumber 		= AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    Double_t p[5]; // Track Params for the Debug Stream
    track.GetExternalParameters(params[0], p);
    TTreeSRedirector &cs = *fgDebugStreamer;
    cs << "MakeTrack"
    << "EventNumber="     << eventNumber
    << "CandidateNumber=" << candidateNumber
    << "nc="     << nc
    << "X="      << params[0]
    << "Y="      << p[0]
    << "Z="      << p[1]
    << "snp="    << p[2]
    << "tnd="    << p[3]
    << "crv="    << p[4]
    << "Yin="    << params[1]
    << "Zin="    << params[2]
    << "snpin="  << params[3]
    << "tndin="  << params[4]
    << "crvin="  << params[5]
    << "track.=" << &track
    << "\n";
  }
  if (nc < 30) return 0x0;

  AliTRDtrackV1 *ptrTrack = SetTrack(&track);
  ptrTrack->SetReconstructor(fReconstructor);
  ptrTrack->CookLabel(.9);
  
  // computes PID for track
  ptrTrack->CookPID();
  // update calibration references using this track
  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if (!calibra){ 
    AliInfo("Could not get Calibra instance\n");
    if(calibra->GetHisto2d()) calibra->UpdateHistogramsV1(ptrTrack);
  }
  return ptrTrack;
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
  // debug level: 7
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
      squality[jLayer]  = bseed[jLayer].IsOK() ? bseed[jLayer].GetQuality(kTRUE) : 1000.;
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
    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) >= 7){
      Int_t eventNumber 		= AliTRDtrackerDebug::GetEventNumber();
      Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
      TLinearFitter *tiltedRieman = GetTiltedRiemanFitter();
      TTreeSRedirector &cstreamer = *fgDebugStreamer;
      cstreamer << "ImproveSeedQuality"
    << "EventNumber=" 		<< eventNumber
    << "CandidateNumber="	<< candidateNumber
    << "Iteration="				<< iter
    << "S0.="							<< &bseed[0]
    << "S1.="							<< &bseed[1]
    << "S2.="							<< &bseed[2]
    << "S3.="							<< &bseed[3]
    << "S4.="							<< &bseed[4]
    << "S5.="							<< &bseed[5]
    << "FitterT.="				<< tiltedRieman
    << "\n";
    }
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
  Double_t likeChi2TC = (fReconstructor->GetRecoParam() ->IsVertexConstrained()) ? 
  											TMath::Exp(-chi2[1] * 0.677) : 1;			// Constrained Tilted Riemann
  Double_t likeChi2TR = TMath::Exp(-chi2[0] * 0.78);			// Non-constrained Tilted Riemann
  Double_t likeAF     = TMath::Exp(-sumdaf * 3.23);
  Double_t trackLikelihood     = likeChi2Z * likeChi2TR * likeAF;

  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) >= 2){
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    TTreeSRedirector &cstreamer = *fgDebugStreamer;
    cstreamer << "CalculateTrackLikelihood0"
        << "EventNumber="			<< eventNumber
        << "CandidateNumber="	<< candidateNumber
        << "LikeChi2Z="				<< likeChi2Z
        << "LikeChi2TR="			<< likeChi2TR
        << "LikeChi2TC="			<< likeChi2TC
        << "LikeAF="					<< likeAF
        << "TrackLikelihood=" << trackLikelihood
        << "\n";
  }

  return trackLikelihood;
}

//____________________________________________________________________
Double_t AliTRDtrackerV1::CookLikelihood(AliTRDseedV1 *cseed, Int_t planes[4])
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

  // ratio of the total number of clusters/track which are expected to be found by the tracker.
  const AliTRDrecoParam *fRecoPars = fReconstructor->GetRecoParam();
  
 	Double_t chi2y = GetChi2Y(&cseed[0]);
  Double_t chi2z = GetChi2Z(&cseed[0]);

  Float_t nclusters = 0.;
  Double_t sumda = 0.;
  for(UChar_t ilayer = 0; ilayer < 4; ilayer++){
    Int_t jlayer = planes[ilayer];
    nclusters += cseed[jlayer].GetN2();
    sumda += TMath::Abs(cseed[jlayer].GetYfitR(1) - cseed[jlayer].GetYref(1));
  }
  nclusters *= .25;

  Double_t likea     = TMath::Exp(-sumda * fRecoPars->GetPhiSlope());
  Double_t likechi2y  = 0.0000000001;
  if (fReconstructor->IsCosmic() || chi2y < fRecoPars->GetChi2YCut()) likechi2y += TMath::Exp(-TMath::Sqrt(chi2y) * fRecoPars->GetChi2YSlope());
  Double_t likechi2z = TMath::Exp(-chi2z * fRecoPars->GetChi2ZSlope());
  Double_t likeN     = TMath::Exp(-(fRecoPars->GetNMeanClusters() - nclusters) / fRecoPars->GetNSigmaClusters());
  Double_t like      = likea * likechi2y * likechi2z * likeN;

  //	AliInfo(Form("sumda(%f) chi2[0](%f) chi2[1](%f) likea(%f) likechi2y(%f) likechi2z(%f) nclusters(%d) likeN(%f)", sumda, chi2[0], chi2[1], likea, likechi2y, likechi2z, nclusters, likeN));
  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) >= 2){
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    Int_t nTracklets = 0; Float_t mean_ncls = 0;
    for(Int_t iseed=0; iseed < kNPlanes; iseed++){
    	if(!cseed[iseed].IsOK()) continue;
    	nTracklets++;
    	mean_ncls += cseed[iseed].GetN2();
    }
    if(nTracklets) mean_ncls /= nTracklets;
    // The Debug Stream contains the seed 
    TTreeSRedirector &cstreamer = *fgDebugStreamer;
    cstreamer << "CookLikelihood"
        << "EventNumber="			<< eventNumber
        << "CandidateNumber=" << candidateNumber
        << "tracklet0.="			<< &cseed[0]
        << "tracklet1.="			<< &cseed[1]
        << "tracklet2.="			<< &cseed[2]
        << "tracklet3.="			<< &cseed[3]
        << "tracklet4.="			<< &cseed[4]
        << "tracklet5.="			<< &cseed[5]
        << "sumda="						<< sumda
        << "chi2y="						<< chi2y
        << "chi2z="						<< chi2z
        << "likea="						<< likea
        << "likechi2y="				<< likechi2y
        << "likechi2z="				<< likechi2z
        << "nclusters="				<< nclusters
        << "likeN="						<< likeN
        << "like="						<< like
        << "meanncls="        << mean_ncls
        << "\n";
  }

  return like;
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
  return idx >= 0 && idx < ncls ? (AliCluster*)fClusters->UncheckedAt(idx) : 0x0;
}

//____________________________________________________________________
AliTRDseedV1* AliTRDtrackerV1::GetTracklet(Int_t idx) const
{
  Int_t ntrklt = fTracklets->GetEntriesFast();
  return idx >= 0 && idx < ntrklt ? (AliTRDseedV1*)fTracklets->UncheckedAt(idx) : 0x0;
}

//____________________________________________________________________
AliKalmanTrack* AliTRDtrackerV1::GetTrack(Int_t idx) const
{
  Int_t ntrk = fTracks->GetEntriesFast();
  return idx >= 0 && idx < ntrk ? (AliKalmanTrack*)fTracks->UncheckedAt(idx) : 0x0;
}

//____________________________________________________________________
Float_t AliTRDtrackerV1::CalculateReferenceX(AliTRDseedV1 *tracklets){
  //
  // Calculates the reference x-position for the tilted Rieman fit defined as middle
  // of the stack (middle between layers 2 and 3). For the calculation all the tracklets
  // are taken into account
  // 
  // Parameters:	- Array of tracklets(AliTRDseedV1)
  //
  // Output:		- The reference x-position(Float_t)
  //
  Int_t nDistances = 0;
  Float_t meanDistance = 0.;
  Int_t startIndex = 5;
  for(Int_t il =5; il > 0; il--){
    if(tracklets[il].IsOK() && tracklets[il -1].IsOK()){
      Float_t xdiff = tracklets[il].GetX0() - tracklets[il -1].GetX0();
      meanDistance += xdiff;
      nDistances++;
    }
    if(tracklets[il].IsOK()) startIndex = il;
  }
  if(tracklets[0].IsOK()) startIndex = 0;
  if(!nDistances){
    // We should normally never get here
    Float_t xpos[2]; memset(xpos, 0, sizeof(Float_t) * 2);
    Int_t iok = 0, idiff = 0;
    // This attempt is worse and should be avoided:
    // check for two chambers which are OK and repeat this without taking the mean value
    // Strategy avoids a division by 0;
    for(Int_t il = 5; il >= 0; il--){
      if(tracklets[il].IsOK()){
  xpos[iok] = tracklets[il].GetX0();
  iok++;
  startIndex = il;
      }
      if(iok) idiff++;	// to get the right difference;
      if(iok > 1) break;
    }
    if(iok > 1){
      meanDistance = (xpos[0] - xpos[1])/idiff;
    }
    else{
      // we have do not even have 2 layers which are OK? The we do not need to fit at all
      return 331.;
    }
  }
  else{
    meanDistance /= nDistances;
  }
  return tracklets[startIndex].GetX0() + (2.5 - startIndex) * meanDistance - 0.5 * (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());
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


//____________________________________________________________________
void AliTRDtrackerV1::SetReconstructor(const AliTRDReconstructor *rec)
{
  fReconstructor = rec;
  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker) > 1){
    if(!fgDebugStreamer){
      TDirectory *savedir = gDirectory;
      fgDebugStreamer = new TTreeSRedirector("TRD.TrackerDebug.root");
      savedir->cd();
    }
  }	
}

//_____________________________________________________________________________
Float_t AliTRDtrackerV1::GetChi2Y(AliTRDseedV1 *tracklets) const
{
  //	Chi2 definition on y-direction

  Float_t chi2 = 0;
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    Double_t distLayer = (tracklets[ipl].GetYfit(0) - tracklets[ipl].GetYref(0));// /tracklets[ipl].GetSigmaY(); 
    chi2 += distLayer * distLayer;
  }
  return chi2;
}

//____________________________________________________________________
void AliTRDtrackerV1::ResetSeedTB()
{
// reset buffer for seeding time bin layers. If the time bin 
// layers are not allocated this function allocates them  

  for(Int_t isl=0; isl<kNSeedPlanes; isl++){
    if(!fSeedTB[isl]) fSeedTB[isl] = new AliTRDchamberTimeBin();
    else fSeedTB[isl]->Clear();
  }
}

//_____________________________________________________________________________
Float_t AliTRDtrackerV1::GetChi2Z(AliTRDseedV1 *tracklets) const 
{
  //	Calculates normalized chi2 in z-direction

  Float_t chi2 = 0;
  // chi2 = Sum ((z - zmu)/sigma)^2
  // Sigma for the z direction is defined as half of the padlength
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    Double_t distLayer = (tracklets[ipl].GetMeanz() - tracklets[ipl].GetZref(0)); // /(tracklets[ipl].GetPadLength()/2); 
    chi2 += distLayer * distLayer;
  }
  return chi2;
}

///////////////////////////////////////////////////////
//                                                   //
// Resources of class AliTRDLeastSquare              //
//                                                   //
///////////////////////////////////////////////////////

//_____________________________________________________________________________
AliTRDtrackerV1::AliTRDLeastSquare::AliTRDLeastSquare(){
  //
  // Constructor of the nested class AliTRDtrackFitterLeastSquare
  //
  memset(fParams, 0, sizeof(Double_t) * 2);
  memset(fSums, 0, sizeof(Double_t) * 5);
  memset(fCovarianceMatrix, 0, sizeof(Double_t) * 3);

}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::AddPoint(Double_t *x, Double_t y, Double_t sigmaY){
  //
  // Adding Point to the fitter
  //
  Double_t weight = 1/(sigmaY * sigmaY);
  Double_t &xpt = *x;
  //	printf("Adding point x = %f, y = %f, sigma = %f\n", xpt, y, sigmaY);
  fSums[0] += weight;
  fSums[1] += weight * xpt;
  fSums[2] += weight * y;
  fSums[3] += weight * xpt * y;
  fSums[4] += weight * xpt * xpt;
  fSums[5] += weight * y * y;
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::RemovePoint(Double_t *x, Double_t y, Double_t sigmaY){
  //
  // Remove Point from the sample
  //
  Double_t weight = 1/(sigmaY * sigmaY);
  Double_t &xpt = *x; 
  fSums[0] -= weight;
  fSums[1] -= weight * xpt;
  fSums[2] -= weight * y;
  fSums[3] -= weight * xpt * y;
  fSums[4] -= weight * xpt * xpt;
  fSums[5] -= weight * y * y;
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::Eval(){
  //
  // Evaluation of the fit:
  // Calculation of the parameters
  // Calculation of the covariance matrix
  //
  
  Double_t denominator = fSums[0] * fSums[4] - fSums[1] *fSums[1];
  if(denominator==0) return;

  //	for(Int_t isum = 0; isum < 5; isum++)
  //		printf("fSums[%d] = %f\n", isum, fSums[isum]);
  //	printf("denominator = %f\n", denominator);
  fParams[0] = (fSums[2] * fSums[4] - fSums[1] * fSums[3])/ denominator;
  fParams[1] = (fSums[0] * fSums[3] - fSums[1] * fSums[2]) / denominator;
  //	printf("fParams[0] = %f, fParams[1] = %f\n", fParams[0], fParams[1]);
  
  // Covariance matrix
  fCovarianceMatrix[0] = fSums[4] - fSums[1] * fSums[1] / fSums[0];
  fCovarianceMatrix[1] = fSums[5] - fSums[2] * fSums[2] / fSums[0];
  fCovarianceMatrix[2] = fSums[3] - fSums[1] * fSums[2] / fSums[0];
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDLeastSquare::GetFunctionValue(Double_t *xpos) const {
  //
  // Returns the Function value of the fitted function at a given x-position
  //
  return fParams[0] + fParams[1] * (*xpos);
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::GetCovarianceMatrix(Double_t *storage) const {
  //
  // Copies the values of the covariance matrix into the storage
  //
  memcpy(storage, fCovarianceMatrix, sizeof(Double_t) * 3);
}

