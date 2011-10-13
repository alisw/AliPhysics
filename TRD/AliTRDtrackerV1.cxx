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

#include <TBranch.h>
#include <TDirectory.h>
#include <TLinearFitter.h>
#include <TTree.h>  
#include <TClonesArray.h>
#include <TTreeStream.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>

#include "AliLog.h"
#include "AliMathBase.h"
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
#include "AliTRDdigitsParam.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDtrackerDebug.h"
#include "AliTRDtrackingChamber.h"
#include "AliTRDchamberTimeBin.h"

ClassImp(AliTRDtrackerV1)
ClassImp(AliTRDtrackerV1::AliTRDLeastSquare)
ClassImp(AliTRDtrackerV1::AliTRDtrackFitterRieman)

AliTRDtrackerV1::ETRDtrackerV1BetheBloch AliTRDtrackerV1::fgBB = AliTRDtrackerV1::kGeant;
Double_t AliTRDtrackerV1::fgTopologicQA[kNConfigs] = {
  0.5112, 0.5112, 0.5112, 0.0786, 0.0786,
  0.0786, 0.0786, 0.0579, 0.0579, 0.0474,
  0.0474, 0.0408, 0.0335, 0.0335, 0.0335
};  
const Double_t AliTRDtrackerV1::fgkX0[kNPlanes]    = {
  300.2, 312.8, 325.4, 338.0, 350.6, 363.2};
// Number of Time Bins/chamber should be also stored independently by the traker
// (also in AliTRDReconstructor) in oder to be able to run HLT. Fix TODO
Int_t AliTRDtrackerV1::fgNTimeBins = 0;
AliRieman* AliTRDtrackerV1::fgRieman = NULL;
TLinearFitter* AliTRDtrackerV1::fgTiltedRieman = NULL;
TLinearFitter* AliTRDtrackerV1::fgTiltedRiemanConstrained = NULL;

//____________________________________________________________________
AliTRDtrackerV1::AliTRDtrackerV1(AliTRDReconstructor *rec) 
  :AliTracker()
  ,fkReconstructor(NULL)
  ,fkRecoParam(NULL)
  ,fGeom(NULL)
  ,fClusters(NULL)
  ,fTracklets(NULL)
  ,fTracks(NULL)
  ,fTracksESD(NULL)
  ,fSieveSeeding(0)
  ,fEventInFile(-1)
{
  //
  // Default constructor.
  // 
  
  SetReconstructor(rec); // initialize reconstructor

  // initialize geometry
  if(!AliGeomManager::GetGeometry()){
    AliFatal("Could not get geometry.");
  }
  fGeom = new AliTRDgeometry();
  fGeom->CreateClusterMatrixArray();
  TGeoHMatrix *matrix = NULL;
  Double_t loc[] = {0., 0., 0.};
  Double_t glb[] = {0., 0., 0.};
  for(Int_t ily=kNPlanes; ily--;){
    Int_t ism = 0;
    while(!(matrix = fGeom->GetClusterMatrix(AliTRDgeometry::GetDetector(ily, 2, ism)))) ism++;
    if(!matrix){
      AliError(Form("Could not get transformation matrix for layer %d. Use default.", ily));
      fR[ily] = fgkX0[ily];
      continue;
    }
    matrix->LocalToMaster(loc, glb);
    fR[ily] = glb[0]+ AliTRDgeometry::AnodePos()-.5*AliTRDgeometry::AmThick() - AliTRDgeometry::DrThick();
  }

  // initialize cluster containers
  for (Int_t isector = 0; isector < AliTRDgeometry::kNsector; isector++) new(&fTrSec[isector]) AliTRDtrackingSector(fGeom, isector);
  
  // initialize arrays
  memset(fTrackQuality, 0, kMaxTracksStack*sizeof(Double_t));
  memset(fSeedLayer, 0, kMaxTracksStack*sizeof(Int_t));
  memset(fSeedTB, 0, kNSeedPlanes*sizeof(AliTRDchamberTimeBin*));
  fTracksESD = new TClonesArray("AliESDtrack", 2*kMaxTracksStack);
  fTracksESD->SetOwner();
}

//____________________________________________________________________
AliTRDtrackerV1::~AliTRDtrackerV1()
{ 
  //
  // Destructor
  //
  
  if(fgRieman) delete fgRieman; fgRieman = NULL;
  if(fgTiltedRieman) delete fgTiltedRieman; fgTiltedRieman = NULL;
  if(fgTiltedRiemanConstrained) delete fgTiltedRiemanConstrained; fgTiltedRiemanConstrained = NULL;
  for(Int_t isl =0; isl<kNSeedPlanes; isl++) if(fSeedTB[isl]) delete fSeedTB[isl];
  if(fTracksESD){ fTracksESD->Delete(); delete fTracksESD; }
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

  if(!fkRecoParam){
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
  AliInfo(Form("Number of tracks: !TRDin[%d]", ntracks));
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
  Int_t det = tracklet->GetDetector();
  Int_t sec = fGeom->GetSector(det);
  Double_t alpha = (sec+.5)*AliTRDgeometry::GetAlpha(),
           sinA  = TMath::Sin(alpha),
           cosA  = TMath::Cos(alpha);
  Double_t local[3];
  local[0] = tracklet->GetX(); 
  local[1] = tracklet->GetY();
  local[2] = tracklet->GetZ();
  Double_t global[3];
  fGeom->RotateBack(det, local, global);

  Double_t cov2D[3]; Float_t cov[6];
  tracklet->GetCovAt(local[0], cov2D);
  cov[0] = cov2D[0]*sinA*sinA;
  cov[1] =-cov2D[0]*sinA*cosA;
  cov[2] =-cov2D[1]*sinA;
  cov[3] = cov2D[0]*cosA*cosA;
  cov[4] = cov2D[1]*cosA;
  cov[5] = cov2D[2];
  // store the global position of the tracklet and its covariance matrix in the track point 
  p.SetXYZ(global[0],global[1],global[2], cov);
  
  // setting volume id
  AliGeomManager::ELayerID iLayer = AliGeomManager::ELayerID(AliGeomManager::kTRD1+fGeom->GetLayer(det));
  Int_t    modId = fGeom->GetSector(det) * AliTRDgeometry::kNstack + fGeom->GetStack(det);
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
  if(!fgRieman) fgRieman = new AliRieman(AliTRDseedV1::kNtb * AliTRDgeometry::kNlayer);
  return fgRieman;
}
  
//_____________________________________________________________________________
Int_t AliTRDtrackerV1::PropagateBack(AliESDEvent *event) 
{
// Propagation of ESD tracks from TPC to TOF detectors and building of the TRD track. For building
// a TRD track an ESD track is used as seed. The informations obtained on the TRD track (measured points,
// covariance, PID, etc.) are than used to update the corresponding ESD track.
// Each track seed is first propagated to the geometrical limit of the TRD detector. 
// Its prolongation is searched in the TRD and if corresponding clusters are found tracklets are 
// constructed out of them (see AliTRDseedV1::AttachClusters()) and the track is updated. 
// Otherwise the ESD track is left unchanged.
// 
// The following steps are performed:
// 1. Selection of tracks based on the variance in the y-z plane.
// 2. Propagation to the geometrical limit of the TRD volume. If track propagation fails the AliESDtrack::kTRDStop is set.
// 3. Prolongation inside the fiducial volume (see AliTRDtrackerV1::FollowBackProlongation()) and marking
// the following status bits:
//   - AliESDtrack::kTRDin - if the tracks enters the TRD fiducial volume
//   - AliESDtrack::kTRDStop - if the tracks fails propagation
//   - AliESDtrack::kTRDbackup - if the tracks fulfills chi2 conditions and qualify for refitting
// 4. Writting to friends, PID, MC label, quality etc. Setting status bit AliESDtrack::kTRDout.
// 5. Propagation to TOF. If track propagation fails the AliESDtrack::kTRDStop is set.
//  

  if(!fClusters || !fClusters->GetEntriesFast()){ 
    AliInfo("No TRD clusters");
    return 0;
  }
  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance(); // Calibration monitor
  if (!calibra) AliInfo("Could not get Calibra instance");
  if (!fgNTimeBins) fgNTimeBins = fkReconstructor->GetNTimeBins(); 

  // Define scalers
  Int_t nFound   = 0, // number of tracks found
        nBacked  = 0, // number of tracks backed up for refit
        nSeeds   = 0, // total number of ESD seeds
        nTRDseeds= 0, // number of seeds in the TRD acceptance
        nTPCseeds= 0; // number of TPC seeds
  Float_t foundMin = 20.0;
  
  Float_t *quality = NULL;
  Int_t   *index   = NULL;
  fEventInFile  = event->GetEventNumberInFile();
  nSeeds   = event->GetNumberOfTracks();
  // Sort tracks according to quality 
  // (covariance in the yz plane)
  if(nSeeds){  
    quality = new Float_t[nSeeds];
    index   = new Int_t[4*nSeeds];
    for (Int_t iSeed = nSeeds; iSeed--;) {
      AliESDtrack *seed = event->GetTrack(iSeed);
      Double_t covariance[15];
      seed->GetExternalCovariance(covariance);
      quality[iSeed] = covariance[0] + covariance[2];
    }
    TMath::Sort(nSeeds, quality, index,kFALSE);
  }
  
  // Propagate all seeds
  Int_t   expectedClr;
  AliTRDtrackV1 track;
  for (Int_t iSeed = 0; iSeed < nSeeds; iSeed++) {
  
    // Get the seeds in sorted sequence
    AliESDtrack *seed = event->GetTrack(index[iSeed]);
    Float_t p4  = seed->GetC(seed->GetBz());
  
    // Check the seed status
    ULong_t status = seed->GetStatus();
    if ((status & AliESDtrack::kTPCout) == 0) continue;
    if ((status & AliESDtrack::kTRDout) != 0) continue;

    // Propagate to the entrance in the TRD mother volume
    track.~AliTRDtrackV1();
    new(&track) AliTRDtrackV1(*seed);
    if(AliTRDgeometry::GetXtrdBeg() > (AliTRDReconstructor::GetMaxStep() + track.GetX()) && !PropagateToX(track, AliTRDgeometry::GetXtrdBeg(), AliTRDReconstructor::GetMaxStep())){
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
      continue;
    }    
    if(!AdjustSector(&track)){
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
      continue;
    }
    if(TMath::Abs(track.GetSnp()) > AliTRDReconstructor::GetMaxSnp()) {
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
      continue;
    }
    nTPCseeds++;
    AliDebug(2, Form("TRD propagate TPC seed[%d] = %d.", iSeed, index[iSeed]));
    // store track status at TRD entrance
    seed->UpdateTrackParams(&track, AliESDtrack::kTRDbackup);

    // prepare track and do propagation in the TRD
    track.SetReconstructor(fkReconstructor);
    track.SetKink(Bool_t(seed->GetKinkIndex(0)));
    track.SetPrimary(status & AliESDtrack::kTPCin);
    expectedClr = FollowBackProlongation(track);
    // check if track entered the TRD fiducial volume
    if(track.GetTrackIn()){ 
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDin);
      nTRDseeds++;
    }
    // check if track was stopped in the TRD
    if (expectedClr<0){      
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
      continue;
    }

    if(expectedClr){
      nFound++;  
      // computes PID for track
      track.CookPID();
      // update calibration references using this track
      if(calibra->GetHisto2d()) calibra->UpdateHistogramsV1(&track);
      // save calibration object
      if (fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 0) { 
        AliTRDtrackV1 *calibTrack = new AliTRDtrackV1(track);
        calibTrack->SetOwner();
        seed->AddCalibObject(calibTrack);
      }
      //update ESD track
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDout);
      track.UpdateESDtrack(seed);
    }

    if ((TMath::Abs(track.GetC(track.GetBz()) - p4) / TMath::Abs(p4) < 0.2) ||(track.Pt() > 0.8)) {

      // Make backup for back propagation
      Int_t foundClr = track.GetNumberOfClusters();
      if (foundClr >= foundMin) {
        track.CookLabel(1. - AliTRDReconstructor::GetLabelFraction());
        //if(track.GetBackupTrack()) UseClusters(track.GetBackupTrack());

        // Sign only gold tracks
        if (track.GetChi2() / track.GetNumberOfClusters() < 4) {
          //if ((seed->GetKinkIndex(0)      ==   0) && (track.Pt() <  1.5)) UseClusters(&track);
        }
        Bool_t isGold = kFALSE;
  
        // Full gold track
        if (track.GetChi2() / track.GetNumberOfClusters() < 5) {
          if (track.GetBackupTrack()) seed->UpdateTrackParams(track.GetBackupTrack(),AliESDtrack::kTRDbackup);
          nBacked++;
          isGold = kTRUE;
        }
  
        // Almost gold track
        if ((!isGold)  && (track.GetNCross() == 0) &&	(track.GetChi2() / track.GetNumberOfClusters()  < 7)) {
          //seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
          if (track.GetBackupTrack()) seed->UpdateTrackParams(track.GetBackupTrack(),AliESDtrack::kTRDbackup);
          nBacked++;
          isGold = kTRUE;
        }
        
        if ((!isGold) && (track.GetBackupTrack())) {
          if ((track.GetBackupTrack()->GetNumberOfClusters() > foundMin) && ((track.GetBackupTrack()->GetChi2()/(track.GetBackupTrack()->GetNumberOfClusters()+1)) < 7)) {
            seed->UpdateTrackParams(track.GetBackupTrack(),AliESDtrack::kTRDbackup);
            nBacked++;
            isGold = kTRUE;
          }
        }
      }
    }
    
    // Propagation to the TOF
    if(!(seed->GetStatus()&AliESDtrack::kTRDStop)) {
      Int_t sm = track.GetSector();
      // default value in case we have problems with the geometry.
      Double_t xtof  = 371.; 
      //Calculate radial position of the beginning of the TOF
      //mother volume. In order to avoid mixing of the TRD 
      //and TOF modules some hard values are needed. This are:
      //1. The path to the TOF module.
      //2. The width of the TOF (29.05 cm)
      //(with the help of Annalisa de Caro Mar-17-2009)
      if(gGeoManager){
        gGeoManager->cd(Form("/ALIC_1/B077_1/BSEGMO%d_1/BTOF%d_1", sm, sm));
        TGeoHMatrix *m = NULL;
        Double_t loc[]={0., 0., -.5*29.05}, glob[3];
        
        if((m=gGeoManager->GetCurrentMatrix())){
          m->LocalToMaster(loc, glob);
          xtof = TMath::Sqrt(glob[0]*glob[0]+glob[1]*glob[1]);
        }
      }
      if(xtof > (AliTRDReconstructor::GetMaxStep() + track.GetX()) && !PropagateToX(track, xtof, AliTRDReconstructor::GetMaxStep())){
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
        continue;
      }
      if(!AdjustSector(&track)){ 
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
        continue;
      }
      if(TMath::Abs(track.GetSnp()) > AliTRDReconstructor::GetMaxSnp()){
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDStop);
        continue;
      }
      //seed->UpdateTrackParams(&track, AliESDtrack::kTRDout);
      // TODO obsolete - delete
      seed->SetTRDQuality(track.StatusForTOF()); 
    }
    seed->SetTRDBudget(track.GetBudget(0));
  }
  if(index) delete [] index;
  if(quality) delete [] quality;

  AliInfo(Form("Number of seeds: TPCout[%d] TRDin[%d]", nTPCseeds, nTRDseeds));
  AliInfo(Form("Number of tracks: TRDout[%d] TRDbackup[%d]", nFound, nBacked));

  // run stand alone tracking
  if (fkReconstructor->IsSeeding()) Clusters2Tracks(event);
  
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
  
  
  if(!fClusters || !fClusters->GetEntriesFast()){ 
    AliInfo("No TRD clusters");
    return 0;
  }
  AliTRDtrackV1 track;
  for (Int_t itrack = 0; itrack < event->GetNumberOfTracks(); itrack++) {
    AliESDtrack *seed = event->GetTrack(itrack);
    ULong_t status = seed->GetStatus();

    new(&track) AliTRDtrackV1(*seed);
    if (track.GetX() < 270.0) {
      seed->UpdateTrackParams(&track, AliESDtrack::kTRDbackup);
      continue;
    }

    // reject tracks which failed propagation in the TRD or
    // are produced by the TRD stand alone tracker
    if(!(status & AliESDtrack::kTRDout)) continue;
    if(!(status & AliESDtrack::kTRDin)) continue;
    nseed++; 

    track.ResetCovariance(50.0);

    // do the propagation and processing
    Bool_t kUPDATE = kFALSE;
    Double_t xTPC = 250.0;
    if(FollowProlongation(track)){	
      // Update the friend track
      if (fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 0){ 
        TObject *o = NULL; Int_t ic = 0;
        AliTRDtrackV1 *calibTrack = NULL; 
        while((o = seed->GetCalibObject(ic++))){
          if(!(calibTrack = dynamic_cast<AliTRDtrackV1*>(o))) continue;
          calibTrack->SetTrackOut(&track);
        }
      }

      // Prolongate to TPC
      if (PropagateToX(track, xTPC, AliTRDReconstructor::GetMaxStep())) { //  -with update
        seed->UpdateTrackParams(&track, AliESDtrack::kTRDrefit);
        found++;
        kUPDATE = kTRUE;
      }
    }
    
    // Prolongate to TPC without update
    if(!kUPDATE) {
      AliTRDtrackV1 tt(*seed);
      if (PropagateToX(tt, xTPC, AliTRDReconstructor::GetMaxStep())) seed->UpdateTrackParams(&tt, AliESDtrack::kTRDbackup);
    }
  }
  AliInfo(Form("Number of seeds: TRDout[%d]", nseed));
  AliInfo(Form("Number of tracks: TRDrefit[%d]", found));
  
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
  for (Int_t iplane = kNPlanes; iplane--;) {
    Int_t   index(-1);
    AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
    AliDebug(2, Form("Tracklet[%p] ly[%d] idx[%d]", (void*)tracklet, iplane, index));
    if(!tracklet) continue;
    if(!tracklet->IsOK()){ 
      AliDebug(1, Form("Tracklet Det[%d] !OK", tracklet->GetDetector()));
      continue;
    }
    Double_t x  = tracklet->GetX();//GetX0();
    // reject tracklets which are not considered for inward refit
    if(x > t.GetX()+AliTRDReconstructor::GetMaxStep()) continue;

    // append tracklet to track
    t.SetTracklet(tracklet, index);
    
    if (x < (t.GetX()-AliTRDReconstructor::GetMaxStep()) && !PropagateToX(t, x+AliTRDReconstructor::GetMaxStep(), AliTRDReconstructor::GetMaxStep())) break;
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

    Double_t cov[3]; tracklet->GetCovAt(x, cov);
    Double_t p[2] = { tracklet->GetY(), tracklet->GetZ()};
    Double_t chi2 = ((AliExternalTrackParam)t).GetPredictedChi2(p, cov);
    if (chi2 < 1e+10 && ((AliExternalTrackParam&)t).Update(p, cov)){ 
      // Register info to track
      t.SetNumberOfClusters();
      t.UpdateChi2(chi2);
      nClustersExpected += tracklet->GetN();
    }
  }

  if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 1){
    Int_t index;
    for(int iplane=0; iplane<AliTRDgeometry::kNlayer; iplane++){
      AliTRDseedV1 *tracklet = GetTracklet(&t, iplane, index);
      if(!tracklet) continue;
      t.SetTracklet(tracklet, index);
    }

    if(fkReconstructor->IsDebugStreaming()){
      Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
      TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
      AliTRDtrackV1 track(t);
      track.SetOwner();
      cstreamer << "FollowProlongation"
          << "EventNumber="	<< eventNumber
          << "ncl="					<< nClustersExpected
          << "track.="			<< &track
          << "\n";
    }
  }
  return nClustersExpected;

}

//_____________________________________________________________________________
Int_t AliTRDtrackerV1::FollowBackProlongation(AliTRDtrackV1 &t)
{
// Extrapolates/Build the TRD track in the TOF direction.
//
// Parameters
//   t : the TRD track which has to be extrapolated
// 
// Output
//   number of clusters attached to the track
//
// Starting from current radial position of track <t> this function
// extrapolates the track through the 6 TRD layers. The following steps
// are being performed for each plane:
// 1. Propagate track to the entrance of the next chamber:
//   - get chamber limits in the radial direction
//   - check crossing sectors 
//   - check track inclination
//   - check track prolongation against boundary conditions (see exclusion boundaries on AliTRDgeometry::IsOnBoundary())
// 2. Build tracklet (see AliTRDseed::AttachClusters() for details) for this layer if needed. If only 
//    Kalman filter is needed and tracklets are already linked to the track this step is skipped.
// 3. Fit tracklet using the information from the Kalman filter.
// 4. Propagate and update track at reference radial position of the tracklet.
// 5. Register tracklet with the tracker and track; update pulls monitoring.
//
// Observation
//   1. During the propagation a bit map is filled detailing the status of the track in each TRD chamber. The following errors are being registered for each tracklet:
// - AliTRDtrackV1::kProlongation : track prolongation failed
// - AliTRDtrackV1::kPropagation : track prolongation failed
// - AliTRDtrackV1::kAdjustSector : failed during sector crossing
// - AliTRDtrackV1::kSnp : too large bending
// - AliTRDtrackV1::kTrackletInit : fail to initialize tracklet
// - AliTRDtrackV1::kUpdate : fail to attach clusters or fit the tracklet
// - AliTRDtrackV1::kUnknown : anything which is not covered before
//   2. By default the status of the track before first TRD update is saved. 
// 
// Debug level 2
//
// Author
//   Alexandru Bercuci <A.Bercuci@gsi.de>
//

  Int_t n = 0;
  Double_t driftLength = .5*AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick();
  AliTRDtrackingChamber *chamber = NULL;
  
  Int_t debugLevel = fkReconstructor->IsDebugStreaming() ? fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) : 0;
  TTreeSRedirector *cstreamer = fkReconstructor->IsDebugStreaming() ? fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker) : 0x0;

  Bool_t kStoreIn(kTRUE),     // toggel store track params. at TRD entry
         kStandAlone(kFALSE), // toggle tracker awarness of stand alone seeding 
         kUseTRD(fkRecoParam->IsOverPtThreshold(t.Pt()));// use TRD measurment to update Kalman

  Int_t startLayer(0);
  AliTRDseedV1 tracklet, *ptrTracklet = NULL;
  // Special case for stand alone tracking
  // - store all tracklets found by seeding
  // - start propagation from first tracklet found
  AliTRDseedV1 *tracklets[kNPlanes];
  memset(tracklets, 0, sizeof(AliTRDseedV1 *) * kNPlanes);
  for(Int_t ip(kNPlanes); ip--;){
    if(!(tracklets[ip] = t.GetTracklet(ip))) continue;
    t.UnsetTracklet(ip);
    if(tracklets[ip]->IsOK()) startLayer=ip;
    kStandAlone = kTRUE;
    kUseTRD = kTRUE;
  } 
  AliDebug(4, Form("SA[%c] Start[%d]\n"
    "  [0]idx[%d] traklet[%p]\n"
    "  [1]idx[%d] traklet[%p]\n"
    "  [2]idx[%d] traklet[%p]\n"
    "  [3]idx[%d] traklet[%p]\n"
    "  [4]idx[%d] traklet[%p]\n"
    "  [5]idx[%d] traklet[%p]"
    , kStandAlone?'y':'n', startLayer
    , t.GetTrackletIndex(0), (void*)tracklets[0]
    , t.GetTrackletIndex(1), (void*)tracklets[1]
    , t.GetTrackletIndex(2), (void*)tracklets[2]
    , t.GetTrackletIndex(3), (void*)tracklets[3]
    , t.GetTrackletIndex(4), (void*)tracklets[4]
    , t.GetTrackletIndex(5), (void*)tracklets[5]));

  // Loop through the TRD layers
  TGeoHMatrix *matrix = NULL;
  Double_t x(0.), y(0.), z(0.);
  for (Int_t ily=startLayer, sm=-1, stk=-1, det=-1; ily < AliTRDgeometry::kNlayer; ily++) {
    AliDebug(2, Form("Propagate to x[%d] = %7.2f", ily, fR[ily]));

    // rough estimate of the entry point
    if (!t.GetProlongation(fR[ily], y, z)){
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kProlongation);
      AliDebug(4, Form("Failed Rough Prolongation to ly[%d] x[%7.2f] y[%7.2f] z[%7.2f]", ily, fR[ily], y, z));
      break;
    }

    // find sector / stack / detector
    sm = t.GetSector();
    // TODO cross check with y value !
    stk = fGeom->GetStack(z, ily);
    det = stk>=0 ? AliTRDgeometry::GetDetector(ily, stk, sm) : -1;
    matrix = det>=0 ? fGeom->GetClusterMatrix(det) : NULL;

    // check if supermodule/chamber is installed
    if( !fGeom->GetSMstatus(sm) ||
        stk<0. ||
        fGeom->IsHole(ily, stk, sm) ||
        !matrix ){ 
      AliDebug(4, Form("Missing Geometry ly[%d]. Guess radial position", ily));
      // propagate to the default radial position
      if(fR[ily] > (AliTRDReconstructor::GetMaxStep() + t.GetX()) && !PropagateToX(t, fR[ily], AliTRDReconstructor::GetMaxStep())){
        n=-1; 
        t.SetStatus(AliTRDtrackV1::kPropagation);
        AliDebug(4, "Failed Propagation [Missing Geometry]");
        break;
      }
      if(!AdjustSector(&t)){
        n=-1; 
        t.SetStatus(AliTRDtrackV1::kAdjustSector);
        AliDebug(4, "Failed Adjust Sector [Missing Geometry]");
        break;
      }
      if(TMath::Abs(t.GetSnp()) > AliTRDReconstructor::GetMaxSnp()){
        n=-1; 
        t.SetStatus(AliTRDtrackV1::kSnp);
        AliDebug(4, "Failed Max Snp [Missing Geometry]");
        break;
      }
      t.SetStatus(AliTRDtrackV1::kGeometry, ily);
      continue;
    }

    // retrieve rotation matrix for the current chamber
    Double_t loc[] = {AliTRDgeometry::AnodePos()- driftLength, 0., 0.};
    Double_t glb[] = {0., 0., 0.};
    matrix->LocalToMaster(loc, glb);
    AliDebug(3, Form("Propagate to det[%3d] x_anode[%7.2f] (%f %f)", det, glb[0]+driftLength, glb[1], glb[2]));

    // Propagate to the radial distance of the current layer
    x = glb[0] - AliTRDReconstructor::GetMaxStep();
    if(x > (AliTRDReconstructor::GetMaxStep() + t.GetX()) && !PropagateToX(t, x, AliTRDReconstructor::GetMaxStep())){
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kPropagation);
      AliDebug(4, Form("Failed Initial Propagation to x[%7.2f]", x));
      break;
    }
    if(!AdjustSector(&t)){
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kAdjustSector);
      AliDebug(4, "Failed Adjust Sector Start");
      break;
    }
    if(TMath::Abs(t.GetSnp()) > AliTRDReconstructor::GetMaxSnp()) {
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kSnp);
      AliDebug(4, Form("Failed Max Snp[%f] MaxSnp[%f]", t.GetSnp(), AliTRDReconstructor::GetMaxSnp()));
      break;
    }
    Bool_t doRecalculate = kFALSE;
    if(sm != t.GetSector()){
      sm = t.GetSector(); 
      doRecalculate = kTRUE;
    }
    if(stk != fGeom->GetStack(z, ily)){
      stk = fGeom->GetStack(z, ily);
      doRecalculate = kTRUE;
    }
    if(doRecalculate){
      det = AliTRDgeometry::GetDetector(ily, stk, sm);
      if(!(matrix = fGeom->GetClusterMatrix(det))){ 
        t.SetStatus(AliTRDtrackV1::kGeometry, ily);
        AliDebug(4, Form("Failed Geometry Matrix ly[%d]", ily));
        continue;
      }
      matrix->LocalToMaster(loc, glb);
      x = glb[0] - AliTRDReconstructor::GetMaxStep();
    }

    // check if track is well inside fiducial volume 
    if (!t.GetProlongation(x+AliTRDReconstructor::GetMaxStep(), y, z)) {
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kProlongation);
      AliDebug(4, Form("Failed Prolongation to x[%7.2f] y[%7.2f] z[%7.2f]", x+AliTRDReconstructor::GetMaxStep(), y, z));
      break;
    }
    if(fGeom->IsOnBoundary(det, y, z, .5)){ 
      t.SetStatus(AliTRDtrackV1::kBoundary, ily);
      AliDebug(4, "Failed Track on Boundary");
      continue;
    }

    ptrTracklet  = tracklets[ily];
    if(!ptrTracklet){ // BUILD TRACKLET
      AliDebug(3, Form("Building tracklet det[%d]", det));
      // check data in supermodule
      if(!fTrSec[sm].GetNChambers()){ 
        t.SetStatus(AliTRDtrackV1::kNoClusters, ily);
        AliDebug(4, "Failed NoClusters");
        continue;
      }
      if(fTrSec[sm].GetX(ily) < 1.){ 
        t.SetStatus(AliTRDtrackV1::kNoClusters, ily);
        AliDebug(4, "Failed NoX");
        continue;
      }
      
      // check data in chamber
      if(!(chamber = fTrSec[sm].GetChamber(stk, ily))){ 
        t.SetStatus(AliTRDtrackV1::kNoClusters, ily);
        AliDebug(4, "Failed No Detector");
        continue;
      }
      if(chamber->GetNClusters() < fgNTimeBins*fkRecoParam ->GetFindableClusters()){ 
        t.SetStatus(AliTRDtrackV1::kNoClusters, ily);
        AliDebug(4, "Failed Not Enough Clusters in Detector");
        continue;
      }      
      // build tracklet
      tracklet.~AliTRDseedV1();
      ptrTracklet = new(&tracklet) AliTRDseedV1(det);
      ptrTracklet->SetReconstructor(fkReconstructor);
      ptrTracklet->SetKink(t.IsKink());
      ptrTracklet->SetPrimary(t.IsPrimary());
      ptrTracklet->SetPadPlane(fGeom->GetPadPlane(ily, stk));
      ptrTracklet->SetX0(glb[0]+driftLength);
      if(!ptrTracklet->Init(&t)){
        n=-1; 
        t.SetStatus(AliTRDtrackV1::kTrackletInit);
        AliDebug(4, "Failed Tracklet Init");
        break;
      }
      if(!ptrTracklet->AttachClusters(chamber, kTRUE, t.Charge()>0?kTRUE:kFALSE, fEventInFile)){
        t.SetStatus(AliTRDtrackV1::kNoAttach, ily);
        if(debugLevel>3){
          AliTRDseedV1 trackletCp(*ptrTracklet);
          UChar_t status(t.GetStatusTRD(ily));
          (*cstreamer)   << "FollowBackProlongation4"
          <<"status="    << status
          <<"tracklet.=" << &trackletCp
          << "\n";
        }
        AliDebug(4, "Failed Attach Clusters");
        continue;
      }
      AliDebug(3, Form("Number of Clusters in Tracklet: %d", ptrTracklet->GetN()));
      if(ptrTracklet->GetN() < fgNTimeBins*fkRecoParam->GetFindableClusters()){
        t.SetStatus(AliTRDtrackV1::kNoClustersTracklet, ily);
        if(debugLevel>3){
          AliTRDseedV1 trackletCp(*ptrTracklet);
          UChar_t status(t.GetStatusTRD(ily));
          (*cstreamer)   << "FollowBackProlongation4"
          <<"status="    << status
          <<"tracklet.=" << &trackletCp
          << "\n";
        }
        AliDebug(4, "Failed N Clusters Attached");
        continue;
      }
      ptrTracklet->UpdateUsed();
    } else AliDebug(2, Form("Use external tracklet ly[%d]", ily));
    // propagate track to the radial position of the tracklet

    // fit tracklet 
    // tilt correction options
    // 0 : no correction
    // 2 : pseudo tilt correction
    if(!ptrTracklet->FitRobust(t.Charge()>0?kTRUE:kFALSE)){
      t.SetStatus(AliTRDtrackV1::kNoFit, ily);
      AliDebug(4, "Failed Tracklet Fit");
      continue;
    } 
    x = ptrTracklet->GetX(); //GetX0();
    if(x > (AliTRDReconstructor::GetMaxStep() + t.GetX()) && !PropagateToX(t, x, AliTRDReconstructor::GetMaxStep())) {
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kPropagation);
      AliDebug(4, Form("Failed Propagation to Tracklet x[%7.2f]", x));
      break;
    }
    if(!AdjustSector(&t)) {
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kAdjustSector);
      AliDebug(4, "Failed Adjust Sector");
      break;
    }
    if(TMath::Abs(t.GetSnp()) > AliTRDReconstructor::GetMaxSnp()) {
      n=-1; 
      t.SetStatus(AliTRDtrackV1::kSnp);
      AliDebug(4, Form("Failed Max Snp[%f] MaxSnp[%f]", t.GetSnp(), AliTRDReconstructor::GetMaxSnp()));
      break;
    }
    Double_t cov[3]; ptrTracklet->GetCovAt(x, cov);
    Double_t p[2] = { ptrTracklet->GetY(), ptrTracklet->GetZ()};
    Double_t chi2 = ((AliExternalTrackParam)t).GetPredictedChi2(p, cov);
    // update Kalman with the TRD measurement
    if(chi2>1e+10){ // TODO
      t.SetStatus(AliTRDtrackV1::kChi2, ily);
      if(debugLevel > 2){
        UChar_t status(t.GetStatusTRD());
        AliTRDseedV1  trackletCp(*ptrTracklet);
        AliTRDtrackV1 trackCp(t);
        trackCp.SetOwner();
        (*cstreamer) << "FollowBackProlongation3"
            << "status="      << status
            << "tracklet.="   << &trackletCp
            << "track.="      << &trackCp
            << "\n";
      }
      AliDebug(4, Form("Failed Chi2[%f]", chi2));
      continue; 
    }
    // mark track as entering the FIDUCIAL volume of TRD
    if(kStoreIn){
      t.SetTrackIn();
      kStoreIn = kFALSE;
    }
    if(kUseTRD){
      if(!((AliExternalTrackParam&)t).Update(p, cov)) {
        n=-1; 
        t.SetStatus(AliTRDtrackV1::kUpdate);
        if(debugLevel > 2){
          UChar_t status(t.GetStatusTRD());
          AliTRDseedV1  trackletCp(*ptrTracklet);
          AliTRDtrackV1 trackCp(t);
          trackCp.SetOwner();
          (*cstreamer) << "FollowBackProlongation3"
              << "status="      << status
              << "tracklet.="   << &trackletCp
              << "track.="      << &trackCp
              << "\n";
        }
        AliDebug(4, Form("Failed Track Update @ y[%7.2f] z[%7.2f] s2y[%f] s2z[%f] covyz[%f]", p[0], p[1], cov[0], cov[2], cov[1]));
        break;
      }
    }
    if(!kStandAlone) ptrTracklet->UseClusters();
    // fill residuals ?!
    AliTracker::FillResiduals(&t, p, cov, ptrTracklet->GetVolumeId());
  

    // register tracklet with the tracker and track
    ptrTracklet->Update(&t);
    ptrTracklet = SetTracklet(ptrTracklet);
    Int_t index(fTracklets->GetEntriesFast()-1);
    t.SetTracklet(ptrTracklet, index);
    // Register info to track
    t.SetNumberOfClusters();
    t.UpdateChi2(chi2);

    n += ptrTracklet->GetN();
    AliDebug(2, Form("Setting Tracklet[%d] @ Idx[%d]", ily, index));

    // Reset material budget if 2 consecutive gold
//     if(ilayer>0 && t.GetTracklet(ilayer-1) && ptrTracklet->GetN() + t.GetTracklet(ilayer-1)->GetN() > 20) t.SetBudget(2, 0.);

    // Make backup of the track until is gold
    Int_t failed(0);
    if(!kStandAlone && (failed = t.MakeBackupTrack())) AliDebug(2, Form("Failed backup on cut[%d]", failed));

  } // end layers loop
  //printf("clusters[%d] chi2[%f] x[%f] status[%d ", n, t.GetChi2(), t.GetX(), t.GetStatusTRD());
  //for(int i=0; i<6; i++) printf("%d ", t.GetStatusTRD(i)); printf("]\n");

  if(n && debugLevel > 1){
    //Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    AliTRDtrackV1 track(t);
    track.SetOwner();
    (*cstreamer) << "FollowBackProlongation2"
        << "EventNumber=" << fEventInFile
        << "track.="      << &track
        << "\n";
  }
  
  return n;
}

//_________________________________________________________________________
Float_t AliTRDtrackerV1::FitRieman(AliTRDseedV1 *tracklets, Double_t *chi2, Int_t *const planes){
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
    fitter->AddPoint(tracklets[ppl[il]].GetX0(), tracklets[ppl[il]].GetYfit(0), tracklets[ppl[il]].GetZfit(0),1,10);
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
  for(Int_t i = 0; i < 4; i++){
    fitter->AddPoint(seedcl[i]->GetX(), seedcl[i]->GetY(), seedcl[i]->GetZ(), 1., 10.);
  }
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
  AliTRDcluster *cl = NULL;
  
  Float_t x, y, z, w, t, error, tilt;
  Double_t uvt[2];
  Int_t nPoints = 0;
  for(Int_t ilr = 0; ilr < AliTRDgeometry::kNlayer; ilr++){
    if(!tracklets[ilr].IsOK()) continue;
    for(Int_t itb = 0; itb < AliTRDseedV1::kNclusters; itb++){
      if(!tracklets[ilr].IsUsable(itb)) continue;
      if(!(cl = tracklets[ilr].GetClusters(itb))) continue;
      if(!cl->IsInChamber()) continue;
      x = cl->GetX();
      y = cl->GetY();
      z = cl->GetZ();
      tilt = tracklets[ilr].GetTilt();
      // Transformation
      t = 1./(x * x + y * y);
      uvt[0] = 2. * x * t;
      uvt[1] = 2. * x * t * tilt ;
      w = 2. * (y + tilt * (z - zVertex)) * t;
      error = 2. * TMath::Sqrt(cl->GetSigmaY2()+tilt*tilt*cl->GetSigmaZ2()) * t;
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
    tracklets[ip].SetC(curvature, 1);

  if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>3) printf("D-AliTRDtrackerV1::FitTiltedRiemanConstraint: Chi2[%f] C[%5.2e] pt[%8.3f]\n", chi2track, curvature, GetBz()*kB2C/curvature);

/*  if(fkReconstructor->GetRecoParam()->GetStreamLevel(AliTRDrecoParam::kTracker()) >= 5){
    //Linear Model on z-direction
    Double_t xref = CalculateReferenceX(tracklets);		// Relative to the middle of the stack
    Double_t slope = fitter->GetParameter(2);
    Double_t zref = slope * xref;
    Float_t chi2Z = CalculateChi2Z(tracklets, zref, slope, xref);
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    TTreeSRedirector &treeStreamer = *fkReconstructor->GetDebugStream(AliTRDReconstructor::kTracker);
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
  AliTRDcluster *cl = NULL;

  Double_t xref = CalculateReferenceX(tracklets);
  Double_t x, y, z, t, tilt, dx, w, we, erry, errz;
  Double_t uvt[4], sumPolY[5], sumPolZ[3];
  memset(sumPolY, 0, sizeof(Double_t) * 5);
  memset(sumPolZ, 0, sizeof(Double_t) * 3);
  Int_t nPoints = 0;
  // Containers for Least-square fitter
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    tilt = tracklets[ipl].GetTilt();
    for(Int_t itb = 0; itb < AliTRDseedV1::kNclusters; itb++){
      if(!(cl = tracklets[ipl].GetClusters(itb))) continue;
      if(!cl->IsInChamber()) continue;
      if (!tracklets[ipl].IsUsable(itb)) continue;
      x = cl->GetX();
      y = cl->GetY();
      z = cl->GetZ();
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
      we *= sigError ? TMath::Sqrt(cl->GetSigmaY2()+tilt*tilt*cl->GetSigmaZ2()) : 0.2;
      fitter->AddPoint(uvt, w, we);
      zfitter.AddPoint(&x, z, static_cast<Double_t>(TMath::Sqrt(cl->GetSigmaZ2())));
      // adding points for covariance matrix estimation
      erry = 1./(TMath::Sqrt(cl->GetSigmaY2()) + 0.1);  // 0.1 is a systematic error (due to misalignment and miscalibration)
      erry *= erry;
      errz = 1./cl->GetSigmaZ2();
      for(Int_t ipol = 0; ipol < 5; ipol++){
        sumPolY[ipol] += erry;
        erry *= x;
        if(ipol < 3){
          sumPolZ[ipol] += errz;
          errz *= x;
        }
      }
      nPoints++;
    }
  }
  if (fitter->Eval()) return 1.e10;
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
    if (TMath::Abs(tracklets[iLayer].GetZfit(0) - zref) > tracklets[iLayer].GetPadLength() * 0.5 + 1.0) 
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
  if (curvature > 0.0) curvature  =  a / TMath::Sqrt(curvature);

  Double_t chi2track = fitter->GetChisquare()/Double_t(nPoints);

  // Prepare error calculation
  TMatrixD covarPolY(3,3);
  covarPolY(0,0) = sumPolY[0]; covarPolY(1,1) = sumPolY[2]; covarPolY(2,2) = sumPolY[4];
  covarPolY(0,1) = covarPolY(1,0) = sumPolY[1];
  covarPolY(0,2) = covarPolY(2,0) = sumPolY[2];
  covarPolY(2,1) = covarPolY(1,2) = sumPolY[3];
  covarPolY.Invert();
  TMatrixD covarPolZ(2,2);
  covarPolZ(0,0) = sumPolZ[0]; covarPolZ(1,1) = sumPolZ[2];
  covarPolZ(1,0) = covarPolZ(0,1) = sumPolZ[1];
  covarPolZ.Invert();

  // Update the tracklets
  Double_t x1, dy, dz;
  Double_t cov[15];
  memset(cov, 0, sizeof(Double_t) * 15);
  for(Int_t iLayer = 0; iLayer < AliTRDtrackerV1::kNPlanes; iLayer++) {

    x  = tracklets[iLayer].GetX0();
    x1 = x - xref;
    y  = 0;
    z  = 0;
    dy = 0;
    dz = 0;
    memset(cov, 0, sizeof(Double_t) * 3);
    TMatrixD transform(3,3);
    transform(0,0) = 1;
    transform(0,1) = x;
    transform(0,2) = x*x;
    transform(1,1) = 1;
    transform(1,2) = x;
    transform(2,2) = 1;
    TMatrixD covariance(transform, TMatrixD::kMult, covarPolY);
    covariance *= transform.T();
    TMatrixD transformZ(2,2);
    transformZ(0,0) = transformZ(1,1) = 1;
    transformZ(0,1) = x;
    TMatrixD covarZ(transformZ, TMatrixD::kMult, covarPolZ);
    covarZ *= transformZ.T();
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
    cov[0] = covariance(0,0);
    cov[2] = covarZ(0,0);
    cov[1] = 0.;

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
    tracklets[iLayer].SetCovRef(cov);
    tracklets[iLayer].SetChi2(chi2track);
  }
  if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>3) printf("D-AliTRDtrackerV1::FitTiltedRieman: Chi2[%f] C[%5.2e] pt[%8.3f]\n", chi2track, curvature, GetBz()*kB2C/curvature);
  
/*  if(fkReconstructor->GetRecoParam()->GetStreamLevel(AliTRDrecoParam::kTracker) >=5){
    TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
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
  //
  // Fit track with a staight line
  // Fills an AliTrackPoint array with np points
  // Function should be used to refit tracks when no magnetic field was on
  //
  AliTRDLeastSquare yfitter, zfitter;
  AliTRDcluster *cl = NULL;

  AliTRDseedV1 work[kNPlanes], *tracklet = NULL;
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
//
// Paramters:   - Array of tracklets (connected to the track candidate)
//              - Flag selecting the error definition
// Output:      - Chi2 values of the track (in Parameter list)
//
// The equations which has to be solved simultaneously are:
// BEGIN_LATEX
// R^{2} = (x-x_{0})^{2} + (y^{*}-y_{0})^{2}
// y^{*} = y - tg(h)(z - z_{t})
// z_{t} = z_{0}+dzdx*(x-x_{r})
// END_LATEX
// with (x, y, z) the coordinate of the cluster, (x_0, y_0, z_0) the coordinate of the center of the Riemann circle,
// R its radius, x_r a constant refrence radial position in the middle of the TRD stack  and dzdx the slope of the 
// track in the x-z plane. Using the following transformations
// BEGIN_LATEX
// t = 1 / (x^{2} + y^{2})
// u = 2 * x * t
// v = 2 * tan(h) * t
// w = 2 * tan(h) * (x - x_{r}) * t
// END_LATEX
// One gets the following linear equation
// BEGIN_LATEX
// a + b * u + c * t + d * v  + e * w = 2 * (y + tg(h) * z) * t
// END_LATEX
// where the coefficients have the following meaning 
// BEGIN_LATEX
// a = -1/y_{0}
// b = x_{0}/y_{0}
// c = (R^{2} -x_{0}^{2} - y_{0}^{2})/y_{0}
// d = z_{0}
// e = dz/dx
// END_LATEX
// The error calculation for the free term is thus
// BEGIN_LATEX
// #sigma = 2 * #sqrt{#sigma^{2}_{y} + (tilt corr ...) + tg^{2}(h) * #sigma^{2}_{z}} * t
// END_LATEX
//
// From this simple model one can compute chi^2 estimates and a rough approximation of pt from the curvature according 
// to the formula:
// BEGIN_LATEX
// C = 1/R = a/(1 + b^{2} + c*a)
// END_LATEX
//
// Authors
//   M.Ivanov <M.Ivanov@gsi.de>
//   A.Bercuci <A.Bercuci@gsi.de>
//   M.Fasel <M.Fasel@gsi.de>

  TLinearFitter *fitter = GetTiltedRiemanFitter();
  fitter->StoreData(kTRUE);
  fitter->ClearPoints();
  AliTRDLeastSquare zfitter;
  AliTRDcluster *cl = NULL;

  AliTRDseedV1 work[kNPlanes], *tracklet = NULL;
  if(!tracklets){
    for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
      if(!(tracklet = track->GetTracklet(ipl))) continue;
      if(!tracklet->IsOK()) continue;
      new(&work[ipl]) AliTRDseedV1(*tracklet);
    }
    tracklets = &work[0];
  }

  Double_t xref = CalculateReferenceX(tracklets);
  if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>3) printf("D-AliTRDtrackerV1::FitRiemanTilt:\nx0[(0)%6.2f (1)%6.2f (2)%6.2f (3)%6.2f (4)%6.2f (5)%6.2f] xref[%6.2f]", tracklets[0].GetX0(), tracklets[1].GetX0(), tracklets[2].GetX0(), tracklets[3].GetX0(), tracklets[4].GetX0(), tracklets[5].GetX0(), xref);
  Double_t x, y, z, t, tilt, dx, w, we;
  Double_t uvt[4];
  Int_t nPoints = 0;
  // Containers for Least-square fitter
  for(Int_t ipl = 0; ipl < kNPlanes; ipl++){
    if(!tracklets[ipl].IsOK()) continue;
    for(Int_t itb = 0; itb < AliTRDseedV1::kNclusters; itb++){
      if(!(cl = tracklets[ipl].GetClusters(itb))) continue;
      //if (!tracklets[ipl].IsUsable(itb)) continue;
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
      we *= sigError ? TMath::Sqrt(cl->GetSigmaY2()) : 0.2;
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
    if (TMath::Abs(tracklets[iLayer].GetZfit(0) - zref) > tracklets[iLayer].GetPadLength() * 0.5 + 1.0) 
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
  Double_t radius    = TMath::Sqrt(tmp);
  Double_t curvature    =  1.0 + b*b - c*a;
  if (curvature > 0.0)  curvature  =  a / TMath::Sqrt(curvature);

  // Calculate chi2 of the fit 
  Double_t chi2 = fitter->GetChisquare()/Double_t(nPoints);
  if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>3) printf("D-AliTRDtrackerV1::FitRiemanTilt:x0[%6.2f] y0[%6.2f] R[%6.2f] chi2[%f]\n", x0, y0, radius, chi2);

  // Update the tracklets
  if(!track){
    for(Int_t ip = 0; ip < kNPlanes; ip++) {
      x = tracklets[ip].GetX0();
      tmp = radius*radius-(x-x0)*(x-x0);  
      if(tmp <= 0.) continue;
      tmp = TMath::Sqrt(tmp);  

      // y:     R^2 = (x - x0)^2 + (y - y0)^2
      //     =>   y = y0 +/- Sqrt(R^2 - (x - x0)^2)
      tracklets[ip].SetYref(0, y0 - (y0>0.?1.:-1)*tmp);
      //     => dy/dx = (x - x0)/Sqrt(R^2 - (x - x0)^2) 
      tracklets[ip].SetYref(1, (x - x0) / tmp);
      tracklets[ip].SetZref(0, z0 + dzdx * (x - xref));
      tracklets[ip].SetZref(1, dzdx);
      tracklets[ip].SetC(curvature);
      tracklets[ip].SetChi2(chi2);
    }
  }
  //update track points array
  if(np && points){
    Float_t xyz[3];
    for(int ip=0; ip<np; ip++){
      points[ip].GetXYZ(xyz);
      xyz[1] = TMath::Abs(xyz[0] - x0) > radius ? 100. : y0 - (y0>0.?1.:-1.)*TMath::Sqrt((radius-(xyz[0]-x0))*(radius+(xyz[0]-x0)));
      xyz[2] = z0 + dzdx * (xyz[0] - xref);
      points[ip].SetXYZ(xyz);
    }
  }
  
  return chi2;
}


//____________________________________________________________________
Double_t AliTRDtrackerV1::FitKalman(AliTRDtrackV1 *track, AliTRDseedV1 * const tracklets, Bool_t up, Int_t np, AliTrackPoint *points)
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


  AliTRDseedV1 tracklet;
  AliTRDseedV1 *ptrTracklet = NULL;

  //Loop through the TRD planes
  for (Int_t jplane = 0; jplane < kNPlanes; jplane++) {
    // GET TRACKLET OR BUILT IT		
    Int_t iplane = up ? jplane : kNPlanes - 1 - jplane;
    if(tracklets){ 
      if(!(ptrTracklet = &tracklets[iplane])) continue;
    }else{
      if(!(ptrTracklet  = track->GetTracklet(iplane))){ 
      /*AliTRDtrackerV1 *tracker = NULL;
        if(!(tracker = dynamic_cast<AliTRDtrackerV1*>( AliTRDrecoParam:Tracker()))) continue;
        ptrTracklet = new(&tracklet) AliTRDseedV1(iplane);
        if(!tracker->MakeTracklet(ptrTracklet, track)) */
        continue;
      }
    }
    if(!ptrTracklet->IsOK()) continue;

    Double_t x = ptrTracklet->GetX0();

    while(ip < np){
      //don't do anything if next marker is after next update point.
      if((up?-1:1) * (points[ip].GetX() - x) - AliTRDReconstructor::GetMaxStep() < 0) break;
      if(((up?-1:1) * (points[ip].GetX() - track->GetX()) < 0) && !PropagateToX(*track, points[ip].GetX(), AliTRDReconstructor::GetMaxStep())) return -1.;
      
      Double_t xyz[3]; // should also get the covariance
      track->GetXYZ(xyz);
      track->Global2LocalPosition(xyz, track->GetAlpha());
      points[ip].SetXYZ(xyz[0], xyz[1], xyz[2]);
      ip++;
    }
    // printf("plane[%d] tracklet[%p] x[%f]\n", iplane, ptrTracklet, x);

    // Propagate closer to the next update point 
    if(((up?-1:1) * (x - track->GetX()) + AliTRDReconstructor::GetMaxStep() < 0) && !PropagateToX(*track, x + (up?-1:1)*AliTRDReconstructor::GetMaxStep(), AliTRDReconstructor::GetMaxStep())) return -1.;

    if(!AdjustSector(track)) return -1;
    if(TMath::Abs(track->GetSnp()) > AliTRDReconstructor::GetMaxSnp()) return -1;
    
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
    if(TMath::Abs(xyz0[0] - xyz1[0]) < 1e-3 && TMath::Abs(xyz0[1] - xyz1[1]) < 1e-3) continue; // check wheter we are at the same global x position
    Double_t param[7];
    if(AliTracker::MeanMaterialBudget(xyz0, xyz1, param) <=0.) break;	
    Double_t xrho = param[0]*param[4]; // density*length
    Double_t xx0  = param[1]; // radiation length
    
    //Propagate the track
    track->PropagateTo(x, xx0, xrho);
    if (!AdjustSector(track)) break;
  
    //Update track
    Double_t cov[3]; ptrTracklet->GetCovAt(x, cov);
    Double_t p[2] = { ptrTracklet->GetY(), ptrTracklet->GetZ()};
    Double_t chi2 = ((AliExternalTrackParam*)track)->GetPredictedChi2(p, cov);
    if(chi2<1e+10) ((AliExternalTrackParam*)track)->Update(p, cov);
    if(!up) continue;

		//Reset material budget if 2 consecutive gold
		if(iplane>0 && track->GetTracklet(iplane-1) && ptrTracklet->GetN() + track->GetTracklet(iplane-1)->GetN() > 20) track->SetBudget(2, 0.);
	} // end planes loop

  // extrapolation
  while(ip < np){
    if(((up?-1:1) * (points[ip].GetX() - track->GetX()) < 0) && !PropagateToX(*track, points[ip].GetX(), AliTRDReconstructor::GetMaxStep())) return -1.;
    
    Double_t xyz[3]; // should also get the covariance
    track->GetXYZ(xyz); 
    track->Global2LocalPosition(xyz, track->GetAlpha());
    points[ip].SetXYZ(xyz[0], xyz[1], xyz[2]);
    ip++;
  }

	return track->GetChi2();
}

//_________________________________________________________________________
Float_t AliTRDtrackerV1::CalculateChi2Z(const AliTRDseedV1 *tracklets, Double_t offset, Double_t slope, Double_t xref)
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
    chi2Z += TMath::Abs(tracklets[iLayer].GetZfit(0) - z);
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

  // Current track X-position
  Double_t xpos = t.GetX()/*,
           mass = t.GetMass()*/;

  // Direction: inward or outward
  Double_t dir  = (xpos < xToGo) ? 1.0 : -1.0;

  while (((xToGo - xpos) * dir) > AliTRDReconstructor::GetEpsilon()) {
//    printf("to go %f\n", (xToGo - xpos) * dir);
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
    if(t.GetProlongation(x,y,z)<0) return 0; // No prolongation possible

    // The global position of the end point of this prolongation step
    xyz1[0] =  x * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha()); 
    xyz1[1] = +x * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;

    // Calculate the mean material budget between start and
    // end point of this prolongation step
    if(AliTracker::MeanMaterialBudget(xyz0, xyz1, param)<=0.) return 0;
    
    // Propagate the track to the X-position after the next step
    if (!t.PropagateTo(x, param[1], param[0]*param[4])) return 0;

/*    // Correct for mean material budget
    Double_t dEdx(0.),
             bg(t.GetP()/mass);
    if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>=3){
      const char *pn[] = {"rho", "x/X0", "<A>", "<Z>", "L", "<Z/A>", "Nb"};
      printf("D-AliTRDtrackerV1::PropagateTo(): x[%6.2f] bg[%6.2f]\n", xpos, bg);
      printf("     param :: %s[%e] %s[%e] %s[%e] %s[%e] %s[%e] %s[%e] %s[%e]\n"
          , pn[0], param[0]
          , pn[1], param[1]
          , pn[2], param[2]
          , pn[3], param[3]
          , pn[4], param[4]
          , pn[5], param[5]
          , pn[6], param[6]);
    }  
    switch(fgBB){
    case kSolid:
      dEdx = AliExternalTrackParam::BetheBlochSolid(bg);
      break;
    case kGas:
      dEdx = AliExternalTrackParam::BetheBlochGas(bg);
      break;
    case kGeant:
      { // mean exitation energy (GeV)
        Double_t mee = ((param[3] < 13.) ? (12. * param[3] + 7.) : (9.76 * param[3] + 58.8 * TMath::Power(param[3],-0.19))) * 1.e-9;
        Double_t mZA = param[5]>1.e-5?param[5]:(param[3]/param[2]);
        if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>=3) printf("D-AliTRDtrackerV1::PropagateTo(): Mee[%e] <Z/A>[%e]\n", mee, mZA);
        // protect against failed calculation of rho in MeanMaterialBudget()
        dEdx = AliExternalTrackParam::BetheBlochGeant(bg, param[0]>1.e-6?param[0]:2.33, 0.2, 3., mee, mZA);
      }
      break;
    }
    if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")>=2) printf("D-AliTRDtrackerV1::PropagateTo(): dEdx(bg=%e, m=%e)= %e[GeV/cm]\n", bg, mass, dEdx);
    if (!t.CorrectForMeanMaterialdEdx(param[1], dir*param[0]*param[4], mass, dEdx)) return 0;
*/
    // Rotate the track if necessary
    if(!AdjustSector(&t)) return 0;

    // New track X-position
    xpos = t.GetX();

  }

  return 1;

}

//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::ReadClusters(TTree *clusterTree)
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
    return kFALSE;
  }
  branch->SetAddress(&clusterArray); 
  
  if(!fClusters){ 
    Float_t nclusters =  fkRecoParam->GetNClusters();
    if(fkReconstructor->IsHLT()) nclusters /= AliTRDgeometry::kNsector;
    fClusters = new TClonesArray("AliTRDcluster", Int_t(nclusters));
    fClusters->SetOwner(kTRUE);
  }
  
  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t) clusterTree->GetEntries();
  Int_t nbytes     = 0;
  Int_t ncl        = 0;
  AliTRDcluster *c = NULL;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);  
    
    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      if(!(c = (AliTRDcluster *) clusterArray->UncheckedAt(iCluster))) continue;
      new((*fClusters)[ncl++]) AliTRDcluster(*c);
      delete (clusterArray->RemoveAt(iCluster)); 
    }
  }
  delete clusterArray;

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliTRDtrackerV1::LoadClusters(TTree *cTree)
{
  //
  // Fills clusters into TRD tracking sectors
  //
  
  fkRecoParam = fkReconstructor->GetRecoParam(); // load reco param for this event

  if(!fkReconstructor->IsWritingClusters()){ 
    fClusters = AliTRDReconstructor::GetClusters();
  } else {
    if(!ReadClusters(cTree)) {
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
Int_t AliTRDtrackerV1::LoadClusters(TClonesArray * const clusters)
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

  fkRecoParam = fkReconstructor->GetRecoParam(); // load reco param for this event
  BuildTrackingContainers();  

  //Int_t ncl  = fClusters->GetEntriesFast();
  //AliInfo(Form("Clusters %d [%6.2f %% in the active volume]", ncl, 100.*float(nin)/ncl));

  return 0;
}


//____________________________________________________________________
Int_t AliTRDtrackerV1::BuildTrackingContainers()
{
// Building tracking containers for clusters

  Int_t nin(0), ncl(fClusters->GetEntriesFast());
  while (ncl--) {
    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(ncl);
    if(c->IsInChamber()) nin++;
    if(fkReconstructor->IsHLT()) c->SetRPhiMethod(AliTRDcluster::kCOG);
    Int_t detector       = c->GetDetector();
    Int_t sector         = fGeom->GetSector(detector);
    Int_t stack          = fGeom->GetStack(detector);
    Int_t layer          = fGeom->GetLayer(detector);
    
    fTrSec[sector].GetChamber(stack, layer, kTRUE)->InsertCluster(c, ncl);
  }

  for(int isector =0; isector<AliTRDgeometry::kNsector; isector++){ 
    if(!fTrSec[isector].GetNChambers()) continue;
    fTrSec[isector].Init(fkReconstructor);
  }

  return nin;
}



//____________________________________________________________________
void AliTRDtrackerV1::UnloadClusters() 
{ 
//
// Clears the arrays of clusters and tracks. Resets sectors and timebins 
// If option "force" is also set the containers are also deleted. This is useful 
// in case of HLT

  if(fTracks){ 
    fTracks->Delete(); 
    if(HasRemoveContainers()){delete fTracks; fTracks = NULL;}
  }
  if(fTracklets){ 
    fTracklets->Delete();
    if(HasRemoveContainers()){delete fTracklets; fTracklets = NULL;}
  }
  if(fClusters){ 
    if(IsClustersOwner()) fClusters->Delete();
    
    // save clusters array in the reconstructor for further use.
    if(!fkReconstructor->IsWritingClusters()){
      AliTRDReconstructor::SetClusters(fClusters);
      SetClustersOwner(kFALSE);
    } else AliTRDReconstructor::SetClusters(NULL);
  }

  for (int i = 0; i < AliTRDgeometry::kNsector; i++) fTrSec[i].Clear();

  // Increment the Event Number
  AliTRDtrackerDebug::SetEventNumber(AliTRDtrackerDebug::GetEventNumber()  + 1);
}

// //____________________________________________________________________
// void AliTRDtrackerV1::UseClusters(const AliKalmanTrack *t, Int_t) const
// {
//   const AliTRDtrackV1 *track = dynamic_cast<const AliTRDtrackV1*>(t);
//   if(!track) return;
// 
//   AliTRDseedV1 *tracklet = NULL;
//   for(Int_t ily=AliTRDgeometry::kNlayer; ily--;){
//     if(!(tracklet = track->GetTracklet(ily))) continue;
//     AliTRDcluster *c = NULL;
//     for(Int_t ic=AliTRDseed::kNclusters; ic--;){
//       if(!(c=tracklet->GetClusters(ic))) continue;
//       c->Use();
//     }
//   }
// }
// 

//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::AdjustSector(AliTRDtrackV1 *const track) 
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
AliTRDseedV1* AliTRDtrackerV1::GetTracklet(AliTRDtrackV1 *const track, Int_t p, Int_t &idx)
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
  AliTRDseedV1 *tracklet = (idx<0) ? NULL : (AliTRDseedV1*)fTracklets->UncheckedAt(idx);

  return tracklet;
}

//____________________________________________________________________
AliTRDseedV1* AliTRDtrackerV1::SetTracklet(const AliTRDseedV1 * const tracklet)
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
AliTRDtrackV1* AliTRDtrackerV1::SetTrack(const AliTRDtrackV1 * const track)
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
  
  Int_t nTracks   = 0;
  Int_t nChambers = 0;
  AliTRDtrackingChamber **stack = NULL, *chamber = NULL;
  for(int istack = 0; istack<AliTRDgeometry::kNstack; istack++){
    if(!(stack = fTrSec[sector].GetStack(istack))) continue;
    nChambers = 0;
    for(int ilayer=0; ilayer<AliTRDgeometry::kNlayer; ilayer++){
      if(!(chamber = stack[ilayer])) continue;
      if(chamber->GetNClusters() < fgNTimeBins * fkRecoParam->GetFindableClusters()) continue;
      nChambers++;
      //AliInfo(Form("sector %d stack %d layer %d clusters %d", sector, istack, ilayer, chamber->GetNClusters()));
    }
    if(nChambers < 4) continue;
    //AliInfo(Form("Doing stack %d", istack));
    nTracks += Clusters2TracksStack(stack, fTracksESD);
  }
  if(nTracks) AliDebug(2, Form("Number of tracks: SM_%02d[%d]", sector, nTracks));

  for(int itrack=0; itrack<nTracks; itrack++){
    AliESDtrack *esdTrack((AliESDtrack*)(fTracksESD->operator[](itrack)));
    Int_t id = esd->AddTrack(esdTrack);

    // set ESD id to stand alone TRD tracks
    if (fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 0){ 
      esdTrack=esd->GetTrack(id);
      TObject *o(NULL); Int_t ic(0);
      AliTRDtrackV1 *calibTrack(NULL); 
      while((o = esdTrack->GetCalibObject(ic++))){
        if(!(calibTrack = dynamic_cast<AliTRDtrackV1*>(o))) continue;
        calibTrack->SetESDid(esdTrack->GetID());
        break;
      }
    }
  }

  // Reset Track and Candidate Number
  AliTRDtrackerDebug::SetCandidateNumber(0);
  AliTRDtrackerDebug::SetTrackNumber(0);

  // delete ESD tracks in the array
  fTracksESD->Delete();
  return nTracks;
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::Clusters2TracksStack(AliTRDtrackingChamber **stack, TClonesArray * const esdTrackList)
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

  AliTRDtrackingChamber *chamber = NULL;
  AliTRDtrackingChamber **ci = NULL;
  AliTRDseedV1 sseed[kMaxTracksStack*6]; // to be initialized
  Int_t pars[4]; // MakeSeeds parameters

  //Double_t alpha = AliTRDgeometry::GetAlpha();
  //Double_t shift = .5 * alpha;
  Int_t configs[kNConfigs];
  
  // Purge used clusters from the containers
  ci = &stack[0];
  for(Int_t ic = kNPlanes; ic--; ci++){
    if(!(*ci)) continue;
    (*ci)->Update();
  }

  // Build initial seeding configurations
  Double_t quality = BuildSeedingConfigs(stack, configs);
  if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 10){
    AliInfo(Form("Plane config %d %d %d Quality %f"
    , configs[0], configs[1], configs[2], quality));
  }

  
  // Initialize contors
  Int_t ntracks,      // number of TRD track candidates
    ntracks1,     // number of registered TRD tracks/iter
    ntracks2 = 0; // number of all registered TRD tracks in stack
  fSieveSeeding = 0;

  // Get stack index
  Int_t ic = 0; ci = &stack[0];
  while(ic<kNPlanes && !(*ci)){ic++; ci++;}
  if(!(*ci)) return ntracks2;
  Int_t istack = fGeom->GetStack((*ci)->GetDetector());

  do{
    // Loop over seeding configurations
    ntracks = 0; ntracks1 = 0;
    for (Int_t iconf = 0; iconf<fkRecoParam->GetNumberOfSeedConfigs(); iconf++) {
      pars[0] = configs[iconf];
      pars[1] = ntracks;
      pars[2] = istack;
      ntracks = MakeSeeds(stack, &sseed[6*ntracks], pars);
      //AliInfo(Form("Number of Tracks after iteration step %d: %d\n", iconf, ntracks));
      if(ntracks == kMaxTracksStack) break;
    }
    AliDebug(2, Form("Candidate TRD tracks %d in iteration %d.", ntracks, fSieveSeeding));
    if(!ntracks) break;
    
    // Sort the seeds according to their quality
    Int_t sort[kMaxTracksStack+1];
    TMath::Sort(ntracks, fTrackQuality, sort, kTRUE);
    if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1") > 2){
      AliDebug(3, "Track candidates classification:");
      for (Int_t it(0); it < ntracks; it++) {
        Int_t jt(sort[it]);
        printf("   %2d idx[%d] Quality[%e]\n", it, jt, fTrackQuality[jt]);
      }
    }
  
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
    Int_t jSieve(0), rejectedCandidates(0);
    do{
      // Check track candidates
      rejectedCandidates=0;
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
          sseed[jseed].UpdateUsed();
          if(!sseed[jseed].IsOK()) continue;
          // check if primary candidate
          if (TMath::Abs(sseed[jseed].GetYref(0) / sseed[jseed].GetX0()) < 0.158) findable++;
          ncl   += sseed[jseed].GetN();
          nused += sseed[jseed].GetNUsed();
          nlayers++;
        }

        // Filter duplicated tracks
        if (nused > 30){
          AliDebug(4, Form("REJECTED : %d idx[%d] quality[%e] tracklets[%d] usedClusters[%d]", itrack, trackIndex, fTrackQuality[trackIndex], nlayers, nused));
          fakeTrack[trackIndex] = kTRUE;
          continue;
        }
        if (ncl>0 && Float_t(nused)/ncl >= .25){
          AliDebug(4, Form("REJECTED : %d idx[%d] quality[%e] tracklets[%d] usedClusters[%d] used/ncl[%f]", itrack, trackIndex, fTrackQuality[trackIndex], nlayers, nused, Float_t(nused)/ncl));
          fakeTrack[trackIndex] = kTRUE;
          continue;
        }

        AliDebug(4, Form("Candidate[%d] Quality[%e] Tracklets[%d] Findable[%d] Ncl[%d] Nused[%d]", trackIndex, fTrackQuality[trackIndex], nlayers, findable, ncl, nused));

        // Classify tracks
        Bool_t skip = kFALSE;
        switch(jSieve){
          case 0: // select 6 tracklets primary tracks, good quality
            if(nlayers > findable || nlayers < kNPlanes) {skip = kTRUE; break;}
            if(TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -5.){skip = kTRUE; break;}
            break;

          case 1: // select shorter primary tracks, good quality
            //if(findable<4){skip = kTRUE; break;}
            if(nlayers < findable){skip = kTRUE; break;}
            if(TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -4.){skip = kTRUE; break;}
            break;

          case 2: // select 6 tracklets secondary tracks
            if(nlayers < kNPlanes) { skip = kTRUE; break;}
            if (TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -6.0){skip = kTRUE; break;}
            break;

          case 3: // select shorter tracks, good quality
            if (nlayers<4){skip = kTRUE; break;}
            if (TMath::Log(1.E-9+fTrackQuality[trackIndex]) < -5.){skip = kTRUE; break;}
            break;

          case 4: // select anything with at least 4 tracklets
            if (nlayers<4){skip = kTRUE; break;}
            //if (TMath::Log(1.E-9+fTrackQuality[trackIndex]) - nused/(nlayers-3.0) < -15.0){skip = kTRUE; break;}
            break;
        }
        if(skip){
          rejectedCandidates++;
          AliDebug(4, Form("REJECTED : %d idx[%d] quality[%e] tracklets[%d] usedClusters[%d]", itrack, trackIndex, fTrackQuality[trackIndex], nlayers, nused));
          continue;
        } else AliDebug(4, Form("ACCEPTED : %d idx[%d] quality[%e] tracklets[%d] usedClusters[%d]", itrack, trackIndex, fTrackQuality[trackIndex], nlayers, nused));

        signedTrack[trackIndex] = kTRUE;

        AliTRDseedV1 *lseed =&sseed[trackIndex*kNPlanes];
        AliTRDtrackV1 *track = MakeTrack(lseed);
        if(!track){
          AliDebug(1, "Track building failed.");
          continue;
        } else { 
          if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1") > 1){
            Int_t ich = 0; while(!(chamber = stack[ich])) ich++;
            AliDebug(2, Form("Track pt=%7.2fGeV/c SM[%2d] Done.", track->Pt(), fGeom->GetSector(chamber->GetDetector())));
          }
        }

        if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 1 && fkReconstructor->IsDebugStreaming()){
          //AliInfo(Form("Track %d [%d] nlayers %d trackQuality = %e nused %d, yref = %3.3f", itrack, trackIndex, nlayers, fTrackQuality[trackIndex], nused, trackParams[1]));

          AliTRDseedV1 *dseed[6];
          for(Int_t iseed = AliTRDgeometry::kNlayer; iseed--;) dseed[iseed] = new AliTRDseedV1(lseed[iseed]);

          //Int_t eventNrInFile = esd->GetEventNumberInFile();
          Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
          Int_t trackNumber = AliTRDtrackerDebug::GetTrackNumber();
          Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
          TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
          cstreamer << "Clusters2TracksStack"
              << "EventNumber="   << eventNumber
              << "TrackNumber="   << trackNumber
              << "CandidateNumber=" << candidateNumber
              << "Iter="        << fSieveSeeding
              << "Like="        << fTrackQuality[trackIndex]
              << "S0.="       << dseed[0]
              << "S1.="       << dseed[1]
              << "S2.="       << dseed[2]
              << "S3.="       << dseed[3]
              << "S4.="       << dseed[4]
              << "S5.="       << dseed[5]
              << "Ncl="       << ncl
              << "NLayers="   << nlayers
              << "Findable="  << findable
              << "NUsed="     << nused
              << "\n";
        }


        AliESDtrack *esdTrack = new ((*esdTrackList)[ntracks0++]) AliESDtrack();
        esdTrack->UpdateTrackParams(track, AliESDtrack::kTRDout);
        esdTrack->SetLabel(track->GetLabel());
        track->UpdateESDtrack(esdTrack);
        // write ESD-friends if neccessary
        if (fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 0){
          AliTRDtrackV1 *calibTrack = new AliTRDtrackV1(*track);
          calibTrack->SetOwner();
          esdTrack->AddCalibObject(calibTrack);
        }
        ntracks1++;
        AliTRDtrackerDebug::SetTrackNumber(AliTRDtrackerDebug::GetTrackNumber() + 1);
      }

      jSieve++;
    } while(jSieve<5 && rejectedCandidates); // end track candidates sieve
    if(!ntracks1) break;

    // increment counters
    ntracks2 += ntracks1;

    if(fkReconstructor->IsHLT()) break;
    fSieveSeeding++;

    // Rebuild plane configurations and indices taking only unused clusters into account
    quality = BuildSeedingConfigs(stack, configs);
    if(quality < 1.E-7) break; //fkReconstructor->GetRecoParam() ->GetPlaneQualityThreshold()) break;
    
    for(Int_t ip = 0; ip < kNPlanes; ip++){ 
      if(!(chamber = stack[ip])) continue;
      chamber->Build(fGeom);//Indices(fSieveSeeding);
    }

    if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 10){ 
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

  Double_t chamberQ[kNPlanes];memset(chamberQ, 0, kNPlanes*sizeof(Double_t));
  AliTRDtrackingChamber *chamber = NULL;
  for(int iplane=0; iplane<kNPlanes; iplane++){
    if(!(chamber = stack[iplane])) continue;
    chamberQ[iplane] = (chamber = stack[iplane]) ?  chamber->GetQuality() : 0.;
  }

  Double_t tconfig[kNConfigs];memset(tconfig, 0, kNConfigs*sizeof(Double_t));
  Int_t planes[] = {0, 0, 0, 0};
  for(int iconf=0; iconf<kNConfigs; iconf++){
    GetSeedingConfig(iconf, planes);
    tconfig[iconf] = fgTopologicQA[iconf];
    for(int iplane=0; iplane<4; iplane++) tconfig[iconf] *= chamberQ[planes[iplane]]; 
  }
  
  TMath::Sort((Int_t)kNConfigs, tconfig, configs, kTRUE);
  //	AliInfo(Form("q[%d] = %f", configs[0], tconfig[configs[0]]));
  // 	AliInfo(Form("q[%d] = %f", configs[1], tconfig[configs[1]]));
  // 	AliInfo(Form("q[%d] = %f", configs[2], tconfig[configs[2]]));
  
  return tconfig[configs[0]];
}

//____________________________________________________________________
Int_t AliTRDtrackerV1::MakeSeeds(AliTRDtrackingChamber **stack, AliTRDseedV1 * const sseed, const Int_t * const ipar)
{
//
// Seed tracklets and build candidate TRD tracks. The procedure is used during barrel tracking to account for tracks which are 
// either missed by TPC prolongation or conversions inside the TRD volume. 
// For stand alone tracking the procedure is used to estimate all tracks measured by TRD. 
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
// The following steps are performed:
// 1. Build seeding layers by collapsing all time bins from each of the four seeding chambers along the 
// radial coordinate. See AliTRDtrackingChamber::GetSeedingLayer() for details. The chambers selection for seeding
// is described in AliTRDtrackerV1::Clusters2TracksStack().
// 2. Using the seeding clusters from the seeding layer (step 1) build combinatorics using the following algorithm:
// - for each seeding cluster in the lower seeding layer find
// - all seeding clusters in the upper seeding layer inside a road defined by a given phi angle. The angle 
//   is calculated on the minimum pt of tracks from vertex accesible to the stand alone tracker.
// - for each pair of two extreme seeding clusters select middle upper cluster using roads defined externally by the 
//   reco params
// - select last seeding cluster as the nearest to the linear approximation of the track described by the first three
//   seeding clusters.
//   The implementation of road calculation and cluster selection can be found in the functions AliTRDchamberTimeBin::BuildCond()
//   and AliTRDchamberTimeBin::GetClusters().   
// 3. Helix fit of the seeding clusters set. (see AliTRDtrackerFitter::FitRieman(AliTRDcluster**)). No tilt correction is 
//    performed at this level 
// 4. Initialize seeding tracklets in the seeding chambers.
// 5. *Filter 0* Chi2 cut on the Y and Z directions. The threshold is set externally by the reco params.
// 6. Attach (true) clusters to seeding tracklets (see AliTRDseedV1::AttachClusters()) and fit tracklet (see 
//    AliTRDseedV1::Fit()). The number of used clusters used by current seeds should not exceed ... (25).
// 7. *Filter 1* Check if all 4 seeding tracklets are correctly constructed.
// 8. Helix fit of the clusters from the seeding tracklets with tilt correction. Refit tracklets using the new 
//    approximation of the track.
// 9. *Filter 2* Calculate likelihood of the track. (See AliTRDtrackerV1::CookLikelihood()). The following quantities are
//    checked against the Riemann fit:
//      - position resolution in y
//      - angular resolution in the bending plane
//      - likelihood of the number of clusters attached to the tracklet
// 10. Extrapolation of the helix fit to the other 2 chambers *non seeding* chambers:
//      - Initialization of extrapolation tracklets with the fit parameters
//      - Attach clusters to extrapolated tracklets
//      - Helix fit of tracklets
// 11. Improve seeding tracklets quality by reassigning clusters based on the last parameters of the track
//      See AliTRDtrackerV1::ImproveSeedQuality() for details.
// 12. Helix fit of all 6 seeding tracklets and chi2 calculation
// 13. Hyperplane fit and track quality calculation. See AliTRDtrackerFitter::FitHyperplane() for details.
// 14. Cooking labels for tracklets. Should be done only for MC
// 15. Register seeds.
//
// Authors:
//   Marian Ivanov <M.Ivanov@gsi.de>
//   Alexandru Bercuci <A.Bercuci@gsi.de>
//   Markus Fasel <M.Fasel@gsi.de>

  AliTRDtrackingChamber *chamber = NULL;
  AliTRDcluster *c[kNSeedPlanes] = {NULL, NULL, NULL, NULL}; // initilize seeding clusters
  AliTRDseedV1 *cseed = &sseed[0]; // initialize tracklets for first track
  Int_t ncl, mcl; // working variable for looping over clusters
  Int_t index[AliTRDchamberTimeBin::kMaxClustersLayer], jndex[AliTRDchamberTimeBin::kMaxClustersLayer];
  // chi2 storage
  // chi2[0] = tracklet chi2 on the Z direction
  // chi2[1] = tracklet chi2 on the R direction
  Double_t chi2[4];

  // this should be data member of AliTRDtrack TODO
  Double_t seedQuality[kMaxTracksStack];
  
  // unpack control parameters
  Int_t config  = ipar[0];
  Int_t ntracks = ipar[1];
  Int_t istack  = ipar[2];
  Int_t planes[kNSeedPlanes]; GetSeedingConfig(config, planes);	
  Int_t planesExt[kNPlanes-kNSeedPlanes]; GetExtrapolationConfig(config, planesExt);


  // Init chambers geometry
  Double_t hL[kNPlanes];       // Tilting angle
  Float_t padlength[kNPlanes]; // pad lenghts
  Float_t padwidth[kNPlanes];  // pad widths
  AliTRDpadPlane *pp = NULL;
  for(int iplane=0; iplane<kNPlanes; iplane++){
    pp                = fGeom->GetPadPlane(iplane, istack);
    hL[iplane]        = TMath::Tan(TMath::DegToRad()*pp->GetTiltingAngle());
    padlength[iplane] = pp->GetLengthIPad();
    padwidth[iplane] = pp->GetWidthIPad();
  }
  
  // Init anode wire position for chambers
  Double_t x0[kNPlanes],       // anode wire position
           driftLength = .5*AliTRDgeometry::AmThick() - AliTRDgeometry::DrThick(); // drift length
  TGeoHMatrix *matrix = NULL;
  Double_t loc[] = {AliTRDgeometry::AnodePos(), 0., 0.};
  Double_t glb[] = {0., 0., 0.};
  AliTRDtrackingChamber **cIter = &stack[0];
  for(int iLayer=0; iLayer<kNPlanes; iLayer++,cIter++){
    if(!(*cIter)) continue;
    if(!(matrix = fGeom->GetClusterMatrix((*cIter)->GetDetector()))){ 
      x0[iLayer] = fgkX0[iLayer];
      continue;
    }
    matrix->LocalToMaster(loc, glb);
    x0[iLayer] = glb[0];
  }

  AliDebug(2, Form("Making seeds Stack[%d] Config[%d] Tracks[%d]...", istack, config, ntracks));

  // Build seeding layers
  ResetSeedTB();
  Int_t nlayers = 0;
  for(int isl=0; isl<kNSeedPlanes; isl++){ 
    if(!(chamber = stack[planes[isl]])) continue;
    if(!chamber->GetSeedingLayer(fSeedTB[isl], fGeom, fkReconstructor)) continue;
    nlayers++;
  }
  if(nlayers < kNSeedPlanes) return ntracks;
  
  
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
      Double_t dzdx = (c[3]->GetZ() - c[0]->GetZ())/dx;
      Double_t dydx   = (c[3]->GetY() - c[0]->GetY())/dx;
      fSeedTB[1]->BuildCond(c[0], cond1, 1, dzdx, dydx);
      fSeedTB[1]->GetClusters(cond1, jndex, mcl);
      //printf("Found c[0] candidates 1 %d\n", mcl);

      Int_t kcl = 0;
      while(kcl<mcl) {
        c[1] = (*fSeedTB[1])[jndex[kcl++]];
        if(!c[1]) continue;
        fSeedTB[2]->BuildCond(c[1], cond2, 2, dzdx, dydx);
        c[2] = fSeedTB[2]->GetNearestCluster(cond2);
        //printf("Found c[1] candidate 2 %p\n", c[2]);
        if(!c[2]) continue;

       	AliDebug(3, Form("Seeding clusters\n 0[%6.3f %6.3f %6.3f]\n 1[%6.3f %6.3f %6.3f]\n 2[%6.3f %6.3f %6.3f]\n 3[%6.3f %6.3f %6.3f].",
          c[0]->GetX(), c[0]->GetY(), c[0]->GetZ(),
          c[1]->GetX(), c[1]->GetY(), c[1]->GetZ(),
          c[2]->GetX(), c[2]->GetY(), c[2]->GetZ(),
          c[3]->GetX(), c[3]->GetY(), c[3]->GetZ()));
              
        for (Int_t il = 0; il < kNPlanes; il++) cseed[il].Reset();
      
        FitRieman(c, chi2);
      
        AliTRDseedV1 *tseed = &cseed[0];
        cIter = &stack[0];
        for(int iLayer=0; iLayer<kNPlanes; iLayer++, tseed++, cIter++){
          Int_t det = (*cIter) ? (*cIter)->GetDetector() : -1;
          tseed->SetDetector(det);
          tseed->SetTilt(hL[iLayer]);
          tseed->SetPadLength(padlength[iLayer]);
          tseed->SetPadWidth(padwidth[iLayer]);
          tseed->SetReconstructor(fkReconstructor);
          tseed->SetX0(det<0 ? fR[iLayer]+driftLength : x0[iLayer]);
          tseed->Init(GetRiemanFitter());
          tseed->SetStandAlone(kTRUE);
        }
      
        Bool_t isFake = kFALSE;
        if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 2 && fkReconstructor->IsDebugStreaming()){
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
          TTreeSRedirector &cs0 = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
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
        if(chi2[0] > fkRecoParam->GetChi2Z()/*7./(3. - sLayer)*//*iter*/){
          AliDebug(3, Form("Filter on chi2Z [%f].", chi2[0]));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        if(chi2[1] > fkRecoParam->GetChi2Y()/*1./(3. - sLayer)*//*iter*/){
          AliDebug(3, Form("Filter on chi2Y [%f].", chi2[1]));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        //AliInfo("Passed chi2 filter.");
      
        // try attaching clusters to tracklets
        Int_t mlayers = 0; 
        AliTRDcluster *cl = NULL;
        for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
          Int_t jLayer = planes[iLayer];
          Int_t nNotInChamber = 0;
          if(!cseed[jLayer].AttachClusters(stack[jLayer], kTRUE)) continue;
          if(/*fkReconstructor->IsHLT()*/kFALSE){ 
            cseed[jLayer].UpdateUsed();
            if(!cseed[jLayer].IsOK()) continue;
          }else{
            cseed[jLayer].Fit();
            cseed[jLayer].UpdateUsed();
            cseed[jLayer].ResetClusterIter();
            while((cl = cseed[jLayer].NextCluster())){
              if(!cl->IsInChamber()) nNotInChamber++;
            }
            //printf("clusters[%d], used[%d], not in chamber[%d]\n", cseed[jLayer].GetN(), cseed[jLayer].GetNUsed(), nNotInChamber);
            if(cseed[jLayer].GetN() - (cseed[jLayer].GetNUsed() + nNotInChamber) < 5) continue; // checking for Cluster which are not in chamber is a much stronger restriction on real data
          }
          mlayers++;
        }

        if(mlayers < kNSeedPlanes){ 
          AliDebug(2, Form("Found only %d tracklets out of %d. Skip.", mlayers, kNSeedPlanes));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }

        // temporary exit door for the HLT
        if(fkReconstructor->IsHLT()){ 
          // attach clusters to extrapolation chambers
          for(int iLayer=0; iLayer<kNPlanes-kNSeedPlanes; iLayer++){
            Int_t jLayer = planesExt[iLayer];
            if(!(chamber = stack[jLayer])) continue;
            if(!cseed[jLayer].AttachClusters(chamber, kTRUE)) continue;
            cseed[jLayer].Fit();
          }
          //FitTiltedRiemanConstraint(&cseed[0], GetZ());
          fTrackQuality[ntracks] = 1.; // dummy value
          ntracks++;
          if(ntracks == kMaxTracksStack) return ntracks;
          cseed += 6; 
          continue;
        }


        // Update Seeds and calculate Likelihood
        // fit tracklets and cook likelihood
        Double_t chi2Vals[4];
        chi2Vals[0] = FitTiltedRieman(&cseed[0], kTRUE);
        for(int iLayer=0; iLayer<kNSeedPlanes; iLayer++){
          Int_t jLayer = planes[iLayer];
          cseed[jLayer].Fit(1);
        }
        Double_t like = CookLikelihood(&cseed[0], planes); // to be checked
      
        if (TMath::Log(1.E-9 + like) < fkRecoParam->GetTrackLikelihood()){
          AliDebug(3, Form("Filter on likelihood %f[%e].", TMath::Log(1.E-9 + like), like));
          AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
          continue;
        }
        //AliInfo(Form("Passed likelihood %f[%e].", TMath::Log(1.E-9 + like), like));
      
        // book preliminary results
        seedQuality[ntracks] = like;
        fSeedLayer[ntracks]  = config;/*sLayer;*/
      
        // attach clusters to the extrapolation seeds
        Int_t elayers(0);
        for(int iLayer=0; iLayer<kNPlanes-kNSeedPlanes; iLayer++){
          Int_t jLayer = planesExt[iLayer];
          if(!(chamber = stack[jLayer])) continue;
      
          // fit extrapolated seed
          if ((jLayer == 0) && !(cseed[1].IsOK())) continue;
          if ((jLayer == 5) && !(cseed[4].IsOK())) continue;
          AliTRDseedV1 pseed = cseed[jLayer];
          if(!pseed.AttachClusters(chamber, kTRUE)) continue;
          pseed.Fit(1);
          cseed[jLayer] = pseed;
          chi2Vals[0] = FitTiltedRieman(cseed,  kTRUE);
          cseed[jLayer].Fit(1);
          elayers++;
        }
      
        // AliInfo("Extrapolation done.");
        // Debug Stream containing all the 6 tracklets
        if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 2 && fkReconstructor->IsDebugStreaming()){
          TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
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
              
        if(fkRecoParam->HasImproveTracklets()){ 
          if(!ImproveSeedQuality(stack, cseed, chi2Vals[0])){
            AliTRDtrackerDebug::SetCandidateNumber(AliTRDtrackerDebug::GetCandidateNumber() + 1);
            AliDebug(3, "ImproveSeedQuality() failed.");
          }
        }
      
        // do track fitting with vertex constraint
        if(fkRecoParam->IsVertexConstrained()) chi2Vals[1] = FitTiltedRiemanConstraint(&cseed[0], GetZ());
        else chi2Vals[1] = -1.;
        chi2Vals[2] = GetChi2Z(&cseed[0]);
        chi2Vals[3] = GetChi2Phi(&cseed[0]);

        // calculate track quality
        fTrackQuality[ntracks] = CalculateTrackLikelihood(&chi2Vals[0]);
                  
        if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 2 && fkReconstructor->IsDebugStreaming()){
          TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
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
              << "Like="				<< like
              << "S0.="				<< &cseed[0]
              << "S1.="				<< &cseed[1]
              << "S2.="				<< &cseed[2]
              << "S3.="				<< &cseed[3]
              << "S4.="				<< &cseed[4]
              << "S5.="				<< &cseed[5]
              << "FitterT.="			<< fitterT
              << "FitterTC.="			<< fitterTC
              << "\n";
        }
        if(AliLog::GetDebugLevel("TRD", "AliTRDtrackerV1")){  
          Double_t pt[]={0., 0.};
          for(Int_t il(0); il<kNPlanes; il++){
            if(!cseed[il].IsOK()) continue;
            pt[0] = GetBz()*kB2C/cseed[il].GetC();
            pt[1] = GetBz()*kB2C/cseed[il].GetC(1);
            break;
          }
          AliDebug(2, Form("Candidate[%2d] pt[%7.3f %7.3f] Q[%e]\n"
            "  [0] x[%6.2f] n[%2d] nu[%d] OK[%c]\n"
            "  [1] x[%6.2f] n[%2d] nu[%d] OK[%c]\n"
            "  [2] x[%6.2f] n[%2d] nu[%d] OK[%c]\n"
            "  [3] x[%6.2f] n[%2d] nu[%d] OK[%c]\n"
            "  [4] x[%6.2f] n[%2d] nu[%d] OK[%c]\n"
            "  [5] x[%6.2f] n[%2d] nu[%d] OK[%c]"
            , ntracks, pt[0], pt[1], fTrackQuality[ntracks]
            ,cseed[0].GetX(), cseed[0].GetN(), cseed[0].GetNUsed(), cseed[0].IsOK()?'y':'n'
            ,cseed[1].GetX(), cseed[1].GetN(), cseed[1].GetNUsed(), cseed[1].IsOK()?'y':'n'
            ,cseed[2].GetX(), cseed[2].GetN(), cseed[2].GetNUsed(), cseed[2].IsOK()?'y':'n'
            ,cseed[3].GetX(), cseed[3].GetN(), cseed[3].GetNUsed(), cseed[3].IsOK()?'y':'n'
            ,cseed[4].GetX(), cseed[4].GetN(), cseed[4].GetNUsed(), cseed[4].IsOK()?'y':'n'
            ,cseed[5].GetX(), cseed[5].GetN(), cseed[5].GetNUsed(), cseed[5].IsOK()?'y':'n'));
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
AliTRDtrackV1* AliTRDtrackerV1::MakeTrack(AliTRDseedV1 * const tracklet)
{
//
// Build a TRD track out of tracklet candidates
//
// Parameters :
//   seeds  : array of tracklets
//   params : array of track parameters as they are estimated by stand alone tracker. 7 elements.
//     [0] - radial position of the track at reference point
//     [1] - y position of the fit at [0]
//     [2] - z position of the fit at [0]
//     [3] - snp of the first tracklet
//     [4] - tgl of the first tracklet
//     [5] - curvature of the Riemann fit - 1/pt
//     [6] - sector rotation angle
//
// Output :
//   The TRD track.
//
// Initialize the TRD track based on the parameters of the fit and a parametric covariance matrix 
// (diagonal with constant variance terms TODO - correct parameterization) 
// 
// In case of HLT just register the tracklets in the tracker and return values of the Riemann fit. For the
// offline case perform a full Kalman filter on the already found tracklets (see AliTRDtrackerV1::FollowBackProlongation() 
// for details). Do also MC label calculation and PID if propagation successfully.

  if(fkReconstructor->IsHLT()) FitTiltedRiemanConstraint(tracklet, 0);
  Double_t alpha = AliTRDgeometry::GetAlpha();
  Double_t shift = AliTRDgeometry::GetAlpha()/2.0;

  // find first good tracklet
  Int_t idx(0); while(idx<kNPlanes && !tracklet[idx].IsOK()) idx++;
  if(idx>2){ AliDebug(1, Form("Found suspect track start @ layer idx[%d]\n"
    "  %c[0] x0[%f] n[%d] nu[%d] OK[%c]\n"
    "  %c[1] x0[%f] n[%d] nu[%d] OK[%c]\n"
    "  %c[2] x0[%f] n[%d] nu[%d] OK[%c]\n"
    "  %c[3] x0[%f] n[%d] nu[%d] OK[%c]\n"
    "  %c[4] x0[%f] n[%d] nu[%d] OK[%c]\n"
    "  %c[5] x0[%f] n[%d] nu[%d] OK[%c]"
    ,idx
    ,idx==0?'*':' ', tracklet[0].GetX0(), tracklet[0].GetN(), tracklet[0].GetNUsed(), tracklet[0].IsOK()?'y':'n'
    ,idx==1?'*':' ', tracklet[1].GetX0(), tracklet[1].GetN(), tracklet[1].GetNUsed(), tracklet[1].IsOK()?'y':'n'
    ,idx==2?'*':' ', tracklet[2].GetX0(), tracklet[2].GetN(), tracklet[2].GetNUsed(), tracklet[2].IsOK()?'y':'n'
    ,idx==3?'*':' ', tracklet[3].GetX0(), tracklet[3].GetN(), tracklet[3].GetNUsed(), tracklet[3].IsOK()?'y':'n'
    ,idx==4?'*':' ', tracklet[4].GetX0(), tracklet[4].GetN(), tracklet[4].GetNUsed(), tracklet[4].IsOK()?'y':'n'
    ,idx==5?'*':' ', tracklet[5].GetX0(), tracklet[5].GetN(), tracklet[5].GetNUsed(), tracklet[5].IsOK()?'y':'n'));
    return NULL;
  }

  Double_t dx(5.);
  Double_t x(tracklet[idx].GetX0() - dx);
  // Build track parameters
  Double_t params[] = {
    tracklet[idx].GetYref(0) - dx*tracklet[idx].GetYref(1) // y
   ,tracklet[idx].GetZref(0) - dx*tracklet[idx].GetZref(1) // z
   ,TMath::Sin(TMath::ATan(tracklet[idx].GetYref(1)))      // snp
   ,tracklet[idx].GetZref(1) / TMath::Sqrt(1. + tracklet[idx].GetYref(1) * tracklet[idx].GetYref(1))   // tgl
   ,tracklet[idx].GetC(fkReconstructor->IsHLT()?1:0)                                   // curvature -> 1/pt
  };
  Int_t sector(fGeom->GetSector(tracklet[idx].GetDetector()));

  Double_t c[15];
  c[ 0] = 0.2; // s^2_y
  c[ 1] = 0.0; c[ 2] = 2.0; // s^2_z
  c[ 3] = 0.0; c[ 4] = 0.0; c[ 5] = 0.02; // s^2_snp
  c[ 6] = 0.0; c[ 7] = 0.0; c[ 8] = 0.0;  c[ 9] = 0.1; // s^2_tgl
  c[10] = 0.0; c[11] = 0.0; c[12] = 0.0;  c[13] = 0.0; c[14] = params[4]*params[4]*0.01; // s^2_1/pt

  AliTRDtrackV1 track(tracklet, params, c, x, sector*alpha+shift);

  AliTRDseedV1 *ptrTracklet = NULL;

  // skip Kalman filter for HLT
  if(/*fkReconstructor->IsHLT()*/kFALSE){ 
    for (Int_t jLayer = 0; jLayer < AliTRDgeometry::kNlayer; jLayer++) {
      track.UnsetTracklet(jLayer);
      ptrTracklet = &tracklet[jLayer];
      if(!ptrTracklet->IsOK()) continue;
      if(TMath::Abs(ptrTracklet->GetYref(1) - ptrTracklet->GetYfit(1)) >= .2) continue; // check this condition with Marian
      ptrTracklet = SetTracklet(ptrTracklet);
      ptrTracklet->UseClusters();
      track.SetTracklet(ptrTracklet, fTracklets->GetEntriesFast()-1);
    }
    AliTRDtrackV1 *ptrTrack = SetTrack(&track);
    ptrTrack->CookPID();
    ptrTrack->CookLabel(.9);
    ptrTrack->SetReconstructor(fkReconstructor);
    return ptrTrack;
  }

  // prevent the error message in AliTracker::MeanMaterialBudget: "start point out of geometry"
  if(TMath::Abs(track.GetX()) + TMath::Abs(track.GetY()) + TMath::Abs(track.GetZ()) > 10000) return NULL;

  track.ResetCovariance(1);
  Int_t nc = TMath::Abs(FollowBackProlongation(track));
  if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 5 && fkReconstructor->IsDebugStreaming()){
    Int_t eventNumber 		= AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    Double_t p[5]; // Track Params for the Debug Stream
    track.GetExternalParameters(x, p);
    TTreeSRedirector &cs = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
    cs << "MakeTrack"
    << "EventNumber="     << eventNumber
    << "CandidateNumber=" << candidateNumber
    << "nc="     << nc
    << "X="      << x
    << "Y="      << p[0]
    << "Z="      << p[1]
    << "snp="    << p[2]
    << "tnd="    << p[3]
    << "crv="    << p[4]
    << "Yin="    << params[0]
    << "Zin="    << params[1]
    << "snpin="  << params[2]
    << "tndin="  << params[3]
    << "crvin="  << params[4]
    << "track.=" << &track
    << "\n";
  }
  if (nc < 30){ 
    UnsetTrackletsTrack(&track);
    return NULL;
  }
  AliTRDtrackV1 *ptrTrack = SetTrack(&track);
  ptrTrack->SetReconstructor(fkReconstructor);
  ptrTrack->CookLabel(.9);
  for(Int_t il(kNPlanes); il--;){
    if(!(ptrTracklet = ptrTrack->GetTracklet(il))) continue;
    ptrTracklet->UseClusters();
  }

  // computes PID for track
  ptrTrack->CookPID();
  // update calibration references using this track
  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if(!calibra){
    AliInfo("Could not get Calibra instance.");
  } else if(calibra->GetHisto2d()){
    calibra->UpdateHistogramsV1(ptrTrack);
  }
  return ptrTrack;
}


//____________________________________________________________________
Bool_t AliTRDtrackerV1::ImproveSeedQuality(AliTRDtrackingChamber **stack, AliTRDseedV1 *cseed, Double_t &chi2)
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
  AliTRDtrackingChamber *chamber = NULL;
  AliTRDseedV1 bseed[AliTRDgeometry::kNlayer];

  Float_t quality(1.e3), 
          lQuality[AliTRDgeometry::kNlayer] = {1.e3, 1.e3, 1.e3, 1.e3, 1.e3, 1.e3};
  Int_t rLayers(0);
  for(Int_t jLayer=AliTRDgeometry::kNlayer; jLayer--;){ 
    bseed[jLayer] = cseed[jLayer];
    if(!bseed[jLayer].IsOK()) continue;
    rLayers++;
    lQuality[jLayer] = bseed[jLayer].GetQuality(kTRUE);
    quality    += lQuality[jLayer];
  }
  quality /= rLayers;
  AliDebug(2, Form("Start N[%d] Q[%f] chi2[%f]", rLayers, quality, chi2));

  for (Int_t iter = 0; iter < 4; iter++) {
    // Try better cluster set
    Int_t nLayers(0); Float_t qualitynew(0.);
    Int_t  indexes[4*AliTRDgeometry::kNlayer];
    TMath::Sort(Int_t(AliTRDgeometry::kNlayer), lQuality, indexes, kFALSE);
    for(Int_t jLayer=AliTRDgeometry::kNlayer; jLayer--;) {
      Int_t bLayer = indexes[jLayer];
      bseed[bLayer].Reset("c");
      if(!(chamber = stack[bLayer])) continue;
      if(!bseed[bLayer].AttachClusters(chamber, kTRUE)) continue;
      bseed[bLayer].Fit(1);
      if(!bseed[bLayer].IsOK()) continue;
      nLayers++;
      lQuality[jLayer] = bseed[jLayer].GetQuality(kTRUE);
      qualitynew    += lQuality[jLayer];
    }
    if(rLayers > nLayers){
      AliDebug(1, Form("Lost %d tracklets while improving.", rLayers-nLayers));
      return iter>0?kTRUE:kFALSE;
    } else rLayers=nLayers;
    qualitynew /= rLayers;

    if(qualitynew > quality){ 
      AliDebug(4, Form("Quality[%f] worsen in iter[%d] to ref[%f].", qualitynew, iter, quality));
      return iter>0?kTRUE:kFALSE;
    } else quality = qualitynew;

    // try improve track parameters
    Float_t chi2new = FitTiltedRieman(bseed, kTRUE);
    if(chi2new > chi2){ 
      AliDebug(4, Form("Chi2[%f] worsen in iter[%d] to ref[%f].", chi2new, iter, chi2));
      return iter>0?kTRUE:kFALSE;
    } else chi2 = chi2new;

    // store better tracklets
    for(Int_t jLayer=AliTRDgeometry::kNlayer; jLayer--;) cseed[jLayer]=bseed[jLayer];
    AliDebug(2, Form("Iter[%d] Q[%f] chi2[%f]", iter, quality, chi2));


    if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 7 && fkReconstructor->IsDebugStreaming()){
      Int_t eventNumber 		= AliTRDtrackerDebug::GetEventNumber();
      Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
      TLinearFitter *tiltedRieman = GetTiltedRiemanFitter();
      TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
      cstreamer << "ImproveSeedQuality"
        << "EventNumber=" 		<< eventNumber
        << "CandidateNumber="	<< candidateNumber
        << "Iteration="				<< iter
        << "S0.="							<< &cseed[0]
        << "S1.="							<< &cseed[1]
        << "S2.="							<< &cseed[2]
        << "S3.="							<< &cseed[3]
        << "S4.="							<< &cseed[4]
        << "S5.="							<< &cseed[5]
        << "FitterT.="				<< tiltedRieman
        << "\n";
    }
  } // Loop: iter

  // we are sure that at least 4 tracklets are OK !
  return kTRUE;
}

//_________________________________________________________________________
Double_t AliTRDtrackerV1::CalculateTrackLikelihood(Double_t *chi2){
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
  
  // Non-constrained Tilted Riemann
  Double_t likeChi2TR = TMath::Exp(-chi2[0] * 0.0078);
  // Constrained Tilted Riemann
  Double_t likeChi2TC(1.);
  if(chi2[1]>0.){
    likeChi2TC = TMath::Exp(-chi2[1] * 0.677);
    Double_t r = likeChi2TC/likeChi2TR;
    if(r>1.e2){;}   // -> a primary track use TC
    else if(r<1.e2) // -> a secondary track use TR
      likeChi2TC =1.;
    else{;}         // -> test not conclusive
  }
  // Chi2 only on Z direction
  Double_t likeChi2Z  = TMath::Exp(-chi2[2] * 0.14);
  // Chi2 angular resolution
  Double_t likeChi2Phi= TMath::Exp(-chi2[3] * 3.23);

  Double_t trackLikelihood     = likeChi2Z * likeChi2TR * likeChi2TC * likeChi2Phi;

  AliDebug(2, Form("Likelihood [%e]\n"
    "  Rieman : chi2[%f] likelihood[%6.2e]\n"
    "  Vertex : chi2[%f] likelihood[%6.2e]\n"
    "  Z      : chi2[%f] likelihood[%6.2e]\n"
    "  Phi    : chi2[%f] likelihood[%6.2e]"
    , trackLikelihood
    , chi2[0], likeChi2TR
    , chi2[1], likeChi2TC
    , chi2[2], likeChi2Z
    , chi2[3], likeChi2Phi
  ));

  if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 2 && fkReconstructor->IsDebugStreaming()){
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
    cstreamer << "CalculateTrackLikelihood0"
        << "EventNumber="			<< eventNumber
        << "CandidateNumber="	<< candidateNumber
        << "LikeChi2Z="				<< likeChi2Z
        << "LikeChi2TR="			<< likeChi2TR
        << "LikeChi2TC="			<< likeChi2TC
        << "LikeChi2Phi=" 		<< likeChi2Phi
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
 	Double_t chi2y = GetChi2Y(&cseed[0]);
  Double_t chi2z = GetChi2Z(&cseed[0]);

  Float_t nclusters = 0.;
  Double_t sumda = 0.;
  for(UChar_t ilayer = 0; ilayer < 4; ilayer++){
    Int_t jlayer = planes[ilayer];
    nclusters += cseed[jlayer].GetN2();
    sumda += TMath::Abs(cseed[jlayer].GetYfit(1) - cseed[jlayer].GetYref(1));
  }
  nclusters *= .25;

  Double_t likea     = TMath::Exp(-sumda * fkRecoParam->GetPhiSlope());
  Double_t likechi2y  = 0.0000000001;
  if (fkReconstructor->IsCosmic() || chi2y < fkRecoParam->GetChi2YCut()) likechi2y += TMath::Exp(-TMath::Sqrt(chi2y) * fkRecoParam->GetChi2YSlope());
  Double_t likechi2z = TMath::Exp(-chi2z * fkRecoParam->GetChi2ZSlope());
  Double_t likeN     = TMath::Exp(-(fkRecoParam->GetNMeanClusters() - nclusters) / fkRecoParam->GetNSigmaClusters());
  Double_t like      = likea * likechi2y * likechi2z * likeN;

  if(fkRecoParam->GetStreamLevel(AliTRDrecoParam::kTracker) >= 2 && fkReconstructor->IsDebugStreaming()){
    Int_t eventNumber = AliTRDtrackerDebug::GetEventNumber();
    Int_t candidateNumber = AliTRDtrackerDebug::GetCandidateNumber();
    Int_t nTracklets = 0; Float_t meanNcls = 0;
    for(Int_t iseed=0; iseed < kNPlanes; iseed++){
    	if(!cseed[iseed].IsOK()) continue;
    	nTracklets++;
    	meanNcls += cseed[iseed].GetN2();
    }
    if(nTracklets) meanNcls /= nTracklets;
    // The Debug Stream contains the seed 
    TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
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
        << "meanncls="        << meanNcls
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
  if(!fClusters) return NULL;
  Int_t ncls = fClusters->GetEntriesFast();
  return idx >= 0 && idx < ncls ? (AliCluster*)fClusters->UncheckedAt(idx) : NULL;
}

//____________________________________________________________________
AliTRDseedV1* AliTRDtrackerV1::GetTracklet(Int_t idx) const
{
  if(!fTracklets) return NULL;
  Int_t ntrklt = fTracklets->GetEntriesFast();
  return idx >= 0 && idx < ntrklt ? (AliTRDseedV1*)fTracklets->UncheckedAt(idx) : NULL;
}

//____________________________________________________________________
AliKalmanTrack* AliTRDtrackerV1::GetTrack(Int_t idx) const
{
  if(!fTracks) return NULL;
  Int_t ntrk = fTracks->GetEntriesFast();
  return idx >= 0 && idx < ntrk ? (AliKalmanTrack*)fTracks->UncheckedAt(idx) : NULL;
}



// //_____________________________________________________________________________
// Int_t AliTRDtrackerV1::Freq(Int_t n, const Int_t *inlist
//           , Int_t *outlist, Bool_t down)
// {    
//   //
//   // Sort eleements according occurancy 
//   // The size of output array has is 2*n 
//   //
// 
//   if (n <= 0) {
//     return 0;
//   }
// 
//   Int_t *sindexS = new Int_t[n];   // Temporary array for sorting
//   Int_t *sindexF = new Int_t[2*n];   
//   for (Int_t i = 0; i < n; i++) {
//     sindexF[i] = 0;
//   }
// 
//   TMath::Sort(n,inlist,sindexS,down); 
// 
//   Int_t last     = inlist[sindexS[0]];
//   Int_t val      = last;
//   sindexF[0]     = 1;
//   sindexF[0+n]   = last;
//   Int_t countPos = 0;
// 
//   // Find frequency
//   for (Int_t i = 1; i < n; i++) {
//     val = inlist[sindexS[i]];
//     if (last == val) {
//       sindexF[countPos]++;
//     }
//     else {      
//       countPos++;
//       sindexF[countPos+n] = val;
//       sindexF[countPos]++;
//       last                = val;
//     }
//   }
//   if (last == val) {
//     countPos++;
//   }
// 
//   // Sort according frequency
//   TMath::Sort(countPos,sindexF,sindexS,kTRUE);
// 
//   for (Int_t i = 0; i < countPos; i++) {
//     outlist[2*i  ] = sindexF[sindexS[i]+n];
//     outlist[2*i+1] = sindexF[sindexS[i]];
//   }
// 
//   delete [] sindexS;
//   delete [] sindexF;
//   
//   return countPos;
// 
// }


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
Float_t AliTRDtrackerV1::GetChi2Y(const AliTRDseedV1 * const tracklets) const
{
  //	Calculates normalized chi2 in y-direction
  // chi2 = Sum chi2 / n_tracklets

  Double_t chi2 = 0.; Int_t n = 0;
  for(Int_t ipl = kNPlanes; ipl--;){
    if(!tracklets[ipl].IsOK()) continue;
    chi2 += tracklets[ipl].GetChi2Y();
    n++;
  }
  return n ? chi2/n : 0.;
}

//_____________________________________________________________________________
Float_t AliTRDtrackerV1::GetChi2Z(const AliTRDseedV1 *const tracklets) const 
{
  //	Calculates normalized chi2 in z-direction
  // chi2 = Sum chi2 / n_tracklets

  Double_t chi2 = 0; Int_t n = 0;
  for(Int_t ipl = kNPlanes; ipl--;){
    if(!tracklets[ipl].IsOK()) continue;
    chi2 += tracklets[ipl].GetChi2Z();
    n++;
  }
  return n ? chi2/n : 0.;
}

//_____________________________________________________________________________
Float_t AliTRDtrackerV1::GetChi2Phi(const AliTRDseedV1 *const tracklets) const 
{
  //  Calculates normalized chi2 for angular resolution
  // chi2 = Sum chi2 / n_tracklets

  Double_t chi2 = 0; Int_t n = 0;
  for (Int_t iLayer = 0; iLayer < kNPlanes; iLayer++) {
    if(!tracklets[iLayer].IsOK()) continue;
    chi2 += tracklets[iLayer].GetChi2Phi();
    n++;
  }
  return n ? chi2/n: 0.;
}

//____________________________________________________________________
Float_t AliTRDtrackerV1::CalculateReferenceX(const AliTRDseedV1 *const tracklets){
 	//
 	// Calculates the reference x-position for the tilted Rieman fit defined as middle
 	// of the stack (middle between layers 2 and 3). For the calculation all the tracklets
 	// are taken into account
 	//
 	// Parameters: - Array of tracklets(AliTRDseedV1)
 	//
 	// Output: - The reference x-position(Float_t)
  // Only kept for compatibility with the old code
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
 	    if(iok) idiff++; // to get the right difference;
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
Double_t AliTRDtrackerV1::FitTiltedRiemanV1(AliTRDseedV1 *const tracklets){
  //
  // Track Fitter Function using the new class implementation of 
  // the Rieman fit
  //
  AliTRDtrackFitterRieman fitter;
  fitter.SetRiemanFitter(GetTiltedRiemanFitter());
  fitter.Reset();
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++) fitter.SetTracklet(il, &tracklets[il]);
  Double_t chi2 = fitter.Eval();
  // Update the tracklets
  Double_t cov[15]; Double_t x0;
  memset(cov, 0, sizeof(Double_t) * 15);
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
    x0 = tracklets[il].GetX0();
    tracklets[il].SetYref(0, fitter.GetYat(x0));
    tracklets[il].SetZref(0, fitter.GetZat(x0));
    tracklets[il].SetYref(1, fitter.GetDyDxAt(x0));
    tracklets[il].SetZref(1, fitter.GetDzDx());
    tracklets[il].SetC(fitter.GetCurvature());
    fitter.GetCovAt(x0, cov);
    tracklets[il].SetCovRef(cov);
    tracklets[il].SetChi2(chi2);
  }
  return chi2;
}

//____________________________________________________________________
void AliTRDtrackerV1::UnsetTrackletsTrack(const AliTRDtrackV1 * const track)
{
//  Remove tracklets from tracker list attached to "track"
  Int_t idx(-1);
  for(Int_t il(0); il<kNPlanes; il++){
    if((idx = track->GetTrackletIndex(il)) < 0) continue;
    delete (fTracklets->RemoveAt(idx));
  }
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
// Fast solving linear regresion in 2D
//         y=a + bx
// The data members have the following meaning
// fParams[0] : a
// fParams[1] : b
// 
// fSums[0] : S
// fSums[1] : Sx
// fSums[2] : Sy
// fSums[3] : Sxy
// fSums[4] : Sxx
// fSums[5] : Syy
// 
// fCovarianceMatrix[0] : s2a
// fCovarianceMatrix[1] : s2b
// fCovarianceMatrix[2] : cov(ab)

  memset(fParams, 0, sizeof(Double_t) * 2);
  memset(fSums, 0, sizeof(Double_t) * 6);
  memset(fCovarianceMatrix, 0, sizeof(Double_t) * 3);

}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::AddPoint(const Double_t *const x, Double_t y, Double_t sigmaY){
  //
  // Adding Point to the fitter
  //
  
  Double_t weight = 1/(sigmaY > 1e-9 ? sigmaY : 1e-9);
  weight *= weight;
  const Double_t &xpt = *x;
  //	printf("Adding point x = %f, y = %f, sigma = %f\n", xpt, y, sigmaY);
  fSums[0] += weight;
  fSums[1] += weight * xpt;
  fSums[2] += weight * y;
  fSums[3] += weight * xpt * y;
  fSums[4] += weight * xpt * xpt;
  fSums[5] += weight * y * y;
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::RemovePoint(const Double_t *const x, Double_t y, Double_t sigmaY){
  //
  // Remove Point from the sample
  //

  Double_t weight = 1/(sigmaY > 1e-9 ? sigmaY : 1e-9);
  weight *= weight;
  const Double_t &xpt = *x; 
  fSums[0] -= weight;
  fSums[1] -= weight * xpt;
  fSums[2] -= weight * y;
  fSums[3] -= weight * xpt * y;
  fSums[4] -= weight * xpt * xpt;
  fSums[5] -= weight * y * y;
}

//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::AliTRDLeastSquare::Eval(){
  //
  // Evaluation of the fit:
  // Calculation of the parameters
  // Calculation of the covariance matrix
  //
  
  Double_t det = fSums[0] * fSums[4] - fSums[1] *fSums[1];
  if(TMath::Abs(det)<1.e-30) return kFALSE;

  //	for(Int_t isum = 0; isum < 5; isum++)
  //		printf("fSums[%d] = %f\n", isum, fSums[isum]);
  //	printf("denominator = %f\n", denominator);
  fParams[0] = (fSums[2] * fSums[4] - fSums[1] * fSums[3])/det;
  fParams[1] = (fSums[0] * fSums[3] - fSums[1] * fSums[2])/det;
  //	printf("fParams[0] = %f, fParams[1] = %f\n", fParams[0], fParams[1]);
  
  // Covariance matrix
  Double_t den = fSums[0]*fSums[4] - fSums[1]*fSums[1];
  fCovarianceMatrix[0] = fSums[4] / den;
  fCovarianceMatrix[1] = fSums[0] / den;
  fCovarianceMatrix[2] = -fSums[1] / den;
/*  fCovarianceMatrix[0] = fSums[4] / fSums[0] - fSums[1] * fSums[1] / (fSums[0] * fSums[0]);
  fCovarianceMatrix[1] = fSums[5] / fSums[0] - fSums[2] * fSums[2] / (fSums[0] * fSums[0]);
  fCovarianceMatrix[2] = fSums[3] / fSums[0] - fSums[1] * fSums[2] / (fSums[0] * fSums[0]);*/



  return kTRUE;
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDLeastSquare::GetFunctionValue(const Double_t *const xpos) const {
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

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDLeastSquare::Reset(){
  //
  // Reset the fitter
  //
  memset(fParams, 0, sizeof(Double_t) * 2);
  memset(fCovarianceMatrix, 0, sizeof(Double_t) * 3);
  memset(fSums, 0, sizeof(Double_t) * 6);
}

///////////////////////////////////////////////////////
//                                                   //
// Resources of class AliTRDtrackFitterRieman        //
//                                                   //
///////////////////////////////////////////////////////

//_____________________________________________________________________________
AliTRDtrackerV1::AliTRDtrackFitterRieman::AliTRDtrackFitterRieman():
  fTrackFitter(NULL),
  fZfitter(NULL),
  fCovarPolY(NULL),
  fCovarPolZ(NULL),
  fXref(0.),
  fSysClusterError(0.)
{
  //
  // Default constructor
  //
  fZfitter = new AliTRDLeastSquare;
  fCovarPolY = new TMatrixD(3,3);
  fCovarPolZ = new TMatrixD(2,2);
  memset(fTracklets, 0, sizeof(AliTRDseedV1 *) * 6);
  memset(fParameters, 0, sizeof(Double_t) * 5);
  memset(fSumPolY, 0, sizeof(Double_t) * 5);
  memset(fSumPolZ, 0, sizeof(Double_t) * 2);
}

//_____________________________________________________________________________
AliTRDtrackerV1::AliTRDtrackFitterRieman::~AliTRDtrackFitterRieman(){
  //
  // Destructor
  //
  if(fZfitter) delete fZfitter;
  if(fCovarPolY) delete fCovarPolY;
  if(fCovarPolZ) delete fCovarPolZ;
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDtrackFitterRieman::Reset(){
  //
  // Reset the Fitter
  //
  if(fTrackFitter){
    fTrackFitter->StoreData(kTRUE);
    fTrackFitter->ClearPoints();
  }
  if(fZfitter){
    fZfitter->Reset();
  }
  fXref = 0.;
  memset(fTracklets, 0, sizeof(AliTRDseedV1 *) * AliTRDgeometry::kNlayer);
  memset(fParameters, 0, sizeof(Double_t) * 5);
  memset(fSumPolY, 0, sizeof(Double_t) * 5);
  memset(fSumPolZ, 0, sizeof(Double_t) * 2);
  for(Int_t irow = 0; irow < fCovarPolY->GetNrows(); irow++)
    for(Int_t icol = 0; icol < fCovarPolY->GetNcols(); icol++){
      (*fCovarPolY)(irow, icol) = 0.;
      if(irow < 2 && icol < 2)
        (*fCovarPolZ)(irow, icol) = 0.;
    }
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDtrackFitterRieman::SetTracklet(Int_t itr, AliTRDseedV1 *tracklet){ 
  //
  // Add tracklet into the fitter
  //
  if(itr >= AliTRDgeometry::kNlayer) return;
  fTracklets[itr] = tracklet; 
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDtrackFitterRieman::Eval(){
  //
  // Perform the fit
  // 1. Apply linear transformation and store points in the fitter
  // 2. Evaluate the fit
  // 3. Check if the result of the fit in z-direction is reasonable
  // if not
  // 3a. Fix the parameters 3 and 4 with the results of a simple least
  //     square fit
  // 3b. Redo the fit with the fixed parameters
  // 4. Store fit results (parameters and errors)
  //
  if(!fTrackFitter){
    return 1e10;
  }
  fXref = CalculateReferenceX();
  for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++) UpdateFitters(fTracklets[il]);
  if(!fTrackFitter->GetNpoints()) return 1e10;
  // perform the fit
  fTrackFitter->Eval();
  fZfitter->Eval();
  fParameters[3] = fTrackFitter->GetParameter(3);
  fParameters[4] = fTrackFitter->GetParameter(4);
  if(!CheckAcceptable(fParameters[3], fParameters[4])) {
    fTrackFitter->FixParameter(3, fZfitter->GetFunctionValue(&fXref));
    fTrackFitter->FixParameter(4, fZfitter->GetFunctionParameter(1));
    fTrackFitter->Eval();
    fTrackFitter->ReleaseParameter(3);
    fTrackFitter->ReleaseParameter(4);
    fParameters[3] = fTrackFitter->GetParameter(3);
    fParameters[4] = fTrackFitter->GetParameter(4);
  }
  // Update the Fit Parameters and the errors
  fParameters[0] = fTrackFitter->GetParameter(0);
  fParameters[1] = fTrackFitter->GetParameter(1);
  fParameters[2] = fTrackFitter->GetParameter(2);

  // Prepare Covariance estimation
  (*fCovarPolY)(0,0) = fSumPolY[0]; (*fCovarPolY)(1,1) = fSumPolY[2]; (*fCovarPolY)(2,2) = fSumPolY[4];
  (*fCovarPolY)(1,0) = (*fCovarPolY)(0,1) = fSumPolY[1];
  (*fCovarPolY)(2,0) = (*fCovarPolY)(0,2) = fSumPolY[2];
  (*fCovarPolY)(2,1) = (*fCovarPolY)(1,2) = fSumPolY[3];
  fCovarPolY->Invert();
  (*fCovarPolZ)(0,0) = fSumPolZ[0]; (*fCovarPolZ)(1,1) = fSumPolZ[2];
  (*fCovarPolZ)(1,0) = (*fCovarPolZ)(0,1) = fSumPolZ[1];
  fCovarPolZ->Invert();
  return fTrackFitter->GetChisquare() / fTrackFitter->GetNpoints();
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDtrackFitterRieman::UpdateFitters(AliTRDseedV1 * const tracklet){
  //
  // Does the transformations and updates the fitters
  // The following transformation is applied
  //
  AliTRDcluster *cl = NULL;
  Double_t x, y, z, dx, t, w, we, yerr, zerr;
  Double_t uvt[4];
  if(!tracklet || !tracklet->IsOK()) return; 
  Double_t tilt = tracklet->GetTilt();
  for(Int_t itb = 0; itb < AliTRDseedV1::kNclusters; itb++){
    if(!(cl = tracklet->GetClusters(itb))) continue;
    if(!cl->IsInChamber()) continue;
    if (!tracklet->IsUsable(itb)) continue;
    x = cl->GetX();
    y = cl->GetY();
    z = cl->GetZ();
    dx = x - fXref;
    // Transformation
    t = 1./(x*x + y*y);
    uvt[0] = 2. * x * t;
    uvt[1] = t;
    uvt[2] = 2. * tilt * t;
    uvt[3] = 2. * tilt * dx * t;
    w = 2. * (y + tilt*z) * t;
    // error definition changes for the different calls
    we = 2. * t;
    we *= TMath::Sqrt(cl->GetSigmaY2()+tilt*tilt*cl->GetSigmaZ2());
    // Update sums for error calculation
    yerr = 1./(TMath::Sqrt(cl->GetSigmaY2()) + fSysClusterError);
    yerr *= yerr;
    zerr = 1./cl->GetSigmaZ2();
    for(Int_t ipol = 0; ipol < 5; ipol++){
      fSumPolY[ipol] += yerr;
      yerr *= x;
      if(ipol < 3){
        fSumPolZ[ipol] += zerr;
        zerr *= x;
      }
    }
    fTrackFitter->AddPoint(uvt, w, we);
    fZfitter->AddPoint(&x, z, static_cast<Double_t>(TMath::Sqrt(cl->GetSigmaZ2())));
  }
}

//_____________________________________________________________________________
Bool_t AliTRDtrackerV1::AliTRDtrackFitterRieman::CheckAcceptable(Double_t offset, Double_t slope){
  // 
  // Check whether z-results are acceptable
  // Definition: Distance between tracklet fit and track fit has to be
  // less then half a padlength
  // Point of comparision is at the anode wire
  //
  Bool_t acceptablez = kTRUE;
  Double_t zref = 0.0;
  for (Int_t iLayer = 0; iLayer < kNPlanes; iLayer++) {
    if(!fTracklets[iLayer]->IsOK()) continue;
    zref = offset + slope * (fTracklets[iLayer]->GetX0() - fXref);
    if (TMath::Abs(fTracklets[iLayer]->GetZfit(0) - zref) > fTracklets[iLayer]->GetPadLength() * 0.5 + 1.0) 
      acceptablez = kFALSE;
  }
  return acceptablez;
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDtrackFitterRieman::GetYat(Double_t x) const {
  //
  // Calculate y position out of the track parameters
  // y:     R^2 = (x - x0)^2 + (y - y0)^2
  //     =>   y = y0 +/- Sqrt(R^2 - (x - x0)^2)
  //          R = Sqrt() = 1/Curvature
  //     =>   y = y0 +/- Sqrt(1/Curvature^2 - (x - x0)^2)
  //
  Double_t y = 0;
  Double_t disc = (x * fParameters[0] + fParameters[1]);
  disc = 1 - fParameters[0]*fParameters[2] + fParameters[1]*fParameters[1] - disc*disc;
  if (disc >= 0) {
    disc = TMath::Sqrt(disc);
    y    = (1.0 - disc) / fParameters[0];
  }
  return y;
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDtrackFitterRieman::GetZat(Double_t x) const {
  //
  // Return z position for a given x position
  // Simple linear function
  //
  return fParameters[3] + fParameters[4] * (x - fXref);
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDtrackFitterRieman::GetDyDxAt(Double_t x) const {
  //
  // Calculate dydx at a given radial position out of the track parameters
  // dy:      R^2 = (x - x0)^2 + (y - y0)^2
  //     =>     y = +/- Sqrt(R^2 - (x - x0)^2) + y0
  //     => dy/dx = (x - x0)/Sqrt(R^2 - (x - x0)^2) 
  // Curvature: cr = 1/R = a/Sqrt(1 + b^2 - c*a)
  //     => dy/dx =  (x - x0)/(1/(cr^2) - (x - x0)^2) 
  //
  Double_t x0 = -fParameters[1] / fParameters[0];
  Double_t curvature = GetCurvature();
  Double_t dy = 0;
  if (-fParameters[2] * fParameters[0] + fParameters[1] * fParameters[1] + 1 > 0) {
    if (1.0/(curvature * curvature) - (x - x0) * (x - x0) > 0.0) {
     Double_t yderiv = (x - x0) / TMath::Sqrt(1.0/(curvature * curvature) - (x - x0) * (x - x0));
      if (fParameters[0] < 0) yderiv *= -1.0;
      dy = yderiv;
    }
  }
  return dy;
}

//_____________________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDtrackFitterRieman::GetCurvature() const {
  //
  // Calculate track curvature
  //
  //
  Double_t curvature =  1.0 + fParameters[1]*fParameters[1] - fParameters[2]*fParameters[0];
  if (curvature > 0.0) 
    curvature  =  fParameters[0] / TMath::Sqrt(curvature);
  return curvature;
}

//_____________________________________________________________________________
void AliTRDtrackerV1::AliTRDtrackFitterRieman::GetCovAt(Double_t x, Double_t *cov) const {
  //
  // Error Definition according to gauss error propagation
  //  
  TMatrixD transform(3,3);
  transform(0,0) = transform(1,1) = transform(2,2) = 1;
  transform(0,1) = transform(1,2) = x;
  transform(0,2) = x*x;
  TMatrixD covariance(transform, TMatrixD::kMult, *fCovarPolY);
  covariance *= transform.T();
  cov[0] = covariance(0,0);
  TMatrixD transformZ(2,2);
  transformZ(0,0) = transformZ(1,1) = 1;
  transformZ(0,1) = x;
  TMatrixD covarZ(transformZ, TMatrixD::kMult, *fCovarPolZ);
  covarZ *= transformZ.T();
  cov[1] = covarZ(0,0);
  cov[2] = 0;
}

//____________________________________________________________________
Double_t AliTRDtrackerV1::AliTRDtrackFitterRieman::CalculateReferenceX(){
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
    if(fTracklets[il]->IsOK() && fTracklets[il -1]->IsOK()){
      Float_t xdiff = fTracklets[il]->GetX0() - fTracklets[il -1]->GetX0();
      meanDistance += xdiff;
      nDistances++;
    }
    if(fTracklets[il]->IsOK()) startIndex = il;
  }
  if(fTracklets[0]->IsOK()) startIndex = 0;
  if(!nDistances){
    // We should normally never get here
    Float_t xpos[2]; memset(xpos, 0, sizeof(Float_t) * 2);
    Int_t iok = 0, idiff = 0;
    // This attempt is worse and should be avoided:
    // check for two chambers which are OK and repeat this without taking the mean value
    // Strategy avoids a division by 0;
    for(Int_t il = 5; il >= 0; il--){
      if(fTracklets[il]->IsOK()){
        xpos[iok] = fTracklets[il]->GetX0();
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
  return fTracklets[startIndex]->GetX0() + (2.5 - startIndex) * meanDistance - 0.5 * (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());
}
