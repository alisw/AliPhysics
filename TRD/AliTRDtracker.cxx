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
//  The standard TRD tracker                                                 //  
//  Based on Kalman filltering approach                                      //
//                                                                           //
//  Authors:                                                                 //
//     M. Ivanov  (Marian.Ivanov@cern.ch)                                    //
//     Y. Belikov (Jouri.Belikov@cern.ch)                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TBranch.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLinearFitter.h>
#include <TObjArray.h> 
#include <TROOT.h>
#include <TTree.h>  
#include <TTreeStream.h>

#include "AliESD.h"
#include "AliAlignObj.h"
#include "AliRieman.h"
#include "AliTrackPointArray.h"

#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliTRDseed.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDtracker.h"
#include "AliTRDReconstructor.h"
#include "AliTRDCalibra.h"
ClassImp(AliTRDtracker)

const  Float_t  AliTRDtracker::fgkMinClustersInTrack =  0.5;  
const  Float_t  AliTRDtracker::fgkLabelFraction      =  0.8;  // ??
const  Double_t AliTRDtracker::fgkMaxChi2            = 12.0; 
const  Double_t AliTRDtracker::fgkMaxSnp             =  0.95; // Corresponds to tan = 3
const  Double_t AliTRDtracker::fgkMaxStep            =  2.0;  // Maximal step size in propagation 

//_____________________________________________________________________________
AliTRDtracker::AliTRDtracker()
  :AliTracker()
  ,fHBackfit(0x0)
  ,fHClSearch(0x0)
  ,fHRefit(0x0)
  ,fHX(0x0)
  ,fHNCl(0x0)
  ,fHNClTrack(0x0)
  ,fHMinYPos(0x0)
  ,fHMinYNeg(0x0)
  ,fHMinZ(0x0)
  ,fHMinD(0x0)
  ,fHDeltaX(0x0)
  ,fHXCl(0x0)
  ,fGeom(0)
  ,fNclusters(0)
  ,fClusters(0)
  ,fNseeds(0)
  ,fSeeds(0)
  ,fNtracks(0)
  ,fTracks(0)
  ,fTimeBinsPerPlane(0)
  ,fAddTRDseeds(kFALSE)
  ,fNoTilt(kFALSE)
  ,fDebugStreamer(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < kTrackingSectors; i++) {
    fTrSec[i] = 0;
  }
  for (Int_t j = 0; j <  5; j++) {
    for (Int_t k = 0; k < 18; k++) {
      fHoles[j][k] = kFALSE;
    }
  }

  InitLogHists();

} 

//_____________________________________________________________________________
AliTRDtracker::AliTRDtracker(const AliTRDtracker &t)
  :AliTracker(t)
  ,fHBackfit(0x0)
  ,fHClSearch(0x0)
  ,fHRefit(0x0)
  ,fHX(0x0)
  ,fHNCl(0x0)
  ,fHNClTrack(0x0)
  ,fHMinYPos(0x0)
  ,fHMinYNeg(0x0)
  ,fHMinZ(0x0)
  ,fHMinD(0x0)
  ,fHDeltaX(0x0)
  ,fHXCl(0x0)
  ,fGeom(0)
  ,fNclusters(0)
  ,fClusters(0)
  ,fNseeds(0)
  ,fSeeds(0)
  ,fNtracks(0)
  ,fTracks(0)
  ,fTimeBinsPerPlane(0)
  ,fAddTRDseeds(kFALSE)
  ,fNoTilt(kFALSE)
  ,fDebugStreamer(0)
{
  //
  // Copy constructor
  //

}

//_____________________________________________________________________________
AliTRDtracker::AliTRDtracker(const TFile *geomfile)
  :AliTracker()
  ,fHBackfit(0x0)
  ,fHClSearch(0x0)
  ,fHRefit(0x0)
  ,fHX(0x0)
  ,fHNCl(0x0)
  ,fHNClTrack(0x0)
  ,fHMinYPos(0x0)
  ,fHMinYNeg(0x0)
  ,fHMinZ(0x0)
  ,fHMinD(0x0)
  ,fHDeltaX(0x0)
  ,fHXCl(0x0)
  ,fGeom(0)
  ,fNclusters(0)
  ,fClusters(new TObjArray(2000))
  ,fNseeds(0)
  ,fSeeds(new TObjArray(2000))
  ,fNtracks(0)
  ,fTracks(new TObjArray(1000))
  ,fTimeBinsPerPlane(0)
  ,fAddTRDseeds(kFALSE)
  ,fNoTilt(kFALSE)
  ,fDebugStreamer(0)
{
  // 
  //  Main constructor
  //  
   
  TDirectory *savedir = gDirectory; 
  TFile      *in      = (TFile *) geomfile; 
 
  if (!in->IsOpen()) {
    AliWarning("geometry file is not open!\n");
    AliWarning("FULL TRD geometry and DEFAULT TRD parameter will be used\n");
  }
  else {
    in->cd();  
    fGeom = (AliTRDgeometry *) in->Get("TRDgeometry");
  }

  if (!fGeom) {
    AliWarning("Cannot find TRD geometry!\n");
    fGeom = new AliTRDgeometry();
  } 
  fGeom->ReadGeoMatrices();

  savedir->cd();  

  for (Int_t geomS = 0; geomS < kTrackingSectors; geomS++) {
    Int_t trS   = CookSectorIndex(geomS);
    fTrSec[trS] = new AliTRDtrackingSector(fGeom,geomS);
    for (Int_t icham = 0; icham < AliTRDgeometry::kNcham; icham++) {
      fHoles[icham][trS] = fGeom->IsHole(0,icham,geomS);
    }
  }

  AliTRDpadPlane *padPlane = AliTRDCommonParam::Instance()->GetPadPlane(0,0);
  Float_t tiltAngle = TMath::Abs(padPlane->GetTiltingAngle());
  if (tiltAngle < 0.1) {
    fNoTilt = kTRUE;
  }

  fTimeBinsPerPlane =  AliTRDcalibDB::Instance()->GetNumberOfTimeBins();

  fDebugStreamer    = new TTreeSRedirector("TRDdebug.root");

  savedir->cd();
  
  InitLogHists();

}   

//_____________________________________________________________________________
AliTRDtracker::~AliTRDtracker()
{
  //
  // Destructor of AliTRDtracker 
  //

  if (fClusters) {
    fClusters->Delete();
    delete fClusters;
  }

  if (fTracks) {
    fTracks->Delete();
    delete fTracks;
  }

  if (fSeeds) {
    fSeeds->Delete();
    delete fSeeds;
  }

  delete fGeom;  

  for (Int_t geomS = 0; geomS < kTrackingSectors; geomS++) {
    delete fTrSec[geomS];
  }

  if (fDebugStreamer) {    
    delete fDebugStreamer;
  }

}   

//_____________________________________________________________________________
Int_t  AliTRDtracker::LocalToGlobalID(Int_t lid)
{
  //
  // Transform internal TRD ID to global detector ID
  //

  Int_t  isector  = fGeom->GetSector(lid);
  Int_t  ichamber = fGeom->GetChamber(lid);
  Int_t  iplan    = fGeom->GetPlane(lid);

  AliAlignObj::ELayerID iLayer = AliAlignObj::kTRD1;
  switch (iplan) {
  case 0:
    iLayer = AliAlignObj::kTRD1;
    break;
  case 1:
    iLayer = AliAlignObj::kTRD2;
    break;
  case 2:
    iLayer = AliAlignObj::kTRD3;
    break;
  case 3:
    iLayer = AliAlignObj::kTRD4;
    break;
  case 4:
    iLayer = AliAlignObj::kTRD5;
    break;
  case 5:
    iLayer = AliAlignObj::kTRD6;
    break;
  };

  Int_t    modId = isector * fGeom->Ncham() + ichamber;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,modId);

  return volid;

}

//_____________________________________________________________________________
Int_t  AliTRDtracker::GlobalToLocalID(Int_t gid)
{
  //
  // Transform global detector ID to local detector ID
  // 

  Int_t modId    = 0;
  AliAlignObj::ELayerID layerId = AliAlignObj::VolUIDToLayer(gid,modId);

  Int_t isector  = modId / fGeom->Ncham();
  Int_t ichamber = modId % fGeom->Ncham();
  Int_t iLayer   = -1;

  switch (layerId) {
  case AliAlignObj::kTRD1:
    iLayer = 0;
    break;
  case AliAlignObj::kTRD2:
    iLayer = 1;
    break;
  case AliAlignObj::kTRD3:
    iLayer = 2;
    break;
  case AliAlignObj::kTRD4:
    iLayer = 3;
    break;
  case AliAlignObj::kTRD5:
    iLayer = 4;
    break;
  case AliAlignObj::kTRD6:
    iLayer = 5;
    break;
  default:
    iLayer =-1;
  }

  if (iLayer < 0) {
    return -1;
  }

  Int_t lid = fGeom->GetDetector(iLayer,ichamber,isector);

  return lid;

}

//_____________________________________________________________________________
Bool_t  AliTRDtracker::Transform(AliTRDcluster *cluster)
{
  //
  // Transform something ... whatever ...
  //

  // Magic constants for geo manager transformation
  const Double_t kX0shift  = 2.52;
  const Double_t kX0shift5 = 3.05;

  //
  // Apply alignment and calibration to transform cluster
  //
  Int_t detector = cluster->GetDetector();
  Int_t plane    = fGeom->GetPlane(cluster->GetDetector());
  Int_t chamber  = fGeom->GetChamber(cluster->GetDetector());
  Int_t sector   = fGeom->GetSector(cluster->GetDetector());

  Double_t dxAmp  = (Double_t) fGeom->CamHght();               // Amplification region
  Double_t driftX = TMath::Max(cluster->GetX()-dxAmp*0.5,0.0); // Drift distance

  //
  // ExB correction
  //
  Double_t vdrift = AliTRDcalibDB::Instance()->GetVdrift(cluster->GetDetector(),0,0);
  Double_t exB    = AliTRDcalibDB::Instance()->GetOmegaTau(vdrift,-AliTracker::GetBz()*0.1);

  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();  
  AliTRDpadPlane    *padPlane    = commonParam->GetPadPlane(plane,chamber);
  Double_t zshiftIdeal = 0.5*(padPlane->GetRow0()+padPlane->GetRowEnd());
  Double_t localPos[3];
  Double_t localPosTracker[3];
  localPos[0] = -cluster->GetX();
  localPos[1] =  cluster->GetY() - driftX * exB;
  localPos[2] =  cluster->GetZ() - zshiftIdeal;

  cluster->SetY(cluster->GetY() - driftX*exB);
  Double_t xplane = (Double_t) AliTRDgeometry::GetTime0(plane); 
  cluster->SetX(xplane- cluster->GetX());

  TGeoHMatrix *matrix =  fGeom->GetCorrectionMatrix(cluster->GetDetector());
  if (!matrix) {
    // No matrix found - if somebody used geometry with holes
    AliError("Invalid Geometry - Default Geometry used\n");
    return kTRUE;   
  }
  matrix->LocalToMaster(localPos,localPosTracker);  

  if (AliTRDReconstructor::StreamLevel() > 1) {
    (* fDebugStreamer) << "Transform"
                       << "Cl.="      << cluster
                       << "matrix.="  << matrix
                       << "Detector=" << detector
                       << "Sector="   << sector
                       << "Plane="    << plane
                       << "Chamber="  << chamber
                       << "lx0="      << localPosTracker[0]
                       << "ly0="      << localPosTracker[1]
                       << "lz0="      << localPosTracker[2]
                       << "\n";
  }

  if (plane == 5) {
     cluster->SetX(localPosTracker[0]+kX0shift5);
  }
  else {
    cluster->SetX(localPosTracker[0]+kX0shift);
  }    

  cluster->SetY(localPosTracker[1]);
  cluster->SetZ(localPosTracker[2]);

  return kTRUE;

}

//_____________________________________________________________________________
// Bool_t  AliTRDtracker::Transform(AliTRDcluster *cluster)
//{
//   //
//   // Is this still needed ????
//   //
//   const Double_t kDriftCorrection  = 1.01;                 // drift coeficient correction
//   const Double_t kTime0Cor         = 0.32;                 // time0 correction
//   //
//   const Double_t kX0shift           = 2.52; 
//   const Double_t kX0shift5          = 3.05; 

//   //
//   // apply alignment and calibration to transform cluster
//   //
//   //
//   Int_t detector = cluster->GetDetector();
//   Int_t plane   = fGeom->GetPlane(cluster->GetDetector());
//   Int_t chamber = fGeom->GetChamber(cluster->GetDetector());
//   Int_t sector  = fGeom->GetSector(cluster->GetDetector());

//   Double_t dxAmp  = (Double_t) fGeom->CamHght();          // Amplification region
//   Double_t driftX = TMath::Max(cluster->GetX()-dxAmp*0.5,0.);  // drift distance
//   //
//   // ExB correction
//   //
//   Double_t vdrift = AliTRDcalibDB::Instance()->GetVdrift(cluster->GetDetector(),0,0);
//   Double_t exB =   AliTRDcalibDB::Instance()->GetOmegaTau(vdrift,-AliTracker::GetBz()*0.1);
//   //

//   AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();  
//   AliTRDpadPlane * padPlane = commonParam->GetPadPlane(plane,chamber);
//   Double_t zshiftIdeal  = 0.5*(padPlane->GetRow0()+padPlane->GetRowEnd());
//   Double_t localPos[3], globalPos[3], localPosTracker[3], localPosTracker2[3];
//   localPos[2] = -cluster->GetX();
//   localPos[0] =  cluster->GetY() - driftX*exB;
//   localPos[1] =  cluster->GetZ() -zshiftIdeal;
//   TGeoHMatrix * matrix =  fGeom->GetGeoMatrix(cluster->GetDetector());
//   matrix->LocalToMaster(localPos, globalPos);
  
//   Double_t sectorAngle = 20.*(sector%18)+10;
//   TGeoHMatrix  rotSector;
//   rotSector.RotateZ(sectorAngle);
//   rotSector.LocalToMaster(globalPos, localPosTracker);
//   //
//   //
//   TGeoHMatrix  matrix2(*matrix);
//   matrix2.MultiplyLeft(&rotSector);
//   matrix2.LocalToMaster(localPos,localPosTracker2);
//   //
//   //
//   //
//   cluster->SetY(cluster->GetY() - driftX*exB);
//   Double_t xplane = (Double_t) AliTRDgeometry::GetTime0(plane); 
//   cluster->SetX(xplane- kDriftCorrection*(cluster->GetX()-kTime0Cor));
//   (*fDebugStreamer)<<"Transform"<<
//     "Cl.="<<cluster<<
//     "matrix.="<<matrix<<
//     "matrix2.="<<&matrix2<<
//     "Detector="<<detector<<
//     "Sector="<<sector<<
//     "Plane="<<plane<<
//     "Chamber="<<chamber<<
//     "lx0="<<localPosTracker[0]<<
//     "ly0="<<localPosTracker[1]<<
//     "lz0="<<localPosTracker[2]<<
//     "lx2="<<localPosTracker2[0]<<
//     "ly2="<<localPosTracker2[1]<<
//     "lz2="<<localPosTracker2[2]<<
//     "\n";
//   //
//   if (plane==5)
//      cluster->SetX(localPosTracker[0]+kX0shift5);
//   else
//     cluster->SetX(localPosTracker[0]+kX0shift);
    
//   cluster->SetY(localPosTracker[1]);
//   cluster->SetZ(localPosTracker[2]);
//   return kTRUE;
// }

//_____________________________________________________________________________
Bool_t AliTRDtracker::AdjustSector(AliTRDtrack *track) 
{
  //
  // Rotates the track when necessary
  //

  Double_t alpha = AliTRDgeometry::GetAlpha(); 
  Double_t y     = track->GetY();
  Double_t ymax  = track->GetX()*TMath::Tan(0.5*alpha);

  // Is this still needed ????
  //Int_t ns = AliTRDgeometry::kNsect;
  //Int_t s=Int_t(track->GetAlpha()/alpha)%ns; 

  if      (y >  ymax) {
    //s = (s+1) % ns;
    if (!track->Rotate( alpha)) {
      return kFALSE;
    }
  } 
  else if (y < -ymax) {
    //s = (s-1+ns) % ns;                           
    if (!track->Rotate(-alpha)) {
      return kFALSE;   
    }
  } 

  return kTRUE;

}

//_____________________________________________________________________________
AliTRDcluster *AliTRDtracker::GetCluster(AliTRDtrack *track, Int_t plane
                                       , Int_t timebin, UInt_t &index)
{
  //
  // Try to find cluster in the backup list
  //

  AliTRDcluster *cl =0;
  Int_t *indexes = track->GetBackupIndexes();

  for (UInt_t i = 0; i < kMaxTimeBinIndex; i++) {
    if (indexes[i] == 0) {
      break;  
    }
    AliTRDcluster *cli = (AliTRDcluster *) fClusters->UncheckedAt(indexes[i]);
    if (!cli) {
      break;
    }
    if (cli->GetLocalTimeBin() != timebin) {
      continue;
    }
    Int_t iplane = fGeom->GetPlane(cli->GetDetector());
    if (iplane == plane) {
      cl    = cli;
      index = indexes[i];
      break;
    }
  }

  return cl;

}

//_____________________________________________________________________________
Int_t  AliTRDtracker::GetLastPlane(AliTRDtrack *track)
{
  //
  // Return last updated plane
  //

  Int_t  lastplane = 0;
  Int_t *indexes   = track->GetBackupIndexes();

  for (UInt_t i = 0; i < kMaxTimeBinIndex; i++) {
    AliTRDcluster *cli = (AliTRDcluster *) fClusters->UncheckedAt(indexes[i]);
    if (!cli) {
      break;
    }
    Int_t iplane = fGeom->GetPlane(cli->GetDetector());
    if (iplane > lastplane) {
      lastplane = iplane;
    }
  }

  return lastplane;

}

//_____________________________________________________________________________
Int_t AliTRDtracker::Clusters2Tracks(AliESD *event)
{
  //
  // Finds tracks within the TRD. The ESD event is expected to contain seeds 
  // at the outer part of the TRD. The seeds
  // are found within the TRD if fAddTRDseeds is TRUE. 
  // The tracks are propagated to the innermost time bin 
  // of the TRD and the ESD event is updated
  //

  Int_t   timeBins = fTrSec[0]->GetNumberOfTimeBins();
  Float_t foundMin = fgkMinClustersInTrack * timeBins; 
  Int_t   nseed    = 0;
  Int_t   found    = 0;
  //Int_t innerTB = fTrSec[0]->GetInnerTimeBin();

  Int_t n = event->GetNumberOfTracks();
  for (Int_t i = 0; i < n; i++) {

    AliESDtrack *seed = event->GetTrack(i);
    ULong_t status = seed->GetStatus();
    if ((status & AliESDtrack::kTRDout) == 0) {
      continue;
    }
    if ((status & AliESDtrack::kTRDin)  != 0) {
      continue;
    }
    nseed++;
    
    AliTRDtrack *seed2 = new AliTRDtrack(*seed);
    //seed2->ResetCovariance(); 
    AliTRDtrack *pt    = new AliTRDtrack(*seed2,seed2->GetAlpha());
    AliTRDtrack &t     = *pt; 
    FollowProlongation(t); 
    if (t.GetNumberOfClusters() >= foundMin) {
      UseClusters(&t);
      CookLabel(pt,1 - fgkLabelFraction);
      //t.CookdEdx();
    }
    found++;

    Double_t xTPC = 250.0;
    if (PropagateToX(t,xTPC,fgkMaxStep)) {
      seed->UpdateTrackParams(pt, AliESDtrack::kTRDin);
    }  
    delete seed2;
    delete pt;

  }

  AliInfo(Form("Number of loaded seeds: %d",nseed));
  AliInfo(Form("Number of found tracks from loaded seeds: %d",found));
  AliInfo(Form("Total number of found tracks: %d",found));
    
  return 0;    

}     
     
//_____________________________________________________________________________
Int_t AliTRDtracker::PropagateBack(AliESD *event) 
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
  Int_t   n        = event->GetNumberOfTracks();

  // Sort tracks
  Float_t *quality = new Float_t[n];
  Int_t   *index   = new Int_t[n];
  for (Int_t i = 0; i < n; i++) {
    AliESDtrack *seed = event->GetTrack(i);
    Double_t covariance[15];
    seed->GetExternalCovariance(covariance);
    quality[i] = covariance[0]+covariance[2];      
    //quality[i] = covariance[0];
  }
  TMath::Sort(n,quality,index,kFALSE);

  for (Int_t i = 0; i < n; i++) {

    //AliESDtrack *seed = event->GetTrack(i);
    AliESDtrack *seed = event->GetTrack(index[i]);
    fHBackfit->Fill(0);

    ULong_t status = seed->GetStatus();
    if ((status & AliESDtrack::kTPCout) == 0) {
      fHBackfit->Fill(1);
      continue;
    }

    if ((status & AliESDtrack::kTRDout) != 0) {
      fHBackfit->Fill(2);
      continue;
    }

    Int_t   lbl = seed->GetLabel();
    AliTRDtrack *track = new AliTRDtrack(*seed);
    track->SetSeedLabel(lbl);
    seed->UpdateTrackParams(track,AliESDtrack::kTRDbackup); // Make backup
    fNseeds++;
    Float_t p4  = track->GetC();
    Int_t expectedClr = FollowBackProlongation(*track);
    
    fHBackfit->Fill(3);
    fHX->Fill(track->GetX());


    // store the last measurement
    /*
    fHNClTrack->Fill(track->GetNumberOfClusters());
    if (track->GetNumberOfClusters() >= foundMin) {

      fHBackfit->Fill(4);
      track->CookdEdx(); 
      CookdEdxTimBin(*track);
      CookLabel(track,1 - fgkLabelFraction);
      if (track->GetBackupTrack()) {
	//fHBackfit->Fill(5);
	UseClusters(track->GetBackupTrack());
	seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);
      }
    }
    */

    /**/
    // inter-tracks competition ???
    if ((TMath::Abs(track->GetC() - p4) / TMath::Abs(p4) < 0.2) || 
        (TMath::Abs(track->GetPt())                      > 0.8)) {
      
      fHBackfit->Fill(4);
      
      // 
      // Make backup for back propagation 
      //
      
      Int_t foundClr = track->GetNumberOfClusters();
      if (foundClr >= foundMin) {
	track->CookdEdx(); 
	CookdEdxTimBin(*track);
	CookLabel(track,1 - fgkLabelFraction);
	if (track->GetBackupTrack()) {
          UseClusters(track->GetBackupTrack());
	}

        // Sign only gold tracks
	if (track->GetChi2() / track->GetNumberOfClusters() < 4) {   
 	  if ((seed->GetKinkIndex(0)      ==   0) &&
              (TMath::Abs(track->GetPt()) <  1.5)) {
            UseClusters(track);
	  }
	}
	Bool_t isGold = kFALSE;
	
        // Full gold track
	if (track->GetChi2() / track->GetNumberOfClusters() < 5) {  
	  //seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	  if (track->GetBackupTrack()) {
            seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);
	    
	  }
	  isGold = kTRUE;
	  //fHBackfit->Fill()
	}

        // Almost gold track
	if ((!isGold)                                              && 
            (track->GetNCross()                              == 0) && 
            (track->GetChi2() / track->GetNumberOfClusters()  < 7)) { 
	  //seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	  if (track->GetBackupTrack()) {
            seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);
	  }
	  isGold = kTRUE;
	}

	if ((!isGold) && 
            (track->GetBackupTrack())) {
	  if ((track->GetBackupTrack()->GetNumberOfClusters() > foundMin) &&
	      ((track->GetBackupTrack()->GetChi2()/(track->GetBackupTrack()->GetNumberOfClusters()+1)) < 7)) {	  
	    seed->UpdateTrackParams(track->GetBackupTrack(),AliESDtrack::kTRDbackup);
	    isGold = kTRUE;
	  }
	}

	if ((track->StatusForTOF() > 0) &&
            (track->GetNCross() == 0) && 
            (Float_t(track->GetNumberOfClusters()) / Float_t(track->GetNExpected())  > 0.4)) {
	  //seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
	}

      }
    }
    /**/

    /**/
    // Debug part of tracking
    TTreeSRedirector &cstream = *fDebugStreamer;
    Int_t eventNr = event->GetEventNumber();
    if (AliTRDReconstructor::StreamLevel() > 0) {
      if (track->GetBackupTrack()) {
	cstream << "Tracks"
	        << "EventNr="  << eventNr
	        << "ESD.="     << seed
	        << "trd.="     << track
	        << "trdback.=" << track->GetBackupTrack()
	        << "\n";
      }
      else {
	cstream << "Tracks"
	        << "EventNr="  << eventNr
	        << "ESD.="     << seed
	        << "trd.="     << track
	        << "trdback.=" << track
	        << "\n";
      }
    }
    /**/

    // Propagation to the TOF (I.Belikov)    
    if (track->GetStop() == kFALSE) {
      fHBackfit->Fill(5);

      Double_t xtof  = 371.0;
      Double_t xTOF0 = 370.0;
    
      Double_t c2    = track->GetSnp() + track->GetC() * (xtof - track->GetX());
      if (TMath::Abs(c2) >= 0.99) {
	
	fHBackfit->Fill(6);
	delete track;
	continue;
      }
      
      PropagateToX(*track,xTOF0,fgkMaxStep);

      // Energy losses taken to the account - check one more time
      c2 = track->GetSnp() + track->GetC() * (xtof - track->GetX());
      if (TMath::Abs(c2) >= 0.99) {
	
	fHBackfit->Fill(7);
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
	  fHBackfit->Fill(8);
	  delete track;
	  continue;
	}
      } 
      else if (y < -ymax) {
	if (!track->Rotate(-AliTRDgeometry::GetAlpha())) {
	  fHBackfit->Fill(9);
	  delete track;
	  continue;
	}
      }
         
      if (track->PropagateTo(xtof)) {
	seed->UpdateTrackParams(track,AliESDtrack::kTRDout);
	fHBackfit->Fill(10);

        for (Int_t i = 0; i < AliESDtrack::kNPlane; i++) {
          for (Int_t j = 0; j < AliESDtrack::kNSlice; j++) {
            seed->SetTRDsignals(track->GetPIDsignals(i,j),i,j);
	  }
          seed->SetTRDTimBin(track->GetPIDTimBin(i),i);
        }
	//seed->SetTRDtrack(new AliTRDtrack(*track));
	if (track->GetNumberOfClusters() > foundMin) {
	  fHBackfit->Fill(11);
          found++;
	}
      }

    }
    else {
      
      fHBackfit->Fill(12);
      
      if ((track->GetNumberOfClusters() >              15) &&
          (track->GetNumberOfClusters() > 0.5*expectedClr)) {
	
	seed->UpdateTrackParams(track,AliESDtrack::kTRDout);
	fHBackfit->Fill(13);

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

    seed->SetTRDQuality(track->StatusForTOF());    
    seed->SetTRDBudget(track->GetBudget(0));    
  
    fHBackfit->Fill(14);
    delete track;
  }

  AliInfo(Form("Number of seeds: %d",fNseeds));
  AliInfo(Form("Number of back propagated TRD tracks: %d",found));
  
  // New seeding
  if (AliTRDReconstructor::SeedingOn()) {
    MakeSeedsMI(3,5,event);
  }

  fSeeds->Clear(); 
  fNseeds = 0;

  delete [] index;
  delete [] quality;
  
  SaveLogHists();
  
  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDtracker::RefitInward(AliESD *event)
{
  //
  // Refits tracks within the TRD. The ESD event is expected to contain seeds 
  // at the outer part of the TRD. 
  // The tracks are propagated to the innermost time bin 
  // of the TRD and the ESD event is updated
  // Origin: Thomas KUHR (Thomas.Kuhr@cern.ch)
  //

  Int_t   timeBins = fTrSec[0]->GetNumberOfTimeBins();
  Float_t foundMin = fgkMinClustersInTrack * timeBins; 
  Int_t   nseed    = 0;
  Int_t   found    = 0;
  //Int_t innerTB    = fTrSec[0]->GetInnerTimeBin();
  AliTRDtrack seed2;
  
  Int_t n = event->GetNumberOfTracks();
  for (Int_t i = 0; i < n; i++) {

    AliESDtrack *seed = event->GetTrack(i);
    new (&seed2) AliTRDtrack(*seed);
    fHRefit->Fill(0);

    if (seed2.GetX() < 270.0) {
      seed->UpdateTrackParams(&seed2,AliESDtrack::kTRDbackup); // Backup TPC track - only update
      fHRefit->Fill(1);
      continue;
    }

    ULong_t status = seed->GetStatus();
    if ((status & AliESDtrack::kTRDout) == 0) {
      fHRefit->Fill(2);
      continue;
    }
    if ((status & AliESDtrack::kTRDin)  != 0) {
      fHRefit->Fill(3);
      continue;
    }
    
    nseed++; 
    fHRefit->Fill(4);

    seed2.ResetCovariance(50.0); 

    AliTRDtrack *pt = new AliTRDtrack(seed2,seed2.GetAlpha());
    Int_t *indexes2 = seed2.GetIndexes();
    for (Int_t i = 0; i < AliESDtrack::kNPlane;i++) {
      for (Int_t j = 0; j < AliESDtrack::kNSlice;j++) {
        pt->SetPIDsignals(seed2.GetPIDsignals(i,j),i,j);
      }
      pt->SetPIDTimBin(seed2.GetPIDTimBin(i),i);
    }

    Int_t *indexes3 = pt->GetBackupIndexes();
    for (Int_t i = 0; i < 200;i++) {
      if (indexes2[i] == 0) {
        break;
      }
      indexes3[i] = indexes2[i];
    }  
        
    //AliTRDtrack *pt = seed2;
    AliTRDtrack &t = *pt; 
    FollowProlongation(t); 
    if (t.GetNumberOfClusters() >= foundMin) {
      //UseClusters(&t);
      //CookLabel(pt, 1-fgkLabelFraction);
      t.CookdEdx();
      CookdEdxTimBin(t);
    }
    found++;

    Double_t xTPC = 250.0;
    if (PropagateToX(t,xTPC,fgkMaxStep)) {

      seed->UpdateTrackParams(pt,AliESDtrack::kTRDrefit);
      fHRefit->Fill(5);

      for (Int_t i = 0; i < AliESDtrack::kNPlane; i++) {
        for (Int_t j = 0; j < AliESDtrack::kNSlice; j++) {
          seed->SetTRDsignals(pt->GetPIDsignals(i,j),i,j);
	}
        seed->SetTRDTimBin(pt->GetPIDTimBin(i),i);
      }
    }
    else {
      // If not prolongation to TPC - propagate without update
      fHRefit->Fill(5);
      AliTRDtrack *seed2 = new AliTRDtrack(*seed);
      seed2->ResetCovariance(5.0); 
      AliTRDtrack *pt2   = new AliTRDtrack(*seed2,seed2->GetAlpha());
      delete seed2;
      if (PropagateToX(*pt2,xTPC,fgkMaxStep)) { 
        pt2->CookdEdx( ); 
        CookdEdxTimBin(*pt2);
	seed->UpdateTrackParams(pt2,AliESDtrack::kTRDrefit);
	fHRefit->Fill(6);

        for (Int_t i = 0; i < AliESDtrack::kNPlane; i++) {
          for (Int_t j = 0; j < AliESDtrack::kNSlice; j++) {
            seed->SetTRDsignals(pt2->GetPIDsignals(i,j),i,j);
	  }
          seed->SetTRDTimBin(pt2->GetPIDTimBin(i),i);
        }
      }
      delete pt2;
    }

    delete pt;

  }   

  AliInfo(Form("Number of loaded seeds: %d",nseed));
  AliInfo(Form("Number of found tracks from loaded seeds: %d",found));

  SaveLogHists();
  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDtracker::FollowProlongation(AliTRDtrack &t)
{
  //
  // Starting from current position on track=t this function tries
  // to extrapolate the track up to timeBin=0 and to confirm prolongation
  // if a close cluster is found. Returns the number of clusters
  // expected to be found in sensitive layers
  // GeoManager used to estimate mean density
  //

  Int_t    sector;
  Int_t    lastplane = GetLastPlane(&t);
  Double_t radLength = 0.0;
  Double_t rho       = 0.0;
  Int_t    expectedNumberOfClusters = 0;

  for (Int_t iplane = lastplane; iplane >= 0; iplane--) {

    Int_t row0    = GetGlobalTimeBin(0,iplane,GetTimeBinsPerPlane()-1);
    Int_t rowlast = GetGlobalTimeBin(0,iplane,0);

    //
    // Propagate track close to the plane if neccessary
    //
    Double_t currentx  = fTrSec[0]->GetLayer(rowlast)->GetX();
    if (currentx < (-fgkMaxStep + t.GetX())) {
      // Propagate closer to chamber - safety space fgkMaxStep      
      if (!PropagateToX(t,currentx+fgkMaxStep,fgkMaxStep)) {
        break;
      }
    }

    if (!AdjustSector(&t)) {
      break;
    }

    //
    // Get material budget
    //
    Double_t xyz0[3];
    Double_t xyz1[3];
    Double_t param[7];
    Double_t x;
    Double_t y;
    Double_t z;

    // Starting global position
    t.GetXYZ(xyz0);   
    // End global position
    x = fTrSec[0]->GetLayer(row0)->GetX();
    if (!t.GetProlongation(x,y,z)) {
      break;
    }
    xyz1[0] =  x * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha()); 
    xyz1[1] = +x * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
    rho       = param[0];
    radLength = param[1]; // Get mean propagation parameters

    //
    // Propagate and update
    //
    sector = t.GetSector();
    //for (Int_t itime=GetTimeBinsPerPlane()-1;itime>=0;itime--) {
    for (Int_t itime = 0 ; itime < GetTimeBinsPerPlane(); itime++) {

      Int_t ilayer = GetGlobalTimeBin(0,iplane,itime);
      expectedNumberOfClusters++;       
      t.SetNExpected(t.GetNExpected() + 1);
      if (t.GetX() > 345.0) {
        t.SetNExpectedLast(t.GetNExpectedLast() + 1);
      }
      AliTRDpropagationLayer &timeBin = *(fTrSec[sector]->GetLayer(ilayer));
      AliTRDcluster *cl = 0;
      UInt_t   index   = 0;
      Double_t maxChi2 = fgkMaxChi2;
      x = timeBin.GetX();

      if (timeBin) {

	AliTRDcluster *cl0 = timeBin[0];
	if (!cl0) {
          // No clusters in given time bin
          continue;         
	}

	Int_t plane   = fGeom->GetPlane(cl0->GetDetector());
	if (plane > lastplane) {
          continue;
	}

	Int_t timebin = cl0->GetLocalTimeBin();
	AliTRDcluster *cl2 = GetCluster(&t,plane,timebin,index);

	if (cl2) {

	  cl = cl2;	
	  //Double_t h01 = GetTiltFactor(cl);    //I.B's fix
	  //maxChi2=t.GetPredictedChi2(cl,h01);

	}	
        if (cl) {

	  //if (cl->GetNPads()<5) 
	  Double_t dxsample = timeBin.GetdX();
	  t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 
          Double_t h01      = GetTiltFactor(cl);
	  Int_t    det      = cl->GetDetector();    
	  Int_t    plane    = fGeom->GetPlane(det);
	  if (t.GetX() > 345.0) {
	    t.SetNLast(t.GetNLast() + 1);
	    t.SetChi2Last(t.GetChi2Last() + maxChi2);
	  }

	  Double_t xcluster = cl->GetX();
	  t.PropagateTo(xcluster,radLength,rho);

          if (!AdjustSector(&t)) {
            break;           //I.B's fix
	  }
	  maxChi2 = t.GetPredictedChi2(cl,h01);
          
	  if (maxChi2<1e+10)
	    if (!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
	      // ????
	    }

	}			

      }

    } 

  }

  return expectedNumberOfClusters;  

}                

//_____________________________________________________________________________
Int_t AliTRDtracker::FollowBackProlongation(AliTRDtrack &t)
{
  //  
  // Starting from current radial position of track <t> this function
  // extrapolates the track up to outer timebin and in the sensitive
  // layers confirms prolongation if a close cluster is found. 
  // Returns the number of clusters expected to be found in sensitive layers
  // Use GEO manager for material Description
  //
  // return number of assigned clusters ? 
  //


  Int_t    sector;
  Int_t    clusters[1000];
  Double_t radLength = 0.0;
  Double_t rho       = 0.0;
  Int_t    expectedNumberOfClusters = 0;
  Float_t  ratio0    = 0.0;
  AliTRDtracklet tracklet;

  // Calibration fill 2D
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();
  if (!calibra) {
    AliInfo("Could not get Calibra instance\n");
  }
  if (calibra->GetMITracking()) {
    calibra->ResetTrack();
  }

  for (Int_t i = 0; i < 1000; i++) {
    clusters[i] = -1;
  }

  for (Int_t iplane = 0; iplane < AliESDtrack::kNPlane; iplane++) {

    int hb = iplane * 10;
    fHClSearch->Fill(hb);

    Int_t    row0     = GetGlobalTimeBin(0,iplane,GetTimeBinsPerPlane()-1);
    Int_t    rowlast  = GetGlobalTimeBin(0,iplane,0);
    Double_t currentx = fTrSec[0]->GetLayer(row0)->GetX();
    if (currentx < t.GetX()) {
      fHClSearch->Fill(hb+1);
      continue;
    }

    //
    // Propagate closer to chamber if neccessary 
    //
    if (currentx > (fgkMaxStep + t.GetX())) {
      if (!PropagateToX(t,currentx-fgkMaxStep,fgkMaxStep)) {
	fHClSearch->Fill(hb+2);
        break;
      }
    }
    if (!AdjustSector(&t)) {
      fHClSearch->Fill(hb+3);
      break;
    }
    if (TMath::Abs(t.GetSnp()) > fgkMaxSnp) {
      fHClSearch->Fill(hb+4);
      break;
    }

    //
    // Get material budget inside of chamber
    //
    Double_t xyz0[3];
    Double_t xyz1[3];
    Double_t param[7];
    Double_t x;
    Double_t y;
    Double_t z;
    // Starting global position
    t.GetXYZ(xyz0);   
    // End global position
    x = fTrSec[0]->GetLayer(rowlast)->GetX();
    if (!t.GetProlongation(x,y,z)) {
      fHClSearch->Fill(hb+5);
      break;
    }
    xyz1[0] =  x * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha()); 
    xyz1[1] = +x * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
    rho       = param[0];
    radLength = param[1]; // Get mean propagation parameters

    //
    // Find clusters
    //
    sector = t.GetSector();
    Float_t ncl = FindClusters(sector,row0,rowlast,&t,clusters,tracklet);
    fHNCl->Fill(tracklet.GetN());

    if (tracklet.GetN() < GetTimeBinsPerPlane()/3) {
      fHClSearch->Fill(hb+6);
      continue;
    }

    //
    // Propagate and update track
    //
    for (Int_t itime = GetTimeBinsPerPlane()-1; itime >= 0; itime--) {

      Int_t ilayer = GetGlobalTimeBin(0, iplane,itime);
      expectedNumberOfClusters++;       
      t.SetNExpected(t.GetNExpected() + 1);
      if (t.GetX() > 345.0) {
        t.SetNExpectedLast(t.GetNExpectedLast() + 1);
      }
      AliTRDpropagationLayer &timeBin = *(fTrSec[sector]->GetLayer(ilayer));
      AliTRDcluster *cl = 0;
      UInt_t   index   = 0;
      Double_t maxChi2 = fgkMaxChi2;
      x = timeBin.GetX();

      if (timeBin) {	

	if (clusters[ilayer] > 0) {
	  index = clusters[ilayer];
	  cl    = (AliTRDcluster *)GetCluster(index);
	  //Double_t h01 = GetTiltFactor(cl);   // I.B's fix
          //maxChi2=t.GetPredictedChi2(cl,h01); //
	}
	
        if (cl) {

	  //if (cl->GetNPads() < 5) 
	  Double_t dxsample = timeBin.GetdX();
	  t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 
          Double_t h01      = GetTiltFactor(cl);
	  Int_t    det      = cl->GetDetector();    
	  Int_t    plane    = fGeom->GetPlane(det);
	  if (t.GetX() > 345.0) {
	    t.SetNLast(t.GetNLast() + 1);
	    t.SetChi2Last(t.GetChi2Last() + maxChi2);
	  }
	  Double_t xcluster = cl->GetX();
	  t.PropagateTo(xcluster,radLength,rho);
	  maxChi2 = t.GetPredictedChi2(cl,h01);

	  if (maxChi2<1e+10)
	    if (!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
	      if (!t.Update(cl,maxChi2,index,h01)) {
		// ????
	      }
	    }  

          if (calibra->GetMITracking()) {
            calibra->UpdateHistograms(cl,&t);
          }

	  // Reset material budget if 2 consecutive gold
	  if (plane > 0) { 
	    if ((t.GetTracklets(plane).GetN() + t.GetTracklets(plane-1).GetN()) > 20) {
	      t.SetBudget(2,0.0);
	    }
	  }

	}

      }

    }

    ratio0 = ncl / Float_t(fTimeBinsPerPlane);
    Float_t ratio1 = Float_t(t.GetNumberOfClusters()+1) / Float_t(t.GetNExpected()+1);	
    if ((tracklet.GetChi2()      <  18.0) && 
        (ratio0                  >   0.8) && 
        (ratio1                  >   0.6) && 
        (ratio0+ratio1           >   1.5) && 
        (t.GetNCross()           ==    0) && 
        (TMath::Abs(t.GetSnp())  <  0.85) &&
        (t.GetNumberOfClusters() >    20)){
      //if (ratio0 > 0.8) {
      t.MakeBackupTrack(); // Make backup of the track until is gold
    }
    
  }

  return expectedNumberOfClusters;  
}         

//_____________________________________________________________________________
Int_t AliTRDtracker::PropagateToX(AliTRDtrack &t, Double_t xToGo, Double_t maxStep)
{
  //
  // Starting from current radial position of track <t> this function
  // extrapolates the track up to radial position <xToGo>. 
  // Returns 1 if track reaches the plane, and 0 otherwise 
  //

  const Double_t kEpsilon = 0.00001;
  //Double_t tanmax = TMath::Tan(0.5*AliTRDgeometry::GetAlpha()); 
  Double_t xpos = t.GetX();
  Double_t dir  = (xpos<xToGo) ? 1.0 : -1.0;

  while (((xToGo-xpos)*dir) > kEpsilon) {

    Double_t step = dir * TMath::Min(TMath::Abs(xToGo-xpos),maxStep);

    Double_t xyz0[3];
    Double_t xyz1[3];
    Double_t param[7];
    Double_t x;
    Double_t y;
    Double_t z;
    // Starting global position
    t.GetXYZ(xyz0);   
    x = xpos + step;

    if (!t.GetProlongation(x,y,z)) {
      return 0; // No prolongation
    }

    xyz1[0] =  x * TMath::Cos(t.GetAlpha()) - y * TMath::Sin(t.GetAlpha()); 
    xyz1[1] = +x * TMath::Sin(t.GetAlpha()) + y * TMath::Cos(t.GetAlpha());
    xyz1[2] =  z;

    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
    if (!t.PropagateTo(x,param[1],param[0])) {
      return 0;
    }
    AdjustSector(&t);
    xpos = t.GetX();

  }

  return 1;

}

//_____________________________________________________________________________
Int_t AliTRDtracker::LoadClusters(TTree *cTree)
{
  //
  // Fills clusters into TRD tracking_sectors 
  // Note that the numbering scheme for the TRD tracking_sectors 
  // differs from that of TRD sectors
  //

  if (ReadClusters(fClusters,cTree)) {
     AliError("Problem with reading the clusters !");
     return 1;
  }
  Int_t ncl = fClusters->GetEntriesFast();
  fNclusters = ncl;
  AliInfo(Form("Sorting %d clusters",ncl));
              
  UInt_t index;
  for (Int_t ichamber = 0; ichamber <  5; ichamber++) {
    for (Int_t isector  = 0; isector  < 18;  isector++) {
      fHoles[ichamber][isector] = kTRUE;
    }
  }

  while (ncl--) {

    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(ncl);
    Int_t detector       = c->GetDetector();
    Int_t localTimeBin   = c->GetLocalTimeBin();
    Int_t sector         = fGeom->GetSector(detector);
    Int_t plane          = fGeom->GetPlane(detector);
    Int_t trackingSector = CookSectorIndex(sector);

    //if (c->GetLabel(0) > 0) {
    if (c->GetQ() > 10) {
      Int_t chamber = fGeom->GetChamber(detector);
      fHoles[chamber][trackingSector] = kFALSE;
    }

    Int_t gtb = fTrSec[trackingSector]->CookTimeBinIndex(plane,localTimeBin);
    if (gtb < 0) {
      continue; 
    }
    Int_t layer = fTrSec[trackingSector]->GetLayerNumber(gtb);

    index = ncl;

    // Apply pos correction
    Transform(c);  
    fHXCl->Fill(c->GetX());

    fTrSec[trackingSector]->GetLayer(layer)->SetX(c->GetX());
    fTrSec[trackingSector]->GetLayer(layer)->InsertCluster(c,index);
  }

  return 0;

}

//_____________________________________________________________________________
void AliTRDtracker::UnloadClusters() 
{ 
  //
  // Clears the arrays of clusters and tracks. Resets sectors and timebins 
  //

  Int_t i;
  Int_t nentr;

  nentr = fClusters->GetEntriesFast();
  for (i = 0; i < nentr; i++) {
    delete fClusters->RemoveAt(i);
  }
  fNclusters = 0;

  nentr = fSeeds->GetEntriesFast();
  for (i = 0; i < nentr; i++) {
    delete fSeeds->RemoveAt(i);
  }

  nentr = fTracks->GetEntriesFast();
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

//_____________________________________________________________________________
void AliTRDtracker::MakeSeedsMI(Int_t /*inner*/, Int_t /*outer*/, AliESD *esd)
{
  //
  // Creates  seeds using clusters between  position inner plane  and outer plane 
  //

  const Double_t kMaxTheta =  1.0;
  const Double_t kMaxPhi   =  2.0;

  const Double_t kRoad0y   =  6.0; // Road for middle cluster 
  const Double_t kRoad0z   =  8.5; // Road for middle cluster 

  const Double_t kRoad1y   =  2.0; // Road in y for seeded cluster
  const Double_t kRoad1z   = 20.0; // Road in z for seeded cluster

  const Double_t kRoad2y   =  3.0; // Road in y for extrapolated cluster
  const Double_t kRoad2z   = 20.0; // Road in z for extrapolated cluster
  const Int_t    kMaxSeed  = 3000;

  Int_t maxSec = AliTRDgeometry::kNsect;  

  // Linear fitters in planes
  TLinearFitter fitterTC(2,"hyp2");  // Fitting with tilting pads - kz fixed - kz= Z/x, + vertex const
  TLinearFitter fitterT2(4,"hyp4");  // Fitting with tilting pads - kz not fixed
  fitterTC.StoreData(kTRUE);
  fitterT2.StoreData(kTRUE);
  AliRieman rieman(1000);            // Rieman fitter
  AliRieman rieman2(1000);           // Rieman fitter

  // Find the maximal and minimal layer for the planes
  Int_t layers[6][2];
  AliTRDpropagationLayer *reflayers[6];
  for (Int_t i = 0; i < 6; i++) {
    layers[i][0] = 10000; 
    layers[i][1] = 0;
  }
  for (Int_t ns = 0; ns < maxSec; ns++) {
    for (Int_t ilayer = 0; ilayer < fTrSec[ns]->GetNumberOfLayers(); ilayer++) {
      AliTRDpropagationLayer &layer = *(fTrSec[ns]->GetLayer(ilayer));
      if (layer == 0) {
        continue;
      }
      Int_t det   = layer[0]->GetDetector();    
      Int_t plane = fGeom->GetPlane(det);
      if (ilayer < layers[plane][0]) {
        layers[plane][0] = ilayer;
      }
      if (ilayer > layers[plane][1]) {
        layers[plane][1] = ilayer;
      }
    }
  }

  AliTRDpadPlane *padPlane = AliTRDCommonParam::Instance()->GetPadPlane(0,0);
  Double_t h01 = TMath::Tan(-TMath::Pi() / 180.0 * padPlane->GetTiltingAngle());
  Double_t hL[6];                                                // Tilting angle
  Double_t xcl[6];                                               // X - position of reference cluster
  Double_t ycl[6];                                               // Y - position of reference cluster
  Double_t zcl[6];                                               // Z - position of reference cluster

  AliTRDcluster *cl[6] = {    0,    0,    0,    0,    0,    0 }; // Seeding clusters
  Float_t padlength[6] = { 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 }; // Current pad-length 

  Double_t chi2R  = 0.0;
  Double_t chi2Z  = 0.0;
  Double_t chi2RF = 0.0;
  Double_t chi2ZF = 0.0;

  Int_t nclusters;                                               // Total number of clusters
  for (Int_t i = 0; i < 6; i++) {
    hL[i] = h01; 
    if (i%2==1) {
      hL[i]*=-1.0;
    }
  }

  // Registered seed
  AliTRDseed *pseed = new AliTRDseed[kMaxSeed*6];
  AliTRDseed *seed[kMaxSeed];
  for (Int_t iseed = 0; iseed < kMaxSeed; iseed++) {
    seed[iseed]= &pseed[iseed*6];
  }
  AliTRDseed *cseed = seed[0];

  Double_t seedquality[kMaxSeed];  
  Double_t seedquality2[kMaxSeed];  
  Double_t seedparams[kMaxSeed][7];
  Int_t    seedlayer[kMaxSeed];
  Int_t    registered = 0;
  Int_t    sort[kMaxSeed];

  //
  // Seeding part
  //
  for (Int_t ns = 0; ns < maxSec; ns++) { // Loop over sectors
  //for (Int_t ns = 0; ns < 5; ns++) { // Loop over sectors

    registered   = 0; // Reset registerd seed counter
    cseed        = seed[registered];
    Float_t iter = 0.0;

    for (Int_t sLayer = 2; sLayer >= 0; sLayer--) {
      //for (Int_t dseed = 5; dseed < 15; dseed += 3) {

      iter += 1.0;
      Int_t dseed = 5 + Int_t(iter) * 3;

      // Initialize seeding layers
      for (Int_t ilayer = 0; ilayer < 6; ilayer++) {
	reflayers[ilayer] = fTrSec[ns]->GetLayer(layers[ilayer][1]-dseed);
	xcl[ilayer]       = reflayers[ilayer]->GetX();
      }

      Double_t xref = (xcl[sLayer+1] + xcl[sLayer+2]) * 0.5;      
      AliTRDpropagationLayer &layer0 = *reflayers[sLayer+0];
      AliTRDpropagationLayer &layer1 = *reflayers[sLayer+1];
      AliTRDpropagationLayer &layer2 = *reflayers[sLayer+2];
      AliTRDpropagationLayer &layer3 = *reflayers[sLayer+3];

      Int_t maxn3 = layer3;
      for (Int_t icl3 = 0; icl3 < maxn3; icl3++) {

	AliTRDcluster *cl3 = layer3[icl3];
	if (!cl3) {
          continue;
	}
	padlength[sLayer+3] = TMath::Sqrt(cl3->GetSigmaZ2() * 12.0);
	ycl[sLayer+3] = cl3->GetY();
	zcl[sLayer+3] = cl3->GetZ();
	Float_t yymin0 = ycl[sLayer+3] - 1.0 - kMaxPhi * (xcl[sLayer+3]-xcl[sLayer+0]);
	Float_t yymax0 = ycl[sLayer+3] + 1.0 + kMaxPhi * (xcl[sLayer+3]-xcl[sLayer+0]);
	Int_t   maxn0  = layer0;

	for (Int_t icl0 = layer0.Find(yymin0); icl0 < maxn0; icl0++) {

	  AliTRDcluster *cl0 = layer0[icl0];
	  if (!cl0) {
            continue;
	  }
	  if (cl3->IsUsed() && cl0->IsUsed()) {
            continue;
	  }
	  ycl[sLayer+0] = cl0->GetY();
	  zcl[sLayer+0] = cl0->GetZ();
	  if (ycl[sLayer+0] > yymax0) {
            break;
	  }
	  Double_t tanphi     = (ycl[sLayer+3]-ycl[sLayer+0]) / (xcl[sLayer+3]-xcl[sLayer+0]); 
	  if (TMath::Abs(tanphi)   > kMaxPhi) {
            continue;
	  }
	  Double_t tantheta   = (zcl[sLayer+3]-zcl[sLayer+0]) / (xcl[sLayer+3]-xcl[sLayer+0]); 
	  if (TMath::Abs(tantheta) > kMaxTheta) {
            continue; 
	  }
	  padlength[sLayer+0] = TMath::Sqrt(cl0->GetSigmaZ2() * 12.0);

	  // Expected position in 1 layer
	  Double_t y1exp = ycl[sLayer+0] + (tanphi)   * (xcl[sLayer+1]-xcl[sLayer+0]);	  
	  Double_t z1exp = zcl[sLayer+0] + (tantheta) * (xcl[sLayer+1]-xcl[sLayer+0]);	  
	  Float_t yymin1 = y1exp - kRoad0y - tanphi;
	  Float_t yymax1 = y1exp + kRoad0y + tanphi;
	  Int_t   maxn1  = layer1;

	  for (Int_t icl1 = layer1.Find(yymin1); icl1 < maxn1; icl1++) {

	    AliTRDcluster *cl1 = layer1[icl1];
	    if (!cl1) {
              continue;
	    }
	    Int_t nusedCl = 0;
	    if (cl3->IsUsed()) nusedCl++;
	    if (cl0->IsUsed()) nusedCl++;
	    if (cl1->IsUsed()) nusedCl++;
	    if (nusedCl > 1) {
              continue;
	    }
	    ycl[sLayer+1] = cl1->GetY();
	    zcl[sLayer+1] = cl1->GetZ();
	    if (ycl[sLayer+1] > yymax1) {
              break;
	    }
	    if (TMath::Abs(ycl[sLayer+1]-y1exp) > kRoad0y+tanphi) {
              continue;
	    }
	    if (TMath::Abs(zcl[sLayer+1]-z1exp) > kRoad0z) {
              continue;
	    }
	    padlength[sLayer+1] = TMath::Sqrt(cl1->GetSigmaZ2() * 12.0);

	    Double_t y2exp  = ycl[sLayer+0]+(tanphi)   * (xcl[sLayer+2]-xcl[sLayer+0]) + (ycl[sLayer+1]-y1exp);	  
	    Double_t z2exp  = zcl[sLayer+0]+(tantheta) * (xcl[sLayer+2]-xcl[sLayer+0]);
	    Int_t    index2 = layer2.FindNearestCluster(y2exp,z2exp,kRoad1y,kRoad1z);
	    if (index2 <= 0) {
              continue; 
	    }
	    AliTRDcluster *cl2 = (AliTRDcluster *) GetCluster(index2);
	    padlength[sLayer+2] = TMath::Sqrt(cl2->GetSigmaZ2() * 12.0);
	    ycl[sLayer+2]       = cl2->GetY();
	    zcl[sLayer+2]       = cl2->GetZ();
	    if (TMath::Abs(cl2->GetZ()-z2exp) > kRoad0z) {
              continue;
	    }

 	    rieman.Reset();
 	    rieman.AddPoint(xcl[sLayer+0],ycl[sLayer+0],zcl[sLayer+0],1,10);
 	    rieman.AddPoint(xcl[sLayer+1],ycl[sLayer+1],zcl[sLayer+1],1,10);
 	    rieman.AddPoint(xcl[sLayer+3],ycl[sLayer+3],zcl[sLayer+3],1,10);	    
	    rieman.AddPoint(xcl[sLayer+2],ycl[sLayer+2],zcl[sLayer+2],1,10);
	    rieman.Update();

	    // Reset fitter
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      cseed[iLayer].Reset();
	    }	  
	    chi2Z = 0.0; 
            chi2R = 0.0;

	    for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
	      cseed[sLayer+iLayer].SetZref(0,rieman.GetZat(xcl[sLayer+iLayer]));
	      chi2Z += (cseed[sLayer+iLayer].GetZref(0)- zcl[sLayer+iLayer])
		     * (cseed[sLayer+iLayer].GetZref(0)- zcl[sLayer+iLayer]);
	      cseed[sLayer+iLayer].SetZref(1,rieman.GetDZat(xcl[sLayer+iLayer]));	      
	      cseed[sLayer+iLayer].SetYref(0,rieman.GetYat(xcl[sLayer+iLayer]));
	      chi2R += (cseed[sLayer+iLayer].GetYref(0)- ycl[sLayer+iLayer])
		     * (cseed[sLayer+iLayer].GetYref(0)- ycl[sLayer+iLayer]);
	      cseed[sLayer+iLayer].SetYref(1,rieman.GetDYat(xcl[sLayer+iLayer]));
	    }
	    if (TMath::Sqrt(chi2R) > 1.0/iter) {
              continue;
	    }
	    if (TMath::Sqrt(chi2Z) > 7.0/iter) {
              continue;
	    }

	    Float_t minmax[2] = { -100.0,  100.0 };
	    for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
	      Float_t max = zcl[sLayer+iLayer]+padlength[sLayer+iLayer] * 0.5
                          + 1.0 - cseed[sLayer+iLayer].GetZref(0);
	      if (max < minmax[1]) {
                minmax[1] = max; 
	      }
	      Float_t min = zcl[sLayer+iLayer]-padlength[sLayer+iLayer] * 0.5
                          - 1.0 - cseed[sLayer+iLayer].GetZref(0);
	      if (min > minmax[0]) {
                minmax[0] = min; 
	      }
	    }

	    Bool_t isFake = kFALSE; 
	    if (cl0->GetLabel(0) != cl3->GetLabel(0)) {
              isFake = kTRUE;
	    }
	    if (cl1->GetLabel(0) != cl3->GetLabel(0)) {
              isFake = kTRUE;
	    }
	    if (cl2->GetLabel(0) != cl3->GetLabel(0)) {
              isFake = kTRUE;
	    }

	    if (AliTRDReconstructor::StreamLevel() > 0) {
	      if ((!isFake) || ((icl3%10) == 0)) {  // Debugging print
		TTreeSRedirector &cstream = *fDebugStreamer;
		cstream << "Seeds0"
		        << "isFake=" << isFake
			<< "Cl0.="   << cl0
			<< "Cl1.="   << cl1
			<< "Cl2.="   << cl2
			<< "Cl3.="   << cl3
			<< "Xref="   << xref
			<< "X0="     << xcl[sLayer+0]
			<< "X1="     << xcl[sLayer+1]
			<< "X2="     << xcl[sLayer+2]
			<< "X3="     << xcl[sLayer+3]
			<< "Y2exp="  << y2exp
			<< "Z2exp="  << z2exp
			<< "Chi2R="  << chi2R
			<< "Chi2Z="  << chi2Z
			<< "Seed0.=" << &cseed[sLayer+0]
			<< "Seed1.=" << &cseed[sLayer+1]
			<< "Seed2.=" << &cseed[sLayer+2]
			<< "Seed3.=" << &cseed[sLayer+3]
			<< "Zmin="   << minmax[0]
			<< "Zmax="   << minmax[1]
			<< "\n";
	      }
	    }

	    ////////////////////////////////////////////////////////////////////////////////////
	    //
	    // Fit seeding part
	    //
	    ////////////////////////////////////////////////////////////////////////////////////

	    cl[sLayer+0] = cl0;
	    cl[sLayer+1] = cl1;
	    cl[sLayer+2] = cl2;
	    cl[sLayer+3] = cl3;
	    Bool_t isOK = kTRUE;

	    for (Int_t jLayer = 0; jLayer < 4; jLayer++) {

	      cseed[sLayer+jLayer].SetTilt(hL[sLayer+jLayer]);
	      cseed[sLayer+jLayer].SetPadLength(padlength[sLayer+jLayer]);
	      cseed[sLayer+jLayer].SetX0(xcl[sLayer+jLayer]);

	      for (Int_t iter = 0; iter < 2; iter++) {

		//
		// In iteration 0 we try only one pad-row
		// If quality not sufficient we try 2 pad-rows - about 5% of tracks cross 2 pad-rows
		//
		AliTRDseed tseed = cseed[sLayer+jLayer];
		Float_t    roadz = padlength[sLayer+jLayer] * 0.5;
		if (iter > 0) {
                  roadz = padlength[sLayer+jLayer];
		}

		Float_t quality = 10000.0;

		for (Int_t iTime = 2; iTime < 20; iTime++) { 

		  AliTRDpropagationLayer &layer = *(fTrSec[ns]->GetLayer(layers[sLayer+jLayer][1]-iTime));
		  Double_t dxlayer = layer.GetX() - xcl[sLayer+jLayer];		 
		  Double_t zexp    = cl[sLayer+jLayer]->GetZ();

		  if (iter > 0) {
		    // Try 2 pad-rows in second iteration
		    zexp = tseed.GetZref(0) + tseed.GetZref(1) * dxlayer;
		    if (zexp > cl[sLayer+jLayer]->GetZ()) {
                      zexp = cl[sLayer+jLayer]->GetZ() + padlength[sLayer+jLayer]*0.5;
		    }
		    if (zexp < cl[sLayer+jLayer]->GetZ()) {
                      zexp = cl[sLayer+jLayer]->GetZ() - padlength[sLayer+jLayer]*0.5;
		    }
		  }

		  Double_t yexp  = tseed.GetYref(0) + tseed.GetYref(1) * dxlayer;
		  Int_t    index = layer.FindNearestCluster(yexp,zexp,kRoad1y,roadz);
		  if (index <= 0) {
                    continue; 
		  }
		  AliTRDcluster *cl = (AliTRDcluster *) GetCluster(index);	      

		  tseed.SetIndexes(iTime,index);
		  tseed.SetClusters(iTime,cl);  // Register cluster
		  tseed.SetX(iTime,dxlayer);    // Register cluster
		  tseed.SetY(iTime,cl->GetY()); // Register cluster
		  tseed.SetZ(iTime,cl->GetZ()); // Register cluster

		}

		tseed.Update();

		// Count the number of clusters and distortions into quality
		Float_t dangle   = tseed.GetYfit(1) - tseed.GetYref(1);
		Float_t tquality = (18.0 - tseed.GetN2()) / 2.0 + TMath::Abs(dangle) / 0.1
                                 + TMath::Abs(tseed.GetYfit(0) - tseed.GetYref(0))   / 0.2
                                 + 2.0 * TMath::Abs(tseed.GetMeanz() - tseed.GetZref(0)) / padlength[jLayer];
		if ((iter == 0) && tseed.IsOK()) {
		  cseed[sLayer+jLayer] = tseed;
		  quality              = tquality;
		  if (tquality < 5) {
                    break;  
		  }
		}
		if (tseed.IsOK() && (tquality < quality)) {
		  cseed[sLayer+jLayer] = tseed;
		}
	
	      } // Loop: iter

	      if (!cseed[sLayer+jLayer].IsOK()) {
		isOK = kFALSE;
		break;
	      }

	      cseed[sLayer+jLayer].CookLabels();
	      cseed[sLayer+jLayer].UpdateUsed();
	      nusedCl += cseed[sLayer+jLayer].GetNUsed();
	      if (nusedCl > 25) {
		isOK = kFALSE;
		break;
	      }

	    } // Loop: jLayer

	    if (!isOK) {
              continue;
	    }
	    nclusters = 0;
	    for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
	      if (cseed[sLayer+iLayer].IsOK()) {
		nclusters += cseed[sLayer+iLayer].GetN2();	    
	      }
	    }

	    // Iteration 0
	    rieman.Reset();
	    for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
	      rieman.AddPoint(xcl[sLayer+iLayer]
			     ,cseed[sLayer+iLayer].GetYfitR(0)
			     ,cseed[sLayer+iLayer].GetZProb()
                             ,1
                             ,10);
	    }
	    rieman.Update();

	    chi2R = 0.0; 
            chi2Z = 0.0;

	    for (Int_t iLayer = 0; iLayer < 4; iLayer++) {
	      cseed[sLayer+iLayer].SetYref(0,rieman.GetYat(xcl[sLayer+iLayer]));
	      chi2R += (cseed[sLayer+iLayer].GetYref(0) - cseed[sLayer+iLayer].GetYfitR(0))
                     * (cseed[sLayer+iLayer].GetYref(0) - cseed[sLayer+iLayer].GetYfitR(0));
	      cseed[sLayer+iLayer].SetYref(1,rieman.GetDYat(xcl[sLayer+iLayer]));
	      cseed[sLayer+iLayer].SetZref(0,rieman.GetZat(xcl[sLayer+iLayer]));
	      chi2Z += (cseed[sLayer+iLayer].GetZref(0) - cseed[sLayer+iLayer].GetMeanz())
                     * (cseed[sLayer+iLayer].GetZref(0) - cseed[sLayer+iLayer].GetMeanz());
	      cseed[sLayer+iLayer].SetZref(1,rieman.GetDZat(xcl[sLayer+iLayer]));
	    }
	    Double_t curv = rieman.GetC();

	    //
	    // Likelihoods
	    //
	    Double_t sumda     = TMath::Abs(cseed[sLayer+0].GetYfitR(1) - cseed[sLayer+0].GetYref(1))
                               + TMath::Abs(cseed[sLayer+1].GetYfitR(1) - cseed[sLayer+1].GetYref(1))
                               + TMath::Abs(cseed[sLayer+2].GetYfitR(1) - cseed[sLayer+2].GetYref(1))
                               + TMath::Abs(cseed[sLayer+3].GetYfitR(1) - cseed[sLayer+3].GetYref(1));
	    Double_t likea     = TMath::Exp(-sumda*10.6);
	    Double_t likechi2  = 0.0000000001;
	    if (chi2R < 0.5) {
              likechi2 += TMath::Exp(-TMath::Sqrt(chi2R) * 7.73);
	    }
	    Double_t likechi2z = TMath::Exp(-chi2Z * 0.088) / TMath::Exp(-chi2Z * 0.019);
	    Double_t likeN     = TMath::Exp(-(72 - nclusters) * 0.19);
	    Double_t like      = likea * likechi2 * likechi2z * likeN;
	    Double_t likePrimY = TMath::Exp(-TMath::Abs(cseed[sLayer+0].GetYref(1) - 130.0*curv) * 1.9);
	    Double_t likePrimZ = TMath::Exp(-TMath::Abs(cseed[sLayer+0].GetZref(1)
                                                      - cseed[sLayer+0].GetZref(0) / xcl[sLayer+0]) * 5.9);
	    Double_t likePrim  = TMath::Max(likePrimY*likePrimZ,0.0005);
					    
	    seedquality[registered] = like; 
	    seedlayer[registered]   = sLayer;
	    if (TMath::Log(0.000000000000001 + like) < -15) {
              continue;
	    }
	    AliTRDseed seedb[6];
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      seedb[iLayer] = cseed[iLayer]; 
	    }

	    ////////////////////////////////////////////////////////////////////////////////////
	    //
	    // Full track fit part
	    //
	    ////////////////////////////////////////////////////////////////////////////////////

	    Int_t nlayers  = 0;
	    Int_t nusedf   = 0;
	    Int_t findable = 0;

	    //
	    // Add new layers  - avoid long extrapolation
	    //
	    Int_t tLayer[2] = { 0, 0 };
	    if (sLayer == 2) {
              tLayer[0] = 1; 
              tLayer[1] = 0;
            }
	    if (sLayer == 1) {
              tLayer[0] = 5; 
              tLayer[1] = 0;
            }
	    if (sLayer == 0) {
              tLayer[0] = 4; 
              tLayer[1] = 5;
            }

	    for (Int_t iLayer = 0; iLayer < 2; iLayer++) {
	      Int_t jLayer = tLayer[iLayer]; // Set tracking layer	      
	      cseed[jLayer].Reset();
	      cseed[jLayer].SetTilt(hL[jLayer]);
	      cseed[jLayer].SetPadLength(padlength[jLayer]);
	      cseed[jLayer].SetX0(xcl[jLayer]);
	      // Get pad length and rough cluster
	      Int_t indexdummy = reflayers[jLayer]->FindNearestCluster(cseed[jLayer].GetYref(0)
								      ,cseed[jLayer].GetZref(0)
                                                                      ,kRoad2y
                                                                      ,kRoad2z);
	      if (indexdummy <= 0) {
                continue; 
	      }
	      AliTRDcluster *cldummy = (AliTRDcluster *) GetCluster(indexdummy);
	      padlength[jLayer]      = TMath::Sqrt(cldummy->GetSigmaZ2() * 12.0);
	    }
	    AliTRDseed::FitRiemanTilt(cseed,kTRUE);

	    for (Int_t iLayer = 0; iLayer < 2; iLayer++) {

	      Int_t    jLayer = tLayer[iLayer];      // set tracking layer	
	      if ((jLayer == 0) && !(cseed[1].IsOK())) {
                continue;  // break not allowed
	      }
	      if ((jLayer == 5) && !(cseed[4].IsOK())) {
                continue;  // break not allowed
	      }
	      Float_t  zexp   = cseed[jLayer].GetZref(0);
	      Double_t zroad  = padlength[jLayer] * 0.5 + 1.0;

	      for (Int_t iter = 0; iter < 2; iter++) {

		AliTRDseed tseed = cseed[jLayer];
		Float_t quality = 10000.0;

		for (Int_t iTime = 2; iTime < 20; iTime++) { 
		  AliTRDpropagationLayer &layer = *(fTrSec[ns]->GetLayer(layers[jLayer][1]-iTime));
		  Double_t dxlayer = layer.GetX()-xcl[jLayer];
		  Double_t yexp    = tseed.GetYref(0) + tseed.GetYref(1) * dxlayer;
		  Float_t  yroad   = kRoad1y;
		  Int_t    index   = layer.FindNearestCluster(yexp,zexp,yroad,zroad);
		  if (index <= 0) {
                    continue; 
		  }
		  AliTRDcluster *cl = (AliTRDcluster *) GetCluster(index);	      
		  tseed.SetIndexes(iTime,index);
		  tseed.SetClusters(iTime,cl);  // Register cluster
		  tseed.SetX(iTime,dxlayer);    // Register cluster
		  tseed.SetY(iTime,cl->GetY()); // Register cluster
		  tseed.SetZ(iTime,cl->GetZ()); // Register cluster
		}

		tseed.Update();
		if (tseed.IsOK()) {
		  Float_t dangle   = tseed.GetYfit(1) - tseed.GetYref(1);
		  Float_t tquality = (18.0 - tseed.GetN2())/2.0 + TMath::Abs(dangle) / 0.1
                                   + TMath::Abs(tseed.GetYfit(0) - tseed.GetYref(0)) / 0.2
                                   + 2.0 * TMath::Abs(tseed.GetMeanz() - tseed.GetZref(0)) / padlength[jLayer];
		  if (tquality < quality) {
		    cseed[jLayer] = tseed;
		    quality       = tquality;
		  }
		}

		zroad *= 2.0;

	      } // Loop: iter

	      if ( cseed[jLayer].IsOK()) {
		cseed[jLayer].CookLabels();
		cseed[jLayer].UpdateUsed();
		nusedf += cseed[jLayer].GetNUsed();
		AliTRDseed::FitRiemanTilt(cseed,kTRUE);
	      }

	    } // Loop: iLayer

	    // Make copy
	    AliTRDseed bseed[6];
	    for (Int_t jLayer = 0; jLayer < 6; jLayer++) {
	      bseed[jLayer] = cseed[jLayer];
	    }
	    Float_t lastquality = 10000.0;
	    Float_t lastchi2    = 10000.0;
	    Float_t chi2        =  1000.0;

	    for (Int_t iter = 0; iter < 4; iter++) {

	      // Sort tracklets according "quality", try to "improve" 4 worst 
	      Float_t sumquality = 0.0;
	      Float_t squality[6];
	      Int_t   sortindexes[6];

	      for (Int_t jLayer = 0; jLayer < 6; jLayer++) {
		if (bseed[jLayer].IsOK()) { 
		  AliTRDseed &tseed = bseed[jLayer];
		  Double_t zcor     = tseed.GetTilt() * (tseed.GetZProb() - tseed.GetZref(0));
		  Float_t  dangle   = tseed.GetYfit(1) - tseed.GetYref(1);
		  Float_t  tquality = (18.0 - tseed.GetN2()) / 2.0 + TMath::Abs(dangle)        / 0.1
                                    + TMath::Abs(tseed.GetYfit(0) - (tseed.GetYref(0) - zcor)) / 0.2
                                    + 2.0 * TMath::Abs(tseed.GetMeanz() - tseed.GetZref(0))    / padlength[jLayer];
		  squality[jLayer]  = tquality;
		}
		else {
                  squality[jLayer] = -1.0;
		}
		sumquality +=squality[jLayer];
	      }

	      if ((sumquality >= lastquality) ||  
                  (chi2       >     lastchi2)) {
                break;
	      }
	      lastquality = sumquality;	 
	      lastchi2    = chi2;
	      if (iter > 0) {
		for (Int_t jLayer = 0; jLayer < 6; jLayer++) {
		  cseed[jLayer] = bseed[jLayer];
		}
	      }
	      TMath::Sort(6,squality,sortindexes,kFALSE);

	      for (Int_t jLayer = 5; jLayer > 1; jLayer--) {

		Int_t bLayer = sortindexes[jLayer];
		AliTRDseed tseed = bseed[bLayer];

		for (Int_t iTime = 2; iTime < 20; iTime++) { 

		  AliTRDpropagationLayer &layer = *(fTrSec[ns]->GetLayer(layers[bLayer][1]-iTime));
		  Double_t dxlayer = layer.GetX() - xcl[bLayer];
		  Double_t zexp    = tseed.GetZref(0);
		  Double_t zcor    = tseed.GetTilt() * (tseed.GetZProb() - tseed.GetZref(0));
		  Float_t  roadz   = padlength[bLayer] + 1;
 		  if (TMath::Abs(tseed.GetZProb() - zexp) > 0.5*padlength[bLayer]) {
                    roadz = padlength[bLayer] * 0.5;
                  }
		  if (tseed.GetZfit(1)*tseed.GetZref(1)   < 0.0) {
                    roadz = padlength[bLayer] * 0.5;
                  }
		  if (TMath::Abs(tseed.GetZProb() - zexp) < 0.1*padlength[bLayer]) {
		    zexp  = tseed.GetZProb(); 
		    roadz = padlength[bLayer] * 0.5;
		  }

		  Double_t yexp  = tseed.GetYref(0) + tseed.GetYref(1) * dxlayer - zcor;
		  Int_t    index = layer.FindNearestCluster(yexp,zexp,kRoad1y,roadz);
		  if (index <= 0) {
                    continue; 
		  }
		  AliTRDcluster *cl = (AliTRDcluster *) GetCluster(index);	      

		  tseed.SetIndexes(iTime,index);
		  tseed.SetClusters(iTime,cl);  // Register cluster
		  tseed.SetX(iTime,dxlayer);    // Register cluster
		  tseed.SetY(iTime,cl->GetY()); // Register cluster
		  tseed.SetZ(iTime,cl->GetZ()); // Register cluster

		}

		tseed.Update();
		if (tseed.IsOK()) {
		  Float_t  dangle   = tseed.GetYfit(1) - tseed.GetYref(1);
		  Double_t zcor     = tseed.GetTilt() * (tseed.GetZProb() - tseed.GetZref(0));
		  Float_t  tquality = (18.0 - tseed.GetN2()) / 2.0 
                                    + TMath::Abs(dangle)     / 0.1
                                    + TMath::Abs(tseed.GetYfit(0) - (tseed.GetYref(0) - zcor)) / 0.2
                                    + 2.0 * TMath::Abs(tseed.GetMeanz() - tseed.GetZref(0))    / padlength[jLayer];
		  if (tquality<squality[bLayer]) {
		    bseed[bLayer] = tseed;
		  }
		}

	      } // Loop: jLayer

	      chi2 = AliTRDseed::FitRiemanTilt(bseed,kTRUE);

	    } // Loop: iter

	    nclusters = 0;
	    nlayers   = 0;
	    findable  = 0;
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (TMath::Abs(cseed[iLayer].GetYref(0) / cseed[iLayer].GetX0()) < 0.15) {
		findable++;
	      }
	      if (cseed[iLayer].IsOK()) {
		nclusters += cseed[iLayer].GetN2();	    
		nlayers++;
	      }
	    }
	    if (nlayers < 3) {
              continue;
	    }
	    rieman.Reset();
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (cseed[iLayer].IsOK()) {
                rieman.AddPoint(xcl[iLayer]
			       ,cseed[iLayer].GetYfitR(0)
			       ,cseed[iLayer].GetZProb()
                               ,1
                               ,10);
	      }
	    }
	    rieman.Update();

	    chi2RF = 0.0;
	    chi2ZF = 0.0;
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (cseed[iLayer].IsOK()) {
		cseed[iLayer].SetYref(0,rieman.GetYat(xcl[iLayer]));
		chi2RF += (cseed[iLayer].GetYref(0) - cseed[iLayer].GetYfitR(0))
                        * (cseed[iLayer].GetYref(0) - cseed[iLayer].GetYfitR(0));
		cseed[iLayer].SetYref(1,rieman.GetDYat(xcl[iLayer]));
		cseed[iLayer].SetZref(0,rieman.GetZat(xcl[iLayer]));
		chi2ZF += (cseed[iLayer].GetZref(0) - cseed[iLayer].GetMeanz())
                        * (cseed[iLayer].GetZref(0) - cseed[iLayer].GetMeanz());
		cseed[iLayer].SetZref(1,rieman.GetDZat(xcl[iLayer]));
	      }
	    }
	    chi2RF /= TMath::Max((nlayers - 3.0),1.0);
	    chi2ZF /= TMath::Max((nlayers - 3.0),1.0);
	    curv    = rieman.GetC();

	    Double_t xref2 = (xcl[2] + xcl[3]) * 0.5; // Middle of the chamber
	    Double_t dzmf  = rieman.GetDZat(xref2);
	    Double_t zmf   = rieman.GetZat(xref2);

	    //
	    // Fit hyperplane
	    //
	    Int_t npointsT = 0;
	    fitterTC.ClearPoints();
	    fitterT2.ClearPoints();
	    rieman2.Reset();

	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {

	      if (!cseed[iLayer].IsOK()) {
                continue;
	      }

	      for (Int_t itime = 0; itime < 25; itime++) {

		if (!cseed[iLayer].IsUsable(itime)) {
                  continue;
		}
                // X relative to the middle chamber
		Double_t x  = cseed[iLayer].GetX(itime) + cseed[iLayer].GetX0() - xref2;  
		Double_t y  = cseed[iLayer].GetY(itime);
		Double_t z  = cseed[iLayer].GetZ(itime);
		// ExB correction to the correction
		// Tilted rieman
		Double_t uvt[6];
                // Global x
		Double_t x2 = cseed[iLayer].GetX(itime) + cseed[iLayer].GetX0();      
		Double_t t  = 1.0 / (x2*x2 + y*y);
		uvt[1] = t;                 // t
		uvt[0] = 2.0 * x2 * uvt[1]; // u 
		uvt[2] = 2.0 * hL[iLayer] * uvt[1];
		uvt[3] = 2.0 * hL[iLayer] * x * uvt[1];	      
		uvt[4] = 2.0 * (y + hL[iLayer]*z) * uvt[1];
		Double_t error = 2.0 * 0.2 * uvt[1];
		fitterT2.AddPoint(uvt,uvt[4],error);

		//
		// Constrained rieman
		// 
		z = cseed[iLayer].GetZ(itime);
		uvt[0] = 2.0 * x2 * t; // u 
		uvt[1] = 2.0 * hL[iLayer] * x2 * uvt[1];	      
		uvt[2] = 2.0 * (y + hL[iLayer] * (z - GetZ())) * t;
		fitterTC.AddPoint(uvt,uvt[2],error);
		rieman2.AddPoint(x2,y,z,1,10);
		npointsT++;

	      }

	    } // Loop: iLayer

	    rieman2.Update();
	    fitterTC.Eval();
	    fitterT2.Eval();
	    Double_t rpolz0 = fitterT2.GetParameter(3);
	    Double_t rpolz1 = fitterT2.GetParameter(4);	    

	    //
	    // Linear fitter  - not possible to make boundaries
	    // Do not accept non possible z and dzdx combinations
	    // 	    
	    Bool_t acceptablez = kTRUE;
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (cseed[iLayer].IsOK()) {
		Double_t zT2 = rpolz0 + rpolz1 * (xcl[iLayer] - xref2);
		if (TMath::Abs(cseed[iLayer].GetZProb() - zT2) > padlength[iLayer] * 0.5 + 1.0) {
		  acceptablez = kFALSE;
		}
	      }
	    }
	    if (!acceptablez) {
	      fitterT2.FixParameter(3,zmf);
	      fitterT2.FixParameter(4,dzmf);
	      fitterT2.Eval();
	      fitterT2.ReleaseParameter(3);
	      fitterT2.ReleaseParameter(4);
	      rpolz0 = fitterT2.GetParameter(3);
	      rpolz1 = fitterT2.GetParameter(4);
	    }

	    Double_t chi2TR = fitterT2.GetChisquare() / Float_t(npointsT);
	    Double_t chi2TC = fitterTC.GetChisquare() / Float_t(npointsT);
	    Double_t polz1c = fitterTC.GetParameter(2);
	    Double_t polz0c = polz1c * xref2;
	    Double_t aC     =  fitterTC.GetParameter(0);
	    Double_t bC     =  fitterTC.GetParameter(1);
	    Double_t cC     =  aC / TMath::Sqrt(bC * bC + 1.0); // Curvature
	    Double_t aR     =  fitterT2.GetParameter(0);
	    Double_t bR     =  fitterT2.GetParameter(1);
	    Double_t dR     =  fitterT2.GetParameter(2);	    
	    Double_t cR     =  1.0 + bR*bR - dR*aR;
	    Double_t dca    =  0.0;	    
	    if (cR > 0.0) {
	      dca = -dR / (TMath::Sqrt(1.0 + bR*bR - dR*aR) + TMath::Sqrt(1.0 + bR*bR)); 
	      cR  =  aR / TMath::Sqrt(cR);
	    }

	    Double_t chi2ZT2 = 0.0;
            Double_t chi2ZTC = 0.0;
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (cseed[iLayer].IsOK()) {
		Double_t zT2 = rpolz0 + rpolz1 * (xcl[iLayer] - xref2);
		Double_t zTC = polz0c + polz1c * (xcl[iLayer] - xref2);
		chi2ZT2 += TMath::Abs(cseed[iLayer].GetMeanz() - zT2);
		chi2ZTC += TMath::Abs(cseed[iLayer].GetMeanz() - zTC);
	      }
	    }
	    chi2ZT2 /= TMath::Max((nlayers - 3.0),1.0);
	    chi2ZTC /= TMath::Max((nlayers - 3.0),1.0);	    

	    AliTRDseed::FitRiemanTilt(cseed,kTRUE);
	    Float_t sumdaf = 0.0;
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (cseed[iLayer].IsOK()) {
		sumdaf += TMath::Abs((cseed[iLayer].GetYfit(1) - cseed[iLayer].GetYref(1))
                                    / cseed[iLayer].GetSigmaY2());
	      }
	    }
	    sumdaf /= Float_t (nlayers - 2.0);

	    //
	    // Likelihoods for full track
	    //
	    Double_t likezf     = TMath::Exp(-chi2ZF * 0.14);
	    Double_t likechi2C  = TMath::Exp(-chi2TC * 0.677);
	    Double_t likechi2TR = TMath::Exp(-chi2TR * 0.78);
	    Double_t likeaf     = TMath::Exp(-sumdaf * 3.23);
	    seedquality2[registered] = likezf * likechi2TR * likeaf; 

	    // Still needed ????
//	    Bool_t isGold = kFALSE;
//	    
// 	    if (nlayers == 6        && TMath::Log(0.000000001+seedquality2[index])<-5.) isGold =kTRUE;   // gold
// 	    if (nlayers == findable && TMath::Log(0.000000001+seedquality2[index])<-4.) isGold =kTRUE;   // gold
// 	    if (isGold &&nusedf<10){
// 	      for (Int_t jLayer=0;jLayer<6;jLayer++){
// 		if ( seed[index][jLayer].IsOK()&&TMath::Abs(seed[index][jLayer].fYfit[1]-seed[index][jLayer].fYfit[1])<0.1)
// 		  seed[index][jLayer].UseClusters();  //sign gold
// 	      }
// 	    }

	    Int_t index0 = 0;
	    if (!cseed[0].IsOK()) {
	      index0 = 1;
	      if (!cseed[1].IsOK()) {
                index0 = 2;
	      }
	    }
	    seedparams[registered][0] = cseed[index0].GetX0();
	    seedparams[registered][1] = cseed[index0].GetYref(0);
	    seedparams[registered][2] = cseed[index0].GetZref(0);
	    seedparams[registered][5] = cR;
	    seedparams[registered][3] = cseed[index0].GetX0() * cR - TMath::Sin(TMath::ATan(cseed[0].GetYref(1)));
	    seedparams[registered][4] = cseed[index0].GetZref(1)
                                      /	TMath::Sqrt(1.0 + cseed[index0].GetYref(1) * cseed[index0].GetYref(1));
	    seedparams[registered][6] = ns;

	    Int_t labels[12];
            Int_t outlab[24];
	    Int_t nlab = 0;
	    for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	      if (!cseed[iLayer].IsOK()) {
                continue;
	      }
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
	      cseed[iLayer].SetC(cR);
	      cseed[iLayer].SetCC(cC);
	      cseed[iLayer].SetChi2(chi2TR);
	      cseed[iLayer].SetChi2Z(chi2ZF);
	    }

	    // Debugging print
	    if (1 || (!isFake)) {
	      Float_t zvertex = GetZ();
	      TTreeSRedirector &cstream = *fDebugStreamer;
	      if (AliTRDReconstructor::StreamLevel() > 0) {
	        cstream << "Seeds1"
			<< "isFake="     << isFake
			<< "Vertex="     << zvertex
			<< "Rieman2.="   << &rieman2
			<< "Rieman.="    << &rieman
			<< "Xref="       << xref
			<< "X0="         << xcl[0]
			<< "X1="         << xcl[1]
			<< "X2="         << xcl[2]
			<< "X3="         << xcl[3]
			<< "X4="         << xcl[4]
			<< "X5="         << xcl[5]
			<< "Chi2R="      << chi2R
			<< "Chi2Z="      << chi2Z
			<< "Chi2RF="     << chi2RF  // Chi2 of trackletes on full track
			<< "Chi2ZF="     << chi2ZF  // Chi2 z on tracklets on full track
			<< "Chi2ZT2="    << chi2ZT2 // Chi2 z on tracklets on full track  - rieman tilt
			<< "Chi2ZTC="    << chi2ZTC // Chi2 z on tracklets on full track  - rieman tilt const
			<< "Chi2TR="     << chi2TR  // Chi2 without vertex constrain
			<< "Chi2TC="     << chi2TC  // Chi2 with    vertex constrain
			<< "C="          << curv    // Non constrained - no tilt correction
			<< "DR="         << dR      // DR parameter    - tilt correction
			<< "DCA="        << dca     // DCA             - tilt correction
			<< "CR="         << cR      // Non constrained curvature - tilt correction
			<< "CC="         << cC      // Constrained curvature
			<< "Polz0="      << polz0c
			<< "Polz1="      << polz1c
			<< "RPolz0="     << rpolz0
			<< "RPolz1="     << rpolz1
			<< "Ncl="        << nclusters
			<< "Nlayers="    << nlayers
			<< "NUsedS="     << nusedCl
			<< "NUsed="      << nusedf
			<< "Findable="   << findable
			<< "Like="       << like
			<< "LikePrim="   << likePrim
			<< "Likechi2C="  << likechi2C
			<< "Likechi2TR=" << likechi2TR
			<< "Likezf="     << likezf
			<< "LikeF="      << seedquality2[registered]
			<< "S0.="        << &cseed[0]
			<< "S1.="        << &cseed[1]
			<< "S2.="        << &cseed[2]
			<< "S3.="        << &cseed[3]
			<< "S4.="        << &cseed[4]
			<< "S5.="        << &cseed[5]
			<< "SB0.="       << &seedb[0]
			<< "SB1.="       << &seedb[1]
			<< "SB2.="       << &seedb[2]
			<< "SB3.="       << &seedb[3]
			<< "SB4.="       << &seedb[4]
			<< "SB5.="       << &seedb[5]
			<< "Label="      << label
			<< "Freq="       << frequency
			<< "sLayer="     << sLayer
			<< "\n";
	      }
	    }

	    if (registered<kMaxSeed - 1) {
	      registered++;
	      cseed = seed[registered];
	    }

	  } // End of loop over layer 1

	} // End of loop over layer 0 

      } // End of loop over layer 3     

    } // End of loop over seeding time bins 

    //
    // Choose best
    //

    TMath::Sort(registered,seedquality2,sort,kTRUE);
    Bool_t signedseed[kMaxSeed];
    for (Int_t i = 0; i < registered; i++) {
      signedseed[i] = kFALSE;
    }

    for (Int_t iter = 0; iter < 5; iter++) {

      for (Int_t iseed = 0; iseed < registered; iseed++) {
      
	Int_t index = sort[iseed];
	if (signedseed[index]) {
          continue;
	}
	Int_t labelsall[1000];
	Int_t nlabelsall = 0;
	Int_t naccepted  = 0;;
	Int_t sLayer     = seedlayer[index];
	Int_t ncl        = 0;
	Int_t nused      = 0;
	Int_t nlayers    = 0;
	Int_t findable   = 0;

	for (Int_t jLayer = 0; jLayer < 6; jLayer++) {

	  if (TMath::Abs(seed[index][jLayer].GetYref(0) / xcl[jLayer]) < 0.15) {
	    findable++;
	  }
	  if (seed[index][jLayer].IsOK()) {
	    seed[index][jLayer].UpdateUsed();
	    ncl   +=seed[index][jLayer].GetN2();
	    nused +=seed[index][jLayer].GetNUsed();
	    nlayers++;
	    // Cooking label
	    for (Int_t itime = 0; itime < 25; itime++) {
	      if (seed[index][jLayer].IsUsable(itime)) {
		naccepted++;
		for (Int_t ilab = 0; ilab < 3; ilab++) {
		  Int_t tindex = seed[index][jLayer].GetClusters(itime)->GetLabel(ilab);
		  if (tindex >= 0) {
		    labelsall[nlabelsall] = tindex;
		    nlabelsall++;
		  }
		}
	      }
	    }
	  }

	}

	if (nused > 30) {
          continue;
	}

 	if (iter == 0) {
	  if (nlayers < 6) {
            continue;
	  }
	  if (TMath::Log(0.000000001+seedquality2[index]) < -5.0) {
            continue; // Gold
	  }
	}

	if (iter == 1) {
	  if (nlayers < findable) {
            continue;
	  }
	  if (TMath::Log(0.000000001+seedquality2[index]) < -4.0) {
            continue;
	  }
	}

	if (iter == 2) {
	  if ((nlayers == findable) || 
              (nlayers ==        6)) {
            continue;
	  }
	  if (TMath::Log(0.000000001+seedquality2[index]) < -6.0) {
            continue;
	  }
	}

	if (iter == 3) {
	  if (TMath::Log(0.000000001+seedquality2[index]) < -5.0) {
            continue;
	  }
	}

	if (iter == 4) {
	  if (TMath::Log(0.000000001+seedquality2[index]) - nused/(nlayers-3.0) < -15.0) {
            continue;
	  }
	}

	signedseed[index] = kTRUE;

	Int_t labels[1000];
        Int_t outlab[1000];
	Int_t nlab = 0;
	for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	  if (seed[index][iLayer].IsOK()) {
	    if (seed[index][iLayer].GetLabels(0) >= 0) {
	      labels[nlab] = seed[index][iLayer].GetLabels(0);
	      nlab++;
	    }
	    if (seed[index][iLayer].GetLabels(1) >= 0) {
	      labels[nlab] = seed[index][iLayer].GetLabels(1);
	      nlab++;
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
	Float_t ratio     = Float_t(nused) / Float_t(ncl);
	if (ratio < 0.25) {
	  for (Int_t jLayer = 0; jLayer < 6; jLayer++) {
	    if ((seed[index][jLayer].IsOK()) && 
                (TMath::Abs(seed[index][jLayer].GetYfit(1) - seed[index][jLayer].GetYfit(1)) < 0.2)) {
	      seed[index][jLayer].UseClusters(); // Sign gold
	    }
	  }
	}

	Int_t eventNr = esd->GetEventNumber();
	TTreeSRedirector &cstream = *fDebugStreamer;

	//
	// Register seed
	//
	AliTRDtrack *track = RegisterSeed(seed[index],seedparams[index]);
	AliTRDtrack  dummy;
	if (!track) {
          track = &dummy;
	}
	else {
	  AliESDtrack esdtrack;
	  esdtrack.UpdateTrackParams(track,AliESDtrack::kTRDout);
	  esdtrack.SetLabel(label);
	  esd->AddTrack(&esdtrack);	
	  TTreeSRedirector &cstream = *fDebugStreamer;
	  if (AliTRDReconstructor::StreamLevel() > 0) {
	    cstream << "Tracks"
		    << "EventNr="  << eventNr
		    << "ESD.="     << &esdtrack
		    << "trd.="     << track
		    << "trdback.=" << track
		    << "\n";
	  }
	}

	if (AliTRDReconstructor::StreamLevel() > 0) {
	  cstream << "Seeds2"
		  << "Iter="      << iter
		  << "Track.="    << track
		  << "Like="      << seedquality[index]
		  << "LikeF="     << seedquality2[index]
		  << "S0.="       << &seed[index][0]
		  << "S1.="       << &seed[index][1]
		  << "S2.="       << &seed[index][2]
		  << "S3.="       << &seed[index][3]
		  << "S4.="       << &seed[index][4]
		  << "S5.="       << &seed[index][5]
		  << "Label="     << label
		  << "Label1="    << label1
		  << "Label2="    << label2
		  << "FakeRatio=" << fakeratio
		  << "Freq="      << frequency
		  << "Ncl="       << ncl
		  << "Nlayers="   << nlayers
		  << "Findable="  << findable
		  << "NUsed="     << nused
		  << "sLayer="    << sLayer
		  << "EventNr="   << eventNr
		  << "\n";
	}

      } // Loop: iseed

    } // Loop: iter

  } // End of loop over sectors

  delete [] pseed;

}
          
//_____________________________________________________________________________
Int_t AliTRDtracker::ReadClusters(TObjArray *array, TTree *clusterTree) const
{
  //
  // Reads AliTRDclusters (option >= 0) or AliTRDrecPoints (option < 0) 
  // from the file. The names of the cluster tree and branches 
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
  
  // Loop through all entries in the tree
  Int_t nEntries   = (Int_t) clusterTree->GetEntries();
  Int_t nbytes     = 0;
  AliTRDcluster *c = 0;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    
    
    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);  
    
    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
    
    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (AliTRDcluster *) clusterArray->UncheckedAt(iCluster);
      AliTRDcluster *co = c;
      array->AddLast(co);
      clusterArray->RemoveAt(iCluster); 
    }

  }

  delete clusterArray;

  return 0;

}

//_____________________________________________________________________________
Bool_t AliTRDtracker::GetTrackPoint(Int_t index, AliTrackPoint &p) const
{
  //
  // Get track space point with index i
  // Origin: C.Cheshkov
  //

  AliTRDcluster *cl = (AliTRDcluster *) fClusters->UncheckedAt(index);
  Int_t  idet     = cl->GetDetector();
  Int_t  isector  = fGeom->GetSector(idet);
  Int_t  ichamber = fGeom->GetChamber(idet);
  Int_t  iplan    = fGeom->GetPlane(idet);
  Double_t local[3];
  local[0] = GetX(isector,iplan,cl->GetLocalTimeBin());
  local[1] = cl->GetY();
  local[2] = cl->GetZ();
  Double_t global[3];
  fGeom->RotateBack(idet,local,global);
  p.SetXYZ(global[0],global[1],global[2]);
  AliAlignObj::ELayerID iLayer = AliAlignObj::kTRD1;
  switch (iplan) {
  case 0:
    iLayer = AliAlignObj::kTRD1;
    break;
  case 1:
    iLayer = AliAlignObj::kTRD2;
    break;
  case 2:
    iLayer = AliAlignObj::kTRD3;
    break;
  case 3:
    iLayer = AliAlignObj::kTRD4;
    break;
  case 4:
    iLayer = AliAlignObj::kTRD5;
    break;
  case 5:
    iLayer = AliAlignObj::kTRD6;
    break;
  };
  Int_t    modId = isector * fGeom->Ncham() + ichamber;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,modId);
  p.SetVolumeID(volid);

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDtracker::CookLabel(AliKalmanTrack *pt, Float_t wrong) const 
{
  //
  // This cooks a label. Mmmmh, smells good...
  //

  Int_t label = 123456789;
  Int_t index;
  Int_t i;
  Int_t j;
  Int_t ncl   = pt->GetNumberOfClusters();

  const Int_t kRange = fTrSec[0]->GetOuterTimeBin() + 1;

  Bool_t labelAdded;

  Int_t **s = new Int_t* [kRange];
  for (i = 0; i < kRange; i++) {
    s[i] = new Int_t[2];
  }
  for (i = 0; i < kRange; i++) {
    s[i][0] = -1;
    s[i][1] =  0;
  }

  Int_t t0;
  Int_t t1;
  Int_t t2;

  for (i = 0; i < ncl; i++) {
    index = pt->GetClusterIndex(i);
    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(index);
    t0=c->GetLabel(0);
    t1=c->GetLabel(1);
    t2=c->GetLabel(2);
  }

  for (i = 0; i < ncl; i++) {
    index = pt->GetClusterIndex(i);
    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(index);
    for (Int_t k = 0; k < 3; k++) { 
      label      = c->GetLabel(k);
      labelAdded = kFALSE; 
      j = 0;
      if (label >= 0) {
        while ((!labelAdded) && (j < kRange)) {
          if ((s[j][0] == label) || 
              (s[j][1] ==     0)) {
            s[j][0] = label; 
            s[j][1] = s[j][1] + 1; 
            labelAdded = kTRUE;
          }
          j++;
        }
      }
    }
  }

  Int_t max = 0;
  label = -123456789;

  for (i = 0; i < kRange; i++) {
    if (s[i][1] > max) {
      max   = s[i][1]; 
      label = s[i][0];
    }
  }

  for (i = 0; i < kRange; i++) {
    delete []s[i];
  }        

  delete []s;

  if ((1.0 - Float_t(max)/ncl) > wrong) {
    label = -label;   
  }

  pt->SetLabel(label); 

}

//_____________________________________________________________________________
void AliTRDtracker::UseClusters(const AliKalmanTrack *t, Int_t from) const 
{
  //
  // Use clusters, but don't abuse them!
  //

  const Float_t kmaxchi2 = 18;
  const Float_t kmincl   = 10;

  AliTRDtrack *track = (AliTRDtrack *) t;

  Int_t ncl = t->GetNumberOfClusters();
  for (Int_t i = from; i < ncl; i++) {
    Int_t index  = t->GetClusterIndex(i);
    AliTRDcluster *c = (AliTRDcluster *) fClusters->UncheckedAt(index);
    Int_t iplane = fGeom->GetPlane(c->GetDetector());
    if (track->GetTracklets(iplane).GetChi2() > kmaxchi2) {
      continue; 
    }
    if (track->GetTracklets(iplane).GetN()    <   kmincl) {
      continue; 
    }
    if (!(c->IsUsed())) {
      c->Use();
    }
  }

}

//_____________________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaY2(Double_t , Double_t , Double_t ) const
{
  //
  // Parametrised "expected" error of the cluster reconstruction in Y 
  //

  Double_t s = 0.08 * 0.08;    
  return s;

}

//_____________________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaZ2(Double_t , Double_t ) const
{
  //
  // Parametrised "expected" error of the cluster reconstruction in Z 
  //

  Double_t s = 9.0 * 9.0 / 12.0;  
  return s;

}                  

//_____________________________________________________________________________
Double_t AliTRDtracker::GetX(Int_t sector, Int_t plane, Int_t localTB) const 
{
  //
  // Returns radial position which corresponds to time bin <localTB>
  // in tracking sector <sector> and plane <plane>
  //

  Int_t index = fTrSec[sector]->CookTimeBinIndex(plane, localTB); 
  Int_t pl    = fTrSec[sector]->GetLayerNumber(index);

  return fTrSec[sector]->GetLayer(pl)->GetX();

}

//_____________________________________________________________________________
AliTRDtracker::AliTRDpropagationLayer
             ::AliTRDpropagationLayer(Double_t x, Double_t dx, Double_t rho
                                    , Double_t radLength, Int_t tbIndex, Int_t plane)
  :fN(0)
  ,fSec(0)
  ,fClusters(NULL)
  ,fIndex(NULL)
  ,fX(x)
  ,fdX(dx)
  ,fRho(rho)
  ,fX0(radLength)
  ,fTimeBinIndex(tbIndex)
  ,fPlane(plane)
  ,fYmax(0)
  ,fYmaxSensitive(0)
  ,fHole(kFALSE)
  ,fHoleZc(0)
  ,fHoleZmax(0)
  ,fHoleYc(0)
  ,fHoleYmax(0)
  ,fHoleRho(0)
  ,fHoleX0(0)
{ 
  //
  // AliTRDpropagationLayer constructor
  //

  for (Int_t i = 0; i < (Int_t) kZones; i++) {
    fZc[i]   = 0; 
    fZmax[i] = 0;
  }

  if (fTimeBinIndex >= 0) { 
    fClusters = new AliTRDcluster*[kMaxClusterPerTimeBin];
    fIndex    = new UInt_t[kMaxClusterPerTimeBin];
  }

  for (Int_t i = 0; i < 5; i++) {
    fIsHole[i] = kFALSE;
  }

}

//_____________________________________________________________________________
void AliTRDtracker::AliTRDpropagationLayer
                  ::SetHole(Double_t Zmax, Double_t Ymax, Double_t rho
                          , Double_t radLength, Double_t Yc, Double_t Zc) 
{
  //
  // Sets hole in the layer 
  //

  fHole     = kTRUE;
  fHoleZc   = Zc;
  fHoleZmax = Zmax;
  fHoleYc   = Yc;
  fHoleYmax = Ymax;
  fHoleRho  = rho;
  fHoleX0   = radLength;

}

//_____________________________________________________________________________
AliTRDtracker::AliTRDtrackingSector
             ::AliTRDtrackingSector(AliTRDgeometry *geo, Int_t gs)
  :fN(0)
  ,fGeom(geo)
  ,fGeomSector(gs)
{
  //
  // AliTRDtrackingSector Constructor
  //

  AliTRDpadPlane *padPlane = 0;
  AliTRDpropagationLayer *ppl = 0;

  // Get holes description from geometry
  Bool_t holes[AliTRDgeometry::kNcham];
  for (Int_t icham = 0; icham < AliTRDgeometry::kNcham; icham++) {
    holes[icham] = fGeom->IsHole(0,icham,gs);
  } 
  
  for (UInt_t i = 0; i < kMaxTimeBinIndex; i++) {
    fTimeBinIndex[i] = -1;
  }

  Double_t x;
  Double_t dx;
  Double_t rho;
  Double_t radLength;

  // Add layers for each of the planes
  Double_t dxAmp   = (Double_t) fGeom->CamHght(); // Amplification region
  //Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region  

  const Int_t kNchambers  = AliTRDgeometry::Ncham();
  Int_t    tbIndex;
  Double_t ymax           = 0;
  Double_t ymaxsensitive  = 0;
  Double_t *zc            = new Double_t[kNchambers];
  Double_t *zmax          = new Double_t[kNchambers];
  Double_t *zmaxsensitive = new Double_t[kNchambers];  

  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliErrorGeneral("AliTRDtrackingSector::Ctor"
                   ,"Could not get common parameters\n");
    return;
  }
    
  for (Int_t plane = 0; plane < AliTRDgeometry::Nplan(); plane++) {

    ymax          = fGeom->GetChamberWidth(plane) / 2.0;
    padPlane      = commonParam->GetPadPlane(plane,0);
    ymaxsensitive = (padPlane->GetColSize(1) * padPlane->GetNcols() - 4.0) / 2.0;    

    for (Int_t ch = 0; ch < kNchambers; ch++) {
      zmax[ch]          = fGeom->GetChamberLength(plane,ch) / 2.0;
      Float_t pad       = padPlane->GetRowSize(1);
      Float_t row0      = commonParam->GetRow0(plane,ch,0);
      Int_t   nPads     = commonParam->GetRowMax(plane,ch,0);
      zmaxsensitive[ch] = Float_t(nPads) * pad / 2.0;      
      zc[ch]            = -(pad * nPads) / 2.0 + row0;
    }

    dx        = AliTRDcalibDB::Instance()->GetVdrift(0,0,0)
              / commonParam->GetSamplingFrequency();
    rho       = 0.00295 * 0.85; //????
    radLength = 11.0;  

    Double_t x0 = (Double_t) AliTRDgeometry::GetTime0(plane);
    //Double_t xbottom = x0 - dxDrift;
    //Double_t xtop    = x0 + dxAmp;

    Int_t nTimeBins =  AliTRDcalibDB::Instance()->GetNumberOfTimeBins();    
    for (Int_t iTime = 0; iTime < nTimeBins; iTime++) {

      Double_t xlayer = iTime * dx - dxAmp;
      //if (xlayer<0) xlayer = dxAmp / 2.0;
      x = x0 - xlayer;

      tbIndex = CookTimeBinIndex(plane,iTime);
      ppl     = new AliTRDpropagationLayer(x,dx,rho,radLength,tbIndex,plane);
      ppl->SetYmax(ymax,ymaxsensitive);
      ppl->SetZ(zc,zmax,zmaxsensitive);
      ppl->SetHoles(holes);
      InsertLayer(ppl);      

    }

  }    

  MapTimeBinLayers();

  delete [] zc;
  delete [] zmax;
  delete [] zmaxsensitive;

}

//_____________________________________________________________________________
AliTRDtracker::AliTRDtrackingSector
             ::AliTRDtrackingSector(const AliTRDtrackingSector &/*t*/)
  :fN(0)
  ,fGeom(0)
  ,fGeomSector(0)
{
  //
  // Copy constructor
  //

}

//_____________________________________________________________________________
Int_t  AliTRDtracker::AliTRDtrackingSector
                    ::CookTimeBinIndex(Int_t plane, Int_t localTB) const
{
  //
  // depending on the digitization parameters calculates "global"
  // time bin index for timebin <localTB> in plane <plane>
  //
  //

  Int_t tbPerPlane = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t gtb        = (plane+1) * tbPerPlane - localTB - 1;
  if (localTB < 0) {
    return -1;
  }
  if (gtb     < 0) {
    return -1;
  }

  return gtb;

}

//_____________________________________________________________________________
void AliTRDtracker::AliTRDtrackingSector
                  ::MapTimeBinLayers() 
{
  //
  // For all sensitive time bins sets corresponding layer index
  // in the array fTimeBins 
  //

  Int_t index;

  for (Int_t i = 0; i < fN; i++) {

    index = fLayers[i]->GetTimeBinIndex();
    
    if (index < 0) {
      continue;
    }
    if (index >= (Int_t) kMaxTimeBinIndex) {
      //AliWarning(Form("Index %d exceeds allowed maximum of %d!\n"
      //               ,index,kMaxTimeBinIndex-1));
      continue;
    }

    fTimeBinIndex[index] = i;

  }

}

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDtrackingSector
                   ::GetLayerNumber(Double_t x) const
{
  // 
  // Returns the number of time bin which in radial position is closest to <x>
  //

  if (x >= fLayers[fN-1]->GetX()) {
    return fN - 1; 
  }
  if (x <= fLayers[   0]->GetX()) {
    return 0; 
  }

  Int_t b = 0;
  Int_t e = fN - 1;
  Int_t m = (b + e) / 2;

  for ( ; b < e; m = (b + e) / 2) {
    if (x > fLayers[m]->GetX()) {
      b = m + 1;
    }
    else {
      e = m;
    }
  }

  if (TMath::Abs(x - fLayers[m]->GetX()) > TMath::Abs(x - fLayers[m+1]->GetX())) {
    return m + 1;
  }
  else {
    return m;
  }

}

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDtrackingSector
                   ::GetInnerTimeBin() const 
{
  // 
  // Returns number of the innermost SENSITIVE propagation layer
  //

  return GetLayerNumber(0);

}

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDtrackingSector
                   ::GetOuterTimeBin() const 
{
  // 
  // Returns number of the outermost SENSITIVE time bin
  //

  return GetLayerNumber(GetNumberOfTimeBins() - 1);

}

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDtrackingSector
                   ::GetNumberOfTimeBins() const 
{
  // 
  // Returns number of SENSITIVE time bins
  //

  Int_t tb;
  Int_t layer;

  for (tb = kMaxTimeBinIndex - 1; tb >= 0; tb--) {
    layer = GetLayerNumber(tb);
    if (layer >= 0) {
      break;
    }
  }

  return tb + 1;

}

//_____________________________________________________________________________
void AliTRDtracker::AliTRDtrackingSector
                  ::InsertLayer(AliTRDpropagationLayer *pl)
{ 
  //
  // Insert layer <pl> in fLayers array.
  // Layers are sorted according to X coordinate.
  //

  if (fN == ((Int_t) kMaxLayersPerSector)) {
    //AliWarning("Too many layers !\n");
    return;
  }

  if (fN == 0) {
    fLayers[fN++] = pl; 
    return;
  }

  Int_t i = Find(pl->GetX());

  memmove(fLayers+i+1,fLayers+i,(fN-i)*sizeof(AliTRDpropagationLayer*));

  fLayers[i] = pl; 
  fN++;

}              

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDtrackingSector
                   ::Find(Double_t x) const 
{
  //
  // Returns index of the propagation layer nearest to X 
  //

  if (x <= fLayers[0]->GetX()) {
    return 0;
  }

  if (x > fLayers[fN-1]->GetX()) {
    return fN;
  }

  Int_t b = 0;
  Int_t e = fN-1;
  Int_t m = (b + e) / 2;

  for (; b < e; m = (b + e) / 2) {
    if (x > fLayers[m]->GetX()) {
      b = m + 1;
    }
    else {
      e = m;
    }
  }

  return m;

}             

//_____________________________________________________________________________
void AliTRDtracker::AliTRDpropagationLayer
                  ::SetZ(Double_t *center, Double_t *w, Double_t *wsensitive )
{
  //
  // set centers and the width of sectors
  //

  for (Int_t icham = 0; icham < AliTRDgeometry::kNcham; icham++) {
    fZc[icham]            = center[icham];  
    fZmax[icham]          = w[icham];
    fZmaxSensitive[icham] = wsensitive[icham];
  }  

}

//_____________________________________________________________________________
void AliTRDtracker::AliTRDpropagationLayer::SetHoles(Bool_t *holes)
{
  //
  // set centers and the width of sectors
  //

  fHole = kFALSE;

  for (Int_t icham = 0; icham < AliTRDgeometry::kNcham; icham++) {
    fIsHole[icham] = holes[icham]; 
    if (holes[icham]) {
      fHole = kTRUE;
    }
  }  

}

//_____________________________________________________________________________
void AliTRDtracker::AliTRDpropagationLayer
                  ::InsertCluster(AliTRDcluster *c, UInt_t index) 
{
  //
  // Insert cluster in cluster array.
  // Clusters are sorted according to Y coordinate.  
  //

  if (fTimeBinIndex < 0) { 
    //AliWarning("Attempt to insert cluster into non-sensitive time bin!\n");
    return;
  }

  if (fN == (Int_t) kMaxClusterPerTimeBin) {
    //AliWarning("Too many clusters !\n"); 
    return;
  }

  if (fN == 0) {
    fIndex[0]       = index; 
    fClusters[fN++] = c; 
    return;
  }

  Int_t i = Find(c->GetY());
  memmove(fClusters+i+1,fClusters+i,(fN-i)*sizeof(AliTRDcluster*));
  memmove(fIndex   +i+1,fIndex   +i,(fN-i)*sizeof(UInt_t)); 
  fIndex[i]    = index; 
  fClusters[i] = c; 
  fN++;

}  

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDpropagationLayer::Find(Float_t y) const 
{
  //
  // Returns index of the cluster nearest in Y    
  //

  if (fN <= 0) {
    return 0;
  }
  if (y <= fClusters[0]->GetY()) {
    return 0;
  }
  if (y >  fClusters[fN-1]->GetY()) {
    return fN;
  }

  Int_t b = 0;
  Int_t e = fN - 1;
  Int_t m = (b + e) / 2;

  for ( ; b < e; m = (b + e) / 2) {
    if (y > fClusters[m]->GetY()) {
      b = m + 1;
    }
    else {
      e = m;
    }
  }

  return m;

}    

//_____________________________________________________________________________
Int_t AliTRDtracker::AliTRDpropagationLayer
                   ::FindNearestCluster(Float_t y, Float_t z, Float_t maxroad
                                      , Float_t maxroadz) const 
{
  //
  // Returns index of the cluster nearest to the given y,z
  //

  Int_t   index   = -1;
  Int_t   maxn    = fN;
  Float_t mindist = maxroad;			

  for (Int_t i = Find(y-maxroad); i < maxn; i++) {
    AliTRDcluster *c = (AliTRDcluster *) (fClusters[i]);
    Float_t ycl = c->GetY();
    if (ycl > (y + maxroad)) {
      break;
    }
    if (TMath::Abs(c->GetZ() - z) > maxroadz) {
      continue;
    }
    if (TMath::Abs(ycl - y)       < mindist) {
      mindist = TMath::Abs(ycl - y);
      index   = fIndex[i];
    }
  }						

  return index;

}             

//_____________________________________________________________________________
Double_t AliTRDtracker::GetTiltFactor(const AliTRDcluster *c) 
{
  //
  // Returns correction factor for tilted pads geometry 
  //

  Int_t    det   = c->GetDetector();    
  Int_t    plane = fGeom->GetPlane(det);
  AliTRDpadPlane *padPlane = AliTRDCommonParam::Instance()->GetPadPlane(plane,0);
  Double_t h01   = TMath::Tan(-TMath::Pi() / 180.0 * padPlane->GetTiltingAngle());

  if (fNoTilt) {
    h01 = 0;
  }

  return h01;

}

//_____________________________________________________________________________
void AliTRDtracker::CookdEdxTimBin(AliTRDtrack &TRDtrack)
{
  //
  // This is setting fdEdxPlane and fTimBinPlane
  // Sums up the charge in each plane for track TRDtrack and also get the 
  // Time bin for Max. Cluster
  // Prashant Shukla (shukla@physi.uni-heidelberg.de)
  //

  Double_t  clscharge[AliESDtrack::kNPlane][AliESDtrack::kNSlice];
  Double_t  maxclscharge[AliESDtrack::kNPlane];
  Int_t     nCluster[AliESDtrack::kNPlane][AliESDtrack::kNSlice];
  Int_t     timebin[AliESDtrack::kNPlane];

  // Initialization of cluster charge per plane.  
  for (Int_t iPlane = 0; iPlane < AliESDtrack::kNPlane; iPlane++) {
    for (Int_t iSlice = 0; iSlice < AliESDtrack::kNSlice; iSlice++) {
      clscharge[iPlane][iSlice] = 0.0;
      nCluster[iPlane][iSlice]  = 0;
    }
  }

  // Initialization of cluster charge per plane.  
  for (Int_t iPlane = 0; iPlane < AliESDtrack::kNPlane; iPlane++) {
    timebin[iPlane]      =  -1;
    maxclscharge[iPlane] = 0.0;
  }

  // Loop through all clusters associated to track TRDtrack
  Int_t nClus = TRDtrack.GetNumberOfClusters();  // from Kalmantrack
  for (Int_t iClus = 0; iClus < nClus; iClus++) {
    Double_t charge = TRDtrack.GetClusterdQdl(iClus);
    Int_t    index  = TRDtrack.GetClusterIndex(iClus);
    AliTRDcluster *pTRDcluster = (AliTRDcluster *) GetCluster(index); 
    if (!pTRDcluster) {
      continue;
    }
    Int_t    tb     = pTRDcluster->GetLocalTimeBin();
    if (!tb) {
      continue;
    }
    Int_t detector  = pTRDcluster->GetDetector();
    Int_t iPlane    = fGeom->GetPlane(detector);
    Int_t iSlice    = tb * AliESDtrack::kNSlice / AliTRDtrack::kNtimeBins;
    clscharge[iPlane][iSlice] = clscharge[iPlane][iSlice] + charge;
    if (charge > maxclscharge[iPlane]) {
      maxclscharge[iPlane] = charge;
      timebin[iPlane]      = tb;
    }
    nCluster[iPlane][iSlice]++;
  } // End of loop over cluster

  // Setting the fdEdxPlane and fTimBinPlane variabales 
  Double_t totalCharge = 0.0;

  for (Int_t iPlane = 0; iPlane < AliESDtrack::kNPlane; iPlane++) {
    for (Int_t iSlice = 0; iSlice < AliESDtrack::kNSlice; iSlice++) {
      if (nCluster[iPlane][iSlice]) {
        clscharge[iPlane][iSlice] /= nCluster[iPlane][iSlice];
      }
      TRDtrack.SetPIDsignals(clscharge[iPlane][iSlice],iPlane,iSlice);
      totalCharge = totalCharge+clscharge[iPlane][iSlice];
    }
    TRDtrack.SetPIDTimBin(timebin[iPlane],iPlane);     
  }

  // Still needed ????
  //  Int_t i;
  //  Int_t nc=TRDtrack.GetNumberOfClusters(); 
  //  Float_t dedx=0;
  //  for (i=0; i<nc; i++) dedx += TRDtrack.GetClusterdQdl(i);
  //  dedx /= nc;
  //  for (Int_t iPlane = 0; iPlane < kNPlane; iPlane++) {
  //    TRDtrack.SetPIDsignals(dedx, iPlane);
  //    TRDtrack.SetPIDTimBin(timbin[iPlane], iPlane);
  //  }

}

//_____________________________________________________________________________
Int_t AliTRDtracker::FindClusters(Int_t sector, Int_t t0, Int_t t1
                                , AliTRDtrack *track
                                , Int_t *clusters, AliTRDtracklet &tracklet)
{
  //
  //
  // Try to find nearest clusters to the track in timebins from t0 to t1 
  //  
  // Correction coeficients - depend on TRD parameters - to be changed accordingly
  //

  Double_t x[100];
  Double_t yt[100];
  Double_t zt[100];
  Double_t xmean = 0.0;       // Reference x
  Double_t dz[10][100];
  Double_t dy[10][100];
  Float_t  zmean[100];
  Float_t  nmean[100];
  Int_t    clfound = 0;
  Int_t    indexes[10][100];  // Indexes of the clusters in the road
  Int_t    best[10][100];     // Index of best matching cluster 
  AliTRDcluster *cl[10][100]; // Pointers to the clusters in the road

  for (Int_t it = 0; it < 100; it++) {
    x[it]        = 0.0;
    yt[it]       = 0.0;
    zt[it]       = 0.0;
    clusters[it] = -2;
    zmean[it]    = 0.0;
    nmean[it]    = 0.0;
    for (Int_t ih = 0; ih < 10;ih++) {
      indexes[ih][it] = -2;   // Reset indexes1
      cl[ih][it]      =  0;
      dz[ih][it]      = -100.0;
      dy[ih][it]      = -100.0;
      best[ih][it]    =  0;
    }
  }  

  Double_t x0        = track->GetX();
  Double_t sigmaz    = TMath::Sqrt(TMath::Abs(track->GetSigmaZ2()));
  Int_t    nall      = 0;
  Int_t    nfound    = 0;
  Double_t h01       = 0.0;
  Int_t    plane     = -1;
  Int_t    detector  = -1;
  Float_t  padlength = 0.0;
  AliTRDtrack track2(* track);
  Float_t  snpy      = track->GetSnp();
  Float_t  tany      = TMath::Sqrt(snpy*snpy / (1.0 - snpy*snpy)); 
  if (snpy < 0.0) {
    tany *= -1.0;
  }

  Double_t sy2       = ExpectedSigmaY2(x0,track->GetTgl(),track->GetPt());
  Double_t sz2       = ExpectedSigmaZ2(x0,track->GetTgl());
  Double_t road      = 15.0 * TMath::Sqrt(track->GetSigmaY2() + sy2);
  if (road > 6.0) {
    road = 6.0;
  }
  //road = 20.0;

  for (Int_t it = 0; it < t1-t0; it++) {

    Double_t maxChi2[2] = { fgkMaxChi2, fgkMaxChi2 };      
    AliTRDpropagationLayer &timeBin = *(fTrSec[sector]->GetLayer(it+t0));
    if (timeBin == 0) {
      continue; // No indexes1
    }

    Int_t    maxn = timeBin;
    x[it]  = timeBin.GetX();
    track2.PropagateTo(x[it]);
    yt[it] = track2.GetY();
    zt[it] = track2.GetZ();
    
    Double_t y    = yt[it];
    Double_t z    = zt[it];
    Double_t chi2 = 1000000.0;
    nall++;

    //
    // Find 2 nearest cluster at given time bin
    // 
    int checkPoint[4] = {0,0,0,0};
    double minY = 123456789;
    double minD[2] = {1,1};

    for (Int_t i = timeBin.Find(y - road); i < maxn; i++) {
      //for (Int_t i = 0; i < maxn; i++) {

      AliTRDcluster *c = (AliTRDcluster *) (timeBin[i]);
      h01 = GetTiltFactor(c);
      if (plane < 0) {
	Int_t det = c->GetDetector();
	plane     = fGeom->GetPlane(det);
	padlength = TMath::Sqrt(c->GetSigmaZ2() * 12.0);
      }

      //if (c->GetLocalTimeBin()==0) continue;
      if (c->GetY() > (y + road)) {
        break;
      }

      fHDeltaX->Fill(c->GetX() - x[it]);
      //printf("%f\t%f\t%f \n", c->GetX(),  x[it], c->GetX()-x[it]);

      if (TMath::Abs(c->GetY()-y) < TMath::Abs(minY)) {
	minY = c->GetY()-y;
	minD[0] = c->GetY()-y;
	minD[1] = c->GetZ()-z;
      }

      checkPoint[0]++;

      fHMinZ->Fill(c->GetZ() - z);
      if ((c->GetZ() - z) * (c->GetZ() - z) > 2 * (12.0 * sz2)) {
        continue;
      }
      checkPoint[1]++;

      Double_t dist = TMath::Abs(c->GetZ() - z);
      if (dist > (0.5 * padlength + 6.0 * sigmaz)) { // 0.5
        continue;   // 6 sigma boundary cut
      }
      checkPoint[2]++;

      Double_t cost = 0.0;
      // Sigma boundary cost function
      if (dist> (0.5 * padlength - sigmaz)){ 
	cost =  (dist - 0.5*padlength) / (2.0 * sigmaz);
	if (cost > -1) {
          cost = (cost + 1.0) * (cost + 1.0);
	}
	else {
          cost = 0.0;
	}
      }
      //Int_t label = TMath::Abs(track->GetLabel());
      //if (c->GetLabel(0)!=label && c->GetLabel(1)!=label&&c->GetLabel(2)!=label) continue;
      chi2 = track2.GetPredictedChi2(c,h01) + cost;
      clfound++;      

      if (chi2 > maxChi2[1]) {
        continue;
      }
      checkPoint[3]++;

      detector = c->GetDetector();
      // Store the clusters in the road
      for (Int_t ih = 2; ih < 9; ih++) {  
	if (cl[ih][it] == 0) {
	  cl[ih][it]      = c;
	  indexes[ih][it] = timeBin.GetIndex(i); // Index - 9 - reserved for outliers
	  break;
	}
      }

      if (chi2 < maxChi2[0]) {
	maxChi2[1]     = maxChi2[0];
	maxChi2[0]     = chi2;
	indexes[1][it] = indexes[0][it];
	cl[1][it]      = cl[0][it];
	indexes[0][it] = timeBin.GetIndex(i);
	cl[0][it]      = c;
	continue;
      }
      maxChi2[1]     = chi2;
      cl[1][it]      = c;
      indexes[1][it] = timeBin.GetIndex(i); 

    }
    
    for(int iCheckPoint = 0; iCheckPoint<4; iCheckPoint++)
      fHFindCl[iCheckPoint]->Fill(checkPoint[iCheckPoint]);

    if (checkPoint[3]) {
      if (track->GetPt() > 0) fHMinYPos->Fill(minY);
      else fHMinYNeg->Fill(minY);

    fHMinD->Fill(minD[0], minD[1]);
     }

    if (cl[0][it]) {
      nfound++;
      xmean += x[it];
    }

  }

  if (nfound < 4) {
    return 0;  
  }
  xmean /= Float_t(nfound);  // Middle x
  track2.PropagateTo(xmean); // Propagate track to the center

  //
  // Choose one of the variants
  //
  Int_t    changes[10];
  Float_t  sumz   = 0.0;
  Float_t  sum    = 0.0;
  Double_t sumdy  = 0.0;
  Double_t sumdy2 = 0.0;
  Double_t sumx   = 0.0;
  Double_t sumxy  = 0.0;
  Double_t sumx2  = 0.0;
  Double_t mpads  = 0.0;

  Int_t    ngood[10];
  Int_t    nbad[10];

  Double_t meanz[10];
  Double_t moffset[10];    // Mean offset
  Double_t mean[10];       // Mean value
  Double_t angle[10];      // Angle

  Double_t smoffset[10];   // Sigma of mean offset
  Double_t smean[10];      // Sigma of mean value
  Double_t sangle[10];     // Sigma of angle
  Double_t smeanangle[10]; // Correlation

  Double_t sigmas[10];     
  Double_t tchi2s[10];     // Chi2s for tracklet

  for (Int_t it = 0; it < 10; it++) {

    ngood[it]      = 0;
    nbad[it]       = 0;

    meanz[it]      = 0.0;
    moffset[it]    = 0.0;    // Mean offset
    mean[it]       = 0.0;    // Mean value
    angle[it]      = 0.0;    // Angle

    smoffset[it]   = 1.0e5; // Sigma of mean offset
    smean[it]      = 1.0e5; // Sigma of mean value
    sangle[it]     = 1.0e5; // Sigma of angle
    smeanangle[it] = 0.0;    // Correlation

    sigmas[it]     = 1.0e5;     
    tchi2s[it]     = 1.0e5; // Chi2s for tracklet

  }

  //
  // Calculate zmean
  //
  for (Int_t it = 0; it < t1 - t0; it++) {
    if (!cl[0][it]) {
      continue;
    }
    for (Int_t dt = -3; dt <= 3; dt++) {
      if (it+dt <     0) {
        continue;
      }
      if (it+dt > t1-t0) {
        continue;
      }
      if (!cl[0][it+dt]) {
        continue;
      }
      zmean[it] += cl[0][it+dt]->GetZ();
      nmean[it] += 1.0;
    }
    zmean[it] /= nmean[it]; 
  }

  for (Int_t it = 0; it < t1 - t0; it++) {

    best[0][it] = 0;

    for (Int_t ih = 0; ih < 10; ih++) {
      dz[ih][it] = -100.0;
      dy[ih][it] = -100.0;
      if (!cl[ih][it]) {
        continue;
      }
      Double_t xcluster = cl[ih][it]->GetX();
      Double_t ytrack;
      Double_t ztrack;
      track2.GetProlongation(xcluster,ytrack,ztrack );
      dz[ih][it] = cl[ih][it]->GetZ()- ztrack;                   // Calculate distance from track in z
      dy[ih][it] = cl[ih][it]->GetY() + dz[ih][it]*h01 - ytrack; // and in y
    }

    // Minimize changes
    if (!cl[0][it]) {
      continue;
    }
    if ((TMath::Abs(cl[0][it]->GetZ()-zmean[it]) > padlength * 0.8) &&
        (cl[1][it])) {
      if (TMath::Abs(cl[1][it]->GetZ()-zmean[it]) < padlength * 0.5) {
	best[0][it] = 1;
      }
    }

  }

  //
  // Iterative choice of "best path"
  //
  Int_t label    = TMath::Abs(track->GetLabel());
  Int_t bestiter = 0;

  for (Int_t iter = 0; iter < 9; iter++) {

    changes[iter] = 0;
    sumz          = 0; 
    sum           = 0; 
    sumdy         = 0;
    sumdy2        = 0;
    sumx          = 0;
    sumx2         = 0;
    sumxy         = 0;
    mpads         = 0; 
    ngood[iter]   = 0; 
    nbad[iter]    = 0;
 
    // Linear fit
    for (Int_t it = 0; it < t1 - t0; it++) {

      if (!cl[best[iter][it]][it]) {
        continue;
      }

      // Calculates pad-row changes
      Double_t zbefore = cl[best[iter][it]][it]->GetZ();
      Double_t zafter  = cl[best[iter][it]][it]->GetZ();
      for (Int_t itd = it - 1; itd >= 0; itd--) {
	if (cl[best[iter][itd]][itd]) {
	  zbefore = cl[best[iter][itd]][itd]->GetZ();
	  break;
	}
      }
      for (Int_t itd = it + 1; itd < t1 - t0; itd++) {
	if (cl[best[iter][itd]][itd]) {
	  zafter  = cl[best[iter][itd]][itd]->GetZ();
	  break;
	}
      }
      if ((TMath::Abs(cl[best[iter][it]][it]->GetZ()-zbefore) > 0.1) &&
          (TMath::Abs(cl[best[iter][it]][it]->GetZ()- zafter) > 0.1)) {
        changes[iter]++;
      }

      Double_t dx = x[it]-xmean; // Distance to reference x
      sumz   += cl[best[iter][it]][it]->GetZ();      
      sum++;
      sumdy  += dy[best[iter][it]][it];
      sumdy2 += dy[best[iter][it]][it]*dy[best[iter][it]][it];
      sumx   += dx;
      sumx2  += dx*dx;
      sumxy  += dx*dy[best[iter][it]][it];
      mpads  += cl[best[iter][it]][it]->GetNPads();
      if ((cl[best[iter][it]][it]->GetLabel(0) == label) || 
          (cl[best[iter][it]][it]->GetLabel(1) == label) || 
          (cl[best[iter][it]][it]->GetLabel(2) == label)) {
	ngood[iter]++;
      }
      else {
	nbad[iter]++;
      }

    }

    //
    // calculates line parameters
    //
    Double_t det  = sum*sumx2 - sumx*sumx;
    angle[iter]   = (sum*sumxy - sumx*sumdy) / det;
    mean[iter]    = (sumx2*sumdy - sumx*sumxy) / det;
    meanz[iter]   = sumz / sum;    
    moffset[iter] = sumdy / sum;
    mpads        /= sum;   // Mean number of pads

    Double_t sigma2 = 0.0; // Normalized residuals - for line fit
    Double_t sigma1 = 0.0; // Normalized residuals - constant fit

    for (Int_t it = 0; it < t1 - t0; it++) {
      if (!cl[best[iter][it]][it]) {
        continue;
      }
      Double_t dx  = x[it] - xmean;
      Double_t ytr = mean[iter] + angle[iter] * dx;
      sigma2 += (dy[best[iter][it]][it] - ytr) 
              * (dy[best[iter][it]][it] - ytr);
      sigma1 += (dy[best[iter][it]][it] - moffset[iter]) 
              * (dy[best[iter][it]][it] - moffset[iter]);
      sum++;
    }
    sigma2          /= (sum - 2);                  // Normalized residuals
    sigma1          /= (sum - 1);                  // Normalized residuals
    smean[iter]      = sigma2 * (sumx2 / det);     // Estimated error2 of mean
    sangle[iter]     = sigma2 * ( sum  / det);     // Estimated error2 of angle
    smeanangle[iter] = sigma2 * (-sumx / det);     // Correlation
    sigmas[iter]     = TMath::Sqrt(sigma1);
    smoffset[iter]   = (sigma1 / sum) + 0.01*0.01; // Sigma of mean offset + unisochronity sigma 

    //
    // Iterative choice of "better path"
    //
    for (Int_t it = 0; it < t1 - t0; it++) {

      if (!cl[best[iter][it]][it]) {
        continue;
      }

      // Add unisochronity + angular effect contribution
      Double_t sigmatr2 = smoffset[iter] + 0.5*tany*tany;                   
      Double_t sweight  = 1.0/sigmatr2 + 1.0/track->GetSigmaY2();
      Double_t weighty  = (moffset[iter] / sigmatr2) / sweight;             // Weighted mean
      Double_t sigmacl  = TMath::Sqrt(sigma1*sigma1 + track->GetSigmaY2());
      Double_t mindist  = 100000.0; 
      Int_t    ihbest   = 0;

      for (Int_t ih = 0; ih < 10; ih++) {
	if (!cl[ih][it]) {
          break;
	}
	Double_t dist2 = (dy[ih][it] - weighty) / sigmacl;
	dist2 *= dist2; // Chi2 distance
	if (dist2 < mindist) {
	  mindist = dist2;
	  ihbest  = ih;
	}
      }
      best[iter+1][it] = ihbest;
    }

    //
    // Update best hypothesy if better chi2 according tracklet position and angle
    //
    Double_t sy2 = smean[iter]  + track->GetSigmaY2();
    Double_t sa2 = sangle[iter] + track->GetSigmaSnp2(); // track->fCee;
    Double_t say = track->GetSigmaSnpY();                // track->fCey;
    //Double_t chi20 = mean[bestiter]*mean[bestiter ] / sy2+angle[bestiter]*angle[bestiter]/sa2;
    //Double_t chi21 = mean[iter]*mean[iter]          / sy2+angle[iter]*angle[iter]/sa2;

    Double_t detchi    = sy2*sa2 - say*say;
    Double_t invers[3] = {sa2/detchi,sy2/detchi,-say/detchi}; // Inverse value of covariance matrix  
    
    Double_t chi20 = mean[bestiter]  * mean[bestiter]       * invers[0]
                   + angle[bestiter] * angle[bestiter]      * invers[1]
                   + 2.0 * mean[bestiter] * angle[bestiter] * invers[2];
    Double_t chi21 = mean[iter]  * mean[iter]               * invers[0]
                   + angle[iter] * angle[iter]              * invers[1]
                   + 2.0 * mean[iter] * angle[iter]         * invers[2];
    tchi2s[iter]   = chi21;

    if ((changes[iter] <= changes[bestiter]) && 
        (chi21 < chi20)) {
      bestiter = iter;
    }

  }

  //
  // Set clusters 
  //
  Double_t sigma2     = sigmas[0]; // Choose as sigma from 0 iteration
  Short_t  maxpos     = -1;
  Float_t  maxcharge  =  0.0;
  Short_t  maxpos4    = -1;
  Float_t  maxcharge4 =  0.0;
  Short_t  maxpos5    = -1;
  Float_t  maxcharge5 =  0.0;

  //if (tchi2s[bestiter]>25.) sigma2*=tchi2s[bestiter]/25.;
  //if (tchi2s[bestiter]>25.) sigma2=1000.;  // dont'accept

  Double_t exB         = AliTRDcalibDB::Instance()->GetOmegaTau(AliTRDcalibDB::Instance()->GetVdrift(0,0,0)
                                                               ,-AliTracker::GetBz()*0.1);
  Double_t expectederr = sigma2*sigma2 + 0.01*0.01;
  if (mpads > 3.5) {
    expectederr += (mpads - 3.5) * 0.04;
  }
  if (changes[bestiter] > 1) {
    expectederr += changes[bestiter] * 0.01; 
  }
  expectederr += (0.03 * (tany-exB)*(tany-exB)) * 15.0;
  //if (tchi2s[bestiter]>18.) expectederr*= tchi2s[bestiter]/18.;
  //expectederr+=10000;

  for (Int_t it = 0; it < t1 - t0; it++) {

    if (!cl[best[bestiter][it]][it]) {
      continue;
    }

    cl[best[bestiter][it]][it]->SetSigmaY2(expectederr); // Set cluster error
    if (!cl[best[bestiter][it]][it]->IsUsed()) {
      cl[best[bestiter][it]][it]->SetY(cl[best[bestiter][it]][it]->GetY()); 
      //cl[best[bestiter][it]][it]->Use();
    }

    // Time bins with maximal charge
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ()) > maxcharge) {
      maxcharge = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
      maxpos    = cl[best[bestiter][it]][it]->GetLocalTimeBin();
    }
    
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ()) > maxcharge4) {
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin() >= 4) {
	maxcharge4 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos4    = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }

    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ()) > maxcharge5) {
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin() >= 5) {
	maxcharge5 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos5    = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }

    // Time bins with maximal charge
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ()) > maxcharge) {
      maxcharge = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
      maxpos    = cl[best[bestiter][it]][it]->GetLocalTimeBin();
    }
    
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ()) > maxcharge4) {
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin() >= 4) {
	maxcharge4 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos4    = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }

    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ()) > maxcharge5) {
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin() >= 5) {
	maxcharge5 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos5    = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }

    clusters[it+t0] = indexes[best[bestiter][it]][it];    

    // Still needed ????
    //if (cl[best[bestiter][it]][it]->GetLocalTimeBin()>4 && 
    //cl[best[bestiter][it]][it]->GetLocalTimeBin()<18) clusters[it+t0]
    // = indexes[best[bestiter][it]][it];    //Test

  } 

  //
  // Set tracklet parameters
  //
  Double_t trackleterr2 = smoffset[bestiter] + 0.01*0.01;
  if (mpads > 3.5) {
    trackleterr2 += (mpads - 3.5) * 0.04;
  }
  trackleterr2 += changes[bestiter] * 0.01;
  trackleterr2 *= TMath::Max(14.0 - nfound,1.0);
  trackleterr2 += 0.2 * (tany-exB)*(tany-exB); 

  // Set tracklet parameters
  tracklet.Set(xmean
              ,track2.GetY() + moffset[bestiter]
              ,meanz[bestiter]
              ,track2.GetAlpha()
              ,trackleterr2);
  tracklet.SetTilt(h01);
  tracklet.SetP0(mean[bestiter]);
  tracklet.SetP1(angle[bestiter]);
  tracklet.SetN(nfound);
  tracklet.SetNCross(changes[bestiter]);
  tracklet.SetPlane(plane);
  tracklet.SetSigma2(expectederr);
  tracklet.SetChi2(tchi2s[bestiter]);
  tracklet.SetMaxPos(maxpos,maxpos4,maxpos5);
  track->SetTracklets(plane,tracklet);
  track->SetNWrong(track->GetNWrong() + nbad[0]);

  //
  // Debuging part
  //
  TClonesArray array0("AliTRDcluster");
  TClonesArray array1("AliTRDcluster");
  array0.ExpandCreateFast(t1 - t0 + 1);
  array1.ExpandCreateFast(t1 - t0 + 1);
  TTreeSRedirector &cstream = *fDebugStreamer;
  AliTRDcluster dummy;
  Double_t dy0[100];
  Double_t dyb[100]; 

  for (Int_t it = 0; it < t1 - t0; it++) {
    dy0[it] = dy[0][it];
    dyb[it] = dy[best[bestiter][it]][it];
    if (cl[0][it]) {
      new(array0[it]) AliTRDcluster(*cl[0][it]);
    }
    else {
      new(array0[it]) AliTRDcluster(dummy);
    }
    if(cl[best[bestiter][it]][it]) {
      new(array1[it]) AliTRDcluster(*cl[best[bestiter][it]][it]);
    }
    else{
      new(array1[it]) AliTRDcluster(dummy);
    }
  }
  TGraph graph0(t1-t0,x,dy0);
  TGraph graph1(t1-t0,x,dyb);
  TGraph graphy(t1-t0,x,yt);
  TGraph graphz(t1-t0,x,zt);

  if (AliTRDReconstructor::StreamLevel() > 0) {
    cstream << "tracklet"
	    << "track.="      << track              // Track parameters
	    << "tany="        << tany               // Tangent of the local track angle 
	    << "xmean="       << xmean              // Xmean - reference x of tracklet  
	    << "tilt="        << h01                // Tilt angle
	    << "nall="        << nall               // Number of foundable clusters 
	    << "nfound="      << nfound             // Number of found clusters
	    << "clfound="     << clfound            // Total number of found clusters in road 
	    << "mpads="       << mpads              // Mean number of pads per cluster
	    << "plane="       << plane              // Plane number 
	    << "detector="    << detector           // Detector number
	    << "road="        << road               // The width of the used road
	    << "graph0.="     << &graph0            // x - y = dy for closest cluster
	    << "graph1.="     << &graph1            // x - y = dy for second closest cluster    
	    << "graphy.="     << &graphy            // y position of the track
	    << "graphz.="     << &graphz            // z position of the track
            //<< "fCl.="        << &array0            // closest cluster
            //<< "fCl2.="       << &array1            // second closest cluster
	    << "maxpos="      << maxpos             // Maximal charge postion
	    << "maxcharge="   << maxcharge          // Maximal charge 
	    << "maxpos4="     << maxpos4            // Maximal charge postion - after bin 4
	    << "maxcharge4="  << maxcharge4         // Maximal charge         - after bin 4
	    << "maxpos5="     << maxpos5            // Maximal charge postion - after bin 5
	    << "maxcharge5="  << maxcharge5         // Maximal charge         - after bin 5
	    << "bestiter="    << bestiter           // Best iteration number 
	    << "tracklet.="   << &tracklet          // Corrspond to the best iteration
	    << "tchi20="      << tchi2s[0]          // Chi2 of cluster in the 0 iteration
	    << "tchi2b="      << tchi2s[bestiter]   // Chi2 of cluster in the best  iteration
	    << "sigmas0="     << sigmas[0]          // Residuals sigma 
	    << "sigmasb="     << sigmas[bestiter]   // Residulas sigma
	    << "ngood0="      << ngood[0]           // Number of good clusters in 0 iteration
	    << "nbad0="       << nbad[0]            // Number of bad clusters in 0 iteration
	    << "ngoodb="      << ngood[bestiter]    //                        in best iteration    
	    << "nbadb="       << nbad[bestiter]     //                        in best iteration
	    << "changes0="    << changes[0]         // Changes of pardrows in iteration number 0 
	    << "changesb="    << changes[bestiter]  // Changes of pardrows in best iteration
	    << "moffset0="    << moffset[0]         // Offset fixing angle in iter=0
	    << "smoffset0="   << smoffset[0]        // Sigma of offset fixing angle in iter=0
	    << "moffsetb="    << moffset[bestiter]  // Offset fixing angle in iter=best
	    << "smoffsetb="   << smoffset[bestiter] // Sigma of offset fixing angle in iter=best
	    << "mean0="       << mean[0]            // Mean dy in iter=0;
	    << "smean0="      << smean[0]           // Sigma of mean dy in iter=0
	    << "meanb="       << mean[bestiter]     // Mean dy in iter=best
	    << "smeanb="      << smean[bestiter]    // Sigma of mean dy in iter=best
	    << "angle0="      << angle[0]           // Angle deviation in the iteration number 0 
	    << "sangle0="     << sangle[0]          // Sigma of angular deviation in iteration number 0
	    << "angleb="      << angle[bestiter]    // Angle deviation in the best iteration   
	    << "sangleb="     << sangle[bestiter]   // Sigma of angle deviation in the best iteration   
	    << "expectederr=" << expectederr        // Expected error of cluster position
	    << "\n";
  }

  return nfound;

}

//_____________________________________________________________________________
Int_t AliTRDtracker::Freq(Int_t n, const Int_t *inlist
                        , Int_t *outlist, Bool_t down)
{    
  //
  // Sort eleements according occurancy 
  // The size of output array has is 2*n 
  //

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

//_____________________________________________________________________________
AliTRDtrack *AliTRDtracker::RegisterSeed(AliTRDseed *seeds, Double_t *params)
{
  //
  // Register a seed
  //

  Double_t alpha = AliTRDgeometry::GetAlpha();
  Double_t shift = AliTRDgeometry::GetAlpha()/2.0;
  Double_t c[15];

  c[ 0] = 0.2;
  c[ 1] = 0.0; c[ 2] = 2.0;
  c[ 3] = 0.0; c[ 4] = 0.0; c[ 5] = 0.02;
  c[ 6] = 0.0; c[ 7] = 0.0; c[ 8] = 0.0;  c[ 9] = 0.1;
  c[10] = 0.0; c[11] = 0.0; c[12] = 0.0;  c[13] = 0.0; c[14] = params[5]*params[5]*0.01;

  Int_t index = 0;
  AliTRDcluster *cl = 0;

  for (Int_t ilayer = 0; ilayer < 6; ilayer++) {
    if (seeds[ilayer].IsOK()) {
      for (Int_t itime = 22; itime > 0; itime--) {
	if (seeds[ilayer].GetIndexes(itime) > 0) {
	  index = seeds[ilayer].GetIndexes(itime);
	  cl    = seeds[ilayer].GetClusters(itime);
	  break;
	}
      }
    }
    if (index > 0) {
      break;
    }
  }
  if (cl == 0) {
    return 0;
  }

  AliTRDtrack *track = new AliTRDtrack(cl
                                      ,index
                                      ,&params[1]
                                      ,c
                                      ,params[0]
                                      ,params[6]*alpha+shift);
  track->PropagateTo(params[0]-5.0);
  track->ResetCovariance(1);

  Int_t rc = FollowBackProlongation(*track);
  if (rc < 30) {
    delete track;
    track = 0;
  }
  else {
    track->CookdEdx();
    CookdEdxTimBin(*track);
    CookLabel(track,0.9);
  }

  return track;
}

//////////////////////////////////////////////////////////////////////////////////////////

void AliTRDtracker::InitLogHists() {
 
  fHBackfit = new TH1D("logTRD_backfit", "", 40, -0.5, 39.5);  
  fHRefit = new TH1D("logTRD_refit", "", 40, -0.5, 39.5);
  fHClSearch = new TH1D("logTRD_clSearch", "", 60, -0.5, 59.5); 
  
  fHX = new TH1D("logTRD_X", ";x (cm)", 200, 50, 400);
  fHNCl = new TH1D("logTRD_ncl", "", 40, -0.5, 39.5);
  fHNClTrack = new TH1D("logTRD_nclTrack", "", 180, -0.5, 179.5);
  
  fHMinYPos = new TH1D("logTRD_minYPos", ";#delta Y (cm)", 400, -6, 6);
  fHMinYNeg = new TH1D("logTRD_minYNeg", ";#delta Y (cm)", 400, -6, 6);
  fHMinZ = new TH1D("logTRD_minZ", ";#delta Z (cm)", 400, -20, 20);
  fHMinD = new TH2D("logTRD_minD", ";#delta Y (cm);#delta Z (cm)", 100, -6, 6, 100, -50, 50);
  
  fHDeltaX = new TH1D("logTRD_deltaX", ";#delta X (cm)", 100, -5, 5);
  fHXCl = new TH1D("logTRD_xCl", ";cluster x position (cm)", 1000, 280, 380);
  
  const char *nameFindCl[4] = {"logTRD_clY", "logTRD_clZ", "logTRD_clB", "logTRD_clG"};
  
  for(int i=0; i<4; i++) {
    fHFindCl[i] = new TH1D(nameFindCl[i], "", 30, -0.5, 29.5);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

void AliTRDtracker::SaveLogHists() {

  TDirectory *sav = gDirectory;
  TFile *logFile = 0;

  TSeqCollection *col = gROOT->GetListOfFiles();
  int N = col->GetEntries();
  for(int i=0; i<N; i++) {
    logFile = (TFile*)col->At(i);
    if (strstr(logFile->GetName(), "AliESDs.root")) break;
  }
   
  logFile->cd();
  fHBackfit->Write(fHBackfit->GetName(), TObject::kOverwrite);
  fHRefit->Write(fHRefit->GetName(), TObject::kOverwrite);
  fHClSearch->Write(fHClSearch->GetName(), TObject::kOverwrite);
  fHX->Write(fHX->GetName(), TObject::kOverwrite);
  fHNCl->Write(fHNCl->GetName(), TObject::kOverwrite);
  fHNClTrack->Write(fHNClTrack->GetName(), TObject::kOverwrite);

  fHMinYPos->Write(fHMinYPos->GetName(), TObject::kOverwrite);
  fHMinYNeg->Write(fHMinYNeg->GetName(), TObject::kOverwrite);
  fHMinD->Write(fHMinD->GetName(), TObject::kOverwrite);
  fHMinZ->Write(fHMinZ->GetName(), TObject::kOverwrite);
  
  fHDeltaX->Write(fHDeltaX->GetName(), TObject::kOverwrite);
  fHXCl->Write(fHXCl->GetName(), TObject::kOverwrite);


  for(int i=0; i<4; i++)
    fHFindCl[i]->Write(fHFindCl[i]->GetName(), TObject::kOverwrite);

  logFile->Flush();
  
  sav->cd();
}

//////////////////////////////////////////////////////////////////////////////////////////
