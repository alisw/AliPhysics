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
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>  
#include <TObjArray.h> 

#include "AliTRDgeometry.h"
#include "AliTRDparameter.h"
#include "AliTRDgeometryHole.h"
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliBarrelTrack.h"
#include "AliESD.h"

#include "AliTRDtracker.h"

ClassImp(AliTRDtracker) 

  const  Float_t     AliTRDtracker::fgkSeedDepth          = 0.5; 
  const  Float_t     AliTRDtracker::fgkSeedStep           = 0.10;   
  const  Float_t     AliTRDtracker::fgkSeedGap            = 0.25;  

  const  Float_t     AliTRDtracker::fgkMaxSeedDeltaZ12    = 40.;  
  const  Float_t     AliTRDtracker::fgkMaxSeedDeltaZ      = 25.;  
  const  Float_t     AliTRDtracker::fgkMaxSeedC           = 0.0052; 
  const  Float_t     AliTRDtracker::fgkMaxSeedTan         = 1.2;  
  const  Float_t     AliTRDtracker::fgkMaxSeedVertexZ     = 150.; 

  const  Double_t    AliTRDtracker::fgkSeedErrorSY        = 0.2;
  const  Double_t    AliTRDtracker::fgkSeedErrorSY3       = 2.5;
  const  Double_t    AliTRDtracker::fgkSeedErrorSZ        = 0.1;

  const  Float_t     AliTRDtracker::fgkMinClustersInSeed  = 0.7;  

  const  Float_t     AliTRDtracker::fgkMinClustersInTrack = 0.5;  
  const  Float_t     AliTRDtracker::fgkMinFractionOfFoundClusters = 0.8;  

  const  Float_t     AliTRDtracker::fgkSkipDepth          = 0.3;
  const  Float_t     AliTRDtracker::fgkLabelFraction      = 0.8;  
  const  Float_t     AliTRDtracker::fgkWideRoad           = 20.;

  const  Double_t    AliTRDtracker::fgkMaxChi2            = 12.; 

const Int_t AliTRDtracker::fgkFirstPlane = 5;
const Int_t AliTRDtracker::fgkLastPlane = 17;


//____________________________________________________________________
AliTRDtracker::AliTRDtracker():AliTracker(),
			       fGeom(0),
			       fPar(0),
			       fNclusters(0),
			       fClusters(0),
			       fNseeds(0),
			       fSeeds(0),
			       fNtracks(0),
			       fTracks(0),
			       fSY2corr(0),
			       fSZ2corr(0),
			       fTimeBinsPerPlane(0),
			       fMaxGap(0),
			       fVocal(kFALSE),
			       fAddTRDseeds(kFALSE),
			       fNoTilt(kFALSE)
{
  // Default constructor

  for(Int_t i=0;i<kTrackingSectors;i++) fTrSec[i]=0;
  for(Int_t j=0;j<5;j++)
    for(Int_t k=0;k<18;k++) fHoles[j][k]=kFALSE;
} 
//____________________________________________________________________
AliTRDtracker::AliTRDtracker(const TFile *geomfile):AliTracker()
{
  // 
  //  Main constructor
  //  

  //Float_t fTzero = 0;
   
  fAddTRDseeds = kFALSE;
  fGeom = NULL;
  fNoTilt = kFALSE;
  
  TDirectory *savedir=gDirectory; 
  TFile *in=(TFile*)geomfile;  
  if (!in->IsOpen()) {
    printf("AliTRDtracker::AliTRDtracker(): geometry file is not open!\n");
    printf("    DETAIL TRD geometry and DEFAULT TRD parameter will be used\n");
  }
  else {
    in->cd();  
//    in->ls();
    fGeom = (AliTRDgeometry*) in->Get("TRDgeometry");
    fPar  = (AliTRDparameter*) in->Get("TRDparameter");
//    fGeom->Dump();
  }

  if(fGeom) {
    //    fTzero = geo->GetT0();
    printf("Found geometry version %d on file \n", fGeom->IsVersion());
  }
  else { 
    printf("AliTRDtracker::AliTRDtracker(): can't find TRD geometry!\n");
    //printf("The DETAIL TRD geometry will be used\n");
    //fGeom = new AliTRDgeometryDetail();
    fGeom = new AliTRDgeometryHole();
    fGeom->SetPHOShole();
    fGeom->SetRICHhole();    
  } 

  if (!fPar) {  
    printf("AliTRDtracker::AliTRDtracker(): can't find TRD parameter!\n");
    printf("The DEFAULT TRD parameter will be used\n");
    fPar = new AliTRDparameter();
  }
  fPar->ReInit();

  savedir->cd();  


  //  fGeom->SetT0(fTzero);

  fNclusters = 0;
  fClusters  = new TObjArray(2000); 
  fNseeds    = 0;
  fSeeds     = new TObjArray(2000);
  fNtracks   = 0;
  fTracks    = new TObjArray(1000);

  for(Int_t geomS = 0; geomS < kTrackingSectors; geomS++) {
    Int_t trS = CookSectorIndex(geomS);
    fTrSec[trS] = new AliTRDtrackingSector(fGeom, geomS, fPar);
    for (Int_t icham=0;icham<AliTRDgeometry::kNcham; icham++){
      fHoles[icham][trS]=fGeom->IsHole(0,icham,geomS);
    }
  }

  Float_t tiltAngle = TMath::Abs(fPar->GetTiltingAngle()); 
  if(tiltAngle < 0.1) {
    fNoTilt = kTRUE;
  }

  fSY2corr = 0.2;
  fSZ2corr = 120.;      

  if(fNoTilt && (tiltAngle > 0.1)) fSY2corr = fSY2corr + tiltAngle * 0.05; 


  // calculate max gap on track

  Double_t dxAmp = (Double_t) fGeom->CamHght();   // Amplification region
  Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region

  Double_t dx = (Double_t) fPar->GetTimeBinSize();   
  Int_t tbAmp = fPar->GetTimeBefore();
  Int_t maxAmp = (Int_t) ((dxAmp+0.000001)/dx);
  if(kTRUE) maxAmp = 0;  // intentional until we change the parameter class 
  Int_t tbDrift = fPar->GetTimeMax();
  Int_t maxDrift = (Int_t) ((dxDrift+0.000001)/dx);

  tbDrift = TMath::Min(tbDrift,maxDrift);
  tbAmp = TMath::Min(tbAmp,maxAmp);

  fTimeBinsPerPlane = tbAmp + tbDrift;
  fMaxGap = (Int_t) (fTimeBinsPerPlane * fGeom->Nplan() * fgkSkipDepth);

  fVocal = kFALSE;


  // Barrel Tracks [SR, 03.04.2003]

  fBarrelFile = 0;
  fBarrelTree = 0;
  fBarrelArray = 0;
  fBarrelTrack = 0;

  savedir->cd();
}   

//___________________________________________________________________
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
  delete fPar;  

  for(Int_t geomS = 0; geomS < kTrackingSectors; geomS++) {
    delete fTrSec[geomS];
  }
}   

//_____________________________________________________________________

void AliTRDtracker::SetBarrelTree(const char *mode) {
  //
  //
  //

  if (!IsStoringBarrel()) return;

  TDirectory *sav = gDirectory;
  if (!fBarrelFile) fBarrelFile = new TFile("AliBarrelTracks.root", "UPDATE");

  char buff[40];
  sprintf(buff,  "BarrelTRD_%d_%s", GetEventNumber(), mode);

  fBarrelFile->cd();
  fBarrelTree = new TTree(buff, "Barrel TPC tracks");
  
  Int_t nRefs = fgkLastPlane - fgkFirstPlane + 1;

  if (!fBarrelArray) fBarrelArray = new TClonesArray("AliBarrelTrack", nRefs);
  for(Int_t i=0; i<nRefs; i++) new((*fBarrelArray)[i]) AliBarrelTrack();
  
  fBarrelTree->Branch("tracks", &fBarrelArray);
  sav->cd();
}
  
//_____________________________________________________________________

void AliTRDtracker::StoreBarrelTrack(AliTRDtrack *ps, Int_t refPlane, Int_t isIn) {
  //
  //
  //
  
  if (!IsStoringBarrel()) return;
  
  static Int_t nClusters;
  static Int_t nWrong;
  static Double_t chi2;
  static Int_t index;
  static Bool_t wasLast = kTRUE;
  
  Int_t newClusters, newWrong;
  Double_t newChi2;
  
  if (wasLast) {   
 
    fBarrelArray->Clear();
    nClusters = nWrong = 0;
    chi2 = 0.0;
    index = 0;
    wasLast = kFALSE;
  }
  
  fBarrelTrack = (AliBarrelTrack*)(*fBarrelArray)[index++];
  ps->GetBarrelTrack(fBarrelTrack);
  
  newClusters = ps->GetNumberOfClusters() - nClusters; 
  newWrong = ps->GetNWrong() - nWrong;
  newChi2 = ps->GetChi2() - chi2;
  
  nClusters =  ps->GetNumberOfClusters();
  nWrong = ps->GetNWrong();
  chi2 = ps->GetChi2();  

  if (refPlane != fgkLastPlane) {
    fBarrelTrack->SetNClusters(newClusters, newChi2);
    fBarrelTrack->SetNWrongClusters(newWrong);
  } else {
    wasLast = kTRUE;
  } 

  fBarrelTrack->SetRefPlane(refPlane, isIn);
}

//_____________________________________________________________________

Bool_t AliTRDtracker::AdjustSector(AliTRDtrack *track) {
  //
  // Rotates the track when necessary
  //

  Double_t alpha = AliTRDgeometry::GetAlpha(); 
  Double_t y = track->GetY();
  Double_t ymax = track->GetX()*TMath::Tan(0.5*alpha);

  //Int_t ns = AliTRDgeometry::kNsect;
  //Int_t s=Int_t(track->GetAlpha()/alpha)%ns; 

  if (y > ymax) {
    //s = (s+1) % ns;
    if (!track->Rotate(alpha)) return kFALSE;
  } else if (y <-ymax) {
    //s = (s-1+ns) % ns;                           
    if (!track->Rotate(-alpha)) return kFALSE;   
  } 

  return kTRUE;
}

//_____________________________________________________________________
inline Double_t f1trd(Double_t x1,Double_t y1,
                      Double_t x2,Double_t y2,
                      Double_t x3,Double_t y3)
{
  //
  // Initial approximation of the track curvature
  //
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr);
}          

//_____________________________________________________________________
inline Double_t f2trd(Double_t x1,Double_t y1,
                      Double_t x2,Double_t y2,
                      Double_t x3,Double_t y3)
{
  //
  // Initial approximation of the track curvature times X coordinate
  // of the center of curvature
  //

  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}          

//_____________________________________________________________________
inline Double_t f3trd(Double_t x1,Double_t y1,
                      Double_t x2,Double_t y2,
                      Double_t z1,Double_t z2)
{
  //
  // Initial approximation of the tangent of the track dip angle
  //

  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}            


AliTRDcluster * AliTRDtracker::GetCluster(AliTRDtrack * track, Int_t plane, Int_t timebin){
  //
  //try to find cluster in the backup list
  //
  AliTRDcluster * cl =0;
  UInt_t *indexes = track->GetBackupIndexes();
  for (UInt_t i=0;i<kMaxTimeBinIndex;i++){
    if (indexes[i]==0) break;  
    AliTRDcluster * cli = (AliTRDcluster*)fClusters->UncheckedAt(indexes[i]);
    if (!cli) break;
    if (cli->GetLocalTimeBin()!=timebin) continue;
    Int_t iplane = fGeom->GetPlane(cli->GetDetector());
    if (iplane==plane) {
      cl = cli;
      break;
    }
  }
  return cl;
}


Int_t  AliTRDtracker::GetLastPlane(AliTRDtrack * track){
  //
  //return last updated plane
  Int_t lastplane=0;
  UInt_t *indexes = track->GetBackupIndexes();
  for (UInt_t i=0;i<kMaxTimeBinIndex;i++){
    AliTRDcluster * cli = (AliTRDcluster*)fClusters->UncheckedAt(indexes[i]);
    if (!cli) break;
    Int_t iplane = fGeom->GetPlane(cli->GetDetector());
    if (iplane>lastplane) {
      lastplane = iplane;
    }
  }
  return lastplane;
}
//___________________________________________________________________
Int_t AliTRDtracker::Clusters2Tracks(AliESD* event)
{
  //
  // Finds tracks within the TRD. The ESD event is expected to contain seeds 
  // at the outer part of the TRD. The seeds
  // are found within the TRD if fAddTRDseeds is TRUE. 
  // The tracks are propagated to the innermost time bin 
  // of the TRD and the ESD event is updated
  //

  Int_t timeBins = fTrSec[0]->GetNumberOfTimeBins();
  Float_t foundMin = fgkMinClustersInTrack * timeBins; 
  Int_t nseed = 0;
  Int_t found = 0;
  Int_t innerTB = fTrSec[0]->GetInnerTimeBin();

  Int_t n = event->GetNumberOfTracks();
  for (Int_t i=0; i<n; i++) {
    AliESDtrack* seed=event->GetTrack(i);
    ULong_t status=seed->GetStatus();
    if ( (status & AliESDtrack::kTRDout ) == 0 ) continue;
    if ( (status & AliESDtrack::kTRDin) != 0 ) continue;
    nseed++;

    AliTRDtrack* seed2 = new AliTRDtrack(*seed);
    //seed2->ResetCovariance(); 
    AliTRDtrack *pt = new AliTRDtrack(*seed2,seed2->GetAlpha());
    AliTRDtrack &t=*pt; 
    FollowProlongation(t, innerTB); 
    if (t.GetNumberOfClusters() >= foundMin) {
      UseClusters(&t);
      CookLabel(pt, 1-fgkLabelFraction);
      //      t.CookdEdx();
    }
    found++;
//    cout<<found<<'\r';     

    if(PropagateToTPC(t)) {
      seed->UpdateTrackParams(pt, AliESDtrack::kTRDin);
    }  
    delete seed2;
    delete pt;
  }     

  cout<<"Number of loaded seeds: "<<nseed<<endl;  
  cout<<"Number of found tracks from loaded seeds: "<<found<<endl;

  // after tracks from loaded seeds are found and the corresponding 
  // clusters are used, look for additional seeds from TRD

  if(fAddTRDseeds) { 
    // Find tracks for the seeds in the TRD
    Int_t timeBins = fTrSec[0]->GetNumberOfTimeBins();
  
    Int_t nSteps = (Int_t) (fgkSeedDepth / fgkSeedStep);
    Int_t gap = (Int_t) (timeBins * fgkSeedGap);
    Int_t step = (Int_t) (timeBins * fgkSeedStep);
  
    // make a first turn with tight cut on initial curvature
    for(Int_t turn = 1; turn <= 2; turn++) {
      if(turn == 2) {
        nSteps = (Int_t) (fgkSeedDepth / (3*fgkSeedStep));
        step = (Int_t) (timeBins * (3*fgkSeedStep));
      }
      for(Int_t i=0; i<nSteps; i++) {
        Int_t outer=timeBins-1-i*step; 
        Int_t inner=outer-gap;

        nseed=fSeeds->GetEntriesFast();
      
        MakeSeeds(inner, outer, turn);
      
        nseed=fSeeds->GetEntriesFast();
	//        printf("\n turn %d, step %d: number of seeds for TRD inward %d\n", 
	//               turn, i, nseed); 
              
        for (Int_t i=0; i<nseed; i++) {   
          AliTRDtrack *pt=(AliTRDtrack*)fSeeds->UncheckedAt(i), &t=*pt; 
          FollowProlongation(t,innerTB); 
          if (t.GetNumberOfClusters() >= foundMin) {
            UseClusters(&t);
            CookLabel(pt, 1-fgkLabelFraction);
            t.CookdEdx();
	    found++;
//            cout<<found<<'\r';     
            if(PropagateToTPC(t)) {
	      AliESDtrack track;
	      track.UpdateTrackParams(pt,AliESDtrack::kTRDin);
	      event->AddTrack(&track);
	      //	      track.SetTRDtrack(new AliTRDtrack(*pt));
            }        
          }
          delete fSeeds->RemoveAt(i);
          fNseeds--;
        }
      }
    }
  }
  
  cout<<"Total number of found tracks: "<<found<<endl;
    
  return 0;    
}     
     
  

//_____________________________________________________________________________
Int_t AliTRDtracker::PropagateBack(AliESD* event) {
  //
  // Gets seeds from ESD event. The seeds are AliTPCtrack's found and
  // backpropagated by the TPC tracker. Each seed is first propagated 
  // to the TRD, and then its prolongation is searched in the TRD.
  // If sufficiently long continuation of the track is found in the TRD
  // the track is updated, otherwise it's stored as originaly defined 
  // by the TPC tracker.   
  //  

  Int_t found=0;  
  Float_t foundMin = 20;

  Int_t n = event->GetNumberOfTracks();
  for (Int_t i=0; i<n; i++) {
    AliESDtrack* seed=event->GetTrack(i);
    ULong_t status=seed->GetStatus();
    if ( (status & AliESDtrack::kTPCout ) == 0 ) continue;
    if ( (status & AliESDtrack::kTRDout) != 0 ) continue;

    Int_t lbl = seed->GetLabel();
    AliTRDtrack *track = new AliTRDtrack(*seed);
    track->SetSeedLabel(lbl);
    seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup); //make backup
    fNseeds++;
    Float_t p4 = track->GetC();
    //
    Int_t expectedClr = FollowBackProlongation(*track);
    /*
      // only debug purpose
    if (track->GetNumberOfClusters()<expectedClr/3){
      AliTRDtrack *track1 = new AliTRDtrack(*seed);
      track1->SetSeedLabel(lbl);
      FollowBackProlongation(*track1);
      AliTRDtrack *track2= new AliTRDtrack(*seed);
      track->SetSeedLabel(lbl);
      FollowBackProlongation(*track2);      
      delete track1;
      delete track2;
    }
    */
     if (TMath::Abs(track->GetC()-p4)/TMath::Abs(p4)>0.2) {
      delete track;
      continue; //too big change of curvature - to be checked
    }

    Int_t foundClr = track->GetNumberOfClusters();
    if (foundClr >= foundMin) {
      track->CookdEdx(0.,1.);
      CookdEdxTimBin(*track);
 
      CookLabel(track, 1-fgkLabelFraction);
      if(track->GetChi2()/track->GetNumberOfClusters()<6) {   // sign only gold tracks
	UseClusters(track);
      }
      Bool_t isGold = kFALSE;
     
      if (track->GetChi2()/track->GetNumberOfClusters()<5) {  //full gold track
	seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	isGold = kTRUE;
      }
      if (!isGold && track->GetNCross()==0&&track->GetChi2()/track->GetNumberOfClusters()<7){ //almost gold track
	seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	isGold = kTRUE;
      }
      if (!isGold && track->GetBackupTrack()){
	if (track->GetBackupTrack()->GetNumberOfClusters()>foundMin&&
	    (track->GetBackupTrack()->GetChi2()/(track->GetBackupTrack()->GetNumberOfClusters()+1))<7){	  
	  seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
	  isGold = kTRUE;
	}
      }
    }
    else{
      delete track;
      continue;
    }

    
    if (track->GetStop()==kFALSE){

      Double_t xtof=371.;
      Double_t c2=track->GetC()*xtof - track->GetEta();
      if (TMath::Abs(c2)>=0.85) {
	delete track;
	continue;
      }
      Double_t xTOF0 = 371. ;          
      PropagateToOuterPlane(*track,xTOF0); 
      //      
      Double_t ymax=xtof*TMath::Tan(0.5*AliTRDgeometry::GetAlpha());
      Double_t y=track->GetYat(xtof);
      if (y > ymax) {
	if (!track->Rotate(AliTRDgeometry::GetAlpha())) {
	  delete track;
	  continue;
	}
      } else if (y <-ymax) {
	if (!track->Rotate(-AliTRDgeometry::GetAlpha())) {
	  delete track;
	  continue;
	}
      }
      
      if (track->PropagateTo(xtof)) {
	seed->UpdateTrackParams(track, AliESDtrack::kTRDout);
        for (Int_t i=0;i<kNPlane;i++) {
           seed->SetTRDsignals(track->GetPIDsignals(i),i);
           seed->SetTRDTimBin(track->GetPIDTimBin(i),i);
        }
	seed->SetTRDtrack(new AliTRDtrack(*track));
	if (track->GetNumberOfClusters()>foundMin) found++;
      }
    }else{
      if (track->GetNumberOfClusters()>15&&track->GetNumberOfClusters()>0.5*expectedClr){
	seed->UpdateTrackParams(track, AliESDtrack::kTRDout);
	//seed->SetStatus(AliESDtrack::kTRDStop);    
        for (Int_t i=0;i<kNPlane;i++) {
           seed->SetTRDsignals(track->GetPIDsignals(i),i);
           seed->SetTRDTimBin(track->GetPIDTimBin(i),i);
        }
	seed->SetTRDtrack(new AliTRDtrack(*track));
	found++;
      }
    }

    delete track;
    
    //End of propagation to the TOF
    //if (foundClr>foundMin)
    //  seed->UpdateTrackParams(track, AliESDtrack::kTRDout);
    

  }
  
  cerr<<"Number of seeds: "<<fNseeds<<endl;  
  cerr<<"Number of back propagated TRD tracks: "<<found<<endl;

  fSeeds->Clear(); fNseeds=0;

  return 0;

}

//_____________________________________________________________________________
Int_t AliTRDtracker::RefitInward(AliESD* event)
{
  //
  // Refits tracks within the TRD. The ESD event is expected to contain seeds 
  // at the outer part of the TRD. 
  // The tracks are propagated to the innermost time bin 
  // of the TRD and the ESD event is updated
  // Origin: Thomas KUHR (Thomas.Kuhr@cern.ch)
  //

  Int_t timeBins = fTrSec[0]->GetNumberOfTimeBins();
  Float_t foundMin = fgkMinClustersInTrack * timeBins; 
  Int_t nseed = 0;
  Int_t found = 0;
  Int_t innerTB = fTrSec[0]->GetInnerTimeBin();

  Int_t n = event->GetNumberOfTracks();
  for (Int_t i=0; i<n; i++) {
    AliESDtrack* seed=event->GetTrack(i);
    AliTRDtrack* seed2 = new AliTRDtrack(*seed);
    if (seed2->GetX()<270){
      seed->UpdateTrackParams(seed2, AliESDtrack::kTRDbackup); // backup TPC track - only update
      delete seed2;
      continue;
    }

    ULong_t status=seed->GetStatus();
    if ( (status & AliESDtrack::kTRDout ) == 0 ) {
      delete seed2;
      continue;
    }
    if ( (status & AliESDtrack::kTRDin) != 0 ) {
      delete seed2;
      continue;
    }
    nseed++;    
    seed2->ResetCovariance(5.); 
    AliTRDtrack *pt = new AliTRDtrack(*seed2,seed2->GetAlpha());
    for (Int_t i=0;i<kNPlane;i++) {
        pt->SetPIDsignals(seed2->GetPIDsignals(i),i);
        pt->SetPIDTimBin(seed2->GetPIDTimBin(i),i);
    }

    UInt_t * indexes2 = seed2->GetIndexes();
    UInt_t * indexes3 = pt->GetBackupIndexes();
    for (Int_t i=0;i<200;i++) {
      if (indexes2[i]==0) break;
      indexes3[i] = indexes2[i];
    }          
    //AliTRDtrack *pt = seed2;
    AliTRDtrack &t=*pt; 
    FollowProlongation(t, innerTB); 
    /*
    if (t.GetNumberOfClusters()<seed->GetTRDclusters(indexes3)*0.5){
      // debug  - why we dont go back?
      AliTRDtrack *pt2 = new AliTRDtrack(*seed2,seed2->GetAlpha());
      UInt_t * indexes2 = seed2->GetIndexes();
      UInt_t * indexes3 = pt2->GetBackupIndexes();
      for (Int_t i=0;i<200;i++) {
	if (indexes2[i]==0) break;
	indexes3[i] = indexes2[i];
      }  
      FollowProlongation(*pt2, innerTB);
      delete pt2;
    }
    */
    if (t.GetNumberOfClusters() >= foundMin) {
      //      UseClusters(&t);
      //CookLabel(pt, 1-fgkLabelFraction);
      //      t.CookdEdx();
    }
    found++;
//    cout<<found<<'\r';     

    if(PropagateToTPC(t)) {
      seed->UpdateTrackParams(pt, AliESDtrack::kTRDrefit);
      for (Int_t i=0;i<kNPlane;i++) {
        seed->SetTRDsignals(pt->GetPIDsignals(i),i);
        seed->SetTRDTimBin(pt->GetPIDTimBin(i),i);
      }
    }else{
      //if not prolongation to TPC - propagate without update
      AliTRDtrack* seed2 = new AliTRDtrack(*seed);
      seed2->ResetCovariance(5.); 
      AliTRDtrack *pt2 = new AliTRDtrack(*seed2,seed2->GetAlpha());
      delete seed2;
      if (PropagateToTPC(*pt2)) { 
        pt2->CookdEdx(0.,1.);
        CookdEdxTimBin(*pt2);
	seed->UpdateTrackParams(pt2, AliESDtrack::kTRDrefit);
        for (Int_t i=0;i<kNPlane;i++) {
          seed->SetTRDsignals(seed2->GetPIDsignals(i),i);
          seed->SetTRDTimBin(seed2->GetPIDTimBin(i),i);
        }
      }
      delete pt2;
    }  

    delete seed2;
    delete pt;
  }   

  cout<<"Number of loaded seeds: "<<nseed<<endl;  
  cout<<"Number of found tracks from loaded seeds: "<<found<<endl;

  return 0;

}


//---------------------------------------------------------------------------
Int_t AliTRDtracker::FollowProlongation(AliTRDtrack& t, Int_t rf)
{
  // Starting from current position on track=t this function tries
  // to extrapolate the track up to timeBin=0 and to confirm prolongation
  // if a close cluster is found. Returns the number of clusters
  // expected to be found in sensitive layers

  Float_t  wIndex, wTB, wChi2;
  Float_t  wYrt, wYclosest, wYcorrect, wYwindow;
  Float_t  wZrt, wZclosest, wZcorrect, wZwindow;
  Float_t  wPx, wPy, wPz, wC;
  Double_t px, py, pz;
  Float_t  wSigmaC2, wSigmaTgl2, wSigmaY2, wSigmaZ2;
  Int_t lastplane = GetLastPlane(&t);

  Int_t trackIndex = t.GetLabel();  

  Int_t ns=Int_t(2*TMath::Pi()/AliTRDgeometry::GetAlpha()+0.5);     

  Int_t tryAgain=fMaxGap;

  Double_t alpha=t.GetAlpha();
  alpha = TVector2::Phi_0_2pi(alpha);

  Int_t s=Int_t(alpha/AliTRDgeometry::GetAlpha())%AliTRDgeometry::kNsect;  
  Double_t radLength, rho, x, dx, y, ymax, z;

  Int_t expectedNumberOfClusters = 0;
  Bool_t lookForCluster;

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning

 
  for (Int_t nr=fTrSec[0]->GetLayerNumber(t.GetX()); nr>rf; nr--) { 

    y = t.GetY(); z = t.GetZ();

    // first propagate to the inner surface of the current time bin 
    fTrSec[s]->GetLayer(nr)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr)->GetX()-dx/2; y = t.GetY(); z = t.GetZ();
    if(!t.PropagateTo(x,radLength,rho)) break;
    y = t.GetY();
    ymax = x*TMath::Tan(0.5*alpha);
    if (y > ymax) {
      s = (s+1) % ns;
      if (!t.Rotate(alpha)) break;
      if(!t.PropagateTo(x,radLength,rho)) break;
    } else if (y <-ymax) {
      s = (s-1+ns) % ns;                           
      if (!t.Rotate(-alpha)) break;   
      if(!t.PropagateTo(x,radLength,rho)) break;
    } 

    y = t.GetY(); z = t.GetZ();

    // now propagate to the middle plane of the next time bin 
    fTrSec[s]->GetLayer(nr-1)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr-1)->GetX(); y = t.GetY(); z = t.GetZ();
    if(!t.PropagateTo(x,radLength,rho)) break;
    y = t.GetY();
    ymax = x*TMath::Tan(0.5*alpha);
    if (y > ymax) {
      s = (s+1) % ns;
      if (!t.Rotate(alpha)) break;
      if(!t.PropagateTo(x,radLength,rho)) break;
    } else if (y <-ymax) {
      s = (s-1+ns) % ns;                           
      if (!t.Rotate(-alpha)) break;   
      if(!t.PropagateTo(x,radLength,rho)) break;
    } 


    if(lookForCluster) {

      expectedNumberOfClusters++;       
      wIndex = (Float_t) t.GetLabel();
      wTB = nr;

      AliTRDpropagationLayer& timeBin=*(fTrSec[s]->GetLayer(nr-1));

      Double_t sy2=ExpectedSigmaY2(x,t.GetTgl(),t.GetPt());
      Double_t sz2=ExpectedSigmaZ2(x,t.GetTgl());

      Double_t road;
      if((t.GetSigmaY2() + sy2) > 0) road=10.*sqrt(t.GetSigmaY2() + sy2);
      else return expectedNumberOfClusters;
      
      wYrt = (Float_t) y;
      wZrt = (Float_t) z;
      wYwindow = (Float_t) road;
      t.GetPxPyPz(px,py,pz);
      wPx = (Float_t) px;
      wPy = (Float_t) py;
      wPz = (Float_t) pz;
      wC  = (Float_t) t.GetC();
      wSigmaC2 = (Float_t) t.GetSigmaC2();
      wSigmaTgl2    = (Float_t) t.GetSigmaTgl2();
      wSigmaY2 = (Float_t) t.GetSigmaY2();
      wSigmaZ2 = (Float_t) t.GetSigmaZ2();
      wChi2 = -1;            
      

      AliTRDcluster *cl=0;
      UInt_t index=0;

      Double_t maxChi2=fgkMaxChi2;

      wYclosest = 12345678;
      wYcorrect = 12345678;
      wZclosest = 12345678;
      wZcorrect = 12345678;
      wZwindow  = TMath::Sqrt(2.25 * 12 * sz2);   

      // Find the closest correct cluster for debugging purposes
      if (timeBin) {
        Float_t minDY = 1000000;
        for (Int_t i=0; i<timeBin; i++) {
          AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
          if((c->GetLabel(0) != trackIndex) &&
             (c->GetLabel(1) != trackIndex) &&
             (c->GetLabel(2) != trackIndex)) continue;
          if(TMath::Abs(c->GetY() - y) > minDY) continue;
          minDY = TMath::Abs(c->GetY() - y);
          wYcorrect = c->GetY();
          wZcorrect = c->GetZ();

          Double_t h01 = GetTiltFactor(c);
          wChi2 = t.GetPredictedChi2(c, h01);
        }
      }                    

      // Now go for the real cluster search

      if (timeBin) {
	//
	//find cluster in history
	cl =0;
	
	AliTRDcluster * cl0 = timeBin[0];
	if (!cl0) {
	  continue;
	}
	Int_t plane = fGeom->GetPlane(cl0->GetDetector());
	if (plane>lastplane) continue;
	Int_t timebin = cl0->GetLocalTimeBin();
	AliTRDcluster * cl2= GetCluster(&t,plane, timebin);
	if (cl2) {
	  cl =cl2;	
	  Double_t h01 = GetTiltFactor(cl);
	  maxChi2=t.GetPredictedChi2(cl,h01);
	}
	if ((!cl) && road>fgkWideRoad) {
	  //if (t.GetNumberOfClusters()>4)
	  //  cerr<<t.GetNumberOfClusters()
	  //	<<"FindProlongation warning: Too broad road !\n";
	  continue;
	}             


	if(!cl){

	  for (Int_t i=timeBin.Find(y-road); i<timeBin; i++) {
	    AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
	    if (c->GetY() > y+road) break;
	    if (c->IsUsed() > 0) continue;
	    if((c->GetZ()-z)*(c->GetZ()-z) > 3 * sz2) continue;
	    
	    Double_t h01 = GetTiltFactor(c);
	    Double_t chi2=t.GetPredictedChi2(c,h01);
	    
	    if (chi2 > maxChi2) continue;
	    maxChi2=chi2;
	    cl=c;
	    index=timeBin.GetIndex(i);
	  }               
	}

        if(!cl) {

          for (Int_t i=timeBin.Find(y-road); i<timeBin; i++) {
            AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
            
            if (c->GetY() > y+road) break;
            if (c->IsUsed() > 0) continue;
            if((c->GetZ()-z)*(c->GetZ()-z) > 12 * sz2) continue;
            
            Double_t h01 = GetTiltFactor(c);
            Double_t chi2=t.GetPredictedChi2(c, h01);
            
            if (chi2 > maxChi2) continue;
            maxChi2=chi2;
            cl=c;
            index=timeBin.GetIndex(i);
          }
        }        
        if (cl) {
          wYclosest = cl->GetY();
          wZclosest = cl->GetZ();
          Double_t h01 = GetTiltFactor(cl);

          t.SetSampledEdx(cl->GetQ()/dx,t.GetNumberOfClusters()); 
	  //printf("Track   position\t%f\t%f\t%f\n",t.GetX(),t.GetY(),t.GetZ());
	  //printf("Cluster position\t%d\t%f\t%f\n",cl->GetLocalTimeBin(),cl->GetY(),cl->GetZ());
	  Int_t det = cl->GetDetector();    
	  Int_t plane = fGeom->GetPlane(det);

	  if(!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
	    //if(!t.Update(cl,maxChi2,index,h01)) {
            //if(!tryAgain--) return 0;
          }  
          else tryAgain=fMaxGap;
        }
        else {
          //if (tryAgain==0) break; 
          tryAgain--;
        }

        /*
        if((((Int_t) wTB)%15 == 0) || (((Int_t) wTB)%15 == 14)) {
          
          printf(" %f", wIndex);       //1
          printf(" %f", wTB);          //2
          printf(" %f", wYrt);         //3
          printf(" %f", wYclosest);    //4
          printf(" %f", wYcorrect);    //5
          printf(" %f", wYwindow);     //6
          printf(" %f", wZrt);         //7
          printf(" %f", wZclosest);    //8
          printf(" %f", wZcorrect);    //9
          printf(" %f", wZwindow);     //10
          printf(" %f", wPx);          //11
          printf(" %f", wPy);          //12
          printf(" %f", wPz);          //13
          printf(" %f", wSigmaC2*1000000);  //14
          printf(" %f", wSigmaTgl2*1000);   //15
          printf(" %f", wSigmaY2);     //16
          //      printf(" %f", wSigmaZ2);     //17
          printf(" %f", wChi2);     //17
          printf(" %f", wC);           //18
          printf("\n");
        } 
        */                        
      }
    }  
  }
  return expectedNumberOfClusters;
  
  
}                

//___________________________________________________________________

Int_t AliTRDtracker::FollowBackProlongation(AliTRDtrack& t)
{
  // Starting from current radial position of track <t> this function
  // extrapolates the track up to outer timebin and in the sensitive
  // layers confirms prolongation if a close cluster is found. 
  // Returns the number of clusters expected to be found in sensitive layers


  Float_t  wIndex, wTB, wChi2;
  Float_t  wYrt, wYclosest, wYcorrect, wYwindow;
  Float_t  wZrt, wZclosest, wZcorrect, wZwindow;
  Float_t  wPx, wPy, wPz, wC;
  Double_t px, py, pz;
  Float_t  wSigmaC2, wSigmaTgl2, wSigmaY2, wSigmaZ2;

  Int_t trackIndex = t.GetLabel();  
  Int_t tryAgain=fMaxGap;

  Double_t alpha=t.GetAlpha();
  TVector2::Phi_0_2pi(alpha);

  Int_t s;

  Int_t outerTB = fTrSec[0]->GetOuterTimeBin();
  Double_t radLength, rho, x, dx, y, ymax = 0, z;
  Bool_t lookForCluster;

  Int_t expectedNumberOfClusters = 0;
  x = t.GetX();

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning

  Int_t nRefPlane = fgkFirstPlane;
  Bool_t isNewLayer = kFALSE; 

  Double_t chi2;
  Double_t minDY;
  Int_t zone =-10;
  Int_t nr;
  for (nr=fTrSec[0]->GetLayerNumber(t.GetX()); nr<outerTB+1; nr++) { 
    
    y = t.GetY(); 
    z = t.GetZ();

    // first propagate to the outer surface of the current time bin 

    s = t.GetSector();
    fTrSec[s]->GetLayer(nr)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr)->GetX()+dx/2; 
    y = t.GetY(); 
    z = t.GetZ();

    if(!t.PropagateTo(x,radLength,rho)) break;
    //    if (!AdjustSector(&t)) break;
    //
    // MI -fix untill correct material desription will be implemented
    //
    Float_t angle =  t.GetAlpha();  // MI - if rotation - we go through the material 
    if (!AdjustSector(&t)) break;
    Int_t cross = kFALSE;
    
    if (TMath::Abs(angle -  t.GetAlpha())>0.000001) cross = kTRUE; //better to stop track
    Int_t currentzone = fTrSec[s]->GetLayer(nr)->GetZone(z);
    if (currentzone==-10) cross = kTRUE;  // we are in the frame
    if (currentzone>-10){   // layer knows where we are
      if (zone==-10) zone = currentzone;
      if (zone!=currentzone) cross=kTRUE;  
    }
    if (cross) {
      t.IncCross();
      if (t.GetNCross()==1) t.MakeBackupTrack();
      if (t.GetNCross()>2) break;
    }
    
    //
    //
    s = t.GetSector();
    if (!t.PropagateTo(x,radLength,rho)) break;

    y = t.GetY();
    z = t.GetZ();

    // Barrel Tracks [SR, 04.04.2003]

    s = t.GetSector();
    if (fTrSec[s]->GetLayer(nr)->IsSensitive() != 
        fTrSec[s]->GetLayer(nr+1)->IsSensitive() ) {

//      if (IsStoringBarrel()) StoreBarrelTrack(&t, nRefPlane++, kTrackBack);
    }

    if (fTrSec[s]->GetLayer(nr-1)->IsSensitive() && 
          ! fTrSec[s]->GetLayer(nr)->IsSensitive()) {
      isNewLayer = kTRUE;
    } else {isNewLayer = kFALSE;}

    y = t.GetY();
    z = t.GetZ();

    // now propagate to the middle plane of the next time bin 
    fTrSec[s]->GetLayer(nr+1)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);

    x = fTrSec[s]->GetLayer(nr+1)->GetX(); 
      if(!t.PropagateTo(x,radLength,rho)) break;
    if (!AdjustSector(&t)) break;
    s = t.GetSector();
      if(!t.PropagateTo(x,radLength,rho)) break;

    y = t.GetY();
    z = t.GetZ();

    if(fVocal) printf("nr+1=%d, x %f, z %f, y %f, ymax %f\n",nr+1,x,z,y,ymax);
    //    printf("label %d, pl %d, lookForCluster %d \n",
    //     trackIndex, nr+1, lookForCluster);

    if(lookForCluster) {
      expectedNumberOfClusters++;       

      wIndex = (Float_t) t.GetLabel();
      wTB = fTrSec[s]->GetLayer(nr+1)->GetTimeBinIndex();

      AliTRDpropagationLayer& timeBin=*(fTrSec[s]->GetLayer(nr+1));
      Double_t sy2=ExpectedSigmaY2(t.GetX(),t.GetTgl(),t.GetPt());
      Double_t sz2=ExpectedSigmaZ2(t.GetX(),t.GetTgl());
      if((t.GetSigmaY2() + sy2) < 0) break;
      Double_t road = 10.*sqrt(t.GetSigmaY2() + sy2); 
      Double_t y=t.GetY(), z=t.GetZ();

      wYrt = (Float_t) y;
      wZrt = (Float_t) z;
      wYwindow = (Float_t) road;
      t.GetPxPyPz(px,py,pz);
      wPx = (Float_t) px;
      wPy = (Float_t) py;
      wPz = (Float_t) pz;
      wC  = (Float_t) t.GetC();
      wSigmaC2 = (Float_t) t.GetSigmaC2();
      wSigmaTgl2    = (Float_t) t.GetSigmaTgl2();
      wSigmaY2 = (Float_t) t.GetSigmaY2();
      wSigmaZ2 = (Float_t) t.GetSigmaZ2();
      wChi2 = -1;            
      
      if (road>fgkWideRoad) {
        if (t.GetNumberOfClusters()>4)
          cerr<<t.GetNumberOfClusters()
              <<"FindProlongation warning: Too broad road !\n";
        return 0;
      }      

      AliTRDcluster *cl=0;
      UInt_t index=0;

      Double_t maxChi2=fgkMaxChi2;

      if (isNewLayer) { 
        road = 3 * road;
        //sz2 = 3 * sz2;
        maxChi2 = 10 * fgkMaxChi2;
      }
      
      if (nRefPlane == fgkFirstPlane) maxChi2 = 20 * fgkMaxChi2; 
      if (nRefPlane == fgkFirstPlane+2) maxChi2 = 15 * fgkMaxChi2;
      if (t.GetNRotate() > 0) maxChi2 = 3 * maxChi2;
      

      wYclosest = 12345678;
      wYcorrect = 12345678;
      wZclosest = 12345678;
      wZcorrect = 12345678;
      wZwindow  = TMath::Sqrt(2.25 * 12 * sz2);   

      // Find the closest correct cluster for debugging purposes
      if (timeBin) {
        minDY = 1000000;
        for (Int_t i=0; i<timeBin; i++) {
          AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
          if((c->GetLabel(0) != trackIndex) &&
             (c->GetLabel(1) != trackIndex) &&
             (c->GetLabel(2) != trackIndex)) continue;
          if(TMath::Abs(c->GetY() - y) > minDY) continue;
          //minDY = TMath::Abs(c->GetY() - y);
          minDY = c->GetY() - y;
          wYcorrect = c->GetY();
          wZcorrect = c->GetZ();

          Double_t h01 = GetTiltFactor(c);
          wChi2 = t.GetPredictedChi2(c, h01);
        }
      }                    

      // Now go for the real cluster search

      if (timeBin) {

        for (Int_t i=timeBin.Find(y-road); i<timeBin; i++) {
          AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
          if (c->GetY() > y+road) break;
          if (c->IsUsed() > 0) continue;
          if((c->GetZ()-z)*(c->GetZ()-z) > 3 * sz2) continue;

          Double_t h01 = GetTiltFactor(c);
          chi2=t.GetPredictedChi2(c,h01);
          
          if (chi2 > maxChi2) continue;
          maxChi2=chi2;
          cl=c;
          index=timeBin.GetIndex(i);

          //check is correct
          if((c->GetLabel(0) != trackIndex) &&
             (c->GetLabel(1) != trackIndex) &&
             (c->GetLabel(2) != trackIndex)) t.AddNWrong();
        }               
	
        if(!cl) {

          for (Int_t i=timeBin.Find(y-road); i<timeBin; i++) {
            AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
            
            if (c->GetY() > y+road) break;
            if (c->IsUsed() > 0) continue;
            if((c->GetZ()-z)*(c->GetZ()-z) > 2.25 * 12 * sz2) continue;
            
            Double_t h01 = GetTiltFactor(c);
            chi2=t.GetPredictedChi2(c,h01);
            
            if (chi2 > maxChi2) continue;
            maxChi2=chi2;
            cl=c;
            index=timeBin.GetIndex(i);
          }
        }        
        
        if (cl) {
          wYclosest = cl->GetY();
          wZclosest = cl->GetZ();

          t.SetSampledEdx(cl->GetQ()/dx,t.GetNumberOfClusters()); 
          Double_t h01 = GetTiltFactor(cl);
	  Int_t det = cl->GetDetector();    
	  Int_t plane = fGeom->GetPlane(det);

	  if(!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
          //if(!t.Update(cl,maxChi2,index,h01)) {
            if(!tryAgain--) return 0;
          }  
          else tryAgain=fMaxGap;
        }
        else {
          if (tryAgain==0) break; 
          tryAgain--;
          
          //if (minDY < 1000000 && isNewLayer) 
            //cout << "\t" << nRefPlane << "\t" << "\t" << t.GetNRotate() <<  "\t" << 
            //  road << "\t" << minDY << "\t" << chi2 << "\t" << wChi2 << "\t" << maxChi2 << endl;
                                                                     
        }

        isNewLayer = kFALSE;

        /*
        if((((Int_t) wTB)%15 == 0) || (((Int_t) wTB)%15 == 14)) {
          
          printf(" %f", wIndex);       //1
          printf(" %f", wTB);          //2
          printf(" %f", wYrt);         //3
          printf(" %f", wYclosest);    //4
          printf(" %f", wYcorrect);    //5
          printf(" %f", wYwindow);     //6
          printf(" %f", wZrt);         //7
          printf(" %f", wZclosest);    //8
          printf(" %f", wZcorrect);    //9
          printf(" %f", wZwindow);     //10
          printf(" %f", wPx);          //11
          printf(" %f", wPy);          //12
          printf(" %f", wPz);          //13
          printf(" %f", wSigmaC2*1000000);  //14
          printf(" %f", wSigmaTgl2*1000);   //15
          printf(" %f", wSigmaY2);     //16
          //      printf(" %f", wSigmaZ2);     //17
          printf(" %f", wChi2);     //17
          printf(" %f", wC);           //18
          printf("\n");
        } 
        */                        
      }
    }  
  }
  if (nr<outerTB) 
    t.SetStop(kTRUE);
  else
    t.SetStop(kFALSE);
  return expectedNumberOfClusters;


}         

//---------------------------------------------------------------------------
Int_t AliTRDtracker::Refit(AliTRDtrack& t, Int_t rf)
{
  // Starting from current position on track=t this function tries
  // to extrapolate the track up to timeBin=0 and to reuse already
  // assigned clusters. Returns the number of clusters
  // expected to be found in sensitive layers
  // get indices of assigned clusters for each layer
  // Origin: Thomas KUHR (Thomas.Kuhr@cern.ch)

  Int_t iCluster[90];
  for (Int_t i = 0; i < 90; i++) iCluster[i] = 0;
  for (Int_t i = 0; i < t.GetNumberOfClusters(); i++) {
    Int_t index = t.GetClusterIndex(i);
    AliTRDcluster *cl=(AliTRDcluster*) GetCluster(index);
    if (!cl) continue;
    Int_t detector=cl->GetDetector();
    Int_t localTimeBin=cl->GetLocalTimeBin();
    Int_t sector=fGeom->GetSector(detector);
    Int_t plane=fGeom->GetPlane(detector);

    Int_t trackingSector = CookSectorIndex(sector);

    Int_t gtb = fTrSec[trackingSector]->CookTimeBinIndex(plane,localTimeBin);
    if(gtb < 0) continue; 
    Int_t layer = fTrSec[trackingSector]->GetLayerNumber(gtb);
    iCluster[layer] = index;
  }
  t.ResetClusters();

  Int_t ns=Int_t(2*TMath::Pi()/AliTRDgeometry::GetAlpha()+0.5);     

  Double_t alpha=t.GetAlpha();
  alpha = TVector2::Phi_0_2pi(alpha);

  Int_t s=Int_t(alpha/AliTRDgeometry::GetAlpha())%AliTRDgeometry::kNsect;  
  Double_t radLength, rho, x, dx, y, ymax, z;

  Int_t expectedNumberOfClusters = 0;
  Bool_t lookForCluster;

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning

 
  for (Int_t nr=fTrSec[0]->GetLayerNumber(t.GetX()); nr>rf; nr--) { 

    y = t.GetY(); z = t.GetZ();

    // first propagate to the inner surface of the current time bin 
    fTrSec[s]->GetLayer(nr)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr)->GetX()-dx/2; y = t.GetY(); z = t.GetZ();
    if(!t.PropagateTo(x,radLength,rho)) break;
    y = t.GetY();
    ymax = x*TMath::Tan(0.5*alpha);
    if (y > ymax) {
      s = (s+1) % ns;
      if (!t.Rotate(alpha)) break;
      if(!t.PropagateTo(x,radLength,rho)) break;
    } else if (y <-ymax) {
      s = (s-1+ns) % ns;                           
      if (!t.Rotate(-alpha)) break;   
      if(!t.PropagateTo(x,radLength,rho)) break;
    } 

    y = t.GetY(); z = t.GetZ();

    // now propagate to the middle plane of the next time bin 
    fTrSec[s]->GetLayer(nr-1)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr-1)->GetX(); y = t.GetY(); z = t.GetZ();
    if(!t.PropagateTo(x,radLength,rho)) break;
    y = t.GetY();
    ymax = x*TMath::Tan(0.5*alpha);
    if (y > ymax) {
      s = (s+1) % ns;
      if (!t.Rotate(alpha)) break;
      if(!t.PropagateTo(x,radLength,rho)) break;
    } else if (y <-ymax) {
      s = (s-1+ns) % ns;                           
      if (!t.Rotate(-alpha)) break;   
      if(!t.PropagateTo(x,radLength,rho)) break;
    } 

    if(lookForCluster) expectedNumberOfClusters++;       

    // use assigned cluster
    if (!iCluster[nr-1]) continue;
    AliTRDcluster *cl=(AliTRDcluster*)GetCluster(iCluster[nr-1]);
    Double_t h01 = GetTiltFactor(cl);
    Double_t chi2=t.GetPredictedChi2(cl, h01);
    t.SetSampledEdx(cl->GetQ()/dx,t.GetNumberOfClusters()); 
    t.Update(cl,chi2,iCluster[nr-1],h01);
  }

  return expectedNumberOfClusters;
}                

//___________________________________________________________________

Int_t AliTRDtracker::PropagateToOuterPlane(AliTRDtrack& t, Double_t xToGo)
{
  // Starting from current radial position of track <t> this function
  // extrapolates the track up to radial position <xToGo>. 
  // Returns 1 if track reaches the plane, and 0 otherwise 

  Int_t ns=Int_t(2*TMath::Pi()/AliTRDgeometry::GetAlpha()+0.5);     

  Double_t alpha=t.GetAlpha();

  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();

  Int_t s=Int_t(alpha/AliTRDgeometry::GetAlpha())%AliTRDgeometry::kNsect;  

  Bool_t lookForCluster;
  Double_t radLength, rho, x, dx, y, ymax, z;

  x = t.GetX();

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning

  Int_t plToGo = fTrSec[0]->GetLayerNumber(xToGo);

  for (Int_t nr=fTrSec[0]->GetLayerNumber(x); nr<plToGo; nr++) { 

    y = t.GetY(); z = t.GetZ();

    // first propagate to the outer surface of the current time bin 
    fTrSec[s]->GetLayer(nr)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr)->GetX()+dx/2; y = t.GetY(); z = t.GetZ();
    if(!t.PropagateTo(x,radLength,rho)) return 0;
    y = t.GetY();
    ymax = x*TMath::Tan(0.5*alpha);
    if (y > ymax) {
      s = (s+1) % ns;
      if (!t.Rotate(alpha)) return 0;
    } else if (y <-ymax) {
      s = (s-1+ns) % ns;                           
      if (!t.Rotate(-alpha)) return 0;   
    } 
    if(!t.PropagateTo(x,radLength,rho)) return 0;

    y = t.GetY(); z = t.GetZ();

    // now propagate to the middle plane of the next time bin 
    fTrSec[s]->GetLayer(nr+1)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr+1)->GetX(); y = t.GetY(); z = t.GetZ();
    if(!t.PropagateTo(x,radLength,rho)) return 0;
    y = t.GetY();
    ymax = x*TMath::Tan(0.5*alpha);
    if (y > ymax) {
      s = (s+1) % ns;
      if (!t.Rotate(alpha)) return 0;
    } else if (y <-ymax) {
      s = (s-1+ns) % ns;                           
      if (!t.Rotate(-alpha)) return 0;   
    } 
    if(!t.PropagateTo(x,radLength,rho)) return 0;
  }
  return 1;
}         

//___________________________________________________________________

Int_t AliTRDtracker::PropagateToTPC(AliTRDtrack& t)
{
  // Starting from current radial position of track <t> this function
  // extrapolates the track up to radial position of the outermost
  // padrow of the TPC. 
  // Returns 1 if track reaches the TPC, and 0 otherwise 

  //Int_t ns=Int_t(2*TMath::Pi()/AliTRDgeometry::GetAlpha()+0.5);     

  Double_t alpha=t.GetAlpha();
  alpha = TVector2::Phi_0_2pi(alpha);

  Int_t s=Int_t(alpha/AliTRDgeometry::GetAlpha())%AliTRDgeometry::kNsect;  

  Bool_t lookForCluster;
  Double_t radLength, rho, x, dx, y, /*ymax,*/ z;

  x = t.GetX();

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning
  Int_t plTPC = fTrSec[0]->GetLayerNumber(246.055);

  for (Int_t nr=fTrSec[0]->GetLayerNumber(x); nr>plTPC; nr--) { 

    y = t.GetY(); 
    z = t.GetZ();

    // first propagate to the outer surface of the current time bin 
    fTrSec[s]->GetLayer(nr)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr)->GetX()-dx/2; 
    
    if(!t.PropagateTo(x,radLength,rho)) return 0;
    AdjustSector(&t);
    if(!t.PropagateTo(x,radLength,rho)) return 0;

    y = t.GetY(); 
    z = t.GetZ();

    // now propagate to the middle plane of the next time bin 
    fTrSec[s]->GetLayer(nr-1)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
    x = fTrSec[s]->GetLayer(nr-1)->GetX(); 
    
    if(!t.PropagateTo(x,radLength,rho)) return 0;
    AdjustSector(&t);
    if(!t.PropagateTo(x,radLength,rho)) return 0;
  } 
  return 1;
}         

//_____________________________________________________________________________
Int_t AliTRDtracker::LoadClusters(TTree *cTree)
{
  // Fills clusters into TRD tracking_sectors 
  // Note that the numbering scheme for the TRD tracking_sectors 
  // differs from that of TRD sectors

  if (ReadClusters(fClusters,cTree)) {
     Error("LoadClusters","Problem with reading the clusters !");
     return 1;
  }
  Int_t ncl=fClusters->GetEntriesFast();
  fNclusters=ncl;
  cout<<"\n LoadSectors: sorting "<<ncl<<" clusters"<<endl;
              
  UInt_t index;
  for (Int_t ichamber=0;ichamber<5;ichamber++)
    for (Int_t isector=0;isector<18;isector++){
      fHoles[ichamber][isector]=kTRUE;
    }


  while (ncl--) {
//    printf("\r %d left  ",ncl); 
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(ncl);
    Int_t detector=c->GetDetector();
    Int_t localTimeBin=c->GetLocalTimeBin();
    Int_t sector=fGeom->GetSector(detector);
    Int_t plane=fGeom->GetPlane(detector);
      
    Int_t trackingSector = CookSectorIndex(sector);
    if (c->GetLabel(0)>0){
      Int_t chamber = fGeom->GetChamber(detector);
      fHoles[chamber][trackingSector]=kFALSE;
    }

    Int_t gtb = fTrSec[trackingSector]->CookTimeBinIndex(plane,localTimeBin);
    if(gtb < 0) continue; 
    Int_t layer = fTrSec[trackingSector]->GetLayerNumber(gtb);

    index=ncl;
    fTrSec[trackingSector]->GetLayer(layer)->InsertCluster(c,index);
  }    
  //  printf("\r\n");
  //
  //
  /*
  for (Int_t isector=0;isector<18;isector++){
    for (Int_t ichamber=0;ichamber<5;ichamber++)      
      if (fHoles[ichamber][isector]!=fGeom->IsHole(0,ichamber,17-isector)) 
	printf("Problem \t%d\t%d\t%d\t%d\n",isector,ichamber,fHoles[ichamber][isector],
	     fGeom->IsHole(0,ichamber,17-isector));
  }
  */
  return 0;
}

//_____________________________________________________________________________
void AliTRDtracker::UnloadClusters() 
{ 
  //
  // Clears the arrays of clusters and tracks. Resets sectors and timebins 
  //

  Int_t i, nentr;

  nentr = fClusters->GetEntriesFast();
  for (i = 0; i < nentr; i++) delete fClusters->RemoveAt(i);
  fNclusters = 0;

  nentr = fSeeds->GetEntriesFast();
  for (i = 0; i < nentr; i++) delete fSeeds->RemoveAt(i);

  nentr = fTracks->GetEntriesFast();
  for (i = 0; i < nentr; i++) delete fTracks->RemoveAt(i);

  Int_t nsec = AliTRDgeometry::kNsect;

  for (i = 0; i < nsec; i++) {    
    for(Int_t pl = 0; pl < fTrSec[i]->GetNumberOfLayers(); pl++) {
      fTrSec[i]->GetLayer(pl)->Clear();
    }
  }

}

//__________________________________________________________________________
void AliTRDtracker::MakeSeeds(Int_t inner, Int_t outer, Int_t turn)
{
  // Creates track seeds using clusters in timeBins=i1,i2

  if(turn > 2) {
    cerr<<"MakeSeeds: turn "<<turn<<" exceeds the limit of 2"<<endl;
    return;
  }

  Double_t x[5], c[15];
  Int_t maxSec=AliTRDgeometry::kNsect;
  
  Double_t alpha=AliTRDgeometry::GetAlpha();
  Double_t shift=AliTRDgeometry::GetAlpha()/2.;
  Double_t cs=cos(alpha), sn=sin(alpha);
  Double_t cs2=cos(2.*alpha), sn2=sin(2.*alpha);
    
      
  Int_t i2 = fTrSec[0]->GetLayerNumber(inner);
  Int_t i1 = fTrSec[0]->GetLayerNumber(outer);
      
  Double_t x1 =fTrSec[0]->GetX(i1);
  Double_t xx2=fTrSec[0]->GetX(i2);
      
  for (Int_t ns=0; ns<maxSec; ns++) {
    
    Int_t nl2 = *(fTrSec[(ns-2+maxSec)%maxSec]->GetLayer(i2));
    Int_t nl=(*fTrSec[(ns-1+maxSec)%maxSec]->GetLayer(i2));
    Int_t nm=(*fTrSec[ns]->GetLayer(i2));
    Int_t nu=(*fTrSec[(ns+1)%maxSec]->GetLayer(i2));
    Int_t nu2=(*fTrSec[(ns+2)%maxSec]->GetLayer(i2));
    
    AliTRDpropagationLayer& r1=*(fTrSec[ns]->GetLayer(i1));
    
    for (Int_t is=0; is < r1; is++) {
      Double_t y1=r1[is]->GetY(), z1=r1[is]->GetZ();
      
      for (Int_t js=0; js < nl2+nl+nm+nu+nu2; js++) {
        
        const AliTRDcluster *cl;
        Double_t x2,   y2,   z2;
        Double_t x3=0., y3=0.;   
        
        if (js<nl2) {
          if(turn != 2) continue;
          AliTRDpropagationLayer& r2=*(fTrSec[(ns-2+maxSec)%maxSec]->GetLayer(i2));
          cl=r2[js];
          y2=cl->GetY(); z2=cl->GetZ();
          
          x2= xx2*cs2+y2*sn2;
          y2=-xx2*sn2+y2*cs2;
        }
        else if (js<nl2+nl) {
          if(turn != 1) continue;
          AliTRDpropagationLayer& r2=*(fTrSec[(ns-1+maxSec)%maxSec]->GetLayer(i2));
          cl=r2[js-nl2];
          y2=cl->GetY(); z2=cl->GetZ();
          
          x2= xx2*cs+y2*sn;
          y2=-xx2*sn+y2*cs;
        }                                
        else if (js<nl2+nl+nm) {
          if(turn != 1) continue;
          AliTRDpropagationLayer& r2=*(fTrSec[ns]->GetLayer(i2));
          cl=r2[js-nl2-nl];
          x2=xx2; y2=cl->GetY(); z2=cl->GetZ();
        }
        else if (js<nl2+nl+nm+nu) {
          if(turn != 1) continue;
          AliTRDpropagationLayer& r2=*(fTrSec[(ns+1)%maxSec]->GetLayer(i2));
          cl=r2[js-nl2-nl-nm];
          y2=cl->GetY(); z2=cl->GetZ();
          
          x2=xx2*cs-y2*sn;
          y2=xx2*sn+y2*cs;
        }              
        else {
          if(turn != 2) continue;
          AliTRDpropagationLayer& r2=*(fTrSec[(ns+2)%maxSec]->GetLayer(i2));
          cl=r2[js-nl2-nl-nm-nu];
          y2=cl->GetY(); z2=cl->GetZ();
          
          x2=xx2*cs2-y2*sn2;
          y2=xx2*sn2+y2*cs2;
        }
        
        if(TMath::Abs(z1-z2) > fgkMaxSeedDeltaZ12) continue;
        
        Double_t zz=z1 - z1/x1*(x1-x2);
        
        if (TMath::Abs(zz-z2)>fgkMaxSeedDeltaZ) continue;
        
        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {cerr<<"TRD MakeSeeds: Straight seed !\n"; continue;}
        
        x[0]=y1;
        x[1]=z1;
        x[4]=f1trd(x1,y1,x2,y2,x3,y3);
        
        if (TMath::Abs(x[4]) > fgkMaxSeedC) continue;      
        
        x[2]=f2trd(x1,y1,x2,y2,x3,y3);
        
        if (TMath::Abs(x[4]*x1-x[2]) >= 0.99999) continue;
        
        x[3]=f3trd(x1,y1,x2,y2,z1,z2);
        
        if (TMath::Abs(x[3]) > fgkMaxSeedTan) continue;
        
        Double_t a=asin(x[2]);
        Double_t zv=z1 - x[3]/x[4]*(a+asin(x[4]*x1-x[2]));
        
        if (TMath::Abs(zv)>fgkMaxSeedVertexZ) continue;
        
        Double_t sy1=r1[is]->GetSigmaY2(), sz1=r1[is]->GetSigmaZ2();
        Double_t sy2=cl->GetSigmaY2(),     sz2=cl->GetSigmaZ2();
        Double_t sy3=fgkSeedErrorSY3, sy=fgkSeedErrorSY, sz=fgkSeedErrorSZ;  

        // Tilt changes
        Double_t h01 = GetTiltFactor(r1[is]);
	Double_t xuFactor = 100.;
	if(fNoTilt) { 
	  h01 = 0;
	  xuFactor = 1;
	}

        sy1=sy1+sz1*h01*h01;
        Double_t syz=sz1*(-h01);
        // end of tilt changes
        
        Double_t f40=(f1trd(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
        Double_t f42=(f1trd(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
        Double_t f43=(f1trd(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
        Double_t f20=(f2trd(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
        Double_t f22=(f2trd(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
        Double_t f23=(f2trd(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
        Double_t f30=(f3trd(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
        Double_t f31=(f3trd(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
        Double_t f32=(f3trd(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
        Double_t f34=(f3trd(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;    

        
        c[0]=sy1;
        //        c[1]=0.;       c[2]=sz1;
        c[1]=syz;       c[2]=sz1*xuFactor;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
        c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
                       c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
        c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
        c[13]=f30*sy1*f40+f32*sy2*f42;
        c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;      
        
        UInt_t index=r1.GetIndex(is);
        
        AliTRDtrack *track=new AliTRDtrack(r1[is],index,x,c,x1,ns*alpha+shift);

        Int_t rc=FollowProlongation(*track, i2);     
        
        if ((rc < 1) ||
            (track->GetNumberOfClusters() < 
             (outer-inner)*fgkMinClustersInSeed)) delete track;
        else {
          fSeeds->AddLast(track); fNseeds++;
//          cerr<<"\r found seed "<<fNseeds;
        }
      }
    }
  }
}            

//_____________________________________________________________________________
Int_t AliTRDtracker::ReadClusters(TObjArray *array, TTree *ClusterTree) const
{
  //
  // Reads AliTRDclusters (option >= 0) or AliTRDrecPoints (option < 0) 
  // from the file. The names of the cluster tree and branches 
  // should match the ones used in AliTRDclusterizer::WriteClusters()
  //
  TObjArray *clusterArray = new TObjArray(400); 
  
  TBranch *branch=ClusterTree->GetBranch("TRDcluster");
  if (!branch) {
    Error("ReadClusters","Can't get the branch !");
    return 1;
  }
  branch->SetAddress(&clusterArray); 
  
  Int_t nEntries = (Int_t) ClusterTree->GetEntries();
  //  printf("found %d entries in %s.\n",nEntries,ClusterTree->GetName());
  
  // Loop through all entries in the tree
  Int_t nbytes = 0;
  AliTRDcluster *c = 0;
  //  printf("\n");

  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    
    
    // Import the tree
    nbytes += ClusterTree->GetEvent(iEntry);  
    
    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
//    printf("\r Read %d clusters from entry %d", nCluster, iEntry);
    
    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (AliTRDcluster*)clusterArray->UncheckedAt(iCluster);
      AliTRDcluster *co = new AliTRDcluster(*c);
      co->SetSigmaY2(c->GetSigmaY2() * fSY2corr);
      Int_t ltb = co->GetLocalTimeBin();
      if(ltb == 19) co->SetSigmaZ2(c->GetSigmaZ2());
      else if(fNoTilt) co->SetSigmaZ2(c->GetSigmaZ2() * fSZ2corr);
      array->AddLast(co);
      delete clusterArray->RemoveAt(iCluster); 
    }
  }

  delete clusterArray;

  return 0;
}

//__________________________________________________________________
void AliTRDtracker::CookLabel(AliKalmanTrack* pt, Float_t wrong) const 
{
  //
  // This cooks a label. Mmmmh, smells good...
  //

  Int_t label=123456789, index, i, j;
  Int_t ncl=pt->GetNumberOfClusters();
  const Int_t kRange = fTrSec[0]->GetOuterTimeBin()+1;

  Bool_t labelAdded;

  //  Int_t s[kRange][2];
  Int_t **s = new Int_t* [kRange];
  for (i=0; i<kRange; i++) {
    s[i] = new Int_t[2];
  }
  for (i=0; i<kRange; i++) {
    s[i][0]=-1;
    s[i][1]=0;
  }

  Int_t t0,t1,t2;
  for (i=0; i<ncl; i++) {
    index=pt->GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    t0=c->GetLabel(0);
    t1=c->GetLabel(1);
    t2=c->GetLabel(2);
  }

  for (i=0; i<ncl; i++) {
    index=pt->GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    for (Int_t k=0; k<3; k++) { 
      label=c->GetLabel(k);
      labelAdded=kFALSE; j=0;
      if (label >= 0) {
        while ( (!labelAdded) && ( j < kRange ) ) {
          if (s[j][0]==label || s[j][1]==0) {
            s[j][0]=label; 
            s[j][1]=s[j][1]+1; 
            labelAdded=kTRUE;
          }
          j++;
        }
      }
    }
  }

  Int_t max=0;
  label = -123456789;

  for (i=0; i<kRange; i++) {
    if (s[i][1]>max) {
      max=s[i][1]; label=s[i][0];
    }
  }

  for (i=0; i<kRange; i++) {
    delete []s[i];
  }        

  delete []s;

  if ((1.- Float_t(max)/ncl) > wrong) label=-label;   

  pt->SetLabel(label); 

}


//__________________________________________________________________
void AliTRDtracker::UseClusters(const AliKalmanTrack* t, Int_t from) const 
{
  //
  // Use clusters, but don't abuse them!
  //

  Int_t ncl=t->GetNumberOfClusters();
  for (Int_t i=from; i<ncl; i++) {
    Int_t index = t->GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    c->Use();
  }
}


//_____________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaY2(Double_t , Double_t , Double_t ) const
{
  // Parametrised "expected" error of the cluster reconstruction in Y 

  Double_t s = 0.08 * 0.08;    
  return s;
}

//_____________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaZ2(Double_t , Double_t ) const
{
  // Parametrised "expected" error of the cluster reconstruction in Z 

  Double_t s = 9 * 9 /12.;  
  return s;
}                  

//_____________________________________________________________________
Double_t AliTRDtracker::GetX(Int_t sector, Int_t plane, Int_t localTB) const 
{
  //
  // Returns radial position which corresponds to time bin <localTB>
  // in tracking sector <sector> and plane <plane>
  //

  Int_t index = fTrSec[sector]->CookTimeBinIndex(plane, localTB); 
  Int_t pl = fTrSec[sector]->GetLayerNumber(index);
  return fTrSec[sector]->GetLayer(pl)->GetX();

}


//_______________________________________________________
AliTRDtracker::AliTRDpropagationLayer::AliTRDpropagationLayer(Double_t x, 
               Double_t dx, Double_t rho, Double_t radLength, Int_t tbIndex)
{ 
  //
  // AliTRDpropagationLayer constructor
  //

  fN = 0; fX = x; fdX = dx; fRho = rho; fX0 = radLength;
  fClusters = NULL; fIndex = NULL; fTimeBinIndex = tbIndex;


  for(Int_t i=0; i < (Int_t) kZones; i++) {
    fZc[i]=0; fZmax[i] = 0;
  }

  fYmax = 0;

  if(fTimeBinIndex >= 0) { 
    fClusters = new AliTRDcluster*[kMaxClusterPerTimeBin];
    fIndex = new UInt_t[kMaxClusterPerTimeBin];
  }

  for (Int_t i=0;i<5;i++) fIsHole[i] = kFALSE;
  fHole = kFALSE;
  fHoleZc = 0;
  fHoleZmax = 0;
  fHoleYc = 0;
  fHoleYmax = 0;
  fHoleRho = 0;
  fHoleX0 = 0;

}

//_______________________________________________________
void AliTRDtracker::AliTRDpropagationLayer::SetHole(
          Double_t Zmax, Double_t Ymax, Double_t rho, 
          Double_t radLength, Double_t Yc, Double_t Zc) 
{
  //
  // Sets hole in the layer 
  //
  fHole = kTRUE;
  fHoleZc = Zc;
  fHoleZmax = Zmax;
  fHoleYc = Yc;
  fHoleYmax = Ymax;
  fHoleRho = rho;
  fHoleX0 = radLength;
}
  

//_______________________________________________________
AliTRDtracker::AliTRDtrackingSector::AliTRDtrackingSector(AliTRDgeometry* geo, Int_t gs, AliTRDparameter* par)
{
  //
  // AliTRDtrackingSector Constructor
  //

  fGeom = geo;
  fPar = par;
  fGeomSector = gs;
  fTzeroShift = 0.13;
  fN = 0;
  //
  // get holes description from geometry
  Bool_t holes[AliTRDgeometry::kNcham];
  //printf("sector\t%d\t",gs);
  for (Int_t icham=0; icham<AliTRDgeometry::kNcham;icham++){
    holes[icham] = fGeom->IsHole(0,icham,gs);
    //printf("%d",holes[icham]);
  } 
  //printf("\n");
  
  for(UInt_t i=0; i < kMaxTimeBinIndex; i++) fTimeBinIndex[i] = -1;


  AliTRDpropagationLayer* ppl;

  Double_t x, xin, xout, dx, rho, radLength;
  Int_t    steps;

  // set time bins in the gas of the TPC

  xin = 246.055; xout = 254.055; steps = 20; dx = (xout-xin)/steps;
  rho = 0.9e-3;  radLength = 28.94;

  for(Int_t i=0; i<steps; i++) {
    x = xin + i*dx + dx/2;
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
    InsertLayer(ppl);
  }

  // set time bins in the outer field cage vessel

  dx = 50e-4; xin = xout; xout = xin + dx; rho = 1.71; radLength = 44.77; // Tedlar
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 0.02; xin = xout; xout = xin + dx; rho = 1.45; radLength = 44.86; // prepreg
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 2.; xin = xout; xout = xin + dx; rho = 1.45*0.02; radLength = 41.28; // Nomex
  steps = 5; dx = (xout - xin)/steps;
  for(Int_t i=0; i<steps; i++) {
    x = xin + i*dx + dx/2;
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
    InsertLayer(ppl);
  }

  dx = 0.02; xin = xout; xout = xin + dx; rho = 1.45; radLength = 44.86; // prepreg
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 50e-4; xin = xout; xout = xin + dx; rho = 1.71; radLength = 44.77; // Tedlar
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);


  // set time bins in CO2

  xin = xout; xout = 275.0; 
  steps = 50; dx = (xout - xin)/steps;
  rho = 1.977e-3;  radLength = 36.2;
  
  for(Int_t i=0; i<steps; i++) {
    x = xin + i*dx + dx/2;
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
    InsertLayer(ppl);
  }

  // set time bins in the outer containment vessel

  dx = 50e-4; xin = xout; xout = xin + dx; rho = 2.7; radLength = 24.01; // Al
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 50e-4; xin = xout; xout = xin + dx; rho = 1.71; radLength = 44.77; // Tedlar
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 0.06; xin = xout; xout = xin + dx; rho = 1.45; radLength = 44.86; // prepreg
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 3.; xin = xout; xout = xin + dx; rho = 1.45*0.02; radLength = 41.28; // Nomex
  steps = 10; dx = (xout - xin)/steps;
  for(Int_t i=0; i<steps; i++) {
    x = xin + i*dx + dx/2;
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
    InsertLayer(ppl);
  }

  dx = 0.06; xin = xout; xout = xin + dx; rho = 1.45; radLength = 44.86; // prepreg
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  dx = 50e-4; xin = xout; xout = xin + dx; rho = 1.71; radLength = 44.77; // Tedlar
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);
  
  dx = 50e-4; xin = xout; xout = xin + dx; rho = 2.7; radLength = 24.01; // Al
  ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
  InsertLayer(ppl);

  Double_t xtrd = (Double_t) fGeom->Rmin();  

  // add layers between TPC and TRD (Air temporarily)
  xin = xout; xout = xtrd;
  steps = 50; dx = (xout - xin)/steps;
  rho = 1.2e-3;  radLength = 36.66;
  
  for(Int_t i=0; i<steps; i++) {
    x = xin + i*dx + dx/2;
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
    InsertLayer(ppl);
  }


  //  Double_t alpha=AliTRDgeometry::GetAlpha();

  // add layers for each of the planes

  Double_t dxRo = (Double_t) fGeom->CroHght();    // Rohacell 
  Double_t dxSpace = (Double_t) fGeom->Cspace();  // Spacing between planes
  Double_t dxAmp = (Double_t) fGeom->CamHght();   // Amplification region
  Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region  
  Double_t dxRad = (Double_t) fGeom->CraHght();   // Radiator
  Double_t dxTEC = dxRad + dxDrift + dxAmp + dxRo; 
  Double_t dxPlane = dxTEC + dxSpace; 

  Int_t tb, tbIndex;
  const Int_t  kNchambers = AliTRDgeometry::Ncham();
  Double_t  ymax = 0;
  //, holeYmax = 0;
  Double_t ymaxsensitive=0;
  Double_t *zc = new Double_t[kNchambers];
  Double_t *zmax = new Double_t[kNchambers];
  Double_t *zmaxsensitive = new Double_t[kNchambers];  
  //  Double_t  holeZmax = 1000.;   // the whole sector is missing

  for(Int_t plane = 0; plane < AliTRDgeometry::Nplan(); plane++) {
    //
    // Radiator 
    xin = xtrd + plane * dxPlane; xout = xin + dxRad;
    steps = 12; dx = (xout - xin)/steps; rho = 0.074; radLength = 40.6; 
    for(Int_t i=0; i<steps; i++) {
      x = xin + i*dx + dx/2;
      ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);      
      InsertLayer(ppl);
    }

    ymax          = fGeom->GetChamberWidth(plane)/2.;
    ymaxsensitive = (fPar->GetColPadSize(plane)*fPar->GetColMax(plane)-4)/2.;
    
    for(Int_t ch = 0; ch < kNchambers; ch++) {
      zmax[ch] = fGeom->GetChamberLength(plane,ch)/2;
      Float_t pad = fPar->GetRowPadSize(plane,ch,0);
      Float_t row0 = fPar->GetRow0(plane,ch,0);
      Int_t nPads = fPar->GetRowMax(plane,ch,0);
      zmaxsensitive[ch] = Float_t(nPads)*pad/2.;      
      //      zc[ch] = (pad * nPads)/2 + row0 - pad/2;
      zc[ch] = (pad * nPads)/2 + row0;
      //zc[ch] = row0+zmax[ch]-AliTRDgeometry::RpadW();

    }

    dx = fPar->GetTimeBinSize(); 
    rho = 0.00295 * 0.85; radLength = 11.0;  

    Double_t x0 = (Double_t) fPar->GetTime0(plane);
    Double_t xbottom = x0 - dxDrift;
    Double_t xtop = x0 + dxAmp;
    //
    // Amplification region
    steps = (Int_t) (dxAmp/dx);

    for(tb = 0; tb < steps; tb++) {
      x = x0 + tb * dx + dx/2;
      tbIndex = CookTimeBinIndex(plane, -tb-1);
      ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,tbIndex);
      ppl->SetYmax(ymax,ymaxsensitive);
      ppl->SetZ(zc, zmax, zmaxsensitive);
      ppl->SetHoles(holes);
      InsertLayer(ppl);
    }
    tbIndex = CookTimeBinIndex(plane, -steps);
    x = (x + dx/2 + xtop)/2;
    dx = 2*(xtop-x);
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,tbIndex);
    ppl->SetYmax(ymax,ymaxsensitive);
    ppl->SetZ(zc, zmax,zmaxsensitive);
    ppl->SetHoles(holes);
    InsertLayer(ppl);

    // Drift region
    dx = fPar->GetTimeBinSize();
    steps = (Int_t) (dxDrift/dx);

    for(tb = 0; tb < steps; tb++) {
      x = x0 - tb * dx - dx/2;
      tbIndex = CookTimeBinIndex(plane, tb);

      ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,tbIndex);
      ppl->SetYmax(ymax,ymaxsensitive);
      ppl->SetZ(zc, zmax, zmaxsensitive);
      ppl->SetHoles(holes);
      InsertLayer(ppl);
    }
    tbIndex = CookTimeBinIndex(plane, steps);
    x = (x - dx/2 + xbottom)/2;
    dx = 2*(x-xbottom);
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,tbIndex);
    ppl->SetYmax(ymax,ymaxsensitive);
    ppl->SetZ(zc, zmax, zmaxsensitive);
    ppl->SetHoles(holes);    
    InsertLayer(ppl);

    // Pad Plane
    xin = xtop; dx = 0.025; xout = xin + dx; rho = 1.7; radLength = 33.0;
    ppl = new AliTRDpropagationLayer(xin+dx/2,dx,rho,radLength,-1);
    ppl->SetYmax(ymax,ymaxsensitive);
    ppl->SetZ(zc, zmax,zmax);
    ppl->SetHoles(holes);	  
    InsertLayer(ppl);

    // Rohacell
    xin = xout; xout = xtrd + (plane + 1) * dxPlane - dxSpace;
    steps = 5; dx = (xout - xin)/steps; rho = 0.074; radLength = 40.6; 
    for(Int_t i=0; i<steps; i++) {
      x = xin + i*dx + dx/2;
      ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
      ppl->SetYmax(ymax,ymaxsensitive);
      ppl->SetZ(zc, zmax,zmax);
      ppl->SetHoles(holes);
      InsertLayer(ppl);
    }

    // Space between the chambers, air
    xin = xout; xout = xtrd + (plane + 1) * dxPlane;
    steps = 5; dx = (xout - xin)/steps; rho = 1.29e-3; radLength = 36.66; 
    for(Int_t i=0; i<steps; i++) {
      x = xin + i*dx + dx/2;
      ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
      InsertLayer(ppl);
    }
  }    

  // Space between the TRD and RICH
  Double_t xRICH = 500.;
  xin = xout; xout = xRICH;
  steps = 200; dx = (xout - xin)/steps; rho = 1.29e-3; radLength = 36.66; 
  for(Int_t i=0; i<steps; i++) {
    x = xin + i*dx + dx/2;
    ppl = new AliTRDpropagationLayer(x,dx,rho,radLength,-1);
    InsertLayer(ppl);
  }

  MapTimeBinLayers();
  delete [] zc;
  delete [] zmax;

}

//______________________________________________________

Int_t  AliTRDtracker::AliTRDtrackingSector::CookTimeBinIndex(Int_t plane, Int_t localTB) const
{
  //
  // depending on the digitization parameters calculates "global"
  // time bin index for timebin <localTB> in plane <plane>
  //

  Double_t dxAmp = (Double_t) fGeom->CamHght();   // Amplification region
  Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region  
  Double_t dx = (Double_t) fPar->GetTimeBinSize();  

  Int_t tbAmp = fPar->GetTimeBefore();
  Int_t maxAmp = (Int_t) ((dxAmp+0.000001)/dx);
  if(kTRUE) maxAmp = 0;   // intentional until we change parameter class 
  Int_t tbDrift = fPar->GetTimeMax();
  Int_t maxDrift = (Int_t) ((dxDrift+0.000001)/dx);

  Int_t tbPerPlane = TMath::Min(tbAmp,maxAmp) + TMath::Min(tbDrift,maxDrift);

  Int_t gtb = (plane+1) * tbPerPlane - localTB - 1 - TMath::Min(tbAmp,maxAmp);

  if((localTB < 0) && 
     (TMath::Abs(localTB) > TMath::Min(tbAmp,maxAmp))) return -1;
  if(localTB >= TMath::Min(tbDrift,maxDrift)) return -1;

  return gtb;


}

//______________________________________________________

void AliTRDtracker::AliTRDtrackingSector::MapTimeBinLayers() 
{
  //
  // For all sensitive time bins sets corresponding layer index
  // in the array fTimeBins 
  //

  Int_t index;

  for(Int_t i = 0; i < fN; i++) {
    index = fLayers[i]->GetTimeBinIndex();
    
    //    printf("gtb %d -> pl %d -> x %f \n", index, i, fLayers[i]->GetX());

    if(index < 0) continue;
    if(index >= (Int_t) kMaxTimeBinIndex) {
      printf("*** AliTRDtracker::MapTimeBinLayers: \n");
      printf("    index %d exceeds allowed maximum of %d!\n",
             index, kMaxTimeBinIndex-1);
      continue;
    }
    fTimeBinIndex[index] = i;
  }

  Double_t x1, dx1, x2, dx2, gap;

  for(Int_t i = 0; i < fN-1; i++) {
    x1 = fLayers[i]->GetX();
    dx1 = fLayers[i]->GetdX();
    x2 = fLayers[i+1]->GetX();
    dx2 = fLayers[i+1]->GetdX();
    gap = (x2 - dx2/2) - (x1 + dx1/2);
    if(gap < -0.01) {
      printf("*** warning: layers %d and %d are overlayed:\n",i,i+1);
      printf("             %f + %f + %f > %f\n", x1, dx1/2, dx2/2, x2);
    }
    if(gap > 0.01) { 
      printf("*** warning: layers %d and %d have a large gap:\n",i,i+1);
      printf("             (%f - %f) - (%f + %f) = %f\n", 
             x2, dx2/2, x1, dx1, gap);
    }
  }
}
  

//______________________________________________________


Int_t AliTRDtracker::AliTRDtrackingSector::GetLayerNumber(Double_t x) const
{
  // 
  // Returns the number of time bin which in radial position is closest to <x>
  //

  if(x >= fLayers[fN-1]->GetX()) return fN-1; 
  if(x <= fLayers[0]->GetX()) return 0; 

  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (x > fLayers[m]->GetX()) b=m+1;
    else e=m;
  }
  if(TMath::Abs(x - fLayers[m]->GetX()) > 
     TMath::Abs(x - fLayers[m+1]->GetX())) return m+1;
  else return m;

}

//______________________________________________________

Int_t AliTRDtracker::AliTRDtrackingSector::GetInnerTimeBin() const 
{
  // 
  // Returns number of the innermost SENSITIVE propagation layer
  //

  return GetLayerNumber(0);
}

//______________________________________________________

Int_t AliTRDtracker::AliTRDtrackingSector::GetOuterTimeBin() const 
{
  // 
  // Returns number of the outermost SENSITIVE time bin
  //

  return GetLayerNumber(GetNumberOfTimeBins() - 1);
}

//______________________________________________________

Int_t AliTRDtracker::AliTRDtrackingSector::GetNumberOfTimeBins() const 
{
  // 
  // Returns number of SENSITIVE time bins
  //

  Int_t tb, layer;
  for(tb = kMaxTimeBinIndex-1; tb >=0; tb--) {
    layer = GetLayerNumber(tb);
    if(layer>=0) break;
  }
  return tb+1;
}

//______________________________________________________

void AliTRDtracker::AliTRDtrackingSector::InsertLayer(AliTRDpropagationLayer* pl)
{ 
  //
  // Insert layer <pl> in fLayers array.
  // Layers are sorted according to X coordinate.

  if ( fN == ((Int_t) kMaxLayersPerSector)) {
    printf("AliTRDtrackingSector::InsertLayer(): Too many layers !\n");
    return;
  }
  if (fN==0) {fLayers[fN++] = pl; return;}
  Int_t i=Find(pl->GetX());

  memmove(fLayers+i+1 ,fLayers+i,(fN-i)*sizeof(AliTRDpropagationLayer*));
  fLayers[i]=pl; fN++;

}              

//______________________________________________________

Int_t AliTRDtracker::AliTRDtrackingSector::Find(Double_t x) const 
{
  //
  // Returns index of the propagation layer nearest to X 
  //

  if (x <= fLayers[0]->GetX()) return 0;
  if (x > fLayers[fN-1]->GetX()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (x > fLayers[m]->GetX()) b=m+1;
    else e=m;
  }
  return m;
}             

//______________________________________________________
void AliTRDtracker::AliTRDpropagationLayer::SetZ(Double_t* center, Double_t *w, Double_t *wsensitive )
{
  //
  // set centers and the width of sectors
  for (Int_t icham=0;icham< AliTRDgeometry::kNcham;icham++){
    fZc[icham] = center[icham];  
    fZmax[icham] = w[icham];
    fZmaxSensitive[icham] = wsensitive[icham];
    //   printf("chamber\t%d\tzc\t%f\tzmax\t%f\tzsens\t%f\n",icham,fZc[icham],fZmax[icham],fZmaxSensitive[icham]);
  }  
}
//______________________________________________________

void AliTRDtracker::AliTRDpropagationLayer::SetHoles(Bool_t *holes)
{
  //
  // set centers and the width of sectors
  fHole = kFALSE;
  for (Int_t icham=0;icham< AliTRDgeometry::kNcham;icham++){
    fIsHole[icham] = holes[icham]; 
    if (holes[icham]) fHole = kTRUE;
  }  
}



void AliTRDtracker::AliTRDpropagationLayer::GetPropagationParameters(
        Double_t y, Double_t z, Double_t &dx, Double_t &rho, Double_t &radLength, 
        Bool_t &lookForCluster) const
{
  //
  // Returns radial step <dx>, density <rho>, rad. length <radLength>,
  // and sensitivity <lookForCluster> in point <y,z>  
  //

  dx  = fdX;
  rho = fRho;
  radLength  = fX0;
  lookForCluster = kFALSE;
  //
  // check dead regions in sensitive volume 
  if(fTimeBinIndex >= 0) {
    //
    Int_t zone=-1;
    for(Int_t ch = 0; ch < (Int_t) kZones; ch++) {
      if  (TMath::Abs(z - fZc[ch]) < fZmaxSensitive[ch]){ 
	zone = ch;
	lookForCluster = !(fIsHole[zone]);
	if(TMath::Abs(y) > fYmaxSensitive){  
	  lookForCluster = kFALSE;
	}
	if (fIsHole[zone]) {
	  //if hole
	  rho = 1.29e-3;
	  radLength = 36.66;
	}
      }    
    }
    return;
  }
  //
  //
  // check hole
  if (fHole==kFALSE) return;
  //
  for(Int_t ch = 0; ch < (Int_t) kZones; ch++) {
    if  (TMath::Abs(z - fZc[ch]) < fZmax[ch]){ 
      if (fIsHole[ch]) {
	//if hole
	rho = 1.29e-3;
	radLength = 36.66;
      }
    }
  }
  return;
}

Int_t  AliTRDtracker::AliTRDpropagationLayer::GetZone( Double_t z) const
{
  //
  //
  if (fTimeBinIndex < 0) return -20;  //unknown 
  Int_t zone=-10;   // dead zone
  for(Int_t ch = 0; ch < (Int_t) kZones; ch++) {
    if(TMath::Abs(z - fZc[ch]) < fZmax[ch]) 
      zone = ch;
  }
  return zone;
}


//______________________________________________________

void AliTRDtracker::AliTRDpropagationLayer::InsertCluster(AliTRDcluster* c, 
                                                          UInt_t index) {

// Insert cluster in cluster array.
// Clusters are sorted according to Y coordinate.  

  if(fTimeBinIndex < 0) { 
    printf("*** attempt to insert cluster into non-sensitive time bin!\n");
    return;
  }

  if (fN== (Int_t) kMaxClusterPerTimeBin) {
    printf("AliTRDpropagationLayer::InsertCluster(): Too many clusters !\n"); 
    return;
  }
  if (fN==0) {fIndex[0]=index; fClusters[fN++]=c; return;}
  Int_t i=Find(c->GetY());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTRDcluster*));
  memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t)); 
  fIndex[i]=index; fClusters[i]=c; fN++;
}  

//______________________________________________________

Int_t AliTRDtracker::AliTRDpropagationLayer::Find(Double_t y) const {

// Returns index of the cluster nearest in Y    

  if (y <= fClusters[0]->GetY()) return 0;
  if (y > fClusters[fN-1]->GetY()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (y > fClusters[m]->GetY()) b=m+1;
    else e=m;
  }
  return m;
}    

//---------------------------------------------------------

Double_t AliTRDtracker::GetTiltFactor(const AliTRDcluster* c) {
//
//  Returns correction factor for tilted pads geometry 
//

  Double_t h01 = sin(TMath::Pi() / 180.0 * fPar->GetTiltingAngle());
  Int_t det = c->GetDetector();    
  Int_t plane = fGeom->GetPlane(det);

  //if((plane == 1) || (plane == 3) || (plane == 5)) h01=-h01;
  if((plane == 0) || (plane == 2) || (plane == 4)) h01=-h01;

  if(fNoTilt) h01 = 0;
  
  return h01;
}


void AliTRDtracker::CookdEdxTimBin(AliTRDtrack& TRDtrack)
{
  // *** ADDED TO GET MORE INFORMATION FOR TRD PID  ---- PS
  // This is setting fdEdxPlane and fTimBinPlane
  // Sums up the charge in each plane for track TRDtrack and also get the 
  // Time bin for Max. Cluster
  // Prashant Shukla (shukla@physi.uni-heidelberg.de)

  //  const Int_t kNPlane = AliTRDgeometry::Nplan();
  //  const Int_t kNPlane = 6;
  Double_t  clscharge[kNPlane], maxclscharge[kNPlane];
  Int_t  nCluster[kNPlane], timebin[kNPlane];

  //Initialization of cluster charge per plane.  
  for (Int_t iPlane = 0; iPlane < kNPlane; iPlane++) {
    clscharge[iPlane] = 0.0;
    nCluster[iPlane] = 0;
    timebin[iPlane] = -1;
    maxclscharge[iPlane] = 0.0;
  }

  // Loop through all clusters associated to track TRDtrack
  Int_t nClus = TRDtrack.GetNumberOfClusters();  // from Kalmantrack
  for (Int_t iClus = 0; iClus < nClus; iClus++) {
    Double_t charge = TRDtrack.GetClusterdQdl(iClus);
    Int_t index = TRDtrack.GetClusterIndex(iClus);
    AliTRDcluster *TRDcluster = (AliTRDcluster *) GetCluster(index); 
    if (!TRDcluster) continue;
    Int_t tb = TRDcluster->GetLocalTimeBin();
    if (!tb) continue;
    Int_t detector = TRDcluster->GetDetector();
    Int_t iPlane   = fGeom->GetPlane(detector);
    clscharge[iPlane] = clscharge[iPlane]+charge;
    if(charge > maxclscharge[iPlane]) {
      maxclscharge[iPlane] = charge;
      timebin[iPlane] = tb;
    }
    nCluster[iPlane]++;
  } // end of loop over cluster

  // Setting the fdEdxPlane and fTimBinPlane variabales 
  Double_t Total_ch = 0;
  for (Int_t iPlane = 0; iPlane < kNPlane; iPlane++) {
    if (nCluster[iPlane]) clscharge[iPlane] /= nCluster[iPlane];
    TRDtrack.SetPIDsignals(clscharge[iPlane], iPlane);
    TRDtrack.SetPIDTimBin(timebin[iPlane], iPlane);
    Total_ch= Total_ch+clscharge[iPlane];
  }
  //  Int_t i;
  //  Int_t nc=TRDtrack.GetNumberOfClusters(); 
  //  Float_t dedx=0;
  //  for (i=0; i<nc; i++) dedx += TRDtrack.GetClusterdQdl(i);
  //  dedx /= nc;
  //  for (Int_t iPlane = 0; iPlane < kNPlane; iPlane++) {
  //    TRDtrack.SetPIDsignals(dedx, iPlane);
  //    TRDtrack.SetPIDTimBin(timbin[iPlane], iPlane);
  //  }

} // end of function





