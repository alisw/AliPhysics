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
#include "AliTRDpadPlane.h"
#include "AliTRDgeometryDetail.h"
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliESD.h"

#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"

#include "TTreeStream.h"
#include "TGraph.h"
#include "AliTRDtracker.h"
#include "TLinearFitter.h"
#include "AliRieman.h"
#include "AliTrackPointArray.h"
#include "AliAlignObj.h"

//

ClassImp(AliTRDtracker) 
ClassImp(AliTRDseed)

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

//   const  Double_t    AliTRDtracker::fgkOffset             = -0.012;
//   const  Double_t    AliTRDtracker::fgkOffsetX            = 0.35;
//   const  Double_t    AliTRDtracker::fgkCoef               = 0.00;
//   const  Double_t    AliTRDtracker::fgkMean               = 8.;
//   const  Double_t    AliTRDtracker::fgkDriftCorrection    = 1.07;
//   const  Double_t    AliTRDtracker::fgkExB                = 0.072;

  const  Double_t    AliTRDtracker::fgkOffset             = -0.019;
  const  Double_t    AliTRDtracker::fgkOffsetX            = 0.26;       // "time offset"  
//  const  Double_t    AliTRDtracker::fgkCoef               = 0.0096;   // angular shift 
  const  Double_t    AliTRDtracker::fgkCoef               = 0.0106;   // angular shift 
  const  Double_t    AliTRDtracker::fgkMean               = 0.;
  const  Double_t    AliTRDtracker::fgkDriftCorrection    = 1.055;   // drift coefficient correction
  const  Double_t    AliTRDtracker::fgkExB                = 0.072;  // ExB angle - for error parameterization


//     poscorrection =  fgkCoef*(GetLocalTimeBin() - fgkMean)+fgkOffset; 

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
  fDebugStreamer = 0;
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
    //    printf("Found geometry version %d on file \n", fGeom->IsVersion());
  }
  else { 
    printf("AliTRDtracker::AliTRDtracker(): can't find TRD geometry!\n");
    //printf("The DETAIL TRD geometry will be used\n");
    //fGeom = new AliTRDgeometryDetail();
    fGeom = new AliTRDgeometryDetail();
    fGeom->SetPHOShole();
    fGeom->SetRICHhole();    
  } 

  if (!fPar) {  
    printf("AliTRDtracker::AliTRDtracker(): can't find TRD parameter!\n");
    printf("The DEFAULT TRD parameter will be used\n");
    fPar = new AliTRDparameter("Pica","Vyjebana");
  }
  fPar = new AliTRDparameter("Pica","Vyjebana");
  fPar->Init();

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
  AliTRDpadPlane *padPlane = AliTRDCommonParam::Instance()->GetPadPlane(0,0);
  Float_t tiltAngle = TMath::Abs(padPlane->GetTiltingAngle());
  //  Float_t tiltAngle = TMath::Abs(fPar->GetTiltingAngle()); 
  if(tiltAngle < 0.1) {
    fNoTilt = kTRUE;
  }

  fSY2corr = 0.2;
  fSZ2corr = 120.;      

  if(fNoTilt && (tiltAngle > 0.1)) fSY2corr = fSY2corr + tiltAngle * 0.05; 


  // calculate max gap on track

  Double_t dxAmp = (Double_t) fGeom->CamHght();   // Amplification region
  Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region

  Double_t dx = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
                         / AliTRDcalibDB::Instance()->GetSamplingFrequency();

  Int_t tbAmp = fPar->GetTimeBefore();
  Int_t maxAmp = (Int_t) ((dxAmp+0.000001)/dx);
  if(kTRUE) maxAmp = 0;  // intentional until we change the parameter class 
  Int_t tbDrift = fPar->GetTimeMax();
  Int_t maxDrift = (Int_t) ((dxDrift+0.000001)/dx)+4; // MI change  - take also last time bins

  tbDrift = TMath::Min(tbDrift,maxDrift);  
  tbAmp = TMath::Min(tbAmp,maxAmp);          

  fTimeBinsPerPlane = tbAmp + tbDrift;
  fMaxGap = (Int_t) (fTimeBinsPerPlane * fGeom->Nplan() * fgkSkipDepth);

  fVocal = kFALSE;
  
  fDebugStreamer = new TTreeSRedirector("TRDdebug.root");

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
  if (fDebugStreamer) {    
    //fDebugStreamer->Close();
    delete fDebugStreamer;
  }
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


AliTRDcluster * AliTRDtracker::GetCluster(AliTRDtrack * track, Int_t plane, Int_t timebin, UInt_t &index){
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
      index = indexes[i];
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
  //
  //Sort tracks
  Float_t *quality =new Float_t[n];
  Int_t *index   =new Int_t[n];
  for (Int_t i=0; i<n; i++) {
    AliESDtrack* seed=event->GetTrack(i);
    Double_t covariance[15];
    seed->GetExternalCovariance(covariance);
    quality[i] = covariance[0]+covariance[2];      
  }
  TMath::Sort(n,quality,index,kFALSE);
  //
  for (Int_t i=0; i<n; i++) {
    //    AliESDtrack* seed=event->GetTrack(i);
    AliESDtrack* seed=event->GetTrack(index[i]);

    ULong_t status=seed->GetStatus();
    if ( (status & AliESDtrack::kTPCout ) == 0 ) continue;
    if ( (status & AliESDtrack::kTRDout) != 0 ) continue;

    Int_t lbl = seed->GetLabel();
    AliTRDtrack *track = new AliTRDtrack(*seed);
    track->SetSeedLabel(lbl);
    seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup); //make backup
    fNseeds++;
    Float_t p4     = track->GetC();
    //
    Int_t expectedClr = FollowBackProlongationG(*track);
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
    if (TMath::Abs(track->GetC()-p4)/TMath::Abs(p4)<0.2 || TMath::Abs(track->GetPt())>0.8 ) {
      // 
      //make backup for back propagation 
      //
      Int_t foundClr = track->GetNumberOfClusters();
      if (foundClr >= foundMin) {
	track->CookdEdx(); 
	CookdEdxTimBin(*track);
	CookLabel(track, 1-fgkLabelFraction);
	if (track->GetBackupTrack()) UseClusters(track->GetBackupTrack());
	if(track->GetChi2()/track->GetNumberOfClusters()<4) {   // sign only gold tracks
	  if (seed->GetKinkIndex(0)==0&&TMath::Abs(track->GetPt())<1.5 ) UseClusters(track);
	}
	Bool_t isGold = kFALSE;
	
	if (track->GetChi2()/track->GetNumberOfClusters()<5) {  //full gold track
	  // seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	   if (track->GetBackupTrack()) seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
	  isGold = kTRUE;
	}
	if (!isGold && track->GetNCross()==0&&track->GetChi2()/track->GetNumberOfClusters()<7){ //almost gold track
	  //	  seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	  if (track->GetBackupTrack()) seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
	  isGold = kTRUE;
	}
	if (!isGold && track->GetBackupTrack()){
	  if (track->GetBackupTrack()->GetNumberOfClusters()>foundMin&&
	      (track->GetBackupTrack()->GetChi2()/(track->GetBackupTrack()->GetNumberOfClusters()+1))<7){	  
	    seed->UpdateTrackParams(track->GetBackupTrack(), AliESDtrack::kTRDbackup);
	    isGold = kTRUE;
	  }
	}
	if (track->StatusForTOF()>0 &&track->fNCross==0 && Float_t(track->fN)/Float_t(track->fNExpected)>0.4){
	  seed->UpdateTrackParams(track, AliESDtrack::kTRDbackup);
	}
      }
    }
    // Debug part of tracking
    TTreeSRedirector& cstream = *fDebugStreamer;
    Int_t eventNr = event->GetEventNumber();
    if (track->GetBackupTrack()){
      cstream<<"Tracks"<<
	"EventNr="<<eventNr<<
	"ESD.="<<seed<<
	"trd.="<<track<<
	"trdback.="<<track->GetBackupTrack()<<	
	"\n";
    }else{
      cstream<<"Tracks"<<
	"EventNr="<<eventNr<<
	"ESD.="<<seed<<
	"trd.="<<track<<
	"trdback.="<<track<<
	"\n";
    }
    //
    //Propagation to the TOF (I.Belikov)    
    if (track->GetStop()==kFALSE){
      
      Double_t xtof=371.;
      Double_t c2=track->GetC()*xtof - track->GetEta();
      if (TMath::Abs(c2)>=0.99) {
	delete track;
	continue;
      }
      Double_t xTOF0 = 365. ;          
      PropagateToOuterPlane(*track,xTOF0); 
      //
      //energy losses taken to the account - check one more time
      c2=track->GetC()*xtof - track->GetEta();
      if (TMath::Abs(c2)>=0.99) {
	delete track;
	continue;
      }

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
	//	seed->SetTRDtrack(new AliTRDtrack(*track));
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
	//seed->SetTRDtrack(new AliTRDtrack(*track));
	found++;
      }
    }
    seed->SetTRDQuality(track->StatusForTOF());    
    seed->SetTRDBudget(track->fBudget[0]);    
  
    delete track;
    //
    //End of propagation to the TOF
    //if (foundClr>foundMin)
    //  seed->UpdateTrackParams(track, AliESDtrack::kTRDout);
    

  }
  
  cerr<<"Number of seeds: "<<fNseeds<<endl;  
  cerr<<"Number of back propagated TRD tracks: "<<found<<endl;
  
  //  MakeSeedsMI(3,5,event); //new seeding


  fSeeds->Clear(); fNseeds=0;
  delete [] index;
  delete [] quality;
  
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
  AliTRDtrack seed2;

  Int_t n = event->GetNumberOfTracks();
  for (Int_t i=0; i<n; i++) {
    AliESDtrack* seed=event->GetTrack(i);
    new(&seed2) AliTRDtrack(*seed);
    if (seed2.GetX()<270){
      seed->UpdateTrackParams(&seed2, AliESDtrack::kTRDbackup); // backup TPC track - only update
      continue;
    }

    ULong_t status=seed->GetStatus();
    if ( (status & AliESDtrack::kTRDout ) == 0 ) {
      continue;
    }
    if ( (status & AliESDtrack::kTRDin) != 0 ) {
      continue;
    }
    nseed++;    
//     if (1/seed2.Get1Pt()>1.5&& seed2.GetX()>260.) {
//       Double_t oldx = seed2.GetX();
//       seed2.PropagateTo(500.);
//       seed2.ResetCovariance(1.);
//       seed2.PropagateTo(oldx);
//     }
//     else{
//       seed2.ResetCovariance(5.); 
//     }

    AliTRDtrack *pt = new AliTRDtrack(seed2,seed2.GetAlpha());
    UInt_t * indexes2 = seed2.GetIndexes();
    for (Int_t i=0;i<kNPlane;i++) {
      pt->SetPIDsignals(seed2.GetPIDsignals(i),i);
      pt->SetPIDTimBin(seed2.GetPIDTimBin(i),i);
    }

    UInt_t * indexes3 = pt->GetBackupIndexes();
    for (Int_t i=0;i<200;i++) {
      if (indexes2[i]==0) break;
      indexes3[i] = indexes2[i];
    }          
    //AliTRDtrack *pt = seed2;
    AliTRDtrack &t=*pt; 
    FollowProlongationG(t, innerTB); 
    if (t.GetNumberOfClusters() >= foundMin) {
      //      UseClusters(&t);
      //CookLabel(pt, 1-fgkLabelFraction);
      t.CookdEdx();
      CookdEdxTimBin(t);
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
        //pt2->CookdEdx(0.,1.);
        pt2->CookdEdx( ); // Modification by PS
        CookdEdxTimBin(*pt2);
	seed->UpdateTrackParams(pt2, AliESDtrack::kTRDrefit);
        for (Int_t i=0;i<kNPlane;i++) {
          seed->SetTRDsignals(pt2->GetPIDsignals(i),i);
          seed->SetTRDTimBin(pt2->GetPIDTimBin(i),i);
        }
      }
      delete pt2;
    }  
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
  Double_t dxsample = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
                         / AliTRDcalibDB::Instance()->GetSamplingFrequency();
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
      if (timeBin&&fVocal) {
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
	AliTRDcluster * cl2= GetCluster(&t,plane, timebin,index);
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

        if (cl) {
	  
          wYclosest = cl->GetY();
          wZclosest = cl->GetZ();
          Double_t h01 = GetTiltFactor(cl);

          //if (cl->GetNPads()<5) 
	    t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 
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
      }
    }  
  }
  return expectedNumberOfClusters;
  
  
}                




//---------------------------------------------------------------------------
Int_t AliTRDtracker::FollowProlongationG(AliTRDtrack& t, Int_t rf)
{
  // Starting from current position on track=t this function tries
  // to extrapolate the track up to timeBin=0 and to confirm prolongation
  // if a close cluster is found. Returns the number of clusters
  // expected to be found in sensitive layers
  // GeoManager used to estimate mean density
  Int_t sector;
  Int_t lastplane = GetLastPlane(&t);
  Int_t tryAgain=fMaxGap;
  Double_t alpha=t.GetAlpha();
  alpha = TVector2::Phi_0_2pi(alpha);
  Double_t radLength = 0.0;
  Double_t rho = 0.0;
  Double_t x, dx;
  //, y, ymax, z;
  Int_t expectedNumberOfClusters = 0;
  Double_t dxsample = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
                         / AliTRDcalibDB::Instance()->GetSamplingFrequency();
  //
  //
  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning
  Double_t tanmax = TMath::Tan(0.5*alpha);  
 
  for (Int_t nr=fTrSec[0]->GetLayerNumber(t.GetX()); nr>rf; nr--) { 
    //
    // propagate track in  non active layers
    //
    if (!(fTrSec[0]->GetLayer(nr)->IsSensitive())){      
      Double_t xyz0[3],xyz1[3],param[7],x,y,z;
      t.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);   //starting global position
      while (nr >rf && (!(fTrSec[0]->GetLayer(nr)->IsSensitive()))){
	x = fTrSec[0]->GetLayer(nr)->GetX();
	nr--;
	if (!t.GetProlongation(x,y,z)) break;
	if (TMath::Abs(y)>x*tanmax){
	  nr--;
	  break;	  
	}
      }
      nr++;
      x = fTrSec[0]->GetLayer(nr)->GetX();
      if (!t.GetProlongation(x,y,z)) break;      
      //
      // minimal mean and maximal budget scan
      Float_t minbudget  =10000;
      Float_t meanbudget =0;
      Float_t maxbudget  =-1;
//      Float_t normbudget =0;
//       for (Int_t idy=-1;idy<=1;idy++)
// 	for (Int_t idz=-1;idz<=1;idz++){
      for (Int_t idy=0;idy<1;idy++)
	for (Int_t idz=0;idz<1;idz++){
	  Double_t y2 = y+idy*TMath::Min(TMath::Sqrt(t.GetSigmaY2()),1.);
	  Double_t z2 = z+idz*TMath::Min(TMath::Sqrt(t.GetSigmaY2()),1.);
	  xyz1[0] = x*TMath::Cos(t.GetAlpha())-y2*TMath::Sin(t.GetAlpha()); 
	  xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y2*TMath::Cos(t.GetAlpha());
	  xyz1[2] = z2;
	  AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
	  Float_t budget = param[0]*param[4];
	  meanbudget+=budget;
	  if (budget<minbudget) minbudget=budget;
	  if (budget>maxbudget) maxbudget=budget;
	}
      t.fBudget[0]+=minbudget;
      t.fBudget[1]+=meanbudget/9.;
      t.fBudget[2]+=minbudget;
      //
      //
      xyz1[0] = x*TMath::Cos(t.GetAlpha())-y*TMath::Sin(t.GetAlpha()); 
      xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y*TMath::Cos(t.GetAlpha());
      xyz1[2] = z;
      AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	      

      t.PropagateTo(x,param[1],param[0]);	
      t.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);   //end  global position	
      AdjustSector(&t);
      continue;
    }
     //
    //
    // stop tracking for highly inclined tracks 
    if (!AdjustSector(&t)) break;
    if (TMath::Abs(t.GetSnp())>0.95) break;
    //
    // propagate and update track in active layers
    //
    Int_t nr0 = nr;  //first active layer
    if (nr >rf  && (fTrSec[0]->GetLayer(nr)->IsSensitive())){      
      Double_t xyz0[3],xyz1[3],param[7],x,y,z;
      t.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);   //starting global position
      while (nr >rf && ((fTrSec[0]->GetLayer(nr)->IsSensitive()))){
	x = fTrSec[0]->GetLayer(nr)->GetX();
	nr--;
	if (!t.GetProlongation(x,y,z)) break;
	if (TMath::Abs(y)>x*tanmax){
	  nr--;
	  break;	  
	}
      }
      //      nr++;
      x = fTrSec[0]->GetLayer(nr)->GetX();
      if (!t.GetProlongation(x,y,z)) break;
      xyz1[0] = x*TMath::Cos(t.GetAlpha())-y*TMath::Sin(t.GetAlpha()); 
      xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y*TMath::Cos(t.GetAlpha());
      xyz1[2] = z;
      // end global position
      AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
      rho = param[0];
      radLength = param[1];   // get mean propagation parameters
    }
    //
    // propagate and update
    if (nr0-nr<5){
      // short tracklet - do not update - edge effect
      x = fTrSec[0]->GetLayer(nr)->GetX();
      t.PropagateTo(x,radLength,rho);
      AdjustSector(&t);
      continue; 
    }
    sector = t.GetSector();
    //
    //    
    for (Int_t ilayer=nr0;ilayer>=nr;ilayer--) {
      expectedNumberOfClusters++;       
      t.fNExpected++;
      if (t.fX>345) t.fNExpectedLast++;
      AliTRDpropagationLayer& timeBin=*(fTrSec[sector]->GetLayer(ilayer));
      AliTRDcluster *cl=0;
      UInt_t index=0;
      Double_t maxChi2=fgkMaxChi2;
      dx = (fTrSec[sector]->GetLayer(ilayer+1))->GetX()-timeBin.GetX();
      x = timeBin.GetX();
      t.PropagateTo(x,radLength,rho);
      // Now go for the real cluster search
      if (timeBin) {
	AliTRDcluster * cl0 = timeBin[0];
	if (!cl0) continue;         // no clusters in given time bin
	Int_t plane = fGeom->GetPlane(cl0->GetDetector());
	if (plane>lastplane) continue;
	Int_t timebin = cl0->GetLocalTimeBin();
	AliTRDcluster * cl2= GetCluster(&t,plane, timebin,index);
	//
	if (cl2) {
	  cl =cl2;	
	  Double_t h01 = GetTiltFactor(cl);
	  maxChi2=t.GetPredictedChi2(cl,h01);
	}
	
        if (cl) {
	  //	  if (cl->GetNPads()<5) 
	    t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 
          Double_t h01 = GetTiltFactor(cl);
	  Int_t det = cl->GetDetector();    
	  Int_t plane = fGeom->GetPlane(det);
	  if (t.fX>345){
	    t.fNLast++;
	    t.fChi2Last+=maxChi2;
	  }
	  if(!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
	    if(!t.Update(cl,maxChi2,index,h01)) {
	      //if(!tryAgain--) return 0;
	    }
          }  
          else tryAgain=fMaxGap;
	  //	  	 	  
	}			
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

  Int_t tryAgain=fMaxGap;

  Double_t alpha=t.GetAlpha();
  TVector2::Phi_0_2pi(alpha);

  Int_t s;
  
  Int_t clusters[1000];
  for (Int_t i=0;i<1000;i++) clusters[i]=-1;

  Int_t outerTB = fTrSec[0]->GetOuterTimeBin();
  //Double_t radLength, rho, x, dx, y, ymax = 0, z;
  Double_t radLength, rho, x, dx, y, z;
  Bool_t lookForCluster;
  Double_t dxsample = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
                         / AliTRDcalibDB::Instance()->GetSamplingFrequency();

  Int_t expectedNumberOfClusters = 0;
  x = t.GetX();

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning
  
  //  Int_t zone =0;
  Int_t nr;
  Float_t ratio0=0;
  AliTRDtracklet tracklet;
  //
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
    //
    // MI -fix untill correct material desription will be implemented
    //
    //Int_t nrotate = t.GetNRotate();
    if (!AdjustSector(&t)) break;    
    //
    //
    y = t.GetY();
    z = t.GetZ();
    s = t.GetSector();

    // now propagate to the middle plane of the next time bin 
    fTrSec[s]->GetLayer(nr+1)->GetPropagationParameters(y,z,dx,rho,radLength,lookForCluster);
//     if (nrotate!=t.GetNRotate()){
//       rho = 1000*2.7; radLength = 24.01;  //TEMPORARY - aluminium in between z - will be detected using GeoModeler in future versions
//     }
    x = fTrSec[s]->GetLayer(nr+1)->GetX(); 
    if(!t.PropagateTo(x,radLength,rho)) break;
    if (!AdjustSector(&t)) break;
    s = t.GetSector();
    //    if(!t.PropagateTo(x,radLength,rho)) break;
    
    if (TMath::Abs(t.GetSnp())>0.95) break;

    y = t.GetY();
    z = t.GetZ();

    if(lookForCluster) {
      if (clusters[nr]==-1) {
	Float_t  ncl   = FindClusters(s,nr,nr+30,&t,clusters,tracklet);
	ratio0 = ncl/Float_t(fTimeBinsPerPlane);
	Float_t  ratio1 = Float_t(t.fN+1)/Float_t(t.fNExpected+1.);	
	if (tracklet.GetChi2()<18.&&ratio0>0.8 && ratio1>0.6 && ratio0+ratio1>1.5 && t.GetNCross()==0 && TMath::Abs(t.GetSnp())<0.85&&t.fN>20){
	  t.MakeBackupTrack();                            // make backup of the track until is gold
	}
// 	if (ncl>4){
//  	  t.PropagateTo(tracklet.GetX());
//  	  t.UpdateMI(tracklet);
//  	  nr = fTrSec[0]->GetLayerNumber(t.GetX())+1;
//  	  continue;
// 	}
      }

      expectedNumberOfClusters++;       
      t.fNExpected++;
      if (t.fX>345) t.fNExpectedLast++;

      AliTRDpropagationLayer& timeBin=*(fTrSec[s]->GetLayer(nr+1));
      Double_t sy2=ExpectedSigmaY2(t.GetX(),t.GetTgl(),t.GetPt());
      if((t.GetSigmaY2() + sy2) < 0) {
	printf("problem\n");
	break;
      }
      Double_t road = 10.*sqrt(t.GetSigmaY2() + sy2); 
      
      if (road>fgkWideRoad) {
        return 0;
      }      

      AliTRDcluster *cl=0;
      UInt_t index=0;
      Double_t maxChi2=fgkMaxChi2;
      
      // Now go for the real cluster search

      if (timeBin) {
	
	if (clusters[nr+1]>0) {
	  index = clusters[nr+1];
	  cl    = (AliTRDcluster*)GetCluster(index);
	  Double_t h01 = GetTiltFactor(cl);
          maxChi2=t.GetPredictedChi2(cl,h01);          
	}
	
        if (cl) {
	  //	  if (cl->GetNPads()<5) 
	    t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 
          Double_t h01 = GetTiltFactor(cl);
	  Int_t det = cl->GetDetector();    
	  Int_t plane = fGeom->GetPlane(det);
	  if (t.fX>345){
	    t.fNLast++;
	    t.fChi2Last+=maxChi2;
	  }
	  if(!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
	    if(!t.Update(cl,maxChi2,index,h01)) {
	      //if(!tryAgain--) return 0;
	    }
          }  
          else tryAgain=fMaxGap;
	  //
	  
	  if (cl->GetLocalTimeBin()==1&&t.fN>20 && float(t.fChi2)/float(t.fN)<5){
	    Float_t  ratio1 = Float_t(t.fN)/Float_t(t.fNExpected);	
	    if (tracklet.GetChi2()<18&&ratio0>0.8&&ratio1>0.6 &&ratio0+ratio1>1.5 && t.GetNCross()==0 && TMath::Abs(t.GetSnp())<0.85){
	      t.MakeBackupTrack();                            // make backup of the track until is gold
	    }
	  }
	  
	}			
        else {
	  // if (tryAgain==0) break; 
          //tryAgain--;                                                                               
        }		
      }
    }  
  }
  if (nr<outerTB) 
    t.SetStop(kTRUE);
  else
    t.SetStop(kFALSE);
  return expectedNumberOfClusters;  
}         



//___________________________________________________________________
Int_t AliTRDtracker::FollowBackProlongationG(AliTRDtrack& t)
{
  
  // Starting from current radial position of track <t> this function
  // extrapolates the track up to outer timebin and in the sensitive
  // layers confirms prolongation if a close cluster is found. 
  // Returns the number of clusters expected to be found in sensitive layers
  // Use GEO manager for material Description
  Int_t tryAgain=fMaxGap;

  Double_t alpha=t.GetAlpha();
  TVector2::Phi_0_2pi(alpha);
  Int_t sector;
  Int_t clusters[1000];
  for (Int_t i=0;i<1000;i++) clusters[i]=-1;
  Double_t dxsample = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
    / AliTRDcalibDB::Instance()->GetSamplingFrequency();
  Int_t outerTB = fTrSec[0]->GetOuterTimeBin();
  Double_t radLength = 0.0;
  Double_t rho = 0.0;
  Double_t x, dx; //y, z;

  Int_t expectedNumberOfClusters = 0;
  x = t.GetX();

  alpha=AliTRDgeometry::GetAlpha();  // note: change in meaning
  Double_t tanmax = TMath::Tan(0.5*alpha);  
  Int_t nr;
  Float_t ratio0=0;
  AliTRDtracklet tracklet;
  //
  //
  //
  for (nr=fTrSec[0]->GetLayerNumber(t.GetX()); nr<outerTB+1; nr++) { 
    //
    // propagate track in  non active layers
    //
    if (!(fTrSec[0]->GetLayer(nr)->IsSensitive())){      
      Double_t xyz0[3],xyz1[3],param[7],x,y,z;
      t.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);   //starting global position
      while (nr <outerTB && (!(fTrSec[0]->GetLayer(nr)->IsSensitive()))){
	x = fTrSec[0]->GetLayer(nr)->GetX();
	nr++;
	if (!t.GetProlongation(x,y,z)) break;
	if (TMath::Abs(y)>x*tanmax){
	  nr++;
	  break;	  
	}
      }
      nr--;
      x = fTrSec[0]->GetLayer(nr)->GetX();
      if (!t.GetProlongation(x,y,z)) break;
      // minimal mean and maximal budget scan
      Float_t minbudget  =10000;
      Float_t meanbudget =0;
      Float_t maxbudget  =-1;
//      Float_t normbudget =0;
//       for (Int_t idy=-1;idy<=1;idy++)
// 	for (Int_t idz=-1;idz<=1;idz++){
      for (Int_t idy=0;idy<1;idy++)
	for (Int_t idz=0;idz<1;idz++){
	  Double_t y2 = y+idy*TMath::Min(TMath::Sqrt(t.GetSigmaY2()),1.);
	  Double_t z2 = z+idz*TMath::Min(TMath::Sqrt(t.GetSigmaZ2()),1.);

	  xyz1[0] = x*TMath::Cos(t.GetAlpha())-y2*TMath::Sin(t.GetAlpha()); 
	  xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y2*TMath::Cos(t.GetAlpha());
	  xyz1[2] = z2;
	  AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
	  Float_t budget = param[0]*param[4];
	  meanbudget+=budget;
	  if (budget<minbudget) minbudget=budget;
	  if (budget>maxbudget) maxbudget=budget;
	}
      t.fBudget[0]+=minbudget;
      t.fBudget[1]+=meanbudget/9.;
      t.fBudget[2]+=minbudget;
      
      xyz1[0] = x*TMath::Cos(t.GetAlpha())-y*TMath::Sin(t.GetAlpha()); 
      xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y*TMath::Cos(t.GetAlpha());
      xyz1[2] = z;
      AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
      t.PropagateTo(x,param[1],param[0]);	
      t.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);   //end  global position	
      AdjustSector(&t);
      continue;
    }
    //
    //
    // stop tracking for highly inclined tracks 
    if (!AdjustSector(&t)) break;
    if (TMath::Abs(t.GetSnp())>0.95) break;
    //
    // propagate and update track in active layers
    //
    Int_t nr0 = nr;  //first active layer
    if (nr <outerTB && (fTrSec[0]->GetLayer(nr)->IsSensitive())){      
      Double_t xyz0[3],xyz1[3],param[7],x,y,z;
      t.GetGlobalXYZ(xyz0[0],xyz0[1],xyz0[2]);   //starting global position
      while (nr <outerTB && ((fTrSec[0]->GetLayer(nr)->IsSensitive()))){
	x = fTrSec[0]->GetLayer(nr)->GetX();
	nr++;
	if (!t.GetProlongation(x,y,z)) break;
	if (TMath::Abs(y)>(x*tanmax)){ 
	  nr++;
	  break;	  
	}
      }
      x = fTrSec[0]->GetLayer(nr)->GetX();
      if (!t.GetProlongation(x,y,z)) break;
       // minimal mean and maximal budget scan
      Float_t minbudget  =10000;
      Float_t meanbudget =0;
      Float_t maxbudget  =-1;
      //      Float_t normbudget =0;
      //     for (Int_t idy=-1;idy<=1;idy++)
      // 	for (Int_t idz=-1;idz<=1;idz++){
      for (Int_t idy=0;idy<1;idy++)
	for (Int_t idz=0;idz<1;idz++){
	  Double_t y2 = y+idy*TMath::Min(TMath::Sqrt(t.GetSigmaY2()),1.);
	  Double_t z2 = z+idz*TMath::Min(TMath::Sqrt(t.GetSigmaZ2()),1.);

	  xyz1[0] = x*TMath::Cos(t.GetAlpha())-y2*TMath::Sin(t.GetAlpha()); 
	  xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y2*TMath::Cos(t.GetAlpha());
	  xyz1[2] = z2;
	  AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);
	  Float_t budget = param[0]*param[4];
	  meanbudget+=budget;
	  if (budget<minbudget) minbudget=budget;
	  if (budget>maxbudget) maxbudget=budget;
	}
      t.fBudget[0]+=minbudget;
      t.fBudget[1]+=meanbudget/9.;
      t.fBudget[2]+=minbudget;
      //
      xyz1[0] = x*TMath::Cos(t.GetAlpha())-y*TMath::Sin(t.GetAlpha()); 
      xyz1[1] = +x*TMath::Sin(t.GetAlpha())+y*TMath::Cos(t.GetAlpha());
      xyz1[2] = z;
      // end global position
      AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
      rho = param[0];
      radLength = param[1];   // get mean propagation parameters
    }
    //
    //
    if (nr-nr0<5){
      // short tracklet - do not update - edge effect
      x = fTrSec[0]->GetLayer(nr+1)->GetX();
      t.PropagateTo(x,radLength,rho);
      AdjustSector(&t);
      continue; 
    }
    //
    //
    sector = t.GetSector();
    Float_t  ncl   = FindClusters(sector,nr0,nr,&t,clusters,tracklet);
    if (tracklet.GetN()-2*tracklet.GetNCross()<10) continue;
    ratio0 = ncl/Float_t(fTimeBinsPerPlane);
    Float_t  ratio1 = Float_t(t.fN+1)/Float_t(t.fNExpected+1.);	
    if (tracklet.GetChi2()<18.&&ratio0>0.8 && ratio1>0.6 && ratio0+ratio1>1.5 && t.GetNCross()==0 && TMath::Abs(t.GetSnp())<0.85&&t.fN>20){
      t.MakeBackupTrack();                            // make backup of the track until is gold
    }
    //
    //
    for (Int_t ilayer=nr0;ilayer<=nr;ilayer++) {
      expectedNumberOfClusters++;       
      t.fNExpected++;
      if (t.fX>345) t.fNExpectedLast++;
      AliTRDpropagationLayer& timeBin=*(fTrSec[sector]->GetLayer(ilayer));
      AliTRDcluster *cl=0;
      UInt_t index=0;
      Double_t maxChi2=fgkMaxChi2;
      dx = (fTrSec[sector]->GetLayer(ilayer-1))->GetX()-timeBin.GetX();
      x = timeBin.GetX();
      t.PropagateTo(x,radLength,rho);
      // Now go for the real cluster search
      if (timeBin) {	
	if (clusters[ilayer]>0) {
	  index = clusters[ilayer];
	  cl    = (AliTRDcluster*)GetCluster(index);
	  Double_t h01 = GetTiltFactor(cl);
          maxChi2=t.GetPredictedChi2(cl,h01);          
	}
	
        if (cl) {
	  //	  if (cl->GetNPads()<5) 
	    t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 
          Double_t h01 = GetTiltFactor(cl);
	  Int_t det = cl->GetDetector();    
	  Int_t plane = fGeom->GetPlane(det);
	  if (t.fX>345){
	    t.fNLast++;
	    t.fChi2Last+=maxChi2;
	  }
	  if(!t.UpdateMI(cl,maxChi2,index,h01,plane)) {
	    if(!t.Update(cl,maxChi2,index,h01)) {
	      //if(!tryAgain--) return 0;
	    }
          }  
          else tryAgain=fMaxGap;
	  //
	  
	  if (cl->GetLocalTimeBin()==1&&t.fN>20 && float(t.fChi2)/float(t.fN)<5){
	    Float_t  ratio1 = Float_t(t.fN)/Float_t(t.fNExpected);	
	    if (tracklet.GetChi2()<18&&ratio0>0.8&&ratio1>0.6 &&ratio0+ratio1>1.5 && t.GetNCross()==0 && TMath::Abs(t.GetSnp())<0.85){
	      t.MakeBackupTrack();                            // make backup of the track until is gold
	    }
	  }
	  // reset material budget if 2 consecutive gold
	  if (plane>0) 
	    if (t.fTracklets[plane].GetN()+t.fTracklets[plane-1].GetN()>20){
	      t.fBudget[2] = 0;
	    }	  
	}			
      }
    }  
  }
  //
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
  Double_t dxsample = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
    / AliTRDcalibDB::Instance()->GetSamplingFrequency();
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
    //if (cl->GetNPads()<5) 
      t.SetSampledEdx(TMath::Abs(cl->GetQ()/dxsample)); 

      //t.SetSampledEdx(cl->GetQ()/dx,t.GetNumberOfClusters()); 
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
  cout<<"\n Read Sectors  clusters"<<endl;
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
    //
    // apply pos correction
    Float_t poscor =  fgkCoef*(c->GetLocalTimeBin() - fgkMean)+fgkOffset; 
    c->SetY(c->GetY()-poscor);
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
//__________________________________________________________________________
void AliTRDtracker::MakeSeedsMI(Int_t /*inner*/, Int_t /*outer*/, AliESD * esd)
{
  //
  // Creates  seeds using clusters between  position inner plane  and outer plane 
  //
  const Double_t maxtheta = 1;
  const Double_t maxphi   = 2.0;
  //
  const Double_t kRoad0y  =  6;     // road for middle cluster 
  const Double_t kRoad0z  =  8.5;   // road for middle cluster 
  //
  const Double_t kRoad1y  =  2;    // road in y for seeded cluster
  const Double_t kRoad1z  =  20;    // road in z for seeded cluster
  //
  const Double_t kRoad2y  =  3;    // road in y for extrapolated cluster
  const Double_t kRoad2z  =  20;   // road in z for extrapolated cluster
  const Int_t    maxseed  = 3000;
  Int_t maxSec=AliTRDgeometry::kNsect;  

  //
  // linear fitters in planes
  TLinearFitter fitterTC(2,"hyp2");  // fitting with tilting pads - kz fixed - kz= Z/x, + vertex const
  TLinearFitter fitterT2(4,"hyp4");  // fitting with tilting pads - kz not fixed
  fitterTC.StoreData(kTRUE);
  fitterT2.StoreData(kTRUE);
  AliRieman rieman(1000);   // rieman fitter
  AliRieman rieman2(1000);   // rieman fitter
  //  
  // find the maximal and minimal layer for the planes
  //
  Int_t layers[6][2];
  AliTRDpropagationLayer* reflayers[6];
  for (Int_t i=0;i<6;i++){layers[i][0]=10000; layers[i][1]=0;}
  for (Int_t ns=0;ns<maxSec;ns++){
    for (Int_t ilayer=0;ilayer<fTrSec[ns]->GetNumberOfLayers();ilayer++){
      AliTRDpropagationLayer& layer=*(fTrSec[ns]->GetLayer(ilayer));
      if (layer==0) continue;
      Int_t det   = layer[0]->GetDetector();    
      Int_t plane = fGeom->GetPlane(det);
      if (ilayer<layers[plane][0]) layers[plane][0] = ilayer;
      if (ilayer>layers[plane][1]) layers[plane][1] = ilayer;
    }
  }
  //
  AliTRDpadPlane *padPlane = AliTRDCommonParam::Instance()->GetPadPlane(0,0);
  Double_t h01 = TMath::Tan(-TMath::Pi() / 180.0 * padPlane->GetTiltingAngle());
  Double_t hL[6];         // tilting angle
  Double_t xcl[6];        // x - position of reference cluster
  Double_t ycl[6];        // y - position of reference cluster
  Double_t zcl[6];        // z - position of reference cluster
  AliTRDcluster *cl[6]={0,0,0,0,0,0};    // seeding clusters
  Float_t padlength[6]={10,10,10,10,10,10};   //current pad-length 
  Double_t chi2R =0, chi2Z=0;
  Double_t chi2RF =0, chi2ZF=0;
  //
  Int_t nclusters;     // total number of clusters
  for (Int_t i=0;i<6;i++) {hL[i]=h01; if (i%2==1) hL[i]*=-1.;}
  //
  //
  //         registered seed
  AliTRDseed *pseed = new AliTRDseed[maxseed*6];
  AliTRDseed *seed[maxseed];
  for (Int_t iseed=0;iseed<maxseed;iseed++) seed[iseed]= &pseed[iseed*6];
  AliTRDseed *cseed = seed[0];
  // 
  Double_t   seedquality[maxseed];  
  Double_t   seedquality2[maxseed];  
  Double_t   seedparams[maxseed][7];
  Int_t      seedlayer[maxseed];
  Int_t      registered =0;
  Int_t      sort[maxseed];
  //
  // seeding part
  //
  for (Int_t ns = 0; ns<maxSec; ns++){         //loop over sectors
  //for (Int_t ns = 0; ns<5; ns++){         //loop over sectors
    registered = 0;   // reset registerd seed counter
    cseed      = seed[registered];
    Float_t iter=0;
    for (Int_t sLayer=2; sLayer>=0;sLayer--){
      //for (Int_t dseed=5;dseed<15; dseed+=3){  //loop over central seeding time bins 
      iter+=1.;
      Int_t dseed = 5+Int_t(iter)*3;
      // Initialize seeding layers
      for (Int_t ilayer=0;ilayer<6;ilayer++){
	reflayers[ilayer] = fTrSec[ns]->GetLayer(layers[ilayer][1]-dseed);
	xcl[ilayer]       = reflayers[ilayer]->GetX();
      }      
      //
      Double_t xref                 = (xcl[sLayer+1] + xcl[sLayer+2])*0.5;      
      AliTRDpropagationLayer& layer0=*reflayers[sLayer+0];
      AliTRDpropagationLayer& layer1=*reflayers[sLayer+1];
      AliTRDpropagationLayer& layer2=*reflayers[sLayer+2];
      AliTRDpropagationLayer& layer3=*reflayers[sLayer+3];
      //
      Int_t maxn3  = layer3;
      for (Int_t icl3=0;icl3<maxn3;icl3++){
	AliTRDcluster *cl3 = layer3[icl3];
	if (!cl3) continue;	
	padlength[sLayer+3] = TMath::Sqrt(cl3->GetSigmaZ2()*12.);
	ycl[sLayer+3] = cl3->GetY();
	zcl[sLayer+3] = cl3->GetZ();
	Float_t yymin0 = ycl[sLayer+3] - 1- maxphi *(xcl[sLayer+3]-xcl[sLayer+0]);
	Float_t yymax0 = ycl[sLayer+3] + 1+ maxphi *(xcl[sLayer+3]-xcl[sLayer+0]);
	Int_t   maxn0 = layer0;  // 
	for (Int_t icl0=layer0.Find(yymin0);icl0<maxn0;icl0++){
	  AliTRDcluster *cl0 = layer0[icl0];
	  if (!cl0) continue;
	  if (cl3->IsUsed()&&cl0->IsUsed()) continue;
	  ycl[sLayer+0] = cl0->GetY();
	  zcl[sLayer+0] = cl0->GetZ();
	  if ( ycl[sLayer+0]>yymax0) break;
	  Double_t tanphi   = (ycl[sLayer+3]-ycl[sLayer+0])/(xcl[sLayer+3]-xcl[sLayer+0]); 
	  if (TMath::Abs(tanphi)>maxphi) continue;
	  Double_t tantheta = (zcl[sLayer+3]-zcl[sLayer+0])/(xcl[sLayer+3]-xcl[sLayer+0]); 
	  if (TMath::Abs(tantheta)>maxtheta) continue; 
	  padlength[sLayer+0] = TMath::Sqrt(cl0->GetSigmaZ2()*12.);
	  //
	  // expected position in 1 layer
	  Double_t y1exp = ycl[sLayer+0]+(tanphi)  *(xcl[sLayer+1]-xcl[sLayer+0]);	  
	  Double_t z1exp = zcl[sLayer+0]+(tantheta)*(xcl[sLayer+1]-xcl[sLayer+0]);	  
	  Float_t yymin1 = y1exp - kRoad0y-tanphi;
	  Float_t yymax1 = y1exp + kRoad0y+tanphi;
	  Int_t   maxn1  = layer1;  // 
	  //
	  for (Int_t icl1=layer1.Find(yymin1);icl1<maxn1;icl1++){
	    AliTRDcluster *cl1 = layer1[icl1];
	    if (!cl1) continue;
	    Int_t nusedCl = 0;
	    if (cl3->IsUsed()) nusedCl++;
	    if (cl0->IsUsed()) nusedCl++;
	    if (cl1->IsUsed()) nusedCl++;
	    if (nusedCl>1) continue;
	    ycl[sLayer+1] = cl1->GetY();
	    zcl[sLayer+1] = cl1->GetZ();
	    if ( ycl[sLayer+1]>yymax1) break;
	    if (TMath::Abs(ycl[sLayer+1]-y1exp)>kRoad0y+tanphi) continue;
	    if (TMath::Abs(zcl[sLayer+1]-z1exp)>kRoad0z)        continue;
	    padlength[sLayer+1] = TMath::Sqrt(cl1->GetSigmaZ2()*12.);
	    //
	    Double_t y2exp  = ycl[sLayer+0]+(tanphi)  *(xcl[sLayer+2]-xcl[sLayer+0])+(ycl[sLayer+1]-y1exp);	  
	    Double_t z2exp  = zcl[sLayer+0]+(tantheta)*(xcl[sLayer+2]-xcl[sLayer+0]);
	    Int_t    index2 = layer2.FindNearestCluster(y2exp,z2exp,kRoad1y,  kRoad1z);
	    if (index2<=0) continue; 
	    AliTRDcluster *cl2 = (AliTRDcluster*)GetCluster(index2);
	    padlength[sLayer+2] = TMath::Sqrt(cl2->GetSigmaZ2()*12.);
	    ycl[sLayer+2] = cl2->GetY();
	    zcl[sLayer+2] = cl2->GetZ();
	    if (TMath::Abs(cl2->GetZ()-z2exp)>kRoad0z)        continue;
	    //
 	    rieman.Reset();
 	    rieman.AddPoint(xcl[sLayer+0],ycl[sLayer+0],zcl[sLayer+0],1,10);
 	    rieman.AddPoint(xcl[sLayer+1],ycl[sLayer+1],zcl[sLayer+1],1,10);
 	    rieman.AddPoint(xcl[sLayer+3],ycl[sLayer+3],zcl[sLayer+3],1,10);	    
	    rieman.AddPoint(xcl[sLayer+2],ycl[sLayer+2],zcl[sLayer+2],1,10);
	    rieman.Update();
	    //
	    // reset fitter
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      cseed[iLayer].Reset();
	    }	  
	    chi2Z =0.; chi2R=0.;
	    for (Int_t iLayer=0;iLayer<4;iLayer++){
	      cseed[sLayer+iLayer].fZref[0] = rieman.GetZat(xcl[sLayer+iLayer]);
	      chi2Z += (cseed[sLayer+iLayer].fZref[0]- zcl[sLayer+iLayer])*
		(cseed[sLayer+iLayer].fZref[0]- zcl[sLayer+iLayer]);
	      cseed[sLayer+iLayer].fZref[1] = rieman.GetDZat(xcl[sLayer+iLayer]);	      
	      cseed[sLayer+iLayer].fYref[0] = rieman.GetYat(xcl[sLayer+iLayer]);
	      chi2R += (cseed[sLayer+iLayer].fYref[0]- ycl[sLayer+iLayer])*
		(cseed[sLayer+iLayer].fYref[0]- ycl[sLayer+iLayer]);
	      cseed[sLayer+iLayer].fYref[1] = rieman.GetDYat(xcl[sLayer+iLayer]);
	    }
	    if (TMath::Sqrt(chi2R)>1./iter) continue;
	    if (TMath::Sqrt(chi2Z)>7./iter) continue;
	    //
	    //
	    //
	    Float_t minmax[2]={-100,100};
	    for (Int_t iLayer=0;iLayer<4;iLayer++){
	      Float_t max = zcl[sLayer+iLayer]+padlength[sLayer+iLayer]*0.5+1 -cseed[sLayer+iLayer].fZref[0];
	      if (max<minmax[1]) minmax[1]=max; 
	      Float_t min = zcl[sLayer+iLayer]-padlength[sLayer+iLayer]*0.5-1 -cseed[sLayer+iLayer].fZref[0];
	      if (min>minmax[0]) minmax[0]=min; 
	    }
	    Bool_t isFake = kFALSE; 
	    if (cl0->GetLabel(0)!=cl3->GetLabel(0)) isFake = kTRUE;
	    if (cl1->GetLabel(0)!=cl3->GetLabel(0)) isFake = kTRUE;
	    if (cl2->GetLabel(0)!=cl3->GetLabel(0)) isFake = kTRUE;
	    if ((!isFake) || (icl3%10)==0 ){  //debugging print
	      TTreeSRedirector& cstream = *fDebugStreamer;
	      cstream<<"Seeds0"<<
		"isFake="<<isFake<<
		"Cl0.="<<cl0<<
		"Cl1.="<<cl1<<
		"Cl2.="<<cl2<<
		"Cl3.="<<cl3<<
		"Xref="<<xref<<
		"X0="<<xcl[sLayer+0]<<
		"X1="<<xcl[sLayer+1]<<
		"X2="<<xcl[sLayer+2]<<
		"X3="<<xcl[sLayer+3]<<
		"Y2exp="<<y2exp<<
		"Z2exp="<<z2exp<<
		"Chi2R="<<chi2R<<
		"Chi2Z="<<chi2Z<<		
		"Seed0.="<<&cseed[sLayer+0]<<
		"Seed1.="<<&cseed[sLayer+1]<<
		"Seed2.="<<&cseed[sLayer+2]<<
		"Seed3.="<<&cseed[sLayer+3]<<
		"Zmin="<<minmax[0]<<
		"Zmax="<<minmax[1]<<
		"\n";
	    }
	    
	    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //<<<<<<<<<<<<<<<<<<    FIT SEEDING PART                  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    cl[sLayer+0] = cl0;
	    cl[sLayer+1] = cl1;
	    cl[sLayer+2] = cl2;
	    cl[sLayer+3] = cl3;
	    Bool_t isOK=kTRUE;
	    for (Int_t jLayer=0;jLayer<4;jLayer++){
	      cseed[sLayer+jLayer].fTilt = hL[sLayer+jLayer];
	      cseed[sLayer+jLayer].fPadLength = padlength[sLayer+jLayer];
	      cseed[sLayer+jLayer].fX0   = xcl[sLayer+jLayer];
	      for (Int_t iter=0; iter<2; iter++){
		//
		// in iteration 0 we try only one pad-row
		// if quality not sufficient we try 2 pad-rows - about 5% of tracks cross 2 pad-rows
		//
		AliTRDseed tseed = cseed[sLayer+jLayer];
		Float_t    roadz  = padlength[sLayer+jLayer]*0.5;
		if (iter>0) roadz = padlength[sLayer+jLayer];
		//
		Float_t quality =10000;
		for (Int_t iTime=2;iTime<20;iTime++){ 
		  AliTRDpropagationLayer& layer = *(fTrSec[ns]->GetLayer(layers[sLayer+jLayer][1]-iTime));
		  Double_t dxlayer= layer.GetX()-xcl[sLayer+jLayer];		 
		  Double_t zexp   = cl[sLayer+jLayer]->GetZ() ;
		  if (iter>0){
		    // try 2 pad-rows in second iteration
		    zexp  = tseed.fZref[0]+ tseed.fZref[1]*dxlayer;
		    if (zexp>cl[sLayer+jLayer]->GetZ()) zexp = cl[sLayer+jLayer]->GetZ()+padlength[sLayer+jLayer]*0.5;
		    if (zexp<cl[sLayer+jLayer]->GetZ()) zexp = cl[sLayer+jLayer]->GetZ()-padlength[sLayer+jLayer]*0.5;
		  }
		  //
		  Double_t yexp  =  tseed.fYref[0]+ 
		    tseed.fYref[1]*dxlayer;
		  Int_t    index = layer.FindNearestCluster(yexp,zexp,kRoad1y, roadz);
		  if (index<=0) continue; 
		  AliTRDcluster *cl = (AliTRDcluster*)GetCluster(index);	      
		  //
		  tseed.fIndexes[iTime]  = index;
		  tseed.fClusters[iTime] = cl;   // register cluster
		  tseed.fX[iTime] = dxlayer;     // register cluster
		  tseed.fY[iTime] = cl->GetY();  // register cluster
		  tseed.fZ[iTime] = cl->GetZ();  // register cluster
		} 
		tseed.Update();
		//count the number of clusters and distortions into quality
		Float_t dangle = tseed.fYfit[1]-tseed.fYref[1];
		Float_t tquality   = (18-tseed.fN2)/2. + TMath::Abs(dangle)/0.1+
		  TMath::Abs(tseed.fYfit[0]-tseed.fYref[0])/0.2+
		  2.*TMath::Abs(tseed.fMeanz-tseed.fZref[0])/padlength[jLayer];
		if (iter==0 && tseed.isOK()) {
		  cseed[sLayer+jLayer] = tseed;
		  quality = tquality;
		  if (tquality<5) break;  
		}
		if (tseed.isOK() && tquality<quality)
		  cseed[sLayer+jLayer] = tseed;				
	      }
	      if (!cseed[sLayer+jLayer].isOK()){
		isOK = kFALSE;
		break;
	      }	     	          
	      cseed[sLayer+jLayer].CookLabels();
	      cseed[sLayer+jLayer].UpdateUsed();
	      nusedCl+= cseed[sLayer+jLayer].fNUsed;
	      if (nusedCl>25){
		isOK = kFALSE;
		break;
	      }	    
	    }
	    //
	    if (!isOK) continue;
	    nclusters=0;
	    for (Int_t iLayer=0;iLayer<4;iLayer++){
	      if (cseed[sLayer+iLayer].isOK()){
		nclusters+=cseed[sLayer+iLayer].fN2;	    
	      }
	    }
	    // 
	    // iteration 0
	    rieman.Reset();
	    for (Int_t iLayer=0;iLayer<4;iLayer++){
	      rieman.AddPoint(xcl[sLayer+iLayer],cseed[sLayer+iLayer].fYfitR[0],
			      cseed[sLayer+iLayer].fZProb,1,10);
	    }
	    rieman.Update();
	    //
	    //
	    chi2R =0; chi2Z=0;
	    for (Int_t iLayer=0;iLayer<4;iLayer++){
	      cseed[sLayer+iLayer].fYref[0] = rieman.GetYat(xcl[sLayer+iLayer]);
	      chi2R += (cseed[sLayer+iLayer].fYref[0]-cseed[sLayer+iLayer].fYfitR[0])*
		(cseed[sLayer+iLayer].fYref[0]-cseed[sLayer+iLayer].fYfitR[0]);
	      cseed[sLayer+iLayer].fYref[1] = rieman.GetDYat(xcl[sLayer+iLayer]);
	      cseed[sLayer+iLayer].fZref[0] = rieman.GetZat(xcl[sLayer+iLayer]);
	      chi2Z += (cseed[sLayer+iLayer].fZref[0]- cseed[sLayer+iLayer].fMeanz)*
		(cseed[sLayer+iLayer].fZref[0]- cseed[sLayer+iLayer].fMeanz);
	      cseed[sLayer+iLayer].fZref[1] = rieman.GetDZat(xcl[sLayer+iLayer]);
	    }
	    Double_t curv = rieman.GetC();
	    //
	    // likelihoods
	    //
	    Double_t sumda = 
	      TMath::Abs(cseed[sLayer+0].fYfitR[1]- cseed[sLayer+0].fYref[1])+
	      TMath::Abs(cseed[sLayer+1].fYfitR[1]- cseed[sLayer+1].fYref[1])+
	      TMath::Abs(cseed[sLayer+2].fYfitR[1]- cseed[sLayer+2].fYref[1])+
	      TMath::Abs(cseed[sLayer+3].fYfitR[1]- cseed[sLayer+3].fYref[1]);
	    Double_t likea = TMath::Exp(-sumda*10.6);
	    Double_t likechi2 = 0.0000000001;
	    if (chi2R<0.5) likechi2+=TMath::Exp(-TMath::Sqrt(chi2R)*7.73);
	    Double_t likechi2z = TMath::Exp(-chi2Z*0.088)/TMath::Exp(-chi2Z*0.019);
	    Double_t likeN    = TMath::Exp(-(72-nclusters)*0.19);
	    Double_t like     = likea*likechi2*likechi2z*likeN;
	    //
	    Double_t likePrimY = TMath::Exp(-TMath::Abs(cseed[sLayer+0].fYref[1]-130*curv)*1.9);
	    Double_t likePrimZ = TMath::Exp(-TMath::Abs(cseed[sLayer+0].fZref[1]-
							cseed[sLayer+0].fZref[0]/xcl[sLayer+0])*5.9);
	    Double_t likePrim  = TMath::Max(likePrimY*likePrimZ,0.0005);
					    
	    seedquality[registered]  = like; 
	    seedlayer[registered]    = sLayer;
	    if (TMath::Log(0.000000000000001+like)<-15) continue;
	    AliTRDseed seedb[6];
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      seedb[iLayer] = cseed[iLayer]; 
	    }
	    //
	    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //<<<<<<<<<<<<<<<   FULL TRACK FIT PART         <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	    //
	    Int_t nlayers            = 0;
	    Int_t nusedf             = 0;
	    Int_t findable           = 0;
	    //
	    // add new layers  - avoid long extrapolation
	    //
	    Int_t tLayer[2]={0,0};
	    if (sLayer==2) {tLayer[0]=1; tLayer[1]=0;}
	    if (sLayer==1) {tLayer[0]=5; tLayer[1]=0;}
	    if (sLayer==0) {tLayer[0]=4; tLayer[1]=5;}
	    //
	    for (Int_t iLayer=0;iLayer<2;iLayer++){
	      Int_t jLayer = tLayer[iLayer];      // set tracking layer	      
	      cseed[jLayer].Reset();
	      cseed[jLayer].fTilt    = hL[jLayer];
	      cseed[jLayer].fPadLength = padlength[jLayer];
	      cseed[jLayer].fX0      = xcl[jLayer];
	      // get pad length and rough cluster
	      Int_t    indexdummy = reflayers[jLayer]->FindNearestCluster(cseed[jLayer].fYref[0], 
									  cseed[jLayer].fZref[0],kRoad2y,kRoad2z);
	      if (indexdummy<=0) continue; 
	      AliTRDcluster *cldummy = (AliTRDcluster*)GetCluster(indexdummy);
	      padlength[jLayer]      = TMath::Sqrt(cldummy->GetSigmaZ2()*12.);
	    }
	    AliTRDseed::FitRiemanTilt(cseed, kTRUE);
	    //
	    for (Int_t iLayer=0;iLayer<2;iLayer++){
	      Int_t jLayer = tLayer[iLayer];      // set tracking layer	
	      if ( (jLayer==0) && !(cseed[1].isOK())) continue;  // break not allowed
	      if ( (jLayer==5) && !(cseed[4].isOK())) continue;  // break not allowed
	      Float_t  zexp  = cseed[jLayer].fZref[0];
	      Double_t zroad =  padlength[jLayer]*0.5+1.;
	      //
	      // 
	      for (Int_t iter=0;iter<2;iter++){
		AliTRDseed tseed = cseed[jLayer];
		Float_t quality = 10000;
		for (Int_t iTime=2;iTime<20;iTime++){ 
		  AliTRDpropagationLayer& layer = *(fTrSec[ns]->GetLayer(layers[jLayer][1]-iTime));
		  Double_t dxlayer     = layer.GetX()-xcl[jLayer];
		  Double_t yexp        = tseed.fYref[0]+tseed.fYref[1]*dxlayer;
		  Float_t  yroad       = kRoad1y;
		  Int_t    index = layer.FindNearestCluster(yexp,zexp, yroad, zroad);
		  if (index<=0) continue; 
		  AliTRDcluster *cl = (AliTRDcluster*)GetCluster(index);	      
		  //
		  tseed.fIndexes[iTime]  = index;
		  tseed.fClusters[iTime] = cl;   // register cluster
		  tseed.fX[iTime] = dxlayer;     // register cluster
		  tseed.fY[iTime] = cl->GetY();  // register cluster
		  tseed.fZ[iTime] = cl->GetZ();  // register cluster
		} 	      
		tseed.Update();
		if (tseed.isOK()){
		  Float_t dangle = tseed.fYfit[1]-tseed.fYref[1];
		  Float_t tquality   = (18-tseed.fN2)/2. + TMath::Abs(dangle)/0.1+		  
		    TMath::Abs(tseed.fYfit[0]-tseed.fYref[0])/0.2+ 
		    2.*TMath::Abs(tseed.fMeanz-tseed.fZref[0])/padlength[jLayer];
		  //
		  if (tquality<quality){
		    cseed[jLayer]=tseed;
		    quality = tquality;
		  }
		}
		zroad*=2.;
	      }
	      if ( cseed[jLayer].isOK()){
		cseed[jLayer].CookLabels();
		cseed[jLayer].UpdateUsed();
		nusedf+= cseed[jLayer].fNUsed;
		AliTRDseed::FitRiemanTilt(cseed, kTRUE);
	      }
	    }
	    //
	    //
	    // make copy
	    AliTRDseed bseed[6];
	    for (Int_t jLayer=0;jLayer<6;jLayer++){
	      bseed[jLayer] = cseed[jLayer];
	    }	    
	    Float_t lastquality = 10000;
	    Float_t lastchi2    = 10000;
	    Float_t chi2        = 1000;

	    //
	    for (Int_t iter =0; iter<4;iter++){
	      //
	      // sort tracklets according "quality", try to "improve" 4 worst 
	      //
	      Float_t sumquality = 0;
	      Float_t squality[6];
	      Int_t   sortindexes[6];
	      for (Int_t jLayer=0;jLayer<6;jLayer++){
		if (bseed[jLayer].isOK()){ 
		  AliTRDseed &tseed = bseed[jLayer];
		  Double_t zcor  =  tseed.fTilt*(tseed.fZProb-tseed.fZref[0]);
		  Float_t dangle = tseed.fYfit[1]-tseed.fYref[1];
		  Float_t tquality  = (18-tseed.fN2)/2. + TMath::Abs(dangle)/0.1+		  
		    TMath::Abs(tseed.fYfit[0]-(tseed.fYref[0]-zcor))/0.2+ 
		    2.*TMath::Abs(tseed.fMeanz-tseed.fZref[0])/padlength[jLayer];
		  squality[jLayer] = tquality;
		}
		else  squality[jLayer]=-1;
		sumquality +=squality[jLayer];
	      }

	      if (sumquality>=lastquality ||  chi2>lastchi2) break;
	      lastquality = sumquality;	 
	      lastchi2    = chi2;
	      if (iter>0){
		for (Int_t jLayer=0;jLayer<6;jLayer++){
		  cseed[jLayer] = bseed[jLayer];
		}
	      }
	      TMath::Sort(6,squality,sortindexes,kFALSE);
	      //
	      //
	      for (Int_t jLayer=5;jLayer>1;jLayer--){
		Int_t bLayer = sortindexes[jLayer];
		AliTRDseed tseed = bseed[bLayer];
		for (Int_t iTime=2;iTime<20;iTime++){ 
		  AliTRDpropagationLayer& layer = *(fTrSec[ns]->GetLayer(layers[bLayer][1]-iTime));
		  Double_t dxlayer= layer.GetX()-xcl[bLayer];
		  //
		  Double_t zexp  =  tseed.fZref[0];
		  Double_t zcor  =  tseed.fTilt*(tseed.fZProb-tseed.fZref[0]);
		  //
		  Float_t  roadz = padlength[bLayer]+1;
 		  if (TMath::Abs(tseed.fZProb-zexp)> padlength[bLayer]*0.5) {roadz = padlength[bLayer]*0.5;}
		  if (tseed.fZfit[1]*tseed.fZref[1]<0) {roadz = padlength[bLayer]*0.5;}
		  if (TMath::Abs(tseed.fZProb-zexp)<0.1*padlength[bLayer]) {
		    zexp = tseed.fZProb; 
		    roadz = padlength[bLayer]*0.5;
		  }
		  //
		  Double_t yexp  =  tseed.fYref[0]+ 
		    tseed.fYref[1]*dxlayer-zcor;
		  Int_t    index = layer.FindNearestCluster(yexp,zexp,kRoad1y, roadz);
		  if (index<=0) continue; 
		  AliTRDcluster *cl = (AliTRDcluster*)GetCluster(index);	      
		  //
		  tseed.fIndexes[iTime]  = index;
		  tseed.fClusters[iTime] = cl;   // register cluster
		  tseed.fX[iTime] = dxlayer;     // register cluster
		  tseed.fY[iTime] = cl->GetY();  // register cluster
		  tseed.fZ[iTime] = cl->GetZ();  // register cluster
		} 
		tseed.Update();
		if (tseed.isOK()) {
		  Float_t dangle = tseed.fYfit[1]-tseed.fYref[1];
		  Double_t zcor  =  tseed.fTilt*(tseed.fZProb-tseed.fZref[0]);
		  //
		  Float_t tquality   = (18-tseed.fN2)/2. + TMath::Abs(dangle)/0.1+		  
		    TMath::Abs(tseed.fYfit[0]-(tseed.fYref[0]-zcor))/0.2+ 
		    2.*TMath::Abs(tseed.fMeanz-tseed.fZref[0])/padlength[jLayer];
		  //
		  if (tquality<squality[bLayer])
		    bseed[bLayer] = tseed;
		}
	      }
	      chi2 = AliTRDseed::FitRiemanTilt(bseed, kTRUE);
	    }
	    //
	    //
	    //
	    nclusters  = 0;
	    nlayers    = 0;
	    findable   = 0;
	    for (Int_t iLayer=0;iLayer<6;iLayer++) {
	      if (TMath::Abs(cseed[iLayer].fYref[0]/cseed[iLayer].fX0)<0.15)
		findable++;
	      if (cseed[iLayer].isOK()){
		nclusters+=cseed[iLayer].fN2;	    
		nlayers++;
	      }
	    }
	    if (nlayers<3) continue;
	    rieman.Reset();
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      if (cseed[iLayer].isOK()) rieman.AddPoint(xcl[iLayer],cseed[iLayer].fYfitR[0],
								   cseed[iLayer].fZProb,1,10);
	    }
	    rieman.Update();
	    //
	    chi2RF =0;
	    chi2ZF =0;
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      if (cseed[iLayer].isOK()){
		cseed[iLayer].fYref[0] = rieman.GetYat(xcl[iLayer]);
		chi2RF += (cseed[iLayer].fYref[0]-cseed[iLayer].fYfitR[0])*
		  (cseed[iLayer].fYref[0]-cseed[iLayer].fYfitR[0]);
		cseed[iLayer].fYref[1] = rieman.GetDYat(xcl[iLayer]);
		cseed[iLayer].fZref[0] = rieman.GetZat(xcl[iLayer]);
		chi2ZF += (cseed[iLayer].fZref[0]- cseed[iLayer].fMeanz)*
		  (cseed[iLayer].fZref[0]- cseed[iLayer].fMeanz);
		cseed[iLayer].fZref[1] = rieman.GetDZat(xcl[iLayer]);
	      }
	    }
	    chi2RF/=TMath::Max((nlayers-3.),1.);
	    chi2ZF/=TMath::Max((nlayers-3.),1.);
	    curv = rieman.GetC();
	    
	    //

	    Double_t xref2    = (xcl[2]+xcl[3])*0.5;  // middle of the chamber
	    Double_t dzmf     = rieman.GetDZat(xref2);
	    Double_t zmf      = rieman.GetZat(xref2);
	    //
	    // fit hyperplane
	    //
	    Int_t npointsT =0;
	    fitterTC.ClearPoints();
	    fitterT2.ClearPoints();
	    rieman2.Reset();
	    for (Int_t iLayer=0; iLayer<6;iLayer++){
	      if (!cseed[iLayer].isOK()) continue;
	      for (Int_t itime=0;itime<25;itime++){
		if (!cseed[iLayer].fUsable[itime]) continue;
		Double_t x   = cseed[iLayer].fX[itime]+cseed[iLayer].fX0-xref2;  // x relative to the midle chamber
		Double_t y   = cseed[iLayer].fY[itime];
		Double_t z   = cseed[iLayer].fZ[itime];
		// ExB correction to the correction
		// tilted rieman
		//
		Double_t uvt[6];
		Double_t x2 = cseed[iLayer].fX[itime]+cseed[iLayer].fX0;      // global x
		//		
		Double_t t = 1./(x2*x2+y*y);
		uvt[1]  = t;    // t
		uvt[0]  = 2.*x2*uvt[1];      // u 
		//
		uvt[2]  = 2.0*hL[iLayer]*uvt[1];
		uvt[3]  = 2.0*hL[iLayer]*x*uvt[1];	      
		uvt[4]  = 2.0*(y+hL[iLayer]*z)*uvt[1];
		//
		Double_t error = 2*0.2*uvt[1];
		fitterT2.AddPoint(uvt,uvt[4],error);
		//
		// constrained rieman
		// 
		z =cseed[iLayer].fZ[itime];
		uvt[0]  = 2.*x2*t;           // u 
		uvt[1]  = 2*hL[iLayer]*x2*uvt[1];	      
		uvt[2]  = 2*(y+hL[iLayer]*(z-GetZ()))*t;
		fitterTC.AddPoint(uvt,uvt[2],error);
		//		
		rieman2.AddPoint(x2,y,z,1,10);
		npointsT++;
	      }
	    }
	    rieman2.Update();
	    fitterTC.Eval();
	    fitterT2.Eval();
	    Double_t rpolz0 = fitterT2.GetParameter(3);
	    Double_t rpolz1 = fitterT2.GetParameter(4);	    
	    //
	    // linear fitter  - not possible to make boundaries
	    // non accept non possible z and dzdx combination
	    // 	    
	    Bool_t   acceptablez =kTRUE;
	    for (Int_t iLayer=0; iLayer<6;iLayer++){
	      if (cseed[iLayer].isOK()){
		Double_t zT2 =  rpolz0+rpolz1*(xcl[iLayer] - xref2);
		if (TMath::Abs(cseed[iLayer].fZProb-zT2)>padlength[iLayer]*0.5+1)
		  acceptablez = kFALSE;
	      }
	    }
	    if (!acceptablez){
	      fitterT2.FixParameter(3,zmf);
	      fitterT2.FixParameter(4,dzmf);
	      fitterT2.Eval();
	      fitterT2.ReleaseParameter(3);
	      fitterT2.ReleaseParameter(4);
	      rpolz0 = fitterT2.GetParameter(3);
	      rpolz1 = fitterT2.GetParameter(4);
	    }
	    //
	    Double_t chi2TR = fitterT2.GetChisquare()/Float_t(npointsT);
	    Double_t chi2TC = fitterTC.GetChisquare()/Float_t(npointsT);
	    //
	    Double_t polz1c = fitterTC.GetParameter(2);
	    Double_t polz0c = polz1c*xref2;
	    //
	    Double_t aC     =  fitterTC.GetParameter(0);
	    Double_t bC     =  fitterTC.GetParameter(1);
	    Double_t CC     =  aC/TMath::Sqrt(bC*bC+1.);     // curvature
	    //
	    Double_t aR     =  fitterT2.GetParameter(0);
	    Double_t bR     =  fitterT2.GetParameter(1);
	    Double_t dR     =  fitterT2.GetParameter(2);	    
	    Double_t CR     =  1+bR*bR-dR*aR;
	    Double_t dca    =  0.;	    
	    if (CR>0){
	      dca = -dR/(TMath::Sqrt(1+bR*bR-dR*aR)+TMath::Sqrt(1+bR*bR)); 
	      CR  = aR/TMath::Sqrt(CR);
	    }
	    //
	    Double_t chi2ZT2=0, chi2ZTC=0;
	    for (Int_t iLayer=0; iLayer<6;iLayer++){
	      if (cseed[iLayer].isOK()){
		Double_t zT2 =  rpolz0+rpolz1*(xcl[iLayer] - xref2);
		Double_t zTC =  polz0c+polz1c*(xcl[iLayer] - xref2);
		chi2ZT2 += TMath::Abs(cseed[iLayer].fMeanz-zT2);
		chi2ZTC += TMath::Abs(cseed[iLayer].fMeanz-zTC);
	      }
	    }
	    chi2ZT2/=TMath::Max((nlayers-3.),1.);
	    chi2ZTC/=TMath::Max((nlayers-3.),1.);	    
	    //
	    //
	    //
	    AliTRDseed::FitRiemanTilt(cseed, kTRUE);
	    Float_t sumdaf = 0;
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      if (cseed[iLayer].isOK())
		sumdaf += TMath::Abs((cseed[iLayer].fYfit[1]-cseed[iLayer].fYref[1])/cseed[iLayer].fSigmaY2);
	    }  
	    sumdaf /= Float_t (nlayers-2.);
	    //
	    // likelihoods for full track
	    //
	    Double_t likezf      = TMath::Exp(-chi2ZF*0.14);
	    Double_t likechi2C   = TMath::Exp(-chi2TC*0.677);
	    Double_t likechi2TR  = TMath::Exp(-chi2TR*0.78);
	    Double_t likeaf      = TMath::Exp(-sumdaf*3.23);
	    seedquality2[registered] = likezf*likechi2TR*likeaf; 
//	    Bool_t isGold = kFALSE;
//	    
// 	    if (nlayers == 6        && TMath::Log(0.000000001+seedquality2[index])<-5.) isGold =kTRUE;   // gold
// 	    if (nlayers == findable && TMath::Log(0.000000001+seedquality2[index])<-4.) isGold =kTRUE;   // gold
// 	    if (isGold &&nusedf<10){
// 	      for (Int_t jLayer=0;jLayer<6;jLayer++){
// 		if ( seed[index][jLayer].isOK()&&TMath::Abs(seed[index][jLayer].fYfit[1]-seed[index][jLayer].fYfit[1])<0.1)
// 		  seed[index][jLayer].UseClusters();  //sign gold
// 	      }
// 	    }
	    //
	    //
	    //
	    Int_t index0=0;
	    if (!cseed[0].isOK()){
	      index0 = 1;
	      if (!cseed[1].isOK()) index0 = 2;
	    }
	    seedparams[registered][0] = cseed[index0].fX0;
	    seedparams[registered][1] = cseed[index0].fYref[0];
	    seedparams[registered][2] = cseed[index0].fZref[0];
	    seedparams[registered][5] = CR;
	    seedparams[registered][3] = cseed[index0].fX0*CR - TMath::Sin(TMath::ATan(cseed[0].fYref[1]));
	    seedparams[registered][4] = cseed[index0].fZref[1]/	      
	      TMath::Sqrt(1+cseed[index0].fYref[1]*cseed[index0].fYref[1]);
	    seedparams[registered][6] = ns;
	    //
	    //
	    Int_t labels[12], outlab[24];
	    Int_t nlab=0;
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      if (!cseed[iLayer].isOK()) continue;
	      if (cseed[iLayer].fLabels[0]>=0) {
		labels[nlab] = cseed[iLayer].fLabels[0];
		nlab++;
	      }
	      if (cseed[iLayer].fLabels[1]>=0) {
		labels[nlab] = cseed[iLayer].fLabels[1];
		nlab++;
	      }	      
	    }
	    Freq(nlab,labels,outlab,kFALSE);
	    Int_t label = outlab[0];
	    Int_t frequency  = outlab[1];
	    for (Int_t iLayer=0;iLayer<6;iLayer++){
	      cseed[iLayer].fFreq  = frequency;
	      cseed[iLayer].fC     = CR;
	      cseed[iLayer].fCC     = CC;
	      cseed[iLayer].fChi2  = chi2TR;
	      cseed[iLayer].fChi2Z = chi2ZF;
	    }
	    //
	    if (1||(!isFake)){  //debugging print
	      Float_t zvertex = GetZ();
	      TTreeSRedirector& cstream = *fDebugStreamer;
	      cstream<<"Seeds1"<<
		"isFake="<<isFake<<
		"Vertex="<<zvertex<<
		"Rieman2.="<<&rieman2<<
		"Rieman.="<<&rieman<<
		"Xref="<<xref<<
		"X0="<<xcl[0]<<
		"X1="<<xcl[1]<<
		"X2="<<xcl[2]<<
		"X3="<<xcl[3]<<
		"X4="<<xcl[4]<<
		"X5="<<xcl[5]<<
		"Chi2R="<<chi2R<<
		"Chi2Z="<<chi2Z<<
		"Chi2RF="<<chi2RF<<                          //chi2 of trackletes on full track
		"Chi2ZF="<<chi2ZF<<                          //chi2 z on tracklets on full track
		"Chi2ZT2="<<chi2ZT2<<                        //chi2 z on tracklets on full track  - rieman tilt
		"Chi2ZTC="<<chi2ZTC<<                        //chi2 z on tracklets on full track  - rieman tilt const
		//
		"Chi2TR="<<chi2TR<<                           //chi2 without vertex constrain
		"Chi2TC="<<chi2TC<<                           //chi2 with    vertex constrain
		"C="<<curv<<                                  // non constrained - no tilt correction
		"DR="<<dR<<                                   // DR parameter          - tilt correction
		"DCA="<<dca<<                                 // DCA                   - tilt correction
		"CR="<<CR<<                                   // non constrained curvature - tilt correction
		"CC="<<CC<<                                   // constrained curvature
		"Polz0="<<polz0c<<
		"Polz1="<<polz1c<<
		"RPolz0="<<rpolz0<<
		"RPolz1="<<rpolz1<<
		"Ncl="<<nclusters<<
		"Nlayers="<<nlayers<<
		"NUsedS="<<nusedCl<<
		"NUsed="<<nusedf<<
		"Findable="<<findable<<
		"Like="<<like<<
		"LikePrim="<<likePrim<<
		"Likechi2C="<<likechi2C<<
		"Likechi2TR="<<likechi2TR<<
		"Likezf="<<likezf<<
		"LikeF="<<seedquality2[registered]<<
		"S0.="<<&cseed[0]<<
		"S1.="<<&cseed[1]<<
		"S2.="<<&cseed[2]<<
		"S3.="<<&cseed[3]<<
		"S4.="<<&cseed[4]<<
		"S5.="<<&cseed[5]<<
		"SB0.="<<&seedb[0]<<
		"SB1.="<<&seedb[1]<<
		"SB2.="<<&seedb[2]<<
		"SB3.="<<&seedb[3]<<
		"SB4.="<<&seedb[4]<<
		"SB5.="<<&seedb[5]<<
		"Label="<<label<<
		"Freq="<<frequency<<
		"sLayer="<<sLayer<<
		"\n";
	    }
	    if (registered<maxseed-1) {
	      registered++;
	      cseed = seed[registered];
	    }
	  }// end of loop over layer 1
	}  // end of loop over layer 0 
      }    // end of loop over layer 3     
    }      // end of loop over seeding time bins 
    //
    // choos best
    //
    TMath::Sort(registered,seedquality2,sort,kTRUE);
    Bool_t signedseed[maxseed];
    for (Int_t i=0;i<registered;i++){
      signedseed[i]= kFALSE;
    }
    for (Int_t iter=0; iter<5; iter++){
      for (Int_t iseed=0;iseed<registered;iseed++){      
	Int_t index = sort[iseed];
	if (signedseed[index]) continue;
	Int_t labelsall[1000];
	Int_t nlabelsall=0;
	Int_t naccepted=0;;
	Int_t sLayer = seedlayer[index];
	Int_t ncl   = 0;
	Int_t nused = 0;
	Int_t nlayers =0;
	Int_t findable   = 0;
	for (Int_t jLayer=0;jLayer<6;jLayer++){
	  if (TMath::Abs(seed[index][jLayer].fYref[0]/xcl[jLayer])<0.15)
	    findable++;
	  if (seed[index][jLayer].isOK()){
	    seed[index][jLayer].UpdateUsed();
	    ncl   +=seed[index][jLayer].fN2;
	    nused +=seed[index][jLayer].fNUsed;
	    nlayers++;
	    //cooking label
	    for (Int_t itime=0;itime<25;itime++){
	      if (seed[index][jLayer].fUsable[itime]){
		naccepted++;
		for (Int_t ilab=0;ilab<3;ilab++){
		  Int_t tindex = seed[index][jLayer].fClusters[itime]->GetLabel(ilab);
		  if (tindex>=0){
		    labelsall[nlabelsall] = tindex;
		    nlabelsall++;
		  }
		}
	      }
	    }
	  }
	}
	//
	if (nused>30) continue;
	//
 	if (iter==0){
	  if (nlayers<6) continue;
	  if (TMath::Log(0.000000001+seedquality2[index])<-5.) continue;   // gold
	}
	//
	if (iter==1){
	  if (nlayers<findable) continue;
	  if (TMath::Log(0.000000001+seedquality2[index])<-4.) continue;  //
	}
	//
	//
	if (iter==2){
	  if (nlayers==findable || nlayers==6) continue;
	  if (TMath::Log(0.000000001+seedquality2[index])<-6.) continue;
	}
	//
	if (iter==3){
	  if (TMath::Log(0.000000001+seedquality2[index])<-5.) continue;
	}
	//
	if (iter==4){
	  if (TMath::Log(0.000000001+seedquality2[index])-nused/(nlayers-3.)<-15.) continue;
	}
	//
	signedseed[index] = kTRUE;
	//
	Int_t labels[1000], outlab[1000];
	Int_t nlab=0;
	for (Int_t iLayer=0;iLayer<6;iLayer++){
	  if (seed[index][iLayer].isOK()){
	    if (seed[index][iLayer].fLabels[0]>=0) {
	      labels[nlab] = seed[index][iLayer].fLabels[0];
	      nlab++;
	    }
	    if (seed[index][iLayer].fLabels[1]>=0) {
	      labels[nlab] = seed[index][iLayer].fLabels[1];
	      nlab++;
	    }	 
	  }     
	}
	Freq(nlab,labels,outlab,kFALSE);
	Int_t label  = outlab[0];
	Int_t frequency  = outlab[1];
	Freq(nlabelsall,labelsall,outlab,kFALSE);
	Int_t label1 = outlab[0];
	Int_t label2 = outlab[2];
	Float_t fakeratio = (naccepted-outlab[1])/Float_t(naccepted);
	Float_t ratio = Float_t(nused)/Float_t(ncl);
	if (ratio<0.25){
	  for (Int_t jLayer=0;jLayer<6;jLayer++){
	    if ( seed[index][jLayer].isOK()&&TMath::Abs(seed[index][jLayer].fYfit[1]-seed[index][jLayer].fYfit[1])<0.2 )
	      seed[index][jLayer].UseClusters();  //sign gold
	  }
	}
	//
	Int_t eventNr = esd->GetEventNumber();
	TTreeSRedirector& cstream = *fDebugStreamer;
	//
	// register seed
	//
	AliTRDtrack * track = RegisterSeed(seed[index],seedparams[index]);
	AliTRDtrack dummy;
	if (!track) track=&dummy;
	else{
	  AliESDtrack esdtrack;
	  esdtrack.UpdateTrackParams(track, AliESDtrack::kTRDout);
	  esdtrack.SetLabel(label);
	  esd->AddTrack(&esdtrack);	
	  TTreeSRedirector& cstream = *fDebugStreamer;
	  cstream<<"Tracks"<<
	    "EventNr="<<eventNr<<
	    "ESD.="<<&esdtrack<<
	    "trd.="<<track<<
	    "trdback.="<<track<<
	    "\n";
	}

	cstream<<"Seeds2"<<
	  "Iter="<<iter<<
	  "Track.="<<track<<
	  "Like="<<seedquality[index]<<
	  "LikeF="<<seedquality2[index]<<
	  "S0.="<<&seed[index][0]<<
	  "S1.="<<&seed[index][1]<<
	  "S2.="<<&seed[index][2]<<
	  "S3.="<<&seed[index][3]<<
	  "S4.="<<&seed[index][4]<<
	  "S5.="<<&seed[index][5]<<
	  "Label="<<label<<
	  "Label1="<<label1<<
	  "Label2="<<label2<<
	  "FakeRatio="<<fakeratio<<
	  "Freq="<<frequency<<
	  "Ncl="<<ncl<<	
	  "Nlayers="<<nlayers<<
	  "Findable="<<findable<<
	  "NUsed="<<nused<<
	  "sLayer="<<sLayer<<
	  "EventNr="<<eventNr<<
	  "\n";
      }
    }
  }        // end of loop over sectors
  delete [] pseed;
}
          
//_____________________________________________________________________________
Int_t AliTRDtracker::ReadClusters(TObjArray *array, TTree *ClusterTree) const
{
  //
  // Reads AliTRDclusters (option >= 0) or AliTRDrecPoints (option < 0) 
  // from the file. The names of the cluster tree and branches 
  // should match the ones used in AliTRDclusterizer::WriteClusters()
  //
  Int_t nsize = Int_t(ClusterTree->GetTotBytes()/(sizeof(AliTRDcluster))); 
  TObjArray *clusterArray = new TObjArray(nsize+1000); 
  
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
//       if (c->GetNPads()>3&&(iCluster%3>0)) {
// 	delete clusterArray->RemoveAt(iCluster);
// 	continue;
//       }
      //      AliTRDcluster *co = new AliTRDcluster(*c);  //remove unnecesary coping - + clusters are together in memory
      AliTRDcluster *co = c;
      co->SetSigmaY2(c->GetSigmaY2() * fSY2corr);
      Int_t ltb = co->GetLocalTimeBin();
      if(ltb == 19) co->SetSigmaZ2(c->GetSigmaZ2());
      else if(fNoTilt) co->SetSigmaZ2(c->GetSigmaZ2() * fSZ2corr);
      array->AddLast(co);
      //      delete clusterArray->RemoveAt(iCluster); 
      clusterArray->RemoveAt(iCluster); 
    }
  }
//   cout<<"Allocated"<<nsize<<"\tLoaded"<<array->GetEntriesFast()<<"\n";

  delete clusterArray;

  return 0;
}

//__________________________________________________________________
Bool_t AliTRDtracker::GetTrackPoint(Int_t index, AliTrackPoint& p) const
{
  //
  // Get track space point with index i
  // Origin: C.Cheshkov
  //

  AliTRDcluster *cl = (AliTRDcluster*)fClusters->UncheckedAt(index);
  Int_t  idet = cl->GetDetector();
  Int_t  isector = fGeom->GetSector(idet);
  Int_t  ichamber= fGeom->GetChamber(idet);
  Int_t  iplan   = fGeom->GetPlane(idet);
  Double_t local[3];
  local[0]=GetX(isector,iplan,cl->GetLocalTimeBin());
  local[1]=cl->GetY();
  local[2]=cl->GetZ();
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
  Int_t modId = isector*fGeom->Ncham()+ichamber;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,modId);
  p.SetVolumeID(volid);

  return kTRUE;

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
  const Float_t kmaxchi2 =18;
  const Float_t kmincl   =10;
  AliTRDtrack * track  = (AliTRDtrack*)t;
  //
  Int_t ncl=t->GetNumberOfClusters();
  for (Int_t i=from; i<ncl; i++) {
    Int_t index = t->GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    //
    Int_t iplane = fGeom->GetPlane(c->GetDetector());
    if (track->fTracklets[iplane].GetChi2()>kmaxchi2) continue; 
    if (track->fTracklets[iplane].GetN()<kmincl) continue; 
    if (!(c->IsUsed())) c->Use();
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
  AliTRDpadPlane *padPlane = 0;

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

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    printf("<AliTRDtracker::AliTRDtrackingSector::AliTRDtrackingSector> ");
    printf("Could not get common params\n");
    return;
  }
    
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
    // Modidified for new pad plane class, 22.04.05 (C.B.)
    // ymaxsensitive = (fPar->GetColPadSize(plane)*fPar->GetColMax(plane)-4)/2.;
    padPlane = commonParam->GetPadPlane(plane,0);
    ymaxsensitive = (padPlane->GetColSize(1)*padPlane->GetNcols()-4)/2.;

    //    ymaxsensitive = (fPar->GetColPadSize(plane)*fPar->GetColMax(plane)-4)/2.;
    
    for(Int_t ch = 0; ch < kNchambers; ch++) {
      zmax[ch] = fGeom->GetChamberLength(plane,ch)/2;
      //
      // Modidified for new pad plane class, 22.04.05 (C.B.)
      //Float_t pad = fPar->GetRowPadSize(plane,ch,0);
      Float_t pad = padPlane->GetRowSize(1);
      //Float_t pad = fPar->GetRowPadSize(plane,ch,0);
      Float_t row0 = commonParam->GetRow0(plane,ch,0);
      Int_t nPads = commonParam->GetRowMax(plane,ch,0);
      zmaxsensitive[ch] = Float_t(nPads)*pad/2.;      
      //      zc[ch] = (pad * nPads)/2 + row0 - pad/2;
      //      zc[ch] = (pad * nPads)/2 + row0;
      zc[ch] = -(pad * nPads)/2 + row0;
      //zc[ch] = row0+zmax[ch]-AliTRDgeometry::RpadW();

    }

    dx  = fgkDriftCorrection*fPar->GetDriftVelocity()
        / AliTRDcalibDB::Instance()->GetSamplingFrequency();
    rho = 0.00295 * 0.85; radLength = 11.0;  

    Double_t x0 = (Double_t) AliTRDgeometry::GetTime0(plane);
    Double_t xbottom = x0 - dxDrift;
    Double_t xtop = x0 + dxAmp;
    //
    // Amplification region
    steps = (Int_t) (dxAmp/dx);

    for(tb = 0; tb < steps; tb++) {
      x = x0 + tb * dx + dx/2+ fgkOffsetX;
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

    dx = fgkDriftCorrection*fPar->GetDriftVelocity()
       / AliTRDcalibDB::Instance()->GetSamplingFrequency();
    steps = (Int_t) (dxDrift/dx)+3;

    for(tb = 0; tb < steps; tb++) {
      x = x0 - tb * dx - dx/2 + fgkOffsetX;                      //temporary fix - fix it the parameters
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
  delete [] zmaxsensitive;

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
  
  Double_t dx = fgkDriftCorrection*(Double_t) fPar->GetDriftVelocity()
                         / AliTRDcalibDB::Instance()->GetSamplingFrequency();

  Int_t tbAmp = fPar->GetTimeBefore();
  Int_t maxAmp = (Int_t) ((dxAmp+0.000001)/dx);
  if(kTRUE) maxAmp = 0;   // intentional until we change parameter class 
  Int_t tbDrift = fPar->GetTimeMax();
  Int_t maxDrift = (Int_t) ((dxDrift+0.000001)/dx)+4; // MI change  - take also last time bins

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
//     if(gap < -0.01) {
//       printf("*** warning: layers %d and %d are overlayed:\n",i,i+1);
//       printf("             %f + %f + %f > %f\n", x1, dx1/2, dx2/2, x2);
//     }
//     if(gap > 0.01) { 
//       printf("*** warning: layers %d and %d have a large gap:\n",i,i+1);
//       printf("             (%f - %f) - (%f + %f) = %f\n", 
//              x2, dx2/2, x1, dx1, gap);
//     }
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



Bool_t AliTRDtracker::AliTRDpropagationLayer::GetPropagationParameters(
        Double_t y, Double_t z, Double_t &dx, Double_t &rho, Double_t &radLength, 
        Bool_t &lookForCluster) const
{
  //
  // Returns radial step <dx>, density <rho>, rad. length <radLength>,
  // and sensitivity <lookForCluster> in point <y,z>  
  //

  Double_t alpha =  AliTRDgeometry::GetAlpha(); 
  Double_t ymax  =  fX*TMath::Tan(0.5*alpha);


  dx  = fdX;
  rho = fRho;
  radLength  = fX0;
  lookForCluster = kFALSE;
  Bool_t cross =kFALSE;
  //
  //
  if ( (ymax-TMath::Abs(y))<3.){   //cross material
    rho*=40.;
    radLength*=40.;
    cross=kTRUE;
  }
  //
  // check dead regions in sensitive volume 
    //
  Int_t zone=-1;
  for(Int_t ch = 0; ch < (Int_t) kZones; ch++) {
    if (TMath::Abs(z - fZc[ch]) > fZmax[ch]) continue;  //not in given zone
    //
    if  (TMath::Abs(z - fZc[ch]) < fZmaxSensitive[ch]){ 
      if (fTimeBinIndex>=0) lookForCluster = !(fIsHole[zone]);
      if(TMath::Abs(y) > fYmaxSensitive){  
	lookForCluster = kFALSE;	
      }
      if (fIsHole[zone]) {
	//if hole
	rho = 1.29e-3;
	radLength = 36.66;
      }
    }else{
      cross = kTRUE; rho = 2.7; radLength = 24.01;  //aluminium in between
    }        
  }
  //
  if (fTimeBinIndex>=0) return cross;
  //
  //
  // check hole
  if (fHole==kFALSE) return cross;
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
  return cross;
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

Int_t AliTRDtracker::AliTRDpropagationLayer::Find(Float_t y) const {

// Returns index of the cluster nearest in Y    

  if (fN<=0) return 0;
  if (y <= fClusters[0]->GetY()) return 0;
  if (y > fClusters[fN-1]->GetY()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (y > fClusters[m]->GetY()) b=m+1;
    else e=m;
  }
  return m;
}    

Int_t AliTRDtracker::AliTRDpropagationLayer::FindNearestCluster(Float_t y, Float_t z, Float_t maxroad, Float_t maxroadz) const 
{
  //
  // Returns index of the cluster nearest to the given y,z
  //
  Int_t index = -1;
  Int_t maxn = fN;
  Float_t mindist = maxroad;			
  //
  for (Int_t i=Find(y-maxroad); i<maxn; i++) {
    AliTRDcluster* c=(AliTRDcluster*)(fClusters[i]);
    Float_t ycl = c->GetY();
    //
    if (ycl > y+maxroad) break;
    if (TMath::Abs(c->GetZ()-z) > maxroadz) continue;      
    if (TMath::Abs(ycl-y)<mindist){
      mindist = TMath::Abs(ycl-y);
      index = fIndex[i];
    }        
  }						
  return index;
}             


//---------------------------------------------------------

Double_t AliTRDtracker::GetTiltFactor(const AliTRDcluster* c) {
//
//  Returns correction factor for tilted pads geometry 
//
  Int_t det = c->GetDetector();    
  Int_t plane = fGeom->GetPlane(det);
  AliTRDpadPlane *padPlane = AliTRDCommonParam::Instance()->GetPadPlane(plane,0);
  Double_t h01 = TMath::Tan(-TMath::Pi() / 180.0 * padPlane->GetTiltingAngle());

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
    // Quality control of TRD track.
    if (nCluster[iPlane]<= 5) {
      clscharge[iPlane]=0.0;
      timebin[iPlane]=-1;
    }
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


Int_t AliTRDtracker::FindClusters(Int_t sector, Int_t t0, Int_t t1, AliTRDtrack * track, Int_t *clusters,AliTRDtracklet&tracklet)
{
  //
  //
  //  try to find nearest clusters to the track in timebins from t0 to t1 
  //  
  //
  //  
  // correction coeficients   - depends on TRD parameters  - to be changed according it
  //

  Double_t x[100],yt[100],zt[100];
  Double_t xmean=0;   //reference x
  Double_t dz[10][100],dy[10][100];
  Float_t zmean[100], nmean[100];
  Int_t    clfound=0;
  Int_t    indexes[10][100];    // indexes of the clusters in the road
  AliTRDcluster *cl[10][100];   // pointers to the clusters in the road
  Int_t    best[10][100];       // index of best matching cluster 
  //
  //

  for (Int_t it=0;it<=t1-t0; it++){
    x[it]=0;
    yt[it]=0;
    zt[it]=0;
    clusters[it+t0]=-2;
    zmean[it]=0;
    nmean[it]=0;
    //
    for (Int_t ih=0;ih<10;ih++){
      indexes[ih][it]=-2;              //reset indexes1
      cl[ih][it]=0;
      dz[ih][it]=-100;
      dy[ih][it]=-100;
      best[ih][it]=0;
    }
  }  
  //
  Double_t x0 = track->GetX();
  Double_t sigmaz = TMath::Sqrt(TMath::Abs(track->GetSigmaZ2()));
  Int_t nall=0;
  Int_t nfound=0;
  Double_t h01 =0;
  Int_t plane =-1;
  Float_t padlength=0;
  AliTRDtrack track2(*track);
  Float_t snpy = track->GetSnp();
  Float_t tany = TMath::Sqrt(snpy*snpy/(1.-snpy*snpy)); 
  if (snpy<0) tany*=-1;
  //
  Double_t sy2=ExpectedSigmaY2(x0,track->GetTgl(),track->GetPt());
  Double_t sz2=ExpectedSigmaZ2(x0,track->GetTgl());
  Double_t road = 15.*sqrt(track->GetSigmaY2() + sy2);
  if (road>6.) road=6.;

  //
  for (Int_t it=0;it<t1-t0;it++){
    Double_t maxChi2[2]={fgkMaxChi2,fgkMaxChi2};      
    AliTRDpropagationLayer& timeBin=*(fTrSec[sector]->GetLayer(it+t0));
    if (timeBin==0) continue;  // no indexes1
    Int_t maxn = timeBin;
    x[it] = timeBin.GetX();
    track2.PropagateTo(x[it]);
    yt[it] = track2.GetY();
    zt[it] = track2.GetZ();
    
    Double_t  y=yt[it],z=zt[it];
    Double_t chi2 =1000000;
    nall++;
    //
    // find 2 nearest cluster at given time bin
    // 
    // 
    for (Int_t i=timeBin.Find(y-road); i<maxn; i++) {
      AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);
      h01 = GetTiltFactor(c);
      if (plane<0){
	Int_t det = c->GetDetector();    
	plane = fGeom->GetPlane(det);
	padlength = TMath::Sqrt(c->GetSigmaZ2()*12.);
      }
      //      if (c->GetLocalTimeBin()==0) continue;
      if (c->GetY() > y+road) break;
      if((c->GetZ()-z)*(c->GetZ()-z) > 12. * sz2) continue;      

      Double_t dist = TMath::Abs(c->GetZ()-z);
      if (dist> (0.5*padlength+6.*sigmaz)) continue;   // 6 sigma boundary cut
      Double_t cost = 0;
      //
      if (dist> (0.5*padlength-sigmaz)){   //  sigma boundary cost function
	cost =  (dist-0.5*padlength)/(2.*sigmaz);
	if (cost>-1) cost= (cost+1.)*(cost+1.);
	else cost=0;
      }      
      //      Int_t label = TMath::Abs(track->GetLabel());
      //      if (c->GetLabel(0)!=label && c->GetLabel(1)!=label&&c->GetLabel(2)!=label) continue;
      chi2=track2.GetPredictedChi2(c,h01)+cost;
      //
      clfound++;      
      if (chi2 > maxChi2[1]) continue;
      
      for (Int_t ih=2;ih<9; ih++){  //store the clusters in the road
	if (cl[ih][it]==0){
	  cl[ih][it] = c;
	  indexes[ih][it] =timeBin.GetIndex(i);   // index - 9 - reserved for outliers
	  break;
	}
      }
      //
      if (chi2 <maxChi2[0]){
	maxChi2[1]     = maxChi2[0];
	maxChi2[0]     = chi2;
	indexes[1][it] = indexes[0][it];
	cl[1][it]      = cl[0][it];
	indexes[0][it] = timeBin.GetIndex(i);
	cl[0][it]      = c;
	continue;
      }
      maxChi2[1]=chi2;
      cl[1][it] = c;
      indexes[1][it] =timeBin.GetIndex(i); 
    }         
    if (cl[0][it]){
      nfound++;
      xmean += x[it];
    }
  }
  //
  if (nfound<4) return 0;  
  xmean /=Float_t(nfound);     // middle x
  track2.PropagateTo(xmean);   // propagate track to the center
  //
  // choose one of the variants
  //
  Int_t changes[10];
  Float_t sumz      = 0;
  Float_t sum       = 0;
  Double_t sumdy    = 0;
  Double_t sumdy2   = 0;
  Double_t sumx     = 0;
  Double_t sumxy    = 0;
  Double_t sumx2    = 0;
  Double_t mpads    = 0;
  //
  Int_t   ngood[10];
  Int_t   nbad[10];
  //
  Double_t meanz[10];
  Double_t moffset[10];    // mean offset
  Double_t mean[10];       // mean value
  Double_t angle[10];      // angle
  //
  Double_t smoffset[10];   // sigma of mean offset
  Double_t smean[10];      // sigma of mean value
  Double_t sangle[10];     // sigma of angle
  Double_t smeanangle[10]; // correlation
  //
  Double_t sigmas[10];     
  Double_t tchi2s[10];      // chi2s for tracklet
  //
  // calculate zmean
  //
  for (Int_t it=0;it<t1-t0;it++){
    if (!cl[0][it]) continue;
    for (Int_t dt=-3;dt<=3;dt++){
      if (it+dt<0) continue;
      if (it+dt>t1-t0) continue;
      if (!cl[0][it+dt]) continue;
      zmean[it]+=cl[0][it+dt]->GetZ();
      nmean[it]+=1.;
    }
    zmean[it]/=nmean[it]; 
  }
  //
  for (Int_t it=0; it<t1-t0;it++){
    best[0][it]=0;
    for (Int_t ih=0;ih<10;ih++){
      dz[ih][it]=-100;
      dy[ih][it]=-100;
      if (!cl[ih][it]) continue;
      //Float_t poscor =  fgkCoef*(cl[ih][it]->GetLocalTimeBin() - fgkMean)+fgkOffset;      
      Float_t poscor =  0;   // applied during loading of clusters      
      if (cl[ih][it]->IsUsed()) poscor=0;  // correction already applied
      dz[ih][it]  = cl[ih][it]->GetZ()- zt[it];                               // calculate distance from track  in z
      dy[ih][it]  = cl[ih][it]->GetY()+ dz[ih][it]*h01 - poscor  -yt[it];     //                                in y
    }
    // minimize changes
    if (!cl[0][it]) continue;
    if (TMath::Abs(cl[0][it]->GetZ()-zmean[it])> padlength*0.8 &&cl[1][it])
      if (TMath::Abs(cl[1][it]->GetZ()-zmean[it])< padlength*0.5){
	best[0][it]=1;
      }
  }
  //
  // iterative choosing of "best path"
  //
  //
  Int_t label = TMath::Abs(track->GetLabel());
  Int_t bestiter=0;
  //
  for (Int_t iter=0;iter<9;iter++){
    //
    changes[iter]= 0;
    sumz      = 0; sum=0; sumdy=0;sumdy2=0;sumx=0;sumx2=0;sumxy=0;mpads=0; ngood[iter]=0; nbad[iter]=0; 
    // linear fit
    for (Int_t it=0;it<t1-t0;it++){
      if (!cl[best[iter][it]][it]) continue;
      //calculates pad-row changes
      Double_t zbefore= cl[best[iter][it]][it]->GetZ();
      Double_t zafter = cl[best[iter][it]][it]->GetZ();
      for (Int_t itd = it-1; itd>=0;itd--) {
	if (cl[best[iter][itd]][itd]) {
	  zbefore= cl[best[iter][itd]][itd]->GetZ();
	  break;
	}
      }
      for (Int_t itd = it+1; itd<t1-t0;itd++) {
	if (cl[best[iter][itd]][itd]) {
	  zafter= cl[best[iter][itd]][itd]->GetZ();
	  break;
	}
      }
      if (TMath::Abs(cl[best[iter][it]][it]->GetZ()-zbefore)>0.1&&TMath::Abs(cl[best[iter][it]][it]->GetZ()-zafter)>0.1) changes[iter]++;
      //
      Double_t dx = x[it]-xmean;  // distance to reference x
      sumz += cl[best[iter][it]][it]->GetZ();      
      sum++;
      sumdy += dy[best[iter][it]][it];
      sumdy2+= dy[best[iter][it]][it]*dy[best[iter][it]][it];
      sumx  += dx;
      sumx2 += dx*dx;
      sumxy  += dx*dy[best[iter][it]][it];
      mpads += cl[best[iter][it]][it]->GetNPads();
      if (cl[best[iter][it]][it]->GetLabel(0)==label || cl[best[iter][it]][it]->GetLabel(1)==label||cl[best[iter][it]][it]->GetLabel(2)==label){
	ngood[iter]++;
      }
      else{
	nbad[iter]++;
      }
    }
    //
    // calculates line parameters
    //
    Double_t det  = sum*sumx2-sumx*sumx;
    angle[iter]   = (sum*sumxy-sumx*sumdy)/det;
    mean[iter]    = (sumx2*sumdy-sumx*sumxy)/det;
    meanz[iter]   = sumz/sum;    
    moffset[iter] = sumdy/sum;
    mpads        /= sum;                         // mean number of pads
    //
    //
    Double_t  sigma2 = 0;   // normalized residuals - for line fit
    Double_t  sigma1 = 0;   // normalized residuals - constant fit
    //
    for (Int_t it=0;it<t1-t0;it++){
      if (!cl[best[iter][it]][it]) continue;
      Double_t dx = x[it]-xmean;
      Double_t ytr = mean[iter]+angle[iter]*dx;
      sigma2 += (dy[best[iter][it]][it]-ytr)*(dy[best[iter][it]][it]-ytr);
      sigma1 +=  (dy[best[iter][it]][it]-moffset[iter])*(dy[best[iter][it]][it]-moffset[iter]);
      sum++;
    }
    sigma2      /=(sum-2);                    // normalized residuals
    sigma1      /=(sum-1);                    // normalized residuals
    //
    smean[iter]       = sigma2*(sumx2/det);   // estimated error2 of mean
    sangle[iter]      = sigma2*(sum/det);     // estimated error2 of angle
    smeanangle[iter]  = sigma2*(-sumx/det);   // correlation
    //
    //
    sigmas[iter]  = TMath::Sqrt(sigma1);      //
    smoffset[iter]= (sigma1/sum)+0.01*0.01;             // sigma of mean offset + unisochronity sigma 
    //
    // iterative choosing of "better path"
    //
    for (Int_t it=0;it<t1-t0;it++){
      if (!cl[best[iter][it]][it]) continue;
      //
      Double_t sigmatr2 = smoffset[iter]+0.5*tany*tany;             //add unisochronity + angular effect contribution
      Double_t sweight  = 1./sigmatr2+1./track->GetSigmaY2();
      Double_t weighty  = (moffset[iter]/sigmatr2)/sweight;         // weighted mean
      Double_t sigmacl  = TMath::Sqrt(sigma1*sigma1+track->GetSigmaY2());   //
      Double_t mindist=100000; 
      Int_t ihbest=0;
      for (Int_t ih=0;ih<10;ih++){
	if (!cl[ih][it]) break;
	Double_t dist2 = (dy[ih][it]-weighty)/sigmacl;
	dist2*=dist2;    //chi2 distance
	if (dist2<mindist){
	  mindist = dist2;
	  ihbest =ih;
	}
      }
      best[iter+1][it]=ihbest;
    }
    //
    //  update best hypothesy if better chi2 according tracklet position and angle
    //
    Double_t sy2 = smean[iter]  + track->GetSigmaY2();
    Double_t sa2 = sangle[iter] + track->fCee;
    Double_t say = track->fCey;
    //    Double_t chi20 = mean[bestiter]*mean[bestiter]/sy2+angle[bestiter]*angle[bestiter]/sa2;
    // Double_t chi21 = mean[iter]*mean[iter]/sy2+angle[iter]*angle[iter]/sa2;

    Double_t detchi    = sy2*sa2-say*say;
    Double_t invers[3] = {sa2/detchi, sy2/detchi, -say/detchi};   //inverse value of covariance matrix  
    
    Double_t chi20 = mean[bestiter]*mean[bestiter]*invers[0]+angle[bestiter]*angle[bestiter]*invers[1]+
      2.*mean[bestiter]*angle[bestiter]*invers[2];
    Double_t chi21 = mean[iter]*mean[iter]*invers[0]+angle[iter]*angle[iter]*invers[1]+
      2*mean[iter]*angle[iter]*invers[2];
    tchi2s[iter] =chi21;
    //
    if (changes[iter]<=changes[bestiter] && chi21<chi20) {
      bestiter =iter;      
    }
  }
  //
  //set clusters 
  //
  Double_t sigma2 = sigmas[0];   // choose as sigma  from 0 iteration
  Short_t maxpos    = -1;
  Float_t maxcharge =  0;
  Short_t maxpos4    = -1;
  Float_t maxcharge4 =  0;
  Short_t maxpos5    = -1;
  Float_t maxcharge5 =  0;

  //if (tchi2s[bestiter]>25.) sigma2*=tchi2s[bestiter]/25.;
  //if (tchi2s[bestiter]>25.) sigma2=1000.;  // dont'accept

  Double_t expectederr = sigma2*sigma2+0.01*0.01;
  if (mpads>3.5) expectederr  +=   (mpads-3.5)*0.04;
  if (changes[bestiter]>1) expectederr+=   changes[bestiter]*0.01; 
  expectederr+=(0.03*(tany-fgkExB)*(tany-fgkExB))*15;
  //  if (tchi2s[bestiter]>18.) expectederr*= tchi2s[bestiter]/18.;
  //expectederr+=10000;
  for (Int_t it=0;it<t1-t0;it++){
    if (!cl[best[bestiter][it]][it]) continue;
    //    Float_t poscor =  fgkCoef*(cl[best[bestiter][it]][it]->GetLocalTimeBin() - fgkMean)+fgkOffset;          
    Float_t poscor =  0;           //applied during loading of cluster
    cl[best[bestiter][it]][it]->SetSigmaY2(expectederr);  // set cluster error
    if (!cl[best[bestiter][it]][it]->IsUsed()){
      cl[best[bestiter][it]][it]->SetY( cl[best[bestiter][it]][it]->GetY()-poscor);  // ExB corrction correction    
      //      cl[best[bestiter][it]][it]->Use();
    }
    //
    //  time bins with maximal charge
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ())> maxcharge){
      maxcharge = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
      maxpos = cl[best[bestiter][it]][it]->GetLocalTimeBin();
    }
    
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ())> maxcharge4){
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin()>=4){
	maxcharge4 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos4 = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ())> maxcharge5){
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin()>=5){
	maxcharge5 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos5 = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }
    //
    //  time bins with maximal charge
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ())> maxcharge){
      maxcharge = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
      maxpos = cl[best[bestiter][it]][it]->GetLocalTimeBin();
    }
    
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ())> maxcharge4){
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin()>=4){
	maxcharge4 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos4 = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }
    if (TMath::Abs(cl[best[bestiter][it]][it]->GetQ())> maxcharge5){
      if (cl[best[bestiter][it]][it]->GetLocalTimeBin()>=5){
	maxcharge5 = TMath::Abs(cl[best[bestiter][it]][it]->GetQ());
	maxpos5 = cl[best[bestiter][it]][it]->GetLocalTimeBin();
      }
    }
    clusters[it+t0] = indexes[best[bestiter][it]][it];    
    //if (cl[best[bestiter][it]][it]->GetLocalTimeBin()>4 && cl[best[bestiter][it]][it]->GetLocalTimeBin()<18) clusters[it+t0] = indexes[best[bestiter][it]][it];    //Test
  } 
  //
  // set tracklet parameters
  //
  Double_t trackleterr2 = smoffset[bestiter]+0.01*0.01;
  if (mpads>3.5) trackleterr2  +=   (mpads-3.5)*0.04;
  trackleterr2+=   changes[bestiter]*0.01;
  trackleterr2*=   TMath::Max(14.-nfound,1.);
  trackleterr2+=   0.2*(tany-fgkExB)*(tany-fgkExB); 
  //
  tracklet.Set(xmean, track2.GetY()+moffset[bestiter], meanz[bestiter], track2.GetAlpha(), trackleterr2);  //set tracklet parameters
  tracklet.SetTilt(h01);
  tracklet.SetP0(mean[bestiter]);
  tracklet.SetP1(angle[bestiter]);
  tracklet.SetN(nfound);
  tracklet.SetNCross(changes[bestiter]);
  tracklet.SetPlane(plane);
  tracklet.SetSigma2(expectederr);
  tracklet.SetChi2(tchi2s[bestiter]);
  tracklet.SetMaxPos(maxpos,maxpos4,maxpos5);
  track->fTracklets[plane] = tracklet;
  track->fNWrong+=nbad[0];
  //
  // Debuging part
  //
  TClonesArray array0("AliTRDcluster");
  TClonesArray array1("AliTRDcluster");
  array0.ExpandCreateFast(t1-t0+1);
  array1.ExpandCreateFast(t1-t0+1);
  TTreeSRedirector& cstream = *fDebugStreamer;
  AliTRDcluster dummy;
  Double_t dy0[100];
  Double_t dyb[100]; 

  for (Int_t it=0;it<t1-t0;it++){
    dy0[it] = dy[0][it];
    dyb[it] = dy[best[bestiter][it]][it];
    if(cl[0][it]) {
      new(array0[it]) AliTRDcluster(*cl[0][it]);
    }
    else{
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
  //
  //
  cstream<<"tracklet"<<
    "track.="<<track<<                                       // track parameters
    "tany="<<tany<<                                          // tangent of the local track angle 
    "xmean="<<xmean<<                                        // xmean - reference x of tracklet  
    "tilt="<<h01<<                                           // tilt angle
    "nall="<<nall<<                                          // number of foundable clusters 
    "nfound="<<nfound<<                                      // number of found clusters
    "clfound="<<clfound<<                                    // total number of found clusters in road 
    "mpads="<<mpads<<                                        // mean number of pads per cluster
    "plane="<<plane<<                                        // plane number 
    "road="<<road<<                                          // the width of the used road
    "graph0.="<<&graph0<<                                    // x - y = dy for closest cluster
    "graph1.="<<&graph1<<                                    // x - y = dy for second closest cluster    
    "graphy.="<<&graphy<<                                    // y position of the track
    "graphz.="<<&graphz<<                                    // z position of the track
    //    "fCl.="<<&array0<<                                       // closest cluster
    // "fCl2.="<<&array1<<                                      // second closest cluster
    "maxpos="<<maxpos<<                                      // maximal charge postion
    "maxcharge="<<maxcharge<<                                // maximal charge 
    "maxpos4="<<maxpos4<<                                    // maximal charge postion - after bin 4
    "maxcharge4="<<maxcharge4<<                              // maximal charge         - after bin 4
    "maxpos5="<<maxpos5<<                                    // maximal charge postion - after bin 5
    "maxcharge5="<<maxcharge5<<                              // maximal charge         - after bin 5
    //
    "bestiter="<<bestiter<<                                  // best iteration number 
    "tracklet.="<<&tracklet<<                                // corrspond to the best iteration
    "tchi20="<<tchi2s[0]<<                                   // chi2 of cluster in the 0 iteration
    "tchi2b="<<tchi2s[bestiter]<<                            // chi2 of cluster in the best  iteration
    "sigmas0="<<sigmas[0]<<                                  // residuals sigma 
    "sigmasb="<<sigmas[bestiter]<<                           // residulas sigma
    //
    "ngood0="<<ngood[0]<<                                    // number of good clusters in 0 iteration
    "nbad0="<<nbad[0]<<                                      // number of bad clusters in 0 iteration
    "ngoodb="<<ngood[bestiter]<<                             //                        in  best iteration    
    "nbadb="<<nbad[bestiter]<<                               //                        in  best iteration
    //
    "changes0="<<changes[0]<<                                // changes of pardrows in iteration number 0 
    "changesb="<<changes[bestiter]<<                         // changes of pardrows in best iteration
    //
    "moffset0="<<moffset[0]<<                                // offset fixing angle in iter=0
    "smoffset0="<<smoffset[0]<<                              // sigma of offset fixing angle in iter=0
    "moffsetb="<<moffset[bestiter]<<                         // offset fixing angle in iter=best
    "smoffsetb="<<smoffset[bestiter]<<                       // sigma of offset fixing angle in iter=best
    //
    "mean0="<<mean[0]<<                                      // mean dy in iter=0;
    "smean0="<<smean[0]<<                                    // sigma of mean dy in iter=0
    "meanb="<<mean[bestiter]<<                               // mean dy in iter=best
    "smeanb="<<smean[bestiter]<<                             // sigma of mean dy in iter=best
    //
    "angle0="<<angle[0]<<                                    // angle deviation in the iteration number 0 
    "sangle0="<<sangle[0]<<                                  // sigma of angular deviation in iteration number 0
    "angleb="<<angle[bestiter]<<                             // angle deviation in the best iteration   
    "sangleb="<<sangle[bestiter]<<                           // sigma of angle deviation in the best iteration   
    //
    "expectederr="<<expectederr<<                            // expected error of cluster position
    "\n";
  //
  //
  return nfound;
}


Int_t  AliTRDtracker::Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down)
{    
  //
  //  Sort eleements according occurancy 
  //  The size of output array has is 2*n 
  //
  Int_t * sindexS = new Int_t[n];     // temp array for sorting
  Int_t * sindexF = new Int_t[2*n];   
  for (Int_t i=0;i<n;i++) sindexF[i]=0;
  //
  TMath::Sort(n,inlist, sindexS, down);  
  Int_t last      = inlist[sindexS[0]];
  Int_t val       = last;
  sindexF[0]      = 1;
  sindexF[0+n]    = last;
  Int_t countPos  = 0;
  //
  //  find frequency
  for(Int_t i=1;i<n; i++){
    val = inlist[sindexS[i]];
    if (last == val)   sindexF[countPos]++;
    else{      
      countPos++;
      sindexF[countPos+n] = val;
      sindexF[countPos]++;
      last =val;
    }
  }
  if (last==val) countPos++;
  // sort according frequency
  TMath::Sort(countPos, sindexF, sindexS, kTRUE);
  for (Int_t i=0;i<countPos;i++){
    outlist[2*i  ] = sindexF[sindexS[i]+n];
    outlist[2*i+1] = sindexF[sindexS[i]];
  }
  delete [] sindexS;
  delete [] sindexF;
  
  return countPos;
}

AliTRDtrack * AliTRDtracker::RegisterSeed(AliTRDseed * seeds, Double_t * params)
{
  //
  //
  //
  Double_t alpha=AliTRDgeometry::GetAlpha();
  Double_t shift=AliTRDgeometry::GetAlpha()/2.;
  Double_t c[15];
  c[0] = 0.2;
  c[1] = 0  ; c[2] = 2;
  c[3] = 0  ; c[4] = 0; c[5] = 0.02;
  c[6] = 0  ; c[7] = 0; c[8] = 0;      c[9] = 0.1;
  c[10] = 0  ; c[11] = 0; c[12] = 0;   c[13] = 0.0; c[14] = params[5]*params[5]*0.01;
  //
  Int_t index =0;
  AliTRDcluster *cl =0;
  for (Int_t ilayer=0;ilayer<6;ilayer++){
    if (seeds[ilayer].isOK()){
      for (Int_t itime=22;itime>0;itime--){
	if (seeds[ilayer].fIndexes[itime]>0){
	  index = seeds[ilayer].fIndexes[itime];
	  cl = seeds[ilayer].fClusters[itime];
	  break;
	}
      }
    }
    if (index>0) break;
  }
  if (cl==0) return 0;
  AliTRDtrack * track  = new AliTRDtrack(cl,index,&params[1],c, params[0],params[6]*alpha+shift);
  track->PropagateTo(params[0]-5.);
  track->ResetCovariance(1);
  //
  Int_t rc=FollowBackProlongationG(*track);
  if (rc<30) {
    delete track;
    track =0;
  }else{
    track->CookdEdx();
    CookdEdxTimBin(*track);
    CookLabel(track, 0.9);
  }
  return track;
}






AliTRDseed::AliTRDseed()
{
  //
  //  
  fTilt =0;         // tilting angle
  fPadLength = 0;   // pad length
  fX0 = 0;           // x0 position
  for (Int_t i=0;i<25;i++){
    fX[i]=0;        // !x position
    fY[i]=0;        // !y position
    fZ[i]=0;        // !z position
    fIndexes[i]=0;  // !indexes
    fClusters[i]=0; // !clusters
  }
  for (Int_t i=0;i<2;i++){
    fYref[i]=0;      // reference y
    fZref[i]=0;      // reference z
    fYfit[i]=0;      // y fit position +derivation
    fYfitR[i]=0;      // y fit position +derivation
    fZfit[i]=0;      // z fit position
    fZfitR[i]=0;      // z fit position
    fLabels[i]=0;    // labels
  }
  fSigmaY  = 0;       
  fSigmaY2 = 0;       
  fMeanz=0;         // mean vaue of z
  fZProb=0;         // max probbable z
  fMPads=0;
  //
  fN=0;            // number of associated clusters
  fN2=0;            // number of not crossed
  fNUsed=0;        // number of used clusters
  fNChange=0;      // change z counter
}

void AliTRDseed::Reset(){
  //
  // reset seed
  //
  for (Int_t i=0;i<25;i++){
    fX[i]=0;        // !x position
    fY[i]=0;        // !y position
    fZ[i]=0;        // !z position
    fIndexes[i]=0;  // !indexes
    fClusters[i]=0; // !clusters
    fUsable[i]  = kFALSE;    
  }
  for (Int_t i=0;i<2;i++){
    fYref[i]=0;      // reference y
    fZref[i]=0;      // reference z
    fYfit[i]=0;      // y fit position +derivation
    fYfitR[i]=0;      // y fit position +derivation
    fZfit[i]=0;      // z fit position
    fZfitR[i]=0;      // z fit position
    fLabels[i]=-1;    // labels
  }
  fSigmaY =0;         //"robust" sigma in y
  fSigmaY2=0;         //"robust" sigma in y
  fMeanz =0;         // mean vaue of z
  fZProb =0;         // max probbable z
  fMPads =0;
  //
  fN=0;            // number of associated clusters
  fN2=0;            // number of not crossed
  fNUsed=0;        // number of used clusters
  fNChange=0;      // change z counter
}

void AliTRDseed::CookLabels(){
  //
  // cook 2 labels for seed
  //
  Int_t labels[200];
  Int_t out[200];
  Int_t nlab =0;
  for (Int_t i=0;i<25;i++){
    if (!fClusters[i]) continue;
    for (Int_t ilab=0;ilab<3;ilab++){
      if (fClusters[i]->GetLabel(ilab)>=0){
	labels[nlab] = fClusters[i]->GetLabel(ilab);
	nlab++;
      }
    }
  }
  Int_t nlab2 = AliTRDtracker::Freq(nlab,labels,out,kTRUE);
  fLabels[0] = out[0];
  if (nlab2>1 && out[3]>1) fLabels[1] =out[2];
}

void   AliTRDseed::UseClusters()
{
  //
  // use clusters
  //
   for (Int_t i=0;i<25;i++){
     if (!fClusters[i]) continue;
     if (!(fClusters[i]->IsUsed())) fClusters[i]->Use();
   }
}


void        AliTRDseed::Update(){
  //
  //
  //
  const Float_t ratio = 0.8;
  const Int_t   kClmin        = 6;
  const Float_t kmaxtan  = 2;
  if (TMath::Abs(fYref[1])>kmaxtan) return;             // too much inclined track
  //
  Float_t  sigmaexp = 0.05+TMath::Abs(fYref[1]*0.25);   // expected r.m.s in y direction
  Float_t  ycrosscor = fPadLength*fTilt*0.5;             // y correction for crossing 
  fNChange =0;
  //
  Double_t sumw, sumwx,sumwx2;
  Double_t sumwy, sumwxy, sumwz,sumwxz;
  Int_t    zints[25];        // histograming of the z coordinate - get 1 and second max probable coodinates in z
  Int_t    zouts[50];        //
  Float_t  allowedz[25];     // allowed z for given time bin
  Float_t  yres[25];         // residuals from reference
  Float_t  anglecor = fTilt*fZref[1];  //correction to the angle
  //
  //
  fN=0; fN2 =0;
  for (Int_t i=0;i<25;i++){
    yres[i] =10000;
    if (!fClusters[i]) continue;
    yres[i] = fY[i]-fYref[0]-(fYref[1]+anglecor)*fX[i];   // residual y
    zints[fN] = Int_t(fZ[i]);
    fN++;    
  }
  if (fN<kClmin) return;
  Int_t nz = AliTRDtracker::Freq(fN,zints,zouts,kFALSE);
  fZProb   = zouts[0];
  if (nz<=1) zouts[3]=0;
  if (zouts[1]+zouts[3]<kClmin) return;
  //
  if (TMath::Abs(zouts[0]-zouts[2])>12.) zouts[3]=0;   // z distance bigger than pad - length
  //
  Int_t  breaktime = -1;
  Bool_t mbefore   = kFALSE;
  Int_t  cumul[25][2];
  Int_t  counts[2]={0,0};
  //
  if (zouts[3]>=3){
    //
    // find the break time allowing one chage on pad-rows with maximal numebr of accepted clusters
    //
    fNChange=1;
    for (Int_t i=0;i<25;i++){
      cumul[i][0] = counts[0];
      cumul[i][1] = counts[1];
      if (TMath::Abs(fZ[i]-zouts[0])<2) counts[0]++;
      if (TMath::Abs(fZ[i]-zouts[2])<2) counts[1]++;
    }
    Int_t  maxcount  = 0;
    for (Int_t i=0;i<24;i++) {
      Int_t after  = cumul[24][0]-cumul[i][0];
      Int_t before = cumul[i][1];
      if (after+before>maxcount) { 
	maxcount=after+before; 
	breaktime=i;
	mbefore=kFALSE;
      }
      after  = cumul[24][1]-cumul[i][1];
      before = cumul[i][0];
      if (after+before>maxcount) { 
	maxcount=after+before; 
	breaktime=i;
	mbefore=kTRUE;
      }
    }
    breaktime-=1;
  }
  for (Int_t i=0;i<25;i++){
    if (i>breaktime)  allowedz[i] =   mbefore  ? zouts[2]:zouts[0];
    if (i<=breaktime) allowedz[i] = (!mbefore) ? zouts[2]:zouts[0];
  }  
  if ( (allowedz[0]>allowedz[24] && fZref[1]<0) || (allowedz[0]<allowedz[24] &&  fZref[1]>0)){
    //
    // tracklet z-direction not in correspondance with track z direction 
    //
    fNChange =0;
    for (Int_t i=0;i<25;i++){
      allowedz[i] =  zouts[0];  //only longest taken
    } 
  }
  //
  if (fNChange>0){
    //
    // cross pad -row tracklet  - take the step change into account
    //
    for (Int_t i=0;i<25;i++){
      if (!fClusters[i]) continue; 
      if (TMath::Abs(fZ[i]-allowedz[i])>2) continue;
      yres[i] = fY[i]-fYref[0]-(fYref[1]+anglecor)*fX[i];   // residual y
      if (TMath::Abs(fZ[i]-fZProb)>2){
	if (fZ[i]>fZProb) yres[i]+=fTilt*fPadLength;
	if (fZ[i]<fZProb) yres[i]-=fTilt*fPadLength;
      }
    }
  }
  //
  Double_t yres2[25];
  Double_t mean,sigma;
  for (Int_t i=0;i<25;i++){
    if (!fClusters[i]) continue;
    if (TMath::Abs(fZ[i]-allowedz[i])>2) continue;
    yres2[fN2] =  yres[i];
    fN2++;
  }
  if (fN2<kClmin){
    fN2 = 0;
    return;
  }
  EvaluateUni(fN2,yres2,mean,sigma,Int_t(fN2*ratio-2));
  if (sigma<sigmaexp*0.8) sigma=sigmaexp;
  fSigmaY = sigma;
  //
  //
  // reset sums
  sumw=0; sumwx=0; sumwx2=0;
  sumwy=0; sumwxy=0; sumwz=0;sumwxz=0;
  fN2 =0;
  fMeanz =0;
  fMPads =0;
  //
  for (Int_t i=0;i<25;i++){
    fUsable[i]=kFALSE;
    if (!fClusters[i]) continue;
    if (TMath::Abs(fZ[i]-allowedz[i])>2)  continue;
    if (TMath::Abs(yres[i]-mean)>4.*sigma) continue;
    fUsable[i] = kTRUE;
    fN2++;
    fMPads+=fClusters[i]->GetNPads();
    Float_t weight =1;
    if (fClusters[i]->GetNPads()>4) weight=0.5;
    if (fClusters[i]->GetNPads()>5) weight=0.2;
    //
    Double_t x = fX[i];
    sumw+=weight; sumwx+=x*weight; sumwx2+=x*x*weight;
    sumwy+=weight*yres[i];  sumwxy+=weight*(yres[i])*x;
    sumwz+=weight*fZ[i];    sumwxz+=weight*fZ[i]*x;
  }
  if (fN2<kClmin){
    fN2 = 0;
    return;
  }
  fMeanz       = sumwz/sumw;
  Float_t correction =0;
  if (fNChange>0){
    // tracklet on boundary
    if (fMeanz<fZProb) correction =   ycrosscor;
    if (fMeanz>fZProb) correction =  -ycrosscor;
  }
  Double_t det = sumw*sumwx2-sumwx*sumwx;
  fYfitR[0]    = (sumwx2*sumwy-sumwx*sumwxy)/det;
  fYfitR[1]    = (sumw*sumwxy-sumwx*sumwy)/det;
  //
  fSigmaY2     =0;
  for (Int_t i=0;i<25;i++){    
    if (!fUsable[i]) continue;
    Float_t delta = yres[i]-fYfitR[0]-fYfitR[1]*fX[i];
    fSigmaY2+=delta*delta;
  }
  fSigmaY2 = TMath::Sqrt(fSigmaY2/Float_t(fN2-2));
  //
  fZfitR[0]    = (sumwx2*sumwz-sumwx*sumwxz)/det;
  fZfitR[1]    = (sumw*sumwxz-sumwx*sumwz)/det;
  fZfit[0]     = (sumwx2*sumwz-sumwx*sumwxz)/det;
  fZfit[1]     = (sumw*sumwxz-sumwx*sumwz)/det;
  fYfitR[0]   += fYref[0]+correction;
  fYfitR[1]   += fYref[1];
  fYfit[0]     = fYfitR[0];
  fYfit[1]     = fYfitR[1];
  //
  //  
  UpdateUsed();
}






void AliTRDseed::UpdateUsed(){
  //
  fNUsed =0;
  for (Int_t i=0;i<25;i++){
     if (!fClusters[i]) continue;
     if ((fClusters[i]->IsUsed())) fNUsed++;
  }
}


void AliTRDseed::EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh)
{
  //
  // robust estimator in 1D case MI version
  //
  //for the univariate case
  //estimates of location and scatter are returned in mean and sigma parameters
  //the algorithm works on the same principle as in multivariate case -
  //it finds a subset of size hh with smallest sigma, and then returns mean and
  //sigma of this subset

  if (hh==0)
    hh=(nvectors+2)/2;
  Double_t faclts[]={2.6477,2.5092,2.3826,2.2662,2.1587,2.0589,1.9660,1.879,1.7973,1.7203,1.6473};
  Int_t *index=new Int_t[nvectors];
  TMath::Sort(nvectors, data, index, kFALSE);
  //
  Int_t    nquant = TMath::Min(Int_t(Double_t(((hh*1./nvectors)-0.5)*40))+1, 11);
  Double_t factor = faclts[nquant-1];
  //
  //
  Double_t sumx  =0;
  Double_t sumx2 =0;
  Int_t    bestindex = -1;
  Double_t bestmean  = 0; 
  Double_t bestsigma = data[index[nvectors-1]]-data[index[0]];   // maximal possible sigma
  for (Int_t i=0; i<hh; i++){
    sumx  += data[index[i]];
    sumx2 += data[index[i]]*data[index[i]];
  }
  //
  Double_t norm = 1./Double_t(hh);
  Double_t norm2 = 1./Double_t(hh-1);
  for (Int_t i=hh; i<nvectors; i++){
    Double_t cmean  = sumx*norm;
    Double_t csigma = (sumx2 - hh*cmean*cmean)*norm2;
    if (csigma<bestsigma){
      bestmean  = cmean;
      bestsigma = csigma;
      bestindex = i-hh;
    }
    //
    //
    sumx  += data[index[i]]-data[index[i-hh]];
    sumx2 += data[index[i]]*data[index[i]]-data[index[i-hh]]*data[index[i-hh]];
  }
  
  Double_t bstd=factor*TMath::Sqrt(TMath::Abs(bestsigma));
  mean  = bestmean;
  sigma = bstd;
  delete [] index;
}


Float_t   AliTRDseed::FitRiemanTilt(AliTRDseed * cseed, Bool_t terror){
  //
  //
  //
  TLinearFitter fitterT2(4,"hyp4");  // fitting with tilting pads - kz not fixed
  fitterT2.StoreData(kTRUE);
  Float_t xref2 = (cseed[2].fX0+cseed[3].fX0)*0.5; // reference x0 for z
  //
  Int_t npointsT =0;
  fitterT2.ClearPoints();
  for (Int_t iLayer=0; iLayer<6;iLayer++){
    if (!cseed[iLayer].isOK()) continue;
    Double_t tilt = cseed[iLayer].fTilt;

    for (Int_t itime=0;itime<25;itime++){
      if (!cseed[iLayer].fUsable[itime]) continue;
      Double_t x   = cseed[iLayer].fX[itime]+cseed[iLayer].fX0-xref2;  // x relative to the midle chamber
      Double_t y   = cseed[iLayer].fY[itime];
      Double_t z   = cseed[iLayer].fZ[itime];
      // tilted rieman
      //
      Double_t uvt[6];
      Double_t x2 = cseed[iLayer].fX[itime]+cseed[iLayer].fX0;      // global x
      Double_t t = 1./(x2*x2+y*y);
      uvt[1]  = t;    // t
      uvt[0]  = 2.*x2*uvt[1];      // u 
      uvt[2]  = 2.0*tilt*uvt[1];
      uvt[3]  = 2.0*tilt*x*uvt[1];	      
      uvt[4]  = 2.0*(y+tilt*z)*uvt[1];
      //
      Double_t error = 2*uvt[1];
      if (terror) error*=cseed[iLayer].fSigmaY;
      else {error *=0.2;} //default error
      fitterT2.AddPoint(uvt,uvt[4],error);
      npointsT++;
    }
  }
  fitterT2.Eval();
  Double_t rpolz0 = fitterT2.GetParameter(3);
  Double_t rpolz1 = fitterT2.GetParameter(4);	    
  //
  // linear fitter  - not possible to make boundaries
  // non accept non possible z and dzdx combination
  // 	    
  Bool_t   acceptablez =kTRUE;
  for (Int_t iLayer=0; iLayer<6;iLayer++){
    if (cseed[iLayer].isOK()){
      Double_t zT2 =  rpolz0+rpolz1*(cseed[iLayer].fX0 - xref2);
      if (TMath::Abs(cseed[iLayer].fZProb-zT2)>cseed[iLayer].fPadLength*0.5+1)
	acceptablez = kFALSE;
    }
  }
  if (!acceptablez){
    Double_t zmf  = cseed[2].fZref[0]+cseed[2].fZref[1]*(xref2-cseed[2].fX0);
    Double_t dzmf = (cseed[2].fZref[1]+ cseed[3].fZref[1])*0.5;
    fitterT2.FixParameter(3,zmf);
    fitterT2.FixParameter(4,dzmf);
    fitterT2.Eval();
    fitterT2.ReleaseParameter(3);
    fitterT2.ReleaseParameter(4);
    rpolz0 = fitterT2.GetParameter(3);
    rpolz1 = fitterT2.GetParameter(4);
  }
  //
  Double_t chi2TR = fitterT2.GetChisquare()/Float_t(npointsT);  
  Double_t params[3];
  params[0]     =  fitterT2.GetParameter(0);
  params[1]     =  fitterT2.GetParameter(1);
  params[2]     =  fitterT2.GetParameter(2);	    
  Double_t CR     =  1+params[1]*params[1]-params[2]*params[0];
  for (Int_t iLayer = 0; iLayer<6;iLayer++){
    Double_t  x = cseed[iLayer].fX0;
    Double_t  y=0,dy=0, z=0, dz=0;
    // y
    Double_t res2 = (x*params[0]+params[1]);
    res2*=res2;
    res2 = 1.-params[2]*params[0]+params[1]*params[1]-res2;
    if (res2>=0){
      res2 = TMath::Sqrt(res2);
      y    = (1-res2)/params[0];
    }
    //dy
    Double_t x0 = -params[1]/params[0];
    if (-params[2]*params[0]+params[1]*params[1]+1>0){
      Double_t Rm1  = params[0]/TMath::Sqrt(-params[2]*params[0]+params[1]*params[1]+1); 
      if ( 1./(Rm1*Rm1)-(x-x0)*(x-x0)>0){
	Double_t res = (x-x0)/TMath::Sqrt(1./(Rm1*Rm1)-(x-x0)*(x-x0));
	if (params[0]<0) res*=-1.;
	dy = res;
      }
    }
    z  = rpolz0+rpolz1*(x-xref2);
    dz = rpolz1;
    cseed[iLayer].fYref[0] = y;
    cseed[iLayer].fYref[1] = dy;
    cseed[iLayer].fZref[0] = z;
    cseed[iLayer].fZref[1] = dz;
    cseed[iLayer].fC  = CR;
    //
  }
  return chi2TR;
}
