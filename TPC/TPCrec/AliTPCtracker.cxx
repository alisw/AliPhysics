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


//-------------------------------------------------------
//          Implementation of the TPC tracker
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
// 
//  AliTPC parallel tracker
//
//  The track fitting is based on Kalman filtering approach
//  The track finding steps:
//      1. Seeding - with and without vertex constraint
//                 - seeding with vertex constain done at first n^2 proble
//                 - seeding without vertex constraint n^3 problem
//      2. Tracking - follow prolongation road - find cluster - update kalman track
//  The seeding and tracking is repeated several times, in different seeding region.
//  This approach enables to find the track which cannot be seeded in some region of TPC
//  This can happen because of low momenta (track do not reach outer radius), or track is currently in the ded region between sectors, or the track is for the moment overlapped with other track (seed quality is poor) ...

//  With this approach we reach almost 100 % efficiency also for high occupancy events.
//  (If the seeding efficiency in a region is about 90 % than with logical or of several 
//  regions we will reach 100% (in theory - supposing independence) 

//  Repeating several seeding - tracking procedures some of the tracks can be find 
//  several times. 

//  The procedures to remove multi find tacks are impremented:
//  RemoveUsed2 - fast procedure n problem - 
//                 Algorithm - Sorting tracks according quality
//                             remove tracks with some shared fraction 
//                             Sharing in respect to all tacks 
//                             Signing clusters in gold region
//  FindSplitted - slower algorithm n^2
//                 Sort the tracks according quality
//                 Loop over pair of tracks
//                 If overlap with other track bigger than threshold - remove track
//  
//  FindCurling  - Finds the pair of tracks which are curling
//               - About 10% of tracks can be find with this procedure
//                 The combinatorial background is too big to be used in High 
//                  multiplicity environment 
//               - n^2 problem - Slow procedure - currently it is disabled because of 
//                  low efficiency
//                 
//  The number of splitted tracks can be reduced disabling the sharing of the cluster.
//  tpcRecoParam-> SetClusterSharing(kFALSE);
//  IT IS HIGHLY non recomended to use it in high flux enviroonment
//  Even using this switch some tracks can be found more than once 
//  (because of multiple seeding and low quality tracks which will not cross full chamber)
//                          
//
// The tracker itself can be debugged  - the information about tracks can be stored in several // phases of the reconstruction
// To enable storage of the TPC tracks in the ESD friend track
// use AliTPCReconstructor::SetStreamLevel(n); 
//
// The debug level -  different procedure produce tree for numerical debugging of code and data (see comments foEStreamFlags in AliTPCtracker.h  )
//

//
// Adding systematic errors to the covariance:
// 
// The systematic errors due to the misalignment and miscalibration are added to the covariance matrix
// of the tracks (not to the clusters as they are dependent):
// The parameters form AliTPCRecoParam are used AliTPCRecoParam::GetSystematicError
// The systematic errors are expressed there in RMS - position (cm), angle (rad), curvature (1/GeV)
// The default values are 0. 
//
// The systematic errors are added to the covariance matrix in following places:
//
// 1. During fisrt itteration - AliTPCtracker::FillESD
// 2. Second iteration - 
//      2.a ITS->TPC   - AliTPCtracker::ReadSeeds 
//      2.b TPC->TRD   - AliTPCtracker::PropagateBack
// 3. Third iteration  -
//      3.a TRD->TPC   - AliTPCtracker::ReadSeeds
//      3.b TPC->ITS   - AliTPCtracker::RefitInward
//

/* $Id$ */

#include "Riostream.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TTimeStamp.h>
#include "AliLog.h"
#include "AliComplexCluster.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliKink.h"
#include "AliV0.h"
#include "AliVParticle.h"
#include "AliHelix.h"
#include "AliRunLoader.h"
#include "AliTPCClustersRow.h"
#include "AliTPCParam.h"
#include "AliTPCReconstructor.h"
#include "AliTPCpolyTrack.h"
#include "AliTPCreco.h"
//#include "AliTPCseed.h"

#include "AliTPCtrackerSector.h" 
#include "AliTPCtracker.h"
#include "TStopwatch.h"
#include "AliTPCReconstructor.h"
#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "TRandom.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "AliTPCTransform.h"
#include "AliTPCClusterParam.h"
#include "AliTPCdEdxInfo.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"
#include "AliDAQ.h"
#include "AliCosmicTracker.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include <math.h>
//
#include "AliESDfriendTrack.h"

using std::cout;
using std::cerr;
using std::endl;
ClassImp(AliTPCtracker)

//__________________________________________________________________
AliTPCtracker::AliTPCtracker()
  :AliTracker(),
  fkNIS(0),
  fInnerSec(0),
  fkNOS(0),
  fOuterSec(0),
  fN(0),
  fSectors(0),
  fInput(0),
  fOutput(0),
  fSeedTree(0),
  fTreeDebug(0),
  fEvent(0),
  fEventHLT(0),
  fDebug(0),
  fNewIO(kFALSE),
  fNtracks(0),
  fSeeds(0),
  fIteration(0),
  fkParam(0),
  fDebugStreamer(0),
  fUseHLTClusters(4),
  fClExtraRoadY(0.),
  fClExtraRoadZ(0.), 
  fExtraClErrYZ2(0), 
  fExtraClErrY2(0),
  fExtraClErrZ2(0),
  fPrimaryDCAZCut(-1),
  fPrimaryDCAYCut(-1),
  fDisableSecondaries(kFALSE),
  fCrossTalkSignalArray(0),
  fClPointersPool(0),
  fClPointersPoolPtr(0),
  fClPointersPoolSize(0),
  fSeedsPool(0),
  fHelixPool(0),
  fETPPool(0),
  fFreeSeedsID(500),
  fNFreeSeeds(0),
  fLastSeedID(-1),
  fAccountDistortions(0)
{
  //
  // default constructor
  //
  for (Int_t irow=0; irow<200; irow++){
    fXRow[irow]=0;
    fYMax[irow]=0;
    fPadLength[irow]=0;
  }
}
//_____________________________________________________________________



Int_t AliTPCtracker::UpdateTrack(AliTPCseed * track, Int_t accept){
  //
  //update track information using current cluster - track->fCurrentCluster


  AliTPCclusterMI* c =track->GetCurrentCluster();
  if (accept > 0) //sign not accepted clusters
    track->SetCurrentClusterIndex1(track->GetCurrentClusterIndex1() | 0x8000);  
  else // unsign accpeted clusters
    track->SetCurrentClusterIndex1(track->GetCurrentClusterIndex1() & 0xffff7fff);  
  UInt_t i = track->GetCurrentClusterIndex1();

  Int_t sec=(i&0xff000000)>>24; 
  Int_t row = (i&0x00ff0000)>>16; 
  if (sec>=fkParam->GetNInnerSector()) row += fkParam->GetNRowLow(); 
  track->SetRow(row);
  track->SetSector(sec);
  track->SetClusterIndex2(row, i);  
  //track->fFirstPoint = row;
  //if ( track->fLastPoint<row) track->fLastPoint =row;
  //  if (track->fRow<0 || track->fRow>=kMaxRow) {
  //  printf("problem\n");
  //}
  if (track->GetFirstPoint()>row) 
    track->SetFirstPoint(row);
  if (track->GetLastPoint()<row) 
    track->SetLastPoint(row);
  
  Double_t angle2 = track->GetSnp()*track->GetSnp();
  //
  //SET NEW Track Point
  //
  if (angle2<1) //PH sometimes angle2 is very big. To be investigated...
  {
    angle2 = TMath::Sqrt(angle2/(1-angle2));
    AliTPCTrackerPoints::Point   &point =*((AliTPCTrackerPoints::Point*)track->GetTrackPoint(row));
    //
    point.SetSigmaY(c->GetSigmaY2()/track->GetCurrentSigmaY2());
    point.SetSigmaZ(c->GetSigmaZ2()/track->GetCurrentSigmaZ2());
    float raty = track->GetErrorY2()>0 ? sqrt(track->GetErrorY2Syst())/track->GetErrorY2() : 0;
    float ratz = track->GetErrorZ2()>0 ? sqrt(track->GetErrorZ2Syst())/track->GetErrorZ2() : 0;
    point.SetErrYSys2TotSq(raty);
    point.SetErrZSys2TotSq(ratz);
    //    point.SetErrY(sqrt(track->GetErrorY2()));
    //    point.SetErrZ(sqrt(track->GetErrorZ2()));
    //
    point.SetX(track->GetX());
    point.SetY(track->GetY());
    point.SetZ(track->GetZ());
    point.SetAngleY(angle2);
    point.SetAngleZ(track->GetTgl());
    if (track->IsShared(row)){
      track->SetErrorY2(track->GetErrorY2()*4);
      track->SetErrorZ2(track->GetErrorZ2()*4);
    }
  }  

  if (accept>0) return 0;
  Double_t chi2 = track->GetPredictedChi2(track->GetCurrentCluster());
  //
//   track->SetErrorY2(track->GetErrorY2()*1.3);
//   track->SetErrorY2(track->GetErrorY2()+0.01);    
//   track->SetErrorZ2(track->GetErrorZ2()*1.3);   
//   track->SetErrorZ2(track->GetErrorZ2()+0.005);      
    //}
  if (track->GetNumberOfClusters()%20==0){
    //    if (track->fHelixIn){
    //  TClonesArray & larr = *(track->fHelixIn);    
    //  Int_t ihelix = larr.GetEntriesFast();
    //  new(larr[ihelix]) AliHelix(*track) ;    
    //}
  }
  if (AliTPCReconstructor::StreamLevel()&kStreamUpdateTrack) {
    Int_t event = (fEvent==NULL)? 0: fEvent->GetEventNumberInFile();
    AliExternalTrackParam param(*track);
    TTreeSRedirector &cstream = *fDebugStreamer;
    cstream<<"UpdateTrack"<<
      "cl.="<<c<<
      "event="<<event<<
      "track.="<<&param<<
      "\n";
  }
  track->SetNoCluster(0);
  return track->Update(c,chi2,i);
}



Int_t AliTPCtracker::AcceptCluster(AliTPCseed * seed, AliTPCclusterMI * cluster)
{
  //
  // decide according desired precision to accept given 
  // cluster for tracking
  Double_t  yt = seed->GetY(),zt = seed->GetZ();
  // RS: use propagation only if the seed in far from the cluster
  const double kTolerance = 10e-4; // assume track is at cluster X if X-distance below this
  if (TMath::Abs(seed->GetX()-cluster->GetX())>kTolerance) seed->GetProlongation(cluster->GetX(),yt,zt); 
  Double_t sy2=0;//ErrY2(seed,cluster);
  Double_t sz2=0;//ErrZ2(seed,cluster);
  ErrY2Z2(seed,cluster,sy2,sz2);
  
  Double_t sdistancey2 = sy2+seed->GetSigmaY2();
  Double_t sdistancez2 = sz2+seed->GetSigmaZ2();
  Double_t dy=seed->GetCurrentCluster()->GetY()-yt;
  Double_t dz=seed->GetCurrentCluster()->GetZ()-zt;
  Double_t rdistancey2 = dy*dy/sdistancey2;
  Double_t rdistancez2 = dz*dz/sdistancez2;
  
  Double_t rdistance2  = rdistancey2+rdistancez2;
  //Int_t  accept =0;
  
  if (AliTPCReconstructor::StreamLevel()>2 && ( (fIteration>0)|| (seed->GetNumberOfClusters()>20))) {
    //  if (AliTPCReconstructor::StreamLevel()>2 && seed->GetNumberOfClusters()>20) {
    Float_t rmsy2 = seed->GetCurrentSigmaY2();
    Float_t rmsz2 = seed->GetCurrentSigmaZ2();
    Float_t rmsy2p30 = seed->GetCMeanSigmaY2p30();
    Float_t rmsz2p30 = seed->GetCMeanSigmaZ2p30();
    Float_t rmsy2p30R  = seed->GetCMeanSigmaY2p30R();
    Float_t rmsz2p30R  = seed->GetCMeanSigmaZ2p30R();
    AliExternalTrackParam param(*seed);
    static TVectorD gcl(3),gtr(3);
    Float_t gclf[3];
    param.GetXYZ(gcl.GetMatrixArray());
    cluster->GetGlobalXYZ(gclf);
    gcl[0]=gclf[0];    gcl[1]=gclf[1];    gcl[2]=gclf[2];
    Int_t nclSeed=seed->GetNumberOfClusters();
    int seedType = seed->GetSeedType();
    if (AliTPCReconstructor::StreamLevel()&kStreamErrParam) { // flag:stream in debug mode cluster and track extrapolation at given row together with error nad shape estimate
      Int_t eventNr = fEvent->GetEventNumberInFile();
	
    (*fDebugStreamer)<<"ErrParam"<<
      "iter="<<fIteration<<
      "eventNr="<<eventNr<<
      "Cl.="<<cluster<<
      "nclSeed="<<nclSeed<<
      "seedType="<<seedType<<
      "T.="<<&param<<
      "dy="<<dy<<
      "dz="<<dz<<
      "yt="<<yt<<
      "zt="<<zt<<
      "gcl.="<<&gcl<<
      "gtr.="<<&gtr<<
      "erry2="<<sy2<<
      "errz2="<<sz2<<
      "rmsy2="<<rmsy2<<
      "rmsz2="<<rmsz2<<	
      "rmsy2p30="<<rmsy2p30<<
      "rmsz2p30="<<rmsz2p30<<	
      "rmsy2p30R="<<rmsy2p30R<<
      "rmsz2p30R="<<rmsz2p30R<<	
      // normalize distance - 
      "rdisty="<<rdistancey2<<
      "rdistz="<<rdistancez2<<
      "rdist="<<rdistance2<< //       
      "\n";
    }
  }
  //return 0;  // temporary
  if (rdistance2>32) return 3;
  
  
  if ((rdistancey2>9. || rdistancez2>9.) && cluster->GetType()==0)  
    return 2;  //suspisiouce - will be changed
  
  if ((rdistancey2>6.25 || rdistancez2>6.25) && cluster->GetType()>0)  
    // strict cut on overlaped cluster
    return  2;  //suspisiouce - will be changed
  
  if ( (rdistancey2>1. || rdistancez2>6.25 ) 
       && cluster->GetType()<0){
    seed->SetNFoundable(seed->GetNFoundable()-1);
    return 2;    
  }

  if (fUseHLTClusters == 3 || fUseHLTClusters == 4) {
    if (fIteration==2){
      if(!AliTPCReconstructor::GetRecoParam()->GetUseHLTOnePadCluster()) {
	if (TMath::Abs(cluster->GetSigmaY2()) < kAlmost0)
	  return 2;
      }
    }
  }

  return 0;
}





//_____________________________________________________________________________
AliTPCtracker::AliTPCtracker(const AliTPCParam *par): 
  AliTracker(), 
  fkNIS(par->GetNInnerSector()/2),
  fInnerSec(0),
  fkNOS(par->GetNOuterSector()/2),
  fOuterSec(0),
  fN(0),
  fSectors(0),
  fInput(0),
  fOutput(0),
  fSeedTree(0),
  fTreeDebug(0),
  fEvent(0),
  fEventHLT(0),
  fDebug(0),
  fNewIO(0),
  fNtracks(0),
  fSeeds(0),
  fIteration(0),
  fkParam(0), 
  fDebugStreamer(0),
  fUseHLTClusters(4),
  fClExtraRoadY(0.),
  fClExtraRoadZ(0.), 
  fExtraClErrYZ2(0), 
  fExtraClErrY2(0),
  fExtraClErrZ2(0),
  fPrimaryDCAZCut(-1),
  fPrimaryDCAYCut(-1),
  fDisableSecondaries(kFALSE),
  fCrossTalkSignalArray(0),
  fClPointersPool(0),
  fClPointersPoolPtr(0),
  fClPointersPoolSize(0),
  fSeedsPool(0),
  fHelixPool(0),
  fETPPool(0),
  fFreeSeedsID(500),
  fNFreeSeeds(0),
  fLastSeedID(-1),
  fAccountDistortions(0)
{
  //---------------------------------------------------------------------
  // The main TPC tracker constructor
  //---------------------------------------------------------------------
  fInnerSec=new AliTPCtrackerSector[fkNIS];         
  fOuterSec=new AliTPCtrackerSector[fkNOS];
 
  Int_t i;
  for (i=0; i<fkNIS; i++) fInnerSec[i].Setup(par,0);
  for (i=0; i<fkNOS; i++) fOuterSec[i].Setup(par,1);

  fkParam = par;  
  Int_t nrowlow = par->GetNRowLow();
  Int_t nrowup = par->GetNRowUp();

  
  for (i=0;i<nrowlow;i++){
    fXRow[i]     = par->GetPadRowRadiiLow(i);
    fPadLength[i]= par->GetPadPitchLength(0,i);
    fYMax[i]     = fXRow[i]*TMath::Tan(0.5*par->GetInnerAngle());
  }

  
  for (i=0;i<nrowup;i++){
    fXRow[i+nrowlow]      = par->GetPadRowRadiiUp(i);
    fPadLength[i+nrowlow] = par->GetPadPitchLength(60,i);
    fYMax[i+nrowlow]      = fXRow[i+nrowlow]*TMath::Tan(0.5*par->GetOuterAngle());
  }

  if (AliTPCReconstructor::StreamLevel()>0) {
    fDebugStreamer = new TTreeSRedirector("TPCdebug.root","recreate");
    AliTPCReconstructor::SetDebugStreamer(fDebugStreamer);
  }
  //
  fSeedsPool = new TClonesArray("AliTPCseed",1000);
  fClPointersPool = new AliTPCclusterMI*[kMaxFriendTracks*kMaxRow];
  memset(fClPointersPool,0,kMaxFriendTracks*kMaxRow*sizeof(AliTPCclusterMI*));
  fClPointersPoolPtr = fClPointersPool;
  fClPointersPoolSize = kMaxFriendTracks;
  //
  // crosstalk array and matrix initialization
  Int_t nROCs   = 72;
  Int_t nTimeBinsAll  = AliTPCcalibDB::Instance()->GetMaxTimeBinAllPads() ;
  Int_t nWireSegments = 11;
  fCrossTalkSignalArray = new TObjArray(nROCs*4);  //  
  fCrossTalkSignalArray->SetOwner(kTRUE);
  for (Int_t isector=0; isector<4*nROCs; isector++){
    TMatrixD * crossTalkSignal = new TMatrixD(nWireSegments,nTimeBinsAll);
    for (Int_t imatrix = 0; imatrix<11; imatrix++)
      for (Int_t jmatrix = 0; jmatrix<nTimeBinsAll; jmatrix++){
        (*crossTalkSignal)[imatrix][jmatrix]=0.;
      }
    fCrossTalkSignalArray->AddAt(crossTalkSignal,isector);
  }

}
//________________________________________________________________________
AliTPCtracker::AliTPCtracker(const AliTPCtracker &t):
  AliTracker(t),
  fkNIS(t.fkNIS),
  fInnerSec(0),
  fkNOS(t.fkNOS),
  fOuterSec(0),
  fN(0),
  fSectors(0),
  fInput(0),
  fOutput(0),
  fSeedTree(0),
  fTreeDebug(0),
  fEvent(0),
  fEventHLT(0),
  fDebug(0),
  fNewIO(kFALSE),
  fNtracks(0),
  fSeeds(0),
  fIteration(0),
  fkParam(0),
  fDebugStreamer(0),
  fUseHLTClusters(4),
  fClExtraRoadY(0.),
  fClExtraRoadZ(0.), 
  fExtraClErrYZ2(0), 
  fExtraClErrY2(0),
  fExtraClErrZ2(0),
  fPrimaryDCAZCut(-1),
  fPrimaryDCAYCut(-1),
  fDisableSecondaries(kFALSE),
  fCrossTalkSignalArray(0),
  fClPointersPool(0),
  fClPointersPoolPtr(0),
  fClPointersPoolSize(0),
  fSeedsPool(0),
  fHelixPool(0),
  fETPPool(0),
  fFreeSeedsID(500),
  fNFreeSeeds(0),
  fLastSeedID(-1),
  fAccountDistortions(0)
{
  //------------------------------------
  // dummy copy constructor
  //------------------------------------------------------------------
  fOutput=t.fOutput;
  for (Int_t irow=0; irow<200; irow++){
    fXRow[irow]=0;
    fYMax[irow]=0;
    fPadLength[irow]=0;
  }

}
AliTPCtracker & AliTPCtracker::operator=(const AliTPCtracker& /*r*/)
{
  //------------------------------
  // dummy 
  //--------------------------------------------------------------
  return *this;
}
//_____________________________________________________________________________
AliTPCtracker::~AliTPCtracker() {
  //------------------------------------------------------------------
  // TPC tracker destructor
  //------------------------------------------------------------------
  delete[] fInnerSec;
  delete[] fOuterSec;
  if (fSeeds) {
    fSeeds->Clear(); 
    delete fSeeds;
  }
  delete[] fClPointersPool;
  if (fCrossTalkSignalArray) delete fCrossTalkSignalArray;
  if (fDebugStreamer) delete fDebugStreamer;
  if (fSeedsPool) {
    for (int isd=fSeedsPool->GetEntriesFast();isd--;) {
      AliTPCseed* seed = (AliTPCseed*)fSeedsPool->At(isd);
      if (seed) {
	seed->SetClusterOwner(kFALSE);
	seed->SetClustersArrayTMP(0);
      }
    }
    delete fSeedsPool;
  }
  if (fHelixPool) fHelixPool->Delete();
  delete fHelixPool;
  if (fETPPool) fETPPool->Delete();
  delete fETPPool;
}


void AliTPCtracker::FillESD(const TObjArray* arr)
{
  //
  //
  //fill esds using updated tracks

  if (!fEvent) return;

  AliESDtrack iotrack;
  
    // write tracks to the event
    // store index of the track
    Int_t nseed=arr->GetEntriesFast();
    //FindKinks(arr,fEvent);
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
      if (!pt) continue; 
      pt->UpdatePoints();

      AddSystCovariance(pt); // correct covariance matrix for clusters syst. error

      AddCovariance(pt);
      if (AliTPCReconstructor::StreamLevel()&kStreamFillESD) {
	(*fDebugStreamer)<<"FillESD"<<  // flag: stream track information in FillESD function (after track Iteration 0)
	  "Tr0.="<<pt<<
	  "\n";       
      }
      //      pt->PropagateTo(fkParam->GetInnerRadiusLow());
      if (pt->GetKinkIndex(0)<=0){  //don't propagate daughter tracks 
	pt->PropagateTo(fkParam->GetInnerRadiusLow());
	pt->SetRow(0); // RS: memorise row
      }
      
      if (( pt->GetPoints()[2]- pt->GetPoints()[0])>5 && pt->GetPoints()[3]>0.8){
	iotrack.~AliESDtrack();
	new(&iotrack) AliESDtrack;
	iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);
        iotrack.SetTPCsignal(pt->GetdEdx(), pt->GetSDEDX(0), pt->GetNCDEDX(0)); 
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	//	iotrack.SetTPCpid(pt->fTPCr);
	//iotrack.SetTPCindex(i); 
	MakeESDBitmaps(pt, &iotrack);
	fEvent->AddTrack(&iotrack);
	continue;
      }
       
      if ( (pt->GetNumberOfClusters()>70)&& (Float_t(pt->GetNumberOfClusters())/Float_t(pt->GetNFoundable()))>0.55) {
	iotrack.~AliESDtrack();
	new(&iotrack) AliESDtrack;
	iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);
        iotrack.SetTPCsignal(pt->GetdEdx(), pt->GetSDEDX(0), pt->GetNCDEDX(0)); 
	iotrack.SetTPCPoints(pt->GetPoints());
	//iotrack.SetTPCindex(i);
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	MakeESDBitmaps(pt, &iotrack);
	//	iotrack.SetTPCpid(pt->fTPCr);
	fEvent->AddTrack(&iotrack);
	continue;
      } 
      //
      // short tracks  - maybe decays

      //RS Seed don't keep their cluster pointers, cache cluster usage stat. for fast evaluation
      // FillSeedClusterStatCache(pt); // RS use this slow method only if info on shared statistics is needed

      if ( (pt->GetNumberOfClusters()>30) && (Float_t(pt->GetNumberOfClusters())/Float_t(pt->GetNFoundable()))>0.70) {
	Int_t found,foundable; //,shared;
	//GetCachedSeedClusterStatistic(0,60,found, foundable,shared,kFALSE); // RS make sure FillSeedClusterStatCache is called
	//pt->GetClusterStatistic(0,60,found, foundable,shared,kFALSE); //RS: avoid this method: seed does not keep clusters
	pt->GetClusterStatistic(0,60,found, foundable);
	if ( (found>20) && (pt->GetNShared()/float(pt->GetNumberOfClusters())<0.2)){
	  iotrack.~AliESDtrack();
	  new(&iotrack) AliESDtrack;
	  iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);	
	  iotrack.SetTPCsignal(pt->GetdEdx(), pt->GetSDEDX(0), pt->GetNCDEDX(0)); 
	  //iotrack.SetTPCindex(i);
	  iotrack.SetTPCPoints(pt->GetPoints());
	  iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	  iotrack.SetV0Indexes(pt->GetV0Indexes());
	  MakeESDBitmaps(pt, &iotrack);
	  //iotrack.SetTPCpid(pt->fTPCr);
	  fEvent->AddTrack(&iotrack);
	  continue;
	}
      }       
      
      if ( (pt->GetNumberOfClusters()>20) && (Float_t(pt->GetNumberOfClusters())/Float_t(pt->GetNFoundable()))>0.8) {
	Int_t found,foundable;//,shared;
	//RS GetCachedSeedClusterStatistic(0,60,found, foundable,shared,kFALSE); // RS make sure FillSeedClusterStatCache is called
	//pt->GetClusterStatistic(0,60,found, foundable,shared,kFALSE); //RS: avoid this method: seed does not keep clusters
	pt->GetClusterStatistic(0,60,found,foundable);
	if (found<20) continue;
	if (pt->GetNShared()/float(pt->GetNumberOfClusters())>0.2) continue;
	//
	iotrack.~AliESDtrack();
	new(&iotrack) AliESDtrack;
	iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);	
        iotrack.SetTPCsignal(pt->GetdEdx(), pt->GetSDEDX(0), pt->GetNCDEDX(0)); 
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	MakeESDBitmaps(pt, &iotrack);
	//iotrack.SetTPCpid(pt->fTPCr);
	//iotrack.SetTPCindex(i);
	fEvent->AddTrack(&iotrack);
	continue;
      }   
      // short tracks  - secondaties
      //
      if ( (pt->GetNumberOfClusters()>30) ) {
	Int_t found,foundable;//,shared;
	//GetCachedSeedClusterStatistic(128,158,found, foundable,shared,kFALSE); // RS make sure FillSeedClusterStatCache is called
	//pt->GetClusterStatistic(128,158,found, foundable,shared,kFALSE);  //RS: avoid this method: seed does not keep clusters
	pt->GetClusterStatistic(128,158,found, foundable);
	if ( (found>20) && (pt->GetNShared()/float(pt->GetNumberOfClusters())<0.2) &&float(found)/float(foundable)>0.8){
	  iotrack.~AliESDtrack();
	  new(&iotrack) AliESDtrack;
	  iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);	
	  iotrack.SetTPCsignal(pt->GetdEdx(), pt->GetSDEDX(0), pt->GetNCDEDX(0)); 
	  iotrack.SetTPCPoints(pt->GetPoints());
	  iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	  iotrack.SetV0Indexes(pt->GetV0Indexes());
	  MakeESDBitmaps(pt, &iotrack);
	  //iotrack.SetTPCpid(pt->fTPCr);	
	  //iotrack.SetTPCindex(i);
	  fEvent->AddTrack(&iotrack);
	  continue;
	}
      }       
      
      if ( (pt->GetNumberOfClusters()>15)) {
	Int_t found,foundable,shared;
	//GetCachedSeedClusterStatistic(138,158,found, foundable,shared,kFALSE); // RS make sure FillSeedClusterStatCache is called
	//RS pt->GetClusterStatistic(138,158,found, foundable,shared,kFALSE);  //RS: avoid this method: seed does not keep clusters
	pt->GetClusterStatistic(138,158,found, foundable);
	if (found<15) continue;
	if (foundable<=0) continue;
	if (pt->GetNShared()/float(pt->GetNumberOfClusters())>0.2) continue;
	if (float(found)/float(foundable)<0.8) continue;
	//
	iotrack.~AliESDtrack();
	new(&iotrack) AliESDtrack;
	iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);	
        iotrack.SetTPCsignal(pt->GetdEdx(), pt->GetSDEDX(0), pt->GetNCDEDX(0)); 
	iotrack.SetTPCPoints(pt->GetPoints());
	iotrack.SetKinkIndexes(pt->GetKinkIndexes());
	iotrack.SetV0Indexes(pt->GetV0Indexes());
	MakeESDBitmaps(pt, &iotrack);
	//	iotrack.SetTPCpid(pt->fTPCr);
	//iotrack.SetTPCindex(i);
	fEvent->AddTrack(&iotrack);
	continue;
      }   
    }
    // >> account for suppressed tracks in the kink indices (RS)
    int nESDtracks = fEvent->GetNumberOfTracks();
    for (int it=nESDtracks;it--;) {
      AliESDtrack* esdTr = fEvent->GetTrack(it);
      if (!esdTr || !esdTr->GetKinkIndex(0)) continue;
      for (int ik=0;ik<3;ik++) {
	int knkId=0;
	if (!(knkId=esdTr->GetKinkIndex(ik))) break; // no more kinks for this track
	AliESDkink* kink = fEvent->GetKink(TMath::Abs(knkId)-1);
	if (!kink) {
	  AliError(Form("ESDTrack%d refers to non-existing kink %d",it,TMath::Abs(knkId)-1));
	  continue;
	}
	kink->SetIndex(it, knkId<0 ? 0:1); // update track index of the kink: mother at 0, daughter at 1
      }
    }

    // << account for suppressed tracks in the kink indices (RS)  
    AliInfo(Form("Number of filled ESDs-\t%d\n",fEvent->GetNumberOfTracks()));
  
}

//_____________________________________________________________________________________
void AliTPCtracker::ErrY2Z2(AliTPCseed* seed, const AliTPCclusterMI * cl, double &erry2, double &errz2)
{
  //
  //
  // Use calibrated cluster error from OCDB
  //
  const AliTPCRecoParam* rp = AliTPCReconstructor::GetRecoParam();
  AliTPCClusterParam * clparam = AliTPCcalibDB::Instance()->GetClusterParam();
  //
  Float_t z = TMath::Abs(fkParam->GetZLength(0)-TMath::Abs(seed->GetZ()));
  Int_t ctype = cl->GetType();  
  Int_t    type = (cl->GetRow()<63) ? 0: (cl->GetRow()>126) ? 1:2;
  double   snp2 = seed->GetSnp()*seed->GetSnp();
  if (snp2>1.-1e-6) snp2 = 1.-1e-6;
  double   tgp2 = snp2/(1.f-snp2);
  double   tgp  = TMath::Sqrt(tgp2);
  double   tgl2m = seed->GetTgl()*seed->GetTgl()*(1+tgp2); 
  Double_t tglm = TMath::Sqrt(tgl2m);
  erry2 = clparam->GetError0Par(0,type, z,tgp);
  errz2 = clparam->GetError0Par(1,type, z,tglm);
  if (ctype<0) { // edge cluster
    erry2+=0.5;  
    errz2+=0.5; 
  }
  erry2 *= erry2;
  errz2 *= errz2;
  // additional systematic error on the cluster
  double serry2=0,serrz2=0;
  AliTPCcalibDB::Instance()->GetTransform()->ErrY2Z2Syst(cl, tgp, seed->GetTgl(), serry2,serrz2);
  erry2 += serry2;
  errz2 += serrz2;
  seed->SetErrorY2(erry2);
  seed->SetErrorZ2(errz2);
  seed->SetErrorY2Syst(serry2);
  seed->SetErrorZ2Syst(serrz2);
 //
}



//_____________________________________________________________________________________
Double_t AliTPCtracker::ErrY2(AliTPCseed* seed, const AliTPCclusterMI * cl){
  //
  //
  // Use calibrated cluster error from OCDB
  //
  const AliTPCRecoParam* rp = AliTPCReconstructor::GetRecoParam();
  AliTPCClusterParam * clparam = AliTPCcalibDB::Instance()->GetClusterParam();
  //
  Float_t z = TMath::Abs(fkParam->GetZLength(0)-TMath::Abs(seed->GetZ()));
  Int_t ctype = cl->GetType();  
  Int_t    type = (cl->GetRow()<63) ? 0: (cl->GetRow()>126) ? 1:2;
  double   angle2 = seed->GetSnp()*seed->GetSnp();
  double   angle = TMath::Sqrt(TMath::Abs(angle2/(1.f-angle2)));
  Double_t erry2 = clparam->GetError0Par(0,type, z,angle);
  if (ctype<0) {
    erry2+=0.5;  // edge cluster
  }
  erry2 *= erry2;
  // additional systematic error on the cluster
  erry2 += AliTPCcalibDB::Instance()->GetTransform()->ErrY2Syst(cl, angle);
  seed->SetErrorY2(erry2);
  //
  return erry2;

//calculate look-up table at the beginning
//   static Bool_t  ginit = kFALSE;
//   static Float_t gnoise1,gnoise2,gnoise3;
//   static Float_t ggg1[10000];
//   static Float_t ggg2[10000];
//   static Float_t ggg3[10000];
//   static Float_t glandau1[10000];
//   static Float_t glandau2[10000];
//   static Float_t glandau3[10000];
//   //
//   static Float_t gcor01[500];
//   static Float_t gcor02[500];
//   static Float_t gcorp[500];
//   //

//   //
//   if (ginit==kFALSE){
//     for (Int_t i=1;i<500;i++){
//       Float_t rsigma = float(i)/100.;
//       gcor02[i] = TMath::Max(0.78 +TMath::Exp(7.4*(rsigma-1.2)),0.6);
//       gcor01[i] = TMath::Max(0.72 +TMath::Exp(3.36*(rsigma-1.2)),0.6);
//       gcorp[i]  = TMath::Max(TMath::Power((rsigma+0.5),1.5),1.2);
//     }

//     //
//     for (Int_t i=3;i<10000;i++){
//       //
//       //
//       // inner sector
//       Float_t amp = float(i);
//       Float_t padlength =0.75;
//       gnoise1 = 0.0004/padlength;
//       Float_t nel     = 0.268*amp;
//       Float_t nprim   = 0.155*amp;
//       ggg1[i]          = fkParam->GetDiffT()*fkParam->GetDiffT()*(2+0.001*nel/(padlength*padlength))/nel;
//       glandau1[i]      = (2.+0.12*nprim)*0.5* (2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
//       if (glandau1[i]>1) glandau1[i]=1;
//       glandau1[i]*=padlength*padlength/12.;      
//       //
//       // outer short
//       padlength =1.;
//       gnoise2   = 0.0004/padlength;
//       nel       = 0.3*amp;
//       nprim     = 0.133*amp;
//       ggg2[i]      = fkParam->GetDiffT()*fkParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
//       glandau2[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
//       if (glandau2[i]>1) glandau2[i]=1;
//       glandau2[i]*=padlength*padlength/12.;
//       //
//       //
//       // outer long
//       padlength =1.5;
//       gnoise3   = 0.0004/padlength;
//       nel       = 0.3*amp;
//       nprim     = 0.133*amp;
//       ggg3[i]      = fkParam->GetDiffT()*fkParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
//       glandau3[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
//       if (glandau3[i]>1) glandau3[i]=1;
//       glandau3[i]*=padlength*padlength/12.;
//       //
//     }
//     ginit = kTRUE;
//   }
//   //
//   //
//   //
//   Int_t amp = int(TMath::Abs(cl->GetQ()));  
//   if (amp>9999) {
//     seed->SetErrorY2(1.);
//     return 1.;
//   }
//   Float_t snoise2;
//   Float_t z = TMath::Abs(fkParam->GetZLength(0)-TMath::Abs(seed->GetZ()));
//   Int_t ctype = cl->GetType();  
//   Float_t padlength= GetPadPitchLength(seed->GetRow());
//   Double_t angle2 = seed->GetSnp()*seed->GetSnp();
//   angle2 = angle2/(1-angle2); 
//   //
//   //cluster "quality"
//   Int_t rsigmay = int(100.*cl->GetSigmaY2()/(seed->GetCurrentSigmaY2()));
//   Float_t res;
//   //
//   if (fSectors==fInnerSec){
//     snoise2 = gnoise1;
//     res     = ggg1[amp]*z+glandau1[amp]*angle2;     
//     if (ctype==0) res *= gcor01[rsigmay];
//     if ((ctype>0)){
//       res+=0.002;
//       res*= gcorp[rsigmay];
//     }
//   }
//   else {
//     if (padlength<1.1){
//       snoise2 = gnoise2;
//       res     = ggg2[amp]*z+glandau2[amp]*angle2; 
//       if (ctype==0) res *= gcor02[rsigmay];      
//       if ((ctype>0)){
// 	res+=0.002;
// 	res*= gcorp[rsigmay];
//       }
//     }
//     else{
//       snoise2 = gnoise3;      
//       res     = ggg3[amp]*z+glandau3[amp]*angle2; 
//       if (ctype==0) res *= gcor02[rsigmay];
//       if ((ctype>0)){
// 	res+=0.002;
// 	res*= gcorp[rsigmay];
//       }
//     }
//   }  

//   if (ctype<0){
//     res+=0.005;
//     res*=2.4;  // overestimate error 2 times
//   }
//   res+= snoise2;
 
//   if (res<2*snoise2)
//     res = 2*snoise2;
  
//   seed->SetErrorY2(res);
//   return res;


}


//__________________________________________________________________________________
Double_t AliTPCtracker::ErrZ2(AliTPCseed* seed, const AliTPCclusterMI * cl){
  //
  //
  // Use calibrated cluster error from OCDB
  //
  const AliTPCRecoParam* rp = AliTPCReconstructor::GetRecoParam();
  AliTPCClusterParam * clparam = AliTPCcalibDB::Instance()->GetClusterParam();
  //
  Float_t z = TMath::Abs(fkParam->GetZLength(0)-TMath::Abs(seed->GetZ()));
  Int_t ctype = cl->GetType();  
  Int_t    type = (cl->GetRow()<63) ? 0: (cl->GetRow()>126) ? 1:2;
  //
  double angle2 = seed->GetSnp()*seed->GetSnp();
  angle2 = seed->GetTgl()*seed->GetTgl()*(1+angle2/(1-angle2)); 
  Double_t angle = TMath::Sqrt(TMath::Abs(angle2));
  Double_t errz2 = clparam->GetError0Par(1,type, z,angle);
  if (ctype<0) {
    errz2+=0.5;  // edge cluster
  }
  errz2*=errz2;
  // additional systematic error on the cluster
  errz2 += AliTPCcalibDB::Instance()->GetTransform()->ErrZ2Syst(cl, seed->GetTgl());
  seed->SetErrorZ2(errz2);
  //
  return errz2;



//   //seed->SetErrorY2(0.1);
//   //return 0.1;
//   //calculate look-up table at the beginning
//   static Bool_t  ginit = kFALSE;
//   static Float_t gnoise1,gnoise2,gnoise3;
//   static Float_t ggg1[10000];
//   static Float_t ggg2[10000];
//   static Float_t ggg3[10000];
//   static Float_t glandau1[10000];
//   static Float_t glandau2[10000];
//   static Float_t glandau3[10000];
//   //
//   static Float_t gcor01[1000];
//   static Float_t gcor02[1000];
//   static Float_t gcorp[1000];
//   //

//   //
//   if (ginit==kFALSE){
//     for (Int_t i=1;i<1000;i++){
//       Float_t rsigma = float(i)/100.;
//       gcor02[i] = TMath::Max(0.81 +TMath::Exp(6.8*(rsigma-1.2)),0.6);
//       gcor01[i] = TMath::Max(0.72 +TMath::Exp(2.04*(rsigma-1.2)),0.6);
//       gcorp[i]  = TMath::Max(TMath::Power((rsigma+0.5),1.5),1.2);
//     }

//     //
//     for (Int_t i=3;i<10000;i++){
//       //
//       //
//       // inner sector
//       Float_t amp = float(i);
//       Float_t padlength =0.75;
//       gnoise1 = 0.0004/padlength;
//       Float_t nel     = 0.268*amp;
//       Float_t nprim   = 0.155*amp;
//       ggg1[i]          = fkParam->GetDiffT()*fkParam->GetDiffT()*(2+0.001*nel/(padlength*padlength))/nel;
//       glandau1[i]      = (2.+0.12*nprim)*0.5* (2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
//       if (glandau1[i]>1) glandau1[i]=1;
//       glandau1[i]*=padlength*padlength/12.;      
//       //
//       // outer short
//       padlength =1.;
//       gnoise2   = 0.0004/padlength;
//       nel       = 0.3*amp;
//       nprim     = 0.133*amp;
//       ggg2[i]      = fkParam->GetDiffT()*fkParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
//       glandau2[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
//       if (glandau2[i]>1) glandau2[i]=1;
//       glandau2[i]*=padlength*padlength/12.;
//       //
//       //
//       // outer long
//       padlength =1.5;
//       gnoise3   = 0.0004/padlength;
//       nel       = 0.3*amp;
//       nprim     = 0.133*amp;
//       ggg3[i]      = fkParam->GetDiffT()*fkParam->GetDiffT()*(2+0.0008*nel/(padlength*padlength))/nel;
//       glandau3[i]  = (2.+0.12*nprim)*0.5*(2.+nprim*nprim*0.001/(padlength*padlength))/nprim;
//       if (glandau3[i]>1) glandau3[i]=1;
//       glandau3[i]*=padlength*padlength/12.;
//       //
//     }
//     ginit = kTRUE;
//   }
//   //
//   //
//   //
//   Int_t amp = int(TMath::Abs(cl->GetQ()));  
//   if (amp>9999) {
//     seed->SetErrorY2(1.);
//     return 1.;
//   }
//   Float_t snoise2;
//   Float_t z = TMath::Abs(fkParam->GetZLength(0)-TMath::Abs(seed->GetZ()));
//   Int_t ctype = cl->GetType();  
//   Float_t padlength= GetPadPitchLength(seed->GetRow());
//   //
//   Double_t angle2 = seed->GetSnp()*seed->GetSnp();
//   //  if (angle2<0.6) angle2 = 0.6;
//   angle2 = seed->GetTgl()*seed->GetTgl()*(1+angle2/(1-angle2)); 
//   //
//   //cluster "quality"
//   Int_t rsigmaz = int(100.*cl->GetSigmaZ2()/(seed->GetCurrentSigmaZ2()));
//   Float_t res;
//   //
//   if (fSectors==fInnerSec){
//     snoise2 = gnoise1;
//     res     = ggg1[amp]*z+glandau1[amp]*angle2;     
//     if (ctype==0) res *= gcor01[rsigmaz];
//     if ((ctype>0)){
//       res+=0.002;
//       res*= gcorp[rsigmaz];
//     }
//   }
//   else {
//     if (padlength<1.1){
//       snoise2 = gnoise2;
//       res     = ggg2[amp]*z+glandau2[amp]*angle2; 
//       if (ctype==0) res *= gcor02[rsigmaz];      
//       if ((ctype>0)){
// 	res+=0.002;
// 	res*= gcorp[rsigmaz];
//       }
//     }
//     else{
//       snoise2 = gnoise3;      
//       res     = ggg3[amp]*z+glandau3[amp]*angle2; 
//       if (ctype==0) res *= gcor02[rsigmaz];
//       if ((ctype>0)){
// 	res+=0.002;
// 	res*= gcorp[rsigmaz];
//       }
//     }
//   }  

//   if (ctype<0){
//     res+=0.002;
//     res*=1.3;
//   }
//   if ((ctype<0) &&amp<70){
//     res+=0.002;
//     res*=1.3;  
//   }
//   res += snoise2;
//   if (res<2*snoise2)
//      res = 2*snoise2;
//   if (res>3) res =3;
//   seed->SetErrorZ2(res);
//   return res;
}





void AliTPCtracker::RotateToLocal(AliTPCseed *seed)
{
  //rotate to track "local coordinata
  Float_t x = seed->GetX();
  Float_t y = seed->GetY();
  Float_t ymax = x*TMath::Tan(0.5*fSectors->GetAlpha());
  
  if (y > ymax) {
    seed->SetRelativeSector((seed->GetRelativeSector()+1) % fN);
    if (!seed->Rotate(fSectors->GetAlpha())) 
      return;
  } else if (y <-ymax) {
    seed->SetRelativeSector((seed->GetRelativeSector()-1+fN) % fN);
    if (!seed->Rotate(-fSectors->GetAlpha())) 
      return;
  }   

}



//_____________________________________________________________________________
Double_t AliTPCtracker::F1old(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) const
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  if ( xr*xr+yr*yr<=0.00000000000001) return 100;
  return -xr*yr/sqrt(xr*xr+yr*yr); 
}



//_____________________________________________________________________________
Double_t AliTPCtracker::F1(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) const
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u;
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  return c2;
}


Double_t AliTPCtracker::F2(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) const 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //  
  Double_t det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10) {
    return 100;
  }
  //
  Double_t u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  Double_t x0 = x3*0.5-y3*u; 
  Double_t y0 = y3*0.5+x3*u;
  Double_t c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  if (det<0) c2*=-1;
  x0+=x1;
  x0*=c2;  
  return x0;
}



//_____________________________________________________________________________
Double_t AliTPCtracker::F2old(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) const
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
Double_t AliTPCtracker::F3(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) const
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}


Double_t AliTPCtracker::F3n(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2, Double_t c) const
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------

  //  Double_t angle1;
  
  //angle1    =  (z1-z2)*c/(TMath::ASin(c*x1-ni)-TMath::ASin(c*x2-ni));
  //
  Double_t d  =  TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  if (TMath::Abs(d*c*0.5)>1) return 0;
  Double_t   angle2    = asinf(d*c*0.5);

  angle2  = (z1-z2)*c/(angle2*2.);
  return angle2;
}

Bool_t   AliTPCtracker::GetProlongation(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z) const
{//-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=x2.
  //-----------------------------------------------------------------
  
  Double_t dx=x2-x1;
  Double_t c1=x[4]*x1 - x[2];
  if (TMath::Abs(c1) >= 0.999) {
    return kFALSE;
  }
  Double_t c2=x[4]*x2 - x[2];
  if (TMath::Abs(c2) >= 0.999) {
    return kFALSE;
  }
  Double_t r1=TMath::Sqrt((1.-c1)*(1.+c1)),r2=TMath::Sqrt((1.-c2)*(1.+c2));  
  y = x[0];
  z = x[1];
  
  Double_t dy = dx*(c1+c2)/(r1+r2);
  //
  Double_t delta = x[4]*dx*(c1+c2)/(c1*r2 + c2*r1);
  Double_t dz = x[3]*asinf(delta)/x[4];
  
  y+=dy;
  z+=dz;
  
  return kTRUE;  
}

Bool_t   AliTPCtracker::GetProlongationLine(Double_t x1, Double_t x2, Double_t x[5], Double_t &y, Double_t &z) const
{//-----------------------------------------------------------------
  // This function find straight line prolongation of a track to a reference plane x=x2.
  //-----------------------------------------------------------------
  
  if (TMath::Abs(x[2]) >= 0.999) return kFALSE;
  Double_t c1=- x[2], dx2r = (x2-x1)/TMath::Sqrt((1.-c1)*(1.+c1));
  y = x[0] + dx2r*c1;
  z = x[1] + dx2r*x[3];
  return kTRUE;  
}


Int_t  AliTPCtracker::LoadClusters (TTree *const tree)
{
  // load clusters
  //
  fInput = tree;
  return LoadClusters();
}


Int_t  AliTPCtracker::LoadClusters(const TObjArray *arr)
{
  //
  // load clusters to the memory
  AliTPCClustersRow *clrow = 0; //RS: why this new? new AliTPCClustersRow("AliTPCclusterMI");
  Int_t lower   = arr->LowerBound();
  Int_t entries = arr->GetEntriesFast();

  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform() ;
  transform->SetCurrentRecoParam((AliTPCRecoParam*)AliTPCReconstructor::GetRecoParam());
  transform->SetCurrentTimeStamp( GetTimeStamp());
  transform->SetCurrentRun( GetRunNumber() );

  AliWarning("Sector Change ins not checked in LoadClusters(const TObjArray *arr)");

  for (Int_t i=lower; i<entries; i++) {
    clrow = (AliTPCClustersRow*) arr->At(i);
    if(!clrow) continue;
    TClonesArray* arr = clrow->GetArray();
    if(!arr) continue;
    int ncl = arr->GetEntriesFast();
    if (ncl<1) continue;
    //  
    Int_t sec,row;
    fkParam->AdjustSectorRow(clrow->GetID(),sec,row);
    
    for (Int_t icl=ncl; icl--;) {
      Transform((AliTPCclusterMI*)(arr->At(icl)));
    }
    //
    // RS: Check for possible sector change due to the distortions: TODO
    //
    AliTPCtrackerRow * tpcrow=0;
    Int_t left=0;
    if (sec<fkNIS*2){
      tpcrow = &(fInnerSec[sec%fkNIS][row]);    
      left = sec/fkNIS;
    }
    else{
      tpcrow = &(fOuterSec[(sec-fkNIS*2)%fkNOS][row]);
      left = (sec-fkNIS*2)/fkNOS;
    }
    if (left ==0){
      tpcrow->SetN1(ncl);
      for (Int_t j=0;j<ncl;++j) 
	tpcrow->SetCluster1(j, *(AliTPCclusterMI*)(arr->At(j)));
    }
    if (left ==1){
      tpcrow->SetN2(ncl);
      for (Int_t j=0;j<ncl;++j) 
	tpcrow->SetCluster2(j, *(AliTPCclusterMI*)(arr->At(j)));
    }
    clrow->GetArray()->Clear(); // RS AliTPCclusterMI does not allocate memory
  }
  //
  //  delete clrow;
  LoadOuterSectors();
  LoadInnerSectors();
  return 0;
}

Int_t  AliTPCtracker::LoadClusters(const TClonesArray *arr)
{
  //
  // load clusters to the memory from one 
  // TClonesArray
  //
  // RS: Check for possible sector change due to the distortions: TODO
  AliWarning("Sector Change ins not checked in LoadClusters(const TClonesArray *arr)");
  //

  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform() ;
  transform->SetCurrentRecoParam((AliTPCRecoParam*)AliTPCReconstructor::GetRecoParam());
  transform->SetCurrentTimeStamp( GetTimeStamp());
  transform->SetCurrentRun( GetRunNumber() );
  //
  AliTPCclusterMI *clust=0;
  Int_t count[72][96] = { {0} , {0} }; 

  // loop over clusters
  for (Int_t icl=0; icl<arr->GetEntriesFast(); icl++) {
    clust = (AliTPCclusterMI*)arr->At(icl);
    if(!clust) continue;
    //printf("cluster: det %d, row %d \n", clust->GetDetector(),clust->GetRow());

    // transform clusters
    Transform(clust);

    // count clusters per pad row
    count[clust->GetDetector()][clust->GetRow()]++;
  }

  // insert clusters to sectors
  for (Int_t icl=0; icl<arr->GetEntriesFast(); icl++) {
    clust = (AliTPCclusterMI*)arr->At(icl);
    if(!clust) continue;

    Int_t sec = clust->GetDetector();
    Int_t row = clust->GetRow();

    // filter overlapping pad rows needed by HLT
    if(sec<fkNIS*2) { //IROCs
     if(row == 30) continue;
    }
    else { // OROCs
      if(row == 27 || row == 76) continue;
    }

    //    Int_t left=0;
    if (sec<fkNIS*2){
      //      left = sec/fkNIS;
      fInnerSec[sec%fkNIS].InsertCluster(clust, count[sec][row], fkParam);    
    }
    else{
      //      left = (sec-fkNIS*2)/fkNOS;
      fOuterSec[(sec-fkNIS*2)%fkNOS].InsertCluster(clust, count[sec][row], fkParam);
    }
  }

  // Load functions must be called behind LoadCluster(TClonesArray*)
  // needed by HLT
  //LoadOuterSectors();
  //LoadInnerSectors();

  return 0;
}

Int_t  AliTPCtracker::LoadClusters()
{
 //
  // load clusters to the memory
  static AliTPCClustersRow *clrow= new AliTPCClustersRow("AliTPCclusterMI");
  //

  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform() ;
  transform->SetCurrentRecoParam((AliTPCRecoParam*)AliTPCReconstructor::GetRecoParam());
  transform->SetCurrentTimeStamp( GetTimeStamp());
  transform->SetCurrentRun( GetRunNumber() );

  //  TTree * tree = fClustersArray.GetTree();
  AliInfo("LoadClusters()\n");

  fNClusters = 0;

  TTree * tree = fInput;
  TBranch * br = tree->GetBranch("Segment");
  br->SetAddress(&clrow);

  // Conversion of pad, row coordinates in local tracking coords.
  // Could be skipped here; is already done in clusterfinder

  double cutZ2X = AliTPCReconstructor::GetPrimaryZ2XCut();
  double cutZOutSector = AliTPCReconstructor::GetZOutSectorCut();
  if (cutZOutSector>0 && AliTPCReconstructor::GetExtendedRoads()) 
    cutZOutSector += AliTPCReconstructor::GetExtendedRoads()[1];
  //
  if (cutZ2X>0 || cutZOutSector>0) {
    AliInfoF("Cut on cluster |Z/X| : %s, on cluster Z on wrong CE side: %s",
	     cutZ2X>0        ? Form("%.3f",cutZ2X) : "N/A",
	     cutZOutSector>0 ? Form("%.3f",cutZOutSector) : "N/A");
  }
  Int_t j=Int_t(tree->GetEntries());
  for (Int_t i=0; i<j; i++) {
    br->GetEntry(i);
    //
    TClonesArray* clArr = clrow->GetArray();
    int nClus = clArr->GetEntriesFast();
    Int_t sec,row;
    fkParam->AdjustSectorRow(clrow->GetID(),sec,row);
    for (Int_t icl=nClus; icl--;) Transform((AliTPCclusterMI*)(clArr->At(icl)));
    //
    AliTPCtrackerRow * tpcrow=0;
    Int_t left=0;
    if (sec<fkNIS*2){
      tpcrow = &(fInnerSec[sec%fkNIS][row]);    
      left = sec/fkNIS;
    }
    else{
      tpcrow = &(fOuterSec[(sec-fkNIS*2)%fkNOS][row]);
      left = (sec-fkNIS*2)/fkNOS;
    }
    int nClusAdd = 0;
    if (left ==0){  // A side     
      for (int k=0;k<nClus;k++) {
	const AliTPCclusterMI& cl = *((AliTPCclusterMI*)clArr->At(k));
	if (cutZOutSector>0 && cl.GetZ()<-cutZOutSector) continue;
	if (cutZ2X>0 && cl.GetZ()/cl.GetX() > cutZ2X) continue;
	tpcrow->SetCluster1(nClusAdd++, cl);
      }
      tpcrow->SetN1(nClusAdd);
    }
    if (left ==1){ // C side
      for (int k=0;k<nClus;k++) {
	const AliTPCclusterMI& cl = *((AliTPCclusterMI*)clArr->At(k));
	if (cutZOutSector>0 && cl.GetZ()>cutZOutSector) continue;
	if (cutZ2X>0 && cl.GetZ()/cl.GetX() < -cutZ2X) continue;
	tpcrow->SetCluster2(nClusAdd++, cl);
      }
      tpcrow->SetN2(nClusAdd);
    }
    fNClusters += nClusAdd;
  }
  //
  //clrow->Clear("C");
  clrow->Clear(); // RS AliTPCclusterMI does not allocate memory
  LoadOuterSectors();
  LoadInnerSectors();

  cout << " =================================================================================================================================== " << endl;
  cout << " AliTPCReconstructor::GetRecoParam()->GetUseIonTailCorrection() =  " << AliTPCReconstructor::GetRecoParam()->GetUseIonTailCorrection() << endl;
  cout << " AliTPCReconstructor::GetRecoParam()->GetCrosstalkCorrection()  =  " << AliTPCReconstructor::GetRecoParam()->GetCrosstalkCorrection()  << endl;
  cout << " =================================================================================================================================== " << endl;

  if (AliTPCReconstructor::GetRecoParam()->GetUseIonTailCorrection()) ApplyTailCancellation();
  if (AliTPCReconstructor::GetRecoParam()->GetCrosstalkCorrection()!=0.) CalculateXtalkCorrection();
  if (AliTPCReconstructor::GetRecoParam()->GetCrosstalkCorrection()!=0.) ApplyXtalkCorrection();
  //if (AliTPCReconstructor::GetRecoParam()->GetUseOulierClusterFilter()) FilterOutlierClusters();  

  /*
  static int maxClus[18][2][kMaxRow]={0};
  int maxAcc=0,nclEv=0, capacity=0;
  for (int isec=0;isec<18;isec++) {
    for (int irow=0;irow<kMaxRow;irow++) {
      AliTPCtrackerRow * tpcrow = irow>62 ?  &(fOuterSec[isec][irow-63]) : &(fInnerSec[isec][irow]);
      maxClus[isec][0][irow] = TMath::Max(maxClus[isec][0][irow], tpcrow->GetN1());
      maxClus[isec][1][irow] = TMath::Max(maxClus[isec][1][irow], tpcrow->GetN2());
      maxAcc += maxClus[isec][0][irow]+maxClus[isec][1][irow];
      nclEv += tpcrow->GetN();
      capacity += tpcrow->GetClusters1()->Capacity();
      capacity += tpcrow->GetClusters2()->Capacity();
    }
  }
  printf("RS:AccumulatedSpace: %d for %d | pointers: %d\n",maxAcc,nclEv,capacity);
  */
  return 0;
}

void  AliTPCtracker::CalculateXtalkCorrection(){
  //
  // Calculate crosstalk estimate
  //
  TStopwatch sw;
  sw.Start();
  const Int_t nROCs   = 72;
  const Int_t   nIterations=3;  // 
  // 0.) reset crosstalk matrix 
  //
  for (Int_t isector=0; isector<nROCs*4; isector++){  //set all ellemts of crosstalk matrix to 0 
    TMatrixD * crossTalkMatrix = (TMatrixD*)fCrossTalkSignalArray->At(isector);
    if (crossTalkMatrix)(*crossTalkMatrix)*=0;
  }
  
  //
  // 1.) Filling part -- loop over clusters
  //
  Double_t missingChargeFactor= AliTPCReconstructor::GetRecoParam()->GetCrosstalkCorrectionMissingCharge();
  for (Int_t iter=0; iter<nIterations;iter++){
    for (Int_t isector=0; isector<36; isector++){      // loop over sectors
      for (Int_t iside=0; iside<2; iside++){           // loop over sides A/C
	AliTPCtrackerSector &sector= (isector<18)?fInnerSec[isector%18]:fOuterSec[isector%18];
	Int_t nrows = sector.GetNRows();
	Int_t sec=0;
	if (isector<18) sec=isector+18*iside;
	if (isector>=18) sec=18+isector+18*iside;
	for (Int_t row = 0;row<nrows;row++){           // loop over rows       
	  //
	  //
	  Int_t wireSegmentID     = fkParam->GetWireSegment(sec,row);
	  Float_t nPadsPerSegment = (Float_t)(fkParam->GetNPadsPerSegment(wireSegmentID));
	  TMatrixD &crossTalkSignal =  *((TMatrixD*)fCrossTalkSignalArray->At(sec));
	  TMatrixD &crossTalkSignalCache =  *((TMatrixD*)fCrossTalkSignalArray->At(sec+nROCs*2));   // this is the cache value of the crosstalk from previous iteration
	  TMatrixD &crossTalkSignalBelow =  *((TMatrixD*)fCrossTalkSignalArray->At(sec+nROCs));
	  Int_t nCols=crossTalkSignal.GetNcols();
	  //
	  AliTPCtrackerRow&  tpcrow = sector[row];       
	  Int_t ncl = tpcrow.GetN1();                  // number of clusters in the row
	  if (iside>0) ncl=tpcrow.GetN2();
	  for (Int_t i=0;i<ncl;i++) {  // loop over clusters
	    AliTPCclusterMI *clXtalk= (iside>0)?(tpcrow.GetCluster2(i)):(tpcrow.GetCluster1(i));
	    
	    Int_t timeBinXtalk = clXtalk->GetTimeBin();      
	    Double_t rmsPadMin2=0.5*0.5+(fkParam->GetDiffT()*fkParam->GetDiffT())*(TMath::Abs((clXtalk->GetZ()-fkParam->GetZLength())))/(fkParam->GetPadPitchWidth(sec)*fkParam->GetPadPitchWidth(sec)); // minimal PRF width - 0.5 is the PRF in the pad units - we should et it from fkparam getters 
	    Double_t rmsTimeMin2=1+(fkParam->GetDiffL()*fkParam->GetDiffL())*(TMath::Abs((clXtalk->GetZ()-fkParam->GetZLength())))/(fkParam->GetZWidth()*fkParam->GetZWidth()); // minimal PRF width - 1 is the TRF in the time bin units - we should et it from fkParam getters
	    Double_t rmsTime2   = clXtalk->GetSigmaZ2()/(fkParam->GetZWidth()*fkParam->GetZWidth()); 
	    Double_t rmsPad2    = clXtalk->GetSigmaY2()/(fkParam->GetPadPitchWidth(sec)*fkParam->GetPadPitchWidth(sec)); 
	    if (rmsPadMin2>rmsPad2){
	      rmsPad2=rmsPadMin2;
	    }
	    if (rmsTimeMin2>rmsTime2){
	      rmsTime2=rmsTimeMin2;
	    }
	    
	    Double_t norm= 2.*TMath::Exp(-1.0/(2.*rmsTime2))+2.*TMath::Exp(-4.0/(2.*rmsTime2))+1.;
	    Double_t qTotXtalk = 0.;   
	    Double_t qTotXtalkMissing = 0.;   
	    for (Int_t itb=timeBinXtalk-2, idelta=-2; itb<=timeBinXtalk+2; itb++,idelta++) {        
	      if (itb<0 || itb>=nCols) continue;
	      Double_t missingCharge=0;
	      Double_t trf= TMath::Exp(-idelta*idelta/(2.*rmsTime2));
	      if (missingChargeFactor>0) {
		for (Int_t dpad=-2; dpad<=2; dpad++){
		  Double_t qPad =   clXtalk->GetMax()*TMath::Exp(-dpad*dpad/(2.*rmsPad2))*trf;
		  if (TMath::Nint(qPad-crossTalkSignalCache[wireSegmentID][itb])<=fkParam->GetZeroSup()){
		    missingCharge+=qPad+crossTalkSignalCache[wireSegmentID][itb];
		  }else{
		    missingCharge+=crossTalkSignalCache[wireSegmentID][itb];
		  }
		}
	      }
	      qTotXtalk = clXtalk->GetQ()*trf/norm+missingCharge*missingChargeFactor;
	      qTotXtalkMissing = missingCharge;
	      crossTalkSignal[wireSegmentID][itb]+= qTotXtalk/nPadsPerSegment; 
	      crossTalkSignalBelow[wireSegmentID][itb]+= qTotXtalkMissing/nPadsPerSegment; 
	    } // end of time bin loop
	  } // end of cluster loop	
	} // end of rows loop
      }  // end of side loop
    }    // end of sector loop
    //
    // copy crosstalk matrix to cached used for next itteration
    //
    //
    // 2.) dump the crosstalk matrices to tree for further investigation
    //     a.) to estimate fluctuation of pedestal in indiviula wire segments
    //     b.) to check correlation between regions
    //     c.) to check relative conribution of signal below threshold to crosstalk
    
    if (AliTPCReconstructor::StreamLevel()&kStreamCrosstalkMatrix) {
      for (Int_t isector=0; isector<nROCs; isector++){  //set all ellemts of crosstalk matrix to 0
	TMatrixD * crossTalkMatrix = (TMatrixD*)fCrossTalkSignalArray->At(isector);
	TMatrixD * crossTalkMatrixBelow = (TMatrixD*)fCrossTalkSignalArray->At(isector+nROCs);
	TMatrixD * crossTalkMatrixCache = (TMatrixD*)fCrossTalkSignalArray->At(isector+nROCs*2);
	TVectorD vecAll(crossTalkMatrix->GetNrows());
	TVectorD vecBelow(crossTalkMatrix->GetNrows());
	TVectorD vecCache(crossTalkMatrixCache->GetNrows());
	//
	for (Int_t itime=0; itime<crossTalkMatrix->GetNcols(); itime++){
	  for (Int_t iwire=0; iwire<crossTalkMatrix->GetNrows(); iwire++){
	    vecAll[iwire]=(*crossTalkMatrix)(iwire,itime);
	    vecBelow[iwire]=(*crossTalkMatrixBelow)(iwire,itime);
	    vecCache[iwire]=(*crossTalkMatrixCache)(iwire,itime);
	  }
	  (*fDebugStreamer)<<"crosstalkMatrix"<<
	    "iter="<<iter<<                      //iteration
	    "sector="<<isector<<                 // sector
	    "itime="<<itime<<                    // time bin index
	    "vecAll.="<<&vecAll<<                // crosstalk charge + charge below threshold
	    "vecCache.="<<&vecCache<<                // crosstalk charge + charge below threshold	  
	    "vecBelow.="<<&vecBelow<<            // crosstalk contribution from signal below threshold
	    "\n";
	}
      }
    }
    if (iter<nIterations-1) for (Int_t isector=0; isector<nROCs*2; isector++){  //set all ellemts of crosstalk matrix to 0 
      TMatrixD * crossTalkMatrix = (TMatrixD*)fCrossTalkSignalArray->At(isector);
      TMatrixD * crossTalkMatrixCache = (TMatrixD*)fCrossTalkSignalArray->At(isector+nROCs*2);
      if (crossTalkMatrix){
	(*crossTalkMatrixCache)*=0;
	(*crossTalkMatrixCache)+=(*crossTalkMatrix);
	(*crossTalkMatrix)*=0;
      }
    }      
  }

  sw.Stop();
  AliInfoF("timing: %e/%e real/cpu",sw.RealTime(),sw.CpuTime());
  //
}




void    AliTPCtracker::FilterOutlierClusters(){
  //
  // filter outlier clusters  
  //
  /*
    1.)..... booking part
    nSectors=72;
    nTimeBins=fParam->Get....
    TH2F hisTime("","", sector,0,sector, nTimeBins,0,nTimeBins);
    TH2F hisPadRow("","", sector,0,sector, nPadRows,0,nPadRows);
    2.) .... filling part
    .... cluster loop { hisTime.Fill(cluster->GetDetector(),cluster->GetTimeBin()); }
    
    3.) ...filtering part 
    sector loop { calculate median,mean80 and rms80 of the nclusters per time bin; calculate median,mean80 and rms80 of the nclusters per par row; .... export values to the debug streamers - to decide which threshold to be used... }
    
    sector loop
    { disable clusters in time bins > mean+ n rms80+ offsetTime disable clusters in padRow > mean+ n rms80+ offsetPadRow // how to dislable clusters? - new bit to introduce } 
    //
    4. Disabling clusters
    
  */
  
  //
  // 1.) booking part 
  // 
  //  AliTPCcalibDB *db=AliTPCcalibDB::Instance();
  Int_t nSectors=AliTPCROC::Instance()->GetNSectors(); 
  Int_t nTimeBins= 1100; // *Bug here - we should get NTimeBins from ALTRO - Parameters not relyable
  Int_t nPadRows=(AliTPCROC::Instance()->GetNRows(0) + AliTPCROC::Instance()->GetNRows(36));
  // parameters for filtering
  const Double_t nSigmaCut=9.;           // should be in recoParam ?
  const Double_t offsetTime=100;         // should be in RecoParam ?  -
  const Double_t offsetPadRow=300;       // should be in RecoParam ?
  const Double_t offsetTimeAccept=8;     // should be in RecoParam ?  - obtained as mean +1 rms in high IR pp
  TH2F hisTime("hisSectorTime","hisSectorTime", nSectors,0,nSectors, nTimeBins,0,nTimeBins);
  TH2F hisPadRow("hisSectorRow","hisSectorRow", nSectors,0,nSectors, nPadRows,0,nPadRows);
  //
  // 2.) Filling part -- loop over clusters
  //
  for (Int_t isector=0; isector<36; isector++){      // loop over sectors
    for (Int_t iside=0; iside<2; iside++){           // loop over sides A/C
      AliTPCtrackerSector &sector= (isector<18)?fInnerSec[isector%18]:fOuterSec[isector%18];
      Int_t nrows = sector.GetNRows();
      for (Int_t row = 0;row<nrows;row++){           // loop over rows       
        AliTPCtrackerRow&  tpcrow = sector[row];       
        Int_t ncl = tpcrow.GetN1();                  // number of clusters in the row
        if (iside>0) ncl=tpcrow.GetN2();
        for (Int_t i=0;i<ncl;i++) {  // loop over clusters
          AliTPCclusterMI *cluster= (iside>0)?(tpcrow.GetCluster2(i)):(tpcrow.GetCluster1(i));
          hisTime.Fill(cluster->GetDetector(),cluster->GetTimeBin()); 
          hisPadRow.Fill(cluster->GetDetector(),cluster->GetRow()); 
	}
      } 
    } 
  } 

  //
  // 3. Filtering part
  //
  TVectorD vecTime(nTimeBins);
  TVectorD vecPadRow(nPadRows);
  TVectorD vecMedianSectorTime(nSectors);
  TVectorD vecRMSSectorTime(nSectors);
  TVectorD vecMedianSectorTimeOut6(nSectors);
  TVectorD vecMedianSectorTimeOut9(nSectors);//
  TVectorD vecMedianSectorTimeOut(nSectors);//
  TVectorD vecMedianSectorPadRow(nSectors);
  TVectorD vecRMSSectorPadRow(nSectors);
  TVectorD vecMedianSectorPadRowOut6(nSectors);
  TVectorD vecMedianSectorPadRowOut9(nSectors);
  TVectorD vecMedianSectorPadRowOut(nSectors);
  TVectorD vecSectorOut6(nSectors);
  TVectorD vecSectorOut9(nSectors);
  TMatrixD matSectorCluster(nSectors,2);
  //
  // 3.a)  median, rms calculations for hisTime 
  //
  for (Int_t isec=0; isec<nSectors; isec++){
    vecMedianSectorTimeOut6[isec]=0;
    vecMedianSectorTimeOut9[isec]=0;
    for (Int_t itime=0; itime<nTimeBins; itime++){
      vecTime[itime]=hisTime.GetBinContent(isec+1, itime+1);      
    }
    Double_t median= TMath::Mean(nTimeBins,vecTime.GetMatrixArray());
    Double_t rms= TMath::RMS(nTimeBins,vecTime.GetMatrixArray()); 
    vecMedianSectorTime[isec]=median;
    vecRMSSectorTime[isec]=rms;
    if ((AliTPCReconstructor::StreamLevel()&kStreamFilterClusterInfo)>0) AliInfo(TString::Format("Sector TimeStat: %d\t%8.0f\t%8.0f",isec,median,rms).Data());
    //
    // declare outliers
    for (Int_t itime=0; itime<nTimeBins; itime++){
      Double_t entries= hisTime.GetBinContent(isec+1, itime+1); 
      if (entries>median+6.*rms+offsetTime) {
	vecMedianSectorTimeOut6[isec]+=1;
      }
      if (entries>median+9.*rms+offsetTime) {
	vecMedianSectorTimeOut9[isec]+=1;
      }
    }
  }    
  //
  // 3.b) median, rms calculations for hisPadRow
  // 
  for (Int_t isec=0; isec<nSectors; isec++){
    vecMedianSectorPadRowOut6[isec]=0;
    vecMedianSectorPadRowOut9[isec]=0;
    for (Int_t ipadrow=0; ipadrow<nPadRows; ipadrow++){
      vecPadRow[ipadrow]=hisPadRow.GetBinContent(isec+1, ipadrow+1);      
    }
    Int_t nPadRowsSector= AliTPCROC::Instance()->GetNRows(isec);
    Double_t median= TMath::Mean(nPadRowsSector,vecPadRow.GetMatrixArray());
    Double_t rms= TMath::RMS(nPadRowsSector,vecPadRow.GetMatrixArray());
    vecMedianSectorPadRow[isec]=median;
    vecRMSSectorPadRow[isec]=rms;
    if ((AliTPCReconstructor::StreamLevel()&kStreamFilterClusterInfo)>0) AliInfo(TString::Format("Sector PadRowStat: %d\t%8.0f\t%8.0f",isec,median,rms).Data());
    //
    // declare outliers
    for (Int_t ipadrow=0; ipadrow<nPadRows; ipadrow++){
      Double_t entries= hisPadRow.GetBinContent(isec+1, ipadrow+1);
      if (entries>median+6.*rms+offsetPadRow) {
        vecMedianSectorPadRowOut6[isec]+=1;
      }
      if (entries>median+9.*rms+offsetPadRow) {
        vecMedianSectorPadRowOut9[isec]+=1;
      }
    }
  }
  //
  // 3.c) filter outlier sectors
  //
  Double_t medianSectorTime = TMath::Median(nSectors, vecTime.GetMatrixArray());
  Double_t mean69SectorTime, rms69SectorTime=0;
  AliMathBase::EvaluateUni(nSectors,  vecTime.GetMatrixArray(), mean69SectorTime,rms69SectorTime,69);
  for (Int_t isec=0; isec<nSectors; isec++){
    vecSectorOut6[isec]=0;
    vecSectorOut9[isec]=0;
    matSectorCluster(isec,0)=0;
    matSectorCluster(isec,1)=0;
    if (TMath::Abs(vecMedianSectorTime[isec])>(mean69SectorTime+6.*(rms69SectorTime+ offsetTimeAccept))) {
      vecSectorOut6[isec]=1;
    }
    if (TMath::Abs(vecMedianSectorTime[isec])>(mean69SectorTime+9.*(rms69SectorTime+ offsetTimeAccept))){
      vecSectorOut9[isec]=1;
    }
  }
  // light version of export variable
  Int_t filteredSector= vecSectorOut9.Sum();                  // light version of export variable
  Int_t filteredSectorTime= vecMedianSectorTimeOut9.Sum();
  Int_t filteredSectorPadRow= vecMedianSectorPadRowOut9.Sum();
  if (fEvent) if (fEvent->GetHeader()){
    fEvent->GetHeader()->SetTPCNoiseFilterCounter(0,TMath::Min(filteredSector,255));
    fEvent->GetHeader()->SetTPCNoiseFilterCounter(1,TMath::Min(filteredSectorTime,255));
    fEvent->GetHeader()->SetTPCNoiseFilterCounter(2,TMath::Min(filteredSectorPadRow,255));
  }
 
  //
  // 4. Disabling clusters in outlier layers
  //
  Int_t counterAll=0;
  Int_t counterOut=0;
  for (Int_t isector=0; isector<36; isector++){      // loop over sectors
    for (Int_t iside=0; iside<2; iside++){           // loop over sides A/C
      AliTPCtrackerSector &sector= (isector<18)?fInnerSec[isector%18]:fOuterSec[isector%18];
      Int_t nrows = sector.GetNRows();
      for (Int_t row = 0;row<nrows;row++){           // loop over rows       
        AliTPCtrackerRow&  tpcrow = sector[row];       
        Int_t ncl = tpcrow.GetN1();                  // number of clusters in the row
        if (iside>0) ncl=tpcrow.GetN2();
        for (Int_t i=0;i<ncl;i++) {  // loop over clusters
          AliTPCclusterMI *cluster= (iside>0)?(tpcrow.GetCluster2(i)):(tpcrow.GetCluster1(i));
	  Double_t medianTime=vecMedianSectorTime[cluster->GetDetector()];
	  Double_t medianPadRow=vecMedianSectorPadRow[cluster->GetDetector()];
	  Double_t rmsTime=vecRMSSectorTime[cluster->GetDetector()];
	  Double_t rmsPadRow=vecRMSSectorPadRow[cluster->GetDetector()];
	  Int_t entriesPadRow=hisPadRow.GetBinContent(cluster->GetDetector()+1, cluster->GetRow()+1);
	  Int_t entriesTime=hisTime.GetBinContent(cluster->GetDetector()+1, cluster->GetTimeBin()+1);
	  Bool_t isOut=kFALSE;
	  if (vecSectorOut9[cluster->GetDetector()]>0.5) {
	    isOut=kTRUE;
	  }
	  
	  if (entriesTime>medianTime+nSigmaCut*rmsTime+offsetTime) {
	    isOut=kTRUE;
	    vecMedianSectorTimeOut[cluster->GetDetector()]++;
	  }
	  if (entriesPadRow>medianPadRow+nSigmaCut*rmsPadRow+offsetPadRow) {
	    isOut=kTRUE;
	    vecMedianSectorPadRowOut[cluster->GetDetector()]++;
	  }
	  counterAll++;
	  matSectorCluster(cluster->GetDetector(),0)+=1;
	  if (isOut){
	    cluster->Disable();
	    counterOut++;
	    matSectorCluster(cluster->GetDetector(),1)+=1;
	  }
	}
      }
    }
  }
  for (Int_t isec=0; isec<nSectors; isec++){
    if ((AliTPCReconstructor::StreamLevel()&kStreamFilterClusterInfo)>0) AliInfo(TString::Format("Sector Stat: %d\t%8.0f\t%8.0f",isec,matSectorCluster(isec,1),matSectorCluster(isec,0)).Data());
  }
  //
  // dump info to streamer - for later tuning of cuts
  //
  if ((AliTPCReconstructor::StreamLevel()&kStreamFilterClusterInfo)>0) {  // stream TPC data ouliers filtering infomation
    AliLog::Flush();
    AliInfo(TString::Format("Cluster counter: (%d/%d) (Filtered/All)",counterOut,counterAll).Data());
    for (Int_t iSec=0; iSec<nSectors; iSec++){
      if (vecSectorOut9[iSec]>0 ||  matSectorCluster(iSec,1)>0) {
	AliInfo(TString::Format("Filtered sector\t%d",iSec).Data());
	Double_t vecMedTime =TMath::Median(72,vecMedianSectorTime.GetMatrixArray());
	Double_t vecMedPadRow =TMath::Median(72,vecMedianSectorPadRow.GetMatrixArray());
	Double_t vecMedCluster=(counterAll-counterOut)/72;
	AliInfo(TString::Format("VecMedianSectorTime\t(%4.4f/%4.4f/%4.4f)",       vecMedianSectorTimeOut[iSec],vecMedianSectorTime[iSec],vecMedTime).Data());
	AliInfo(TString::Format("VecMedianSectorPadRow\t(%4.4f/%4.4f/%4.4f)",     vecMedianSectorPadRowOut[iSec],vecMedianSectorPadRow[iSec],vecMedPadRow).Data());
	AliInfo(TString::Format("MatSectorCluster\t(%4.4f/%4.4f/%4.4f)\n",          matSectorCluster(iSec,1), matSectorCluster(iSec,0),  vecMedCluster).Data());
	AliLog::Flush();
      }
    }
    AliLog::Flush();
    Int_t eventNr = fEvent->GetEventNumberInFile();
    (*fDebugStreamer)<<"filterClusterInfo"<<
      // minimal set variables for the ESDevent
      "eventNr="<<eventNr<<
      "counterAll="<<counterAll<<
      "counterOut="<<counterOut<<
      "matSectotCluster.="<<&matSectorCluster<<                   // 
      //
      "filteredSector="<<filteredSector<<                        //  counter filtered sectors                   
      "filteredSectorTime="<<filteredSectorTime<<                //  counter filtered time bins
      "filteredSectorPadRow="<<filteredSectorPadRow<<            //  counter filtered pad-rows
      // per sector outlier information
      "medianSectorTime="<<medianSectorTime<<                    // median number of clusters per sector/timebin
      "mean69SectorTime="<<mean69SectorTime<<                    // LTM statistic  mean of clusters per sector/timebin
      "rms69SectorTime="<<rms69SectorTime<<                      // LTM statistic  RMS of clusters per sector/timebin
      "vecSectorOut6.="<<&vecSectorOut6<<                        // flag array sector  - 6 sigma +accept margin outlier
      "vecSectorOut9.="<<&vecSectorOut9<<                        // flag array sector  - 9 sigma + accept margin outlier
      // per sector/timebin outlier detection
      "vecMedianSectorTime.="<<&vecMedianSectorTime<<
      "vecRMSSectorTime.="<<&vecRMSSectorTime<<
      "vecMedianSectorTimeOut6.="<<&vecMedianSectorTimeOut6<<
      "vecMedianSectorTimeOut9.="<<&vecMedianSectorTimeOut9<<
      "vecMedianSectorTimeOut0.="<<&vecMedianSectorTimeOut<<
      // per sector/pad-row outlier detection
      "vecMedianSectorPadRow.="<<&vecMedianSectorPadRow<<
      "vecRMSSectorPadRow.="<<&vecRMSSectorPadRow<<
      "vecMedianSectorPadRowOut6.="<<&vecMedianSectorPadRowOut6<<
      "vecMedianSectorPadRowOut9.="<<&vecMedianSectorPadRowOut9<<
      "vecMedianSectorPadRowOut9.="<<&vecMedianSectorPadRowOut<<
      "\n";
    ((*fDebugStreamer)<<"filterClusterInfo").GetTree()->Write();
    fDebugStreamer->GetFile()->Flush();
  }
}

void AliTPCtracker::UnloadClusters()
{
  //
  // unload clusters from the memory
  //
  Int_t nrows = fOuterSec->GetNRows();
  for (Int_t sec = 0;sec<fkNOS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCtrackerRow*  tpcrow = &(fOuterSec[sec%fkNOS][row]);
      //      if (tpcrow){
      //	if (tpcrow->fClusters1) delete []tpcrow->fClusters1; 
      //	if (tpcrow->fClusters2) delete []tpcrow->fClusters2; 
      //}
      tpcrow->ResetClusters();
    }
  //
  nrows = fInnerSec->GetNRows();
  for (Int_t sec = 0;sec<fkNIS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCtrackerRow*  tpcrow = &(fInnerSec[sec%fkNIS][row]);
      //if (tpcrow){
      //	if (tpcrow->fClusters1) delete []tpcrow->fClusters1; 
      //if (tpcrow->fClusters2) delete []tpcrow->fClusters2; 
      //}
      tpcrow->ResetClusters();
    }

  fNClusters = 0;
  return ;
}

void AliTPCtracker::FillClusterArray(TObjArray* array) const{
  //
  // Filling cluster to the array - For visualization purposes
  //
  Int_t nrows=0;
  nrows = fOuterSec->GetNRows();
  for (Int_t sec = 0;sec<fkNOS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCtrackerRow*  tpcrow = &(fOuterSec[sec%fkNOS][row]);
      if (!tpcrow) continue;
      for (Int_t icl = 0;icl<tpcrow->GetN();icl++){
	array->AddLast((TObject*)((*tpcrow)[icl]));
      }
    } 
  nrows = fInnerSec->GetNRows();
  for (Int_t sec = 0;sec<fkNIS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCtrackerRow*  tpcrow = &(fInnerSec[sec%fkNIS][row]);
      if (!tpcrow) continue;
      for (Int_t icl = 0;icl<tpcrow->GetN();icl++){
	array->AddLast((TObject*)(*tpcrow)[icl]);
      }
    }
}


void AliTPCtracker::Transform(AliTPCclusterMI * cluster){
  //
  // transformation
  //
  const double kMaxY2X = AliTPCTransform::GetMaxY2X();  // tg of sector angular span
  const double kSinSect = TMath::Sin(TMath::Pi()/9), kCosSect = TMath::Cos(TMath::Pi()/9);
  //
  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform() ;
  if (!transform) {
    AliFatal("Tranformations not in calibDB");
    return;
  }
  //  if (!transform->GetCurrentRecoParam()) transform->SetCurrentRecoParam((AliTPCRecoParam*)AliTPCReconstructor::GetRecoParam());
  Double_t x[3]={static_cast<Double_t>(cluster->GetRow()),static_cast<Double_t>(cluster->GetPad()),static_cast<Double_t>(cluster->GetTimeBin())};
  Int_t idROC = cluster->GetDetector();
  //  transform->Transform(x,&idROC,0,1);
  transform->Transform(x,&idROC,0,cluster->GetLabel(0));
  const float* clCorr = transform->GetLastMapCorrection();
  const float* clCorrRef = transform->GetLastMapCorrectionRef();
  //
  cluster->SetDistortions(clCorr[0]-clCorrRef[0],
			  clCorr[1]-clCorrRef[1],
			  clCorr[2]-clCorrRef[2]); // memorize distortions (difference to reference one)
  // store the dispersion difference
  cluster->SetDistortionDispersion(clCorr[3]); // ref error is already subtracted
  //
  cluster->SetX(x[0]);
  cluster->SetY(x[1]);
  cluster->SetZ(x[2]);
  // in debug mode  check the transformation
  //
  if ((AliTPCReconstructor::StreamLevel()&kStreamTransform)>0) { 
    Float_t gx[3];
    cluster->GetGlobalXYZ(gx);
    Int_t event = (fEvent==NULL)? 0: fEvent->GetEventNumberInFile();
    TTreeSRedirector &cstream = *fDebugStreamer;
    Int_t timeStamp=transform->GetCurrentTimeStamp();
    float* nCclCorr = (float*)transform->GetLastMapCorrection();  
    float* nCclCorrRef = (float*)transform->GetLastMapCorrectionRef();
    transform->SetDebugStreamer(fDebugStreamer);

    cstream<<"Transform"<<  // needed for debugging of the cluster transformation, resp. used for later visualization 
      "event="<<event<<
      "timeStamp="<<timeStamp<<
      "x0="<<x[0]<<
      "x1="<<x[1]<<
      "x2="<<x[2]<<
      "gx0="<<gx[0]<<
      "gx1="<<gx[1]<<
      "gx2="<<gx[2]<<
      "dx="<<nCclCorr[0]<<
      "dy="<<nCclCorr[1]<<
      "dz="<<nCclCorr[2]<<
      "dxRef="<<nCclCorrRef[0]<<
      "dyRef="<<nCclCorrRef[1]<<
      "dzRef="<<nCclCorrRef[2]<<
      "Cl.="<<cluster<<
      "\n"; 
  }
  // The old stuff:
  //
  // 
  //
  //if (!fkParam->IsGeoRead()) fkParam->ReadGeoMatrices();
  if (AliTPCReconstructor::GetRecoParam()->GetUseSectorAlignment() && (!calibDB->HasAlignmentOCDB())){
    TGeoHMatrix  *mat = fkParam->GetClusterMatrix(cluster->GetDetector());
    //TGeoHMatrix  mat;
    Double_t pos[3]= {cluster->GetX(),cluster->GetY(),cluster->GetZ()};
    Double_t posC[3]={cluster->GetX(),cluster->GetY(),cluster->GetZ()};
    if (mat) mat->LocalToMaster(pos,posC);
    else{
      // chack Loading of Geo matrices from GeoManager - TEMPORARY FIX
    }
    cluster->SetX(posC[0]);
    cluster->SetY(posC[1]);
    cluster->SetZ(posC[2]);
  }
}

void  AliTPCtracker::ApplyXtalkCorrection(){
  //
  // ApplyXtalk correction 
  // Loop over all clusters
  //      add to each cluster signal corresponding to common Xtalk mode for given time bin at given wire segment
  // cluster loop
  TStopwatch sw;
  sw.Start();
  for (Int_t isector=0; isector<36; isector++){  //loop tracking sectors
    for (Int_t iside=0; iside<2; iside++){       // loop over sides A/C
      AliTPCtrackerSector &sector= (isector<18)?fInnerSec[isector%18]:fOuterSec[isector%18];
      Int_t nrows     = sector.GetNRows();       
      for (Int_t row = 0;row<nrows;row++){           // loop over rows       
	AliTPCtrackerRow&  tpcrow = sector[row];     // row object   
	Int_t ncl = tpcrow.GetN1();                  // number of clusters in the row
	if (iside>0) ncl=tpcrow.GetN2();
	Int_t xSector=0;    // sector number in the TPC convention 0-72
	if (isector<18){  //if IROC
	  xSector=isector+(iside>0)*18;
	}else{
	  xSector=isector+18;  // isector -18 +36   
	  if (iside>0) xSector+=18;
	}	
	TMatrixD &crossTalkMatrix= *((TMatrixD*)fCrossTalkSignalArray->At(xSector));
	Int_t wireSegmentID     = fkParam->GetWireSegment(xSector,row);
	for (Int_t i=0;i<ncl;i++) {
	  AliTPCclusterMI *cluster= (iside>0)?(tpcrow.GetCluster2(i)):(tpcrow.GetCluster1(i));
	  Int_t iTimeBin=TMath::Nint(cluster->GetTimeBin());
	  Double_t xTalk= crossTalkMatrix[wireSegmentID][iTimeBin];
	  cluster->SetMax(cluster->GetMax()+xTalk);
	  const Double_t kDummy=4;
	  Double_t sumxTalk=xTalk*kDummy; // should be calculated via time response function
	  cluster->SetQ(cluster->GetQ()+sumxTalk);


          if ((AliTPCReconstructor::StreamLevel()&kStreamXtalk)>0) {  // flag: stream crosstalk correctio as applied to cluster
            TTreeSRedirector &cstream = *fDebugStreamer;
            if (gRandom->Rndm() > 0.){
              cstream<<"Xtalk"<<
                "isector=" << isector <<               // sector [0,36]
                "iside=" << iside <<                   // side A or C
                "row=" << row <<                       // padrow
                "i=" << i <<                           // index of the cluster 
                "xSector=" << xSector <<               // sector [0,72] 
                "wireSegmentID=" << wireSegmentID <<   // anode wire segment id [0,10] 
                "iTimeBin=" << iTimeBin <<             // timebin of the corrected cluster 
                "xTalk=" << xTalk <<                   // Xtalk contribution added to Qmax
                "sumxTalk=" << sumxTalk <<             // Xtalk contribution added to Qtot (roughly 3*Xtalk) 
                "cluster.=" << cluster <<              // corrected cluster object 
                "\n";
            }
          }// dump the results to the debug streamer if in debug mode
	}
      }
    }
  }
  sw.Stop();
  AliInfoF("timing: %e/%e real/cpu",sw.RealTime(),sw.CpuTime());
  //
}

void  AliTPCtracker::ApplyTailCancellation(){
  //
  // Correct the cluster charge for the ion tail effect 
  // The TimeResponse function accessed via  AliTPCcalibDB (TPC/Calib/IonTail)
  //
  TStopwatch sw;
  sw.Start();
  // Retrieve
  TObjArray *ionTailArr = (TObjArray*)AliTPCcalibDB::Instance()->GetIonTailArray();
  if (!ionTailArr) {AliFatal("TPC - Missing IonTail OCDB object");}
  TObject *rocFactorIROC  = ionTailArr->FindObject("factorIROC");
  TObject *rocFactorOROC  = ionTailArr->FindObject("factorOROC");   
  Float_t factorIROC      = (atof(rocFactorIROC->GetTitle()));
  Float_t factorOROC      = (atof(rocFactorOROC->GetTitle()));

  // find the number of clusters for the whole TPC (nclALL)
  Int_t nclALL=0;
  for (Int_t isector=0; isector<36; isector++){
    AliTPCtrackerSector &sector= (isector<18)?fInnerSec[isector%18]:fOuterSec[isector%18];
    nclALL += sector.GetNClInSector(0);
    nclALL += sector.GetNClInSector(1);
  }

  //RS changed from heap allocation in the loop to stack allocation
  TGraphErrors * graphRes[20]; 
  Float_t        indexAmpGraphs[20];      
  //  AliTPCclusterMI* rowClusterArray[kMaxClusterPerRow]; // caches clusters for each row  // RS avoid trashing the heap

  // start looping over all clusters 
  for (Int_t iside=0; iside<2; iside++){    // loop over sides
    //
    //
    for (Int_t secType=0; secType<2; secType++){  //loop over inner or outer sector
      // cache experimantal tuning factor for the different chamber type 
      const Float_t ampfactor = (secType==0)?factorIROC:factorOROC;
//       std::cout << " ampfactor = " << ampfactor << std::endl;
      //
      for (Int_t sec = 0;sec<fkNOS;sec++){        //loop overs sectors
        //
        //
        // Cache time response functions and their positons to COG of the cluster       
	// TGraphErrors ** graphRes   = new TGraphErrors *[20]; // RS avoid heap allocations if stack can be used
        // Float_t * indexAmpGraphs   = new Float_t[20];        // RS Moved outside of the loop
	memset(graphRes,0,20*sizeof(TGraphErrors*));
	memset(indexAmpGraphs,0,20*sizeof(float));
        //for (Int_t icache=0; icache<20; icache++)  //RS
        //{
        //  graphRes[icache]       = NULL;
        //  indexAmpGraphs[icache] = 0;
	// }
        /////////////////////////////  --> position fo sie loop
        if (!AliTPCcalibDB::Instance()->GetTailcancelationGraphs(sec+36*secType+18*iside,graphRes,indexAmpGraphs))
        {
          continue;
        }

        // set time range from graph
        const Int_t timeRangeMax=graphRes[0]?graphRes[0]->GetN():600;
//         std::cout << " timeRangeMax = " << timeRangeMax << std::endl;

        AliTPCtrackerSector &sector= (secType==0)?fInnerSec[sec]:fOuterSec[sec];  
        Int_t nrows     = sector.GetNRows();                                       // number of rows
        Int_t nclSector = sector.GetNClInSector(iside);                            // ncl per sector to be used for debugging

        for (Int_t row = 0;row<nrows;row++){           // loop over rows

          AliTPCtrackerRow&  tpcrow = sector[row];     // row object   
          Int_t ncl = tpcrow.GetN1();                  // number of clusters in the row
          if (iside>0) ncl=tpcrow.GetN2();

          // Order clusters in time for the proper correction of ion tail
          Float_t qTotArray[ncl];          // arrays to be filled with modified Qtot and Qmax values in order to avoid float->int conversion  
          Float_t qMaxArray[ncl];
          Int_t sortedClusterIndex[ncl];
          Float_t sortedClusterTimeBin[ncl];
          //TObjArray *rowClusterArray = new TObjArray(ncl);  // cache clusters for each row  // RS avoid trashing the heap
          AliTPCclusterMI* rowClusterArray[ncl]; // caches clusters for each row  // RS avoid trashing the heap 
          //  memset(rowClusterArray,0,sizeof(AliTPCclusterMI*)*ncl);  //.Clear();
          //if (rowClusterArray.GetSize()<ncl) rowClusterArray.Expand(ncl);
          for (Int_t i=0;i<ncl;i++) 
          {
            qTotArray[i]=0;
            qMaxArray[i]=0;
            sortedClusterIndex[i]=i;
            AliTPCclusterMI *rowcl= (iside>0)?(tpcrow.GetCluster2(i)):(tpcrow.GetCluster1(i));
            rowClusterArray[i] = rowcl;
	    //if (rowcl) {
            //  rowClusterArray.AddAt(rowcl,i);
            //} else {
            //  rowClusterArray.RemoveAt(i);
            //}
            // Fill the timebin info to the array in order to sort wrt tb
            if (!rowcl) {
	      sortedClusterTimeBin[i]=0.0;
	    } else {
	      sortedClusterTimeBin[i] = rowcl->GetTimeBin();
	    }

          } 
	  TMath::Sort(ncl,sortedClusterTimeBin,sortedClusterIndex,kFALSE);       // sort clusters in time
     
          // Main cluster correction loops over clusters
          for (Int_t icl0=0; icl0<ncl;icl0++){    // first loop over clusters

            AliTPCclusterMI *cl0= rowClusterArray[sortedClusterIndex[icl0]]; //RS static_cast<AliTPCclusterMI*>(rowClusterArray.At(sortedClusterIndex[icl0]));
            
            if (!cl0) continue;
            Int_t nclPad=0;
            //for (Int_t icl1=0; icl1<ncl;icl1++){  // second loop over clusters	   
	    // RS: time increases with index since sorted -> cl0->GetTimeBin()>cl1->GetTimeBin() means that icl0>icl1
            for (Int_t icl1=0; icl1<icl0;icl1++){  // second loop over clusters
              AliTPCclusterMI *cl1= rowClusterArray[sortedClusterIndex[icl1]];//RS static_cast<AliTPCclusterMI*>(rowClusterArray.At(sortedClusterIndex[icl1]));
	      if (!cl1) continue;
	      // RS no needed with proper loop organization
	      //if (cl0->GetTimeBin()<= cl1->GetTimeBin()) continue;               // no contibution to the tail if later

	      int dpad = TMath::Abs(cl0->GetPad()-cl1->GetPad());
	      if (dpad>4) continue;           // no contribution if far away in pad direction

	      // RS no point in iterating further with sorted clusters once large distance reached
              if (cl0->GetTimeBin()-cl1->GetTimeBin()>=timeRangeMax) continue; // out of the range of response function

	      // RS: what about dpad=4?
              if (dpad<4) nclPad++;           // count ncl for every pad for debugging
            
              // Get the correction values for Qmax and Qtot and find total correction for a given cluster
              Double_t ionTailMax=0.;  
              Double_t ionTailTotal=0.;  
              GetTailValue(ampfactor,ionTailMax,ionTailTotal,graphRes,indexAmpGraphs,cl0,cl1);
              ionTailMax=TMath::Abs(ionTailMax);
              ionTailTotal=TMath::Abs(ionTailTotal);
              qTotArray[icl0]+=ionTailTotal;
              qMaxArray[icl0]+=ionTailMax;

              // Dump some info for debugging while clusters are being corrected
              if ((AliTPCReconstructor::StreamLevel()&kStreamIonTail)>0) {  // flag: stream ion tail correction  as applied to cluster
                TTreeSRedirector &cstream = *fDebugStreamer;
                if (gRandom->Rndm() > 0.999){
                  cstream<<"IonTail"<<
                      "cl0.="         <<cl0          <<   // cluster 0 (to be corrected)
                      "cl1.="         <<cl1          <<   // cluster 1 (previous cluster)
                      "ionTailTotal=" <<ionTailTotal <<   // ion Tail from cluster 1 contribution to cluster0
                      "ionTailMax="   <<ionTailMax   <<   // ion Tail from cluster 1 contribution to cluster0 
                      "\n";
                }
              }// dump the results to the debug streamer if in debug mode
            
            }//end of second loop over clusters
            
            // Set corrected values of the corrected cluster          
            cl0->SetQ(TMath::Nint(Float_t(cl0->GetQ())+Float_t(qTotArray[icl0])));
            cl0->SetMax(TMath::Nint(Float_t(cl0->GetMax())+qMaxArray[icl0]));
          
            // Dump some info for debugging after clusters are corrected 
            if ((AliTPCReconstructor::StreamLevel()&kStreamIonTail)>0) {
              TTreeSRedirector &cstream = *fDebugStreamer;
              if (gRandom->Rndm() > 0.999){
              cstream<<"IonTailCorrected"<<
                  "cl0.="                     << cl0              <<   // cluster 0 with huge Qmax
                  "ionTailTotalPerCluster="   << qTotArray[icl0]  <<
                  "ionTailMaxPerCluster="     << qMaxArray[icl0]  <<
                  "nclALL="                   << nclALL           <<
                  "nclSector="                << nclSector        <<
                  "nclRow="                   << ncl              <<
                  "nclPad="                   << nclPad           <<
                  "row="                      << row              <<
                  "sector="                   << sec              <<
                  "icl0="                     << icl0             <<
                  "\n";
              }
            }// dump the results to the debug streamer if in debug mode
          
          }//end of first loop over cluster
          // delete rowClusterArray; // RS was moved to stack allocation
        }//end of loop over rows
        for (int i=0; i<20; i++) delete graphRes[i];
	//        delete [] graphRes; //RS was changed to stack allocation
	//        delete [] indexAmpGraphs;
      
      }//end of loop over sectors
    }//end of loop over IROC/OROC
  }// end of side loop
  sw.Stop();
  AliInfoF("timing: %e/%e real/cpu",sw.RealTime(),sw.CpuTime());
  //
}
//_____________________________________________________________________________
void AliTPCtracker::GetTailValue(Float_t ampfactor,Double_t &ionTailMax, Double_t &ionTailTotal,TGraphErrors **graphRes,Float_t *indexAmpGraphs,AliTPCclusterMI *cl0,AliTPCclusterMI *cl1){

  //
  // Function in order to calculate the amount of the correction to be added for a given cluster, return values are ionTailTaoltal and ionTailMax
  // Parameters:
  // cl0 -  cluster to be modified
  // cl1 -  source cluster ion tail of this cluster will be added to the cl0 (accroding time and pad response function)
  // 
  const float kMinPRF       = 0.5f;                          // minimal PRF width
  ionTailTotal              = 0.;                            // correction value to be added to Qtot of cl0
  ionTailMax                = 0.;                            // correction value to be added to Qmax of cl0

  Float_t qTot0             =  cl0->GetQ();                  // cl0 Qtot info
  Float_t qTot1             =  cl1->GetQ();                  // cl1 Qtot info
  Int_t sectorPad           =  cl1->GetDetector();           // sector number
  Int_t padcl0              =  TMath::Nint(cl0->GetPad());   // pad0
  Int_t padcl1              =  TMath::Nint(cl1->GetPad());   // pad1
  Float_t padWidth          = (sectorPad < 36)?0.4:0.6;      // pad width in cm
  const Int_t deltaTimebin  =  TMath::Nint(TMath::Abs(cl1->GetTimeBin()-cl0->GetTimeBin()))+12;  //distance between pads of cl1 and cl0 increased by 12 bins
  float rmsPad1I            = (cl1->GetSigmaY2()==0)?0.5f/kMinPRF:(0.5f*padWidth/sqrtf(cl1->GetSigmaY2()));
  float rmsPad0I            = (cl0->GetSigmaY2()==0)?0.5f/kMinPRF:(0.5f*padWidth/sqrtf(cl0->GetSigmaY2()));
  
  // RS avoid useless calculations
  //Double_t sumAmp1=0.;
  //for (Int_t idelta =-2; idelta<=2;idelta++){
  //  sumAmp1+=TMath::Exp(-idelta*idelta*rmsPad1I);
  //}
  // Double_t sumAmp0=0.;
  //for (Int_t idelta =-2; idelta<=2;idelta++){
  //  sumAmp0+=TMath::Exp(-idelta*idelta*rmsPad0I));
  //}

  float tmp = expf(-rmsPad1I);
  float sumAmp1 = 1.f/(1.f+2.f*tmp*(1.f+tmp*tmp*tmp));
  tmp = expf(-rmsPad0I);
  float sumAmp0 = 1.f/(1.f+2.f*tmp*(1.f+tmp*tmp*tmp));

  // Apply the correction  -->   cl1 corrects cl0 (loop over cl1's pads and find which pads of cl0 are going to be corrected)
  Int_t padScan=2;      // +-2 pad-timebin window will be scanned
  for (Int_t ipad1=padcl1-padScan; ipad1<=padcl1+padScan; ipad1++) {
    //
    //
    Float_t deltaPad1  = TMath::Abs(cl1->GetPad()-(Float_t)ipad1);
    Float_t amp1       = TMath::Exp(-(deltaPad1*deltaPad1)*rmsPad1I)*sumAmp1;  // normalized pad response function
    Float_t qTotPad1   = amp1*qTot1;                                           // used as a factor to multipliy the response function
      
    // find closest value of cl1 to COG (among the time response functions' amplitude array --> to select proper t.r.f.)
    Int_t ampIndex = 0;
    Float_t diffAmp  = TMath::Abs(deltaPad1-indexAmpGraphs[0]);
    for (Int_t j=0;j<20;j++) {
      if (diffAmp > TMath::Abs(deltaPad1-indexAmpGraphs[j]) && indexAmpGraphs[j]!=0)
        {
          diffAmp  = TMath::Abs(deltaPad1-indexAmpGraphs[j]);
          ampIndex = j;
        }
    }
    if (!graphRes[ampIndex]) continue;
    if (deltaTimebin+2 >= graphRes[ampIndex]->GetN()) continue;
    if (graphRes[ampIndex]->GetY()[deltaTimebin+2]>=0) continue;
     
    for (Int_t ipad0=padcl0-padScan; ipad0<=padcl0+padScan; ipad0++) {
      //
      //
      if (ipad1!=ipad0) continue;                                     // check if ipad1 channel sees ipad0 channel, if not no correction to be applied.
      
      Float_t deltaPad0  = TMath::Abs(cl0->GetPad()-(Float_t)ipad0);
      Float_t amp0       = expf(-(deltaPad0*deltaPad0)*rmsPad0I)*sumAmp0;  // normalized pad resp function
      Float_t qMaxPad0   = amp0*qTot0;
           
      // Add 5 timebin range contribution around the max peak (-+2 tb window)
      for (Int_t itb=deltaTimebin-2; itb<=deltaTimebin+2; itb++) {

        if (itb<0) continue; 
        if (itb>=graphRes[ampIndex]->GetN()) continue;
       
        // calculate contribution to qTot
        Float_t tailCorr =  TMath::Abs((qTotPad1*ampfactor)*(graphRes[ampIndex])->GetY()[itb]);
        if (ipad1!=padcl0) { 
          ionTailTotal += TMath::Min(qMaxPad0,tailCorr);   // for side pad
        } else {             
          ionTailTotal += tailCorr;                        // for center pad
        }
        // calculate contribution to qMax
        if (itb == deltaTimebin && ipad1 == padcl0) ionTailMax += tailCorr;   
        
      } // end of tb correction loop which is applied over 5 tb range

    } // end of cl0 loop
  } // end of cl1 loop
  
}

//_____________________________________________________________________________
Int_t AliTPCtracker::LoadOuterSectors() {
  //-----------------------------------------------------------------
  // This function fills outer TPC sectors with clusters.
  //-----------------------------------------------------------------
  Int_t nrows = fOuterSec->GetNRows();
  UInt_t index=0;
  for (Int_t sec = 0;sec<fkNOS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCtrackerRow*  tpcrow = &(fOuterSec[sec%fkNOS][row]);  
      Int_t sec2 = sec+2*fkNIS;
      //left
      Int_t ncl = tpcrow->GetN1();
      while (ncl--) {
	AliTPCclusterMI *c= (tpcrow->GetCluster1(ncl));
	index=(((sec2<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //right
      ncl = tpcrow->GetN2();
      while (ncl--) {
	AliTPCclusterMI *c= (tpcrow->GetCluster2(ncl));
	index=((((sec2+fkNOS)<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //
      // write indexes for fast acces
      //
      for (Int_t i=510;i--;) tpcrow->SetFastCluster(i,-1);
      for (Int_t i=0;i<tpcrow->GetN();i++){
        Int_t zi = Int_t((*tpcrow)[i]->GetZ()+255.);
	tpcrow->SetFastCluster(zi,i);  // write index
      }
      Int_t last = 0;
      for (Int_t i=0;i<510;i++){
	if (tpcrow->GetFastCluster(i)<0)
	  tpcrow->SetFastCluster(i,last);
	else
	  last = tpcrow->GetFastCluster(i);
      }
    }  
  fN=fkNOS;
  fSectors=fOuterSec;
  return 0;
}


//_____________________________________________________________________________
Int_t  AliTPCtracker::LoadInnerSectors() {
  //-----------------------------------------------------------------
  // This function fills inner TPC sectors with clusters.
  //-----------------------------------------------------------------
  Int_t nrows = fInnerSec->GetNRows();
  UInt_t index=0;
  for (Int_t sec = 0;sec<fkNIS;sec++)
    for (Int_t row = 0;row<nrows;row++){
      AliTPCtrackerRow*  tpcrow = &(fInnerSec[sec%fkNIS][row]);
      //
      //left
      Int_t ncl = tpcrow->GetN1();
      while (ncl--) {
	AliTPCclusterMI *c= (tpcrow->GetCluster1(ncl));
	index=(((sec<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //right
      ncl = tpcrow->GetN2();
      while (ncl--) {
	AliTPCclusterMI *c= (tpcrow->GetCluster2(ncl));
	index=((((sec+fkNIS)<<8)+row)<<16)+ncl;
	tpcrow->InsertCluster(c,index);
      }
      //
      // write indexes for fast acces
      //
      for (Int_t i=0;i<510;i++)
	tpcrow->SetFastCluster(i,-1);
      for (Int_t i=0;i<tpcrow->GetN();i++){
        Int_t zi = Int_t((*tpcrow)[i]->GetZ()+255.);
	tpcrow->SetFastCluster(zi,i);  // write index
      }
      Int_t last = 0;
      for (Int_t i=0;i<510;i++){
	if (tpcrow->GetFastCluster(i)<0)
	  tpcrow->SetFastCluster(i,last);
	else
	  last = tpcrow->GetFastCluster(i);
      }

    }  
   
  fN=fkNIS;
  fSectors=fInnerSec;
  return 0;
}



//_________________________________________________________________________
AliTPCclusterMI *AliTPCtracker::GetClusterMI(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  if (index<0) return 0; // no cluster
  Int_t sec=(index&0xff000000)>>24; 
  Int_t row=(index&0x00ff0000)>>16; 
  Int_t ncl=(index&0x00007fff)>>00;

  const AliTPCtrackerRow * tpcrow=0;
  TClonesArray * clrow =0;

  if (sec<0 || sec>=fkNIS*4) {
    AliWarning(Form("Wrong sector %d",sec));
    return 0x0;
  }

  if (sec<fkNIS*2){
    AliTPCtrackerSector& tracksec = fInnerSec[sec%fkNIS];
    if (tracksec.GetNRows()<=row) return 0;
    tpcrow = &(tracksec[row]);
    if (tpcrow==0) return 0;

    if (sec<fkNIS) {
      if (tpcrow->GetN1()<=ncl) return 0;
      clrow = tpcrow->GetClusters1();
    }
    else {
      if (tpcrow->GetN2()<=ncl) return 0;
      clrow = tpcrow->GetClusters2();
    }
  }
  else {
    AliTPCtrackerSector& tracksec = fOuterSec[(sec-fkNIS*2)%fkNOS];
    if (tracksec.GetNRows()<=row) return 0;
    tpcrow = &(tracksec[row]);
    if (tpcrow==0) return 0;

    if (sec-2*fkNIS<fkNOS) {
      if (tpcrow->GetN1()<=ncl) return 0;
      clrow = tpcrow->GetClusters1();
    }
    else {
      if (tpcrow->GetN2()<=ncl) return 0;
      clrow = tpcrow->GetClusters2();
    }
  }

  return (AliTPCclusterMI*)clrow->At(ncl);
  
}



Int_t AliTPCtracker::FollowToNext(AliTPCseed& t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  //
  const double kRoadY = 1., kRoadZ = 1.;
  Double_t  x= GetXrow(nr), ymax=GetMaxY(nr);
  //
  //
  AliTPCclusterMI *cl=0;
  Int_t tpcindex= t.GetClusterIndex2(nr);
  //
  // update current shape info every 5 pad-row
  //  if ( (nr%5==0) || t.GetNumberOfClusters()<2 || (t.fCurrentSigmaY2<0.0001) ){
    GetShape(&t,nr);    
    //}
  //  
  if (fIteration>0 && tpcindex>=-1){  //if we have already clusters 
    //        
    if (tpcindex==-1) return 0; //track in dead zone
    if (tpcindex >= 0){     //
      cl = t.GetClusterPointer(nr);         // cluster might not be there during track finding, but is attached in refits
      if (cl==0) cl = GetClusterMI(tpcindex); 
      t.SetCurrentClusterIndex1(tpcindex); 
    }
    if (cl){      
      Int_t relativesector = ((tpcindex&0xff000000)>>24)%18;  // if previously accepted cluster in different sector
      Float_t angle = relativesector*fSectors->GetAlpha()+fSectors->GetAlphaShift();
      //
      if (angle<-TMath::Pi()) angle += 2*TMath::Pi();
      if (angle>=TMath::Pi()) angle -= 2*TMath::Pi();
      
      if (TMath::Abs(angle-t.GetAlpha())>0.001){
	Double_t rotation = angle-t.GetAlpha();
	t.SetRelativeSector(relativesector);
	if (!t.Rotate(rotation)) {
          t.SetClusterIndex(nr, t.GetClusterIndex(nr) | 0x8000);
          return 0;
        }	
      }
      if (!t.PropagateTo(cl->GetX())) { // RS: go directly to cluster X
        t.SetClusterIndex(nr, t.GetClusterIndex(nr) | 0x8000);
        return 0;
      }
      //
      t.SetCurrentCluster(cl); 
      t.SetRow(nr); // RS: memorise row
      Int_t accept = AcceptCluster(&t,t.GetCurrentCluster());
      if ((tpcindex&0x8000)==0) accept =0;
      if (accept<3) { 
	//if founded cluster is acceptible
	if (cl->IsUsed(11)) {  // id cluster is shared inrease uncertainty
	  t.SetErrorY2(t.GetErrorY2()+0.03);
	  t.SetErrorZ2(t.GetErrorZ2()+0.03); 
	  t.SetErrorY2(t.GetErrorY2()*3);
	  t.SetErrorZ2(t.GetErrorZ2()*3); 
	}
	t.SetNFoundable(t.GetNFoundable()+1);
	UpdateTrack(&t,accept);
	return 1;
      }
      else { // Remove old cluster from track
	t.SetClusterIndex(nr, -3);
	//RS t.SetClusterPointer(nr, 0);
      }
    }
  }
  if (TMath::Abs(t.GetSnp())>AliTPCReconstructor::GetMaxSnpTracker()) return 0;  // cut on angle
  if (fIteration>1 && IsFindable(t)){
    // not look for new cluster during refitting    
    t.SetNFoundable(t.GetNFoundable()+1);
    return 0;
  }
  //
  UInt_t index=0;
  //  if (TMath::Abs(t.GetSnp())>0.95 || TMath::Abs(x*t.GetC()-t.GetEta())>0.95) return 0;// patch 28 fev 06
  //
  if (fAccountDistortions && !DistortX(&t,x,nr)) {if (fIteration==0) t.SetRemoval(10); return 0;}
  if (!t.PropagateTo(x))                         {if (fIteration==0) t.SetRemoval(10); return 0;}
  t.SetRow(nr); //RS:? memorise reached row?
  //
  double y=t.GetY(), z=t.GetZ();
  double yEdgeDist =  y;
  if  (fAccountDistortions) yEdgeDist -= GetYSectEdgeDist(t.GetRelativeSector(),nr,y,z);
  Bool_t rot = kFALSE;
  if (y>0 && yEdgeDist>ymax) { //RS y sign here is used to deduce which edge is used
    t.SetRelativeSector((t.GetRelativeSector() + 1) % fN);
    if (!t.Rotate(fSectors->GetAlpha())) return 0;
    rot = kTRUE;
  }
  else if (y<0 && yEdgeDist<-ymax) {
    t.SetRelativeSector((t.GetRelativeSector() + fN - 1) % fN);
    if (!t.Rotate(-fSectors->GetAlpha())) return 0;
    rot = kTRUE;
  }
  if (rot) {
    x = GetXrow(nr);
    if (fAccountDistortions && !DistortX(&t,x,nr)) {if (fIteration==0) t.SetRemoval(10); return 0;}
    if (!t.PropagateTo(x))                         {if (fIteration==0) t.SetRemoval(10); return 0; }
    y = t.GetY();
    yEdgeDist =  y;
    if (fAccountDistortions) yEdgeDist -= GetYSectEdgeDist(t.GetRelativeSector(),nr,y,z);
  }
  //
  if (!IsActive(t.GetRelativeSector(),nr)) { // RS:? How distortions affect this
    t.SetInDead(kTRUE);
    t.SetClusterIndex2(nr,-1); 
    return 0;
  }
  //AliInfo(Form("A - Sector%d phi %f - alpha %f", t.fRelativeSector,y/x, t.GetAlpha()));
  Bool_t isActive  = IsActive(t.GetRelativeSector(),nr); //RS:? Why do we check this again?
  Bool_t isActive2 = (nr>=fInnerSec->GetNRows()) ? 
    fOuterSec[t.GetRelativeSector()][nr-fInnerSec->GetNRows()].GetN()>0 : 
    fInnerSec[t.GetRelativeSector()][nr].GetN()>0;
  
  if (!isActive || !isActive2) return 0;

  const AliTPCtrackerRow &krow=GetRow(t.GetRelativeSector(),nr);
  if ( (t.GetSigmaY2()<0) || t.GetSigmaZ2()<0) return 0;
  //
  // RS: account for eventual modifications in dead zone definition due to distortions
  double margin = (y>0 ? ymax-yEdgeDist : ymax + yEdgeDist);
  if (margin<krow.GetDeadZone()){
    t.SetInDead(kTRUE);
    t.SetClusterIndex2(nr,-1); 
    return 0;
  } 
  else {
    if (IsFindable(t))  t.SetNFoundable(t.GetNFoundable()+1);
    else                return 0;
  }   
  //calculate 
  if (krow && (cl=krow.FindNearest2(y,z,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ,index)) ) t.SetCurrentClusterIndex1(krow.GetIndex(index));
  if (cl) {
    t.SetCurrentCluster(cl); 
    //    t.SetRow(nr); //RS: memorise row | already set at propagation
    if (fIteration==2&&cl->IsUsed(10)) return 0; 
    Int_t accept = AcceptCluster(&t,t.GetCurrentCluster());
    if (fIteration==2&&cl->IsUsed(11)) {
      t.SetErrorY2(t.GetErrorY2()+0.03);
      t.SetErrorZ2(t.GetErrorZ2()+0.03); 
      t.SetErrorY2(t.GetErrorY2()*3);
      t.SetErrorZ2(t.GetErrorZ2()*3); 
    }
    /*    
    if (t.fCurrentCluster->IsUsed(10)){
      //
      //     

      t.fNShared++;
      if (t.fNShared>0.7*t.GetNumberOfClusters()) {
	t.fRemoval =10;
	return 0;
      }
    }
    */
    if (accept<3) UpdateTrack(&t,accept);

  } else {  
    if ( fIteration==0 && t.GetNFoundable()*0.5 > t.GetNumberOfClusters()) t.SetRemoval(10);
    
  }
  return 1;
}



//_________________________________________________________________________
Bool_t AliTPCtracker::GetTrackPoint(Int_t index, AliTrackPoint &p ) const
{
  // Get track space point by index
  // return false in case the cluster doesn't exist
  AliTPCclusterMI *cl = GetClusterMI(index);
  if (!cl) return kFALSE;
  Int_t sector = (index&0xff000000)>>24;
  //  Int_t row = (index&0x00ff0000)>>16;
  Float_t xyz[3];
  //  xyz[0] = fkParam->GetPadRowRadii(sector,row);
  xyz[0] = cl->GetX();
  xyz[1] = cl->GetY();
  xyz[2] = cl->GetZ();
  Float_t sin,cos;
  fkParam->AdjustCosSin(sector,cos,sin);
  Float_t x = cos*xyz[0]-sin*xyz[1];
  Float_t y = cos*xyz[1]+sin*xyz[0];
  Float_t cov[6];
  Float_t sigmaY2 = 0.027*cl->GetSigmaY2();
  if (sector < fkParam->GetNInnerSector()) sigmaY2 *= 2.07;
  Float_t sigmaZ2 = 0.066*cl->GetSigmaZ2();
  if (sector < fkParam->GetNInnerSector()) sigmaZ2 *= 1.77;
  cov[0] = sin*sin*sigmaY2;
  cov[1] = -sin*cos*sigmaY2;
  cov[2] = 0.;
  cov[3] = cos*cos*sigmaY2;
  cov[4] = 0.;
  cov[5] = sigmaZ2;
  p.SetXYZ(x,y,xyz[2],cov);
  AliGeomManager::ELayerID iLayer;
  Int_t idet;
  if (sector < fkParam->GetNInnerSector()) {
    iLayer = AliGeomManager::kTPC1;
    idet = sector;
  }
  else {
    iLayer = AliGeomManager::kTPC2;
    idet = sector - fkParam->GetNInnerSector();
  }
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,idet);
  p.SetVolumeID(volid);
  return kTRUE;
}



Int_t AliTPCtracker::UpdateClusters(AliTPCseed& t,  Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------
  t.SetCurrentCluster(0);
  t.SetCurrentClusterIndex1(-3);   
   
  Int_t     row = (fAccountDistortions ? GetRowNumber(&t):GetRowNumber(t.GetX()))-1;  //RS: with distortions cannot rely on X
  Double_t  ymax= GetMaxY(nr);

  if (row < nr) return 1; // don't prolongate if not information until now -

  const double kRoadY = 1., kRoadZ = 1.;

  Double_t x= GetXrow(nr);
  if (fAccountDistortions && !DistortX(&t,x,nr)) return 0; // RS: if needed, account distortion
  //
  if (!t.PropagateTo(x)) return 0;
  t.SetRow(nr); //RS: memorise row
  //
  double y=t.GetY(), z=t.GetZ();
  double yEdgeDist =  y;
  if  (fAccountDistortions) yEdgeDist -= GetYSectEdgeDist(t.GetRelativeSector(),nr,y,z);
  if (y>0 && yEdgeDist>ymax) { //RS y sign is used to determine which edge is used
    t.SetRelativeSector((t.GetRelativeSector()+1) % fN);
    if (!t.Rotate(fSectors->GetAlpha())) return 0;
    t.SetRow(-1); //RS: after rotation the row is not known
    return 1;
  }
  else if (y<0 && yEdgeDist<-ymax) {
    t.SetRelativeSector((t.GetRelativeSector()+fN-1) % fN);
    if (!t.Rotate(-fSectors->GetAlpha())) return 0;
    t.SetRow(-1); //RS: after rotation the row is not known
    return 1;
  }
  //
  if (TMath::Abs(t.GetSnp())>AliTPCReconstructor::GetMaxSnpTracker()) return 0;

  if (!IsActive(t.GetRelativeSector(),nr)) {
    t.SetInDead(kTRUE);
    t.SetClusterIndex2(nr,-1); 
    return 0;
  }
  //AliInfo(Form("A - Sector%d phi %f - alpha %f", t.fRelativeSector,y/x, t.GetAlpha()));

  AliTPCtrackerRow &krow=GetRow(t.GetRelativeSector(),nr);
  double margin = (y>0 ? ymax-yEdgeDist : ymax + yEdgeDist);
  if (margin<krow.GetDeadZone()) {
    t.SetInDead(kTRUE);
    t.SetClusterIndex2(nr,-1); 
    return 0;
  } 
  else {
    //      if (TMath::Abs(t.GetZ())<(AliTPCReconstructor::GetCtgRange()*t.GetX()+10) && (TMath::Abs(t.GetSnp())<AliTPCReconstructor::GetMaxSnpTracker())) 
    if (IsFindable(t)) t.SetNFoundable(t.GetNFoundable()+1);
    else return 0;      
  }

  // update current
  if ( (nr%2==0) || t.GetNumberOfClusters()<2){
    //    t.fCurrentSigmaY = GetSigmaY(&t);
    //t.fCurrentSigmaZ = GetSigmaZ(&t);
    GetShape(&t,nr);
  }
    
  AliTPCclusterMI *cl=0;
  Int_t index=0;
  //
  Double_t roady = 1.;
  Double_t roadz = 1.;
  //

  if (!cl){
    index = t.GetClusterIndex2(nr);    
    if ( (index >= 0) && (index&0x8000)==0){
      //RS cl = t.GetClusterPointer(nr);
      //RS if ( (cl==0) && (index >= 0)) cl = GetClusterMI(index);
      if ( index >= 0 ) cl = GetClusterMI(index);
      t.SetCurrentClusterIndex1(index);
      if (cl) {
	t.SetCurrentCluster(cl);
	return 1;
      }
    }
  }

  //  if (index<0) return 0;
  UInt_t uindex = TMath::Abs(index);

  if (krow) {
    //cl = krow.FindNearest2(y+10,z,roady,roadz,uindex);      
    cl = krow.FindNearest2(y,z,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ,uindex);      
  }

  if (cl) t.SetCurrentClusterIndex1(krow.GetIndex(uindex));   
  t.SetCurrentCluster(cl);

  return 1;
}


Int_t AliTPCtracker::FollowToNextCluster(AliTPCseed & t, Int_t nr) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation to next pad row
  //-----------------------------------------------------------------

  //update error according neighborhoud

  if (t.GetCurrentCluster()) {
    t.SetRow(nr); 
    Int_t accept = AcceptCluster(&t,t.GetCurrentCluster());
    
    if (t.GetCurrentCluster()->IsUsed(10)){
      //
      //
      //  t.fErrorZ2*=2;
      //  t.fErrorY2*=2;
      t.SetNShared(t.GetNShared()+1);
      if (t.GetNShared()>0.7*t.GetNumberOfClusters()) {
	t.SetRemoval(10);
	return 0;
      }
    }   
    if (fIteration>0) accept = 0;
    if (accept<3)  UpdateTrack(&t,accept);  
 
  } else {
    if (fIteration==0){
      //RS: with distortions related cluster errors the track error may grow, don't use this cut
      //if ( t.GetNumberOfClusters()>18 && ( (t.GetSigmaY2()+t.GetSigmaZ2())>0.16 + fExtraClErrYZ2 )) t.SetRemoval(10); 
      if ( t.GetNumberOfClusters()>18 && t.GetChi2()/t.GetNumberOfClusters()>6 ) t.SetRemoval(10);      

      if (( (t.GetNFoundable()*0.5 > t.GetNumberOfClusters()) || t.GetNoCluster()>15)) t.SetRemoval(10);
    }
  }
  return 1;
}



//_____________________________________________________________________________
Int_t AliTPCtracker::FollowProlongation(AliTPCseed& t, Int_t rf, Int_t step, Bool_t fromSeeds) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  // Double_t xt=t.GetX(); // RS: with distortions we cannot rely on rowID from X
  //
  Double_t alpha=t.GetAlpha();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  //
  t.SetRelativeSector(Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN);
    
  Int_t first = fAccountDistortions ? GetRowNumber(&t):GetRowNumber(t.GetX());  //RS: with distortions cannot rely on X
  //Int_t first = GetRowNumber(xt);
  if (!fromSeeds)
    first -= step;
  if (first < 0)
    first = 0;
  for (Int_t nr= first; nr>=rf; nr-=step) {
    // update kink info
    if (t.GetKinkIndexes()[0]>0){
      for (Int_t i=0;i<3;i++){
	Int_t index = t.GetKinkIndexes()[i];
	if (index==0) break;
	if (index<0) continue;
	//
	AliKink * kink = (AliKink*)fEvent->GetKink(index-1);
	if (!kink){
	  printf("PROBLEM\n");
	}
	else{
	  Int_t kinkrow = kink->GetTPCRow0()+2+Int_t(0.5/(0.05+kink->GetAngle(2)));
	  if (kinkrow==nr){
	    AliExternalTrackParam paramd(t);
	    kink->SetDaughter(paramd);
	    kink->SetStatus(2,5);
	    kink->Update();
	  }
	}
      }
    }

    if (nr==80) t.UpdateReference();
    if (nr<fInnerSec->GetNRows()) 
      fSectors = fInnerSec;
    else
      fSectors = fOuterSec;
    if (FollowToNext(t,nr)==0) 
      if (!t.IsActive()) 
	return 0;
    
  }   
  return 1;
}






Int_t AliTPCtracker::FollowBackProlongation(AliTPCseed& t, Int_t rf, Bool_t fromSeeds) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  //
  //  Double_t xt=t.GetX();   //RS: with distortions cannot rely on row from X
  Double_t alpha=t.GetAlpha();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  t.SetRelativeSector(Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN);
    
  Int_t first = t.GetFirstPoint();
  Int_t ri = fAccountDistortions ? GetRowNumber(&t) :  GetRowNumber(t.GetX());  //RS: with distortions cannot rely on X
  //  Int_t ri = GetRowNumber(xt); 
  if (!fromSeeds)
    ri += 1;

  if (first<ri) first = ri;
  //
  if (first<0) first=0;
  for (Int_t nr=first; nr<=rf; nr++) {
    //    if ( (TMath::Abs(t.GetSnp())>0.95)) break;//patch 28 fev 06
    if (t.GetKinkIndexes()[0]<0){
      for (Int_t i=0;i<3;i++){
	Int_t index = t.GetKinkIndexes()[i];
	if (index==0) break;
	if (index>0) continue;
	index = TMath::Abs(index);
	AliKink * kink = (AliKink*)fEvent->GetKink(index-1);
	if (!kink){
	  printf("PROBLEM\n");
	}
	else{
	  Int_t kinkrow = kink->GetTPCRow0()-2-Int_t(0.5/(0.05+kink->GetAngle(2)));
	  if (kinkrow==nr){
	    AliExternalTrackParam paramm(t);
	    kink->SetMother(paramm);
	    kink->SetStatus(2,1);
	    kink->Update();
	  }
	}
      }      
    }
    //
    if (nr<fInnerSec->GetNRows()) 
      fSectors = fInnerSec;
    else
      fSectors = fOuterSec;

    FollowToNext(t,nr);                                                             
  }   
  return 1;
}




   
Float_t AliTPCtracker::OverlapFactor(AliTPCseed * s1, AliTPCseed * s2, Int_t &sum1, Int_t & sum2)
{
  // overlapping factor
  //
  sum1=0;
  sum2=0;
  Int_t sum=0;
  //
  Float_t dz2 =(s1->GetZ() - s2->GetZ());
  dz2*=dz2;  

  Float_t dy2 =TMath::Abs((s1->GetY() - s2->GetY()));
  dy2*=dy2;
  Float_t distance = TMath::Sqrt(dz2+dy2);
  if (distance>4.) return 0; // if there are far away  - not overlap - to reduce combinatorics
 
  //  Int_t offset =0;
  Int_t firstpoint = TMath::Min(s1->GetFirstPoint(),s2->GetFirstPoint());
  Int_t lastpoint = TMath::Max(s1->GetLastPoint(),s2->GetLastPoint());
  if (lastpoint>=kMaxRow) 
    lastpoint = kMaxRow-1;
  if (firstpoint<0) 
    firstpoint = 0;
  if (firstpoint>lastpoint) {
    firstpoint =lastpoint;
    //    lastpoint  = kMaxRow-1;
  }
    
  
  // for (Int_t i=firstpoint-1;i<lastpoint+1;i++){
  //   if (s1->GetClusterIndex2(i)>0) sum1++;
  //   if (s2->GetClusterIndex2(i)>0) sum2++;
  //   if (s1->GetClusterIndex2(i)==s2->GetClusterIndex2(i) && s1->GetClusterIndex2(i)>0) {
  //     sum++;
  //   }
  // }
  // RS: faster version + fix(?): clusterindex starts from 0 
  for (Int_t i=firstpoint-1;i<lastpoint+1;i++){
    int ind1=s1->GetClusterIndex2(i), ind2=s2->GetClusterIndex2(i);
    if (ind1>=0) sum1++;
    if (ind2>=0) sum2++;
    if (ind1==ind2 && ind1>=0) sum++;
  }
  if (sum<5) return 0;
  
  Float_t summin = TMath::Min(sum1+1,sum2+1);
  Float_t ratio = (sum+1)/Float_t(summin);
  return ratio;
}

void  AliTPCtracker::SignShared(AliTPCseed * s1, AliTPCseed * s2)
{
  // shared clusters
  //
  Float_t thetaCut = 0.2;//+10.*TMath::Sqrt(s1->GetSigmaTglZ()+ s2->GetSigmaTglZ());
  if (TMath::Abs(s1->GetTgl()-s2->GetTgl())>thetaCut) return;
  Float_t minCl = TMath::Min(s1->GetNumberOfClusters(),s2->GetNumberOfClusters());
  Int_t cutN0 = TMath::Max(5,TMath::Nint(0.1*minCl));
  
  //
  Int_t sumshared=0;
  //
  //Int_t firstpoint = TMath::Max(s1->GetFirstPoint(),s2->GetFirstPoint());
  //Int_t lastpoint = TMath::Min(s1->GetLastPoint(),s2->GetLastPoint());
  Int_t firstpoint = 0;
  Int_t lastpoint = kMaxRow;
  //
  //  if (firstpoint>=lastpoint-5) return;;

  for (Int_t i=firstpoint;i<lastpoint;i++){
    //    if ( (s1->GetClusterIndex2(i)&0xFFFF8FFF)==(s2->GetClusterIndex2(i)&0xFFFF8FFF) && s1->GetClusterIndex2(i)>0) {
    if ( (s1->GetClusterIndex2(i))==(s2->GetClusterIndex2(i)) && s1->GetClusterIndex2(i)>=0) {
      sumshared++;
    }
  }
  if (sumshared>cutN0){
    // sign clusters
    //
    for (Int_t i=firstpoint;i<lastpoint;i++){
      //      if ( (s1->GetClusterIndex2(i)&0xFFFF8FFF)==(s2->GetClusterIndex2(i)&0xFFFF8FFF) && s1->GetClusterIndex2(i)>0) {
      if ( (s1->GetClusterIndex2(i))==(s2->GetClusterIndex2(i)) && s1->GetClusterIndex2(i)>=0) {
	const AliTPCTrackerPoints::Point *p1  = s1->GetTrackPoint(i);
	const AliTPCTrackerPoints::Point *p2  = s2->GetTrackPoint(i);; 
	if (s1->IsActive()&&s2->IsActive()){
	  s1->SetShared(i);
	  s2->SetShared(i);
	}	
      }
    }
  }
  //  
  if (sumshared>cutN0){
    for (Int_t i=0;i<4;i++){
      if (s1->GetOverlapLabel(3*i)==0){
	s1->SetOverlapLabel(3*i,  s2->GetLabel());
	s1->SetOverlapLabel(3*i+1,sumshared);
	s1->SetOverlapLabel(3*i+2,s2->GetUniqueID());
	break;
      }	
    }
    for (Int_t i=0;i<4;i++){
      if (s2->GetOverlapLabel(3*i)==0){
	s2->SetOverlapLabel(3*i,  s1->GetLabel());
	s2->SetOverlapLabel(3*i+1,sumshared);
	s2->SetOverlapLabel(3*i+2,s1->GetUniqueID());
	break;
      }	
    }    
  }
}

void  AliTPCtracker::SignShared(TObjArray * arr)
{
  //
  //sort trackss according sectors
  //  
  for (Int_t i=0; i<arr->GetEntriesFast(); i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    //if (pt) RotateToLocal(pt);
    pt->SetSort(0);
  }
  arr->UnSort();
  arr->Sort();  // sorting according relative sectors
  arr->Expand(arr->GetEntries());
  //
  //
  Int_t nseed=arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    for (Int_t j=0;j<12;j++){
      pt->SetOverlapLabel(j,0);
    }
  }
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) continue;
    if (pt->GetRemoval()>10) continue;
    for (Int_t j=i+1; j<nseed; j++){
      AliTPCseed *pt2=(AliTPCseed*)arr->UncheckedAt(j);
      if (TMath::Abs(pt->GetRelativeSector()-pt2->GetRelativeSector())>1) continue;
      //      if (pt2){
      if (pt2->GetRemoval()<=10) {
	//if ( TMath::Abs(pt->GetRelativeSector()-pt2->GetRelativeSector())>0) break;
	SignShared(pt,pt2);
      }
    }  
  }
}


void AliTPCtracker::SortTracks(TObjArray * arr, Int_t mode) const
{
  //
  //sort tracks in array according mode criteria
  Int_t nseed = arr->GetEntriesFast();    
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    pt->SetSort(mode);
  }
  arr->UnSort();
  arr->Sort();
}


void AliTPCtracker::RemoveUsed2(TObjArray * arr, Float_t factor1,  Float_t factor2, Int_t minimal)
{
  //
  // Loop over all tracks and remove overlaped tracks (with lower quality)
  // Algorithm:
  //    1. Unsign clusters
  //    2. Sort tracks according quality
  //       Quality is defined by the number of cluster between first and last points 
  //       
  //    3. Loop over tracks - decreasing quality order
  //       a.) remove - If the fraction of shared cluster less than factor (1- n or 2) 
  //       b.) remove - If the minimal number of clusters less than minimal and not ITS
  //       c.) if track accepted - sign clusters
  //
  //Called in - AliTPCtracker::Clusters2Tracks()
  //          - AliTPCtracker::PropagateBack()
  //          - AliTPCtracker::RefitInward()
  //
  // Arguments:
  //   factor1 - factor for constrained
  //   factor2 -        for non constrained tracks 
  //            if (Float_t(shared+1)/Float_t(found+1)>factor) - DELETE
  //
  UnsignClusters();
  //
  Int_t nseed = arr->GetEntriesFast();  
  //  Float_t quality = new Float_t[nseed];
  //  Int_t   * indexes = new Int_t[nseed];
  Float_t quality[nseed]; //RS Use stack allocations
  Int_t   indexes[nseed];
  Int_t good =0;
  //
  //
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt){
      quality[i]=-1;
      continue;
    }
    pt->UpdatePoints();    //select first last max dens points
    Float_t * points = pt->GetPoints();
    if (points[3]<0.8) quality[i] =-1;
    quality[i] = (points[2]-points[0])+pt->GetNumberOfClusters();
    //prefer high momenta tracks if overlaps
    quality[i] *= TMath::Sqrt(TMath::Abs(pt->Pt())+0.5); 
  }
  TMath::Sort(nseed,quality,indexes);
  //
  //
  for (Int_t itrack=0; itrack<nseed; itrack++) {
    Int_t trackindex = indexes[itrack];
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(trackindex);   
    if (!pt) continue;
    //
    if (quality[trackindex]<0){
      MarkSeedFree( arr->RemoveAt(trackindex) );
      continue;
    }
    //
    //
    Int_t first = Int_t(pt->GetPoints()[0]);
    Int_t last  = Int_t(pt->GetPoints()[2]);
    Double_t factor = (pt->GetBConstrain()) ? factor1: factor2;
    //
    Int_t found,foundable,shared;
    GetSeedClusterStatistic(pt, first,last, found, foundable,shared,kFALSE); //RS: seeds don't keep their clusters
    //RS pt->GetClusterStatistic(first,last, found, foundable,shared,kFALSE); // better to get statistic in "high-dens"  region do't use full track as in line bellow
    //MI  pt->GetClusterStatistic(0,kMaxRow, found, foundable,shared,kFALSE);
    Bool_t itsgold =kFALSE;
    if (pt->GetESD()){
      Int_t dummy[12];
      if (pt->GetESD()->GetITSclusters(dummy)>4) itsgold= kTRUE;
    }
    if (!itsgold){
      //
      if (Float_t(shared+1)/Float_t(found+1)>factor){
	if (pt->GetKinkIndexes()[0]!=0) continue;  //don't remove tracks  - part of the kinks
	if( (AliTPCReconstructor::StreamLevel()&kStreamRemoveUsed)>0){ // flag:stream  information about TPC tracks which were descarded (double track removal)
	  TTreeSRedirector &cstream = *fDebugStreamer;
	  cstream<<"RemoveUsed"<<
	    "iter="<<fIteration<<
	    "pt.="<<pt<<
	    "\n";
	}
	MarkSeedFree( arr->RemoveAt(trackindex) );
	continue;
      }      
      if (pt->GetNumberOfClusters()<50&&(found-0.5*shared)<minimal){  //remove short tracks
	if (pt->GetKinkIndexes()[0]!=0) continue;  //don't remove tracks  - part of the kinks
	if( (AliTPCReconstructor::StreamLevel()&kStreamRemoveShort)>0){ // flag:stream  information about TPC tracks which were discarded (short track removal)
	  TTreeSRedirector &cstream = *fDebugStreamer;
	  cstream<<"RemoveShort"<<
	    "iter="<<fIteration<<
	    "pt.="<<pt<<
	    "\n";
	}
	MarkSeedFree( arr->RemoveAt(trackindex) );
	continue;
      }
    }

    good++;
    //if (sharedfactor>0.4) continue;
    if (pt->GetKinkIndexes()[0]>0) continue;
    //Remove tracks with undefined properties - seems
    if (pt->GetSigmaY2()<kAlmost0) continue; // ? what is the origin ? 
    //
    for (Int_t i=first; i<last; i++) {
      Int_t index=pt->GetClusterIndex2(i);
      // if (index<0 || index&0x8000 ) continue;
      if (index<0 || index&0x8000 ) continue;
      AliTPCclusterMI *c= GetClusterMI(index); //RS pt->GetClusterPointer(i);        
      if (!c) continue;
      c->Use(10);  
    }    
  }
  fNtracks = good;
  if (fDebug>0){
    Info("RemoveUsed","\n*****\nNumber of good tracks after shared removal\t%d\n",fNtracks);
  }
  //  delete []quality; // RS was moved to stack allocation
  //  delete []indexes;
}

void AliTPCtracker::DumpClusters(Int_t iter, TObjArray *trackArray) 
{
  //
  // Dump clusters after reco
  // signed and unsigned cluster can be visualized   
  // 1. Unsign all cluster
  // 2. Sign all used clusters
  // 3. Dump clusters
  UnsignClusters();
  Int_t nseed = trackArray->GetEntries();
  for (Int_t i=0; i<nseed; i++){
    AliTPCseed *pt=(AliTPCseed*)trackArray->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    Bool_t isKink=pt->GetKinkIndex(0)!=0;
    for (Int_t j=0; j<kMaxRow; ++j) {
      Int_t index=pt->GetClusterIndex2(j);
      if (index<0) continue;
      AliTPCclusterMI *c= GetClusterMI(index); //RS pt->GetClusterPointer(j);
      if (!c) continue;
      if (isKink) c->Use(100);   // kink
      c->Use(10);                // by default usage 10
    }
  }
  //
  Int_t eventNr = fEvent->GetEventNumberInFile();

  for (Int_t sec=0;sec<fkNIS;sec++){
    for (Int_t row=0;row<fInnerSec->GetNRows();row++){
      TClonesArray *cla = fInnerSec[sec][row].GetClusters1();
      for (Int_t icl =0;icl< fInnerSec[sec][row].GetN1();icl++){    
	AliTPCclusterMI* cl = (AliTPCclusterMI*)cla->At(icl);
	Float_t gx[3];	cl->GetGlobalXYZ(gx);
	(*fDebugStreamer)<<"clDump"<< 
	  "eventNr="<<eventNr<<
	  "iter="<<iter<<
	  "cl.="<<cl<<      
	  "gx0="<<gx[0]<<
	  "gx1="<<gx[1]<<
	  "gx2="<<gx[2]<<
	  "\n";
      }
      cla = fInnerSec[sec][row].GetClusters2();
      for (Int_t icl =0;icl< fInnerSec[sec][row].GetN2();icl++){
	AliTPCclusterMI* cl = (AliTPCclusterMI*)cla->At(icl);
	Float_t gx[3];	cl->GetGlobalXYZ(gx);
	(*fDebugStreamer)<<"clDump"<< 
	  "eventNr="<<eventNr<<
	  "iter="<<iter<<
	  "cl.="<<cl<<
	  "gx0="<<gx[0]<<
	  "gx1="<<gx[1]<<
	  "gx2="<<gx[2]<<
	  "\n";
      }
    }
  }
  
  for (Int_t sec=0;sec<fkNOS;sec++){
    for (Int_t row=0;row<fOuterSec->GetNRows();row++){
      TClonesArray *cla = fOuterSec[sec][row].GetClusters1();
      for (Int_t icl =0;icl< fOuterSec[sec][row].GetN1();icl++){
	Float_t gx[3];	
	AliTPCclusterMI* cl = (AliTPCclusterMI*) cla->At(icl);
	cl->GetGlobalXYZ(gx);
	(*fDebugStreamer)<<"clDump"<< 
	  "eventNr="<<eventNr<<
	  "iter="<<iter<<
	  "cl.="<<cl<<
	  "gx0="<<gx[0]<<
	  "gx1="<<gx[1]<<
	  "gx2="<<gx[2]<<
	  "\n";      
      }
      cla = fOuterSec[sec][row].GetClusters2();
      for (Int_t icl =0;icl< fOuterSec[sec][row].GetN2();icl++){
	Float_t gx[3];	
	AliTPCclusterMI* cl = (AliTPCclusterMI*) cla->At(icl);
	cl->GetGlobalXYZ(gx);
	(*fDebugStreamer)<<"clDump"<< 
	  "eventNr="<<eventNr<<
	  "iter="<<iter<<
	  "cl.="<<cl<<
	  "gx0="<<gx[0]<<
	  "gx1="<<gx[1]<<
	  "gx2="<<gx[2]<<
	  "\n";      
      }
    }
  }
  
}
void AliTPCtracker::UnsignClusters() 
{
  //
  // loop over all clusters and unsign them
  //
  
  for (Int_t sec=0;sec<fkNIS;sec++){
    for (Int_t row=0;row<fInnerSec->GetNRows();row++){
      TClonesArray *cla = fInnerSec[sec][row].GetClusters1();
      for (Int_t icl =0;icl< fInnerSec[sec][row].GetN1();icl++)
	//	if (cl[icl].IsUsed(10)) 	
	((AliTPCclusterMI*) cla->At(icl))->Use(-1);
      cla = fInnerSec[sec][row].GetClusters2();
      for (Int_t icl =0;icl< fInnerSec[sec][row].GetN2();icl++)
	//if (cl[icl].IsUsed(10)) 	
	((AliTPCclusterMI*) cla->At(icl))->Use(-1);
    }
  }
  
  for (Int_t sec=0;sec<fkNOS;sec++){
    for (Int_t row=0;row<fOuterSec->GetNRows();row++){
      TClonesArray *cla = fOuterSec[sec][row].GetClusters1();
      for (Int_t icl =0;icl< fOuterSec[sec][row].GetN1();icl++)
	//if (cl[icl].IsUsed(10)) 	
	((AliTPCclusterMI*) cla->At(icl))->Use(-1);
      cla = fOuterSec[sec][row].GetClusters2();
      for (Int_t icl =0;icl< fOuterSec[sec][row].GetN2();icl++)
	//if (cl[icl].IsUsed(10)) 	
	((AliTPCclusterMI*) cla->At(icl))->Use(-1);
    }
  }
  
}



void AliTPCtracker::SignClusters(const TObjArray * arr, Float_t fnumber, Float_t fdensity)
{
  //
  //sign clusters to be "used"
  //
  // snumber and sdensity sign number of sigmas - bellow mean value to be accepted
  // loop over "primaries"
  
  Float_t sumdens=0;
  Float_t sumdens2=0;
  Float_t sumn   =0;
  Float_t sumn2  =0;
  Float_t sumchi =0;
  Float_t sumchi2 =0;

  Float_t sum    =0;

  TStopwatch timer;
  timer.Start();

  Int_t nseed = arr->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    if (!(pt->IsActive())) continue;
    Float_t dens = pt->GetNumberOfClusters()/Float_t(pt->GetNFoundable());
    if ( (dens>0.7) && (pt->GetNumberOfClusters()>70)){
      sumdens += dens;
      sumdens2+= dens*dens;
      sumn    += pt->GetNumberOfClusters();
      sumn2   += pt->GetNumberOfClusters()*pt->GetNumberOfClusters();
      Float_t chi2 = pt->GetChi2()/pt->GetNumberOfClusters();
      if (chi2>5) chi2=5;
      sumchi  +=chi2;
      sumchi2 +=chi2*chi2;
      sum++;
    }
  }

  Float_t mdensity = 0.9;
  Float_t meann    = 130;
  Float_t meanchi  = 1;
  Float_t sdensity = 0.1;
  Float_t smeann    = 10;
  Float_t smeanchi  =0.4;
  

  if (sum>20){
    mdensity = sumdens/sum;
    meann    = sumn/sum;
    meanchi  = sumchi/sum;
    //
    sdensity = sumdens2/sum-mdensity*mdensity;
    if (sdensity >= 0)
       sdensity = TMath::Sqrt(sdensity);
    else
       sdensity = 0.1;
    //
    smeann   = sumn2/sum-meann*meann;
    if (smeann >= 0)
      smeann   = TMath::Sqrt(smeann);
    else 
      smeann   = 10;
    //
    smeanchi = sumchi2/sum - meanchi*meanchi;
    if (smeanchi >= 0)
      smeanchi = TMath::Sqrt(smeanchi);
    else
      smeanchi = 0.4;
  }


  //REMOVE  SHORT DELTAS or tracks going out of sensitive volume of TPC
  //
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }
    if (pt->GetBSigned()) continue;
    if (pt->GetBConstrain()) continue;    
    //if (!(pt->IsActive())) continue;
    /*
    Int_t found,foundable,shared;    
    pt->GetClusterStatistic(0,kMaxRow,found, foundable,shared);
    if (shared/float(found)>0.3) {
      if (shared/float(found)>0.9 ){
	//MarkSeedFree( arr->RemoveAt(i) );
      }
      continue;
    }
    */
    Bool_t isok =kFALSE;
    if ( (pt->GetNShared()/pt->GetNumberOfClusters()<0.5) &&pt->GetNumberOfClusters()>60)
      isok = kTRUE;
    if ((TMath::Abs(1/pt->GetC())<100.) && (pt->GetNShared()/pt->GetNumberOfClusters()<0.7))
      isok =kTRUE;
    if  (TMath::Abs(pt->GetZ()/pt->GetX())>1.1)
      isok =kTRUE;
    if ( (TMath::Abs(pt->GetSnp()>0.7) && pt->GetD(0,0)>60.))
      isok =kTRUE;
    
    if (isok)     
      for (Int_t j=0; j<kMaxRow; ++j) {	
	Int_t index=pt->GetClusterIndex2(j);
	if (index<0) continue;
	AliTPCclusterMI *c= GetClusterMI(index);//pt->GetClusterPointer(j);
	if (!c) continue;
	//if (!(c->IsUsed(10))) c->Use();  
	c->Use(10);  
      }
  }
  
  
  //
  Double_t maxchi  = meanchi+2.*smeanchi;

  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i);    
    if (!pt) {
      continue;
    }    
    //if (!(pt->IsActive())) continue;
    if (pt->GetBSigned()) continue;
    Double_t chi     = pt->GetChi2()/pt->GetNumberOfClusters();
    if (chi>maxchi) continue;

    Float_t bfactor=1;
    Float_t dens = pt->GetNumberOfClusters()/Float_t(pt->GetNFoundable());
   
    //sign only tracks with enoug big density at the beginning
    
    if ((pt->GetDensityFirst(40)<0.75) && pt->GetNumberOfClusters()<meann) continue; 
    
    
    Double_t mindens = TMath::Max(double(mdensity-sdensity*fdensity*bfactor),0.65);
    Double_t minn    = TMath::Max(Int_t(meann-fnumber*smeann*bfactor),50);
   
    //    if (pt->fBConstrain) mindens = TMath::Max(mdensity-sdensity*fdensity*bfactor,0.65);
    if ( (pt->GetRemoval()==10) && (pt->GetSnp()>0.8)&&(dens>mindens))
      minn=0;

    if ((dens>mindens && pt->GetNumberOfClusters()>minn) && chi<maxchi ){
      //Int_t noc=pt->GetNumberOfClusters();
      pt->SetBSigned(kTRUE);
      for (Int_t j=0; j<kMaxRow; ++j) {

	Int_t index=pt->GetClusterIndex2(j);
	if (index<0) continue;
	AliTPCclusterMI *c= GetClusterMI(index); //RS pt->GetClusterPointer(j);
	if (!c) continue;
	//	if (!(c->IsUsed(10))) c->Use();  
	c->Use(10);  
      }
    }
  }
  //  gLastCheck = nseed;
  //  arr->Compress();
  if (fDebug>0){
    timer.Print();
  }
}



Int_t AliTPCtracker::RefitInward(AliESDEvent *event)
{
  //
  // back propagation of ESD tracks
  //
  //return 0;
  if (!event) return 0;
  fEvent = event;
  fEventHLT = 0;
  // extract correction object for multiplicity dependence of dEdx
  TObjArray * gainCalibArray = AliTPCcalibDB::Instance()->GetTimeGainSplinesRun(event->GetRunNumber());

  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  if (!transform) {
    AliFatal("Tranformations not in RefitInward");
    return 0;
  }
  transform->SetCurrentRecoParam((AliTPCRecoParam*)AliTPCReconstructor::GetRecoParam());
  const AliTPCRecoParam * recoParam = AliTPCcalibDB::Instance()->GetTransform()->GetCurrentRecoParam();
  Int_t nContribut = event->GetNumberOfTracks();
  TGraphErrors * graphMultDependenceDeDx = 0x0;
  if (recoParam && recoParam->GetUseMultiplicityCorrectionDedx() && gainCalibArray) {
    if (recoParam->GetUseTotCharge()) {
      graphMultDependenceDeDx = (TGraphErrors *) gainCalibArray->FindObject("TGRAPHERRORS_MEANQTOT_MULTIPLICITYDEPENDENCE_BEAM_ALL");
    } else {
      graphMultDependenceDeDx = (TGraphErrors *) gainCalibArray->FindObject("TGRAPHERRORS_MEANQMAX_MULTIPLICITYDEPENDENCE_BEAM_ALL");
    }
  }
  //
  ReadSeeds(event,2);
  fIteration=2;
  //PrepareForProlongation(fSeeds,1);
  PropagateForward2(fSeeds);
  RemoveUsed2(fSeeds,0.4,0.4,20);

  Int_t entriesSeed=fSeeds->GetEntries();
  TObjArray arraySeed(entriesSeed);
  for (Int_t i=0;i<entriesSeed;i++) {
    arraySeed.AddAt(fSeeds->At(i),i);    
  }
  SignShared(&arraySeed);
  //  FindCurling(fSeeds, event,2); // find multi found tracks
  FindSplitted(fSeeds, event,2); // find multi found tracks
  if ((AliTPCReconstructor::StreamLevel()&kStreamFindMultiMC)>0)  FindMultiMC(fSeeds, fEvent,2); // flag: stream MC infomation about the multiple find track (ONLY for MC data)

  Int_t ntracks=0;
  Int_t nseed = fSeeds->GetEntriesFast();
  //
  // RS: the cluster pointers are not permanently attached to the seed during the tracking, need to attach temporarily
  AliTPCclusterMI* seedClusters[kMaxRow];
  int seedsInFriends = 0;
  int seedsInFriendsNorm = event->GetNTPCFriend2Store();
  if (seedsInFriendsNorm>nseed) seedsInFriendsNorm = nseed; // all friends are stored
  //
  fClPointersPoolPtr = fClPointersPool;
  //
  for (Int_t i=0;i<nseed;i++) {
    AliTPCseed * seed = (AliTPCseed*) fSeeds->UncheckedAt(i);
    if (!seed) continue;
    if (seed->GetKinkIndex(0)>0)  UpdateKinkQualityD(seed);  // update quality informations for kinks
    AliESDtrack *esd=event->GetTrack(i);
    //
    //RS: if needed, attach temporary cluster array
    const AliTPCclusterMI** seedClustersSave = seed->GetClusters();
    if (!seedClustersSave) { //RS: temporary attach clusters
      for (int ir=kMaxRow;ir--;) {
	int idx = seed->GetClusterIndex2(ir);
	seedClusters[ir] = idx<0 ? 0 : GetClusterMI(idx);
      }
      seed->SetClustersArrayTMP(seedClusters);
    }
    //
    //
    if (seed->GetNumberOfClusters()<60 && seed->GetNumberOfClusters()<(esd->GetTPCclusters(0) -5)*0.8) {
      AliExternalTrackParam paramIn;
      AliExternalTrackParam paramOut;
      //
      Int_t ncl = seed->RefitTrack(seed,&paramIn,&paramOut);
      //
      if ((AliTPCReconstructor::StreamLevel() & kStreamRecoverIn)>0) {  // flag: stream track information for track  failing in RefitInward function and recovered back
	(*fDebugStreamer)<<"RecoverIn"<< 
	  "seed.="<<seed<<
	  "esd.="<<esd<<
	  "pin.="<<&paramIn<<
	  "pout.="<<&paramOut<<
	  "ncl="<<ncl<<
	  "\n";
      }
      if (ncl>15) {
	seed->Set(paramIn.GetX(),paramIn.GetAlpha(),paramIn.GetParameter(),paramIn.GetCovariance());
	seed->SetNumberOfClusters(ncl);
      }
    }

    seed->PropagateTo(fkParam->GetInnerRadiusLow());
    seed->SetRow(0);
    seed->UpdatePoints();

    AddSystCovariance(seed); // correct covariance matrix for clusters syst. error

    AddCovariance(seed);

    MakeESDBitmaps(seed, esd);
    seed->CookdEdx(0.02,0.6);
    CookLabel(seed,0.1); //For comparison only
    //
    if (((AliTPCReconstructor::StreamLevel()&kStreamRefitInward)>0) && seed!=0) {
      TTreeSRedirector &cstream = *fDebugStreamer;
      cstream<<"RefitInward"<<  // flag: stream track information in RefitInward function (after tracking Iteration 2)
	"Esd.="<<esd<<
	"Track.="<<seed<<
	"\n"; 
    }
    if (seed->GetNumberOfClusters()>15) {
      esd->UpdateTrackParams(seed,AliESDtrack::kTPCrefit); 
      esd->SetTPCPoints(seed->GetPoints());
      esd->SetTPCPointsF(seed->GetNFoundable());
      Int_t   ndedx = seed->GetNCDEDX(0);
      Float_t sdedx = seed->GetSDEDX(0);
      Float_t dedx  = seed->GetdEdx();
      // apply mutliplicity dependent dEdx correction if available
      if (graphMultDependenceDeDx) {
	Double_t corrGain =  AliTPCcalibDButil::EvalGraphConst(graphMultDependenceDeDx, nContribut);
	dedx += (1 - corrGain)*50.; // MIP is normalized to 50
      }
      esd->SetTPCsignal(dedx, sdedx, ndedx);
      //
      // fill new dEdx information
      //
      Double32_t signal[4]; 
      Double32_t signalMax[4]; 
      Char_t ncl[3]; 
      Char_t nrows[3];
      //
      for(Int_t iarr=0;iarr<3;iarr++) {
	signal[iarr] = seed->GetDEDXregion(iarr+1);
	signalMax[iarr] = seed->GetDEDXregion(iarr+5);
	ncl[iarr] = seed->GetNCDEDX(iarr+1);
	nrows[iarr] = seed->GetNCDEDXInclThres(iarr+1);
      }
      signal[3] = seed->GetDEDXregion(4);
      signalMax[3] = seed->GetDEDXregion(8);
      
      //
      AliTPCdEdxInfo * infoTpcPid = new AliTPCdEdxInfo();
      infoTpcPid->SetTPCSignalRegionInfo(signal, ncl, nrows);
      infoTpcPid->SetTPCSignalsQmax(signalMax);
      esd->SetTPCdEdxInfo(infoTpcPid);
      //
      // add seed to the esd track in Calib level
      //
      Bool_t storeFriend = seedsInFriendsNorm>0 && (!esd->GetFriendNotStored()) 
	&& seedsInFriends<(kMaxFriendTracks-1) 
	&& gRandom->Rndm()<(kMaxFriendTracks)/Float_t(seedsInFriendsNorm);
      //      if (AliTPCReconstructor::StreamLevel()>0 &&storeFriend){
      //AliInfoF("Store: %d Stored %d / %d FrOff: %d",storeFriend,seedsInFriends,seedsInFriendsNorm,esd->GetFriendNotStored());
      if (storeFriend){ // RS: seed is needed for calibration, regardless on streamlevel
	seedsInFriends++;
	// RS: this is the only place where the seed is created not in the pool, 
	// since it should belong to ESDevent
	//AliTPCseed * seedCopy = new AliTPCseed(*seed, kTRUE); 
	//esd->AddCalibObject(seedCopy);
	//
	//RS to avoid the cloning the seeds and clusters we will declare the seed to own its
	// clusters and reattach the clusters pointers from the pool, so they are saved in the friends
	seed->SetClusterOwner(kTRUE);
	memcpy(fClPointersPoolPtr,seedClusters,kMaxRow*sizeof(AliTPCclusterMI*));
	seed->SetClustersArrayTMP(fClPointersPoolPtr);
	esd->AddCalibObject(seed);
	fClPointersPoolPtr += kMaxRow;
      }
      else seed->SetClustersArrayTMP((AliTPCclusterMI**)seedClustersSave);
      //
      ntracks++;
    }
    else{
      //printf("problem\n");
    }
    //
    //RS if seed does not own clusters, then it was not added to friends: detach temporary clusters !!!
    if (!seedClustersSave && !seed->GetClusterOwner()) seed->SetClustersArrayTMP(0); 
    //
  }
  //FindKinks(fSeeds,event);
  if ((AliTPCReconstructor::StreamLevel()&kStreamClDump)>0)  DumpClusters(2,fSeeds);  // dump clusters at the end of process (signed with useage flags)
  Info("RefitInward","Number of refitted tracks %d",ntracks);

  AliCosmicTracker::FindCosmic(event, kTRUE);

  FillClusterOccupancyInfo();

  return 0;
}


Int_t AliTPCtracker::PropagateBack(AliESDEvent *event)
{
  //
  // back propagation of ESD tracks
  //
  if (!event) return 0;
  fEvent = event;
  fEventHLT = 0;
  fIteration = 1;
  ReadSeeds(event,1);
  PropagateBack(fSeeds); 
  RemoveUsed2(fSeeds,0.4,0.4,20);
  //FindCurling(fSeeds, fEvent,1);  
  FindSplitted(fSeeds, event,1); // find multi found tracks
  if ((AliTPCReconstructor::StreamLevel()&kStreamFindMultiMC)>0)  FindMultiMC(fSeeds, fEvent,1); // find multi found tracks
  //
  // RS: the cluster pointers are not permanently attached to the seed during the tracking, need to attach temporarily
  AliTPCclusterMI* seedClusters[kMaxRow];
  //
  Int_t nseed = fSeeds->GetEntriesFast();
  Int_t ntracks=0;
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed * seed = (AliTPCseed*) fSeeds->UncheckedAt(i);
    if (!seed) continue;
    AliESDtrack *esd=event->GetTrack(i);
    if (!esd) continue; //never happen
    //
    //RS: if needed, attach temporary cluster array
    const AliTPCclusterMI** seedClustersSave = seed->GetClusters();
    if (!seedClustersSave) { //RS: temporary attach clusters
      for (int ir=kMaxRow;ir--;) {
	int idx = seed->GetClusterIndex2(ir);
	seedClusters[ir] = idx<0 ? 0 : GetClusterMI(idx);
      }
      seed->SetClustersArrayTMP(seedClusters);
    }
    //
    if (seed->GetKinkIndex(0)<0)  UpdateKinkQualityM(seed);  // update quality informations for kinks
    seed->UpdatePoints();
    
    AddSystCovariance(seed); // correct covariance matrix for clusters syst. error

    AddCovariance(seed);
    if (seed->GetNumberOfClusters()<60 && seed->GetNumberOfClusters()<(esd->GetTPCclusters(0) -5)*0.8){
      AliExternalTrackParam paramIn;
      AliExternalTrackParam paramOut;
      Int_t ncl = seed->RefitTrack(seed,&paramIn,&paramOut);
      if ((AliTPCReconstructor::StreamLevel()&kStreamRecoverBack)>0) { // flag: stream track information for track  faling PropagateBack function and recovered back (RefitTrack)
	(*fDebugStreamer)<<"RecoverBack"<<
	  "seed.="<<seed<<
	  "esd.="<<esd<<
	  "pin.="<<&paramIn<<
	  "pout.="<<&paramOut<<
	  "ncl="<<ncl<<
	  "\n";
      }
      if (ncl>15) {
	seed->Set(paramOut.GetX(),paramOut.GetAlpha(),paramOut.GetParameter(),paramOut.GetCovariance());
	seed->SetNumberOfClusters(ncl);
      }
    }
    seed->CookdEdx(0.02,0.6);
    CookLabel(seed,0.1); //For comparison only
    if (seed->GetNumberOfClusters()>15){
      esd->UpdateTrackParams(seed,AliESDtrack::kTPCout);
      esd->SetTPCPoints(seed->GetPoints());
      esd->SetTPCPointsF(seed->GetNFoundable());
      Int_t   ndedx = seed->GetNCDEDX(0);
      Float_t sdedx = seed->GetSDEDX(0);
      Float_t dedx  = seed->GetdEdx();
      esd->SetTPCsignal(dedx, sdedx, ndedx);
      ntracks++;
      Int_t eventnumber = event->GetEventNumberInFile();// patch 28 fev 06
      // This is most likely NOT the event number you'd like to use. It has nothing to do with the 'real' event number      
      if (((AliTPCReconstructor::StreamLevel()&kStreamCPropagateBack)>0) && esd) {
	(*fDebugStreamer)<<"PropagateBack"<<  // flag: stream track information in PropagateBack function (after tracking Iteration 1)
	  "Tr0.="<<seed<<
	  "esd.="<<esd<<
	  "EventNrInFile="<<eventnumber<<
	  "\n";       
      }
    }
    //
    if (!seedClustersSave) seed->SetClustersArrayTMP(0); //RS detach temporary clusters !!!
    //
  }
  if (AliTPCReconstructor::StreamLevel()&kStreamClDump)  DumpClusters(1,fSeeds); // flag:stream clusters at the end of process (signed with usage flag)
  //FindKinks(fSeeds,event);
  Info("PropagateBack","Number of back propagated tracks %d",ntracks);
  fEvent =0;
  fEventHLT = 0;

  return 0;
}


Int_t AliTPCtracker::PostProcess(AliESDEvent *event)
{
  //
  // Post process events 
  //
  if (!event) return 0;

  //
  // Set TPC event status
  // 

  // event affected by HV dip
  // reset TPC status
  if(IsTPCHVDipEvent(event)) { 
    event->ResetDetectorStatus(AliDAQ::kTPC);
  }
 
  //printf("Status %d \n", event->IsDetectorOn(AliDAQ::kTPC));

  return 0;
}


 void AliTPCtracker::DeleteSeeds()
{
  //
  fSeeds->Clear();
  ResetSeedsPool();
  //  delete fSeeds; // RS avoid unnecessary deletes
  // fSeeds =0;
}

void AliTPCtracker::ReadSeeds(const AliESDEvent *const event, Int_t direction)
{
  //
  //read seeds from the event
  
  Int_t nentr=event->GetNumberOfTracks();
  if (fDebug>0){
    Info("ReadSeeds", "Number of ESD tracks: %d\n", nentr);
  }
  if (fSeeds) 
    DeleteSeeds();
  if (!fSeeds) fSeeds = new TObjArray(nentr);
  else if (fSeeds->GetSize()<nentr) fSeeds->Expand(nentr);  //RS reuse array
  //
  UnsignClusters();
  //  Int_t ntrk=0;  
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);
    ULong_t status=esd->GetStatus();
    if (!(status&AliESDtrack::kTPCin)) continue;
    AliTPCtrack t(*esd);
    t.SetNumberOfClusters(0);
    //    AliTPCseed *seed = new AliTPCseed(t,t.GetAlpha());
    AliTPCseed *seed = new( NextFreeSeed() ) AliTPCseed(t/*,t.GetAlpha()*/);
    seed->SetPoolID(fLastSeedID);
    seed->SetUniqueID(esd->GetID());
    AddCovariance(seed);   //add systematic ucertainty
    for (Int_t ikink=0;ikink<3;ikink++) {
      Int_t index = esd->GetKinkIndex(ikink);
      seed->GetKinkIndexes()[ikink] = index;
      if (index==0) continue;
      index = TMath::Abs(index);
      AliESDkink * kink = fEvent->GetKink(index-1);
      if (kink&&esd->GetKinkIndex(ikink)<0){
	if ((status & AliESDtrack::kTRDrefit) != 0) kink->SetStatus(1,2);
	if ((status & AliESDtrack::kITSout) != 0)   kink->SetStatus(1,0);
      }
      if (kink&&esd->GetKinkIndex(ikink)>0){
	if ((status & AliESDtrack::kTRDrefit) != 0) kink->SetStatus(1,6);
	if ((status & AliESDtrack::kITSout) != 0)   kink->SetStatus(1,4);
      }

    }
    // RS: resetting is done in the AliTPCtrack constructor 
    // if (((status&AliESDtrack::kITSout)==0)&&(direction==1)) seed->ResetCovariance(10.); 
    // if ( direction ==2 &&(status & AliESDtrack::kTRDrefit) == 0 ) seed->ResetCovariance(10.);
    //if ( direction ==2 && ((status & AliESDtrack::kTPCout) == 0) ) {
    //  fSeeds->AddAt(0,i);
    //  MarkSeedFree( seed );
    //  continue;    
    //}
    
    //
    //
    // rotate to the local coordinate system
    //   
    fSectors=fInnerSec; fN=fkNIS;    
    Double_t alpha=seed->GetAlpha();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN;
    alpha =ns*fSectors->GetAlpha() + fSectors->GetAlphaShift();
    alpha-=seed->GetAlpha();  
    if (alpha<-TMath::Pi()) alpha += 2*TMath::Pi();
    if (alpha>=TMath::Pi()) alpha -= 2*TMath::Pi();
    if (TMath::Abs(alpha) > 0.001) { // This should not happen normally
      AliWarning(Form("Rotating track over %f",alpha));
      if (!seed->Rotate(alpha)) {
	MarkSeedFree( seed );
	continue;
      }
    }
    seed->SetESD(esd);
    // sign clusters
    if (esd->GetKinkIndex(0)<=0){
      for (Int_t irow=0;irow<kMaxRow;irow++){
	Int_t index = seed->GetClusterIndex2(irow);    
	if (index >= 0){ 
	  //
	  AliTPCclusterMI * cl = GetClusterMI(index);
	  //RS seed->SetClusterPointer(irow,cl);
	  if (cl){
	    if ((index & 0x8000)==0){
	      cl->Use(10);  // accepted cluster	  
	    }else{
	      cl->Use(6);   // close cluster not accepted
	    }	
	  }else{
	    Info("ReadSeeds","Not found cluster");
	  }
	}
      }
    }
    fSeeds->AddAt(seed,i);
  }
}

//_____________________________________________________________________________
void AliTPCtracker::MakeSeeds3(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2,  Float_t cuts[4],
				 Float_t deltay, Int_t ddsec) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  // SEEDING WITH VERTEX CONSTRAIN 
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut

  const double kRoadY = 1., kRoadZ = 0.6;

  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;

  Double_t x[5], c[15];
  //  Int_t di = i1-i2;
  //
  AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed();
  seed->SetPoolID(fLastSeedID);
  Double_t alpha=fSectors->GetAlpha(), shift=fSectors->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);
  //
  //  Double_t x1 =fOuterSec->GetX(i1);
  //Double_t xx2=fOuterSec->GetX(i2);
  
  Double_t x1 =GetXrow(i1);
  Double_t xx2=GetXrow(i2);

  Double_t x3=GetX(), y3=GetY(), z3=GetZ();

  Int_t imiddle = (i2+i1)/2;    //middle pad row index
  Double_t xm = GetXrow(imiddle); // radius of middle pad-row
  const AliTPCtrackerRow& krm=GetRow(sec,imiddle); //middle pad -row
  //
  Int_t ns =sec;   

  const AliTPCtrackerRow& kr1=GetRow(ns,i1);
  Double_t ymax  = GetMaxY(i1)-kr1.GetDeadZone()-1.5;  
  Double_t ymaxm = GetMaxY(imiddle)-kr1.GetDeadZone()-1.5;  

  //
  // change cut on curvature if it can't reach this layer
  // maximal curvature set to reach it
  Double_t dvertexmax  = TMath::Sqrt((x1-x3)*(x1-x3)+(ymax+5-y3)*(ymax+5-y3));
  if (dvertexmax*0.5*cuts[0]>0.85){
    cuts[0] = 0.85/(dvertexmax*0.5+1.);
  }
  Double_t r2min = 1/(cuts[0]*cuts[0]);  //minimal square of radius given by cut

  //  Int_t ddsec = 1;
  if (deltay>0) ddsec = 0; 
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    if (kr1[is]->IsUsed(10)) continue;
    if (kr1[is]->IsDisabled()) {
      continue;
    }

    Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();    
    //if (TMath::Abs(y1)>ymax) continue;

    if (deltay>0 && TMath::Abs(ymax-TMath::Abs(y1))> deltay ) continue;  // seed only at the edge

    // find possible directions    
    Float_t anglez = (z1-z3)/(x1-x3); 
    Float_t extraz = z1 - anglez*(x1-xx2);  // extrapolated z      
    //
    //
    //find   rotation angles relative to line given by vertex and point 1
    Double_t dvertex2 = (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3);
    Double_t dvertex  = TMath::Sqrt(dvertex2);
    Double_t angle13  = TMath::ATan((y1-y3)/(x1-x3));
    Double_t cs13     = cos(-angle13), sn13 = sin(-angle13);            
    
    //
    // loop over 2 sectors
    Int_t dsec1=-ddsec;
    Int_t dsec2= ddsec;
    if (y1<0)  dsec2= 0;
    if (y1>0)  dsec1= 0;
    
    Double_t dddz1=0;  // direction of delta inclination in z axis
    Double_t dddz2=0;
    if ( (z1-z3)>0)
      dddz1 =1;    
    else
      dddz2 =1;
    //
    for (Int_t dsec = dsec1; dsec<=dsec2;dsec++){
      Int_t sec2 = sec + dsec;
      // 
      //      AliTPCtrackerRow&  kr2  = fOuterSec[(sec2+fkNOS)%fkNOS][i2];
      //AliTPCtrackerRow&  kr2m = fOuterSec[(sec2+fkNOS)%fkNOS][imiddle];
      AliTPCtrackerRow&  kr2  = GetRow((sec2+fkNOS)%fkNOS,i2);
      AliTPCtrackerRow&  kr2m = GetRow((sec2+fkNOS)%fkNOS,imiddle);
      Int_t  index1 = TMath::Max(kr2.Find(extraz-0.6-dddz1*TMath::Abs(z1)*0.05)-1,0);
      Int_t  index2 = TMath::Min(kr2.Find(extraz+0.6+dddz2*TMath::Abs(z1)*0.05)+1,kr2);

      // rotation angles to p1-p3
      Double_t cs13r     = cos(-angle13+dsec*alpha)/dvertex, sn13r = sin(-angle13+dsec*alpha)/dvertex;            
      Double_t x2,   y2,   z2; 
      //
      //      Double_t dymax = maxangle*TMath::Abs(x1-xx2);

      //
      Double_t dxx0 =  (xx2-x3)*cs13r;
      Double_t dyy0 =  (xx2-x3)*sn13r;
      for (Int_t js=index1; js < index2; js++) {
	const AliTPCclusterMI *kcl = kr2[js];
	if (kcl->IsUsed(10)) continue; 	
	if (kcl->IsDisabled()) {
	  continue;
	}
	//
	//calcutate parameters
	//	
	Double_t yy0 =  dyy0 +(kcl->GetY()-y3)*cs13r;
	// stright track
	if (TMath::Abs(yy0)<0.000001) continue;
	Double_t xx0 =  dxx0 -(kcl->GetY()-y3)*sn13r;
	Double_t y0  =  0.5*(xx0*xx0+yy0*yy0-xx0)/yy0;
	Double_t r02 = (0.25+y0*y0)*dvertex2;	
	//curvature (radius) cut
	if (r02<r2min) continue;		
       
	nin0++;	
	//
	Double_t c0  = 1/TMath::Sqrt(r02);
	if (yy0>0) c0*=-1.;	
	       
       
	//Double_t dfi0   = 2.*TMath::ASin(dvertex*c0*0.5);
	//Double_t dfi1   = 2.*TMath::ASin(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);
	Double_t dfi0   = 2.*asinf(dvertex*c0*0.5);
	Double_t dfi1   = 2.*asinf(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);  
	//
	//
	Double_t z0  =  kcl->GetZ();  
	Double_t zzzz2    = z1-(z1-z3)*dfi1/dfi0;
	if (TMath::Abs(zzzz2-z0)>0.5) continue;       
	nin1++;              
	//	
	Double_t dip    = (z1-z0)*c0/dfi1;        
	Double_t x0 = (0.5*cs13+y0*sn13)*dvertex*c0;
	//
	y2 = kcl->GetY(); 
	if (dsec==0){
	  x2 = xx2; 
	  z2 = kcl->GetZ();	  
	}
	else
	  {
	    // rotation	
	    z2 = kcl->GetZ();  
	    x2= xx2*cs-y2*sn*dsec;
	    y2=+xx2*sn*dsec+y2*cs;
	  }
	
	x[0] = y1;
	x[1] = z1;
	x[2] = x0;
	x[3] = dip;
	x[4] = c0;
	//
	//
	// do we have cluster at the middle ?
	Double_t ym=0,zm=0;
	if (!GetProlongation(x1,xm,x,ym,zm)) continue;
	UInt_t dummy; 
	AliTPCclusterMI * cm=0;
	if (TMath::Abs(ym)-ymaxm<0){	  
	  cm = krm.FindNearest2(ym,zm,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ,dummy);
	  if ((!cm) || (cm->IsUsed(10))) {	  
	    continue;
	  }
	}
	else{	  
	  // rotate y1 to system 0
	  // get state vector in rotated system 
	  Double_t yr1  = (-0.5*sn13+y0*cs13)*dvertex*c0;
	  Double_t xr2  =  x0*cs+yr1*sn*dsec;
	  Double_t xr[5]={kcl->GetY(),kcl->GetZ(), xr2, dip, c0};
	  //
	  if (!GetProlongation(xx2,xm,xr,ym,zm)) continue;
	  if (TMath::Abs(ym)-ymaxm<0){
	    cm = kr2m.FindNearest2(ym,zm,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ,dummy);
	    if ((!cm) || (cm->IsUsed(10))) {	  
	      continue;
	    }
	  }
	}
       

	// Double_t dym = 0;
	// Double_t dzm = 0;
	// if (cm){
	//   dym = ym - cm->GetY();
	//   dzm = zm - cm->GetZ();
	// }
	nin2++;


	//
	//
        Double_t sy1=kr1[is]->GetSigmaY2()*2., sz1=kr1[is]->GetSigmaZ2()*2.;
        Double_t sy2=kcl->GetSigmaY2()*2.,     sz2=kcl->GetSigmaZ2()*2.;
	//Double_t sy3=400*3./12., sy=0.1, sz=0.1;
	Double_t sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	//Double_t sy3=25000*x[4]*x[4]*60+0.5, sy=0.1, sz=0.1;

	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
	
	Double_t f30=(F3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
	
        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
        c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
                       c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
        c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
        c[13]=f30*sy1*f40+f32*sy2*f42;
        c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	//	if (!BuildSeed(kr1[is],kcl,0,x1,x2,x3,x,c)) continue;
	
        UInt_t index=kr1.GetIndex(is);
	if (seed) {MarkSeedFree(seed); seed = 0;}
	AliTPCseed *track = seed = new( NextFreeSeed() ) AliTPCseed(x1, ns*alpha+shift, x, c, index);
	seed->SetPoolID(fLastSeedID);
	seed->SetRow(i1); //RS: memorise current row
	track->SetIsSeeding(kTRUE);
	track->SetSeed1(i1);
	track->SetSeed2(i2);
	track->SetSeedType(3);

       
	//if (dsec==0) {
	  FollowProlongation(*track, (i1+i2)/2,1);
	  Int_t foundable,found,shared;	  
	  GetSeedClusterStatistic(track,(i1+i2)/2,i1, found, foundable, shared, kTRUE); //RS: seeds don't keep their clusters
	  //RS track->GetClusterStatistic((i1+i2)/2,i1, found, foundable, shared, kTRUE);
	  if ((found<0.55*foundable)  || shared>0.5*found || (track->GetSigmaY2()+track->GetSigmaZ2())>0.5){
	    MarkSeedFree(seed); seed = 0;
	    continue;
	  }
	  //}
	
	nin++;
	FollowProlongation(*track, i2,1);
	
	
	//Int_t rc = 1;
	track->SetBConstrain(1);
	//	track->fLastPoint = i1+fInnerSec->GetNRows();  // first cluster in track position
	track->SetLastPoint(i1);  // first cluster in track position
	track->SetFirstPoint(track->GetLastPoint());
	
	if (track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	    track->GetNumberOfClusters() < track->GetNFoundable()*0.6 || 
	    track->GetNShared()>0.4*track->GetNumberOfClusters() ) {
	  MarkSeedFree(seed); seed = 0;
	  continue;
	}
	nout1++;
	Double_t zv, bz=GetBz();
        if ( !track->GetZAt(0.,bz,zv) )       { MarkSeedFree( seed ); seed = 0; continue; }
	//
	if (fDisableSecondaries) {
	  if (TMath::Abs(zv)>fPrimaryDCAZCut) { MarkSeedFree( seed ); seed = 0; continue; }
	  double yv; 
	  if ( !track->GetZAt(0.,bz,yv) )     { MarkSeedFree( seed ); seed = 0; continue; }
	  if (TMath::Abs(zv)>fPrimaryDCAZCut) { MarkSeedFree( seed ); seed = 0; continue; }
	}
	
	//
        // Z VERTEX CONDITION
	if (TMath::Abs(zv-z3)>cuts[2]) {
	  FollowProlongation(*track, TMath::Max(i2-20,0));
          if ( !track->GetZAt(0.,bz,zv) ) continue;
	  if (TMath::Abs(zv-z3)>cuts[2]){
	    FollowProlongation(*track, TMath::Max(i2-40,0));
            if ( !track->GetZAt(0.,bz,zv) ) continue;
	    if (TMath::Abs(zv-z3)>cuts[2] &&(track->GetNumberOfClusters() > track->GetNFoundable()*0.7)){
	      // make seed without constrain
	      AliTPCseed * track2 = MakeSeed(track,0.2,0.5,1.);
	      FollowProlongation(*track2, i2,1);
	      track2->SetBConstrain(kFALSE);
	      track2->SetSeedType(1);
	      arr->AddLast(track2); 
	      MarkSeedFree( seed ); seed = 0;
	      continue;		
	    }
	    else{
	      MarkSeedFree( seed ); seed = 0;
	      continue;
	    
	    }
	  }
	}
      
	track->SetSeedType(0);
	arr->AddLast(track); // note, track is seed, don't free the seed
	seed = new( NextFreeSeed() ) AliTPCseed; 	
	seed->SetPoolID(fLastSeedID);
	nout2++;
	// don't consider other combinations
	if (track->GetNumberOfClusters() > track->GetNFoundable()*0.8)
	  break;
      }
    }
  }
  if (fDebug>3){
    Info("MakeSeeds3","\nSeeding statistic:\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2);
  }
  if (seed) MarkSeedFree( seed );
}

//_____________________________________________________________________________
void AliTPCtracker::MakeSeeds3Dist(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2,  Float_t cuts[4],
				   Float_t deltay, Int_t ddsec) {
  //-----------------------------------------------------------------
  // This function creates track seeds, accounting for distortions
  // SEEDING WITH VERTEX CONSTRAIN 
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut
  const double kRoadY = 1., kRoadZ = 0.6;

  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;

  Double_t x[5], c[15];
  //  Int_t di = i1-i2;
  //
  AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed();
  seed->SetPoolID(fLastSeedID);
  Double_t alpha=fSectors->GetAlpha(), shift=fSectors->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);
  //
  Double_t x1def = GetXrow(i1);
  Double_t xx2def = GetXrow(i2);

  Double_t x3=GetX(), y3=GetY(), z3=GetZ();

  Int_t imiddle = (i2+i1)/2;    //middle pad row index
  const AliTPCtrackerRow& krm=GetRow(sec,imiddle); //middle pad -row
  //
  Int_t ns =sec;   

  const AliTPCtrackerRow& kr1=GetRow(ns,i1);
  Double_t ymax  = GetMaxY(i1)-kr1.GetDeadZone()-1.5;  
  Double_t ymaxm = GetMaxY(imiddle)-kr1.GetDeadZone()-1.5; 

  //
  // change cut on curvature if it can't reach this layer
  // maximal curvature set to reach it
  Double_t dvertexmax  = TMath::Sqrt((x1def-x3)*(x1def-x3)+(ymax+5-y3)*(ymax+5-y3));
  if (dvertexmax*0.5*cuts[0]>0.85){
    cuts[0] = 0.85/(dvertexmax*0.5+1.);
  }
  Double_t r2min = 1/(cuts[0]*cuts[0]);  //minimal square of radius given by cut

  //  Int_t ddsec = 1;
  if (deltay>0) ddsec = 0; 
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    if (kr1[is]->IsUsed(10)) continue;
    if (kr1[is]->IsDisabled()) {
      continue;
    }
    const AliTPCclusterMI* clkr1 = kr1[is];
    Double_t x1=clkr1->GetX(), y1=clkr1->GetY(), z1=clkr1->GetZ();    
    //if (TMath::Abs(y1)>ymax) continue;
    double y1EdgeDist = y1;
    if  (fAccountDistortions) y1EdgeDist -= GetYSectEdgeDist(sec,ns,y1,z1);
    if (deltay>0) {
      double margin = (y1>0 ? ymax-y1EdgeDist : ymax + y1EdgeDist);
      if (margin<deltay ) continue;  // seed only at the edge
    }
    //
    // find possible directions
    double dx13 = x1-x3, dy13 = y1-y3, dz13 = z1-z3;
    double anglez = dz13/dx13;
    double extraz = z1 - anglez*(x1-xx2def);  // extrapolated z      
    //
    //
    //find   rotation angles relative to line given by vertex and point 1
    Double_t dvertex2 = dx13*dx13+dy13*dy13;
    Double_t dvertex  = TMath::Sqrt(dvertex2);
    Double_t angle13  = TMath::ATan(dy13/dx13);
    Double_t cs13     = cos(-angle13), sn13 = sin(-angle13);            
    
    //
    // loop over 2 sectors
    Int_t dsec1=-ddsec;
    Int_t dsec2= ddsec;
    if (y1<0)  dsec2= 0;
    if (y1>0)  dsec1= 0;
    
    Double_t dddz1=0;  // direction of delta inclination in z axis
    Double_t dddz2=0;
    if ( dz13>0 ) dddz1 =1;    
    else          dddz2 =1;
    //
    for (Int_t dsec = dsec1; dsec<=dsec2;dsec++){
      Int_t sec2 = sec + dsec;
      // 
      //      AliTPCtrackerRow&  kr2  = fOuterSec[(sec2+fkNOS)%fkNOS][i2];
      //AliTPCtrackerRow&  kr2m = fOuterSec[(sec2+fkNOS)%fkNOS][imiddle];
      AliTPCtrackerRow&  kr2  = GetRow((sec2+fkNOS)%fkNOS,i2);
      AliTPCtrackerRow&  kr2m = GetRow((sec2+fkNOS)%fkNOS,imiddle);
      Int_t  index1 = TMath::Max(kr2.Find(extraz-0.6-dddz1*TMath::Abs(z1)*0.05)-1,0);
      Int_t  index2 = TMath::Min(kr2.Find(extraz+0.6+dddz2*TMath::Abs(z1)*0.05)+1,kr2);

      // rotation angles to p1-p3
      Double_t cs13r     = cos(-angle13+dsec*alpha)/dvertex, sn13r = sin(-angle13+dsec*alpha)/dvertex;            
      //
      //      Double_t dymax = maxangle*TMath::Abs(x1-xx2);
      //
      for (Int_t js=index1; js < index2; js++) {
	const AliTPCclusterMI *kcl = kr2[js];
	if (kcl->IsUsed(10)) continue; 	
	if (kcl->IsDisabled()) {
	  continue;
	}
	//
	double x2=kcl->GetX(), y2=kcl->GetY(), z2=kcl->GetZ();
	double dy23 = y2-y3, dx23 = x2-x3;
	//calcutate parameters
	Double_t dxx0 =  dx23*cs13r;
	Double_t dyy0 =  dx23*sn13r;
	//	
	Double_t yy0 =  dyy0 +dy23*cs13r;
	// stright track
	if (TMath::Abs(yy0)<0.000001) continue;
	Double_t xx0 =  dxx0 -dy23*sn13r;
	Double_t y0  =  0.5*(xx0*xx0+yy0*yy0-xx0)/yy0;
	Double_t r02 = (0.25+y0*y0)*dvertex2;	
	//curvature (radius) cut
	if (r02<r2min) continue;		
       
	nin0++;	
	//
	Double_t c0  = 1/TMath::Sqrt(r02);
	if (yy0>0) c0*=-1.;	
       
	Double_t dfi0   = 2.*TMath::ASin(dvertex*c0*0.5);
	Double_t dfi1   = 2.*TMath::ASin(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);
	//	Double_t dfi0   = 2.*AliTPCFastMath::FastAsin(dvertex*c0*0.5);
	//	Double_t dfi1   = 2.*AliTPCFastMath::FastAsin(TMath::Sqrt(yy0*yy0+(1-xx0)*(1-xx0))*dvertex*c0*0.5);  
	//
	Double_t zzzz2    = z1-dz13*dfi1/dfi0;
	if (TMath::Abs(zzzz2-z2)>0.5) continue;       
	nin1++;              
	//	
	Double_t dip    = (z1-z2)*c0/dfi1;        
	Double_t x0 = (0.5*cs13+y0*sn13)*dvertex*c0;
	//
	if (dsec!=0) {
	  // rotation
	  double xt2 = x2;
	  x2= xt2*cs-y2*sn*dsec;
	  y2=+xt2*sn*dsec+y2*cs;
	}
	
	x[0] = y1;
	x[1] = z1;
	x[2] = x0;
	x[3] = dip;
	x[4] = c0;
	//
	//
	// do we have cluster at the middle ?
	Double_t xm = GetXrow(imiddle), ym=0, zm=0; // radius of middle pad-row
	if (!GetProlongation(x1,xm,x,ym,zm)) continue;
	// account for distortion
	double dxDist = GetDistortionX(xm,ym,zm,sec,imiddle);
	if (TMath::Abs(dxDist)>0.05 && !GetProlongation(x1,xm+dxDist,x,ym,zm)) continue; //RS:? can we use straight line here?
	UInt_t dummy; 
	AliTPCclusterMI * cm=0;
	double ymEdgeDist = ym;
	if (fAccountDistortions) ymEdgeDist -= GetYSectEdgeDist(sec,imiddle,ym,zm); // ym shifted by edge distortion
	if ( (ym>0&&ymEdgeDist<ymaxm) || (ym<=0&&ymEdgeDist>-ymaxm) ) { //RS  //	if (TMath::Abs(ym)-ymaxm<0){
	  cm = krm.FindNearest2(ym,zm,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ,dummy);
	  if ((!cm) || (cm->IsUsed(10))) continue;
	}
	else{	  
	  // rotate y1 to system 0
	  // get state vector in rotated system 
	  Double_t yr1  = (-0.5*sn13+y0*cs13)*dvertex*c0;
	  Double_t xr2  =  x0*cs+yr1*sn*dsec;
	  Double_t xr[5]={kcl->GetY(),kcl->GetZ(), xr2, dip, c0};
	  //
	  if (!GetProlongation(kcl->GetX(),xm,xr,ym,zm)) continue;
	  double dxDist = GetDistortionX(xm,ym,zm,sec2,imiddle);
	  if (TMath::Abs(dxDist)>0.05 && !GetProlongation(x1,xm+dxDist,x,ym,zm)) continue; //RS:? can we use straight line here?
	  //
	  ymEdgeDist = ym;
	  if (fAccountDistortions) ymEdgeDist -= GetYSectEdgeDist(sec2,imiddle,ym,zm); // ym shifted by edge distortion
	  if ( (ym>0&&ymEdgeDist<ymaxm) || (ym<=0&&ymEdgeDist>-ymaxm) ) { //RS //if (TMath::Abs(ym)-ymaxm<0){
	    cm = kr2m.FindNearest2(ym,zm,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ,dummy);
	    if ((!cm) || (cm->IsUsed(10))) continue;
	  }
	}
       
	// Double_t dym = 0;
	// Double_t dzm = 0;
	// if (cm){
	//   dym = ym - cm->GetY();
	//   dzm = zm - cm->GetZ();
	// }
	nin2++;


	//
	//
        Double_t sy1=kr1[is]->GetSigmaY2()*2., sz1=kr1[is]->GetSigmaZ2()*2.;
        Double_t sy2=kcl->GetSigmaY2()*2.,     sz2=kcl->GetSigmaZ2()*2.;
	//Double_t sy3=400*3./12., sy=0.1, sz=0.1;
	Double_t sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	//Double_t sy3=25000*x[4]*x[4]*60+0.5, sy=0.1, sz=0.1;

	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
	
	Double_t f30=(F3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
	
        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
        c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
                       c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
        c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
        c[13]=f30*sy1*f40+f32*sy2*f42;
        c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	//	if (!BuildSeed(kr1[is],kcl,0,x1,x2,x3,x,c)) continue;
	
        UInt_t index=kr1.GetIndex(is);
	if (seed) {MarkSeedFree(seed); seed = 0;}
	AliTPCseed *track = seed = new( NextFreeSeed() ) AliTPCseed(x1, ns*alpha+shift, x, c, index);
	seed->SetPoolID(fLastSeedID);
	seed->SetRow(i1); //RS: memorise current row
	track->SetIsSeeding(kTRUE);
	track->SetSeed1(i1);
	track->SetSeed2(i2);
	track->SetSeedType(3);

       
	//if (dsec==0) {
	  FollowProlongation(*track, (i1+i2)/2,1);
	  Int_t foundable,found,shared;
	  GetSeedClusterStatistic(track,(i1+i2)/2,i1, found, foundable, shared, kTRUE); //RS: seeds don't keep their clusters
	  //RS track->GetClusterStatistic((i1+i2)/2,i1, found, foundable, shared, kTRUE);
	  if ((found<0.55*foundable)  || shared>0.5*found) {
	    //RS: with distortions related cluster errors the track error may grow, don't use this cut
	    //|| (track->GetSigmaY2()+track->GetSigmaZ2())>0.5){
	    MarkSeedFree(seed); seed = 0;
	    continue;
	  }
	  //}
	
	nin++;
	FollowProlongation(*track, i2,1);
	
	
	//Int_t rc = 1;
	track->SetBConstrain(1);
	//	track->fLastPoint = i1+fInnerSec->GetNRows();  // first cluster in track position
	track->SetLastPoint(i1);  // first cluster in track position
	track->SetFirstPoint(track->GetLastPoint());
	
	if (track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	    track->GetNumberOfClusters() < track->GetNFoundable()*0.6 || 
	    track->GetNShared()>0.4*track->GetNumberOfClusters() ) {
	  MarkSeedFree(seed); seed = 0;
	  continue;
	}
	nout1++;
	Double_t zv, bz=GetBz();
        if ( !track->GetZAt(0.,bz,zv) )       { MarkSeedFree( seed ); seed = 0; continue; }
	//
	if (fDisableSecondaries) {
	  if (TMath::Abs(zv)>fPrimaryDCAZCut) { MarkSeedFree( seed ); seed = 0; continue; }
	  double yv; 
	  if ( !track->GetZAt(0.,bz,yv) )     { MarkSeedFree( seed ); seed = 0; continue; }
	  if (TMath::Abs(zv)>fPrimaryDCAZCut) { MarkSeedFree( seed ); seed = 0; continue; }
	}
	
	//
        // Z VERTEX CONDITION
	if (TMath::Abs(zv-z3)>cuts[2]) {
	  FollowProlongation(*track, TMath::Max(i2-20,0));
          if ( !track->GetZAt(0.,bz,zv) ) continue;
	  if (TMath::Abs(zv-z3)>cuts[2]){
	    FollowProlongation(*track, TMath::Max(i2-40,0));
            if ( !track->GetZAt(0.,bz,zv) ) continue;
	    if (TMath::Abs(zv-z3)>cuts[2] &&(track->GetNumberOfClusters() > track->GetNFoundable()*0.7)){
	      // make seed without constrain
	      AliTPCseed * track2 = MakeSeed(track,0.2,0.5,1.);
	      FollowProlongation(*track2, i2,1);
	      track2->SetBConstrain(kFALSE);
	      track2->SetSeedType(1);
	      arr->AddLast(track2); 
	      MarkSeedFree( seed ); seed = 0;
	      continue;		
	    }
	    else{
	      MarkSeedFree( seed ); seed = 0;
	      continue;
	    
	    }
	  }
	}
      
	track->SetSeedType(0);
	arr->AddLast(track); // note, track is seed, don't free the seed
	seed = new( NextFreeSeed() ) AliTPCseed; 	
	seed->SetPoolID(fLastSeedID);
	nout2++;
	// don't consider other combinations
	if (track->GetNumberOfClusters() > track->GetNFoundable()*0.8)
	  break;
      }
    }
  }
  if (fDebug>3){
    Info("MakeSeeds3Dist","\nSeeding statistic:\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2);
  }
  if (seed) MarkSeedFree( seed );
}


void AliTPCtracker::MakeSeeds5(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2,  Float_t cuts[4],
				 Float_t deltay) {
  


  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut


  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;
  Int_t nout3 =0;
  Double_t x[5], c[15];
  //
  // make temporary seed
  AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed;
  seed->SetPoolID(fLastSeedID);
  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  //
  //

  // first 3 padrows
  Double_t x1 = GetXrow(i1-1);
  const    AliTPCtrackerRow& kr1=GetRow(sec,i1-1);
  Double_t y1max  = GetMaxY(i1-1)-kr1.GetDeadZone()-1.5;  
  //
  Double_t x1p = GetXrow(i1);
  const    AliTPCtrackerRow& kr1p=GetRow(sec,i1);
  //
  Double_t x1m = GetXrow(i1-2);
  const    AliTPCtrackerRow& kr1m=GetRow(sec,i1-2);

  //
  //last 3 padrow for seeding
  AliTPCtrackerRow&  kr3  = GetRow((sec+fkNOS)%fkNOS,i1-7);
  Double_t    x3   =  GetXrow(i1-7);
  //  Double_t    y3max= GetMaxY(i1-7)-kr3.fDeadZone-1.5;  
  //
  AliTPCtrackerRow&  kr3p  = GetRow((sec+fkNOS)%fkNOS,i1-6);
  Double_t    x3p   = GetXrow(i1-6);
  //
  AliTPCtrackerRow&  kr3m  = GetRow((sec+fkNOS)%fkNOS,i1-8);
  Double_t    x3m   = GetXrow(i1-8);

  //
  //
  // middle padrow
  Int_t im = i1-4;                           //middle pad row index
  Double_t xm         = GetXrow(im);         // radius of middle pad-row
  const AliTPCtrackerRow& krm=GetRow(sec,im);   //middle pad -row
  //  Double_t ymmax = GetMaxY(im)-kr1.fDeadZone-1.5;  
  //
  //
  Double_t deltax  = x1-x3;
  Double_t dymax   = deltax*cuts[1];
  Double_t dzmax   = deltax*cuts[3];
  //
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    if (kr1[is]->IsUsed(10)) continue;	
    if (kr1[is]->IsDisabled()) {
      continue;
    }

    Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();    
    //
    if (deltay>0 && TMath::Abs(y1max-TMath::Abs(y1))> deltay ) continue;  // seed only at the edge    
    // 
    Int_t  index1 = TMath::Max(kr3.Find(z1-dzmax)-1,0);
    Int_t  index2 = TMath::Min(kr3.Find(z1+dzmax)+1,kr3);
    //    
    Double_t y3,   z3;
    //
    //
    UInt_t index;
    for (Int_t js=index1; js < index2; js++) {
      const AliTPCclusterMI *kcl = kr3[js];
      if (kcl->IsDisabled()) {
	continue;
      }

      if (kcl->IsUsed(10)) continue;
      y3 = kcl->GetY(); 
      // apply angular cuts
      if (TMath::Abs(y1-y3)>dymax) continue;
      //x3 = x3; 
      z3 = kcl->GetZ();	
      if (TMath::Abs(z1-z3)>dzmax) continue;
      //
      Double_t angley = (y1-y3)/(x1-x3);
      Double_t anglez = (z1-z3)/(x1-x3);
      //
      Double_t erry = TMath::Abs(angley)*(x1-x1m)*0.5+0.5;
      Double_t errz = TMath::Abs(anglez)*(x1-x1m)*0.5+0.5;
      //
      Double_t yyym = angley*(xm-x1)+y1;
      Double_t zzzm = anglez*(xm-x1)+z1;

      const AliTPCclusterMI *kcm = krm.FindNearest2(yyym,zzzm,erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      if (!kcm) continue;
      if (kcm->IsUsed(10)) continue;
      if (kcm->IsDisabled()) {
	continue;
      }

      erry = TMath::Abs(angley)*(x1-x1m)*0.4+0.5;
      errz = TMath::Abs(anglez)*(x1-x1m)*0.4+0.5;
      //
      //
      //
      Int_t used  =0;
      Int_t found =0;
      //
      // look around first
      const AliTPCclusterMI *kc1m = kr1m.FindNearest2(angley*(x1m-x1)+y1,
						      anglez*(x1m-x1)+z1,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc1m){
	found++;
	if (kc1m->IsUsed(10)) used++;
      }
      const AliTPCclusterMI *kc1p = kr1p.FindNearest2(angley*(x1p-x1)+y1,
						      anglez*(x1p-x1)+z1,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc1p){
	found++;
	if (kc1p->IsUsed(10)) used++;
      }
      if (used>1)  continue;
      if (found<1) continue; 

      //
      // look around last
      const AliTPCclusterMI *kc3m = kr3m.FindNearest2(angley*(x3m-x3)+y3,
						      anglez*(x3m-x3)+z3,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc3m){
	found++;
	if (kc3m->IsUsed(10)) used++;
      }
      else 
	continue;
      const AliTPCclusterMI *kc3p = kr3p.FindNearest2(angley*(x3p-x3)+y3,
						      anglez*(x3p-x3)+z3,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc3p){
	found++;
	if (kc3p->IsUsed(10)) used++;
      }
      else 
	continue;
      if (used>1)  continue;
      if (found<3) continue;       
      //
      Double_t x2,y2,z2;
      x2 = xm;
      y2 = kcm->GetY();
      z2 = kcm->GetZ();
      //
                  	
      x[0]=y1;
      x[1]=z1;
      x[4]=F1(x1,y1,x2,y2,x3,y3);
      //if (TMath::Abs(x[4]) >= cuts[0]) continue;
      nin0++;
      //
      x[2]=F2(x1,y1,x2,y2,x3,y3);
      nin1++;
      //
      x[3]=F3n(x1,y1,x2,y2,z1,z2,x[4]);
      //if (TMath::Abs(x[3]) > cuts[3]) continue;
      nin2++;
      //
      //
      Double_t sy1=0.1,  sz1=0.1;
      Double_t sy2=0.1,  sz2=0.1;
      Double_t sy3=0.1,  sy=0.1, sz=0.1;
      
      Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
      Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
      Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
      Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
      Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
      Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
      
      Double_t f30=(F3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
      Double_t f31=(F3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
      Double_t f32=(F3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
      Double_t f34=(F3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
      
      c[0]=sy1;
      c[1]=0.;       c[2]=sz1; 
      c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
      c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
      c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
      c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
      c[13]=f30*sy1*f40+f32*sy2*f42;
      c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
      
      //	if (!BuildSeed(kr1[is],kcl,0,x1,x2,x3,x,c)) continue;
      
      index=kr1.GetIndex(is);
      if (seed) {MarkSeedFree( seed ); seed = 0;}
      AliTPCseed *track = seed = new( NextFreeSeed() ) AliTPCseed(x1, sec*alpha+shift, x, c, index);
      seed->SetPoolID(fLastSeedID);
      seed->SetRow(i1-1); //RS: memorise current row      
      track->SetIsSeeding(kTRUE);

      nin++;      
      FollowProlongation(*track, i1-7,1);
      if (track->GetNumberOfClusters() < track->GetNFoundable()*0.75 || 
	  track->GetNShared()>0.6*track->GetNumberOfClusters() || ( track->GetSigmaY2()+ track->GetSigmaZ2())>0.6+fExtraClErrYZ2){
	MarkSeedFree( seed ); seed = 0;
	continue;
      }
      nout1++;
      nout2++;	
      //Int_t rc = 1;
      FollowProlongation(*track, i2,1);
      track->SetBConstrain(0);
      track->SetLastPoint(i1+fInnerSec->GetNRows());  // first cluster in track position
      track->SetFirstPoint(track->GetLastPoint());
      
      if (track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	  track->GetNumberOfClusters()<track->GetNFoundable()*0.7 || 
	  track->GetNShared()>2. || track->GetChi2()/track->GetNumberOfClusters()>6 || ( track->GetSigmaY2()+ track->GetSigmaZ2())>0.5+fExtraClErrYZ2) {
	MarkSeedFree( seed ); seed = 0;
	continue;
      }
   
      {
	FollowProlongation(*track, TMath::Max(i2-10,0),1);
	AliTPCseed * track2 = MakeSeed(track,0.2,0.5,0.9);
	FollowProlongation(*track2, i2,1);
	track2->SetBConstrain(kFALSE);
	track2->SetSeedType(4);
	arr->AddLast(track2);
	MarkSeedFree( seed ); seed = 0;
      }
      
   
      //arr->AddLast(track); 
      //seed = new AliTPCseed; 	
      nout3++;
    }
  }
  
  if (fDebug>3){
    Info("MakeSeeds5","\nSeeding statiistic:\t%d\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2,nout3);
  }
  if (seed) MarkSeedFree(seed);
}

void AliTPCtracker::MakeSeeds5Dist(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2,  Float_t cuts[4],
				 Float_t deltay) {
  


  //-----------------------------------------------------------------
  // This function creates track seeds, accounting for distortions
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut
  // cuts[3]   - fP3 cut


  Int_t nin0  = 0;
  Int_t nin1  = 0;
  Int_t nin2  = 0;
  Int_t nin   = 0;
  Int_t nout1 = 0;
  Int_t nout2 = 0;
  Int_t nout3 =0;
  Double_t x[5], c[15];
  //
  // make temporary seed
  AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed;
  seed->SetPoolID(fLastSeedID);
  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  //
  //

  // first 3 padrows
  Double_t x1Def = GetXrow(i1-1);
  const    AliTPCtrackerRow& kr1=GetRow(sec,i1-1);
  Double_t y1max  = GetMaxY(i1-1)-kr1.GetDeadZone()-1.5;  
  //
  Double_t x1pDef = GetXrow(i1);
  const    AliTPCtrackerRow& kr1p=GetRow(sec,i1);
  //
  Double_t x1mDef = GetXrow(i1-2);
  const    AliTPCtrackerRow& kr1m=GetRow(sec,i1-2);

  double dx11mDef = x1Def-x1mDef;
  double dx11pDef = x1Def-x1pDef;
  //
  //last 3 padrow for seeding
  AliTPCtrackerRow&  kr3  = GetRow((sec+fkNOS)%fkNOS,i1-7);
  Double_t    x3Def   =  GetXrow(i1-7);
  //
  AliTPCtrackerRow&  kr3p  = GetRow((sec+fkNOS)%fkNOS,i1-6);
  Double_t    x3pDef   = GetXrow(i1-6);
  //
  AliTPCtrackerRow&  kr3m  = GetRow((sec+fkNOS)%fkNOS,i1-8);
  Double_t    x3mDef   = GetXrow(i1-8);

  //
  double dx33mDef = x3Def-x3mDef;
  double dx33pDef = x3Def-x3pDef;

  //
  // middle padrow
  Int_t im = i1-4;                           //middle pad row index
  Double_t xmDef         = GetXrow(im);         // radius of middle pad-row
  const AliTPCtrackerRow& krm=GetRow(sec,im);   //middle pad -row
  //
  //
  Double_t deltax  = x1Def-x3Def;
  Double_t dymax   = deltax*cuts[1];
  Double_t dzmax   = deltax*cuts[3];
  //
  // loop over clusters  
  for (Int_t is=0; is < kr1; is++) {
    //
    if (kr1[is]->IsUsed(10)) continue;	
    if (kr1[is]->IsDisabled()) {
      continue;
    }
    const AliTPCclusterMI* clkr1 = kr1[is];
    Double_t x1=clkr1->GetX(), y1=clkr1->GetY(), z1=clkr1->GetZ();    
    //
    double y1EdgeDist =  y1;
    if  (fAccountDistortions) y1EdgeDist -= GetYSectEdgeDist(sec,i1-1,y1,z1);
    if (deltay>0) {
      double margin = (y1>0 ? y1max-y1EdgeDist : y1max + y1EdgeDist);
      if (margin<deltay ) continue;  // seed only at the edge
    }

    // 
    Int_t  index1 = TMath::Max(kr3.Find(z1-dzmax)-1,0);
    Int_t  index2 = TMath::Min(kr3.Find(z1+dzmax)+1,kr3);
    //    
    Double_t x3,y3,z3;
    //
    //
    UInt_t index;
    for (Int_t js=index1; js < index2; js++) {
      const AliTPCclusterMI *kcl = kr3[js];
      if (kcl->IsDisabled()) {
	continue;
      }

      if (kcl->IsUsed(10)) continue;
      y3 = kcl->GetY(); 
      // apply angular cuts
      if (TMath::Abs(y1-y3)>dymax) continue;
      //x3 = x3; 
      z3 = kcl->GetZ();	
      if (TMath::Abs(z1-z3)>dzmax) continue;
      
      x3 = kcl->GetX();
      //
      double dx13 = x1-x3;
      if (TMath::Abs(dx13)<0.1) {
	//AliErrorF("Wrong X correction? Sec%d : row%d@X=%.2f->%.2f (z=%.2f) row%d@X=%.2f->%.2f (z=%.2f)\n",sec,
	//	  i1-1,x1Def,x1,z1, i1-7,x3Def,x3,z3);
	continue; // distortions should not make distance so small
      }

      Double_t angley = (y1-y3)/dx13;
      Double_t anglez = (z1-z3)/dx13;
      //
      Double_t erry = TMath::Abs(angley)*dx11mDef*0.5+0.5; // RS: use ideal X differences assuming
      Double_t errz = TMath::Abs(anglez)*dx11mDef*0.5+0.5; // that the distortions are ~same on close rows
      //
      Double_t yyym = angley*(xmDef-x1Def)+y1; // RS: idem, assume distortions cancels in the difference
      Double_t zzzm = anglez*(xmDef-x1Def)+z1;

      const AliTPCclusterMI *kcm = krm.FindNearest2(yyym,zzzm,erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      if (!kcm) continue;
      if (kcm->IsUsed(10)) continue;
      if (kcm->IsDisabled()) {
	continue;
      }

      erry = TMath::Abs(angley)*dx11mDef*0.4+0.5;
      errz = TMath::Abs(anglez)*dx11mDef*0.4+0.5;
      //
      //
      //
      Int_t used  =0;
      Int_t found =0;
      //
      // look around first
      const AliTPCclusterMI *kc1m = kr1m.FindNearest2(-angley*dx11mDef+y1, // RS: idem, assume distortions cancels in the difference
						      -anglez*dx11mDef+z1,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc1m){
	found++;
	if (kc1m->IsUsed(10)) used++;
      }
      const AliTPCclusterMI *kc1p = kr1p.FindNearest2(-angley*dx11pDef+y1, // RS: idem, assume distortions cancels in the difference
						      -anglez*dx11pDef+z1,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc1p){
	found++;
	if (kc1p->IsUsed(10)) used++;
      }
      if (used>1)  continue;
      if (found<1) continue; 

      //
      // look around last
      const AliTPCclusterMI *kc3m = kr3m.FindNearest2(-angley*dx33mDef+y3, // RS: idem, assume distortions cancels in the difference
						      -anglez*dx33mDef+z3,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc3m){
	found++;
	if (kc3m->IsUsed(10)) used++;
      }
      else 
	continue;
      const AliTPCclusterMI *kc3p = kr3p.FindNearest2(-angley*dx33pDef+y3, // RS: idem, assume distortions cancels in the difference
						      -anglez*dx33pDef+z3,
						      erry+fClExtraRoadY,errz+fClExtraRoadZ,index);
      //
      if (kc3p){
	found++;
	if (kc3p->IsUsed(10)) used++;
      }
      else 
	continue;
      if (used>1)  continue;
      if (found<3) continue;       
      //
      Double_t x2,y2,z2;
      x2 = kcm->GetX();
      y2 = kcm->GetY();
      z2 = kcm->GetZ();
      //
                  	
      x[0]=y1;
      x[1]=z1;
      x[4]=F1(x1,y1,x2,y2,x3,y3);
      //if (TMath::Abs(x[4]) >= cuts[0]) continue;
      nin0++;
      //
      x[2]=F2(x1,y1,x2,y2,x3,y3);
      nin1++;
      //
      x[3]=F3n(x1,y1,x2,y2,z1,z2,x[4]);
      //if (TMath::Abs(x[3]) > cuts[3]) continue;
      nin2++;
      //
      //
      Double_t sy1=0.1,  sz1=0.1;
      Double_t sy2=0.1,  sz2=0.1;
      Double_t sy3=0.1,  sy=0.1, sz=0.1;
      
      Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
      Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
      Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
      Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
      Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
      Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
      
      Double_t f30=(F3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
      Double_t f31=(F3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
      Double_t f32=(F3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
      Double_t f34=(F3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;
      
      c[0]=sy1;
      c[1]=0.;       c[2]=sz1; 
      c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
      c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
      c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
      c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
      c[13]=f30*sy1*f40+f32*sy2*f42;
      c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
      
      //	if (!BuildSeed(kr1[is],kcl,0,x1,x2,x3,x,c)) continue;
      
      index=kr1.GetIndex(is);
      if (seed) {MarkSeedFree( seed ); seed = 0;}
      AliTPCseed *track = seed = new( NextFreeSeed() ) AliTPCseed(x1, sec*alpha+shift, x, c, index);
      seed->SetPoolID(fLastSeedID);
      seed->SetRow(i1-1); // RS: memorise current row      
      track->SetIsSeeding(kTRUE);

      nin++;      
      FollowProlongation(*track, i1-7,1);
      if (track->GetNumberOfClusters() < track->GetNFoundable()*0.75 || 
	  track->GetNShared()>0.6*track->GetNumberOfClusters()) {
	//RS: with distortions related cluster errors the track error may grow, don't use this cut
	//|| ( track->GetSigmaY2()+ track->GetSigmaZ2())>0.6){
	MarkSeedFree( seed ); seed = 0;
	continue;
      }
      nout1++;
      nout2++;	
      //Int_t rc = 1;
      FollowProlongation(*track, i2,1);
      track->SetBConstrain(0);
      track->SetLastPoint(i1+fInnerSec->GetNRows());  // first cluster in track position
      track->SetFirstPoint(track->GetLastPoint());
      
      if (track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	  track->GetNumberOfClusters()<track->GetNFoundable()*0.7 || 
	  track->GetNShared()>2. || track->GetChi2()/track->GetNumberOfClusters()>6) {
	//RS: with distortions related cluster errors the track error may grow, don't use this cut
	//|| ( track->GetSigmaY2()+ track->GetSigmaZ2())>0.5 ) {
	MarkSeedFree( seed ); seed = 0;
	continue;
      }
   
      {
	FollowProlongation(*track, TMath::Max(i2-10,0),1);
	AliTPCseed * track2 = MakeSeed(track,0.2,0.5,0.9);
	FollowProlongation(*track2, i2,1);
	track2->SetBConstrain(kFALSE);
	track2->SetSeedType(4);
	arr->AddLast(track2);
	MarkSeedFree( seed ); seed = 0;
      }
      
   
      //arr->AddLast(track); 
      //seed = new AliTPCseed; 	
      nout3++;
    }
  }
  
  if (fDebug>3){
    Info("MakeSeeds5Dist","\nSeeding statiistic:\t%d\t%d\t%d\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin,nout1,nout2,nout3);
  }
  if (seed) MarkSeedFree(seed);
}


//_____________________________________________________________________________
void AliTPCtracker::MakeSeeds2(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t */*cuts[4]*/,
				 Float_t deltay, Bool_t /*bconstrain*/) {
  //-----------------------------------------------------------------
  // This function creates track seeds - without vertex constraint
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut        - not applied
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut    - not applied 
  // cuts[3]   - fP3 cut
  const double kRoadZ = 1.2, kRoadY = 1.2;

  Int_t nin0=0;
  Int_t nin1=0;
  Int_t nin2=0;
  Int_t nin3=0;
  //  Int_t nin4=0;
  //Int_t nin5=0;

  

  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  Int_t row0 = (i1+i2)/2;
  Int_t drow = (i1-i2)/2;
  const AliTPCtrackerRow& kr0=fSectors[sec][row0];
  AliTPCtrackerRow * kr=0;

  AliTPCpolyTrack polytrack;
  Int_t nclusters=fSectors[sec][row0];
  AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed;
  seed->SetPoolID(fLastSeedID);

  Int_t sumused=0;
  Int_t cused=0;
  Int_t cnused=0;
  for (Int_t is=0; is < nclusters; is++) {  //LOOP over clusters
    Int_t nfound =0;
    Int_t nfoundable =0;
    for (Int_t iter =1; iter<2; iter++){   //iterations
      const AliTPCtrackerRow& krm=fSectors[sec][row0-iter];
      const AliTPCtrackerRow& krp=fSectors[sec][row0+iter];      
      const AliTPCclusterMI * cl= kr0[is];
      if (cl->IsDisabled()) {
	continue;
      }

      if (cl->IsUsed(10)) {
	cused++;
      }
      else{
	cnused++;
      }
      Double_t x = kr0.GetX();
      // Initialization of the polytrack
      nfound =0;
      nfoundable =0;
      polytrack.Reset();
      //
      Double_t y0= cl->GetY();
      Double_t z0= cl->GetZ();
      Float_t erry = 0;
      Float_t errz = 0;
      
      Double_t ymax = fSectors->GetMaxY(row0)-kr0.GetDeadZone()-1.5;
      if (deltay>0 && TMath::Abs(ymax-TMath::Abs(y0))> deltay ) continue;  // seed only at the edge
      
      erry = (0.5)*cl->GetSigmaY2()/TMath::Sqrt(cl->GetQ())*6;	    
      errz = (0.5)*cl->GetSigmaZ2()/TMath::Sqrt(cl->GetQ())*6;      
      polytrack.AddPoint(x,y0,z0,erry, errz);

      sumused=0;
      if (cl->IsUsed(10)) sumused++;


      Float_t roady = (5*TMath::Sqrt(cl->GetSigmaY2()+0.2)+1.)*iter;
      Float_t roadz = (5*TMath::Sqrt(cl->GetSigmaZ2()+0.2)+1.)*iter;
      //
      x = krm.GetX();
      AliTPCclusterMI * cl1 = krm.FindNearest(y0,z0,roady+fClExtraRoadY,roadz+fClExtraRoadZ);
      if (cl1 && TMath::Abs(ymax-TMath::Abs(y0))) {
	erry = (0.5)*cl1->GetSigmaY2()/TMath::Sqrt(cl1->GetQ())*3;	    
	errz = (0.5)*cl1->GetSigmaZ2()/TMath::Sqrt(cl1->GetQ())*3;
	if (cl1->IsUsed(10))  sumused++;
	polytrack.AddPoint(x,cl1->GetY(),cl1->GetZ(),erry,errz);
      }
      //
      x = krp.GetX();
      AliTPCclusterMI * cl2 = krp.FindNearest(y0,z0,roady+fClExtraRoadY,roadz+fClExtraRoadZ);
      if (cl2) {
	erry = (0.5)*cl2->GetSigmaY2()/TMath::Sqrt(cl2->GetQ())*3;	    
	errz = (0.5)*cl2->GetSigmaZ2()/TMath::Sqrt(cl2->GetQ())*3;
	if (cl2->IsUsed(10)) sumused++;	 
	polytrack.AddPoint(x,cl2->GetY(),cl2->GetZ(),erry,errz);
      }
      //
      if (sumused>0) continue;
      nin0++;
      polytrack.UpdateParameters();
      // follow polytrack
      //
      Double_t yn,zn;
      nfoundable = polytrack.GetN();
      nfound     = nfoundable; 
      //
      for (Int_t ddrow = iter+1; ddrow<drow;ddrow++){
	Float_t maxdist = 0.8*(1.+3./(ddrow));
	for (Int_t delta = -1;delta<=1;delta+=2){
	  Int_t row = row0+ddrow*delta;
	  kr = &(fSectors[sec][row]);
	  Double_t xn = kr->GetX();
	  Double_t ymax1 = fSectors->GetMaxY(row)-kr->GetDeadZone()-1.5;
	  polytrack.GetFitPoint(xn,yn,zn);
	  if (TMath::Abs(yn)>ymax1) continue;
	  nfoundable++;
	  AliTPCclusterMI * cln = kr->FindNearest(yn,zn,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ);
	  if (cln) {
	    Float_t dist =  TMath::Sqrt(  (yn-cln->GetY())*(yn-cln->GetY())+(zn-cln->GetZ())*(zn-cln->GetZ()));
	    if (dist<maxdist){
	      /*
	      erry = (dist+0.3)*cln->GetSigmaY2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddrow));	    
	      errz = (dist+0.3)*cln->GetSigmaZ2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddrow));
	      if (cln->IsUsed(10)) {
		//	printf("used\n");
		sumused++;
		erry*=2;
		errz*=2;
	      }
	      */
	      erry=0.1;
	      errz=0.1;
	      polytrack.AddPoint(xn,cln->GetY(),cln->GetZ(),erry, errz);
	      nfound++;
	    }
	  }
	}
	if ( (sumused>3) || (sumused>0.5*nfound) || (nfound<0.6*nfoundable))  break;     
	polytrack.UpdateParameters();
      }           
    }
    if ( (sumused>3) || (sumused>0.5*nfound))  {
      //printf("sumused   %d\n",sumused);
      continue;
    }
    nin1++;
    Double_t dy,dz;
    polytrack.GetFitDerivation(kr0.GetX(),dy,dz);
    AliTPCpolyTrack track2;
    
    polytrack.Refit(track2,0.5+TMath::Abs(dy)*0.3,0.4+TMath::Abs(dz)*0.3);
    if (track2.GetN()<0.5*nfoundable) continue;
    nin2++;

    if ((nfound>0.6*nfoundable) &&( nfoundable>0.4*(i1-i2))) {
      //
      // test seed with and without constrain
      for (Int_t constrain=0; constrain<=0;constrain++){
	// add polytrack candidate

	Double_t x[5], c[15];
	Double_t x1,x2,x3,y1,y2,y3,z1,z2,z3;
	track2.GetBoundaries(x3,x1);	
	x2 = (x1+x3)/2.;
	track2.GetFitPoint(x1,y1,z1);
	track2.GetFitPoint(x2,y2,z2);
	track2.GetFitPoint(x3,y3,z3);
	//
	//is track pointing to the vertex ?
	Double_t x0,y0,z0;
	x0=0;
	polytrack.GetFitPoint(x0,y0,z0);

	if (constrain) {
	  x2 = x3;
	  y2 = y3;
	  z2 = z3;
	  
	  x3 = 0;
	  y3 = 0;
	  z3 = 0;
	}
	x[0]=y1;
	x[1]=z1;
	x[4]=F1(x1,y1,x2,y2,x3,y3);
		
	//	if (TMath::Abs(x[4]) >= cuts[0]) continue;  //
	x[2]=F2(x1,y1,x2,y2,x3,y3);
	
	//if (TMath::Abs(x[4]*x1-x[2]) >= cuts[1]) continue;
	//x[3]=F3(x1,y1,x2,y2,z1,z2);
	x[3]=F3n(x1,y1,x3,y3,z1,z3,x[4]);
	//if (TMath::Abs(x[3]) > cuts[3]) continue;

	
	Double_t sy =0.1, sz =0.1;
	Double_t sy1=0.02, sz1=0.02;
	Double_t sy2=0.02, sz2=0.02;
	Double_t sy3=0.02;

	if (constrain){
	  sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	}
	
	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;

	Double_t f30=(F3(x1,y1+sy,x3,y3,z1,z3)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x3,y3,z1+sz,z3)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x3,y3+sy,z1,z3)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x3,y3,z1,z3+sz)-x[3])/sz;

	
	c[0]=sy1;
	c[1]=0.;       c[2]=sz1;
	c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
	c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
	c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
	c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
	c[13]=f30*sy1*f40+f32*sy2*f42;
	c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	//Int_t row1 = fSectors->GetRowNumber(x1);
	Int_t row1 = GetRowNumber(x1);

	UInt_t index=0;
	//kr0.GetIndex(is);
	if (seed) {MarkSeedFree( seed ); seed = 0;}
	AliTPCseed *track = seed = new( NextFreeSeed() ) AliTPCseed(x1,sec*alpha+shift,x,c,index);
	seed->SetPoolID(fLastSeedID);
	track->SetIsSeeding(kTRUE);
	Int_t rc=FollowProlongation(*track, i2);	
	if (constrain) track->SetBConstrain(1);
	else
	  track->SetBConstrain(0);
	track->SetLastPoint(row1+fInnerSec->GetNRows());  // first cluster in track position
	track->SetFirstPoint(track->GetLastPoint());

	if (rc==0 || track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	    track->GetNumberOfClusters() < track->GetNFoundable()*0.6 || 
	    track->GetNShared()>0.4*track->GetNumberOfClusters()) {
	  MarkSeedFree( seed ); seed = 0;
	}
	else {
	  arr->AddLast(track); // track IS seed, don't free seed
	  seed = new( NextFreeSeed() ) AliTPCseed;
	  seed->SetPoolID(fLastSeedID);
	}
	nin3++;
      }
    }  // if accepted seed
  }
  if (fDebug>3){
    Info("MakeSeeds2","\nSeeding statiistic:\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin3);
  }
  if (seed) MarkSeedFree( seed );
}

//_____________________________________________________________________________
void AliTPCtracker::MakeSeeds2Dist(TObjArray * arr, Int_t sec, Int_t i1, Int_t i2, Float_t */*cuts[4]*/,
				 Float_t deltay, Bool_t /*bconstrain*/) {
  //-----------------------------------------------------------------
  // This function creates track seeds, accounting for distortions - without vertex constraint
  //-----------------------------------------------------------------
  // cuts[0]   - fP4 cut        - not applied
  // cuts[1]   - tan(phi)  cut
  // cuts[2]   - zvertex cut    - not applied 
  // cuts[3]   - fP3 cut
  const double kRoadZ = 1.2, kRoadY = 1.2;

  Int_t nin0=0;
  Int_t nin1=0;
  Int_t nin2=0;
  Int_t nin3=0;
  //  Int_t nin4=0;
  //Int_t nin5=0;

  AliFatal("This method is still not fully aware of distortions, should not be used");
  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  //  Double_t cs=cos(alpha), sn=sin(alpha);
  Int_t row0 = (i1+i2)/2;
  Int_t drow = (i1-i2)/2;
  const AliTPCtrackerRow& kr0=fSectors[sec][row0];
  AliTPCtrackerRow * kr=0;

  AliTPCpolyTrack polytrack;
  Int_t nclusters=fSectors[sec][row0];
  AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed;
  seed->SetPoolID(fLastSeedID);

  Int_t sumused=0;
  Int_t cused=0;
  Int_t cnused=0;
  Int_t rowMax = -1;
  //
  for (Int_t is=0; is < nclusters; is++) {  //LOOP over clusters
    Int_t nfound =0;
    Int_t nfoundable =0;
    for (Int_t iter =1; iter<2; iter++){   //iterations
      const AliTPCtrackerRow& krm=fSectors[sec][row0-iter];
      const AliTPCtrackerRow& krp=fSectors[sec][row0+iter];      
      const AliTPCclusterMI * cl= kr0[is];
      if (cl->IsDisabled()) {
	continue;
      }

      if (cl->IsUsed(10)) {
	cused++;
      }
      else{
	cnused++;
      }
      
      Double_t x = cl->GetX(), xDist = x - kr0.GetX(); // approximate X distortion in proximity of cl
      // Initialization of the polytrack
      nfound =0;
      nfoundable =0;
      polytrack.Reset();
      //
      Double_t y0= cl->GetY();
      Double_t z0= cl->GetZ();
      Float_t erry = 0;
      Float_t errz = 0;
      
      Double_t ymax = fSectors->GetMaxY(row0)-kr0.GetDeadZone()-1.5; // RS: watch dead zones
      if (deltay>0 && TMath::Abs(ymax-TMath::Abs(y0))> deltay ) continue;  // seed only at the edge
      
      erry = (0.5)*cl->GetSigmaY2()/TMath::Sqrt(cl->GetQ())*6;	    
      errz = (0.5)*cl->GetSigmaZ2()/TMath::Sqrt(cl->GetQ())*6;      
      polytrack.AddPoint(x,y0,z0,erry, errz);
      rowMax = row0;      
      sumused=0;
      if (cl->IsUsed(10)) sumused++;


      Float_t roady = (5*TMath::Sqrt(cl->GetSigmaY2()+0.2)+1.)*iter;
      Float_t roadz = (5*TMath::Sqrt(cl->GetSigmaZ2()+0.2)+1.)*iter;
      //
      x = krm.GetX() + xDist; // RS: assume distortion at krm is similar to that at kr
      AliTPCclusterMI * cl1 = krm.FindNearest(y0,z0,roady+fClExtraRoadY,roadz+fClExtraRoadZ);
      if (cl1 && TMath::Abs(ymax-TMath::Abs(y0))) {
	erry = (0.5)*cl1->GetSigmaY2()/TMath::Sqrt(cl1->GetQ())*3;	    
	errz = (0.5)*cl1->GetSigmaZ2()/TMath::Sqrt(cl1->GetQ())*3;
	if (cl1->IsUsed(10))  sumused++;
	//RS: use real cluster X instead of approximately distorted
	polytrack.AddPoint(cl1->GetX(),cl1->GetY(),cl1->GetZ(),erry,errz); 
      }
      //
      x = krp.GetX() + xDist; // RS: assume distortion at krp is similar to that at kr
      AliTPCclusterMI * cl2 = krp.FindNearest(y0,z0,roady+fClExtraRoadY,roadz+fClExtraRoadZ);
      if (cl2) {
	erry = (0.5)*cl2->GetSigmaY2()/TMath::Sqrt(cl2->GetQ())*3;	    
	errz = (0.5)*cl2->GetSigmaZ2()/TMath::Sqrt(cl2->GetQ())*3;
	if (cl2->IsUsed(10)) sumused++;	 
	//RS: use real cluster X instead of approximately distorted
	polytrack.AddPoint(cl2->GetX(),cl2->GetY(),cl2->GetZ(),erry,errz); 
	rowMax = row0+iter;
      }
      //
      if (sumused>0) continue;
      nin0++;
      polytrack.UpdateParameters();
      // follow polytrack
      //
      Double_t yn,zn;
      nfoundable = polytrack.GetN();
      nfound     = nfoundable; 
      //
      for (Int_t ddrow = iter+1; ddrow<drow;ddrow++){
	Float_t maxdist = 0.8*(1.+3./(ddrow));
	for (Int_t delta = -1;delta<=1;delta+=2){
	  Int_t row = row0+ddrow*delta;
	  kr = &(fSectors[sec][row]);
	  Double_t xn = kr->GetX(); //RS: use row X as a first guess about distorted cluster x
	  Double_t ymax1 = fSectors->GetMaxY(row)-kr->GetDeadZone()-1.5; // RS: watch dead zones
	  polytrack.GetFitPoint(xn,yn,zn);
	  double dxn = GetDistortionX(xn, yn, zn, sec, row);
	  if (TMath::Abs(dxn)>0.1) {    //RS: account for distortion
	    xn += dxn;	
	    polytrack.GetFitPoint(xn,yn,zn);    
	  }
	  //
	  if (TMath::Abs(yn)>ymax1) continue; // RS:? watch dead zones
	  nfoundable++;
	  AliTPCclusterMI * cln = kr->FindNearest(yn,zn,kRoadY+fClExtraRoadY,kRoadZ+fClExtraRoadZ);
	  if (cln) {
	    Float_t dist =  TMath::Sqrt(  (yn-cln->GetY())*(yn-cln->GetY())+(zn-cln->GetZ())*(zn-cln->GetZ()));
	    if (dist<maxdist){
	      /*
	      erry = (dist+0.3)*cln->GetSigmaY2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddrow));	    
	      errz = (dist+0.3)*cln->GetSigmaZ2()/TMath::Sqrt(cln->GetQ())*(1.+1./(ddrow));
	      if (cln->IsUsed(10)) {
		//	printf("used\n");
		sumused++;
		erry*=2;
		errz*=2;
	      }
	      */
	      erry=0.1;
	      errz=0.1;
	      polytrack.AddPoint(xn,cln->GetY(),cln->GetZ(),erry, errz);
	      if (row>rowMax) rowMax = row;
	      nfound++;
	    }
	  }
	}
	if ( (sumused>3) || (sumused>0.5*nfound) || (nfound<0.6*nfoundable))  break;     
	polytrack.UpdateParameters();
      }           
    }
    if ( (sumused>3) || (sumused>0.5*nfound))  {
      //printf("sumused   %d\n",sumused);
      continue;
    }
    nin1++;
    Double_t dy,dz;
    polytrack.GetFitDerivation(kr0.GetX(),dy,dz); // RS: Note: derivative is at ideal row X
    AliTPCpolyTrack track2;
    
    polytrack.Refit(track2,0.5+TMath::Abs(dy)*0.3,0.4+TMath::Abs(dz)*0.3);
    if (track2.GetN()<0.5*nfoundable) continue;
    nin2++;

    if ((nfound>0.6*nfoundable) &&( nfoundable>0.4*(i1-i2))) {
      //
      // test seed with and without constrain
      for (Int_t constrain=0; constrain<=0;constrain++){
	// add polytrack candidate

	Double_t x[5], c[15];
	Double_t x1,x2,x3,y1,y2,y3,z1,z2,z3;
	track2.GetBoundaries(x3,x1);	
	x2 = (x1+x3)/2.;
	track2.GetFitPoint(x1,y1,z1);
	track2.GetFitPoint(x2,y2,z2);
	track2.GetFitPoint(x3,y3,z3);
	//
	//is track pointing to the vertex ?
	Double_t x0,y0,z0;
	x0=0;
	polytrack.GetFitPoint(x0,y0,z0);

	if (constrain) {
	  x2 = x3;
	  y2 = y3;
	  z2 = z3;
	  
	  x3 = 0;
	  y3 = 0;
	  z3 = 0;
	}
	x[0]=y1;
	x[1]=z1;
	x[4]=F1(x1,y1,x2,y2,x3,y3);
		
	//	if (TMath::Abs(x[4]) >= cuts[0]) continue;  //
	x[2]=F2(x1,y1,x2,y2,x3,y3);
	
	//if (TMath::Abs(x[4]*x1-x[2]) >= cuts[1]) continue;
	//x[3]=F3(x1,y1,x2,y2,z1,z2);
	x[3]=F3n(x1,y1,x3,y3,z1,z3,x[4]);
	//if (TMath::Abs(x[3]) > cuts[3]) continue;

	
	Double_t sy =0.1, sz =0.1;
	Double_t sy1=0.02, sz1=0.02;
	Double_t sy2=0.02, sz2=0.02;
	Double_t sy3=0.02;

	if (constrain){
	  sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;
	}
	
	Double_t f40=(F1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
	Double_t f42=(F1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
	Double_t f43=(F1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
	Double_t f20=(F2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(F2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f23=(F2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;

	Double_t f30=(F3(x1,y1+sy,x3,y3,z1,z3)-x[3])/sy;
	Double_t f31=(F3(x1,y1,x3,y3,z1+sz,z3)-x[3])/sz;
	Double_t f32=(F3(x1,y1,x3,y3+sy,z1,z3)-x[3])/sy;
	Double_t f34=(F3(x1,y1,x3,y3,z1,z3+sz)-x[3])/sz;

	
	c[0]=sy1;
	c[1]=0.;       c[2]=sz1;
	c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
	c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
	c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
	c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
	c[13]=f30*sy1*f40+f32*sy2*f42;
	c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
	
	//Int_t row1 = GetRowNumber(x1); // RS: this is now substituted by rowMax

	UInt_t index=0;
	//kr0.GetIndex(is);
	if (seed) {MarkSeedFree( seed ); seed = 0;}
	AliTPCseed *track = seed = new( NextFreeSeed() ) AliTPCseed(x1,sec*alpha+shift,x,c,index);
	seed->SetPoolID(fLastSeedID);
	seed->SetRow(rowMax); //RS: memorise row of x1
	track->SetIsSeeding(kTRUE);
	Int_t rc=FollowProlongation(*track, i2);	
	if (constrain) track->SetBConstrain(1);
	else
	  track->SetBConstrain(0);
	track->SetLastPoint(rowMax); //row1+fInnerSec->GetNRows());  // first cluster in track position
	track->SetFirstPoint(track->GetLastPoint());

	if (rc==0 || track->GetNumberOfClusters()<(i1-i2)*0.5 || 
	    track->GetNumberOfClusters() < track->GetNFoundable()*0.6 || 
	    track->GetNShared()>0.4*track->GetNumberOfClusters()) {
	  MarkSeedFree( seed ); seed = 0;
	}
	else {
	  arr->AddLast(track); // track IS seed, don't free seed
	  seed = new( NextFreeSeed() ) AliTPCseed;
	  seed->SetPoolID(fLastSeedID);
	}
	nin3++;
      }
    }  // if accepted seed
  }
  if (fDebug>3){
    Info("MakeSeeds2","\nSeeding statiistic:\t%d\t%d\t%d\t%d",nin0,nin1,nin2,nin3);
  }
  if (seed) MarkSeedFree( seed );
}


AliTPCseed *AliTPCtracker::MakeSeed(AliTPCseed *const track, Float_t r0, Float_t r1, Float_t r2)
{
  //
  //
  //reseed using track points
  Int_t p0 = int(r0*track->GetNumberOfClusters());     // point 0 
  Int_t p1 = int(r1*track->GetNumberOfClusters());
  Int_t p2 = int(r2*track->GetNumberOfClusters());   // last point
  Int_t pp2=0;
  Double_t  x0[3],x1[3],x2[3];
  for (Int_t i=0;i<3;i++){
    x0[i]=-1;
    x1[i]=-1;
    x2[i]=-1;
  }

  // find track position at given ratio of the length
  Int_t  sec0=0, sec1=0, sec2=0;
  Int_t index=-1;
  Int_t clindex;
  for (Int_t i=0;i<kMaxRow;i++){
    if (track->GetClusterIndex2(i)>=0){
      index++;
      const AliTPCTrackerPoints::Point *trpoint =track->GetTrackPoint(i);
      if ( (index<p0) || x0[0]<0 ){
	if (trpoint->GetX()>1){
	  clindex = track->GetClusterIndex2(i);
	  if (clindex >= 0){	
	    x0[0] = trpoint->GetX();
	    x0[1] = trpoint->GetY();
	    x0[2] = trpoint->GetZ();
	    sec0  = ((clindex&0xff000000)>>24)%18;
	  }
	}
      }

      if ( (index<p1) &&(trpoint->GetX()>1)){
	clindex = track->GetClusterIndex2(i);
	if (clindex >= 0){
	  x1[0] = trpoint->GetX();
	  x1[1] = trpoint->GetY();
	  x1[2] = trpoint->GetZ();
	  sec1  = ((clindex&0xff000000)>>24)%18;
	}
      }
      if ( (index<p2) &&(trpoint->GetX()>1)){
	clindex = track->GetClusterIndex2(i);
	if (clindex >= 0){
	  x2[0] = trpoint->GetX();
	  x2[1] = trpoint->GetY();
	  x2[2] = trpoint->GetZ(); 
	  sec2  = ((clindex&0xff000000)>>24)%18;
	  pp2 = i;
	}
      }
    }
  }
  
  Double_t alpha, cs,sn, xx2,yy2;
  //
  alpha = (sec1-sec2)*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= x1[0]*cs-x1[1]*sn;
  yy2= x1[0]*sn+x1[1]*cs;
  x1[0] = xx2;
  x1[1] = yy2;
  //
  alpha = (sec0-sec2)*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= x0[0]*cs-x0[1]*sn;
  yy2= x0[0]*sn+x0[1]*cs;
  x0[0] = xx2;
  x0[1] = yy2;
  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=x2[1];
  x[1]=x2[2];
  x[4]=F1(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]);
  //  if (x[4]>1) return 0;
  x[2]=F2(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]);
  x[3]=F3n(x2[0],x2[1],x0[0],x0[1],x2[2],x0[2],x[4]);
  //if (TMath::Abs(x[3]) > 2.2)  return 0;
  //if (TMath::Abs(x[2]) > 1.99) return 0;
  //  
  Double_t sy =0.1,  sz =0.1;
  //
  Double_t sy1=0.02+track->GetSigmaY2(), sz1=0.02+track->GetSigmaZ2();
  Double_t sy2=0.01+track->GetSigmaY2(), sz2=0.01+track->GetSigmaZ2();
  Double_t sy3=0.01+track->GetSigmaY2();
  //
  Double_t f40=(F1(x2[0],x2[1]+sy,x1[0],x1[1],x0[0],x0[1])-x[4])/sy;
  Double_t f42=(F1(x2[0],x2[1],x1[0],x1[1]+sy,x0[0],x0[1])-x[4])/sy;
  Double_t f43=(F1(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]+sy)-x[4])/sy;
  Double_t f20=(F2(x2[0],x2[1]+sy,x1[0],x1[1],x0[0],x0[1])-x[2])/sy;
  Double_t f22=(F2(x2[0],x2[1],x1[0],x1[1]+sy,x0[0],x0[1])-x[2])/sy;
  Double_t f23=(F2(x2[0],x2[1],x1[0],x1[1],x0[0],x0[1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(x2[0],x2[1]+sy,x0[0],x0[1],x2[2],x0[2])-x[3])/sy;
  Double_t f31=(F3(x2[0],x2[1],x0[0],x0[1],x2[2]+sz,x0[2])-x[3])/sz;
  Double_t f32=(F3(x2[0],x2[1],x0[0],x0[1]+sy,x2[2],x0[2])-x[3])/sy;
  Double_t f34=(F3(x2[0],x2[1],x0[0],x0[1],x2[2],x0[2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  //  Int_t row1 = fSectors->GetRowNumber(x2[0]);
  AliTPCseed *seed = new( NextFreeSeed() )  AliTPCseed(x2[0], sec2*fSectors->GetAlpha()+fSectors->GetAlphaShift(), x, c, 0);
  seed->SetPoolID(fLastSeedID);
  seed->SetRow(pp2); //RS: memorise current row
  //  Double_t y0,z0,y1,z1, y2,z2;
  //seed->GetProlongation(x0[0],y0,z0);
  // seed->GetProlongation(x1[0],y1,z1);
  //seed->GetProlongation(x2[0],y2,z2);
  //  seed =0;
  seed->SetLastPoint(pp2);
  seed->SetFirstPoint(pp2);
  

  return seed;
}


AliTPCseed *AliTPCtracker::ReSeed(const AliTPCseed *track, Float_t r0, Float_t r1, Float_t r2)
{
  //
  //
  //reseed using founded clusters 
  //
  // Find the number of clusters
  Int_t nclusters = 0;
  for (Int_t irow=0;irow<kMaxRow;irow++){
    if (track->GetClusterIndex(irow)>0) nclusters++;
  }
  //
  Int_t ipos[3];
  ipos[0] = TMath::Max(int(r0*nclusters),0);             // point 0 cluster
  ipos[1] = TMath::Min(int(r1*nclusters),nclusters-1);   // 
  ipos[2] = TMath::Min(int(r2*nclusters),nclusters-1);   // last point
  //
  //
  Double_t  xyz[3][3]={{0}};
  Int_t     row[3]={0},sec[3]={0,0,0};
  //
  // find track row position at given ratio of the length
  Int_t index=-1;
  for (Int_t irow=0;irow<kMaxRow;irow++){    
    if (track->GetClusterIndex2(irow)<0) continue;
    index++;
    for (Int_t ipoint=0;ipoint<3;ipoint++){
      if (index<=ipos[ipoint]) row[ipoint] = irow;
    }        
  }
  //
  //Get cluster and sector position
  for (Int_t ipoint=0;ipoint<3;ipoint++){
    Int_t clindex = track->GetClusterIndex2(row[ipoint]);    
    AliTPCclusterMI * cl = clindex<0 ? 0:GetClusterMI(clindex);
    if (cl==0) {
      //Error("Bug\n");
      //      AliTPCclusterMI * cl = GetClusterMI(clindex);
      return 0;
    }
    sec[ipoint]     = ((clindex&0xff000000)>>24)%18;
    xyz[ipoint][0]  = GetXrow(row[ipoint]);
    xyz[ipoint][1]  = cl->GetY();
    xyz[ipoint][2]  = cl->GetZ();
  }
  //
  //
  // Calculate seed state vector and covariance matrix

  Double_t alpha, cs,sn, xx2,yy2;
  //
  alpha = (sec[1]-sec[2])*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= xyz[1][0]*cs-xyz[1][1]*sn;
  yy2= xyz[1][0]*sn+xyz[1][1]*cs;
  xyz[1][0] = xx2;
  xyz[1][1] = yy2;
  //
  alpha = (sec[0]-sec[2])*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= xyz[0][0]*cs-xyz[0][1]*sn;
  yy2= xyz[0][0]*sn+xyz[0][1]*cs;
  xyz[0][0] = xx2;
  xyz[0][1] = yy2;
  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=xyz[2][1];
  x[1]=xyz[2][2];
  x[4]=F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[2]=F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[3]=F3n(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2],x[4]);
  //  
  Double_t sy =0.1,  sz =0.1;
  //
  Double_t sy1=0.2, sz1=0.2;
  Double_t sy2=0.2, sz2=0.2;
  Double_t sy3=0.2;
  //
  Double_t f40=(F1(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f42=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f43=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[4])/sy;
  Double_t f20=(F2(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f22=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f23=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(xyz[2][0],xyz[2][1]+sy,xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f31=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2]+sz,xyz[0][2])-x[3])/sz;
  Double_t f32=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1]+sy,xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f34=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  //  Int_t row1 = fSectors->GetRowNumber(xyz[2][0]);
  AliTPCseed *seed=new( NextFreeSeed() ) AliTPCseed(xyz[2][0], sec[2]*fSectors->GetAlpha()+fSectors->GetAlphaShift(), x, c, 0);
  seed->SetPoolID(fLastSeedID);
  seed->SetRow(row[2]); //RS: memorise row
  seed->SetLastPoint(row[2]);
  seed->SetFirstPoint(row[2]);  
  return seed;
}


AliTPCseed *AliTPCtracker::ReSeed(AliTPCseed *track,Int_t r0, Bool_t forward)
{
  //
  //
  //reseed using founded clusters 
  //
  Double_t  xyz[3][3];
  Int_t     row[3]={0,0,0};
  Int_t     sec[3]={0,0,0};
  //
  // forward direction
  if (forward){
    for (Int_t irow=r0;irow<kMaxRow;irow++){
      if (track->GetClusterIndex(irow)>0){
	row[0] = irow;
	break;
      }
    }
    for (Int_t irow=kMaxRow;irow>r0;irow--){
      if (track->GetClusterIndex(irow)>0){
	row[2] = irow;
	break;
      }
    }
    for (Int_t irow=row[2]-15;irow>row[0];irow--){
      if (track->GetClusterIndex(irow)>0){
	row[1] = irow;
	break;
      }
    }
    //
  }
  if (!forward){
    for (Int_t irow=0;irow<r0;irow++){
      if (track->GetClusterIndex(irow)>0){
	row[0] = irow;
	break;
      }
    }
    for (Int_t irow=r0;irow>0;irow--){
      if (track->GetClusterIndex(irow)>0){
	row[2] = irow;
	break;
      }
    }    
    for (Int_t irow=row[2]-15;irow>row[0];irow--){
      if (track->GetClusterIndex(irow)>0){
	row[1] = irow;
	break;
      }
    } 
  }
  //
  if ((row[2]-row[0])<20) return 0;
  if (row[1]==0) return 0;
  //
  //
  //Get cluster and sector position
  for (Int_t ipoint=0;ipoint<3;ipoint++){
    Int_t clindex = track->GetClusterIndex2(row[ipoint]);
    AliTPCclusterMI * cl = clindex<0 ? 0:GetClusterMI(clindex);
    if (cl==0) {
      //Error("Bug\n");
      //      AliTPCclusterMI * cl = GetClusterMI(clindex);
      return 0;
    }
    sec[ipoint]     = ((clindex&0xff000000)>>24)%18;
    xyz[ipoint][0]  = GetXrow(row[ipoint]);
    const AliTPCTrackerPoints::Point * point = track->GetTrackPoint(row[ipoint]);    
    if (point&&ipoint<2){
      //
       xyz[ipoint][1]  = point->GetY();
       xyz[ipoint][2]  = point->GetZ();
    }
    else{
      xyz[ipoint][1]  = cl->GetY();
      xyz[ipoint][2]  = cl->GetZ();
    }
  }
  //
  //
  //
  //
  // Calculate seed state vector and covariance matrix

  Double_t alpha, cs,sn, xx2,yy2;
  //
  alpha = (sec[1]-sec[2])*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= xyz[1][0]*cs-xyz[1][1]*sn;
  yy2= xyz[1][0]*sn+xyz[1][1]*cs;
  xyz[1][0] = xx2;
  xyz[1][1] = yy2;
  //
  alpha = (sec[0]-sec[2])*fSectors->GetAlpha();
  cs = TMath::Cos(alpha);
  sn = TMath::Sin(alpha); 
  xx2= xyz[0][0]*cs-xyz[0][1]*sn;
  yy2= xyz[0][0]*sn+xyz[0][1]*cs;
  xyz[0][0] = xx2;
  xyz[0][1] = yy2;
  //
  //
  //
  Double_t x[5],c[15];
  //
  x[0]=xyz[2][1];
  x[1]=xyz[2][2];
  x[4]=F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[2]=F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]);
  x[3]=F3n(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2],x[4]);
  //  
  Double_t sy =0.1,  sz =0.1;
  //
  Double_t sy1=0.2, sz1=0.2;
  Double_t sy2=0.2, sz2=0.2;
  Double_t sy3=0.2;
  //
  Double_t f40=(F1(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f42=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[4])/sy;
  Double_t f43=(F1(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[4])/sy;
  Double_t f20=(F2(xyz[2][0],xyz[2][1]+sy,xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f22=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1]+sy,xyz[0][0],xyz[0][1])-x[2])/sy;
  Double_t f23=(F2(xyz[2][0],xyz[2][1],xyz[1][0],xyz[1][1],xyz[0][0],xyz[0][1]+sy)-x[2])/sy;
  //
  Double_t f30=(F3(xyz[2][0],xyz[2][1]+sy,xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f31=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2]+sz,xyz[0][2])-x[3])/sz;
  Double_t f32=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1]+sy,xyz[2][2],xyz[0][2])-x[3])/sy;
  Double_t f34=(F3(xyz[2][0],xyz[2][1],xyz[0][0],xyz[0][1],xyz[2][2],xyz[0][2]+sz)-x[3])/sz;
  
  
  c[0]=sy1;
  c[1]=0.;       c[2]=sz1;
  c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
  c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
  c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
  c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
  c[13]=f30*sy1*f40+f32*sy2*f42;
  c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;
  
  //  Int_t row1 = fSectors->GetRowNumber(xyz[2][0]);
  AliTPCseed *seed=new( NextFreeSeed() )  AliTPCseed(xyz[2][0], sec[2]*fSectors->GetAlpha()+fSectors->GetAlphaShift(), x, c, 0);
  seed->SetPoolID(fLastSeedID);
  seed->SetRow(row[2]); //RS: memorise row
  seed->SetLastPoint(row[2]);
  seed->SetFirstPoint(row[2]);  
  for (Int_t i=row[0];i<row[2];i++){
    seed->SetClusterIndex(i, track->GetClusterIndex(i));
  }

  return seed;
}



void  AliTPCtracker::FindMultiMC(const TObjArray * array, AliESDEvent */*esd*/, Int_t iter)
{
  //
  //  find multi tracks - THIS FUNCTION IS ONLY FOR DEBUG PURPOSES
  //                      USES MC LABELS
  //  Use AliTPCReconstructor::StreamLevel()& kStreamFindMultiMC if you want to tune parameters - cuts
  //
  //  Two reasons to have multiple find tracks
  //  1. Curling tracks can be find more than once
  //  2. Splitted tracks 
  //     a.) Multiple seeding to increase tracking efficiency - (~ 100% reached)        
  //     b.) Edge effect on the sector boundaries
  //
  //
  //  Algorithm done in 2 phases - because of CPU consumption
  //  it is n^2 algorithm - for lead-lead 20000x20000 combination are investigated                           
  //
  //  Algorihm for curling tracks sign:
  //    1 phase -makes a very rough fast cuts to minimize combinatorics
  //                 a.) opposite sign
  //                 b.) one of the tracks - not pointing to the primary vertex - 
  //                 c.) delta tan(theta)
  //                 d.) delta phi
  //    2 phase - calculates DCA between tracks  - time consument

  //
  //    fast cuts 
  //
  //    General cuts    - for splitted tracks and for curling tracks
  //
  const Float_t kMaxdPhi      = 0.2;  // maximal distance in phi
  //
  //    Curling tracks cuts
  //
  //
  //
  //
  Int_t nentries = array->GetEntriesFast();  
  if (!fHelixPool) fHelixPool = new TClonesArray("AliHelix",nentries+50);
  fHelixPool->Clear();
  TClonesArray& helixes = *fHelixPool;
  Float_t  xm[nentries];
  Float_t  dz0[nentries];
  Float_t  dz1[nentries];
  //
  //
  TStopwatch timer;
  timer.Start();
  //
  // Find track COG in x direction - point with best defined parameters
  //
  for (Int_t i=0;i<nentries;i++){
    AliTPCseed* track = (AliTPCseed*)array->At(i);    
    if (!track) continue;
    track->SetCircular(0);
    new (helixes[i]) AliHelix(*track);
    Int_t ncl=0;
    xm[i]=0;
    Float_t dz[2];
    track->GetDZ(GetX(),GetY(),GetZ(),GetBz(),dz);
    dz0[i]=dz[0];
    dz1[i]=dz[1];
    for (Int_t icl=0; icl<kMaxRow; icl++){
      int tpcindex= track->GetClusterIndex2(icl);
      const AliTPCclusterMI * cl = (tpcindex<0) ? 0:GetClusterMI(tpcindex);
      //RS      AliTPCclusterMI * cl =  track->GetClusterPointer(icl);
      if (cl) {
	xm[i]+=cl->GetX();
	ncl++;
      }
    }
    if (ncl>0) xm[i]/=Float_t(ncl);
  }  
  //
  for (Int_t i0=0;i0<nentries;i0++){
    AliTPCseed * track0 = (AliTPCseed*)array->At(i0);
    if (!track0) continue;    
    AliHelix* hlxi0 = (AliHelix*)helixes[i0];
    Float_t xc0 = hlxi0->GetHelix(6);
    Float_t yc0 = hlxi0->GetHelix(7);
    Float_t r0  = hlxi0->GetHelix(8);
    Float_t rc0 = TMath::Sqrt(xc0*xc0+yc0*yc0);
    Float_t fi0 = TMath::ATan2(yc0,xc0);
    
    for (Int_t i1=i0+1;i1<nentries;i1++){
      AliTPCseed * track1 = (AliTPCseed*)array->At(i1);
      if (!track1) continue;      
      Int_t lab0=track0->GetLabel();
      Int_t lab1=track1->GetLabel();
      if (TMath::Abs(lab0)!=TMath::Abs(lab1)) continue;
      //
      AliHelix* hlxi1 = (AliHelix*)helixes[i1];
      Float_t xc1 = hlxi1->GetHelix(6);
      Float_t yc1 = hlxi1->GetHelix(7);
      Float_t r1  = hlxi1->GetHelix(8);
      Float_t rc1 = TMath::Sqrt(xc1*xc1+yc1*yc1);
      Float_t fi1 = TMath::ATan2(yc1,xc1);
      //
      Float_t dfi = fi0-fi1;
      //
      //
      if (dfi>1.5*TMath::Pi())  dfi-=TMath::Pi();  // take care about edge effect 
      if (dfi<-1.5*TMath::Pi()) dfi+=TMath::Pi();  // 
      if (TMath::Abs(dfi)>kMaxdPhi&&hlxi0->GetHelix(4)*hlxi1->GetHelix(4)<0){
	//
	// if short tracks with undefined sign 
	fi1 =  -TMath::ATan2(yc1,-xc1);
	dfi = fi0-fi1;
      }
      Float_t dtheta = TMath::Abs(track0->GetTgl()-track1->GetTgl())<TMath::Abs(track0->GetTgl()+track1->GetTgl())? track0->GetTgl()-track1->GetTgl():track0->GetTgl()+track1->GetTgl();
      
      //
      // debug stream to tune "fast cuts" 
      //
      Double_t dist[3];   // distance at X 
      Double_t mdist[3]={0,0,0};  // mean distance X+-40cm
      track0->GetDistance(track1,0.5*(xm[i0]+xm[i1])-40.,dist,AliTracker::GetBz());
      for (Int_t i=0;i<3;i++) mdist[i]+=TMath::Abs(dist[i]);
      track0->GetDistance(track1,0.5*(xm[i0]+xm[i1])+40.,dist,AliTracker::GetBz());
      for (Int_t i=0;i<3;i++) mdist[i]+=TMath::Abs(dist[i]);
      track0->GetDistance(track1,0.5*(xm[i0]+xm[i1]),dist,AliTracker::GetBz());
      for (Int_t i=0;i<3;i++) mdist[i]+=TMath::Abs(dist[i]);
      for (Int_t i=0;i<3;i++) mdist[i]*=0.33333;
      
      Float_t sum =0;
      Float_t sums=0;
      for (Int_t icl=0; icl<kMaxRow; icl++){
	Int_t tpcindex0 = track0->GetClusterIndex2(icl);
	Int_t tpcindex1 = track1->GetClusterIndex2(icl);
	//RS AliTPCclusterMI * cl0 =  track0->GetClusterPointer(icl);
	//RS AliTPCclusterMI * cl1 =  track1->GetClusterPointer(icl);
	//RS if (cl0&&cl1) {
	if (tpcindex0>=0 && tpcindex1>=0) {
	  sum++;
	  if (tpcindex0==tpcindex1) sums++; //RS if (cl0==cl1) sums++;
	}
      }
      //
      if ((AliTPCReconstructor::StreamLevel()&kStreamFindMultiMC)>0) {  // flag: stream MC infomation about the multiple find track (ONLY for MC data)
      TTreeSRedirector &cstream = *fDebugStreamer;
      cstream<<"Multi"<<
	"iter="<<iter<<
	"lab0="<<lab0<<
	"lab1="<<lab1<<   
	"Tr0.="<<track0<<       // seed0
	"Tr1.="<<track1<<       // seed1
	"h0.="<<hlxi0<<
	"h1.="<<hlxi1<<
	//
	"sum="<<sum<<           //the sum of rows with cl in both
	"sums="<<sums<<         //the sum of shared clusters
	"xm0="<<xm[i0]<<        // the center of track
	"xm1="<<xm[i1]<<        // the x center of track
	// General cut variables                   
	"dfi="<<dfi<<           // distance in fi angle
	"dtheta="<<dtheta<<     // distance int theta angle
	//
	"dz00="<<dz0[i0]<<
	"dz01="<<dz0[i1]<<
	"dz10="<<dz1[i1]<<
	"dz11="<<dz1[i1]<<
	"dist0="<<dist[0]<<     //distance x
	"dist1="<<dist[1]<<     //distance y
	"dist2="<<dist[2]<<     //distance z
	"mdist0="<<mdist[0]<<   //distance x
	"mdist1="<<mdist[1]<<   //distance y
	"mdist2="<<mdist[2]<<   //distance z
	//
	"r0="<<r0<<
	"rc0="<<rc0<<
	"fi0="<<fi0<<
	"fi1="<<fi1<<
	"r1="<<r1<<
	"rc1="<<rc1<<
	"\n";
	}
    }
  }    
  if (fHelixPool) fHelixPool->Clear();
  //  delete [] helixes; // RS moved to stack
  //  delete [] xm;
  //  delete [] dz0;
  //  delete [] dz1;
  if (AliTPCReconstructor::StreamLevel()>0) {
    AliInfo("Time for curling tracks removal DEBUGGING MC");
    timer.Print();
  }
}



void  AliTPCtracker::FindSplitted(TObjArray * array, AliESDEvent */*esd*/, Int_t /*iter*/){
  //
  // Find Splitted tracks and remove the one with worst quality  
  // Corresponding debug streamer to tune selections - "Splitted2"
  // Algorithm:
  // 0. Sort tracks according quality
  // 1. Propagate the tracks to the reference radius
  // 2. Double_t loop to select close tracks (only to speed up process)
  // 3. Calculate cluster overlap ratio - and remove the track if bigger than a threshold
  // 4. Delete temporary parameters
  // 
  const Double_t xref=GetXrow(63);  // reference radius -IROC/OROC boundary
  //    rough cuts
  const Double_t kCutP1=10;       // delta Z cut 10 cm 
  const Double_t kCutP2=0.15;     // delta snp(fi) cut 0.15 
  const Double_t kCutP3=0.15;     // delta tgl(theta) cut 0.15
  const Double_t kCutAlpha=0.15;  // delta alpha cut
  Int_t firstpoint = 0;
  Int_t lastpoint = kMaxRow;
  //
  Int_t nentries = array->GetEntriesFast();  
  if (!fETPPool) fETPPool = new TClonesArray("AliExternalTrackParam",nentries+50);
  else fETPPool->Clear();
  TClonesArray &params = *fETPPool;
  //
  //
  TStopwatch timer;
  timer.Start();
  //
  //0.  Sort tracks according quality
  //1.  Propagate the ext. param to reference radius
  Int_t nseed = array->GetEntriesFast();  
  if (nseed<=0) return;
  Float_t quality[nseed];
  Int_t   indexes[nseed];
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)array->UncheckedAt(i);    
    if (!pt){
      quality[i]=-1;
      continue;
    }
    pt->UpdatePoints();    //select first last max dens points
    Float_t * points = pt->GetPoints();
    if (points[3]<0.8) quality[i] =-1;
    quality[i] = (points[2]-points[0])+pt->GetNumberOfClusters();
    //prefer high momenta tracks if overlaps
    quality[i] *= TMath::Sqrt(TMath::Abs(pt->Pt())+0.5);
    AliExternalTrackParam* parI = new (params[i]) AliExternalTrackParam(*pt);
    // params[i]=(*pt);
    AliTracker::PropagateTrackParamOnlyToBxByBz(parI,xref,5.,kTRUE);
    //    AliTracker::PropagateTrackToBxByBz(parI,xref,pt->GetMass(),1.,kTRUE); //RS What is the point of 2nd propagation
  }
  TMath::Sort(nseed,quality,indexes);
  //
  // 3. Loop over pair of tracks
  //
  for (Int_t i0=0; i0<nseed; i0++) {
    Int_t index0=indexes[i0];
    if (!(array->UncheckedAt(index0))) continue;
    AliTPCseed *s1 = (AliTPCseed*)array->UncheckedAt(index0);  
    if (!s1->IsActive()) continue;
    AliExternalTrackParam &par0=*(AliExternalTrackParam*)params[index0];
    for (Int_t i1=i0+1; i1<nseed; i1++) {
      Int_t index1=indexes[i1];
      if (!(array->UncheckedAt(index1))) continue;
      AliTPCseed *s2 = (AliTPCseed*)array->UncheckedAt(index1);  
      if (!s2->IsActive()) continue;
      if (s2->GetKinkIndexes()[0]!=0)
	if (s2->GetKinkIndexes()[0] == -s1->GetKinkIndexes()[0]) continue;
      AliExternalTrackParam &par1=*(AliExternalTrackParam*)params[index1];
      if (TMath::Abs(par0.GetParameter()[3]-par1.GetParameter()[3])>kCutP3) continue;
      if (TMath::Abs(par0.GetParameter()[1]-par1.GetParameter()[1])>kCutP1) continue;
      if (TMath::Abs(par0.GetParameter()[2]-par1.GetParameter()[2])>kCutP2) continue;
      Double_t dAlpha= TMath::Abs(par0.GetAlpha()-par1.GetAlpha());
      if (dAlpha>TMath::Pi()) dAlpha-=TMath::Pi();
      if (TMath::Abs(dAlpha)>kCutAlpha) continue;
      //
      Int_t sumShared=0;
      Int_t nall0=0;
      Int_t nall1=0;
      Int_t firstShared=lastpoint, lastShared=firstpoint;
      Int_t firstRow=lastpoint, lastRow=firstpoint;
      //
      // for (Int_t i=firstpoint;i<lastpoint;i++){
      // 	if (s1->GetClusterIndex2(i)>0) nall0++;
      // 	if (s2->GetClusterIndex2(i)>0) nall1++;
      // 	if  (s1->GetClusterIndex2(i)>0 && s2->GetClusterIndex2(i)>0) {
      // 	  if (i<firstRow) firstRow=i;
      // 	  if (i>lastRow)  lastRow=i;
      // 	}
      // 	if ( (s1->GetClusterIndex2(i))==(s2->GetClusterIndex2(i)) && s1->GetClusterIndex2(i)>0) {
      // 	  if (i<firstShared) firstShared=i;
      // 	  if (i>lastShared)  lastShared=i;
      // 	  sumShared++;
      // 	}
      // }
      //
      // RS: faster version + fix(?)  ">" -> ">="
      for (Int_t i=firstpoint;i<lastpoint;i++){
	int ind1=s1->GetClusterIndex2(i),ind2=s2->GetClusterIndex2(i);
	if (ind1>=0) nall0++; // RS: ">" -> ">="
	if (ind2>=0) nall1++;
	if  (ind1>=0 && ind2>=0) {
	  if (i<firstRow) firstRow=i;
	  if (i>lastRow)  lastRow=i;
	}
	if ( (ind1==ind2) && ind1>=0) {
	  if (i<firstShared) firstShared=i;
	  if (i>lastShared)  lastShared=i;
	  sumShared++;
	}
      }

      Double_t ratio0 = Float_t(sumShared)/Float_t(TMath::Min(nall0+1,nall1+1));
      Double_t ratio1 = Float_t(sumShared)/Float_t(TMath::Max(nall0+1,nall1+1));
      
      if ((AliTPCReconstructor::StreamLevel()&kStreamSplitted2)>0){ // flag:stream information about discarded TPC tracks pair algorithm 
	TTreeSRedirector &cstream = *fDebugStreamer;
	Int_t n0=s1->GetNumberOfClusters();
	Int_t n1=s2->GetNumberOfClusters();
	Int_t n0F=s1->GetNFoundable();
	Int_t n1F=s2->GetNFoundable();
	Int_t lab0=s1->GetLabel();
	Int_t lab1=s2->GetLabel();

	cstream<<"Splitted2"<<    // flag:stream information about discarded TPC tracks pair algorithm 
	  "iter="<<fIteration<<
	  "lab0="<<lab0<<        // MC label if exist
	  "lab1="<<lab1<<        // MC label if exist
	  "index0="<<index0<<
	  "index1="<<index1<<
 	  "ratio0="<<ratio0<<      // shared ratio
 	  "ratio1="<<ratio1<<      // shared ratio
	  "p0.="<<&par0<<        // track parameters
	  "p1.="<<&par1<<
	  "s0.="<<s1<<           // full seed
	  "s1.="<<s2<<
	  "n0="<<n0<<     // number of clusters track 0
	  "n1="<<n1<<     // number of clusters track 1
	  "nall0="<<nall0<<     // number of clusters track 0
	  "nall1="<<nall1<<     // number of clusters track 1
	  "n0F="<<n0F<<   // number of findable
	  "n1F="<<n1F<<   // number of findable
	  "shared="<<sumShared<<    // number of shared clusters
	  "firstS="<<firstShared<<  // first and the last shared row
	  "lastS="<<lastShared<<
	  "firstRow="<<firstRow<<   // first and the last row with cluster
	  "lastRow="<<lastRow<<     //
	  "\n";
      }
      //
      // remove track with lower quality
      //
      if (ratio0>AliTPCReconstructor::GetRecoParam()->GetCutSharedClusters(0) ||
	  ratio1>AliTPCReconstructor::GetRecoParam()->GetCutSharedClusters(1)){
	//
	//
	//
	MarkSeedFree( array->RemoveAt(index1) );
      }
    }
  }
  //
  // 4. Delete temporary array
  //
  if (fETPPool) fETPPool->Clear();
  // delete [] params; // RS moved to stack
  //  delete [] quality;
  //  delete [] indexes;

}



void  AliTPCtracker::FindCurling(const TObjArray * array, AliESDEvent */*esd*/, Int_t iter)
{
  //
  //  find Curling tracks
  //  Use AliTPCReconstructor::StreamLevel()&kStreamFindCurling if you want to tune parameters - cuts
  //
  //
  //  Algorithm done in 2 phases - because of CPU consumption
  //  it is n^2 algorithm - for lead-lead 20000x20000 combination are investigated                           
  //  see detal in MC part what can be used to cut
  //
  //    
  //
  const Float_t kMaxC         = 400;  // maximal curvature to of the track
  const Float_t kMaxdTheta    = 0.15;  // maximal distance in theta
  const Float_t kMaxdPhi      = 0.15;  // maximal distance in phi
  const Float_t kPtRatio      = 0.3; // ratio between pt
  const Float_t kMinDCAR      = 2.;   // distance to the primary vertex in r - see cpipe cut      

  //
  //    Curling tracks cuts
  //
  //
  const Float_t kMaxDeltaRMax = 40;   // distance in outer radius
  const Float_t kMaxDeltaRMin = 5.;   // distance in lower radius - see cpipe cut
  const Float_t kMinAngle     = 2.9;  // angle between tracks
  const Float_t kMaxDist      = 5;    // biggest distance 
  //
  // The cuts can be tuned using the "MC information stored in Multi tree ==> see FindMultiMC
  /* 
     Fast cuts:
     TCut csign("csign","Tr0.fP[4]*Tr1.fP[4]<0"); //opposite sign
     TCut cmax("cmax","abs(Tr0.GetC())>1/400");
     TCut cda("cda","sqrt(dtheta^2+dfi^2)<0.15");
     TCut ccratio("ccratio","abs((Tr0.fP[4]+Tr1.fP[4])/(abs(Tr0.fP[4])+abs(Tr1.fP[4])))<0.3");
     TCut cpipe("cpipe", "min(abs(r0-rc0),abs(r1-rc1))>5");    
     //
     TCut cdrmax("cdrmax","abs(abs(rc0+r0)-abs(rc1+r1))<40")
     TCut cdrmin("cdrmin","abs(abs(rc0+r0)-abs(rc1+r1))<10")
     //
     Multi->Draw("dfi","iter==0"+csign+cmax+cda+ccratio); ~94% of curling tracks fulfill 
     Multi->Draw("min(abs(r0-rc0),abs(r1-rc1))","iter==0&&abs(lab1)==abs(lab0)"+csign+cmax+cda+ccratio+cpipe+cdrmin+cdrmax); //80%
     //
     Curling2->Draw("dfi","iter==0&&abs(lab0)==abs(lab1)"+csign+cmax+cdtheta+cdfi+ccratio)

  */
  //  
  //
  //
  Int_t nentries = array->GetEntriesFast();
  if (!fHelixPool) fHelixPool = new TClonesArray("AliHelix",nentries+100);
  fHelixPool->Clear();
  TClonesArray& helixes = *fHelixPool;

  for (Int_t i=0;i<nentries;i++){
    AliTPCseed* track = (AliTPCseed*)array->At(i);    
    if (!track) continue;
    track->SetCircular(0);
    new (helixes[i]) AliHelix(*track);
  }
  //
  //
  TStopwatch timer;
  timer.Start();
  Double_t phase[2][2]={{0,0},{0,0}},radius[2]={0,0};

  //
  // Find tracks
  //
  //
  for (Int_t i0=0;i0<nentries;i0++){
    AliTPCseed * track0 = (AliTPCseed*)array->At(i0);
    if (!track0) continue;    
    if (TMath::Abs(track0->GetC())<1/kMaxC) continue;
    AliHelix* hlxi0 = (AliHelix*)helixes[i0];
    Float_t xc0 = hlxi0->GetHelix(6);
    Float_t yc0 = hlxi0->GetHelix(7);
    Float_t r0  = hlxi0->GetHelix(8);
    Float_t rc0 = TMath::Sqrt(xc0*xc0+yc0*yc0);
    Float_t fi0 = TMath::ATan2(yc0,xc0);
    
    for (Int_t i1=i0+1;i1<nentries;i1++){
      AliTPCseed * track1 = (AliTPCseed*)array->At(i1);
      if (!track1) continue;      
      if (TMath::Abs(track1->GetC())<1/kMaxC) continue;    
      AliHelix* hlxi1 = (AliHelix*)helixes[i1];
      Float_t xc1 = hlxi1->GetHelix(6);
      Float_t yc1 = hlxi1->GetHelix(7);
      Float_t r1  = hlxi1->GetHelix(8);
      Float_t rc1 = TMath::Sqrt(xc1*xc1+yc1*yc1);
      Float_t fi1 = TMath::ATan2(yc1,xc1);
      //
      Float_t dfi = fi0-fi1;
      //
      //
      if (dfi>1.5*TMath::Pi())  dfi-=TMath::Pi();  // take care about edge effect 
      if (dfi<-1.5*TMath::Pi()) dfi+=TMath::Pi();  // 
      Float_t dtheta = TMath::Abs(track0->GetTgl()-track1->GetTgl())<TMath::Abs(track0->GetTgl()+track1->GetTgl())? track0->GetTgl()-track1->GetTgl():track0->GetTgl()+track1->GetTgl();
      //
      //
      // FIRST fast cuts
      if (track0->GetBConstrain()&&track1->GetBConstrain())  continue;  // not constrained
      if (track1->GetSigned1Pt()*track0->GetSigned1Pt()>0)   continue; // not the same sign
      if ( TMath::Abs(track1->GetTgl()+track0->GetTgl())>kMaxdTheta) continue; //distance in the Theta
      if ( TMath::Abs(dfi)>kMaxdPhi) continue;  //distance in phi
      if ( TMath::Sqrt(dfi*dfi+dtheta*dtheta)>kMaxdPhi) continue; //common angular offset
      //
      Float_t pt0 = track0->GetSignedPt();
      Float_t pt1 = track1->GetSignedPt();
      if ((TMath::Abs(pt0+pt1)/(TMath::Abs(pt0)+TMath::Abs(pt1)))>kPtRatio) continue;      
      if ((iter==1) && TMath::Abs(TMath::Abs(rc0+r0)-TMath::Abs(rc1+r1))>kMaxDeltaRMax) continue;
      if ((iter!=1) &&TMath::Abs(TMath::Abs(rc0-r0)-TMath::Abs(rc1-r1))>kMaxDeltaRMin) continue;
      if (TMath::Min(TMath::Abs(rc0-r0),TMath::Abs(rc1-r1))<kMinDCAR) continue;
      //
      //
      // Now find closest approach
      //
      //
      //
      Int_t npoints = hlxi0->GetRPHIintersections(*hlxi1, phase, radius,10);
      if (npoints==0) continue;
      hlxi0->GetClosestPhases(*hlxi1, phase);
      //
      Double_t xyz0[3];
      Double_t xyz1[3];
      Double_t hangles[3];
      hlxi0->Evaluate(phase[0][0],xyz0);
      hlxi1->Evaluate(phase[0][1],xyz1);

      hlxi0->GetAngle(phase[0][0],*hlxi1,phase[0][1],hangles);
      Double_t deltah[2],deltabest;
      if (TMath::Abs(hangles[2])<kMinAngle) continue;

      if (npoints>0){
	Int_t ibest=0;
	hlxi0->ParabolicDCA(*hlxi1,phase[0][0],phase[0][1],radius[0],deltah[0],2);
	if (npoints==2){
	  hlxi0->ParabolicDCA(*hlxi1,phase[1][0],phase[1][1],radius[1],deltah[1],2);
	  if (deltah[1]<deltah[0]) ibest=1;
	}
	deltabest  = TMath::Sqrt(deltah[ibest]);
	hlxi0->Evaluate(phase[ibest][0],xyz0);
	hlxi1->Evaluate(phase[ibest][1],xyz1);
	hlxi0->GetAngle(phase[ibest][0],*hlxi1,phase[ibest][1],hangles);
	Double_t radiusbest = TMath::Sqrt(radius[ibest]);
	//
	if (deltabest>kMaxDist) continue;
	//	if (mindcar+mindcaz<40 && (TMath::Abs(hangles[2])<kMinAngle ||deltabest>3)) continue;
	Bool_t sign =kFALSE;
	if (hangles[2]>kMinAngle) sign =kTRUE;
	//
	if (sign){
	  //	  circular[i0] = kTRUE;
	  //      circular[i1] = kTRUE;
	  if (track0->OneOverPt()<track1->OneOverPt()){
	    track0->SetCircular(track0->GetCircular()+1);
	    track1->SetCircular(track1->GetCircular()+2);
	  }
	  else{
	    track1->SetCircular(track1->GetCircular()+1);
	    track0->SetCircular(track0->GetCircular()+2);
	  }
	}		
	if ((AliTPCReconstructor::StreamLevel()&kStreamFindCurling)>0){  // flag: stream track infroamtion if the FindCurling tracks method	  
	  //
	  //debug stream to tune "fine" cuts	  
	  Int_t lab0=track0->GetLabel();
	  Int_t lab1=track1->GetLabel();
          TTreeSRedirector &cstream = *fDebugStreamer;
	  cstream<<"Curling2"<<
	    "iter="<<iter<<
	    "lab0="<<lab0<<
	    "lab1="<<lab1<<   
	    "Tr0.="<<track0<<
	    "Tr1.="<<track1<<
	    //
	    "r0="<<r0<<
	    "rc0="<<rc0<<
	    "fi0="<<fi0<<
	    "r1="<<r1<<
	    "rc1="<<rc1<<
	    "fi1="<<fi1<<
	    "dfi="<<dfi<<
	    "dtheta="<<dtheta<<
	    //
	    "npoints="<<npoints<<                      
	    "hangles0="<<hangles[0]<<
	    "hangles1="<<hangles[1]<<
	    "hangles2="<<hangles[2]<<                    
	    "xyz0="<<xyz0[2]<<
	    "xyzz1="<<xyz1[2]<<
	    "radius="<<radiusbest<<
	    "deltabest="<<deltabest<< 
	    "phase0="<<phase[ibest][0]<<
	    "phase1="<<phase[ibest][1]<<
	    "\n"; 	  	  

	}
      }
    }
  }
  if (fHelixPool) fHelixPool->Clear();
  //  delete [] helixes; //RS moved to stack
  if (AliTPCReconstructor::StreamLevel()>1) {
    AliInfo("Time for curling tracks removal");
    timer.Print();
  }
}


void  AliTPCtracker::FindKinks(TObjArray * array, AliESDEvent *esd)
{
  //
  //  find kinks
  //
  //
  // RS something is wrong in this routine: not all seeds are assigned to daughters and mothers array, but they all are queried
  // to check later

  TObjArray kinks(10000);
  Int_t nentries = array->GetEntriesFast();
  Char_t   sign[nentries];
  UChar_t  nclusters[nentries];
  Float_t  alpha[nentries];
  UChar_t  usage[nentries];
  Float_t  zm[nentries];
  Float_t  z0[nentries]; 
  Float_t  fim[nentries];
  Bool_t   shared[nentries];
  Bool_t   circular[nentries];
  Float_t dca[nentries];
  AliKink  *kink         = new AliKink();
  //
  if (!fHelixPool) fHelixPool = new TClonesArray("AliHelix",nentries+100);
  fHelixPool->Clear();
  TClonesArray& helixes = *fHelixPool;

  //const AliESDVertex * primvertex = esd->GetVertex();
  //
  //  nentries = array->GetEntriesFast();
  //
  
  //
  //
  for (Int_t i=0;i<nentries;i++){
    sign[i]=0;
    usage[i]=0;
    AliTPCseed* track = (AliTPCseed*)array->At(i);    
    if (!track) continue;
    track->SetCircular(0);
    shared[i] = kFALSE;
    track->UpdatePoints();
    if (( track->GetPoints()[2]- track->GetPoints()[0])>5 && track->GetPoints()[3]>0.8){
    }
    nclusters[i]=track->GetNumberOfClusters();
    alpha[i] = track->GetAlpha();
    AliHelix* hlxi = new (helixes[i]) AliHelix(*track);
    Double_t xyz[3];
    hlxi->Evaluate(0,xyz);
    sign[i] = (track->GetC()>0) ? -1:1;
    Double_t x,y,z;
    x=160;
    if (track->GetProlongation(x,y,z)){
      zm[i]  = z;
      fim[i] = alpha[i]+TMath::ATan2(y,x);
    }
    else{
      zm[i]  = track->GetZ();
      fim[i] = alpha[i];
    }   
    z0[i]=1000;
    circular[i]= kFALSE;
    if (track->GetProlongation(0,y,z))  z0[i] = z;
    dca[i] = track->GetD(0,0);    
  }
  //
  //
  TStopwatch timer;
  timer.Start();
  Int_t ncandidates =0;
  Int_t nall =0;
  Int_t ntracks=0; 
  Double_t phase[2][2]={{0,0},{0,0}},radius[2]={0,0};

  //
  // Find circling track
  //
  for (Int_t i0=0;i0<nentries;i0++){
    AliTPCseed * track0 = (AliTPCseed*)array->At(i0);
    if (!track0) continue;    
    if (track0->GetNumberOfClusters()<40) continue;
    if (TMath::Abs(1./track0->GetC())>200) continue;
    AliHelix* hlxi0 = (AliHelix*)helixes[i0];
    //
    for (Int_t i1=i0+1;i1<nentries;i1++){
      AliTPCseed * track1 = (AliTPCseed*)array->At(i1);
      if (!track1) continue;
      if (track1->GetNumberOfClusters()<40)                  continue;
      if ( TMath::Abs(track1->GetTgl()+track0->GetTgl())>0.1) continue;
      if (track0->GetBConstrain()&&track1->GetBConstrain()) continue;
      if (TMath::Abs(1./track1->GetC())>200) continue;
      if (track1->GetSigned1Pt()*track0->GetSigned1Pt()>0)      continue;
      if (track1->GetTgl()*track0->GetTgl()>0)      continue;
      if (TMath::Max(TMath::Abs(1./track0->GetC()),TMath::Abs(1./track1->GetC()))>190) continue;
      if (track0->GetBConstrain()&&track1->OneOverPt()<track0->OneOverPt()) continue; //returning - lower momenta
      if (track1->GetBConstrain()&&track0->OneOverPt()<track1->OneOverPt()) continue; //returning - lower momenta
      //
      Float_t mindcar = TMath::Min(TMath::Abs(dca[i0]),TMath::Abs(dca[i1]));
      if (mindcar<5)   continue;
      Float_t mindcaz = TMath::Min(TMath::Abs(z0[i0]-GetZ()),TMath::Abs(z0[i1]-GetZ()));
      if (mindcaz<5) continue;
      if (mindcar+mindcaz<20) continue;
      //
      AliHelix* hlxi1 = (AliHelix*)helixes[i1];
      //
      Float_t xc0 = hlxi0->GetHelix(6);
      Float_t yc0 = hlxi0->GetHelix(7);
      Float_t r0  = hlxi0->GetHelix(8);
      Float_t xc1 = hlxi1->GetHelix(6);
      Float_t yc1 = hlxi1->GetHelix(7);
      Float_t r1  = hlxi1->GetHelix(8);
	
      Float_t rmean = (r0+r1)*0.5;
      Float_t delta =TMath::Sqrt((xc1-xc0)*(xc1-xc0)+(yc1-yc0)*(yc1-yc0));
      //if (delta>30) continue;
      if (delta>rmean*0.25) continue;
      if (TMath::Abs(r0-r1)/rmean>0.3) continue; 
      //
      Int_t npoints = hlxi0->GetRPHIintersections(*hlxi1, phase, radius,10);
      if (npoints==0) continue;
      hlxi0->GetClosestPhases(*hlxi1, phase);
      //
      Double_t xyz0[3];
      Double_t xyz1[3];
      Double_t hangles[3];
      hlxi0->Evaluate(phase[0][0],xyz0);
      hlxi1->Evaluate(phase[0][1],xyz1);

      hlxi0->GetAngle(phase[0][0],*hlxi1,phase[0][1],hangles);
      Double_t deltah[2],deltabest;
      if (hangles[2]<2.8) continue;
      if (npoints>0){
	Int_t ibest=0;
	hlxi0->ParabolicDCA(*hlxi1,phase[0][0],phase[0][1],radius[0],deltah[0],2);
	if (npoints==2){
	  hlxi0->ParabolicDCA(*hlxi1,phase[1][0],phase[1][1],radius[1],deltah[1],2);
	  if (deltah[1]<deltah[0]) ibest=1;
	}
	deltabest  = TMath::Sqrt(deltah[ibest]);
	hlxi0->Evaluate(phase[ibest][0],xyz0);
	hlxi1->Evaluate(phase[ibest][1],xyz1);
	hlxi0->GetAngle(phase[ibest][0],*hlxi1,phase[ibest][1],hangles);
	Double_t radiusbest = TMath::Sqrt(radius[ibest]);
	//
	if (deltabest>6) continue;
	if (mindcar+mindcaz<40 && (hangles[2]<3.12||deltabest>3)) continue;
	Bool_t lsign =kFALSE;
	if (hangles[2]>3.06) lsign =kTRUE;
	//
	if (lsign){
	  circular[i0] = kTRUE;
	  circular[i1] = kTRUE;
	  if (track0->OneOverPt()<track1->OneOverPt()){
	    track0->SetCircular(track0->GetCircular()+1);
	    track1->SetCircular(track1->GetCircular()+2);
	  }
	  else{
	    track1->SetCircular(track1->GetCircular()+1);
	    track0->SetCircular(track0->GetCircular()+2);
	  }
	}		
	if (lsign&&((AliTPCReconstructor::StreamLevel()&kStreamFindCurling)>0)){	  
	  //debug stream	  
	  Int_t lab0=track0->GetLabel();
	  Int_t lab1=track1->GetLabel();
          TTreeSRedirector &cstream = *fDebugStreamer;
	  cstream<<"Curling"<<
	    "lab0="<<lab0<<
	    "lab1="<<lab1<<   
	    "Tr0.="<<track0<<
	    "Tr1.="<<track1<<	   
	    "dca0="<<dca[i0]<<
	    "dca1="<<dca[i1]<<
	    "mindcar="<<mindcar<<
	    "mindcaz="<<mindcaz<<
	    "delta="<<delta<<
	    "rmean="<<rmean<<
	    "npoints="<<npoints<<                      
	    "hangles0="<<hangles[0]<<
	    "hangles2="<<hangles[2]<<                    
	    "xyz0="<<xyz0[2]<<
	    "xyzz1="<<xyz1[2]<<
	    "z0="<<z0[i0]<<
	    "z1="<<z0[i1]<<
	    "radius="<<radiusbest<<
	    "deltabest="<<deltabest<< 
	    "phase0="<<phase[ibest][0]<<
	    "phase1="<<phase[ibest][1]<<
	    "\n"; 	  	  
	}
      }
    }
  }
  //
  //  Finf kinks loop
  // 
  //
  for (Int_t i =0;i<nentries;i++){
    if (sign[i]==0) continue;
    AliTPCseed * track0 = (AliTPCseed*)array->At(i);
    if (track0==0) {
      AliInfo("seed==0");
      continue;
    }
    ntracks++;
    //
    Double_t cradius0 = 40*40;
    Double_t cradius1 = 270*270;
    Double_t cdist1=8.;
    Double_t cdist2=8.;
    Double_t cdist3=0.55; 
    AliHelix* hlxi = (AliHelix*)helixes[i];
    for (Int_t j =i+1;j<nentries;j++){
      nall++;
      if (sign[j]*sign[i]<1) continue;
      int ncltot = nclusters[i];
      ncltot += nclusters[j];
      if ( ncltot>200) continue;
      if ( ncltot<80) continue;
      if ( TMath::Abs(zm[i]-zm[j])>60.) continue;
      if ( TMath::Abs(fim[i]-fim[j])>0.6 && TMath::Abs(fim[i]-fim[j])<5.7 ) continue;
      //AliTPCseed * track1 = (AliTPCseed*)array->At(j);  Double_t phase[2][2],radius[2];    
      AliHelix* hlxj = (AliHelix*)helixes[j];
      Int_t npoints = hlxi->GetRPHIintersections(*hlxj, phase, radius,20);
      if (npoints<1) continue;
      // cuts on radius      
      if (npoints==1){
	if (radius[0]<cradius0||radius[0]>cradius1) continue;
      }
      else{
	if ( (radius[0]<cradius0||radius[0]>cradius1) && (radius[1]<cradius0||radius[1]>cradius1) ) continue;
      }
      //      
      Double_t delta1=10000,delta2=10000;
      // cuts on the intersection radius
      hlxi->LinearDCA(*hlxj,phase[0][0],phase[0][1],radius[0],delta1);
      if (radius[0]<20&&delta1<1) continue; //intersection at vertex
      if (radius[0]<10&&delta1<3) continue; //intersection at vertex
      if (npoints==2){ 
	hlxi->LinearDCA(*hlxj,phase[1][0],phase[1][1],radius[1],delta2);
	if (radius[1]<20&&delta2<1) continue;  //intersection at vertex
	if (radius[1]<10&&delta2<3) continue;  //intersection at vertex	
      }
      //
      Double_t distance1 = TMath::Min(delta1,delta2);
      if (distance1>cdist1) continue;  // cut on DCA linear approximation
      //
      npoints = hlxi->GetRPHIintersections(*hlxj, phase, radius,20);
      hlxi->ParabolicDCA(*hlxj,phase[0][0],phase[0][1],radius[0],delta1);
      if (radius[0]<20&&delta1<1) continue; //intersection at vertex
      if (radius[0]<10&&delta1<3) continue; //intersection at vertex
      //
      if (npoints==2){ 
	hlxi->ParabolicDCA(*hlxj,phase[1][0],phase[1][1],radius[1],delta2);	
	if (radius[1]<20&&delta2<1) continue;  //intersection at vertex
	if (radius[1]<10&&delta2<3) continue;  //intersection at vertex	
      }            
      distance1 = TMath::Min(delta1,delta2);
      Float_t rkink =0;
      if (delta1<delta2){
	rkink = TMath::Sqrt(radius[0]);
      }
      else{
	rkink = TMath::Sqrt(radius[1]);
      }
      if (distance1>cdist2) continue;
      //
      //
      AliTPCseed * track1 = (AliTPCseed*)array->At(j);
      //
      //
      Int_t row0 = GetRowNumber(rkink); 
      if (row0<10)  continue;
      if (row0>150) continue;
      //
      //
      Float_t dens00=-1,dens01=-1;
      Float_t dens10=-1,dens11=-1;
      //
      Int_t found,foundable;//,ishared;
      //RS Seed don't keep their cluster pointers, cache cluster usage stat. for fast evaluation
      // FillSeedClusterStatCache(track0); // RS: Use this slow method only if sharing stat. is needed and used
      //
      //GetCachedSeedClusterStatistic(0,row0-5, found, foundable,ishared,kFALSE); // RS make sure FillSeedClusterStatCache is called
      //RS track0->GetClusterStatistic(0,row0-5, found, foundable,ishared,kFALSE);
      track0->GetClusterStatistic(0,row0-5, found, foundable);
      if (foundable>5) dens00 = Float_t(found)/Float_t(foundable);
      //
      //GetCachedSeedClusterStatistic(row0+5,155, found, foundable,ishared,kFALSE); // RS make sure FillSeedClusterStatCache is called
      //RS track0->GetClusterStatistic(row0+5,155, found, foundable,ishared,kFALSE);
      track0->GetClusterStatistic(row0+5,155, found, foundable);
      if (foundable>5) dens01 = Float_t(found)/Float_t(foundable);
      //
      //
      //RS Seed don't keep their cluster pointers, cache cluster usage stat. for fast evaluation
      // FillSeedClusterStatCache(track1); // RS: Use this slow method only if sharing stat. is needed and used
      //
      //GetCachedSeedClusterStatistic(0,row0-5, found, foundable,ishared,kFALSE); // RS make sure FillSeedClusterStatCache is called
      //RS track1->GetClusterStatistic(0,row0-5, found, foundable,ishared,kFALSE);
      track1->GetClusterStatistic(0,row0-5, found, foundable);
      if (foundable>10) dens10 = Float_t(found)/Float_t(foundable);
      //
      //GetCachedSeedClusterStatistic(row0+5,155, found, foundable,ishared,kFALSE); // RS make sure FillSeedClusterStatCache is called
      //RS track1->GetClusterStatistic(row0+5,155, found, foundable,ishared,kFALSE);
      track1->GetClusterStatistic(row0+5,155, found, foundable);
      if (foundable>10) dens11 = Float_t(found)/Float_t(foundable);
      //
      if (dens00<dens10 && dens01<dens11) continue;
      if (dens00>dens10 && dens01>dens11) continue;
      if (TMath::Max(dens00,dens10)<0.1)  continue;
      if (TMath::Max(dens01,dens11)<0.3)  continue;
      //
      if (TMath::Min(dens00,dens10)>0.6)  continue;
      if (TMath::Min(dens01,dens11)>0.6)  continue;

      //
      AliTPCseed * ktrack0, *ktrack1;
      if (dens00>dens10){
	ktrack0 = track0;
	ktrack1 = track1;
      }
      else{
	ktrack0 = track1;
	ktrack1 = track0;
      }
      if (TMath::Abs(ktrack0->GetC())>5) continue; // cut on the curvature for mother particle
      AliExternalTrackParam paramm(*ktrack0);
      AliExternalTrackParam paramd(*ktrack1);
      if (row0>60&&ktrack1->GetReference().GetX()>90.)new (&paramd) AliExternalTrackParam(ktrack1->GetReference()); 
      //
      //
      kink->SetMother(paramm);
      kink->SetDaughter(paramd);
      kink->Update();

      Float_t x[3] = { static_cast<Float_t>(kink->GetPosition()[0]),static_cast<Float_t>(kink->GetPosition()[1]),static_cast<Float_t>(kink->GetPosition()[2])};
      Int_t index[4];
      fkParam->Transform0to1(x,index);
      fkParam->Transform1to2(x,index);
      row0 = GetRowNumber(x[0]); 

      if (kink->GetR()<100) continue;
      if (kink->GetR()>240) continue;
      if (kink->GetPosition()[2]/kink->GetR()>AliTPCReconstructor::GetCtgRange()) continue;  //out of fiducial volume
      if (kink->GetDistance()>cdist3) continue;
      Float_t dird = kink->GetDaughterP()[0]*kink->GetPosition()[0]+kink->GetDaughterP()[1]*kink->GetPosition()[1];  // rough direction estimate
      if (dird<0) continue;

      Float_t dirm = kink->GetMotherP()[0]*kink->GetPosition()[0]+kink->GetMotherP()[1]*kink->GetPosition()[1];  // rough direction estimate
      if (dirm<0) continue;
      Float_t mpt = TMath::Sqrt(kink->GetMotherP()[0]*kink->GetMotherP()[0]+kink->GetMotherP()[1]*kink->GetMotherP()[1]);
      if (mpt<0.2) continue;

      if (mpt<1){
	//for high momenta momentum not defined well in first iteration
	Double_t qt   =  TMath::Sin(kink->GetAngle(2))*ktrack1->GetP();
	if (qt>0.35) continue; 
      }
      
      kink->SetLabel(CookLabel(ktrack0,0.4,0,row0),0);
      kink->SetLabel(CookLabel(ktrack1,0.4,row0,kMaxRow),1);
      if (dens00>dens10){
	kink->SetTPCDensity(dens00,0,0);
	kink->SetTPCDensity(dens01,0,1);
	kink->SetTPCDensity(dens10,1,0);
	kink->SetTPCDensity(dens11,1,1);
	kink->SetIndex(i,0);
	kink->SetIndex(j,1);
      }
      else{
	kink->SetTPCDensity(dens10,0,0);
	kink->SetTPCDensity(dens11,0,1);
	kink->SetTPCDensity(dens00,1,0);
	kink->SetTPCDensity(dens01,1,1);
	kink->SetIndex(j,0);
	kink->SetIndex(i,1);
      }

      if (mpt<1||kink->GetAngle(2)>0.1){
	//	angle and densities  not defined yet
	if (kink->GetTPCDensityFactor()<0.8) continue;
	if ((2-kink->GetTPCDensityFactor())*kink->GetDistance() >0.25) continue;
	if (kink->GetAngle(2)*ktrack0->GetP()<0.003) continue; //too small angle
	if (kink->GetAngle(2)>0.2&&kink->GetTPCDensityFactor()<1.15) continue;
	if (kink->GetAngle(2)>0.2&&kink->GetTPCDensity(0,1)>0.05) continue;

	Float_t criticalangle = track0->GetSigmaSnp2()+track0->GetSigmaTgl2();
	criticalangle+= track1->GetSigmaSnp2()+track1->GetSigmaTgl2();
	criticalangle= 3*TMath::Sqrt(criticalangle);
	if (criticalangle>0.02) criticalangle=0.02;
	if (kink->GetAngle(2)<criticalangle) continue;
      }
      //
      Int_t drow = Int_t(2.+0.5/(0.05+kink->GetAngle(2)));  // overlap region defined
      Float_t shapesum =0;
      Float_t sum = 0;
      for ( Int_t row = row0-drow; row<row0+drow;row++){
	if (row<0) continue;
	if (row>155) continue;
	//RS if (ktrack0->GetClusterPointer(row)) {
	if (ktrack0->GetClusterIndex2(row)>=0) {
	  const AliTPCTrackerPoints::Point *point = ktrack0->GetTrackPoint(row);
	  shapesum+=point->GetSigmaY()+point->GetSigmaZ();
	  sum++;
	}
	//RS if (ktrack1->GetClusterPointer(row)){
	if (ktrack1->GetClusterIndex2(row)>=0) {
	  const AliTPCTrackerPoints::Point *point =ktrack1->GetTrackPoint(row);
	  shapesum+=point->GetSigmaY()+point->GetSigmaZ();
	  sum++;
	}	
      }
      if (sum<4){
	kink->SetShapeFactor(-1.);
      }
      else{
	kink->SetShapeFactor(shapesum/sum);
      }      
      //      esd->AddKink(kink);
      //
      //      kink->SetMother(paramm);
      //kink->SetDaughter(paramd);
      
      Double_t chi2P2 = paramm.GetParameter()[2]-paramd.GetParameter()[2];
      chi2P2*=chi2P2;
      chi2P2/=paramm.GetCovariance()[5]+paramd.GetCovariance()[5];
      Double_t chi2P3 = paramm.GetParameter()[3]-paramd.GetParameter()[3];
      chi2P3*=chi2P3;
      chi2P3/=paramm.GetCovariance()[9]+paramd.GetCovariance()[9];
      //
      if ((AliTPCReconstructor::StreamLevel()&kStreamFindKinks)>0) {   // flag: stream track infroamtion in the FindKinks method
	(*fDebugStreamer)<<"kinkLpt"<<
	  "chi2P2="<<chi2P2<<
	  "chi2P3="<<chi2P3<<
	  "p0.="<<&paramm<<
	  "p1.="<<&paramd<<
	  "k.="<<kink<<
	  "\n";
      }
      if ( chi2P2+chi2P3<AliTPCReconstructor::GetRecoParam()->GetKinkAngleCutChi2(0)){
	continue;
      }
      //
      kinks.AddLast(kink);
      kink = new AliKink;
      ncandidates++;
    }
  }
  //
  // sort the kinks according quality - and refit them towards vertex
  //
  Int_t       nkinks    = kinks.GetEntriesFast();
  Float_t    quality[nkinks];
  Int_t      indexes[nkinks];
  AliTPCseed *mothers[nkinks];
  AliTPCseed *daughters[nkinks];
  memset(mothers,0,nkinks*sizeof(AliTPCseed*));
  memset(daughters,0,nkinks*sizeof(AliTPCseed*));
  //
  //
  for (Int_t i=0;i<nkinks;i++){
    quality[i] =100000;
    AliKink *kinkl = (AliKink*)kinks.At(i);
    //
    // refit kinks towards vertex
    // 
    Int_t index0 = kinkl->GetIndex(0);
    Int_t index1 = kinkl->GetIndex(1);
    AliTPCseed * ktrack0 = (AliTPCseed*)array->At(index0);
    AliTPCseed * ktrack1 = (AliTPCseed*)array->At(index1);
    //
    Int_t sumn=ktrack0->GetNumberOfClusters()+ktrack1->GetNumberOfClusters();
    //
    // Refit Kink under if too small angle
    //
    if (kinkl->GetAngle(2)<0.05){
      //
      // RS: if possible, remove kink before reseeding
      if (kinkl->GetDistance()>0.5 || kinkl->GetR()<110 || kinkl->GetR()>240) {
	delete kinks.RemoveAt(i);
	continue;
      }
      //
      kinkl->SetTPCRow0(GetRowNumber(kinkl->GetR()));
      Int_t row0 = kinkl->GetTPCRow0();
      Int_t drow = Int_t(2.+0.5/(0.05+kinkl->GetAngle(2)));
      //
      Int_t last  = row0-drow;
      if (last<40) last=40;
      if (last<ktrack0->GetFirstPoint()+25) last = ktrack0->GetFirstPoint()+25;
      AliTPCseed* seed0 = ReSeed(ktrack0,last,kFALSE);
      //
      Int_t first = row0+drow;
      if (first>130) first=130;
      if (first>ktrack1->GetLastPoint()-25) first = TMath::Max(ktrack1->GetLastPoint()-25,30);
      AliTPCseed* seed1 = ReSeed(ktrack1,first,kTRUE);
      //
      if (seed0 && seed1) {
	kinkl->SetStatus(1,8);
	if (RefitKink(*seed0,*seed1,*kinkl)) kinkl->SetStatus(1,9);
	row0 = GetRowNumber(kinkl->GetR());
	sumn = seed0->GetNumberOfClusters()+seed1->GetNumberOfClusters();
	mothers[i]   = seed0;
	daughters[i] = seed1;
      }
      else {
	delete kinks.RemoveAt(i);
	if (seed0) MarkSeedFree( seed0 );
	if (seed1) MarkSeedFree( seed1 );
	continue;
      }
    }
    //
    if (kinkl) quality[i] = 160*((0.1+kinkl->GetDistance())*(2.-kinkl->GetTPCDensityFactor()))/(sumn+40.);  //the longest -clossest will win
  }
  TMath::Sort(nkinks,quality,indexes,kFALSE);
  //
  //remove double find kinks
  //
  for (Int_t ikink0=1;ikink0<nkinks;ikink0++){
    AliKink * kink0 = (AliKink*) kinks.At(indexes[ikink0]);
    if (!kink0) continue;
    //
    for (Int_t ikink1=0;ikink1<ikink0;ikink1++){ 
      kink0 = (AliKink*) kinks.At(indexes[ikink0]);
      if (!kink0) continue;
      AliKink * kink1 = (AliKink*) kinks.At(indexes[ikink1]);
      if (!kink1) continue;
      // if not close kink continue
      if (TMath::Abs(kink1->GetPosition()[2]-kink0->GetPosition()[2])>10) continue;
      if (TMath::Abs(kink1->GetPosition()[1]-kink0->GetPosition()[1])>10) continue;
      if (TMath::Abs(kink1->GetPosition()[0]-kink0->GetPosition()[0])>10) continue;
      //
      AliTPCseed &mother0   =  mothers[indexes[ikink0]]   ? *mothers[indexes[ikink0]]   : *((AliTPCseed*)array->At(kink0->GetIndex(0)));
      AliTPCseed &daughter0 =  daughters[indexes[ikink0]] ? *daughters[indexes[ikink0]] : *((AliTPCseed*)array->At(kink0->GetIndex(1)));
      AliTPCseed &mother1   =  mothers[indexes[ikink1]]   ? *mothers[indexes[ikink1]]   : *((AliTPCseed*)array->At(kink1->GetIndex(0)));
      AliTPCseed &daughter1 =  daughters[indexes[ikink1]] ? *daughters[indexes[ikink1]] : *((AliTPCseed*)array->At(kink1->GetIndex(1)));

      Int_t row0 = (kink0->GetTPCRow0()+kink1->GetTPCRow0())/2;
      //
      Int_t same  = 0;
      Int_t both  = 0;
      Int_t samem = 0;
      Int_t bothm = 0;
      Int_t samed = 0;
      Int_t bothd = 0;
      //
      for (Int_t i=0;i<row0;i++){
	if (mother0.GetClusterIndex(i)>0 && mother1.GetClusterIndex(i)>0){
	  both++;
	  bothm++;
	  if (mother0.GetClusterIndex(i)==mother1.GetClusterIndex(i)){
	    same++;
	    samem++;
	  }
	}
      }

      for (Int_t i=row0;i<158;i++){
	//if (daughter0.GetClusterIndex(i)>0 && daughter0.GetClusterIndex(i)>0){ // RS: Bug? 
	if (daughter0.GetClusterIndex(i)>0 && daughter1.GetClusterIndex(i)>0){
	  both++;
	  bothd++;
	  if (mother0.GetClusterIndex(i)==mother1.GetClusterIndex(i)){
	    same++;
	    samed++;
	  }
	}
      }
      Float_t ratio = Float_t(same+1)/Float_t(both+1);
      Float_t ratiom = Float_t(samem+1)/Float_t(bothm+1);
      Float_t ratiod = Float_t(samed+1)/Float_t(bothd+1);
      if (ratio>0.3 && ratiom>0.5 &&ratiod>0.5) {
	Int_t sum0 = mother0.GetNumberOfClusters()+daughter0.GetNumberOfClusters();
	Int_t sum1 = mother1.GetNumberOfClusters()+daughter1.GetNumberOfClusters();
	if (sum1>sum0){
	  shared[kink0->GetIndex(0)]= kTRUE;
	  shared[kink0->GetIndex(1)]= kTRUE;	  
	  delete kinks.RemoveAt(indexes[ikink0]);
	  break;
	}
	else{
	  shared[kink1->GetIndex(0)]= kTRUE;
	  shared[kink1->GetIndex(1)]= kTRUE;	  
	  delete kinks.RemoveAt(indexes[ikink1]);
	}
      }
    }
  }


  for (Int_t i=0;i<nkinks;i++){
    AliKink * kinkl = (AliKink*) kinks.At(indexes[i]);
    if (!kinkl) continue;
    kinkl->SetTPCRow0(GetRowNumber(kinkl->GetR()));
    Int_t index0 = kinkl->GetIndex(0);
    Int_t index1 = kinkl->GetIndex(1);
    if (circular[index0]||(circular[index1]&&kinkl->GetDistance()>0.2)) continue;
    kinkl->SetMultiple(usage[index0],0);
    kinkl->SetMultiple(usage[index1],1);
    if (kinkl->GetMultiple()[0]+kinkl->GetMultiple()[1]>2) continue;
    if (kinkl->GetMultiple()[0]+kinkl->GetMultiple()[1]>0 && quality[indexes[i]]>0.2) continue;
    if (kinkl->GetMultiple()[0]+kinkl->GetMultiple()[1]>0 && kinkl->GetDistance()>0.2) continue;
    if (circular[index0]||(circular[index1]&&kinkl->GetDistance()>0.1)) continue;

    AliTPCseed * ktrack0 = (AliTPCseed*)array->At(index0);
    AliTPCseed * ktrack1 = (AliTPCseed*)array->At(index1);
    if (!ktrack0 || !ktrack1) continue;
    Int_t index = esd->AddKink(kinkl);
    //
    //
    if ( ktrack0->GetKinkIndex(0)==0 && ktrack1->GetKinkIndex(0)==0) {  //best kink
      if ( (mothers[indexes[i]]) && daughters[indexes[i]] && 
	   mothers[indexes[i]]->GetNumberOfClusters()>20 &&  // where they reseeded?
	   daughters[indexes[i]]->GetNumberOfClusters()>20 && 
	  (mothers[indexes[i]]->GetNumberOfClusters()+daughters[indexes[i]]->GetNumberOfClusters())>100){
	*ktrack0 = *mothers[indexes[i]];
	*ktrack1 = *daughters[indexes[i]];
      }
    }
    //
    ktrack0->SetKinkIndex(usage[index0],-(index+1));
    ktrack1->SetKinkIndex(usage[index1], (index+1));
    usage[index0]++;
    usage[index1]++;
  }
  //
  // Remove tracks corresponding to shared kink's
  //
  for (Int_t i=0;i<nentries;i++){
    AliTPCseed * track0 = (AliTPCseed*)array->At(i);
    if (!track0) continue;
    if (track0->GetKinkIndex(0)!=0) continue;
    if (shared[i]) MarkSeedFree( array->RemoveAt(i) );
  }

  //
  //
  RemoveUsed2(array,0.5,0.4,30);
  UnsignClusters();
  // RS: the cluster pointers are not permanently attached to the seed during the tracking, need to attach temporarily   
  AliTPCclusterMI* seedClusters[kMaxRow] = {0};                                                                          
  //
  for (Int_t i=0;i<nentries;i++){
    AliTPCseed * track0 = (AliTPCseed*)array->At(i);
    if (!track0) continue;
    //RS: if needed, attach temporary cluster array    
    const AliTPCclusterMI** seedClustersSave = track0->GetClusters();
    if (!seedClustersSave) { //RS: temporary attach clusters
      for (int ir=kMaxRow;ir--;) {
	int idx = track0->GetClusterIndex2(ir);
	seedClusters[ir] = idx<0 ? 0 : GetClusterMI(idx);
      }
      track0->SetClustersArrayTMP(seedClusters);
    }
    track0->CookdEdx(0.02,0.6);
    track0->CookPID();
    if (!seedClustersSave) track0->SetClustersArrayTMP(0);
  }
  //
  // RS use stack allocation instead of the heap
  AliTPCseed mother,daughter;
  AliKink kinkl;
  //
  for (Int_t i=0;i<nentries;i++){
    AliTPCseed * track0 = (AliTPCseed*)array->At(i);
    if (!track0) continue;
    if (track0->Pt()<1.4) continue;
    //remove double high momenta tracks - overlapped with kink candidates
    Int_t ishared=0;
    Int_t all   =0;
    for (Int_t icl=track0->GetFirstPoint();icl<track0->GetLastPoint(); icl++){
      Int_t tpcindex = track0->GetClusterIndex2(icl);
      if (tpcindex<0) continue;
      AliTPCclusterMI *cl = GetClusterMI(tpcindex);
      all++;
      if (cl->IsUsed(10)) ishared++;
    }
    if (Float_t(ishared+1)/Float_t(all+1)>0.5) {  
      MarkSeedFree( array->RemoveAt(i) );
      continue;
    }
    //
    if (track0->GetKinkIndex(0)!=0) continue;
    if (track0->GetNumberOfClusters()<80) continue;

    // AliTPCseed *pmother = new AliTPCseed(); // RS use stack allocation
    // AliTPCseed *pdaughter = new AliTPCseed();
    // AliKink *pkink = new AliKink;
    //
    if (CheckKinkPoint(track0,mother,daughter, kinkl)){
      if (mother.GetNumberOfClusters()<30||daughter.GetNumberOfClusters()<20) {
	//	delete pmother; // RS used stack allocation
	//	delete pdaughter;
	//	delete pkink;
	continue;  //too short tracks
      }
      if (mother.Pt()<1.4) {
	//	delete pmother; // RS used stack allocation
	//	delete pdaughter;
	//	delete pkink;
	continue;
      }
      Int_t row0= kinkl.GetTPCRow0();
      if (kinkl.GetDistance()>0.5 || kinkl.GetR()<110. || kinkl.GetR()>240.) {
	//	delete pmother; // RS used stack allocation
	//	delete pdaughter;
	//	delete pkink;
	continue;
      }
      //
      Int_t index = esd->AddKink(&kinkl);      
      mother.SetKinkIndex(0,-(index+1));
      daughter.SetKinkIndex(0,index+1);
      if (mother.GetNumberOfClusters()>50) {
	MarkSeedFree( array->RemoveAt(i) );
	AliTPCseed* mtc = new( NextFreeSeed() ) AliTPCseed(mother);
	mtc->SetPoolID(fLastSeedID);
	array->AddAt(mtc,i);
      }
      else{
	AliTPCseed* mtc = new( NextFreeSeed() ) AliTPCseed(mother);
	mtc->SetPoolID(fLastSeedID);
	array->AddLast(mtc);
      }
      AliTPCseed* dtc = new( NextFreeSeed() ) AliTPCseed(daughter);
      dtc->SetPoolID(fLastSeedID);
      array->AddLast(dtc);      
      for (Int_t icl=0;icl<row0;icl++) {
	Int_t tpcindex= mother.GetClusterIndex2(icl);
	if (tpcindex<0) continue;
	AliTPCclusterMI *cl = GetClusterMI(tpcindex);
	if (cl) cl->Use(20);
      }
      //
      for (Int_t icl=row0;icl<158;icl++) {
	Int_t tpcindex= mother.GetClusterIndex2(icl);
	if (tpcindex<0) continue;
	AliTPCclusterMI *cl = GetClusterMI(tpcindex);
	if (cl) cl->Use(20);
      }
      //
    }
    //delete pmother; // RS used stack allocation
    //delete pdaughter;
    //delete pkink;
  }

  for (int i=nkinks;i--;) {
    if (mothers[i]) MarkSeedFree( mothers[i] );
    if (daughters[i]) MarkSeedFree( daughters[i] );
  }
  //  delete [] daughters;
  //  delete [] mothers;
  //
  // RS: most of heap array are converted to stack arrays
  //  delete [] dca;
  //  delete []circular;
  //  delete []shared;
  //  delete []quality;
  //  delete []indexes;
  //
  delete kink;
  //  delete[]fim;
  //  delete[] zm;
  //  delete[] z0;
  //  delete [] usage;
  //  delete[] alpha;
  //  delete[] nclusters;
  //  delete[] sign;
  // delete[] helixes;
  if (fHelixPool) fHelixPool->Clear();
  kinks.Delete();
  //delete kinks;

  AliInfo(Form("Ncandidates=\t%d\t%d\t%d\t%d\n",esd->GetNumberOfKinks(),ncandidates,ntracks,nall));
  timer.Print();
}

Int_t AliTPCtracker::RefitKink(AliTPCseed &mother, AliTPCseed &daughter, const AliESDkink &knk)
{
  //
  // refit kink towards to the vertex
  //
  //
  AliKink &kink=(AliKink &)knk;

  Int_t row0 = GetRowNumber(kink.GetR());
  FollowProlongation(mother,0);
  mother.Reset(kFALSE);
  //
  FollowProlongation(daughter,row0);
  daughter.Reset(kFALSE);
  FollowBackProlongation(daughter,158);
  daughter.Reset(kFALSE);
  Int_t first = TMath::Max(row0-20,30); 
  Int_t last  = TMath::Min(row0+20,140);
  //
  const Int_t kNdiv =5;
  AliTPCseed  param0[kNdiv];  // parameters along the track
  AliTPCseed  param1[kNdiv];  // parameters along the track
  AliKink     kinks[kNdiv];   // corresponding kink parameters
  //
  Int_t rows[kNdiv];
  for (Int_t irow=0; irow<kNdiv;irow++){
    rows[irow] = first +((last-first)*irow)/(kNdiv-1);
  }
  // store parameters along the track
  //
  for (Int_t irow=0;irow<kNdiv;irow++){
    FollowBackProlongation(mother, rows[irow]);
    FollowProlongation(daughter,rows[kNdiv-1-irow]);       
    param0[irow] = mother;
    param1[kNdiv-1-irow] = daughter;
  }
  //
  // define kinks 
  for (Int_t irow=0; irow<kNdiv-1;irow++){
    if (param0[irow].GetNumberOfClusters()<kNdiv||param1[irow].GetNumberOfClusters()<kNdiv) continue;
    kinks[irow].SetMother(param0[irow]);
    kinks[irow].SetDaughter(param1[irow]);
    kinks[irow].Update();
  }
  //
  // choose kink with best "quality"
  Int_t index =-1;
  Double_t mindist = 10000;
  for (Int_t irow=0;irow<kNdiv;irow++){
    if (param0[irow].GetNumberOfClusters()<20||param1[irow].GetNumberOfClusters()<20) continue;
    if (TMath::Abs(kinks[irow].GetR())>240.) continue;
    if (TMath::Abs(kinks[irow].GetR())<100.) continue;
    //
    Float_t normdist = TMath::Abs(param0[irow].GetX()-kinks[irow].GetR())*(0.1+kink.GetDistance());
    normdist/= (param0[irow].GetNumberOfClusters()+param1[irow].GetNumberOfClusters()+40.);
    if (normdist < mindist){
      mindist = normdist;
      index = irow;
    }
  }
  //
  if (index==-1) return 0;
  //
  //
  param0[index].Reset(kTRUE);
  FollowProlongation(param0[index],0);
  //
  mother = param0[index];
  daughter = param1[index];  // daughter in vertex
  //
  kink.SetMother(mother);
  kink.SetDaughter(daughter);
  kink.Update();
  kink.SetTPCRow0(GetRowNumber(kink.GetR()));
  kink.SetTPCncls(param0[index].GetNumberOfClusters(),0);
  kink.SetTPCncls(param1[index].GetNumberOfClusters(),1);
  kink.SetLabel(CookLabel(&mother,0.4, 0,kink.GetTPCRow0()),0);
  kink.SetLabel(CookLabel(&daughter,0.4, kink.GetTPCRow0(),kMaxRow),1);
  mother.SetLabel(kink.GetLabel(0));
  daughter.SetLabel(kink.GetLabel(1));

  return 1;
}


void AliTPCtracker::UpdateKinkQualityM(AliTPCseed * seed){
  //
  // update Kink quality information for mother after back propagation
  //
  if (seed->GetKinkIndex(0)>=0) return; 
  for (Int_t ikink=0;ikink<3;ikink++){
    Int_t index = seed->GetKinkIndex(ikink);
    if (index>=0) break;
    index = TMath::Abs(index)-1;
    AliESDkink * kink = fEvent->GetKink(index);
    kink->SetTPCDensity(-1,0,0);
    kink->SetTPCDensity(1,0,1);
    //
    Int_t row0 = kink->GetTPCRow0() - 2 - Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (row0<15) row0=15;
    //
    Int_t row1 = kink->GetTPCRow0() + 2 +  Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (row1>145) row1=145;
    //
    Int_t found,foundable;//,shared;
    //GetSeedClusterStatistic(seed,0,row0, found, foundable,shared,kFALSE); //RS: seeds don't keep their clusters
    //RS seed->GetClusterStatistic(0,row0, found, foundable,shared,kFALSE);
    seed->GetClusterStatistic(0,row0, found, foundable);
    if (foundable>5)   kink->SetTPCDensity(Float_t(found)/Float_t(foundable),0,0);
    //
    //GetSeedClusterStatistic(seed,row1,155, found, foundable,shared,kFALSE); //RS: seeds don't keep their clusters
    //RS seed->GetClusterStatistic(row1,155, found, foundable,shared,kFALSE);
    seed->GetClusterStatistic(row1,155, found, foundable);
    if (foundable>5)   kink->SetTPCDensity(Float_t(found)/Float_t(foundable),0,1);
  }
    
}

void AliTPCtracker::UpdateKinkQualityD(AliTPCseed * seed){
  //
  // update Kink quality information for daughter after refit
  //
  if (seed->GetKinkIndex(0)<=0) return; 
  for (Int_t ikink=0;ikink<3;ikink++){
    Int_t index = seed->GetKinkIndex(ikink);
    if (index<=0) break;
    index = TMath::Abs(index)-1;
    AliESDkink * kink = fEvent->GetKink(index);
    kink->SetTPCDensity(-1,1,0);
    kink->SetTPCDensity(-1,1,1);
    //
    Int_t row0 = kink->GetTPCRow0() -2 - Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (row0<15) row0=15;
    //
    Int_t row1 = kink->GetTPCRow0() +2 +  Int_t( 0.5/ (0.05+kink->GetAngle(2)));
    if (row1>145) row1=145;
    //
    Int_t found,foundable;//,shared;
    //GetSeedClusterStatistic(0,row0, found, foundable,shared,kFALS); //RS: seeds don't keep their clusters
    //seed->GetClusterStatistic(0,row0, found, foundable,shared,kFALSE);
    seed->GetClusterStatistic(0,row0, found, foundable);
    if (foundable>5)   kink->SetTPCDensity(Float_t(found)/Float_t(foundable),1,0);
    //
    //GetSeedClusterStatistic(row1,155, found, foundable,shared,kFALSE); //RS: seeds don't keep their clusters
    // seed->GetClusterStatistic(row1,155, found, foundable,shared,kFALSE);
    seed->GetClusterStatistic(row1,155, found, foundable);
    if (foundable>5)   kink->SetTPCDensity(Float_t(found)/Float_t(foundable),1,1);
  }
    
}


Int_t  AliTPCtracker::CheckKinkPoint(AliTPCseed*seed,AliTPCseed &mother, AliTPCseed &daughter, const AliESDkink &knk)
{
  //
  // check kink point for given track
  // if return value=0 kink point not found
  // otherwise seed0 correspond to mother particle
  //           seed1 correspond to daughter particle
  //           kink  parameter of kink point
  AliKink &kink=(AliKink &)knk;

  Int_t middlerow = (seed->GetFirstPoint()+seed->GetLastPoint())/2;
  Int_t first = seed->GetFirstPoint(); 
  Int_t last  = seed->GetLastPoint();
  if (last-first<20) return 0;          // shortest length - 2*30 = 60 pad-rows

  
  AliTPCseed *seed1 = ReSeed(seed,middlerow+20, kTRUE);  //middle of chamber
  if (!seed1) return 0;
  FollowProlongation(*seed1,seed->GetLastPoint()-20);
  seed1->Reset(kTRUE);
  FollowProlongation(*seed1,158);
  seed1->Reset(kTRUE);  
  last = seed1->GetLastPoint();
  //
  AliTPCseed *seed0 = new( NextFreeSeed() ) AliTPCseed(*seed);
  seed0->SetPoolID(fLastSeedID);
  seed0->Reset(kFALSE);
  seed0->Reset();
  //
  AliTPCseed  param0[20];  // parameters along the track
  AliTPCseed  param1[20];  // parameters along the track
  AliKink     kinks[20];   // corresponding kink parameters
  Int_t rows[20];
  for (Int_t irow=0; irow<20;irow++){
    rows[irow] = first +((last-first)*irow)/19;
  }
  // store parameters along the track
  //
  for (Int_t irow=0;irow<20;irow++){
    FollowBackProlongation(*seed0, rows[irow]);
    FollowProlongation(*seed1,rows[19-irow]);       
    param0[irow] = *seed0;
    param1[19-irow] = *seed1;
  }
  //
  // define kinks 
  for (Int_t irow=0; irow<19;irow++){
    kinks[irow].SetMother(param0[irow]);
    kinks[irow].SetDaughter(param1[irow]);
    kinks[irow].Update();
  }
  //
  // choose kink with biggest change of angle
  Int_t index =-1;
  Double_t maxchange= 0;
  for (Int_t irow=1;irow<19;irow++){
    if (TMath::Abs(kinks[irow].GetR())>240.) continue;
    if (TMath::Abs(kinks[irow].GetR())<110.) continue;
    Float_t quality = TMath::Abs(kinks[irow].GetAngle(2))/(3.+TMath::Abs(kinks[irow].GetR()-param0[irow].GetX()));
    if ( quality > maxchange){
      maxchange = quality;
      index = irow;
      //
    }
  }
  MarkSeedFree( seed0 );
  MarkSeedFree( seed1 );
  if (index<0) return 0;
  //
  Int_t row0    = GetRowNumber(kinks[index].GetR());   //row 0 estimate
  seed0 = new( NextFreeSeed() ) AliTPCseed(param0[index]);
  seed0->SetPoolID(fLastSeedID);
  seed1 = new( NextFreeSeed() ) AliTPCseed(param1[index]);
  seed1->SetPoolID(fLastSeedID);
  seed0->Reset(kFALSE);
  seed1->Reset(kFALSE);
  seed0->ResetCovariance(10.);
  seed1->ResetCovariance(10.);
  FollowProlongation(*seed0,0);
  FollowBackProlongation(*seed1,158);
  mother = *seed0; // backup mother at position 0
  seed0->Reset(kFALSE);  
  seed1->Reset(kFALSE);
  seed0->ResetCovariance(10.);
  seed1->ResetCovariance(10.);
  //
  first = TMath::Max(row0-20,0);
  last  = TMath::Min(row0+20,158);
  //
  for (Int_t irow=0; irow<20;irow++){
    rows[irow] = first +((last-first)*irow)/19;
  }
  // store parameters along the track
  //
  for (Int_t irow=0;irow<20;irow++){
    FollowBackProlongation(*seed0, rows[irow]);
    FollowProlongation(*seed1,rows[19-irow]);       
    param0[irow] = *seed0;
    param1[19-irow] = *seed1;
  }
  //
  // define kinks 
  for (Int_t irow=0; irow<19;irow++){
    kinks[irow].SetMother(param0[irow]);
    kinks[irow].SetDaughter(param1[irow]);
    //    param0[irow].Dump();
    //param1[irow].Dump();
    kinks[irow].Update();
  }
  //
  // choose kink with biggest change of angle
  index =-1;
  maxchange= 0;
  for (Int_t irow=0;irow<20;irow++){
    if (TMath::Abs(kinks[irow].GetR())>250.) continue;
    if (TMath::Abs(kinks[irow].GetR())<90.) continue;
    Float_t quality = TMath::Abs(kinks[irow].GetAngle(2))/(3.+TMath::Abs(kinks[irow].GetR()-param0[irow].GetX()));
    if ( quality > maxchange){
      maxchange = quality;
      index = irow;
      //
    }
  }
  //
  //
  if (index==-1 || param0[index].GetNumberOfClusters()+param1[index].GetNumberOfClusters()<100){
    MarkSeedFree( seed0 );
    MarkSeedFree( seed1 );
    return 0;
  }

  //  Float_t anglesigma = TMath::Sqrt(param0[index].fC22+param0[index].fC33+param1[index].fC22+param1[index].fC33);
  
  kink.SetMother(param0[index]);
  kink.SetDaughter(param1[index]);
  kink.Update();

  Double_t chi2P2 = param0[index].GetParameter()[2]-param1[index].GetParameter()[2];
  chi2P2*=chi2P2;
  chi2P2/=param0[index].GetCovariance()[5]+param1[index].GetCovariance()[5];
  Double_t chi2P3 = param0[index].GetParameter()[3]-param1[index].GetParameter()[3];
  chi2P3*=chi2P3;
  chi2P3/=param0[index].GetCovariance()[9]+param1[index].GetCovariance()[9];
  //
  if (AliTPCReconstructor::StreamLevel()&kStreamFindKinks) { // flag: stream track infroamtion in the FindKinks method
    (*fDebugStreamer)<<"kinkHpt"<<
      "chi2P2="<<chi2P2<<
      "chi2P3="<<chi2P3<<
      "p0.="<<&param0[index]<<
      "p1.="<<&param1[index]<<
      "k.="<<&kink<<
      "\n";
  }
  if ( chi2P2+chi2P3<AliTPCReconstructor::GetRecoParam()->GetKinkAngleCutChi2(0)){
    MarkSeedFree( seed0 );
    MarkSeedFree( seed1 );
    return 0; 
  }


  row0    = GetRowNumber(kink.GetR());   
  kink.SetTPCRow0(row0);
  kink.SetLabel(CookLabel(seed0,0.5,0,row0),0);
  kink.SetLabel(CookLabel(seed1,0.5,row0,158),1);
  kink.SetIndex(-10,0);
  kink.SetIndex(int(param0[index].GetNumberOfClusters()+param1[index].GetNumberOfClusters()),1);
  kink.SetTPCncls(param0[index].GetNumberOfClusters(),0);
  kink.SetTPCncls(param1[index].GetNumberOfClusters(),1);
  //
  //
  //  new (&mother) AliTPCseed(param0[index]);
  daughter = param1[index];
  daughter.SetLabel(kink.GetLabel(1));  
  param0[index].Reset(kTRUE);
  FollowProlongation(param0[index],0);    
  mother = param0[index];
  mother.SetLabel(kink.GetLabel(0));
  if ( chi2P2+chi2P3<AliTPCReconstructor::GetRecoParam()->GetKinkAngleCutChi2(1)){
    mother=*seed;
  }
  MarkSeedFree( seed0 );
  MarkSeedFree( seed1 );
  //
  return 1;
}




AliTPCseed*  AliTPCtracker::ReSeed(AliTPCseed *t)
{
  //
  // reseed - refit -  track
  //
  Int_t first = 0;
  //  Int_t last  = fSectors->GetNRows()-1;
  //
  if (fSectors == fOuterSec){
    first = TMath::Max(first, t->GetFirstPoint()-fInnerSec->GetNRows());
    //last  = 
  }
  else
    first = t->GetFirstPoint();
  //
  AliTPCseed * seed = MakeSeed(t,0.1,0.5,0.9);
  FollowBackProlongation(*t,fSectors->GetNRows()-1);
  t->Reset(kFALSE);
  FollowProlongation(*t,first);
  return seed;
}







//_____________________________________________________________________________
Int_t AliTPCtracker::ReadSeeds(const TFile *inp) {
  //-----------------------------------------------------------------
  // This function reades track seeds.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 

  TFile *in=(TFile*)inp;
  if (!in->IsOpen()) {
     cerr<<"AliTPCtracker::ReadSeeds(): input file is not open !\n";
     return 1;
  }

  in->cd();
  TTree *seedTree=(TTree*)in->Get("Seeds");
  if (!seedTree) {
     cerr<<"AliTPCtracker::ReadSeeds(): ";
     cerr<<"can't get a tree with track seeds !\n";
     return 2;
  }
  AliTPCtrack *seed=new AliTPCtrack; 
  seedTree->SetBranchAddress("tracks",&seed);
  
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     AliTPCseed* sdc = new( NextFreeSeed() ) AliTPCseed(*seed/*,seed->GetAlpha()*/);
     sdc->SetPoolID(fLastSeedID);
     fSeeds->AddLast(sdc);
  }
  
  delete seed; // RS: this seed is not from the pool, delete it !!!
  delete seedTree; 
  savedir->cd();
  return 0;
}

Int_t AliTPCtracker::Clusters2TracksHLT (AliESDEvent *const esd, const AliESDEvent *hltEvent)
{
  //
  // clusters to tracks
  if (fSeeds) DeleteSeeds();
  else ResetSeedsPool();
  fEvent = esd; 
  fEventHLT = hltEvent;
  fAccountDistortions = AliTPCReconstructor::GetRecoParam()->GetAccountDistortions();
  if (AliTPCReconstructor::GetRecoParam()->GetUseOulierClusterFilter()) FilterOutlierClusters();  
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;  
  transform->SetCurrentRecoParam((AliTPCRecoParam*)AliTPCReconstructor::GetRecoParam());
  transform->SetCurrentTimeStamp( GetTimeStamp());
  transform->SetCurrentRun( GetRunNumber());

  //transform->SetCurrentTimeStamp( esd->GetTimeStamp());
  //transform->SetCurrentRun(esd->GetRunNumber());
  //
  if (AliTPCReconstructor::GetExtendedRoads()){
    fClExtraRoadY = AliTPCReconstructor::GetExtendedRoads()[0];
    fClExtraRoadZ = AliTPCReconstructor::GetExtendedRoads()[1];
    AliInfoF("Additional errors for roads: Y:%f Z:%f",fClExtraRoadY,fClExtraRoadZ);
  }
  const Double_t *errCluster = (AliTPCReconstructor::GetSystematicErrorCluster()) ?  
    AliTPCReconstructor::GetSystematicErrorCluster() : 
    AliTPCReconstructor::GetRecoParam()->GetSystematicErrorCluster();
  //
  fExtraClErrY2 = errCluster[0]*errCluster[0];
  fExtraClErrZ2 = errCluster[1]*errCluster[1];
  fExtraClErrYZ2 = fExtraClErrY2 + fExtraClErrZ2;
  AliInfoF("Additional errors for clusters: Y:%f Z:%f",errCluster[0],errCluster[1]);
  //
  if (AliTPCReconstructor::GetPrimaryDCACut()) {
    fPrimaryDCAYCut = AliTPCReconstructor::GetPrimaryDCACut()[0];
    fPrimaryDCAZCut = AliTPCReconstructor::GetPrimaryDCACut()[1];
    fDisableSecondaries = kTRUE;
    AliInfoF("Only primaries will be tracked with DCAY=%f and DCAZ=%f cuts",fPrimaryDCAYCut,fPrimaryDCAZCut);
  }
  //
  Clusters2Tracks();
  fEventHLT = 0;
  if (!fSeeds) return 1;
  FillESD(fSeeds);
  if ((AliTPCReconstructor::StreamLevel()&kStreamClDump)>0)  DumpClusters(0,fSeeds);
  return 0;
  //
}

Int_t AliTPCtracker::Clusters2Tracks(AliESDEvent *const esd)
{
  //
  // clusters to tracks
  return Clusters2TracksHLT( esd, 0);
}

//_____________________________________________________________________________
Int_t AliTPCtracker::Clusters2Tracks() {
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 
  TStopwatch timer;

  fIteration = 0;
  fSeeds = Tracking();

  if (fDebug>0){
    Info("Clusters2Tracks","Time for tracking: \t");timer.Print();timer.Start();
  }
  //activate again some tracks
  for (Int_t i=0; i<fSeeds->GetEntriesFast(); i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<20) {
      MarkSeedFree( fSeeds->RemoveAt(i) );
      continue;
    } 
    CookLabel(pt,0.1);
    if (pt->GetRemoval()==10) {
      if (pt->GetDensityFirst(20)>0.8 || pt->GetDensityFirst(30)>0.8 || pt->GetDensityFirst(40)>0.7)
	pt->Desactivate(10);  // make track again active  // MvL: should be 0 ?
      else{
	pt->Desactivate(20); 	
	MarkSeedFree( fSeeds->RemoveAt(i) );
      }
    } 
  }
  //
  RemoveUsed2(fSeeds,0.85,0.85,0);
  if (AliTPCReconstructor::GetRecoParam()->GetDoKinks()) FindKinks(fSeeds,fEvent);
  //FindCurling(fSeeds, fEvent,0);  
  if (AliTPCReconstructor::StreamLevel()&kStreamFindMultiMC)  FindMultiMC(fSeeds, fEvent,-1); // find multi found tracks
  RemoveUsed2(fSeeds,0.5,0.4,20);
  FindSplitted(fSeeds, fEvent,0); // find multi found tracks
  if (AliTPCReconstructor::StreamLevel()&kStreamFindMultiMC)  FindMultiMC(fSeeds, fEvent,0); // find multi found tracks

 //  //
//   // refit short tracks
//   //
  Int_t nseed=fSeeds->GetEntriesFast();
  //
  Int_t found = 0;
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;    
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<15) {
      MarkSeedFree( fSeeds->RemoveAt(i) );
      continue;
    }
    CookLabel(pt,0.1); //For comparison only
    //if ((pt->IsActive() || (pt->fRemoval==10) )&& nc>50 &&pt->GetNumberOfClusters()>0.4*pt->fNFoundable){
    if ((pt->IsActive() || (pt->GetRemoval()==10) )){
      found++;      
      if (fDebug>0) cerr<<found<<'\r';      
      pt->SetLab2(i);
    }
    else
      MarkSeedFree( fSeeds->RemoveAt(i) );
  }

  
  //RemoveOverlap(fSeeds,0.99,7,kTRUE);  
  SignShared(fSeeds);  
  //RemoveUsed(fSeeds,0.9,0.9,6);
  // 
  nseed=fSeeds->GetEntriesFast();
  found = 0;
  // RS: the cluster pointers are not permanently attached to the seed during the tracking, need to attach temporarily
  AliTPCclusterMI* seedClusters[kMaxRow];
  //
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t = *pt;
    if (!pt) continue;    
    Int_t nc=t.GetNumberOfClusters();
    if (nc<15) {
      MarkSeedFree( fSeeds->RemoveAt(i) );
      continue;
    }
    t.SetUniqueID(i);
    
    //    CheckKinkPoint(&t,0.05);
    if ((pt->IsActive() || (pt->GetRemoval()==10) )){
      found++;
      if (fDebug>0) cerr<<found<<'\r';
      pt->SetLab2(i);
    }
    else {MarkSeedFree( fSeeds->RemoveAt(i) );}
    //
    // temporarilly attach clusters and cook dedx
    const AliTPCclusterMI** seedClustersSave = t.GetClusters();
    if (!seedClustersSave) { //RS: temporary attach clusters
      for (int ir=kMaxRow;ir--;) {
	int idx = t.GetClusterIndex2(ir);
	seedClusters[ir] = idx<0 ? 0 : GetClusterMI(idx);
      }
      t.SetClustersArrayTMP(seedClusters);
    }
    t.CookdEdx(0.02,0.6);
    if (!seedClustersSave) t.SetClustersArrayTMP(0);
    //
  }

  SortTracks(fSeeds, 1);
   
  //  fNTracks = found;
  if (fDebug>0){
    Info("Clusters2Tracks","Time for overlap removal, track writing and dedx cooking: \t"); timer.Print();timer.Start();
  }
  //
  //  cerr<<"Number of found tracks : "<<"\t"<<found<<endl;  
  Info("Clusters2Tracks","Number of found tracks %d",found);  
  savedir->cd();
  //  UnloadClusters();
  //  
  return 0;
}

void AliTPCtracker::Tracking(TObjArray * arr)
{
  //
  // tracking of the seeds
  //
  fSectors = fOuterSec;
  ParallelTracking(arr,150,63);
  fSectors = fOuterSec;
  ParallelTracking(arr,63,0);

}

TObjArray * AliTPCtracker::Tracking(Int_t seedtype, Int_t i1, Int_t i2, Float_t cuts[4], Float_t dy, Int_t dsec)
{
  //
  //
  //tracking routine
  static TObjArray arrTracks;
  TObjArray * arr = &arrTracks;
  // 
  //  AliInfoF(">> %d | s:%d i1:%d i2:%d dy:%f dsec:%d",arr->GetEntriesFast(),seedtype,i1,i2,dy,dsec); // RSTMP


  fSectors = fOuterSec;
  TStopwatch timer;
  timer.Start();
  for (Int_t sec=0;sec<fkNOS;sec++){
    if (fAccountDistortions) {
      if (seedtype==3) MakeSeeds3Dist(arr,sec,i1,i2,cuts,dy, dsec);		     
      if (seedtype==4) MakeSeeds5Dist(arr,sec,i1,i2,cuts,dy);    
      if (seedtype==2) MakeSeeds2Dist(arr,sec,i1,i2,cuts,dy); //RS
    }
    else {
      if (seedtype==3) MakeSeeds3(arr,sec,i1,i2,cuts,dy, dsec);		     
      if (seedtype==4) MakeSeeds5(arr,sec,i1,i2,cuts,dy);    
      if (seedtype==2) MakeSeeds2(arr,sec,i1,i2,cuts,dy); //RS
    }
    //
  }
  if (fDebug>0){
    Info("Tracking","\nSeeding - %d\t%d\t%d\t%d\n",seedtype,i1,i2,arr->GetEntriesFast());
    timer.Print();
    timer.Start();
  }
  Tracking(arr);  
  if (fDebug>0){
    timer.Print();
  }
  //  AliInfoF("<< %d | s:%d i1:%d i2:%d dy:%f dsec:%d",arr->GetEntriesFast(),seedtype,i1,i2,dy,dsec); // RSTMP
  return arr;
}

TObjArray * AliTPCtracker::Tracking()
{
  // tracking
  //
  if (AliTPCReconstructor::GetRecoParam()->GetSpecialSeeding()) return TrackingSpecial();
  TStopwatch timer;
  timer.Start();
  Int_t nup=fOuterSec->GetNRows()+fInnerSec->GetNRows();

  TObjArray * seeds = fSeeds;
  if (!seeds) seeds = new TObjArray;
  else seeds->Clear();
  TObjArray * arr=0;
  Int_t fLastSeedRowSec=AliTPCReconstructor::GetRecoParam()->GetLastSeedRowSec();
  Int_t gapPrim = AliTPCReconstructor::GetRecoParam()->GetSeedGapPrim();
  Int_t gapSec = AliTPCReconstructor::GetRecoParam()->GetSeedGapSec();
  
  Int_t gap =20;
  Float_t cuts[4];
  cuts[0] = 0.002;
  cuts[1] = 1.5;
  cuts[2] = 3.;
  cuts[3] = 3.;
  Float_t fnumber  = 3.0;
  Float_t fdensity = 3.0;

  // make HLT seeds
  if (AliTPCReconstructor::GetRecoParam()->GetUseHLTPreSeeding()) {
    arr = MakeSeedsHLT( fEventHLT );
    if( arr ){
      SumTracks(seeds,arr);     
      delete arr;
      arr=0;
      //cout<<"HLT tracks left after sorting: "<<seeds->GetEntriesFast()<<endl;
      //SignClusters(seeds,fnumber,fdensity);    
    }
  }
  
  //  
  //find primaries  
  cuts[0]=0.0066;
  for (Int_t delta = 0; delta<18; delta+=gapPrim){
    //
    cuts[0]=0.0070;
    cuts[1] = 1.5;
    arr = Tracking(3,nup-1-delta,nup-1-delta-gap,cuts,-1,1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity); 
    //
    for (Int_t i=2;i<6;i+=2){
      // seed high pt tracks
      cuts[0]=0.0022;
      cuts[1]=0.3;
      arr = Tracking(3,nup-i-delta,nup-i-delta-gap,cuts,-1,0);
      SumTracks(seeds,arr);   
      SignClusters(seeds,fnumber,fdensity);        
    }
  }
  fnumber  = 4;
  fdensity = 4.;
  //  RemoveUsed(seeds,0.9,0.9,1);
  //  UnsignClusters();
  //  SignClusters(seeds,fnumber,fdensity);    

  //find primaries  
  cuts[0]=0.0077;
  for (Int_t delta = 20; delta<120; delta+=gapPrim){
    //
    // seed high pt tracks
    cuts[0]=0.0060;
    cuts[1]=0.3;
    cuts[2]=6.;
    arr = Tracking(3,nup-delta,nup-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);            

    cuts[0]=0.003;
    cuts[1]=0.3;
    cuts[2]=6.;
    arr = Tracking(3,nup-delta-5,nup-delta-5-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);            
  }

  cuts[0] = 0.01;
  cuts[1] = 2.0;
  cuts[2] = 3.;
  cuts[3] = 2.0;
  fnumber  = 2.;
  fdensity = 2.;
  
  if (fDebug>0){
    Info("Tracking()","\n\nPrimary seeding\t%d\n\n",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }
  //  RemoveUsed(seeds,0.75,0.75,1);
  //UnsignClusters();
  //SignClusters(seeds,fnumber,fdensity);
  
  if (fDisableSecondaries) return seeds;

  // find secondaries

  cuts[0] = 0.3;
  cuts[1] = 1.5;
  cuts[2] = 3.;
  cuts[3] = 1.5;

  arr = Tracking(4,nup-1,nup-1-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-2,nup-2-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-3,nup-3-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-5,nup-5-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  arr = Tracking(4,nup-7,nup-7-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //
  //
  arr = Tracking(4,nup-9,nup-9-gap,cuts,-1);
  SumTracks(seeds,arr);   
  SignClusters(seeds,fnumber,fdensity);   
  //


  for (Int_t delta = 9; delta<30; delta+=gapSec){
    //
    cuts[0] = 0.3;
    cuts[1] = 1.5;
    cuts[2] = 3.;
    cuts[3] = 1.5;
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
    //
    arr = Tracking(4,nup-3-delta,nup-5-delta-gap,cuts,4);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity); 
    //
  } 
  fnumber  = 1;
  fdensity = 1;
  //
  // change cuts
  fnumber  = 2.;
  fdensity = 2.;
  cuts[0]=0.0080;


  // find secondaries
  for (Int_t delta = 30; delta<fLastSeedRowSec; delta+=gapSec){
    //
    cuts[0] = 0.3;
    cuts[1] = 3.5;
    cuts[2] = 3.;
    cuts[3] = 3.5;
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
    //
    arr = Tracking(4,nup-5-delta,nup-5-delta-gap,cuts,5 );
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
  }
 
  if (fDebug>0){
    Info("Tracking()","\n\nSecondary seeding\t%d\n\n",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }

  return seeds;
  //
      
}


TObjArray * AliTPCtracker::TrackingSpecial()
{
  //
  // seeding adjusted for laser and cosmic tests - short tracks with big inclination angle
  // no primary vertex seeding tried
  //
  TStopwatch timer;
  timer.Start();
  Int_t nup=fOuterSec->GetNRows()+fInnerSec->GetNRows();

  TObjArray * seeds = fSeeds;
  if (!seeds) seeds = new TObjArray;
  else seeds->Clear();
  TObjArray * arr=0;
  
  Int_t   gap  = 15;
  Float_t cuts[4];
  Float_t fnumber  = 3.0;
  Float_t fdensity = 3.0;
  
  // find secondaries
  cuts[0] = AliTPCReconstructor::GetRecoParam()->GetMaxC();   // max curvature
  cuts[1] = 3.5;    // max tan(phi) angle for seeding
  cuts[2] = 3.;     // not used (cut on z primary vertex)     
  cuts[3] = 3.5;    // max tan(theta) angle for seeding

  for (Int_t delta = 0; nup-delta-gap-1>0; delta+=3){
    //
    arr = Tracking(4,nup-1-delta,nup-1-delta-gap,cuts,-1);
    SumTracks(seeds,arr);   
    SignClusters(seeds,fnumber,fdensity);   
  } 
 
  if (fDebug>0){
    Info("Tracking()","\n\nSecondary seeding\t%d\n\n",seeds->GetEntriesFast());
    timer.Print();
    timer.Start();
  }

  return seeds;
  //
      
}


void AliTPCtracker::SumTracks(TObjArray *arr1,TObjArray *&arr2)
{
  //
  //sum tracks to common container
  //remove suspicious tracks
  // RS: Attention: supplied tracks come in the static array, don't delete them
  Int_t nseed = arr2->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt=(AliTPCseed*)arr2->UncheckedAt(i);    
    if (pt){
      //
      // remove tracks with too big curvature
      //
      if (TMath::Abs(pt->GetC())>AliTPCReconstructor::GetRecoParam()->GetMaxC()){
	MarkSeedFree( arr2->RemoveAt(i) );
	continue;
      }
       // REMOVE VERY SHORT  TRACKS
      if (pt->GetNumberOfClusters()<20){ 
	MarkSeedFree( arr2->RemoveAt(i) );
	continue;
      }// patch 28 fev06
      // NORMAL ACTIVE TRACK
      if (pt->IsActive()){
	arr1->AddLast(arr2->RemoveAt(i));
	continue;
      }
      //remove not usable tracks
      if (pt->GetRemoval()!=10){
	MarkSeedFree( arr2->RemoveAt(i) );
	continue;
      }
     
      // ENABLE ONLY ENOUGH GOOD STOPPED TRACKS
      if (pt->GetDensityFirst(20)>0.8 || pt->GetDensityFirst(30)>0.8 || pt->GetDensityFirst(40)>0.7)
	arr1->AddLast(arr2->RemoveAt(i));
      else{      
	MarkSeedFree( arr2->RemoveAt(i) );
      }
    }
  }
  // delete arr2;  arr2 = 0; // RS: this is static array, don't delete it
}



void  AliTPCtracker::ParallelTracking(TObjArray *const arr, Int_t rfirst, Int_t rlast)
{
  //
  // try to track in parralel

  Int_t nseed=arr->GetEntriesFast();
  //prepare seeds for tracking
  for (Int_t i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i), &t=*pt; 
    if (!pt) continue;
    if (!t.IsActive()) continue;
    // follow prolongation to the first layer
    if ( (fSectors ==fInnerSec) || (t.GetFirstPoint()-fkParam->GetNRowLow()>rfirst+1) )  
      FollowProlongation(t, rfirst+1);
  }


  //
  for (Int_t nr=rfirst; nr>=rlast; nr--){ 
    if (nr<fInnerSec->GetNRows()) 
      fSectors = fInnerSec;
    else
      fSectors = fOuterSec;
    // make indexes with the cluster tracks for given       

    // find nearest cluster
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i), &t=*pt;       
      if (!pt) continue;
      if (nr==80) pt->UpdateReference();
      if (!pt->IsActive()) continue;
      //      if ( (fSectors ==fOuterSec) && (pt->fFirstPoint-fkParam->GetNRowLow())<nr) continue;
      if (pt->GetRelativeSector()>17) {
	continue;
      }
      UpdateClusters(t,nr);
    }
    // prolonagate to the nearest cluster - if founded
    for (Int_t i=0; i<nseed; i++) {
      AliTPCseed *pt=(AliTPCseed*)arr->UncheckedAt(i); 
      if (!pt) continue;
      if (!pt->IsActive()) continue; 
      // if ((fSectors ==fOuterSec) && (pt->fFirstPoint-fkParam->GetNRowLow())<nr) continue;
      if (pt->GetRelativeSector()>17) {
	continue;
      }
      FollowToNextCluster(*pt,nr);
    }
  }    
}

void AliTPCtracker::PrepareForBackProlongation(const TObjArray *const arr,Float_t fac) const
{
  //
  //
  // if we use TPC track itself we have to "update" covariance
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) {
      pt->Modify(fac);
      //
      //rotate to current local system at first accepted  point    
      Int_t index  = pt->GetClusterIndex2(pt->GetFirstPoint()); 
      Int_t sec    = (index&0xff000000)>>24;
      sec = sec%18;
      Float_t angle1 = fInnerSec->GetAlpha()*sec+fInnerSec->GetAlphaShift();
      if (angle1>TMath::Pi()) 
	angle1-=2.*TMath::Pi();
      Float_t angle2 = pt->GetAlpha();
      
      if (TMath::Abs(angle1-angle2)>0.001){
	if (!pt->Rotate(angle1-angle2)) return;
	//angle2 = pt->GetAlpha();
	//pt->fRelativeSector = pt->GetAlpha()/fInnerSec->GetAlpha();
	//if (pt->GetAlpha()<0) 
	//  pt->fRelativeSector+=18;
	//sec = pt->fRelativeSector;
      }
	
    }
    
  }


}
void AliTPCtracker::PrepareForProlongation(TObjArray *const arr, Float_t fac) const
{
  //
  //
  // if we use TPC track itself we have to "update" covariance
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) {
      pt->Modify(fac);
      pt->SetFirstPoint(pt->GetLastPoint()); 
    }
    
  }


}

Int_t AliTPCtracker::PropagateBack(const TObjArray *const arr)
{
  //
  // make back propagation
  //
  Int_t nseed= arr->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt&& pt->GetKinkIndex(0)<=0) { 
      //AliTPCseed *pt2 = new AliTPCseed(*pt);
      fSectors = fInnerSec;
      //FollowBackProlongation(*pt,fInnerSec->GetNRows()-1);
      //fSectors = fOuterSec;
      FollowBackProlongation(*pt,fInnerSec->GetNRows()+fOuterSec->GetNRows()-1,1);     
      //if (pt->GetNumberOfClusters()<(pt->fEsd->GetTPCclusters(0)) ){
      //	Error("PropagateBack","Not prolonged track %d",pt->GetLabel());
      //	FollowBackProlongation(*pt2,fInnerSec->GetNRows()+fOuterSec->GetNRows()-1);
      //}
    }
    if (pt&& pt->GetKinkIndex(0)>0) {
      AliESDkink * kink = fEvent->GetKink(pt->GetKinkIndex(0)-1);
      pt->SetFirstPoint(kink->GetTPCRow0());
      fSectors = fInnerSec;
      FollowBackProlongation(*pt,fInnerSec->GetNRows()+fOuterSec->GetNRows()-1,1);  
    }
    CookLabel(pt,0.3);
  }
  return 0;
}


Int_t AliTPCtracker::PropagateForward2(const TObjArray *const arr)
{
  //
  // make forward propagation
  //
  Int_t nseed= arr->GetEntriesFast();
  //
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)arr->UncheckedAt(i);
    if (pt) { 
      FollowProlongation(*pt,0,1,1);
      CookLabel(pt,0.3);
    }
    
  }
 return 0;
}


Int_t AliTPCtracker::PropagateForward()
{
  //
  // propagate track forward
  //UnsignClusters();
  Int_t nseed = fSeeds->GetEntriesFast();
  for (Int_t i=0;i<nseed;i++){
    AliTPCseed *pt = (AliTPCseed*)fSeeds->UncheckedAt(i);
    if (pt){
      AliTPCseed &t = *pt;
      Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
      if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
      if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
      t.SetRelativeSector(Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN);
    }
  }
  
  fSectors = fOuterSec;
  ParallelTracking(fSeeds,fOuterSec->GetNRows()+fInnerSec->GetNRows()-1,fInnerSec->GetNRows());
  fSectors = fInnerSec;
  ParallelTracking(fSeeds,fInnerSec->GetNRows()-1,0);
  return 1;
}






Int_t AliTPCtracker::PropagateBack(AliTPCseed *const pt, Int_t row0, Int_t row1)
{
  //
  // make back propagation, in between row0 and row1
  //
  
  if (pt) { 
    fSectors = fInnerSec;
    Int_t  r1;
    //
    if (row1<fSectors->GetNRows()) 
      r1 = row1;
    else 
      r1 = fSectors->GetNRows()-1;

    if (row0<fSectors->GetNRows()&& r1>0 )
      FollowBackProlongation(*pt,r1);
    if (row1<=fSectors->GetNRows())
      return 0;
    //
    r1 = row1 - fSectors->GetNRows();
    if (r1<=0) return 0;
    if (r1>=fOuterSec->GetNRows()) return 0;
    fSectors = fOuterSec;
    return FollowBackProlongation(*pt,r1);
  }        
  return 0;
}




void  AliTPCtracker::GetShape(AliTPCseed * seed, Int_t row)
{
  // gets cluster shape
  // 
  AliTPCClusterParam * clparam = AliTPCcalibDB::Instance()->GetClusterParam();
  Float_t zdrift = TMath::Abs((fkParam->GetZLength(0)-TMath::Abs(seed->GetZ())));
  Int_t type = (seed->GetSector() < fkParam->GetNSector()/2) ? 0: (row>126) ? 1:2;
  Double_t angulary  = seed->GetSnp();

  if (TMath::Abs(angulary)>AliTPCReconstructor::GetMaxSnpTracker()) {
    angulary = TMath::Sign(AliTPCReconstructor::GetMaxSnpTracker(),angulary);
  }

  angulary = angulary*angulary/((1.-angulary)*(1.+angulary));
  Double_t angularz  = seed->GetTgl()*seed->GetTgl()*(1.+angulary);
  
  Double_t sigmay =  clparam->GetRMS0(0,type,zdrift,TMath::Sqrt(TMath::Abs(angulary)));
  Double_t sigmaz =  clparam->GetRMS0(1,type,zdrift,TMath::Sqrt(TMath::Abs(angularz)));
  seed->SetCurrentSigmaY2(sigmay*sigmay);
  seed->SetCurrentSigmaZ2(sigmaz*sigmaz);
  // Float_t sd2 = TMath::Abs((fkParam->GetZLength(0)-TMath::Abs(seed->GetZ())))*fkParam->GetDiffL()*fkParam->GetDiffL();
//   //  Float_t padlength =  fkParam->GetPadPitchLength(seed->fSector);
//   Float_t padlength =  GetPadPitchLength(row);
//   //
//   Float_t sresy = (seed->GetSector() < fkParam->GetNSector()/2) ? 0.2 :0.3;
//   seed->SetCurrentSigmaY2(sd2+padlength*padlength*angulary/12.+sresy*sresy);  
//   //
//   Float_t sresz = fkParam->GetZSigma();
//   seed->SetCurrentSigmaZ2(sd2+padlength*padlength*angularz*angularz*(1+angulary)/12.+sresz*sresz);
  /*
  Float_t wy = GetSigmaY(seed);
  Float_t wz = GetSigmaZ(seed);
  wy*=wy;
  wz*=wz;
  if (TMath::Abs(wy/seed->fCurrentSigmaY2-1)>0.0001 || TMath::Abs(wz/seed->fCurrentSigmaZ2-1)>0.0001 ){
    printf("problem\n");
  }
  */
}



//__________________________________________________________________________
void AliTPCtracker::CookLabel(AliKalmanTrack *tk, Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  //  AliTPCseed * t = dynamic_cast<AliTPCseed*>(tk);
  //  if(!t){
  //  printf("%s:%d wrong type \n",(char*)__FILE__,__LINE__);
  //  return;
  //  }
  AliTPCseed * t = (AliTPCseed*)tk; // RS avoid slow dynamic cast

  Int_t noc=t->GetNumberOfClusters();
  if (noc<10){
    //printf("\nnot founded prolongation\n\n\n");
    //t->Dump();
    return ;
  }
  Int_t lb[kMaxRow];
  Int_t mx[kMaxRow];
  AliTPCclusterMI *clusters[kMaxRow];
  //
  for (Int_t i=0;i<kMaxRow;i++) {
    clusters[i]=0;
    lb[i]=mx[i]=0;
  }

  Int_t i;
  Int_t current=0;
  for (i=0; i<kMaxRow && current<noc; i++) {
     
     Int_t index=t->GetClusterIndex2(i);
     if (index<=0) continue; 
     if (index&0x8000) continue;
     //     
     clusters[current++]=GetClusterMI(index);
     // if (t->GetClusterPointer(i)){
     //   clusters[current]=t->GetClusterPointer(i);     
     //   current++;
     // }
  }
  noc = current;

  Int_t lab=123456789;
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i];
    if (!c) continue;
    lab=TMath::Abs(c->GetLabel(0));
    Int_t j;
    for (j=0; j<noc; j++) if (lb[j]==lab || mx[j]==0) break;
    lb[j]=lab;
    (mx[j])++;
  }

  Int_t max=0;
  for (i=0; i<noc; i++) if (mx[i]>max) {max=mx[i]; lab=lb[i];}
    
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i]; 
    if (!c) continue;
    if (TMath::Abs(c->GetLabel(1)) == lab ||
        TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }
  if (noc<=0) { lab=-1; return;}
  if ((1.- Float_t(max)/(noc)) > wrong) lab=-lab;

  else {
     Int_t tail=Int_t(0.10*noc);
     max=0;
     Int_t ind=0;
     for (i=1; i<kMaxRow&&ind<tail; i++) {
       //       AliTPCclusterMI *c=clusters[noc-i];
       AliTPCclusterMI *c=clusters[i];
       if (!c) continue;
       if (lab == TMath::Abs(c->GetLabel(0)) ||
           lab == TMath::Abs(c->GetLabel(1)) ||
           lab == TMath::Abs(c->GetLabel(2))) max++;
       ind++;
     }
     if (max < Int_t(0.5*tail)) lab=-lab;
  }

  t->SetLabel(lab);

  //  delete[] lb;
  //delete[] mx;
  //delete[] clusters;
}


//__________________________________________________________________________
Int_t AliTPCtracker::CookLabel(AliTPCseed *const t, Float_t wrong,Int_t first, Int_t last) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  if (noc<10){
    //printf("\nnot founded prolongation\n\n\n");
    //t->Dump();
    return -1;
  }
  Int_t lb[kMaxRow];
  Int_t mx[kMaxRow];
  AliTPCclusterMI *clusters[kMaxRow];
  //
  for (Int_t i=0;i<kMaxRow;i++) {
    clusters[i]=0;
    lb[i]=mx[i]=0;
  }

  Int_t i;
  Int_t current=0;
  for (i=0; i<kMaxRow && current<noc; i++) {
    if (i<first) continue;
    if (i>last)  continue;
     Int_t index=t->GetClusterIndex2(i);
     if (index<=0) continue; 
     if (index&0x8000) continue;
     //     
     clusters[current++]=GetClusterMI(index);
     // if (t->GetClusterPointer(i)){
     //   clusters[current]=t->GetClusterPointer(i);     
     //   current++;
     // }
  }
  noc = current;
  //if (noc<5) return -1;
  Int_t lab=123456789;
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i];
    if (!c) continue;
    lab=TMath::Abs(c->GetLabel(0));
    Int_t j;
    for (j=0; j<noc; j++) if (lb[j]==lab || mx[j]==0) break;
    lb[j]=lab;
    (mx[j])++;
  }

  Int_t max=0;
  for (i=0; i<noc; i++) if (mx[i]>max) {max=mx[i]; lab=lb[i];}
    
  for (i=0; i<noc; i++) {
    AliTPCclusterMI *c=clusters[i]; 
    if (!c) continue;
    if (TMath::Abs(c->GetLabel(1)) == lab ||
        TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }
  if (noc<=0) { lab=-1; return -1;}
  if ((1.- Float_t(max)/(noc)) > wrong) lab=-lab;

  else {
     Int_t tail=Int_t(0.10*noc);
     max=0;
     Int_t ind=0;
     for (i=1; i<kMaxRow&&ind<tail; i++) {
       //       AliTPCclusterMI *c=clusters[noc-i];
       AliTPCclusterMI *c=clusters[i];
       if (!c) continue;
       if (lab == TMath::Abs(c->GetLabel(0)) ||
           lab == TMath::Abs(c->GetLabel(1)) ||
           lab == TMath::Abs(c->GetLabel(2))) max++;
       ind++;
     }
     if (max < Int_t(0.5*tail)) lab=-lab;
  }

  //  t->SetLabel(lab);
  return lab;
  //  delete[] lb;
  //delete[] mx;
  //delete[] clusters;
}


Int_t  AliTPCtracker::GetRowNumber(Double_t x[3]) const 
{
  //return pad row number for given x vector
  Float_t phi = TMath::ATan2(x[1],x[0]);
  if(phi<0) phi=2.*TMath::Pi()+phi;
  //  Get the local angle in the sector philoc
  const Float_t kRaddeg = 180/3.14159265358979312;
  Float_t phiangle   = (Int_t (phi*kRaddeg/20.) + 0.5)*20./kRaddeg;
  Double_t localx    = x[0]*TMath::Cos(phiangle)-x[1]*TMath::Sin(phiangle);
  return GetRowNumber(localx);
}



void AliTPCtracker::MakeESDBitmaps(AliTPCseed *t, AliESDtrack *esd)
{
  //-----------------------------------------------------------------------
  // Fill the cluster and sharing bitmaps of the track
  //-----------------------------------------------------------------------

  Int_t firstpoint = 0;
  Int_t lastpoint = kMaxRow;
  //  AliTPCclusterMI    *cluster;
  
  Int_t nclsf = 0;
  TBits clusterMap(kMaxRow);
  TBits sharedMap(kMaxRow);
  TBits fitMap(kMaxRow);
  for (int iter=firstpoint; iter<lastpoint; iter++) {
    // Change to cluster pointers to see if we have a cluster at given padrow
    Int_t tpcindex= t->GetClusterIndex2(iter);
    if (tpcindex>=0) {
      clusterMap.SetBitNumber(iter, kTRUE);
      if (t->IsShared(iter)) sharedMap.SetBitNumber(iter,kTRUE); // RS shared flag moved to seed
      //
      if ( (tpcindex&0x8000) == 0)  {
	fitMap.SetBitNumber(iter, kTRUE);
	nclsf++;
      }
    }
    // if (t->GetClusterIndex(iter) >= 0 && (t->GetClusterIndex(iter) & 0x8000) == 0) {
    //   fitMap.SetBitNumber(iter, kTRUE);
    //   nclsf++;
    // }
  }
  esd->SetTPCClusterMap(clusterMap);
  esd->SetTPCSharedMap(sharedMap);
  esd->SetTPCFitMap(fitMap);
  if (nclsf != t->GetNumberOfClusters())
    AliDebug(3,Form("Inconsistency between ncls %d and indices %d (found %d)",t->GetNumberOfClusters(),nclsf,esd->GetTPCClusterMap().CountBits()));
}

Bool_t AliTPCtracker::IsFindable(AliTPCseed & track){
  //
  // return flag if there is findable cluster at given position
  //
  Float_t kDeltaZ=10;
  Float_t z = track.GetZ();
  
  if (TMath::Abs(z)<(AliTPCReconstructor::GetCtgRange()*track.GetX()+kDeltaZ) && 
      TMath::Abs(z)<fkParam->GetZLength(0) && 
      (TMath::Abs(track.GetSnp())<AliTPCReconstructor::GetMaxSnpTracker()))
    return kTRUE;
  return kFALSE;      
}


void AliTPCtracker::AddCovariance(AliTPCseed * seed){
  //
  // Adding systematic error estimate to the covariance matrix
  //                !!!! the systematic error for element 4 is in 1/GeV 
  // 03.03.2012     MI changed in respect to the previous versions
  // in case systemtic errors defiend statically no correlation added (needed for CPass0)
  const Double_t *param = (AliTPCReconstructor::GetSystematicError()!=NULL) ? AliTPCReconstructor::GetSystematicError():AliTPCReconstructor::GetRecoParam()->GetSystematicError();
  //
  // use only the diagonal part if not specified otherwise
  if (!AliTPCReconstructor::GetRecoParam()->GetUseSystematicCorrelation() || AliTPCReconstructor::GetSystematicError()) return AddCovarianceAdd(seed);
  //
  Double_t *covarS= (Double_t*)seed->GetCovariance();
  Double_t factor[5]={1,1,1,1,1};
  factor[0]= TMath::Sqrt(TMath::Abs((covarS[0] + param[0]*param[0])/covarS[0]));
  factor[1]= TMath::Sqrt(TMath::Abs((covarS[2] + param[1]*param[1])/covarS[2]));
  factor[2]= TMath::Sqrt(TMath::Abs((covarS[5] + param[2]*param[2])/covarS[5]));
  factor[3]= TMath::Sqrt(TMath::Abs((covarS[9] + param[3]*param[3])/covarS[9]));
  factor[4]= TMath::Sqrt(TMath::Abs((covarS[14] +param[4]*param[4])/covarS[14]));
  //
  factor[0]=factor[2];
  factor[4]=factor[2];
  // 0
  // 1    2
  // 3    4    5
  // 6    7    8    9 
  // 10   11   12   13   14
  for (Int_t i=0; i<5; i++){
    for (Int_t j=i; j<5; j++){
      Int_t index=seed->GetIndex(i,j);
      covarS[index]*=factor[i]*factor[j];
    }
  }
}


void AliTPCtracker::AddCovarianceAdd(AliTPCseed * seed){
  //
  // Adding systematic error - as additive factor without correlation
  //
  //                !!!! the systematic error for element 4 is in 1/GeV 
  // 03.03.2012     MI changed in respect to the previous versions

  const Double_t *param = (AliTPCReconstructor::GetSystematicError()!=NULL) ? AliTPCReconstructor::GetSystematicError():AliTPCReconstructor::GetRecoParam()->GetSystematicError();
  Double_t *covarIn= (Double_t*)seed->GetCovariance();
  Double_t covar[15];
  for (Int_t i=0;i<15;i++) covar[i]=0;
  // 0
  // 1    2
  // 3    4    5
  // 6    7    8    9 
  // 10   11   12   13   14
  covar[0] = param[0]*param[0];
  covar[2] = param[1]*param[1];
  covar[5] = param[2]*param[2];
  covar[9] = param[3]*param[3];
  covar[14]= param[4]*param[4];
  //
  covar[1]=TMath::Sqrt((covar[0]*covar[2]))*covarIn[1]/TMath::Sqrt((covarIn[0]*covarIn[2]));
  //
  covar[3]=TMath::Sqrt((covar[0]*covar[5]))*covarIn[3]/TMath::Sqrt((covarIn[0]*covarIn[5]));
  covar[4]=TMath::Sqrt((covar[2]*covar[5]))*covarIn[4]/TMath::Sqrt((covarIn[2]*covarIn[5]));
  //
  covar[6]=TMath::Sqrt((covar[0]*covar[9]))*covarIn[6]/TMath::Sqrt((covarIn[0]*covarIn[9]));
  covar[7]=TMath::Sqrt((covar[2]*covar[9]))*covarIn[7]/TMath::Sqrt((covarIn[2]*covarIn[9]));
  covar[8]=TMath::Sqrt((covar[5]*covar[9]))*covarIn[8]/TMath::Sqrt((covarIn[5]*covarIn[9]));
  //
  covar[10]=TMath::Sqrt((covar[0]*covar[14]))*covarIn[10]/TMath::Sqrt((covarIn[0]*covarIn[14]));
  covar[11]=TMath::Sqrt((covar[2]*covar[14]))*covarIn[11]/TMath::Sqrt((covarIn[2]*covarIn[14]));
  covar[12]=TMath::Sqrt((covar[5]*covar[14]))*covarIn[12]/TMath::Sqrt((covarIn[5]*covarIn[14]));
  covar[13]=TMath::Sqrt((covar[9]*covar[14]))*covarIn[13]/TMath::Sqrt((covarIn[9]*covarIn[14]));
  //
  seed->AddCovariance(covar);
}

//_____________________________________________________________________________
Bool_t  AliTPCtracker::IsTPCHVDipEvent(AliESDEvent const *esdEvent)
{
  //
  // check events affected by TPC HV dip
  //
  if(!esdEvent) return kFALSE;

  // Init TPC OCDB
  AliTPCcalibDB *db=AliTPCcalibDB::Instance();
  if(!db) return kFALSE;
  db->SetRun(esdEvent->GetRunNumber());

  // maximum allowed voltage before an event is identified as a dip event
  // and scanning period
  const Double_t kTPCHVdip          = db->GetParameters()->GetMaxDipVoltage(); 
  const Double_t dipEventScanPeriod = db->GetParameters()->GetVoltageDipScanPeriod();
  const Double_t tevSec             = esdEvent->GetTimeStamp();
  
  for(Int_t sector=0; sector<72; sector++)
  {
    // don't use excluded chambers, since the state is not defined at all
    if (!db->GetChamberHVStatus(sector)) continue;
    
    // get hv sensor of the chamber
    AliDCSSensor *sensor = db->GetChamberHVSensor(sector);
    if (!sensor) continue;
    TGraph *grSensor=sensor->GetGraph();
    if (!grSensor) continue;
    if (grSensor->GetN()<1) continue;
    
    // get median
    const Double_t median = db->GetChamberHighVoltageMedian(sector);
    if(median < 1.) continue;
    
    for (Int_t ipoint=0; ipoint<grSensor->GetN()-1; ++ipoint){
      Double_t nextTime=grSensor->GetX()[ipoint+1]*3600+sensor->GetStartTime();
      if (tevSec-dipEventScanPeriod>nextTime) continue;
      const Float_t deltaV=TMath::Abs(grSensor->GetY()[ipoint]-median);
      if (deltaV>kTPCHVdip) {
        AliDebug(3,Form("HV dip detected in ROC '%02d' with dV '%.2f' at time stamp '%.0f'",sector,deltaV,tevSec));
        return kTRUE;
      }
      if (nextTime>tevSec+dipEventScanPeriod) break;
    }
  }
  
  return kFALSE;
}

//________________________________________
void AliTPCtracker::MarkSeedFree(TObject *sd) 
{
  // account that this seed is "deleted" 
  //  AliTPCseed* seed = dynamic_cast<AliTPCseed*>(sd);
  //  if (!seed) {
  //    AliError(Form("Freeing of non-AliTPCseed %p from the pool is requested",sd)); 
  //    return;
  //  }
  AliTPCseed* seed = (AliTPCseed*)sd;
  int id = seed->GetPoolID();
  if (id<0) {
    AliError(Form("Freeing of seed %p NOT from the pool is requested",sd)); 
    return;
  }
  //  AliInfo(Form("%d %p",id, seed));
  fSeedsPool->RemoveAt(id);
  if (fFreeSeedsID.GetSize()<=fNFreeSeeds) fFreeSeedsID.Set( 2*fNFreeSeeds + 100 );
  fFreeSeedsID.GetArray()[fNFreeSeeds++] = id;
}

//________________________________________
TObject *&AliTPCtracker::NextFreeSeed()
{
  // return next free slot where the seed can be created
  fLastSeedID = fNFreeSeeds ? fFreeSeedsID.GetArray()[--fNFreeSeeds] : fSeedsPool->GetEntriesFast();
  //  AliInfo(Form("%d",fLastSeedID));
  return (*fSeedsPool)[ fLastSeedID ];
  //
}

//________________________________________
void AliTPCtracker::ResetSeedsPool()
{
  // mark all seeds in the pool as unused
  AliInfo(Form("CurrentSize: %d, BookedUpTo: %d, free: %d",fSeedsPool->GetSize(),fSeedsPool->GetEntriesFast(),fNFreeSeeds));
  fNFreeSeeds = 0;
  fSeedsPool->Clear(); // RS: nominally the seeds may allocate memory...
  
}

Int_t  AliTPCtracker::PropagateToRowHLT(AliTPCseed *pt, int nrow)
{
  AliTPCseed &t=*pt;
  Double_t x= GetXrow(nrow);
  Double_t  ymax= GetMaxY(nrow);
  Int_t rotate = 0;
  Int_t nRotations=0;
  int ret = 1;
  do{
    rotate = 0;
    if (!t.PropagateTo(x) ){
      //cout<<"can't propagate to row "<<nrow<<", x="<<t.GetX()<<" -> "<<x<<endl; 
      //t.Print();
      ret = 0;
      break;
    }
    t.SetRow(nrow);
    Double_t y = t.GetY();
    double yEdgeDist = y;
    if  (fAccountDistortions) yEdgeDist -= GetYSectEdgeDist(t.GetRelativeSector(),nrow,y,t.GetZ());
    if( y>0 && yEdgeDist>ymax) {
      if( rotate!=-1 ) rotate=1;
    } else if  (y<0 && yEdgeDist<-ymax) {	
      if( rotate!=1 ) rotate = -1;
    }
    if( rotate==0 ) break;
    //cout<<"rotate at row "<<nrow<<": "<<rotate<<endl;
    if (!t.Rotate( rotate==1 ?fSectors->GetAlpha() :-fSectors->GetAlpha())) {
      //cout<<"can't rotate "<<endl; 
      ret=0;
      break;
    }
    nRotations+= rotate;
  }while(rotate!=0);
  if( nRotations!=0 ){
    int newSec= t.GetRelativeSector()+nRotations;
    if( newSec>=fN ) newSec-=fN;
    else if( newSec<0 ) newSec +=fN; 
    //cout<<"rotate at row "<<nrow<<": "<<nRotations<<" times "<<" sec "
    //<< t.GetRelativeSector()<<" -> "<<newSec<<endl;
    t.SetRelativeSector(newSec);
  }
  return ret;
}

void  AliTPCtracker::TrackFollowingHLT(TObjArray *const arr )
{
  //
  // try to track in parralel

  AliFatal("RS: This method is not yet aware of cluster pointers no present in the in-memory tracks");

  Int_t nRows=fOuterSec->GetNRows()+fInnerSec->GetNRows();
  fSectors=fInnerSec;

  Int_t nseed=arr->GetEntriesFast();
  //cout<<"Parallel tracking My.."<<endl;
  double shapeY2[kMaxRow], shapeZ2[kMaxRow];
  Int_t clusterIndex[kMaxRow];
 
  for (Int_t iSeed=0; iSeed<nseed; iSeed++) {
    //if( iSeed!=1 ) continue;
    AliTPCseed *pt=(AliTPCseed*) (arr->UncheckedAt(iSeed));
    if (!pt) continue;
    AliTPCseed &t=*pt;    
    
    //cout <<"Pt "<<t.GetSigned1Pt()<<endl;
    
    // t.Print();
    
    for( int iter=0; iter<3; iter++ ){
            
      t.Reset();
      t.SetLastPoint(0);  // first cluster in track position
      t.SetFirstPoint(nRows-1);
      t.ResetCovariance(.1);
      t.SetNumberOfClusters(0);
      for( int i=0; i<nRows; i++ ){
	shapeY2[i]=1.;
	shapeZ2[i]=1.;
	clusterIndex[i]=-1;
	t.SetClusterIndex2(i,-1); 
	t.SetClusterIndex(i,-1); 
      }

     // pick up the clusters
      
      Double_t roady = 20.;
      Double_t roadz = 20.;
      double roadr = 5;

      AliTPCseed t0(t);
      t0.Reset();
      int nClusters = 0;      
      {
	t0.SetRelativeSector(t.GetRelativeSector());
	t0.SetLastPoint(0);  // first cluster in track position
	t0.SetFirstPoint(kMaxRow);
	for (Int_t nr=0; nr<nRows; nr++){ 
	  if( nr<fInnerSec->GetNRows() ) fSectors=fInnerSec;
	  else fSectors=fOuterSec;

	  if( !PropagateToRowHLT(&t0, nr ) ){ break; }
	  if (TMath::Abs(t0.GetSnp())>AliTPCReconstructor::GetMaxSnpTracker()){
	    //cout<<"Snp is too big: "<<t0.GetSnp()<<endl; 
	    continue;
	  }
	  if (!IsActive(t0.GetRelativeSector(),nr)) {
	    continue;
	  }
	  
	  if( iter==0 ){
	    GetShape(&t0,nr); 
	    shapeY2[nr]=t0.GetCurrentSigmaY2();
	    shapeZ2[nr]=t0.GetCurrentSigmaZ2();
	  }

	  AliTPCtrackerRow &krow=GetRow(t0.GetRelativeSector(),nr);
	  if( !krow ) continue; 

	  t.SetClusterIndex2(nr,-3); // foundable
	  t.SetClusterIndex(nr,-3); 

	  AliTPCclusterMI *cl=0;
	  UInt_t uindex = 0;
	  cl = krow.FindNearest2(t0.GetY(),t0.GetZ(),roady+fClExtraRoadY,roadz+fClExtraRoadZ,uindex); 
	  if (!cl ) continue;
	  double dy = cl->GetY()-t0.GetY();
	  double dz = cl->GetZ()-t0.GetZ();
	  double dr = sqrt(dy*dy+dz*dz);
	  if( dr>roadr ){ 
	    //cout<<"row "<<nr<<", best cluster r= "<<dr<<" y,z = "<<dy<<" "<<dz<<endl;
	    continue;
	  }
	  //cout<<"row "<<nr<<", found cluster r= "<<dr<<" y,z = "<<dy<<" "<<dz<<endl;

	  t0.SetClusterPointer(nr, cl);	  
	  clusterIndex[nr] = krow.GetIndex(uindex);	  
	  if( t0.GetFirstPoint()>nr ) t0.SetFirstPoint(nr);
	  t0.SetLastPoint(nr);
	  nClusters++;
	}
      }

      if( nClusters <3 ){
	//cout<<"NOT ENOUGTH CLUSTERS: "<<nClusters<<endl;
	break;
      }
      Int_t basePoints[3] = {t0.GetFirstPoint(),t0.GetFirstPoint(),t0.GetLastPoint()};

      // find midpoint
      {
	Int_t midRow = (t0.GetLastPoint()-t0.GetFirstPoint())/2;
	int dist=200;
	for( int nr=t0.GetFirstPoint()+1; nr< t0.GetLastPoint(); nr++){
	  //if (t0.GetClusterIndex2(nr)<0) continue;
	  if( !t0.GetClusterPointer(nr) ) continue;	  
	  int d = TMath::Abs(nr-midRow);
	  if( d < dist ){
	    dist = d;
	    basePoints[1] = nr;
	  }
	}
      }

      // first fit 3 base points
      if( 1||iter<2 ){
	//cout<<"Fit3: "<<endl;
 	for( int icl=0; icl<3; icl++){
	  int nr = basePoints[icl];
	  int lr=nr;
	  if( nr>=fInnerSec->GetNRows()){ 
	    lr = nr - fInnerSec->GetNRows();
	    fSectors=fOuterSec;
	  } else fSectors=fInnerSec;
	  
	  AliTPCclusterMI *cl=t0.GetClusterPointer(nr);
	  if(!cl){
	    //cout<<"WRONG!!!!"<<endl; 
	    continue;
	  }
	  int iSec = cl->GetDetector() %fkNIS;
	  int rotate = iSec - t.GetRelativeSector();	
	  if( rotate!=0 ){
	    //cout<<"Rotate at row"<<nr<<" to "<<rotate<<" sectors"<<endl;
	    if (!t.Rotate( rotate*fSectors->GetAlpha()) ) {
	      //cout<<"can't rotate "<<endl; 
	      break;
	    }
	    t.SetRelativeSector(iSec);
	  }
	  Double_t x= cl->GetX();
	  if (!t.PropagateTo(x)){
	    //cout<<"can't propagate to x="<<x<<endl; 
	    break;
	  }    
	  t.SetRow(nr); // RS: memorise row
	  if (TMath::Abs(t.GetSnp())>AliTPCReconstructor::GetMaxSnpTracker()){
	    //cout<<"Snp is too big: "<<t.GetSnp()<<endl; 
	    break;
	  }
	  //cout<<"fit3 : row "<<nr<<" ind = "<<clusterIndex[nr]<<endl;
	  
	  t.SetCurrentClusterIndex1(clusterIndex[nr]);
	  t.SetCurrentCluster(cl);
	  t.SetRow(lr);	
	  
	  t.SetErrorY2(shapeY2[nr]);
	  t.SetErrorZ2(shapeZ2[nr]);
	  if( icl==0 ){
	    double cov[15];
	    for( int j=0; j<15; j++ ) cov[j]=0;
	    cov[0]=10;
	    cov[2]=10;
	    cov[5]=.5;
	    cov[9]=.5;
	    cov[14]=1.;
	    t.AliExternalTrackParam::AddCovariance(cov);
	  }
	  if( !UpdateTrack(&t,0) ){
	    //cout<<"Can not update"<<endl;
	    //t.Print();
	    t.SetClusterIndex2(nr,-1); 
	    t.SetClusterIndex(nr,-1); 
	    t.SetClusterPointer(nr,0); 
	    break;
	  }      	
	  //t.SetClusterPointer(nr, cl);
	}
	
	//t.SetLastPoint(t0.GetLastPoint());
	//t.SetFirstPoint(t0.GetFirstPoint());

	//cout<<"Fit: "<<endl;
 	for (Int_t nr=t0.GetLastPoint(); nr>=t0.GetFirstPoint(); nr-- ){
	  int lr=nr;
	  if( nr>=fInnerSec->GetNRows()){ 
	    lr = nr - fInnerSec->GetNRows();
	    fSectors=fOuterSec;
	  } else fSectors=fInnerSec;
	  
	  if(1|| iter<2 ){
	    if( nr == basePoints[0] ) continue;
	    if( nr == basePoints[1] ) continue;
	    if( nr == basePoints[2] ) continue;
	  }
	  AliTPCclusterMI *cl=t0.GetClusterPointer(nr);
	  if(!cl) continue;
	  
	  int iSec = cl->GetDetector() %fkNIS;
	  int rotate = iSec - t.GetRelativeSector();	
	  if( rotate!=0 ){
	    //cout<<"Rotate at row"<<nr<<" to "<<rotate<<" sectors"<<endl;
	    if (!t.Rotate( rotate*fSectors->GetAlpha()) ) {
	      //cout<<"can't rotate "<<endl; 
	      break;
	    }
	    t.SetRelativeSector(iSec);
	  }
	  Double_t x= cl->GetX();
	  if (!t.PropagateTo(x)){
	    //cout<<"can't propagate to x="<<x<<endl; 
	    break;
	  }    
	  if (TMath::Abs(t.GetSnp())>AliTPCReconstructor::GetMaxSnpTracker()){
	    //cout<<"Snp is too big: "<<t.GetSnp()<<endl; 
	    break;
	  } 	
	  t.SetRow(nr); //RS: memorise row
	  //cout<<"fit: row "<<nr<<" ind = "<<clusterIndex[nr]<<endl;

	  t.SetCurrentClusterIndex1(clusterIndex[nr]);
	  t.SetCurrentCluster(cl);
	  t.SetRow(lr);	
	  t.SetErrorY2(shapeY2[nr]);
	  t.SetErrorZ2(shapeZ2[nr]);
	  
	  if( !UpdateTrack(&t,0) ){
	    //cout<<"Can not update"<<endl;
	    //t.Print();
	    t.SetClusterIndex2(nr,-1); 
	    t.SetClusterIndex(nr,-1); 
	    break;
	  }      	
	  //t.SetClusterPointer(nr, cl);
	}
      }
      //cout<<"After iter "<<iter<<": N clusters="<<t.GetNumberOfClusters()<<" : "<<nClusters<<endl;
    }

    //cout<<"fitted track"<<iSeed<<endl;
    //t.Print();
    //cout<<"Statistics: "<<endl;    
    Int_t foundable,found,shared;
    //t.GetClusterStatistic(0,nRows, found, foundable, shared, kTRUE);
    t.GetClusterStatistic(0,nRows, found, foundable); //RS shared info is not used, use faster method
    t.SetNFoundable(foundable);
    //cout<<"found "<<found<<" foundable "<<foundable<<" shared "<<shared<<endl;
    
  }
}


TObjArray * AliTPCtracker::MakeSeedsHLT(const AliESDEvent *hltEvent)
{
  // tracking
  //  
  AliFatal("RS: This method is not yet aware of cluster pointers no present in the in-memory tracks");
  if( !hltEvent ) return 0;  
 

  Int_t nentr=hltEvent->GetNumberOfTracks();
 
  AliInfo(Form("Using %d HLT tracks for seeding",nentr));
 
  TObjArray * seeds = new TObjArray(nentr);

  Int_t nup=fOuterSec->GetNRows()+fInnerSec->GetNRows();
  Int_t index = 0;
  
  Int_t nTr=hltEvent->GetNumberOfTracks();      
    
  for( int itr=0; itr<nTr; itr++ ){
    //if( itr!=97 ) continue;
    const AliExternalTrackParam *param = hltEvent->GetTrack(itr)->GetTPCInnerParam();
    if( !param ) continue;
    //if( TMath::Abs(esdTr->GetSigned1Pt())>1 ) continue;
    //if( TMath::Abs(esdTr->GetTgl())>1. ) continue;
    AliTPCtrack tr;    
    tr.Set(param->GetX(),param->GetAlpha(),param->GetParameter(),param->GetCovariance());
    tr.SetNumberOfClusters(0);
    AliTPCseed * seed = new( NextFreeSeed() ) AliTPCseed(tr);

    Double_t alpha=seed->GetAlpha();// - fSectors->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    //
    seed->SetRelativeSector(Int_t(alpha/fSectors->GetAlpha()+0.0001)%fN);
    Double_t alphaSec = fSectors->GetAlpha() * seed->GetRelativeSector() + fSectors->GetAlphaShift();
 
    if (alphaSec >= TMath::Pi()) alphaSec -= 2.*TMath::Pi();
    if (alphaSec < -TMath::Pi()) alphaSec += 2.*TMath::Pi();

    seed->Rotate(alphaSec - alpha);
    
    seed->SetPoolID(fLastSeedID);
    seed->SetIsSeeding(kTRUE);
    seed->SetSeed1(nup-1);
    seed->SetSeed2(nup-2);
    seed->SetSeedType(0);
    seed->SetFirstPoint(-1);
    seed->SetLastPoint(-1);
    seeds->AddLast(seed); // note, track is seed, don't free the seed
    index++;
    //if( index>3 ) break;
  }
 

  fSectors = fOuterSec;
  
  TrackFollowingHLT(seeds );

  nTr = seeds->GetEntriesFast();
  for( int itr=0; itr<nTr; itr++ ){
    AliTPCseed * seed = (AliTPCseed*) seeds->UncheckedAt(itr);
    if( !seed ) continue;
    //FollowBackProlongation(*seed,0);
    // cout<<seed->GetNumberOfClusters()<<endl;
    Int_t foundable,found;//,shared;
    //    seed->GetClusterStatistic(0,nup, found, foundable, shared, kTRUE); //RS shared info is not used, use faster method
    seed->GetClusterStatistic(0,nup, found, foundable);
    seed->SetNFoundable(foundable);
    //cout<<"found "<<found<<" foundable "<<foundable<<" shared "<<shared<<endl;
    //if ((found<0.55*foundable)  || shared>0.5*found ){// || (seed->GetSigmaY2()+seed->GetSigmaZ2())>0.5){
    //MarkSeedFree(seeds->RemoveAt(itr)); 
    //continue;
    //}
    if (seed->GetNumberOfClusters()<30 || 
	seed->GetNumberOfClusters() < seed->GetNFoundable()*0.6 || 
	seed->GetNShared()>0.4*seed->GetNumberOfClusters() ) {
      MarkSeedFree(seeds->RemoveAt(itr)); 
      continue;
    }      

    for( int ir=0; ir<nup; ir++){
      AliTPCclusterMI *c = seed->GetClusterPointer(ir);      
      if( c ) c->Use(10);
    }
  }
  std::cout<<"\n\nHLT tracks left: "<<seeds->GetEntries()<<" out of "<<hltEvent->GetNumberOfTracks()<<endl<<endl;
  return seeds;    
}

void AliTPCtracker::FillClusterOccupancyInfo()
{
  //fill the cluster occupancy info into the ESD friend
  AliESDfriend* esdFriend = static_cast<AliESDfriend*>(fEvent->FindListObject("AliESDfriend"));
  if (!esdFriend) return;

  for (Int_t isector=0; isector<18; isector++){
    AliTPCtrackerSector &iroc = fInnerSec[isector];
    AliTPCtrackerSector &oroc = fOuterSec[isector];
    //all clusters
    esdFriend->SetNclustersTPC(isector,   iroc.GetNClInSector(0));
    esdFriend->SetNclustersTPC(isector+18,iroc.GetNClInSector(1));
    esdFriend->SetNclustersTPC(isector+36,oroc.GetNClInSector(0));
    esdFriend->SetNclustersTPC(isector+54,oroc.GetNClInSector(1));
    //clusters used in tracking
    esdFriend->SetNclustersTPCused(isector,    iroc.GetNClUsedInSector(0));
    esdFriend->SetNclustersTPCused(isector+18, iroc.GetNClUsedInSector(1));
    esdFriend->SetNclustersTPCused(isector+36, oroc.GetNClUsedInSector(0));
    esdFriend->SetNclustersTPCused(isector+54, oroc.GetNClUsedInSector(1));
  }
}

void AliTPCtracker::FillSeedClusterStatCache(const AliTPCseed* seed)
{
  //RS fill cache for seed's cluster statistics evaluation buffer
  for (int i=kMaxRow;i--;) {
    Int_t index = seed->GetClusterIndex2(i);
    fClStatFoundable[i] = (index!=-1);
    fClStatFound[i] = fClStatShared[i] = kFALSE;
    if (index<0 || (index&0x8000)) continue;
    const AliTPCclusterMI* cl = GetClusterMI(index);
    if (cl->IsUsed(10)) fClStatShared[i] = kTRUE;
  }
}

void AliTPCtracker::GetCachedSeedClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2) const
{
  //RS calculation of cached seed statistics
  found       = 0;
  foundable   = 0;
  shared      = 0;
  //
  for (Int_t i=first;i<last; i++){
    if (fClStatFoundable[i]) foundable++;
    else continue;
    //
    if (fClStatFound[i]) found++;
    else continue;
    //
    if (fClStatShared[i]) {
      shared++;
      continue;
    }
    if (!plus2) continue; //take also neighborhoud
    //
    if ( i>0 && fClStatShared[i-1]) {
      shared++;
      continue;
    }
    if ( i<(kMaxRow-1) && fClStatShared[i+1]) {
      shared++;
      continue;
    }    
  }
}

void AliTPCtracker::GetSeedClusterStatistic(const AliTPCseed* seed, Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2) const
{
  // RS: get cluster stat. on given region, faster for small regions than cachin the whole seed stat.
  //
  found       = 0;
  foundable   = 0;
  shared      =0;
  for (Int_t i=first;i<last; i++){
    Int_t index = seed->GetClusterIndex2(i);
    if (index!=-1) foundable++;
    if (index&0x8000) continue;
    if (index>=0) found++;
    else continue;
    const AliTPCclusterMI* cl = GetClusterMI(index);
    if (cl->IsUsed(10)) {
      shared++;
      continue;
    }
    if (!plus2) continue; //take also neighborhoud
    //
    if (i>0) {
      index = seed->GetClusterIndex2(i-1);
      cl = index<0 ? 0:GetClusterMI(index);
      if (cl && cl->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    if (i<(kMaxRow-1)) {
      index = seed->GetClusterIndex2(i+1);
      cl = index<0 ? 0:GetClusterMI(index);      
      if (cl && cl->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    
  }
}

Bool_t AliTPCtracker::DistortX(const AliTPCseed* seed, double& x, int row)
{
  // distort X by the distorion at track location on given row
  //
  //RS:? think to unset fAccountDistortions if no maps are used
  if (!AliTPCReconstructor::GetRecoParam()->GetUseCorrectionMap()) return kTRUE;
  double xyz[3]; 
  int rowInp = row;//RS
  if (!seed->GetYZAt(x,AliTracker::GetBz(),&xyz[1])) return kFALSE; //RS:? Think to cache fBz?
  xyz[0] = x;
  int roc = seed->GetRelativeSector();
  if (seed->GetZ()<0) roc += 18;
  if (row>62) {
    roc += 36;
    row -= 63;
  }
  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform();
  if (!transform) AliFatal("Tranformations not in calibDB");
  x += transform->GetCorrMapComponent(roc,row,xyz,0);
  return kTRUE;
}

Double_t AliTPCtracker::GetDistortionX(double x, double y, double z, int sec, int row)
{
  // get X distortion at location on given row
  //
  if (!AliTPCReconstructor::GetRecoParam()->GetUseCorrectionMap()) return 0;
  double xyz[3] = {x,y,z}; 
  int rowInp = row;
  int secInp = sec;
  if (z<0) sec += 18;
  if (row>62) {
    sec += 36;
    row -= 63;
  }
  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform() ;
  if (!transform) AliFatal("Tranformations not in calibDB");
  return transform->GetCorrMapComponent(sec,row,xyz,0);
}

Double_t AliTPCtracker::GetYSectEdgeDist(int sec, int row, double y, double z) 
{
  // get the signed  shift for maxY of the sector/row accounting for distortion
  // Slow way, to speed up
  if (!AliTPCReconstructor::GetRecoParam()->GetUseCorrectionMap()) return 0;
  double ymax = 0.9*GetMaxY(row); // evaluate distortions at 5% of pad lenght from the edge
  if (y<0) ymax = -ymax;
  double xyz[3] = {GetXrow(row),ymax,z}; 
  if (z<0) sec += 18;
  if (row>62) {
    sec += 36;
    row -= 63;
  }
  AliTPCcalibDB * calibDB = AliTPCcalibDB::Instance();
  AliTPCTransform *transform = calibDB->GetTransform();
  if (!transform) AliFatal("Tranformations not in calibDB");
  // change of distance from the edge due to the X shift 
  double dxtg = transform->GetCorrMapComponent(sec,row,xyz,0)*AliTPCTransform::GetMaxY2X();
  double dy = transform->GetCorrMapComponent(sec,row,xyz,1);
  return dy + (y>0?dxtg:-dxtg);
  //
}

Int_t AliTPCtracker::GetTrackSector(double alpha)
{
  //convert alpha to sector
  if (alpha<0) alpha += TMath::Pi()*2;
  int sec = alpha/(TMath::Pi()/9);
  return sec;
}

void AliTPCtracker::CleanESDFriendsObjects(AliESDEvent* esd)
{
  // RS: remove seeds stored in friend's calib object contained w/o changing its ownership
  //
  AliInfo("Removing own seeds from friend tracks");
  AliESDfriend* esdF = esd->FindFriend();
  int ntr = esd->GetNumberOfTracks();
  for (int itr=ntr;itr--;) {
    AliESDtrack* trc = esd->GetTrack(itr);
    AliESDfriendTrack* trcF = (AliESDfriendTrack*)trc->GetFriendTrack();
    if (!trcF) continue;
    AliTPCseed* seed = (AliTPCseed*)trcF->GetTPCseed();
    if (seed) {
      trcF->RemoveCalibObject((TObject*)seed);
      seed->SetClusterOwner(kFALSE);
      seed->SetClustersArrayTMP(0);
    }
  }
  //
}

//__________________________________________________________________
void AliTPCtracker::CleanESDTracksObjects(TObjArray* trcList)
{
  // RS: remove seeds stored in friend's calib object contained w/o changing its ownership
  //
  int ntr = trcList->GetEntriesFast();
  for (int itr=0;itr<ntr;itr++) {
    TObject *obj = trcList->At(itr);
    AliESDfriendTrack* trcF =  (obj->IsA()==AliESDtrack::Class()) ? 
      (AliESDfriendTrack*)((AliESDtrack*)obj)->GetFriendTrack() : (AliESDfriendTrack*)obj;
    if (!trcF) continue;
    AliTPCseed* seed = (AliTPCseed*)trcF->GetTPCseed();
    if (seed) {
      trcF->RemoveCalibObject((TObject*)seed);
      seed->SetClusterOwner(kFALSE);
      seed->SetClustersArrayTMP(0);
    }
  }
  //
}

//__________________________________________________________________
void AliTPCtracker::AddSystCovariance(AliTPCseed* t)
{
  // calculate correction to covariance matrix at given point
  //
  float amplCorr  = AliTPCReconstructor::GetRecoParam()->GetSystCovAmplitude(); // amplitude of correction
  if (amplCorr<1e-6) return;                                                    // no correction
  float fluctCorr = AliTPCReconstructor::GetRecoParam()->GetDistFluctCorrelation(); // distortion fluctuations correlation 

  double corY[6]={0},corZ[3]={0},jacobC = kB2C*GetBz()/2;
  float xArr[kMaxRow][3], xRef = t->GetX();
  float gammaY[kMaxRow],gammaZ[kMaxRow]; // sysErr^2 / totErr^4
  int nclFit = 0;
  for (int ip=kMaxRow;ip--;) {
    int index = t->GetClusterIndex2(ip);
    if ( index<0 || (index&0x8000)) continue; // missing or discared cluster
    const AliTPCTrackerPoints::Point *trpoint =t->GetTrackPoint(ip);
    xArr[nclFit][0] = 1.;
    double dx = trpoint->GetX()-xRef;
    xArr[nclFit][1] = dx;
    xArr[nclFit][2] = dx*dx;
    gammaY[nclFit] = trpoint->GetErrYSys2TotSq(); // sysE/totE^2
    gammaZ[nclFit] = trpoint->GetErrZSys2TotSq(); // sysE/totE^2
    nclFit++;
  }
  // X^{T} W M W^{T} X matrix
  int cnt=0; 
  for (int row=0;row<3;row++) {
    for (int col=0;col<=row;col++) { // matrix is symmetric
      for (int ip=0;ip<nclFit;ip++) {
	for (int jp=0;jp<nclFit;jp++) {
	  double corr = (ip==jp ? 1.0:fluctCorr) * amplCorr;
	  double xprod = xArr[ip][row]*xArr[jp][col]*corr;
	  corY[cnt] += xprod*gammaY[ip]*gammaY[jp];
	  if (row<2) corZ[cnt] += xprod*gammaZ[ip]*gammaZ[jp];
	}
      }
      cnt++;
    }
  }
  //
  // account for jacobian C->P[4]
  corY[3] *= jacobC;
  corY[4] *= jacobC;
  corY[5] *= jacobC*jacobC;
  //
  //
  double *covP = (double*)t->GetCovariance(),
    &C00=covP[0],
    &C10=covP[1] ,&C11=covP[2],
    &C20=covP[3] ,&C21=covP[4] ,&C22=covP[5],
    &C30=covP[6] ,&C31=covP[7] ,&C32=covP[8] ,&C33=covP[9],
    &C40=covP[10],&C41=covP[11],&C42=covP[12],&C43=covP[13],&C44=covP[14];
  double &g00=corY[0],&g11=corZ[0],&g20=corY[1],&g22=corY[2],&g31=corZ[1],&g33=corZ[2],&g40=corY[3],&g42=corY[4],&g44=corY[5];

  double
    cg00=C00*g00 + C20*g20 + C40*g40, cg01=C10*g11 + C30*g31, cg02=C00*g20 + C20*g22 + C40*g42, cg03=C10*g31 + C30*g33, cg04=C00*g40 + C20*g42 + C40*g44,
    cg10=C10*g00 + C21*g20 + C41*g40, cg11=C11*g11 + C31*g31, cg12=C10*g20 + C21*g22 + C41*g42, cg13=C11*g31 + C31*g33, cg14=C10*g40 + C21*g42 + C41*g44,
    cg20=C20*g00 + C22*g20 + C42*g40, cg21=C21*g11 + C32*g31, cg22=C20*g20 + C22*g22 + C42*g42, cg23=C21*g31 + C32*g33, cg24=C20*g40 + C22*g42 + C42*g44,
    cg30=C30*g00 + C32*g20 + C43*g40, cg31=C31*g11 + C33*g31, cg32=C30*g20 + C32*g22 + C43*g42, cg33=C31*g31 + C33*g33, cg34=C30*g40 + C32*g42 + C43*g44,
    cg40=C40*g00 + C42*g20 + C44*g40, cg41=C41*g11 + C43*g31, cg42=C40*g20 + C42*g22 + C44*g42, cg43=C41*g31 + C43*g33, cg44=C40*g40 + C42*g42 + C44*g44;
  //
  double 
    a00 = C00*cg00 + C10*cg01 + C20*cg02 + C30*cg03 + C40*cg04, //  row0
    a10 = C00*cg10 + C10*cg11 + C20*cg12 + C30*cg13 + C40*cg14, //  row1 
    a11 = C10*cg10 + C11*cg11 + C21*cg12 + C31*cg13 + C41*cg14,
    a20 = C00*cg20 + C10*cg21 + C20*cg22 + C30*cg23 + C40*cg24, //  row2 
    a21 = C10*cg20 + C11*cg21 + C21*cg22 + C31*cg23 + C41*cg24, 
    a22 = C20*cg20 + C21*cg21 + C22*cg22 + C32*cg23 + C42*cg24,
    a30 = C00*cg30 + C10*cg31 + C20*cg32 + C30*cg33 + C40*cg34, //  row3
    a31 = C10*cg30 + C11*cg31 + C21*cg32 + C31*cg33 + C41*cg34, 
    a32 = C20*cg30 + C21*cg31 + C22*cg32 + C32*cg33 + C42*cg34, 
    a33 = C30*cg30 + C31*cg31 + C32*cg32 + C33*cg33 + C43*cg34,
    a40 = C00*cg40 + C10*cg41 + C20*cg42 + C30*cg43 + C40*cg44, //  row4 
    a41 = C10*cg40 + C11*cg41 + C21*cg42 + C31*cg43 + C41*cg44, 
    a42 = C20*cg40 + C21*cg41 + C22*cg42 + C32*cg43 + C42*cg44, 
    a43 = C30*cg40 + C31*cg41 + C32*cg42 + C33*cg43 + C43*cg44, 
    a44 = C40*cg40 + C41*cg41 + C42*cg42 + C43*cg43 + C44*cg44;
  //
  // make sure diagonal elements are positive
  if (a00<0) a00 = 0;
  if (a11<0) a11 = 0;
  if (a22<0) a22 = 0;
  if (a33<0) a33 = 0;
  if (a44<0) a44 = 0;
  //    
  C00 += a00;
  C10 += a10;  C11 += a11;
  C20 += a20;  C21 += a21;  C22 += a22;
  C30 += a30;  C31 += a31;  C32 += a32;  C33 += a33;
  C40 += a40;  C41 += a41;  C42 += a42;  C43 += a43;  C44 += a44;
  
}
