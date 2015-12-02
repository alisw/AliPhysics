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


/// \class AliTPCcalibAlignInterpolation
/// Class to produce TPC time dependent space point distortion maps using the ITS, TRD and TOF 
/// as a reference detector
///  
/// Related to task https://alice.its.cern.ch/jira/browse/ATO-108
///  - code created addopting compiled macro for open gating grid analysis form TPC git:
///    $NOTES/SpaceChargeDistortion/code/spaceChargeDistortions.C
/// 
/// \author Marian Ivanov,  marian.ivanov@cern.ch

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "TTreeStream.h"
#include "AliESDfriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTrackPointArray.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliTrackerBase.h"
#include "AliGeomManager.h"
#include "TVectorF.h"
#include "TStopwatch.h"
#include "TProfile.h"
#include "TGraphErrors.h"
//#include "THnBase.h"
#include "THn.h"
#include "AliSysInfo.h"
#include "TMatrixD.h"
 #include "TF1.h"
#include "TDatabasePDG.h"
#include "TTreeStream.h"
#include "TStatToolkit.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"
#include "AliTPCRecoParam.h"
#include "AliTPCreco.h"
#include "AliTPCcalibAlignInterpolation.h"
#include "AliPID.h"
#include "TSystem.h"

const Int_t AliTPCcalibAlignInterpolation_kMaxPoints=500;

ClassImp(AliTPCcalibAlignInterpolation)


AliTPCcalibAlignInterpolation::AliTPCcalibAlignInterpolation() : 
  AliTPCcalibBase(),
  fOnTheFlyFill(0),  // flag - on the fly filling of histograms
  fHisITSDRPhi(0),      // TPC-ITS residual histograms
  fHisITSTRDDRPhi(0),   // TPC-ITS+TRD residual histograms
  fHisITSTOFDRPhi(0),   // TPC-ITS_TOF residual histograms
  fHisITSDZ(0),      // TPC-ITS residual histograms
  fHisITSTRDDZ(0),   // TPC-ITS+TRD residual histograms
  fHisITSTOFDZ(0),   // TPC-ITS_TOF residual histograms
  fRhoTPC(0.9e-3),
  fX0TPC(28.94),
  fStreamer(0),         // calibration streamer 
  fStreamLevelTrack(0),      // stream level - In mode 0 only basic information needed for calibration  stored (see EStream
  fSyswatchStep(100),      // dump system resource information after  fSyswatchStep tracks
  fTrackCounter(0)           // processed track counter
{
  
}   
AliTPCcalibAlignInterpolation::AliTPCcalibAlignInterpolation(const Text_t *name, const Text_t *title, Bool_t onTheFlyFill):
  AliTPCcalibBase(),
  fOnTheFlyFill(onTheFlyFill),  // flag - on the fly filling of histograms
  fHisITSDRPhi(0),      // TPC-ITS residual histograms
  fHisITSTRDDRPhi(0),   // TPC-ITS+TRD residual histograms
  fHisITSTOFDRPhi(0),   // TPC-ITS_TOF residual histograms  
  fHisITSDZ(0),      // TPC-ITS residual histograms
  fHisITSTRDDZ(0),   // TPC-ITS+TRD residual histograms
  fHisITSTOFDZ(0),   // TPC-ITS_TOF residual histograms
  fRhoTPC(0.9e-3),
  fX0TPC(28.94),
  fStreamer(0),         // calibration streamer 
  fStreamLevelTrack(0),      // stream level - In mode 0 only basic information needed for calibration  stored (see EStream
  fSyswatchStep(100),      // dump system resource information after  fSyswatchStep tracks  
  fTrackCounter(0)           // processed track counter
{
  // create output histograms
  SetName(name);
  SetTitle(title);
  if (onTheFlyFill) CreateResidualHistosInterpolation();
}   

AliTPCcalibAlignInterpolation::~AliTPCcalibAlignInterpolation(){
  //
  //
  //
  if (fStreamer){
    // fStreamer->GetFile()->Close();
    fStreamer->GetFile()->cd();
    if (fHisITSDRPhi) fHisITSDRPhi->Write();
    if (fHisITSTRDDRPhi) fHisITSTRDDRPhi->Write();
    if (fHisITSTOFDRPhi) fHisITSTOFDRPhi->Write();
  }
  delete fStreamer;
  fStreamer=0;
  delete fHisITSDRPhi;
  delete fHisITSTRDDRPhi;
  delete fHisITSTOFDRPhi;
}


void AliTPCcalibAlignInterpolation::Terminate(){
  //
  // Terminate function
  // call base terminate + Eval of fitters
  //
  Info("AliTPCcalibAlignInterpolation","Terminate");
  if (fStreamer){
    // fStreamer->GetFile()->Close();
    fStreamer->GetFile()->cd();
    if (fHisITSDRPhi) fHisITSDRPhi->Write();
    if (fHisITSTRDDRPhi) fHisITSTRDDRPhi->Write();
    if (fHisITSTOFDRPhi) fHisITSTOFDRPhi->Write();
  }
  delete fStreamer;
  fStreamer=0;

  AliTPCcalibBase::Terminate();
}


Bool_t  AliTPCcalibAlignInterpolation::RefitITStrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackITS, 
						     Double_t &chi2, Double_t &npoints, int* sortInd){
  //
  // Makes a refit of the ITS track
  // Input: AliESDfriendTrack, particle mass, outer ITS TrackParam 
  // Output: AliExternalTrackParam of the ITS track refit at the last layer of ITS
  //
  const Double_t sigma2[6][3] = {
    {0.002*0.002, 0,  0.013*0.013},
    {0.002*0.002, 0,  0.013*0.013},
    {0.050*0.050, 0,  0.050*0.050},
    {0.050*0.050, 0,  0.050*0.050},
    {0.003*0.003, 0,  0.100*0.100},
    {0.003*0.003, 0,  0.100*0.100},
  };    // ITS intrincsic resolution in (y,z)  - error from the points can be used SD layer 2-3 sighnificantly bigger error
  // !!!! We should set ITS error parameterization form outside !!!!
  const Double_t kMaxRadius=50;
  static Int_t sortedIndexLoc[AliTPCcalibAlignInterpolation_kMaxPoints]={0};
  chi2=0;
  npoints=0; 
  //
  if (friendTrack->GetITSOut()==NULL) return kFALSE;  
  //
  trackITS = *((AliExternalTrackParam*)friendTrack->GetITSOut());
  // Reset track to the vertex
  //if (!AliTrackerBase::PropagateTrackToBxByBz(&trackITS,0,mass,1,kFALSE)) return kFALSE;
  if (!AliTrackerBase::PropagateTrackParamOnlyTo(&trackITS,0.,5,kFALSE)) return kFALSE;
  trackITS.ResetCovariance(1000.);
  
  // Get space points
  AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
  if (!pointarray){
    printf("Space points are not stored in the friendTrack!\n");
    return kFALSE;
  };
  Int_t nPoints = pointarray->GetNPoints();  // # space points of all available detectors                                            
                                             // Sort space points first
  int *sortedIndex = sortInd;
  if (!sortedIndex) {
    SortPointArray(pointarray, sortedIndexLoc);  // space point indices sorted by radius in increasing order
    sortedIndex = sortedIndexLoc;
  }
  //
  // Propagate track through ITS space points
  AliTrackPoint spacepoint;
  Int_t volId=0,modId=0,layerId=0;
  
  for (Int_t iPoint=0;iPoint<nPoints;iPoint++){
    pointarray->GetPoint(spacepoint,sortedIndex[iPoint]);
    int lr = AliGeomManager::VolUIDToLayer(spacepoint.GetVolumeID())-1;
    if (lr<0||lr>=6) continue;
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    trackITS.Global2LocalPosition(xyz,alpha);
    if (xyz[0]>kMaxRadius) break;  // use only ITS points - maybe we should indexes of elements
    if (!trackITS.Rotate(alpha)) return kFALSE;
    if (!AliTrackerBase::PropagateTrackToBxByBz(&trackITS,xyz[0],mass,1,kFALSE)) return kFALSE;
    Double_t pos[2] = {xyz[1], xyz[2]};
    const Double_t* cov = sigma2[lr];
    volId = spacepoint.GetVolumeID();
    //    layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    //     if (layerId ==AliGeomManager::kSDD1 || layerId ==AliGeomManager::kSDD2) cov[0]*=16.; cov[2]*=16.;}      
    double chi2cl = trackITS.GetPredictedChi2(pos,cov);
    chi2 += chi2cl;
    npoints++;
    if (!trackITS.Update(pos,cov)) return kFALSE;
  }
  return npoints>0;
}


Bool_t AliTPCcalibAlignInterpolation::RefitTRDtrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackTRD, 
						    Double_t &chi2, Double_t &npoints, Int_t* sortInd){
  //
  // Makes a refit of the TRD track using TOF and TRD points
  // Input: AliESDfriendTrack, particle mass, inner TRD TrackParam 
  // Output: AliExternalTrackParam of the TRD track refit - in the first layer of TRD
  // Here we forgot about the tiliting pads of TRD - we assume delta Z and delta y are uncorelated
  //      given approximation is in average tru - in case avearaging of significantly bigger than pad length
  //  
  const Double_t sigmaTRD2[2] = {0.04*0.04, 5*5};
  const Double_t sigmaTOF2[2] = {1, 1};
  static Int_t sortedIndexLoc[AliTPCcalibAlignInterpolation_kMaxPoints]={0};
  const Double_t kMaxRadius=390;
  const Double_t kMinRadius=280;
  //
  chi2=0; 
  npoints=0;  
  //
  if (friendTrack->GetTRDIn() == NULL) return kFALSE;
  trackTRD = *((AliExternalTrackParam*)friendTrack->GetTRDIn());
  
  
  // Reset track outside TRD
  if (!AliTrackerBase::PropagateTrackParamOnlyTo(&trackTRD,kMaxRadius,5,kFALSE)) return kFALSE;
  //if (!AliTrackerBase::PropagateTrackToBxByBz(&trackTRD,kMaxRadius,mass,1,kFALSE)) return kFALSE;
  trackTRD.ResetCovariance(1000.);
      
  // Get space points
  AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
  if (!pointarray){
    printf("Space points are not stored in the friendTrack!\n");
    return kFALSE;
  };
  Int_t nPoints = pointarray->GetNPoints();  // # space points of all available detectors
                                             // Sort space points first
  int *sortedIndex = sortInd;
  if (!sortedIndex) {
    SortPointArray(pointarray, sortedIndexLoc);  // space point indices sorted by radius in increasing order
    sortedIndex = sortedIndexLoc;
  }
  
  // Propagate track through TRD space points
  AliTrackPoint spacepoint;
  Int_t volId,modId,layerId, npfit=0;
  for (Int_t iPoint=nPoints;iPoint--;) {
    pointarray->GetPoint(spacepoint,sortedIndex[iPoint]);
    volId = spacepoint.GetVolumeID();
    layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if (layerId <AliGeomManager::kTRD1) break;
    if (layerId>AliGeomManager::kTOF) continue;
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    trackTRD.Global2LocalPosition(xyz,alpha);
    if (xyz[0]<kMinRadius) break;  // use only TRD points
    if (!trackTRD.Rotate(alpha)) return kFALSE;
    if (!AliTrackerBase::PropagateTrackToBxByBz(&trackTRD,xyz[0],mass,1,kFALSE)) return kFALSE;
    Double_t pos[2] = {xyz[1], xyz[2]};
    Double_t cov[3] = {sigmaTRD2[0],0,sigmaTRD2[1]};
    if (layerId==AliGeomManager::kTOF) {cov[0]=sigmaTOF2[0]; cov[2]=sigmaTOF2[1];};
    double chi2cl = trackTRD.GetPredictedChi2(pos,cov);
    chi2 += chi2cl;
    if (!trackTRD.Update(pos,cov)) return kFALSE;
    npfit++;
  }
  npoints = npfit;
  return npoints>0;
}


Bool_t  AliTPCcalibAlignInterpolation::RefitTOFtrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackTOF, 
						     Double_t &chi2, Double_t &npoints, Int_t* sortInd){
  //
  // Makes a refit of the TRD track
  // Input: AliESDfriendTrack, particle mass, OUTER ITS track - propagated to the TOF point and updated by TOF point 
  // Output: AliExternalTrackParam of the TOF track refit - at the TOF point
  Double_t sigma2[2] = {1., 1.};      // should be parameterized
  const Double_t kTOFRadius = 370;
  //
  chi2=0; 
  npoints=0;
  //
  static Int_t sortedIndexLoc[AliTPCcalibAlignInterpolation_kMaxPoints]={0};
  //  if (!AliTrackerBase::PropagateTrackParamOnlyTo(&trackTOF,kTOFRadius,15,kTRUE)) return kFALSE;
  if (!AliTrackerBase::PropagateTrackToBxByBz(&trackTOF,kTOFRadius,mass,10,kTRUE)) return kFALSE;
  // RS why don't we reset the cov. matrix here?
  Int_t volId,modId,layerId;
      
  // Get space points
  AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
  if (!pointarray){
    printf("Space points are not stored in the friendTrack!\n");
    return kFALSE;
  };
  Int_t nPoints = pointarray->GetNPoints();  // # space points of all available detectors
                                             // Sort space points first
  int *sortedIndex = sortInd;
  if (!sortedIndex) {
    SortPointArray(pointarray, sortedIndexLoc);  // space point indices sorted by radius in increasing order
    sortedIndex = sortedIndexLoc;
  }

  // Propagate track through TRD space points
  AliTrackPoint spacepoint;
  int npfit = 0;
  for (Int_t iPoint=nPoints;iPoint--;){  
    pointarray->GetPoint(spacepoint,sortedIndex[iPoint]);
    volId = spacepoint.GetVolumeID();
    layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if (layerId !=AliGeomManager::kTOF) continue;
    
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    trackTOF.Global2LocalPosition(xyz,alpha);
    if (!trackTOF.Rotate(alpha)) return kFALSE;
    if (!AliTrackerBase::PropagateTrackToBxByBz(&trackTOF,xyz[0],mass,1,kFALSE)) return kFALSE;
    Double_t pos[2] = {xyz[1], xyz[2]};
    Double_t cov[3] = {sigma2[0],0,sigma2[1]};
    double chi2cl = trackTOF.GetPredictedChi2(pos,cov);
    chi2 += chi2cl;
    if (!trackTOF.Update(pos,cov)) return kFALSE;
    npfit++;
    break; // there is just 1 TOF poitn
  }
  npoints = npfit;
  return npoints>0;
}

Bool_t  AliTPCcalibAlignInterpolation::SortPointArray(AliTrackPointArray *pointarray, Int_t * sortedIndex){
  //
  // Fill array of indexes to the pointArray (array sorted in increasing order)
  //
  if (sortedIndex==NULL) return kFALSE;
  Int_t nPoints = pointarray->GetNPoints();
  if (!nPoints) return kFALSE;
  Double_t rp[nPoints];
  const float* x = pointarray->GetX();
  const float* y = pointarray->GetY();
  for (Int_t iPoint=nPoints;iPoint--;) rp[iPoint] = x[iPoint]*x[iPoint]+y[iPoint]*y[iPoint];
  TMath::Sort(nPoints,rp,sortedIndex,kFALSE);
  return kTRUE;
}



void AliTPCcalibAlignInterpolation::ProcessStandalone(const char * inputList){
  //
  // Process ESD information standalone without full colibration train
  // Input:
  //   inputList - list of the input ESD files
  //
  // code from test macro to be set here

}



void  AliTPCcalibAlignInterpolation::Process(AliESDEvent *esdEvent){
  //
  // Create distortion maps out of residual histograms of ITS-TRD interpolation and TPC space points
  // JIRA ticket: https://alice.its.cern.ch/jira/browse/ATO-108
  //
  const Double_t kMaxChi2=10;         // max track/track chi2 
  const Double_t kMaxAlignTolerance=0.1;   // max track/track chi2 
  const Double_t kMaxSNP = 0.8; // max snp tolerated
  //
  static Bool_t firstCall = kTRUE;
  if (firstCall) {
    firstCall = kFALSE;
    ExtractTPCGasData();
  }
  //
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(esdEvent->FindListObject("AliESDfriend"));
  if (!esdFriend) return;
  if (esdFriend->TestSkipBit()) return;
  Int_t nPrimTracks= (esdEvent->GetPrimaryVertex()!=NULL)? esdEvent->GetPrimaryVertex()->GetNContributors():0;
  Int_t nPrimTracksSPD= (esdEvent->GetPrimaryVertexSPD()!=NULL)? esdEvent->GetPrimaryVertexSPD()->GetNContributors():0;
  Int_t nTracks = esdEvent->GetNumberOfTracks();  // Get number of tracks in ESD
  Int_t nSPD=  esdEvent->GetMultiplicity()->GetNumberOfITSClusters(0,1);
  Int_t nSDD=  esdEvent->GetMultiplicity()->GetNumberOfITSClusters(2,3);
  Int_t nSSD=  esdEvent->GetMultiplicity()->GetNumberOfITSClusters(4,5);
  if (nTracks==0) return;
  if (!fStreamer) fStreamer = new TTreeSRedirector("ResidualHistos.root","recreate");
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  TVectorF vecNClTPC(72);
  TVectorF vecNClTPCused(72);
  for (Int_t isec=0; isec<72;isec++){
    vecNClTPC[isec]=esdFriend->GetNclustersTPC(isec);
    vecNClTPCused[isec]=esdFriend->GetNclustersTPCused(isec);
  }
  Long64_t gid = esdEvent->GetHeader()->GetEventIdAsLong(); 
  Int_t timeStamp= esdEvent->GetTimeStamp();
  (*fStreamer)<<"eventInfo"<< // store event info - used to calculate per sector currents
    "gid="<<gid<<
    "timeStamp="<<timeStamp<<
    "nSPD="<<nSPD<<
    "nSDD="<<nSDD<<
    "nSSD="<<nSSD<<
    "nPrimTrack="<<nPrimTracks<<
    "nPrimTrackSPD="<<nPrimTracksSPD<<
    "nTracks="<<nTracks<<
    "vecNClTPC.="<<&vecNClTPC<<
    "vecNClTPCused.="<<&vecNClTPCused<<
    "\n";
  

  //
  const Int_t nPointsAlloc=AliTPCcalibAlignInterpolation_kMaxPoints; 
  const Int_t kMaxLayer=kMaxRow;
  AliExternalTrackParam trackArrayITS[kMaxLayer];
  AliExternalTrackParam trackArrayTRD[kMaxLayer];
  AliExternalTrackParam trackArrayTOF[kMaxLayer];
  AliExternalTrackParam trackArrayITSTRD[kMaxLayer];
  AliExternalTrackParam trackArrayITSTOF[kMaxLayer];
  AliTPCclusterMI clusterArray[kMaxLayer];
  //
  //MakeResidualHistosInterpolation();
  //
  Int_t sortedIndex[AliTPCcalibAlignInterpolation_kMaxPoints];
  TVectorF deltaITS0(kMaxLayer), deltaITS1(kMaxLayer); 
  TVectorF deltaTRD0(kMaxLayer), deltaTRD1(kMaxLayer); 
  TVectorF deltaTOF0(kMaxLayer), deltaTOF1(kMaxLayer); 
  TVectorF vecR(kMaxLayer), vecPhi(kMaxLayer), vecZ(kMaxLayer), vecSec(kMaxLayer);
  static int evCnt=0;
  Bool_t backupUseComposedCorrection = transform->GetCurrentRecoParamNonConst()->GetUseComposedCorrection();
  transform->GetCurrentRecoParamNonConst()->SetUseComposedCorrection(kFALSE);
  
  for (Int_t iTrack=0;iTrack<nTracks;iTrack++){ // Track loop
    // 0.) For each track in each event, get the AliESDfriendTrack
    AliESDtrack *esdTrack = esdEvent->GetTrack(iTrack);
    AliESDfriendTrack *friendTrack = (AliESDfriendTrack*)esdTrack->GetFriendTrack();
    if (!friendTrack) continue;      
    if (esdTrack->GetITSNcls()<4 || esdTrack->GetTPCNcls()<15) continue;
    Double_t mass = esdTrack->GetMass();  // particle mass    
    Double_t tofDiff=esdTrack->GetTOFExpTDiffSpec(AliPID::kPion);
    // Get TPC seed
    TObject *calibObject=0;
    AliTPCseed *seed = (AliTPCseed*)friendTrack->GetTPCseed();
    if (!seed) continue;
    //
    // 1.) Start with AliExternalTrackParam *ITSOut and *TRDIn 
    //
    AliExternalTrackParam paramITS;
    Double_t itsChi2=0, itsNCl=0;
    AliExternalTrackParam paramTRD;
    Double_t trdChi2=0, trdNCl=0;
    AliExternalTrackParam paramTOF;
    Double_t tofChi2=0, tofNCl=0;            
    //
    // prepare sorted points
    AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
    if (!pointarray) continue;
    Int_t nPointsAll = pointarray->GetNPoints();  // # space points of all available detectors
    SortPointArray(pointarray, sortedIndex);  // space point indices sorted by radius in increasing order

    // 2.) ITS, TRD and ITS-TRD refits
    //
    if (!RefitITStrack(friendTrack,mass,paramITS, itsChi2, itsNCl, sortedIndex)) continue;
    if (itsNCl<4) continue; 
    //
    RefitTRDtrack(friendTrack,mass,paramTRD, trdChi2, trdNCl, sortedIndex); 
    paramTOF=paramITS;
    RefitTOFtrack(friendTrack,mass,paramTOF, tofChi2, tofNCl, sortedIndex);
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("Refitting",fTrackCounter,1,0,0);
    if ((trdNCl+tofNCl+itsNCl)<5) continue; // - use ITS only tracks also  -should it be option?
    //
    // 3.) Propagate to TPC volume, histogram residuals to TPC clusters and dump all information to TTree
    //
    AliTrackPoint spacepoint;  
    Int_t volId,modId,layerId;      
    fTrackCounter++; // increase total track number
    //
    // 3.a) Make a local copy of clusters and apply transformation
    //
    //
    for (Int_t iPoint=kMaxLayer;iPoint--;){
      //
      // reset track interpolation statuses
      trackArrayITS[iPoint].SetUniqueID(0);
      trackArrayTRD[iPoint].SetUniqueID(0);
      trackArrayITSTRD[iPoint].SetUniqueID(0);
      trackArrayTOF[iPoint].SetUniqueID(0);
      trackArrayITSTOF[iPoint].SetUniqueID(0);
      //
      const AliTPCclusterMI *cluster=seed->GetClusterPointer(iPoint);
      if (!cluster){
	clusterArray[iPoint].SetVolumeId(0);
	continue;
      }
      clusterArray[iPoint]=*cluster;
      Int_t i[1]={cluster->GetDetector()};
      Double_t x[3]={static_cast<Double_t>(cluster->GetRow()),cluster->GetPad(),cluster->GetTimeBin()};
      transform->Transform(x,i,0,1);
      clusterArray[iPoint].SetX(x[0]);
      clusterArray[iPoint].SetY(x[1]);
      clusterArray[iPoint].SetZ(x[2]);
    }
    transform->GetCurrentRecoParamNonConst()->SetUseComposedCorrection(backupUseComposedCorrection);
    //
    // 4.) Propagate  ITS tracks outward
    // 
    Bool_t itsOK=kTRUE;
    int npUpdITS = 0;
    for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++) {
      //trackArrayITS[iPoint].SetUniqueID(0);
      AliTPCclusterMI &cluster=clusterArray[iPoint];
      if (cluster.GetVolumeId()==0) continue;
      Float_t fxyz[3] = {0};
      cluster.GetGlobalXYZ(fxyz);
      Double_t xyz[3]={fxyz[0],fxyz[1],fxyz[2]};
      Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
      paramITS.Global2LocalPosition(xyz,alpha);	
      if (!(itsOK=paramITS.Rotate(alpha))) break;
      // full material correction makes sense only when crossing the boundary of the TPC
      itsOK = (++npUpdITS)==1 ? 
	AliTrackerBase::PropagateTrackToBxByBz(&paramITS,xyz[0],mass,1,kFALSE) :
	PropagateInTPCTo(&paramITS,xyz[0],fRhoTPC,fX0TPC,mass) &&
	TMath::Abs(paramITS.GetSnp())<kMaxSNP;
      if (itsOK){
	trackArrayITS[iPoint]=paramITS;
	trackArrayITS[iPoint].SetUniqueID(1);
      }
      else break; // no sense to propagate farther
    }
    if (!itsOK) continue; // no point in continuing if ITS failed
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("ExtrapolateITS",fTrackCounter,2,0,0);  
    //
    // 5.) Propagate  TRD/TOF tracks inwards
    //
    Bool_t trdOK=(trdNCl>0);
    Bool_t tofOK=(tofNCl>0);
    int npUpdTRD = 0, npUpdTOF = 0;
    //
    for (Int_t iPoint=kMaxLayer;iPoint--;){  
      if (!trdOK && !tofOK) break; // no sense to continue;
      AliTPCclusterMI &cluster=clusterArray[iPoint];
      //      if (cluster==NULL) continue;
      if (cluster.GetVolumeId()==0) continue;
      Float_t fxyz[3] = {0};
      cluster.GetGlobalXYZ(fxyz);
      Double_t alpha = TMath::ATan2(fxyz[1],fxyz[0]);            

      if (trdOK){
	Double_t xyz[3]={fxyz[0],fxyz[1],fxyz[2]};
	paramTRD.Global2LocalPosition(xyz,alpha);	
	// material correction makes sense only when crossing the boundary of the TPC
	trdOK = paramTRD.Rotate(alpha) && ((++npUpdTRD)==1 ? 
					   AliTrackerBase::PropagateTrackToBxByBz(&paramTRD,xyz[0],mass,1,kFALSE) :
					   PropagateInTPCTo(&paramTRD,xyz[0],fRhoTPC,fX0TPC,mass)) 
	  &&                               TMath::Abs(paramTRD.GetSnp())<kMaxSNP;
	if (trdOK){
	  trackArrayTRD[iPoint]=paramTRD;
	  trackArrayTRD[iPoint].SetUniqueID(1);
	  //
	  trackArrayITSTRD[iPoint]=paramTRD;
	  if (trackArrayITS[iPoint].GetUniqueID()) { // update by ITS only if the latter is OK
	    AliTrackerBase::UpdateTrack(trackArrayITSTRD[iPoint], trackArrayITS[iPoint]);
	    Double_t chi2=trackArrayITSTRD[iPoint].GetY()-trackArrayITS[iPoint].GetY();
	    chi2*=chi2;
	    chi2/=trackArrayITSTRD[iPoint].GetSigmaY2()+trackArrayITS[iPoint].GetSigmaY2()+kMaxAlignTolerance;
	    if (chi2<kMaxChi2) trackArrayITSTRD[iPoint].SetUniqueID(1);
	  }
	}
      }
      if (tofOK){
	Double_t xyz[3]={fxyz[0],fxyz[1],fxyz[2]};
	paramTOF.Global2LocalPosition(xyz,alpha);	
	// material correction makes sense only when crossing the boundary of the TPC
	tofOK = paramTOF.Rotate(alpha) && ((++npUpdTOF)==1 ?
					   AliTrackerBase::PropagateTrackToBxByBz(&paramTOF,xyz[0],mass,1,kFALSE) :
					   PropagateInTPCTo(&paramTOF,xyz[0],fRhoTPC,fX0TPC,mass))
	  &&                               TMath::Abs(paramTOF.GetSnp())<kMaxSNP;
	if (tofOK){
	  trackArrayTOF[iPoint]=paramTOF;
	  trackArrayTOF[iPoint].SetUniqueID(1);

	  trackArrayITSTOF[iPoint]=paramTOF;
	  if (trackArrayITS[iPoint].GetUniqueID()) {  // update by ITS only if the latter is OK
	    AliTrackerBase::UpdateTrack(trackArrayITSTOF[iPoint], trackArrayITS[iPoint]);
	    Double_t chi2=trackArrayITSTOF[iPoint].GetY()-trackArrayITS[iPoint].GetY();
	    chi2*=chi2;
	    chi2/=trackArrayITSTOF[iPoint].GetSigmaY2()+trackArrayITS[iPoint].GetSigmaY2()+kMaxAlignTolerance;
	    if (chi2<kMaxChi2)  trackArrayITSTOF[iPoint].SetUniqueID(1);
	  }
	}
      }
      //
    }
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("InterpolateTRD",fTrackCounter,3,0,0);  

    if ( ((fStreamLevelTrack&kStremInterpolation)>0)&&(fTrackCounter%fSyswatchStep==0)){
      //if ((fTrackCounter%fSyswatchStep==0)){
      for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++){
	AliTPCclusterMI &cluster=clusterArray[iPoint];
	//if (cluster==NULL) continue;
	if (cluster.GetVolumeId()==0) continue;

	(*fStreamer)<<"interpolation"<<
          "itrack="<<fTrackCounter<<  // total track #
          "cluster.="<<&cluster<<  // space points                                    //
          "itsNCl="<<itsNCl<<
          "trdNCl="<<trdNCl<<
          "tofNCl="<<tofNCl<<
	  "itsOK="<<itsOK<<
	  "trdOK="<<trdOK<<
	  "tofOK="<<tofOK<<
          //
          "itsChi2="<<itsChi2<<
          "trdChi2="<<trdChi2<<
          "tofChi2="<<tofChi2<<
	  "tofBC="<<tofDiff<<
          //
          "trackITS.="<<&trackArrayITS[iPoint]<<  // ITS fit
          "trackTRD.="<<&trackArrayTRD[iPoint]<<  // TRD fit
          "trackTOF.="<<&trackArrayTOF[iPoint]<<  // TOF fit
          "trackITSTRD.="<<&trackArrayITSTRD[iPoint]<<  // ITS-TRD fit
          "trackITSTOF.="<<&trackArrayITSTOF[iPoint]<<  // ITS-TOF fit
          "\n";	
      }
    }
    UShort_t counter=0;
    Double_t rounding=200;
    //    
    memset( deltaITS0.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaITS1.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTRD0.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTRD1.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTOF0.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTOF1.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    //
    memset( vecR.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( vecPhi.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( vecZ.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( vecSec.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    //
    for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++){      
      //
      deltaITS0[counter]=deltaITS1[counter]=deltaTRD0[counter]=deltaTRD1[counter]=deltaTOF0[counter]=deltaTOF1[counter]=-999;
      vecR[counter] = -999;
      //
      AliTPCclusterMI &cluster=clusterArray[iPoint];
      if (cluster.GetVolumeId()==0) continue;
      Double_t   zsignSector=((cluster.GetDetector()%36)<18) ? 1.:-1.;
      //if (zsignSector*cluster.GetZ()<0.) continue;
      //
      if (trackArrayITS[iPoint].GetUniqueID()>0) { // deltas make sense only if ITS was ok
	deltaITS0[counter]= TMath::Nint(trackArrayITS[iPoint].GetY()*rounding)/rounding;
	deltaITS1[counter]= TMath::Nint((trackArrayITS[iPoint].GetZ()-cluster.GetZ())*rounding)/rounding;
	//
	if (trackArrayITSTRD[iPoint].GetUniqueID()>0){
	  deltaTRD0[counter]= TMath::Nint(trackArrayITSTRD[iPoint].GetY()*rounding)/rounding;
	  deltaTRD1[counter]= TMath::Nint((trackArrayITSTRD[iPoint].GetZ()-cluster.GetZ())*rounding)/rounding;
	}
	if (trackArrayITSTOF[iPoint].GetUniqueID()>0){
	  deltaTOF0[counter]= TMath::Nint(trackArrayITSTOF[iPoint].GetY()*rounding)/rounding;
	  deltaTOF1[counter]= TMath::Nint((trackArrayITSTOF[iPoint].GetZ()-cluster.GetZ())*rounding)/rounding;
	}
	// vecR(kMaxLayer), vecPhi(kMaxLayer), vecZ(kMaxLayer);
	vecR[counter]=trackArrayITS[iPoint].GetX();
	vecPhi[counter]=trackArrayITS[iPoint].GetAlpha();
	vecZ[counter]=trackArrayITS[iPoint].GetZ();
	vecSec[counter]=cluster.GetDetector();
	counter++;
      }
    }
    AliExternalTrackParam * ip = (AliExternalTrackParam *)esdTrack->GetInnerParam();
    Int_t timeStamp= esdEvent->GetTimeStamp();
    (*fStreamer)<<"delta"<<
      "nTracks="<<nTracks<<               // number of tracks in event (pileup indicator)
      "nPrimTracks="<<nPrimTracks<<       // number of tracks pointed to primary vertes of selected event
      "timeStamp="<<timeStamp<<           // time stamp
      "itrack="<<fTrackCounter<<          // total track #
      "itsNCl="<<itsNCl<<
      "trdNCl="<<trdNCl<<
      "tofNCl="<<tofNCl<<
      "itsOK="<<itsOK<<
      "trdOK="<<trdOK<<
      "tofOK="<<tofOK<<
      "itsChi2="<<itsChi2<<
      "trdChi2="<<trdChi2<<
      "tofChi2="<<tofChi2<<
      "tofBC="<<tofDiff<<
      //
      "track.="<<ip<<                    // track parameters at inner wal of TPC
      "npValid="<<counter<<
      "vecR.="<<&vecR<<          
      "vecPhi.="<<&vecPhi<<
      "vecSec.="<<&vecSec<<              // sector number
      "vecZ.="<<&vecZ<<
      "its0.="<<&deltaITS0<<
      "its1.="<<&deltaITS1<<
      "trd0.="<<&deltaTRD0<<
      "trd1.="<<&deltaTRD1<<
      "tof0.="<<&deltaTOF0<<
      "tof1.="<<&deltaTOF1<<
      "\n";    
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("FittTree",fTrackCounter,4,0,0);  
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("FillHistos",fTrackCounter,5,0,0);  
  }
  transform->GetCurrentRecoParamNonConst()->SetUseComposedCorrection( backupUseComposedCorrection);
  //
 // end of track loop
}

void AliTPCcalibAlignInterpolation::CreateResidualHistosInterpolation(Double_t dy, Double_t dz, Int_t selHis){
  //
  // Make cluster residual histograms
  //
  Double_t xminTrack[9], xmaxTrack[9];
  Double_t xminTrackITS[9], xmaxTrackITS[9];
  Int_t    binsTrack[9], binsTrackITS[9];
  TString  axisName[9],axisTitle[9];
  //
  // 0 - local   q/pt
  // 1 - global  phi in sector number  as float
  // 2 - local   r
  // 3 - local   kz
  // 4 - delta   of interest

  // 
  // gx,gy,gz - will be taken from the TPC
  //
  //
  axisName[0]="qpt";    axisTitle[0]="q/pt (c/GeV)";                         // to fill : track.GetSigned1Pt() 
  binsTrack[0]=5;        xminTrack[0]=-2.5;        xmaxTrack[0]=2.5; 
  binsTrackITS[0]=5;     xminTrackITS[0]=-2.5;     xmaxTrackITS[0]=2.5; 

  //
  axisName[1]="sector";  axisTitle[1]="Sector Number";              // to fill:   9*atan2(gy,gx)/pi+ if (sector>0) sector+18
  binsTrack[1]=180;      xminTrack[1]=0;           xmaxTrack[1]=18; 
  binsTrackITS[1]=180;   xminTrackITS[1]=0;        xmaxTrackITS[1]=18; 
  //
  axisName[2]="R";       axisTitle[2]="r (cm)";                          // to fill:    gr=sqrt(gy**2+gx**2)
  binsTrack[2]=53;       xminTrack[2]=85.;         xmaxTrack[2]=245.; 
  binsTrackITS[2]=53;    xminTrackITS[2]=85.;      xmaxTrackITS[2]=245.; 
  //
  //
  axisName[3]="kZ";      axisTitle[3]="z/r";                          // to fill : gz/gr 
  binsTrack[3]=20;       xminTrack[3]=-1.0;        xmaxTrack[3]=1.0;  // +-1 for ITS+TRD and ITS+TOF 
  binsTrackITS[3]=20;    xminTrackITS[3]=-1.8;     xmaxTrackITS[3]=1.8;  // +-1.8 for the ITS 
  //
  axisName[4]="delta";   axisTitle[4]="#Delta (cm)";                 // to fill    local(clusterY-track.y)
  binsTrack[4]=100;       xminTrack[4]=-dy;        xmaxTrack[4]=dy; 
  binsTrackITS[4]=100;    xminTrackITS[4]=-dy;     xmaxTrackITS[4]=dy; 

  // 
  binsTrack[4]=TMath::Min(Int_t(20.+2.*dy/0.05),120); // buffer should be smaller than 1 GBy
  if (selHis==0 ||selHis<0) fHisITSDRPhi = new THnF("deltaRPhiTPCITS","#Delta_{Y} (cm)", 5, binsTrackITS,xminTrackITS, xmaxTrackITS);
  if (selHis==1 ||selHis<0) fHisITSTRDDRPhi = new THnF("deltaRPhiTPCITSTRD","#Delta_{Y} (cm) TPC-(ITS+TRD)", 5, binsTrack,xminTrack, xmaxTrack);
  if (selHis==2 ||selHis<0) fHisITSTOFDRPhi = new THnF("deltaRPhiTPCITSTOF","#Delta_{Y} (cm) TPC-(ITS+TOF)", 5, binsTrack,xminTrack, xmaxTrack);
  //
  binsTrack[4]=TMath::Min(Int_t(20.+2.*dz/0.05),120); // buffer should be smaller than 1 GBy
  if (selHis==3 ||selHis<0) fHisITSDZ = new THnF("deltaZTPCITS","#Delta_{Z} (cm)", 5, binsTrackITS,xminTrackITS, xmaxTrackITS);
  if (selHis==4 ||selHis<0) fHisITSTRDDZ = new THnF("deltaZTPCITSTRD","#Delta_{Z} (cm) TPC-(ITS+TRD)", 5, binsTrack,xminTrack, xmaxTrack);
  if (selHis==5 ||selHis<0) fHisITSTOFDZ = new THnF("deltaZTPCITSTOF","#Delta_{Z} (cm) TPC-(ITS+TOF)", 5, binsTrack,xminTrack, xmaxTrack);
  //
  //
  //
  THn *hisToFill[6]={GetHisITSDRPhi(), GetHisITSTRDDRPhi(), GetHisITSTOFDRPhi(), GetHisITSDZ(), GetHisITSTRDDZ(), GetHisITSTOFDZ()};
  for (Int_t ihis=0; ihis<6; ihis++){
    if (hisToFill[ihis]) for (Int_t ivar2=0;ivar2<5;ivar2++){ 
      hisToFill[ihis]->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
      hisToFill[ihis]->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());      
    }
  }

}



void  AliTPCcalibAlignInterpolation::CreateDistortionMapsFromFile(const char * inputFile, const char *outputFile){
  //
  // Create distortion maps from residual histograms
  // TPC cluster to ITS, ITS-TRD and ITS-TOF track fits
  //
  TFile *fHistos  = TFile::Open(inputFile);
  
  THnF *histoITS = (THnF*) fHistos->Get("deltaRPhiTPCITS");
  THnF *histoITSTRD= (THnF*) fHistos->Get("deltaRPhiTPCITSTRD");
  THnF *histoITSTOF = (THnF*) fHistos->Get("deltaRPhiTPCITSTOF");
  THnF *histoITSZ = (THnF*) fHistos->Get("deltaZTPCITS");
  THnF *histoITSTRDZ= (THnF*) fHistos->Get("deltaZTPCITSTRD");
  THnF *histoITSTOFZ = (THnF*) fHistos->Get("deltaZTPCITSTOF");
  
  TTreeSRedirector * pcstream = new TTreeSRedirector(outputFile,"recreate");
  
  TMatrixD projectionInfo(5,5);
  projectionInfo(0,0)=0;  projectionInfo(0,1)=0;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=4;  projectionInfo(1,1)=0;  projectionInfo(1,2)=1; 
  projectionInfo(2,0)=3;  projectionInfo(2,1)=0;  projectionInfo(2,2)=1;
  projectionInfo(3,0)=2;  projectionInfo(3,1)=0;  projectionInfo(3,2)=1;
  projectionInfo(4,0)=1;  projectionInfo(4,1)=0;  projectionInfo(4,2)=1;
  
  TStatToolkit::MakeDistortionMap(4, histoITS,    pcstream, projectionInfo); 
  TStatToolkit::MakeDistortionMap(4, histoITSTRD, pcstream, projectionInfo); 
  TStatToolkit::MakeDistortionMap(4, histoITSTOF, pcstream, projectionInfo); 
  TStatToolkit::MakeDistortionMap(4, histoITSZ,    pcstream, projectionInfo); 
  TStatToolkit::MakeDistortionMap(4, histoITSTRDZ, pcstream, projectionInfo); 
  TStatToolkit::MakeDistortionMap(4, histoITSTOFZ, pcstream, projectionInfo); 
  delete pcstream;
  //
}

void    AliTPCcalibAlignInterpolation::FillHistogramsFromChain(const char * residualList, Double_t dy, Double_t dz, Int_t startTime, Int_t stopTime, Int_t maxStat, Int_t selHis,const char * residualInfoFile ){
  /**
   * Trees with point-track residuals to residual histogram
   * @param residualList  text file with tree list
   * Output - ResidualHistograms.root file with hitogram within AliTPCcalibAlignInterpolation object
   */
  //
  //
  // 
  ::Info(" AliTPCcalibAlignInterpolation::FillHistogramsFromChain","Start %s\n", residualList);
  //
  // 0.) Load current information file and bookd variables
  // 
  const Int_t nSec=81;         // 72 sector +5 sumarry info+ 4 medians +
  TVectorF meanNcl(nSec);      // mean current estator ncl per sector
  TVectorF meanNclUsed(nSec);  // mean current estator ncl per sector
  Double_t meanTime=0, maxTime=startTime, minTime=stopTime;
  Int_t currentTrack=0;  

  TFile *finfo = TFile::Open(residualInfoFile);
  TTree *treeInfo=0;
  if (finfo) treeInfo=(TTree*)finfo->Get("sumaryTime"); 
  TGraphErrors * nclArray[nSec]={0};
  TGraphErrors * nclArrayUsed[nSec]={0};
  
  if (treeInfo) {
    for (Int_t iSec=0; iSec<nSec; iSec++){
      nclArray[iSec]=0;
      nclArrayUsed[iSec]=0;
      treeInfo->SetBranchAddress(TString::Format("grNcl%d.",iSec).Data(),&nclArray[iSec]);
      treeInfo->SetBranchAddress(TString::Format("grNclUsed%d.",iSec).Data(),&nclArrayUsed[iSec]);
    }
    treeInfo->GetEntry(0);
  }
  //
  // 1.) Fill histograms and mean informations
  //
  const Int_t knPoints=kMaxRow;
  AliTPCcalibAlignInterpolation * calibInterpolation = new  AliTPCcalibAlignInterpolation("calibInterpolation","calibInterpolation",kFALSE);
  calibInterpolation->CreateResidualHistosInterpolation(dy,dz,selHis);
  TString branches[6]={"its0.","trd0.","tof0.", "its1.","trd1.","tof1."};
  //
  TVectorF *vecDelta= 0;
  TVectorF *vecR=0;
  TVectorF *vecSec=0;
  TVectorF *vecPhi=0;
  TVectorF *vecZ=0;
  Int_t timeStamp=0;
  AliExternalTrackParam *param = 0;
  //
  TString  esdList0 = gSystem->GetFromPipe(TString::Format("cat %s",residualList).Data());
  TObjArray *esdArray= esdList0.Tokenize("\n");  
  Int_t nesd = esdArray->GetEntriesFast();  
  //
  THn *hisToFill[6]={calibInterpolation->GetHisITSDRPhi(), calibInterpolation->GetHisITSTRDDRPhi(),  calibInterpolation->GetHisITSTOFDRPhi(), calibInterpolation->GetHisITSDZ(), calibInterpolation->GetHisITSTRDDZ(),  calibInterpolation->GetHisITSTOFDZ()};
  TTreeSRedirector * fout = 0;
  if (selHis<0)  {
    if (startTime<=0) fout=new TTreeSRedirector("ResidualHistograms.root","recreate");
    if (startTime>0) fout=new TTreeSRedirector(TString::Format("ResidualHistograms_Time%d.root",startTime).Data(),"recreate");
  }
  if (selHis>=0) {
    if (startTime<=0)  fout=new TTreeSRedirector(TString::Format("ResidualHistograms_His%d.root",selHis).Data(),"recreate");
    if (startTime>0)   fout=new TTreeSRedirector(TString::Format("ResidualHistograms_His%d_Time%d.root",selHis,startTime).Data(),"recreate");
  }
  TH1 * hisTime=0;
  if (startTime>0) hisTime=new TH1F("hisTrackTime","hisTrackTime",(stopTime-startTime)/20,startTime,stopTime);
  TStopwatch timerAll;
  UShort_t npValid=knPoints;
  for (Int_t ihis=0; ihis<6; ihis++){    
    if (selHis>=0 && ihis!=selHis) continue;
    for (Int_t iesd=0; iesd<nesd; iesd++){
      TStopwatch timerFile;
      TFile *esdFile = TFile::Open(esdArray->At(iesd)->GetName(),"read");
      if (!esdFile) continue;
      TTree *tree = (TTree*)esdFile->Get("delta");
      tree->SetBranchStatus("*",kFALSE);
      if (!tree) continue;
      ::Info(" AliTPCcalibAlignInterpolation::FillHistogramsFromChain", "Processing file \t %s\n",esdArray->At(iesd)->GetName());
      AliSysInfo::AddStamp(esdArray->At(iesd)->GetName(),ihis,iesd,currentTrack);
      tree->SetBranchStatus("timeStamp",kTRUE);
      TBranch *br = tree->GetBranch("timeStamp");
      tree->SetBranchStatus("vecR.",kTRUE);
      tree->SetBranchStatus("vecSec.",kTRUE);
      tree->SetBranchStatus("vecPhi.",kTRUE);
      tree->SetBranchStatus("vecZ.",kTRUE);
      tree->SetBranchStatus("track.",kTRUE);      
      tree->SetBranchAddress("vecR.",&vecR);
      tree->SetBranchAddress("vecSec.",&vecSec);
      tree->SetBranchAddress("vecPhi.",&vecPhi);
      tree->SetBranchAddress("vecZ.",&vecZ);
      tree->SetBranchAddress("track.",&param);
      br->SetAddress(&timeStamp);
      if (tree->GetBranch("npValid")!=NULL) {
	tree->SetBranchStatus("npValid",kTRUE);
	tree->SetBranchAddress("npValid",&npValid);
      }
      tree->SetBranchStatus(branches[ihis],kTRUE);
      tree->SetBranchAddress(branches[ihis],&vecDelta);
      
      Int_t ntracks=tree->GetEntries();
      //
      for (Int_t itrack=0; itrack<ntracks; itrack++){
	if (startTime>0){
	  br->GetEntry(itrack);
	  if (timeStamp<startTime  || timeStamp>stopTime) continue;
	  hisTime->Fill(timeStamp);
	}
	tree->GetEntry(itrack);
	const Float_t *vSec= vecSec->GetMatrixArray();
	const Float_t *vPhi= vecPhi->GetMatrixArray();
	const Float_t *vR  = vecR->GetMatrixArray();
	const Float_t *vZ  = vecZ->GetMatrixArray();
	const Float_t *vDelta  = vecDelta->GetMatrixArray();
	//
	currentTrack++;
	if (timeStamp<minTime) minTime=0;
	if (timeStamp>maxTime) maxTime=0;
	meanTime+=timeStamp;
	if (treeInfo) for (Int_t iSec=0; iSec<nSec; iSec++){
	  meanNcl[iSec]+=nclArray[iSec]->Eval(timeStamp);
	  meanNclUsed[iSec]+=nclArrayUsed[iSec]->Eval(timeStamp);
	}

	if (maxStat>0 &&currentTrack>maxStat) break;
	//for (Int_t ipoint=0; ipoint<knPoints; ipoint++){
	for (Int_t ipoint=0; ipoint<npValid; ipoint++){
	  if (vR[ipoint]<=0 || vDelta[ipoint]<-990.) continue;
	  Double_t sector=9.*vPhi[ipoint]/TMath::Pi();
	  if (sector<0) sector+=18;
	  Double_t deltaPhi=vPhi[ipoint]-TMath::Pi()*(sector+0.5)/9.;
	  Double_t localX = TMath::Cos(deltaPhi)*vR[ipoint];
	  Double_t xxx[5]={ param->GetParameter()[4], sector, localX,   vZ[ipoint]/localX, vDelta[ipoint]};
	  if (xxx[4]==0) continue;
	  Double_t side=-1.+2.*((TMath::Nint(vSec[ipoint])%36)<18);
	  if ((vZ[ipoint]*side)<-1) xxx[3]=side*0.001; // do not mix z on A side and C side 	  
	  hisToFill[ihis]->Fill(xxx);	  
	}
      }
      timerFile.Print();
      delete tree;
      delete esdFile;
      
    }    
    fout->GetFile()->cd();
    hisToFill[ihis]->Write();
  }
  if (hisTime) hisTime->Write();
  ::Info(" AliTPCcalibAlignInterpolation::FillHistogramsFromChain","End of processing\n");
  timerAll.Print();
  //
  // 2.) Fill metadata information
  //
  if (currentTrack>0){
    meanTime/=currentTrack;
    if (treeInfo) for (Int_t iSec=0; iSec<nSec; iSec++){
      meanNcl[iSec]/=currentTrack;
      meanNclUsed[iSec]/=currentTrack;
    }
  }
  (*fout)<<"metaData"<<
    "startTime="<<startTime<<        // start time  as requested
    "stopTime="<<stopTime<<          // stop time as requested
    "meanTime="<<meanTime<<          // mean time 
    "minTime="<<minTime<<            // minimal time stamp in data sample
    "maxTime="<<maxTime<<            // maximal time stamp in data sample
    "ntracksUsed="<<currentTrack<<   // number of tracks acumulated in time interval
    "meanNcl.="<<&meanNcl<<          // current estimator - mean number of clusters
    "meanNclUsed.="<<&meanNclUsed;   // current estimator - mean number of clusters
  
  for (Int_t iSec=0; iSec<nSec; iSec++){
    (*fout)<<"metaData"<<
      TString::Format("grNcl%d.=",iSec).Data()<< nclArray[iSec]<<
      TString::Format("grNclUsed%d.=",iSec).Data()<< nclArrayUsed[iSec];
  }
  (*fout)<<"metaData"<<"\n";

  delete fout;
}


void     AliTPCcalibAlignInterpolation::FillHistogramsFromStreamers(const char * residualList, Double_t dy, Double_t dz, Int_t downscale){
  /**
   * Input list of ErrParam trees as defined in the AliTPCtracker in debug mode 
   * @param residualList  text file with tree list
   * Output - ResidualHistograms.root file with hitogram within AliTPCcalibAlignInterpolation object
   residualList="residual.list"
   dy=1; dz=1
   */
  //
  //
  // 
  AliTPCcalibAlignInterpolation * calibInterpolation = new  AliTPCcalibAlignInterpolation("calibInterpolation","calibInterpolation",kFALSE);
  calibInterpolation->CreateResidualHistosInterpolation(dy,dz);
  TString  esdList0 = gSystem->GetFromPipe(TString::Format("cat %s",residualList).Data());
  TObjArray *esdArray= esdList0.Tokenize("\n");  
  Int_t nesd = esdArray->GetEntriesFast();  
  //
  THn *hisToFill[6]={calibInterpolation->GetHisITSDRPhi(), calibInterpolation->GetHisITSTRDDRPhi(),  calibInterpolation->GetHisITSTOFDRPhi(), calibInterpolation->GetHisITSDZ(), calibInterpolation->GetHisITSTRDDZ(),  calibInterpolation->GetHisITSTOFDZ()};
  //
  //
  AliExternalTrackParam * param=0;
  AliTPCclusterMI * cl=0;
  Int_t iter=0;
  Int_t currentCl=0;
  for (Int_t iesd=0; iesd<nesd; iesd++){
    TFile *esdFile = TFile::Open(esdArray->At(iesd)->GetName(),"read");
    if (!esdFile) continue;
    TTree *tree = (TTree*)esdFile->Get("ErrParam"); 
    if (!tree) continue;
    tree->SetBranchAddress("Cl.",&cl);
    tree->SetBranchAddress("T.",&param);    
    tree->SetBranchAddress("iter",&iter);    
    Int_t nCl=tree->GetEntries();
    for (Int_t iCl=0; iCl<nCl; iCl+=downscale){
      tree->GetEntry(iCl);
      if (iCl%100000==0) printf("%d\n",iCl);
      currentCl++;
      Float_t xyz[3]={0};
      cl->GetGlobalXYZ(xyz);
      Double_t phi = TMath::ATan2(xyz[1],xyz[0]);
      Double_t radius=TMath::Sqrt(xyz[1]*xyz[1]+xyz[0]*xyz[0]);
      param->Rotate(phi);
      param->PropagateTo(radius,0.); // for big distortion we should query field, for small deltas we are using straight approximtion 
      Double_t sector=9*phi/TMath::Pi();
      if (sector<0) sector+=18;
      Double_t deltaY=param->GetY();
      Double_t deltaZ=param->GetZ()-cl->GetZ();
      Double_t localX = cl->GetX();
      Double_t   zsignSector=((cl->GetDetector()%36)<18) ? 1.:-1.;
      if (zsignSector*cl->GetZ()<0.) continue;
      Double_t xxx[5]={ param->GetParameter()[4], sector, localX,   cl->GetZ()/cl->GetX(),  deltaY};
      hisToFill[iter]->Fill(xxx);	  
      xxx[4]=deltaZ;
      hisToFill[3+iter]->Fill(xxx);	  
    }
  }
  TFile * fout = TFile::Open("ResidualHistograms.root","recreate");
  calibInterpolation->GetHisITSDRPhi()->Write("deltaYIter0");
  calibInterpolation->GetHisITSTRDDRPhi()->Write("deltaYIter1");
  calibInterpolation->GetHisITSTOFDRPhi()->Write("deltaYIter2");
  calibInterpolation->GetHisITSDZ()->Write("deltaZIter0");
  calibInterpolation->GetHisITSTRDDZ()->Write("deltaZIter1");
  calibInterpolation->GetHisITSTOFDZ()->Write("deltaZIter2");
  delete fout;
}




TTree*  AliTPCcalibAlignInterpolation::AddFriendDistortionTree(TTree * tree, const char * fname,  const char *treeName, const char *friendAlias){
  //
  //
  //
  TFile * fin = TFile::Open(fname);
  if (fin==NULL) {
    ::Error("AliTPCcalibAlignInterpolation::AddFriendDistotionTree",TString::Format("file %s not readable", fname).Data());
    return 0;
  }
  TTree * treeFriend = (TTree*) fin->Get(treeName);
  
  if (treeFriend==NULL){
    ::Error("AliTPCcalibAlignInterpolation::AddFriendDistotionTree",TString::Format("file %s not readable", fname).Data());
    return 0;
  }
  if (tree==NULL) {
    tree = treeFriend;
  }else{
    tree->AddFriend(treeFriend,TString::Format("%s",friendAlias).Data());
    tree->SetAlias(TString::Format("%sOK",friendAlias).Data(),TString::Format("%s.rms>0&&abs(%s.mean-%s.meanG)<2&&%s.chi2G>0&&%s.rmsG<2&&%s.rmsG/%s.rms<2",friendAlias,friendAlias,friendAlias,friendAlias,friendAlias,friendAlias,friendAlias).Data());
    tree->SetAlias(TString::Format("%sDrawOK",friendAlias).Data(),TString::Format("%s.rms>0&&abs(%s.mean-%s.meanG)<4&&%s.chi2G>0",friendAlias,friendAlias,friendAlias,friendAlias).Data()); 
  }
  return tree;
}

//_____________________________________________________________________________
Bool_t AliTPCcalibAlignInterpolation::PropagateInTPCTo(AliExternalTrackParam* t, Double_t xk, Double_t rho,Double_t x0, Double_t mass) 
{
  //-----------------------------------------------------------------
  //  This function propagates a track to a reference plane x=xk.
  //  rho - density of the crossed matrial (g/cm^3)
  //  x0  - radiation length of the crossed material (g/cm^2) 
  //-----------------------------------------------------------------
  //
  Double_t old[3]={t->GetX(),t->GetY(),t->GetZ()};
  Double_t b[3]; AliTrackerBase::GetBxByBz(old,b);
  if (!t->PropagateToBxByBz(xk,b)) return kFALSE;

  Double_t d = TMath::Sqrt((t->GetX()-old[0])*(t->GetX()-old[0]) + 
                           (t->GetY()-old[1])*(t->GetY()-old[1]) + 
                           (t->GetZ()-old[2])*(t->GetZ()-old[2]));
  if (old[0] < xk) d = -d;
  if (!t->CorrectForMeanMaterial(d*rho/x0,d*rho,mass,
				 kFALSE,AliExternalTrackParam::BetheBlochGas)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliTPCcalibAlignInterpolation::ExtractTPCGasData()
{
  // get TPC gas rho and X0
  double p0[3] = {90,1,45};
  double p1[3] = {240,1,120};
  double par[10];
  AliTrackerBase::MeanMaterialBudget(p0,p1,par);
  fRhoTPC = par[0]>0 ? par[0] : 0.9e-3;
  double l = par[4];
  fX0TPC  = par[1]>0 ? par[4]/par[1] : 28.94;
  //
  AliInfoF("Propagation in TPC will use rho=%.2f X0=%.2f",fRhoTPC,fX0TPC);
}


void AliTPCcalibAlignInterpolation::MakeEventStatInfo(const char * inputList, Int_t timeInterval, Int_t id, Int_t skip){
  //
  /// Code to query statistical event information from the ResidualTrees.root file 
  /// output written to file residualInfo.root
  ///   \param const char * inputList - ascii file with input list
  ///   \param Int_t timeInterval     - length of time interval (beginning of time intervals rounded)
  ///   \param id                     - additional ID added to the tree
  ///   \param skip                   - parameter skip file
  /// Algorithm:
  ///   1.) Cache information per files - beginTime and endTime for file
  ///   2.) Cache information per time interval

  /*
    run=240204;
    GetResidualStatInfo("cat residual.list",300,run,1);
  */
  TObjArray *array = TString(gSystem->GetFromPipe(TString::Format("%s",inputList).Data())).Tokenize("\n");
  Int_t nFiles=array->GetEntries();
  if (nFiles<=0) {
    ::Error("GetResidualStatInfo. Wrong input list",inputList);
    return;
  }
  TStopwatch timer;
  //
  // 1.) Cache information per files - beginTime and endTime for file
  //
  TStopwatch timer1;
  TTreeSRedirector * pcstream = new TTreeSRedirector("residualInfo.root", "recreate");
  for (Int_t iFile=0; iFile<nFiles; iFile+=skip){
    timer.Start();
    printf("%d\t%s\n",iFile,array->At(iFile)->GetName());
    TFile * f = TFile::Open(array->At(iFile)->GetName());
    if (f==NULL) continue;
    TTree * treeInfo = (TTree*)f->Get("eventInfo");
    if (treeInfo==NULL) continue;
    Int_t entriesInfo=treeInfo->GetEntries();
    Int_t entries=treeInfo->Draw("B1","1","goff");
    Double_t maxTime=TMath::MaxElement(entries,treeInfo->GetV1());
    Double_t minTime=TMath::MinElement(entries,treeInfo->GetV1());
    Double_t meanTime=TMath::Mean(entries,treeInfo->GetV1());
    TObjString fname(array->At(iFile)->GetName());
    (*pcstream)<<"summary1"<<
      "iFile="<<iFile<<
      "fname.="<<&fname<<
      "events="<<entriesInfo<<
      "minTime="<<minTime<<
      "maxTime="<<maxTime<<
      "meanTime="<<meanTime<<
      "\n";
    timer.Print();
  }
  delete pcstream;
  ::Info("GetResidualStatInfo","Total time");
  timer1.Print();
  //
  // 2.) Cache information per time interval
  //
  TStopwatch timer2;
  pcstream = new TTreeSRedirector("residualInfo.root", "update");
  TTree * treeSummary1=(TTree*)(pcstream->GetFile()->Get("summary1"));
  Int_t entries = treeSummary1->Draw("minTime","1","goff");
  Long64_t minTime = TMath::MinElement(entries, treeSummary1->GetV1());
  entries = treeSummary1->Draw("maxTime","1","goff");
  Long64_t maxTime = TMath::MaxElement(entries, treeSummary1->GetV1());
  minTime=timeInterval*(minTime/timeInterval);
  maxTime=timeInterval*(1+(maxTime/timeInterval));
  Int_t nIntervals=(maxTime-minTime)/timeInterval;
  Int_t nIntervalsQA=(maxTime-minTime)/15;
  //
  TH1F  * hisEvent= new TH1F("hisEvent","hisEvent",nIntervalsQA,minTime,maxTime);
  const Int_t nSec=81; // 72 sector +5 sumarry info+ 4 medians
  TProfile * profArrayNcl[nSec]={0};
  TProfile * profArrayNclUsed[nSec]={0};
  TGraphErrors * grArrayNcl[nSec]={0};
  TGraphErrors * grArrayNclUsed[nSec]={0};
  TProfile * profArrayITSNcl[3]={0};
  TGraphErrors * grArrayITSNcl[3]={0};
  
  for (Int_t isec=0; isec<nSec; isec++){
    profArrayNcl[isec]=new TProfile(TString::Format("TPCnclSec%d",isec).Data(), TString::Format("TPCnclSec%d",isec).Data(), nIntervalsQA,minTime,maxTime);
    profArrayNclUsed[isec]=new TProfile(TString::Format("TPCnclUsedSec%d",isec).Data(), TString::Format("TPCnclUsedSec%d",isec).Data(), nIntervalsQA,minTime,maxTime);
  }
   for (Int_t iits=0; iits<3; iits++){
    profArrayITSNcl[iits]=new TProfile(TString::Format("ITSnclSec%d",iits).Data(), TString::Format("ITSnclSec%d",iits).Data(), nIntervalsQA,minTime,maxTime);    
  }

  TVectorF *vecNClTPC=0;
  TVectorF *vecNClTPCused=0;
  Int_t nITS[3]={0};
  Int_t timeStamp=0;
  for (Int_t iFile=0; iFile<nFiles; iFile+=skip){
    timer.Start();
    printf("%d\t%s\n",iFile,array->At(iFile)->GetName());    
    TFile * f = TFile::Open(array->At(iFile)->GetName());
    if (f==NULL) continue;
    TTree * treeInfo = (TTree*)f->Get("eventInfo"); 
    if (treeInfo==NULL) continue;
    treeInfo->SetBranchAddress("vecNClTPC.",&vecNClTPC);
    treeInfo->SetBranchAddress("vecNClTPCused.",&vecNClTPCused);
    treeInfo->SetBranchAddress("nSPD",&nITS[0]);
    treeInfo->SetBranchAddress("nSDD",&nITS[1]);
    treeInfo->SetBranchAddress("nSSD",&nITS[2]);
    Bool_t hasTimeStamp=(treeInfo->GetBranch("timeStamp")!=NULL);
    if (hasTimeStamp) treeInfo->SetBranchAddress("timeStamp",&timeStamp);
    if (!hasTimeStamp) ((TBranch*)(treeInfo->GetListOfBranches()->At(1)))->SetAddress(&timeStamp);
    Int_t treeEntries=treeInfo->GetEntries();
    for (Int_t iEntry=0; iEntry<treeEntries; iEntry++){
      treeInfo->GetEntry(iEntry);
      hisEvent->Fill(timeStamp);
      for (Int_t isec=0; isec<72; isec++){
	profArrayNcl[isec]->Fill(timeStamp, (*vecNClTPC)[isec]);
	profArrayNclUsed[isec]->Fill(timeStamp, (*vecNClTPC)[isec]);
	if (isec<36){
	  if (isec<18) 	profArrayNcl[72]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if (isec>=18) profArrayNcl[73]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if (isec<18) 	profArrayNclUsed[72]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	  if (isec>=18) profArrayNclUsed[73]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	}else{
	  if ((isec%36)<18)  profArrayNcl[74]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if ((isec%36)>=18) profArrayNcl[75]->Fill(timeStamp, (*vecNClTPC)[isec]);
	  if ((isec%36)<18)  profArrayNclUsed[74]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	  if ((isec%36)>=18) profArrayNclUsed[75]->Fill(timeStamp, (*vecNClTPCused)[isec]);
	}
	profArrayNcl[76]->Fill(timeStamp, (*vecNClTPC)[isec]);
	profArrayNclUsed[76]->Fill(timeStamp, (*vecNClTPCused)[isec]);
      }
      profArrayNcl[77]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[0])));
      profArrayNcl[78]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[18])));
      profArrayNcl[79]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[36])));
      profArrayNcl[80]->Fill(timeStamp, TMath::Median(18, &((vecNClTPC->GetMatrixArray())[54])));
      //
      profArrayNclUsed[77]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[0])));
      profArrayNclUsed[78]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[18])));
      profArrayNclUsed[79]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[36])));
      profArrayNclUsed[80]->Fill(timeStamp, TMath::Median(18, &((vecNClTPCused->GetMatrixArray())[54])));
      for (Int_t iits=0; iits<3; iits++){
	profArrayITSNcl[iits]->Fill(timeStamp,nITS[iits]);
      }
    }
    timer.Print();
  }
  timer2.Print();
  TGraphErrors grEvent(hisEvent);
  (*pcstream)<<"sumaryTime"<<
    "id="<<id<<
    "grEvent.="<<&grEvent;
  for (Int_t isec=0; isec<nSec; isec++){
    grArrayNcl[isec] = new TGraphErrors((profArrayNcl[isec]));
    grArrayNclUsed[isec] = new TGraphErrors((profArrayNclUsed[isec]));
    (*pcstream)<<"sumaryTime"<<
      TString::Format("grNcl%d.=",isec).Data()<<grArrayNcl[isec]<<
      TString::Format("grNclUsed%d.=",isec).Data()<<grArrayNclUsed[isec];
  }
  for (Int_t iits=0; iits<3; iits++){
    grArrayITSNcl[iits] = new TGraphErrors((profArrayITSNcl[iits]));
    (*pcstream)<<"sumaryTime"<<
      TString::Format("grITSNcl%d.=",iits).Data()<<grArrayITSNcl[iits];
  }
  
  
  (*pcstream)<<"sumaryTime"<<"\n";
  for (Int_t isec=0; isec<nSec; isec++){
    delete 	profArrayNcl[isec];
    delete 	profArrayNclUsed[isec];
    delete 	grArrayNcl[isec];
    delete 	grArrayNclUsed[isec];
  }
  delete hisEvent;
  delete pcstream;

  printf("StatInfo.minTime\t%d\n",minTime);
  printf("StatInfo.maxTime\t%d\n",maxTime);
  delete array;
}
