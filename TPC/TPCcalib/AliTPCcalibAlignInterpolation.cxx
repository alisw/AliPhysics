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
#include "TVectorD.h"
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
#include "AliCDBManager.h"
#include "AliMagF.h"
#include "AliGRPManager.h"
#include <TGeoGlobalMagField.h>
#include "TSystem.h"
#include "TGrid.h"
#include "TCut.h"
#include "AliNDLocalRegression.h"
#include "AliMathBase.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliTPCParam.h"
#include "AliTPCClusterParam.h"
#include "AliDAQ.h"
#include "AliGRPObject.h"

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
  ULong64_t gid = esdEvent->GetHeader()->GetEventIdAsLong(); 
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
  Bool_t backupUseCorrectionMaps = transform->GetCurrentRecoParamNonConst()->GetUseCorrectionMap();
  transform->GetCurrentRecoParamNonConst()->SetUseCorrectionMap(kFALSE);
  Bool_t backupAccountDistortions = transform->GetCurrentRecoParamNonConst()->GetAccountDistortions();
  transform->GetCurrentRecoParamNonConst()->SetAccountDistortions(kFALSE);

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
    //
    // 4.) Propagate  ITS tracks outward
    // 
    Bool_t itsOK=kTRUE;
    int npUpdITS = 0;
    for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++) {
      //trackArrayITS[iPoint].SetUniqueID(0);
      AliTPCclusterMI &cluster=clusterArray[iPoint];
      if (cluster.GetVolumeId()==0) continue;
      double alpSect = ((cluster.GetDetector()%18)+0.5)*20*TMath::DegToRad();
      double xyz[3] = {cluster.GetX(),cluster.GetY(),cluster.GetZ()}; // sector tracking frame
      paramITS.Local2GlobalPosition(xyz,alpSect);	
      Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
      paramITS.Global2LocalPosition(xyz,alpha);	 // cluster frame
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
      double alpSect = ((cluster.GetDetector()%18)+0.5)*20*TMath::DegToRad();
      double xyz[3] = {cluster.GetX(),cluster.GetY(),cluster.GetZ()}; // sector tracking frame
      paramTRD.Local2GlobalPosition(xyz,alpSect);	
      Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
      paramTRD.Global2LocalPosition(xyz,alpha); // cluster frame

      if (trdOK){
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
    Bool_t saveBit = ip->TestBit(kAlignmentBugFixedBit);
    ip->SetBit(kAlignmentBugFixedBit,kTRUE);    // Flag that the alignment bug is fixed
    Int_t timeStamp= esdEvent->GetTimeStamp();
    (*fStreamer)<<"delta"<<
      "nTracks="<<nTracks<<               // number of tracks in event (pileup indicator)
      "nPrimTracks="<<nPrimTracks<<       // number of tracks pointed to primary vertes of selected event
      "timeStamp="<<timeStamp<<           // time stamp
      "itrack="<<fTrackCounter<<          // total track #
      "gid="<<gid<<                       // global ID of the event
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
    //
    ip->SetBit(kAlignmentBugFixedBit,saveBit);
  }
  transform->GetCurrentRecoParamNonConst()->SetUseComposedCorrection( backupUseComposedCorrection);
  transform->GetCurrentRecoParamNonConst()->SetUseCorrectionMap(backupUseCorrectionMaps);
  transform->GetCurrentRecoParamNonConst()->SetAccountDistortions(backupAccountDistortions);

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
  axisName[kQ2PT]="qpt";     axisTitle[kQ2PT]="q/pt (c/GeV)";                         // to fill : track.GetSigned1Pt() 
  binsTrack[kQ2PT]=7;        xminTrack[kQ2PT]=-3;        xmaxTrack[kQ2PT]=3; 
  binsTrackITS[kQ2PT]=7;     xminTrackITS[kQ2PT]=-3;     xmaxTrackITS[kQ2PT]=3; 
  //
  axisName[kSect]="sector";  axisTitle[kSect]="Sector Number";              // to fill:   9*atan2(gy,gx)/pi+ if (sector>0) sector+18
  binsTrack[kSect]=181;      xminTrack[kSect]=-0.05;           xmaxTrack[kSect]=18.05; 
  binsTrackITS[kSect]=181;   xminTrackITS[kSect]=-0.05;        xmaxTrackITS[kSect]=18.05; 
  //
  axisName[kLocX]="R";       axisTitle[kLocX]="r (cm)";                          // to fill:    gr=sqrt(gy**2+gx**2)
  binsTrack[kLocX]=53;       xminTrack[kLocX]=85.;         xmaxTrack[kLocX]=245.; 
  binsTrackITS[kLocX]=53;    xminTrackITS[kLocX]=85.;      xmaxTrackITS[kLocX]=245.; 
  //
  axisName[kZ2X]="kZ";       axisTitle[kZ2X]="z/r";                          // to fill : gz/gr 
  binsTrack[kZ2X]=20;        xminTrack[kZ2X]=-1.0;         xmaxTrack[kZ2X]=1.0;  // +-1 for ITS+TRD and ITS+TOF 
  binsTrackITS[kZ2X]=20;     xminTrackITS[kZ2X]=-1.8;      xmaxTrackITS[kZ2X]=1.8;  // +-1.8 for the ITS 
  //
  axisName[kDelt]="delta";   axisTitle[kDelt]="#Delta (cm)";                 // to fill    local(clusterY-track.y)
  binsTrack[kDelt]=100;      xminTrack[kDelt]=-dy;        xmaxTrack[kDelt]=dy; 
  binsTrackITS[kDelt]=100;   xminTrackITS[kDelt]=-dy;     xmaxTrackITS[kDelt]=dy; 
  // 
  binsTrack[kDelt]=TMath::Min(Int_t(20.+2.*dy/0.05),100); // buffer should be smaller than 1 GBy
  if (selHis==0 ||selHis<0) fHisITSDRPhi =    new THnF("deltaRPhiTPCITS","#Delta_{Y} (cm)", kNDim, binsTrackITS,xminTrackITS, xmaxTrackITS);
  if (selHis==1 ||selHis<0) fHisITSTRDDRPhi = new THnF("deltaRPhiTPCITSTRD","#Delta_{Y} (cm) TPC-(ITS+TRD)", kNDim, binsTrack,xminTrack, xmaxTrack);
  if (selHis==2 ||selHis<0) fHisITSTOFDRPhi = new THnF("deltaRPhiTPCITSTOF","#Delta_{Y} (cm) TPC-(ITS+TOF)", kNDim, binsTrack,xminTrack, xmaxTrack);
  //
  binsTrack[kDelt]=TMath::Min(Int_t(20.+2.*dz/0.05),100); // buffer should be smaller than 1 GBy
  xminTrack[kDelt]=-dz;        xmaxTrack[kDelt]=dz; 
  xminTrackITS[kDelt]=-dz;        xmaxTrackITS[kDelt]=dz; 
  if (selHis==3 ||selHis<0) fHisITSDZ = new THnF("deltaZTPCITS","#Delta_{Z} (cm)", kNDim, binsTrackITS,xminTrackITS, xmaxTrackITS);
  if (selHis==4 ||selHis<0) fHisITSTRDDZ = new THnF("deltaZTPCITSTRD","#Delta_{Z} (cm) TPC-(ITS+TRD)", kNDim, binsTrack,xminTrack, xmaxTrack);
  if (selHis==5 ||selHis<0) fHisITSTOFDZ = new THnF("deltaZTPCITSTOF","#Delta_{Z} (cm) TPC-(ITS+TOF)", kNDim, binsTrack,xminTrack, xmaxTrack);
  //
  //
  //
  THn *hisToFill[6]={GetHisITSDRPhi(), GetHisITSTRDDRPhi(), GetHisITSTOFDRPhi(), GetHisITSDZ(), GetHisITSTRDDZ(), GetHisITSTOFDZ()};
  for (Int_t ihis=0; ihis<6; ihis++){
    if (hisToFill[ihis]) for (Int_t ivar2=0;ivar2<kNDim;ivar2++){ 
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

void    AliTPCcalibAlignInterpolation::FillHistogramsFromChain(const char * residualList, Double_t dy, Double_t dz, 
							       Int_t startTime, Int_t stopTime, Int_t maxStat, 
							       Int_t selHis,const char * residualInfoFile,
							       Bool_t fixAlignmentBug){
  /**
   * Trees with point-track residuals to residual histogram
   * @param residualList  text file with tree list
   * Output - ResidualHistograms.root file with hitogram within AliTPCcalibAlignInterpolation object
   */
  //
  //
  // 
  ::Info(" AliTPCcalibAlignInterpolation::FillHistogramsFromChain","Start %s\n", residualList);
  Int_t cacheSize= 200000000;
  if (gSystem->Getenv("treeCacheSize")) cacheSize=TString(gSystem->Getenv("treeCacheSize")).Atoi();
  Bool_t autoCache = kFALSE;
  if (gSystem->Getenv("autoCacheSize")) {
    autoCache = Bool_t(TString(gSystem->Getenv("autoCacheSize")).Atoi());
  }
  Int_t cacheLearnEntries = 1;
  if (gSystem->Getenv("cacheLearnEntriesProjection")) {
    cacheLearnEntries = TString(gSystem->Getenv("cacheLearnEntriesProjection")).Atoi();
  }
  Printf("************* cacheSize = %d, autoCache = %d, cacheLearnEntries = %d", cacheSize, (Int_t)autoCache, cacheLearnEntries);

  const Int_t kNDim1 = kNDim-1;
  const Double_t kernelSigma2I[4]={1./0.25,1./0.25,1./0.25,1./0.25};  // inverse kernel sigma in bin width units
  const Double_t kFillGap=0.02  ;  // weight for the "non primary distortion info" - 
  //                                used to fill the gap without measurement (PHOS hole)
  const Double_t kFillGapITS=0.01;
  // track and cluster quality cuts - see also AliTPCcalibAlignInterpolation::CalculateDistance
  const Int_t   kMaxSkippedCluster=10;  // 10 cluster
  const Float_t kMaxRMSTrackCut=2.0;    // maximal RMS (cm) between the tracks 
  const Float_t kMaxRMSClusterCut=0.3;    // maximal RMS (cm) between the cluster and local mean
  const Float_t kMaxDeltaClusterCut=0.5;    // maximal delta(cm) between the cluster and local mean
  //
  // gap weight is kFillGap + Exp(-dist): don't calculate exponent if dist is >kMaxExpArg
  const Double_t kMaxExpArg = -TMath::Log(TMath::Max(kFillGap*0.1, 1e-3)); 
  //
  Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
  float bz = fixAlignmentBug ? InitForAlignmentBugFix(runNumber) : 0;

  // 0.) Load current information file and bookd variables
  // 
  const Int_t nSec=81;         // 72 sector +5 sumarry info+ 4 medians +
  const Double_t kMaxZSect[2]={2.49725e+02,2.49698e+02};
  TVectorF meanNcl(nSec);      // mean current estator ncl per sector
  TVectorF meanNclUsed(nSec);  // mean current estator ncl per sector
  Double_t meanTime=0, maxTime=startTime, minTime=stopTime;
  Int_t currentTrack=0;  

  TFile *finfo = TFile::Open(residualInfoFile);
  TTree *treeInfo=0;
  if (finfo) {
    treeInfo=(TTree*)finfo->Get("summaryTime"); 
  }else{
    ::Fatal("AliTPCcalibAlignInterpolation::FillHistogramsFromChain","residualInfoFile %s does not exist",residualInfoFile);
  }  
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
  }else{
    ::Fatal("AliTPCcalibAlignInterpolation::FillHistogramsFromChain","residualInfoFile %s does not contain tree summaryTime",residualInfoFile);
  }
  //
  // 0.a) Load drift velocity calibration in case availbel
  //
  TVectorD     *vdriftParam=0;
  TGraphErrors *vdriftGraph=0;  
  TFile *fdrift = TFile::Open("fitDrift.root"); 
  AliSysInfo::AddStamp("FillHistogramsFromChain.LoadDriftBegin",1,0);
  if (fdrift){
    TTree * tree = (TTree*)fdrift->Get("fitTimeStat");
    if (tree==NULL){
      ::Error("LoadDriftCalibration FAILED", "tree fitTimeStat not avaliable in file fitDrift.root");
    }else{      
      tree->SetBranchAddress("grTRDReg.",&vdriftGraph);
      tree->SetBranchAddress("paramRobust.",&vdriftParam);
      tree->GetEntry(0);
      if (vdriftGraph==NULL || vdriftGraph->GetN()<=0){
	::Info("LoadDriftCalibration FAILED", "ITS/TRD drift calibration not availalble. Trying ITS/TOF");
	tree->SetBranchAddress("grTOFReg.",&vdriftGraph);
	tree->GetEntry(0);
      }else{
	::Info("LoadDriftCalibration", "tree fitTimeStat loaded from the tree");
      }
    }
    
  }else{
    ::Error("LoadDriftCalibration FAILED", "fitDrift.root not present");
  }
  AliSysInfo::AddStamp("FillHistogramsFromChain.LoadDriftEND",1,1);
  //
  // 1.) Fill histograms and mean informations
  //
  const Int_t knPoints=kMaxRow;
  AliTPCcalibAlignInterpolation * calibInterpolation = new  AliTPCcalibAlignInterpolation("calibInterpolation","calibInterpolation",kFALSE);
  calibInterpolation->CreateResidualHistosInterpolation(dy,dz,selHis);
  TString branches[6]={"its0.","trd0.","tof0.", "its1.","trd1.","tof1."};
  //
  TVectorF *vecDeltaOther= 0;
  TVectorF *vecDeltaOtherITS= 0;
  TVectorF *vecDelta= 0;
  TVectorF *vecDeltaITS= 0;
  TVectorF *vecR=0;
  TVectorF *vecSec=0;
  TVectorF *vecPhi=0;
  TVectorF *vecZ=0;
  TVectorF *vecLocalDelta = new TVectorF(knPoints);
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
  Long64_t fillCounter=0;
  Long64_t clusterCounter=0;
  AliSysInfo::AddStamp("FillHistogramsFromChain.BEGIN",2,0);
  for (Int_t ihis=0; ihis<6; ihis++){    
    if (selHis>=0 && ihis!=selHis) continue;
    Double_t binWidth[4]={0};
    for (Int_t idim=0; idim<4; idim++) binWidth[idim]=hisToFill[ihis]->GetAxis(idim)->GetBinWidth(1);
    for (Int_t iesd=0; iesd<nesd; iesd++){
      TStopwatch timerFile;
      TString fileNameString(esdArray->At(iesd)->GetName());
      if (fileNameString.Contains("alien://") && (!gGrid || (gGrid && !gGrid->IsConnected()))) TGrid::Connect("alien://");
      TFile *esdFile = TFile::Open(fileNameString.Data(),"read");
      if (!esdFile) continue;
      TTree *tree = (TTree*)esdFile->Get("delta");
      if (!autoCache) {
	tree->SetCacheSize(cacheSize);
	tree->SetCacheLearnEntries(cacheLearnEntries);
      }
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
      tree->SetBranchStatus("track.*",kTRUE);      
      tree->SetBranchAddress("vecR.",&vecR);
      tree->SetBranchAddress("vecSec.",&vecSec);
      tree->SetBranchAddress("vecPhi.",&vecPhi);
      tree->SetBranchAddress("vecZ.",&vecZ);
      tree->SetBranchAddress("track.",&param);
      //
      // before doing anything else, check if alignment bug indeed must be fixed
      Bool_t fixAlignmentBugForFile = fixAlignmentBug;
      if (fixAlignmentBugForFile) {
	TBranch *brParam = tree->GetBranch("track.");
	brParam->GetEntry(0);
	if (param->TestBit(kAlignmentBugFixedBit)) {
	  ::Info(" AliTPCcalibAlignInterpolation::FillHistogramsFromChain",
		 "AlignmentBugFix is requested but is not needed for this chunk\n");
	  fixAlignmentBugForFile = kFALSE;
	}
      }
      //
      br->SetAddress(&timeStamp);
      if (tree->GetBranch("npValid")!=NULL) {
	tree->SetBranchStatus("npValid",kTRUE);
	tree->SetBranchAddress("npValid",&npValid);
      }
      tree->SetBranchStatus(branches[ihis],kTRUE);
      tree->SetBranchAddress(branches[ihis],&vecDelta);
      // 
      // if aligment bug fix is needed, we need also other delta
      if (fixAlignmentBugForFile) {
	int ihisOther = ihis<=2 ? ihis+3 : ihis-3;
	tree->SetBranchStatus(branches[ihisOther],kTRUE);
	tree->SetBranchAddress(branches[ihisOther],&vecDeltaOther);
      }

      if (ihis<=2 &&ihis!=0){
	tree->SetBranchStatus(branches[0],kTRUE);
	tree->SetBranchAddress(branches[0],&vecDeltaITS);
	if (fixAlignmentBugForFile) {
	  tree->SetBranchStatus(branches[3],kTRUE);
	  tree->SetBranchAddress(branches[3],&vecDeltaOtherITS);
	}
      }
      else if (ihis>2 && ihis!=3){
	tree->SetBranchStatus(branches[3],kTRUE);
	tree->SetBranchAddress(branches[3],&vecDeltaITS);
	if (fixAlignmentBugForFile) {
	  tree->SetBranchStatus(branches[0],kTRUE);
	  tree->SetBranchAddress(branches[0],&vecDeltaOtherITS);
	}
      }
      
      // prepare aux info for histo bin calculation
      Long64_t nBProd[kNDim] = {0};
      int nBinDim[kNDim];
      double bsize[kNDim],bsizeI[kNDim],limMin[kNDim],limMax[kNDim];
      nBProd[kNDim1] = 1;
      THn* curHis = hisToFill[ihis];
      TNDArrayT<float>& arrND = (TNDArrayT<float>&)curHis->GetArray(); // for direct access
      for (int i=kNDim;i--;) {
	TAxis* ax = curHis->GetAxis(i);
	limMin[i] = ax->GetXmin();
	limMax[i] = ax->GetXmax();
	nBinDim[i] = ax->GetNbins();
	bsize[i]  = (limMax[i]-limMin[i])/nBinDim[i];
	bsizeI[i] = 1./bsize[i];
	if (i<kNDim1) nBProd[i] = nBProd[i+1]*(nBinDim[i+1]+2); // +2 to account for under/over-flows
      }
      //
      Int_t ntracks=tree->GetEntries();
      if (!autoCache) {
	tree->SetCacheSize(0);
	tree->SetCacheSize(cacheSize);
	tree->SetCacheLearnEntries(cacheLearnEntries);
	tree->AddBranchToCache("*", kTRUE);
      }
      //
      for (Int_t itrack=0; itrack<ntracks; itrack++){
	if (startTime>0){
	  br->GetEntry(itrack);
	  if (timeStamp<startTime  || timeStamp>stopTime) continue;
	  hisTime->Fill(timeStamp);
	}
	tree->GetEntry(itrack);
	Double_t corrTime = (vdriftGraph!=NULL) ? vdriftGraph->Eval(timeStamp):0;
	const Float_t *vSec= vecSec->GetMatrixArray();
	const Float_t *vPhi= vecPhi->GetMatrixArray();
	const Float_t *vR  = vecR->GetMatrixArray();
	const Float_t *vZ  = vecZ->GetMatrixArray();
	const Float_t *vDelta  = vecDelta->GetMatrixArray();
	const Float_t *vDeltaITS  = (vecDeltaITS!=NULL) ? vecDeltaITS->GetMatrixArray():0;
	const Float_t *vDeltaOther = vecDeltaOther ? vecDeltaOther->GetMatrixArray():0;
	const Float_t *vDeltaOtherITS = vecDeltaOtherITS ? vecDeltaOtherITS->GetMatrixArray():0;
	const Float_t *vLocalDelta=vecLocalDelta->GetMatrixArray();
	//
	// calculate distance and aplly track and cluster quality cut
	Float_t rmsTrack=3, rmsCluster=1;
	Int_t nSkippedCluster=AliTPCcalibAlignInterpolation::CalculateDistance(*vecDelta,*vecDeltaITS, *vecSec, *vecLocalDelta, npValid, rmsTrack, rmsCluster,1.5);
	if (nSkippedCluster>kMaxSkippedCluster) continue;
	if (rmsTrack>kMaxRMSTrackCut) continue;
	if (rmsCluster>kMaxRMSClusterCut) continue;	
	//
	currentTrack++;

	if (timeStamp<minTime) minTime=timeStamp;
	if (timeStamp>maxTime) maxTime=timeStamp;
	meanTime+=timeStamp;
	if (treeInfo) for (Int_t iSec=0; iSec<nSec; iSec++){
	  meanNcl[iSec]+=nclArray[iSec]->Eval(timeStamp);
	  meanNclUsed[iSec]+=nclArrayUsed[iSec]->Eval(timeStamp);
	}

	if (maxStat>0 &&currentTrack>maxStat) break;
	double xxx[kNDim] = {0.};
	//for (Int_t ipoint=0; ipoint<knPoints; ipoint++){
	for (Int_t ipoint=0; ipoint<npValid; ipoint++){
	  if (vR[ipoint]<=0 || vDelta[ipoint]<-990.) continue;
	  if (TMath::Abs(vDelta[ipoint])<0.000001) continue; // RS Do we need this?
	  if (TMath::Abs(vLocalDelta[ipoint])> kMaxDeltaClusterCut) continue;
	  float phiUse = vPhi[ipoint], rUse = vR[ipoint], zUse = vZ[ipoint];
	  float deltaITSUse = (vDeltaITS) ? vDeltaITS[ipoint]:0;
	  float deltaRefUse = vDelta[ipoint];
	  int rocID = TMath::Nint(vSec[ipoint]);

	  if (fixAlignmentBugForFile) { // was it produced w/o bug fix
	    float dAux = vDeltaOther[ipoint];
	    if (ihis<3) FixAlignmentBug(rocID, param->GetParameter()[4], bz, phiUse, rUse, zUse, deltaRefUse, dAux);
	    else        FixAlignmentBug(rocID, param->GetParameter()[4], bz, phiUse, rUse, zUse, dAux, deltaRefUse);
	    //
	    if (vDeltaITS) {
	      dAux = vDeltaOtherITS[ipoint];
	      float phiAux = vPhi[ipoint], rAux = vR[ipoint], zAux = vZ[ipoint];
	      if (ihis<3) FixAlignmentBug(rocID, param->GetParameter()[4], bz, phiAux, rAux, zAux, deltaITSUse, dAux);
	      else        FixAlignmentBug(rocID, param->GetParameter()[4], bz, phiAux, rAux, zAux, dAux, deltaITSUse);
	    }
	    //
	  }
	  //
	  Double_t sector=9.*phiUse/TMath::Pi();
	  if (sector<0) sector+=18;
	  Double_t deltaPhi=phiUse-TMath::Pi()*(Int_t(sector)+0.5)/9.;
	  Double_t localX = TMath::Cos(deltaPhi)*rUse;

	  xxx[kQ2PT] = param->GetParameter()[4];
	  xxx[kSect] = sector;
	  xxx[kLocX] = localX;
	  Double_t side=-1.+2.*((rocID%36)<18);
	  double maxZ = kMaxZSect[side<0];
	  xxx[kZ2X] = (zUse*side)<-1 ? side*0.001 : zUse/localX; // do not mix z on A side and C side ?? RS
	  // apply drift velocity calibration if available
	  
	  if (ihis>2){  // if z residuals and vdrift calibration existing
	    Double_t drift = (side>0) ? maxZ-zUse : zUse+maxZ;
	    Double_t gy    = TMath::Sin(phiUse)*localX;
	    Double_t pvecFit[3];
	    pvecFit[0]= side;             // z shift (cm)
	    pvecFit[1]= drift*gy/maxZ;   // global y gradient
	    pvecFit[2]= drift;            // drift length
	    Double_t expected = (vdriftParam!=NULL) ? (*vdriftParam)[0]+(*vdriftParam)[1]*pvecFit[0]+(*vdriftParam)[2]*pvecFit[1]+(*vdriftParam)[3]*pvecFit[2]:0;
	    deltaRefUse= side*(deltaRefUse*side-(expected+corrTime*drift));
	    deltaITSUse= side*(deltaITSUse*side-(expected+corrTime*drift));
	  }
	  xxx[kDelt] = deltaRefUse;
	  clusterCounter++;
	  //
	  // calculate axis bins and global bin
	  Int_t binIndex[kNDim]={0};
	  //
	  if (xxx[kNDim1]<limMin[kNDim1]) binIndex[kNDim1] = 0; // underflow
	  else if (xxx[kNDim1]<limMax[kNDim1]) binIndex[kNDim1] = 1+int((xxx[kNDim1] - limMin[kNDim1])*bsizeI[kNDim1]); // range
	  else binIndex[kNDim1] = nBinDim[kNDim1]+1; // oveflow
	  //
	  ULong64_t binToFill = binIndex[kNDim1]; // global bin
	  for (int id=kNDim1;id--;) {
	    if (xxx[id]<limMin[id]) binIndex[id] = 0; // underflow
	    else if (xxx[id]<limMax[id]) binIndex[id] = 1+int((xxx[id] - limMin[id])*bsizeI[id]); // range
	    else binIndex[id] = nBinDim[id]+1; // oveflow
	    binToFill += binIndex[id]*nBProd[id];
	  }
	  //
	  if (vDeltaITS){
	    xxx[kDelt] = deltaITSUse;
	    int binDeltITS;
	    if (deltaITSUse<limMin[kDelt]) binDeltITS = 0; // underflow
	    else if (deltaITSUse<limMax[kDelt]) binDeltITS = 1+int((deltaITSUse - limMin[kDelt])*bsizeI[kDelt]); // range
	    else binDeltITS = nBinDim[kDelt]+1; // oveflow
	    Long64_t binToFillITS = binToFill + (binDeltITS-binIndex[kDelt])*nBProd[kDelt]; // global bin for ITS
	    arrND.At(binToFillITS) += kFillGapITS;	    // curHis->Fill(xxx,kFillGapITS);
	    xxx[kDelt] = deltaRefUse;
	  }
	  arrND.At(binToFill) += 1.; // curHis->Fill(xxx,1.);

	  Double_t dbinCenter[kNDim];	  
	  for (Int_t idim=kNDim;idim--;) {
	    // fractional distance to the center of the bin
	    // double bincenter = limMin[idim]+bsize[idim]*(binIndex[idim]-0.5); // binindex=0 is underflow!
	    // dbinCenter[idim] = (xxx[idim]-bincenter)*bsizeI[idim];
	    dbinCenter[idim] = (xxx[idim]-limMin[idim])*bsize[idim] - (binIndex[idim]-0.5); // binindex=0 is underflow!
	  }
	  // dbinCenter[kDelt] = 0;
	  //
	  for (Int_t iqpt=-1; iqpt<=1; iqpt++){  //qpt
	    int qptBin = binIndex[kQ2PT]+iqpt;
	    if (qptBin<0||qptBin>nBinDim[kQ2PT]) continue;
	    double dqpt = dbinCenter[kQ2PT] + iqpt;
	    dqpt *= dqpt*kernelSigma2I[kQ2PT];
	    binToFill += iqpt*nBProd[kQ2PT];
	    //for (Int_t isec=0; isec<=0; isec++){  //sector - (Not defined yet if we should make bin respone functio and unfold later) 
	    //  int secBin = binIndex[kSect]+isec;
	    //  if (secBin<0||secBin>nBinDim[kSect]) continue;
	    //  double dsec = dbinCenter[kSect] + isec;
	    //  dsec *= dsec*kernelSigma2I[kSect];
	    //  binToFill += isec*nBProd[kSect];
	    //  for (Int_t ilocx=0; ilocx<=0; ilocx++){   //local X
	    //    int locxBin = binIndex[kLocX]+ilocx;
	    //    if (locxBin<0||locxBin>nBinDim[kLocX]) continue;
	    //    double dlocx = dbinCenter[kLocX] + ilocx;
	    //    dlocx *= dlocx*kernelSigma2I[kLocX];
	    //    binToFill += ilocx*nBProd[kLocX];
	    for (Int_t iz2x=-2; iz2x<=2; iz2x++){ // Z/x
	      if ( xxx[kZ2X]*(xxx[kZ2X]+bsize[kZ2X]*iz2x) < 0 ) continue; // do not mix a and C side
	      int z2xBin = binIndex[kZ2X]+iz2x;
	      if (z2xBin<0||z2xBin>nBinDim[kZ2X]) continue; 
	      double dz2x = dbinCenter[kZ2X] + iz2x;
	      dz2x *= dz2x*kernelSigma2I[kZ2X];
	      binToFill += iz2x*nBProd[kZ2X];
	      //
	      {
		// Looping is over, fill histo
		Double_t weightAll= 2.*(dz2x+dqpt); // 2.*(dz2x+dqpt+dsec+dlocx);
		weightAll = weightAll>kMaxExpArg ? kFillGap : kFillGap+TMath::Exp(-weightAll);
		arrND.At(binToFill) += weightAll; // curHis->Fill(...)
		if (fillCounter==0) printf("Start to Fill\n");
		fillCounter++;
	      }
	      //
	      binToFill -= iz2x*nBProd[kZ2X];
	    } // z2x fill loop
	    //   binToFill -= ilocx*nBProd[kLocX];
	    // }   // xloc fill loop	      
	    //    binToFill -= isec*nBProd[kSect];
	    //  }     // sector fill loop 
	    binToFill -= iqpt*nBProd[kQ2PT]; // restore
	  } // qpt fill loop
	}
      }
      tree->PrintCacheStats();

      timerFile.Print();
      delete tree;
      delete esdFile;
      
    }    
    fout->GetFile()->cd();
    hisToFill[ihis]->Write();
  }
  AliSysInfo::AddStamp("FillHistogramsFromChain.END",2,1);
  if (hisTime) hisTime->Write();
  ::Info(" AliTPCcalibAlignInterpolation::FillHistogramsFromChain","End of processing\n");
  timerAll.Print();
  printf("StatInfo.fillCounter:\t%lld\n",fillCounter);
  printf("StatInfo.clusterCounter:\t%lld\n",clusterCounter);
  printf("StatInfo.trackCounter:\t%d\n",currentTrack);
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
    "runNumber="<<runNumber<<        // runNumber
    "selHis="<<selHis<<              // selected histogram type
    "fillCounter="<<fillCounter<<    // number of histogram fills
    "clusterCounter="<<clusterCounter<<    // number of clusters used for fill
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
    TString fileNameString(esdArray->At(iesd)->GetName());
    if (fileNameString.Contains("alien://") && (!gGrid || (gGrid && !gGrid->IsConnected()))) TGrid::Connect("alien://");
    TFile *esdFile = TFile::Open(fileNameString.Data(),"read");
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
      double alpSect = ((cl->GetDetector()%18)+0.5)*20*TMath::DegToRad();
      double xyz[3] = {cl->GetX(),cl->GetY(),cl->GetZ()}; // sector tracking frame
      param->Local2GlobalPosition(xyz,alpSect);	
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
    ::Error("AliTPCcalibAlignInterpolation::AddFriendDistortionTree", "file %s not readable", fname);
    return 0;
  }
  TTree * treeFriend = (TTree*) fin->Get(treeName);
  
  if (treeFriend==NULL){
    ::Error("AliTPCcalibAlignInterpolation::AddFriendDistortionTree", "file %s not readable", fname);
    return 0;
  }
  if (tree==NULL) {
    tree = treeFriend;
    tree->SetName("Default");
    tree->SetTitle("Default");
    treeFriend = (TTree*) fin->Get(treeName);
  }
  {
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
    ::Error("GetResidualStatInfo", "Wrong input list %s", inputList);
    return;
  }
  TStopwatch timer;
  //
  // 1.) Cache information per files - beginTime and endTime for file
  //
  TStopwatch timer1;
  TTreeSRedirector * pcstream = new TTreeSRedirector("residualInfo.root", "recreate");
  Int_t cacheSize=100000000; // 100 MBy cache
  if (gSystem->Getenv("treeCacheSize")) cacheSize=TString(gSystem->Getenv("treeCacheSize")).Atoi();
  Bool_t autoCache = kFALSE;
  if (gSystem->Getenv("autoCacheSize")) autoCache = Bool_t(TString(gSystem->Getenv("autoCacheSize")).Atoi());
  Int_t cacheLearnEntries = 1;
  if (gSystem->Getenv("cacheLearnEntriesStatInfo")) {
    cacheLearnEntries = TString(gSystem->Getenv("cacheLearnEntriesStatInfo")).Atoi();
  }
  Printf("************* cacheSize = %d, autoCache = %d, cacheLearnEntries = %d", cacheSize, (Int_t)autoCache, cacheLearnEntries);

  TChain * chainInfo  = AliXRDPROOFtoolkit::MakeChain("residual.list","eventInfo",0,-1);
  if (!autoCache) {
    chainInfo->SetCacheSize(cacheSize);
    chainInfo->SetCacheLearnEntries(cacheLearnEntries);
  }
  TChain * chainTracks  = AliXRDPROOFtoolkit::MakeChain("residual.list","delta",0,-1);
  if (!autoCache) {
    chainTracks->SetCacheSize(cacheSize);
    chainTracks->SetCacheLearnEntries(cacheLearnEntries);
  }
  //
  Int_t gidRounding=128;                        // git has to be rounded
  Int_t neventsAll=chainInfo->GetEntries();     // total amount of events
  Int_t ntracksAll=chainTracks->GetEntries();   // total amount of tracks

  //chainInfo->SetEstimate(-1);                 // using -1  make a problem - too much memory allocate with autimatic switch . crash in case of 4MBy limits in the next Draw 
  chainInfo->SetEstimate(neventsAll);   
  chainInfo->Draw("timeStamp:gid/128","timeStamp>0","goff");          
  //
  //
  Long64_t minTime=0,maxTime=0;
  double minGID=0,maxGID=0,meanGID=0,meanTime=0;
  if (neventsAll) {
    double minTimeD=0,maxTimeD=0;
    TStatToolkit::GetMinMaxMean(chainInfo->GetV1(),neventsAll,minTimeD,maxTimeD, meanTime);
    minTime = minTimeD;
    maxTime = maxTimeD;
    TStatToolkit::GetMinMaxMean(chainInfo->GetV2(),neventsAll,minGID,maxGID, meanGID);
    minGID*=128;
    maxGID*=128;
    meanGID*=128;
  }
  (*pcstream)<<"summary1"<<
    "id="<<id<<                // chain id - usually should be run number
    "nevents="<<neventsAll<<   // total number of events
    "ntracks="<<ntracksAll<<   // total number of tracks
    "minTime="<<minTime<<      // minimal time stamp in sample
    "maxTime="<<maxTime<<      // max time stamp in sample
    "meanTime="<<meanTime<<    // mean time
    "minGID="<<minGID<<        // minimal event gid in sample (rounded to 128)
    "maxGID="<<maxGID<<        // max  event gid in sample (rounded to 128)
    "meanGID="<<meanGID<<      // mean event gid in sample (rounded to 128)
    "\n";
  delete pcstream;
  ::Info("GetResidualStatInfo","Total time");
  timer1.Print();
  //
  // 2.) Cache information per time interval
  //
  TStopwatch timer2;
  pcstream = new TTreeSRedirector("residualInfo.root", "update");
  //  Int_t entries = neventsAll;

  Long64_t minTimeQA = timeInterval*(minTime/timeInterval);
  Long64_t maxTimeQA = timeInterval*(1+(maxTime/timeInterval));
  Int_t nIntervals=(maxTimeQA-minTimeQA)/timeInterval;
  Int_t nIntervalsQA=(maxTimeQA-minTimeQA)/15;
  //
  TH1F  * hisEvent= new TH1F("hisEvent","hisEvent",nIntervalsQA,minTimeQA,maxTimeQA);
  const Int_t nSec=81; // 72 sector +5 sumarry info+ 4 medians
  TProfile * profArrayNcl[nSec]={0};
  TProfile * profArrayNclUsed[nSec]={0};
  TGraphErrors * grArrayNcl[nSec]={0};
  TGraphErrors * grArrayNclUsed[nSec]={0};
  TProfile * profArrayITSNcl[3]={0};
  TGraphErrors * grArrayITSNcl[3]={0};
  
  for (Int_t isec=0; isec<nSec; isec++){
    profArrayNcl[isec]=new TProfile(TString::Format("TPCnclSec%d",isec).Data(), TString::Format("TPCnclSec%d",isec).Data(), nIntervalsQA,minTimeQA,maxTimeQA);
    profArrayNclUsed[isec]=new TProfile(TString::Format("TPCnclUsedSec%d",isec).Data(), TString::Format("TPCnclUsedSec%d",isec).Data(), nIntervalsQA,minTimeQA,maxTimeQA);
  }
   for (Int_t iits=0; iits<3; iits++){
    profArrayITSNcl[iits]=new TProfile(TString::Format("ITSnclSec%d",iits).Data(), TString::Format("ITSnclSec%d",iits).Data(), nIntervalsQA,minTimeQA,maxTimeQA);    
  }

  TVectorF *vecNClTPC=0;
  TVectorF *vecNClTPCused=0;
  Int_t nITS[3]={0};
  Int_t timeStamp=0;
  for (Int_t iFile=0; iFile<nFiles; iFile+=skip){
    timer.Start();
    TString fileName = array->At(iFile)->GetName();
    if (fileName.Contains("alien://") && (!gGrid || (gGrid && !gGrid->IsConnected()))) TGrid::Connect("alien://");
    printf("%d\t%s\n",iFile,fileName.Data());    
    AliSysInfo::AddStamp(fileName.Data(),1,iFile);
    TFile * f = TFile::Open(fileName.Data());
    if (f==NULL) continue;
    TTree * treeInfo = (TTree*)f->Get("eventInfo"); 
    if (treeInfo==NULL) continue;
    if (!autoCache) {
      treeInfo->SetCacheSize(cacheSize);
      treeInfo->SetCacheLearnEntries(cacheLearnEntries);
    }
    treeInfo->SetBranchAddress("vecNClTPC.",&vecNClTPC);
    treeInfo->SetBranchAddress("vecNClTPCused.",&vecNClTPCused);
    treeInfo->SetBranchAddress("nSPD",&nITS[0]);
    treeInfo->SetBranchAddress("nSDD",&nITS[1]);
    treeInfo->SetBranchAddress("nSSD",&nITS[2]);
    treeInfo->AddBranchToCache(treeInfo->GetBranch("vecNClTPC."), kTRUE);
    treeInfo->AddBranchToCache(treeInfo->GetBranch("vecNClTPCused."), kTRUE);
    treeInfo->AddBranchToCache(treeInfo->GetBranch("nSPD"), kTRUE);
    treeInfo->AddBranchToCache(treeInfo->GetBranch("nSDD"), kTRUE);
    treeInfo->AddBranchToCache(treeInfo->GetBranch("nSSD"), kTRUE);
    Bool_t hasTimeStamp=(treeInfo->GetBranch("timeStamp")!=NULL);
    if (hasTimeStamp) treeInfo->SetBranchAddress("timeStamp",&timeStamp);
    if (!hasTimeStamp) ((TBranch*)(treeInfo->GetListOfBranches()->At(1)))->SetAddress(&timeStamp);
    treeInfo->AddBranchToCache(treeInfo->GetBranch("timeStamp"), kTRUE);
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
    treeInfo->PrintCacheStats();
    AliSysInfo::AddStamp(fileName.Data(),2,iFile);
    delete treeInfo;
    delete f;
    AliSysInfo::AddStamp(fileName.Data(),3,iFile);
    
  }
  timer2.Print();
  TGraphErrors grEvent(hisEvent);
  (*pcstream)<<"summaryTime"<<
    "id="<<id<<
    "grEvent.="<<&grEvent;
  for (Int_t isec=0; isec<nSec; isec++){
    grArrayNcl[isec] = new TGraphErrors((profArrayNcl[isec]));
    grArrayNclUsed[isec] = new TGraphErrors((profArrayNclUsed[isec]));
    (*pcstream)<<"summaryTime"<<
      TString::Format("grNcl%d.=",isec).Data()<<grArrayNcl[isec]<<
      TString::Format("grNclUsed%d.=",isec).Data()<<grArrayNclUsed[isec];
  }
  for (Int_t iits=0; iits<3; iits++){
    grArrayITSNcl[iits] = new TGraphErrors((profArrayITSNcl[iits]));
    (*pcstream)<<"summaryTime"<<
      TString::Format("grITSNcl%d.=",iits).Data()<<grArrayITSNcl[iits];
  }
  
  
  (*pcstream)<<"summaryTime"<<"\n";
  for (Int_t isec=0; isec<nSec; isec++){
    delete 	profArrayNcl[isec];
    delete 	profArrayNclUsed[isec];
    delete 	grArrayNcl[isec];
    delete 	grArrayNclUsed[isec];
  }
  delete hisEvent;
  delete pcstream;

  printf("StatInfo.minTime\t%lld\n",minTime); //this formatting does not work on my Debian why it was changed ?
  printf("StatInfo.maxTime\t%lld\n",maxTime);
  //printf("StatInfo.minTime\t%f\n",Double_t(minTime)); //this formatting does not work on my Debian why it was changed ?
  //printf("StatInfo.maxTime\t%f\n",Double_t(maxTime));

  
  delete array;
}



Bool_t AliTPCcalibAlignInterpolation::FitDrift(double deltaT, double sigmaT, double time0, double_t time1,Bool_t fixAlignmentBug, Bool_t tofBCValidation){
  //
  //  Fit time dependence of the drift velocity for ITS-TRD and ITS-TOF scenario
  /*  Intput:
        "residual.list"   - ascii file with files assumed to be in working directory
        double deltaT     - time binning for drift velocity
        double sigmaT     - kernel width for time smoothing
        double time0      - starting time for drift velocity calculation
        double time1      - stop time for drift velocty calculation
        * in case time0 and time1 not specified - full time range in the selected sample used (time0=minTime, time1=maxTime)
      Output:
        "fitDrift.root"   - small file with the drift velocity calibration created  
         robustFit tree   - drift vlocity, time0, z shift and gy gradiend calibration using TLinearFitter::EvalRobust estimator
                          -   20 statistically independent values to QA procedure
                          -   resulting values choosen using robust median estimator
         fitTime          - set of graphs drift velocity as function of time
                          - n different graph values 
                          - resulting time dependent graph used as an logal regression with kernel sigmaT  per interval 
         fitTimeStat      - comparing median and local regression statistic 
  */  
  /*
    To do:
       - 1.) make version:
       - 1.a) with- without bunch 0 crossing 
       - 1.b) disentagle fit for the TRD and for the TOF
       - 2.) use time bin sigma estimate for the outlier rejection in the AliNDLocalFit fit 

   */
 
  /*
    fileList="residual.list"
    time0=0; time1=0;
    deltaT=60; sigmaT=600
    fitDrift(60,00,0);
  */  
  const Double_t kMinEntries=1000;
  const Double_t kInvalidR = 70;
  const Double_t kInvalidRes = -900;
  const Double_t kMaxDist0=20;
  const Double_t kMaxDist1=5;
  const Double_t kDumpSample=0.01;
  const Double_t kBCcutMin=-5;
  const Double_t kBCcutMax=20;
  const Double_t robFraction=0.95;
  //const Double_t regRobust=0.0000000001;
  
  Int_t maxEntries=1000000;
  Int_t maxPointsRobust=4000000;
  //
  const Double_t kMaxZSect[2]={2.49725e+02,2.49698e+02};
  TCut  selection="";
  Int_t entriesAll=0;
  Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();

  Int_t cacheSize=100000000; // 100 MBy cache
  if (gSystem->Getenv("treeCacheSize")) cacheSize=TString(gSystem->Getenv("treeCacheSize")).Atoi();
  Bool_t autoCache = kFALSE;
  if (gSystem->Getenv("autoCacheSize")) autoCache = Bool_t(TString(gSystem->Getenv("autoCacheSize")).Atoi());
  Int_t cacheLearnEntries = 1;
  if (gSystem->Getenv("cacheLearnEntriesVDrift")) {
    cacheLearnEntries = TString(gSystem->Getenv("cacheLearnEntriesVDrift")).Atoi();
  }
  Printf("************* cacheSize = %d, autoCache = %d, cacheLearnEntries = %d", cacheSize, (Int_t)autoCache, cacheLearnEntries);

  //
  //
  if (deltaT<=0 || sigmaT<=0){
    ::Error("AliTPCcalibAlignInterpolation::FitDrift FAILED ","Invalid parameter value for the deltaT %.1f and sigmaT %.1f", deltaT, sigmaT);
    return kFALSE;
  }
  if (TString(gSystem->GetFromPipe("cat residual.list | grep -c alien://")).Atoi()>0) TGrid::Connect("alien");

  TChain * chainDelta = AliXRDPROOFtoolkit::MakeChain("residual.list","delta",0,-1);
  //
  // before doing anything else, check if alignment bug indeed must be fixed
  Bool_t fixAlignmentBugForFile = fixAlignmentBug;
  AliExternalTrackParam* param = 0;
  if (fixAlignmentBugForFile) {
    TBranch *brParam = chainDelta->GetBranch("track.");
    brParam->SetAddress(&param);
    brParam->GetEntry(0);
    if (param->TestBit(kAlignmentBugFixedBit)) {
      ::Info(" AliTPCcalibAlignInterpolation::FitDrift",
	     "AlignmentBugFix is requested but is not needed for this chunk\n");
      fixAlignmentBugForFile = kFALSE;
    }
  }
  //
  if (!autoCache) {
    chainDelta->SetCacheSize(0);
    chainDelta->SetCacheSize(cacheSize);
    chainDelta->SetCacheLearnEntries(cacheLearnEntries);
  }
  AliSysInfo::AddStamp("FitDrift.chainDeltaGetEntriesBegin",1,0,0);
  entriesAll = chainDelta->GetEntries();
  AliSysInfo::AddStamp("FitDrift.chainDeltaGetEntriesEnd",1,1,0);

  if (entriesAll<kMinEntries) {
    ::Error("fitDrift FAILED","Not enough tracks in the chain.  Ntracks=%d",entriesAll); 
    return kFALSE;
  }
  maxEntries=TMath::Min(maxEntries, entriesAll);
  TTreeSRedirector *pcstream = new TTreeSRedirector("fitDrift.root","recreate");
  if (time0==time1){ 
    AliSysInfo::AddStamp("FitDrift.chainInfoGetTimeBegin",2,0,0);
    TChain * chainInfo=  AliXRDPROOFtoolkit::MakeChain("residual.list","eventInfo",0,-1); 
    if (!autoCache) {
      chainInfo->SetCacheSize(cacheSize);
      chainInfo->SetCacheLearnEntries(cacheLearnEntries);
    }
    //    chainInfo->SetEstimate(-1); // i cashed here -1 does not work
    chainInfo->SetEstimate(maxEntries); // i cashed here -1 does not work
    Int_t entries = chainInfo->Draw("timeStamp","","goff",maxEntries);
    chainInfo->PrintCacheStats();
    AliSysInfo::AddStamp("FitDrift.chainInfoGetTimeEND",2,1,0);
    if (entries) TStatToolkit::GetMinMax(chainInfo->GetV1(),entries,time0,time1);
  }
  // 0.) Cache variables:  to be done using loop
  //     Variables to cache:
  //       tof1.fElements
  //       trd1.fElements
  //       vecZ.fElements
  //       vecR.fElements
  //       vecSec.fElements
  //       vecPhi.fElements
  //       timeStamp
  //       npValid
  //
  AliSysInfo::AddStamp("FitDrift.StartCache",3,0,0);

  float bz = fixAlignmentBugForFile ? InitForAlignmentBugFix(runNumber) : 0;

  // check if TOF was present
  const float tofBCMin = -25.f;
  const float tofBCMax = 50.f;
  //
  if (tofBCValidation) {
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man->IsDefaultStorageSet()) man->SetDefaultStorage("raw://");
    if (man->GetRun()<0) man->SetRun(runNumber);  
    AliGRPObject* grp = (AliGRPObject*)man->Get(AliCDBPath("GRP/GRP/Data"))->GetObject();
    Int_t activeDetectors = grp->GetDetectorMask();
    if (!(activeDetectors&AliDAQ::DetectorPattern("TOF"))) {
      AliWarningClass("Disabling TOF BC validation since TOF is not in the GRP");
      tofBCValidation = kFALSE;
    }
    else {
      AliInfoClassF("TOF BC validation is enabled with window %.2f : %.2f ns",tofBCMin,tofBCMax);
    }
  }
  //
  chainDelta->SetEstimate(maxEntries*160/5.);
  TString toRead = "tof1.fElements:trd1.fElements:vecZ.fElements:vecR.fElements:vecSec.fElements:vecPhi.fElements:timeStamp:tofBC:its1.fElements";
  TString cutEv = "Entry$%5==0 && vecR.fElements>10";
  if (tofBCValidation) cutEv += Form(" && tofOK && !(tofBC<%.3f || tofBC>%.3f)",tofBCMin,tofBCMax);
  //
  if (fixAlignmentBugForFile) toRead += ":tof0.fElements:trd0.fElements:track.fP[4]"; // needed only in case we need to fix alignment bug
  //
  Int_t entriesFit0 = chainDelta->Draw(toRead.Data(),cutEv.Data(),"goffpara",maxEntries); 
  chainDelta->PrintCacheStats();
  AliInfoClassF("Selected %d points from first %d entries",entriesFit0,maxEntries);
  
  float lossFactorOK = 0.5f; // allowed loss factor
  float lossFactor = entriesFit0*5./160./maxEntries;
  int prescaling = 10;
  if (lossFactor < lossFactorOK) {
    prescaling = TMath::Max(lossFactor/lossFactorOK,2.0f);
  }
  Int_t entriesFit=entriesFit0/prescaling;
  AliInfoClassF("Selection loss factor: %f -> random prescaling: %d -> entriesFit: %d",lossFactor,prescaling,entriesFit);

  AliSysInfo::AddStamp("FitDrift.EndCache",3,1,0);

  AliSysInfo::AddStamp("FitDrift.BeginFill",4,1,0);

  TVectorD * deltaZTOF  = new TVectorD(entriesFit); //
  TVectorD * deltaZTRD  = new TVectorD(entriesFit); //
  TVectorD * deltaYTOF  = new TVectorD(entriesFit); //
  TVectorD * deltaYTRD  = new TVectorD(entriesFit); //
  TVectorD * vecZ      = new TVectorD(entriesFit); //
  TVectorD * vecR      = new TVectorD(entriesFit); //
  TVectorD * vecSec    = new TVectorD(entriesFit); //
  TVectorD * vecPhi    = new TVectorD(entriesFit); //
  TVectorD * vecTime   = new TVectorD(entriesFit); //
  TVectorD * vecTOFBC   = new TVectorD(entriesFit); //
  TVectorD * deltaZITS  = new TVectorD(entriesFit); //
  
  for (Int_t i=0; i<entriesFit; i++){
    Int_t index=gRandom->Rndm()*entriesFit0;
    (*deltaZTOF)[i]= chainDelta->GetVal(0)[index];
    (*deltaZTRD)[i]= chainDelta->GetVal(1)[index];
    (*vecZ)[i]= chainDelta->GetVal(2)[index];
    (*vecR)[i]= chainDelta->GetVal(3)[index];
    (*vecSec)[i]= chainDelta->GetVal(4)[index];
    (*vecPhi)[i]= chainDelta->GetVal(5)[index];
    (*vecTime)[i]= chainDelta->GetVal(6)[index];    
    (*vecTOFBC)[i]= chainDelta->GetVal(7)[index];    
    (*deltaZITS)[i]= chainDelta->GetVal(8)[index];
    //
    if ( (*vecR)[i] < kInvalidR ) continue; // there was no cluster at this point
    // if needed, fix alignment bug 
    if (fixAlignmentBugForFile) {
      float dyTOF = chainDelta->GetVal(9)[index];
      float dyTRD = chainDelta->GetVal(10)[index];
      float q2pt = chainDelta->GetVal(11)[index];
      //
      float dzITS = (*deltaZITS)[i];
      float dzTOF = (*deltaZTOF)[i];
      float dzTRD = (*deltaZTRD)[i];
      float rTOF  = (*vecR)[i], rTRD = rTOF;
      float zTOF  = (*vecZ)[i] + dzTOF-dzITS , zTRD = (*vecZ)[i] + dzTRD-dzITS;
      float phiTOF= (*vecPhi)[i], phiTRD = phiTOF;
      
      int rocID = TMath::Nint((*vecSec)[i]);

      if (dyTOF>kInvalidRes) {
	FixAlignmentBug(rocID, q2pt, bz, phiTOF, rTOF, zTOF, dyTOF, dzTOF);
	(*deltaZTOF)[i] = dzTOF;
	(*vecR)[i]      = rTOF;
	(*vecZ)[i]      = zTOF;
	(*vecPhi)[i]    = phiTOF;
      }
      if (dyTRD>kInvalidRes) {
	FixAlignmentBug(rocID, q2pt, bz, phiTRD, rTRD, zTRD, dyTRD, dzTRD);
	//
	(*deltaZTRD)[i] = dzTRD;
	(*vecR)[i]      = rTRD;
	(*vecZ)[i]      = zTRD;
	(*vecPhi)[i]    = phiTRD;
      }
      //
    }    
  }
  AliSysInfo::AddStamp("FitDrift.EndFill",4,1,0);


  //
  // 1.) Make robust first estimate
  //
  TVectorD paramRobust(5);  
  TVectorD paramTRD(5);  
  TVectorD paramTOF(5);  
  TVectorD paramRobustBC(5);  
  TVectorD paramTRDBC(5);  
  TVectorD paramTOFBC(5);
  //  
  AliSysInfo::AddStamp("FitDrift.StartRobust",2,0,0);
  TF1 * fpol1 = new TF1("f1","[0]+[1]*x",0,250);
  //
  if (pcstream->GetFile()->Get("robustFit")==0){
    for (Int_t iter=0; iter<20; iter++){
      AliSysInfo::AddStamp("FitDrift.RobustIter",3,iter,0);
      TLinearFitter *fitterRobust= new TLinearFitter(4,TString::Format("hyp%d",3).Data());
      TLinearFitter *fitterTRD= new TLinearFitter(4,TString::Format("hyp%d",3).Data());
      TLinearFitter *fitterTOF= new TLinearFitter(4,TString::Format("hyp%d",3).Data());
      TLinearFitter *fitterRobustBC= new TLinearFitter(4,TString::Format("hyp%d",3).Data());
      TLinearFitter *fitterTRDBC= new TLinearFitter(4,TString::Format("hyp%d",3).Data());
      TLinearFitter *fitterTOFBC= new TLinearFitter(4,TString::Format("hyp%d",3).Data());
      TH2F * hisDeltaZ[12];
      for (Int_t ihis=0; ihis<12; ihis++) {
	hisDeltaZ[ihis] = new TH2F("hisZ","hisZ",21,40,250,100,-5,5);
      }

      Int_t maxPoints= TMath::Min(maxPointsRobust,entriesFit); 
      for (Int_t ipoint=0; ipoint<maxPoints; ipoint++){
	Double_t  radius= (*vecR)[ipoint];
	
	if (radius<kInvalidR) continue;   // bad point

	Double_t pvecFit[10]={0};
	if (gRandom->Rndm()>0.1) continue;
	Int_t sector   = TMath::Nint((*vecSec)[ipoint]);
	Double_t side  = -1.+((sector%36)<18)*2.;
	Double_t z= (*vecZ)[ipoint];
	double maxZ = kMaxZSect[side<0];
	Double_t drift = (side>0) ? maxZ-(*vecZ)[ipoint] : (*vecZ)[ipoint]+maxZ;
	if (drift>maxZ) continue;
	Double_t phi   = (*vecPhi)[ipoint];
	Double_t gy    = TMath::Sin(phi)*radius;
	Float_t tofBC=(*vecTOFBC)[ipoint];
	pvecFit[0]= side;             // z shift (cm)
	pvecFit[1]= drift*gy/maxZ;   // global y gradient
	pvecFit[2]= drift;            // drift length
	Double_t dZTOF=(*deltaZTOF)[ipoint];
	Double_t dZTRD=(*deltaZTRD)[ipoint];
	Int_t hisIndex=(((1-side)/2)*6);

	if (iter==0){
	  if ( dZTOF>kInvalidRes &&TMath::Abs(dZTOF)<kMaxDist0) fitterRobust->AddPoint(pvecFit,dZTOF*side,1);
	  if ( dZTRD>kInvalidRes &&TMath::Abs(dZTRD)<kMaxDist0) fitterRobust->AddPoint(pvecFit,dZTRD*side,1);
	}else{
	  Double_t expected = paramRobust[0]+paramRobust[1]*pvecFit[0]+paramRobust[2]*pvecFit[1]+paramRobust[3]*pvecFit[2];
	  if ( dZTRD>kInvalidRes &&TMath::Abs(dZTRD*side-expected)<kMaxDist1) {
	    fitterRobust->AddPoint(pvecFit,dZTRD*side,1);	
	    fitterTRD->AddPoint(pvecFit,dZTRD*side,1);	
	    hisDeltaZ[hisIndex+0]->Fill(drift, dZTRD*side-expected);
	    hisDeltaZ[hisIndex+1]->Fill(drift, dZTRD*side-expected);
	    if (tofBC>kBCcutMin&&tofBC<kBCcutMax){
	      fitterRobustBC->AddPoint(pvecFit,dZTRD*side,1);	
	      fitterTRDBC->AddPoint(pvecFit,dZTRD*side,1);	
	      hisDeltaZ[hisIndex+3]->Fill(drift, dZTRD*side-expected);
	      hisDeltaZ[hisIndex+4]->Fill(drift, dZTRD*side-expected);
	    }
	  }
	  if ( dZTOF>kInvalidRes &&TMath::Abs(dZTOF*side-expected)<kMaxDist1) {
	    hisDeltaZ[hisIndex+0]->Fill(drift, dZTOF*side-expected);
	    hisDeltaZ[hisIndex+2]->Fill(drift, dZTOF*side-expected);
	    fitterRobust->AddPoint(pvecFit,dZTOF*side,1);
	    fitterTOF->AddPoint(pvecFit,dZTOF*side,1);
	    if (tofBC>kBCcutMin&&tofBC<kBCcutMax){
	      fitterRobustBC->AddPoint(pvecFit,dZTOF*side,1);
	      fitterTOFBC->AddPoint(pvecFit,dZTOF*side,1);
	      hisDeltaZ[hisIndex+3]->Fill(drift, dZTOF*side-expected);
	      hisDeltaZ[hisIndex+5]->Fill(drift, dZTOF*side-expected);
	    }
	  }

	  if (gRandom->Rndm()<kDumpSample){
	    Double_t cTime = (*vecTime)[ipoint];
	    (*pcstream)<<"dumpSample"<<
	      "iter="<<iter<<
	      "run="<<runNumber<<
	      "ctime="<<cTime<<
	      "sector="<<sector<<
	      "side="<<side<<
	      "drift="<<drift<<
	      "radius="<<radius<<
	      "phi="<<phi<<
	      "tofBC="<<tofBC<<
	      "gy="<<gy<<
	      "z="<<z<<
	      "expected="<<expected<<
	      "paramRobust.="<<&paramRobust<<      //  drift fit using all tracks
	      "dZTOF="<<dZTOF<<
	      "dZTRD="<<dZTRD<<
	      "\n";
	  }
	}
      }
      printf("iter=%d\n",iter);
      if (fitterRobust->GetNpoints()<kMinEntries*0.1){ 	
	::Error("fitDrift FAILED","Not enough points in the chain. Iter=%d\tN=%d" , iter, fitterRobust->GetNpoints());
	delete pcstream;
	return kFALSE;
      }
      //fitterRobust->EvalRobust(robFraction, regRobust);   // ROOT has to be patched to use regularization
      fitterRobust->EvalRobust(robFraction);
      //fitterRobust->Eval(); // EvalRobust sometimes failed - looks like related to the random selection of subset of data - can be only one side 
      fitterRobust->GetParameters(paramRobust);
      Double_t npoints= fitterRobust->GetNpoints();
      Double_t chi2= TMath::Sqrt(fitterRobust->GetChisquare()/npoints);
      paramRobust.Print();

      Int_t nTRD = fitterTRD->GetNpoints();
      Int_t nTOF = fitterTOF->GetNpoints();
      Int_t nRobustBC = fitterRobustBC->GetNpoints();
      Int_t nTRDBC = fitterTRDBC->GetNpoints();
      Int_t nTOFBC = fitterTOFBC->GetNpoints();

      if (nRobustBC>kMinEntries*0.1){	fitterRobustBC->EvalRobust(robFraction);  fitterRobustBC->GetParameters(paramRobustBC);}
      if (nTRD>kMinEntries*0.1){	fitterTRD->EvalRobust(robFraction);	   fitterTRD->GetParameters(paramTRD);}
      if (nTOF>kMinEntries*0.1){	fitterTOF->EvalRobust(robFraction);	   fitterTOF->GetParameters(paramTOF);}
      if (nTRDBC>kMinEntries*0.1){	fitterTRDBC->EvalRobust(robFraction);	   fitterTRDBC->GetParameters(paramTRDBC);}
      if (nTOFBC>kMinEntries*0.1){	fitterTOFBC->EvalRobust(robFraction);	   fitterTOFBC->GetParameters(paramTOFBC); }

   //    if (nRobustBC>kMinEntries*0.1){	fitterRobustBC->EvalRobust(robFraction, regRobust);  fitterRobustBC->GetParameters(paramRobustBC);}
//       if (nTRD>kMinEntries*0.1){	fitterTRD->EvalRobust(robFraction, regRobust);	   fitterTRD->GetParameters(paramTRD);}
//       if (nTOF>kMinEntries*0.1){	fitterTOF->EvalRobust(robFraction, regRobust);	   fitterTOF->GetParameters(paramTOF);}
//       if (nTRDBC>kMinEntries*0.1){	fitterTRDBC->EvalRobust(robFraction, regRobust);	   fitterTRDBC->GetParameters(paramTRDBC);}
//       if (nTOFBC>kMinEntries*0.1){	fitterTOFBC->EvalRobust(robFraction, regRobust);	   fitterTOFBC->GetParameters(paramTOFBC); }
      //
      TGraphErrors * grDeltaZ[12] = {0};
      TGraphErrors * grRMSZ[12] = {0};
      TVectorD * fitDeltaZ[12] = {0};
      TObjArray fitArray(3);
      for (Int_t ihis=0; ihis<12; ihis++){
	hisDeltaZ[ihis]->FitSlicesY(0,0,-1,0,"QNR",&fitArray);
	grDeltaZ[ihis] = new TGraphErrors(((TH1D*)fitArray.At(1)));
	grRMSZ[ihis]   = new TGraphErrors(((TH1D*)fitArray.At(2)));
	fitDeltaZ[ihis] = new TVectorD(2);
	grDeltaZ[ihis]->Fit(fpol1,"q","q");
	fpol1->GetParameters(fitDeltaZ[ihis]->GetMatrixArray());
	fitArray.Delete();
	delete hisDeltaZ[ihis];
	hisDeltaZ[ihis] = 0;
      }

      delete fitterRobust;    
      delete fitterTRD;    
      delete fitterTOF;    
      delete fitterRobustBC;    
      delete fitterTRDBC;    
      delete fitterTOFBC;    

      (*pcstream)<<"robustFit"<<
	"iter="<<iter<<
	"time0="<<time0<<
	"time1="<<time1<<
	"run="<<runNumber<<
	"npoints="<<npoints<<
	"chi2="<<chi2<<                   
	"nTRD="<<nTRD<<                      // number of TRD points
	"nTOF="<<nTOF<<                      // number of TRD points
	"nTRDBC="<<nTRDBC<<                  // number of TRD points BC
	"nTOFBC="<<nTOFBC<<                  // number of TRD points BC
	//
	"paramRobust.="<<&paramRobust<<      //  drift fit using all tracks
	"paramTRD.="<<&paramTRD<<            //  drift fit using TRD tracks
	"paramTOF.="<<&paramTOF<<	     //  drift fit using TOF tracks
	"paramRobustBC.="<<&paramRobustBC<<  //  drift fit using all tracks with BC
	"paramTRDBC.="<<&paramTRDBC<<        //  drift fit using TRD tracks with BC
	"paramTOFBC.="<<&paramTOFBC;         //  drift fit using TOF tracks with BC
      
      for (Int_t ihis=0; ihis<12; ihis++){
	(*pcstream)<<"robustFit"<<
	  TString::Format("grDeltaZ%d.=",ihis).Data()<<grDeltaZ[ihis]<<     // residual histogram drift fit - mean 
	  TString::Format("grRMSZ%d.=",ihis).Data()<<grDeltaZ[ihis]<<       //  residual histogram drift fit - rms
	  TString::Format("fitDeltaZ%d.=",ihis).Data()<<fitDeltaZ[ihis];     //  residual histogram drift fit - linear fit
      }
      (*pcstream)<<"robustFit"<<	
	"\n";    
      for (Int_t ihis=0; ihis<12; ihis++){
	delete grDeltaZ[ihis];
	delete grRMSZ[ihis];
	delete fitDeltaZ[ihis];	
      }
    }
    // delete pcstream;  
    //pcstream = new TTreeSRedirector("fitDrift.root","update");
  }
  delete fpol1;
  //
  TTree * treeRobust= (TTree*)(pcstream->GetFile()->Get("robustFit"));  
  Int_t entriesR= treeRobust->Draw("paramRobust.fElements[0]:paramRobust.fElements[1]:paramRobust.fElements[2]:paramRobust.fElements[3]","","goffPara");
  for (Int_t ipar=0; ipar<4; ipar++) {paramRobust[ipar]=TMath::Median(entriesR-5, &(treeRobust->GetVal(ipar)[5]));}
  //
  treeRobust->Draw("paramRobustBC.fElements[0]:paramRobustBC.fElements[1]:paramRobustBC.fElements[2]:paramRobustBC.fElements[3]","","goffPara");
  for (Int_t ipar=0; ipar<4; ipar++) {paramRobustBC[ipar]=TMath::Median(entriesR-5, &(treeRobust->GetVal(ipar)[5]));}
  treeRobust->Draw("paramTRD.fElements[0]:paramTRD.fElements[1]:paramTRD.fElements[2]:paramTRD.fElements[3]","","goffPara");
  for (Int_t ipar=0; ipar<4; ipar++) {paramTRD[ipar]=TMath::Median(entriesR-5, &(treeRobust->GetVal(ipar)[5]));}
  treeRobust->Draw("paramTOF.fElements[0]:paramTOF.fElements[1]:paramTOF.fElements[2]:paramTOF.fElements[3]","","goffPara");
  for (Int_t ipar=0; ipar<4; ipar++) {paramTOF[ipar]=TMath::Median(entriesR-5, &(treeRobust->GetVal(ipar)[5]));}
  treeRobust->Draw("paramTRDBC.fElements[0]:paramTRDBC.fElements[1]:paramTRDBC.fElements[2]:paramTRDBC.fElements[3]","","goffPara");
  for (Int_t ipar=0; ipar<4; ipar++) {paramTRDBC[ipar]=TMath::Median(entriesR-5, &(treeRobust->GetVal(ipar)[5]));}
  treeRobust->Draw("paramTOFBC.fElements[0]:paramTOFBC.fElements[1]:paramTOFBC.fElements[2]:paramTOFBC.fElements[3]","","goffPara");
  for (Int_t ipar=0; ipar<4; ipar++) {paramTOFBC[ipar]=TMath::Median(entriesR-5, &(treeRobust->GetVal(ipar)[5]));}
  paramRobust.Print();
  //
  //
  // 3.) Make drift fit per time interval
  //
  Int_t nTimeBins= Int_t((time1-time0)/deltaT)+1;
  Int_t nParam   = nTimeBins+1;
  TVectorD vecFit(nParam);
  Double_t *pvecFit=vecFit.GetMatrixArray();
  //TLinearFitter *fitterTRD= new TLinearFitter(nParam,TString::Format("hyp%d",nParam+1).Data());
  //TLinearFitter *fitterTOF= new TLinearFitter(nParam,TString::Format("hyp%d",nParam+1).Data());
  //
  TObjArray arrayFit(3);
  TH2F hisTRD("hisTRD","hisTRD",nTimeBins,time0,time1,100,-0.02,0.02);
  TH2F hisTOF("hisTOF","hisTOF",nTimeBins,time0,time1,100,-0.02,0.02);

  for (Int_t iter=0; iter<10; iter++){
    hisTRD.Reset();
    hisTOF.Reset();
    for (Int_t ipoint=0; ipoint<entriesFit; ipoint++){
      if (ipoint%10!=iter) continue;  // points correlated - can be skipped
      Double_t  radius= (*vecR)[ipoint];

      if (radius<kInvalidR) continue;  // no cluster

      Int_t sector   = TMath::Nint((*vecSec)[ipoint]);
      Double_t side  = -1.+((sector%36)<18)*2.;
      Double_t z= (*vecZ)[ipoint];
      double maxZ = kMaxZSect[side<0];
      Double_t drift = (side>0) ? maxZ-(*vecZ)[ipoint] : (*vecZ)[ipoint]+maxZ;
      if (drift>maxZ) continue;
      Double_t phi   = (*vecPhi)[ipoint];
      Double_t gy    = TMath::Sin(phi)*radius;
      pvecFit[0]= side;             // z shift (cm)
      pvecFit[1]= drift*gy/maxZ;   // global y gradient
      pvecFit[2]= drift;            // drift length
      Double_t expected = paramRobust[0]+paramRobust[1]*pvecFit[0]+paramRobust[2]*pvecFit[1]+paramRobust[3]*pvecFit[2];
      Int_t time=(*vecTime)[ipoint];
      //
      Double_t dZTOF=(*deltaZTOF)[ipoint];
      if (dZTOF>kInvalidRes) {
	if (TMath::Abs(dZTOF*side-expected)<kMaxDist1) {
	  dZTOF *= side;
	  //      fitterTOF->AddPoint(pvecFit,dZTOF,1);
	  hisTOF.Fill(time,(dZTOF-expected)/drift,drift/maxZ); 
	}
      }
      Double_t dZTRD=(*deltaZTRD)[ipoint];
      if ( dZTRD>kInvalidRes) {
	if (TMath::Abs(dZTRD*side-expected)<kMaxDist1) {
	  dZTRD*=side;
	  //fitterTOF->AddPoint(pvecFit,dZTRD,1);	
	  hisTRD.Fill(time,(dZTRD-expected)/drift,drift/maxZ); 
	}
      }  

    }
    printf("iter=%d\n",iter);
    TGraphErrors *grTRD=NULL, *grTOF=NULL;
    TGraphErrors *grTRD2=NULL, *grTOF2=NULL;
    Int_t nclTRD=hisTRD.GetEntries();
    Int_t nclTOF=hisTOF.GetEntries();
    if (nclTRD>kMinEntries){
      hisTRD.FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
      grTRD = new TGraphErrors((TH1D*)arrayFit.At(1));
      grTRD2 = new TGraphErrors((TH1D*)arrayFit.At(2));
      arrayFit.Delete();
    }
    if (nclTOF>kMinEntries){
      hisTOF.FitSlicesY(0,0,-1,0,"QNR",&arrayFit);
      grTOF  = new TGraphErrors((TH1D*)arrayFit.At(1));
      grTOF2 = new TGraphErrors((TH1D*)arrayFit.At(2));
      arrayFit.Delete();
    }
    if (grTRD==NULL) {
      grTRD=grTOF; // we should have at minimum one of the histograms not empty
      grTRD2=grTOF2;
    }
    if (grTOF==NULL) {
      grTOF=grTRD;
      grTOF2=grTRD2;
    }
    (*pcstream)<<"fitTime"<<
      "iter="<<iter<<      
      "nclTRD="<<nclTRD<<             // 
      "nclTOF="<<nclTOF<<             // 
      "grTRD.="<<grTRD<<              // time dependent drift correction TRD - mean in time bin
      "grTOF.="<<grTOF<<	      // time dependent drift correction TOF - mean in time bin
      "grTRD2.="<<grTRD2<<            // time dependent drift correction TRD - sigma in time bin
      "grTOF2.="<<grTOF2<<	      // time dependent drift correction TOF - sigma in time bin
      "time0="<<time0<<
      "time1="<<time1<<
      "run="<<runNumber<<
      "paramRobust.="<<&paramRobust<< // time independent parameters
      "paramTRD.="<<&paramTRD<<            //  drift fit using TRD tracks
      "paramTOF.="<<&paramTOF<<	     //  drift fit using TOF tracks
      "paramRobustBC.="<<&paramRobustBC<<  //  drift fit using all tracks with BC
      "paramTRDBC.="<<&paramTRDBC<<        //  drift fit using TRD tracks with BC
      "paramTOFBC.="<<&paramTOFBC<<	     //  drift fit using TOF tracks with BC      
      "\n";
  }
  delete pcstream;
  pcstream = new TTreeSRedirector("fitDrift.root","update"); 
  //
  // 3.) Make local regression and median fit
  //
  TTree * treeFit= (TTree*)(pcstream->GetFile()->Get("fitTime"));  
  Int_t entriesGr = treeFit->Draw("grTRD.fY:grTOF.fY:grTRD.fX:Iteration$","1","goffpara");
  Int_t nbins = TMath::MaxElement(entriesGr, treeFit->GetV4())+1;  

  Double_t dtime0=0,dtime1=0;
  if (entriesGr) TStatToolkit::GetMinMax(treeFit->GetV3(),entriesGr,dtime0,dtime1);
  Int_t ngraphs =entriesGr/nbins;
  //
  // 3.a) local regression fit
  //
  treeFit->Draw("grTRD.fY-grTRD.fY[Iteration$+1]:grTRD.fX","1","goff");
  Double_t deltaRMS,deltaMean;
  AliMathBase::EvaluateUni(entriesGr,treeFit->GetV1(),deltaMean, deltaRMS,entriesGr*0.8);
  TCut cutTRD=TString::Format("abs(grTRD.fY-grTOF.fY)<0.01&&abs(grTRD.fY-grTRD.fY[Iteration$+1]-%f)<%f",deltaMean, 3*deltaRMS).Data();
  TCut cutTOF=TString::Format("abs(grTRD.fY-grTOF.fY)<0.01&&abs(grTOF.fY-grTOF.fY[Iteration$+1]-%f)<%f",deltaMean, 3*deltaRMS).Data();

  THnD *hN = new THnD("hN","hN", 1, &nbins, &dtime0, &dtime1);
  AliNDLocalRegression * pfitTRD = new  AliNDLocalRegression("pfitTRD","pfitTRD");
  AliNDLocalRegression * pfitTOF = new  AliNDLocalRegression("pfitTOF","pfitTOF");
  pfitTRD->SetHistogram((THn*)(hN->Clone()));
  pfitTOF->SetHistogram((THn*)(hN->Clone()));
  pfitTRD->MakeFit(treeFit, "grTRD.fY:grTRD.fEY+0.01", "grTRD.fX",cutTRD, TString::Format("(grTRD.fX/grTRD.fX)+%f",sigmaT),"2:2",0.0001);
  pfitTOF->MakeFit(treeFit, "grTOF.fY:grTOF.fEY+0.01", "grTOF.fX",cutTOF, TString::Format("(grTRD.fX/grTRD.fX)+%f",sigmaT),"2:2",0.0001);
  AliNDLocalRegression::AddVisualCorrection(pfitTRD,104);
  AliNDLocalRegression::AddVisualCorrection(pfitTOF,204);  
  //
  //
  TVectorD medianTRD(nbins);
  TVectorD medianTOF(nbins);
  TVectorD medianQA(nbins);
  TVectorD regTRD(nbins);
  TVectorD regTOF(nbins);
  TVectorD regQA(nbins);
  //
  TVectorD rmsTRD(nbins);
  TVectorD rmsTOF(nbins);
  TVectorD rmsQA(nbins);
  TVectorD vecTimeg(nbins);
  TVectorD vecWorking(entriesGr);
  treeFit->Draw("grTRD.fY:grTOF.fY:grTRD.fX:Iteration$","1","goffpara");
  for (Int_t itype=0; itype<2; itype++){
    for (Int_t itime=0; itime<nbins; itime++){
      for (Int_t igr=0; igr<ngraphs; igr++){
	Int_t index=igr*nbins+itime;
	vecWorking[igr]=treeFit->GetVal(itype)[index];
      }
      Double_t grtime = treeFit->GetVal(2)[itime];
      vecTimeg[itime]=treeFit->GetVal(2)[itime];
      if (itype==0) {
	medianTRD[itime]=TMath::Median(ngraphs,vecWorking.GetMatrixArray()); 
	regTRD[itime]= pfitTRD->Eval(&grtime);
	rmsTRD[itime]=TMath::RMS(ngraphs,vecWorking.GetMatrixArray())/TMath::Sqrt(ngraphs); 
      }
      if (itype==1) {
	medianTOF[itime]=TMath::Median(ngraphs,vecWorking.GetMatrixArray());       
	regTOF[itime]= pfitTOF->Eval(&grtime);
	rmsTOF[itime]=TMath::RMS(ngraphs,vecWorking.GetMatrixArray())/TMath::Sqrt(ngraphs);       
      }
    }
  }
  for (Int_t itime=0; itime<nbins; itime++){  
    medianQA[itime]= medianTRD[itime]-medianTOF[itime];
    regQA[itime]= regTRD[itime]-regTOF[itime];
    rmsQA[itime]=TMath::Sqrt(rmsTRD[itime]*rmsTRD[itime]+rmsTOF[itime]*rmsTOF[itime]);
  }
  TGraphErrors *grmedTRD = new TGraphErrors(nbins, vecTimeg.GetMatrixArray(), medianTRD.GetMatrixArray(),0, rmsTRD.GetMatrixArray());
  TGraphErrors *grmedTOF =  new TGraphErrors(nbins, vecTimeg.GetMatrixArray(), medianTOF.GetMatrixArray(),0, rmsTOF.GetMatrixArray());
  TGraphErrors *grmedQA =  new TGraphErrors(nbins, vecTimeg.GetMatrixArray(), medianQA.GetMatrixArray(),0, rmsQA.GetMatrixArray());
  TGraphErrors *grregTRD = new TGraphErrors(nbins, vecTimeg.GetMatrixArray(), regTRD.GetMatrixArray(),0, rmsTRD.GetMatrixArray());
  TGraphErrors *grregTOF =  new TGraphErrors(nbins, vecTimeg.GetMatrixArray(), regTOF.GetMatrixArray(),0, rmsTOF.GetMatrixArray());
  TGraphErrors *grregQA =  new TGraphErrors(nbins, vecTimeg.GetMatrixArray(), regQA.GetMatrixArray(),0, rmsQA.GetMatrixArray());
  //
  (*pcstream)<<"fitTimeStat"<<
    "grTRDReg.="<<grregTRD<<      // time dependent drift correction TRD - regression estimator
    "grTOFReg.="<<grregTOF<<      // time dependent drift correction TOF - regression estimator
    "grQAReg.="<<grregQA<<        // time dependent drift correction TOF - regression estimator
    "grTRDMed.="<<grmedTRD<<      // time dependent drift correction TRD - median estimator
    "grTOFMed.="<<grmedTOF<<      // time dependent drift correction TOF - median estimator
    "grQAMed.="<<grmedQA<<        // time dependent drift correction TOF - median estimator
    "time0="<<time0<<             
    "time1="<<time1<<
    "run="<<runNumber<<
    "paramRobust.="<<&paramRobust<< // time independent parameters
    "\n";
  delete grmedTRD;
  delete grmedTOF;
  delete grmedQA;
  delete grregTRD;
  delete grregTOF;
  delete grregQA;
  
  delete deltaZTOF;
  delete deltaZTRD;
  delete vecZ;
  delete vecR;
  delete vecSec;
  delete vecPhi;
  delete vecTime;
  delete vecTOFBC;

  delete pcstream;   
  return kTRUE;
}



void  AliTPCcalibAlignInterpolation::MakeNDFit(const char * inputFile, const char * inputTree, Float_t sector0,  Float_t sector1,  Float_t theta0, Float_t theta1){
  //
  /// 
  /// Make ND local regression, QA  for later usage
  /// Parameters:

  // Algorithm:
  //  1.) Make NDLocal regression fits
  //  2.) Make regularization (smoothing)
  //  3.) Make QA trending variable
  //  4.) Make QA plots
  //  5.) Export QA trending variables into trending tree
  //
  /*
    Example usage:
    const char * inputFile="ResidualMapFull_1.root"
    const char * inputTree="deltaRPhiTPCITSTRDDist"
    Float_t sector0=3, sector1 =5;
    Float_t theta0=0, theta1=1;
    AliTPCcalibAlignInterpolation::MakeNDFit(inputFile,inputTree, sector0,sector1, theta0,theta1);
  */

  TTreeSRedirector * pcstream = new TTreeSRedirector(TString::Format("%sFit_sec%d_%d_theta%d_%d.root",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data(),"recreate");
  TTreeSRedirector * pcstreamFit = new TTreeSRedirector(TString::Format("fitTree_%sFit_sec%d_%d_theta%d_%d.root",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data(),"recreate");
  Int_t runNumber=TString(gSystem->Getenv("runNumber")).Atoi();
  TFile * fdist = TFile::Open(inputFile);
  if (!fdist){
    ::Error("MakeNDFit","Intput file %s not accessible\n",inputFile);
    return;
  }
  TTree *treeDist = (TTree*)fdist->Get(inputTree);
  if (!treeDist){
    ::Error("MakeNDFit","Intput tree %s not accessible\n",inputTree);
    return;    
  }
  TTree *treeMeta = (TTree*)fdist->Get("metaData");
  pcstream->GetFile()->cd();
  TTree *treeMetaCopy =  treeMeta->CopyTree("1");
  treeMetaCopy->Write("metaData");
  delete treeMetaCopy;
  //
  // 1.) Make NDLocal regression fits
  //
  const Double_t pxmin=8.48499984741210938e+01; //param.GetPadRowRadii(0,0)-param.GetPadPitchLength(0,0)/2
  const Double_t pxmax=2.46600006103515625e+02; //2.46600006103515625e+02param.GetPadRowRadii(36,param.GetNRow(36)-1)+param.GetPadPitchLength(36,param.GetNRow(36)-1)/2.
  Int_t     ndim=4;
  Int_t     nbins[4]= {30,  (Int_t)((sector1-sector0-0.1)*15),     int(abs(theta1-theta0)*10),        3};  // {radius, phi bin, }
  Double_t  xmin[4] = {pxmin,  sector0+0.05,   theta0,                            -2.0};
  Double_t  xmax[4] = {pxmax, sector1-0.05,   theta1,               2.0};
  //

  THnF* hN= new THnF("exampleFit","exampleFit", ndim, nbins, xmin,xmax);
  treeDist->SetAlias("isOK","rms>0&&vecLTM.fElements[1]*binMedian>0");
  treeDist->SetAlias("delta","vecLTM.fElements[1]");
  TCut cutFit="isOK";
  TCut cutAcceptFit=TString::Format("sectorCenter>%f&&sectorCenter<%f&&kZCenter>%f&&kZCenter<%f", sector0-0.5,sector1+0.5,theta0,theta1).Data();
  TCut cutAcceptDraw=TString::Format("sectorCenter>%f&&sectorCenter<%f&&kZCenter>%f&&kZCenter<%f", sector0,sector1,theta0,theta1).Data();
  
  AliNDLocalRegression *fitCorrs[6]={0};
  for (Int_t icorr=0; icorr<1; icorr++){
    fitCorrs[icorr]= new  AliNDLocalRegression();
    fitCorrs[icorr]->SetName(TString::Format("%sFit%d_sec%d_%d_theta%d_%d",inputTree,icorr, Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data());  
    Int_t hashIndex=fitCorrs[icorr]->GetVisualCorrectionIndex();
    fitCorrs[icorr]->SetHistogram((THn*)(hN->Clone()));  
    TStopwatch timer;
    fitCorrs[0]->SetStreamer(pcstream);
    if (icorr==0) fitCorrs[icorr]->MakeFit(treeDist,"delta:(1+sqrt(abs(qptCenter)))/sqrt(entries)", "RCenter:sectorCenter:kZCenter:qptCenter",cutFit+cutAcceptFit,"3:0.05:0.1:18","2:2:2:2",0.0001);
    timer.Print();
    AliNDLocalRegression::AddVisualCorrection(fitCorrs[icorr]);
    treeDist->SetAlias(TString::Format("delta_Fit%d",icorr).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+0)",hashIndex).Data());
  }
  //
  // 2.) Make regularization (smoothing)
  //    
  Int_t nDims=4;
  Int_t indexes[4]={0,1,2,3};
  Double_t relWeight0[12]={1,4,16,   1,4,16, 1,4,16, 1,4,16};
  Double_t relWeightC[12]={0.5,4,16,   0.5,4,16, 0.5,4,16, 0.5,4,16};
  fitCorrs[1]=(AliNDLocalRegression *)fitCorrs[0]->Clone();
  fitCorrs[1]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight0, 0);
  fitCorrs[2]=(AliNDLocalRegression *)fitCorrs[1]->Clone();
  fitCorrs[2]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight0, 0);
  fitCorrs[3]=(AliNDLocalRegression *)fitCorrs[1]->Clone();
  fitCorrs[4]=(AliNDLocalRegression *)fitCorrs[2]->Clone();
  fitCorrs[3]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeightC, 0, kTRUE);
  fitCorrs[4]->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeightC, 0, kTRUE);
  //
  fitCorrs[1]->SetName(TString::Format("%s_Smooth1",fitCorrs[0]->GetName()).Data());
  fitCorrs[2]->SetName(TString::Format("%s_Smooth2",fitCorrs[1]->GetName()).Data());
  fitCorrs[3]->SetName(TString::Format("%s_SmoothConst1",fitCorrs[0]->GetName()).Data());
  fitCorrs[4]->SetName(TString::Format("%s_SmoothConst2",fitCorrs[1]->GetName()).Data()); 
  for (Int_t icorr=1; icorr<5; icorr++){
    Int_t hashIndex=fitCorrs[icorr]->GetVisualCorrectionIndex();
    AliNDLocalRegression::AddVisualCorrection(fitCorrs[icorr]);
    treeDist->SetAlias(TString::Format("delta_Fit%d",icorr).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+0)",hashIndex).Data());
  }  
  //
  // 3.) Make QA variables and store fits
  //
  gStyle->SetOptFit(1);
  treeDist->Draw(">>drawlist3",cutFit+cutAcceptDraw,"entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("drawlist3");
  treeDist->SetEntryList(elist);

  Double_t quantiles[10]={0.001,0.01,0.05,0.1,0.2, 0.8,0.9,0.99,0.999};
  Double_t deltas[10]={0};  
  const char * pchVar[9]={"delta-delta_Fit0",
			   "delta-delta_Fit1",
			   "delta-delta_Fit2",
			   "delta-delta_Fit3",
			   "delta-delta_Fit4",
			   "delta_Fit0-delta_Fit1",
			   "delta_Fit0-delta_Fit2",
			   "delta_Fit0-delta_Fit3",
			   "delta_Fit0-delta_Fit4"};
  const char * pchTittle[9]={"map-fit_{0} (cm)",
			      "map-fit_{0Reg1} (cm)",
			      "map-fit_{0Reg2} (cm)",
			      "map-fit_{0Reg1Const} (cm)",
			      "map-fit_{0Reg2Const} (cm)",
			      "fit_{O}-fit_{0Reg1} (cm)",
			      "fit_{0}-fit_{0Reg2} (cm)",
			      "fit_{0}-fit_{0Reg1Const} (cm)",
			      "fit_{0}-fit_{0Reg2Const} (cm)"};
  TGraph *    grQuantiles[18]={0};  
  TVectorD    vecRMS1(18);
  for (Int_t idiff=0; idiff<18; idiff++){
    Int_t chEntries=0;
    if (idiff<9){
      chEntries=treeDist->Draw(pchVar[idiff%9],"1","goff");
    }
    if (idiff>=9){
      chEntries=treeDist->Draw(TString::Format("(%s)/(rms/sqrt(entries))",pchVar[idiff%9]).Data(),"1","goff");
    }
    for (Int_t iq=0; iq<10; iq++){
      deltas[iq]=TMath::KOrdStat(chEntries,treeDist->GetV1(), Int_t(chEntries*quantiles[iq]));
    }
    grQuantiles[idiff]=new TGraph(10,quantiles,deltas);
    grQuantiles[idiff]->SetTitle(pchTittle[idiff%9]);
    Double_t mean, rms=0;
    AliMathBase::EvaluateUni(chEntries,treeDist->GetV1(),mean, rms,0.80*chEntries);
    vecRMS1[idiff]=rms;
  }
  vecRMS1.Print();
  TH1* his=0;
  TVectorD slopePols1(5);
  for (Int_t i=0; i<5; i++){
    treeDist->Draw(TString::Format("delta:delta_Fit%d>>hisQA2D%d",i,i).Data(),cutFit+cutAcceptDraw,"colzgoff"); 
    his=treeDist->GetHistogram();
    his->Fit("pol1");
    slopePols1[i]= his->GetFunction("pol1")->GetParameter(1);
  }

  TFile * fout = pcstream->GetFile();
  pcstream->GetFile()->cd();
  for (Int_t iter=0; iter<5; iter++){
    if (TString(gSystem->Getenv("AliTPCcalibAlignInterpolation_MakeNDFit_KeepCovariance")).Atoi()==0){
      if (iter!=0) fitCorrs[iter]->CleanCovariance();
    }
    fitCorrs[iter]->Write();
    if (TString(gSystem->Getenv("AliTPCcalibAlignInterpolation_MakeNDFit_DumpFitTree")).Atoi()>0){
      fitCorrs[iter]->DumpToTree(4, (*pcstreamFit)<<TString::Format("tree%s", fitCorrs[iter]->GetName()).Data());
    }
  }
  //
  // 4.) Make standard QA plot   
  //
  gStyle->SetLabelSize(0.06,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1200,1000);
  canvasQA->Divide(1,3);
  //
  // 4.1) delta:fit correaltion
  canvasQA->cd(1)->SetLogz();
  treeDist->Draw("delta:delta_Fit0>>hisQA2D",cutFit+cutAcceptDraw,"colzgoff");
  his=treeDist->GetHistogram();
  his->GetXaxis()->SetTitle("#Delta_{fit} (cm)");
  his->GetYaxis()->SetTitle("#Delta_{map} (cm)");
  his->Fit("pol1");
  Double_t slopePol1= his->GetFunction("pol1")->GetParameter(1);
  his->Write("hisQA2D");
  his->Draw("colz");
  // 
  // 4.2) delta-fitReg0 and fitReg0-fitReg1
  canvasQA->cd(2)->SetLogy();
  treeDist->SetLineColor(1); treeDist->SetMarkerColor(1);
  treeDist->Draw("delta-delta_Fit0>>hisQA1D(300)",cutFit+cutAcceptDraw,"goff");
  his=treeDist->GetHistogram();  
  his->GetXaxis()->SetTitle("#Delta_{map}-#Delta_{fit} (cm)");
  his->Fit("gaus");
  Double_t mean=his->GetMean();
  Double_t rms=his->GetRMS();
  Double_t rmsG= his->GetFunction("gaus")->GetParameter(2);
  //
  treeDist->GetHistogram()->Write("hisQA1D");
  treeDist->SetLineColor(2); treeDist->SetMarkerColor(2);
  treeDist->Draw("delta_Fit0-delta_Fit1>>hisQA1DifFit(300)",cutFit+cutAcceptDraw,"same");
  his=treeDist->GetHistogram();
  Double_t meanDiffFit=his->GetMean();
  Double_t rmsDiffFit=his->GetRMS();
  treeDist->GetHistogram()->Write("hisQA1DifFit");
  //
  // 4.3) (delta-fitReg0) "pull" 
  treeDist->SetLineColor(1); treeDist->SetMarkerColor(1);
  canvasQA->cd(3)->SetLogy();
  treeDist->Draw("(delta_Fit0-delta)/(rms/sqrt(entries))>>hisQAPull",cutFit+cutAcceptDraw,"goff");
  his=treeDist->GetHistogram();
  his->Fit("gaus");
  Double_t meanPull=his->GetMean();
  Double_t rmsPull=his->GetRMS();
  Double_t rmsPullG= his->GetFunction("gaus")->GetParameter(2);
  treeDist->GetHistogram()->Write("hisQAPull");
  treeDist->SetLineColor(2); treeDist->SetMarkerColor(2);
  treeDist->Draw("(delta_Fit0-delta_Fit1)/(rms/sqrt(entries))>>hisQAPullDifFit",cutFit+cutAcceptDraw,"same");
  his=treeDist->GetHistogram();
  Double_t meanPullDiffFit=his->GetMean();
  Double_t rmsPullDiffFit=his->GetRMS();
  treeDist->GetHistogram()->Write("hisQAPullDiffFit");

  canvasQA->SaveAs((TString::Format("%sFit_sec%d_%d_theta%d_%dQA.png",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data()));

  TCanvas *canvasQAFit = new TCanvas("canvasQAFit","canvasQAFit",1200,1000); 
  canvasQAFit->SetRightMargin(0.01);
  canvasQAFit->Divide(1,5,0,0); 
  treeDist->SetMarkerStyle(25);
  treeDist->SetMarkerSize(0.5);
  TCut defaultCut="abs(qptCenter)<0.1";
  //
  {
    canvasQAFit->cd(1)->SetRightMargin(0.1);
    treeDist->Draw("delta:sectorCenter:RCenter",defaultCut+"abs(abs(kZCenter)-0.1)<0.06","colz");
    canvasQAFit->cd(2)->SetRightMargin(0.1);
    treeDist->Draw("delta-delta_Fit0:sectorCenter:RCenter",defaultCut+"abs(abs(kZCenter)-0.1)<0.06","colz");
    canvasQAFit->cd(3)->SetRightMargin(0.1);
    treeDist->Draw("delta-delta_Fit1:sectorCenter:RCenter",defaultCut+"abs(abs(kZCenter)-0.1)<0.06","colz");
    canvasQAFit->cd(4)->SetRightMargin(0.1);
    treeDist->Draw("delta-delta_Fit2:sectorCenter:RCenter",defaultCut+"abs(abs(kZCenter)-0.1)<0.06","colz");
    canvasQAFit->cd(5)->SetRightMargin(0.1);
    treeDist->Draw("delta_Fit0-delta_Fit2:sectorCenter:RCenter",defaultCut+"abs(abs(kZCenter)-0.1)<0.06","colz");
  }
  canvasQAFit->SaveAs((TString::Format("%sFit_sec%d_%d_theta%d_%dQAFit.png",inputTree,Int_t(sector0),Int_t(sector1),Int_t(theta0),Int_t(theta1)).Data()));
  TObjString input=inputTree;
  //
  // 5.) Export trending variables - used for validation of the fit
  //
  treeMeta->Draw("runNumber:selHis:fillCounter:clusterCounter:meanTime:ntracksUsed:startTime:stopTime","1","goffpara");
  (*pcstream)<<"NDFitTrending"<<               // cp of subset of info from meta data (rest accessible in the metadata tree also avaialble in file)
    "run="<<treeMeta->GetVal(0)[0]<<      
    "selHis="<<treeMeta->GetVal(1)[0]<<
    "fillCounter="<<treeMeta->GetVal(2)[0]<<
    "clusterCounter="<<treeMeta->GetVal(3)[0]<<
    "meanTime="<<treeMeta->GetVal(4)[0]<<
    "ntracksUsed="<<treeMeta->GetVal(5)[0]<<
    "startTime="<<treeMeta->GetVal(6)[0]<<
    "stopTime="<<treeMeta->GetVal(7)[0];
  Int_t entriesCl= treeMeta->Draw("grNcl72.fY:grNcl73.fY:grNcl74.fY:grNcl75.fY","(grNcl72.fX>startTime&&grNcl72.fX<stopTime&&grNcl72.fY!=0)","goffpara");
  TVectorF vecNcl(8);
  if (entriesCl>0) {
    for (Int_t icl=0; icl<4; icl++){
      vecNcl[icl]=TMath::Median(entriesCl, treeMeta->GetVal(icl));
    }
    entriesCl= treeMeta->Draw("grNclUsed72.fY:grNclUsed73.fY:grNclUsed74.fY:grNclUsed75.fY","(grNcl72.fX>startTime&&grNcl72.fX<stopTime&&grNcl72.fY!=0)","goffpara");
    for (Int_t icl=0; icl<4; icl++){
      vecNcl[icl+4]=TMath::Median(entriesCl, treeMeta->GetVal(icl));
    }
  }
  (*pcstream)<<"NDFitTrending"<<               // cp of subset of info from meta data (rest accessible in the metadata tree also avaialble in file)
    "vecNclCounter.="<<&vecNcl;
  //
  (*pcstream)<<"NDFitTrending"<<
    "input.="<<&input<<                // name of the input file
    "inputTree="<<inputTree<<          // name of the input tree
    "runNumber="<<runNumber<<          // run number ID
    "sec0="<<sector0<<                 // range: sector0 
    "sec1="<<sector1<<                 // range: sector1
    "theta0="<<theta0<<                // range: theta0
    "theta1="<<theta1<<                // range: theta1 
    //                                 // QA plots statistical info
    "slopePols1.="<<&slopePols1<<      // data:fit - pol1 fit for differnt Regularization
    "slopePol1="<<slopePol1<<
    //
    "mean="<<mean<<                    // mean <value-fitND0>
    "rms="<<rms<<                      // rms  (value-fitND0>
    "rmsG="<<rmsG<<                    // gaus fit rms  (value-fitND0>
    "meanPull="<<meanPull<<            // normalized mean <value-fitND0>
    "rmsPull="<<rmsPull<<              // normalized error<value-fitND0>
    "rmsPullG="<<rmsPull<<             // gaus fit normalized error<value-fitND0>
    "meanDiffFit="<<meanDiffFit<<      //  
    "rmsDiffFit="<<rmsDiffFit<<
    "meanPullDiffFit="<<meanPullDiffFit<<
    "rmsPullDiffFit="<<rmsPullDiffFit<<
    "vecRMS1.="<<&vecRMS1;  // rms of the diffs of different estimators (Kernles, Regularization)
  for (Int_t iq=0; iq<18; iq++){
    (*pcstream)<<"NDFitTrending"<<
      TString::Format("grQuantiles%d.=",iq).Data()<<grQuantiles[iq];
  }
  (*pcstream)<<"NDFitTrending"<<"\n";
  delete treeMeta;
  delete pcstream;
  delete pcstreamFit;  
}


TTree* AliTPCcalibAlignInterpolation::LoadDistortionTrees(const char * maplist, Int_t cacheSize, Int_t markerStyle, Float_t markerSize){
  //
  // Load distortion trees specified in the maplist 
  // Loading distortion maps as a friend trees used for
  //    - correction for reference distortions (e.g map at low IR)
  //    - distortion maps correlation studies
  //    - distortion maps scaling fitting
  // 
  // To obtain run number and TimeBin ID we assume naming convention as used in the calibration procedure. 
  // This naming convention is hardwired in the code. 
  //
  TTree * treeReturn=0;
  TTree * tree=0;
  TObjArray* array = TString(gSystem->GetFromPipe(TString::Format("cat %s",maplist).Data())).Tokenize("\n");  
  TObjArray*arrayOK=new TObjArray(array->GetEntries());
  for (Int_t i=0; i<array->GetEntries(); i++){
    printf("%s\n",array->At(i)->GetName());
    TString fname(array->At(i)->GetName());
    Int_t index=fname.Index("/000");
    TString runName(&(fname[index+1]),9);
    index=fname.Index("/Time");
    TString timeString(&(fname[index+9]),4);    
    if (TString(array->At(i)->GetName()).Contains("_0.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(treeReturn,array->At(i)->GetName(),"deltaRPhiTPCITSDist",runName+"_"+timeString+"ITSY");
    }
    if (TString(array->At(i)->GetName()).Contains("_1.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(treeReturn,array->At(i)->GetName(),"deltaRPhiTPCITSTRDDist",runName+"_"+timeString+"TRDY");
    }
    if (TString(array->At(i)->GetName()).Contains("_2.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(treeReturn,array->At(i)->GetName(),"deltaRPhiTPCITSTOFDist",runName+"_"+timeString+"TOFY");
    }
    if (TString(array->At(i)->GetName()).Contains("_3.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(treeReturn,array->At(i)->GetName(),"deltaZTPCITSDist",runName+"_"+timeString+"ITSZ");
    }
    if (TString(array->At(i)->GetName()).Contains("_4.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(treeReturn,array->At(i)->GetName(),"deltaZTPCITSTRDDist",runName+"_"+timeString+"TRDZ");
    }
    if (TString(array->At(i)->GetName()).Contains("_5.root")){
      tree = AliTPCcalibAlignInterpolation::AddFriendDistortionTree(treeReturn,array->At(i)->GetName(),"deltaZTPCITSTOFDist",runName+"_"+timeString+"TOFZ");
    }
    if (tree) {
      arrayOK->AddLast(array->At(i));
      treeReturn=tree;
    }    
  } 
  treeReturn->SetCacheSize(cacheSize);
  treeReturn->SetMarkerStyle(markerStyle); 
  treeReturn->SetMarkerSize(markerSize);
  return treeReturn;
}


Bool_t  AliTPCcalibAlignInterpolation::LoadNDLocalFit(TTree * tree, const char *chTree){
  ///
  ///  Load ND local fits. We assume data are organized in particular directory structure, in directories together with input maps
  ///  
  ///   
  /// Input:    
  ///   \param TTree * tree        - input tree with distortion maps and "friend trees" per time bins
  ///   \param const char *chTree  - name of the "distortion branch"

  if ( tree->GetListOfFriends()->FindObject(chTree)==NULL){
    ::Error("AliTPCcalibAlignInterpolation::LoadNDLocal","Tree %s does not exist",chTree);
    return kFALSE;
  }
  TString floc = tree->GetListOfFriends()->FindObject(chTree)->GetTitle();
  TString fdir = gSystem->DirName(floc);
  //
  TObjArray *ndFileList = ( gSystem->GetFromPipe(TString::Format("ls %s/delta*root",fdir.Data()).Data())).Tokenize("\n");
  
  if (ndFileList->GetEntries()==0){
    ::Error(" AliTPCcalibAlignInterpolation::LoadNDLocal","File with NDLocal %s",chTree);
    return kFALSE;
  }
  for (Int_t ind=0; ind<ndFileList->GetEntries(); ind++){
    //
    TString  fname=ndFileList->At(ind)->GetName();
    TFile * f = TFile::Open(fname.Data());
    TList * arrKey = f->GetListOfKeys();
    for (Int_t ikey=0; ikey<arrKey->GetEntries(); ikey++){
      TString keyName=arrKey->At(ikey)->GetName();
      if (keyName.Contains("delta")==0) continue;
      TObject * o = f->Get(keyName);
      AliNDLocalRegression * reg  = dynamic_cast<AliNDLocalRegression*>(o);
      if (reg==NULL){
	delete o;
	continue;
      }
      TString aliasName = TString::Format("%s.%s",chTree,reg->GetName());
      aliasName.ReplaceAll("-","Min");
      ::Info("AliTPCcalibAlignInterpolation::LoadNDLocal","Loaded ND local regression %s/%s as alias %s", fname.Data(), reg->GetName(),aliasName.Data());
      reg->SetName(aliasName);
      Int_t hashIndex=reg->GetVisualCorrectionIndex();
      reg->AddVisualCorrection(reg, hashIndex);      
      tree->SetAlias(aliasName, TString::Format("AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+0)",hashIndex).Data());
      tree->SetAlias(aliasName+"L", TString::Format("AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+0)-(AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter-2)+AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+2))*0.5",hashIndex,hashIndex,hashIndex).Data());
      tree->SetAlias(aliasName+"M", TString::Format("(AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter-2)+AliNDLocalRegression::GetCorrND(%d,RCenter,sectorCenter,kZCenter,qptCenter+2))*0.5",hashIndex, hashIndex).Data());
    }
  }
  return kTRUE;

}


void  AliTPCcalibAlignInterpolation::DrawMapEstimatorComparison(TTree * tree, const char* chtree,  Float_t radius, Float_t kZ, TCut &  selection, const char *figType){
  // Predefined plot: 
  //    Draw distortion map comparison
  //    Compare median and LTM estimator of the mean value of distortion in the bin
  //
  
  // Example usage:
  /* 
     const char* chtree="000244918_his1TRDY";    // Low IR one polarity
     const char* chtree="000246391_his1TRDY";    // Low IR another polarity
     Float_t radius=100;
     Float_t kZ=0.1;
     figType="png"; 
     AliTPCcalibAlignInterpolation::DrawMapEstimatorComparison(tree, chtree, radius,kZ,"abs(qptCenter)<0.1",figType);
  */
  if (!tree) {
    ::Error("DrawEstimatorComparison","Tree not available");
    return;
  }
  if (chtree && tree->GetListOfFriends()->FindObject(chtree)==NULL){
    ::Error("DrawEstimatorComparison","Ttree %s not available",chtree);
    return;
  }
  gStyle->SetLabelSize(0.07,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(0.6,"Y");
  gStyle->SetTitleOffset(0.4,"Z");
  gStyle->SetOptTitle(1);


  TCanvas * canvasC = new TCanvas("canvasC","canvasC",1400,1000);
  TPad * pad=0;
  TPad *pads[4]={0};
  for (Int_t ipad=0; ipad<4; ipad++){
    pads[ipad] = new TPad("pad1","This is pad1",0.00,ipad/4.,  1,(ipad+1.)/4.);
    pads[ipad]->SetBottomMargin(0);
    pads[ipad]->SetTopMargin(0);
  }
  pads[0]->SetBottomMargin(0.15);
  pads[3]->SetTopMargin(0.05);
  for (Int_t ipad=0; ipad<4; ipad++){
    pads[ipad]->Draw();
    pads[ipad]->SetGrid(1,1);
  }

  if (chtree){
    pads[3]->cd();
    tree->Draw(TString::Format("%s.binMedian:sectorCenter:RCenter",chtree).Data(),selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&%s.rms>0",kZ,radius,chtree).Data(),"colz");
    pads[2]->cd();
    tree->Draw(TString::Format("%s.vecLTM.fElements[1]:sectorCenter:RCenter",chtree).Data(),selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&%s.rms>0",kZ,radius,chtree).Data(),"colz");
    pads[1]->cd();
    tree->Draw(TString::Format("%s.meanG:sectorCenter:RCenter",chtree).Data(),selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&%s.rms>0&&abs(%s.meanG-%s.binMedian)<1",kZ,radius,chtree,chtree,chtree).Data(),"colz");
    pads[0]->cd();
    tree->Draw(TString::Format("%s.vecLTM.fElements[1]-%s.binMedian:sectorCenter:RCenter",chtree,chtree).Data(),selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&%s.rms>0",kZ,radius,chtree).Data(),"colz");
  }else{
    pads[3]->cd();
    tree->Draw("binMedian:sectorCenter:RCenter",selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&rms>0",kZ,radius).Data(),"colz");
    pads[2]->cd();
    tree->Draw("vecLTM.fElements[1]:sectorCenter:RCenter",selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&rms>0",kZ,radius).Data(),"colz");
    pads[1]->cd();
    tree->Draw("meanG:sectorCenter:RCenter",selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&rms>0&&abs(meanG-binMedian)<1",kZ,radius).Data(),"colz");
    pads[0]->cd();
    tree->Draw("binMedian-vecLTM.fElements[1]:sectorCenter:RCenter",selection+TString::Format("abs(kZCenter+(%2.2f))<0.06&&abs(RCenter-%2.2f)<20&&rms>0",kZ,radius).Data(),"colz");
  }
  
  if (figType) canvasC->SaveAs(TString::Format("distortionMap_%s_R%2.0f_kZ%2.2f.%s",chtree,radius,kZ,figType).Data());
}



Bool_t  AliTPCcalibAlignInterpolation::DrawScalingComparison(TTree * tree, const char* chRef, const char *chBin0, const char *chBin1,  Float_t R0, Float_t R1, Float_t kZ, const char *figType){
  ///
  /// Make predefined plot: 
  ///    Draw distortion maps comparison and save figure (figType) in current working directory.
  ///       Fig 1.): Map run/timeBin chBin0 corrected for reference map (chRef)offset 
  ///       Fig 2.): Map run/timeBin chBin1 corrected for reference map (chRef)offset  
  ///       Fig 3.): Map (chBin1-chRef)-scale*(chBin0-chRef)
  /// 
  /// Input:    
  ///   \param TTree * tree - input tree with distortion maps and "friend trees" per time bins
  ///   \param chRef        - reference distortion map
  ///   \param chBin0       - run or time bin shown in firs row
  ///   \param chBin1       - run or time bin 
  /// \return kTRUE if comparison of distortion map possible and figure saved
  ///    TString::Format("distortionMapScaling_%s_%s_%s_R0%2.0f_R1%2.0f_kZ%2.2f.%s",chBin0,chBin1,chRef,R0,R1,kZ,figType).Data();
  /// Example usage:
  /// To be used for systematic studies of TPC distortion. 
  /// E.g:   
  /*
    TTree * tree = AliTPCcalibAlignInterpolation::LoadDistortionTrees("map.list");
    const char* chRef="000246391_his1TRDY";       // Low IR refernce run for given B-field polarity
    const char* chBin0="000246980_his1TRDY";      // time bin at the beginnig of the fill
    const char* chBin1="000246994_his1TRDY";      // time bin at the end of the fill
    AliTPCcalibAlignInterpolation::DrawScalingComparison(tree, chRef, chBin0, chBin1, 85, 130, 0.1, "png");
  */
  //
  // 0.) Check variables
  //
  TCut defaultCut="abs(qptCenter)<0.1";
  if (tree==NULL) {
    ::Error("AliTPCcalibAlignInterpolation::DrawScalingComparison","Tree not set");
    return kFALSE;
  }
  Int_t tentries[3]={0};
  const char *vars[3]={chRef, chBin0, chBin1};
  for (Int_t i=0; i<3;i++){
    tentries[i]=tree->Draw(TString::Format("%s.binMedian",vars[i]),defaultCut,"goff", 10000);
    if (tentries[i]<=0){
      ::Error("AliTPCcalibAlignInterpolation::DrawScalingComparison","Expression %s  or tree %s not valid ", vars[i], tree->GetName());
      return kFALSE;
    }
  }
  //
  // 1.) get scaling factor using A side scaling (OROC scaling on C side more problematic)
  //
  Int_t entries = tree->Draw(TString::Format("%s.binMedian-%s.binMedian:%s.binMedian-%s.binMedian", chBin0,chRef, chBin1,chRef).Data(),defaultCut+"kZCenter>0","goff");
  TGraph * gr0 = new TGraph(entries,tree->GetV1(),tree->GetV2());
  TGraph * gr1 = new TGraph(entries,tree->GetV2(),tree->GetV1());
  gr0->Fit("pol1");
  gr1->Fit("pol1");
  Double_t slope=(gr0->GetFunction("pol1")->GetParameter(1)+ 1/gr1->GetFunction("pol1")->GetParameter(1))*0.5;  
  tree->SetAlias("norm0",TString::Format("%f*(%s.binMedian-%s.binMedian)",slope,chBin0,chRef).Data());
  //
  // 2. Make plots 
  //
  gStyle->SetLabelSize(0.07,"XYZ");
  gStyle->SetLabelSize(0.03,"Y");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(0.5,"Y");
  gStyle->SetTitleOffset(0.4,"Z");
  gStyle->SetOptTitle(1);
  //
  TCanvas * canvasC = new TCanvas("canvasC","canvasC",1400,1000);
  TPad * pad=0;
  TPad *pad3 = new TPad("pad1","This is pad1",0.00,0.0,  1,0.33);
  TPad *pad2 = new TPad("pad2","This is pad2",0.00,0.33, 1,0.66);
  TPad *pad1 = new TPad("pad3","This is pad3",0.00,0.66, 1,1);
  pad1->SetBottomMargin(0);
  pad2->SetBottomMargin(0);
  pad3->SetBottomMargin(0.15);
  pad2->SetTopMargin(0);
  pad3->SetTopMargin(0);
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad1->SetGrid(1,1);
  pad2->SetGrid(1,1);
  pad3->SetGrid(1,1);
  TCut cutAccept=defaultCut+TString::Format("abs(kZCenter+(%2.2f))<0.06&&RCenter>%2.0f&&RCenter<%2.0f&&%s.rms>0",kZ,R0,R1,chRef).Data();
  Int_t isOK=0;
  const Int_t kMinEntries=100;
  {
    pad1->cd();
    isOK+=tree->Draw(TString::Format("(%s.binMedian-%s.binMedian):sectorCenter:RCenter", chBin0,chRef),cutAccept,"colz")>kMinEntries;
    pad2->cd();
    isOK+=tree->Draw(TString::Format("(%s.binMedian-%s.binMedian):sectorCenter:RCenter", chBin1,chRef),cutAccept,"colz")>kMinEntries;
    pad3->cd();
    isOK+=tree->Draw(TString::Format("(%s.binMedian-%s.binMedian)-norm0:sectorCenter:RCenter", chBin1,chRef),cutAccept,"colz")>kMinEntries;
  }
  if (isOK<3){
    ::Error("AliTPCcalibAlignInterpolation::DrawScalingComparison","Not enough points in selected region R<%2.0f,%2.0f>", R0,R1);
    return kFALSE;
  }
  if (figType) canvasC->SaveAs(TString::Format("distortionMapScaling_%s_%s_%s_R0%2.0f_R1%2.0f_kZ%2.2f.%s",chBin0,chBin1,chRef,R0,R1,kZ,figType).Data());
  return kTRUE;
}



//_______________________________________________
double AliTPCcalibAlignInterpolation::GetTgPhi(double x, double y2x, double q2p, double b)
{
  // calculate tangent of primary track at any frame at given x,y
  double y = y2x*x;
  double c = q2p*b*(-0.299792458e-3);
  if (TMath::Abs(c)<1e-9) return y2x;
  double r2 = y*y+x*x;
  double det = 4./r2 - c*c;
  double snp;
  if (det<0) {
    snp = TMath::Sign(-0.8,c);
    //printf("track of q2p=%f cannot reach x:%f y:%f, define snp as %f \n",q2p,x,y,snp);
  }
  else {
    snp = 0.5*(y*TMath::Sqrt(det)-c*x); // snp at vertex
    snp += x*c;  // snp at x,y
  }
  return snp/TMath::Sqrt((1-snp)*(1+snp));
}

//______________________________________________
void AliTPCcalibAlignInterpolation::FixAlignmentBug(int sect, float q2pt, float bz,
						    float& alp, float& x, float &z, float &deltaY, float &deltaZ)
{
  // fix alignment bug: https://alice.its.cern.ch/jira/browse/ATO-339?focusedCommentId=170850&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-170850
  //
  // NOTE: deltaZ in the buggy code is calculated as Ztrack_with_bug - Zcluster_w/o_bug
  static TGeoHMatrix *mCache[72] = {0};
  if (sect<0||sect>=72) {
    AliErrorClassF("Invalid sector %d",sect);
    return;
  }
  int lr = sect/36 ? (AliGeomManager::kTPC2) : (AliGeomManager::kTPC1);
  TGeoHMatrix* mgt = mCache[sect];
  if (!mgt) {
    int volID = AliGeomManager::LayerToVolUIDSafe(lr,sect%36);
    mgt = new TGeoHMatrix(*AliGeomManager::GetTracking2LocalMatrix(volID));
    mgt->MultiplyLeft(AliGeomManager::GetMatrix(volID));
    mCache[sect] = mgt;
    printf("Caching matrix for sector %d\n",sect);
  }  
  double alpSect = ((sect%18)+0.5)*20.*TMath::DegToRad();

  // cluster in its proper alpha frame with alignment bug
  double xyzClUse[3] = {x,0,z}; // this is what we read from the residual tree
  double xyzTrUse[3] = {x, deltaY, z}; // track in bad cluster frame
  //
  // recover cluster Z position by adding deltaZ
  double zClSave = xyzClUse[2] -= deltaZ;  // here the cluster is not affected by Z alignment component of the bug!
  static AliExternalTrackParam trDummy;
  trDummy.Local2GlobalPosition(xyzClUse,alp); // misaligned cluster in global frame
  double xyz0[3]={xyzClUse[0],xyzClUse[1],xyzClUse[2]};
  mgt->MasterToLocal(xyz0,xyzClUse);
  // we got ideal cluster in the sector tracking frame,  but now the Z is wrong, since it was not affected by the bug!!!
  //
  xyzClUse[2] = zClSave;

  // go to ideal cluster frame
  trDummy.Local2GlobalPosition(xyzClUse,alpSect); // ideal global
  double alpFix = TMath::ATan2(xyzClUse[1],xyzClUse[0]);    // fixed cluster phi
  trDummy.Global2LocalPosition(xyzClUse,alpFix);     // fixed cluster in in its frame
  //
  trDummy.Local2GlobalPosition(xyzTrUse,alp); // track in global frame
  trDummy.Global2LocalPosition(xyzTrUse,alpFix); // track in cluster frame
  alp = alpFix;
  //
  double dx = xyzTrUse[0] - xyzClUse[0]; // x might not be the same after alignment fix
  // deduce track slopes assuming it comes from the vertex
  double tgphi = GetTgPhi(xyzClUse[0],xyzTrUse[1]/xyzClUse[0],q2pt,bz);
  xyzTrUse[1] -= dx*tgphi;
  xyzTrUse[2] -= dx*xyzClUse[2]/xyzClUse[0]; // z2x
  //
  x = xyzClUse[0];
  z = xyzTrUse[2]; // we still use track Z as a reference ...
  deltaY = xyzTrUse[1]-xyzClUse[1];
  deltaZ = xyzTrUse[2]-xyzClUse[2];
  //
}


void  AliTPCcalibAlignInterpolation::MakeVDriftOCDB(const char *inputFile, Int_t run, TString  targetOCDBstorage, const char * testDiffCDB/*=0*/){
  ///
  /// Make OCDB entry using information using the vdrift calibration obtained in the  AliTPCcalibAlignInterpolation::FitDrift function (suing ResidualTrees.root)
  /// 
  ///
  /* 
     char * inputFile= "fitDrift.root"
     char * testDiffCDB="/cvmfs/alice.cern.ch/calibration/data/2015/OCDB/TPC/Calib/TimeDrift/Run244918_244918_v3_s0.root";
     Int_t  run=244918;

     .x  $ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C(run,"local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB/")
     AliTPCcalibAlignInterpolation::MakeVDriftOCDB( inputFile,run,"", testDiffCDB)
  */

  /*
   The same calibration was done previoulsy using the AliTPCcalibTime class
   Many graphs in previously done RUN1 calibration used for monitoring purposes and do not have equivalent in the new calibration 
   using ResidualTrees  
   Old CPass0 calibration () - 114 graphs exported and monitored later using trending  (O (100 kBy) per run)
      * 108 graphs for Align+drift = 9 params x 3 detectors x 4 version ) 
      ** vdrift calibration (3 param)  and alingment graphs (6)
      ** each pair of detector ITS->TPC, TRD->TPC, TOF->TPC
      ** Kalman filter storted in 4 version
      *** raw input per time bin, kalman forward propagation, Kalaman back propagation, smoothed Kalman  
      ** usage 
      *** ITS->TPC smoothed version in case exist
      *** TOF->TPC smoothed in case ITS->TPC not availale
        
      * 6 graphs for the Laser CE backup graphs
      * backup of ClusterParam and RecoParam as it was used in the process of  calibration 
      ** calibration depends
  
   New calibration input data ( O(3 MBy) for 2 hours measurement):
   
   Export variables:
     only vdrift calibration graphs (_DELTAZ,DRIFTVD, TO and VDGY)  will  be exported
     smoothed version equivalent (grTRDReg) and raw calibration equivalent (grTRDmed)
     ITS  - mean (TRD+TOF)*0.5
     TRD  - TRD
     TOF  - TOF
  */  
  TObjArray * driftArray = new TObjArray();

  //
  // 0. Initialize OCDB if not done before
  //
  AliTPCParam *param =0;
  if (AliCDBManager::Instance()->GetDefaultStorage()!=NULL){
    //
    //
    AliTPCClusterParam *clParam =   AliTPCcalibDB::Instance()->GetClusterParam();
    param= AliTPCcalibDB::Instance()->GetParameters();
    TObjArray *recoParams = new TObjArray(4) ;
    for (Int_t i=0;i<4;i++) recoParams->AddAt(AliTPCcalibDB::Instance()->GetRecoParam(i),i);
    driftArray->AddLast(clParam);
    driftArray->AddLast(recoParams);    
    ::Info("AliTPCcalibAlignInterpolation::MakeVDriftOCDB","Residual calibration used. We have to trust  - your OCDB is correct, as we can nt check.");
  }else{
    ::Fatal("AliTPCcalibAlignInterpolation::MakeVDriftOCDB","Im sorry. OCDB has to be intilaized before. We need to get Parameters objects, which were used in previous calibration. In ideal case the OCDB entry to be setup to the entries, indeed used to produce Residual trees ... But this we can not CHECK");
  }

  //
  // 1. Read calibratio data from the tree
  //
  TVectorD     *vdriftParam=0;
  TGraphErrors *grTRDReg=0;  
  TGraphErrors *grTOFReg=0;  
  TGraphErrors *grTRDMed=0;  
  TGraphErrors *grTOFMed=0;  
  TFile *fdrift = TFile::Open(inputFile);
  if (fdrift){
    TTree * tree = (TTree*)fdrift->Get("fitTimeStat");
    if (tree==NULL){
      ::Fatal("MakeVDriftOCDB.LoadDriftCalibration FAILED", "tree fitTimeStat not avaliable in file fitDrift.root");
    }else{
      if (tree->GetBranch("grTRDReg."))  tree->SetBranchAddress("grTRDReg.",&grTRDReg);
      if (tree->GetBranch("grTOFReg."))  tree->SetBranchAddress("grTOFReg.",&grTOFReg);
      if (tree->GetBranch("grTRDMed."))  tree->SetBranchAddress("grTRDMed.",&grTRDMed);
      if (tree->GetBranch("grTOFMed."))  tree->SetBranchAddress("grTOFMed.",&grTOFMed);
      if (tree->GetBranch("paramRobust.")) tree->SetBranchAddress("paramRobust.",&vdriftParam);
      tree->GetEntry(0);
    }    
  }else{
    ::Fatal("MakeVDriftOCDB.LoadDriftCalibration FAILED", "Input file %s not accesible", inputFile);
  }  
  if (grTRDReg==NULL && grTOFReg==NULL){
    ::Fatal("MakeVDriftOCDB.LoadDriftCalibration FAILED", "drift calibration enot present");
  }
  //
  // 2.)  Transform v drift calibration object into the format used in previous calibration
  //
  TGraphErrors   *graphDELTAZ=0;
  TGraphErrors   *graphT0=0;
  TGraphErrors   *graphVDGY=0;
  TGraphErrors   *graphDRIFTVD=0;
  TGraphErrors   *graph=0;

  const char *chDet[3]={"ITS","TRD","TOF"};
  Double_t atime[2]={0,0};
  Double_t deltaZ[2]={0,0};
  Double_t t0[2]={0,0};
  Double_t vdgy[2]={0,0};
  graphDRIFTVD=(grTRDMed!=NULL) ? grTRDMed:grTOFMed;
  //
  atime[0]=graphDRIFTVD->GetX()[0];
  atime[1]=graphDRIFTVD->GetX()[graphDRIFTVD->GetN()-1];
  for (Int_t ipoint=0; ipoint<=1; ipoint++){
    deltaZ[ipoint]=-(*vdriftParam)[1];  // unit OK
    vdgy[ipoint]=-(*vdriftParam)[2];   // units OK
    t0[ipoint]=-(*vdriftParam)[0]/(1+(*vdriftParam)[3]);       // t0 to be normalized to the ms
    t0[ipoint]/=(param->GetDriftV()/1000000.); 
  }
  

  Double_t vdrifCorrRun=1./(1.+(*vdriftParam)[3]);

  //
  for (Int_t idet=0; idet<3; idet++){
    for (Int_t itype=0; itype<2; itype++){  //0. original version 1.)smoothed version
      TString grPrefix="ALIGN_";
      grPrefix+=chDet[idet];
      if (itype==0) grPrefix+="_TPC_";
      if (itype==1) grPrefix+="B_TPC_";
      //
      graphDELTAZ=new TGraphErrors(2, atime, deltaZ);
      graphDELTAZ->SetName(grPrefix+"DELTAZ");
      driftArray->AddLast(graphDELTAZ);
      graphT0=new TGraphErrors(2, atime, t0);
      graphT0->SetName(grPrefix+"T0");
      driftArray->AddLast(graphT0);
      graphVDGY=new TGraphErrors(2, atime, vdgy);
      graphVDGY->SetName(grPrefix+"VDGY");
      driftArray->AddLast(graphVDGY);
      //
      // drift velocity
      //graph=(grTRDMed!=NULL) ? new TGraphErrors(*grTRDMed): new TGraphErrors(*grTOFMed); // normal constructor give us seg fault
      graph=(TGraphErrors*)((grTRDMed!=NULL) ? grTRDMed->Clone(): grTOFMed->Clone()); // normal constructor give us 0
      Int_t npoints = graph->GetN();
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
	graph->GetY()[ipoint]=-1;
	if (idet==1 && grTRDMed) { // TRD
	  if (itype==0) graph->GetY()[ipoint]=(vdrifCorrRun*(1-grTRDMed->GetY()[ipoint]))-1.;
	  if (itype==1) graph->GetY()[ipoint]=(vdrifCorrRun*(1-grTRDReg->GetY()[ipoint]))-1.;
	}
	if (idet==2 && grTOFMed) { // TOF
	  if (itype==0) graph->GetY()[ipoint]=(vdrifCorrRun*(1-grTOFMed->GetY()[ipoint]))-1.;
	  if (itype==1) graph->GetY()[ipoint]=(vdrifCorrRun*(1-grTOFReg->GetY()[ipoint]))-1.;
	}
	if (idet==0 &&  (grTRDMed&&grTOFMed)) { // (TRD+TOF)*0.5
	  if (itype==0) graph->GetY()[ipoint]=(vdrifCorrRun*(1-(grTRDMed->GetY()[ipoint]+grTOFMed->GetY()[ipoint])*0.5))-1.;
	  if (itype==1) graph->GetY()[ipoint]=(vdrifCorrRun*(1-(grTRDReg->GetY()[ipoint]+grTOFReg->GetY()[ipoint])*0.5))-1.;
	}
      }
      graph->SetName(grPrefix+"DRIFTVD");
      driftArray->AddLast(graph);
    }
  }
  //
  // 3. Store selected graphs in OCDB
  //
  AliCDBStorage* targetStorage = 0x0;
  if (targetOCDBstorage.Length()==0) {
    targetOCDBstorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    targetStorage = AliCDBManager::Instance()->GetStorage(targetOCDBstorage.Data());
  }
  else if (targetOCDBstorage.CompareTo("same",TString::kIgnoreCase) == 0 ){
    targetStorage = AliCDBManager::Instance()->GetDefaultStorage();
  }
  else {
    targetStorage = AliCDBManager::Instance()->GetStorage(targetOCDBstorage.Data());
  }

  AliCDBMetaData* metaData = new AliCDBMetaData;
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion(">v5-07-20"); //AliRoot version
  metaData->SetComment("AliTPCcalibAlignInterpolation Calibration of the time dependence of the drift velocity using Residual trees");
  AliCDBId id1("TPC/Calib/TimeDrift", run, run);
  
  //now the entry owns the metadata, but NOT the data
  AliCDBEntry *driftCDBentry=new AliCDBEntry(driftArray,id1,metaData,kFALSE);
  targetStorage->Put(driftCDBentry); 
  //
  //
  // 4. Make diff to reference CDB  if sepecfied
  // 
  if (testDiffCDB){    
    TFile *f = TFile::Open(testDiffCDB);
    AliCDBEntry * entry=(AliCDBEntry*)f->Get("AliCDBEntry");
    TObjArray * array = (TObjArray*)entry->GetObject();
    //array->ls();
    TGraphErrors * grRef= (TGraphErrors * )array->FindObject("ALIGN_TRDB_TPC_DRIFTVD");
    TGraphErrors * gr= (TGraphErrors * )driftArray->FindObject("ALIGN_TRDB_TPC_DRIFTVD");
    grRef->SetMarkerStyle(25);
    gr->SetMarkerStyle(21);
    Double_t minimum=TMath::Min(grRef->GetMinimum(), gr->GetMinimum());
    Double_t maximum=TMath::Max(grRef->GetMaximum(), gr->GetMaximum());
    grRef->SetMinimum(minimum);
    grRef->SetMaximum(maximum);

    grRef->Draw("ap");
    gr->Draw("p");
  }

  delete driftCDBentry;


}


Float_t  AliTPCcalibAlignInterpolation::CalculateDistance(const TVectorF &track0, const TVectorF &track1, const TVectorF &vecSec, TVectorF &vecDelta, Int_t npValid, Float_t &rmsTrack,  Float_t &rmsCluster, Float_t lpNorm){
  ///  Calculate track and cluster RMS distance. RMS distance is than used in toder to
  ///    1.) Calcute rms distance between the track0 and reference track1    
  ///    2.) Calculate rms distance between the cluster and local median cluster position in region +-2 padrows (wihing one sector)
  ///        and store delta in exported vecDelta array
  /// \param[in] track0             vector: cluster-track residuals for track of interest
  /// \param[in] track1             vector: cluster-track residuals for reference track
  /// \param[in] vecSec             vector: cluster sector index
  /// \param[in] lpNorm             the p value for the p-type norm (https://en.wikipedia.org/wiki/Norm_(mathematics))
  /// \param[out] rmsTrack          calculated 2 track distance (lpNorm)
  /// \param[out] rmsCluster        calculated cluster local mean distance (lpNorm)
  ///
  ///  Return value - number of the clusters with too big RMS
  /// 
  // In the following routine AliTPCcalibAlignInterpolation::FillHistogramsFromChain cut applied on quality of track and custer infromation
  // Suggested cuts (based on the resolution parametrrization studies)
  //    - nOut <10 - less than 15 clusters rejected   - kMaxSkippedCluster=10;
  //    - rmsTrack<2 cm    kMaxRMSTrackCut=2.0;
  //    - rmsCluster<0.3 cm kMaxRMSClusterCut=0.3;  
  //    - |vecDelta|<0.5 cm  kMaxDeltaClusterCut=0.5; 
  //
  // Parameters of algorithm for the moment set as a constant 
  const Float_t   kMinFractionPoints=0.5;
  const Float_t kMaxDist=20;
  const Float_t kMaxDistTrack=5;
  Float_t maxRMSTrack=rmsTrack;
  Float_t maxRMSCluster=rmsCluster;

  //
  Float_t distance2=0;
  Int_t npoints=0;
   // cache matrix pointers
  const Float_t *delta0=track0.GetMatrixArray();
  const Float_t *delta1=track1.GetMatrixArray();
  const Float_t *sec=vecSec.GetMatrixArray();
  //
  // 1.) Calcute rms distance between the track0 and reference track1    
  //
  for (Int_t irow=0; irow<npValid; irow++){
    vecDelta[irow]=1000;
    if (TMath::Abs(delta0[irow])==0) continue;  // 
    if (TMath::Abs(delta1[irow])==0) continue;  // 
    if (TMath::Abs(delta0[irow])>kMaxDist) continue;  // 
    if (TMath::Abs(delta1[irow])>kMaxDist) continue;  //
    if (TMath::Abs(delta0[irow]-delta1[irow])>kMaxDistTrack) continue;  //
    npoints++;
    distance2+=(delta0[irow]-delta1[irow])*(delta0[irow]-delta1[irow]);
  }
  if (npoints<=80) return -1;
  if (npoints<npValid*kMinFractionPoints) return -1; // not enough points
  rmsTrack = TMath::Sqrt(distance2/npoints);
  if (rmsTrack>maxRMSTrack) {
    return -1;
  }
  //
  //    2.) Calculate rms distance between the cluster and local median cluster posotion in region +-2 padrows (wihing one sector)
  //         and store delta in exported vecDelta array
  Int_t nOutCounter=0;  // counter of clusters 
  rmsCluster=0;
  Float_t rmsSum=0;
  Float_t deltaArray[10]={0};
  for (Int_t irow=0; irow<npValid; irow++){
    Int_t idelta=1;
    vecDelta[irow]=1000;
    if (TMath::Abs(delta0[irow])>kMaxDist) continue;
    deltaArray[0]=delta0[irow];
    Int_t nforMedian=1;
    Int_t sector36=Int_t(sec[irow])%36;
    for (idelta=1; idelta<3; idelta++) {
      if ((irow-idelta)<0) break;
      if ((irow+idelta)>=npValid) break;
      if (TMath::Abs(sec[irow-idelta])>72) break;
      if (TMath::Abs(sec[irow+idelta])>72) break;
      if (Int_t(sec[irow-idelta])%36!=sector36) break;  // sometimes looks like sector not initialized
      if (Int_t(sec[irow+idelta])%36!=sector36) break;
      if (TMath::Abs(delta0[irow-idelta])<kMaxDist) deltaArray[nforMedian++]=delta0[irow-idelta];
      if (TMath::Abs(delta0[irow+idelta])<kMaxDist) deltaArray[nforMedian++]=delta0[irow+idelta];
    }
    if (irow==0) vecDelta[irow]=delta0[0]-delta0[1];
    if (irow==npValid-1) vecDelta[irow]=delta0[irow]-delta0[irow-1];
    if (nforMedian>1){
      vecDelta[irow]=delta0[irow]-TMath::Mean(nforMedian, deltaArray);
    }
    if (vecDelta[irow]>maxRMSCluster){
      nOutCounter++;
    }else{
      rmsCluster+=TMath::Power(TMath::Abs(vecDelta[irow]),lpNorm);
      rmsSum++;
    }
  }
  if (rmsSum<=0) return 0;
  rmsCluster=TMath::Power(rmsCluster/rmsSum,1./lpNorm);
  return nOutCounter;
}


Float_t AliTPCcalibAlignInterpolation::InitForAlignmentBugFix(int run, const char* ocdb)
{
  ::Info(" AliTPCcalibAlignInterpolation::InitForAlignmentBugFix","Alignment bug fix is requested\n");
  //
  // this requires the field and the geometry ...
  if (run<1) ::Fatal("tstw","InitForBugFix: Run number is not provided");
  Bool_t geomOK = AliGeomManager::GetGeometry() != 0;
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!geomOK || !fld) { // need to setup ocdb?
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man->IsDefaultStorageSet()) man->SetDefaultStorage(ocdb);
    if (man->GetRun()!=run) man->SetRun(run);
  }
  if (!geomOK) {
    AliGeomManager::LoadGeometry();
    AliGeomManager::ApplyAlignObjsFromCDB("TPC");
  }
  if (!fld) {
    AliGRPManager grpMan;
    grpMan.ReadGRPEntry();
    grpMan.SetMagField();
    fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  }
  return fld->SolenoidField();
}
