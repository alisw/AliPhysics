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
//#include "THnSparse.h"
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

#include "AliTPCcalibAlignInterpolation.h"

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
  fStreamer(0),         // calibration streamer 
  fStreamLevel(0),      // stream level - In mode 0 only basic information needed for calibration  stored (see EStream
  fSyswatchStep(1),      // dump system resource information after  fSyswatchStep tracks
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
  fStreamer(0),         // calibration streamer 
  fStreamLevel(0),      // stream level - In mode 0 only basic information needed for calibration  stored (see EStream
  fSyswatchStep(1),      // dump system resource information after  fSyswatchStep tracks  
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


Bool_t  AliTPCcalibAlignInterpolation::RefitITStrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackITS, Double_t &chi2, Double_t &npoints){
  //
  // Makes a refit of the ITS track
  // Input: AliESDfriendTrack, particle mass, outer ITS TrackParam 
  // Output: AliExternalTrackParam of the ITS track refit at the last layer of ITS
  //
  const Double_t kMaxRadius=50;
  Int_t sortedIndex[AliTPCcalibAlignInterpolation_kMaxPoints]={0};
  if (friendTrack->GetITSOut()==NULL) return kFALSE;  
  Double_t sigma[2] = {0.003,0.01};    // ITS intrincsic resolution in (y,z)  - error from the points can be used SD layer 2-3 sighnificantly bigger error
  // !!!! We should set ITS error parameterization form outside !!!!
  //
  trackITS = *((AliExternalTrackParam*)friendTrack->GetITSOut());
  // Reset track to the vertex
  AliTrackerBase::PropagateTrackToBxByBz(&trackITS,0.,mass,1,kFALSE);
  trackITS.ResetCovariance(1000.);
  
  // Get space points
  AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
  if (!pointarray){
    printf("Space points are not stored in the friendTrack!\n");
    return kFALSE;
  };
  Int_t nPoints = pointarray->GetNPoints();  // # space points of all available detectors                                            
                                             // Sort space points first
  SortPointArray(pointarray, sortedIndex);  // space point indices sorted by radius in increasing order
  
  // Propagate track through ITS space points
  AliTrackPoint spacepoint;
  chi2=0;
  npoints=0; 
  Int_t volId=0,modId=0,layerId=0;
  
  for (Int_t iPoint=0;iPoint<nPoints;iPoint++){
    pointarray->GetPoint(spacepoint,sortedIndex[iPoint]);
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    trackITS.Global2LocalPosition(xyz,alpha);
    if (xyz[0]>kMaxRadius) break;  // use only ITS points - maybe we should indexes of elements
    trackITS.Rotate(alpha);
    Bool_t status = AliTrackerBase::PropagateTrackToBxByBz(&trackITS,xyz[0],mass,1,kFALSE);
    if (status){
      Double_t pos[2] = {xyz[1], xyz[2]};
      Double_t cov[3] = {sigma[0]*sigma[0],0, sigma[1]*sigma[1]};   
      volId = spacepoint.GetVolumeID();
      layerId = AliGeomManager::VolUIDToLayer(volId,modId);
      if (layerId ==AliGeomManager::kSDD1 || layerId ==AliGeomManager::kSDD2) {
        cov[0]*=16.; cov[2]*=16.;
      }      
      double chi2cl = trackITS.GetPredictedChi2(pos,cov);
      chi2 += chi2cl;
      npoints++;
      trackITS.Update(pos,cov);      
    }
  }
  return kTRUE;
}


Bool_t AliTPCcalibAlignInterpolation::RefitTRDtrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackTRD, Double_t &chi2, Double_t &npoints){
  //
  // Makes a refit of the TRD track
  // Input: AliESDfriendTrack, particle mass, inner TRD TrackParam 
  // Output: AliExternalTrackParam of the TRD track refit - in the first layer of TRD
  // Here we forgot about the tiliting pads of TRD - we assume delta Z and delta y are uncorelated
  //      given approximation is in average tru - in case avearaging of significantly bigger than pad length
  //  
  Double_t sigma[2] = {0.03, 5};
  Int_t sortedIndex[AliTPCcalibAlignInterpolation_kMaxPoints]={0};
  const Double_t kMaxRadius=370;
  const Double_t kMinRadius=280;
  if (friendTrack->GetTRDIn() == NULL) return kFALSE;
  trackTRD = *((AliExternalTrackParam*)friendTrack->GetTRDIn());
  
  
  // Reset track outside TRD
  AliTrackerBase::PropagateTrackToBxByBz(&trackTRD,370.,mass,1,kFALSE);
  trackTRD.ResetCovariance(1000.);
  
  // Get space points
  AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
  if (!pointarray){
    printf("Space points are not stored in the friendTrack!\n");
    return kFALSE;
  };
  Int_t nPoints = pointarray->GetNPoints();  // # space points of all available detectors
                                             // Sort space points first
  SortPointArray(pointarray, sortedIndex);  // space point indices sorted by radius in increasing order

  // Propagate track through TRD space points
  AliTrackPoint spacepoint;
  chi2=0; 
  npoints=0;
  for (Int_t iPoint=nPoints-1;iPoint>=0;iPoint--){
    pointarray->GetPoint(spacepoint,sortedIndex[iPoint]);
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    trackTRD.Global2LocalPosition(xyz,alpha);
    if (xyz[0]>kMaxRadius) continue;  // use only TRD points
    if (xyz[0]<kMinRadius) break;  // use only TRD points
    trackTRD.Rotate(alpha);
    Bool_t status = AliTrackerBase::PropagateTrackToBxByBz(&trackTRD,xyz[0],mass,1,kFALSE);
    if (status){
      Double_t pos[2] = {xyz[1], xyz[2]};
      Double_t cov[3] = {sigma[0]*sigma[0],0,sigma[1]*sigma[1]};
      double chi2cl = trackTRD.GetPredictedChi2(pos,cov);
      chi2 += chi2cl;
      npoints++;
      trackTRD.Update(pos,cov);
    }
  }
  return kTRUE;
}


Bool_t  AliTPCcalibAlignInterpolation::RefitTOFtrack(AliESDfriendTrack *friendTrack, Double_t mass, AliExternalTrackParam &trackTOF, Double_t &chi2, Double_t &npoints){
  //
  // Makes a refit of the TRD track
  // Input: AliESDfriendTrack, particle mass, OUTER ITS track - propagated to the TOF point and updated by TOF point 
  // Output: AliExternalTrackParam of the TOF track refit - at the TOF point
  Double_t sigma[2] = {1., 1.};      // should be parameterized
  const Double_t kTOFRadius = 370;
  Int_t sortedIndex[AliTPCcalibAlignInterpolation_kMaxPoints]={0};
  AliTrackerBase::PropagateTrackToBxByBz(&trackTOF,kTOFRadius,mass,10,kTRUE);
  Int_t volId,modId,layerId;
  
  // Get space points
  AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
  if (!pointarray){
    printf("Space points are not stored in the friendTrack!\n");
    return kFALSE;
  };
  Int_t nPoints = pointarray->GetNPoints();  // # space points of all available detectors
                                             // Sort space points first
  SortPointArray(pointarray, sortedIndex);  // space point indices sorted by radius in increasing order
  // Propagate track through TRD space points
  AliTrackPoint spacepoint;
  chi2=0; 
  npoints=0;
  for (Int_t iPoint=nPoints-1;iPoint>=0;iPoint--){  
    pointarray->GetPoint(spacepoint,sortedIndex[iPoint]);
    volId = spacepoint.GetVolumeID();
    layerId = AliGeomManager::VolUIDToLayer(volId,modId);
    if (layerId !=AliGeomManager::kTOF) continue;
    
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    trackTOF.Global2LocalPosition(xyz,alpha);
    trackTOF.Rotate(alpha);
    Bool_t status = AliTrackerBase::PropagateTrackToBxByBz(&trackTOF,xyz[0],mass,1,kFALSE);
    if (status){
      Double_t pos[2] = {xyz[1], xyz[2]};
      Double_t cov[3] = {sigma[0]*sigma[0],0,sigma[1]*sigma[1]};
      double chi2cl = trackTOF.GetPredictedChi2(pos,cov);
      chi2 += chi2cl;
      npoints++;
      trackTOF.Update(pos,cov);
    }
  }
  return kTRUE;  
}




Bool_t  AliTPCcalibAlignInterpolation::SortPointArray(AliTrackPointArray *pointarray, Int_t * sortedIndex){
  //
  // Fill array of indexes to the pointArray (array sorted in increasing order)
  //
  if (sortedIndex==NULL) return kFALSE;
  Int_t nPoints = pointarray->GetNPoints();
  Double_t xunsorted[nPoints];
  AliTrackPoint spacepoint;
  AliExternalTrackParam param;
  for (Int_t iPoint=0;iPoint<nPoints;iPoint++){
    pointarray->GetPoint(spacepoint,iPoint);
    Double_t xyz[3] = {(Double_t)spacepoint.GetX(),(Double_t)spacepoint.GetY(),(Double_t)spacepoint.GetZ()};
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    param.Global2LocalPosition(xyz,alpha);
    xunsorted[iPoint] = xyz[0];
  }
  if (nPoints==0) {
    return kFALSE;
  }
  TMath::Sort(nPoints,xunsorted,sortedIndex,kFALSE);
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
  AliESDfriend *esdFriend=static_cast<AliESDfriend*>(esdEvent->FindListObject("AliESDfriend"));
  if (!esdFriend) return;
  if (esdFriend->TestSkipBit()) return;
  Int_t nPrimTracks= (esdEvent->GetPrimaryVertex()!=NULL)? esdEvent->GetPrimaryVertex()->GetNContributors():0;
  Int_t nTracks = esdEvent->GetNumberOfTracks();  // Get number of tracks in ESD
  if (nTracks==0) return;
  if (!fStreamer) fStreamer = new TTreeSRedirector("ResidualHistos.root","recreate");
  //
  const Int_t nPointsAlloc=AliTPCcalibAlignInterpolation_kMaxPoints; 
  const Int_t kMaxLayer=159;
  AliExternalTrackParam trackArrayITS[kMaxLayer];
  AliExternalTrackParam trackArrayTRD[kMaxLayer];
  AliExternalTrackParam trackArrayTOF[kMaxLayer];
  AliExternalTrackParam trackArrayITSTRD[kMaxLayer];
  AliExternalTrackParam trackArrayITSTOF[kMaxLayer];
  //
  //MakeResidualHistosInterpolation();
  //
  Int_t sortedIndex[AliTPCcalibAlignInterpolation_kMaxPoints];
  TVectorF deltaITS0(kMaxLayer), deltaITS1(kMaxLayer); 
  TVectorF deltaTRD0(kMaxLayer), deltaTRD1(kMaxLayer); 
  TVectorF deltaTOF0(kMaxLayer), deltaTOF1(kMaxLayer); 
  TVectorF vecR(kMaxLayer), vecPhi(kMaxLayer), vecZ(kMaxLayer);
  
  for (Int_t iTrack=0;iTrack<nTracks;iTrack++){ // Track loop
    // 0.) For each track in each event, get the AliESDfriendTrack
    AliESDtrack *esdTrack = esdEvent->GetTrack(iTrack);
    AliESDfriendTrack *friendTrack = esdFriend->GetTrack(iTrack);
    if (!friendTrack) continue;      
    if (esdTrack->GetITSNcls()<4) continue;
    Double_t mass = esdTrack->GetMass();  // particle mass    
    Double_t tofBC=esdTrack->GetTOFBunchCrossing();
    // Get TPC seed
    TObject *calibObject=0;
    AliTPCseed *seed = 0;
    for (Int_t l=0;(calibObject=friendTrack->GetCalibObject(l));++l) {
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) break;
    }
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
    // 2.) ITS, TRD and ITS-TRD refits
    //
    RefitITStrack(friendTrack,mass,paramITS, itsChi2, itsNCl );
    if (itsNCl<4) continue; 
    RefitTRDtrack(friendTrack,mass,paramTRD, trdChi2, trdNCl); 
    paramTOF=paramITS;
    RefitTOFtrack(friendTrack,mass,paramTOF, tofChi2, tofNCl );
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("Refitting",fTrackCounter,1,0,0);        
    //if ((trdNCl+tofNCl)==0) continue; // - use ITS only tracks also  -should it be option?
    //
    // 3.) Propagate to TPC volume, histogram residuals to TPC clusters and dump all information to TTree
    //
    AliTrackPoint spacepoint;  
    Int_t volId,modId,layerId;      
    AliTrackPointArray *pointarray = (AliTrackPointArray*)friendTrack->GetTrackPointArray();
    if (!pointarray) continue;
    Int_t nPointsAll = pointarray->GetNPoints();  // # space points of all available detectors
    SortPointArray(pointarray, sortedIndex);  // space point indices sorted by radius in increasing order
    fTrackCounter++; // increase total track number
    //
    // 4.) Propagate  ITS tracks outward
    // 
    //  
    for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++){
      AliTPCclusterMI *cluster=seed->GetClusterPointer(iPoint);
      if (cluster==NULL) continue;
      Float_t fxyz[3] = {0};
      cluster->GetGlobalXYZ(fxyz);
      Double_t xyz[3]={fxyz[0],fxyz[1],fxyz[2]};
      Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
      paramITS.Global2LocalPosition(xyz,alpha);	
      paramITS.Rotate(alpha);
      Bool_t itsOK = AliTrackerBase::PropagateTrackToBxByBz(&paramITS,xyz[0],mass,1,kFALSE);
      if (itsOK){
	trackArrayITS[iPoint]=paramITS;
	trackArrayITS[iPoint].SetUniqueID(1);
      }
    }
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("ExtrapolateITS",fTrackCounter,2,0,0);  
    //
    // 5.) Propagate  TRD/TOF tracks inwards
    //
    for (Int_t iPoint=kMaxLayer-1;iPoint>=0;iPoint--){
      AliTPCclusterMI *cluster=seed->GetClusterPointer(iPoint);
      if (cluster==NULL) continue;
      Float_t fxyz[3] = {0};
      cluster->GetGlobalXYZ(fxyz);
      Double_t xyz[3]={fxyz[0],fxyz[1],fxyz[2]};
      Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);            
      if (trdNCl>0){
	paramTRD.Global2LocalPosition(xyz,alpha);	
	paramTRD.Rotate(alpha);
	Bool_t trdOK = AliTrackerBase::PropagateTrackToBxByBz(&paramTRD,xyz[0],mass,1,kFALSE);
	if (trdOK){
	  trackArrayTRD[iPoint]=paramTRD;
	  trackArrayITSTRD[iPoint]=paramTRD;
	  AliTrackerBase::UpdateTrack(trackArrayITSTRD[iPoint], trackArrayITS[iPoint]);	  
	}
      }
      if (tofNCl>0){
	paramTOF.Global2LocalPosition(xyz,alpha);	
	paramTOF.Rotate(alpha);
	Bool_t tofOK = AliTrackerBase::PropagateTrackToBxByBz(&paramTOF,xyz[0],mass,1,kFALSE);
	if (tofOK){
	  trackArrayTOF[iPoint]=paramTOF;
	  trackArrayITSTOF[iPoint]=paramTOF;
	  trackArrayTOF[iPoint].SetUniqueID(1);
	  trackArrayITSTOF[iPoint].SetUniqueID(1);
	  AliTrackerBase::UpdateTrack(trackArrayITSTOF[iPoint], trackArrayITS[iPoint]);	  
	}
      }
    }
    if (fTrackCounter%fSyswatchStep==0) AliSysInfo::AddStamp("InterpolateTRD",fTrackCounter,3,0,0);  

    if ( ((fStreamLevel&kStremInterpolation)>0)&&(fTrackCounter%fSyswatchStep==0)){
      for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++){
	AliTPCclusterMI *cluster=seed->GetClusterPointer(iPoint);
	if (cluster==NULL) continue;
	(*fStreamer)<<"interpolation"<<
          "itrack="<<fTrackCounter<<  // total track #
          "cluster.="<<cluster<<  // space points                                    //
          "itsNCl="<<itsNCl<<
          "trdNCl="<<trdNCl<<
          "tofNCl="<<tofNCl<<
          //
          "itsChi2="<<itsChi2<<
          "trdChi2="<<trdChi2<<
          "tofChi2="<<tofChi2<<
	  "tofBC="<<tofBC<<
          //
          "trackITS.="<<&trackArrayITS[iPoint]<<  // ITS fit
          "trackTRD.="<<&trackArrayTRD[iPoint]<<  // TRD fit
          "trackTOF.="<<&trackArrayTOF[iPoint]<<  // TOF fit
          "trackITSTRD.="<<&trackArrayITSTRD[iPoint]<<  // ITS-TRD fit
          "trackITSTOF.="<<&trackArrayITSTOF[iPoint]<<  // ITS-TOF fit
          "\n";	
      }
    }
    Int_t counter=0;
    Double_t rounding=100;
    memset( deltaITS0.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaITS1.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTRD0.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTRD1.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTOF0.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    memset( deltaTOF1.GetMatrixArray(), 0,kMaxLayer*sizeof(Float_t));
    
    for (Int_t iPoint=0;iPoint<kMaxLayer;iPoint++){      
      AliTPCclusterMI *cluster=seed->GetClusterPointer(iPoint);
      if (cluster==NULL) continue;
      deltaITS0[counter]= Int_t(trackArrayITS[iPoint].GetY()*rounding)/rounding;
      deltaITS1[counter]= Int_t((cluster->GetZ()-trackArrayITS[iPoint].GetZ())*rounding)/rounding;
      deltaTRD0[counter]=0;
      deltaTRD1[counter]=0;
      deltaTOF0[counter]=0;
      deltaTOF1[counter]=0;
      if (trackArrayTRD[iPoint].GetUniqueID()>0){
	deltaTRD0[counter]= Int_t(trackArrayITSTRD[iPoint].GetY()*rounding)/rounding;
	deltaTRD1[counter]= Int_t((cluster->GetZ()-trackArrayITSTRD[iPoint].GetZ())*rounding)/rounding;
      }
      if (trackArrayTOF[iPoint].GetUniqueID()>0){
	deltaTOF0[counter]= Int_t(trackArrayITSTOF[iPoint].GetY()*rounding)/rounding;
	deltaTOF1[counter]= Int_t((cluster->GetZ()-trackArrayITSTOF[iPoint].GetZ())*rounding)/rounding;
      }
      //      vecR(kMaxLayer), vecPhi(kMaxLayer), vecZ(kMaxLayer);
      vecR[counter]=trackArrayITS[iPoint].GetX();
      vecPhi[counter]=trackArrayITS[iPoint].GetAlpha();
      vecZ[counter]=trackArrayITS[iPoint].GetZ();
      counter++;
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
      "itsChi2="<<itsChi2<<
      "trdChi2="<<trdChi2<<
      "tofChi2="<<tofChi2<<
      "tofBC="<<tofBC<<
      //
      "track.="<<ip<<                    // track parameters at inner wal of TPC
      "vecR.="<<&vecR<<
      "vecPhi.="<<&vecPhi<<
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
  //
 // end of track loop
}

void AliTPCcalibAlignInterpolation::CreateResidualHistosInterpolation(Double_t dy, Double_t dz){
  //
  // Make cluster residual histograms
  //
  Double_t xminTrack[9], xmaxTrack[9];
  Double_t xminTrackITS[9], xmaxTrackITS[9];
  Int_t    binsTrack[9], binsTrackITS[9];
  TString  axisName[9],axisTitle[9];
  //
  // 0 - delta   of interest
  // 1 - global  phi in sector number  as float
  // 2 - local   r
  // 3 - local   kz
  // 4 - local   q/pt
  // 
  // gx,gy,gz - will be taken from the TPC
  //
  axisName[0]="delta";   axisTitle[0]="#Delta (cm)";                 // to fill    local(clusterY-track.y)
  binsTrack[0]=90;       xminTrack[0]=-dy;        xmaxTrack[0]=dy; 
  binsTrackITS[0]=90;    xminTrackITS[0]=-dy;     xmaxTrackITS[0]=dy; 
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
  binsTrackITS[3]=36;    xminTrackITS[3]=-1.8;     xmaxTrackITS[3]=1.8;  // +-1.8 for the ITS 
  //
  axisName[4]="qpt";    axisTitle[4]="q/pt (c/GeV)";                         // to fill : track.GetSigned1Pt() 
  binsTrack[4]=5;        xminTrack[4]=-2.5;        xmaxTrack[4]=2.5; 
  binsTrackITS[4]=5;     xminTrackITS[4]=-2.5;     xmaxTrackITS[4]=2.5; 


  //
  fHisITSDRPhi = new THnF("deltaRPhiTPCITS","#Delta_{Y} (cm)", 5, binsTrackITS,xminTrackITS, xmaxTrackITS);
  fHisITSTRDDRPhi = new THnF("deltaRPhiTPCITSTRD","#Delta_{Y} (cm) TPC-(ITS+TRD)", 5, binsTrack,xminTrack, xmaxTrack);
  fHisITSTOFDRPhi = new THnF("deltaRPhiTPCITSTOF","#Delta_{Y} (cm) TPC-(ITS+TOF)", 5, binsTrack,xminTrack, xmaxTrack);
  //
  fHisITSDZ = new THnF("deltaZTPCITS","#Delta_{Z} (cm)", 5, binsTrackITS,xminTrackITS, xmaxTrackITS);
  fHisITSTRDDZ = new THnF("deltaZTPCITSTRD","#Delta_{Z} (cm) TPC-(ITS+TRD)", 5, binsTrack,xminTrack, xmaxTrack);
  fHisITSTOFDZ = new THnF("deltaZTPCITSTOF","#Delta_{Z} (cm) TPC-(ITS+TOF)", 5, binsTrack,xminTrack, xmaxTrack);
  //
  //
  //
  for (Int_t ivar2=0;ivar2<5;ivar2++){
    fHisITSDRPhi->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    fHisITSDRPhi->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    fHisITSTRDDRPhi->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    fHisITSTRDDRPhi->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    fHisITSTOFDRPhi->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    fHisITSTOFDRPhi->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    fHisITSDZ->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    fHisITSDZ->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    fHisITSTRDDZ->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    fHisITSTRDDZ->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
    fHisITSTOFDZ->GetAxis(ivar2)->SetName(axisName[ivar2].Data());
    fHisITSTOFDZ->GetAxis(ivar2)->SetTitle(axisName[ivar2].Data());
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

void    AliTPCcalibAlignInterpolation::FillHistogramsFromChain(const char * residualList, Double_t dy, Double_t dz, Int_t downscale){
  /**
   * Trees with point-track residuals to residual histogram
   * @param residualList  text file with tree list
   * Output - ResidualHistograms.root file with hitogram within AliTPCcalibAlignInterpolation object
   */
  //
  //
  // 
  const Int_t knPoints=159;
  AliTPCcalibAlignInterpolation * calibInterpolation = new  AliTPCcalibAlignInterpolation("calibInterpolation","calibInterpolation",kFALSE);
  calibInterpolation->CreateResidualHistosInterpolation(dy,dz);
  TString branches[6]={"its0.","trd0.","tof0.", "its1.","trd1.","tof1."};
  //
  TVectorF *vecDelta= 0;
  TVectorF *vecR=0;
  TVectorF *vecPhi=0;
  TVectorF *vecZ=0;
  AliExternalTrackParam *param = 0;
  //
  TString  esdList0 = gSystem->GetFromPipe(TString::Format("cat %s",residualList).Data());
  TObjArray *esdArray= esdList0.Tokenize("\n");  
  Int_t nesd = esdArray->GetEntriesFast();  
  //
  THn *hisToFill[6]={calibInterpolation->GetHisITSDRPhi(), calibInterpolation->GetHisITSTRDDRPhi(),  calibInterpolation->GetHisITSTOFDRPhi(), calibInterpolation->GetHisITSDZ(), calibInterpolation->GetHisITSTRDDZ(),  calibInterpolation->GetHisITSTOFDZ()};
  Int_t currentTrack=0;
  for (Int_t ihis=0; ihis<6; ihis++){    
    for (Int_t iesd=0; iesd<nesd; iesd+=downscale){
      TFile *esdFile = TFile::Open(esdArray->At(iesd)->GetName(),"read");
      if (!esdFile) continue;
      TTree *tree = (TTree*)esdFile->Get("delta"); 
      if (!tree) continue;
      AliSysInfo::AddStamp("xxx",ihis,iesd,currentTrack);
      tree->SetBranchAddress("vecR.",&vecR);
      tree->SetBranchAddress("vecPhi.",&vecPhi);
      tree->SetBranchAddress("vecZ.",&vecZ);
      tree->SetBranchAddress("track.",&param);
      tree->SetBranchAddress(branches[ihis],&vecDelta);
      Int_t ntracks=tree->GetEntries();
      for (Int_t itrack=0; itrack<ntracks; itrack++){
	tree->GetEntry(itrack);
	currentTrack++;
	for (Int_t ipoint=0; ipoint<knPoints; ipoint++){
	  if ((*vecR)[ipoint]<=0) continue;
	  Double_t sector=9.*(*vecPhi)[ipoint]/TMath::Pi();
	  if (sector<0) sector+=18;
	  Double_t deltaPhi=(*vecPhi)[ipoint]-TMath::Pi()*(sector+0.5)/9.;
	  Double_t localX = TMath::Cos(deltaPhi)*(*vecR)[ipoint];
	  Double_t xxx[5]={(*vecDelta)[ipoint], sector, localX,   (*vecZ)[ipoint]/(*vecR)[ipoint], param->GetParameter()[4] };
	  if (xxx[0]==0) continue;
	  hisToFill[ihis]->Fill(xxx);	  
	}	
      }
    }
  }
  TFile * fout = TFile::Open("ResidualHistograms.root","recreate");
  calibInterpolation->GetHisITSDRPhi()->Write();
  calibInterpolation->GetHisITSTRDDRPhi()->Write();
  calibInterpolation->GetHisITSTOFDRPhi()->Write();
  calibInterpolation->GetHisITSDZ()->Write();
  calibInterpolation->GetHisITSTRDDZ()->Write();
  calibInterpolation->GetHisITSTOFDZ()->Write();
  //calibInterpolation->Write();
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
  }
  return tree;
}
