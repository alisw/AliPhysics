/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Author: Daniel Mühlheim													*
 * Version 1.0																*
 *																			*
 * Permission to use, copy, modify and distribute this software and its	 	*
 * documentation strictly for non-commercial purposes is hereby granted	 	*
 * without fee, provided that the above copyright notice appears in all	 	*
 * copies and that both the copyright notice and this permission notice	 	*
 * appear in the supporting documentation. The authors make no claims		*
 * about the suitability of this software for any purpose. It is			*
 * provided "as is" without express or implied warranty.					*
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Basic Track Matching Class
//---------------------------------------------
////////////////////////////////////////////////


#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliCaloTrackMatcher.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliTrackerBase.h"
#include "AliV0ReaderV1.h"

#include "TChain.h"

#include <vector>
#include <map>
#include <utility>

class iostream;

using namespace std;


ClassImp(AliCaloTrackMatcher)

//________________________________________________________________________
AliCaloTrackMatcher::AliCaloTrackMatcher(const char *name, Int_t clusterType) : AliAnalysisTaskSE(name),
  fClusterType(clusterType),
  fV0ReaderName(""),
  fMatchingWindow(200),
  fMatchingResidual(0.2),
  fRunNumber(-1),
  fGeomEMCAL(NULL),
  fGeomPHOS(NULL),
  fMapTrackToCluster(),
  fMapClusterToTrack(),
  fNEntries(1),
  fVectorDeltaEtaDeltaPhi(0),
  fMap_TrID_ClID_ToIndex()

{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliCaloTrackMatcher::~AliCaloTrackMatcher()
{
    // default deconstructor
    fMapTrackToCluster.clear();
    fMapClusterToTrack.clear();
    fVectorDeltaEtaDeltaPhi.clear();
    fMap_TrID_ClID_ToIndex.clear();
}

//________________________________________________________________________
void AliCaloTrackMatcher::Terminate(Option_t *)
{
  fMapTrackToCluster.clear();
  fMapClusterToTrack.clear();
  fVectorDeltaEtaDeltaPhi.clear();
  fMap_TrID_ClID_ToIndex.clear();
}

//________________________________________________________________________
void AliCaloTrackMatcher::UserCreateOutputObjects()
{
    // Create User Output Objects
}

//________________________________________________________________________
void AliCaloTrackMatcher::Initialize(Int_t runNumber)
{
  // Initialize function to be called once before analysis
  fMapTrackToCluster.clear();
  fMapClusterToTrack.clear();
  fNEntries = 1;
  fVectorDeltaEtaDeltaPhi.clear();
  fMap_TrID_ClID_ToIndex.clear();

  if(fRunNumber == -1 || fRunNumber != runNumber){
    if(fClusterType == 1){
      fGeomEMCAL = AliEMCALGeometry::GetInstanceFromRunNumber(runNumber);
      if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
    }else if(fClusterType == 2){
      fGeomPHOS = AliPHOSGeometry::GetInstance();
      if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
    }
    fRunNumber = runNumber;
  }
}

//________________________________________________________________________
void AliCaloTrackMatcher::UserExec(Option_t *){

  // process only for EMCal (1) or PHOS (2) clusters
  if(fClusterType == 1 || fClusterType == 2){
    Initialize(fInputEvent->GetRunNumber());
    ProcessEvent(fInputEvent);
  }

//  DEBUGGING SECTION
//  cout << "******************************" << endl;
//  cout << "******************************" << endl;
//  cout << "NEW EVENT !" << endl;
//  cout << "vector etaphi:" << endl;
//  cout << fVectorDeltaEtaDeltaPhi.size() << endl;
//  cout << "multimap" << endl;
//  mapT::iterator iter = fMap_TrID_ClID_ToIndex.begin();
//  for (iter = fMap_TrID_ClID_ToIndex.begin(); iter != fMap_TrID_ClID_ToIndex.end(); ++iter){
//       Float_t dEta, dPhi = 0;
//       GetTrackClusterMatchingResidual(iter->first.first,iter->first.second,dEta,dPhi);
//       cout << "  [" << iter->first.first << "/" << iter->first.second << ", " << iter->second << "] - (" << dEta << "/" << dPhi << ")" << endl;
//   }
//  cout << "mapTrackToCluster" << endl;
//  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(fInputEvent);
//  for (Int_t itr=0;itr<esdev->GetNumberOfTracks();itr++){
//    AliVTrack *inTrack = esdev->GetTrack(itr);
//    if(!inTrack) continue;
//    cout << itr << " - " << GetNMatchedClusterIDsForTrack(inTrack->GetID(),5,-5,5,-5) << "\t\t";
//  }
//  cout << endl;
//  multimap<Int_t,Int_t>::iterator it;
//  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it) cout << it->first << " => " << it->second << '\n';
//  cout << "mapClusterToTrack" << endl;
//  Int_t tempClus = it->second;
//  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it) cout << it->first << " => " << it->second << '\n';
//  vector<Int_t> tempTracks = GetMatchedTrackIDsForCluster(tempClus, 5, -5, 5, -5);
//  for(Int_t iJ=0; iJ<tempTracks.size();iJ++){
//    cout << tempClus << " - " << tempTracks.at(iJ) << endl;
//  }
}

//________________________________________________________________________
void AliCaloTrackMatcher::ProcessEvent(AliVEvent *event)
{
  Int_t nClus = event->GetNumberOfCaloClusters();
  Int_t nModules = 0;
  if(fClusterType == 1) nModules = fGeomEMCAL->GetNumberOfSuperModules();
  else if(fClusterType == 2) nModules = fGeomPHOS->GetNModules();

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }

  for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){
    AliExternalTrackParam *trackParam = 0;
    AliVTrack *inTrack = 0x0;
    if(esdev){
      inTrack = esdev->GetTrack(itr);
      if(!inTrack) continue;
      AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);

      const AliExternalTrackParam *in = esdt->GetInnerParam();
      if (!in){AliDebug(2, "Could not get InnerParam of Track, continue");continue;}
      trackParam = new AliExternalTrackParam(*in);
    } else if(aodev) {
      if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
        inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
      } else {
        for(Int_t ii=0;ii<aodev->GetNumberOfTracks();ii++) {
          inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(ii));
          if(inTrack){
            if(inTrack->GetID() == itr) {
              break;
            }
          }
        }
      }

      if(!inTrack) continue;
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);

      Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
      aodt->GetPxPyPz(pxpypz);
      aodt->GetXYZ(xyz);
      aodt->GetCovarianceXYZPxPyPz(cv);
      trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
    }

    if (!trackParam) {AliError("Could not get TrackParameters, continue");continue;}
    AliExternalTrackParam emcParam(*trackParam);
    Float_t eta, phi, pt;

    //propagate tracks to emc surfaces
    if(fClusterType == 1){
      if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
        delete trackParam;
        continue;
      }

      if( TMath::Abs(eta) > 0.75 ) {
        delete trackParam;
        continue;
      }
      // Save some time and memory in case of no DCal present
      if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
        delete trackParam;
        continue;
      }
    }else if(fClusterType == 2){
      if( !AliTrackerBase::PropagateTrackToBxByBz(&emcParam, 460., 0.139, 20, kTRUE, 0.8, -1)){
        delete trackParam;
        continue;
      }
//to do: implement of distance checks
    }

    Float_t dEta=-999, dPhi=-999;
    Float_t clsPos[3] = {0.,0.,0.};
    Double_t exPos[3] = {0.,0.,0.};
    if (!emcParam.GetXYZ(exPos)) continue;

//cout << inTrack->GetID() << " - " << trackParam << endl;
//cout << "eta/phi: " << eta << ", " << phi << endl;
//cout << "nClus: " << nClus << endl;
    for(Int_t iclus=0;iclus < nClus;iclus++){
      AliVCluster* cluster = event->GetCaloCluster(iclus);
      if (!cluster) continue;
//cout << "-------------------------LOOPING: " << iclus << ", " << cluster->GetID() << endl;
      cluster->GetPosition(clsPos);
      Double_t dR = TMath::Sqrt(TMath::Power(exPos[0]-clsPos[0],2)+TMath::Power(exPos[1]-clsPos[1],2)+TMath::Power(exPos[2]-clsPos[2],2));
//cout << "dR: " << dR << endl;
      if (dR > fMatchingWindow) continue;
      Double_t clusterR = TMath::Sqrt( clsPos[0]*clsPos[0] + clsPos[1]*clsPos[1] );

      AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
      if(fClusterType == 1){
        if (!cluster->IsEMCAL()) continue;
        if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, cluster, 0.139, 5., dEta, dPhi)) continue;
      }else if(fClusterType == 2){
        if (!cluster->IsPHOS()) continue;
        if(!AliTrackerBase::PropagateTrackToBxByBz(&trackParamTmp, clusterR, 0.139, 5., kTRUE, 0.8, -1)) continue;
        Double_t trkPos[3] = {0,0,0};
        trackParamTmp.GetXYZ(trkPos);
        TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
        TVector3 clsPosVec(clsPos);
        dPhi = clsPosVec.DeltaPhi(trkPosVec);
        dEta = clsPosVec.Eta()-trkPosVec.Eta();
      }


      Float_t dR2 = dPhi*dPhi + dEta*dEta;
//cout << dEta << " - " << dPhi << " - " << dR2 << endl;
      if(dR2 > fMatchingResidual) continue;
//cout << "MATCHED!!!!!!!" << endl;

      fMapTrackToCluster.insert(make_pair(inTrack->GetID(),cluster->GetID()));
      fMapClusterToTrack.insert(make_pair(cluster->GetID(),inTrack->GetID()));
      fVectorDeltaEtaDeltaPhi.push_back(make_pair(dEta,dPhi));
      fMap_TrID_ClID_ToIndex[make_pair(inTrack->GetID(),cluster->GetID())] = fNEntries++;
      if( (Int_t)fVectorDeltaEtaDeltaPhi.size() != (fNEntries-1)) AliFatal("Fatal error in AliCaloTrackMatcher, vector and map are not in sync!");
    }

    delete trackParam;
  }

  return;
}

//________________________________________________________________________
Bool_t AliCaloTrackMatcher::GetTrackClusterMatchingResidual(Int_t trackID, Int_t clusterID, Float_t &dEta, Float_t &dPhi)
{
  Int_t position = fMap_TrID_ClID_ToIndex[make_pair(trackID,clusterID)];
  if(position == 0) return kFALSE;

  pairFloat tempEtaPhi = fVectorDeltaEtaDeltaPhi.at(position-1);
  dEta = tempEtaPhi.first;
  dPhi = tempEtaPhi.second;
  return kTRUE;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedTrackIDsForCluster(Int_t clusterID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg)
{
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(it->second,it->first,tempDEta,tempDPhi)){
        if( (dEtaNeg < tempDEta) && (tempDEta < dEtaPos) && (dPhiNeg < tempDPhi) && (tempDPhi < dPhiPos) ) matched++;
      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForTrack(Int_t trackID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg)
{
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == trackID){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(it->first,it->second,tempDEta,tempDPhi)){
        if( (dEtaNeg < tempDEta) && (tempDEta < dEtaPos) && (dPhiNeg < tempDPhi) && (tempDPhi < dPhiPos) ) matched++;
      }
    }
  }

  return matched;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedTrackIDsForCluster(Int_t clusterID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg)
{
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(it->second,it->first,tempDEta,tempDPhi)){
        if( (dEtaNeg < tempDEta) && (tempDEta < dEtaPos) && (dPhiNeg < tempDPhi) && (tempDPhi < dPhiPos) ) tempMatchedTracks.push_back(it->second);
      }
    }
  }

  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForTrack(Int_t trackID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg)
{
  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == trackID){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(it->first,it->second,tempDEta,tempDPhi)){
        if( (dEtaNeg < tempDEta) && (tempDEta < dEtaPos) && (dPhiNeg < tempDPhi) && (tempDPhi < dPhiPos) ) tempMatchedClusters.push_back(it->second);
      }
    }
  }

  return tempMatchedClusters;
}
