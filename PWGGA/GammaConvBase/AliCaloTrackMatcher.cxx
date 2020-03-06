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

#include "TAxis.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <utility>

class iostream;

using namespace std;


ClassImp(AliCaloTrackMatcher)

//________________________________________________________________________
AliCaloTrackMatcher::AliCaloTrackMatcher(const char *name, Int_t clusterType, Int_t runningMode) : AliAnalysisTaskSE(name),
  fClusterType(clusterType),
  fRunningMode(runningMode),
  fV0ReaderName(""),
  fCorrTaskSetting(""),
  fAnalysisTrainMode("Grid"),
  fMatchingWindow(200),
  fMatchingResidual(0.2),
  fRunNumber(-1),
  fGeomEMCAL(NULL),
  fGeomPHOS(NULL),
  fMapTrackToCluster(),
  fMapClusterToTrack(),
  fNEntries(1),
  fVectorDeltaEtaDeltaPhi(0),
  fMap_TrID_ClID_ToIndex(),
  fSecMapTrackToCluster(),
  fSecMapClusterToTrack(),
  fSecNEntries(1),
  fSecVectorDeltaEtaDeltaPhi(0),
  fSecMap_TrID_ClID_ToIndex(),
  fSecMap_TrID_ClID_AlreadyTried(),
  fListHistos(NULL),
  fHistControlMatches(NULL),
  fSecHistControlMatches(NULL),
  fDoLightOutput(kFALSE)
{
    // Default constructor
    DefineInput(0, TChain::Class());

    if( !(fRunningMode == 0 || fRunningMode == 1 || fRunningMode == 2 || fRunningMode == 5 || fRunningMode == 6) )
      AliFatal(Form("AliCaloTrackMatcher: Running mode not defined: '%d' ! Please set a proper mode in the AddTask, returning...",fRunningMode));
}

//________________________________________________________________________
AliCaloTrackMatcher::~AliCaloTrackMatcher(){
    // default deconstructor
    fMapTrackToCluster.clear();
    fMapClusterToTrack.clear();
    fVectorDeltaEtaDeltaPhi.clear();
    fMap_TrID_ClID_ToIndex.clear();

    fSecMapTrackToCluster.clear();
    fSecMapClusterToTrack.clear();
    fSecVectorDeltaEtaDeltaPhi.clear();
    fSecMap_TrID_ClID_ToIndex.clear();
    fSecMap_TrID_ClID_AlreadyTried.clear();

    if(fHistControlMatches) delete fHistControlMatches;
    if(fSecHistControlMatches) delete fSecHistControlMatches;
    if(fAnalysisTrainMode.EqualTo("Grid")){
        if(fListHistos != NULL){
            delete fListHistos;
        }
    }
}

//________________________________________________________________________
void AliCaloTrackMatcher::Terminate(Option_t *){
  fMapTrackToCluster.clear();
  fMapClusterToTrack.clear();
  fVectorDeltaEtaDeltaPhi.clear();
  fMap_TrID_ClID_ToIndex.clear();

  fSecMapTrackToCluster.clear();
  fSecMapClusterToTrack.clear();
  fSecVectorDeltaEtaDeltaPhi.clear();
  fSecMap_TrID_ClID_ToIndex.clear();
  fSecMap_TrID_ClID_AlreadyTried.clear();
}

//________________________________________________________________________
void AliCaloTrackMatcher::UserCreateOutputObjects(){
  if(fListHistos != NULL){
    delete fListHistos;
    fListHistos = NULL;
  }

  if(fDoLightOutput) return;

  if(fListHistos == NULL){
    fListHistos = new TList();
    fListHistos->SetOwner(kTRUE);
    if(!fCorrTaskSetting.CompareTo(""))
      fListHistos->SetName(Form("CaloTrackMatcher_%i_%i",fClusterType,fRunningMode));
    else
      fListHistos->SetName(Form("CaloTrackMatcher_%i_%i_%s",fClusterType,fRunningMode,fCorrTaskSetting.Data()));
  }

  // Create User Output Objects
  fHistControlMatches  = new TH2F(Form("ControlMatches_%i_%i",fClusterType,fRunningMode),Form("ControlMatches_%i_%i",fClusterType,fRunningMode),7,-0.5,6.5,50,0.05,200.0);
  SetLogBinningYTH2(fHistControlMatches);
  fHistControlMatches->GetYaxis()->SetTitle("track pT (GeV/c)");
  fHistControlMatches->GetXaxis()->SetBinLabel(1,"nTr in");
  fHistControlMatches->GetXaxis()->SetBinLabel(2,"no inner Tr params || track not in right direction");
  fHistControlMatches->GetXaxis()->SetBinLabel(3,"failed 1st pro-step");
  fHistControlMatches->GetXaxis()->SetBinLabel(4,"Tr not in EMCal acc");
  fHistControlMatches->GetXaxis()->SetBinLabel(5,"failed 2nd pro-step");
  fHistControlMatches->GetXaxis()->SetBinLabel(6,"w/o match to cluster");
  fHistControlMatches->GetXaxis()->SetBinLabel(7,"nTr out, w/ match");
  fListHistos->Add(fHistControlMatches);

  fSecHistControlMatches  = new TH2F(Form("ControlSecMatches_%i_%i",fClusterType,fRunningMode),Form("ControlSecMatches_%i_%i",fClusterType,fRunningMode),7,-0.5,6.5,50,0.05,200.0);
  SetLogBinningYTH2(fSecHistControlMatches);
  fSecHistControlMatches->GetYaxis()->SetTitle("track pT (GeV/c)");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(1,"nTr in");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(2,"no inner Tr params || track not in right direction");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(3,"failed 1st pro-step");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(4,"Tr not in EMCal acc");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(5,"failed 2nd pro-step");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(6,"w/o match to cluster");
  fSecHistControlMatches->GetXaxis()->SetBinLabel(7,"nTr out, w/ match");
  fListHistos->Add(fSecHistControlMatches);
}

//________________________________________________________________________
void AliCaloTrackMatcher::Initialize(Int_t runNumber){
  // Initialize function to be called once before analysis
  fMapTrackToCluster.clear();
  fMapClusterToTrack.clear();
  fNEntries = 1;
  fVectorDeltaEtaDeltaPhi.clear();
  fMap_TrID_ClID_ToIndex.clear();

  fSecMapTrackToCluster.clear();
  fSecMapClusterToTrack.clear();
  fSecNEntries = 1;
  fSecVectorDeltaEtaDeltaPhi.clear();
  fSecMap_TrID_ClID_ToIndex.clear();
  fSecMap_TrID_ClID_AlreadyTried.clear();

  if(fRunNumber == -1 || fRunNumber != runNumber){
    if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
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

  // main method of AliCaloTrackMatcher, first initialize and then process event

  //DebugV0Matching();

  // do processing only for EMCal (1), DCal (3) or PHOS (2) clusters, otherwise do nothing
  if(fClusterType == 1 || fClusterType == 2 || fClusterType == 3 || fClusterType == 4){
    Initialize(fInputEvent->GetRunNumber());
    ProcessEvent(fInputEvent);
  }

  //DebugMatching();

  return;
}

//________________________________________________________________________
void AliCaloTrackMatcher::ProcessEvent(AliVEvent *event){
  Int_t nClus = 0;
  TClonesArray * arrClusters = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nClus = event->GetNumberOfCaloClusters();
  } else {
    arrClusters = dynamic_cast<TClonesArray*>(event->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(arrClusters){
      nClus = arrClusters->GetEntries();
    }else{
      AliError(Form("Could not find %sClustersBranch despite correction framework being used!",fCorrTaskSetting.Data()));
      return;
    }
  }
  Int_t nModules = 0;
  if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4) nModules = fGeomEMCAL->GetNumberOfSuperModules();
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
  static AliESDtrackCuts *EsdTrackCuts = 0x0;
  static int prevRun = -1;
  // Using standard function for setting Cuts
  if (esdev && (fRunningMode == 0 || fRunningMode == 6)){
    Int_t runNumber = fInputEvent->GetRunNumber();
    if (prevRun!=runNumber) {
      delete EsdTrackCuts;
      EsdTrackCuts = 0;
      prevRun = runNumber;
    }
    if (!EsdTrackCuts) {
      // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
      if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
        if(fRunningMode == 6) EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
        else EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
        // else if run2 data use 2015 PbPb cuts
      } else if (runNumber>=209122){
        //EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
        // hard coded track cuts for the moment, because AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() gives spams warnings
        EsdTrackCuts = new AliESDtrackCuts();
        // TPC; clusterCut = 1, cutAcceptanceEdges = kTRUE, removeDistortedRegions = kFALSE
        EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
        EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
        EsdTrackCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);  // only dead zone and not clusters per length
        //EsdTrackCuts->AliESDtrackCuts::SetCutOutDistortedRegionsTPC(kTRUE);
        EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
        EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
        EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
        // ITS; selPrimaries = 1
        EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
        EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                                                AliESDtrackCuts::kAny);
        if(fRunningMode != 6){
          EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
          EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
        }
        EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
        EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
        EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
        EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);
        // else use 2011 version of track cuts
      }else{
        if(fRunningMode == 6) EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
        else EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      }
      EsdTrackCuts->SetMaxDCAToVertexZ(2);
    }
  }

  for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){
    AliExternalTrackParam *trackParam = 0;
    AliVTrack *inTrack = 0x0;
    if(esdev){
      inTrack = esdev->GetTrack(itr);
      if(!inTrack) continue;
      FillfHistControlMatches(0.,inTrack->Pt());
      AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);

      if(fRunningMode == 0){
        if(TMath::Abs(inTrack->Eta())>0.8 && (fClusterType == 1 || fClusterType == 3  || fClusterType == 4)){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(TMath::Abs(inTrack->Eta())>0.3 &&  fClusterType == 2 ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(inTrack->Pt()<0.5){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(!EsdTrackCuts->AcceptTrack(esdt)) {FillfHistControlMatches(1.,inTrack->Pt()); continue;}
      }else if(fRunningMode == 5 || fRunningMode == 6){
        if(TMath::Abs(inTrack->Eta())>0.8 && (fClusterType == 1 || fClusterType == 3  || fClusterType == 4)){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(TMath::Abs(inTrack->Eta())>0.3 &&  fClusterType == 2 ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(inTrack->Pt()<0.3){FillfHistControlMatches(1.,inTrack->Pt()); continue;}

        if(fRunningMode == 6){
          if(!EsdTrackCuts->AcceptTrack(esdt)) {FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        }
      }

      const AliExternalTrackParam *in = esdt->GetInnerParam();
      if (!in){AliDebug(2, "Could not get InnerParam of Track, continue"); FillfHistControlMatches(1.,inTrack->Pt()); continue;}
      trackParam = new AliExternalTrackParam(*in);
    } else if(aodev) {
      inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
      if(!inTrack) continue;
      FillfHistControlMatches(0.,inTrack->Pt());
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);

      if(fRunningMode == 0){
        if(inTrack->GetID()<0){FillfHistControlMatches(1.,inTrack->Pt()); continue;} // Avoid double counting of tracks
        if(TMath::Abs(inTrack->Eta())>0.8 && (fClusterType == 1 || fClusterType == 3  || fClusterType == 4)){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(TMath::Abs(inTrack->Eta())>0.3 &&  fClusterType == 2 ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(inTrack->Pt()<0.5){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(!aodt->IsHybridGlobalConstrainedGlobal()){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
      }else if(fRunningMode == 5 || fRunningMode == 6){
        if(inTrack->GetID()<0){FillfHistControlMatches(1.,inTrack->Pt()); continue;} // Avoid double counting of tracks
        if(TMath::Abs(inTrack->Eta())>0.8 && (fClusterType == 1 || fClusterType == 3  || fClusterType == 4)){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(TMath::Abs(inTrack->Eta())>0.3 &&  fClusterType == 2 ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        if(inTrack->Pt()<0.3){FillfHistControlMatches(1.,inTrack->Pt()); continue;}

        if(fRunningMode == 6){
          if(!aodt->IsHybridGlobalConstrainedGlobal()){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        }
      }

      Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
      aodt->GetPxPyPz(pxpypz);
      aodt->GetXYZ(xyz);
      aodt->GetCovarianceXYZPxPyPz(cv);

      // check for EMC tracks already propagated tracks are out of bounds
      if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        if( TMath::Abs(aodt->GetTrackEtaOnEMCal()) > 0.75 ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}

        // conditions for run1
        if( fClusterType == 1 && nModules < 13 && ( aodt->GetTrackPhiOnEMCal() < 70*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 190*TMath::DegToRad())){FillfHistControlMatches(1.,inTrack->Pt()); continue;}

        // conditions for run2
        if( nModules > 12 ){
          if( fClusterType == 3 && ( aodt->GetTrackPhiOnEMCal() < 250*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 340*TMath::DegToRad()) ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
          if( fClusterType == 1 && ( aodt->GetTrackPhiOnEMCal() < 70*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 190*TMath::DegToRad()) ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
          if( fClusterType == 4 && ( aodt->GetTrackPhiOnEMCal() < 70*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 190*TMath::DegToRad()) && ( aodt->GetTrackPhiOnEMCal() < 250*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 340*TMath::DegToRad()) ){FillfHistControlMatches(1.,inTrack->Pt()); continue;}
        }
      }
      trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
    }

    if (!trackParam) {AliError("Could not get TrackParameters, continue"); FillfHistControlMatches(1.,inTrack->Pt()); continue;}
    AliExternalTrackParam emcParam(*trackParam);
    Float_t eta, phi, pt;

    //propagate tracks to emc surfaces
    if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
      if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
        delete trackParam;
        FillfHistControlMatches(2.,inTrack->Pt());
        continue;
      }
      if( TMath::Abs(eta) > 0.75 ) {
        delete trackParam;
        FillfHistControlMatches(3.,inTrack->Pt());
        continue;
      }
      // Save some time and memory in case of no DCal present
      if( fClusterType == 1 && nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
        delete trackParam;
        FillfHistControlMatches(3.,inTrack->Pt());
        continue;
      }
      // Save some time and memory in case of run2
      if( nModules > 12 ){
        if (fClusterType == 3 && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
          delete trackParam;
          FillfHistControlMatches(3.,inTrack->Pt());
          continue;
        }
        if( fClusterType == 1 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
          delete trackParam;
          FillfHistControlMatches(3.,inTrack->Pt());
          continue;
        }
        if( fClusterType == 4 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())
                              && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
          delete trackParam;
          FillfHistControlMatches(3.,inTrack->Pt());
          continue;
        }
      }
    }else if(fClusterType == 2){
      if( !AliTrackerBase::PropagateTrackToBxByBz(&emcParam, 460., 0.139, 20, kTRUE, 0.8, -1)){
        delete trackParam;
        FillfHistControlMatches(2.,inTrack->Pt());
        continue;
      }
      Double_t trkPos[3] = {0,0,0};
      if(emcParam.GetXYZ(trkPos)){
        TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
        if (TMath::Abs(trkPosVec.Eta()) > 0.25 ){
          delete trackParam;
          FillfHistControlMatches(3.,inTrack->Pt());
          continue;
        }
        if (trkPosVec.Phi() < 230*TMath::DegToRad() || trkPosVec.Phi() > 350*TMath::DegToRad() ){
          delete trackParam;
          FillfHistControlMatches(3.,inTrack->Pt());
          continue;
        }
      }
    }

    Float_t dEta=-999, dPhi=-999;
    Float_t clsPos[3] = {0.,0.,0.};
    Double_t exPos[3] = {0.,0.,0.};
    if (!emcParam.GetXYZ(exPos)){
      delete trackParam;
      FillfHistControlMatches(2.,inTrack->Pt());
      continue;
    }

    // cout << inTrack->GetID() << " - " << trackParam << endl;
    // cout << "eta/phi: " << eta << ", " << phi << endl;
    // cout << "nClus: " << nClus << endl;
    Int_t nClusterMatchesToTrack = 0;
    for(Int_t iclus=0;iclus < nClus;iclus++){
      AliVCluster* cluster = NULL;
      if(arrClusters){
        if(esdev){
          if(arrClusters)
            cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClusters->At(iclus));
        } else if(aodev){
          if(arrClusters)
            cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClusters->At(iclus));
        }
      }
      else
        cluster = event->GetCaloCluster(iclus);
      if (!cluster){
        if(arrClusters) delete cluster;
        continue;
      }
      // cout << "-------------------------LOOPING: " << iclus << ", " << cluster->GetID() << endl;
      cluster->GetPosition(clsPos);
      Double_t dR = TMath::Sqrt(TMath::Power(exPos[0]-clsPos[0],2)+TMath::Power(exPos[1]-clsPos[1],2)+TMath::Power(exPos[2]-clsPos[2],2));
      //cout << "dR: " << dR << endl;
      if (dR > fMatchingWindow){
        if(arrClusters) delete cluster;
        continue;
      }
      Double_t clusterR = TMath::Sqrt( clsPos[0]*clsPos[0] + clsPos[1]*clsPos[1] );
      AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
      if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        if (!cluster->IsEMCAL()){
          if(arrClusters) delete cluster;
          continue;
        }
        if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, cluster, 0.139, 5., dEta, dPhi)){
          FillfHistControlMatches(4.,inTrack->Pt());
          if(arrClusters) delete cluster;
          continue;
        }
      }else if(fClusterType == 2){
        if (!cluster->IsPHOS()){
          if(arrClusters) delete cluster;
          continue;
        }
        if(!AliTrackerBase::PropagateTrackToBxByBz(&trackParamTmp, clusterR, 0.139, 5., kTRUE, 0.8, -1)){
          FillfHistControlMatches(4.,inTrack->Pt());
          if(arrClusters) delete cluster;
          continue;
        }
        Double_t trkPos[3] = {0,0,0};
        trackParamTmp.GetXYZ(trkPos);
        TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
        TVector3 clsPosVec(clsPos);
        dPhi = clsPosVec.DeltaPhi(trkPosVec);
        dEta = clsPosVec.Eta()-trkPosVec.Eta();
      }

      Float_t dR2 = dPhi*dPhi + dEta*dEta;

      //cout << dEta << " - " << dPhi << " - " << dR2 << endl;
      if(dR2 > fMatchingResidual){
        if(arrClusters) delete cluster;
        continue;
      }
      nClusterMatchesToTrack++;
      if(aodev){
        fMapTrackToCluster.insert(make_pair(itr,cluster->GetID()));
        fMapClusterToTrack.insert(make_pair(cluster->GetID(),itr));
      }else{
        fMapTrackToCluster.insert(make_pair(inTrack->GetID(),cluster->GetID()));
        fMapClusterToTrack.insert(make_pair(cluster->GetID(),inTrack->GetID()));
      }
      fVectorDeltaEtaDeltaPhi.push_back(make_pair(dEta,dPhi));
      fMap_TrID_ClID_ToIndex[make_pair(inTrack->GetID(),cluster->GetID())] = fNEntries++;
      if( (Int_t)fVectorDeltaEtaDeltaPhi.size() != (fNEntries-1)) AliFatal("Fatal error in AliCaloTrackMatcher, vector and map are not in sync!");
      if(arrClusters) delete cluster;
    }
    if(nClusterMatchesToTrack == 0) FillfHistControlMatches(5.,inTrack->Pt());
    else FillfHistControlMatches(6.,inTrack->Pt());
    delete trackParam;
  }

  return;
}

//________________________________________________________________________
Bool_t AliCaloTrackMatcher::PropagateV0TrackToClusterAndGetMatchingResidual(AliVTrack* inSecTrack, AliVCluster* cluster, AliVEvent* event, Float_t &dEta, Float_t &dPhi){

  //if V0-track to cluster match is already available return stored residuals
//   cout << "matching sec track: " << inSecTrack->GetID() << "\t to cluster" << cluster->GetID() << endl;
  if(GetSecTrackClusterMatchingResidual(inSecTrack->GetID(),cluster->GetID(), dEta, dPhi)){
//     cout << "RESIDUAL ALREADY AVAILABLE! - " << dEta << "/" << dPhi << endl;
    return kTRUE;
  }

  if(IsSecTrackClusterAlreadyTried(inSecTrack->GetID(),cluster->GetID())){
  //cout << "PROPAGATION ALREADY FAILED! - " << inSecTrack->GetID() << "/" << cluster->GetID() << endl;
    return kFALSE;
  }

  //cout << "running matching! - " << inSecTrack->GetID() << "/" << cluster->GetID() << endl;
  //if match has not yet been computed, go on:
  Int_t nModules = 0;
  if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4) nModules = fGeomEMCAL->GetNumberOfSuperModules();
  else if(fClusterType == 2) nModules = fGeomPHOS->GetNModules();

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

  if(!cluster->IsEMCAL() && !cluster->IsPHOS()){AliError("Cluster is neither EMCAL nor PHOS, returning");return kFALSE;}

  Float_t clusterPosition[3] = {0,0,0};
  cluster->GetPosition(clusterPosition);
  Double_t clusterR = TMath::Sqrt( clusterPosition[0]*clusterPosition[0] + clusterPosition[1]*clusterPosition[1] );

  if(!inSecTrack) return kFALSE;
  FillfSecHistControlMatches(0.,inSecTrack->Pt());

  if (inSecTrack->Pt() < 0.3 ) {
    FillfSecHistControlMatches(1.,inSecTrack->Pt());
    fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
    return kFALSE;
  }

  AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inSecTrack);
  AliAODTrack *aodt = 0;
  if (!esdt) {
    aodt = dynamic_cast<AliAODTrack*>(inSecTrack);
    if (!aodt){
      AliError("Track is neither ESD nor AOD, continue");
      fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
      return kFALSE;
    }
  }

  AliExternalTrackParam *trackParam = 0;
  if (esdt) {
    const AliExternalTrackParam *in = esdt->GetInnerParam();
    if (!in){
      AliDebug(2, "Could not get InnerParam of Track, continue");
      FillfSecHistControlMatches(1.,inSecTrack->Pt());
      fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
      return kFALSE;
    }
    trackParam = new AliExternalTrackParam(*in);
  } else {
    // check if tracks should be propagated at all
    if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
      if (TMath::Abs(aodt->GetTrackEtaOnEMCal()) > 0.8){
        FillfSecHistControlMatches(1.,inSecTrack->Pt());
        fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
        return kFALSE;
      }
      if( nModules < 13 ){
        if (( aodt->GetTrackPhiOnEMCal() < 60*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 200*TMath::DegToRad())){
          FillfSecHistControlMatches(1.,inSecTrack->Pt());
          fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
          return kFALSE;
        }
      } else if( nModules > 12 ){
        if (fClusterType == 3 && ( aodt->GetTrackPhiOnEMCal() < 250*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 340*TMath::DegToRad())){
          FillfSecHistControlMatches(1.,inSecTrack->Pt());
          fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
          return kFALSE;
        }
        if( fClusterType == 1 && ( aodt->GetTrackPhiOnEMCal() < 60*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 200*TMath::DegToRad())){
          FillfSecHistControlMatches(1.,inSecTrack->Pt());
          fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
          return kFALSE;
        }
        if( fClusterType == 4 && ( aodt->GetTrackPhiOnEMCal() < 60*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 200*TMath::DegToRad())
                              && ( aodt->GetTrackPhiOnEMCal() < 250*TMath::DegToRad() || aodt->GetTrackPhiOnEMCal() > 340*TMath::DegToRad()) ){
          FillfSecHistControlMatches(1.,inSecTrack->Pt());
          fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
          return kFALSE;
        }
      }
    } else {
      if ( aodt->Phi() < 230*TMath::DegToRad() || aodt->Phi() > 350*TMath::DegToRad()){
        FillfSecHistControlMatches(1.,inSecTrack->Pt());
        fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
        return kFALSE;
      }
      if (TMath::Abs(aodt->Eta()) > 0.3 ){
        FillfSecHistControlMatches(1.,inSecTrack->Pt());
        fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
        return kFALSE;
      }
    }

    Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
    aodt->GetPxPyPz(pxpypz);
    aodt->GetXYZ(xyz);
    aodt->GetCovarianceXYZPxPyPz(cv);
    trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
  }
  if(!trackParam){
    AliError("Could not get TrackParameters, continue");
    FillfSecHistControlMatches(1.,inSecTrack->Pt());
    fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
    return kFALSE;
  }

  Bool_t propagated = kFALSE;
  AliExternalTrackParam emcParam(*trackParam);
  Float_t dPhiTemp = 0;
  Float_t dEtaTemp = 0;

  if(cluster->IsEMCAL()){
    Float_t eta = 0;Float_t phi = 0;Float_t pt = 0;
    propagated = AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 430, 0.000510999, 20, eta, phi, pt);
    if(propagated){
      if( TMath::Abs(eta) > 0.8 ) {
        delete trackParam;
        FillfSecHistControlMatches(3.,inSecTrack->Pt());
        fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
        return kFALSE;
      }
      // Save some time and memory in case of no DCal present
      if( nModules < 13 && ( phi < 60*TMath::DegToRad() || phi > 200*TMath::DegToRad())){
        delete trackParam;
        FillfSecHistControlMatches(3.,inSecTrack->Pt());
        fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
        return kFALSE;
      }

      propagated = AliEMCALRecoUtils::ExtrapolateTrackToCluster(&emcParam, cluster, 0.000510999, 5, dEtaTemp, dPhiTemp);
      if(!propagated){
        delete trackParam;
        FillfSecHistControlMatches(4.,inSecTrack->Pt());
        fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
        return kFALSE;
      }
    }else{
      delete trackParam;
      FillfSecHistControlMatches(2.,inSecTrack->Pt());
      fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
      return kFALSE;
    }

  }else if(cluster->IsPHOS()){
    propagated = AliTrackerBase::PropagateTrackToBxByBz(&emcParam, clusterR, 0.000510999, 20, kTRUE, 0.8, -1);
    if (propagated){
      Double_t trkPos[3] = {0,0,0};
      emcParam.GetXYZ(trkPos);
      TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
      TVector3 clsPosVec(clusterPosition);
      dPhiTemp = clsPosVec.DeltaPhi(trkPosVec);
      dEtaTemp = clsPosVec.Eta()-trkPosVec.Eta();
    }else{
      delete trackParam;
      FillfSecHistControlMatches(2.,inSecTrack->Pt());
      fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
      return kFALSE;}
  }

  if (propagated){
    Float_t dR2 = dPhiTemp*dPhiTemp + dEtaTemp*dEtaTemp;
    //cout << dEtaTemp << " - " << dPhiTemp << " - " << dR2 << endl;
    if(dR2 > fMatchingResidual){
      FillfSecHistControlMatches(5.,inSecTrack->Pt());
      fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
      //cout << "NO MATCH! - " << inSecTrack->GetID() << "/" << cluster->GetID() << endl;
      delete trackParam;
      return kFALSE;
    }
    //cout << "MATCHED!!!!!!!" << endl;

    if(aodev){
      //need to search for position in case of AOD
      Int_t TrackPos = -1;
      for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
        AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
        if(currTrack->GetID() == inSecTrack->GetID()){
          TrackPos = iTrack;
          break;
        }
      }
      if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: PropagateV0TrackToClusterAndGetMatchingResidual - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",inSecTrack->GetID()));
      fSecMapTrackToCluster.insert(make_pair(TrackPos,cluster->GetID()));
      fSecMapClusterToTrack.insert(make_pair(cluster->GetID(),TrackPos));
    }else{
      fSecMapTrackToCluster.insert(make_pair(inSecTrack->GetID(),cluster->GetID()));
      fSecMapClusterToTrack.insert(make_pair(cluster->GetID(),inSecTrack->GetID()));
    }
    fSecVectorDeltaEtaDeltaPhi.push_back(make_pair(dEtaTemp,dPhiTemp));
    fSecMap_TrID_ClID_ToIndex[make_pair(inSecTrack->GetID(),cluster->GetID())] = fSecNEntries++;
    if( (Int_t)fSecVectorDeltaEtaDeltaPhi.size() != (fSecNEntries-1)) AliFatal("Fatal error in AliCaloTrackMatcher, vector and map are not in sync!");

    FillfSecHistControlMatches(6.,inSecTrack->Pt());
    dEta = dEtaTemp;
    dPhi = dPhiTemp;
    delete trackParam;
    return kTRUE;
  }else AliFatal("Fatal error in AliCaloTrackMatcher, track is labeled as sucessfully propagated although this should be impossible!");

  fSecMap_TrID_ClID_AlreadyTried[make_pair(inSecTrack->GetID(),cluster->GetID())] = 1.;
  delete trackParam;
  return kFALSE;
}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
Bool_t AliCaloTrackMatcher::GetTrackClusterMatchingResidual(Int_t trackID, Int_t clusterID, Float_t &dEta, Float_t &dPhi){
  Int_t position = fMap_TrID_ClID_ToIndex[make_pair(trackID,clusterID)];
  if(position == 0) return kFALSE;

  pairFloat tempEtaPhi = fVectorDeltaEtaDeltaPhi.at(position-1);
  dEta = tempEtaPhi.first;
  dPhi = tempEtaPhi.second;
  return kTRUE;
}
//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) matched++;
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) matched++;
        }
      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedTrackIDsForCluster(AliVEvent *event, Int_t clusterID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )matched++;
      }
    }
  }
  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dR){
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) matched++;
      }
    }
  }
  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForTrack(AliVEvent *event, Int_t trackID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){

  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return matched;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) matched++;
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) matched++;
        }
      }
    }
  }
  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForTrack(AliVEvent *event, Int_t trackID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return matched;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )matched++;

      }
    }
  }
  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForTrack(AliVEvent *event, Int_t trackID, Float_t dR){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return matched;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) matched++;
      }
    }
  }
  return matched;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) tempMatchedTracks.push_back(it->second);
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) tempMatchedTracks.push_back(it->second);
        }
      }
    }
  }
  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedTrackIDsForCluster(AliVEvent *event, Int_t clusterID,  TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )tempMatchedTracks.push_back(it->second);

      }
    }
  }
  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedTrackIDsForCluster(AliVEvent *event, Int_t clusterID,  Float_t dR){
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) tempMatchedTracks.push_back(it->second);
      }
    }
  }
  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForTrack(AliVEvent *event, Int_t trackID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return tempMatchedClusters;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) tempMatchedClusters.push_back(it->second);
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) tempMatchedClusters.push_back(it->second);
        }
      }
    }
  }

  return tempMatchedClusters;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForTrack(AliVEvent *event, Int_t trackID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return tempMatchedClusters;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )tempMatchedClusters.push_back(it->second);
      }
    }
  }
  return tempMatchedClusters;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForTrack(AliVEvent *event, Int_t trackID, Float_t dR){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return tempMatchedClusters;
  for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) tempMatchedClusters.push_back(it->second);
      }
    }
  }
  return tempMatchedClusters;
}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
Bool_t AliCaloTrackMatcher::GetSecTrackClusterMatchingResidual(Int_t trackID, Int_t clusterID, Float_t &dEta, Float_t &dPhi){
  Int_t position = fSecMap_TrID_ClID_ToIndex[make_pair(trackID,clusterID)];
  if(position == 0) return kFALSE;

  pairFloat tempEtaPhi = fSecVectorDeltaEtaDeltaPhi.at(position-1);
  dEta = tempEtaPhi.first;
  dPhi = tempEtaPhi.second;
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloTrackMatcher::IsSecTrackClusterAlreadyTried(Int_t trackID, Int_t clusterID){
  Int_t position = fSecMap_TrID_ClID_AlreadyTried[make_pair(trackID,clusterID)];
  if(position == 0) return kFALSE;
  else return kTRUE;
}
//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedSecTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) matched++;
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) matched++;
        }
      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedSecTrackIDsForCluster(AliVEvent *event, Int_t clusterID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )matched++;
      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedSecTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dR){
  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) matched++;
      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForSecTrack(AliVEvent *event, Int_t trackID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return matched;
  for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) matched++;
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) matched++;
        }
      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForSecTrack(AliVEvent *event, Int_t trackID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return matched;
  for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )matched++;

      }
    }
  }

  return matched;
}

//________________________________________________________________________
Int_t AliCaloTrackMatcher::GetNMatchedClusterIDsForSecTrack(AliVEvent *event, Int_t trackID, Float_t dR){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  Int_t matched = 0;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return matched;
  for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) matched++;
      }
    }
  }

  return matched;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedSecTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) tempMatchedTracks.push_back(it->second);
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) tempMatchedTracks.push_back(it->second);
        }
      }
    }
  }

  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedSecTrackIDsForCluster(AliVEvent *event, Int_t clusterID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )tempMatchedTracks.push_back(it->second);
      }
    }
  }

  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedSecTrackIDsForCluster(AliVEvent *event, Int_t clusterID, Float_t dR){
  vector<Int_t> tempMatchedTracks;
  multimap<Int_t,Int_t>::iterator it;
  for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it){
    if(it->first == clusterID){
      Float_t tempDEta, tempDPhi;
      AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(it->second));
      if(!tempTrack) continue;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->first,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) tempMatchedTracks.push_back(it->second);
      }
    }
  }

  return tempMatchedTracks;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForSecTrack(AliVEvent *event, Int_t trackID, Float_t dEtaMax, Float_t dEtaMin, Float_t dPhiMax, Float_t dPhiMin){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return tempMatchedClusters;
  for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if(tempTrack->Charge()>0){
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin < tempDPhi) && (tempDPhi < dPhiMax) ) tempMatchedClusters.push_back(it->second);
        }else if(tempTrack->Charge()<0){
          dPhiMin*=-1;
          dPhiMax*=-1;
          if( (dEtaMin < tempDEta) && (tempDEta < dEtaMax) && (dPhiMin > tempDPhi) && (tempDPhi > dPhiMax) ) tempMatchedClusters.push_back(it->second);
        }
      }
    }
  }

  return tempMatchedClusters;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForSecTrack(AliVEvent *event, Int_t trackID, TF1* fFuncPtDepEta, TF1* fFuncPtDepPhi){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return tempMatchedClusters;
  for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        Bool_t match_dEta = kFALSE;
        Bool_t match_dPhi = kFALSE;
        if( TMath::Abs(tempDEta) < fFuncPtDepEta->Eval(tempTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(tempDPhi) < fFuncPtDepPhi->Eval(tempTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;

        if (match_dPhi && match_dEta )tempMatchedClusters.push_back(it->second);
      }
    }
  }

  return tempMatchedClusters;
}

//________________________________________________________________________
vector<Int_t> AliCaloTrackMatcher::GetMatchedClusterIDsForSecTrack(AliVEvent *event, Int_t trackID, Float_t dR){
  Int_t TrackPos = -1;
  if(event->IsA()==AliAODEvent::Class()){ // for AOD, we have to look for position of track in the event
    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++){
      AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(currTrack->GetID() == trackID){
        TrackPos = iTrack;
        break;
      }
    }
    if(TrackPos == -1) AliFatal(Form("AliCaloTrackMatcher: GetNMatchedClusterIDsForTrack - track (ID: '%i') cannot be retrieved from event, should be impossible as it has been used in maim task before!",trackID));
  }else TrackPos = trackID; // for ESD just take trackID

  vector<Int_t> tempMatchedClusters;
  multimap<Int_t,Int_t>::iterator it;
  AliVTrack* tempTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(TrackPos));
  if(!tempTrack) return tempMatchedClusters;
  for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it){
    if(it->first == TrackPos){
      Float_t tempDEta, tempDPhi;
      if(GetTrackClusterMatchingResidual(tempTrack->GetID(),it->second,tempDEta,tempDPhi)){
        if (TMath::Sqrt(tempDEta*tempDEta + tempDPhi*tempDPhi) < dR ) tempMatchedClusters.push_back(it->second);
      }
    }
  }

  return tempMatchedClusters;
}

//________________________________________________________________________
Float_t AliCaloTrackMatcher::SumTrackEtAroundCluster(AliVEvent* event, Int_t clusterID, Float_t dR){
  Float_t sumTrackEt = 0.;
  vector<Int_t> labelsMatched = GetMatchedTrackIDsForCluster(event, clusterID, dR);
  if((Int_t) labelsMatched.size()<1) return sumTrackEt;

  TLorentzVector vecTrack;
  for (UInt_t i = 0; i < labelsMatched.size(); i++){
    AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatched.at(i)));
    if(!currTrack) continue;
    vecTrack.SetPxPyPzE(currTrack->Px(),currTrack->Py(),currTrack->Pz(),currTrack->E());
    sumTrackEt += vecTrack.Et();
  }

  return sumTrackEt;
}

//________________________________________________________________________
void AliCaloTrackMatcher::SetLogBinningYTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetYaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}

//________________________________________________________________________
void AliCaloTrackMatcher::DebugV0Matching(){
  if(fSecVectorDeltaEtaDeltaPhi.size()>0){
    cout << "******************************" << endl;
    cout << "******************************" << endl;
    cout << "NEW EVENT !" << endl;
    cout << "vector etaphi:" << endl;
    cout << fSecVectorDeltaEtaDeltaPhi.size() << endl;
    cout << "multimap" << endl;
    mapT::iterator iter = fSecMap_TrID_ClID_ToIndex.begin();
    for (iter = fSecMap_TrID_ClID_ToIndex.begin(); iter != fSecMap_TrID_ClID_ToIndex.end(); ++iter){
      Float_t dEta, dPhi = 0;
      if(!GetSecTrackClusterMatchingResidual(iter->first.first,iter->first.second,dEta,dPhi)) continue;
      cout << "  [" << iter->first.first << "/" << iter->first.second << ", " << iter->second << "] - (" << dEta << "/" << dPhi << ")" << endl;
    }
    cout << "mapTrackToCluster" << endl;
    AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(fInputEvent);
    for (Int_t itr=0;itr<esdev->GetNumberOfTracks();itr++){
      AliVTrack *inTrack = esdev->GetTrack(itr);
      if(!inTrack) continue;
      TString tCharge;
      if(inTrack->Charge()>0) tCharge = "+";
      else if(inTrack->Charge()<0) tCharge = "-";
      cout << itr << " (" << tCharge << ") - " << GetNMatchedClusterIDsForSecTrack(fInputEvent,inTrack->GetID(),5,-5,0.2,-0.4) << "\t\t";
    }
    cout << endl;
    multimap<Int_t,Int_t>::iterator it;
    for (it=fSecMapTrackToCluster.begin(); it!=fSecMapTrackToCluster.end(); ++it) cout << it->first << " => " << it->second << '\n';
    cout << "mapClusterToTrack" << endl;
    Int_t tempClus = it->second;
    for (it=fSecMapClusterToTrack.begin(); it!=fSecMapClusterToTrack.end(); ++it) cout << it->first << " => " << it->second << '\n';
    vector<Int_t> tempTracks = GetMatchedSecTrackIDsForCluster(fInputEvent,tempClus, 5, -5, 0.2, -0.4);
    for(UInt_t iJ=0; iJ<tempTracks.size();iJ++){
      cout << tempClus << " - " << tempTracks.at(iJ) << endl;
    }
  }
  return;
}

//________________________________________________________________________
void AliCaloTrackMatcher::DebugMatching(){
  if(fVectorDeltaEtaDeltaPhi.size()>0){
    cout << "******************************" << endl;
    cout << "******************************" << endl;
    cout << "NEW EVENT !" << endl;
    cout << "vector etaphi:" << endl;
    cout << fVectorDeltaEtaDeltaPhi.size() << endl;
    cout << "multimap" << endl;
    mapT::iterator iter = fMap_TrID_ClID_ToIndex.begin();
    for (iter = fMap_TrID_ClID_ToIndex.begin(); iter != fMap_TrID_ClID_ToIndex.end(); ++iter){
      Float_t dEta, dPhi = 0;
      if(!GetTrackClusterMatchingResidual(iter->first.first,iter->first.second,dEta,dPhi)) continue;
      cout << "  [" << iter->first.first << "/" << iter->first.second << ", " << iter->second << "] - (" << dEta << "/" << dPhi << ")" << endl;
    }
    cout << "mapTrackToCluster" << endl;
    AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(fInputEvent);
    for (Int_t itr=0;itr<esdev->GetNumberOfTracks();itr++){
      AliVTrack *inTrack = esdev->GetTrack(itr);
      if(!inTrack) continue;
      TString tCharge;
      if(inTrack->Charge()>0) tCharge = "+";
      else if(inTrack->Charge()<0) tCharge = "-";
      cout << itr << " (" << tCharge << ") - " << GetNMatchedClusterIDsForTrack(fInputEvent,inTrack->GetID(),5,-5,0.2,-0.4) << "\t\t";
    }
    cout << endl;
    multimap<Int_t,Int_t>::iterator it;
    for (it=fMapTrackToCluster.begin(); it!=fMapTrackToCluster.end(); ++it) cout << it->first << " => " << it->second << '\n';
    cout << "mapClusterToTrack" << endl;
    Int_t tempClus = it->second;
    for (it=fMapClusterToTrack.begin(); it!=fMapClusterToTrack.end(); ++it) cout << it->first << " => " << it->second << '\n';
    vector<Int_t> tempTracks = GetMatchedTrackIDsForCluster(fInputEvent,tempClus, 5, -5, 0.2, -0.4);
    for(UInt_t iJ=0; iJ<tempTracks.size();iJ++){
      cout << tempClus << " - " << tempTracks.at(iJ) << endl;
    }
  }
  return;
}
