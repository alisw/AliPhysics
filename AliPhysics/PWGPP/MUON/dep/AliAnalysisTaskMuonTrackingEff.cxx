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

//Class to calculate the intrinsic efficiency of the detection elements of the
//MUON tracking chambers in function of the position in the detection element.
//Work on ESD only
//Author:  Nicolas LE BRIS - SUBATECH Nantes
// Modified by Matthieu LENHARDT - SUBATECH Nantes
// Modified by Antoine LARDEUX - SUBATECH Nantes
// Modified by Philippe PILLOT - SUBATECH Nantes

// ROOT includes
#include <TROOT.h>
#include <TList.h>
#include <THnSparse.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoGlobalMagField.h>
#include <TVectorD.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TFile.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliInputEventHandler.h"
#include "AliCounterCollection.h"

// ANALYSIS includes
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskMuonTrackingEff.h"

//MUON includes
#include "AliMUONCDB.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONRecoParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONConstants.h"
#include "AliMUON2DMap.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONTrackerData.h"

//include MUON/mapping:
#include "AliMpDEManager.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
#include "AliMpArea.h"
#include "AliMpDEIterator.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"

ClassImp(AliAnalysisTaskMuonTrackingEff)

const Int_t AliAnalysisTaskMuonTrackingEff::fgkNofDE[11] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26, 156};
const Int_t AliAnalysisTaskMuonTrackingEff::fgkNofBusPath = 888;
const Int_t AliAnalysisTaskMuonTrackingEff::fgkNofManu = 16828;

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff() :
  AliAnalysisTaskSE(),
  fOCDBLoaded(kFALSE),
  fOCDBpath(""),
  fAlignOCDBpath(""),
  fRecoParamOCDBpath(""),
  fCentMin(-FLT_MAX),
  fCentMax(FLT_MAX),
  fMuonTrackCuts(0x0),
  fPtCut(-1.),
  fUseMCLabel(kFALSE),
  fEnableDisplay(kFALSE),
  fTransformer(0x0),
  fDEPlanes(0x0),
  fClusters(0x0),
  fEvents(0x0),
  fChamberTDHistList(0x0),
  fChamberTTHistList(0x0),
  fChamberSDHistList(0x0),
  fExtraHistList(0x0)
{
  /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::AliAnalysisTaskMuonTrackingEff(TString name) :
  AliAnalysisTaskSE(name),
  fOCDBLoaded(kFALSE),
  fOCDBpath("raw://"),
  fAlignOCDBpath(""),
  fRecoParamOCDBpath(""),
  fCentMin(-FLT_MAX),
  fCentMax(FLT_MAX),
  fMuonTrackCuts(0x0),
  fPtCut(-1.),
  fUseMCLabel(kFALSE),
  fEnableDisplay(kFALSE),
  fTransformer(0x0),
  fDEPlanes(0x0),
  fClusters(0x0),
  fEvents(0x0),
  fChamberTDHistList(0x0),
  fChamberTTHistList(0x0),
  fChamberSDHistList(0x0),
  fExtraHistList(0x0)
{
  /// Constructor
  
  // Output slots
  DefineOutput(1, AliCounterCollection::Class());
  DefineOutput(2, AliCounterCollection::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMuonTrackingEff::~AliAnalysisTaskMuonTrackingEff()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fMuonTrackCuts;
    delete fClusters;
    delete fEvents;
    delete fChamberTDHistList;
    delete fChamberTTHistList;
    delete fChamberSDHistList;
    delete fExtraHistList;
  }
  delete fTransformer;
  delete fDEPlanes;
}

//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::NotifyRun()
{
  /// Load the OCDB and the Geometry
  
  // Load it only once
  if (fOCDBLoaded) return;
  
  // OCDB
  AliCDBManager* man = AliCDBManager::Instance();
  if (man->IsDefaultStorageSet()) printf("EfficiencyTask: CDB default storage already set!\n");
  else {
    man->SetDefaultStorage(fOCDBpath.Data());
    if (!fAlignOCDBpath.IsNull()) man->SetSpecificStorage("MUON/Align/Data",fAlignOCDBpath.Data());
    if (!fRecoParamOCDBpath.IsNull()) man->SetSpecificStorage("MUON/Calib/RecoParam",fRecoParamOCDBpath.Data());
  }
  if (man->GetRun() > -1) printf("EfficiencyTask: run number already set!\n");
  else man->SetRun(fCurrentRunNumber);
  
  // Geometry
  if (!AliGeomManager::GetGeometry()) {
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;  
    if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return;
  }
  fTransformer = new AliMUONGeometryTransformer();
  fTransformer->LoadGeometryData();
  
  // Mapping
  if (!AliMpSegmentation::Instance(kFALSE) || !AliMpDDLStore::Instance(kFALSE)) {
    if (!AliMUONCDB::LoadMapping()) return;
  }
  
  // vectors (x0, y0, z0, a, b, c) defining the plane of each DE in the global frame
  Double_t pl0[3] = {0., 0., 0.};
  Double_t pl1[3] = {0., 0., 1.};
  Double_t pg0[3], pg1[3];
  fDEPlanes = new TObjArray(1026);
  fDEPlanes->SetOwner(kTRUE);
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    AliMpDEIterator it;
    it.First(i);
    while (!it.IsDone()) {
      Int_t deId = it.CurrentDEId();
      fTransformer->Local2Global(deId, pl0[0], pl0[1], pl0[2], pg0[0], pg0[1], pg0[2]);
      fTransformer->Local2Global(deId, pl1[0], pl1[1], pl1[2], pg1[0], pg1[1], pg1[2]);
      TVectorD *plane = new TVectorD(6);
      (*plane)[0] = pg0[0];
      (*plane)[1] = pg0[1];
      (*plane)[2] = pg0[2];
      (*plane)[3] = pg1[0] - pg0[0];
      (*plane)[4] = pg1[1] - pg0[1];
      (*plane)[5] = pg1[2] - pg0[2];
      fDEPlanes->AddAt(plane, deId);
      it.Next();
    }
  }
  
  // Prepare the tracker (will load RecoParam and magnetic field if not already done)
  if (!AliMUONESDInterface::GetTracker()) AliMUONESDInterface::ResetTracker();
  
  // get the trackCuts for this run
  if (!fMuonTrackCuts) AliFatal("You must specify the requested selections (AliMuonTrackCut obj is missing)");
  fMuonTrackCuts->SetRun(fInputHandler);
  fMuonTrackCuts->Print();
  
  fOCDBLoaded = kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::UserCreateOutputObjects()
{
  /// Define output objects
  
  // count detected, accepted and expected clusters
  fClusters = new AliCounterCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fClusters->AddRubric("Cluster", "Detected/Accepted/Expected");
  fClusters->AddRubric("Chamber", AliMpConstants::NofTrackingChambers());
  fClusters->AddRubric("DE", fgkNofDE[10]);
  fClusters->AddRubric("BusPatch", fgkNofBusPath);
  fClusters->AddRubric("Manu", fgkNofManu);
  fClusters->AddRubric("Channel", AliMpConstants::ManuNofChannels());
  fClusters->Init();
  
  // events analyzed
  fEvents = new AliCounterCollection(GetOutputSlot(2)->GetContainer()->GetName());
  fEvents->AddRubric("event", "any");
  fEvents->AddRubric("run", 100000);
  fEvents->Init();
  
  fChamberTDHistList = new TList();
  fChamberTDHistList->SetOwner();
  fChamberTTHistList = new TList();
  fChamberTTHistList->SetOwner();
  fChamberSDHistList = new TList();
  fChamberSDHistList->SetOwner();
  fExtraHistList = new TList();
  fExtraHistList->SetOwner();
  
  THnSparse *hn = 0x0;
  TString histName, histTitle;
  
  // centrality bins
  Int_t nCentBins = 22;
  Double_t centRange[2] = {-5., 105.};
  
  // prepare binning for THnSparse
  // 1: Ch or DE Id
  // 2: centrality
  // 3: pt
  // 4: y
  // 5: phi
  // 6: sign
  const Int_t nDims = 6;
  Int_t nBins[nDims] = {0, nCentBins, 20, 15, 15, 2};
  Double_t xMin[nDims] = {0., centRange[0], 0., -4., 0., -2.};
  Double_t xMax[nDims] = {0., centRange[1], 20., -2.5, TMath::TwoPi(), 2.};
  
  for (Int_t iCh = 0; iCh < 10; iCh++)
  {
    // histograms per chamber
    nBins[0] = fgkNofDE[iCh];
    xMin[0] = 0.; xMax[0] = static_cast<Double_t>(fgkNofDE[iCh]);
    histTitle.Form("ChamberNbr %d", iCh+1);
    histName.Form("TD_ChamberNbr%d", iCh+1);
    hn = new THnSparseT<TArrayF>(histName, histTitle, nDims, nBins, xMin, xMax);
    fChamberTDHistList->AddAt(hn, iCh);
    histName.Form("TT_ChamberNbr%d",iCh+1);
    hn = new THnSparseT<TArrayF>(histName, histTitle, nDims, nBins, xMin, xMax);
    fChamberTTHistList->AddAt(hn, iCh);
    histName.Form("SD_ChamberNbr%d", iCh+1);
    hn = new THnSparseT<TArrayF>(histName, histTitle, nDims, nBins, xMin, xMax);
    fChamberSDHistList->AddAt(hn, iCh);
    
  }
  
  // global histograms per chamber
  nBins[0] = 10;
  xMin[0] = 0.5; xMax[0] = 10.5;
  hn = new THnSparseT<TArrayF>("TD_Chambers 11", "Chambers 11", nDims, nBins, xMin, xMax);
  fChamberTDHistList->AddAt(hn, 10);
  hn = new THnSparseT<TArrayF>("TT_Chambers 11", "Chambers 11", nDims, nBins, xMin, xMax);
  fChamberTTHistList->AddAt(hn, 10);
  hn = new THnSparseT<TArrayF>("SD_Chambers 11", "Chambers 11", nDims, nBins, xMin, xMax);
  fChamberSDHistList->AddAt(hn, 10);

  //Extra histograms
  TH1F *fHistCent = new TH1F("fHistCent", "centrality distribution", nCentBins, centRange[0], centRange[1]);
  fExtraHistList->AddAt(fHistCent,0);
  TH1F *fHistPt = new TH1F("fHistPt", "pt distribution", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPt,1);
  TH1F *fHistY = new TH1F("fHistY", "y distribution", 60, -4., -2.5);
  fExtraHistList->AddAt(fHistY,2);
  TH1F *fHistTheta = new TH1F("fHistTheta", "theta distribution", 120, 2.8, 3.2);
  fExtraHistList->AddAt(fHistTheta,3);
  TH1F *fHistP = new TH1F("fHistP", "momentum distribution", 250, 0., 500.);
  fExtraHistList->AddAt(fHistP,4);
  TH1F *fHistZVtx = new TH1F("fHistZVtx", "Z vertex distribution", 200, -50., 50.);
  fExtraHistList->AddAt(fHistZVtx,5);
  TH1F *fHistPhi = new TH1F("fHistPhi", "phi distribution", 60, 0., TMath::TwoPi());
  fExtraHistList->AddAt(fHistPhi,6);
  TH1F *fHistPtRap2p5To2p75 = new TH1F("fHistPtRap2p5To2p75", "2.5 < y < 2.75", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPtRap2p5To2p75,7);
  TH1F *fHistPtRap2p75To3p0 = new TH1F("fHistPtRap2p75To3p0", "2.75 < y < 3.0", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPtRap2p75To3p0,8);
  TH1F *fHistPtRap3p0To3p25 = new TH1F("fHistPtRap3p0To3p25", "3.0 < y < 3.25", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPtRap3p0To3p25,9);
  TH1F *fHistPtRap3p25To3p5 = new TH1F("fHistPtRap3p25To3p5", "3.25 < y < 3.5", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPtRap3p25To3p5,10);
  TH1F *fHistPtRap3p5To3p75 = new TH1F("fHistPtRap3p5To3p75", "3.5 < y < 3.75", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPtRap3p5To3p75,11);
  TH1F *fHistPtRap3p75To4p0 = new TH1F("fHistPtRap3p75To4p0", "3.75 < y < 4.0", 250, 0., 50.);
  fExtraHistList->AddAt(fHistPtRap3p75To4p0,12);
  TH1F *fHistRapPt1p0To2p0 = new TH1F("fHistRapPt1p0To2p0", "1.0 < pT < 2.0", 60, -4., -2.5);
  fExtraHistList->AddAt(fHistRapPt1p0To2p0,13);
  TH1F *fHistRapPt2p0To5p0 = new TH1F("fHistRapPt2p0To5p0", "2.0 < pT < 5.0", 60, -4., -2.5);
  fExtraHistList->AddAt(fHistRapPt2p0To5p0,14);
  TH1F *fHistRapPt5p0To8p0 = new TH1F("fHistRapPt5p0To8p0", "5.0 < pT < 8.0", 60, -4., -2.5);
  fExtraHistList->AddAt(fHistRapPt5p0To8p0,15);
  TH1F *fHistRapPt2p0To4p0 = new TH1F("fHistRapPt2p0To4p0", "2.0 < pT < 4.0", 60, -4., -2.5);
  fExtraHistList->AddAt(fHistRapPt2p0To4p0,16);
  TH1F *fHistRapPt4p0To8p0 = new TH1F("fHistRapPt4p0To8p0", "4.0 < pT < 8.0", 60, -4., -2.5);
  fExtraHistList->AddAt(fHistRapPt4p0To8p0,17);
  
  TH2D *hDXYOverDXYMax = new TH2D("hDXYOverDXYMax", "DXY / DXYMax;DX / DXMax;DY / DYMax", 100, -1., 1., 100, -1., 1.);
  fExtraHistList->AddAt(hDXYOverDXYMax,18);
  TH2D *hDXOverDXMax = new TH2D("hDXOverDXMax", "DX / DXMax vs pXZ;pXZ;DX / DXMax", 50, 0., 500., 100, -1., 1.);
  fExtraHistList->AddAt(hDXOverDXMax,19);
  TH2D *hDYOverDYMax = new TH2D("hDYOverDYMax", "DY / DYMax vs pYZ;pYZ;DY / DYMax", 50, 0., 500., 100, -1., 1.);
  fExtraHistList->AddAt(hDYOverDYMax,20);
  
  THnSparse *hKine = new THnSparseT<TArrayF>("hKine", "kinematics distribution", nDims-1, &nBins[1], &xMin[1], &xMax[1]);
  fExtraHistList->AddAt(hKine,21);
  
  TH1F *fHistXVtx = new TH1F("fHistXVtx", "X vertex distribution", 200, -1., 1.);
  fExtraHistList->AddAt(fHistXVtx,22);
  TH1F *fHistYVtx = new TH1F("fHistYVtx", "Y vertex distribution", 200, -1., 1.);
  fExtraHistList->AddAt(fHistYVtx,23);
  
  // post the output data at least once
  PostData(1, fClusters);  
  PostData(2, fEvents);
  PostData(3, fChamberTDHistList);
  PostData(4, fChamberTTHistList);
  PostData(5, fChamberSDHistList);
  PostData(6, fExtraHistList);
}

//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::UserExec(Option_t *)
{
  /// Main event loop
  
  // check the OCDB has been loaded properly
  if (!fOCDBLoaded) AliFatal("Problem occur while loading OCDB objects");
  
  // get the current event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  // get the centrality
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(esd->FindListObject("MultSelection"));
  Double_t cent = multSelection ? multSelection->GetMultiplicityPercentile("V0M") : -1.;
  if (cent < fCentMin || cent > fCentMax) return;
  static_cast<TH1F*>(fExtraHistList->At(0))->Fill(cent);
  
  // total number of events analyzed
  fEvents->Count(Form("event:any/run:%d",fCurrentRunNumber));
  
  // loop over tracks
  AliMUONTrack track;
  Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack *esdTrack = esd->GetMuonTrack(iTrack);
    
    // apply track selections
    Double_t pT = esdTrack->Pt();
    if (!esdTrack->ContainTrackerData() || !fMuonTrackCuts->IsSelected(esdTrack) || (fPtCut > 0. && pT < fPtCut) ||
	(fUseMCLabel && (esdTrack->GetLabel() < 0 || esdTrack->TestBit(BIT(22))))) continue;
    
    // fill histograms
    Double_t y = esdTrack->Y();
    Double_t phi = esdTrack->Phi();
    static_cast<TH1F*>(fExtraHistList->At(1))->Fill(pT);
    static_cast<TH1F*>(fExtraHistList->At(2))->Fill(y);
    static_cast<TH1F*>(fExtraHistList->At(3))->Fill(esdTrack->Theta());
    static_cast<TH1F*>(fExtraHistList->At(4))->Fill(esdTrack->P());
    static_cast<TH1F*>(fExtraHistList->At(5))->Fill(esdTrack->GetZ());
    static_cast<TH1F*>(fExtraHistList->At(22))->Fill(esdTrack->GetNonBendingCoor());
    static_cast<TH1F*>(fExtraHistList->At(23))->Fill(esdTrack->GetBendingCoor());
    static_cast<TH1F*>(fExtraHistList->At(6))->Fill(phi);
    if ( (y < -2.5) && (y > -2.75) ) static_cast<TH1F*>(fExtraHistList->At(7))->Fill(pT);
    if ( (y < -2.75) && (y > -3.0) ) static_cast<TH1F*>(fExtraHistList->At(8))->Fill(pT);
    if ( (y < -3.0) && (y > -3.25) ) static_cast<TH1F*>(fExtraHistList->At(9))->Fill(pT);
    if ( (y < -3.25) && (y > -3.5) ) static_cast<TH1F*>(fExtraHistList->At(10))->Fill(pT);
    if ( (y < -3.5) && (y > -3.75) ) static_cast<TH1F*>(fExtraHistList->At(11))->Fill(pT);
    if ( (y < -3.75) && (y > -4.0) ) static_cast<TH1F*>(fExtraHistList->At(12))->Fill(pT);
    if ( (pT > 1.0) && (pT < 2.0) ) static_cast<TH1F*>(fExtraHistList->At(13))->Fill(y);
    if ( (pT > 2.0) && (pT < 5.0) ) static_cast<TH1F*>(fExtraHistList->At(14))->Fill(y);
    if ( (pT > 5.0) && (pT < 8.0) ) static_cast<TH1F*>(fExtraHistList->At(15))->Fill(y);
    if ( (pT > 2.0) && (pT < 4.0) ) static_cast<TH1F*>(fExtraHistList->At(16))->Fill(y);
    if ( (pT > 4.0) && (pT < 8.0) ) static_cast<TH1F*>(fExtraHistList->At(17))->Fill(y);
    
    // convert to MUON track
    AliMUONESDInterface::ESDToMUON(*esdTrack, track);
    Double_t trackInfo[6] = {0., cent, pT, y, phi, static_cast<Double_t>(esdTrack->Charge())};
    static_cast<THnSparse*>(fExtraHistList->At(21))->Fill(&trackInfo[1]);
    
    // tag the removable clusters/chambers, i.e. not needed to fulfill the tracking conditions
    Bool_t removableChambers[10];
    Bool_t isValidTrack = TagRemovableClusters(track, removableChambers);
    
    // loop over clusters
    Int_t previousCh = -1;
    TObjArray *trackParams = track.GetTrackParamAtCluster();
    AliMUONTrackParam *trackParam = static_cast<AliMUONTrackParam*>(trackParams->First());
    while (trackParam) {
      
      AliMUONVCluster* cluster = trackParam->GetClusterPtr();
      Int_t currentCh = cluster->GetChamberId();
      Int_t currentDE = cluster->GetDetElemId();
      
      // find the pads at the position of the cluster
      Double_t pos[3] = {cluster->GetX(), cluster->GetY(), cluster->GetZ()};
      AliMpPad pad[2];
      FindPads(currentDE, pos, pad);
      
      AliMUONTrackParam *nextTrackParam = static_cast<AliMUONTrackParam*>(trackParams->After(trackParam));
      Int_t nextCh = nextTrackParam ? nextTrackParam->GetClusterPtr()->GetChamberId() : 10;
      
      // record all clusters/chambers
      RecordCluster(currentCh, currentDE, pad, trackInfo, "Detected", fChamberSDHistList, (currentCh != nextCh));
      
      // record removable clusters/chambers
      if (trackParam->IsRemovable()) {
	Bool_t recordChamber = (removableChambers[currentCh] && currentCh != nextCh);
	RecordCluster(currentCh, currentDE, pad, trackInfo, "Accepted", fChamberTDHistList, recordChamber);
	RecordCluster(currentCh, currentDE, pad, trackInfo, "Expected", fChamberTTHistList, recordChamber);
      }
      
      // record missing clusters/chambers prior to current one
      while (previousCh < currentCh-1 && (isValidTrack || previousCh < 5))
	FindAndRecordMissingClusters(*trackParam, ++previousCh, trackInfo);
      
      // stop if we reached station 4 or 5 and the track is not valid
      if (!isValidTrack && currentCh > 5) break;
      
      // record missing cluster on the same chamber
      if (currentCh != previousCh && currentCh != nextCh)
	FindAndRecordMissingClusters(*trackParam, currentCh, trackInfo);
      
      if (nextTrackParam) {
	
	// prepare next step
	previousCh = currentCh;
	trackParam = nextTrackParam;
	
      } else {
	
	// record missing clusters/chambers next to the last chamber
	while (++currentCh < 10 && (isValidTrack || currentCh < 6))
	  FindAndRecordMissingClusters(*trackParam, currentCh, trackInfo);
	
	break;
      }
      
    }
    
  }
  
  // post the output data:
  PostData(1, fClusters);
  PostData(2, fEvents);
  PostData(3, fChamberTDHistList);
  PostData(4, fChamberTTHistList);
  PostData(5, fChamberSDHistList);
  PostData(6, fExtraHistList);
}
    
//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Terminate(Option_t *)
{
  /// final plots
  
  if (!fEnableDisplay) return;
  
  fClusters = static_cast<AliCounterCollection*>(GetOutputData(1));
  if (!fClusters) return;
  
  // load mapping locally if not already done
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    if (gROOT->IsBatch()) man->SetDefaultStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
    else man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }
  if (man->GetRun() < 0) man->SetRun(0);
  if (!AliMpSegmentation::Instance(kFALSE) || !AliMpDDLStore::Instance(kFALSE)) {
    if (!AliMUONCDB::LoadMapping()) return;
  }
  
  TString clusterKey[3] = {"Detected", "Accepted", "Expected"};
  
  // list of fired DEs
  TObjArray *deKeys = fClusters->GetKeyWords("DE").Tokenize(",");
  Int_t nDEs = deKeys->GetEntriesFast();
  
  // loop over cluster types
  for (Int_t iKey = 0; iKey < 3; iKey++) {
    
    // create the cluster store
    AliMUON2DMap clustersStore(kTRUE);
    
    // loop over fired DEs
    for (Int_t iDE = 0; iDE < nDEs; iDE++) {
      
      // get DE Id
      TString deKey = static_cast<TObjString*>(deKeys->UncheckedAt(iDE))->GetString();
      Int_t deId = deKey.Atoi();
      
      // get numbers of clusters in this DE for each manu/channel combination
      TH2D *channelVsManu = fClusters->Get("channel", "Manu",
					   Form("Cluster:%s/DE:%s", clusterKey[iKey].Data(), deKey.Data()));
      Int_t nManus = channelVsManu->GetNbinsX();
      Int_t nChannels = channelVsManu->GetNbinsY();
      
      // loop over fired manus
      for (Int_t iManu = 1; iManu <= nManus; iManu++) {
	
	// get manu Id
	TString manuKey = channelVsManu->GetXaxis()->GetBinLabel(iManu);
	Int_t manuId = manuKey.Atoi();
	
	// loop over fired channels
	for (Int_t iChannel = 1; iChannel <= nChannels; iChannel++) {
	  
	  // get channel Id
	  TString channelKey = channelVsManu->GetYaxis()->GetBinLabel(iChannel);
	  Int_t channelId = channelKey.Atoi();
	  
	  // get the number of clusters in this pad
	  Int_t nClusters = static_cast<Int_t>(channelVsManu->GetBinContent(iManu, iChannel));
	  if (nClusters < 1) continue;
	  
	  // register the clusters
	  AliMUONVCalibParam* c = static_cast<AliMUONVCalibParam*>(clustersStore.FindObject(deId, manuId));
	  if (!c) {
	    c = new AliMUONCalibParamNI(1, AliMpConstants::ManuNofChannels(), deId, manuId);
	    clustersStore.Add(c);
	  }
	  c->SetValueAsInt(channelId, 0, nClusters);
	  
	}
	
      }
      
      // clean memory
      delete channelVsManu;
      
    }
    
    // create the tracker data
    TString suffix = GetName();
    suffix.ReplaceAll("MuonTrackingEfficiency","");
    AliMUONTrackerData clustersData(Form("%sClusters%s", clusterKey[iKey].Data(), suffix.Data()),
				    Form("%s clusters %s", clusterKey[iKey].Data(), suffix.Data()), 1, kFALSE);
    clustersData.SetDimensionName(0, "count");
    clustersData.Add(clustersStore);
    
    // save it to a file
    TFile *outFile = TFile::Open("DisplayResults.root", "UPDATE");
    if (outFile && outFile->IsOpen()) {
      clustersData.Write(0x0, TObject::kOverwrite);
      outFile->Close();
    }
    
  }
  
  // clean memory
  delete deKeys;
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonTrackingEff::TagRemovableClusters(AliMUONTrack &track, Bool_t removableChambers[10])
{
  /// Identify clusters/chambers that can be removed from the track
  /// return kTRUE if the track as it is satisfies the tracking conditions
  
  for (Int_t i = 0; i < 10; i++) removableChambers[i] = kFALSE;
  
  // check if track is valid as it is
  UInt_t requestedStationMask = AliMUONESDInterface::GetTracker()->GetRecoParam()->RequestedStationMask();
  Bool_t request2ChInSameSt45 = !AliMUONESDInterface::GetTracker()->GetRecoParam()->MakeMoreTrackCandidates();
  Bool_t isValidTrack = track.IsValid(requestedStationMask, request2ChInSameSt45);
  
  // count the number of clusters per chamber and the number of chambers hit per station
  Int_t nClInCh[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Int_t nChInSt[5] = {0, 0, 0, 0, 0};
  Int_t previousCh = -1;
  Int_t nClusters = track.GetNClusters();
  for (Int_t i = 0; i < nClusters; i++) {
    
    AliMUONTrackParam *trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(i));
    Int_t currentCh = trackParam->GetClusterPtr()->GetChamberId();
    Int_t currentSt = currentCh/2;
    
    nClInCh[currentCh]++;
    
    if (currentCh != previousCh) {
      previousCh = currentCh;
      nChInSt[currentSt]++;
    }
    
  }
  
  // tag removable clusters/chambers
  for (Int_t i = 0; i < nClusters; i++) {
    
    AliMUONTrackParam *trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(i));
    Int_t currentCh = trackParam->GetClusterPtr()->GetChamberId();
    Int_t currentSt = currentCh/2;
    
    // for stations 1, 2 and 3
    if (currentSt < 3) {
      
      if (!((1 << currentSt) & requestedStationMask) || // station not required
	  nChInSt[currentSt] == 2) { // or both chamber hit in the station
	
	removableChambers[currentCh] = kTRUE;
	trackParam->SetRemovable(kTRUE);
	
      } else if (nClInCh[currentCh] == 2) { // 2 clusters in the chamber
	
	trackParam->SetRemovable(kTRUE);
	
      }
      
    } else { // for stations 4 and 5
      
      // if the track is already not valid we can certainly not remove more cluster on station 4 or 5
      if (!isValidTrack) continue;
      
      if (((request2ChInSameSt45 && // tracking starts with 2 chamber hits in the same station:
	    nChInSt[4-currentSt/4] == 2) || // --> 2 hits in the other stations
	   (!request2ChInSameSt45 && // or tracking starts with 2 chamber hits in stations 4&5:
	    nChInSt[3]+nChInSt[4] >= 3)) && // --> at least 3 hits in stations 4&5
	  (!((1 << currentSt) & requestedStationMask) || // + this station not requested
	   nChInSt[currentSt] == 2)) { // or 2 hits in it
        
        removableChambers[currentCh] = kTRUE;
        trackParam->SetRemovable(kTRUE);
        
      } else if (nClInCh[currentCh] == 2) { // 2 clusters in the chamber
	
	trackParam->SetRemovable(kTRUE);
	
      }
      
    }
    
  }
  
  return isValidTrack;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::FindAndRecordMissingClusters(AliMUONTrackParam &param, Int_t chamber,
								  Double_t trackInfo[6])
{
  /// Find which detection elements should have been hit and record the missing clusters
  
  static const Double_t maxDZ = 0.01; // max distance between extrapolated track and DE to stop to extrapolate
  static const Double_t maxDevX = 0.01; // max X-deviation per dZ(cm) between the linear and the correct extrapolation
  static const Double_t maxDevY = 0.1; // max Y-deviation per dZ(cm) between the linear and the correct extrapolation
  static const Double_t minDXY = 1.; // min half-size of track area to intersect with DE to account for bad DE area
  
  Bool_t missingChamber = (chamber != param.GetClusterPtr()->GetChamberId());
  Int_t startDE = param.GetClusterPtr()->GetDetElemId();
  Double_t pos[3], maxDX, maxDY, dZ = 0.;
  AliMUONTrackParam param1, param2;
  AliMpPad pad[2];
  
  // extrapolate the track to the missing chamber if needed
  param1.SetParameters(param.GetParameters());
  param1.SetZ(param.GetZ());
  if (missingChamber && !AliMUONTrackExtrap::ExtrapToZ(&param1, AliMUONConstants::DefaultChamberZ(chamber))) return;
  
  Double_t pX = param1.Px();
  Double_t pY = param1.Py();
  Double_t pZ = param1.Pz();
  Double_t pXZ = TMath::Sqrt(pX*pX + pZ*pZ);
  Double_t pYZ = TMath::Sqrt(pY*pY + pZ*pZ);
  
  // loop over DEs
  AliMpDEIterator it;
  it.First(chamber);
  while (!it.IsDone()) {
    Int_t deId = it.CurrentDEId();
    
    // skip current cluster
    if (deId == startDE) {
      it.Next();
      continue;
    }
    
    // reset parameters
    param2.SetParameters(param1.GetParameters());
    param2.SetZ(param1.GetZ());
    
    // check if the track can cross this DE
    Int_t nStep = 0;
    Bool_t crossDE = kFALSE;
    Bool_t extrapOk = kTRUE;
    do {
      
      // plots to check that the correct extrapolation is within the area opened around the linear extrapolation
      if (nStep >= 1) {
	Double_t dX = pos[0]-param2.GetNonBendingCoor();
	Double_t dY = pos[1]-param2.GetBendingCoor();
	Double_t dXOverDXMax = (dZ > 0.) ? dX*pXZ/(maxDevX*dZ) : 0.;
	Double_t dYOverDYMax = (dZ > 0.) ? dY*pYZ/(maxDevY*dZ) : 0.;
	static_cast<TH2D*>(fExtraHistList->At(18))->Fill(dXOverDXMax, dYOverDYMax);
	static_cast<TH2D*>(fExtraHistList->At(19))->Fill(pXZ, dXOverDXMax);
	static_cast<TH2D*>(fExtraHistList->At(20))->Fill(pYZ, dYOverDYMax);
	if (dXOverDXMax >= 1.) printf("st = %d; pXZ = %f; dZ = %f; dX = %f; dX*PXZ/(dZ*maxDevX) = %f\n",
				      chamber/2, pXZ, dZ, dX, dXOverDXMax);
	if (dYOverDYMax >= 1.) printf("st = %d; pYZ = %f; dZ = %f; dY = %f; dY*PYZ/(dZ*maxDevY) = %f\n",
				      chamber/2, pYZ, dZ, dY, dYOverDYMax);
      }
      nStep++;
      
      // build the area in which the track can eventually cross this DE and check if it overlaps with the DE area
      Intersect(param2, deId, pos);
      dZ = TMath::Abs(pos[2]-param2.GetZ());
      maxDX = minDXY + maxDevX*dZ/pXZ;
      maxDY = minDXY + maxDevY*dZ/pYZ;
      AliMpArea area(pos[0], pos[1], maxDX, maxDY);
      crossDE = OverlapDE(area, deId);
      
    } while(crossDE && dZ > maxDZ && (extrapOk = AliMUONTrackExtrap::ExtrapToZ(&param2, pos[2])));
    
    // find the pads (if any) at the position of the missing cluster and register it
    if (crossDE && extrapOk && FindPads(deId, pos, pad)) {
      RecordCluster(chamber, deId, pad, trackInfo, "Expected", fChamberTTHistList, missingChamber);
      missingChamber = kFALSE;
    }
    
    it.Next();
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::Intersect(AliMUONTrackParam &param, Int_t deId, Double_t p[3])
{
  /// Find the intersection point between the track (assuming straight line) and the DE in the global frame
  
  Double_t pos[3] = {param.GetNonBendingCoor(), param.GetBendingCoor(), param.GetZ()};
  Double_t slope[2] = {param.GetNonBendingSlope(), param.GetBendingSlope()};
  TVectorD &plane = *(static_cast<TVectorD*>(fDEPlanes->UncheckedAt(deId)));
  
  p[2] = (plane[3]*(slope[0]*pos[2]-pos[0]+plane[0]) + plane[4]*(slope[1]*pos[2]-pos[1]+plane[1]) + plane[5]*plane[2]) /
  (plane[3]*slope[0] + plane[4]*slope[1] + plane[5]);
  p[0] = slope[0]*(p[2]-pos[2]) + pos[0];
  p[1] = slope[1]*(p[2]-pos[2]) + pos[1];
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonTrackingEff::OverlapDE(AliMpArea &area, Int_t deId)
{
  /// Check whether (global) area overlaps with the given DE
  
  AliMpArea* globalDEArea = fTransformer->GetDEArea(deId);
  if (!globalDEArea) return kFALSE;
  
  return area.Overlap(*globalDEArea);
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonTrackingEff::RecordCluster(Int_t chamber, Int_t deId, AliMpPad pad[2], Double_t trackInfo[6],
						   TString clusterKey, TList *chamberHistList, Bool_t recordChamber)
{
  /// Register the cluster in the given stores
  
  // register the pads
  for (Int_t iCath = 0; iCath < 2; iCath++) if (pad[iCath].IsValid()) {
    Int_t manuId = pad[iCath].GetManuId();
    Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(deId,manuId);
    fClusters->Count(Form("Cluster:%s/Chamber:%d/DE:%d/BusPatch:%d/Manu:%d/Channel:%d",
		     clusterKey.Data(), chamber, deId, busPatchId, manuId, pad[iCath].GetManuChannel()));
  }
  
  // register the DE
  trackInfo[0] = static_cast<Double_t>(deId%100);
  static_cast<THnSparse*>(chamberHistList->At(chamber))->Fill(trackInfo);
  
  // register the chamber
  if (recordChamber) {
    trackInfo[0] = static_cast<Double_t>(chamber+1);
    static_cast<THnSparse*>(chamberHistList->At(10))->Fill(trackInfo);
  }
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonTrackingEff::FindPads(Int_t deId, Double_t pos[3], AliMpPad pad[2])
{
  /// Look for pads at the cluster's location
  
  static const AliMpPad emptyPad;
  
  // compute the cluster position in the DE frame
  Double_t localPos[3];
  fTransformer->Global2Local(deId, pos[0], pos[1], pos[2], localPos[0], localPos[1], localPos[2]);
  
  // find pads at this position
  for (Int_t iCath = 0; iCath < 2; iCath++) {
    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::GetCathodType(iCath));
    if (seg) pad[iCath] = seg->PadByPosition(localPos[0], localPos[1], kFALSE);
    else pad[iCath] = emptyPad;
  }
  
  return (pad[0].IsValid() || pad[1].IsValid());
  
}

