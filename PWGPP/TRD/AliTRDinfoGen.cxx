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

/* $Id: AliTRDinfoGen.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//
//  Tender wagon for TRD performance/calibration train
//
// In this wagon the information from
//   - ESD
//   - Friends [if available]
//   - MC [if available]
// are grouped into AliTRDtrackInfo objects and fed to worker tasks
//
//  Authors:
//    Markus Fasel <M.Fasel@gsi.de>
//    Alexandru Bercuci <A.Bercuci@gsi.de>
//
////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TString.h>
#include <TVectorT.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliMCEventHandler.h"

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDHeader.h"
#include "AliESDRun.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0KineCuts.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "TTreeStream.h"

#include <cstdio>
#include <climits>
#include <cstring>
#include <iostream>
#include <memory>

#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"
#include "AliTRDcalibDB.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDinfoGen.h"
#include "AliTRDpwgppHelper.h"
#include "info/AliTRDtrackInfo.h"
#include "info/AliTRDeventInfo.h"
#include "info/AliTRDv0Info.h"
#include "info/AliTRDchmbInfo.h"
#include "info/AliTRDtriggerInfo.h"
#include "info/AliTRDeventCuts.h"

ClassImp(AliTRDinfoGen)

const Float_t AliTRDinfoGen::fgkITS = 100.; // to be checked
const Float_t AliTRDinfoGen::fgkTPC = 290.;
const Float_t AliTRDinfoGen::fgkTRD = 365.;

const Float_t AliTRDinfoGen::fgkTrkDCAxy  = 1.;
const Float_t AliTRDinfoGen::fgkTrkDCAz   = 1.;
const Int_t   AliTRDinfoGen::fgkNclTPC    = 70;
const Float_t AliTRDinfoGen::fgkPt        = 0.2;
const Float_t AliTRDinfoGen::fgkEta       = 0.9;
AliTRDReconstructor* AliTRDinfoGen::fgReconstructor(NULL);
AliTRDgeometry* AliTRDinfoGen::fgGeo(NULL);
//____________________________________________________________________
AliTRDinfoGen::AliTRDinfoGen()
  :AliAnalysisTaskSE()
  ,fESDev(NULL)
  ,fMCev(NULL)
  ,fEventCut(NULL)
  ,fTrackCut(NULL)
  ,fV0Identifier(NULL)
  ,fV0Cut(NULL)
  ,fOCDB("local://$ALICE_ROOT/OCDB")
  ,fTrackInfo(NULL)
  ,fEventInfo(NULL)
  ,fV0Info(NULL)
  ,fTracksBarrel(NULL)
  ,fTracksITS(NULL)
  ,fTracksSA(NULL)
  ,fTracksKink(NULL)
  ,fV0List(NULL)
  ,fClusters(NULL)
  ,fContainer(NULL)
  ,fRecos(NULL)
  ,fDebugStream(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDinfoGen", "MC-REC TRD-track list generator");
}

//____________________________________________________________________
AliTRDinfoGen::AliTRDinfoGen(char* name)
  :AliAnalysisTaskSE(name)
  ,fESDev(NULL)
  ,fMCev(NULL)
  ,fEventCut(NULL)
  ,fTrackCut(NULL)
  ,fV0Identifier(NULL)
  ,fV0Cut(NULL)
  ,fOCDB("local://$ALICE_ROOT/OCDB")
  ,fTrackInfo(NULL)
  ,fEventInfo(NULL)
  ,fV0Info(NULL)
  ,fTracksBarrel(NULL)
  ,fTracksITS(NULL)
  ,fTracksSA(NULL)
  ,fTracksKink(NULL)
  ,fV0List(NULL)
  ,fClusters(NULL)
  ,fContainer(NULL)
  ,fRecos(NULL)
  ,fDebugStream(NULL)
{
  //
  // Default constructor
  //
  SetTitle("MC-REC TRD-track list generator");
  DefineOutput(AliTRDpwgppHelper::kTracksBarrel, TObjArray::Class());
  DefineOutput(AliTRDpwgppHelper::kTracksITS,    TObjArray::Class());
  DefineOutput(AliTRDpwgppHelper::kTracksSA,     TObjArray::Class());
  DefineOutput(AliTRDpwgppHelper::kTracksKink,   TObjArray::Class());
  DefineOutput(AliTRDpwgppHelper::kEventInfo,    AliTRDeventInfo::Class());
  DefineOutput(AliTRDpwgppHelper::kV0List,       TObjArray::Class());
  DefineOutput(AliTRDpwgppHelper::kClusters,     TObjArray::Class());
  DefineOutput(AliTRDpwgppHelper::kMonitor,      TObjArray::Class()); // histogram list
}

//____________________________________________________________________
AliTRDinfoGen::~AliTRDinfoGen()
{
// Destructor
  if(fgGeo) delete fgGeo;
  if(fgReconstructor) delete fgReconstructor;
  if(fDebugStream) delete fDebugStream;
  if(fV0Cut) delete fV0Cut;
  if(fV0Identifier) delete fV0Identifier;
  if(fTrackCut) delete fTrackCut;
  if(fEventCut) delete fEventCut;
  if(fTrackInfo) delete fTrackInfo; fTrackInfo = NULL;
  if(fEventInfo) delete fEventInfo; fEventInfo = NULL;
  if(fV0Info) delete fV0Info; fV0Info = NULL;
  if(fTracksBarrel){
    fTracksBarrel->Delete(); delete fTracksBarrel;
    fTracksBarrel = NULL;
  }
  if(fTracksITS){
    fTracksITS->Delete(); delete fTracksITS;
    fTracksITS = NULL;
  }
  if(fTracksSA){
    fTracksSA->Delete(); delete fTracksSA;
    fTracksSA = NULL;
  }
  if(fTracksKink){ 
    fTracksKink->Delete(); delete fTracksKink;
    fTracksKink = NULL;
  }
  if(fV0List){ 
    fV0List->Delete(); 
    delete fV0List;
    fV0List = NULL;
  }
  if(fClusters){
    fClusters->Delete(); delete fClusters;
    fClusters = NULL;
  }
  if(fContainer && !(AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode())){
    fContainer->Delete(); 
    delete fContainer;
    fContainer = NULL;
  }
}

//____________________________________________________________________
Bool_t AliTRDinfoGen::GetRefFigure(Int_t)
{
// General graphs for PWGPP/TRD train
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }
  fContainer->At(kStatTrk)->Draw("bar");
  return kTRUE;
}

//____________________________________________________________________
void AliTRDinfoGen::UserCreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
 
  fTrackInfo = new AliTRDtrackInfo();
  fEventInfo = new AliTRDeventInfo();
  fV0Info    = new AliTRDv0Info();
  fTracksBarrel = new TObjArray(200); fTracksBarrel->SetOwner(kTRUE);
  fTracksITS    = new TObjArray(20); fTracksITS->SetOwner(kTRUE);
  fTracksSA     = new TObjArray(20); fTracksSA->SetOwner(kTRUE);
  fTracksKink   = new TObjArray(20); fTracksKink->SetOwner(kTRUE);
  fV0List       = new TObjArray(10); fV0List->SetOwner(kTRUE);
  fClusters     = new TObjArray(AliTRDgeometry::kNdet); fClusters->SetOwner(kTRUE);

  // define general monitor
  fContainer = new TObjArray(kNclasses); fContainer->SetOwner(kTRUE);
  TH1 *h=new TH1I("hStat", "Run statistics;Observable;Entries", Int_t(kNObjects), -0.5, Float_t(kNObjects)-0.5);
  TAxis *ax(h->GetXaxis());
  ax->SetBinLabel(Int_t(kTracksESD) + 1, "ESD");
  ax->SetBinLabel(Int_t(kTracksMC) + 1, "MC");
  ax->SetBinLabel(Int_t(kV0) + 1, "V0");
  ax->SetBinLabel(Int_t(kTPC) + 1, "TPC");
  ax->SetBinLabel(Int_t(kITS) + 1, "ITS");
  ax->SetBinLabel(Int_t(kTRDin) + 1, "TRDin");
  ax->SetBinLabel(Int_t(kTRDout) + 1, "TRDout");
  ax->SetBinLabel(Int_t(kBarrel) + 1, "Barrel");
  ax->SetBinLabel(Int_t(kBarrelMC) + 1, "BarrelMC");
  ax->SetBinLabel(Int_t(kSA) + 1, "SA");
  ax->SetBinLabel(Int_t(kSAMC) + 1, "SAMC");
  ax->SetBinLabel(Int_t(kKink) + 1, "Kink");
  ax->SetBinLabel(Int_t(kKinkMC) + 1, "KinkMC");
  ax->SetBinLabel(Int_t(kBarrelFriend) + 1, "BFriend");
  ax->SetBinLabel(Int_t(kSAFriend) + 1, "SFriend");
  fContainer->AddAt(h, kStatTrk);
  h=new TH1I("hEv", "Run statistics;Event Class;Entries", 4, -0.5, 3.5);
  ax = h->GetXaxis();
  ax->SetBinLabel(1, "Low");
  ax->SetBinLabel(2, "High");
  ax->SetBinLabel(3, "Cosmic");
  ax->SetBinLabel(4, "Calib");
  fContainer->AddAt(h, kEvType);
  TH2I* h2=new TH2I("hBCtrack", "Track Statistics;Fill Bunch;TOF BC;Entries", 3500, -0.5, 3499.5, 31, -10.5, 20.5);
  fContainer->AddAt(h2, kBC);
  fContainer->AddAt(new AliTRDtriggerInfo(), kTrigger);
  TObjArray *chmb = new TObjArray(AliTRDgeometry::kNdet);
  chmb->SetName("Chambers Status"); chmb->SetOwner(kTRUE);
  fContainer->AddAt(chmb, kChmb);

  PostData(AliTRDpwgppHelper::kTracksBarrel, fTracksBarrel);
  PostData(AliTRDpwgppHelper::kTracksITS, fTracksITS);
  PostData(AliTRDpwgppHelper::kTracksSA,     fTracksSA);
  PostData(AliTRDpwgppHelper::kTracksKink,   fTracksKink);
  PostData(AliTRDpwgppHelper::kEventInfo,    fEventInfo);
  PostData(AliTRDpwgppHelper::kV0List,       fV0List);
  PostData(AliTRDpwgppHelper::kClusters,     fClusters);
  PostData(AliTRDpwgppHelper::kMonitor,      fContainer);
}

//____________________________________________________________________
Bool_t AliTRDinfoGen::Load(const Char_t *file, const Char_t *dir, const Char_t *name)
{
// Load data from performance file

  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(!gFile->cd(dir)){
    AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
    gFile->Close();
    return kFALSE;
  }
  const Char_t *tn=(name ? name : GetName());
  if(!(fContainer = (TObjArray*)gDirectory->Get(tn))){
    AliWarning(Form("Missing histogram container %s.", tn));
    gFile->Close();
    return kFALSE;
  }
  gFile->Close();
  return kTRUE;
}

//____________________________________________________________________
void AliTRDinfoGen::UserExec(Option_t *){
  //
  // Run the Analysis
  //

  fTracksBarrel->Delete();
  fTracksITS->Delete();
  fTracksSA->Delete();
  fTracksKink->Delete();
  fV0List->Delete();
  fClusters->Delete();
  fEventInfo->Delete("");

  fESDev = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESDev){
    AliError("Failed retrieving ESD event");
    return;
  }
  // WARNING
  // This part may conflict with other detectors !!
  if(!IsInitOCDB()){
    AliCDBEntry* obj(NULL);
    AliCDBManager* ocdb = AliCDBManager::Instance();
    if(ocdb->IsDefaultStorageSet()){
      AliInfo("OCDB :: initialized externally.");
    } else {
      AliInfo("OCDB :: initializing locally ...");
      // prepare OCDB access
      ocdb->SetDefaultStorage(fOCDB.Data());
      //ocdb->SetSpecificStorage("TRD/Align/Data", "local:///home/niham/abercuci/local");
      ocdb->SetRun(fESDev->GetRunNumber());
      // create geo manager
      if(!(obj = ocdb->Get(AliCDBPath("GRP", "Geometry", "Data")))){
        AliError("GEOMETRY failed initialization.");
      } else {
        AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
        AliGeomManager::GetNalignable("TRD");
        AliGeomManager::GetNalignable("TPC");
        AliGeomManager::GetNalignable("ITS");
        AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD");
      }
      //init magnetic field
      if(!TGeoGlobalMagField::Instance()->IsLocked() &&
        !fESDev->InitMagneticField()){
        AliError("MAGNETIC FIELD failed initialization.");
      }
      // set no of time bins
      AliTRDtrackerV1::SetNTimeBins(AliTRDcalibDB::Instance()->GetNumberOfTimeBinsDCS());
    }
    AliInfo(Form("OCDB :  Loc[%s] Run[%d] TB[%d]", ocdb->GetDefaultStorage()->GetURI().Data(), ocdb->GetRun(), AliTRDtrackerV1::GetNTimeBins()));

    // load misalignment
    fgGeo = new AliTRDgeometry;
    fgGeo->CreateClusterMatrixArray();
    MakeChambers();
    // load reco param list from OCDB
    AliInfo("Initializing TRD reco params ...");
    fgReconstructor = new AliTRDReconstructor();
    if(!(obj = ocdb->Get(AliCDBPath("TRD", "Calib", "RecoParam")))){
      AliError("RECO PARAM failed initialization.");
    } else {
      obj->PrintMetaData();
      fRecos = (TObjArray*)obj->GetObject();
    }
    SetInitOCDB();
  }

  AliDebug(2, "+++++++++++++++++++++++++++++++++++++++");
  AliDebug(2, Form("Analyse Ev[%d]", fESDev->GetEventNumberInFile()));

  // set reco param valid for this event
  TH1 *h = (TH1I*)fContainer->At(kEvType);
  if(!fRecos){
    fgReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
    h->Fill(0);
  } else {
    for(Int_t ireco(0); ireco<fRecos->GetEntriesFast(); ireco++){
      AliTRDrecoParam *reco((AliTRDrecoParam*)fRecos->At(ireco));
      Int_t es(reco->GetEventSpecie());
      if(!(es&fESDev->GetEventSpecie())) continue;
      fgReconstructor->SetRecoParam(reco);
      if(AliLog::GetDebugLevel("PWGPP/TRD", "AliTRDinfoGen")>5) reco->Dump();
      TString s;
      if(es&AliRecoParam::kLowMult){ s="LowMult"; h->Fill(0);}
      else if(es&AliRecoParam::kHighMult){ s="HighMult"; h->Fill(1);}
      else if(es&AliRecoParam::kCosmic){ s="Cosmic"; h->Fill(2);}
      else if(es&AliRecoParam::kCalib){ s="Calib"; h->Fill(3);}
      else s="Unknown";
      AliDebug(2, Form("  Using reco param \"%s\" for event %d.", s.Data(), fESDev->GetEventNumberInFile()));
      break;
    }
  }

  // link MC if available
  fMCev = MCEvent();
  
  // trigger monitor
  AliTRDtriggerInfo *ti = (AliTRDtriggerInfo*)fContainer->At(kTrigger);
  TObjArray *evTriggers = fESDev->GetFiredTriggerClasses().Tokenize(" ");
  //printf("Ev[%03d] Triggers[%s]\n", fESDev->GetEventNumberInFile(), fESDev->GetFiredTriggerClasses().Data());
  for(Int_t iet(evTriggers->GetEntriesFast()); iet--;) ti->Add((*evTriggers)[iet]->GetName());
  evTriggers->Delete(); delete evTriggers;

  // event selection based on vertex cuts and trigger
  if(UseLocalEvSelection() && !fEventCut->IsSelected(fESDev, IsCollision())){
    AliDebug(2, "Event failed selection on vertex and trigger");
    return;
  }

  if(!fESDfriend){
    AliError("Failed retrieving ESD friend event");
    return;
  }
  if(HasMCdata() && !fMCev){
    AliError("Failed retrieving MC event");
    return;
  }

  new(fEventInfo)AliTRDeventInfo(fESDev->GetHeader(), const_cast<AliESDRun *>(fESDev->GetESDRun()));
  // Determine centrality
  // Author: Ionut Arsene <I.C.Arsene@gsi.de>
  TString beamtype = fESDev->GetESDRun()->GetBeamType();
  if(beamtype.Contains("Pb-Pb") || beamtype.Contains("A-A")){
    const AliMultiplicity *mult = fESDev->GetMultiplicity();
    fEventInfo->SetMultiplicity(mult?mult->GetNumberOfTracklets():0);
    const AliCentrality *cent = fESDev->GetCentrality();
    // centrality for different options V0 = "V0M", ITS = "TKL" etc
    fEventInfo->SetCentrality(cent?cent->GetCentralityPercentile("TKL"):-1.);
    AliDebug(2, Form("  Beam Type[%s] Mult[%d/%4d] Cent[%d/%5.2f]",
        beamtype.Data(),
        fEventInfo->GetMultiplicity(), mult?mult->GetNumberOfTracklets():0,
        fEventInfo->GetCentrality(), cent?cent->GetCentralityPercentile("TKL"):-1.));
  } else {
    fEventInfo->SetMultiplicity(0);
    fEventInfo->SetCentrality(-1.);
    AliDebug(2, Form("  Beam Type[%s]", beamtype.Data()));
  }
  UShort_t evBC(fESDev->GetBunchCrossNumber());

  // electron identifier from conversions
  if(!fV0Identifier) fV0Identifier = new AliESDv0KineCuts();
  fV0Identifier->SetEvent(fESDev);

  Bool_t *trackMap(NULL);
  AliStack * mStack(NULL);
  Int_t nTracksMC = HasMCdata() ? fMCev->GetNumberOfTracks() : 0, nTracksESD = fESDev->GetNumberOfTracks();
  if(HasMCdata()){
    mStack = fMCev->Stack();
    if(!mStack){
      AliError("Failed retrieving MC Stack");
      return;
    }
    trackMap = new Bool_t[nTracksMC];
    memset(trackMap, 0, sizeof(Bool_t) * nTracksMC);
  }
  
  Double32_t dedx[100]; Int_t nSlices(0);
  Int_t nTRDout(0), nTRDin(0), nTPC(0), nITS(0)
//        ,nclsTrklt
       ,nBarrel(0), nBarrelITS(0), nSA(0), nKink(0)
       ,nBarrelFriend(0), nBarrelITSFriend(0), nSAFriend(0)
       ,nBarrelMC(0), nSAMC(0), nKinkMC(0);
  AliESDtrack *esdTrack = NULL;
  AliESDfriendTrack *esdFriendTrack = NULL;
  TObject *calObject = NULL;
  AliTRDtrackV1 *track = NULL;
  AliTRDseedV1 *tracklet = NULL;
  AliTRDcluster *cl = NULL;


  // LOOP 0 - over ESD v0s
  Float_t bField(fESDev->GetMagneticField());
  AliESDv0 *v0(NULL);
  Int_t v0pid[AliPID::kSPECIES];
  for(Int_t iv0(0); iv0<fESDev->GetNumberOfV0s(); iv0++){
    // Take only V0s from the On-the-fly v0 finder
    if(!(v0 = fESDev->GetV0(iv0))) continue;
    if(!v0->GetOnFlyStatus()) continue;
    // register v0
    if(fV0Cut) new(fV0Info) AliTRDv0Info(*fV0Cut);
    else new(fV0Info) AliTRDv0Info();
    fV0Info->SetMagField(bField);
    fV0Info->SetV0tracks(fESDev->GetTrack(v0->GetPindex()), fESDev->GetTrack(v0->GetNindex()));
    fV0Info->SetV0Info(v0);
    // tag conversion electrons
    Int_t pdgV0, pdgP, pdgN;
    if(fV0Identifier->ProcessV0(v0, pdgV0, pdgP, pdgN)) {
      switch(pdgV0){
      case kGamma: fV0Info->SetDecay(AliTRDv0Info::kGamma); break;
      case kK0Short: fV0Info->SetDecay(AliTRDv0Info::kK0s); break;
      case kLambda0: fV0Info->SetDecay(AliTRDv0Info::kLambda); break;
      case kLambda0Bar: fV0Info->SetDecay(AliTRDv0Info::kAntiLambda); break;
      default: AliDebug(1, Form("V0[%+4d] -> +[%+4d] -[%+4d]. Decay not mapped.", pdgV0, pdgP, pdgN));
      }
    }
    fV0List->Add(new AliTRDv0Info(*fV0Info));//  kFOUND=kFALSE;
  }


  // LOOP 1 - over ESD tracks
  AliTRDv0Info *v0info=NULL;
  for(Int_t itrk = 0; itrk < nTracksESD; itrk++){
    new(fTrackInfo) AliTRDtrackInfo();
    esdTrack = fESDev->GetTrack(itrk);
    AliDebug(3, Form("\n%3d ITS[%d] TPC[%d] TRD[%d] TOF-BC[%d]\n", itrk, esdTrack->GetNcls(0), esdTrack->GetNcls(1), esdTrack->GetNcls(2), esdTrack->GetTOFBunchCrossing()));
    if(esdTrack->GetStatus()&AliESDtrack::kITSout) nITS++;
    if(esdTrack->GetStatus()&AliESDtrack::kTPCout) nTPC++;
    if(esdTrack->GetStatus()&AliESDtrack::kTRDout) nTRDout++;
    if(esdTrack->GetStatus()&AliESDtrack::kTRDin) nTRDin++;

/*    Int_t ns(esdTrack->GetNumberOfTRDslices());
    printf("  %3d ITS[%c] TPC[%c] TRDin[%c] TRDout[%c] TRDStop[%c] ns[%d]\n", itrk,
      (esdTrack->GetStatus()&AliESDtrack::kITSout)?'y':'n',
      (esdTrack->GetStatus()&AliESDtrack::kTPCout)?'y':'n',
      (esdTrack->GetStatus()&AliESDtrack::kTRDin)?'y':'n',
      (esdTrack->GetStatus()&AliESDtrack::kTRDout)?'y':'n',
      (esdTrack->GetStatus()&AliESDtrack::kTRDStop)?'y':'n', ns);
    if(ns){
      for(Int_t ipl(0); ipl<AliTRDgeometry::kNlayer; ipl++){
        Double_t sp, p(esdTrack->GetTRDmomentum(ipl, &sp));
        printf("    [%d] p[%6.3f+-%6.3f] dEdx={", ipl, p, sp);
        for(Int_t is(0); is<8; is++) printf("%7.2f ", esdTrack->GetTRDslice(ipl, is)); printf("}\n");
      }
    }
*/
    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    Double_t xyz[3];
    if(op){
      op->GetXYZ(xyz);
      op->Global2LocalPosition(xyz, op->GetAlpha());
      AliDebug(3, Form("op @ X[%7.3f]\n", xyz[0]));
    }

    // read MC info
    Int_t label = -1; UInt_t alab=UINT_MAX;
    if(HasMCdata()){
      label = esdTrack->GetLabel(); 
      alab = TMath::Abs(label);
      // register the track
      if(alab < UInt_t(nTracksMC)){ 
        trackMap[alab] = kTRUE; 
      } else { 
        AliError(Form("MC label[%d] outside scope for Ev[%d] Trk[%d].", label, (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), itrk));
        continue; 
      }
      AliMCParticle *mcParticle(NULL);
      if(!(mcParticle = (AliMCParticle*) fMCev->GetTrack(alab))){
        AliError(Form("MC particle label[%d] missing for Ev[%d] Trk[%d].", label, (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), itrk));
        continue;
      }
      fTrackInfo->SetMC();
      fTrackInfo->SetMCeta(mcParticle->Eta());
      fTrackInfo->SetMCphi(mcParticle->Phi());
      fTrackInfo->SetMCpt(mcParticle->Pt());
      fTrackInfo->SetPDG(mcParticle->Particle()->GetPdgCode());
      fTrackInfo->SetPrimary(mcParticle->Particle()->IsPrimary());
      fTrackInfo->SetLabel(label);
      fTrackInfo->SetTRDlabel(esdTrack->GetTRDLabel());
      AliTrackReference *ref(NULL);
      for(Int_t iref(0); iref<mcParticle->GetNumberOfTrackReferences(); iref++){
        if(!(ref = mcParticle->GetTrackReference(iref))) continue;
        if(ref->DetectorId() != AliTrackReference::kTRD) continue;
        AliDebug(4, Form("  TRD trackRef[%2d] @ r[%7.3f] [REC]", iref, ref->LocalX()));
        fTrackInfo->AddTrackRef(ref);
      }
      AliDebug(3, Form("Lab[%4d] pdg[%+4d] NtrackRefs[%d(%d)]", alab, mcParticle->Particle()->GetPdgCode(), fTrackInfo->GetNTrackRefs(), mcParticle->GetNumberOfTrackReferences()));
    }

    // copy some relevant info to TRD track info
    fTrackInfo->SetStatus(esdTrack->GetStatus());
    fTrackInfo->SetTrackId(esdTrack->GetID());
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    fTrackInfo->SetESDpid(p);
    fTrackInfo->SetESDpidQuality(esdTrack->GetTRDntrackletsPID());
    fTrackInfo->SetESDeta(esdTrack->Eta());
    Double_t loc[3];
    if(esdTrack->GetXYZAt(298., fESDev->GetMagneticField(), loc)) fTrackInfo->SetESDphi(TMath::ATan2(loc[1], loc[0]));
    fTrackInfo->SetESDpt(esdTrack->Pt());
    if(!nSlices) nSlices = esdTrack->GetNumberOfTRDslices();
    memset(dedx, 0, 100*sizeof(Double32_t));
    Int_t in(0);
    for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++)
      for(Int_t is=0; is<nSlices; is++) 
        dedx[in++]=esdTrack->GetTRDslice(il, is);
    for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++) dedx[in++]=esdTrack->GetTRDmomentum(il);
    fTrackInfo->SetSlices(in, dedx);
    fTrackInfo->SetNumberOfClustersRefit(esdTrack->GetNcls(2));
    // some other Informations which we may wish to store in order to find problematic cases
    fTrackInfo->SetKinkIndex(esdTrack->GetKinkIndex(0));
    fTrackInfo->SetTPCncls(static_cast<UShort_t>(esdTrack->GetNcls(1)));
    fTrackInfo->SetTPCdedx(esdTrack->GetTPCsignal());
    Float_t tofTime = esdTrack->GetTOFsignal() - fESDev->GetT0TOF(0);
    fTrackInfo->SetTOFbeta(tofTime>0.?((esdTrack->GetIntegratedLength()/(tofTime*TMath::C()))*10e9):-999.);
    fTrackInfo->SetTOFbc(esdTrack->GetTOFBunchCrossing()==AliVTrack::kTOFBCNA?0:esdTrack->GetTOFBunchCrossing());
//    nclsTrklt = 0;
  
    // set V0pid info
    //printf("%4d Looking for V0s...\n" , fTrackInfo->GetTrackId());
    for(Int_t iv(0); iv<fV0List->GetEntriesFast(); iv++){
      if(!(v0info = (AliTRDv0Info*)fV0List->At(iv))) continue;
      if(!v0info->GetV0Daughter(1) && !v0info->GetV0Daughter(-1)) continue;
      if(!v0info->HasTrack(fTrackInfo)) continue;
      //v0info->Print();
      memset(v0pid, 0, AliPID::kSPECIES*sizeof(Int_t));
      fTrackInfo->SetV0();
      for(Int_t is=AliPID::kSPECIES; is--;) v0pid[is] = v0info->GetPID(is, fTrackInfo); fTrackInfo->SetV0pid(v0pid);
      if(v0info->IsDecay(AliTRDv0Info::kGamma)) fTrackInfo->SetElectron();
      else if(v0info->IsDecay(AliTRDv0Info::kK0s)) fTrackInfo->SetPion();
      else if(v0info->IsDecay(AliTRDv0Info::kLambda)) esdTrack->Charge()>0?fTrackInfo->SetProton():fTrackInfo->SetPion();
      else if(v0info->IsDecay(AliTRDv0Info::kAntiLambda)) esdTrack->Charge()<0?fTrackInfo->SetProton():fTrackInfo->SetPion();
      
      //TODO one track can be attached to more than one v0. Ideally one would need a list of v0 attached to the track info
    }

    // read track REC info
    if((esdFriendTrack = (fESDfriend->GetNumberOfTracks() > itrk) ? fESDfriend->GetTrack(itrk): NULL)) {
      fTrackInfo->SetTPCoutParam(esdFriendTrack->GetTPCOut());
      fTrackInfo->SetITSoutParam(esdFriendTrack->GetITSOut());
      const AliTrackPointArray *tps(NULL);
      if((tps=esdFriendTrack->GetTrackPointArray()) && HasTrackPoints()) fTrackInfo->SetTrackPointArray(tps);
      Int_t icalib = 0;
      while((calObject = esdFriendTrack->GetCalibObject(icalib++))){
        if(strcmp(calObject->IsA()->GetName(),"AliTRDtrackV1") != 0) continue; // Look for the TRDtrack
        if(!(track = dynamic_cast<AliTRDtrackV1*>(calObject))) break;
        AliDebug(4, Form("TRD track OK"));
        // Set the clusters to unused
        for(Int_t ipl = 0; ipl < AliTRDgeometry::kNlayer; ipl++){
          if(!(tracklet = track->GetTracklet(ipl))) continue;
          tracklet->ResetClusterIter();
          while((cl = tracklet->NextCluster())) cl->Use(0);
        }
        fTrackInfo->SetTrack(track);
        ((TH2I*)fContainer->At(kBC))->Fill(evBC, esdTrack->GetTOFBunchCrossing());
        break;
      }
      AliDebug(3, Form("Ntracklets[%d]", fTrackInfo->GetNTracklets()));
    } else AliDebug(3, Form("No Friends for trk[%3d] Ntrk[%3d]", itrk, fESDfriend->GetNumberOfTracks()));
    if(op) fTrackInfo->SetOuterParam(op);

    if(DebugLevel() >= 1){
      AliTRDtrackInfo info(*fTrackInfo);
      (*DebugStream()) << "trackInfo"
      << "TrackInfo.=" << &info
      << "\n";
      info.Delete("");
    }

    ULong_t status(esdTrack->GetStatus());
    if((status&AliESDtrack::kTPCout)){  // TPC prolongation
      if(!esdTrack->GetKinkIndex(0)){ // Barrel  Track Selection
        Bool_t selected(kTRUE);
        if(UseLocalTrkSelection()){
          if(esdTrack->Pt() < fgkPt){ 
            AliDebug(3, Form("Reject TPC Trk[%3d] Ev[%4d] Pt[%5.2f]", itrk, fESDev->GetEventNumberInFile(), esdTrack->Pt()));
            selected = kFALSE;
          }
          if(selected && TMath::Abs(esdTrack->Eta()) > fgkEta){
            AliDebug(3, Form("Reject TPC Trk[%3d] Ev[%4d] Eta[%5.2f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(esdTrack->Eta())));
            selected = kFALSE;
          }
          if(selected && esdTrack->GetTPCNcls() < fgkNclTPC){ 
            AliDebug(3, Form("Reject TPC Trk[%3d] Ev[%4d] NclTPC[%d]", itrk, fESDev->GetEventNumberInFile(), esdTrack->GetTPCNcls()));
            selected = kFALSE;
          }
          if(selected && !(esdTrack->GetStatus()&AliESDtrack::kITSrefit)){ //SPD refit flag (Ionut)
            AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] ITSrefit", itrk, fESDev->GetEventNumberInFile()));
            selected = kFALSE;
          }
          UChar_t itsMap = esdTrack->GetITSClusterMap();
          Bool_t firstSPDlayerHit = (itsMap & (UChar_t(1)<<0)),
                 secondSPDlayerHit = (itsMap & (UChar_t(1)<<1));
          if(selected && !(firstSPDlayerHit || secondSPDlayerHit)){ // request at least 1 SPD cluster (Ionut)
            AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] no SPD cluster", itrk, fESDev->GetEventNumberInFile()));
            selected = kFALSE;
          }

          Float_t par[2], cov[3];
          esdTrack->GetImpactParameters(par, cov);
          if(IsCollision()){ // cuts on DCA
            if(selected && TMath::Abs(par[0]) > fgkTrkDCAxy){ 
              AliDebug(3, Form("Reject TPC Trk[%3d] Ev[%4d] DCAxy[%f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(par[0])));
              selected = kFALSE;
            }
            if(selected && TMath::Abs(par[1]) > fgkTrkDCAz){ 
              AliDebug(3, Form("Reject TPC Trk[%3d] Ev[%4d] DCAz[%f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(par[1])));
              selected = kFALSE;
            }
          } else if(selected && fMCev && !fMCev->IsPhysicalPrimary(alab)){;
            AliDebug(3, Form("Reject TPC Trk[%3d] Ev[%4d] Primary", itrk, fESDev->GetEventNumberInFile()));
            selected = kFALSE;
          }
        }
        if(fTrackCut && !fTrackCut->IsSelected(esdTrack)) selected = kFALSE;
        if(selected){ 
          fTracksBarrel->Add(new AliTRDtrackInfo(*fTrackInfo));
          nBarrel++;
          if(fTrackInfo->GetTrack()) nBarrelFriend++;
        }
      } else {
        fTracksKink->Add(new AliTRDtrackInfo(*fTrackInfo));
        nKink++;
      }
    } else if((status&AliESDtrack::kITSout)) { // ITS prolongation 
      Bool_t selected(kTRUE);
      if(UseLocalTrkSelection()){
        if(esdTrack->Pt() < fgkPt){
          AliDebug(3, Form("Reject ITS Trk[%3d] Ev[%4d] Pt[%5.2f]", itrk, fESDev->GetEventNumberInFile(), esdTrack->Pt()));
          selected = kFALSE;
        }
        if(selected && TMath::Abs(esdTrack->Eta()) > fgkEta){
          AliDebug(3, Form("Reject ITS Trk[%3d] Ev[%4d] Eta[%5.2f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(esdTrack->Eta())));
          selected = kFALSE;
        }
        Float_t par[2], cov[3];
        esdTrack->GetImpactParameters(par, cov);
        if(IsCollision()){ // cuts on DCA
          if(selected && TMath::Abs(par[0]) > fgkTrkDCAxy){
            AliDebug(3, Form("Reject ITS Trk[%3d] Ev[%4d] DCAxy[%f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(par[0])));
            selected = kFALSE;
          }
          if(selected && TMath::Abs(par[1]) > fgkTrkDCAz){
            AliDebug(3, Form("Reject ITS Trk[%3d] Ev[%4d] DCAz[%f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(par[1])));
            selected = kFALSE;
          }
        } else if(selected && fMCev && !fMCev->IsPhysicalPrimary(alab)){;
          AliDebug(3, Form("Reject ITS Trk[%3d] Ev[%4d] Primary", itrk, fESDev->GetEventNumberInFile()));
          selected = kFALSE;
        }
      }
      if(fTrackCut && !fTrackCut->IsSelected(esdTrack)) selected = kFALSE;
      if(selected){
        fTracksITS->Add(new AliTRDtrackInfo(*fTrackInfo));
        nBarrelITS++;
        if(fTrackInfo->GetTrack()) nBarrelITSFriend++;
      }
    } else if((status&AliESDtrack::kTRDout) && !(status&AliESDtrack::kTRDin)){ // TRD SA tracking
      fTracksSA->Add(new AliTRDtrackInfo(*fTrackInfo));
      nSA++;
      if(fTrackInfo->GetTrack()) nSAFriend++;
    }
    fTrackInfo->Delete("");
  }
  // read clusters REC info
  TTree * treeR(NULL); AliESDInputHandlerRP *esdRPhandler(NULL);
  if(fInputHandler) esdRPhandler = dynamic_cast<AliESDInputHandlerRP*>(fInputHandler);
  if(esdRPhandler){
    if((treeR = esdRPhandler->GetTreeR("TRD"))) {
      TObjArray *recPoints(NULL);
      if((treeR->GetBranch("TRDcluster"))){
        treeR->SetBranchAddress("TRDcluster", &recPoints);
        for(Int_t idet(0); idet<treeR->GetEntries(); idet++){
          treeR->GetEntry(idet);
          if(!recPoints->GetEntries()){
            AliDebug(1, Form("Missing entry %d from TreeR", idet));
            continue;
          }
          AliTRDcluster *c = (AliTRDcluster*)(*recPoints)[0];
          if(!c){
            AliDebug(1, Form("Missing first cluster in entry %d from TreeR", idet));
            continue;
          }
          fClusters->AddAt(recPoints->Clone(Form("%03d", c->GetDetector())), c->GetDetector());
        }
      } else AliDebug(3, "No TRDcluster branch");
    } else AliDebug(3, "No RecPoints");
  } else AliDebug(3, "No AliESDInputHandlerRP");


  // LOOP 2 - over MC tracks which are passing TRD where the track is not reconstructed
  if(HasMCdata()){
    AliDebug(10, "Output of the MC track map:");
    for(Int_t itk = 0; itk < nTracksMC;  itk++) AliDebug(10, Form("trackMap[%d] = %s", itk, trackMap[itk] == kTRUE ? "TRUE" : "kFALSE"));
  
    AliTrackReference *ref(NULL);
    for(Int_t itk = 0; itk < nTracksMC; itk++){
      if(trackMap[itk]) continue;
      AliMCParticle *mcParticle =  (AliMCParticle*) fMCev->GetTrack(TMath::Abs(itk));

      Int_t nRefsTRD(0);
      for(Int_t iref(0); iref<mcParticle->GetNumberOfTrackReferences(); iref++){ // count TRD TR
        if(!(ref = mcParticle->GetTrackReference(iref))) continue;
        if(ref->DetectorId() != AliTrackReference::kTRD) continue;
        if(!nRefsTRD){ // build track info for this pure MC track
          new(fTrackInfo) AliTRDtrackInfo();
          fTrackInfo->SetMC();
          fTrackInfo->SetPDG(mcParticle->Particle()->GetPdgCode());
        }
        AliDebug(4, Form("  TRD trackRef[%2d] @ r[%7.3f] [MC]", iref, ref->LocalX()));
        fTrackInfo->AddTrackRef(ref);
      }
      // In this stage we at least require 1 hit inside TRD. What will be done with this tracks is a task for the
      // analysis job
      if(!nRefsTRD) continue;
      AliDebug(3, Form("Lab[%4d] pdg[%+4d] NtrackRefs[%d(%d)]", itk, mcParticle->Particle()->GetPdgCode(), fTrackInfo->GetNTrackRefs(), mcParticle->GetNumberOfTrackReferences()));
      fTrackInfo->SetPrimary(mcParticle->Particle()->IsPrimary());
      fTrackInfo->SetLabel(itk);
      if(DebugLevel() >= 1){
        AliTRDtrackInfo info(*fTrackInfo);
        (*DebugStream()) << "trackInfo"
        << "TrackInfo.=" << &info
        << "\n";
        info.Delete("");
      }
      AliDebug(3, Form("Add MC track @ label[%d] nTRDrefs[%d].", itk, nRefsTRD));
      // check where the track starts
      ref = mcParticle->GetTrackReference(0);
      if(ref->LocalX() < fgkITS){ 
        fTracksBarrel->Add(new AliTRDtrackInfo(*fTrackInfo));
        //fTracksITS->Add(new AliTRDtrackInfo(*fTrackInfo));
        nBarrelMC++;
      } else if(ref->LocalX() < fgkTPC) {
        fTracksKink->Add(new AliTRDtrackInfo(*fTrackInfo));
        nKinkMC++;
      } else if(nRefsTRD>6){
        fTracksSA->Add(new AliTRDtrackInfo(*fTrackInfo));
        nSAMC++;
      }
      fTrackInfo->Delete("");
    }
    delete[] trackMap;
  }
  AliDebug(1, Form(
    "\nEv[%3d] Tracks: ESD[%d] MC[%d] V0[%d]\n"
    "        TPCout[%d] ITSout[%d] TRDin[%d] TRDout[%d]\n"
    "        Barrel[%3d+%3d=%3d] ITS[%3d=%3d] SA[%2d+%2d=%2d] Kink[%2d+%2d=%2d]"
    ,(Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nTracksESD, nTracksMC, fV0List->GetEntries()
    , nTPC, nITS, nTRDin, nTRDout
    ,nBarrel, nBarrelMC, fTracksBarrel->GetEntries()
    ,nBarrelITS, fTracksITS->GetEntries()
    ,nSA, nSAMC, fTracksSA->GetEntries()
    ,nKink, nKinkMC, fTracksKink->GetEntries()
  ));
  // save track statistics
  h = (TH1I*)fContainer->At(kStatTrk);
  h->Fill(Float_t(kTracksESD), nTracksESD);
  h->Fill(Float_t(kTracksMC), nTracksMC);
  h->Fill(Float_t(kV0), fV0List->GetEntries());
  h->Fill(Float_t(kTPC), nTPC);
  h->Fill(Float_t(kITS), nITS);
  h->Fill(Float_t(kTRDin), nTRDin);
  h->Fill(Float_t(kTRDout), nTRDout);
  h->Fill(Float_t(kBarrel), nBarrel);
  h->Fill(Float_t(kBarrelMC), nBarrelMC);
  h->Fill(Float_t(kSA), nSA);
  h->Fill(Float_t(kSAMC), nSAMC);
  h->Fill(Float_t(kKink), nKink);
  h->Fill(Float_t(kKinkMC), nKinkMC);
  h->Fill(Float_t(kBarrelFriend), nBarrelFriend);
  h->Fill(Float_t(kSAFriend), nSAFriend);
}


//____________________________________________________________________
void AliTRDinfoGen::MakeChambers()
{
// Build chamber position and status
  if(!fContainer){
    AliError("Missing container");
    return;
  }
  AliTRDcalibDB *calib(AliTRDcalibDB::Instance());
  if(!calib){
    AliError("No access to calibration data");
    return;
  }
  TObjArray *chmb = (TObjArray*)fContainer->At(kChmb);
  Int_t stat(0);
  Double_t alpha(0.), cs(-2.), sn(0.), pos[4];
  for(Int_t isec(0); isec<AliTRDgeometry::kNsector; isec++){
    alpha = (0.5+isec)*AliTRDgeometry::GetAlpha();
    cs    = TMath::Cos(alpha);
    sn    = TMath::Sin(alpha);

    for(Int_t istk(0); istk<AliTRDgeometry::kNstack; istk++){
      for(Int_t ilyr(0); ilyr<AliTRDgeometry::kNlayer; ilyr++){
        Int_t idet(AliTRDgeometry::GetDetector(ilyr, istk, isec));
        TGeoHMatrix *matrix(fgGeo->GetClusterMatrix(idet));
        if(!matrix){
          AliDebug(2, Form("Missing matrix for %03d [%02d_%d_%d]", idet, isec, istk, ilyr));
          continue;
        }
        AliDebug(2, Form("Read info for %03d [%02d_%d_%d]", idet, isec, istk, ilyr));
        AliTRDpadPlane *pp(fgGeo->GetPadPlane(ilyr, istk));
        Double_t zm(0.5 * (pp->GetRow0() + pp->GetRowEnd())),
                 loc0[] = {AliTRDgeometry::AnodePos(), pp->GetCol0(), zm-pp->GetRow0()},
                 loc1[] = {AliTRDgeometry::AnodePos(), pp->GetColEnd(), zm-pp->GetRowEnd()},
                 glb[3] = {1,1,1};
        matrix->LocalToMaster(loc0, glb);
        Float_t phi = TMath::ATan2(glb[0]*sn + glb[1]*cs, glb[0]*cs - glb[1]*sn),
                tgl = glb[2]/glb[0]/TMath::Sqrt(1.+glb[1]*glb[1]/glb[0]/glb[0]),
                eta = -TMath::Log(TMath::Tan(0.5 *  (0.5*TMath::Pi() - TMath::ATan(tgl))));
        pos[0] = eta; pos[1] = phi;
        matrix->LocalToMaster(loc1, glb);
        phi = TMath::ATan2(glb[0]*sn + glb[1]*cs, glb[0]*cs - glb[1]*sn);
        tgl = glb[2]/glb[0]/TMath::Sqrt(1.+glb[1]*glb[1]/glb[0]/glb[0]);
        eta = -TMath::Log(TMath::Tan(0.5 *  (0.5*TMath::Pi() - TMath::ATan(tgl))));
        pos[2] = eta; pos[3] = phi;
        stat = 0;
        if(calib->IsChamberGood(idet)){
          if(calib->IsHalfChamberNoData(idet, 0)) stat += 2;
          if(calib->IsHalfChamberNoData(idet, 1)) stat += 3;
        } else stat = 1;
        chmb->Add(new AliTRDchmbInfo(idet, stat, pos));
      }
    }
  }
}

//____________________________________________________________________
void AliTRDinfoGen::MakeSummary()
{
// Build summary plots
  if(!fContainer){
    AliError("Missing results");
    return;
  }
  TH1 *h1(NULL); TVirtualPad *p(NULL); TCanvas *cOut(NULL);

  const Int_t nx(1024), ny(1024);
  cOut = new TCanvas("infoGenSummary", "Run Statistics", nx, ny);
  cOut->Divide(2,2, 1.e-5, 1.e-5);
  //=========
  p=cOut->cd(1);p->SetRightMargin(0.025);p->SetTopMargin(0.01);p->SetLogy();
  h1 = (TH1*)fContainer->At(kStatTrk);
  h1->SetBarOffset(0.06); h1->SetBarWidth(0.88); h1->SetFillColor(kGreen); h1->SetFillStyle(3001);
  h1->Draw("bar1");
  //=========
  p=cOut->cd(2);p->SetRightMargin(0.025);p->SetTopMargin(0.01);p->SetLogy();
  h1 = (TH1*)fContainer->At(kEvType);
  h1->SetBarOffset(0.04); h1->SetBarWidth(0.92);h1->SetFillColor(kGreen); h1->SetFillStyle(3001);
  h1->Draw("bar1");
  //=========
  p=cOut->cd(3);p->SetRightMargin(0.025);p->SetTopMargin(0.01);p->SetLogz();
  TH2 *h2((TH2*)fContainer->At(kBC));
  h1 = h2->ProjectionX("_px");
  Int_t n(0); Int_t bins[3500];
  for(Int_t ib(1); ib<=3500; ib++){
    if(h1->GetBinContent(ib) < 1) continue;
    bins[n++] = ib;
  }
  delete h1;

  TAxis *ay(h2->GetYaxis());
  TH2 *hs = new TH2I("hBC_Summary", Form("%s;%s;%s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()),
                     n, -0.5, n-0.5, ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
  hs->SetLineColor(kBlack);hs->SetLineWidth(1);
  hs->SetMarkerColor(kRed);
  TAxis *ax(hs->GetXaxis()); ax->CenterTitle(); ax->SetTitleOffset(1.4);
  for(Int_t ib(0); ib<n; ib++){
    ax->SetBinLabel(ib+1, Form("%d", bins[ib]));
    for(Int_t iy(1); iy<=ay->GetNbins(); iy++){
      hs->SetBinContent(ib+1, iy, h2->GetBinContent(bins[ib], iy));
    }
  }
  hs->Draw("textbox");

  //=========
  p=cOut->cd(4); p->SetRightMargin(0.0215);p->SetLeftMargin(0.414);//p->SetLogz();
  TObject *o = fContainer->At(kTrigger);
  if(o){
    if((h1 = dynamic_cast<TH1I*>(o))) {
      h1->GetXaxis()->SetTitleOffset(6.5); h1->GetXaxis()->CenterTitle();
      h1->SetFillStyle(3001);h1->SetFillColor(kGreen);
      h1->SetBarWidth(0.8);h1->SetBarOffset(0.1);
      ((TH1I*)o)->Draw("hbar2");
    } else o->Draw();
  }
  cOut->SaveAs(Form("%s.gif", cOut->GetName()));
}

//____________________________________________________________________
void AliTRDinfoGen::SetLocalEvSelection(const AliTRDeventCuts &ec)
{
// Set event cuts from outside
  if(!fEventCut) fEventCut = new AliTRDeventCuts(ec);
  else new(fEventCut) AliTRDeventCuts(ec);
  fEventCut->Print();
}

//____________________________________________________________________
void AliTRDinfoGen::SetLocalV0Selection(const AliTRDv0Info &v0)
{
// Set V0 cuts from outside

  if(!fV0Cut) fV0Cut = new AliTRDv0Info(v0);
  else new(fV0Cut) AliTRDv0Info(v0);
  fV0Cut->Print();
}

//____________________________________________________________________
TTreeSRedirector* AliTRDinfoGen::DebugStream()
{
// Manage debug stream for task
  if(!fDebugStream){
    TDirectory *savedir = gDirectory;
    fDebugStream = new TTreeSRedirector("TRD.DebugInfoGen.root");
    savedir->cd();
  }
  return fDebugStream;
}

//____________________________________________________________________
void AliTRDinfoGen::Terminate(Option_t* /*option*/)
{
// Process run information
  AliInfo("");
  if(!(fContainer = dynamic_cast<TObjArray *>(GetOutputData(AliTRDpwgppHelper::kMonitor)))) return;
  AliInfo(Form("fContainer(%p)", (void*)fContainer));

  AliTRDtriggerInfo* ti(NULL);
  if(UseLocalEvSelection()){
    if(!(ti = (AliTRDtriggerInfo*)fContainer->At(kTrigger))) return;
    for(Int_t ix(0); ix<ti->GetNTriggers(); ix++){
      if(fEventCut->CheckTrigger(ti->GetTrigger(ix))) ti->SetSelectTrigger(ix);
      //ax->SetBinLabel(ix, Form("#color[2]{%s}", ax->GetBinLabel(ix)));
    }
  }
}

