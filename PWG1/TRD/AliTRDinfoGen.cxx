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
#include <TH1S.h>
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
#include "AliMCEventHandler.h"

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDHeader.h"
#include "AliESDRun.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"
#include "AliMCParticle.h"
#include "AliMultiplicity.h"
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
#include "AliTRDgeometry.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTRDinfoGen.h"
#include "AliTRDpwg1Helper.h"
#include "info/AliTRDtrackInfo.h"
#include "info/AliTRDeventInfo.h"
#include "info/AliTRDv0Info.h"
#include "info/AliTRDeventCuts.h"

ClassImp(AliTRDinfoGen)

const Float_t AliTRDinfoGen::fgkITS = 100.; // to be checked
const Float_t AliTRDinfoGen::fgkTPC = 290.;
const Float_t AliTRDinfoGen::fgkTRD = 365.;

const Float_t AliTRDinfoGen::fgkEvVertexZ = 15.;
const Int_t   AliTRDinfoGen::fgkEvVertexN = 1;

const Float_t AliTRDinfoGen::fgkTrkDCAxy  = 3.;
const Float_t AliTRDinfoGen::fgkTrkDCAz   = 10.;
const Int_t   AliTRDinfoGen::fgkNclTPC    = 70;
const Float_t AliTRDinfoGen::fgkPt        = 0.2;
const Float_t AliTRDinfoGen::fgkEta       = 0.9;
AliTRDReconstructor* AliTRDinfoGen::fgReconstructor(NULL);
AliTRDgeometry* AliTRDinfoGen::fgGeo(NULL);
//____________________________________________________________________
AliTRDinfoGen::AliTRDinfoGen()
  :AliAnalysisTaskSE()
  ,fEvTrigger(NULL)
  ,fESDev(NULL)
  ,fMCev(NULL)
  ,fEventCut(NULL)
  ,fTrackCut(NULL)
  ,fV0Cut(NULL)
  ,fOCDB("local://$ALICE_ROOT/OCDB")
  ,fTrackInfo(NULL)
  ,fEventInfo(NULL)
  ,fV0Info(NULL)
  ,fTracksBarrel(NULL)
  ,fTracksSA(NULL)
  ,fTracksKink(NULL)
  ,fV0List(NULL)
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
  ,fEvTrigger(NULL)
  ,fESDev(NULL)
  ,fMCev(NULL)
  ,fEventCut(NULL)
  ,fTrackCut(NULL)
  ,fV0Cut(NULL)
  ,fOCDB("local://$ALICE_ROOT/OCDB")
  ,fTrackInfo(NULL)
  ,fEventInfo(NULL)
  ,fV0Info(NULL)
  ,fTracksBarrel(NULL)
  ,fTracksSA(NULL)
  ,fTracksKink(NULL)
  ,fV0List(NULL)
  ,fContainer(NULL)
  ,fRecos(NULL)
  ,fDebugStream(NULL)
{
  //
  // Default constructor
  //
  SetTitle("MC-REC TRD-track list generator");
  DefineOutput(AliTRDpwg1Helper::kTracksBarrel, TObjArray::Class());
  DefineOutput(AliTRDpwg1Helper::kTracksSA, TObjArray::Class());
  DefineOutput(AliTRDpwg1Helper::kTracksKink, TObjArray::Class());
  DefineOutput(AliTRDpwg1Helper::kEventInfo, AliTRDeventInfo::Class());
  DefineOutput(AliTRDpwg1Helper::kV0List, TObjArray::Class());
  DefineOutput(AliTRDpwg1Helper::kMonitor, TObjArray::Class()); // histogram list
}

//____________________________________________________________________
AliTRDinfoGen::~AliTRDinfoGen()
{
// Destructor
  if(fgGeo) delete fgGeo;
  if(fgReconstructor) delete fgReconstructor;
  if(fDebugStream) delete fDebugStream;
  if(fEvTrigger) delete fEvTrigger;
  if(fV0Cut) delete fV0Cut;
  if(fTrackCut) delete fTrackCut;
  if(fEventCut) delete fEventCut;
  if(fTrackInfo) delete fTrackInfo; fTrackInfo = NULL;
  if(fEventInfo) delete fEventInfo; fEventInfo = NULL;
  if(fV0Info) delete fV0Info; fV0Info = NULL;
  if(fTracksBarrel){ 
    fTracksBarrel->Delete(); delete fTracksBarrel;
    fTracksBarrel = NULL;
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
  if(fContainer){ 
    fContainer->Delete(); 
    delete fContainer;
    fContainer = NULL;
  }
}

//____________________________________________________________________
Bool_t AliTRDinfoGen::GetRefFigure(Int_t)
{
// General graphs for PWG1/TRD train
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
  fTracksSA = new TObjArray(20); fTracksSA->SetOwner(kTRUE);
  fTracksKink = new TObjArray(20); fTracksKink->SetOwner(kTRUE);
  fV0List = new TObjArray(10); fV0List->SetOwner(kTRUE);

  // define general monitor
  fContainer = new TObjArray(kNclasses); fContainer->SetOwner(kTRUE);
  TH1 *h=new TH1I("hStat", "Run statistics;Observable;Entries", Int_t(kNObjects), -0.5, Float_t(kNObjects)-0.5);
  TAxis *ax(h->GetXaxis());
  ax->SetBinLabel(Int_t(kTracksESD) + 1, "ESD");
  ax->SetBinLabel(Int_t(kTracksMC) + 1, "MC");
  ax->SetBinLabel(Int_t(kV0) + 1, "V0");
  ax->SetBinLabel(Int_t(kTPC) + 1, "TPC");
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
  h=new TH1I("hBC", "TOF Bunch Cross statistics;bc index;Entries", 21, -0.5, 20.5);
  fContainer->AddAt(h, kBunchCross);
  PostData(AliTRDpwg1Helper::kMonitor, fContainer);
}

//____________________________________________________________________
Bool_t AliTRDinfoGen::Load(const Char_t *file, const Char_t *dir, const Char_t *name)
{
// Load data from performance file

  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(dir){
    if(!gFile->cd(dir)){
      AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
      return kFALSE;
    }
  }
  TObjArray *o(NULL);
  const Char_t *tn=(name ? name : GetName());
  if(!(o = (TObjArray*)gDirectory->Get(tn))){
    AliWarning(Form("Missing histogram container %s.", tn));
    return kFALSE;
  }
  fContainer = (TObjArray*)o->Clone(GetName());
  gFile->Close();
  return kTRUE;
}

//____________________________________________________________________
void AliTRDinfoGen::UserExec(Option_t *){
  //
  // Run the Analysis
  //

  fTracksBarrel->Delete();
  fTracksSA->Delete();
  fTracksKink->Delete();
  fV0List->Delete();
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
      AliInfo("OCDB :: initialized.");
    } else {
      AliInfo("OCDB :: initializing locally ...");
      // prepare OCDB access
      ocdb->SetDefaultStorage(fOCDB.Data());
      ocdb->SetRun(fESDev->GetRunNumber());
      // create geo manager
      if(!(obj = ocdb->Get(AliCDBPath("GRP", "Geometry", "Data")))){
        AliError("GEOMETRY failed initialization.");
      } else {
        AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
        AliGeomManager::GetNalignable("TRD");
        AliGeomManager::ApplyAlignObjsFromCDB("TRD");
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
      if(AliLog::GetDebugLevel("PWG1/TRD", "AliTRDinfoGen")>5) reco->Dump();
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

  // event selection : trigger cut
  if(UseLocalEvSelection() && fEvTrigger){ 
    Bool_t kTRIGGERED(kFALSE);
    std::auto_ptr<TObjArray> trig(fEvTrigger->Tokenize(" "));
    for(Int_t itrig=trig->GetEntriesFast(); itrig--;){
      const Char_t *trigClass(((TObjString*)(*trig)[itrig])->GetName());
      if(fESDev->IsTriggerClassFired(trigClass)) {
        AliDebug(2, Form("Ev[%4d] Trigger[%s]", fESDev->GetEventNumberInFile(), trigClass));
        kTRIGGERED = kTRUE;
        break; 
      }
    }
    if(!kTRIGGERED){ 
      AliDebug(2, Form("Reject Ev[%4d] Trigger", fESDev->GetEventNumberInFile()));
      return;
    }
    // select only physical events
    if(fESDev->GetEventType() != 7){ 
      AliDebug(2, Form("Reject Ev[%4d] EvType[%d]", fESDev->GetEventNumberInFile(), fESDev->GetEventType()));
      return;
    }
  }

  // if the required trigger is a collision trigger then apply event vertex cut
  if(UseLocalEvSelection() && IsCollision()){
    const AliESDVertex *vertex = fESDev->GetPrimaryVertex();
    if(TMath::Abs(vertex->GetZv())<1.e-10 || 
       TMath::Abs(vertex->GetZv())>fgkEvVertexZ || 
       vertex->GetNContributors()<fgkEvVertexN) {
      AliDebug(2, Form("Reject Ev[%4d] Vertex Zv[%f] Nv[%d]", fESDev->GetEventNumberInFile(), TMath::Abs(vertex->GetZv()), vertex->GetNContributors()));
      return;
    }
  }

  if(fEventCut && !fEventCut->IsSelected(fESDev, IsCollision())) return;

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
  Int_t centralityBin = -1;
  AliDebug(2, Form("  Beam Type: %s", fESDev->GetESDRun()->GetBeamType()));
  if(TString(fESDev->GetESDRun()->GetBeamType()).Contains("Pb-Pb")){
    centralityBin = 4;
    const AliMultiplicity *mult = fESDev->GetMultiplicity();
    Double_t zdcNeutronEnergy = fESDev->GetZDCN1Energy()+fESDev->GetZDCN2Energy();
    Double_t itsNTracklets = mult->GetNumberOfTracklets();
    Double_t centralitySlopes[6] = {0.0, 4.0, 8.0, 20.0, 50.0, 1000000.};
    AliDebug(1, Form("zdcNeutronEnergy: %f, itsNTracklets: %f\n", zdcNeutronEnergy, itsNTracklets));
    for(Int_t iCent=1; iCent<=5; ++iCent) {
      if(zdcNeutronEnergy>centralitySlopes[iCent-1]*itsNTracklets && zdcNeutronEnergy<centralitySlopes[iCent]*itsNTracklets)
        centralityBin=iCent - 1;
    }
    AliDebug(2, Form("  Centrality Class: %d", centralityBin));
  }
  fEventInfo->SetCentrality(centralityBin);
  AliDebug(2, Form("  Bunch Crosses: %d", fESDev->GetBunchCrossNumber()));

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
  Int_t nTRDout(0), nTRDin(0), nTPC(0)
       ,nclsTrklt
       ,nBarrel(0), nSA(0), nKink(0)
       ,nBarrelFriend(0), nSAFriend(0)
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
    if(!(v0 = fESDev->GetV0(iv0))) continue;
    // register v0
    if(fV0Cut) new(fV0Info) AliTRDv0Info(*fV0Cut);
    else new(fV0Info) AliTRDv0Info();
    fV0Info->SetMagField(bField);
    fV0Info->SetV0tracks(fESDev->GetTrack(v0->GetPindex()), fESDev->GetTrack(v0->GetNindex()));
    fV0Info->SetV0Info(v0);
    fV0List->Add(new AliTRDv0Info(*fV0Info));//  kFOUND=kFALSE;
  }


  // LOOP 1 - over ESD tracks
  AliTRDv0Info *v0info=NULL;
  for(Int_t itrk = 0; itrk < nTracksESD; itrk++){
    new(fTrackInfo) AliTRDtrackInfo();
    esdTrack = fESDev->GetTrack(itrk);
    AliDebug(3, Form("\n%3d ITS[%d] TPC[%d] TRD[%d] TOF-BC[%d]\n", itrk, esdTrack->GetNcls(0), esdTrack->GetNcls(1), esdTrack->GetNcls(2), esdTrack->GetTOFBunchCrossing()));

    if(esdTrack->GetStatus()&AliESDtrack::kTPCout) nTPC++;
    if(esdTrack->GetStatus()&AliESDtrack::kTRDout) nTRDout++;
    if(esdTrack->GetStatus()&AliESDtrack::kTRDin) nTRDin++;

    // look at external track param
    const AliExternalTrackParam *op = esdTrack->GetOuterParam();
    Double_t xyz[3];
    if(op){
      op->GetXYZ(xyz);
      op->Global2LocalPosition(xyz, op->GetAlpha());
      AliDebug(3, Form("op @ X[%7.3f]\n", xyz[0]));
    }

    // read MC info
    Int_t fPdg = -1;
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
      AliMCParticle *mcParticle = NULL; 
      if(!(mcParticle = (AliMCParticle*) fMCev->GetTrack(alab))){
        AliError(Form("MC particle label[%d] missing for Ev[%d] Trk[%d].", label, (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), itrk));
        continue;
      }
      fPdg = mcParticle->Particle()->GetPdgCode();
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      Int_t iref = 0; AliTrackReference *ref = NULL; 
      while(iref<nRefs){
        ref = mcParticle->GetTrackReference(iref);
        if(ref->LocalX() > fgkTPC) break;
        iref++;
      }

      fTrackInfo->SetMC();
      fTrackInfo->SetPDG(fPdg);
      fTrackInfo->SetPrimary(mcParticle->Particle()->IsPrimary());
      fTrackInfo->SetLabel(label);
      Int_t jref = iref;//, kref = 0;
      while(jref<nRefs){
        ref = mcParticle->GetTrackReference(jref);
        if(ref->LocalX() > fgkTRD) break;
        AliDebug(4, Form("  trackRef[%2d (%2d)] @ %7.3f OK", jref-iref, jref, ref->LocalX()));
        fTrackInfo->AddTrackRef(ref);
        jref++;
      }
      AliDebug(3, Form("NtrackRefs[%d(%d)]", fTrackInfo->GetNTrackRefs(), nRefs));
    }

    // copy some relevant info to TRD track info
    fTrackInfo->SetStatus(esdTrack->GetStatus());
    fTrackInfo->SetTrackId(esdTrack->GetID());
    Double_t p[AliPID::kSPECIES]; esdTrack->GetTRDpid(p);
    fTrackInfo->SetESDpid(p);
    fTrackInfo->SetESDpidQuality(esdTrack->GetTRDntrackletsPID());
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
    fTrackInfo->SetTOFbc(static_cast<UShort_t>(esdTrack->GetTOFBunchCrossing()));
    nclsTrklt = 0;
  
    // set V0pid info
    for(Int_t iv(0); iv<fV0List->GetEntriesFast(); iv++){
      if(!(v0info = (AliTRDv0Info*)fV0List->At(iv))) continue;
      if(!v0info->GetV0Daughter(1) && !v0info->GetV0Daughter(-1)) continue;
      if(!v0info->HasTrack(fTrackInfo)) continue;
      memset(v0pid, 0, AliPID::kSPECIES*sizeof(Int_t));
      fTrackInfo->SetV0();
      for(Int_t is=AliPID::kSPECIES; is--;){v0pid[is] = v0info->GetPID(is, fTrackInfo);}
      fTrackInfo->SetV0pid(v0pid);
      fTrackInfo->SetV0();
      //const AliTRDtrackInfo::AliESDinfo *ei = fTrackInfo->GetESDinfo();
      break;
    }

    // read REC info
    esdFriendTrack = (fESDfriend->GetNumberOfTracks() > itrk) ? fESDfriend->GetTrack(itrk): NULL;

    if(esdFriendTrack){
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
        h = (TH1I*)fContainer->At(kBunchCross); h->Fill(esdTrack->GetTOFBunchCrossing());
        break;
      }
      AliDebug(3, Form("Ntracklets[%d]\n", fTrackInfo->GetNTracklets()));
    } else AliDebug(3, "No ESD friends");
    if(op) fTrackInfo->SetOuterParam(op);

    if(DebugLevel() >= 1){
      AliTRDtrackInfo info(*fTrackInfo);
      (*DebugStream()) << "trackInfo"
      << "TrackInfo.=" << &info
      << "\n";
      info.Delete("");
    }

    ULong_t status(esdTrack->GetStatus());
    if((status&AliESDtrack::kTPCout)){
      if(!esdTrack->GetKinkIndex(0)){ // Barrel  Track Selection
        Bool_t selected(kTRUE);
        if(UseLocalTrkSelection()){
          if(esdTrack->Pt() < fgkPt){ 
            AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] Pt[%5.2f]", itrk, fESDev->GetEventNumberInFile(), esdTrack->Pt()));
            selected = kFALSE;
          }
          if(selected && TMath::Abs(esdTrack->Eta()) > fgkEta){
            AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] Eta[%5.2f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(esdTrack->Eta())));
            selected = kFALSE;
          }
          if(selected && esdTrack->GetTPCNcls() < fgkNclTPC){ 
            AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] NclTPC[%d]", itrk, fESDev->GetEventNumberInFile(), esdTrack->GetTPCNcls()));
            selected = kFALSE;
          }
          Float_t par[2], cov[3];
          esdTrack->GetImpactParameters(par, cov);
          if(IsCollision()){ // cuts on DCA
            if(selected && TMath::Abs(par[0]) > fgkTrkDCAxy){ 
              AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] DCAxy[%f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(par[0])));
              selected = kFALSE;
            }
            if(selected && TMath::Abs(par[1]) > fgkTrkDCAz){ 
              AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] DCAz[%f]", itrk, fESDev->GetEventNumberInFile(), TMath::Abs(par[1])));
              selected = kFALSE;
            }
          } else if(selected && fMCev && !fMCev->IsPhysicalPrimary(alab)){;
            AliDebug(3, Form("Reject Trk[%3d] Ev[%4d] Primary", itrk, fESDev->GetEventNumberInFile()));
            selected = kFALSE;
          }
        }
        if(fTrackCut && !fTrackCut->IsSelected(esdTrack)) selected = kFALSE;
        if(selected){ 
          fTracksBarrel->Add(new AliTRDtrackInfo(*fTrackInfo));
          nBarrel++;
          if(fTrackInfo->GetTrack()) 
            nBarrelFriend++;
        }
      } else {
        fTracksKink->Add(new AliTRDtrackInfo(*fTrackInfo));
        nKink++;
      }
    } else if((status&AliESDtrack::kTRDout) && !(status&AliESDtrack::kTRDin)){ 
      fTracksSA->Add(new AliTRDtrackInfo(*fTrackInfo));
      nSA++;
      if(fTrackInfo->GetTrack()) 
        nSAFriend++;
    }
    fTrackInfo->Delete("");
  }

  // LOOP 2 - over MC tracks which are passing TRD where the track is not reconstructed
  if(HasMCdata()){
    AliDebug(10, "Output of the MC track map:");
    for(Int_t itk = 0; itk < nTracksMC;  itk++) AliDebug(10, Form("trackMap[%d] = %s", itk, trackMap[itk] == kTRUE ? "TRUE" : "kFALSE"));
  
    for(Int_t itk = 0; itk < nTracksMC; itk++){
      if(trackMap[itk]) continue;
      AliMCParticle *mcParticle =  (AliMCParticle*) fMCev->GetTrack(TMath::Abs(itk));
      Int_t fPdg = mcParticle->Particle()->GetPdgCode();
      Int_t nRefs = mcParticle->GetNumberOfTrackReferences();
      Int_t iref = 0; AliTrackReference *ref = NULL; 
      Int_t nRefsTRD = 0;
      new(fTrackInfo) AliTRDtrackInfo();
      fTrackInfo->SetMC();
      fTrackInfo->SetPDG(fPdg);
      while(iref<nRefs){ // count TRD TR
        Bool_t kIN(kFALSE);
        ref = mcParticle->GetTrackReference(iref);
        if(ref->LocalX() > fgkTPC && ref->LocalX() < fgkTRD){
          fTrackInfo->AddTrackRef(ref);
          nRefsTRD++;kIN=kTRUE;
        }
        AliDebug(4, Form("  trackRef[%2d] @ x[%7.3f] %s", iref, ref->LocalX(), kIN?"IN":"OUT"));
        iref++;
      }
      if(!nRefsTRD){
        // In this stage we at least require 1 hit inside TRD. What will be done with this tracks is a task for the 
        // analysis job
        fTrackInfo->Delete("");
        continue;
      }
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
    "        TPCout[%d] TRDin[%d] TRDout[%d]\n"
    "        Barrel[%3d+%3d=%3d] SA[%2d+%2d=%2d] Kink[%2d+%2d=%2d]"
    ,(Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nTracksESD, nTracksMC, fV0List->GetEntries()
    , nTPC, nTRDin, nTRDout
    ,nBarrel, nBarrelMC, fTracksBarrel->GetEntries()
    ,nSA, nSAMC, fTracksSA->GetEntries()
    ,nKink, nKinkMC, fTracksKink->GetEntries()
  ));
  // save track statistics
  h = (TH1I*)fContainer->At(kStatTrk);
  h->Fill(Float_t(kTracksESD), nTracksESD);
  h->Fill(Float_t(kTracksMC), nTracksMC);
  h->Fill(Float_t(kV0), fV0List->GetEntries());
  h->Fill(Float_t(kTPC), nTPC);
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

  PostData(AliTRDpwg1Helper::kTracksBarrel, fTracksBarrel);
  PostData(AliTRDpwg1Helper::kTracksSA, fTracksSA);
  PostData(AliTRDpwg1Helper::kTracksKink, fTracksKink);
  PostData(AliTRDpwg1Helper::kEventInfo, fEventInfo);
  PostData(AliTRDpwg1Helper::kV0List, fV0List);
}

//____________________________________________________________________
void AliTRDinfoGen::SetLocalV0Selection(const AliTRDv0Info *v0)
{
// Set V0 cuts from outside

  if(!fV0Cut) fV0Cut = new AliTRDv0Info(*v0);
  else new(fV0Cut) AliTRDv0Info(*v0);
}

//____________________________________________________________________
void AliTRDinfoGen::SetTrigger(const Char_t *trigger)
{
  if(!fEvTrigger) fEvTrigger = new TString(trigger);
  else (*fEvTrigger) = trigger;
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



