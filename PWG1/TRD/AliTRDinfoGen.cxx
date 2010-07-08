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
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTRDtrackV1.h"
#include "AliTrackReference.h"
#include "AliTRDgeometry.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "TTreeStream.h"

#include <cstdio>
#include <climits>
#include <cstring>
#include <iostream>

#include "AliTRDinfoGen.h"
#include "info/AliTRDtrackInfo.h"
#include "info/AliTRDeventInfo.h"
#include "info/AliTRDv0Info.h"
#include "info/AliTRDeventCuts.h"
#include "macros/AliTRDperformanceTrain.h"

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

//____________________________________________________________________
AliTRDinfoGen::AliTRDinfoGen()
  :AliAnalysisTaskSE()
  ,fEvTrigger(NULL)
  ,fESDev(NULL)
  ,fMCev(NULL)
  ,fEventCut(NULL)
  ,fTrackCut(NULL)
  ,fV0Cut(NULL)
  ,fTrackInfo(NULL)
  ,fEventInfo(NULL)
  ,fV0Info(NULL)
  ,fTracksBarrel(NULL)
  ,fTracksSA(NULL)
  ,fTracksKink(NULL)
  ,fV0List(NULL)
  ,fDebugStream(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("infoGen", "MC-REC TRD-track list generator");
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
  ,fTrackInfo(NULL)
  ,fEventInfo(NULL)
  ,fV0Info(NULL)
  ,fTracksBarrel(NULL)
  ,fTracksSA(NULL)
  ,fTracksKink(NULL)
  ,fV0List(NULL)
  ,fDebugStream(NULL)
{
  //
  // Default constructor
  //
  SetTitle("MC-REC TRD-track list generator");
  DefineOutput(kTracksBarrel, TObjArray::Class());
  DefineOutput(kTracksSA, TObjArray::Class());
  DefineOutput(kTracksKink, TObjArray::Class());
  DefineOutput(kEventInfo, AliTRDeventInfo::Class());
  DefineOutput(kV0List, TObjArray::Class());
}

//____________________________________________________________________
AliTRDinfoGen::~AliTRDinfoGen()
{
// Destructor
  
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
}

//____________________________________________________________________
void AliTRDinfoGen::UserExec(Option_t *){
  //
  // Run the Analysis
  //

  fESDev = dynamic_cast<AliESDEvent*>(InputEvent());
  fMCev = MCEvent();

  if(!fESDev){
    AliError("Failed retrieving ESD event");
    return;
  }

  // event selection : trigger cut
  if(UseLocalEvSelection() && fEvTrigger){ 
    Bool_t kTRIGGERED(kFALSE);
    const TObjArray *trig = fEvTrigger->Tokenize(" ");
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

  fTracksBarrel->Delete();
  fTracksSA->Delete();
  fTracksKink->Delete();
  fV0List->Delete();
  fEventInfo->Delete("");
  new(fEventInfo)AliTRDeventInfo(fESDev->GetHeader(), const_cast<AliESDRun *>(fESDev->GetESDRun()));
  
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
       ,nBarrelMC(0), nSAMC(0), nKinkMC(0);
  AliESDtrack *esdTrack = NULL;
  AliESDfriendTrack *esdFriendTrack = NULL;
  TObject *calObject = NULL;
  AliTRDtrackV1 *track = NULL;
  AliTRDseedV1 *tracklet = NULL;
  AliTRDcluster *cl = NULL;
  // LOOP 0 - over ESD tracks
  for(Int_t itrk = 0; itrk < nTracksESD; itrk++){
    new(fTrackInfo) AliTRDtrackInfo();
    esdTrack = fESDev->GetTrack(itrk);
    AliDebug(3, Form("\n%3d ITS[%d] TPC[%d] TRD[%d]\n", itrk, esdTrack->GetNcls(0), esdTrack->GetNcls(1), esdTrack->GetNcls(2)));

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
    nclsTrklt = 0;
  

    // read REC info
    esdFriendTrack = fESDfriend->GetTrack(itrk);
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
            AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] Pt[%5.2f]", fESDev->GetEventNumberInFile(), itrk, esdTrack->Pt()));
            selected = kFALSE;
          }
          if(selected && TMath::Abs(esdTrack->Eta()) > fgkEta){
            AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] Eta[%5.2f]", fESDev->GetEventNumberInFile(), itrk, TMath::Abs(esdTrack->Eta())));
            selected = kFALSE;
          }
          if(selected && esdTrack->GetTPCNcls() < fgkNclTPC){ 
            AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] NclTPC[%d]", fESDev->GetEventNumberInFile(), itrk, esdTrack->GetTPCNcls()));
            selected = kFALSE;
          }
          Float_t par[2], cov[3];
          esdTrack->GetImpactParameters(par, cov);
          if(IsCollision()){ // cuts on DCA
            if(selected && TMath::Abs(par[0]) > fgkTrkDCAxy){ 
              AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] DCAxy[%f]", fESDev->GetEventNumberInFile(), itrk, TMath::Abs(par[0])));
              selected = kFALSE;
            }
            if(selected && TMath::Abs(par[1]) > fgkTrkDCAz){ 
              AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] DCAz[%f]", fESDev->GetEventNumberInFile(), itrk, TMath::Abs(par[1])));
              selected = kFALSE;
            }
          } else if(selected && fMCev && !fMCev->IsPhysicalPrimary(alab)){;
            AliDebug(2, Form("Reject Ev[%4d] Trk[%3d] Primary", fESDev->GetEventNumberInFile(), itrk));
            selected = kFALSE;
          }
        }
        if(fTrackCut && !fTrackCut->IsSelected(esdTrack)) selected = kFALSE;
        if(selected){ 
          fTracksBarrel->Add(new AliTRDtrackInfo(*fTrackInfo));
          nBarrel++;
        }
      } else {
        fTracksKink->Add(new AliTRDtrackInfo(*fTrackInfo));
        nKink++;
      }
    } else if((status&AliESDtrack::kTRDout) && !(status&AliESDtrack::kTRDin)){ 
      fTracksSA->Add(new AliTRDtrackInfo(*fTrackInfo));
      nSA++;
    }
    fTrackInfo->Delete("");
  }

  // LOOP 1 - over ESD v0s
  Bool_t kFOUND(kFALSE);
  Float_t bField(fESDev->GetMagneticField());
  AliESDv0 *v0(NULL);
  for(Int_t iv0(0); iv0<fESDev->GetNumberOfV0s(); iv0++){
    if(!(v0 = fESDev->GetV0(iv0))) continue;

    // purge V0 irelevant for TRD
    Int_t jb(0); AliTRDtrackInfo *ti(NULL);
    for(Int_t ib(0); ib<nBarrel; ib++){ 
      ti =   (AliTRDtrackInfo*)fTracksBarrel->At(ib);
      Int_t id(ti->GetTrackId());
      if(id!=v0->GetPindex() && id!=v0->GetNindex()) continue; 
      kFOUND=kTRUE; ti->SetV0(); jb=ib+1;
      break;
    }
    if(!kFOUND) continue;

    // register v0
    if(fV0Cut) new(fV0Info) AliTRDv0Info(*fV0Cut);
    else new(fV0Info) AliTRDv0Info();
    fV0Info->SetMagField(bField);
    fV0Info->SetV0tracks(fESDev->GetTrack(v0->GetPindex()), fESDev->GetTrack(v0->GetNindex()));
    fV0Info->SetV0Info(v0);
    // mark the other leg too 
    for(Int_t ib(jb); ib<nBarrel; ib++){ 
      ti =   (AliTRDtrackInfo*)fTracksBarrel->At(ib);
      if(!fV0Info->HasTrack(ti)) continue;
      ti->SetV0();
      break;
    }
    //fV0Info->Print("a");
    fV0List->Add(new AliTRDv0Info(*fV0Info)); kFOUND=kFALSE;
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

  PostData(kTracksBarrel, fTracksBarrel);
  PostData(kTracksSA, fTracksSA);
  PostData(kTracksKink, fTracksKink);
  PostData(kEventInfo, fEventInfo);
  PostData(kV0List, fV0List);
}

//____________________________________________________________________
void AliTRDinfoGen::SetLocalV0Selection(AliTRDv0Info *v0)
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
  if(!fDebugStream){
    TDirectory *savedir = gDirectory;
    fDebugStream = new TTreeSRedirector("TRD.DebugInfoGen.root");
    savedir->cd();
  }
  return fDebugStream;
}


