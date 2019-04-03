/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

// Task to create a tree for TOF trigger efficiency studies
// evgeny.kryshen@cern.ch

#include "AliTOFTriggerEfficiencyTask.h"

#include "AliAnalysisTaskSE.h"
#include "TChain.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliVHeader.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TObjString.h"
#include "TH2D.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliUpcParticle.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliTOFTriggerMask.h"
#include "AliTriggerIR.h"
#include "AliESDtrackCuts.h"
#include "AliTOFGeometry.h"
#include "AliTrackerBase.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "TArrayI.h"

ClassImp(AliTOFTriggerEfficiencyTask)

//-----------------------------------------------------------------------------
AliTOFTriggerEfficiencyTask::AliTOFTriggerEfficiencyTask(const char* name) :
AliAnalysisTaskSE(name),
fTrackCutsBit0(NULL),
fTrackCutsBit5(NULL),
fMatOrig(), 
fMatCurr(),
fIsMC(0),
fListOfHistos(NULL),
fTree(NULL),
fTriggersVsRun(NULL),
fTracks(NULL),
fClassesFired(),
fRunNumber(0),
fPeriod(-1),
fOrbit(-1),
fBC(-1),
fL0inputs(0),
fVtxX(-1000),
fVtxY(-1000),
fVtxZ(-1000),
fVtxTPC(kFALSE),
fVtxContributors(0),
fNofTracklets(0),
fIR1(),
fIR2(),
fFOmap(),
fFiredChipMap(),
fTriggerMask(),
fTOFhits(),
fTOFhitTimes(),
fTrackIndices(),
fNofTOFtrgPads()

{
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliTOFTriggerEfficiencyTask::NotifyRun(){
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(fCurrentRunNumber);
  AliGRPManager* fGRPManager = new AliGRPManager();
  if(!fGRPManager->ReadGRPEntry()) AliFatal("Cannot get GRP entry"); 
  if(!fGRPManager->SetMagField())  AliFatal("Problem with magnetic field setup");
  if (!AliGeomManager::GetGeometry()) AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("ITS TRD TOF");
  for (int i=0;i<18;i++) {
    AliGeomManager::GetOrigGlobalMatrix( Form("TOF/sm%02d",i) ,fMatOrig[i]);
    fMatCurr[i] = *AliGeomManager::GetMatrix( Form("TOF/sm%02d",i) );
  }
}

//-----------------------------------------------------------------------------
void AliTOFTriggerEfficiencyTask::UserCreateOutputObjects(){
  fTrackCutsBit5 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fTrackCutsBit0 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fTriggersVsRun = new TH2D("fTriggersVsRun","",8,0,8,30000,270000,300000);
  fListOfHistos->Add(fTriggersVsRun);
  fTracks      = new TClonesArray("AliUpcParticle",10);
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fClassesFired",&fClassesFired);
  fTree->Branch("fTracks",&fTracks);
  fTree->Branch("fPeriod",&fPeriod);
  fTree->Branch("fOrbit",&fOrbit);
  fTree->Branch("fBC",&fBC);
  fTree->Branch("fRunNumber",&fRunNumber);
  fTree->Branch("fNofTracklets",&fNofTracklets);
  fTree->Branch("fVtxX",&fVtxX);
  fTree->Branch("fVtxY",&fVtxY);
  fTree->Branch("fVtxZ",&fVtxZ);
  fTree->Branch("fVtxTPC",&fVtxTPC);
  fTree->Branch("fVtxContributors",&fVtxContributors);
  fTree->Branch("fNofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("fL0inputs",&fL0inputs);
  fTree->Branch("fIR1",&fIR1);
  fTree->Branch("fIR2",&fIR2);
  fTree->Branch("fFOmap",&fFOmap);
  fTree->Branch("fFiredChipMap",&fFiredChipMap);
  fTree->Branch("fTriggerMask",&fTriggerMask,"fTriggerMask[72]/i");
  fTree->Branch("fTOFhits",&fTOFhits);
  fTree->Branch("fTOFhitTimes",&fTOFhitTimes);
  fTree->Branch("fTrackIndices",&fTrackIndices);
  fTree->Branch("fNofTOFtrgPads",&fNofTOFtrgPads);

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliTOFTriggerEfficiencyTask::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  if (!fIsMC) {
    Bool_t accept = 0;
    if (fClassesFired.String().Contains("CINT7-B"))    { accept = 1; fTriggersVsRun->Fill(0.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP25-B"))   { accept = 1; fTriggersVsRun->Fill(1.5,fRunNumber); }
    if (fClassesFired.String().Contains("CINT7ZAC-B")) { accept = 1; fTriggersVsRun->Fill(2.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP30-B"))   { accept = 1; fTriggersVsRun->Fill(3.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP31-B"))   { accept = 1; fTriggersVsRun->Fill(4.5,fRunNumber); }
    if (!accept) { PostData(1,fListOfHistos); return; }
  } else {
    fTriggersVsRun->Fill(0.5,fRunNumber); 
  }

  Bool_t int7selected = fInputHandler->IsEventSelected() & AliVEvent::kINT7;
  if (fClassesFired.String().Contains("CINT7-B") && !int7selected) { PostData(1,fListOfHistos); return; }
  
  fPeriod       = fInputEvent->GetPeriodNumber();
  fOrbit        = fInputEvent->GetOrbitNumber();
  fBC           = fInputEvent->GetBunchCrossNumber();
  fL0inputs     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();

  AliVMultiplicity* mult = fInputEvent->GetMultiplicity();
  fNofTracklets = mult->GetNumberOfTracklets();
  fFOmap        = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  // avoid high multiplicity events
  if (fNofTracklets>100) { PostData(1,fListOfHistos); return; }
  
  for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fInputEvent->GetNumberOfITSClusters(i);

  for (UInt_t k=0;k<72;k++){
    fTriggerMask[k] = fInputEvent->GetTOFHeader()->GetTriggerMask()->GetTriggerMask(k);
  }
  fNofTOFtrgPads = fInputEvent->GetTOFHeader()->GetNumberOfTOFtrgPads();

  const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(fInputEvent);

  TClonesArray* tofClusters = esd->GetESDTOFClusters();

  Int_t nTOFhits = 0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
    AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
    if (cl->GetNMatchableTracks()!=1) continue;
    nTOFhits+=cl->GetNTOFhits();
  }

  fTOFhits.Reset();
  fTOFhitTimes.Reset();
  fTrackIndices.Reset();
  fTOFhits.Set(nTOFhits);
  fTOFhitTimes.Set(nTOFhits);
  fTrackIndices.Set(nTOFhits);
  Int_t hitCounts=0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
    AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
    if (cl->GetNMatchableTracks()!=1) continue;
    Int_t trackIndex = cl->GetTrackIndex(0);
    for (Int_t ihit=0;ihit<cl->GetNTOFhits();ihit++){
      AliESDTOFHit* hit = (AliESDTOFHit*) cl->GetTOFHit(ihit);
      Float_t t = hit->GetTime();
      Int_t channel = hit->GetTOFchannel();
      fTOFhits.AddAt(channel,hitCounts);
      fTOFhitTimes.AddAt(t,hitCounts);
      fTrackIndices.AddAt(trackIndex,hitCounts);
      hitCounts++;
    }
  }

  const AliVVertex* vertex  = fInputEvent->GetPrimaryVertex();
  fVtxX   = vertex->GetX();
  fVtxY   = vertex->GetY();
  fVtxZ   = vertex->GetZ();
  fVtxTPC = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");
  fVtxContributors = vertex->GetNContributors();
  
  // vertex selection
  if (fVtxX<-0.5 || fVtxX>0.5 || fVtxY<-0.5 || fVtxY>1.0 || fVtxZ<-10 || fVtxZ>10 || fVtxContributors<1) 
  { PostData(1,fListOfHistos); return; }
  
  fTracks->Clear();
  for(Int_t iTrack=0; iTrack<fInputEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(iTrack);
    if (!track) continue;
    if(!fTrackCutsBit5->AcceptTrack(track))continue;
    Float_t pt     = track->Pt();
    if (pt<0.6) continue;
    ULong_t status = track->GetStatus();
    if (!(status & AliVTrack::kTOFin)) continue;
    Float_t eta    = track->Eta();
    Float_t phi0   = track->Phi();
    Short_t charge = track->Charge();

    AliESDtrack* esdTrack = (AliESDtrack*) track->Clone();
    AliExternalTrackParam* trc = (AliExternalTrackParam*) esdTrack->GetOuterParam();
    if (!trc) { delete esdTrack; continue; }
    // assuming trc is the copy of the track at the vertex or TPC outer param, w/o any preliminary rotations.
    // propagate to approximate X of TOF in arbitrary frame, rotating (kTRUE in PropagateTrackToBxByBz below)
    // the frame at each step in such a way that the cos(phi) is ~1.
    Float_t rTOF = AliTOFGeometry::RinTOF();
    if (!AliTrackerBase::PropagateTrackToBxByBz(trc,rTOF,esdTrack->GetMassForTracking(), 1., kTRUE, 0.8)) continue; // propagation failed

    // go to specific frame sector and reach target X in this frame:
    double phi = trc->PhiPos()*TMath::RadToDeg();
    if (phi<0) phi += 360;
    int sectOld = -1, sect = int(phi/20.);
    Bool_t failed = kFALSE;
    while (sectOld!=sect) {
      sectOld = sect;
      double alpha = (sect*20.+10)*TMath::DegToRad();
      if (!trc->Rotate(alpha) || !AliTrackerBase::PropagateTrackToBxByBz(trc,rTOF,esdTrack->GetMassForTracking(), 1., kFALSE)) { // don't rotate at every step anymore
        failed = kTRUE; break;
      }
      // make sure the propagation did not change the sector
      phi = trc->PhiPos()*TMath::RadToDeg();
      if (phi<0) phi += 360;
      sect = int(phi/20.);
    }
    if (failed) { delete esdTrack; continue;}
    
    Double_t cv[21];
    trc->GetCovarianceXYZPxPyPz(cv);
//    Double_t* cov = trc->GetCovariance();
    Double_t sigma_y2 = trc->GetSigmaY2();
    Double_t sigma_z2 = trc->GetSigmaZ2();
//    printf("%f\n",sigma_y2);
//    trc->
    Double_t pos[3]={0};
    Double_t locTmp[3]={0};
    Float_t posm[3]={0};
    Float_t posf[3]={0};
    Int_t detId[5] = {-1,-1,-1,-1,-1};
    Int_t inCounts = 0;
    // search for TOF strips in 1mm steps
    while (rTOF<AliTOFGeometry::Rmax()){
      rTOF += 0.1; // 1mm
      if(!AliTrackerBase::PropagateTrackParamOnlyToBxByBz(trc,rTOF,1,kFALSE)) break;
      phi = trc->PhiPos()*TMath::RadToDeg();
      if (phi<0) phi += 360;
      sect = int(phi/20.);
      if (sect!=sectOld) {
        sectOld = sect;
        double alpha = (sect*20.+10)*TMath::DegToRad();
        if (!trc->Rotate(alpha)) break;
      }
      trc->GetXYZ(pos);
      posm[0]=pos[0];
      posm[1]=pos[1];
      posm[2]=pos[2];
      fMatCurr[sect].MasterToLocal(pos,locTmp); // go to sector local frame, accounting for misalignment
      fMatOrig[sect].LocalToMaster(locTmp,pos); // go back to ideal lab frame
      posf[0]=pos[0];
      posf[1]=pos[1];
      posf[2]=pos[2];
      
      AliTOFGeometry::GetDetID(posf,detId);
      Float_t padDx = AliTOFGeometry::GetPadDx(posf);
      Float_t padDy = AliTOFGeometry::GetPadDy(posf);
      Float_t padDz = AliTOFGeometry::GetPadDz(posf);
      if (detId[4]==-1 || padDy<-0.1 || padDy>0.1) {
        inCounts = 0; 
        continue;
      }
      inCounts++;
      if (inCounts>1) continue;
      if (padDx<-1.25 || padDx>1.25) continue;
      if (padDz<-1.75 || padDz>1.75) continue;
      const Int_t index = AliTOFGeometry::GetIndex(detId); // [3]=padz,[4]=padx
      if (index<0) continue;
      float* fstatus = (float *) &status;
      AliUpcParticle* part = new ((*fTracks)[fTracks->GetEntriesFast()]) AliUpcParticle(pt,eta,phi0,charge,index,18);
      part->SetAt(posf[0],0);
      part->SetAt(posf[1],1);
      part->SetAt(posf[2],2);
      part->SetAt(cv[0],3);
      part->SetAt(cv[2],4);
      part->SetAt(cv[5],5);
      part->SetAt(padDx,6);
      part->SetAt(padDy,7);
      part->SetAt(padDz,8);
      part->SetAt(posm[0],9);
      part->SetAt(posm[1],10);
      part->SetAt(posm[2],11);
      part->SetAt(track->GetTPCCrossedRows(),12);
      part->SetAt(track->GetTOFsignal(),13);
      part->SetAt(iTrack,14);
      part->SetAt(*(fstatus),15);
      part->SetAt(sigma_y2,16);
      part->SetAt(sigma_z2,17);
    }
    delete esdTrack;
  }

  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------
