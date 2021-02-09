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

// Task to create a tree for SPD trigger efficiency studies
// evgeny.kryshen@cern.ch

#include "AliSPDTriggerEfficiencyTask.h"

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

ClassImp(AliSPDTriggerEfficiencyTask)

//-----------------------------------------------------------------------------
AliSPDTriggerEfficiencyTask::AliSPDTriggerEfficiencyTask(const char* name) :
AliAnalysisTaskSE(name),
fTrackCutsBit0(NULL),
fTrackCutsBit5(NULL),
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
fNofITSClusters(),
fIR1(),
fIR2(),
fFOmap(),
fFiredChipMap()
{
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliSPDTriggerEfficiencyTask::UserCreateOutputObjects(){
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

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliSPDTriggerEfficiencyTask::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  if (!fMCEvent) {
    Bool_t accept = 0;
    if (fClassesFired.String().Contains("CINT7-B"))     { accept = 1; fTriggersVsRun->Fill(0.5,fRunNumber); }
    if (fClassesFired.String().Contains("CINT7ZAC-B"))  { accept = 1; fTriggersVsRun->Fill(1.5,fRunNumber); }
    if (!accept) { PostData(1,fListOfHistos); return; }
  } else {
    fTriggersVsRun->Fill(0.5,fRunNumber); 
  }

  Bool_t int7selected = fInputHandler->IsEventSelected() & AliVEvent::kINT7;
  if (fClassesFired.String().Contains("CINT7") && !int7selected) { PostData(1,fListOfHistos); return; }
  
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

  const AliVVertex* vertex  = fInputEvent->GetPrimaryVertex();
  fVtxX   = vertex->GetX();
  fVtxY   = vertex->GetY();
  fVtxZ   = vertex->GetZ();
  fVtxTPC = TString(vertex->GetName()).CompareTo("PrimaryVertex") && TString(vertex->GetName()).CompareTo("SPDVertex");
  fVtxContributors = vertex->GetNContributors();
  
  // vertex selection
  if (fVtxX<-0.5 || fVtxX>0.5 || fVtxY<-0.5 || fVtxY>1.0 || fVtxContributors<2) 
  { PostData(1,fListOfHistos); return; }
  
  fTracks->Clear();
  for (Int_t ipart=0;ipart<fInputEvent->GetNumberOfTracks();ipart++){
    AliESDtrack* track = (AliESDtrack*) fInputEvent->GetTrack(ipart);
    if (!track) continue;
    Bool_t isBit0 = fTrackCutsBit0->AcceptTrack(track);
    Bool_t isBit5 = fTrackCutsBit5->AcceptTrack(track);
    if (!isBit0) continue;
    UChar_t itsMap = track->GetITSClusterMap();
    if (!TESTBIT(itsMap,0) && !TESTBIT(itsMap,1)) continue;
    UInt_t filterMap = 0;
    filterMap |= isBit0 << 0;
    filterMap |= isBit5 << 5;
    Double_t p[3];
    track->GetPxPyPz(p);
    AliUpcParticle* part = new ((*fTracks)[fTracks->GetEntriesFast()]) AliUpcParticle(10,3);
    part->SetI(0,track->Charge());
    part->SetI(1,track->GetLabel());
    part->SetI(2,track->GetITSClusterMap());
    part->SetI(3,track->GetStatus());
    part->SetI(4,filterMap);
    part->SetI(5,track->GetTPCCrossedRows());
    part->SetI(6,track->GetITSModuleIndex(0));
    part->SetI(7,track->GetITSModuleIndex(1));
    part->SetI(8,track->GetITSModuleIndex(6));
    part->SetI(9,track->GetITSModuleIndex(7));
    part->SetF(0,p[0]);
    part->SetF(1,p[1]);
    part->SetF(2,p[2]);
  }

  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------
