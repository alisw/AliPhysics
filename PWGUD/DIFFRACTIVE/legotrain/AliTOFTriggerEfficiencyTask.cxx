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
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TObjString.h"
#include "TH2D.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliTOFTriggerMask.h"
#include "AliTriggerIR.h"
#include "TArrayI.h"

ClassImp(AliTOFTriggerEfficiencyTask)

//-----------------------------------------------------------------------------
AliTOFTriggerEfficiencyTask::AliTOFTriggerEfficiencyTask(const char* name) :
AliAnalysisTaskSE(name),
fListOfHistos(NULL),
fTree(NULL),
fTriggersVsRun(NULL),
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
fTriggerMask(),
fTOFhits(),
fTOFhitTimes()
{
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
void AliTOFTriggerEfficiencyTask::UserCreateOutputObjects(){
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();
  fTriggersVsRun = new TH2D("fTriggersVsRun","",10,0,10,30000,270000,300000);
  fListOfHistos->Add(fTriggersVsRun);
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fClassesFired",&fClassesFired);
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
  fTree->Branch("fL0inputs",&fL0inputs);
  fTree->Branch("fIR1",&fIR1);
  fTree->Branch("fIR2",&fIR2);
  fTree->Branch("fTriggerMask",&fTriggerMask,"fTriggerMask[72]/i");
  fTree->Branch("fTOFhits",&fTOFhits);
  fTree->Branch("fTOFhitTimes",&fTOFhitTimes);

  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliTOFTriggerEfficiencyTask::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  fRunNumber  = fInputEvent->GetRunNumber();
  Bool_t isMC = (fMCEvent!=NULL);

  if (!isMC) {
    Bool_t accept = 0;
    if (fClassesFired.String().Contains("CINT7-B"))    { accept = 1; fTriggersVsRun->Fill(0.5,fRunNumber); }
    if (fClassesFired.String().Contains("CINT7ZAC-B")) { accept = 1; fTriggersVsRun->Fill(1.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP13-B"))   { accept = 1; fTriggersVsRun->Fill(2.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP25-B"))   { accept = 1; fTriggersVsRun->Fill(3.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP26-B"))   { accept = 1; fTriggersVsRun->Fill(4.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP29-B"))   { accept = 1; fTriggersVsRun->Fill(5.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP29-U"))   { accept = 1; fTriggersVsRun->Fill(6.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP30-B"))   { accept = 1; fTriggersVsRun->Fill(7.5,fRunNumber); }
    if (fClassesFired.String().Contains("CCUP31-B"))   { accept = 1; fTriggersVsRun->Fill(8.5,fRunNumber); }
    if (fClassesFired.String().Contains("CTRUE-B"))    { accept = 1; fTriggersVsRun->Fill(9.5,fRunNumber); }
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

  // avoid high multiplicity events
  if (fNofTracklets>100) { PostData(1,fListOfHistos); return; }
  
  for (UInt_t k=0;k<72;k++){
    fTriggerMask[k] = fInputEvent->GetTOFHeader()->GetTriggerMask()->GetTriggerMask(k);
  }

  const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(fInputEvent);
  if (!esd) return;
  
  TClonesArray* tofClusters = esd->GetESDTOFClusters();
  Int_t nTOFhits = 0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
    AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
    if (cl->GetNMatchableTracks()!=1) continue;
    nTOFhits+=cl->GetNTOFhits();
  }

  fTOFhits.Reset();
  fTOFhitTimes.Reset();
  fTOFhits.Set(nTOFhits);
  fTOFhitTimes.Set(nTOFhits);
  Int_t hitCounts=0;
  for (Int_t icl=0;icl<tofClusters->GetEntriesFast();icl++){
    AliESDTOFCluster* cl = (AliESDTOFCluster*) tofClusters->At(icl);
    if (cl->GetNMatchableTracks()!=1) continue;
    for (Int_t ihit=0;ihit<cl->GetNTOFhits();ihit++){
      AliESDTOFHit* hit = (AliESDTOFHit*) cl->GetTOFHit(ihit);
      Float_t t = hit->GetTime();
      Int_t channel = hit->GetTOFchannel();
      fTOFhits.AddAt(channel,hitCounts);
      fTOFhitTimes.AddAt(t,hitCounts);
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
  
  fTree->Fill();
  PostData(1,fListOfHistos);
  PostData(2,fTree);
}
