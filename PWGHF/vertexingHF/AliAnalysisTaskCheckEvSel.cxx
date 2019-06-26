/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

//*************************************************************************
// Class AliAnalysisTaskCheckEvSel
// Task to check event selection with D2H code
// Author: Francesco Prino
//*************************************************************************


#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliVertexingHFUtils.h"
#include "AliNormalizationCounter.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAnalysisTaskCheckEvSel.h"

ClassImp(AliAnalysisTaskCheckEvSel)


//________________________________________________________________________
AliAnalysisTaskCheckEvSel::AliAnalysisTaskCheckEvSel():
  AliAnalysisTaskSE(),
  fOutput(0x0),
  fReadMC(kFALSE),
  fSystem(0),
  fCutOnzVertexSPD(0),
  fAODProtection(1),
  fHistNEvents(0x0),
  fHistNEventsVsCent(0x0),
  fHistNEventsVsCL1(0x0),
  fHistWhyRej(0x0),
  fHistNEventsVsWhyRej(0x0),
  fHistNEventsVsTime(0x0),
  fHistNTrackletsBeforePileup(0x0),
  fHistNTrackletsAfterPileup(0x0),
  fHistNCL1BeforePileup(0x0),
  fHistNCL1AfterPileup(0x0),
  fHistCentrality(0x0),
  fHistCL0vsV0MCentrality(0x0),
  fHistNTracksTPCoutVsV0Cent(0x0),
  fHistNTracksFB4VsV0Cent(0x0),
  fHistNTracksFB4EtaPosVsV0Cent(0x0),
  fHistNTracksFB4EtaNegVsV0Cent(0x0),
  fHistNTracksBC0VsV0Cent(0x0),
  fHistNTrackletsVsV0Cent(0x0),
  fHistNTrackletsGoldenVsV0Cent(0x0),
  fHistNTrackletsGoldenVsV0CentVsZvert(0x0),  
  fHistPhiTrackelts(0x0),
  fHistNCL0VsV0Cent(0x0),
  fHistNCL1VsV0Cent(0x0),
  fHistT0AmplVsV0Ampl(0x0),
  fHistT0AmplVsV0Cent(0x0),
  fHistT0AmplVsNCL0(0x0),
  fHistT0AmplVsCL0Cent(0x0),
  fHistNTracksTPCoutVsNTracklets(0x0),
  fHistNTracksFB4VsNTracklets(0x0),
  fHistNTracksBC0VsNTracksFB4(0x0),
  fHistZVertexSPDBeforeCuts(0x0),
  fHistZVertexSPDBeforeSPDCut(0x0),
  fHistZVertexSPDAfterCuts(0x0),
  fHistZVertexSPDBadTrackVert(0x0),
  fEventProp(0x0),
  fNtupleEvProp(0x0),
  fEnableEvPropNtuple(kFALSE),
  fNtupleZvtxDistVsWhyRej(0x0),
  fEnableVertexNtuple(kFALSE),
  fUseAliEventCuts(kFALSE),
  fCounter(0),
  fAnalysisCuts(0x0)
{
  // default constructor
  fAnalysisCuts=new AliRDHFCutsD0toKpi("EvSelCuts");
  fAnalysisCuts->SetTriggerMask(AliVEvent::kMB | AliVEvent::kINT7);
  fAnalysisCuts->SetTriggerClass("");
  fAnalysisCuts->SetOptPileup(AliRDHFCuts::kNoPileupSelection);
}
//________________________________________________________________________
AliAnalysisTaskCheckEvSel::AliAnalysisTaskCheckEvSel(Bool_t readMC, Int_t system, AliRDHFCutsD0toKpi* cuts):
  AliAnalysisTaskSE("EvSelTask"),
  fOutput(0x0),
  fReadMC(readMC),
  fSystem(system),
  fCutOnzVertexSPD(0),
  fAODProtection(1),
  fHistNEvents(0x0),
  fHistNEventsVsCent(0x0),
  fHistNEventsVsCL1(0x0),
  fHistWhyRej(0x0),
  fHistNEventsVsWhyRej(0x0),
  fHistNEventsVsTime(0x0),
  fHistNTrackletsBeforePileup(0x0),
  fHistNTrackletsAfterPileup(0x0),
  fHistNCL1BeforePileup(0x0),
  fHistNCL1AfterPileup(0x0),
  fHistCentrality(0x0),
  fHistCL0vsV0MCentrality(0x0),
  fHistNTracksTPCoutVsV0Cent(0x0),
  fHistNTracksFB4VsV0Cent(0x0),
  fHistNTracksFB4EtaPosVsV0Cent(0x0),
  fHistNTracksFB4EtaNegVsV0Cent(0x0),
  fHistNTracksBC0VsV0Cent(0x0),
  fHistNTrackletsVsV0Cent(0x0),
  fHistNTrackletsGoldenVsV0Cent(0x0),
  fHistNTrackletsGoldenVsV0CentVsZvert(0x0),
  fHistPhiTrackelts(0x0),
  fHistNCL0VsV0Cent(0x0),
  fHistNCL1VsV0Cent(0x0),
  fHistT0AmplVsV0Ampl(0x0),
  fHistT0AmplVsV0Cent(0x0),
  fHistT0AmplVsNCL0(0x0),
  fHistT0AmplVsCL0Cent(0x0),
  fHistNTracksTPCoutVsNTracklets(0x0),
  fHistNTracksFB4VsNTracklets(0x0),
  fHistNTracksBC0VsNTracksFB4(0x0),
  fHistZVertexSPDBeforeCuts(0x0),
  fHistZVertexSPDBeforeSPDCut(0x0),
  fHistZVertexSPDAfterCuts(0x0),
  fHistZVertexSPDBadTrackVert(0x0),
  fEventProp(0x0),
  fNtupleEvProp(0x0),
  fEnableEvPropNtuple(kFALSE),
  fNtupleZvtxDistVsWhyRej(0x0),
  fEnableVertexNtuple(kFALSE),
  fUseAliEventCuts(kFALSE),
  fCounter(0),
  fAnalysisCuts(cuts)
{
  // default constructor

  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,AliNormalizationCounter::Class());
  DefineOutput(3,TNtuple::Class());
  DefineOutput(4,TNtuple::Class());
}

//________________________________________________________________________
AliAnalysisTaskCheckEvSel::~AliAnalysisTaskCheckEvSel()
{
  //
  // Destructor
  //
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNEventsVsCent;
    delete fHistNEventsVsCL1;
    delete fHistWhyRej;
    delete fHistNEventsVsWhyRej;
    delete fHistNEventsVsTime;
    delete fHistNTrackletsBeforePileup;
    delete fHistNTrackletsAfterPileup;
    delete fHistNCL1BeforePileup;
    delete fHistNCL1AfterPileup;
    delete fHistCentrality;
    delete fHistCL0vsV0MCentrality;
    delete fHistNTracksTPCoutVsV0Cent;
    delete fHistNTracksFB4VsV0Cent;
    delete fHistNTracksFB4EtaPosVsV0Cent;
    delete fHistNTracksFB4EtaNegVsV0Cent;
    delete fHistNTracksBC0VsV0Cent;
    delete fHistNTrackletsVsV0Cent;
    delete fHistNTrackletsGoldenVsV0Cent;
    delete fHistNTrackletsGoldenVsV0CentVsZvert;
    delete fHistPhiTrackelts;
    delete fHistNCL0VsV0Cent;
    delete fHistNCL1VsV0Cent;
    delete fHistT0AmplVsV0Ampl;
    delete fHistT0AmplVsV0Cent;
    delete fHistT0AmplVsNCL0;
    delete fHistT0AmplVsCL0Cent;
    delete fHistNTracksTPCoutVsNTracklets;
    delete fHistNTracksFB4VsNTracklets;
    delete fHistNTracksBC0VsNTracksFB4;
    delete fHistZVertexSPDBeforeCuts;
    delete fHistZVertexSPDBeforeSPDCut;
    delete fHistZVertexSPDAfterCuts;  
    delete fHistZVertexSPDBadTrackVert;
    delete fEventProp;
  }
  delete fOutput;
  delete fCounter;
  if(fNtupleZvtxDistVsWhyRej) delete fNtupleZvtxDistVsWhyRej;
  if(fNtupleEvProp) delete fNtupleEvProp;

}

//________________________________________________________________________
void AliAnalysisTaskCheckEvSel::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskCheckEvSel::UserCreateOutputObjects() \n");

  Double_t maxMult=200.;
  if(fSystem==1) maxMult=5000.;
  else if(fSystem==2) maxMult=500.;

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");
  
  fHistNEvents = new TH1F("hNEvents", "number of events ",21,-0.5,20.5);
  ConfigureEvSelAxis(fHistNEvents->GetXaxis());
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  fHistNEventsVsCent = new TH2F("hNEventsVsCent", " ; ; Centrality ",21,-0.5,20.5,101,0.,101.);
  ConfigureEvSelAxis(fHistNEventsVsCent->GetXaxis());
  fOutput->Add(fHistNEventsVsCent);

  fHistNEventsVsCL1 = new TH2F("hNEventsVsCL1", " ; ; N_{CL1}",21,-0.5,20.5,200,-0.5,2*maxMult-0.5);
  ConfigureEvSelAxis(fHistNEventsVsCL1->GetXaxis());
  fOutput->Add(fHistNEventsVsCL1);

  fHistNEventsVsTime = new TH1F("hNEventsVsTime", " ; Timestamp",44640,1541462400,1544140800);
  fOutput->Add(fHistNEventsVsTime);
  
  fHistWhyRej = new TH1F("hWhyRej"," ; WhyRej",11,-0.5,10.5);
  fOutput->Add(fHistWhyRej);
  fHistNEventsVsWhyRej = new TH2F("hNEventsVsWhyRej", " ; ; WhyRej ",21,-0.5,20.5,11,-0.5,10.5);
  ConfigureEvSelAxis(fHistNEventsVsWhyRej->GetXaxis());
  fOutput->Add(fHistNEventsVsWhyRej);

  fHistNTrackletsBeforePileup = new TH1F("hNTrackletsBeforePileup"," ; N_{tracklets}",200,-0.5,maxMult-0.5);
  fHistNTrackletsAfterPileup = new TH1F("hNTrackletsAfterPileup"," ; N_{tracklets}",200,-0.5,maxMult-0.5);
  fHistNCL1BeforePileup = new TH1F("hNCL1BeforePileup"," ; N_{CL1}",200,-0.5,maxMult-0.5);
  fHistNCL1AfterPileup = new TH1F("hNCL1AfterPileup"," ; N_{CL1}",200,-0.5,maxMult-0.5);
  fOutput->Add(fHistNTrackletsBeforePileup);
  fOutput->Add(fHistNTrackletsAfterPileup);
  fOutput->Add(fHistNCL1BeforePileup);
  fOutput->Add(fHistNCL1AfterPileup);

  fHistCentrality = new TH1F("hCentrality"," ; Centrality ; ",105,0.,105.);
  fHistCL0vsV0MCentrality = new TH2F("hCL0vsV0MCentrality"," ; Centrality V0M ; Centrality CL0",105,0.,105.,105,0.,105.);
  fOutput->Add(fHistCentrality);
  fOutput->Add(fHistCL0vsV0MCentrality);

  fHistNTracksTPCoutVsV0Cent = new TH2F("hNTracksTPCoutVsV0Cent"," ; Centrality ; N_{tracks, TPCout}",105,0.,105.,200,-0.5,2*maxMult-0.5);
  fHistNTracksFB4VsV0Cent = new TH2F("hNTracksFB4VsV0Cent"," ; Centrality ; N_{tracks, FiltBit4}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTracksFB4EtaPosVsV0Cent = new TH2F("hNTracksFB4EtaPosVsV0Cent"," ; Centrality ; N_{tracks, FiltBit4, |#eta|>0}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTracksFB4EtaNegVsV0Cent = new TH2F("hNTracksFB4EtaNegVsV0Cent"," ; Centrality ; N_{tracks, FiltBit4, |#eta|<0}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTracksBC0VsV0Cent = new TH2F("hNTracksBC0VsV0Cent"," ; Centrality ; N_{tracks, TOFBC=0}",105,0.,105.,200,-0.5,maxMult-0.5);  
  fHistNTrackletsVsV0Cent = new TH2F("hNTrackletsVsV0Cent"," ; Centrality ; N_{tracklets}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTrackletsGoldenVsV0Cent = new TH2F("hNTrackletsGoldenVsV0Cent"," ; Centrality ; N_{tracklets}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTrackletsGoldenVsV0CentVsZvert = new TH3F("hNTrackletsGoldenVsV0CentVsZvert"," ; Centrality ; N_{tracklets} ; z_{vertex} (cm)",105,0.,105.,200,-0.5,maxMult-0.5,30,-15.,15.);
  fHistNTracksTPCoutVsNTracklets = new TH2F("hNTracksTPCoutVsNTracklets"," ; N_{tracklets} ; N_{tracks, TPCout}",200,-0.5,maxMult-0.5,200,-0.5,2*maxMult-0.5);
  fHistNTracksFB4VsNTracklets = new TH2F("hNTracksFB4VsNTracklets"," ; N_{tracklets} ; N_{tracks, FiltBit4}",200,-0.5,maxMult-0.5,200,-0.5,maxMult-0.5);
  fHistNTracksBC0VsNTracksFB4 = new TH2F("hNTracksBC0VsNTracksFB4"," ; N_{tracks, FiltBit4}; N_{tracks, TOFBC=0}",200,-0.5,maxMult-0.5,200,-0.5,maxMult-0.5);
  fOutput->Add(fHistNTracksTPCoutVsV0Cent);
  fOutput->Add(fHistNTracksFB4VsV0Cent);
  fOutput->Add(fHistNTracksFB4EtaPosVsV0Cent);
  fOutput->Add(fHistNTracksFB4EtaNegVsV0Cent);
  fOutput->Add(fHistNTracksBC0VsV0Cent);
  fOutput->Add(fHistNTrackletsVsV0Cent);
  fOutput->Add(fHistNTrackletsGoldenVsV0Cent);
  fOutput->Add(fHistNTrackletsGoldenVsV0CentVsZvert);
  fOutput->Add(fHistNTracksTPCoutVsNTracklets);
  fOutput->Add(fHistNTracksFB4VsNTracklets);
  fOutput->Add(fHistNTracksBC0VsNTracksFB4);

  fHistPhiTrackelts = new TH1F("hPhiTrackelts"," ; #varphi (rad)",200,0.,2.*TMath::Pi());
  fOutput->Add(fHistPhiTrackelts);

  fHistNCL0VsV0Cent = new TH2F("hNCL0VsV0Cent"," ; Centrality ; N_{CL0}",105,0.,105.,200,-0.5,2*maxMult-0.5);
  fHistNCL1VsV0Cent = new TH2F("hNCL1VsV0Cent"," ; Centrality ; N_{CL1}",105,0.,105.,200,-0.5,2*maxMult-0.5);
  fHistT0AmplVsV0Ampl = new TH2F("hT0AmplVsV0Ampl"," ; V0 amplitude ; T0 amplitude",200,0.,60000,200,0.,2000.);
  fHistT0AmplVsV0Cent = new TH2F("hT0AmplVsV0Cent"," ; Centrality ; T0 amplitude",105,0.,105.,200,0.,2000.);
  fHistT0AmplVsNCL0 = new TH2F("hT0AmplVsNCL0"," ; N_{CL0} ; T0 amplitude",200,-0.5,2*maxMult-0.5,200,0.,2000.);
  fHistT0AmplVsCL0Cent = new TH2F("hT0AmplVsCL0Cent"," ; Centrality CL0 ; T0 amplitude",105,0.,105.,200,0.,2000.);
  fOutput->Add(fHistNCL0VsV0Cent);
  fOutput->Add(fHistNCL1VsV0Cent);
  fOutput->Add(fHistT0AmplVsV0Ampl);
  fOutput->Add(fHistT0AmplVsV0Cent);
  fOutput->Add(fHistT0AmplVsNCL0);
  fOutput->Add(fHistT0AmplVsCL0Cent);
  
  fHistZVertexSPDBeforeCuts = new TH2F("hZVertexSPDBeforeCuts"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDBeforeCuts);
  fHistZVertexSPDBeforeSPDCut = new TH2F("hZVertexSPDBeforeSPDCut"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDBeforeSPDCut);
  fHistZVertexSPDAfterCuts = new TH2F("hZVertexSPDAfterCuts"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDAfterCuts);
  fHistZVertexSPDBadTrackVert =new TH2F("hZVertexSPDBadTrackVert"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDBadTrackVert);

  const Int_t nVarForSp=8;
  Int_t nBinsForSp[nVarForSp]={105,30,200,200,200,200,200,200};
  Double_t minForSparse[nVarForSp]={0.,-15.,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5};
  Double_t maxForSparse[nVarForSp]={105.,15.,maxMult-0.5,maxMult-0.5,maxMult-0.5,maxMult-0.5,maxMult-0.5,maxMult-0.5};
  
  fEventProp = new THnSparseF("hEventProps","",nVarForSp,nBinsForSp,minForSparse,maxForSparse);
  fEventProp->GetAxis(0)->SetTitle("Centrality");
  fEventProp->GetAxis(1)->SetTitle("zSPDvertex (cm)");
  fEventProp->GetAxis(2)->SetTitle("nTracklets eta<1");
  fEventProp->GetAxis(3)->SetTitle("nTracklets golden SPD eta<1.4");
  fEventProp->GetAxis(4)->SetTitle("nTracklets golden SPD eta<1");
  fEventProp->GetAxis(5)->SetTitle("nGener eta<1");
  fEventProp->GetAxis(6)->SetTitle("nGener golden SPD eta<1.4");
  fEventProp->GetAxis(7)->SetTitle("nGener golden SPD eta<1");
  fOutput->Add(fEventProp);
  
  PostData(1,fOutput);

  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(2,fCounter);

  fNtupleEvProp = new TNtuple("fNtupleEvProp","ntupleEvProp","zSPDvert:centrV0:centrCL0:V0ampl:T0ampl:nTracklets:nClSPD0:nClSPD1:nTracksFB4");
  PostData(3,fNtupleEvProp);
  fNtupleZvtxDistVsWhyRej = new TNtuple("ntupleZvtxDistVsWhyRej","fNtupleZvtxDistVsWhyRej","zSPDvertex:zTRKvertex:NcontributorszSPD:NcontributorszTRK:whyrejection:vtxtype");
  PostData(4,fNtupleZvtxDistVsWhyRej);

}

//________________________________________________________________________
void AliAnalysisTaskCheckEvSel::ConfigureEvSelAxis(TAxis* ax){
  /// configure bin labels for event selection bits
  ///
  ax->SetBinLabel(1,"nEvents read");
  ax->SetBinLabel(2,"Rejected due to mismatch in trees");
  ax->SetBinLabel(3,"nEvents with good AOD");
  ax->SetBinLabel(4,"Rejected due to trigger");
  ax->SetBinLabel(5,"Rejected due to phys sel");
  ax->SetBinLabel(6,"Rejected due time range cut");
  ax->SetBinLabel(7,"Rejected due to bad centrality estimator");
  ax->SetBinLabel(8,"Rejected due to centrality flattening");
  ax->SetBinLabel(9,"Rejected due to centrality out of range");
  ax->SetBinLabel(10,"Rejected due to not reco vertex");
  ax->SetBinLabel(11,"Rejected for contr vertex");
  ax->SetBinLabel(12,"Rejected for bad track vertex");
  ax->SetBinLabel(13,"Rejected for vertex out of accept");
  ax->SetBinLabel(14,"Rejected for pileup events");
  ax->SetBinLabel(15,"Rejected due to centrality correlations");
  ax->SetBinLabel(16,"Rejected due to mult vs. V0 cut");
  ax->SetBinLabel(17,"Passing Phys sel + trigger");
  ax->SetBinLabel(18,"Passing Phys sel + trigger + Pileup");
  ax->SetBinLabel(19,"Passing Phys sel + trigger + Pileup + zVertex");
  ax->SetBinLabel(20,"Passing Phys sel + trigger + Pileup + zVertex + Centrality");
  ax->SetBinLabel(21,"Passing IsEventSelected");
  ax->SetNdivisions(1,kFALSE);
}
//________________________________________________________________________
void AliAnalysisTaskCheckEvSel::UserExec(Option_t */*option*/){
  //Build the 3-track combinatorics (+-+ and -+-) for D+->Kpipi decays
  
  if(fDebug > 1) printf("AnalysisTaskCheckEvSel::UserExec() \n");
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  fHistNEvents->Fill(0); // count event

  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(1);
      return;
    }
  }

  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if(!aod){
    printf("AliAnalysisTaskCheckEvSel::UserExec: AOD not found!\n");
    return;
  }
  
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  Int_t nGener=0;
  Int_t nGenerEta14=0;
  Int_t nGenerEta1=0;
  Int_t nGenerGoldenSPD=0;
  Int_t nGenerGoldenSPDEta14=0;
  Int_t nGenerGoldenSPDEta1=0;
  if(fReadMC && fMCEvent){
    for (int iMC = 0; iMC < fMCEvent->GetNumberOfTracks(); ++iMC) {
      AliVParticle *part = (AliVParticle*)fMCEvent->GetTrack(iMC);
      if(!part->IsPhysicalPrimary()) continue;
      if(part->Charge()==0) continue;
      nGener+=1;
      Double_t eta=part->Eta();
      Double_t phi=part->Phi();
      if(TMath::Abs(eta)<1.4) nGenerEta14+=1.;
      if(TMath::Abs(eta)<1) nGenerEta1+=1.;
      if(phi>3.9){
	nGenerGoldenSPD+=1.;
	if(TMath::Abs(eta)<1.4) nGenerGoldenSPDEta14+=1.;
	if(TMath::Abs(eta)<1) nGenerGoldenSPDEta1+=1.;
      }
    }
  }

  // Post the data already here

  Bool_t isEvSel=fAnalysisCuts->IsEventSelected(aod);
  if(fUseAliEventCuts) isEvSel=fAnalysisCuts->IsEventSelectedWithAliEventCuts(aod);
  Int_t wrej=fAnalysisCuts->GetWhyRejection();
  Double_t centr=fAnalysisCuts->GetCentrality(aod,AliRDHFCuts::kCentV0M);
  if(fSystem==2) centr=fAnalysisCuts->GetCentrality(aod,AliRDHFCuts::kCentZNA);// p-Pb
  Double_t centrCL0=fAnalysisCuts->GetCentrality(aod,AliRDHFCuts::kCentCL0);
  const AliVVertex *vertex = aod->GetPrimaryVertex();
  const AliVVertex *vertexSPD = aod->GetPrimaryVertexSPD();
  Double_t ntrkl = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.);
  Double_t nTrackletsGolden=0;
  Double_t nTrackletsGoldenEta1=0;
  Double_t nTrackletsGoldenEta14=0;
  Double_t nTrackletsEta1=0;
  Double_t nTrackletsEta14=0;
  AliAODTracklets* tracklets=aod->GetTracklets();
  Int_t nTr=tracklets->GetNumberOfTracklets();
  for(Int_t iTr=0; iTr<nTr; iTr++){
    Double_t theta=tracklets->GetTheta(iTr);
    Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
    Double_t phi=tracklets->GetPhi(iTr);
    fHistPhiTrackelts->Fill(phi);
    if(TMath::Abs(eta)<1.4) nTrackletsEta14+=1.;
    if(TMath::Abs(eta)<1) nTrackletsEta1+=1.;
    if(phi>3.9){
      nTrackletsGolden+=1.;
      if(TMath::Abs(eta)<1.4) nTrackletsGoldenEta14+=1.;
      if(TMath::Abs(eta)<1) nTrackletsGoldenEta1+=1.;
    }
  }

  Double_t ncl0 = aod->GetNumberOfITSClusters(0);
  Double_t ncl1 = aod->GetNumberOfITSClusters(1);
  Double_t zvSPD=vertexSPD->GetZ();
  Double_t zvTRK=vertex->GetZ();
  AliAODVZERO* vzer=(AliAODVZERO*)aod->GetVZEROData();
  Float_t v0ampl=vzer->GetMTotV0A()+vzer->GetMTotV0C();
  AliAODTZERO* tz=aod->GetTZEROData();
  Double_t t0ampl=0;
  for(Int_t j=0; j<24; j++) t0ampl+=tz->GetAmp(j);

  if(!isEvSel)fHistWhyRej->Fill(wrej); 
    
  if(fEnableVertexNtuple) {
    Float_t wrej4ntuple = (Float_t)wrej;
    if(isEvSel) wrej4ntuple = -1;
    Float_t vertextype=0;
    if(vertex) {
      TString vertextitle = vertex->GetTitle();
      //Vtx tracks
      if(vertextitle.Contains("VertexerTracks")) vertextype=3;
      //vtx SPD
      else if(vertextitle.Contains("Z")) vertextype=1;
      else if(vertextitle.Contains("3D")) vertextype=2;
    }
      
    Float_t vec4ntuple[6] = {(Float_t)zvSPD,(Float_t)zvTRK,(Float_t)vertexSPD->GetNContributors(),(Float_t)vertex->GetNContributors(),wrej4ntuple,vertextype};
    fNtupleZvtxDistVsWhyRej->Fill(vec4ntuple);
    PostData(4,fNtupleZvtxDistVsWhyRej);
  }
  if(fAnalysisCuts->IsEventRejectedDueToTrigger()) fHistNEventsVsWhyRej->Fill(3,wrej);
  if(fAnalysisCuts->IsEventRejectedDuePhysicsSelection()) fHistNEventsVsWhyRej->Fill(4,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToTimeRangeCut()) fHistNEventsVsWhyRej->Fill(5,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality()) fHistNEventsVsWhyRej->Fill(6,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToCentralityFlattening()) fHistNEventsVsWhyRej->Fill(7,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex()) fHistNEventsVsWhyRej->Fill(9,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors()) fHistNEventsVsWhyRej->Fill(10,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToBadTrackVertex()) fHistNEventsVsWhyRej->Fill(11,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) fHistNEventsVsWhyRej->Fill(12,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToPileup()) fHistNEventsVsWhyRej->Fill(13,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToCentralityEstimCorrel()) fHistNEventsVsWhyRej->Fill(14,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToTRKV0CentralityCorrel()) fHistNEventsVsWhyRej->Fill(15,wrej);

  fHistZVertexSPDBeforeCuts->Fill(zvSPD,zvTRK);
  
  if(!fAnalysisCuts->IsEventRejectedDueToTrigger() &&
     !fAnalysisCuts->IsEventRejectedDuePhysicsSelection() &&
     !fAnalysisCuts->IsEventRejectedDueToNotRecoVertex() &&
     !fAnalysisCuts->IsEventRejectedDueToVertexContributors() &&
     vertexSPD->GetNContributors()>=1){
    fHistZVertexSPDBeforeSPDCut->Fill(zvSPD,zvTRK);
    Double_t dz = vertexSPD->GetZ()-vertex->GetZ();
    Bool_t okSpdTrk=kTRUE;
    if(fCutOnzVertexSPD==2 && TMath::Abs(dz)>0.5) okSpdTrk=kFALSE;
    if(okSpdTrk && fCutOnzVertexSPD==3){
      double covTrc[6],covSPD[6];
      vertex->GetCovarianceMatrix(covTrc);
      vertexSPD->GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) okSpdTrk=kFALSE;
    }
    if(!okSpdTrk) fHistZVertexSPDBadTrackVert->Fill(zvSPD,zvTRK);
  }
  if(isEvSel) fHistZVertexSPDAfterCuts->Fill(zvSPD,zvTRK);

  fHistNEvents->Fill(2);
  fHistNEventsVsCent->Fill(2,centr);
  fHistNEventsVsCL1->Fill(2,ncl1);

  Int_t binToFill=-1;
  if(fAnalysisCuts->IsEventRejectedDueToTrigger() || fAnalysisCuts->IsEventRejectedDuePhysicsSelection() || fAnalysisCuts->IsEventRejectedDueToTimeRangeCut() ){
    if(fAnalysisCuts->IsEventRejectedDueToTrigger()) binToFill=3;
    if(fAnalysisCuts->IsEventRejectedDuePhysicsSelection()) binToFill=4;
    if(fAnalysisCuts->IsEventRejectedDueToTimeRangeCut()) binToFill=5;
  }else{
    if(fAnalysisCuts->IsEventRejectedDueToCentrality() || fAnalysisCuts->IsEventRejectedDueToCentralityFlattening()){
      if(wrej==3) binToFill=6;
      else if(wrej==4) binToFill=7;
      else binToFill=8;
    }else{
      if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex() || fAnalysisCuts->IsEventRejectedDueToVertexContributors() || fAnalysisCuts->IsEventRejectedDueToBadTrackVertex()){
	if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex()) binToFill=9;
	if(fAnalysisCuts->IsEventRejectedDueToVertexContributors()) binToFill=10;
	if(fAnalysisCuts->IsEventRejectedDueToBadTrackVertex()) binToFill=11;
      }else{
	if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) binToFill=12;
	else if(fAnalysisCuts->IsEventRejectedDueToPileup()) binToFill=13;
	else if(fAnalysisCuts->IsEventRejectedDueToCentralityEstimCorrel()) binToFill=14;
	else if(fAnalysisCuts->IsEventRejectedDueToTRKV0CentralityCorrel())  binToFill=15;
      }
    }
  }
  if(binToFill>0){
    fHistNEvents->Fill(binToFill);
    fHistNEventsVsCent->Fill(binToFill,centr);
    fHistNEventsVsCL1->Fill(binToFill,ncl1);
  }
  if(fAnalysisCuts->IsEventRejectedDueToTrigger()==0 && fAnalysisCuts->IsEventRejectedDuePhysicsSelection()==0){
    fHistNEvents->Fill(16);
    fHistNEventsVsCent->Fill(16,centr);
    fHistNEventsVsCL1->Fill(16,ncl1);
    fHistNEventsVsWhyRej->Fill(16,wrej);    
    if(fAnalysisCuts->IsEventRejectedDueToPileup()==0){
      fHistNEvents->Fill(17);
      fHistNEventsVsCent->Fill(17,centr);
      fHistNEventsVsCL1->Fill(17,ncl1);
      fHistNEventsVsWhyRej->Fill(17,wrej);
      if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()==0){
	fHistNEvents->Fill(18);
	fHistNEventsVsCent->Fill(18,centr);
	fHistNEventsVsCL1->Fill(18,ncl1);
	fHistNEventsVsWhyRej->Fill(18,wrej);
	if(fAnalysisCuts->IsEventRejectedDueToCentrality()==0 && fAnalysisCuts->IsEventRejectedDueToCentralityFlattening()==0){
	  fHistNEvents->Fill(19);
	  fHistNEventsVsCent->Fill(19,centr);
	  fHistNEventsVsCL1->Fill(19,ncl1);
	  fHistNEventsVsWhyRej->Fill(19,wrej);
	}
      }
    }
  }
  if(isEvSel || (!isEvSel && wrej==1)){
    fHistNTrackletsBeforePileup->Fill(ntrkl);
    fHistNCL1BeforePileup->Fill(ncl1);
  }
  if(isEvSel){
    fHistNEvents->Fill(20);
    fHistNEventsVsCent->Fill(20,centr);
    fHistNEventsVsCL1->Fill(20,ncl1);
    fHistNEventsVsWhyRej->Fill(20,wrej);
    fHistNTrackletsAfterPileup->Fill(ntrkl);
    fHistNCL1AfterPileup->Fill(ncl1);
    Int_t runNumb=aod->GetRunNumber();
    if(fAnalysisCuts->GetUseTimeRangeCutForPbPb2018() && runNumb >= 295369 && runNumb <= 297624){
      fHistNEventsVsTime->Fill(aod->GetTimeStamp());
    }
  }

  if(fAnalysisCuts->GetUseCentrality()>0 && fAnalysisCuts->IsEventSelectedInCentrality(aod)!=0) return;
  // events not passing the centrality selection can be removed immediately. For the others we must count the generated D mesons
  fCounter->StoreEvent(aod,fAnalysisCuts,fReadMC);

  if(isEvSel){

    Int_t ntracksTPCout=0;
    Int_t ntracksFB4=0;
    Int_t ntracksFB4EtaPos=0;
    Int_t ntracksFB4EtaNeg=0;
    Int_t ntracksBC0=0;
    
    Double_t magField  = aod->GetMagneticField();

    for(Int_t iTr=0; iTr<aod->GetNumberOfTracks(); iTr++){
      AliAODTrack* track=dynamic_cast<AliAODTrack*>(aod->GetTrack(iTr));
      
      if(track->GetStatus() & AliESDtrack::kTPCout) ntracksTPCout++;
      if(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){
	ntracksFB4++;
	Double_t etatrack=track->Eta();
	if(etatrack>0) ntracksFB4EtaPos++;
	else ntracksFB4EtaNeg++;
	Int_t tofBC=track->GetTOFBunchCrossing(magField);
	if(tofBC==0) ntracksBC0++;
      }
    }
    fHistCentrality->Fill(centr);
    fHistCL0vsV0MCentrality->Fill(centr,centrCL0);
    fHistNTracksTPCoutVsV0Cent->Fill(centr,ntracksTPCout);
    fHistNTracksFB4VsV0Cent->Fill(centr,ntracksFB4);
    fHistNTracksFB4EtaPosVsV0Cent->Fill(centr,ntracksFB4EtaPos);
    fHistNTracksFB4EtaNegVsV0Cent->Fill(centr,ntracksFB4EtaNeg);
    fHistNTracksBC0VsV0Cent->Fill(centr,ntracksBC0);
    fHistNTrackletsVsV0Cent->Fill(centr,ntrkl);
    fHistNTrackletsGoldenVsV0Cent->Fill(centr,nTrackletsGolden);
    fHistNTrackletsGoldenVsV0CentVsZvert->Fill(centr,nTrackletsGolden,zvSPD);
    fHistNTracksTPCoutVsNTracklets->Fill(ntracksTPCout,ntrkl);
    fHistNTracksFB4VsNTracklets->Fill(ntracksFB4,ntrkl);
    fHistNTracksBC0VsNTracksFB4->Fill(ntracksBC0,ntracksFB4);
    fHistNCL0VsV0Cent->Fill(centr,ncl0);
    fHistNCL1VsV0Cent->Fill(centr,ncl1);
    fHistT0AmplVsV0Cent->Fill(centr,t0ampl);
    fHistT0AmplVsV0Ampl->Fill(v0ampl,t0ampl);
    fHistT0AmplVsNCL0->Fill(ncl0,t0ampl);
    fHistT0AmplVsCL0Cent->Fill(centrCL0,t0ampl);
    Double_t vec4Sp[8];
    vec4Sp[0]=centr;
    vec4Sp[1]=zvSPD;
    vec4Sp[2]=nTrackletsEta1;
    vec4Sp[3]=nTrackletsGoldenEta14;
    vec4Sp[4]=nTrackletsGoldenEta1;
    vec4Sp[5]=nGenerEta1;
    vec4Sp[6]=nGenerGoldenSPDEta14;
    vec4Sp[7]=nGenerGoldenSPDEta1;
    fEventProp->Fill(vec4Sp);
    if(fEnableEvPropNtuple){
      Float_t vec4ep[9]={(Float_t)zvSPD,(Float_t)centr,(Float_t)centrCL0,
			 (Float_t)v0ampl,(Float_t)t0ampl,
			 (Float_t)ntrkl,(Float_t)ncl0,(Float_t)ncl1,
			 (Float_t)ntracksFB4};
      fNtupleEvProp->Fill(vec4ep);
      PostData(3,fNtupleEvProp);
    }
  }

  PostData(1,fOutput);
  PostData(2,fCounter);
  
  return;
}

//_________________________________________________________________
void AliAnalysisTaskCheckEvSel::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AliAnalysisTaskCheckEvSel: Terminate() \n");
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  if(fHistNEvents){
    printf("Number of analyzed events = %d\n",(Int_t)fHistNEvents->GetBinContent(2));
  }else{
    printf("ERROR: fHistNEvents not available\n");
    return;
  }
  return;
}

