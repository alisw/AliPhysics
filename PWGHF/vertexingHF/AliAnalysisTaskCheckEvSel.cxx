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

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
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
  fHistNEventsVsWhyRej(0x0),
  fHistNTrackletsBeforePileup(0x0),
  fHistNTrackletsAfterPileup(0x0),
  fHistNCL1BeforePileup(0x0),
  fHistNCL1AfterPileup(0x0),
  fHistCentrality(0x0),
  fHistNTracksTPCoutVsV0Cent(0x0),
  fHistNTracksFB4VsV0Cent(0x0),
  fHistNTracksBC0VsV0Cent(0x0),
  fHistNTrackletsVsV0Cent(0x0),
  fHistNTracksTPCoutVsNTracklets(0x0),
  fHistNTracksFB4VsNTracklets(0x0),
  fHistNTracksBC0VsNTracksFB4(0x0),
  fHistZVertexSPDBeforeCuts(0x0),
  fHistZVertexSPDBeforeSPDCut(0x0),
  fHistZVertexSPDAfterCuts(0x0),
  fHistZVertexSPDBadTrackVert(0x0),
  fNtupleZvtxDistVsWhyRej(0x0),
  fEnableVertexNtuple(kFALSE),
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
  fHistNEventsVsWhyRej(0x0),
  fHistNTrackletsBeforePileup(0x0),
  fHistNTrackletsAfterPileup(0x0),
  fHistNCL1BeforePileup(0x0),
  fHistNCL1AfterPileup(0x0),
  fHistCentrality(0x0),
  fHistNTracksTPCoutVsV0Cent(0x0),
  fHistNTracksFB4VsV0Cent(0x0),
  fHistNTracksBC0VsV0Cent(0x0),
  fHistNTrackletsVsV0Cent(0x0),
  fHistNTracksTPCoutVsNTracklets(0x0),
  fHistNTracksFB4VsNTracklets(0x0),
  fHistNTracksBC0VsNTracksFB4(0x0),
  fHistZVertexSPDBeforeCuts(0x0),
  fHistZVertexSPDBeforeSPDCut(0x0),
  fHistZVertexSPDAfterCuts(0x0),
  fHistZVertexSPDBadTrackVert(0x0),
  fNtupleZvtxDistVsWhyRej(0x0),
  fEnableVertexNtuple(kFALSE),
  fCounter(0),
  fAnalysisCuts(cuts)
{
  // default constructor

  DefineOutput(1,TList::Class());  //My private output
  DefineOutput(2,AliNormalizationCounter::Class());
  if(fEnableVertexNtuple){
    // Output slot #3 writes into a TNtuple container
    DefineOutput(3,TNtuple::Class());  //My private output
  }
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
    delete fHistNEventsVsWhyRej;
    delete fHistNTrackletsBeforePileup;
    delete fHistNTrackletsAfterPileup;
    delete fHistNCL1BeforePileup;
    delete fHistNCL1AfterPileup;
    delete fHistCentrality;
    delete fHistNTracksTPCoutVsV0Cent;
    delete fHistNTracksFB4VsV0Cent;
    delete fHistNTracksBC0VsV0Cent;
    delete fHistNTrackletsVsV0Cent;
    delete fHistNTracksTPCoutVsNTracklets;
    delete fHistNTracksFB4VsNTracklets;
    delete fHistNTracksBC0VsNTracksFB4;
    delete fHistZVertexSPDBeforeCuts;
    delete fHistZVertexSPDBeforeSPDCut;
    delete fHistZVertexSPDAfterCuts;  
    delete fHistZVertexSPDBadTrackVert;
  }
  delete fOutput;
  delete fCounter;
  if(fNtupleZvtxDistVsWhyRej) delete fNtupleZvtxDistVsWhyRej;
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
  
  fHistNEvents = new TH1F("hNEvents", "number of events ",18,-0.5,17.5);
  ConfigureEvSelAxis(fHistNEvents->GetXaxis());
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  fHistNEventsVsCent = new TH2F("hNEventsVsCent", " ; ; Centrality ",18,-0.5,17.5,101,0.,101.);
  ConfigureEvSelAxis(fHistNEventsVsCent->GetXaxis());
  fOutput->Add(fHistNEventsVsCent);

  fHistNEventsVsCL1 = new TH2F("hNEventsVsCL1", " ; ; N_{CL1}",18,-0.5,17.5,200,-0.5,2*maxMult-0.5);
  ConfigureEvSelAxis(fHistNEventsVsCL1->GetXaxis());
  fOutput->Add(fHistNEventsVsCL1);

  fHistNEventsVsWhyRej = new TH2F("hNEventsVsWhyRej", " ; ; WhyRej ",18,-0.5,17.5,11,-0.5,10.5);
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
  fOutput->Add(fHistCentrality);

  fHistNTracksTPCoutVsV0Cent = new TH2F("hNTracksTPCoutVsV0Cent"," ; Centrality ; N_{tracks, TPCout}",105,0.,105.,200,-0.5,2*maxMult-0.5);
  fHistNTracksFB4VsV0Cent = new TH2F("hNTracksFB4VsV0Cent"," ; Centrality ; N_{tracks, FiltBit4}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTracksBC0VsV0Cent = new TH2F("hNTracksBC0VsV0Cent"," ; Centrality ; N_{tracks, TOFBC=0}",105,0.,105.,200,-0.5,maxMult-0.5);  
  fHistNTrackletsVsV0Cent = new TH2F("hNTrackletsVsV0Cent"," ; Centrality ; N_{tracklets}",105,0.,105.,200,-0.5,maxMult-0.5);
  fHistNTracksTPCoutVsNTracklets = new TH2F("hNTracksTPCoutVsNTracklets"," ; N_{tracklets} ; N_{tracks, TPCout}",200,-0.5,maxMult-0.5,200,-0.5,2*maxMult-0.5);
  fHistNTracksFB4VsNTracklets = new TH2F("hNTracksFB4VsNTracklets"," ; N_{tracklets} ; N_{tracks, FiltBit4}",200,-0.5,maxMult-0.5,200,-0.5,maxMult-0.5);
  fHistNTracksBC0VsNTracksFB4 = new TH2F("hNTracksBC0VsNTracksFB4"," ; N_{tracks, FiltBit4}; N_{tracks, TOFBC=0}",200,-0.5,maxMult-0.5,200,-0.5,maxMult-0.5);
  fOutput->Add(fHistNTracksTPCoutVsV0Cent);
  fOutput->Add(fHistNTracksFB4VsV0Cent);
  fOutput->Add(fHistNTracksBC0VsV0Cent);
  fOutput->Add(fHistNTrackletsVsV0Cent);
  fOutput->Add(fHistNTracksTPCoutVsNTracklets);
  fOutput->Add(fHistNTracksFB4VsNTracklets);
  fOutput->Add(fHistNTracksBC0VsNTracksFB4);

  fHistZVertexSPDBeforeCuts = new TH2F("hZVertexSPDBeforeCuts"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDBeforeCuts);
  fHistZVertexSPDBeforeSPDCut = new TH2F("hZVertexSPDBeforeSPDCut"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDBeforeSPDCut);
  fHistZVertexSPDAfterCuts = new TH2F("hZVertexSPDAfterCuts"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDAfterCuts);
  fHistZVertexSPDBadTrackVert =new TH2F("hZVertexSPDBadTrackVert"," ; z_{SPDvertex} ; z_{TRKvertex}",400,-20.,20.,400,-20.,20.);
  fOutput->Add(fHistZVertexSPDBadTrackVert);
  
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();

  PostData(1,fOutput);
  
  if(fHistZVertexSPDBadTrackVert) {
    OpenFile(3); // 3 is the slot number of the ntuple
    
    fNtupleZvtxDistVsWhyRej = new TNtuple("fNtupleZvtxDistVsWhyRej","fNtupleZvtxDistVsWhyRej","zSPDvertex:zTRKvertex:NcontributorszSPD:NcontributorszTRK:whyrejection:vtxtype");
  }
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
  ax->SetBinLabel(6,"Rejected due to bad centrality estimator");
  ax->SetBinLabel(7,"Rejected due to centrality flattening");
  ax->SetBinLabel(8,"Rejected due to centrality out of range");
  ax->SetBinLabel(9,"Rejected due to not reco vertex");
  ax->SetBinLabel(10,"Rejected for contr vertex");
  ax->SetBinLabel(11,"Rejected for bad track vertex");
  ax->SetBinLabel(12,"Rejected for vertex out of accept");
  ax->SetBinLabel(13,"Rejected for pileup events");
  ax->SetBinLabel(14,"Passing Phys sel + trigger");
  ax->SetBinLabel(15,"Passing Phys sel + trigger + Pileup");
  ax->SetBinLabel(16,"Passing Phys sel + trigger + Pileup + zVertex");
  ax->SetBinLabel(17,"Passing Phys sel + trigger + Pileup + zVertex + Centrality");
  ax->SetBinLabel(18,"Passing IsEventSelected");
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
  

  // Post the data already here

  Bool_t isEvSel=fAnalysisCuts->IsEventSelected(aod);
  Int_t wrej=fAnalysisCuts->GetWhyRejection();
  Double_t centr=fAnalysisCuts->GetCentrality(aod,AliRDHFCuts::kCentZNA);
  const AliVVertex *vertex = aod->GetPrimaryVertex();
  const AliVVertex *vertexSPD = aod->GetPrimaryVertexSPD();
  Double_t ntrkl = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.); 
  Double_t ncl1 = aod->GetNumberOfITSClusters(1);
  Double_t zvSPD=vertexSPD->GetZ();
  Double_t zvTRK=vertex->GetZ();

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
    PostData(3,fNtupleZvtxDistVsWhyRej);
  }
  
  if(fAnalysisCuts->IsEventRejectedDueToTrigger()) fHistNEventsVsWhyRej->Fill(3,wrej);
  if(fAnalysisCuts->IsEventRejectedDuePhysicsSelection()) fHistNEventsVsWhyRej->Fill(4,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality()) fHistNEventsVsWhyRej->Fill(5,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToCentralityFlattening()) fHistNEventsVsWhyRej->Fill(6,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex()) fHistNEventsVsWhyRej->Fill(8,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors()) fHistNEventsVsWhyRej->Fill(9,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToBadTrackVertex()) fHistNEventsVsWhyRej->Fill(10,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) fHistNEventsVsWhyRej->Fill(11,wrej);
  if(fAnalysisCuts->IsEventRejectedDueToPileup()) fHistNEventsVsWhyRej->Fill(12,wrej);

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
  if(fAnalysisCuts->IsEventRejectedDueToTrigger() || fAnalysisCuts->IsEventRejectedDuePhysicsSelection()){
    if(fAnalysisCuts->IsEventRejectedDueToTrigger()) binToFill=3;
    if(fAnalysisCuts->IsEventRejectedDuePhysicsSelection()) binToFill=4;
  }else{
    if(fAnalysisCuts->IsEventRejectedDueToCentrality() || fAnalysisCuts->IsEventRejectedDueToCentralityFlattening()){
      if(wrej==3) binToFill=5;
      else if(wrej==4) binToFill=6;
      else binToFill=7;
    }else{
      if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex() || fAnalysisCuts->IsEventRejectedDueToVertexContributors() || fAnalysisCuts->IsEventRejectedDueToBadTrackVertex()){
	if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex()) binToFill=9;
	if(fAnalysisCuts->IsEventRejectedDueToVertexContributors()) binToFill=9;
	if(fAnalysisCuts->IsEventRejectedDueToBadTrackVertex()) binToFill=10;
      }else{
	if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) binToFill=11;
	else if(fAnalysisCuts->IsEventRejectedDueToPileup()) binToFill=12;
      }
    }
  }
  if(binToFill>0){
    fHistNEvents->Fill(binToFill);
    fHistNEventsVsCent->Fill(binToFill,centr);
    fHistNEventsVsCL1->Fill(binToFill,ncl1);
  }
  if(fAnalysisCuts->IsEventRejectedDueToTrigger()==0 && fAnalysisCuts->IsEventRejectedDuePhysicsSelection()==0){
    fHistNEvents->Fill(13);
    fHistNEventsVsCent->Fill(13,centr);
    fHistNEventsVsCL1->Fill(13,ncl1);
    fHistNEventsVsWhyRej->Fill(13,wrej);    
    if(fAnalysisCuts->IsEventRejectedDueToPileup()==0){
      fHistNEvents->Fill(14);
      fHistNEventsVsCent->Fill(14,centr);
      fHistNEventsVsCL1->Fill(14,ncl1);
      fHistNEventsVsWhyRej->Fill(14,wrej);
      if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()==0){
	fHistNEvents->Fill(15);
	fHistNEventsVsCent->Fill(15,centr);
	fHistNEventsVsCL1->Fill(15,ncl1);
	fHistNEventsVsWhyRej->Fill(15,wrej);
	if(fAnalysisCuts->IsEventRejectedDueToCentrality()==0 && fAnalysisCuts->IsEventRejectedDueToCentralityFlattening()==0){
	  fHistNEvents->Fill(16);
	  fHistNEventsVsCent->Fill(16,centr);
	  fHistNEventsVsCL1->Fill(16,ncl1);
	  fHistNEventsVsWhyRej->Fill(16,wrej);
	}
      }
    }
  }
  if(isEvSel || (!isEvSel && wrej==1)){
    fHistNTrackletsBeforePileup->Fill(ntrkl);
    fHistNCL1BeforePileup->Fill(ncl1);
  }
  if(isEvSel){
    fHistNEvents->Fill(17);
    fHistNEventsVsCent->Fill(17,centr);
    fHistNEventsVsCL1->Fill(17,ncl1);
    fHistNEventsVsWhyRej->Fill(17,wrej);
    fHistNTrackletsAfterPileup->Fill(ntrkl);
    fHistNCL1AfterPileup->Fill(ncl1);
  }

  if(fAnalysisCuts->GetUseCentrality()>0 && fAnalysisCuts->IsEventSelectedInCentrality(aod)!=0) return;
  // events not passing the centrality selection can be removed immediately. For the others we must count the generated D mesons
  fCounter->StoreEvent(aod,fAnalysisCuts,fReadMC);

  if(isEvSel){

    Int_t ntracksTPCout=0;
    Int_t ntracksFB4=0;
    Int_t ntracksBC0=0;
    
    Double_t magField  = aod->GetMagneticField();

    for(Int_t iTr=0; iTr<aod->GetNumberOfTracks(); iTr++){
      AliAODTrack* track=dynamic_cast<AliAODTrack*>(aod->GetTrack(iTr));
      
      if(track->GetStatus() & AliESDtrack::kTPCout) ntracksTPCout++;
      if(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){
	ntracksFB4++;
	Int_t tofBC=track->GetTOFBunchCrossing(magField);
	if(tofBC==0) ntracksBC0++;
      }
    }
    fHistCentrality->Fill(centr);
    fHistNTracksTPCoutVsV0Cent->Fill(centr,ntracksTPCout);
    fHistNTracksFB4VsV0Cent->Fill(centr,ntracksFB4);
    fHistNTracksBC0VsV0Cent->Fill(centr,ntracksBC0);
    fHistNTrackletsVsV0Cent->Fill(centr,ntrkl);
    fHistNTracksTPCoutVsNTracklets->Fill(ntracksTPCout,ntrkl);
    fHistNTracksFB4VsNTracklets->Fill(ntracksFB4,ntrkl);
    fHistNTracksBC0VsNTracksFB4->Fill(ntracksBC0,ntracksFB4);
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

