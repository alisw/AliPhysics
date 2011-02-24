/**************************************************************************
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

/* $Id$ */

//*************************************************************************
// Class AliNormalizationCounter
// Class to store the informations relevant for the normalization in the 
// barrel for each run
// Authors: G. Ortona, ortona@to.infn.it
// D. Caffarri, davide.caffarri@pd.to.infn.it
// with many thanks to P. Pillot
/////////////////////////////////////////////////////////////

#include "AliNormalizationCounter.h"
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODEvent.h>
#include <AliVParticle.h>
#include <AliTriggerAnalysis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TString.h>
#include <TCanvas.h>
#include <AliPhysicsSelection.h>
#include <AliMultiplicity.h>

ClassImp(AliNormalizationCounter)

//____________________________________________
AliNormalizationCounter::AliNormalizationCounter(): 
TNamed(),
fCounters(),
fESD(kFALSE),
fRejectPileUp(kFALSE),
fHistTrackFilterEvMult(0),
fHistTrackAnaEvMult(0),
fHistTrackFilterSpdMult(0),
fHistTrackAnaSpdMult(0)
{
  // empty constructor
}

//__________________________________________________				
AliNormalizationCounter::AliNormalizationCounter(const char *name): 
TNamed(name,name),
fCounters(name),
fESD(kFALSE),
fRejectPileUp(kFALSE),
fHistTrackFilterEvMult(0),
fHistTrackAnaEvMult(0),
fHistTrackFilterSpdMult(0),
fHistTrackAnaSpdMult(0)
{
  //default constructor
  fCounters.AddRubric("Event","triggered/V0AND/PileUp/PbPbC0SMH-B-NOPF-ALLNOTRD/Candles0.2/Candles0.2spd1/Candles0.3/2xCandles0.2/CandleITSsa/Candle35clsTPC/PrimaryVTracks/PrimaryV/PrimaryVSPD/!V0A&Candle02/!V0A&Candle025/!V0A&Candle03/!V0A&PrimaryVTracks/!V0A&PrimaryV/!V0A&2xCandles02/Candid(Filter)/Candid(Analysis)/NCandid(Filter)/NCandid(Analysis)");
  fCounters.AddRubric("Run", 1000000);
  fCounters.Init();
  fHistTrackFilterEvMult=new TH2F("FiltCandidvsTracksinEv","FiltCandidvsTracksinEv",10000,-0.5,9999.5,200,-0.5,199.5);
  fHistTrackFilterEvMult->GetYaxis()->SetTitle("NCandidates");
  fHistTrackFilterEvMult->GetXaxis()->SetTitle("NTracksinEvent");
  fHistTrackAnaEvMult=new TH2F("AnaCandidvsTracksinEv","AnaCandidvsTracksinEv",10000,-0.5,9999.5,100,-0.5,99.5);
  fHistTrackAnaEvMult->GetYaxis()->SetTitle("NCandidates");
  fHistTrackAnaEvMult->GetXaxis()->SetTitle("NTracksinEvent");
  fHistTrackFilterSpdMult=new TH2F("FilterCandidvsSpdMult","FilterCandidvsSpdMult",5000,-0.5,4999.5,200,-0.5,199.5);
  fHistTrackFilterSpdMult->GetYaxis()->SetTitle("NCandidates");
  fHistTrackFilterSpdMult->GetXaxis()->SetTitle("NSPDTracklets");
  fHistTrackAnaSpdMult=new TH2F("AnaCandidvsSpdMult","AnaCandidvsSpdMult",5000,-0.5,4999.5,100,-0.5,99.5);
  fHistTrackAnaSpdMult->GetYaxis()->SetTitle("NCandidates");
  fHistTrackAnaSpdMult->GetXaxis()->SetTitle("NSPDTracklets");
}
//______________________________________________
AliNormalizationCounter::~AliNormalizationCounter()
{
  //destructor
  if(fHistTrackFilterEvMult){
    delete fHistTrackFilterEvMult;
    fHistTrackFilterEvMult =0;
  }
  if(fHistTrackAnaEvMult){
    delete fHistTrackAnaEvMult;
    fHistTrackAnaEvMult=0;
  }
  if(fHistTrackFilterSpdMult){
    delete fHistTrackFilterSpdMult;
    fHistTrackFilterSpdMult=0;
  }
  if(fHistTrackAnaSpdMult){
    delete fHistTrackAnaSpdMult;
    fHistTrackAnaSpdMult=0;
  }
}

//______________________________________________
Long64_t AliNormalizationCounter::Merge(TCollection* list){
  if (!list) return 0;
  if (list->IsEmpty()) return 0;//(Long64_t)fCounters.Merge(list);

  TIter next(list);
  const TObject* obj = 0x0;
  while ((obj = next())) {
    
    // check that "obj" is an object of the class AliNormalizationCounter
    const AliNormalizationCounter* counter = dynamic_cast<const AliNormalizationCounter*>(obj);
    if (!counter) {
      AliError(Form("object named %s is not AliNormalizationCounter! Skipping it.", counter->GetName()));
      continue;
    }

    Add(counter);

  }
  
  return (Long64_t)1;//(Long64_t)fCounters->GetEntries();
}
//_______________________________________
void AliNormalizationCounter::Add(const AliNormalizationCounter *norm){
  fCounters.Add(&(norm->fCounters));
  fHistTrackFilterEvMult->Add(norm->fHistTrackFilterEvMult);
  fHistTrackAnaEvMult->Add(norm->fHistTrackAnaEvMult);
  fHistTrackFilterSpdMult->Add(norm->fHistTrackFilterSpdMult);
  fHistTrackAnaSpdMult->Add(norm->fHistTrackAnaSpdMult);
}
//_______________________________________
void AliNormalizationCounter::StoreEvent(AliVEvent *event,Bool_t mc){
  //

  Bool_t v0A=kFALSE; 
  Bool_t v0B=kFALSE;
  Bool_t flag02=kFALSE;
  Bool_t flag03=kFALSE;
  Int_t flag0202=0;
  Bool_t flagPV=kFALSE;
  Bool_t flagPVT=kFALSE; 
  Bool_t flag35cls=kFALSE;
  Bool_t flagITSsa=kFALSE;
  Bool_t flag02spd=kFALSE;

  //Run Number
  Int_t runNumber = event->GetRunNumber();
 
  //Find CINT1B
  AliESDEvent *eventESD = (AliESDEvent*)event;
  if(!eventESD){AliError("ESD event not available");return;}
  if(mc&&event->GetEventType() != 0)return;
  //event must be either physics or MC
  if(!(event->GetEventType() == 7||event->GetEventType() == 0))return;
  
  fCounters.Count(Form("Event:triggered/Run:%d",runNumber));
      
  //Find V0AND
  AliTriggerAnalysis trAn; /// Trigger Analysis
  v0B = trAn.IsOfflineTriggerFired(eventESD , AliTriggerAnalysis::kV0C);
  v0A = trAn.IsOfflineTriggerFired(eventESD , AliTriggerAnalysis::kV0A);
  if(v0A&&v0B){fCounters.Count(Form("Event:V0AND/Run:%d",runNumber));}
  
 //FindPrimary vertex  
  AliAODEvent *eventAOD = (AliAODEvent*)event;
  AliVVertex *vtrc =  (AliVVertex*)event->GetPrimaryVertex();
  if(vtrc && vtrc->GetNContributors()>0){
    fCounters.Count(Form("Event:PrimaryV/Run:%d",runNumber));
    flagPV=kTRUE;
  }

  AliAODVertex *vrtcSPD = (AliAODVertex*)eventAOD->GetPrimaryVertexSPD();
  if(vrtcSPD){
    if(vrtcSPD->GetNContributors()>0)fCounters.Count(Form("Event:PrimaryVSPD/Run:%d",runNumber));
  }
  
  if(fESD){
    const AliESDVertex *vtrc1 =  eventESD->GetPrimaryVertexTracks();
    if(vtrc1 && vtrc1->GetNContributors()>0){
      fCounters.Count(Form("Event:PrimaryVTracks/Run:%d",runNumber));
      flagPVT=kTRUE;
    }
  }

  //trigger
  TString trigclass=eventAOD->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL"))fCounters.Count(Form("Event:PbPbC0SMH-B-NOPF-ALLNOTRD/Run:%d",runNumber));

  //PileUp
  if(eventAOD->IsPileupFromSPD()){
    fCounters.Count(Form("Event:PileUp/Run:%d",runNumber));
    if(fRejectPileUp==1)return;
  }

  //Find Candle
  Int_t trkEntries = (Int_t)event->GetNumberOfTracks();
  
  for(Int_t i=0;i<trkEntries;i++){
    AliAODTrack *track=(AliAODTrack*)event->GetTrack(i);
    UShort_t nClusTPC=track->GetTPCNcls();
    Int_t nSPD=0;
    if(TESTBIT(track->GetITSClusterMap(),1))nSPD++;
    if(TESTBIT(track->GetITSClusterMap(),0))nSPD++;
    if(nClusTPC==0&&(track->GetStatus()&AliESDtrack::kITSrefit)&&!flagITSsa){
      flagITSsa=kTRUE;
      fCounters.Count(Form("Event:CandleITSsa/Run:%d",runNumber));
    }
    if((nClusTPC>=35)&&(track->GetStatus()&AliESDtrack::kITSrefit)&&(track->GetStatus()&AliESDtrack::kTPCrefit)&&(track->Pt()>0.2)&&(!flag35cls)){
      fCounters.Count(Form("Event:Candle35clsTPC/Run:%d",runNumber));
      flag35cls=kTRUE;
    }
    if((nClusTPC>=70)&&(track->GetStatus()&AliESDtrack::kITSrefit)&&(track->GetStatus()&AliESDtrack::kTPCrefit)){
      
      if((track->Pt()>0.2)&&flag0202<2){
	if(!flag02)fCounters.Count(Form("Event:Candles0.2/Run:%d",runNumber));
	flag02=kTRUE;
	flag0202++;
      }
      if((track->Pt()>0.2)&&!flag02spd&&nSPD>=1){
	fCounters.Count(Form("Event:Candles0.2spd1/Run:%d",runNumber));
	flag02spd=kTRUE;
      }
      if((track->Pt()>0.3)&&(!flag03)){
	fCounters.Count(Form("Event:Candles0.3/Run:%d",runNumber));
	flag03=kTRUE;
      }
    }
    if((flag02)&&(flag03)&&flag0202>=2&&flag35cls&&flagITSsa) break; 
  }
  
  if(!(v0A&&v0B)&&(flag02))fCounters.Count(Form("Event:!V0A&Candle02/Run:%d",runNumber));
  if(!(v0A&&v0B)&&(flag03))fCounters.Count(Form("Event:!V0A&Candle03/Run:%d",runNumber));
  if(!(v0A&&v0B)&&flagPVT)fCounters.Count(Form("Event:!V0A&PrimaryVTracks/Run:%d",runNumber));
  if(!(v0A&&v0B)&&flagPV)fCounters.Count(Form("Event:!V0A&PrimaryV/Run:%d",runNumber));
  if(flag0202>1)fCounters.Count(Form("Event:2xCandles0.2/Run:%d",runNumber));
  if(!(v0A&&v0B)&&flag0202>1)fCounters.Count(Form("Event:!V0A&2xCandles02/Run:%d",runNumber));
  
  //delete eventESD;

  return;
}
//_____________________________________________________________________
void AliNormalizationCounter::StoreCandidates(AliVEvent *event,Int_t nCand,Bool_t flagFilter){
  
  Int_t ntracks=event->GetNumberOfTracks();
  if(flagFilter)fHistTrackFilterEvMult->Fill(ntracks,nCand);
  else fHistTrackAnaEvMult->Fill(ntracks,nCand);
  Int_t nSPD=0;
  if(fESD){
    AliESDEvent *ESDevent=(AliESDEvent*)event;
    const AliMultiplicity *alimult = ESDevent->GetMultiplicity();
    nSPD = alimult->GetNumberOfTracklets();

  }else{
    AliAODEvent *aodEvent =(AliAODEvent*)event;
    AliAODTracklets *trklets=aodEvent->GetTracklets();
    nSPD = trklets->GetNumberOfTracklets();
  }
  if(flagFilter)fHistTrackFilterSpdMult->Fill(nSPD,nCand);
  else fHistTrackAnaSpdMult->Fill(nSPD,nCand);
  
  Int_t runNumber = event->GetRunNumber();
  if(nCand==0)return;
  if(flagFilter){
    fCounters.Count(Form("Event:Candid(Filter)/Run:%d",runNumber));
    for(Int_t i=0;i<nCand;i++)fCounters.Count(Form("Event:NCandid(Filter)/Run:%d",runNumber));
  }else{
    fCounters.Count(Form("Event:Candid(Analysis)/Run:%d",runNumber));
    for(Int_t i=0;i<nCand;i++)fCounters.Count(Form("Event:NCandid(Analysis)/Run:%d",runNumber));
  }
  return;
}
//_______________________________________________________________________
TH1D* AliNormalizationCounter::DrawAgainstRuns(TString candle,Bool_t drawHist){
  //
  fCounters.SortRubric("Run");
  TString selection;
  selection.Form("event:%s",candle.Data());
  TH1D* histoneD = fCounters.Get("run",selection.Data());

  histoneD->Sumw2();
  if(drawHist)histoneD->DrawClone();
  return histoneD;
}
//___________________________________________________________________________
TH1D* AliNormalizationCounter::DrawRatio(TString candle1,TString candle2){
  //
  fCounters.SortRubric("Run");
  TString name;

  name.Form("%s/%s",candle1.Data(),candle2.Data());
  TH1D* num=DrawAgainstRuns(candle1.Data(),kFALSE);
  TH1D* den=DrawAgainstRuns(candle2.Data(),kFALSE);

  den->SetTitle(candle2.Data());
  den->SetName(candle2.Data());
  num->Divide(num,den,1,1,"B");
  num->SetTitle(name.Data());
  num->SetName(name.Data());
  num->DrawClone();
  return num;
}
//___________________________________________________________________________
void AliNormalizationCounter::PrintRubrics(){
  fCounters.PrintKeyWords();
}
//___________________________________________________________________________
Double_t AliNormalizationCounter::GetSum(TString candle){
  TString selection="event:";
  selection.Append(candle);
  return fCounters.GetSum(selection.Data());
}
//___________________________________________________________________________
TH2F* AliNormalizationCounter::GetHist(Bool_t filtercuts,Bool_t spdtracklets,Bool_t drawHist){
  if(filtercuts){
    if(spdtracklets){
      if(drawHist)fHistTrackFilterSpdMult->DrawCopy("LEGO2Z 0");
      return fHistTrackFilterSpdMult;
    }else{
      if(drawHist)fHistTrackFilterEvMult->DrawCopy("LEGO2Z 0");
      return fHistTrackFilterEvMult;
    }
  }else{
    if(spdtracklets){
      if(drawHist)fHistTrackAnaSpdMult->DrawCopy("LEGO2Z 0");
      return fHistTrackAnaSpdMult;
    }else{
      if(drawHist)fHistTrackAnaEvMult->DrawCopy("LEGO2Z 0");
      return fHistTrackAnaEvMult;
    }
  }
}
