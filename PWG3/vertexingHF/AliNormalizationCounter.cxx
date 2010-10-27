#include "AliNormalizationCounter.h"
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODEvent.h>
#include <AliVParticle.h>
#include <AliTriggerAnalysis.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TList.h>
#include <TString.h>
#include <TCanvas.h>
#include <AliPhysicsSelection.h>

ClassImp(AliNormalizationCounter)

//____________________________________________
AliNormalizationCounter::AliNormalizationCounter(): 
AliCounterCollection(),
fESD(kFALSE)
{
  // empty constructor
}

//__________________________________________________				
AliNormalizationCounter::AliNormalizationCounter(const char *name): 
AliCounterCollection(name),
fESD(kFALSE)
{
  //default constructor
  AddRubric("Event","triggered/V0AND/Candles0.2/Candles0.25/Candles0.3/2xCandles0.2/PrimaryVTracks/PrimaryV/!V0A&Candle02/!V0A&Candle025/!V0A&Candle03/!V0A&PrimaryVTracks/!V0A&PrimaryV/!V0A&2xCandles02/Candid(Filter)/Candid(Analysis)/NCandid(Filter)/NCandid(Analysis)");//new line
  AddRubric("Run", 10000000);//new line
  Init();//new line
}
//______________________________________________
AliNormalizationCounter::~AliNormalizationCounter()
{
  //destructor
}


//_______________________________________
void AliNormalizationCounter::StoreEvent(AliVEvent *event,Bool_t mc){
  //

  Bool_t v0A=kFALSE; 
  Bool_t v0B=kFALSE;
  Bool_t flag02=kFALSE;
  Bool_t flag025=kFALSE;
  Bool_t flag03=kFALSE;
  Int_t flag0202=0;
  Bool_t flagPV=kFALSE;
  Bool_t flagPVT=kFALSE; 

  //Run Number
  Int_t runNumber = event->GetRunNumber();
 
  //Find CINT1B
  
  AliESDEvent *eventESD = (AliESDEvent*)event;
  if(!eventESD){AliError("ESD event not available");return;}
  if(mc&&event->GetEventType() != 0)return;
  //event must be either physics or MC
  if(!(event->GetEventType() == 7||event->GetEventType() == 0))return;
  
  Count(Form("Event:triggered/Run:%d",runNumber));
      
  //Find V0AND
  AliTriggerAnalysis trAn; /// Trigger Analysis
  v0B = trAn.IsOfflineTriggerFired(eventESD , AliTriggerAnalysis::kV0C);
  v0A = trAn.IsOfflineTriggerFired(eventESD , AliTriggerAnalysis::kV0A);
  if(v0A&&v0B){Count(Form("Event:V0AND/Run:%d",runNumber));}
  
  //Find Candle
  Int_t trkEntries = (Int_t)event->GetNumberOfTracks();
  
  for(Int_t i=0;i<trkEntries;i++){
    AliAODTrack *track=(AliAODTrack*)event->GetTrack(i);
    UShort_t nClusTPC=track->GetTPCNcls();
    if((nClusTPC>=70)&&(track->GetStatus()&AliESDtrack::kITSrefit)&&(track->GetStatus()&AliESDtrack::kTPCrefit)){
      
      if((track->Pt()>0.2)&&flag0202<2){
	if(!flag02)Count(Form("Event:Candles0.2/Run:%d",runNumber));
	flag02=kTRUE;
	flag0202++;
      }
      if((track->Pt()>0.25)&&(!flag025)){
	Count(Form("Event:Candles0.25/Run:%d",runNumber));
	flag025=kTRUE;
      }
      if((track->Pt()>0.3)&&(!flag03)){
	Count(Form("Event:Candles0.3/Run:%d",runNumber));
	flag03=kTRUE;
      }
    }
    if((flag02)&&(flag025)&&(flag03)&&flag0202>=2) break; //i=trkEntries+1;
  }
  
  //FindPrimary vertex  
  AliVVertex *vtrc =  (AliVVertex*)event->GetPrimaryVertex();
  if(vtrc && vtrc->GetNContributors()>0){
    Count(Form("Event:PrimaryV/Run:%d",runNumber));
    flagPV=kTRUE;
  }
  
  if(fESD){
    const AliESDVertex *vtrc1 =  eventESD->GetPrimaryVertexTracks();
    if(vtrc1 && vtrc1->GetNContributors()>0){
      Count(Form("Event:PrimaryVTracks/Run:%d",runNumber));
      flagPVT=kTRUE;
    }
  }
  
  if(!(v0A&&v0B)&&(flag02))Count(Form("Event:!V0A&Candle02/Run:%d",runNumber));
  if(!(v0A&&v0B)&&(flag025))Count(Form("Event:!V0A&Candle025/Run:%d",runNumber));
  if(!(v0A&&v0B)&&(flag03))Count(Form("Event:!V0A&Candle03/Run:%d",runNumber));
  if(!(v0A&&v0B)&&flagPVT)Count(Form("Event:!V0A&PrimaryVTracks/Run:%d",runNumber));
  if(!(v0A&&v0B)&&flagPV)Count(Form("Event:!V0A&PrimaryV/Run:%d",runNumber));
  if(flag0202>1)Count(Form("Event:2xCandles0.2/Run:%d",runNumber));
  if(!(v0A&&v0B)&&flag0202>1)Count(Form("Event:!V0A&2xCandles02/Run:%d",runNumber));
  
  //delete eventESD;

  return;
}
//_____________________________________________________________________
void AliNormalizationCounter::StoreCandidates(AliVEvent *event,Int_t nCand,Bool_t flagFilter){
  //

  Int_t runNumber = event->GetRunNumber();
  if(nCand==0)return;
  if(flagFilter){
    Count(Form("Event:Candid(Filter)/Run:%d",runNumber));
    for(Int_t i=0;i<nCand;i++)Count(Form("Event:NCandid(Filter)/Run:%d",runNumber));
  }else{
    Count(Form("Event:Candid(Analysis)/Run:%d",runNumber));
    for(Int_t i=0;i<nCand;i++)Count(Form("Event:NCandid(Analysis)/Run:%d",runNumber));
  }
  return;
}
//_______________________________________________________________________
TH1D* AliNormalizationCounter::DrawAgainstRuns(TString candle="candid(filter)"){
  //
  TString selection;
  selection.Form("event:%s",candle.Data());
  TH2D* hist = Draw("event","run",selection.Data());
  TH1D* histoneD =(TH1D*)(hist->ProjectionX()->Clone());
  histoneD->DrawClone();
  return histoneD;
}
//___________________________________________________________________________
TH1D* AliNormalizationCounter::DrawRatio(TString candle1="candid(filter)",TString candle2="triggered"){
  //
  TString name;
  name.Form("%s/%s",candle1.Data(),candle2.Data());
  TH1D* num=(TH1D*)(DrawAgainstRuns(candle1.Data())->Clone(name.Data()));
  TH1D* den=DrawAgainstRuns(candle2.Data());
  den->SetTitle(candle2.Data());
  den->SetName(candle2.Data());
  printf("%f %f",num->GetEntries(),den->GetEntries());
  num->Divide(num,den,1,1,"B");
  num->DrawClone();
  return num;
}

