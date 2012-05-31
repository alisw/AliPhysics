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

/* $Id:  $ */

////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//     Implementation of the TPC Raw drift velocity and Altro L1 Phase  calibration   //
//                                                                                    //
//               Origin: Jens Wiechula, J.Wiechula@gsi.de                             //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////
//
//
// *************************************************************************************
// *                                Class Description                                  *
// *************************************************************************************
/*

----example---
TFile f("CalibAltro.root");
AliTPCCalibRaw *al=(AliTPCCalibRaw*)f.Get(f.GetListOfKeys()->At(0)->GetName())
{
TCanvas *c1=(TCanvas*)gROOT->FindObject("c1");
if (!c1) c1=new TCanvas("c1","c1");
c1->Clear();

TH2F h2f("h2","h2;RCU;fail",216,0,216,al->GetNevents(),0,al->GetNevents());
Bool_t first=kTRUE;
Int_t i,iev;
for (i=0;i<216;++i) {
  TVectorF *v=al->GetALTROL1PhaseFailEventsRCU(i);
  if (!v) continue;
  for (iev=0;iev<al->GetNevents();++iev) {
    h2f->SetBinContent(i+1,iev+1,(*v)(iev));
  }
//   TH1F h(*v);
//   h.SetLineColor(i/216.*50+50);
//   ((TH1F*)h.Clone(Form("h%d",i)))->Draw(first?"":"same");
//   c1->Modified();
//   c1->Update();
  first=kFALSE;
}
h2f->Draw("col");
}

*/



//Root includes
#include <TROOT.h>
#include <TH2C.h>
#include <TH1F.h>
#include <TMap.h>
#include <TGraph.h>
#include <TObjString.h>
#include <TTimeStamp.h>
#include <TCanvas.h>
//AliRoot includes
#include "AliTPCCalROC.h"
#include "AliAltroRawStream.h"
#include "AliLog.h"
//class header
#include "AliTPCCalibRaw.h"

ClassImp(AliTPCCalibRaw)

AliTPCCalibRaw::AliTPCCalibRaw() :
  AliTPCCalibRawBase(),
  fPeakDetMinus(1),
  fPeakDetPlus(2),
  fNFailL1Phase(0),
  fNFailL1PhaseEvent(0),
  fNSecTime(600), //default 10 minutes
  fNBinsTime(60), //default 60*10 minutes = 10 hours
  fPadProcessed(kFALSE),
  fCurrentChannel(-1),
  fCurrentSector(-1),
  fLastSector(-2),
  fCurrentRow(-1),
  fCurrentPad(-1),
  fLastTimeBinProc(0),
  fPeakTimeBin(0),
  fLastSignal(0),
  fNOkPlus(0),
  fNOkMinus(0),
  fNanoSec(0),
  fArrCurrentPhaseDist(4),
  fArrCurrentPhase(kNRCU),
  fArrFailEventNumber(100),
  fArrALTROL1Phase(100000),
  fArrALTROL1PhaseEvent(kNRCU),
  fArrALTROL1PhaseFailEvent(kNRCU),
  fHnDrift(0x0),
  fVOccupancyEvent(100000),
  fVSignalSumEvent(100000),
  fVOccupancySenEvent(100000),
  fVSignalSumSenEvent(100000),
  fVNfiredPadsSenEvent(100000),
  fVTimeStampEvent(100000)
{
  //
  // Default ctor
  //
  SetNameTitle("AliTPCCalibRaw","AliTPCCalibRaw");
  CreateDVhist();
  for (Int_t ircu=0;ircu<kNRCU;++ircu) fArrCurrentPhase.GetMatrixArray()[ircu]=-1;
  fFirstTimeBin=850;
  fLastTimeBin=1020;
}
//_____________________________________________________________________
AliTPCCalibRaw::AliTPCCalibRaw(const TMap *config) :
AliTPCCalibRawBase(),
fPeakDetMinus(1),
fPeakDetPlus(2),
fNFailL1Phase(0),
fNFailL1PhaseEvent(0),
fNSecTime(600), //default 10 minutes
fNBinsTime(60), //default 60*10 minutes = 10 hours
fPadProcessed(kFALSE),
fCurrentChannel(-1),
fCurrentSector(-1),
fLastSector(-2),
fCurrentRow(-1),
fCurrentPad(-1),
fLastTimeBinProc(0),
fPeakTimeBin(0),
fLastSignal(0),
fNOkPlus(0),
fNOkMinus(0),
fNanoSec(0),
fArrCurrentPhaseDist(4),
fArrCurrentPhase(kNRCU),
fArrFailEventNumber(100),
fArrALTROL1Phase(100000),
fArrALTROL1PhaseEvent(kNRCU),
fArrALTROL1PhaseFailEvent(kNRCU),
fHnDrift(0x0),
fVOccupancyEvent(100000),
fVSignalSumEvent(100000),
fVOccupancySenEvent(100000),
fVSignalSumSenEvent(100000),
fVNfiredPadsSenEvent(100000),
fVTimeStampEvent(100000)
{
  //
  // Default ctor
  //
  SetNameTitle("AliTPCCalibRaw","AliTPCCalibRaw");
  CreateDVhist();
  for (Int_t ircu=0;ircu<kNRCU;++ircu) fArrCurrentPhase.GetMatrixArray()[ircu]=-1;
  fFirstTimeBin=850;
  fLastTimeBin=1020;
  if (config->GetValue("FirstTimeBin")) fFirstTimeBin = ((TObjString*)config->GetValue("FirstTimeBin"))->GetString().Atoi();
  if (config->GetValue("LastTimeBin")) fLastTimeBin = ((TObjString*)config->GetValue("LastTimeBin"))->GetString().Atoi();
  if (config->GetValue("DebugLevel")) fDebugLevel = ((TObjString*)config->GetValue("DebugLevel"))->GetString().Atoi();
}

//_____________________________________________________________________
AliTPCCalibRaw::~AliTPCCalibRaw()
{
  //
  // dtor
  //
  delete fHnDrift;  
}
//_____________________________________________________________________
// AliTPCCalibRaw& AliTPCCalibRaw::operator = (const  AliTPCCalibRaw &source)
// {
//   //
//   // assignment operator
//   //
//   if (&source == this) return *this;
//   new (this) AliTPCCalibRaw(source);
//   
//   return *this;
// }

//_____________________________________________________________________
Int_t AliTPCCalibRaw::Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
             const Int_t iTimeBin, const Float_t signal)
{
  //
  // Data filling method
  //
  if (iRow<0) return 0;
  if (iPad<0) return 0;
  if (iTimeBin<0) return 0;
  //
  Int_t iChannel  = fROC->GetRowIndexes(isector)[iRow]+iPad; //  global pad position in sector
  //occupancy
  fVOccupancyEvent.GetMatrixArray()[GetNevents()]++;
  fVSignalSumEvent.GetMatrixArray()[GetNevents()]+=signal;
  //occupancy in sensitive regions
  Int_t npads=(Int_t)fROC->GetNPads(isector,iRow);
  Int_t cpad=iPad-npads/2;
  if (isector<(Int_t)fROC->GetNInnerSector()){
    //IROC case (spot)
    if ( iRow>19 && iRow<46 ){
      if ( TMath::Abs(cpad)<7 ){
        fVOccupancySenEvent.GetMatrixArray()[GetNevents()]++;
        fVSignalSumSenEvent.GetMatrixArray()[GetNevents()]+=signal;
        if (iChannel!=fCurrentChannel) fVNfiredPadsSenEvent.GetMatrixArray()[GetNevents()]++;
      }
    }
  } else if ( iRow>75 ){
    //OROC case (outer corners and last three rows are sensitive)
    Int_t padEdge=(Int_t)TMath::Min(iPad,npads-iPad);
    Int_t nrows=(Int_t)fROC->GetNRows(isector);
    if ((nrows-iRow-1)<3 || padEdge<((((Int_t)iRow-76)/4+1))*2){
      fVOccupancySenEvent.GetMatrixArray()[GetNevents()]++;
      fVSignalSumSenEvent.GetMatrixArray()[GetNevents()]+=signal;
      if (iChannel!=fCurrentChannel) fVNfiredPadsSenEvent.GetMatrixArray()[GetNevents()]++;
    }
  }
  //
  if ( (iTimeBin>fLastTimeBin) || (iTimeBin<fFirstTimeBin)   ) return 0;
  //don't process edge pads
  if (IsEdgePad(isector,iRow,iPad)) return 0;
//   Double_t x[kHnBinsDV]={1,isector,0};
//   fHnDrift->Fill(x);
  if (fCurrentChannel==iChannel){
    if (fPadProcessed) return 0;
  } else {
    fPadProcessed=kFALSE;
    fNOkPlus=0;
    fNOkMinus=0;
    fPeakTimeBin=0;
    fLastSignal=0;
  }
//   Double_t x2[kHnBinsDV]={2,isector,0};
//   fHnDrift->Fill(x2);
  

  if (signal>fLastSignal) ++fNOkPlus;
  else if(signal<fLastSignal && fNOkPlus>=fPeakDetPlus){
    ++fNOkMinus;
    if (!fPeakTimeBin) fPeakTimeBin=fLastTimeBinProc;
    if ( fNOkMinus>=fPeakDetMinus ) {
      Double_t x[kHnBinsDV]={fPeakTimeBin,isector,(fTimeStamp-fFirstTimeStamp)/fNSecTime};
      fHnDrift->Fill(x);
    } 
  } else {
    fNOkPlus=0;
    fNOkMinus=0;
    fPeakTimeBin=0;
    fLastSignal=0;
  }

  fLastTimeBinProc=iTimeBin;
  fLastSignal=TMath::Nint(signal);
  fCurrentChannel = iChannel;
  return 0;
}
//_____________________________________________________________________
void AliTPCCalibRaw::UpdateDDL(){
  //
  // fill ALTRO L1 information
  //
  
  //set nanoseconds
  if (!fNanoSec) {
    TTimeStamp s;
    fNanoSec=s.GetNanoSec();
  }
  // current phase
  Int_t phase=(Int_t)(GetL1PhaseTB()*4.);
  //Fill pahse information of current rcu and event
  fArrCurrentPhase.GetMatrixArray()[fCurrDDLNum]=phase;
  //increase phase counter
  ++((fArrCurrentPhaseDist.GetMatrixArray())[phase]);
  
}
//_____________________________________________________________________
void AliTPCCalibRaw::ResetEvent()
{
  //
  // Reset event counters
  //

  fCurrentChannel=-1;
  fCurrentRow=-1;
  fCurrentPad=-1;
  fArrCurrentPhaseDist.Zero();
}
//_____________________________________________________________________
void AliTPCCalibRaw::EndEvent()
{
  //
  // End event analysis
  //

  
  //find phase of the current event
  Int_t phaseMaxEntries=-1;
  Int_t maxEntries=0;
  for (Int_t i=0;i<fArrCurrentPhaseDist.GetNrows();++i){
    Int_t entries=(Int_t)fArrCurrentPhaseDist[i];
    if (maxEntries<entries) {
      maxEntries=entries;
      phaseMaxEntries=i;
    }
  }
  // store phase of current event
  if (fArrALTROL1Phase.GetNrows()-1<=GetNevents())
    fArrALTROL1Phase.ResizeTo(GetNevents()+10000);
  (fArrALTROL1Phase.GetMatrixArray())[GetNevents()]=phaseMaxEntries;
  
  //loop over RCUs and test failures
  UInt_t fail=0;
  for (Int_t ircu=0;ircu<kNRCU;++ircu){
    Int_t phase=(Int_t)fArrCurrentPhase[ircu];
    if (phase<0) continue;
    if (phase!=phaseMaxEntries){
      TVectorF *arr=MakeArrL1PhaseRCU(fCurrDDLNum,kTRUE);
      if (arr->GetNrows()-1<=(Int_t)fNFailL1PhaseEvent) arr->ResizeTo(arr->GetNrows()+100);
      (arr->GetMatrixArray())[fNFailL1PhaseEvent]=phase;
      ++fNFailL1Phase;
      fail=1;
      }
    //reset current phase information
    fArrCurrentPhase[ircu]=-1;
  }
  if (fail){
    if (fArrFailEventNumber.GetNrows()-1<=(Int_t)fNFailL1PhaseEvent) fArrFailEventNumber.ResizeTo(fArrFailEventNumber.GetNrows()+100);
    fArrFailEventNumber.GetMatrixArray()[fNFailL1PhaseEvent]=GetNevents();
  }
  fNFailL1PhaseEvent+=fail;
  //time stamps
  fVTimeStampEvent.GetMatrixArray()[GetNevents()]=GetTimeStamp()-fFirstTimeStamp+fNanoSec*1e-9;
  fNanoSec=0;
  //occupance related
  if (fVOccupancyEvent.GetNrows()-1<=GetNevents()){
    fVOccupancyEvent.ResizeTo(GetNevents()+10000);
    fVSignalSumEvent.ResizeTo(GetNevents()+10000);
    fVOccupancySenEvent.ResizeTo(GetNevents()+10000);
    fVSignalSumSenEvent.ResizeTo(GetNevents()+10000);
    fVTimeStampEvent.ResizeTo(GetNevents()+10000);
    fVNfiredPadsSenEvent.ResizeTo(GetNevents()+10000);
  }
  IncrementNevents();
}
//_____________________________________________________________________
TH2C *AliTPCCalibRaw::MakeHistL1RCUEvents(Int_t type)
{
  // Create a 2D histo RCU:Events indicating the there was a deviation
  // from the mean L1 phase of the event
  //
  //type: 0=Failures, 1=Phases

  //number of relavant events, depending on version
  Int_t nevents=GetNevents();
  //check version
  Bool_t newVersion=kFALSE;
  for (Int_t ircu=0; ircu<kNRCU; ++ircu){
    const TVectorF *v=GetALTROL1PhaseEventsRCU(ircu);
    if (!v) continue;
    if ((UInt_t)(v->GetNrows())==fNFailL1PhaseEvent){
      newVersion=kTRUE;
      nevents=fNFailL1PhaseEvent;
    }
    break;
  }
  TH2C *h2 = new TH2C("hL1FailRCUEvents","L1 Failures;RCU;Event",kNRCU,0,kNRCU,nevents,0,nevents);
  Int_t add=0;
  for (Int_t ircu=0;ircu<kNRCU;++ircu) {
    const TVectorF *v=GetALTROL1PhaseEventsRCU(ircu);
    if (type==0){
      add=1;
      h2->SetMinimum(0);
      h2->SetMaximum(2);
    } else if (type==1) {
      add=0;
      h2->SetMinimum(0);
      h2->SetMaximum(4);
    }
    if (!v) continue;
    for (Int_t iev=0;iev<nevents;++iev) {
      Float_t val=(*v)(iev);
      Float_t phase=fArrALTROL1Phase.GetMatrixArray()[iev];
      if (newVersion) {
        Int_t event=(Int_t)fArrFailEventNumber.GetMatrixArray()[iev];
        phase=fArrALTROL1Phase.GetMatrixArray()[event];
      }
      if (type==0) val=(val!=phase);
      h2->SetBinContent(ircu+1,iev+1,val+add);
    }
  }
  return h2;
}
//_____________________________________________________________________
TH1F *AliTPCCalibRaw::MakeHistL1PhaseDist()
{
  //
  // L1 phase distribution. Should be flat in ideal case
  //
  TH1F *h=new TH1F("L1phaseDist","Normalized L1 phase distribution;phase;fraction of events",4,0,4);
  h->Sumw2();
  for (Int_t iev=0;iev<GetNevents();++iev) h->Fill(fArrALTROL1Phase.GetMatrixArray()[iev]);
  if (GetNevents()>0) h->Scale(1./GetNevents());
  h->SetMinimum(0);
  h->SetMaximum(1);
  return h;
}
//_____________________________________________________________________
TVectorF *AliTPCCalibRaw::MakeVectL1PhaseDist()
{
  //
  // L1 phase distribution. Should be flat in ideal case
  //
  TVectorF *v=new TVectorF(4);
  for (Int_t iev=0;iev<GetNevents();++iev) {
    Int_t phase=(Int_t)fArrALTROL1Phase.GetMatrixArray()[iev];
    ((v->GetMatrixArray())[phase])+=1./GetNevents();
  }
  return v;
}
//_____________________________________________________________________
TH2C *AliTPCCalibRaw::MakeHistL1RCUEventsIROC(Int_t type)
{
  //
  // Create a 2D histo RCU:Events indicating the there was a deviation
  // from the mean L1 phase of the event
  //
  TH2C *h2 = new TH2C("hL1FailRCUEventsIROC","L1 Failures IROCs;RCU;Event",72,0,36,GetNevents(),0,GetNevents());
  for (Int_t ircu=0;ircu<72;++ircu) {
    const TVectorF *v=0;
    if (type==0)      v=GetALTROL1PhaseFailEventsRCU(ircu);
    else if (type==1) v=GetALTROL1PhaseEventsRCU(ircu);
    if (!v) continue;
    for (Int_t iev=0;iev<GetNevents();++iev) {
      h2->SetBinContent(ircu+1,iev+1,(*v)(iev));
    }
  }
  return h2;
}
//_____________________________________________________________________
TH2C *AliTPCCalibRaw::MakeHistL1RCUEventsOROC(Int_t type)
{
  //
  // Create a 2D histo RCU:Events indicating the there was a deviation
  // from the mean L1 phase of the event
  //
  TH2C *h2 = new TH2C("hL1FailRCUEventsOROC","L1 Failures OROCs;RCU;Event",144,0,36,GetNevents(),0,GetNevents());
  for (Int_t ircu=72;ircu<kNRCU;++ircu) {
    const TVectorF *v=0;
    if (type==0)      v=GetALTROL1PhaseFailEventsRCU(ircu);
    else if (type==1) v=GetALTROL1PhaseEventsRCU(ircu);
    if (!v) continue;
    for (Int_t iev=0;iev<GetNevents();++iev) {
      h2->SetBinContent(ircu-72+1,iev+1,(*v)(iev));
    }
  }
  return h2;
}
//_____________________________________________________________________
void AliTPCCalibRaw::CreateDVhist()
{
  //
  // Setup the HnSparse for the drift velocity determination
  //
  if (fHnDrift) return;
  //HnSparse bins
  //time bin, roc, time
  Int_t    bins[kHnBinsDV] = {fLastTimeBin-fFirstTimeBin, 72, fNBinsTime};
  Double_t xmin[kHnBinsDV] = {fFirstTimeBin,0,0};
  Double_t xmax[kHnBinsDV] = {fLastTimeBin,72,fNBinsTime};
  fHnDrift=new THnSparseI("fHnDrift",Form("Drift velocity using last time bin;time bin[#times 100ns];ROC;Time bin [#times %us]",fNSecTime),kHnBinsDV, bins, xmin, xmax);
    
}
//_____________________________________________________________________
void AliTPCCalibRaw::Analyse()
{
  //
  // Analyse Data
  //

  //resize arrays
  fArrALTROL1Phase.ResizeTo(GetNevents());
  for (Int_t ircu=0;ircu<kNRCU;++ircu){
    TVectorF *arr=MakeArrL1PhaseRCU(ircu);//MakeArrL1PhaseRCU(ircu);
    if (!arr) continue;
    arr->ResizeTo(fNFailL1PhaseEvent);
    fArrFailEventNumber.ResizeTo(fNFailL1PhaseEvent);
//    TVectorF *arrF=MakeArrL1PhaseFailRCU(ircu);
//     arrF->ResizeTo(1);
  }
  //resize occupancy arrays only save event occupancy in sensitive regions by default
  //save the rest in debub mode
  fVOccupancySenEvent.ResizeTo(GetNevents());
  if (fDebugLevel>0){
    fVOccupancyEvent.ResizeTo(GetNevents());
    fVSignalSumEvent.ResizeTo(GetNevents());
    fVSignalSumSenEvent.ResizeTo(GetNevents());
    fVNfiredPadsSenEvent.ResizeTo(GetNevents());
    fVTimeStampEvent.ResizeTo(GetNevents());
  } else {
    fVOccupancyEvent.ResizeTo(0);
    fVSignalSumEvent.ResizeTo(0);
    fVSignalSumSenEvent.ResizeTo(0);
    fVNfiredPadsSenEvent.ResizeTo(0);
    fVTimeStampEvent.ResizeTo(0);
  }
  //Analyse drift velocity TODO
  
}
//_____________________________________________________________________
TGraph* AliTPCCalibRaw::MakeGraphOccupancy(const Int_t type, const Int_t xType)
{
  //
  // create occupancy graph (of samples abouve threshold)
  // type=0:   number of samples
  // type=1:   mean data volume (ADC counts/sample)
  // type=2:   data volume (ADC counts)
  // type=3:   samples per ADC count
  // type=4:   sample occupancy
  //
  // type=5: number of sample sensitive / number of samples
  //
  // same in sensitive regions:
  // type=10:  number of samples
  // type=11:  mean data volume (ADC counts/sample)
  // type=12:  data volume (ADC counts)
  // type=13:  samples per ADC count
  // type=14:  sample occupancy
  //
  // type=16: number of samples sensitive / number of pads sensitive
  // type=17: pad occupancy in sensitive regions
  // xType=0:  vs. time stamp
  // xType=1:  vs. event counter
  //

  TString title("Event occupancy");
  TString xTitle("Time");
  TString yTitle("number of samples");
  TGraph *gr=new TGraph(GetNevents());
  if (fVSignalSumEvent.GetNrows()==0&&!(type==10||type==14)) return 0;
  TVectorF *vOcc=&fVOccupancyEvent;
  TVectorF *vSum=&fVSignalSumEvent;
  TVectorF *vPads=&fVNfiredPadsSenEvent;
  Double_t norm=557568.;
  if (type!=14&&fVOccupancyEvent.GetNrows()==0){
    AliWarning("In non debug mode only occupancy in sensitive regions vs. event awailable!!!");
    return 0;
  }
  if (type>=10){
    vOcc=&fVOccupancySenEvent;
    vSum=&fVSignalSumSenEvent;
    vPads=&fVNfiredPadsSenEvent;
    norm=33012.;
  }
  for (Int_t i=0;i<GetNevents(); ++i){
    Double_t nAboveThreshold=vOcc->GetMatrixArray()[i];
    
    Double_t nSumADC        =1;
    Double_t timestamp      =1;
    Double_t nPads          =1;

    if (fVOccupancyEvent.GetNrows()>0){
      nSumADC        =vSum->GetMatrixArray()[i];
      timestamp      =fVTimeStampEvent.GetMatrixArray()[i]+fFirstTimeStamp;
      nPads          =vPads->GetMatrixArray()[i];
    }
    Double_t x=timestamp;
    Double_t y=0;
    //
    if (xType==1)     x=i;
    //
    if (type%10==0)            y=nAboveThreshold;
    if (type%10==1&&nAboveThreshold>0)      y=nSumADC/nAboveThreshold;
    if (type%10==2)            y=nSumADC;
    if (type%10==3&&nSumADC>0) y=nAboveThreshold/nSumADC;
    if (type%10==4)            y=nAboveThreshold/(norm*(fLastTimeBin-fFirstTimeBin));
    if (type==5)               y=fVOccupancySenEvent.GetMatrixArray()[i]/fVOccupancyEvent.GetMatrixArray()[i];
    if (type==16&&nPads>0)     y=nAboveThreshold/nPads;
    if (type==17)              y=nPads/norm;
    //
    gr->SetPoint(i,x,y);
  }
  if (xType==1) xTitle="Event";
  if (type%10==1)      yTitle="Mean ADC counts/sample";
  else if (type%10==2) yTitle="Data volume [ADC counts]";
  else if (type%10==3) yTitle="samples per ADC count";
  else if (type%10==4) yTitle="sample occupancy";
  if (type==5)         yTitle="N samples (sensitive) / N samples";
  if (type%10==6)      yTitle="N samples / N pads";
  if (type==17)        yTitle="Pad Occupancy";
  if (type>=10)        yTitle+=" (sensitive)";
  title=yTitle+":"+xTitle;
  title+=";"+xTitle+";"+yTitle;
  gr->SetTitle(title.Data());
  gr->SetEditable(kFALSE);
  return gr;
}
//_____________________________________________________________________
// TGraph* AliTPCCalibRaw::MakeGraphNoiseEvents()
// {
  //
  // Not implemented for the moment
  //
//   return 0;  
// }
//_____________________________________________________________________
TCanvas* AliTPCCalibRaw::MakeCanvasOccupancy(const Int_t xType, Bool_t sen)
{
  //
  // Create a canvas with occupancy information of all 'type's (see MakeGraphOccupancy)
  // xType=0: vs. timestamp
  // xType=1: vs. event number
  //
  // sen=kTRUE: for sensitive regions
  //

  TString name("RawOccupancy_");
  TString title("Raw Occupancy vs. ");
  if (xType==0){
    name+="Time";
    title+="time";
  } else if (xType==1){
    name+="Event";
    title+="event";
  }
  if (sen){
    name+="Sen";
    title+=" (sensitive)";
  }
  TCanvas *c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(name.Data());
  if (!c) c=new TCanvas(name.Data(),title.Data());
  c->Clear();
  c->Divide(2,2);
  for (Int_t i=0;i<4;++i){
    c->cd(i+1);
    TGraph *gr=MakeGraphOccupancy(i+10*(Int_t)sen,xType);
    gr->Draw("alp");
  }
  return c;
}

//_____________________________________________________________________
void AliTPCCalibRaw::Merge(AliTPCCalibRaw * const sig)
{
  //
  // Merge sig with this instance
  //

  if (!sig) return;
  MergeBase(sig);
  //Add last time bin distribution histogram
  fHnDrift->Add(sig->fHnDrift);

  //Add occupancy data
  
}

//_____________________________________________________________________
Long64_t AliTPCCalibRaw::Merge(TCollection * const list)
{
  //
  // Merge all objects of this type in list
  //
  
  Long64_t nmerged=1;
  
  TIter next(list);
  AliTPCCalibRaw *ce=0;
  TObject *o=0;
  
  while ( (o=next()) ){
    ce=dynamic_cast<AliTPCCalibRaw*>(o);
    if (ce){
      Merge(ce);
      ++nmerged;
    }
  }
  
  return nmerged;
}

