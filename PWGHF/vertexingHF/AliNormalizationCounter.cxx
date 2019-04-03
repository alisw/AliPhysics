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

#include "AliLog.h"
#include "AliNormalizationCounter.h"
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliAODEvent.h>
#include <AliAODMCHeader.h>
#include <AliAODVZERO.h>
#include <AliVParticle.h>
#include <AliTriggerAnalysis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>
#include <TCanvas.h>
#include <AliPhysicsSelection.h>
#include <AliMultiplicity.h>

/// \cond CLASSIMP
ClassImp(AliNormalizationCounter);
/// \endcond

//____________________________________________
AliNormalizationCounter::AliNormalizationCounter(): 
TNamed(),
fCounters(),
fESD(kFALSE),
fMultiplicity(kFALSE),
fMultiplicityEtaRange(1.0),
fSpherocity(kFALSE),
fSpherocitySteps(100.),
fHistTrackFilterEvMult(0),
fHistTrackAnaEvMult(0),
fHistTrackFilterSpdMult(0),
fHistTrackAnaSpdMult(0),
fHistGenVertexZ(0),
fHistGenVertexZRecoPV(0),
fHistRecoVertexZ(0)
{
  // empty constructor
}

//__________________________________________________				
AliNormalizationCounter::AliNormalizationCounter(const char *name): 
TNamed(name,name),
fCounters(name),
fESD(kFALSE),
fMultiplicity(kFALSE),
fMultiplicityEtaRange(1.0),
fSpherocity(kFALSE),
fSpherocitySteps(100.),
fHistTrackFilterEvMult(0),
fHistTrackAnaEvMult(0),
fHistTrackFilterSpdMult(0),
fHistTrackAnaSpdMult(0),
fHistGenVertexZ(0),
fHistGenVertexZRecoPV(0),
fHistRecoVertexZ(0)
{
  ;
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
  delete fHistGenVertexZ;
  delete fHistGenVertexZRecoPV;
  delete fHistRecoVertexZ;
}

//______________________________________________
void AliNormalizationCounter::Init()
{
  //variables initialization
  fCounters.AddRubric("Event","triggered/V0AND/PileUp/PbPbC0SMH-B-NOPF-ALLNOTRD/Candles0.3/PrimaryV/countForNorm/noPrimaryV/zvtxGT10/!V0A&Candle03/!V0A&PrimaryV/Candid(Filter)/Candid(Analysis)/NCandid(Filter)/NCandid(Analysis)");
  if(fMultiplicity)  fCounters.AddRubric("Multiplicity", 5000);
  if(fSpherocity)  fCounters.AddRubric("Spherocity", (Int_t)fSpherocitySteps+1);
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
  fHistGenVertexZ=new TH1F("hGenVertexZ","generated z vertex ; z_{vertex}^{MC} (cm) ; counts",600,-30,30);
  fHistGenVertexZRecoPV=new TH1F("hGenVertexZRecoPV","generated z vertex (events with reconstructed vertex); z_{vertex}^{MC} (cm) ; counts",600,-30,30);
  fHistRecoVertexZ=new TH1F("hRecoVertexZ","reconstructed z vertex; z_{vertex}^{rec} (cm) ; counts",600,-30,30);
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
  fHistGenVertexZ->Add(norm->fHistGenVertexZ);
  fHistGenVertexZRecoPV->Add(norm->fHistGenVertexZRecoPV);
  fHistRecoVertexZ->Add(norm->fHistRecoVertexZ);
}
//_______________________________________
/*
Stores the variables used for normalization as function of run number
returns kTRUE if the event is to be counted for normalization
(pass event selection cuts OR has no primary vertex)
 */
void AliNormalizationCounter::StoreEvent(AliVEvent *event,AliRDHFCuts *rdCut,Bool_t mc, Int_t multiplicity, Double_t spherocity){
  //

  Bool_t isEventSelected = rdCut->IsEventSelected(event);

  // events not passing physics selection. do nothing
  if(rdCut->IsEventRejectedDuePhysicsSelection()) return;

  Bool_t v0A=kFALSE; 
  Bool_t v0B=kFALSE;
  Bool_t flag03=kFALSE;
  Bool_t flagPV=kFALSE;

  //Run Number
  Int_t runNumber = event->GetRunNumber();
 
  // Evaluate the multiplicity
  if(fMultiplicity && multiplicity==-9999) Multiplicity(event);

  //Find CINT1B
  AliESDEvent *eventESD = (AliESDEvent*)event;
  if(!eventESD){AliError("ESD event not available");return;}
  if(mc&&event->GetEventType() != 0)return;
  //event must be either physics or MC
  if(!(event->GetEventType() == 7||event->GetEventType() == 0))return;
  
  FillCounters("triggered",runNumber,multiplicity,spherocity);

  //Find V0AND
  AliTriggerAnalysis trAn; /// Trigger Analysis
  AliAODVZERO* aodV0 = (AliAODVZERO*)event->GetVZEROData();
  Bool_t isPP2012 = kFALSE;
  if(runNumber>=176326 && runNumber<=193766) isPP2012=kTRUE;
  if(aodV0 && !isPP2012){
    v0B = trAn.IsOfflineTriggerFired(eventESD , AliTriggerAnalysis::kV0C);
    v0A = trAn.IsOfflineTriggerFired(eventESD , AliTriggerAnalysis::kV0A);
  }
  if(v0A&&v0B) FillCounters("V0AND",runNumber,multiplicity,spherocity);
  
  //FindPrimary vertex  
  // AliVVertex *vtrc =  (AliVVertex*)event->GetPrimaryVertex();
  // if(vtrc && vtrc->GetNContributors()>0){
  //   fCounters.Count(Form("Event:PrimaryV/Run:%d",runNumber));
  //   flagPV=kTRUE;
  // }

  //trigger
  AliAODEvent *eventAOD = (AliAODEvent*)event;
  TString trigclass=eventAOD->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL")){
    FillCounters("PbPbC0SMH-B-NOPF-ALLNOTRD",runNumber,multiplicity,spherocity);
  }

  //FindPrimary vertex  
  if(isEventSelected){
    FillCounters("PrimaryV",runNumber,multiplicity,spherocity);
    flagPV=kTRUE;
  }else{
    if(rdCut->GetWhyRejection()==0){
      FillCounters("noPrimaryV",runNumber,multiplicity,spherocity);
    }
    //find good vtx outside range
    if(rdCut->GetWhyRejection()==6){
      FillCounters("zvtxGT10",runNumber,multiplicity,spherocity);
      FillCounters("PrimaryV",runNumber,multiplicity,spherocity);
      flagPV=kTRUE;
    }
    if(rdCut->GetWhyRejection()==1){
      FillCounters("PileUp",runNumber,multiplicity,spherocity);
    }
  }
  //to be counted for normalization
  if(rdCut->CountEventForNormalization()){
    FillCounters("countForNorm",runNumber,multiplicity,spherocity);
  }
  // fill histograms of vertex position
  if(mc){
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(eventAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if(mcHeader){
      Double_t zMCVertex=mcHeader->GetVtxZ();
      fHistGenVertexZ->Fill(zMCVertex);
      if(flagPV) fHistGenVertexZRecoPV->Fill(zMCVertex);
    }
  }
  if(flagPV){
    AliVVertex *vtrc =  (AliVVertex*)event->GetPrimaryVertex();
    if(vtrc && vtrc->GetNContributors()>0) fHistRecoVertexZ->Fill(vtrc->GetZ()); 
  }

  //Find Candle
  Int_t trkEntries = (Int_t)event->GetNumberOfTracks();
  for(Int_t i=0;i<trkEntries&&!flag03;i++){
    AliAODTrack *track=(AliAODTrack*)event->GetTrack(i);
    if((track->Pt()>0.3)&&(!flag03)){
      FillCounters("Candles0.3",runNumber,multiplicity,spherocity);
      flag03=kTRUE;
      break;
    }
  }
  
  if(!(v0A&&v0B)&&(flag03)){ 
    FillCounters("!V0A&Candle03",runNumber,multiplicity,spherocity);
  }
  if(!(v0A&&v0B)&&flagPV){
    FillCounters("!V0A&PrimaryV",runNumber,multiplicity,spherocity);
  }
  
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
  Int_t multiplicity = Multiplicity(event);
  if(nCand==0)return;
  if(flagFilter){
    if(fMultiplicity) 
      fCounters.Count(Form("Event:Candid(Filter)/Run:%d/Multiplicity:%d",runNumber,multiplicity));
    else 
      fCounters.Count(Form("Event:Candid(Filter)/Run:%d",runNumber));
    for(Int_t i=0;i<nCand;i++){ 
      if(fMultiplicity) 
	fCounters.Count(Form("Event:NCandid(Filter)/Run:%d/Multiplicity:%d",runNumber,multiplicity));
      else 
	fCounters.Count(Form("Event:NCandid(Filter)/Run:%d",runNumber));
    }
  }else{
    if(fMultiplicity) 
      fCounters.Count(Form("Event:Candid(Analysis)/Run:%d/Multiplicity:%d",runNumber,multiplicity));
    else
      fCounters.Count(Form("Event:Candid(Analysis)/Run:%d",runNumber));
    for(Int_t i=0;i<nCand;i++){ 
      if(fMultiplicity) 
	fCounters.Count(Form("Event:NCandid(Analysis)/Run:%d/Multiplicity:%d",runNumber,multiplicity));
      else
	fCounters.Count(Form("Event:NCandid(Analysis)/Run:%d",runNumber));
    }
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
//___________________________________________________________________________
Double_t AliNormalizationCounter::GetNEventsForNorm(){
  Double_t noVtxzGT10=GetSum("noPrimaryV")*GetSum("zvtxGT10")/GetSum("PrimaryV");
  return GetSum("countForNorm")-noVtxzGT10;
}
//___________________________________________________________________________
Double_t AliNormalizationCounter::GetNEventsForNorm(Int_t runnumber){
  TString listofruns = fCounters.GetKeyWords("RUN");
  if(!listofruns.Contains(Form("%d",runnumber))){
    printf("WARNING: %d is not a valid run number\n",runnumber);
    fCounters.Print("Run","",kTRUE);
    return 0.;
  }
  TString suffix;suffix.Form("/RUN:%d",runnumber);
  TString zvtx;zvtx.Form("zvtxGT10%s",suffix.Data());
  TString noPV;noPV.Form("noPrimaryV%s",suffix.Data());
  TString pV;pV.Form("PrimaryV%s",suffix.Data());
  TString tbc;tbc.Form("countForNorm%s",suffix.Data());
  Double_t noVtxzGT10=GetSum(noPV.Data())*GetSum(zvtx.Data())/GetSum(pV.Data());
  return GetSum(tbc.Data())-noVtxzGT10;
}

//___________________________________________________________________________
Double_t AliNormalizationCounter::GetNEventsForNorm(Int_t minmultiplicity, Int_t maxmultiplicity){

  if(!fMultiplicity) {
    AliInfo("Sorry, you didn't activate the multiplicity in the counter!");
    return 0.;
  }

  TString listofruns = fCounters.GetKeyWords("Multiplicity");

  Int_t nmultbins = maxmultiplicity - minmultiplicity;
  Double_t sumnoPV=0., sumZvtx=0., sumPv=0., sumEvtNorm=0.;
  for (Int_t ibin=0; ibin<=nmultbins; ibin++) {
    //    cout << " Looking at bin "<< ibin+minmultiplicity<<endl;
    if(!listofruns.Contains(Form("%d",ibin+minmultiplicity))){
      //      AliInfo(Form("WARNING: %d is not a valid multiplicity number. \n",ibin+minmultiplicity));
      continue;
    }
    TString suffix;suffix.Form("/Multiplicity:%d",ibin+minmultiplicity);
    TString zvtx;zvtx.Form("zvtxGT10%s",suffix.Data());
    TString noPV;noPV.Form("noPrimaryV%s",suffix.Data());
    TString pV;pV.Form("PrimaryV%s",suffix.Data());
    TString tbc;tbc.Form("countForNorm%s",suffix.Data());
    sumnoPV += GetSum(noPV.Data());
    sumZvtx += GetSum(zvtx.Data());
    sumPv += GetSum(pV.Data());
    sumEvtNorm += GetSum(tbc.Data());
  }
  Double_t noVtxzGT10 = sumPv>0. ? sumnoPV * sumZvtx / sumPv : 0.;
  return sumEvtNorm - noVtxzGT10;
}
//___________________________________________________________________________
Double_t AliNormalizationCounter::GetNEventsForNorm(Int_t minmultiplicity, Int_t maxmultiplicity, Double_t minspherocity, Double_t maxspherocity){

  if(!fMultiplicity || !fSpherocity) {
    AliInfo("You must activate both multiplicity and spherocity in the counters to use this method!");
    return 0.;
  }

  TString listofruns = fCounters.GetKeyWords("Multiplicity");
  TString listofruns2 = fCounters.GetKeyWords("Spherocity");
  TObjArray* arr=listofruns2.Tokenize(",");
  Int_t nSphVals=arr->GetEntries();

  Int_t nmultbins = maxmultiplicity - minmultiplicity;
  Int_t minSphToInteger=minspherocity*fSpherocitySteps;
  Int_t maxSphToInteger=maxspherocity*fSpherocitySteps;
  Int_t nstbins = maxSphToInteger - minSphToInteger;
  Double_t sumnoPV=0., sumZvtx=0., sumPv=0., sumEvtNorm=0.;
  for (Int_t ibin=0; ibin<=nmultbins; ibin++) {
    if(!listofruns.Contains(Form("%d",ibin+minmultiplicity))) continue;
    for (Int_t ibins=0; ibins<nstbins; ibins++) {
      Bool_t analyze=kFALSE;
      for(Int_t j=0; j<nSphVals; j++) if((((TObjString*)arr->At(j))->String()).Atoi()==(ibins+minSphToInteger)) analyze=kTRUE;
      if(!analyze) continue;
      if(listofruns2.Contains(",") && !listofruns2.Contains(Form("%d",ibins+minSphToInteger))) continue;
      TString suffix;suffix.Form("/Multiplicity:%d/Spherocity:%d",ibin+minmultiplicity,ibins+minSphToInteger);
      TString zvtx;zvtx.Form("zvtxGT10%s",suffix.Data());
      TString noPV;noPV.Form("noPrimaryV%s",suffix.Data());
      TString pV;pV.Form("PrimaryV%s",suffix.Data());
      TString tbc;tbc.Form("countForNorm%s",suffix.Data());
      sumnoPV += GetSum(noPV.Data());
      sumZvtx += GetSum(zvtx.Data());
      sumPv += GetSum(pV.Data());
      sumEvtNorm += GetSum(tbc.Data());
    }
  }
  delete arr;
  Double_t noVtxzGT10 = sumPv>0. ? sumnoPV * sumZvtx / sumPv : 0.;
  return sumEvtNorm - noVtxzGT10;
}

//___________________________________________________________________________
Double_t AliNormalizationCounter::GetNEventsForNormSpheroOnly(Double_t minspherocity, Double_t maxspherocity){

  if(!fSpherocity) {
    AliInfo("Sorry, you didn't activate the sphericity in the counter!");
    return 0.;
  }

  TString listofruns = fCounters.GetKeyWords("Spherocity");
  TObjArray* arr=listofruns.Tokenize(",");
  Int_t nSphVals=arr->GetEntries();

  Int_t minSphToInteger=minspherocity*fSpherocitySteps;
  Int_t maxSphToInteger=maxspherocity*fSpherocitySteps;
  Int_t nstbins = maxSphToInteger - minSphToInteger;
  Double_t sumnoPV=0., sumZvtx=0., sumPv=0., sumEvtNorm=0.;
  for (Int_t ibin=0; ibin<nstbins; ibin++) {
    Bool_t analyze=kFALSE;
    for(Int_t j=0; j<nSphVals; j++) if((((TObjString*)arr->At(j))->String()).Atoi()==(ibin+minSphToInteger)) analyze=kTRUE;
    if(!analyze) continue;
   
    TString suffix;suffix.Form("/Spherocity:%d",ibin+minSphToInteger);
    TString zvtx;zvtx.Form("zvtxGT10%s",suffix.Data());
    TString noPV;noPV.Form("noPrimaryV%s",suffix.Data());
    TString pV;pV.Form("PrimaryV%s",suffix.Data());
    TString tbc;tbc.Form("countForNorm%s",suffix.Data());
    sumnoPV += GetSum(noPV.Data());
    sumZvtx += GetSum(zvtx.Data());
    sumPv += GetSum(pV.Data());
    sumEvtNorm += GetSum(tbc.Data());
  }
  delete arr;
  Double_t noVtxzGT10 = sumPv>0. ? sumnoPV * sumZvtx / sumPv : 0.;
  return sumEvtNorm - noVtxzGT10;
}
//___________________________________________________________________________
Double_t AliNormalizationCounter::GetSum(TString candle,Int_t minmultiplicity, Int_t maxmultiplicity){
  // counts events of given type in a given multiplicity range

  if(!fMultiplicity) {
    AliInfo("Sorry, you didn't activate the multiplicity in the counter!");
    return 0.;
  }

  TString listofruns = fCounters.GetKeyWords("Multiplicity");
  Double_t sum=0.;
  for (Int_t ibin=minmultiplicity; ibin<=maxmultiplicity; ibin++) {
    //    cout << " Looking at bin "<< ibin+minmultiplicity<<endl;
    if(!listofruns.Contains(Form("%d",ibin))){
      //      AliInfo(Form("WARNING: %d is not a valid multiplicity number. \n",ibin));
      continue;
    }
    TString suffix=Form("/Multiplicity:%d",ibin);
    TString name=Form("%s%s",candle.Data(),suffix.Data());
    sum += GetSum(name.Data());
  }
  return sum;
}

//___________________________________________________________________________
TH1D* AliNormalizationCounter::DrawNEventsForNorm(Bool_t drawRatio){
  //usare algebra histos
  fCounters.SortRubric("Run");
  TString selection;

  selection.Form("event:noPrimaryV");
  TH1D* hnoPrimV = fCounters.Get("run",selection.Data());
  hnoPrimV->Sumw2();

  selection.Form("event:zvtxGT10");
  TH1D*  hzvtx= fCounters.Get("run",selection.Data());
  hzvtx->Sumw2();

  selection.Form("event:PrimaryV");
  TH1D* hPrimV = fCounters.Get("run",selection.Data());
  hPrimV->Sumw2();

  hzvtx->Multiply(hnoPrimV);
  hzvtx->Divide(hPrimV);

  selection.Form("event:countForNorm");
  TH1D* hCountForNorm = fCounters.Get("run",selection.Data());
  hCountForNorm->Sumw2();

  hCountForNorm->Add(hzvtx,-1.);

  if(drawRatio){
    selection.Form("event:triggered");
    TH1D* htriggered = fCounters.Get("run",selection.Data());
    htriggered->Sumw2();
    hCountForNorm->Divide(htriggered);
  }

  hCountForNorm->DrawClone();
  return hCountForNorm;
}

//___________________________________________________________________________
Int_t AliNormalizationCounter::Multiplicity(AliVEvent* event){

  Int_t multiplicity = 0;
  AliAODEvent *eventAOD = (AliAODEvent*)event;
  AliAODTracklets * aodTracklets = (AliAODTracklets*)eventAOD->GetTracklets();
  Int_t ntracklets = (Int_t)aodTracklets->GetNumberOfTracklets();
  for(Int_t i=0;i<ntracklets; i++){
    Double_t theta = aodTracklets->GetTheta(i);
    Double_t eta = -TMath::Log( TMath::Tan(theta/2.) ); // check the formula
    if(TMath::Abs(eta)<fMultiplicityEtaRange){ // set the proper cut on eta
      multiplicity++;
    }
  }

  return multiplicity;
}

//___________________________________________________________________________
void AliNormalizationCounter::FillCounters(TString name, Int_t runNumber, Int_t multiplicity, Double_t spherocity){


  Int_t sphToInteger=spherocity*fSpherocitySteps;
  if(fMultiplicity  && !fSpherocity) 
    fCounters.Count(Form("Event:%s/Run:%d/Multiplicity:%d",name.Data(),runNumber,multiplicity));
  else if(fMultiplicity  && fSpherocity) 
    fCounters.Count(Form("Event:%s/Run:%d/Multiplicity:%d/Spherocity:%d",name.Data(),runNumber,multiplicity,sphToInteger));
  else if(!fMultiplicity  && fSpherocity) 
    fCounters.Count(Form("Event:%s/Run:%d/Spherocity:%d",name.Data(),runNumber,sphToInteger));
  else 
    fCounters.Count(Form("Event:%s/Run:%d",name.Data(),runNumber));
  return;
}
