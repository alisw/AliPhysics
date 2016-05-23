/*************************************************************************
* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskCorrelation3p.h"
#include "AliThreeParticleCorrelator.h"
#include "AliCorrelation3p.h"
#include "AliCorrelation3p_noQA.h"
#include "AliCFPI0.h"
#include "AliFilteredTrack.h"
#include "AliFilteredEvent.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliOADBContainer.h"
#include "AliCentrality.h"
#include "AliVVZERO.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "THashList.h"
#include "TMath.h"
#include "AliMCParticle.h"
#include "TCanvas.h"
#include <TRandom3.h>
#include <sstream>
#include <memory>
#include<iostream>

ClassImp(AliAnalysisTaskCorrelation3p)
//
//Task to create three particle correlations.
//Authors:
// Matthias Richter
// Paul Baetzing || pbatzing@cern.ch
//

// using std::string;
// using std::cin;
// using std::stringstream;

AliAnalysisTaskCorrelation3p::AliAnalysisTaskCorrelation3p()
  : AliAnalysisTaskSE("AliAnalysisTaskCorrelation3p")
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption("")
  , fstarttime(-1)
  , fCentrality(NULL)
  , fVertexobj(NULL)
  , fRun(130848)
  , fNEventsProcessed(0)
  , fNEventsParsed(0)
  , fNEventsToProcess(-1)
  , fStartAtEvent(1)
  , fVertex()
  , fperiod(AliAnalysisTaskCorrelation3p::P10h)
  , fCollisionType(AliAnalysisTaskCorrelation3p::PbPb)
  , fisESD(kFALSE)
  , fisAOD(kFALSE)
  , fisDstTree(kFALSE)
  , fgenerate(kFALSE)
  , fQA(kFALSE)
  , fqatask(kFALSE)
  , fMoreOutputs(kFALSE)
  , fWeights(NULL)
  , fWeightshpt(NULL)
  , fpTfunction(NULL)
  , fRandom(NULL)
  , fMcArray(NULL)
  , fTrackCuts(NULL)
  , fCorrelator(NULL)
  , fRunNumberList(NULL)
  , fNruns(1000)
  , fRunFillValue(0.0)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fMaxTracksperEvent(-1)
  , fCentralityEstimator("V0M")
  , ftrigger(AliAnalysisTaskCorrelation3p::tracks)
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fBinVer(0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinTriggerPt(3.0)
  , fMaxTriggerPt(8.0)
  , fMinAssociatedPt(1.0)
  , fMaxAssociatedPt(3.0)
  , fMinNClustersTPC(70)
  , fCutMask(0)
  , fMinClusterEnergy(0.3)
  , fMinBCDistance(0.0)
  , fMinNCells(3)
  , fMinM02(0.2)
  , fTOFCutEnabled(1)
  , fTOFCut(100.e-9)
  , fMassInvMean(0.135)
  , fMassInvSigma(0.01)
  , fphospions(kTRUE)
  , femcalpions(kTRUE)
{
  // default constructor
  // 
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-5,-2.5,2.5,5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
  DefineSlots();
}

AliAnalysisTaskCorrelation3p::AliAnalysisTaskCorrelation3p(const char *name, const char* opt)
  : AliAnalysisTaskSE(name)
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption(opt)
  , fstarttime(-1)
  , fCentrality(NULL)
  , fVertexobj(NULL)
  , fRun(130848)
  , fNEventsProcessed(0)
  , fNEventsParsed(0)
  , fNEventsToProcess(-1)
  , fStartAtEvent(1)
  , fVertex()
  , fperiod(AliAnalysisTaskCorrelation3p::P10h)
  , fCollisionType(AliAnalysisTaskCorrelation3p::PbPb)
  , fisESD(kFALSE)
  , fisAOD(kFALSE)
  , fisDstTree(kFALSE)
  , fgenerate(kFALSE)
  , fQA(kFALSE)  
  , fqatask(kFALSE)
  , fMoreOutputs(kFALSE)
  , fWeights(NULL)
  , fWeightshpt(NULL)
  , fpTfunction(NULL)
  , fRandom(NULL)
  , fMcArray(NULL)
  , fTrackCuts(NULL)
  , fCorrelator(NULL)
  , fRunNumberList(NULL)
  , fNruns(1000)
  , fRunFillValue(0.0)
  , fMBinEdges(TArrayD())
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fMaxTracksperEvent(-1)
  , fCentralityEstimator("V0M")
  , ftrigger(AliAnalysisTaskCorrelation3p::tracks)
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fBinVer(0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)
  , fAcceptancecut(0.8)
  , fMinTriggerPt(3.0)
  , fMaxTriggerPt(8.0)
  , fMinAssociatedPt(1.0)
  , fMaxAssociatedPt(3.0)
  , fMinNClustersTPC(70)
  , fCutMask(0)
  , fMinClusterEnergy(0.3)
  , fMinBCDistance(0.0)
  , fMinNCells(3)
  , fMinM02(0.2)
  , fTOFCutEnabled(1)
  , fTOFCut(100.e-9)
  , fMassInvMean(0.135)
  , fMassInvSigma(0.01)
  , fphospions(kTRUE)
  , femcalpions(kTRUE)
  {
  // constructor with options
  //
  //
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-7.5,-2.5,2.5,7.5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
  DefineSlots();
  }

int AliAnalysisTaskCorrelation3p::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  return 0;
}

AliAnalysisTaskCorrelation3p::~AliAnalysisTaskCorrelation3p()
{
  // destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }
  if (fTrackCuts) delete fTrackCuts;
  fTrackCuts=NULL;
  if(fRunNumberList) delete fRunNumberList;
  fRunNumberList=NULL;
  if(fRandom) delete fRandom;
  fRandom = NULL;
}

void AliAnalysisTaskCorrelation3p::UserCreateOutputObjects()
{
  // create result objects and add to output list
  MakeRunNumbers();//Needs to be done once in the beginning
  
  TTimeStamp now;
  fRandom = new TRandom3();
  fRandom->SetSeed(now.GetNanoSec());
  fstarttime = now.GetSec();
  
  TH1::SetDefaultSumw2(kTRUE);//want the collection of weights on all histograms.
  TString collisiontype;
  if(fCollisionType==pp) collisiontype.Append("pp");
  if(fCollisionType==PbPb) collisiontype.Append("PbPb");
  fOutput = new THashList;
  fOutput->SetOwner();
  if(!fQA&&!fqatask){
    if(fMoreOutputs) AliWarning("Creating User Output Objects.");
    //Create the appropriate ThreeParticleCorrelators and add the used one to be fCorrelator.
    AliThreeParticleCorrelator<AliCorrelation3p>* correlator=new AliThreeParticleCorrelator<AliCorrelation3p>;
    fCorrelator=correlator;

    //Initialize QA histograms and add them to fOutput
    InitializeQAhistograms();
    //Intitialize the Multiplicity and ZVertex bins.
    const Int_t    MaxNofEvents=fMaxNEventMix;
    const Int_t    MinNofTracks=fMinNofTracksMix;
    const Int_t    nofMBins=fMBinEdges.GetSize()-1;
    Double_t 	 MBinsTemp[nofMBins+1];
    for(int i=0; i<=nofMBins; ++i) MBinsTemp[i] = fMBinEdges.At(i);
    const Int_t    nofZBins=fZBinEdges.GetSize()-1;//5;
    Double_t 	 ZBinsTemp[nofZBins+1];
    for(int i=0; i<=nofZBins; ++i) ZBinsTemp[i] = fZBinEdges.At(i);
    //Create the AliEventPoolManager 
    AliEventPoolManager* poolMgr = new AliEventPoolManager(MaxNofEvents, MinNofTracks, nofMBins, (Double_t*)MBinsTemp, nofZBins, (Double_t*)ZBinsTemp);
    poolMgr->SetTargetValues(MinNofTracks,1.0E-4,1.0);
    correlator->InitEventMixing(poolMgr);
    correlator->SetNMaxMixed(fMaxTracksperEvent);
    
    //initialize track worker and add to the output if appropriate
    TString tracksname;
    if(ftrigger == AliAnalysisTaskCorrelation3p::tracks) tracksname = Form("tracktrigger_correlation_%.0f_%.0f", fMinTriggerPt, fMaxTriggerPt);
    if(ftrigger == AliAnalysisTaskCorrelation3p::pi0)    tracksname = Form("pi0trigger_correlation_%.0f_%.0f", fMinTriggerPt, fMaxTriggerPt);
    TString triggertype;
    if(ftrigger == AliAnalysisTaskCorrelation3p::tracks) triggertype = "tracks";
    if(ftrigger == AliAnalysisTaskCorrelation3p::pi0) triggertype = "pi0";
    
    TString triggerinit = Form("minTriggerPt=%.1f maxTriggerPt=%.1f minAssociatedPt=%.1f maxAssociatedPt=%.1f collisiontype=%s triggertype=%s", fMinTriggerPt, fMaxTriggerPt, fMinAssociatedPt, fMaxAssociatedPt,collisiontype.Data(),triggertype.Data());
    AliCorrelation3p* workertracktrigger =new AliCorrelation3p(tracksname, fMBinEdges, fZBinEdges);
    workertracktrigger->SetAcceptanceCut(fAcceptancecut);
    workertracktrigger->SetBinningVersion(fBinVer);
    workertracktrigger->Init(triggerinit);
    correlator->Add(workertracktrigger);
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 0));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 1));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 2));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 3));
    fOutput->Add(workertracktrigger);
    }
  if(!fQA&&fqatask){
    //Create the appropriate ThreeParticleCorrelators and add the used one to be fCorrelator.
    AliThreeParticleCorrelator<AliCorrelation3p_noQA>* correlator=new AliThreeParticleCorrelator<AliCorrelation3p_noQA>;
    fCorrelator=correlator;
    //Initialize QA histograms and add them to fOutput
    InitializeQAhistograms();
    //Intitialize the Multiplicity and ZVertex bins.
    const Int_t    MaxNofEvents=fMaxNEventMix;
    const Int_t    MinNofTracks=fMinNofTracksMix;
    const Int_t    nofMBins=fMBinEdges.GetSize()-1;
    Double_t 	 MBinsTemp[nofMBins+1];
    for(int i=0; i<=nofMBins; ++i) MBinsTemp[i] = fMBinEdges.At(i);
    const Int_t    nofZBins=fZBinEdges.GetSize()-1;//5;
    Double_t 	 ZBinsTemp[nofZBins+1];
    for(int i=0; i<=nofZBins; ++i) ZBinsTemp[i] = fZBinEdges.At(i);
    //Create the AliEventPoolManager 
    AliEventPoolManager* poolMgr = new AliEventPoolManager(MaxNofEvents, MinNofTracks, nofMBins, (Double_t*)MBinsTemp, nofZBins, (Double_t*)ZBinsTemp);
    poolMgr->SetTargetValues(MinNofTracks,1.0E-4,1.0);
    correlator->InitEventMixing(poolMgr);
    
    //initialize track worker and add to the output if appropriate
    TString tracksname;
    if(ftrigger == AliAnalysisTaskCorrelation3p::tracks) tracksname = Form("tracktrigger_correlation_%.0f_%.0f", fMinTriggerPt, fMaxTriggerPt);
    if(ftrigger == AliAnalysisTaskCorrelation3p::pi0)    tracksname = Form("pi0trigger_correlation_%.0f_%.0f", fMinTriggerPt, fMaxTriggerPt);
    TString triggertype;
    if(ftrigger == AliAnalysisTaskCorrelation3p::tracks) triggertype = "tracks";
    if(ftrigger == AliAnalysisTaskCorrelation3p::pi0) triggertype = "pi0";
    TString triggerinit = Form("minTriggerPt=%.1f maxTriggerPt=%.1f minAssociatedPt=%.1f maxAssociatedPt=%.1f collisiontype=%s triggertype=%s", fMinTriggerPt, fMaxTriggerPt, fMinAssociatedPt, fMaxAssociatedPt,collisiontype.Data(),triggertype.Data());
    AliCorrelation3p_noQA* workertracktrigger =new AliCorrelation3p_noQA(tracksname, fMBinEdges, fZBinEdges);
    workertracktrigger->SetAcceptanceCut(fAcceptancecut);
    workertracktrigger->SetBinningVersion(fBinVer);
    workertracktrigger->Init(triggerinit);
    correlator->Add(workertracktrigger);
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 0));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 1));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 2));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 3));
    fOutput->Add(workertracktrigger);
    }
  
  if(fQA){
    InitializeQAhistograms();
  }
  // all tasks must post data once for all outputs
  if(fMoreOutputs) AliWarning("Posting once for all outputs.");

  PostData(1, fOutput);
}

void AliAnalysisTaskCorrelation3p::UserExec(Option_t* /*option*/)
{
  if(fgenerate){execgenerate();return;}//Toy MC generator, skips all data processing.
  fNEventsParsed +=1;
//   if(fNEventsParsed<fStartAtEvent) return;
//   if(fNEventsToProcess>0&&fNEventsProcessed>=fNEventsToProcess)return;
  // process the event
  TObject* pInput=InputEvent();
  if (!pInput) {AliError("failed to get input");return;}
  AliVEvent *pEvent = dynamic_cast<AliVEvent*>(pInput);
  if(!pEvent){AliError(Form("input of wrong class type %s, expecting AliVEvent", pInput->ClassName()));return;}
  //In the MC case, get the Array.
//   if(fefficiencies)GetMCArray();
  //if it is not found, return without doing anything:
//   if(fefficiencies&&!fMcArray) return;
  //Find out if it is AOD or ESD.  
  if(!fisDstTree){
    fisESD=pEvent->IsA()==AliESDEvent::Class();
    fisAOD=pEvent->IsA()==AliAODEvent::Class();
  }
  //Get the runnumber and find which bin and fill value this corresponds to.
  fRun = pEvent->GetRunNumber();
  if(fQA){
    TAxis* runnumberaxis= dynamic_cast<TH1D*>(fOutput->FindObject("EventsperRun"))->GetXaxis();
    if (runnumberaxis){double RunBin = runnumberaxis->FindBin(Form("%i",fRun));fRunFillValue = runnumberaxis->GetBinCenter(RunBin);}   
  }
  if(!fisDstTree)GetCentralityAndVertex();
  if(fisDstTree)GetCentralityAndVertex(pEvent);
  if(!SelectEvent()) return;//events are rejected.
  if(fCollisionType==AliAnalysisTaskCorrelation3p::pp) FillHistogram("centVsZVertex",fMultiplicity,fVertex[2]);//only fill with selected events.
  if(fCollisionType==AliAnalysisTaskCorrelation3p::PbPb) FillHistogram("centVsZVertex",fCentralityPercentile,fVertex[2]);
  //initialize the period dependent cuts.

  UsePeriod();  
  if(fQA){
    //Fill Events/run histogram.
    FillHistogram("EventsperRun", fRunFillValue);
    FillHistogram("NEventsVertex",fRunFillValue,fVertex[2]);
    FillHistogram("NEventsCent",fRunFillValue,fCentralityPercentile);
  }
  //To fill with tracks and pions:
  TObjArray* allrelevantParticles = NULL;
  
  if(!fQA){
    allrelevantParticles = new TObjArray();
    allrelevantParticles->SetOwner();//In order for it to work in the event pool.
  }
  fNTriggers=0.0;//Reset fNTriggers
  fNAssociated=0.0;//Reset fNAssociated
  //Fill all the tracks
   GetTracks(allrelevantParticles, pEvent);
//   FillHistogram("centVsNofTracks",fCentralityPercentile,);
  //Fill all the pi0 candidates if appropriate
  if (ftrigger == AliAnalysisTaskCorrelation3p::pi0) GetPi0s(allrelevantParticles, pEvent);
  
  FillHistogram("Ntriggers",fNTriggers);
  FillHistogram("NAssociated",fNAssociated);
  if(fQA){
    FillHistogram("NTriggersperRun",fRunFillValue,fNTriggers);
    FillHistogram("NAssociatedperRun",fRunFillValue,fNAssociated);
  }
  if(fNTriggers>=1)FillHistogram("NAssociatedETriggered",fNAssociated);
  //if fQA the correlations are not build.
  if(!fQA){
    if(fCollisionType==AliAnalysisTaskCorrelation3p::PbPb&&!fqatask){dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p>*>(fCorrelator)->SetEventVzM(fVertex[2],fCentralityPercentile);}
    if(fCollisionType==AliAnalysisTaskCorrelation3p::pp&&!fqatask){dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p>*>(fCorrelator)->SetEventVzM(fVertex[2],fMultiplicity);}
    if(fCollisionType==AliAnalysisTaskCorrelation3p::PbPb&&fqatask){dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p_noQA>*>(fCorrelator)->SetEventVzM(fVertex[2],fCentralityPercentile);}
    if(fCollisionType==AliAnalysisTaskCorrelation3p::pp&&fqatask){dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p_noQA>*>(fCorrelator)->SetEventVzM(fVertex[2],fMultiplicity);}
    //Do the actual correlations.
    if(!((fNTriggers+fNAssociated)==0))fCorrelator->Execute(NULL, allrelevantParticles);//correlate for events that contain at least one trigger or associated.
    else delete allrelevantParticles;
    //Post the output
  }
  fNEventsProcessed +=1;
  if(fMoreOutputs&&((fNEventsProcessed%10000 ==0)||fNEventsProcessed==1)){
    TTimeStamp now;
    int timediff = (now.GetSec()-fstarttime);
    int sec = timediff%60;
    int min = ((timediff-sec)/60)%60;
    int hour = (timediff-sec-60*min)/(60*60);
    if(hour ==0&&min==0){
    AliWarning(Form("Number of Events finished:%i. Number of Events parsed:%i. Time: %i. Time spent: %i seconds",fNEventsProcessed,fNEventsParsed,now.GetTime(),sec));
    }
    if(hour ==0&&min!=0){
    AliWarning(Form("Number of Events finished:%i. Number of Events parsed:%i. Time: %i. Time spent: %i minutes %i seconds",fNEventsProcessed,fNEventsParsed,now.GetTime(),min,sec));
    }
    if(hour !=0&&min!=0){
    AliWarning(Form("Number of Events finished:%i. Number of Events parsed:%i. Time: %i. Time spent: %i hours %i minutes %i seconds",fNEventsProcessed,fNEventsParsed,now.GetTime(),hour,min,sec));
    }
  }
  PostData(1, fOutput);
}

void AliAnalysisTaskCorrelation3p::FinishTaskOutput()
{
  // end of the processing
    TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject("trackCount")) ;
    if (hist) AliWarning(Form("FinishTaskOutput: %i events(s)",(int)hist->GetEntries()));
    TTimeStamp now;
    int timediff = (now.GetSec()-fstarttime);
    int sec = timediff%60;
    int min = ((timediff-sec)/60)%60;
    int hour = (timediff-sec-60*min)/(60*60);    
    if(hour ==0&&min==0){
    AliWarning(Form("Time spent since UserCreateOutputObjects: %i seconds",sec));
    }
    if(hour ==0&&min!=0){
    AliWarning(Form("Time spent since UserCreateOutputObjects: %i minutes %i seconds",min,sec));
    }
    if(hour !=0&&min!=0){
    AliWarning(Form("Time spent since UserCreateOutputObjects: %i hours %i minutes %i seconds",hour,min,sec));}
}

void AliAnalysisTaskCorrelation3p::Terminate(Option_t *)
{
  // last action on the client

}

Int_t AliAnalysisTaskCorrelation3p::GetTracks(TObjArray* allrelevantParticles, AliVEvent *pEvent)
{
  Int_t nofTracks = 0;
  Int_t VZbin;
  if(fWeights){
//     MultBin = fWeights->GetXaxis()->FindBin(fCentralityPercentile);
    VZbin   = fWeights->GetYaxis()->FindBin(fVertex[2]);
  }
  nofTracks=pEvent->GetNumberOfTracks();
  FillHistogram("trackCount",nofTracks);
  for (int i=0; i<nofTracks; i++) {
    Double_t Weight = 1.0;
    AliVParticle* t=pEvent->GetTrack(i);
    if (!t) continue;
    if(fQA){
      FillHistogram("TracksperRun",fRunFillValue);
    }
    FillHistogram("trackUnselectedPt",t->Pt(),1.0);
    FillHistogram("trackUnselectedPhi",t->Phi(),1.0);
    FillHistogram("trackUnselectedTheta",t->Theta(),1.0);
    if (!IsSelected(t)) continue;
    if(fWeights){
      Int_t etabin = fWeights->GetXaxis()->FindBin(t->Eta());
      if(t->Pt()<2.0){
	Int_t pTbin  = fWeights->GetZaxis()->FindBin(t->Pt());
	Weight = fWeights->GetBinContent(etabin,VZbin,pTbin);
      }
      else{
	Weight = fWeightshpt->GetBinContent(etabin,VZbin)*fpTfunction->Eval(t->Pt());
      }
    }
    if(allrelevantParticles){
	AliFilteredTrack * filp = new AliFilteredTrack(*t);
	filp->SetEff(Weight);
      allrelevantParticles->Add(filp);
    }
    if(fQA){
      FillHistogram("selectedTracksperRun",fRunFillValue);
      FillHistogram("NTracksVertex",fRunFillValue,fVertex[2]);
      FillHistogram("NTracksCent",fRunFillValue,fCentralityPercentile);
    }
    FillHistogram("trackPt",t->Pt(),Weight);
    if(fQA&&fWeights&&fWeightshpt){
      if(t->Pt()<4.0){
	FillHistogram("Track_Cent_Vertex_lpT",fCentralityPercentile,fVertex[2],t->Pt());
      }
      else{
	FillHistogram("Track_Cent_Vertex_eta",fCentralityPercentile,fVertex[2]);
      }
    }
    if(dynamic_cast<AliAODTrack*>(t)){
      if(dynamic_cast<AliAODTrack*>(t)->IsGlobalConstrained())FillHistogram("trackPtconstrained",t->Pt(),Weight);
      if(!dynamic_cast<AliAODTrack*>(t)->IsGlobalConstrained())FillHistogram("trackPtnotconstrained",t->Pt(),Weight);

    }
    FillHistogram("trackPhi",t->Phi(),Weight);
    FillHistogram("trackTheta",t->Theta(),Weight);
    if(IsSelectedTrigger(t)){
      fNTriggers+=1;
      FillHistogram("trackTriggerPt",t->Pt(),Weight);
      FillHistogram("trackTriggerPhi",t->Phi(),Weight);
      FillHistogram("trackTriggerTheta",t->Theta(),Weight);
    }
    if(IsSelectedAssociated(t)){
      fNAssociated+=1;
      FillHistogram("trackAssociatedPt",t->Pt(),Weight);
      FillHistogram("trackAssociatedPhi",t->Phi(),Weight);
      FillHistogram("trackAssociatedTheta",t->Theta(),Weight);
    }
  }

  return nofTracks;
}
void AliAnalysisTaskCorrelation3p::GetPi0s(TObjArray* allrelevantParticles, AliVEvent* pEvent)
{
  TObjArray allCaloClusters;
  Int_t nclu = pEvent->GetNumberOfCaloClusters();
  FillHistogram("clusterCount",nclu);
  Int_t ncluPHOS = 0;
  Int_t ncluEMCAL = 0;
  for (Int_t i=0;  i<nclu;  i++) {
    AliVCluster *clu = pEvent->GetCaloCluster(i);
    if (!clu) continue;
    if (!clu->IsPHOS() && !clu->IsEMCAL()) continue;
    if (!GoodCluster(clu)) continue;
    allCaloClusters.Add(clu);
    TLorentzVector lorentzMomentum;
    clu->GetMomentum(lorentzMomentum, fVertex);
    if(clu->IsPHOS()){
      ncluPHOS += 1;
      FillHistogram("PhosClustersperRun",fRunFillValue);
      FillHistogram("clusterPtPHOS",lorentzMomentum.Pt());
      FillHistogram("clusterPhiPHOS",lorentzMomentum.Phi() + 2*TMath::Pi());
      FillHistogram("clusterThetaPHOS",lorentzMomentum.Theta());
    }
    if(clu->IsEMCAL()){
      ncluEMCAL += 1;
      FillHistogram("EmcalClustersperRun",fRunFillValue);
      FillHistogram("clusterPtEMCAL",lorentzMomentum.Pt());
      FillHistogram("clusterPhiEMCAL",lorentzMomentum.Phi());
      FillHistogram("clusterThetaEMCAL",lorentzMomentum.Theta());
    }
  }
  FillHistogram("clusterCountemcal",ncluEMCAL);
  FillHistogram("clusterCountphos",ncluPHOS);
  
  const Int_t nCalo=allCaloClusters.GetEntriesFast() ;
  Int_t npi = 0;
  Int_t npiPHOS = 0;
  Int_t npiEMCAL = 0;
    for(Int_t i1=0; i1 < nCalo-1; i1++){
      AliVCluster *clu1 = (AliVCluster*)allCaloClusters.At(i1);
     for (Int_t i2=i1+1; i2<nCalo; i2++){
      AliVCluster *clu2 = (AliVCluster*)allCaloClusters.At(i2);
	if(!(clu1->IsPHOS() == clu2->IsPHOS())) continue;//only look at pairs from the same detector.
	TLorentzVector ph1, ph2;
	clu1->GetMomentum(ph1, fVertex);
	clu2->GetMomentum(ph2, fVertex);
	TLorentzVector p12  = ph1  + ph2;
	
	Double_t m=p12.M() ;
	if (!(fMassInvMean-fMassInvSigma<m && m<fMassInvMean+fMassInvSigma)) continue;//select only pi0 candidates
	
	npi +=1;
	AliCFPI0 * pi0loc =new AliCFPI0(p12);
	pi0loc->SetIsPhos(clu1->IsPHOS());
	pi0loc->SetIsEmcal(clu1->IsEMCAL());
 	if(pi0loc->IsPHOS()) FillHistogram("PhosPionsperRun",fRunFillValue);
 	if(pi0loc->IsEMCAL()) FillHistogram("EmcalPionsperRun",fRunFillValue);
	
	if (! IsSelected(pi0loc)) continue;
	fNTriggers+=1;
	if(allrelevantParticles) allrelevantParticles->Add(pi0loc);
	if(IsSelectedTrigger(pi0loc))FillHistogram("pi0TriggerPt",pi0loc->Pt());

	if (clu1->IsPHOS()) {
	  npiPHOS +=1;
 	  FillHistogram("PhosSelectedPionsperRun",fRunFillValue);
	  FillHistogram("pi0ptphos",p12.Pt());
	  FillHistogram("pi0phiphos",p12.Phi()+ 2*TMath::Pi());
	  FillHistogram("pi0ThetaPHOS",p12.Theta());
	}
	if (clu1->IsEMCAL()) {
	  npiEMCAL +=1;
	  FillHistogram("EmcalSelectedPionsperRun",fRunFillValue);
	  FillHistogram("pi0ptemcal",p12.Pt());
	  FillHistogram("pi0phiemcal",p12.Phi());
	  FillHistogram("pi0thetaemcal",p12.Theta());
	}
     }
    }
  FillHistogram("pi0count",npi);
  FillHistogram("pi0countEMCAL",npiEMCAL);
  FillHistogram("pi0countPHOS",npiPHOS);
}


void AliAnalysisTaskCorrelation3p::UsePeriod()
{
  bool is2010 = false;
  bool is2011 = false;
  if(fperiod==AliAnalysisTaskCorrelation3p::P10b||fperiod==AliAnalysisTaskCorrelation3p::P10c||fperiod==AliAnalysisTaskCorrelation3p::P10d||fperiod==AliAnalysisTaskCorrelation3p::P10e||fperiod==AliAnalysisTaskCorrelation3p::P10h)is2010 =true;
  if(fperiod==AliAnalysisTaskCorrelation3p::P11a||fperiod==AliAnalysisTaskCorrelation3p::P11h) is2011 = true;
  //initializes track cuts. Checks if the cuts exist, and does not move the pointer otherwise.
  if(!fTrackCuts&&is2010&& fisESD){
      fTrackCuts = AliESDtrackCuts::AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(); 
  }
  if(fTrackCuts&&is2010&& fisESD){
      delete fTrackCuts;
      fTrackCuts = AliESDtrackCuts::AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(); 
  }
  if(!fTrackCuts&&is2011&& fisESD){
      fTrackCuts = AliESDtrackCuts::AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  }
  if(fTrackCuts&&is2011&& fisESD){
    delete fTrackCuts;  
    fTrackCuts = AliESDtrackCuts::AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  }

}

Bool_t AliAnalysisTaskCorrelation3p::IsSelected(AliVParticle* p)
{
  //Performs selection cuts for tracks and triggers
  if (p->IsA()==AliCFPI0::Class()) return IsSelectedTrigger(p);
  if (p->IsA()==AliESDtrack::Class() && IsSelectedTrackESD(p)) return IsSelectedTrigger(p)||IsSelectedAssociated(p);
  if (p->IsA()==AliAODTrack::Class() && IsSelectedTrackAOD(p)) return IsSelectedTrigger(p)||IsSelectedAssociated(p);
  if (p->IsA()==AliFilteredTrack::Class()&& IsSelectedTrackFiltered(p)){dynamic_cast<AliFilteredTrack*>(p)->Calculate();return IsSelectedTrigger(p)||IsSelectedAssociated(p);}
  if (p->IsA()==AliAODMCParticle::Class() && dynamic_cast<AliAODMCParticle*>(p)->IsPhysicalPrimary()) return IsSelectedTrigger(p)||IsSelectedAssociated(p);
  return kFALSE;
}

Bool_t AliAnalysisTaskCorrelation3p::IsSelectedTrigger(AliVParticle* p)
{
  if (dynamic_cast<AliCFPI0*>(p)){if((dynamic_cast<AliCFPI0*>(p)->IsEMCAL()&&!femcalpions)|| (dynamic_cast<AliCFPI0*>(p)->IsPHOS()&&!fphospions))return kFALSE;}
  if (p->Pt()<=fMinTriggerPt) return kFALSE;
  if (fMaxTriggerPt>fMinTriggerPt && p->Pt()>fMaxTriggerPt) return kFALSE;
  float etatrigger=p->Eta();
  if (etatrigger<=-fAcceptancecut || etatrigger>=fAcceptancecut) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrelation3p::IsSelectedAssociated(AliVParticle* p)
{
  if (p->Pt()<=fMinAssociatedPt) return kFALSE;
  if (fMaxAssociatedPt>fMinAssociatedPt && p->Pt()>fMaxAssociatedPt) return kFALSE;
  float etaAssociated=p->Eta();
  if (etaAssociated<=-fAcceptancecut || etaAssociated>=fAcceptancecut) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrelation3p::IsSelectedTrackAOD(AliVParticle* t)
{
  AliAODTrack *AODt = dynamic_cast<AliAODTrack*>(t);
  Bool_t isselected = kTRUE;
  Double_t DCAtang=-999.0;
  Double_t DCAlong=-999.0;
  GetDCA(DCAtang,DCAlong,AODt);
  //Hybrid tracks give flat distributions
  if(fCutMask == 0) isselected = AODt->IsHybridGlobalConstrainedGlobal();
  else if(fCutMask == 1) isselected = AODt->TestFilterBit(BIT(4));
  else if(fCutMask == 2) isselected = AODt->TestFilterBit(BIT(5));
  else if(fCutMask == 3) isselected = AODt->TestFilterBit(BIT(6));
  else isselected = AODt->IsHybridGlobalConstrainedGlobal(); // defaults to global hybrid.
//   if( (AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2))&&isselected)   FillHistogram("TrackDCAandonITSselected",DCAtang,DCAlong,1);
//   if(!(AODt->HasPointOnITSLayer(1)||AODt->HasPointOnITSLayer(2))&&isselected)   FillHistogram("TrackDCAandonITSselected",DCAtang,DCAlong,0);
  return isselected; 
}

Bool_t AliAnalysisTaskCorrelation3p::IsSelectedTrackESD(AliVParticle* t)
{
  return fTrackCuts->IsSelected(t);
}

Bool_t AliAnalysisTaskCorrelation3p::IsSelectedTrackFiltered(AliVParticle* t)
{
  if(dynamic_cast<AliFilteredTrack*>(t)->IsGlobalHybrid()&&(fCutMask==0||fCutMask>5))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(t)->IsBIT4()&&(fCutMask==1))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(t)->IsBIT5()&&(fCutMask==2))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(t)->IsBIT6()&&(fCutMask==3))return kTRUE;
  if((dynamic_cast<AliFilteredTrack*>(t)->IsBIT6()|dynamic_cast<AliFilteredTrack*>(t)->IsBIT5())&&(fCutMask==4))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(t)->IsGlobalHybrid()&&!(dynamic_cast<AliFilteredTrack*>(t)->IsBIT6()|dynamic_cast<AliFilteredTrack*>(t)->IsBIT5())&&(fCutMask==5))return kTRUE;
  return kFALSE;
}

void AliAnalysisTaskCorrelation3p::GetMCArray()
{//Gets the MCarray if we are in AOD.
    fMcArray = 0;
    AliAODInputHandler* aodHandler=dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (aodHandler){
      AliAODEvent *aod=aodHandler->GetEvent();
//       if (aod && fefficiencies) {
// 	fMcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
// 	if (!fMcArray) AliError("Could not retrieve MC array!");
//       }
//       else AliError("Could not retrieve AOD event! MC is only supported in AOD.");
    }
}

void AliAnalysisTaskCorrelation3p::GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt)
{
if(AODt->TestBit(AliAODTrack::kIsDCA)){
  DCAtang = AODt->DCA();
  DCAlong = AODt->ZAtDCA();
}
else{
  if(fVertexobj){
    Double_t fBzkg = dynamic_cast<AliAODEvent*>(InputEvent())->GetMagneticField();
    Double_t* dca = new Double_t[2];
    Double_t* dcacov = new Double_t[3];

    Double_t kVeryBigno = 1000000;
    if(AODt->PropagateToDCA(fVertexobj,fBzkg,kVeryBigno,dca,dcacov)){DCAtang=dca[0];DCAlong = dca[1];}
    else{DCAtang = -999;DCAlong=-999;}
    delete[] dca;
    delete[] dcacov;
    }
  }
} 


void AliAnalysisTaskCorrelation3p::GetCentralityAndVertex()
{
  fVertexobj = InputEvent()->GetPrimaryVertex();
  if(fisAOD)
  {
    // Fill AliAODEvent interface specific information
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(dynamic_cast<AliAODEvent*>(InputEvent())->GetHeader());
    fCentrality =  dynamic_cast<AliCentrality*>(header->GetCentralityP());
    fMultiplicity = header->GetRefMultiplicity();
  }  
  if(fisESD)
  {
    fCentrality = dynamic_cast<AliESDEvent*>(InputEvent())->GetCentrality();
    fMultiplicity = 0.0 ;//not implemented atm
  }
  if(fCentrality)fCentralityPercentile = fCentrality->GetCentralityPercentile(fCentralityEstimator);

  //Get the primary Vertex
  if( fVertexobj ) {
    fVertex[0] = fVertexobj->GetX();
    fVertex[1] = fVertexobj->GetY();
    fVertex[2] = fVertexobj->GetZ();}
  else return;
}

void AliAnalysisTaskCorrelation3p::GetCentralityAndVertex(AliVEvent* pEvent)
{
  //for DstTree:
  if(fCollisionType == pp){fMultiplicity = dynamic_cast<AliFilteredEvent*>(pEvent)->GetCentralityP();}
  if(fCollisionType == PbPb){fCentralityPercentile = dynamic_cast<AliFilteredEvent*>(pEvent)->GetCentralityP();}

  //Get the primary Vertex
  fVertex[0] = dynamic_cast<AliFilteredEvent*>(pEvent)->GetfVertexX();
  fVertex[1] = dynamic_cast<AliFilteredEvent*>(pEvent)->GetfVertexY();
  fVertex[2] = dynamic_cast<AliFilteredEvent*>(pEvent)->GetfVertexZ();
}

Bool_t AliAnalysisTaskCorrelation3p::SelectEvent()
{//This function provides the event cuts for this class.
  if(fisDstTree){
    if(fCollisionType == pp){
    FillHistogram("multiplicity",fMultiplicity,0.75);
    }
    if(fCollisionType == PbPb){
      FillHistogram("centrality",fCentralityPercentile,0.75);
    }
    FillHistogram("vertex",fVertex[2],0.75);
  }
  if(fCollisionType==pp&&!fisDstTree){//With pp, the following cuts are applied:
    FillHistogram("multiplicity",fMultiplicity,0.25);
    FillHistogram("vertex",fVertex[2],0.25);    
    if(!fVertexobj){AliError("Vertex object not found.");return kFALSE;}//Runs only after GetCentralityAndVertex().
    if(fVertexobj->GetNContributors()<1) return kFALSE; // no tracks go into reconstructed vertex
    if(abs(fVertex[2])>fMaxVz) return kFALSE;//Vertex is too far out
    if(fMultiplicity>fMaxNumberOfTracksInPPConsidered) return kFALSE;//Out of multiplicity bounds in pp, no histograms will be filled.
    if(InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.))return kFALSE;  //reject for pileup.
    FillHistogram("multiplicity",fMultiplicity,0.75);
    FillHistogram("vertex",fVertex[2],0.75);    
  }
  if(fCollisionType==PbPb&&!fisDstTree){
    FillHistogram("centrality",fCentralityPercentile,0.25);
    FillHistogram("multiplicity",fMultiplicity,0.25);
    FillHistogram("vertex",fVertex[2],0.25);    
    if(!fVertexobj){AliError("Vertex object not found.");return kFALSE;}//Runs only after GetCentralityAndVertex().
    if(fVertexobj->GetNContributors()<1) return kFALSE; // no tracks go into reconstructed vertex
    if(abs(fVertex[2])>fMaxVz) return kFALSE;//Vertex is too far out
    if(!fCentrality){AliError("Centrality object not found.");return kFALSE;}//Centrality must be defined in the PbPb case.
    if(fCentrality->GetQuality()!=0)return kFALSE;//bad centrality.
    if(fCentralityPercentile<0) return kFALSE;//centrality is not defined
    if(fCentralityPercentile>fMaxMult) return kFALSE;//Out of centrality bounds in PbPb, will not fill any histogram.
    FillHistogram("centrality",fCentralityPercentile,0.75);
    FillHistogram("multiplicity",fMultiplicity,0.75);
    FillHistogram("vertex",fVertex[2],0.75);
  }
  
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrelation3p::GoodCluster(AliVCluster* clu)
{
  //Function to do cuts on clusters.
  Bool_t isgood = kTRUE;
  
  if(clu->IsPHOS()){
    isgood = isgood || clu->E()<fMinClusterEnergy;
    isgood = isgood || clu->GetDistanceToBadChannel()<fMinBCDistance;
    isgood = isgood || clu->GetNCells()<fMinNCells;
    isgood = isgood || clu->GetM02()<fMinM02;
    if(fTOFCutEnabled) isgood = isgood || clu->GetTOF()<fTOFCut;
  }
  
  if(clu->IsEMCAL()){//TODO: Find required cuts on EMCAL photons
  }
  return isgood;
}

void AliAnalysisTaskCorrelation3p::InitializeQAhistograms()
{
  //Function that initializes the QA histograms 
  
  if (!fOutput) return;
  //QA histograms
//   fOutput->Add(new TH3D("TrackDCAandonITS","DCA tangential vs DCA longitudinal vs Is in the first two ITS layers",50,-2,2,50,-5,5,2,-0.5,1.5));
//   fOutput->Add(new TH3D("TrackDCAandonITSselected","DCA tangential vs DCA longitudinal vs Is in the first two ITS layers for selected events",50,-2,2,50,-5,5,2,-0.5,1.5));
//   fOutput->Add(new TH3D("Eventbeforeselection","Vertex vs Multiplicity vs Centrality before event selection.", 50,-15,15,50,0,4000,50,0,100));
//   fOutput->Add(new TH3D("Eventafterselection","Vertex vs Multiplicity vs Centrality after event selection.", 50,-15,15,50,0,4000,50,0,100));
  fOutput->Add(new TH1D("trackCount", "trackCount", 1000,  0, 15000));
  fOutput->Add(new TH1D("trackUnselectedPt"   , "trackPt"   , 1000,  0, 20));
  fOutput->Add(new TH1D("trackPt"   			, "trackPt"   				, 1000,  0, 20));
  fOutput->Add(new TH1D("trackPtconstrained"	   	, "trackPt for tracks constrained to the vertex"   , 1000,  0, 20));
  fOutput->Add(new TH1D("trackPtnotconstrained"  	, "trackPt for tracks not constrained to the vertex"   , 1000,  0, 20));
  fOutput->Add(new TH1D("trackAssociatedPt" , "Pt of associated Tracks", 1000, fMinAssociatedPt, fMaxAssociatedPt));
  fOutput->Add(new TH1D("trackTriggerPt" , "Pt of Trigger Tracks", 1000, fMinTriggerPt, fMaxTriggerPt));
  fOutput->Add(new TH1D("trackUnselectedPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackTriggerPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackAssociatedPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackUnselectedTheta", "trackTheta",  180, 0, TMath::Pi()));
  fOutput->Add(new TH1D("trackTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("trackTriggerTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("trackAssociatedTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("Ntriggers","Number of triggers per event",50,-0.5,49.5));
  fOutput->Add(new TH1D("NAssociated","Number of Associated per event",200,-0.5,199.5));
  fOutput->Add(new TH1D("NAssociatedETriggered","Number of Associated per event that contains a trigger.",200,-0.5,199.5));
  if(fWeights)fOutput->Add(fWeights);
  if(fWeightshpt)fOutput->Add(fWeightshpt);
  if(fpTfunction)fOutput->Add(fpTfunction);

  //   if (ftrigger == AliAnalysisTaskCorrelation3p::pi0 || ftrigger == AliAnalysisTaskCorrelation3p::pi0MC){
//     fOutput->Add(new TH1D("clusterCount", "clusterCount", 1000,  0, 2000));
//     fOutput->Add(new TH1D("clusterCountphos", "clusterCountphos", 1000,  0, 2000));
//     fOutput->Add(new TH1D("clusterCountemcal", "clusterCountemcal", 1000,  0, 2000));
//     fOutput->Add(new TH1D("clusterPtPHOS"   , "clusterPtPHOS"   , 1000,  0, 2000));
//     fOutput->Add(new TH1D("clusterPtEMCAL"   , "clusterPtEMCAL"   , 1000,  0, 2000));
//     fOutput->Add(new TH1D("clusterPhiPHOS"  , "clusterPhiPHOS"  ,  180,  0., 2*TMath::Pi()));
//     fOutput->Add(new TH1D("clusterPhiEMCAL"  , "clusterPhiEMCAL"  ,  180,  0., 2*TMath::Pi()));
//     fOutput->Add(new TH1D("clusterThetaPHOS", "clusterThetaPHOS",  180, -1.*TMath::Pi(), TMath::Pi()));
//     fOutput->Add(new TH1D("clusterThetaEMCAL", "clusterThetaEMCAL",  180, -1.*TMath::Pi(), TMath::Pi()));
//     fOutput->Add(new TH1D("pi0count", "pi0count", 1000, 0, 2000));
//     fOutput->Add(new TH1D("pi0countPHOS", "pi0countPHOS", 1000, 0, 2000));
//     fOutput->Add(new TH1D("pi0countEMCAL", "pi0countEMCAL", 1000, 0, 2000));
//     fOutput->Add(new TH1D("pi0ptphos","pi0ptphos", 1000, 0, 2000));  
//     fOutput->Add(new TH1D("pi0phiphos","pi0phiphos",180, 0.,2*TMath::Pi()));
//     fOutput->Add(new TH1D("pi0ThetaPHOS", "pi0ThetaPHOS",180, -1.*TMath::Pi(), TMath::Pi()));
//     fOutput->Add(new TH1D("pi0ptemcal", "pi0ptemcal", 1000,0,2000));  
//     fOutput->Add(new TH1D("pi0phiemcal", "pi0phiemcal", 180, 0., 2*TMath::Pi()));
//     fOutput->Add(new TH1D("pi0thetaemcal", "pi0thetaemcal",180,-1.*TMath::Pi(), TMath::Pi()));
//     fOutput->Add(new TH1D("pi0TriggerPt","Pt of Trigger pi0s", 1000, fMinTriggerPt,fMaxTriggerPt));
//   }
//   fOutput->Add(new TH1D("vzeroMult" , "V0 Multiplicity",  200,  0, 30000));
  if(fCollisionType==PbPb)fOutput->Add(new TH2D("centrality", "Centrality before and after selection",  100,  0, 100,2,0,1));
  if(!fisDstTree||(fCollisionType==pp))fOutput->Add(new TH2D("multiplicity", "Multiplicity of tracks,  before and after selection",  100,  0, fMaxNumberOfTracksInPPConsidered,2,0,1));
  fOutput->Add(new TH2D("vertex", "Vertex of tracks,selected vs unselected",  100,  -15.0, 15.0,2,0,1));
  if(fCollisionType==PbPb)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, 100, 100, -10, 10));
  if(fCollisionType==pp)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, fMaxNumberOfTracksInPPConsidered, 100, -10, 10));

  if(fQA){
    //QA per run histograms:
    TH1D * eventsperrun 	= new TH1D("EventsperRun", "# Events per Run", fNruns, 0, 1);
    TH1D * TracksperRun 	= new TH1D("TracksperRun", "# tracks per Run", fNruns, 0,1);
    TH1D * selectedTracksperRun = new TH1D("selectedTracksperRun", "# selected tracks per Run", fNruns, 0,1);
    TH2D * NTriggersperRun 	= new TH2D("NTriggersperRun","# triggers per event per Run",fNruns, 0,1,50,-0.5,49.5);
    TH2D * NAssociatedperRun 	= new TH2D("NAssociatedperRun","# associated per event per Run",fNruns, 0,1,100,-0.5,99.5);
    TH2D * NTracksVertex	= new TH2D("NTracksVertex","#selected tracks per run and vertex",fNruns,0,1,100,-10.0,10.0);
    TH2D * NEventsVertex	= new TH2D("NEventsVertex","Events per run and vertex",fNruns,0,1,100,-10.0,10.0);
    TH2D * NTracksCent	   	= new TH2D("NTracksCent","#selected tracks per run and vertex",fNruns,0,1,100,0.0,100.0);
    TH2D * NEventsCent	   	= new TH2D("NEventsCent","Events per run and vertex",fNruns,0,1,100,0.0,100.0);
    
    for(int i=0; i<fNruns; i++){
      TString lable = Form("%i",fRunNumberList[i]);
      eventsperrun->GetXaxis()->SetBinLabel(i+1, lable);
      eventsperrun->GetXaxis()->LabelsOption("v");
      TracksperRun->GetXaxis()->SetBinLabel(i+1, lable);
      TracksperRun->GetXaxis()->LabelsOption("v");
      selectedTracksperRun->GetXaxis()->SetBinLabel(i+1, lable);
      selectedTracksperRun->GetXaxis()->LabelsOption("v");
      NTriggersperRun->GetXaxis()->SetBinLabel(i+1,lable);
      NTriggersperRun->GetXaxis()->LabelsOption("v");
      NAssociatedperRun->GetXaxis()->SetBinLabel(i+1,lable);
      NAssociatedperRun->GetXaxis()->LabelsOption("v");    
      NTracksVertex->GetXaxis()->SetBinLabel(i+1,lable);
      NTracksVertex->GetXaxis()->LabelsOption("v"); 
      NEventsVertex->GetXaxis()->SetBinLabel(i+1,lable);
      NEventsVertex->GetXaxis()->LabelsOption("v"); 
      NTracksCent->GetXaxis()->SetBinLabel(i+1,lable);
      NTracksCent->GetXaxis()->LabelsOption("v"); 
      NEventsCent->GetXaxis()->SetBinLabel(i+1,lable);
      NEventsCent->GetXaxis()->LabelsOption("v"); 
    }
    fOutput->Add(eventsperrun);
    fOutput->Add(TracksperRun);
    fOutput->Add(selectedTracksperRun);
    fOutput->Add(NTriggersperRun);
    fOutput->Add(NAssociatedperRun);
    fOutput->Add(NTracksVertex);
    fOutput->Add(NEventsVertex);
    fOutput->Add(NTracksCent);
    fOutput->Add(NEventsCent);
  }
  if(fQA&&fWeights&&fWeightshpt){
    TH3D * histtrackslpt = (TH3D*)(fWeights->Clone("Track_Cent_Vertex_lpT"));
    histtrackslpt->Reset();
    histtrackslpt->SetTitle("Tracks in Centrality vs Vertex vs pT");
    fOutput->Add(histtrackslpt);
    TH2D * histtrackshpt = (TH2D*)(fWeightshpt->Clone("Track_Cent_Vertex_eta"));
    histtrackshpt->Reset();    
    histtrackshpt->SetTitle("Tracks in Centrality vs Vertex");
    fOutput->Add(histtrackshpt);    
  }
//   if (ftrigger == AliAnalysisTaskCorrelation3p::pi0 || ftrigger == AliAnalysisTaskCorrelation3p::pi0MC){
//     TH1D * PhosClustersperRun = new TH1D("PhosClustersperRun", "# clusters in Phos per Run", fNruns, 0, 1);
//     TH1D * PhosPionsperRun = new TH1D("PhosPionsperRun", "# Pi0s in Phos per Run.", fNruns, 0, 1);
//     TH1D * PhosSelectedPionsperRun = new TH1D("PhosSelectedPionsperRun", "# selected Pi0 triggers in Phos per Run.", fNruns, 0, 1);
//     TH1D * EmcalClustersperRun = new TH1D("EmcalClustersperRun", "# clusters in Emcal per Run.", fNruns, 0,1);
//     TH1D * EmcalPionsperRun = new TH1D("EmcalPionsperRun", "# Pi0s in Emcal per Run.", fNruns, 0,1);
//     TH1D * EmcalSelectedPionsperRun = new TH1D("EmcalSelectedPionsperRun", "# selected Pi0 triggers in Emcal per Run.", fNruns, 0,1);
//     for(int i=0; i<fNruns; i++){
//       TString lable = Form("%i",fRunNumberList[i]);
// 
//       PhosClustersperRun->GetXaxis()->SetBinLabel(i+1, lable);
//       PhosClustersperRun->GetXaxis()->LabelsOption("v");
//       PhosPionsperRun->GetXaxis()->SetBinLabel(i+1, lable);
//       PhosPionsperRun->GetXaxis()->LabelsOption("v");
//       PhosSelectedPionsperRun->GetXaxis()->SetBinLabel(i+1, lable);
//       PhosSelectedPionsperRun->GetXaxis()->LabelsOption("v");
//       EmcalClustersperRun->GetXaxis()->SetBinLabel(i+1, lable);
//       EmcalClustersperRun->GetXaxis()->LabelsOption("v");
//       EmcalPionsperRun->GetXaxis()->SetBinLabel(i+1, lable);
//       EmcalPionsperRun->GetXaxis()->LabelsOption("v");
//       EmcalSelectedPionsperRun->GetXaxis()->SetBinLabel(i+1, lable);
//       EmcalSelectedPionsperRun->GetXaxis()->LabelsOption("v");
//       }
//     fOutput->Add(PhosClustersperRun);
//     fOutput->Add(PhosPionsperRun);
//     fOutput->Add(PhosSelectedPionsperRun);
//     fOutput->Add(EmcalClustersperRun);
//     fOutput->Add(EmcalPionsperRun);
//     fOutput->Add(EmcalSelectedPionsperRun);}
}

void AliAnalysisTaskCorrelation3p::InitializeEffHistograms()
{
  //Function that initializes the Efficiency histogram
  Int_t    nofMBins=2*(fMBinEdges.GetSize()-1);
  Double_t MBinmined=fMBinEdges.At(0);
  Double_t MBinmaxed=fMBinEdges.At(fMBinEdges.GetSize()-1);
  Int_t    nofZBins=2*(fZBinEdges.GetSize()-1);
  Double_t ZBinmined=fZBinEdges.At(0);
  Double_t ZBinmaxed=fZBinEdges.At(fZBinEdges.GetSize()-1);
  Int_t nphi = 72;
  Double_t phimin = 0.0;
  Double_t phimax = 2.0*TMath::Pi();
  Int_t nEta = 126;
  Double_t EtaMin = -0.9;
  Double_t EtaMax =  0.9;
  Double_t pTmin = fMinAssociatedPt;
  Double_t pTmax = fMaxTriggerPt;
  Int_t npT = (pTmax-pTmin)/0.1 - 0.5;//rounding up

  Int_t  bins[5]   = {nofMBins, nofZBins,nphi,nEta,npT};
  Double_t xmin[5] = {MBinmined,ZBinmined,phimin,EtaMin,pTmin};
  Double_t xmax[5] = {MBinmaxed,ZBinmaxed,phimax,EtaMax,pTmax};

  fOutput->Add(new THnD("hnTracksinBins","Tracks in different bins.",5,bins,xmin,xmax));
//  dynamic_cast<THnD*>(fOutput->FindObject("hnTracksinBins"))->Sumw2();
  fOutput->Add(new THnD("hnTracksinBinsRecPP","Tracks in different bins from reconstruction that originate from a PhysicalPrimary MC particle.",5,bins,xmin,xmax));
//  dynamic_cast<THnD*>(fOutput->FindObject("hnTracksinBinsRecPP"))->Sumw2();  
  fOutput->Add(new THnD("hnTracksinBinsMC","Tracks in different bins, MC truth for charged particles.",5,bins,xmin,xmax));
//  dynamic_cast<THnD*>(fOutput->FindObject("hnTracksinBinsMC"))->Sumw2();  
}


void AliAnalysisTaskCorrelation3p::MakeRunNumbers()
{
  //Initialize array for run numbers.
  Int_t runnumbersP10b[fNRunsP10b] = {117222, 117220, 117116, 117112, 117109, 117099, 117092, 117063, 117060, 117059, 117053, 117052, 117050, 117048, 116787, 116645, 116643, 116574, 116571, 116562, 116432, 116431, 116429,116403, 116402, 116372, 116360, 116358, 116288, 116102, 116081, 116079, 115521, 115414, 115406, 115401, 115399, 115393, 115369,115345, 115335, 115328,  115327, 115322, 115318, 115312, 115310,  115193, 115186, 115056, 114931, 114930, 114924, 114920, 114918, 114798, 114786};
  Int_t runnumbersP10c[fNRunsP10c] = {121040, 121039, 120829, 120825, 120824, 120823, 120822, 120821, 120820, 120758, 120750, 120741, 120671, 120617, 120616, 120505, 120504, 120503, 120244,120079, 120076, 120073, 120072, 120069, 120067, 119862, 119859, 119856, 119853, 119849, 119846, 119845, 119844, 119842, 119841, 119163, 119161, 119159, 118561, 118560, 118558, 118556, 118518, 118512, 118507, 118506};
  Int_t runnumbersP10d[fNRunsP10d] = {126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126168,126167, 126160, 126158, 126097, 126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843,  125842, 125633, 125632, 125630, 125628, 125296, 125295, 125186, 125156, 125140, 125139, 125134, 125133, 125101, 125100, 125097, 125085, 125083, 125023, 124751, 122375, 122374};
  Int_t runnumbersP10e[fNRunsP10e] = {130850, 130848, 130847, 130844, 130842, 130840, 130834, 130804, 130803, 130802, 130799, 130798, 130795, 130793, 130704, 130696, 130628, 130623, 130621, 130620, 130609, 130608, 130601, 130526, 130524, 130520, 130519, 130517, 130481, 130480, 130479, 130375, 130360, 130358, 130356, 130354, 130343, 130342, 130178, 130172, 130168, 130158, 130157, 130151, 130149, 129983, 129966, 129962, 129961, 129960, 129959, 129744, 129742, 129738, 129736, 129735, 129734, 129729, 129726, 129725, 129723, 129666, 129659, 129653, 129652, 129651, 129650, 129647, 129641, 129639, 129599, 129587, 129586, 129540, 129536, 129528, 129527, 129525, 129524, 129523, 129521, 129520, 129519, 129516, 129515, 129514, 129513, 129512, 129042, 128913, 128855, 128853, 128850, 128843, 128836, 128835, 128834, 128833, 128824, 128823, 128820, 128819, 128778, 128777, 128678, 128677, 128621, 128615, 128611, 128609, 128605, 128596, 128594, 128592, 128590, 128582, 128506, 128505, 128504, 128503, 128498, 128495, 128494, 128486, 128452, 128366}; 
  Int_t runnumbersP10h[fNRunsP10h] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135};
  Int_t runnumbersP11a[fNRunsP11a] = {146860, 146859, 146858, 146856, 146824, 146817, 146807, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746, 146402, 146369,146292, 146287, 146282, 146277, 146273, 146272, 146223, 146220, 146208, 146158, 146156, 146153, 146152, 146148, 146147, 146141, 146099, 146079, 146072, 146071, 146027, 146026, 146025, 146024, 146023, 145674, 145455, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145300, 145292, 145290, 145289, 145288};
  Int_t runnumbersP11h[fNRunsP11h] = {170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915};

  if (fperiod == AliAnalysisTaskCorrelation3p::P10b){
    fNruns = fNRunsP10b;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10b[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::pp;
  }
  else if (fperiod == AliAnalysisTaskCorrelation3p::P10c){
    fNruns = fNRunsP10c;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10c[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::pp;
  }
  else if (fperiod == AliAnalysisTaskCorrelation3p::P10d){
    fNruns = fNRunsP10d;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10d[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::pp;
  }
  else if (fperiod == AliAnalysisTaskCorrelation3p::P10e){
    fNruns = fNRunsP10e;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10e[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::pp;
  }
  else if (fperiod == AliAnalysisTaskCorrelation3p::P11a){
    fNruns = fNRunsP11a;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP11a[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::pp;
  }
  else if (fperiod ==  AliAnalysisTaskCorrelation3p::P10h){
    fNruns = fNRunsP10h;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10h[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::PbPb;
  }
  else if (fperiod == AliAnalysisTaskCorrelation3p::P11h){
    fNruns = fNRunsP11h;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP11h[i];
    //Set the correct collision type
    fCollisionType = AliAnalysisTaskCorrelation3p::PbPb;
  }
  if(fCollisionType==pp)for(int i=0;i<fMBinEdges.GetSize();i++){
      fMBinEdges.AddAt(fMaxNumberOfTracksInPPConsidered*fMBinEdges.At(i)/fMBinEdges.At(fMBinEdges.GetSize()-1),i);
    }
  return ;
}

void AliAnalysisTaskCorrelation3p::SetMixingScheme(Int_t MaxNEventMix, Int_t MinNofTracksMix, TArrayD MBinEdges, TArrayD ZBinEdges)
{
  fMaxNEventMix= MaxNEventMix;
  fMinNofTracksMix = MinNofTracksMix;
  for(int i=0; i<MBinEdges.GetSize()-1; ++i)
    if(MBinEdges.At(i) > MBinEdges.At(i+1)) AliFatal("edges are not sorted");
  for(int i=0; i<ZBinEdges.GetSize()-1; ++i)
    if(ZBinEdges.At(i) > ZBinEdges.At(i+1)) AliFatal("edges are not sorted");  
  fMBinEdges = MBinEdges;
  fZBinEdges = ZBinEdges;
  fMaxMult = fMBinEdges.At(fMBinEdges.GetSize()-1);
}

void AliAnalysisTaskCorrelation3p::FillHistogram(const char* key, Double_t x)
{
  TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p::FillHistogram(const char* key, Double_t x, Double_t y)
{
  TH2 * hist = dynamic_cast<TH2*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y) ;
  else if(dynamic_cast<TH1*>(fOutput->FindObject(key)))dynamic_cast<TH1*>(fOutput->FindObject(key))->Fill(x,y);
  else AliError(Form("can not find histogram (of instance TH2 or TH1) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z)
{
  TH3 * hist = dynamic_cast<TH3*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y,z) ;
  else if(dynamic_cast<TH2*>(fOutput->FindObject(key)))dynamic_cast<TH2*>(fOutput->FindObject(key))->Fill(x,y,z);
  else AliError(Form("can not find histogram (of instance TH3 or TH2) <%s> ",key)) ;
}

void AliAnalysisTaskCorrelation3p::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z,Double_t a, Double_t b)
{
  THnD * hist = dynamic_cast<THnD*>(fOutput->FindObject(key)) ;
  if(hist){
    Double_t s[5] = {x,y,z,a,b};
    hist->Fill(s) ;
  }
  else AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}
void AliAnalysisTaskCorrelation3p::execgenerate()
{//This function is called once for every event, and generates the particle array, and then executes the 3p correlations. Should be run locally.
  delete gRandom;
  gRandom=fRandom;
  Int_t nevents;
  Int_t NJetPairs; //number of jet pairs per event
  Int_t Njetpart_n; // number of particles in near-side jet
  Double_t nearwidth;//gaussian width of the near side jet in radians.
  Int_t Njetpart_a; // number of particles in away-side jet
  Double_t awaywidth;//gaussian width of the away side jet in radians.  
  Int_t NJetTriplets; // number of jet triplets per event
  Int_t Njetpart_n_3; //number of particles in near-side jet
  Double_t nearwidth_3;//gaussian width of the near side jet in radians.
  Int_t Njetpart_a_3;// number of particles in away-side jet
  Double_t awaywidth_3;//gaussian width of the away side jet in radians.
  Double_t Trijetawaysplitting;//splitting of the away side jets in jet triplets in radians.
  Double_t NFLOWparticles;//number of particles from flow.
  Double_t v2; // v2
  Double_t v3; // v3
  Int_t NFlat; // background particles.
  

  
  TFile* settingstest =TFile::Open("generatorsetting.root","READ");
  if(!settingstest){
    //File does not exist yet, make one
    Askforgensettings();
  }
  else settingstest->Close();
  //file exists now, read existing settings.
  TFile* settings =TFile::Open("generatorsetting.root","READ");
  TDirectory * dir = settings->GetDirectory("");
  TTree* tree;
  dir->GetObject("Parameters",tree);
  tree->SetBranchAddress("nevents",&nevents);
  tree->SetBranchAddress("NJetPairs",&NJetPairs);
  tree->SetBranchAddress("Njetpart_n",&Njetpart_n);
  tree->SetBranchAddress("nearwidth",&nearwidth);
  tree->SetBranchAddress("Njetpart_a",&Njetpart_a);
  tree->SetBranchAddress("awaywidth",&awaywidth);
  tree->SetBranchAddress("NJetTriplets",&NJetTriplets);
  tree->SetBranchAddress("Njetpart_n_3",&Njetpart_n_3);
  tree->SetBranchAddress("nearwidth_3",&nearwidth_3);
  tree->SetBranchAddress("Njetpart_a_3",&Njetpart_a_3);
  tree->SetBranchAddress("awaywidth_3",&awaywidth_3);
  tree->SetBranchAddress("Trijetawaysplitting",&Trijetawaysplitting);
  tree->SetBranchAddress("NFLOWparticles",&NFLOWparticles);
  tree->SetBranchAddress("v2",&v2);
  tree->SetBranchAddress("v3",&v3);
  tree->SetBranchAddress("NFlat",&NFlat);
  tree->GetEntry(0);
  TF1* triggerpT = dynamic_cast<TF1*>(settings->Get("triggerpT"));
  TF1* associatedpT = dynamic_cast<TF1*>(settings->Get("associatedpT"));
  TGraphAsymmErrors * v2graph = dynamic_cast<TGraphAsymmErrors*>(settings->Get("v2ptdist"));
  TGraphAsymmErrors * v3graph = dynamic_cast<TGraphAsymmErrors*>(settings->Get("v3ptdist"));
  settings->Close();
  if(dynamic_cast<TH1*>(fOutput->FindObject("trackCount"))->GetEntries()==0&&!(NJetPairs==0&&NJetTriplets==0&&NFLOWparticles==0&&NFlat==0)){
    AliWarning(Form("This will now generate %i events where each event contains the following:",11*nevents));
    if((NJetPairs!=0&&Njetpart_a==0)||(NJetTriplets!=0&&Njetpart_a_3==0)){
      if(NFLOWparticles ==0 && NFlat!=0){
	AliWarning(Form("    1 trigger and %i associated particles, uniformly distributed in the phase space covered.",NFlat));
      }
      if(NFLOWparticles !=0 && NFlat==0 && (v2!=0||v3!=0)){
	if(v2!=0)AliWarning(Form("    1 trigger with a v2 of %f and %i associated particles from flow ",v2*20.0/100.0,(int)NFLOWparticles));
	if(v2==0)AliWarning(Form("    1 uniformly distributed trigger and %i associated particles from flow ",(int)NFLOWparticles));
	if(v2!=0&&v3!=0) AliWarning(Form("with a v2 of %f and a v3 of %f",v2,v3));
	if(v2!=0&&v3==0) AliWarning(Form("with a v2 of %f", v2 ));
	if(v2==0&&v3!=0) AliWarning(Form("with a v3 of %f",v3));
      }
      if(NFLOWparticles !=0 && NFlat!=0){
	if(v2!=0)AliWarning(Form("    1 trigger with a v2 of %f,%i uniformly distributed (background) associated particles and %i associated particles from flow ",v2*20.0/100.0,NFlat,(int)NFLOWparticles));
	if(v2==0)AliWarning(Form("    1 uniformly distributed trigger ,%i uniformly distributed (background) associated particles and %i associated particles from flow ",NFlat,(int)NFLOWparticles));
	if(v2!=0&&v3!=0) AliWarning(Form("with a v2 of %f and a v3 of %f",v2,v3));
	if(v2!=0&&v3==0) AliWarning(Form("with a v2 of %f",v2));
	if(v2==0&&v3!=0) AliWarning(Form("with a v3 of %f",v3));	
      }
    }
    if(NJetPairs!=0&&Njetpart_a!=0){
      AliWarning(Form("    %i correlated Dijets, where the nearside contains 1 trigger and %i associated, spread around the jet axis with a width of %f radians,",NJetPairs,Njetpart_n-1,nearwidth));
      AliWarning(Form("    and the away side contains %i associated particles, spread around the jet axis with a width of %f radians.",Njetpart_a, awaywidth));
    }
    if(NJetTriplets!=0&&Njetpart_a_3!=0){
      AliWarning(Form("    %i correlated Jet triplets, where the nearside contains 1 trigger and %i associated, spread around the jet axis with a width of %f radians,",NJetTriplets,Njetpart_n_3-1,nearwidth_3));
      AliWarning(Form("    and the away side jets contain %i associated particles each, spread around the jet axis with a width of %f radians.",Njetpart_a_3,awaywidth));
      AliWarning(Form("    the away side jets are split in angular distance by %f radians.",Trijetawaysplitting));
    }
    if(NFLOWparticles !=0 &&!((NJetPairs!=0&&Njetpart_a==0)||(NJetTriplets!=0&&Njetpart_a_3==0))){
      AliWarning(Form("   %f associated particles from flow ",Trijetawaysplitting));
      if(v2!=0&&v3!=0) AliWarning(Form("with a v2 of %f and a v3 of %f",v2 ,v3));
      if(v2!=0&&v3==0) AliWarning(Form("with a v2 of %f", v2));
      if(v2==0&&v3!=0) AliWarning(Form("with a v3 of %f", v3));	
    }
    if(NFlat !=0 &&!((NJetPairs!=0&&Njetpart_a==0)||(NJetTriplets!=0&&Njetpart_a_3==0))){
      AliWarning(Form("    %i associated particles, uniformly distributed in the phase space covered.",NFlat));
    }
  }
  Double_t Pii = TMath::Pi();
  TObjArray allrelevantParticles;
  int nTrig = 0;
  int nAss = 0;
  int ndijets=0;
  int ntrijets=0;
  int nFlatP=0;
  int nFlowP=0;
  int nTotal = 0;
  for(Int_t iev=0;iev<nevents;iev++){
    int evnr = dynamic_cast<TH1*>(fOutput->FindObject("trackCount"))->GetEntries();
    if(evnr%100==0) AliWarning(Form("event no: %i",evnr));
    allrelevantParticles.Clear();
    nTrig = 0;
    nAss = 0;
    nTotal = 0;
    //Pick a reaction plane. If there is flow, this correlates directions of dijets and the flow a bit, if not everything is flat.
    Double_t reacplane = TMath::TwoPi()*gRandom->Rndm();
    TF1 *PhiDist = new TF1("PhiDist","1+2*[1]*cos(2*(x-[0]))+2*[2]*cos(3*(x-[0]))",0,TMath::TwoPi());
    PhiDist->SetParameter(0,reacplane);    
    //generate NJetPairs dijets
    ndijets = 0;
    while(ndijets<NJetPairs){
      if(v2graph&&!triggerpT) v2 = v2graph->Eval(7.0);//No pT dist, all triggers are 7.0
      if(v2graph&&triggerpT) v2 = v2graph->Eval(0.5*triggerpT->GetXmax()+0.5*triggerpT->GetXmin());
      if(v2 > 0) PhiDist->SetParameter(1,v2);
      else PhiDist->SetParameter(1,0.0);
      if(v3graph&&!triggerpT) v3 = v3graph->Eval(7.0);//No pT dist, all triggers are 7.0
      if(v3graph&&triggerpT)  v3 = v3graph->Eval(0.5*triggerpT->GetXmax()+0.5*triggerpT->GetXmin());
      if(v3 > 0) PhiDist->SetParameter(2,v3);
      else PhiDist->SetParameter(2,0.0);
      // Choosing random jet axis azimuthal angles and near- and away-side
      // pseudorapidities (flat distributions)
      Double_t jetphi = PhiDist->GetRandom(); //slight correlation between reaction plane and jet direction.
//       jetphi = reacplane;
      Double_t jeteta_near = 2.0*(gRandom->Rndm()-0.5);
      Double_t jeteta_away = 2.0*(gRandom->Rndm()-0.5);
      Short_t charge = 1;
      if(gRandom->Rndm()>0.5) charge *=-1;
      Short_t nearcharge = charge;
      Short_t awaycharge = -1*charge;//total charge of the away side jet 0 or opposite of the near side jet.

      // Generate dijet.  No transverse momentum effects considered,
      // only direction (azimuthal angle and pseudorapidity)
      for (Int_t ipart=0; ipart<Njetpart_n; ipart++){
	Double_t PT;//pT of the particle, if there is a distribution, use it.
	if(triggerpT&&associatedpT){
	  if(ipart == 0)PT = triggerpT->GetRandom();
	  else PT = associatedpT->GetRandom();
	}
	else{
	  if(ipart == 0)PT = 7.0;
	  else PT = 3.0;
	}
	Double_t pxpart, pypart, pzpart;
	Double_t incldif = gRandom->Gaus(0.0,nearwidth); // Gaussian width
	Double_t angledif = TMath::Pi()*(gRandom->Rndm()-0.5);
	Double_t near_phidif = incldif*TMath::Cos(angledif);// gRandom->Gaus(0.0,nearwidth); // Gaussian width
	Double_t near_thetadif =incldif*TMath::Sin(angledif);// gRandom->Gaus(0.0,nearwidth);
	Double_t phipart = jetphi + near_phidif;
	Double_t thetapart = 2.0*TMath::ATan(TMath::Exp(-1.0*jeteta_near))+ near_thetadif;
	while((phipart<(-0.5*Pii))||(phipart>(1.5*Pii))){
	  if (phipart<(-0.5*Pii)) phipart += TMath::TwoPi();
	  if (phipart>(1.5*Pii)) phipart -= TMath::TwoPi();
	}
	while((thetapart<0)||(thetapart>Pii)){
	  if (thetapart>(Pii)) thetapart-=TMath::TwoPi();//if it is in the range pi->2pi, it is moved to -pi->0
	  if (thetapart<(0.0)&&thetapart>-Pii) thetapart*=-1;//if it is in the range -pi->0, it is moved into 0->pi
	  if (thetapart<(0.0)&&thetapart<-Pii) thetapart+=TMath::TwoPi();	//if in range -2pi->-pi, it is moved into 0->pi
	}
	Double_t etapart = -TMath::Log(TMath::Tan(thetapart/2.0));
	Double_t P = PT*TMath::CosH(etapart);//pT*cosh(eta)	
	TParticle* mypart = new TParticle();
	nTotal +=1;
	pxpart = P*TMath::Sin(thetapart)*TMath::Cos(phipart);
	pypart = P*TMath::Sin(thetapart)*TMath::Sin(phipart);
	pzpart = P*TMath::Cos(thetapart);
	mypart->SetMomentum(pxpart,pypart,pzpart,pxpart*pxpart+pypart*pypart+pzpart*pzpart);
	mypart->SetPdgCode(nearcharge*11);//electron if charge is positive, positron if negative.
	nearcharge *= -1; //Flip charge for each particle (neutral or one charge in the jet).
	AliMCParticle* myalivpart= new AliMCParticle(mypart);
	FillHistogram("trackPt",myalivpart->Pt());
	FillHistogram("trackPhi",myalivpart->Phi());
	FillHistogram("trackTheta",myalivpart->Theta());
	if(IsSelectedTrigger(myalivpart)){
	  FillHistogram("trackTriggerPt",myalivpart->Pt());
	  FillHistogram("trackTriggerPhi",myalivpart->Phi());
	  FillHistogram("trackTriggerTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nTrig+=1;
	}
	if(IsSelectedAssociated(myalivpart)){
	  FillHistogram("trackAssociatedPt",myalivpart->Pt());
	  FillHistogram("trackAssociatedPhi",myalivpart->Phi());
	  FillHistogram("trackAssociatedTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nAss+=1;
	}      
	if(!(IsSelectedAssociated(myalivpart)||IsSelectedTrigger(myalivpart))){
	  //if neither, delete the thing.
	  delete myalivpart->Particle();delete myalivpart;
	}
      }
      for (Int_t ipart=0; ipart<Njetpart_a; ipart++){
	Double_t PT;//pT of the particle, if there is a distribution, use it.
	if(triggerpT&&associatedpT) 	PT = associatedpT->GetRandom();
	else				PT = 3.0;	
	Double_t pxpart, pypart, pzpart;
	Double_t incldifa = gRandom->Gaus(0.0,awaywidth); // Gaussian width
	Double_t angledifa = TMath::Pi()*(gRandom->Rndm()-0.5);	
	Double_t away_phidif = incldifa*TMath::Cos(angledifa);
	Double_t away_thetadif = incldifa*TMath::Sin(angledifa);
	Double_t phipart = jetphi + TMath::Pi() + away_phidif;
	Double_t thetapart = 2.0*TMath::ATan(TMath::Exp(-1.0*jeteta_away))+ away_thetadif;
	while((phipart<(-0.5*Pii))||(phipart>(1.5*Pii))){
	  if (phipart<(-0.5*Pii)) phipart += TMath::TwoPi();
	  if (phipart>(1.5*Pii)) phipart -= TMath::TwoPi();
	}
	while((thetapart<0)||(thetapart>Pii)){
	  if (thetapart>(Pii)) thetapart-=TMath::TwoPi();//if it is in the range pi->2pi, it is moved to -pi->0
	  if (thetapart<(0.0)&&thetapart>-Pii) thetapart*=-1;//if it is in the range -pi->0, it is moved into 0->pi
	  if (thetapart<(0.0)&&thetapart<-Pii) thetapart+=TMath::TwoPi();	//if in range -2pi->-pi, it is moved into 0->pi
	}
	Double_t etapart = -TMath::Log(TMath::Tan(thetapart/2.0));
	Double_t P = PT*TMath::CosH(etapart);//pT*cosh(eta)
	TParticle * mypart = new TParticle();
	nTotal +=1;
	pxpart = P*TMath::Sin(thetapart)*TMath::Cos(phipart);
	pypart = P*TMath::Sin(thetapart)*TMath::Sin(phipart);
	pzpart = P*TMath::Cos(thetapart);
	mypart->SetMomentum(pxpart,pypart,pzpart,pxpart*pxpart+pypart*pypart+pzpart*pzpart);
	mypart->SetPdgCode(awaycharge*11);//electron if charge is positive, positron if negative.
	awaycharge *= -1; //Flip charge for each particle (neutral or one charge in the jet).
	AliMCParticle* myalivpart= new AliMCParticle(mypart);
	FillHistogram("trackPt",myalivpart->Pt());
	FillHistogram("trackPhi",myalivpart->Phi());
	FillHistogram("trackTheta",myalivpart->Theta());
	if(IsSelectedTrigger(myalivpart)){
	  FillHistogram("trackTriggerPt",myalivpart->Pt());
	  FillHistogram("trackTriggerPhi",myalivpart->Phi());
	  FillHistogram("trackTriggerTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nTrig+=1;
	}
	if(IsSelectedAssociated(myalivpart)){
	  FillHistogram("trackAssociatedPt",myalivpart->Pt());
	  FillHistogram("trackAssociatedPhi",myalivpart->Phi());
	  FillHistogram("trackAssociatedTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nAss+=1;
	}
	if(!(IsSelectedAssociated(myalivpart)||IsSelectedTrigger(myalivpart))){
	  //if neither, delete the thing.
	  delete myalivpart->Particle();delete myalivpart;
	}
      }
      ndijets+=1;
    }
    //generate NJetTriplets tri-jets:
    ntrijets = 0;
    while(ntrijets<NJetTriplets){
      if(v2graph&&!triggerpT) v2 = v2graph->Eval(7.0);//No pT dist, all triggers are 7.0
      if(v2graph&&triggerpT) v2 = v2graph->Eval(0.5*triggerpT->GetXmax()+0.5*triggerpT->GetXmin());
      if(v2 > 0) PhiDist->SetParameter(1,v2);
      else PhiDist->SetParameter(1,0.0);
      if(v3graph&&!triggerpT) v3 = v3graph->Eval(7.0);//No pT dist, all triggers are 7.0
      if(v3graph&&triggerpT)  v3 = v3graph->Eval(0.5*triggerpT->GetXmax()+0.5*triggerpT->GetXmin());
      if(v3 > 0) PhiDist->SetParameter(2,v3);
      else PhiDist->SetParameter(2,0.0);
      // Choosing random jet axis azimuthal angles and near- and away-side
      // pseudorapidities (flat distributions)
      Double_t jetphi =  PhiDist->GetRandom(); //slight correlation between reaction plane and jet direction.
      Double_t jetdeltaphiaw = Trijetawaysplitting*(gRandom->Rndm());//Flat distribution in delta phi with width as given.
      Double_t jeteta_near = 2.0*(gRandom->Rndm()-0.5);
      Double_t jeteta_away = 2.0*(gRandom->Rndm()-0.5);//for now both away jets have the same width and number of particles.
      Short_t charge = 1;
      if(gRandom->Rndm()>0.5) charge *=-1;
      Short_t nearcharge = charge;
      Short_t awaycharge = -1*charge;//first is opposite charge of nearside, second is same.
      Short_t awaycharge2 = charge;

      // Generate trijet.  No transverse momentum effects considered,
      // only direction (azimuthal angle and pseudorapidity)
      for (Int_t ipart=0; ipart<Njetpart_n_3; ipart++){
	Double_t PT;//pT of the particle, if there is a distribution, use it.
	if(triggerpT&&associatedpT){
	  if(ipart == 0)PT = triggerpT->GetRandom();
	  else PT = associatedpT->GetRandom();
	}
	else{
	  if(ipart == 0)PT = 7.0;
	  else PT = 3.0;
	}
	
	Double_t pxpart, pypart, pzpart;
	Double_t incldif = gRandom->Gaus(0.0,nearwidth_3); // Gaussian width
	Double_t angledif = TMath::Pi()*(gRandom->Rndm()-0.5);
	Double_t near_phidif = incldif*TMath::Cos(angledif);// gRandom->Gaus(0.0,nearwidth); // Gaussian width
	Double_t near_thetadif =incldif*TMath::Sin(angledif);// gRandom->Gaus(0.0,nearwidth);
	Double_t phipart = jetphi + near_phidif;
	Double_t thetapart = 2.0*TMath::ATan(TMath::Exp(-1.0*jeteta_near))+ near_thetadif;
	while((phipart<(-0.5*Pii))||(phipart>(1.5*Pii))){
	  if (phipart<(-0.5*Pii)) phipart += TMath::TwoPi();
	  if (phipart>(1.5*Pii)) phipart -= TMath::TwoPi();
	}
	while((thetapart<0)||(thetapart>Pii)){
	  if (thetapart>(Pii)) thetapart-=TMath::TwoPi();//if it is in the range pi->2pi, it is moved to -pi->0
	  if (thetapart<(0.0)&&thetapart>-Pii) thetapart*=-1;//if it is in the range -pi->0, it is moved into 0->pi
	  if (thetapart<(0.0)&&thetapart<-Pii) thetapart+=TMath::TwoPi();	//if in range -2pi->-pi, it is moved into 0->pi
	}
	Double_t etapart = -TMath::Log(TMath::Tan(thetapart/2.0));
	Double_t P = PT*TMath::CosH(etapart);//pT*cosh(eta)	
	TParticle* mypart = new TParticle();
	nTotal +=1;
	pxpart = P*TMath::Sin(thetapart)*TMath::Cos(phipart);
	pypart = P*TMath::Sin(thetapart)*TMath::Sin(phipart);
	pzpart = P*TMath::Cos(thetapart);
	mypart->SetMomentum(pxpart,pypart,pzpart,pxpart*pxpart+pypart*pypart+pzpart*pzpart);
	mypart->SetPdgCode(nearcharge*11);//electron if charge is positive, positron if negative.
	nearcharge *= -1; //Flip charge for each particle (neutral or one charge in the jet).
	AliMCParticle* myalivpart= new AliMCParticle(mypart);
	FillHistogram("trackPt",myalivpart->Pt());
	FillHistogram("trackPhi",myalivpart->Phi());
	FillHistogram("trackTheta",myalivpart->Theta());
	if(IsSelectedTrigger(myalivpart)){
	  FillHistogram("trackTriggerPt",myalivpart->Pt());
	  FillHistogram("trackTriggerPhi",myalivpart->Phi());
	  FillHistogram("trackTriggerTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nTrig+=1;
	}
	if(IsSelectedAssociated(myalivpart)){
	  FillHistogram("trackAssociatedPt",myalivpart->Pt());
	  FillHistogram("trackAssociatedPhi",myalivpart->Phi());
	  FillHistogram("trackAssociatedTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nAss+=1;
	}      
	if(!(IsSelectedAssociated(myalivpart)||IsSelectedTrigger(myalivpart))){
	  //if neither, delete the thing.
	  delete myalivpart->Particle();delete myalivpart;
	}
      }
      //First away side jet.
      for (Int_t ipart=0; ipart<Njetpart_a_3; ipart++){
	Double_t PT;//pT of the particle, if there is a distribution, use it.
	if(triggerpT&&associatedpT)PT = associatedpT->GetRandom();
	else PT = 3.0;
	Double_t pxpart, pypart, pzpart;	
	Double_t incldif = gRandom->Gaus(0.0,awaywidth_3); // Gaussian width
	Double_t angledif = TMath::Pi()*(gRandom->Rndm()-0.5);
	Double_t away_phidif =incldif*TMath::Cos(angledif);// gRandom->Gaus(0.0,awaywidth_3);
	Double_t away_thetadif = incldif*TMath::Sin(angledif);// gRandom->Gaus(0.0,awaywidth_3);
	Double_t phipart = jetphi + TMath::Pi()+jetdeltaphiaw + away_phidif;
	Double_t thetapart = 2.0*TMath::ATan(TMath::Exp(-1.0*jeteta_away))+ away_thetadif;
	while((phipart<(-0.5*Pii))||(phipart>(1.5*Pii))){
	  if (phipart<(-0.5*Pii)) phipart += TMath::TwoPi();
	  if (phipart>(1.5*Pii)) phipart -= TMath::TwoPi();
	}
	while((thetapart<0)||(thetapart>Pii)){
	  if (thetapart>(Pii)) thetapart-=TMath::TwoPi();//if it is in the range pi->2pi, it is moved to -pi->0
	  if (thetapart<(0.0)&&thetapart>-Pii) thetapart*=-1;//if it is in the range -pi->0, it is moved into 0->pi
	  if (thetapart<(0.0)&&thetapart<-Pii) thetapart+=TMath::TwoPi();	//if in range -2pi->-pi, it is moved into 0->pi
	}
	Double_t etapart = -TMath::Log(TMath::Tan(thetapart/2.0));
	Double_t P = PT*TMath::CosH(etapart);//pT*cosh(eta)
	TParticle * mypart = new TParticle();
	nTotal +=1;
	pxpart = P*TMath::Sin(thetapart)*TMath::Cos(phipart);
	pypart = P*TMath::Sin(thetapart)*TMath::Sin(phipart);
	pzpart = P*TMath::Cos(thetapart);
	mypart->SetMomentum(pxpart,pypart,pzpart,pxpart*pxpart+pypart*pypart+pzpart*pzpart);
	mypart->SetPdgCode(awaycharge*11);//electron if charge is positive, positron if negative.
	awaycharge *= -1; //Flip charge for each particle (neutral or one charge in the jet).
	AliMCParticle* myalivpart= new AliMCParticle(mypart);
	FillHistogram("trackPt",myalivpart->Pt());
	FillHistogram("trackPhi",myalivpart->Phi());
	FillHistogram("trackTheta",myalivpart->Theta());
	if(IsSelectedTrigger(myalivpart)){
	  FillHistogram("trackTriggerPt",myalivpart->Pt());
	  FillHistogram("trackTriggerPhi",myalivpart->Phi());
	  FillHistogram("trackTriggerTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nTrig+=1;
	}
	if(IsSelectedAssociated(myalivpart)){
	  FillHistogram("trackAssociatedPt",myalivpart->Pt());
	  FillHistogram("trackAssociatedPhi",myalivpart->Phi());
	  FillHistogram("trackAssociatedTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nAss+=1;
	}      
	if(!(IsSelectedAssociated(myalivpart)||IsSelectedTrigger(myalivpart))){
	  //if neither, delete the thing.
	  delete myalivpart->Particle();delete myalivpart;
	}
      }
      //second away side jet.
      for (Int_t ipart=0; ipart<Njetpart_a_3; ipart++){
	Double_t PT;//pT of the particle, if there is a distribution, use it.
	if(triggerpT&&associatedpT) PT = associatedpT->GetRandom();
	else PT = 3.0;	
	Double_t pxpart, pypart, pzpart;	
	Double_t incldif = gRandom->Gaus(0.0,awaywidth_3); // Gaussian width
	Double_t angledif = TMath::Pi()*(gRandom->Rndm()-0.5);
	Double_t away_phidif =incldif*TMath::Cos(angledif);// gRandom->Gaus(0.0,awaywidth_3);
	Double_t away_thetadif = incldif*TMath::Sin(angledif);// gRandom->Gaus(0.0,awaywidth_3);
	Double_t phipart = jetphi + TMath::Pi()-jetdeltaphiaw + away_phidif;
	Double_t thetapart = 2.0*TMath::ATan(TMath::Exp(-1.0*jeteta_away))+ away_thetadif;
	while((phipart<(-0.5*Pii))||(phipart>(1.5*Pii))){
	  if (phipart<(-0.5*Pii)) phipart += TMath::TwoPi();
	  if (phipart>(1.5*Pii)) phipart -= TMath::TwoPi();
	}
	while((thetapart<0)||(thetapart>Pii)){
	  if (thetapart>(Pii)) thetapart-=TMath::TwoPi();//if it is in the range pi->2pi, it is moved to -pi->0
	  if (thetapart<(0.0)&&thetapart>-Pii) thetapart*=-1;//if it is in the range -pi->0, it is moved into 0->pi
	  if (thetapart<(0.0)&&thetapart<-Pii) thetapart+=TMath::TwoPi();	//if in range -2pi->-pi, it is moved into 0->pi
	}
	Double_t etapart = -TMath::Log(TMath::Tan(thetapart/2.0));
	Double_t P = PT*TMath::CosH(etapart);//pT*cosh(eta)
	TParticle * mypart = new TParticle();
	nTotal +=1;
	pxpart = P*TMath::Sin(thetapart)*TMath::Cos(phipart);
	pypart = P*TMath::Sin(thetapart)*TMath::Sin(phipart);
	pzpart = P*TMath::Cos(thetapart);
	mypart->SetMomentum(pxpart,pypart,pzpart,pxpart*pxpart+pypart*pypart+pzpart*pzpart);
	mypart->SetPdgCode(awaycharge2*11);//electron if charge is positive, positron if negative.
	awaycharge2 *= -1; //Flip charge for each particle (neutral or one charge in the jet).
	AliMCParticle* myalivpart= new AliMCParticle(mypart);
	FillHistogram("trackPt",myalivpart->Pt());
	FillHistogram("trackPhi",myalivpart->Phi());
	FillHistogram("trackTheta",myalivpart->Theta());
	if(IsSelectedTrigger(myalivpart)){
	  FillHistogram("trackTriggerPt",myalivpart->Pt());
	  FillHistogram("trackTriggerPhi",myalivpart->Phi());
	  FillHistogram("trackTriggerTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nTrig+=1;
	}
	if(IsSelectedAssociated(myalivpart)){
	  FillHistogram("trackAssociatedPt",myalivpart->Pt());
	  FillHistogram("trackAssociatedPhi",myalivpart->Phi());
	  FillHistogram("trackAssociatedTheta",myalivpart->Theta());
	  allrelevantParticles.Add(myalivpart);
	  nAss+=1;
	}      
	if(!(IsSelectedAssociated(myalivpart)||IsSelectedTrigger(myalivpart))){
	  //if neither, delete the thing.
	  delete myalivpart->Particle();delete myalivpart;
	}
      }
      ntrijets+=1;
    }
    //Generate flat and flow background particles.
    nFlatP = 0;
    nFlowP = 0;
    while(nFlatP<NFlat||nFlowP<NFLOWparticles){
      Bool_t isflow = ((NFlat-nFlatP)<1)&&(NFLOWparticles>0);//is a flow particle, if all bg particles are made.
      Double_t PT;//pT of the particle, if there is a distribution, use it.
      if(triggerpT&&associatedpT){
	PT = associatedpT->GetRandom();
	if(nFlatP==0&&gRandom->Rndm()>0.5&&!isflow)PT = triggerpT->GetRandom();
	if(nFlowP==0&&isflow)PT = triggerpT->GetRandom();
      }
      else{
	PT = 3.0;
	if(nFlatP==0&&gRandom->Rndm()>0.5&&!isflow)PT = 7.0;
	if(nFlowP==0&&isflow)PT = 7.0;
      }      
      if(isflow){
	if(v2graph) v2 = v2graph->Eval(PT);
	PhiDist->SetParameter(1,v2);
	if(v3graph)  v3 = v3graph->Eval(PT);
	PhiDist->SetParameter(2,v3);
      }
      Double_t pxpart,pypart,pzpart;
      Double_t pphi;
      if(isflow)  pphi = PhiDist->GetRandom();
      if(!isflow) pphi = TMath::TwoPi()*(gRandom->Rndm()-0.25);
      Double_t peta = 2.0*(gRandom->Rndm()-0.5);
      Double_t ptheta = 2.0*TMath::ATan(TMath::Exp(-1.0*peta));
      Short_t  charge = 1;
      if(gRandom->Rndm()>0.5)charge*=-1;
      while((pphi<(-0.5*Pii))||(pphi>(1.5*Pii))){
	if (pphi<(-0.5*Pii)) pphi += TMath::TwoPi();
	if (pphi>(1.5*Pii)) pphi -= TMath::TwoPi();
      }
      if (ptheta>(Pii)) ptheta-=TMath::TwoPi();//if it is in the range pi->2pi, it is moved to -pi->0
      if (ptheta<(0.0)&&ptheta>-Pii) ptheta*=-1;//if it is in the range -pi->0, it is moved into 0->pi
      if (ptheta<(0.0)&&ptheta<-Pii) ptheta+=TMath::TwoPi();	//if in range -2pi->-pi, it is moved into 0->pi
      Double_t P = PT*TMath::CosH(peta);//pT*cosh(eta)	
      TParticle* mypart = new TParticle();
      nTotal +=1;
      pxpart = P*TMath::Sin(ptheta)*TMath::Cos(pphi);
      pypart = P*TMath::Sin(ptheta)*TMath::Sin(pphi);
      pzpart = P*TMath::Cos(ptheta);
      mypart->SetMomentum(pxpart,pypart,pzpart,pxpart*pxpart+pypart*pypart+pzpart*pzpart);
      mypart->SetPdgCode(charge*11);
      AliMCParticle* myalivpart= new AliMCParticle(mypart);

      FillHistogram("trackPt",myalivpart->Pt());
      FillHistogram("trackPhi",myalivpart->Phi());
      FillHistogram("trackTheta",myalivpart->Theta());
      if(IsSelectedTrigger(myalivpart)){
	FillHistogram("trackTriggerPt",myalivpart->Pt());
	FillHistogram("trackTriggerPhi",myalivpart->Phi());
	FillHistogram("trackTriggerTheta",myalivpart->Theta());
	allrelevantParticles.Add(myalivpart);
	nTrig+=1;
      }
      if(IsSelectedAssociated(myalivpart)){
	FillHistogram("trackAssociatedPt",myalivpart->Pt());
	FillHistogram("trackAssociatedPhi",myalivpart->Phi());
	FillHistogram("trackAssociatedTheta",myalivpart->Theta());
	allrelevantParticles.Add(myalivpart);
	nAss+=1;
      }      
      if(!(IsSelectedAssociated(myalivpart)||IsSelectedTrigger(myalivpart))){
	//if neither, delete the thing.
	delete myalivpart->Particle();delete myalivpart;
      }
      if(!isflow)nFlatP+=1;
      if(isflow)nFlowP+=1;
    }
    //After this the generated event is analyzed
    FillHistogram("Ntriggers",nTrig);
    FillHistogram("NAssociated",nAss);
    FillHistogram("trackCount",nTotal);

    dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p>*>(fCorrelator)->SetEventVzM(0.0,50);
    fCorrelator->Execute(NULL,&allrelevantParticles);
    PostData(1,fOutput);
    AliMCParticle * deletepart;
    for(int i=0;i<allrelevantParticles.GetEntriesFast();i++){
      deletepart = (AliMCParticle*)allrelevantParticles.At(i);
      if(!IsSelectedAssociated(deletepart)){delete deletepart->Particle();delete deletepart;} //delete triggers and non involved particles.
    }
    delete PhiDist;
  }
}

void AliAnalysisTaskCorrelation3p::Askforgensettings()
{//Ask the user for generator settings, write them in generatorsetting.root file.
  TFile* settings = TFile::Open("generatorsetting.root","RECREATE");
  settings->cd();
  TTree *tree = new TTree("Parameters","parameters for particle generation");
  
//   int nints = 5;//number of setters with integer value.
//   int ndoubles = 2;//number of setters with double value.
//   TArrayD doublearray = TArrayD(nints);
//   TArrayI intarray = TArrayI(ndoubles);
  

  Int_t nevents=1000; //intarray:0
  tree->Branch("nevents",&nevents,"nevents/I");
  Int_t NJetPairs    =  1; //intarray:3 number of jet pairs per event
  tree->Branch("NJetPairs",&NJetPairs,"NJetPairs/I");
  Int_t Njetpart_n   =  5; //number of particles in near-side jet
  tree->Branch("Njetpart_n",&Njetpart_n,"Njetpart_n/I");
  Double_t nearwidth = 0.1;//gaussian width of the near side jet in radians.
  tree->Branch("nearwidth",&nearwidth,"nearwidth/D");  
  Int_t Njetpart_a   =  3; // number of particles in away-side jet
  tree->Branch("Njetpart_a",&Njetpart_a,"Njetpart_a/I");
  Double_t awaywidth = 0.2;//gaussian width of the away side jet in radians.  
  tree->Branch("awaywidth",&awaywidth,"awaywidth/D");
  Int_t NJetTriplets =  0; // number of jet triplets per event
  tree->Branch("NJetTriplets",&NJetTriplets,"NJetTriplets/I");
  Int_t Njetpart_n_3  =  5; //number of particles in near-side jet
  tree->Branch("Njetpart_n_3",&Njetpart_n_3,"Njetpart_n_3/I");
  Double_t nearwidth_3 = 0.1;//gaussian width of the near side jet in radians.
  tree->Branch("nearwidth_3",&nearwidth_3,"nearwidth_3/D");
  Int_t Njetpart_a_3   =  3; // number of particles in away-side jet
  tree->Branch("Njetpart_a_3",&Njetpart_a_3,"Njetpart_a_3/I");
  Double_t awaywidth_3 = 0.2;//gaussian width of the away side jet in radians.
  tree->Branch("awaywidth_3",&awaywidth_3,"awaywidth_3/D");
  Double_t Trijetawaysplitting = 0.5;//splitting of the Phi of the away side dijet in radians.
  tree->Branch("Trijetawaysplitting",&Trijetawaysplitting,"Trijetawaysplitting/D");
  Double_t NFLOWparticles = 0;//number of particles from flow.
  tree->Branch("NFLOWparticles",&NFLOWparticles,"NFLOWparticles/D");
  Double_t v2        =0.0; //doublearray:0 v2
  tree->Branch("v2",&v2,"v2/D");
  Double_t v3        =0.0; //doublearray:1 v3
  tree->Branch("v3",&v3,"v3/D");
  Int_t NFlat        =0;//Number of background particles (flat angular distribution.
  tree->Branch("NFlat",&NFlat,"NFlat/I");

  
  AliWarning("This parser reads in parameters for particle generation.");
  AliWarning("The defaults is one Phi back to back jet pair with the following parameters:");
  AliWarning(Form("   # events = %i", nevents));
  AliWarning(Form("   # particles in near side jet = %i",Njetpart_n)) ;
  AliWarning(Form("   width of the near side jet   = %f [radians]",awaywidth));
  AliWarning(Form("   # particles in away side jet = %i",Njetpart_a));
  AliWarning(Form("   width of the away side jet   = %f [radians]",awaywidth));
  AliWarning("In addition you can choose to add three-jet objects and v2 and v3.");
  AliWarning("The toy MC can be initialized with a realistic pT distribution,");
  AliWarning("parametrized from real Data. The production must have");
  AliWarning("an AnalysisResults.root file located in ../../LHC1xx relative to this folder.");
  AliWarning("Standard is no pT distribution.");
  AliWarning("To use standard values, please type 0, to use other values type 1:");
  bool worked = false;
  string input;
  while(!worked){
    int defaultv;
    getline(cin, input);
    stringstream ss(input);
    ss>>defaultv;
    if(ss.fail()) {AliWarning("You need to enter either 1 or 0:"); continue;}
    if(defaultv==0){//just write the defaults to file.
      tree->Fill();
      tree->Write();
      settings->Close();
      return;}
    if(defaultv==1){worked=true;}
    else AliWarning("You need to enter either 1 or 0:");
  }
  worked = false;
  while(!worked){
    AliWarning("Please enter the number of events per (data) event:");
    int nev;
    getline(cin, input);
    stringstream ss(input);
    ss>>nev;
    if(ss.fail()) { AliWarning("You need to enter an integer:"); continue;}
    else {nevents=nev;worked = true;}
  }
  worked = false;
  if(nevents==0){
    AliWarning("You entered 0 events, finishing");
    tree->Fill();
    tree->Write();
    settings->Close();
    return;}
  while(!worked){
    AliWarning("Please provide the real data period you want to extrapolate the pT distribution from: LHC1");
    getline(cin, input);
    const char* period = input.c_str();
    AliWarning(Form("You choose LHC %s.",period));
    AliWarning(Form("Will now try to open the file ../../LHC1%s/results.root",period));
    TFile * AResults = TFile::Open(Form("../../LHC1%s/results.root",period),"READ");
    if(!AResults){AliWarning("File not found. Defaulting to no pT distribution. Please restart and provide the file if the parametrization is nessecary.");worked = true;}
    else{
      if(AResults->IsZombie()){AliWarning("File not found. Defaulting to no pT distribution. Please restart and provide the file if the parametrization is nessecary.");worked =true;}
      else{
	TString cent = TString("");
	//In the PbPb case, one needs to provide the correct centrality bins:
	if(TString(period).CompareTo("0h")==0||TString(period).CompareTo("1h")==0){
	  TList * list = AResults->GetListOfKeys();
	  TObjArray * CentBins = new TObjArray(20);
	  CentBins->SetOwner(true);
	  AliWarning("Please provide the centrality window from the following alternatives:");
	  int j = 1;
	  for(int i = 0;i<list->GetEntries();i++){
	    TObjString * name = new TObjString(list->At(i)->GetName());
	    if(name->String().Contains("BinM")&&!name->String().Contains("Z")){AliWarning(Form("%i: %s",j,name->String().Data())); CentBins->AddAt(name,j);j++;}
	    else delete name;
	  }
	  bool worked2 = false;
	  while(!worked2){
	    AliWarning("Please type in the number of the preferred multiplicity bin: ");
	    getline(cin, input);
	    stringstream ss(input);
	    ss>>j;	  
	    if(ss.fail()){AliWarning("You need to enter an integer from the list above.");continue;}
	    else{
	      TObjString * tmp = dynamic_cast<TObjString*>(CentBins->At(j));
	      if(!tmp){AliWarning("You need to enter an integer from the list above:");continue;}
	      else {cent.Append(tmp->GetString().Data());worked2 = true;}
	    }
	  }
	  delete CentBins;
	}
	TH1D * triggerpT;
	TH1D * associatedpT;
	if(cent.CompareTo("") == 0){
	  TDirectory * statsdir =  AResults->GetDirectory("bin_stats");
	  if(!statsdir){AliWarning("Could not locate directory bin_stats in results.root. Defaulting to no pT distribution. Please provide the correct file if the parametrization is nessecary."); worked = true; continue;}
	  triggerpT = dynamic_cast<TH1D*>(statsdir->Get("trigger_pT"));
	  associatedpT = dynamic_cast<TH1D*>(statsdir->Get("associated_pT"));
	}
	else{
	  TDirectory * statsdir = AResults->GetDirectory(Form("%s/bin_stats",cent.Data()));
	  if(!statsdir){AliWarning(Form("Could not locate directory %s/bin_stats in results.root. Defaulting to no pT distribution. Please provide the correct file if the parametrization is nessecary.",cent.Data())); worked = true; continue;}
	  triggerpT = dynamic_cast<TH1D*>(statsdir->Get("trigger_pT"));
	  associatedpT = dynamic_cast<TH1D*>(statsdir->Get("associated_pT"));	  
	}
	if(!triggerpT||!associatedpT){AliWarning(Form("Could not locate histograms in directory %s/bin_stats in results.root. Defaulting to no pT distribution. Please provide the correct file if the parametrization is nessecary.",cent.Data())); worked = true; continue;}
	TF1* funct = pTdistribution(triggerpT,"triggerpT");
	TF1* funca = pTdistribution(associatedpT,"associatedpT");
	settings->cd();
	funct->Write();
	funca->Write();
	AResults->Close();
	delete funct;delete funca;
	worked = true;
      }
    }
  }
  worked = false;
    while(!worked){
    AliWarning("Please enter the number of dijets per event:");
    int njets;
    getline(cin, input);
    stringstream ss(input);
    ss>>njets;
    if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
    else {NJetPairs=njets;worked = true;}
  }  
  worked = false;
  if(NJetPairs !=0){//or there is no point
    while(!worked){
      AliWarning("Please enter the number of particles in the nearside in dijets:");
      int npart;
      getline(cin, input);
      stringstream ss(input);
      ss>>npart;
      if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
      if(Njetpart_n<1.0){AliWarning("The number you enter must be an integer bigger then 1"); continue;}
      else {Njetpart_n=npart;worked = true;}
    }  
    worked = false;
    while(!worked){
      AliWarning("Please enter the width of the near side jet of the dijet in radians (double):");
      double width;
      getline(cin, input);
      stringstream ss(input);
      ss>>width;
      if(ss.fail()) {AliWarning("You need to enter a double:"); continue;}
      else {nearwidth=width;worked = true;}
    }  
    worked = false;
    while(!worked){
      AliWarning("Please enter the number of particles in the awayside jet of the dijet :");
      int npart;
      getline(cin, input);
      stringstream ss(input);
      ss>>npart;
      if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
      else {Njetpart_a=npart;worked = true;}
    }  
    worked = false;
      while(!worked){
      AliWarning("Please enter the width of the away side jet in radians (double):");
      double width;
      getline(cin, input);
      stringstream ss(input);
      ss>>width;
      if(ss.fail()) {AliWarning("You need to enter a double:"); continue;}
      else {awaywidth=width;worked = true;}
    }  
    worked = false;
  }
  while(!worked){
    AliWarning("Please enter the number of trijets per event:");
    int njets;
    getline(cin, input);
    stringstream ss(input);
    ss>>njets;
    if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
    else {NJetTriplets=njets;worked = true;}
  }  
  worked = false;
  if(NJetTriplets!=0){//or there is no point
    while(!worked){
      AliWarning("Please enter the number of particles in the nearside in jet triplet:");
      int npart;
      getline(cin, input);
      stringstream ss(input);
      ss>>npart;
      if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
      if(Njetpart_n<1.0){AliWarning("The number you enter must be an integer bigger then 1"); continue;}
      else {Njetpart_n=npart;worked = true;}
    }      
    worked = false;
    while(!worked){
      AliWarning("Please enter the width of the near side jet of the jet triplet in radians (double):");
      double width;
      getline(cin, input);
      stringstream ss(input);
      ss>>width;
      if(ss.fail()) {AliWarning( "You need to enter a double:"); continue;}
      else {nearwidth_3=width;worked = true;}
    }  
    worked = false;
    while(!worked){
      AliWarning("Please enter the number of particles in the awayside jets of the jet triplet :");
      int npart;
      getline(cin, input);
      stringstream ss(input);
      ss>>npart;
      if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
      else {Njetpart_a_3=npart;worked = true;}
    }  
    worked = false;
      while(!worked){
      AliWarning("Please enter the width of the away side jets in radians (double):");
      double width;
      getline(cin, input);
      stringstream ss(input);
      ss>>width;
      if(ss.fail()) {AliWarning("You need to enter a double:"); continue;}
      else {awaywidth_3=width;worked = true;}
    }  
    worked = false;
    while(!worked){
      AliWarning("Please enter the maximum value of the away side dijet splitting in radians (double):");
      double width;
      getline(cin, input);
      stringstream ss(input);
      ss>>width;
      if(ss.fail()) {AliWarning("You need to enter a double:"); continue;}
      else {Trijetawaysplitting=width;worked = true;}
    }  
    worked = false;
  }  
  
  while(!worked){
    AliWarning("Please enter the number of particles from flow per event:");
    int nflow;
    getline(cin, input);
    stringstream ss(input);
    ss>>nflow;
    if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
    else {NFLOWparticles=nflow;worked = true;}
  }  
  worked = false;
  if(NFLOWparticles!=0){
    while(!worked){
      AliWarning("If you want to use a realistic pT distributed flow from ../../flowpt.root, please enter 1. If you want to manually insert v2 and v3, please enter 0:");
      int yesno;
      getline(cin,input);
      stringstream ss(input);
      ss>>yesno;
      if(ss.fail()||(yesno!=0&&yesno!=1)){AliWarning("You need to enter either 1 or 0:");continue;}
      else{
	if(yesno == 0){
	  bool inworked =false;
	  while(!inworked){
	    AliWarning("Please enter v2 (double) to add flow:");
	    double v2l;
	    getline(cin, input);
	    stringstream ss1(input);
	    ss1>>v2l;
	    if(ss1.fail()) {AliWarning("You need to enter a double:"); continue;}
	    else {v2=v2l;inworked = true;}
	    }  
	  inworked = false;
	  while(!inworked){
	    AliWarning("Please enter v3 (double) to add flow:");
	    double v3l;
	    getline(cin, input);
	    stringstream ss1(input);
	    ss1>>v3l;
	    if(ss1.fail()) {AliWarning("You need to enter a double:"); continue;}
	    else {v3=v3l;inworked = true;}
	    }  
	  worked = true;
	}
	if(yesno == 1){
	  AliWarning("Trying to open ../../flowpt.root");
	  TFile * flowpt = TFile::Open("../../flowpt.root","READ");
	  if(!flowpt){AliWarning("File not found. Defaulting to no flow distribution. Please restart and provide the file if the parametrization is nessecary.");}
	  else{
	    if(flowpt->IsZombie()){AliWarning("File not found. Defaulting to no flow distribution. Please restart and provide the file if the parametrization is nessecary.");}
	    else{
	      bool worked3 = false;
	      AliWarning("Please choose a centrality range for the flow:");
	      AliWarning("1   -    0%- 5%");
	      AliWarning("2   -    5%-10%");
	      AliWarning("3   -   10%-20%");
	      AliWarning("4   -   20%-30%");
	      AliWarning("5   -   30%-40%");
	      AliWarning("6   -   40%-50%");
	      AliWarning("Your choice:");
	      int choice;
	      while(!worked3){
		getline(cin, input);
		stringstream ss1(input);
		ss1>>choice;
		if(ss1.fail()) {AliWarning("You need to enter an integer of the ones above:"); continue;}
		else {worked3 = true;}
	      }
	      TGraphAsymmErrors * v2graph;
	      TGraphAsymmErrors * v3graph;
	      if(choice == 1){
		v2graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v20t5prc"));
		v3graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v30t5prc"));}
	      if(choice == 2){
		v2graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v25t10prc"));
		v3graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v35t10prc"));}
	      if(choice == 3){
		v2graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v210t20prc"));
		v3graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v310t20prc"));}
	      if(choice == 4){
		v2graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v220t30prc"));
		v3graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v320t30prc"));}
	      if(choice == 5){
		v2graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v230t40prc"));
		v3graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v330t40prc"));}	     
	      if(choice == 6){
		v2graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v240t50prc"));
		v3graph = dynamic_cast<TGraphAsymmErrors*>(flowpt->Get("v340t50prc"));}	
	      if(!v2graph||!v3graph){AliWarning("The file does not contain the expected distributions. Please provide the correct file. Defaulting to no flow.");worked = true; continue;}
	      settings->cd();
	      v2graph->Write("v2ptdist");
	      v3graph->Write("v3ptdist");
	    }
	  }
	  worked = true;
	  flowpt->Close();
	}
      }
    }
  }
  
  worked=kFALSE;
  while(!worked){
   AliWarning("Please enter the number of background associated particles with a flat distribution:");
    int nflat;
    getline(cin, input);
    stringstream ss(input);
    ss>>nflat;
    if(ss.fail()) {AliWarning("You need to enter an integer:"); continue;}
    else {NFlat=nflat;worked = true;}
  }  
  worked = false;  
  worked = false;
  settings->cd();
  tree->Fill();
  tree->Write();
  settings->Close();
//   delete tree;
  return;
 
  
}
Double_t pTdist(Double_t *x, Double_t *par)
{
//Mixture of power-law and exponential function. For now only power-law.
  Double_t B = par[0];//B is the scale of the power law region.
  Double_t k = par[1];//k is the power law exponent.  
//   Double_t A = par[2];//A is the scale of the exponential.
//   Double_t b = par[3];//b is the factor of the exponent.
  if(x[0]<1.0e-4)return 0.0;
  else return B*TMath::Power(x[0],k);
}

TF1* AliAnalysisTaskCorrelation3p::pTdistribution(TH1D* hist, const char* name)
{
  Double_t minpT = hist->GetXaxis()->GetBinCenter(1);
  Double_t maxpT=hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
  bool minis = false;
  for(int i =1;i<=hist->GetNbinsX();i++)
  {
    if(!minis && hist->GetBinContent(i)>0.0){minpT = hist->GetXaxis()->GetBinLowEdge(i);minis = true;}
    if( minis && hist->GetBinContent(i)>0.0){maxpT = hist->GetXaxis()->GetBinUpEdge(i);}
  }
  TF1 * pTfunc = new TF1(name, pTdist,minpT,maxpT,2);
  pTfunc->SetParNames("B","k");
  pTfunc->SetParLimits(0,1.0,1.0e10);
  pTfunc->SetParLimits(1,-1.0e5,-1.0e-10);
  pTfunc->SetParameters(1.0,-1.2);
  hist->Fit(pTfunc,"MSQ","",minpT,maxpT);
  return pTfunc;
}











