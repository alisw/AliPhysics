#ifndef __CINT__
#include "AliAnalysisTaskCorrelation3p.h"
#endif
AliAnalysisTaskCorrelation3p* AddTaskThreePartTracksTree ()
{
  //Add a task AliAnalysisTaskCorrelation3p to the analysis train in charged track analysis, for PbPb data 
  //Defaults to 10h data with MB trigger.
  //Author: Paul Baetzing
  /* $Id$ */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr){
    ::Error("AddTaskThreePartTracks", "No analysis manager to connect to");
    return NULL;
  }
    if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskThreePartTracks", "This task requires an input event handler");
    return NULL;
  }
  gROOT->LoadMacro("TaskConfig.C");
  const char* fname = Form("%s_%1.0f_%1.0f",name,MinTriggerPt,MaxTriggerPt,MinAssociatedPt,MaxAssociatedPt);
  const char* tname = Form("%s_%1.0f_%1.0f_%1.0f_%1.0f",name,MinTriggerPt,MaxTriggerPt,MinAssociatedPt,MaxAssociatedPt);
  AliAnalysisTaskCorrelation3p* task = new AliAnalysisTaskCorrelation3p(Form("%sTask", tname), "");

  task->SetTrigger(AliAnalysisTaskCorrelation3p::tracks);
  task->SetMinTriggerPt(MinTriggerPt);
  task->SetMaxTriggerPt(MaxTriggerPt);
  task->SetMinAssociatedPt(MinAssociatedPt);
  task->SetMaxAssociatedPt(MaxAssociatedPt);
  task->SetAcceptanceCut(Acceptancecut);
  task->SetTrackCut(cutmask);
  task->SetBinVer(binver);
  task->SetMaxTracksPerEvent(maxntracksmix);
  task->SetMoreOutputs(MoreOutput);
  task->SetLeading(Leading);

  if(QAonly)  task->SetQA();
  if(QAtask)  task->SetQAtask(true);
  task->SetDstTree(true);
  if(TString(file).CompareTo("")!=0&&grid)   task->SetWeights(Form("alien:///alice/cern.ch/user/p/pbatzing/efficiencies/%s",file));
  if(TString(file).CompareTo("")!=0&&!grid)  task->SetWeights(Form("%s",file));

  //Mixing scheme:
  Double_t *Mbin = new Double_t[NMBins+1];
  Mbin[0] = Mbin0;
  Mbin[1] = Mbin1;
  if(NMBins>1) Mbin[2] = Mbin2;
  if(NMBins>2) Mbin[3] = Mbin3;
  if(NMBins>3) Mbin[4] = Mbin4;
  if(NMBins>4) Mbin[5] = Mbin5;
  if(NMBins>5) Mbin[6] = Mbin6;
  if(NMBins>6) Mbin[7] = Mbin7;
  TArrayD tMbin(NMBins+1, Mbin);
  //z vertex binning.
  Double_t *Zbin = new Double_t[NZBins+1];
  Zbin[0] = Zbin0;
  Zbin[1] = Zbin1;
  if(NZBins>1) Zbin[2] = Zbin2;
  if(NZBins>2) Zbin[3] = Zbin3;
  if(NZBins>3) Zbin[4] = Zbin4;
  if(NZBins>4) Zbin[5] = Zbin5;
  if(NZBins>5) Zbin[6] = Zbin6;
  if(NZBins>6) Zbin[7] = Zbin7;
  if(NZBins>7) Zbin[8] = Zbin8;
  if(NZBins>8) Zbin[9] = Zbin9;
  if(NZBins>9) Zbin[10] = Zbin10;
  if(NZBins>10) Zbin[11] = Zbin11;
  if(NZBins>11) Zbin[12] = Zbin12;
  if(NZBins>12) Zbin[13] = Zbin13;
  if(NZBins>13) Zbin[14] = Zbin14;
  if(NZBins>14) Zbin[15] = Zbin15;
  if(NZBins>15) Zbin[16] = Zbin16;
  if(NZBins>16) Zbin[17] = Zbin17;
  if(NZBins>17) Zbin[18] = Zbin18;
  if(NZBins>18) Zbin[19] = Zbin19;
  TArrayD tZbin(NZBins+1, Zbin);  
  task->SetMixingScheme(MaxNEventsMix,MinNTracksMix,tMbin,tZbin);
  
  if( TString(period).Contains("10b") )P10b,P10c,P10d,P10e,P10h,P11a,P11h,
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10b);
  if( TString(period).Contains("10c") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10c);
  if( TString(period).Contains("10d") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10d);
  if( TString(period).Contains("10e") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10e);
  if( TString(period).Contains("11a") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P11a);  
  
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", tname));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), fname));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  return task;
};
