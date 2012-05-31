// $Id: AliAnalysisTaskHLTCalo.h 40285 2010-04-09 14:04:51Z kkanaki $

#ifndef ALIANALYSISTASKHLTCALO_H
#define ALIANALYSISTASKHLTCALO_H

//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

/** @file AliAnalysisTaskHLTTPC.h
    @author Zhongbao Yin, Kalliopi Kanaki
    @date
    @brief An analysis task to compare the offline and HLT esd trees
*/


// forward declarations
class TH1F;
class TH2F;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliESDRun;
class TObjArray;
class TRefArray;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHLTCalo : public AliAnalysisTaskSE {
 
public:  

  AliAnalysisTaskHLTCalo();
  AliAnalysisTaskHLTCalo(const char *name);
  virtual ~AliAnalysisTaskHLTCalo() {}
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //virtual Bool_t Notify();
  virtual void NotifyRun();

  virtual void DoSpecificStuff(const AliESDEvent * /*evESD*/, const AliESDEvent * /*evHLTESD*/) {return;}
  virtual void CreateSpecificStuff(const TList * /*fOutputList*/) {return;}
  virtual Bool_t IsThisDetector(AliESDCaloCluster * cluster) { return cluster->IsPHOS(); }
  virtual Int_t GetClusters(AliESDEvent * event, TRefArray * clusters) { return event->GetPHOSClusters(clusters); }
 
  //Use only triggered events
  void SetUseHLTTriggerDecision(Bool_t useHLT = kFALSE) {fUseHLTTrigger = useHLT;}

private:
  Bool_t fUseHLTTrigger; // boolean to enable the HLT triggered events
  AliESDRun *fESDRun;  //!Transient, pointer to esds
  TList *fOutputList;  //List of histograms to be stored
  
  TH1F *fHistOfflResiduals; //histogram
  TH1F *fHistOnlResiduals; //histogram
  TH1F *fHistOfflDz; //histogram
  TH1F *fHistOnlDz; //histogram
  TH1F *fHistOfflDxy; //histogram
  TH1F *fHistOnlDxy; //histogram


  TH1F *fHistOfflResidualsPos; //histogram
  TH1F *fHistOnlResidualsPos; //histogram
  TH1F *fHistOfflDzPos; //histogram
  TH1F *fHistOnlDzPos; //histogram
  TH1F *fHistOfflDxyPos; //histogram
  TH1F *fHistOnlDxyPos; //histogram


  TH1F *fHistOfflResidualsNeg; //histogram
  TH1F *fHistOnlResidualsNeg; //histogram
  TH1F *fHistOfflDzNeg; //histogram
  TH1F *fHistOnlDzNeg; //histogram
  TH1F *fHistOfflDxyNeg; //histogram
  TH1F *fHistOnlDxyNeg; //histogram

  TH2F * fHistNclvsNcl; //histogram
  TH2F * fHistTotEVsTotE;//histogram

  
  Int_t fNevt; //Number of events
  TObjArray *fTrgClsArray; //Trigger cluster array
  
  TObjArray * fGlobalHistoProdArrOff; //!transient array of histogram producer classes 
  TObjArray * fGlobalHistoProdArrOn; //!transient array of histogram producer classes 

  TRefArray * fClustersArray; //!transient Array to contain calo clusters

  TString fCaloName; //!transient PHOS or EMCAL

  /** copy constructor */
  AliAnalysisTaskHLTCalo(const AliAnalysisTaskHLTCalo&); 
  /** assignment operator */
  AliAnalysisTaskHLTCalo& operator=(const AliAnalysisTaskHLTCalo&); 

  ClassDef(AliAnalysisTaskHLTCalo, 1);
};

#endif
