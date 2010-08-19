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

  virtual void DoSpecificStuff(AliESDEvent * evESD, AliESDEvent * evHLTESD) {return;}
  virtual void CreateSpecificStuff(TList  * fOutputList) {return;}
  virtual Bool_t IsThisDetector(AliESDCaloCluster * cluster) { return cluster->IsPHOS(); }
  virtual Int_t GetClusters(AliESDEvent * event, TRefArray * clusters) { return event->GetPHOSClusters(clusters); }

private:
  
  AliESDRun *fESDRun;  //!
  TList *fOutputList;
  
  TH1F *fHistOfflResiduals; 
  TH1F *fHistOnlResiduals; 
  TH1F *fHistOfflDz; 
  TH1F *fHistOnlDz; 
  TH1F *fHistOfflDxy; 
  TH1F *fHistOnlDxy; 


  TH1F *fHistOfflResidualsPos; 
  TH1F *fHistOnlResidualsPos; 
  TH1F *fHistOfflDzPos; 
  TH1F *fHistOnlDzPos; 
  TH1F *fHistOfflDxyPos; 
  TH1F *fHistOnlDxyPos; 


  TH1F *fHistOfflResidualsNeg; 
  TH1F *fHistOnlResidualsNeg; 
  TH1F *fHistOfflDzNeg; 
  TH1F *fHistOnlDzNeg; 
  TH1F *fHistOfflDxyNeg; 
  TH1F *fHistOnlDxyNeg; 

  TH2F * fHistNclvsNcl; 
  TH2F * fHistTotEVsTotE;

  
  Int_t fNevt;
  TObjArray *fTrgClsArray;
  
  TObjArray * fGlobalHistoProdArrOff; //!transient 
  TObjArray * fGlobalHistoProdArrOn; //!transient 

  TRefArray * fClustersArray; //!transient

  TString fName; //!transient

  /** copy constructor */
  AliAnalysisTaskHLTCalo(const AliAnalysisTaskHLTCalo&); 
  /** assignment operator */
  AliAnalysisTaskHLTCalo& operator=(const AliAnalysisTaskHLTCalo&); 

  ClassDef(AliAnalysisTaskHLTCalo, 1);
};

#endif
