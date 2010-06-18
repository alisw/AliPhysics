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

  AliAnalysisTaskHLTCalo(const char *name);
  virtual ~AliAnalysisTaskHLTCalo() {}
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //virtual Bool_t Notify();
  virtual void NotifyRun();

  virtual void DoSpecificStuff(AliESDEvent * evESD, AliESDEvent * evHLTESD) = 0;
  virtual void CreateSpecificStuff(TList  * fOutputList) = 0;
  virtual Bool_t IsThisDetector(AliESDCaloCluster * cluster) = 0;
  virtual Int_t GetClusters(AliESDEvent * event, TRefArray * clusters) = 0;

private:
  
  AliESDRun *fESDRun;  //!
  TList *fOutputList;
  
  TH1F *fHistOfflResiduals; //!
  TH1F *fHistOnlResiduals; //!

  TH1F *fHistOfflDz; //!
  TH1F *fHistOnlDz; //!

  TH1F *fHistOfflDxy; //!
  TH1F *fHistOnlDxy; //!
  
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

  ClassDef(AliAnalysisTaskHLTCalo, 0);
};

#endif
