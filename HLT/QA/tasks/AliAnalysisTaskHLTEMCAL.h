// $Id: AliAnalysisTaskHLTEMCAL.h 40285 2010-04-09 14:04:51Z kkanaki $

#ifndef ALIANALYSISTASKHLTEMCAL_H
#define ALIANALYSISTASKHLTEMCAL_H

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

#include "AliAnalysisTaskHLTCalo.h"

class AliAnalysisTaskHLTEMCAL : public AliAnalysisTaskHLTCalo {
 
public:  
  AliAnalysisTaskHLTEMCAL(const char *name);
  virtual ~AliAnalysisTaskHLTEMCAL() {}
  
private:
  

  /** copy constructor */
  AliAnalysisTaskHLTEMCAL(const AliAnalysisTaskHLTEMCAL&); 
  /** assignment operator */
  AliAnalysisTaskHLTEMCAL& operator=(const AliAnalysisTaskHLTEMCAL&); 
  
  void CreateSpecificStuff(TList * fOutputList);
  void DoSpecificStuff(AliESDEvent * evESD, AliESDEvent * evHLTESD);
  
  Int_t GetClusters(AliESDEvent * event, TRefArray * clusters);
  Bool_t IsThisDetector(AliESDCaloCluster * cluster);


  ClassDef(AliAnalysisTaskHLTEMCAL, 0);
};

#endif
