/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/******************************** 
 * integrated flow estimate by  *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
 *******************************/ 
 
#ifndef AliFittingQDistribution_H
#define AliFittingQDistribution_H

#include "AliFlowCommonConstants.h"

class TObjArray;
class TList;
class TFile;

class TH1;
class TProfile;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

//================================================================================================================

class AliFittingQDistribution{
 public:
  AliFittingQDistribution();
  virtual ~AliFittingQDistribution(); 
  
  virtual void CreateOutputObjects();
  virtual void Make(AliFlowEventSimple* anEvent);
  virtual void Finish();
  virtual void WriteHistograms(TString* outputFileName);
  
  TList* GetHistList() const {return this->fHistList;}      //output histogram list
 
 private:
  AliFittingQDistribution(const AliFittingQDistribution& afqd);
  AliFittingQDistribution& operator=(const AliFittingQDistribution& afqd);
  AliFlowTrackSimple* fTrack;                               //track
  
  TList*             fHistList;         //list to hold all output histograms
  TProfile*          fAvMultIntFlowFQD;    //avarage selected multiplicity
  TH1D*              fIntFlowResultsFQD;   //integrated flow final results
  AliFlowCommonHist* fCommonHists;      //common control histograms
  TH1D*              fQDistributionFQD;            //q-distribution
      
  ClassDef(AliFittingQDistribution, 0);
};

//================================================================================================================

#endif





