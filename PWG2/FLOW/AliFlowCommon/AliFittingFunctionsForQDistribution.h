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
 *                              *
 * based on the macro written   *
 *     by Sergei Voloshin       *                        
 *******************************/ 

#ifndef AliFittingFunctionsForQDistribution_H
#define AliFittingFunctionsForQDistribution_H

#include "AliFlowCommonConstants.h"

class TH1;
class TProfile;

class TObjArray;
class TList;
class TFile;

//================================================================================================================

class AliFittingFunctionsForQDistribution{
 public:
  AliFittingFunctionsForQDistribution();
  virtual ~AliFittingFunctionsForQDistribution();
  AliFittingFunctionsForQDistribution(TProfile *AvMult, TH1D *QDistribution, TH1D *intFlowRes, TH1D *sigma2, AliFlowCommonHistResults *chr);
 
  void Calculate();

 private:
  AliFittingFunctionsForQDistribution(const AliFittingFunctionsForQDistribution& fun);
  AliFittingFunctionsForQDistribution& operator=(const AliFittingFunctionsForQDistribution& fun);
  
  TProfile                 *fAvMultFQD;           //avarage selected multiplicity for int. flow
  TH1D                     *fQDistributionFQD;    //q-distribution
  TH1D                     *fIntFlowResFQD;       //integrated flow final result
  TH1D                     *fSigma2;              //sigma^2
  AliFlowCommonHistResults *fchrFQD;              //final results for integrated flow stored in the common histograms

  ClassDef(AliFittingFunctionsForQDistribution, 0);
};

//================================================================================================================

#endif





