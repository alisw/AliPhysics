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
  
  virtual void Init();
  virtual void Make(AliFlowEventSimple* anEvent);
  virtual void Finish();
  virtual void WriteHistograms(TString* outputFileName);
  virtual void WriteHistograms(TString outputFileName);

//----------------------------------------------------------------------------------------------------------------
//                                            setters and getters                                                 
//----------------------------------------------------------------------------------------------------------------      
  TList* GetHistList() const {return this->fHistList;}
  
  void SetWeightsList(TList* wlist) {this->fWeightsList = wlist;}
  TList* GetWeightsList() const {return this->fWeightsList;}   
  
  void SetIntFlowResults(TH1D* ifr)  {this->fIntFlowResultsFQD = ifr;};
  TH1D* GetIntFlowResults() const    {return this->fIntFlowResultsFQD;};
  
  void SetCommonHistsResults(AliFlowCommonHistResults* chr)  {this->fCommonHistsResults = chr;};
  AliFlowCommonHistResults* GetCommonHistsResults() const    {return this->fCommonHistsResults;};
  
  void SetAverageMultiplicity(TProfile* am)  {this->fAvMultIntFlowFQD = am;};
  TProfile* GetAverageMultiplicity() const   {return this->fAvMultIntFlowFQD;};
  
  void SetQDistribution(TH1D* qd)  {this->fQDistributionFQD = qd;};
  TH1D* GetQDistribution() const   {return this->fQDistributionFQD;};
  
  void SetSigma2(TH1D* s2)  {this->fSigma2 = s2;};
  TH1D* GetSigma2() const   {return this->fSigma2;};
  
  void SetUsePhiWeights(Bool_t const uPhiW) {this->fUsePhiWeights = uPhiW;};
  Bool_t GetUsePhiWeights() const {return this->fUsePhiWeights;};
  
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  
  void SetUseEtaWeights(Bool_t const uEtaW) {this->fUseEtaWeights = uEtaW;};
  Bool_t GetUseEtaWeights() const {return this->fUseEtaWeights;};
//----------------------------------------------------------------------------------------------------------------
 
 private:
  AliFittingQDistribution(const AliFittingQDistribution& afqd);
  AliFittingQDistribution& operator=(const AliFittingQDistribution& afqd);
  
  AliFlowTrackSimple*        fTrack;                   //track
   
  TList*                     fHistList;                //list to hold all output histograms
  TList*                     fWeightsList;             //list to hold all histograms with weights
  
  TProfile*                  fAvMultIntFlowFQD;        //avarage selected multiplicity
  TH1D*                      fIntFlowResultsFQD;       //integrated flow final results
  TH1D*                      fSigma2;                  //sigma^2
  AliFlowCommonHist*         fCommonHists;             //common control histograms
  AliFlowCommonHistResults*  fCommonHistsResults;      //final results for integrated flow stored in the common histograms 
  TH1D*                      fQDistributionFQD;        //q-distribution 
  
  Bool_t                     fUsePhiWeights;           //phi weights
  Bool_t                     fUsePtWeights;            //v_2(pt) weights
  Bool_t                     fUseEtaWeights;           //v_2(eta) weights
         
  ClassDef(AliFittingQDistribution, 0);
};

//================================================================================================================

#endif





