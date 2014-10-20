#ifndef ALITPCPIDETAQA_H
#define ALITPCPIDETAQA_H

/*
This task determines the eta dependence of the TPC signal.
For this purpose, only tracks fulfilling some standard quality cuts are taken into account.
The obtained data can be used to derive the functional behaviour of the eta dependence.
Such a function can be plugged into this task to correct for the eta dependence and to see
if there is then really no eta dependence left.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

class TH1F;
class TF1;
class THnSparse;
class AliPIDResponse;
class AliPID;

#include "AliTPCPIDBase.h"

class AliTPCPIDEtaQA : public AliTPCPIDBase {
 public:
  AliTPCPIDEtaQA();
  AliTPCPIDEtaQA(const char *name);
  virtual ~AliTPCPIDEtaQA();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
  Double_t GetPtThresholdForPhiCut() const { return fPtThresholdForPhiCut; };     
  virtual void  SetPtThresholdForPhiCut(Double_t threshold) { fPtThresholdForPhiCut = threshold; };
  
 protected:
  virtual void   SetUpHist(THnSparse* hist, Double_t* binsPt) const;
  
 private:
  Double_t fPtThresholdForPhiCut; // Only apply phi cut on tracks with pT larger equal this threshold
  
  TF1* fPhiCutSecondBranchLow; // phi prime cut, low, second branch (very low pT)
  TF1* fPhiCutSecondBranchHigh; // phi prime cut, high, second branch (very low pT)
    
  THnSparseI* fhPIDdataAll; //! data histogram
  
  //OLD clusterQA THnSparseI* fhNumClustersPhiPrimePtBeforeCut; //! QA histogra - before phi prime cut
  //OLD clusterQA THnSparseI* fhNumClustersPhiPrimePtAfterCut; //! QA histogra - after phi prime cut
  
  TH1F* fhPhiPrimeCutEfficiency; //! Effeciency of phiPrime cut as a functio of pT
  
  TObjArray* fOutputContainer; //! output data container
   
  AliTPCPIDEtaQA(const AliTPCPIDEtaQA&); // not implemented
  AliTPCPIDEtaQA& operator=(const AliTPCPIDEtaQA&); // not implemented
  
  ClassDef(AliTPCPIDEtaQA, 1); // example of analysis
};

#endif
