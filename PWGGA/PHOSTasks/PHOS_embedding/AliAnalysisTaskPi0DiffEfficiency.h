#ifndef ALIANALYSISTASKPI0DIFFEFFICIENCY_H
#define ALIANALYSISTASKPI0DIFFEFFICIENCY_H

// Task calculating effciency of pi0 registration 
// as a difference between spectrum after and before 
// embedding

class TClonesArray;
class AliAODCaloCluster ;


#include "AliAnalysisTaskPi0Efficiency.h"

class AliAnalysisTaskPi0DiffEfficiency : public AliAnalysisTaskPi0Efficiency {
public:
  AliAnalysisTaskPi0DiffEfficiency(const char *name = "AliAnalysisTaskPi0DiffEfficiency");
  virtual ~AliAnalysisTaskPi0DiffEfficiency() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
private:
  AliAnalysisTaskPi0DiffEfficiency(const AliAnalysisTaskPi0DiffEfficiency&); // not implemented
  AliAnalysisTaskPi0DiffEfficiency& operator=(const AliAnalysisTaskPi0DiffEfficiency&); // not implemented
  Bool_t IsSameCluster(AliAODCaloCluster * c1,AliAODCaloCluster * c2)const ;
 
private:
  TClonesArray * fPHOSEvent1 ;      //!PHOS photons in current event
  TClonesArray * fPHOSEvent2 ;      //!PHOS photons in current event
 
  ClassDef(AliAnalysisTaskPi0DiffEfficiency, 2); // PHOS analysis task
};

#endif
