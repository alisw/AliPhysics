#ifndef AliAnalysisTaskggMC_cxx
#define AliAnalysisTaskggMC_cxx

// task to analyze background correlations in gg interferometry
// Author: D.Peresunko

class THashList ;
class AliPHOSGeometry;
class AliCaloPhoton ;
class AliAODMCParticle ;

#include "AliAnalysisTaskgg.h"

class AliAnalysisTaskggMC : public AliAnalysisTaskgg {
public:
  
  enum CutList {kDefault, kDisp, kCPV, kBoth, kDistance1, kDistance2, kDistance3} ;  
  
  
  AliAnalysisTaskggMC(const char *name = "AliAnalysisTaskggMC");
  virtual ~AliAnalysisTaskggMC() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
protected:
  void ProcessMC() ;  
  Int_t FindAODLabel(Int_t esdLabel)const ;
 
  
private:
  AliAnalysisTaskggMC(const AliAnalysisTaskggMC&); // not implemented
  AliAnalysisTaskggMC& operator=(const AliAnalysisTaskggMC&); // not implemented

  Double_t PionHBTWeight(const AliAODMCParticle * ph1, const AliAODMCParticle * ph2) const ;
  Double_t FlowWeight(const AliAODMCParticle * ph1, const AliAODMCParticle * ph2) const ;
  Double_t EtaPhiWeight(const AliAODMCParticle * ph1, const AliAODMCParticle * ph2) const ;
  AliAODMCParticle* TestCommonParent(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2) const ;
 
private:
  TClonesArray * fStack ;
  TList *       fMCEvents ;   //Containers for events with primary photons
  TClonesArray* fMCEvent ;    //primary photons in current event
  TList *       fMCEventsH ;   //Containers for events with primary hadrons
  TClonesArray* fMCEventH ;    //primary hadrons in current event

  
  
  ClassDef(AliAnalysisTaskggMC, 1); // PHOS analysis task
};

#endif
