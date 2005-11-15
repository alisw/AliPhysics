#ifndef ALIFLOWANALYSIS_H
#define ALIFLOWANALYSIS_H

//*********************************************************
// class AliFlowAnalysis
// Flow Analysis
// S.Radomski@gsi.de
// Piotr.Skowronski@cern.ch
//*********************************************************

#include "AliAnalysis.h"

class AliESD;
class AliAOD;
class AliStack;
class AliAODParticleCut;
    
class AliFlowAnalysis: public AliAnalysis
{ 
  public: 
     AliFlowAnalysis();
     virtual ~AliFlowAnalysis();

    Int_t Init();
    Int_t ProcessEvent(AliAOD* aodrec, AliAOD* aodsim = 0x0);
    Int_t Finish();
   
    void SetParticleCut(AliAODParticleCut* pcut){fPartCut = pcut;}
    static Double_t GetEventPlane(AliESD* esd);
    static void     GetFlow(AliESD* esd,Double_t& v2,Double_t& psi);

    Double_t        GetEventPlane(AliAOD* aod);
    void            GetFlow(AliAOD* aod,Double_t& v2,Double_t& psi);
  protected:
    
  private:
    AliAODParticleCut* fPartCut;//Particle Cut
    ClassDef(AliFlowAnalysis,1)
};

#endif
