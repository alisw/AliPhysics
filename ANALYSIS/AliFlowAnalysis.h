#ifndef ALIFLOWANALYSIS_H
#define ALIFLOWANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliFlowAnalysis
//
// Flow Analysis
//
//
// S.Radomski@gsi.de
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include "AliAnalysis.h"

class AliESD;
class AliStack;

class AliFlowAnalysis: public AliAnalysis
{ 
  public: 
     AliFlowAnalysis(){}
     ~AliFlowAnalysis(){}

    Int_t Init();
    Int_t ProcessEvent(AliESD* esd, AliStack* stack = 0x0);
    Int_t Finish();
   
    static Double_t GetEventPlane(AliESD* esd);
       
  protected:

  private:
 
    ClassDef(AliFlowAnalysis,1)
};

#endif
