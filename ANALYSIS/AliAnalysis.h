#ifndef ALIANALYSIS_H
#define ALIANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliAnalysis
//
// Base class for analysis
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <TTask.h>

class AliAOD;
class AliStack;
 
class AliAnalysis: public TTask
{
  public: 
    AliAnalysis();
    AliAnalysis(const char* name,const char* title);
    virtual ~AliAnalysis();
    
    virtual Int_t Init() = 0;
    virtual Int_t ProcessEvent(AliAOD* aodrec, AliAOD* aodsim = 0x0) = 0;
    virtual Int_t Finish() = 0;
    
  protected:
    
  private:
    ClassDef(AliAnalysis,1)
};

#endif
