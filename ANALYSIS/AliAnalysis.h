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
class AliEventCut;
 
class AliAnalysis: public TTask
{
  public: 
    AliAnalysis();
    AliAnalysis(const char* name,const char* title);
    virtual ~AliAnalysis();
    
    virtual Int_t Init() = 0;
    virtual Int_t ProcessEvent(AliAOD* aodrec, AliAOD* aodsim = 0x0) = 0;
    virtual Int_t Finish() = 0;

    void          EventCutOnRec(Bool_t flag){fCutOnRec = flag;}
    void          EventCutOnSim(Bool_t flag){fCutOnSim = flag;}
    void          SetEventCut(AliEventCut* evcut);
    
  protected:
    Bool_t        Pass(AliAOD* recevent, AliAOD* simevent);
    AliEventCut*  fEventCut;//event cut

    Bool_t        fCutOnSim;//flag indicating that event cut is performed on simulated particles 
    Bool_t        fCutOnRec;//flag indicating that event cut is performed on reconstructed tracks
    
  private:
    ClassDef(AliAnalysis,1)
};

#endif
