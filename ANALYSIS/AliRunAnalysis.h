#ifndef ALIRUNANALYSIS_H
#define ALIRUNANALYSIS_H

///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
// Analysis manager
// Author: Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <TTask.h>
#include <TObjArray.h>

class AliEventCut;
class AliReader;
class AliAnalysis;
class AliAOD;

class AliRunAnalysis: public TTask
{
  public: 
    AliRunAnalysis();
    virtual ~AliRunAnalysis();
    
    Int_t         Run();
    void          Add(TTask *t){TTask::Add(t);}
    void          Add(AliAnalysis* a);
    void          SetReader(AliReader* reader){fReader = reader;}
    
    const char*   GetName() const {return "RunAnalysis";}
    void          EventCutOnRec(Bool_t flag){fCutOnRec = flag;}
    void          EventCutOnSim(Bool_t flag){fCutOnSim = flag;}
    void          SetEventCut(AliEventCut* evcut);
    void          SetOwner(Bool_t owner=kTRUE){fAnalysies.SetOwner(owner);}
    
  protected:
    TObjArray     fAnalysies;//arry with analysies
    AliReader*    fReader;//arry with directories to read data from
    
    AliEventCut*  fEventCut;//event cut    
    
    Bool_t        fCutOnSim;//flag indicating that event cut is performed on simulated particles 
    Bool_t        fCutOnRec;//flag indicating that event cut is performed on reconstructed tracks
    
    Bool_t        Rejected(AliAOD* recevent, AliAOD* simevent);
    
  private:
    AliRunAnalysis(const AliRunAnalysis & src);
    AliRunAnalysis & operator=(const AliRunAnalysis & src);

    void SetName(const char *){}//change SetName to be private
    
    ClassDef(AliRunAnalysis,1)
};

#endif
