#ifndef ALIRUNANALYSIS_H
#define ALIRUNANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliRunAnalysis
//
// Analysis manager
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <TTask.h>
#include <TObjArray.h>
#include "AliAnalysis.h"

class AliEventCut;
class TFile;
class AliReader;

class AliRunAnalysis: public TTask
{
  public: 
    AliRunAnalysis();
    virtual ~AliRunAnalysis();
    
    Int_t         Run();
    void          Add(AliAnalysis* a);
    void          SetReader(AliReader* reader){fReader = reader;}
    
    const char*   GetName(){return "RunAnalysis";}
    
  protected:
    TObjArray     fAnalysies;//arry with analysies
    AliReader*    fReader;//arry with directories to read data from
    
    AliEventCut*  fEventCut;//event cut    
    
    Bool_t        fCutOnSim;//flag indicating that event cut is performed on simulated particles 
    Bool_t        fCutOnRec;//flag indicating that event cut is performed on reconstructed tracks
    
    Bool_t        Pass(AliAOD* recevent, AliAOD* simevent);
    
  private:
    void SetName(const char *){}//change SetName to be private
    
    ClassDef(AliRunAnalysis,1)
};

#endif
