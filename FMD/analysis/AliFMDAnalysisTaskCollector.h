#ifndef ALIFMDANALYSISTASKCOLLECTOR_H
#define ALIFMDANALYSISTASKCOLLECTOR_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"
#include "TH1F.h"
#include "TObjArray.h"

class AliESDEvent;
class TChain;
class AliAODEvent;




class AliFMDAnalysisTaskCollector : public AliAnalysisTaskSE
{
 public:
    AliFMDAnalysisTaskCollector();
    AliFMDAnalysisTaskCollector(const char* name);
 AliFMDAnalysisTaskCollector(const AliFMDAnalysisTaskCollector& o) : AliAnalysisTaskSE(),
      fDebug(o.fDebug),
      fOutputList(o.fOutputList),
      fArray(o.fArray),
      fZvtxDist(o.fZvtxDist)  {}
    
    AliFMDAnalysisTaskCollector& operator=(const AliFMDAnalysisTaskCollector&) { return *this; }
    virtual ~AliFMDAnalysisTaskCollector() {;}
    // Implementation of interface methods
   
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    void ReadFromFile(const Char_t* filename, Bool_t store=kFALSE, Int_t speciesOption = 0);
   static Double_t  TripleLandau(Double_t *x, Double_t *par);
 private:
    void          GetVertex(Double_t* vertexXYZ); 
    Int_t         fDebug;        //  Debug flag
    TList*        fOutputList;
    TObjArray*    fArray;
    TH1F*         fZvtxDist;
   
    ClassDef(AliFMDAnalysisTaskCollector, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
