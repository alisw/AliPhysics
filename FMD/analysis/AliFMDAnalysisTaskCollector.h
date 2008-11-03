#ifndef ALIFMDANALYSISTASKCOLLECTOR_H
#define ALIFMDANALYSISTASKCOLLECTOR_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"
#include "TH1F.h"
#include "TObjArray.h"

class AliESDEvent;
class TChain;
class AliAODEvent;



class AliFMDAnalysisTaskCollector : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskCollector();
    AliFMDAnalysisTaskCollector(const char* name);
    virtual ~AliFMDAnalysisTaskCollector() {;}
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
 private:
    Int_t         fDebug;        //  Debug flag
    TChain*       fChain;        //! chained files
    AliESDEvent*  fESD;          //! ESD
    AliAODEvent*  fAOD;          //! AOD
    TList*        fOutputList;
    TObjArray*    fArray;
    TH1F*         fEdistHist;
    TH1F*         fZvtxDist;
   
    ClassDef(AliFMDAnalysisTaskCollector, 0); // Analysis task for FMD analysis
};
 
#endif
