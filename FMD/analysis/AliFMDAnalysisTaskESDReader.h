#ifndef ALIFMDANALYSISTASKESDREADER_H
#define ALIFMDANALYSISTASKESDREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"
#include "TH1F.h"
#include "TObjString.h"
#include "AliESDFMD.h"
#include "TTree.h"
#include "AliESDVertex.h"
class AliESDEvent;
class TChain;
class AliAODEvent;



class AliFMDAnalysisTaskESDReader : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskESDReader();
    AliFMDAnalysisTaskESDReader(const char* name);
 AliFMDAnalysisTaskESDReader(const AliFMDAnalysisTaskESDReader& o) : AliAnalysisTask(),
      fDebug(o.fDebug),fChain(o.fChain), fESD(o.fESD),fOutputESD(o.fOutputESD) {}
    
    virtual ~AliFMDAnalysisTaskESDReader() {;}
    AliFMDAnalysisTaskESDReader& operator=(const AliFMDAnalysisTaskESDReader&) { return *this; }
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option );
    virtual void CreateOutputObjects() {};
    virtual void Init() {}
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t* /* option*/) {}
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
 private:
    Int_t         fDebug;        //  Debug flag
    TChain*       fChain;        //! chained files
    AliESDEvent*  fESD;          //! ESD
    AliESDEvent*  fOutputESD;
    
    ClassDef(AliFMDAnalysisTaskESDReader, 0); // Analysis task for FMD analysis
};
 
#endif
