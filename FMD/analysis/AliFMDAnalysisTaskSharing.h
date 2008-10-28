#ifndef ALIFMDANALYSISTASKSHARING_H
#define ALIFMDANALYSISTASKSHARING_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"
#include "TH1F.h"
#include "TObjString.h"
#include "AliESDFMD.h"
#include "TTree.h"
class AliESDEvent;
class TChain;
class AliAODEvent;



class AliFMDAnalysisTaskSharing : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskSharing();
    AliFMDAnalysisTaskSharing(const char* name);
    virtual ~AliFMDAnalysisTaskSharing() {;}
 AliFMDAnalysisTaskSharing(const AliFMDAnalysisTaskSharing& o) : AliAnalysisTask(),
      fDebug(o.fDebug),
      fESD(o.fESD),
      fOutputESD(o.fOutputESD),
      foutputESDFMD(o.foutputESDFMD),
      fSharedThis(o.fSharedThis),
      fSharedPrev(o.fSharedPrev) {}
    AliFMDAnalysisTaskSharing& operator=(const AliFMDAnalysisTaskSharing&) { return *this; }
    
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init() {}
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t* /* option*/) {}
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    Float_t GetMultiplicityOfStrip(Float_t mult, Float_t Eprev, Float_t Enext, Int_t   det, Char_t  ring);
 private:
    Int_t         fDebug;        //  Debug flag
    AliESDEvent*  fESD;          //! ESD
    AliESDEvent*  fOutputESD;
    AliESDFMD*    foutputESDFMD;
    Bool_t        fSharedThis;
    Bool_t        fSharedPrev;
    ClassDef(AliFMDAnalysisTaskSharing, 0); // Analysis task for FMD analysis
};
 
#endif
