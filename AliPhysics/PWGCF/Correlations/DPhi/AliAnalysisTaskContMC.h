#ifndef ALIANALYSISTASKCONTMC_H
#define ALIANALYSISTASKCONTMC_H

class TH1F;
class TH2F;
class TH3F;
class AliAODEvent;
class AliHelperPID;

#include "AliAnalysisTaskSE.h"
#include "AliHelperPID.h"

class AliAnalysisTaskContMC : public AliAnalysisTaskSE
{
 public:
  
  // constructors
 AliAnalysisTaskContMC() : AliAnalysisTaskSE(), fAOD(0), fNSigmaPID(0), fIsMC(0), fOutput(0), fHistID(0)
    {}
  AliAnalysisTaskContMC(const char *name);
  virtual ~AliAnalysisTaskContMC() {}
  
  AliHelperPID * GetPID()         {  return fNSigmaPID; }
  void SetPID      (AliHelperPID      * pid)   {   fNSigmaPID  = pid;}
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; }
  Bool_t GetIsMC()           const           { return fIsMC;}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  
  AliAODEvent           *fAOD;         //! AOD object
  AliHelperPID           *fNSigmaPID;  // NSigmaPID object
  Bool_t          fIsMC;// true if processing MC
  TList *fOutput; //! tlist with output
  TH3F *fHistID; //! histo 
 
  AliAnalysisTaskContMC(const AliAnalysisTaskContMC&);
  AliAnalysisTaskContMC& operator=(const AliAnalysisTaskContMC&);
  
  ClassDef(AliAnalysisTaskContMC, 1);
};

#endif
