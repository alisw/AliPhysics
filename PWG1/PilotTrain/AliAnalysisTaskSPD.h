#ifndef ALIANALYSISTASKSPD_H
#define ALIANALYSISTASKSPD_H


#include "AliAnalysisTaskSE.h"

class AliITSsegmentationSPD;
class AliAnalysisTaskSPD : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskSPD();
  AliAnalysisTaskSPD(const char *name);
  virtual ~AliAnalysisTaskSPD();

  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  

 private:

  UInt_t GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip); //see AliITSRawStreamSPD
  UInt_t GetOfflineChipKeyFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip); // see AliITSRawStreamSPD
  UInt_t GetOnlineEqIdFromOffline(UInt_t module); // see AliITSRawStreamSPD
  UInt_t GetOnlineHSFromOffline(UInt_t module); // see AliITSRawStreamSPD
  UInt_t GetOnlineChipFromOffline(UInt_t module, UInt_t colM); // see AliITSRawStreamSPD
  Int_t GetModuleNumber(UInt_t iDDL, UInt_t iModule); // see AliITSRawStreamSPD

  AliAnalysisTaskSPD(const AliAnalysisTaskSPD &source);
  AliAnalysisTaskSPD& operator=(const AliAnalysisTaskSPD &source);

  static const Int_t fgkDDLModuleMap[20][12];  // mapping DDL/module -> module number

  TList   *fOutput; //! list of histos
  AliITSsegmentationSPD *fSegSPD;
  ClassDef(AliAnalysisTaskSPD,1);  
};


#endif
