#ifndef AliAnalysisTaskEmcalTriggerInfoQA_H
#define AliAnalysisTaskEmcalTriggerInfoQA_H

class TH1;
class TList;
class TClonesArray;
class TString;
class AliEMCALTriggerPatchInfo;
class AliEmcalTriggerSetupInfo;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskEmcalTriggerInfoQA : public AliAnalysisTaskEmcal
{
  public:
    AliAnalysisTaskEmcalTriggerInfoQA();
    AliAnalysisTaskEmcalTriggerInfoQA(const char *name);
    virtual ~AliAnalysisTaskEmcalTriggerInfoQA();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExecOnce();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

  void SetCaloTriggerPatchInfoName(const char *name) { fCaloTriggerPatchInfoName      = name; }
  void SetCaloTriggerSetupInfoName(const char *name) { fCaloTriggerSetupInfoName      = name; }

  private:
    TList *fOutput;          //! Output list
    TH1 **fHistos;           //! histos
    TClonesArray *fTriggersInfo;     //! jet array
    AliEmcalTriggerSetupInfo *fTriggerSetup;   //! tracks array
    
    Bool_t fIsInitialized;  //! init flag
    
    TString    fCaloTriggerPatchInfoName;      // trigger array name
    TString    fCaloTriggerSetupInfoName;  // track bins
    

    void FillPatch( AliEMCALTriggerPatchInfo *patch, Int_t type );

    AliAnalysisTaskEmcalTriggerInfoQA(const AliAnalysisTaskEmcalTriggerInfoQA&); // not implemented
    AliAnalysisTaskEmcalTriggerInfoQA& operator=(const AliAnalysisTaskEmcalTriggerInfoQA&); // not implemented
    
    ClassDef(AliAnalysisTaskEmcalTriggerInfoQA, 2); // example of analysis
};
#endif
