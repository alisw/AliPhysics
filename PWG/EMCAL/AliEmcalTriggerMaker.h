#ifndef ALIEMCALTRIGGERMAKER_H
#define ALIEMCALTRIGGERMAKER_H

// $Id$

class TClonesArray;
class AliEmcalTriggerSetupInfo;
class AliAODCaloTrigger;
class AliVVZERO;

#include "AliAnalysisTaskEmcal.h"

class AliEmcalTriggerMaker : public AliAnalysisTaskEmcal {
 public:
  AliEmcalTriggerMaker();
  AliEmcalTriggerMaker(const char *name);
  virtual ~AliEmcalTriggerMaker();

  void ExecOnce();
  Bool_t Run();
	void RunSimpleOfflineTrigger();

  void SetCaloTriggersOutName(const char *name) { fCaloTriggersOutName      = name; }
  void SetCaloTriggerSetupOutName(const char *name) { fCaloTriggerSetupOutName      = name; }

  void SetTriggerThresholdJetLow( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  
  void SetV0InName(const char *name) { fV0InName      = name; }

  Int_t IsEGA(Int_t level) {if (level > 1 || level < 0) AliError("EGA: check the requested threshold"); return fEGA[level];}
  Int_t IsEJE(Int_t level) {if (level > 1 || level < 0) AliError("EJE: check the requested threshold"); return fEJE[level];}
  
 protected:  
  TString            fCaloTriggersOutName;    // name of output track array
  TString            fCaloTriggerSetupOutName;    // name of output track array
  TString            fV0InName;    // name of output track array
  TClonesArray      *fCaloTriggersOut;        //!trigger array out
  AliEmcalTriggerSetupInfo  *fCaloTriggerSetupOut;        //!trigger setup
  AliAODCaloTrigger *fSimpleOfflineTriggers; //! simple offline trigger
  AliVVZERO         *fV0;                    //! V0 object
  
  Double_t           fPatchADCSimple[48][64]; //! patch map for simple offline trigger
  Int_t              fThresholdConstants[4][3]; // simple offline trigger thresholds constants

  Int_t              fEGA[2];
  Int_t              fEJE[2];
  
 private:
  AliEmcalTriggerMaker(const AliEmcalTriggerMaker&);            // not implemented
  AliEmcalTriggerMaker &operator=(const AliEmcalTriggerMaker&); // not implemented

  Bool_t NextTrigger( Bool_t &isOfflineSimple );

  ClassDef(AliEmcalTriggerMaker, 3); // Task to make array of EMCAL particle
};
#endif
