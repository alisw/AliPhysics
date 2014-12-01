#ifndef ALIEMCALTRIGGERMAKER_H
#define ALIEMCALTRIGGERMAKER_H

// $Id: AliEmcalTriggerMaker.h 64593 2013-10-18 10:23:58Z loizides $

class TClonesArray;
class AliEmcalTriggerSetupInfo;
class AliAODCaloTrigger;
class AliVVZERO;
class THistManager;

#include "AliLog.h"
#include "AliEmcalTriggerBitConfig.h"
#include "AliAnalysisTaskEmcal.h"


class AliEmcalTriggerMaker : public AliAnalysisTaskEmcal {
 public:
  enum TriggerMakerTriggerType_t {
    kTMEMCalJet = 0,
    kTMEMCalGamma = 1,
    kTMEMCalLevel0 = 2
  };
  enum TriggerMakerBitConfig_t {
    kOldConfig = 0,
    kNewConfig = 1
  };
  AliEmcalTriggerMaker();
  AliEmcalTriggerMaker(const char *name, Bool_t doQA = kFALSE);
  virtual ~AliEmcalTriggerMaker();

  void SetRunQA(Bool_t doQA = kTRUE) { fDoQA = doQA; }
  void SetCaloTriggersOutName(const char *name)     { fCaloTriggersOutName      = name; }
  void SetCaloTriggerSetupOutName(const char *name) { fCaloTriggerSetupOutName  = name; }
  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[2][0] = a; fThresholdConstants[2][1] = b; fThresholdConstants[2][2] = c; }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[0][0] = a; fThresholdConstants[0][1] = b; fThresholdConstants[0][2] = c; }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[3][1] = b; fThresholdConstants[3][2] = c; }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) { fThresholdConstants[3][0] = a; fThresholdConstants[1][1] = b; fThresholdConstants[1][2] = c; }
  void SetV0InName(const char *name) { fV0InName      = name; }

  void SetRunTriggerType(TriggerMakerTriggerType_t type, Bool_t doTrigger = kTRUE) { fRunTriggerType[type] = doTrigger; }
  void SetUseTriggerBitConfig(TriggerMakerBitConfig_t bitConfig) { fUseTriggerBitConfig = bitConfig; }
  void SetTriggerBitConfig(const AliEmcalTriggerBitConfig *conf) { fTriggerBitConfig = conf; }

  Bool_t IsEJE(Int_t tBits) const { if( tBits & ( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetHighBit()) | 1 << (fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetJetHighBit()) )) return kTRUE; else return kFALSE; }
  Bool_t IsEGA(Int_t tBits) const { if( tBits & ( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaHighBit()) | 1 << (fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetGammaHighBit()) )) return kTRUE; else return kFALSE; }
  Bool_t IsLevel0(Int_t tBits) const { if( tBits & (1 << (fTriggerBitConfig->GetLevel0Bit() + fTriggerBitConfig->GetLevel0Bit()) | (1 << fTriggerBitConfig->GetLevel0Bit()))) return kTRUE; return kFALSE; }

 protected:  
  enum{
	  kPatchCols = 48,
	  kPatchRows = 64
  };
  void                       UserCreateOutputObjects();
  void                       ExecOnce();
  Bool_t                     Run();
  void                       RunSimpleOfflineTrigger();
  Bool_t                     NextTrigger( Bool_t &isOfflineSimple );
  AliEmcalTriggerPatchInfo*  ProcessPatch(TriggerMakerTriggerType_t type, Bool_t isOfflineSimple);
  Bool_t 					           CheckForL0(const AliVCaloTrigger &trg) const;

  TString                    fCaloTriggersOutName;      // name of output track array
  TString                    fCaloTriggerSetupOutName;  // name of output track array
  TString                    fV0InName;                 // name of output track array
  TriggerMakerBitConfig_t    fUseTriggerBitConfig;      // type of trigger config
  Int_t                      fThresholdConstants[4][3]; // simple offline trigger thresholds constants
  const AliEmcalTriggerBitConfig  *fTriggerBitConfig;         // Trigger bit configuration, aliroot-dependent
  TClonesArray              *fCaloTriggersOut;          //!trigger array out
  AliEmcalTriggerSetupInfo  *fCaloTriggerSetupOut;      //!trigger setup
  AliAODCaloTrigger         *fSimpleOfflineTriggers;    //!simple offline trigger
  AliVVZERO                 *fV0;                       //!V0 object
  Double_t                   fPatchADCSimple[kPatchCols][kPatchRows];   //!patch map for simple offline trigger
  Int_t                      fPatchADC[kPatchCols][kPatchRows];         //!ADC values map
  Int_t                      fITrigger;                 //!trigger counter
  Bool_t                     fRunTriggerType[3];        // Run patch maker for a given trigger type
  Bool_t                     fDoQA;                     // Fill QA histograms
  THistManager              *fQAHistos;                 //! Histograms for QA

 private:
  AliEmcalTriggerMaker(const AliEmcalTriggerMaker&);            // not implemented
  AliEmcalTriggerMaker &operator=(const AliEmcalTriggerMaker&); // not implemented

  ClassDef(AliEmcalTriggerMaker, 4) // Task to make array of EMCAL trigger patches
};
#endif
