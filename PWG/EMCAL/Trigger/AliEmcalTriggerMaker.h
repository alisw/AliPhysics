/**
 * \file AliEmcalTriggerMaker.h
 * \brief Class to make array of trigger patch objects in AOD/ESD events.
 *
 * Class to make array of trigger patch objects in AOD/ESD events.
 * The input for the process are:
 *   - AliCaloTrigger objects from ESS/AOD, which contain raw trigger information
 *   - the CaloCells, which contain offline/FEE information
 *
 * The output is a list of AliEmcalTriggerPatchInfo objects which is stored in ESD/AOD (Use event->FindListObject to get them) with three types of trigger patches:
 *  -# Online trigger info
 *  -# Trigger info based on the offline FEE ADCs (SimpleOffline)
 *  -# The highest gamma and jet patch in the event, even if it does
 *     not pass the threshold (RecalcJet and RecalcGamma); with two versions
 *     -# based on the online trigger information
 *     -# based offline FEE information
 * The different types of patches are distinguished by bitflags according
 * to the enum AliEmcalTriggerPatchInfo::TriggerMakerBits and
 *   EMCAL/AliEmcalTriggerTypes.h
 *
 * \author Jiri Kral <>, University of Jyv&aumlskul&auml
 * \date Jun 26, 2013
 */
#ifndef ALIEMCALTRIGGERMAKER_H
#define ALIEMCALTRIGGERMAKER_H

/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class AliEmcalTriggerSetupInfo;
class AliAODCaloTrigger;
class AliVVZERO;
class THistManager;

template<class T> class AliEmcalTriggerDataGridAP;

#include "AliLog.h"
#include "AliEmcalTriggerBitConfig.h"
#include "AliEmcalTriggerChannelContainerAP.h"
#include "AliAnalysisTaskEmcal.h"

/**
 * \class AliEmcalTriggerMaker
 * \brief EMCAL trigger patch maker
 *
 * This class creates EMCAL trigger patches found in reconstructed events.
 * Two methods can be applied.
 * -# Trigger patches can be created from reconstructed online trigger patch
 * information using the online TRU amplitude.
 * -# A simple offline trigger estimates a TRU amplitude from calibrated
 * EMCAL cells, and runs the same trigger algorithm on the offline estimated
 * amplitudes.
 *
 * For reconstructed trigger patches, an object AliEmcalTriggerPatchInfo which
 * stores all relevant information of the trigger patch is created and stored
 * in a TClonesArray. This TClonesArray is appended to the ESD or AOD event,
 * and can be retrieved in other analysis tasks in the same train by the container
 * name, defined in the function SetCaloTriggerOutName.
 */
class AliEmcalTriggerMaker : public AliAnalysisTaskEmcal {
 public:
  /**
   * \enum TriggerMakerTriggerType_t
   * \brief Definition of different trigger patch types
   *
   * This enumeration defines the different trigger patch types
   * processed by the trigger maker. Each trigger patch type has
   * a certain patch size and therefore a certain length and
   * geometric center
   */
  enum TriggerMakerTriggerType_t {
    kTMEMCalJet = 0,       ///< EMCAL Jet trigger
    kTMEMCalGamma = 1,     ///< EMCAL Gamma trigger
    kTMEMCalLevel0 = 2,    ///< EMCAL Level0 patches
    kTMEMCalRecalcJet = 3, ///< EMCAL Jet patches, recalculated
    kTMEMCalRecalcGamma = 4///< EMCAL Gamma patches, recalculated
  };

  /***
   * \enum TriggerMakerTriggerBitConfig_t
   * \brief Definition of trigger bit configurations
   *
   * This enumeration handles different trigger bit configurations for the
   * EMCAL Level1 triggers (with and without different thresholds) applied
   * in the reconstruction of different samples.
   */
  enum TriggerMakerBitConfig_t {
    kOldConfig = 0,///< Old configuration, no distinction between high and low threshold
    kNewConfig = 1 ///< New configuration, distiction between high and low threshold
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
  void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel; }
  void SetV0InName(const char *name) { fV0InName      = name; }
  void AddHotFastor(int col, int row) { fBadChannels.AddChannel(col, row); }

  void SetRunTriggerType(TriggerMakerTriggerType_t type, Bool_t doTrigger = kTRUE) { fRunTriggerType[type] = doTrigger; }
  void SetUseTriggerBitConfig(TriggerMakerBitConfig_t bitConfig) { fUseTriggerBitConfig = bitConfig; }
  void SetTriggerBitConfig(const AliEmcalTriggerBitConfig *conf) { fTriggerBitConfig = conf; }

  /**
   * Switch on rejection of patches which leave the EMCAL acceptance in \f$ \eta \f$ and \f$ \phi \f$
   * \param doReject If true we reject patches outside the EMCAL acceptance
   */
  void SetRejectOffAcceptancePatches(Bool_t doReject = kTRUE) { fRejectOffAcceptancePatches = doReject; }

  inline Bool_t IsEJE(Int_t tBits) const;
  inline Bool_t IsEGA(Int_t tBits) const;
  inline Bool_t IsLevel0(Int_t tBits) const;

 protected:
  static const TString                      fgkTriggerTypeNames[5];       ///< Histogram name tags
  static const int                          kColsEta;                     ///< Number of columns in eta direction

  void                                      UserCreateOutputObjects();
  void                                      ExecOnce();
  Bool_t                                    Run();
  void                                      RunSimpleOfflineTrigger();
  Bool_t                                    NextTrigger( Bool_t &isOfflineSimple );
  AliEmcalTriggerPatchInfo*                 ProcessPatch(TriggerMakerTriggerType_t type, Bool_t isOfflineSimple);
  Bool_t 					                          CheckForL0(const AliVCaloTrigger &trg) const;

  AliEmcalTriggerChannelContainerAP         fBadChannels;                 ///< Container of bad channels
  TString                                   fCaloTriggersOutName;         ///< name of output track array
  TString                                   fCaloTriggerSetupOutName;     ///< name of output track array
  TString                                   fV0InName;                    ///< name of output track array
  TriggerMakerBitConfig_t                   fUseTriggerBitConfig;         ///< type of trigger config
  Int_t                                     fThresholdConstants[4][3];    ///< simple offline trigger thresholds constants
  const AliEmcalTriggerBitConfig            *fTriggerBitConfig;           ///< Trigger bit configuration, aliroot-dependent
  TClonesArray                              *fCaloTriggersOut;            //!<! trigger array out
  AliEmcalTriggerSetupInfo                  *fCaloTriggerSetupOut;        //!<! trigger setup
  AliAODCaloTrigger                         *fSimpleOfflineTriggers;      //!<! simple offline trigger
  AliVVZERO                                 *fV0;                         //!<! V0 object
  AliEmcalTriggerDataGridAP<float>          *fPatchAmplitudes;            //!<! TRU Amplitudes (for L0)
  AliEmcalTriggerDataGridAP<double>         *fPatchADCSimple;             //!<! patch map for simple offline trigger
  AliEmcalTriggerDataGridAP<int>            *fPatchADC;                   //!<! ADC values map
  AliEmcalTriggerDataGridAP<char>           *fLevel0TimeMap;              //!<! Map needed to store the level0 times
  Int_t                                     fITrigger;                    //!<! trigger counter
  Bool_t                                    fRunTriggerType[5];           ///< Run patch maker for a given trigger type
  Bool_t                                    fDoQA;                        ///< Fill QA histograms
  Bool_t                                    fRejectOffAcceptancePatches;  ///< Switch for rejection of patches outside the acceptance
  THistManager                              *fQAHistos;                   //!<! Histograms for QA

  Int_t                                     fDebugLevel;                  ///< Debug lebel;

 private:
  AliEmcalTriggerMaker(const AliEmcalTriggerMaker&);            // not implemented
  AliEmcalTriggerMaker &operator=(const AliEmcalTriggerMaker&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMaker, 6) // Task to make array of EMCAL trigger patches
  /// \endcond
};

/**
 * Check whehter trigger is jet patch trigger according to trigger bits
 * @param tBits Trigger bits of the fastor
 * @return True if fastor belongs to a jet patch trigger, false otherwise
 */
Bool_t AliEmcalTriggerMaker::IsEJE(Int_t tBits) const {
  if( tBits & ( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetJetHighBit()) | 1 << (fTriggerBitConfig->GetJetLowBit()) | 1 << (fTriggerBitConfig->GetJetHighBit()) ))
    return kTRUE;
  else
    return kFALSE;
}

/**
 * Check whether trigger is gamma patch trigger according to trigger bits
 * @param tBits Trigger bits of the fastor
 * @return True if fastor belongs to a gamma trigger, false otherwise
 */
Bool_t AliEmcalTriggerMaker::IsEGA(Int_t tBits) const {
  if( tBits & ( 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetTriggerTypesEnd() + fTriggerBitConfig->GetGammaHighBit()) | 1 << (fTriggerBitConfig->GetGammaLowBit()) | 1 << (fTriggerBitConfig->GetGammaHighBit()) ))
    return kTRUE;
  else
    return kFALSE;
}

/**
 * Check whether trigger is level0 trigger according to trigger bits
 * @param tBits Trigger bits of the fastor
 * @return True if fastor belongs to a level0 trigger, false otherwise
 */
Bool_t AliEmcalTriggerMaker::IsLevel0(Int_t tBits) const {
  if( tBits & (1 << (fTriggerBitConfig->GetLevel0Bit() + fTriggerBitConfig->GetLevel0Bit()) | (1 << fTriggerBitConfig->GetLevel0Bit())))
    return kTRUE;
  return kFALSE;
}

#endif
