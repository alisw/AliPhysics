/**
 * \file AliEmcalTriggerSetupInfo.h
 * \brief Manager for constants used in the trigger maker.
 *
 * \author Jiri Kral <>, University of Jyv&aumlskul&auml
 * \date Jun 26, 2013
 */
#ifndef ALIEMCALTRIGGERSETUPINFO_H
#define ALIEMCALTRIGGERSETUPINFO_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "AliEMCALTriggerConstants.h"

static const Double_t kEMCL1ADCtoADCSum = 3.40;             ///<

/**
 * \class AliEmcalTriggerSetupInfo
 * \brief Settings manager for the trigger patch algorithm
 * \ingroup EMCALTGRFW
 *
 * This class contains the main settings (trigger threshold) for the different EMCAL
 * Level1 trigger classes.
 */
class AliEmcalTriggerSetupInfo: public TNamed {
 public:
  AliEmcalTriggerSetupInfo();
  AliEmcalTriggerSetupInfo(const AliEmcalTriggerSetupInfo &p); 
  AliEmcalTriggerSetupInfo &operator=(const AliEmcalTriggerSetupInfo &p);
  virtual ~AliEmcalTriggerSetupInfo();

  /**
   * Get lower trigger threshold of the jet trigger (online)
   * @return Trigger threshold
   */
  Int_t GetThresholdJetLow() const { return fThresholds[2]; }
  /**
   * Get higher trigger threshold of the jet trigger (online)
   * @return Trigger threshold
   */
  Int_t GetThresholdJetHigh() const { return fThresholds[0]; }
  /**
   * Get lower trigger threshold of the jet trigger (offline)
   * @return Trigger threshold
   */
  Int_t GetThresholdJetLowSimple() const { return fThresholdsSimple[2]; }
  /**
   * Get higher trigger threshold of the jet trigger (offline)
   * @return Trigger threshold
   */
  Int_t GetThresholdJetHighSimple() const { return fThresholdsSimple[0]; }
  
  /**
   * Get lower trigger threshold of the gamma trigger (online)
   * @return Trigger threshold
   */
  Int_t GetThresholdGammaLow() const { return fThresholds[3]; }
  /**
   * Get higher trigger threshold of the gamma trigger (online)
   * @return Trigger threshold
   */
  Int_t GetThresholdGammaHigh() const { return fThresholds[1]; }
  /**
   * Get lower trigger threshold of the gamma trigger (offline)
   * @return Trigger threshold
   */
  Int_t GetThresholdGammaLowSimple() const { return fThresholdsSimple[3]; }
  /**
   * Get higher trigger threshold of the gamma trigger (offline)
   * @return Trigger threshold
   */
  Int_t GetThresholdGammaHighSimple() const { return fThresholdsSimple[1]; }

  /**
   * Get lower online trigger threshold of the jet trigger, converted to energy
   * @return Trigger threshold, converted to energy
   */
  Double_t GetThresholdGeVRoughJetLow() const { return ((Double_t)fThresholds[2])*EMCALTrigger::kEMCL1ADCtoGeV; }
  /**
   * Get higher online trigger threshold of the jet trigger, converted to energy
   * @return Trigger threshold, converted to energy
   */
  Double_t GetThresholdGeVRoughJetHigh() const { return ((Double_t)fThresholds[0])*EMCALTrigger::kEMCL1ADCtoGeV; }
  /**
   * Get lower offline trigger threshold of the jet trigger, converted to energy
   * @return Trigger threshold, converted to energy
   */
  Double_t GetThresholdGeVRoughJetLowSimple() const { return ((Double_t)fThresholdsSimple[2])*EMCALTrigger::kEMCL1ADCtoGeV; }
  /**
   * Get higher offline trigger threshold of the jet trigger, converted to energy
   * @return Trigger threshold, converted to energy
   */
  Double_t GetThresholdGeVRoughJetHighSimple() const { return ((Double_t)fThresholdsSimple[0])*EMCALTrigger::kEMCL1ADCtoGeV; }
  
  /**
   * Set trigger thresholds for the online trigger for the different trigger classes
   * @param i0 Jet trigger, high threshold
   * @param i1 Gamma trigger, high threshold
   * @param i2 Jet trigger, low threshold
   * @param i3 Gamma trigger, low threshold
   */
  void SetThresholds( Int_t i0, Int_t i1, Int_t i2, Int_t i3 ) {
            fThresholds[0] = i0; fThresholds[1] = i1; fThresholds[2] = i2; fThresholds[3] = i3;}
  /**
   * Set trigger thresholds for the simple offline trigger for the different trigger classes
   * @param i0 Jet trigger, high threshold
   * @param i1 Gamma trigger, high threshold
   * @param i2 Jet trigger, low threshold
   * @param i3 Gamma trigger, low threshold
   */
  void SetThresholdsSimple( Int_t i0, Int_t i1, Int_t i2, Int_t i3 ) {
            fThresholdsSimple[0] = i0; fThresholdsSimple[1] = i1; fThresholdsSimple[2] = i2; fThresholdsSimple[3] = i3;}
            
  void Clean();


 protected:
  Int_t             fThresholds[4];                 ///< per event L1 online thresholds in ADC counts
  Int_t             fThresholdsSimple[4];           ///< per event L1 simple offline thresholds

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerSetupInfo, 2) // Emcal trigger setup class
  /// \endcond
};
#endif
