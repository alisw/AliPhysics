/**
 * \file AliEmcalTriggerBitConfig.h
 * \brief Definition of EMCAL trigger bit configurations
 *
 * Two configuation are defined:
 * - 2 triggers, 1 threshold
 * - 2 triggers, 2 thresholds
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGERBITCONFIG_H
#define ALIEMCALTRIGGERBITCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>

/**
 * @class AliEmcalTriggerBitConfig
 * @brief Definition of trigger bit configurations
 *
 * Trigger bit configuration used in the trigger patch maker and by the trigger patches
 * themselves in order to identify of which type the trigger patch is. Can be adapted to different
 * trigger bit configurations use in different reconstructions
 */
class AliEMCALTriggerBitConfig  : public TNamed {
public:

  /**
   * Dummy constructor for the configuraiton base classes, not to be callled
   */
  AliEMCALTriggerBitConfig();

  /**
   * Constructor initialising the configurations. Used by the inheriting classes
   * @param l0bit Trigger bit for Level0 triggers
   * @param jhighbit Trigger bit for jet trigger, high threshold
   * @param jlowbit Trigger bit for jet trigger, low threshold
   * @param ghighbit Trigger bit for gamma trigger, high threshold
   * @param glowbit Trigger bit for gamma trigger, low threshold
   * @param mcoffset Offset for MC
   */
  AliEMCALTriggerBitConfig(Int_t l0bit, Int_t j1bit, Int_t j2bit, Int_t g1bit, Int_t g2bit, Int_t mcoffset);

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerBitConfig() {}

  /**
   * Initialise from other object
   * @param ref Reference to initialize this trigger bit configuaration from
   */
  void Initialise(const AliEMCALTriggerBitConfig &ref);

  /**
   * Get the trigger bit for Level0 triggers
   * @return Trigger bit number
   */
  Int_t GetLevel0Bit() const { if(fL0Bit < 0) Fatal("AliEMCALTriggerBitConfig", "Invalid trigger configuration: Level0 bit < 0"); return fL0Bit; }

  /**
   * Get trigger bit for the jet trigger, high threshold
   * @return Trigger bit number
   */
  Int_t GetJetHighBit() const { if(fJHighBit < 0) Fatal("AliEMCALTriggerBitConfig", "Invalid trigger configuration: Jet high bit < 0"); return fJHighBit; }

  /**
   * Get trigger bit for the jet trigger, low threshold
   * @return Trigger bit number
   */
  Int_t GetJetLowBit() const { if(fJLowBit < 0) Fatal("AliEMCALTriggerBitConfig", "Invalid trigger configuration: Jet low bit < 0"); return fJLowBit; }

  /**
   * Get trigger bit for the gamma trigger, high treshold
   * @return Trigger bit number
   */
  Int_t GetGammaHighBit() const { if(fGHighBit < 0) Fatal("AliEMCALTriggerBitConfig", "Invalid trigger configuration: Gamma high bit < 0"); return fGHighBit; }

  /**
   * Get trigger bit for the gamma trigger, low threshold
   * @return Trigger bit number
   */
  Int_t GetGammaLowBit() const { if(fGLowBit < 0) Fatal("AliEMCALTriggerBitConfig", "Invalid trigger configuration: Gamma low bit < 0"); return fGLowBit; }

  /**
   * Get the Monte-Carlo offset
   * @return MC offset
   */
  Int_t GetTriggerTypesEnd() const {if(fTriggerTypesEnd < 0) Fatal("AliEMCALTriggerBitConfig", "Invalid trigger configuration: MC Offset bit < 0"); return fTriggerTypesEnd; }

protected:
 Int_t fL0Bit;                      ///< Level0 bit
 Int_t fJHighBit;                   ///< Jet High bit
 Int_t fJLowBit;                    ///< Jet Low bit
 Int_t fGHighBit;                   ///< Gamma High bit
 Int_t fGLowBit;                    ///< Gamma Low bit
 Int_t fTriggerTypesEnd;            ///< Monte-Carlo offset

 /// \cond CLASSIMP
 ClassDef(AliEMCALTriggerBitConfig, 1);
 /// \endcond
};

/**
 * @class AliEmcalTriggerBitConfigOld
 * @brief Definition of old trigger bit configuration
 *
 * Definition of old trigger bit configuration:
 * - 2 triggers, 1 threshold
 */
class AliEMCALTriggerBitConfigOld : public AliEMCALTriggerBitConfig{
public:

  /**
   * Constructor, initializing the configuration
   */
  AliEMCALTriggerBitConfigOld();

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerBitConfigOld() {}

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerBitConfigOld, 1);
  /// \endcond
};

/**
 * @class AliEmcalTriggerBitConfigNew
 * @brief Definition of new trigger bit configuration
 *
 * Definition of new trigger bit configuration:
 * - 2 triggers, 2 thresholds
 */
class AliEMCALTriggerBitConfigNew : public AliEMCALTriggerBitConfig{
public:

  /**
   * Constructor, initializing the configuration
   */
  AliEMCALTriggerBitConfigNew();

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerBitConfigNew() {}

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerBitConfigNew, 1);
  /// \endcond
};

#endif /* ALIEMCALTRIGGERBITCONFIG_H */
