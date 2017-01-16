#ifndef ALIEMCALTRIGGERPATCHINFOV1_H
#define ALIEMCALTRIGGERPATCHINFOV1_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCALTriggerPatchInfo.h"

/**
 * @class AliEMCALTriggerPatchInfoV1
 * @brief Temporary class handling
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Jan. 16, 2017
 *
 * This class temporarily extends AliEMCALTriggerPatchInfo with
 * a field for the smeared energy until the additional data member
 * in the additional class becomes available for AliPhysics
 */
class AliEMCALTriggerPatchInfoV1 : public AliEMCALTriggerPatchInfo {
public:

  /**
   * Constructor
   */
  AliEMCALTriggerPatchInfoV1();

  /**
   * Destructor
   */
  virtual ~AliEMCALTriggerPatchInfoV1() {}

  /**
   * Set the smeared energy
   * @param e Smeared energy
   */
  void SetSmearedEnergyV1(double e) { fSmearedEnergyV1 = e;}

  /**
   * Get the smeared energy
   * @return Smeared energy
   */
  Double_t GetSmearedEnergyV1() const { return fSmearedEnergyV1; }

  Double_t GetSmearedETV1() const;


  /**
   * Allocate a new AliEMCALTriggerPatchInfo object and initialize it
   * @param col0        Start column of the patch
   * @param row0        Start row of the patch
   * @param size        Size of the patch
   * @param adc         ADC signal of the patch
   * @param offlineAdc  Offline ADC signal of the patch
   * @param patchE      Energy of the patch (sum of cell amplitudes)
   * @param bitmask     Trigger bit mask of the patch
   * @param vertex      Primary vertex of the event
   * @param geom        Pointer to the EMCal geometry object
   * @return            Pointer to a new and initialized AliEMCALTriggerPatchInfo object (caller is responsible for releasing memory)
   */
  static AliEMCALTriggerPatchInfoV1* CreateAndInitializeV1(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom);

  /**
   * Allocate a new AliEMCALTriggerPatchInfo object and initialize it
   * @param col0        Start column of the patch
   * @param row0        Start row of the patch
   * @param size        Size of the patch
   * @param adc         ADC signal of the patch
   * @param offlineAdc  Offline ADC signal of the patch
   * @param bitmask     Trigger bit mask of the patch
   * @param geom        Pointer to the EMCal geometry object
   * @return            Pointer to a new and initialized AliEMCALTriggerPatchInfo object (caller is responsible for releasing memory)
   */
  static AliEMCALTriggerPatchInfoV1* CreateAndInitializeV1(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, UInt_t bitmask, const AliEMCALGeometry* geom);


private:

  Double_t                fSmearedEnergyV1;             ///< Smeared energy

  ClassDef(AliEMCALTriggerPatchInfoV1, 1)
};

#endif /* ALIEMCALTRIGGERPATCHINFOV1_H */
