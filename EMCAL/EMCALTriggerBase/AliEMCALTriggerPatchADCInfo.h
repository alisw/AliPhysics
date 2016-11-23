#ifndef ALIEMCALTRIGGERPATCHADCINFO_H
#define ALIEMCALTRIGGERPATCHADCINFO_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliEMCALTriggerDataGrid.h"

/**
 * @class AliEMCALTriggerPatchADCInfo
 * @brief Helper class containing all ADC infos for a given trigger patch
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Oct 13, 2016
 *
 * This class is a simple helper class providing access to the single FastOR L1
 * ADC values within a trigger elementary patch. Single ADC values can be obtained
 * via the position in the patch using column and row.
 *
 * ~~~{.cxx}
 * AliEMCALTriggerPatchInfo recpatch;
 * AliEMCALTriggerPatchADCInfo patch
 * for(int icol = 0; icol < recpatch->GetPatchSize(); icol++){
 *   for(int irow = 0; irow < recpatch->GetPatchSize(); irow++){
 *     std::cout << "ADC[" << icol << "," << irow << "] = " << patch.GetADC(icol, irow) << std::endl;
 *   }
 * }
 * ~~~
 */
class AliEMCALTriggerPatchADCInfo : public TObject {
public:

  /**
   * Dummy constructor
   */
  AliEMCALTriggerPatchADCInfo();

  /**
   * Main constructor, initializing patch size and data container
   * for the ADC values
   * @param[in] patchsize Size of the trigger patch
   */
  AliEMCALTriggerPatchADCInfo(UChar_t patchsize);

  /**
   * Copy constructor
   * @param [in]ref Reference for the copy
   */
  AliEMCALTriggerPatchADCInfo(const AliEMCALTriggerPatchADCInfo &ref);

  /**
   * Assignment operator
   * @param[in] ref Reference for the assignment
   * @return
   */
  AliEMCALTriggerPatchADCInfo &operator=(const AliEMCALTriggerPatchADCInfo &ref);

  /**
   * Destructor, deleting underlying data container
   */
  virtual ~AliEMCALTriggerPatchADCInfo() {}

  /**
   * Get the ADC value at the given position in the patch. The
   * position (col, row) is supposed to be relative to the starting
   * position of the trigger patch.
   * @param[in] col Coloumn of the fastor postion relative to the patch starting position
   * @param[in] row Row of the fastor position relative to the patch starting postition
   * @return Fastor ADC value at the given position of the patch
   */
  Int_t GetADC(UChar_t col, UChar_t row) const;

  /**
   * Calculate sum of all FastOR ADC values connected
   * @return Sum of all FastOR ADC values
   */
  Int_t GetSumADC() const;

  /**
   * Get the maximum FastOR ADC value within the trigger
   * patch
   * @return Maximum FastOR ADC value in trigger patch
   */
  Int_t GetMaxADC() const;

  /**
   * Get number of non-zero FastOR ADC values within a trigger
   * patch (equivalent to contributing fastors).
   * @return Number of non-zero (contributing) FastORs
   */
  Int_t GetNFastorsContrib() const;

  /**
   * Get the size of the patch
   * @return Size of the patch
   */
  UChar_t GetPatchSize() const { return fPatchSize; }

  /**
   * Set the size of the patch. Attention: Allocates space for
   * the ADC values
   * @param[in] patchsize Size of the patch
   */
  void SetPatchSize(UChar_t patchsize);

  /**
   * Set ADC value for a given fastor position inside the trigger
   * patch. The position (col, row) is supposed to be relative to
   * the starting position of the trigger patch.
   * @param ADC Fastor ADC value for the given fastor position
   * @param col Coloumn of the fastor postion relative to the patch starting position
   * @param row Row of the fastor position relative to the patch starting postition
   */
  void SetADC(Int_t adc, UChar_t col, UChar_t row);

private:
  UChar_t                               fPatchSize;       ///< Size of the patch
  AliEMCALTriggerDataGrid<Int_t>        fADCValues;       ///< underlying container with ADC data

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerPatchADCInfo, 1)
  /// \endcond
};



#endif /* PWG_EMCAL_EMCALTRIGGER_ALIEMCALTRIGGERPATCHADCINFO_H_ */
