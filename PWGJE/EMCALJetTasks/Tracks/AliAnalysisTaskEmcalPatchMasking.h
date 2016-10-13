#ifndef ALIANALYSISTASKEMCALPATCHMASKING_H
#define ALIANALYSISTASKEMCALPATCHMASKING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliEMCALTriggerDataGrid.h"
#include <TString.h>
#include <vector>

class AliOADBContainer;
class THistManager;

namespace EMCalTriggerPtAnalysis {

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
  UShort_t GetADC(UChar_t col, UChar_t row) const;

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
  void SetADC(UShort_t ADC, UChar_t col, UChar_t row);

private:
  UChar_t                               fPatchSize;       ///< Size of the patch
  AliEMCALTriggerDataGrid<UShort_t>     fADCValues;       ///< underlying container with ADC data

  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerPatchADCInfo, 1)
  /// \endcond
};

/**
 * @class AliAnalysisTaskEmcalPatchMasking
 * @brief Investigation of observables related to masked fastors within trigger patches
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Oct 13, 2016
 */
class AliAnalysisTaskEmcalPatchMasking : public AliAnalysisTaskEmcal {
public:

  /**
   * Dummy constructor
   */
  AliAnalysisTaskEmcalPatchMasking();

  /**
   * Constructor specifying the name of the task. Also prepares
   * the output container, steered by AliAnalysisTaskEmcal
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalPatchMasking(const char *name);

  /**
   * Destructor, cleaning up stuff belonging to this task
   */
  virtual ~AliAnalysisTaskEmcalPatchMasking();

  /**
   * Define trigger class the task runs over, consisting of a trigger bit and (optionally)
   * a pattern in the trigger string
   * @param[in] triggerbits
   * @param[in] triggerstring
   */
  void SetRequireTrigger(UInt_t triggerbits, const TString & triggerstring = "") {
    fTriggerBits = triggerbits; fTriggerPattern = triggerstring;
  }

  /**
   * Define OADB container with masked fastors
   * @param[in] oadbcont Name of the OADB container with the masked fastor abs IDs
   */
  void SetMaskedFastorOADB(const TString &oadbcont) { fNameMaskedFastorOADB = oadbcont; }

protected:

  virtual void UserCreateOutputObjects();
  virtual bool IsEventSelected();
  virtual bool Run();

  /**
   * Implementation of framework method Run: Initialize OADB handler for
   * masked fastors in case a masked fastor OADB container name is provided.
   */
  virtual void ExecOnce();

  /**
   * Implementation of framework method RunChanged: Loading of new
   * fastor bad channel map for the given run
   * @param[in] newrun new run number
   */
  virtual void RunChanged(int newrun);

  /**
   * Helper function filling interal ADC data grid.
   */
  void PrepareL1FastorADC();

  /**
   * Create Fastor-by-fastor ADC table for the input patch. The indexing
   * of the position is in col-row space, relative to the starting position
   * of the patch.
   * @param[in] patch Trigger patch for which ADC values are checked
   * @return Table structure with ADC values
   */
  AliEMCALTriggerPatchADCInfo *MakeFastorADCValuesForPatch(const AliEMCALTriggerPatchInfo &patch) const;

  /**
   * Perform analysis for given trigger patch: Collect single fastor
   * information and fill histograms
   * @param[in] patch Patch to process
   */
  void ProcessPatch(const AliEMCALTriggerPatchInfo &patch);

  /**
   * Perform analysis similar to ProcessPatch, but only for the max patch. Histograms are
   * different
   * @param[in] patch Patch to process
   * @param[in] maxtype String defining which type of max patch
   */
  void ProcessMaxPatch(const AliEMCALTriggerPatchInfo &patch, const TString &maxtype);

  THistManager                                *fHistos;                 //!<! Histogram manager
  AliEMCALTriggerDataGrid<int>                fL1ADC;                   ///< L1 ADC values
  ULong_t                                     fTriggerBits;             ///< Trigger bit selection
  TString                                     fTriggerPattern;          ///< Trigger pattern

  TString                                     fNameMaskedFastorOADB;    ///< Name of the masked fastor OADB container
  AliOADBContainer                            *fMaskedFastorOADB;       //!<! OADB container of the masked fastors
  std::vector<UShort_t>                       fListMaskedFastors;       ///< List of masked fastors

private:
  AliAnalysisTaskEmcalPatchMasking(const AliAnalysisTaskEmcalPatchMasking &);
  AliAnalysisTaskEmcalPatchMasking &operator=(const AliAnalysisTaskEmcalPatchMasking &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalPatchMasking, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALPATCHMASKING_H */
