#ifndef ALIANALYSISTASKEMCALPATCHMASKING_H
#define ALIANALYSISTASKEMCALPATCHMASKING_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "AliEMCALTriggerDataGrid.h"
#include <TString.h>
#include <vector>

class AliEMCALTriggerPatchADCInfoAP;
class AliOADBContainer;
class THistManager;

namespace PWGJE {
  
namespace EMCALJetTasks {

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
  AliEMCALTriggerPatchADCInfoAP *MakeFastorADCValuesForPatch(const AliEMCALTriggerPatchInfo &patch) const;

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
  AliEMCALTriggerDataGrid<Int_t>              fL1ADC;                   ///< L1 ADC values
  ULong_t                                     fTriggerBits;             ///< Trigger bit selection
  TString                                     fTriggerPattern;          ///< Trigger pattern

  TString                                     fNameMaskedFastorOADB;    ///< Name of the masked fastor OADB container
  AliOADBContainer                            *fMaskedFastorOADB;       //!<! OADB container of the masked fastors
  std::vector<UShort_t>                       fListMaskedFastors;       ///< List of masked fastors

private:
  AliAnalysisTaskEmcalPatchMasking(const AliAnalysisTaskEmcalPatchMasking &);
  AliAnalysisTaskEmcalPatchMasking &operator=(const AliAnalysisTaskEmcalPatchMasking &);

  ClassDef(AliAnalysisTaskEmcalPatchMasking, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALPATCHMASKING_H */
