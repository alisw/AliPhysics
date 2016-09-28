#ifndef ALIEMCALFASTORMONITORTASK_H
#define ALIEMCALFASTORMONITORTASK_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include <TString.h>

class AliEMCALGeometry;
class THistManager;

/**
 * @class AliEmcalFastOrMonitorTask
 * @brief Simlple monitoring of EMCAL FastOr quantities
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Sept 8, 2016
 * @ingroup EMCALFWTASKS
 *
 * Simple monitoring task for FastOr related quantities. Filling the following
 * distributions:
 * - Frequency
 * - Amplitude
 * - L0timeSum
 * - Number of LO times
 * - Position in the col-row space
 *
 * The wagon can be added to the train using the corresponding AddMacro (AddEmcalFastOrMonitorTask):
 *
 *  ~~~{.cxx}
 *  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddEmcalFastOrMonitorTask.C");
 *  AddEmcalFastOrMonitorTask();
 *  ~~~
 */
class AliEmcalFastOrMonitorTask : public AliAnalysisTaskSE {
public:

  /**
   * Default constructor. For ROOT I/O
   */
  AliEmcalFastOrMonitorTask();

  /**
   * Named constructor
   * @param name Name of the task
   */
  AliEmcalFastOrMonitorTask(const char *name);

  /**
   * Destructor
   */
  virtual ~AliEmcalFastOrMonitorTask();

  /**
   * Define trigger selection. It can consist of trigger bits and a trigger string
   * @param triggerbits Trigger bit selection
   * @param triggerstring Trigger string (i.e. EG1, EG2, DG1, ...)
   */
  void SetRequestTrigger(ULong_t triggerbits, TString triggerstring = "") {
    fRequestTrigger = triggerbits;
    fTriggerPattern = triggerstring;
  }

  /**
   * Add masked fastor to the list of masked fastors. Masked fastors will be
   * ignored in ADC spectrum.
   * @param fastorID Abs ID of the fastor to be masked
   */
  void AddMaskedFastor(int fastorID){ fMaskedFastors.push_back(fastorID); }

protected:

  /**
   * Creating output objects. In this case only the histogram handler is created.
   * Histograms are done inside the function ExecOnce.
   */
  virtual void UserCreateOutputObjects();

  /**
   * Event loop. Filling the monitoring histograms for each FastOR:
   * - Frequency
   * - Amplitude
   * - L0timeSum
   * - Number of LO times
   * - Position in the col-row space
   * @param
   */
  virtual void UserExec(Option_t *);

  /**
   * Performing initial initializations. In contrast to UserCreateOutputObjects,
   * which is called before the event loop, ExecOnce is called for the first event
   * within the event loop. At that step some basic event information is already
   * available,
   */
  virtual void ExecOnce();

  /**
   * Performing run-dependent initializations. This function is useful
   * i.e. to load parameters from the OCDB/OADB
   */
  virtual void RunChanged();

  THistManager                            *fHistos;           //!<! Histogram handler
  AliEMCALGeometry                        *fGeom;             //!<! EMCAL Geometry object
  Bool_t                                  fLocalInitialized;  ///< Switch whether task is initialized (for ExecOnce)
  Int_t                                   fOldRun;            ///< Old Run (for RunChanged())

  ULong_t                                 fRequestTrigger;    ///< Trigger selection bits
  TString                                 fTriggerPattern;    ///< Trigger string pattern used in addition to the trigger selection bits

  std::vector<int>                        fMaskedFastors;     ///< List of masked fastors

  /// \cond CLASSIMP
  ClassDef(AliEmcalFastOrMonitorTask, 1);
  /// \endcond
};

#endif /* ALIEMCALFASTORMONITORTASK_H */
