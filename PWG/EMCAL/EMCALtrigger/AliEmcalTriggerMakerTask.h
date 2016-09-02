#ifndef ALIEMCALTRIGGERMAKERTASK_H
#define ALIEMCALTRIGGERMAKERTASK_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEmcalTriggerMakerKernel.h"
#include "AliAnalysisTaskEmcal.h"
#include <TString.h>
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <functional>
#endif

class TClonesArray;
class THistManager;
class AliVVZERO;
class AliEMCALTriggerPatchInfo;

/**
 * @class AliEmcalTriggerMakerTask
 * @brief EMCAL trigger maker task
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @author Salvatore Aiola, Yale University
 * @since Oct. 22nd, 2015
 * @ingroup EMCALTRGFW
 *
 * The EMCAL trigger maker task steers the process building
 * trigger patches, performed via the trigger maker kernel,
 * and provides the interface to the user. As this, it forwards
 * necessary data like FASTor amplitudes and time sums as well
 * as cell data to the trigger maker kernel and reads out the
 * patches. The patches as type of raw patches are converted
 * into full patches (AliEMCALTriggerPatchInfo) and stored
 * in a TClonesArray which is added to the input event. On user
 * request it can also fill QA histograms which are added to
 * the common root file.
 */
class AliEmcalTriggerMakerTask : public AliAnalysisTaskEmcal {
public:
  /***
   * @enum TriggerMakerTriggerBitConfig_t
   * @brief Definition of trigger bit configurations
   *
   * This enumeration handles different trigger bit configurations for the
   * EMCAL Level1 triggers (with and without different thresholds) applied
   * in the reconstruction of different samples.
   */
  enum TriggerMakerBitConfig_t {
    kOldConfig = 0,///< Old configuration, no distinction between high and low threshold
    kNewConfig = 1 ///< New configuration, distinction between high and low threshold
  };

  /**
   * @enum Exception_t
   * @brief Definition of various exception codes used in the trigger maker
   *
   * Collection of all possible error codes issued in the trigger maker task
   * class.
   */
  enum Exception_t {
    kInvalidChannelException = 1      ///< kInvalidChannelException
  };

  /**
   * Dummy constructor
   */
  AliEmcalTriggerMakerTask();

  /**
   * Main constructor, initializing also AliAnalysisTaskEmcal and Create
   * the Kernel.
   * @param name Name of the task
   * @param doQA If true we switch on QA
   */
  AliEmcalTriggerMakerTask(const char *name, Bool_t doQA = kFALSE);

  /**
   * Destructor
   */
  virtual ~AliEmcalTriggerMakerTask();

  /**
   * Initializing output objects. Kernel is initialized in
   * ExecOnce as geometry - which depends on the run number -
   * is needed for this.
   */
  virtual void      UserCreateOutputObjects();

  /**
   * Initialize the trigger maker kernel. This function is called
   * for the first event in order to obtain the run number from it.
   * Also initialises AliAnalysisTaskEmcal.
   */
  virtual void      ExecOnce();

  /**
   * Run trigger patch finding.
   *
   * As the patch finding is implemented in the class AliEmcalTriggerMakerKernel,
   * the analysis task only takes care about propagating necessary information to
   * the kernel (ADCs, cell energy, V0 amplitudes), steering the patch finder, and
   * appending the output to the ESD event.
   *
   * @return Always true.
   */
  virtual Bool_t    Run();

  /**
   * Switch on basic QA of the trigger maker
   * @param[in] doQA If true QA is switched on.
   */
  void SetRunQA(Bool_t doQA = kTRUE) { fDoQA = doQA; }

  /**
   * Set range for L0 time
   * @param[in] min Minimum L0 time (default is 7)
   * @param[in] max Maximum L0 time (default is 10)
   */
  void SetL0TimeRange(Int_t min, Int_t max) { if (fTriggerMaker) fTriggerMaker->SetL0TimeRange(min, max); }


  /**
   * Set the name of the output container
   * @param[in] name Name of the output container
   */
  void SetCaloTriggersOutName(const char *name)     { fCaloTriggersOutName      = name; }

  /**
   * Set the name of the V0 input object
   * @param[in] name Name of the V0 input object
   */
  void SetV0InName(const char *name) { fV0InName      = name; }

  /**
   * Set the name of the OADB container providing a list of absolute IDs of FEE channels to
   * be masked (TObjArray of TParameter<int>)
   * @param[in] name Name of the container providing masked channel abs IDs
   */
  void SetMaskedFEEChannelOADBContainer(const TString &name) { fBadFEEChannelOADB = name; }

  /**
   * Trigger bit configuration to be used in the trigger patch maker.
   * @param[in] bitConfig Type of the trigger bit config (old - 3 bit, new - 5 bit)
   */
  void SetUseTriggerBitConfig(TriggerMakerBitConfig_t bitConfig);

  /**
   * Defining whether to use the FastOR masking from the OCDB. Attention:
   * In case the option is set to true the CDBconnect task is required to
   * run before the trigger maker task.
   * @param[in] doLoad If true the FastOR masking is obtained from the OCDB
   */
  void SetLoadFastORMaskingFromOCDB(Bool_t doLoad = kTRUE) { fLoadFastORMaskingFromOCDB = doLoad; }

  /**
   * Apply online bad channels (FastOR masking) to offline FEE energies. Using this
   * the same acceptance which was available online can be applied to offline patches.
   * The masking is handled inside the trigger maker kernel.
   * @param[in] doApply If true the online FastOR masking is applied to offline energies
   */
  void SetApplyTRUMaskingToFEE(Bool_t doApply = kTRUE) {
    if(fTriggerMaker) fTriggerMaker->SetApplyOnlineBadChannelMaskingToOffline(doApply);
  }

  void SetTriggerThresholdJetLow   ( Int_t a, Int_t b, Int_t c ) {
    if(fTriggerMaker) fTriggerMaker->SetTriggerThresholdJetLow(a, b, c);
  }
  void SetTriggerThresholdJetHigh  ( Int_t a, Int_t b, Int_t c ) {
    if(fTriggerMaker) fTriggerMaker->SetTriggerThresholdJetHigh(a, b, c);
  }
  void SetTriggerThresholdGammaLow ( Int_t a, Int_t b, Int_t c ) {
    if(fTriggerMaker) fTriggerMaker->SetTriggerThresholdGammaLow(a, b, c);
  }
  void SetTriggerThresholdGammaHigh( Int_t a, Int_t b, Int_t c ) {
    if(fTriggerMaker) fTriggerMaker->SetTriggerThresholdGammaHigh(a, b, c);
  }

  AliEmcalTriggerMakerKernel* GetTriggerMaker() const { return fTriggerMaker; }

  void SetUseL0Amplitudes(Bool_t b)        { fUseL0Amplitudes = b           ; }

protected:

#if ! defined(__CINT__) || defined(__MAKECINT__)
  /**
   * Closure producing a handler converting a set of mask / bit number into a channel
   * ID. In case the mask / bit number is invalid the function will return -1
   * @return function that converts a set of mask / bit number into a channel ID
   */
  std::function<int (unsigned int, unsigned int)> GetMaskHandler() const;
#endif

  /**
   * Filling basic QA histograms implemented in the trigger maker task for
   * a given trigger patch which is of a certain patch type. Note that the
   * patchtype provided here defines the histograms to be filled (a trigger
   * patch can be of multiple types).
   * @param[in] patchtype Patchtype of the histograms to be filled
   * @param[in] recpatch Patch to be monitored
   */
  void FillQAHistos(const TString &patchtype, const AliEMCALTriggerPatchInfo &recpatch);

  /**
   * RunChanged method of the trigger maker:
   * In case of a bad channel container provided, initialize bad channel map
   * for the new run.
   */
  virtual void RunChanged();

  /**
   * Initialize bad channels from the OADB container. Called by the RunChanged() method.
   */
  void InitializeBadFEEChannels();

  /**
   * Initialize the FastOR masking from the OCDB
   */
  void InitializeFastORMaskingFromOCDB();

  AliEmcalTriggerMakerKernel              *fTriggerMaker;             ///< The actual trigger maker kernel
  AliVVZERO                               *fV0;                       //!<! VZERO data

  TString                                 fCaloTriggersOutName;       ///< name of output track array
  TString                                 fV0InName;                  ///< name of output track array
  TString                                 fBadFEEChannelOADB;         ///< name of the OADB container containing channels to be masked inside the trigger maker
  Bool_t                                  fUseL0Amplitudes;           ///< Use L0 amplitudes instead of L1 time sum (useful for runs where STU was not read)
  Bool_t                                  fLoadFastORMaskingFromOCDB; ///< Load FastOR masking from the OCDB
  TClonesArray                            *fCaloTriggersOut;          //!<! trigger array out

  Bool_t                                  fDoQA;                      ///< Fill QA histograms
  THistManager                            *fQAHistos;                 //!<! Histograms for QA

private:
  AliEmcalTriggerMakerTask(const AliEmcalTriggerMakerTask &);
  AliEmcalTriggerMakerTask &operator=(const AliEmcalTriggerMakerTask &);

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerMakerTask, 2);
  /// \endcond
};

#endif
