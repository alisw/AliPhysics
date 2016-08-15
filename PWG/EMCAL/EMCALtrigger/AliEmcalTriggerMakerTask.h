#ifndef ALIEMCALTRIGGERMAKERTASK_H
#define ALIEMCALTRIGGERMAKERTASK_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEmcalTriggerMakerKernel.h"
#include "AliAnalysisTaskEmcal.h"
#include <TString.h>

class TClonesArray;
class THistManager;
class AliVVZERO;
class AliEMCALTriggerPatchInfo;

/**
 * \class AliEmcalTriggerMakerTask
 * \brief EMCAL trigger maker task
 * \ingroup EMCALTRGFW
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
   * @param doQA If true QA is switched on.
   */
  void SetRunQA(Bool_t doQA = kTRUE) { fDoQA = doQA; }

  /**
   * Set range for L0 time
   * @param min Minimum L0 time (default is 7)
   * @param max Maximum L0 time (default is 10)
   */
  void SetL0TimeRange(Int_t min, Int_t max) { if (fTriggerMaker) fTriggerMaker->SetL0TimeRange(min, max); }


  /**
   * Set the name of the output container
   * @param name Name of the output container
   */
  void SetCaloTriggersOutName(const char *name)     { fCaloTriggersOutName      = name; }

  /**
   * Set the name of the V0 input object
   * @param name Name of the V0 input object
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
   * @param bitConfig Type of the trigger bit config (old - 3 bit, new - 5 bit)
   */
  void SetUseTriggerBitConfig(TriggerMakerBitConfig_t bitConfig);

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

  AliEmcalTriggerMakerKernel              *fTriggerMaker;             ///< The actual trigger maker kernel
  AliVVZERO                               *fV0;                       //!<! VZERO data

  TString                                 fCaloTriggersOutName;       ///< name of output track array
  TString                                 fV0InName;                  ///< name of output track array
  TString                                 fBadFEEChannelOADB;         ///< name of the OADB container containing channels to be masked inside the trigger maker
  Bool_t                                  fUseL0Amplitudes;           ///< Use L0 amplitudes instead of L1 time sum (useful for runs where STU was not read)
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
