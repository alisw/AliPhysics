/**
 * \file AliEMCalTriggerPatchAnalysisComponent.h
 * \brief Analysis component for EMCAL trigger patches
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>
#include "AliEMCalTriggerTracksAnalysisComponent.h"

class AliEMCALTriggerPatchInfo;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerPatchAnalysisComponent
 * \brief Analysis component for EMCAL trigger patches
 *
 * Analysis components for trigger patches. Fills THnSparses with the different energy definitions
 * (amplitude, estimated (rough) patch energy, calibrated (offline) energy and the patch position on
 * the EMCAL surface for EMCAL trigger patches of different categories created by the EMCAL trigger
 * patch maker. Main patches defined as the patches of a given type with the highest energy.
 */
class AliEMCalTriggerPatchAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerPatchAnalysisComponent();
  AliEMCalTriggerPatchAnalysisComponent(const char *name, Bool_t withEventSelection = kFALSE);
  virtual ~AliEMCalTriggerPatchAnalysisComponent() { }

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  /**
   * Check whether thresholds for online patches are swapped
   *
   * \return true if thresholds are swapped, fales otherwise
   */
  Bool_t IsSwapOnlineThresholds() const { return fSwapOnlineThresholds; }

  /**
   * Check whether thresholds for offline patches are swapped
   *
   * \return true if thresholds are swapped, fales otherwise
   */
  Bool_t IsSwapOfflineThresholds() const { return fSwapOfflineThresholds; }

  /**
   * Swap online thresholds (i.e. low threshold becomes high threshold and vice versa)
   *
   * \param doSwap Do swap the thresholds
   */
  void SetSwapOnlineThresholds(Bool_t doSwap = kTRUE) { fSwapOnlineThresholds = doSwap; }

  /**
   * Swap offline thresholds (i.e. low threshold becomes high threshold and vice versa)
   *
   * \param doSwap Do swap the thresholds
   */
  void SetSwapOfflineThresholds(Bool_t doSwap = kTRUE) { fSwapOfflineThresholds = doSwap; }

protected:

  /**
   * \class AliEmcalTriggerPatchHandlerFactory
   * \brief Factory class handling the identification of trigger patches as a given type
   *
   * This class handles the identification as trigger patch of a given type, neglecting whether the patch is an
   * online or an offline patch. The actual identification method is delegated to handler classes for each trigger
   * patch type inheriting from AliEmcalTriggerPatchHandler.
   */
  class AliEmcalTriggerPatchHandlerFactory{
  public:
    /**
     * Constructor
     * \param swapThresholdsOnline If true we swap thresholds for online patches
     * \param swapThresholdsOffline If true we swap thresholds for offline patches
     */
    AliEmcalTriggerPatchHandlerFactory(Bool_t swapThresholdsOnline, Bool_t swapThresholdsOffline):
      fSwapThresholdsOnline(swapThresholdsOnline),
      fSwapThresholdsOffline(swapThresholdsOffline)
    {}
    /**
     * Destructor
     */
    virtual ~AliEmcalTriggerPatchHandlerFactory() {}
    Bool_t IsPatchOfType(const AliEMCALTriggerPatchInfo *const patch, TString patchtype) const;

  protected:

    /**
     * \class AliEmcalTriggerPatchHandler
     * \brief Base class for trigger patch handlers in the trigger patch analysis component
     *
     * This class is the base class of patch handler classes, which identify given trigger patches
     * based on the trigger class.
     */
    class AliEmcalTriggerPatchHandler{
    public:
      /**
       * Constructor
       */
      AliEmcalTriggerPatchHandler(Bool_t doSwapOnline, Bool_t doSwapOffline):
        fPatchSwapThresholdsOnline(doSwapOnline),
        fPatchSwapThresholdsOffline(doSwapOffline)
      {}
      /**
       * Destructor
       */
      virtual ~AliEmcalTriggerPatchHandler() {}

      virtual Bool_t IsOfType(const AliEMCALTriggerPatchInfo * const patch) const = 0;
    protected:
      Bool_t                      fPatchSwapThresholdsOnline;          ///< Swap thresholds for online patches
      Bool_t                      fPatchSwapThresholdsOffline;         ///< Swap thresholds for offline patches
    };

    /**
     * \class AliEmcalTriggerPatchHandlerJetLow
     * \brief Patch handler implementation of the jet low trigger
     *
     * Helper class implementing the handling of the jet low trigger in the patch analysis component.
     */
    class AliEmcalTriggerPatchHandlerJetLow : public AliEmcalTriggerPatchHandler{
    public:
      /**
       * Constructor
       */
      AliEmcalTriggerPatchHandlerJetLow(Bool_t doSwapOnline, Bool_t doSwapOffline) :
        AliEmcalTriggerPatchHandler(doSwapOnline, doSwapOffline)
      {}
      /**
       * Destructor
       */
      virtual ~AliEmcalTriggerPatchHandlerJetLow() {}
      virtual Bool_t IsOfType(const AliEMCALTriggerPatchInfo * const patch) const;
    };

    /**
     * \class AliEmcalTriggerPatchHandlerJetHigh
     * \brief Patch handler implementation of the jet high trigger
     *
     * Helper class implementing the handling of the jet high trigger in the patch analysis component.
     */
    class AliEmcalTriggerPatchHandlerJetHigh : public AliEmcalTriggerPatchHandler{
    public:
      /**
       * Constructor
       */
      AliEmcalTriggerPatchHandlerJetHigh(Bool_t doSwapOnline, Bool_t doSwapOffline) :
        AliEmcalTriggerPatchHandler(doSwapOnline, doSwapOffline)
      {}
      /**
       * Destructor
       */
      virtual ~AliEmcalTriggerPatchHandlerJetHigh() {}
      virtual Bool_t IsOfType(const AliEMCALTriggerPatchInfo * const patch) const;
    };

    /**
     * \class AliEmcalTriggerPatchHandlerGammaLow
     * \brief Patch handler implementation of the gamma low trigger
     *
     * Helper class implementing the handling of the gamma low trigger in the patch analysis component.
     */
    class AliEmcalTriggerPatchHandlerGammaLow : public AliEmcalTriggerPatchHandler{
    public:
      /**
       * Constructor
       */
      AliEmcalTriggerPatchHandlerGammaLow(Bool_t doSwapOnline, Bool_t doSwapOffline) :
        AliEmcalTriggerPatchHandler(doSwapOnline, doSwapOffline)
      {}
      /**
       * Destructor
       */
      virtual ~AliEmcalTriggerPatchHandlerGammaLow() {}
      virtual Bool_t IsOfType(const AliEMCALTriggerPatchInfo * const patch) const;
    };

    /**
     * \class AliEmcalTriggerPatchHandlerGammaHigh
     * \brief Patch handler implementation of the gamma high trigger
     *
     * Helper class implementing the handling of the gamma high trigger in the patch analysis component.
     */
    class AliEmcalTriggerPatchHandlerGammaHigh : public AliEmcalTriggerPatchHandler{
    public:
      /**
       * Constructor
       */
      AliEmcalTriggerPatchHandlerGammaHigh(Bool_t doSwapOnline, Bool_t doSwapOffline) :
        AliEmcalTriggerPatchHandler(doSwapOnline, doSwapOffline)
      {}
      /**
       * Destructor
       */
      virtual ~AliEmcalTriggerPatchHandlerGammaHigh() {}
      virtual Bool_t IsOfType(const AliEMCALTriggerPatchInfo * const patch) const;
    };

    /**
     * \class AliEmcalTriggerPatchHandlerLevel0
     * \brief Patch handler implementation of the level0 trigger
     *
     * Helper class implementing the handling of the level0 trigger in the patch analysis component.
     */
    class AliEmcalTriggerPatchHandlerLevel0 : public AliEmcalTriggerPatchHandler{
    public:
      /**
       * Constructor
       */
      AliEmcalTriggerPatchHandlerLevel0(Bool_t doSwapOnline, Bool_t doSwapOffline) :
        AliEmcalTriggerPatchHandler(doSwapOnline, doSwapOffline)
      {}
      /**
       * Destructor
       */
      virtual ~AliEmcalTriggerPatchHandlerLevel0() {}
      virtual Bool_t IsOfType(const AliEMCALTriggerPatchInfo * const patch) const;
    };

    Bool_t                      fSwapThresholdsOnline;          ///< Swap thresholds for online patches
    Bool_t                      fSwapThresholdsOffline;         ///< Swap thresholds for offline patches
  };

  void FillStandardMonitoring(const AliEMCALTriggerPatchInfo * const patch, TString eventclass = "");
  void FillTriggerInfoHistogram(TString histo, const AliEMCALTriggerPatchInfo *const patch);
  void FillAmplitudeHistogram(TString histo, const AliEMCALTriggerPatchInfo *const patch);

  Bool_t                        fSwapOnlineThresholds;          ///< Swap trigger thresholds for online patches
  Bool_t                        fSwapOfflineThresholds;         ///< Swap trigger thresholds for offline patches
  Bool_t                        fWithEventSelection;            ///< Define whether patches are analysed with event selection

  ClassDef(AliEMCalTriggerPatchAnalysisComponent, 1);     // Component for trigger patch analysis
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H */
