#ifndef ALIEMCALCORRECTIONEVENTMANAGER_H
#define ALIEMCALCORRECTIONEVENTMANAGER_H

#include "AliAnalysisTaskSE.h"

/**
 * @class AliEmcalCorrectionEventManager
 * @ingroup EMCALCOREFW
 * @brief Class to manage the current event for correction components in the EMCal correction framework
 *
 * Manages the input event given to each correction component. The manager is given access to both the
 * internal and external (embedded) event and then returns the configured event when InputEvent() is called.
 * This allows the event manager to replace any class to AliAnalysisTaskSE::InputEvent() seamlessly.
 *
 * This is particularly useful for a component which only uses (for example) embedded input objects, as it
 * can be supplied the embedded event instead of the internal event. In such a case, the component doesn't
 * need to know anything about the embedding framework. This also allows the integration of the PHOS tender
 * into the EMCal Correction Framework.
 *
 * Note that this class defaults to using the **embedded event** so that any call to a properly designed
 * function will cause the internal event to be used. For an example of such usage, see
 * AliEmcalCorrectionTask::AddContainersToComponent() and the corresponding calls to
 * AliEmcalCorrectionComponent::SetUsingInputEvent(). In that case, if there is a non-embedded container,
 * it will call SetUsingInputEvent(true), which will ensure that the component uses the internal event.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @date Jun 29, 2017
 */

class AliEmcalCorrectionEventManager : public AliAnalysisTaskSE
{
 public:
  AliEmcalCorrectionEventManager();
  AliEmcalCorrectionEventManager(const char * name);
  virtual ~AliEmcalCorrectionEventManager();

  AliVEvent * InputEvent() const;
  /// True if the embedding event should be used.
  bool UseEmbeddingEvent() const { return fUseEmbedding; }

  /**
   * Set the input event of the current event. This event will be returned by
   * InputEvent() unless fUseEmbedding is true.
   *
   * @param inputEvent Input event of the current event
   */
  void SetInputEvent(AliVEvent * inputEvent) { fInputEvent = inputEvent; }

  /**
   * Determines whether the embedded event is returned when InputEvent() is called.
   *
   * @param b If true, return the embedded event when InputEvent() is called.
   */
  void SetUseEmbeddingEvent(bool b = true) { fUseEmbedding = b; }

 protected:
  bool fUseEmbedding;      ///< If true, return the embedded event instead

 private:
  AliEmcalCorrectionEventManager(const AliEmcalCorrectionEventManager &);               // Not implemented
  AliEmcalCorrectionEventManager &operator=(const AliEmcalCorrectionEventManager &);    // Not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionEventManager, 1); // AliEmcalCorrectionEventManager
  /// \endcond
};

#endif /* AliEmcalCorrectionEventManager.h */
