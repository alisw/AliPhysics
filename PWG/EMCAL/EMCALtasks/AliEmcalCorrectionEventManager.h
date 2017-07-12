#ifndef ALIEMCALCORRECTIONEVENTMANAGER_H
#define ALIEMCALCORRECTIONEVENTMANAGER_H

#include "AliAnalysisTaskSE.h"

/**
 * @class AliEmcalCorrectionEventManager
 * @ingroup EMCALCOREFW
 * @brief Class to manage the current event for correction components in the EMCal correction framework
 *
 * TODO: Additional description
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
  /// Whether the embedding should be used.
  bool UseEmbeddingEvent() const { return fUseEmbedding; }

  /**
   * Set the input event of the current event. It will be ignored if fUseEmbedding
   * is true.
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
