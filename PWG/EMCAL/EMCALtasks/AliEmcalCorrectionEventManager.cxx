// AliEmcalCorrectionEventManager
//

#include "AliEmcalCorrectionEventManager.h"

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

/**
 * Standard constructor
 *
 * fUseEmbedding defaults to true because our condition turns it off a standard input
 * event is seen.
 */
AliEmcalCorrectionEventManager::AliEmcalCorrectionEventManager():
  AliAnalysisTaskSE("AliEmcalCorrectionEventManager"),
  fUseEmbedding(true)
{

}

/**
 * Named constructor.
 *
 * fUseEmbedding defaults to true because our condition turns it off a standard input
 * event is seen.
 */
AliEmcalCorrectionEventManager::AliEmcalCorrectionEventManager(const char * name):
  AliAnalysisTaskSE(name),
  fUseEmbedding(true)
{

}

/**
 * Deconstructor
 */
AliEmcalCorrectionEventManager::~AliEmcalCorrectionEventManager()
{
}

/**
 * Return the proper input event depending on the whether the current or embedded event
 * was configured.
 *
 * @return Requested input event
 */
AliVEvent * AliEmcalCorrectionEventManager::InputEvent() const
{
  AliVEvent * event = fInputEvent;
  if (fUseEmbedding) {
    const auto embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (!embeddingHelper) return 0;

    event = embeddingHelper->GetExternalEvent();
  }

  return event;
}
