// AliEmcalCorrectionEventManager
//

#include "AliEmcalCorrectionEventManager.h"
#include "AliEmcalContainerUtils.h"

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
  return AliEmcalContainerUtils::GetEvent(fInputEvent, fUseEmbedding);
}
