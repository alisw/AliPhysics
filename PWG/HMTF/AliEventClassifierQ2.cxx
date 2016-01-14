#include <vector>
#include <iostream>

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"

#include "AliEventClassifierQ2.h"
#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliEventClassifierQ2)

AliEventClassifierQ2::AliEventClassifierQ2(const char* name, const char* title,
					     TList *taskOutputList)
  : AliEventClassifierBase(name, title, taskOutputList)
{
  fExpectedMinValue = 0;
  fExpectedMaxValue = 100;

}

void AliEventClassifierQ2::CalculateClassifierValue(AliMCEvent *event, AliStack *stack) {
  fClassifierValue = 0.0;
  // If it is not a pythia header, this should fail
  AliGenPythiaEventHeader* header = dynamic_cast<AliGenPythiaEventHeader*>(event->GenEventHeader());
  if(!header) {
    return;
  }
  else fClassifierValue = header->GetPtHard();
}
