#include <vector>
#include <iostream>

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"

#include "AliEventClassifierMPI.h"
#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliEventClassifierMPI)

AliEventClassifierMPI::AliEventClassifierMPI(const char* name, const char* title,
					     TList *taskOutputList)
  : AliEventClassifierBase(name, title, taskOutputList)
{
  fExpectedMinValue = 0;
  fExpectedMaxValue = 250;
}

void AliEventClassifierMPI::CalculateClassifierValue(AliMCEvent *event, AliStack *stack) {
  fClassifierValue = 0.0;
  // If it is not a pythia header, this should fail
  AliGenPythiaEventHeader* header = dynamic_cast<AliGenPythiaEventHeader*>(event->GenEventHeader());
  if(!header) {
    return;
  }
  else fClassifierValue = header->GetNMPI();
}
