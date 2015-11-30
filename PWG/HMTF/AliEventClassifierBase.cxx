#include <vector>
#include <iostream>

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliEventClassifierBase.h"

using namespace std;

ClassImp(AliEventClassifierBase)

AliEventClassifierBase::AliEventClassifierBase()
 : TNamed(),
  fClassifierValue(0),
  fExpectedMinValue(0),
  fExpectedMaxValue(0),
  fClassifierValueIsCached(false),
  fClassifierOutputList(0),
  fTaskOutputList(0),
  fSafePdgCodes(0)
{
  
}
AliEventClassifierBase::AliEventClassifierBase(const char* name, const char* title, TList *taskOutputList)
  : TNamed(name, title),
    fClassifierValue(0),
    fExpectedMinValue(0),
    fExpectedMaxValue(0),
    fClassifierValueIsCached(false),
    fClassifierOutputList(0),
    fTaskOutputList(taskOutputList),
    fSafePdgCodes(0)
{
  fClassifierOutputList = new TList();
  fClassifierOutputList->SetName(name);
  fTaskOutputList->Add(fClassifierOutputList);

  // asign valid pdgs; detour for cxx98
  Int_t pdgs[] = {111, 211, 321, 2212, 310, 3122, 3312, 3334, 333, 313};
  fSafePdgCodes.assign(&pdgs[0], &pdgs[0] + (sizeof(pdgs) / sizeof(Int_t)));
}


Float_t AliEventClassifierBase::GetClassifierValue(AliMCEvent *event, AliStack *stack) {
  if(!fClassifierValueIsCached) {
    CalculateClassifierValue(event, stack);
    fClassifierValueIsCached = true;
  }
  return fClassifierValue;
}
