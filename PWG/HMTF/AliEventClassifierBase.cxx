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
  fTaskOutputList(0)
{
  
}
AliEventClassifierBase::AliEventClassifierBase(const char* name, const char* title, TList *taskOutputList)
  : TNamed(name, title),
    fClassifierValue(0),
    fExpectedMinValue(0),
    fExpectedMaxValue(0),
    fClassifierValueIsCached(false),
    fClassifierOutputList(0),
    fTaskOutputList(taskOutputList)
{
  fClassifierOutputList = new TList();
  fClassifierOutputList->SetName(name);
  fTaskOutputList->Add(fClassifierOutputList);
}


Float_t AliEventClassifierBase::GetClassifierValue(AliMCEvent *event, AliStack *stack) {
  if(!fClassifierValueIsCached) {
    CalculateClassifierValue(event, stack);
    fClassifierValueIsCached = true;
  }
  return fClassifierValue;
}
