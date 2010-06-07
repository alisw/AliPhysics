//
// Class AliRsnValueList
//
// This class defines a base classe to implement a function
// which uses the internal RSN package event format (AliRsnEvent).
// It contains some default flags which turn out to be useful:
//  - a flag to select only the "true" pairs (tracks from same resonance)
//  - a flag to know if the computation is done over two events (mixing)
//
// Any kind of analysis object should be implemented as inheriting from this
// because the AliRsnAnalyzer which executes the analysis will accept a collection
// of such objects, in order to have a unique format of processing method
//
// The user who implements a kind of computation type should inherit from
// this class and override the virtual functions defined in it, which
// initialize the final output histogram and define how to process data.
//
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//

#include <TString.h>

#include "AliLog.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"
#include "AliRsnPairParticle.h"
#include "AliRsnValue.h"

#include "AliRsnValueList.h"

ClassImp(AliRsnValueList)

//________________________________________________________________________________________
AliRsnValueList::AliRsnValueList() :
  TNamed(),
  fValueList("AliRsnValue", 0),
  fArray(0)
{
//
// Constructor.
//
}

//________________________________________________________________________________________
AliRsnValueList::AliRsnValueList(const AliRsnValueList &copy) :
  TNamed(copy),
  fValueList(copy.fValueList),
  fArray(copy.fArray)
{
//
// Copy constructor.
//
}

//________________________________________________________________________________________
const AliRsnValueList& AliRsnValueList::operator=(const AliRsnValueList& copy)
{
//
// Assignment operator.
//

  // copy name and title
  SetName(copy.GetName());
  SetTitle(copy.GetTitle());

  // add a copy of each axis to this list
  // this will update the value list accordingly
  Int_t i, n = copy.fValueList.GetEntries();
  for (i = 0; i < n; i++)
  {
    AliRsnValue *val = (AliRsnValue*)copy.fValueList[i];
    AddValue(val);
  }

  return *this;
}

//________________________________________________________________________________________
void AliRsnValueList::AddValue(const AliRsnValue *axis)
{
//
// Adds a new value to the list
// and upgrades the list of values accordingly
//

  Int_t size = fValueList.GetEntries();
  new(fValueList[size]) AliRsnValue(*axis);
  fArray.Set(size + 1);
}

//________________________________________________________________________________________
Bool_t AliRsnValueList::Eval(TObject * const obj, const AliRsnPairDef *pairDef)
{
//
// Fill function histogram using the passed object,
// whose type is dynamically casted when used and checked.
//

  AliDebug(AliLog::kDebug +2,"->");

  Int_t  i, size = fValueList.GetSize();
  for (i = 0; i < size; i++)
  {
    AliRsnValue *val = (AliRsnValue*)fValueList[i];
    fArray[i] = val->Eval(obj, pairDef);
  }
}
