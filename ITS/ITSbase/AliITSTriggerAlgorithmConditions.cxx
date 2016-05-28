////////////////////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                                         //
//                                                                                //
// Class for storing conditions data from Pixel Trigger (PIT) algorithms.         //
// This holds a sub set of the conditions data needed.                            //
// It is used by AliITSTriggerConditions, which holds all the information.        //
// AliITSTriggerConditions contains a TObjArray of this type.                     //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

#include "AliITSTriggerAlgorithmConditions.h"
#include <TObjString.h>

ClassImp(AliITSTriggerAlgorithmConditions)

//__________________________________________________________________________
AliITSTriggerAlgorithmConditions::AliITSTriggerAlgorithmConditions():
TObject(),
fId(0),
fLabel(TString("label")),
fDescription(TString("descr")),
fNumParam(0),
fParamNames(TObjArray(3)),
fParamValues(TArrayI(3))
{
  // default constructor
  fParamNames.SetOwner(kTRUE);
}
//__________________________________________________________________________
AliITSTriggerAlgorithmConditions::AliITSTriggerAlgorithmConditions(UShort_t id, const Char_t* label, const Char_t* descr):
TObject(),
fId(id),
fLabel(label),
fDescription(descr),
fNumParam(0),
fParamNames(TObjArray(3)),
fParamValues(TArrayI(3))
{
  // optional constructor
  fParamNames.SetOwner(kTRUE);
}
//__________________________________________________________________________
AliITSTriggerAlgorithmConditions::AliITSTriggerAlgorithmConditions(const AliITSTriggerAlgorithmConditions& cond):
TObject(),
fId(cond.fId),
fLabel(cond.fLabel),
fDescription(cond.fDescription),
fNumParam(cond.fNumParam),
fParamNames(cond.fParamNames),
fParamValues(cond.fParamValues)
{
  // default constructor
  fParamNames.SetOwner(kTRUE);
}
//__________________________________________________________________________
AliITSTriggerAlgorithmConditions::~AliITSTriggerAlgorithmConditions() 
{
  // destructor
  ClearParams();
}
//__________________________________________________________________________
AliITSTriggerAlgorithmConditions& AliITSTriggerAlgorithmConditions::operator=(const AliITSTriggerAlgorithmConditions& cond) {
  // assignment operator
  if (this!=&cond) {
    fId = cond.fId;
    fLabel = cond.fLabel;
    fDescription = cond.fDescription;
    fNumParam = cond.fNumParam;
    fParamNames = cond.fParamNames;
    fParamValues = cond.fParamValues;
  }
  return *this;
}
//__________________________________________________________________________
void AliITSTriggerAlgorithmConditions::ClearParams() {
  // clears parameter list
  fParamNames.Clear();
  fNumParam=0;
}
//__________________________________________________________________________
void AliITSTriggerAlgorithmConditions::AddParam(const Char_t* name, Int_t value) {
  // adds a new parameter with name 'name' and value 'value'
  // if the name is already present in the list, the parameter value will be over-written
  UShort_t findIndex=fNumParam;
  for (UInt_t i=0; i<fNumParam; i++) {
    if (((TObjString*)fParamNames.At(i))->String().CompareTo(name, TString::kIgnoreCase) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex<fNumParam) {
    fParamValues[findIndex]=value;
  }
  else {
    fParamNames.AddAtAndExpand(new TObjString(name),fNumParam);
    Int_t valSize = fParamValues.GetSize();
    if (valSize<=fNumParam) fParamValues.Set(valSize*2);
    fParamValues[fNumParam]=value;
    fNumParam++;
  }
}
//__________________________________________________________________________
const Char_t* AliITSTriggerAlgorithmConditions::GetParamNameI(UShort_t index) const {
  // returns parameter name for parameter at position index
  if (index>=fNumParam) {
    Error("AliITSTriggerAlgorithmConditions::GetParamNameI", "index %d out of range", index);
    return "dummy";
  }
  return ((TObjString*)fParamNames.At(index))->String().Data();
}
//__________________________________________________________________________
Int_t AliITSTriggerAlgorithmConditions::GetParamValueI(UShort_t index) const {
  // returns paramter value at position index
  if (index>=fNumParam) {
    Error("AliITSTriggerAlgorithmConditions::GetParamValueI", "index %d out of range", index);
    return -1;
  }
  return fParamValues.At(index);
}
//__________________________________________________________________________
Int_t AliITSTriggerAlgorithmConditions::GetParamValueN(const Char_t* name) const {
  // returns parameter value for parameter with name 'name'
  UShort_t findIndex=fNumParam;
  for (UInt_t i=0; i<fNumParam; i++) {
    if (((TObjString*)fParamNames.At(i))->String().CompareTo(name, TString::kIgnoreCase) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumParam) {
    Error("AliITSTriggerAlgorithmConditions::GetParamValueN", "name %s not found", name);
    return -1;
  }
  return fParamValues.At(findIndex);
}
