#include "AliITSUParamList.h"

ClassImp(AliITSUParamList)

//___________________________________________________
AliITSUParamList::AliITSUParamList(Int_t n, const Double_t *parVal) 
:  AliParamList(n,parVal)
  ,fParamObj(0) 
{
  // def-ctor
}

//___________________________________________________
AliITSUParamList::AliITSUParamList(const AliITSUParamList& src)
  : AliParamList(src)
  , fParamObj( src.fParamObj ? (TObjArray*) src.fParamObj->Clone() : 0)
{
  // copy c-tor
}

//_____________________________________________________________________
AliITSUParamList& AliITSUParamList::operator=(const AliITSUParamList& src)
{
  // copy op.
  if (this != &src) {
    this->~AliITSUParamList();
    new(this) AliITSUParamList(src);
  }
  return *this;
  //
}

//_____________________________________________________________________
AliITSUParamList::~AliITSUParamList()
{
  // d-tor
  delete fParamObj;
}

//_____________________________________________________________________
void AliITSUParamList::AddParamObject(TObject* obj)
{
  // add new custom object
  if (!fParamObj) {
    fParamObj = new TObjArray();
    fParamObj->SetOwner();
  }
  fParamObj->AddLast(obj);
  
}

//_____________________________________________________________________
void AliITSUParamList::Print(Option_t *opt) const
{
  // print itself
  AliParamList::Print(opt);
  //
  if (fParamObj) fParamObj->Print();
  //
}
