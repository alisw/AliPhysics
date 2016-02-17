#include "AliITSMFTParamList.h"

ClassImp(AliITSMFTParamList)

//___________________________________________________
AliITSMFTParamList::AliITSMFTParamList(Int_t n, const Double_t *parVal) 
:  AliParamList(n,parVal)
  ,fParamObj(0) 
{
  // def-ctor
}

//___________________________________________________
AliITSMFTParamList::AliITSMFTParamList(const AliITSMFTParamList& src)
  : AliParamList(src)
  , fParamObj( src.fParamObj ? (TObjArray*) src.fParamObj->Clone() : 0)
{
  // copy c-tor
}

//_____________________________________________________________________
AliITSMFTParamList& AliITSMFTParamList::operator=(const AliITSMFTParamList& src)
{
  // copy op.
  if (this != &src) {
    this->~AliITSMFTParamList();
    new(this) AliITSMFTParamList(src);
  }
  return *this;
  //
}

//_____________________________________________________________________
AliITSMFTParamList::~AliITSMFTParamList()
{
  // d-tor
  delete fParamObj;
}

//_____________________________________________________________________
void AliITSMFTParamList::AddParamObject(TObject* obj)
{
  // add new custom object
  if (!fParamObj) {
    fParamObj = new TObjArray();
    fParamObj->SetOwner();
  }
  fParamObj->AddLast(obj);
  
}

//_____________________________________________________________________
void AliITSMFTParamList::Print(Option_t *opt) const
{
  // print itself
  AliParamList::Print(opt);
  //
  if (fParamObj) fParamObj->Print();
  //
}
