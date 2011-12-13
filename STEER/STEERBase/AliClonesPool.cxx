#include <TClass.h>
#include <TString.h>
#include "AliLog.h"
#include "AliClonesPool.h"

ClassImp(AliClonesPool)


//_______________________________________________________________________________
AliClonesPool::AliClonesPool(const char* classname, Int_t size)
: TClonesArray(classname,size), fFreeID(size>1000 ? size/100:10), fNFree(0), fLastID(-1), fCastToTrack(kFALSE)
{
  // ctor, delegate to TClonesArray
  TClass* cl = TClass::GetClass(classname);
  if (cl->InheritsFrom(AliVParticle::Class())) fCastToTrack = kTRUE;
  else AliInfo(Form("Pool will use TObject::GetUniqueID to accout %s clones",classname));
  //
}

//_______________________________________________________________________________
AliClonesPool::AliClonesPool(const TClass* cl, Int_t size)
  : TClonesArray(cl,size), fFreeID(size>1000 ? size/100:10), fNFree(0), fLastID(-1), fCastToTrack(kFALSE)
{
  // ctor, delegate to TClonesArray
  if (cl->InheritsFrom(AliVParticle::Class())) fCastToTrack = kTRUE;
  else AliInfo(Form("Pool will use TObject::GetUniqueID to accout %s clones",cl->ClassName()));
  //
}

//_______________________________________________________________________________
AliClonesPool::AliClonesPool(const AliClonesPool& src)
  : TClonesArray(src), fFreeID(src.fFreeID), fNFree(src.fNFree), fLastID(src.fLastID), fCastToTrack(src.fCastToTrack)
{
  // ctor, delegate to TClonesArray
}

//_______________________________________________________________________________
AliClonesPool& AliClonesPool::operator=(const AliClonesPool& src)
{
  if (this!=&src) {
    TClonesArray::operator=(src);
    fFreeID = src.fFreeID;
    fNFree  = src.fNFree;
    fLastID = src.fLastID;
    fCastToTrack = src.fCastToTrack;
  }
  return *this;
}

//_______________________________________________________________________________
AliClonesPool::~AliClonesPool()
{
  // d-tor
  Delete();
}

//_______________________________________________________________________________
void AliClonesPool::Clear(Option_t* opt)
{
  // reset all
  TClonesArray::Clear(opt);
  fNFree = 0;
  fLastID = -1;
}

//_______________________________________________________________________________
void AliClonesPool::MarkSlotFree(TObject *sd) 
{
  // account that this seed is "deleted" 
  int id;
  if (!sd || (id=GetCloneID(sd))<0 ) return;
  if (!IsReset()) {
    //  if (id<0) {AliError(Form("Freeing of seed %p NOT from the pool is requested",sd)); return;}
    sd->Clear("");
    RemoveAt(id);
    if (fFreeID.GetSize()<=fNFree) fFreeID.Set( 2*fNFree + 100 );
    fFreeID.GetArray()[fNFree++] = id;
  }
  //  if (fCastToTrack) ((AliVParticle*)sd)->SetPoolID(-1);
  //  else              sd->SetUniqueID(0);
}

//_______________________________________________________________________________
void AliClonesPool::PrintSummary(Option_t *opt) const
{
  // print summary
  printf("Pool:%s, Booked %d, Size: %d, Freed: %d\n",GetName(), GetEntriesFast(),GetSize(),fNFree);
  TString optS = opt; optS.ToLower();
  if (optS.Contains("l")) TClonesArray::Print();
}
