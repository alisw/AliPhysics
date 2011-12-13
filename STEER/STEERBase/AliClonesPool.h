#ifndef ALICLONESPOOL_H
#define ALICLONESPOOL_H

#include <TClonesArray.h>
#include <TArrayI.h>
#include "AliVParticle.h"

class AliClonesPool: public TClonesArray 
{
 public:
  AliClonesPool() : fFreeID(),fNFree(0),fLastID(-1),fCastToTrack(kFALSE) {}
  AliClonesPool(const char* classname, Int_t size = 1000);
  AliClonesPool(const TClass* cl, Int_t size = 1000);
  AliClonesPool(const AliClonesPool& src);
  AliClonesPool& operator=(const AliClonesPool& src);
  virtual ~AliClonesPool();
  virtual void Clear(Option_t* option="");
  //
  Int_t  GetLastID()              const {return fLastID;}
  Bool_t IsReset()                const {return fLastID<0;}
  TObject *&NextFreeSlot();
  void MarkSlotFree(TObject *sd);
  void RegisterClone(TObject *sd) const {fCastToTrack ? ((AliVParticle*)sd)->SetPoolID(fLastID) : sd->SetUniqueID(UInt_t(fLastID+1));}
  void PrintSummary(const Option_t* opt="") const;
  //
 protected:
  Int_t GetCloneID(TObject* obj) const {return fCastToTrack ? ((AliVParticle*)obj)->GetPoolID() : (obj->GetUniqueID()-1);}

 protected:
  TArrayI fFreeID;              // array of ID's of free slots
  Int_t   fNFree;               // number of slots freed in the pool
  Int_t   fLastID;              // id of the slot which is returned by the NextFreeSlot method
  Bool_t  fCastToTrack;         // AliVParticle-derived classes has GetPoolID method, use it
  //
  ClassDef(AliClonesPool,0)   // pool of clones
};

inline TObject *& AliClonesPool::NextFreeSlot() {
  // return slot where next object will be created (using new with placement)
  return (*this)[ fLastID = fNFree ? fFreeID.GetArray()[--fNFree] : GetEntriesFast() ];
}


#endif
