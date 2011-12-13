#ifndef ALIPOOLSSET_H
#define ALIPOOLSSET_H

#include <TObject.h>
#include "AliPoolN.h"
#include "AliClonesPool.h"


class AliPoolsSet : public TObject 
{
 public:
  //
  enum {kPoolN,           // special pool for numeric arrays (AliPoolN class)
	kPoolExtTrPar,    // pool for AliExternalTrackParam objects (AliClonesPool)
	kPoolTrPoints,    // pool for AliTrackPointArray objects (AliClonesPool)
	kPoolTrITS,       // pool for ITS tracks (AliClonesPool)
	kPoolTrTRD,       // pool for TRD tracks  (AliClonesPool)
	kPoolTrFriend,    // pool for AliESDfriendTrack (AliClonesPool)
	kPoolTPCdEdx,     // pool for AliTPCdEdx (AliClonesPool)
	kPoolTPCSeed,     // pool for TPC seeds
	kPoolTPCKink,     // pool for TPC kinks
	kMaxPools};       // total number of pools
  //
  AliPoolsSet();
  AliPoolsSet(const AliPoolsSet& src) : TObject(src),fPools(src.fPools) {}
  virtual ~AliPoolsSet() {DeletePools();}
  AliPoolsSet& operator=(const AliPoolsSet& src) { if (this!=&src) {this->TObject::operator=(src);fPools=src.fPools;} return *this;}
  //
  AliPoolN*       GetPoolN()          const {return dynamic_cast<AliPoolN*>(GetPool(kPoolN));}
  AliClonesPool*  GetPoolC(Int_t id)  const {return dynamic_cast<AliClonesPool*>(GetPool(id));}
  TObject*        GetPool(Int_t id)   const {return fPools[id];}
  // detector-specific pools require setters
  Bool_t          SetPool(TObject* pool, Int_t id);
  //
  void            ResetPools();
  void            InitPools();
  void            DeletePools();
  void            Print(Option_t* opt="") const;
  //
  // shortcuts
  AliClonesPool*  GetPoolExtTrPar()      const {return GetPoolC(kPoolExtTrPar);}
  AliClonesPool*  GetPoolTrPoints()      const {return GetPoolC(kPoolTrPoints);}
  AliClonesPool*  GetPoolTrFriend()      const {return GetPoolC(kPoolTrFriend);}  
  AliClonesPool*  GetPoolTPCdEdx()       const {return GetPoolC(kPoolTPCdEdx);}
  AliClonesPool*  GetPoolTrITS()         const {return GetPoolC(kPoolTrITS);}  
  AliClonesPool*  GetPoolTrTRD()         const {return GetPoolC(kPoolTrTRD);}  
  AliClonesPool*  GetPoolTPCSeed()       const {return GetPoolC(kPoolTPCSeed);}  
  //
 protected:
  //
  TObjArray            fPools;          //! array of pools
  //
  ClassDef(AliPoolsSet,0) // set of pools
};


#endif
