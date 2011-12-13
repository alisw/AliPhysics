#include "AliPoolsSet.h"
#include "AliLog.h"

ClassImp(AliPoolsSet)

//__________________________________________
AliPoolsSet::AliPoolsSet() : fPools(kMaxPools)
{
  // def ctor
  fPools.SetOwner();
}

//___________________________________________________________
void AliPoolsSet::InitPools()
{
  // Create obligatory pools (which always needed). 
  // The pools for specific detector objects are initialized by respective trackers of AliReconstruction
  AliInfo("Initializing pools>>");
  if (!GetPoolN()) SetPool(new AliPoolN(100000), kPoolN);
  if (!GetPool(kPoolExtTrPar)) {
    SetPool(new AliClonesPool("AliExternalTrackParam",5000), kPoolExtTrPar);
    GetPoolC(kPoolExtTrPar)->SetName("ExternalTrackParam");
  }
  if (!GetPool(kPoolTrFriend)) {
    SetPool(new AliClonesPool("AliESDfriendTrack",5000),kPoolTrFriend);
    GetPoolC(kPoolTrFriend)->SetName("ESDfriendTrack");
  }
  if (!GetPool(kPoolTPCdEdx)) {
    SetPool(new AliClonesPool("AliTPCdEdxInfo",5000), kPoolTPCdEdx);
    GetPoolC(kPoolTPCdEdx)->SetName("TPCdEdxInfo");
  }
  //
  AliInfo("Initializing pools<<");
}

//___________________________________________________________
void AliPoolsSet::DeletePools()
{
  // delete all global pools
  fPools.Delete();
}

//___________________________________________________________
void AliPoolsSet::ResetPools()
{
  // reset all pools
  for (int i=kMaxPools;i--;) {
    TObject* pool = GetPool(i);
    if (!pool) continue;
    pool->Clear("C");
  }
}

//___________________________________________________________
Bool_t AliPoolsSet::SetPool(TObject* pool, Int_t id)
{
  // attach the pool
  if (!pool) return kFALSE;
  if (id<0 || id>=kMaxPools) AliFatal(Form("Defined pool id's are 0:%d, %d requested",kMaxPools,id));
  if (!pool->InheritsFrom("AliPoolN") && !pool->InheritsFrom("AliClonesPool")) 
    AliFatal(Form("Supported pool types: %s,%s | supplied %s","AliPoolN","AliClonesPool",pool->ClassName()));
  //
  fPools[id] = pool;
  return kTRUE;
}

//___________________________________________________________
void AliPoolsSet::Print(Option_t* opt) const
{
  // print summary
  for (int i=0;i<kMaxPools;i++) {
    TObject* pl = GetPool(i);
    if (!pl) continue;
    if      ( pl->IsA() == AliClonesPool::Class() ) ((AliClonesPool*)pl)->PrintSummary(opt);
    else if ( pl->IsA() == AliPoolN::Class() )      ((AliPoolN*)pl)->PrintSummary(opt);
  }
}
