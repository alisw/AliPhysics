#include"AliPIDperfCut.h"

ClassImp(AliPIDperfCut);

AliPIDperfCut::AliPIDperfCut(const char *name):
  TNamed(name,name)
{}
//--------------------------------------------------------------------
AliPIDperfCut::AliPIDperfCut():
  TNamed("PIDperfCut","PIDperfCut")
{}
//--------------------------------------------------------------------
Bool_t AliPIDperfCut::IsSelected(AliVTrack */*track*/,AliPID::EParticleType /*type*/) const{
  return 1;
}
