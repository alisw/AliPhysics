#include "AliLumiRef.h"
#include "AliLog.h"

ClassImp(AliLumiRef);

AliLumiRef::AliLumiRef(const char* trigger, const char* comment, UInt_t run, float sig, float eff) :
  TNamed(trigger,comment),
  fRefSigma(sig),
  fRefEff(eff)
{
  SetRunStart(run);
}

AliLumiRef::AliLumiRef(const AliLumiRef& src) :
  TNamed(src),
  fRefSigma(src.fRefSigma),
  fRefEff(src.fRefSigma)
{}

AliLumiRef &AliLumiRef::operator=(const AliLumiRef& src)
{
  if(&src == this) return *this;
  TNamed::operator=(src);
  SetRefSigma(src.GetRefSigma());
  SetRefEff(src.GetRefEff());
  return *this;
}

void AliLumiRef::Print(const Option_t *) const
{
  printf("%6d %-30s %8.2f %4.2f | %s\n",
	 GetRunStart(),GetRefTrigger(),GetRefSigma(),GetRefEff(),GetComment());
}

Int_t AliLumiRef::Compare(const TObject *obj) const
{
  // compare according to run number
  AliLumiRef* ref = (AliLumiRef*)obj;
  if (GetRunStart()>ref->GetRunStart()) return 1;
  if (GetRunStart()<ref->GetRunStart()) return -1;
  AliFatalF("Two LumiReferences cannot have the same run %d",GetRunStart());
  return 0;
}
