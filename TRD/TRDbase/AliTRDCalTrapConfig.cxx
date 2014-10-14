#include "TString.h"

#include "AliLog.h"
#include "AliTRDCalTrapConfig.h"

ClassImp(AliTRDCalTrapConfig)

AliTRDCalTrapConfig::AliTRDCalTrapConfig() :
  TObject(),
  fConfigList()
{
  // ctor

}

AliTRDCalTrapConfig::~AliTRDCalTrapConfig()
{
  // dtor

}

AliTRDtrapConfig* AliTRDCalTrapConfig::Get(const TString &name) const
{
  return (AliTRDtrapConfig*) fConfigList.FindObject(name.Data());
}

void AliTRDCalTrapConfig::Print(Option_t * /* option */) const
{
  TIter config(&fConfigList);

  while (AliTRDtrapConfig *cfg = (AliTRDtrapConfig*) config()) {
    AliInfo(Form("found TRAPconfig: %s - %s", cfg->GetName(), cfg->GetTitle()));
  }
}
