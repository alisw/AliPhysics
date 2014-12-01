#include "AliEMCALConfiguration.h"
#include "AliJSONData.h"
#include "AliEMCALConfigurationMatcher.h"

ClassImp(AliEMCALConfigurationMatcher)

AliEMCALConfigurationMatcher::AliEMCALConfigurationMatcher(AliEMCALConfiguration *userConfig, AliEMCALConfiguration *defaultConfig):
  fUserConfiguration(userConfig),
  fDefaultConfiguration(defaultConfig)
{
}

AliJSONValue *AliEMCALConfigurationMatcher::GetValue(const char *key) const {
   if (fUserConfiguration->HasKey(key)) return fUserConfiguration->GetValue(key);
   return fDefaultConfiguration->GetValue(key);
}

