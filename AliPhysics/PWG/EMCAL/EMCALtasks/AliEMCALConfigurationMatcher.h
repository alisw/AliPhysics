/*
 * AliEMCALConfiguration.h
 *
 *  Created on: 13.11.2014
 *      Author: markusfasel
 */


#ifndef PWG_EMCAL_ALIEMCALCONFIGURATIONMATCHER_H_
#define PWG_EMCAL_ALIEMCALCONFIGURATIONMATCHER_H_

#include<TObject.h>

class AliEMCALConfiguration;
class AliJSONValue;

class AliEMCALConfigurationMatcher : public TObject{
public:
  AliEMCALConfigurationMatcher(AliEMCALConfiguration *userConfig, AliEMCALConfiguration *defaultConfig);
  virtual ~AliEMCALConfigurationMatcher() {}

  void SetUserConfiguration(AliEMCALConfiguration *conf) { fUserConfiguration = conf; }
  void SetDefaultConfiguration(AliEMCALConfiguration *conf) { fDefaultConfiguration = conf; } 

  AliJSONValue *GetValue(const char * key) const;

protected:
  AliEMCALConfiguration *fUserConfiguration;
  AliEMCALConfiguration *fDefaultConfiguration;

private:
  AliEMCALConfigurationMatcher(const AliEMCALConfigurationMatcher &);
  AliEMCALConfigurationMatcher &operator=(const AliEMCALConfigurationMatcher &);
  
  ClassDef(AliEMCALConfigurationMatcher, 1);
};

#endif
