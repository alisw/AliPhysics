/*
 * AliEMCALConfigHandler.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef PWG_EMCAL_ALIEMCALCONFIGHANDLER_H_
#define PWG_EMCAL_ALIEMCALCONFIGHANDLER_H_

#include <TObject.h>

class AliEMCALConfiguration;
class AliEMCALConfigurationValue;
class TList;

class AliEMCALConfigHandler : public TObject {
public:
  AliEMCALConfigHandler();
  virtual ~AliEMCALConfigHandler();

  AliEMCALConfiguration *CreateConfiguration(const char *name);
  void AddParam(const char *configName, const char *key, AliEMCALConfigurationValue *value);
  const char *GetConfigurationString(const char *configname) const;
  AliEMCALConfiguration *FindConfiguration(const char *configName);

protected:
  TList *fConfiguration;

  ClassDef(AliEMCALConfigHandler, 1);
};

#endif /* PWG_EMCAL_ALIEMCALCONFIGHANDLER_H_ */
