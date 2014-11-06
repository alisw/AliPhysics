/*
 * AliEMCALConfiguration.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef PWG_EMCAL_ALIEMCALCONFIGURATION_H_
#define PWG_EMCAL_ALIEMCALCONFIGURATION_H_

#include <TNamed.h>

class TList;

class AliEMCALConfiguration : public TNamed {
public:
  AliEMCALConfiguration(const char *name);
  virtual ~AliEMCALConfiguration();

  void AddParam(const char *name, AliEMCALConfigurationValue *value);
  void AddConfiguration(AliEMCALConfiguration * conf);
  void Build(const char * jsonstring);
  void Build(TList *entries);
  const char *CreateJSONString() const;

  TObject *GetValue(const char *key);

protected:
  TList *fParams;

  ClassDef(AliEMCALConfiguration, 1);
};

#endif /* PWG_EMCAL_ALIEMCALCONFIGURATION_H_ */
