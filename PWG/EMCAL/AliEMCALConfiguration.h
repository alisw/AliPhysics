/*
 * AliEMCALConfiguration.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef PWG_EMCAL_ALIEMCALCONFIGURATION_H_
#define PWG_EMCAL_ALIEMCALCONFIGURATION_H_

#include <ostream>
#include <TNamed.h>

class TList;

class AliEMCALConfigurationValue;

class AliEMCALConfiguration : public TNamed {
public:
  AliEMCALConfiguration(const char *name);
  virtual ~AliEMCALConfiguration();

  void AddParam(const char *name, AliEMCALConfigurationValue *value);
  void AddConfiguration(AliEMCALConfiguration * conf);
  void Build(const char * jsonstring);
  void Build(TList *entries);
  const char *CreateJSONString() const;

  Bool_t HasKey(const char *key) const { return GetValue(key) != NULL; }
  AliEMCALConfigurationValue *GetValue(const char *key) const ;
  void Print(Option_t *) const;

protected:
  TList *fParams;

private:
  AliEMCALConfiguration(const AliEMCALConfiguration &ref);
  AliEMCALConfiguration &operator=(const AliEMCALConfiguration &ref);

  ClassDef(AliEMCALConfiguration, 1);
};

std::ostream &operator<<(std::ostream &, const AliEMCALConfiguration &); 

#endif /* PWG_EMCAL_ALIEMCALCONFIGURATION_H_ */
