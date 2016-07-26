/*
 * AliEMCALConfiguration.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef _ALIEMCALCONFIGURATION_H_
#define _ALIEMCALCONFIGURATION_H_

#include <ostream>
#include <string>
#include <TNamed.h>

class TList;

class AliJSONValue;

class AliEMCALConfiguration : public TNamed {
public:
  AliEMCALConfiguration(const char *name);
  virtual ~AliEMCALConfiguration();

  void AddParam(const char *name, AliJSONValue *value);
  void AddConfiguration(AliEMCALConfiguration * conf);
  void Build(const char * jsonstring);
  void Build(TList *entries);
  std::string CreateJSONString() const;

  Bool_t HasKey(const char *key) const { return GetValue(key) != NULL; }
  AliJSONValue *GetValue(const char *key) const ;
  void Print(Option_t *) const;

protected:
  TList *fParams;

private:
  AliEMCALConfiguration(const AliEMCALConfiguration &ref);
  AliEMCALConfiguration &operator=(const AliEMCALConfiguration &ref);

  ClassDef(AliEMCALConfiguration, 1);
};

std::ostream &operator<<(std::ostream &, const AliEMCALConfiguration &); 

#endif /* _ALIEMCALCONFIGURATION_H_ */
