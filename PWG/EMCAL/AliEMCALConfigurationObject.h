/*
 * AliEMCALConfigurationObject.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef PWG_EMCAL_ALIEMCALCONFIGURATIONOBJECT_H_
#define PWG_EMCAL_ALIEMCALCONFIGURATIONOBJECT_H_

#include <ostream>
#include <TNamed.h>

class AliEMCALConfigurationValue : public TObject {
public:
  AliEMCALConfigurationValue() {}
  virtual ~AliEMCALConfigurationValue() {}

  virtual const char *ToString() const = 0;

  ClassDef(AliEMCALConfigurationValue,1);
};

class AliEMCALConfigurationValueInt : public AliEMCALConfigurationValue{
public:
  AliEMCALConfigurationValueInt(Int_t val):
    AliEMCALConfigurationValue(),
    fValue(val)
  {}
  virtual ~AliEMCALConfigurationValueInt() {}

  void SetValue(Int_t value) { fValue = value; }
  Int_t GetValue() const { return fValue; }

  virtual const char *ToString() const ;
private:
  Int_t fValue;

  ClassDef(AliEMCALConfigurationValueInt,1);
};

class AliEMCALConfigurationValueFloat : public AliEMCALConfigurationValue{
public:
  AliEMCALConfigurationValueFloat(Float_t val):
    AliEMCALConfigurationValue(),
    fValue(val)
  {}
  virtual ~AliEMCALConfigurationValueFloat() {}

  void SetValue(Float_t value) { fValue = value; }
  Float_t GetValue() const { return fValue; }

  virtual const char *ToString() const;
private:
  Float_t fValue;

  ClassDef(AliEMCALConfigurationValueFloat,1);
};

class AliEMCALConfigurationValueDouble : public AliEMCALConfigurationValue{
public:
  AliEMCALConfigurationValueDouble(Double_t val):
    AliEMCALConfigurationValue(),
    fValue(val)
  {}
  virtual ~AliEMCALConfigurationValueDouble() {}

  void SetValue(Double_t value) { fValue = value; }
  Double_t GetValue() const { return fValue; }

  virtual const char *ToString() const;
private:
  Double_t fValue;

  ClassDef(AliEMCALConfigurationValueDouble,1);
};

class AliEMCALConfigurationValueBool : public AliEMCALConfigurationValue{
public:
  AliEMCALConfigurationValueBool(Bool_t val):
    AliEMCALConfigurationValue(),
    fValue(val)
  {}
  virtual ~AliEMCALConfigurationValueBool() {}

  void SetValue(Bool_t value) { fValue = value; }
  Bool_t GetValue() const { return fValue; }

  virtual const char *ToString() const { return fValue ? "true" : "false"; }
private:
  Bool_t fValue;

  ClassDef(AliEMCALConfigurationValueBool,1);
};

class AliEMCALConfigurationValueString : public AliEMCALConfigurationValue{
public:
  AliEMCALConfigurationValueString(const char *val):
    AliEMCALConfigurationValue(),
    fValue(val)
  {}
  virtual ~AliEMCALConfigurationValueString() {}

  void SetValue(const char * value) { fValue = value; }
  const char * GetValue() const { return fValue; }

  virtual const char *ToString() const { return fValue; }
private:
  TString fValue;

  ClassDef(AliEMCALConfigurationValueString,1);
};

class AliEMCALConfigurationObject : public TNamed {
public:
  AliEMCALConfigurationObject(const char *name, AliEMCALConfigurationValue *value):
    TNamed(name, ""),
    fValue(value)
  {}
  AliEMCALConfigurationObject(const char *key, const char *value);

  virtual ~AliEMCALConfigurationObject(){
    delete fValue;
  }

  void SetValue(AliEMCALConfigurationValue *val){
    if(fValue) delete fValue;
    fValue = val;
  }

  AliEMCALConfigurationValue *GetValue() const { return fValue; }
  const char *ToString() const;

protected:
  AliEMCALConfigurationObject(const AliEMCALConfigurationObject &ref);
  AliEMCALConfigurationObject &operator=(const AliEMCALConfigurationObject &ref);

  AliEMCALConfigurationValue *fValue;

  ClassDef(AliEMCALConfigurationObject, 1);
};

std::ostream &operator<<(std::ostream &, const AliEMCALConfigurationValue &);
std::ostream &operator<<(std::ostream &, const AliEMCALConfigurationObject &);

#endif /* PWG_EMCAL_ALIEMCALCONFIGURATIONOBJECT_H_ */
