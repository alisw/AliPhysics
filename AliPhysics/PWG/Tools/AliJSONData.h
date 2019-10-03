/*
 * AliEMCALConfigurationObject.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef _ALIJSONDATA_H_
#define _ALIJSONDATA_H_

#include <ostream>
#include <string>
#include <TNamed.h>
#include <TString.h>

/************************************************************
 *
 * Declaration of wrapper types used in the JSON handling
 *
 * **********************************************************/

class AliJSONValue : public TObject {
public:
  AliJSONValue() {}
  virtual ~AliJSONValue() {}

  virtual std::string ToString() const = 0;

  ClassDef(AliJSONValue,1);
};

class AliJSONInt : public AliJSONValue{
public:
  AliJSONInt(Int_t val):
    AliJSONValue(),
    fValue(val)
  {}
  virtual ~AliJSONInt() {}

  void SetValue(Int_t value) { fValue = value; }
  Int_t GetValue() const { return fValue; }

  virtual std::string ToString() const ;
private:
  Int_t fValue;

  ClassDef(AliJSONInt,1);
};

class AliJSONFloat : public AliJSONValue{
public:
  AliJSONFloat(Float_t val):
    AliJSONValue(),
    fValue(val)
  {}
  virtual ~AliJSONFloat() {}

  void SetValue(Float_t value) { fValue = value; }
  Float_t GetValue() const { return fValue; }

  virtual std::string ToString() const;
private:
  Float_t fValue;

  ClassDef(AliJSONFloat,1);
};

class AliJSONDouble : public AliJSONValue{
public:
  AliJSONDouble(Double_t val):
    AliJSONValue(),
    fValue(val)
  {}
  virtual ~AliJSONDouble() {}

  void SetValue(Double_t value) { fValue = value; }
  Double_t GetValue() const { return fValue; }

  virtual std::string ToString() const;
private:
  Double_t fValue;

  ClassDef(AliJSONDouble,1);
};

class AliJSONBool : public AliJSONValue{
public:
  AliJSONBool(Bool_t val):
    AliJSONValue(),
    fValue(val)
  {}
  virtual ~AliJSONBool() {}

  void SetValue(Bool_t value) { fValue = value; }
  Bool_t GetValue() const { return fValue; }

  virtual std::string ToString() const { return fValue ? "true" : "false"; }
private:
  Bool_t fValue;

  ClassDef(AliJSONBool,1);
};

class AliJSONString : public AliJSONValue{
public:
  AliJSONString(const char *val):
    AliJSONValue(),
    fValue(val)
  {}
  virtual ~AliJSONString() {}

  void SetValue(const char * value) { fValue = value; }
  const char * GetValue() const { return fValue; }

  virtual std::string ToString() const { return std::string(fValue.Data()); }
private:
  TString fValue;

  ClassDef(AliJSONString,1);
};

/********************************************************************
 *                                                                  *
 * Declaration of the JSON key-value pair                           *
 *                                                                  *
 ********************************************************************/

class AliJSONData : public TNamed {
public:
  AliJSONData(const char *name, AliJSONValue *value):
    TNamed(name, ""),
    fValue(value)
  {}
  AliJSONData(const char *key, const char *value);

  virtual ~AliJSONData(){
    delete fValue;
  }

  void SetValue(AliJSONValue *val){
    if(fValue) delete fValue;
    fValue = val;
  }

  AliJSONValue *GetValue() const { return fValue; }
  std::string ToString() const;

protected:
  AliJSONData(const AliJSONData &ref);
  AliJSONData &operator=(const AliJSONData &ref);

  AliJSONValue *fValue;

  ClassDef(AliJSONData, 1);
};

std::ostream &operator<<(std::ostream &, const AliJSONValue &);
std::ostream &operator<<(std::ostream &, const AliJSONData &);

#endif /* _ALIJSONDATA_H_ */
