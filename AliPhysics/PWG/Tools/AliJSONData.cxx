/*
 * AliJSONData.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */
#include <sstream>
#include <TString.h>

#include <AliJSONData.h>

ClassImp(AliJSONValue)
ClassImp(AliJSONInt)
ClassImp(AliJSONFloat)
ClassImp(AliJSONDouble)
ClassImp(AliJSONBool)
ClassImp(AliJSONString)
ClassImp(AliJSONData)

std::string AliJSONInt::ToString() const {
  std::stringstream stringbuilder;
  stringbuilder << fValue;
  return stringbuilder.str();
}

std::string AliJSONFloat::ToString() const {
  std::stringstream stringbuilder;
  stringbuilder << fValue;
  return stringbuilder.str();
}

std::string AliJSONDouble::ToString() const {
  std::stringstream stringbuilder;
  stringbuilder << fValue;
  return stringbuilder.str();
}

AliJSONData::AliJSONData(const char* key, const char* value):
  TNamed(key, ""),
  fValue(NULL)
{
   TString valstring(value);
   if(!valstring.CompareTo("true"))
     fValue = new AliJSONBool(kTRUE);
   else if(!valstring.CompareTo("false"))
     fValue = new AliJSONBool(kFALSE);
   else if(valstring.IsDigit()){
     if(valstring.IsFloat())
       fValue = new AliJSONDouble(valstring.Atof());
     else
       fValue = new AliJSONInt(valstring.Atoi());
   } else
     fValue = new AliJSONString(value);
}

std::string AliJSONData::ToString() const {
  std::stringstream jsonbuilder;
  jsonbuilder << "\"" << GetName() << "\":\"" << fValue->ToString() << "\"";
  return jsonbuilder.str().c_str();
}

std::ostream &operator<<(std::ostream &os, const AliJSONValue &val){
  os << val.ToString();
  return os;
}

std::ostream &operator<<(std::ostream &os, const AliJSONData &obj){
  os << obj.ToString();
  return os;
}

