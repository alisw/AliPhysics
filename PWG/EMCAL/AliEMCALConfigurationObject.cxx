/*
 * AliEMCALConfigurationObject.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */
#include <sstream>
#include <TString.h>

#include <AliEMCALConfigurationObject.h>


const char* AliEMCALConfigurationValueInt::ToString() const { 
  std::stringstream stringbuilder;
  stringbuilder << fValue;
  return stringbuilder.str().c_str();
}

const char* AliEMCALConfigurationValueFloat::ToString() const {
  std::stringstream stringbuilder;
  stringbuilder << fValue;
  return stringbuilder.str().c_str();
}

const char* AliEMCALConfigurationValueDouble::ToString() const {
  std::stringstream stringbuilder;
  stringbuilder << fValue;
  return stringbuilder.str().c_str();
}

AliEMCALConfigurationObject::AliEMCALConfigurationObject(const char* key, const char* value):
  TNamed(key, ""),
  fValue(NULL)
{
   TString valstring(value);
   if(!valstring.CompareTo("true"))
     fValue = new AliEMCALConfigurationValueBool(kTRUE);
   else if(!valstring.CompareTo("false"))
     fValue = new AliEMCALConfigurationValueBool(kFALSE);
   else if(valstring.IsDigit()){
     if(valstring.IsFloat())
       fValue = new AliEMCALConfigurationValueDouble(valstring.Atof());
     else
       fValue = new AliEMCALConfigurationValueInt(valstring.Atoi());
   } else
     fValue = new AliEMCALConfigurationValueString(value);
}

const char* AliEMCALConfigurationObject::ToString() const {
  std::stringstream jsonbuilder;
  jsonbuilder << "\"" << GetName() << "\":\"" << fValue->ToString() << "\"";
  return jsonbuilder.str().c_str();
}
