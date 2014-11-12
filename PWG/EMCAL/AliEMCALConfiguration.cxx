/*
 * AliEMCALConfiguration.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */
#include <cstring>
#include <sstream>
#include <iostream>
#include <TList.h>

#include "AliEMCALJSONReader.h"
#include "AliEMCALConfiguration.h"
#include "AliEMCALConfigurationObject.h"


AliEMCALConfiguration::AliEMCALConfiguration(const char* name) :
  TNamed(name, ""),
  fParams(NULL)
{
  fParams = new TList;
  fParams->SetOwner();
}

AliEMCALConfiguration::~AliEMCALConfiguration() {
  if(!fParams) delete fParams;
}

void AliEMCALConfiguration::AddParam(const char* name,
    AliEMCALConfigurationValue* value) {
  AliEMCALConfigurationObject *entry = dynamic_cast<AliEMCALConfigurationObject *>(fParams->FindObject(name));
  if(entry) entry->SetValue(value);
  else fParams->Add(new AliEMCALConfigurationObject(name, value));
}

void AliEMCALConfiguration::Build(const char *jsonstring) {
  AliEMCALJSONReader parser;
  Build(parser.Decode(jsonstring));
}

void AliEMCALConfiguration::Build(TList *entries){
  TIter objects(entries);
  for(TIter entry = objects.Begin(); entry != objects.End(); ++entry){
    AliEMCALConfigurationObject *val = dynamic_cast<AliEMCALConfigurationObject *>(*entry);
    if(val)
      fParams->Add(val);
    else{
      TList *conf = dynamic_cast<TList *>(*entry);
      if(conf){
        AliEMCALConfiguration *daughter = new AliEMCALConfiguration(conf->GetName());
        daughter->Build(conf);
        fParams->Add(daughter);
      }
    }
  }
}

void AliEMCALConfiguration::Print(Option_t * /*value*/) const {
  std::cout << "Configuration " << GetName() << ":" << std::endl;
  std::cout << "=================================================" << std::endl;
  TIter parIter(fParams);
  AliEMCALConfigurationObject *conf(NULL);
  while((conf = dynamic_cast<AliEMCALConfigurationObject *>(parIter()))){
    std::cout << "Key " << conf->GetName() << ", value " << conf->GetValue()->ToString() << std::endl;
  }
  std::cout << "=================================================" << std::endl;
}
void AliEMCALConfiguration::AddConfiguration(AliEMCALConfiguration* conf) {
  fParams->Add(conf);
}

AliEMCALConfigurationValue* AliEMCALConfiguration::GetValue(const char *key) {
  AliEMCALConfigurationObject *val = dynamic_cast<AliEMCALConfigurationObject *>(fParams->FindObject(key));
  if(!val) return NULL;
  return val->GetValue();
}

const char* AliEMCALConfiguration::CreateJSONString() const {
  std::stringstream jsonbuilder;
  jsonbuilder << "{";
  TIter confentries(fParams);
  bool isFirst = true;
  for(TIter it = confentries.Begin(); it != confentries.End(); ++it){
    AliEMCALConfiguration *conf = dynamic_cast<AliEMCALConfiguration *>(*it);
    if(conf){
      if(!isFirst) jsonbuilder << ",";
      jsonbuilder << "\"" << conf->GetName() << "\":" << conf->CreateJSONString();
    } else {
      AliEMCALConfigurationObject *obj = dynamic_cast<AliEMCALConfigurationObject *>(*it);
      if(obj){
        if(!isFirst) jsonbuilder << ",";
        jsonbuilder << obj->ToString();
      }
    }
    if(isFirst) isFirst = false;
  }
  jsonbuilder << "}";

  char * result = new char[jsonbuilder.str().length()];
  strcpy(result, jsonbuilder.str().c_str());  
  return result;
}

std::ostream &operator<<(std::ostream & os, const AliEMCALConfiguration &conf){
  os << conf.CreateJSONString();
  return os;
}
