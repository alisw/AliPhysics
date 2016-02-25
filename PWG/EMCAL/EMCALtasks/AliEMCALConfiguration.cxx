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

#include "AliJSONReader.h"
#include "AliEMCALConfiguration.h"
#include "AliJSONData.h"


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
    AliJSONValue* value) {
  AliJSONData *entry = dynamic_cast<AliJSONData *>(fParams->FindObject(name));
  if(entry) entry->SetValue(value);
  else fParams->Add(new AliJSONData(name, value));
}

void AliEMCALConfiguration::Build(const char *jsonstring) {
  AliJSONReader parser;
  Build(parser.Decode(jsonstring));
}

void AliEMCALConfiguration::Build(TList *entries){
  TIter objects(entries);
  for(TIter entry = objects.Begin(); entry != objects.End(); ++entry){
    AliJSONData *val = dynamic_cast<AliJSONData *>(*entry);
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
  AliJSONData *conf(NULL);
  while((conf = dynamic_cast<AliJSONData *>(parIter()))){
    std::cout << "Key " << conf->GetName() << ", value " << conf->GetValue()->ToString() << std::endl;
  }
  std::cout << "=================================================" << std::endl;
}
void AliEMCALConfiguration::AddConfiguration(AliEMCALConfiguration* conf) {
  fParams->Add(conf);
}

AliJSONValue* AliEMCALConfiguration::GetValue(const char *key) const {
  AliJSONData *val = dynamic_cast<AliJSONData *>(fParams->FindObject(key));
  if(!val) return NULL;
  return val->GetValue();
}

std::string AliEMCALConfiguration::CreateJSONString() const {
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
      AliJSONData *obj = dynamic_cast<AliJSONData *>(*it);
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
