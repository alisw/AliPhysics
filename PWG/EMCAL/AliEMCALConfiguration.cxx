/*
 * AliEMCALConfiguration.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */
#include <sstream>
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

void AliEMCALConfiguration::AddConfiguration(AliEMCALConfiguration* conf) {
  fParams->Add(conf);
}

TObject* AliEMCALConfiguration::GetValue(const char *key) {
  return fParams->FindObject(key);
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
  return jsonbuilder.str().c_str();
}
