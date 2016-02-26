/*
 * AliEMCALConfigHandler.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#include <sstream>
#include <TList.h>

#include "AliEMCALConfigHandler.h"
#include "AliEMCALConfiguration.h"
#include "AliJSONData.h"

ClassImp(AliEMCALConfigHandler)

AliEMCALConfigHandler::AliEMCALConfigHandler() :
  TObject(),
  fConfiguration(NULL)
{
  fConfiguration = new TList;
  fConfiguration->SetOwner();
}

AliEMCALConfigHandler::~AliEMCALConfigHandler() {
  delete fConfiguration;
}

AliEMCALConfiguration *AliEMCALConfigHandler::CreateConfiguration(const char* name) {
  AliEMCALConfiguration *config = FindConfiguration(name);
  if(config)
    this->Error("AliEMCALConfigHandler::CreateConfiguration", "Configuration with name %s already exists, not duplicating", name);
  else{
    config = new AliEMCALConfiguration(name);
    fConfiguration->Add(config);
  }
  return config;
}

void AliEMCALConfigHandler::AddParam(const char* configName, const char* key,
    AliJSONValue* value) {
  AliEMCALConfiguration *config = FindConfiguration(configName);
  if(!config){
    this->Warning("AliEMCALConfigHandler", "Configuration with name %s does not exist, creating it", configName);
    config = CreateConfiguration(configName);
  }
  config->AddParam(key, value);
}

std::string AliEMCALConfigHandler::GetConfigurationString() const {
  std::stringstream jsonbuilder;
  jsonbuilder << "{";
  TIter confentries(fConfiguration);
  bool isFirst = true;
  for(TIter it = confentries.Begin(); it != confentries.End(); ++it){
    AliEMCALConfiguration *conf = static_cast<AliEMCALConfiguration *>(*it);
    if(!isFirst) jsonbuilder << ",";
    jsonbuilder << "\"" << conf->GetName() << "\":" << conf->CreateJSONString();
    if(isFirst) isFirst = false;
  }
  jsonbuilder << "}";
  return jsonbuilder.str();
}

std::string AliEMCALConfigHandler::GetConfigurationString(const char* configname) const {
  AliEMCALConfiguration *conf = FindConfiguration(configname);
  if(!conf) return "";
  return conf->CreateJSONString();
}

AliEMCALConfiguration* AliEMCALConfigHandler::FindConfiguration(const char* configName) const {
  return dynamic_cast<AliEMCALConfiguration *>(fConfiguration->FindObject(configName));
}
