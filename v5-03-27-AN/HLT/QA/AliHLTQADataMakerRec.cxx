// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTQADataMakerRec.cxx
/// @author Matthias Richter
/// @date   2010-03-10
/// @brief  Steering class for the HLT offline QA
///
#include "AliHLTQADataMakerRec.h"
#include "AliHLTMisc.h"
#include "AliHLTModuleAgent.h"
#include "AliRecoParam.h"
#include <iostream>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TDirectory.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTQADataMakerRec)

AliHLTQADataMakerRec::AliHLTQADataMakerRec()
  : AliHLTQADataMakerBase()
  , fPlugins()
  , fFlags(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  LoadAgents();
}

AliHLTQADataMakerRec::~AliHLTQADataMakerRec()
{
  // see header file for class documentation
}

int AliHLTQADataMakerRec::LoadAgents()
{
  // iterate over available agents and query class names of plugins
  TString plugins;
  for (AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent(); 
       pAgent!=NULL;
       pAgent=AliHLTModuleAgent::GetNextAgent()) {
    const char* modulePlugins=pAgent->GetQAPlugins();
    if (!modulePlugins || !modulePlugins[0]) continue;
    if (!plugins.IsNull() && !plugins.EndsWith(" ")) plugins+=" ";
    plugins+=modulePlugins;
  }
  if (!plugins.IsNull()) return LoadPlugins(plugins);
  return 0;
}

int AliHLTQADataMakerRec::LoadPlugins(const char* plugins)
{
  // load plugins from list of blank separated class names
  int iResult=0;
  TString strPlugins=plugins;
  TObjArray* tokens=strPlugins.Tokenize(" ");
  if (tokens) {
    TIter next(tokens);
    TObject* obj=NULL;
    while ((obj=next())) {
      TObjString* objstring=dynamic_cast<TObjString*>(obj);
      if (!objstring) continue;
      AliHLTQADataMakerBase* plugin=AliHLTMisc::LoadInstance((AliHLTQADataMakerBase*)0,
							     objstring->GetString().Data());
      if (!plugin) continue;
      AliInfo(Form("using HLT QA plugin %s", plugin->IsA()->GetName()));
      fPlugins.Add(plugin);
    }
    delete tokens;
  }
  return iResult;
}

void AliHLTQADataMakerRec::StartOfDetectorCycle()
{
  // see header file for class documentation
  
  // this function is called multiple times by the framework, actually for every QA task
  // however, here we don't have any argument for the task
  // this class is initialized right before StartOfDetectorCycle is called, and depending
  // on the availibility of thr histogram arrays one can tell where we are ;-)
  unsigned init=0;
  if (fDigitsQAList!=NULL && (fFlags&kDigitsListInit)==0) {
    init|=kDigitsListInit; // indicate that plugins should be initialized for that task
    fFlags|=kDigitsListInit; // indicate that it was initialized
  }
  if (fESDsQAList!=NULL && (fFlags&kESDsListInit)==0) {
    init|=kESDsListInit; // indicate that plugins should be initialized for that task
    fFlags|=kESDsListInit; // indicate that it was initialized
  }
  if (fRawsQAList!=NULL && (fFlags&kRawsListInit)==0) {
    init|=kRawsListInit; // indicate that plugins should be initialized for that task
    fFlags|=kRawsListInit; // indicate that it was initialized
  }
  if (fRecPointsQAList!=NULL && (fFlags&kRecPointsListInit)==0) {
    init|=kRecPointsListInit; // indicate that plugins should be initialized for that task
    fFlags|=kRecPointsListInit; // indicate that it was initialized
  }

  TIter next(&fPlugins);
  TObject* obj=NULL;
  while ((obj=next())) {
    AliHLTQADataMakerBase* plugin=dynamic_cast<AliHLTQADataMakerBase*>(obj);
    if (!plugin) continue;
    // transfer the properties set in AliQAManager::GetQADataMaker to the plugin
    plugin->SetName(GetName());
    plugin->SetUniqueID(GetUniqueID());
    if (init&kDigitsListInit) plugin->Init(AliQAv1::GetTaskIndex("Digits"), 0);
    if (init&kESDsListInit) plugin->Init(AliQAv1::GetTaskIndex("ESDs"), 0);
    if (init&kRawsListInit) plugin->Init(AliQAv1::GetTaskIndex("Raws"), 0);
    if (init&kRecPointsListInit) plugin->Init(AliQAv1::GetTaskIndex("RecPoints"), 0);
    plugin->StartOfDetectorCycle();
  }
}

void AliHLTQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray** list)
{
  // see header file for class documentation
  TIter next(&fPlugins);
  TObject* obj=NULL;
  TDirectory* dirBackup=gDirectory;
  gDirectory=NULL;
  while ((obj=next())) {
    AliHLTQADataMakerBase* plugin=dynamic_cast<AliHLTQADataMakerBase*>(obj);
    if (!plugin) continue;
    plugin->SetEventSpecie(GetEventSpecie());

    TObjArray** pluginList=NULL;
    if (task==AliQAv1::kESDS) {
      pluginList=plugin->GetESDsQAList();
    } else if (task==AliQAv1::kRECPOINTS) {
      pluginList=plugin->GetRecPointsQAList();
    } else if (task==AliQAv1::kDIGITS) {
      pluginList=plugin->GetDigitsQAList();
    } else if (task==AliQAv1::kRAWS) {
      pluginList=plugin->GetRawsQAList();
    }

    if (pluginList) {
      plugin->EndOfDetectorCycle(task, pluginList);
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
	if (!pluginList[specie]) continue;
	TIter nextentry(pluginList[specie]);
	TObject* entry=NULL;
	while ((entry=nextentry())) {
	  AliInfo(Form("cloning histogram %s for specie %d", entry->GetName(), specie));
	  list[specie]->Add(entry->Clone());
	}
      }
    }
  }
  gDirectory=dirBackup;
}

void AliHLTQADataMakerRec::MakeRaws(AliRawReader * rawReader)
{
  // see header file for class documentation
  if (!rawReader) return;

  TIter next(&fPlugins);
  TObject* obj=NULL;
  while ((obj=next())) {
    AliHLTQADataMakerBase* plugin=dynamic_cast<AliHLTQADataMakerBase*>(obj);
    if (!plugin) continue;
    plugin->SetEventSpecie(GetEventSpecie());
    plugin->MakeRaws(rawReader);
  }
}

void AliHLTQADataMakerRec::MakeESDs(AliESDEvent * esd, AliESDEvent* hltesd)
{
  // HLT QA on ESDs
  TIter next(&fPlugins);
  TObject* obj=NULL;
  while ((obj=next())) {
    AliHLTQADataMakerBase* plugin=dynamic_cast<AliHLTQADataMakerBase*>(obj);
    if (!plugin) continue;
    plugin->SetEventSpecie(GetEventSpecie());
    plugin->MakeESDs(esd, hltesd);
  }
}
