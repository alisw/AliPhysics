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

/// @file   AliHLTGlobalHistoComponent.cxx
/// @author Matthias Richter
/// @date   2010-09-16
/// @brief  A histogramming component for global ESD properties based
///         on the AliHLTTTreeProcessor

#include "AliHLTGlobalHistoComponent.h"
#include "AliESDEvent.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalHistoComponent)

AliHLTGlobalHistoComponent::AliHLTGlobalHistoComponent()
  : AliHLTTTreeProcessor()
  , fEvent(0)
  , fNofTracks(0)
  , fTrackVariables()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTGlobalHistoComponent::~AliHLTGlobalHistoComponent()
{
  // see header file for class documentation
  fTrackVariables.Reset();
}

void AliHLTGlobalHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTAllDataTypes);
}

TTree* AliHLTGlobalHistoComponent::CreateTree(int /*argc*/, const char** /*argv*/)
{
  // create the tree and branches
  int iResult=0;
  TTree* pTree=new TTree("ESDproperties", "HLT ESD properties");
  if (!pTree) return NULL;

  const char* trackVariableNames = {
    "Track_pt "
    "Track_phi "
    "Track_eta "
    // "Track_dcar "
    // "Track_dcaz "
  };
  int maxTrackCount=20000; // FIXME: make configurable
  if ((iResult=fTrackVariables.Init(maxTrackCount, trackVariableNames))<0) {
    HLTError("failed to initialize internal structure for track properties");
  }
  
  if (iResult>=0) {
    pTree->Branch("event",        &fEvent, "event/I");
    pTree->Branch("trackcount",   &fNofTracks, "trackcount/I");
    for (int i=0; i<fTrackVariables.Variables(); i++) {
      TString specifier=fTrackVariables.GetKey(i);
      float* pArray=fTrackVariables.GetArray(specifier);
      specifier+="[trackcount]/f";
      pTree->Branch(fTrackVariables.GetKey(i), pArray, specifier.Data());
    }
  } else {
    delete pTree;
    pTree=NULL;
  }

  return pTree;
}

void AliHLTGlobalHistoComponent::FillHistogramDefinitions()
{
  /// default histogram definitions
}

int AliHLTGlobalHistoComponent::FillTree(TTree* pTree, const AliHLTComponentEventData& /*evtData*/, 
					 AliHLTComponentTriggerData& /*trigData*/ )
{
  /// fill the tree from the ESD
  int iResult=0;
  if (!IsDataEvent()) return 0;

  ResetVariables();

  // fetch ESD from input stream
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  esd->GetStdContent();

  // fill track variables
  fNofTracks=esd->GetNumberOfTracks();
  for (int i=0; i<fNofTracks; i++) {
    AliESDtrack *esdTrack = esd->GetTrack(i);
    if (!esdTrack) continue;
    
    fTrackVariables.Fill("Track_pt", esdTrack->Pt());
    fTrackVariables.Fill("Track_phi", esdTrack->Phi());
    fTrackVariables.Fill("Track_eta", esdTrack->Theta());
  }
  HLTInfo("added parameters for %d tracks", fNofTracks);

  if (iResult<0) {
    // fill an empty event
    ResetVariables();
  }
  fEvent++;

  pTree->Fill();
  return iResult;
}

int AliHLTGlobalHistoComponent::ResetVariables()
{
  /// reset all filling variables
  fNofTracks=0;
  fTrackVariables.ResetCount();
  return 0;
}

AliHLTComponentDataType AliHLTGlobalHistoComponent::GetOriginDataType() const
{
  // get the origin of the output data
  return kAliHLTVoidDataType;
}

// AliHLTUInt32_t AliHLTGlobalHistoComponent::GetDataSpec() const;
// {
//   // get specifications of the output data
//   return 0;
// }

AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::AliHLTGlobalHistoVariables()
 : fCapacity(0)
 , fArrays()
 , fCount()
 , fKeys()
{
  /// default constructor
}

AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::~AliHLTGlobalHistoVariables()
{
  /// destructor
  Reset();
}

int AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::Init(int capacity, const char* names)
{
  /// init the arrays
  int iResult=0;
  TString initializer(names);
  TObjArray* pTokens=initializer.Tokenize(" ");
  fCapacity=capacity;
  if (pTokens) {
    int entries=pTokens->GetEntriesFast();
    fArrays.resize(entries);
    fCount.resize(entries);
    for (int i=0; i<entries; i++) {
      fKeys[pTokens->At(i)->GetName()]=i;
      fArrays[i]=new float[fCapacity];
    }
    delete pTokens;
  }
  assert(fArrays.size()==fCount.size());
  assert(fArrays.size()==fKeys.size());
  if (fArrays.size()!=fCount.size() ||
      fArrays.size()!=fKeys.size()) {
    return -EFAULT;
  }

  ResetCount();
  return iResult;
}

int AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::Reset()
{
  /// reset the arrays
  for (vector<float*>::iterator element=fArrays.begin();
       element!=fArrays.end();
       element++) {
    delete *element;
  }
  fArrays.clear();
  fCount.clear();
  fKeys.clear();

  return 0;
}

int AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::ResetCount()
{
  /// reset the fill counts
  for (vector<int>::iterator element=fCount.begin();
       element!=fCount.end();
       element++) {
    *element=0;
  }
  return 0;
}

int AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::Fill(unsigned index, float value)
{
  /// fill variable at index
  assert(fArrays.size()==fCount.size());
  if (index>=fArrays.size() || index>=fCount.size()) return -ENOENT;
  if (fCount[index]>=fCapacity) return -ENOSPC;

  (fArrays[index])[fCount[index]++]=value;
  return fCount[index];
}

int AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::Fill(const char* key, float value)
{
  /// fill variable at key
  int index=FindKey(key);
  if (index<0) return -ENOENT;

  return Fill(index, value);
}

float* AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::GetArray(const char* key)
{
  /// get array at key
  int index=FindKey(key);
  if (index<0) return NULL;
  assert((unsigned)index<fArrays.size());
  if ((unsigned)index>=fArrays.size()) return NULL;

  return fArrays[index];
}

int AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::FindKey(const char* key) const
{
  /// fill variable at key
  map<string, int>::const_iterator element=fKeys.find(key);
  if (element==fKeys.end()) return -ENOENT;
  return element->second;
}

const char* AliHLTGlobalHistoComponent::AliHLTGlobalHistoVariables::GetKey(int index) const
{
  /// fill variable at key
  for (map<string, int>::const_iterator element=fKeys.begin();
       element!=fKeys.end(); element++) {
    if (element->second==index) return element->first.c_str();
  }
  return NULL;
}
