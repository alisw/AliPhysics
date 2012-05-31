//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTV0FINDERCOMPONENT_H
#define ALIHLTV0FINDERCOMPONENT_H
/// This file is property of and copyright by the ALICE HLT Project         
/// ALICE Experiment at CERN, All rights reserved.                         
/// See cxx source for full Copyright notice                               

/// @file   AliHLTV0FinderComponent.h
/// @author Timur Pocheptsov
/// @date   2010-12-26
/// @brief  V0 finder component
///

#include <vector>

#include "AliHLTVertexFinderBase.h"
#include "AliKFVertex.h"

class AliHLTV0FinderComponent : public AliHLTVertexFinderBase
{
public:
  AliHLTV0FinderComponent();

  //AliHLTComponent's final-overriders.
  const char* GetComponentID();
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn();

  int DoInit(int argc, const char** argv);
  int ScanConfigurationArgument(int argc, const char** argv);
  int DoDeinit();

  using AliHLTProcessor::DoEvent;
  int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

private:
  //Aux. staff.

  //Read AliKFVertex, produced by primary finder.
  bool ReadPrimaryVertex();
  //Read input tracks.
  bool ReadTracks();
  //
  void FindPrimaryDeviations();
  //If track (id) is from primary.
  bool IsPrimaryTrack(int id)const;
  //
  void FindV0s();
  //PushBack the output (v0s pairs' indices).
  int DoOutput();

  AliKFVertex fPrimaryVtx; //Primary vertex.
  std::vector<int> fPrimaryTracks; //"bit mask": 1 - track is prim, 0 - no.
  int fNPrimaryTracks; //Number of primary tracks.
  int fMinPrimID; //Min id of primary tracks.
  int fMaxPrimID; //Max id of primary tracks.

  //V0 finder cuts.
  double fDaughterPrimDeviation; //daughters deviation from prim vertex >= cut
  double fPrimDeviation; //v0 deviation from prim vertex <= cut
  double fChi; //v0 sqrt(chi^2/NDF) <= cut
  double fDecayLengthInSigmas; //v0 decay length/sigma_length >= cut
  int fPosPID; //Pid for a positive track, when constracting v0.
  int fNegPID; //Pid for a negative track, when constracting v0.

  //Output of V0 finder.
  std::vector<int> fV0s; //Indices of track pairs, participating in a V0s.

  //For gammas, special version of AliKFParitcle must be constructed.
  bool fGammaFinder;

  //defaults for cuts.
  static const double fgDaughterPrimDeviation;
  static const double fgPrimDeviation;
  static const double fgChi;
  static const double fgDecayLengthInSigmas;

  //Compiler generated dtor, copy-ctor and copy-assignment operators are ok.

  ClassDef(AliHLTV0FinderComponent, 0);
};

#endif
