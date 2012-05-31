//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTPRIMARYVERTEXFINDERCOMPONENT_H
#define ALIHLTPRIMARYVERTEXFINDERCOMPONENT_H
/// This file is property of and copyright by the ALICE HLT Project         
/// ALICE Experiment at CERN, All rights reserved.                         
/// See cxx source for full Copyright notice                               

/// @file   AliHLTPrimaryVertexFinderComponent.h
/// @author Timur Pocheptsov
/// @date   2010-12-26
/// @brief  Primary vertex finder component
///

#include <vector>

#include "AliHLTVertexFinderBase.h"
#include "AliKFVertex.h"

//Primary vertex finder, developed by Sergey Gorbunov.
//Based on KF package.
//Produces primary vertex (AliKFVertex object) and
//indices of input tracks, which participates
//in primary construction.
//Can be configured with two options:
//-fitTracksToVertex 0/1
//-constrainedTrackDeviation value.

class AliHLTPrimaryVertexFinderComponent : public AliHLTVertexFinderBase
{
public:
  AliHLTPrimaryVertexFinderComponent();

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
  int DoEvent(const AliHLTComponentEventData& evtData,
              AliHLTComponentTriggerData& trigData);

private:
  //Aux. staff.
  void FindPrimaryVertex();
  int DoOutput();

  std::vector<char> fPrimaryOutput; //Ids of primary tracks.
  AliKFVertex fPrimaryVtx; //Reconstructed KF primary vertex.
  bool fFitTracksToVertex; //Flag to store vertex constrained track parameters
  double fConstrainedTrackDeviation; //Deviation of a track from prim.vtx <=cut

  static const double fgDefaultDeviation; //Default value for fConstrainedTrackDeviation.

  //Compiler generated dtor, copy-ctor and copy-assignment operators are ok.

  ClassDef(AliHLTPrimaryVertexFinderComponent, 0);
};

#endif
