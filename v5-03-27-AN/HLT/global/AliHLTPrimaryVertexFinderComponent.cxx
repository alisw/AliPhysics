// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timur Pocheptsov <Timur.Pocheptsov@cern.ch>           *
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

/// @file   AliHLTPrimaryVertexFinderComponent.cxx
/// @author Timur Pocheptsov
/// @date   2010-12-26
/// @brief  Primary vertex finder component
///

#include <algorithm>
#include <cerrno>

#include "TString.h"
#include "TMath.h"

#include "AliHLTPrimaryVertexFinderComponent.h"
#include "AliExternalTrackParam.h"
#include "AliHLTDataTypes.h"

ClassImp(AliHLTPrimaryVertexFinderComponent)

const double AliHLTPrimaryVertexFinderComponent::fgDefaultDeviation = 4.;

//________________________________________________________________________
AliHLTPrimaryVertexFinderComponent::AliHLTPrimaryVertexFinderComponent()
                             : fPrimaryOutput(),
                               fPrimaryVtx(),
                               fFitTracksToVertex(true),
                               fConstrainedTrackDeviation(fgDefaultDeviation)
{
  //Default ctor.
}

//________________________________________________________________________
const char* AliHLTPrimaryVertexFinderComponent::GetComponentID()
{
  //Component's "name".
  return "PrimaryVertexFinder";
}

//________________________________________________________________________
void AliHLTPrimaryVertexFinderComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  //Input to the primary vertex finder:
  //a)tracks from ESD object;
  //b)hlt tracks (ITS)
  //c)hlt tracks (TPC)
  list.clear();

  list.push_back(kAliHLTDataTypeESDObject);
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginITS);
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginTPC);
}

//________________________________________________________________________
AliHLTComponentDataType AliHLTPrimaryVertexFinderComponent::GetOutputDataType()
{
  //Data type of output.
  return kAliHLTMultipleDataType;
}

//________________________________________________________________________
int AliHLTPrimaryVertexFinderComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
  //Types for outputs from primary vertex finder (for V0 finder).
  list.clear();

  //Indices of tracks, participating in a primary.
  list.push_back(kAliHLTDataTypePrimaryFinder | kAliHLTDataOriginOut);
  //KFVertex - primary vertex.
  list.push_back(kAliHLTDataTypeKFVertex | kAliHLTDataOriginOut);

  return list.size();
}

//________________________________________________________________________
void AliHLTPrimaryVertexFinderComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //These numbers are complete crap.
  constBase = 80000;
  inputMultiplier = 2.;
}

//________________________________________________________________________
AliHLTComponent* AliHLTPrimaryVertexFinderComponent::Spawn()
{
  //Create primary vertex finder componet.
  return new AliHLTPrimaryVertexFinderComponent;
}

//________________________________________________________________________
int AliHLTPrimaryVertexFinderComponent::DoInit(int /*argc*/, const char** /*argv*/)
{
  //Process options.
  //1. Default parameters.
  fFitTracksToVertex = true;
  fConstrainedTrackDeviation = fgDefaultDeviation;

  //2. Parameters from OCDB.
  TString cdbPath("HLT/ConfigHLT/");
  cdbPath += GetComponentID();

  //This part will be uncommented as soon as
  //OCDB object is added.
  /*
  int res = ConfigureFromCDBTObjString(cdbPath);

  if (res < 0)
    return res;

  //3. "Command line" parameters.
  if (argc)
    res = ConfigureFromArgumentString(argc, argv);

  return res;*/

  return 0;
}

//________________________________________________________________________
int AliHLTPrimaryVertexFinderComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  //Scan the name of option and its parameters from the list.
  //Return number of processed entries.
  //Possible arguments:
  //-fitTracksToVertex 1/0
  //-constrainedTrackDeviation value
  AliHLTUtility::CmdLineParser parser;
  parser.Add("-fitTrackToVertex", &fFitTracksToVertex);
  parser.Add("-constrainedTrackDeviation", &fConstrainedTrackDeviation);

  const int nParsed = parser.Parse(argc, argv, 0);
  if (nParsed < 0) {
    HLTError(parser.GetError().Data());
    return -EPROTO;
  }

  return nParsed;
}

//________________________________________________________________________
int AliHLTPrimaryVertexFinderComponent::DoDeinit()
{
  //Reset parameters to default.
  fFitTracksToVertex = true;
  fConstrainedTrackDeviation = fgDefaultDeviation;

  return 0;
}

//________________________________________________________________________
int AliHLTPrimaryVertexFinderComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
                                       AliHLTComponentTriggerData& /*trigData*/)
{
  //Find primary vertex.

  if (GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR))
    return 0;


  //Clean all previous track infos and output block.
  fTrackInfos.clear();
  fPrimaryOutput.clear();

  //Initialize KF package.
  AliKFParticle::SetField(GetBz());

  //The logic, how vertex finder reads input tracks,
  //is taken from the original global vertexer.
  //First, try to read tracks from ESD event.
  //Both positive and negative PID are 211 ("default").
  //Other hypotesis are not possible here.
  ReadESDTracks(211, 211);
  //If no good esd tracks or no esd at all:
  if (!fTrackInfos.size())
    ReadHLTTracks(kAliHLTDataTypeTrack | kAliHLTDataOriginITS, 211, 211);
  //If no good its tracks:
  if (!fTrackInfos.size())
    ReadHLTTracks(kAliHLTDataTypeTrack | kAliHLTDataOriginTPC, 211, 211);

  if (!fTrackInfos.size()) {
    HLTWarning("No input tracks found");
    return 0;
  }

  FindPrimaryVertex();

  return DoOutput();
}

namespace
{

struct VertexDeviation
{
  int fI; //Index in fTrackInfos array.
  double fD; //Deviation from primary vertex.

  bool operator < (const VertexDeviation& rhs) const
  {
    return fD < rhs.fD;
  }
};

}

//________________________________________________________________________
void AliHLTPrimaryVertexFinderComponent::FindPrimaryVertex()
{
  //Find event's primary vertex.

  ///////////////////////////////////////////////////////////////////////
  //Some changes must be done here to read the initial guess (?)
  //for primary vertex.
  //Select rough region (in sigmas) in which the vertex could be found,
  //all tracks outside these limits are rejected from the primary vertex finding.
  fPrimaryVtx.Initialize();
  fPrimaryVtx.SetBeamConstraint(0., 0., 0., 3., 3., 5.3);
  ////////////////////////////////////////////////////////////////////////

  std::vector<const AliKFParticle*> vSelected(fTrackInfos.size());
  std::vector<VertexDeviation> devs(fTrackInfos.size());

  int nSelected = 0;
  for (VectorSize_t i = 0; i < fTrackInfos.size(); ++i) {
    const AliKFParticle& p = fTrackInfos[i].fParticle;
    const double chi = p.GetDeviationFromVertex(fPrimaryVtx);
    if (chi > fConstrainedTrackDeviation)
      continue;

    devs[nSelected].fI = i;
    devs[nSelected].fD = chi;
    vSelected[nSelected] = &fTrackInfos[i].fParticle;
    nSelected++;
  }

  //Fit
  while (nSelected > 2) {
    //Primary vertex finder with rejection of outliers
    for (int i = 0; i < nSelected; ++i)
      vSelected[i] = &fTrackInfos[devs[i].fI].fParticle;

    const double xv = fPrimaryVtx.GetX();
    const double yv = fPrimaryVtx.GetY();
    const double zv = fPrimaryVtx.GetZ(); //Values from the previous iteration.

    fPrimaryVtx.Initialize();
    fPrimaryVtx.SetBeamConstraint(0, 0, 0, 3., 3., 5.3);
    fPrimaryVtx.SetVtxGuess(xv, yv, zv);

    // refilled for every iteration
    //0: pointer to production vertex, -1. : mass, true : constrained.
    fPrimaryVtx.Construct(&vSelected[0], nSelected, 0, -1., true);

    for (int it = 0; it < nSelected; ++it) {
      const AliKFParticle& p = fTrackInfos[devs[it].fI].fParticle;
      if (nSelected <= 20) {
        //Exclude the current track from the sample and recalculate the vertex
        AliKFVertex tmp(fPrimaryVtx - p);
        devs[it].fD = p.GetDeviationFromVertex(tmp);
      } else {
        devs[it].fD = p.GetDeviationFromVertex(fPrimaryVtx);
      }
    }

    //Sort tracks with increasing chi2 (used for rejection)
    std::sort(&devs[0], &devs[0] + nSelected);

    //Remove 30% of the tracks (done for performance, only if there are more than 20 tracks)
    int nRemove = int(0.3 * nSelected);
    if (nSelected - nRemove <= 20)
      nRemove = 1;// removal based on the chi2 of every track

    int firstRemove = nSelected - nRemove;
    while (firstRemove < nSelected) {
      if (devs[firstRemove].fD >= fConstrainedTrackDeviation)
        break;
      firstRemove++;
    }

    if (firstRemove >= nSelected)
      break;

    nSelected = firstRemove;
  }

  if (nSelected < 3) {//No vertex for less than 3 contributors.
    fPrimaryVtx.NDF() = -3;
    fPrimaryVtx.Chi2() = 0.;
    nSelected = 0;
  }

  if (nSelected) {
    //Prepare output block.
    fPrimaryOutput.resize(sizeof(PrimaryFinderBlock) + sizeof(int) * (nSelected - 1));
    PrimaryFinderBlock* out = reinterpret_cast<PrimaryFinderBlock*>(&fPrimaryOutput[0]);

    out->fFitTracksFlag = fFitTracksToVertex;
    out->fNPrimaryTracks = nSelected;

    int minID = fTrackInfos[devs[0].fI].fID;
    int maxID = minID;
    for (int i = 0; i < nSelected; ++i) {
      const int id = fTrackInfos[devs[i].fI].fID;
      minID = TMath::Min(minID, id);
      maxID = TMath::Max(maxID, id);
      out->fPrimTrackIds[i] = id;
    }

    out->fMinPrimID = minID;
    out->fMaxPrimID = maxID;
  }
}

//________________________________________________________________________
int AliHLTPrimaryVertexFinderComponent::DoOutput()
{
  //Primary vertex finder output.
  if (!fPrimaryOutput.size()) {
    //Vertex not found.
    //Messages? return values?
    HLTWarning("No primary vertex was found");
    return 0;
  }

  //1. indices of primary tracks;
  //int - type of PushBack's parameter.
  const int iResult = PushBack(&fPrimaryOutput[0], fPrimaryOutput.size(),
                               kAliHLTDataTypePrimaryFinder | kAliHLTDataOriginOut);
  if (iResult < 0)
    return iResult;
  //2. primary vertex (AliKFVertex).
  return PushBack(&fPrimaryVtx, kAliHLTDataTypeKFVertex | kAliHLTDataOriginOut, 0);
}
