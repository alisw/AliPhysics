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

/// @file   AliHLTV0FinderComponent.cxx
/// @author Timur Pocheptsov
/// @date   2010-12-26
/// @brief  V0 finder component
///

#include "TString.h"

#include "AliHLTDataTypes.h"
#include "AliHLTV0FinderComponent.h"

ClassImp(AliHLTV0FinderComponent)

const double AliHLTV0FinderComponent::fgDaughterPrimDeviation = 2.5;
const double AliHLTV0FinderComponent::fgPrimDeviation = 3.5;
const double AliHLTV0FinderComponent::fgChi = 3.5;
const double AliHLTV0FinderComponent::fgDecayLengthInSigmas = 3.;

//________________________________________________________________________
AliHLTV0FinderComponent::AliHLTV0FinderComponent()
                  : fPrimaryVtx(),
                    fPrimaryTracks(),
                    fNPrimaryTracks(0),
                    fMinPrimID(0),
                    fMaxPrimID(0),
                    fDaughterPrimDeviation(fgDaughterPrimDeviation),
                    fPrimDeviation(fgPrimDeviation),
                    fChi(fgChi),
                    fDecayLengthInSigmas(fgDecayLengthInSigmas),
                    fPosPID(211),
                    fNegPID(211),
                    fV0s(),
                    fGammaFinder(false)
{
  //Default ctor.
}

//________________________________________________________________________
const char* AliHLTV0FinderComponent::GetComponentID()
{
  //Component's "name".
  return "V0Finder";
}

//________________________________________________________________________
void AliHLTV0FinderComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  //Input to the primary vertex finder:
  //a)tracks from ESD object;
  //b)hlt tracks (ITS)
  //c)hlt tracks (TPC)
  //d)primary vertex (AliKFVertex)
  //e)indices of primary tracks.
  list.clear();
  //Input tracks.
  list.push_back(kAliHLTDataTypeESDObject);
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginITS);
  list.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginTPC);
  //Input from primary finder:
  //Primary vertex.
  list.push_back(kAliHLTDataTypeKFVertex | kAliHLTDataOriginOut);
  //Primary tracks' indices.
  list.push_back(kAliHLTDataTypePrimaryFinder | kAliHLTDataOriginOut);
}

//________________________________________________________________________
AliHLTComponentDataType AliHLTV0FinderComponent::GetOutputDataType()
{
  //Data type of output.
  return kAliHLTMultipleDataType;
}

//________________________________________________________________________
int AliHLTV0FinderComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  //Output type for V0 finder.
  tgtList.clear();
  //Indices of tracks, participating in V0s.
  tgtList.push_back(kAliHLTDataTypeV0Finder | kAliHLTDataOriginOut);

  return tgtList.size();
}

//________________________________________________________________________
void AliHLTV0FinderComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //These numbers are complete crap.
  constBase = 80000;
  inputMultiplier = 2.;
}

//________________________________________________________________________
AliHLTComponent* AliHLTV0FinderComponent::Spawn()
{
  //Create primary vertex finder componet.
  return new AliHLTV0FinderComponent;
}

//________________________________________________________________________
int AliHLTV0FinderComponent::DoInit(int argc, const char** argv)
{
  //1. Default parameters.
  fDaughterPrimDeviation = fgDaughterPrimDeviation;
  fPrimDeviation = fgPrimDeviation;
  fChi = fgChi;
  fDecayLengthInSigmas = fgDecayLengthInSigmas;
  fPosPID = 211;
  fNegPID = 211;
  fGammaFinder = false;

  //2. Parameters from OCDB.
  TString cdbPath("HLT/ConfigHLT/");
  cdbPath += GetComponentID();

  int res = ConfigureFromCDBTObjString(cdbPath);
  if (res < 0)
    return res;

  //3. "Command line" parameters.
  if (argc)
    res = ConfigureFromArgumentString(argc, argv);

  fV0s.clear();
  fV0s.push_back(0); //Number of v0s.

  return res;
}

//________________________________________________________________________
int AliHLTV0FinderComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  //Scan one argument and its parameters from the list
  //Return number of processed entries.
  //Possible arguments:
  //-daughterPrimDeviation num
  //-primDeviation num
  //-chi num
  //-decayLengthInSigmas num
  //-posPID int_num
  //-negPid int_num
  //-gammaFinder 0/1

  AliHLTUtility::CmdLineParser parser;
  parser.Add("-daughterPrimDeviation", &fDaughterPrimDeviation);
  parser.Add("-primDeviation", &fPrimDeviation);
  parser.Add("-chi", &fChi);
  parser.Add("-decayLengthInSigmas", &fDecayLengthInSigmas);
  parser.Add("-posPID", &fPosPID);
  parser.Add("-negPID", &fNegPID);
  parser.Add("-gammaFinder", &fGammaFinder);

  const int nParsed = parser.Parse(argc, argv, 0);
  if (nParsed < 0) {
    HLTError(parser.GetError().Data());
    return -EPROTO;
  }

  return nParsed;
}

//________________________________________________________________________
int AliHLTV0FinderComponent::DoDeinit()
{
  //Reset parameters to default.
  fDaughterPrimDeviation = fgDaughterPrimDeviation;
  fPrimDeviation = fgPrimDeviation;
  fChi = fgChi;
  fDecayLengthInSigmas = fgDecayLengthInSigmas;
  fPosPID = 211;
  fNegPID = 211;
  fGammaFinder = false;

  return 0;
}

//________________________________________________________________________
int AliHLTV0FinderComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
                            AliHLTComponentTriggerData& /*trigData*/)
{
  //Find primary vertex.
  if (GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR))
    return 0;

  //Clean all previous track infos.
  fTrackInfos.clear();
  fPrimaryTracks.clear();

  //Initialize KF package.
  AliKFParticle::SetField(GetBz());

  // no output if there is no vertex information
  if (!ReadPrimaryVertex())
    return 0;

  // start from clean array
  // the output block contains at least the number of pairs even
  // if this is zero
  fV0s.clear();
  fV0s.push_back(0); //Number of v0s.

  if (ReadTracks()) {
    FindV0s();
  }

  return DoOutput();
}

//________________________________________________________________________
bool AliHLTV0FinderComponent::ReadPrimaryVertex()
{
  //Primary finder produces primary vertex (AliKFVertex) and
  //indices for primary tracks.
  //1. Try to extract primary vertex.
  const TObject* obj = GetFirstInputObject(kAliHLTDataTypeKFVertex | kAliHLTDataOriginOut);
  const AliKFVertex* kfVtx = dynamic_cast<const AliKFVertex*>(obj);

  if (!kfVtx) {
    HLTError("V0 finder requires KF vertex (primary vertex) as input");
    return false;
  }

  //2. Try to read primary track indices.
  const AliHLTComponentBlockData* p = GetFirstInputBlock(kAliHLTDataTypePrimaryFinder
                                                          | kAliHLTDataOriginOut);
  if (!p || !p->fSize || !p->fPtr) {
    HLTError("Array of primary track indices expected");
    return false;
  }

  //Data block from primary finder
  const PrimaryFinderBlock* blk = static_cast<PrimaryFinderBlock*>(p->fPtr);
  //Track ids must be positive integers.
  if (blk->fMinPrimID < 0 || blk->fMaxPrimID < 0) {
    HLTError("Got negative track ID from primary finder, internal HLT error");
    return false;
  }

  //3. Got correct data, modify the component's state.
  //KF vertex.
  fPrimaryVtx = *kfVtx;
  //Primary tracks.
  fNPrimaryTracks = blk->fNPrimaryTracks;
  fMinPrimID = blk->fMinPrimID;
  fMaxPrimID = blk->fMaxPrimID;

  fPrimaryTracks.assign(fMaxPrimID + 1, 0);
  for (int i = 0; i < fNPrimaryTracks; ++i)
    fPrimaryTracks[blk->fPrimTrackIds[i]] = 1;

  return true;
}

//________________________________________________________________________
bool AliHLTV0FinderComponent::ReadTracks()
{
  //The logic, how vertex finder reads input tracks,
  //is taken from the original global vertexer.
  //First, try to read tracks from ESD event.
  ReadESDTracks(fPosPID, fNegPID);
  if (fTrackInfos.size()) {
    FindPrimaryDeviations();
    return true;
  }

  //No good esd tracks, try:
  ReadHLTTracks(kAliHLTDataTypeTrack | kAliHLTDataOriginITS, fPosPID, fNegPID);
  if (fTrackInfos.size()) {
    FindPrimaryDeviations();
    return true;
  }

  //If no good its tracks, try:
  ReadHLTTracks(kAliHLTDataTypeTrack | kAliHLTDataOriginTPC, fPosPID, fNegPID);
  if (fTrackInfos.size()) {
    FindPrimaryDeviations();
    return true;
  }

  HLTError("No input tracks found for V0 finder");

  return false;
}

//________________________________________________________________________
void AliHLTV0FinderComponent::FindPrimaryDeviations()
{
  //Quite a tricky part.
  for (VectorSize_t i = 0; i < fTrackInfos.size(); ++i) {
    AliHLTTrackInfo& info = fTrackInfos[i];
    if (IsPrimaryTrack(info.fID)) {
      info.fPrimUsed = true;
      //The way primary deviation is computed in primary finder:
      if (fNPrimaryTracks <= 20) {
        AliKFVertex tmp(fPrimaryVtx - info.fParticle);
        info.fPrimDeviation = info.fParticle.GetDeviationFromVertex(tmp);
      } else
        info.fPrimDeviation = info.fParticle.GetDeviationFromVertex(fPrimaryVtx);
    } else {
      info.fPrimUsed = false;
      info.fPrimDeviation = info.fParticle.GetDeviationFromVertex(fPrimaryVtx);
    }
  }
}

//________________________________________________________________________
bool AliHLTV0FinderComponent::IsPrimaryTrack(int id)const
{
  if (id < fMinPrimID || id > fMaxPrimID)
    return false;

  return fPrimaryTracks[id];
}

//________________________________________________________________________
void AliHLTV0FinderComponent::FindV0s()
{
  //Here's the core.
  if (fPrimaryVtx.GetNContributors() < 3)
    return;

  for (int iTr = 0, ei = fTrackInfos.size(); iTr < ei; ++iTr) {
    AliHLTTrackInfo& info = fTrackInfos[iTr];
    if (info.fParticle.GetQ() > 0)
      continue;
    if (info.fPrimDeviation < fDaughterPrimDeviation)
      continue;

    for (int jTr = 0; jTr < ei; ++jTr) {
      AliHLTTrackInfo& jnfo = fTrackInfos[jTr];
      if (jnfo.fParticle.GetQ() < 0)
        continue;
      if (jnfo.fPrimDeviation < fDaughterPrimDeviation)
        continue;

      //Check if the particles fit
      if (info.fParticle.GetDeviationFromParticle(jnfo.fParticle) > fChi)
        continue;

      //Construct V0 mother
      AliKFParticle v0(info.fParticle, jnfo.fParticle /*, bGammaFinder*/);
      //Check V0 Chi^2
      if (v0.GetChi2() < 0. || v0.GetChi2() > fChi * fChi * v0.GetNDF())
        continue;

      //Subtruct daughters from primary vertex
      AliKFVertex primVtxCopy(fPrimaryVtx);

      if (info.fPrimUsed) {
        if (primVtxCopy.GetNContributors() <= 2)
          continue;
        primVtxCopy -= info.fParticle;
      }

      if (jnfo.fPrimUsed) {
        if (primVtxCopy.GetNContributors() <= 2)
          continue;
        primVtxCopy -= jnfo.fParticle;
      }

      //Check v0 Chi^2 deviation from primary vertex
      if (v0.GetDeviationFromVertex(primVtxCopy) > fPrimDeviation)
        continue;
      //Add V0 to primary vertex to improve the primary vertex resolution
      primVtxCopy += v0;
      //Set production vertex for V0
      v0.SetProductionVertex(primVtxCopy);
      //Get V0 decay length with estimated error
      double length = 0., sigmaLength = 0.;
      if (v0.GetDecayLength(length, sigmaLength))
        continue;
      //Reject V0 if it decays too close[sigma] to the primary vertex
      if (length  < fDecayLengthInSigmas * sigmaLength)
        continue;
      //Keep v0
      fV0s.push_back(info.fID);
      fV0s.push_back(jnfo.fID);

      fV0s[0] += 1;
    }
  }
}

//________________________________________________________________________
int AliHLTV0FinderComponent::DoOutput()
{
  //Save V0s' track pairs' indices.
  if (!fV0s.size()) {
    HLTError("internal data mismatch of output structure");
    return 0;
  }

  HLTInfo("Number of v0 candidates: %d", fV0s[0]);

  //Indices of primary tracks.
  return PushBack(&fV0s[0], fV0s.size() * sizeof(int),
                  kAliHLTDataTypeV0Finder | kAliHLTDataOriginOut);
}
