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

/// @file   AliHLTVertexFinderBase.cxx
/// @author Timur Pocheptsov
/// @date   2010-12-26
/// @brief  Base class for vertex finder components
///

#include <stdexcept>
#include <cmath>

//AliHLTExternalTrackParam uses typedes from Rtypes.h
//but does not include it.
#include "Rtypes.h"

#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTVertexFinderBase.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliKFVertex.h"
#include "AliESDv0.h"

ClassImp(AliHLTVertexFinderBase)

//_________________________________________________________________
bool AliHLTVertexFinderBase::AliHLTTrackInfo::IsKFFinite()const
{
  //Check KF particle.

  //In C99, 'isfinite' is a macro.
  //But still, I add using directive, in case
  //it's a function in std.
  using namespace std;

  for (int i = 0; i < 8; ++i)
    if (!isfinite(fParticle.GetParameter(i)))
      return false;
  for (int i = 0; i < 36; ++i)//Not sure, probably, 27 is enough.
    if (!isfinite(fParticle.GetCovariance(i)))
      return false;

  return true;
}

//_________________________________________________________________
void AliHLTVertexFinderBase::ReadESDTracks(int posPID, int negPID)
{
  //Try to read input tracks from AliESDEvent.
  //Use esd track's index as "track id" (fID).
  //IMPORTANT: fTrackInfos is _NOT_ cleared here,
  //must be done externally (if you need to).

  const TObject* iter = GetFirstInputObject(kAliHLTDataTypeESDObject);
  for (; iter; iter = GetNextInputObject()) {
    //Cast away constness before dynamic_cast.
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>((TObject*)iter);
    if (!esd)//From original code. Not sure, if this is possible at all.
      continue;

    ReadESDTracks(esd, negPID, posPID);
    //Only one esd event is taken.
    break;
  }
}

//_________________________________________________________________
void AliHLTVertexFinderBase::ReadESDTracks(AliESDEvent* esd, int posPID, int negPID)
{
  //Try to read input tracks from AliESDEvent.
  //Use esd track's index as "track id" (fID).
  //IMPORTANT: fTrackInfos is _NOT_ cleared here,
  //must be done externally (if you need to).
  esd->GetStdContent();

  const int nTracks = esd->GetNumberOfTracks();
  if (nTracks)
    fTrackInfos.reserve(nTracks + fTrackInfos.size());

  for (int i = 0; i < nTracks; ++i) {
    AliESDtrack* pTrack = esd->GetTrack(i);
    //This checks of track parameters are
    //from the original global vertexer.
    if (!pTrack)
      continue;
    if (pTrack->GetKinkIndex(0) > 0)
      continue;
    if (!(pTrack->GetStatus() & AliESDtrack::kTPCin))
      continue;

    //i: track id, 211: pid, false: not used in primary, 0.: prim. deviation.
    AliHLTTrackInfo newTrackInfo(i, *pTrack, pTrack->Charge() > 0 ? posPID : negPID, false, 0.);
    if (!newTrackInfo.IsKFFinite())
      continue;

    fTrackInfos.push_back(newTrackInfo);
  }
}

//________________________________________________________________________
void AliHLTVertexFinderBase::ReadHLTTracks(const AliHLTComponentDataType& blockType,
                                           int posPID, int negPID)
{
  //Read HLT tracks as input.
  //IMPORTANT: fTrackInfos is _NOT_ cleared here,
  //must be done externally (if you need it).

  const AliHLTComponentBlockData* block = GetFirstInputBlock(blockType);
  for (; block; block = GetNextInputBlock()) {
    const AliHLTTracksData* dataPtr = static_cast<AliHLTTracksData*>(block->fPtr);
    const AliHLTExternalTrackParam* hltTrk = dataPtr->fTracklets;

    const int nTracks = dataPtr->fCount;
    if (nTracks)
      fTrackInfos.reserve(fTrackInfos.size() + nTracks);

    for (int i = 0; i < nTracks; ++i) {
      //Ugly conversion from one track format to another
      //and to AliKFParitcle then, to be changed later.
      AliHLTGlobalBarrelTrack tmpTrk(*hltTrk);

      AliHLTTrackInfo newTrackInfo(hltTrk->fTrackID, tmpTrk,
                                   tmpTrk.Charge() > 0 ? posPID : negPID, false, 0.);
      if (!newTrackInfo.IsKFFinite())
        continue;

      fTrackInfos.push_back(newTrackInfo);

      const unsigned dSize = sizeof(AliHLTExternalTrackParam)
                            + hltTrk->fNPoints * sizeof(unsigned);
      hltTrk = (AliHLTExternalTrackParam*)((char *)hltTrk + dSize);
    }
  }
}

//________________________________________________________________________
void AliHLTVertexFinderBase::FillESD(AliESDEvent* esd, AliKFVertex* vtx,
                                     const void* primData, const void* v0Data)
{
  //Put the output of a primary finder and v0 finder to the esd event.
  //Code was taken from the original global vertexer.
  //Code is an absolute mess.

  typedef std::map<int, int> Map_t;
  typedef Map_t::const_iterator MapIter_t;

  //Map esdId -> esdTrackIndex.
  Map_t mapId;

  double params[3];
  for (int i = 0; i < 3; ++i)
    params[i] = vtx->Parameters()[i];
  double cov[6];
  for (int i = 0; i < 6; ++i)
    cov[i] = vtx->CovarianceMatrix()[i];

  AliESDVertex vESD(params, cov, vtx->GetChi2(), vtx->GetNContributors());
  esd->SetPrimaryVertexTPC(&vESD);
  esd->SetPrimaryVertexTracks(&vESD);

  //Relate tracks to the primary vertex
  const PrimaryFinderBlock* prim = static_cast<const PrimaryFinderBlock*>(primData);
  const int nESDTracks = esd->GetNumberOfTracks();
  std::vector<bool> constrainedToVtx(nESDTracks);

  if (prim->fFitTracksFlag) {
    for (int i = 0; i < nESDTracks; ++i) {
      if (!esd->GetTrack(i))
        continue;

      mapId[esd->GetTrack(i)->GetID()] = i;
    }

    for (int i = 0; i < prim->fNPrimaryTracks; ++i) {
      MapIter_t it = mapId.find(prim->fPrimTrackIds[i]);
      if (it == mapId.end())
        continue;
      const int itr = it->second;
      //100. is an argument for parameter maxd in AliESDtrack - cut on impact parameter.
      esd->GetTrack(itr)->RelateToVertex(&vESD, esd->GetMagneticField(), 100.);
      constrainedToVtx[itr] = true;
    }
  }

  //Add v0s.
  if (v0Data) {
  const int* v0s = static_cast<const int*>(v0Data);
  const int nV0s = v0s[0];
  ++v0s;
  for (int i = 0; i < nV0s; ++i) {
    MapIter_t it = mapId.find(v0s[2 * i]);
    if (it==mapId.end())
      continue;
    const int iTr = it->second;

    it = mapId.find(v0s[2 * i + 1]);
    if (it == mapId.end())
      continue;
    const int jTr = it->second;

    AliESDv0 v0(*esd->GetTrack(iTr), iTr, *esd->GetTrack(jTr), jTr);
    esd->AddV0(&v0);
    // relate the tracks to the vertex
    if (prim->fFitTracksFlag) {
      if (constrainedToVtx[iTr] || constrainedToVtx[jTr])
        continue;

      double pos[3];
      double sigma[3] = {.1, .1, .1};
      v0.XvYvZv(pos);
      AliESDVertex v0ESD(pos, sigma);
      esd->GetTrack(iTr)->RelateToVertex(&v0ESD, esd->GetMagneticField(), 100.);
      esd->GetTrack(jTr)->RelateToVertex(&v0ESD, esd->GetMagneticField(), 100.);
      constrainedToVtx[iTr] = true;
      constrainedToVtx[jTr] = true;
    }
  }
  }
}

namespace AliHLTUtility
{

//_______________________________________________________________________
Parameter::Parameter()
             : fWasSet(false),
               fConstraint(none),
               fBool(0),
               fInt(0),
               fDouble(0),
               fCompound(0),
               fConstraintChecker(0)
{
  //Default ctor.
}

//_______________________________________________________________________
Parameter::Parameter(bool* b)
             : fWasSet(false),
               fConstraint(none),
               fBool(b),
               fInt(0),
               fDouble(0),
               fCompound(0),
               fConstraintChecker(0)
{
  //Parameter of type bool.
}

//_______________________________________________________________________
Parameter::Parameter(int* i, Constraint c)
             : fWasSet(false),
               fConstraint(c),
               fBool(0),
               fInt(i),
               fDouble(0),
               fCompound(0),
               fConstraintChecker(0)
{
  //Parameter of type int.
  SetConstraintChecker();
}

//_______________________________________________________________________
Parameter::Parameter(double* d, Constraint c)
             : fWasSet(false),
               fConstraint(c),
               fBool(0),
               fInt(0),
               fDouble(d),
               fCompound(0),
               fConstraintChecker(0)
{
  //Parameter of type double.
  SetConstraintChecker();
}

//_______________________________________________________________________
Parameter::Parameter(CompoundType *ct)
             : fWasSet(false),
               fConstraint(none),
               fBool(0),
               fInt(0),
               fDouble(0),
               fCompound(ct),
               fConstraintChecker(0)
{
  //Parameter of more complex user-defined type.
  //All checks and conversions must be implemented
  //by user of CompoundType.
}

//_______________________________________________________________________
unsigned Parameter::SetParameter(unsigned argc, const char** argv, unsigned currPos)
{
  //Set parameter from command line tokens.

  //It's up to compound parameter to parse.
  if (fCompound)
    return fCompound->SetParameter(argc, argv, currPos);

  //Now, int, bool or double must be set from a string.
  //Primitive checks are done here.
  if (currPos == argc)
    throw std::runtime_error("value expected");
  if (fWasSet)
    throw std::runtime_error("parameter was set already");

  const TString val(argv[currPos]);
  if (!val.Length())
    throw std::runtime_error("value expected");

  if (fBool) {
    if (val.Length() != 1 || (val[0] != '0' && val[0] != '1'))
      throw std::runtime_error("expected 0 or 1 for bool parameter");
    *fBool = val[0] - '0';
  } else if (fInt)
    *fInt = val.Atoi();
  else
    *fDouble = val.Atof();

  if (fConstraintChecker)
    fConstraintChecker(*this);

  fWasSet = true;

  //For double, int and bool - only one element of argv is processed.
  return 1;
}

//_______________________________________________________________________
void Parameter::SetConstraintChecker()
{
  //Set function pointer.
  if (fConstraint == positive)
    fConstraintChecker = CheckPositive;
  else if (fConstraint == nonNegative)
    fConstraintChecker = CheckNonNegative;
  //No constraints for bool or CompoundType.
}

//Aux. functions to check constraints.
//_______________________________________________________________________
void Parameter::CheckPositive(const Parameter& param)
{
  //val > 0
  if (param.fInt) {
    if (*param.fInt <= 0)
      throw std::runtime_error("integer parameter must be positive");
  } else if (param.fDouble) {
     if (*param.fDouble <= 0.)//Hmmmm.
      throw std::runtime_error("double parameter must be positive");
  }
  //For bool or CompoundType there is no check.
}

//_______________________________________________________________________
void Parameter::CheckNonNegative(const Parameter& param)
{
  //val >= 0
  if (param.fInt) {
    if (*param.fInt < 0)
      throw std::runtime_error("integer parameter must be non-negative");
  } else if (param.fDouble) {
    if (*param.fDouble < 0.)
      throw std::runtime_error("double parameter must be non-negative");
  }
  //For bool or CompoundType there is no check.
}


//Command line arguments parser.

//_______________________________________________________________________
void CmdLineParser::Add(const TString& cmd, bool* b)
{
  //Command with bool parameter.
  if (fParameters.find(cmd) == fParameters.end())
    fParameters[cmd] = Parameter(b);
}

//_______________________________________________________________________
void CmdLineParser::Add(const TString& cmd, int* i, Parameter::Constraint c)
{
  //Command with int parameter, with constraint.
  if (fParameters.find(cmd) == fParameters.end())
    fParameters[cmd] = Parameter(i, c);
}

//_______________________________________________________________________
void CmdLineParser::Add(const TString& cmd, double* d, Parameter::Constraint c)
{
  //Coomand with a parameter of type double, possibly with a constraint.
  if (fParameters.find(cmd) == fParameters.end())
    fParameters[cmd] = Parameter(d, c);
}

//_______________________________________________________________________
void CmdLineParser::Add(const TString& cmd, CompoundType *ct)
{
  //Command with a parameter of compound type.
  if (fParameters.find(cmd) == fParameters.end())
    fParameters[cmd] = Parameter(ct);
}

//_______________________________________________________________________
void CmdLineParser::IgnoreCommand(const TString& cmd, unsigned nParams)
{
  //Command and number of parameters to ignore.
  //I do not check, if it's in map already.
  fIgnored[cmd] = nParams;
}

//_______________________________________________________________________
int CmdLineParser::Parse(unsigned argc, const char** argv, unsigned currPos)
{
  //Parse. Returns number of processed arguments.
  //Negative number - error (as it's in HLT)
  fError.Clear();

  unsigned nProcessed = 0;
  while (currPos < argc) {
    //Command.
    const TString command(argv[currPos]);
    //Check, if this command must be ignored.
    SkipMapIter_t skipCmd = fIgnored.find(command);
    if (skipCmd != fIgnored.end()) {
      //Command and its parameters must be skipped
      nProcessed += 1 + skipCmd->second;
      currPos += 1 + skipCmd->second;
      continue;
    }

    //Corresponding parameter.
    MapIter_t it = fParameters.find(command);
    if (it == fParameters.end()) {
      fError = "Unknown command: " + command;
      return -1;
    }

    ++currPos;
    ++nProcessed;

    try {
      const unsigned n = it->second.SetParameter(argc, argv, currPos);
      nProcessed += n;
      currPos += n;
    } catch (const std::runtime_error& e) {
      fError = "Command " + command + " error: " + e.what();
      return -1;
    }
  }

  return nProcessed;
}

}
