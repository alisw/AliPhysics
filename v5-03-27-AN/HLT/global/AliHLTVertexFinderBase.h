//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTVERTEXFINDERBASE_H
#define ALIHLTVERTEXFINDERBASE_H
/// This file is property of and copyright by the ALICE HLT Project         
/// ALICE Experiment at CERN, All rights reserved.                         
/// See cxx source for full Copyright notice                               

/// @file   AliHLTVertexFinderBase.h
/// @author Timur Pocheptsov
/// @date   2010-12-26
/// @brief  Base class for vertex finder components
///


#include <vector>
#include <map>

#include "TString.h"

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
#include "AliKFParticle.h"


//Auxiliary base class, used by primary and V0 vertex finders.
//It defines a nested type AliHLTTrackInfo and some
//functions to read tracks from input blocks -
//this is usefull for both vertexers.
//This class has to inherit AliHLTProcessor to be able
//to call member functions from AliHLTComponent.
//This class is not for end user, it's just a
//common base for vertex finders.

class AliESDEvent;
class AliKFVertex;

class AliHLTVertexFinderBase : public AliHLTProcessor
{
protected:

  class AliHLTTrackInfo
  {
  public:
    AliHLTTrackInfo()
        : fID(0),
          fParticle(),
          fPrimUsed(false),
          fPrimDeviation(0.)
    {
    }

    AliHLTTrackInfo(int id, const AliVTrack& track,
                    int pid, bool primUsed, double dev)
        : fID(id),
          fParticle(track, pid),
          fPrimUsed(primUsed),
          fPrimDeviation(dev)
    {
    }

    //Check AliKFParticle parameters and covariance matrix.
    bool IsKFFinite()const;

    //Track ID taken from input track, just an index in some "collection".
    int fID;

    AliKFParticle fParticle; //Corresponding to track KFParticle.
    bool fPrimUsed; //Particle was used for primary vertex fit.
    double fPrimDeviation; //Primary deviation.
  };

  //These data types are
  //used in FillESD functions,
  //they also used in vertex finders
  //(v0 finder uses primary finder's block).
  struct PrimaryFinderBlock
  {
    int fFitTracksFlag;
    int fNPrimaryTracks;
    int fMinPrimID;
    int fMaxPrimID;
    int fPrimTrackIds[1];
  };

  struct V0FinderBlock
  {
    int fNV0s;
    int fV0s[1];
  };

protected:

  //Compiler-generated def ctor will call vector's def
  //ctor for member, but still, compiler wants mem-init-list.
  AliHLTVertexFinderBase() : fTrackInfos()
  {
  }

  //Read tracks from AliESDEvent
  //into fTrackInfos.
  //IMPORTANT: fTrackInfos is not
  //cleared here before filling,
  //tracks are appended to the end.
  //A user (derived class)
  //has to clear it, if needed.
  void ReadESDTracks(int posPID, int negPID);
  void ReadESDTracks(AliESDEvent* esd, int posPID, int negPID);
  //Read HLT tracks produced by ITS
  //or TPC (somewhere else?).
  //IMPORTANT: fTrackInfos is not
  //cleared here before filling,
  //tracks are appended to the end.
  void ReadHLTTracks(const AliHLTComponentDataType& blockType, int posPID, int negPID);

public:
  //Produce output for ESD collector from primary and v0 finders' outputs.
  static void FillESD(AliESDEvent* esd, AliKFVertex* primVtx, const void* primData,
                      const void* v0Data);

protected:

  typedef std::vector<AliHLTTrackInfo> TrackInfoVector_t;
  typedef TrackInfoVector_t::size_type VectorSize_t;

  TrackInfoVector_t fTrackInfos;

  ClassDef(AliHLTVertexFinderBase, 0);
};

namespace AliHLTUtility
{

//Small utility to automate command line parsing.
//Command line consists of commands and parameters.
//Currently, parameters can be of type int, double,
//bool (so, command with one parameter of built-in type)
//or CompoundType (this is for more complex
//expressions like -command param1 param2 param3).
//If SetParameter failes, it MUST throw std::runtime_error
//with description. CmdLineParser will convert this
//exception into negative integer and will compose
//error message.

class CompoundType
{
public:
  virtual ~CompoundType()
  {
  }

  virtual unsigned SetParameter(unsigned argc, const char** argv, unsigned currPos) = 0;
};

//Parameter for a command. If the type is CompoundType, it's up to CompoundType
//to do all parsing. If it's int, double, bool - checks are done here:
//1. Check if parameter was set already (user specified the same command more than
//   once - this is an error).
//2. Check if parameter for command was specified.
//3. Parameter has correct value: 0 or 1 for bool and
//   for int or double non-empty constraint is satisfied
//
//Object of type Parameter holds a pointer to real object of type int, bool, double
//(or CompoundType) and this real objects are updated
//from command line parameters, using SetParameter function.

class Parameter
{
public:
  enum Constraint
  {
    none,
    positive,
    nonNegative
    //Something else.
  };

  Parameter();
  Parameter(bool* b);
  Parameter(int* i, Constraint c);
  Parameter(double* d, Constraint c);
  Parameter(CompoundType* cp);

  //Returns number of argv elements processed, starting from
  //currPos.
  unsigned SetParameter(unsigned argc, const char** argv, unsigned currPos);

private:
  //If any constraint was set:
  void CheckConstraint();
  //
  void SetConstraintChecker();

  bool fWasSet; //Parameter was set already.
  Constraint fConstraint; //Constraint on parameter.

  bool* fBool; //Real object to set of type bool.
  int* fInt; //Real object to set of type int.
  double* fDouble; //Real object to set of type double.
  CompoundType* fCompound; //"Object" to set of complex type.

  void (*fConstraintChecker)(const Parameter& param); //Constraint checker.

  static void CheckPositive(const Parameter& param);
  static void CheckNonNegative(const Parameter& param);
  //Something else.

  //Compiler generated dtor, copy-ctor and copy-assignment operators are ok.
};

//CmdLineParser parses char** argv, which constains commands and
//options for these commands.
//Returns the number of processed argv elements.

class CmdLineParser
{
public:

  CmdLineParser() : fParameters(), fIgnored(), fError()
  {
  }

  void Add(const TString& cmd, bool* b);
  void Add(const TString& cmd, int* i, Parameter::Constraint c = Parameter::none);
  void Add(const TString& cmd, double* d, Parameter::Constraint c = Parameter::none);
  void Add(const TString& cmd, CompoundType* ct);

  //Skip -command nParams * params
  void IgnoreCommand(const TString& cmd, unsigned nParams);

  int Parse(unsigned argc, const char** argv, unsigned currPos);

  const TString& GetError()const
  {
    return fError;
  }

private:

  typedef std::map<TString, Parameter> ParameterMap_t;
  typedef ParameterMap_t::iterator MapIter_t;
  typedef std::map<TString, unsigned> SkipMap_t;
  typedef SkipMap_t::const_iterator SkipMapIter_t;

  //Map of commands and it's parameters.
  ParameterMap_t fParameters;
  //These commands and its parameters (if any) will be skipped.
  SkipMap_t fIgnored;

  TString fError; //Error message describing parsing errors.
};

}

#endif
