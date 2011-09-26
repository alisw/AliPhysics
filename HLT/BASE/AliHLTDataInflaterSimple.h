//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTDATAINFLATERSIMPLE_H
#define ALIHLTDATAINFLATERSIMPLE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTDataInflaterSimple.h
/// @author Matthias Richter
/// @date   2011-09-01
/// @brief  Data inflater implementation for format of AliHLTDataDeflaterSimple
/// @note   

#include "AliHLTDataInflater.h"
#include "AliHLTDataDeflaterSimple.h"

class AliHLTDataInflaterSimple : public AliHLTDataInflater
{
public:
  /// standard constructor
  AliHLTDataInflaterSimple();
  /// destructor
  ~AliHLTDataInflaterSimple();

  /// add a parameter definition to the configuration, return reference id
  int AddParameterDefinition(const char* name, int bitLength, int reducedBitLength);

  /// overloaded from AliHLTDataInflater
  virtual bool NextValue(AliHLTUInt64_t& value, AliHLTUInt32_t& length);

protected:
private:
  /** copy constructor prohibited */
  AliHLTDataInflaterSimple(const AliHLTDataInflaterSimple&);
  /** assignment operator prohibited */
  AliHLTDataInflaterSimple& operator=(const AliHLTDataInflaterSimple&);

  /// parameter definitions
  vector<AliHLTDataDeflaterSimple::AliHLTDataDeflaterParameter> fParameterDefinitions; //!

  /// current parameter during reading
  int fCurrentParameter; //!

  ClassDef(AliHLTDataInflaterSimple, 0)
};

#endif //ALIHLTDATAINFLATERSIMPLE_H
