//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTESDTRACK_H
#define ALIHLTESDTRACK_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTESDtrack.h
/// @author Matthias Richter
/// @date   2010-10-29
/// @brief  An AliESDtrack child class doing the conversion between 
///         AliHLTESDOptTrack and AliESDtrack
/// @note   

#include "AliESDtrack.h"

class AliHLTOnlineESDtrack;

/**
 * @class AliHLTESDtrack
 * @brief AliESDtrack child implementing the conversion from AliHLTOnlineESDtrack
 *        to AliESDtrack
 *
 * The class has no own members, only specific methods.
 */
class AliHLTESDtrack : public AliESDtrack {
 public:
  /// standard constructor
  AliHLTESDtrack();
  /// copy constructor
  AliHLTESDtrack(const AliHLTESDtrack& t);
  /// destructor
  virtual ~AliHLTESDtrack();

  AliHLTESDtrack& operator=(const AliHLTESDtrack& t);
  AliHLTESDtrack& operator=(const AliHLTOnlineESDtrack& t);

private:

  ClassDef(AliHLTESDtrack, 1); // AliESDtrack child for AliHLTESDEvent
};
#endif
