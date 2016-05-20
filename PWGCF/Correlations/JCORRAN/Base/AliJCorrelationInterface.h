/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Interface that all correlation analysis must fulfill

//===========================================================
// AliJCorrelations.h
//
//   Jussi Viinikainen
//===========================================================

#ifndef ALIJCORRELATIONINTERFACE_H
#define ALIJCORRELATIONINTERFACE_H

#include "AliJConst.h"
#include "AliJBaseTrack.h"

using namespace std;

class AliJCorrelationInterface {
  
public:
  
  AliJCorrelationInterface();  // default constructor
  virtual ~AliJCorrelationInterface(){;} //destructor
  
  virtual void FillHisto(corrFillType cFTyp, fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2) = 0; // virtual histogram filler method needed in AliJEventPool.cxx

};

#endif






















