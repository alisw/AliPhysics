//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSMIPCLUSTER_H
#define ALIHLTPHOSMIPCLUSTER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSConstants.h"
#include "Rtypes.h"

using namespace PhosHLTConst;

class AliHLTPHOSDigit;

#define Z_RANGE 5
#define X_RANGE 5

class  AliHLTPHOSMipCluster
{
public:  
  AliHLTPHOSMipCluster();
  virtual ~AliHLTPHOSMipCluster();
  void AddDigit(AliHLTPHOSDigit *digit);
  Float_t GetCenterAmplitude();
  void ResetMipCluster();
  void SetCenterCoordinate(int z, int x);
  Int_t GetZ();
  Int_t GetX();
  Int_t GetEntries();
  Float_t Get3x3Sum();

  void PrintInfo();
  void Remap(); // maps the digits relative to the center corrdinate

private:
  void ResetDigit(int z, int x);
  void CopyDigit(AliHLTPHOSDigit *inDigit, AliHLTPHOSDigit *copiedDigit);
  int fZ;
  int fX;
  int fEntries;
  AliHLTPHOSDigit *fMipcluster[Z_RANGE][X_RANGE];
};

#endif
