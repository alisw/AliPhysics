// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   testAliHLTGlobalBarrelTrack.C
    @author Matthias Richter
    @date   2010-12-17
    @brief  Test macro/program for the AliHLTGlobalBarrelTrack
 */

#ifndef __CINT__
#include <ostream>
#include <vector>
#include "AliHLTGlobalBarrelTrack.h"
#include "AliExternalTrackParam.h"
#include "TDatime.h"
#include "TRandom.h"
#endif //__CINT__

/**
 * Get a random number in the given range.
 */
int GetRandom(int min, int max)
{
  if (max-min<2) return min;
  static TRandom rand;
  static bool seedSet=false;
  if (!seedSet) {
    TDatime dt;
    rand.SetSeed(dt.Get());
    seedSet=true;
  }
  return min+rand.Integer(max-min);
}

int testAliHLTGlobalBarrelTrack()
{
  // test the AliHLTGlobalBarrelTrack class

  // 1)
  // the 5 track parameters are named in the AliHLTExternalTrackParam
  // while AliExternalTrackParam just uses an array[5]
  // the members have the same order, fY is the first one
  AliHLTExternalTrackParam hltparam;
  memset(&hltparam, 0, sizeof(AliHLTExternalTrackParam));
  hltparam.fY=GetRandom(0, 100000);
  hltparam.fZ=GetRandom(0, 100000);
  hltparam.fSinPsi=GetRandom(0, 100000);
  hltparam.fTgl=GetRandom(0, 100000);
  hltparam.fq1Pt=GetRandom(0, 100000);

  AliHLTGlobalBarrelTrack t(hltparam);
  if (TMath::Abs(hltparam.fY-t.GetY())>=1 ||
      TMath::Abs(hltparam.fZ-t.GetZ())>=1 ||
      TMath::Abs(hltparam.fSinPsi-t.GetSnp())>=1 ||
      TMath::Abs(hltparam.fTgl-t.GetTgl())>=1 ||
      TMath::Abs(hltparam.fq1Pt-t.GetSigned1Pt())) {
	cerr << "parameter sequence of AliExternalTrackParam and AliHLTExternalTrackParam does not match" << endl;
	return -1;
      }
  
  return 0;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  if ((iResult=testAliHLTGlobalBarrelTrack())<0) {
  }
  return iResult;
}
