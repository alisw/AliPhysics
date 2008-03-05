//-*- Mode: C++ -*-
// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTTPCRunStatistics.cxx
    @author Jochen Thaeder
    @date   
    @brief  TPC class for run statistics, derived from @see AliHLTRunStatistics
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCRunStatistics.h"

ClassImp(AliHLTTPCRunStatistics)
    
  AliHLTTPCRunStatistics::AliHLTTPCRunStatistics() :
    fNTotalTracks(0),
    fNTracksAboveClusterThreshold(0),
    fClusterThreshold(0),
    fNTotalCluster(0),
    fNUsedCluster(0),
    fAvgClusterPerTrack(0) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCRunStatistics::~AliHLTTPCRunStatistics() {
  // see header file for class documentation
}


