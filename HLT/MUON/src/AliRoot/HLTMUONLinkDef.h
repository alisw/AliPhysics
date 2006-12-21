/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedefs;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ nestedfunction;

#pragma link C++ function AliHLTMUONVersion;
#pragma link C++ function AliHLTMUONMajorVersion;
#pragma link C++ function AliHLTMUONMinorVersion;
#pragma link C++ function AliHLTMUONBuildNumber;

#pragma link C++ class AliHLTMUONRegion+;
#pragma link C++ class AliHLTMUONPoint+;
#pragma link C++ class AliHLTMUONTriggerRecord+;
#pragma link C++ class AliHLTMUONADCStream+;
#pragma link C++ class AliHLTMUONTrack+;

#pragma link C++ class AliHLTMUONMicrodHLT+;

#pragma link C++ class AliHLTMUONADCStreamSource+;
#pragma link C++ struct AliHLTMUONADCStreamSource::AliDataBlock+;

#pragma link C++ class AliHLTMUONTriggerSource+;
#pragma link C++ enum AliHLTMUONTriggerSource::AreaType+;
#pragma link C++ enum AliHLTMUONTriggerSource::SourceType+;
#pragma link C++ class AliHLTMUONTriggerSource::AliEventData+;

#pragma link C++ class AliHLTMUONClusterSource+;
#pragma link C++ enum AliHLTMUONClusterSource::AreaType+;
#pragma link C++ enum AliHLTMUONClusterSource::SourceType+;
#pragma link C++ class AliHLTMUONClusterSource::AliBlockData+;
#pragma link C++ class AliHLTMUONClusterSource::AliEventData+;

#pragma link C++ class AliHLTMUONTrackSink+;
#pragma link C++ class AliHLTMUONTrackSink::AliEventData+;

#pragma link C++ class AliHLTMUONTrackerInterface+;
#pragma link C++ class AliHLTMUONTrackerCallback+;

#pragma link C++ class AliHLTMUONClusterFinderInterface+;
#pragma link C++ class AliHLTMUONClusterFinderCallback+;

#pragma link C++ class AliHLTMUONTracker+;
#pragma link C++ class AliHLTMUONHitReconstructor+;
#endif // __CINT__

