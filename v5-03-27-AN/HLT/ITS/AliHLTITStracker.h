//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTITSTRACKER_H
#define ALIHLTITSTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          High Level Trigger ITS tracker
//     reads AliITSclusterV2 clusters and HLT tracks as input and creates 
//       AliITStrackV2 tracks which then can be stored to the HLT ESD.
//     The class inheritates all the main features of the off-line
//       AliITStrackerV2 with some changes in order to load the HLT tracks
//       instead of the off-line TPC tracks. For time performance reasons
//       the tracker is supposed to perform only one reconstruction pass
//       and eventually propagate the tracks back to TPC inner border.
//       For detals how to use it, see RunHLTITS.C macro.
//      
//           Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch 
//-------------------------------------------------------------------------

#include "AliITStrackerV2.h"

class AliESDEvent;
class AliHLTITStrack;

//-------------------------------------------------------------------------
class AliHLTITStracker : public AliITStrackerV2 {
public:
  AliHLTITStracker():AliITStrackerV2(){ fConstraint[0]=1; fConstraint[1]=0; }
  AliHLTITStracker(const Char_t *geom) : AliITStrackerV2(geom){ fConstraint[0]=1; fConstraint[1]=0; }

  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);

  ClassDef(AliHLTITStracker,1)   //HLT ITS tracker
};

typedef AliHLTITStracker AliL3ITStracker; // for backward compatibility

#endif
