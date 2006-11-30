#ifndef ALIL3HOUGHKALMANTRACK_H
#define ALIL3HOUGHKALMANTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//             High Level Trigger TPC Hough Kalman Track Class
//
//        Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch 
//-------------------------------------------------------------------------


/*****************************************************************************
 *                          October 11, 2004                                 *
 * The class inherits from the off-line AliTPCtrack class.                   *
 * It is used to transform AliHLTHoughTrack into AliTPCTrack, which is        *
 * then stored as AliESDtrack object in the ESD                              *
 *****************************************************************************/

#include <AliTPCtrack.h>

class AliHLTHoughTrack;
class AliHLTHoughBaseTransformer;

class AliHLTHoughKalmanTrack : public AliTPCtrack {
public:
  AliHLTHoughKalmanTrack(const AliHLTHoughTrack& t) throw (const Char_t *);

  ClassDef(AliHLTHoughKalmanTrack,1)   //HLT TPC Hough track
};

typedef AliHLTHoughKalmanTrack AliL3HoughKalmanTrack; // for backward compatibility

#endif
