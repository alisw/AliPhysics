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
 * It is used to transform AliL3HoughTrack into AliTPCTrack, which is        *
 * then stored as AliESDtrack object in the ESD                              *
 *****************************************************************************/

#include <AliTPCtrack.h>

class AliL3HoughTrack;
class AliL3HoughBaseTransformer;

class AliL3HoughKalmanTrack : public AliTPCtrack {
public:
  AliL3HoughKalmanTrack(const AliL3HoughTrack& t) throw (const Char_t *);

  ClassDef(AliL3HoughKalmanTrack,1)   //HLT TPC Hough track
};

#endif
