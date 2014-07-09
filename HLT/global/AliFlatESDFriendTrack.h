#ifndef ALIFLATESDFRIENDTRACK_H
#define ALIFLATESDFRIENDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */

/*
Cp - Track parameters constrained to the primary vertex
Ip - Track parameters estimated at the inner wall of TPC
TPCInner - Track parameters estimated at the inner wall of TPC using the TPC stand-alone 
Op - Track parameters estimated at the point of maximal radial coordinate reached during the tracking
*/

#include "Rtypes.h"

#include "AliFlatTPCCluster.h"
#include "AliVVfriendTrack.h"
#include "AliFlatESDMisc.h"

class AliESDtrack;
class AliESDfriendTrack;
class AliExternalTrackParam;

class AliFlatESDFriendTrack :public AliVVfriendTrack 
{
 public:
  AliFlatESDFriendTrack();
  ~AliFlatESDFriendTrack() {}

  //implementation of AliVVfriendTrack methods 

  AliVVTPCseed* GetTPCseed() const {return NULL;}
  AliVVTRDseed* GetTRDseed() const {return NULL;}
  const AliTrackPointArray *GetTrackPointArray() const { return NULL; }
  const AliExternalTrackParam * GetITSOut() const { return NULL; } 
  const AliExternalTrackParam * GetTPCOut() const { return  NULL; } 
  const AliExternalTrackParam * GetTRDIn()  const { return NULL; } 

  // own methods

  void Reinitialize()
  {
    new (this) AliFlatESDFriendTrack(AliFlatESDReinitialize);
  }

private: 
  AliFlatESDFriendTrack(const AliFlatESDFriendTrack &);
  AliFlatESDFriendTrack& operator=(const AliFlatESDFriendTrack& );   

  // special constructor, to be called by placement new,
  // when accessing information after reinterpret_cast
  // so that vtable is generated, but values are not overwritten
  AliFlatESDFriendTrack(AliFlatESDSpecialConstructorFlag);

};


#endif
