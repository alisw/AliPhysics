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
#include "AliVVMisc.h"

class AliESDtrack;
class AliESDfriendTrack;
class AliExternalTrackParam;
class AliTrackPointArray;
class AliVVtrack;

class AliFlatESDFriendTrack :public AliVVfriendTrack 
{
 public:
  AliFlatESDFriendTrack();
  ~AliFlatESDFriendTrack() {}
  // constructor for reinitialisation of vtable
  AliFlatESDFriendTrack( AliVVConstructorReinitialisationFlag );

  //implementation of AliVVfriendTrack methods 

  //AliVVTPCseed* GetTPCseed() const {return NULL;}
  AliVVtrack* GetTPCseed() const {return NULL;}
  //AliVVTRDseed* GetTRDseed() const {return NULL;}
  const AliVVtrackPointArray *GetTrackPointArray() const { return NULL; }
  //const AliExternalTrackParam * GetITSOut() const { return NULL; } 
  //const AliExternalTrackParam * GetTPCOut() const { return  NULL; } 
  //const AliExternalTrackParam * GetTRDIn()  const { return NULL; } 
  const AliVVtrack * GetITSOut() const { return NULL; } 
  const AliVVtrack * GetTPCOut() const { return  NULL; } 
  const AliVVtrack * GetTRDIn()  const { return NULL; } 

  // own methods

  void Reinitialize() { new (this) AliFlatESDFriendTrack( AliVVReinitialize ); }
  
 private: 

  AliFlatESDFriendTrack(const AliFlatESDFriendTrack &);
  AliFlatESDFriendTrack& operator=(const AliFlatESDFriendTrack& ); 

};


#endif
