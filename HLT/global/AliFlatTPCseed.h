#ifndef ALIFLATTPCSEED_H
#define ALIFLATTPCSEED_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */


#include "Rtypes.h"

#include "AliFlatTPCCluster.h"
#include "AliVfriendTrack.h"
#include "AliVMisc.h"
#include "AliFlatExternalTrackParam.h"
#include "AliFlatTPCseed.h"

class AliESDtrack;
class AliESDfriendTrack;
class AliExternalTrackParam;
class AliTrackPointArray;
class AliTPCseed;

class AliFlatTPCseed  
{
 public:

  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatTPCseed();
  ~AliFlatTPCseed() {}
 
  // constructor and method for reinitialisation of virtual table  
  AliFlatTPCseed( AliVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatTPCseed( AliVReinitialize ); }

  // -- Set methods
 
  void Reset();

  void SetFromTPCseed( const AliTPCseed *p );
  void GetTPCseed( AliTPCseed *p ) const;

  
  // --------------------------------------------------------------------------------
  // -- Size methods

  static size_t EstimateSize(){
    return sizeof(AliFlatTPCseed) + 6*sizeof(AliFlatExternalTrackParam);
  }

  size_t GetSize() const { return fContent -  reinterpret_cast<const Byte_t*>(this) + fContentSize; }

 private: 

  AliFlatTPCseed(const AliFlatTPCseed &);
  AliFlatTPCseed& operator=(const AliFlatTPCseed& ); 

  // --------------------------------------------------------------------------------

  ULong64_t fContentSize;                      // Size of this object

  // --------------------------------------------------------------------------------
  
  Byte_t fContent[1];                  // Variale size object, which contains all data

};


#endif
