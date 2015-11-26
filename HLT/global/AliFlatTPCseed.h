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
#include "AliTPCseed.h"

class AliESDtrack;
class AliESDfriendTrack;
class AliExternalTrackParam;

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

  // -- Getters 

  void GetTPCseed( AliTPCseed *p ) const;

  Int_t GetLabel() const { return fLabel; }
  Int_t GetNClusters() const { return fNTPCClusters; }

  const AliFlatTPCCluster *GetClusters() const { return reinterpret_cast< const AliFlatTPCCluster* >( fContent ); }
  AliFlatTPCCluster *GetClustersNonConst(){ return reinterpret_cast< AliFlatTPCCluster* >( fContent ); }

  // -- Set methods
 
  void Reset();

  void SetFromTPCseed( const AliTPCseed *p );

  void SetExternalTrackParam( const AliExternalTrackParam* p ){ fParam.SetExternalTrackParam(  p ); }

  void SetLabel( Int_t lab ){ fLabel=lab; }

  void AddCluster( const AliTPCclusterMI *cl, const AliTPCTrackerPoints::Point *p ){
		if(cl){
			GetClustersNonConst()[fNTPCClusters++].SetTPCCluster( cl, p );
			fContentSize+=sizeof(AliFlatTPCCluster);
		}
	}


  // --------------------------------------------------------------------------------
  // -- Size methods

  static size_t EstimateSize(){
    return sizeof(AliFlatTPCseed) + 160*sizeof(AliFlatTPCCluster);
  }

  size_t GetSize() const { return fContent -  reinterpret_cast<const Byte_t*>(this) + fContentSize; }

 private: 

  AliFlatTPCseed(const AliFlatTPCseed &);
  AliFlatTPCseed& operator=(const AliFlatTPCseed& ); 

  // --------------------------------------------------------------------------------

  ULong64_t fContentSize;                      // Size of this object
  
  // --------------------------------------------------------------------------------
  
  AliFlatExternalTrackParam fParam;
  Int_t fLabel;
  Int_t fNTPCClusters;
  Byte_t fContent[1];                  // Variale size object, which contains all data

};

inline AliFlatTPCseed::AliFlatTPCseed()
  :
  fContentSize(0),
  fParam(),
  fLabel(-1),
  fNTPCClusters(0)
{
  // constructor
  fContent[0]=0;
}

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatTPCseed::AliFlatTPCseed( AliVConstructorReinitialisationFlag  )
{
  // reinitialise vtable
  fParam.Reinitialize();
  AliFlatTPCCluster *clusters = reinterpret_cast< AliFlatTPCCluster* >( fContent );  
  for( Int_t ic=0; ic<fNTPCClusters; ic++ ){
    clusters[ic].Reinitialize();
  }
}
#pragma GCC diagnostic warning "-Weffc++" 


inline void AliFlatTPCseed::Reset()
{
  // Reset
  fContentSize = 0;
  fNTPCClusters = 0;
  fLabel=-1;
}

#endif
