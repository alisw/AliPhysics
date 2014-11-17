#ifndef ALIFLATESDFRIENDTRACK_H
#define ALIFLATESDFRIENDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/*
 * See implementation file for documentation
 */


#include "Rtypes.h"

#include "AliVfriendTrack.h"
#include "AliVMisc.h"
#include "AliFlatTPCseed.h"
#include "AliFlatExternalTrackParam.h"

class AliESDtrack;
class AliESDfriendTrack;
class AliExternalTrackParam;
class AliTrackPointArray;
class AliTPCseed;

class AliFlatESDFriendTrack :public AliVfriendTrack 
{
 public:

  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatESDFriendTrack();
  ~AliFlatESDFriendTrack() {}
 
  // constructor and method for reinitialisation of virtual table  
  AliFlatESDFriendTrack( AliVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatESDFriendTrack( AliVReinitialize ); }

  // --------------------   AliVfriendTrack interface    ---------------------------------

  Int_t GetTPCseed( AliTPCseed &) const;

  Int_t GetTrackParamTPCOut( AliExternalTrackParam &p ) const {
    if( !fTPCOutFlag ) return -1;
    fTPCOut.GetExternalTrackParam(p);
    return 0;
  }
  Int_t GetTrackParamITSOut( AliExternalTrackParam &p ) const { return GetTrackParam( fITSOutPointer, p ); }
  Int_t GetTrackParamTRDIn( AliExternalTrackParam &p )  const { return GetTrackParam( fTRDInPointer,  p ); }

  //virtual const AliVtrackPointArray *GetTrackPointArray() const {return NULL;}

  TObject* GetCalibObject(Int_t) const {return NULL;}
  const AliExternalTrackParam* GetTPCOut() const {return NULL;}
  const AliExternalTrackParam* GetITSOut() const {return NULL;}

  
  // bit manipulation for filtering

  void SetSkipBit(Bool_t skip){ fBitFlags = skip; }
  Bool_t TestSkipBit() const { return (fBitFlags!=0); }
  
  // ------------------- Own methods  ---------------------------------------------------------

  // -- Set methods
 
  void Reset();

  Int_t SetFromESDfriendTrack( const AliESDfriendTrack* track, size_t allocatedMemory );

  void SetTrackParamTPCOut( const AliExternalTrackParam *p ){
    if( !p ){ fTPCOutFlag=0; return; }
    fTPCOut.SetExternalTrackParam( p );
    fTPCOutFlag = 1;
  }
  
  void AddTrackParamTPCOut( const AliExternalTrackParam *p ){ SetTrackParamTPCOut( p ); }

  void AddTrackParamITSOut( const AliExternalTrackParam *p ){ AddTrackParam( fITSOutPointer, p ); }
  void AddTrackParamTRDIn ( const AliExternalTrackParam *p ){ AddTrackParam( fTRDInPointer,  p );  }

  void AddTPCseed         ( const AliTPCseed *p );

  // -- 

  AliFlatTPCseed* AddTPCseedStart();
  void AddTPCseedEnd( size_t tpcSeedSize );

	
	
  const AliFlatTPCseed* GetFlatTPCseed( ) const{return reinterpret_cast<const AliFlatTPCseed*>(fContent+fTPCseedPointer); };
  
  const AliFlatESDFriendTrack *GetNextTrack() const { return reinterpret_cast<const AliFlatESDFriendTrack*>(fContent+fContentSize); }
  AliFlatESDFriendTrack *GetNextTrackNonConst() { return reinterpret_cast<AliFlatESDFriendTrack*>(fContent+fContentSize); }
 
  // --------------------------------------------------------------------------------
  // -- Size methods

  static size_t EstimateSize(){
    return sizeof(AliFlatESDFriendTrack) + 3*sizeof(AliFlatExternalTrackParam) + AliFlatTPCseed::EstimateSize();
  }

  size_t GetSize() const { return fContent -  reinterpret_cast<const Byte_t*>(this) + fContentSize; }

 private: 

  AliFlatESDFriendTrack(const AliFlatESDFriendTrack &);
  AliFlatESDFriendTrack& operator=(const AliFlatESDFriendTrack& ); 

  Int_t GetTrackParam( Long64_t ptr, AliExternalTrackParam &param ) const;
  void  AddTrackParam( Long64_t &ptr, const AliExternalTrackParam *p );

  // --------------------------------------------------------------------------------

  ULong64_t fContentSize;                      // Size of this object
  AliFlatExternalTrackParam fTPCOut; // TPC Out track parameters (always in the structure, because it may be reset in calibration)
  Long64_t  fITSOutPointer;        // pointer to ITSOut track param in fContent
  Long64_t  fTRDInPointer;        // pointer to TRDIn track param in fContent
  Long64_t  fTPCseedPointer;       // pointer to TPCseed in fContent
  Bool_t    fBitFlags; // bit flags
  Bool_t    fTPCOutFlag;        // flag if TPCOut is set

  // --------------------------------------------------------------------------------
  
  Byte_t fContent[1];                  // Variale size object, which contains all data

  ClassDef(AliFlatESDFriendTrack, 0)

};

inline Int_t AliFlatESDFriendTrack::GetTrackParam( Long64_t ptr, AliExternalTrackParam &param ) const
{
  if( ptr<0 ) return -1;
  const AliFlatExternalTrackParam *fp = reinterpret_cast< const AliFlatExternalTrackParam* >( fContent + ptr );
  fp->GetExternalTrackParam( param );
	return 0;
}

inline void AliFlatESDFriendTrack::AddTrackParam( Long64_t &ptr, const AliExternalTrackParam *p )
{
  if(!p ) return;
  if( ptr<0 ){
    ptr = fContentSize;
    fContentSize += sizeof(AliFlatExternalTrackParam);
  }
  AliFlatExternalTrackParam *fp = reinterpret_cast< AliFlatExternalTrackParam* >( fContent + ptr );
  fp->SetExternalTrackParam( p );
}

inline void AliFlatESDFriendTrack::AddTPCseed( const AliTPCseed *p )
{
  fTPCseedPointer = -1;
  if(!p ) return;
  fTPCseedPointer = fContentSize;
  AliFlatTPCseed *fp = reinterpret_cast< AliFlatTPCseed* >( fContent + fTPCseedPointer );
  fp->SetFromTPCseed( p );
  fContentSize += fp->GetSize();  
}

inline Int_t AliFlatESDFriendTrack::GetTPCseed( AliTPCseed &s ) const
{
  if( fTPCseedPointer<0 ) return -1;
  const AliFlatTPCseed *fp = reinterpret_cast< const AliFlatTPCseed* >( fContent + fTPCseedPointer );
  fp->GetTPCseed( &s );
  return 0;
}

inline AliFlatTPCseed* AliFlatESDFriendTrack::AddTPCseedStart()
{
  fTPCseedPointer = fContentSize;
  return  reinterpret_cast< AliFlatTPCseed* >( fContent + fTPCseedPointer );
}

inline void AliFlatESDFriendTrack::AddTPCseedEnd( size_t tpcSeedSize ){
  fContentSize += tpcSeedSize;
}

#endif
