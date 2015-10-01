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
 
  Int_t GetTrackParamTPCOut( AliExternalTrackParam &p ) const { return GetTrackParam( fTPCOutPointer, p ); }
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

  void ResetTrackParamTPCOut( const AliExternalTrackParam *p ){ ResetTrackParam( fTPCOutPointer, p ); }
  void ResetTPCseed( const AliTPCseed* s );
  
  void SetTrackParamTPCOut( const AliExternalTrackParam *p ){ SetTrackParam( fTPCOutPointer, p ); }
  void SetTrackParamITSOut( const AliExternalTrackParam *p ){ SetTrackParam( fITSOutPointer, p ); }
  void SetTrackParamTRDIn ( const AliExternalTrackParam *p ){ SetTrackParam( fTRDInPointer,  p );  }

  void SetTPCseed         ( const AliTPCseed *p );

  // -- 

  AliFlatTPCseed* SetTPCseedStart();
  void SetTPCseedEnd( size_t tpcSeedSize );

	
	
  const AliFlatTPCseed* GetFlatTPCseed( ) const{return fTPCseedPointer<0? NULL : reinterpret_cast<const AliFlatTPCseed*>(fContent+fTPCseedPointer); };
  
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
  void  SetTrackParam( Long64_t &ptr, const AliExternalTrackParam *p );
  void  ResetTrackParam( Long64_t &ptr, const AliExternalTrackParam *p );

  // --------------------------------------------------------------------------------

  ULong64_t fContentSize;                      // Size of this object
  Long64_t  fTPCOutPointer;        // pointer to TPCOut track param in fContent
  Long64_t  fITSOutPointer;        // pointer to ITSOut track param in fContent
  Long64_t  fTRDInPointer;        // pointer to TRDIn track param in fContent
  Long64_t  fTPCseedPointer;       // pointer to TPCseed in fContent
  Bool_t    fBitFlags; // bit flags

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

inline void AliFlatESDFriendTrack::SetTrackParam( Long64_t &ptr, const AliExternalTrackParam *p )
{
  if(!p ) return;
  if( ptr<0 ){
    ptr = fContentSize;
    fContentSize += sizeof(AliFlatExternalTrackParam);
  }
  AliFlatExternalTrackParam *fp = reinterpret_cast< AliFlatExternalTrackParam* >( fContent + ptr );
  fp->SetExternalTrackParam( p );
}

inline void AliFlatESDFriendTrack::ResetTrackParam( Long64_t &ptr, const AliExternalTrackParam *p )
{
  if( ptr<0 || !p ) return;
  AliFlatExternalTrackParam *fp = reinterpret_cast< AliFlatExternalTrackParam* >( fContent + ptr );
  fp->SetExternalTrackParam( p );
}

inline void AliFlatESDFriendTrack::ResetTPCseed( const AliTPCseed* s )
{
  if (!s) return;
  AliFlatTPCseed* fs = reinterpret_cast<AliFlatTPCseed*>(fContent+fTPCseedPointer);
  fs->SetFromTPCseed( s );
}

inline void AliFlatESDFriendTrack::SetTPCseed( const AliTPCseed *p )
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

inline AliFlatTPCseed* AliFlatESDFriendTrack::SetTPCseedStart()
{
  fTPCseedPointer = fContentSize;
  return  reinterpret_cast< AliFlatTPCseed* >( fContent + fTPCseedPointer );
}

inline void AliFlatESDFriendTrack::SetTPCseedEnd( size_t tpcSeedSize ){
  fContentSize += tpcSeedSize;
}

#endif
