#ifndef ALIFLATESDTRACK_H
#define ALIFLATESDTRACK_H

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
#include "AliVMisc.h"
#include "AliVVtrack.h"
#include "AliFlatExternalTrackParam.h"

class AliESDtrack;
class AliExternalTrackParam;

class AliFlatESDTrack :public AliVVtrack {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors

  AliFlatESDTrack();
  ~AliFlatESDTrack() {}  

  // constructor and method for reinitialisation of virtual table
  AliFlatESDTrack( AliVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatESDTrack( AliVReinitialize ); }

  // --------------------   AliVVtrack interface    ---------------------------------

  Int_t GetTrackParam         ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x0  ); }
  Int_t GetTrackParamRefitted ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x1  ); }
  Int_t GetTrackParamIp       ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x2  ); }
  Int_t GetTrackParamTPCInner ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x4  ); }
  Int_t GetTrackParamOp       ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x8  ); }
  Int_t GetTrackParamCp       ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x10 ); }
  Int_t GetTrackParamITSOut   ( AliExternalTrackParam &p ) const { return GetExternalTrackParam( p, 0x20 ); }

  UShort_t GetTPCNcls() const {return GetNumberOfTPCClusters(); }
  Double_t GetPt() const {
    const AliFlatExternalTrackParam *f = GetFlatTrackParam();
    return (f) ?f->GetPt() : kVeryBig;
  }
  

  // --------------------------------------------------------------------------------

  // -- Set methods
 
  Int_t SetFromESDTrack( const AliESDtrack* track );

  Int_t SetExternalTrackParam( 
			      const AliExternalTrackParam* refittedParam,
			      const AliExternalTrackParam* innerParam,
			      const AliExternalTrackParam* innerTPC,
			      const AliExternalTrackParam* outerParam,
			      const AliExternalTrackParam* constrainedParam,
			      const AliExternalTrackParam* outerITSParam
			       );

  void SetNumberOfTPCClusters( Int_t nCl ) { fNTPCClusters = nCl; } 
  void SetNumberOfITSClusters( Int_t nCl ) { fNITSClusters = nCl; } 

  
  // --------------------------------------------------------------------------------
  // -- Getter methods

  const AliFlatExternalTrackParam* GetFlatTrackParam()         const { return GetFlatParam( 0x0  ); }
  const AliFlatExternalTrackParam* GetFlatTrackParamRefitted() const { return GetFlatParam( 0x1  ); }
  const AliFlatExternalTrackParam* GetFlatTrackParamIp()       const { return GetFlatParam( 0x2  ); } 
  const AliFlatExternalTrackParam* GetFlatTrackParamTPCInner() const { return GetFlatParam( 0x4  ); } 
  const AliFlatExternalTrackParam* GetFlatTrackParamOp()       const { return GetFlatParam( 0x8  ); }     
  const AliFlatExternalTrackParam* GetFlatTrackParamCp()       const { return GetFlatParam( 0x10 ); }
  const AliFlatExternalTrackParam* GetFlatTrackParamITSOut()   const { return GetFlatParam( 0x20 ); }

  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  

  Int_t GetNumberOfTPCClusters() const { return fNTPCClusters; } 
  Int_t GetNumberOfITSClusters() const { return fNITSClusters; } 
    
  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  
  const AliFlatESDTrack *GetNextTrack() const { return reinterpret_cast<const AliFlatESDTrack*>(fContent+fContentSize); }
  AliFlatESDTrack *GetNextTrackNonConst() { return reinterpret_cast<AliFlatESDTrack*>(fContent+fContentSize); }
 
  // --------------------------------------------------------------------------------
  // -- Size methods

  static size_t EstimateSize(){
    return sizeof(AliFlatESDTrack) + 6*sizeof(AliFlatExternalTrackParam);
  }

  size_t GetSize() { return fContent -  reinterpret_cast<Byte_t*>(this) + fContentSize; }
    
 private:

  AliFlatESDTrack(const AliFlatESDTrack&);
  AliFlatESDTrack& operator=(const AliFlatESDTrack&);

  const AliFlatExternalTrackParam* GetFlatParam( UShort_t flag ) const {
    if( flag==0 ) return ( fTrackParamMask ) ? reinterpret_cast<const AliFlatExternalTrackParam*>(fContent) : NULL;
    else return (fTrackParamMask & flag) ? reinterpret_cast<const AliFlatExternalTrackParam*>(fContent) + CountBits(fTrackParamMask, flag-1) : NULL;
  }

  Int_t GetExternalTrackParam( AliExternalTrackParam &p, UShort_t flag  ) const;

  Int_t FillExternalTrackParam(const AliExternalTrackParam* param, UShort_t flag);

  static UInt_t CountBits(Byte_t field, UInt_t mask = 0xFFFFFFFF);
 
  // --------------------------------------------------------------------------------
  // -- Fixed size member objects
  //    -> Try to align in memory

  Byte_t   fTrackParamMask;            // Bit mask specfifying which ExternalTrackParam are present
  Int_t    fNTPCClusters;                 // Number of TPC clusters in track
  Int_t    fNITSClusters;                 // Number of ITS clusters in track
  // Bool_t   fMCLabels;

  ULong64_t fContentSize;                      // Size of this object
  
  // --------------------------------------------------------------------------------
  // -- Variable Size Object
  Byte_t fContent[1];                  // Variale size object, which contains all data

};

inline UInt_t AliFlatESDTrack::CountBits(Byte_t field, UInt_t mask) {
  // Count bits in field
  UInt_t count = 0, reg = field & mask;
  for (; reg; count++) reg &= reg - 1; 
  return count;
}

inline Int_t AliFlatESDTrack::GetExternalTrackParam( AliExternalTrackParam &p, UShort_t flag) const
{
  // Get external track parameters  
  const AliFlatExternalTrackParam *f = GetFlatParam ( flag );
  if( !f ) return -1;  
  Float_t par[5] = { f->GetY(), f->GetZ(), f->GetSnp(), f->GetTgl(), f->GetSigned1Pt() };
  p.Set( f->GetX(), f->GetAlpha(), par, f->GetCov() );
  return 0;
}


#endif
