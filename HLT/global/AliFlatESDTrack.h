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
#include "AliVVMisc.h"
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
  AliFlatESDTrack( AliVVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatESDTrack( AliVVReinitialize ); }

 // --------------------------------------------------------------------------------
  // -- Set methods
 
  Int_t Set( const AliESDtrack* track );

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
  AliFlatExternalTrackParam* GetTrackParamRefitted(){
    return (fTrackParamMask & 0x1) ? reinterpret_cast<AliFlatExternalTrackParam*>(fContent) : NULL;
  } 

  AliFlatExternalTrackParam* GetTrackParamIp() { 
    return (fTrackParamMask & 0x2) ? reinterpret_cast<AliFlatExternalTrackParam*>(fContent) + CountBits(fTrackParamMask, 0x1) : NULL;
  } 

  AliFlatExternalTrackParam* GetTrackParamTPCInner() { 
    return (fTrackParamMask & 0x4) ? reinterpret_cast<AliFlatExternalTrackParam*>(fContent) + CountBits(fTrackParamMask, 0x3) : NULL;
  } 

  AliFlatExternalTrackParam* GetTrackParamOp() {      
    return (fTrackParamMask & 0x8) ? reinterpret_cast<AliFlatExternalTrackParam*>(fContent) + CountBits(fTrackParamMask, 0x7) : NULL;
  } 

  AliFlatExternalTrackParam* GetTrackParamCp() {
    return (fTrackParamMask & 0x10) ? reinterpret_cast<AliFlatExternalTrackParam*>(fContent) + CountBits(fTrackParamMask, 0xF) : NULL;
  } 

  AliFlatExternalTrackParam* GetTrackParamITSOut() {
    return (fTrackParamMask & 0x20) ? reinterpret_cast<AliFlatExternalTrackParam*>(fContent) + CountBits(fTrackParamMask, 0x1F) : NULL;
  } 

  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  

  Int_t GetNumberOfTPCClusters() {
    return fNTPCClusters;
  } 
  
  Int_t GetNumberOfITSClusters() {
    return fNITSClusters;
  } 
    
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

#endif
