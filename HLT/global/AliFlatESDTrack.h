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

#include "AliFlatTPCCluster.h"
#include "AliFlatExternalTrackParam.h"
#include "AliVVtrack.h"

class AliESDtrack;
class AliESDfriendTrack;
class AliExternalTrackParam;

//class AliFlatESDTrack: public AliVVtrack {
class AliFlatESDTrack {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliFlatESDTrack();   
  AliFlatESDTrack(const AliESDtrack* track, AliESDfriendTrack* friendTrack); 
  ~AliFlatESDTrack();  

  // --------------------------------------------------------------------------------
  // -- Fill / Set methods
  Int_t FillExternalTrackParam( 
			       const AliExternalTrackParam* refittedParam,
			       const AliExternalTrackParam* innerParam,
			       const AliExternalTrackParam* innerTPC,
			       const AliExternalTrackParam* outerParam,
			       const AliExternalTrackParam* constrainedParam,
			       const AliExternalTrackParam* outerITSParam
			     );

  AliFlatTPCCluster *GetNextTPCClusterPointer(){ return GetTPCCluster(fNTPCClusters); }

  void StoreLastTPCCluster(){  
     ++fNTPCClusters;
     fSize += sizeof(AliFlatTPCCluster);
  }

  void SetNumberOfITSClusters( Int_t nCl ) { fNITSClusters = nCl; } 

  Int_t Fill( const AliESDtrack* track, AliESDfriendTrack* friendTrack);
  
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
  
  AliFlatTPCCluster *GetTPCClusters() {
    return reinterpret_cast< AliFlatTPCCluster*>(fContent + sizeof(AliFlatExternalTrackParam)*CountBits(fTrackParamMask));
  } 
  
  AliFlatTPCCluster *GetTPCCluster(Int_t ind) {
    return GetTPCClusters() + ind*sizeof(AliFlatTPCCluster);
  }

  Int_t GetNumberOfITSClusters() {
    return fNITSClusters;
  } 
  
  
  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  
  AliFlatESDTrack *GetNextTrack() {return reinterpret_cast<AliFlatESDTrack*>(fContent+fSize);}
  
  // --------------------------------------------------------------------------------
  // -- Size methods
  static ULong64_t EstimateSize(Bool_t useESDFriends = kTRUE, Int_t nTPCClusters = 160 );
         ULong64_t GetSize()  {return fContent -  reinterpret_cast<Byte_t*>(this) + fSize;}
    
 private:
  AliFlatESDTrack(const AliFlatESDTrack&);
  AliFlatESDTrack& operator=(const AliFlatESDTrack&);

  Int_t FillExternalTrackParam(const AliExternalTrackParam* param, UShort_t flag);

  UInt_t CountBits(Byte_t field, UInt_t mask = 0xFFFFFFFF);
 
  // --------------------------------------------------------------------------------
  // -- Fixed size member objects
  //    -> Try to align in memory

  Byte_t   fTrackParamMask;            // Bit mask specfifying which ExternalTrackParam are present
  Int_t    fNTPCClusters;                 // Number of TPC clusters in track
  Int_t    fNITSClusters;                 // Number of ITS clusters in track
  // Bool_t   fMCLabels;

  ULong64_t fSize;                      // Size of this object
  
  // --------------------------------------------------------------------------------
  // -- Variable Size Object
  Byte_t fContent[1];                  // Variale size object, which contains all data

};
#endif
