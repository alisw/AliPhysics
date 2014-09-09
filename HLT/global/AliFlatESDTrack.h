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
#include "AliVParticle.h"
#include "AliFlatExternalTrackParam.h"

class AliESDtrack;
class AliExternalTrackParam;
class AliFlatESDTrack;

class AliFlatESDTrack :public AliVTrack {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors

  AliFlatESDTrack();
  virtual ~AliFlatESDTrack() {}  

  // constructor and method for reinitialisation of virtual table
  AliFlatESDTrack( AliVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatESDTrack( AliVReinitialize ); }

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

  size_t GetSize() const { return fContent -  reinterpret_cast<const Byte_t*>(this) + fContentSize; }
    
  // ---------------------------------------------------------------------------------
  // AliVParticle interface
  virtual Double_t Pt() const {const AliFlatExternalTrackParam* p=GetFlatTrackParam(); return (p)?p->GetPt():0.;}
  virtual Double_t Px() const {return 0.;}
  virtual Double_t Py() const {return 0.;}
  virtual Double_t Pz() const {return 0.;}
  virtual Double_t P() const {return 0.;}
  virtual Bool_t PxPyPz(Double_t*) const {return kFALSE;}
  virtual Double_t Xv() const {return 0.;}
  virtual Double_t Yv() const {return 0.;}
  virtual Double_t Zv() const {return 0.;}
  virtual Bool_t XvYvZv(Double_t*) const {return 0.;}
  virtual Double_t OneOverPt() const {return 0.;}
  virtual Double_t Phi() const {return 0.;}
  virtual Double_t Theta() const {return 0.;}
  virtual Double_t E() const {return 0.;}
  virtual Double_t M() const {return 0.;}
  virtual Double_t Eta() const {return 0.;}
  virtual Double_t Y() const {return 0.;}
  virtual Short_t Charge() const {return 0.;}
  virtual Int_t GetLabel() const {return 0.;}
  virtual Int_t PdgCode() const {return 0.;}
  virtual const Double_t* PID() const {return NULL;} 
  virtual Int_t    GetID() const {return 0.;}
  virtual UChar_t  GetITSClusterMap() const {return 0.;}
  virtual ULong_t  GetStatus() const {return 0.;}
  virtual Bool_t   GetCovarianceXYZPxPyPz(Double_t cv[21]) const {if (cv[0]); return kFALSE;}
  virtual Bool_t   PropagateToDCA(const AliVVertex* /*vtx*/, Double_t /*b*/, Double_t /*maxd*/, Double_t dz[2], Double_t covar[3]) {if (dz[0]==covar[3]); return kFALSE;}

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

// _______________________________________________________________________________________________________
inline AliFlatESDTrack::AliFlatESDTrack() :
  fTrackParamMask(0),
  fNTPCClusters(0),
  fNITSClusters(0),
  fContentSize(0)
{
  // Default constructor
}

inline AliFlatESDTrack::AliFlatESDTrack( AliVConstructorReinitialisationFlag f )
  :
  fTrackParamMask(fTrackParamMask ),
  fNTPCClusters( fNTPCClusters ),
  fNITSClusters( fNITSClusters ),
  fContentSize( fContentSize )
{
  // Constructor for reinitialisation of vtable
}

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
  f->GetExternalTrackParam( p );
  return 0;
}


#endif
