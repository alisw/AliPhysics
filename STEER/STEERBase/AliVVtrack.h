#ifndef ALIVVTRACK_H
#define ALIVVTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Mikolaj Krzewicki mkrzewic@cern.ch     */

/*
 * See implementation file for documentation
 */

#include "Rtypes.h"
#include "AliVVMisc.h"
class AliExternalTrackParam;

class AliVVtrack {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliVVtrack() {} 
  virtual ~AliVVtrack() {}

  // constructor and method for reinitialisation of virtual table
  AliVVtrack( AliVVConstructorReinitialisationFlag ) {}
  void Reinitialize() { new (this) AliVVtrack( AliVVReinitialize ); }

 // --------------------------------------------------------------------------------

  // --------------------------------------------------------------------------------
  // -- Getter methods
  /*
  virtual Int_t GetTrackParamRefitted( AliExternalTrackParam & ) const = 0 ;
  virtual Int_t GetTrackParamIp( AliExternalTrackParam & ) const = 0 ;
  virtual Int_t GetTrackParamTPCInner( AliExternalTrackParam & ) const = 0 ;
  virtual Int_t GetTrackParamOp( AliExternalTrackParam & ) const = 0 ;
  virtual Int_t GetTrackParamCp( AliExternalTrackParam & ) const = 0 ;
  virtual Int_t GetTrackParamITSOut( AliExternalTrackParam & ) const = 0 ;
  */
  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  virtual UShort_t GetTPCNcls() const {return 0;}//= 0;
  virtual Double_t GetPt() const {return 0;}//= 0;


  // may be for the future

  // virtual Float_t GetTPCClusterInfo(Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159, Int_t bitType=0 ) const ;
  // virtual UShort_t GetTPCncls(Int_t row0=0,Int_t row1=159) const ;
  // virtual Bool_t IsOn(Int_t /*mask*/) const ;
  // virtual void GetImpactParametersTPC(Float_t& /*xy*/,Float_t& /*z*/) const ;
  // virtual ULong_t GetStatus() const ;
  // virtual Int_t GetKinkIndex(Int_t /*i*/) const ;
  // virtual Int_t GetNcls(Int_t /*idet*/) const ;
  // virtual void GetIntegratedTimes(Double_t* /*times*/, Int_t nspec=AliPID::kSPECIES) const ;
  // virtual Char_t GetITSclusters(Int_t* /*idx*/) const ;
  // virtual Float_t GetTPCCrossedRows() const ;
  // virtual Double_t GetTPCsignal() const ;
  // virtual Double_t GetTOFsignal() const ;
  // virtual UChar_t GetTRDclusters(Int_t* /*idx*/) const ;  

  //ClassDef(AliVVtrack, 0)   // base class for track data

};
#endif
