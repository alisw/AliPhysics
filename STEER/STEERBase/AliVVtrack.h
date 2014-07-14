#ifndef ALIVVTRACK_H
#define ALIVVTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Mikolaj Krzewicki mkrzewic@cern.ch     */

/*
 * See implementation file for documentation
 */

#include "AliPID.h"

struct AliFlatTPCCluster;
struct AliFlatExternalTrackParam;

class AliVVtrack {
 public:
  // --------------------------------------------------------------------------------
  // -- Constructor / Destructors
  AliVVtrack() {} 
  virtual ~AliVVtrack() {}

  // --------------------------------------------------------------------------------
  // -- Getter methods
  virtual AliFlatExternalTrackParam* GetTrackParamRefitted() ;
  virtual AliFlatExternalTrackParam* GetTrackParamIp() ;
  virtual AliFlatExternalTrackParam* GetTrackParamTPCInner() ;
  virtual AliFlatExternalTrackParam* GetTrackParamOp() ;
  virtual AliFlatExternalTrackParam* GetTrackParamCp() ;
  virtual AliFlatExternalTrackParam* GetTrackParamITSOut() ;

  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  virtual Int_t GetNumberOfTPCClusters() ;
  virtual AliFlatTPCCluster* GetTPCClusters() ;
  virtual AliFlatTPCCluster* GetTPCCluster(Int_t /*ind*/) ;
  virtual Int_t GetNumberOfITSClusters() ;
  virtual AliVVtrack *GetNextTrack() ;

  virtual Bool_t GetXYZ(Double_t* p) const ;
  virtual Bool_t GetXYZAt(Double_t x, Double_t y, Double_t* z) const ;

  //AliESDv0
  virtual void  GetXYZ(Double_t& x, Double_t& y, Double_t& z) const ;

  virtual Double_t GetTgl()  const ;
  virtual UShort_t GetTPCNclsF() const ;

  virtual Double_t GetTOFsignalDz() const ;
  virtual void GetImpactParameters(Float_t& /*xy*/,Float_t& /*z*/) const ;
  //TODO:
  virtual void GetDZ(Double_t /*x*/,Double_t /*y*/,Double_t /*z*/,Double_t /*b*/, Float_t dz[2]) const ;

  virtual Float_t GetTPCClusterInfo(Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159, Int_t bitType=0 ) const ;
  virtual UShort_t GetTPCncls(Int_t row0=0,Int_t row1=159) const ;
  virtual Bool_t IsOn(Int_t /*mask*/) const ;
  virtual void GetDirection(Double_t d[3]) const ;
  virtual const Double_t *GetParameter() const ;
  virtual void GetImpactParametersTPC(Float_t& /*xy*/,Float_t& /*z*/) const ;
  virtual Int_t GetNumberOfClusters() const ;
  virtual const AliVVtrack* GetTPCInnerParam() const ;
  virtual Double_t Pt() const ;
  virtual Double_t GetP() const ;
  virtual Double_t GetTPCmomentum() const ;
  virtual ULong_t GetStatus() const ;
  virtual const AliVVtrack * GetOuterParam() const ;
  virtual const AliVVtrack * GetInnerParam() const ;
  virtual Int_t GetKinkIndex(Int_t /*i*/) const ;
  virtual Double_t Eta() const ;
  virtual Double_t GetY() const ;
  
  virtual Double_t GetX() const ;
  virtual Double_t GetZ() const ;
  virtual Int_t GetNcls(Int_t /*idet*/) const ;
  virtual void GetIntegratedTimes(Double_t* /*times*/, Int_t nspec=AliPID::kSPECIES) const ;
  virtual Double_t GetSigned1Pt()  const ;
  virtual Double_t GetLinearD(Double_t /*xv*/, Double_t /*yv*/) const ;
  virtual const AliVVtrack *GetConstrainedParam() const ;
  virtual Double_t GetAlpha() const ;
  virtual Char_t GetITSclusters(Int_t* /*idx*/) const ;
  virtual Double_t GetSign() const ;
  virtual UShort_t GetTPCNcls() const ;
  virtual Float_t GetTPCCrossedRows() const ;
  virtual Double_t GetTPCsignal() const ;
  virtual Double_t GetTOFsignal() const ;
  virtual UChar_t GetTRDclusters(Int_t* /*idx*/) const ;
  
  //AliTPCtrack
  virtual Int_t GetNFoundable() const ;
  virtual Double_t GetdEdx()  const ;

  ClassDef(AliVVtrack, 1)   // base class for track data

};
#endif
