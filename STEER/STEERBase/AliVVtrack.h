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
  virtual AliFlatExternalTrackParam* GetTrackParamRefitted() { return NULL; } 
  virtual AliFlatExternalTrackParam* GetTrackParamIp() { return NULL; } 
  virtual AliFlatExternalTrackParam* GetTrackParamTPCInner() { return NULL; } 
  virtual AliFlatExternalTrackParam* GetTrackParamOp() { return NULL; } 
  virtual AliFlatExternalTrackParam* GetTrackParamCp() { return NULL; } 
  virtual AliFlatExternalTrackParam* GetTrackParamITSOut() { return NULL; } 

  // --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  
  virtual Int_t GetNumberOfTPCClusters() { return 0; } 
  virtual AliFlatTPCCluster* GetTPCClusters() { return NULL; } 
  virtual AliFlatTPCCluster* GetTPCCluster(Int_t /*ind*/) { return NULL; } 
  virtual Int_t GetNumberOfITSClusters() { return 0; }
  virtual AliVVtrack *GetNextTrack() {return NULL; }

  virtual Bool_t GetXYZ(Double_t* /*p*/) const {return kFALSE;}
  virtual void  GetXYZ(Double_t& /*x*/, Double_t& /*y*/, Double_t& /*z*/) const {}
  virtual Double_t GetTgl()  const {return 0.;}
  virtual UShort_t GetTPCNclsF() const { return 0;}

  virtual Double_t GetTOFsignalDz() const {return 0.;}
  virtual void GetImpactParameters(Float_t& /*xy*/,Float_t& /*z*/) const {}
  //TODO:
  virtual void GetDZ(Double_t /*x*/,Double_t /*y*/,Double_t /*z*/,Double_t /*b*/, Float_t dz[2]) const {if (dz[0]==0) return;}

  virtual Float_t GetTPCClusterInfo(Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159, Int_t bitType=0 ) const {return 0.*nNeighbours*type*row0*row1*bitType;}
  virtual UShort_t GetTPCncls(Int_t row0=0,Int_t row1=159) const {return 0*row0*row1;}
  virtual Bool_t IsOn(Int_t /*mask*/) const {return kFALSE;}
  virtual void GetDirection(Double_t d[3]) const {if (d[0]==0) return;}
  virtual const Double_t *GetParameter() const {return 0;}
  virtual void GetImpactParametersTPC(Float_t& /*xy*/,Float_t& /*z*/) const {}
  virtual Int_t GetNumberOfClusters() const {return 0;} 
  virtual const AliVVtrack* GetTPCInnerParam() const {return NULL;}
  virtual Double_t Pt() const {return 0.;}
  virtual Double_t GetP() const {return 0.;}
  virtual Double_t GetTPCmomentum() const {return 0.;}
  virtual ULong_t GetStatus() const {return 0;}
  virtual const AliVVtrack * GetOuterParam() const { return NULL;}
  virtual const AliVVtrack * GetInnerParam() const { return NULL;}
  virtual Int_t GetKinkIndex(Int_t /*i*/) const { return 0;}
  virtual Double_t Eta() const {return 0.;}
  virtual Double_t GetY() const {return 0.;}
  
  virtual Double_t GetX() const {return 0.;}
  virtual Double_t GetZ() const {return 0.;}
  virtual Int_t GetNcls(Int_t /*idet*/) const {return 0;}
  virtual void GetIntegratedTimes(Double_t* /*times*/, Int_t nspec=AliPID::kSPECIES) const {if (nspec<0) return;}
  virtual Double_t GetSigned1Pt()  const {return 0.;}  
  virtual Double_t GetLinearD(Double_t /*xv*/, Double_t /*yv*/) const {return 0.;}
  virtual const AliVVtrack *GetConstrainedParam() const {return NULL;}
  virtual Double_t GetAlpha() const {return 0.;}
  virtual Char_t GetITSclusters(Int_t* /*idx*/) const {return 0;}
  virtual Double_t GetSign() const {return 0.;}
  virtual UShort_t GetTPCNcls() const { return 0;}
  virtual Float_t GetTPCCrossedRows() const {return 0.;}
  virtual Double_t GetTPCsignal() const {return 0.;}
  virtual Double_t GetTOFsignal() const {return 0.;}
  virtual UChar_t GetTRDclusters(Int_t* /*idx*/) const {return 0;}
  
  //AliTPCtrack
  virtual Int_t GetNFoundable() const {return 0;} 
  virtual Double_t GetdEdx()  const {return 0.;}

  ClassDef(AliVVtrack, 1)   // base class for track data

};
#endif
