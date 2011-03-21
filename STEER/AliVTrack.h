#ifndef AliVTrack_H
#define AliVTrack_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//     base class for ESD and AOD tracks
//     Author: A. Dainese
//-------------------------------------------------------------------------

#include <TBits.h>

#include "AliVParticle.h"

class AliVVertex;
class AliExternalTrackParam;

class AliVTrack: public AliVParticle {

public:
  AliVTrack() { }
  virtual ~AliVTrack() { }
  AliVTrack(const AliVTrack& vTrack); 
  AliVTrack& operator=(const AliVTrack& vTrack);

  virtual Int_t    GetID() const = 0;
  virtual UChar_t  GetITSClusterMap() const = 0;
  virtual Float_t  GetTPCClusterInfo(Int_t /*nNeighbours*/, Int_t /*type*/, Int_t /*row0*/=0, Int_t /*row1*/=159) const {return 0.;}
  virtual UShort_t GetTPCNcls() const { return 0;}
  virtual UShort_t GetTPCNclsF() const { return 0;}
  
  //pid info
  virtual Double_t  GetITSsignal()       const {return 0.;}
  virtual Double_t  GetTPCsignal()       const {return 0.;}
  virtual UShort_t  GetTPCsignalN()      const {return 0 ;}
  virtual Double_t  GetTPCmomentum()     const {return 0.;}
  virtual Double_t  GetTOFsignal()       const {return 0.;}
  
  virtual ULong_t  GetStatus() const = 0;
  virtual Bool_t   GetXYZ(Double_t *p) const = 0;
  virtual Double_t GetBz() const;
  virtual void     GetBxByBz(Double_t b[3]) const;
  virtual Bool_t   GetCovarianceXYZPxPyPz(Double_t cv[21]) const = 0;
  virtual Bool_t   PropagateToDCA(const AliVVertex *vtx,Double_t b,Double_t maxd,Double_t dz[2],Double_t covar[3]) = 0;
  virtual const    AliExternalTrackParam * GetOuterParam() const { return NULL; }
  virtual Int_t    GetNcls(Int_t /*idet*/) const { return 0; }
  virtual Bool_t   GetPxPyPz(Double_t */*p*/) const { return kFALSE; }
  virtual void     SetID(Short_t /*id*/) {;}


  ClassDef(AliVTrack,1)  // base class for tracks
};

#endif
