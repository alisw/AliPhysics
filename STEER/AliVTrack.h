#ifndef AliVTrack_H
#define AliVTrack_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//     base class for ESD and AOD tracks
//     Author: A. Dainese
//-------------------------------------------------------------------------

#include "AliVParticle.h"


class AliVTrack: public AliVParticle {

public:
  AliVTrack() { }
  virtual ~AliVTrack() { }
  AliVTrack(const AliVTrack& vTrack); 
  AliVTrack& operator=(const AliVTrack& vTrack);

  virtual Int_t    GetID() const = 0;
  virtual UChar_t  GetITSClusterMap() const = 0;
  virtual UShort_t GetTPCNcls() const { return 0;}
  virtual ULong_t  GetStatus() const = 0;
  virtual Bool_t   GetXYZ(Double_t *p) const = 0;
  virtual Double_t GetBz() const;
  virtual void     GetBxByBz(Double_t b[3]) const;
  virtual Bool_t   GetCovarianceXYZPxPyPz(Double_t cv[21]) const = 0;

  ClassDef(AliVTrack,1)  // base class for tracks
};

#endif
