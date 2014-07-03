#ifndef ALIVVEXTERNALTRACKPARAM_H
#define ALIVVEXTERNALTRACKPARAM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Author: Mikolaj Krzewicki, mkrzewic@cern.ch     */

/**
 * >> interface class <<
 */

#include "Rtypes.h"

class AliVVexternalTrackParam
{
  public:
  AliVVexternalTrackParam() {}
  AliVVexternalTrackParam(Bool_t) {}
  virtual ~AliVVexternalTrackParam() {}
  virtual void SetAlpha(Float_t /*alpha*/)                 {}
  virtual void SetX(Float_t /*x*/)                         {}
  virtual void SetY(Float_t /*y*/)                         {}
  virtual void SetZ(Float_t /*z*/)                         {}
  virtual void SetSnp(Float_t /*snp*/)                     {}
  virtual void SetTgl(Float_t /*tgl*/)                     {}
  virtual void SetSigned1Pt(Float_t /*signed1Pt*/)         {}
  virtual void SetCovEntry(Int_t /*idx*/, Float_t /*cov*/) {}

  virtual Float_t  GetAlpha()             const {return 0.;}
  virtual Float_t  GetX()                 const {return 0.;}
  virtual Float_t  GetY()                 const {return 0.;}
  virtual Float_t  GetZ()                 const {return 0.;}
  virtual Float_t  GetSnp()               const {return 0.;}
  virtual Float_t  GetTgl()               const {return 0.;}
  virtual Float_t  GetSigned1Pt()         const {return 0.;}
  virtual Float_t* GetCov()               const {return NULL;}
  virtual Float_t  GetCovEntry(Int_t /*idx*/) const {return 0.;}
};

#endif
