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
