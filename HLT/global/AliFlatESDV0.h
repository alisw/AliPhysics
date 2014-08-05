#ifndef ALIFLATESDV0_H
#define ALIFLATESDV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * >> Flat structure representing a ESD v0 vertex <<
 */

#include "Rtypes.h"

class AliFlatESDV0
{
 public:

  AliFlatESDV0(): fNegTrackID(-1), fPosTrackID(-1) {}
  ~AliFlatESDV0(){}

  static size_t GetSize(){ return sizeof(AliFlatESDV0); }

  void SetNegTrackID( Int_t id ){ fNegTrackID = id; }
  void SetPosTrackID( Int_t id ){ fPosTrackID = id; }

  Int_t GetNegTrackID() const { return fNegTrackID; }
  Int_t GetPosTrackID() const { return fPosTrackID; }

 private:

  AliFlatESDV0(const AliFlatESDV0&);
  AliFlatESDV0& operator=(const AliFlatESDV0&);

  Int_t fNegTrackID;
  Int_t fPosTrackID;
};

#endif
