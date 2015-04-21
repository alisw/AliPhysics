#ifndef ALIPI0MIXEDDICLUSTERSTATSTRUCT_H
#define ALIPI0MIXEDDICLUSTERSTATSTRUCT_H

#include <TObject.h>
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */
/// $Id$
/// \class Alipi0mixedDiClusterStatStruct
/// \brief mixed di cluster information structure pion/eta in PbPb (not used.)
///
/// \author Astrid Morreale astridmorreale@cern.ch subatech
/// \date April 10 2015
class mixedDiClusterStatStruct: public TObject
{

  public:

  /// object name (re-implemented)
  virtual const char*	GetName() const
  { return "mixeddiClusterStatStruct"; }

  /// default contructor
  mixedDiClusterStatStruct( void ):
  isMBmx(kFALSE),
  isAnyINTmx(kFALSE),
  isCentralmx(kFALSE),
  isSemiCentralmx(kFALSE),
  isEgamx(kFALSE),
  kAllMBmx(kFALSE),
  Mxmass(0),
  Mxpt(0),
  centMixedV0(0),
  centMixedSPD(0)

{}



  Bool_t  isMBmx;
  Bool_t  isAnyINTmx;
  Bool_t  isCentralmx;
  Bool_t  isSemiCentralmx;
  Bool_t  isEgamx;
  Bool_t  kAllMBmx;
  Float_t Mxmass;
  Float_t Mxpt;
  Float_t centMixedV0;
  Float_t centMixedSPD;

  ClassDef(mixedDiClusterStatStruct,1)

};

#endif
