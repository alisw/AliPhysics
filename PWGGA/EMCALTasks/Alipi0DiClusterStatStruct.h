#ifndef ALIPI0DICLUSTERSTATSTRUCT_H
#define ALIPI0DICLUSTERSTATSTRUCT_H

#include <TObject.h>
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */
/// $Id$
/// \class Alipi0DiClusterStatStruct
/// \brief di cluster information structure pion/eta in PbPb
///
/// \author Astrid Morreale astridmorreale@cern.ch subatech
/// \date April 10 2015

class diClusterStatStruct: public TObject
{

  public:

  /// object name (re-implemented)
  virtual const char*	GetName() const
  { return "diclusterStatStruct"; }

  /// default contructor
  diClusterStatStruct( void ):
 NCellFlag(kFALSE),  EcFlag(kFALSE), M02Flag(kFALSE), D2BadChFlag(kFALSE),
  piE(0),
  piphi(0),
  pieta(0),
  ptpi(0),
  pipx(0), pipy(0), pipz(0),
  asympi(0), masspi(0)
{}

  Bool_t    genFlag;
  Bool_t    NCellFlag;
  Bool_t    EcFlag;
  Bool_t    M02Flag;
  Bool_t    D2BadChFlag;
  Float_t   piE;
  Float_t   piphi;
  Float_t   pieta;
  Float_t   ptpi;
  Float_t   pipx;
  Float_t   pipy;
  Float_t   pipz;
  Float_t   asympi;
  Float_t   masspi;


  ClassDef(diClusterStatStruct,1)

};

#endif
