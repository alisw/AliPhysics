//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTD0TOKPI_H
#define ALIHLTD0TOKPI_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTD0toKpi.h
/// @author Gaute Ovrebekk
/// @date   2009-10-28
/// @brief  Class for calculating D0->Kpi

#include "TObject.h"

class AliESDtrack;
class AliAODVertex;
class AliESDVertex;
class TObjArray;

class AliHLTD0toKpi : public TObject
{
public:
  AliHLTD0toKpi();
  Double_t InvMass(AliESDtrack* d1, AliESDtrack* d2);
  void cosThetaStar(AliESDtrack* n, AliESDtrack* p,Double_t &D0,Double_t &D0bar);
  Double_t pointingAngle(AliESDtrack* n, AliESDtrack* p, Double_t *pv, Double_t *sv);

  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray, Double_t b, const AliESDVertex *v);

private:
  
  ClassDef(AliHLTD0toKpi, 1)
};
#endif //ALIHLTD0TOKPI_H
