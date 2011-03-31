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

class AliESDVertex;
class TObjArray;
class AliExternalTrackParam;
class AliAODVertex;

class AliHLTD0toKpi : public TObject
{
public:
  AliHLTD0toKpi();

  Double_t InvMass(const AliExternalTrackParam* d1, const AliExternalTrackParam* d2);
  void CosThetaStar(const AliExternalTrackParam* n, const AliExternalTrackParam* p,Double_t &D0,Double_t &D0bar);
  Double_t PointingAngle(const AliExternalTrackParam* n, const AliExternalTrackParam* p, const Double_t *pv, const Double_t *sv);
  Double_t Pt(const AliExternalTrackParam* d1, const AliExternalTrackParam* d2);

  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray, Double_t b, const AliESDVertex *v, bool useKF);

private:
  
  ClassDef(AliHLTD0toKpi, 1)
};
#endif //ALIHLTD0TOKPI_H
