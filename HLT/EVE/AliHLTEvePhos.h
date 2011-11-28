//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEPHOS_H
#define ALIHLTEVEPHOS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEvePhos.h
/// @author Svein Lindal
/// @brief  PHOS Instance of Eve display processor

#include "AliESDEvent.h"
#include "AliHLTEveCalo.h"

class TEveElementList;
class AliPHOSGeoUtils;

class AliHLTEvePhos : public AliHLTEveCalo {

public:
  
  /** Constructor  **/
  AliHLTEvePhos();

  /** Destructor **/
 ~AliHLTEvePhos();

private:

  /** copy constructor prohibited */
  AliHLTEvePhos(const AliHLTEvePhos&);
  /** assignment operator prohibited */
  AliHLTEvePhos& operator = (const AliHLTEvePhos& );

  /** inherited from AliHLTEveCalo */
  void CreateElementList();
  
  /** inherited from AliHLTEveCalo */
  void AddClusters(Float_t * pos, Int_t module, Float_t energy);
  void AddClusters(Float_t * pos, Int_t module, Float_t energy, Int_t nCells);

  /** inherited from AliHLTEveCalo */
  void AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy);

  Int_t GetClusters(AliESDEvent * event, TRefArray * clusters) { return event->GetPHOSClusters(clusters); }

  void ProcessESDCluster(AliESDCaloCluster * cluster);

  AliPHOSGeoUtils * fGeoUtils;  //PHOS geometry

  ClassDef(AliHLTEvePhos, 0);
};

#endif
