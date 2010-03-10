/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTEveCalo.h
/// @author Svein Lindal
/// @brief  EMCAL Instance of Eve display processor

#ifndef ALIHLTEVEEMCAL_H
#define ALIHLTEVEEMCAL_H

#include "AliHLTEveCalo.h"
class TEveElementList;
class AliEMCALGeoUtils;

class AliHLTEveEmcal : public AliHLTEveCalo {

public:
  
  /** Constructor  **/
  AliHLTEveEmcal();

  /** Destructor **/
 ~AliHLTEveEmcal();
  
private:

  /** copy constructor prohibited */
  AliHLTEveEmcal(const AliHLTEveEmcal&);
  /** assignment operator prohibited */
  AliHLTEveEmcal& operator = (const AliHLTEveEmcal );

  void AddClusters(Float_t * pos, Int_t module, Float_t energy);

  void AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy);
  
  TEveElementList * CreateElementList();

  AliEMCALGeoUtils * fGeoUtils;

  ClassDef(AliHLTEveEmcal, 0);
};

#endif
