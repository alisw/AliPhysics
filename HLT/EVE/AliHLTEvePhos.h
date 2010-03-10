/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTEveCalo.h
/// @author Svein Lindal
/// @brief  PHOS Instance of Eve display processor

#ifndef ALIHLTEVEPHOS_H
#define ALIHLTEVEPHOS_H

#include "AliHLTEveCalo.h"
class TEveElementList;

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
  AliHLTEvePhos& operator = (const AliHLTEvePhos );

  /** inherited from AliHLTEveCalo */
  TEveElementList * CreateElementList();
  
  /** inherited from AliHLTEveCalo */
  void AddClusters(Float_t * pos, Int_t module, Float_t energy);

  /** inherited from AliHLTEveCalo */
  void AddDigits(UShort_t fX, UShort_t fZ, Int_t module, Float_t energy);

  ClassDef(AliHLTEvePhos, 0);
};

#endif
