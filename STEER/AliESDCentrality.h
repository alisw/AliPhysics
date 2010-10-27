//-*- Mode: C++ -*-
#ifndef ALIESDCentrality_H
#define ALIESDCentrality_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliCentralitySelectionTask
//   author: Alberica Toia
//*****************************************************

#include "TNamed.h"

class AliESDCentrality : public TNamed
{
 public:

  AliESDCentrality();  /// constructor
  ~AliESDCentrality();  /// destructor
  AliESDCentrality(const AliESDCentrality& cnt); /// copy constructor
  AliESDCentrality& operator=(const AliESDCentrality& cnt);   /// assignment operator

  /// set centrality result
  void SetCentrality(Float_t cent) {fCentrality = cent;}

  /// get centrality result
  Float_t GetCentralityPercentile();
  Int_t   GetCentralityClass10();
  Int_t   GetCentralityClass5();
  Bool_t  IsEventInCentralityClass(Float_t a, Float_t b);

 private:
  Float_t fCentrality;   // Centrality

  ClassDef(AliESDCentrality, 1)
};
#endif //ALIESDCENTRALITY_H
