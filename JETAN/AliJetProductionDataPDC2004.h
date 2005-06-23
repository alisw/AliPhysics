#ifndef ALIJETPRODUCTIONDATAPDC2004_H
#define ALIJETPRODUCTIONDATAPDC2004_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Service class for jet production data
// Physics Data Challenge 2004
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//---------------------------------------------------------------------
 
#include <AliJetProductionData.h>
 
class AliJetProductionDataPDC2004 : public AliJetProductionData
{
 public:
   AliJetProductionDataPDC2004();
  ~AliJetProductionDataPDC2004();
 private:
  ClassDef(AliJetProductionDataPDC2004, 1)
};
 
#endif
