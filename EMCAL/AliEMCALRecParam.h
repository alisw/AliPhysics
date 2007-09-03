#ifndef ALIEMCALRECPARAM_H
#define ALIEMCALRECPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------------
// Container of EMCAL reconstruction parameters
// The purpose of this object is to store it to OCDB
// and retrieve it in AliEMCALClusterizerv1
// Author: Yuri Kharlov
//-----------------------------------------------------------------------------

// --- ROOT system ---

#include "TObject.h" 

class AliEMCALRecParam : public TObject
{
public:
  
  AliEMCALRecParam() ;
  virtual ~AliEMCALRecParam() {}
  Float_t GetClusteringThreshold() const     {return fClusteringThreshold;}
  Float_t GetW0                 () const     {return fW0                 ;}
  Float_t GetMinECut            () const     {return fMinECut            ;}
  void SetClusteringThreshold(Float_t thrsh)   {fClusteringThreshold = thrsh;}
  void SetW0                 (Float_t w0)      {fW0 = w0                    ;}
  void SetMinECut            (Float_t minEcut) {fMinECut = minEcut          ;}
  virtual void Print(Option_t * option="") const ; 

private:
  Float_t fClusteringThreshold ; // minimum energy to seed a EC digit in a cluster
  Float_t fW0 ;                  // logarithmic weight for the cluster center of gravity calculation
  Float_t fMinECut;              // Minimum energy for a digit to be a member of a cluster

  ClassDef(AliEMCALRecParam,1)   // Reconstruction parameters

} ;

#endif //  ALIEMCALRECPARAM_H
