#ifndef ALIPMDUTILITY_H
#define ALIPMDUTILITY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Utility class for PMD                              //
//                                                     //
//-----------------------------------------------------//
#include "Rtypes.h"
class AliPMDUtility
{
 public:
  AliPMDUtility();
  AliPMDUtility(Float_t px, Float_t py, Float_t pz);
  virtual ~AliPMDUtility();

  void RectGeomCellPos(Int_t ism, Int_t xpad, Int_t ypad,
		       Float_t & xpos, Float_t & ypos);
  void RectGeomCellPos(Int_t ism, Float_t xpad, Float_t ypad,
		       Float_t & xpos, Float_t & ypos);
  void SetPxPyPz(Float_t px, Float_t py, Float_t pz);
  void SetXYZ(Float_t xpos, Float_t ypos, Float_t zpos);
  void CalculateEta();
  void CalculatePhi();
  void CalculateEtaPhi();
  Float_t GetTheta() const;
  Float_t GetEta() const;
  Float_t GetPhi() const;
  
 protected:
  Float_t fPx;     // Momentum along x
  Float_t fPy;     // Momentum along y
  Float_t fPz;     // Momentum along z
  Float_t fTheta;  // Polar angle in radian
  Float_t fEta;    // Pseudo-rapidity
  Float_t fPhi;    // Azimuthal angle in radian
  
  ClassDef(AliPMDUtility,3) // Utility class for the detector set:PMD
};

#endif
