#ifndef ALIPHOS_H
#define ALIPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

////////////////////////////////////////////////
//   Abstract Base Class for PHOS             //
//  Version SUBATECH                          //
//  Author  Laurent Aphecetche SUBATECH       //
//   The only provided method here is         // 
//   CreateMaterials, which defines the       // 
//   materials common to all PHOS versions.   // 
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

// --- AliRoot header files ---

#include "AliDetector.h"
#include "AliPHOSGeometry.h" 


class AliPHOS : public AliDetector {

 public:

  AliPHOS(const char* name, const char* title) ;
  AliPHOS() ;
  virtual ~AliPHOS() ; 
 
  virtual void CreateMaterials() ;
  virtual AliPHOSGeometry *  GetGeometry() = 0 ; 

  ClassDef(AliPHOS,2) // Photon Spectrometer Detector

} ;

#endif // ALIPHOS_H
