#ifndef ALIPHOS_H
#define ALIPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for PHOS     
//                  
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- AliRoot header files ---

#include "AliDetector.h"
#include "AliPHOSGeometry.h" 


class AliPHOS : public AliDetector {

 public:

  AliPHOS(const char* name, const char* title): AliDetector(name,title) {} 
  AliPHOS() : AliDetector() {} 
  virtual ~AliPHOS() ; 
 
  virtual void CreateMaterials() ;               // defines the material of the detector
  virtual AliPHOSGeometry *  GetGeometry() = 0 ; // hands the pointer to the unique geometry object

  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
