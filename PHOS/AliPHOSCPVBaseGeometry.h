#ifndef ALIPHOSCPVBASEGEOMETRY_H
#define ALIPHOSCPVBASEGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry base class for PHOS:CVS (Charged particle veto)
// It contains only virtual member functions-getters
// which will be implemented for PPSD and CPV in
// appropriate derived classes
// Author  : Yuri Kharlov (IHEP, Protvino)
// 14 September 2000

#include "TObject.h"

class AliPHOSCPVBaseGeometry : public TObject {

public: 

           AliPHOSCPVBaseGeometry()     {}
  virtual ~AliPHOSCPVBaseGeometry(void) {}

  // Return common for PPSD and CPV geometrical parameters

  virtual Float_t GetCPVBoxSize(Int_t index)        = 0 ;

  // Return PPSD geometrical parameters

  virtual Float_t GetAnodeThickness(void)           = 0 ;
  virtual Float_t GetAvalancheGap(void)             = 0 ;
  virtual Float_t GetCathodeThickness(void)         = 0 ;
  virtual Float_t GetCompositeThickness(void)       = 0 ;
  virtual Float_t GetConversionGap(void)            = 0 ;
  virtual Float_t GetLeadConverterThickness(void)   = 0 ;
  virtual Float_t GetLeadToMicro2Gap(void)          = 0 ;
  virtual Float_t GetLidThickness(void)             = 0 ;
  virtual Float_t GetMicromegas1Thickness(void)     = 0 ;
  virtual Float_t GetMicromegas2Thickness(void)     = 0 ;
  virtual Float_t GetMicromegasWallThickness(void)  = 0 ;
  virtual Float_t GetMicro1ToLeadGap(void)          = 0 ;
  virtual Int_t   GetNumberOfPadsPhi(void)          = 0 ;
  virtual Int_t   GetNumberOfPadsZ(void)            = 0 ;
  virtual Int_t   GetNumberOfModulesPhi(void)       = 0 ;
  virtual Int_t   GetNumberOfModulesZ(void)         = 0 ;
  virtual Float_t GetPCThickness(void)              = 0 ;
  virtual Float_t GetPhiDisplacement(void)          = 0 ;
  virtual Float_t GetPPSDModuleSize(Int_t index)    = 0 ;
  virtual Float_t GetZDisplacement(void)            = 0 ;
  virtual void    SetLeadConverterThickness(Float_t x) = 0 ;

  // Return CPV geometrical parameters

  virtual Int_t   GetNumberOfCPVLayers(void)        = 0 ;
  virtual Bool_t  IsLeadConverterExists(void)       = 0 ;
  virtual Float_t GetCPVActiveSize(Int_t index)     = 0 ;
  virtual Int_t   GetNumberOfCPVChipsPhi(void)      = 0 ;
  virtual Int_t   GetNumberOfCPVChipsZ(void)        = 0 ;
  virtual Float_t GetGassiplexChipSize(Int_t index) = 0 ;
  virtual Float_t GetCPVGasThickness(void)          = 0 ;
  virtual Float_t GetCPVTextoliteThickness(void)    = 0 ;
  virtual Float_t GetCPVCuNiFoilThickness(void)     = 0 ;
  virtual Float_t GetFTPosition(Int_t index)        = 0 ;
  virtual Float_t GetCPVFrameSize(Int_t index)      = 0 ;
 
  ClassDef(AliPHOSCPVBaseGeometry,1)        // CPV base geometry class 

} ;

#endif // ALIPHOSCPVBASEGEOMETRY_H
