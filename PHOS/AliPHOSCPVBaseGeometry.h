#ifndef ALIPHOSCPVBASEGEOMETRY_H
#define ALIPHOSCPVBASEGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry base class for PHOS:CVS (Charged particle veto)
//
//*-- Author  : Yuri Kharlov (IHEP, Protvino)
//    14 September 2000

#include "TObject.h"

class AliPHOSCPVBaseGeometry : public TObject {

public: 

           AliPHOSCPVBaseGeometry()     {}
  virtual ~AliPHOSCPVBaseGeometry(void) {}

  // Return common for PPSD and CPV geometrical parameters

  virtual Float_t GetCPVBoxSize(Int_t index)        { return 0; }

  // Return PPSD geometrical parameters

  virtual Float_t GetAnodeThickness(void)           { return 0; }
  virtual Float_t GetAvalancheGap(void)             { return 0; }
  virtual Float_t GetCathodeThickness(void)         { return 0; }
  virtual Float_t GetCompositeThickness(void)       { return 0; }
  virtual Float_t GetConversionGap(void)            { return 0; }
  virtual Float_t GetLeadConverterThickness(void)   { return 0; }
  virtual Float_t GetLeadToMicro2Gap(void)          { return 0; }
  virtual Float_t GetLidThickness(void)             { return 0; }
  virtual Float_t GetMicromegas1Thickness(void)     { return 0; }
  virtual Float_t GetMicromegas2Thickness(void)     { return 0; }
  virtual Float_t GetMicromegasWallThickness(void)  { return 0; }
  virtual Float_t GetMicro1ToLeadGap(void)          { return 0; }
  virtual Int_t   GetNumberOfPadsPhi(void)          { return 0; }
  virtual Int_t   GetNumberOfPadsZ(void)            { return 0; }
  virtual Int_t   GetNumberOfModulesPhi(void)       { return 0; }
  virtual Int_t   GetNumberOfModulesZ(void)         { return 0; }
  virtual Float_t GetPCThickness(void)              { return 0; }
  virtual Float_t GetPhiDisplacement(void)          { return 0; }
  virtual Float_t GetPPSDModuleSize(Int_t index)    { return 0; }
  virtual Float_t GetZDisplacement(void)            { return 0; }

  // Return CPV geometrical parameters

  virtual Int_t   GetNumberOfCPVLayers(void)        { return 0; }
  virtual Bool_t  IsLeadConverterExists(void)       { return 0; }
  virtual Float_t GetCPVActiveSize(Int_t index)     { return 0; }
  virtual Int_t   GetNumberOfCPVChipsPhi(void)      { return 0; }
  virtual Int_t   GetNumberOfCPVChipsZ(void)        { return 0; }
  virtual Float_t GetGassiplexChipSize(Int_t index) { return 0; }
  virtual Float_t GetCPVGasThickness(void)          { return 0; }
  virtual Float_t GetCPVTextoliteThickness(void)    { return 0; }
  virtual Float_t GetCPVCuNiFoilThickness(void)     { return 0; }
  virtual Float_t GetFTPosition(Int_t index)        { return 0; }
  virtual Float_t GetCPVFrameSize(Int_t index)      { return 0; }
 
  ClassDef(AliPHOSCPVBaseGeometry,1)        // CPV base geometry class 

} ;

#endif // ALIPHOSCPVBASEGEOMETRY_H
