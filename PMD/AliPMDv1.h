#ifndef ALIPMDV1_H
#define ALIPMDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////

class TGeoManager;
 
#include "AliPMD.h"

//___________________________________________
 
class AliPMDv1 : public AliPMD {
  
public:
  AliPMDv1();
  AliPMDv1(const char *name, const char *title);
  virtual      ~AliPMDv1() {}
  virtual void  CreateGeometry();
  virtual void  CreatePMD();
  virtual void  CreateSupermodule();
  virtual void  GetParameters();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
  virtual void  AddAlignableVolumes() const;
  void          SetSectorAlignable() const;
  void          SetCpvOff();
  void          SetPreOff();
  void          SetModuleOff(Int_t imodule);

private:

  static const Int_t   fgkNcolUM1 = 48;  // Number of cols in UM, type 1
  static const Int_t   fgkNcolUM2 = 96;  // Number of cols in UM, type 2
  static const Int_t   fgkNrowUM1 = 96;  // Number of rows in UM, type 1
  static const Int_t   fgkNrowUM2 = 48;  // Number of rows in UM, type 2
  static const Float_t fgkCellRadius = 0.25;   // Radius of a hexagonal cell
  static const Float_t fgkCellWall   = 0.02;   // Thickness of cell Wall
  static const Float_t fgkCellDepth  = 0.50;   // Gas thickness
  static const Float_t fgkThPCB      = 0.16;   // Thickness of PCB
  static const Float_t fgkThLead     = 1.5;    // Thickness of Pb
  static const Float_t fgkThSteel    = 0.5;    // Thickness of Steel
  static const Float_t fgkGap        = 0.025;  // Air Gap
  static const Float_t fgkZdist      = 361.5;   // z-position of the detector
  static const Float_t fgkSqroot3    = 1.7320508;   // Square Root of 3
  static const Float_t fgkSqroot3by2 = 0.8660254;   // Square Root of 3 by 2
  static const Float_t fgkSSBoundary = 0.3;  // Stainless steel boundary
  static const Float_t fgkThSS       = 1.23; // thickness of SS frame
  static const Float_t fgkThTopG10   = 0.33; // top G10 material thickness
  static const Float_t fgkThBotG10   = 0.4;  // bottom G10 material thickness

  Int_t   fModStatus[48];    // To position different modules

  Float_t fSMthick;     // Thickness of the full PMD profile
  Float_t fSMthickpmd;  // Thickness of the PMD detector only
  Float_t fDthick;      // Thickness of the pre/veto module
  Float_t fSMLengthax;  // Supermodule length along X, type A
  Float_t fSMLengthay;  // Supermodule length along Y, type A
  Float_t fSMLengthbx;  // Supermodule length along X, type B
  Float_t fSMLengthby;  // Supermodule length along Y, type A
  Int_t   fMedSens;     // Sensitive Medium Ar+CO2
  Float_t fDboxmm1[3];  // Master MODULE EMPA of aluminum for PMD
  Float_t fDboxmm12[3]; // Master MODULE EMCA of aluminum for CPV
  Float_t fDboxmm2[3];  // Master MODULE EMPB of aluminum for PMD
  Float_t fDboxmm22[3]; // Master MODULE EMCB of aluminum for CPV

 
  ClassDef(AliPMDv1,5)     //Hits manager for set:PMD
};
 
#endif
