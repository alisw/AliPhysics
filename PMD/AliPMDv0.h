#ifndef ALIPMDV0_H
#define ALIPMDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliPMD.h"

//___________________________________________
 
class AliPMDv0 : public AliPMD {

public:
  AliPMDv0();
  AliPMDv0(const char *name, const char *title);
  virtual      ~AliPMDv0() {}
  virtual void  CreateGeometry();
  virtual void  CreatePMD();
  virtual void  CreateSupermodule();
  virtual void  GetParameters();
  virtual void  CreateMaterials();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  StepManager();
  
 private:
  static const Int_t   fgkNcellHole;     // Hole Dimension
  static const Float_t fgkCellRadius;    // Radius of a hexagonal cell
  static const Float_t fgkCellWall;      // Thickness of cell Wall
  static const Float_t fgkCellDepth;     // Gas thickness
  static const Float_t fgkBoundary;      // Thickness of Boundary wall
  static const Float_t fgkThBase;        // Thickness of Base plate
  static const Float_t fgkThAir;         // Thickness of Air
  static const Float_t fgkThPCB;         // Thickness of PCB
  static const Float_t fgkThLead;        // Thickness of Pb
  static const Float_t fgkThSteel;       // Thickness of Steel
  static const Float_t fgkZdist;         // z-position of the detector
  static const Float_t fgkSqroot3;       // Square Root of 3
  static const Float_t fgkSqroot3by2;    // Square Root of 3 by 2 
  static const Float_t fgkPi;            // Value of pi
  
  Float_t fSMthick;     // Thickness of the supermodule
  Float_t fSMLength;    // Supermodule length
  Int_t   fMedSens;     // Sensitive Medium (Ar+C02)
  Int_t   fNcellSM;     // Number of cells in SuperModule
  
  ClassDef(AliPMDv0,1)  //Hits manager for set:PMD
};
 
#endif


