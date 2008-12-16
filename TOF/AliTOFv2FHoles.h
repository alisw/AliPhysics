//_________________________________________________________________________
// Implementation version v2 of TOF Manager class
// HOLES FOR PHOS AND HMPID (RICH) DETECTOR
// NOT Official code!
// Requested by Andreas Morsch in order v2 TOF version
// is compliant with Frame v0 (Frame with holes)  
//*-- 
//*-- Authors: Pierella, Seganti, Vicinanza (Bologna and Salerno University)

#ifndef ALITOFv2FHoles_H
#define ALITOFv2FHoles_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

 
#include "AliTOF.h"
 
 
class AliTOFv2FHoles : public AliTOF {
 
public:
  AliTOFv2FHoles();
  AliTOFv2FHoles(const char *name, const char *title);
  virtual ~AliTOFv2FHoles() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 2;}
  virtual void   TOFpc(Float_t xtof,Float_t ytof,Float_t zlenC,Float_t zlenB,
                       Float_t zlenA,Float_t ztof0);
  virtual void   StepManager();
  virtual void   DrawModule() const;
  virtual void   DrawDetectorModules();
  virtual void   DrawDetectorStrips();
//  virtual void   DrawDetectorModulesinFrame();
//  virtual void   DrawDetectorStripsinFrame();

private:
  Int_t fIdFTOA; // FTOA volume identifier (outer plate A)
  Int_t fIdFTOB; // FTOB volume identifier (outer plate B)
  Int_t fIdFTOC; // FTOC volume identifier (outer plate C)
  Int_t fIdFLTA; // FLTA volume identifier (inner plate A)
  Int_t fIdFLTB; // FLTB volume identifier (inner plate B)
  Int_t fIdFLTC; // FLTC volume identifier (inner plate C)
  
   ClassDef(AliTOFv2FHoles,1)  //Time Of Flight version 2
};
 
#endif /* ALITOFv2FHoles_H  */
