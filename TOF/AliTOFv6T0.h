#ifndef ALITOFv6T0_H
#define ALITOFv6T0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________//
//                                                                         //
// Implementation version v6 of TOF Manager class                          //
// FULL COVERAGE VERSION + OPTION FOR PHOS HOLES                           //
//                                                                         //
// -- Authors: G. Cara Romeo, A. De Caro                                   //
//                                                                         //
//_________________________________________________________________________//

#include "AliTOF.h"
 
 
class AliTOFv6T0 : public AliTOF {

public:
  AliTOFv6T0();
  AliTOFv6T0(const char *name, const char *title);
  virtual ~AliTOFv6T0() {};
  virtual void   BuildGeometry();
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 8;}
  virtual void   AddAlignableVolumes() const;
  virtual void   TOFpc(Float_t xtof,  Float_t ytof, Float_t zlenA);
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t) {};
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) {};
  virtual void   StepManager();
  virtual void   DrawModule() const;
  virtual void   DrawDetectorModules() const;
  virtual void   DrawDetectorStrips() const;
 
 protected:

  void MaterialMixer(Float_t* p,Float_t* a,Float_t* m,Int_t n) const;

private:

  void CreateModules(Float_t xtof, Float_t ytof, Float_t zlenA,
		     Float_t xFLT, Float_t yFLT, Float_t zFLTA) const;
  void MakeStripsInModules(Float_t ytof, Float_t zlenA) const;
  void CreateModuleCovers(Float_t xtof, Float_t zlenA) const;
  void CreateBackZone(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void MakeFrontEndElectronics() const;
  void MakeFEACooling(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void MakeNinoMask(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void MakeSuperModuleCooling(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void MakeSuperModuleServices(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void MakeModulesInBTOFvolumes(Float_t ytof, Float_t zlenA) const;
  void MakeCoversInBTOFvolumes() const;
  void MakeBackInBTOFvolumes(Float_t ytof) const;
  void MakeReadoutCrates(Float_t ytof) const;

  Int_t fIdFTOA; // FTOA volume identifier (outer plate A)
  Int_t fIdFTOB; // FTOB volume identifier (outer plate B)
  Int_t fIdFTOC; // FTOC volume identifier (outer plate C)
  Int_t fIdFLTA; // FLTA volume identifier (inner plate A)
  Int_t fIdFLTB; // FLTB volume identifier (inner plate B)
  Int_t fIdFLTC; // FLTC volume identifier (inner plate C)
  Bool_t fTOFHoles; // Selecting Geometry with and w/o holes
 
  //private:

  static const Float_t fgkModuleWallThickness;  // wall thickness (cm)
  static const Float_t fgkInterCentrModBorder1; // 1st distance of
						// border between
						// central and
						// intermediate
						// modules respect to
						// the central module
						// centre (cm)
  static const Float_t fgkInterCentrModBorder2; // 2nd distance of
						// border between the
						// central and
						// intermediate
						// modules respect to
						// the central module
						// centre (cm)
  static const Float_t fgkExterInterModBorder1; // 1st distance of
						// border between the
						// intermediate and
						// external modules
						// respect to the
						// central module
						// centre (cm)
  static const Float_t fgkExterInterModBorder2; // 2nd distance of
						// border between the
						// intermediate and
						// external
						// modules respect to
						// the central module
						// centre (cm)
  static const Float_t fgkLengthInCeModBorder;  // height of border
						// between the central
						// and intermediate
						// modules (cm)
  static const Float_t fgkLengthExInModBorder;  // height of border
						// between the
						// intermediate and
						// external modules
						// (cm)
  static const Float_t fgkModuleCoverThickness; // thickness of cover
						// modules zone (cm)
  static const Float_t fgkFEAwidth1; // mother volume width of each of
				     // two external FEA in a
				     // supermodule (cm)
  static const Float_t fgkFEAwidth2; // mother volume width of two
				     // internal FEA in a supermodule
				     // (cm)
  static const Float_t fgkSawThickness; // services alluminium wall
					// thickness (cm)
  static const Float_t fgkCBLw;  // cables&tubes block width (cm)
  static const Float_t fgkCBLh1; // min. height of cables&tubes block
				 // (cm)
  static const Float_t fgkCBLh2; // max. height of cables&tubes block
				 // (cm)

  ClassDef(AliTOFv6T0,0)  //Time Of Flight version 6
};
 
#endif /* ALITOFv6T0_H */
