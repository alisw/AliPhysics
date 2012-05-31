#ifndef ALITOFV6T0_H
#define ALITOFV6T0_H
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
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 8;}
  virtual void   AddAlignableVolumes() const;
  virtual void   TOFpc(Float_t xtof,  Float_t ytof, Float_t zlenA);
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t) {};
  virtual void   TOFpc(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) {};
  virtual void   StepManager();
 
 protected:

  void MaterialMixer(Float_t * p, const Float_t * const a,
		     const Float_t * const m, Int_t n) const;

private:

  void CreateModules(Float_t xtof, Float_t ytof, Float_t zlenA,
		     Float_t xFLT, Float_t yFLT, Float_t zFLTA) const;
  void MakeStripsInModules(Float_t ytof, Float_t zlenA) const;
  void CreateModuleCovers(Float_t xtof, Float_t zlenA) const;
  void CreateBackZone(Float_t xtof, Float_t ytof, Float_t zlenA) const;
  void MakeFrontEndElectronics(Float_t xtof) const;
  void MakeFEACooling(Float_t xtof) const;
  void MakeNinoMask(Float_t xtof) const;
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

  //private:

  static const Bool_t fgkFEAwithMasks[18]; // Selecting TOF sectors containing FEA cooling masks

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
  /*
  static const Float_t fgkLengthInCeModBorder;  // height of border
						// between the central
						// and intermediate
						// modules (cm)
						*/
  static const Float_t fgkLengthInCeModBorderU; // height of upper border
						// between the central
						// and intermediate
						// modules (cm)
  static const Float_t fgkLengthInCeModBorderD; // height of lower border
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
  static const Float_t fgkBetweenLandMask; // distance between the L
					   // element and the Nino
					   // mask (cm)
  static const Float_t fgkAl1parameters[3]; // (cm)
  static const Float_t fgkAl2parameters[3]; // (cm)
  static const Float_t fgkAl3parameters[3]; // (cm)
  static const Float_t fgkRoof1parameters[3]; // (cm)
  static const Float_t fgkRoof2parameters[3]; // (cm)
  static const Float_t fgkFEAparameters[3]; // (cm)
  //static const Float_t fgkFCAparameters[3]; // (cm)
  static const Float_t fgkBar[3]; // (cm)
  static const Float_t fgkBar1[3]; // (cm)
  static const Float_t fgkBar2[3]; // (cm)
  static const Float_t fgkBarS[3]; // (cm)
  static const Float_t fgkBarS1[3]; // (cm)
  static const Float_t fgkBarS2[3]; // (cm)

  ClassDef(AliTOFv6T0,1)  //Time Of Flight version 6
};
 
#endif /* ALITOFV6T0_H */
