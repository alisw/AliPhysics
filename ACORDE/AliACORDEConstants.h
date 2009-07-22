#ifndef ALIACORDECONSTANTS_H
#define ALIACORDECONSTANTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliACORDEConstants class
//
// This class serves to group constants needed by ACORDE detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
//
// Author: Arturo Fernandez, Enrique Gamez, Mario Rodriguez Cahuantzi
//         FCFM-UAP, Mexico.
//
// Last update: Nov. 24th 08
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>

enum ECRMode {
  kSingleMuons,
  kMuonBundle,
  kMuonFlux
};

class AliACORDEConstants : public TObject {
public:
  virtual ~AliACORDEConstants();

  static AliACORDEConstants* Instance();

  // constant for geometry
  Float_t ModuleLength() const;
  Float_t ModuleWidth() const;
  Float_t ModuleHeight() const;
  Float_t ModulePositionX(Int_t i) const;
  Float_t ModulePositionY(Int_t i) const;
  Float_t ModulePositionZ(Int_t i) const;
  Float_t SupportModulePositionX(Int_t i) const;
  Float_t SupportModulePositionY(Int_t i) const;
  Float_t SupportModulePositionZ(Int_t i) const;
  Float_t ExtraModulePositionZ(Int_t i) const;
  Int_t ModuleElectronicChannel(Int_t i) const;
  Float_t ExtraModulePositionX() const;
  Float_t ExtraModulePositionY() const;
  Float_t PlasticLength() const;
  Float_t PlasticWidth() const;
  Float_t PlasticHeight() const;
  Float_t ProfileWidth() const;
  Float_t ProfileThickness() const;
  Float_t Depth() const;

  // constant to convert hits in digits
  Float_t HitEnergyThreshold() const { return fgkHitEnergyThreshold;}
  Float_t MaxHitTimeDifference() const { return fgkMaxHitTimeDifference;}
  // constants for trigger
  Int_t MultiMuonThreshold() const { return fgkMultiMuonThreshold;}
  Float_t MultiMuonWindow() const { return fgkMultiMuonWindow;}

protected:

  AliACORDEConstants();

  static AliACORDEConstants* fgInstance; // static instanton

  static const Float_t fgkModuleLength; // Module lenght
  static const Float_t fgkModuleWidth;  // Module width
  static const Float_t fgkModuleHeight; // Module height
  static const Float_t fgkModulePositionX[60]; // position in ALICE
  static const Float_t fgkModulePositionY[60]; // of center of module
  static const Float_t fgkModulePositionZ[60]; 
  static const Float_t fgkSupportModulePositionX[60];
  static const Float_t fgkSupportModulePositionY[60];
  static const Float_t fgkSupportModulePositionZ[60];
  static const Float_t fgkExtraModulePositionZ[4];
  static const Float_t fgkExtraModulePositionX;
  static const Float_t fgkExtraModulePositionY;
  static const Int_t fgkModuleElectronicChannel[60];

  static const Float_t fgkPlasticLength; // Plastic length
  static const Float_t fgkPlasticWidth;  // Plastic width
  static const Float_t fgkPlasticHeight; // Plastic height


  static const Float_t fgkProfileWidth;   // profile of the module
  static const Float_t fgkProfileThickness; 
  static const Float_t fgkDepth; // Alice IP depth from surface

  static const Float_t fgkHitEnergyThreshold;
  static const Float_t fgkMaxHitTimeDifference;
  static const Int_t fgkMultiMuonThreshold;
  static const Float_t fgkMultiMuonWindow;


 private:
  ClassDef(AliACORDEConstants, 0)   // ACORDE(ACORDE) global constants
};

typedef AliACORDEConstants AliCRTConstants; // for backward compatibility

#endif // ALIACORDECONSTANTS_H
