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
  Float_t InsideModulePositionX(Int_t i) const;
  Float_t InsideModulePositionY(Int_t i) const;
  Float_t InsideModulePositionZ(Int_t i) const;
  Float_t CenterModulePositionX(Int_t i) const;
  Float_t CenterModulePositionY(Int_t i) const;
  Float_t CenterModulePositionZ(Int_t i) const;
  Float_t OutsideModulePositionX(Int_t i) const;
  Float_t OutsideModulePositionY(Int_t i) const;
  Float_t OutsideModulePositionZ(Int_t i) const;
  Float_t OldModulePositionX(Int_t i) const;
  Float_t OldModulePositionY(Int_t i) const;
  Float_t OldModulePositionZ(Int_t i) const;
  Float_t SupportModulePositionX(Int_t i) const;
  Float_t SupportModulePositionY(Int_t i) const;
  Float_t SupportModulePositionZ(Int_t i) const;
  Float_t OldExtraModulePositionZ(Int_t i) const;
  Int_t OldModuleElectronicChannel(Int_t i) const;
  Float_t OldExtraModulePositionX() const;
  Float_t OldExtraModulePositionY() const;
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
  static const Float_t fgkInsideModulePositionX[60]; 
  static const Float_t fgkInsideModulePositionY[60];
  static const Float_t fgkInsideModulePositionZ[60];
  static const Float_t fgkCenterModulePositionX[60];
  static const Float_t fgkCenterModulePositionY[60];
  static const Float_t fgkCenterModulePositionZ[60];
  static const Float_t fgkOutsideModulePositionX[60];
  static const Float_t fgkOutsideModulePositionY[60];
  static const Float_t fgkOutsideModulePositionZ[60];
  static const Float_t fgkOldModulePositionX[60]; // OLD position in ALICE
  static const Float_t fgkOldModulePositionY[60]; // of center of module
  static const Float_t fgkOldModulePositionZ[60]; 
  static const Float_t fgkSupportModulePositionX[60];
  static const Float_t fgkSupportModulePositionY[60];
  static const Float_t fgkSupportModulePositionZ[60];
  static const Float_t fgkOldExtraModulePositionZ[4];
  static const Float_t fgkOldExtraModulePositionX;
  static const Float_t fgkOldExtraModulePositionY;
  static const Int_t fgkOldModuleElectronicChannel[60];

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
