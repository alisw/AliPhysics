#ifndef ALITOFGEOMETRYV5_H
#define ALITOFGEOMETRYV5_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF geometry class (new version)                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTOFGeometry.h"

class AliTOFGeometryV5: public AliTOFGeometry {

 public:
  AliTOFGeometryV5();
  virtual ~AliTOFGeometryV5();
  
  void    Init();
  Bool_t  IsInsideThePad(Int_t *det, Float_t *pos); 
  Float_t DistanceToPad(Int_t *det, Float_t *pos, Float_t *dist3d=0);
  Int_t   GetPlate(Float_t *pos);
  Int_t   GetStrip(Float_t *pos);
  Int_t   GetSector(Float_t *pos);
  Int_t   GetPadX(Float_t *pos);
  Int_t   GetPadZ(Float_t *pos);
  Float_t GetX(Int_t *det);
  Float_t GetY(Int_t *det);
  Float_t GetZ(Int_t *det);
  Float_t GetPadDx(Float_t *pos);
  Float_t GetPadDy(Float_t *pos);
  Float_t GetPadDz(Float_t *pos);
  //Float_t GetAngles(Int_t iplate, Int_t istrip)  const {return fAngles[iplate][istrip];};
  //Float_t GetHeights(Int_t iplate, Int_t istrip) const {return fHeights[iplate][istrip];};
  //Float_t GetDistances(Int_t iplate, Int_t istrip) const {return fDistances[iplate][istrip];};

  Float_t NStirpC()     { return kNStripC;};
  Int_t   NMaxNstrip()  { return kMaxNstrip;};
  Int_t   NPadXSector() { return (AliTOFGeometry::kNStripA + 2*AliTOFGeometry::kNStripB +
				  2*kNStripC)*AliTOFGeometry::kNpadX*AliTOFGeometry::kNpadZ;};

  Float_t RinTOF()      { return fgkxTOF;};
  Float_t Rmin()        { return fgkRmin;};
  Float_t Rmax()        { return fgkRmax;};

  Float_t ZlenA()       { return fgkZlenA;};
  Float_t ZlenB()       { return fgkZlenB;};
  Float_t ZlenC()       { return fgkZlenC;};
  Float_t MaxhZtof()    { return fgkMaxhZtof;};
  Float_t StripLength() { return fgkStripLength;};

  void Translation(Float_t *xyz, Float_t translationVector[3]);
  void Rotation(Float_t *xyz, Double_t rotationAngles[6]);
  void InverseRotation(Float_t *xyz, Double_t rotationAngles[6]);

  protected:

  //private:

  static const Int_t kNStripC;         // number of strips in C type module 
  static const Int_t kMaxNstrip;       // Max. number of strips

  static const Float_t fgkZlenA;       // length (cm) of the A module
  static const Float_t fgkZlenB;       // length (cm) of the B module
  static const Float_t fgkZlenC;       // length (cm) of the C module
  static const Float_t fgkMaxhZtof;    // Max half z-size of TOF (cm)
  static const Float_t fgkStripLength; // Strip Length (rho X phi direction) (cm)

  static const Float_t fgkRmin;        // Inner radius of the TOF (cm)
  static const Float_t fgkRmax;        // Outer radius of the TOF (cm)
  static const Float_t fgkxTOF;        // Inner TOF Radius used in Reconstruction (cm)

  Bool_t fHoles; //logical for geometry version (w/wo holes) 

  //Float_t *fAngles[kNPlates];//Strip tilt angles
  //Float_t *fHeights[kNPlates];//Strip heights
  //Float_t *fDistances[kNPlates];//Strip distances

  ClassDef(AliTOFGeometryV5,0) // TOF Geometry class
};

#endif
