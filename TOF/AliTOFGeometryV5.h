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
   
  void    ImportGeometry();
  void    Init();
  Bool_t  IsInsideThePad(Int_t *det, TGeoHMatrix mat, Float_t *pos) const;
  Float_t DistanceToPad(Int_t *det, TGeoHMatrix mat, Float_t *pos, Float_t *dist3d=0) const;
  Bool_t  IsInsideThePadPar(Int_t *det, Float_t *pos) const; 
  Float_t DistanceToPadPar(Int_t *det, Float_t *pos, Float_t *dist3d=0) const;
  void    GetVolumePath(Int_t *ind, Char_t *path );
  Int_t   GetPlate(Float_t *pos) const;
  Int_t   GetStrip(Float_t *pos) const;
  Int_t   GetSector(Float_t *pos) const;
  Int_t   GetPadX(Float_t *pos) const;
  Int_t   GetPadZ(Float_t *pos) const;
  void    GetPos(Int_t *det,Float_t *pos);
  Float_t GetX(Int_t *det) const;
  Float_t GetY(Int_t *det) const;
  Float_t GetZ(Int_t *det) const;
  Float_t GetPadDx(Float_t *pos);
  Float_t GetPadDy(Float_t *pos);
  Float_t GetPadDz(Float_t *pos);

  Int_t   NStripC()  const   { return kNStripC;};
  Int_t   NPadXSector() const { return (AliTOFGeometry::kNStripA + 2*AliTOFGeometry::kNStripB +
				  2*kNStripC)*AliTOFGeometry::kNpadX*AliTOFGeometry::kNpadZ;};

  Float_t RinTOF()  const    { return fgkxTOF;};
  Float_t Rmin()   const     { return fgkRmin;};
  Float_t Rmax()   const     { return fgkRmax;};

  Float_t ZlenA() const      { return fgkZlenA;};
  Float_t ZlenB() const      { return fgkZlenB;};
  Float_t ZlenC() const      { return fgkZlenC;};
  Float_t MaxhZtof() const   { return fgkMaxhZtof;};

  void Translation(Float_t *xyz, Float_t translationVector[3]) const;
  void Rotation(Float_t *xyz, Double_t rotationAngles[6]) const;
  void InverseRotation(Float_t *xyz, Double_t rotationAngles[6]) const;

  protected:

  enum {
    kNStripC    = 19 // number of strips in C type module 
  };

  //private:

  static const Float_t fgkZlenA;       // length (cm) of the A module
  static const Float_t fgkZlenB;       // length (cm) of the B module
  static const Float_t fgkZlenC;       // length (cm) of the C module
  static const Float_t fgkMaxhZtof;    // Max half z-size of TOF (cm)

  static const Float_t fgkRmin;        // Inner radius of the TOF (cm)
  static const Float_t fgkRmax;        // Outer radius of the TOF (cm)
  static const Float_t fgkxTOF;        // Inner TOF Radius used in Reconstruction (cm)

  ClassDef(AliTOFGeometryV5,1) // TOF Geometry class
};

#endif
