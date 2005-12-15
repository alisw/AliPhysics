#ifndef ALITOFGEOMETRY_H
#define ALITOFGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TOF geometry class                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"

class AliTOFGeometry: public TObject{

 public:
  AliTOFGeometry();
  virtual ~AliTOFGeometry();

  static  Int_t NStripA()     { return kNStripA;};
  static  Int_t NStripB()     { return kNStripB;};
  virtual Int_t NStripC()     { return kNStripC;};
  virtual Int_t NMaxNstrip()  { return kMaxNstrip;};
  static  Int_t NpadX()       { return kNpadX;};
  static  Int_t NpadZ()       { return kNpadZ;};
  static  Int_t NpadXStrip()  { return kNpadX*kNpadZ;};
  static  Int_t NSectors()    { return kNSectors;};
  static  Int_t NPlates()     { return kNPlates;};
  virtual Int_t NPadXSector() { return (kNStripA + 2*kNStripB +
					2*kNStripC)*kNpadX*kNpadZ;};

  virtual Float_t RinTOF()    { return fgkxTOF;};
  virtual Float_t Rmin()      { return fgkRmin;};
  virtual Float_t Rmax()      { return fgkRmax;};

  static  Float_t XPad()      { return fgkXPad;};
  static  Float_t ZPad()      { return fgkZPad;};

  static  Int_t TimeDiff()    { return fgkTimeDiff;};
  static  Int_t MaxTOFTree()  { return kMaxTOFTree;};

  static  Int_t NDDL()        { return kNDDL;};
  static  Int_t NTRM()        { return kNTRM;}
  static  Int_t NTdc()        { return kNTdc;};
  static  Int_t NCh()         { return kNCh;};
  static  Int_t NPadXTRM()    { return kNCh*kNTdc;};

  virtual  Float_t ZlenA()       { return kZlenA;};
  virtual  Float_t ZlenB()       { return kZlenB;};
  virtual  Float_t ZlenC()       { return kZlenC;};
  virtual  Float_t MaxhZtof()    { return kMaxhZtof;};
  virtual  Float_t StripLength() { return kStripLength;};

  static  Float_t SigmaForTail1() { return fgkSigmaForTail1;};
  static  Float_t SigmaForTail2() { return fgkSigmaForTail2;};
 
  static  Double_t GetAlpha()  { return 2 * 3.14159265358979323846 / kNSectors; }; 
 
  static Float_t TdcBinWidth() {return fgkTdcBin;};


  virtual void    Init();
  virtual void    SetHoles(Bool_t holes) {fHoles = holes;};
  virtual Bool_t  GetHoles() const {return fHoles;};

  virtual Bool_t  IsInsideThePad(Int_t */*det*/, Float_t */*pos*/) {return kFALSE;};
  virtual Float_t DistanceToPad(Int_t */*det*/, Float_t */*pos*/, Float_t *dist3d=0) {return dist3d[0];};
  virtual void    GetPos(Int_t */*det*/,Float_t */*pos*/);
  virtual void    GetDetID(Float_t */*pos*/,Int_t */*det*/);
  virtual Int_t   GetPlate(Float_t */*pos*/) {return -1;};
  virtual Int_t   GetStrip(Float_t */*pos*/) {return -1;};
  virtual Int_t   GetSector(Float_t */*pos*/) {return -1;};
  virtual Int_t   GetPadX(Float_t */*pos*/) {return -1;};
  virtual Int_t   GetPadZ(Float_t */*pos*/) {return -1;};
  virtual Float_t GetX(Int_t */*det*/) {return -500.;};
  virtual Float_t GetY(Int_t */*det*/) {return -500.;};
  virtual Float_t GetZ(Int_t */*det*/) {return -500.;};

  Float_t GetAngles(Int_t iplate, Int_t istrip)  const {return fAngles[iplate][istrip];};
  Float_t GetHeights(Int_t iplate, Int_t istrip) const {return fHeights[iplate][istrip];};
  Float_t GetDistances(Int_t iplate, Int_t istrip) const {return fDistances[iplate][istrip];};

  //private:
  protected:

  enum {
    kNStripA    = 15, // number of strips in A type module 
    kNStripB    = 19, // number of strips in B type module 
    kNpadX      = 48, // Number of pads along X 
    kNpadZ      = 2,  // Number of pads along Z
    kNSectors   = 18, // Number of Sectors
    kNPlates    = 5,  // Number of Plates
    kMaxTOFTree = 5   // numer of geom. levels: 
  };

  // DAQ characteristics
  // cfr. TOF-TDR pag. 105 for Glossary
  // TARODA : TOF-ALICE Read Out and Data Acquisition system
  enum {
    kNDDL        =    4, // Number of DDL (Detector Data Link) per sector
    kNTRM        =   10, // Number of TRM ( Readout Module) per DDL
    kNTdc        =   30, // Number of Tdc (Time to Digital Converter) per TRM
    kNCh         =    8  // Number of channels per Tdc
  };

  static const Int_t fgkTimeDiff;      // Min signal separation (ps)

  mutable Int_t kNStripC;       // number of strips in C type module 
  mutable Int_t kMaxNstrip;     // Max. number of strips

  mutable Float_t kZlenA;       // length (cm) of the A module
  mutable Float_t kZlenB;       // length (cm) of the B module
  mutable Float_t kZlenC;       // length (cm) of the C module
  mutable Float_t kMaxhZtof;    // Max half z-size of TOF (cm)
  mutable Float_t kStripLength; // Strip Length (rho X phi direction) (cm)

  mutable Float_t fgkRmin;     // Inner radius of the TOF (cm)
  mutable Float_t fgkRmax;     // Outer radius of the TOF (cm)
  mutable Float_t fgkxTOF;     // Inner TOF Radius used in Reconstruction (cm)

  static const Float_t fgkXPad;     // Pad size in the x direction (cm)
  static const Float_t fgkZPad;     // Pad size in the z direction (cm)

  static const Float_t fgkSigmaForTail1;//Sig1 for simulation of TDC tails 
  static const Float_t fgkSigmaForTail2;//Sig2 for simulation of TDC tails

  Bool_t fHoles; //logical for geometry version (w/wo holes)

  Float_t *fAngles[kNPlates]; //Strip Tilt Angles
  Float_t *fHeights[kNPlates];//Strip heights
  Float_t *fDistances[kNPlates];//Strip distances

  Float_t fPhiSec; //sector Phi width (deg)

  static const Float_t fgkTdcBin;   // time-window for the TDC bins [ps]

  ClassDef(AliTOFGeometry,2) // TOF Geometry base class
};

#endif
