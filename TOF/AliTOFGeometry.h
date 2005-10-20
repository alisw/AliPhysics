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
  static  Int_t NStripC()     { return kNStripC;};
  static  Int_t NpadX()       { return kNpadX;};
  static  Int_t NpadZ()       { return kNpadZ;};
  static  Int_t NpadXStrip()  { return kNpadX*kNpadZ;};
  static  Int_t NSectors()    { return kNSectors;};
  static  Int_t NPlates()     { return kNPlates;};
  static  Int_t NPadXSector() { return (kNStripA + 2*kNStripB +
				       2*kNStripC)*kNpadX*kNpadZ;};
  static  Int_t TimeDiff()    { return fgkTimeDiff;};
  static  Int_t MaxTOFTree()  { return kMaxTOFTree;};


  static  Int_t NDDL()        { return kNDDL;};
  static  Int_t NTRM()        { return kNTRM;}
  static  Int_t NTdc()        { return kNTdc;};
  static  Int_t NCh()         { return kNCh;};
  static  Int_t NPadXTRM()    { return kNCh*kNTdc;};


  static  Float_t RinTOF()      { return fgkxTOF;};
  static  Float_t Rmin()        { return fgkRmin;};
  static  Float_t Rmax()        { return fgkRmax;};
  static  Float_t ZlenA()       { return fgkZlenA;};
  static  Float_t ZlenB()       { return fgkZlenB;};
  static  Float_t ZlenC()       { return fgkZlenC;};
  static  Float_t XPad()        { return fgkXPad;};
  static  Float_t ZPad()        { return fgkZPad;};
  static  Float_t MaxhZtof()    { return fgkMaxhZtof;};
  static  Float_t StripLength() { return fgkStripLength;};
  static  Float_t DeadBndX()    { return fgkDeadBndX;};
  static  Float_t DeadBndZ()    { return fgkDeadBndZ;};
  static  Float_t OverSpc()     { return fgkOverSpc;};

  static  Float_t SigmaForTail1() { return fgkSigmaForTail1;};
  static  Float_t SigmaForTail2() { return fgkSigmaForTail2;};
  static  Float_t SpeedOfLight()  { return fgkSpeedOfLight;};
  static  Float_t PionMass()      { return fgkPionMass;};
  static  Float_t KaonMass()      { return fgkKaonMass;};
  static  Float_t ProtonMass()    { return fgkProtonMass;};
  static  Float_t ElectronMass()  { return fgkElectronMass;};
  static  Float_t MuonMass()      { return fgkMuonMass;};
 
  static  Double_t GetAlpha()  { return 2 * 3.14159265358979323846 / kNSectors; }; 
 
  static Float_t TdcBinWidth() {return fgkTdcBin;};


  virtual void    Init();
  virtual void    SetHoles(Bool_t holes) {fHoles = holes;};
  virtual Bool_t  GetHoles() const {return fHoles;};
  virtual Bool_t  IsInsideThePad(Int_t *det, Float_t *pos); 
  virtual Float_t DistanceToPad(Int_t *det, Float_t *pos, Float_t *dist3d=0); 
  virtual void    GetPos(Int_t *det,Float_t *pos);
  virtual void    GetDetID(Float_t *pos,Int_t *det);
  virtual Int_t   GetPlate(Float_t *pos);
  virtual Int_t   GetStrip(Float_t *pos);
  virtual Int_t   GetSector(Float_t *pos);
  virtual Int_t   GetPadX(Float_t *pos);
  virtual Int_t   GetPadZ(Float_t *pos);
  virtual Float_t GetX(Int_t *det);
  virtual Float_t GetY(Int_t *det);
  virtual Float_t GetZ(Int_t *det);
  virtual Float_t GetMinPlateTheta(Int_t iPlate);
  virtual Float_t GetMaxPlateTheta(Int_t iPlate);
  virtual Float_t GetMinStripTheta(Int_t iPlate, Int_t iStrip);
  virtual Float_t GetMaxStripTheta(Int_t iPlate, Int_t iStrip);
  virtual Float_t GetStripTheta(Int_t iPlate, Int_t iStrip);
  virtual Float_t GetAngles(Int_t iplate, Int_t istrip)  const {return fAngles[iplate][istrip];};
  virtual Float_t GetHeights(Int_t iplate, Int_t istrip) const {return fHeights[iplate][istrip];};

  private:

  enum {
    kNStripA    = 15, // number of strips in A type module 
    kNStripB    = 19, // number of strips in B type module 
    kNStripC    = 20, // number of strips in C type module 
    kNpadX      = 48, // Number of pads along X 
    kNpadZ      = 2,  // Number of pads along Z
    kNSectors   = 18, // Number of Sectors
    kNPlates    = 5,  // Number of Plates
    kMaxNstrip  = 20, // Max. number of strips
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

  static const Float_t fgkRmin;        // Inner radius of the TOF (cm)
  static const Float_t fgkRmax;        // Outer radius of the TOF (cm)
  static const Float_t fgkZlenA;       // length (cm) of the A module
  static const Float_t fgkZlenB;       // length (cm) of the B module
  static const Float_t fgkZlenC;       // length (cm) of the C module
  static const Float_t fgkXPad;        // Pad size in the x direction (cm)
  static const Float_t fgkZPad;        // Pad size in the z direction (cm)
  static const Float_t fgkMaxhZtof;    // Max half z-size of TOF (cm)
  static const Float_t fgkStripLength; // Strip Length (rho X phi direction) (cm)
  static const Float_t fgkDeadBndX;    // Dead Boundaries of a Strip along Z direction (width)
  static const Float_t fgkDeadBndZ;    // Dead Boundaries of a Strip along X direction (length)
  static const Float_t fgkOverSpc;     // Space available for sensitive layers in radial direction (cm)

  static const Float_t fgkxTOF;// Inner TOF Radius used in Reconstruction (cm)

  static const Float_t fgkSigmaForTail1;//Sig1 for simulation of TDC tails 
  static const Float_t fgkSigmaForTail2;//Sig2 for simulation of TDC tails
  static const Float_t fgkSpeedOfLight;// c (10^9 m/s)
  static const Float_t fgkPionMass;// pion mass (Gev/c^2)
  static const Float_t fgkKaonMass;// kaon mass (Gev/c^2)
  static const Float_t fgkProtonMass;// proton mass (Gev/c^2)
  static const Float_t fgkElectronMass;// electron mass (Gev/c^2)
  static const Float_t fgkMuonMass;// muon mass (Gev/c^2)


  static const Float_t fgkDprecMin;//num.prec.tolerance on Thmin 
  static const Float_t fgkDprecMax;//num.prec.tolerance on Thma 
  static const Float_t fgkDprecCen;//num.prec.tolerance on <Theta> 
  Bool_t fHoles; //logical for geometry version (w/wo holes) 
  Float_t fAngles[kNPlates][kMaxNstrip]; //Strip Tilt Angles
  Float_t fHeights[kNPlates][kMaxNstrip];//Strip heights
  Float_t fPhiSec; //sector Phi width (deg)

  static const Float_t fgkTdcBin;   // time-window for the TDC bins [ps]

  ClassDef(AliTOFGeometry,1) // TOF Geometry base class
};

#endif
