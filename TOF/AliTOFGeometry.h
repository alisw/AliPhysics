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



class AliTOFGeometry{

 public:
  AliTOFGeometry();
  virtual ~AliTOFGeometry();

  static  Int_t NStripA()     { return fgkNStripA;};
  static  Int_t NStripB()     { return fgkNStripB;};
  static  Int_t NStripC()     { return fgkNStripC;};
  static  Int_t NpadX()       { return fgkNpadX;};
  static  Int_t NpadZ()       { return fgkNpadZ;};
  static  Int_t NSectors()    { return fgkNSectors;};
  static  Int_t NPlates()     { return fgkNPlates;};
  static  Int_t NPadXSector() { return (fgkNStripA + 2*fgkNStripB +
				       2*fgkNStripC)*fgkNpadX*fgkNpadZ;};
  static  Int_t TimeDiff()    { return fgkTimeDiff;};
  static  Int_t MaxTOFTree()  { return fgkMaxTOFTree;};


  static  Float_t Rmin()     { return fgkRmin;};
  static  Float_t Rmax()     { return fgkRmax;};
  static  Float_t ZlenA()    { return fgkZlenA;};
  static  Float_t ZlenB()    { return fgkZlenB;};
  static  Float_t ZlenC()    { return fgkZlenC;};
  static  Float_t XPad()     { return fgkXPad;};
  static  Float_t ZPad()     { return fgkZPad;};
  static  Float_t MaxhZtof() { return fgkMaxhZtof;};


  static  Float_t SigmaForTail1() { return fgkSigmaForTail1;};
  static  Float_t SigmaForTail2() { return fgkSigmaForTail2;};
  static  Float_t SpeedOfLight()  { return fgkSpeedOfLight;};
  static  Float_t PionMass()      { return fgkPionMass;};
  static  Float_t KaonMass()      { return fgkKaonMass;};
  static  Float_t ProtonMass()    { return fgkProtonMass;};
  static  Float_t ElectronMass()  { return fgkElectronMass;};
  static  Float_t MuonMass()      { return fgkMuonMass;};
 
 

  virtual void    Init();
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
    fgkNStripA    = 15, // number of strips in A type module 
    fgkNStripB    = 19, // number of strips in B type module 
    fgkNStripC    = 20, // number of strips in C type module 
    fgkNpadX      = 48, // Number of pads along X 
    fgkNpadZ      = 2,  // Number of pads along Z
    fgkNSectors   = 18, // Number of Sectors
    fgkNPlates    = 5,  // Number of Plates
    fgkMaxNstrip  = 20, // Max. number of strips
    fgkMaxTOFTree = 5   // numer of geom. levels: 
  };

  static const Int_t fgkTimeDiff;  // Min signal separation (ps)

  static const Float_t fgkRmin;    // Inner radius of the TOF (cm)
  static const Float_t fgkRmax;    // Outer radius of the TOF (cm)
  static const Float_t fgkZlenA;   // length (cm) of the A module
  static const Float_t fgkZlenB;   // length (cm) of the B module
  static const Float_t fgkZlenC;   // length (cm) of the C module
  static const Float_t fgkXPad;    // Pad size in the x direction (cm)
  static const Float_t fgkZPad;    // Pad size in the z direction (cm)
  static const Float_t fgkMaxhZtof;// Max half z-size of TOF (cm)


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
 
  Float_t fAngles[fgkNPlates][fgkMaxNstrip]; //Strip Tilt Angles
  Float_t fHeights[fgkNPlates][fgkMaxNstrip];//Strip heights
  Float_t fPhiSec; //sector Phi width (deg)

  ClassDef(AliTOFGeometry,0) // TOF Geometry base class
};

#endif
