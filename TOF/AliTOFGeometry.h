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

#include "TGeoMatrix.h"
#include "TObject.h"

//class TGeoMatrix;

class AliTOFGeometry: public TObject{

 public:
  AliTOFGeometry();
  virtual ~AliTOFGeometry();

  static  Int_t NStripA()     { return kNStripA;};
  static  Int_t NStripB()     { return kNStripB;};
  static  Int_t NStripC()     { return kNStripC;};
  static  Int_t NStrip(Int_t nPlate);
  static  Int_t NMaxNstrip()  { return kMaxNstrip;};
  static  Int_t NpadX()       { return kNpadX;};
  static  Int_t NpadZ()       { return kNpadZ;};
  static  Int_t NpadXStrip()  { return kNpadX*kNpadZ;};
  static  Int_t NSectors()    { return kNSectors;};
  static  Int_t NPlates()     { return kNPlates;};
  static  Int_t NStripXSector() { return (kNStripA + 2*kNStripB +
						2*kNStripC);};
  static  Int_t NPadXSector() { return (kNStripA + 2*kNStripB +
					2*kNStripC)*kNpadX*kNpadZ;};

  static  Float_t RinTOF()  { return fgkxTOF;};
  static  Float_t Rmin()      { return fgkRmin;};
  static  Float_t Rmax()      { return fgkRmax;};

  static  Float_t XPad()     { return fgkXPad;};
  static  Float_t ZPad()     { return fgkZPad;};

  static  Float_t StripLength() { return fgkStripLength;};

  static  Float_t DeadTime()       { return fgkDeadTime;};
  static  Float_t MatchingWindow() { return fgkMatchingWindow;};

  static  Int_t MaxTOFTree()  { return kMaxTOFTree;};

  static  Int_t NDDL()        { return kNDDL;};
  static  Int_t NTRM()        { return kNTRM;}
  static  Int_t NTdc()        { return kNTdc;};
  static  Int_t NChain()      { return kNChain;};
  static  Int_t NCh()         { return kNCh;};
  static  Int_t NPadXTRM()    { return kNCh*kNTdc*kNChain;};

  static  Float_t ZlenA()       { return fgkZlenA;};
  static  Float_t ZlenB()       { return fgkZlenB;};
  static  Float_t ZlenC()       { return fgkZlenC;};
  static  Float_t MaxhZtof()    { return fgkMaxhZtof;};

  static  Float_t SigmaForTail1() { return fgkSigmaForTail1;};
  static  Float_t SigmaForTail2() { return fgkSigmaForTail2;};
 
  static  Double_t GetAlpha()  { return 2 * 3.14159265358979323846 / kNSectors; }; 
 
  static Float_t TdcBinWidth() {return fgkTdcBin;};
  static Float_t ToTBinWidth() {return fgkToTBin;};
  static Float_t BunchCrossingBinWidth() {return fgkBunchCrossingBin;};

  static Float_t SlewTOTMin() {return fgkSlewTOTMin;};
  static Float_t SlewTOTMax() {return fgkSlewTOTMax;};

  virtual void    ImportGeometry();
  virtual void    SetHoles(Bool_t holes) {fgHoles = holes;};
  static  Bool_t  GetHoles() {return fgHoles;};
  static  Float_t DistanceToPadPar(Int_t *det, const Float_t * pos, Float_t *dist3d=0);
  virtual Bool_t  IsInsideThePadPar(Int_t *det, const Float_t * pos) const;
  virtual Bool_t  IsInsideThePad(TGeoHMatrix *mat, const Float_t * pos, Float_t *dist3d=0) const;
  static  void    GetVolumePath(const Int_t * ind, Char_t *path );
  static  void    GetVolumePath(Int_t sector, Char_t *path );
  static  void    GetVolumePath(Int_t sector, Int_t plate, Int_t strip, Char_t *path );
  static  void    GetPos(Int_t *det,Float_t *pos);
  static  void    GetPosPar(Int_t *det,Float_t *pos);
  static  void    GetDetID(Float_t *pos,Int_t *det);
  static  Int_t   GetPlate(const Float_t * pos);
  static  Int_t   GetStrip(const Float_t * pos);
  static  Int_t   GetSector(const Float_t * pos);
  static  Int_t   GetPadX(const Float_t * pos);
  static  Int_t   GetPadZ(const Float_t * pos);
  static  Float_t GetX(const Int_t * det);
  static  Float_t GetY(const Int_t * det);
  static  Float_t GetZ(const Int_t * det);
  virtual void    DetToStripRF(Int_t nPadX, Int_t nPadZ,
			       Float_t &x,  Float_t &z) const;
  virtual void    DetToSectorRF(Int_t vol[5], Double_t coord[4][3]);
  static  Float_t GetPadDx(const Float_t * pos);
  static  Float_t GetPadDy(const Float_t * pos);
  static  Float_t GetPadDz(const Float_t * pos);
  static  void Translation(Float_t *xyz, Float_t translationVector[3]);
  static  void Rotation(Float_t *xyz, Double_t rotationAngles[6]);
  static  void InverseRotation(Float_t *xyz, Double_t rotationAngles[6]);

  static Float_t GetAngles(Int_t iplate, Int_t istrip)  {return fgkAngles[iplate][istrip];};
  static Float_t GetHeights(Int_t iplate, Int_t istrip)  {return fgkHeights[iplate][istrip];};
  static Float_t GetDistances(Int_t iplate, Int_t istrip)  {return fgkDistances[iplate][istrip];};

  static Int_t GetIndex(const Int_t * detId); // Get channel index from det Id (for calibration mainly)
  static void GetVolumeIndices(Int_t index, Int_t *detId); // Get volume index from channel index

  static UShort_t GetAliSensVolIndex(Int_t sec, Int_t pla, Int_t str); // Get the index of the TOF alignable volume in the AliGeomManager order
  static Int_t GetStripNumber(Int_t isector, Int_t iplate, Int_t istrip); // Get the serial number of the TOF alignable volume, i.e. the TOF strip
  static Int_t GetStripNumberPerSM(Int_t iplate, Int_t istrip); // Get the serial number of the TOF strip in a TOF SM
  static void GetStripAndModule(Int_t iStripPerSM, Int_t &iplate, Int_t &istrip); // Return the module and strip per module corresponding to the strip number per SM
  void PadRF2TrackingRF(Float_t *ctrackPos, Float_t *differenceT); // Convert the track coordinates from pad RF to tracking RF

  static Int_t GetTOFsupermodule(Int_t index); // Return the TOF supermodule where TOF channel index is located

  private:

  enum {
    kNStripA    = 15, // number of strips in A type module 
    kNStripB    = 19, // number of strips in B type module 
    kNStripC    = 19, // number of strips in C type module 
    kNpadX      = 48, // Number of pads along X 
    kNpadZ      = 2,  // Number of pads along Z
    kNSectors   = 18, // Number of Sectors
    kNPlates    = 5,  // Number of Plates
    kMaxTOFTree = 5,  // numer of geom. levels: 
    kMaxNstrip  = 20  // Max. number of strips
  };

  // DAQ characteristics
  // cfr. TOF-TDR pag. 105 for Glossary
  // TARODA : TOF-ALICE Read Out and Data Acquisition system
  enum {
    kNDDL        =    4, // Number of DDL (Detector Data Link) per sector
    kNTRM        =   12, // Number of TRM ( Readout Module) per DDL
    kNTdc        =   15, // Number of Tdc (Time to Digital Converter) per TRM
    kNChain      =    2, // Number of chains per TRM
    kNCh         =    8  // Number of channels per Tdc
  };

  static const Float_t fgkDeadTime;       // Single channel dead time (ps)
  static const Float_t fgkMatchingWindow; // Matching window (ps)

  static const Float_t fgkZlenA;       // length (cm) of the A module
  static const Float_t fgkZlenB;       // length (cm) of the B module
  static const Float_t fgkZlenC;       // length (cm) of the C module
  static const Float_t fgkMaxhZtof;    // Max half z-size of TOF (cm)

  static const Float_t fgkRmin;        // Inner radius of the TOF (cm)
  static const Float_t fgkRmax;        // Outer radius of the TOF (cm)
  static const Float_t fgkxTOF;        // Inner TOF Radius used in Reconstruction (cm)

  static const Float_t fgkStripLength; // Strip Length (rho X phi direction) (cm)

  static const Float_t fgkXPad;     // Pad size in the x direction (cm)
  static const Float_t fgkZPad;     // Pad size in the z direction (cm)

  static const Float_t fgkSigmaForTail1;//Sig1 for simulation of TDC tails 
  static const Float_t fgkSigmaForTail2;//Sig2 for simulation of TDC tails

  static const Float_t fgkPhiSec; //sector Phi width (deg)

  static Bool_t fgHoles; //logical for geometry version (w/wo holes)

  static const Float_t fgkAngles[kNPlates][kMaxNstrip];   //Strip Tilt Angles
  static const Float_t fgkHeights[kNPlates][kMaxNstrip];  //Strip heights
  static const Float_t fgkDistances[kNPlates][kMaxNstrip];//Strip distances


  static const Float_t fgkTdcBin;   // time-of-flight bin width [ps]
  static const Float_t fgkToTBin;   // time-over-threshold bin width [ps]
  static const Float_t fgkBunchCrossingBin; // bunch-crossing bin width [ps]

  static const Float_t fgkSlewTOTMin; // min TOT for slewing correction [ns]
  static const Float_t fgkSlewTOTMax; // max TOT for slewing correction [ns]

  ClassDef(AliTOFGeometry,9) // TOF Geometry base class
};

#endif
