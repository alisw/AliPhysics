#ifndef ALIITSV11GEOMETRYSDD_H
#define ALIITSV11GEOMETRYSDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// class AliITSv11GeometrySDD
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

class TGeoVolume;
class TGeoTranslation;
class TGeoCombiTrans;
class TGeoArb8;
class TGeoNode;
class TGeoMedium;
class TGeoMatrix;
class AliITSgeom;
class AliITSv11GeomCableFlat;

#include "AliITSv11Geometry.h"


class AliITSv11GeometrySDD : public AliITSv11Geometry {

 public:
  AliITSv11GeometrySDD();
  AliITSv11GeometrySDD(Int_t debug);
  AliITSv11GeometrySDD(const AliITSv11GeometrySDD &source);
  AliITSv11GeometrySDD& operator=(const AliITSv11GeometrySDD &source);
  virtual ~AliITSv11GeometrySDD();

  // Main functions
  virtual void  Layer3(TGeoVolume *moth);
  virtual void  Layer4(TGeoVolume *moth);
  virtual Int_t ExportSensorGeometry(AliITSgeom *geom, Int_t iLaySDD,
				      Int_t startMod);
  virtual Int_t GetCurrentLayLaddDet(Int_t &lay, Int_t &ladd, Int_t&det) const;

  // Functions for coding, testing, debugging 
  void          AddHybrids(bool b)     {fAddHybrids    = b;};
  void          AddSensors(bool b)     {fAddSensors    = b;};
  void          AddHVcables(bool b)    {fAddHVcables   = b;};
  void          AddCables(bool b)      {fAddCables     = b;};
  void          AddCoolingSyst(bool b) {fAddCoolingSyst= b;};

  void          CheckOverlaps(Double_t precision = 0.01);
  void          AddOnlyLay3Ladder(Int_t min,Int_t max){
                  fAddOnlyLadder3min = min; fAddOnlyLadder3max = max; };
  void          AddOnlyLay4Ladder(Int_t min,Int_t max) {
                  fAddOnlyLadder4min = min; fAddOnlyLadder4max = max;};
  void          ShowOnePiece(TGeoVolume *Moth);

  virtual void  SetParameters();
  TGeoVolume*   GetMotherVolume() const { return fMotherVol;};
  const char*   GetSenstiveVolumeMame() const {return fgSDDsensitiveVolName;};
  Int_t         GetLay3NLadders() const;
  Int_t         GetLay4NLadders() const;

  private:

  // Create ladder virtual volumes and its detectors
  virtual TGeoVolume*  CreateLadder(Int_t iLay);
  virtual TGeoVolume*  CreateDetectors(Int_t iLay);
  // Create virtual volumes inside a ladder volume
  virtual TGeoVolume*  CreateLadderSegment(Int_t iLay, Int_t iSeg);
  virtual TGeoVolume*  CreateEndLadder(Int_t iLay);
  // Create some basic objects  
  virtual TGeoVolume*  CreateHybrid(Int_t iLRSide);
  virtual TGeoVolume*  CreatePinSupport();
  virtual TGeoVolume*  CreateCoolPipeSupportL();
  virtual TGeoVolume*  CreateCoolPipeSupportR();
  virtual TGeoVolume*  CreateSDDsensor();
  virtual TGeoVolume*  CreateBaseThermalBridge();
  void                 CreateBasicObjects();


  // Check that the nedium exists
  virtual TGeoMedium* GetMedium(const char* mediumName);

  // Create a TGeoCombiTrans: general rotation in phi and (dy,dz) translation 
  TGeoCombiTrans* CreateCombiTrans( const char *name,
				    Double_t dy, Double_t dz, Double_t dphi);

  // add (dx,dy,dz) translation to a initial TGeoCombiTrans
  void AddTranslationToCombiTrans( TGeoCombiTrans* ct,
			  Double_t dx=0, Double_t dy=0, Double_t dz=0) const;

  // Create one side of the CF corner of the CF structure
  TGeoArb8* CreateLadderSide( Double_t dz, Double_t angle, Double_t xSign,
			      Double_t L, Double_t H, Double_t l);

  //----------------------------
  TGeoVolume* fPinSupport;        //!  pins glued to sensors
  TGeoVolume* fCoolPipeSupportL;  //!  half of cooling pipe support
  TGeoVolume* fCoolPipeSupportR;  //!  half of cooling pipe support
  TGeoVolume* fSDDsensor;         //!  sensor and HV cables on it
  TGeoVolume* fBaseThermalBridge; //!  Base of hybrid thermal bridge
  TGeoVolume* fHybrid;            //!  hybrid volume

  static const Int_t fgkNladdSegCommonVol = 19;       //  Number of vol.
  TGeoVolume* fLaddSegCommonVol[fgkNladdSegCommonVol];//! volumes in ladder
  TGeoMatrix* fLaddSegCommonTr[fgkNladdSegCommonVol]; //! their transf.

  AliITSv11GeomCableFlat *fDigitCableLay3A; // layer 3 cables, side A
  AliITSv11GeomCableFlat *fDigitCableLay3B; // layer 3 cables, side A
  AliITSv11GeomCableFlat *fDigitCableLay4A; // layer 4 cables, side B
  AliITSv11GeomCableFlat *fDigitCableLay4B; // layer 4 cables, side B

  TGeoVolume *fMotherVol;    //! mother volume given in LayerX() funct.
  bool  fAddHybrids;         //  Insert hybrids ?   (default TRUE)
  bool  fAddSensors;         //  Insert sensors ?   (default TRUE)
  bool  fAddHVcables;        //  Insert HV cables ? (default TRUE)
  bool  fAddCables;          //  Insert cables ?    (default TRUE)
  bool  fAddCoolingSyst;     //  Insert cooling system ? (default TRUE)
  bool  fCoolingOn;          //  Insert cooling fluid ?  (default TRUE)
  Int_t fAddOnlyLadder3min;  //  first ladder index
  Int_t fAddOnlyLadder3max;  //  last  ladder index
  Int_t fAddOnlyLadder4min;  //  first ladder index
  Int_t fAddOnlyLadder4max;  //  last  ladder index
  Int_t fColorCarbonFiber;   //  display colors
  Int_t fColorRyton;         //  ===
  Int_t fColorPhynox;        //  ===
  Int_t fColorSilicon;       //  ===
  Int_t fColorAl;            //  ===
  Int_t fColorPolyhamide;    //  ===
  Int_t fColorGlass;         //  ===
  Int_t fColorSMD;           //  ===
  Int_t fColorSMDweld;       //  ===

  //--------------------------------------  parameters for the SDD geometry

  static const char* fgSDDsensitiveVolName;       // name of sensitive vol

  static const Int_t    fgkLay3Nladd;             // 14
  static const Int_t    fgkLay3Ndet;              //  6
  static const Double_t fgkLay3Rmin;              // min. radius of tube
  static const Double_t fgkLay3Rmax;              // max. radius of tube
  static const Double_t fgkLay3Length;            // length of layer 3 tube
  static const Double_t fgkLay3LadderLength;      // tot. length of ladder
  static const Double_t fgkLay3DetShortRadius;    // radius from beam axis
  static const Double_t fgkLay3DetLongRadius;     // radius from beam axis
  static const Double_t fgkLay3LaddTopCornerEnd;  // Ends of ladder 3
  static const Double_t fgkLay3ZPlusEndLength;    // ===

  static const Int_t    fgkLay4Nladd;             // 22
  static const Int_t    fgkLay4Ndet;              //  8
  static const Double_t fgkLay4Rmin;              // min. radius of tube
  static const Double_t fgkLay4Rmax;              // max. radius of tube
  static const Double_t fgkLay4Length;            // length of layer 4 tube
  static const Double_t fgkLay4LadderLength;      // tot. length of ladder
  static const Double_t fgkLay4DetShortRadius;    // radius from beam axis
  static const Double_t fgkLay4DetLongRadius;     // radius from beam axis
  static const Double_t fgkLay4LaddTopCornerEnd;  // Ends of ladder 3
  static const Double_t fgkLay4ZPlusEndLength;    // ===

  static const Double_t fgkSegmentLength;         // length of 1 ladder seg.
  static const Double_t fgkLadderWidth;           // carbon fiber structure 
  static const Double_t fgkLadderHeight;          // including bottom beam
  static const Double_t fgkLadderSegBoxDW;        // To include hybrids in box
  static const Double_t fgkLadderSegBoxDH;        // To include hybrids in box

  static const Double_t fgkLadderBeamRadius;      // carbon fiber beam radius
  static const Double_t fgkLadderLa;              // parameters defining
  static const Double_t fgkLadderHa;              //   the V side shape
  static const Double_t fgkLadderLb;              //   of the carbon
  static const Double_t fgkLadderHb;              //   fiber ladder
  static const Double_t fgkLadderl;               //   ============

  static const Double_t fgkBottomBeamAngle;       // bottom beam angle
  static const Double_t fgkBeamSidePhi;           // side beam angle

  static const Double_t fgkWaferThickness;        // sensor thickness (Y)
  static const Double_t fgkWaferWidth;            // width (X)
  static const Double_t fgkWaferLength;           // length (Z)
  static const Double_t fgkWaferThickSens;        // sensitive volume thich
  static const Double_t fgkWaferWidthSens;        // sens. volume width
  static const Double_t fgkWaferLengthSens;       // sens. volume length

  static const Double_t fgkSensorGlassLX;         // dimensions of glass
  static const Double_t fgkSensorGlassLZ;         //  (on which pins are
  static const Double_t fgkSensorGlassLY;         //   glued)
  static const Double_t fgkGlassDXOnSensor;       // Position of glass
  static const Double_t fgkGlassDZOnSensor;       //   on sensor

  static const Double_t fgkLadWaferSep;           // ladder-sensor dist.
  static const Double_t fgkPinR;                  // pins radius
  static const Double_t fgkPinSuppWidth;          // ===
  static const Double_t fgkPinSuppHeight;         // ===
  static const Double_t fgkPinSuppRmax;           // Parameters for pin
  static const Double_t fgkPinSuppLength;         //   supports on
  static const Double_t fgkPinSuppThickness;      //   carbon fiber
  static const Double_t fgkPinSuppConeAngle;      //   ladder
  static const Double_t fgkPinDXminOnSensor;      // ===
  static const Double_t fgkPinPinDDXOnSensor;     // ===
  static const Double_t fgkPinDYOnSensor;         // ===

  static const Double_t fgkCoolPipeInnerDiam;     // Water cooling
  static const Double_t fgkCoolPipeOuterDiam;     //   pipe
  static const Double_t fgkLay3CoolPipeSuppH;     // Heights of water
  static const Double_t fgkLay4CoolPipeSuppH;     //   pipes on ladders
  static const Double_t fgkCoolPipeSuppHeight;    // ===
  static const Double_t fgkCoolPipeSuppMaxLength; // ===
  static const Double_t fgkCoolPipeSuppWidthExt;  // Parameters for
  static const Double_t fgkCoolPipeSuppWidthIn;   //  cooling pipes
  static const Double_t fgkCoolPipeSuppHoleDiam;  //  on carbon fiber
  static const Double_t fgkCoolPipeSuppFulWidth;  //  ladder 
  static const Double_t fgkCoolPipeSuppTongW;     // ===
  static const Double_t fgkCoolPipeSuppAngle;     // ===
  static const Double_t fgkCoolPipeSuppSlitL;     // ===
  static const Double_t fgkCoolPipeSuppAxeDist;   // ===

  static const Double_t fgkBTBthick;              // BTB for :
  static const Double_t fgkBTBlength;             // Base of Thermal Bridge
  static const Double_t fgkBTBwidth;              // =====================
  static const Double_t fgkBTBaxisAtoBottom;      // axis A is the same as
  static const Double_t fgkBTBaxisAtoBase;        // the cooling pipe axis
  static const Double_t fgkRadiusAminBTB;         // ===
  static const Double_t fgkRadiusBminBTB;         // ===
  static const Double_t fgkBTBHoleLength;         // ===
  static const Double_t fgkBTBHolewidth;          // ===
  static const Double_t fgkBTBHoleRefX;           // ===
  static const Double_t fgkBTBHoleRefY;           // ===

  static const Double_t fgkHybridLength;          // Hybrid parameters :
  static const Double_t fgkHybridWidth;           // ===
  static const Double_t fgkHybridAngle;           // Hybrid on ladder in phi

  static const Double_t fgkHybRndHoleRad;         // ===
  static const Double_t fgkHybRndHoleZ;           // ===
  static const Double_t fgkHybRndHoleX;           // ===

  static const Double_t fgkHybFLlowHoleDZ;        // FLlow : low flex
  static const Double_t fgkHybFLlowHolePasDX;     // ===
  static const Double_t fgkHybFLlowHoleAmbDX;     // ===
  // (center of ships to the border)
  static const Double_t fgkHybFLlowChipZ4;        // Z1 to Z4 : position
  static const Double_t fgkHybFLlowChipZ3;        //   in z of the chip
  static const Double_t fgkHybFLlowChipZ2;        //   centers
  static const Double_t fgkHybFLlowChipZ1;        // ===
  static const Double_t fgkHybFLlowPasX;          // Pascal center X pos
  static const Double_t fgkHybFLlowAmbX;          // Ambra center X pos
  static const Double_t fgkHybChipsDZ;            // Z dimension of chips
  static const Double_t fgkHybPascalDX;           // X dimension of Pascal
  static const Double_t fgkHybAmbraDX;            // X dimension of Ambra
  static const Double_t fgkHybFLUpperWidth;       // bFLUpper : upper flex
  static const Double_t fgkHybFLUpperLength;      // ===
  static const Double_t fgkHybFLUpperAlDZ;        // ===
  static const Double_t fgkHybFLUpperAldx;        // ===

  static const Double_t fgkHybridThBridgeThick;   // Thicknesses :
  static const Double_t fgkHybAlThick;            // ===
  static const Double_t fgkHybUpThick;            // ===
  static const Double_t fgkHybGlueScrnThick;      // ===
  static const Double_t fgkHybGlueLowThick;       // ===
  static const Double_t fgkHybGlueUpThick;        // ===
  static const Double_t fgkHybAlCCThick;          // ===
  static const Double_t fgkHybUpCCThick;          // ===
  static const Double_t fgkHybChipThick;          // ===
  static const Double_t fgkHybGlueAgThick;        // ===
  static const Double_t fgkHybUnderNiThick;       // ===
  static const Int_t    fgkNHybSMD;               // Number of SMD
  static const Double_t fgkHybSMDposX[25];        // X pos. of SMD
  static const Double_t fgkHybSMDposZ[25];        // Z pos. of SMD
  static const Double_t fgkHybSMDmiddleW;         // SMD width
  static const Double_t fgkHybSMDmiddleL;         // SMD length
  static const Double_t fgkHybSMDendW;            // end SMD witdh
  static const Double_t fgkHybSMDendL;            // end SMD length
  static const Double_t fgkHybSMDheight;          // SMD height

  static const Double_t fgkDigitCablWidth;        // Digital
  static const Double_t fgkDigitCablAlThick;      // cables
  static const Double_t fgkDigitCablPolyThick;    // ===

  //HV cables
  static const Double_t fgkWaHVcableAlThick;      // Wrap-around
  static const Double_t fgkWaHVcablePolyThick;    //   High Voltage
  static const Double_t fgkWaHVcableLength;       //   cables
  static const Double_t fgkWaHVcableWitdh;        //   (on sensor)
  static const Double_t fgkWaHVcableDW;           // ===

  static const Double_t fgkTransitHVAlThick;      // Transition
  static const Double_t fgkTransitHVPolyThick;    //   High Voltage
  static const Double_t fgkTransitHVHeadLX;       //   cables
  static const Double_t fgkTransitHVHeadLZ;       //   (on sensor)
  static const Double_t fgkTransitHVBondingLZ;    // ===
  static const Double_t fgkTransitHVtailLength;   // ===
  static const Double_t fgkTransitHVtailWidth;    // ===
  static const Double_t fgkTransitHVtailXpos;     // ===
  static const Double_t fgkTransitHVsideLZ;       // ===
  static const Double_t fgkTransitHVsideLeftZ;    // ===
  static const Double_t fgkTransitHVsideRightZ;   // ===

  static const Double_t fgkLongHVcablePolyThick;  // Long High
  static const Double_t fgkLongHVcableAlThick;    //   Voltage
  static const Double_t fgkLongHVcableSeparation; //   cables

  static const Double_t fgkmu;  // 1 micron, or more for debugging

  // calculated parameters
  Double_t fLay3LadderUnderSegDH;  // To include HVcables in box
  Double_t fLay4LadderUnderSegDH;  // To include HVcables in box
  Double_t fLay3LaddShortRadius;   // ladder 3 to beam axis radius
  Double_t fLay3LaddLongRadius;    // ladder 3 to beam axis radius
  Double_t fLay4LaddShortRadius;   // ladder 4 to beam axis radius
  Double_t fLay4LaddLongRadius;    // ladder 4 to beam axis radius

  // parameters that be modified
  Double_t fLay3sensorZPos[6];     // Z pos of sensors in layer 3
  Double_t fLay4sensorZPos[8];     // Z pos of sensors in layer 4

  ClassDef(AliITSv11GeometrySDD,1) // ITS v11 SDD geometry
};


#endif
