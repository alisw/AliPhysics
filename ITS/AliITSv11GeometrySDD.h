#ifndef ALIITSV11GEOMETRYSDD_H
#define ALIITSV11GEOMETRYSDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// class AliITSv11GeometrySDD
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

class TGeoVolume;
class TGeoVolumeAssembly;
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
  const char*   GetSenstiveVolumeName3() const {return fgSDDsensitiveVolName3;};
  const char*   GetSenstiveVolumeName4() const {return fgSDDsensitiveVolName4;};
  Int_t         GetLay3NLadders() const;
  Int_t         GetLay4NLadders() const;

  private:

  virtual TGeoVolumeAssembly*  CreateLadder(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateDetectorsAssembly(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateLadderSegment(Int_t iLay, Int_t iSeg);
  virtual TGeoVolumeAssembly*  CreateEndLadder(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateEndLadderCards(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateSupportRing(Int_t iLay);

  // Create some basic objects : 
  virtual void                 CreateSDDsensor();
  virtual TGeoVolume*          CreateHybrid(Int_t iLRSide);
  virtual TGeoVolume*          CreatePinSupport();
  virtual TGeoVolume*          CreateCoolPipeSupportL();
  virtual TGeoVolume*          CreateCoolPipeSupportR();
  virtual TGeoVolume*          CreateBaseThermalBridge();

  virtual TGeoVolumeAssembly*  CreateCarlosCard(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateLVCard(Int_t orientation);
  virtual TGeoVolumeAssembly*  CreateHVCard(Int_t iLay);

  void                         CreateBasicObjects();


  // Check that the nedium exists
  virtual TGeoMedium* GetMedium(const char* mediumName);

  // Create a TGeoCombiTrans: general rotation in phi and (dy,dz) translation 
  TGeoCombiTrans* CreateCombiTrans( const char *name,
				    Double_t dy, Double_t dz, Double_t dphi,
				    Bool_t planeSym=kFALSE);

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
  TGeoVolume* fSDDsensor3;        //!  sensor of lay. 3 and HV cables on it
  TGeoVolume* fSDDsensor4;        //!  sensor of lay. 4 and HV cables on it
  TGeoVolume* fBaseThermalBridge; //!  Base of hybrid thermal bridge
  TGeoVolume* fHybrid;            //!  hybrid volume

  static const Int_t fgkNladdSegCommonVol = 19;       //  Number of vol.
  TGeoVolume* fLaddSegCommonVol[fgkNladdSegCommonVol];//! volumes in ladder
  TGeoMatrix* fLaddSegCommonTr[fgkNladdSegCommonVol]; //! their transf.

  AliITSv11GeomCableFlat *fDigitCableLay3A; // layer 3 cables, side A
  AliITSv11GeomCableFlat *fDigitCableLay3B; // layer 3 cables, side B
  AliITSv11GeomCableFlat *fDigitCableLay4A; // layer 4 cables, side A
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
  Int_t fColorStesalite;     //  ===

  //--------------------------------------  parameters for the SDD geometry

  static const char* fgSDDsensitiveVolName3;       // sens. vol. name for lay. 3
  static const char* fgSDDsensitiveVolName4;       // sens. vol. name for lay. 4

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


  static const Double_t fgkRubyDX; // ruby dx with respect to the middle (to ladder z axis)
  static const Double_t fgkRubyZladd3;
  static const Double_t fgkRubyZladd4;

  static const Double_t fgkLadFoot_X; // Length of ladder foot
  static const Double_t fgkLadFoot_Z; // width 
  static const Double_t fgkLadFoot_Y; // thickness
  static const Double_t fgkLadFootMiddleY; // thickness in the middle part
  static const Double_t fgkLadBox1_X;
  static const Double_t fgkLadFingerPrint_X;
  static const Double_t fgkLadFingerPrint_Y;
  static const Double_t fgkLadFingerPrintBorder;
  static const Double_t fgkRubyCageHoleZ;
  static const Double_t fgkRubyCageHoleX;
  static const Double_t fgkRubyCageHoleY;
  static const Double_t fgkRubyCageAxisShift;
  static const Double_t fgkScrewM4diam;
  static const Double_t fgkRubyScrewShiftToCenterY;
  static const Double_t fgkRubyHoleDiam;

// the end ladder cooling pipe and its heat exchanger
  static const Double_t fgkEndLadPipeUlengthLay3;
  static const Double_t fgkEndLadPipeUlengthLay4;
  static const Double_t fgkEndLadPipeUwidth;
  static const Double_t fgkEndLadPipeRadius;
  static const Double_t fgkEndLadPipeInnerDiam;
  static const Double_t fgkEndLadPipeOuterDiam;

  static const Double_t fgkEndLadPipeArmZLay3;   //
  static const Double_t fgkEndLadPipeArmZLay4;   //
  static const Double_t fgkEndLadPipeArmX;    // the arms of the U cooling tube
  static const Double_t fgkEndLadPipeArmY;
  static const Double_t fgkEndLadPipeArmBoxDY; // shift in Y of the arms from the axis
  static const Double_t fgkEndLadPipeArmBoxDX;  // shift in X of the arms from the axis
  static const Double_t fgkEndLadPipeArmZpos; // 


  // approx dim for now - all of the following has to be checked
  // once Beppe provide the drawing

  // Carlos Card :
  static const Double_t fgkLVcard_X;
  static const Double_t fgkLVcard_Y;
  static const Double_t fgkLVcard_Z;
  static const Double_t fgkLVcard_CuZ;

  static const Double_t fgkLVChip0_X;
  static const Double_t fgkLVChip0_Y;
  static const Double_t fgkLVChip0_Z; // all except si layer
  static const Double_t fgkLVChip0_SiZ; //???????????????????????????????????????????????????
  static const Double_t fgkLVChip0_PosX;
  static const Double_t fgkLVChip0_PosY;

  static const Double_t fgkLVChip1_X;
  static const Double_t fgkLVChip1_Y;
  static const Double_t fgkLVChip1_Z;
  static const Double_t fgkLVChip1_SiZ;
  static const Double_t fgkLVChip1_PosX;
  static const Double_t fgkLVChip1_PosY;

  static const Double_t fgkLVChip2_X;
 static const Double_t fgkLVChip2_Y;
  static const Double_t fgkLVChip2_Z;
  static const Double_t fgkLVChip2_SiZ;
  static const Double_t fgkLVChip2_PosX;
  static const Double_t fgkLVChip2_PosY;

  static const Double_t fgkLVChip3_X;
  static const Double_t fgkLVChip3_Y;
  static const Double_t fgkLVChip3_Z;
  static const Double_t fgkLVChip3_SiZ;
  static const Double_t fgkLVChip3_PosX;
  static const Double_t fgkLVChip3_PosY;

  static const Double_t fgkLVcool_X1;
  static const Double_t fgkLVcool_Y1;
  static const Double_t fgkLVcool_Z1;

  static const Double_t fgkLVcool_X2;
  static const Double_t fgkLVcool_Y2;
  static const Double_t fgkLVcool_Z2;

  static const Double_t fgkLVcool_X3;
  static const Double_t fgkLVcool_Y3;
  static const Double_t fgkLVcoolPosY;

  // HV card :

  static const Double_t fgkHVCardCeramX;
  static const Double_t fgkHVCardCeramY;
  static const Double_t fgkHVCardCeramZ;

  static const Double_t fgkHVCardCapa1X;
  static const Double_t fgkHVCardCapa1Z;
  static const Double_t fgkHVCardCapa1Ymid;
  static const Double_t fgkHVCardCapa1Yend;
  static const Double_t fgkHVCardCapa1PosX;
  static const Double_t fgkHVCardCapa1PosY;

  static const Double_t fgkHVCardCapa2X;
  static const Double_t fgkHVCardCapa2Z;
  static const Double_t fgkHVCardCapa2Ymid;
  static const Double_t fgkHVCardCapa2Yend;
  static const Double_t fgkHVCardCapa2PosX;
  static const Double_t fgkHVCardCapa2PosY;

  static const Double_t fgkHVCardCapa3Xmid;
  static const Double_t fgkHVCardCapa3Xend;
  static const Double_t fgkHVCardCapa3Z;
  static const Double_t fgkHVCardCapa3Y;

  static const Double_t fgkHVCardCapa3PosX1;
  static const Double_t fgkHVCardCapa3PosX2;
  static const Double_t fgkHVCardCapa3PosX3;
  static const Double_t fgkHVCardCapa3PosX4;
  static const Double_t fgkHVCardCapa3PosX5;
  static const Double_t fgkHVCardCapa3PosY1;
  static const Double_t fgkHVCardCapa3PosY2;
  static const Double_t fgkHVCardCapa3PosY3;

  static const Double_t fgkHVCardCool1X;
  static const Double_t fgkHVCardCool1Y;
  static const Double_t fgkHVCardCool1Z;
  static const Double_t fgkHVCardCool2X;
  static const Double_t fgkHVCardCool2Y;
  static const Double_t fgkHVCardCool2Z;
  static const Double_t fgkHVCardCool3X;
  static const Double_t fgkHVCardCool3Y;
  static const Double_t fgkHVCardCool3Z;
  static const Double_t fgkHVCardCoolDY;

  static const Double_t fgkCarlosSuppX1;
  static const Double_t fgkCarlosSuppY1;
  static const Double_t fgkCarlosSuppX2;
  static const Double_t fgkCarlosSuppY2;
  static const Double_t fgkCarlosSuppZ;
  static const Double_t fgkCarlosSuppAngle;
  static const Double_t fgkCarlosSuppX3;
  static const Double_t fgkCarlosSuppY3;
  static const Double_t fgkCarlosSuppZ3;
  static const Double_t fgkCarlosSuppTopLen;

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

  ClassDef(AliITSv11GeometrySDD,2) // ITS v11 SDD geometry
};


#endif
