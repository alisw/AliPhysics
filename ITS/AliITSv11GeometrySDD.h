#ifndef ALIITSV11GEOMETRYSDD_H
#define ALIITSV11GEOMETRYSDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

//*************************************************************************
//
// SDD geometry, based on ROOT geometrical modeler
//
// Its integration to the aliroot framework is done in the AliITSv11Hybrid
// class (AliITSv11 not being functionnal so far)
//
// This geometry has no dependence with aliroot, you can run it with root
// only, provided that the AliITSv11GeomCable classes are also compiled
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************


class TGeoVolume;
class TGeoPcon;
class TGeoVolumeAssembly;
class TGeoTranslation;
class TGeoCombiTrans;
class TGeoArb8;
class TGeoNode;
class TGeoMedium;
class TGeoMatrix;
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
  virtual void  ForwardLayer3(TGeoVolume *moth);
  virtual void  ForwardLayer4(TGeoVolume *moth);
  virtual void  SDDCables(TGeoVolume *moth);

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
  virtual TGeoVolumeAssembly*  CreateDetectorsAssemblyLadd2();
  virtual TGeoVolume*          CreateLadderSegment(Int_t iLay, Int_t iSeg);
  virtual TGeoVolumeAssembly*  CreateEndLadder(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateEndLadderCards(Int_t iLay);
  virtual TGeoVolume*          CreateEndLadderCardsV(Int_t iLay);
  virtual TGeoVolumeAssembly*  CreateSupportRing(Int_t iLay);

  // Create some basic objects : 
  virtual void                 CreateSDDsensor();
  virtual TGeoVolume*          CreateHybrid(Int_t iLRSide);
  virtual TGeoVolume*          CreatePinSupport();
  virtual TGeoVolume*          CreateCoolPipeSupportL();
  virtual TGeoVolume*          CreateCoolPipeSupportR();
  virtual TGeoVolume*          CreateBaseThermalBridge();

  virtual TGeoVolumeAssembly*  CreateLadderFoot();
  virtual TGeoVolumeAssembly*  CreateCarlosCard(Int_t iLay);
  virtual Int_t                CreateLVCard();
  virtual TGeoVolumeAssembly*  CreateHVCard(Int_t iLay);

  void                         CreateBasicObjects();
  Double_t                     GetConeZ(Double_t r, Double_t refR1, Double_t refR2,
					Double_t refZ1, Double_t refZ2) const;
  TGeoPcon*                    CreateConeConstSection(Double_t r1max, Double_t z1,
						      Double_t r2max, Double_t z2,
						      Double_t section, Int_t nDiv=1);
  Int_t CreateAndInsetConeCablePart(TGeoVolume *mother, Double_t angle,
				    Int_t nLay3, Int_t nLay4,
				    Double_t r1, Double_t z1,
				    Double_t r2, Double_t z2);

  // Check that the medium exists
  virtual TGeoMedium* GetMedium(const char* mediumName);

  // Create a TGeoCombiTrans: general rotation in phi and (dy,dz) translation 
  TGeoCombiTrans* CreateCombiTrans( const char *name,
				    Double_t dy, Double_t dz, Double_t dphi,
				    Bool_t planeSym=kFALSE);

  // add (dx,dy,dz) translation to a initial TGeoCombiTrans
  void AddTranslationToCombiTrans( TGeoCombiTrans* ct,
			  Double_t dx=0, Double_t dy=0, Double_t dz=0) const;

  // Create one side of the CF corner of the CF structure
  TGeoArb8* CreateLadderSide( const char *name,
			      Double_t dz, Double_t angle, Double_t xSign,
			      Double_t L, Double_t H, Double_t l);

  //----------------------------
  TGeoVolume* fPinSupport;        //!  pins glued to sensors
  TGeoVolume* fCoolPipeSupportL;  //!  half of cooling pipe support
  TGeoVolume* fCoolPipeSupportR;  //!  half of cooling pipe support
  TGeoVolume* fSDDsensor3;        //!  sensor of lay. 3 and HV cables on it
  TGeoVolume* fSDDsensor4;        //!  sensor of lay. 4 and HV cables on it
  TGeoVolume* fBaseThermalBridge; //!  Base of hybrid thermal bridge
  TGeoVolume* fHybrid;            //!  hybrid volume
  TGeoVolumeAssembly *fLadderFoot;//!  ladder foot in stesalite
  TGeoVolumeAssembly *fCardLVR;   //!  low voltage card, right
  TGeoVolumeAssembly *fCardLVL;   //!  low voltage card, left
  TGeoVolumeAssembly *fCardHV;    //!  high voltage card
  TGeoVolumeAssembly *fCardCarlos;//!  end-ladder CARLOS card
  TGeoVolumeAssembly *fRaccordoL; //!  link between cooling tubes at end ladder
  TGeoVolume* fCommonVol[2];      //!  some common vol. used in several places
  TGeoMatrix* fCommonTr[2];       //!  some common transformations

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

  static const char* fgSDDsensitiveVolName3;      // sens. vol. name for lay. 3
  static const char* fgSDDsensitiveVolName4;      // sens. vol. name for lay. 4

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

  static const Double_t fgkEndLaddCardsShortRadiusLay3; // ref radius of end ladder cards
  static const Double_t fgkEndLaddCardsShortRadiusLay4; // ref radius  of end ladder cards
  static const Double_t fgkDistEndLaddCardsLadd;        // dist. between U cooling tube and ladder

  static const Double_t fgkSegmentLength;         // length of 1 ladder seg.
  static const Double_t fgkLadderWidth;           // carbon fiber structure 
  static const Double_t fgkLadderHeight;          // including bottom beam
  static const Double_t fgkLadderSegBoxDW;        // To include hybrids in box
  static const Double_t fgkLadderSegBoxDH;        // To include hybrids in box
  //  static const Double_t fgkLadderSegBoxDHCorr;    // To include hybrids in box

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
  static const Double_t fgkHybCC2SensorLen;       // ===
  static const Double_t fgkHybCC2SensorWid;       // ===
  static const Double_t fgkHybCC2SensorAng;       // ===
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


  static const Double_t fgkRubyDX;                // ruby dx with respect to the middle (to ladder z axis)
  static const Double_t fgkRubyZladd3;            // Z of ruby, ladder 3
  static const Double_t fgkRubyZladd4;            // Z of ruby, ladder 4

  static const Double_t fgkLadFootX;              // Length of ladder foot
  static const Double_t fgkLadFootZ;              // width 
  static const Double_t fgkLadFootY;              // thickness
  static const Double_t fgkLadFootMiddleY;        // thickness in the middle part
  static const Double_t fgkLadBox1X;              // size in X
  static const Double_t fgkLadFingerPrintX;       // size in X
  static const Double_t fgkLadFingerPrintY ;      // size in Y
  static const Double_t fgkLadFingerPrintBorder;  // size in X
  static const Double_t fgkRubyCageHoleZ;         // size in Z
  static const Double_t fgkRubyCageHoleX;         // size in X
  static const Double_t fgkRubyCageHoleY;         // size in Y
  static const Double_t fgkRubyCageAxisShift;     // shift in X
  static const Double_t fgkScrewM4diam;           // M4 screw standard diameter 
  static const Double_t fgkRubyScrewShiftToCenterY; // screw placement 
  static const Double_t fgkRubyHoleDiam;          // guess what

// the end ladder cooling pipe and its heat exchanger
  static const Double_t fgkEndLadPipeUlengthLay3; // length in Z of the U cooling tube
  static const Double_t fgkEndLadPipeUlengthLay4; // length in Z of the U cooling tube
  static const Double_t fgkEndLadPipeUwidth;      // width
  static const Double_t fgkEndLadPipeRadius;      // radius
  static const Double_t fgkEndLadPipeInnerDiam;   // InnerDiam
  static const Double_t fgkEndLadPipeOuterDiam;   // OuterDiam

  static const Double_t fgkEndLadPipeArmZLay3;    // the arms of the U cooling tube
  static const Double_t fgkEndLadPipeArmZLay4;    // (rectangular part surrounding the)
  static const Double_t fgkEndLadPipeArmX;        // the tube
  static const Double_t fgkEndLadPipeArmY;        // X, Y : size in the correxponding axis
  static const Double_t fgkEndLadPipeArmBoxDY;    // shift in Y of the arms from the axis
  static const Double_t fgkEndLadPipeArmBoxDX;    // shift in X of the arms from the axis
  static const Double_t fgkEndLadPipeArmZpos;     // position with respect to tube  


  // approx dim for now - all of the following has to be checked
  // once Beppe provide the drawing...

  // Carlos Card :
  static const Double_t fgkLVcardX;              // size of the card itself in X
  static const Double_t fgkLVcardY;              // size of the card itself in Y
  static const Double_t fgkLVcardZ;              // size of the card itself in Z
  static const Double_t fgkLVcardCuZ;            // Cu thickness

  static const Double_t fgkLVChip0X;             // chip #0
  static const Double_t fgkLVChip0Y;             // ...
  static const Double_t fgkLVChip0Z;             //  thickness without si layer
  static const Double_t fgkLVChip0SiZ;           // Si layer thickness
  static const Double_t fgkLVChip0PosX;          // Position with respect to the card
  static const Double_t fgkLVChip0PosY;          // Position with respect to the card

  static const Double_t fgkLVChip1X;             // same
  static const Double_t fgkLVChip1Y;             // conventions
  static const Double_t fgkLVChip1Z;             // as
  static const Double_t fgkLVChip1SiZ;           // chip 0
  static const Double_t fgkLVChip1PosX;          // ==
  static const Double_t fgkLVChip1PosY;          // ==

  static const Double_t fgkLVChip2X;             // same
  static const Double_t fgkLVChip2Y;             // conventions
  static const Double_t fgkLVChip2Z;             // as
  static const Double_t fgkLVChip2SiZ;           // chip 0
  static const Double_t fgkLVChip2PosX;          // ==
  static const Double_t fgkLVChip2PosY;          // ==

  static const Double_t fgkLVChip3X;             // same
  static const Double_t fgkLVChip3Y;             // conventions
  static const Double_t fgkLVChip3Z;             // as
  static const Double_t fgkLVChip3SiZ;           // chip 0
  static const Double_t fgkLVChip3PosX;          // ==
  static const Double_t fgkLVChip3PosY;          // ==

  static const Double_t fgkLVcoolX1;             // pieces of alCu12
  static const Double_t fgkLVcoolY1;             // for heat exchange
  static const Double_t fgkLVcoolZ1;             // with the cooling tube

  static const Double_t fgkLVcoolX2;             // X,Y,Z are
  static const Double_t fgkLVcoolY2;             // dimensions
  static const Double_t fgkLVcoolZ2;             // of the pieces

  static const Double_t fgkLVcoolX3;             // ==
  static const Double_t fgkLVcoolY3;             // ==
  static const Double_t fgkLVcoolPosY;           // ==

  // HV card :
  static const Double_t fgkHVCardCeramX;         // size in X of the ceramic card
  static const Double_t fgkHVCardCeramY;         // size in Y
  static const Double_t fgkHVCardCeramZ;         // size in Z

  static const Double_t fgkHVCardCapa1X;         // size in X of the capa 1
  static const Double_t fgkHVCardCapa1Z;         // size in Z
  static const Double_t fgkHVCardCapa1Ymid;      // size of the middle part
  static const Double_t fgkHVCardCapa1Yend;      // ...
  static const Double_t fgkHVCardCapa1PosX;      // position on the card
  static const Double_t fgkHVCardCapa1PosY;      // position on the card

  static const Double_t fgkHVCardCapa2X;         // idem for second type capa
  static const Double_t fgkHVCardCapa2Z;         //   love me
  static const Double_t fgkHVCardCapa2Ymid;      //   ...
  static const Double_t fgkHVCardCapa2Yend;      //   tender,
  static const Double_t fgkHVCardCapa2PosX;      //   ...
  static const Double_t fgkHVCardCapa2PosY;      //   love me true

  static const Double_t fgkHVCardCapa3Xmid;      //  idem for third type capa
  static const Double_t fgkHVCardCapa3Xend;      //  ===
  static const Double_t fgkHVCardCapa3Z;         //  ===
  static const Double_t fgkHVCardCapa3Y;         //  ===

  static const Double_t fgkHVCardCapa3PosX1;     // this capa is placed
  static const Double_t fgkHVCardCapa3PosX2;     // in several positions
  static const Double_t fgkHVCardCapa3PosX3;     // ...
  static const Double_t fgkHVCardCapa3PosX4;     // ===
  static const Double_t fgkHVCardCapa3PosX5;     // ===
  static const Double_t fgkHVCardCapa3PosY1;     // ===
  static const Double_t fgkHVCardCapa3PosY2;     // ===
  static const Double_t fgkHVCardCapa3PosY3;     // ===

  static const Double_t fgkHVCardCool1X;         // cooling
  static const Double_t fgkHVCardCool1Y;         // pieces for
  static const Double_t fgkHVCardCool1Z;         // heat exchange
  static const Double_t fgkHVCardCool2X;         // with
  static const Double_t fgkHVCardCool2Y;         // cooling U tube
  static const Double_t fgkHVCardCool2Z;         // ===
  static const Double_t fgkHVCardCool3X;         // ===
  static const Double_t fgkHVCardCool3Y;         // ===
  static const Double_t fgkHVCardCool3Z;         // ===
  static const Double_t fgkHVCardCoolDY;         // ===

  static const Double_t fgkCarlosSuppX1;         // piece with which
  static const Double_t fgkCarlosSuppY1;         // the carlos card
  static const Double_t fgkCarlosSuppX2;         // is fixed
  static const Double_t fgkCarlosSuppY2;         // ===
  static const Double_t fgkCarlosSuppZ;          // ===
  static const Double_t fgkCarlosSuppAngle;      //  ===
  static const Double_t fgkCarlosSuppX3;         // ===
  static const Double_t fgkCarlosSuppY3;         // ===
  static const Double_t fgkCarlosSuppZ3;         // ===
  static const Double_t fgkCarlosSuppTopLen;     // ===

  // screws fixing the board on the U tube
  static const Double_t fgkLittleScrewHeadR;     // screws fixing boards
  static const Double_t fgkLittleScrewHeadH;     // Value to be checked
  static const Double_t fgkLittleScrewR;         // ===
  static const Double_t fgkShiftLittleScrewLV;   // ===
  static const Double_t fgkLittleLVScrewHeadR;   // ===

  // CARLOS board
  static const Double_t fgkCarlosCardX1;         // length (first part of Carlos card)
  static const Double_t fgkCarlosCardY1;         // thickness
  static const Double_t fgkCarlosCardZ1;         // width 
  static const Double_t fgkCarlosCardCuY;        // thickness of Cu layer (strips)
  static const Double_t fgkCarlosCardX2;         // length (2nd part of Carlos card)
  static const Double_t fgkCarlosCardZ2;         // width 

  static const Double_t fgkCarlosCardChipSiThick;  // Carlos Chip thicknes - value to be checked
  static const Double_t fgkCarlosCardShift;        // (value to be checked) shift in z w.r.t. heat bridge 

  // size and position of various chips on carlos end-ladder board
  static const Double_t fgkCarlosU1X;            // chip size in X
  static const Double_t fgkCarlosU1Y;            // chip size in Y
  static const Double_t fgkCarlosU1Z;            // chip size in Z
  static const Double_t fgkCarlosU1posX;         // position in X
  static const Double_t fgkCarlosU1posZ;         // position in Z

  static const Double_t fgkCarlosU2X;            // chip size in X
  static const Double_t fgkCarlosU2Y;            // chip size in Y
  static const Double_t fgkCarlosU2Z;            // chip size in Z
  static const Double_t fgkCarlosU2posX;         // position in X
  static const Double_t fgkCarlosU2posZ;         // position in Z
  
  static const Double_t fgkCarlosU3X;            // same convention
  static const Double_t fgkCarlosU3Y;            // ===
  static const Double_t fgkCarlosU3Z;            // ===
  static const Double_t fgkCarlosU3posX;         // ===
  static const Double_t fgkCarlosU3posZ;         // ===

  // U4 like U3  
  static const Double_t fgkCarlosU4posX;         // same convention
  static const Double_t fgkCarlosU4posZ;         // ===

  static const Double_t fgkCarlosU17X;           // same convention
  static const Double_t fgkCarlosU17Y;           // ===
  static const Double_t fgkCarlosU17Z;           // ===
  static const Double_t fgkCarlosU17posX;        // ===
  static const Double_t fgkCarlosU17posZ;        // ===
  
  static const Double_t fgkCarlosU35X;           // same convention
  static const Double_t fgkCarlosU35Y;           // ===
  static const Double_t fgkCarlosU35Z;           // ===
  static const Double_t fgkCarlosU35posX;        // ===
  static const Double_t fgkCarlosU35posZ;        // ===

  static const Double_t fgkCarlosU36X;           // same convention
  static const Double_t fgkCarlosU36Y;           // ===
  static const Double_t fgkCarlosU36Z;           // ===
  static const Double_t fgkCarlosU36posX;        // ===
  static const Double_t fgkCarlosU36posZ;        // ===
  
  static const Double_t fgkCarlosQZ1X;           // same convention
  static const Double_t fgkCarlosQZ1Y;           // look more thick than design number (0.7) ! to be checked
  static const Double_t fgkCarlosQZ1Z;           // to be checked
  static const Double_t fgkCarlosQZ1posX;        // to be checked
  static const Double_t fgkCarlosQZ1posZ;        // to be checked

  // some pieces at the end of the carbon fiber ladder
  static const Double_t fgkCoolPipeLay3Len;  // value to be checked
  static const Double_t fgkCoolPipeLay4Len;  // ===
  static const Double_t fgkHVguideX1;    // ===
  static const Double_t fgkHVguideY1;    // ===
  static const Double_t fgkHVguideZ1;    // ===
  static const Double_t fgkHVguideZ2;    // ===
  static const Double_t fgkHVguideDX;    // ===
  static const Double_t fgkHVguideSuppFullZ;    // ===
  
  // Cooling connector between phynox and plastic cooling water tubes
  static const Double_t fgkConnectorCoolTubeRmin; // internal radius
  static const Double_t fgkConnectorCoolTubeR1; // value to be checked
  static const Double_t fgkConnectorCoolTubeL1;  // ===
  static const Double_t fgkConnectorCoolTubeR2;  // ===
  static const Double_t fgkConnectorCoolTubeL2;  // ===
  static const Double_t fgkConnectorCoolTubeR3;  // ===
  static const Double_t fgkConnectorCoolTubeL3;  // ===

  // parameters for coding SDD cables on SDD and SSD cones
  static const Double_t fgkSectionCuPerMod;    // area of copper per mod.
  static const Double_t fgkSectionPlastPerMod; // area of plast per mod.
  static const Double_t fgkSectionGlassPerMod; // area of optical fiber per mod.
  static const Double_t fgkSectionCoolPolyuEL; // area of cooling tubes on End Ladders
  static const Double_t fgkSectionCoolWaterEL; // area of cooling water on End Ladders
  static const Double_t fgkEndLadderEarthCableR; // radius of the earth cable on End Ladders
  static const Double_t fgkCableBendRatio; // ??? this factor account for the bending of cables
  static const Double_t fgkHybridAlFoilThick; // Thickness of Al foil on hybrid side
  static const Double_t fgkHybridAlFoilWide; // Width of Al foil on hybrid side
  static const Double_t fgkHybridAlFoilSide; // Side length of Al foil on hybrid side

  static const Double_t fgkConeSDDr1; // define SDD cone slope and pos
  static const Double_t fgkConeSDDr2; // define SDD cone slope and pos
  static const Double_t fgkConeSDDz1; // define SDD cone slope and pos
  static const Double_t fgkConeSDDz2; // define SDD cone slope and pos

  static const Double_t fgkSDDCableR1; // ??? // part 1 of "cable cone"
  static const Double_t fgkSDDCableR2; // ??? // part 1/2 of "cable cone"
  static const Double_t fgkSDDCableR3; // ??? // part 2 of "cable cone"

  static const Double_t fgkSDDCableDZint; // length of intermediate cylinder
  static const Double_t fgkSDDCableR5; // third part of "cable cone"
  static const Double_t fgkSDDCableZ5; // third part of "cable cone"



  // distance from the heat bridge center to the card center :
  static const Double_t fgkCarlosCard2HeatBridge;// distance from the heat bridge center to the card center

  static const Double_t fgkmu;  // 1 micron, or more for debugging

  // calculated parameters
  Double_t fLay3LadderUnderSegDH;  // To include HVcables in box
  Double_t fLay4LadderUnderSegDH;  // To include HVcables in box
  Double_t fLay3LaddShortRadius;   // ladder 3 to beam axis radius
  Double_t fLay3LaddLongRadius;    // ladder 3 to beam axis radius
  Double_t fLay4LaddShortRadius;   // ladder 4 to beam axis radius
  Double_t fLay4LaddLongRadius;    // ladder 4 to beam axis radius

  // parameters that can be modified
  Double_t fLay3sensorZPos[6];     // Z pos of sensors in layer 3
  Double_t fLay4sensorZPos[8];     // Z pos of sensors in layer 4

  ClassDef(AliITSv11GeometrySDD,0) // ITS v11 SDD geometry
};


#endif
