#ifndef ALIITSV11GEOMETRYSDD_H
#define ALIITSV11GEOMETRYSDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Ludovic Gaudichet (Ludovic.Gaudichet@to.infn.it)

class TGeoVolume;
class TGeoCombiTrans;
class TGeoArb8;
class TGeoNode;

#include <TObjArray.h>

//----------------------------------------------------------------------
class AliITSv11GeomSDDcable: public TObject {
 public:
  AliITSv11GeomSDDcable():TObject(){fInitialNode = 0; fPointArray.SetOwner(); };
  virtual ~AliITSv11GeomSDDcable();
  void AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt, Double_t *coord);
  Int_t GetNCheckPoints() {return fVolumeArray.GetEntriesFast();};
  Int_t GetCheckPoint( Int_t iCheckPt, Int_t nOccur, Int_t motherLevel,
		       Double_t *coord);
  void SetInitialNode(TGeoVolume *vol);
  void ResetInitialNode();

 protected:
  bool CheckDaughter(Int_t i, TGeoNode* node);

  Int_t fNodeInd[50];
  TObjArray fVolumeArray;
  TObjArray fPointArray;
  TGeoVolume *fCurrentVol;
  TGeoNode *fInitialNode;

  ClassDef(AliITSv11GeomSDDcable,1) // ITS v11 SDD cable geometry
};



//----------------------------------------------------------------------
class AliITSv11GeomSDDcableNap : public AliITSv11GeomSDDcable {
 public:
  AliITSv11GeomSDDcableNap(Double_t width, Double_t thick);
  virtual ~AliITSv11GeomSDDcableNap() {};
  Int_t CreateAndInsertCableSegment(Int_t p1, Int_t p2, TGeoVolume *motherVol);
 protected:
  Double_t fWidth;
  Double_t fThick;

  ClassDef(AliITSv11GeomSDDcableNap,1) // ITS v11 SDD cable Nap geometry
};



//----------------------------------------------------------------------
class AliITSv11GeometrySDD : public AliITSv11Geometry {

 public:
  AliITSv11GeometrySDD();
  AliITSv11GeometrySDD(Int_t debug);
  virtual ~AliITSv11GeometrySDD(){};

  // Main functions
  virtual void Layer3(TGeoVolume *Moth);
  virtual void Layer4(TGeoVolume *Moth);

  // Functions for coding, testing, debugging 
  void ShowOnePiece(TGeoVolume *Moth);
  void CheckOverlaps(Double_t precision = 0.01);
  void AddOnlySegment(Int_t i) { fAddOnlySegment=i; };
  void AddOnlyLay3Ladder(Int_t min,Int_t max){
    fAddOnlyLadder3min = min; fAddOnlyLadder3max = max; };
  void AddOnlyLay4Ladder(Int_t min,Int_t max) {
    fAddOnlyLadder4min = min; fAddOnlyLadder4max = max;};

  virtual void SetGeomParameters();


  AliITSv11GeomSDDcable cable;

  private:

  // Create ladder virtual volumes and its detectors (inside layers 3 and 4)
  virtual TGeoVolume* CreateLay3Ladder();
  virtual TGeoVolume* CreateLay3Detectors();
  virtual TGeoVolume* CreateLay4Ladder();
  virtual TGeoVolume* CreateLay4Detectors();

  // Create virtual volumes inside a ladder volume
  virtual TGeoVolume* CreateLadderSegment(Int_t iLay, Int_t iSeg);
  virtual TGeoVolume* CreateEndLadder(Int_t iLay, Int_t);
  virtual TGeoVolume* CreateHybrid(Int_t iSeg);


  // Create a TGeoCombiTrans : general rotation in phi and (dy,dz) translation 
  TGeoCombiTrans* CreateCombiTrans( const char *name,
				    Double_t dy, Double_t dz, Double_t dphi);

  // add (dx,dy,dz) translation to a initial TGeoCombiTrans
  void AddTranslationToCombiTrans( TGeoCombiTrans* ct,
				    Double_t dx=0, Double_t dy=0, Double_t dz=0);

  // Create one side of the CF corner of the CF structure
  TGeoArb8* CreateLadderSide( Double_t dz, Double_t angle, Double_t xSign,
			      Double_t L, Double_t H, Double_t l);

  // Create de CF structure in the ladder segment
  void AddLadderCFstruct(Double_t dy, TGeoVolume* vol);

  // Create a pin support object
  virtual TGeoVolume* CreatePinSupport();

  virtual TGeoVolume* CreateCoolPipeSupportL();
  virtual TGeoVolume* CreateCoolPipeSupportR();

  // Create a SDD detector object
  virtual TGeoVolume* CreateSDDsensor();

  // Create the base of the thermal bridge, supporting hybrids on ladder
  virtual TGeoVolume* CreateBaseThermalBridge();

  //**************************************************

  Int_t fAddOnlyLadder3min;
  Int_t fAddOnlyLadder3max;
  Int_t fAddOnlyLadder4min;
  Int_t fAddOnlyLadder4max;
  Int_t fAddOnlySegment;
  Int_t fColorCarbonFiber;
  Int_t fColorRyton;
  Int_t fColorPhynox;
  Int_t fColorSilicon;

  bool fCoolingOn;

  // parameters of the SDD geometry :
  static const Double_t fLay3Rmin;
  static const Double_t fLay3Rmax;
  static const Double_t fLay3Length;
  static const Double_t fLay3LadderLength;
  static const Double_t fLay3DetLongRadius;
  static const Double_t fLay3LaddTopCornerEnd;
  static const Double_t fLay3ZPlusEndLength;
  static const Int_t    fLay3Nladd;
  static const Int_t    fLay3Ndet;
  static const Double_t fLay3DetShortRadius;

  static const Double_t fLay4Rmin;
  static const Double_t fLay4Rmax;
  static const Double_t fLay4Length;
  static const Double_t fLay4LadderLength;
  static const Double_t fLay4LaddTopCornerEnd;
  static const Double_t fLay4ZPlusEndLength;
  static const Int_t    fLay4Nladd;
  static const Int_t    fLay4Ndet;
  static const Double_t fLay4DetShortRadius;
  static const Double_t fLay4DetLongRadius;

  static const Double_t fSegmentLength;
  static const Double_t fLadderWidth;
  static const Double_t fLadderHeight;
  static const Double_t fLadderBeamRadius;
  static const Double_t fLadderLa;
  static const Double_t fLadderHa;
  static const Double_t fLadderLb;
  static const Double_t fLadderHb;
  static const Double_t fLadderl;

  static const Double_t fBottomBeamAngle;
  static const Double_t fBeamSidePhi;

  static const Double_t fWaferThickness;
  static const Double_t fWaferWidth;
  static const Double_t fWaferLength;
  
  static const Double_t fHybridLength;
  static const Double_t fHybridWidth;
  static const Double_t fHybridThBridgeThick;

  static const Double_t fLadWaferSep;
  static const Double_t fPinSuppWidth;
  static const Double_t fPinSuppHeight;
  static const Double_t fPinSuppRmax;
  static const Double_t fPinR;
  static const Double_t fPinSuppLength;
  static const Double_t fPinSuppThickness;
  static const Double_t fPinSuppConeAngle;
  static const Double_t fPinDXminOnSensor;
  static const Double_t fPinPinDDXOnSensor;
  static const Double_t fPinDYOnSensor;

  static const Double_t fCoolPipeInnerDiam;
  static const Double_t fCoolPipeOuterDiam;
  static const Double_t fCoolPipeSuppHeight;
  static const Double_t fCoolPipeSuppMaxLength;
  static const Double_t fCoolPipeSuppWidthExt;
  static const Double_t fCoolPipeSuppWidthIn;
  static const Double_t fCoolPipeSuppHoleDiam;
  static const Double_t fCoolPipeSuppFulWidth;
  static const Double_t fCoolPipeSuppTongW;
  static const Double_t fCoolPipeSuppAngle;
  static const Double_t fCoolPipeSuppSlitL;
  static const Double_t fCoolPipeSuppAxeDist;
  static const Double_t fLay3CoolPipeSuppH;
  static const Double_t fLay4CoolPipeSuppH;

  static const Double_t fBTBthick;
  static const Double_t fBTBlength;
  static const Double_t fBTBwidth;
  static const Double_t fBTBaxisAtoBottom;
  static const Double_t fBTBaxisAtoBase;
  static const Double_t fRadiusAminBTB;
  static const Double_t fRadiusBminBTB;
  static const Double_t fBTBHoleLength;
  static const Double_t fBTBHolewidth;
  static const Double_t fBTBHoleRefX;
  static const Double_t fBTBHoleRefY;

  Double_t fLay3LaddShortRadius;
  Double_t fLay3LaddLongRadius;

  Double_t fLay4LaddShortRadius;
  Double_t fLay4LaddLongRadius;

  Double_t fLay3sensorZPos[6];
  Double_t fLay4sensorZPos[8];


  ClassDef(AliITSv11GeometrySDD,1) // ITS v11 SDD geometry
};





#endif
