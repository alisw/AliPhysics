#ifndef ALIITSV11GEOMETRYSDD_H
#define ALIITSV11GEOMETRYSDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Ludovic Gaudichet (Ludovic.Gaudichet@to.infn.it)

class TGeoVolume;
class TGeoCombiTrans;
class TGeoArb8;



class AliITSv11GeometrySDD : public AliITSv11Geometry {

 public:
  AliITSv11GeometrySDD();
  AliITSv11GeometrySDD(Int_t debug);
  virtual ~AliITSv11GeometrySDD(){};


  virtual void SetGeomParameters();
  virtual void Layer3(TGeoVolume *Moth);
  virtual void Layer4(TGeoVolume *Moth);

  // Insert only the specified element
  void AddOnlyLadder3(Int_t i) { fAddOnlyLadder3=i; };
  void AddOnlyLadder4(Int_t i) { fAddOnlyLadder4=i; };
  void AddOnlySegment(Int_t i) { fAddOnlySegment=i; };

  private:

  // Create virtual volumes of a ladder and its detectors (inside layers 3 and 4)
  virtual TGeoVolume* CreateLay3Ladder();
  virtual TGeoVolume* CreateLay3Detectors();
  virtual TGeoVolume* CreateLay4Ladder();
  virtual TGeoVolume* CreateLay4Detectors();

  // Create virtual volumes inside a ladder volume
  virtual TGeoVolume* CreateLadderSegment();
  virtual TGeoVolume* CreateEndLadder(Double_t length);
  virtual TGeoVolume* CreateHybrid();


  // Create a TGeoCombiTrans, allows a general rotation in dPhi and dy,dz translations 
  TGeoCombiTrans* CreateCombiTrans(const char *name, Double_t dy, Double_t dz, Double_t dphi);

  // add (dx,dy,dz) translation to a initial TGeoCombiTrans
  void AddTranslationTotCombiTrans(TGeoCombiTrans* ct, Double_t dx=0, Double_t dy=0, Double_t dz=0);

  // Create one side of the CF corner of the CF structure
  TGeoArb8* CreateLadderSide(Double_t dz, Double_t angle, Double_t xSign,
			  Double_t L, Double_t H, Double_t l);

  // Create de CF structure in the ladder segment
  void AddLadderCFstruct(Double_t dy, TGeoVolume* vol);

  // Create a pin support object
  virtual TGeoVolume* CreatePinSupport(Double_t rotY);

  // Create a SDD detector object
  virtual TGeoVolume* CreateSDDsensor();


  //**************************************************

  Int_t fAddOnlyLadder3;
  Int_t fAddOnlyLadder4;
  Int_t fAddOnlySegment;


  Double_t fLay3Rmin;
  Double_t fLay3Rmax;
  Double_t fLay3Length;
  Double_t fLay3LadderLength;
  Double_t fLay3LaddShortRadius;
  Double_t fLay3LaddLongRadius;
  Double_t fLay3DetShortRadius;
  Double_t fLay3DetLongRadius;
  Double_t fLay3LaddTopCornerEnd;
  Int_t    fLay3Nladd;
  Int_t    fLay3Ndet;
  Double_t fLay3ZPlusEndLength;
  Double_t fLay3sensorZPos[6];


  Double_t fLay4Rmin;
  Double_t fLay4Rmax;
  Double_t fLay4Length;
  Double_t fLay4LaddShortRadius;
  Double_t fLay4LaddLongRadius;
  Double_t fLay4DetShortRadius;
  Double_t fLay4DetLongRadius;
  Double_t fLay4LadderLength;
  Double_t fLay4LaddTopCornerEnd;
  Int_t    fLay4Nladd;
  Int_t    fLay4Ndet;
  Double_t fLay4ZPlusEndLength;
  Double_t fLay4sensorZPos[8];

  Double_t fSegmentLength;
  Double_t fLadderWidth;
  Double_t fLadderHeight;
  Double_t fLadderBeamRadius;
  Double_t fLadderLa;
  Double_t fLadderHa;
  Double_t fLadderLb;
  Double_t fLadderHb;
  Double_t fLadderl;

  Double_t fBottomBeamAngle;
  Double_t fBeamSidePhi;

  Double_t fWaferThickness;
  Double_t fWaferWidth;
  Double_t fWaferLength;
  
  Double_t fHybridLength;
  Double_t fHybridWidth;

  Double_t fLadWaferSep;
  Double_t fPinSuppWidth;
  Double_t fPinSuppHeight;
  Double_t fPinSuppRmax;
  Double_t fPinR;
  Double_t fPinSuppLength;
  Double_t fPinSuppThickness;
  Double_t fPinSuppConeAngle;

  ClassDef(AliITSv11GeometrySDD,1) // ITS v11 SDD geometry
};



#endif
