#ifndef ALIITSV11GEOMCABLEFLAT_H
#define ALIITSV11GEOMCABLEFLAT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
//   Class for flat cables
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

class TGeoVolume;
class TGeoTranslation;
//class TGeoCombiTrans;
//class TGeoArb8;
//class TGeoNode;
class TGeoMedium;

#include "AliITSv11GeomCable.h"

class AliITSv11GeomCableFlat : public AliITSv11GeomCable {

 public:
  AliITSv11GeomCableFlat();
  AliITSv11GeomCableFlat(const char* name, Double_t width, Double_t thick);
  AliITSv11GeomCableFlat(const AliITSv11GeomCableFlat &source);
  AliITSv11GeomCableFlat& operator=(const AliITSv11GeomCableFlat &source);
  virtual ~AliITSv11GeomCableFlat() {};

  void        SetNLayers(Int_t nLayers);
  Int_t       SetLayer(Int_t nLayer,Double_t thick,TGeoMedium *medium,Int_t color=0);
  void        AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt,
			     Double_t *coord, Double_t *orthVect);
  TGeoVolume* CreateAndInsertCableSegment(Int_t p2, Double_t rotation=0);
  TGeoVolume* CreateAndInsertBoxCableSegment(Int_t p2, Double_t rotation=0);
  TGeoVolume* CreateAndInsertCableCylSegment(Int_t p2, Double_t rotation=0);

  void        SetWidth(Double_t width) { fWidth = width;};
  void        SetThickness(Double_t thick) {fThick = thick;};
  Double_t    GetWidth() const {return fWidth;};
  Double_t    GetThickness() const {return fThick;};

  virtual void  PrintCheckPoints() const;
  virtual Int_t GetPoint(Int_t iCheckPt, Double_t *coord) const;
  virtual Int_t GetVect(Int_t iCheckPt, Double_t *coord) const;

 protected:
  TGeoVolume *CreateSegment( Double_t *coord1,Double_t *coord2,
			     Double_t *localVect1, Double_t *localVect2 );

  TGeoVolume *CreateBoxSegment( Double_t *coord1,Double_t *coord2);

  TGeoVolume *CreateCylSegment( Double_t &phi, Double_t &r);

  Double_t  fWidth;                                 // width
  Double_t  fThick;                                 // total thickness
  Int_t     fNlayer;                                // number of layers
  Double_t  fPreviousX[3];                          // used internally
  Double_t         fLayThickness[fgkCableMaxLayer]; // layer thicknesses
  TGeoTranslation *fTranslation[fgkCableMaxLayer];  // layer translations
  TGeoMedium      *fLayMedia[fgkCableMaxLayer];     // layer media
  Int_t            fLayColor[fgkCableMaxLayer];     // layer colors

  ClassDef(AliITSv11GeomCableFlat,1)
};


#endif
