#ifndef ALIITSV11GEOMCABLEROUND_H
#define ALIITSV11GEOMCABLEROUND_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//*************************************************************************
//   Class for round cables
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************

class TGeoVolume;
class TGeoMedium;

#include "AliITSv11GeomCable.h"

class AliITSv11GeomCableRound : public AliITSv11GeomCable {

 public:
  AliITSv11GeomCableRound(const char* name, Double_t radius);
  AliITSv11GeomCableRound(const AliITSv11GeomCableRound &source);
  AliITSv11GeomCableRound& operator=(const AliITSv11GeomCableRound &source);
  virtual ~AliITSv11GeomCableRound() {};

  virtual Int_t GetPoint(Int_t iCheckPt, Double_t *coord) const;
  virtual Int_t GetVect(Int_t iCheckPt, Double_t *coord) const;

  void          AddCheckPoint( TGeoVolume *vol, Int_t iCheckPt,
			       Double_t *coord, Double_t *orthVect);
  Int_t         CreateAndInsertCableSegment(Int_t p2);
  Int_t         CreateAndInsertTorusSegment(Int_t p2, Double_t rotation=0);
  void          PrintCheckPoints() const;

  void          SetNLayers(Int_t nLayers);
  Int_t         SetLayer(Int_t nLayer,Double_t thick,TGeoMedium *medium,
			 Int_t color=0);
  void          SetPhi(Double_t phi1, Double_t phi2)
                      {fPhiMin=phi1; fPhiMax=phi2;};

 protected:
  TGeoVolume*   CreateSegment( Double_t *coord1,Double_t *coord2,
			       Double_t *localVect1, Double_t *localVect2, Int_t p);
  TGeoVolume*   CreateTorus(  Double_t &phi, Double_t &r, Int_t p);

  Double_t   fRadius;                         // total radius
  Int_t      fNlayer;                         // number of layers
  Double_t   fPhiMin;                         // minimum phi
  Double_t   fPhiMax;                         // maximum phi
  Double_t   fLayThickness[fgkCableMaxLayer]; // layer thicknesses
  Int_t      fLayColor[fgkCableMaxLayer];     // layer colors
  TGeoMedium *fLayMedia[fgkCableMaxLayer];    // layer media

  ClassDef(AliITSv11GeomCableRound,1)
};


#endif
