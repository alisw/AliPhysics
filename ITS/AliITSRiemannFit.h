#ifndef ALIITSRIEMANNFIT_H
#define ALIITSRIEMANNFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TLorentzVector.h"
#include "TTree.h"
#include "AliITS.h"
#include "TVector3.h"

struct Point_tl{
  Int_t lay,lad,det,track;
  Float_t fx,fy,fz,fr;               // global position of point 
  Float_t fdE,fdx,fdy,fdz;               // Errors
  TLorentzVector fOrigin,fMomentum;  // position and momentum of 
                                     //  particle at its origin
  Int_t fCode;                       // Geant code of particle
  const Char_t *fName;
  Float_t fPt;                       // Pt at the origin
  Float_t phi,eta,vertexPhi;         // phi eta on layer and phi on vertex
};


class AliITSRiemannFit : public TObject{
 public:
  AliITSRiemannFit();
  AliITSRiemannFit(Int_t size,Int_t ntracks);
  ~AliITSRiemannFit();
  
  Int_t GetSize() const {return this->fSizeEvent;}
  Int_t GetPrimaryTracks() const {return this->fPrimaryTracks;}
  Int_t GetPoints() const {return this->fPoints;}
  Int_t GetParticles() const {return this->fParticles;}
  Int_t GetLayPoints(Int_t layer) const {return this->fPLay[layer-1];}
  Point_tl **GetPointRecs() const {return this->fPointRecs;}
  Float_t GetX(Int_t i) const {return this->fPointRecs[i]->fx;}
  Float_t GetY(Int_t i) const {return this->fPointRecs[i]->fy;}
  Float_t GetZ(Int_t i) const {return this->fPointRecs[i]->fz;}
  Float_t GetdX(Int_t i) const {return this->fPointRecs[i]->fdx;}
  Float_t GetdY(Int_t i) const {return this->fPointRecs[i]->fdy;}
  Float_t GetdZ(Int_t i) const {return this->fPointRecs[i]->fdz;}
  
  void     InitPoints(Int_t ntracks,AliITS *ITS,TTree *TR,Int_t nparticles);
  void     WritePoints(void);
  void     ReadPoints(void);
  static Int_t SolveCubic(Double_t a,Double_t b,Double_t c,Double_t&,Double_t&,Double_t&);
  Int_t FitHelix(Int_t tracknumber,Double_t Px,Double_t Py,Double_t Pz,
		 Double_t& fd0,Double_t& fphi,Double_t& u0, Double_t& v0, Double_t& rho,
		 Double_t& omega, Double_t& z0,
		 Double_t& vpar,Double_t& chisql,Double_t& fCorrLin,Double_t& fFit,
		 Int_t first=1,Int_t second=1,Int_t third=1,Int_t fourth=1,Int_t fifth=1,Int_t sixth=1);  
 private:
  Int_t fSizeEvent;      // size of array 
  Int_t fPrimaryTracks;  // number of primary tracks in the event
  Int_t fPoints;         // number of Reconstructed Points in the event
  Int_t fParticles;      // number of particles in the event
  Int_t fPLay[6];           // number of points in each layer
  Point_tl **fPointRecs;
  //
  // test erase
/*    Point_tl **fspdi,**fspdo; // This are for the first two layers and vertex analysis */
  
  ClassDef(AliITSRiemannFit,1)  // Fast fit of helices on ITS RecPoints
    };
#endif

