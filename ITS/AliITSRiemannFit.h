#ifndef ALIITSRIEMANNFIT_H
#define ALIITSRIEMANNFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/////////////////////////////////////////////////////////////////////
// Class for helix fit on the Riemann sphere                       //
/////////////////////////////////////////////////////////////////////
#include<TLorentzVector.h>
class TTree;
class TVector3;

class AliITSRiemannFit : public TObject{
 public:
  AliITSRiemannFit();
  AliITSRiemannFit(Int_t size,Int_t ntracks);
  AliITSRiemannFit(const AliITSRiemannFit& rec);
  AliITSRiemannFit& operator=(const AliITSRiemannFit &source);

  ~AliITSRiemannFit();
  class  AliPointtl {
    public :
      AliPointtl();
    // getters
      Int_t GetLay() const {return fLay;}
      Int_t GetLad() const {return fLad;}
      Int_t GetDet() const {return fDet;}
      Int_t GetTrack() const {return fTrack;}
      Float_t GetX() const {return fx;}
      Float_t GetY() const {return fy;}
      Float_t GetZ() const {return fz;}
      Float_t GetR() const {return fr;}
      Float_t GetdE() const {return fdE;}
      Float_t GetdX() const {return fdx;}
      Float_t GetdY() const {return fdy;}
      Float_t GetdZ() const {return fdz;}
      TLorentzVector* GetOrigin() const {return fOrigin;}
      TLorentzVector* GetMomentum() const {return fMomentum;}
      Int_t GetCode() const {return fCode;}
      Char_t* GetName() {return fName;}
      Float_t GetPt() const {return fPt;}
      Float_t GetPhi() const {return fPhi;}
      Float_t GetEta() const {return fEta;}
      Float_t GetVertexPhi() const {return fVertexPhi;}
      //setters
      void SetLay(Int_t l=0) { fLay = l;}
      void SetLad(Int_t l=0) { fLad = l;}
      void SetDet(Int_t d=0) { fDet = d;}
      void SetTrack(Int_t t=0) { fTrack = t;}
      void SetX(Float_t x=0) { fx = x;}
      void SetY(Float_t y=0) { fy = y;}
      void SetZ(Float_t z=0) { fz = z;}
      void SetR(Float_t r=0) { fr = r;}
      void SetdE(Float_t de=0) { fdE = de;}
      void SetdX(Float_t dx=0) { fdx = dx;}
      void SetdY(Float_t dy=0) { fdy = dy;}
      void SetdZ(Float_t dz=0) { fdz = dz;}
      void SetOrigin(TLorentzVector *ori=0) { fOrigin = ori;}
      void SetMomentum(TLorentzVector *mo=0) { fMomentum = mo;}
      void SetCode(Int_t c=0) { fCode = c;}
      void SetName(Char_t *n=0) { fName = n;}
      void SetPt(Float_t pt=0) { fPt = pt;}
      void SetPhi(Float_t phi=0) { fPhi = phi;}
      void SetEta(Float_t eta=0) { fEta = eta;}
      void SetVertexPhi(Float_t vert=0) { fVertexPhi = vert;}
    private :
      // copy constructor (NO copy ctr. allowed)
      AliPointtl(const AliPointtl& ap);
    // assignment operator (NO assignment allowed)
    AliPointtl& operator=(const AliPointtl& ap){
       this->~AliPointtl(); new(this) AliPointtl(ap);return *this;
    }
      Int_t fLay,fLad,fDet,fTrack;       // layer,ladder,detector and track
      Float_t fx,fy,fz,fr;               // global position of point 
      Float_t fdE,fdx,fdy,fdz;               // Errors
      TLorentzVector* fOrigin;    // position and momentum of 
      TLorentzVector* fMomentum;  //  particle at its origin
      Int_t fCode;                       // Geant code of particle
      Char_t *fName;                    // name
      Float_t fPt;                       // Pt at the origin
      Float_t fPhi,fEta,fVertexPhi;         // phi eta on layer and phi on vertex
  };
  Int_t GetSize() const {return this->fSizeEvent;}
  Int_t GetPrimaryTracks() const {return this->fPrimaryTracks;}
  Int_t GetPoints() const {return this->fPoints;}
  Int_t GetParticles() const {return this->fParticles;}
  Int_t GetLayPoints(Int_t layer) const {return this->fPLay[layer-1];}
  AliPointtl **GetPointRecs() const {return this->fPointRecs;}
  Float_t GetX(Int_t i) const {return this->fPointRecs[i]->GetX();}
  Float_t GetY(Int_t i) const {return this->fPointRecs[i]->GetY();}
  Float_t GetZ(Int_t i) const {return this->fPointRecs[i]->GetZ();}
  Float_t GetdX(Int_t i) const {return this->fPointRecs[i]->GetdX();}
  Float_t GetdY(Int_t i) const {return this->fPointRecs[i]->GetdY();}
  Float_t GetdZ(Int_t i) const {return this->fPointRecs[i]->GetdZ();}
  
  void     InitPoints(Int_t ntracks,TTree *TR,Int_t nparticles);
  void     WritePoints(void);
  void     ReadPoints(void);
  static Int_t SolveCubic(Double_t a,Double_t b,Double_t c,Double_t& x1,Double_t& x2,Double_t& x3);
  Int_t FitHelix(Int_t tracknumber,Double_t Px,Double_t Py,Double_t Pz,
		 Double_t& fd0,Double_t& fphi,Double_t& u0, Double_t& v0, Double_t& rho,
		 Double_t& omega, Double_t& z0,
		 Double_t& vpar,Double_t& chisql,Double_t& fCorrLin,Double_t& fFit,
		 Int_t first=1,Int_t second=1,Int_t third=1,Int_t fourth=1,Int_t fifth=1,Int_t sixth=1);  
 Int_t FitHelix(Int_t NPoints, TVector3** fPointRecs,
		 TVector3** fPointRecErrors,Float_t& f1, 
		 Float_t& f2, Float_t& f3);
 Int_t LinearFit(Int_t npoints, TVector3 **input, 
		  TVector3 **errors, Double_t omega,
		  Double_t &thu0, Double_t &thv0, Double_t &phi,TVector2 &zData, TVector3 &zError, 
		  Double_t &corrLin);

 private:

  static Double_t Fitfunction(Double_t *x, Double_t* par);

  Int_t fSizeEvent;      // size of array 
  Int_t fPrimaryTracks;  // number of primary tracks in the event
  Int_t fPoints;         // number of Reconstructed Points in the event
  Int_t fParticles;      // number of particles in the event
  Int_t fPLay[6];           // number of points in each layer
  AliPointtl **fPointRecs;    //rec points

  
  ClassDef(AliITSRiemannFit,1)  // Fast fit of helices on ITS RecPoints
    };
#endif

