#ifndef ALIL3TRACK_H
#define ALIL3TRACK_H

#include <string.h>

#include "AliL3RootTypes.h"

class AliL3Track {
  
 private:

  Int_t fNHits;

  Int_t fMCid;  //Assigned id from MC data.
  Double_t   fKappa;  // Curvature
  Double_t   fRadius; // Radius of the helix (projected to a circle)
  Double_t   fCenterX; // x coordinate of the center of the helix (projected to a circle)
  Double_t   fCenterY; // y coordinate of the center of the helix (projected to a circle)
  Bool_t   fFromMainVertex; // true if tracks origin is the main vertex, otherwise false
  
  Int_t fRowRange[2]; //Subsector where this track was build
  Int_t fSector;      //Sector # where  this track was build

  //data from momentum fit
  Int_t    fQ;    //charge measured fit
    
  //track parameters:
  Double_t fPhi0; //azimuthal angle of the first point
  Double_t fPsi; //azimuthal angle of the momentum 
  Double_t fR0;  //radius of the first point
  Double_t fTanl; //tan of dipangle at (r,phi,z)
  Double_t fZ0;  //z coordinate of the first point
  Double_t fPt; //transverse momentum
  Double_t fLength;
  
  Double_t fPterr;
  Double_t fPsierr;
  Double_t fZ0err;
  Double_t fTanlerr;

  Double_t fFirstPoint[3];
  Double_t fLastPoint[3];
  Double_t fPoint[3];
  Double_t fPointPsi; //azimuthal angle of the momentum at Point
  Bool_t fIsPoint;    //Helix crosses the X-plane
  Bool_t IsPoint(Bool_t ispoint) {fIsPoint = ispoint;return fIsPoint;}
  
  Bool_t fIsLocal; //Track given in local coordinates.

  UInt_t fHitNumbers[174];  //Array of hit numbers for this track

 protected:

  static Float_t BFACT;
  static Float_t bField;
  static Double_t pi;
 
 public:
  
  AliL3Track();
  virtual ~AliL3Track();
  
  virtual void Set(AliL3Track* track);
  
  Bool_t CalculateReferencePoint(Double_t angle);//Calculate Reference Point
  Bool_t CalculateEdgePoint(Double_t angle);//Calculate crossing point with line
  Bool_t CalculatePoint(Double_t xplane);//Calculate crossing point with X-plane
  Bool_t IsPoint() {return fIsPoint;}
  void CalculateHelix();
  Double_t GetDistance(Double_t x0,Double_t x1){return 0;}

  Bool_t IsLocal() {return fIsLocal;}

  // getter
  Double_t GetFirstPointX() {return fFirstPoint[0];}
  Double_t GetFirstPointY() {return fFirstPoint[1];}
  Double_t GetFirstPointZ() {return fFirstPoint[2];}
  Double_t GetLastPointX() {return fLastPoint[0];}
  Double_t GetLastPointY() {return fLastPoint[1];}
  Double_t GetLastPointZ() {return fLastPoint[2];}

  Double_t GetPointPsi() {return fPointPsi;}
  Double_t GetPointX() {return fPoint[0];}
  Double_t GetPointY() {return fPoint[1];}
  Double_t GetPointZ() {return fPoint[2];}

  Double_t GetPt() const {return fPt;}
  Double_t GetTgl() const {return fTanl;}
  Double_t GetPhi0() const {return fPhi0;}
  Double_t GetPsi() const {return fPsi;}
  Double_t GetR0() const {return fR0;}
  Double_t GetZ0() const {return fZ0;}

  Double_t   GetKappa()            const { return fKappa;}
  Double_t   GetRadius()           const { return fRadius;}
  Double_t   GetCenterX()          const { return fCenterX;}
  Double_t   GetCenterY()          const { return fCenterY;}

  Int_t GetNHits() {return fNHits;}
  Int_t   GetNumberOfPoints()   const {return fNHits;}
  Bool_t  ComesFromMainVertex() const { return fFromMainVertex;}
    
  Double_t   GetPx()               const { return fPt*cos(fPsi);}
  Double_t   GetPy()               const { return fPt*sin(fPsi);}
  Double_t   GetPz()               const { return fPt*fTanl;}
  
  Double_t   GetP() const;
  Double_t   GetPseudoRapidity() const;
  Double_t   GetEta() const; 
  Double_t   GetRapidity() const;
  
  Int_t   GetCharge()           const { return fQ;}
  Int_t GetMCid() const {return fMCid;}
  Double_t GetLength()  const {return fLength;}

  Int_t GetFirstRow()  const {return fRowRange[0];}
  Int_t GetLastRow()  const {return fRowRange[1];}

  UInt_t *GetHitNumbers() {return fHitNumbers;}

  // setter   
  void SetMCid(Int_t f) {fMCid = f;}
  void SetFirstPoint(Double_t f,Double_t g,Double_t h) {fFirstPoint[0]=f; fFirstPoint[1]=g; fFirstPoint[2]=h;}
  void SetLastPoint(Double_t f,Double_t g,Double_t h) {fLastPoint[0]=f; fLastPoint[1]=g; fLastPoint[2]=h;}
  
  void SetHits(Int_t nhits,UInt_t *hits) {memcpy(fHitNumbers,hits,nhits*sizeof(UInt_t));}
  
  void SetPhi0(Double_t f) {fPhi0 = f;}
  void SetPsi(Double_t f) {fPsi = f;}
  void SetR0(Double_t f) {fR0 = f;}
  void SetTgl(Double_t f) {fTanl =f;}
  void SetZ0(Double_t f) {fZ0 = f;}
  void SetPt(Double_t f) {fPt = f;}
  void SetLength(Double_t f) {fLength = f;}
  void SetPterr(Double_t f) {fPterr = f;}
  void SetPsierr(Double_t f) {fPsierr = f;}
  void SetZ0err(Double_t f) {fZ0err = f;}
  void SetTglerr(Double_t f) {fTanlerr = f;}
  void SetKappa(Double_t f) {fKappa = f;}

  void SetNHits(Int_t f) {fNHits = f;}
    
  void  SetRowRange(Int_t f,Int_t g) {fRowRange[0]=f; fRowRange[1]=g;}
  void  SetSector(Int_t f) {fSector = f;}
  
  void   SetRadius(Double_t f)         { fRadius = f; }
  void   SetCenterX(Double_t f)        { fCenterX = f; }
  void   SetCenterY(Double_t f)        { fCenterY = f; }
  
  void   SetCharge(Int_t f)            { fQ = f; }
  
  void   ComesFromMainVertex(Bool_t f) { fFromMainVertex = f; }
  
  ClassDef(AliL3Track,1) //Conformal mapping track class
};
    
#endif
    
