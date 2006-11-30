// @(#) $Id$

#ifndef AliHLTVERTEXFINDER_H
#define AliHLTVERTEXFINDER_H

class AliHLTSpacePointData;
class AliHLTVertexData;
class AliHLTVertex;

#include "AliHLTVertexArray.h"

class AliHLTVertexFinder:public AliHLTVertexArray {
 private:
  Double_t fX;     // x
  Double_t fY;     // y
  Double_t fZ;     // z
  Double_t fPhi;   // phi
  Double_t fR;     // radius
  
  Double_t fXErr;  // x error
  Double_t fYErr;  // y error
  Double_t fZErr;  // z error
  Double_t fMWxy;  // xy weight

 public:
  AliHLTVertexFinder(); 
  AliHLTVertexFinder(AliHLTVertexFinder &vf) : AliHLTVertexArray(vf){;}
  virtual ~AliHLTVertexFinder();

  void Reset();
  void Read(Int_t ncluster, AliHLTSpacePointData* hits);
  void Analyze();
  void Write(AliHLTVertex *vertex) const;
  void Write(AliHLTVertexData *vertex) const;

  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetXErr() const {return fXErr;}
  Double_t GetYErr() const {return fYErr;}
  Double_t GetZErr() const {return fZErr;}
  Double_t GetPhi()  const {return fPhi;}
  Double_t GetR()    const {return fR;}
  Double_t GetXYWeight() const {return fMWxy;}
  void SetX(Double_t f) {fX=f;}
  void SetY(Double_t f) {fY=f;}
  void SetZ(Double_t f) {fZ=f;}
  void SetXErr(Double_t f) {fXErr=f;}
  void SetYErr(Double_t f) {fYErr=f;}
  void SetZErr(Double_t f) {fZErr=f;}

  void SetXYWeight(Double_t f) {fMWxy = f;}
 
  ClassDef(AliHLTVertexFinder,1)  // Vertex finder class
};

typedef AliHLTVertexFinder AliL3VertexFinder; // for backward compatibility

#endif
