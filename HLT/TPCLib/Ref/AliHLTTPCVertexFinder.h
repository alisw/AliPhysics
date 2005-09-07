// @(#) $Id$

#ifndef AliHLTTPCVERTEXFINDER_H
#define AliHLTTPCVERTEXFINDER_H

#include "AliHLTTPCVertexArray.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertexData.h"

class AliHLTTPCVertex;

class AliHLTTPCVertexFinder:public AliHLTTPCVertexArray{
  private:
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fPhi;
    Double_t fR;

    Double_t fXErr;
    Double_t fYErr;
    Double_t fZErr;

    Double_t fMWxy;

  public:
    AliHLTTPCVertexFinder(); 
    AliHLTTPCVertexFinder(AliHLTTPCVertexFinder&){;}
    virtual ~AliHLTTPCVertexFinder();

    void Reset();
    void Read(Int_t ncluster, AliHLTTPCSpacePointData* hits);
    void Analyze();
    void Write(AliHLTTPCVertex *vertex);
    void Write(AliHLTTPCVertexData *vertex);

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
 
    ClassDef(AliHLTTPCVertexFinder,1)  // Vertex finder class
};
#endif
