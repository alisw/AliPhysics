#ifndef AliL3VERTEX_H
#define AliL3VERTEX_H

#include "AliL3RootTypes.h"
#include "AliL3VertexData.h"

class AliL3Vertex {
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
    AliL3Vertex(); 
    AliL3Vertex(AliL3Vertex&){;}
    virtual ~AliL3Vertex();
    void SetZero();
       
    void Read(AliL3VertexData *vertex);

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
 
    ClassDef(AliL3Vertex,1)  // Level3
};
#endif
