// @(#) $Id$

#ifndef ALIL3VERTEX_H
#define ALIL3VERTEX_H

class AliHLTVertexData;

class AliHLTVertex {

  public:
    AliHLTVertex(); 

    virtual ~AliHLTVertex();

    void SetZero();
    void Read(const AliHLTVertexData *vertex);

    Double_t GetX()    const {return fX;}
    Double_t GetY()    const {return fY;}
    Double_t GetZ()    const {return fZ;}
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

  private:
    AliHLTVertex(const AliHLTVertex&){;}
    AliHLTVertex& operator=(const AliHLTVertex&){return *this;}

    Double_t fX;   //x 
    Double_t fY;   //y 
    Double_t fZ;   //z 
    Double_t fPhi; //phi
    Double_t fR;   //R
    Double_t fXErr; //error in x
    Double_t fYErr; //error in z
    Double_t fZErr; //error in y
    Double_t fMWxy; //weight
 
    ClassDef(AliHLTVertex,1)  // Vertex base class
};

typedef AliHLTVertex AliL3Vertex; // for backward compatibility

#endif
