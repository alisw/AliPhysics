#ifndef ALIMUONMATHIESON_H
#define ALIMUONMATHIESON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

class AliSegmentation;
class AliMUONGeometrySegmentation;

class AliMUONMathieson 
{
 public:
    AliMUONMathieson();
    virtual ~AliMUONMathieson(){}
 
    // Get anode cathode Pitch
    Float_t Pitch() const        {return fPitch;}
    // Set anode cathode Pitch
    void    SetPitch(Float_t p1) {fPitch = p1;};

    // Set Mathieson parameters
    // Mathieson \sqrt{Kx3} and derived Kx2 and Kx4
    void SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3);
    // Mathieson \sqrt{Kx3}
    void    SetSqrtKx3(Float_t p1) {fSqrtKx3 = p1;};
    // Mathieson Kx2
    void    SetKx2(Float_t p1)      {fKx2 = p1;};
    // Mathieson Kx4
    void    SetKx4(Float_t p1)      {fKx4 = p1;};
    // Mathieson \sqrt{Ky3} and derived Ky2 and Ky4
    void SetSqrtKy3AndDeriveKy2Ky4(Float_t SqrtKy3);
    // Mathieson \sqrt{Ky3}
    void    SetSqrtKy3(Float_t p1)   {fSqrtKy3 = p1;};
    // Mathieson Ky2
    void    SetKy2(Float_t p1) {fKy2 = p1;};
    // Mathieson Ky4
    void    SetKy4(Float_t p1) {fKy4 = p1;};
    // Charge disintegration
    Float_t  IntXY(AliSegmentation * segmentation);
    Float_t  IntXY(Int_t id, AliMUONGeometrySegmentation* segmentation);

    ClassDef(AliMUONMathieson,1) // Implementation of Mathieson response
 protected:
  
    Float_t fSqrtKx3;                  // Mathieson Sqrt(Kx3)
    Float_t fKx2;                      // Mathieson Kx2
    Float_t fKx4;                      // Mathieson Kx4 = Kx1/Kx2/Sqrt(Kx3)  
    Float_t fSqrtKy3;                  // Mathieson Sqrt(Ky3)
    Float_t fKy2;                      // Mathieson Ky2
    Float_t fKy4;                      // Mathieson Ky4 = Ky1/Ky2/Sqrt(Ky3)
    Float_t fPitch;                    // anode-cathode pitch
};
#endif











