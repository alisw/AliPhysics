#ifndef ALIMUONMATHIESON_H
#define ALIMUONMATHIESON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONMathieson
/// \brief Implementation of Mathieson response

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONGeometrySegmentation;

class AliMUONMathieson : public TObject
{
 public:
    AliMUONMathieson();
    virtual ~AliMUONMathieson();
 
    /// Get anode cathode Pitch
    Float_t Pitch() const        {return fPitch;}
    // Set anode cathode Pitch
    void    SetPitch(Float_t p1);

    // Set Mathieson parameters
    //
    
    /// Mathieson \a sqrt{Kx3} and derived \a Kx2 and \a Kx4
    void    SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3);
    
    /// Mathieson \a sqrt{Kx3}
    void    SetSqrtKx3(Float_t p1) {fSqrtKx3 = p1;};
    
    /// Mathieson \a Kx2
    void    SetKx2(Float_t p1)      {fKx2 = p1;};
    
    /// Mathieson \a Kx4
    void    SetKx4(Float_t p1)      {fKx4 = p1;};
    
    /// Mathieson \a sqrt{Ky3} and derived \a Ky2 and \a Ky4
    void SetSqrtKy3AndDeriveKy2Ky4(Float_t SqrtKy3);
    
    /// Mathieson \a sqrt{Ky3}
    void    SetSqrtKy3(Float_t p1)   {fSqrtKy3 = p1;};
    
    /// Mathieson \a Ky2
    void    SetKy2(Float_t p1) {fKy2 = p1;};
    
    /// Mathieson \a Ky4
    void    SetKy4(Float_t p1) {fKy4 = p1;};
    
    /// \deprecated To be removed when old (s)digitizers go off.
    Float_t  IntXY(Int_t id, AliMUONGeometrySegmentation* segmentation) const;
    
    /// Charge integration on region \a (x1,y1,x2,y2).
    Float_t IntXY(Float_t xi1, Float_t yi1, Float_t xi2, Float_t yi2) const;
    
 private:
  
    Float_t fSqrtKx3;                  ///< Mathieson Sqrt(Kx3)
    Float_t fKx2;                      ///< Mathieson Kx2
    Float_t fKx4;                      ///< Mathieson Kx4 = Kx1/Kx2/Sqrt(Kx3)  
    Float_t fSqrtKy3;                  ///< Mathieson Sqrt(Ky3)
    Float_t fKy2;                      ///< Mathieson Ky2
    Float_t fKy4;                      ///< Mathieson Ky4 = Ky1/Ky2/Sqrt(Ky3)
    Float_t fPitch;                    ///< anode-cathode pitch
    Float_t fInversePitch;             ///< 1/Pitch

  ClassDef(AliMUONMathieson,3) // Implementation of Mathieson response
};
#endif











