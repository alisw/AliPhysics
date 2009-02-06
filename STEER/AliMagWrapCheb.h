
// Author: ruben.shahoyan@cern.ch   20/03/2007

///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  Wrapper for the set of mag.field parameterizations by Chebyshev polinomials  //
//  To obtain the field in cartesian coordinates/components use                  //
//    Field(double* xyz, double* bxyz);                                          //
//  For cylindrical coordinates/components:                                      //
//    FieldCyl(double* rphiz, double* brphiz)                                    //
//                                                                               //
//  The solenoid part is parameterized in the volume  R<500, -550<Z<550 cm       //
//                                                                               //
//  The region R<423 cm,  -343.3<Z<481.3 for 30kA and -343.3<Z<481.3 for 12kA    //
//  is parameterized using measured data while outside the Tosca calculation     //
//  is used (matched to data on the boundary of the measurements)                //
//                                                                               //
//  Two options are possible:                                                    //
//  1) _BRING_TO_BOUNDARY_ is defined in the AliCheb3D:                          //
//     If the querried point is outside of the validity region then the field    //
//     at the closest point on the fitted surface is returned.                   //
//  2) _BRING_TO_BOUNDARY_ is not defined in the AliCheb3D:                      //
//     If the querried point is outside of the validity region the return        //
//     value for the field components are set to 0.                              //
//                                                                               //
//  To obtain the field integral in the TPC region from given point to nearest   //
//  cathod plane (+- 250 cm) use:                                                //
//  GetTPCInt(double* xyz, double* bxyz);  for Cartesian frame                   //
//  or                                                                           //
//  GetTPCIntCyl(Double_t *rphiz, Double_t *b); for Cylindrical frame            //
//                                                                               //
//                                                                               //
//  The units are kiloGauss and cm.                                              //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

#ifndef ALIMAGWRAPCHEB_H
#define ALIMAGWRAPCHEB_H

#include <TMath.h>
#include <TNamed.h>
#include "AliCheb3D.h"

class TSystem;
class TArrayF;
class TArrayI;

class AliMagWrapCheb: public TNamed
{
 public:
  AliMagWrapCheb();
  AliMagWrapCheb(const AliMagWrapCheb& src);
  ~AliMagWrapCheb() {Clear();}
  //
  void       CopyFrom(const AliMagWrapCheb& src);
  AliMagWrapCheb& operator=(const AliMagWrapCheb& rhs);
  virtual void Clear(const Option_t * = "");
  //
  Int_t      GetNParamsSol()                              const {return fNParamsSol;}
  Int_t      GetNSegZSol()                                const {return fNSegZSol;}
  Float_t*   GetSegZSol() const {return fSegZSol;}
  //
  Int_t      GetNParamsTPCInt()                           const {return fNParamsTPCInt;}
  Int_t      GetNSegZTPCInt()                             const {return fNSegZTPCInt;}
  //
  Int_t      GetNParamsDip()                              const {return fNParamsDip;}
  Int_t      GetNSegZDip()                                const {return fNZSegDip;}
  //
  Float_t    GetMaxZ()                                    const {return GetMaxZSol();}
  Float_t    GetMinZ()                                    const {return fParamsDip ? GetMinZDip() : GetMinZSol();}
  //
  Float_t    GetMinZSol()                                 const {return fMinZSol;}
  Float_t    GetMaxZSol()                                 const {return fMaxZSol;}
  Float_t    GetMaxRSol()                                 const {return fMaxRSol;}
  //
  Float_t    GetMinZDip()                                 const {return fMinZDip;}
  Float_t    GetMaxZDip()                                 const {return fMaxZDip;}
  //
  Float_t    GetMinZTPCInt()                              const {return fMinZTPCInt;}
  Float_t    GetMaxZTPCInt()                              const {return fMaxZTPCInt;}
  Float_t    GetMaxRTPCInt()                              const {return fMaxRTPCInt;}
  //
  AliCheb3D* GetParamSol(Int_t ipar)                      const {return (AliCheb3D*)fParamsSol->UncheckedAt(ipar);}
  AliCheb3D* GetParamTPCInt(Int_t ipar)                   const {return (AliCheb3D*)fParamsTPCInt->UncheckedAt(ipar);}
  AliCheb3D* GetParamDip(Int_t ipar)                      const {return (AliCheb3D*)fParamsDip->UncheckedAt(ipar);}
  //
  virtual void Print(Option_t * = "")                     const;
  //
  virtual void Field(const Double_t *xyz, Double_t *b)    const;
  Double_t     GetBz(const Double_t *xyz)                 const;
  //
  void FieldCyl(const Double_t *rphiz, Double_t  *b)      const;
  void GetTPCInt(const Double_t *xyz, Double_t *b)        const;
  void GetTPCIntCyl(const Double_t *rphiz, Double_t *b)   const;
  //
  Int_t       FindDipSegment(const Double_t *xyz)         const; 
  static void CylToCartCylB(const Double_t *rphiz, const Double_t *brphiz,Double_t *bxyz);
  static void CylToCartCartB(const Double_t *xyz,  const Double_t *brphiz,Double_t *bxyz);
  static void CartToCylCartB(const Double_t *xyz,  const Double_t *bxyz,  Double_t *brphiz);
  static void CartToCylCylB(const Double_t *rphiz, const Double_t *bxyz,  Double_t *brphiz);
  static void CartToCyl(const Double_t *xyz,  Double_t *rphiz);
  static void CylToCart(const Double_t *rphiz,Double_t *xyz);
  //
#ifdef  _INC_CREATION_ALICHEB3D_                          // see AliCheb3D.h for explanation
  void         LoadData(const char* inpfile);
  //
  AliMagWrapCheb(const char* inputFile);
  void       SaveData(const char* outfile)                const;
  Int_t      SegmentDipDimension(Float_t** seg,const TObjArray* par,int npar, int dim, 
				 Float_t xmn,Float_t xmx,Float_t ymn,Float_t ymx,Float_t zmn,Float_t zmx);
  //
  void       AddParamSol(const AliCheb3D* param);
  void       AddParamTPCInt(const AliCheb3D* param);
  void       AddParamDip(const AliCheb3D* param);
  void       BuildTableDip();
  void       BuildTableSol();
  void       BuildTableTPCInt();
  void       ResetTPCInt();
  //
  //
#endif
  //
 protected:
  void     FieldCylSol(const Double_t *rphiz, Double_t *b)    const;
  Double_t FieldCylSolBz(const Double_t *rphiz)               const;
  //
 protected:
  //
  Int_t      fNParamsSol;            // Total number of parameterization pieces for Sol 
  Int_t      fNSegZSol;              // Number of segments in Z for Solenoid field
  //
  Int_t      fNParamsTPCInt;         // Total number of parameterization pieces for TPC field integral 
  Int_t      fNSegZTPCInt;           // Number of segments in Z for TPC field integral
  //
  Int_t      fNParamsDip;            // Total number of parameterization pieces for dipole 
  Int_t      fNZSegDip;              // number of distinct Z segments in Dipole
  Int_t      fNYSegDip;              // number of distinct Y segments in Dipole
  Int_t      fNXSegDip;              // number of distinct X segments in Dipole
  //
  Float_t*   fSegZSol;               //[fNSegZSol]      upper boundaries of Z segments
  Float_t*   fSegRSol;               //[fNParamsSol]    upper boundaries of R segments
  //
  Float_t*   fSegZTPCInt;            //[fNSegZTPCInt]    upper boundaries of Z segments
  Float_t*   fSegRTPCInt;            //[fNParamsTPCInt]  upper boundaries of R segments
  //
  Float_t*   fSegZDip;               //[fNZSegDip] coordinates of distinct Z segments in Dipole
  Float_t*   fSegYDip;               //[fNYSegDip] coordinated of Y segments for each Zsegment in Dipole
  Float_t*   fSegXDip;               //[fNXSegDip] coordinated of X segments for each Ysegment in Dipole
  //
  Int_t*     fNSegRSol;              //[fNSegZSol]      number of R segments for each Z segment
  Int_t*     fSegZIdSol;             //[fNSegZSol]      Id of the first R segment of each Z segment in the fSegRSol...
  //
  Int_t*     fNSegRTPCInt;           //[fNSegZTPCInt]   number of R segments for each Z segment
  Int_t*     fSegZIdTPCInt;          //[fNSegZTPCInt]   Id of the first R segment of each Z segment in the fSegRTPCInt...
  //
  Int_t*     fBegSegYDip;            //[fNZSegDip] beginning of Y segments array for each Z segment
  Int_t*     fNSegYDip;              //[fNZSegDip] number of Y segments for each Z segment
  Int_t*     fBegSegXDip;            //[fNYSegDip] beginning of X segments array for each Y segment
  Int_t*     fNSegXDip;              //[fNYSegDip] number of X segments for each Y segment
  Int_t*     fSegIDDip;              //[fNXSegDip] ID of the dipole parameterization for given XYZ segment
  //
  Float_t    fMinZSol;               // Min Z of Sol parameterization (in CYL. coordinates)
  Float_t    fMaxZSol;               // Max Z of Sol parameterization (in CYL. coordinates)
  Float_t    fMaxRSol;               // Max R of Sol parameterization (in CYL. coordinates)
  //
  Float_t    fMinZDip;               // Min Z of Dipole parameterization
  Float_t    fMaxZDip;               // Max Z of Dipole parameterization
  //
  Float_t    fMinZTPCInt;            // Min Z of TPCInt parameterization (in CYL. coordinates)
  Float_t    fMaxZTPCInt;            // Max Z of TPCInt parameterization (in CYL. coordinates)
  Float_t    fMaxRTPCInt;            // Max R of TPCInt parameterization (in CYL. coordinates)
  // 
  TObjArray* fParamsSol;             // Parameterization pieces for Solenoid field
  TObjArray* fParamsDip;             // Parameterization pieces for Dipole field
  TObjArray* fParamsTPCInt;          // Parameterization pieces for Solenoid field integrals in TPC region
  //
  ClassDef(AliMagWrapCheb,4)         // Wrapper class for the set of Chebishev parameterizations of Alice mag.field
  //
 };


//__________________________________________________________________________________________
inline void AliMagWrapCheb::FieldCyl(const Double_t *rphiz, Double_t *b) const
{
  // compute field in Cylindircal coordinates
  //  if (rphiz[2]<GetMinZSol() || rphiz[2]>GetMaxZSol() || rphiz[0]>GetMaxRSol()) {for (int i=3;i--;) b[i]=0; return;}
  FieldCylSol(rphiz,b);
}

//__________________________________________________________________________________________________
inline void AliMagWrapCheb::CylToCartCylB(const Double_t *rphiz, const Double_t *brphiz,Double_t *bxyz)
{
  // convert field in cylindrical coordinates to cartesian system, point is in cyl.system
  Double_t btr = TMath::Sqrt(brphiz[0]*brphiz[0]+brphiz[1]*brphiz[1]);
  Double_t psiPLUSphi = TMath::ATan2(brphiz[1],brphiz[0]) + rphiz[1];
  bxyz[0] = btr*TMath::Cos(psiPLUSphi);
  bxyz[1] = btr*TMath::Sin(psiPLUSphi);
  bxyz[2] = brphiz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagWrapCheb::CylToCartCartB(const Double_t* xyz, const Double_t *brphiz, Double_t *bxyz)
{
  // convert field in cylindrical coordinates to cartesian system, point is in cart.system
  Double_t btr = TMath::Sqrt(brphiz[0]*brphiz[0]+brphiz[1]*brphiz[1]);
  Double_t phiPLUSpsi = TMath::ATan2(xyz[1],xyz[0]) +  TMath::ATan2(brphiz[1],brphiz[0]);
  bxyz[0] = btr*TMath::Cos(phiPLUSpsi);
  bxyz[1] = btr*TMath::Sin(phiPLUSpsi);
  bxyz[2] = brphiz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagWrapCheb::CartToCylCartB(const Double_t *xyz, const Double_t *bxyz, Double_t *brphiz)
{
  // convert field in cylindrical coordinates to cartesian system, poin is in cart.system
  Double_t btr = TMath::Sqrt(bxyz[0]*bxyz[0]+bxyz[1]*bxyz[1]);
  Double_t psiMINphi = TMath::ATan2(bxyz[1],bxyz[0]) - TMath::ATan2(xyz[1],xyz[0]);
  //
  brphiz[0] = btr*TMath::Cos(psiMINphi);
  brphiz[1] = btr*TMath::Sin(psiMINphi);
  brphiz[2] = bxyz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagWrapCheb::CartToCylCylB(const Double_t *rphiz, const Double_t *bxyz, Double_t *brphiz)
{
  // convert field in cylindrical coordinates to cartesian system, point is in cyl.system
  Double_t btr = TMath::Sqrt(bxyz[0]*bxyz[0]+bxyz[1]*bxyz[1]);
  Double_t psiMINphi =  TMath::ATan2(bxyz[1],bxyz[0]) - rphiz[1];
  brphiz[0] = btr*TMath::Cos(psiMINphi);
  brphiz[1] = btr*TMath::Sin(psiMINphi);
  brphiz[2] = bxyz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagWrapCheb::CartToCyl(const Double_t *xyz, Double_t *rphiz)
{
  rphiz[0] = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  rphiz[1] = TMath::ATan2(xyz[1],xyz[0]);
  rphiz[2] = xyz[2];
}

//__________________________________________________________________________________________________
inline void AliMagWrapCheb::CylToCart(const Double_t *rphiz, Double_t *xyz)
{
  xyz[0] = rphiz[0]*TMath::Cos(rphiz[1]);
  xyz[1] = rphiz[0]*TMath::Sin(rphiz[1]);
  xyz[2] = rphiz[2];
}

#endif
