
// Author: ruben.shahoyan@cern.ch   20/03/2007

///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  Wrapper for the set of mag.field parameterizations by Chebyshev polinomials  //
//  To obtain the field in cartesian coordinates/components use                  //
//    Field(float* xyz, float* bxyz);                                            //
//  For cylindrical coordinates/components:                                      //
//    FieldCyl(float* rphiz, float* brphiz)                                      //
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
//  GetTPCInt(float* xyz, float* bxyz);  for Cartesian frame                     //
//  or                                                                           //
//  GetTPCIntCyl(Float_t *rphiz, Float_t *b); for Cylindrical frame              //
//                                                                               //
//                                                                               //
//  The units are kiloGauss and cm.                                              //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

#ifndef ALIMAGFCHEB_H
#define ALIMAGFCHEB_H

#include <TMath.h>
#include <TNamed.h>
#include "AliCheb3D.h"

class TSystem;

class AliMagFCheb: public TNamed
{
 public:
  AliMagFCheb();
  AliMagFCheb(const AliMagFCheb& src);
  ~AliMagFCheb() {Clear();}
  //
  void       CopyFrom(const AliMagFCheb& src);
  AliMagFCheb& operator=(const AliMagFCheb& rhs);
  virtual void Clear(const Option_t * = "");
  //
  Int_t      GetNParamsSol()                              const {return fNParamsSol;}
  Int_t      GetNSegZSol()                                const {return fNSegZSol;}
  float*     GetSegZSol() const {return fSegZSol;}
  //
  Int_t      GetNParamsTPCInt()                           const {return fNParamsTPCInt;}
  Int_t      GetNSegZTPCInt()                             const {return fNSegZTPCInt;}
  //
  Int_t      GetNParamsDip()                              const {return fNParamsDip;}
  Int_t      GetNSegZDip()                                const {return fNZSegDip;}
  //
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
  Int_t      FindDipSegment(const float *xyz)             const;
  AliCheb3D* GetParamSol(Int_t ipar)                      const {return (AliCheb3D*)fParamsSol->UncheckedAt(ipar);}
  AliCheb3D* GetParamTPCInt(Int_t ipar)                   const {return (AliCheb3D*)fParamsTPCInt->UncheckedAt(ipar);}
  AliCheb3D* GetParamDip(Int_t ipar)                      const {return (AliCheb3D*)fParamsDip->UncheckedAt(ipar);}
  //
  virtual void Print(Option_t * = "")                     const;
  //
  virtual void Field(Float_t *xyz, Float_t *b)            const;
  virtual void FieldCyl(const Float_t *rphiz, Float_t *b) const;
  //
  virtual void GetTPCInt(Float_t *xyz, Float_t *b)        const;
  virtual void GetTPCIntCyl(Float_t *rphiz, Float_t *b)   const;
  //
  static void CylToCartCylB(const float *rphiz, const float *brphiz,float *bxyz);
  static void CylToCartCartB(const float *xyz,  const float *brphiz,float *bxyz);
  static void CartToCylCartB(const float *xyz,  const float *bxyz,  float *brphiz);
  static void CartToCylCylB(const float *rphiz, const float *bxyz,  float *brphiz);
  static void CartToCyl(const float *xyz,  float *rphiz);
  static void CylToCart(const float *rphiz,float *xyz);
  //
#ifdef  _INC_CREATION_ALICHEB3D_                          // see AliCheb3D.h for explanation
  void         LoadData(const char* inpfile);
  //
  AliMagFCheb(const char* inputFile);
  void       SaveData(const char* outfile)              const;
  Int_t      SegmentDipDimension(float** seg,const TObjArray* par,int npar, int dim, 
				 float xmn,float xmx,float ymn,float ymx,float zmn,float zmx);
  //
  void       AddParamSol(const AliCheb3D* param);
  void       AddParamTPCInt(const AliCheb3D* param);
  void       AddParamDip(const AliCheb3D* param);
  void       BuildTableDip();
  void       BuildTableSol();
  void       BuildTableTPCInt();
  void       ResetTPCInt();

#endif
  //
 protected:
  virtual void FieldCylSol(const Float_t *rphiz, Float_t *b)    const;
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
  ClassDef(AliMagFCheb,3)            // Wrapper class for the set of Chebishev parameterizations of Alice mag.field
  //
 };


//__________________________________________________________________________________________
inline void AliMagFCheb::FieldCyl(const Float_t *rphiz, Float_t *b) const
{
  // compute field in Cylindircal coordinates
  //  if (rphiz[2]<GetMinZSol() || rphiz[2]>GetMaxZSol() || rphiz[0]>GetMaxRSol()) {for (int i=3;i--;) b[i]=0; return;}
  FieldCylSol(rphiz,b);
}

//__________________________________________________________________________________________________
inline void AliMagFCheb::CylToCartCylB(const float *rphiz, const float *brphiz,float *bxyz)
{
  // convert field in cylindrical coordinates to cartesian system, point is in cyl.system
  float btr = TMath::Sqrt(brphiz[0]*brphiz[0]+brphiz[1]*brphiz[1]);
  float psiPLUSphi = TMath::ATan2(brphiz[1],brphiz[0]) + rphiz[1];
  bxyz[0] = btr*TMath::Cos(psiPLUSphi);
  bxyz[1] = btr*TMath::Sin(psiPLUSphi);
  bxyz[2] = brphiz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagFCheb::CylToCartCartB(const float *xyz, const float *brphiz,float *bxyz)
{
  // convert field in cylindrical coordinates to cartesian system, point is in cart.system
  float btr = TMath::Sqrt(brphiz[0]*brphiz[0]+brphiz[1]*brphiz[1]);
  float phiPLUSpsi = TMath::ATan2(xyz[1],xyz[0]) +  TMath::ATan2(brphiz[1],brphiz[0]);
  bxyz[0] = btr*TMath::Cos(phiPLUSpsi);
  bxyz[1] = btr*TMath::Sin(phiPLUSpsi);
  bxyz[2] = brphiz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagFCheb::CartToCylCartB(const float *xyz, const float *bxyz, float *brphiz)
{
  // convert field in cylindrical coordinates to cartesian system, poin is in cart.system
  float btr = TMath::Sqrt(bxyz[0]*bxyz[0]+bxyz[1]*bxyz[1]);
  float psiMINphi = TMath::ATan2(bxyz[1],bxyz[0]) - TMath::ATan2(xyz[1],xyz[0]);
  //
  brphiz[0] = btr*TMath::Cos(psiMINphi);
  brphiz[1] = btr*TMath::Sin(psiMINphi);
  brphiz[2] = bxyz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagFCheb::CartToCylCylB(const float *rphiz, const float *bxyz, float *brphiz)
{
  // convert field in cylindrical coordinates to cartesian system, point is in cyl.system
  float btr = TMath::Sqrt(bxyz[0]*bxyz[0]+bxyz[1]*bxyz[1]);
  float psiMINphi =  TMath::ATan2(bxyz[1],bxyz[0]) - rphiz[1];
  brphiz[0] = btr*TMath::Cos(psiMINphi);
  brphiz[1] = btr*TMath::Sin(psiMINphi);
  brphiz[2] = bxyz[2];
  //
}

//__________________________________________________________________________________________________
inline void AliMagFCheb::CartToCyl(const float *xyz,float *rphiz)
{
  rphiz[0] = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  rphiz[1] = TMath::ATan2(xyz[1],xyz[0]);
  rphiz[2] = xyz[2];
}

//__________________________________________________________________________________________________
inline void AliMagFCheb::CylToCart(const float *rphiz, float *xyz)
{
  xyz[0] = rphiz[0]*TMath::Cos(rphiz[1]);
  xyz[1] = rphiz[0]*TMath::Sin(rphiz[1]);
  xyz[2] = rphiz[2];
}

#endif
