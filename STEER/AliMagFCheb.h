#ifndef ALIMAGFCHEB_H
#define ALIMAGFCHEB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


// Author: ruben.shahoyan@cern.ch   20/03/2007
//
///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  Wrapper for the set of mag.field parameterizations by Chebyshev polinomials  //
//  To obtain the field in cartesian coordinates/components use                  //
//    Field(float* xyz, float* bxyz);                                            //
//  For cylindrical coordinates/components:                                      //
//    FieldCyl(float* rphiz, float* brphiz)                                      //
//                                                                               //
//  For the moment only the solenoid part is parameterized in the volume defined //
//  by R<500, -550<Z<550 cm                                                      //
//                                                                               //
//  The region R<423 cm,  -343.3<Z<481.3 for 30kA and -343.3<Z<481.3 for 12kA    //
//  is parameterized using measured data while outside the Tosca calculation     //
//  is used (matched to data on the boundary of the measurements)                //
//                                                                               //
//  If the querried point is outside the validity region no the return values    //
//  for the field components are set to 0.                                       //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////


#include <TSystem.h>
#include <TNamed.h>
#include "AliCheb3D.h"
#include "AliCheb3DCalc.h"

class AliMagFCheb: public TNamed
{
 public:
    AliMagFCheb();
    AliMagFCheb(const char* inputFile);
    AliMagFCheb(const AliMagFCheb &src);
    AliMagFCheb& operator= (const AliMagFCheb &rhs);
    
   ~AliMagFCheb();
  //
  void       AddParamSol(AliCheb3D* param);
  void       AddParamDip(AliCheb3D* param);
  void       BuildTableSol();
  //
  Int_t      GetNParamsSol()                              const {return fNParamsSol;}
  Int_t      GetNSegZSol()                                const {return fNSegZSol;}
  Int_t      GetNSegRSol(int iz)                          const {return iz<fNParamsSol ? fNSegRSol[iz]:0;}
  Int_t      GetSegIDSol(int iz,int ir)                   const {return iz<fNParamsSol&&ir<fNSegRSol[iz] ? fSegZIdSol[iz]+ir:-1;}
  //
  Float_t    GetMinZSol()                                 const {return fMinZSol;}
  Float_t    GetMaxZSol()                                 const {return fMaxZSol;}
  Float_t    GetMaxRSol()                                 const {return fMaxRSol;}
  AliCheb3D* GetParamSol(Int_t ipar)                      const {return (AliCheb3D*)fParamsSol->UncheckedAt(ipar);}
  AliCheb3D* GetParamDip(Int_t ipar)                      const {return (AliCheb3D*)fParamsDip->UncheckedAt(ipar);}
  //
  void         LoadData(const char* inpfile);
  //
  virtual void Print(Option_t * = "")                     const;
  //
  virtual void Field(Float_t *xyz, Float_t *b)            const;
  virtual void FieldCyl(Float_t *rphiz, Float_t *b)       const;
  //
  //
#ifdef  _INC_CREATION_ALICHEB3D_                          // see AliCheb3D.h for explanation
  void         SaveData(const char* outfile)              const;
#endif
  //
 protected:
  void         Init0();
  virtual void FieldCylSol(Float_t *rphiz, Float_t *b)    const;
  //
 protected:
  //
  Int_t      fNParamsSol;            // Total number of parameterization pieces for Sol 
  Int_t      fNSegZSol;              // Number of segments is Z
  //
  Int_t      fNParamsDip;            // Total number of parameterization pieces for dipole 
  //
  Float_t*   fSegZSol;               //[fNSegZSol]       upper boundaries of Z segments
  Float_t*   fSegRSol;               //[fNParamsSol]     upper boundaries of R segments
  //
  Int_t*     fNSegRSol;              //[fNSegZSol]       number of R segments for each Z segment
  Int_t*     fSegZIdSol;             //[fNSegZSol]       Id of the first R segment of each Z segment in the fSegRSol...
  //
  Float_t    fMinZSol;               // Min Z of Sol parameterization (in CYL. coordinates)
  Float_t    fMaxZSol;               // Max Z of Sol parameterization (in CYL. coordinates)
  Float_t    fMaxRSol;               // Max R of Sol parameterization (in CYL. coordinates)
  //
  TObjArray* fParamsSol;             // Parameterization pieces for Solenoid field
  TObjArray* fParamsDip;             // Parameterization pieces for Dipole field
  //
  ClassDef(AliMagFCheb,1)            // Wrapper class for the set of Chebishev parameterizations of Alice mag.field
  //
 };


//__________________________________________________________________________________________
inline void AliMagFCheb::FieldCyl(Float_t *rphiz, Float_t *b) const
{
  // compute field in Cylindircal coordinates
  if (rphiz[2]<GetMinZSol() || rphiz[2]>GetMaxZSol() || rphiz[0]>GetMaxRSol()) {for (int i=3;i--;) b[i]=0; return;}
  FieldCylSol(rphiz,b);
}

#endif
