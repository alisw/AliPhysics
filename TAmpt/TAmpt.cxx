/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// This class implements an interface to the Ampt event generator

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h>
#include "Acommon.h"
#include "TAmpt.h"
#include "AliAmptRndm.h"

#ifndef WIN32
# define amptset amptset_
# define amptsetdef amptsetdef_
# define ampt ampt_
# define invflv invflv_
# define profile profile_
# define rluget_ampt rluget_ampt_
# define rluset_ampt rluset_ampt_
# define type_of_call
#else
# define amptset AMPTSET
# define amptsetdef AMPTSETDEF
# define ampt AMPT
# define invflv INVFLV
# define profile PROFILE
# define rluget_ampt RLUGET_AMPT
# define rluset_ampt RLUSET_AMPT
# define type_of_call _stdcall
#endif

#ifndef WIN32
extern "C" void type_of_call amptset(Double_t & , const char *, 
                                     const char *, const char *,
                                     Int_t & , Int_t &, Int_t &,
                                     Int_t &, const int, 
                                     const int, const int);
extern "C" void type_of_call amptsetdef();
extern "C" void type_of_call ampt(const char *, Double_t  &,
                                  Double_t &, const int);
extern "C" Int_t type_of_call invflv(Int_t &);
extern "C" float type_of_call profile(Float_t &);
extern "C" void type_of_call rluget_ampt(Int_t & lfn, Int_t & move);
extern "C" void type_of_call rluset_ampt(Int_t & lfn, Int_t & move);
#else
#endif


ClassImp(TAmpt)

//______________________________________________________________________________
TAmpt::TAmpt()
  : TGenerator("Ampt","Ampt"),
    fEfrm(5500.),
    fFrame("CMS"),
    fProj("A"),
    fTarg("A"),
    fIap(208),
    fIzp(82),
    fIat(208),
    fIzt(82),
    fBmin(0.),
    fBmax(5.)
{
  // Default constructor 
  amptsetdef();
}

//______________________________________________________________________________
TAmpt::TAmpt(Double_t efrm, const char *frame="CMS", 
             const char *proj="A", const char *targ="A", 
             Int_t iap=207, Int_t izp=82, Int_t iat=207, Int_t izt=82,
             Double_t bmin=0, Double_t bmax=20)
  : TGenerator("Ampt","Ampt"),
    fEfrm(efrm),
    fFrame(frame),
    fProj(proj),
    fTarg(targ),
    fIap(iap),
    fIzp(izp),
    fIat(iat),
    fIzt(izt),
    fBmin(bmin),
    fBmax(bmax)
{
  // TAmpt constructor: 
  // Note that there may be only one functional TAmpt object
  // at a time, so it's not useful to create more than one 
  // instance of it.
  amptsetdef();
}

//______________________________________________________________________________
TAmpt::~TAmpt()
{
  // Destructor
}

//______________________________________________________________________________
TObjArray* TAmpt::ImportParticles(Option_t */*option*/)
{
  //  Import created particles.

  fParticles->Clear();
  Int_t numpart = HBT.nlast;
  printf("TAmpt: AMPT stack contains %d particles.\n", numpart);

  for (Int_t i=0; i < numpart; ++i) {
    Int_t status = 1;
    Double_t px = HBT.plast[i][0];//GeV/c
    Double_t py = HBT.plast[i][1];//GeV/c
    Double_t pz = HBT.plast[i][2];//GeV/c
    Double_t ma = HBT.plast[i][3];//GeV/c/c
//    Double_t vx = 0;//HBT.xlast[i][0]*1e-12;//mm
//    Double_t vy = 0;//HBT.xlast[i][1]*1e-12;//mm
//    Double_t vz = 0;//HBT.xlast[i][2]*1e-12;//mm
//    Double_t vt = 0;//HBT.xlast[i][3]*1e-12;//mm/c
    Double_t vx = HBT.xlast[i][0]*1e-12;//mm
    Double_t vy = HBT.xlast[i][1]*1e-12;//mm
    Double_t vz = HBT.xlast[i][2]*1e-12;//mm
    Double_t vt = HBT.xlast[i][3]*1e-12;//mm/c
    Int_t pdg   = invflv(HBT.lblast[i]);
    TParticle *p = new TParticle(pdg,
                                 status,
                                 -1,
                                 -1,
                                 -1,
                                 -1,
                                 px,
                                 py,
                                 pz,
                                 TMath::Sqrt(ma*ma+px*px+py*py+pz*pz),
                                 vx,
                                 vy,
                                 vz,
                                 vt);
    if((px==0)&&(py==0)) {
      if(pz<0)
        p->SetUniqueID(0);
      else 
        p->SetUniqueID(10);
    } else 
      p->SetUniqueID(999);
    fParticles->Add(p);
  }
  return fParticles;
}

//______________________________________________________________________________
Int_t TAmpt::ImportParticles(TClonesArray *particles, Option_t */*option*/)
{
  // Import created particles.

  if (particles == 0) 
    return 0;

  TClonesArray &particlesR = *particles;
  particlesR.Clear();

  Int_t numpart = HBT.nlast;
  printf("TAmpt: AMPT stack contains %d particles.\n", numpart);

  //at this point not clear how to read particle history, just take primaries.
  for (Int_t i=0; i < numpart; ++i) {
    Int_t status = 1;
    Double_t px = HBT.plast[i][0];//GeV/c
    Double_t py = HBT.plast[i][1];//GeV/c
    Double_t pz = HBT.plast[i][2];//GeV/c
    Double_t ma = HBT.plast[i][3];//GeV/c/c
//    Double_t vx = 0;//HBT.xlast[i][0]*1e-12;//mm
//    Double_t vy = 0;//HBT.xlast[i][1]*1e-12;//mm
//    Double_t vz = 0;//HBT.xlast[i][2]*1e-12;//mm
//    Double_t vt = 0;//HBT.xlast[i][3]*1e-12;//mm/c
    Double_t vx = HBT.xlast[i][0]*1e-12;//mm
    Double_t vy = HBT.xlast[i][1]*1e-12;//mm
    Double_t vz = HBT.xlast[i][2]*1e-12;//mm
    Double_t vt = HBT.xlast[i][3]*1e-12;//mm/c
    Int_t pdg  = invflv(HBT.lblast[i]);
    //printf("i %d pdg %d px %f py %f pz %f vx %f vy %f vz %f vt %f\n", i, pdg, px, py, pz, vx, vy, vz, vt);
    new(particlesR[i]) TParticle(pdg,
                                 status,
                                 -1,
                                 -1,
                                 -1,
                                 -1,
                                 px,
                                 py,
                                 pz,
                                 TMath::Sqrt(ma*ma+px*px+py*py+pz*pz),
                                 vx,
                                 vy,
                                 vz,
                                 vt);
    if((px==0)&&(py==0)){
      if(pz<0)
        particlesR[i]->SetUniqueID(0);
      else 
        particlesR[i]->SetUniqueID(10);
    } else 
      particlesR[i]->SetUniqueID(999);
  }
  return numpart;
}

//______________________________________________________________________________
Int_t TAmpt::ImportNucleons(TClonesArray *nucleons, Option_t */*option*/)
{
  // Import created particles.

  if (nucleons == 0) 
    return 0;

  TClonesArray &nucleonsR = *nucleons;
  nucleonsR.Clear();

  Int_t nA = HPARNT.ihnt2[0];
  for (Int_t i=0; i < nA; ++i) {
    Double_t x = HJCRDN.yp[i][0] + 0.5*GetBB();
    Double_t y = HJCRDN.yp[i][1];
    Double_t z = HJCRDN.yp[i][2];
    Int_t    p = HSTRNG.nfp[3][i];
    Int_t    s = HSTRNG.nfp[4][i];
    new(nucleonsR[i]) TParticle(p,
                                s,
                                -1,
                                -1,
                                -1,
                                -1,
                                 0,
                                 0,
                                 0,
                                 0,
                                 x,
                                 y,
                                 z,
                                 0);
    nucleonsR[i]->SetUniqueID(1);
  }
  Int_t nB = HPARNT.ihnt2[2];
  for (Int_t i=0; i < nB; ++i) {
    Double_t x = HJCRDN.yt[i][0] - 0.5*HPARNT.hint1[18];
    Double_t y = HJCRDN.yt[i][1];
    Double_t z = HJCRDN.yt[i][2];
    Int_t    p = HSTRNG.nft[3][i];
    Int_t    s = HSTRNG.nft[4][i];
    new(nucleonsR[nA+i]) TParticle(p,
                                   s,
                                   -1,
                                   -1,
                                   -1,
                                   -1,
                                   0,
                                   0,
                                   0,
                                   0,
                                   x,
                                   y,
                                   z,
                                   0);
    nucleonsR[nA+i]->SetUniqueID(-1);
  }
  return nA+nB;
}

//______________________________________________________________________________
void TAmpt::SetEFRM(Float_t efrm)
{
  // Set the centre of mass (CMS) or lab-energy (LAB)
  fEfrm=efrm;
} 

//______________________________________________________________________________
void TAmpt::SetFRAME(const char* frame)
{
  // Set the frame type ("CMS" or "LAB")
  fFrame=frame;
} 

//______________________________________________________________________________
void TAmpt::SetPROJ(const char* proj)
{
  // Set the projectile type
  fProj=proj;
} 

//______________________________________________________________________________
void TAmpt::SetTARG(const char* targ)
{
  // Set the target type
  fTarg=targ;
} 

//______________________________________________________________________________
void TAmpt::SetIAP(Int_t iap)
{
  // Set the projectile atomic number
  fIap=iap;
} 

//______________________________________________________________________________
void TAmpt::SetIZP(Int_t izp)
{
  // Set the projectile charge number
  fIzp=izp;
} 

//______________________________________________________________________________
void TAmpt::SetIAT(Int_t iat)
{
  // Set the target atomic number
  fIat=iat;
} 

//______________________________________________________________________________
void TAmpt::SetIZT(Int_t izt)
{
  // Set the target charge number
  fIzt=izt;
} 

//______________________________________________________________________________
void TAmpt::SetBMIN(Float_t bmin)
{
  // Set the minimum impact parameter
  fBmin=bmin;
} 

//______________________________________________________________________________
void TAmpt::SetBMAX(Float_t bmax)
{
  // Set the maximum impact parameter
  fBmax=bmax;
} 

//______________________________________________________________________________
Float_t TAmpt::GetEFRM() const
{
  // Get the centre of mass (CMS) or lab-energy (LAB)
  return fEfrm;
} 

//______________________________________________________________________________
const char* TAmpt::GetFRAME() const
{
  // Get the frame type ("CMS" or "LAB")
  return fFrame.Data();
} 

//______________________________________________________________________________
const char* TAmpt::GetPROJ() const
{
  // Get the projectile type
  return fProj;
} 

//______________________________________________________________________________
const char* TAmpt::GetTARG() const
{
  // Set the target type
  return fTarg;
} 

//______________________________________________________________________________
Int_t TAmpt::GetIAP() const
{
  // Get the projectile atomic number
  return fIap;
} 

//______________________________________________________________________________
Int_t TAmpt::GetIZP() const
{
  // Get the projectile charge number
  return fIzp;
} 

//______________________________________________________________________________
Int_t TAmpt::GetIAT() const
{
  // Get the target atomic number
  return fIat;
} 

//______________________________________________________________________________
Int_t TAmpt::GetIZT() const
{
  // Get the target charge number
  return fIzt;
} 

//______________________________________________________________________________
Float_t TAmpt::GetBMIN() const
{
  // Get the minimum impact parameter
  return fBmin;
} 

//______________________________________________________________________________
Float_t TAmpt::GetBMAX() const
{
  // Get the maximum impact parameter
  return fBmax;
} 

//====================== access to common HIPARNT ===============================

//______________________________________________________________________________
void TAmpt::SetHIPR1(Int_t key,Float_t value)
{
  // Set the values of array HIPR1 in common HIPARNT
  if ( key<1 || key>100 ) {
    printf ("ERROR in TAmpt:SetHIPR1(key,value): \n ");
    printf ("      key=%i is out of range [1..100]!\n",key);
    return;
  }
  HPARNT.hipr1[key-1]=value;
}

//______________________________________________________________________________
Float_t TAmpt::GetHIPR1(Int_t key) const
{
  // Get the values of array HIPR1 in common HIPARNT
  if ( key<1 || key>100 ) {
    printf ("ERROR in TAmpt:GetHIPR1(key): \n ");
    printf ("      key=%i is out of range [1..100]!\n",key);
    return 0;
  }
  return HPARNT.hipr1[key-1];
}

//______________________________________________________________________________
void TAmpt::SetIHPR2(Int_t key,Int_t value)
{
  // Set the values of array HIPR2 in common HIPARNT
  if ( key<1 || key>50 ) {
    printf ("ERROR in TAmpt:SetIHPR2(key,value): \n ");
    printf ("      key=%i is out of range [1..50]!\n",key);
    return;
  }
  HPARNT.ihpr2[key-1]=value;
}

//______________________________________________________________________________
Int_t TAmpt::GetIHPR2(Int_t key) const
{
  // Get the values of array HIPR2 in common HIPARNT
  if ( key<1 || key>50 ) {
    printf ("ERROR in TAmpt:GetIHPR2(key): \n ");
    printf ("      key=%i is out of range [1..50]!\n",key);
    return 0;
  }
  return HPARNT.ihpr2[key-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetHINT1(Int_t key) const
{
  // Get the values of array HINT1 in common HIPARNT
  if ( key<1 || key>100 ) {
    printf ("ERROR in TAmpt:GetHINT1(key): \n ");
    printf ("      key=%i is out of range [1..100]!\n",key);
    return 0;
  }
  return HPARNT.hint1[key-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetIHNT2(Int_t key) const
{
  // Get the values of array HINT2 in common HIPARNT
  if ( key<1 || key>50 ) {
    printf ("ERROR in TAmpt:GetIHNT2(key): \n ");
    printf ("      key=%i is out of range [1..50]!\n",key);
    return 0;
  }
  return HPARNT.ihnt2[key-1];
}

//====================== access to common HMAIN1 ===============================

//______________________________________________________________________________
Int_t  TAmpt::GetNATT() const
{
  // Get the number of particles produces
  return HMAIN1.natt;
}

//______________________________________________________________________________
Float_t  TAmpt::GetEATT() const
{
  // Get total energy of particles
  return HMAIN1.eatt;
}

//______________________________________________________________________________
Int_t  TAmpt::GetJATT() const
{
  // Get number of hard scatterings
  return HMAIN1.jatt;
}

//______________________________________________________________________________
Int_t  TAmpt::GetNT() const
{
  // Get number of target participants
  return HMAIN1.nt;
}

//______________________________________________________________________________
Int_t  TAmpt::GetNP() const
{
  // Get number of projectile participants
  return HMAIN1.np;
}

//______________________________________________________________________________
Int_t  TAmpt::GetN0() const
{
  // Get number of N-N collisions
  return HMAIN1.n0;
}

//______________________________________________________________________________
Int_t  TAmpt::GetN01() const
{
  // Get number of N-wounded collisions
  return HMAIN1.n01;
}

//______________________________________________________________________________
Int_t  TAmpt::GetN10() const
{
  // Get number of wounded-N collisions
  return HMAIN1.n10;
}

//______________________________________________________________________________
Int_t  TAmpt::GetN11() const
{
  // Get number of wounded-wounded collisions
  return HMAIN1.n11;
}

//______________________________________________________________________________
Float_t  TAmpt::GetBB() const
{
  // Get impact parameter

  return HPARNT.hint1[18];
}

//====================== access to common HMAIN2 ===============================

//______________________________________________________________________________
Int_t TAmpt::GetKATT(Int_t key1, Int_t key2) const
{
  // Get values of array KATT in common HMAIN2
  if ( key1<1 || key1>200000 ) {
    printf("ERROR in TAmpt::GetKATT(key1,key2):\n");
    printf("      key1=%i is out of range [1..200000]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>4 ) {
    printf("ERROR in TAmpt::GetKATT(key1,key2):\n");
    printf("      key2=%i is out of range [1..4]\n",key2);
    return 0;
  }
  return HMAIN2.katt[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPATT(Int_t key1, Int_t key2) const
{
  // Get values of array PATT in common HMAIN2
  if ( key1<1 || key1>200000 ) {
    printf("ERROR in TAmpt::GetPATT(key1,key2):\n");
    printf("      key1=%i is out of range [1..130000]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>4 ) {
    printf("ERROR in TAmpt::GetPATT(key1,key2):\n");
    printf("      key2=%i is out of range [1..4]\n",key2);
    return 0;
  }
  return HMAIN2.patt[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetVATT(Int_t key1, Int_t key2) const
{
  // Get values of array VATT in common HMAIN2
  if ( key1<1 || key1>200000 ) {
    printf("ERROR in TAmpt::GetVATT(key1,key2):\n");
    printf("      key1=%i is out of range [1..130000]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>4 ) {
    printf("ERROR in TAmpt::GetVATT(key1,key2):\n");
    printf("      key2=%i is out of range [1..4]\n",key2);
    return 0;
  }
  return HMAIN2.vatt[key2-1][key1-1];
}

//====================== access to common HIJJET1 ===============================

//______________________________________________________________________________
Int_t TAmpt::GetNPJ(Int_t key) const
{
  // Get values of array NPJ of common HIJJET1
  if ( key<1 || key>300 ) {
    printf("ERROR in TAmpt::GetNPJ(key):\n");
    printf("      key=%i is out of range [1..300]\n",key);
    return 0;
  }
  return HJJET1.npj[key-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetKFPJ(Int_t key1, Int_t key2) const
{
  // Get values of array KFPJ in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetKFPJ(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetKFPJ(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.kfpj[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJPX(Int_t key1, Int_t key2) const
{
  // Get values of array PJPX in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJPX(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJPX(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjpx[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJPY(Int_t key1, Int_t key2) const
{
  // Get values of array PJPY in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJPY(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJPY(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjpy[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJPZ(Int_t key1, Int_t key2) const
{
  // Get values of array PJPZ in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJPZ(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJPZ(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjpz[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJPE(Int_t key1, Int_t key2) const
{
  // Get values of array PJPE in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJPE(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJPE(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjpe[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJPM(Int_t key1, Int_t key2) const
{
  // Get values of array PJPM in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJPM(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJPM(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjpm[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetNTJ(Int_t key) const
{
  // Get values of array NTJ in common HIJJET1
  if ( key<1 || key>300 ) {
    printf("ERROR in TAmpt::GetNTJ(key):\n");
    printf("      key=%i is out of range [1..300]\n",key);
    return 0;
  }
  return HJJET1.ntj[key-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetKFTJ(Int_t key1, Int_t key2) const
{
  // Get values of array KFTJ in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetKFTJ(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetKFTJ(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.kftj[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJTX(Int_t key1, Int_t key2) const
{
  // Get values of array PJTX in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJTX(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJTX(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjtx[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJTY(Int_t key1, Int_t key2) const
{
  // Get values of array PJTY in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJTY(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJTY(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjty[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJTZ(Int_t key1, Int_t key2) const
{
  // Get values of array PJTZ in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJTZ(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJTZ(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjtz[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJTE(Int_t key1, Int_t key2) const
{
  // Get values of array PJTE in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJTE(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJTE(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjte[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPJTM(Int_t key1, Int_t key2) const
{
  // Get values of array PJTM in common HIJJET1
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPJTM(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>500 ) {
    printf("ERROR in TAmpt::GetPJTM(key1,key2):\n");
    printf("      key2=%i is out of range [1..500]\n",key2);
    return 0;
  }
  return HJJET1.pjtm[key2-1][key1-1];
}

//====================== access to common HIJJET1 ===============================

//______________________________________________________________________________
Int_t TAmpt::GetNSG() const
{
  // Get value of NSG in common HIJJET2
  return HJJET2.nsg;
}

//______________________________________________________________________________
Int_t TAmpt::GetNJSG(Int_t key) const
{
  // Get values of array NJSG in common HIJJET2
  if ( key<1 || key>900 ) {
    printf ("ERROR in TAmpt:GetNJSG(key): \n ");
    printf ("      key=%i is out of range [1..900]!\n",key);
    return 0;
  }
  return HJJET2.njsg[key-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetIASG(Int_t key1, Int_t key2) const
{
  // Get values of IASG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetIASG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>3 ) {
    printf("ERROR in TAmpt::GetIASG(key1,key2):\n");
    printf("      key2=%i is out of range [1..3]\n",key2);
    return 0;
  }
  return HJJET2.iasg[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetK1SG(Int_t key1, Int_t key2) const
{
  // Get values of K1SG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetK1SG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetK1SG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.k1sg[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetK2SG(Int_t key1, Int_t key2) const
{
  // Get values of K2SG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetK2SG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetK2SG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.k2sg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPXSG(Int_t key1, Int_t key2) const
{
  // Get values of PXSG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetPXSG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetPXSG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.pxsg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPYSG(Int_t key1, Int_t key2) const
{
  // Get values of PYSG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetPYSG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetPYSG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.pysg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPZSG(Int_t key1, Int_t key2) const
{
  // Get values of PZSG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetPZSG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetPZSG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.pzsg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPESG(Int_t key1, Int_t key2) const
{
  // Get values of PESG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetPESG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetPESG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.pesg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPMSG(Int_t key1, Int_t key2) const
{
  // Get values of PMSG in common HIJJET2
  if ( key1<1 || key1>900 ) {
    printf("ERROR in TAmpt::GetPMSG(key1):\n");
    printf("      key1=%i is out of range [1..900]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>100 ) {
    printf("ERROR in TAmpt::GetPMSG(key1,key2):\n");
    printf("      key2=%i is out of range [1..100]\n",key2);
    return 0;
  }
  return HJJET2.pmsg[key2-1][key1-1];
}

//====================== access to common HISTRNG ===============================

//______________________________________________________________________________
Int_t TAmpt::GetNFP(Int_t key1, Int_t key2) const
{
  // Get values of array NFP in common HISTRNG
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetNFP(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>15 ) {
    printf("ERROR in TAmpt::GetNFP(key1,key2):\n");
    printf("      key2=%i is out of range [1..15]\n",key2);
    return 0;
  }
  return HSTRNG.nfp[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPP(Int_t key1, Int_t key2) const
{
  // Get values of array PP in common HISTRNG
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPP(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>15 ) {
    printf("ERROR in TAmpt::GetPP(key1,key2):\n");
    printf("      key2=%i is out of range [1..15]\n",key2);
    return 0;
  }
  return HSTRNG.pp[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t TAmpt::GetNFT(Int_t key1, Int_t key2) const
{
  // Get values of array NFT in common HISTRNG
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetNFT(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>15 ) {
    printf("ERROR in TAmpt::GetNFT(key1,key2):\n");
    printf("      key2=%i is out of range [1..15]\n",key2);
    return 0;
  }
  return HSTRNG.nft[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t TAmpt::GetPT(Int_t key1, Int_t key2) const
{
  // Get values of array PT in common HISTRNG
  if ( key1<1 || key1>300 ) {
    printf("ERROR in TAmpt::GetPT(key1):\n");
    printf("      key1=%i is out of range [1..300]\n",key1);
    return 0;
  }
  if ( key2<1 || key2>15 ) {
    printf("ERROR in TAmpt::GetPT(key1,key2):\n");
    printf("      key2=%i is out of range [1..15]\n",key2);
    return 0;
  }
  return HSTRNG.pt[key2-1][key1-1];
}

void TAmpt::SetPARJ(Int_t key, Float_t parm) 
{
  // Set values of array PARJ in common HISTRNG
  if ( key < 1 || key > 200) {
    printf("ERROR in TAmpt::SetPARJ(key,parm):\n");
    printf("      key=%i is out of range [1..200]\n",key);
    return;
  }
  LUDAT1.parj[key-1] = parm;
}


void TAmpt::SetMSTJ(Int_t key, Int_t parm) 
{
  // Set values of array MSTJ in common HISTRNG
  if ( key < 1 || key > 200) {
    printf("ERROR in TAmpt::SetMSTJ(key,parm):\n");
    printf("      key=%i is out of range [1..200]\n",key);
    return;
  }
  LUDAT1.mstj[key-1] = parm;
}

//====================== access to Ampt subroutines =========================

//______________________________________________________________________________
void TAmpt::Initialize()
{
  // Calls Ampset with the either default parameters or the ones set by the user  //
  // via SetEFRM, SetFRAME, SetPROJ, SetTARG, SetIAP, SetIZP, SetIAT, SetIZT      //

  if ( (!strcmp(fFrame.Data(), "CMS     "  )) &&
       (!strcmp(fFrame.Data(), "LAB     "  ))){
    printf("WARNING! In TAmpt:Initialize():\n");
    printf(" specified frame=%s is neither CMS or LAB\n",fFrame.Data());
    printf(" resetting to default \"CMS\" .");
    fFrame="CMS";
  }

  if ( (!strcmp(fProj.Data(), "A       "     )) &&
       (!strcmp(fProj.Data(), "P       "     )) &&
       (!strcmp(fProj.Data(), "PBAR    "  ))){
    printf("WARNING! In TAmpt:Initialize():\n");
    printf(" specified projectile=%s is neither A, P or PBAR\n",fProj.Data());
    printf(" resetting to default \"A\" .");
    fProj="A";
  }

  if ( (!strcmp(fTarg.Data(), "A       "     )) &&
       (!strcmp(fTarg.Data(), "P       "     )) &&
       (!strcmp(fTarg.Data(), "PBAR    "  ))){
    printf("WARNING! In TAmpt:Initialize():\n");
    printf(" specified target=%s is neither A, P or PBAR\n",fTarg.Data());
    printf(" resetting to default \"A\" .");
    fTarg="A";
  }

  //printf(" %s-%s at %f GeV \n",fProj.Data(),fTarg.Data(),fEfrm);
  Amptset(fEfrm,fFrame.Data(),fProj.Data(),fTarg.Data(),fIap,fIzp,fIat,fIzt);
}

//______________________________________________________________________________
void TAmpt::GenerateEvent()
{
  // Generates one event;

  //printf("Next event ------------------------\n");
  Ampt(fFrame.Data(),fBmin,fBmax);
}

//______________________________________________________________________________
void TAmpt::Amptset(double efrm, const char *frame, const char *proj,  
                    const char *targ, int iap, int izp, int iat, int izt)
{
  // Call AMPT routine HIJSET passing the parameters in a way accepted by
  // Fortran routines				   

  int s1 = strlen(frame);
  int s2 = strlen(proj);
  int s3 = strlen(targ);
  //printf("s1 = %d s2 = %d s3 = %d\n",s1,s2,s3);
#ifndef WIN32 
  amptset(efrm, frame, proj, targ, iap, izp, iat, izt, s1, s2, s3);
#else
  amptset(efrm, frame, s1, proj, s2, targ, s3, iap, izp, iat, izt);
#endif
}

//______________________________________________________________________________
void TAmpt::Ampt(const char *frame, double bmin, double bmax)
{
  // Call AMPT routine HIJSET passing the parameters in a way accepted by
  // Fortran routines				   
  
  int s1 = strlen(frame);
  
#ifndef WIN32 
  ampt(frame, bmin, bmax, s1);
#else
  ampt(frame, s1, bmin, bmax);
#endif
}

Float_t TAmpt::Profile(float b)
{
  // Call AMPT routine PROFILE 
  return profile(b);
}

//______________________________________________________________________________
void TAmpt::Rluget(Int_t lfn, Int_t move)
{
  // write seed to file
  rluget_ampt(lfn, move);
}

//______________________________________________________________________________
void TAmpt::Rluset(Int_t lfn, Int_t move)
{
  // read seed from file 
  rluset_ampt(lfn, move);
}

//______________________________________________________________________________
void TAmpt::SetIsoft(Int_t i)     
{
  // set ampt mode.
  ANIM.isoft    = i;  
}

//______________________________________________________________________________
void TAmpt::SetNtMax(Int_t max)   
{
  // set maximum number of timesteps
  INPUT2.ntmax  = max;
}

//______________________________________________________________________________
void TAmpt::SetIpop(Int_t pop)    
{
  // set flag for popcorn mechanism(netbaryon stopping)
  POPCORN.ipop  = pop;
}

//______________________________________________________________________________
void TAmpt::SetXmu(Float_t m)     
{
  // set parton screening mass in fm^-1.
  PARA2.xmu     = m;  
}

//______________________________________________________________________________
void TAmpt::SetAlpha(Float_t alpha)     
{
  // set parton screening mass in fm^-1.
  PARA2.alpha     = alpha;  
}

//______________________________________________________________________________
void TAmpt::SetStringFrag(Float_t a, Float_t b)
{
  // Set string frag parameters, f(z) = 1/z*(1-z)^a*exp(-b*mz^2/z).
   SetPARJ(41, a);
   SetPARJ(42, b);
}
