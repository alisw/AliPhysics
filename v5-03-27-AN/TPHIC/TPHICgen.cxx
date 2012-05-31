/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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
 *                                                                        *
 **************************************************************************/
//------------------------------------------------------------------------
// TPHICgen is an interface class to fortran event generator of
// two-photon processes in ultra-peripheral ion collisions
//%
// Yuri.Kharlov@cern.ch
// 15 April 2003
//------------------------------------------------------------------------

#include "TRandom.h"

#include "TPHICgen.h"
//#include "TClonesArray.h"
//#include "TParticle.h"
#define IN_TPHICGEN_CXX
#include "TPHICcommon.h"

#ifndef WIN32
# define gginit  gginit_
# define ggrun   ggrun_
# define ggexit  ggexit_
# define ggrnd   ggrnd_
#else
# define gginit  gginit
# define ggrun   ggrun
# define ggexit  ggexit
# define ggrnd   ggrnd
#endif

extern "C" {
  void gginit();
  void ggrun() ;
  void ggexit();
//   Double_t pyr_(Int_t*);
  Double_t ggrnd(Int_t*) {
    Double_t r;
    do r=gRandom->Rndm(); while(0 >= r || r >= 1);
    return r;
//     return pyr_(0);
  }
}

ClassImp(TPHICgen)

//------------------------------------------------------------------------------
TPHICgen::TPHICgen() : TGenerator("TPHIC","TPHIC")
{
  // TPHIC constructor
}

//------------------------------------------------------------------------------
TPHICgen::~TPHICgen()
{
  // Destructor
}

//______________________________________________________________________________
void TPHICgen::Initialize()
{
  // Initialize TPHIC
  const Float_t kNucleonMass = 0.931494;
  SetAMAS(GGINI.ia*kNucleonMass);
  gginit();
}

//______________________________________________________________________________
void TPHICgen::GenerateEvent()
{
  //produce one event
  ggrun();
}

//______________________________________________________________________________
void TPHICgen::Finish()
{
  // calculate cross section and print out cross section and related variables
  ggexit();
}

//______________________________________________________________________________

// Setters
void TPHICgen::SetIPROC    (Int_t   iproc   ) const
{
  GGINI.iproc    = iproc;        
}
//______________________________________________________________________________
void TPHICgen::SetNEVENT   (Int_t   nevent  ) const
{
  GGINI.nevent   = nevent;        
}
//______________________________________________________________________________
void TPHICgen::SetILUMF    (Int_t   ilumf   ) const
{
  GGINI.ilumf    = ilumf;        
}
//______________________________________________________________________________
void TPHICgen::SetLUMFIL   (TString lumfil  ) const
{
  for (Int_t i=0; i<lumfil.Length(); i++)
    GGINI.lumfil[i]   = lumfil[i];
}
//______________________________________________________________________________
void TPHICgen::SetEBMN     (Float_t ebmn    ) const
{
  GGINI.ebmn     = ebmn;         
}
//______________________________________________________________________________
void TPHICgen::SetIZ       (Int_t   iz      ) const
{
  GGINI.iz       = iz;           
}
//______________________________________________________________________________
void TPHICgen::SetIA       (Int_t   ia      ) const
{
  GGINI.ia       = ia;           
}
//______________________________________________________________________________
void TPHICgen::SetAMAS     (Float_t amas    ) const
{
  GGINI.amas     = amas;         
}
//______________________________________________________________________________
void TPHICgen::SetAMIN     (Float_t amin    ) const
{
  GGINI.amin     = amin;         
}
//______________________________________________________________________________
void TPHICgen::SetAMAX     (Float_t amax    ) const
{
  GGINI.amax     = amax;         
}
//______________________________________________________________________________
void TPHICgen::SetYMIN     (Float_t ymin    ) const
{
  GGINI.ymin     = ymin;         
}
//______________________________________________________________________________
void TPHICgen::SetYMAX     (Float_t ymax    ) const
{
  GGINI.ymax     = ymax;         
}
//______________________________________________________________________________
void TPHICgen::SetNMAS     (Int_t   nmas    ) const
{
  GGINI.nmas     = nmas;         
}
//______________________________________________________________________________
void TPHICgen::SetNY       (Int_t   ny      ) const
{
  GGINI.ny       = ny;           
}
//______________________________________________________________________________
void TPHICgen::SetKFERM    (Int_t   kferm   ) const
{
  GGINI.kferm    = kferm;        
}
//______________________________________________________________________________
void TPHICgen::SetKFONIUM  (Int_t   kfonium ) const
{
  GGINI.kfonium  = kfonium;      
}
//______________________________________________________________________________
void TPHICgen::SetXMRES    (Float_t xmres   ) const
{
  GGINI.xmres    = xmres;        
}
//______________________________________________________________________________
void TPHICgen::SetXGTRES   (Float_t xgtres  ) const
{
  GGINI.xgtres   = xgtres;       
}
//______________________________________________________________________________
void TPHICgen::SetXGGRES   (Float_t xggres  ) const
{
  GGINI.xggres   = xggres;       
}
//______________________________________________________________________________
void TPHICgen::SetMODDCY   (Int_t   moddcy  ) const
{
  GGINI.moddcy   = moddcy;       
}
//______________________________________________________________________________
void TPHICgen::SetTHETAMIN (Float_t thetamin) const
{
  GGINI.thetamin = thetamin;     
}
//______________________________________________________________________________
void TPHICgen::SetKV1      (Int_t   kv1     ) const
{
  GGINI.kv1      = kv1;
}
//______________________________________________________________________________
void TPHICgen::SetKV2      (Int_t   kv2     ) const
{
  GGINI.kv2      = kv2;
}

//______________________________________________________________________________
// Getters for COMMON /GGEVNT/
Float_t TPHICgen::GetWSQ() const
{
  return GGEVNT.wsq;
}
//______________________________________________________________________________
Float_t TPHICgen::GetYGG() const
{
  return GGEVNT.ygg;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXMG1() const
{
  return GGEVNT.xmg1;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXMG2() const
{
  return GGEVNT.xmg2;
}
//______________________________________________________________________________
Float_t TPHICgen::GetP2G(Int_t i) const
{
  return GGEVNT.p2g[i];
}
//______________________________________________________________________________
Float_t TPHICgen::GetPTAG1(Int_t i) const
{
  return GGEVNT.ptag1[i];
}
//______________________________________________________________________________
Float_t TPHICgen::GetPTAG2(Int_t i) const
{
  return GGEVNT.ptag2[i];
}
//______________________________________________________________________________
Int_t   TPHICgen::GetNGG() const
{
  return GGEVNT.ngg;
}
//______________________________________________________________________________
Int_t   TPHICgen::GetKGG(Int_t i) const
{
  return GGEVNT.kgg[i];
}
//______________________________________________________________________________
Float_t TPHICgen::GetPGG(Int_t i, Int_t j) const
{
  return GGEVNT.pgg[i][j];
}

//______________________________________________________________________________
// Getters for COMMON /GGXS/
Float_t TPHICgen::GetXSMAX0() const
{
  return GGXS.xsmax0;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSCUR0() const
{
  return GGXS.xscur0;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSCUR () const
{
  return GGXS.xscur;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSTOT () const
{
  return GGXS.xstot;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSTOTE() const
{
  return GGXS.xstote;
}
