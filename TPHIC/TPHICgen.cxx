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
#include "TClonesArray.h"
#include "TParticle.h"
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
void TPHICgen::SetIPROC    (const Int_t   iproc   )
{
  GGINI.iproc    = iproc;        
}
//______________________________________________________________________________
void TPHICgen::SetNEVENT   (const Int_t   nevent  )
{
  GGINI.nevent   = nevent;        
}
//______________________________________________________________________________
void TPHICgen::SetILUMF    (const Int_t   ilumf   )
{
  GGINI.ilumf    = ilumf;        
}
//______________________________________________________________________________
void TPHICgen::SetLUMFIL   (const TString lumfil  )
{
  for (Int_t i=0; i<lumfil.Length(); i++)
    GGINI.lumfil[i]   = lumfil[i];
}
//______________________________________________________________________________
void TPHICgen::SetEBMN     (const Float_t ebmn    )
{
  GGINI.ebmn     = ebmn;         
}
//______________________________________________________________________________
void TPHICgen::SetIZ       (const Int_t   iz      )
{
  GGINI.iz       = iz;           
}
//______________________________________________________________________________
void TPHICgen::SetIA       (const Int_t   ia      )
{
  GGINI.ia       = ia;           
}
//______________________________________________________________________________
void TPHICgen::SetAMAS     (const Float_t amas    )
{
  GGINI.amas     = amas;         
}
//______________________________________________________________________________
void TPHICgen::SetAMIN     (const Float_t amin    )
{
  GGINI.amin     = amin;         
}
//______________________________________________________________________________
void TPHICgen::SetAMAX     (const Float_t amax    )
{
  GGINI.amax     = amax;         
}
//______________________________________________________________________________
void TPHICgen::SetYMIN     (const Float_t ymin    )
{
  GGINI.ymin     = ymin;         
}
//______________________________________________________________________________
void TPHICgen::SetYMAX     (const Float_t ymax    )
{
  GGINI.ymax     = ymax;         
}
//______________________________________________________________________________
void TPHICgen::SetNMAS     (const Int_t   nmas    )
{
  GGINI.nmas     = nmas;         
}
//______________________________________________________________________________
void TPHICgen::SetNY       (const Int_t   ny      )
{
  GGINI.ny       = ny;           
}
//______________________________________________________________________________
void TPHICgen::SetKFERM    (const Int_t   kferm   )
{
  GGINI.kferm    = kferm;        
}
//______________________________________________________________________________
void TPHICgen::SetKFONIUM  (const Int_t   kfonium )
{
  GGINI.kfonium  = kfonium;      
}
//______________________________________________________________________________
void TPHICgen::SetXMRES    (const Float_t xmres   )
{
  GGINI.xmres    = xmres;        
}
//______________________________________________________________________________
void TPHICgen::SetXGTRES   (const Float_t xgtres  )
{
  GGINI.xgtres   = xgtres;       
}
//______________________________________________________________________________
void TPHICgen::SetXGGRES   (const Float_t xggres  )
{
  GGINI.xggres   = xggres;       
}
//______________________________________________________________________________
void TPHICgen::SetMODDCY   (const Int_t   moddcy  )
{
  GGINI.moddcy   = moddcy;       
}
//______________________________________________________________________________
void TPHICgen::SetTHETAMIN (const Float_t thetamin)
{
  GGINI.thetamin = thetamin;     
}
//______________________________________________________________________________
void TPHICgen::SetKV1      (const Int_t   kv1     )
{
  GGINI.kv1      = kv1;
}
//______________________________________________________________________________
void TPHICgen::SetKV2      (const Int_t   kv2     )
{
  GGINI.kv2      = kv2;
}

//______________________________________________________________________________
// Getters for COMMON /GGEVNT/
Float_t TPHICgen::GetWSQ()
{
  return GGEVNT.wsq;
}
//______________________________________________________________________________
Float_t TPHICgen::GetYGG()
{
  return GGEVNT.ygg;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXMG1()
{
  return GGEVNT.xmg1;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXMG2()
{
  return GGEVNT.xmg2;
}
//______________________________________________________________________________
Float_t TPHICgen::GetP2G(const Int_t i)
{
  return GGEVNT.p2g[i];
}
//______________________________________________________________________________
Float_t TPHICgen::GetPTAG1(const Int_t i)
{
  return GGEVNT.ptag1[i];
}
//______________________________________________________________________________
Float_t TPHICgen::GetPTAG2(const Int_t i)
{
  return GGEVNT.ptag2[i];
}
//______________________________________________________________________________
Int_t   TPHICgen::GetNGG()
{
  return GGEVNT.ngg;
}
//______________________________________________________________________________
Int_t   TPHICgen::GetKGG(const Int_t i)
{
  return GGEVNT.kgg[i];
}
//______________________________________________________________________________
Float_t TPHICgen::GetPGG(const Int_t i, const Int_t j)
{
  return GGEVNT.pgg[i][j];
}

//______________________________________________________________________________
// Getters for COMMON /GGXS/
Float_t TPHICgen::GetXSMAX0()
{
  return GGXS.xsmax0;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSCUR0()
{
  return GGXS.xscur0;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSCUR ()
{
  return GGXS.xscur;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSTOT ()
{
  return GGXS.xstot;
}
//______________________________________________________________________________
Float_t TPHICgen::GetXSTOTE()
{
  return GGXS.xstote;
}
