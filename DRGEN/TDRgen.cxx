/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
// TDRgen is an interface class to fortran event generator of
// two-pomeron processes in HE PP collisions
//%
// Sergey.Evdokimov@cern.ch
// 9 Oct 2013
//------------------------------------------------------------------------

#include "TRandom.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "TDRgen.h"
#define IN_TDRGEN_CXX
#include "DRGENcommon.h"

# define pp2pxp_gen  pp2pxp_gen_
# define pp2rnd   pp2rnd_

extern "C" {
  void pp2pxp_gen();
  Float_t pp2rnd(Int_t*) {
    Float_t r;
    do r=gRandom->Rndm(); while(0 >= r || r >= 1);
    return r;
  }
}

ClassImp(TDRgen)

//------------------------------------------------------------------------------
TDRgen::TDRgen() : TGenerator("DRgen","DRgen")
{
  // DRgen constructor
  SetIPROC(2,2);
  SetSqrtS(7000.);
  Double_t procXsec[20]={3.0,1.5,1.0,0.3,1.0};
  SetProcXsec(procXsec);
  SetAMIN(0.280);
  SetAMAX(3.0);
}

//------------------------------------------------------------------------------
TDRgen::~TDRgen()
{
  // Destructor
}

//______________________________________________________________________________
void TDRgen::Initialize()
{
  // Initialize DRgen
  pp2pxp_gen();
}

//______________________________________________________________________________
void TDRgen::GenerateEvent()
{
  //produce one event
 pp2pxp_gen() ;
}

//______________________________________________________________________________
void TDRgen::Finish()
{
  //do nothing here  
}

//______________________________________________________________________________

// Setters
void TDRgen::SetIPROC    (Int_t   iproc1, Int_t iproc2   ) const
{
  PP2INIT.Iproc1    = iproc1;
  PP2INIT.Iproc2    = iproc2;

}
//______________________________________________________________________________
void TDRgen::SetSqrtS   (Double_t sqrts  ) const
{
  PP2INIT.sqrtS = sqrts;      
}
//______________________________________________________________________________
void TDRgen::SetAMIN    (Double_t amin   ) const
{
  PP2INIT.aMmin = amin;
}
//______________________________________________________________________________
void TDRgen::SetAMAX   (Double_t amax  ) const
{
  PP2INIT.aMmax = amax;
}
//______________________________________________________________________________
void TDRgen::SetProcXsec     (Double_t* procXsec    ) const
{
  for (Int_t i=0;i<20;i++)
    {
      PP2INIT.ProcXsec[i] = procXsec[i];
    }
}
//______________________________________________________________________________
void TDRgen::SetF2Polarization     (Double_t* polarization    ) const
{
  for (Int_t i=0;i<8;i++)
    {
      PP2INIT.F2Polarization[i] = polarization[i];
    }
}

//______________________________________________________________________________
// Getters for COMMON /PP2EVNT/
Int_t TDRgen::GetIevent() const
{
  return PP2EVNT.ievent;
}
//______________________________________________________________________________
Int_t TDRgen::GetIproc() const
{
  return PP2EVNT.Iproc;
}
//______________________________________________________________________________
Int_t TDRgen::GetNpart() const
{
  return PP2EVNT.Npart;
}
//______________________________________________________________________________
Int_t TDRgen::GetIParticleNumber() const
{
  return PP2EVNT.iParticleNumber;
}
//______________________________________________________________________________
void TDRgen::ImportParticles(TClonesArray* particles)
{

  TClonesArray &clonesParticles = *particles;
  clonesParticles.Clear();
  
  for (Int_t i =0;i<PP2EVNT.iParticleNumber&&i<10;i++)// 10 particles max. in common block now.
    {
      //	TParticle(Int_t pdg, Int_t status, Int_t mother1, Int_t mother2, Int_t daughter1, Int_t daughter2, Double_t px, Double_t py, Double_t pz, Double_t etot, Double_t vx, Double_t vy, Double_t vz, Double_t time)
      new(clonesParticles[i]) TParticle(
					PP2EVNT.iParticleCode[i],
					PP2EVNT.iParticleStatus[i],
					0,
					0,
					0,
					0,
					PP2EVNT.pParticle[i][0],
					PP2EVNT.pParticle[i][1],
					PP2EVNT.pParticle[i][2],
					PP2EVNT.pParticle[i][3],
					0,
					0,
					0,
					0);
      // clonesParticles.ConstructedAt(i)->Print();
    }
}
