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
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Particle generator for the simplescopic TRD simulator                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <stdlib.h>
 
#include <TRandom.h>
#include <TMCProcess.h>

#include "AliRun.h"

#include "AliTRDsimpleGen.h"
#include "AliTRDsimpleMC.h"
#include "AliMC.h"
 
ClassImp(AliTRDsimpleGen)
 
//_____________________________________________________________________________
AliTRDsimpleGen::AliTRDsimpleGen():TObject()
{                       
  //
  // AliTRDsimpleGen default constructor
  //

  fPdg    = 211;
  fMomMin = 1.0;
  fMomMax = 1.0;             
                                            
}                                                                               
 
//_____________________________________________________________________________
AliTRDsimpleGen::AliTRDsimpleGen(const AliTRDsimpleGen &g):TObject(g)
{
  //
  // AliTRDsimpleGen copy constructor
  //
 
  ((AliTRDsimpleGen &) g).Copy(*this);
 
}
 
//_____________________________________________________________________________
AliTRDsimpleGen::~AliTRDsimpleGen()
{
  //
  // AliTRDsimpleGen destructor
  //
 
}                                                                               
 
//_____________________________________________________________________________
AliTRDsimpleGen &AliTRDsimpleGen::operator=(const AliTRDsimpleGen &g)
{
  //
  // Assignment operator
  //
 
  if (this != &g) ((AliTRDsimpleGen &) g).Copy(*this);
  return *this;
 
}
 
//_____________________________________________________________________________
void AliTRDsimpleGen::Copy(TObject &g) const
{
  //
  // Copy function
  //                             

  ((AliTRDsimpleGen &) g).fPdg     = fPdg;                                                 
  ((AliTRDsimpleGen &) g).fMomMin  = fMomMin;                                                 
  ((AliTRDsimpleGen &) g).fMomMax  = fMomMax;                                                 

}

//_____________________________________________________________________________
void AliTRDsimpleGen::NewParticle(Int_t ievent)
{
  //
  // Generate a new particle and initialize the MC object
  // 

  if (ievent == 0) {
    printf("\n");
    printf("<AliTRDsimpleGen> Generate particles with PDG code %d\n",fPdg);
    if (fMomMax > fMomMin) {
      printf("<AliTRDsimpleGen> Momentum range = %4.2f - %4.2f GeV/c\n"
            ,fMomMin,fMomMax);
    }
    else {
      printf("<AliTRDsimpleGen> Fixed momentum = %4.2f GeV/c\n"
            ,fMomMax);
    }
    printf("\n");

    // Add one dummy particle to the stack so that AddHit will work
    Float_t mom[3] = { 0.0 };
    Float_t vtx[3] = { 0.0 };
    Float_t pol[3] = { 0.0 };
    Int_t   ntr    = 0;
    gAlice->GetMCApp()->PushTrack(0,-1,fPdg,mom,vtx,pol,0.0,kPPrimary,ntr);

  }

  Double_t p = fMomMax;
  if (fMomMax > fMomMin) {
    p = (fMomMax - fMomMin) * gRandom->Rndm() + fMomMin;
  }

  Double_t px = p;
  Double_t py = 0.0;      
  Double_t pz = 0.0;      

  ((AliTRDsimpleMC *) gMC)->NewTrack(ievent,fPdg,px,py,pz);

}                                   
