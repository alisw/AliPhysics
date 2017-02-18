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

//-----------------------------------------------------------------------------
//   Class: AliDpmJetRndm
//   Responsibilities: Interface to Root random number generator 
//                     from Fortran (re-implements FINCTION dt_rndm_dpmjet)
//                     Very similar to AliHijingRndm
//   Note: Since AliGenDpmJet belongs to another module (TDPMjet) one cannot
//         pass the ponter to the generator via static variable
//   Collaborators: AliGenDPMjet class
//-----------------------------------------------------------------------------

#include <TRandom.h>

#include "AliDpmJetRndm.h"

TRandom * AliDpmJetRndm::fgDpmJetRandom=0;

ClassImp(AliDpmJetRndm)


//_______________________________________________________________________
void AliDpmJetRndm::SetDpmJetRandom(TRandom *ran) {
  //
  // Sets the pointer to an existing random numbers generator
  //
  if(ran) fgDpmJetRandom=ran;
  else fgDpmJetRandom=gRandom;
}

//_______________________________________________________________________
TRandom * AliDpmJetRndm::GetDpmJetRandom() {
  //
  // Retrieves the pointer to the random numbers generator
  //
  return fgDpmJetRandom;
}

#ifndef WIN32
# define dt_rndm_dpmjet   dt_rndm_dpmjet_
# define dt_rndmst_dpmjet dt_rndmst_dpmjet_
# define dt_rndmin_dpmjet dt_rndmin_dpmjet_
# define dt_rndmou_dpmjet dt_rndmou_dpmjet_
# define rninit_dpmjet 	  rninit_dpmjet_
# define type_of_call
#else
# define dt_rndm_dpmjet   DT_RNDM_DPMJET_
# define dt_rndmst_dpmjet DT_RNDMST_DPMJET
# define dt_rndmin_dpmjet DT_RNDMIN_DPMJET
# define dt_rndmou_dpmjet DT_RNDMOU_DPMJET
# define rninit_dpmjet    RNINIT_DPMJET
# define type_of_call _stdcall
#endif


extern "C" {
  void type_of_call dt_rndmst_(Int_t &, Int_t &, Int_t &, Int_t &)
  {printf("Dummy version of dt_rndmst reached\n");}

  void type_of_call dt_rndmin_(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &, Int_t &)
  {printf("Dummy version of dt_rndmin reached\n");}

  void type_of_call dt_rndmou_(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &, Int_t &)
  {printf("Dummy version of dt_rndmou reached\n");}

  void type_of_call dt_rndmte_(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &, Int_t &)
  {printf("Dummy version of dt_rndmou reached\n");}

  void type_of_call rninit_(Int_t &, Int_t &, Int_t &, Int_t &)
  {printf("Dummy version of rninit reached\n");}

  Double_t type_of_call dt_rndm_(Int_t &) 
  {
    // Wrapper to static method which retrieves the 
    // pointer to the Root (C++) generator
      Float_t r;
      do r = AliDpmJetRndm::GetDpmJetRandom()->Rndm();
      while(0 >= r || r >= 1);
      return r;
  }
}


