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
// 
// Author: Mikael.Mieskolainen@cern.ch


/* $Id:*/

#include <TROOT.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include "DCommon.h"
#include "TDime.h"
#include "AliDimeRndm.h"

#ifndef WIN32
# define dimeinit dimeinit_
# define dimegenerate dimegenerate_
# define type_of_call
#else
# define dimeinit DIMEINIT
# define dimegenerate DIMEGENERATE
# define type_of_call _stdcall
#endif

#ifndef WIN32
extern "C" void type_of_call dimeinit();
extern "C" void type_of_call dimegenerate(int& success);
#else
#endif



ClassImp(TDime)


TDime::TDime():
    TGenerator("Dime","Dime"),
    fEfrm(7000.0),
    fProcess("pipm      "),
    fFormf("orexp     "),
    fFsi("true      "),
    fEcut(0.0),
    fRmin(-2.0),
    fRmax(2.0)
{
// Default constructor

// Setup random number generator (if user doesn't set it, doesn't work without)
AliDimeRndm::SetDimeRandom(new TRandom3(123456789)); // default seed 123456789
}

//______________________________________________________________________________
TDime::TDime(Double_t efrm) :
    TGenerator("Dime","Dime"),
    fEfrm(efrm),
    fProcess("pipm      "),
    fFormf("orexp     "),
    fFsi("true      "),
    fEcut(0.0),
    fRmin(-2.0),
    fRmax(2.0)
{
// TDime constructor:
// Note that there may be only one functional TDime object
// at a time, so it's not use to create more than one instance of it.

// Setup random number generator (if user doesn't set it, doesn't work without)
AliDimeRndm::SetDimeRandom(new TRandom3(123456789)); // default seed 123456789
}

//______________________________________________________________________________
TDime::~TDime()
{
// Destructor
}

void TDime::Initialize() {

    VARS.rts = fEfrm;
    CUTS.rmax = fRmax;
    CUTS.rmin = fRmin;
    CUTS.ecut = fEcut;
    Int_t len = 10; // char array length in DIME

    // Resize, in a case of user not padding the string with spaces
    fProcess.Resize(len);
    fFormf.Resize(len);
    fFsi.Resize(len);

    if (  fProcess.CompareTo("pipm      ") != 0 && 
          fProcess.CompareTo("pi0       ") != 0 &&
          fProcess.CompareTo("kpkm      ") != 0 &&
          fProcess.CompareTo("ks        ") != 0 &&
          fProcess.CompareTo("rho       ") != 0) {

      printf("TDime::Initialize() : Unknown process flag: %s \n", fProcess.Data());
      return;
    }
    if (  fFormf.CompareTo("orexp     ") != 0 && 
          fFormf.CompareTo("exp       ") != 0 &&
          fFormf.CompareTo("power     ") != 0) {

      printf("TDime::Initialize() : Unknown form factor: %s\n", fFormf.Data());
      return;
    }
    if (  fFsi.CompareTo("true      ") != 0 && 
          fFsi.CompareTo("false     ") != 0) {

      printf("TDime::Initialize() : Fsi flag invalid: %s (\"true\" or \"false\") \n", fFsi.Data());
      return;
    }

    // Now copy those to FORTRAN
    strncpy(FLAGS.pflag,     fProcess.Data(), len);
    strncpy(FLAGS.fsi,       fFsi.Data(),     len); // #len chars
    strncpy(FLAGS.ppbar,     "false     ",    len); // #len chars
    strncpy(FLAGS.cuts,      "true      ",    len); // #len chars
    strncpy(FLAGS.unw,       "true      ",    len); // <- Always true!
    strncpy(FF.formf,        fFormf.Data(),   len);
    VARS.iin = 1;
    dimeinit();

    return;
}

void TDime::GenerateEvent()
{
  Int_t ok = 0;
  while(!ok)
    dimegenerate(ok);
   //for (Int_t i = 0; i < HEPEUP.NUP; i++) {
   //printf("%5d %5d %5d %5d %5d %13.3f %13.3f\n", i,
   //	 HEPEUP.IDUP[i], HEPEUP.ISTUP[i], HEPEUP.MOTHUP[i][0],
   //	 HEPEUP.ICOLUP[i][0], HEPEUP.PUP[i][2], HEPEUP.VTIMUP[i]);
   //}

}


TObjArray* TDime::ImportParticles(Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
    fParticles->Clear();
    Int_t nump = 0;

    const Int_t F2C = -1; // Fortran array to C offset (indexing)

    Int_t numpart = HEPEUP.NUP;
    printf("\n TDime:: Full DIME stack contains %d particles \n", numpart);

    if (!strcmp(option,"") || !strcmp(option,"Final")) {
	for (Int_t i = 0; i < numpart; i++) {

	    if (HEPEUP.ISTUP[i] == 1) { // is final
//
//  Use the common block values for the TParticle constructor
//
		TParticle* p = new TParticle(
					     HEPEUP.IDUP[i], HEPEUP.ISTUP[i], HEPEUP.MOTHUP[i][0] + F2C, HEPEUP.MOTHUP[i][1] + F2C,
		    -1, -1,
		    HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
		    0.0, 0.0, 0.0, 0.0);
		fParticles->Add(p);
		++nump;
	    }
	}
    }
    else if (!strcmp(option,"All")) {
	nump = numpart;
	for (Int_t i = 0; i < numpart; i++) {
           Int_t iParent = HEPEUP.MOTHUP[i][0]-1;
	    if (iParent >= 0) {
                TParticle *mother = (TParticle*) (fParticles->UncheckedAt(iParent));
		mother->SetLastDaughter(i);
		if (mother->GetFirstDaughter()==-1)
		    mother->SetFirstDaughter(i);
	    }

	    TParticle* p = new TParticle(
					 HEPEUP.IDUP[i], HEPEUP.ISTUP[i], HEPEUP.MOTHUP[i][0] + F2C, HEPEUP.MOTHUP[i][1] + F2C,
		    -1, -1,
		    HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
		    0.0, 0.0, 0.0, 0.0);
	    fParticles->Add(p);
	}
    }
    return fParticles;
}

Int_t TDime::ImportParticles(TClonesArray *particles, Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  if (particles == 0) return 0;
  TClonesArray &particlesR = *particles;
  particlesR.Clear();
  Int_t nump = 0;

  const Int_t F2C = -1; // Fortran array to C offset (indexing)

  Int_t numpart = HEPEUP.NUP;
  printf("\n TDime:: Full DIME stack contains %d particles\n", numpart);

  if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (Int_t i = 0; i < numpart; i++) {
	if (HEPEUP.ISTUP[i] == 1) { // is final
//
//  Use the common block values for the TParticle constructor
//
	    new(particlesR[i]) // [i], This means there are empty/null elements in TClonesArray for non-final states
	      TParticle(
			HEPEUP.IDUP[i], HEPEUP.ISTUP[i], HEPEUP.MOTHUP[i][0] + F2C, HEPEUP.MOTHUP[i][1] + F2C,
			-1, -1,
			HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
			0.0, 0.0, 0.0, 0.0);
    	    ++nump;
	  }
      }
  }
  else if (!strcmp(option,"All")) {
      nump = numpart;
      for (Int_t i = 0; i < numpart; i++) {

	Int_t iParent = HEPEUP.MOTHUP[i][0]-1;

	if (iParent >= 0) {
	  TParticle *mother = (TParticle*) (particlesR.UncheckedAt(iParent));
	  mother->SetLastDaughter(i);
	  if (mother->GetFirstDaughter()==-1)
	    mother->SetFirstDaughter(i);
	}

	  new(particlesR[i]) TParticle(
				       HEPEUP.IDUP[i], HEPEUP.ISTUP[i], HEPEUP.MOTHUP[i][0] + F2C, HEPEUP.MOTHUP[i][1] + F2C,
				       -1, -1,
				       HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
				       0.0, 0.0, 0.0, 0.0);
      }
  }
  return nump;
}
