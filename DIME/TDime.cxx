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

/* $Id:*/

#include <TROOT.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include "DCommon.h"
#include "TDime.h"

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
  fEfrm(5500.),
  fProcess("rho       "),
  fEcut(0),
  fRmin(-1.8),
  fRmax(1.8)
{
// Default constructor 
}

//______________________________________________________________________________
TDime::TDime(Double_t efrm) :
    TGenerator("Dime","Dime"),
    fEfrm(efrm),
    fProcess("rho      "),
    fEcut(0.),
    fRmin(-1.8),
    fRmax(1.8)
{
// TDime constructor: 
// Note that there may be only one functional TDime object
// at a time, so it's not use to create more than one instance of it.
}

//______________________________________________________________________________
TDime::~TDime()
{
// Destructor
}

  void  TDime::Initialize()
{
    VARS.rts = fEfrm;  
    CUTS.rmax = fRmax; 
    CUTS.rmin = fRmin;
    CUTS.ecut = fEcut;
    Int_t len = 10;
    //    strncpy(FLAGS.pflag,     "rho       ", len);
    strncpy(FLAGS.pflag,     fProcess.Data(), len);
    strncpy(FLAGS.fsi,       "true      ", len);
    strncpy(FLAGS.ppbar,     "false     ", len);
    strncpy(FLAGS.cuts,      "true      ", len);
    strncpy(FLAGS.unw,       "true      ", len);
    strncpy(FF.formf,        "orexp     ", len);
    VARS.iin = 1;
    dimeinit();
}

void  TDime::GenerateEvent()
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

    Int_t numpart = HEPEUP.NUP;
    printf("\n TDime: DIME stack contains %d particles.", numpart);

    if (!strcmp(option,"") || !strcmp(option,"Final")) {
	for (Int_t i = 0; i < numpart; i++) {
	  
	    if (HEPEUP.ISTUP[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
		nump++;
		TParticle* p = new TParticle(
					     HEPEUP.IDUP[i], HEPEUP.MOTHUP[i][0], HEPEUP.MOTHUP[i][1] ,
		    -1, -1, -1,
		    HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
		    0., 0., 0., 0.0);
		fParticles->Add(p);
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
					 HEPEUP.IDUP[i], HEPEUP.MOTHUP[i][0]-1, HEPEUP.MOTHUP[i][1]-1 ,
		    -1, -1, -1,
		    HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
		    0., 0., 0., 0.);
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

  Int_t numpart = HEPEUP.NUP;
  printf("\n TDime: DIME stack contains %d particles.", numpart);
 
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (Int_t i = 0; i < numpart; i++) {
	if (HEPEUP.ISTUP[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
	    nump++;
	    new(particlesR[i]) 
	      TParticle(
			HEPEUP.IDUP[i], HEPEUP.MOTHUP[i][0], HEPEUP.MOTHUP[i][1] ,
			-1, -1, -1,
			HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
			0., 0., 0., 0.0);
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
				       HEPEUP.IDUP[i], HEPEUP.MOTHUP[i][0], HEPEUP.MOTHUP[i][1] ,
				       -1, -1, -1,
				       HEPEUP.PUP[i][0], HEPEUP.PUP[i][1], HEPEUP.PUP[i][2], HEPEUP.PUP[i][3] ,
				       0., 0., 0., 0.0
				       );
      }
  }
  return nump;
}

