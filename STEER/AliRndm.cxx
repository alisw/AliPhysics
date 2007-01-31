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
// Random number class for AliRoot
// This class allows to have different 
// random number generator for different
// elements of AliRoot                                                                          //
// It also allows saving and retrieving of the random number seeds
///////////////////////////////////////////////////////////////////////////////

#include <TClass.h>
#include <TFile.h>
#include <TError.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "AliRndm.h"
#include "AliLog.h"

ClassImp(AliRndm)

//_______________________________________________________________________
AliRndm::AliRndm():
  fRandom(gRandom)
{
  // 
  // Default ctor
  //
}

//_______________________________________________________________________
AliRndm::AliRndm(const AliRndm& rn):
  fRandom(gRandom)
{
  //
  // Copy constructor
  //
  rn.Copy(*this);
}

//_______________________________________________________________________
void AliRndm::Copy(AliRndm&) const
{
  AliFatalClass("Not implemented");
}


//_____________________________________________________________________________
void AliRndm::Rndm(Float_t* array, Int_t size) const
{
  //
  // Return an array of n random numbers uniformly distributed 
  // between 0 and 1 not included
  //
  for(Int_t i=0; i<size; i++) 
#ifdef CKNONE
    array[i]=fRandom->Rndm();
#else
    do array[i]=fRandom->Rndm(); while(0>=array[i] || array[i]>=1);
#endif
}

//_____________________________________________________________________________
void AliRndm::ReadRandom(const char *filename)
{
  //
  // Reads saved random generator status from filename
  //
  char *fntmp = gSystem->ExpandPathName(filename);
  TFile *file = new TFile(fntmp,"r");
  delete [] fntmp;
  if(!file) {
    AliErrorClass(Form("Could not open file %s",filename));
  } else {
    if(!fRandom) fRandom = new TRandom();
    fRandom->Read("Random");
    file->Close();
    delete file;
  }
}

//_____________________________________________________________________________
void AliRndm::WriteRandom(const char *filename) const
{
  //
  // Writes random generator status to filename
  //
  char *fntmp = gSystem->ExpandPathName(filename);
  TFile *file = new TFile(fntmp,"new");
  delete [] fntmp;
  if(!file) {
    AliErrorClass(Form("Could not open file %s",filename));
  } else {
    fRandom->Write();
    file->Close();
    delete file;
  }
}
