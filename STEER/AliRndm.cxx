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

/*
$Log$
Revision 1.1  2000/11/30 07:12:48  alibrary
Introducing new Rndm and QA classes

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include "TFile.h"

#include "AliRndm.h"
#include "TRandom3.h"

ClassImp(AliRndm)


//_____________________________________________________________________________
void AliRndm::Rndm(Float_t* array, const Int_t size) const
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
  char *fntmp = gSystem->ExpandPathName(filename);
  TFile *file = new TFile(fntmp,"r");
  delete [] fntmp;
  if(!file) {
    printf("AliRndm:: Could not open file %s\n",filename);
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
  char *fntmp = gSystem->ExpandPathName(filename);
  TFile *file = new TFile(fntmp,"new");
  delete [] fntmp;
  if(!file) {
    printf("AliRndm:: Could not open file %s\n",filename);
  } else {
    fRandom->Write();
    file->Close();
    delete file;
  }
}
