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

/* $Id: AliADReconstructor.cxx 20956 2007-09-26 14:22:18Z mrodrigu $ */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Class for AD reconstruction                                         //
//////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"

#include "AliADReconstructor.h"
#include "AliESDEvent.h"
#include "AliADdigit.h"
#include "AliESDAD.h"

ClassImp(AliADReconstructor)

AliADReconstructor:: AliADReconstructor():
  AliReconstructor(),
  fESDAD(0x0),
  fDigitsArray(0)
{
  // Default constructor  
  // Get calibration data

}

//_____________________________________________________________________________
AliADReconstructor& AliADReconstructor::operator = 
  (const AliADReconstructor& /*reconstructor*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliADReconstructor::~AliADReconstructor()
{
// destructor
  delete fESDAD;
  delete fDigitsArray;
}

//_____________________________________________________________________________
void AliADReconstructor::Init()
{
// initializer
    fESDAD  = new AliESDAD;
}

void AliADReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{

  printf("Converting digits for AD .. not implemented \n");
}

void AliADReconstructor::FillESD(TTree* digitsTree, TTree* /*clustersTree*/,AliESDEvent* esd) const
{

  printf("Running AD Reconstruction \n");

  // fills ESD with AD Digits

  if (!digitsTree)
    {
      AliError("No digits tree!");
      return;
    }

  TBranch* digitBranch = digitsTree->GetBranch("ADdigit");
  if (!digitBranch) {
    AliError("No AD digits branch found!");
    return;
  }
  digitBranch->SetAddress(&fDigitsArray);

  digitsTree->GetEvent(0);

  Bool_t ADHits[16];
  for(Int_t i = 0; i < 16; i++) { ADHits[i] = kFALSE; }

  Int_t nDigits = fDigitsArray->GetEntriesFast();
    
  for (Int_t d=0; d<nDigits; d++) {    
    AliADdigit* digit = (AliADdigit*) fDigitsArray->At(d);
    Int_t module = digit->GetCell();
 //   printf("AD Module: %d\n",module);
    ADHits[module] = kTRUE;
  }  
  if (!esd) {
	AliError("NO AD ESD branch found!");
	return;
}
  fESDAD->SetADBitCell(ADHits);

  if (esd)
    {
      AliDebug(1, Form("Writing AD data to ESD Tree"));
      esd->SetADData(fESDAD);
    }

  fDigitsArray->Clear();
}


