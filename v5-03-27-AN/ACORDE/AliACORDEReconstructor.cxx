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

/* $Id: AliACORDEReconstructor.cxx 20956 2007-09-26 14:22:18Z mrodrigu $ */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Class for ACORDE reconstruction                                         //
//////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"

#include "AliACORDEReconstructor.h"
#include "AliACORDERawStream.h"
#include "AliESDEvent.h"
#include "AliACORDEdigit.h"
#include "AliACORDERecoParam.h"

ClassImp(AliACORDEReconstructor)

AliACORDEReconstructor:: AliACORDEReconstructor():
  AliReconstructor(),
  fESDACORDE(0x0),
  fAcordeRecoParam(0x0),
  fCalibData(0x0),
  fDigitsArray(0)
{
  // Default constructor  
  // Get calibration data

  fCalibData = GetCalibData();
  fAcordeRecoParam = GetRecoParam();
}

//_______________________________________________________________________
AliACORDECalibData *AliACORDEReconstructor::GetCalibData() const
{
  return 0x0;
}
//____________________________________________________________________________
AliACORDERecoParam *AliACORDEReconstructor::GetRecoParam() const
{
  return 0x0;
}
//_____________________________________________________________________________
AliACORDEReconstructor& AliACORDEReconstructor::operator = 
  (const AliACORDEReconstructor& /*reconstructor*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliACORDEReconstructor::~AliACORDEReconstructor()
{
// destructor
  delete fESDACORDE;
  delete fDigitsArray;
}

//_____________________________________________________________________________
void AliACORDEReconstructor::Init()
{
// initializer
    fESDACORDE  = new AliESDACORDE;
}

void AliACORDEReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  if (!fDigitsArray)
    fDigitsArray = new TClonesArray("AliACORDEdigit", 60);

  digitsTree->Branch("ACORDEdigit", &fDigitsArray);

  rawReader->Reset();
  AliACORDERawStream rawStream(rawReader);
  if (rawStream.Next()) {
    for(Int_t iChannel = 0; iChannel < 60; iChannel++) {
      Int_t  index = iChannel / 30;
      Int_t  bit   = iChannel % 30;
      if (rawStream.GetWord(index) & (1 << bit))
        new ((*fDigitsArray)[fDigitsArray->GetEntriesFast()]) AliACORDEdigit(iChannel,0);
    }
  }

  digitsTree->Fill();

  fDigitsArray->Clear();
}

void AliACORDEReconstructor::FillESD(TTree* digitsTree, TTree* /*clustersTree*/,AliESDEvent* esd) const
{

  // fills ESD with ACORDE Digits

  if (!digitsTree)
    {
      AliError("No digits tree!");
      return;
    }

  TBranch* digitBranch = digitsTree->GetBranch("ACORDEdigit");
  if (!digitBranch) {
    AliError("No ACORDE digits branch found!");
    return;
  }
  digitBranch->SetAddress(&fDigitsArray);

  digitsTree->GetEvent(0);

  Bool_t AcoHitSingle[60],AcoHitMulti[60];
  for(Int_t i = 0; i < 60; i++) { AcoHitSingle[i] = AcoHitMulti[i] = kFALSE; }

  Int_t nDigits = fDigitsArray->GetEntriesFast();
    
  for (Int_t d=0; d<nDigits; d++) {    
    AliACORDEdigit* digit = (AliACORDEdigit*) fDigitsArray->At(d);
    Int_t module = digit->GetModule();

    AcoHitSingle[module] = kTRUE;
    AcoHitMulti[module] = kTRUE;
  }  
  if (!esd) {
	AliError("NO ACORDE ESD branch found!");
	return;
}
  TString ActiveTriggerDetector = esd->GetFiredTriggerClasses();
  if (ActiveTriggerDetector.Contains("ASL")) fESDACORDE->SetACORDEBitPattern(AcoHitSingle);
  else if (ActiveTriggerDetector.Contains("AMU")) fESDACORDE->SetACORDEBitPattern(AcoHitMulti);
	else fESDACORDE->SetACORDEBitPattern(AcoHitSingle);

  if (esd)
    {
      AliDebug(1, Form("Writing ACORDE data to ESD Tree"));
      esd->SetACORDEData(fESDACORDE);
    }

  fDigitsArray->Clear();
}


