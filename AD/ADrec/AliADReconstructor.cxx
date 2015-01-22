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
#include <TParameter.h>

#include "AliRawReader.h"
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"

#include "AliADReconstructor.h"
#include "AliADdigit.h"
#include "AliESDAD.h"
#include "AliADConst.h"
#include "AliADCalibData.h"
#include "AliADRawStream.h"

ClassImp(AliADReconstructor)
//_____________________________________________________________________________
AliADReconstructor:: AliADReconstructor():
  AliReconstructor(),
  fESDAD(0x0),
  fCalibData(NULL),
  fDigitsArray(0)

{
  fCalibData = GetCalibData();
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

//_____________________________________________________________________________
void AliADReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
// converts RAW to digits 

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  if (!fDigitsArray){
    fDigitsArray = new TClonesArray("AliADdigit", 16);
    digitsTree->Branch("ADDigit", &fDigitsArray);
    }

  rawReader->Reset();
  AliADRawStream rawStream(rawReader);
  if (rawStream.Next()) { 

    for(Int_t iChannel=0; iChannel < 16; ++iChannel) {
      Int_t offlineCh = rawStream.GetOfflineChannel(iChannel);
      // ADC charge samples
      Short_t chargeADC[kNClocks];
      for(Int_t iClock=0; iClock < kNClocks; ++iClock) {
	chargeADC[iClock] = rawStream.GetPedestal(iChannel,iClock);
      }
      // Integrator flag
      Bool_t integrator = rawStream.GetIntegratorFlag(iChannel,kNClocks/2);
      Bool_t BBflag = rawStream.GetBBFlag(iChannel,kNClocks/2); 
      Bool_t BGflag = rawStream.GetBGFlag(iChannel,kNClocks/2);
   
      // HPTDC data (leading time and width)
      Int_t board = AliADCalibData::GetBoardNumber(offlineCh);
      Float_t time = rawStream.GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      Float_t width = rawStream.GetWidth(iChannel)*fCalibData->GetWidthResolution(board);
      // Add a digit
      if(!fCalibData->IsChannelDead(iChannel)){
	  new ((*fDigitsArray)[fDigitsArray->GetEntriesFast()]) AliADdigit(offlineCh, time, width,integrator, chargeADC, BBflag, BGflag);
      }
    }

    digitsTree->Fill();
  }

  fDigitsArray->Clear();

}

//_____________________________________________________________________________
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
    Int_t module = digit->PMNumber();
 //   printf("AD Module: %d\n",module);
    ADHits[module] = kTRUE;
  }  
  if (!esd) {
	AliError("NO AD ESD branch found!");
	return;
}
  //fESDAD->SetADBitCell(ADHits);

  if (esd)
    {
      AliDebug(1, Form("Writing AD data to ESD Tree"));
      esd->SetADData(fESDAD);
    }

  fDigitsArray->Clear();
}

//_____________________________________________________________________________
AliADCalibData* AliADReconstructor::GetCalibData() const
{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  //entry = man->Get("AD/Calib/Data");
  //if(!entry){
    //AliWarning("Load of calibration data from default storage failed!");
    //AliWarning("Calibration data will be loaded from local storage ($ALICE_ROOT)");
	
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(1);
    entry = man->Get("AD/Calib/Data");
  //}
  // Retrieval of data in directory AD/Calib/Data:

  AliADCalibData *calibdata = 0;

  if (entry) calibdata = (AliADCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}


