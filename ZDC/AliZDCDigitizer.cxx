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
//  			ZDC digitizer class                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

// --- ROOT system
#include <TTree.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TRandom.h>

// --- AliRoot header files
#include "AliLog.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliZDCSDigit.h"
#include "AliZDCDigit.h"
#include "AliZDCFragment.h"
#include "AliZDCDigitizer.h"

ClassImp(AliZDCDigitizer)


//____________________________________________________________________________
AliZDCDigitizer::AliZDCDigitizer()
{
// Default constructor    

}

//____________________________________________________________________________
AliZDCDigitizer::AliZDCDigitizer(AliRunDigitizer* manager):
  AliDigitizer(manager)
{
// Constructor    

}

//____________________________________________________________________________
AliZDCDigitizer::~AliZDCDigitizer()
{
// Destructor

}


//____________________________________________________________________________
Bool_t AliZDCDigitizer::Init()
{
// Initialize the digitizer

  for (Int_t i = 0; i < 3; i++) {
    for(Int_t j = 0; j < 5; j++) {
      fPMGain[i][j]   = 100000.;
    }
  }

// ADC Caen V965
  fADCRes[0]   = 0.0000008; // ADC Resolution high gain: 200 fC/adcCh
  fADCRes[1]   = 0.0000064; // ADC Resolution low gain:  25  fC/adcCh

  return kTRUE;
}

//____________________________________________________________________________
void AliZDCDigitizer::Exec(Option_t* /*option*/)
{
// Execute digitization

  Float_t pm[3][5];
  for (Int_t iSector1 = 0; iSector1 < 3; iSector1++)
    for (Int_t iSector2 = 0; iSector2 < 5; iSector2++) {
      pm[iSector1][iSector2] = 0;
    }
  // pm[0][...] = light in ZN [C, Q1, Q2, Q3, Q4]
  // pm[1][...] = light in ZP [C, Q1, Q2, Q3, Q4]
  // pm[2][...] = light in ZEM [x, 1, 2, x, x]

  // impact parameter and number of spectators
  Float_t impPar = -1;
  Int_t specN = 0;
  Int_t specP = 0;

  // loop over input streams
  for (Int_t iInput = 0; iInput < fManager->GetNinputs(); iInput++) {

    // get run loader and ZDC loader
    AliRunLoader* runLoader = 
      AliRunLoader::GetRunLoader(fManager->GetInputFolderName(iInput));
    AliLoader* loader = runLoader->GetLoader("ZDCLoader");
    if (!loader) continue;

    // load sdigits
    loader->LoadSDigits();
    TTree* treeS = loader->TreeS();
    if (!treeS) continue;
    AliZDCSDigit sdigit;
    AliZDCSDigit* psdigit = &sdigit;
    treeS->SetBranchAddress("ZDC", &psdigit);

    // loop over sdigits
    for (Int_t iSDigit = 0; iSDigit < treeS->GetEntries(); iSDigit++) {
      treeS->GetEntry(iSDigit);
      if (!psdigit) continue;

      if ((sdigit.GetSector(1) < 0) || (sdigit.GetSector(1) > 4)) {
	AliError(Form("\nsector[0] = %d, sector[1] = %d\n", 
                      sdigit.GetSector(0), sdigit.GetSector(1)));
	continue;
      }
      pm[sdigit.GetSector(0)-1][sdigit.GetSector(1)] += sdigit.GetLightPM();
      //printf("\n\t sector[0] = %d, sector[1] = %d, pm[%d][%d] = %f \n",
      //	  sdigit.GetSector(0), sdigit.GetSector(1),sdigit.GetSector(0)-1,
      //	  sdigit.GetSector(1), pm[sdigit.GetSector(0)-1][sdigit.GetSector(1)]); // Chiara debugging!
    }

    // unload sdigits
    loader->UnloadSDigits();

    // get the impact parameter and the number of spectators in case of hijing
    if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
    AliHeader* header = runLoader->GetAliRun()->GetHeader();
    if (!header) continue;
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if (!genHeader) continue;
    if (!genHeader->InheritsFrom(AliGenHijingEventHeader::Class())) continue;
    impPar = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
    // Until there is only 1 ZDC set the # of spectators must be divided by 2!!!
    specN = ((AliGenHijingEventHeader*) genHeader)->Spectatorsn() / 2;
    specP = ((AliGenHijingEventHeader*) genHeader)->Spectatorsp() / 2;
    AliDebug(2, Form("\n b = %f fm, Nspecn = %d, Nspecp = %d\n",
                     impPar, specN, specP));
  }

  // add spectators
  if (impPar >= 0) {
    Int_t freeSpecN, freeSpecP;
    Fragmentation(impPar, specN, specP, freeSpecN, freeSpecP);
    SpectatorSignal(1, freeSpecN, pm);
    SpectatorSignal(2, freeSpecP, pm);
  }

  // get the output run loader and loader
  AliRunLoader* runLoader = 
    AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  AliLoader* loader = runLoader->GetLoader("ZDCLoader");
  if (!loader) {
    AliError("no ZDC loader found");
    return;
  }

  // create the output tree
  const char* mode = "update";
  if (runLoader->GetEventNumber() == 0) mode = "recreate";
  loader->LoadDigits(mode);
  loader->MakeTree("D");
  TTree* treeD = loader->TreeD();
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;
  const Int_t kBufferSize = 4000;
  treeD->Branch("ZDC", "AliZDCDigit", &pdigit, kBufferSize);

  // Create digits
  Int_t sector[2];
  Int_t digi[2];
  for (sector[0] = 1; sector[0] <= 3; sector[0]++)
    for (sector[1] = 0; sector[1] < 5; sector[1]++)  {
        if ((sector[0] == 3) && ((sector[1] < 1) || (sector[1] > 2))) continue;
        for (Int_t res = 0; res < 2; res++){
          digi[res] = Phe2ADCch(1, 0, pm[sector[0]-1][sector[1]], res) 
	            + Pedestal();
	  //printf("\n\t sector[0] = %d, sector[1] = %d, digi[%d] = %d\n",
	  //	sector[0], sector[1], res, digi[res]); // Chiara debugging!
	}
        new(pdigit) AliZDCDigit(sector, digi);
        treeD->Fill();
    }

  // write the output tree
  loader->WriteDigits("OVERWRITE");
  loader->UnloadDigits();
}


//_____________________________________________________________________________
void AliZDCDigitizer::Fragmentation(Float_t impPar, Int_t specN, Int_t specP,
                                    Int_t &freeSpecN, Int_t &freeSpecP) const
{
// simulate fragmentation of spectators

  Int_t zz[100], nn[100];
  AliZDCFragment frag(impPar);
  for (Int_t j=0; j<=99; j++){
     zz[j] =0;
     nn[j] =0;
  }

  // Fragments generation
  Int_t nAlpha;
  frag.GenerateIMF(zz, nAlpha);

  // Attach neutrons
  Int_t ztot=0;
  Int_t ntot=0;
  frag.AttachNeutrons(zz, nn, ztot, ntot);
  freeSpecN = specN-ntot-2*nAlpha;
  freeSpecP = specP-ztot-2*nAlpha;
  if(freeSpecN<0) freeSpecN=0;
  if(freeSpecP<0) freeSpecP=0;
  AliDebug(2, Form("FreeSpn = %d, FreeSpp = %d", freeSpecN, freeSpecP));
}

//_____________________________________________________________________________
void AliZDCDigitizer::SpectatorSignal(Int_t SpecType, Int_t numEvents, 
                                      Float_t pm[3][5]) const
{
// add signal of the spectators
 
  TFile* file = NULL;
  if (SpecType == 1) {		// --- Signal for spectator neutrons
    file = TFile::Open("$ALICE/$ALICE_LEVEL/ZDC/ZNsignalntu.root");
  } else if (SpecType == 2) {	// --- Signal for spectator protons
    file = TFile::Open("$ALICE/$ALICE_LEVEL/ZDC/ZPsignalntu.root");
  }
  if (!file || !file->IsOpen()) {
    AliError("Opening of file failed");
    return;
  }

  TNtuple* zdcSignal = (TNtuple*) file->Get("ZDCSignal");
  Int_t nentries = (Int_t) zdcSignal->GetEntries();
  
  Float_t *entry, hitsSpec[7];
  Int_t pl, i, j, k, iev=0, rnd[125], volume[2];
  for(pl=0;pl<125;pl++){
     rnd[pl] = 0;
  }
  if (numEvents > 125) {
    AliWarning(Form("numEvents (%d) is larger than 125", numEvents));
    numEvents = 125;
  }
  for(pl=0;pl<numEvents;pl++){
     rnd[pl] = (Int_t) (9999*gRandom->Rndm());
     if(rnd[pl] >= 9998) rnd[pl] = 9997;
     //printf("	rnd[%d] = %d\n",pl,rnd[pl]);     
  }
  // Sorting vector in ascending order with C function QSORT 
  qsort((void*)rnd,numEvents,sizeof(Int_t),comp);
  do{
     for(i=0; i<nentries; i++){  
  	zdcSignal->GetEvent(i);
  	entry = zdcSignal->GetArgs();
  	if(entry[0] == rnd[iev]){
          for(k=0; k<2; k++) volume[k] = (Int_t) entry[k+1];
          for(j=0; j<7; j++) hitsSpec[j] = entry[j+3];

	  Float_t lightQ = hitsSpec[4];
	  Float_t lightC = hitsSpec[5];
	  AliDebug(3, Form("volume = (%d, %d), lightQ = %.3f, lightC = %.3f",
                           volume[0], volume[1], lightQ, lightC));
	  if (volume[0] < 3) {  // ZN or ZP
            pm[volume[0]-1][0] += lightC;
            pm[volume[0]-1][volume[1]] += lightQ;
	  } else {
            if (volume[1] == 1) pm[2][1] += lightC;
            else                pm[2][2] += lightQ;
	}


  	}
  	else if(entry[0] > rnd[iev]){
	  iev++;
	  continue;
	}
     }
  }while(iev<numEvents);
  
  file->Close();
  delete file;
}


//_____________________________________________________________________________
Int_t AliZDCDigitizer::Phe2ADCch(Int_t Det, Int_t Quad, Float_t Light, 
                                 Int_t Res) const
{
// Evaluation of the ADC channel corresponding to the light yield Light

  return (Int_t) (Light * fPMGain[Det-1][Quad] * fADCRes[Res]);
}

//_____________________________________________________________________________
Int_t AliZDCDigitizer::Pedestal() const
{
  // --- Pedestal value -> extracted from a gaussian distribution
  // obtained from the beam test on the ZEM prototype (Aug. 2000)

  Float_t pedMean  = 50.;
  Float_t pedWidth = 5.;

  return (Int_t) gRandom->Gaus(pedMean, pedWidth);
}
