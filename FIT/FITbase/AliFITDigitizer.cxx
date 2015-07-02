/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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


/******************************************************************
 *    Produde digits from hits
 *       digits is TObject and includes
 *	We are writing array if C & A  TDC
 *	C & A  ADC (will need for slow simulation)
 *	TOF first particle C & A
 *	mean time and time difference (vertex position)
 *
 *      Alla.Maevskaya@cern.ch 
 ****************************************************************/


#include <TArrayI.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h> 

#include "AliLog.h"
#include "AliFITDigitizer.h"
#include "AliFIT.h"
#include "AliFITHits.h"
#include "AliFITDigit.h"
#include "AliDigitizationInput.h"
#include "AliRun.h"
#include <AliLoader.h>
#include <AliRunLoader.h>
#include <stdlib.h>

ClassImp(AliFITDigitizer)

//___________________________________________
  AliFITDigitizer::AliFITDigitizer()  :AliDigitizer(),
				     fFIT(0),
				     fHits(0),
				     fDigits(0),
                                     fNdigits(0)
{
// Default ctor - don't use it

}

//___________________________________________
AliFITDigitizer::AliFITDigitizer(AliDigitizationInput* digInput) 
  :AliDigitizer(digInput),
   fFIT(0),
   fHits(0),
   fDigits(0),
   fNdigits(0)
{
// ctor which should be used
 
}


//------------------------------------------------------------------------
AliFITDigitizer::~AliFITDigitizer()
{
// Destructor

  AliDebug(1,"FIT");

 }

//------------------------------------------------------------------------
Bool_t AliFITDigitizer::Init()
{
// Initialization
  AliDebug(1," Init");
 return kTRUE;

}
 
//---------------------------------------------------------------------
void AliFITDigitizer::Digitize(Option_t* /*option*/)
{

  /*
    Produde digits from hits
    digits is TObject and includes
    We are writing array if C & A for each channel CFD, LED, QT0 and QT1
    C & A  ADC (will need for slow simulation)
  */
  
  
  
  //output loader 
  AliDebug(1,"start...");
  //input loader
  //
  // From hits to digits
  //
  
  AliRunLoader *outRL = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());
  AliLoader * outFitLoader = outRL->GetLoader("FITLoader");
  
  fFIT  = static_cast<AliFIT*>(gAlice->GetDetector("FIT"));
  if (!fFIT) {
    AliError("Can not get FIT from gAlice");
    return;
  }  
  fFIT->ResetDigits();
  
  DigitizeHits();
  
  //load digits    
  outFitLoader->LoadDigits("UPDATE");
  TTree *treeD  = outFitLoader->TreeD();
  if (treeD == 0x0) {
    outFitLoader->MakeTree("D");
    treeD = outFitLoader->TreeD();
  }
  treeD->Reset();
  fFIT  = (AliFIT*)outRL ->GetAliRun()->GetDetector("FIT");
  // Make a branch in the tree 
  fFIT->MakeBranch("D");
  treeD->Fill();
  
  outFitLoader->WriteDigits("OVERWRITE");
  fFIT->ResetDigits();
  outFitLoader->UnloadDigits();
}
//____________________________________________________________________________
void AliFITDigitizer::DigitizeHits()
{

  Int_t hit, nhits;
  Float_t time[240], besttime[240];
  Int_t countE[240];
  Int_t timeCFD, timeLED, timeQT1, timeQT0;

  Int_t threshold = 0; //photoelectrons
  Float_t channelWidth = 24.4 ; 
  Int_t pmt, mcp, volume, qt; 
  //eqailized distance from IP
  Float_t zdetC = 84;
  Float_t zdetA = 335.;
  Float_t c = 0.0299792458; // cm/ps
  Float_t eqdistance = (zdetA - zdetC) /c;
  Float_t ph2Mip = 318; 
  
  AliFITHits  *startHit;
  TBranch *brHits=0;
 
  Int_t nFiles=fDigInput->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {
    if (inputFile < nFiles-1) {
      AliWarning(Form("ignoring input stream %d", inputFile));
      continue;
     }

  for (Int_t i0=0; i0<240; i0++)
      {
	time[i0]=besttime[i0]=999999; countE[i0]=0;
      }
  AliRunLoader * inRL = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inputFile));
  AliLoader * fitLoader = inRL->GetLoader("FITLoader");
  if (!inRL->GetAliRun()) inRL->LoadgAlice();
  Int_t numpmt;
  //read Hits 
  fitLoader->LoadHits("READ");//probably it is necessary to load them before
  fHits = fFIT->Hits ();
  TTree *th = fitLoader->TreeH();
  brHits = th->GetBranch("FIT");
  if (brHits) {
    fFIT->SetHitsAddressBranch(brHits);
  }else{
    AliWarning("Branch FIT hit not found for this event");
    continue;      
  } 
  Int_t ntracks    = (Int_t) th->GetEntries();
  if (ntracks<=0) return;
  // Start loop on tracks in the hits containers
  for (Int_t track=0; track<ntracks;track++) {
    brHits->GetEntry(track);
    nhits = fHits->GetEntriesFast();
    for (hit=0;hit<nhits;hit++) 
      {
	startHit   = (AliFITHits*) fHits->UncheckedAt(hit);
	if (!startHit) {
	  AliError("The unchecked hit doesn't exist");
	  break;
	}
	Float_t ipart = startHit->Particle();
	if (ipart<49) continue;
	pmt = startHit->Pmt();
	mcp = startHit->MCP();
	volume = startHit->Volume();
	if(volume==2) continue;
	Float_t z = startHit->Z();
	numpmt= 4*mcp + pmt;
	besttime[numpmt] = startHit->Time();
	if(besttime[numpmt]<time[numpmt]) time[numpmt]=besttime[numpmt];
	countE[numpmt]++;
      } //hits loop
  } //track loop
  
  for (Int_t ipmt=0; ipmt<240; ipmt++)
    {
      if (countE[ipmt]>threshold && time[ipmt]<50000 && time[ipmt]>0 ) {
	//fill ADC
	// QTC procedure:
	// 1MIP ->318phe  ;
	qt= 1000* countE[ipmt] /ph2Mip;  // 318 ph/Mip 
	//  fill TDC
	if (ipmt>95) time[ipmt] = time[ipmt] + eqdistance;
	timeCFD = Int_t (gRandom->Gaus(time[ipmt], 50)/channelWidth ); 
	timeLED =  Int_t (time[ipmt]/channelWidth );
	timeQT0 = 0;
	timeQT1 = qt ;
	AliDebug(1,Form("Digits:::::  numpmt %i  time CFD %i QTC  %i :: counts %i \n ",  ipmt, timeCFD, timeQT1,  countE[ipmt]) );
	fFIT-> AddDigit(ipmt,   timeCFD, timeLED, timeQT0,  timeQT1, 0);
      } //hitted PMTs
    } //pmt loop
  fitLoader->UnloadHits();
  
  }
}

//____________________________________________________________________________
void AliFITDigitizer::AddDigit(Int_t npmt,  
				Int_t timeCFD, Int_t timeLED, Int_t timeQT0, 
				Int_t timeQT1) 
 { 
 
// Adds Digit 
 
  TClonesArray &ldigits = *fDigits;  
	 
  new(ldigits[fNdigits++]) AliFITDigit(npmt,  
				       timeCFD, timeLED, timeQT0, timeQT1);
	 
}
//____________________________________________________________________________
void AliFITDigitizer::ResetDigits()
{

// Clears Digits

  fNdigits = 0;
  if (fDigits) fDigits->Clear();
}
