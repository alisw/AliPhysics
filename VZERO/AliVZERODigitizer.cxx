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
 
//_________________________________________________________________________
//
// This  constructs Digits out of Hits
//

// --- Standard library ---
#include <Riostream.h>
#include <stdlib.h>

// --- ROOT system ---
#include <TFile.h>
#include <TFolder.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualMC.h>

// --- AliRoot header files ---
#include "AliVZEROConst.h"
#include "AliRun.h"
#include "AliVZERO.h"
#include "AliVZEROhit.h"
#include "AliHit.h"
#include "AliDetector.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliConfig.h"
#include "AliRunDigitizer.h"
#include "AliDigitizer.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliVZERODigitizer.h"
#include "AliVZEROdigit.h"
#include "AliMC.h"

ClassImp(AliVZERODigitizer)

 AliVZERODigitizer::AliVZERODigitizer()
{
   if (!fDigits) fDigits = new TClonesArray("AliVZEROdigit", 1000);
   
   fNevents = 0 ;     
   fDigits = 0 ;
   fNdigits = 0;
   fHits = 0 ;
   fRunLoader = 0;
  
   fPhotoCathodeEfficiency =   0.18;
   fPMVoltage              =  768.0;
   fPMGain = TMath::Power((fPMVoltage / 112.5) ,7.04277);     
}

//____________________________________________________________________________ 
  AliVZERODigitizer::AliVZERODigitizer(AliRunDigitizer* manager)
                    :AliDigitizer(manager)
                    
{
  // constructor
  
  fNevents = 0;     
  fDigits = 0;
  fNdigits = 0;
  fHits = 0;
  fRunLoader = 0;
  
  fPhotoCathodeEfficiency =   0.18;
  fPMVoltage              =  768.0;
  fPMGain = TMath::Power( (fPMVoltage / 112.5) ,7.04277 );
}
           
//____________________________________________________________________________ 
  AliVZERODigitizer::~AliVZERODigitizer()
{
  // destructor
  
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits=0; }
}

//____________________________________________________________________________
void AliVZERODigitizer::OpengAliceFile(const char *file)
{
  // Loads galice.root file and corresponding header, kinematics
  // hits and  digits 
  
  fRunLoader = AliRunLoader::Open(file,AliConfig::GetDefaultEventFolderName(),
                                  "UPDATE");
  
  if (!fRunLoader)
   {
     Error("Open","Can not open session for file %s.",file);
   }
  
  fRunLoader->LoadgAlice();
  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();

  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice)
    { printf("<AliVZEROdigitizer::Open> ");
      printf("AliRun object found on file.\n");}
  else
    { printf("<AliVZEROdigitizer::Open> ");
      printf("Could not find AliRun object.\n");}
    
  // Initialise Hit and Digit arrays
  fHits   = new TClonesArray ("AliVZEROhit", 1000);
  fDigits = new TClonesArray ("AliVZEROdigit", 1000);

  fVZERO  = (AliVZERO*) gAlice->GetDetector("VZERO");
  fVZEROLoader = fRunLoader->GetLoader("VZEROLoader");
  
  if (fVZEROLoader == 0x0){
      cerr<<"Hits2Digits : Can not find VZERO or VZEROLoader\n";}
    
  Int_t retval = fVZEROLoader->LoadHits("read");
  if (retval){
     Error("Open","Error occured while loading hits... Exiting.");
     return;}
      
  fVZEROLoader->LoadDigits("recreate");  
}

//____________________________________________________________________________
void AliVZERODigitizer::Exec() 
 { 

  Int_t nbytes;
  fNdigits = 0;
  Int_t N; 
  Int_t map[96];
  Int_t cell = 0;
  Float_t cPM = fPhotoCathodeEfficiency * fPMGain;
	     
  for(Int_t i=0; i<96; i++) map[i] = 0; 
          
  fNevents = (Int_t) fRunLoader->TreeE()->GetEntries ();  
  printf(" Number of events in file =  %d \n", fNevents);
  
  for (Int_t ievent = 0; ievent < fNevents; ievent++){
      
      for(Int_t i=0; i<96; i++) map[i] = 0;     
      
      fRunLoader->GetEvent(ievent);
      
      fTreeH = fVZEROLoader->TreeH();
      if (fTreeH == 0x0)
       { Error("Exec","Cannot get TreeH");
         return; }
             
      fTreeD = fVZEROLoader->TreeD();
      if (fTreeD == 0x0)
      { fVZEROLoader->MakeTree("D");
        fVZEROLoader->MakeDigitsContainer();
        fTreeD = fVZEROLoader->TreeD(); }
	
      Int_t bufsize = 16000;
      fTreeD->Branch("VZERODigit", &fDigits, bufsize); 
      
//    Now make Digits from hits

      if (fVZERO)
       {
         fHits = fVZERO->Hits();
         
         Int_t ntracks = (Int_t) fTreeH->GetEntries ();
//	 printf(" Number of Tracks in the TreeH = %d \n", ntracks);
         for (Int_t track = 0; track < ntracks; track++)
           {
             gAlice->ResetHits ();
             nbytes += fTreeH->GetEvent(track);
             fParticle = fRunLoader->Stack()->Particle(track);
             Int_t nhits = fHits->GetEntriesFast();
             for (Int_t hit = 0; hit < nhits; hit++)
                 {
                   fVZEROHit = (AliVZEROhit *)fHits->UncheckedAt(hit);
		   N    = fVZEROHit->Nphot();
		   cell = fVZEROHit->Cell();			 		
		   map[cell] = map[cell] + N;
                 }           // hit   loop
           }                 // track loop
            
           Int_t icount = 0; 
	   
           for(Int_t i=0; i<96; i++) {
	      Float_t q1 = Float_t ( map[i] )* cPM * kQe;
	      Float_t noise = gRandom->Gaus(10.5,3.22);
              Float_t PMresponse  =  q1/kC*TMath::Power(ktheta/kthau,1/(1-ktheta/kthau)) 
	                          + noise*1e-3;
	      map[i] = Int_t( PMresponse * 200.0);
	      if(map[i] > 3) {
	         icount++;
//	         printf(" Event, cell, adc = %d %d %d\n", ievent, i, map[i]);
	         AddDigit(ievent, i, map[i]);} 
            }

           fTreeD->Reset();
           fTreeD->Fill();
	   ResetDigit();
           fVZEROLoader->WriteDigits("OVERWRITE");  
	   fVZEROLoader->UnloadDigits();     
      }                     // VZERO loop
    }  //event loop
}

//____________________________________________________________________________
void AliVZERODigitizer::AddDigit(Int_t eventnumber, Int_t cellnumber, Int_t adc) 
 { 
 
// Adds Digit 
 
  TClonesArray &ldigits = *fDigits;  
  new(ldigits[fNdigits++]) AliVZEROdigit(eventnumber,cellnumber,adc);
}
//____________________________________________________________________________
void AliVZERODigitizer::ResetDigit()
{
// Clears Digits
  fNdigits = 0;
  if (fDigits) fDigits->Clear();
}
