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
// This is a TTask that constructs SDigits out of Hits
// A Summable Digits is the sum of all hits in a pad
// 
//
//-- Author: F. Pierella
//////////////////////////////////////////////////////////////////////////////


#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"

#include "AliTOFHitMap.h"
#include "AliTOFSDigit.h"
#include "AliTOFhit.h"
#include "AliTOF.h"
#include "AliTOFv1.h"
#include "AliTOFv2.h"
#include "AliTOFv3.h"
#include "AliTOFv4.h"
#include "AliTOFSDigitizer.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "AliMC.h"

#include "TFile.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

ClassImp(AliTOFSDigitizer)

//____________________________________________________________________________ 
  AliTOFSDigitizer::AliTOFSDigitizer():TTask("AliTOFSDigitizer","") 
{
  // ctor
  fNevents = 0 ;     
//  fSDigits = 0 ;
  fHits = 0 ;

}
           
//____________________________________________________________________________ 
  AliTOFSDigitizer::AliTOFSDigitizer(char* HeaderFile,char *SdigitsFile ):TTask("AliTOFSDigitizer","") 
{
  fNevents = 0 ;     // Number of events to digitize, 0 means all evens in current file
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}

//____________________________________________________________________________ 
  AliTOFSDigitizer::~AliTOFSDigitizer()
{
  // dtor
}

//____________________________________________________________________________
void AliTOFSDigitizer::Exec(Option_t *option) { 


  AliTOF *TOF = (AliTOF *) gAlice->GetDetector ("TOF");

  if (!TOF) {
    Error("AliTOFSDigitizer","TOF not found");
    return;
  }

  if (fNevents == 0)
    fNevents = (Int_t) gAlice->TreeE()->GetEntries();

  for (Int_t ievent = 0; ievent < fNevents; ievent++) {
    gAlice->GetEvent(ievent);
    TTree *TH = gAlice->TreeH ();
    if (!TH)
      return;
    if (gAlice->TreeS () == 0)
      gAlice->MakeTree ("S");

      
    //Make branches
    char branchname[20];
    sprintf (branchname, "%s", TOF->GetName ());
    //Make branch for digits
    TOF->MakeBranch ("S");
    
    //Now made SDigits from hits

    Int_t    vol[5];       // location for a digit
    Float_t  digit[2];     // TOF digit variables
    TParticle *particle;
    AliTOFhit *tofHit;
    TClonesArray *TOFhits = TOF->Hits();

// create hit map
    AliTOFHitMap *hitMap = new AliTOFHitMap(TOF->SDigits());

    Int_t ntracks = static_cast<Int_t>(TH->GetEntries());
    for (Int_t track = 0; track < ntracks; track++)
    {
      gAlice->ResetHits();
      TH->GetEvent(track);
      particle = gAlice->Particle(track);
      Int_t nhits = TOFhits->GetEntriesFast();

      for (Int_t hit = 0; hit < nhits; hit++)
      {
	tofHit = (AliTOFhit *) TOFhits->UncheckedAt(hit);
	vol[0] = tofHit->GetSector();
	vol[1] = tofHit->GetPlate();
	vol[2] = tofHit->GetStrip();
	vol[3] = tofHit->GetPadx();
	vol[4] = tofHit->GetPadz();

	// 95% of efficiency to be inserted here
	// edge effect to be inserted here
	// cross talk  to be inserted here

	Float_t idealtime = tofHit->GetTof(); // unit s
	idealtime *= 1.E+12;  // conversion from s to ps
	// fTimeRes is given usually in ps
	Float_t tdctime   = gRandom->Gaus(idealtime, TOF->GetTimeRes());
	digit[0] = tdctime;

	// typical Landau Distribution to be inserted here
	// instead of Gaussian Distribution
	Float_t idealcharge = tofHit->GetEdep();
	Float_t adccharge = gRandom->Gaus(idealcharge, TOF->GetChrgRes());
	digit[1] = adccharge;
	Int_t tracknum = tofHit->GetTrack();

	// check if two digit are on the same pad; in that case we sum
	// the two or more digits
	if (hitMap->TestHit(vol) != kEmpty) {
	  AliTOFSDigit *sdig = static_cast<AliTOFSDigit*>(hitMap->GetHit(vol));
	  sdig->Update(tdctime,adccharge,tracknum);
	} else {
	  TOF->AddSDigit(tracknum, vol, digit);
	  hitMap->SetHit(vol);
	}
      } // end loop on hits for the current track
    } // end loop on ntracks

    delete hitMap;
      
    gAlice->TreeS()->Reset();
    gAlice->TreeS()->Fill();
    gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
  }				//event loop


}
 
//__________________________________________________________________
void AliTOFSDigitizer::SetSDigitsFile(char * file ){
  if(!fSDigitsFile.IsNull())
    cout << "Changing SDigits file from " <<(char *)fSDigitsFile.Data() << " to " << file << endl ;
  fSDigitsFile=file ;
}
//__________________________________________________________________
void AliTOFSDigitizer::Print(Option_t* option)const
{
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  if(fSDigitsFile.IsNull())
    cout << " Writing SDigitis to file galice.root "<< endl ;
  else
    cout << "    Writing SDigitis to file  " << (char*) fSDigitsFile.Data() << endl ;

}
