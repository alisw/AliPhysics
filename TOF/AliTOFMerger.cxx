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
//Piotr.Skowronski@cern.ch :
//Corrections applied in order to compile (only) with new I/O and folder structure
//To be implemented correctly by responsible

#include "AliRunLoader.h"
#include "AliLoader.h"

#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>


#include "AliTOFMerger.h"
#include "AliTOF.h"
#include "AliTOFSDigitizer.h"
#include "AliTOFhit.h"
#include "AliTOFdigit.h"

#include "AliRun.h"
#include "AliPDG.h"

#include <stdlib.h>
#include <Riostream.h>
#include <Riostream.h>

ClassImp(AliTOFMerger)

//___________________________________________
  AliTOFMerger::AliTOFMerger() 
{
// Default ctor    
    fNDigits = 0;
    fEvNrSig = 0;
    fEvNrBgr = 0;
    fMerge =kDigitize;
    fDigits = 0;
    fSDigits =0;
    fFnBgr = 0;
    fFnSig = 0;
    fBgrFile = 0;
    fRunLoader     = 0 ;
}

//------------------------------------------------------------------------
AliTOFMerger::~AliTOFMerger()
{
// Dtor
  if(fSDigits)  {
    fSDigits->Delete();
    delete fSDigits ;
    fSDigits = 0;
  }
  delete fBgrFile;
  fBgrFile = 0;
  if(fFnBgr) {
    delete[] fFnBgr;
    fFnBgr = 0;
  }
  if(fFnSig) {
    delete[] fFnSig;
    fFnSig = 0;
  }
  delete fRunLoader;
}


//------------------------------------------------------------------------
void AliTOFMerger::Init()
{
// Initialisation
    if (fMerge) fBgrFile = InitBgr();
    
}



//------------------------------------------------------------------------
TFile* AliTOFMerger::InitBgr()
{
// Initialise background event
   fRunLoader = AliRunLoader::Open(fFnBgr);//open session and mount on default event folder
   
    TFile *file = new TFile(fFnBgr);
// add error checking later
    printf("\n AliTOFMerger has opened %s file with background event \n", fFnBgr);
    return file;
}

//------------------------------------------------------------------------
void AliTOFMerger::Digitise()
{
// as in FMD
// keep galice.root for signal and name differently the file for 
// background when add! otherwise the track info for signal will be lost !

#ifdef DEBUG
  cout<<"ALiTOFMerger::>SDigits2Digits start...\n";
#endif
  Int_t retval;

  if (fRunLoader == 0x0)
   {
     Error("Exec","Event is not loaded. Exiting");
     return;
   }

  if (fRunLoader->GetAliRun() == 0x0) {
    retval = fRunLoader->LoadgAlice();
    if (retval)
      {
	Error("Exec","Error occured while loading gAlice. Exiting");
	return;
      }
  }

  if (fRunLoader->TreeE() == 0x0) {
    retval = fRunLoader->LoadHeader();
    if (retval)
      {
	Error("Exec","Error occured while loading header. Exiting");
	return;
      }
  }

  if (fRunLoader->TreeK() == 0x0) {
    retval = fRunLoader->LoadKinematics("READ");
    if (retval)
      {
	Error("Exec","Error occured while loading kinematics. Exiting");
	return;
      }
  }

  AliLoader* gime = fRunLoader->GetLoader("TOFLoader");
  if (gime == 0x0)
   {
     Error("Exec","Can not find TOF loader in event. Exiting.");
     return;
   }
  gAlice = fRunLoader->GetAliRun();
  
  AliTOF* TOF = (AliTOF *) gAlice->GetDetector("TOF") ;


  TFile *f1 =0;
  TTree *TK = fRunLoader->TreeK();
  if (TK) f1 = TK->GetCurrentFile();

  fRunLoader->GetEvent(fEvNrSig);
  
  if(gime->TreeD() == 0)    	
    gime->MakeTree("D") ;
  
  gime->TreeD()->Reset();

  // read and write digits for signal
   ReadWriteDigit(fEvNrSig);

   if(fMerge)
    { 
      // bgr file
      fBgrFile->cd();
      // gAlice->TreeS()->Reset();
      gAlice = (AliRun*)fBgrFile->Get("gAlice");
      ReadWriteDigit(fEvNrBgr);
    } //if merge


  f1->cd();
   
  //Make branch for digits
  TOF->MakeBranch("D");

  gime->TreeD()->Reset();
  gime->TreeD()->Fill();
  
  fDigits   = TOF->Digits();
  
  gime->WriteDigits("OVERWRITE");
  
  gAlice->ResetDigits();

}

//---------------------------------------------------------------------

void AliTOFMerger::ReadWriteDigit(Int_t iEvNum)
{
//
// Read Sdigits from the current file and write them in the TreeD
//
  AliTOFdigit* tofsdigit;
  
  AliTOF * tofinfile = (AliTOF *) gAlice->GetDetector("TOF") ;
  
  Int_t retval = fRunLoader->GetEvent(iEvNum);
  if (retval)
   {
     Error("ReadWriteDigit","Error while getting event %d",iEvNum);
     return;
   }

  AliLoader* gime = fRunLoader->GetLoader("TOFLoader");
  if (gime == 0x0)
   {
     Error("Exec","Can not find TOF loader in event. Exiting.");
     return;
   }
  
  
  
  if(gime->TreeS()==0) 
   {
    cout<<" TreeS==0 -> return"<<gime->TreeS()<<endl; 
    return ;
   }
  
  Int_t ndig, k;
  Int_t    tracks[3];    // track info
  Int_t    vol[5];       // location for a digit
  Float_t  digit[2];     // TOF digit variables

  gAlice->ResetDigits();
  gime->TreeS()->GetEvent(iEvNum);
  TClonesArray * TOFdigits   = tofinfile->SDigits();
  
  ndig=TOFdigits->GetEntries();

  for (k=0; k<ndig; k++) {
    tofsdigit= (AliTOFdigit*) TOFdigits->UncheckedAt(k);

	tracks[0] = tofsdigit->GetTrack(0);
	tracks[1] = tofsdigit->GetTrack(1);
	tracks[2] = tofsdigit->GetTrack(2);

	vol[0] = tofsdigit->GetSector();
	vol[1] = tofsdigit->GetPlate();
	vol[2] = tofsdigit->GetStrip();
	vol[3] = tofsdigit->GetPadx();
	vol[4] = tofsdigit->GetPadz();

	digit[0] = tofsdigit->GetTdc();
	digit[1] = tofsdigit->GetAdc();

	new ((*fDigits)[fNDigits++]) AliTOFdigit(tracks, vol, digit);
  } // end loop on sdigits in the current file
}



