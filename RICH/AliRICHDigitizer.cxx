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

/* 
   $Log$
   Revision 1.4  2002/05/28 07:53:10  morsch
   Wrong order of arguments in for-statement corrected.

   Revision 1.3  2001/12/05 14:53:34  hristov
   Destructor corrected

   Revision 1.2  2001/11/07 14:50:31  hristov
   Minor correction of the Log part

   Revision 1.1  2001/11/02 15:37:26  hristov
   Digitizer class created. Code cleaning and bug fixes (J.Chudoba)
*/
#include <iostream.h> 

#include <TTree.h> 
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TParticle.h>

#include "AliRICHDigitizer.h"
#include "AliRICHChamber.h"
#include "AliHitMap.h"
#include "AliRICHHitMapA1.h"
#include "AliRICH.h"
#include "AliRICHHit.h"
#include "AliRICHSDigit.h"
#include "AliRICHDigit.h"
#include "AliRICHTransientDigit.h"
#include "AliRun.h"
#include "AliPDG.h"
#include "AliRunDigitizer.h"

ClassImp(AliRICHDigitizer)

//___________________________________________
  AliRICHDigitizer::AliRICHDigitizer()
{
// Default constructor - don't use it   
  fHits = 0;
  fSDigits = 0;
  fHitMap = 0;
  fTDList = 0;
}

////////////////////////////////////////////////////////////////////////
AliRICHDigitizer::AliRICHDigitizer(AliRunDigitizer* manager) 
    :AliDigitizer(manager)
{
// ctor which should be used
  fHits = 0;
  fSDigits = 0;
  fHitMap = 0;
  fTDList = 0;
  fDebug   = 0; 
  if (GetDebug()>2) 
    cerr<<"AliRICHDigitizer::AliRICHDigitizer"
	<<"(AliRunDigitizer* manager) was processed"<<endl;
}

////////////////////////////////////////////////////////////////////////
AliRICHDigitizer::~AliRICHDigitizer()
{
// Destructor
    if (fHits) {
	fHits->Delete();
	delete fHits;
    }
    if (fSDigits) {
	fSDigits->Delete();
	delete fSDigits;
    }
    for (Int_t i=0; i<kNCH; i++ )
	delete fHitMap[i];
    delete [] fHitMap;
    if (fTDList) {
	fTDList->Delete();
	delete fTDList;
  }
}

//------------------------------------------------------------------------
Bool_t AliRICHDigitizer::Exists(const AliRICHSDigit *padhit)
{
  return (fHitMap[fNch]->TestHit(padhit->PadX(),padhit->PadY()));
}

//------------------------------------------------------------------------
void AliRICHDigitizer::Update(AliRICHSDigit *padhit)
{
  AliRICHTransientDigit *pdigit = 
    static_cast<AliRICHTransientDigit*>(
      fHitMap[fNch]->GetHit(padhit->PadX(),padhit->PadY()));

  // update charge
  //
  Int_t iqpad    = Int_t(padhit->QPad()); // charge per pad
  pdigit->AddSignal(iqpad);
  pdigit->AddPhysicsSignal(iqpad);		

  // update list of tracks
  //
  Int_t track, charge;
  track = fTrack+fMask;
  if (fSignal) {
    charge = iqpad;
  } else {
    charge = kBgTag;
  }
  pdigit->UpdateTrackList(track,charge);
}

//------------------------------------------------------------------------
void AliRICHDigitizer::CreateNew(AliRICHSDigit *padhit)
{
  fTDList->AddAtAndExpand(
    new AliRICHTransientDigit(fNch,fDigits),fCounter);
  fHitMap[fNch]->SetHit(padhit->PadX(),padhit->PadY(),fCounter);

  AliRICHTransientDigit* pdigit = 
    static_cast<AliRICHTransientDigit*>(fTDList->Last());
  // list of tracks
  Int_t track, charge;
  if (fSignal) {
    track = fTrack;
    charge = padhit->QPad();
  } else {
    track = kBgTag;
    charge = kBgTag;
  }
  pdigit->AddToTrackList(track,charge);
  fCounter++;
}


////////////////////////////////////////////////////////////////////////
Bool_t AliRICHDigitizer::Init()
{
// Initialisation
  fHits     = new TClonesArray("AliRICHHit",1000);
  fSDigits  = new TClonesArray("AliRICHSDigit",1000);
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
//void AliRICHDigitizer::Digitise(Int_t nev, Int_t flag)
void AliRICHDigitizer::Exec(Option_t* option)
{
  TString optionString = option;
  if (optionString.Data() == "deb") {
    cout<<"AliMUONDigitizer::Exec: called with option deb "<<endl;
    fDebug = 3;
  }

  AliRICHChamber*       iChamber;
  AliSegmentation*  segmentation;
  
  fTDList = new TObjArray;
    
  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
  pRICH->MakeBranchInTreeD(fManager->GetTreeD());
  fHitMap= new AliHitMap* [kNCH];
        
  for (Int_t i =0; i<kNCH; i++) {
    iChamber= &(pRICH->Chamber(i));
    segmentation=iChamber->GetSegmentationModel(1);
    fHitMap[i] = new AliRICHHitMapA1(segmentation, fTDList);
  }

// Loop over files to digitize
  fSignal = kTRUE;
  fCounter = 0;
  for (Int_t inputFile=0; inputFile<fManager->GetNinputs();
       inputFile++) {

// Connect RICH branches

    if (inputFile > 0 ) fSignal = kFALSE;
    TBranch *branchHits = 0;
    TBranch *branchSDigits = 0;
    TTree *treeH = fManager->GetInputTreeH(inputFile);
    if (GetDebug()>2) {
      cerr<<"   inputFile  "<<inputFile<<endl;
      cerr<<"   treeH, fHits "<<treeH<<" "<<fHits<<endl;
    }
    if (treeH && fHits) {
      branchHits = treeH->GetBranch("RICH");
      if (branchHits) {
	fHits->Clear();
	branchHits->SetAddress(&fHits);
      }
      else
	Error("Exec","branch RICH was not found");
    }
    if (GetDebug()>2) cerr<<"   branchHits = "<<branchHits<<endl;

    if (treeH && fSDigits) {
      branchSDigits = treeH->GetBranch("RICHSDigits");
      if (branchSDigits) 
	branchSDigits->SetAddress(&fSDigits);
      else
	Error("exec","branch RICHSDigits was not found");
    }
    if (GetDebug()>2) cerr<<"   branchSDigits = "<<branchSDigits<<endl;


    //
    //   Loop over tracks
    //
    
    Int_t ntracks =(Int_t) treeH->GetEntries();
    for (fTrack=0; fTrack<ntracks; fTrack++) {
      fHits->Clear();
      fSDigits->Clear();  
      branchHits->GetEntry(fTrack);
      branchSDigits->GetEntry(fTrack);


      //
      //   Loop over hits
      for(Int_t i = 0; i < fHits->GetEntriesFast(); ++i) {
	AliRICHHit* mHit = static_cast<AliRICHHit*>(fHits->At(i));
	fNch = mHit->Chamber()-1;  // chamber number
	if (fNch >= kNCH) {
	  cerr<<"AliRICHDigitizer: chamber nr. fNch out of range: "<<fNch<<endl;
	  cerr<<"               track: "<<fTrack<<endl;
	  continue;
	}
	iChamber = &(pRICH->Chamber(fNch));
	  
	//
	// Loop over pad hits
	for (AliRICHSDigit* mPad=
	       (AliRICHSDigit*)pRICH->FirstPad(mHit,fSDigits);
	     mPad;
	     mPad=(AliRICHSDigit*)pRICH->NextPad(fSDigits))
	{
	  Int_t iqpad    = mPad->QPad();       // charge per pad
	  fDigits[0]=mPad->PadX();
	  fDigits[1]=mPad->PadY();
	  fDigits[2]=iqpad;
	  fDigits[3]=iqpad;
	  fDigits[4]=mPad->HitNumber();
		  
	  // build the list of fired pads and update the info
	  if (Exists(mPad)) {
	    Update(mPad);
	  } else {
	    CreateNew(mPad);
	  }
	} //end loop over clusters
      } // hit loop
    } // track loop
  } // end file loop
  if (GetDebug()>2) cerr<<"END OF FILE LOOP"<<endl;

  Int_t tracks[kMAXTRACKSPERRICHDIGIT];
  Int_t charges[kMAXTRACKSPERRICHDIGIT];
  Int_t nentries=fTDList->GetEntriesFast();
    
  for (Int_t nent=0;nent<nentries;nent++) {
    AliRICHTransientDigit *transDigit=(AliRICHTransientDigit*)fTDList->At(nent);
    if (transDigit==0) continue; 
    Int_t ich=transDigit->GetChamber();
    Int_t q=transDigit->Signal(); 
    iChamber=&(pRICH->Chamber(ich));
    AliRICHResponse * response=iChamber->GetResponseModel();
    Int_t adcmax= (Int_t) response->MaxAdc();
      
      
    // add white noise and do zero-suppression and signal truncation (new electronics,old electronics gaus 1.2,0.2)
    //printf("Treshold: %d\n",iChamber->fTresh->GetHitIndex(transDigit->PadX(),transDigit->PadY()));

// tmp change
//    Int_t pedestal = iChamber->fTresh->GetHitIndex(transDigit->PadX(),transDigit->PadY());
    Int_t pedestal = 0;


    //printf("Pedestal:%d\n",pedestal);
    //Int_t pedestal=0;
    Float_t treshold = (pedestal + 4*2.2);
      
    Float_t meanNoise = gRandom->Gaus(2.2, 0.3);
    Float_t noise     = gRandom->Gaus(0, meanNoise);
    q+=(Int_t)(noise + pedestal);
    //q+=(Int_t)(noise);
    //          magic number to be parametrised !!! 
    if ( q <= treshold) 
    {
      q = q - pedestal;
      continue;
    }
    q = q - pedestal;
    if ( q >= adcmax) q=adcmax;
    fDigits[0]=transDigit->PadX();
    fDigits[1]=transDigit->PadY();
    fDigits[2]=q;
    fDigits[3]=transDigit->Physics();
    fDigits[4]=transDigit->Hit();

    Int_t nptracks = transDigit->GetNTracks();  

    // this was changed to accomodate the real number of tracks
    if (nptracks > kMAXTRACKSPERRICHDIGIT) {
      printf("Attention - tracks > 10 %d\n",nptracks);
      nptracks=kMAXTRACKSPERRICHDIGIT;
    }
    if (nptracks > 2) {
      printf("Attention - tracks > 2  %d \n",nptracks);
    }
    for (Int_t tr=0;tr<nptracks;tr++) {
      tracks[tr]=transDigit->GetTrack(tr);
      charges[tr]=transDigit->GetCharge(tr);
    }      //end loop over list of tracks for one pad
    if (nptracks < kMAXTRACKSPERRICHDIGIT ) {
      for (Int_t t=nptracks; t<kMAXTRACKSPERRICHDIGIT; t++) {
	tracks[t]=0;
	charges[t]=0;
      }
    }
    //write file
    //if (ich==2)
    //fprintf(points,"%4d,      %4d,      %4d\n",digits[0],digits[1],digits[2]);
      
    // fill digits
    pRICH->AddDigits(ich,tracks,charges,fDigits);
  }	
  fManager->GetTreeD()->Fill();

  pRICH->ResetDigits();
  fTDList->Delete();    // or fTDList->Clear(); ???
  for(Int_t ii=0;ii<kNCH;++ii) {
    if (fHitMap[ii]) {
      delete fHitMap[ii];
      fHitMap[ii]=0;
    }
  }
    
  TClonesArray *richDigits;
  for (Int_t k=0;k<kNCH;k++) {
    richDigits = pRICH->DigitsAddress(k);
    Int_t ndigit=richDigits->GetEntriesFast();
    printf ("Chamber %d digits %d \n",k,ndigit);
  }
  pRICH->ResetDigits(); /// ??? should it be here???
  fManager->GetTreeD()->Write(0,TObject::kOverwrite);
  delete [] fHitMap;
  delete fTDList;

  if (fHits)    fHits->Delete();
  if (fSDigits) fSDigits->Delete();

}
