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

/* $Id$ */

#include <Riostream.h> 

#include <TTree.h> 
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TParticle.h>

#include "AliRICHMerger.h"
#include "AliRICHChamber.h"
#include "AliHitMap.h"
#include "AliRICHHitMapA1.h"
#include "AliRICH.h"
#include "AliRICHSDigit.h"
#include "AliRICHDigit.h"
#include "AliRICHTransientDigit.h"
#include "AliRun.h"
#include "AliPDG.h"

ClassImp(AliRICHMerger)

//___________________________________________
  AliRICHMerger::AliRICHMerger()
{
// Default constructor    
  fEvNrSig = 0;
  fEvNrBgr = 0;
  fMerge   = kDigitize;
  fFnBgr   = 0;
  fHitMap  = 0;
  fList    = 0;
  fTrH1    = 0;
  fHitsBgr = 0;
  fSDigitsBgr = 0;
  fHitMap     = 0;
  fList       = 0;
  fBgrFile    = 0;
}

//------------------------------------------------------------------------
AliRICHMerger::~AliRICHMerger()
{
// Destructor
  if (fTrH1)       delete fTrH1;
  if (fHitsBgr)    delete fHitsBgr;
  if (fSDigitsBgr) delete fSDigitsBgr;
  if (fHitMap)     delete fHitMap;
  if (fList)       delete fList;
  if (fBgrFile)    delete fBgrFile;
}

//------------------------------------------------------------------------
Bool_t AliRICHMerger::Exists(const AliRICHSDigit *padhit)
{
  return (fHitMap[fNch]->TestHit(padhit->PadX(),padhit->PadY()));
}

//------------------------------------------------------------------------
void AliRICHMerger::Update(AliRICHSDigit *padhit)
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
  if (fSignal) {
    track = fTrack;
    charge = iqpad;
  } else {
    track = kBgTag;
    charge = kBgTag;
  }
  pdigit->UpdateTrackList(track,charge);
}

//------------------------------------------------------------------------
void AliRICHMerger::CreateNew(AliRICHSDigit *padhit)
{
  fList->AddAtAndExpand(
    new AliRICHTransientDigit(fNch,fDigits),fCounter);
  fHitMap[fNch]->SetHit(padhit->PadX(),padhit->PadY(),fCounter);

  AliRICHTransientDigit* pdigit = 
    static_cast<AliRICHTransientDigit*>(fList->Last());
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


//------------------------------------------------------------------------
void AliRICHMerger::Init()
{
// Initialisation
    
  if (fMerge && !fBgrFile) fBgrFile = InitBgr();
}



//------------------------------------------------------------------------
TFile* AliRICHMerger::InitBgr()
{
// Initialise background event
  TFile *file = new TFile(fFnBgr);
// add error checking later
  printf("\n AliRICHMerger has opened %s file with background event \n", fFnBgr);
  fHitsBgr     = new TClonesArray("AliRICHhit",1000);
  fSDigitsBgr  = new TClonesArray("AliRICHSDigit",1000);
  return file;
}

//------------------------------------------------------------------------
void AliRICHMerger::Digitise(Int_t nev, Int_t flag)
{

  // keep galice.root for signal and name differently the file for 
  // background when add! otherwise the track info for signal will be lost !

  Int_t particle;

  AliRICHChamber*       iChamber;
  AliSegmentation*  segmentation;
  
  fList = new TObjArray;
    
  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
  fHitMap= new AliHitMap* [kNCH];
    
  if (fMerge ) {
    fBgrFile->cd();
    //
    // Get Hits Tree header from file
    if(fHitsBgr) fHitsBgr->Clear();
    if(fSDigitsBgr) fSDigitsBgr->Clear();
    if(fTrH1) delete fTrH1;
    fTrH1 = 0;
      
    char treeName[20];
    sprintf(treeName,"TreeH%d",fEvNrBgr);
    fTrH1 = (TTree*)gDirectory->Get(treeName);
    if (!fTrH1) {
      printf("ERROR: cannot find Hits Tree for event:%d\n",fEvNrBgr);
    }
    //
    // Set branch addresses
    TBranch *branch;
    char branchname[20];
    sprintf(branchname,"%s",pRICH->GetName());
    if (fTrH1 && fHitsBgr) {
      branch = fTrH1->GetBranch(branchname);
      if (branch) branch->SetAddress(&fHitsBgr);
    }
    if (fTrH1 && fSDigitsBgr) {
      branch = fTrH1->GetBranch("RICHSDigits");
      if (branch) branch->SetAddress(&fSDigitsBgr);
    }
  }
    
  for (Int_t i =0; i<kNCH; i++) {
    iChamber= &(pRICH->Chamber(i));
    segmentation=iChamber->GetSegmentationModel(1);
    fHitMap[i] = new AliRICHHitMapA1(segmentation, fList);
  }
  //
  //   Loop over tracks
  //
    
  fSignal = kTRUE;
  fCounter = 0;
  TTree *treeH = pRICH->TreeH();
  Int_t ntracks =(Int_t) treeH->GetEntries();
  for (fTrack=0; fTrack<ntracks; fTrack++) {
    gAlice->ResetHits();
    treeH->GetEvent(fTrack);
    //
    //   Loop over hits
    for(AliRICHhit* mHit=(AliRICHhit*)pRICH->FirstHit(-1); 
	mHit;
	mHit=(AliRICHhit*)pRICH->NextHit()) 
    {
	  
      fNch = mHit->Chamber()-1;  // chamber number
      Int_t trackIndex = mHit->Track();
      if (fNch >= kNCH) {
	cerr<<"AliRICHMerger: chamber nr. fNch out of range: "<<fNch<<endl;
	cerr<<"               track: "<<fTrack<<endl;
	continue;
      }
      iChamber = &(pRICH->Chamber(fNch));
	  

// 
// If flag is set, create digits only for some particles
//
      TParticle *current = (TParticle*)gAlice->Particle(trackIndex);
      if (current->GetPdgCode() >= 50000050)
      {
	TParticle *motherofcurrent = (TParticle*)gAlice->Particle(current->GetFirstMother());
	particle = motherofcurrent->GetPdgCode();
      }
      else
      {
	particle = current->GetPdgCode();
      }

	  
      //printf("Flag:%d\n",flag);
      //printf("Track:%d\n",track);
      //printf("Particle:%d\n",particle);
	  
      Int_t digitise=1;
	  
      if (flag == 1) 
	if(TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
	  digitise=0;
	  
      if (flag == 2)
	if(TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
	   || TMath::Abs(particle)==311)
	  digitise=0;
	  
      if (flag == 3 && TMath::Abs(particle)==2212)
	digitise=0;
	  
      if (flag == 4 && TMath::Abs(particle)==13)
	digitise=0;
	  
      if (flag == 5 && TMath::Abs(particle)==11)
	digitise=0;
	  
      if (flag == 6 && TMath::Abs(particle)==2112)
	digitise=0;
	  
	  
      //printf ("Particle: %d, Flag: %d, Digitise: %d\n",particle,flag,digitise); 
	  
	  
      if (digitise)
      {
	      
	//
	// Loop over pad hits
	for (AliRICHSDigit* mPad=
	       (AliRICHSDigit*)pRICH->FirstPad(mHit,pRICH->SDigits());
	     mPad;
	     mPad=(AliRICHSDigit*)pRICH->NextPad(pRICH->SDigits()))
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
      }// track type condition
    } // hit loop
  } // track loop
    
  // open the file with background
    
  if (fMerge) {
    fSignal = kFALSE;
    ntracks = (Int_t)fTrH1->GetEntries();
    //
    //   Loop over tracks
    //
    for (fTrack = 0; fTrack < ntracks; fTrack++) {
	
      if (fHitsBgr)       fHitsBgr->Clear();
      if (fSDigitsBgr)    fSDigitsBgr->Clear();
	
      fTrH1->GetEvent(fTrack);
      //
      //   Loop over hits
      AliRICHhit* mHit;
      for(Int_t i = 0; i < fHitsBgr->GetEntriesFast(); ++i) 
      {	
	mHit   = (AliRICHhit*) (*fHitsBgr)[i];
	fNch   = mHit->Chamber()-1;  // chamber number
	iChamber = &(pRICH->Chamber(fNch));

	//
	// Loop over pad hits
	for (AliRICHSDigit* mPad =
	       (AliRICHSDigit*)pRICH->FirstPad(mHit,fSDigitsBgr);
	     mPad;
	     mPad = (AliRICHSDigit*)pRICH->NextPad(fSDigitsBgr))
	{
	  fDigits[0] = mPad->PadX(); 
	  fDigits[1] = mPad->PadY();  
	  fDigits[2] = mPad->QPad();
	  fDigits[3] = 0;
	  fDigits[4] = -1; // set hit number to -1 for bgr
		
	  // build the list of fired pads and update the info
	  if (!Exists(mPad)) {
	    CreateNew(mPad);
	  } else {
	    Update(mPad);
	  } //  end if !Exists
	} //end loop over clusters
      } // hit loop
    } // track loop
      
    TTree *treeK = gAlice->TreeK();
    TFile *file = NULL;
      
    if (treeK) file = treeK->GetCurrentFile();
    file->cd();
  } // if fMerge
    
  Int_t tracks[kMAXTRACKSPERRICHDIGIT];
  Int_t charges[kMAXTRACKSPERRICHDIGIT];
  Int_t nentries=fList->GetEntriesFast();
    
  for (Int_t nent=0;nent<nentries;nent++) {
    AliRICHTransientDigit *transDigit=(AliRICHTransientDigit*)fList->At(nent);
    if (transDigit==0) continue; 
    Int_t ich=transDigit->GetChamber();
    Int_t q=transDigit->Signal(); 
    iChamber=&(pRICH->Chamber(ich));
    AliRICHResponse * response=iChamber->GetResponseModel();
    Int_t adcmax= (Int_t) response->MaxAdc();
      
      
    // add white noise and do zero-suppression and signal truncation (new electronics,old electronics gaus 1.2,0.2)
    //printf("Treshold: %d\n",iChamber->fTresh->GetHitIndex(transDigit->PadX(),transDigit->PadY()));
    Int_t pedestal = iChamber->fTresh->GetHitIndex(transDigit->PadX(),transDigit->PadY());

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
      //printf("cat,ich,ix,iy,q %d %d %d %d %d \n",
      //icat,ich,digits[0],digits[1],q);
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
  gAlice->TreeD()->Fill();

  fList->Delete();
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
  gAlice->TreeD()->Write(0,TObject::kOverwrite);
  // reset tree
  //    gAlice->TreeD()->Reset();
  //    delete fList; // deleted in dtor
  // gObjectTable->Print();

}
