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

#include <Riostream.h> 

#include <TTree.h> 
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TParticle.h>

#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliRICHDigitizer.h"
#include "AliRICHChamber.h"
#include "AliHitMap.h"
#include "AliRICHHitMapA1.h"
#include "AliRICH.h"
#include "AliRICHSDigit.h"
#include "AliRICHDigit.h"
#include "AliRICHTransientDigit.h"
#include "AliRun.h"
#include "AliPDG.h"
#include "AliRunDigitizer.h"

ClassImp(AliRICHDigitizer)

//__________________________________________________________________________________________________
AliRICHDigitizer::AliRICHDigitizer()
{//default ctor - don't use it   
  fHits = 0;
  fSDigits = 0;
  fHitMap = 0;
  fTDList = 0;
}//default ctor
//__________________________________________________________________________________________________
AliRICHDigitizer::AliRICHDigitizer(AliRunDigitizer *pManager) 
                 :AliDigitizer(pManager)
{//main ctor which should be used
  if(fManager->GetDebug())Info("main ctor","Start.");
  fHits = 0;
  fSDigits = 0;
  fHitMap = 0;
  fTDList = 0;
  fDebug  = pManager->GetDebug(); 
}//main ctor
//__________________________________________________________________________________________________
AliRICHDigitizer::~AliRICHDigitizer()
{//dtor
  if(fManager->GetDebug())Info("dtor","Start.");
  
  if(fHits)    {fHits->Delete();   delete fHits;}
  if(fSDigits) {fSDigits->Delete();delete fSDigits;}
  if(fTDList)  {fTDList->Delete(); delete fTDList;}
  for(Int_t i=0; i<kNCH; i++) delete fHitMap[i];    delete [] fHitMap;
}//dtor
//__________________________________________________________________________________________________
Bool_t AliRICHDigitizer::Exists(const AliRICHSDigit *p) 
{
  return (fHitMap[fNch]->TestHit(p->PadX(),p->PadY()));
}//Exists
//__________________________________________________________________________________________________
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
}//Update()
//__________________________________________________________________________________________________
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
}//CreateNew()
//__________________________________________________________________________________________________
Bool_t AliRICHDigitizer::Init()
{// Initialisation
  if(GetDebug())Info("Init","Start.");
  fHits     = new TClonesArray("AliRICHhit",1000);
  fSDigits  = new TClonesArray("AliRICHSDigit",1000);
  return kTRUE;
}//Init()
//__________________________________________________________________________________________________
void AliRICHDigitizer::Exec(Option_t* option)
{
  if(GetDebug())Info("Exec","Start with option=%s",option);

  AliRICHChamber*       iChamber;
  AliSegmentation*  segmentation;
 
  AliRunLoader *pInAL, *pOutAL;//in and out Run loaders
  AliLoader    *pInRL, *pOutRL;//in and out RICH loaders
 
  pOutAL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  pOutRL = pOutAL->GetLoader("RICHLoader");

  fTDList = new TObjArray;
  
     
  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
  
  if(!pOutRL->TreeD()) pOutRL->MakeTree("D");  pRICH->MakeBranch("D");
  
  
  fHitMap= new AliHitMap* [kNCH];
        
  for (Int_t i =0; i<kNCH; i++) {
    iChamber= &(pRICH->Chamber(i));
    segmentation=iChamber->GetSegmentationModel();
    fHitMap[i] = new AliRICHHitMapA1(segmentation, fTDList);
  }

// Loop over files to digitize
  fSignal = kTRUE;
  fCounter = 0;
  for (Int_t inputFile=0;inputFile<fManager->GetNinputs();inputFile++){//files loop
    pInAL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    pInRL = pInAL->GetLoader("RICHLoader");
    
// Connect RICH branches

    if (inputFile > 0 ) fSignal = kFALSE;
    TBranch *branchHits = 0;
    TBranch *branchSDigits = 0;

    pInRL->LoadHits("READ");
    
    TTree *treeH = pInRL->TreeH();
    if(pInRL->TreeH()) {
      branchHits = treeH->GetBranch("RICH");
      if(branchHits){
        fHits->Clear();
        branchHits->SetAddress(&fHits);}
      else
        Error("Exec","branch RICH was not found");
    }

    if(treeH) {
      branchSDigits = treeH->GetBranch("RICHSpecials");
      if(branchSDigits) 
        branchSDigits->SetAddress(&fSDigits);
      else
        Error("exec","branch RICHSpecials was not found");
    }
    
    for (fTrack=0; fTrack<treeH->GetEntries(); fTrack++) {//prims loop
      fHits->Clear();                   fSDigits->Clear();  
      fHits->Print();

      branchHits->GetEntry(fTrack);     branchSDigits->GetEntry(fTrack);
      for(Int_t i=0;i<fHits->GetEntriesFast();++i){//hits loop
	AliRICHhit* pHit = static_cast<AliRICHhit*>(fHits->At(i));
	fNch = pHit->Chamber()-1;  // chamber number
	if (fNch >= kNCH) {
	  cerr<<"AliRICHDigitizer: chamber nr. fNch out of range: "<<fNch<<endl;
	  cerr<<"               track: "<<fTrack<<endl;
	  continue;
	}
	iChamber = &(pRICH->Chamber(fNch));
	  
	for (AliRICHSDigit* mPad=FirstPad(pHit,fSDigits);mPad;mPad=NextPad(fSDigits)){//specials loop
	  Int_t iqpad    = mPad->QPad();       // charge per pad
	  fDigits[0]=mPad->PadX();
	  fDigits[1]=mPad->PadY();
	  fDigits[2]=iqpad;
	  fDigits[3]=iqpad;
	  fDigits[4]=mPad->HitNumber();
	  if (Exists(mPad))	  // build the list of fired pads and update the info
            Update(mPad);
	  else
            CreateNew(mPad);
	}//specials loop
      }//hits loop
    }//prims loop
  }//files loop
  if(GetDebug()) Info("Exec","END OF FILE LOOP");

  Int_t tracks[kMAXTRACKSPERRICHDIGIT];
  Int_t charges[kMAXTRACKSPERRICHDIGIT];
  Int_t nentries=fTDList->GetEntriesFast();
    
  for (Int_t nent=0;nent<nentries;nent++){//transient digits loop
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

    
    if(nptracks>kMAXTRACKSPERRICHDIGIT){// this was changed to accomodate the real number of tracks
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
    pRICH->AddDigitOld(ich+1,tracks,charges,fDigits);//OLD
  }//transient digits loop      
  pOutRL->TreeD()->Fill();

  //pRICH->ResetDigits();
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
    printf ("Chamber %d digits %d \n",k+1,ndigit);
  }
  pRICH->ResetDigits(); /// ??? should it be here???
  
  pOutRL->WriteDigits("OVERWRITE");

  delete [] fHitMap;
  delete fTDList;
  pOutRL->UnloadHits();
  pOutRL->UnloadDigits();
//  if (fHits)    fHits->Delete();
//  if (fSDigits) fSDigits->Delete();
  if(GetDebug())Info("Exec","Stop.");
}//Exec()
//__________________________________________________________________________________________________

static Int_t sMaxIterPad=0;    // Static variables for the pad-hit iterator routines
static Int_t sCurIterPad=0;

//__________________________________________________________________________________________________
AliRICHSDigit* AliRICHDigitizer::FirstPad(AliRICHhit*  hit,TClonesArray *clusters ) 
{// Initialise the pad iterator Return the address of the first sdigit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->PHlast() > 0) {
	sMaxIterPad=Int_t(hit->PHlast());
	sCurIterPad=Int_t(hit->PHfirst());
	return (AliRICHSDigit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
    
}//FirstPad
//__________________________________________________________________________________________________
AliRICHSDigit* AliRICHDigitizer::NextPad(TClonesArray *clusters) 
{// Iterates over pads
  
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliRICHSDigit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}//NextPad
//__________________________________________________________________________________________________
