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
Revision 1.1  2001/02/27 22:13:34  jbarbosa
Implementing merger class.

*/

#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TParticle.h>


// #include "AliMerger.h"
// #include "AliMergable.h"
#include "AliRICHMerger.h"
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
    fTrList     = 0;
    fAddress    = 0; 
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
    if (fTrList)     delete fTrList;
    if (fAddress)    delete fAddress; 
}

//------------------------------------------------------------------------
Bool_t AliRICHMerger::Exists(const AliRICHSDigit *mergable)
{
    AliRICHSDigit *padhit = (AliRICHSDigit*) mergable;
    return (fHitMap[fNch]->TestHit(padhit->PadX(),padhit->PadY()));
}

//------------------------------------------------------------------------
void AliRICHMerger::Update(AliRICHSDigit *mergable)
{
    AliRICHSDigit *padhit = (AliRICHSDigit*) mergable;    
    AliRICHTransientDigit* pdigit;
    Int_t ipx      = padhit->PadX();        // pad number on X
    Int_t ipy      = padhit->PadY();        // pad number on Y
    Int_t iqpad    = Int_t(padhit->QPad()); // charge per pad

    pdigit = (AliRICHTransientDigit*) fHitMap[fNch]->GetHit(ipx, ipy);
    // update charge
    //
    (*pdigit).AddSignal(iqpad);
    (*pdigit).AddPhysicsSignal(iqpad);		
    // update list of tracks
    //
    TObjArray* fTrList = (TObjArray*)pdigit->TrackList();
    Int_t lastEntry = fTrList->GetLast();
    TVector *pTrack = (TVector*)fTrList->At(lastEntry);
    TVector &ptrk   = *pTrack;
    TVector &trinfo = *((TVector*) (*fAddress)[fCountadr-1]);
    Int_t lastTrack = Int_t(ptrk(0));

    if (trinfo(0) != kBgTag) {
	if (lastTrack == fTrack) {
	    Int_t lastCharge = Int_t(ptrk(1));
	    lastCharge += iqpad;
	    fTrList->RemoveAt(lastEntry);
	    trinfo(1) = lastCharge;
	    fTrList->AddAt(&trinfo,lastEntry);
	} else {
	    fTrList->Add(&trinfo);
	}
    } else {
	if (lastTrack != -1) fTrList->Add(&trinfo);
    }
}

//------------------------------------------------------------------------
void AliRICHMerger::CreateNew(AliRICHSDigit *mergable)
{
    AliRICHSDigit *padhit = (AliRICHSDigit*) mergable;    
    AliRICHTransientDigit* pdigit;

    Int_t ipx      = padhit->PadX();       // pad number on X
    Int_t ipy      = padhit->PadY();       // pad number on Y
    fList->AddAtAndExpand(
	new AliRICHTransientDigit(fNch,fDigits),fCounter);
    fHitMap[fNch]->SetHit(ipx, ipy, fCounter);
    fCounter++;
    pdigit = (AliRICHTransientDigit*)fList->At(fList->GetLast());
    // list of tracks
    TObjArray *fTrList = (TObjArray*)pdigit->TrackList();
    TVector &trinfo    = *((TVector*) (*fAddress)[fCountadr-1]);
    fTrList->Add(&trinfo);
}


//------------------------------------------------------------------------
void AliRICHMerger::Init()
{
// Initialisation
    
    if (fMerge) fBgrFile = InitBgr();
}



//------------------------------------------------------------------------
TFile* AliRICHMerger::InitBgr()
{
// Initialise background event
    TFile *file = new TFile(fFnBgr);
// add error checking later
    printf("\n AliRICHMerger has opened %s file with background event \n", fFnBgr);
    fHitsBgr     = new TClonesArray("AliRICHHit",1000);
    fSDigitsBgr  = new TClonesArray("AliRICHSDigit",1000);
    return file;
}

//------------------------------------------------------------------------
void AliRICHMerger::Digitise(Int_t nev, Int_t flag)
{

 // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !

  Int_t particle;

  //FILE* points; //these will be the digits...
  
  //points=fopen("points.dat","w");
  
  AliRICHChamber*       iChamber;
  AliSegmentation*  segmentation;
  
    Int_t digitise=0;
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    static TClonesArray *pAddress=0;
    if(!pAddress) pAddress=new TClonesArray("TVector",1000);
    Int_t digits[5]; 
    
    AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
    AliHitMap* pHitMap[10];
    Int_t i;
    for (i=0; i<10; i++) {pHitMap[i]=0;}
    
    if (fMerge ) {
      fBgrFile->cd();
      //fBgrFile->ls();
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
    
    AliHitMap* hm;
    Int_t countadr=0;
    Int_t counter=0;
    for (i =0; i<kNCH; i++) {
      iChamber= &(pRICH->Chamber(i));
      segmentation=iChamber->GetSegmentationModel(1);
      pHitMap[i] = new AliRICHHitMapA1(segmentation, list);
    }
    //
    //   Loop over tracks
    //
    
    TTree *treeH = gAlice->TreeH();
    Int_t ntracks =(Int_t) treeH->GetEntries();
    for (Int_t track=0; track<ntracks; track++) {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      //
      //   Loop over hits
      for(AliRICHHit* mHit=(AliRICHHit*)pRICH->FirstHit(-1); 
	  mHit;
	  mHit=(AliRICHHit*)pRICH->NextHit()) 
	{
	  
	  Int_t   nch   = mHit->fChamber-1;  // chamber number
	  Int_t   index = mHit->Track();
	  if (nch >kNCH) continue;
	  iChamber = &(pRICH->Chamber(nch));
	  
	  TParticle *current = (TParticle*)gAlice->Particle(index);
	  
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
	  
	  digitise=1;
	  
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
		  Int_t ipx      = mPad->fPadX;       // pad number on X
		  Int_t ipy      = mPad->fPadY;       // pad number on Y
		  Int_t iqpad    = mPad->fQpad;       // charge per pad
		  //
		  //
		  //printf("X:%d, Y:%d, Q:%d\n",ipx,ipy,iqpad);
		  
		  Float_t thex, they, thez;
		  segmentation=iChamber->GetSegmentationModel(0);
		  segmentation->GetPadC(ipx,ipy,thex,they,thez);
		  new((*pAddress)[countadr++]) TVector(2);
		  TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
		  trinfo(0)=(Float_t)track;
		  trinfo(1)=(Float_t)iqpad;
		  
		  digits[0]=ipx;
		  digits[1]=ipy;
		  digits[2]=iqpad;
		  
		  AliRICHTransientDigit* pdigit;
		  // build the list of fired pads and update the info
		  if (!pHitMap[nch]->TestHit(ipx, ipy)) {
		    list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
		    pHitMap[nch]->SetHit(ipx, ipy, counter);
		    counter++;
		    pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
		    // list of tracks
		    TObjArray *trlist=(TObjArray*)pdigit->TrackList();
		    trlist->Add(&trinfo);
		  } else {
		    pdigit=(AliRICHTransientDigit*) pHitMap[nch]->GetHit(ipx, ipy);
		    // update charge
		    (*pdigit).fSignal+=iqpad;
		    // update list of tracks
		    TObjArray* trlist=(TObjArray*)pdigit->TrackList();
		    Int_t lastEntry=trlist->GetLast();
		    TVector *ptrkP=(TVector*)trlist->At(lastEntry);
		    TVector &ptrk=*ptrkP;
		    Int_t lastTrack=Int_t(ptrk(0));
		    Int_t lastCharge=Int_t(ptrk(1));
		    if (lastTrack==track) {
		      lastCharge+=iqpad;
		      trlist->RemoveAt(lastEntry);
		      trinfo(0)=lastTrack;
		      trinfo(1)=lastCharge;
		      trlist->AddAt(&trinfo,lastEntry);
		    } else {
		      trlist->Add(&trinfo);
		    }
		    // check the track list
		    Int_t nptracks=trlist->GetEntriesFast();
		    if (nptracks > 2) {
		      printf("Attention - tracks:  %d (>2)\n",nptracks);
		      //printf("cat,nch,ix,iy %d %d %d %d  \n",icat+1,nch,ipx,ipy);
		      for (Int_t tr=0;tr<nptracks;tr++) {
			TVector *pptrkP=(TVector*)trlist->At(tr);
			TVector &pptrk=*pptrkP;
			trk[tr]=Int_t(pptrk(0));
			chtrk[tr]=Int_t(pptrk(1));
		      }
		    } // end if nptracks
		  } //  end if pdigit
		} //end loop over clusters
	    }// track type condition
	} // hit loop
    } // track loop
    
    // open the file with background
    
    if (fMerge) {
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
	AliRICHHit* mHit;
	for(int i = 0; i < fHitsBgr->GetEntriesFast(); ++i) 
	  {	
	    mHit   = (AliRICHHit*) (*fHitsBgr)[i];
	    fNch   = mHit->Chamber()-1;  // chamber number
	    iChamber = &(pRICH->Chamber(fNch));
	    //Float_t xbgr = mHit->X();
	    //Float_t ybgr = mHit->Y();
	    Bool_t cond  = kFALSE;
	    cond  = kTRUE;
	    //
	    // Loop over pad hits
	    for (AliRICHSDigit* mPad =
		   (AliRICHSDigit*)pRICH->FirstPad(mHit,fSDigitsBgr);
		 mPad;
		 mPad = (AliRICHSDigit*)pRICH->NextPad(fSDigitsBgr))
	      {
		Int_t ipx      = mPad->PadX();        // pad number on X
		Int_t ipy      = mPad->PadY();        // pad number on Y
		Int_t iqpad    = Int_t(mPad->QPad()); // charge per pad
		
		new((*fAddress)[fCountadr++]) TVector(2);
		TVector &trinfo = *((TVector*) (*fAddress)[fCountadr-1]);
		trinfo(0) = kBgTag;  // tag background
		trinfo(1) = kBgTag;
		
		fDigits[0] = ipx;
		fDigits[1] = ipy;
		fDigits[3] = iqpad;
		
		// build the list of fired pads and update the info
		if (!Exists(mPad)) {
		  CreateNew(mPad);
		} else {
		  Update(mPad);
		} //  end if !Exists
	      } //end loop over clusters
	  } // hit loop
      } // track loop
      
      TTree *fAli = gAlice->TreeK();
      TFile *file = NULL;
      
      if (fAli) file = fAli->GetCurrentFile();
      file->cd();
    } // if fMerge
    
    Int_t tracks[10];
    Int_t charges[10];
    //cout<<"Start filling digits \n "<<endl;
    Int_t nentries=list->GetEntriesFast();
    //printf(" \n \n nentries %d \n",nentries);
    
    // start filling the digits
    
    for (Int_t nent=0;nent<nentries;nent++) {
      AliRICHTransientDigit *address=(AliRICHTransientDigit*)list->At(nent);
      if (address==0) continue; 
      
      Int_t ich=address->fChamber;
      Int_t q=address->fSignal; 
      iChamber=&(pRICH->Chamber(ich));
      AliRICHResponse * response=iChamber->GetResponseModel();
      Int_t adcmax= (Int_t) response->MaxAdc();
      
      
      // add white noise and do zero-suppression and signal truncation (new electronics,old electronics gaus 1.2,0.2)
      //printf("Treshold: %d\n",iChamber->fTresh->GetHitIndex(address->fPadX,address->fPadY));
      Int_t pedestal = iChamber->fTresh->GetHitIndex(address->fPadX,address->fPadY);

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
      digits[0]=address->fPadX;
      digits[1]=address->fPadY;
      digits[2]=q;
      
      TObjArray* trlist=(TObjArray*)address->TrackList();
      Int_t nptracks=trlist->GetEntriesFast();
      
      // this was changed to accomodate the real number of tracks
      if (nptracks > 10) {
	printf("Attention - tracks > 10 %d\n",nptracks);
	nptracks=10;
      }
      if (nptracks > 2) {
	printf("Attention - tracks > 2  %d \n",nptracks);
	//printf("cat,ich,ix,iy,q %d %d %d %d %d \n",
	//icat,ich,digits[0],digits[1],q);
      }
      for (Int_t tr=0;tr<nptracks;tr++) {
	TVector *ppP=(TVector*)trlist->At(tr);
	TVector &pp  =*ppP;
	tracks[tr]=Int_t(pp(0));
	charges[tr]=Int_t(pp(1));
      }      //end loop over list of tracks for one pad
      if (nptracks < 10 ) {
	for (Int_t t=nptracks; t<10; t++) {
	  tracks[t]=0;
	  charges[t]=0;
	}
      }
      //write file
      //if (ich==2)
	//fprintf(points,"%4d,      %4d,      %4d\n",digits[0],digits[1],digits[2]);
      
      // fill digits
      pRICH->AddDigits(ich,tracks,charges,digits);
    }	
    gAlice->TreeD()->Fill();
    
    list->Delete();
    for(Int_t ii=0;ii<kNCH;++ii) {
      if (pHitMap[ii]) {
	hm=pHitMap[ii];
	delete hm;
	pHitMap[ii]=0;
      }
    }
    
    //TTree *TD=gAlice->TreeD();
    //Stat_t ndig=TD->GetEntries();
    //cout<<"number of digits  "<<ndig<<endl;
    TClonesArray *fDch;
    for (int k=0;k<kNCH;k++) {
      fDch= pRICH->DigitsAddress(k);
      int ndigit=fDch->GetEntriesFast();
      printf ("Chamber %d digits %d \n",k,ndigit);
    }
    pRICH->ResetDigits();
    char hname[30];
    sprintf(hname,"TreeD%d",nev);
    gAlice->TreeD()->Write(hname);
    
    // reset tree
    //    gAlice->TreeD()->Reset();
    delete list;
    pAddress->Clear();
    // gObjectTable->Print();

}





