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
*/

#include <TTree.h> 
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>


// #include "AliMerger.h"
// #include "AliMergable.h"
#include "AliMUONMerger.h"
#include "AliMUONConstants.h"
#include "AliMUONChamber.h"
#include "AliHitMap.h"
#include "AliMUONHitMapA1.h"
#include "AliMUON.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONDigit.h"
#include "AliMUONTransientDigit.h"
#include "AliRun.h"
#include "AliPDG.h"

ClassImp(AliMUONMerger)

//___________________________________________
AliMUONMerger::AliMUONMerger()
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
    fPadHitsBgr = 0;
    fHitMap     = 0;
    fList       = 0;
    fTrList     = 0;
    fAddress    = 0; 
}

//------------------------------------------------------------------------
AliMUONMerger::~AliMUONMerger()
{
// Destructor
    if (fTrH1)       delete fTrH1;
    if (fHitsBgr)    delete fHitsBgr;
    if (fPadHitsBgr) delete fPadHitsBgr;
    if (fHitMap)     delete fHitMap;
    if (fList)       delete fList;
    if (fTrList)     delete fTrList;
    if (fAddress)    delete fAddress; 
}

//------------------------------------------------------------------------
Bool_t AliMUONMerger::Exists(const AliMUONPadHit *mergable)
{
    AliMUONPadHit *padhit = (AliMUONPadHit*) mergable;
    return (fHitMap[fNch]->TestHit(padhit->PadX(),padhit->PadY()));
}

//------------------------------------------------------------------------
void AliMUONMerger::Update(AliMUONPadHit *mergable)
{
    AliMUONPadHit *padhit = (AliMUONPadHit*) mergable;    
    AliMUONTransientDigit* pdigit;
    Int_t ipx      = padhit->PadX();        // pad number on X
    Int_t ipy      = padhit->PadY();        // pad number on Y
    Int_t iqpad    = Int_t(padhit->QPad()); // charge per pad

    pdigit = (AliMUONTransientDigit*) fHitMap[fNch]->GetHit(ipx, ipy);
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
void AliMUONMerger::CreateNew(AliMUONPadHit *mergable)
{
    AliMUONPadHit *padhit = (AliMUONPadHit*) mergable;    
    AliMUONTransientDigit* pdigit;

    Int_t ipx      = padhit->PadX();       // pad number on X
    Int_t ipy      = padhit->PadY();       // pad number on Y
    fList->AddAtAndExpand(
	new AliMUONTransientDigit(fNch,fDigits),fCounter);
    fHitMap[fNch]->SetHit(ipx, ipy, fCounter);
    fCounter++;
    pdigit = (AliMUONTransientDigit*)fList->At(fList->GetLast());
    // list of tracks
    TObjArray *fTrList = (TObjArray*)pdigit->TrackList();
    TVector &trinfo    = *((TVector*) (*fAddress)[fCountadr-1]);
    fTrList->Add(&trinfo);
}


//------------------------------------------------------------------------
void AliMUONMerger::Init()
{
// Initialisation
    
    if (fMerge) fBgrFile = InitBgr();
}



//------------------------------------------------------------------------
TFile* AliMUONMerger::InitBgr()
{
// Initialise background event
    TFile *file = new TFile(fFnBgr);
// add error checking later
    printf("\n AliMUONMerger has opened %s file with background event \n", fFnBgr);
    fHitsBgr     = new TClonesArray("AliMUONHit",1000);
    fPadHitsBgr  = new TClonesArray("AliMUONPadHit",1000);
    return file;
}

//------------------------------------------------------------------------
void AliMUONMerger::Digitise()
{

    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  
    AliMUONChamber*   iChamber;
    AliSegmentation*  segmentation;

    fList = new TObjArray;
    if(!fAddress) fAddress = new TClonesArray("TVector",10000);

    AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
    fHitMap= new AliHitMap* [AliMUONConstants::NCh()];
    for (Int_t i = 0; i < AliMUONConstants::NCh(); i++) {fHitMap[i] = 0;}
    if (fMerge ) {
	fBgrFile->cd();
	fBgrFile->ls();
        //
	// Get Hits Tree header from file
	if(fHitsBgr) fHitsBgr->Clear();
	if(fPadHitsBgr) fPadHitsBgr->Clear();
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
	sprintf(branchname,"%s",pMUON->GetName());
	if (fTrH1 && fHitsBgr) {
	    branch = fTrH1->GetBranch(branchname);
	    if (branch) branch->SetAddress(&fHitsBgr);
	}
	if (fTrH1 && fPadHitsBgr) {
	    branch = fTrH1->GetBranch("MUONCluster");
	    if (branch) branch->SetAddress(&fPadHitsBgr);
	}
    }
    //
    // loop over cathodes
    //
    AliHitMap* hm;
    fCountadr = 0;
    for (int icat = 0; icat < 2; icat++) { 
	fCounter = 0;
	Int_t * nmuon = new Int_t [AliMUONConstants::NCh()];
	for (Int_t i = 0; i < AliMUONConstants::NCh(); i++) {
	    iChamber = &(pMUON->Chamber(i));
	    if (iChamber->Nsec() == 1 && icat == 1) {
		continue;
	    } else {
		segmentation = iChamber->SegmentationModel(icat+1);
	    }
	    fHitMap[i] = new AliMUONHitMapA1(segmentation, fList);
	    nmuon[i] = 0;
	}

//
//   Loop over tracks
//

	TTree *treeH  = gAlice->TreeH();
	Int_t ntracks = (Int_t) treeH->GetEntries();
	Int_t jj;

	Float_t ** xhit = new Float_t * [AliMUONConstants::NCh()];
	for (jj = 0; jj < AliMUONConstants::NCh(); jj++) xhit[jj] = new Float_t[2];
	Float_t ** yhit = new Float_t * [AliMUONConstants::NCh()];
	for (jj = 0; jj < AliMUONConstants::NCh(); jj++) yhit[jj] = new Float_t[2];

	for (fTrack = 0; fTrack < ntracks; fTrack++) {
	    gAlice->ResetHits();
	    treeH->GetEvent(fTrack);
//
//   Loop over hits
	    for(AliMUONHit* mHit = (AliMUONHit*)pMUON->FirstHit(-1); 
		mHit;
		mHit = (AliMUONHit*)pMUON->NextHit()) 
	    {
		fNch = mHit->Chamber()-1;  // chamber number
		if (fNch > AliMUONConstants::NCh()-1) continue;
		iChamber = &(pMUON->Chamber(fNch));
		/*
		if (fMerge) {
		    if (mHit->Particle() == kMuonPlus || 
			mHit->Particle() == kMuonMinus) {
			xhit[fNch][nmuon[fNch]] = mHit->X();
			yhit[fNch][nmuon[fNch]] = mHit->Y();
			nmuon[fNch]++;
			if (nmuon[fNch] > 2) printf("MUON: nmuon %d\n",nmuon[fNch]);
		    }
		}
		*/
		
//
// Loop over pad hits
		for (AliMUONPadHit* mPad =
			 (AliMUONPadHit*)pMUON->FirstPad(mHit,pMUON->PadHits());
		     mPad;
		     mPad = (AliMUONPadHit*)pMUON->NextPad(pMUON->PadHits()))
		{
		    Int_t cathode  = mPad->Cathode();      // cathode number
		    Int_t ipx      = mPad->PadX();         // pad number on X
		    Int_t ipy      = mPad->PadY();         // pad number on Y
		    Int_t iqpad    = Int_t(mPad->QPad());  // charge per pad
		    if (cathode != (icat+1)) continue;

		    segmentation = iChamber->SegmentationModel(cathode);

		    new((*fAddress)[fCountadr++]) TVector(2);

		    TVector &trinfo = *((TVector*) (*fAddress)[fCountadr-1]);
		    trinfo(0) = (Float_t)fTrack;
		    trinfo(1) = (Float_t)iqpad;

		    fDigits[0] = ipx;
		    fDigits[1] = ipy;
		    fDigits[2] = icat;
		    fDigits[3] = iqpad;
		    fDigits[4] = iqpad;
		    if (mHit->Particle() == kMuonPlus ||
			mHit->Particle() == kMuonMinus) {
			fDigits[5] = mPad->HitNumber();
		    } else fDigits[5] = -1;

		    // build the list of fired pads and update the info

		    if (!Exists(mPad)) {
			CreateNew(mPad);
		    } else {
			Update(mPad);
		    } //  end if pdigit
		} //end loop over clusters
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
		if (fPadHitsBgr)    fPadHitsBgr->Clear();

		fTrH1->GetEvent(fTrack);
//
//   Loop over hits
		AliMUONHit* mHit;
		for(int i = 0; i < fHitsBgr->GetEntriesFast(); ++i) 
		{	
		    mHit   = (AliMUONHit*) (*fHitsBgr)[i];
		    fNch   = mHit->Chamber()-1;  // chamber number
		    iChamber = &(pMUON->Chamber(fNch));
                    Float_t xbgr = mHit->X();
		    Float_t ybgr = mHit->Y();
		    Bool_t cond  = kFALSE;
		    /*
		    for (Int_t imuon = 0; imuon < nmuon[fNch]; imuon++) {
			Float_t dist = (xbgr-xhit[fNch][imuon])*(xbgr-xhit[fNch][imuon])
			    +(ybgr-yhit[fNch][imuon])*(ybgr-yhit[fNch][imuon]);
			if (dist < 100.) cond = kTRUE;
		    }
		    */
		    cond  = kTRUE;
//
// Loop over pad hits
		    for (AliMUONPadHit* mPad =
			     (AliMUONPadHit*)pMUON->FirstPad(mHit,fPadHitsBgr);
			 mPad;
			 mPad = (AliMUONPadHit*)pMUON->NextPad(fPadHitsBgr))
		    {
			Int_t cathode  = mPad->Cathode();     // cathode number
			Int_t ipx      = mPad->PadX();        // pad number on X
			Int_t ipy      = mPad->PadY();        // pad number on Y
			Int_t iqpad    = Int_t(mPad->QPad()); // charge per pad

			if (cathode != (icat+1)) continue;
			new((*fAddress)[fCountadr++]) TVector(2);
			TVector &trinfo = *((TVector*) (*fAddress)[fCountadr-1]);
			trinfo(0) = kBgTag;  // tag background
			trinfo(1) = kBgTag;
			
			fDigits[0] = ipx;
			fDigits[1] = ipy;
			fDigits[2] = icat;
			fDigits[3] = iqpad;
			fDigits[4] = 0;
			fDigits[5] = -1;
			
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
	delete [] xhit;
	delete [] yhit;

	Int_t tracks[10];
	Int_t charges[10];
	Int_t nentries = fList->GetEntriesFast();
	
	for (Int_t nent = 0; nent < nentries; nent++) {
	    AliMUONTransientDigit *address = (AliMUONTransientDigit*)fList->At(nent);
	    if (address == 0) continue; 
	    Int_t ich = address->Chamber();
	    Int_t   q = address->Signal(); 
	    iChamber = &(pMUON->Chamber(ich));
//
//  Digit Response (noise, threshold, saturation, ...)
	    AliMUONResponse * response = iChamber->ResponseModel();
	    q = response->DigitResponse(q);
	    
	    if (!q) continue;
	    
	    fDigits[0] = address->PadX();
	    fDigits[1] = address->PadY();
	    fDigits[2] = address->Cathode();
	    fDigits[3] = q;
	    fDigits[4] = address->Physics();
	    fDigits[5] = address->Hit();
	    
	    TObjArray* fTrList = (TObjArray*)address->TrackList();
	    Int_t nptracks = fTrList->GetEntriesFast();

	    // this was changed to accomodate the real number of tracks

	    if (nptracks > 10) {
		printf("\n Attention - nptracks > 10 %d \n", nptracks);
		nptracks = 10;
	    }
	    if (nptracks > 2) {
		printf("Attention - nptracks > 2  %d \n",nptracks);
		printf("cat,ich,ix,iy,q %d %d %d %d %d \n",icat,ich,fDigits[0],fDigits[1],q);
	    }
	    for (Int_t tr = 0; tr < nptracks; tr++) {
		TVector *ppP = (TVector*)fTrList->At(tr);
		if(!ppP ) printf("ppP - %p\n",ppP);
		TVector &pp  = *ppP;
		tracks[tr]   = Int_t(pp(0));
		charges[tr]  = Int_t(pp(1));
	    }      //end loop over list of tracks for one pad
            // Sort list of tracks according to charge
	    if (nptracks > 1) {
		pMUON->SortTracks(tracks,charges,nptracks);
	    }
	    if (nptracks < 10 ) {
		for (Int_t i = nptracks; i < 10; i++) {
		    tracks[i]  = 0;
		    charges[i] = 0;
		}
	    }
	    
	    // fill digits
	    pMUON->AddDigits(ich,tracks,charges,fDigits);
	    // delete fTrList;
	}
	gAlice->TreeD()->Fill();
	pMUON->ResetDigits();
	fList->Delete();

	
	for(Int_t ii = 0; ii < AliMUONConstants::NCh(); ++ii) {
	    if (fHitMap[ii]) {
		hm=fHitMap[ii];
		delete hm;
		fHitMap[ii] = 0;
	    }
	}
	delete [] nmuon;    
    } //end loop over cathodes
    delete [] fHitMap;
    char hname[30];
    sprintf(hname,"TreeD%d",fEvNrSig);
    gAlice->TreeD()->Write(hname);
    // reset tree
    gAlice->TreeD()->Reset();
    delete fList;
    
    if (fAddress)    fAddress->Delete();
    if (fHitsBgr)    fHitsBgr->Delete();
    if (fPadHitsBgr) fPadHitsBgr->Delete();
    // gObjectTable->Print();
}





