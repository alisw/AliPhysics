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
#include <TPDGCode.h>

#include "AliHitMap.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONHit.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONMerger.h"
#include "AliMUONPadHit.h"
#include "AliMUONTransientDigit.h"
#include "AliRun.h"

ClassImp(AliMUONMerger)

//----------------------------------------------------------------------
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
    fBgrFile    = 0;
    fDebug      = 0;
}

//----------------------------------------------------------------------
AliMUONMerger::AliMUONMerger(const AliMUONMerger&)
{
// Protected copy constructor

  Fatal("AliMUONMergerModule", "Not implemented.");
}

//------------------------------------------------------------------------
AliMUONMerger::~AliMUONMerger()
{
// Destructor
    if (fTrH1)       delete fTrH1;
    if (fHitsBgr)    delete fHitsBgr;
    if (fPadHitsBgr) delete fPadHitsBgr;
    if (fHitMap)     delete [] fHitMap;
    if (fList)       delete fList;
    if (fBgrFile)    delete fBgrFile;
}

//----------------------------------------------------------------------
AliMUONMerger&  AliMUONMerger::operator=(const AliMUONMerger& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  Fatal("operator=", "Not implemented.");
    
  return *this;  
}    
          
//------------------------------------------------------------------------
Bool_t AliMUONMerger::Exists(const AliMUONPadHit *padhit) const
{
// test if the given padhit was already fired
    return (fHitMap[fNch]->TestHit(padhit->PadX(),padhit->PadY()));
}

//------------------------------------------------------------------------
void AliMUONMerger::Update(AliMUONPadHit *padhit)
{
// add new contribution to the fired padhit
    AliMUONTransientDigit *pdigit = 
      static_cast<AliMUONTransientDigit*>(
      fHitMap[fNch]->GetHit(padhit->PadX(),padhit->PadY()));

    // update charge
    //
    Int_t iqpad    = padhit->QPad();        // charge per pad
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
void AliMUONMerger::CreateNew(AliMUONPadHit *padhit)
{
// add new transient digit to the list, update hit map
    fList->AddAtAndExpand(
	new AliMUONTransientDigit(fNch,fDigits),fCounter);
    fHitMap[fNch]->SetHit(padhit->PadX(),padhit->PadY(),fCounter);
    AliMUONTransientDigit* pdigit = 
      static_cast<AliMUONTransientDigit*>(fList->Last());
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
void AliMUONMerger::Init()
{
// Initialisation
    // open only once the background file !!
    if (fMerge && !fBgrFile) fBgrFile = InitBgr();
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
  // Obsolete sep 2003 Gines MARTINEZ 
  //    AliMUONChamber*   iChamber;
//      AliSegmentation*  segmentation;

//      fList = new TObjArray;

//      AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
//      fHitMap= new AliHitMap* [AliMUONConstants::NCh()];
//      if (fMerge ) {
//  	fBgrFile->cd();
//          //
//  	// Get Hits Tree header from file
//  	//if(fHitsBgr) fHitsBgr->Clear();     // Useless because line 327
//  	//if(fPadHitsBgr) fPadHitsBgr->Clear(); // Useless because line 328
//  	if(fTrH1) delete fTrH1;
//  	fTrH1 = 0;
	
//  	char treeName[20];
//  	sprintf(treeName,"TreeH%d",fEvNrBgr);
//  	fTrH1 = (TTree*)gDirectory->Get(treeName);
//  	if (!fTrH1) {
//  	    printf("ERROR: cannot find Hits Tree for event:%d\n",fEvNrBgr);
//  	}
//          //
//  	// Set branch addresses
//  	TBranch *branch;
//  	char branchname[20];
//  	sprintf(branchname,"%s",pMUON->GetName());
//  	if (fTrH1 && fHitsBgr) {
//  	    branch = fTrH1->GetBranch(branchname);
//  	    if (branch) branch->SetAddress(&fHitsBgr);
//  	}
//  	if (fTrH1 && fPadHitsBgr) {
//  	    branch = fTrH1->GetBranch("MUONCluster");
//  	    if (branch) branch->SetAddress(&fPadHitsBgr);
//  	}
//      }
//      //
//      // loop over cathodes
//      //
//      fSignal = kTRUE;
//      for (int icat = 0; icat < 2; icat++) { 
//  	fCounter = 0;
//  	for (Int_t i = 0; i < AliMUONConstants::NCh(); i++) {
//  	    iChamber = &(pMUON->Chamber(i));
//  	    if (iChamber->Nsec() == 1 && icat == 1) {
//  		continue;
//  	    } else {
//  		segmentation = iChamber->SegmentationModel(icat+1);
//  	    }
//  	    fHitMap[i] = new AliMUONHitMapA1(segmentation, fList);
//  	}

//  //
//  //   Loop over tracks
//  //

//  /******************************************************************/
//        TTree* treeH = pMUON->TreeH();
//        if (treeH == 0x0)
//         {
//           cerr<<"AliMUONMerger::Exec: Can not get TreeH"<<endl;
//           return;
//         }
//  /******************************************************************/     

	
	
//  	Int_t ntracks = (Int_t) treeH->GetEntries();
//  	treeH->SetBranchStatus("*",0); // switch off all branches
//          treeH->SetBranchStatus("MUON*",1); // switch on only MUON

//  	for (fTrack = 0; fTrack < ntracks; fTrack++) {
//  	    gAlice->ResetHits();
//  	    treeH->GetEntry(fTrack,0);
//  //
//  //   Loop over hits
//  	    for(AliMUONHit* mHit = (AliMUONHit*)pMUON->FirstHit(-1); 
//  		mHit;
//  		mHit = (AliMUONHit*)pMUON->NextHit()) 
//  	    {
//  		fNch = mHit->Chamber()-1;  // chamber number
//  		if (fNch > AliMUONConstants::NCh()-1) continue;
//  		iChamber = &(pMUON->Chamber(fNch));
		
//  //
//  // Loop over pad hits
//  		for (AliMUONPadHit* mPad =
//  			 (AliMUONPadHit*)pMUON->FirstPad(mHit,pMUON->PadHits());
//  		     mPad;
//  		     mPad = (AliMUONPadHit*)pMUON->NextPad(pMUON->PadHits()))
//  		{
//  		    Int_t cathode  = mPad->Cathode();      // cathode number
//  		    if (cathode != (icat+1)) continue;
//  		    Int_t iqpad    = Int_t(mPad->QPad());  // charge per pad
//  //		    segmentation = iChamber->SegmentationModel(cathode);
//  		    fDigits[0] = mPad->PadX();  
//  		    fDigits[1] = mPad->PadY();
//  		    if (!(fHitMap[fNch]->ValidateHit(fDigits[0], fDigits[1]))) continue;
//  		    fDigits[2] = icat;
//  		    fDigits[3] = iqpad;
//  		    fDigits[4] = iqpad;
//  		    if (mHit->Particle() == kMuonPlus ||
//  			mHit->Particle() == kMuonMinus) {
//  			fDigits[5] = mPad->HitNumber();
//  		    } else fDigits[5] = -1;

//  		    // build the list of fired pads and update the info

//  		    if (!Exists(mPad)) {
//  			CreateNew(mPad);
//  		    } else {
//  			Update(mPad);
//  		    } //  end if pdigit
//  		} //end loop over clusters
//  	    } // hit loop
//  	} // track loop

//  	// open the file with background
       
//  	if (fMerge) {
//              fSignal = kFALSE;
//  	    ntracks = (Int_t)fTrH1->GetEntries();
//  //
//  //   Loop over tracks
//  //
//  	    for (fTrack = 0; fTrack < ntracks; fTrack++) {

//  		if (fHitsBgr)       fHitsBgr->Clear();
//  		if (fPadHitsBgr)    fPadHitsBgr->Clear();

//  		fTrH1->GetEvent(fTrack);
//  //
//  //   Loop over hits
//  		AliMUONHit* mHit;
//  		for(Int_t i = 0; i < fHitsBgr->GetEntriesFast(); ++i) 
//  		{	
//  		    mHit   = (AliMUONHit*) (*fHitsBgr)[i];
//  		    fNch   = mHit->Chamber()-1;  // chamber number
//  		    iChamber = &(pMUON->Chamber(fNch));
//  //
//  // Loop over pad hits
//  		    for (AliMUONPadHit* mPad =
//  			     (AliMUONPadHit*)pMUON->FirstPad(mHit,fPadHitsBgr);
//  			 mPad;
//  			 mPad = (AliMUONPadHit*)pMUON->NextPad(fPadHitsBgr))
//  		    {
//  			Int_t cathode  = mPad->Cathode();     // cathode number
//  			Int_t ipx      = mPad->PadX();        // pad number on X
//  			Int_t ipy      = mPad->PadY();        // pad number on Y
//  			Int_t iqpad    = Int_t(mPad->QPad()); // charge per pad
//  			if (!(fHitMap[fNch]->ValidateHit(ipx, ipy))) continue;

//  			if (cathode != (icat+1)) continue;
//  			fDigits[0] = ipx;
//  			fDigits[1] = ipy;
//  			fDigits[2] = icat;
//  			fDigits[3] = iqpad;
//  			fDigits[4] = 0;
//  			fDigits[5] = -1;
			
//  			// build the list of fired pads and update the info
//  			if (!Exists(mPad)) {
//  			    CreateNew(mPad);
//  			} else {
//  			    Update(mPad);
//  			} //  end if !Exists
//  		    } //end loop over clusters
//  		} // hit loop
//  	    } // track loop

//  	    TTree *treeK = gAlice->TreeK();
//              TFile *file = NULL;
	    
//  	    if (treeK) file = treeK->GetCurrentFile();
//  	    file->cd();
//  	} // if fMerge

//  	Int_t tracks[kMAXTRACKS];
//  	Int_t charges[kMAXTRACKS];
//  	Int_t nentries = fList->GetEntriesFast();
	
//  	for (Int_t nent = 0; nent < nentries; nent++) {
//  	    AliMUONTransientDigit *address = (AliMUONTransientDigit*)fList->At(nent);
//  	    if (address == 0) continue; 
//  	    Int_t ich = address->Chamber();
//  	    Int_t   q = address->Signal(); 
//  	    iChamber = &(pMUON->Chamber(ich));
//  //
//  //  Digit Response (noise, threshold, saturation, ...)
//  	    AliMUONResponse * response = iChamber->ResponseModel();
//  	    q = response->DigitResponse(q,address);
	    
//  	    if (!q) continue;
	    
//  	    fDigits[0] = address->PadX();
//  	    fDigits[1] = address->PadY();
//  	    fDigits[2] = address->Cathode();
//  	    fDigits[3] = q;
//  	    fDigits[4] = address->Physics();
//  	    fDigits[5] = address->Hit();
	    
//  	    Int_t nptracks = address->GetNTracks();

//  	    if (nptracks > kMAXTRACKS) {
//  	        if (fDebug>0)
//  		  printf("\n Attention - nptracks > kMAXTRACKS %d \n", nptracks);
//  		nptracks = kMAXTRACKS;
//  	    }
//  	    if (nptracks > 2) {
//  	        if (fDebug>0) {
//  		  printf("Attention - nptracks > 2  %d \n",nptracks);
//  		  printf("cat,ich,ix,iy,q %d %d %d %d %d \n",icat,ich,fDigits[0],fDigits[1],q);
//  		}
//  	    }
//  	    for (Int_t tr = 0; tr < nptracks; tr++) {
//  		tracks[tr]   = address->GetTrack(tr);
//  		charges[tr]  = address->GetCharge(tr);
//  	    }      //end loop over list of tracks for one pad
//              // Sort list of tracks according to charge
//  	    if (nptracks > 1) {
//  		SortTracks(tracks,charges,nptracks);
//  	    }
//  	    if (nptracks < kMAXTRACKS ) {
//  		for (Int_t i = nptracks; i < kMAXTRACKS; i++) {
//  		    tracks[i]  = 0;
//  		    charges[i] = 0;
//  		}
//  	    }
	    
//  	    // fill digits
//  	    pMUON->AddDigits(ich,tracks,charges,fDigits);
//  	}
//  	gAlice->TreeD()->Fill();
//  	pMUON->ResetDigits();
//  	fList->Delete();

	
//  	for(Int_t ii = 0; ii < AliMUONConstants::NCh(); ++ii) {
//  	    if (fHitMap[ii]) {
//  		delete fHitMap[ii];
//  		fHitMap[ii] = 0;
//  	    }
//  	}
//      } //end loop over cathodes
//      delete [] fHitMap;
//      delete fList;
    
//  no need to delete ... and it makes a crash also
//    if (fHitsBgr)    fHitsBgr->Delete();
//    if (fPadHitsBgr) fPadHitsBgr->Delete();
    // gObjectTable->Print();
}



void AliMUONMerger::SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr)
{
  //
  // Sort the list of tracks contributing to a given digit
  // Only the 3 most significant tracks are acctually sorted
  //
  
  //
  //  Loop over signals, only 3 times
  //
  
  Int_t qmax;
  Int_t jmax;
  Int_t idx[3] = {-2,-2,-2};
  Int_t jch[3] = {-2,-2,-2};
  Int_t jtr[3] = {-2,-2,-2};
  Int_t i,j,imax;
  
  if (ntr<3) imax=ntr;
  else imax=3;
  for(i=0;i<imax;i++){
    qmax=0;
    jmax=0;
    
    for(j=0;j<ntr;j++){
      
      if((i == 1 && j == idx[i-1]) 
	 ||(i == 2 && (j == idx[i-1] || j == idx[i-2]))) continue;
      
      if(charges[j] > qmax) {
	qmax = charges[j];
	jmax=j;
      }       
    } 
    
    if(qmax > 0) {
      idx[i]=jmax;
      jch[i]=charges[jmax]; 
      jtr[i]=tracks[jmax]; 
    }
    
  } 
  
  for(i=0;i<3;i++){
    if (jtr[i] == -2) {
         charges[i]=0;
         tracks[i]=0;
    } else {
         charges[i]=jch[i];
         tracks[i]=jtr[i];
    }
  }
}



