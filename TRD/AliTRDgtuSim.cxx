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

/* $Id: AliTRDgtuSim.cxx 28397 2008-09-02 09:33:00Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  GTU simulation                                                        //
//                                                                        //
//  Authors: J. Klein (Jochen.Klein@cern.ch)                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"

#include "AliTRDgtuSim.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDgtuTMU.h"
#include "AliTRDtrackGTU.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDtrackletMCM.h"
#include "AliESDEvent.h"

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"

ClassImp(AliTRDgtuSim)

AliTRDgtuSim::AliTRDgtuSim(AliRunLoader *rl) 
  : TObject(),
  fRunLoader(rl),
  fTMU(0x0),
  fTrackletArray(0x0),
  fTrackTree(0x0),
  fTrackletTree(0x0)
{
  fTrackletTree = new TTree("gtutracklets", "Tree with GTU tracklets");
  fTrackletTree->SetDirectory(0);
}

AliTRDgtuSim::~AliTRDgtuSim() 
{
  if (fTrackletArray)
    fTrackletArray->Delete();
  delete fTrackletArray;
  delete fTrackletTree;
}

Bool_t AliTRDgtuSim::RunGTUFromTrackletFile(TString filename, Int_t event, Int_t noev) 
{
    AliInfo(Form("Running the GTU simulation on file: %s", filename.Data()));
    ifstream input(filename.Data());
    
    std::string str;
    TString string;
    int lineno = 0;
    
    Int_t iEventPrev = -1;
    Int_t iStackPrev = -1;
    Int_t iSecPrev = -1;
    Int_t iSec = -1;
    Int_t iStack = -1;
    Int_t iLink = -1;
    Int_t iEvent = -1;
    Int_t evcnt = -1;
    
    fTMU = 0x0;
    
    AliInfo("--------- Reading from file ----------");
    while (getline(input, str)) {
	lineno++;
	string = str;
	AliInfo(Form("Line %i : %s", lineno, string.Data()));
	
	TObjArray *tokens = string.Tokenize(" ");
	if (tokens->GetEntriesFast() < 7) {
	    AliWarning(Form("Invalid input in line %i, too few parameters", lineno));
	    continue;
	}

	if ( ((TObjString*) tokens->At(0))->GetString().Atoi() < event) 
	    continue;

	iEvent = ((TObjString*) tokens->At(0))->GetString().Atoi();
	iSec = ((TObjString*) tokens->At(1))->GetString().Atoi();
	iStack = ((TObjString*) tokens->At(2))->GetString().Atoi();
	iLink = 2 * ((TObjString*) tokens->At(3))->GetString().Atoi() + ((TObjString*) tokens->At(4))->GetString().Atoi();

	if (iEvent != iEventPrev || iStack != iStackPrev || iSec != iSecPrev) {
	    if(fTMU) {
		TList *ListOfTracks = new TList();
		fTMU->SetStack(iStackPrev);
		fTMU->SetSector(iSecPrev);
		fTMU->RunTMU(ListOfTracks);
		AliInfo(Form("--- There are %i tracks. Writing ...", ListOfTracks->GetEntries()));
		WriteTracksToTree(ListOfTracks);
		fTMU->WriteTrackletsToTree(fTrackletTree);
		WriteTracksToDataFile(ListOfTracks, iEventPrev);
		if (ListOfTracks->GetEntries() > 0) 
		    AliInfo(Form("   %d GeV/c", ((AliTRDtrackGTU*) ListOfTracks->At(0))->GetPt() ));
		delete fTMU;
		fTMU = new AliTRDgtuTMU(); 
		delete ListOfTracks;
		ListOfTracks = 0x0;
	    } else {
		fTMU = new AliTRDgtuTMU();
	    }
	    iStackPrev = iStack;
	    iSecPrev = iSec;
	    iEventPrev = iEvent;
	    evcnt++;
	    if (evcnt == noev)
		break;
	}
	for (Int_t i = 5; i < tokens->GetEntriesFast(); i++) {
	    UInt_t trackletWord = 0;
	    sscanf(((TObjString*) tokens->At(i))->GetString().Data(), "%x", &trackletWord);
	    if (trackletWord == 0x10001000) 
		break;
	    AliInfo(Form("%i. tracklet: %s -> 0x%08x", i-4, ((TObjString*) tokens->At(i))->GetString().Data(), trackletWord));
	    AliTRDtrackletWord *trkl = new AliTRDtrackletWord(trackletWord);
	    fTMU->AddTracklet(trkl, iLink);
	}
    }
    
    if (fTMU && evcnt < noev) {
	TList *ListOfTracks = new TList();
	fTMU->SetStack(iStackPrev);
	fTMU->SetSector(iSecPrev);
	fTMU->RunTMU(ListOfTracks);
	WriteTracksToTree(ListOfTracks);
	fTMU->WriteTrackletsToTree(fTrackletTree);
	WriteTracksToDataFile(ListOfTracks, iEventPrev);
	delete fTMU;
	delete ListOfTracks;
	fTMU = 0x0;
    }

    AliInfo(Form("Analyzed %i events", evcnt));
    return kTRUE; 
}

Bool_t AliTRDgtuSim::RunGTU(AliLoader *loader, AliESDEvent *esd) 
{
    if (!LoadTracklets(loader)) {
	AliError("Could not load the tracklets. Aborting ...");
	return kFALSE;
    }

    Int_t iStackPrev = -1;
    Int_t iSecPrev = -1;
    Int_t iSec = -1;
    Int_t iStack = -1;
    Int_t iLink = -1;

    if (fTMU) {
	delete fTMU;
	fTMU = 0x0;
    }
    
    TList *ListOfTracks = new TList();
    
    TIter next(fTrackletArray);
    AliTRDtrackletBase *trkl;

    while ((trkl = (AliTRDtrackletBase*) next())) {
	iSec = trkl->GetDetector() / 30;
	iStack = (trkl->GetDetector() % 30) / 6;
	iLink = 2 * (trkl->GetDetector() % 6) + (trkl->GetYbin() < 0 ? 0 : 1);

	if (iStack != iStackPrev || iSec != iSecPrev) {
	    if(fTMU) {
		fTMU->SetStack(iStackPrev);
		fTMU->SetSector(iSecPrev);
		fTMU->RunTMU(ListOfTracks);
		WriteTracksToTree(ListOfTracks);
		fTMU->WriteTrackletsToTree(fTrackletTree);
		WriteTracksToESD(ListOfTracks, esd);
		fTMU->Reset();
		ListOfTracks->Delete();
	    } else {
		fTMU = new AliTRDgtuTMU();
	    }
	    iStackPrev = iStack;
	    iSecPrev = iSec;
	}
	fTMU->AddTracklet(trkl, iLink);
    }
    
    if (fTMU) {
	fTMU->SetStack(iStackPrev);
	fTMU->SetSector(iSecPrev);
	fTMU->RunTMU(ListOfTracks);
	WriteTracksToTree(ListOfTracks);
	fTMU->WriteTrackletsToTree(fTrackletTree);
	WriteTracksToESD(ListOfTracks, esd);
	delete fTMU;
	fTMU = 0x0;
	ListOfTracks->Delete();
    }

    delete ListOfTracks;

    return kTRUE;
}

Bool_t AliTRDgtuSim::LoadTracklets(AliLoader *loader) 
{
  AliInfo("Loading tracklets ...");

  if (!loader) {
    AliError("No loader given!");
    return kFALSE;
  }

  AliDataLoader *trackletLoader = loader->GetDataLoader("tracklets");
  if (!trackletLoader) {
      AliError("No tracklet loader found!");
      return kFALSE;
  }

  trackletLoader->Load();
  
  TTree *trackletTree = trackletLoader->Tree();
  if (!trackletTree) {
    AliError("No tracklet tree found");
    return kFALSE;
  }  


  TBranch *trklbranch = trackletTree->GetBranch("mcmtrklbranch");
  if (trklbranch) {
      if (!fTrackletArray)
	  fTrackletArray = new TClonesArray("AliTRDtrackletMCM", 1000);
      else if ((TClass::GetClass("AliTRDtrackletMCM"))->InheritsFrom(fTrackletArray->Class()))
	  fTrackletArray->Delete();
      else {
	  fTrackletArray->Delete();
	  delete fTrackletArray;
	  fTrackletArray = new TClonesArray("AliTRDtrackletMCM", 1000);
      }

      AliTRDtrackletMCM *trkl = new AliTRDtrackletMCM; 
      trklbranch->SetAddress(&trkl);
      for (Int_t iTracklet = 0; iTracklet < trklbranch->GetEntries(); iTracklet++) {
	  trklbranch->GetEntry(iTracklet);
	  new ((*fTrackletArray)[fTrackletArray->GetEntries()]) AliTRDtrackletMCM(*trkl);
      }
      return kTRUE;
  }

  trklbranch = trackletTree->GetBranch("trkbranch");

  if (!trklbranch) {
      AliError("Could not get trkbranch");
      return kFALSE;
  }
  
  if (!fTrackletArray)
      fTrackletArray = new TClonesArray("AliTRDtrackletWord", 1000);
  else if ((TClass::GetClass("AliTRDtrackletWord"))->InheritsFrom(fTrackletArray->Class()))
      fTrackletArray->Delete();
  else {
      fTrackletArray->Delete();
      delete fTrackletArray;
      fTrackletArray = new TClonesArray("AliTRDtrackletWord", 1000);
  }

  Int_t notrkl = 0;
  UInt_t *leaves = new UInt_t[258];
  AliInfo(Form("No. of entries: %i", trklbranch->GetEntries()));
  
  for (Int_t iEntry = 0; iEntry < trklbranch->GetEntries(); iEntry++) {
      trklbranch->SetAddress(leaves);
      trklbranch->GetEntry(iEntry);
      for (Int_t iTracklet = 0; iTracklet < 256; iTracklet++) {
	if (leaves[2 + iTracklet] == 0)
	  break;
	new((*fTrackletArray)[notrkl]) AliTRDtrackletWord(leaves[2 + iTracklet], leaves[0] + leaves[1]);
	notrkl++;
      }
      AliInfo(Form("Entry: %3i: Det: %3i, side: %i, 1st tracklet: 0x%08x, no: %i", iEntry, leaves[0], leaves[1], leaves[2], notrkl));
  }

  return kTRUE;
}

Bool_t AliTRDgtuSim::WriteTracksToDataFile(TList *ListOfTracks, Int_t event) 
{
    Int_t sm = 0;
    Int_t stack = 0;

    FILE *out;
    out = fopen("test.data", "a");

    AliInfo(Form("%i tracks found in event %i", ListOfTracks->GetSize(), event));
    fprintf(out, "0 %5i %2i %i  00000000\n", event, sm, stack);
    for (Int_t i = 0; i < ListOfTracks->GetSize(); i++) {
	AliTRDtrackGTU *trk = (AliTRDtrackGTU*) ListOfTracks->At(i);
	sm = trk->GetSector();
	stack = trk->GetStack();
	fprintf(out, "1 %5i %2i %2i %3i %3i %3i %3i %3i %3i %3i %4i %f\n", event, sm, stack, trk->GetTrackletMask(), 
	       trk->GetTrackletIndex(5), 
	       trk->GetTrackletIndex(4), 
	       trk->GetTrackletIndex(3), 
	       trk->GetTrackletIndex(2), 
	       trk->GetTrackletIndex(1), 
	       trk->GetTrackletIndex(0),
	       trk->GetPtInt(), 
	       trk->GetPt());
    }
    fclose(out);
    return kTRUE;
}

Bool_t AliTRDgtuSim::WriteTracksToTree(TList *ListOfTracks, Int_t /*event*/) 
{
  AliInfo(Form("Writing %i tracks to the tree...", ListOfTracks->GetEntries()));

  if (!ListOfTracks)
    return kFALSE;

  if (ListOfTracks->GetEntries() <= 0) 
    return kTRUE;

  if (!fTrackTree) {
    fTrackTree = new TTree("gtutracks", "GTU tracks");
    fTrackTree->SetDirectory(0);
  }

  AliTRDtrackGTU *trk = 0x0;
  TBranch *branch = fTrackTree->GetBranch("TRDgtuTrack");
  if (!branch) {
      branch = fTrackTree->Branch("TRDgtuTrack", "AliTRDtrackGTU", &trk, 32000, 99);
  }

  TIter next(ListOfTracks);
  while ((trk = (AliTRDtrackGTU*) next())) {
      trk->CookLabel();
      branch->SetAddress(&trk);
      fTrackTree->Fill();   
  }
  fTrackTree->ResetBranchAddress(branch);

  return kTRUE; 
}

Bool_t AliTRDgtuSim::WriteTreesToFile() {
  TFile *f = TFile::Open("TRD.GtuTracking.root", "RECREATE");
  f->cd();
  if (fTrackTree)
    f->WriteTObject(fTrackTree);
  if (fTrackletTree)
    f->WriteTObject(fTrackletTree);
  f->Close();
  return kTRUE;
}

Bool_t AliTRDgtuSim::WriteTracksToESD(TList *ListOfTracks, AliESDEvent *esd) 
{
    if (esd) {
	TIter next(ListOfTracks);
	while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) next()) {
	    AliESDTrdTrack *trdtrack = trk->CreateTrdTrack();
	    esd->AddTrdTrack(trdtrack);
	    delete trdtrack;
	}
    }
    return kTRUE;
}
