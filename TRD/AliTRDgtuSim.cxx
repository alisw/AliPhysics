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

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliTreeLoader.h"
#include "AliLog.h"
#include "AliESDTrdTrack.h"

#include "AliTRDgtuSim.h"
#include "AliTRDgtuTMU.h"
#include "AliTRDtrackGTU.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDtrackletMCM.h"
#include "AliESDEvent.h"

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
  // destructor

  if (fTrackletArray)
    fTrackletArray->Delete();
  delete fTrackletArray;
  delete fTrackletTree;
}

Bool_t AliTRDgtuSim::RunGTUFromTrackletFile(TString filename, Int_t event, Int_t noev) 
{
  // run the GTU from a file of tracklets 
  // used for comparison to VHDL simulation

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
    
    AliDebug(5,"--------- Reading from file ----------");
    while (getline(input, str)) {
	lineno++;
	string = str;
	AliDebug(5,Form("Line %i : %s", lineno, string.Data()));
	
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
		TList *listOfTracks = new TList();
		fTMU->SetStack(iStackPrev);
		fTMU->SetSector(iSecPrev);
		fTMU->RunTMU(listOfTracks);
		AliDebug(1,Form("--- There are %i tracks. Writing ...", listOfTracks->GetEntries()));
		WriteTracksToTree(listOfTracks);
		fTMU->WriteTrackletsToTree(fTrackletTree);
		WriteTracksToDataFile(listOfTracks, iEventPrev);
		if (listOfTracks->GetEntries() > 0) 
		    AliDebug(2,Form("   %4.1f GeV/c", ((AliTRDtrackGTU*) listOfTracks->At(0))->GetPt() ));
		delete fTMU;
		fTMU = new AliTRDgtuTMU(); 
		delete listOfTracks;
		listOfTracks = 0x0;
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
	    AliDebug(2,Form("%i. tracklet: %s -> 0x%08x", i-4, ((TObjString*) tokens->At(i))->GetString().Data(), trackletWord));
	    AliTRDtrackletWord *trkl = new AliTRDtrackletWord(trackletWord);
	    fTMU->AddTracklet(trkl, iLink);
	}
    }
    
    if (fTMU && evcnt < noev) {
	TList *listOfTracks = new TList();
	fTMU->SetStack(iStackPrev);
	fTMU->SetSector(iSecPrev);
	fTMU->RunTMU(listOfTracks);
	WriteTracksToTree(listOfTracks);
	fTMU->WriteTrackletsToTree(fTrackletTree);
	WriteTracksToDataFile(listOfTracks, iEventPrev);
	delete fTMU;
	delete listOfTracks;
	fTMU = 0x0;
    }

    AliInfo(Form("Analyzed %i events", evcnt));
    return kTRUE; 
}

Bool_t AliTRDgtuSim::RunGTU(AliLoader *loader, AliESDEvent *esd) 
{
  // run the GTU on tracklets taken from the loader
  // if specified the GTU tracks are written to the ESD event 

    if (!LoadTracklets(loader)) {
	AliError("Could not load the tracklets. Nothing done ...");
	return kFALSE;
    }

    AliDebug(1, Form("running on %i tracklets", fTrackletArray->GetEntriesFast()));

    Int_t iStackPrev = -1;
    Int_t iSecPrev = -1;
    Int_t iSec = -1;
    Int_t iStack = -1;
    Int_t iLink = -1;

    if (fTMU) {
	delete fTMU;
	fTMU = 0x0;
    }
    
    TList *listOfTracks = new TList();
    
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
		fTMU->RunTMU(listOfTracks);
		WriteTracksToTree(listOfTracks);
		fTMU->WriteTrackletsToTree(fTrackletTree);
		WriteTracksToESD(listOfTracks, esd);
		fTMU->Reset();
		listOfTracks->Delete();
	    } else {
		fTMU = new AliTRDgtuTMU();
	    }
	    iStackPrev = iStack;
	    iSecPrev = iSec;
	}
	AliDebug(1, Form("adding tracklet: 0x%08x", trkl->GetTrackletWord()));
	fTMU->AddTracklet(trkl, iLink);
    }
    
    if (fTMU) {
	fTMU->SetStack(iStackPrev);
	fTMU->SetSector(iSecPrev);
	fTMU->RunTMU(listOfTracks);
	WriteTracksToTree(listOfTracks);
	fTMU->WriteTrackletsToTree(fTrackletTree);
	WriteTracksToESD(listOfTracks, esd);
	delete fTMU;
	fTMU = 0x0;
	listOfTracks->Delete();
    }

    delete listOfTracks;

    return kTRUE;
}

Bool_t AliTRDgtuSim::LoadTracklets(AliLoader *const loader) 
{
  // load the tracklets using the given loader

  AliDebug(1,"Loading tracklets ...");

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
  TTree *trackletTree = 0x0;

  // simulated tracklets
  trackletTree = trackletLoader->Tree();
  if (trackletTree) {
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
  }

  // raw tracklets
  AliTreeLoader *tl = (AliTreeLoader*) trackletLoader->GetBaseLoader("tracklets-raw");
  trackletTree = tl ? tl->Load(), tl->Tree() : 0x0;

  if (trackletTree) {
    if (!fTrackletArray)
      fTrackletArray = new TClonesArray("AliTRDtrackletWord", 1000);
    else if ((TClass::GetClass("AliTRDtrackletWord"))->InheritsFrom(fTrackletArray->Class()))
      fTrackletArray->Delete();
    else {
      fTrackletArray->Delete();
      delete fTrackletArray;
      fTrackletArray = new TClonesArray("AliTRDtrackletWord", 1000);
    }
    
    Int_t hc; 
    TClonesArray *ar = 0x0;
    trackletTree->SetBranchAddress("hc", &hc);
    trackletTree->SetBranchAddress("trkl", &ar);

    for (Int_t iEntry = 0; iEntry < trackletTree->GetEntries(); iEntry++) {
      trackletTree->GetEntry(iEntry);
      printf("%i tracklets in HC %i\n", ar->GetEntriesFast(), hc);
      for (Int_t iTracklet = 0; iTracklet < ar->GetEntriesFast(); iTracklet++) {
	AliTRDtrackletWord *trklWord = (AliTRDtrackletWord*) (*ar)[iTracklet];
	new((*fTrackletArray)[fTrackletArray->GetEntriesFast()]) AliTRDtrackletWord(trklWord->GetTrackletWord(), hc);
      }
    }
    return kTRUE;
  }
  
  AliError("No raw tracklet tree found\n");

  return kFALSE;
}

Bool_t AliTRDgtuSim::WriteTracksToDataFile(TList *listOfTracks, Int_t event) 
{
  // write the found tracks to a data file
  // used for comparison to VHDL simulation

    Int_t sm = 0;
    Int_t stack = 0;

    FILE *out;
    out = fopen("test.data", "a");

    AliDebug(1,Form("%i tracks found in event %i", listOfTracks->GetSize(), event));
    fprintf(out, "0 %5i %2i %i  00000000\n", event, sm, stack);
    for (Int_t i = 0; i < listOfTracks->GetSize(); i++) {
	AliTRDtrackGTU *trk = (AliTRDtrackGTU*) listOfTracks->At(i);
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

Bool_t AliTRDgtuSim::WriteTracksToTree(TList *listOfTracks, Int_t /*event*/) 
{
  // write the tracks to the tree for intermediate storage

  AliDebug(1,Form("Writing %i tracks to the tree...", listOfTracks->GetEntries()));

  if (!listOfTracks)
    return kFALSE;

  if (listOfTracks->GetEntries() <= 0) 
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

  TIter next(listOfTracks);
  while ((trk = (AliTRDtrackGTU*) next())) {
      trk->CookLabel();
      branch->SetAddress(&trk);
      fTrackTree->Fill();   
  }
  fTrackTree->ResetBranchAddress(branch);

  return kTRUE; 
}

Bool_t AliTRDgtuSim::WriteTreesToFile() const {
  // write the trees holding tracklets and tracks to file

  TFile *f = TFile::Open("TRD.GtuTracking.root", "RECREATE");
  f->cd();
  if (fTrackTree)
    f->WriteTObject(fTrackTree);
  if (fTrackletTree)
    f->WriteTObject(fTrackletTree);
  f->Close();
  return kTRUE;
}

Bool_t AliTRDgtuSim::WriteTracksToESD(const TList * const listOfTracks, AliESDEvent *esd) 
{
  // fill the found tracks to the given ESD event

    if (esd) {
	TIter next(listOfTracks);
	while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) next()) {
	    AliESDTrdTrack *trdtrack = trk->CreateTrdTrack();
	    esd->AddTrdTrack(trdtrack);
	    delete trdtrack;
	}
    }
    return kTRUE;
}

Bool_t AliTRDgtuSim::WriteTracksToLoader()
{
  // write the GTU tracks to the dedicated loader
  // these tracks contain more information than the ones in the ESD

  if (!fTrackTree) {
    AliError("No track tree found!");
    return kFALSE;
  }

  AliRunLoader *rl = AliRunLoader::Instance();
  AliDataLoader *dl = 0x0;
  if (rl)
    dl = rl->GetLoader("TRDLoader")->GetDataLoader("gtutracks");
  if (!dl) {
    AliError("Could not get the GTU-track data loader!");
    return kFALSE;
  }

  TTree *trackTree = dl->Tree();
  if (!trackTree) {
    dl->MakeTree();
    trackTree = dl->Tree();
  }
  
  AliTRDtrackGTU *trk = 0x0;
  if (!trackTree->GetBranch("TRDtrackGTU"))
    trackTree->Branch("TRDtrackGTU", "AliTRDtrackGTU", &trk, 32000);
  
  AliDebug(1,Form("Found %d tracks", fTrackTree->GetEntries()));

  for (Int_t iTrack = 0; iTrack < fTrackTree->GetEntries(); iTrack++) {
    fTrackTree->SetBranchAddress("TRDgtuTrack", &trk);
    fTrackTree->GetEntry(iTrack);
    trackTree->SetBranchAddress("TRDtrackGTU", &trk);
    trackTree->Fill();
  }

  dl->WriteData("OVERWRITE");

  return kTRUE;
}
