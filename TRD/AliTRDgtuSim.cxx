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
#include "TObject.h"
#include "TClonesArray.h"

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliTreeLoader.h"
#include "AliLog.h"
#include "AliESDTrdTrack.h"
#include "AliESDTrdTracklet.h"

#include "AliTRDgtuSim.h"
#include "AliTRDfeeParam.h"
#include "AliTRDgtuTMU.h"
#include "AliTRDtrackGTU.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDtrackletMCM.h"
#include "AliESDEvent.h"

ClassImp(AliTRDgtuSim)

AliTRDgtuSim::AliTRDgtuSim(AliRunLoader *rl)
  : TObject(),
  fRunLoader(rl),
  fFeeParam(AliTRDfeeParam::Instance()),
  fTMU(0x0),
  fTrackletArray(0x0)
{

}

AliTRDgtuSim::~AliTRDgtuSim()
{
  // destructor

  if (fTrackletArray)
    fTrackletArray->Clear();
  delete fTrackletArray;
}

Bool_t AliTRDgtuSim::RunGTUFromTrackletFile(TString filename, Int_t event, Int_t noev)
{
  // run the GTU from a file of tracklets
  // used for comparison to VHDL simulation

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

  TClonesArray trklArray("AliTRDtrackletWord", 100);
  TClonesArray trklArrayGTU("AliTRDtrackletGTU", 100);

  AliDebug(1, Form("--------- Reading from %s ----------", filename.Data()));
  while (getline(input, str)) {
    lineno++;
    string = str;

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

    if ((iEvent != iEventPrev) ||
	(iStack != iStackPrev) ||
	(iSec != iSecPrev)) {
      if(fTMU) {
	TList *listOfTracks = new TList();
	fTMU->SetStack(iStackPrev);
	fTMU->SetSector(iSecPrev);
	fTMU->RunTMU(listOfTracks);
	AliDebug(1,Form("--- There are %i tracks. Writing ...", listOfTracks->GetEntries()));
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
      sscanf(((TObjString*) tokens->At(i))->GetString().Data(), "%i", &trackletWord);
      if (trackletWord == 0x10001000)
	break;
      AliDebug(2, Form("link: %2i trkl: %2i - %s -> 0x%08x",
		       iLink, i-4, ((TObjString*) tokens->At(i))->GetString().Data(), trackletWord));
      AliTRDtrackletWord *tracklet = new (trklArray[trklArray.GetEntriesFast()])       AliTRDtrackletWord(trackletWord);
      AliTRDtrackletGTU   *trkl    = new (trklArrayGTU[trklArrayGTU.GetEntriesFast()]) AliTRDtrackletGTU(tracklet);
      if (fTMU)
	fTMU->AddTracklet(trkl, iLink);
    }
  }

  if (fTMU && evcnt < noev) {
    TList *listOfTracks = new TList();
    fTMU->SetStack(iStackPrev);
    fTMU->SetSector(iSecPrev);
    fTMU->RunTMU(listOfTracks);
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

  if (!fFeeParam->GetTracklet())
    return kFALSE;

  if (fTrackletArray)
    fTrackletArray->Clear();

  if (loader) {
    if (!LoadTracklets(loader)) {
	AliError("Could not load the tracklets. Nothing done ...");
	return kFALSE;
    }
  }
  else {
    LoadTracklets(esd);
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

    while (AliTRDtrackletGTU *trkl = (AliTRDtrackletGTU*) next()) {
	iSec = trkl->GetDetector() / 30;
	iStack = (trkl->GetDetector() % 30) / 6;
	iLink = trkl->GetHCId() % 12;

	if (iStack != iStackPrev || iSec != iSecPrev) {
	    if(fTMU) {
		fTMU->SetStack(iStackPrev);
		fTMU->SetSector(iSecPrev);
		fTMU->RunTMU(listOfTracks);
		WriteTracksToLoader(listOfTracks);
		WriteTracksToESD(listOfTracks, esd);
		fTMU->Reset();
		listOfTracks->Delete();
	    } else {
		fTMU = new AliTRDgtuTMU();
	    }
	    iStackPrev = iStack;
	    iSecPrev = iSec;
	    AliDebug(1, Form("now in sec %i, stack %i", iSec, iStack));
	}
	AliDebug(1, Form("adding tracklet: 0x%08x in sec %i stack %i link %i",
			 trkl->GetTrackletWord(), trkl->GetDetector() / 30, (trkl->GetDetector() % 30) / 6, trkl->GetHCId() % 12));
	if (fTMU) {
	  fTMU->AddTracklet(trkl, iLink);
	}
    }

    if (fTMU) {
	fTMU->SetStack(iStackPrev);
	fTMU->SetSector(iSecPrev);
	fTMU->RunTMU(listOfTracks);
	WriteTracksToLoader(listOfTracks);
	WriteTracksToESD(listOfTracks, esd);
	delete fTMU;
	fTMU = 0x0;
	listOfTracks->Delete();
    }

    delete listOfTracks;

    return kTRUE;
}

Bool_t AliTRDgtuSim::LoadTracklets(const AliESDEvent *const esd)
{
  AliDebug(1,"Loading tracklets from ESD event ...");

  if (!fTrackletArray)
    fTrackletArray = new TClonesArray("AliTRDtrackletGTU", 1000);

  for (Int_t iTracklet = 0; iTracklet < esd->GetNumberOfTrdTracklets(); iTracklet++) {
    AliESDTrdTracklet *trkl = esd->GetTrdTracklet(iTracklet);
    new ((*fTrackletArray)[fTrackletArray->GetEntries()]) AliTRDtrackletGTU(trkl);
  }

  return kTRUE;
}

Bool_t AliTRDgtuSim::LoadTracklets(AliLoader *const loader)
{
  // load the tracklets using the given loader

  AliDebug(1,"Loading tracklets ...");

  if (!fFeeParam->GetTracklet())
    return kFALSE;

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
	fTrackletArray = new TClonesArray("AliTRDtrackletGTU", 1000);

      AliTRDtrackletMCM *trkl = 0x0;
      trklbranch->SetAddress(&trkl);
      for (Int_t iTracklet = 0; iTracklet < trklbranch->GetEntries(); iTracklet++) {
	trklbranch->GetEntry(iTracklet);
	new ((*fTrackletArray)[fTrackletArray->GetEntries()]) AliTRDtrackletGTU(new AliTRDtrackletMCM(*trkl));
	((AliTRDtrackletGTU *)((*fTrackletArray)[fTrackletArray->GetEntries()-1]))->SetMCMtrackletIndex(iTracklet);
      }
      return kTRUE;
    }
  }

  // raw tracklets
  AliTreeLoader *tl = (AliTreeLoader*) trackletLoader->GetBaseLoader("tracklets-raw");
  trackletTree = tl ? tl->Load(), tl->Tree() : 0x0;

  if (trackletTree) {
    if (!fTrackletArray)
      fTrackletArray = new TClonesArray("AliTRDtrackletGTU", 1000);

    Int_t hc;
    TClonesArray *ar = 0x0;
    trackletTree->SetBranchAddress("hc", &hc);
    trackletTree->SetBranchAddress("trkl", &ar);

    for (Int_t iEntry = 0; iEntry < trackletTree->GetEntries(); iEntry++) {
      trackletTree->GetEntry(iEntry);
      AliDebug(2, Form("%i tracklets in HC %i", ar->GetEntriesFast(), hc));
      for (Int_t iTracklet = 0; iTracklet < ar->GetEntriesFast(); iTracklet++) {
	AliTRDtrackletWord *trklWord = (AliTRDtrackletWord*) (*ar)[iTracklet];
	new((*fTrackletArray)[fTrackletArray->GetEntriesFast()]) AliTRDtrackletGTU(new AliTRDtrackletWord(trklWord->GetTrackletWord(), hc));
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
    // fprintf(out, "0 %5i %2i %i  00000000\n", event, sm, stack);
    for (Int_t i = 0; i < listOfTracks->GetSize(); i++) {
	AliTRDtrackGTU *trk = (AliTRDtrackGTU*) listOfTracks->At(i);
	sm = trk->GetSector();
	stack = trk->GetStack();

	ULong64_t trackWord = 1;
	AppendBits(trackWord,   1, 0);
	AppendBits(trackWord,   6, trk->GetTrackletMask());
	AppendBits(trackWord,  18, (Int_t) trk->GetA());
	AppendBits(trackWord,  18, (Int_t) trk->GetB());
	AppendBits(trackWord,  12, (Int_t) trk->GetC());
	AppendBits(trackWord,   8, trk->GetPID());
	fprintf(out, "ev. %i sec. %i stack %i - track word: 0x%016llx, ",
		event, sm, stack, trackWord);

	trackWord = 0;
	AppendBits(trackWord, 11, 0); // flags
	AppendBits(trackWord,  3, 0);
	AppendBits(trackWord, 13, trk->GetYapprox());
	AppendBits(trackWord,  6, trk->GetTrackletIndex(5));
	AppendBits(trackWord,  6, trk->GetTrackletIndex(4));
	AppendBits(trackWord,  6, trk->GetTrackletIndex(3));
	AppendBits(trackWord,  6, trk->GetTrackletIndex(2));
	AppendBits(trackWord,  6, trk->GetTrackletIndex(1));
	AppendBits(trackWord,  6, trk->GetTrackletIndex(0));
	fprintf(out, "extended track word: 0x%016llx\n", trackWord);

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

Bool_t AliTRDgtuSim::WriteTracksToESD(const TList * const listOfTracks, AliESDEvent *esd)
{
  // fill the found tracks to the given ESD event

    if (esd) {
	TIter next(listOfTracks);
	while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) next()) {
	    AliESDTrdTrack *trdtrack = trk->CreateTrdTrack();
	    if (trdtrack->GetLabel() < 0)
	      trdtrack->SetLabel(-2);
	    esd->AddTrdTrack(trdtrack);
	    delete trdtrack;
	}
    }
    return kTRUE;
}

Bool_t AliTRDgtuSim::WriteTracksToLoader(const TList * const listOfTracks)
{
  // write the GTU tracks to the dedicated loader
  // these tracks contain more information than the ones in the ESD

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

  AliDebug(1, Form("Writing %i tracks to loader", listOfTracks->GetEntries()));
  TIter next(listOfTracks);
  while ((trk = (AliTRDtrackGTU*) next())) {
    trackTree->SetBranchAddress("TRDtrackGTU", &trk);
    trackTree->Fill();
  }

  dl->WriteData("OVERWRITE");

  return kTRUE;
}
