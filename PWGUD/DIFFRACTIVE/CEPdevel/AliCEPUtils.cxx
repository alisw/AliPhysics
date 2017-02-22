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
//
// AliCEPUtils
// for
// AliAnalysisTaskCEP
//
//
//  Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  rewritten by
//  Paul Buehler <paul.buehler@oeaw.ac.at>


#include "AliCEPUtils.h"

//------------------------------------------------------------------------------
AliCEPUtils::AliCEPUtils():
    fTPCnclsS()
  , fTrackDCAz()
  , fTrackEtaMin()
  , fTrackEtaMax()
  , fTrackCutListPrim(0x0)
{
  
  // initialize the track cuts
  fTrackCutListPrim = new TList();
	fTrackCutListPrim->SetOwner();
  
}

//------------------------------------------------------------------------------
AliCEPUtils::~AliCEPUtils()
{

  if (fTrackCutListPrim) {
    fTrackCutListPrim->Delete();
    delete fTrackCutListPrim;
    fTrackCutListPrim = 0x0;
  }

}

//------------------------------------------------------------------------------
TH1F* AliCEPUtils::GetHistStatsFlow()
{
	// setup the stats flow histogram

	TH1F *hist = new TH1F("statsFlow","",
    AliCEPBase::kBinLastValue,0,AliCEPBase::kBinLastValue);
	TAxis* axis = hist->GetXaxis();

	axis->SetBinLabel(AliCEPBase::kBinTotalInput+1,   "total input");
	axis->SetBinLabel(AliCEPBase::kBinGoodInput+1,    "good input");
	axis->SetBinLabel(AliCEPBase::kBinMCEvent+1,      "MC");
	axis->SetBinLabel(AliCEPBase::kBinPhysEvent+1,    "physics event");
	axis->SetBinLabel(AliCEPBase::kBinEventCut+1,     "passed event cut");
	axis->SetBinLabel(AliCEPBase::kBinPhysel+1,       "Physics selected");
	axis->SetBinLabel(AliCEPBase::kBinPileup+1,       "pileup");
	axis->SetBinLabel(AliCEPBase::kBinClusterCut+1,   "passed cluster cut");
	axis->SetBinLabel(AliCEPBase::kBinDGTrigger+1,    "DG trigger");
	axis->SetBinLabel(AliCEPBase::kBinSharedCluster+1,"passed shared cluster test");
	axis->SetBinLabel(AliCEPBase::kBinVtx+1,          "Vtx ok");
	axis->SetBinLabel(AliCEPBase::kBinMBOR+1,         "MBOR");
	axis->SetBinLabel(AliCEPBase::kBinMBAND+1,        "MBAND");
	axis->SetBinLabel(AliCEPBase::kBinnoV0+1,         "!V0");
	axis->SetBinLabel(AliCEPBase::kBinnoFMD+1,        "!FMD");
	axis->SetBinLabel(AliCEPBase::kBinnoAD+1,         "!AD");
	axis->SetBinLabel(AliCEPBase::kBinSaved+1,        "saved");

	return hist;

}

//------------------------------------------------------------------------------
Int_t AliCEPUtils::GetEventType(const AliVEvent *Event)
{
	// checks of which type a event is:
  // A:  CINT1A-ABCE-NOPF-ALL         : beam from A-side
  // C:  CINT1C-ABCE-NOPF-ALL         : beam from C-side
  // E:  CINT1-E-NOPF-ALL, CDG5-E     : no beam
  // I:  CINT1B-ABCE-NOPF-ALL, CDG5-I : beam from both sides
  // AC: CDG5-AC                      : beam from one side, side not specified

	TString firedTriggerClasses = Event->GetFiredTriggerClasses();
  // printf("<I - AliCEPUtils::GetEventType> firedTriggerClasses: %s\n",
  //   firedTriggerClasses.Data());

	if (firedTriggerClasses.Contains("CINT1A-ABCE-NOPF-ALL")) { // A
		return AliCEPBase::kBinEventA;
	}
	if (firedTriggerClasses.Contains("CINT1C-ABCE-NOPF-ALL")) { // C
		return AliCEPBase::kBinEventC;
	}
	if (firedTriggerClasses.Contains("CINT1B-ABCE-NOPF-ALL")) { // I
		return AliCEPBase::kBinEventI;
	}
	if (firedTriggerClasses.Contains("CINT1-E-NOPF-ALL")) { // E
		return AliCEPBase::kBinEventE;
	}
	if (firedTriggerClasses.Contains("CDG5-E")) { // E
		return AliCEPBase::kBinEventE;
	}
	if (firedTriggerClasses.Contains("CDG5-I")) { // I
		return AliCEPBase::kBinEventI;
	}
	if (firedTriggerClasses.Contains("CDG5-AC")) { // AC
		return AliCEPBase::kBinEventAC;
	}
	return AliCEPBase::kBinEventUnknown;
  
}

//------------------------------------------------------------------------------
// this is a copy of the method described on
// https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
UInt_t AliCEPUtils::GetVtxPos(AliVEvent *Event, TVector3 *fVtxPos)
{
    
  // initialize
  UInt_t fVtxType = AliCEPBase::kVtxUnknown;
  fVtxPos->SetXYZ(-999.9,-999.9,-999.9);

  // On AOD, only one primary vertex is stored: aodEv->GetPrimaryVertex()
  if (Event->GetDataLayoutType()==AliVEvent::kAOD) {
    
    fVtxType |= AliCEPBase::kVtxAOD;
  
  } else {
  
    const AliESDVertex *trkVertex = ((AliESDEvent*)Event)->GetPrimaryVertexTracks();
    const AliESDVertex *spdVertex = ((AliESDEvent*)Event)->GetPrimaryVertexSPD();
    Bool_t hasSPD = spdVertex->GetStatus();
    Bool_t hasTrk = trkVertex->GetStatus();
  
    //Note that AliVertex::GetStatus checks that N_contributors is > 0
    //reject events if both are explicitly requested and none is available
    if (!(hasSPD && hasTrk)) return AliCEPBase::kVtxUnknown;
  
    // check the spd vertex resolution and reject if not satisfied
    if (hasSPD) {
      if (spdVertex->IsFromVertexerZ() &&
        !(spdVertex->GetDispersion()<0.04 &&
        spdVertex->GetZRes()<0.25)) return AliCEPBase::kVtxUnknown;
    }
  
    // reject events if none between the SPD or track verteces are available
    // if no trk vertex, try to fall back to SPD vertex;
    if (hasTrk) {
      fVtxType |= AliCEPBase::kVtxTracks;
      if (hasSPD) {
        fVtxType |= AliCEPBase::kVtxSPD;
        // check the proximity between the spd vertex and trak vertex, and reject if not satisfied
        if (TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5) return AliCEPBase::kVtxUnknown;
      }
      
    } else {
      fVtxType |= AliCEPBase::kVtxSPD;
    }
  
  }
  
  // Cut on the vertex z position
  const AliVVertex *vertex = Event->GetPrimaryVertex();
  if (TMath::Abs(vertex->GetZ())>10) return AliCEPBase::kVtxUnknown;
  
  // set the vertex position fVtxPos
  fVtxPos->SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());
  
  return fVtxType;
  
}

//------------------------------------------------------------------------------
// This routine takes an ESD event as input
// loops over all tracks and assigns them a status word (fTrackStatus)
// according to the result of various tests
// fTracks is at output an array of AliESDTracks
Int_t AliCEPUtils::AnalyzeTracks(AliESDEvent* fESDEvent,
  TObjArray* fTracks,TArrayI* fTrackStatus)
{

  // Initialisations
  UInt_t trackstat;
  Int_t nTracks = fESDEvent->GetNumberOfTracks();

  // sets event of all tracks to esd
  fESDEvent->ConnectTracks();

  // reset fTracks and fTrackStatus
  fTracks->Clear();
  fTrackStatus->Reset();
  fTrackStatus->Set(nTracks);
  if (nTracks == 0) return nTracks;

  // prepare for pplication of track cuts
  AliESDtrackCuts *cut;
  
  // prepare for the FiredChips test
  const AliMultiplicity *mult = fESDEvent->GetMultiplicity();
  Int_t statusLay;
  Int_t idet = -1;
  Float_t xloc,zloc;
  UInt_t Modules[2],eq,hs,chip;
  Bool_t goodtrack,ktmp;
  
  // get the fired chips
  TArrayI *fFiredChips = new TArrayI(1200);
  Int_t nFiredChips = 0;
  for (Int_t ii=0; ii<1200; ii++) {
    // was chip fired?
    if (!mult->TestFiredChipMap(ii)) continue;
    
    // get the module of the fired chip
    AliSPDUtils::GetOnlineFromOfflineChipKey(ii,eq,hs,chip);
    fFiredChips->SetAt(AliSPDUtils::GetOfflineModuleFromOnline(eq,hs,chip),nFiredChips);
    nFiredChips++;
  }
  fFiredChips->Set(nFiredChips);
  
  // retrieve main vertex
  const AliVVertex *vtx = fESDEvent->GetPrimaryVertex();
  
  // prepare magnetic field for propagation to DCA
  Double_t bfield = fESDEvent->GetMagneticField();
  Double_t dca[2], cov[3];

  // initialize V0 daughters    
  TArrayI *v0daughters = new TArrayI(nTracks);
  v0daughters->Reset(kFALSE);
  Int_t nV0  = fESDEvent->GetNumberOfV0s();
  for (Int_t ii=0; ii<nV0; ii++) {
    AliESDv0 *V0 = fESDEvent->GetV0(ii);
    if (!V0) continue;

    // positive and negative daughters
    v0daughters->SetAt(kTRUE,V0->GetPindex());
    v0daughters->SetAt(kTRUE,V0->GetNindex());
  }
  
  // loop over all tracks of fESDEvent
  // determine the TrackStatus
  AliESDtrack* track = NULL;
  for (Int_t ii=0; ii<nTracks; ii++) {
    track = (AliESDtrack*) fESDEvent->GetTrack(ii);
    track->SetESDEvent(fESDEvent);
    if (!track) continue;
    
    // add track to buffer
    fTracks->Add(track);
    
    // go through the list of selection tests
    // and update the TrackStatus word trackstat accordingly
    // see AliCEPBase.h for a definition of the TrackStatus word bits
    // initialize trackstat
    trackstat = AliCEPBase::kTTBaseLine;
    
    // TOFBunchCrossing
    if (track->GetTOFBunchCrossing() ==0) ;
      trackstat |=  AliCEPBase::kTTTOFBunchCrossing;

    // number of TPC shared clusters <= fTPCnclsS(3)    
    if (track->GetTPCnclsS() <= fTPCnclsS)
      trackstat |=  AliCEPBase::kTTTPCScluster;
    
    // DCA to vertex is < 500
    if (vtx) {
      Bool_t DCAok =
        track->PropagateToDCA(vtx,bfield,fTrackDCA,dca,cov); 
       if (DCAok)
        trackstat |=  AliCEPBase::kTTDCA;
    }
    
    // is a daughter of a V0
    if (v0daughters->At(ii))
      trackstat |= AliCEPBase::kTTV0;

    // is an ITS pure track
    if (track->GetStatus() & AliESDtrack::kITSpureSA)
      trackstat |=  AliCEPBase::kTTITSpure;

    // |Zv-VtxZ| <= fTrackDCAz(6)
    if (vtx) {
      Bool_t DCAzok = TMath::Abs(track->Zv() - vtx->GetZ()) <= fTrackDCAz;
      if (DCAzok) 
        trackstat |= AliCEPBase::kTTZv;
    }
    
    // fTrackEtaMin(-0.9)<eta<fTrackEtaMax(0.9)    
    if (track->Eta() >= fTrackEtaMin && track->Eta() <= fTrackEtaMax)
      trackstat |= AliCEPBase::kTTeta;

    // accepted by ITSTPC and ITSSA criteria
    if (cut = (AliESDtrackCuts*)fTrackCutListPrim->At(0))
    {
      if (cut->AcceptTrack(track)) {
        trackstat |= AliCEPBase::kTTAccITSTPC;
	    }
    }
    if (cut = (AliESDtrackCuts*)fTrackCutListPrim->At(1))
    {
      if (cut->AcceptTrack(track))
        trackstat |= AliCEPBase::kTTAccITSSA;
	  }
    
    // FiredChips test
    // test whether both modules associated with the track are modules
    // with fired chips
    // get the 2 modules associated with the track
    Bool_t retc = track->GetITSModuleIndexInfo(0,idet,statusLay,xloc,zloc);
    if (retc && statusLay!=5) Modules[0] = idet;
    retc = track->GetITSModuleIndexInfo(1,idet,statusLay,xloc,zloc);
    if (retc && statusLay!=5) Modules[1] = idet;
    
    // check whether both modules are modules of fired chips
    goodtrack = kTRUE;
    for (Int_t iM = 0; iM<2; iM++) {
      ktmp = kFALSE;
      for (Int_t jj=0; jj<nFiredChips; jj++) {
        if (fFiredChips->At(jj) == Modules[iM]) {
          ktmp = kTRUE;
          continue;
        }
      }

      if (!ktmp) {
        goodtrack = kFALSE;
        // leave the for lopp
        break;
      }
    }
    if (goodtrack)
      trackstat |= AliCEPBase::kTTFiredChips;
    
    // save trackstat
    fTrackStatus->AddAt(trackstat,ii);
  }
  
  fFiredChips->Reset();
  delete fFiredChips;
  v0daughters->Reset();
  delete v0daughters;

  return nTracks;

}

// ------------------------------------------------------------------------------
Bool_t AliCEPUtils::checkstatus(UInt_t stat, UInt_t mask, UInt_t pattern)
{
  
  // mask defines which bits of stat and pattern must agree
  return (stat & mask) == pattern;

}

// ------------------------------------------------------------------------------
Int_t AliCEPUtils::countstatus(TArrayI *stats, UInt_t mask, UInt_t pattern)
{
  
  // initialisations
  Int_t ngood = 0;
  
  for (Int_t ii=0; ii<stats->GetSize(); ii++) {
    if (checkstatus(stats->At(ii),mask,pattern)) ngood++;
  }
  
  return ngood;

}

// ------------------------------------------------------------------------------
Int_t AliCEPUtils::countstatus(TArrayI *stats, UInt_t mask, UInt_t pattern, TArrayI* indices)
{

  // initialisations
  Int_t ngood = 0;
  indices->Reset(0);
  indices->Set(stats->GetSize());
  
  for (Int_t ii=0; ii<stats->GetSize(); ii++) {
    if (checkstatus(stats->At(ii),mask,pattern)) {
      
      // increment the counter
      ngood++;
    
      // update the indices
      indices->SetAt(ii,ngood-1);
      // printf("%i ",ii);
    }
  }
  indices->Set(ngood);
  // printf("\n");
  
  return ngood;

}

// ------------------------------------------------------------------------------
Int_t AliCEPUtils::countstatus(TArrayI *stats,
  TArrayI *masks, TArrayI *patterns)
{

  // checks status words with masks and patterns
  // makes logical OR of all tests
  
  // initialisations
  Bool_t ktmp;
  Int_t ngood = 0;

  // number of masks and patterns must agree
  Int_t ntests = masks->GetSize();
  if (ntests != patterns->GetSize()) return ngood;
  Int_t nstats = stats->GetSize();

  for (Int_t ii=0; ii<nstats; ii++) {
    
    ktmp = kFALSE;
    for (Int_t jj=0; jj<ntests;jj++) {
      
      if (checkstatus(stats->At(ii),masks->At(jj),patterns->At(jj))) {
        ktmp = kTRUE;
        break;
      }
    }
      
    if (ktmp) ngood++;

  }
  
  return ngood;

}

// ------------------------------------------------------------------------------
Int_t AliCEPUtils::countstatus(TArrayI *stats,
  TArrayI *masks, TArrayI *patterns, TArrayI* indices)
{

  // checks status words with masks and patterns
  // makes logical OR of all tests
  
  // initialisations
  Bool_t ktmp;
  Int_t ngood = 0;
  indices->Reset(0);

  // number of masks and patterns must agree
  Int_t ntests = masks->GetSize();
  if (ntests != patterns->GetSize()) return ngood;
  Int_t nstats = stats->GetSize();
  indices->Set(nstats);

  for (Int_t ii=0; ii<nstats; ii++) {
    
    ktmp = kFALSE;
    for (Int_t jj=0; jj<ntests;jj++) {
      
      if (checkstatus(stats->At(ii),masks->At(jj),patterns->At(jj))) {
        ktmp = kTRUE;
        break;
      }
    }
      
    if (ktmp) {
      // increment the counter
      ngood++;
    
      // update the indices
      indices->SetAt(ii,ngood-1);
      // printf("%i ",ii);
      
    }

  }
  indices->Set(ngood);
  // printf("\n");
  
  return ngood;

}

// ------------------------------------------------------------------------------
Int_t AliCEPUtils::GetCEPTracks(
  AliESDEvent *ESDEvent, TArrayI *stats, TArrayI* indices)
{

  // now use the stats to scrutinize the event and tracks
  // e.g. Martin's selection
  // 1. no track with !kTTTPCScluster
  // 2. kTTITSpure -> npureITSTracks
  // 3. kTTDCA && !kTTV0 && !kTTITSpure && kTTZv -> nTrackSel
  // 4. 3. && kTTeta && (kTTAccITSTPC || kTTAccITSSA) -> nTrackAccept
  // 5. nTrackSel>=npureITSTracks && nTrackSel>=nTracklets && nTrackAccept==nTrackSel
  // 6. pass fpassedFiredChipsTest
  //  
  // if all criteria are met the event is accepted and the number of
  // good tracks is nTrackSel = nMartinSel
  //
  // possible return values:
  //  >0: number of selected CEP tracks
  //  -1: tracks with !kTTTPCScluster
  //  -2: more tracklets or ITSpure tracks than selected tracks
  //  -3: one of the selected tracks has missed a further tests
  //  -4: not all fired chips are associated with a selected track
  
  // initialisations
  UInt_t mask, pattern;
  TArrayI *masks    = new TArrayI();
  TArrayI *patterns = new TArrayI();
  indices->Reset();

  // nbad track with !kTTTPCScluster
  mask = AliCEPBase::kTTTPCScluster;
  pattern = 0;
  Int_t nbad = countstatus(stats,mask,pattern);
  if (nbad > 0) return -1;
    
  // the number of tracks which are considered for further testing
  mask = AliCEPBase::kTTDCA+AliCEPBase::kTTV0+AliCEPBase::kTTITSpure+
    AliCEPBase::kTTZv;
  pattern =  AliCEPBase::kTTDCA+AliCEPBase::kTTZv;
  Int_t nTrackSel = countstatus(stats,mask,pattern,indices);

  // the number of ITSpure tracks
  mask = AliCEPBase::kTTITSpure;
  pattern = AliCEPBase::kTTITSpure;
  Int_t npureITSTracks = countstatus(stats,mask,pattern);
    
  // nTrackSel>=npureITSTracks && nTrackSel>=nTracklets
  const AliMultiplicity *mult = ESDEvent->GetMultiplicity();
  Int_t nTracklets = mult->GetNumberOfTracklets();
  if (nTrackSel<npureITSTracks || nTrackSel<nTracklets) return -2;
  
  // further track tests
  Int_t nTrackAccept;
  mask += AliCEPBase::kTTeta;
  pattern += AliCEPBase::kTTeta;
  masks->Set(2);
  patterns->Set(2);
  masks->AddAt(mask+AliCEPBase::kTTAccITSTPC,0);
  masks->AddAt(mask+AliCEPBase::kTTAccITSSA, 1);
  patterns->AddAt(pattern+AliCEPBase::kTTAccITSTPC,0);
  patterns->AddAt(pattern+AliCEPBase::kTTAccITSSA, 1);
  nTrackAccept = countstatus(stats,masks,patterns,indices);
  if (nTrackAccept<nTrackSel) return -3;
  
  // FiredChips test
  Bool_t fpassedFiredChipsTest = TestFiredChips(ESDEvent,indices);
  if (!fpassedFiredChipsTest) return -4;

  printf("<I - UserExec> nBad            : %i\n",nbad);
  printf("<I - UserExec> nTracklets      : %i\n",nTracklets);
  printf("<I - UserExec> npureITSTracks  : %i\n",npureITSTracks);
  printf("<I - UserExec> nTrackSel       : %i\n",nTrackSel);
  printf("<I - UserExec> FiredChipsTest  : %i\n",fpassedFiredChipsTest);
  printf("<I - UserExec> nTrackAccept    : %i\n",nTrackAccept);

  // clean up
  if (masks) {
    delete masks;
    masks = 0x0;
  }
  if (patterns) {
    delete patterns;
    patterns = 0x0;
  }
  if (indices) {
    delete indices;
    indices = 0x0;
  }
  
}

// ------------------------------------------------------------------------------
Int_t AliCEPUtils::GetResiduals(AliESDEvent* fESDEvent)
{

	// determines the number of tracklets in an event, which are not 
  // associated with a track

	// initialisations
  Int_t nResiduals = 0;
  Int_t id1 = -1, id2 = -1;

	const AliMultiplicity *mult = fESDEvent->GetMultiplicity();
	if (mult) {
  
		for (Int_t ii = 0; ii < mult->GetNumberOfTracklets(); ii++) {
			if (!mult->GetTrackletTrackIDs(ii,0,id1,id2))
        nResiduals++;
		}
	}
  
  return nResiduals;

}

// ------------------------------------------------------------------------------
Bool_t AliCEPUtils::TestFiredChips(AliESDEvent *esd, TArrayI *indices)
{
  
  // All fired chips must be associated with a track listed in indices
  // If there is at least one fired chip without associated track, then
  // return kFALSE 
  // If no chips are fired then return kFALSE
  
  // initialisations
  Bool_t retc;
  UInt_t eq,hs,chip,module;
  Int_t statusLay,idet = -1;
  Float_t xloc,zloc;
  Int_t nFiredChips = 0;
  
  const AliMultiplicity *mult = esd->GetMultiplicity();
  Int_t Ntracks = indices->GetSize();
  UInt_t *Modules = new UInt_t[2*Ntracks];

  // get the modules associated with a selected track
  for (Int_t ii=0; ii<Ntracks; ii++)
  {
    idet = -1;
    AliESDtrack* track = esd->GetTrack(indices->At(ii));
    retc = track->GetITSModuleIndexInfo(0,idet,statusLay,xloc,zloc);
    if (retc && statusLay!=5) Modules[2*ii] = idet;
    retc = track->GetITSModuleIndexInfo(1,idet,statusLay,xloc,zloc);
    if (retc && statusLay!=5) Modules[2*ii+1] = idet;
  }

  // check that all fired chips are assosiated with one selected track
  for (Int_t ii=0; ii<1200; ii++)
  {
    if (!mult->TestFiredChipMap(ii)) continue;
    nFiredChips++;
    
    AliSPDUtils::GetOnlineFromOfflineChipKey(ii,eq,hs,chip);
    module = AliSPDUtils::GetOfflineModuleFromOnline(eq,hs,chip);

    Bool_t ktmp = kFALSE;
    for (Int_t iM = 0; iM<2*Ntracks; iM++)
    {
      if (Modules[iM]==module) ktmp = kTRUE;
    }
    if(!ktmp) 
    {
      // printf("module without track %i\n",module);
      delete[] Modules;
      return kFALSE;
    }
  }

  delete[] Modules;
  if (nFiredChips>0) return kTRUE;
  return kFALSE;

}

//------------------------------------------------------------------------------
void AliCEPUtils::SPDLoadGeom(Int_t run)
{
	// method to get the gGeomanager
	// it is called at the CreatedOutputObject stage
	// to comply with the CAF environment

	AliCDBManager *man = AliCDBManager::Instance();

	TString cdbpath;
	if (man->IsDefaultStorageSet()) {
		const AliCDBStorage *dsto = man->GetDefaultStorage();
		cdbpath = TString(dsto->GetBaseFolder());
		//printf("default was set to: %s\n", cdbpath.Data());
	}
	else { //should not be used!
		// man->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
		// would be needed on grid
		man->SetDefaultStorage(gSystem->Getenv("TRAIN_CDB_PATH"));
		cdbpath = TString(gSystem->Getenv("TRAIN_CDB_PATH"));
	}

	man->SetSpecificStorage("ITS/Align/Data",cdbpath);
	man->SetSpecificStorage("GRP/Geometry/Data",cdbpath);
	man->SetRun(run);

	AliCDBEntry* obj = man->Get(AliCDBPath("GRP", "Geometry", "Data"));
	if (!obj) {
		printf("AliCEPUtils failed loading geometry object for run %i\n",run);
		return;
	}
	AliGeomManager::SetGeometry((TGeoManager*)obj->GetObject());
	AliGeomManager::ApplyAlignObjsFromCDB("ITS");
}

//------------------------------------------------------------------------------
void AliCEPUtils::DetermineMCprocessType (
  AliMCEvent *fMCEvent, TString fMCGenerator, Int_t &fMCProcess)
{
	//
	// retrieves the MC process type from the AliGenEventHeader and classifies
	// them
	//

	// get MC information
	fMCProcess = AliCEPBase::kProctypeUnknown;
	Int_t fMCProcessType = AliCEPBase::kProctypeUnknown;

	if (fMCEvent) {
		AliGenEventHeader* header = fMCEvent->GenEventHeader();
    
    if (header) {
			// cover all possible generators
      //
      // this is the list of generators with a specific header
      // status 10/07/2016
      // Generator  Name in header
      // cocktail
      // DPMjet
      // Epos3
      // Epos
      // GeVSim
      // HepMC
      // Herwig
      // Hijing
      // Pythia     Pythia
      // Toy
      // 
      // this is the list of generators without specific header
      // DIME       Dime
      // Starlight
		  
      // get the name of this generator
      fMCGenerator = TString(header->GetName());
      //printf("MC generator name: %s\n",fMCGenerator.Data());
      Int_t nprod = header->NProduced();
			// printf("Number of produced particles: %i\n",nprod);

      // Pythia
			if (fMCGenerator == "Pythia") {
				fMCProcess = ((AliGenPythiaEventHeader*)header)->ProcessType();
				//printf("Pythia process type: %i\n",fMCProcess);
        switch(fMCProcess) {
				case 92:
				case 93:
				case 94:
				case 104: fMCProcessType = AliCEPBase::kProctypeSD; break;
				case 105: fMCProcessType = AliCEPBase::kProctypeDD; break;
				case 106: fMCProcessType = AliCEPBase::kProctypeCD; break;
				default:  fMCProcessType = AliCEPBase::kProctypeND; break;
				}
			}
			
      // DIME
			else if (fMCGenerator == "Dime") {
				fMCProcessType = AliCEPBase::kProctypeCD;
			}
			
      // DPMjet = Phojet
			else if (fMCGenerator == "DPMJET") {
				// see TDPMjet.h for definition of process codes
        
        fMCProcess = ((AliGenDPMjetEventHeader*)header)->ProcessType();
				//printf("DPMjet process type: %i\n",fMCProcess);
				switch(fMCProcess) {
				case 1:  fMCProcessType = AliCEPBase::kProctypeND; break;
				case 3:  fMCProcessType = AliCEPBase::kProctypeSD; break;
				case 4:  fMCProcessType = AliCEPBase::kProctypeDD; break;
				case 5:  fMCProcessType = AliCEPBase::kProctypeCD; break;
				default: fMCProcessType = AliCEPBase::kProctypeND; break;
				}
			}
      
		}
    
    // if (fMCProcessType == AliCEPBase::kProctypeCD)
    //   printf("Central Diffractive Event detected!\n");
    printf("MC process ID %i\n",fMCProcess);
	}
  
}

// ------------------------------------------------------------------------------
void AliCEPUtils::InitTrackCuts(Bool_t IsRun1, Int_t clusterCut)
{
	// TrackCuts
	fTrackCutListPrim->Clear();

	// Important message for 7TeV analysis (LHC10b,c,d,e)
	/*
	   Alexander Kalweit 
	   Email to PWG conveners on 22 Apr 2014
	   LHC10b&c (pass2):
	   ==================

     Default cut which is currently recommended:
     AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE/kFALSE, 1)
     Important is the second argument (=1) which replaces the cut on 70 clusters
     with a crossed rows cuts.
     !!! Please note, that a cut on 70 clusters is strongly discouraged in
     LHC10b&c pass2 data analysis !!!
     Changing to number of clusters (=0) and variations of the cut to 60 or 80
     should be included in the systematic studies

	   LHC10deh (pass2):
	   ==================
     Default cut which is currently recommended:
     AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE/kFALSE, 0) In this
     period, a cut on 70 clusters should be okay, however, changing to a crossed
     rows cut and lowering the cut to 60 clusters should be included in the
     systematic error study.
	   */

	// Will be used as standard trackcuts
  // distinguish between run1 and run 2 data
  
  Int_t    minclsTPC = 70;
  Int_t    minXrowsTPC = 70;
  Float_t  minNCrossedRowsTPC = 120; 
  Float_t  minRatioCrossedRowsOverFindableClustersTPC = 0.8; 
  Float_t  maxFractionSharedTPCCluster = 0.4;
  Double_t maxchi2perTPCcl = 4.;
  Double_t maxdcazITSTPC = 2.0;
  
  // Run2
  if (IsRun1) {

    // ITS+TPC
    AliESDtrackCuts *fcutITSTPC_P = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, clusterCut);
    
    fcutITSTPC_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
		if (clusterCut == 1) {
			fcutITSTPC_P->SetMinNCrossedRowsTPC(minXrowsTPC);
		} else {
			fcutITSTPC_P->SetMinNClustersTPC(minclsTPC);
		}
		fcutITSTPC_P->SetName("ITSTPC");
		AddTrackCut(fcutITSTPC_P);
		
    // ITS
    AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
		fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
		fcutITSSA_P->SetName("ITSSA");
		AddTrackCut(fcutITSSA_P);

  } else { //Run2
  
    // ITS+TPC
		AliESDtrackCuts *fcutITSTPC_P = new AliESDtrackCuts;
    //fcutITSTPC_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
		fcutITSTPC_P -> SetRequireTPCRefit(kTRUE);
		fcutITSTPC_P -> SetAcceptKinkDaughters(kFALSE);

		fcutITSTPC_P -> SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
		fcutITSTPC_P -> SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
		fcutITSTPC_P -> SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
		fcutITSTPC_P -> SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);
		fcutITSTPC_P -> SetMaxChi2TPCConstrainedGlobal(36);
    
		fcutITSTPC_P -> SetRequireITSRefit(kTRUE);
		fcutITSTPC_P -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		fcutITSTPC_P -> SetMaxChi2PerClusterITS(36);
    
		fcutITSTPC_P -> SetMaxDCAToVertexXYPtDep("(0.0182+0.0350/pt^1.01)");
		fcutITSTPC_P -> SetMaxDCAToVertexZ(maxdcazITSTPC);
    
		fcutITSTPC_P -> SetEtaRange(-2.0,2.0);
		fcutITSTPC_P -> SetPtRange(0.15);
    
		fcutITSTPC_P->SetName("ITSTPC");
		AddTrackCut(fcutITSTPC_P);
		
    // ITS
    AliESDtrackCuts *fcutITSSA_P = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kTRUE, 0);
		fcutITSSA_P->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
		fcutITSSA_P->SetName("ITSSA");
		AddTrackCut(fcutITSSA_P);
  
	}

}

// ------------------------------------------------------------------------------
