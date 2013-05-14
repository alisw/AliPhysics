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

/* $Id: AliTRDgtuTMU.cxx 28397 2008-09-02 09:33:00Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Track Matching Unit (TMU) simulation                                  //
//                                                                        //
//  Author: J. Klein (Jochen.Klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TList.h"
#include "TVectorD.h"
#include "TMath.h"

#include "AliESDEvent.h"
#include "AliESDTrdTrack.h"

#include "AliLog.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

#include "AliTRDgtuParam.h"
#include "AliTRDgtuTMU.h"
#include "AliTRDtrackGTU.h"

ClassImp(AliTRDgtuTMU)

AliTRDgtuTMU::AliTRDgtuTMU(Int_t stack, Int_t sector) :
  TObject(),
  fTracklets(0x0),
  fTrackletsPostInput(0x0),
  fZChannelTracklets(0x0),
  fTracks(0x0),
  fGtuParam(0x0),
  fStack(-1),
  fSector(-1)
{
  // constructor which initializes the position information of the TMU

  fGtuParam = AliTRDgtuParam::Instance();

  // store tracklets per link
  fTracklets = new TObjArray*[fGtuParam->GetNLinks()];
  for (Int_t iLink = 0; iLink < fGtuParam->GetNLinks(); iLink++) {
    fTracklets[iLink] = new TObjArray();
  }

  // tracklets after input units per layer
  fTrackletsPostInput = new TObjArray*[fGtuParam->GetNLayers()];
  fZChannelTracklets = new TList*[fGtuParam->GetNLayers()];
  for (Int_t layer = 0;  layer <  fGtuParam->GetNLayers(); layer++) {
    fZChannelTracklets[layer] = new TList[fGtuParam->GetNZChannels()];
    fTrackletsPostInput[layer] = new TObjArray();
  }

  fTracks = new TList*[fGtuParam->GetNZChannels()];
  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
      fTracks[zch] = new TList[fGtuParam->GetNRefLayers()];
  }

  if (stack > -1)
      SetStack(stack);
  if (sector > -1)
      SetSector(sector);
}

AliTRDgtuTMU::~AliTRDgtuTMU()
{
  // destructor

  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
    delete [] fTracks[zch];
  }
  delete [] fTracks;
  for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
    delete [] fZChannelTracklets[layer];
    //    delete [] fTrackletsPostInput[layer];
  }
  delete [] fZChannelTracklets;
  //  delete [] fTrackletsPostInput;

  for (Int_t iLink = 0; iLink < fGtuParam->GetNLinks(); iLink++) {
    delete fTracklets[iLink];
  }
  delete [] fTracklets;
}

Bool_t AliTRDgtuTMU::SetSector(Int_t sector)
{
  // set the sector

  if (sector > -1 && sector < fGtuParam->GetGeo()->Nsector() ) {
    fSector = sector;
    return kTRUE;
  }

  AliError(Form("Invalid sector given: %i", sector));
  return kFALSE;
}

Bool_t AliTRDgtuTMU::SetStack(Int_t stack)
{
  // Set the stack (necessary for tracking)

  if (stack > -1 && stack < fGtuParam->GetGeo()->Nstack() ) {
    fStack = stack;
    return kTRUE;
  }

  AliError(Form("Invalid stack given: %i", stack));
  return kFALSE;
}

Bool_t AliTRDgtuTMU::Reset()
{
  // delete all tracks

  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
    for (Int_t reflayeridx = 0; reflayeridx < fGtuParam->GetNRefLayers(); reflayeridx++) {
      fTracks[zch][reflayeridx].Clear();
    }
  }

  // delete all added tracklets
  for (Int_t iLink = 0; iLink < fGtuParam->GetNLinks(); iLink++) {
    fTracklets[iLink]->Clear();
  }
  for (Int_t layer = 0; layer < fGtuParam->GetNLinks()/2; layer++) {
    fTrackletsPostInput[layer]->Clear();
    for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++)
      fZChannelTracklets[layer][zch].Clear();
  }

  fSector = -1;
  fStack = -1;

  return kTRUE;
}

Bool_t AliTRDgtuTMU::AddTracklet(AliTRDtrackletGTU *tracklet, Int_t link)
{
  // add a tracklet on the given link

  fTracklets[link]->Add(tracklet);
  return kTRUE;
}


Bool_t AliTRDgtuTMU::RunTMU(TList *ListOfTracks, AliESDEvent *esd)
{
  // performs the analysis as in a TMU module of the GTU, i. e.
  // track matching
  // calculation of track parameteres (pt, deflection, ???)

  if (fStack < 0 || fSector < 0) {
    AliError("No valid stack/sector set for this TMU! No tracking!");
    return kFALSE;
  }

  // ----- Input units -----
  AliDebug(1,"--------- Running Input units ----------");
  for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
    if (!RunInputUnit(layer)) {
      AliError(Form("Input unit in layer %i failed", layer));
      return kFALSE;
    }
  }

  // ----- Z-channel units -----
  AliDebug(1,"--------- Running Z-channel units ----------");
  for (Int_t layer = 0;  layer <  fGtuParam->GetNLayers(); layer++) {
    fZChannelTracklets[layer] = new TList[fGtuParam->GetNZChannels()];
    if (!RunZChannelUnit(layer)) {
      AliError(Form("Z-Channel unit in layer %i failed", layer));
      return kFALSE;
    }
  }

  // ----- track finding -----
  AliDebug(1,"--------- Running tracking units ----------");
  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
    if (!RunTrackFinder(zch, ListOfTracks)) {
      AliError(Form("Track Finder in z-channel %i failed", zch));
      return kFALSE;
    }
  }

  // ----- Track Merging -----
  if (!RunTrackMerging(ListOfTracks)) {
    AliError("Track merging failed");
    return kFALSE;
  }

  // ----- track reconstruction -----
  if (!RunTrackReconstruction(ListOfTracks)) {
    AliError("Track reconstruction failed");
    return kFALSE;
  }

  // ----- label calculation and ESD storage -----
  TIter next(ListOfTracks);
  while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) next()) {
    trk->CookLabel();
    if (esd) {
      AliESDTrdTrack *trdtrack = trk->CreateTrdTrack();
      esd->AddTrdTrack(trdtrack);
      delete trdtrack;
    }
  }

  return kTRUE;
}

Bool_t AliTRDgtuTMU::RunInputUnit(Int_t layer)
{
  // precalculations for the tracking and reconstruction

  Int_t iTrkl0 = 0; // A-side tracklet
  Int_t iTrkl1 = 0; // B-side tracklet

  while ((iTrkl0 < fTracklets[2*layer + 0]->GetEntriesFast()) || (iTrkl1 < fTracklets[2*layer + 1]->GetEntriesFast())) {
    if (iTrkl1 >= fTracklets[2*layer + 1]->GetEntriesFast()) {
      fTrackletsPostInput[layer]->AddLast(fTracklets[2*layer + 0]->At(iTrkl0));
      iTrkl0++;
    }
    else {
      if (iTrkl0 >= fTracklets[2*layer + 0]->GetEntriesFast()) {
	fTrackletsPostInput[layer]->AddLast(fTracklets[2*layer + 1]->At(iTrkl1));
	iTrkl1++;
      }
      else {
	if (((AliTRDtrackletGTU*) fTracklets[2*layer + 1]->At(iTrkl1))->GetZbin() <
	    ((AliTRDtrackletGTU*) fTracklets[2*layer + 0]->At(iTrkl0))->GetZbin()) {
	  fTrackletsPostInput[layer]->AddLast(fTracklets[2*layer + 1]->At(iTrkl1));
	  iTrkl1++;

	}
	else {
	  fTrackletsPostInput[layer]->AddLast(fTracklets[2*layer + 0]->At(iTrkl0));
	  iTrkl0++;
	}
      }
    }
  }

  TIter next(fTrackletsPostInput[layer]);

  while ( AliTRDtrackletGTU *trk = (AliTRDtrackletGTU*) next() ) {
    trk->SetIndex(fTrackletsPostInput[layer]->IndexOf(trk));

    Int_t alpha = (trk->GetYbin() >> fGtuParam->GetBitExcessY()) * fGtuParam->GetCiAlpha(layer);
    alpha = ( 2 * trk->GetdY() - (alpha >> fGtuParam->GetBitExcessAlpha()) + 1 ) >> 1;
    // wrapping projected alpha as in hardware
    if ((alpha < -64) || (alpha > 63))
      AliError(Form("alpha out of range: %i", alpha));
    alpha = alpha & 0x7f;
    if (alpha & 0x40)
      trk->SetAlpha(0xffffffc0 | alpha);
    else
      trk->SetAlpha(alpha);

    Int_t yproj = trk->GetdY() * (fGtuParam->GetCiYProj(layer));
    yproj = ( ( ( (yproj >> fGtuParam->GetBitExcessYProj()) + trk->GetYbin() ) >> 2) + 1) >> 1;
    trk->SetYProj(yproj);

    trk->SetYPrime(trk->GetYbin() + fGtuParam->GetYt(fStack, layer, trk->GetZbin()));

    AliDebug(10, Form("0x%08x: idx: %3i, z: %2i, y: %5i, dy: %3i, y': %5i, y_proj: %5i, alpha: %3i, pid: %3i, c: %5i, yt: %5i",
		      trk->GetTrackletWord(), trk->GetIndex(), trk->GetZbin(), trk->GetYbin(), trk->GetdY(), trk->GetYPrime(),
		      trk->GetYProj(), trk->GetAlpha(), trk->GetPID(),
		      fGtuParam->GetCiYProj(layer), fGtuParam->GetYt(fStack, layer, trk->GetZbin()) ));
  }
  return kTRUE;
}

Bool_t AliTRDgtuTMU::RunZChannelUnit(Int_t layer)
{
  // run the z-channel unit

  TIter next(fTrackletsPostInput[layer]);

  while (AliTRDtrackletGTU *trk = (AliTRDtrackletGTU*) next()) {
    for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
      if (fGtuParam->IsInZChannel(fStack, layer, zch, trk->GetZbin()) ) {
	trk->SetSubChannel(zch, fGtuParam->GetZSubchannel(fStack, layer, zch, trk->GetZbin()) );

	TIter nexttrkl(&fZChannelTracklets[layer][zch], kIterBackward);
	AliTRDtrackletGTU *t = 0x0;
	while ((t = (AliTRDtrackletGTU*) nexttrkl.Next())) {
	  if (t->GetSubChannel(zch) < trk->GetSubChannel(zch) ||
	      ((t->GetSubChannel(zch) == trk->GetSubChannel(zch)) && (t->GetYProj() < trk->GetYProj())) ) {
	    break;
	  }
	}
	if (t)
	  fZChannelTracklets[layer][zch].AddAfter(t, trk);
	else
	  fZChannelTracklets[layer][zch].AddFirst(trk);
      }
    }
    AliDebug(10, Form("stack %d, layer %2d: 0x%08x Z(2..0)=%i/%i/%i",
		      fStack, layer, trk->GetTrackletWord(),
		      fGtuParam->IsInZChannel(fStack, layer, 2, trk->GetZbin()) ? trk->GetSubChannel(2) : -1,
		      fGtuParam->IsInZChannel(fStack, layer, 1, trk->GetZbin()) ? trk->GetSubChannel(1) : -1,
		      fGtuParam->IsInZChannel(fStack, layer, 0, trk->GetZbin()) ? trk->GetSubChannel(0) : -1
		      ));
  }
  return kTRUE;
}

Bool_t AliTRDgtuTMU::RunTrackFinder(Int_t zch, TList* /* ListOfTracks */)
{
  // run the track finding

   Int_t 	*notr = new Int_t[fGtuParam->GetNLayers()];
   Int_t 	*ptrA = new Int_t[fGtuParam->GetNLayers()];
   Int_t 	*ptrB = new Int_t[fGtuParam->GetNLayers()];

   Bool_t 	*bHitA 	   = new Bool_t[fGtuParam->GetNLayers()];
   Bool_t 	*bHitB 	   = new Bool_t[fGtuParam->GetNLayers()];
   Bool_t 	*bAligned  = new Bool_t[fGtuParam->GetNLayers()];
   Bool_t 	*bAlignedA = new Bool_t[fGtuParam->GetNLayers()];
   Bool_t 	*bAlignedB = new Bool_t[fGtuParam->GetNLayers()];
   Bool_t 	*bDone 	   = new Bool_t[fGtuParam->GetNLayers()];
   Bool_t 	 ready;

   Int_t 	*inc 	  = new Int_t[fGtuParam->GetNLayers()];
   Int_t 	*incprime = new Int_t[fGtuParam->GetNLayers()];

// ----- signals within current layer -----
   Int_t 	yPlus;
   Int_t 	yMinus;
   Int_t 	ybPlus;
   Int_t 	ybMinus;
   Int_t 	alphaPlus;
   Int_t 	alphaMinus;
   Int_t 	nHits;
   Int_t 	nUnc;
   Int_t 	nWayBeyond;

   AliTRDtrackletGTU 	*trkRA 	= 0x0;	// reference tracklet A
   AliTRDtrackletGTU 	*trkRB 	= 0x0;	// reference tracklet B
   AliTRDtrackletGTU 	*trkA  	= 0x0;	// tracklet A in current layer
   AliTRDtrackletGTU 	*trkB  	= 0x0;	// tracklet B in current layer
/*
   AliTRDtrackletGTU 	*trkEnd = new AliTRDtrackletGTU();
   for (Int_t i = 0; i < fGtuParam->GetNZChannels(); i++)
       trkEnd->SetSubChannel(i, 7);
   trkEnd->SetYProj(0);
   trkEnd->SetAlpha(0);
*/

   for (Int_t refLayerIdx = 0; refLayerIdx < fGtuParam->GetNRefLayers(); refLayerIdx++) {
     Int_t reflayer = fGtuParam->GetRefLayer(refLayerIdx);
     AliDebug(5,Form("Tracking for z-channel: %i, reflayer: %i", zch, reflayer));

     ready = kFALSE; // ready if all channels done

//      for (Int_t iLayer = 0; iLayer < fGtuParam->GetNLayers(); iLayer++) {
//        for (Int_t iTrkl = 0; iTrkl < fZChannelTracklets[iLayer][zch].GetEntries(); iTrkl++) {
// 	 printf("%4i/%4i(%i,%i) ",
// 		((AliTRDtrackletGTU*) fZChannelTracklets[iLayer][zch].At(iTrkl))->GetYProj(),
// 		((AliTRDtrackletGTU*) fZChannelTracklets[iLayer][zch].At(iTrkl))->GetAlpha(),
// 		((AliTRDtrackletGTU*) fZChannelTracklets[iLayer][zch].At(iTrkl))->GetIndex(),
// 		((AliTRDtrackletGTU*) fZChannelTracklets[iLayer][zch].At(iTrkl))->GetSubChannel(zch)
// 		);
//        }
//        printf("\n");
//      }

     for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
       notr[layer] = fZChannelTracklets[layer][zch].GetEntries();
       ptrA[layer] = 0; // notr[layer] > 0 ? 0 : -1;
       ptrB[layer] = 1; // notr[layer] > 1 ? 1 : -1;

// not necessary here
//       bDone[layer] = (ptrB[layer] >= notr[layer] - 1); // trkB is last one
//       bDone[layer] = (notr[layer] <= 0);
//       ready = ready && bDone[layer];

       if (reflayer == 1)
	 AliDebug(5,Form("in layer: %i (zchannel = %i) there are: %i tracklets", layer, zch, notr[layer]));
     }

     if (ptrA[reflayer] < 0 && ptrB[reflayer] < 0)
       continue;

     while (!ready) {
       // ----- reference tracklets -----
       trkRA = 0x0;
       trkRB = 0x0;
       if (0 <= ptrA[reflayer] && ptrA[reflayer] < notr[reflayer])
	 trkRA = (AliTRDtrackletGTU*) fZChannelTracklets[reflayer][zch].At(ptrA[reflayer]);
       else  {
	 AliDebug(10,Form("No valid tracklet in the reference at ptr: %i! Nothing done!", ptrA[reflayer]));
	 break;
       }

       if (0 <= ptrB[reflayer] && ptrB[reflayer] < notr[reflayer])
	 trkRB = (AliTRDtrackletGTU*) fZChannelTracklets[reflayer][zch].At(ptrB[reflayer]);

       yPlus  	  = trkRA->GetYProj() + fGtuParam->GetDeltaY();
       yMinus 	  = trkRA->GetYProj() - fGtuParam->GetDeltaY();
       alphaPlus  = trkRA->GetAlpha() + fGtuParam->GetDeltaAlpha();
       alphaMinus = trkRA->GetAlpha() - fGtuParam->GetDeltaAlpha();
       if (trkRB) {
	 ybPlus 	  = trkRB->GetYProj() + fGtuParam->GetDeltaY();
	 ybMinus 	  = trkRB->GetYProj() - fGtuParam->GetDeltaY();
       }
       else { // irrelevant (should be, is it?)
	 ybPlus 	  = trkRA->GetYProj() + fGtuParam->GetDeltaY();
	 ybMinus 	  = trkRA->GetYProj() - fGtuParam->GetDeltaY();
       }

       nHits 	  = 0;
       nUnc 	  = 0;
       nWayBeyond = 0;

       for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
	 bHitA[layer] = bHitB[layer] = bAligned[layer] = kFALSE;
	 inc[layer] = incprime[layer] = 0;

	 if (layer == reflayer) {
	   bHitA[layer] = kTRUE;
	   bAligned[layer] = kTRUE;
	   nHits++;
	   continue;
	 }

	 trkA = 0x0;
	 trkB = 0x0;
	 if (0 <= ptrA[layer] && ptrA[layer] < notr[layer])
	   trkA = (AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrA[layer]);
	 if (0 <= ptrB[layer] && ptrB[layer] < notr[layer])
	   trkB = (AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrB[layer]);

	 bAlignedA[layer] = kFALSE;
	 bAlignedB[layer] = kFALSE;

	 if (trkA) {
	   bHitA[layer] = ( !(trkA->GetSubChannel(zch) < trkRA->GetSubChannel(zch) || (trkA->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkA->GetYProj() < yMinus) ) &&
			    !(trkA->GetSubChannel(zch) > trkRA->GetSubChannel(zch) || (trkA->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkA->GetYProj() > yPlus)  ) &&
			    !(trkA->GetAlpha() < alphaMinus) &&
			    !(trkA->GetAlpha() > alphaPlus) );
	   bAlignedA[layer] = !(trkA->GetSubChannel(zch) < trkRA->GetSubChannel(zch) || (trkA->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkA->GetYProj() < yMinus) );
	 }
	 else {
	   bHitA[layer] = 0;
	   bAlignedA[layer] = kTRUE;
	 }

	 if (trkB) {
	   bHitB[layer] = ( !(trkB->GetSubChannel(zch) < trkRA->GetSubChannel(zch) || (trkB->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkB->GetYProj() < yMinus) ) &&
			    !(trkB->GetSubChannel(zch) > trkRA->GetSubChannel(zch) || (trkB->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkB->GetYProj() > yPlus) ) &&
			    !(trkB->GetAlpha() < alphaMinus) &&
			    !(trkB->GetAlpha() > alphaPlus) );
	   bAlignedB[layer] = (trkB->GetSubChannel(zch) > trkRA->GetSubChannel(zch) || (trkB->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkB->GetYProj() > yPlus) );
	 }
	 else {
	   bHitB[layer] = 0;
	   bAlignedB[layer] = kTRUE;
	 }

	 bAligned[layer] = bAlignedA[layer] || bAlignedB[layer]; //???
//	 bAligned[layer] = bAlignedA[layer]; //???

	 if (bAligned[layer] && (bHitA[layer] || bHitB[layer]) )
	   nHits++;
	 else if (!bAligned[layer] )
	   nUnc++;
	 if (trkRB) {
	   if (trkA) {
	     if ((trkA->GetSubChannel(zch) > trkRB->GetSubChannel(zch)) || (trkA->GetSubChannel(zch) == trkRB->GetSubChannel(zch) && trkA->GetYProj() > ybPlus) )
	       nWayBeyond++;
	   }
	   else
	     nWayBeyond++;
	 }

	 //  pre-calculation for the layer shifting (alignment w. r. t. trkRB)
	 if (trkA) {
	     if(trkRB) {
		 if ((trkA->GetSubChannel(zch) < trkRB->GetSubChannel(zch)) || (trkA->GetSubChannel(zch) == trkRB->GetSubChannel(zch) && trkA->GetYProj() < ybMinus )) // could trkA be aligned for trkRB
		     incprime[layer] = 1;
		 else
		     incprime[layer] = 0;
	     }
	     else
		 incprime[layer] = 1;

	     if (trkB) {
		 if (trkRB) {
		     if ((trkB->GetSubChannel(zch) < trkRB->GetSubChannel(zch)) || (trkB->GetSubChannel(zch) == trkRB->GetSubChannel(zch) && trkB->GetYProj() < ybMinus )) // could trkB be aligned for trkRB
			 incprime[layer] = 2;
		 }
		 else
		     incprime[layer] = 2;
	     }
	 }
       } // end of loop over layers

       AliDebug(5,Form("logic calculation finished, Nhits: %i %s",
		       nHits, (nHits >= 4) ? "-> track found" : ""));

       if (nHits >= 4) {
	 // ----- track registration -----
	 AliTRDtrackGTU *track = new AliTRDtrackGTU();
	 track->SetSector(fSector);
	 track->SetStack(fStack);
	 for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++ ) {
	   if (bHitA[layer] || layer == reflayer)
	     track->AddTracklet((AliTRDtrackletGTU* ) fZChannelTracklets[layer][zch].At(ptrA[layer]), layer );
	   else if (bHitB[layer])
	     track->AddTracklet((AliTRDtrackletGTU* ) fZChannelTracklets[layer][zch].At(ptrB[layer]), layer );
	 }

	 Bool_t registerTrack = kTRUE;
	 for (Int_t layerIdx = refLayerIdx-1; layerIdx >= 0; layerIdx--) {
           if (track->IsTrackletInLayer(fGtuParam->GetRefLayer(layerIdx))) {
             if ((track->GetTracklet(fGtuParam->GetRefLayer(layerIdx)))->GetSubChannel(zch) > 0) {
               AliDebug(1,"Not registered");
		 registerTrack = kFALSE;
             }
           }
	 }
	 if (registerTrack) {
	     track->SetZChannel(zch);
	     track->SetRefLayerIdx(refLayerIdx);
	     fTracks[zch][refLayerIdx].Add(track);
	 }
       }

       if ( (nUnc != 0) && (nUnc + nHits >= 4) ) // could this position of the reference layer give some track //??? special check in case of hit?
	  inc[reflayer] = 0;
       else if (nWayBeyond > 2) // no track possible for both reference tracklets
	 inc[reflayer] = 2;
       else
	 inc[reflayer] = 1;

       if (inc[reflayer] != 0) // reflayer is shifted
	 for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
	   if (layer == reflayer)
	     continue;
	   inc[layer] = incprime[layer];
	 }
       else { // reflayer is not shifted
	 for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
	   if (layer == reflayer || notr[layer] == 0)
	     continue;
	   inc[layer] = 0;
	   trkA = 0x0;
	   trkB = 0x0;
	   if (0 <= ptrA[layer] && ptrA[layer] < notr[layer])
	     trkA = (AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrA[layer]);

	   if (0 <= ptrB[layer] && ptrB[layer] < notr[layer])
	     trkB = (AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrB[layer]);

	   if (trkA) {
	     if ( !(trkA->GetSubChannel(zch) < trkRA->GetSubChannel(zch) || (trkA->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkA->GetYProj() < yMinus) ) &&
		  !(trkA->GetSubChannel(zch) > trkRA->GetSubChannel(zch) || (trkA->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkA->GetYProj() > yPlus ) ) )  // trkA could hit trkRA
	     // if ( !(trkA->GetSubChannel(zch) < trkRA->GetSubChannel(zch) || (trkA->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkA->GetYProj() < yMinus) ) )
		   inc[layer] = 0;
	       else if (trkB) {
		   if ( trkB->GetSubChannel(zch) < trkRA->GetSubChannel(zch) || (trkB->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkB->GetYProj() < yMinus) )
		       inc[layer] = 2;
		   else if ( !(trkB->GetSubChannel(zch) > trkRA->GetSubChannel(zch) || (trkB->GetSubChannel(zch) == trkRA->GetSubChannel(zch) && trkB->GetYProj() > yPlus) ) )
		       inc[layer] = 1;
		   else
		       inc[layer] = incprime[layer];
	       }
	       else
		   inc[layer] = incprime[layer];
	   }
	 }
       }

       ready = kTRUE;
       for (Int_t layer = fGtuParam->GetNLayers()-1; layer >= 0; layer--) {

	 bDone[layer] = ptrB[layer] < 0 || ptrB[layer] >= notr[layer];
	 ready = ready && bDone[layer];

	 if (inc[layer] != 0 && ptrA[layer] >= notr[layer])
	   AliError(Form("Invalid increment: %i at ptrA: %i, notr: %i", inc[layer], ptrA[layer], notr[layer]));

//	 AliInfo(Form("Shifting layer: %i, notr: %i, ptrA: %i, ptrB: %i, inc: %i", layer, notr[layer], ptrA[layer], ptrB[layer], inc[layer]));
	 AliDebug(10,Form("    Layer: %i   %2i(%2i, %2i, %4i, %3i)%s%s   %2i(%2i, %2i, %4i, %3i)%s%s   +%i   %s  (no: %i)",
			  layer,
			  ptrA[layer],
			  (0 <= ptrA[layer] && ptrA[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrA[layer]))->GetIndex() : -1,
			  (0 <= ptrA[layer] && ptrA[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrA[layer]))->GetSubChannel(zch) : -1,
			  (0 <= ptrA[layer] && ptrA[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrA[layer]))->GetYProj() : -1,
			  (0 <= ptrA[layer] && ptrA[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrA[layer]))->GetAlpha() : -1,
			  bHitA[layer] ? "*" : " ", bAlignedA[layer] ? "+" : " ",
			  ptrB[layer],
			  (0 <= ptrB[layer] && ptrB[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrB[layer]))->GetIndex() : -1,
			  (0 <= ptrB[layer] && ptrB[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrB[layer]))->GetSubChannel(zch) : -1,
			  (0 <= ptrB[layer] && ptrB[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrB[layer]))->GetYProj() : -1,
			  (0 <= ptrB[layer] && ptrB[layer] < notr[layer]) ? ((AliTRDtrackletGTU*) fZChannelTracklets[layer][zch].At(ptrB[layer]))->GetAlpha() : -1,
			  bHitB[layer] ? "*" : " ", bAlignedB[layer] ? "+" : " ",
			  inc[layer], bDone[layer] ? "done" : "    ", notr[layer]));
 	 ptrA[layer] += inc[layer];
	 ptrB[layer] += inc[layer];
       }

     } // end of while
   } // end of loop over reflayer

   delete [] bHitA;
   delete [] bHitB;
   delete [] bAligned;
   delete [] bDone;
   delete [] inc;
   delete [] incprime;
   delete [] bAlignedA;
   delete [] bAlignedB;
   delete [] notr;
   delete [] ptrA;
   delete [] ptrB;

   return kTRUE;
}

Bool_t AliTRDgtuTMU::RunTrackMerging(TList* ListOfTracks)
{
    TList **tracksRefMerged = new TList*[fGtuParam->GetNZChannels()];
    TList **tracksRefUnique = new TList*[fGtuParam->GetNZChannels()];

    for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
	tracksRefMerged[zch] = new TList;
	tracksRefUnique[zch] = new TList;
    }

    TList *tracksZMergedStage0 = new TList;
    TList *tracksZUniqueStage0 = new TList;

    TList **tracksZSplitted = new TList*[2];
    for (Int_t i = 0; i < 2; i++)
	tracksZSplitted[i] = new TList;

    TList *tracksZMergedStage1 = new TList;

    AliTRDtrackGTU **trkInRefLayer = new AliTRDtrackGTU*[fGtuParam->GetNRefLayers()];

    // Bool_t done = kFALSE;
    Int_t minIdx = 0;
    AliTRDtrackGTU *trkStage0 = 0x0;

    for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
        // ----- Merging and Unification in Reflayers (seed_merger) -----
	do {
	  // done = kTRUE;
	    trkStage0 = 0x0;
	    for (Int_t refLayerIdx = 0; refLayerIdx < fGtuParam->GetNRefLayers(); refLayerIdx++) {
		trkInRefLayer[refLayerIdx] = (AliTRDtrackGTU*) fTracks[zch][refLayerIdx].First();
		if (trkInRefLayer[refLayerIdx] == 0) {
		    continue;
		}
		else if (trkStage0 == 0x0 ) {
		    trkStage0 = trkInRefLayer[refLayerIdx];
		    minIdx = refLayerIdx;
		    // done = kFALSE;
		}
		else if ( (trkInRefLayer[refLayerIdx]->GetZSubChannel() < trkStage0->GetZSubChannel()) ||
			  ((trkInRefLayer[refLayerIdx]->GetZSubChannel() == trkStage0->GetZSubChannel()) &&
			   ((trkInRefLayer[refLayerIdx]->GetYapprox() >> 3) < (trkStage0->GetYapprox() >> 3)) ) ) {
		    minIdx = refLayerIdx;
		    trkStage0 = trkInRefLayer[refLayerIdx];
		    // done = kFALSE;
		}
	    }
            if (!trkStage0)
              break;
	    tracksRefMerged[zch]->Add(trkStage0);
	    fTracks[zch][minIdx].RemoveFirst();
	} while (trkStage0 != 0);

	Uniquifier(tracksRefMerged[zch], tracksRefUnique[zch]);

	AliDebug(2, Form("zchannel %i", zch));
	TIter trackRefMerged(tracksRefMerged[zch]);
	while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) trackRefMerged())
	  AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, z_idx=%i",
			   AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
			   trk->GetTrackletIndex(5),
			   trk->GetTrackletIndex(4),
			   trk->GetTrackletIndex(3),
			   trk->GetTrackletIndex(2),
			   trk->GetTrackletIndex(1),
			   trk->GetTrackletIndex(0),
			   trk->GetYapprox() >> 3,
			   trk->GetZSubChannel()));
	AliDebug(2, "uniquified:");
	TIter trackRefUnique(tracksRefUnique[zch]);
	while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) trackRefUnique())
	  AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, z_idx=%i",
			   AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
			   trk->GetTrackletIndex(5),
			   trk->GetTrackletIndex(4),
			   trk->GetTrackletIndex(3),
			   trk->GetTrackletIndex(2),
			   trk->GetTrackletIndex(1),
			   trk->GetTrackletIndex(0),
			   trk->GetYapprox() >> 3,
			   trk->GetZSubChannel()));
    }

// ----- Merging in zchannels - 1st stage -----

    if (AliTRDgtuParam::GetUseGTUmerge()) {
      Int_t notEmpty;
      do {
	Bool_t lowerThan[3] = { kFALSE, kFALSE, kFALSE };
	AliTRDtrackGTU *trk[3] = { 0x0, 0x0, 0x0 };
	for (Int_t iChannel = 0; iChannel < fGtuParam->GetNZChannels(); ++iChannel)
	  trk[iChannel] = (AliTRDtrackGTU*) tracksRefUnique[iChannel]->First();

	for (Int_t iChannel = 0; iChannel < fGtuParam->GetNZChannels(); ++iChannel) {
	  AliTRDtrackGTU *trk1 = trk[iChannel];
	  AliTRDtrackGTU *trk2 = trk[(iChannel + 1) % 3];
	  if (trk1 && trk2) {
	    Int_t sortnum1 = (trk1->GetZChannel() + 3 * trk1->GetZSubChannel()) / 2 - 1;
	    Int_t sortnum2 = (trk2->GetZChannel() + 3 * trk2->GetZSubChannel()) / 2 - 1;
	    AliDebug(5, Form("comparing tracks %i - %i: %i - %i",
			     trk1->GetZChannel(), trk2->GetZChannel(),
			     sortnum1, sortnum2));
	    if ( (sortnum1 < sortnum2) ||
		 ((sortnum1 == sortnum2) &&
		  ((trk1->GetYapprox() >> 3) < (trk2->GetYapprox() >> 3)) ) ) {
	      lowerThan[iChannel] = kTRUE;
	    }
	  }
	}

	notEmpty = (trk[2] ? (1 << 2) : 0) |
	  (trk[1] ? (1 << 1) : 0) |
	  (trk[0] ? (1 << 0) : 0);
	Int_t pop = -1;

	switch (notEmpty) {
	  // one track only
	case 0x1:
	  pop = 0;
	  break;
	case 0x2:
	  pop = 1;
	  break;
	case 0x4:
	  pop = 2;
	  break;

	  // two tracks
	case 0x3:
	  if (lowerThan[0])
	    pop = 0;
	  else
	    pop = 1;
	  break;
	case 0x5:
	  if (lowerThan[2])
	    pop = 2;
	  else
	    pop = 0;
	  break;
	case 0x6:
	  if (lowerThan[1])
	    pop = 1;
	  else
	    pop = 2;
	  break;

	  // three tracks
	case 0x7:
	  if (lowerThan[0]) {
	    if (lowerThan[2])
	      pop = 2;
	    else
	      pop = 0;
	  } else {
	    if (lowerThan[1])
	      pop = 1;
	    else
	      pop = 2;
	  }
	  break;

	  // no tracks
	default:
	  // nop
	  ;
	}

	if (pop > -1) {
	  tracksZMergedStage0->Add(trk[pop]);
	  tracksRefUnique[pop]->RemoveFirst();
	}
      } while (notEmpty);
    }
    else {
      // there is still a difference between this implementation and
      // the hardware algorithm, only for expert use

      do {
	// done = kTRUE;
	trkStage0 = 0x0;
	// compare tracks from all adjacent zchannels
	// (including comparison of channels 2 and 0)
        for (Int_t i = fGtuParam->GetNZChannels() - 1; i > -1; i--) {
	  Int_t zch = i % 3;
	  AliTRDtrackGTU *trk = (AliTRDtrackGTU*) tracksRefUnique[zch]->First();
	  if (trk == 0) {
	    continue;
	  }
	  else if (trkStage0 == 0x0 ) {
	    trkStage0 = trk;
	    minIdx = zch;
	    // done = kFALSE;
	  }
	  else {
	    Int_t sortnum1 = (trk->GetZChannel() + 3 * trk->GetZSubChannel()) / 2 - 1;
	    Int_t sortnum2 = (trkStage0->GetZChannel() + 3 * trkStage0->GetZSubChannel()) / 2 - 1;
	    AliDebug(5, Form("comparing tracks %i - %i: %i - %i",
			     trk->GetZChannel(), trkStage0->GetZChannel(),
			     sortnum1, sortnum2));
	    if ( (sortnum1 < sortnum2) ||
		 ((sortnum1 == sortnum2) &&
		  ((trk->GetYapprox() >> 3) < (trkStage0->GetYapprox() >> 3)) ) ) {
	      minIdx = zch;
	      trkStage0 = trk;
	      // done = kFALSE;
	    }
	  }
	}

	if (!trkStage0)
          break;
	tracksZMergedStage0->Add(trkStage0);
	tracksRefUnique[minIdx]->RemoveFirst();
      } while (trkStage0 != 0);
    }

    Uniquifier(tracksZMergedStage0, tracksZUniqueStage0);

    AliDebug(2, "stage 0:");
    TIter trackZMergedStage0(tracksZMergedStage0);
    while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) trackZMergedStage0())
      AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, zch=%i idx=%i",
		       AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
		       trk->GetTrackletIndex(5),
		       trk->GetTrackletIndex(4),
		       trk->GetTrackletIndex(3),
		       trk->GetTrackletIndex(2),
		       trk->GetTrackletIndex(1),
		       trk->GetTrackletIndex(0),
		       trk->GetYapprox() >> 3,
		       trk->GetZChannel(),
		       trk->GetZSubChannel()));

    AliDebug(2, "uniquified:");
    TIter trackZUniqueStage0(tracksZUniqueStage0);
    while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) trackZUniqueStage0())
      AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, zch=%i idx=%i",
		       AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
		       trk->GetTrackletIndex(5),
		       trk->GetTrackletIndex(4),
		       trk->GetTrackletIndex(3),
		       trk->GetTrackletIndex(2),
		       trk->GetTrackletIndex(1),
		       trk->GetTrackletIndex(0),
		       trk->GetYapprox() >> 3,
		       trk->GetZChannel(),
		       trk->GetZSubChannel()));

// ----- Splitting in z -----

    TIter next(tracksZUniqueStage0);
    while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) next()) {
	tracksZSplitted[(trk->GetZChannel() + 3 * (trk->GetZSubChannel() - 1)) % 2]->Add(trk);
    }

    for (Int_t i = 0; i < 2; i++) {
      AliDebug(2, Form("split %i:", i));
      TIter trackZSplit(tracksZSplitted[i]);
      while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) trackZSplit())
	AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, zch=%i idx=%i",
			 AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
			 trk->GetTrackletIndex(5),
			 trk->GetTrackletIndex(4),
			 trk->GetTrackletIndex(3),
			 trk->GetTrackletIndex(2),
			 trk->GetTrackletIndex(1),
			 trk->GetTrackletIndex(0),
			 trk->GetYapprox() >> 3,
			 trk->GetZChannel(),
			 trk->GetZSubChannel()));
    }

// ----- Merging in zchanels - 2nd stage -----

    do {
      // done = kTRUE;
	trkStage0 = 0x0;
	for (Int_t i = 1; i >= 0; i--) {
	    AliTRDtrackGTU *trk = (AliTRDtrackGTU*) tracksZSplitted[i]->First();
	    if (trk == 0) {
		continue;
	    }
	    else if (trkStage0 == 0x0 ) {
		trkStage0 = trk;
		minIdx = i;
		// done = kFALSE;
	    }
	    else if (  (((trk->GetZChannel() + 3 * (trk->GetZSubChannel() - 1)) / 2) <  ((trkStage0->GetZChannel() + 3 * (trkStage0->GetZSubChannel() - 1)) / 2)) ||
		      ((((trk->GetZChannel() + 3 * (trk->GetZSubChannel() - 1)) / 2) == ((trkStage0->GetZChannel() + 3 * (trkStage0->GetZSubChannel() - 1)) / 2)) &&
		       ((trk->GetYapprox() >> 3) < (trkStage0->GetYapprox() >> 3)) ) ) {
		minIdx = i;
		trkStage0 = trk;
		// done = kFALSE;
	    }
	}

        if (!trkStage0)
          break;
	tracksZMergedStage1->Add(trkStage0);
	tracksZSplitted[minIdx]->RemoveFirst();
    } while (trkStage0 != 0);

    Uniquifier(tracksZMergedStage1, ListOfTracks);

    AliDebug(2, "merged:");
    TIter trackZMergedStage1(tracksZMergedStage1);
    while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) trackZMergedStage1())
      AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, zch=%i idx=%i",
		       AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
		       trk->GetTrackletIndex(5),
		       trk->GetTrackletIndex(4),
		       trk->GetTrackletIndex(3),
		       trk->GetTrackletIndex(2),
		       trk->GetTrackletIndex(1),
		       trk->GetTrackletIndex(0),
		       trk->GetYapprox() >> 3,
		       trk->GetZChannel(),
		       trk->GetZSubChannel()));

    AliDebug(2, "uniquified:");
    TIter track(ListOfTracks);
    while (AliTRDtrackGTU *trk = (AliTRDtrackGTU*) track())
      AliDebug(2, Form("track ref layer %i : %i %i %i %i %i %i, y=%i, zch=%i idx=%i",
		       AliTRDgtuParam::GetRefLayer(trk->GetRefLayerIdx()),
		       trk->GetTrackletIndex(5),
		       trk->GetTrackletIndex(4),
		       trk->GetTrackletIndex(3),
		       trk->GetTrackletIndex(2),
		       trk->GetTrackletIndex(1),
		       trk->GetTrackletIndex(0),
		       trk->GetYapprox() >> 3,
		       trk->GetZChannel(),
		       trk->GetZSubChannel()));

    // cleaning up
    for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) {
      delete tracksRefMerged[zch];
      delete tracksRefUnique[zch];
    }
    delete [] tracksRefMerged;
    delete [] tracksRefUnique;


    delete tracksZMergedStage0;
    delete tracksZUniqueStage0;

    for (Int_t i = 0; i < 2; i++)
      delete tracksZSplitted[i];
    delete [] tracksZSplitted;

    delete tracksZMergedStage1;

    delete [] trkInRefLayer;

    return kTRUE;
}

Bool_t AliTRDgtuTMU::RunTrackReconstruction(TList* ListOfTracks)
{
  // run the track reconstruction for all tracks in the list

  TIter next(ListOfTracks);
  while (AliTRDtrackGTU *track = (AliTRDtrackGTU*) next()) {
    CalculateTrackParams(track);
    CalculatePID(track);
    AliDebug(1, Form("track with pid = %i", track->GetPID()));
  }
  return kTRUE;
}

Bool_t AliTRDgtuTMU::CalculatePID(AliTRDtrackGTU *track)
{
  // calculate PID for the given track
  if (!track) {
    AliError("No track to calculate!");
    return kFALSE;
  }

  if (AliTRDgtuParam::GetUseGTUconst()) {
    // averaging as in GTU
    AliDebug(1, "using GTU constants for PID calculation");
    ULong64_t coeff;

    // selection of coefficient for averaging
    switch(track->GetTrackletMask()){
    case 0x3f:
      // 6 tracklets
      coeff=0x5558; // approx. 1/6
      break;

    case 0x3e:
    case 0x3d:
    case 0x3b:
    case 0x37:
    case 0x2f:
    case 0x1f:
      // 5 tracklets
      coeff=0x6666; // approx. 1/5
      break;

    default:
      // 4 tracklets
      coeff=0x8000; // = 0.25
    }
    coeff &= 0x1ffff; // 17-bit constant

    ULong64_t sum = 0;
    Int_t i = 0;
    for (Int_t iLayer = 0; iLayer < fGtuParam->GetNLayers(); iLayer++) {
      if ((track->GetTrackletMask() >> iLayer) & 1) {
	sum += track->GetTracklet(iLayer)->GetPID();
	++i;
      }
    }

    Float_t av = 1./i * sum;
    sum = sum & 0x7ff;
    ULong64_t prod   = (sum * coeff) & 0xfffffffffull;
    ULong64_t prodFinal = ((prod >> 17) + ((prod >> 16) & 1)) & 0xff;

    if (TMath::Abs((prodFinal & 0xff) - av) > 0.5)
      AliError(Form("problem with PID averaging (hw <-> ar): %3lli <-> %4.1f", prodFinal & 0xff, av));
    track->SetPID(prodFinal & 0xff);

    return kTRUE;
  }
  else {

    // simple averaging
    Int_t nTracklets = 0;
    Int_t pidSum = 0;
    for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
      if (!track->IsTrackletInLayer(layer)) {
	continue;
      }
      AliTRDtrackletGTU *trk = track->GetTracklet(layer);
      pidSum += trk->GetPID();
      nTracklets++;
    }

    if (nTracklets>0)
      track->SetPID(pidSum/nTracklets);
    else
      AliError("Track without contributing tracklets, no PID assigned");

    return kTRUE;
  }
}

Bool_t AliTRDgtuTMU::CalculateTrackParams(AliTRDtrackGTU *track)
{
  // calculate the track parameters

  if (!track) {
    AliError("No track to calculate!");
    return kFALSE;
  }

  Int_t a = 0;
  Float_t b = 0;
  Float_t c = 0;
  Float_t x1;
  Float_t x2;

  AliDebug(5,Form("There are %i tracklets in this track.", track->GetNTracklets()));

  for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
    if (!track->IsTrackletInLayer(layer)) {
      continue;
    }
    AliTRDtrackletGTU *trk = track->GetTracklet(layer);
    if (!trk) {
      AliError(Form("Could not get tracklet in layer %i\n", layer));
      continue;
    }
    AliDebug(10,Form("  layer %i trk yprime: %6i, aki: %6i", layer, trk->GetYPrime(),
		     fGtuParam->GetAki(track->GetTrackletMask(), layer)));
    a += (((fGtuParam->GetAki(track->GetTrackletMask(), layer) * trk->GetYPrime()) >> 7) + 1) >> 1;
    b += fGtuParam->GetBki(track->GetTrackletMask(), layer) * trk->GetYPrime() * fGtuParam->GetBinWidthY();
    c += fGtuParam->GetCki(track->GetTrackletMask(), layer) * trk->GetYPrime() * fGtuParam->GetBinWidthY();
  }
  if (a < 0)
      a += 3;
  a = a >> 2;

  track->SetFitParams(a << 2, b, c);

  fGtuParam->GetIntersectionPoints(track->GetTrackletMask(), x1, x2);

  AliDebug(5,Form("  Track parameters: a16 = %i, a18 = %i, b = %f, c = %f, x1 = %f, x2 = %f, pt = %f (trkl mask: %i)",
		  a, a << 2, b, c, x1, x2, track->GetPt(), track->GetTrackletMask()));

  return kTRUE;
}


Bool_t AliTRDgtuTMU::Uniquifier(const TList *inlist, TList *outlist)
{
  // remove multiple occurences of the same track

    TIter next(inlist);
    AliTRDtrackGTU *trkStage0 = 0x0;
    AliTRDtrackGTU *trkStage1 = 0x0;

    do {
	trkStage0 = (AliTRDtrackGTU*) next();

	Bool_t tracksEqual = kFALSE;
	if (trkStage0 != 0 && trkStage1 != 0) {
	    for (Int_t layer = 0; layer < fGtuParam->GetNLayers(); layer++) {
		if (trkStage0->GetTrackletIndex(layer) != -1 && trkStage0->GetTrackletIndex(layer) == trkStage1->GetTrackletIndex(layer)) {
		    tracksEqual = kTRUE;
		    break;
		}
	    }
	}

	if (tracksEqual) {
	    if (trkStage0->GetNTracklets() >= trkStage1->GetNTracklets())
		trkStage1 = trkStage0;
	}
	else {
	    if (trkStage1 != 0x0)
		outlist->Add(trkStage1);
	    trkStage1 = trkStage0;
	}

    } while (trkStage1 != 0x0);
    return kTRUE;
}
