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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  A GTU track                                                           //
//                                                                        //
//  Author: J. Klein (Jochen.Klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

/* $Id: AliTRDtrackGTU.cxx 27566 2008-07-24 15:31:08Z cblume $ */

#include "TObject.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliTRDgtuParam.h"
#include "AliTRDtrackGTU.h"
#include "AliTRDtrackletGTU.h"
#include "AliESDTrdTrack.h"

ClassImp(AliTRDtrackGTU)
    
AliTRDtrackGTU::AliTRDtrackGTU() :
  TObject(),
  fStack(-1),
  fSector(-1),
  fPt(0),
  fPID(0),
  fTracklets(0x0),
  fTrackletMask(0),
  fNTracklets(0),
  fRefLayerIdx(-1),
  fZChannel(-1),
  fZSubChannel(-1),
  fA(0),
  fB(0),
  fC(0)
{
  fTracklets = new TClonesArray("AliTRDtrackletGTU", 6);
  for (Int_t iTracklet = 0; iTracklet < 6; iTracklet++)
      new ((*fTracklets)[iTracklet]) AliTRDtrackletGTU();
//  fTracklets->BypassStreamer(kFALSE);
}

AliTRDtrackGTU::~AliTRDtrackGTU()
{
  fTracklets->Delete();
  delete fTracklets;
}

void AliTRDtrackGTU::AddTracklet(AliTRDtrackletGTU *tracklet, Int_t layer) 
{
  if ( (fTrackletMask & (1 << layer)) != 0 ) {
    AliError(Form("Only one tracklet per layer (%i) possible! Mask: 0x%02x", layer, fTrackletMask));
    return;
  }

  new ((*fTracklets)[layer]) AliTRDtrackletGTU(*tracklet);
  fNTracklets++;
  fTrackletMask |= (1 << layer);
}

AliTRDtrackletGTU* AliTRDtrackGTU::GetTracklet(Int_t layer) 
{
  return ((AliTRDtrackletGTU*) (*fTracklets)[layer]);
}

Int_t AliTRDtrackGTU::GetNTracklets() const
{
  return fNTracklets;
}

Bool_t AliTRDtrackGTU::IsTrackletInLayer(Int_t layer) const 
{
  if ( (GetTrackletMask() & (1 << layer)) != 0)
    return kTRUE;
  else 
    return kFALSE;
}

void AliTRDtrackGTU::SetFitParams(Float_t a, Float_t b, Float_t c) 
{
  fA = a; 
  fB = b;
  fC = c;
}

Int_t AliTRDtrackGTU::GetZSubChannel() 
{
  if (fZSubChannel < 0) {
    for (Int_t layer = 0; layer < AliTRDgtuParam::GetNLayers(); layer++)
    {
      if (IsTrackletInLayer(layer))
	fZSubChannel = ((AliTRDtrackletGTU*) (*fTracklets)[layer])->GetSubChannel(GetZChannel());
    }
  }
  return fZSubChannel;
}

Int_t AliTRDtrackGTU::GetYapprox() 
{
  for (Int_t layer = 0; layer < AliTRDgtuParam::GetNLayers(); layer++) 
  {
    if (IsTrackletInLayer(layer))
      return ((AliTRDtrackletGTU*) (*fTracklets)[layer])->GetYProj();
  }
  return 0;
}

AliESDTrdTrack* AliTRDtrackGTU::CreateTrdTrack() const
{
    AliESDTrdTrack *trk = new AliESDTrdTrack();
    trk->SetPt(1./128. * fPt);
    trk->SetPID(fPID);
    trk->SetDetector(fSector * 30 + fStack * 6);
    return trk;
}
