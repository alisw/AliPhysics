// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSector2D.h"
#include "AliEveTPCSector3D.h"

#include <EveDet/AliEveTPCSectorData.h>

#include <TEveTrans.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TH1S.h>
#include <TH2S.h>
#include <TVirtualPad.h>

//==============================================================================
//==============================================================================
// AliEveTPCSector2D
//==============================================================================

//______________________________________________________________________________
//
// Displays TPC raw-data in 2D.
//
// fShowMax: true  - display maximum value for given time interval
//           false - display integral
// fAverage: only available when fShowMax = false; divide by time window width
//
// fUseTexture: use OpenGL textures to display data (fast rendering,
//   updates take the same time)
//

ClassImp(AliEveTPCSector2D)

/******************************************************************************/

AliEveTPCSector2D::AliEveTPCSector2D(const Text_t* n, const Text_t* t) :
  AliEveTPCSectorViz(n,t),

  fShowMax (kTRUE),
  fAverage (kFALSE),

  fUseTexture (kTRUE),
  fPickEmpty  (kFALSE),
  fPickMode   (1)
{
  // Constructor.
}

/******************************************************************************/

void AliEveTPCSector2D::MakeSector3D()
{
  // Make a 3D sector with same setting as this one.
  // It is added as a child ot this object.

  AliEveTPCSector3D* s = new AliEveTPCSector3D;
  s->SetDataSource(fTPCData);
  s->SetSectorID(fSectorID);
  s->SetAutoTrans(fAutoTrans);
  s->CopyVizParams(this);
  AddElement(s);
  ElementChanged(kFALSE, kTRUE);
}

/******************************************************************************/

void AliEveTPCSector2D::ComputeBBox()
{
  // Compute boundig-box.

  const AliEveTPCSectorData::SegmentInfo&  iSeg = AliEveTPCSectorData::GetInnSeg();
  const AliEveTPCSectorData::SegmentInfo& o2Seg = AliEveTPCSectorData::GetOut2Seg();

  BBoxInit();
  Float_t w = o2Seg.GetNMaxPads()*o2Seg.GetPadWidth()/2;
  fBBox[0] = -w;
  fBBox[1] =  w;
  fBBox[2] =  iSeg.GetRLow();
  fBBox[3] =  o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  fBBox[4] = -0.5; // Fake z-width to 1 cm.
  fBBox[5] =  0.5;
}

/******************************************************************************/

void AliEveTPCSector2D::PadSelected(Int_t row, Int_t pad)
{
  // Called when ctrl-mouse-left-click registered over a pad.

  // EVE -> Std convention
  Int_t sseg = fSectorID, srow = row;
  if (row >= AliEveTPCSectorData::GetInnSeg().GetNRows()) {
    sseg += 36;
    srow -= AliEveTPCSectorData::GetInnSeg().GetNRows();
  }
  switch (fPickMode)
    {
    case 0: {
      printf("AliEveTPCSector2DGL::ProcessSelection segment=%d, row=%d, pad=%d\n",
	     sseg, srow, pad);
      break;
    }
    case 1: {
      AliEveTPCSectorData* sectorData = GetSectorData();
      if (sectorData == 0) return;
      Int_t mint = fMinTime;
      Int_t maxt = fMaxTime;
      TH1S* h = new TH1S(Form("Seg%d_Row%d_TEvePad%d", sseg, srow, pad),
			 Form("Segment %d, Row %d, TEvePad %d", sseg, srow, pad),
			 maxt - mint +1 , mint, maxt);
      h->SetXTitle("Time");
      h->SetYTitle("ADC");
      AliEveTPCSectorData::PadIterator i = sectorData->MakePadIterator(row, pad);
      while (i.Next())
	h->Fill(i.Time(), i.Signal());
      h->Draw();
      gPad->Modified();
      gPad->Update();
      break;
    }
    case 2: {
      AliEveTPCSectorData* sectorData = GetSectorData();
      if (sectorData == 0) return;
      Int_t mint = fMinTime;
      Int_t maxt = fMaxTime;
      Int_t npad = AliEveTPCSectorData::GetNPadsInRow(row);
      TH2S* h = new TH2S(Form("Seg%d_Row%d", sseg, srow),
			 Form("Segment %d, Row %d", sseg, srow),
			 maxt - mint +1 , mint, maxt,
			 npad, 0, npad - 1);
      h->SetXTitle("Time");
      h->SetYTitle("TEvePad");
      h->SetZTitle("ADC");
      AliEveTPCSectorData::RowIterator i = sectorData->MakeRowIterator(row);
      while (i.NextPad())
	while (i.Next())
	  h->Fill(i.Time(), i.TEvePad(), i.Signal());
      h->Draw();
      gPad->Modified();
      gPad->Update();
      break;
    }
    } // switch
}

/******************************************************************************/

void AliEveTPCSector2D::Paint(Option_t* )
{
  // Paint object.

  if(fRnrSelf == kFALSE)
    return;

  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = GetMainColor();
  buffer.fTransparency = GetMainTransparency();
  if (HasMainTrans()) RefMainTrans().SetBuffer3D(buffer);
  buffer.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    return;
  }

  Error("AliEveTPCSector2D::Paint", "only direct OpenGL rendering supported.");
}
