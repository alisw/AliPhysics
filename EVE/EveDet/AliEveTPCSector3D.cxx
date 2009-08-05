// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSector3D.h"
#include <EveDet/AliEveTPCSectorData.h>

#include <TEveTrans.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TStyle.h>

//==============================================================================
//==============================================================================
// AliEveTPCSector3D
//==============================================================================

//______________________________________________________________________________
//
// Visualization of TPC raw-data in 3D.

ClassImp(AliEveTPCSector3D)

AliEveTPCSector3D::AliEveTPCSector3D(const Text_t* n, const Text_t* t) :
  AliEveTPCSectorViz(n, t),

  fBoxSet       (n, t),
  fPointSetArray(n, t),
  fPointFrac    (0.25),
  fPointSize    (3),

  fPointSetOn     (0),
  fPointSetMaxVal (0),

  fDriftVel  (1.07),
  fZStep     (250.0/900)
{
  // Constructor.

  fRnrFrame = kFALSE;
  ComputeBBox();
}

/******************************************************************************/

void AliEveTPCSector3D::SetRnrFrame(Bool_t rf)
{
  // Setter for fRnrFrame.

  if(fRnrFrame != rf) {
    fRnrFrame = rf;
    IncRTS();
  }
}

/******************************************************************************/

void AliEveTPCSector3D::ComputeBBox()
{
  // Compute bounding box containing whole sector.

  const AliEveTPCSectorData::SegmentInfo&  iSeg = AliEveTPCSectorData::GetInnSeg();
  const AliEveTPCSectorData::SegmentInfo& o2Seg = AliEveTPCSectorData::GetOut2Seg();

  BBoxInit();

  Float_t w = 0.5*o2Seg.GetNMaxPads()*o2Seg.GetPadWidth();
  fBBox[0] = -w;
  fBBox[1] =  w;
  fBBox[2] =  iSeg.GetRLow();
  fBBox[3] =  o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  fBBox[4] =  0;
  fBBox[5] =  AliEveTPCSectorData::GetZLength();
  Float_t* b = fBoxSet.AssertBBox();
  for (Int_t i=0; i<6; ++i) { b[i] = fBBox[i]; }
}

void AliEveTPCSector3D::Paint(Option_t* /*option*/)
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

  Error("AliEveTPCSector3D::Paint", "only direct OpenGL rendering supported.");
  return;
}

/******************************************************************************/

void AliEveTPCSector3D::LoadPadrow(AliEveTPCSectorData::RowIterator& iter,
                                   Float_t xs, Float_t ys,
                                   Float_t pw, Float_t ph)
{
  // Load data of one padrow. Fill internal boxset and pointset objects.

  Short_t pad, time, val;
  Float_t x0, z0;
  Float_t ym = ys + 0.5*ph;
  Float_t zs = fZStep/fDriftVel;

  while (iter.NextPad())
  {
    pad = iter.TEvePad();
    while (iter.Next())
    {
      time = iter.Time();
      val  = iter.Signal();

      if (val <= fThreshold || time < fMinTime || time > fMaxTime)
	continue;

      if (fPointSetOn && val <= fPointSetMaxVal)
      {
	fPointSetArray.Fill(xs + (pad+0.5)*pw, ym, (time+0.5)*zs, val);
      }
      else
      {
	x0 = xs + pad*pw;
	z0 = time*zs;
	fBoxSet.AddBox(x0, ys, z0, pw, ph, zs);
        fBoxSet.DigitColor(ColorFromArray(val));
      }
    }
  }
}

void AliEveTPCSector3D::UpdateBoxesAndPoints()
{
  // Populate BoxSet and PointSet with digit information.

  // printf("AliEveTPCSector3D update boxes\n");

  fBoxSet.Reset(TEveBoxSet::kBT_AABox, kTRUE, 16384);
  // Brutally delete sub-pointsets so that destruction via TEveManager is
  // avoided. This only works because fPointSetArray is never published.
  for (Int_t i = 0; i < fPointSetArray.GetNBins(); ++i)
    delete fPointSetArray.GetBin(i);
  fPointSetArray.RemoveElementsLocal();

  AliEveTPCSectorData* data = GetSectorData();
  if (data != 0) {
    Bool_t isOn[3];
    isOn[0] = fRnrInn;
    isOn[1] = fRnrOut1;
    isOn[2] = fRnrOut2;

    SetupColorArray();
    SetupPointSetArray();

    // Loop over 3 main segments
    for (Int_t sId = 0; sId <= 2; ++sId)
    {
      if (isOn[sId] == kFALSE)
        continue;
      const AliEveTPCSectorData::SegmentInfo& sInfo = AliEveTPCSectorData::GetSeg(sId);
      Float_t sy = sInfo.GetRLow();
      for (Int_t row=sInfo.GetFirstRow(); row<=sInfo.GetLastRow(); ++row)
      {
        AliEveTPCSectorData::RowIterator i = data->MakeRowIterator(row);
        Float_t sx = -0.5*AliEveTPCSectorData::GetNPadsInRow(row)*sInfo.GetPadWidth();
        LoadPadrow(i, sx, sy, sInfo.GetPadWidth(), sInfo.GetPadHeight());
        sy += sInfo.GetPadHeight();
      }
    }

    fBoxSet.RefitPlex();
    if (fPointSetOn)
      fPointSetArray.CloseBins();
  }
}

void AliEveTPCSector3D::SetupPointSetArray()
{
  // Setup fPointSetArray for current settings.

  Int_t nBins = (Int_t) TMath::Nint(fPointFrac*gStyle->GetNumberOfColors());
  if (nBins > 0) {
    fPointSetOn = kTRUE;
    fPointSetMaxVal = fThreshold + (Int_t) TMath::Nint(fPointFrac*(fMaxVal - fThreshold));
    // printf("SetupPointSetArray frac=%f nbins=%d psmv=%d (%d,%d)\n", fPointFrac, nBins, fPointSetMaxVal, fThreshold, fMaxVal);
    fPointSetArray.InitBins("", nBins, fThreshold, fPointSetMaxVal);
    for (Int_t b=0; b<nBins; ++b) {
      fPointSetArray.GetBin(b)->SetMarkerColor(gStyle->GetColorPalette(b));
    }
  } else {
    fPointSetOn = kFALSE;
  }
}
