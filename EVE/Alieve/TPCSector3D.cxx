// $Header$

#include "TPCSector3D.h"
#include <Alieve/TPCSectorData.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TStyle.h>
#include <TColor.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCSector3D
//

ClassImp(TPCSector3D)

TPCSector3D::TPCSector3D(const Text_t* n, const Text_t* t) :
  TPCSectorViz(n, t),

  fBoxSet       (n, t),
  fPointSetArray(n, t),
  fPointFrac    (0.25),
  fPointSize    (3),

  fPointSetOn     (0),
  fPointSetMaxVal (0),

  fDriftVel  (1),
  fZStep     (250.0/450)
{
  fRnrFrame = kFALSE;
  ComputeBBox();
}

TPCSector3D::~TPCSector3D()
{}

/**************************************************************************/

void TPCSector3D::SetRnrFrame(Bool_t rf)
{
  if(fRnrFrame != rf) {
    fRnrFrame = rf;
    IncRTS();
  }
}

/**************************************************************************/

void TPCSector3D::ComputeBBox()
{
  const TPCSectorData::SegmentInfo&  iSeg = TPCSectorData::GetInnSeg();
  const TPCSectorData::SegmentInfo& o2Seg = TPCSectorData::GetOut2Seg();

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  bbox_init();
#else
  BBoxInit();
#endif
  Float_t w = 0.5*o2Seg.GetNMaxPads()*o2Seg.GetPadWidth();
  fBBox[0] = -w;
  fBBox[1] =  w;
  fBBox[2] =  iSeg.GetRLow();
  fBBox[3] =  o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  fBBox[4] =  0;
  fBBox[5] =  TPCSectorData::GetZLength();
  Float_t* b = fBoxSet.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b[i] = fBBox[i]; }

}

void TPCSector3D::Paint(Option_t* /*option*/)
{
  if(fRnrElement == kFALSE)
    return;

  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = 1;
  buffer.fTransparency = 0;
  fHMTrans.SetBuffer3D(buffer);
  buffer.SetSectionsValid(TBuffer3D::kCore);
   
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    return;
  }

  Error("TPCSector3D::Paint", "only direct OpenGL rendering supported.");
  return;
}

/**************************************************************************/

void TPCSector3D::LoadPadrow(TPCSectorData::RowIterator& iter,
                             Float_t xs, Float_t ys, Float_t pw, Float_t ph) 
{
  Short_t pad, time, val;
  Float_t x0, x1, z0, z1;
  Float_t ym = ys + 0.5*ph;
  Float_t ye = ys + ph;
  Float_t zs = fZStep/fDriftVel;

  while (iter.NextPad()) {
    pad = iter.Pad();
    while (iter.Next()) {
      time = iter.Time();
      val  = iter.Signal();

      if(val <= fThreshold || time < fMinTime || time > fMaxTime)
	continue;

      if(fPointSetOn && val <= fPointSetMaxVal) {
	fPointSetArray.Fill(xs + (pad+0.5)*pw, ym, (time+0.5)*zs, val);
      } else {
	fBoxSet.fBoxes.push_back(Reve::Box());
	ColorFromArray(val, fBoxSet.fBoxes.back().color);
	x0 = xs + pad*pw;
	x1 = x0 + pw;
	z0 = time*zs;
	z1 = z0 + zs;
	Float_t* p = fBoxSet.fBoxes.back().vertices; 
	// front
	p[0] = x0;  p[1] = ys;  p[2] = z0;  p += 3;
	p[0] = x1;  p[1] = ys;  p[2] = z0;  p += 3;
	p[0] = x1;  p[1] = ye;  p[2] = z0;  p += 3;
	p[0] = x0;  p[1] = ye;  p[2] = z0;  p += 3;
	// back
	p[0] = x0;  p[1] = ys;  p[2] = z1;  p += 3;
	p[0] = x1;  p[1] = ys;  p[2] = z1;  p += 3;
	p[0] = x1;  p[1] = ye;  p[2] = z1;  p += 3;
	p[0] = x0;  p[1] = ye;  p[2] = z1;
      }
    }
  }
}

void TPCSector3D::UpdateBoxes()
{
  // Populate parent class Reve::BoxSet with digit information.

  // printf("TPCSector3D update boxes\n");

  fBoxSet.ClearSet();
  fPointSetArray.RemoveElements();

  TPCSectorData* data = GetSectorData();
  if (data != 0) {
    Bool_t isOn[3];
    isOn[0] = fRnrInn;
    isOn[1] = fRnrOut1;
    isOn[2] = fRnrOut2;

    SetupColorArray();
    SetupPointSetArray();

    // Loop over 3 main segments
    for (Int_t sId = 0; sId <= 2; ++sId) {
      if(isOn[sId] == kFALSE)
        continue;
      const TPCSectorData::SegmentInfo& sInfo = TPCSectorData::GetSeg(sId);
      Float_t sy = sInfo.GetRLow();
      for (Int_t row=sInfo.GetFirstRow(); row<=sInfo.GetLastRow(); ++row) {
        TPCSectorData::RowIterator i = data->MakeRowIterator(row);
        Float_t sx = -0.5*TPCSectorData::GetNPadsInRow(row)*sInfo.GetPadWidth();
        LoadPadrow(i, sx, sy, sInfo.GetPadWidth(), sInfo.GetPadHeight());
        sy += sInfo.GetPadHeight();
      }
    }

    if(fPointSetOn)
      fPointSetArray.CloseBins();
  }
}

void TPCSector3D::SetupPointSetArray()
{
  Int_t   nBins = (Int_t) TMath::Nint(fPointFrac*gStyle->GetNumberOfColors());
  if(nBins > 0) {
    fPointSetOn = kTRUE;
    fPointSetMaxVal = fThreshold + (Int_t) TMath::Nint(fPointFrac*(fMaxVal - fThreshold));
    // printf("SetupPointSetArray frac=%f nbins=%d psmv=%d (%d,%d)\n", fPointFrac, nBins, fPointSetMaxVal, fThreshold, fMaxVal);
    fPointSetArray.InitBins("", nBins, fThreshold, fPointSetMaxVal, kFALSE);
    for(Int_t b=0; b<nBins; ++b) {
      fPointSetArray.GetBin(b)->SetMarkerColor(gStyle->GetColorPalette(b));
    }
  } else {
    fPointSetOn = kFALSE;
  }
}
