// $Header$

#include "TPCSector3D.h"
#include <Alieve/TPCSectorData.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCSector3D
//

ClassImp(TPCSector3D)

TPCSector3D::TPCSector3D(const Text_t* n, const Text_t* t) :
  TPCSectorViz(n, t),

  fBoxSet    (n, t),
  fDriftVel  (1),
  fZStep     (250.0/450)
{
  fRnrFrame = false;
  MakeFrameBox();
  ComputeBBox();
}

TPCSector3D::~TPCSector3D()
{}

/**************************************************************************/

UInt_t TPCSector3D::IncRTS()
{
  UpdateBoxes();
  return ++fRTS;
}

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
  Float_t w = o2Seg.GetNMaxPads()*o2Seg.GetPadWidth()/2;
  fBBox[0] = -w;
  fBBox[1] =  w;
  fBBox[2] =  iSeg.GetRLow();
  fBBox[3] =  o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  fBBox[4] = -0.5;
  fBBox[5] =  250.5;
  Float_t* b = fBoxSet.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b[i] = fBBox[i]; }

}

void TPCSector3D::Paint(Option_t* option)
{
  if(fRnrElement) {
    fBoxSet.SetTrans(fTrans);
    memcpy(fBoxSet.ArrTrans(), fMatrix, 16*sizeof(Double_t));
    fBoxSet.Paint(option);
  }
}

/**************************************************************************/

void TPCSector3D::MakeFrameBox()
{
  const TPCSectorData::SegmentInfo&  iSeg = TPCSectorData::GetInnSeg();
  const TPCSectorData::SegmentInfo& o2Seg = TPCSectorData::GetOut2Seg();

  Float_t x1, x2, y1, y2, D;
  x1 = 0.5*TPCSectorData::GetNPadsInRow(0)*iSeg.GetPadWidth();
  x2 = 0.5*o2Seg.GetNMaxPads()*o2Seg.GetPadWidth();
  y1 = iSeg.GetRLow();
  y2 = o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  D  = 0.5*fZStep;

  ColorFromIdx(fFrameColor, fFrameBox.color);
  Float_t* p = fFrameBox.vertices;
  // bottom
  p[0] = -x1;  p[1] = y1;  p[2] = -D;  p += 3;
  p[0] =  x1;  p[1] = y1;  p[2] = -D;  p += 3;
  p[0] =  x2;  p[1] = y2;  p[2] = -D;  p += 3;
  p[0] = -x2;  p[1] = y2;  p[2] = -D;  p += 3;
  // top
  p[0] = -x1;  p[1] = y1;  p[2] =  D;  p += 3;
  p[0] =  x1;  p[1] = y1;  p[2] =  D;  p += 3;
  p[0] =  x2;  p[1] = y2;  p[2] =  D;  p += 3;
  p[0] = -x2;  p[1] = y2;  p[2] =  D;  p += 3;
}

/**************************************************************************/

void TPCSector3D::LoadPadrow(TPCSectorData::RowIterator& iter,
                             Float_t xs, Float_t ys, Float_t pw, Float_t ph) 
{
  Short_t pad, time, val;
  Float_t x0, x1, z0, z1;
  Float_t ye = ys + ph;

  while (iter.NextPad()) {
    pad = iter.Pad();
    while (iter.Next()) {
      time = iter.Time();
      val  = iter.Signal();

      if(val <= fThreshold || time < fMinTime || time > fMaxTime)
	continue;

      fBoxSet.fBoxes.push_back(Reve::Box());
      SetupColor(val, fBoxSet.fBoxes.back().color);
      x0 = xs + pad*pw;
      x1 = x0 + pw;
      z0 = fZStep*time/fDriftVel;
      z1 = z0 + fZStep/fDriftVel;
      Float_t* p = fBoxSet.fBoxes.back().vertices; 
      // front
      p[0] = x0;  p[1] = ys;  p[2] = z0;  p += 3;
      p[0] = x1;  p[1] = ys;  p[2] = z0;  p += 3;
      p[0] = x1;  p[1] = ye;  p[2] = z0;  p += 3;
      p[0] = x0;  p[1] = ye;  p[2] = z0;  p += 3;
      // back
      z0 += fZStep;
      p[0] = x0;  p[1] = ys;  p[2] = z0;  p += 3;
      p[0] = x1;  p[1] = ys;  p[2] = z0;  p += 3;
      p[0] = x1;  p[1] = ye;  p[2] = z0;  p += 3;
      p[0] = x0;  p[1] = ye;  p[2] = z0;
    }
  }
}

void TPCSector3D::UpdateBoxes()
{
  // Populate parent class Reve::BoxSet with digit information.

  // printf("TPCSector3D update boxes\n");

  fBoxSet.ClearSet();

  TPCSectorData* data = GetSectorData();
  if (data != 0) {
    Bool_t isOn[3];
    isOn[0] = fRnrInn;
    isOn[1] = fRnrOut1;
    isOn[2] = fRnrOut2;

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
  }

  if(fRnrFrame)
    fBoxSet.fBoxes.push_back(fFrameBox);
}
