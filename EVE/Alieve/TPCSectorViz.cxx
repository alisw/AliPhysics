// $Header$

#include "TPCSectorViz.h"

#include <Alieve/TPCData.h>
#include <Alieve/TPCSectorData.h>
#include <AliTPCParam.h>

#include <TStyle.h>
#include <TColor.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCSectorViz
//
// Base class for TPC raw-data visualization.
// See TPCSector2D and TPCSector3D for concrete implementations.

ClassImp(TPCSectorViz)

/**************************************************************************/

TPCSectorViz::TPCSectorViz(const Text_t* n, const Text_t* t) :
  Reve::RenderElement(fFrameColor),
  TNamed(n, t),

  fTPCData  (0),
  fSectorID (0),

  fMinTime   (0),
  fMaxTime   (450),
  fThreshold (5),
  fMaxVal    (80),

  fRnrInn   (kTRUE),
  fRnrOut1  (kTRUE),
  fRnrOut2  (kTRUE),

  fFrameColor ((Color_t) 4),
  fRnrFrame (kTRUE),
  fAutoTrans(kFALSE),
  fRTS      (1),

  fColorArray (0)
{}

TPCSectorViz::~TPCSectorViz()
{
  if(fTPCData) fTPCData->DecRefCount();
  delete [] fColorArray;
}

void TPCSectorViz::CopyVizParams(const TPCSectorViz& v)
{
  fMinTime   = v.fMinTime;
  fMaxTime   = v.fMaxTime;
  fThreshold = v.fThreshold;
  fMaxVal    = v.fMaxVal;

  fRnrInn    = v.fRnrInn;
  fRnrOut1   = v.fRnrOut1;
  fRnrOut2   = v.fRnrOut2;
}

/**************************************************************************/

void TPCSectorViz::SetDataSource(TPCData* data)
{
  if(data == fTPCData) return;
  if(fTPCData) fTPCData->DecRefCount();
  fTPCData = data;
  if(fTPCData) fTPCData->IncRefCount();
  IncRTS();
}

void TPCSectorViz::SetSectorID(Int_t id)
{
  if(id < 0)  id = 0;
  if(id > 35) id = 35;
  fSectorID = id;
  if(fAutoTrans)
    SetAutoTrans(kTRUE); // Force repositioning.
  IncRTS();
}

TPCSectorData* TPCSectorViz::GetSectorData() const
{
  return fTPCData ? fTPCData->GetSectorData(fSectorID) : 0;
}

/**************************************************************************/

void TPCSectorViz::SetThreshold(Short_t t)
{
  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  ClearColorArray();
  IncRTS();
}

void TPCSectorViz::SetMaxVal(Int_t mv)
{
  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  ClearColorArray();
  IncRTS();
}

/**************************************************************************/

void TPCSectorViz::SetAutoTrans(Bool_t trans) 
{
  fAutoTrans = trans;
  if(fAutoTrans) {
    fHMTrans.UnitTrans();

    using namespace TMath;
    Float_t c = Cos((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t s = Sin((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t z = TPCSectorData::GetZLength();
    Float_t d = -1;
    if(fSectorID >= 18) {
      z = -z;
      d = -d;
    }

    // column major
    fHMTrans[0]  = -c;
    fHMTrans[1]  = -s;
    fHMTrans[4]  = -s;
    fHMTrans[5]  =  c;
    fHMTrans[10] =  d;
    fHMTrans[14] =  z;
    fHMTrans[15] =  1;
  }
}

/**************************************************************************/

void TPCSectorViz::SetupColor(Int_t val, UChar_t* pixel) const
{
  using namespace TMath;
  Float_t div  = Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) Nint(nCol*(val - fThreshold)/div);

  ColorFromIdx(gStyle->GetColorPalette(Min(nCol - 1, cBin)), pixel);
}

void TPCSectorViz::ClearColorArray()
{
  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

void TPCSectorViz::SetupColorArray() const
{
  if(fColorArray)
    return;

  fColorArray = new UChar_t [4 * (fMaxVal - fThreshold + 1)];
  UChar_t* p = fColorArray;
  for(Int_t v=fThreshold; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);
}
