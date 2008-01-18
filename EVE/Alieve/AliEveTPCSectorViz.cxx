// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSectorViz.h"

#include <Alieve/AliEveTPCData.h>
#include <Alieve/AliEveTPCSectorData.h>
#include <AliTPCParam.h>

#include <TStyle.h>
#include <TColor.h>


//______________________________________________________________________
// AliEveTPCSectorViz
//
// Base class for TPC raw-data visualization.
// See AliEveTPCSector2D and AliEveTPCSector3D for concrete implementations.

ClassImp(AliEveTPCSectorViz)

/**************************************************************************/

AliEveTPCSectorViz::AliEveTPCSectorViz(const Text_t* n, const Text_t* t) :
  TEveElement(fFrameColor),
  TNamed(n, t),

  fTPCData  (0),
  fSectorID (0),

  fMinTime   (0),
  fMaxTime   (450),
  fThreshold (5),
  fMaxVal    (128),

  fRnrInn   (kTRUE),
  fRnrOut1  (kTRUE),
  fRnrOut2  (kTRUE),

  fFrameColor ((Color_t) 4),
  fRnrFrame   (kTRUE),
  fAutoTrans  (kFALSE),
  fRTS        (1),

  fColorArray (0)
{}

AliEveTPCSectorViz::~AliEveTPCSectorViz()
{
  if(fTPCData) fTPCData->DecRefCount();
  delete [] fColorArray;
}

void AliEveTPCSectorViz::CopyVizParams(const AliEveTPCSectorViz& v)
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

void AliEveTPCSectorViz::SetDataSource(AliEveTPCData* data)
{
  if(data == fTPCData) return;
  if(fTPCData) fTPCData->DecRefCount();
  fTPCData = data;
  if(fTPCData) fTPCData->IncRefCount();
  IncRTS();
}

void AliEveTPCSectorViz::SetSectorID(Int_t id)
{
  if(id < 0)  id = 0;
  if(id > 35) id = 35;
  fSectorID = id;
  if(fAutoTrans)
    SetAutoTrans(kTRUE); // Force repositioning.
  IncRTS();
}

AliEveTPCSectorData* AliEveTPCSectorViz::GetSectorData() const
{
  return fTPCData ? fTPCData->GetSectorData(fSectorID) : 0;
}

/**************************************************************************/

void AliEveTPCSectorViz::SetThreshold(Short_t t)
{
  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  ClearColorArray();
  IncRTS();
}

void AliEveTPCSectorViz::SetMaxVal(Int_t mv)
{
  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  ClearColorArray();
  IncRTS();
}

/**************************************************************************/

void AliEveTPCSectorViz::SetAutoTrans(Bool_t trans)
{
  fAutoTrans = trans;
  if(fAutoTrans) {
    fHMTrans.UnitTrans();

    using namespace TMath;
    Float_t c = Cos((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t s = Sin((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t z = AliEveTPCSectorData::GetZLength();
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

void AliEveTPCSectorViz::SetupColor(Int_t val, UChar_t* pixel) const
{
  using namespace TMath;
  Float_t div  = Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) Nint(nCol*(val - fThreshold)/div);

  TEveUtil::TEveUtil::ColorFromIdx(gStyle->GetColorPalette(Min(nCol - 1, cBin)), pixel);
}

void AliEveTPCSectorViz::ClearColorArray()
{
  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

void AliEveTPCSectorViz::SetupColorArray() const
{
  if(fColorArray)
    return;

  fColorArray = new UChar_t [4 * (fMaxVal - fThreshold + 1)];
  UChar_t* p = fColorArray;
  for(Int_t v=fThreshold; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);
}
