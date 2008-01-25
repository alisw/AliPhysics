// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSectorViz.h"

#include <EveDet/AliEveTPCData.h>
#include <EveDet/AliEveTPCSectorData.h>
#include <AliTPCParam.h>

#include <TStyle.h>
#include <TColor.h>


//______________________________________________________________________________
//
// Base class for TPC raw-data visualization.
// See AliEveTPCSector2D and AliEveTPCSector3D for concrete implementations.

ClassImp(AliEveTPCSectorViz)

/******************************************************************************/

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
  fHMTrans    (),
  fAutoTrans  (kFALSE),
  fRTS        (1),

  fColorArray (0)
{
  // Constructor.
}

AliEveTPCSectorViz::~AliEveTPCSectorViz()
{
  // Destructor.

  if (fTPCData) fTPCData->DecRefCount();
  delete [] fColorArray;
}

void AliEveTPCSectorViz::CopyVizParams(const AliEveTPCSectorViz& v)
{
  // Copy basic viualization parameters from another TPCSectorViz.

  fMinTime   = v.fMinTime;
  fMaxTime   = v.fMaxTime;
  fThreshold = v.fThreshold;
  fMaxVal    = v.fMaxVal;

  fRnrInn    = v.fRnrInn;
  fRnrOut1   = v.fRnrOut1;
  fRnrOut2   = v.fRnrOut2;
}

/******************************************************************************/

void AliEveTPCSectorViz::SetDataSource(AliEveTPCData* data)
{
  // Set the data source.

  if (data == fTPCData) return;
  if (fTPCData) fTPCData->DecRefCount();
  fTPCData = data;
  if (fTPCData) fTPCData->IncRefCount();
  IncRTS();
}

void AliEveTPCSectorViz::SetSectorID(Int_t id)
{
  // Set sector id.

  if (id < 0)  id = 0;
  if (id > 35) id = 35;
  fSectorID = id;
  if (fAutoTrans)
    SetAutoTrans(kTRUE); // Force repositioning.
  IncRTS();
}

AliEveTPCSectorData* AliEveTPCSectorViz::GetSectorData() const
{
  // Get sector-data.

  return fTPCData ? fTPCData->GetSectorData(fSectorID) : 0;
}

/******************************************************************************/

void AliEveTPCSectorViz::SetThreshold(Short_t t)
{
  // Set visualization threshold.

  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  ClearColorArray();
  IncRTS();
}

void AliEveTPCSectorViz::SetMaxVal(Int_t mv)
{
  // Set visualization max signal value.
  // Signals above this will have the same color.

  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  ClearColorArray();
  IncRTS();
}

/******************************************************************************/

void AliEveTPCSectorViz::SetAutoTrans(Bool_t trans)
{
  // Set automatic update of transformation matrix.
  // The position is calculated immediately.

  fAutoTrans = trans;
  if (fAutoTrans) {
    fHMTrans.UnitTrans();

    using namespace TMath;
    Float_t c = Cos((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t s = Sin((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t z = AliEveTPCSectorData::GetZLength();
    Float_t d = -1;
    if (fSectorID >= 18) {
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

/******************************************************************************/

void AliEveTPCSectorViz::SetupColor(Int_t val, UChar_t* pixel) const
{
  // Set pixel color to represent signal val from the current palette.

  using namespace TMath;
  Float_t div  = Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) Nint(nCol*(val - fThreshold)/div);

  TEveUtil::ColorFromIdx(gStyle->GetColorPalette(Min(nCol - 1, cBin)), pixel);
}

void AliEveTPCSectorViz::ClearColorArray()
{
  // Clear cached color array.

  if (fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

void AliEveTPCSectorViz::SetupColorArray() const
{
  // Initialize cached color array.

  if (fColorArray)
    return;

  fColorArray = new UChar_t [4 * (fMaxVal - fThreshold + 1)];
  UChar_t* p = fColorArray;
  for (Int_t v=fThreshold; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);
}
