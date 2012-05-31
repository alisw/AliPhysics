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

#include <TEveTrans.h>
#include <TStyle.h>
#include <TMath.h>

//==============================================================================
//==============================================================================
// AliEveTPCSectorViz
//==============================================================================

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
  fMaxTime   (1000),
  fThreshold (5),
  fMaxVal    (128),

  fRnrInn   (kTRUE),
  fRnrOut1  (kTRUE),
  fRnrOut2  (kTRUE),

  fFrameColor (4),
  fRnrFrame   (kTRUE),
  fAutoTrans  (kFALSE),
  fRTS        (1),

  fColorArray (0)
{
  // Constructor.

   InitMainTrans();
}

AliEveTPCSectorViz::~AliEveTPCSectorViz()
{
  // Destructor.

  if (fTPCData) fTPCData->DecRefCount();
  delete [] fColorArray;
}

void AliEveTPCSectorViz::CopyVizParams(const TEveElement* el)
{
  // Copy basic viualization parameters from another TPCSectorViz.

  const AliEveTPCSectorViz* v = dynamic_cast<const AliEveTPCSectorViz*>(el);
  if (v) {
    fMinTime   = v->fMinTime;
    fMaxTime   = v->fMaxTime;
    fThreshold = v->fThreshold;
    fMaxVal    = v->fMaxVal;

    fRnrInn    = v->fRnrInn;
    fRnrOut1   = v->fRnrOut1;
    fRnrOut2   = v->fRnrOut2;
  }
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
  if (fAutoTrans)
  {
    using namespace TMath;
    Float_t c = Cos((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t s = Sin((fSectorID + 0.5)*20*Pi()/180 - PiOver2());
    Float_t z = AliEveTPCSectorData::GetZLength();

    InitMainTrans();
    TEveTrans& t = RefMainTrans();

    if (fSectorID < 18) {
      // column major
      t[0]  = -c;  t[1]  = -s;
      t[4]  = -s;  t[5]  =  c;
      t[10] = -1;
      t[14] =  z;
    } else {
      t[0]  =  c;  t[1]  =  s;
      t[4]  = -s;  t[5]  =  c;
      t[10] =  1;
      t[14] = -z;
    }
  }
}
void AliEveTPCSectorViz::SetUseTrans(Bool_t t)
{
  // Set flag spcifying if transformation matrix should be applied.

  RefMainTrans().SetUseTrans(t);
}

/******************************************************************************/

void AliEveTPCSectorViz::SetupColor(Int_t val, UChar_t* pixel) const
{
  // Set pixel color to represent signal val from the current palette.

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = TMath::Nint(nCol*(val - fThreshold)/div);

  TEveUtil::ColorFromIdx(gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin)), pixel);
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
