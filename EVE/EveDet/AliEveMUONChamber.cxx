// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include "AliEveMUONChamber.h"

#include <EveDet/AliEveMUONData.h>
#include <EveDet/AliEveMUONChamberData.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>
#include <TStyle.h>
#include <TMath.h>


//______________________________________________________________________________
// AliEveMUONChamber
//

ClassImp(AliEveMUONChamber)

//______________________________________________________________________________
AliEveMUONChamber::AliEveMUONChamber(Int_t id, const Text_t* n, const Text_t* t) :
  TEveElement(fFrameColor),
  TNamed(n,t),
  fMUONData(0),
  fFrameColor(2),
  fRTS(1),
  fChamberID(0),
  fQuadSet1(n,t),
  fQuadSet2(n,t),
  fPointSet1(n),
  fPointSet2(n),
  fThreshold(0),
  fMaxVal(4096),
  fClusterSize(5),
  fHitSize(5),
  fColorArray(0)
{
  //
  // constructor
  //

  Char_t name[256];
  if (id < 10) {
    snprintf(name,256,"Chamber %02d (trac)",id);
  } else {
    snprintf(name,256,"Chamber %02d (trig)",id);
  }
  SetName(name);

  ComputeBBox();

}

//______________________________________________________________________________
AliEveMUONChamber::~AliEveMUONChamber()
{
  //
  // destructor
  //

  if(fMUONData) fMUONData->DecRefCount();

}

//______________________________________________________________________________
void AliEveMUONChamber::ComputeBBox()
{
  //
  // bounding box
  //

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  bbox_init();
#else
  BBoxInit();
#endif

  fBBox[0] = - 400.0;
  fBBox[1] = + 400.0;
  fBBox[2] = - 400.0;
  fBBox[3] = + 400.0;
  fBBox[4] = -1800.0;
  fBBox[5] = + 500.0;

  Float_t* b1 = fQuadSet1.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b1[i] = fBBox[i]; }
  Float_t* b2 = fQuadSet2.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b2[i] = fBBox[i]; }
  Float_t* b3 = fPointSet1.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b3[i] = fBBox[i]; }
  Float_t* b4 = fPointSet2.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b4[i] = fBBox[i]; }

}

//______________________________________________________________________________
void AliEveMUONChamber::Paint(Option_t*)
{
  //
  // draw...
  //

  if(fRnrSelf == kFALSE)
    return;

  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  buffer.fID           = this;
  buffer.fColor        = 2;
  buffer.fTransparency = 0;
  buffer.fLocalFrame   = 0;

  buffer.SetSectionsValid(TBuffer3D::kCore);
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    //printf("AliEveMUONChamber::Paint viewer was happy with Core buff3d.\n");
    return;
  }

  printf("AliEveMUONChamber::Paint only GL supported.\n");
  return;

}

//______________________________________________________________________________
void AliEveMUONChamber::SetThreshold(Short_t t)
{
  //
  // digits amplitude threshold
  //

  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  ClearColorArray();
  IncRTS();

}

//______________________________________________________________________________
void AliEveMUONChamber::SetMaxVal(Int_t mv)
{
  //
  // digits amplitude maximum value
  //

  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  ClearColorArray();
  IncRTS();

}

//______________________________________________________________________________
void AliEveMUONChamber::SetClusterSize(Int_t size)
{
  //
  // cluster point size
  //

  fClusterSize = TMath::Max(1, size);
  IncRTS();

}

//______________________________________________________________________________
void AliEveMUONChamber::SetHitSize(Int_t size)
{
  //
  // hit point size
  //

  fHitSize = TMath::Max(1, size);
  IncRTS();

}

//______________________________________________________________________________
void AliEveMUONChamber::SetupColor(Int_t val, UChar_t* pixel) const
{
  //
  // RGBA color for amplitude "val"
  //

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) TMath::Nint(nCol*(val - fThreshold)/div);

  TEveUtil::TEveUtil::ColorFromIdx(gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin)), pixel);

}

//______________________________________________________________________________
Int_t AliEveMUONChamber::ColorIndex(Int_t val) const
{
  //
  // index color
  //

  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) TMath::Nint(nCol*(val - fThreshold)/div);

  return gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin));

}

//______________________________________________________________________________
void AliEveMUONChamber::SetupColorArray() const
{
  //
  // build array of colors
  //

  if(fColorArray)
    return;

  fColorArray = new UChar_t [4 * (fMaxVal - fThreshold + 1)];
  UChar_t* p = fColorArray;
  for(Int_t v=fThreshold; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);

}

//______________________________________________________________________________
void AliEveMUONChamber::ClearColorArray()
{
  //
  // delete array of colors
  //

  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

//______________________________________________________________________________
void AliEveMUONChamber::SetDataSource(AliEveMUONData* data)
{
  // Set source of data.

  if (data == fMUONData) return;
  if(fMUONData) fMUONData->DecRefCount();
  fMUONData = data;
  if(fMUONData) fMUONData->IncRefCount();
  IncRTS();
}

//______________________________________________________________________________
AliEveMUONChamberData* AliEveMUONChamber::GetChamberData() const
{
  // Return source of data.

  return fMUONData ? fMUONData->GetChamberData(fChamberID) : 0;
}

//______________________________________________________________________________
void AliEveMUONChamber::UpdateQuads()
{
  // Update digit representation.

  fQuadSet1.Reset(TEveQuadSet::kQT_RectangleXY, kTRUE, 32);
  fQuadSet2.Reset(TEveQuadSet::kQT_RectangleXY, kTRUE, 32);
  fPointSet1.Reset();
  fPointSet2.Reset();

  AliEveMUONChamberData* data = GetChamberData();

  Float_t *buffer;
  Float_t x0, y0, z, w, h, clsq;
  Int_t charge, cathode, nDigits, nClusters, nHits, oldSize, ic1, ic2;
  Double_t clsX, clsY, clsZ;
  Float_t hitX, hitY, hitZ;

  if (data != 0) {

    SetupColorArray();

    // digits

    nDigits = data->GetNDigits();

    for (Int_t id = 0; id < nDigits; id++) {

      buffer = data->GetDigitBuffer(id);

      x0 = buffer[0]-buffer[2];
      y0 = buffer[1]-buffer[3];
      w  = 2*buffer[2];
      h  = 2*buffer[3];
      z  = buffer[4];
      charge = (Int_t)buffer[5];
      cathode = (Int_t)buffer[6];

      if (charge <= fThreshold) continue;

      if (cathode == 0) {

	fQuadSet1.AddQuad(x0, y0, z, w, h);
	fQuadSet1.QuadColor(ColorIndex(charge));

      }

      if (cathode == 1) {

	fQuadSet2.AddQuad(x0, y0, z, w, h);
	fQuadSet2.QuadColor(ColorIndex(charge));

      }

    } // end digits loop

    // clusters

    nClusters = data->GetNClusters()/2;  // only one cathode plane
    oldSize = fPointSet1.GrowFor(nClusters);
    ic1 = ic2 = 0;
    for (Int_t ic = 0; ic < (nClusters*2); ic++) {

      buffer = data->GetClusterBuffer(ic);

      clsX    = (Double_t)buffer[0];
      clsY    = (Double_t)buffer[1];
      clsZ    = (Double_t)buffer[2];
      clsq    = buffer[3];
      cathode = (Int_t)buffer[4];

      if (cathode == 0) {
	fPointSet1.SetPoint(ic1,clsX,clsY,clsZ);
	ic1++;
      }

    } // end clusters loop

    // hits

    nHits = data->GetNHits();
    oldSize = fPointSet2.GrowFor(nHits);
    for (Int_t ih = 0; ih < nHits; ih++) {
      buffer = data->GetHitBuffer(ih);
      hitX = buffer[0];
      hitY = buffer[1];
      hitZ = buffer[2];
      fPointSet2.SetPoint(ih,hitX,hitY,hitZ);
    }

  } // end data

}

//______________________________________________________________________________
void AliEveMUONChamber::SetChamberID(Int_t id)
{
  // Set id of chamber to display.

  if (id <  0) id = 0;
  if (id > 13) id = 13;

  fChamberID = id;
  IncRTS();
}

