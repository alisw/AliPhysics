#include "TPCSector2D.h"

#include <Alieve/TPCData.h>
#include <Alieve/TPCSectorData.h>

#include <AliTPCParam.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(TPCSector2D)

/**************************************************************************/

void TPCSector2D::Init()
{
  fTPCData   = 0;

  fSectorID  = 0;
  fShowMax   = true;
  fMinTime   = 0;
  fMaxTime   = 1;
  fthreshold  = 5;
  fMaxVal    = 80;

  fRnrFrame   = true;
  fUseTexture = true;

  fTrans      = false;
}

TPCSector2D::~TPCSector2D()
{
  if(fTPCData) fTPCData->DecRefCount();
}

/**************************************************************************/

void TPCSector2D::SetDataSource(TPCData* data)
{
  if(data == fTPCData) return;
  if(fTPCData) fTPCData->DecRefCount();
  fTPCData = data;
  if(fTPCData) fTPCData->IncRefCount();
  ++fRTS;
}

void TPCSector2D::SetSectorID(Int_t segment)
{
  if(segment < 0 ) segment = 0;
  if(segment > 35) segment = 35;
  fSectorID = segment;
  SetName(Form("TPCSector2D %d", fSectorID));
  ++fRTS;
}

/**************************************************************************/

void TPCSector2D::ComputeBBox()
{
  const TPCSectorData::SegmentInfo&  iSeg = TPCSectorData::GetInnSeg();
  const TPCSectorData::SegmentInfo& o2Seg = TPCSectorData::GetOut2Seg();

  bbox_init();
  Float_t w = o2Seg.GetNMaxPads()*o2Seg.GetPadWidth()/2;
  fBBox[0] = -w;
  fBBox[1] =  w;
  fBBox[2] =  iSeg.GetRLow();
  fBBox[3] =  o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  fBBox[4] = -0.5; // Fake z-width to 1 cm.
  fBBox[5] =  0.5;
}

/**************************************************************************/

void TPCSector2D::SetTrans(Bool_t trans) 
{
  fTrans = trans;
  if(fTrans) {
    for (Int_t k = 0; k< 16; k++)
      fMatrix[k] = 0.;

    Float_t c = TMath::Cos((fSectorID + 0.5)*20*TMath::Pi()/180 - TMath::Pi()/2);
    Float_t s = TMath::Sin((fSectorID + 0.5)*20*TMath::Pi()/180 - TMath::Pi()/2);
    Float_t z = TPCSectorData::GetParam().GetZLength();
    if(fSectorID >= 18) z = -z;
  
    // column major
    fMatrix[0]  = -c;
    fMatrix[1]  = -s;
    fMatrix[4]  = -s;
    fMatrix[5]  =  c;
    fMatrix[10] = -1;
    fMatrix[14] =  z;
    fMatrix[15] =  1;
  }
}

/**************************************************************************/

void TPCSector2D::Paint(Option_t* )
{
  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = 1;
  buffer.fTransparency = 0;
  buffer.fLocalFrame   = fTrans; 
  if (fTrans)
    memcpy(buffer.fLocalMaster, fMatrix, 16*sizeof(Double_t));
  buffer.SetSectionsValid(TBuffer3D::kCore);
   
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    // printf("TPCSector2D::Paint viewer was happy with Core buff3d.\n");
    return;
  }

  printf("TPCSector2D::Paint only GL supported.\n");
  return;

  /*
    if (reqSections & TBuffer3D::kRawSizes) {
    Int_t nbPnts = fQuads.size()*4;
    Int_t nbSegs = nbPnts;
    if (!buffer.SetRawSizes(nbPnts, 3*nbPnts, nbSegs, 3*nbSegs, fQuads.size(), fQuads.size()*6)) {
    return;
    }
    buffer.SetSectionsValid(TBuffer3D::kRawSizes); 
    }

    if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
    // Points
    Int_t pidx = 0;
    for (std::vector<Quad>::iterator i=fQuads.begin(); i!=fQuads.end(); ++i) {
    for (Int_t k = 0; k < 12; k++ ){
    buffer.fPnts[pidx] = (*i).vertices[k]; 
    pidx++;
    }
    }

    // Segments
    Int_t sidx = 0;
    for (Int_t q = 0; q < fQuads.size(); ++q) {
    for (Int_t s = 0; s < 4; ++s ) {
    buffer.fSegs[3*sidx ] = 4; 
    buffer.fSegs[3*sidx+1] = sidx;
    if (s == 3)
    buffer.fSegs[3*sidx+2] = q*4;
    else
    buffer.fSegs[3*sidx+2] = sidx + 1;
    sidx ++;
    }
    }

    // Polygons
    for (Int_t q = 0; q < fQuads.size(); ++q) {
    buffer.fPols[6*q] = fQuads[q].color;   
    buffer.fPols[6*q +1] = 4;
    buffer.fPols[6*q +2] = 4*q +0;
    buffer.fPols[6*q +3] = 4*q +1;
    buffer.fPols[6*q +4] = 4*q +2;
    buffer.fPols[6*q +5] = 4*q +3;
    }

    buffer.SetSectionsValid(TBuffer3D::kRaw);
    buffer.fColor = 5;
    }
   
  */
  // gPad->GetViewer3D()->AddObject(buffer);
}
