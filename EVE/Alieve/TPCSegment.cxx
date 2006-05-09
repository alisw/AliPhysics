#include "TPCSegment.h"


using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(TPCSegment)

/**************************************************************************/

void TPCSegment::Init()
{
  fID = 0;
  fInfo = 0;

  fTrans = false;

  fRnrFrame = true;
  fUseTexture = true;
  fTreshold = 1;
  fShowMax = true;

  fMinTime   = 0;
  fMaxTime   = 1;
  fTreshold  = 5;
  fMaxVal    = 80;
}

TPCSegment::~TPCSegment()
{
  if(fInfo) fInfo->DecRefCount();
}

/**************************************************************************/

void TPCSegment::SetInfo(TPCDigitsInfo* info)
{
  if(fInfo) fInfo->DecRefCount();
  fInfo = info;
  if(fInfo) fInfo->IncRefCount();
}

void TPCSegment::SetSegmentID(Int_t segment)
{
  if(segment < 0 ) segment = 0;
  if(segment > 36) segment = 36;
  fID = segment;
  SetName(Form("TPCSegment %d", fID));
  ++fRTS;
}

/**************************************************************************/

void TPCSegment::ComputeBBox()
{
  Float_t b = fInfo->fInnSeg.fRlow;
  Float_t w = fInfo->fOut2Seg.fNMaxPads* fInfo->fOut2Seg.fPadWidth/2;
  Float_t h = fInfo->fOut2Seg.fRlow +
    fInfo->fOut2Seg.fNRows* fInfo->fOut2Seg.fPadLength - fInfo->fInnSeg.fRlow;

  bbox_init();
  fBBox[0] = -w;   fBBox[1] = w;
  fBBox[2] =  b;   fBBox[3] = b + h;
  fBBox[4] = -0.5; fBBox[5] = 0.5;   // Fake z-width to 1 cm.
}

/**************************************************************************/

void TPCSegment::SetTrans(Bool_t trans) 
{
  fTrans = trans;
  if(fTrans) {
    for (Int_t k = 0; k< 16; k++)
      fMatrix[k] = 0.;

    Float_t z, s, c;
    if(fID < 18) {
      z =  fInfo->fParameter->GetZLength();
    } else {
      z = -fInfo->fParameter->GetZLength();
    } 
  
    // column major ii
    fMatrix[14] = z;
    fMatrix[15] = 1;

    c = TMath::Cos((fID + 0.5)*20*TMath::Pi()/180 - TMath::Pi()/2);
    s = TMath::Sin((fID + 0.5)*20*TMath::Pi()/180 - TMath::Pi()/2);
  
    fMatrix[0] = -c;
    fMatrix[1] = -s;
    fMatrix[4] = -s;
    fMatrix[5] =  c;
    fMatrix[10] = -1;
  }
}

/**************************************************************************/
void TPCSegment::Paint(Option_t* )
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
   
  // We fill kCore on first pass and try with viewer
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    // printf("TPCSegment::Paint viewer was happy with Core buff3d.\n");
    return;
  }
  printf("TPCSegment::Paint only GL supported.\n");
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
