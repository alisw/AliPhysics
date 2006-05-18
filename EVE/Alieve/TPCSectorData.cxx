// $Header$

#include "TPCSectorData.h"
#include <AliTPCParamSR.h>

#include <string.h>

//______________________________________________________________________
// TPCSectorData
//
// Stores data from a fiven TPC sector.
//
// Row addresses grow linearly by radius, there is no separation on
// inner/outer segments. The SegmentInfo objects can be used to get
// information about low-level segments.
//
// A lot of TPC-sector information is stored as static data.
//
// For accessing data, see for example TPCSector2DGL::CreateTexture()
// and LoadPadrow().
//

using namespace Reve;
using namespace Alieve;

ClassImp(TPCSectorData);

AliTPCParam* TPCSectorData::fgParam    = 0;
Int_t        TPCSectorData::fgNAllRows = 0;
Int_t        TPCSectorData::fgNAllPads = 0;
Int_t*       TPCSectorData::fgRowBegs  = 0;

TPCSectorData::SegmentInfo TPCSectorData::fgInnSeg;
TPCSectorData::SegmentInfo TPCSectorData::fgOut1Seg;
TPCSectorData::SegmentInfo TPCSectorData::fgOut2Seg;

TPCSectorData::SegmentInfo* TPCSectorData::fgSegInfoPtrs[3] = {0};

/**************************************************************************/

void TPCSectorData::InitStatics()
{
  if(fgParam != 0) return;

  fgParam    = new AliTPCParamSR;
  fgNAllRows = fgParam->GetNRowLow() + fgParam->GetNRowUp();
  fgNAllPads = 0;
  fgRowBegs  = new Int_t[fgNAllRows + 1];

  Int_t row = 0;
  for(Int_t i=0; i<fgParam->GetNRowLow(); ++i, ++row) {
    fgRowBegs[row] = fgNAllPads;
    fgNAllPads += fgParam->GetNPadsLow(i);
  }
  for(Int_t i=0; i<fgParam->GetNRowUp(); ++i, ++row) {
    fgRowBegs[row] = fgNAllPads;
    fgNAllPads += fgParam->GetNPadsUp(i);
  }
  fgRowBegs[fgNAllRows] = fgNAllPads;


  // Fill SegmentInfos, used by rendering classes.

  // General paramameters
  fgInnSeg.fPadWidth   = fgParam->GetInnerPadPitchWidth();
  fgInnSeg.fPadHeight  = fgParam->GetInnerPadPitchLength();
  fgInnSeg.fRLow       = fgParam->GetPadRowRadiiLow(0);
  fgInnSeg.fNRows      = fgParam->GetNRowLow();
  fgInnSeg.fFirstRow   = 0;
  fgInnSeg.fLastRow    = fgInnSeg.fNRows - 1;
  fgInnSeg.fNMaxPads   = fgParam->GetNPadsLow(fgInnSeg.fNRows - 1);
  fgSegInfoPtrs[0]     = &fgInnSeg;

  fgOut1Seg.fPadWidth  = fgParam->GetOuterPadPitchWidth();
  fgOut1Seg.fPadHeight = fgParam->GetOuter1PadPitchLength();
  fgOut1Seg.fRLow      = fgParam->GetPadRowRadiiUp(0);
  fgOut1Seg.fNRows     = fgParam->GetNRowUp1();
  fgOut1Seg.fFirstRow  = fgInnSeg.fNRows;
  fgOut1Seg.fLastRow   = fgOut1Seg.fFirstRow + fgOut1Seg.fNRows - 1;
  fgOut1Seg.fNMaxPads  = fgParam->GetNPadsUp(fgOut1Seg.fNRows - 1);
  fgSegInfoPtrs[1]     = &fgOut1Seg;
 
  fgOut2Seg.fPadWidth  = fgParam->GetOuterPadPitchWidth();
  fgOut2Seg.fPadHeight = fgParam->GetOuter2PadPitchLength();
  fgOut2Seg.fRLow      = fgParam->GetPadRowRadiiUp(fgOut1Seg.fNRows);
  fgOut2Seg.fNRows     = fgParam->GetNRowUp() - fgOut1Seg.fNRows;
  fgOut2Seg.fFirstRow  = fgOut1Seg.fLastRow + 1;
  fgOut2Seg.fLastRow   = fgOut2Seg.fFirstRow + fgOut2Seg.fNRows - 1;
  fgOut2Seg.fNMaxPads  = fgParam->GetNPadsUp(fgParam->GetNRowUp() - 1);
  fgSegInfoPtrs[2]     = &fgOut2Seg;

  // Set stepsize array
  Int_t k, npads;
  // Inn
  k=0, npads = fgParam->GetNPadsLow(0);
  for (int row = 0; row < fgInnSeg.fNRows; ++row) {
    if (fgParam->GetNPadsLow(row) > npads) {
      npads = fgParam->GetNPadsLow(row);
      fgInnSeg.fYStep[k] = row*fgInnSeg.fPadHeight + fgInnSeg.fRLow;
      k++;
    }
  }
  fgInnSeg.fNYSteps = k;
  // Out1 seg 
  k=0; npads = fgParam->GetNPadsUp(0);
  for (int row = 0; row < fgOut1Seg.fNRows; ++row) {
    if (fgParam->GetNPadsUp(row) > npads) {
      npads = fgParam->GetNPadsUp(row);
      fgOut1Seg.fYStep[k] = row*fgOut1Seg.fPadHeight + fgOut1Seg.fRLow ;
      k++;
    }
  }
  fgOut1Seg.fNYSteps = k;
  // Out2 seg
  k=0; npads = fgParam->GetNPadsUp(fgOut1Seg.fNRows);
  for (int row = fgOut1Seg.fNRows; row < fgParam->GetNRowUp() ;row++ ) {
    if (fgParam->GetNPadsUp(row) > npads) {
      npads = fgParam->GetNPadsUp(row);
      fgOut2Seg.fYStep[k] = (row - fgOut1Seg.fNRows)*fgOut2Seg.fPadHeight + fgOut2Seg.fRLow ;
      k++;
    }
  }
  fgOut2Seg.fNYSteps = k;
}

Int_t TPCSectorData::GetNPadsInRow(Int_t row)
{
  if(row < 0 || row >= fgNAllRows) return 0;
  return fgRowBegs[row + 1] - fgRowBegs[row];
}

const TPCSectorData::SegmentInfo& TPCSectorData::GetSeg(Int_t seg)
{
  static const SegmentInfo null;

  if(seg < 0 || seg > 2)
    return null;
  else
    return *fgSegInfoPtrs[seg];
}

/**************************************************************************/
// True member functions start here.
/**************************************************************************/

void TPCSectorData::NewBlock()
{
  fBlocks.push_back(new Short_t[fBlockSize]);
  fBlockPos = 0;
}

/**************************************************************************/

TPCSectorData::TPCSectorData(Int_t sector, Int_t bsize) :
  fSectorID(sector), fNPadsFilled(0),
  fBlockSize(bsize), fBlockPos(0),
  fCurrentRow(0), fCurrentPad(0), fCurrentPos(0)
{
  if(fgParam == 0) InitStatics();

  fPads.assign(fgNAllPads, PadData());
  fBlocks.reserve(16);
  fBlockPos = fBlockSize; // Enforce creation of a new block.
}


TPCSectorData::~TPCSectorData()
{
  for(std::vector<Short_t*>::iterator b=fBlocks.begin(); b!=fBlocks.end(); ++b)
    delete [] *b;
}

/**************************************************************************/

void TPCSectorData::Print(Option_t* /*opt*/) const
{
  printf("TPCSectorData sector=%d, NPadsFilled=%d, NBlocks=%d, BlockPos=%d\n",
	 fSectorID, fNPadsFilled, fBlocks.size(), fBlockPos);
}

/**************************************************************************/

void TPCSectorData::BeginPad(Int_t row, Int_t pad, Bool_t reverseTime)
{
  fCurrentRow = row;
  fCurrentPad = pad;
  if(reverseTime) {
    fCurrentPos  = 2046;
    fCurrentStep = -2;
  } else {
    fCurrentPos  = 0;
    fCurrentStep = 2;
  }
}

void TPCSectorData::EndPad()
{
  Short_t *beg, *end;
  if(fCurrentStep > 0) {
    beg = fPadBuffer;
    end = fPadBuffer + fCurrentPos;
  } else {
    beg = fPadBuffer + fCurrentPos + 2;
    end = fPadBuffer + 2048;
  }
  Short_t* wpos = beg;
  Short_t* rpos = beg;

  // Compress pad buffer
  while(rpos < end) {
    Short_t* spos = rpos;
    Short_t  t    = spos[0];
    while(true) {
      rpos += 2;
      if(rpos >= end || *rpos > t + 1 || t == 0)
	break;
      ++t;
    }
    Short_t n = t - spos[0] + 1;
    if(n == 1) {
      wpos[0] = -spos[0];
      wpos[1] =  spos[1];
      wpos += 2;
    } else {
      wpos[0] = spos[0];
      wpos[2] = spos[1];
      wpos[1] = n;
      wpos += 3; spos += 3;
      while(--n) {
	*wpos = *spos;
	++wpos;	spos += 2;
      }
    }
  }

  // Copy buffer to storage, set PadData
  if(wpos > beg) {
    Short_t len = wpos - beg;
    if(len > fBlockSize - fBlockPos)
      NewBlock();
    Short_t *dest = fBlocks.back() + fBlockPos;
    memcpy(dest, beg, len*sizeof(Short_t));
    fBlockPos += len;

    PadData& pad = fPads[PadIndex(fCurrentRow, fCurrentPad)];
    pad.SetDataLength(dest, len);
  }

  ++fNPadsFilled;
}

/**************************************************************************/

const TPCSectorData::PadData& TPCSectorData::GetPadData(Int_t padAddr)
{
  static const PadData null;

  if(padAddr < 0 || padAddr >= fgNAllPads) return null;
  return fPads[padAddr];
}

const TPCSectorData::PadData& TPCSectorData::GetPadData(Int_t row, Int_t pad)
{
  static const PadData null;

  Int_t np = GetNPadsInRow(row);
  if(np == 0 || pad < 0 || pad >= np) return null;
  return GetPadData(fgRowBegs[row] + pad);
}

TPCSectorData::PadIterator TPCSectorData::MakePadIterator(Int_t padAddr, Short_t thr)
{
  return PadIterator(GetPadData(padAddr), thr);
}

TPCSectorData::PadIterator TPCSectorData::MakePadIterator(Int_t row, Int_t pad, Short_t thr)
{
  return PadIterator(GetPadData(row, pad), thr);
}

TPCSectorData::RowIterator TPCSectorData::MakeRowIterator(Int_t row, Short_t thr)
{
  Short_t npads = GetNPadsInRow(row);
  if(npads > 0)
    return RowIterator(&fPads[fgRowBegs[row]], npads, thr);
  else
    return RowIterator(0, 0);
}

/**************************************************************************/
// TPCSectorData::PadData
/**************************************************************************/

void TPCSectorData::PadData::Print(Option_t* /*opt*/)
{
  printf("addr=%p, len=%hd>\n", (void*)fData, fLength);
  for(Int_t i=0; i<fLength; ++i)
    printf("  %3d %hd\n", i, fData[i]);
}

/**************************************************************************/
// TPCSectorData::PadIterator
/**************************************************************************/

Bool_t TPCSectorData::PadIterator::Next()
{
  if(fPos >= fEnd) return kFALSE;
  if(fNChunk > 0) {
    ++fTime;
    --fNChunk;
    fSignal = *fPos; ++fPos;
  } else {
    fTime = fPos[0];
    if(fTime <= 0) {
      fTime   = -fTime;
      fSignal = fPos[1];
      fPos += 2;
    } else {
      fNChunk = fPos[1] - 1;
      fSignal = fPos[2];
      fPos += 3;
    }
  }
  return (fSignal > fThreshold) ? kTRUE : Next();
}

void TPCSectorData::PadIterator::Reset()
{
  // Return to the beginning of the pad-data. Must call Next() to get to
  // the first stored signal.

  fPos    = fBeg;
  fTime   = -1;
  fSignal = -1;
  fNChunk = 0;
}

void TPCSectorData::PadIterator::Reset(const PadData& pd)
{
  // Reinitialize to new pad-data. Must call Next() to get to the
  // first stored signal.

  fBeg = pd.Data();
  fEnd = pd.Data() + pd.Length();
  fPos = pd.Data();
  Reset();
}

void TPCSectorData::PadIterator::Test()
{
  while(Next())
    printf("  %3d %d\n", fTime, fSignal);
}

/**************************************************************************/
// TPCSectorData::RowIterator
/**************************************************************************/

Bool_t TPCSectorData::RowIterator::NextPad()
{
  ++fPad;
  if(fPad >= fNPads) return kFALSE;
  Reset(fPadArray[fPad]);
  return kTRUE;
}

void TPCSectorData::RowIterator::ResetRow()
{
  // Return to the beginning of the row. Must call NextPad() to get to
  // the zeroth pad.

  fPad = -1;
}

void TPCSectorData::RowIterator::ResetRow(const PadData* first, Short_t npads)
{
  // Reinitialize to another pad-data array. Must call NextPad() to
  // get to the zeroth pad.

  fPadArray =  first;
  fNPads    =  npads;
  fPad      = -1;
}

void TPCSectorData::RowIterator::Test()
{
  while(NextPad()) {
    printf("Pad %d\n", fPad);
    PadIterator::Test();
  }
}

/**************************************************************************/
// TPCSectorData::SegmentInfo
/**************************************************************************/

ClassImp(TPCSectorData::SegmentInfo);

TPCSectorData::SegmentInfo::SegmentInfo()
{
  memset(this, sizeof(SegmentInfo), 0);
}
