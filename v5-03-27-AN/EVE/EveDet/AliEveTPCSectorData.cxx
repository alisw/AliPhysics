// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTPCSectorData.h"
#include <AliTPCParamSR.h>

#include <set>
#include <string.h>

//==============================================================================
//==============================================================================
// AliEveTPCSectorData
//==============================================================================

//______________________________________________________________________________
//
// Stores data from a fiven TPC sector.
//
// Row addresses grow linearly by radius, there is no separation on
// inner/outer segments. The SegmentInfo objects can be used to get
// information about low-level segments.
//
// A lot of TPC-sector geometry information is stored as static data.
//
// For accessing data, see for example AliEveTPCSector2DGL::CreateTexture()
// and LoadPadrow().
//

ClassImp(AliEveTPCSectorData)

AliTPCParam* AliEveTPCSectorData::fgParam    = 0;
Float_t      AliEveTPCSectorData::fgZLength  = 0;
Int_t        AliEveTPCSectorData::fgNAllRows = 0;
Int_t        AliEveTPCSectorData::fgNAllPads = 0;
Int_t*       AliEveTPCSectorData::fgRowBegs  = 0;

AliEveTPCSectorData::SegmentInfo AliEveTPCSectorData::fgInnSeg;
AliEveTPCSectorData::SegmentInfo AliEveTPCSectorData::fgOut1Seg;
AliEveTPCSectorData::SegmentInfo AliEveTPCSectorData::fgOut2Seg;

AliEveTPCSectorData::SegmentInfo* AliEveTPCSectorData::fgSegInfoPtrs[3] = {0};

/******************************************************************************/

void AliEveTPCSectorData::InitStatics()
{
  // Initialize static variables.

  if (fgParam != 0) return;

  fgParam    = new AliTPCParamSR;
  fgZLength  = fgParam->GetZLength(0) + 0.275;
  fgNAllRows = fgParam->GetNRowLow()  + fgParam->GetNRowUp();
  fgNAllPads = 0;
  fgRowBegs  = new Int_t[fgNAllRows + 1];

  {
    Int_t row = 0;
    for (Int_t i=0; i<fgParam->GetNRowLow(); ++i, ++row)
    {
      fgRowBegs[row] = fgNAllPads;
      fgNAllPads += fgParam->GetNPadsLow(i);
    }
    for (Int_t i=0; i<fgParam->GetNRowUp(); ++i, ++row)
    {
      fgRowBegs[row] = fgNAllPads;
      fgNAllPads += fgParam->GetNPadsUp(i);
    }
    fgRowBegs[fgNAllRows] = fgNAllPads;
  }

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

  // Set stepsize arrays

  { // Inner
    Int_t k=0, npads = fgParam->GetNPadsLow(0);
    for (Int_t row = 0; row < fgInnSeg.fNRows; ++row)
    {
      if (fgParam->GetNPadsLow(row) > npads)
      {
        npads = fgParam->GetNPadsLow(row);
        fgInnSeg.fYStep[k] = row*fgInnSeg.fPadHeight + fgInnSeg.fRLow;
        k++;
      }
    }
    fgInnSeg.fNYSteps = k;
  }

  {  // Outer 1 seg
    Int_t k=0, npads = fgParam->GetNPadsUp(0);
    for (Int_t row = 0; row < fgOut1Seg.fNRows; ++row)
    {
      if (fgParam->GetNPadsUp(row) > npads)
      {
        npads = fgParam->GetNPadsUp(row);
        fgOut1Seg.fYStep[k] = row*fgOut1Seg.fPadHeight + fgOut1Seg.fRLow ;
        k++;
      }
    }
    fgOut1Seg.fNYSteps = k;
  }

  {  // Outer 2 seg
    Int_t k=0, npads = fgParam->GetNPadsUp(fgOut1Seg.fNRows);
    for (Int_t row = fgOut1Seg.fNRows; row < fgParam->GetNRowUp(); ++row)
    {
      if (fgParam->GetNPadsUp(row) > npads)
      {
        npads = fgParam->GetNPadsUp(row);
        fgOut2Seg.fYStep[k] = (row - fgOut1Seg.fNRows)*fgOut2Seg.fPadHeight + fgOut2Seg.fRLow ;
        k++;
      }
    }
    fgOut2Seg.fNYSteps = k;
  }
}

Int_t AliEveTPCSectorData::GetNPadsInRow(Int_t row)
{
  // Return number of pads in given row.

  if (row < 0 || row >= fgNAllRows) return 0;
  return fgRowBegs[row + 1] - fgRowBegs[row];
}

const AliEveTPCSectorData::SegmentInfo& AliEveTPCSectorData::GetSeg(Int_t seg)
{
  // Return reference to segment geometry information.
  // 0 ~ inner, 1 ~ middle, 2 ~ outer.

  static const SegmentInfo null;

  if (seg < 0 || seg > 2)
    return null;
  else
    return *fgSegInfoPtrs[seg];
}

/******************************************************************************/
// True member functions start here.
/******************************************************************************/

void AliEveTPCSectorData::NewBlock()
{
  // Create new data-block. Position is set to the beginning.

  fBlocks.push_back(new Short_t[fkBlockSize]);
  fBlockPos = 0;
}

/******************************************************************************/

AliEveTPCSectorData::AliEveTPCSectorData(Int_t sector, Int_t bsize) :
  fSectorID(sector),  fNPadsFilled(0), fPads(),
  fkBlockSize(bsize), fBlockPos(0),    fBlocks(),
  fCurrentRow(0), fCurrentPad(0), fCurrentPos(0), fCurrentStep(0)
{
  // Constructor.
	
  memset(fPadBuffer,0,2048*sizeof(Short_t));

  if (fgParam == 0) InitStatics();

  fPads.assign(fgNAllPads, PadData());
  fBlocks.reserve(16);
  fBlockPos = fkBlockSize; // Enforce creation of a new block.
}


AliEveTPCSectorData::~AliEveTPCSectorData()
{
  // Destructor.

  for (std::vector<Short_t*>::iterator b=fBlocks.begin(); b!=fBlocks.end(); ++b)
    delete [] *b;
}

void AliEveTPCSectorData::DropData()
{
  // Drop data, deallocate data-blocks.

  fPads.assign(fgNAllPads, PadData());
  for (std::vector<Short_t*>::iterator b=fBlocks.begin(); b!=fBlocks.end(); ++b)
    delete [] *b;
  fBlocks.clear();
  fBlockPos = fkBlockSize; // Enforce creation of a new block.
}

/******************************************************************************/

void AliEveTPCSectorData::Print(Option_t* /*opt*/) const
{
  // Print summary information.

  printf("AliEveTPCSectorData sector=%d, NPadsFilled=%d, NBlocks=%d, BlockPos=%d\n",
	 fSectorID, fNPadsFilled, (Int_t) fBlocks.size(), fBlockPos);
}

/******************************************************************************/

void AliEveTPCSectorData::BeginPad(Int_t row, Int_t pad, Bool_t reverseTime)
{
  // Begin filling of pad-data as specified with arguments.

  fCurrentRow = row;
  fCurrentPad = pad;
  if (reverseTime) {
    fCurrentPos  = 2046;
    fCurrentStep = -2;
  } else {
    fCurrentPos  = 0;
    fCurrentStep = 2;
  }
  //printf("begpad for row=%d pad=%d\n  buf=%p pos=%d step=%d\n",
  //     fCurrentRow, fCurrentPad,
  //     fPadBuffer, fCurrentPos, fCurrentStep);
}

void AliEveTPCSectorData::EndPad(Bool_t autoPedestal, Short_t threshold)
{
  // End filling of pad-data. At this point data is compressed and moved
  // into the cuurent position in memory block.

  Short_t *beg, *end;
  if (fCurrentStep > 0) {
    beg = fPadBuffer;
    end = fPadBuffer + fCurrentPos;
  } else {
    beg = fPadBuffer + fCurrentPos + 2;
    end = fPadBuffer + 2048;
  }

  //printf("endpad for row=%d pad=%d\n  buf=%p beg=%p end=%p pos=%d step=%d\n",
  //     fCurrentRow, fCurrentPad,
  //     fPadBuffer, beg, end, fCurrentPos, fCurrentStep);
  if (beg >= end)
    return;

  if (autoPedestal) {
    Short_t array[1024];
    Short_t* val;
    val = beg + 1;
    while (val <= end) {
      array[(val-beg)/2] = *val;
      val += 2;
    }
    Short_t pedestal = TMath::Nint(TMath::Median((end-beg)/2, array));
    val = beg + 1;
    while (val <= end) {
      *val -= pedestal;
      val += 2;
    }
    Short_t* wpos = beg;
    Short_t* rpos = beg;
    while (rpos < end) {
      if (rpos[1] >= threshold) {
	wpos[0] = rpos[0];
	wpos[1] = rpos[1];
	wpos += 2;
      }
      rpos += 2;
    }
    end = wpos;
  }

  Short_t* wpos = beg;
  Short_t* rpos = beg;

  // Compress pad buffer
  while (rpos < end) {
    Short_t* spos = rpos;
    Short_t  t    = spos[0];
    while (true) {
      rpos += 2;
      if (rpos >= end || *rpos > t + 1 || t == 0)
	break;
      ++t;
    }
    Short_t n = t - spos[0] + 1;
    if (n == 1) {
      wpos[0] = -spos[0];
      wpos[1] =  spos[1];
      wpos += 2;
    } else {
      wpos[0] = spos[0];
      wpos[2] = spos[1];
      wpos[1] = n;
      wpos += 3; spos += 3;
      while (--n) {
	*wpos = *spos;
	++wpos;	spos += 2;
      }
    }
  }

  // Copy buffer to storage, set PadData
  if (wpos > beg) {
    Short_t len = wpos - beg;
    if (len > fkBlockSize - fBlockPos)
      NewBlock();
    Short_t *dest = fBlocks.back() + fBlockPos;
    memcpy(dest, beg, len*sizeof(Short_t));
    fBlockPos += len;

    PadData& pad = fPads[PadIndex(fCurrentRow, fCurrentPad)];
    pad.SetDataLength(dest, len);
  }

  ++fNPadsFilled;
}

/******************************************************************************/

const AliEveTPCSectorData::PadData& AliEveTPCSectorData::GetPadData(Int_t padAddr) const
{
  // Get pad-data reference by absolute index.

  static const PadData kNull;

  if (padAddr < 0 || padAddr >= fgNAllPads) return kNull;
  return fPads[padAddr];
}

const AliEveTPCSectorData::PadData& AliEveTPCSectorData::GetPadData(Int_t row, Int_t pad) const
{
  // Get pad-data reference by row and pad number.

  static const PadData kNull;

  Int_t np = GetNPadsInRow(row);
  if (np == 0 || pad < 0 || pad >= np) return kNull;
  return GetPadData(fgRowBegs[row] + pad);
}

AliEveTPCSectorData::PadIterator AliEveTPCSectorData::MakePadIterator(Int_t padAddr, Short_t thr)
{
  // Get pad-data iterator by absolute index.

  return PadIterator(GetPadData(padAddr), thr);
}

AliEveTPCSectorData::PadIterator AliEveTPCSectorData::MakePadIterator(Int_t row, Int_t pad, Short_t thr)
{
  // Get pad-data iterator by row and pad number.

  return PadIterator(GetPadData(row, pad), thr);
}

AliEveTPCSectorData::RowIterator AliEveTPCSectorData::MakeRowIterator(Int_t row, Short_t thr)
{
  // Get row iterator.

  Short_t npads = GetNPadsInRow(row);
  if (npads > 0)
    return RowIterator(&fPads[fgRowBegs[row]], npads, thr);
  else
    return RowIterator(0, 0);
}

/******************************************************************************/
// AliEveTPCSectorData::PadData
/******************************************************************************/

void AliEveTPCSectorData::PadData::Print(Option_t* /*opt*/)
{
  // Print summary information.

  printf("addr=%p, len=%hd>\n", (void*)fData, fLength);
  for (Int_t i=0; i<fLength; ++i)
    printf("  %3d %hd\n", i, fData[i]);
}

/******************************************************************************/
// AliEveTPCSectorData::PadIterator
/******************************************************************************/

Bool_t AliEveTPCSectorData::PadIterator::Next()
{
  // Move iterator to next signal above the iteration threshold.
  // Returns false when the end of data is reached.

  if (fPos >= fEnd) return kFALSE;
  if (fNChunk > 0) {
    ++fTime;
    --fNChunk;
    fSignal = *fPos; ++fPos;
  } else {
    fTime = fPos[0];
    if (fTime <= 0) {
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

void AliEveTPCSectorData::PadIterator::Reset()
{
  // Return to the beginning of the pad-data. Must call Next() to get to
  // the first stored signal.

  fPos    = fBeg;
  fTime   = -1;
  fSignal = -1;
  fNChunk = 0;
}

void AliEveTPCSectorData::PadIterator::Reset(const PadData& pd)
{
  // Reinitialize to new pad-data. Must call Next() to get to the
  // first stored signal.

  fBeg = pd.Data();
  fEnd = pd.Data() + pd.Length();
  fPos = pd.Data();
  Reset();
}

void AliEveTPCSectorData::PadIterator::Test()
{
  while (Next())
    printf("  %3d %d\n", fTime, fSignal);
}

/******************************************************************************/
// AliEveTPCSectorData::RowIterator
/******************************************************************************/

Bool_t AliEveTPCSectorData::RowIterator::NextPad()
{
  // Move iterator to next pad.

  ++fPad;
  if(fPad >= fNPads) return kFALSE;
  Reset(fPadArray[fPad]);
  return kTRUE;
}

void AliEveTPCSectorData::RowIterator::ResetRow()
{
  // Return to the beginning of the row. Must call NextPad() to get to
  // the zeroth pad.

  fPad = -1;
}

void AliEveTPCSectorData::RowIterator::ResetRow(const PadData* first, Short_t npads)
{
  // Reinitialize to another pad-data array. Must call NextPad() to
  // get to the zeroth pad.

  fPadArray =  first;
  fNPads    =  npads;
  fPad      = -1;
}

/******************************************************************************/
// AliEveTPCSectorData::SegmentInfo
/******************************************************************************/

//______________________________________________________________________________
//
// Stores geometry data about a segment needed for fast data-access
// and rendering

ClassImp(AliEveTPCSectorData::SegmentInfo)

AliEveTPCSectorData::SegmentInfo::SegmentInfo() :
  TObject(),

  fPadWidth(0), fPadHeight(0),
  fRLow(0), fNRows(0), fFirstRow(0), fLastRow(0),
  fNMaxPads(0),
  fNYSteps(0)
{
  // Constructor.

  memset(fYStep, 0, sizeof(fYStep));
}
