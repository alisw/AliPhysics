// $Header$

#ifndef ALIEVE_TPCSectorData_H
#define ALIEVE_TPCSectorData_H

#include <Reve/Reve.h>

#include <TObject.h>

#include <vector>

class AliTPCParam;

namespace Alieve {

class TPCSectorData : public TObject
{
  TPCSectorData(const TPCSectorData&);            // Not implemented
  TPCSectorData& operator=(const TPCSectorData&); // Not implemented

public:

  class PadData
  {
  protected:
    Short_t* fData;
    Short_t  fLength;

  public:
    PadData(Short_t* d=0, Short_t l=0) : fData(d), fLength(l) {}

    PadData(const PadData& p) : fData(p.fData), fLength(p.fLength) {}
    PadData& operator=(const PadData& p)
    { fData = p.fData; fLength = p.fLength; return *this; }

    Short_t* Data()   const { return fData; }
    Short_t  Length() const { return fLength; }

    void SetDataLength(Short_t* d, Short_t l) { fData = d; fLength = l; }

    void Print(Option_t* opt="");
  };

  class PadIterator
  {
  protected:
    Short_t *fBeg, *fEnd, *fPos;
    Short_t  fTime, fSignal;
    Short_t  fThreshold;
    Short_t  fNChunk;

  public:
    PadIterator(const PadData& pd, Short_t thr=0) :
      fBeg(pd.Data()), fEnd(pd.Data() + pd.Length()), fPos(pd.Data()),
      fTime(-1), fSignal(-1), fThreshold(thr), fNChunk(0)
    {}
    PadIterator(const PadIterator& i) :
      fBeg(i.fBeg), fEnd(i.fEnd), fPos(i.fPos),
      fTime(i.fTime), fSignal(i.fSignal), fThreshold(i.fThreshold), fNChunk(i.fNChunk)
    {}
    virtual ~PadIterator() {}

    PadIterator& operator=(const PadIterator& i) {
      fBeg = i.fBeg; fEnd = i.fEnd; fPos = i.fPos;
      fTime = i.fTime; fSignal = i.fSignal; fThreshold = i.fThreshold; fNChunk = i.fNChunk;
      return *this;
    }

    Bool_t Next();
    void   Reset();
    void   Reset(const PadData& pd);

    Short_t Time()   const { return fTime; }
    Short_t Signal() const { return fSignal; }

    Short_t Threshold()    const { return fThreshold; }
    void SetThreshold(Short_t t) { fThreshold = t; }

    void Test();
  };

  class RowIterator : public PadIterator
  {
  protected:
    const PadData* fPadArray;
    Short_t fNPads;
    Short_t fPad;

  public:
    RowIterator(const PadData* first, Short_t npads, Short_t thr=0) :
      PadIterator(*first, thr),
      fPadArray(first), fNPads(npads),
      fPad(-1)
    {}
    RowIterator(const RowIterator& i) :
      PadIterator(i),
      fPadArray(i.fPadArray), fNPads(i.fNPads), fPad(i.fPad)
    {}

    RowIterator& operator=(const RowIterator& i) {
      fPadArray = i.fPadArray; fNPads = i.fNPads; fPad = i.fPad;
      return *this;
    }

    Bool_t NextPad();
    void   ResetRow();
    void   ResetRow(const PadData* first, Short_t npads);
    
    Short_t Pad() const { return fPad; }

    void Test();
  };

  class SegmentInfo : public TObject
  {
    friend class TPCSectorData;

  private:
    Float_t   fPadWidth;
    Float_t   fPadHeight;
    Float_t   fRLow;      // Radius at the bottom of first row
    Int_t     fNRows;     // Number of rows in this segment
    Int_t     fFirstRow;  // First row index within sector
    Int_t     fLastRow;   // Last row index within sector
    Int_t     fNMaxPads;  // Maximum number of pads in a row
    Int_t     fNYSteps;   // Number of steps in pad-count
    Float_t   fYStep[64]; // Y coords where pad-count changes

  public:
    SegmentInfo();

    Float_t GetPadWidth()  const { return fPadWidth; }
    Float_t GetPadHeight() const { return fPadHeight; }
    Float_t GetRLow()      const { return fRLow; }
    Int_t   GetNRows()     const { return fNRows; }
    Int_t   GetFirstRow()  const { return fFirstRow; }
    Int_t   GetLastRow()   const { return fLastRow; }
    Int_t   GetNMaxPads()  const { return fNMaxPads; }
    Int_t   GetNYSteps()   const { return fNYSteps; }
    Float_t GetYStep(Int_t step) const { return fYStep[step]; }

    ClassDef(SegmentInfo, 0);
  };

private:
  static AliTPCParam *fgParam;
  static Float_t      fgZLength;
  static Int_t        fgNAllRows;
  static Int_t        fgNAllPads;
  static Int_t       *fgRowBegs;

  static SegmentInfo  fgInnSeg;
  static SegmentInfo  fgOut1Seg;
  static SegmentInfo  fgOut2Seg;

  static SegmentInfo* fgSegInfoPtrs[3];

protected:
  Int_t                 fSectorID;
  Int_t                 fNPadsFilled;
  std::vector<PadData>  fPads;

  // Blocks of pad-data.
  const Int_t           fBlockSize;
  Int_t                 fBlockPos;
  std::vector<Short_t*> fBlocks;

  void NewBlock();

  // Intermediate buffer/vars for pad-data.
  Short_t fPadBuffer[2048];
  Int_t   fCurrentRow, fCurrentPad, fCurrentPos, fCurrentStep;

  Int_t PadIndex(Int_t row, Int_t pad) { return fgRowBegs[row] + pad; }

public:
  TPCSectorData(Int_t sector, Int_t bsize=65536);
  virtual ~TPCSectorData();

  void DropData();

  virtual void Print(Option_t* opt="") const;

  void BeginPad(Int_t row, Int_t pad, Bool_t reverseTime=kFALSE);
  void RegisterData(Short_t time, Short_t signal);
  void EndPad(Bool_t autoPedestal=kFALSE, Short_t threshold=0);

  const PadData& GetPadData(Int_t padAddr);
  const PadData& GetPadData(Int_t row, Int_t pad);

  PadIterator MakePadIterator(Int_t padAddr, Short_t thr=0);
  PadIterator MakePadIterator(Int_t row, Int_t pad, Short_t thr=0);

  RowIterator MakeRowIterator(Int_t row, Short_t thr=0);

  // --- Static functions

  static const AliTPCParam& GetParam() { return *fgParam; }
  static Float_t GetZLength()  { return fgZLength;  }
  static Int_t   GetNAllRows() { return fgNAllRows; }
  static Int_t   GetNAllPads() { return fgNAllPads; }

  static Int_t GetNPadsInRow(Int_t row);

  static const SegmentInfo& GetInnSeg()  { return fgInnSeg;  }
  static const SegmentInfo& GetOut1Seg() { return fgOut1Seg; }
  static const SegmentInfo& GetOut2Seg() { return fgOut2Seg; }

  static const SegmentInfo& GetSeg(Int_t seg);
  
  static void InitStatics();


  //----------------------------------------------------------------
  // Hack for noisy pad-row removal
  //----------------------------------------------------------------

  class PadRowHack
  {
  public:
    Int_t   fRow, fPad;
    Int_t   fThrExt;
    Float_t fThrFac; // Actual threshold = fThrExt + fThrFac*thr

    PadRowHack(Int_t r, Int_t p, Int_t te=0, Float_t tf=1) :
      fRow(r), fPad(p), fThrExt(te), fThrFac(tf) {}
    bool operator<(const PadRowHack& a) const
    { return (fRow == a.fRow) ? fPad < a.fPad : fRow < a.fRow; }
  };

  PadRowHack* GetPadRowHack(Int_t r, Int_t p);
  void AddPadRowHack(Int_t r, Int_t p, Int_t te=0, Float_t tf=1);
  void RemovePadRowHack(Int_t r, Int_t p);
  void DeletePadRowHack();

protected:
  void* fPadRowHackSet;
  

  ClassDef(TPCSectorData, 0);
}; // endclass TPCSectorData


inline void TPCSectorData::RegisterData(Short_t time, Short_t signal)
{
  fPadBuffer[fCurrentPos]   = time;
  fPadBuffer[fCurrentPos+1] = signal;
  fCurrentPos += fCurrentStep;
}

}

#endif
