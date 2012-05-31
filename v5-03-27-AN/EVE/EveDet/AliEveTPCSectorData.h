// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSectorData_H
#define AliEveTPCSectorData_H

#include <TObject.h>

#include <vector>

class AliTPCParam;

//------------------------------------------------------------------------------
// AliEveTPCSectorData
//
// Constainer for pad-data of a single TPC sector.
//  Also stores relevant geometry information in static data-members.
//

class AliEveTPCSectorData : public TObject
{
public:
  // --- Inner classes ---

  class PadData
  {
  public:
    PadData(Short_t* d=0, Short_t l=0) : fData(d), fLength(l) {}

    PadData(const PadData& p) : fData(p.fData), fLength(p.fLength) {}
    PadData& operator=(const PadData& p)
      { if(this!=&p){fData = p.fData; fLength = p.fLength;} return *this; }

    Short_t* Data()   const { return fData; }
    Short_t  Length() const { return fLength; }

    void SetDataLength(Short_t* d, Short_t l) { fData = d; fLength = l; }

    void Print(Option_t* opt="");

  protected:
    Short_t* fData;   // Data for given pad.
    Short_t  fLength; // Length of pad-data.
  };

  class PadIterator
  {
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

    PadIterator& operator=(const PadIterator& i) {if(this!=&i){
      fBeg = i.fBeg; fEnd = i.fEnd; fPos = i.fPos;
      fTime = i.fTime; fSignal = i.fSignal; fThreshold = i.fThreshold; fNChunk = i.fNChunk;
      }return *this;
    }

    Bool_t Next();
    void   Reset();
    void   Reset(const PadData& pd);

    Short_t Time()   const { return fTime; }
    Short_t Signal() const { return fSignal; }

    Short_t Threshold()    const { return fThreshold; }
    void SetThreshold(Short_t t) { fThreshold = t; }

    void Test();

  protected:
    Short_t *fBeg, *fEnd;    // Begin and end of data.
    Short_t *fPos;           // Current position.
    Short_t  fTime, fSignal; // Current time and signal.
    Short_t  fThreshold;     // Threshold for data iteration. 
    Short_t  fNChunk;        // Number of contiguous signals still to read.
  };

  class RowIterator : public PadIterator
  {
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

    RowIterator& operator=(const RowIterator& i) {if(this!=&i){
      fPadArray = i.fPadArray; fNPads = i.fNPads; fPad = i.fPad;
      }return *this;
    }

    Bool_t NextPad();
    void   ResetRow();
    void   ResetRow(const PadData* first, Short_t npads);

    Short_t TEvePad() const { return fPad; }

  protected:
    const PadData *fPadArray; // Pointer to array of pad-data.
    Short_t        fNPads;    // Number of pads in row.
    Short_t        fPad;      // Current pad.
  };

  class SegmentInfo : public TObject
  {
    friend class AliEveTPCSectorData;

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

  private:
    Float_t   fPadWidth;  // Width of pad in this segment.
    Float_t   fPadHeight; // Height of pad in this segment.
    Float_t   fRLow;      // Radius at the bottom of first row.
    Int_t     fNRows;     // Number of rows in this segment.
    Int_t     fFirstRow;  // First row index within sector.
    Int_t     fLastRow;   // Last row index within sector.
    Int_t     fNMaxPads;  // Maximum number of pads in a row.
    Int_t     fNYSteps;   // Number of steps in pad-count.
    Float_t   fYStep[64]; // Y coords where pad-count changes.

    ClassDef(SegmentInfo, 0);
  };

  // --- Interface ---

  AliEveTPCSectorData(Int_t sector, Int_t bsize=65536);
  virtual ~AliEveTPCSectorData();

  void DropData();

  virtual void Print(Option_t* opt="") const;

  void BeginPad(Int_t row, Int_t pad, Bool_t reverseTime=kFALSE);
  void RegisterData(Short_t time, Short_t signal);
  void EndPad(Bool_t autoPedestal=kFALSE, Short_t threshold=0);

  const PadData& GetPadData(Int_t padAddr) const;
  const PadData& GetPadData(Int_t row, Int_t pad) const;

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


protected:
  Int_t                 fSectorID;      // Sector id.
  Int_t                 fNPadsFilled;   // Number of filled pads.
  std::vector<PadData>  fPads;          // Vector of pad-data.

  // Blocks of pad-data.
  const Int_t           fkBlockSize;    // Size of pad-data block.
  Int_t                 fBlockPos;      // Position in current block.
  std::vector<Short_t*> fBlocks;        // Vector of blocks.

  void NewBlock();


  // Intermediate buffer/vars used during filling of pad-data.
  Short_t fPadBuffer[2048];     // Buffer for current pad.
  Int_t   fCurrentRow;          // Current row.
  Int_t   fCurrentPad;          // Current pad.
  Int_t   fCurrentPos;          // Current position in pad-buffer.
  Int_t   fCurrentStep;         // Step, can be -2 or +2, depending on fill direction.

  Int_t PadIndex(Int_t row, Int_t pad) const { return fgRowBegs[row] + pad; }


private:
  static AliTPCParam *fgParam;     // Global TPC parameters.
  static Float_t      fgZLength;   // Z-length of a sector.
  static Int_t        fgNAllRows;  // Number of rows in all segments.
  static Int_t        fgNAllPads;  // Number of pads in all segments.
  static Int_t       *fgRowBegs;   // Ids of pads at row-beginnings.

  static SegmentInfo  fgInnSeg;    // Geometry information for inner segment.
  static SegmentInfo  fgOut1Seg;   // Geometry information for middle segment.
  static SegmentInfo  fgOut2Seg;   // Geometry information for outer segment.

  static SegmentInfo* fgSegInfoPtrs[3]; // Array of geometry information objects, for access by segment id.

  AliEveTPCSectorData(const AliEveTPCSectorData&);            // Not implemented
  AliEveTPCSectorData& operator=(const AliEveTPCSectorData&); // Not implemented

  ClassDef(AliEveTPCSectorData, 0); // Holds pad-data of a single TPC sector. Also stores geometry information in static data-members.
};


inline void AliEveTPCSectorData::RegisterData(Short_t time, Short_t signal)
{
  // Register data for given time.

  fPadBuffer[fCurrentPos]   = time;
  fPadBuffer[fCurrentPos+1] = signal;
  fCurrentPos += fCurrentStep;
}

#endif
