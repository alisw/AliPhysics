#ifndef AliTRDSIGNALINDEX_H
#define AliTRDSIGNALINDEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
#include "TObject.h"

/////////////////////////////////////////////////////////////
//  General container for data from TRD detector segments  //
//  Adapted from AliDigits, origin M.Ivanov                //
/////////////////////////////////////////////////////////////

//class TArrayI;
#include "TArrayI.h"

class AliTRDSignalIndex : public TObject
{
 public:

  AliTRDSignalIndex(); 
  AliTRDSignalIndex(Int_t nrow, Int_t ncol,Int_t ntime);
  AliTRDSignalIndex(const AliTRDSignalIndex &d);
  virtual ~AliTRDSignalIndex(); // destructor
  AliTRDSignalIndex &operator=(const AliTRDSignalIndex &d); 
  virtual void   Copy(TObject &d) const; 
  virtual void   Allocate(Int_t nrow, Int_t ncol,Int_t ntime);
  virtual void   Reset();

  virtual void   ResetCounters()
  {
    // reset the counters/iterators
    fPositionRow = 0;
    fPositionCol = fPositionRow + 1;
    fPositionTbin = 1;
    
    fLastRow = -1;
    fLastCol = -1;
    fLastTbin = -1;

    fResetCounters = kTRUE;
  }

  virtual void   ResetTbinCounter()
  {
    // reset the time bin counter
    
    fPositionTbin = 1;
  }

  virtual void   AddIndexTBin(Int_t row, Int_t col, Int_t tbin);

  Bool_t  NextRCIndex(Int_t &row, Int_t &col); // get the next pad (row and column) and return kTRUE on success
  Bool_t  NextRCTbinIndex(Int_t &row, Int_t &col, Int_t &tbin); // get the next timebin of a pad (row and column) and return kTRUE on success
  Bool_t  NextTbinIndex(Int_t &tbin); // get the next active timebin and return kTRUE on success

  Int_t   GetCurrentRow() {return (*fIndex)[fPositionRow];} // current row
  Int_t   GetCurrentCol() {return (*fIndex)[fPositionCol];} // current col
  Int_t   GetCurtentTbin() {return (*fIndex)[fPositionCol + fPositionTbin];} //current tbin
  
  void    ClearAll(); // clear the array, actually destroy and recreate w/o allocating
  
  //void    Dump(); // printf content - one way of iterating demo
  //void    Dump2(); // printf content - another way of iterating demo

  //Bool_t  IsAllocated() const {if (fIndex) return kTRUE; else return kFALSE;}
  Bool_t  IsAllocated() const 
  {
    // return kTRUE if array allocated and there is no need to call allocate
    if (!fIndex) 
      return kFALSE; 
    if (fIndex->GetSize() <= 0) 
      return kFALSE; 
    else return kTRUE;
  }

  void SetSM(Int_t ix) 
  { 
    // Set which SM
    fSM = ix;
  };

  void SetStack(Int_t ix) 
  { 
    // Set which stack
    fStack = ix;
  };

  void SetChamber(Int_t ix) 
  {
    // aka set stack
    SetStack(ix);
  }

  void SetLayer(Int_t ix) 
  { 
    // Set which layer
    fLayer = ix;
  };

  void SetPlane(Int_t ix) 
  {
    // aka set plane
    SetLayer(ix);
  }

  void SetDetNumber(Int_t ix)
  {
    // Set Det Number
    fDet = ix;
  }
  
  const Int_t GetDetNumber()     {return fDet;} // Get Det number
  const Int_t GetLayer()   {return fLayer;} // Layer = Plane = position of the chamber in TRD
  const Int_t GetPlane()   {return fLayer;} // Layer = Plane = position of the chamber in TRD
  const Int_t GetStack()   {return fStack;} // Stack = Chameber = position of the chamber in TRD
  const Int_t GetChamber() {return fStack;} // Stack = Chameber = position of the chamber in TRD
  const Int_t GetSM()      {return fSM;} // Super module of the TRD

  const Bool_t HasEntry() {return fHasEntry;} // Return status if has an entry

  TArrayI *GetArray() {return fIndex;} // Get the tarrayi pointer for god knows what reason
  
 private:

  Int_t   fDet;  // det number
  Int_t   fLayer; // aka plane - position in the full TRD
  Int_t   fStack; // aka chamber - position in the full TRD
  Int_t   fSM; // super module - position in the full TRD

 protected:

  TArrayI  *fIndex; //! monitor active pads and tbins

  Int_t    fPositionRow; // position in the index - jumps by 1 + 1 + fNtbins
  Int_t    fPositionCol; // position in the index - jumps by 1 + 1 + fNtbins
  Int_t    fPositionTbin; // position in the tbin - goes from 0 to fNtbins

  Int_t    fLastRow; // to keep track what is the RC combination
  Int_t    fLastCol; // to keep track what is the RC combination
  Int_t    fLastTbin; // to keep track what is the Tbin - will catch if raw data bogus

  Int_t    fNrows; // number of rows in the chamber
  Int_t    fNcols; // number of cols in the chamber
  Int_t    fNtbins; // number of tbins in the chamber

  Int_t    fMaxLimit; // max number of things in the array  = nrow * ncol * ntime + nrow * ncol * 2

  Bool_t   fResetCounters; // reset counter status

  Bool_t   fHasEntry; // kTRUE flag if we have an entry 

  ClassDef(AliTRDSignalIndex,1)      // Data container for one TRD detector segment

};
 
#endif
