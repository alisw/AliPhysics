#ifndef AliTRDSIGNALINDEX_H
#define AliTRDSIGNALINDEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
#include "TObject.h"

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  General container for data from TRD detector segments                 //
//  Adapted from AliDigits, origin M.Ivanov                               //
//                                                                        //
//  Author:                                                               //
//    Mateusz Ploskon (ploskon@ikf.uni-frankfurt.de)                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TArrayI.h"

class AliTRDSignalIndex : public TObject
{

 public:

  AliTRDSignalIndex(); 
  AliTRDSignalIndex(Int_t nrow, Int_t ncol,Int_t ntime);
  AliTRDSignalIndex(const AliTRDSignalIndex &d);
  virtual ~AliTRDSignalIndex();
  AliTRDSignalIndex &operator=(const AliTRDSignalIndex &d); 

  virtual void     Copy(TObject &d) const; 
  virtual void     Allocate(Int_t nrow, Int_t ncol,Int_t ntime);

  virtual void     Reset();
  virtual void     ResetContentConditional(Int_t nrow, Int_t ncol,Int_t ntime);
  virtual void     ResetContent();
  virtual void     ResetCounters();
  virtual void     ResetTbinCounter() { fPositionTbin  = 1; }

  virtual void     AddIndexTBin(Int_t row, Int_t col, Int_t tbin);

          // Get the next pad (row and column) and return kTRUE on success
          Bool_t   NextRCIndex(Int_t &row, Int_t &col); 
          // Get the next timebin of a pad (row and column) and return kTRUE on success
          Bool_t   NextRCTbinIndex(Int_t &row, Int_t &col, Int_t &tbin); 
          // Get the next active timebin and return kTRUE on success
          Bool_t   NextTbinIndex(Int_t &tbin); 

          Int_t    GetCurrentRow() const  { return (*fIndex)[fPositionRow];                 }
          Int_t    GetCurrentCol() const  { return (*fIndex)[fPositionCol];                 }
          Int_t    GetCurrentTbin() const { return (*fIndex)[fPositionCol + fPositionTbin]; }
  
          // Clear the array, actually destroy and recreate w/o allocating
          void     ClearAll(); 
          // Return kTRUE if array allocated and there is no need to call allocate
          Bool_t   IsAllocated() const    { if (!fIndex)                return kFALSE; 
                                            if (fIndex->GetSize() <= 0) return kFALSE; 
                                            else                        return kTRUE;       }

          void     SetSM(Int_t ix)        { fSM      =    ix; }
          void     SetStack(Int_t ix)     { fStack   =    ix; }
          void     SetLayer(Int_t ix)     { fLayer   =    ix; }
          void     SetDetNumber(Int_t ix) { fDet     =    ix; }
  
  virtual Int_t    GetDetNumber() const   { return fDet;      } // Get Det number
  virtual Int_t    GetLayer() const       { return fLayer;    } // Layer position of the chamber in TRD
  virtual Int_t    GetStack() const       { return fStack;    } // Stack position of the chamber in TRD
  virtual Int_t    GetSM() const          { return fSM;       } // Super module of the TRD
          TArrayI *GetArray() const       { return fIndex;    } // Get the tarrayi pointer for god knows what reason

  virtual Bool_t   HasEntry() const       { return fHasEntry; } // Return status if has an entry

  virtual Int_t    GetNrow() const        { return fNrows;    } // Get Nrows
  virtual Int_t    GetNcol() const        { return fNcols;    } // Get Ncols
  virtual Int_t    GetNtime() const       { return fNtbins;   } // Get Ntbins

 private:

  Int_t     fDet;                //  Detector number
  Int_t     fLayer;              //  Layer position in the full TRD
  Int_t     fStack;              //  Stack position in the full TRD
  Int_t     fSM;                 //  Super module - position in the full TRD

  TArrayI  *fIndex;              //! Monitor active pads and tbins

  Int_t     fPositionRow;        //  Position in the index - jumps by 1 + 1 + fNtbins
  Int_t     fPositionCol;        //  Position in the index - jumps by 1 + 1 + fNtbins
  Int_t     fPositionTbin;       //  Position in the tbin - goes from 0 to fNtbins

  Int_t     fLastRow;            //  To keep track what is the RC combination
  Int_t     fLastCol;            //  To keep track what is the RC combination
  Int_t     fLastTbin;           //  To keep track what is the Tbin - will catch if raw data bogus

  Int_t     fNrows;              //  Number of rows in the chamber
  Int_t     fNcols;              //  Number of cols in the chamber
  Int_t     fNtbins;             //  Number of tbins in the chamber

  Int_t     fMaxLimit;           //  Max number of things in the array  = nrow*ncol*ntime + nrow*ncol*2
  Bool_t    fResetCounters;      //  Reset counter status
  Bool_t    fHasEntry;           //  kTRUE flag if we have an entry 

  ClassDef(AliTRDSignalIndex,2)  //  Data container for one TRD detector segment

};
#endif
