/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  General container for data from TRD detector segments                    //
//  Adapted from AliDigits, origin M.Ivanov                                  //
//                                                                           //
//  Author:                                                                  //
//    Mateusz Ploskon (ploskon@ikf.uni-frankfurt.de)                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TArrayI.h"

#include "AliLog.h"

#include "AliTRDSignalIndex.h"

ClassImp(AliTRDSignalIndex)

//_____________________________________________________________________________
AliTRDSignalIndex::AliTRDSignalIndex()
  :TObject()
  ,fDet(-1)
  ,fLayer(-1)
  ,fStack(-1)
  ,fSM(-1)
  ,fIndex(NULL)
  ,fPositionRow(0)
  ,fPositionCol(0)
  ,fPositionTbin(0)
  ,fLastRow(0)
  ,fLastCol(0)
  ,fLastTbin(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fNtbins(0)
  ,fMaxLimit(0)
  ,fResetCounters(kTRUE)
  ,fHasEntry(kFALSE)
{
  //
  // Default contructor
  //

  ResetCounters();

}

//_____________________________________________________________________________
AliTRDSignalIndex::AliTRDSignalIndex(Int_t nrow, Int_t ncol,Int_t ntime)
  :TObject()
  ,fDet(-1)
  ,fLayer(-1)
  ,fStack(-1)
  ,fSM(-1)
  ,fIndex(NULL)
  ,fPositionRow(0)
  ,fPositionCol(0)
  ,fPositionTbin(0)
  ,fLastRow(0)
  ,fLastCol(0)
  ,fLastTbin(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fNtbins(0)
  ,fMaxLimit(0)
  ,fResetCounters(kTRUE)
  ,fHasEntry(kFALSE)
{
  //
  // Not the default contructor... hmmm...
  //

  Allocate(nrow, ncol, ntime);  

}

//_____________________________________________________________________________
AliTRDSignalIndex::AliTRDSignalIndex(const AliTRDSignalIndex &a)
  :TObject(a)
  ,fDet(a.fDet)
  ,fLayer(a.fLayer)
  ,fStack(a.fStack)
  ,fSM(a.fSM)
  ,fIndex(a.fIndex)
  ,fPositionRow(a.fPositionRow)
  ,fPositionCol(a.fPositionCol)
  ,fPositionTbin(a.fPositionTbin)
  ,fLastRow(a.fLastRow)
  ,fLastCol(a.fLastCol)
  ,fLastTbin(a.fLastTbin)
  ,fNrows(a.fNrows)
  ,fNcols(a.fNcols)
  ,fNtbins(a.fNtbins)
  ,fMaxLimit(a.fMaxLimit)
  ,fResetCounters(a.fResetCounters)
  ,fHasEntry(a.fHasEntry)
{
  //
  // Copy constructor
  //

}

//_____________________________________________________________________________
AliTRDSignalIndex::~AliTRDSignalIndex()
{
  //
  // Destructor
  //

  delete fIndex;
  fIndex = NULL;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::Copy(TObject &a) const
{
  //
  // Copy function
  //

  ((AliTRDSignalIndex &)a).fDet           = fDet;
  ((AliTRDSignalIndex &)a).fLayer         = fLayer;
  ((AliTRDSignalIndex &)a).fStack         = fStack;
  ((AliTRDSignalIndex &)a).fSM            = fSM;
  ((AliTRDSignalIndex &)a).fIndex         = fIndex;
  ((AliTRDSignalIndex &)a).fPositionRow   = fPositionRow;
  ((AliTRDSignalIndex &)a).fPositionTbin  = fPositionTbin;
  ((AliTRDSignalIndex &)a).fLastRow       = fLastRow;
  ((AliTRDSignalIndex &)a).fLastCol       = fLastCol;
  ((AliTRDSignalIndex &)a).fLastTbin      = fLastTbin;
  ((AliTRDSignalIndex &)a).fNrows         = fNrows;
  ((AliTRDSignalIndex &)a).fNcols         = fNcols;
  ((AliTRDSignalIndex &)a).fNtbins        = fNtbins;
  ((AliTRDSignalIndex &)a).fMaxLimit      = fMaxLimit;
  ((AliTRDSignalIndex &)a).fResetCounters = fResetCounters;
  ((AliTRDSignalIndex &)a).fHasEntry      = fHasEntry;

}

//_____________________________________________________________________________
AliTRDSignalIndex& AliTRDSignalIndex::operator = (const AliTRDSignalIndex& a)
{
  //
  // Assignment operator
  //

  if (this != &a) ((AliTRDSignalIndex &) a).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::Allocate(Int_t nrow, Int_t ncol,Int_t ntime)
{
  //
  // Create the arrays
  //

  fNrows = nrow;
  fNcols = ncol;
  fNtbins = ntime;

  fMaxLimit = nrow * ncol * ntime + nrow * ncol * 2;
  if (fIndex) {
    delete fIndex;
    fIndex = NULL;
  }

  fIndex = new TArrayI(fMaxLimit);
  fIndex->Reset(-1);

  ResetCounters();

  fHasEntry = kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::Reset()
{
  //
  // Reset the array but keep the size - realloc
  //

  fDet   = -1;
  fLayer = -1;
  fStack = -1;
  fSM    = -1;

  // All will be lost 
  Allocate(fNrows, fNcols, fNtbins);

}

//_____________________________________________________________________________
void AliTRDSignalIndex::ResetContent()
{
  //
  // Reset the array but keep the size - no realloc
  //

  fIndex->Reset(-1);
  ResetCounters();

  fHasEntry = kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::ResetContentConditional(Int_t nrow, Int_t ncol,Int_t ntime)
{
  //
  // Reset the array but keep the size if no need to enlarge - no realloc
  //

  fDet   = -1;
  fLayer = -1;
  fStack = -1;
  fSM    = -1;

  if ((nrow  > fNrows) || 
      (ncol  > fNcols) || 
      (ntime > fNtbins)) {
    Allocate(nrow, ncol, ntime);
  }
  else {
    fIndex->Reset(-1);
    ResetCounters();
    fHasEntry = kFALSE;
  }

}

//_____________________________________________________________________________
void AliTRDSignalIndex::ClearAll()
{
  //
  // Reset the values - clear all!
  //

  fDet    = -1;
  fLayer  = -1;
  fStack  = -1;
  fSM     = -1;

  fNrows  = -1;
  fNcols  = -1;
  fNtbins = -1;

  if (fIndex) {
    delete fIndex;
    fIndex = NULL;
  }
  fIndex = new TArrayI();
  ResetCounters();

  fHasEntry = kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::AddIndexTBin(Int_t row, Int_t col, Int_t tbin)
{
  //
  // Store the index row-column-tbin as an interesting one
  // The RC index is updated to!!!
  // This is to be used in the TRD clusterizer!
  //

  if (row * col * tbin + row * col * 2 >= fMaxLimit) {
    AliError(Form("Out-of-limits fPositionCol + fNtbins %d. Limit is: %d"
                 ,fPositionCol + fNtbins
                 ,fMaxLimit));
    return;
  }

  if ((row != fLastRow) || 
      (col != fLastCol)) {

    // New RC combination
    if (fResetCounters == kTRUE) {
      fPositionRow    = 0;
      fPositionCol   = 1;  
      fResetCounters = kFALSE;
    }
    else {
      fPositionRow += fNtbins + 2;
      fPositionCol += fNtbins + 2;
    }

    fPositionTbin = 1;

    (*fIndex)[fPositionRow] = row;
    (*fIndex)[fPositionCol] = col;
    (*fIndex)[fPositionCol + fPositionTbin] = tbin;
    ++fPositionTbin;

  }
  else {

    // Same RCT combination ?
      
    (*fIndex)[fPositionCol + fPositionTbin] = tbin;
    ++fPositionTbin;      
  
  }
  
  fLastRow  = row;
  fLastCol  = col;
  fLastTbin = tbin;

  fHasEntry = kTRUE;

}

//_____________________________________________________________________________
Bool_t  AliTRDSignalIndex::NextRCIndex(Int_t &row, Int_t &col)
{
  //
  // Return the position (index in the data array) of the next available pad
  //

  if (fResetCounters == kTRUE) {
    fPositionRow   = 0;
    fPositionCol   = 1;
    fResetCounters = kFALSE;
  }
  else {
    fPositionRow  += fNtbins + 2;
    fPositionCol  += fNtbins + 2;
  }

  if (fPositionRow >= fMaxLimit) {
    return kFALSE;
  }

  fPositionTbin = 1;

  row = (*fIndex)[fPositionRow];
  col = (*fIndex)[fPositionCol];

  if ((row > -1) && 
      (col > -1)) { 
    return kTRUE;
  }
  
  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDSignalIndex::NextRCTbinIndex(Int_t &row, Int_t &col, Int_t &tbin)
{
  //
  // Return the position (index in the data array) of the next available tbin 
  // within the current pad
  //

  if (fPositionRow >= fMaxLimit) {
    return kFALSE;
  }

  if (NextTbinIndex(tbin)) {
    row  = (*fIndex)[fPositionRow];
    col  = (*fIndex)[fPositionCol];
    fResetCounters = kFALSE;
    return kTRUE;
  }
  else {
    if (NextRCIndex(row, col)) {
      return NextRCTbinIndex(row, col, tbin);
    }
  }

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDSignalIndex::NextTbinIndex(Int_t &tbin)
{
  //
  // Return the position (index in the data array) of the next available tbin 
  // within the current pad
  //

  if ((fPositionCol + fPositionTbin >= fMaxLimit) || 
      (fPositionTbin                >  fNtbins  )) {
    return kFALSE;
  }

  tbin = (*fIndex)[fPositionCol + fPositionTbin];

  if (tbin > -1) {
    ++fPositionTbin;
    return kTRUE;
  }

  return kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::ResetCounters()    
{ 
  //
  // Reset the counters/iterators
  //

  fPositionRow   =  0;
  fPositionCol   =  fPositionRow + 1;
  fPositionTbin  =  1;
  fLastRow       = -1;
  fLastCol       = -1;
  fLastTbin      = -1;
  fResetCounters = kTRUE; 

}
