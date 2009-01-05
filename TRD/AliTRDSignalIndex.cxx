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

#include "TObject.h"

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
  ,fBoolIndex(NULL)
  ,fSortedIndex(NULL)
  ,fMaxLimit(0)
  ,fPositionRC(0)
  ,fSortedWasInit(kFALSE)
  ,fCurrRow(0)
  ,fCurrCol(0)
  ,fCurrTbin(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fNtbins(0)
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
  ,fBoolIndex(NULL)
  ,fSortedIndex(NULL)
  ,fMaxLimit(0)
  ,fPositionRC(0)
  ,fSortedWasInit(kFALSE)
  ,fCurrRow(0)
  ,fCurrCol(0)
  ,fCurrTbin(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fNtbins(0)
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
  ,fBoolIndex(NULL)
  ,fSortedIndex(NULL)
  ,fMaxLimit(a.fMaxLimit)
  ,fPositionRC(a.fPositionRC)
  ,fSortedWasInit(a.fSortedWasInit)
  ,fCurrRow(a.fCurrRow)
  ,fCurrCol(a.fCurrCol)
  ,fCurrTbin(a.fCurrTbin)
  ,fNrows(a.fNrows)
  ,fNcols(a.fNcols)
  ,fNtbins(a.fNtbins)
  ,fHasEntry(a.fHasEntry)
{
  //
  // Copy constructor
  //

  fBoolIndex = new Bool_t[fMaxLimit];
  memcpy(fBoolIndex, a.fBoolIndex, fMaxLimit*sizeof(Bool_t));

  fSortedIndex = new Short_t[2*fMaxLimit];
  memcpy(fSortedIndex, a.fSortedIndex, 2*fMaxLimit*sizeof(Short_t));
}

//_____________________________________________________________________________
AliTRDSignalIndex::~AliTRDSignalIndex()
{
  //
  // Destructor
  //

  if (fBoolIndex) {
    delete [] fBoolIndex;
    fBoolIndex = NULL;
  }

if (fSortedIndex) {
    delete [] fSortedIndex;
    fSortedIndex = NULL;
  }

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
  ((AliTRDSignalIndex &)a).fMaxLimit    = fMaxLimit;
  ((AliTRDSignalIndex &)a).fPositionRC    = fPositionRC;
  ((AliTRDSignalIndex &)a).fSortedWasInit = fSortedWasInit;
  ((AliTRDSignalIndex &)a).fCurrRow       = fCurrRow;
  ((AliTRDSignalIndex &)a).fCurrCol       = fCurrCol;
  ((AliTRDSignalIndex &)a).fCurrTbin      = fCurrTbin;
  ((AliTRDSignalIndex &)a).fNrows         = fNrows;
  ((AliTRDSignalIndex &)a).fNcols         = fNcols;
  ((AliTRDSignalIndex &)a).fNtbins        = fNtbins;
  ((AliTRDSignalIndex &)a).fHasEntry      = fHasEntry;

  if(((AliTRDSignalIndex &)a).fBoolIndex)
    {
      delete [] ((AliTRDSignalIndex &)a).fBoolIndex;
    }
  ((AliTRDSignalIndex &)a).fBoolIndex = new Bool_t[fMaxLimit];
  memcpy(((AliTRDSignalIndex &)a).fBoolIndex, fBoolIndex, fMaxLimit*sizeof(Bool_t));

  if(((AliTRDSignalIndex &)a).fSortedIndex)
    {
      delete [] ((AliTRDSignalIndex &)a).fSortedIndex;
    }
  ((AliTRDSignalIndex &)a).fSortedIndex = new Short_t[2*fMaxLimit];
  memcpy(((AliTRDSignalIndex &)a).fSortedIndex, fSortedIndex, 2*fMaxLimit*sizeof(Short_t));

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
void AliTRDSignalIndex::Allocate(const Int_t nrow, const Int_t ncol, const Int_t ntime)
{
  //
  // Create the arrays
  //

  fNrows = nrow;
  fNcols = ncol;
  fNtbins = ntime;

  fMaxLimit = nrow * ncol + 1;

  if (fBoolIndex) {
    delete [] fBoolIndex;
    fBoolIndex = NULL;
  }
  if (fSortedIndex) {
    delete [] fSortedIndex;
    fSortedIndex = NULL;
  }

  fBoolIndex = new Bool_t[fMaxLimit];
  fSortedIndex = new Short_t[2*fMaxLimit];

  ResetArrays();
 
  ResetCounters();

  fHasEntry = kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::ResetArrays()
{
  memset(fBoolIndex,0x00,sizeof(Bool_t)*fMaxLimit);
  memset(fSortedIndex,0xFF,2*sizeof(Short_t)*fMaxLimit);
  fSortedWasInit = kFALSE;
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

  ResetArrays();
  ResetCounters();

  fHasEntry = kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::ResetContentConditional(const Int_t nrow, const Int_t ncol, const Int_t ntime)
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
    ResetArrays();
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

  if (fBoolIndex) {
    delete [] fBoolIndex;
    fBoolIndex = NULL;
  }

  if (fSortedIndex) {
    delete [] fSortedIndex;
    fSortedIndex = NULL;
  }
  
  ResetCounters();

  fHasEntry = kFALSE;
  fSortedWasInit = kFALSE;
  fMaxLimit = 0;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::AddIndexTBin(Int_t row, Int_t col, Int_t /*tbin*/)
{
  //
  // This function is now obsolate, it will be deleted in future. 
  //
  // Store the index row-column-tbin as an interesting one
  // The RC index is updated to!!!
  // This is to be used in the TRD clusterizer!
  //
 
  AddIndexRC(row, col);

}

//_____________________________________________________________________________
void AliTRDSignalIndex::AddIndexRC(const Int_t row, const Int_t col)
{
  //
  // Store the index row-column as an interesting one
  // The RC index is updated to!!!
  // This is to be used in the TRD clusterizer!
  //

  if (row * col + 1 >= fMaxLimit) {
    AliError(Form("Out-of-limits row * col %d. Limit is: %d"
                 ,row * col
                 ,fMaxLimit));
    return;
  }
 
  fBoolIndex[row*fNcols+col]=kTRUE;

  fHasEntry = kTRUE;

}

//_____________________________________________________________________________
Bool_t  AliTRDSignalIndex::NextRCIndex(Int_t &row, Int_t &col)
{
  //
  // Returns next used RC combination
  //

  if(fSortedIndex[fPositionRC]>-1){
    row = fCurrRow = fSortedIndex[fPositionRC];
    fPositionRC++;
    col = fCurrCol = fSortedIndex[fPositionRC];
    fPositionRC++;
    return kTRUE;
  }
  else {
    if(fSortedWasInit || !fHasEntry)
      {
        ResetCounters();
	row = fCurrRow;
	col = fCurrCol;
	return kFALSE;
      }
    else
      {
	InitSortedIndex();
	return NextRCIndex(row, col);
      }
  }

}

//_____________________________________________________________________________
Bool_t AliTRDSignalIndex::NextRCTbinIndex(Int_t &row, Int_t &col, Int_t &tbin)
{
  //
  // Returns the next tbin, or if there is no next time bin, it returns the
  // next used RC combination.
  //  

  if (NextTbinIndex(tbin)) {
    row = fCurrRow;
    col = fCurrCol;
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
  // Returns the next tbin of the current RC combination
  //
  
  if(fCurrTbin<fNtbins)
    {
      tbin = fCurrTbin++;
      return kTRUE;
    }

  return kFALSE;

}

//_____________________________________________________________________________
void AliTRDSignalIndex::InitSortedIndex()
{
  //
  // Creates the SortedIndex
  //

  fSortedWasInit = kTRUE;
  int pos=0;
  for(int row = 0; row < fNrows; row++)
    for(int col = 0; col < fNcols; col++)
      if(IsBoolIndex(row, col)){
	fSortedIndex[pos] = row;
	pos++;
	fSortedIndex[pos] = col;
	pos++;
      }
}

//_____________________________________________________________________________
void AliTRDSignalIndex::ResetCounters()    
{ 
  //
  // Reset the counters/iterators
  //

  fCurrRow    = -1;
  fCurrCol    = -1;
  fCurrTbin   = -1;
  fPositionRC =  0;
}
