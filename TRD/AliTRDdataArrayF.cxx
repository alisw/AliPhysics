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
//  General container for integer data of a TRD detector segment.            //
//  Adapted from AliDigits (origin: M.Ivanov).                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliTRDdataArrayF.h"
#include "AliTRDarrayI.h"
#include "AliTRDarrayF.h"

ClassImp(AliTRDdataArrayF)

//_____________________________________________________________________________
AliTRDdataArrayF::AliTRDdataArrayF()
   :AliTRDdataArray()
   ,fElements(0)
   ,fThreshold(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDdataArrayF::AliTRDdataArrayF(Int_t nrow, Int_t ncol, Int_t ntime)
   :AliTRDdataArray(nrow,ncol,ntime)
   ,fElements(0)
   ,fThreshold(0)
{
  //
  // Creates a AliTRDdataArrayF with the dimensions <nrow>, <ncol>, and <ntime>.
  // The row- and column dimensions are compressible.
  //

  Allocate(nrow,ncol,ntime);
  
}

//_____________________________________________________________________________
AliTRDdataArrayF::AliTRDdataArrayF(const AliTRDdataArrayF &a)
   :AliTRDdataArray(a)
   ,fElements(a.fElements)
   ,fThreshold(a.fThreshold)
{
  //
  // AliTRDdataArrayF copy constructor
  //

}

//_____________________________________________________________________________
AliTRDdataArrayF::~AliTRDdataArrayF()
{
  //
  // Destructor
  //

  if (fElements) {
    delete fElements;
    fElements = 0;
  }
  
}

//_____________________________________________________________________________
void AliTRDdataArrayF::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocates memory for a AliTRDdataArrayF with the dimensions 
  // <nrow>, <ncol>, and <ntime>.
  // The row- and column dimensions are compressible.
  //

  if (fNelems < 0) {
    AliTRDdataArray::Allocate(nrow,ncol,ntime);
  }

  if (fElements) {
    delete fElements;
  }

  fElements = new AliTRDarrayF();
  fElements->Set(fNelems);

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Copy(TObject &a) const
{
  //
  // Copy function
  //

  fElements->Copy(*((AliTRDdataArrayF &) a).fElements);

  ((AliTRDdataArrayF &) a).fThreshold = fThreshold;

  AliTRDdataArray::Copy(a);

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Reset() 
{ 
  //
  // Reset the array (old content gets deleted)
  //
  
  if (fElements) {
    delete fElements;
  }

  fElements = new AliTRDarrayF();
  fElements->Set(0); 

  AliTRDdataArray::Reset();

}

//_____________________________________________________________________________
Int_t AliTRDdataArrayF::GetSize() const
{
  //
  // Returns the size of the complete object
  //

  Int_t size = sizeof(this);

  if (fIndex)    size += sizeof(fIndex)    
                         + fIndex->GetSize()    * sizeof(Int_t);
  if (fElements) size += sizeof(fElements) 
                         + fElements->GetSize() * sizeof(Float_t);

  return size;

}

//_____________________________________________________________________________
Int_t AliTRDdataArrayF::GetDataSize() const 
{
  //
  // Returns the size of only the data part
  //

  if (fElements == 0) 
    return 0;
  else 
    return sizeof(fElements) + fElements->GetSize() * sizeof(Float_t);

}

//_____________________________________________________________________________
Int_t AliTRDdataArrayF::GetOverThreshold(Float_t threshold) 
{
  //
  // Returns the number of entries over threshold
  //
 
  if ((fElements == 0) || (fElements->GetSize() <= 0))
    return 0;
 
  Int_t over = 0;

  for (Bool_t cont = First(); cont == kTRUE; cont = Next()) {
    if ((fCurrentIdx1 < 0) || (fCurrentIdx1 >= fNdim1)) continue;
    if ((fCurrentIdx2 < 0) || (fCurrentIdx2 >= fNdim2)) continue;
    if (fElements->At(fCurrentIndex) > threshold) over++;
  }

  return over;

}

//_____________________________________________________________________________
Float_t AliTRDdataArrayF::GetData(Int_t row, Int_t col, Int_t time) const
{
  //
  // Returns the data value at a given position of the array
  // Includes boundary checking
  //

  if ((row >= 0) && (col >= 0) && (time >= 0)) {
    Int_t idx1 = GetIdx1(row,col);
    if ((idx1 >= 0) && (time < fNdim2)) {
      if (fBufType == 0) return GetDataFast(idx1,time);
      if (fBufType == 1) return GetData1(idx1,time);
    }
    else {
      if (idx1 >= 0) {
        AliError(Form("time %d out of bounds (size: %d, this: 0x%08x)"
                     ,time,fNdim2,this));
      }
    }
  }

  return -1;

}

//_____________________________________________________________________________
Float_t AliTRDdataArrayF::GetDataFast(Int_t idx1, Int_t idx2) const
{
  //
  // Returns the data value at a given position of the array
  // No boundary checking
  //

  return fElements->At(fIndex->At(idx2)+idx1);

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Compress(Int_t bufferType, Float_t threshold)
{
  //
  // Compresses the buffer
  //

  fThreshold = threshold;
  Compress(bufferType);

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Compress(Int_t bufferType)
{
  //
  // Compresses the buffer
  //

  if (fBufType  < 0) {
    AliError("Buffer does not exist");
    return;
  }
  if (fBufType == bufferType) {
    return;
  }  
  if (fBufType > 0) {
    Expand();
  }
  if (fBufType !=0)  {
    AliError("Buffer does not exist");
    return;
  }

  // Compress a buffer of type 1
  if (bufferType == 1) {
    Compress1();
  }

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Expand()
{  
  //
  // Expands the compressed buffer
  //

  if (fBufType  < 0) {
    AliError("Buffer does not exist");
    return;
  }
  if (fBufType == 0) {  
    return;
  } 
 
  // Expand a buffer of type 1
  if (fBufType == 1) Expand1();
  
  fBufType = 0;

}

//_____________________________________________________________________________
Bool_t AliTRDdataArrayF::First() 
{
  //
  // Returns the position of the first valid data value
  //

  if (fBufType == 0) return First0();
  if (fBufType == 1) return First1();
  return kFALSE;

}

//_____________________________________________________________________________
Bool_t  AliTRDdataArrayF::Next()
{
  //
  // Returns the position of the next valid data value
  //

  if (fBufType == 0) return Next0();
  if (fBufType == 1) return Next1();
  return kFALSE;

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Expand1()
{
  //
  // Expands a buffer of type 1
  //

  Int_t i, k;

  fNelems = fNdim1 * fNdim2;

  Float_t *buf = new Float_t[fNelems];
  memset(buf,0,fNelems*sizeof(Float_t)); 

  fIndex->Set(fNdim2);

  for (i = 0, k = 0; i < fNdim2; i++, k += fNdim1) (*fIndex)[i] = k;

  Int_t idx1 = 0;
  Int_t idx2 = 0;
  Int_t n    = fElements->fN;

  for (i = 0; i < n; i++){

    // Negative sign counts the unwritten values (under threshold)
    if ((*fElements)[i] < 0) {
      idx1 -= TMath::Nint(fElements->At(i));
    } 
    else {
      buf[(*fIndex)[idx2] + idx1] = (*fElements)[i];
      idx1++;
    }
    if (idx1 == fNdim1) {
      idx1 = 0;
      idx2++;
    }
    else { 
      if (idx1 > fNdim1){
	Reset();
	return;
      }      
    }
  }

  fElements->Adopt(fNelems,buf); 
   
}

//_____________________________________________________________________________
void AliTRDdataArrayF::Compress1()
{
  //
  // Compress a buffer of type 1
  //

  AliTRDarrayF *buf   = new AliTRDarrayF();  
  buf->Set(fNelems);
  AliTRDarrayI *index = new AliTRDarrayI();
  index->Set(fNdim2);

  Int_t icurrent = -1;
  Int_t izero;
  for (Int_t idx2 = 0; idx2 < fNdim2; idx2++){      

    // Set the idx2 pointer
    (*index)[idx2] = icurrent + 1;

    // Reset the zero counter 
    izero = 0;  

    for (Int_t idx1 = 0; idx1 < fNdim1; idx1++){
      // If below threshold
      if (GetDataFast(idx1,idx2) <= fThreshold) {
        izero++;
      }
      else {
	if (izero > 0) {
	  // If we have currently izero counts under threshold
	  icurrent++;	  
	  if (icurrent >= buf->fN) buf->Expand(icurrent*2);
          // Store the number of entries below zero
	  (*buf)[icurrent] = -izero;  
	  izero = 0;
	} 
	icurrent++;
	if (icurrent >= buf->fN) buf->Expand(icurrent*2);
	(*buf)[icurrent] = GetDataFast(idx1,idx2);	    
      } // If signal larger than threshold	  	
    } // End of loop over idx1

    if (izero > 0) {
      icurrent++;	  
      if (icurrent >= buf->fN) buf->Expand(icurrent*2);
      // Store the number of entries below zero
      (*buf)[icurrent] = -izero;  
    }

  }

  buf->Expand(icurrent+1);
  if (fElements) delete fElements;
  fElements = buf;
  fNelems   = fElements->fN;
  fBufType  = 1;
  if (fIndex) delete fIndex;
  fIndex    = index;

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Expand2()
{
  //
  // Expands a buffer of type 2 
  //

  Int_t i, k;

  Float_t *buf = new Float_t[fNelems];
  memset(buf,0,fNelems*sizeof(Float_t)); 

  fNelems = fNdim1 * fNdim2;
  fIndex->Set(fNdim2);

  for (i = 0, k = 0; i < fNdim2; i++, k += fNdim1) (*fIndex)[i] = k;

  Int_t idx1 = 0;
  Int_t idx2 = 0;
  Int_t n    = fElements->fN;
  for (i = 0; i < n; i++){
    // Negative sign counts the unwritten values (under threshold)
    if ((*fElements)[i] < 0) {
      //idx1 -= (Int_t) fElements->At(i); 
      idx1 -= TMath::Nint(fElements->At(i)); 
    }
    else {
      buf[(*fIndex)[idx2]+idx1] = fElements->At(i);
      idx1++;
    }
    if (idx1 == fNdim1) {
      idx1 = 0;
      idx2++;
    }
    else { 
      if (idx1 > fNdim1){
	Reset();
	return;
      }      
    }
  }

  fElements->Adopt(fNelems,buf);    

}

//_____________________________________________________________________________
void AliTRDdataArrayF::Compress2()
{
  //
  // Compress a buffer of type 2 - not implemented!
  //

}

//_____________________________________________________________________________
Bool_t AliTRDdataArrayF::First0() 
{
  //
  // Returns the first entry for a buffer of type 0
  //

  fCurrentIdx1  = -1;
  fCurrentIdx2  = -1;
  fCurrentIndex = -1;

  Int_t i;
  for (i = 0; ((i < fNelems) && (fElements->At(i) <= fThreshold)); i++)
  if (i == fNelems) return kFALSE;

  fCurrentIdx1  = i % fNdim1;
  fCurrentIdx2  = i / fNdim1;
  fCurrentIndex = i;
  return kTRUE;	

}

//_____________________________________________________________________________
Bool_t AliTRDdataArrayF::Next0()
{
  //
  // Returns the next entry for a buffer of type 0
  //

  if (fCurrentIndex < 0) return kFALSE; 

  Int_t i;
  for (i = fCurrentIndex + 1; 
       ((i < fNelems) && (fElements->At(i) <= fThreshold)); 
       i++) {}
  if (i >= fNelems)  {
    fCurrentIndex = -1;
    return kFALSE;
  }

  fCurrentIdx1  = i % fNdim1;
  fCurrentIdx2  = i / fNdim1;
  fCurrentIndex = i;
  return kTRUE;	

}

//_____________________________________________________________________________
Bool_t AliTRDdataArrayF::First1() 
{
  //
  // Returns the first entry for a buffer of type 1
  //

  fCurrentIdx1  = -1;
  fCurrentIdx2  =  0;
  fCurrentIndex = -1;

  Int_t i;
  for (i = 0; i < fNelems; i++){
    if (fElements->At(i) < 0) {
      //fCurrentIdx1 -= (Int_t) fElements->At(i);
      fCurrentIdx1 -= TMath::Nint(fElements->At(i));
    }
    else {     
      fCurrentIdx1++;
    }
    if (fCurrentIdx1 >= fNdim1) {
      fCurrentIdx2++;
      fCurrentIdx1 -= fNdim1;
    }
    if (fElements->At(i) > fThreshold) break;
  }

  fCurrentIndex = i;
  if (fCurrentIndex >= 0) return kTRUE;
  fCurrentIdx1  = -1;
  fCurrentIdx2  = -1;
  return kFALSE;	

}

//_____________________________________________________________________________
Bool_t AliTRDdataArrayF::Next1() 
{
  //
  // Returns the next entry for a buffer of type 1
  //

  if (fCurrentIndex < 0) return kFALSE;

  Int_t i;
  for (i = fCurrentIndex + 1; i < fNelems; i++){
    if (fElements->At(i) < 0) {
      //fCurrentIdx1 -= (Int_t) fElements->At(i);
      fCurrentIdx1 -= TMath::Nint(fElements->At(i));
    }
    else {      
      fCurrentIdx1++;
    }
    if (fCurrentIdx1 >= fNdim1) {
      fCurrentIdx2++;
      fCurrentIdx1 -= fNdim1;
    }
    if (fElements->At(i) > fThreshold) break;
  }

  fCurrentIndex =  i;
  if ((i >= 0) && (i < fNelems)) return kTRUE;
  fCurrentIdx1  = -1;
  fCurrentIdx2  = -1;
  return kFALSE;

}

//_____________________________________________________________________________
Float_t AliTRDdataArrayF::GetData1(Int_t idx1, Int_t idx2) const
{
  //
  // Returns the value at a given position of the array
  //
  
  Int_t i, n2;

  if ((idx2 + 1) >= fNdim2) {
    n2 = fNelems;
  }
  else {
    n2 = fIndex->At(idx2 + 1);
  }

  // Current idx1    
  Int_t curidx1 = 0; 
 
  for (i = fIndex->At(idx2); ((i < n2) && (curidx1 < idx1)); i++){
    if (fElements->At(i) < 0) {
      //curidx1 -= (Int_t) fElements->At(i);
      curidx1 -= TMath::Nint(fElements->At(i));
    }
    else {      
      curidx1++;
    }
  }

  if ((curidx1 == idx1) && (fElements->At(i) > 0)) {
    return fElements->At(i);
  }
  else {
    return 0;
  }

}

//_____________________________________________________________________________
void AliTRDdataArrayF::SetData(Int_t row, Int_t col, Int_t time, Float_t value)
{
  //
  // Sets the data value at a given position of the array
  // Includes boundary checking
  //

  if ((row >= 0) && (col >= 0) && (time >= 0)) {
    Int_t idx1 = GetIdx1(row,col);
    if ((idx1 >= 0) && (time < fNdim2)) {
      SetDataFast(idx1,time,value);
    }
    else {
      if (idx1 >= 0) {
        AliError(Form("time %d out of bounds (size: %d, this: 0x%08x)"
                     ,time,fNdim2,this));
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDdataArrayF::SetDataFast(Int_t idx1, Int_t idx2, Float_t value)
{
  //
  // Sets the data value at a given position of the array
  // No boundary checking
  //

  (*fElements)[fIndex->fArray[idx2]+idx1] = value;

}

//_____________________________________________________________________________
AliTRDdataArrayF &AliTRDdataArrayF::operator=(const AliTRDdataArrayF &a)
{
  //
  // Assignment operator
  //

  if (this != &a) ((AliTRDdataArrayF &) a).Copy(*this);
  return *this;

}

