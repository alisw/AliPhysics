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

/*
$Log$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  General container for data of a TRD detector segment.                    //
//  Adapted from AliDigits (origin: M.Ivanov).                               //
//  The main difference is that we used 4 byte integer, so that this class   //
//  can also be used as a dictionary between digits and MC particles.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TClass.h"
#include "TError.h"
#include "AliTRDsegmentID.h"
#include "AliTRDarrayI.h"
#include "AliTRDdataArray.h"

ClassImp(AliTRDdataArray)

//_____________________________________________________________________________
AliTRDdataArray::AliTRDdataArray()
{
  //
  // Default constructor
  //

  fIndex     = 0;
  fElements  = 0;
  fThreshold = 0;
  Reset();

}

//_____________________________________________________________________________
AliTRDdataArray::~AliTRDdataArray()
{
  //
  // Destructor
  //

  if (fIndex)    fIndex->Delete();;
  if (fElements) fElements->Delete();
  
}


//_____________________________________________________________________________
void AliTRDdataArray::Reset() 
{ 
  //
  // Reset the array (old content gets deleted)
  //

  if (fIndex)    delete fIndex;
  fIndex    = new AliTRDarrayI;
  
  if (fElements) delete fElements;
  fElements = new AliTRDarrayI;

  fNdim1 = fNdim2 = fNelems = -1; 
  fElements->Set(0); 
  fIndex->Set(0); 
  fBufType = -1;

  fNrow  = 0;
  fNcol  = 0;
  fNtime = 0;

}

//_____________________________________________________________________________
void AliTRDdataArray::Allocate(Int_t nrow, Int_t ncol, Int_t ntime)
{
  //
  // Allocate an empty buffer of the size <nrow> x <ncol> x <ntime>
  //

  Reset();

  if (nrow  <= 0) {
    Error("AliTRDdataArray::Allocate","The number of rows has to be positive");
    return;
  }
  if (ncol  <= 0) {
    Error("AliTRDdataArray::Allocate","The number of columns has to be positive");
    return;
  }
  if (ntime <= 0) {
    Error("AliTRDdataArray::Allocate","The number of timebins has to be positive");
    return;
  }

  // The two-dimensional array row/column gets mapped into the first 
  // dimension of the array. The second array dimension, which is not compressible,
  // corresponds to the time direction
  fNdim1  = nrow * ncol;
  fNdim2  = ntime;
  fNelems = fNdim1 * fNdim2;

  fNrow   = nrow;
  fNcol   = ncol;
  fNtime  = ntime;

  fElements->Set(fNelems);
  fIndex->Set(fNdim2);

  for (Int_t i = 0, k = 0; i < fNdim2; i++, k += fNdim1) { 
    (*fIndex)[i] = k;
  }

  fBufType = 0;

}

//_____________________________________________________________________________
Int_t AliTRDdataArray::GetSize()
{
  //
  // Returns the size of the complete object
  //

  Int_t size = sizeof(this);

  if (fIndex)    size += sizeof(fIndex)    
                         + fIndex->GetSize()    * sizeof(Int_t);
  if (fElements) size += sizeof(fElements) 
                         + fElements->GetSize() * sizeof(Int_t);

  return size;

}

//_____________________________________________________________________________
Int_t AliTRDdataArray::GetDataSize()
{
  //
  // Returns the size of only the data part
  //

  if (fElements == 0) 
    return 0;
  else 
    return sizeof(fElements) + fElements->GetSize() * sizeof(Int_t);

}

//_____________________________________________________________________________
Int_t AliTRDdataArray::GetOverThreshold(Float_t threshold)
{
  //
  // Returns the number of entries over threshold
  //
 
  if ((fElements == 0) || (fElements->GetSize() <= 0))
    return 0;
 
  Int_t  over = 0;

  for (Bool_t cont = First(); cont == kTRUE; cont = Next()) {
    if ((fCurrentIdx1 < 0) || (fCurrentIdx1 > fNdim1)) continue;
    if ((fCurrentIdx2 < 0) || (fCurrentIdx2 > fNdim2)) continue;
    if (fElements->At(fCurrentIndex) > threshold) over++;
  }

  return over;

}

//_____________________________________________________________________________
Int_t AliTRDdataArray::GetData(Int_t row, Int_t col, Int_t time)
{
  //
  // Returns the data value at a given position of the array
  //

  if (fBufType == 0) return GetDataFast(GetIdx1(row,col),time);
  if (fBufType == 1) return GetData1(GetIdx1(row,col),time);

  return 0;

}

//_____________________________________________________________________________
void AliTRDdataArray::SetData(Int_t row, Int_t col, Int_t time, Int_t value)
{
  //
  // Sets the data value at a given position of the array
  //

  SetDataFast(GetIdx1(row,col),time,value);

}

//_____________________________________________________________________________
void AliTRDdataArray::Expand()
{  
  //
  // Expands the compressed buffer
  //

  if (fBufType  < 0) {
    Error("AliTRDdataArray::Expand","Buffer does not exist");
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
void AliTRDdataArray::Compress(Int_t bufferType)
{
  //
  // Compresses the buffer
  //

  if (fBufType  < 0) {
    Error("AliTRDdataArray::Compress","Buffer does not exist");
    return;
  }
  if (fBufType == bufferType) {
    return;
  }  
  if (fBufType > 0) {
    Expand();
  }
  if (fBufType !=0)  {
    Error("AliTRDdataArray::Compress","Buffer does not exist");
    return;
  }

  // Compress a buffer of type 1
  if (bufferType == 1) {
    Compress1();
  }

}

//_____________________________________________________________________________
void AliTRDdataArray::Compress(Int_t bufferType, Int_t threshold)
{
  //
  // Compresses the buffer
  //

  fThreshold = threshold;
  Compress(bufferType);

}

//_____________________________________________________________________________
Bool_t AliTRDdataArray::First()
{
  //
  // Returns the position of the first valid data value
  //

  if (fBufType == 0) return First0();
  if (fBufType == 1) return First1();
  return kFALSE;

}

//_____________________________________________________________________________
Bool_t  AliTRDdataArray::Next()
{
  //
  // Returns the position of the next valid data value
  //

  if (fBufType == 0) return Next0();
  if (fBufType == 1) return Next1();
  return kFALSE;

}
 
//_____________________________________________________________________________
void AliTRDdataArray::Expand1()
{
  //
  // Expands a buffer of type 1
  //

  Int_t i, k;

  fNelems = fNdim1 * fNdim2;

  Int_t *buf = new Int_t[fNelems];

  fIndex->Set(fNdim2);

  for (i = 0, k = 0; i < fNdim2; i++, k += fNdim1) (*fIndex)[i] = k;

  Int_t idx1 = 0;
  Int_t idx2 = 0;
  Int_t N    = fElements->fN;

  for (i = 0; i < N; i++){

    // Negative sign counts the unwritten values (under threshold)
    if ((*fElements)[i] < 0) {
      idx1 -= fElements->At(i);
    } 
    else {
      buf[(*fIndex)[idx2] + idx1] = fElements->At(i);
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
void AliTRDdataArray::Compress1()
{
  //
  // Compress a buffer of type 1
  //

  AliTRDarrayI  buf;  
  buf.Set(fNelems);
  AliTRDarrayI  index;
  index.Set(fNdim2);

  Int_t icurrent = -1;
  Int_t izero;
  for (Int_t idx2 = 0; idx2 < fNdim2; idx2++){      

    // Set the idx2 pointer
    index[idx2] = icurrent + 1;

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
	  if (icurrent >= buf.fN) buf.Expand(icurrent*2);
          // Store the number of entries below zero
	  buf[icurrent] = -izero;  
	  izero = 0;
	} 
	icurrent++;
	if (icurrent >= buf.fN) buf.Expand(icurrent*2);
	buf[icurrent] = GetDataFast(idx1,idx2);	    
      } // If signal larger than threshold	  	
    } // End of loop over idx1

    if (izero > 0) {
      icurrent++;	  
      if (icurrent >= buf.fN) buf.Expand(icurrent*2);
      // Store the number of entries below zero
      buf[icurrent] = -izero; 
    }

  }

  buf.Expand(icurrent+1);
  (*fElements) = buf;
  fNelems   = fElements->fN;
  fBufType  = 1;
  (*fIndex) = index;

}

//_____________________________________________________________________________
void AliTRDdataArray::Expand2()
{
  //
  // Expands a buffer of type 2 
  //

  Int_t i, k;
  Int_t *buf = new Int_t[fNelems];

  fNelems = fNdim1 * fNdim2;
  fIndex->Set(fNdim2);

  for (i = 0, k = 0; i < fNdim2; i++, k += fNdim1) (*fIndex)[i] = k;

  Int_t idx1 = 0;
  Int_t idx2 = 0;
  Int_t N    = fElements->fN;
  for (i = 0; i < N; i++){
    // Negative sign counts the unwritten values (under threshold)
    if ((*fElements)[i] < 0) {
      idx1 -= fElements->At(i); 
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
void AliTRDdataArray::Compress2()
{
  /*

  //
  // Compress a buffer of type 2
  //

  AliArrayS  buf;  //lets have the nearly the "worst case"
  buf.Set(fNelems);
  AliTRDarrayI  index;
  index.Set(fNdim2);
  Int_t icurrent=-1;
  Int_t izero;

  for (Int_t col = 0; col<fNdim2; col++){      
    index[col]=icurrent+1;//set collumn pointer
    izero = 0;  //reset zer counter at the begining of the column
    Int_t lastrow=0;
    Int_t lastrowval=GetDigitFast(row,0);

    for (Int_t row = 1; row< fNdim1;row++){
      //if under threshold
      Int_t val = GetDigitFast(row,col);
      Int_t dif = val -lastrowval;
      
      if (TMath::Abs(dif)<fThreshold)  izero++;
      else{
	if (izero>0) {
	  //if we have currently izero count under threshold
	  icurrent++;	  
	  if (icurrent>=buf.fN) buf.Expand(icurrent*2);
	  buf[icurrent]= -izero;  //write how many under zero
	  izero = 0;
	} //end of reseting izero
	icurrent++;
	if (icurrent>=buf.fN) buf.Expand(icurrent*2);
	buf[icurrent] = GetDigitFast(row,col);	    
      }//if signal bigger then threshold	  	
    } //end of loop over rows
    
    if (izero>0) {
      icurrent++;	  
      if (icurrent>=buf.fN) buf.Expand(icurrent*2);
      buf[icurrent]= -izero;  //write how many under zero
    }
  }//end of lopping over digits
  buf.Expand(icurrent+1);
  (*fElements)=buf;
  fNelems = fElements->fN;
  fBufType = 1;
  (*fIndex) =index;
  //end of compresing bufer of type 1 

  */

}

//_____________________________________________________________________________
Bool_t AliTRDdataArray::First0()
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
Bool_t AliTRDdataArray::Next0()
{
  //
  // Returns the next entry for a buffer of type 0
  //

  if (fCurrentIndex < 0) return kFALSE; 

  Int_t i;
  for (i = fCurrentIndex + 1; 
       ((i < fNelems) && (fElements->At(i) <= fThreshold)); 
       i++);
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
Bool_t AliTRDdataArray::First1()
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
      fCurrentIdx1-=fElements->At(i);
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
Bool_t AliTRDdataArray::Next1()
{
  //
  // Returns the next entry for a buffer of type 1
  //

  if (fCurrentIndex < 0) return kFALSE;

  Int_t i;
  for (i = fCurrentIndex + 1; i < fNelems; i++){
    if (fElements->At(i) < 0) {
      fCurrentIdx1 -= fElements->At(i);
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
Int_t AliTRDdataArray::GetData1(Int_t idx1, Int_t idx2)
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
      curidx1 -= fElements->At(i);
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

