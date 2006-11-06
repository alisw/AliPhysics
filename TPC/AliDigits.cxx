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

/*MI change -- for Rule checker
          -- added copy constructor and assignmet operator 
	  -- new GetSize return size of object in Bytes
          -- added GetDigitSize and GetOverTh function
	  -- added GetNRows, GetNCols function
          -- for Marek -I had it in my code  
*/ 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Alice  digits array  object  AliDigits                                  //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TClass.h"
#include <Riostream.h>
#include "TError.h"
#include "AliSegmentID.h"
#include "AliH2F.h"
#include <TArrayI.h>
#include <TArrayS.h>
#include "AliDigits.h"



//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
ClassImp(AliDigits)


 AliDigits::AliDigits()
           :AliSegmentID(),
            fNrows(0),
            fNcols(0),
	    fElements(0),
            fIndex(0),
            fBufType(0),
            fThreshold(0),
            fNelems(0),
            fCurrentRow(0),
            fCurrentCol(0),
            fCurrentIndex(0) 
{
  // 
  //default constructor
  //
  Invalidate();
}

AliDigits::AliDigits(const AliDigits& digits)
          :AliSegmentID(digits),
            fNrows(0),
            fNcols(0),
	    fElements(0),
            fIndex(0),
            fBufType(0),
            fThreshold(0),
            fNelems(0),
            fCurrentRow(0),
            fCurrentCol(0),
            fCurrentIndex(0)
{
  //
  //copy constructor
  //
  fNrows = digits.fNrows;
  fNcols = digits.fNcols;
  fElements = new TArrayS(*(digits.fElements));
  fIndex = new TArrayI(*(digits.fIndex));
  fBufType = digits.fBufType;
  fThreshold = digits.fThreshold;
  fNelems    = digits.fNelems;
}

AliDigits & AliDigits::operator =(const AliDigits & digits)
{
 //assignment operator
  fNrows = digits.fNrows;
  fNcols = digits.fNcols;
  if (fElements) delete fElements;
  fElements = new TArrayS(*(digits.fElements));
  if (fIndex) delete fIndex;
  fIndex = new TArrayI(*(digits.fIndex));
  fBufType = digits.fBufType;
  fThreshold = digits.fThreshold;
  fNelems    = digits.fNelems; 
  return (*this);
}

AliDigits::~AliDigits()
{
  //
  //default destructor
  if (fIndex !=0 ) {
    delete fIndex;
  }
  if (fElements != 0) {
    delete fElements;
  }
  
  
}


Bool_t AliDigits::OutOfBoundsError(const char *where, Int_t row, Int_t column) 
{
   // Generate an out-of-bounds error. Always returns false.
   ::Error(where, "row %d  col %d out of bounds (size: %d x %d, this: 0x%08x)", 
	   row, column, fNrows, fNcols, this);
   return kFALSE;
}


void AliDigits::Invalidate() 
{ 
  //
  //set default (invalid parameters)
  if (fIndex != 0)  delete  fIndex;
  fIndex = new TArrayI;
  
  if (fElements!= 0)     delete  fElements;
  
  fElements = new TArrayS;
  
  fNrows = fNcols =fNelems= -1; 
  fElements->Set(0); 
  fIndex->Set(0); 
  fBufType = -1;
}

void AliDigits::Allocate(Int_t rows, Int_t columns)
{
  //
  //construct empty buffer fDigits with size rows x columns
  Invalidate();
  if (rows <= 0) {
      Error("Allocate", "no of rows has to be positive");
      return;
   }
   if (columns <= 0) {
      Error("Allocate", "no of columns has to be positive");
      return;
   }
  fNrows = rows;
  fNcols=columns;
  fNelems = fNrows * fNcols;
  fElements->Set(fNelems);
  fIndex->Set(fNcols);
  for (Int_t i =0,k=0; i<fNcols;i++,k+=fNrows) 
  (*fIndex)[i]=k;
  fBufType =0;
}


Int_t AliDigits::GetSize()
{
  //
  //return size of object
  //
  Int_t size = sizeof(this);
  if (fIndex!=0) size+= sizeof(fIndex)+fIndex->GetSize()*sizeof(Int_t);
  if (fElements!=0) size+= sizeof(fElements)+fElements->GetSize()*sizeof(Short_t);
  return size;
}

Int_t AliDigits::GetDigitSize() //return total size of pure digit
{
  //
  //return size of PURE DIGITS
  //
  if (fElements==0) return 0;
  else return sizeof(fElements)+fElements->GetSize()*sizeof(Short_t);
}

Int_t AliDigits::GetOverTh(Float_t threshold,Float_t x1, Float_t x2, Float_t y1, Float_t y2)
{
  //
  //return number of digits over threshold
  // 
 if ( (fElements==0) || (fElements->GetSize()<=0)) return 0;
 
 if (x1<=x2) {
    x1=0;
    x2=fNrows;
  }
  if (y1<=y2) {
     y1=0;
     y2=fNcols;
  }
  Int_t over=0;

  Bool_t cont=First();
  for ( cont=First(); cont==kTRUE;cont=Next()) {
    if ( (CurrentRow()<x1) || (CurrentRow()>x2)) continue;
    if ( (CurrentColumn()<y1) || (CurrentColumn()>y2)) continue;
    if (CurrentDigit()>threshold) over++;
  }
  return over;
}


Short_t AliDigits::GetDigit(Int_t row, Int_t column)
{
  //
  // return digit for given row and collumn
  if (fBufType ==0) return GetDigitFast(row,column);
  if (fBufType ==1) return GetDigit1(row,column);

  return 0;
}


void AliDigits::ExpandBuffer()
{  
  //
  //expand buffer to two dimensional array
  if (fBufType<0)  {
    Error("ExpandBuffer", "buffer doesn't exist");
    return;
  }
  if (fBufType==0)      return;
  
  //expanding of buffer type 1
  if (fBufType==1) ExpandBuffer1();
  
  fBufType = 0;
}

void AliDigits::CompresBuffer(Int_t bufferType,Int_t threshold)
{
  //
  //compres buffer according buffertype algorithm
  if (fBufType<0)  {
    Error("CompressBuffer", "buffer doesn't exist");
    return;
  }
  if (fBufType == bufferType) return;
  //
  if (fBufType>0) ExpandBuffer();
  if (fBufType !=0)  {
    Error("CompressBuffer", "buffer doesn't exist");
    return;
  }
  fThreshold = threshold;
  //compress buffer of type 1
  if ( bufferType == 1) CompresBuffer1();//end of compresing bufer of type 1 
}

Bool_t AliDigits::First()
{
  //adjust  first valid current digit
  if (fBufType ==0) return First0();
  if (fBufType ==1) return First1();
  return kFALSE;
}

Bool_t  AliDigits::Next()
{
  //addjust next valid current digit
  if (fBufType ==0) return Next0();
  if (fBufType ==1) return Next1();
  return kFALSE;
}
 
void AliDigits::AcceptHisto(AliH2F * his)
{
  //
  //make digits buffer with value according histograms values
  //for testing purpose  
  Int_t idim =his->GetNbinsX();
  Int_t jdim =his->GetNbinsY();
  if ( (idim<1)|| (jdim<1)) {
    return;
  }
  //allocate proper buffer size
  Allocate(idim,jdim);
  //set digits values
  for (Int_t i = 0; i<idim;i++)    
    for (Int_t j = 0; j<jdim;j++)
      {
        Int_t index = his->GetBin(i+1,j+1);     
        SetDigitFast((Short_t)his->GetBinContent(index),i,j);
      }   
}

AliH2F *  AliDigits::GenerHisto()
{
  //
  //make digits histo 
  char ch[30];
  sprintf(ch,"Segment_%d ",GetID());
  if ( (fNrows<1)|| (fNcols<1)) {
    return 0;
  }
  AliH2F * his  = new AliH2F("Digit histo",ch,fNrows,0,fNrows,fNcols,0,fNcols);
  ExpandBuffer();
  //set histogram  values
  for (Int_t i = 0; i<fNrows;i++)    
    for (Int_t j = 0; j<fNcols;j++)
        his->Fill(i,j,GetDigitFast(i,j));
  return his;
}

AliH2F *AliDigits::DrawDigits(const char *option,Float_t x1, Float_t x2, Float_t y1, Float_t y2)
{
  //
  //draw digits in given array
  //
  AliH2F *h2f = GenerHisto();
  if (x1>=0) {
      AliH2F *h2fsub = h2f->GetSubrange2d(x1,x2,y1,y2);
      delete h2f;
      h2f=h2fsub;
  }
  if (h2f==0) return 0;
  if (option!=0) h2f->Draw(option);
  else h2f->Draw();
  return h2f;  
}

void AliDigits::ExpandBuffer1()
{
  //
  //expand buffer of type to twodimensional array
  Int_t i,k;
  fNelems = fNrows*fNcols;
  Short_t * buf = new Short_t[fNelems];
  memset(buf,0,fNelems*sizeof(Short_t)); //MI change - 4.12.2000
  fIndex->Set(fNcols);
  for (i =0,k=0 ;i<fNcols;i++,k+=fNrows) (*fIndex)[i]=k;
  Int_t col=0;
  Int_t row = 0;
  Int_t n=fElements->fN;
  for (i=0;i<n;i++){
    //oposite signa means how many unwrited (under threshold) values
    if ((*fElements)[i]<0) row-=fElements->At(i); 
    else {
      buf[(*fIndex)[col]+row]=fElements->At(i);
      row++;
    }
    if (row==fNrows) {
      row=0;
      col++;
    }else 
      if (row>fNrows){
	Invalidate();
	return;
      }      
  }
  fElements->Adopt(fNelems,buf);    
}

void AliDigits::CompresBuffer1()
{
  //
  //compres buffer according  algorithm 1
  //
  TArrayS  buf;  //lets have the nearly the "worst case"
  buf.Set(fNelems);
  TArrayI  index;
  index.Set(fNcols);
  Int_t icurrent=-1;
  Int_t izero;
  Short_t * cbuff = fElements->GetArray();  //MI change

  for (Int_t col = 0; col<fNcols; col++){      
    index[col]=icurrent+1;//set collumn pointer
    izero = 0;  //reset zer counter at the begining of the column
    for (Int_t row = 0; row< fNrows;row++){
      //if under threshold
      //if (GetDigitFast(row,col)<=fThreshold)  izero++;
      if (*cbuff<=fThreshold)  izero++;

      else{
	if (izero>0) {
	  //if we have currently izero count under threshold
	  icurrent++;	  
	  if (icurrent>=buf.fN) buf.Set(icurrent*2);
	  buf[icurrent]= -izero;  //write how many under zero
	  izero = 0;
	} //end of reseting izero
	icurrent++;
	if (icurrent>=buf.fN) buf.Set(icurrent*2);
	//buf[icurrent] = GetDigitFast(row,col);	    
	buf[icurrent] = *cbuff;	    
      }//if signal bigger then threshold	
       cbuff++;
    } //end of loop over rows
    if (izero>0) {
      icurrent++;	  
      if (icurrent>=buf.fN) buf.Set(icurrent*2);
      buf[icurrent]= -izero;  //write how many under zero
    }
  }//end of lopping over digits
  buf.Set(icurrent+1);
  (*fElements)=buf;
  fNelems = fElements->fN;
  fBufType = 1;
  (*fIndex) =index;
  //end of compresing bufer of type 1 
}



Bool_t AliDigits::First0()
{
  //
  //first for the buffer type 0
  fCurrentRow = -1;
  fCurrentCol = -1;
  fCurrentIndex = -1;
  Int_t i;
  for (i=0; (( i<fNelems) && (fElements->At(i)<=fThreshold));i++);  //MI1211
  if (i == fNelems) return kFALSE;
  fCurrentCol =i/fNrows;
  fCurrentRow =i%fNrows;
  fCurrentIndex = i;
  return kTRUE;	
}

Bool_t AliDigits::Next0()
{
  //
  //next for the buffer type 0
  //
  if (fCurrentIndex<0) return kFALSE;  // if we didn't adjust first 
  Int_t i;
  for (i=fCurrentIndex+1; ( (i<fNelems) && (fElements->At(i)<=fThreshold) ) ;i++);
  if (i >= fNelems)  {
    fCurrentIndex = -1;
    return kFALSE;
  }
  fCurrentCol =i/fNrows;
  fCurrentRow =i%fNrows;
  fCurrentIndex = i;
  return kTRUE;	
}

Bool_t AliDigits::First1()
{
  //
  //first for the buffer type 1
  fCurrentRow = -1;
  fCurrentCol = 0;
  fCurrentIndex = -1;
  Int_t i;
  for (i=0; i<fNelems; i++){
    if (fElements->At(i) < 0) fCurrentRow-=fElements->At(i);
    else      
      fCurrentRow++;
    if (fCurrentRow>=fNrows) {
       fCurrentCol++;
       fCurrentRow-=fNrows;
    }
    if (fElements->At(i)>fThreshold) break;
  }
  fCurrentIndex = i;
  if (fCurrentIndex>=0&&i<fNelems) return kTRUE;
  fCurrentRow =-1;
  fCurrentCol =-1;
  return kFALSE;	
}

Bool_t AliDigits::Next1()
{
  //
  //next for the buffer type 1
  if (fCurrentIndex<0) return kFALSE;  // if we didn't adjust first 
  Int_t i;
  for (i=fCurrentIndex+1; i<fNelems;i++){
    if (fElements->At(i) < 0) fCurrentRow-=fElements->At(i);
    else      
      fCurrentRow++;
    if (fCurrentRow>=fNrows) {
      fCurrentCol++;
      fCurrentRow-=fNrows;
    }
    if (fElements->At(i)>fThreshold) break;
  }
  fCurrentIndex = i;
  if ( (i>=0) && (i<fNelems) ) return kTRUE;
  fCurrentRow =-1;
  fCurrentCol =-1;
  return kFALSE;
}

Short_t AliDigits::GetDigit1(Int_t row, Int_t column)
{
  //
  //return digit for given row and column  the buffer type 1
  //no control performed
  
  Int_t i,n2;
  if ( (column+1)>=fNcols) n2 = fNelems;
  else
    n2 = fIndex->At(column+1);
  Int_t irow = 0; //current row    
 
  for (i=fIndex->At(column); ( (i<n2) && (irow<row) );i++){
    if (fElements->At(i) < 0) irow-=fElements->At(i);
    else      
      irow++;
  }
  if ( irow == row ) return fElements->At(i);
  return -1;
}

