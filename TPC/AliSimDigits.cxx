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
Revision 1.1.4.2  2000/04/10 11:37:42  kowal2

Digits handling in a new data structure

*/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Alice segment manager object                                             //
//                                
//  AliSimDigits object   (derived from AliDigits)                            //
//  provide additional track information to digit                            //
//  
//   Origin: Marian Ivanov  GSI Darmstadt                                     //
//
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "TClass.h"
#include <iostream.h>
#include "TError.h"
#include "AliSegmentID.h"
#include "AliH2F.h"
#include "AliArrayI.h"
#include "AliArrayS.h"
#include "AliDigits.h"
#include "AliSimDigits.h"
#include "AliTPC.h"



//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
ClassImp(AliSimDigits)

AliSimDigits::AliSimDigits()
{
  //  AliDigits::Invalite();
  fTracks = 0;
  fTrIndex = 0;  
  InvalidateTrack();
}
AliSimDigits::~AliSimDigits()
{
  //  AliDigits::Invalidate();
  if (fTracks != 0) fTracks->Delete();
  if (fTrIndex != 0) fTrIndex->Delete();  

}

void AliSimDigits::InvalidateTrack() 
{ 
  //
  //set default (invalid parameters)
  if ( fTracks != 0) delete fTracks;
  fTracks = new AliArrayI;
  if ( fTrIndex  != 0) delete fTrIndex;
  fTrIndex = new AliArrayI;

  for (Int_t i = 0; i<3; i++){
    fTracks->Set(0);
    fTrIndex->Set(0);
  }
}

void  AliSimDigits::AllocateTrack(Int_t length)
{
  //
  //construct empty buffer fElements and fTracks with size fNrows x fNcols x
  //length 
  InvalidateTrack();
  fNlevel = length;
  fTracks->Set(fNcols*fNrows*fNlevel);
  fTrIndex->Set(0); 
  fTrBufType =0;
}

Int_t AliSimDigits::GetTrackID(Int_t row, Int_t column, Int_t level) 
{
  //
  //Get track ID 
  if (fTrBufType == 0) return  GetTrackIDFast(row, column,level);
  if (fTrBufType == 1) return  GetTrackID1(row, column,level); 
  if (fTrBufType == 2) return  GetTrackID2(row, column,level); 
  return 0;
}

void AliSimDigits::ExpandTrackBuffer()
{  
  //
  //expand buffer to two dimensional array 
  if (fTrBufType<0)  {
    Error("ExpandBuffer", "buffer doesn't exist");
    return;
  }
  if (fTrBufType==0)      return;  //buffer expanded
  if (fTrBufType==1)  {ExpandTrackBuffer1(); return;}
  if (fTrBufType==2)  ExpandTrackBuffer2();
 
}

void AliSimDigits::CompresTrackBuffer(Int_t bufType)
{
  //
  //compres buffer according buffertype algorithm
  //
  if (fTrBufType<0)  {
    Error("CompressBuffer", "buffer doesn't exist");
    return;
  }
  if (fTrBufType == bufType) return;
  //
  if (fTrBufType>0) ExpandTrackBuffer();
  if (fTrBufType !=0)  {
    Error("CompressBuffer", "buffer doesn't exist");
    return;
  }
  //compress buffer of type 1
  
  if (bufType==1)      {CompresTrackBuffer1();return;}
  if (bufType==2)      CompresTrackBuffer2();
   
}

Int_t  AliSimDigits::GetTrackID1(Int_t row, Int_t column, Int_t level)
{
  //return  track ID of digits - for buffer compresion 2
  Int_t i,n1,n2;
  i = level*fNcols+column;
  if ( (i+1)>=fTrIndex->fN) n2 = fTracks->fN;
  else 
    n2 = fTrIndex->At(i+1);
  n1 = fTrIndex->At(i);
  Int_t rownew = 0;
  Int_t rowold=0;
  Int_t ID;
  for (i = n1;(i<n2);i++){
    ID = 0;
    Int_t num = fTracks->At(i);
    if (num<0) {
      rownew-=num;  
      rowold = rownew;
      i++;
      if (i<n2){
	num = fTracks->At(i);
	rownew+=num;
	i++;
	ID = fTracks->At(i);
      }
    }
    else {
      rowold = rownew;
      rownew+=num;
      i++;
      ID = fTracks->At(i);
    }
ID-=2;
    if ( (row>=rowold) && (row<=rownew) ) return ID;
    if (row < rownew ) return -2; //empty track
  }
  return -2;
}

void  AliSimDigits::ExpandTrackBuffer1()
{
  //
  //expand  track compressed according algorithm 1 (track ID comression independent to the digit compression)
  // !!in expanded tracks we don't use fTrIndex array
  //  
  fTrBufType = 0;
  Int_t i,j;
  Int_t all   = fNrows*fNcols;  //total number of digits
  Int_t elems = all*fNlevel;  //length of the buffer

  AliArrayI * buf = new AliArrayI;
  buf->Set(elems);
  fTrIndex->Set(0);
  //
  Int_t level = 0;
  Int_t col=0;
  Int_t row = 0;
  Int_t N=fTracks->fN;
  //
  for (i=0;i<N;i++){
    //oposite signa means how many unwrited (under threshold) values
    Int_t num = fTracks->At(i);
    if (num<0) row-=num;   //negative number mean number of zeroes (no tracks of gibven level no need to write to array) 
    else {
      i++;
      Int_t ID =  fTracks->At(i);
      for (j = 0; j<num; j++,row++) (*buf)[level*all+col*fNrows+row]=ID;       
    }
    if (row>=fNrows) {
      row=0;
      col++;
    }
    if (col>=fNcols) {
      col=0;
      level++;
    }    
  }//end of loop over digits
  delete fTracks;
  fTracks = buf;
}

void  AliSimDigits::CompresTrackBuffer1()
{
  //
  //comress track according algorithm 1 (track ID comression independent to the digit compression)
  //
  fTrBufType = 1;  

  AliArrayI *  buf = new AliArrayI;   //create  new buffer 
  buf->Set(fNrows*fNcols*fNlevel); //lets have the nearly the "worst case"
  AliArrayI *  index = new AliArrayI;
  index->Set(fNcols*fNlevel);

  Int_t icurrent=-1;  //current index
  Int_t izero;      //number of zero
  Int_t inum;      //number of digits  with the same current track ID  
  Int_t lastID =0;  //last track ID
  for (Int_t lev =0; lev<fNlevel; lev++){    //loop over levels 
    for (Int_t col = 0; col<fNcols; col++){    //loop over columns
      izero = 0;
      inum =  0;
      lastID = 0;
      (*index)[lev*fNcols+col]=icurrent+1;//set collumn pointer
      Int_t ID;  //current ID
      for (Int_t row = 0; row< fNrows;row++){ //loop over rows
	ID = GetTrackIDFast(row,col,lev);
	if (ID <= 0) {
	  if ( inum> 0 ) { //if we have some tracks in buffer
	    icurrent++;
	    if ((icurrent+1)>=buf->fN) buf->Expand(icurrent*2);
	    (*buf)[icurrent] = inum;
	    icurrent++;
	    (*buf)[icurrent] = lastID;	
	    inum = 0;      
	    lastID = 0;
	  }
	  izero++;
	}
	else
          if (ID != lastID) 
	    if ( izero > 0 ) { 
	      //if we have currently izero count of non tracks digits
	      icurrent++;	  
	      if (icurrent>=buf->fN) buf->Expand(icurrent*2);
	      (*buf)[icurrent]= -izero;  //write how many under zero
	      inum++;
	      izero = 0;	     
	      lastID = ID;
	    }
	    else{ 
	      //if we change track ID from another track ID	    
	      icurrent++;	  
	      if ((icurrent+1)>=buf->fN) buf->Expand(icurrent*2);
	      (*buf)[icurrent] = inum;
	      icurrent++;
	      (*buf)[icurrent] = lastID;	
	      lastID = ID;
	      inum = 1;      
	      izero = 0;
	    }	
	  else {	  
	    inum++;
	  }
	
      }//end of loop over rows
      if ( izero > 0 ) { 
	//if we have currently izero count of non tracks digits
	icurrent++;	  
	if (icurrent>=buf->fN) buf->Expand(icurrent*2);
	(*buf)[icurrent]= -izero;  //write how many under zero	
      }
      if ( inum> 0 ) { //if we have some tracks in buffer
	icurrent++;
	if ((icurrent+1)>=buf->fN) buf->Expand(icurrent*2);
	(*buf)[icurrent] = inum;
	icurrent++;
	(*buf)[icurrent] = ID;	
      }      
    }//end of loop over columns
  }//end of loop over differnet track level  
  buf->Expand(icurrent+1);
  delete fTracks;
  fTracks = buf;
  delete fTrIndex;
  fTrIndex = index;
}



void  AliSimDigits::ExpandTrackBuffer2()
{
  //
  //comress track according algorithm 2 (track ID comression according  digit compression)
  fTrBufType = 0;
}

void  AliSimDigits::CompresTrackBuffer2()
{
  //
  //comress track according algorithm 2 (track ID comression according  digit compression)
  fTrBufType = 2;
}


Int_t  AliSimDigits::GetTrackID2(Int_t row, Int_t column, Int_t level)
{
  //returnb track ID of digits - for buffer compresion 2
  return -2;
}



AliH2F *  AliSimDigits::DrawTracks( const char *option,Int_t level, 
			      Float_t x1, Float_t x2, Float_t y1, Float_t y2)
{
  //
  //draw digits in given array
  //  
  //make digits histo 
  char ch[30];
  sprintf(ch,"Track Segment_%d level %d ",GetID(),level );
  if ( (fNrows<1)|| (fNcols<1)) {
    return 0;
  }
  AliH2F * his  = new AliH2F("Track histo",ch,fNrows,0,fNrows,fNcols,0,fNcols);
  ExpandTrackBuffer();
  //set histogram  values
  for (Int_t i = 0; i<fNrows;i++)    
    for (Int_t j = 0; j<fNcols;j++)
        his->Fill(i,j,GetTrackIDFast(i,j,level));
  if (x1>=0) {
      AliH2F *h2fsub = his->GetSubrange2d(x1,x2,y1,y2);
      delete his;
      his=h2fsub;
  }
  if (his==0) return 0;
  if (option!=0) his->Draw(option);
  else his->Draw();
  return his;  
}

TClonesArray *  AliSimDigits::GenerTPCClonesArray(TClonesArray * arr)
{
  //
  //generate TClonnesArray of digits
  //
  TClonesArray * digits;
  if (arr==0)  digits=new TClonesArray("AliTPCdigit",300);
  else digits = arr; 
  Int_t index = digits->GetEntriesFast();
  for (Int_t row =0; row<fNrows; row++)
    for (Int_t col =0; col<fNcols; col++){
      Int_t amp = GetDigit(row,col);
      if (amp>GetThreshold()){
	AliTPCdigit dig;
	dig.fPad = col;
	dig.fTime = row;
	dig.fSignal= amp;
	dig.fPadRow =fSegmentID;
	dig.fSector =fSegmentID;
	dig.fTracks[0]= GetTrackID(row,col,0);
	dig.fTracks[1]= GetTrackID(row,col,1);
	dig.fTracks[2]= GetTrackID(row,col,2);
	TClonesArray &ldigits = *digits;
	new(ldigits[index++]) AliTPCdigit(dig);
      }
    }    
  return digits;
}

