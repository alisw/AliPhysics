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

/* $Id:  */

//____________________________________________________________________
//                                                                          
// FIT 
// Class for reading FIT RAW data in TOF data format
//
#include "AliFITRawReader.h"
#include "AliBitPacking.h"
#include "TBits.h"

#include <Riostream.h>
#include "TMath.h"
#include "TH1F.h"
#include "AliLog.h"

using std::cout;
using std::endl;
using std::ios_base;
ClassImp(AliFITRawReader)
  
  AliFITRawReader::AliFITRawReader (AliRawReader *rawReader)
    :  TTask("FITRawReader","read FIT raw  data"),
       fRawReader(rawReader),
       fData(NULL),
       fPosition(0),
       fBunchID(0),
       fPrintout(kFALSE)
    
{
  //
  // create an object to read FIT raw digits
  
  fRawReader->Reset();
  fRawReader->Select("FIT"); 
  for ( Int_t k=0; k<1000; k++)   fAllData[k] = -1;
 
}
//_____________________________________________________________________________
AliFITRawReader::~AliFITRawReader ()
{
  // 
}

//_____________________________________________________________________________

Bool_t  AliFITRawReader::Next()
{
  // read the next raw digit
  // returns kFALSE if there is no digit left
// allData array collect data from all channels in one :
// allData[0] - allData[79]   CFD channels side C
// allData[80] - allData[159] CFD channels side A
// allData[160] -   allData[239]  QT0 channels side C
// allData[240] -   allData[319]  QT0 channels side A
// allData[320] -   allData[399]  QT1 channels side C
// allData[400] -   allData[479]  QT1 channels side A

  UInt_t word;
  Int_t time=0,  itdc=0, ichannel=0, uu; 
  Int_t numberOfWordsInTRM=0, iTRM=0;
  Int_t tdcTime, koef;
  Int_t trm_chain_header =  0x00000000;
  Int_t  trm_chain_trailer[2] = { 0x10000000, 0x20000000};
  UInt_t  filler =  0x70000000;
  Bool_t correct=kTRUE;
  Int_t header;
  
  for ( Int_t k=0; k<1000; k++)   fAllData[k] = -1;
  
  
  do {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);
  
  fPosition = 0;
  cout.setf( ios_base::hex, ios_base::basefield );
  if(fPrintout)
    cout<<" CDH :: BC ID "<< (fRawReader->GetBCID())<<
      " Event size"<<fRawReader->GetDataSize()<<
      " orbit ID "<<fRawReader->GetOrbitID()<< 
      " event index "<<fRawReader->GetEventIndex()<<
      " event type " <<fRawReader->GetType()<<endl;
  //DRM header
  for (Int_t i=0; i<6; i++) {
    word = GetNextWord();
    if(fPrintout && i==0) cout<<" DRM header:: event words "<<AliBitPacking::UnpackWord(word,4, 20)<<endl;;
    if (fPrintout && i==4 ) cout<<" L0BC ID "<< AliBitPacking::UnpackWord(word,4, 15)<<endl;
    header = AliBitPacking::UnpackWord(word,28,31);
    if( header !=4 ) 
      {
	AliWarning(Form(" !!!! wrong  DRM header  %x!!!!", word));
	fRawReader->AddFatalErrorLog(kWrongDRMHeader,Form("w=%x",word));
	break;
      }
  }
  for (Int_t ntrm=0; ntrm<4; ntrm++)
    {
      //TRMheader  
      word = GetNextWord();
      header = AliBitPacking::UnpackWord(word,28,31);
      if ( header != 4 )
	{
	  AliWarning(Form(" !!!! wrong TRM header  %x!!!!", word));
	  fRawReader->AddMajorErrorLog(kWrongTRMHeader,Form("w=%x",word));
	  break;
	}
      numberOfWordsInTRM=AliBitPacking::UnpackWord(word,4,16);
      if(fPrintout) {
	cout<<" TRM header :: event words "<<numberOfWordsInTRM;
	cout<<" ACQ bits "<<AliBitPacking::UnpackWord(word,17,18);
	cout<<" L bit "<<AliBitPacking::UnpackWord(word,19,19)<<endl;
      }
      iTRM=AliBitPacking::UnpackWord(word,0,3);
      for( Int_t ichain=0; ichain<2; ichain++)
	{
	  //chain header
	  word = GetNextWord();
	  uu = word & trm_chain_header;
	  if(uu != trm_chain_header) 
	    {
	      AliWarning(Form(" !!!! wrong CHAIN  0  header %x!!!!", word));
	      fRawReader->AddMajorErrorLog(kWrongChain0Header,Form("w=%x",word));
	      break;
	    }
	  fBunchID=AliBitPacking::UnpackWord(word,4,15);
	  if(fPrintout)
	    cout<<" chain "<< ichain<<" header:: BunchID  "<<fBunchID<<endl;
	  word = GetNextWord();
	  tdcTime =  AliBitPacking::UnpackWord(word,31,31);   
	  while(tdcTime==1)
	    {
	      correct = kTRUE;
	      itdc=AliBitPacking::UnpackWord(word,24,27);
	      ichannel=AliBitPacking::UnpackWord(word,21,23);
	      time=AliBitPacking::UnpackWord(word,0,20);
	      
	      koef = GetChannel(iTRM,itdc,ichain,ichannel);
	      if (koef != 0 && fPrintout) 
		cout<<"RawReader>> "<<"koef "<<koef<<" trm "<<iTRM<<" tdc "<<itdc<<" chain "<<ichain<<		    " channel "<<ichannel<<" time "<<time<<endl;
	      if (koef ==-1 ){
		AliWarning(Form("Incorrect lookup table ! "));
		fRawReader->AddMajorErrorLog(kIncorrectLUT);
		correct=kFALSE;
	      }
	      if(correct)   fAllData[koef]=time; 
	      word = GetNextWord();
		tdcTime =  AliBitPacking::UnpackWord(word,31,31);   
	    }
	    
	  uu = word&trm_chain_trailer[ichain];
	  if(uu != trm_chain_trailer[ichain] )
	    {
	      AliWarning(Form(" !!!! wrong CHAIN %i trailer %x !!!!", ichain, word));
	      fRawReader->AddMajorErrorLog(kWrongChain0Trailer,Form("w=%x",word));
	      break;
	    }
	  if(fPrintout)
	    cout<<" trailer:: event counter "<< AliBitPacking::UnpackWord(word,16,27)<<endl;
	}
      
	word = GetNextWord(); //TRM trailer
	header = AliBitPacking::UnpackWord(word,28,31);
	if( header !=5 )
	  {
	    AliWarning(Form(" !!!! wrong TRM GLOBAL trailer  %x!!!!", word));
	    fRawReader->AddMajorErrorLog(kWrongTRMTrailer,Form("w=%x",word));
	    break;
	  }
	if(fPrintout)
	  cout<<"  TRM trailer :: event counter "<< AliBitPacking::UnpackWord(word,16,27)<<endl;
    } //TRM loop
  word = GetNextWord(); //
  if (word == filler )  word = GetNextWord(); 
  header = AliBitPacking::UnpackWord(word,28,31);
  if( header !=5 )
    {
      AliWarning(Form(" !!!! wrong DRM GLOBAL trailer  %x!!!!", word));
      // fRawReader->AddFatalErrorLog(kWrongDRMTrailer,Form("w=%x",word));
    }
  if(fPrintout)
    cout<<" DRM trailer ::event counter "<< AliBitPacking::UnpackWord(word,4,15)<<endl;
  cout.setf( ios_base::dec, ios_base::basefield );
  
  return kTRUE;
}
//_____________________________________________________________________________
Int_t AliFITRawReader::GetPosition()
{
  // Sets the position in the
  // input stream
  if (((fRawReader->GetDataSize() * 8) % 32) != 0)
    AliFatal(Form("Incorrect raw data size ! %d words are found !",fRawReader->GetDataSize()));
  return (fRawReader->GetDataSize() * 8) / 32;
}
//_____________________________________________________________________________
UInt_t AliFITRawReader::GetNextWord()
{
  // Read the next 32 bit word in backward direction
  // The input stream access is given by fData and fPosition
  
  
  //   fPosition--;
  Int_t iBit = fPosition * 32;
  Int_t iByte = iBit / 8;
  
  UInt_t word = 0;
  word  = fData[iByte+3]<<24;
  word |= fData[iByte+2]<<16;
  word |= fData[iByte+1]<<8;
  word |= fData[iByte];
  fPosition++;
  
  return word;   
}

//_____________________________________________________________________________
UInt_t AliFITRawReader::GetChannel(Int_t iTRM, Int_t iTDC, Int_t iChain, Int_t ichannel)
{
  return   iTRM*120 + iChain*60 + iTDC*4 +ichannel;
}
