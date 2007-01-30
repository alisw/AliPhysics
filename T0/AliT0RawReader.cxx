#include "AliT0RawReader.h"
#include "AliT0Parameters.h"
#include "AliBitPacking.h"
#include "TBits.h"

#include <Riostream.h>
#include "TMath.h"
#include "TH1F.h"
#include "TArrayI.h"
#include "AliLog.h"
 
ClassImp(AliT0RawReader)
  
  AliT0RawReader::AliT0RawReader (AliRawReader *rawReader)
    :  TTask("T0RawReader","read raw T0 data"),
       fRawReader(rawReader),
       fData(NULL),
       fPosition(0)
{
  //
// create an object to read T0raw digits
  AliDebug(1,"Start ");
 
  fRawReader->Reset();
  fRawReader->Select("T0");
 
 
}
 AliT0RawReader::~AliT0RawReader ()
{
  // 
}

Bool_t  AliT0RawReader::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left
//"LookUpTable":
// Amplitude LED TRM=0; chain=0; TDC 0 -5        channel 0,2,4,6
// Time CFD      TRM=0; chain=0; TDC 6 - 11      channel 0,2,4,6
// mean time     TRM=0; chain=0; TDC 12          channel 0
// T0A           TRM=0; chain=0; TDC 12          channel 2
// T0C           TRM=0; chain=0; TDC 12          channel 4
// vertex        TRM=0; chain=0; TDC 12          channel 6
// mult QTC0        TRM=0; chain=0; TDC 13          channel 0
// mult QTC1        TRM=0; chain=0; TDC 13          channel 2

// Charge QTC0   TRM=1; chain=0; TDC 0 -5        channel 0,2,4,6
// Charge QTC1   TRM=1; chain=0; TDC 6 - 11      channel 0,2,4,6
// T0A trigger          TRM=1; chain=0; TDC 12          channel 0
// T0C trigger          TRM=1; chain=0; TDC 12          channel 2
// vertex trigger       TRM=1; chain=0; TDC 12          channel 4
// trigger central      TRM=1; chain=0; TDC 13          channel 0
// tigger semicenral    TRM=1; chain=0; TDC 13          channel 2
//
// allData array collect data from all channels in one :
// allData[0] - allData[23] 24 CFD channels
// allData[24] -   allData[47] 24 LED channels
//  allData[48]  mean (T0) signal  
// allData[49]   time difference (vertex)


 //  if (fDigits == 0x0) fDigits = new AliT0digit(); 
  // fTree->Branch("T0","AliT0digit",&fDigits,405,1);
 
  UInt_t word;
  Int_t time=0,  itdc=0, ichannel=0; 
  Int_t numberOfWordsInTRM=0, iTRM=0;
  Int_t tdcTime, koef,hit, meanTime, timeDiff ;



  AliT0Parameters* param = AliT0Parameters::Instance();   //-->Zhenya
  param->Init();
 
 for ( Int_t k=0; k<110; k++) {
    for ( Int_t jj=0; jj<5; jj++) {
      fAllData[k][jj]=0;
    }
  }
    do {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);
  
  //  fPosition = GetPosition();
  fPosition = 0;

  //DRM header
  for (Int_t i=0; i<4; i++) {
    word = GetNextWord();
  }
  //TRMheader  
   word = GetNextWord();
   numberOfWordsInTRM=AliBitPacking::UnpackWord(word,4,16);
   iTRM=AliBitPacking::UnpackWord(word,0,3);

   //chain header
   Int_t ichain=0;
   word = GetNextWord();
  
   for (Int_t i=0; i<numberOfWordsInTRM; i++) {
     word = GetNextWord();
     tdcTime =  AliBitPacking::UnpackWord(word,31,31);   

     if ( tdcTime == 1)
       {
	 itdc=AliBitPacking::UnpackWord(word,24,27);
	 ichannel=AliBitPacking::UnpackWord(word,21,23);
	 time=AliBitPacking::UnpackWord(word,0,20);
	 //  koef = itdc*4 + ichannel/2;
	 koef = param->GetChannel(iTRM,itdc,ichain,ichannel);
	 //	 cout<<" RawReader::Next ::"<<iTRM<<"  "<<itdc<<" "<<ichain<<" "<<ichannel<<" "<<  koef<<" "<<time<<endl;
	 if(fAllData[koef][0] == 0)  fAllData[koef][0]=time;  // yield only 1st particle
	  
       }
   }
   word = GetNextWord(); //chain trailer
   word = GetNextWord(); //TRM trailer
     
  //TRMheader  
   word = GetNextWord();
   numberOfWordsInTRM=AliBitPacking::UnpackWord(word,4,16);
   iTRM=AliBitPacking::UnpackWord(word,0,3);
   
   //chain header
   word = GetNextWord();
   
   for (Int_t iword=0; iword<numberOfWordsInTRM; iword++) {
     word = GetNextWord();
     tdcTime =  AliBitPacking::UnpackWord(word,31,31);   

     if ( tdcTime == 1)
       {
	 itdc=AliBitPacking::UnpackWord(word,24,27);
	 ichannel=AliBitPacking::UnpackWord(word,21,23);
	 time=AliBitPacking::UnpackWord(word,0,20);
	 //	 koef = itdc*4 + ichannel/2;
	 koef = param->GetChannel(iTRM,itdc,ichain,ichannel);
 
	   if(fAllData[koef][0] == 0)	 fAllData[koef][0]=time;
		 //	 if(allData[koef+55] == 0) allData[koef+55]=time; // yield only 1st particle
       }
   }
   meanTime = fAllData[49][0];  // T0 !!!!!!
   timeDiff = fAllData[50][0];

   word = GetNextWord();
   word = GetNextWord();
   return kTRUE;
}
//_____________________________________________________________________________
Int_t AliT0RawReader::GetPosition()
{
  // Sets the position in the
  // input stream
  if (((fRawReader->GetDataSize() * 8) % 32) != 0)
    AliFatal(Form("Incorrect raw data size ! %d words are found !",fRawReader->GetDataSize()));
  return (fRawReader->GetDataSize() * 8) / 32;
}
//_____________________________________________________________________________
UInt_t AliT0RawReader::GetNextWord()
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

