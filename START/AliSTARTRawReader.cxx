#include "AliSTARTRawReader.h"
#include "AliSTARTRawData.h"
#include "AliSTARTdigit.h"
#include "AliBitPacking.h"
#include "TBits.h"

#include <Riostream.h>
#include "TMath.h"
#include "TH1F.h"
#include "TArrayI.h"
#include "AliLog.h"
 
ClassImp(AliSTARTRawReader)
  
  AliSTARTRawReader::AliSTARTRawReader (AliRawReader *rawReader, TTree* tree)
    :  TTask("STARTRawReader","read raw START data"),
       fDigits(NULL),
       fTree(tree),
       fRawReader(rawReader),
       fData(NULL),
       fPosition(0)
{
  //
// create an object to read STARTraw digits
  AliDebug(1,"Start ");
  if (fDigits == 0x0) fDigits = new AliSTARTdigit(); 
  fTree->Branch("START","AliSTARTdigit",&fDigits,405,1);
 
  fRawReader->Reset();
  fRawReader->Select("START");
 
 
}
 AliSTARTRawReader::~AliSTARTRawReader ()
{
  // 
}

Bool_t  AliSTARTRawReader::Next()
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

  UInt_t word;
  Int_t time=0,  itdc=0, ichannel=0; 
  Int_t numberOfWordsInTRM=0, iTRM=0;
  Int_t tdcTime, koef, meanTime, timeDiff ;
  Int_t allData[107];

  TArrayI *timeTDC1 = new TArrayI(24);
  TArrayI * chargeTDC1 = new TArrayI(24);
  TArrayI *timeTDC2 = new TArrayI(24);
  TArrayI *chargeTDC2 = new TArrayI(24);
   
  for ( Int_t k=0; k<107; k++)  allData[k]=0;
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
   word = GetNextWord();
  
   for (Int_t i=0; i<numberOfWordsInTRM; i++) {
     word = GetNextWord();
     tdcTime =  AliBitPacking::UnpackWord(word,31,31);   

     if ( tdcTime == 1)
       {
	 itdc=AliBitPacking::UnpackWord(word,24,27);
	 ichannel=AliBitPacking::UnpackWord(word,21,23);
	 time=AliBitPacking::UnpackWord(word,0,20);
	 koef = itdc*4 + ichannel/2;
	 allData[koef]=time;
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
	 koef = itdc*4 + ichannel/2;
	 allData[koef+54]=time;
       }
   }
      
   for (Int_t in=0; in<24; in++)
     {
       timeTDC1->AddAt(allData[in],in);
       timeTDC2->AddAt(allData[in+24],in);
       chargeTDC1->AddAt(allData[in+54],in);
       chargeTDC2->AddAt(allData[in+78],in);
     }      

   meanTime = allData[48];  // T0 !!!!!!
   timeDiff = allData[49];

   word = GetNextWord();
   word = GetNextWord();
   
   fDigits->SetTime(*timeTDC2);
   fDigits->SetADC(*chargeTDC1);
   
   fDigits->SetTimeAmp(*timeTDC1);
   fDigits->SetADCAmp(*chargeTDC2);

   fDigits->SetMeanTime(meanTime);
   fDigits->SetDiffTime(timeDiff);
   fTree->Fill();
   
   delete timeTDC1 ;
   delete chargeTDC1;
   delete timeTDC2 ;
   delete chargeTDC2;
   
   return kTRUE;
}
//_____________________________________________________________________________
/*
void AliSTARTRawReader::UnpackTime(Int_t outTime, Int_t outCh)
{
      UInt_t word=0;
      UInt_t unpackword=0;
    
      word = GetNextWord();
      unpackword=AliBitPacking::UnpackWord(word,0,12);
      outTime=unpackword;
      unpackword= AliBitPacking::UnpackWord(word,21,27);
      outCh=unpackword;  
 }
 */
//_____________________________________________________________________________
Int_t AliSTARTRawReader::GetPosition()
{
  // Sets the position in the
  // input stream
  if (((fRawReader->GetDataSize() * 8) % 32) != 0)
    AliFatal(Form("Incorrect raw data size ! %d words are found !",fRawReader->GetDataSize()));
  return (fRawReader->GetDataSize() * 8) / 32;
}
//_____________________________________________________________________________
UInt_t AliSTARTRawReader::GetNextWord()
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

