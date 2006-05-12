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
       fRawReader(rawReader)
{
  //
// create an object to read STARTraw digits
  AliDebug(1,"Start ");
  if (fDigits == 0x0) fDigits = new AliSTARTdigit(); 
  fTree->Branch("START","AliSTARTdigit",&fDigits,405,1);
 
  fRawReader->Reset();
  fRawReader->Select(13,0,1);
 
 
}
 AliSTARTRawReader::~AliSTARTRawReader ()
{
  // 
}

Bool_t  AliSTARTRawReader::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left
  AliBitPacking *pack ;

  UInt_t word, unpackword;
  Int_t time,  pmt;
  TArrayI *timeTDC1 = new TArrayI(24);
  TArrayI * chargeTDC1 = new TArrayI(24);
  TArrayI *timeTDC2 = new TArrayI(24);
  TArrayI *chargeTDC2 = new TArrayI(24);

  do {
    if (!fRawReader->ReadNextData(fData)) return kFALSE;
  } while (fRawReader->GetDataSize() == 0);
  
  fPosition = GetPosition();
 //  trigger 
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;

   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;

   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
    unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;

   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;

   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;


//multiplicity
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;
   fDigits->SetSumMult(time);   

   // best time differece  
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;

   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;
   fDigits->SetDiffTime(time);   
  // best time left 
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;
   fDigits->SetTimeBestLeft(time);   
   
  // Best time right &left  
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
 
   word=0;
   unpackword=0;
   
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;
   fDigits->SetTimeBestRight(time);  
  // mean 
   word=0;
   unpackword=0;
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,8,31);
   time=unpackword;
 
   word=0;
   unpackword=0;
   
   word = GetNextWord();
   unpackword=pack->UnpackWord(word,0,8);
   pmt=unpackword;
   fDigits->SetMeanTime(time);  

   for (Int_t i=0; i<24; i++)   //QTC amplified
     {
      word=0;
      unpackword=0;
    
      word = GetNextWord();
      unpackword=pack->UnpackWord(word,8,31);
      time=unpackword;
      word=0;
      unpackword=0;
      word = GetNextWord();
      unpackword= pack->UnpackWord(word,0,8);
      pmt=unpackword;  
      chargeTDC2->AddAt(time,pmt-72);
     }

  for (Int_t i=0; i<24; i++)
    {
      //  QTC
      word=0;
      unpackword=0;
      word = GetNextWord();
      unpackword=pack->UnpackWord(word,8,31);
      time=unpackword;
      word=0;
      unpackword=0;
      word = GetNextWord();
      unpackword= pack->UnpackWord(word,0,8);
      pmt=unpackword; //
      chargeTDC1->AddAt(time,pmt-48);
    }
  
  for (Int_t i=0; i<24; i++) //time CFD
    {
      word=0;
      unpackword=0;
      word = GetNextWord();
      unpackword=pack->UnpackWord(word,8,31);
      time=unpackword;
      word=0;
      unpackword=0;
      word = GetNextWord();
       unpackword=pack->UnpackWord(word,0,8);
      pmt=unpackword;
      timeTDC2->AddAt(time,pmt-24);
    } 

  
  for (Int_t i=0; i<24; i++) //time LED
    {
      word=0;
      unpackword=0;
      word = GetNextWord();
      unpackword=pack->UnpackWord(word,8,31);
      time=unpackword; 
      
      word=0;
      unpackword=0;
      word = GetNextWord();
      unpackword=pack->UnpackWord(word,0,8);
      pmt=unpackword;
      timeTDC1->AddAt(time,pmt);
    }
 
  
  fDigits->SetTime(*timeTDC2);
  fDigits->SetADC(*chargeTDC1);
  
  fDigits->SetTimeAmp(*timeTDC1);
  fDigits->SetADCAmp(*chargeTDC2);
    fTree->Fill();

    delete timeTDC1 ;
    delete chargeTDC1;
    delete timeTDC2 ;
    delete chargeTDC2;
  
    return kTRUE;
}
 
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

   fPosition--;
  Int_t iBit = fPosition * 32;
  Int_t iByte = iBit / 8;

  UInt_t word = 0;
  word  = fData[iByte+3]<<24;
  word |= fData[iByte+2]<<16;
  word |= fData[iByte+1]<<8;
  word |= fData[iByte];
// fPosition--;
  return word;

}

