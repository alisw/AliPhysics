#include "AliSTARTRawReader.h"
#include "AliSTARTRawData.h"
#include "AliRawReaderFile.h" 


#include <Riostream.h>
#include "TMath.h"
#include "TH1F.h"
#include "TArrayI.h"
#include "AliLog.h"
 
ClassImp(AliSTARTRawReader)

  AliSTARTRawReader::AliSTARTRawReader (): TTask("STARTRawReader","read raw data"),
					   fPMTId(-1),
					   fTimeTDC1(0),
					   fChargeADC1(0),
					   fTimeTDC2(0),
					   fChargeADC2(0)
{
  //
// create an object to read STARTraw digits
  AliDebug(1,"Start ");
  fTimeTDC1 = new TArrayI(24); 
  fChargeADC1 = new TArrayI(24); 
  fTimeTDC2 = new TArrayI(24); 
  fChargeADC2 = new TArrayI(24); 
 
}
 AliSTARTRawReader::~AliSTARTRawReader ()
{
  delete fTimeTDC1;
  delete fTimeTDC2;
  delete fChargeADC1 ; 
  delete fChargeADC2; 
 
}
//------------------------------------------------------------------------------------------------

UInt_t AliSTARTRawReader::UnpackWord(UInt_t packedWord, Int_t startBit, Int_t stopBit)
{
  //This method unpacks a words of StopBit-StartBit+1 bits starting from "StopBits"  
  UInt_t word;
  UInt_t offSet;
  Int_t length;
  length=stopBit-startBit+1;
  offSet=(UInt_t)TMath::Power(2,length)-1;
  offSet<<=startBit;
  word=packedWord&offSet;
  word>>=startBit;
  return word;
}
//---------------------------------------------------------------------------------------
Bool_t AliSTARTRawReader::NextThing( AliRawReader *fRawReader)
{
// read the next raw digit
// returns kFALSE if there is no digit left

  UInt_t word, unpackword; 
  UInt_t fADC, fTime;
  fRawReader->Select(13);
 
  if (!fRawReader->ReadNextInt(fData)) return kFALSE;
   

  Int_t size=fRawReader->GetDataSize();
  for (Int_t i=0; i<size/32; i++)
    {
      word=0;
      unpackword=0;
      fRawReader->ReadNextInt(word);
      unpackword=UnpackWord(word,0,5);
      fPMTId=unpackword;
      word=0;
      unpackword=0;

      fRawReader->ReadNextInt(word);
      unpackword=UnpackWord(word,8,31);
      fTime=unpackword;
      fTimeTDC1->AddAt(fTime,fPMTId);
      word=0;
      unpackword=0;

      fRawReader->ReadNextInt(word);
      unpackword=UnpackWord(word,0,5);
      fPMTId=unpackword;
      word=0;
      unpackword=0;

      fRawReader->ReadNextInt(word);
  
      unpackword=UnpackWord(word,8,31);
      fTime=unpackword;
      fTimeTDC2->AddAt(fTime,fPMTId);

      word=0;
      unpackword=0;
      fRawReader->ReadNextInt(word);
      unpackword=UnpackWord(word,0,5);
      fPMTId=unpackword;
      word=0;
      unpackword=0;

      fRawReader->ReadNextInt(word);
      unpackword= UnpackWord(word,8,31);
      fADC=unpackword;
      fChargeADC1 -> AddAt(fADC, fPMTId); 
 
      word=0;
      unpackword=0;

      fRawReader->ReadNextInt(word);
      unpackword=UnpackWord(word,0,5);
      fPMTId=unpackword;

      word=0;
      unpackword=0;
      fRawReader->ReadNextInt(word);
  
      unpackword= UnpackWord(word,8,31);
      fADC=unpackword;
      fChargeADC2 -> AddAt(fADC, fPMTId); 
      }
  return kTRUE;

}
 
//--------------------------------------------
void AliSTARTRawReader::GetTime (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTimeTDC1->At(i);
    }
}
//--------------------------------------------
//--------------------------------------------
void AliSTARTRawReader::GetADC (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fChargeADC1->At(i);
    }
}
