/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
Revision 0.01  2004/6/11 A.De Caro, S.B.Sellitto, R.Silvestri
        First implementation: global methods RawDataTOF
                                             GetDigits
                                             (Un)PackWord
*/
//
// This class contains the methods to create the Raw Data files
// for the TOF detector starting from the Digits.
// In this preliminary implementation, we defined the structure
// of the ALICE-TOF raw data starting from the current format
// for the TOF digits and the TOF raw data.
//

#include <stdlib.h>
#include <Riostream.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TMath.h>
#include "AliTOF.h"
#include "AliTOFGeometry.h"
#include "AliTOFdigit.h"
#include "AliTOFDDLRawData.h"
#include "AliRawDataHeader.h"

ClassImp(AliTOFDDLRawData)

//----------------------------------------------------------------------------------------
AliTOFDDLRawData::AliTOFDDLRawData()
{
  //Default constructor
  fIndex=-1;
  fVerbose=0;
}

//----------------------------------------------------------------------------------------

AliTOFDDLRawData::AliTOFDDLRawData(const AliTOFDDLRawData &source) : 
    TObject(source){
  //Copy Constructor
  this->fIndex=source.fIndex;
  this->fVerbose=source.fVerbose;
  return;
}

//---------------------------------------------------------------------------------------

AliTOFDDLRawData& AliTOFDDLRawData::operator=(const AliTOFDDLRawData &source){
  //Assigment operator
  this->fIndex=source.fIndex;
  this->fVerbose=source.fVerbose;
  return *this;
}

//---------------------------------------------------------------------------------------

void AliTOFDDLRawData::GetDigits(TClonesArray *TOFdigits,Int_t nDDL,UInt_t *buf)
{
  //This method packs the TOF digits in a proper 32 bits structure
  Int_t iDDL=(Int_t)((nDDL/4.-(Int_t)(nDDL/4.))*4);
  Int_t iSector=(Int_t)(nDDL/4.);
  Int_t iTRM=0;
  Int_t iTDC=0;
  Int_t iCH=-1;
  Int_t sector; 
  Int_t plate;
  Int_t strip;
  Int_t padx;
  Int_t padz;
  Int_t totCharge;
  Int_t timeOfFlight;
  Int_t error=0;
  Int_t eureka;
  UInt_t word;
  UInt_t baseWord;
  Int_t ndigits = TOFdigits->GetEntries();
  AliTOFdigit *digs;
  ofstream ftxt;
  if(!ndigits) 
    {
      Error("GetDigits", "No found TOF digits\n");      
      return;
    }

  if (fVerbose==2) ftxt.open("TOFdigits.txt",ios::app);
  for (Int_t digit=0;digit<ndigits;digit++) {
    digs = (AliTOFdigit*)TOFdigits->UncheckedAt(digit);
    sector=digs->GetSector(); // Sector Number (0-17)
    plate=digs->GetPlate();   // Plate Number (0-4)
    strip=digs->GetStrip();   // Strip Number (0-14/18/19)
    padx=digs->GetPadx();     // Pad Number in x direction (0-47)
    padz=digs->GetPadz();     // Pad Number in z direction (0-1)
    eureka=digs->GetTotPad(); // Global Pad Number inside a Sector
    totCharge = (Int_t)digs->GetAdc();
    timeOfFlight = (Int_t)digs->GetTdc();
    /*
    Int_t istriPlate=0;
    switch (plate)
      {
      case 0:
	break;
      case 1:
	istriPlate = AliTOFGeometry::NStripC();
	break;
      case 2:
	istriPlate = AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB();
	break;
      case 3:
	istriPlate = AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA();
	break;
      case 4:
	istriPlate = AliTOFGeometry::NStripC()+2*AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA();
	break;
      }

    eureka=2*padx+padz+AliTOFGeometry::NpadXStrip()*(strip+istriPlate);
    
    if (eureka!=digs->GetTotPad()) printf(" eureka = %d AND digs->GetTotPad() = %d",eureka,digs->GetTotPad());
    */
    if (sector!=iSector || (Int_t)((Float_t)eureka/AliTOF::NPadXTRM()/AliTOF::NTRM())!=iDDL) continue;
    
    if (fVerbose==2) ftxt <<" Sector: "<<sector<<" plate: "<<plate<<" strip "<<strip<<" padx "<<padx<<" padz "<<padz<<" eureka "<<eureka<<endl;
    
    iTRM = (Int_t)((Float_t)eureka/AliTOF::NPadXTRM() - AliTOF::NTRM()*iDDL);
    
    iTDC = (Int_t)(AliTOF::NTdc()* 
		   (
		    (Float_t)eureka/AliTOF::NPadXTRM() -
		    (Int_t)((Float_t)eureka/AliTOF::NPadXTRM())
		    )
		   );
    
    iCH  = (Int_t)(AliTOF::NCh() * 
		   (
		    (Float_t)eureka/AliTOF::NPadXTRM()*AliTOF::NTdc() - (Int_t)((Float_t)eureka/AliTOF::NPadXTRM()*AliTOF::NTdc()) -
		    (Int_t)((Float_t)eureka/AliTOF::NPadXTRM()*AliTOF::NTdc() - (Int_t)((Float_t)eureka/AliTOF::NPadXTRM()*AliTOF::NTdc()))
		    )
		   );
    
    if (fVerbose==2) ftxt << "DDL: "<<iDDL<<" Sector: "<<sector<<" TRM: "<<iTRM<<" TDC: "<<iTDC<<" Channel: "<<iCH<<" totCharge: "<<totCharge<<" tof: "<<timeOfFlight<<endl;
    
    baseWord=0;
    
    word=iTRM;
    PackWord(baseWord,word, 0, 3); // TRM ID
    word=iTDC;
    PackWord(baseWord,word, 4, 8); // TDC ID
    word=iCH;
    PackWord(baseWord,word, 9,11); // CH ID

    // temporary control
    if (totCharge<0) word=TMath::Abs(totCharge);
    else word=totCharge;
    PackWord(baseWord,word,12,31); // Charge (TOT)
    
    fIndex++;
    buf[fIndex]=baseWord;
    
    baseWord=0;
    
    word=error;
    PackWord(baseWord,word,0, 7); // Error flag
    word=timeOfFlight;
    PackWord(baseWord,word,8,31); // time-of-flight
    
    fIndex++;
    buf[fIndex]=baseWord;
    word=0;
    baseWord=0;
    
  }//end for
  
  if (fVerbose==2) ftxt.close();

  return;

}//end GetDigits

//-------------------------------------------------------------------------------------

void AliTOFDDLRawData::PackWord(UInt_t &BaseWord, UInt_t Word, Int_t StartBit, Int_t StopBit)
{
  //This method packs a word into the Baseword buffer starting form the "StartBit" 
  //and tacking StopBit-StartBit+1 bits
  UInt_t dummyWord,offSet;
  Int_t  length;
  UInt_t sum;
  //The BaseWord is being filled with 1 from StartBit to StopBit
  length=StopBit-StartBit+1;
  sum=(UInt_t)TMath::Power(2,length)-1;

  if(Word > sum){
    Error("PackWord", "Word to be filled is not within desired length\n"
	  "Word:%d Start bit:%d Stop Bit:%d",Word,StartBit,StopBit);
    return;
  }
  offSet=sum;
  offSet<<=StartBit;
  BaseWord=BaseWord|offSet;

  //The Word to be filled is shifted to the position StartBit
  //and the remaining  Left and Right bits are filled with 1
  sum=(UInt_t)TMath::Power(2,StartBit)-1;
  dummyWord=0xFFFFFFFF<<length;
  dummyWord +=Word;
  dummyWord<<=StartBit;
  dummyWord+=sum;
  BaseWord=BaseWord&dummyWord;

  return;
}

//------------------------------------------------------------------------------------------------

void AliTOFDDLRawData::UnpackWord(UInt_t PackedWord, Int_t StartBit, Int_t StopBit, UInt_t &Word)
{
  //This method unpacks a words of StopBit-StartBit+1 bits starting from "StopBits"  
  UInt_t offSet;
  Int_t length;
  length=StopBit-StartBit+1;
  offSet=(UInt_t)TMath::Power(2,length)-1;
  offSet<<=StartBit;
  Word=PackedWord&offSet;
  Word>>=StartBit;
  return;
}

//---------------------------------------------------------------------------------------

Int_t AliTOFDDLRawData::RawDataTOF(TBranch* branch){
  //This method creates the Raw data files for TOF detector
  const Int_t kSize=5000; //2*AliTOFGeometry::NpadXSector() 
                          //max number of digits per DDL file times 2
  UInt_t buf[kSize];
  UInt_t baseWord;
  UInt_t word;

  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  char fileName[15];
  ofstream outfile;         // logical name of the output file 
  AliRawDataHeader header;

  //loop over TOF DDL files
  for(Int_t i=0;i<72;i++){
    sprintf(fileName,"TOF_%d.ddl",i+kDDLOffset); //The name of the output file
#ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile.tellp();
    outfile.write((char*)(&header),sizeof(header));

    baseWord=0;
    word=i;
    PackWord(baseWord,word,0, 31); // Number of DDL file

    fIndex++;
    buf[fIndex]=baseWord;

    branch->GetEvent();

    //For each DDL file, buf contains the array of data words in Binary format
    //fIndex gives the number of 32 bits words in the buffer for each DDL file
    GetDigits(digits,i,buf);
    outfile.write((char *)buf,((fIndex+1)*sizeof(UInt_t)));

    for(Int_t ii=0;ii<(fIndex+1);ii++) buf[ii]=0;
    fIndex=-1;
    
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile.tellp();
    outfile.seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile.write((char*)(&header),sizeof(header));
    outfile.close();

  }//end for
  
  return 0;  
}

//-----------------------------------------------------------------------------------------
/*
void AliTOFDDLRawData::WriteChipHeader(Int_t ChipAddr,Int_t EventCnt,UInt_t &BaseWord)
{
  //This method writes a chip header 
  BaseWord=0;
  PackWord(BaseWord,ChipAddr,0,3);
  PackWord(BaseWord,EventCnt,4,10);
  PackWord(BaseWord,0x7,11,13);
  PackWord(BaseWord,0x1,14,15);
  return;
}//end WriteChipHeader
*/
//----------------------------------------------------------------------------------------
/*
void AliTOFDDLRawData::ReadChipHeader(Int_t &ChipAddr,Int_t &EventCnt,UInt_t BaseWord)
{
  //This method reads a chip header
  UInt_t temp=0;
  UnpackWord(BaseWord,0,3,temp);
  ChipAddr=(Int_t)temp;
  UnpackWord(BaseWord,4,10,temp);
  EventCnt=(Int_t)temp;
  if(fVerbose)
    Info("ReadChipHeader", "Chip:&d Half Stave module:%d",ChipAddr,EventCnt);
  return;
}//end ReadChipHeader
*/
//----------------------------------------------------------------------------------------
/*
void  AliTOFDDLRawData::WriteChipTrailer(UInt_t *buf,Int_t ChipHitCount,UInt_t &BaseWord)
{
  //This method writes a chip trailer
  //pixel fill word
  if((ChipHitCount%2)!=0){
    PackWord(BaseWord,0xFEDC,0,15);
  }
  PackWord(BaseWord,ChipHitCount,16,28);
  PackWord(BaseWord,0x0,30,31);
  fIndex++;
  buf[fIndex]=BaseWord;
  BaseWord=0;
  return;
}//end WriteChipTrailer
*/
//------------------------------------------------------------------------------------------
/*
void  AliTOFDDLRawData::ReadChipTrailer(Int_t &ChipHitCount,UInt_t BaseWord)
{
  //This method reads a chip trailer
  UInt_t temp=0;
  UnpackWord(BaseWord,16,28,temp);
  ChipHitCount=(Int_t)temp;
  return;
}//end ReadChipTrailer
*/
//------------------------------------------------------------------------------------------
