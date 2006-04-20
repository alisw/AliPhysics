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
Revision 0.02  2005/7/25 A.De Caro
        Update number of bits allocated for time-of-flight
               and 'charge' measurements

Revision 0.01  2004/6/11 A.De Caro, S.B.Sellitto, R.Silvestri
        First implementation: global methods RawDataTOF
                                             GetDigits
*/

//////////////////////////////////////////////////////////////////
//
// This class contains the methods to create the Raw Data files
// for the TOF detector starting from the Digits.
// In this preliminary implementation, we defined the structure
// of the ALICE-TOF raw data starting from the current format
// for the TOF digits and the TOF raw data.
//
//////////////////////////////////////////////////////////////////

#include "Riostream.h"

#include "TBranch.h"
#include "TClonesArray.h"
#include "TMath.h"

#include "AliBitPacking.h"
#include "AliLog.h"
#include "AliRawDataHeader.h"

#include "AliTOFDDLRawData.h"
#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"
#include "AliTOFRawStream.h"

ClassImp(AliTOFDDLRawData)

//----------------------------------------------------------------------------------------
AliTOFDDLRawData::AliTOFDDLRawData()
{
  //Default constructor
  fIndex=-1;
  fVerbose=0;
  fTOFgeometry = 0;
}

//----------------------------------------------------------------------------------------
AliTOFDDLRawData::AliTOFDDLRawData(AliTOFGeometry *tofGeom)
{
  //Constructor
  fIndex=-1;
  fVerbose=0;
  fTOFgeometry = tofGeom;
}

//----------------------------------------------------------------------------------------

AliTOFDDLRawData::AliTOFDDLRawData(const AliTOFDDLRawData &source) : 
    TObject(source){
  //Copy Constructor
  this->fIndex=source.fIndex;
  this->fVerbose=source.fVerbose;
  this->fTOFgeometry=source.fTOFgeometry;
  return;
}

//---------------------------------------------------------------------------------------

AliTOFDDLRawData& AliTOFDDLRawData::operator=(const AliTOFDDLRawData &source){
  //Assigment operator
  this->fIndex=source.fIndex;
  this->fVerbose=source.fVerbose;
  this->fTOFgeometry=source.fTOFgeometry;
  return *this;
}

//---------------------------------------------------------------------------------------

void AliTOFDDLRawData::GetDigits(TClonesArray *TOFdigits,Int_t nDDL,UInt_t *buf)
{

  //This method packs the TOF digits in a proper 32 bits structure

  Int_t iDDL    = nDDL%AliTOFGeometry::NDDL();
  Int_t iSector = (Int_t)((Float_t)nDDL/AliTOFGeometry::NDDL());
  Int_t iTRM = 0;
  Int_t iTDC = 0;
  Int_t iCH  =-1;
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
      AliError("No found TOF digits");
      return;
    }

  if (fVerbose==2) ftxt.open("TOFdigits.txt",ios::app);

  for (Int_t digit=0;digit<ndigits;digit++) {
    digs = (AliTOFdigit*)TOFdigits->UncheckedAt(digit);
    sector = digs->GetSector(); // Sector Number (0-17)
    plate  = digs->GetPlate();  // Plate Number (0-4)
    strip  = digs->GetStrip();  // Strip Number (0-14/18/19)
    padx   = digs->GetPadx();   // Pad Number in x direction (0-47)
    padz   = digs->GetPadz();   // Pad Number in z direction (0-1)
    eureka = digs->GetTotPad(fTOFgeometry); // Global Pad Number inside a Sector
    //   totCharge = (Int_t)digs->GetAdc(); //Use realistic ToT, for Standard  production with no miscalibration/Slewing it == fAdC in digit (see AliTOFDigitizer)
    totCharge = (Int_t)digs->GetToT();
    timeOfFlight = (Int_t)digs->GetTdc();

    if (sector!=iSector || (Int_t)((Float_t)eureka/AliTOFGeometry::NPadXTRM()/AliTOFGeometry::NTRM())!=iDDL) continue;
    
    if (fVerbose==2) ftxt << " Sector: " << sector << " Plate: " << plate << " Strip: " << strip << " PadZ: " << padz << " PadX: " << padx << " totPadNumber " << eureka << endl;
    
    iTRM = (Int_t)((Float_t)eureka/AliTOFGeometry::NPadXTRM()) - iDDL*AliTOFGeometry::NTRM();

    
    iTDC = (Int_t)((Float_t)eureka/AliTOFGeometry::NCh()) - (iDDL*AliTOFGeometry::NTRM() + iTRM) * AliTOFGeometry::NTdc();
    /*
    iTDC = (Int_t)(AliTOFGeometry::NTdc()* 
		   (
		    (Float_t)eureka/AliTOFGeometry::NPadXTRM() -
		    (Int_t)((Float_t)eureka/AliTOFGeometry::NPadXTRM())
		    )
		   );
    */

    iCH  = eureka - ((iDDL*AliTOFGeometry::NTRM() + iTRM) * AliTOFGeometry::NTdc() + iTDC) * AliTOFGeometry::NCh();
    /*
    iCH  = (Int_t)(AliTOFGeometry::NCh() * 
		   (
		    (Float_t)eureka/AliTOFGeometry::NPadXTRM()*AliTOFGeometry::NTdc() - (Int_t)((Float_t)eureka/AliTOFGeometry::NPadXTRM()*AliTOFGeometry::NTdc()) -
		    (Int_t)((Float_t)eureka/AliTOFGeometry::NPadXTRM()*AliTOFGeometry::NTdc() - (Int_t)((Float_t)eureka/AliTOFGeometry::NPadXTRM()*AliTOFGeometry::NTdc()))
		    )
		   );
    */

    if (fVerbose==2) ftxt << "DDL: "<<nDDL<<" TRM: "<<iTRM<<" TDC: "<<iTDC<<" Channel: "<<iCH<<" totCharge: "<<totCharge<<" tof: "<<timeOfFlight<<endl;

    AliDebug(2,Form("%2i %2i %2i %2i   %2i %2i %2i %2i %2i %7i %8i",nDDL,iTRM,iTDC,iCH,sector,plate,strip,padz,padx,totCharge,timeOfFlight));
    
    baseWord=0;
    
    word=iTRM;
    AliBitPacking::PackWord(word,baseWord, 0, 3); // TRM ID
    word=iTDC;
    AliBitPacking::PackWord(word,baseWord, 4, 8); // TDC ID
    word=iCH;
    AliBitPacking::PackWord(word,baseWord, 9,11); // CH ID

    // temporary control
    if (totCharge<0) word=TMath::Abs(totCharge);
    else word=totCharge;
    AliBitPacking::PackWord(word,baseWord,12,31); // Charge (TOT) // v0.01
    //AliBitPacking::PackWord(word,baseWord,12,19); // Charge (TOT) // v0.02
    //AliBitPacking::PackWord(0,baseWord,20,31); // v0.02

    fIndex++;
    buf[fIndex]=baseWord;
    
    baseWord=0;
    
    word=error;
    AliBitPacking::PackWord(word,baseWord,0, 7); // Error flag
    word=timeOfFlight;
    AliBitPacking::PackWord(word,baseWord,8,31); // time-of-flight // v0.01
    //AliBitPacking::PackWord(word,baseWord,8,19); // time-of-flight // v0.02
    //AliBitPacking::PackWord(0,baseWord,20,30); // v0.02
    //AliBitPacking::PackWord(1,baseWord,31,31); // v0.02
    
    fIndex++;
    buf[fIndex]=baseWord;
    word=0;
    baseWord=0;
    
  }//end for
  
  if (fVerbose==2) ftxt.close();

  return;

}//end GetDigits

//---------------------------------------------------------------------------------------

Int_t AliTOFDDLRawData::RawDataTOF(TBranch* branch){
  //
  // This method creates the Raw data files for TOF detector
  //

  const Int_t kSize = 5000; //max number of digits per DDL file times 2

  UInt_t buf[kSize];
  //UInt_t baseWord; // v0.01
  //UInt_t word; // v0.01

  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  char fileName[15];
  ofstream outfile;         // logical name of the output file 
  AliRawDataHeader header;
  UInt_t sizeRawData = 0;

  //loop over TOF DDL files
  for(Int_t i = 0; i<AliTOFGeometry::NDDL()*AliTOFGeometry::NSectors(); i++){

    sprintf(fileName,"TOF_%d.ddl",i+AliTOFRawStream::kDDLOffset); //The name of the output file
#ifndef __DECCXX
    outfile.open(fileName,ios::binary);
#else
    outfile.open(fileName);
#endif

    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile.tellp();
    outfile.write((char*)(&header),sizeof(header));

    /*
    // v0.01
    baseWord=0;
    word=i;
    //AliBitPacking::PackWord(word,baseWord,0, 31); // Number of DDL file
    AliBitPacking::PackWord(word,baseWord,0, 6); // Number of DDL file
    AliBitPacking::PackWord(0,baseWord,7,31);

    fIndex++;
    buf[fIndex]=baseWord;
    */

    branch->GetEvent();

    //For each DDL file, buf contains the array of data words in Binary format
    //fIndex gives the number of 32 bits words in the buffer for each DDL file
    GetDigits(digits,i,buf);
    outfile.write((char *)buf,((fIndex+1)*sizeof(UInt_t)));

    for(Int_t ii=0;ii<(fIndex+1);ii++) buf[ii]=0;
    fIndex=-1;
    
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile.tellp();
    sizeRawData = currentFilePosition - dataHeaderPosition - sizeof(header);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile.seekp(dataHeaderPosition);
    outfile.write((char*)(&header),sizeof(header));
    outfile.seekp(currentFilePosition);

    outfile.close();

  }//end for
  
  return 0;  
}

//---------------------------------------------------------------------------
