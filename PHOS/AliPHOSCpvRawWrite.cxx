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

/* $Id: $ */

/* History:
 *
 * $Log$
 */

//_________________________________________________________________________
//  Create a raw data stream for the CPV detector
//  Input:  AliPHOSDigit or TClonesArray of AliPHOSDigits
//  Output: AliFstream, a raw data stream in DDL format
//  Number of CPV modules          :  5
//  Number of DDL per CPV module   :  1
//  Number of rows per DDL         : 16
//  Number of Dilogic chips per row: 10
//  Number of pads per Dilogic     : 48 (nX,nZ)=(8,6)
//--
//  Author: Yuri Kharlov
//  Reimplemented from AliHMPIDRawStream
//  6 July 2008
//_________________________________________________________________________

// --- ROOT system ---

#include "TObject.h"
#include "TClonesArray.h"

// --- Standard library ---
#include <assert.h>

// --- AliRoot header files ---
#include "AliBitPacking.h"
#include "AliDAQ.h"
#include "AliFstream.h"
#include "AliRawDataHeaderSim.h"
#include "AliPHOSCpvRawWrite.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGeometry.h"

ClassImp(AliPHOSCpvRawWrite)

//____________________________________________________________________________
AliPHOSCpvRawWrite::AliPHOSCpvRawWrite() : 
  TObject(),
  fNDDL(5),
  fNRow(16),
  fNDilogic(10),
  fNPad(48)
{
  // ctor
}

//____________________________________________________________________________
AliPHOSCpvRawWrite::~AliPHOSCpvRawWrite()
{
  // dtor
}

//____________________________________________________________________________
void AliPHOSCpvRawWrite::WriteRaw(const TObjArray *pDigAll)
{
// Write a list of digits for a given chamber in raw data stream
// Arguments: pDigAll- list of digits as TObjArray of 5 TClonesArray's
// Returns: none

  AliPHOSGeometry *geom;
  if (!(geom = AliPHOSGeometry::GetInstance())) 
        geom = AliPHOSGeometry::GetInstance("IHEP","");

  Int_t  ddl,row,dilogic,address; //32-bits data word 
  Int_t  cntLpad;
  Int_t  cntLrow;
  Int_t  cntL=0;        //data words counters for DDLs
  Int_t  cntLeoe;
  UInt_t posL;
  UInt_t cntLseg;
  UInt_t cntwInLseg=0;
  Int_t  cntLdig=0;
  
  UInt_t posLmarker;
  Int_t  digcnt=0;

  Int_t isDigThere[5][16][10][48];
  
  Int_t relId[4];

  for(Int_t iCh=0;iCh<fNDDL;iCh++){//chambers loop
    cntL=0;
    for(Int_t irow=0; irow<fNRow; irow++){
      for(Int_t idil=0; idil<fNDilogic; idil++){
	for(Int_t ipad=0; ipad<fNPad; ipad++){
	  isDigThere[iCh][irow][idil][ipad]=-1;
	}
      }
    }
    
    AliFstream* ddlL;                                 //output stream
    
    AliRawDataHeaderSim header; header.SetAttribute(0);  //empty DDL header
    
    ddlL = new AliFstream(AliDAQ::DdlFileName("CPV",iCh));

    //write dummy header as place holder, actual header
    //will be rewritten later when total size of DDL is known

    ddlL->WriteBuffer((char*)&header,sizeof(header));
    
    UInt_t w32=0;                 //32 bits data word 
    digcnt=0;
    
    TClonesArray *pDigCh=(TClonesArray *)pDigAll->At(iCh); //list of digits for current chamber 
   
    for(Int_t iDig=0;iDig<pDigCh->GetEntriesFast();iDig++){//digits loop
      AliPHOSDigit *pDig1=(AliPHOSDigit*)pDigCh->At(iDig);
      HWaddress(pDig1,w32,ddl,row,dilogic,address);
      isDigThere[ddl][row][dilogic][address]=iDig;
    }  
    
    Int_t cntRrow, cntLrow, cntLseg, cntLeoe, cntL, cntwInLseg;
    for(Int_t row = 0; row < fNRow; row++){ //AliHMPIDRawStream::kNRows=25!
      cntRrow=0;cntLrow=0;cntLseg=0,cntLeoe=0;
      posLmarker=ddlL->Tellp(); WriteRowMarker(ddlL,(UInt_t)1);   cntL++; cntLrow++; cntwInLseg++;
      for(Int_t dil = 0; dil < fNDilogic; dil++){
	cntLpad=0;
        for(Int_t pad = 0; pad < fNPad; pad++){
	  if (isDigThere[iCh][row][dil][pad]!=-1) {
	    AliPHOSDigit *pDig=(AliPHOSDigit*)pDigCh->At(isDigThere[iCh][row][dil][pad]);             
	    HWaddress(pDig,w32,iCh,row,dil,pad);
	    if(pDig->GetAmp() < 0 ) continue;  //We can turn of the zero sup for pedestal simulation
	    ddlL->WriteBuffer((char*)&w32,sizeof(w32));   cntL++; cntLpad++; cntLrow++;  cntLdig++; cntwInLseg++;
          }//isDig
	}//pad
        WriteEoE(ddlL,row,dil,cntLpad); cntL++;  cntLrow++;    cntLeoe++;   cntwInLseg++; // write EoE markers
      }//dil
      if(row%8==0){  // Why %8 ???
        WriteSegMarker(ddlL,row,cntwInLseg); cntL++;  cntLseg++; cntwInLseg=0;
      }
      // Find the marker position write and  go back to the actual position to continue writing                    
      posL=ddlL->Tellp();   ddlL->Seekp(posLmarker);    WriteRowMarker(ddlL,(UInt_t)(cntLrow-1)); ddlL->Seekp(posL);
    }//row
    //rewrite header with size set to
    header.fSize=sizeof(header)+cntL*sizeof(w32); ddlL->Seekp(0); ddlL->WriteBuffer((char*)&header,sizeof(header)); delete ddlL;
    
  }//chambers loop
}//WriteRaw()

//____________________________________________________________________________
void AliPHOSCpvRawWrite::HWaddress(const AliPHOSDigit *digit, 
				   UInt_t &w32, Int_t &ddl, Int_t &row, Int_t &dilogic, Int_t &address)
{
// Convert AliPHOSDigit to raw word format
// Arguments: w32,ddl,row,dilogic,address where to write the results

  AliPHOSGeometry *geom;
  if (!(geom = AliPHOSGeometry::GetInstance())) 
        geom = AliPHOSGeometry::GetInstance("IHEP","");

  Int_t relid[4] ;
  geom->AbsToRelNumbering(digit->GetId(), relid) ;
  if (relid[0]<geom->GetNModules()) Error("HWaddress","Digit does not belong to CPV!");

  Int_t module = relid[0];
  Int_t cellX  = relid[2];
  Int_t cellZ  = relid[3];
  UInt_t charge = (UInt_t)digit->GetAmp();

  ddl     = module-1;                      // DDL# 0..4
  row     = (cellX-1)/8;                   // row# 0..16
  dilogic = (cellZ-1)/6;                   // Dilogic# 0..10
  address = (cellX-1)%6 + 6*((cellZ-1)%8); // Address 0..47
  
  w32=0;
  AliBitPacking::PackWord(charge ,w32, 0,11); // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq  Qdc          bits (00..11) counts (0..4095)
  assert(0<=address && address<fNPad);
  AliBitPacking::PackWord(address,w32,12,17); // 3322 2222 2222 1111 1111 1000 0000 0000  DILOGIC addr bits (12..17) counts (0..47)
  assert(0<=dilogic && dilogic<fNDilogic);
  AliBitPacking::PackWord(dilogic,w32,18,21); // 1098 7654 3210 9876 5432 1098 7654 3210  DILOGIC No.  bits (18..21) counts (1..10)
  assert(0<=row && row<fNRow);
  AliBitPacking::PackWord(row ,w32,22,26); //                                          Row No.   bits (22..26) counts (1..24)
  AliBitPacking::PackWord((UInt_t)0, w32,27,27);  //To make sure set the 27th bit to Zero so we can distinguish it from the EoE
}

//____________________________________________________________________________
void AliPHOSCpvRawWrite::WriteEoE(AliFstream *ddl,UInt_t row,UInt_t dil,UInt_t wordCnt  )
{
  //Writes the EoE word from real data and pedestals into the ddl stream
  //Arguments:  ddl stream, row number, dilogic number and the number of words before the EoE
  //Returns:   nothing
  UInt_t e=1;
  UInt_t w32=0;
  assert(0<=row&&row<=fNRow);
  AliBitPacking::PackWord((UInt_t)row     ,w32,22,26);  // row number (1...24)
  assert(0<=dil&&dil<=fNDilogic);
  AliBitPacking::PackWord((UInt_t)dil     ,w32,18,21);  // DILOGIC number (1...10)
  AliBitPacking::PackWord(          e     ,w32, 7,17);  // event number -- not used
  AliBitPacking::PackWord((UInt_t)wordCnt ,w32, 0, 6);  // word counter (0...47)
  AliBitPacking::PackWord((UInt_t)1       ,w32,27,27);  // bit 27 is always 1 by definition of EoE
  ddl->WriteBuffer((char*)&w32,sizeof(w32));      
}

//____________________________________________________________________________
void AliPHOSCpvRawWrite::WriteRowMarker(AliFstream *ddl,UInt_t size)
{
  //Writes the row marker for real data and pedestal into the ddl stream
  //Arguments: ddl stream and the size of the block of the given row, the siye is at least the 10 EoE words!
  //Returns:   nothing

  UInt_t w32=0;
  UInt_t marker=13992;                         //for pedestal=12968  ==  32a8 for zero suppressed 36a8
  AliBitPacking::PackWord(size,  w32, 16,31);  //number of roaw written after row marker (digits and EoE)
  AliBitPacking::PackWord(marker,w32,0,15);    //the marker word
  ddl->WriteBuffer((char*)&w32,sizeof(w32));              
}

//____________________________________________________________________________
void AliPHOSCpvRawWrite::WriteSegMarker(AliFstream *ddl,UInt_t row, Int_t nwInSeg)
{
  //Writes the segment marker (after 8 rows) into the ddl stream
  //Arguments: ddl stream and the segment: row 8 -> 0x5800, row 16 -> 5801, row 24 -> 5802 for pedestal
  //Returns:   nothing
  UInt_t w32=0;
  
  //Segment marker: 2736 == ab0
  AliBitPacking::PackWord((UInt_t)2736   ,w32,20,31); //ab0 the segment marker word
  AliBitPacking::PackWord((UInt_t)nwInSeg,w32, 8,19); //number of words in the segment
  AliBitPacking::PackWord((UInt_t)(row/8),w32, 0, 7); //segment 0,1,2    
  ddl->WriteBuffer((char*)&w32,sizeof(w32)); 
  //Printf("Segment word created is: %x",w32);
}
