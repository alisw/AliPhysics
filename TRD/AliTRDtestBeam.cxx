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

#include "AliTRDtestBeam.h"

#include "AliTRDRawStreamTB.h"
#include "AliRawReaderMemory.h"

#include <iostream>
#include <fstream>

/*
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
*/

ClassImp(AliTRDtestBeam)

const Long_t AliTRDtestBeam::fgkFileHeadSize = 544; // ?
const Long_t AliTRDtestBeam::fgkEventHeadSize = 68; //?
const Long_t AliTRDtestBeam::fgkLdcHeadSize = 68; //?
const Long_t AliTRDtestBeam::fgkEquipHeadSize = 28; //
const Int_t AliTRDtestBeam::fgkVmeIn =1; //VME event in
const Int_t AliTRDtestBeam::fgkSimIn =1; //Si-strips in

//typedef char byte;

//offsets in bytes
const Int_t AliTRDtestBeam::fgkPosRun = 20; //run nr. (in file and event header)
const Int_t AliTRDtestBeam::fgkPosLength = 0; //event/equip. length
const Int_t AliTRDtestBeam::fgkEqId = 8;  //equipment id.
const Int_t AliTRDtestBeam::fgkPosSiOff = 12;  //Si data size offset (3 extra words!!!)      
     
using namespace std;

//____________________________________________________________________________ 
AliTRDtestBeam::AliTRDtestBeam() :
  fDataStream(0),
  fHeaderIsRead(0),
  fEventCount(0),
  fLimit(4), 
  fCurrent(0),
  fDdlOff(0),
  fSiOff(0),
  fQdcOff(0),
  fDdlSize(0),
  fFileHeader(0),
  fEventHeader(0),
  fEventData(0),
  fNSi1(0),
  fNSi2(0),
  fCher(0),
  fPb(0)
{
  //
  // Standard construction
  //

}
//____________________________________________________________________________ 
AliTRDtestBeam::AliTRDtestBeam(const char *filename) :
  fDataStream(0),
  fHeaderIsRead(0),
  fEventCount(0),
  fLimit(4), 
  fCurrent(0),
  fDdlOff(0),
  fSiOff(0),
  fQdcOff(0),
  fDdlSize(0),
  fFileHeader(0),
  fEventHeader(0),
  fEventData(0),
  fNSi1(0),
  fNSi2(0),
  fCher(0),
  fPb(0)
{
  //
  // AliTRDtestBeam constructor
  //

  fDataStream = new ifstream(filename, ifstream::in | ifstream::binary );
  cout << fDataStream->is_open() << endl;
  //fHeaderIsRead = kTRUE;
  fHeaderIsRead = kTRUE;

  fFileHeader = new Char_t[fgkFileHeadSize];
  fEventHeader = new Char_t[fgkEventHeadSize];
  fEventData = new Char_t[fLimit];

}

//____________________________________________________________________________
AliTRDtestBeam::AliTRDtestBeam(const AliTRDtestBeam &tb)
 :TObject(tb),
  fDataStream(0),
  fHeaderIsRead(0),
  fEventCount(0),
  fLimit(4),
  fCurrent(0),
  fDdlOff(0),
  fSiOff(0),
  fQdcOff(0),
  fDdlSize(0),
  fFileHeader(0),
  fEventHeader(0),
  fEventData(0),
  fNSi1(0),
  fNSi2(0),
  fCher(0),
  fPb(0)
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________ 
AliTRDtestBeam::~AliTRDtestBeam() 
{
  //
  // Destructor
  //

  if (fDataStream)  delete fDataStream;
  if (fEventHeader) delete fEventHeader;
  if (fFileHeader)  delete fFileHeader;
  if (fEventData)   delete fEventData;

}

//____________________________________________________________________________ 
Int_t AliTRDtestBeam::NextEvent() 
{
  //
  // Read the next event
  //
  
  Long_t dataSize=0,ldcOff; //,ldc_id,ldc2_id;
  Long_t ldcSize,eqId; //,ev_l2;
  Long_t eventNr,evL1;
  Long_t word;
  
  if ( !fHeaderIsRead ) {
    fDataStream->read(fFileHeader, fgkFileHeadSize);
    if(fDataStream->fail()) {
      cerr << "Error reading file header! " << endl;	
      return false;
    }
    cout  << " Run nr.  " << Int(fgkPosRun, fFileHeader) << endl;    
    fHeaderIsRead=kTRUE;
  }

  fDataStream->read(fEventHeader, fgkEventHeadSize);
  if(fDataStream->fail()) {
    cerr << "End of file, Event " << fEventCount  << endl;	
    return false;
  }
  
  dataSize = Int(fgkPosLength, fEventHeader)-fgkEventHeadSize; //?
  eventNr = Int((4+fgkPosRun), fEventHeader); //ev.nr.
  //cout << " Event " << eventNr <<" size "<< dataSize <<endl;
  
    if (eventNr <= fEventCount-1) { //watch-out ...event counter starts at 1?
      cout << fEventCount << " End of file?, Event " << fEventCount << endl;	
      return false;
    }
    //cout <<  "Run " << Int(fgkPosRun, header)<< " , Event " <<eventNr <<endl;
    
    // enough space for data?
    if (fLimit < dataSize) {
      delete[] fEventData;
      fEventData = new Char_t[dataSize];
      fLimit = dataSize;
    }
    
    fDataStream->read(fEventData, dataSize);
    
    if(fDataStream->fail()) {
      cerr << "End of file, Event " << fEventCount; // << endl;	
	return false;
    }
    
    //cout  << " ...IDs (size) : ";
    
    ldcOff=0; // size of data from one DDL link
    
    for ( size_t k = 0; k < 2; k++ ) { // 2 LDCs (DDL & VME)
      
      ldcSize = Int(ldcOff+fgkPosLength, fEventData); //
      //ldcSize1=(ldcSize-fgkLdcHeadSize);
      eqId = Int(ldcOff+fgkLdcHeadSize+fgkEqId, fEventData);	    
      //cout  << eqId <<" ("<<ldcSize<<") ";	    
      
      evL1 = Int((4+ldcOff+fgkPosRun), fEventData); //ev.nr.
      if ( evL1 != eventNr ){
	//cerr << "eqId " <<eqId<<" event nr. mismatch? " << eventNr <<" / "<< evL1 <<" ...LDC data size (header:68) " <<ldcSize<<endl;
      }
      
      if (eqId == 1024) {  //DDL data
	fDdlOff = ldcOff; //+fgkLdcHeadSize+fgkEquipHeadSize + 32;
	fDdlSize = ldcSize;
      }
      
      if (eqId == 550) {  //Si-strip data (+QDC)
	//cout << "550" << endl;
	fSiOff=ldcOff+fgkLdcHeadSize+fgkEquipHeadSize+fgkPosSiOff;
	word = Int(fSiOff, fEventData);
	Short_t lenSi1 = (word >> 16) & 0xffff;
	Short_t lenSi2 = word & 0xffff;
	fQdcOff=fSiOff+4*(lenSi1+lenSi2+1)+fgkEquipHeadSize+4; 
      } 
      else if (eqId == 1182) {  //QDC first...
	//cout << "1182" << endl;
	fQdcOff=ldcOff+fgkLdcHeadSize+fgkEquipHeadSize+fgkPosSiOff;
	fSiOff=fQdcOff+fgkEquipHeadSize+4;
      }
      
      ldcOff=ldcSize;
      
    }
    //cout << endl;

    //cout << "DDL = " << fDdlOff << endl;
    // cout << "Si  = " << fSiOff << endl;
    //cout << "QDC = " << fQdcOff << endl;
    
    DecodeSi();

    fEventCount++; //event counter
    return true;
}

//____________________________________________________________________________
Int_t AliTRDtestBeam::DecodeSi() 
{
  //
  // Decode the silicon detector
  //
  
  if (fSiOff < 0) return 0;
  
  // cout << "decoding Si data" << endl;

  Long_t word;
  
  word=Int(fSiOff, fEventData);
  fNSi1 = (word >> 16) & 0xffff;
  fNSi2 = word & 0xffff;
  
  Int_t cSi=fSiOff; //   
  for (int i = 0; i < fNSi1; i++) {
    fSi1Address[i] =  ( Int(cSi, fEventData) >> 12 ) & 0x7ff;
    fSi1Charge[i] = Int(cSi, fEventData)  & 0xfff;
    cSi+=4;
  }
    
  for (int i = 0; i < fNSi2; i++) {  //1,for Date!
    fSi2Address[i] =  ( Int(cSi, fEventData) >> 12 ) & 0x7ff;
    fSi2Charge[i] = Int(cSi, fEventData)  & 0xfff;
    cSi+=4;
  }  
  
  // reconstruction

  int aLenSiX = 640;

  int amaxX=0;
  int amaxY=0;

  Int_t q, a;  
  Int_t Nst1=0,Nst2=0;
  Int_t QclX=0,QclY=0, NclX=0,NclY=0, NstX=0,NstY=0;
  const Int_t Thr = 20;

  Nst1=0;
  NstX=0;
  NstY=0;
  NclX=0;
  NclY=0;
  QclX=0;
  QclY=0;
  
  for( int i = 0; i < GetNSi1(); i++ ) {
 
    if (fSi1Address[i] == 0) continue; // noise

    q = fSi1Charge[i];
    a = fSi1Address[i];

    if ( q > Thr ) 
    {
	if ( i > 0 && i < (GetNSi1()-1) ) {

	    if ( (a-fSi1Address[i+1]) == -1 &&
		 (a-fSi1Address[i-1]) == 1) 
	    {  
		Nst1++;	  
		if (a < aLenSiX) {
		    QclX = q+fSi1Charge[i+1]+fSi1Charge[i-1];
		    NclX++;
		    NstX+=3;
		    amaxX = a;
		}
		else {
		    QclY = q+fSi1Charge[i+1]+fSi1Charge[i-1];
		    NclY++;
		    NstY+=3;
		    amaxY = a;
		}
		i+=1;
	    }
	    else if ( (a-fSi1Address[i-1]) == 1)
	    {  
		Nst1++;	  
		if (a < aLenSiX) {
		    QclX = q+fSi1Charge[i-1];
		    NclX++;
		    NstX+=2;
		    amaxX = a;
		}
		else {
		    QclY = q+fSi1Charge[i-1];
		    NclY++;
		    NstY+=2;
		    amaxY = a;
		}
	    }
	    else if ( (a-fSi1Address[i+1]) == -1)
	    {  
		Nst1++;	  
		if (a < aLenSiX) {
		    QclX = q+fSi1Charge[i+1];
		    NclX++;
		    NstX+=2;
		    amaxX = a;
		}
		else {
		    QclY = q+fSi1Charge[i+1];
		    NclY++;
		    NstY+=2;
		    amaxY = a;
		}
		i+=1;
	    }
	}
    }
  }
  if (Nst1==2 && NstX<4 && NstY<4 ) {
      fX[0] = (float)(amaxX*0.05);  // [mm]
      fY[0] = (float)((amaxY-aLenSiX)*0.05);
      fQx[0] = (float)QclX;
      fQy[0] = (float)QclY;
  }
  else {
      fX[0] = -1.;
      fY[0] = -1.;
      fQx[0] = 0.;
      fQy[0] = 0.;
  }
  
  // ...and Si2

  Nst2=0;
  NstX=0;
  NstY=0;
  NclX=0;
  NclY=0;
  QclX=0;
  QclY=0;

  for( int i = 0; i < GetNSi2(); i++ ) {
    
    if (fSi2Address[i] == 1279) continue; // noise
    if (fSi2Address[i] == 0) continue;    // noise
    
    q = fSi2Charge[i];
    a = fSi2Address[i];

    if ( q > Thr/2 ) //...as Si2 has 1/2 gain! 
    {
	if ( i > 0 && i < (GetNSi2()-1) ) {

	    if ( (a-fSi2Address[i+1]) == -1 &&
		 (a-fSi2Address[i-1]) == 1) 
	    {  
		Nst2++;	  
		if (a < aLenSiX) {
		    QclX = q+fSi2Charge[i+1]+fSi2Charge[i-1];
		    NclX++;
		    NstX+=3;
		    amaxX = a;
		}
		else {
		    QclY = q+fSi2Charge[i+1]+fSi2Charge[i-1];
		    NclY++;
		    NstY+=3;
		    amaxY = a;
		}
		i+=1;
	    }
	    else if ( (a-fSi2Address[i-1]) == 1)
	    {  
		Nst2++;	  
		if (a < aLenSiX) {
		    QclX = q+fSi2Charge[i-1];
		    NclX++;
		    NstX+=2;
		    amaxX = a;
		}
		else {
		    QclY = q+fSi2Charge[i-1];
		    NclY++;
		    NstY+=2;
		    amaxY = a;
		}
	    }
	    else if ( (a-fSi2Address[i+1]) == -1)
	    {  
		Nst2++;	  
		if (a < aLenSiX) {
		    QclX = q+fSi2Charge[i+1];
		    NclX++;
		    NstX+=2;
		    amaxX = a;
		}
		else {
		    QclY = q+fSi2Charge[i+1];
		    NclY++;
		    NstY+=2;
		    amaxY = a;
		}
		i+=1;
	    }
	}
    }
  }
  
  if (Nst2==2 && NstX<4 && NstY<4 ) {
      fX[1] = (float)(amaxX*0.05);  // [mm]
      fY[1] = (float)((amaxY-aLenSiX)*0.05);
      fQx[1] = (float)QclX;
      fQy[1] = (float)QclY;
  }
  else {
      fX[1] = -1.;
      fY[1] = -1.;
      fQx[1] = 0.;
      fQy[1] = 0.;
  }
  
  if (fQdcOff < 0) return 0;
 
  word=Int(fQdcOff, fEventData);
  fPb   = (Double_t)((word >> 16) & 0xFFF);
  fCher = (Double_t)((word ) & 0xFFF);

  //cout << fCher << " " << fPb << endl;
  return 1;

}
//____________________________________________________________________________ 
AliTRDRawStreamTB *AliTRDtestBeam::GetTRDrawStream() 
{
  //
  // Get the TRD raw stream
  //
  
  // needs AliTRDRawStreamTB  
  //cout << "Chamber reader:" << (Int_t)(fEventData+fDdlOff) << " " << fDdlSize << endl;
  //int ifout = open("dump.dat", O_WRONLY | O_TRUNC | O_CREAT);
  //write(ifout, (void*)(fEventData+fDdlOff+16), fDdlSize);
  //close(ifout);

  AliRawReaderMemory *reader = new AliRawReaderMemory((UChar_t*)(fEventData+fDdlOff), (UInt_t)fDdlSize);
  reader->SetEquipmentID(1024);
  reader->ReadHeader();
  AliTRDRawStreamTB::RawBufferMissAligned(kTRUE);
  AliTRDRawStreamTB::SupressWarnings(kTRUE);
 
  AliTRDRawStreamTB *tb = new AliTRDRawStreamTB(reader); 
  tb->Init();
  return tb;
  /*
    return 

    AliRawReaderMemory *rmem = data->GetRawReader();
    rmem->ReadHeader();
    
    AliTRDRawStreamTB tb(rmem);
    tb.Init();
    AliTRDRawStreamTB::SupressWarnings(kTRUE);
    
  */
}

//____________________________________________________________________________ 
Int_t AliTRDtestBeam::Int(Int_t i, Char_t *start) 
{
  //
  // ?????
  //
  
  bool swap = kFALSE;

  if(swap) {
    char *q=(char*)(start+i); 
    char p[] = {q[3], q[2], q[1], q[0]};
    return *((int*) p);
  } else return *((int*)(start+i));

}

//____________________________________________________________________________ 
