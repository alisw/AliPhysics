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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

//#include <>

ClassImp(AliTRDtestBeam)

const Long_t AliTRDtestBeam::file_head_size = 544; // ?
const Long_t AliTRDtestBeam::event_head_size = 68; //?
const Long_t AliTRDtestBeam::ldc_head_size = 68; //?
const Long_t AliTRDtestBeam::equip_head_size = 28; //
const Int_t AliTRDtestBeam::vme_in =1; //VME event in
const Int_t AliTRDtestBeam::sim_in =1; //Si-strips in

//typedef char byte;

//offsets in bytes
const Int_t AliTRDtestBeam::pos_run = 20; //run nr. (in file and event header)
const Int_t AliTRDtestBeam::pos_length = 0; //event/equip. length
const Int_t AliTRDtestBeam::pos_eqid = 8;  //equipment id.
const Int_t AliTRDtestBeam::pos_sioff = 12;  //Si data size offset (3 extra words!!!)      
     
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

  fDataStream = new ifstream(filename, ifstream::in | ifstream::binary );
  cout << fDataStream->is_open() << endl;
  //fHeaderIsRead = kTRUE;
  fHeaderIsRead = kTRUE;

  fFileHeader = new Char_t[file_head_size];
  fEventHeader = new Char_t[event_head_size];
  fEventData = new Char_t[fLimit];
}
 
//____________________________________________________________________________ 

Int_t AliTRDtestBeam::NextEvent() {
  
  Long_t data_size=0,ldc_off; //,ldc_id,ldc2_id;
  Long_t ldc_size,eq_id; //,ev_l2;
  Long_t event_nr,ev_l1;
  Long_t word;
  
  if ( !fHeaderIsRead ) {
    fDataStream->read(fFileHeader, file_head_size);
    if(fDataStream->fail()) {
      cerr << "Error reading file header! " << endl;	
      return false;
    }
    cout  << " Run nr.  " << Int(pos_run, fFileHeader) << endl;    
    fHeaderIsRead=kTRUE;
  }

  fDataStream->read(fEventHeader, event_head_size);
  if(fDataStream->fail()) {
    cerr << "End of file, Event " << fEventCount  << endl;	
    return false;
  }
  
  data_size = Int(pos_length, fEventHeader)-event_head_size; //?
  event_nr = Int((4+pos_run), fEventHeader); //ev.nr.
  //cout << " Event " << event_nr <<" size "<< data_size <<endl;
  
    if (event_nr <= fEventCount-1) { //watch-out ...event counter starts at 1?
      cout << fEventCount << " End of file?, Event " << fEventCount << endl;	
      return false;
    }
    //cout <<  "Run " << Int(pos_run, header)<< " , Event " <<event_nr <<endl;
    
    // enough space for data?
    if (fLimit < data_size) {
      delete[] fEventData;
      fEventData = new Char_t[data_size];
      fLimit = data_size;
    }
    
    fDataStream->read(fEventData, data_size);
    
    if(fDataStream->fail()) {
      cerr << "End of file, Event " << fEventCount; // << endl;	
	return false;
    }
    
    //cout  << " ...IDs (size) : ";
    
    ldc_off=0; // size of data from one DDL link
    
    for ( size_t k = 0; k < 2; k++ ) { // 2 LDCs (DDL & VME)
      
      ldc_size = Int(ldc_off+pos_length, fEventData); //
      //ldc_size1=(ldc_size-ldc_head_size);
      eq_id = Int(ldc_off+ldc_head_size+pos_eqid, fEventData);	    
      //cout  << eq_id <<" ("<<ldc_size<<") ";	    
      
      ev_l1 = Int((4+ldc_off+pos_run), fEventData); //ev.nr.
      if ( ev_l1 != event_nr ){
	//cerr << "Eq_id " <<eq_id<<" event nr. mismatch? " << event_nr <<" / "<< ev_l1 <<" ...LDC data size (header:68) " <<ldc_size<<endl;
      }
      
      if (eq_id == 1024) {  //DDL data
	fDdlOff = ldc_off; //+ldc_head_size+equip_head_size + 32;
	fDdlSize = ldc_size;
      }
      
      if (eq_id == 550) {  //Si-strip data (+QDC)
	//cout << "550" << endl;
	fSiOff=ldc_off+ldc_head_size+equip_head_size+pos_sioff;
	word = Int(fSiOff, fEventData);
	Short_t LenSi1 = (word >> 16) & 0xffff;
	Short_t LenSi2 = word & 0xffff;
	fQdcOff=fSiOff+4*(LenSi1+LenSi2+1)+equip_head_size+4; 
      } 
      else if (eq_id == 1182) {  //QDC first...
	//cout << "1182" << endl;
	fQdcOff=ldc_off+ldc_head_size+equip_head_size+pos_sioff;
	fSiOff=fQdcOff+equip_head_size+4;
      }
      
      ldc_off=ldc_size;
      
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
Int_t AliTRDtestBeam::DecodeSi() {
  
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

  int LenSiX = 640;

  int qmaxX; int amaxX;
  int qmaxY; int amaxY;
  
  qmaxX = 5;
  qmaxY = 5;
  amaxX = -1;
  amaxY = -1+LenSiX;
 
  for( int i = 0; i < GetNSi1(); i++ ) {
 
    if (fSi1Address[i] == 0) continue; // noise
   
    if (fSi1Address[i] < LenSiX ) {
      if( fSi1Charge[i] > qmaxX ) {
	qmaxX = fSi1Charge[i];
	amaxX = fSi1Address[i];
      }
    } else  {
      if( fSi1Charge[i] > qmaxY ) {
	qmaxY = fSi1Charge[i];
	amaxY = fSi1Address[i];
      }
    }
  }
  
  fX[0] = (float)(amaxX*0.05);  // [mm]
  fY[0] = (float)((amaxY-LenSiX)*0.05);
  fQx[0] = (float)qmaxX;
  fQy[0] = (float)qmaxY;
  
  // 
  qmaxX = 5;
  qmaxY = 5;
  amaxX = -1;
  amaxY = -1+LenSiX;

  for( int i = 0; i < GetNSi2(); i++ ) {
    
    if (fSi2Address[i] == 1279) continue; // noise
    if (fSi2Address[i] == 0) continue;    // noise
    
    if(fSi2Address[i] < LenSiX) {
      if( fSi2Charge[i] > qmaxX ) {
	qmaxX = fSi2Charge[i];
	amaxX = fSi2Address[i];
      }
    } else {
      if( fSi2Charge[i] > qmaxY ) {
	//if (fSi2Charge[i] > 50) cout << fSi2Charge[i] << " " << i << " " <<  fSi2Address[i] << endl;
	qmaxY = fSi2Charge[i];
	amaxY = fSi2Address[i];
      }
    }
  }
  
  fX[1] = (float)(amaxX*0.05);  // [mm]
  fY[1] = (float)((amaxY-LenSiX)*0.05);
  fQx[1] = (float)qmaxX;
  fQy[1] = (float)qmaxY;
  
  if (fQdcOff < 0) return 0;
 
  word=Int(fQdcOff, fEventData);
  fPb   = (Double_t)((word >> 16) & 0xFFF);
  fCher = (Double_t)((word ) & 0xFFF);

  //cout << fCher << " " << fPb << endl;
  return 1;
}
//____________________________________________________________________________ 
/**/
AliTRDRawStreamTB *AliTRDtestBeam::GetTRDrawStream() {
  
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
/**/
//____________________________________________________________________________ 

Int_t AliTRDtestBeam::Int(Int_t i, Char_t *start) {
  
  bool swap = kFALSE;

  if(swap) {
    char *q=(char*)(start+i); 
    char p[] = {q[3], q[2], q[1], q[0]};
    return *((int*) p);
  } else return *((int*)(start+i));
}

//____________________________________________________________________________ 
