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

/* $Id$ */


//This class contains all the necessary methods to create the Raw Data
//files (slides) for the ITS data challenges for:
//SPD 
//SDD
//SSD

#include <stdlib.h>
//#include <Riostream.h>
#include <TClonesArray.h>
#include <TTree.h>
#include "AliITSdigit.h"
#include "AliITSDDLRawData.h"
#include "AliRawDataHeaderSim.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSDDLModuleMapSDD.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSIntMap.h"
#include "AliBitPacking.h"
#include "AliDAQ.h"
#include "AliFstream.h"
#include "AliITSFOSignalsSPD.h"

ClassImp(AliITSDDLRawData)

////////////////////////////////////////////////////////////////////////////////////////
AliITSDDLRawData::AliITSDDLRawData():
fVerbose(0),
fIndex(-1),
fHalfStaveModule(-1),
fSDDRawFormat(7){
  //Default constructor

}

////////////////////////////////////////////////////////////////////////////////////////

AliITSDDLRawData::AliITSDDLRawData(const AliITSDDLRawData &source) : 
    TObject(source),
fVerbose(source.fVerbose),
fIndex(source.fIndex),
fHalfStaveModule(source.fHalfStaveModule),
fSDDRawFormat(source.fSDDRawFormat){
  //Copy Constructor
}

////////////////////////////////////////////////////////////////////////////////////////

AliITSDDLRawData& AliITSDDLRawData::operator=(const AliITSDDLRawData &source){
  //Assigment operator
  if(this==&source) return *this;
  fIndex=source.fIndex;
  fHalfStaveModule=source.fHalfStaveModule;
  fVerbose=source.fVerbose;
  fSDDRawFormat=source.fSDDRawFormat;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////
//STRIP 
//

void AliITSDDLRawData::GetDigitsSSD(TClonesArray *ITSdigits,Int_t mod,Int_t modR,Int_t ddl,UInt_t *buf){
  //This method packs the SSD digits in a proper 32 bits structure
  // Revised by Enrico Fragiacomo
  Int_t ix;
  Int_t iz;
  Int_t is;
  UInt_t word;
  UInt_t baseWord;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  ofstream ftxt;
  if(ndigits){
    if (fVerbose==2){
      ftxt.open("SSDdigits.txt",ios::app);
    }
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      iz=digs->GetCoord1();  // If iz==0, O side and if iz=1 N side
      ix=digs->GetCoord2();  // Strip Number
      is=digs->GetCompressedSignal();  // ADC Signal
      // cout<<" Module:"<<mod-500<<" N/P side:"<<iz<<" Strip Number:"<<ix<<" Amplidute:"<<is-1<<endl;
      if(is<0) is = 0;
      if(is>4095) is = 4095;
      if (fVerbose==2)
	ftxt<<"DDL:"<<ddl<<" Mod: "<<modR<<" N/P: "<<iz<<" Strip: "<<ix<<" Value: "<<is-1<<endl;

      baseWord=0;

      word=is;
      AliBitPacking::PackWord(word,baseWord,0,11);//ADC data

      word = (iz==0) ? ix : 1535-ix ; // on N-side 1535-768 -> 0-767
      AliBitPacking::PackWord(word,baseWord,12,22);//Strip Number

      word = mod%12; // ADC-number (12 ADCs per AD module)
      word += ( word<6 ) ? 0 : 2; // ADC range 0-5 and 8-13
      AliBitPacking::PackWord(word,baseWord,24,27);//ADC Channel

      word = mod/12+1; // AD-number (AD module index ranges 1-9)
      AliBitPacking::PackWord(word,baseWord,28,31);//AD slot
      fIndex++;
      buf[fIndex]=baseWord;
    }//end for
  }//end if
  if (fVerbose==2)
    ftxt.close();
  return;
}//end GetDigitsSSD

////////////////////////////////////////////////////////////////////////////////////////
//Silicon Drift Detector
//

void AliITSDDLRawData::GetDigitsSDDCompressed(TClonesArray *ITSdigits, Int_t mod, UInt_t *buf){ 
//This method packs the SDD digits in the compressed format (32 bit per digit)
// see AliITSRawStreamSDDCompressed for details on the dta format

  UInt_t dataWord=0;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  if(ndigits){
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      Int_t iz=digs->GetCoord1();  // Anode
      Int_t ix=digs->GetCoord2();  // Time
      Int_t is=digs->GetCompressedSignal();  // ADC Signal - 8 bit
      dataWord=mod<<27;
      Int_t sid=0;
      if(iz>=256){
	sid=1;
	iz-=256;
      }
      dataWord+=sid<<26;
      dataWord+=iz<<18;
      dataWord+=ix<<10;
      UInt_t adcEncoded=0;
      Int_t shift=0;
      if(is < 8) shift=2;
      else if(is<16) shift=3;
      else if(is<32) shift=4;
      else if(is<64) shift=5;
      else if(is<128) shift=6;
      else shift=7;
      adcEncoded=shift+((is-(1<<shift))<<3);
      dataWord+=adcEncoded;
      fIndex++;
      buf[fIndex]=dataWord;
    }
  }
  UInt_t finalWord=UInt_t(15)<<28;
  finalWord+=mod;
  fIndex++;
  buf[fIndex]=finalWord;  
}

//______________________________________________________________________

void AliITSDDLRawData::GetDigitsSDD(TClonesArray *ITSdigits,Int_t mod,Int_t modR,Int_t ddl,UInt_t *buf){  
  //This method packs the SDD digits in a proper 32 bits structure
  Int_t ix;
  Int_t iz;
  Int_t is;
  UInt_t word=0;
  UInt_t baseWord=0;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  ofstream ftxt;
  Int_t digarr[512][256];
  for(Int_t i=0;i<512;i++){
    for(Int_t j=0;j<256;j++){
      digarr[i][j]=0;
    }
  }
  //word to select the 12 carlos for the 12 modules
  UInt_t carlosid=0x30000000+mod;
  
  fIndex++;
  buf[fIndex]=carlosid;
  Int_t first=0;
  Int_t last=0;
  Int_t diff=0;
  Int_t nbit=0;
  UInt_t word2=0;
  Bool_t flag = kFALSE;
  baseWord=0;
  Int_t bitinfo1[4] = {3,8,3,7}; //vector with info on bit for timebin info 
  Int_t wordinfo1[4]= {0,0,0,0}; //vector with word info for timebin info 
  Int_t bitinfo2[2] = {3,18};    //vector with info on bit for EOR (end of row) info
  Int_t wordinfo2[3]= {1,65593}; //vector with word info for anode info

  /* for time bin info: word          n bits   meaning
                         0               3      next info is timebin 
                         8               3      next word is 8 bit long
                       tb value          8      timebin value
		       n (2->7)          3      next info is n bit long
		        signal           n      signal value

     for anode info:     1               3      next 18 bits are for EOR 
                                                increments the anode value

                         EOR             18     error codes + other info
  */
             
  if(ndigits){
    if (fVerbose==2)
      ftxt.open("SDDdigits.txt",ios::app);
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      iz=digs->GetCoord1();  // Anode
      ix=digs->GetCoord2();  // Time
      is=digs->GetCompressedSignal();  // ADC Signal
      digarr[iz][ix]=is;
      if (fVerbose==2)
	ftxt<<"DDL:"<<ddl<<" MID:"<<modR<<" An:"<<iz<<" T:"<<ix<<" A:"<<is<<endl;
      if (is>255){Error("GetDigitsSDD", "bits words is needed)!!!");}
    }
      
    for(Int_t anode=0;anode<512;anode++){
      if(flag){
	last = first+diff-1;
	AliBitPacking::PackWord(word2,baseWord,first,last);
	flag = kFALSE;
	first = last+1;
	diff=0;
      }
      if(anode == 256){
	last = 0;
	first = 0;
	flag = kFALSE;
	diff = 0;
	word2=0;
      }
      
      for(Int_t tb=0;tb<256;tb++){
	if(digarr[anode][tb]!=0){
	  if(flag){      
	    last = first+diff-1;
	    AliBitPacking::PackWord(word2,baseWord,first,last);
	    flag = kFALSE;
	    first = last+1;
	    diff=0;
	  }
	  wordinfo1[1] = tb;
	  //non lossy compression as it is done in Carlos 
	  //(data are already 10to8bit compressed by AMBRA

	  /* if value < 8  value = value - (1 << 2) (word is 2 bit long) 
             if value < 16 value = value - (1 << 3) (word is 3 bit long)
             if value < 32 value = value - (1 << 4) (word is 4 bit long)
	     if value < 64 value = value - (1 << 5) (word is 5 bit long)
	     if value <128 value = value - (1 << 6) (word is 6 bit long)
	     if value >=128value = value - (1 << 7) (word is 7 bit long)

	  */
	  //if(digarr[anode][tb]<4) continue; // not write <4 cnts above tL
	  if(digarr[anode][tb]<8){
	    bitinfo1[3] = 2;
	    wordinfo1[2] = 2;
	    wordinfo1[3] = digarr[anode][tb]-(1 << bitinfo1[3]);
	  }	  	  
	  if(digarr[anode][tb]>=8 && digarr[anode][tb]<16){
	    bitinfo1[3] = 3;
	    wordinfo1[2] = 3;
	    wordinfo1[3] = digarr[anode][tb]-(1 << bitinfo1[3]);
	  }
	  if(digarr[anode][tb]>=16 && digarr[anode][tb]<32){
	    bitinfo1[3] = 4;
	    wordinfo1[2] = 4;
	    wordinfo1[3] = digarr[anode][tb]-(1 << bitinfo1[3]);
	  }
	  if(digarr[anode][tb]>=32 && digarr[anode][tb]<64){
	    bitinfo1[3] = 5;
	    wordinfo1[2] = 5;
	    wordinfo1[3] = digarr[anode][tb]-(1 << bitinfo1[3]);
	  }
	  if(digarr[anode][tb]>=64 && digarr[anode][tb]<128){
	    bitinfo1[3] = 6;
	    wordinfo1[2] = 6;
	    wordinfo1[3] = digarr[anode][tb]-(1 << bitinfo1[3]);
	  }
	  if(digarr[anode][tb]>=128){
	    bitinfo1[3] = 7;
	    wordinfo1[2] = 7;
	    wordinfo1[3] = digarr[anode][tb]-(1 << bitinfo1[3]);
	  }
	  
	  for(Int_t ie=0;ie<4;ie++){
	    
	    if(flag){      
	      last = first+diff-1;
	      AliBitPacking::PackWord(word2,baseWord,first,last);
	      flag = kFALSE;
	      first = last+1;
	      diff=0;
	    }
	    last = first+bitinfo1[ie]-1;
	    if(first < 30 && last < 30){	  	  
	      AliBitPacking::PackWord(wordinfo1[ie],baseWord,first,last); 
	      first = last+1;
	    }
	    else{
	      if(first<=29){
		UInt_t w = AliBitPacking::UnpackWord(wordinfo1[ie],0,29-first);
		AliBitPacking::PackWord(w,baseWord,first,29);
		Int_t lb = 29-first+1;
		diff = bitinfo1[ie]-lb;
		word2 = AliBitPacking::UnpackWord(wordinfo1[ie],lb,lb+diff-1);
		flag = kTRUE;
		if(anode<256) word = 2;//channel 0 of carlos
		else word = 3; //channel 1 of carlos
		AliBitPacking::PackWord(word,baseWord,30,31);
		fIndex++;
		buf[fIndex]=baseWord;
		first=0;
		last = 0;
		baseWord=0;
		word = 0;
	      }
	      else{
		word2 = wordinfo1[ie];
		diff = bitinfo1[ie];
		flag = kTRUE;
		if(anode<256) word = 2; //channel 0 of carlos
		else word = 3; //channel 1 of carlos
		AliBitPacking::PackWord(word,baseWord,30,31);
		fIndex++;
		buf[fIndex]=baseWord;
		first=0;
		last=0;
		baseWord=0;
		word = 0;
	      }
	    }
	  }
	  
	}//END IF
	
      }//end loop on tb
    
      for(Int_t i=0;i<2;i++){
	if(flag){      
	  last = first+diff-1;
	  AliBitPacking::PackWord(word2,baseWord,first,last);
	  flag = kFALSE;
	  first = last+1;
	  diff=0;
	}
	
	word = wordinfo2[i];
	nbit = bitinfo2[i];
	last = first+nbit-1;
	if(first < 30 && last < 30){	  	  
	  AliBitPacking::PackWord(word,baseWord,first,last); //3 bit code =1 -> next 18 bits for EOR
	  first = last+1;
	}
	
	else{
	  if(first<=29){
	    UInt_t w = AliBitPacking::UnpackWord(word,0,29-first);
	    AliBitPacking::PackWord(w,baseWord,first,29);
	    Int_t lb = 29-first+1;
	    diff = nbit-lb;	   
	    word2 = AliBitPacking::UnpackWord(word,lb,lb+diff-1);
	    flag = kTRUE;
	    if(anode<256) word = 2;
	    else word = 3;
	    AliBitPacking::PackWord(word,baseWord,30,31);
	    fIndex++;
	    buf[fIndex]=baseWord;
	    first=0;
	    last = 0;
	    baseWord=0;
	    if(anode==255){
	      flag=kFALSE;
	      word2=0;
	    }
	  }
	  else{
	    word2 = word;
	    diff = nbit;
	    flag = kTRUE;
	    if(anode<256) word = 2;
	    else word = 3;
	    AliBitPacking::PackWord(word,baseWord,30,31);
	    fIndex++;
	    buf[fIndex]=baseWord;
	    first=0;
	    last=0;
	    baseWord=0;
	    if(anode==255){
	      flag=kFALSE;
	      word2=0;
	    }
	  }
	}
      }
    } //end for
    
  }
  if(fVerbose==2)
    ftxt.close();
  return;
}//end GetDigitsSDD

////////////////////////////////////////////////////////////////////////////////////////
//PIXEL 
//

void AliITSDDLRawData::GetDigitsSPD(TClonesArray *ITSdigits,Int_t mod,Int_t ddl, UInt_t *buf, AliITSFOSignalsSPD* foSignals){
  //This method packs the SPD digits in a proper 32 structure
  //Since data is zero suppressed,the coordinates for the chip having zero digits 
  //doesn't get listed in the galice.root file. However the SPD format requires 
  //the empty chip to be written with chip header and chip trailer.

  Int_t chipLow  = AliITSRawStreamSPD::GetOnlineChipFromOffline(mod,0);
  Int_t chipHigh = AliITSRawStreamSPD::GetOnlineChipFromOffline(mod,159);

  if (chipLow>chipHigh) {chipLow  -= 4; chipHigh += 4;}
  UInt_t eq = AliITSRawStreamSPD::GetOnlineEqIdFromOffline(mod);
  UInt_t hs = AliITSRawStreamSPD::GetOnlineHSFromOffline(mod);

  // create int map to later hold all digits sorted
  AliITSIntMap* digMap = new AliITSIntMap();

  UInt_t baseWord=0;
 
  Int_t ndigits = ITSdigits->GetEntries(); //number of digits in the current module
  //cout<<"      Number of digits in the current module:"<<ndigits<<" module:"<<mod<<endl;

  // _______________________________________________________________________
  // Preprocess the digits - sort them in integer map (Henrik Tydesjo)
  // Needed to have exact same order as in real raw data
  AliITSdigit *digs;
 ofstream ftxt;
  if (ndigits) {
    //loop over digits
    if (fVerbose==2) ftxt.open("SPDdigits.txt",ios::app);
    for (Int_t digit=0; digit<ndigits; digit++) {
      digs = (AliITSdigit*) ITSdigits->UncheckedAt(digit);
      /*---------------------------------------------------------------------------
       *     Each module contains 5 read out chips of 256 rows and 32 columns.
       *     So, the cell number in Z direction varies from 0 to 159.
       *     ---------------------------------------------------------------------*/
      Int_t iz=digs->GetCoord1();  // Cell number in Z direction 
      Int_t ix=digs->GetCoord2();  // Cell number in X direction

      if(fVerbose==2) ftxt<<"DDL:"<<ddl<<" Mod:"<<mod<<" Row:"<<ix<<" Col:"<<iz<<endl;
      UInt_t dummyDDL, dummyHS, chip, col, row;
      AliITSRawStreamSPD::OfflineToOnline(mod,iz,ix,dummyDDL,dummyHS,chip,col,row);

      //  insert digit into map...
      // (reverse order of cols and rows as in real raw data)
      digMap->Insert(chip*256*32+(31-col)*256+(255-row),row);
    }
  }

   // _______________________________________________________________________
  // Procedure for writing raw data (Henrik Tydesjo)
  // Reimplemented because of unreadability (5 Mar 2009)
  // Now also with fast-or signals
  Int_t  previousChip = chipLow-1;
  Int_t  chip         = chipLow-1;
  UInt_t chipHitCount = 0;


  UInt_t nrHits = digMap->GetNrEntries();
  for (UInt_t nHit=0; nHit<nrHits; nHit++) {

    Int_t key  = digMap->GetKeyIndex(nHit);
    chip = key/(256*32);
    Int_t col  = 31 - (key%(256*32))/256;
    Int_t row  = digMap->GetValIndex(nHit);

    // add trailer for previous chip (if there was one...)
    if (chip>previousChip && previousChip>chipLow-1) {
      WriteChipTrailer(buf, chipHitCount, foSignals->GetSignal(eq,hs,previousChip), baseWord);
    }

    // add headers/trailers for chips without hits (if any)
    for (Int_t ch=previousChip+1; ch<chip; ch++) {
      WriteChipHeader(ch, hs, baseWord);
      WriteChipTrailer(buf, 0, foSignals->GetSignal(eq,hs,ch), baseWord);
    }

    // if new chip, add header
    if (chip>previousChip) {
      WriteChipHeader(chip, hs, baseWord);
      chipHitCount = 0;
      previousChip = chip;
    }

    chipHitCount++;

    // add pixel hit
    WriteHit(buf,row,col,baseWord);

  }

  // add trailer for last chip (if there was one...)
  if (chip>chipLow-1) {
    WriteChipTrailer(buf, chipHitCount, foSignals->GetSignal(eq,hs,chip), baseWord);
  }

  // add REMAINING headers/trailers for chips without hits (if any)
  for (Int_t ch=chip+1; ch<=chipHigh; ch++) {
    WriteChipHeader(ch, hs, baseWord);
    WriteChipTrailer(buf, 0, foSignals->GetSignal(eq,hs,ch), baseWord);
  }
  // _______________________________________________________________________


  delete digMap;

  if(fVerbose==2)
    ftxt.close();
  return;
}//end GetDigitsSPD

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSPD(TBranch* branch, AliITSFOSignalsSPD* foSignals){
  //This method creates the Raw data files for SPD detectors
  const Int_t kSize=21000; //256*32*5=40960 max number of digits per module
  UInt_t buf[kSize];      //One buffer cell can contain 2 digits 
  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  TString fileName;
  AliFstream* outfile;         // logical name of the output file 
  AliRawDataHeaderSim header;

  //loop over DDLs
  for(Int_t ddl=0;ddl<AliDAQ::NumberOfDdls("ITSSPD");ddl++){
    fileName.Form("%s",AliDAQ::DdlFileName("ITSSPD",ddl)); //The name of the output file.
    outfile = new AliFstream(fileName.Data());
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile->Tellp();
    outfile->WriteBuffer((char*)(&header),sizeof(header));
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<AliITSRawStreamSPD::kModulesPerDDL; mod++){
      Int_t moduleNumber = AliITSRawStreamSPD::GetModuleNumber(ddl, mod);
      digits->Clear();
      branch->GetEvent(moduleNumber);
      //For each Module, buf contains the array of data words in Binary format	  
      //fIndex gives the number of 32 bits words in the buffer for each module
      GetDigitsSPD(digits, moduleNumber, ddl, buf, foSignals);
      outfile->WriteBuffer((char *)buf,((fIndex+1)*sizeof(UInt_t)));
      for(Int_t i=0;i<(fIndex+1);i++){
	buf[i]=0;
      }//end for
      fIndex=-1;
    }//end for
    
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile->Tellp();
    outfile->Seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    outfile->WriteBuffer((char*)(&header),sizeof(header));
    delete outfile;
  }//end for

  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSSD(TBranch* branch){

    //This method creates the Raw data files for SSD detectors
  const Int_t kSize=1536;//768*2 Number of stripe * number of sides(N and P)
  UInt_t buf[kSize];      
  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  TString fileName;
  AliFstream* outfile;         // logical name of the output file 
  AliRawDataHeaderSim header;

  //loop over DDLs  
  for(Int_t i=0;i<AliDAQ::NumberOfDdls("ITSSSD");i++){
    fileName.Form("%s",AliDAQ::DdlFileName("ITSSSD",i)); //The name of the output file.
    outfile = new AliFstream(fileName.Data());
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile->Tellp();
    outfile->WriteBuffer((char*)(&header),sizeof(header));
    
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<AliITSRawStreamSSD::kModulesPerDDL; mod++){
      Int_t moduleNumber = AliITSRawStreamSSD::GetModuleNumber(i, mod);
      if(moduleNumber!=-1){
	digits->Clear();
	branch->GetEvent(moduleNumber);
	//For each Module, buf contains the array of data words in Binary format	  
	//fIndex gives the number of 32 bits words in the buffer for each module
	GetDigitsSSD(digits,mod,moduleNumber,i,buf);
	outfile->WriteBuffer((char *)buf,((fIndex+1)*sizeof(UInt_t)));
	fIndex=-1;
      }//end if
    }//end for

    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile->Tellp();
    outfile->Seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile->WriteBuffer((char*)(&header),sizeof(header));
    delete outfile;
  }//end for

  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSDD(TBranch* branch, AliITSDDLModuleMapSDD* ddlsdd){
    //This method creates the Raw data files for SDD detectors
  const Int_t kSize=131072; //256*512
  UInt_t buf[kSize];      
  fIndex=-1;

  TClonesArray*& digits = * (TClonesArray**) branch->GetAddress();
  TString fileName;
  AliFstream* outfile;             // logical name of the output file 
  AliRawDataHeaderSim header;
  
  if(fSDDRawFormat!=0){ 
    for(Int_t ibit=0; ibit<8; ibit++) header.SetAttribute(ibit);
  }else{
    for(Int_t ibit=0; ibit<5; ibit++) header.SetAttribute(ibit);
    for(Int_t ibit=5; ibit<8; ibit++) header.ResetAttribute(ibit);  
  }
  UInt_t skippedword=0; 
  UInt_t carlosFooterWord=0;
  UInt_t fifoFooterWord=0;
  UInt_t jitterWord=0;
  Bool_t retcode;
  retcode = AliBitPacking::PackWord(0x3FFFFFFF,carlosFooterWord,0,31);
  retcode = AliBitPacking::PackWord(0x3F1F1F1F,fifoFooterWord,0,31);
  if(fSDDRawFormat!=0) retcode = AliBitPacking::PackWord(0x7F000000,jitterWord,0,31);
  else retcode = AliBitPacking::PackWord(0x80000000,jitterWord,0,31);
 
  //loop over DDLs  
  for(Int_t i=0;i<AliDAQ::NumberOfDdls("ITSSDD");i++){
    fileName.Form("%s",AliDAQ::DdlFileName("ITSSDD",i)); //The name of the output file.
    outfile = new AliFstream(fileName.Data());
    //write Dummy DATA HEADER
    UInt_t dataHeaderPosition=outfile->Tellp();
    outfile->WriteBuffer((char*)(&header),sizeof(header));


    //first 1 "dummy" word to be skipped
    if(fSDDRawFormat!=0){
      retcode = AliBitPacking::PackWord(0xFFFFFFFF,skippedword,0,31);
      outfile->WriteBuffer((char*)(&skippedword),sizeof(skippedword));
    }

    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<AliITSRawStreamSDD::kModulesPerDDL; mod++){
      Int_t moduleNumber = ddlsdd->GetModuleNumber(i, mod);
      if(moduleNumber!=-1){
	digits->Clear();
	branch->GetEvent(moduleNumber);

	//For each Module, buf contains the array of data words in Binary format	  
	//fIndex gives the number of 32 bits words in the buffer for each module
	//	cout<<"MODULE NUMBER:"<<mapSDD[i][mod]<<endl;
	if(fSDDRawFormat==0){
	  GetDigitsSDDCompressed(digits,mod,buf);
	  outfile->WriteBuffer((char *)buf,((fIndex+1)*sizeof(UInt_t)));
	}else{
	  GetDigitsSDD(digits,mod,moduleNumber,i,buf);
	  outfile->WriteBuffer((char *)buf,((fIndex+1)*sizeof(UInt_t)));
	  for(Int_t iw=0;iw<3;iw++) outfile->WriteBuffer((char*)(&carlosFooterWord),sizeof(carlosFooterWord));
	}
	fIndex=-1;
      }//end if
    }//end for
    // 12 words with FIFO footers (=4 FIFO x 3 3F1F1F1F words per DDL)
    if(fSDDRawFormat!=0){
      for(Int_t iw=0;iw<12;iw++) outfile->WriteBuffer((char*)(&fifoFooterWord),sizeof(fifoFooterWord));
    }
    outfile->WriteBuffer((char*)(&jitterWord),sizeof(jitterWord));      
    //Write REAL DATA HEADER
    UInt_t currentFilePosition=outfile->Tellp();
    outfile->Seekp(dataHeaderPosition);
    header.fSize=currentFilePosition-dataHeaderPosition;
    header.SetAttribute(0);  // valid data
    outfile->WriteBuffer((char*)(&header),sizeof(header));
    delete outfile;
  }//end for

  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::WriteChipHeader(Int_t ChipAddr,Int_t halfStave,UInt_t &BaseWord){
  //This method writes a chip header 
  //cout<<"Chip: "<<ChipAddr<<" Half Stave module:"<<halfStave<<endl;
  BaseWord=0;
  AliBitPacking::PackWord(ChipAddr,BaseWord,16,19);
  //  At the moment the event count is always 0 (bits 20-26)
  AliBitPacking::PackWord(0,BaseWord,20,26);
  AliBitPacking::PackWord(halfStave,BaseWord,27,29);
  AliBitPacking::PackWord(0x1,BaseWord,30,31);
  return;
}//end WriteChipHeader

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::WriteChipTrailer(UInt_t *buf, Int_t ChipHitCount, Bool_t foBit, UInt_t &BaseWord){
  //This method writes a chip trailer
  //pixel fill word
  if((ChipHitCount%2)!=0){
    AliBitPacking::PackWord(0xC000,BaseWord,16,31);
  }
  AliBitPacking::PackWord(ChipHitCount,BaseWord,0,11);
  AliBitPacking::PackWord(0x0,BaseWord,12,12);
  AliBitPacking::PackWord(foBit,BaseWord,13,13);
  AliBitPacking::PackWord(0x0,BaseWord,14,15);
  fIndex++;
  buf[fIndex]=BaseWord;
  BaseWord=0;
  return;
}//end WriteChipTrailer

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::WriteHit(UInt_t *buf,Int_t RowAddr,Int_t HitAddr,UInt_t &BaseWord){
  //This method writs an hit
  if(!BaseWord){
    AliBitPacking::PackWord(HitAddr,BaseWord,16,20);
    AliBitPacking::PackWord(RowAddr,BaseWord,21,28);
    AliBitPacking::PackWord(2,BaseWord,30,31);
  }//end if
  else{
    AliBitPacking::PackWord(HitAddr,BaseWord,0,4);
    AliBitPacking::PackWord(RowAddr,BaseWord,5,12);
    AliBitPacking::PackWord(2,BaseWord,14,15);
    fIndex++;
    buf[fIndex]=BaseWord;
    BaseWord=0;
  }//end else
  return;
}//end WriteHit

