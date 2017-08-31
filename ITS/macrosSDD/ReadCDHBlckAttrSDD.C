#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH2F.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSDDCompressed.h"
#endif

// Macro to easily print the content of the CDH Block Attributes
// Author: M. Sitta,  sitta@to.infn.it

void PrintOutRecordInfo(Char_t *recordType,
			Int_t evNum,
			AliRawReader *rd,
			UChar_t firstCDH,
			Bool_t printDetails=kTRUE){

  printf("\n%s record at iev %d: GDC %d\n", recordType, evNum, rd->GetGDCId());
  if (printDetails) {
    printf(" -> DDL %2d on LDC %2d -- ", rd->GetDDLID(), rd->GetLDCId());
    printf("CDH Attribute is 0x%02x\n", firstCDH);

    UChar_t *data;
    UChar_t cdhAttr;

    while(1) {
      Bool_t nodata = kFALSE;
      do{
	if(!rd->ReadNextData(data)) {
	  nodata = kTRUE;
	  break;
	}
      }while(rd->GetDataSize()==0);
      if (nodata) break;
      cdhAttr = rd->GetBlockAttributes();

      printf(" -> DDL %2d on LDC %2d -- ", rd->GetDDLID(), rd->GetLDCId());
      printf("CDH Attribute is 0x%02x\n", cdhAttr);
    }
  }

}

void ReadCDHBlckAttrSDD(Char_t datafil[100],
			Int_t  firstEv=0,
			Int_t  lastEv=20,
			Bool_t printPhysics=kFALSE){

  // Main function

  Int_t iev=firstEv;
  AliRawReader *rd; 
  if(strstr(datafil,".root")!=0){
    rd=new AliRawReaderRoot(datafil,iev);
  }else{
    rd=new AliRawReaderDate(datafil,iev);
  }
  TStopwatch *evtime=new TStopwatch();

  UChar_t *data;
  UChar_t cdhAttr;
  do{
    evtime->Start();

    rd->Reset();
    rd->Select("ITSSDD");
    Bool_t nodata = kFALSE;
    do{
      if(!rd->ReadNextData(data)) {
	nodata = kTRUE;
	break;
      }
    }while(rd->GetDataSize()==0);
    if (nodata) continue;
    cdhAttr = rd->GetBlockAttributes();

    if (rd->GetType() == AliRawEventHeaderBase::kStartOfData)
      PrintOutRecordInfo("SOD",iev,rd,cdhAttr);

    if (rd->GetType() == AliRawEventHeaderBase::kEndOfData)
      PrintOutRecordInfo("EOD",iev,rd,cdhAttr);

    if (rd->GetType() == AliRawEventHeaderBase::kPhysicsEvent) {
      if (printPhysics)
	PrintOutRecordInfo("PHYSICS",iev,rd,cdhAttr);
      else
	PrintOutRecordInfo("PHYSICS",iev,rd,cdhAttr,kFALSE);
    }

    iev++;
  }while(rd->NextEvent()&&iev<=lastEv);

}

void ReadCDHBlckAttrSDD(Int_t  nrun,
			Int_t  n1,
			Int_t  n2,
			Int_t  year=2017,
			Char_t* dir="LHC17a",
			Int_t  firstEv=0,
			Int_t  lastEv=20,
			Bool_t printPhysics=kFALSE){

  // Get file directly from alien

  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/%d/%s/%09d/raw/%02d%09d%03d.%d.root",year,dir,nrun,year-2000,nrun,n1,n2);
  printf("Open file %s\n",filnam);
  ReadCDHBlckAttrSDD(filnam,firstEv,lastEv,printPhysics);
}
