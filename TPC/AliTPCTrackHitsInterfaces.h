#include "AliClassInfo.h"
#include "AliTPCTrackHits.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id */

/************************************************/
/* Automaticaly generated interface for class     
                 AliTrackHitsInfo                                
**************************************************/


class AliClassAliTrackHitsInfo : public AliClassInfo {
public:
	AliClassAliTrackHitsInfo(){
	  SetName("AliTrackHitsInfo");
	  SetTitle("Interface for AliTrackHitsInfo class ");
	  fgList.Add(this);
	  fSize = sizeof(AliTrackHitsInfo);
	}
	const char * GetClassName(){ return "AliTrackHitsInfo";}
	virtual TClass* GetClass(){return AliTrackHitsInfo::Class();}
	void CTORBuffer(void * pointer, UInt_t size=1)
	{
	  AliTrackHitsInfo * last = &((AliTrackHitsInfo*)pointer)[size];
	  AliTrackHitsInfo * p = (AliTrackHitsInfo*)pointer;
	  while (p!=last) new (p++)AliTrackHitsInfo;
	}
	void DTORBuffer(void * pointer, UInt_t size=1)
	{
	  AliTrackHitsInfo * last = &((AliTrackHitsInfo*)pointer)[size];
	  AliTrackHitsInfo * p = (AliTrackHitsInfo*)pointer;
	  while (p!=last) (p++)->~AliTrackHitsInfo();
	}
	void StreamBuffer(TBuffer &b,const void * pointer, UInt_t size=1)
	{
	  for (UInt_t i=0;i<size;i++) ((AliTrackHitsInfo*)pointer)[i].Streamer(b);
	}
	  void ObjectDump(void *p) {((AliTrackHitsInfo*)p)->Dump();}
};
AliClassAliTrackHitsInfo galiclass____AliTrackHitsInfo; 

/************************************************/
/* Automaticaly generated interface for class     
                 AliTrackHitsParam                                
**************************************************/


class AliClassAliTrackHitsParam : public AliClassInfo {
public:
	AliClassAliTrackHitsParam(){
	  SetName("AliTrackHitsParam");
	  SetTitle("Interface for AliTrackHitsParam class ");
	  fgList.Add(this);
	  fSize = sizeof(AliTrackHitsParam);
	}
	const char * GetClassName(){ return "AliTrackHitsParam";}
	virtual TClass* GetClass(){return AliTrackHitsParam::Class();}
	void CTORBuffer(void * pointer, UInt_t size=1)
	{
	  AliTrackHitsParam * last = &((AliTrackHitsParam*)pointer)[size];
	  AliTrackHitsParam * p = (AliTrackHitsParam*)pointer;
	  while (p!=last) new (p++)AliTrackHitsParam;
	}
	void DTORBuffer(void * pointer, UInt_t size=1)
	{
	  AliTrackHitsParam * last = &((AliTrackHitsParam*)pointer)[size];
	  AliTrackHitsParam * p = (AliTrackHitsParam*)pointer;
	  while (p!=last) (p++)->~AliTrackHitsParam();
	}
	void StreamBuffer(TBuffer &b,const void * pointer, UInt_t size=1)
	{
	  for (UInt_t i=0;i<size;i++) ((AliTrackHitsParam*)pointer)[i].Streamer(b);
	}
	  void ObjectDump(void *p) {((AliTrackHitsParam*)p)->Dump();}
};
AliClassAliTrackHitsParam galiclass____AliTrackHitsParam; 

/************************************************/
/* Automaticaly generated interface for class     
                 AliHitInfo                                
**************************************************/


class AliClassAliHitInfo : public AliClassInfo {
public:
	AliClassAliHitInfo(){
	  SetName("AliHitInfo");
	  SetTitle("Interface for AliHitInfo class ");
	  fgList.Add(this);
	  fSize = sizeof(AliHitInfo);
	}
	const char * GetClassName(){ return "AliHitInfo";}
	virtual TClass* GetClass(){return AliHitInfo::Class();}
	void CTORBuffer(void * pointer, UInt_t size=1)
	{
	  AliHitInfo * last = &((AliHitInfo*)pointer)[size];
	  AliHitInfo * p = (AliHitInfo*)pointer;
	  while (p!=last) new (p++)AliHitInfo;
	}
	void DTORBuffer(void * pointer, UInt_t size=1)
	{
	  AliHitInfo * last = &((AliHitInfo*)pointer)[size];
	  AliHitInfo * p = (AliHitInfo*)pointer;
	  while (p!=last) (p++)->~AliHitInfo();
	}
	void StreamBuffer(TBuffer &b,const void * pointer, UInt_t size=1)
	{
	  for (UInt_t i=0;i<size;i++) ((AliHitInfo*)pointer)[i].Streamer(b);
	}
	  void ObjectDump(void *p) {((AliHitInfo*)p)->Dump();}
};
AliClassAliHitInfo galiclass____AliHitInfo; 
