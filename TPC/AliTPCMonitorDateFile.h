#ifndef ALITPCMONITORDATEFILE_H
#define ALITPCMONITORDATEFILE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <string>
#include <fstream>
#include "AliTPCMonitorDateFormat.h"
#include "TNamed.h"

using namespace std;

class AliTPCMonitorDateFile : public TNamed {
	public:
	AliTPCMonitorDateFile();
	~AliTPCMonitorDateFile();
	
	void    AllocateArray(int size);
	void    CloseDateFile();
	
	Int_t   GetAllocatedSizeofArray();
	Int_t   GetFileSize(); 
	Int_t   GetFilePosition();
	Char_t* GetMemoryPointer();
	
	Bool_t  IsLastEvent();
	Bool_t  IsEventValid();
	
	Bool_t  IsDateFileOpen();
	void    OpenDateFile(string name);
	
	void    ReadEvent();
	void    ResetFilePos();




	private:
	void SetFileSize();
	ifstream *     fin;                            // file to be read
	Int_t          ffilePos;                       // position in file
	Char_t         fmem[512];                      // array for event header
	Char_t*        fbigMem;                        // array for event data
        UInt_t         fbigMemsize;                    // size of data array
	Bool_t         fisBigMemAllocated;             // flag for already allocated array
	Int_t          ffileSize;                      // size of DATE file
	string         ffilename;                      // name of DATE file
	Bool_t         finitFile;                      // flag for opened file
	Bool_t         freadPosOverflow;               // data position overflow flag 
	
	ClassDef(AliTPCMonitorDateFile,1);
};

#endif

