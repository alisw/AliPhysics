#ifndef ALITPCMONITORDATEFILE_H
#define ALITPCMONITORDATEFILE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDateFile class
//// 
//// Class for handling the data structure in a DATE file
//// 
//// Authors: Roland Bramm
////          Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////


#include <string>
#include "TNamed.h"

using namespace std;

class AliTPCMonitorDateFile : public TNamed {
	public: 
	AliTPCMonitorDateFile();
	AliTPCMonitorDateFile(const  AliTPCMonitorDateFile &datefile);
	AliTPCMonitorDateFile& operator= (const AliTPCMonitorDateFile& datefile);
	~AliTPCMonitorDateFile();
	
	void    AllocateArray(int size);
	void    CloseDateFile();
	
	Int_t   GetAllocatedSizeofArray() const;
	Int_t   GetFileSize() const; 
	Int_t   GetFilePosition() const;
	Char_t* GetMemoryPointer();
	
	Bool_t  IsLastEvent() const;
	Bool_t  IsEventValid() const;
	
	Bool_t  IsDateFileOpen();
	void    OpenDateFile(string name);
	
	void    ReadEvent();
	void    ResetFilePos();




	private:
	
	void SetFileSize();

	Int_t          ffilePos;                       // position in file
	Char_t*        fbigMem;                        // array for event data
        UInt_t         fbigMemsize;                    // size of data array
	Bool_t         fisBigMemAllocated;             // flag for already allocated array
	Int_t          ffileSize;                      // size of DATE file
	string         ffilename;                      // name of DATE file
	Bool_t         finitFile;                      // flag for opened file
	Bool_t         freadPosOverflow;               // data position overflow flag 
	ifstream *     fin;                            // file to be read
	
	ClassDef(AliTPCMonitorDateFile,1);
};

#endif

