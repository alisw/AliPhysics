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

/*
$Log$
*/ 

#include "AliTPCMonitorDateFile.h"

ClassImp(AliTPCMonitorDateFile)

//____________________________________________________________________________
AliTPCMonitorDateFile::AliTPCMonitorDateFile() 
{
  // Constructor
  fin = new ifstream();
  fisBigMemAllocated = false;
  freadPosOverflow = false;
  ffilePos = 0;
  ffileSize = 0;
  fbigMemsize = 0;
  finitFile = false;
  ffilename = "";
}

//____________________________________________________________________________
AliTPCMonitorDateFile::~AliTPCMonitorDateFile() 
{
  // Destructor
  delete fin;
  if(fisBigMemAllocated == true){
    delete[] fbigMem;
  }
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::OpenDateFile(string name) 
{
  // Open DATE file 
  if( (ffilename != name) && (finitFile == true) ) {
    CloseDateFile();
    ResetFilePos();
  }
  fin->open(name.c_str());
  ffilename = name;
  SetFileSize();
  finitFile = true;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFile::IsDateFileOpen() 
{
  // Check if DATE file is open
  return fin->is_open();
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::CloseDateFile() 
{
  // Close DATE file
  fin->close();
  freadPosOverflow = false;
  delete fin;
  fin = new ifstream();
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::ResetFilePos() 
{
  // Reset file position 
  ffilePos = 0;
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::SetFileSize() 
{
  // Set size of DATE file
  fin->seekg (0, ios::end);
  ffileSize = fin->tellg();
  fin->seekg (0, ios::beg);
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFile::GetFileSize() 
{
  // Return size of DATE file
  return ffileSize;
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::ReadEvent() 
{
  // Read in event from file
  Int_t size;
  Char_t swapcarry[4];
	Bool_t toSwapEndian = false;
	//Fetch some bytes to get headers to know how much
	fin->seekg(ffilePos);
	//if( (getFilePosition() + sizeof(fmem)) <= (unsigned int)getFileSize()){
	UInt_t sfmem = 272;
	if( (GetFilePosition() + sfmem) <= (UInt_t)GetFileSize()){

	  
	  fin->read((Char_t*)fmem,sfmem);
		fin->seekg(ffilePos);
		AliTPCMonitorDateFormat *DateForm;
		DateForm = new AliTPCMonitorDateFormat((Char_t *)&fmem);
		

	  
		toSwapEndian = DateForm->IsEventWrongEndian();
		if(toSwapEndian == true){
			delete DateForm;
			for(Int_t i = 0; i < 68; i++) {
				swapcarry[0] = fmem[(i*4)+0];
				swapcarry[1] = fmem[(i*4)+1];
				swapcarry[2] = fmem[(i*4)+2];
				swapcarry[3] = fmem[(i*4)+3];
				fmem[(i*4)+0] = swapcarry[3];
				fmem[(i*4)+1] = swapcarry[2];
				fmem[(i*4)+2] = swapcarry[1];
				fmem[(i*4)+3] = swapcarry[0];
			}
			DateForm = new AliTPCMonitorDateFormat((Char_t *)&fmem);
		}
		size = DateForm->GetEventSize();
		if(size > GetAllocatedSizeofArray()) {
		  AllocateArray((Int_t)(size*1.1));
		}
		if( (GetFilePosition() + size) <= GetFileSize()){
		  fin->read((Char_t*)fbigMem,size);
		  if(toSwapEndian == true){
				for(Int_t i = 0; i < (size/4); i++) {
					swapcarry[0] = fbigMem[(i*4)+0];
					swapcarry[1] = fbigMem[(i*4)+1];
					swapcarry[2] = fbigMem[(i*4)+2];
					swapcarry[3] = fbigMem[(i*4)+3];
					fbigMem[(i*4)+0] = swapcarry[3];
					fbigMem[(i*4)+1] = swapcarry[2];
					fbigMem[(i*4)+2] = swapcarry[1];
					fbigMem[(i*4)+3] = swapcarry[0];
				}
			}
			ffilePos += size;
		}else{
		  freadPosOverflow = true;
		}
		delete DateForm;
	}else{
	  freadPosOverflow = true;
	}
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFile::IsLastEvent() 
{
  // Check if event is last event in file
  Bool_t retval;
  if( (ffilePos < ffileSize) && (freadPosOverflow == false) )
    retval = false;
  else
    retval = true;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFile::IsEventValid()
{
  // Check Over flow flag 
  Bool_t retval;
  if(freadPosOverflow == false)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFile::GetFilePosition() 
{
  // Return current position in file
  return ffilePos;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFile::GetAllocatedSizeofArray() 
{
  // Return size of allocated data array
  return fbigMemsize;
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::AllocateArray(Int_t size) 
{
  // allocate data array
  if(fisBigMemAllocated == false) {
    fbigMem = new Char_t[size];
    fbigMemsize = size;
    fisBigMemAllocated = true;
  } else {
    delete fbigMem;
    fbigMem = new Char_t[size];
    fbigMemsize = size;
    fisBigMemAllocated = true;
  }
}

//____________________________________________________________________________
Char_t *AliTPCMonitorDateFile::GetMemoryPointer() 
{
  // Return pointer to data array
  return fbigMem;
}

