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
Revision 1.2  2007/10/12 13:36:27  cvetan
Coding convention fixes from Stefan

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDateFile class
//// 
//// Class for handling the data structure in a DATE file
//// Used to read DATE files for the TPC raw data Monitor 
////
//// Author: Roland Bramm
////         Stefan Kniege, IKF, Frankfurt
////       
////
/////////////////////////////////////////////////////////////////////////



#include "AliTPCMonitorDateFile.h"
#include "AliTPCMonitorDateFormat.h"
#include <Riostream.h>
ClassImp(AliTPCMonitorDateFile)

//____________________________________________________________________________
AliTPCMonitorDateFile::AliTPCMonitorDateFile() :
  ffilePos(0),
  fbigMem(0),
  fbigMemsize(0),
  fisBigMemAllocated(false),
  ffileSize(0),
  ffilename(""),
  finitFile(false),
  freadPosOverflow(false),
  fin(new ifstream())
{
  // Constructor
}

//____________________________________________________________________________
AliTPCMonitorDateFile::AliTPCMonitorDateFile(const AliTPCMonitorDateFile &datefile) :
  TNamed(datefile.GetName(),datefile.GetTitle()),
  ffilePos(datefile.ffilePos),
  fbigMem(datefile.fbigMem),
  fbigMemsize(datefile.fbigMemsize),
  fisBigMemAllocated(datefile.fisBigMemAllocated),
  ffileSize(datefile.ffileSize),
  ffilename(datefile.ffilename),
  finitFile(datefile.finitFile),
  freadPosOverflow(datefile.freadPosOverflow),
  fin(new ifstream())
{
  // copy constructor 
}

//____________________________________________________________________________
AliTPCMonitorDateFile &AliTPCMonitorDateFile:: operator= (const AliTPCMonitorDateFile& datefile)
{

  // assignment operator 
  if(this!=&datefile)
    {
      ffilePos=datefile.ffilePos;
      fbigMem=datefile.fbigMem;
      fbigMemsize=datefile.fbigMemsize;
      fisBigMemAllocated=datefile.fisBigMemAllocated;
      ffileSize=datefile.ffileSize;
      ffilename=datefile.ffilename;
      finitFile=datefile.finitFile;
      freadPosOverflow=datefile.freadPosOverflow;
      fin = new ifstream();
    }
  return *this;
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
Int_t AliTPCMonitorDateFile::GetFileSize() const
{
  // Return size of DATE file
  return ffileSize;
}

//____________________________________________________________________________
void AliTPCMonitorDateFile::ReadEvent() 
{
  // Read in event from file
  Int_t size;
  Char_t         fmem[512];                      // array for event header
  Char_t swapcarry[4];
	Bool_t toSwapEndian = false;
	//Fetch some bytes to get headers to know how much
	fin->seekg(ffilePos);
	//if( (getFilePosition() + sizeof(fmem)) <= (unsigned int)getFileSize()){
	UInt_t sfmem = 272;
	if( (GetFilePosition() + sfmem) <= (UInt_t)GetFileSize()){

	  
	  fin->read((Char_t*)fmem,sfmem);
		fin->seekg(ffilePos);
		AliTPCMonitorDateFormat *dateform;
		dateform = new AliTPCMonitorDateFormat((Char_t *)&fmem);
		

	  
		toSwapEndian = dateform->IsEventWrongEndian();
		if(toSwapEndian == true){
			delete dateform;
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
			dateform = new AliTPCMonitorDateFormat((Char_t *)&fmem);
		}
		size = dateform->GetEventSize();
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
		delete dateform;
	}else{
	  freadPosOverflow = true;
	}
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFile::IsLastEvent() const
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
Bool_t AliTPCMonitorDateFile::IsEventValid() const
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
Int_t AliTPCMonitorDateFile::GetFilePosition() const
{
  // Return current position in file
  return ffilePos;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFile::GetAllocatedSizeofArray() const
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

