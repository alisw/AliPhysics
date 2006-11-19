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

// Storing digits in a binary file
// according to the DDL mapping
// To be used in Alice Data Challenges
// This class is used by AliVZERODDL.C macro
// Author: B. Cheynis

#include <Riostream.h>
#include <TObjArray.h>
#include "AliRawDataHeader.h"
#include "AliVZEROBuffer.h"

//#include "TFile.h"
//#include "TTree.h"

ClassImp(AliVZEROBuffer)

//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer():TObject(),
    fVerbose(0),
    f(),
    fNumberOfDigits(0)
{
  //
  // default constructor
  //
}
//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer(const char* fileName):TObject(),
    fVerbose(0),
    f(),
    fNumberOfDigits(0)
{
  // Constructor
#ifndef __DECCXX
  f.open(fileName,ios::binary|ios::out);
#else
  f.open(fileName,ios::out);
#endif
  // fout=new TFile(fileName,"recreate");
  // tree=new TTree("tree","Values");
  AliRawDataHeader header;
  f.write((char*)(&header), sizeof(header));

}

//_____________________________________________________________________________
AliVZEROBuffer::~AliVZEROBuffer(){
  // Destructor, it closes the IO stream
  AliRawDataHeader header;
  header.fSize = f.tellp();
  header.SetAttribute(0);  // valid data
  f.seekp(0);
  f.write((char*)(&header), sizeof(header));
  f.close();
  //delete tree;
  //delete fout;
}

//_____________________________________________________________________________
AliVZEROBuffer::AliVZEROBuffer(const AliVZEROBuffer &source):TObject(source),
   fVerbose(0),
   f(),
   fNumberOfDigits(0)

{
  // Copy Constructor
  this->fVerbose=source.fVerbose;
  return;
}

//_____________________________________________________________________________
AliVZEROBuffer& AliVZEROBuffer::operator=(const AliVZEROBuffer &source)

{
  //Assigment operator
  this->fVerbose=source.fVerbose;
  return *this;
}

//_____________________________________________________________________________
void AliVZEROBuffer::WriteBinary(Int_t cell,Int_t ADC, Int_t Time){
  // It writes VZERO digits as a raw data file. 
  // Being called by AliVZERODDL.C

  struct DataFile{
    Int_t cell;
    Int_t ADC;
    Int_t Time;
  };
  
  DataFile  data;
  data.cell = cell;
  data.ADC  = ADC;
  data.Time = Time;

  fNumberOfDigits++;
  f.write((char*)(&data),sizeof(data));

}

