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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD raw data conversion class                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometryFull.h"
#include "AliTRDparameter.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDarrayI.h"

ClassImp(AliTRDrawData)

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData():TObject()
{

  fDebug         = 0;
  fDigitsManager = NULL;

}

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData(const AliTRDrawData &r)
{
  //
  // AliTRDrawData copy constructor
  //

  ((AliTRDrawData &) r).Copy(*this);

}

//_____________________________________________________________________________
AliTRDrawData::~AliTRDrawData()
{
  //
  // Destructor
  //

  if (fDigitsManager) {
    delete fDigitsManager;
    fDigitsManager = NULL;
  }

}

//_____________________________________________________________________________
AliTRDrawData &AliTRDrawData::operator=(const AliTRDrawData &r)
{
  //
  // Assignment operator
  //

  if (this != &r) ((AliTRDrawData &) r).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDrawData::Copy(TObject &r)
{
  //
  // Copy function
  //

  ((AliTRDrawData &) r).fDebug         = fDebug;
  ((AliTRDrawData &) r).fDigitsManager = NULL;

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::OpenInput(const Char_t *name)
{
  //
  // Opens a ROOT-file with the TRD digits 
  //

  // Create the digits manager
  if (fDigitsManager) {
    delete fDigitsManager;
  }
  fDigitsManager = new AliTRDdigitsManager();
  fDigitsManager->SetDebug(fDebug);

  // Open the input file
  return fDigitsManager->Open(name);

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digit2Raw(const Char_t *name1, const Char_t *name2)
{
  //
  // Convert the digits to raw data byte stream. The output is written
  // into the the binary files <name1> and <name2>.
  //
  // The pseudo raw data format is currently defined like this:
  //
  //          LDC header (8 bytes)
  //                  FLAG
  //                  LDC no.
  //                  Number of detectors with data (not yet implemented)
  //                  5 empty bytes
  //
  //          Subevent (= single chamber) header (8 bytes)
  //                  FLAG
  //                  Detector number (2 bytes)
  //                  Number of data bytes (2 bytes)
  //                  Number of pads with data (2 bytes)
  //                  1 empty byte
  //
  //          Data bank
  //

  const Int_t kLDCHeaderLength      = 8;
  const Int_t kSubeventHeaderLength = 8;

  const Int_t kLDCDummyFlag         = 0xAA;
  const Int_t kSubeventDummyFlag    = 0xBB;

  int headerLDC[2];
  int headerSubevent[2];

  Int_t          ntotalbyte[2] = { 0 };
  Int_t          nbyte = 0;
  Int_t          npads = 0;
  unsigned char *byte_p;
  unsigned char *header_p;

  AliTRDgeometryFull *geo = new AliTRDgeometryFull();
  AliTRDparameter    *par = new AliTRDparameter("TRDparameter"
                                               ,"TRD parameter class");
  AliTRDdataArrayI   *digits;

  if (fDebug) {
    Info("Digit2Raw","Open the LDC output files %s, %s"
        ,name1,name2);
  }
  ofstream *outputFile1 = new ofstream(name1, ios::out | ios::binary);
  ofstream *outputFile2 = new ofstream(name2, ios::out | ios::binary);
  ofstream *outputFile;

  if (!fDigitsManager) {
    Error("Digit2Raw","No input file open\n");
    return kFALSE;
  }

  // Read in the digit arrays
  if (!fDigitsManager->ReadDigits()) {
    return kFALSE;
  }

  // Count the number of chambers with data
  Int_t ndetLDC0 = 0;
  Int_t ndetLDC1 = 0;

  if (fDebug > 1) {
    Info("Digit2Raw","Write the LDC headers");
  }

  // Write the LDC header 1
  byte_p    = (unsigned char *) headerLDC;
  header_p  = byte_p;
  *byte_p++ = kLDCDummyFlag;
  *byte_p++ = 1;
  *byte_p++ = ndetLDC0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  outputFile1->write(header_p,kLDCHeaderLength);
  ntotalbyte[0] += kLDCHeaderLength;

  if (fDebug > 1) {
    Info("Digit2Raw","LDC header 0 = %d, %d",headerLDC[0],headerLDC[1]);
  }

  // Write the LDC header 1
  byte_p    = (unsigned char *) headerLDC;
  header_p  = byte_p;
  *byte_p++ = kLDCDummyFlag;
  *byte_p++ = 2;
  *byte_p++ = ndetLDC1;
  *byte_p++ = 0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  *byte_p++ = 0;
  outputFile2->write(header_p,kLDCHeaderLength);
  ntotalbyte[1] += kLDCHeaderLength;

  if (fDebug > 1) {
    Info("Digit2Raw","LDC header 1 = %d, %d",headerLDC[0],headerLDC[1]);
  }

  // Loop through all detectors
  for (Int_t det = 0; det < AliTRDgeometry::Ndet(); det++) {

    Int_t cham      = geo->GetChamber(det);
    Int_t plan      = geo->GetPlane(det);
    Int_t sect      = geo->GetSector(det);
    Int_t rowMax    = par->GetRowMax(plan,cham,sect);
    Int_t colMax    = par->GetColMax(plan);
    Int_t timeMax   = par->GetTimeMax();
    Int_t bufferMax = rowMax*colMax*timeMax;
    int  *buffer    = new int[bufferMax];

    npads  = 0;
    nbyte  = 0;
    byte_p = (unsigned char *) buffer;

    // Determine the LDC (resp. output file)
    Int_t ldc;
    if (sect < 9) {
      outputFile = outputFile1;
      ldc = 0;
    }
    else {
      outputFile = outputFile2;
      ldc = 1;
    }

    // Get the digits array
    digits = fDigitsManager->GetDigits(det);
    digits->Expand();

    // Loop through the detector pixel
    for (Int_t col = 0; col < colMax; col++) {
      for (Int_t row = 0; row < rowMax; row++) {

	// Check whether data exists for this pad
        Bool_t dataflag = kFALSE;
        for (Int_t time = 0; time < timeMax; time++) {
          Int_t data = digits->GetDataUnchecked(row,col,time);
          if (data) {
            dataflag = kTRUE;
            break;
	  }
	}

        if (dataflag) {

          npads++;

	  // The pad row number
          *byte_p++ = row + 1;
	  // The pad column number
          *byte_p++ = col + 1;
          nbyte += 2;

          Int_t nzero = 0;
          for (Int_t time = 0; time < timeMax; time++) {

            Int_t data = digits->GetDataUnchecked(row,col,time);

            if (!data) {
              nzero++;
              if ((nzero ==       256) || 
                  (time  == timeMax-1)) {
                *byte_p++ = 0;
                *byte_p++ = nzero-1;
                nbyte += 2;
                nzero  = 0;
      	      }
	    }
            else {
              if (nzero) {
                *byte_p++ = 0;
                *byte_p++ = nzero-1;
                nbyte += 2;
                nzero  = 0;
	      }
              // High byte (MSB always set)
              *byte_p++ = ((data >> 8) | 128);
              // Low byte
              *byte_p++ = (data & 0xff);
              nbyte += 2;
	    }

	  }

	}

      }

    }

    // Fill the end of the buffer with zeros
    while (nbyte % 4) {  
      *byte_p++ = 0;
      nbyte++;
    }

    if (fDebug > 1) {
      Info("Digit2Raw","LDC = %d, det = %d, nbyte = %d (%d)",ldc,det,nbyte,bufferMax);
    }

    // Write the subevent header
    byte_p    = (unsigned char *) headerSubevent;
    header_p  = byte_p;
    *byte_p++ = kSubeventDummyFlag;
    *byte_p++ = (det   & 0xff);
    *byte_p++ = (det   >> 8);
    *byte_p++ = (nbyte & 0xff);
    *byte_p++ = (nbyte >> 8);
    *byte_p++ = (npads & 0xff);
    *byte_p++ = (npads >> 8);
    *byte_p++ = 0;
    outputFile->write(header_p,kSubeventHeaderLength);

    // Write the buffer to the file
    byte_p = (unsigned char *) buffer;
    outputFile->write(byte_p,nbyte);

    ntotalbyte[ldc] += nbyte + kSubeventHeaderLength;

    delete buffer;

  }

  if (fDebug) {
    Info("Digit2Raw","Total size: LDC0 = %d, LDC1 = %d",ntotalbyte[0],ntotalbyte[1]);
  }

  outputFile1->close();
  outputFile2->close();

  delete geo;
  delete par;
  delete outputFile1;
  delete outputFile2;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Raw2Digit(const Char_t *name1, const Char_t *name2)
{

  const Int_t  kLDCHeaderLength      = 8;
  const Int_t  kSubeventHeaderLength = 8;

  const Int_t  kNldc = 2;
  const Char_t *name = 0;

  int headerLDC[2];
  int headerSubevent[2];

  Int_t             npads     = 0;
  Int_t             nbyte     = 0;
  unsigned char    *byte_p;
  ifstream         *inputFile = 0;
  AliTRDdataArrayI *digits    = 0;

  AliTRDgeometryFull *geo = new AliTRDgeometryFull();
  AliTRDparameter    *par = new AliTRDparameter("TRDparameter"
                                               ,"TRD parameter class");

  // Create the digits manager
  if (fDigitsManager) {
    delete fDigitsManager;
  }
  fDigitsManager = new AliTRDdigitsManager();
  fDigitsManager->SetDebug(fDebug);
  fDigitsManager->CreateArrays();

  for (Int_t ldc = 0; ldc < kNldc; ldc++) {

    if      (ldc == 0) {
      name = name1;
    }
    else if (ldc == 1) {
      name = name2;
    }
    if (fDebug) {
      Info("Raw2Digit","Open the LDC input file %s",name);
    }
    inputFile = new ifstream(name, ios::in | ios::binary);

    // Read the LDC header
    byte_p = (unsigned char *) headerLDC;
    inputFile->read(byte_p,kLDCHeaderLength);

    if (fDebug > 1) {
      Info("Raw2Digit","LDC header no. %d:",ldc);
      Info("Raw2Digit","\tflag   = %d",*byte_p++);
      Info("Raw2Digit","\tldc no = %d",*byte_p++);
      Info("Raw2Digit","\tndet   = %d",*byte_p++);
      Info("Raw2Digit","\tempty  = %d",*byte_p++);
      Info("Raw2Digit","\tempty  = %d",*byte_p++);
      Info("Raw2Digit","\tempty  = %d",*byte_p++);
      Info("Raw2Digit","\tempty  = %d",*byte_p++);
      Info("Raw2Digit","\tempty  = %d",*byte_p++);
    }

    // Loop through the subevents
    byte_p = (unsigned char *) headerSubevent;
    while (inputFile->read(byte_p,kSubeventHeaderLength)) {

      Int_t flag   = *byte_p++;
      Int_t detl   = *byte_p++;
      Int_t deth   = *byte_p++;
      Int_t det    = detl   + (deth   << 8);
      Int_t nbytel = *byte_p++;
      Int_t nbyteh = *byte_p++;
            nbyte  = nbytel + (nbyteh << 8);
      Int_t npadsl = *byte_p++;
      Int_t npadsh = *byte_p++;
            npads  = npadsl + (npadsh << 8);
      if (fDebug > 2) {
        Info("Raw2Digit","Subevent header:");
        Info("Raw2Digit","\tflag  = %d",flag);
        Info("Raw2Digit","\tdet   = %d",det);
        Info("Raw2Digit","\tnbyte = %d",nbyte);
        Info("Raw2Digit","\tnpads = %d",npads);
        Info("Raw2Digit","\tempty = %d",*byte_p++);
      }      

      // Create the data buffer
      Int_t cham      = geo->GetChamber(det);
      Int_t plan      = geo->GetPlane(det);
      Int_t sect      = geo->GetSector(det);
      Int_t rowMax    = par->GetRowMax(plan,cham,sect);
      Int_t colMax    = par->GetColMax(plan);
      Int_t timeMax   = par->GetTimeMax();
      Int_t bufferMax = rowMax*colMax*timeMax;
      int   *buffer   = new int[bufferMax];
      byte_p          = (unsigned char *) buffer;      
      memset(buffer,0,bufferMax*sizeof(int));

      // Add a container for the digits of this detector
      digits = fDigitsManager->GetDigits(det);
      // Allocate memory space for the digits buffer
      if (digits->GetNtime() == 0) {
        digits->Allocate(rowMax,colMax,timeMax);
      }

      // Read the data   
      inputFile->read(byte_p,nbyte);

      Int_t time;
      Int_t nzero;
      Int_t data;
      Int_t low;
      Int_t high;
      Int_t signal;

      // Decompress the data
      while (nbyte > 0) {

        // The pad row number
        Int_t row = (*byte_p++) - 1;
        // The pad column number
        Int_t col = (*byte_p++) - 1;
        nbyte -= 2;

        time = nzero = 0;

        while ((time  < timeMax) &&
               (nbyte >       0)) {

          data = *byte_p++;
          nbyte--;

          if (data) {
            if (!nzero) {
              // signal for given timebim
              low    = *byte_p++;
              high   = data & 127;
              signal = low + (high << 8);
              if ((row <       0) || (col <       0) || (time <        0) ||
                  (row >= rowMax) || (col >= colMax) || (time >= timeMax)) {
                Error("Raw2Digit"
                     ,"row=%d(%d) col=%d(%d) time=%d(%d)"
		     ,row,rowMax,col,colMax,time,timeMax);
	      }
              else {
                digits->SetDataUnchecked(row,col,time,signal);
	      }
              nbyte--;
              time++;
	    }
            else {
              time += data + 1;
              nzero = 0;
	    }
	  }
          else {
            if (!nzero) {
              nzero = 1;
	    }
            else {
              time++;
              nzero = 0;
	    }
	  }

	}

      }

      digits->Compress(1,0);

      delete buffer;

      byte_p = (unsigned char *) headerSubevent;

    } 

    inputFile->close();
    delete inputFile;
    inputFile = 0;

  }

  delete geo;
  delete par;

  return kTRUE;

}
