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
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDRawStream.h"
#include "AliRawDataHeader.h"
#include "AliRawReader.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"
#include "AliDAQ.h"

ClassImp(AliTRDrawData)

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData():TObject()
{
  //
  // Default constructor
  //

  fDebug         = 0;

}

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData(const AliTRDrawData &r):TObject()
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
void AliTRDrawData::Copy(TObject &r) const
{
  //
  // Copy function
  //

  ((AliTRDrawData &) r).fDebug         = fDebug;

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digits2Raw(TTree *digitsTree)
{
  //
  // Convert the digits to raw data byte stream. The output is written
  // into the the binary files TRD_<DDL number>.ddl.
  //
  // The pseudo raw data format is currently defined like this:
  //
  //          DDL data header
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

  const Int_t kNumberOfDDLs         = AliDAQ::NumberOfDdls("TRD");
  const Int_t kSubeventHeaderLength = 8;
  const Int_t kSubeventDummyFlag    = 0xBB;
  Int_t       headerSubevent[3];

  ofstream     **outputFile = new ofstream* [kNumberOfDDLs];
  UInt_t        *bHPosition = new UInt_t    [kNumberOfDDLs];
  Int_t         *ntotalbyte = new Int_t     [kNumberOfDDLs];
  Int_t          nbyte = 0;
  Int_t          npads = 0;
  unsigned char *bytePtr;
  unsigned char *headerPtr;

  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->SetDebug(fDebug);

  // Read in the digit arrays
  if (!digitsManager->ReadDigits(digitsTree)) {
    delete digitsManager;
    return kFALSE;
  }

  AliTRDgeometry   *geo = new AliTRDgeometry();
  AliTRDdataArrayI *digits;

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    printf("<AliTRDrawData::Digits2Raw> ");
    printf("Could not get common params\n");
    return 0;
  }
  
  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    printf("<AliTRDdigitizer::Digits2Raw> ");
    printf("Could not get calibration object\n");
    return kFALSE;
  }
    
  // the event header
  AliRawDataHeader header;

  // Open the output files
  for (Int_t iDDL = 0; iDDL < kNumberOfDDLs; iDDL++) {
    char name[20];
    strcpy(name,AliDAQ::DdlFileName("TRD",iDDL));
#ifndef __DECCXX
    outputFile[iDDL] = new ofstream(name, ios::binary);
#else
    outputFile[iDDL] = new ofstream(name);
#endif

    // Write a dummy data header
    bHPosition[iDDL] = outputFile[iDDL]->tellp();
    outputFile[iDDL]->write((char*)(&header),sizeof(header));
    ntotalbyte[iDDL] = 0;
  }

  // Loop through all detectors
  for (Int_t det = 0; det < AliTRDgeometry::Ndet(); det++) {

    Int_t cham      = geo->GetChamber(det);
    Int_t plan      = geo->GetPlane(det);
    Int_t sect      = geo->GetSector(det);
    Int_t rowMax    = commonParam->GetRowMax(plan,cham,sect);
    Int_t colMax    = commonParam->GetColMax(plan);
    Int_t timeTotal = calibration->GetNumberOfTimeBins();
    Int_t bufferMax = rowMax*colMax*timeTotal;
    Int_t *buffer   = new Int_t[bufferMax];

    npads   = 0;
    nbyte   = 0;
    bytePtr = (unsigned char *) buffer;

    Int_t iDDL = sect;

    // Get the digits array
    digits = digitsManager->GetDigits(det);
    digits->Expand();

    // Loop through the detector pixel
    for (Int_t col = 0; col < colMax; col++) {
      for (Int_t row = 0; row < rowMax; row++) {

	// Check whether data exists for this pad
        Bool_t dataflag = kFALSE;
        for (Int_t time = 0; time < timeTotal; time++) {
          Int_t data = digits->GetDataUnchecked(row,col,time);
          if (data) {
            dataflag = kTRUE;
            break;
	  }
	}

        if (dataflag) {

          npads++;

	  // The pad row number
          *bytePtr++ = row + 1;
	  // The pad column number
          *bytePtr++ = col + 1;
          nbyte += 2;

          Int_t nzero = 0;
          for (Int_t time = 0; time < timeTotal; time++) {

            Int_t data = digits->GetDataUnchecked(row,col,time);

            if (!data) {
              nzero++;
              if ((nzero ==       256) || 
                  (time  == timeTotal-1)) {
                *bytePtr++ = 0;
                *bytePtr++ = nzero-1;
                nbyte += 2;
                nzero  = 0;
      	      }
	    }
            else {
              if (nzero) {
                *bytePtr++ = 0;
                *bytePtr++ = nzero-1;
                nbyte += 2;
                nzero  = 0;
	      }
              // High byte (MSB always set)
              *bytePtr++ = ((data >> 8) | 128);
              // Low byte
              *bytePtr++ = (data & 0xff);
              nbyte += 2;
	    }

	  }

	}

      }

    }

    // Fill the end of the buffer with zeros
    while (nbyte % 4) {  
      *bytePtr++ = 0;
      nbyte++;
    }

    if (fDebug > 1) {
      Info("Digits2Raw","det = %d, nbyte = %d (%d)",det,nbyte,bufferMax);
    }

    // Write the subevent header
    bytePtr    = (unsigned char *) headerSubevent;
    headerPtr  = bytePtr;
    *bytePtr++ = kSubeventDummyFlag;
    *bytePtr++ = (det   & 0xff);
    *bytePtr++ = (det   >> 8);
    *bytePtr++ = (nbyte & 0xff);
    *bytePtr++ = (nbyte >> 8);
    *bytePtr++ = (nbyte >> 16);
    *bytePtr++ = (npads & 0xff);
    *bytePtr++ = (npads >> 8);
    outputFile[iDDL]->write((char*)headerPtr,kSubeventHeaderLength);

    // Write the buffer to the file
    bytePtr = (unsigned char *) buffer;
    outputFile[iDDL]->write((char*)bytePtr,nbyte);

    ntotalbyte[iDDL] += nbyte + kSubeventHeaderLength;

    delete buffer;

  }

  if (fDebug) {
    for (Int_t iDDL = 0; iDDL < kNumberOfDDLs; iDDL++) {
      Info("Digits2Raw","Total size: DDL %d = %d",iDDL,ntotalbyte[iDDL]);
    }
  }

  // Update the data headers and close the output files
  for (Int_t iDDL = 0; iDDL < kNumberOfDDLs; iDDL++) {
    header.fSize = UInt_t(outputFile[iDDL]->tellp()) - bHPosition[iDDL];
    header.SetAttribute(0);  // valid data
    outputFile[iDDL]->seekp(bHPosition[iDDL]);
    outputFile[iDDL]->write((char*)(&header),sizeof(header));

    outputFile[iDDL]->close();
    delete outputFile[iDDL];
  }

  delete geo;
  delete digitsManager;

  delete [] outputFile;
  delete [] bHPosition;
  delete [] ntotalbyte;




  return kTRUE;

}

//_____________________________________________________________________________
AliTRDdigitsManager* AliTRDrawData::Raw2Digits(AliRawReader* rawReader)
{
  //
  // Read the raw data digits and put them into the returned digits manager
  //

  AliTRDdataArrayI *digits    = 0;
  AliTRDdataArrayI *track0    = 0;
  AliTRDdataArrayI *track1    = 0;
  AliTRDdataArrayI *track2    = 0; 

  AliTRDgeometry *geo = new AliTRDgeometry();

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    printf("<AliTRDrawData::Raw2Digits> ");
    printf("Could not get common params\n");
    return 0;
  }
    
  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    printf("<AliTRDdigitizer::Raw2Digits> ");
    printf("Could not get calibration object\n");
    return 0;
  }

  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->SetDebug(fDebug);
  digitsManager->CreateArrays();

  AliTRDRawStream input(rawReader);

  // Loop through the digits
  while (input.Next()) {

    Int_t det    = input.GetDetector();
    Int_t npads  = input.GetNPads();

    if (input.IsNewDetector()) {

      if (digits) digits->Compress(1,0);
      if (track0) track0->Compress(1,0);
      if (track1) track1->Compress(1,0);
      if (track2) track2->Compress(1,0);

      if (fDebug > 2) {
	Info("Raw2Digits","Subevent header:");
	Info("Raw2Digits","\tdet   = %d",det);
	Info("Raw2Digits","\tnpads = %d",npads);
      }      

      // Create the data buffer
      Int_t cham      = geo->GetChamber(det);
      Int_t plan      = geo->GetPlane(det);
      Int_t sect      = geo->GetSector(det);
      Int_t rowMax    = commonParam->GetRowMax(plan,cham,sect);
      Int_t colMax    = commonParam->GetColMax(plan);
      Int_t timeTotal = calibration->GetNumberOfTimeBins();

      // Add a container for the digits of this detector
      digits = digitsManager->GetDigits(det);
      track0 = digitsManager->GetDictionary(det,0);
      track1 = digitsManager->GetDictionary(det,1);
      track2 = digitsManager->GetDictionary(det,2);
      // Allocate memory space for the digits buffer
      if (digits->GetNtime() == 0) {
        digits->Allocate(rowMax,colMax,timeTotal);
        track0->Allocate(rowMax,colMax,timeTotal);
        track1->Allocate(rowMax,colMax,timeTotal);
        track2->Allocate(rowMax,colMax,timeTotal);
      }

    } 

    digits->SetDataUnchecked(input.GetRow(),input.GetColumn(),
			     input.GetTime(),input.GetSignal());
    track0->SetDataUnchecked(input.GetRow(),input.GetColumn(),
                             input.GetTime(),               -1);
    track1->SetDataUnchecked(input.GetRow(),input.GetColumn(),
                             input.GetTime(),               -1);
    track2->SetDataUnchecked(input.GetRow(),input.GetColumn(),
                             input.GetTime(),               -1);
  }

  if (digits) digits->Compress(1,0);
  if (track0) track0->Compress(1,0);
  if (track1) track1->Compress(1,0);
  if (track2) track2->Compress(1,0);

  delete geo;

  return digitsManager;

}
