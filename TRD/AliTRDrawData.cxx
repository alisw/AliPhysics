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
#include <TMath.h>

#include "AliDAQ.h"
#include "AliRawDataHeader.h"
#include "AliRawReader.h"
#include "AliLog.h"

#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDRawStream.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDrawData)

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData()
  :TObject()
  ,fRawVersion(1)    // Default Raw Data version set here
  ,fCommonParam(0)
  ,fCalibration(0)
  ,fGeo(0)
  ,fNumberOfDDLs(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData(const AliTRDrawData &r)
  :TObject(r)
  ,fRawVersion(1)    // Default Raw Data version set here
  ,fCommonParam(0)
  ,fCalibration(0)
  ,fGeo(0)
  ,fNumberOfDDLs(0)
{
  //
  // Copy constructor
  //

}

//_____________________________________________________________________________
AliTRDrawData::~AliTRDrawData()
{
  //
  // Destructor
  //

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::SetRawVersion(Int_t v)
{
  //
  // Set the raw data version (Currently only version 0, 1 and 2 are available)
  //

  if ( (v >= 0) && (v <= 2) ) {
    fRawVersion = v;
    return kTRUE;
  }

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digits2Raw(TTree *digitsTree, TTree *tracks )
{
  //
  // Initialize necessary parameters and call one
  // of the raw data simulator selected by SetRawVersion.
  //
  // Currently tracklet output is not spported yet and it
  // will be supported in higher version simulator.
  //

  fNumberOfDDLs = AliDAQ::NumberOfDdls("TRD");

  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();

  if (!digitsManager->ReadDigits(digitsTree)) {
    delete digitsManager;
    return kFALSE;
  }

  if (tracks != NULL) {
    delete digitsManager;
    AliError("Tracklet input is not supported yet.");
    return kFALSE;
  }

  fGeo = new AliTRDgeometry();

  fCommonParam = AliTRDCommonParam::Instance();
  if (!fCommonParam) {
    AliError("Could not get common params");
    delete fGeo;
    delete digitsManager;
    return kFALSE;
  }

  fCalibration = AliTRDcalibDB::Instance();
  if (!fCalibration) {
    AliError("Could not get calibration object");
    delete fGeo;
    delete digitsManager;
    return kFALSE;
  }

  Int_t retval = kTRUE;

  // Call appropriate Raw Simulator
  if      ( fRawVersion == 0 )                    retval = Digits2RawV0(digitsManager); 
  else if ( fRawVersion > 0 && fRawVersion <= 2 ) retval = Digits2RawVx(digitsManager); 
  else {
    retval = kFALSE;
    AliWarning(Form("Unsupported raw version (fRawVersion=%d).",fRawVersion));
  }

  // Cleanup
  delete fGeo;
  delete digitsManager;

  return retval;

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digits2RawV0(AliTRDdigitsManager* digitsManager)
{
  //
  // Bogdan's raw simulator (offline use only)
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

  const Int_t kSubeventHeaderLength = 8;
  const Int_t kSubeventDummyFlag    = 0xBB;
  Int_t       headerSubevent[3];

  ofstream     **outputFile = new ofstream* [fNumberOfDDLs];
  UInt_t        *bHPosition = new UInt_t    [fNumberOfDDLs];
  Int_t         *ntotalbyte = new Int_t     [fNumberOfDDLs];
  Int_t          nbyte = 0;
  Int_t          npads = 0;
  unsigned char *bytePtr;
  unsigned char *headerPtr;

  AliTRDdataArrayI *digits;
  AliRawDataHeader  header;   // The event header

  // Open the output files
  for (Int_t iDDL = 0; iDDL < fNumberOfDDLs; iDDL++) {

    char name[20];
    sprintf(name, "TRD_%d.ddl", iDDL + AliTRDRawStream::kDDLOffset);
#ifndef __DECCXX
    outputFile[iDDL] = new ofstream(name, ios::binary);
#else
    outputFile[iDDL] = new ofstream(name);
#endif

    // Write a dummy data header
    bHPosition[iDDL] = outputFile[iDDL]->tellp();
    outputFile[iDDL]->write((char *) (& header),sizeof(header));
    ntotalbyte[iDDL] = 0;

  }

  // Loop through all detectors
  for (Int_t det = 0; det < AliTRDgeometry::Ndet(); det++) {

    Int_t cham      = fGeo->GetChamber(det);
    Int_t plan      = fGeo->GetPlane(det);
    Int_t sect      = fGeo->GetSector(det);
    Int_t rowMax    = fCommonParam->GetRowMax(plan,cham,sect);
    Int_t colMax    = fCommonParam->GetColMax(plan);
    Int_t timeTotal = fCalibration->GetNumberOfTimeBins();
    Int_t bufferMax = rowMax * colMax * timeTotal;
    Int_t *buffer   = new Int_t[bufferMax];

    npads   = 0;
    nbyte   = 0;
    bytePtr = (unsigned char *) buffer;

    Int_t iDDL = sect;

    // Get the digits array
    digits = digitsManager->GetDigits(det);
    digits->Expand();
    // This is to take care of switched off super modules
    if (digits->GetNtime() == 0) {
      continue;
    }

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

    AliDebug(1,Form("det = %d, nbyte = %d (%d)",det,nbyte,bufferMax));

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
    outputFile[iDDL]->write((char *) headerPtr,kSubeventHeaderLength);

    // Write the buffer to the file
    bytePtr = (unsigned char *) buffer;
    outputFile[iDDL]->write((char *) bytePtr,nbyte);

    ntotalbyte[iDDL] += nbyte + kSubeventHeaderLength;

    delete buffer;

  }

  // Update the data headers and close the output files
  for (Int_t iDDL = 0; iDDL < fNumberOfDDLs; iDDL++) {

    header.fSize = UInt_t(outputFile[iDDL]->tellp()) - bHPosition[iDDL];
    header.SetAttribute(0);  // valid data
    outputFile[iDDL]->seekp(bHPosition[iDDL]);
    outputFile[iDDL]->write((char *) (&header),sizeof(header));

    outputFile[iDDL]->close();
    delete outputFile[iDDL];

  }

  delete [] outputFile;
  delete [] bHPosition;
  delete [] ntotalbyte;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digits2RawVx(AliTRDdigitsManager *digitsManager)
{
  //
  // Raw data simulator for all versions > 0. This is prepared for real data.
  // This version simulate only raw data with ADC data and not with tracklet.
  //

  const Int_t kMaxHcWords = (fGeo->TBmax()/3)*fGeo->ADCmax()*fGeo->MCMmax()*fGeo->ROBmaxC1()/2 + 100 + 20;

  // Buffer to temporary store half chamber data
  UInt_t     *hc_buffer   = new UInt_t[kMaxHcWords];

  // sect is same as iDDL, so I use only sect here.
  for (Int_t sect = 0; sect < fGeo->Nsect(); sect++) { 

    char name[1024];
    sprintf(name,"TRD_%d.ddl",sect + AliTRDRawStream::kDDLOffset);

#ifndef __DECCXX
    ofstream *of = new ofstream(name, ios::binary);
#else
    ofstream *of = new ofstream(name);
#endif

    // Write a dummy data header
    AliRawDataHeader  header;  // the event header
    UInt_t hpos = of->tellp();
    of->write((char *) (& header), sizeof(header));

    // Reset payload byte size (payload does not include header).
    Int_t npayloadbyte = 0;

    // GTU common data header (5x4 bytes per super module, shows link mask)
    for( Int_t cham = 0; cham < fGeo->Ncham(); cham++ ) {
      UInt_t GtuCdh = (0xe << 28);
      for( Int_t plan = 0; plan < fGeo->Nplan(); plan++) {
	Int_t iDet = fGeo->GetDetector(plan, cham, sect);
	// If chamber status is ok, we assume that the optical link is also OK.
        // This is shown in the GTU link mask.
	if ( fCalibration->GetChamberStatus(iDet) )
	  GtuCdh = GtuCdh | (3 << (2*plan));
      }
      of->write((char *) (& GtuCdh), sizeof(GtuCdh));
      npayloadbyte += 4;
    }

    // Prepare chamber data
    for( Int_t cham = 0; cham < fGeo->Ncham(); cham++) {
      for( Int_t plan = 0; plan < fGeo->Nplan(); plan++) {

        Int_t iDet = fGeo->GetDetector(plan,cham,sect);

        // Get the digits array
        AliTRDdataArrayI *digits = digitsManager->GetDigits(iDet);
        digits->Expand();

        Int_t hcwords = 0;

        // Process A side of the chamber
	if ( fRawVersion >= 1 && fRawVersion <= 2 ) hcwords = ProduceHcDataV1andV2(digits,0,iDet,hc_buffer,kMaxHcWords);
        of->write((char *) hc_buffer, hcwords*4);
        npayloadbyte += hcwords*4;

        // Process B side of the chamber
	if ( fRawVersion >= 1 && fRawVersion <= 2 ) hcwords = ProduceHcDataV1andV2(digits,1,iDet,hc_buffer,kMaxHcWords);
        of->write((char *) hc_buffer, hcwords*4);
        npayloadbyte += hcwords*4;

      }
    }

    // Complete header
    header.fSize = UInt_t(of->tellp()) - hpos;
    header.SetAttribute(0);  // Valid data
    of->seekp(hpos);         // Rewind to header position
    of->write((char *) (& header), sizeof(header));
    of->close();
    delete of;

  }

  delete hc_buffer;
  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDrawData::ProduceHcDataV1andV2(AliTRDdataArrayI *digits, Int_t side
                                        , Int_t det, UInt_t *buf, Int_t maxSize)
{
  //
  // This function simulates: 1) SM-I commissiong data Oct. 06 (fRawVersion == 1).
  //                          2) Full Raw Production Version   (fRawVersion == 2)
  //
  // Produce half chamber data (= an ORI data) for the given chamber (det) and side (side)
  // where
  //
  //   side=0 means A side with ROB positions 0, 2, 4, 6.
  //   side=1 means B side with ROB positions 1, 3, 5, 7.
  //
  // Chamber type (C0 orC1) is determined by "det" automatically.
  // Appropriate size of buffer (*buf) must be prepared prior to calling this function.
  // Pointer to the buffer and its size must be given to "buf" and "maxSize".
  // Return value is the number of valid data filled in the buffer in unit of 32 bits
  // UInt_t words.
  // If buffer size if too small, the data is truncated with the buffer size however
  // the function will finish without crash (this behaviour is similar to the MCM).
  //

  Int_t          nw = 0;                       // Number of written    words
  Int_t          of = 0;                       // Number of overflowed words
  Int_t        plan = fGeo->GetPlane( det );   // Plane
  Int_t        cham = fGeo->GetChamber( det ); // Chamber
  Int_t        sect = fGeo->GetSector( det );  // Sector (=iDDL)
  Int_t        nRow = fCommonParam->GetRowMax( plan, cham, sect );
  Int_t        nCol = fCommonParam->GetColMax( plan );
  const Int_t nTBin = fCalibration->GetNumberOfTimeBins();
  Int_t      kCtype = 0;                       // Chamber type (0:C0, 1:C1)
  Int_t         iEv = 0xA;                     // Event ID. Now fixed to 10, how do I get event id?
  UInt_t          x = 0;                       // General used number

  // Check the nCol and nRow.
  if ((nCol == 144) && 
      (nRow == 16 || nRow == 12)) {
    kCtype = (nRow-12) / 4;
  } 
  else {
    AliError(Form("This type of chamber is not supported (nRow=%d, nCol=%d)."
                 ,nRow,nCol));
    return 0;
  }

  AliDebug(1,Form("Producing raw data for sect=%d plan=%d cham=%d side=%d"
                 ,sect,plan,cham,side));

  // Tracklet should be processed here but not implemented yet

  // Write end of tracklet marker
  if (nw < maxSize) {
    buf[nw++] = end_of_tracklet_marker;
  } 
  else {
    of++;
  }

  // Half Chamber header
  if      ( fRawVersion == 1 ) {
    // Now it is the same version as used in SM-I commissioning.
    Int_t  dcs = det+100;      // DCS Serial (in simulation, it is meaningless
    x = (dcs<<20) | (sect<<15) | (plan<<12) | (cham<<9) | (side<<8) | 1;
    if (nw < maxSize) {
      buf[nw++] = x; 
    }
    else {
      of++;
    }
  } 
  else if ( fRawVersion == 2 ) {
    // h[0] (there are 2 HC header)
    Int_t minorv = 0;      // The minor version number
    Int_t add    = 1;      // The number of additional header words to follow
    x = (1<<31) | (fRawVersion<<24) | (minorv<<17) | (add<<14) | (sect<<9) | (plan<<6) | (cham<<3) | (side<<2) | 1;
    if (nw < maxSize) {
      buf[nw++] = x; 
    }
    else {
      of++;
    }
    // h[1]
    Int_t bc_ctr   = 99; // bunch crossing counter. Here it is set to 99 always for no reason
    Int_t pt_ctr   = 15; // pretrigger counter. Here it is set to 15 always for no reason
    Int_t pt_phase = 11; // pretrigger phase. Here it is set to 11 always for no reason
    x = (bc_ctr<<16) | (pt_ctr<<12) | (pt_phase<<8) | ((nTBin-1)<<2) | 1;
    if (nw < maxSize) {
      buf[nw++] = x; 
    }
    else {
      of++;
    }
  }

  // Scan for ROB and MCM
  for (Int_t iRobRow = 0; iRobRow < (kCtype + 3); iRobRow++ ) {
    Int_t iRob = iRobRow * 2 + side;
    for (Int_t iMcm = 0; iMcm < fGeo->MCMmax(); iMcm++ ) {
      Int_t padrow = iRobRow * 4 + iMcm / 4;

      // MCM header
      x = ((iRob * fGeo->MCMmax() + iMcm) << 24) | ((iEv % 0x100000) << 4) | 0xC;
      if (nw < maxSize) {
        buf[nw++] = x; 
      }
      else {
        of++;
      }

      // ADC data
      for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
	Int_t padcol = fGeo->GetPadCol(iRob, iMcm, iAdc);
	UInt_t aa = !(iAdc & 1) + 2;
        UInt_t *a = new UInt_t[nTBin+2];
        // 3 timebins are packed into one 32 bits word
        for (Int_t iT = 0; iT < nTBin; iT+=3) { 
          if ((padcol >=    0) && (padcol <  nCol)) {
	    a[iT  ] = ((iT    ) < nTBin ) ? digits->GetDataUnchecked(padrow,padcol,iT    ) : 0;
	    a[iT+1] = ((iT + 1) < nTBin ) ? digits->GetDataUnchecked(padrow,padcol,iT + 1) : 0;
	    a[iT+2] = ((iT + 2) < nTBin ) ? digits->GetDataUnchecked(padrow,padcol,iT + 2) : 0; 
	  } 
	  else {
	    a[iT] = a[iT+1] = a[iT+2] = 0; // This happenes at the edge of chamber (should be pedestal! How?)
	  }
	  x = (a[iT+2] << 22) | (a[iT+1] << 12) | (a[iT] << 2) | aa;
	  if (nw < maxSize) {
	    buf[nw++] = x; 
	  }
	  else {
	    of++;
	  }
	}
        // Diagnostics
        Float_t avg = 0;
        Float_t rms = 0;
        for (Int_t iT = 0; iT < nTBin; iT++) {
          avg += (Float_t) (a[iT]);
	}
        avg /= (Float_t) nTBin;
        for (Int_t iT = 0; iT < nTBin; iT++) {
          rms += ((Float_t) (a[iT]) - avg) * ((Float_t) (a[iT]) - avg);
	}
        rms = TMath::Sqrt(rms / (Float_t) nTBin);
        if (rms > 1.7) {
          AliDebug(2,Form("Large RMS (>1.7)  (ROB,MCM,ADC)=(%02d,%02d,%02d), avg=%03.1f, rms=%03.1f"
			  ,iRob,iMcm,iAdc,avg,rms));
	}
        delete a;
      }
    }
  }

  // Write end of raw data marker
  if (nw < maxSize) {
    buf[nw++] = end_of_event_marker; 
  }
  else {
    of++;
  }
  if (of != 0) {
    AliWarning("Buffer overflow. Data is truncated. Please increase buffer size and recompile.");
  }

  return nw;

}

//_____________________________________________________________________________
AliTRDdigitsManager *AliTRDrawData::Raw2Digits(AliRawReader *rawReader)
{
  //
  // Read raw data and convert to digits
  //

  if ( fRawVersion == 0 ) {
    return Raw2DigitsV0(rawReader);  // fRawVersion == 0
  }
  else {
    return Raw2DigitsVx(rawReader);  // fRawVersion > 0
  }

}

//_____________________________________________________________________________
AliTRDdigitsManager *AliTRDrawData::Raw2DigitsV0(AliRawReader *rawReader)
{
  //
  // Bogdan's raw data reader (for offline only).
  //

  AliTRDdataArrayI *digits = 0;
  AliTRDdataArrayI *track0 = 0;
  AliTRDdataArrayI *track1 = 0;
  AliTRDdataArrayI *track2 = 0; 

  AliTRDgeometry *geo = new AliTRDgeometry();

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliError("Could not get common parameters");
    return 0;
  }

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("Could not get calibration object");
    return 0;
  }

  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
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

      AliDebug(2,"Subevent header:");
      AliDebug(2,Form("\tdet   = %d",det));
      AliDebug(2,Form("\tnpads = %d",npads));

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
                             input.GetTime(),                0);
    track1->SetDataUnchecked(input.GetRow(),input.GetColumn(),
                             input.GetTime(),                0);
    track2->SetDataUnchecked(input.GetRow(),input.GetColumn(),
                             input.GetTime(),                0);

  }

  if (digits) digits->Compress(1,0);
  if (track0) track0->Compress(1,0);
  if (track1) track1->Compress(1,0);
  if (track2) track2->Compress(1,0);

  delete geo;

  return digitsManager;

}

//_____________________________________________________________________________
AliTRDdigitsManager *AliTRDrawData::Raw2DigitsVx(AliRawReader *rawReader)
{

  //
  // This is executed for all Raw Data Versions > 0. Raw data is read and filled
  // into digits array. Next function is not used.
  //

  AliTRDdataArrayI  *digits      = 0;

  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliError("Could not get common params");
    return 0;
  }

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("Could not get calibration object");
    return 0;
  }

  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->CreateArrays();

  AliTRDRawStream input(rawReader, digitsManager, digits);
  input.SetRawVersion( fRawVersion );
  input.ReadAll();     // Loop through the digits

  delete digits;

  return digitsManager;

}
