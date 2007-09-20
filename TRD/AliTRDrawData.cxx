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
#include "TClass.h"

#include "AliDAQ.h"
#include "AliRawDataHeader.h"
#include "AliRawReader.h"
#include "AliLog.h"

#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDRawStream.h"
#include "AliTRDRawStreamV2.h"

#include "AliTRDcalibDB.h"
#include "AliFstream.h"

#include "AliTRDSignalIndex.h"

#include "AliTRDfeeParam.h"
#include "AliTRDmcmSim.h"

ClassImp(AliTRDrawData)

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData()
  :TObject()
  ,fGeo(NULL)
  ,fFee(NULL)
  ,fNumberOfDDLs(0)
{
  //
  // Default constructor
  //
  fFee = AliTRDfeeParam::Instance();
  fNumberOfDDLs = AliDAQ::NumberOfDdls("TRD");
}

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData(const AliTRDrawData &r)
  :TObject(r)
  ,fGeo(NULL)
  ,fFee(NULL)
  ,fNumberOfDDLs(0)
{
  //
  // Copy constructor
  //
  fFee = AliTRDfeeParam::Instance();
  fNumberOfDDLs = AliDAQ::NumberOfDdls("TRD");
}

//_____________________________________________________________________________
AliTRDrawData::~AliTRDrawData()
{
  //
  // Destructor
  //
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

  if (!AliTRDcalibDB::Instance()) {
    AliError("Could not get calibration object");
    delete fGeo;
    delete digitsManager;
    return kFALSE;
  }

  Int_t retval = kTRUE;
  Int_t rv     = fFee->GetRAWversion();

  // Call appropriate Raw Simulator
  if ( rv > 0 && rv <= 3 ) retval = Digits2Raw(digitsManager); 
  else {
    retval = kFALSE;
    AliWarning(Form("Unsupported raw version (%d).", rv));
  }

  // Cleanup
  delete fGeo;
  delete digitsManager;

  return retval;

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digits2Raw(AliTRDdigitsManager *digitsManager)
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

    AliFstream* of = new AliFstream(name);

    // Write a dummy data header
    AliRawDataHeader  header;  // the event header
    UInt_t hpos = of->Tellp();
    of->WriteBuffer((char *) (& header), sizeof(header));

    // Reset payload byte size (payload does not include header).
    Int_t npayloadbyte = 0;

    // GTU common data header (5x4 bytes per super module, shows link mask)
    for( Int_t cham = 0; cham < fGeo->Ncham(); cham++ ) {
      UInt_t GtuCdh = (UInt_t)(0xe << 28);
      for( Int_t plan = 0; plan < fGeo->Nplan(); plan++) {
	Int_t iDet = fGeo->GetDetector(plan, cham, sect);
	// If chamber status is ok, we assume that the optical link is also OK.
        // This is shown in the GTU link mask.
	if ( AliTRDcalibDB::Instance()->GetChamberStatus(iDet) )
	  GtuCdh = GtuCdh | (3 << (2*plan));
      }
      of->WriteBuffer((char *) (& GtuCdh), sizeof(GtuCdh));
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
	Int_t rv = fFee->GetRAWversion();

        // Process A side of the chamber
	if ( rv >= 1 && rv <= 2 ) hcwords = ProduceHcDataV1andV2(digits,0,iDet,hc_buffer,kMaxHcWords);
	if ( rv == 3 )            hcwords = ProduceHcDataV3     (digits,0,iDet,hc_buffer,kMaxHcWords);

        of->WriteBuffer((char *) hc_buffer, hcwords*4);
        npayloadbyte += hcwords*4;

        // Process B side of the chamber
	if ( rv >= 1 && rv <= 2 ) hcwords = ProduceHcDataV1andV2(digits,1,iDet,hc_buffer,kMaxHcWords);
	if ( rv >= 3 )            hcwords = ProduceHcDataV3     (digits,1,iDet,hc_buffer,kMaxHcWords);

        of->WriteBuffer((char *) hc_buffer, hcwords*4);
        npayloadbyte += hcwords*4;
      }
    }

    // Complete header
    header.fSize = UInt_t(of->Tellp()) - hpos;
    header.SetAttribute(0);  // Valid data
    of->Seekp(hpos);         // Rewind to header position
    of->WriteBuffer((char *) (& header), sizeof(header));
    delete of;
  }

  delete [] hc_buffer;
  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDrawData::ProduceHcDataV1andV2(AliTRDdataArrayI *digits, Int_t side
                                        , Int_t det, UInt_t *buf, Int_t maxSize)
{
  //
  // This function simulates: 1) SM-I commissiong data Oct. 06 (Raw Version == 1).
  //                          2) Full Raw Production Version   (Raw Version == 2)
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
  Int_t        nRow = fGeo->GetRowMax( plan, cham, sect );
  Int_t        nCol = fGeo->GetColMax( plan );
  const Int_t nTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t      kCtype = 0;                       // Chamber type (0:C0, 1:C1)
  Int_t         iEv = 0xA;                     // Event ID. Now fixed to 10, how do I get event id?
  UInt_t          x = 0;                       // General used number
  Int_t          rv = fFee->GetRAWversion();

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
    buf[nw++] = kEndoftrackletmarker;
  } 
  else {
    of++;
  }

  // Half Chamber header
  if      ( rv == 1 ) {
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
  else if ( rv == 2 ) {
    // h[0] (there are 3 HC header)
    Int_t minorv = 0;      // The minor version number
    Int_t add    = 2;      // The number of additional header words to follow
    x = (1<<31) | (rv<<24) | (minorv<<17) | (add<<14) | (sect<<9) | (plan<<6) | (cham<<3) | (side<<2) | 1;
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
    // h[2]
    Int_t ped_setup       = 1;    // Pedestal filter setup (0:1). Here it is always 1 for no reason
    Int_t gain_setup      = 1;    // Gain filter setup (0:1). Here it is always 1 for no reason
    Int_t tail_setup      = 1;    // Tail filter setup (0:1). Here it is always 1 for no reason
    Int_t xt_setup        = 0;    // Cross talk filter setup (0:1). Here it is always 0 for no reason
    Int_t nonlin_setup    = 0;    // Nonlinearity filter setup (0:1). Here it is always 0 for no reason
    Int_t bypass_setup    = 0;    // Filter bypass (for raw data) setup (0:1). Here it is always 0 for no reason
    Int_t common_additive = 10;   // Digital filter common additive (0:63). Here it is always 10 for no reason
    x = (ped_setup<<31) | (gain_setup<<30) | (tail_setup<<29) | (xt_setup<<28) | (nonlin_setup<<27)
      | (bypass_setup<<26) | (common_additive<<20) | 1;
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
	Int_t padcol = fFee->GetPadColFromADC(iRob, iMcm, iAdc);
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
        delete [] a;
      }
    }
  }

  // Write end of raw data marker
  if (nw < maxSize) {
    buf[nw++] = kEndofrawdatamarker; 
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
Int_t AliTRDrawData::ProduceHcDataV3(AliTRDdataArrayI *digits, Int_t side
				     , Int_t det, UInt_t *buf, Int_t maxSize)
{
  //
  // This function simulates: Raw Version == 3 (Zero Suppression Prototype)
  //

  Int_t          nw = 0;                       // Number of written    words
  Int_t          of = 0;                       // Number of overflowed words
  Int_t        plan = fGeo->GetPlane( det );   // Plane
  Int_t        cham = fGeo->GetChamber( det ); // Chamber
  Int_t        sect = fGeo->GetSector( det );  // Sector (=iDDL)
  Int_t        nRow = fGeo->GetRowMax( plan, cham, sect );
  Int_t        nCol = fGeo->GetColMax( plan );
  const Int_t nTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t      kCtype = 0;                       // Chamber type (0:C0, 1:C1)
  //Int_t         iEv = 0xA;                     // Event ID. Now fixed to 10, how do I get event id?
  UInt_t          x = 0;                       // General used number
  Int_t          rv = fFee->GetRAWversion();

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
    buf[nw++] = kEndoftrackletmarker;
  } 
  else {
    of++;
  }

  // Half Chamber header
  // h[0] (there are 3 HC header)
  Int_t minorv = 0;    // The minor version number
  Int_t add    = 2;    // The number of additional header words to follow
  x = (1<<31) | (rv<<24) | (minorv<<17) | (add<<14) | (sect<<9) | (plan<<6) | (cham<<3) | (side<<2) | 1;
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
  // h[2]
  Int_t ped_setup       = 1;  // Pedestal filter setup (0:1). Here it is always 1 for no reason
  Int_t gain_setup      = 1;  // Gain filter setup (0:1). Here it is always 1 for no reason
  Int_t tail_setup      = 1;  // Tail filter setup (0:1). Here it is always 1 for no reason
  Int_t xt_setup        = 0;  // Cross talk filter setup (0:1). Here it is always 0 for no reason
  Int_t nonlin_setup    = 0;  // Nonlinearity filter setup (0:1). Here it is always 0 for no reason
  Int_t bypass_setup    = 0;  // Filter bypass (for raw data) setup (0:1). Here it is always 0 for no reason
  Int_t common_additive = 10; // Digital filter common additive (0:63). Here it is always 10 for no reason
  x = (ped_setup<<31) | (gain_setup<<30) | (tail_setup<<29) | (xt_setup<<28) | (nonlin_setup<<27)
    | (bypass_setup<<26) | (common_additive<<20) | 1;
  if (nw < maxSize) {
    buf[nw++] = x; 
  }
  else {
    of++;
  }


  // Scan for ROB and MCM
  for (Int_t iRobRow = 0; iRobRow < (kCtype + 3); iRobRow++ ) {
    Int_t iRob = iRobRow * 2 + side;
    for (Int_t iMcm = 0; iMcm < fGeo->MCMmax(); iMcm++ ) {

      AliTRDmcmSim *mcm = new AliTRDmcmSim();
      mcm->Init( det, iRob, iMcm );
      Int_t padrow = mcm->GetRow();

      // Copy ADC data to MCM simulator
      for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
	Int_t padcol = mcm->GetCol( iAdc );
	if ((padcol >=    0) && (padcol <  nCol)) {
	  for (Int_t iT = 0; iT < nTBin; iT++) { 
	    mcm->SetData( iAdc, iT, digits->GetDataUnchecked( padrow, padcol, iT) );
	  } 
	} else {  // this means it is out of chamber, and masked ADC
	  mcm->SetDataPedestal( iAdc );
	}
      }
      // Simulate process in MCM
      mcm->Filter();     // Apply filter
      mcm->ZSMapping();  // Calculate zero suppression mapping
      //mcm->DumpData( "trdmcmdata.txt", "RFZS" ); // debugging purpose

      // Write MCM data to buffer
      Int_t temp_nw =  mcm->ProduceRawStream( &buf[nw], maxSize - nw );
      if( temp_nw < 0 ) {
	of += temp_nw;
	nw += maxSize - nw;
	AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
      } else {
	nw += temp_nw;
      }

      delete mcm;

    }
  }

  // Write end of raw data marker
  if (nw < maxSize) {
    buf[nw++] = kEndofrawdatamarker; 
  }
  else {
    of++;
  }
  if (of != 0) {
    AliError("Buffer overflow. Data is truncated. Please increase buffer size and recompile.");
  }

  return nw;

}

//_____________________________________________________________________________
AliTRDdigitsManager *AliTRDrawData::Raw2Digits(AliRawReader *rawReader)
{
  //
  // Vx of the raw data reading
  //

  AliTRDdataArrayI *digits = 0;
  AliTRDdataArrayI *track0 = 0;
  AliTRDdataArrayI *track1 = 0;
  AliTRDdataArrayI *track2 = 0; 

  AliTRDSignalIndex *indexes = 0;
  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->CreateArrays();

  //AliTRDRawStream input(rawReader);
  AliTRDRawStreamV2 input(rawReader);
  input.SetRawVersion( fFee->GetRAWversion() );
  input.Init();

  AliInfo(Form("Stream version: %s", input.IsA()->GetName()));

  // Loop through the digits
  Int_t lastdet = -1;
  Int_t det    = 0;
  Int_t it = 0;
  while (input.Next()) {

      det    = input.GetDet();

      if (det != lastdet) { // If new detector found
	
	  lastdet = det;

	  if (digits) digits->Compress(1,0);
	  if (track0) track0->Compress(1,0);
	  if (track1) track1->Compress(1,0);
	  if (track2) track2->Compress(1,0);
	
	  // Add a container for the digits of this detector
	  digits = digitsManager->GetDigits(det);
	  track0 = digitsManager->GetDictionary(det,0);
	  track1 = digitsManager->GetDictionary(det,1);
	  track2 = digitsManager->GetDictionary(det,2);

	  // Allocate memory space for the digits buffer
	  if (digits->GetNtime() == 0) 
	    {
	      digits->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	      track0->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	      track1->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	      track2->Allocate(input.GetMaxRow(),input.GetMaxCol(), input.GetNumberOfTimeBins());
	    }

	  indexes = digitsManager->GetIndexes(det);
	  indexes->SetSM(input.GetSM());
	  indexes->SetStack(input.GetStack());
	  indexes->SetLayer(input.GetLayer());
	  indexes->SetDetNumber(det);
	  if (indexes->IsAllocated() == kFALSE)
	    indexes->Allocate(input.GetMaxRow(), input.GetMaxCol(), input.GetNumberOfTimeBins());
	}
    
      // 3 timebin data are stored per word
      for (it = 0; it < 3; it++)
	{
	  if ( input.GetTimeBin() + it < input.GetNumberOfTimeBins() )
	    {
	      if (input.GetSignals()[it] > 0)
		{
		  digits->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, input.GetSignals()[it]);

		  indexes->AddIndexTBin(input.GetRow(), input.GetCol(),
					input.GetTimeBin() + it);
		  track0->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, 0);
		  track1->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, 0);
		  track2->SetDataUnchecked(input.GetRow(), input.GetCol(),
					   input.GetTimeBin() + it, 0);
		}
	    }
	}
  }

  if (digits) digits->Compress(1,0);
  if (track0) track0->Compress(1,0);
  if (track1) track1->Compress(1,0);
  if (track2) track2->Compress(1,0);

  return digitsManager;

}
