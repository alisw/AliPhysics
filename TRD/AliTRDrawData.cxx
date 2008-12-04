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

#include <TMath.h>
#include "TClass.h"

#include "AliDAQ.h"
#include "AliRawDataHeaderSim.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliFstream.h"

#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDrawStreamBase.h"
#include "AliTRDrawOldStream.h"
#include "AliTRDRawStreamV2.h"
#include "AliTRDcalibDB.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDfeeParam.h"
#include "AliTRDmcmSim.h"

ClassImp(AliTRDrawData)

Int_t AliTRDrawData::fgRawFormatVersion = AliTRDrawData::kRawOldFormat;

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

  const Int_t kMaxHcWords = (fGeo->TBmax()/3)
                          * fGeo->ADCmax()
                          * fGeo->MCMmax()
                          * fGeo->ROBmaxC1()/2 
                          + 100 + 20;

  // Buffer to temporary store half chamber data
  UInt_t     *hcBuffer    = new UInt_t[kMaxHcWords];
  
  Bool_t newEvent = kFALSE;  // only for correct readout tree

  // sect is same as iDDL, so I use only sect here.
  for (Int_t sect = 0; sect < fGeo->Nsector(); sect++) { 

    char name[1024];
    sprintf(name,"TRD_%d.ddl",sect + AliTRDrawOldStream::kDDLOffset);

    AliFstream* of = new AliFstream(name);

    // Write a dummy data header
    AliRawDataHeaderSim  header;  // the event header
    UInt_t hpos = of->Tellp();
    of->WriteBuffer((char *) (& header), sizeof(header));

    // Reset payload byte size (payload does not include header).
    Int_t npayloadbyte = 0;

    
	if ( fgRawFormatVersion == 0 ){
    // GTU common data header (5x4 bytes per super module, shows link mask)
    for( Int_t stack = 0; stack < fGeo->Nstack(); stack++ ) {
      UInt_t gtuCdh = (UInt_t)(0xe << 28);
      for( Int_t layer = 0; layer < fGeo->Nlayer(); layer++) {
	Int_t iDet = fGeo->GetDetector(layer, stack, sect);
	
	// If chamber status is ok, we assume that the optical link is also OK.
        // This is shown in the GTU link mask.
	if ( AliTRDcalibDB::Instance()->GetChamberStatus(iDet) )
	  gtuCdh = gtuCdh | (3 << (2*layer));
      }
      of->WriteBuffer((char *) (& gtuCdh), sizeof(gtuCdh));
      npayloadbyte += 4;
    }
	}

    // Prepare chamber data
    for( Int_t stack = 0; stack < fGeo->Nstack(); stack++) {
      for( Int_t layer = 0; layer < fGeo->Nlayer(); layer++) {

        Int_t iDet = fGeo->GetDetector(layer,stack,sect);
	if (iDet == 0) newEvent = kTRUE; // it is expected that each event has at least one tracklet; this is only needed for correct readout tree
	// Get the digits array
	AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(iDet);
        if (digits->HasData() ) {  // second part is new!! and is for indicating a new event

          digits->Expand();

          Int_t hcwords = 0;
	  Int_t rv = fFee->GetRAWversion();


	if ( fgRawFormatVersion == 0 ){
      // Process A side of the chamber
	  if ( rv >= 1 && rv <= 2 ) {
            hcwords = ProduceHcDataV1andV2(digits,0,iDet,hcBuffer,kMaxHcWords);
	  }
	  if ( rv == 3 ) { 
   
	      hcwords = ProduceHcDataV3     (digits,0,iDet,hcBuffer,kMaxHcWords,newEvent);
	      //hcwords = ProduceHcDataV3     (digits,0,iDet,hcBuffer,kMaxHcWords);
	    if(newEvent == kTRUE) newEvent = kFALSE;
	  }

          of->WriteBuffer((char *) hcBuffer, hcwords*4);
          npayloadbyte += hcwords*4;

      // Process B side of the chamber
	  if ( rv >= 1 && rv <= 2 ) {
            hcwords = ProduceHcDataV1andV2(digits,1,iDet,hcBuffer,kMaxHcWords);
	  }
	  if ( rv >= 3 ) {
	   
	      hcwords = ProduceHcDataV3     (digits,1,iDet,hcBuffer,kMaxHcWords,newEvent);
	      //hcwords = ProduceHcDataV3     (digits,1,iDet,hcBuffer,kMaxHcWords);
	  }

          of->WriteBuffer((char *) hcBuffer, hcwords*4);
          npayloadbyte += hcwords*4;

	} else { // real data format
		// Process A side of the chamber
		hcwords = ProduceHcData(digits,0,iDet,hcBuffer,kMaxHcWords,newEvent);
		if(newEvent == kTRUE) newEvent = kFALSE;
        of->WriteBuffer((char *) hcBuffer, hcwords*4);
        npayloadbyte += hcwords*4;
		//for ( Int_t i=0; i<hcwords; i++ ) AliInfo(Form("Buf : %X",hcBuffer[i]));

		// Process B side of the chamber
		hcwords = ProduceHcData(digits,1,iDet,hcBuffer,kMaxHcWords,newEvent);
        of->WriteBuffer((char *) hcBuffer, hcwords*4);
        npayloadbyte += hcwords*4;
	}

	}

      }
    }

    // Complete header
    header.fSize = UInt_t(of->Tellp()) - hpos;
    header.SetAttribute(0);  // Valid data
    of->Seekp(hpos);         // Rewind to header position
    of->WriteBuffer((char *) (& header), sizeof(header));
    delete of;
  }

  delete [] hcBuffer;

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDrawData::ProduceSMIndexData(UInt_t *buf, Int_t& nw){
	// 
	// This function generates 
	// 	 1) SM index words : ssssssss ssssssss vvvv rrrr r d t mmmmm
	// 	    - s : size of SM header (number of header, default = 0x0001)
	// 	    - v : SM header version (default = 0xa)
	// 	    - r : reserved for future use (default = 00000)
	// 	    - d : track data enabled bit (default = 0)
	// 	    - t : tracklet data enabled bit (default = 1)
	// 	    - m : stack mask (each bit corresponds a stack, default = 11111)
	//
	// 	 2) SM header : rrr c vvvv vvvvvvvv vvvv rrrr bbbbbbbb
	// 	    - r : reserved for future use (default = 000)
	// 	    - c : clean check out flag (default = 1)
	// 	    - v : hardware design revision (default = 0x0404)
	// 	    - r : reserved for future use (default = 0x0)
	// 	    - b : physical board ID (default = 0x71)
	//
	// 	 3) stack index words : ssssssss ssssssss vvvv mmmm mmmmmmmm
	// 	    - s : size of stack header (number of header, (default = 0x0007)
	// 	    - v : header version (default = 0xa)
	// 	    - m : link mask (default = 0xfff)
	//
	// 	 4) stack header : vvvvvvvv vvvvvvvv bbbbbbbb rrrr rrr c
	// 	    - v : hardware design revision (default = 0x0404)
	// 	    - b : physical board ID (default = 0x5b)
	// 	    - r : reserved for future use (default = 0000 000)
	// 	    - c : clean checkout flag (default = 1)
	// 	
	// 	 and 6 dummy words(0x00000000)
	//
	
	buf[nw++] = 0x0001a03f;	// SM index words
	buf[nw++] = 0x10404071;	// SM header

	for (Int_t istack=0; istack<5; istack++){
		buf[nw++] = 0x0007afff; // stack index words
		buf[nw++] = 0x04045b01;	// stack header
		for (Int_t i=0;i<6;i++) buf[nw++] = 0x00000000; // 6 dummy words
	} // loop over 5 stacks

}

//_____________________________________________________________________________
Int_t AliTRDrawData::ProduceHcData(AliTRDarrayADC *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize, Bool_t newEvent = kFALSE){
	//
	// This function can be used for both ZS and NZS data
	//

  	Int_t           nw = 0;                       // Number of written    words
  	Int_t           of = 0;                       // Number of overflowed words
  	Int_t        layer = fGeo->GetLayer( det );   // Layer
  	Int_t        stack = fGeo->GetStack( det );   // Stack
  	Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
  	Int_t         nRow = fGeo->GetRowMax( layer, stack, sect );
  	Int_t         nCol = fGeo->GetColMax( layer );
  	const Int_t kNTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
	Int_t       kCtype = 0;                       // Chamber type (0:C0, 1:C1)

	Bool_t tracklet_on = fFee->GetTracklet();     // **new**   

	// Check the nCol and nRow.
  	if ((nCol == 144) && (nRow == 16 || nRow == 12)) {
    	kCtype = (nRow-12) / 4;
  	} else {
    	AliError(Form("This type of chamber is not supported (nRow=%d, nCol=%d).",nRow,nCol));
    	return 0;
  	}

  	AliDebug(1,Form("Producing raw data for sect=%d layer=%d stack=%d side=%d",sect,layer,stack,side));

  	AliTRDmcmSim** mcm = new AliTRDmcmSim*[(kCtype + 3)*(fGeo->MCMmax())];

	if (newEvent) ProduceSMIndexData(buf, nw);		// SM index words , Stack index words

	if (!tracklet_on) WriteIntermediateWordsV2(buf,nw,of,maxSize,det,side);	// no tracklet or NZS
 
  	
  	// scanning direction such, that tracklet-words are sorted in ascending z and then in ascending y order
  	// ROB numbering on chamber and MCM numbering on ROB increase with decreasing z and increasing y
  	for (Int_t iRobRow =  (kCtype + 3)-1; iRobRow >= 0; iRobRow-- ) {
    	Int_t iRob = iRobRow * 2 + side;
    	// MCM on one ROB
    	for (Int_t iMcmRB = 0; iMcmRB < fGeo->MCMmax(); iMcmRB++ ) {
		Int_t iMcm = 16 - 4*(iMcmRB/4 + 1) + (iMcmRB%4);
		Int_t entry = iRobRow*(fGeo->MCMmax()) + iMcm;
	
		mcm[entry] = new AliTRDmcmSim();
		mcm[entry]->Init( det, iRob, iMcm , newEvent);
		//mcm[entry]->Init( det, iRob, iMcm);
		if (newEvent == kTRUE) newEvent = kFALSE; // only one mcm is concerned with new event
		Int_t padrow = mcm[entry]->GetRow();

		// Copy ADC data to MCM simulator
		for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
		    Int_t padcol = mcm[entry]->GetCol( iAdc );
		    if ((padcol >=    0) && (padcol <  nCol)) {
				for (Int_t iT = 0; iT < kNTBin; iT++) { 
				  mcm[entry]->SetData( iAdc, iT, digits->GetData( padrow, padcol, iT) );
				} 
		    } 
	    	else {  // this means it is out of chamber, and masked ADC
				mcm[entry]->SetDataPedestal( iAdc );
	   	 	}
		}

		// Simulate process in MCM
		//		mcm[entry]->Filter();     // Apply filter
		//		mcm[entry]->ZSMapping();  // Calculate zero suppression mapping
		mcm[entry]->CopyArrays();
		mcm[entry]->GeneratefZSM1Dim();
		mcm[entry]->RestoreZeros();

		if (tracklet_on) {
		    mcm[entry]->Tracklet(); 
		    Int_t tempNw =  mcm[entry]->ProduceTrackletStream( &buf[nw], maxSize - nw );
		    //Int_t tempNw = 0;
		    if( tempNw < 0 ) {
				of += tempNw;
				nw += maxSize - nw;
				AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
		    } else {
				nw += tempNw;
		    }
		} else { // no tracklets: write raw-data already in this loop 
		    // Write MCM data to buffer
		    Int_t tempNw =  mcm[entry]->ProduceRawStreamV2( &buf[nw], maxSize - nw );
		    if( tempNw < 0 ) {
				of += tempNw;
				nw += maxSize - nw;
				AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
		    } else {
				nw += tempNw;
		    }
		    
		    delete mcm[entry];
		}

		//mcm->DumpData( "trdmcmdata.txt", "RFZS" ); // debugging purpose
    	}
	}

  	// if tracklets are switched on, raw-data can be written only after all tracklets
	if (tracklet_on) {
      	WriteIntermediateWordsV2(buf,nw,of,maxSize,det,side); 
  
      	// Scan for ROB and MCM
      	for (Int_t iRobRow =  (kCtype + 3)-1; iRobRow >= 0; iRobRow-- ) {
	  	//Int_t iRob = iRobRow * 2 + side;
	  	// MCM on one ROB
	  		for (Int_t iMcmRB = 0; iMcmRB < fGeo->MCMmax(); iMcmRB++ ) {
	      		Int_t iMcm = 16 - 4*(iMcmRB/4 + 1) + (iMcmRB%4);
	      
		      	Int_t entry = iRobRow*(fGeo->MCMmax()) + iMcm; 
		      
		      	// Write MCM data to buffer
		      	Int_t tempNw =  mcm[entry]->ProduceRawStreamV2( &buf[nw], maxSize - nw );
		      	if( tempNw < 0 ) {
			  		of += tempNw;
			  		nw += maxSize - nw;
			  		AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
		      	} else {
			  		nw += tempNw;
		      	}
	      
		      delete mcm[entry];
	  
	  		}	
      	}
	}

  	delete [] mcm;
  
  	// Write end of raw data marker
  	if (nw < maxSize) {
    	buf[nw++] = kEndofrawdatamarker; 
  	} else {
    	of++;
  	}
  	
	if (of != 0) {
    	AliError("Buffer overflow. Data is truncated. Please increase buffer size and recompile.");
  	}

  	return nw;
}


//_____________________________________________________________________________
Int_t AliTRDrawData::ProduceHcDataV1andV2(AliTRDarrayADC *digits, Int_t side
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

  Int_t           nw = 0;                       // Number of written    words
  Int_t           of = 0;                       // Number of overflowed words
  Int_t        layer = fGeo->GetLayer( det );   // Layer
  Int_t        stack = fGeo->GetStack( det );   // Stack
  Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
  Int_t         nRow = fGeo->GetRowMax( layer, stack, sect );
  Int_t         nCol = fGeo->GetColMax( layer );
  const Int_t kNTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t       kCtype = 0;                       // Chamber type (0:C0, 1:C1)
  Int_t          iEv = 0xA;                     // Event ID. Now fixed to 10, how do I get event id?
  UInt_t           x = 0;                       // General used number
  Int_t           rv = fFee->GetRAWversion();

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

  AliDebug(1,Form("Producing raw data for sect=%d layer=%d stack=%d side=%d"
                 ,sect,layer,stack,side));

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
    x = (dcs<<20) | (sect<<15) | (layer<<12) | (stack<<9) | (side<<8) | 1;
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
    x = (1<<31) | (rv<<24) | (minorv<<17) | (add<<14) | (sect<<9) | (layer<<6) | (stack<<3) | (side<<2) | 1;
    if (nw < maxSize) {
      buf[nw++] = x; 
    }
    else {
      of++;
    }
    // h[1]
    Int_t bcCtr   = 99; // bunch crossing counter. Here it is set to 99 always for no reason
    Int_t ptCtr   = 15; // pretrigger counter. Here it is set to 15 always for no reason
    Int_t ptPhase = 11; // pretrigger phase. Here it is set to 11 always for no reason
    x = (bcCtr<<16) | (ptCtr<<12) | (ptPhase<<8) | ((kNTBin-1)<<2) | 1;
    if (nw < maxSize) {
      buf[nw++] = x; 
    }
    else {
      of++;
    }
    // h[2]
    Int_t pedSetup       = 1;    // Pedestal filter setup (0:1). Here it is always 1 for no reason
    Int_t gainSetup      = 1;    // Gain filter setup (0:1). Here it is always 1 for no reason
    Int_t tailSetup      = 1;    // Tail filter setup (0:1). Here it is always 1 for no reason
    Int_t xtSetup        = 0;    // Cross talk filter setup (0:1). Here it is always 0 for no reason
    Int_t nonlinSetup    = 0;    // Nonlinearity filter setup (0:1). Here it is always 0 for no reason
    Int_t bypassSetup    = 0;    // Filter bypass (for raw data) setup (0:1). Here it is always 0 for no reason
    Int_t commonAdditive = 10;   // Digital filter common additive (0:63). Here it is always 10 for no reason
    x = (pedSetup<<31) | (gainSetup<<30) | (tailSetup<<29) | (xtSetup<<28) | (nonlinSetup<<27)
      | (bypassSetup<<26) | (commonAdditive<<20) | 1;
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
        UInt_t *a = new UInt_t[kNTBin+2];
        // 3 timebins are packed into one 32 bits word
        for (Int_t iT = 0; iT < kNTBin; iT+=3) { 
          if ((padcol >=    0) && (padcol <  nCol)) {
	    a[iT  ] = ((iT    ) < kNTBin ) ? digits->GetData(padrow,padcol,iT    ) : 0;
	    a[iT+1] = ((iT + 1) < kNTBin ) ? digits->GetData(padrow,padcol,iT + 1) : 0;
	    a[iT+2] = ((iT + 2) < kNTBin ) ? digits->GetData(padrow,padcol,iT + 2) : 0; 
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
        for (Int_t iT = 0; iT < kNTBin; iT++) {
          avg += (Float_t) (a[iT]);
	}
        avg /= (Float_t) kNTBin;
        for (Int_t iT = 0; iT < kNTBin; iT++) {
          rms += ((Float_t) (a[iT]) - avg) * ((Float_t) (a[iT]) - avg);
	}
        rms = TMath::Sqrt(rms / (Float_t) kNTBin);
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

//Int_t AliTRDrawData::ProduceHcDataV3(AliTRDarrayADC *digits, Int_t side , Int_t det, UInt_t *buf, Int_t maxSize)
Int_t AliTRDrawData::ProduceHcDataV3(AliTRDarrayADC *digits, Int_t side , Int_t det, UInt_t *buf, Int_t maxSize, Bool_t newEvent = kFALSE)
{
  //
  // This function simulates: Raw Version == 3 (Zero Suppression Prototype)
  //

  Int_t           nw = 0;                       // Number of written    words
  Int_t           of = 0;                       // Number of overflowed words
  Int_t        layer = fGeo->GetLayer( det );   // Layer
  Int_t        stack = fGeo->GetStack( det );   // Stack
  Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
  Int_t         nRow = fGeo->GetRowMax( layer, stack, sect );
  Int_t         nCol = fGeo->GetColMax( layer );
  const Int_t kNTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t       kCtype = 0;                       // Chamber type (0:C0, 1:C1)
  //Int_t          iEv = 0xA;                     // Event ID. Now fixed to 10, how do I get event id?

 

  Bool_t tracklet_on = fFee->GetTracklet();     // **new**

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

  AliDebug(1,Form("Producing raw data for sect=%d layer=%d stack=%d side=%d"
                 ,sect,layer,stack,side));

  AliTRDmcmSim** mcm = new AliTRDmcmSim*[(kCtype + 3)*(fGeo->MCMmax())];

  // in case no tracklet-words are processed: write the tracklet-endmarker as well as all additional words immediately and write 
  // raw-data in one go; if tracklet-processing is enabled, first all tracklet-words of a half-chamber have to be processed before the
  // additional words (tracklet-endmarker,headers,...)are written. Raw-data is written in a second loop;
  
  if (!tracklet_on) {
      WriteIntermediateWords(buf,nw,of,maxSize,det,side); 
  }
  
  // Scan for ROB and MCM
  // scanning direction such, that tracklet-words are sorted in ascending z and then in ascending y order
  // ROB numbering on chamber and MCM numbering on ROB increase with decreasing z and increasing y
  for (Int_t iRobRow =  (kCtype + 3)-1; iRobRow >= 0; iRobRow-- ) {
    Int_t iRob = iRobRow * 2 + side;
    // MCM on one ROB
    for (Int_t iMcmRB = 0; iMcmRB < fGeo->MCMmax(); iMcmRB++ ) {
	Int_t iMcm = 16 - 4*(iMcmRB/4 + 1) + (iMcmRB%4);
	Int_t entry = iRobRow*(fGeo->MCMmax()) + iMcm;
	
	mcm[entry] = new AliTRDmcmSim();
	mcm[entry]->Init( det, iRob, iMcm , newEvent);
	//mcm[entry]->Init( det, iRob, iMcm);
	if (newEvent == kTRUE) newEvent = kFALSE; // only one mcm is concerned with new event
	Int_t padrow = mcm[entry]->GetRow();

	// Copy ADC data to MCM simulator
	for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
	    Int_t padcol = mcm[entry]->GetCol( iAdc );
	    if ((padcol >=    0) && (padcol <  nCol)) {
		for (Int_t iT = 0; iT < kNTBin; iT++) { 
		  mcm[entry]->SetData( iAdc, iT, digits->GetData( padrow, padcol, iT) );
		} 
	    } 
	    else {  // this means it is out of chamber, and masked ADC
		mcm[entry]->SetDataPedestal( iAdc );
	    }
	}

	// Simulate process in MCM
	//	mcm[entry]->Filter();     // Apply filter
	//	mcm[entry]->ZSMapping();  // Calculate zero suppression mapping
	mcm[entry]->CopyArrays();
	mcm[entry]->GeneratefZSM1Dim();
	mcm[entry]->RestoreZeros();

	if (tracklet_on) {
	    mcm[entry]->Tracklet(); 
	    Int_t tempNw =  mcm[entry]->ProduceTrackletStream( &buf[nw], maxSize - nw );
	    //Int_t tempNw = 0;
	    if( tempNw < 0 ) {
		of += tempNw;
		nw += maxSize - nw;
		AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
	    } else {
		nw += tempNw;
	    }
	}
	// no tracklets: write raw-data already in this loop
	else {
	     // Write MCM data to buffer
	    Int_t tempNw =  mcm[entry]->ProduceRawStream( &buf[nw], maxSize - nw );
	    if( tempNw < 0 ) {
		of += tempNw;
		nw += maxSize - nw;
		AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
	    } else {
		nw += tempNw;
	    }
	    
	    delete mcm[entry];
	}
	    
	


	//mcm->DumpData( "trdmcmdata.txt", "RFZS" ); // debugging purpose
    }
  }

  // if tracklets are switched on, raw-data can be written only after all tracklets
  if (tracklet_on) {
      WriteIntermediateWords(buf,nw,of,maxSize,det,side); 
  
  
      // Scan for ROB and MCM
      for (Int_t iRobRow =  (kCtype + 3)-1; iRobRow >= 0; iRobRow-- ) {
	  //Int_t iRob = iRobRow * 2 + side;
	  // MCM on one ROB
	  for (Int_t iMcmRB = 0; iMcmRB < fGeo->MCMmax(); iMcmRB++ ) {
	      Int_t iMcm = 16 - 4*(iMcmRB/4 + 1) + (iMcmRB%4);
	      
	      Int_t entry = iRobRow*(fGeo->MCMmax()) + iMcm; 
	      
	      // Write MCM data to buffer
	      Int_t tempNw =  mcm[entry]->ProduceRawStream( &buf[nw], maxSize - nw );
	      if( tempNw < 0 ) {
		  of += tempNw;
		  nw += maxSize - nw;
		  AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
	      } else {
		  nw += tempNw;
	      }
	      
	      delete mcm[entry];
	  
	  }
      }
  }

  delete [] mcm;
  
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

  AliTRDarrayADC *digits = 0;
  AliTRDarrayDictionary *track0 = 0;
  AliTRDarrayDictionary *track1 = 0;
  AliTRDarrayDictionary *track2 = 0;  

  //AliTRDSignalIndex *indexes = 0;
  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->CreateArrays();

  //AliTRDrawOldStream input(rawReader);
  //   AliTRDRawStreamV2 input(rawReader);
  //   input.SetRawVersion( fFee->GetRAWversion() );
  //   input.Init();

  AliTRDrawStreamBase *pinput = AliTRDrawStreamBase::GetRawStream(rawReader);
  AliTRDrawStreamBase &input = *pinput;
  input.SetRawVersion( fFee->GetRAWversion() ); //<= ADDED by MinJung

  AliInfo(Form("Stream version: %s", input.IsA()->GetName()));

  // Loop through the digits
  Int_t det    = 0;

  while (det >= 0)
    {
      det = input.NextChamber(digitsManager);
      if (det >= 0)
	{
	  // get...
	  digits = (AliTRDarrayADC *) digitsManager->GetDigits(det);
	  track0 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,0);
	  track1 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,1);
	  track2 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,2);
	  // and compress
 	  if (digits) digits->Compress();  
 	  if (track0) track0->Compress();   
 	  if (track1) track1->Compress();     
 	  if (track2) track2->Compress();
	}
    }

  delete pinput;
  pinput = NULL;

  return digitsManager;
}


//_____________________________________________________________________________
void AliTRDrawData::WriteIntermediateWords(UInt_t* buf, Int_t& nw, Int_t& of, const Int_t& maxSize, const Int_t& det, const Int_t& side) {
    
    Int_t        layer = fGeo->GetLayer( det );   // Layer
    Int_t        stack = fGeo->GetStack( det );   // Stack
    Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
    Int_t           rv = fFee->GetRAWversion();
    const Int_t kNTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
    UInt_t           x = 0;

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
    x = (1<<31) | (rv<<24) | (minorv<<17) | (add<<14) | (sect<<9) | (layer<<6) | (stack<<3) | (side<<2) | 1;
    if (nw < maxSize) {
	buf[nw++] = x; 
    }
    else {
	of++;
    }
    // h[1]
    Int_t bcCtr   = 99; // bunch crossing counter. Here it is set to 99 always for no reason
    Int_t ptCtr   = 15; // pretrigger counter. Here it is set to 15 always for no reason
    Int_t ptPhase = 11; // pretrigger phase. Here it is set to 11 always for no reason
    x = (bcCtr<<16) | (ptCtr<<12) | (ptPhase<<8) | ((kNTBin-1)<<2) | 1;
    if (nw < maxSize) {
	buf[nw++] = x; 
    }
    else {
	of++;
    }
    // h[2]
    Int_t pedSetup       = 1;  // Pedestal filter setup (0:1). Here it is always 1 for no reason
    Int_t gainSetup      = 1;  // Gain filter setup (0:1). Here it is always 1 for no reason
    Int_t tailSetup      = 1;  // Tail filter setup (0:1). Here it is always 1 for no reason
    Int_t xtSetup        = 0;  // Cross talk filter setup (0:1). Here it is always 0 for no reason
    Int_t nonlinSetup    = 0;  // Nonlinearity filter setup (0:1). Here it is always 0 for no reason
    Int_t bypassSetup    = 0;  // Filter bypass (for raw data) setup (0:1). Here it is always 0 for no reason
    Int_t commonAdditive = 10; // Digital filter common additive (0:63). Here it is always 10 for no reason
    x = (pedSetup<<31) | (gainSetup<<30) | (tailSetup<<29) | (xtSetup<<28) | (nonlinSetup<<27)
	| (bypassSetup<<26) | (commonAdditive<<20) | 1;
    if (nw < maxSize) {
	buf[nw++] = x; 
    }
    else {
	of++;
    } 
}

//_____________________________________________________________________________
void AliTRDrawData::WriteIntermediateWordsV2(UInt_t* buf, Int_t& nw, Int_t& of, const Int_t& maxSize, const Int_t& det, const Int_t& side) {
    
    Int_t        layer = fGeo->GetLayer( det );   // Layer
    Int_t        stack = fGeo->GetStack( det );   // Stack
    Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
    Int_t           rv = fFee->GetRAWversion();
    const Int_t kNTBin = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
	Bool_t tracklet_on = fFee->GetTracklet();
    UInt_t           x = 0;

    // Write end of tracklet marker
    if (nw < maxSize) buf[nw++] = fgkEndOfTrackletMarker; else of++;
   
  	// Half Chamber header
  	// h[0] (there are 2 HC headers) xmmm mmmm nnnn nnnq qqss sssp ppcc ci01
	// 	, where	 x : Raw version speacial number (=1)
	// 		     m : Raw version major number (test pattern, ZS, disable tracklet, 0, options)
	// 		     n : Raw version minor number
	// 		     q : number of addtional header words (default = 1)
	// 		     s : SM sector number (ALICE numbering)
	// 		     p : plane(layer) number
	//			 c : chamber(stack) number
	//			 i : side number (0:A, 1:B)
	Int_t majorv = 0;	// The major version number 
    Int_t minorv = 0;	// The minor version number
    Int_t add    = 1;	// The number of additional header words to follow : now 1, previous 2
	Int_t TP	 = 0;	// test pattern (default=0)
	Int_t ZS	 = (rv==3) ? 1 : 0;			// zero suppression
	Int_t DT	 = (tracklet_on) ? 0 : 1; 	// disable tracklet 

	majorv = (TP<<6) | (ZS<<5) | (DT<<4) | 1;	// major version

    x = (1<<31) | (majorv<<24) | (minorv<<17) | (add<<14) | (sect<<9) | (layer<<6) | (stack<<3) | (side<<2) | 1;
    if (nw < maxSize) buf[nw++] = x; else of++;
   
    // h[1]		tttt ttbb bbbb bbbb bbbb bbpp pphh hh01
	// , where  t : number of time bins
	//          b : bunch crossing number
	//          p : pretrigger counter
	//          h : pretrigger phase
    Int_t bcCtr   = 99; // bunch crossing counter. Here it is set to 99 always for no reason
    Int_t ptCtr   = 15; // pretrigger counter. Here it is set to 15 always for no reason
    Int_t ptPhase = 11; // pretrigger phase. Here it is set to 11 always for no reason
    //x = (bcCtr<<16) | (ptCtr<<12) | (ptPhase<<8) | ((kNTBin-1)<<2) | 1; 	// old format
    x = ((kNTBin-1)<<26) | (bcCtr<<10) | (ptCtr<<6) | (ptPhase<<2) | 1;
    if (nw < maxSize) buf[nw++] = x; else of++;
  
}

//_____________________________________________________________________________
AliTRDdigitsManager *AliTRDrawData::Raw2DigitsOLD(AliRawReader *rawReader)
{
  //
  // Vx of the raw data reading
  //

  AliTRDarrayADC *digits = 0;
  AliTRDarrayDictionary *track0 = 0;
  AliTRDarrayDictionary *track1 = 0;
  AliTRDarrayDictionary *track2 = 0; 

  AliTRDSignalIndex *indexes = 0;
  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->CreateArrays();

  //AliTRDrawOldStream input(rawReader);
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

 	  if (digits) digits->Compress();
 	  if (track0) track0->Compress();       
 	  if (track1) track1->Compress();       
 	  if (track2) track2->Compress();
	
	  // Add a container for the digits of this detector
	  digits = (AliTRDarrayADC *) digitsManager->GetDigits(det);
	  track0 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,0);
	  track1 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,1);
	  track2 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,2);

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
		  digits->SetData(input.GetRow(), input.GetCol(),input.GetTimeBin() + it, input.GetSignals()[it]);

		  indexes->AddIndexTBin(input.GetRow(), input.GetCol(),
					input.GetTimeBin() + it);
		  track0->SetData(input.GetRow(), input.GetCol(), input.GetTimeBin() + it, 0);
		  track1->SetData(input.GetRow(), input.GetCol(), input.GetTimeBin() + it, 0);
		  track2->SetData(input.GetRow(), input.GetCol(), input.GetTimeBin() + it, 0);
		}
	    }
	}
  }

  if (digits) digits->Compress();
  if (track0) track0->Compress();        
  if (track1) track1->Compress();       
  if (track2) track2->Compress();

  return digitsManager;

}
