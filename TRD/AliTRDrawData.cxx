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
#include "AliTreeLoader.h"

#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDrawStreamBase.h"
#include "AliTRDcalibDB.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDfeeParam.h"
#include "AliTRDmcmSim.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDdigitsParam.h"

ClassImp(AliTRDrawData)

Int_t AliTRDrawData::fgDataSuppressionLevel = 0;

//_____________________________________________________________________________
AliTRDrawData::AliTRDrawData()
  :TObject()
  ,fRunLoader(NULL)
  ,fGeo(NULL)
  ,fFee(NULL)
  ,fNumberOfDDLs(0)
  ,fTrackletTree(NULL)
  ,fTrackletContainer(NULL)
  ,fSMindexPos(0)
  ,fStackindexPos(0)
  ,fEventCounter(0)
  ,fDigitsParam(NULL)
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
  ,fRunLoader(NULL)
  ,fGeo(NULL)
  ,fFee(NULL)
  ,fNumberOfDDLs(0)
  ,fTrackletTree(NULL)
  ,fTrackletContainer(NULL)
  ,fSMindexPos(0)
  ,fStackindexPos(0)
  ,fEventCounter(0)
  ,fDigitsParam(NULL)
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

  if (fTrackletContainer){
    delete fTrackletContainer;
    fTrackletContainer = NULL;
  }

}

//_____________________________________________________________________________
Bool_t AliTRDrawData::Digits2Raw(TTree *digitsTree, const TTree *tracks )
{
  //
  // Initialize necessary parameters and call one
  // of the raw data simulator selected by SetRawVersion.
  //
  // Currently tracklet output is not spported yet and it
  // will be supported in higher version simulator.
  //

  AliTRDdigitsManager* const digitsManager = new AliTRDdigitsManager();

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
  Bool_t newSM    = kFALSE;  // new SM flag, for writing SM index words
  Bool_t newStack = kFALSE;  // new stack flag, for writing stack index words

  // Get digits parameter
  fDigitsParam = digitsManager->GetDigitsParam();

  // sect is same as iDDL, so I use only sect here.
  for (Int_t sect = 0; sect < fGeo->Nsector(); sect++) { 

    char name[1024];
    sprintf(name,"TRD_%d.ddl",sect + AliTRDrawStreamBase::kDDLOffset);

    AliFstream* of = new AliFstream(name);

    // Write a dummy data header
    AliRawDataHeaderSim  header;  // the event header
    UInt_t hpos = of->Tellp();
    of->WriteBuffer((char *) (& header), sizeof(header));
    
    // Reset payload byte size (payload does not include header).
    Int_t npayloadbyte = 0;

    // check the existance of the data
    // SM index word and Stack index word
    UInt_t *iwbuffer = new UInt_t[109]; // index word buffer; max 109 = 2 SM headers + 67 dummy headers + 5*8 stack headers
    Int_t nheader = 0;
    UInt_t bStackMask = 0x0;
    Bool_t bStackHasData = kFALSE;
    Bool_t bSMHasData = kFALSE;
    
    //iwbuffer[nheader++] = 0x0001a020;   // SM index words 
    iwbuffer[nheader++] = 0x0044a020;   // SM index words | additional SM header:48 = 1 SM header + 47 dummy words(for future use)
    iwbuffer[nheader++] = 0x10404071;   // SM header
    for ( Int_t i=0; i<66; i++ ) iwbuffer[nheader++] = 0x00000000;  // dummy words 
    iwbuffer[nheader++] = 0x10000000;   // end of dummy words
    
    for ( Int_t stack= 0; stack < fGeo->Nstack(); stack++) {
      UInt_t linkMask = 0x0;
      for( Int_t layer = 0; layer < fGeo->Nlayer(); layer++) {
	Int_t iDet = fGeo->GetDetector(layer,stack,sect);
	AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(iDet);
	if ( fgDataSuppressionLevel==0 || digits->HasData() ) {
	  bStackMask = bStackMask | ( 1 << stack ); // active stack mask for new stack
	  linkMask = linkMask | ( 3 << (2*layer) );    // 3 = 0011
	  bStackHasData = kTRUE;
	  bSMHasData = kTRUE;
	} // has data
      } // loop over layer
      
      if ( fgDataSuppressionLevel==0 || bStackHasData ){
	iwbuffer[nheader++] = 0x0007a000 | linkMask;    // stack index word + link masks
	if (fgDataSuppressionLevel==0) iwbuffer[nheader-1] = 0x0007afff;  // no suppression
	iwbuffer[nheader++] = 0x04045b01;               // stack header
	for (Int_t i=0;i<6;i++) iwbuffer[nheader++] = 0x00000000; // 6 dummy words
	bStackHasData = kFALSE;
      }
    } // loop over stack
    
    if ( fgDataSuppressionLevel==0 || bSMHasData ){
      iwbuffer[0] = iwbuffer[0] | bStackMask;  // add stack masks to SM index word
      if (fgDataSuppressionLevel==0) iwbuffer[0] = 0x0044a03f;    // no suppression : all stacks are active
      of->WriteBuffer((char *) iwbuffer, nheader*4);
      AliDebug(11, Form("SM %d index word: %08x", iwbuffer[0]));
      AliDebug(11, Form("SM %d header: %08x", iwbuffer[1]));
    }
    // end of SM & stack header ------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------
    
    // Prepare chamber data
    for( Int_t stack = 0; stack < fGeo->Nstack(); stack++) {
      for( Int_t layer = 0; layer < fGeo->Nlayer(); layer++) {
	
	Int_t iDet = fGeo->GetDetector(layer,stack,sect);
	if (iDet == 0){
	  newEvent = kTRUE; // it is expected that each event has at least one tracklet; 
	  // this is only needed for correct readout tree
	  fEventCounter++;
	  AliDebug(11, Form("New event!! Event counter: %d",fEventCounter));
	}
	
	if ( stack==0 && layer==0 ) newSM = kTRUE;  // new SM flag
	if ( layer==0 ) newStack = kTRUE;           // new stack flag
	AliDebug(15, Form("stack : %d, layer : %d, iDec : %d\n",stack,layer,iDet));
	// Get the digits array
	AliTRDarrayADC *digits = (AliTRDarrayADC *) digitsManager->GetDigits(iDet);
	if (fgDataSuppressionLevel==0 || digits->HasData() ) {  // second part is new!! and is for indicating a new event
	  
	  if (digits->HasData()) digits->Expand();
	  
	  Int_t hcwords = 0;
	  
	  // Process A side of the chamber
	  hcwords = ProduceHcData(digits,0,iDet,hcBuffer,kMaxHcWords,newEvent,newSM);
	  if ( newEvent ) newEvent = kFALSE;
	  //AssignLinkMask(hcBuffer, layer);  // active link mask for this layer(2*HC)
	  of->WriteBuffer((char *) hcBuffer, hcwords*4);
	  npayloadbyte += hcwords*4;
	  
	  // Process B side of the chamber
	  hcwords = ProduceHcData(digits,1,iDet,hcBuffer,kMaxHcWords,newEvent,newSM);
	  of->WriteBuffer((char *) hcBuffer, hcwords*4);
	  npayloadbyte += hcwords*4;
	} // has data
	
      } // loop over layer
    } // loop over stack
    
    // Complete header
    header.fSize = UInt_t(of->Tellp()) - hpos;
    header.SetAttribute(0);  // Valid data
    of->Seekp(hpos);         // Rewind to header position
    of->WriteBuffer((char *) (& header), sizeof(header));
    delete of;
  } // loop over sector(SM)
  
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
	
    //buf[nw++] = 0x0001a03f;   // SM index words
    fSMindexPos = nw;       // memorize position of the SM index word for re-allocating stack mask
    buf[nw++] = 0x0001a020; // SM index words
    buf[nw++] = 0x10404071; // SM header

    fStackindexPos = nw;    // memorize position of the stack index word for future adding
	/*  
    for (Int_t istack=0; istack<5; istack++){
        buf[nw++] = 0x0007afff; // stack index words
        buf[nw++] = 0x04045b01; // stack header
        for (Int_t i=0;i<6;i++) buf[nw++] = 0x00000000; // 6 dummy words
    } // loop over 5 stacks
	*/
}

//_____________________________________________________________________________
void AliTRDrawData::AssignStackMask(UInt_t *buf, Int_t nStack){
    //
    // This function re-assign stack mask active(from 0 to 1) in the SM index word
    //   
    buf[fSMindexPos] = buf[fSMindexPos] | ( 1 << nStack );
}

//_____________________________________________________________________________  
Int_t AliTRDrawData::AddStackIndexWords(UInt_t *buf, Int_t /*nStack*/, Int_t nMax){
    // 
    // This function add stack index words and stack header when there is data for the stack
    //
    //   1) stack index words : ssssssss ssssssss vvvv mmmm mmmmmmmm 
    //      - s : size of stack header (number of header, (default = 0x0007)       
    //      - v : header version (default = 0xa)
    //      - m : link mask (default = 0xfff)
    //      - m : link mask (starting value = 0x000)
    //
    //   2) stack header : vvvvvvvv vvvvvvvv bbbbbbbb rrrr rrr c
    //      - v : hardware design revision (default = 0x0404)
    //      - b : physical board ID (default = 0x5b)
    //      - r : reserved for future use (default = 0000 000)
    //      - c : clean checkout flag (default = 1)
    //  
    //   and 6 dummy words(0x00000000)
    //

    Int_t nAddedWords = 0;  // Number of added words
    if ( ShiftWords(buf, fStackindexPos, 8, nMax)== kFALSE ){
        AliError("Adding stack header failed.");
        return 0;
    }

    buf[fStackindexPos++] = 0x0007a000; // stack index words
    buf[fStackindexPos++] = 0x04045b01; // stack header
    for (Int_t i=0;i<6;i++) buf[fStackindexPos++] = 0x00000000; // 6 dummy words 
    nAddedWords += 8;

    return nAddedWords;
}

//_____________________________________________________________________________
void AliTRDrawData::AssignLinkMask(UInt_t *buf, Int_t nLayer){
    //
    // This function re-assign link mask active(from 0 to 1) in the stack index word
    //   
    buf[fStackindexPos-8] = buf[fStackindexPos-8] | ( 3 << (2*nLayer) );    // 3 = 0011 
}

//_____________________________________________________________________________ 
Bool_t AliTRDrawData::ShiftWords(UInt_t *buf, Int_t nStart, Int_t nWords, Int_t nMax){
    //  
    // This function shifts n words
    //
    //if ( nStart+nWords > sizeof(buf)/sizeof(UInt_t) ){
    //  AliError("Words shift failed. No more buffer space.");
    //  return kFALSE;
    //}

    for ( Int_t iw=nMax; iw>nStart-1; iw--){
        buf[iw+nWords] = buf[iw];
    }
    return kTRUE;
}

//_____________________________________________________________________________
Int_t AliTRDrawData::ProduceHcData(AliTRDarrayADC *digits, Int_t side, Int_t det, UInt_t *buf, Int_t maxSize, Bool_t /*newEvent = kFALSE*/, Bool_t /*newSM = kFALSE*/){
	//
        // This function produces the raw data for one HC, i.e. tracklets, tracklet endmarkers, 
        // raw data, raw data endmarkers. 
	// This function can be used for both ZS and NZS data
	//

  	Int_t           nw = 0;                       // Number of written    words
  	Int_t           of = 0;                       // Number of overflowed words
  	Int_t      *tempnw = &nw;                     // Number of written    words for temp. buffer
  	Int_t      *tempof = &of;                     // Number of overflowed words for temp. buffer
  	Int_t        layer = fGeo->GetLayer( det );   // Layer
  	Int_t        stack = fGeo->GetStack( det );   // Stack
  	Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
	const Int_t kCtype = fGeo->GetStack(det) == 2 ? 0 : 1;                       // Chamber type (0:C0, 1:C1)

	Bool_t trackletOn = fFee->GetTracklet();     // tracklet simulation active?

  	AliDebug(1,Form("Producing raw data for sect=%d layer=%d stack=%d side=%d",sect,layer,stack,side));
        
  	AliTRDmcmSim* mcm = new AliTRDmcmSim();

	UInt_t *tempBuffer = buf; // tempBuffer used to write ADC data
	                          // different in case of tracklet writing
	
	if (trackletOn) {
	  tempBuffer = new UInt_t[maxSize];
	  tempnw = new Int_t(0);
	  tempof = new Int_t(0);
	}
	  
	WriteIntermediateWords(tempBuffer,*tempnw,*tempof,maxSize,det,side);

	if (digits->HasData()) {
	  // scanning direction such, that tracklet-words are sorted in ascending z and then in ascending y order
	  // ROB numbering on chamber and MCM numbering on ROB increase with decreasing z and increasing y
	  for (Int_t iRobRow = 0; iRobRow <= (kCtype + 3)-1; iRobRow++ ) {
	    // ROB number should be increasing
	    Int_t iRob = iRobRow * 2 + side;
	    // MCM on one ROB
	    for (Int_t iMcmRB = 0; iMcmRB < fGeo->MCMmax(); iMcmRB++ ) {
	      Int_t iMcm = 16 - 4*(iMcmRB/4 + 1) + (iMcmRB%4);
	      
	      mcm->Init(det, iRob, iMcm);
	      mcm->SetData(digits);     // no filtering done here (already done in digitizer)
	      if (trackletOn) {
		mcm->Tracklet();
		Int_t tempNw = mcm->ProduceTrackletStream(&buf[nw], maxSize - nw);
		if(  tempNw < 0 ) {
		  of += tempNw;
		  nw += maxSize - nw;
		  AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
		} else {
		  nw += tempNw;
		}
	      }
	      mcm->ZSMapping();  // Calculate zero suppression mapping
	      // at the moment it has to be rerun here
	      // Write MCM data to temp. buffer
	      Int_t tempNw = mcm->ProduceRawStream( &tempBuffer[*tempnw], maxSize - *tempnw, fEventCounter );
	      if ( tempNw < 0 ) {
		*tempof += tempNw;
		*tempnw += maxSize - nw;
		AliError(Form("Buffer overflow detected. Please increase the buffer size and recompile."));
	      } else {
		*tempnw += tempNw;
	      }
	    }
	  }
	  
	  delete mcm;

	  // in case of tracklet writing copy temp data to final buffer
	  if (trackletOn) {
	    if (nw + *tempnw < maxSize) {
	      memcpy(&buf[nw], tempBuffer, *tempnw * sizeof(UInt_t));
	      nw += *tempnw;
	    }
	    else {
	      AliError("Buffer overflow detected");
	    }
	    delete [] tempBuffer;
	    delete tempof;
	    delete tempnw;
	  }
	}

  	// Write end of raw data marker
  	if (nw+3 < maxSize) {
          buf[nw++] = fgkEndOfDataMarker;
          buf[nw++] = fgkEndOfDataMarker;
          buf[nw++] = fgkEndOfDataMarker;
          buf[nw++] = fgkEndOfDataMarker;
  	} else {
          of += 4;
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

  rawReader->Select("TRD"); //[mj]

  AliTRDarrayADC *digits = 0;
  AliTRDarrayDictionary *track0 = 0;
  AliTRDarrayDictionary *track1 = 0;
  AliTRDarrayDictionary *track2 = 0;  

  //AliTRDSignalIndex *indexes = 0;
  // Create the digits manager
  AliTRDdigitsManager* digitsManager = new AliTRDdigitsManager();
  digitsManager->CreateArrays();

  if (!fTrackletContainer) {
    const Int_t kTrackletChmb=256;
    fTrackletContainer = new UInt_t *[2];
    fTrackletContainer[0] = new UInt_t[kTrackletChmb];
    fTrackletContainer[1] = new UInt_t[kTrackletChmb];
    memset(fTrackletContainer[0], 0, kTrackletChmb*sizeof(UInt_t)); //jkl
    memset(fTrackletContainer[1], 0, kTrackletChmb*sizeof(UInt_t)); //jkl
  }

  AliTRDrawStreamBase *pinput = AliTRDrawStreamBase::GetRawStream(rawReader);
  AliTRDrawStreamBase &input = *pinput;
  input.SetRawVersion( fFee->GetRAWversion() ); //<= ADDED by MinJung

  AliInfo(Form("Stream version: %s", input.IsA()->GetName()));

  // ----- preparing tracklet output -----
  AliDataLoader *trklLoader = AliRunLoader::Instance()->GetLoader("TRDLoader")->GetDataLoader("tracklets");
  if (!trklLoader) {
    //AliInfo("Could not get the tracklets data loader, adding it now!");
    trklLoader = new AliDataLoader("TRD.Tracklets.root","tracklets", "tracklets");
    AliRunLoader::Instance()->GetLoader("TRDLoader")->AddDataLoader(trklLoader);
  }
  AliTreeLoader *trklTreeLoader = dynamic_cast<AliTreeLoader*> (trklLoader->GetBaseLoader("tracklets-raw"));
  if (!trklTreeLoader) {
    trklTreeLoader = new AliTreeLoader("tracklets-raw", trklLoader);
    trklLoader->AddBaseLoader(trklTreeLoader);
  }

  if (!trklTreeLoader->Tree())
    trklTreeLoader->MakeTree();

  // Loop through the digits
  Int_t det    = 0;

  while (det >= 0)
    {
      det = input.NextChamber(digitsManager,fTrackletContainer);

      if (*(fTrackletContainer[0]) > 0 || *(fTrackletContainer[1]) > 0) WriteTracklets(det);

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

  if (trklTreeLoader)
    trklTreeLoader->WriteData("OVERWRITE");
  if (trklLoader) 
    trklLoader->UnloadAll();

  if (fTrackletContainer){
    delete [] fTrackletContainer[0];
    delete [] fTrackletContainer[1];
    delete [] fTrackletContainer;
    fTrackletContainer = NULL;
  }

  delete pinput;
  pinput = NULL;

  return digitsManager;
}

//_____________________________________________________________________________
void AliTRDrawData::WriteIntermediateWords(UInt_t* buf, Int_t& nw, Int_t& of, const Int_t& maxSize, const Int_t& det, const Int_t& side) {
    // 
    // write tracklet end marker(0x10001000) 
    // and half chamber headers(H[0] and H[1])
    //
    
    Int_t        layer = fGeo->GetLayer( det );   // Layer
    Int_t        stack = fGeo->GetStack( det );   // Stack
    Int_t         sect = fGeo->GetSector( det );  // Sector (=iDDL)
    Int_t           rv = fFee->GetRAWversion();
    const Int_t kNTBin = fDigitsParam->GetNTimeBins(det);
    Bool_t  trackletOn = fFee->GetTracklet();
    UInt_t           x = 0;

    // Write end of tracklet marker
    if (nw < maxSize){
      buf[nw++] = fgkEndOfTrackletMarker;
      buf[nw++] = fgkEndOfTrackletMarker;     // the number of tracklet end marker should be more than 2
    }
    else {
      of++;
    }
   
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
    Int_t tp	 = 0;	// test pattern (default=0)
    Int_t zs	 = (rv==3) ? 1 : 0;			// zero suppression
    Int_t dt	 = (trackletOn) ? 0 : 1; 	// disable tracklet 
    
    majorv = (tp<<6) | (zs<<5) | (dt<<4) | 1;	// major version
    
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
    x = ((kNTBin)<<26) | (bcCtr<<10) | (ptCtr<<6) | (ptPhase<<2) | 1;
    if (nw < maxSize) buf[nw++] = x; else of++;
}

//_____________________________________________________________________________
Bool_t AliTRDrawData::WriteTracklets(Int_t det)
{
  //
  // Write the raw data tracklets into seperate file
  //

  UInt_t **leaves = new UInt_t *[2];
  for (Int_t i=0; i<2 ;i++){
    leaves[i] = new UInt_t[258];
    leaves[i][0] = det; // det
    leaves[i][1] = i;   // side
    memcpy(leaves[i]+2, fTrackletContainer[i], sizeof(UInt_t) * 256);
  }

  if (!fTrackletTree){
    AliDataLoader *dl = fRunLoader->GetLoader("TRDLoader")->GetDataLoader("tracklets");
    dl->MakeTree();
    fTrackletTree = dl->Tree();
  }

  TBranch *trkbranch = fTrackletTree->GetBranch("trkbranch");
  if (!trkbranch) {
    trkbranch = fTrackletTree->Branch("trkbranch",leaves[0],"det/i:side/i:tracklets[256]/i");
  }

  for (Int_t i=0; i<2; i++){
    if (leaves[i][2]>0) {
      trkbranch->SetAddress(leaves[i]);
      fTrackletTree->Fill();
    }
  }

  //  AliDataLoader *dl = fRunLoader->GetLoader("TRDLoader")->GetDataLoader("tracklets"); //jkl: wrong
  //  dl->WriteData("OVERWRITE"); //jkl: wrong
  //dl->Unload();
  delete [] leaves;

  return kTRUE;
}
