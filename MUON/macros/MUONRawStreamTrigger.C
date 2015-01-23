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

// $Id$

/// \ingroup macros
/// \file MUONRawStreamTrigger.C
/// \brief Macro for reading trigger raw data
///
/// \author Ch. Finck, Subatech, April 2006
///
/// Added example routines to show how to use the interface of the high
/// performance decoder AliMUONRawStreamTriggerrHP.
///    --  Artur Szostak <artursz@iafrica.com>
///
/// Implement "digits" iterator.
/// This macro is interface with AliRawReader for RAW.
/// The different stucture of the patload are readout and stored in TClonesArray
/// with AliMUONRawStreamTrigger class.
/// The macro just simply reads again the TClonesArray contents.
/// The parameter of each structure could be seen in the container classes
/// AliMUONDarcHeader, AliMUONRegHeader, AliMUONLocalStruct.
/// The class AliMUONDDLTrigger manages the structure containers.
/// The number of structures in the rawdata file could be set.
/// The DATE format reading is no more supported please use the MUONTRGda code.


#if !defined(__CINT__) || defined(__MAKECINT__)

// RAW includes
#include "AliRawReader.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"

// MUON includes
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
#include "AliMpTriggerCrate.h"
#include "AliMpDDLStore.h"
#include "AliMpCDB.h"

#include "TStopwatch.h"


#endif

void MUONRawStreamTrigger(Int_t maxEvent = 1, Int_t minDDL = 0, Int_t maxDDL = 1, TString fileName = "./")
{
   
   TStopwatch timer;
   timer.Start(kTRUE);

   AliRawReader* rawReader = AliRawReader::Create(fileName.Data());

   // Load mapping
     if ( ! AliMpCDB::LoadDDLStore() ) {
       printf("Could not access mapping from OCDB !\n");
     }

   // raw stream
   AliMUONRawStreamTrigger* rawStream   = new AliMUONRawStreamTrigger(rawReader);

   // set the number ofreg & local that are PRESENT in the rawdata file
   // it's NOT the number to be read.
   // default wise set to 8, 16 respectively.
   //    rawStream->SetMaxReg(2);
   //    rawStream->SetMaxLoc(xx);

   // containers
   AliMUONDDLTrigger*       ddlTrigger  = 0x0;
   AliMUONDarcHeader*       darcHeader  = 0x0;
   AliMUONRegHeader*        regHeader   = 0x0;
   AliMUONLocalStruct*      localStruct = 0x0;


   // Loop over events  
   Int_t iEvent = 0;

   while (rawReader->NextEvent()) {

     if (iEvent == maxEvent)
       break;

     printf("Event %d\n",iEvent++);

     // read DDL while < 2 DDL
     while(rawStream->NextDDL()) {

      if (rawStream->GetDDL() < minDDL || rawStream->GetDDL() > maxDDL)
	 continue;

       printf("\niDDL %d\n", rawStream->GetDDL());

       ddlTrigger = rawStream->GetDDLTrigger();
       darcHeader = ddlTrigger->GetDarcHeader();

       printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());

       // loop over regional structures
       Int_t nReg = darcHeader->GetRegHeaderEntries();
       for(Int_t iReg = 0; iReg < nReg ;iReg++){   //REG loop

//	 printf("RegionalId %d\n", iReg);

	 regHeader =  darcHeader->GetRegHeaderEntry(iReg);
	 //  printf("Reg length %d\n",regHeader->GetHeaderLength());

	 // crate info  
	 AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
	                            GetTriggerCrate(rawStream->GetDDL(), iReg);

	 // loop over local structures
	 Int_t nLocal = regHeader->GetLocalEntries();
	 for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) {  

	   localStruct = regHeader->GetLocalEntry(iLocal);

	   Int_t iLocCard = crate->GetLocalBoardId(localStruct->GetId());

	   if ( !iLocCard ) continue; // empty slot

	   // check if trigger 
 	   if (localStruct->GetTriggerX() 
	       || localStruct->GetTriggerY()) { // no empty data

	     printf("LocalId %d\n", localStruct->GetId());

	     Int_t loStripX  = (Int_t)localStruct->GetXPos();
	     Int_t loStripY  = (Int_t)localStruct->GetYPos();
	     Int_t loDev     = (Int_t)localStruct->GetXDev();
	       
	     printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", iLocCard, loStripX, loStripY, loDev);

	   }
	 } // iLocal
       } // iReg
     } // NextDDL
   }// NextEvent

   delete rawReader;
   delete rawStream;

   timer.Print();

}


void MUONRawStreamTriggerHPExpert(Int_t maxEvent = 1, Int_t minDDL = 0, Int_t maxDDL = 1, TString fileName = "./")
{
   /// Decodes raw MUON trigger data using the fast decoder.
   
   TStopwatch timer;
   timer.Start(kTRUE);

   AliRawReader* rawReader = AliRawReader::Create(fileName.Data());

   // Load mapping
   if ( ! AliMpCDB::LoadDDLStore() ) {
     printf("Could not access mapping from OCDB !\n");
   }

   // raw stream
   AliMUONRawStreamTriggerHP* rawStream = new AliMUONRawStreamTriggerHP(rawReader);

   // Light weight interfaces to the raw data.
   const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
   const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
   const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

   // Loop over events.
   Int_t iEvent = 0;

   while (rawReader->NextEvent()) {

     if (iEvent == maxEvent)
       break;

     printf("Event %d\n",iEvent++);

     // read DDL while < 2 DDL
     while(rawStream->NextDDL()) {

      if (rawStream->GetDDL() < minDDL || rawStream->GetDDL() > maxDDL)
	 continue;

       printf("\niDDL %d\n", rawStream->GetDDL());

       darcHeader = rawStream->GetHeaders();

       printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());

       // loop over regional structures
       Int_t nReg = (Int_t) rawStream->GetRegionalHeaderCount();
       for(Int_t iReg = 0; iReg < nReg ;iReg++){   //REG loop

//	 printf("RegionalId %d\n", iReg);

	 regHeader =  rawStream->GetRegionalHeader(iReg);

	 // crate info
	 AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
	                            GetTriggerCrate(rawStream->GetDDL(), iReg);

	 // loop over local structures
	 Int_t nLocal = rawStream->GetLocalStructCount(iReg);
	 for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) {

	   localStruct = rawStream->GetLocalStruct(iReg, iLocal);

	   Int_t iLocCard = crate->GetLocalBoardId(localStruct->GetId());

	   if ( !iLocCard ) continue; // empty slot
	 	
	   // check if trigger 
 	   if (localStruct->GetTriggerX() 
	       || localStruct->GetTriggerY()) { // no empty data

	     printf("LocalId %d\n", localStruct->GetId());

	     Int_t loStripX  = (Int_t)localStruct->GetXPos();
	     Int_t loStripY  = (Int_t)localStruct->GetYPos();
	     Int_t loDev     = (Int_t)localStruct->GetXDev();

	     printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", iLocCard, loStripX, loStripY, loDev);

	   }
	 } // iLocal
       } // iReg
     } // NextDDL
   }// NextEvent

   delete rawReader;
   delete rawStream;

   timer.Print();

}


void MUONRawStreamTriggerHPExpert2(Int_t maxEvent = 1, Int_t minDDL = 0, Int_t maxDDL = 1, TString fileName = "./")
{
   /// Decodes raw MUON trigger data using the fast decoder.
   /// This method shows a slightly different method of iteration over the
   /// data using the AliMUONRawStreamTriggerHP interface.
   
   TStopwatch timer;
   timer.Start(kTRUE);

   AliRawReader* rawReader = AliRawReader::Create(fileName.Data());

   // Load mapping
   if ( ! AliMpCDB::LoadDDLStore() ) {
     printf("Could not access mapping from OCDB !\n");
   }

   // raw stream
   AliMUONRawStreamTriggerHP* rawStream = new AliMUONRawStreamTriggerHP(rawReader);

   // Light weight interfaces to the raw data.
   const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
   const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
   const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

   // Loop over events.
   Int_t iEvent = 0;

   while (rawReader->NextEvent()) {

     if (iEvent == maxEvent)
       break;

     printf("Event %d\n",iEvent++);

     // read DDL while < 2 DDL
     while(rawStream->NextDDL()) {

      if (rawStream->GetDDL() < minDDL || rawStream->GetDDL() > maxDDL)
	 continue;

       printf("\niDDL %d\n", rawStream->GetDDL());

       darcHeader = rawStream->GetHeaders();

       printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());

       // loop over regional structures
       Int_t nReg = (Int_t) rawStream->GetRegionalHeaderCount();
       for(Int_t iReg = 0; iReg < nReg ;iReg++){   //REG loop

//	 printf("RegionalId %d\n", iReg);

	 regHeader =  rawStream->GetRegionalHeader(iReg);

	 // crate info
	 AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
	                            GetTriggerCrate(rawStream->GetDDL(), iReg);

	 // loop over local structures
	 Int_t nLocal = regHeader->GetLocalStructCount();
	 for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) {

	   localStruct = regHeader->GetLocalStruct(iLocal);

	   Int_t iLocCard = crate->GetLocalBoardId(localStruct->GetId());

	   if ( !iLocCard ) continue; // empty slot
	 	
	   // check if trigger 
 	   if (localStruct->GetTriggerX() 
	       || localStruct->GetTriggerY()) { // no empty data

	     printf("LocalId %d\n", localStruct->GetId());

	     Int_t loStripX  = (Int_t)localStruct->GetXPos();
	     Int_t loStripY  = (Int_t)localStruct->GetYPos();
	     Int_t loDev     = (Int_t)localStruct->GetXDev();

	     printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", iLocCard, loStripX, loStripY, loDev);

	   }
	 } // iLocal
       } // iReg
     } // NextDDL
   }// NextEvent

   delete rawReader;
   delete rawStream;

   timer.Print();

}


void MUONRawStreamTriggerHPExpert3(Int_t maxEvent = 1, Int_t minDDL = 0, Int_t maxDDL = 1, TString fileName = "./")
{
   /// Decodes raw MUON trigger data using the fast decoder.
   /// This method shows a jet another method of iteration over the
   /// data using the AliMUONRawStreamTriggerHP interface.
   
   TStopwatch timer;
   timer.Start(kTRUE);

   AliRawReader* rawReader = AliRawReader::Create(fileName.Data());

   // Load mapping
   if ( ! AliMpCDB::LoadDDLStore() ) {
     printf("Could not access mapping from OCDB !\n");
   }

   // raw stream
   AliMUONRawStreamTriggerHP* rawStream = new AliMUONRawStreamTriggerHP(rawReader);

   // Light weight interfaces to the raw data.
   const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
   const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
   const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

   // Loop over events.
   Int_t iEvent = 0;

   while (rawReader->NextEvent()) {

     if (iEvent == maxEvent)
       break;

     printf("Event %d\n",iEvent++);

     // read DDL while < 2 DDL
     while(rawStream->NextDDL()) {

      if (rawStream->GetDDL() < minDDL || rawStream->GetDDL() > maxDDL)
	 continue;

       printf("\niDDL %d\n", rawStream->GetDDL());

       darcHeader = rawStream->GetHeaders();

       printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());

       // loop over regional structures
       Int_t iReg = 0;
       regHeader = rawStream->GetFirstRegionalHeader();
       while (regHeader != NULL)
       {
//	 printf("RegionalId %d\n", iReg);

	 regHeader =  rawStream->GetRegionalHeader(iReg);

	 // crate info
	 AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
	                            GetTriggerCrate(rawStream->GetDDL(), iReg);

	 // loop over local structures
	 Int_t iLocal = 0;
	 localStruct = regHeader->GetFirstLocalStruct();
	 while (localStruct != NULL)
	 {
	   Int_t iLocCard = crate->GetLocalBoardId(localStruct->GetId());
	   if ( !iLocCard ) continue; // empty slot
	 	
	   // check if trigger 
 	   if (localStruct->GetTriggerX() 
	       || localStruct->GetTriggerY()) { // no empty data

	     printf("LocalId %d\n", localStruct->GetId());

	     Int_t loStripX  = (Int_t)localStruct->GetXPos();
	     Int_t loStripY  = (Int_t)localStruct->GetYPos();
	     Int_t loDev     = (Int_t)localStruct->GetXDev();

	     printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", iLocCard, loStripX, loStripY, loDev);

	   }
	   iLocal++;
	   localStruct = localStruct->Next();
	 } // iLocal
	 iReg++;
	 regHeader = regHeader->Next();
       } // iReg
     } // NextDDL
   }// NextEvent

   delete rawReader;
   delete rawStream;

   timer.Print();

}


void MUONRawStreamTriggerSimple(Int_t maxEvent = 1, TString fileName = "./")
{
  /// Reads the raw data in fileName, using a simplified interface (iterator
  /// over local structure response).

  TStopwatch timer;
  timer.Start(kTRUE);

   AliRawReader* rawReader = 0x0;

   if (fileName.EndsWith("/")) {
     rawReader = new AliRawReaderFile(fileName);// DDL files
   } else if (fileName.EndsWith(".root")) {
     rawReader = new AliRawReaderRoot(fileName);
   } 

   // raw stream
   AliMUONRawStreamTrigger* rawStream   = new AliMUONRawStreamTrigger(rawReader);

   // set the number of reg & local that are PRESENT in the rawdata file
   // it's NOT the number to be read.
   // default wise set to 8, 16 respectively.
   //    rawStream->SetMaxReg(2);
   //    rawStream->SetMaxLoc(xx);

   UChar_t id;   
   UChar_t dec;
   Bool_t trigY; 
   UChar_t yPos; 
   UChar_t sXDev; 
   UChar_t xDev;
   UChar_t xPos;

   Bool_t triggerX; 
   Bool_t triggerY; 

 
   TArrayS xPattern; 
   TArrayS yPattern;

   // Loop over events  
   Int_t iEvent = 0;

   while (rawReader->NextEvent()) {

     if (iEvent == maxEvent)
	 break;

     printf("Event %d\n",iEvent++);

     rawStream->First();

     // read while there are digits
     while( rawStream->Next(id, dec, trigY, yPos, sXDev, xDev, xPos,
			    triggerX, triggerY, xPattern, yPattern) ) 
     {
       if ( triggerX || triggerY )  // no empty data
	   printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", id, xPos, yPos, xDev);

     }// Next
  
   }// NextEvent

   delete rawReader;
   delete rawStream;

   timer.Print();

}


void MUONRawStreamTriggerHPSimple(Int_t maxEvent = 1, TString fileName = "./")
{
  /// Reads the raw data in fileName, using the simplified interface of the
  /// high performance decoder.

  TStopwatch timer;
  timer.Start(kTRUE);

  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  AliMUONRawStreamTriggerHP* rawStream   = new AliMUONRawStreamTriggerHP(rawReader);

  UChar_t id;
  UChar_t dec;
  Bool_t trigY;
  UChar_t yPos;
  UChar_t sXDev;
  UChar_t xDev;
  UChar_t xPos;
  Bool_t triggerX;
  Bool_t triggerY;
  TArrayS xPattern;
  TArrayS yPattern;

  // Loop over events
  Int_t iEvent = 0;

  while (rawReader->NextEvent())
  {

    if (iEvent == maxEvent)
       break;

     printf("Event %d\n", iEvent++);

     rawStream->First();

     // read while there are digits
     while ( rawStream->Next(id, dec, trigY, yPos, sXDev, xDev, xPos,
                            triggerX, triggerY, xPattern, yPattern) )
     {
       if ( triggerX || triggerY )  // no empty data
          printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", id, xPos, yPos, xDev);
     }// Next
  }// NextEvent

  delete rawReader;
  delete rawStream;

  timer.Print();
}


void MUONRawStreamTriggerHPSimple2(Int_t maxEvent = 1, TString fileName = "./")
{
  /// Reads the raw data in fileName, using the simplified interface of the
  /// high performance decoder.

  TStopwatch timer;
  timer.Start(kTRUE);

  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  AliMUONRawStreamTriggerHP* rawStream   = new AliMUONRawStreamTriggerHP(rawReader);

  const AliMUONRawStreamTriggerHP::AliLocalStruct* localStruct;

  // Loop over events
  Int_t iEvent = 0;

  while (rawReader->NextEvent())
  {

    if (iEvent == maxEvent)
       break;

     printf("Event %d\n", iEvent++);

     rawStream->First();

     // read while there are digits
     while ( (localStruct = rawStream->Next()) != NULL )
     {
       if ( localStruct->GetTriggerX() || localStruct->GetTriggerY() )  // no empty data
          printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n",
          	localStruct->GetId(),
          	localStruct->GetXPos(),
          	localStruct->GetYPos(),
          	localStruct->GetXDev()
          );
     }// Next
  }// NextEvent

  delete rawReader;
  delete rawStream;

  timer.Print();
}

