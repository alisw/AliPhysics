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

//
// Macro for reading tracker raw data
// Ch. Finck, Subatech Febuary
//

#if !defined(__CINT__) || defined(__MAKECINT__)

// RAW includes
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"

// MUON includes
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCrate.h"
#include "AliMUONLocalTriggerBoard.h"

#endif

// Macro to read rawdata for trigger
// Ch. Finck, Subatech, April. 2006
//
// This macro is interface with AliRawReader for RAW
// The different stucture of the patload are readout and stored in TClonesArray
// with AliMUONRawStreamTrigger classe.
// The macro just simpy read again the TClonesArray contents.
// The parameter of each structure could be seen in the container classes
// AliMUONDarcHeader, AliMUONRegHeader, AliMUONLocalStruct.
// The class AliMUONDDLTrigger manages the structure containers.
// The number of structures in the rawdata file could be set.

void MUONRawStreamTrigger(Int_t maxEvent = 1, Int_t minDDL = 0, Int_t maxDDL = 1, TString fileName = "./")
{
   AliRawReader* rawReader = 0x0;

   if (fileName.EndsWith("/")) {
     rawReader = new AliRawReaderFile(fileName);// DDL files
   } else if (fileName.EndsWith(".root")) {
     rawReader = new AliRawReaderRoot(fileName);
   } else if (!fileName.IsNull()) {
     rawReader = new AliRawReaderDate(fileName);// DATE file
   }

   // raw stream
   AliMUONRawStreamTrigger* rawStream   = new AliMUONRawStreamTrigger(rawReader);

   // set the number of DDL reg & local that are PRESENT in the rawdata file
   // it's NOT the number to be read.
   // default wise set to 2, 8, 16 respectively.
   //    rawStream->SetMaxDDL(xx);
   //    rawStream->SetMaxReg(xx);
   //    rawStream->SetMaxLoc(xx);

   // containers
   AliMUONDDLTrigger*       ddlTrigger  = 0x0;
   AliMUONDarcHeader*       darcHeader  = 0x0;
   AliMUONRegHeader*        regHeader   = 0x0;
   AliMUONLocalStruct*      localStruct = 0x0;

   // crate manager
   AliMUONTriggerCrateStore* crateManager = new AliMUONTriggerCrateStore();   
   crateManager->ReadFromFile();

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

	 printf("RegionalId %d\n", iReg);

	 regHeader =  darcHeader->GetRegHeaderEntry(iReg);
	 //  printf("Reg length %d\n",regHeader->GetHeaderLength());

	 // crate info
	 AliMUONTriggerCrate* crate = crateManager->Crate(rawStream->GetDDL(), iReg);
	 TObjArray *boards = crate->Boards();

	 // loop over local structures
	 Int_t nLocal = regHeader->GetLocalEntries();
	 for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) {  

	   localStruct = regHeader->GetLocalEntry(iLocal);

	   // check if trigger 
	   if (localStruct->GetTriggerY() == 0) { // no empty data

	       // local trigger circuit number
	       AliMUONLocalTriggerBoard* localBoard = (AliMUONLocalTriggerBoard*)boards->At(iLocal+1);

	       printf("LocalId %d\n", localStruct->GetId());

	       Int_t iLocCard = localBoard->GetNumber();
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
}
