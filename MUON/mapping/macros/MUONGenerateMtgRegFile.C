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
//
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include "TString.h"
#include "AliMpDDLStore.h"
#include "AliMpDDL.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpCDB.h"

#endif


// Macro to generate regionnal file for Mtg package (online) from AliRoot trigger mapping
// Ch. Finck (20.08.07)

void MUONGenerateMtgRegFile(Int_t mode = 2, Int_t coinc = 0, UInt_t defaultMask = 0xFFFF, 
			    Int_t version = 20, TString name =  "MtgRegionalCrate" )
{
    // mode = 0,1; 2; 3 -> writing mode no filter; mask; overwrite
    // coinc = 0,1 -> coincidence 3/4, 4/4

    UInt_t mask;
    
    TString fileName = name + Form("-%d", version) + ".dat";
    // Load DDL store
    if ( !AliMpCDB::LoadDDLStore() ) {
      printf("Could not access DDL Store from OCDB !\n");
    }

    FILE* fp = fopen(fileName.Data(), "w");
    // instanciate the elec. mapping class
    AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();

    // loop over the trigger DDL (Right: 20, Left: 21)
    for (Int_t iDDL = 20; iDDL <= 21; ++iDDL) {

      // get ddl object
      AliMpDDL* ddl = ddlStore->GetDDL(iDDL);

      Int_t nCrate = ddl->GetNofTriggerCrates();
    
      // loop over the number of crate in DDL
      for (Int_t index = 0; index < nCrate; ++index) {

	// get crate object
	AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(iDDL, index);
	fprintf(fp, "%s\n", crate->GetName());
	fprintf(fp, "%02x\n", index + (iDDL-20)*8);
	fprintf(fp, "%d\n", mode);
	fprintf(fp, "%d\n", coinc);

	if (index == 7)
	    mask = 0x1ff;
	else
	    mask = defaultMask;
	fprintf(fp, "%04x\n", mask);

	Int_t nLocal = crate->GetNofLocalBoards();

	for (Int_t iLocal = 0; iLocal  < nLocal; ++iLocal) {

	  // get local board Id from crate object
	  Int_t localId = crate->GetLocalBoardId(iLocal);

	  // get local board object
	  AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(localId);

	  // print out some info from local board
          fprintf(fp, "%02d %s %03d %03x\n", localBoard->GetSlot(), localBoard->GetName(), 
                  localId, localBoard->GetSwitch());
          
          // print DE id
          fprintf(fp, " ");
          for (Int_t i = 0; i < localBoard->GetNofDEs(); ++i)
            fprintf(fp, "%4d ", localBoard->GetDEId(i));
          fprintf(fp, "\n");
          
          // print copy card numbers
          fprintf(fp, " %4d %4d", localBoard->GetInputXfrom(), localBoard->GetInputXto());
          fprintf(fp, " %4d %4d", localBoard->GetInputYfrom(), localBoard->GetInputYto());
          fprintf(fp, " %4d\n",   localBoard->GetTC());
	}
      }
    }
    fclose(fp);
}

