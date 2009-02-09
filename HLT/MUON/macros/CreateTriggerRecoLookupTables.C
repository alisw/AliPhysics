/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
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

/**
 * \ingroup macros
 * \file CreateTriggerRecoLookupTables.C
 * \brief Macro used to generate lookup tables for the trigger reconstructor components.
 *
 * This macro is used to generate the lookup tables for the trigger reconstructor
 * component. All alignment and geometry data is taken from the CDB.
 *
 * \note The LUT files must be generated on the same platform / machine on which
 * they will be used, since they may not be binary compatible across platforms.
 *
 * To run this macro copy "rootlogon.C" from $ALICE_ROOT/HLT/MUON/macros
 * into the current directory, then from the shell command prompt run one of
 * the following commands:
 * \code
 *   > aliroot $ALICE_ROOT/HLT/MUON/macros/CreateTriggerRecoLookupTables.C
 * \endcode
 * or
 * \code
 *   > aliroot -b -q -l $ALICE_ROOT/HLT/MUON/macros/CreateTriggerRecoLookupTables.C+
 * \endcode
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "TSystem.h"
#endif

/**
 * Generates the lookup tables for the AliHLTMUONTriggerReconstructorComponent
 * components. The tables are generated from the CDB database information.
 * \param CDBPath  This is the CDB path to use as the DB storage.
 *                 (Default = local://$ALICE_ROOT/OCDB)
 * \param run  This is the run number to use for the CDB (Default = 0).
 * \param useCrateId  Indicates if the crate ID should be used for the lookup table
 *            indexing rather than just a sequencial number (Default = true).
 */
void CreateTriggerRecoLookupTables(
		const char* CDBPath = "local://$ALICE_ROOT/OCDB",
		Int_t run = 0,
		bool useCrateId = true
	)
{
	gSystem->Load("libAliHLTMUON.so");

	for (Int_t ddl = 20; ddl < 22; ddl++)
	{
		Char_t filename[64];
		sprintf(filename, "Lut%d.dat", ddl+1);
		cout << "Generating LUT for DDL " << ddl+1
			<< " and writing output to file " << filename << endl;
		bool ok = AliHLTMUONTriggerReconstructorComponent::GenerateLookupTable(
				ddl, filename, CDBPath, run, useCrateId
			);
		if (! ok) return;
	}
	
	cout << "Lookup tables have been generated successfully." << endl;
}
