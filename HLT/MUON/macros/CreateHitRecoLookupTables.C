/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
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
 * \file CreateHitRecoLookupTables.C
 * \brief Macro used to generate lookup tables for the hit reconstructor components.
 *
 * This macro is used to generate the lookup tables for the hit reconstructor
 * component. All alignment and geometry data is taken from the CDB.
 *
 * \note The LUT files must be generated on the same platform / machine on which
 * they will be used, since they may not be binary compatible across platforms.
 *
 * To run this macro copy "rootlogon.C" from $ALICE_ROOT/HLT/MUON/macros
 * into the current directory, then from the shell command prompt run one of
 * the following commands:
 * \code
 *   > aliroot $ALICE_ROOT/HLT/MUON/macros/CreateHitRecoLookupTables.C
 * \endcode
 * or
 * \code
 *   > aliroot -b -q -l $ALICE_ROOT/HLT/MUON/macros/CreateHitRecoLookupTables.C+
 * \endcode
 *
 * \author Indranil Das <indra.das@saha.ac.in>
 */

/*
Purpose:  A macro to generate LookupTable
          in the following form
buspatchId+manuid+channelId  buspatchId Ix  IY  X  Y  B/NB

Created:  7/10/2005
Modified: 22/12/2005
Modified: 09/02/2006
Modified: 09/04/2007
Modified: 24/08/2007 (To adopt to AliRoot v4-06-Release)

Run Info: To run this macro copy "rootlogon.C" from $ALICE_ROOT/HLT/MUON/macros
  into the current directory then compile and run from inside AliRoot using
      root [0] .x CreateHitRecoLookupTables.C+

Author:   Indranil Das, HEP, SINP, Kolkata
Email:    indra.das@saha.ac.in

18 Apr 2008: Moved lookup table generation code to the AliHLTMUONHitReconstructorComponent
component since it is required there anyway and we want to avoid code duplication.
This also makes this macro much cleaner.
    -- Artur Szostak <artursz@iafrica.com>
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "AliHLTMUONHitReconstructorComponent.h"
#include "TSystem.h"
#endif

/**
 * Generates the ASCII lookup tables for the AliHLTMUONHitReconstructorComponent
 * components. The tables are generated from the CDB database information.
 * \param CDBPath  This is the CDB path to use as the DB storage.
 *                 (Default = local://$ALICE_ROOT/OCDB)
 * \param run  This is the run number to use for the CDB (Default = 0).
 */
void CreateHitRecoLookupTables(const char* CDBPath = "local://$ALICE_ROOT/OCDB", Int_t run = 0)
{
	gSystem->Load("libAliHLTMUON.so");

	for (Int_t ddl = 12; ddl < 20; ddl++)
	{
		Char_t filename[64];
		sprintf(filename, "Lut%d.dat", ddl+1);
		cout << "Generating LUT for DDL " << ddl+1
			<< " and writing output to file " << filename << endl;
		bool ok = AliHLTMUONHitReconstructorComponent::GenerateLookupTable(
				ddl, filename, CDBPath, run
			);
		if (! ok) return;
	}
	
	cout << "Lookup tables have been generated successfully." << endl;
}
