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

// $Id: $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliHLTPluginBase.h"
#include "AliHLTConfiguration.h"
#include "AliReconstruction.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

/**
 * @file HLToutputTodHLTRootObjects.C
 * @brief Macro for converting dHLT output in HLTOUT into ROOT objects or dumping to file.
 *
 * This macro converts dHLT output data blocks in HLTOUT to dHLT ROOT
 * objects which can be used for dHLT specific analysis or debugging.
 * Alternatively the dHLT data blocks are dumped to file in raw binary format.
 *
 * @note The macro must be run in the directory with the raw0, raw1 etc
 * directories. Also no reconstructed ROOT files should be in this directory
 * either. Backup any reconstructed AliESDs.root files if you want, and then
 * delete all root files in the directory where you will run this macro.
 *
 * @param dataSource  Indicates where to find the raw data. If this is a path
 *      then the raw data is expected to be in DDL file format stored in directories
 *      called raw0, raw1 and so on.
 *      If a file name is given instead then it is assumed to be a DATE file unless
 *      it ends in ".root", in which case it is assumed to be rootified raw data.
 * @param dumpBinary  Indicates if the data should be dumped in raw binary format
 *      to files with the file name prefix "output-". If set to false then the
 *      output is rootified and written to a ROOT file "output.root".
 *
 * @author Artur Szostak <artursz@iafrica.com>
 * @ingroup alihlt_dimuon_macros
 */
void HLToutputTodHLTRootObjects(const char* dataSource = "./", bool dumpBinary = false)
{
	// setup of the HLT system
	gSystem->Load("libHLTrec.so");
	AliHLTSystem* sys = AliHLTPluginBase::GetInstance();
	if (sys == NULL)
	{
		cerr << "FATAL ERROR: Cannot get HLT system instance." << endl;
		return;
	}

	// Configure the chain.
	AliHLTConfiguration pub("Pub", "AliHLTOUTPublisher" , NULL, "-origin 'MUON'");
	TString sources = "Pub";
	if (dumpBinary)
	{
		AliHLTConfiguration sink("dhlt-output-dump", "FileWriter", sources, "-datafile output.dat -specfmt");
	}
	else
	{
		AliHLTConfiguration convert("convert", "MUONRootifier", sources, "");
		AliHLTConfiguration sink("dhlt-output-dump", "ROOTFileWriter", "convert", "-concatenate-events -datafile output.root -specfmt");
	}
	
	// Setup the reconstruction and run.
	AliReconstruction rec;
	rec.SetInput(dataSource);
	rec.SetLoadAlignData("");
	rec.SetRunLocalReconstruction("HLT");
	rec.SetRunTracking("");
	rec.SetFillESD("");
	rec.SetRunQA(":");
	rec.SetFillTriggerESD(kFALSE);
	rec.SetRunVertexFinder(kFALSE);
	rec.SetOption("HLT", "libAliHLTMUON.so loglevel=0x78 chains=dhlt-output-dump");
	rec.Run();
}
