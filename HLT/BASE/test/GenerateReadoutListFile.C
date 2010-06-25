// $Id: $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   GenerateReadoutListFile.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   25 June 2010
/// @brief  Macro for generating AliHLTReadoutList readout list objects for testing.
///
/// This macro generates AliHLTReadoutList objects and writes them to a ROOT file
/// directly as an object and also in a TTree to check the behaviour of the streamers.
/// The generated file should be used by testAliHLTEventDDLBackwardCompatibility.C
/// to test the backward compatibility of AliHLTReadoutList objects.
///
/// The macro can be run as follows:
/// \code
///   aliroot -b -q GenerateReadoutListFile.C\(\"myoutputfile.root\"\)
/// \endcode
/// This will generate the file myoutputfile.root in the current directory.
/// By omitting the file name the default file name "output.root" will be used.

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliHLTReadoutList.h"
#include "TFile.h"
#include "TTree.h"
#include "Riostream.h"
#endif

void GenerateReadoutListFile(const char* filename = "output.root")
{
	/// Generates the test file with readout list objects.
	
	TFile* file = new TFile(filename, "RECREATE");
	if (file == NULL)
	{
		cerr << "ERROR: Could not create file: " << filename << endl;
		return;
	}
	
	AliHLTReadoutList* r = new AliHLTReadoutList();
	TTree* tree = new TTree("rltree","Tree containing readout list objects.");
	tree->Branch("readoutlist", "AliHLTReadoutList", &r);
	
	r->Enable(AliHLTReadoutList::kTRG);
	r->Write("readoutlist1", TObject::kOverwrite);
	tree->Fill();
	r->Enable(AliHLTReadoutList::kEMCAL);
	r->Write("readoutlist2", TObject::kOverwrite);
	tree->Fill();
	r->Enable(AliHLTReadoutList::kDAQTEST);
	r->Write("readoutlist3", TObject::kOverwrite);
	tree->Fill();
	r->Enable(AliHLTReadoutList::kHLT);
	r->Write("readoutlist4", TObject::kOverwrite);
	tree->Fill();
	r->Disable(AliHLTReadoutList::kEMCAL);
	r->Write("readoutlist5", TObject::kOverwrite);
	tree->Fill();
	
	tree->Write("rltree", TObject::kOverwrite);
	delete file;
}
