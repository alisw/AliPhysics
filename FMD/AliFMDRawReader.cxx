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

//____________________________________________________________________
//
// Class to read ADC values from a AliRawReader object. 
//
// This class uses the AliFMDRawStreamer class to read the ALTRO
// formatted data. 
// 
//          +-------+
//          | TTask |
//          +-------+
//              ^
//              |
//      +-----------------+  <<references>>  +--------------+
//      | AliFMDRawReader |<>----------------| AliRawReader |
//      +-----------------+                  +--------------+
//              |                                  ^
//              | <<uses>>                         |
//              V                                  |
//      +-----------------+      <<uses>>          |
//      | AliFMDRawStream |------------------------+
//      +-----------------+
//              |
//              V
//      +----------------+
//      | AliAltroStream |
//      +----------------+
//
#include <AliLog.h>		// ALILOG_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDRawStream.h"	// ALIFMDRAWSTREAM_H 
#include "AliRawReader.h"	// ALIRAWREADER_H 
#include "AliFMDRawReader.h"	// ALIFMDRAWREADER_H 
#include <TArrayI.h>		// ROOT_TArrayI
#include <TTree.h>		// ROOT_TTree
#include <TClonesArray.h>	// ROOT_TClonesArray

//____________________________________________________________________
ClassImp(AliFMDRawReader)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRawReader::AliFMDRawReader(AliRawReader* reader, TTree* tree) 
  : TTask("FMDRawReader", "Reader of Raw ADC values from the FMD"),
    fTree(tree),
    fReader(reader), 
    fSampleRate(1)
{
  // Default CTOR
  AliFMDParameters* pars = AliFMDParameters::Instance();
  fSampleRate = pars->GetSampleRate();
}


//____________________________________________________________________
void
AliFMDRawReader::Exec(Option_t*) 
{
  // Read raw data into the digits array
  if (!fReader->ReadHeader()) {
    Error("ReadAdcs", "Couldn't read header");
    return;
  }

  Int_t n = 0;
  TClonesArray* array = new TClonesArray("AliFMDDigit");
  fTree->Branch("FMD", &array);

  // Use AliAltroRawStream to read the ALTRO format.  No need to
  // reinvent the wheel :-) 
  AliFMDRawStream input(fReader, fSampleRate);
  // Select FMD DDL's 
  fReader->Select(AliFMDParameters::kBaseDDL >> 8);
  
  Int_t    oldDDL      = -1;
  Int_t    count       = 0;
  UShort_t detector    = 1; // Must be one here
  UShort_t oldDetector = 0;
  Bool_t   next        = kTRUE;

  // local Cache 
  TArrayI counts(10);
  counts.Reset(-1);
  
  // Loop over data in file 
  while (next) {
    next = input.Next();

    count++; 
    Int_t ddl = fReader->GetDDLID();
    if (ddl != oldDDL || input.IsNewStrip() || !next) {
      // Make a new digit, if we have some data (oldDetector == 0,
      // means that we haven't really read anything yet - that is,
      // it's the first time we get here). 
      if (oldDetector > 0) {
	// Got a new strip. 
	AliDebug(10, Form("Add a new strip: FMD%d%c[%2d,%3d] "
			  "(current: FMD%d%c[%2d,%3d])", 
			  oldDetector, input.PrevRing(), 
			  input.PrevSector() , input.PrevStrip(),
			  detector , input.Ring(), input.Sector(), 
			  input.Strip()));
	new ((*array)[n]) AliFMDDigit(oldDetector, 
				      input.PrevRing(), 
				      input.PrevSector(), 
				      input.PrevStrip(), 
				      counts[0], counts[1], counts[2]);
	n++;
#if 0
	AliFMDDigit* digit = 
	  static_cast<AliFMDDigit*>(fFMD->Digits()->
				    UncheckedAt(fFMD->GetNdigits()-1));
#endif 
      }
	
      if (!next) { 
	AliDebug(10, Form("Read %d channels for FMD%d", 
			  count + 1, detector));
	break;
      }
    
    
      // If we got a new DDL, it means we have a new detector. 
      if (ddl != oldDDL) {
	if (detector != 0) 
	  AliDebug(10, Form("Read %d channels for FMD%d", 
			    count + 1, detector));
	// Reset counts, and update the DDL cache 
	count       = 0;
	oldDDL      = ddl;
	// Check that we're processing a FMD detector 
	Int_t detId = fReader->GetDetectorID();
	if (detId != (AliFMDParameters::kBaseDDL >> 8)) {
	  Error("ReadAdcs", "Detector ID %d != %d",
		detId, (AliFMDParameters::kBaseDDL >> 8));
	  break;
	}
	// Figure out what detector we're deling with 
	oldDetector = detector;
	switch (ddl) {
	case 0: detector = 1; break;
	case 1: detector = 2; break;
	case 2: detector = 3; break;
	default:
	  Error("ReadAdcs", "Unknown DDL 0x%x for FMD", ddl);
	  return;
	}
	AliDebug(10, Form("Reading ADCs for 0x%x  - That is FMD%d",
			  fReader->GetEquipmentId(), detector));
      }
      counts.Reset(-1);
    }
    
    counts[input.Sample()] = input.Count();
    
    AliDebug(10, Form("ADC of FMD%d%c[%2d,%3d] += %d",
		      detector, input.Ring(), input.Sector(), 
		      input.Strip(), input.Count()));
    oldDetector = detector;
  }
  fTree->Fill();
  return;

}

//____________________________________________________________________
// 
// EOF
//
