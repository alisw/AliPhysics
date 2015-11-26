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

/* $Id: $ */

///////////////////////////////////////////////////////////////////////////////
//
// This class provides access to PHOS/EMCAL digits in raw data.
//
// It loops over all PHOS/EMCAL digits in the raw data given by the AliRawReader.
// The Next method goes to the next digit. If there are no digits left
// it returns kFALSE.
// Several getters provide information about the current digit.
// usage: 
//    AliRawReader *reader = AliRawReader::Create(fileName);
//    AliCaloRawStreamV3 *stream = new AliCaloRawStreamV3(reader,calo);
//    while (reader->NextEvent())
//      while (stream->NextDDL())
//        while (stream->NextChannel()) ...
///
/// Yuri Kharlov. 23 June 2009
///////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliCaloRawStreamV3.h"
#include "AliRawReader.h"
#include "AliCaloAltroMapping.h"

ClassImp(AliCaloRawStreamV3)


//_____________________________________________________________________________
AliCaloRawStreamV3::AliCaloRawStreamV3(AliRawReader* rawReader, TString calo, AliAltroMapping **mapping) :
  AliAltroRawStreamV3(rawReader),
  fModule(-1),
  fRow(-1),
  fColumn(-1),
  fCaloFlag(0),
  fNModules(0),
  fNRCU(0),
  fNSides(0),
  fCalo(calo),
  fExternalMapping(kFALSE)
{
  // create an object to read PHOS/EMCAL raw digits
  SelectRawData(calo);

  // PHOS and EMCAL have different number of RCU per module
  //For PHOS (different mappings for different modules)
  if(fCalo == "PHOS")  {
    fNModules = 5;
    fNRCU     = 4;
    fNSides   = 1;
  }
  //For EMCAL (the same mapping for all modules)
  if(fCalo == "EMCAL")  {
    fNModules = 1;
    fNRCU     = 2;
    fNSides   = 2;
  }
  TString sides[]={"A","C"};

  if (mapping == NULL) {
    // Read mapping files from $ALICE_ROOT/CALO/mapping/*.data
    TString path = gSystem->Getenv("ALICE_ROOT");
    path += "/"+fCalo+"/mapping/";
    TString path1, path2;
    for(Int_t m = 0; m < fNModules; m++) {
      path1 = path;
      if     (fCalo == "EMCAL") {
	path1 += "RCU";
      }
      else if(fCalo == "PHOS" ) {
	path1 += "Mod";
	path1 += m;
	path1 += "RCU";
      }
      for(Int_t j = 0; j < fNSides; j++){
	for(Int_t i = 0; i < fNRCU; i++) {
	  path2 = path1;
	  path2 += i;
	  if(fCalo == "EMCAL") path2 += sides[j];
	  path2 += ".data";
	  AliDebug(2,Form("Mapping file: %s",path2.Data()));
	  fMapping[m*fNSides*fNRCU + j*fNRCU + i] = new AliCaloAltroMapping(path2.Data());
	}
      }
    }
  }
  else {
    // Mapping is supplied by reconstruction
    fExternalMapping = kTRUE;
    for(Int_t i = 0; i < fNModules*fNRCU*fNSides; i++)
      fMapping[i] = mapping[i];
  }
}

//_____________________________________________________________________________
AliCaloRawStreamV3::AliCaloRawStreamV3(const AliCaloRawStreamV3& stream) :
  AliAltroRawStreamV3(stream),
  fModule(-1),
  fRow(-1),
  fColumn(-1),
  fCaloFlag(0),
  fNModules(0),
  fNRCU(0),
  fNSides(0),
  fCalo(""),
  fExternalMapping(kFALSE)
{  
  // Dummy copy constructor
  Fatal("AliCaloRawStreamV3", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliCaloRawStreamV3& AliCaloRawStreamV3::operator = (const AliCaloRawStreamV3& 
					      /* stream */)
{
  // Dummy assignment operator
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliCaloRawStreamV3::~AliCaloRawStreamV3()
{
// destructor

  if (!fExternalMapping)
    for(Int_t i = 0; i < fNModules*fNRCU*fNSides; i++)
      delete fMapping[i];
}

//_____________________________________________________________________________
void AliCaloRawStreamV3::Reset()
{
  // reset PHOS/EMCAL raw stream params
  AliAltroRawStreamV3::Reset();
  fModule = fRow = fColumn = -1;
  fCaloFlag = 0;
  fCalo="";
}

//_____________________________________________________________________________
Bool_t AliCaloRawStreamV3::NextChannel()
{
  // Read next PHOS/EMCAL signal
  // Apply the PHOS/EMCAL altro mapping to get
  // the module,row and column indeces

  if (AliAltroRawStreamV3::NextChannel()) {
    ApplyAltroMapping();
    return kTRUE;
  }
  else
    return kFALSE;
}

//_____________________________________________________________________________
void AliCaloRawStreamV3::ApplyAltroMapping()
{
  // Take the DDL index, load
  // the corresponding altro mapping
  // object and fill the sector,row and pad indeces

  Int_t ddlNumber = GetDDLNumber();
  fModule = ddlNumber / fNRCU;

  Int_t rcuIndex = ddlNumber % fNRCU;

  if( fNModules > 1) rcuIndex += fModule*fNRCU*fNSides;
  if( fNRCU == 2 ){ // EMCAL may need to increase RCU index for the maps
    if (fModule%2 == 1) { rcuIndex += 2; } // other='C' side maps
  }

  Short_t hwAddress = GetHWAddress();
  if(rcuIndex > -1 && rcuIndex < 20 && hwAddress > -1) {
    fRow      = fMapping[rcuIndex]->GetPadRow(hwAddress);
    fColumn   = fMapping[rcuIndex]->GetPad(hwAddress);
    fCaloFlag = fMapping[rcuIndex]->GetSector(hwAddress);
  }
}
