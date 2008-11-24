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
///
/// This class provides access to PHOS/EMCAL digits in raw data.
///
/// It loops over all PHOS/EMCAL digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
/// usage: 
/// root > AliRawReaderFile rawReader ; 
/// root > AliCaloRawStream input(&rawReader) ; 
/// root > while (input.Next()) .....
///
///Modification: Class exported from PHOS to be used by EMCAL and PHOS
///November 2006 Gustavo Conesa Balbastre 
///////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TSystem.h>

#include "AliCaloRawStream.h"
#include "AliRawReader.h"
#include "AliCaloAltroMapping.h"

ClassImp(AliCaloRawStream)


//_____________________________________________________________________________
AliCaloRawStream::AliCaloRawStream(AliRawReader* rawReader, TString calo, AliAltroMapping **mapping) :
  AliAltroRawStream(rawReader),
  fModule(-1),
  fPrevModule(-1),
  fRow(-1),
  fPrevRow(-1),
  fColumn(-1),
  fPrevColumn(-1),
  fCaloFlag(0),
  fFilter(0),
  fNRCU(0),
  fNSides(0),
  fCalo(calo),
  fExternalMapping(kFALSE)
{
  // create an object to read PHOS/EMCAL raw digits
  SelectRawData(calo);

  // PHOS and EMCAL have differen number of RCU per module
  //For PHOS
  fNRCU = 4;
  fNSides = 1;
  //For EMCAL
  TString sides[]={"A","C"};
  if(fCalo == "EMCAL")  {
    fNRCU = 2;
    fNSides = 2;
  }

  if (mapping == NULL) {
    TString path = gSystem->Getenv("ALICE_ROOT");
    path += "/"+fCalo+"/mapping/RCU";
    TString path2;
    for(Int_t j = 0; j < fNSides; j++){
      for(Int_t i = 0; i < fNRCU; i++) {
	path2 = path;
	path2 += i;
	if(fCalo == "EMCAL") path2 += sides[j];
	path2 += ".data";
	//printf("AliCaloRawStream::RCU:  %s\n",path2.Data());
	fMapping[j*fNRCU+ i] = new AliCaloAltroMapping(path2.Data());
      }
    }
  }
  else {
    fExternalMapping = kTRUE;
    //printf("AliCaloRawStream::External mapping N: RCU %d, sides %d \n",  fNRCU,fNSides);
    for(Int_t i = 0; i < fNRCU*fNSides; i++)
      fMapping[i] = mapping[i];
    
  }
  
  SetNoAltroMapping(kFALSE);
}

//_____________________________________________________________________________
AliCaloRawStream::AliCaloRawStream(const AliCaloRawStream& stream) :
  AliAltroRawStream(stream),
  fModule(-1),
  fPrevModule(-1),
  fRow(-1),
  fPrevRow(-1),
  fColumn(-1),
  fPrevColumn(-1),
  fCaloFlag(0),
  fFilter(0),
  fNRCU(0),
  fNSides(0),
  fCalo(""),
  fExternalMapping(kFALSE)
{  
  Fatal("AliCaloRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliCaloRawStream& AliCaloRawStream::operator = (const AliCaloRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliCaloRawStream::~AliCaloRawStream()
{
// destructor

  if (!fExternalMapping)
    for(Int_t i = 0; i < fNRCU*fNSides; i++)
      delete fMapping[i];
}

//_____________________________________________________________________________
void AliCaloRawStream::Reset()
{
  // reset phos/emcal raw stream params
  AliAltroRawStream::Reset();
  fModule = fPrevModule = fRow = fPrevRow = fColumn = fPrevColumn = -1;
  fFilter = fCaloFlag = 0;
  fCalo="";
}

//_____________________________________________________________________________
Bool_t AliCaloRawStream::Next()
{
  // Read next PHOS/EMCAL signal
  // Apply the PHOS/EMCAL altro mapping to get
  // the module,row and column indeces
  fPrevModule = fModule;
  fPrevRow = fRow;
  fPrevColumn = fColumn;
  if (AliAltroRawStream::Next()) {
    if (IsNewHWAddress()) {
      ApplyAltroMapping();
      if ( fFilter > 0 ) { // some data should be filtered out
	if ( (fFilter & (1<<fCaloFlag)) != 0) {  
	  // this particular data should be filtered out
	  Next(); // go to the next address instead
	}
      }
    }
    return kTRUE;
  }
  else
    return kFALSE;
}

//_____________________________________________________________________________
void AliCaloRawStream::ApplyAltroMapping()
{
  // Take the DDL index, load
  // the corresponding altro mapping
  // object and fill the sector,row and pad indeces
  Int_t ddlNumber = GetDDLNumber();
  fModule = ddlNumber / fNRCU;

  Int_t rcuIndex = ddlNumber % fNRCU;

  if(fCalo=="EMCAL"){ // EMCAL may need to increase RCU index for the maps
    if (fModule%2 == 1) { rcuIndex += 2; } // other='C' side maps
  }

  Short_t hwAddress = GetHWAddress();
  fRow = fMapping[rcuIndex]->GetPadRow(hwAddress);
  fColumn = fMapping[rcuIndex]->GetPad(hwAddress);
  fCaloFlag = fMapping[rcuIndex]->GetSector(hwAddress);

}
