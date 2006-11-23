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
  AliCaloRawStream::AliCaloRawStream(AliRawReader* rawReader, TString calo) :
  AliAltroRawStream(rawReader),
  fModule(-1),
  fPrevModule(-1),
  fRow(-1),
  fPrevRow(-1),
  fColumn(-1),
  fPrevColumn(-1),
  fGain(0)
{
// create an object to read PHOS/EMCAL raw digits

  SelectRawData(calo);

  // PHOS and EMCAL have differen number of RCU per module
  fNRCU = 4;
  if(calo == "EMCAL")  fNRCU = 2;

  TString path = gSystem->Getenv("ALICE_ROOT/");
  path += calo+"/mapping/RCU";
  TString path2;
  for(Int_t i = 0; i < fNRCU; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    fMapping[i] = new AliCaloAltroMapping(path2.Data());
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
  fGain(0),
  fNRCU(0)
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

  for(Int_t i = 0; i < fNRCU; i++) delete fMapping[i];
}

//_____________________________________________________________________________
void AliCaloRawStream::Reset()
{
  // reset phos/emcal raw stream params
  AliAltroRawStream::Reset();
  fModule = fPrevModule = fRow = fPrevRow = fColumn = fPrevColumn = -1;
  fGain = 0;
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
    if (IsNewHWAddress())
      ApplyAltroMapping();
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

  Short_t hwAddress = GetHWAddress();
  fRow = fMapping[rcuIndex]->GetPadRow(hwAddress);
  fColumn = fMapping[rcuIndex]->GetPad(hwAddress);
  fGain = fMapping[rcuIndex]->GetSector(hwAddress);

}
