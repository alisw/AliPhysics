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

// *
// *
// *
// * this class defines the TOF object to be stored
// * in Reference data a run-by-run basis in order to have
// * info about readout electronics
// * 
// *
// *
// *

#include "AliTOFReadoutInfo.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliTOFReadoutInfo)

//_________________________________________________________

AliTOFReadoutInfo::AliTOFReadoutInfo() :
  TObject(),
  fChainEfficiency(NULL),
  fTRMData(NULL),
  fTRMEmptyEvent(NULL),
  fTRMBadEventCounter(NULL),
  fTRMBadCRC(NULL),
  fChainData(NULL),
  fChainBadStatus(NULL),
  fChainBadEventCounter(NULL),
  fTDCError(NULL),
  fTDCErrorFlags(NULL)
{
  /*
   * default constructor
   */
}

//_________________________________________________________

AliTOFReadoutInfo::~AliTOFReadoutInfo()
{
  /*
   * default destructor
   */

}

//_________________________________________________________

AliTOFReadoutInfo::AliTOFReadoutInfo(const AliTOFReadoutInfo &source) :
  TObject(source),
  fChainEfficiency(source.fChainEfficiency),
  fTRMData(source.fTRMData),
  fTRMEmptyEvent(source.fTRMEmptyEvent),
  fTRMBadEventCounter(source.fTRMBadEventCounter),
  fTRMBadCRC(source.fTRMBadCRC),
  fChainData(source.fChainData),
  fChainBadStatus(source.fChainBadStatus),
  fChainBadEventCounter(source.fChainBadEventCounter),
  fTDCError(source.fTDCError),
  fTDCErrorFlags(source.fTDCErrorFlags)
{
  /*
   * copy constructor
   */

}

//_________________________________________________________

AliTOFReadoutInfo &
AliTOFReadoutInfo::operator=(const AliTOFReadoutInfo &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  
  fChainEfficiency = source.fChainEfficiency;
  fTRMData = source.fTRMData;
  fTRMEmptyEvent = source.fTRMEmptyEvent;
  fTRMBadEventCounter = source.fTRMBadEventCounter;
  fTRMBadCRC = source.fTRMBadCRC;
  fChainData = source.fChainData;
  fChainBadStatus = source.fChainBadStatus;
  fChainBadEventCounter = source.fChainBadEventCounter;
  fTDCError = source.fTDCError;
  fTDCErrorFlags = source.fTDCErrorFlags;

  return *this;
}

