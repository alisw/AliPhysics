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

//-------------------------------------------------------------------------
//                      Implementation of   Class AliESDHeader
//   Header data
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

#include "AliESDHeader.h"


ClassImp(AliESDHeader)

//______________________________________________________________________________
AliESDHeader::AliESDHeader() :
  AliVHeader(),
  fTriggerMask(0),
  fOrbitNumber(0),
  fTimeStamp(0),
  fEventType(0),
  fPeriodNumber(0),
  fEventNumberInFile(0),
  fBunchCrossNumber(0),
  fTriggerCluster(0)
{
  // default constructor

  SetName("AliESDHeader");
}


AliESDHeader::~AliESDHeader() 
{
  // destructor
}


AliESDHeader::AliESDHeader(const AliESDHeader &header) :
  AliVHeader(header),
  fTriggerMask(header.fTriggerMask),
  fOrbitNumber(header.fOrbitNumber),
  fTimeStamp(header.fTimeStamp),
  fEventType(header.fEventType),
  fPeriodNumber(header.fPeriodNumber),
  fEventNumberInFile(header.fEventNumberInFile),
  fBunchCrossNumber(header.fBunchCrossNumber),
  fTriggerCluster(header.fTriggerCluster)
{
  // copy constructor

  SetName(header.fName);
  SetTitle(header.fTitle);
}

AliESDHeader& AliESDHeader::operator=(const AliESDHeader &header)
{ 
  // assigment operator
  if(this!=&header) {
    AliVHeader::operator=(header);
    fTriggerMask = header.fTriggerMask;
    fOrbitNumber = header.fOrbitNumber;
    fTimeStamp = header.fTimeStamp;
    fEventType = header.fEventType;
    fPeriodNumber = header.fPeriodNumber;
    fEventNumberInFile = header.fEventNumberInFile;
    fBunchCrossNumber = header.fBunchCrossNumber;
    fTriggerCluster = header.fTriggerCluster;
    SetName(header.fName);
    SetTitle(header.fTitle);

  } 
  return *this;
}



//______________________________________________________________________________
void AliESDHeader::Reset()
{
  // reset all data members
  fTriggerMask       = 0;
  fOrbitNumber       = 0;
  fTimeStamp         = 0;
  fEventType         = 0;
  fPeriodNumber      = 0;
  fEventNumberInFile = 0;
  fBunchCrossNumber  = 0;
  fTriggerCluster    = 0;
}

//______________________________________________________________________________
void AliESDHeader::Print(const Option_t *) const
{
  // Print some data members
  printf("Event # %d in file Bunch crossing # %d Orbit # %d Trigger %lld \n",
	 GetEventNumberInFile(),
	 GetBunchCrossNumber(),
	 GetOrbitNumber(),
	 GetTriggerMask());
}

