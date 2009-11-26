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
#include "AliTriggerScalersESD.h"
#include "AliTriggerScalersRecordESD.h"
#include "AliTriggerIR.h"
#include "AliLog.h" 

ClassImp(AliESDHeader)

//______________________________________________________________________________
AliESDHeader::AliESDHeader() :
  AliVHeader(),
  fTriggerMask(0),
  fOrbitNumber(0),
  fTimeStamp(0),
  fEventType(0),
  fEventSpecie(0),
  fPeriodNumber(0),
  fEventNumberInFile(0),
  fBunchCrossNumber(0),
  fTriggerCluster(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fTriggerScalers()
{
  // default constructor

  SetName("AliESDHeader");
  for(Int_t i = 0; i<kNMaxIR ; i++) fIRArray[i] = 0;

}
AliESDHeader::~AliESDHeader() 
{
  // destructor
  for(Int_t i=0;i<kNMaxIR;i++)if(fIRArray[i])delete fIRArray[i];
}


AliESDHeader::AliESDHeader(const AliESDHeader &header) :
  AliVHeader(header),
  fTriggerMask(header.fTriggerMask),
  fOrbitNumber(header.fOrbitNumber),
  fTimeStamp(header.fTimeStamp),
  fEventType(header.fEventType),
  fEventSpecie(header.fEventSpecie),
  fPeriodNumber(header.fPeriodNumber),
  fEventNumberInFile(header.fEventNumberInFile),
  fBunchCrossNumber(header.fBunchCrossNumber),
  fTriggerCluster(header.fTriggerCluster),
  fL0TriggerInputs(header.fL0TriggerInputs),
  fL1TriggerInputs(header.fL1TriggerInputs),
  fL2TriggerInputs(header.fL2TriggerInputs),
  fTriggerScalers(header.fTriggerScalers)
{
  // copy constructor
  SetName(header.fName);
  SetTitle(header.fTitle);
  for(Int_t i = 0; i<kNMaxIR ; i++) {
    if(header.fIRArray[i])fIRArray[i] = new AliTriggerIR(*header.fIRArray[i]);
    else fIRArray[i]=0;
  }
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
    fEventSpecie = header.fEventSpecie;
    fPeriodNumber = header.fPeriodNumber;
    fEventNumberInFile = header.fEventNumberInFile;
    fBunchCrossNumber = header.fBunchCrossNumber;
    fTriggerCluster = header.fTriggerCluster;
    fL0TriggerInputs = header.fL0TriggerInputs;
    fL1TriggerInputs = header.fL1TriggerInputs;
    fL2TriggerInputs = header.fL2TriggerInputs;
    fTriggerScalers = header.fTriggerScalers;
    for(Int_t i = 0; i<kNMaxIR ; i++) {
       if(header.fIRArray[i])fIRArray[i] = new AliTriggerIR(*header.fIRArray[i]);
       else fIRArray[i]=0;
    }
    SetName(header.fName);
    SetTitle(header.fTitle);

  } 
  return *this;
}

void AliESDHeader::Copy(TObject &obj) const 
{  
  // this overwrites the virtual TOBject::Copy()
  // to allow run time copying without casting
  // in AliESDEvent

  if(this==&obj)return;
  AliESDHeader *robj = dynamic_cast<AliESDHeader*>(&obj);
  if(!robj)return; // not an AliESDHeader
  *robj = *this;

}
//______________________________________________________________________________
void AliESDHeader::Reset()
{
  // reset all data members
  fTriggerMask       = 0;
  fOrbitNumber       = 0;
  fTimeStamp         = 0;
  fEventType         = 0;
  fEventSpecie       = 0;
  fPeriodNumber      = 0;
  fEventNumberInFile = 0;
  fBunchCrossNumber  = 0;
  fTriggerCluster    = 0;
  fL0TriggerInputs   = 0;
  fL1TriggerInputs   = 0;
  fL2TriggerInputs   = 0;
  fTriggerScalers.Reset();
  for(Int_t i=0;i<kNMaxIR;i++)if(fIRArray[i]){
   delete fIRArray[i];
   fIRArray[i]=0;
  }
}
//______________________________________________________________________________
Bool_t AliESDHeader::AddTriggerIR(const AliTriggerIR* ir)
{
 // Adds trigger interaction record to array
 for(Int_t i=0;i<kNMaxIR;i++){
  if(!fIRArray[i]){
    fIRArray[i]=const_cast<AliTriggerIR*>(ir);
    return 0;
  }
 }
 //AliErrorClass("Attempt to add # of IRs > kNMaxIR \n");
 return 1;
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

