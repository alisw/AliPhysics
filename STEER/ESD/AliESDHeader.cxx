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
#include "AliTriggerConfiguration.h"
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
  fTriggerScalers(),
  fTriggerScalersDeltaEvent(),
  fTriggerScalersDeltaRun(),
  fTriggerInputsNames(kNTriggerInputs),
  fCTPConfig(NULL),
  fIRBufferArray()
{
  // default constructor

  SetName("AliESDHeader");
  for(Int_t i = 0; i<kNMaxIR ; i++) fIRArray[i] = 0;
  fTriggerInputsNames.SetOwner(kTRUE);
}

AliESDHeader::~AliESDHeader() 
{
  // destructor
  for(Int_t i=0;i<kNMaxIR;i++)if(fIRArray[i])delete fIRArray[i];
  delete fCTPConfig;

  fIRBufferArray.Delete();
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
  fTriggerScalers(header.fTriggerScalers),
  fTriggerScalersDeltaEvent(header.fTriggerScalersDeltaEvent),
  fTriggerScalersDeltaRun(header.fTriggerScalersDeltaRun),
  fTriggerInputsNames(TObjArray(kNTriggerInputs)),
  fCTPConfig(header.fCTPConfig),
  fIRBufferArray()
{
  // copy constructor
  SetName(header.fName);
  SetTitle(header.fTitle);
  for(Int_t i = 0; i<kNMaxIR ; i++) {
    if(header.fIRArray[i])fIRArray[i] = new AliTriggerIR(*header.fIRArray[i]);
    else fIRArray[i]=0;
  }
  for(Int_t i = 0; i < kNTriggerInputs; i++) {
    TNamed *str = (TNamed *)((header.fTriggerInputsNames).At(i));
    if (str) fTriggerInputsNames.AddAt(new TNamed(*str),i);
  }

  for(Int_t i = 0; i < (header.fIRBufferArray).GetEntries(); ++i) {
    AliTriggerIR *ir = (AliTriggerIR*)((header.fIRBufferArray).At(i));
    if (ir) fIRBufferArray.Add(new AliTriggerIR(*ir));
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
    fTriggerScalersDeltaEvent = header.fTriggerScalersDeltaEvent;
    fTriggerScalersDeltaRun = header.fTriggerScalersDeltaRun;
    fCTPConfig = header.fCTPConfig;

    fTriggerInputsNames.Clear();
    for(Int_t i = 0; i < kNTriggerInputs; i++) {
      TNamed *str = (TNamed *)((header.fTriggerInputsNames).At(i));
      if (str) fTriggerInputsNames.AddAt(new TNamed(*str),i);
    }
    for(Int_t i = 0; i<kNMaxIR ; i++) {
       if(header.fIRArray[i])fIRArray[i] = new AliTriggerIR(*header.fIRArray[i]);
       else fIRArray[i]=0;
    }
    SetName(header.fName);
    SetTitle(header.fTitle);

    fIRBufferArray.Delete();
    for(Int_t i = 0; i < (header.fIRBufferArray).GetEntries(); ++i) {
      AliTriggerIR *ir = (AliTriggerIR*)((header.fIRBufferArray).At(i));
      if (ir) fIRBufferArray.Add(new AliTriggerIR(*ir));
    }
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
  fTriggerScalersDeltaEvent.Reset();
  fTriggerScalersDeltaRun.Reset();
  fTriggerInputsNames.Clear();
  delete fCTPConfig;
  fCTPConfig = 0;
  for(Int_t i=0;i<kNMaxIR;i++)if(fIRArray[i]){
   delete fIRArray[i];
   fIRArray[i]=0;
  }

  fIRBufferArray.Delete();
}
//______________________________________________________________________________
Bool_t AliESDHeader::AddTriggerIR(const AliTriggerIR* ir)
{
  // Add an IR object into the array
  // of IRs in the ESD header

 fIRBufferArray.Add(new AliTriggerIR(*ir));

 return kTRUE;
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
         printf("List of the active trigger inputs: ");
  	 for(Int_t i = 0; i < kNTriggerInputs; i++) {
    	   TNamed *str = (TNamed *)((fTriggerInputsNames).At(i));
    	   if (str) printf("%s ",str->GetName());
         }
         printf("\n");
}

//______________________________________________________________________________
void AliESDHeader::SetActiveTriggerInputs(const char*name, Int_t index)
{
  // Fill the active trigger inputs names
  // into the corresponding fTriggerInputsNames (TObjArray of TNamed)
  if (index >= kNTriggerInputs || index < 0) {
    AliError(Form("Index (%d) is outside the allowed range (0,59)!",index));
    return;
  }

  fTriggerInputsNames.AddAt(new TNamed(name,NULL),index);
}
//______________________________________________________________________________
const char* AliESDHeader::GetTriggerInputName(Int_t index, Int_t trglevel) const
{
  // Get the trigger input name
  // at the specified position in the trigger mask and trigger level (0,1,2)
  TNamed *trginput = 0;
  if (trglevel == 0) trginput = (TNamed *)fTriggerInputsNames.At(index);
  if (trglevel == 1) trginput = (TNamed *)fTriggerInputsNames.At(index+24);  
  if (trglevel == 2) trginput = (TNamed *)fTriggerInputsNames.At(index+48); 
  if (trginput) return trginput->GetName();
  else return "";
}
//______________________________________________________________________________
TString AliESDHeader::GetActiveTriggerInputs() const
{
  // Returns the list with the names of the active trigger inputs
  TString trginputs;
  for(Int_t i = 0; i < kNTriggerInputs; i++) {
    TNamed *str = (TNamed *)((fTriggerInputsNames).At(i));
    if (str) {
      trginputs += " ";
      trginputs += str->GetName();
      trginputs += " ";
    }
  }

  return trginputs;
}
//______________________________________________________________________________
TString AliESDHeader::GetFiredTriggerInputs() const
{
  // Returns the list with the names of the fired trigger inputs
  TString trginputs;
  for(Int_t i = 0; i < kNTriggerInputs; i++) {
      TNamed *str = (TNamed *)((fTriggerInputsNames.At(i)));
      if (i < 24 && (fL0TriggerInputs & (1 << i))) {
        if (str) {
	  trginputs += " ";
	  trginputs += str->GetName();
          trginputs += " ";
        }
      }
      if (i >= 24 && i < 48 && (fL1TriggerInputs & (1 << (i-24)))) {
        if (str) {
	  trginputs += " ";
	  trginputs += str->GetName();
          trginputs += " ";
        }
      }
      if (i >= 48 && (fL2TriggerInputs & (1 << (i-48)))) {
        if (str) {
	  trginputs += " ";
	  trginputs += str->GetName();
          trginputs += " ";
        }
      }

  }
  return trginputs;
}
//______________________________________________________________________________
Bool_t AliESDHeader::IsTriggerInputFired(const char *name) const
{
  // Checks if the trigger input is fired 
 
  TNamed *trginput = (TNamed *)fTriggerInputsNames.FindObject(name);
  if (!trginput) return kFALSE;

  Int_t inputIndex = fTriggerInputsNames.IndexOf(trginput);
  if (inputIndex < 0) return kFALSE;
  
  if (fL0TriggerInputs & (1 << inputIndex)) return kTRUE;
  else if (fL1TriggerInputs & (1 << (inputIndex-24))) return kTRUE;
  else if (fL2TriggerInputs & (1 << (inputIndex-48))) return kTRUE;
  else return kFALSE;
}
