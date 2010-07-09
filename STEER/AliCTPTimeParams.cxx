/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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


#include <Riostream.h>

#include <TObjArray.h>
#include <TObjString.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>

#include "AliLog.h"
#include "AliCTPTimeParams.h"
#include "AliCTPInputTimeParams.h"

ClassImp(AliCTPTimeParams)

//______________________________________________________________________________
AliCTPTimeParams::AliCTPTimeParams():
TNamed(),
fDelayL1L0(0),
fDelayL2L0(0),
fCTPInputTimeParams()
{
//Default constructor
}

//______________________________________________________________________________
AliCTPTimeParams::AliCTPTimeParams(const AliCTPTimeParams &timeparams):
 TNamed(),
 fDelayL1L0(timeparams.fDelayL1L0),
 fDelayL2L0(timeparams.fDelayL2L0),
 fCTPInputTimeParams()
{
 for (Int_t i = 0; i < timeparams.fCTPInputTimeParams.GetSize(); i++) {
 if ( timeparams.fCTPInputTimeParams[i] ) fCTPInputTimeParams.Add(timeparams.fCTPInputTimeParams[i]->Clone());
 }
 // copy constructor
}


//______________________________________________________________________________
AliCTPTimeParams& AliCTPTimeParams::operator=(const AliCTPTimeParams &timeparams)
{
 // assignment operator
 if(this==&timeparams) return *this;
 ((TNamed *)this)->operator=(timeparams);
 fDelayL1L0=timeparams.fDelayL1L0;
 fDelayL2L0=timeparams.fDelayL2L0;
 
 for (Int_t i = 0; i < timeparams.fCTPInputTimeParams.GetSize(); i++) {
 if ( timeparams.fCTPInputTimeParams[i] ) fCTPInputTimeParams.Add(timeparams.fCTPInputTimeParams[i]->Clone());
 }

 return *this;
}

//______________________________________________________________________________
AliCTPTimeParams::~AliCTPTimeParams()
{
 //Destructor

 fCTPInputTimeParams.SetOwner();
 fCTPInputTimeParams.Delete();
}

//______________________________________________________________________________
void AliCTPTimeParams::AddInput( TString& inputName, UInt_t& inputLevel, UInt_t inputDelay, TString inputEdge, UInt_t deltamin, UInt_t deltamax )
{
 fCTPInputTimeParams.AddLast( new AliCTPInputTimeParams(inputName, inputLevel, inputDelay, inputEdge, deltamin, deltamax ));
}

//______________________________________________________________________________
void AliCTPTimeParams::AddDelayL0L1L2(Int_t delayL1L0, UInt_t delayL2L0)
{
 fDelayL1L0 = delayL1L0;
 fDelayL2L0 = delayL2L0;
}
//______________________________________________________________________________
AliCTPInputTimeParams* AliCTPTimeParams::GetTimeParamsForInput( TString inputname)
{
 // Get AliCTPInputTimeParams for input name
Int_t ninputs = fCTPInputTimeParams.GetEntriesFast();
for ( Int_t i=0; i < ninputs; i++ )
 {
  AliCTPInputTimeParams* ctpinputtime = (AliCTPInputTimeParams*)fCTPInputTimeParams.At(i);
  if (inputname == ctpinputtime->GetInputName() ) return ctpinputtime;
 }
cout << "Input: " << inputname << " not found." << endl;
return NULL;
}
//______________________________________________________________________________
Int_t AliCTPTimeParams::GetDeltasforClass(TString classname,Int_t& deltamin,Int_t& deltamax)
{
 // Get deltamin and deltamax for given class
 // Assumes that descriptor = DINPU
 // ret: 0=OK, 1= input doea not exist
 TString input(classname(1,4)); 
 AliCTPInputTimeParams* tprm = GetTimeParamsForInput(input);
 if(tprm){
  deltamin=tprm->GetDeltaMin();
  deltamax=tprm->GetDeltaMax();
  return 0;
 }
 return 1;
}
//______________________________________________________________________________
AliCTPTimeParams* AliCTPTimeParams::LoadCTPTimeParams(TString filename)
{
 // Load pre-created CTP time parameters from database/file
 // By default files are stored in GRP/CTP folder
 // The filename is constructed as GRP/CTP/<file>.cfg
  if( gSystem->AccessPathName( filename.Data() )) {
 // AliError( Form( "File (%s) not found!", filename.Data()));
  return NULL;
 }

 ifstream *file = new ifstream( filename.Data() );
 if(!*file) {
  //AliErrorClass( Form( "Error opening file (%s) !", filename.Data()));
  file->close();
  delete file;
  return NULL;
 }

 AliCTPTimeParams *ctptime = new AliCTPTimeParams();

 TString strline;
 
 while (strline.ReadLine(*file)) {
    if (ctptime->ProcessCTPTimeParamsLine(strline) == kFALSE) {
     delete ctptime;
     break;
    }
 }

 file->close();
 delete file;

 return ctptime;
}

//______________________________________________________________________________
AliCTPTimeParams* AliCTPTimeParams::LoadCTPTimeParamsFromString(const char* timeparams)
{

 // Loads configuration from string

   if (!timeparams)
     return 0;

   AliCTPTimeParams *ctptime = new AliCTPTimeParams();

   TObjArray* tokens = TString(timeparams).Tokenize("\n");
   for (Int_t i=0; i<tokens->GetEntries(); i++)
   {
     TObjString* string = dynamic_cast<TObjString*>(tokens->At(i));
     if (!string)
       continue;

     if (ctptime->ProcessCTPTimeParamsLine(string->String()) == kFALSE)
     {
        ctptime  = 0x0;
        break;
     }
   }

   delete tokens;
   if (ctptime) return ctptime;
   else return NULL;
}

//______________________________________________________________________________
Bool_t AliCTPTimeParams::ProcessCTPTimeParamsLine(const char* line)
{
 UInt_t level = 0;
 TString strline(line);
  if (strline.BeginsWith("L012")) {
   strline.ReplaceAll("L012", "");
   TObjArray *tokens = strline.Tokenize(" \t");
      

   AddDelayL0L1L2(((TObjString*)tokens->At(0))->String().Atoi(),((TObjString*)tokens->At(1))->String().Atoi());
   delete tokens;
  }
  else {
    if (strline.BeginsWith("0")) { level = 0; }   // determine the input level (0, 1 or 2)
    else if (strline.BeginsWith("1")) { level = 1; }
    else if (strline.BeginsWith("2")) { level = 2; } 
    else return kFALSE; // file not in the right format!
   
    TObjArray *tokens = strline.Tokenize(" \t");
    Int_t ntokens = tokens->GetEntriesFast();
    if (ntokens == 5)  AddInput(((TObjString*)tokens->At(0))->String(), level, ((TObjString*)tokens->At(2))->String().Atoi(), ((TObjString*)tokens->At(1))->String(),  ((TObjString*)tokens->At(3))->String().Atoi(), ((TObjString*)tokens->At(4))->String().Atoi());
    else if (ntokens == 3) AddInput(((TObjString*)tokens->At(0))->String(), level, ((TObjString*)tokens->At(2))->String().Atoi(), ((TObjString*)tokens->At(1))->String(),  0, 0); //the old format is used - no DeltaMin & DeltaMax!
         else return kFALSE; // file not in the right format!
  }
return kTRUE;
}

//______________________________________________________________________________
void AliCTPTimeParams::Print(const Option_t*) const
{
  //Print
 cout << "Delay L0 - L1 = " << fDelayL1L0 << endl;
 cout << "Delay L0 - L2 = " << fDelayL2L0 << endl;
 fCTPInputTimeParams.Print();
 cout << ""<<endl;
}
//______________________________________________________________________________
